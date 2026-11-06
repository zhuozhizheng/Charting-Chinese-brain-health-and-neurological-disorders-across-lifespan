# -*- coding: utf-8 -*-
"""Serial run orchestrator with JSON state, resume support and a global lock
(the pipeline scratch directory is shared and overwrite-style)."""
import json
import os
import subprocess
import threading
import time
import uuid

from . import discovery, export, llm, pipeline, reports, snapshot

TERMINAL = ('done', 'failed', 'blocked')


def external_pipeline_running():
    try:
        out = subprocess.run(['pgrep', '-f', r'freesurf_luhao\.py --input='],
                             capture_output=True, text=True).stdout.split()
    except OSError:
        return False
    return bool(out)


class Runner:
    def __init__(self, cfg):
        self.cfg = cfg
        self.out_root = cfg['out_root']
        os.makedirs(self.out_root, exist_ok=True)

    def run_dir(self, run_id):
        return os.path.join(self.out_root, run_id)

    def state_path(self, run_id):
        return os.path.join(self.run_dir(run_id), 'state.json')

    def load_state(self, run_id):
        with open(self.state_path(run_id), encoding='utf-8') as fh:
            return json.load(fh)

    def save_state(self, state):
        path = self.state_path(state['run_id'])
        tmp = path + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(state, fh, ensure_ascii=False, indent=1)
        os.replace(tmp, path)

    def list_runs(self):
        out = []
        for d in sorted(os.listdir(self.out_root)):
            if os.path.exists(os.path.join(self.out_root, d, 'state.json')):
                out.append(d)
        return out

    def start_run(self, root, subjects=None, dry_run=False, run_id=None):
        pairs, skipped = discovery.discover(root)
        if subjects:
            want = set(subjects)
            pairs = [p for p in pairs if p[0] in want]
        run_id = run_id or time.strftime('%Y%m%d-%H%M%S-') + uuid.uuid4().hex[:6]
        state = {'run_id': run_id, 'root': os.path.abspath(root),
                 'created': time.time(), 'dry_run': bool(dry_run),
                 'skipped': skipped,
                 'subjects': [{'subject': s, 'series': d, 'status': 'planned'}
                              for s, d in pairs]}
        os.makedirs(self.run_dir(run_id), exist_ok=True)
        self.save_state(state)
        if not dry_run:
            threading.Thread(target=self._execute, args=(run_id,), daemon=True).start()
        return run_id

    def finish_subject(self, series_dir, run_id=None):
        """Report-only path for a series already processed by the external
        pipeline (run_freesurf_mulu.sh): skips the pipeline stage and accepts
        scratch files regardless of age."""
        series_dir = os.path.abspath(series_dir.rstrip('/'))
        parent = os.path.dirname(series_dir)
        name = discovery._safe(os.path.basename(parent) or os.path.basename(series_dir))
        run_id = run_id or time.strftime('finish-%Y%m%d-%H%M%S-') + uuid.uuid4().hex[:6]
        state = {'run_id': run_id, 'root': parent,
                 'created': time.time(), 'dry_run': False, 'skipped': [],
                 'subjects': [{'subject': name, 'series': series_dir,
                               'status': 'planned', 'skip_pipeline': True}]}
        os.makedirs(self.run_dir(run_id), exist_ok=True)
        self.save_state(state)
        threading.Thread(target=self._execute, args=(run_id,), daemon=True).start()
        return run_id

    def _acquire_lock(self):
        path = os.path.join(self.out_root, 'run.lock')
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
        return path

    def _release_lock(self, path):
        try:
            os.remove(path)
        except OSError:
            pass

    def _execute(self, run_id):
        try:
            lock = self._acquire_lock()
        except FileExistsError:
            state = self.load_state(run_id)
            state['error'] = 'another run holds %s/run.lock' % self.out_root
            for s in state['subjects']:
                if s['status'] == 'planned':
                    s['status'] = 'blocked'
            self.save_state(state)
            return
        try:
            state = self.load_state(run_id)
            for entry in state['subjects']:
                if entry['status'] == 'done':
                    continue
                self._process(state, entry)
                self.save_state(state)
        finally:
            self._release_lock(lock)

    def _process(self, state, entry):
        work = os.path.join(self.run_dir(state['run_id']), entry['subject'])
        os.makedirs(work, exist_ok=True)
        entry.update(status='running', dir=work, started=time.time())
        self.save_state(state)
        if entry.get('skip_pipeline'):
            start_ts = 0
        else:
            if external_pipeline_running():
                entry.update(status='blocked', finished=time.time(),
                             error='external freesurf_luhao.py run in progress; '
                                   'wait for it to finish, then resume this run')
                return
            start_ts = time.time()
            log = os.path.join(work, 'pipeline.log')
            entry['pipeline_log'] = log
            rc = pipeline.run_pipeline(self.cfg, entry['series'], log)
            entry['pipeline_rc'] = rc
            if rc != 0:
                entry.update(status='failed', error='pipeline exit %d' % rc,
                             finished=time.time())
                return
        ok, missing, _copied = snapshot.snapshot(self.cfg, work, start_ts)
        if not ok:
            entry.update(status='failed',
                         error='snapshot missing: %s' % ','.join(missing),
                         finished=time.time())
            return
        snapshot.write_subject_meta(work, entry['series'])
        oki, lerr = llm.write_interpretations(work)
        if not oki:
            entry.update(status='failed', error=lerr, finished=time.time())
            return
        ok2, missing2 = reports.build_reports(self.cfg, work,
                                              os.path.join(work, 'build_reports.log'))
        entry['reports'] = [os.path.join(work, o) for o in reports.OUTPUTS] if ok2 else []
        entry['status'] = 'done' if ok2 else 'failed'
        if not ok2:
            entry['error'] = 'report build missing: %s' % ','.join(missing2)
        else:
            try:
                dst, copied = export.export_to_sub(work, entry['series'])
                entry['export_dir'] = dst
                entry['exported'] = copied
            except OSError as exc:
                entry['export_error'] = str(exc)
        entry['finished'] = time.time()
