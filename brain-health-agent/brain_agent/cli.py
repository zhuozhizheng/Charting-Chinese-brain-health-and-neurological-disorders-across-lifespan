# -*- coding: utf-8 -*-
import argparse
import json
import os
import time

from . import discovery, reports
from .config import load_config
from .runner import Runner, TERMINAL


def _wait(runner, run_id):
    while True:
        state = runner.load_state(run_id)
        if all(s['status'] in TERMINAL for s in state['subjects']):
            return state
        time.sleep(5)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='brain-agent',
                                 description='brain-health report agent')
    ap.add_argument('--config', help='path to config json')
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('discover', help='list DICOM series under a root (read-only)')
    p.add_argument('root')

    p = sub.add_parser('run', help='process a root folder end to end')
    p.add_argument('root')
    p.add_argument('--run-id')
    p.add_argument('--subjects', help='comma separated subject subset')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--wait', action='store_true', help='block until the run finishes')

    p = sub.add_parser('finish', help='build reports for a series already processed '
                                      'by run_freesurf_mulu.sh (skips the pipeline stage)')
    p.add_argument('series_dir')
    p.add_argument('--run-id')
    p.add_argument('--wait', action='store_true', help='block until finished')

    p = sub.add_parser('status', help='print run state json')
    p.add_argument('run_id')

    sub.add_parser('runs', help='list run ids')

    p = sub.add_parser('report-only', help='rebuild the two HTML reports in an existing snapshot dir')
    p.add_argument('work_dir')

    p = sub.add_parser('serve', help='start the Flask service')
    p.add_argument('--host', default='127.0.0.1')
    p.add_argument('--port', type=int, default=8765)

    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    runner = Runner(cfg)

    if a.cmd == 'discover':
        pairs, skipped = discovery.discover(a.root)
        print(json.dumps({'subjects': pairs, 'skipped': skipped},
                         ensure_ascii=False, indent=1))
    elif a.cmd == 'run':
        subjects = a.subjects.split(',') if a.subjects else None
        run_id = runner.start_run(a.root, subjects, a.dry_run, a.run_id)
        print(run_id)
        if a.wait and not a.dry_run:
            print(json.dumps(_wait(runner, run_id), ensure_ascii=False, indent=1))
    elif a.cmd == 'finish':
        run_id = runner.finish_subject(a.series_dir, a.run_id)
        print(run_id)
        if a.wait:
            print(json.dumps(_wait(runner, run_id), ensure_ascii=False, indent=1))
    elif a.cmd == 'status':
        print(json.dumps(runner.load_state(a.run_id), ensure_ascii=False, indent=1))
    elif a.cmd == 'runs':
        print('\n'.join(runner.list_runs()))
    elif a.cmd == 'report-only':
        work = os.path.abspath(a.work_dir)
        ok, missing = reports.build_reports(cfg, work,
                                            os.path.join(work, 'build_reports.log'))
        print('OK' if ok else 'MISSING: %s' % ','.join(missing))
        raise SystemExit(0 if ok else 1)
    elif a.cmd == 'serve':
        from .server import create_app
        create_app(cfg).run(host=a.host, port=a.port, threaded=True)


if __name__ == '__main__':
    main()
