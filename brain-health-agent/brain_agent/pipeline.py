# -*- coding: utf-8 -*-
import shlex
import subprocess


def run_pipeline(cfg, series_dir, log_path):
    cmd = ('source {sh} && conda activate {env} && python {script} --input={inp} '
           '--out_start_freesurfer={fs} --out_csv={sc}').format(
        sh=shlex.quote(cfg['conda_sh']), env=cfg['conda_env'],
        script=shlex.quote(cfg['pipeline_script']), inp=shlex.quote(series_dir),
        fs=shlex.quote(cfg['fs_out_dir']), sc=shlex.quote(cfg['scratch_dir']))
    with open(log_path, 'wb') as fh:
        proc = subprocess.run(['bash', '-c', cmd], stdout=fh, stderr=subprocess.STDOUT)
    return proc.returncode
