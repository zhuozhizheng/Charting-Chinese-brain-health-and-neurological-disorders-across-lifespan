# -*- coding: utf-8 -*-
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(HERE)

DEFAULTS = {
    'conda_sh': os.path.expanduser('~/anaconda3/etc/profile.d/conda.sh'),
    'conda_env': 'py38',
    'pipeline_script': '',
    'scratch_dir': '',
    'fs_out_dir': '',
    'out_root': os.path.join(WORKSPACE, 'agent_out'),
    'build_script': os.path.join(WORKSPACE, 'report_builder', 'build_reports.py'),
    'python_bin': '/usr/bin/python3',
}


def load_config(path=None):
    cfg = dict(DEFAULTS)
    p = path or os.environ.get('BRAIN_AGENT_CONFIG') or os.path.join(HERE, 'config.json')
    if p and os.path.exists(p):
        with open(p, encoding='utf-8') as fh:
            cfg.update(json.load(fh))
    cfg['conda_sh'] = os.path.expanduser(cfg['conda_sh'])
    return cfg
