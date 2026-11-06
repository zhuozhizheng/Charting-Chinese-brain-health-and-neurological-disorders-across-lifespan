# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

OUTPUTS = ['Individualized-brain-health-report_population-normative.html',
           'Individualized-brain-health-report_subtype-specific.html',
           'Individualized-brain-health-report_population-normative_zh.html',
           'Individualized-brain-health-report_population-normative_en.html',
           'Individualized-brain-health-report_subtype-specific_zh.html',
           'Individualized-brain-health-report_subtype-specific_en.html']


def build_reports(cfg, work_dir, log_path):
    shutil.copy2(cfg['build_script'], os.path.join(work_dir, 'build_reports.py'))
    with open(log_path, 'wb') as fh:
        proc = subprocess.run([cfg['python_bin'], 'build_reports.py'], cwd=work_dir,
                              stdout=fh, stderr=subprocess.STDOUT)
    missing = [o for o in OUTPUTS if not os.path.exists(os.path.join(work_dir, o))]
    return proc.returncode == 0 and not missing, missing
