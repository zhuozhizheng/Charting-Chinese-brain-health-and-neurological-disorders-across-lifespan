# -*- coding: utf-8 -*-
"""Copy built reports and interpretation artifacts into the series' sibling
_sub folder (the folder freesurf_luhao.py creates next to the DICOM series,
which also receives the full tem_data scratch copy during the pipeline run)."""
import os
import shutil

from . import reports

ARTIFACTS = ['interpretations.json', 'subject_meta.json']


def series_sub_dir(series_dir):
    return os.path.abspath(series_dir).rstrip('/') + '_sub'


def export_to_sub(work_dir, series_dir):
    dst = series_sub_dir(series_dir)
    os.makedirs(dst, exist_ok=True)
    copied = []
    for name in reports.OUTPUTS + ARTIFACTS:
        src = os.path.join(work_dir, name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dst, name))
            copied.append(name)
    return dst, copied
