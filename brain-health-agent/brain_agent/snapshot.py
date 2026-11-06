# -*- coding: utf-8 -*-
"""Snapshot the shared scratch directory (tem_data) into a per-subject folder
before the next subject overwrites it."""
import csv
import datetime
import json
import os
import shutil

CURVES = ['GMV', 'sGMV', 'WMV', 'Ventricles', 'mean_thickness',
          'total_surface_area', 'Cerebellum', 'BrainStem']
MAPS = ['cortical_volume', 'cortical_thickness', 'cortical_area', 'subcortical_volume']
TABLES = ['global_summary', 'local_summary']


REQUIRED_CSVS = ['individual_raw_data.csv', 'individual_centile_score.csv']
# Per-norm individual centile exports: population-general and subtype-specific.
# Copied when the pipeline provides them; the report builder falls back to the
# legacy individual_centile_score.csv when a norm-specific export is absent.
OPTIONAL_CSVS = ['individual_centile_score_overall_population.csv',
                 'individual_centile_score_subtype_specific.csv']


def required_files():
    out = list(REQUIRED_CSVS)
    for stem in CURVES + MAPS + TABLES:
        out.append(stem + '.png')
        out.append(stem + '_Subtype_specific.png')
    return out


def snapshot(cfg, dst_dir, start_ts):
    os.makedirs(dst_dir, exist_ok=True)
    missing, copied = [], []
    for name in required_files():
        src = os.path.join(cfg['scratch_dir'], name)
        if not os.path.exists(src) or os.path.getmtime(src) < start_ts:
            missing.append(name)
            continue
        shutil.copy2(src, os.path.join(dst_dir, name))
        copied.append(name)
    for name in OPTIONAL_CSVS:
        src = os.path.join(cfg['scratch_dir'], name)
        if os.path.exists(src) and os.path.getmtime(src) >= start_ts:
            shutil.copy2(src, os.path.join(dst_dir, name))
            copied.append(name)
    return (not missing), missing, copied


def _parse_date(s):
    s = (s or '').strip()
    if len(s) < 8:
        return None
    try:
        return datetime.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    except ValueError:
        return None


def write_subject_meta(work_dir, series_dir):
    meta = {}
    t1 = os.path.join(os.path.dirname(series_dir),
                      os.path.basename(series_dir) + '_sub', 'T1.csv')
    if os.path.exists(t1):
        with open(t1, newline='', encoding='utf-8', errors='replace') as fh:
            row = next(csv.DictReader(fh), None) or {}
        meta['id'] = (row.get('PatientName') or '').strip()
        sex = (row.get('PatientSex') or '').strip().upper()
        if sex in ('1', 'F'):
            meta['sex'] = 'female'
        elif sex in ('0', 'M'):
            meta['sex'] = 'male'
        study = _parse_date(row.get('StudyDate'))
        birth = _parse_date(row.get('PatientBirthDate'))
        if study:
            meta['study_date'] = study.isoformat()
        if study and birth:
            meta['age'] = str((study - birth).days // 365)
        meta['site'] = (row.get('InstitutionName') or '').strip()
    with open(os.path.join(work_dir, 'subject_meta.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    return meta
