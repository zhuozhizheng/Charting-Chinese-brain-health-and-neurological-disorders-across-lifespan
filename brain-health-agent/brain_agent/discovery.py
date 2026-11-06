# -*- coding: utf-8 -*-
"""Locate DICOM series folders under a root, mirroring run_freesurf_mulu.sh
(find *.IMA/*.dcm -> dirname -> dedup) plus the pipeline's own naming rules
(basename must contain "t1" and must not end with "_sub")."""
import os
import re

SERIES_EXTS = ('.ima', '.dcm')


def _safe(name):
    return re.sub(r'[^\w.\-]+', '_', name).strip('_') or 'subject'


def discover(root):
    hits = set()
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith(SERIES_EXTS):
                hits.add(dirpath)
                break
    pairs, skipped, seen = [], [], {}
    for d in sorted(hits):
        base = os.path.basename(d)
        if base.endswith('_sub') or 't1' not in base.lower():
            skipped.append([d, 'name filter (needs "t1", not "_sub")'])
            continue
        parent = os.path.dirname(d)
        subject = os.path.basename(parent)
        if not subject or os.path.samefile(parent, root):
            subject = base
        subject = _safe(subject)
        seen[subject] = seen.get(subject, 0) + 1
        if seen[subject] > 1:
            subject = '%s_%d' % (subject, seen[subject])
        pairs.append((subject, d))
    return pairs, skipped
