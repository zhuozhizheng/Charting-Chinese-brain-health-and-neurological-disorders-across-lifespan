# -*- coding: utf-8 -*-
"""Optional LLM step: per-subject bilingual medical/neuroscience interpretations.

Config lives in brain_agent/llm_config.json (server-side only):
{"base_url": "...", "api_key": "...", "model": "...", "timeout": 120}
When configured, the runner asks the model for a strict-JSON interpretation set
before building the HTML reports; build_reports.py then renders it.
"""
import csv
import json
import os
import shutil
import subprocess
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(HERE, 'llm_config.json')
QODER_TIMEOUT = 900

GLOBAL = [
    ('TCV', 'intracranial total volume (ml)', 1000.0),
    ('GMV', 'whole-brain grey-matter volume (ml)', 1000.0),
    ('sGMV', 'subcortical grey-matter volume (ml)', 1000.0),
    ('WMV', 'white-matter volume (ml)', 1000.0),
    ('Ventricles', 'ventricular system volume (ml)', 1000.0),
    ('mean_thickness', 'mean cortical thickness (mm)', 1.0),
    ('total_surface_arrea', 'total cortical surface area (cm2)', 100.0),
    ('cerebellum_total', 'total cerebellar volume (ml)', 1000.0),
    ('Brain.Stem', 'brainstem volume (ml)', 1000.0),
]
SUBCORT = ['Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala',
           'Accumbens.area', 'VentralDC']

KEYS_CURVES = ['GMV', 'sGMV', 'WMV', 'Ventricles', 'mean_thickness',
               'total_surface_area', 'Cerebellum', 'BrainStem']
KEYS_SUB = [k + '_sub' for k in KEYS_CURVES]
KEYS_LISTS = {'discuss64': 3, 'sub62': 4}
KEYS_TEXT = KEYS_CURVES + KEYS_SUB + ['discuss61', 'discuss62', 'discuss63',
                                      'sub61', 'sub63'] + list(KEYS_LISTS)


def load_config():
    if not os.path.exists(CONFIG_PATH):
        return None
    with open(CONFIG_PATH, encoding='utf-8') as fh:
        cfg = json.load(fh)
    if not (cfg.get('base_url') and cfg.get('api_key') and cfg.get('model')):
        return None
    cfg.setdefault('timeout', 180)
    return cfg


def save_config(cfg):
    with open(CONFIG_PATH, 'w', encoding='utf-8') as fh:
        json.dump(cfg, fh, ensure_ascii=False, indent=1)
    os.chmod(CONFIG_PATH, 0o600)


def masked(cfg):
    out = dict(cfg or {})
    key = out.get('api_key') or ''
    out['api_key_set'] = bool(key)
    out['api_key_masked'] = (key[:4] + '***' + key[-4:]) if len(key) > 8 else ('***' if key else '')
    out.pop('api_key', None)
    return out


def _flag(c):
    if c < 0.05:
        return 'markedly low'
    if c < 0.15:
        return 'low'
    if c > 0.95:
        return 'markedly high'
    if c > 0.85:
        return 'high'
    return 'mid'


def build_summary(work_dir):
    def row(path):
        with open(path, newline='', encoding='utf-8') as fh:
            rows = list(csv.reader(fh))
        return {h: v for h, v in zip(rows[0], rows[1])}

    def clean(d):
        return {k: v for k, v in d.items()
                if v is not None and str(v).strip() not in ('', 'NA', 'NaN')}

    legacy_p = os.path.join(work_dir, 'individual_centile_score.csv')
    base = row(legacy_p) if os.path.exists(legacy_p) else {}

    def load(primary):
        """Norm-specific centiles merged over the legacy export (which fills
        metrics the norm-specific pivot omits); NA/empty values are ignored."""
        p = os.path.join(work_dir, primary)
        if os.path.exists(p):
            merged = dict(base)
            merged.update(clean(row(p)))
            return merged
        return dict(base)

    raw = row(os.path.join(work_dir, 'individual_raw_data.csv'))
    cen = load('individual_centile_score_overall_population.csv')
    sub_p = os.path.join(work_dir, 'individual_centile_score_subtype_specific.csv')
    has_sub = os.path.exists(sub_p)
    cen_sub = load('individual_centile_score_subtype_specific.csv') if has_sub else None

    def block(cen):
        out = {'global': [], 'subcortical': [], 'cortical': {}}
        for key, label, div in GLOBAL:
            v = float(raw[key]) / div
            c = float(cen[key])
            out['global'].append({'metric': label, 'value': round(v, 3),
                                  'centile': round(c, 4), 'flag': _flag(c)})
        for side in ('Left', 'Right'):
            for r in SUBCORT:
                k = '%s.%s' % (side, r)
                c = float(cen[k])
                out['subcortical'].append({'structure': k,
                                           'volume_ml': round(float(raw[k]) / 1000.0, 2),
                                           'centile': round(c, 4), 'flag': _flag(c)})
        th = [float(cen[h]) for h in cen if h.endswith('_thickness') and h[:3] in ('lh_', 'rh_')]
        ar = [float(cen[h]) for h in cen if h.endswith('_area') and h[:3] in ('lh_', 'rh_')]
        vo = [float(cen[h]) for h in cen if h.endswith('_volume') and h[:3] in ('lh_', 'rh_')]
        low_th = sorted([h for h in cen if h.endswith('_thickness') and h[:3] in ('lh_', 'rh_')],
                        key=lambda h: float(cen[h]))[:6]
        out['cortical'] = {
            'n_parcels': len(th),
            'thickness_below_0.05': sum(1 for v in th if v < 0.05),
            'thickness_below_0.10': sum(1 for v in th if v < 0.10),
            'thickness_mean_centile': round(sum(th) / len(th), 3),
            'area_mean_centile': round(sum(ar) / len(ar), 3),
            'volume_below_0.10': sum(1 for v in vo if v < 0.10),
            'lowest_thickness_parcels': low_th,
        }
        return out

    summary = block(cen)
    if cen_sub is not None:
        summary['subtype_specific'] = block(cen_sub)
    return summary


SYSTEM_PROMPT = (
    'You are a clinical neuroscience report writer for structural-MRI morphometry '
    '(FreeSurfer-style centiles, 0-1). You receive one subject\'s quantitative summary as JSON: '
    'the top-level global/subcortical/cortical centiles use the population-general norm, and when '
    'a "subtype_specific" object is present it repeats the same metrics with centiles under the '
    'subtype-specific normative model. '
    'Return STRICT JSON only, with exactly two top-level objects "zh" and "en". Each contains '
    'these keys: %s. Keys discuss64 and sub62 are arrays of strings (3 and 4 items). '
    'All other keys are single strings. zh in Simplified Chinese, en in English, same content. '
    'Curve keys describe that metric for THIS subject (value, centile, flag, trajectory position, '
    'medical and neuroscience significance) using the population-general norm; *_sub keys compare '
    'the subtype-specific norm (summary.subtype_specific centiles) with the general norm '
    '(top-level centiles), quoting both; discuss61 phenotype summary; discuss62 neurodevelopmental '
    'perspective; discuss63 '
    'technical/measurement caveats; discuss64 clinical associations and recommendations; '
    'sub61 what the subtype comparison adds; sub62 mechanism and clinical implications; '
    'sub63 interpretation caveats. Quote the given numbers exactly. No diagnosis, '
    'research-use framing only.' % ', '.join(KEYS_TEXT))


def _call(cfg, summary):
    body = json.dumps({'model': cfg['model'],
                       'messages': [{'role': 'system', 'content': SYSTEM_PROMPT},
                                    {'role': 'user', 'content': json.dumps(summary, ensure_ascii=False)}],
                       'temperature': 0.3}).encode()
    req = urllib.request.Request(cfg['base_url'].rstrip('/') + '/chat/completions',
                                 data=body,
                                 headers={'Content-Type': 'application/json',
                                          'Authorization': 'Bearer ' + cfg['api_key']})
    with urllib.request.urlopen(req, timeout=cfg['timeout']) as resp:
        out = json.load(resp)
    return out['choices'][0]['message']['content']


def validate(text):
    data = json.loads(_extract_json(text))
    for lang in ('zh', 'en'):
        node = data[lang]
        for key in KEYS_TEXT:
            val = node[key]
            if key in KEYS_LISTS:
                if not isinstance(val, list) or len(val) != KEYS_LISTS[key] or not all(
                        isinstance(x, str) and x.strip() for x in val):
                    raise ValueError('bad list for %s/%s' % (lang, key))
            elif not isinstance(val, str) or not val.strip():
                raise ValueError('missing %s/%s' % (lang, key))
    return data


def _extract_json(text):
    text = text.strip()
    if text.startswith('```'):
        text = text.strip('`')
        if text.startswith('json'):
            text = text[4:]
    start, end = text.find('{'), text.rfind('}')
    if start == -1 or end <= start:
        raise ValueError('no JSON object in model output')
    return text[start:end + 1]


def qoder_binary():
    path = shutil.which('qoder')
    if path:
        return path
    home = os.path.expanduser('~/.qoder/entry/qoder')
    return home if os.path.exists(home) else None


def _call_qoder(summary):
    prompt = ('Subject quantitative summary JSON:\n'
              + json.dumps(summary, ensure_ascii=False)
              + '\nReturn STRICT JSON only, following the system instructions.')
    cmd = [qoder_binary(), '-p', prompt,
           '--system-prompt', SYSTEM_PROMPT,
           '--tools', '',
           '--no-session-persistence',
           '--permission-mode', 'dont_ask']
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=QODER_TIMEOUT)
    if proc.returncode != 0:
        raise RuntimeError('qoder cli exit %d: %s'
                           % (proc.returncode, (proc.stderr or proc.stdout)[-400:]))
    return proc.stdout


def write_interpretations(work_dir):
    """Returns (ok, error).

    Engine priority: configured OpenAI-compatible endpoint first, then the
    Qoder CLI (default, no key needed). When neither is available the step is
    skipped and the report builder falls back to its frozen demo text with a
    banner. When an engine is available but fails, the subject fails.
    """
    cfg = load_config()
    engines = []
    if cfg:
        engines.append(('llm', lambda s: _call(cfg, s)))
    if qoder_binary():
        engines.append(('qoder', _call_qoder))
    if not engines:
        return True, None
    try:
        summary = build_summary(work_dir)
    except Exception as exc:
        return False, 'summary build failed: %s' % exc
    last = None
    for name, call in engines:
        for _attempt in range(2):
            try:
                data = validate(call(summary))
                with open(os.path.join(work_dir, 'interpretations.json'), 'w',
                          encoding='utf-8') as fh:
                    json.dump(data, fh, ensure_ascii=False, indent=1)
                return True, None
            except Exception as exc:
                last = '%s: %s' % (name, exc)
    return False, 'interpretation generation failed: %s' % last


def ping(cfg):
    body = json.dumps({'model': cfg['model'],
                       'messages': [{'role': 'user', 'content': 'ping'}],
                       'max_tokens': 5}).encode()
    req = urllib.request.Request(cfg['base_url'].rstrip('/') + '/chat/completions',
                                 data=body,
                                 headers={'Content-Type': 'application/json',
                                          'Authorization': 'Bearer ' + cfg['api_key']})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)['choices'][0]['message']['content']
