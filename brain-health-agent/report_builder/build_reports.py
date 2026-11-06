# -*- coding: utf-8 -*-
"""Build two self-contained individualized brain-health reports (population-normative
and subtype-specific normative) from the pipeline PNG outputs and individual CSVs."""
import base64
import csv
import io
import json
import os
import re
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None
BASE = os.path.dirname(os.path.abspath(__file__))
ASSET = os.path.join(BASE, 'report_assets')
os.makedirs(ASSET, exist_ok=True)

# ---------------------------------------------------------------- data loading
def read_one_row(path):
    with open(path, newline='', encoding='utf-8') as fh:
        rows = list(csv.reader(fh))
    hdr = rows[0]
    val = rows[1]
    return {h: v for h, v in zip(hdr, val)}

raw = read_one_row(os.path.join(BASE, 'individual_raw_data.csv'))

_legacy_path = os.path.join(BASE, 'individual_centile_score.csv')
cen_leg = read_one_row(_legacy_path) if os.path.exists(_legacy_path) else {}

def _clean(d):
    return {k: v for k, v in d.items()
            if v is not None and str(v).strip() not in ('', 'NA', 'NaN')}

def load_centile(primary):
    """Norm-specific individual centiles merged over the legacy export.

    The per-norm pivot exports omit metrics not modelled for that norm (e.g.
    cerebellum_GM/cerebellum_WM are excluded from the deviation models), so the
    legacy per-subject centiles fill those gaps and NA values are ignored. Falls
    back entirely to the legacy file when the norm-specific export is absent.
    Returns (row_dict, source_filename)."""
    path = os.path.join(BASE, primary)
    if os.path.exists(path):
        merged = dict(cen_leg)
        merged.update(_clean(read_one_row(path)))
        return merged, primary
    return dict(cen_leg), 'individual_centile_score.csv'

cen_pop, cen_pop_src = load_centile('individual_centile_score_overall_population.csv')
cen_sub, cen_sub_src = load_centile('individual_centile_score_subtype_specific.csv')
if not cen_pop or not cen_sub:
    raise SystemExit('no individual centile CSV found in %s' % BASE)

def f(x):
    return float(x)

def _build_global_metrics(cen):
    return [
    ('TCV',                '颅内总容积 TCV',        'ml',   f(raw['TCV']) / 1000.0,              f(cen['TCV']), 'Intracranial total volume TCV'),
    ('GMV',                '全脑灰质体积 GMV',      'ml',   f(raw['GMV']) / 1000.0,              f(cen['GMV']), 'Whole-brain grey-matter volume GMV'),
    ('sGMV',               '皮层下灰质体积 sGMV',   'ml',   f(raw['sGMV']) / 1000.0,             f(cen['sGMV']), 'Subcortical grey-matter volume sGMV'),
    ('WMV',                '白质体积 WMV',          'ml',   f(raw['WMV']) / 1000.0,              f(cen['WMV']), 'White-matter volume WMV'),
    ('Ventricles',         '脑室系统体积',          'ml',   f(raw['Ventricles']) / 1000.0,       f(cen['Ventricles']), 'Ventricular system volume'),
    ('mean_thickness',     '平均皮层厚度',          'mm',   f(raw['mean_thickness']),            f(cen['mean_thickness']), 'Mean cortical thickness'),
    ('total_surface',      '皮层总表面积',          'cm2',  f(raw['total_surface_arrea']) / 100.0, f(cen['total_surface_arrea']), 'Total cortical surface area'),
    ('cerebellum_GM',      '小脑灰质体积',          'ml',   f(raw['cerebellum_GM']) / 1000.0,    f(cen['cerebellum_GM']), 'Cerebellar grey-matter volume'),
    ('cerebellum_WM',      '小脑白质体积',          'ml',   f(raw['cerebellum_WM']) / 1000.0,    f(cen['cerebellum_WM']), 'Cerebellar white-matter volume'),
    ('cerebellum_total',   '小脑总体积',            'ml',   f(raw['cerebellum_total']) / 1000.0, f(cen['cerebellum_total']), 'Total cerebellar volume'),
    ('Brain.Stem',         '脑干体积',              'ml',   f(raw['Brain.Stem']) / 1000.0,       f(cen['Brain.Stem']), 'Brainstem volume'),
]

SUBCORT = ['Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala',
           'Accumbens.area', 'VentralDC']
SUBCORT_CN = {'Thalamus': '丘脑', 'Caudate': '尾状核', 'Putamen': '壳核', 'Pallidum': '苍白球',
              'Hippocampus': '海马', 'Amygdala': '杏仁核', 'Accumbens.area': '伏隔核',
              'VentralDC': '腹侧间脑(含丘脑前核/乳头体)'}
SUBCORT_EN = {'Thalamus': 'Thalamus', 'Caudate': 'Caudate', 'Putamen': 'Putamen', 'Pallidum': 'Pallidum',
              'Hippocampus': 'Hippocampus', 'Amygdala': 'Amygdala', 'Accumbens.area': 'Nucleus accumbens',
              'VentralDC': 'Ventral diencephalon (incl. anterior thalamic nuclei/mammillary bodies)'}

def flag(c):
    if c < 0.05:
        return (bi('显著偏低', 'Markedly low'), 'low')
    if c < 0.15:
        return (bi('偏低', 'Low'), 'lowish')
    if c > 0.95:
        return (bi('显著偏高', 'Markedly high'), 'high')
    if c > 0.85:
        return (bi('偏高', 'High'), 'highish')
    return (bi('居中', 'Mid'), 'mid')

def set_centile(new_cen):
    """Populate the centile-derived globals that section 3 renders from, for one
    norm (population-general or subtype-specific). Called before each report build.
    Cortical dicts keep lh_/rh_ columns distinct so all 68 parcels are counted."""
    global cen, GLOBAL_METRICS, thick_c, vol_c, area_c
    global n_th05, n_th10, mean_th, mean_ar, n_ar90, n_vol10, low_th_names, low_vol_names
    cen = new_cen
    GLOBAL_METRICS = _build_global_metrics(cen)
    thick_c = {h: f(cen[h]) for h in cen if h.endswith('_thickness') and h[:3] in ('lh_', 'rh_')}
    vol_c = {h: f(cen[h]) for h in cen if h.endswith('_volume') and h[:3] in ('lh_', 'rh_')}
    area_c = {h: f(cen[h]) for h in cen if h.endswith('_area') and h[:3] in ('lh_', 'rh_')}

    def _disp(kv):
        return (kv[0].replace('_thickness', '').replace('_volume', '').replace('_area', ''), kv[1])

    all_th = sorted(thick_c.items(), key=lambda kv: kv[1])
    n_th05 = sum(1 for v in thick_c.values() if v < 0.05)
    n_th10 = sum(1 for v in thick_c.values() if v < 0.10)
    mean_th = sum(thick_c.values()) / len(thick_c)
    mean_ar = sum(area_c.values()) / len(area_c)
    n_ar90 = sum(1 for v in area_c.values() if v > 0.90)
    n_vol10 = sum(1 for v in vol_c.values() if v < 0.10)
    low_th_names = '、'.join('%s(%.3f)' % _disp(kv) for kv in all_th[:6])
    low_vol = sorted(vol_c.items(), key=lambda kv: kv[1])[:6]
    low_vol_names = '、'.join('%s(%.2f)' % _disp(kv) for kv in low_vol)

# ---------------------------------------------------------------- image assets
FONT_R = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
FONT_B = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

def _rot_label(parts, color):
    tmp = Image.new('RGBA', (4000, 500), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp)
    x = 50.0
    for txt, font, dy in parts:
        d.text((x, 380 + dy), txt, font=font, fill=color, anchor='ls')
        x += d.textlength(txt, font=font)
    a = np.array(tmp)
    ys, xs = np.where(a[:, :, 3] > 8)
    tmp = tmp.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1))
    return tmp.rotate(90, expand=True, resample=Image.BICUBIC)

def _paste(im, tile, x, y):
    im.paste(tile, (int(x), int(y)), tile)

def _text_tile(text, font, color):
    tmp = Image.new('RGBA', (1200, 300), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp)
    d.text((50, 200), text, font=font, fill=color, anchor='ls')
    a = np.array(tmp)
    ys, xs = np.where(a[:, :, 3] > 8)
    return tmp.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1))

def edit_ytitle_mt(im):
    d = ImageDraw.Draw(im)
    d.rectangle([0, 147, 76, 646], fill=(255, 255, 255))
    tile = _rot_label([('mean thickness(mm)', ImageFont.truetype(FONT_R, 46), 0)], (0, 0, 0, 255))
    _paste(im, tile, 44 - tile.width / 2.0, 396 - tile.height / 2.0)
    return im

def edit_ytitle_tsa(im):
    d = ImageDraw.Draw(im)
    d.rectangle([0, 107, 76, 690], fill=(255, 255, 255))
    f = ImageFont.truetype(FONT_R, 46)
    fs = ImageFont.truetype(FONT_R, 29)
    tile = _rot_label([('total surface area(mm', f, 0), ('2', fs, -17), (')', f, 0)], (0, 0, 0, 255))
    _paste(im, tile, 44 - tile.width / 2.0, 398 - tile.height / 2.0)
    return im

def edit_colorbar(im):
    a = np.array(im)
    a[650:687, 1898:2000] = 255
    a[687:705, 1919:2000] = 255
    im = Image.fromarray(a)
    f = ImageFont.truetype(FONT_R, 29)
    t1 = _text_tile('Centile', f, (0, 0, 0, 255))
    _paste(im, t1, 1951 - t1.width / 2.0, 664 - t1.height)
    t2 = _text_tile('1.00', f, (0, 0, 0, 255))
    _paste(im, t2, 1982 - t2.width, 676)
    return im

def edit_global_rows(im):
    a = np.array(im)
    return Image.fromarray(np.concatenate([a[:248], a[310:]]))

def _header_tile(parts, color=(51, 51, 51, 255)):
    tmp = Image.new('RGBA', (900, 120), (255, 255, 255, 0))
    dd = ImageDraw.Draw(tmp)
    x = 30.0
    for txt, font, dy in parts:
        dd.text((x, 80 + dy), txt, font=font, fill=color, anchor='ls')
        x += dd.textlength(txt, font=font)
    a = np.array(tmp)
    ys, xs = np.where(a[:, :, 3] > 8)
    return tmp.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1))

def edit_global_header(im):
    d = ImageDraw.Draw(im)
    d.rectangle([680, 0, 992, 30], fill=(255, 255, 255))
    f = ImageFont.truetype(FONT_B, 14)
    for right, txt in [(779, 'Volume(ml)'), (986, 'Volume(centile)')]:
        tile = _header_tile([(txt, f, 0)])
        _paste(im, tile, right + 1 - tile.width, 19 - tile.height)
    return im

def edit_local_header(im):
    d = ImageDraw.Draw(im)
    d.rectangle([210, 0, 992, 30], fill=(255, 255, 255))
    f = ImageFont.truetype(FONT_B, 14)
    fs = ImageFont.truetype(FONT_B, 9)
    specs = [(297, [('Volume(ml)', f, 0)]),
             (432, [('Volume(centile)', f, 0)]),
             (560, [('Thickness(mm)', f, 0)]),
             (716, [('Thickness(centile)', f, 0)]),
             (872, [('Area(cm', f, 0), ('2', fs, -5), (')', f, 0)]),
             (986, [('Area(centile)', f, 0)])]
    for right, parts in specs:
        tile = _header_tile(parts)
        _paste(im, tile, right + 1 - tile.width, 19 - tile.height)
    return im

EDITS = {'mean_thickness': edit_ytitle_mt,
         'total_surface_area': edit_ytitle_tsa,
         'cortical_volume': edit_colorbar, 'cortical_thickness': edit_colorbar,
         'cortical_area': edit_colorbar, 'subcortical_volume': edit_colorbar,
         'global_summary': lambda im: edit_global_header(edit_global_rows(im)),
         'local_summary': edit_local_header}

def optimize(src_name, max_dim=None, white_bg=True):
    src = os.path.join(BASE, src_name)
    dst = os.path.join(ASSET, src_name.replace('.png', '_opt.png'))
    im = Image.open(src)
    if max_dim and max(im.size) > max_dim:
        im.thumbnail((max_dim, max_dim), Image.LANCZOS)
    if im.mode in ('RGBA', 'LA', 'PA') or 'transparency' in im.info:
        im = im.convert('RGBA')
        white = Image.new('RGBA', im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(white, im)
    im = im.convert('RGB')
    if white_bg:
        arr = np.array(im).astype(np.int16)
        bg = arr[2, 2]
        mask = (np.abs(arr - bg).max(axis=2) <= 12)
        arr[mask] = 255
        im = Image.fromarray(arr.astype(np.uint8))
    stem = src_name.replace('.png', '').replace('_Subtype_specific', '')
    if stem in EDITS:
        im = EDITS[stem](im.convert('RGB'))
    im = im.quantize(colors=256, dither=Image.NONE)
    im.save(dst, optimize=True)
    return dst

CURVES = ['GMV.png', 'sGMV.png', 'WMV.png', 'Ventricles.png', 'mean_thickness.png',
          'total_surface_area.png', 'Cerebellum.png', 'BrainStem.png']
MAPS = ['cortical_volume.png', 'cortical_thickness.png', 'cortical_area.png', 'subcortical_volume.png']
TABLES = ['global_summary.png', 'local_summary.png']

CACHE = {}
def asset_b64(name, max_dim=None, white_bg=True):
    key = (name, max_dim, white_bg)
    if key not in CACHE:
        path = optimize(name, max_dim, white_bg)
        with open(path, 'rb') as fh:
            CACHE[key] = base64.b64encode(fh.read()).decode()
    return CACHE[key]

def img_tag(name, max_dim=None, white_bg=True, alt=''):
    return ('<img src="data:image/png;base64,%s" alt="%s" loading="lazy" '
            'onclick="openLB(this)" />' % (asset_b64(name, max_dim, white_bg), alt))

# ---------------------------------------------------------------- svg lollipop
def lollipop_svg():
    items = sorted(thick_c.items(), key=lambda kv: kv[1])
    n = len(items)
    bw, gap, left, top, H = 12, 5, 46, 18, 240
    W = left + n * (bw + gap) + 20
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" role="img">' % (W, top + H + 70)]
    parts.append('<text x="%d" y="14" class="svg-title" data-lang="zh">68 个皮层分区厚度百分位(升序,受检者个体)</text>' % left)
    parts.append('<text x="%d" y="14" class="svg-title" data-lang="en">Thickness centiles of 68 cortical parcels (ascending, subject)</text>' % left)
    for yv in (0.0, 0.05, 0.5, 1.0):
        y = top + H - yv * H
        parts.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="1" %s/>'
                     % (left - 6, y, W - 8, y, '#c0392b' if yv == 0.05 else '#d8dee6',
                        'stroke-dasharray="4 3"' if yv == 0.05 else ''))
        parts.append('<text x="%d" y="%.1f" text-anchor="end" class="svg-tick">%.2f</text>' % (left - 10, y + 4, yv))
    for i, (name, c) in enumerate(items):
        x = left + i * (bw + gap)
        h = max(2.0, c * H)
        col = '#c0392b' if c < 0.05 else ('#e67e22' if c < 0.15 else '#27ae60')
        parts.append('<rect x="%d" y="%.1f" width="%d" height="%.1f" rx="2" fill="%s"><title>%s: %.4f</title></rect>'
                     % (x, top + H - h, bw, h, col, name, c))
    parts.append('<text x="%d" y="%d" class="svg-note" data-lang="zh">红=&lt;第5百分位 橙=&lt;第15百分位 绿=其余;虚线=第5百分位界;悬停查看分区名</text>'
                 % (left, top + H + 34))
    parts.append('<text x="%d" y="%d" class="svg-note" data-lang="en">red = &lt;5th centile, orange = &lt;15th, green = rest; dashed = 5th centile; hover for parcel name</text>'
                 % (left, top + H + 34))
    parts.append('</svg>')
    return ''.join(parts)

# ---------------------------------------------------------------- html helpers
def bi(cn, en):
    return '<span data-lang="zh">%s</span><span data-lang="en">%s</span>' % (cn, en)

def plain_cn(s):
    s = re.sub(r'<span data-lang="en">.*?</span>', '', s)
    return s.replace('<span data-lang="zh">', '').replace('</span>', '')

INTERP = {}
_interp_path = os.path.join(BASE, 'interpretations.json')
if os.path.exists(_interp_path):
    with open(_interp_path, encoding='utf-8') as fh:
        INTERP = json.load(fh)

def it(key, cn, en):
    idx = None
    if '[' in key:
        key, i = key.split('[')
        idx = int(i.rstrip(']'))
    out = []
    for lang, fb in (('zh', cn), ('en', en)):
        val = INTERP.get(lang, {}).get(key)
        if idx is not None:
            val = val[idx] if isinstance(val, list) and idx < len(val) else None
        out.append(val if isinstance(val, str) and val.strip() else fb)
    return bi(out[0], out[1])

def chip(c):
    label, cls = flag(c)
    return '<span class="chip %s">%s</span>' % (cls, label), cls

def table_global():
    rows = []
    for _k, cn, unit, v, c, en in GLOBAL_METRICS:
        lab, cls = chip(c)
        rows.append('<tr><td>%s(%s)</td><td class="num">%.2f</td><td class="num">%.4f</td>'
                    '<td class="num">%.1f</td><td>%s</td></tr>'
                    % (bi(cn, en), {'ml': 'ml', 'mm': 'mm', 'cm2': 'cm<sup>2</sup>'}[unit], v, c, c * 100, lab))
    return ('<table class="dt"><thead><tr><th>%s</th><th>%s</th><th>%s</th>'
            '<th>%s</th><th>%s</th></tr></thead><tbody>%s</tbody></table>'
            % (bi('指标', 'Metric'), bi('实测值', 'Measured'), bi('百分位', 'Centile'),
               bi('百分位(%)', 'Centile (%)'), bi('判读', 'Reading'), ''.join(rows)))

def table_subcort():
    rows = []
    for side, cn_side, en_side in (('Left', '左', 'Left '), ('Right', '右', 'Right ')):
        for r in SUBCORT:
            key = '%s.%s' % (side, r)
            v = f(raw[key]) / 1000.0
            c = f(cen[key])
            lab, cls = chip(c)
            rows.append('<tr><td>%s</td><td class="num">%.2f</td><td class="num">%.4f</td>'
                        '<td class="num">%.1f</td><td>%s</td></tr>'
                        % (bi(cn_side + SUBCORT_CN[r], en_side + SUBCORT_EN[r]), v, c, c * 100, lab))
    return ('<table class="dt"><thead><tr><th>%s</th><th>%s</th><th>%s</th>'
            '<th>%s</th><th>%s</th></tr></thead><tbody>%s</tbody></table>'
            % (bi('皮层下核团', 'Subcortical nucleus'), bi('体积(ml)', 'Volume (ml)'),
               bi('百分位', 'Centile'), bi('百分位(%)', 'Centile (%)'),
               bi('判读', 'Reading'), ''.join(rows)))

FIG = {}
FIG['GMV.png'] = (bi('全脑灰质体积(GMV)年龄轨迹常模曲线', 'Whole-brain grey-matter volume (GMV) age-trajectory normative curves'),
                  bi('实线为男(浅蓝)/女(红)性别人群中位数轨迹,点线为 95% 参考区间;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                     'Solid lines = male (light blue) / female (red) population median trajectories; dotted lines = 95% reference interval; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['sGMV.png'] = (bi('皮层下灰质体积(sGMV)年龄轨迹常模曲线', 'Subcortical grey-matter volume (sGMV) age-trajectory normative curves'),
                   bi('sGMV 为丘脑、纹状体、海马、杏仁核等皮层下核团灰质合计;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                      'sGMV sums the grey matter of subcortical nuclei (thalamus, striatum, hippocampus, amygdala); dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['WMV.png'] = (bi('白质体积(WMV)年龄轨迹常模曲线', 'White-matter volume (WMV) age-trajectory normative curves'),
                  bi('白质体积反映髓鞘化纤维束总量,青春期后持续上升至中年平台;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                     'White-matter volume reflects the total myelinated fibre bundles, rising after puberty to a mid-life plateau; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['Ventricles.png'] = (bi('脑室系统体积年龄轨迹常模曲线', 'Ventricular system volume age-trajectory normative curves'),
                         bi('脑室体积随龄呈加速增长,与实质体积互为消长;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                            'Ventricular volume grows with accelerating age dependence, in inverse proportion to parenchymal volume; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['mean_thickness.png'] = (bi('平均皮层厚度年龄轨迹常模曲线', 'Mean cortical thickness age-trajectory normative curves'),
                             bi('平均皮层厚度为全脑皮层厚度均值,是皮层完整性最敏感的宏观指标之一;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                                'Mean cortical thickness averages cortical thickness across the brain and is one of the most sensitive macroscopic markers of cortical integrity; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['total_surface_area.png'] = (bi('皮层总表面积年龄轨迹常模曲线', 'Total cortical surface area age-trajectory normative curves'),
                                 bi('皮层总面积反映皮层柱数量的切向扩张,儿童期达峰后轻度下降;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                                    'Total cortical area reflects the tangential expansion of cortical column numbers, peaking in childhood then declining slightly; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['Cerebellum.png'] = (bi('小脑总体积年龄轨迹常模曲线', 'Total cerebellar volume age-trajectory normative curves'),
                         bi('小脑体积含小脑皮层与白质;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                            'Cerebellar volume includes cerebellar cortex and white matter; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
FIG['BrainStem.png'] = (bi('脑干体积年龄轨迹常模曲线', 'Brainstem volume age-trajectory normative curves'),
                        bi('脑干体积反映生命中枢与上/下行传导束的形态保有;暗红点=受检者实测值;绿点=受检者个体伪正常态数值。',
                           'Brainstem volume reflects the preservation of vital centres and ascending/descending tracts; dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'))
for m in MAPS:
    FIG[m] = (bi('皮层/皮层下分区百分位映射图', 'Cortical/subcortical regional centile mapping'),
              bi('inflated 皮层表面四视图(双侧外侧面与内侧面),色标为百分位 0–1:蓝=高百分位、红=低百分位。鼠标悬停于表面视图可按 DKT atlas 显示对应脑区名称(英文,示意性区域定位,仅供读图定向)。',
                 'Inflated cortical surface, four views (bilateral lateral and medial); colour scale = centile 0-1: blue = high centile, red = low centile. Hover over the surface views to see DKT atlas region names (English; schematic localisation, for reading orientation only).'))
FIG['subcortical_volume.png'] = (bi('皮层下分区百分位映射图', 'Subcortical regional centile mapping'),
                                 bi('皮层下核团冠状面断层视图,色标为百分位 0–1:蓝=高百分位、红=低百分位。鼠标悬停可按 FreeSurfer aseg 模板显示对应核团名称(英文,示意性区域定位,仅供读图定向)。',
                                    'Coronal slice through the subcortical nuclei; colour scale = centile 0-1: blue = high centile, red = low centile. Hover to see FreeSurfer aseg structure names (English; schematic localisation, for reading orientation only).'))
FIG['global_summary.png'] = (bi('局部指标汇总表 — 整体与皮下核团体积', 'Regional summary table - whole-brain and subcortical nuclear volumes'),
                             bi('常模输出的整体体积与皮下核团分区体积百分位汇总表;表头单位统一以括号标注(如 Volume(ml)、Volume(centile))。',
                                'Normative-model output table of whole-brain and subcortical regional volumes with centiles; header units are uniformly parenthesised (e.g. Volume(ml), Volume(centile)).'))
FIG['local_summary.png'] = (bi('局部指标汇总表 — 皮层分区', 'Regional summary table - cortical parcels'),
            bi('常模输出的 68 个皮层分区体积/厚度/面积及百分位汇总表;表头单位统一以括号标注;高亮条为高百分位分区。',
               'Normative-model output table of volume/thickness/area and centiles for the 68 cortical parcels; header units are uniformly parenthesised; highlighted rows are high-centile parcels.'))
for _k in ['GMV', 'sGMV', 'WMV', 'Ventricles', 'mean_thickness', 'total_surface_area', 'Cerebellum', 'BrainStem']:
    _t, _c = FIG[_k + '.png']
    FIG[_k + '_Subtype_specific.png'] = (_t + bi('(亚型特异性常模)', ' (subtype-specific normative model)'), _c)
for _m in MAPS:
    if _m == 'subcortical_volume.png':
        FIG[_m.replace('.png', '_Subtype_specific.png')] = (
            bi('皮层下分区百分位映射图', 'Subcortical regional centile mapping'),
            bi('皮层下核团冠状面断层视图,色标为百分位 0–1:蓝=高百分位、红=低百分位。鼠标悬停可按 FreeSurfer aseg 模板显示对应核团名称(英文,示意性区域定位,仅供读图定向)。',
               'Coronal slice through the subcortical nuclei; colour scale = centile 0-1: blue = high centile, red = low centile. Hover to see FreeSurfer aseg structure names (English; schematic localisation, for reading orientation only).'))
    else:
        FIG[_m.replace('.png', '_Subtype_specific.png')] = (
            bi('皮层/皮层下分区百分位映射图', 'Cortical/subcortical regional centile mapping'),
            bi('inflated 皮层表面四视图(双侧外侧面与内侧面),色标为百分位 0–1:蓝=高百分位、红=低百分位。鼠标悬停于表面视图可按 DKT atlas 显示对应脑区名称(英文,示意性区域定位,仅供读图定向)。',
               'Inflated cortical surface, four views (bilateral lateral and medial); colour scale = centile 0-1: blue = high centile, red = low centile. Hover over the surface views to see DKT atlas region names (English; schematic localisation, for reading orientation only).'))
for _t in TABLES:
    if _t == 'global_summary.png':
        FIG[_t.replace('.png', '_Subtype_specific.png')] = (
            bi('局部指标汇总表 — 整体与皮下核团体积', 'Regional summary table - whole-brain and subcortical nuclear volumes'),
            bi('常模输出的整体体积与皮下核团分区体积百分位汇总表;表头单位统一以括号标注(如 Volume(ml)、Volume(centile))。',
               'Normative-model output table of whole-brain and subcortical regional volumes with centiles; header units are uniformly parenthesised (e.g. Volume(ml), Volume(centile)).'))
    else:
        FIG[_t.replace('.png', '_Subtype_specific.png')] = (
            bi('局部指标汇总表 — 皮层分区', 'Regional summary table - cortical parcels'),
            bi('常模输出的 68 个皮层分区体积/厚度/面积及百分位汇总表;表头单位统一以括号标注;高亮条为高百分位分区。',
               'Normative-model output table of volume/thickness/area and centiles for the 68 cortical parcels; header units are uniformly parenthesised; highlighted rows are high-centile parcels.'))

CSS = """
:root{--ink:#1c2733;--mut:#5b6b7c;--bg:#f2f5f9;--card:#ffffff;--line:#e3e9f0;
--acc:#2563eb;--acc2:#0e7490;--low:#c0392b;--lowish:#e67e22;--high:#2980b9;--highish:#5dade2;--mid:#27ae60;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.75 -apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;}
.wrap{max-width:1080px;margin:0 auto;padding:0 20px 80px}
header.hero{background:linear-gradient(135deg,#0f2b46 0%,#155e94 55%,#0e7490 100%);color:#fff;
padding:52px 20px 44px;margin-bottom:34px;border-radius:0 0 22px 22px}
header.hero .wrap{padding-bottom:0}
.hero h1{margin:0 0 10px;font-size:30px;letter-spacing:.5px}
.hero p.sub{margin:0;opacity:.92;font-size:15px}
.hero .meta{margin-top:18px;display:flex;flex-wrap:wrap;gap:8px}
.hero .meta span{background:rgba(255,255,255,.14);border:1px solid rgba(255,255,255,.25);
border-radius:999px;padding:4px 14px;font-size:13px}
nav.toc{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:14px 20px;margin:0 0 26px}
nav.toc b{margin-right:10px}
nav.toc a{color:var(--acc);text-decoration:none;margin-right:16px;font-size:14px;white-space:nowrap}
nav.toc a:hover{text-decoration:underline}
section{background:var(--card);border:1px solid var(--line);border-radius:16px;
padding:26px 30px;margin-bottom:26px;box-shadow:0 1px 3px rgba(16,42,67,.06)}
h2{font-size:21px;margin:0 0 6px;padding-bottom:10px;border-bottom:2px solid var(--line);
display:flex;align-items:center;gap:10px}
h2 .no{background:var(--acc);color:#fff;border-radius:8px;font-size:14px;padding:2px 10px;font-weight:600}
h3{font-size:16.5px;margin:26px 0 8px;color:#123c63}
p{margin:10px 0}
.mut{color:var(--mut);font-size:13.5px}
figure{margin:18px 0;padding:14px;border:1px solid var(--line);border-radius:14px;background:#fbfdff}
figure img{width:100%;border-radius:10px;cursor:zoom-in;display:block;background:#fff;
border:1px solid #eef2f7}
figure.map img{max-width:760px;margin:0 auto}
.mapwrap{position:relative;max-width:760px;margin:0 auto}
.maptip{position:absolute;left:0;top:0;pointer-events:none;display:none;z-index:6;
background:rgba(16,28,40,.92);color:#fff;font:12px/1.4 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
padding:4px 9px;border-radius:7px;white-space:nowrap;transform:translate(-50%,-140%);
box-shadow:0 2px 8px rgba(0,0,0,.25)}
figcaption{margin-top:10px;font-size:13.5px;color:var(--mut)}
figcaption b.t{display:block;color:var(--ink);font-size:14.5px;margin-bottom:2px}
.interp{background:#f4f9ff;border-left:4px solid var(--acc);border-radius:0 10px 10px 0;
padding:12px 16px;margin:10px 0;font-size:14.5px}
.interp.sub{background:#fff7f2;border-left-color:#e67e22}
table.dt{border-collapse:collapse;width:100%;margin:14px 0;font-size:14px;background:#fff}
table.dt th,table.dt td{border:1px solid var(--line);padding:7px 10px;text-align:left}
table.dt th{background:#eef4fb;font-weight:600}
table.dt tr:nth-child(even) td{background:#f8fbfd}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.chip{border-radius:999px;padding:1px 10px;font-size:12.5px;color:#fff;white-space:nowrap}
.chip.low{background:var(--low)}.chip.lowish{background:var(--lowish)}
.chip.high{background:var(--high)}.chip.highish{background:var(--highish)}.chip.mid{background:var(--mid)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}
@media(max-width:860px){.grid2{grid-template-columns:1fr}}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:10px;margin:14px 0}
.kv div{background:#f6f9fc;border:1px solid var(--line);border-radius:10px;padding:10px 14px;font-size:14px}
.kv b{display:block;color:var(--mut);font-size:12.5px;font-weight:600;margin-bottom:2px}
.callout{background:#fff8e6;border:1px solid #f0dca8;border-radius:12px;padding:14px 18px;margin:16px 0;font-size:14.5px}
.svgbox{overflow-x:auto;padding:6px 0}
.svgbox svg{min-width:900px;width:100%}
.svg-title{font-size:15px;font-weight:600;fill:#123c63}
.svg-tick{font-size:11px;fill:#5b6b7c}
.svg-note{font-size:12px;fill:#5b6b7c}
footer{color:var(--mut);font-size:13px;text-align:center;padding:10px 0 30px}
#lb{position:fixed;inset:0;background:rgba(10,20,32,.88);display:none;align-items:center;
justify-content:center;z-index:99;cursor:zoom-out;padding:24px}
#lb img{max-width:96vw;max-height:92vh;border-radius:8px;background:#fff}
ul.tight{margin:8px 0;padding-left:22px}
ul.tight li{margin:5px 0}
html[lang="en"] [data-lang="zh"]{display:none}
html:not([lang="en"]) [data-lang="en"]{display:none}
.langbar{position:fixed;top:14px;right:16px;z-index:120;display:flex;gap:6px}
.langbar button{cursor:pointer;border:1px solid rgba(255,255,255,.55);background:rgba(16,28,40,.55);
color:#fff;border-radius:999px;padding:4px 13px;font-size:12.5px}
.langbar button.on{background:#fff;color:#123c63;font-weight:600}
.langbar a.langsw{display:inline-block;background:rgba(16,28,40,.55);border:1px solid rgba(255,255,255,.55);
color:#fff;border-radius:999px;padding:4px 13px;font-size:12.5px;text-decoration:none}
"""

JS = """
function openLB(el){var lb=document.getElementById('lb');
lb.firstChild.src=el.src;lb.style.display='flex';}
document.addEventListener('click',function(e){if(e.target.id==='lb')e.target.style.display='none';});
document.addEventListener('keydown',function(e){if(e.key==='Escape')document.getElementById('lb').style.display='none';});
var MAPVIEWS=[
{x0:103,y0:431,x1:842,y1:974,t:'lat',ant:1},
{x0:1026,y0:431,x1:1795,y1:974,t:'med',ant:0},
{x0:103,y0:1041,x1:842,y1:1568,t:'lat',ant:0},
{x0:1026,y0:1041,x1:1795,y1:1568,t:'med',ant:1}];
var MAPVN={lat:'lateral view',med:'medial view'};
var DKT={lat:[
['frontalpole',0.05,0.52],['lateralorbitofrontal',0.12,0.70],['parsorbitalis',0.17,0.62],
['rostralmiddlefrontal',0.20,0.45],['superiorfrontal',0.30,0.20],['caudalmiddlefrontal',0.33,0.40],
['parstriangularis',0.30,0.55],['parsopercularis',0.36,0.58],['precentral',0.43,0.30],
['postcentral',0.51,0.28],['supramarginal',0.56,0.42],['angular',0.64,0.45],
['superiorparietal',0.56,0.14],['inferiorparietal',0.65,0.30],['lateraloccipital',0.83,0.38],
['superiortemporal',0.48,0.60],['bankssts',0.55,0.55],['transversetemporal',0.42,0.66],
['middletemporal',0.40,0.76],['inferiortemporal',0.30,0.84],['temporalpole',0.16,0.78],
['fusiform',0.66,0.78],['lingual',0.78,0.72]],
med:[
['frontalpole',0.04,0.50],['medialorbitofrontal',0.11,0.68],['rostralanteriorcingulate',0.30,0.50],
['caudalanteriorcingulate',0.42,0.42],['posteriorcingulate',0.55,0.40],['isthmuscingulate',0.66,0.46],
['superiorfrontal',0.26,0.22],['paracentral',0.47,0.16],['precuneus',0.63,0.22],['cuneus',0.79,0.32],
['pericalcarine',0.87,0.48],['lingual',0.77,0.66],['fusiform',0.63,0.76],['parahippocampal',0.52,0.82],
['entorhinal',0.40,0.80],['temporalpole',0.14,0.76],['middletemporal',0.28,0.84],
['inferiortemporal',0.44,0.86],['lateraloccipital',0.90,0.60]]};
var DKT_LAT_MED=1;
var ASEG=[
['caudate',717,905,32,62],['caudate',1183,905,32,62],
['lateral-ventricle',820,880,120,85],['lateral-ventricle',1080,890,125,95],
['pallidum',637,1180,45,73],['pallidum',1260,1187,48,73],
['putamen',550,1173,55,145],['putamen',1337,1173,62,140],
['thalamus',823,1080,110,115],['thalamus',1070,1087,110,135],
['amygdala',783,1280,130,75],['amygdala',1090,1280,125,75],
['hippocampus',637,1433,110,125],['hippocampus',1237,1433,105,130]];
function mapSetup(img){
  var wrap=img.parentNode,tip=wrap.querySelector('.maptip'),cv=null,ctx=null;
  var hov=img.getAttribute('data-hov');
  function mask(){if(cv)return;cv=document.createElement('canvas');cv.width=400;cv.height=400;
    ctx=cv.getContext('2d',{willReadFrequently:true});ctx.drawImage(img,0,0,400,400);}
  img.addEventListener('mousemove',function(e){
    var r=img.getBoundingClientRect();
    var nx=(e.clientX-r.left)/r.width*2000,ny=(e.clientY-r.top)/r.height*2000;
    var label=null;
    if(hov==='aseg'){
      for(var k=0;k<ASEG.length;k++){var el=ASEG[k];
        var ux=(nx-el[1])/el[3],uy=(ny-el[2])/el[4];
        if(ux*ux+uy*uy<=1){label=el[0]+' · aseg';break;}}
    }else{
      var v=null;
      for(var i=0;i<4;i++){var w=MAPVIEWS[i];if(nx>=w.x0&&nx<=w.x1&&ny>=w.y0&&ny<=w.y1){v=w;break;}}
      if(v){
        try{mask();}catch(err){v=null;}
        if(v){
          var d=ctx.getImageData(Math.round(nx/2000*399),Math.round(ny/2000*399),1,1).data;
          if(d[0]>242&&d[1]>242&&d[2]>242)v=null;
          else if(Math.abs(d[0]-d[1])<10&&Math.abs(d[1]-d[2])<10&&d[0]>105&&d[0]<175)v=null;
        }
      }
      if(v){
        var lx=(nx-v.x0)/(v.x1-v.x0),ly=(ny-v.y0)/(v.y1-v.y0);
        if(v.ant)lx=1-lx;
        var A=DKT[v.t],best=A[0],bd=1e9;
        for(var j=0;j<A.length;j++){var dx=A[j][1]-lx,dy=A[j][2]-ly,dd=dx*dx+dy*dy;if(dd<bd){bd=dd;best=A[j];}}
        label=best[0]+' · '+MAPVN[v.t];
      }
    }
    if(!label){tip.style.display='none';return;}
    tip.textContent=label;
    tip.style.display='block';
    tip.style.left=(e.clientX-r.left)+'px';tip.style.top=(e.clientY-r.top)+'px';});
  img.addEventListener('mouseleave',function(){tip.style.display='none';});
}
function setLang(l){document.documentElement.lang=l;
try{localStorage.setItem('bh_lang',l);}catch(e){}
var bz=document.getElementById('langZh'),be=document.getElementById('langEn');
if(bz)bz.className=(l==='zh'?'on':'');if(be)be.className=(l==='en'?'on':'');}
var _pm=location.pathname.match(/_(zh|en)\.html$/);
var _bar=document.querySelector('.langbar');
if(_pm&&_bar){var _cur=_pm[1];
var _oh=location.pathname.replace(/_(zh|en)\.html$/,'_'+(_cur==='zh'?'en':'zh')+'.html');
_bar.innerHTML='<a class="langsw" href="'+_oh+'">'+(_cur==='zh'?'EN':'中文')+'</a>';}
window.addEventListener('load',function(){
  var l='zh';try{l=localStorage.getItem('bh_lang')||'zh';}catch(e){}
  setLang(l==='en'?'en':'zh');
  var bz=document.getElementById('langZh'),be=document.getElementById('langEn');
  if(bz)bz.onclick=function(){setLang('zh');};
  if(be)be.onclick=function(){setLang('en');};
  document.querySelectorAll('img.mapimg').forEach(function(im){mapSetup(im);});});
"""

def page(title, hero_sub, body):
    return ('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"/>'
            '<meta name="viewport" content="width=device-width,initial-scale=1"/>'
            '<title>%s</title><style>%s</style></head><body>'
            '<div id="lb"><img alt=""/></div>'
            '<div class="langbar"><button id="langZh" type="button">中文</button>'
            '<button id="langEn" type="button">EN</button></div>'
            '%s<div class="wrap"><footer>%s</footer></div>'
            '<script>%s</script></body></html>'
            % (title, CSS, body,
               bi('本报告由结构化 MRI 自由表面形态计量管线自动生成,仅供科学研究使用,不构成任何诊断或临床决策依据。',
                  'This report is auto-generated by a structural-MRI free-surface morphometry pipeline, for scientific research only; '
                  'it does not constitute any diagnosis or clinical decision support.'),
               JS))

def mono(html, lang):
    drop, keep = ('en', 'zh') if lang == 'zh' else ('zh', 'en')
    html = re.sub(r'<span data-lang="%s">.*?</span>' % drop, '', html, flags=re.S)
    html = re.sub(r'<text[^>]*? data-lang="%s">.*?</text>' % drop, '', html, flags=re.S)
    html = re.sub(r'<span data-lang="%s">(.*?)</span>' % keep, r'\1', html, flags=re.S)
    html = re.sub(r'<text([^>]*?) data-lang="%s">(.*?)</text>' % keep,
                  r'<text\1>\2</text>', html, flags=re.S)
    if lang == 'en':
        html = html.replace('<html lang="zh-CN">', '<html lang="en">', 1)
    return html

def fig(name, extra_class='', max_dim=None, white_bg=True, hover=''):
    t, cap = FIG[name]
    tag = img_tag(name, max_dim, white_bg, plain_cn(t).replace('"', '&quot;'))
    if hover:
        tag = ('<div class="mapwrap">%s<div class="maptip"></div></div>'
               % tag.replace('<img ', '<img class="mapimg" data-hov="%s" ' % hover, 1))
    return ('<figure class="%s">%s<figcaption><b class="t">%s</b>%s</figcaption></figure>'
            % (extra_class, tag, t, cap))

def interp(text, sub=False):
    return '<div class="interp%s">%s</div>' % (' sub' if sub else '', text)

_REPORT_DATE = time.strftime('%Y-%m-%d')
SCANNER = '3.0 T · UIH uMR 880'
META = {'id': 'DEMO-S001', 'age': '24', 'sex_zh': '女', 'sex_en': 'female',
        'site_zh': '—', 'site_en': '—', 'moca_zh': '—', 'moca_en': '—',
        'edu_zh': '—', 'edu_en': '—',
        'dx_zh': 'HC(健康对照,研究登记)', 'dx_en': 'HC (healthy control, research registration)'}
_meta_path = os.path.join(BASE, 'subject_meta.json')
if os.path.exists(_meta_path):
    with open(_meta_path, encoding='utf-8') as fh:
        _m = json.load(fh)
    if _m.get('id'):
        META['id'] = _m['id']
    if _m.get('sex'):
        META['sex_zh'] = '女' if _m['sex'] == 'female' else '男'
        META['sex_en'] = _m['sex']
    if _m.get('age'):
        META['age'] = _m['age']
    if _m.get('site'):
        META['site_zh'] = META['site_en'] = _m['site']

SUBJ_BLOCK = """
<section id="subj"><h2><span class="no">1</span>%s</h2>
<div class="kv">
<div><b>%s</b>%s</div>
<div><b>%s</b>%s</div>
<div><b>%s</b>%s</div>
<div><b>%s</b>%s</div>
<div><b>%s</b>%s</div>
<div><b>MoCA / MMSE</b>%s</div>
<div><b>%s</b>%s</div>
<div><b>%s</b>%s</div>
</div>
<p class="mut">%s</p>
</section>
""" % (bi('受检者与数据质量概览', 'Subject and data-quality overview'),
       bi('受检者 ID', 'Subject ID'), META['id'],
       bi('年龄', 'Age'), bi('%s 岁' % META['age'], '%s y' % META['age']),
       bi('性别', 'Sex'), bi(META['sex_zh'], META['sex_en']),
       bi('登记诊断', 'Registered diagnosis'), bi(META['dx_zh'], META['dx_en']),
       bi('扫描设备', 'Scanner'), SCANNER,
       bi(META['moca_zh'], META['moca_en']),
       bi('受教育年限', 'Education (years)'), bi(META['edu_zh'], META['edu_en']),
       bi('采集中心', 'Site'), bi(META['site_zh'], META['site_en']),
       bi('形态计量来自 FreeSurfer 风格自由表面管线(aseg + aparc 分区);百分位为个体实测值在同龄、同性别参考人群分布中的位置(0–1,越低表示越低于参考人群)。',
          'Morphometry comes from a FreeSurfer-style free-surface pipeline (aseg + aparc parcels); centiles place the individual measured value within the age- and sex-matched reference distribution (0-1; lower = further below the reference population).'))


def quant_section(no):
    stats_zh = ('统计:68 个皮层分区中,厚度百分位 &lt; 0.05 者 %d 个、&lt; 0.10 者 %d 个,均值 %.3f;'
                '面积百分位均值 %.3f(&gt;0.90 者 %d 个);体积百分位 &lt; 0.10 者 %d 个。厚度最低分区:%s;体积最低分区:%s。'
                % (n_th05, n_th10, mean_th, mean_ar, n_ar90, n_vol10, low_th_names, low_vol_names))
    stats_en = ('Statistics: of the 68 cortical parcels, %d have thickness centile &lt; 0.05 and %d &lt; 0.10 (mean %.3f); '
                'mean area centile %.3f (%d parcels &gt; 0.90); %d parcels with volume centile &lt; 0.10. '
                'Lowest thickness parcels: %s; lowest volume parcels: %s.'
                % (n_th05, n_th10, mean_th, mean_ar, n_ar90, n_vol10, low_th_names, low_vol_names))
    return """
<section id="quant"><h2><span class="no">%d</span>%s</h2>
<h3>%s</h3>
%s
<h3>%s</h3>
%s
<h3>%s</h3>
<div class="svgbox">%s</div>
<p class="mut">%s</p>
</section>
""" % (no, bi('受检者定量百分位汇总', 'Subject quantitative centile summary'),
       bi('全局指标', 'Global metrics'), table_global(),
       bi('皮层下核团(双侧)', 'Subcortical nuclei (bilateral)'), table_subcort(),
       bi('皮层分区厚度百分位分布(受检者个体)', 'Cortical parcel thickness centile distribution (subject)'),
       lollipop_svg(), bi(stats_zh, stats_en))


# ---------------------------------------------------------------- conclusion
_DKT_CN = {
    'bankssts': '颞上沟岸区', 'caudalanteriorcingulate': '尾状前扣带', 'caudalmiddlefrontal': '尾状中额回',
    'cuneus': '楔叶', 'entorhinal': '内嗅皮层', 'frontalpole': '额极', 'fusiform': '梭状回',
    'inferiorparietal': '下顶叶', 'inferiortemporal': '下颞回', 'insula': '岛叶',
    'isthmuscingulate': '峡部扣带', 'lateraloccipital': '枕外侧回', 'lateralorbitofrontal': '眶额外侧回',
    'lingual': '舌回', 'medialorbitofrontal': '眶额内侧回', 'middletemporal': '中颞回',
    'paracentral': '中央旁小叶', 'parahippocampal': '海马旁回', 'parsopercularis': '额下回盖部',
    'parsorbitalis': '额下回眶部', 'parstriangularis': '额下回三角部', 'pericalcarine': '距状周皮层',
    'postcentral': '中央后回', 'posteriorcingulate': '后扣带', 'precentral': '中央前回',
    'precuneus': '楔前叶', 'rostralanteriorcingulate': '嘴侧前扣带', 'rostralmiddlefrontal': '嘴侧中额回',
    'superiorfrontal': '额上回', 'superiorparietal': '顶上小叶', 'superiortemporal': '颞上回',
    'supramarginal': '缘上回', 'temporalpole': '颞极', 'transversetemporal': '颞横回',
}
_GLOBAL_CN = {'GMV': '全脑灰质', 'sGMV': '皮层下灰质', 'WMV': '白质', 'TCV': '颅内总容积',
              'mean_thickness': '平均皮层厚度', 'total_surface_arrea': '皮层总表面积',
              'cerebellum_total': '小脑总体积', 'Brain.Stem': '脑干', 'Ventricles': '脑室'}
_MET_CN = {'area': '面积', 'thickness': '厚度', 'volume': '体积'}
_CONCL_SKIP_COLS = {'', 'Age', 'Sex'}


def _metric_label(k):
    if k in _GLOBAL_CN:
        return _GLOBAL_CN[k], k
    m = re.match(r'(lh|rh)_([a-z]+)_(area|thickness|volume)$', k)
    if m:
        side, reg, met = m.groups()
        zh = ('左' if side == 'lh' else '右') + _DKT_CN.get(reg, reg) + _MET_CN[met]
        en = ('left ' if side == 'lh' else 'right ') + reg + ' ' + met
        return zh, en
    m = re.match(r'(Left|Right)\.(.+)$', k)
    if m:
        zh = ('左' if m.group(1) == 'Left' else '右') + SUBCORT_CN.get(m.group(2), m.group(2))
        en = ('left ' if m.group(1) == 'Left' else 'right ') + SUBCORT_EN.get(m.group(2), m.group(2))
        return zh, en
    return k, k


def _risk_grade(n):
    if n >= 10:
        return '高风险', 'high risk'
    if n >= 5:
        return '中风险', 'medium risk'
    if n >= 2:
        return '中低风险', 'medium-low risk'
    if n == 1:
        return '低风险', 'low risk'
    return '较低神经系统疾病风险', 'lower neurological-disease risk'


def _norm_metric_keys(primary):
    """Metric columns of the norm-specific centile export (the modelled set:
    228 regional/global metrics + Ventricles), excluding id/age/sex columns."""
    path = os.path.join(BASE, primary)
    if os.path.exists(path):
        with open(path, newline='', encoding='utf-8') as fh:
            hdr = next(csv.reader(fh))
        return {h for h in hdr
                if h not in _CONCL_SKIP_COLS and not h.startswith('Freesufer_Path')}
    return {h for h in cen_leg
            if h not in _CONCL_SKIP_COLS and not h.startswith('Freesufer_Path')}


def abnormal_indicators(cen_dict, keys=None):
    """Ventricles abnormal when centile > 0.975; every other modelled metric
    abnormal when centile < 0.025. Returns [(metric, centile)] sorted by severity."""
    ks = keys if keys is not None else set(cen_dict)
    out = []
    for k in ks:
        try:
            v = float(cen_dict.get(k))
        except (TypeError, ValueError):
            continue
        if (k == 'Ventricles' and v > 0.975) or (k != 'Ventricles' and v < 0.025):
            out.append((k, v))
    out.sort(key=lambda kv: (kv[0] != 'Ventricles', -kv[1] if kv[0] == 'Ventricles' else kv[1]))
    return out


def conclusion_section(cen_dict, keys=None):
    abn = abnormal_indicators(cen_dict, keys)
    n_total = len(keys) if keys is not None else len(cen_dict)
    grade_zh, grade_en = _risk_grade(len(abn))
    parts_zh, parts_en = [], []
    for k, v in abn[:6]:
        zh, en = _metric_label(k)
        parts_zh.append('%s%.3f' % (zh, v))
        parts_en.append('%s %.3f' % (en, v))
    names_zh = '、'.join(parts_zh)
    names_en = ', '.join(parts_en)
    if len(abn) > 6:
        names_zh += '等(共%d项)' % len(abn)
        names_en += ', etc. (%d in total)' % len(abn)
    if abn:
        sum_zh = '异常指标%d项(%s);神经系统疾病风险评分:%s。' % (len(abn), names_zh, grade_zh)
        sum_en = '%d abnormal indicators (%s); neurological-disease risk: %s.' % (len(abn), names_en, grade_en)
    else:
        sum_zh = '异常指标0项;神经系统疾病风险评分:%s。' % grade_zh
        sum_en = 'No abnormal indicators (0); neurological-disease risk: %s.' % grade_en
    rule = bi('判定规则:全部%d项形态计量指标中,脑室以百分位&gt;0.975判为异常,其余%d项以百分位&lt;0.025判为异常;'
              '神经系统疾病风险分级:0项=较低神经系统疾病风险、1项=低风险、2–4项=中低风险、5–9项=中风险、≥10项=高风险。'
              '本结论为研究性概述,不构成诊断或临床决策依据。'
              % (n_total, n_total - 1),
              'Rule: among the %d morphometric metrics, ventricles are abnormal when centile &gt; 0.975 and the other %d '
              'metrics when centile &lt; 0.025; neurological-disease risk grades: 0 = lower, 1 = low, 2-4 = medium-low, '
              '5-9 = medium, &gt;=10 = high. Research summary only; not a diagnosis or clinical decision support.'
              % (n_total, n_total - 1))
    return ('<section id="conclusion"><h2><span class="no">7</span>%s</h2>' % bi('结论总结', 'Conclusions summary')
            + '<div class="callout">%s</div>' % it('concl', sum_zh, sum_en)
            + '<p class="mut">%s</p>' % rule
            + '</section>')


# ================================================================ report 1
_pop_keys = _norm_metric_keys('individual_centile_score_overall_population.csv')
_sub_keys = _norm_metric_keys('individual_centile_score_subtype_specific.csv')


def build_population():
    hero = ('<header class="hero"><div class="wrap"><h1>%s</h1>'
            '<p class="sub">%s</p>'
            '<div class="meta"><span>%s</span><span>%s</span>'
            '<span>%s</span><span>%s</span>'
            '<span>%s</span></div></div></header>'
            % (bi('个体化脑健康结构影像报告(I)', 'Individualized Brain-Health Structural Imaging Report (I)'),
               bi('全人群通用常模参照 · 年龄-性别特异性轨迹对照 · 百分位定量解读',
                  'Population-general normative reference · age- and sex-specific trajectory comparison · quantitative centile interpretation'),
               bi('受检者:' + META['id'], 'Subject: ' + META['id']),
               bi('%s 岁 / %s' % (META['age'], META['sex_zh']),
                  '%s y / %s' % (META['age'], META['sex_en'])),
               SCANNER,
               bi('报告日期:' + _REPORT_DATE, 'Report date: ' + _REPORT_DATE),
               bi('常模类型:全人群通用常模', 'Norm type: population-general')))
    toc = ('<div class="wrap"><nav class="toc"><b>%s</b>'
           '<a href="#legend">%s</a><a href="#subj">%s</a>'
           '<a href="#global">%s</a>'
           '<a href="#quant">%s</a><a href="#maps">%s</a>'
           '<a href="#tables">%s</a><a href="#discuss">%s</a>'
           '<a href="#conclusion">%s</a></nav></div>'
           % (bi('目录', 'Contents'), bi('阅读说明', 'Reading guide'), bi('受检者概览', 'Subject overview'),
              bi('全局指标曲线', 'Global metric curves'), bi('定量百分位汇总', 'Quantitative centile summary'),
              bi('分区百分位映射图', 'Regional centile maps'), bi('汇总表', 'Summary tables'),
              bi('综合解读', 'Integrated interpretation'), bi('结论总结', 'Conclusions summary')))

    s0 = """
<section id="legend"><h2><span class="no">0</span>%s</h2>
<ul class="tight">
<li>%s</li>
<li>%s</li>
<li>%s</li>
<li>%s</li>
<li>%s</li>
</ul>
</section>""" % (
        bi('阅读说明与图例', 'Reading guide and legend'),
        bi('<b>曲线:</b>浅蓝 = 男性参考人群,红 = 女性参考人群;实线 = 中位数轨迹,点线 = 95% 参考区间(第 2.5–97.5 百分位)。',
           '<b>Curves:</b> light blue = male reference population, red = female; solid = median trajectory, dotted = 95% reference interval (2.5th-97.5th centile).'),
        bi('<b>散点:</b>暗红点 = 受检者实测值;绿点 = 受检者个体伪正常态数值。',
           '<b>Scatter:</b> dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'),
        bi('<b>百分位:</b>0–1 取值,表示个体实测值在同龄同性别参考人群中的累积位置;&lt;0.05 判为显著偏低,&gt;0.95 判为显著偏高。',
           '<b>Centile:</b> 0-1 value, the cumulative position of the individual measured value within the age- and sex-matched reference population; &lt;0.05 = markedly low, &gt;0.95 = markedly high.'),
        bi('<b>映射图色标:</b>蓝 = 高百分位,红 = 低百分位;皮层表面视图悬停按 DKT atlas 显示脑区名称,皮层下冠状面映射图悬停按 FreeSurfer aseg 模板显示核团名称(均为英文)。',
           '<b>Map colour scale:</b> blue = high centile, red = low centile; hovering cortical surface views shows DKT atlas region names and hovering the subcortical coronal map shows FreeSurfer aseg structure names (both in English).'),
        bi('<b>单位标注:</b>全文单位一律以括号形式标注于指标名后(如 Volume(ml)、Thickness(mm)、Surface.area(cm²))。',
           '<b>Unit notation:</b> units are uniformly parenthesised after the metric name (e.g. Volume(ml), Thickness(mm), Surface.area(cm2)).'))


    if not INTERP:
        s0 += '<div class="callout">%s</div>' % bi(
            '注意:解读引擎(Qoder CLI 或已配置大模型)不可用,以下解读为冻结演示模板文本,不适用于其他受检者。',
            'Note: no interpretation engine (Qoder CLI or a configured model) was available; the interpretations below are frozen demo template text and must not be applied to other subjects.')
    g = ['<section id="global"><h2><span class="no">2</span>%s</h2>' % bi('全局脑结构指标(含小脑与脑干):常模曲线对照', 'Global brain structural metrics (incl. cerebellum and brainstem): normative curve comparison'),
         '<p class="mut">%s</p>' % bi('以下八幅曲线图将受检者实测值(暗红点)置于同龄、同性别人群轨迹中对照(含小脑与脑干),并标出受检者个体伪正常态数值(绿点);每图附医学与神经科学解读。',
                                     'The eight curve plots below place the subject measured values (dark-red dot) within age- and sex-matched population trajectories (including cerebellum and brainstem), and mark the subject individual pseudo-normal value (green dot); each figure carries a medical and neuroscience interpretation.')]
    g.append(fig('GMV.png') + interp(it('GMV',
        '受检者 GMV 为 374.7 ml,位于同龄女性第 9.5 百分位,低于女性中位数轨迹并接近 95% 区间下界。'
        '灰质体积是神经元胞体、树突、突触与胶质含量的全局代理指标,与计算容量和认知储备相关。'
        '值得注意的是其颅内总容积 TCV 处于第 61.5 百分位(全脑尺寸正常),故灰质偏低并非全脑等比例缩小的结果,'
        '而是灰质成分相对于个体自身脑容量的相对不足——这种"脑容量正常而灰质占比低"的模式更支持成分特异性(而非全局性)偏差。',
        'The subject\'s GMV is 374.7 ml, at the 9.5th centile of age-matched females, below the female median trajectory and close to the lower bound of the 95% interval. '
        'Grey-matter volume is a global proxy for neuronal somata, dendrites, synapses and glial content, and relates to computational capacity and cognitive reserve. '
        'Notably, intracranial total volume (TCV) sits at the 61.5th centile (normal whole-brain size), so the grey-matter deficit is not a proportional whole-brain shrinkage '
        'but a relative shortage of grey matter with respect to the individual\'s own brain volume - a "normal brain volume with low grey-matter fraction" pattern that favours component-specific (rather than global) deviation.')))
    g.append(fig('sGMV.png') + interp(it('sGMV',
        '皮层下灰质体积 49.6 ml,第 14.5 百分位,轻度偏低,散点位于女性中位数轨迹下方。'
        '皮层下核团(丘脑-纹状体-海马-杏仁核)是皮层环路的枢纽与中继,其灰质合计偏低与后文海马、丘脑的单核团低百分位相互印证,'
        '提示灰质不足同时累及皮层与皮层下两级结构,符合弥散性(而非单灶性)受累模式。',
        'Subcortical grey-matter volume 49.6 ml, 14.5th centile, mildly low, with the scatter point below the female median trajectory. '
        'Subcortical nuclei (thalamus-striatum-hippocampus-amygdala) are hubs and relays of cortical circuits; this combined grey-matter deficit corroborates the low single-nucleus centiles of hippocampus and thalamus below, '
        'indicating that grey-matter insufficiency involves both cortical and subcortical levels, consistent with a diffuse (not single-focus) pattern.')))
    g.append(fig('WMV.png') + interp(it('WMV',
        '白质体积 433.7 ml,第 54.1 百分位,恰居同龄中位数水平。白质体积主要反映髓鞘化轴突总量,'
        '在 20–40 岁仍处平台期。"灰质低、白质正常"的解离提示:髓鞘发育与纤维束构建过程相对保有,'
        '偏差集中于灰质(神经元胞体/突触)成分;该解离本身即是有鉴别价值的表型特征,'
        '可与以白质脱失为主的疾病模式(如脱髓鞘病)相区分。',
        'White-matter volume 433.7 ml, 54.1st centile, exactly at the age-matched median. White-matter volume mainly reflects total myelinated axons '
        'and remains on a plateau between 20-40 years. The "low grey matter, normal white matter" dissociation indicates relatively preserved myelination and fibre-tract construction, '
        'with deviation confined to the grey-matter (neuronal soma/synaptic) component; this dissociation is itself a discriminating phenotypic feature, '
        'separating it from disease patterns dominated by white-matter loss (e.g. demyelinating disease).')))
    g.append(fig('Ventricles.png') + interp(it('Ventricles',
        '脑室系统体积 27.0 ml,第 88.4 百分位,偏高。颅腔为封闭容积,脑室与实质互为消长(ex vacuo 机制);'
        '在年轻受检者中,该水平更宜解读为实质(尤其灰质)体积偏低的代偿性镜像,而非原发性脑积水;'
        '但脑室百分位是纵向监测最敏感的指标之一,建议作为随访锚点:若随时间上移跨越区间,则提示进行性实质丢失。',
        'Ventricular system volume 27.0 ml, 88.4th centile, high. The cranial cavity is a closed volume: ventricles and parenchyma trade off (ex vacuo mechanism). '
        'In a young subject this level is best read as a compensatory mirror of low parenchymal (especially grey-matter) volume rather than primary hydrocephalus; '
        'ventricular centile is however one of the most sensitive longitudinal monitoring metrics and is recommended as a follow-up anchor: upward crossing of the interval over time would indicate progressive parenchymal loss.')))
    g.append(fig('mean_thickness.png') + interp(it('mean_thickness',
        '平均皮层厚度 2.156 mm,第 0.03 百分位,显著低于同龄女性参考(中位数约 2.4–2.5 mm),为本报告最突出的单指标偏差。'
        '皮层厚度反映神经元分层、树突棘密度、胶质与微血管含量,是皮层完整性最敏感的宏观指标。'
        '24 岁处于青春期修剪完成后的稳定期,如此低值在神经退行性框架下不典型(退行性病程多于中年后起病);'
        '应优先系统排除技术因素(序列/SNR、厚度算法与版本、表面重建质量),若在跨管线、跨时间点可复现,'
        '则指向弥散性皮层发育性变薄,并与认知筛查所提示的偏低表现构成可检验的临床-影像关联。',
        'Mean cortical thickness 2.156 mm, 0.03rd centile, markedly below age-matched female reference (median about 2.4-2.5 mm) - the most prominent single-metric deviation in this report. '
        'Cortical thickness reflects neuronal lamination, dendritic spine density, glial and microvascular content, and is the most sensitive macroscopic marker of cortical integrity. '
        'At 24 years, after completion of adolescent pruning, such a low value is atypical in a neurodegenerative frame (degenerative courses usually begin after midlife); '
        'technical factors (sequence/SNR, thickness algorithm and version, surface reconstruction quality) should be excluded first; if reproducible across pipelines and time points, '
        'it points to diffuse developmental cortical thinning and forms a testable clinical-imaging association with the low cognitive-screening performance.')))
    g.append(fig('total_surface_area.png') + interp(it('total_surface_area',
        '皮层总表面积 1557.3 cm<sup>2</sup>,第 62.7 百分位,居中偏高。按放射单元假说,皮层面积取决于皮层柱数量的切向扩张,'
        '与厚度(柱内放射结构、修剪与髓鞘平衡)在发育与遗传上相对独立。'
        '"面积正常-厚度显著低"的解离模式将偏差定位到放射(纵向)维度,而切向扩张与折叠生成相对保有;'
        '该模式在神经发育障碍鉴别中具有指向意义,也提示不应将本例简单归类为"弥漫性脑萎缩"。',
        'Total cortical surface area 1557.3 cm<sup>2</sup>, 62.7th centile, mid-to-high. Under the radial unit hypothesis, cortical area depends on tangential expansion of cortical column numbers, '
        'developmentally and genetically relatively independent of thickness (intra-column radial structure, pruning and myelination balance). '
        'The "normal area - markedly low thickness" dissociation localises the deviation to the radial (longitudinal) dimension while tangential expansion and folding generation are relatively preserved; '
        'this pattern is informative in the differential diagnosis of neurodevelopmental disorders and argues against simply classifying this case as "diffuse cerebral atrophy".')))
    g.append(fig('Cerebellum.png') + interp(it('Cerebellum',
        '小脑总体积 109.7 ml(第 8.1 百分位),其中小脑灰质 86.8 ml(10.8)、小脑白质 22.9 ml(9.1),整体偏低,'
        '暗红点贴近女性 95% 区间下界。小脑不仅是运动协调与平衡中枢,还经小脑-丘脑-皮层环路参与工作记忆、语言、'
        '执行与情绪调节("小脑认知情感综合征"即源于此)。小脑与大脑皮层在发育时间窗与遗传程序上耦合,'
        '小脑偏低与弥散性皮层变薄共现,提示二者可能共享同一发育性背景,而非独立事件。',
        'Total cerebellar volume 109.7 ml (8.1st centile), including cerebellar grey matter 86.8 ml (10.8) and white matter 22.9 ml (9.1), overall low, '
        'with the dark-red point near the lower bound of the female 95% interval. Beyond motor coordination and balance, the cerebellum participates in working memory, language, '
        'executive and emotional regulation via the cerebello-thalamo-cortical loop ("cerebellar cognitive affective syndrome" arises from this). '
        'Cerebellum and cerebral cortex are coupled in developmental timing and genetic programmes; the co-occurrence of low cerebellar volume with diffuse cortical thinning suggests a shared developmental background rather than independent events.')))
    g.append(fig('BrainStem.png') + interp(it('BrainStem',
        '脑干体积 19.4 ml,第 55.5 百分位,居中。脑干承载呼吸、心血管中枢与全部上/下行传导束,'
        '其形态保有提示生命功能相关结构未受累,也为"偏差选择性累及灰质成分"的整体模式补充了一枚阴性证据。',
        'Brainstem volume 19.4 ml, 55.5th centile, mid-range. The brainstem carries respiratory and cardiovascular centres and all ascending/descending tracts; '
        'its preserved morphology indicates sparing of life-function structures and adds a piece of negative evidence to the overall pattern of "deviation selectively involving grey-matter components".')))
    g.append('</section>')

    q = quant_section(3)

    m = ['<section id="maps"><h2><span class="no">4</span>%s</h2>' % bi('皮层/皮层下分区百分位映射图', 'Cortical/subcortical regional centile mapping'),
         '<p class="mut">%s</p>' % bi('以下映射图为<b>全人群通用常模</b>下的分区百分位映射,与第 3 节定量汇总同一常模口径;'
         '鼠标悬停时,皮层表面四视图按 DKT atlas 显示脑区名称(英文),皮层下冠状面映射图按 FreeSurfer aseg 模板显示核团名称(英文)。',
         'The maps below show regional centile mapping under the <b>population-general normative model</b>, the same calibre as the section 3 quantitative summary; '
         'on hover, the four cortical surface views show DKT atlas region names (English) and the subcortical coronal map shows FreeSurfer aseg structure names (English).')]
    for name in MAPS:
        m.append(fig(name, extra_class='map', max_dim=2000,
                     hover='aseg' if name.startswith('subcortical') else 'dkt'))
    m.append('</section>')

    t = ['<section id="tables"><h2><span class="no">5</span>%s</h2>' % bi('局部指标汇总表', 'Regional summary tables'),
         '<p class="mut">%s</p>' % bi('常模输出的局部指标汇总表:第一张给出 TCV、GMV、sGMV、WMV、脑室、小脑、脑干等整体体积'
         '与皮下核团分区体积及百分位;第二张给出 68 个皮层分区的体积/厚度/面积及百分位。'
         '两表表头单位统一以括号标注(Volume(ml)、Volume(centile)、Thickness(mm)、Thickness(centile)、Area(cm²)、Area(centile))。',
         'Normative-model regional summary tables: the first lists whole-brain volumes (TCV, GMV, sGMV, WMV, ventricles, cerebellum, brainstem) '
         'and subcortical nuclear volumes with centiles; the second lists volume/thickness/area and centiles for the 68 cortical parcels. '
         'Header units are uniformly parenthesised (Volume(ml), Volume(centile), Thickness(mm), Thickness(centile), Area(cm2), Area(centile)).')]
    t.append(fig('global_summary.png', white_bg=False))
    t.append(fig('local_summary.png', extra_class='map', white_bg=False))
    t.append('</section>')

    p61_zh = ('受检者(24 岁,女)的结构表型可归纳为:<b>弥散性皮层变薄</b>(68 区中 %d 区厚度百分位 &lt; 0.05、%d 区 &lt; 0.10,均值 %.3f),'
              '<b>全脑与皮层下灰质体积偏低</b>(GMV 9.5、sGMV 14.5 百分位),<b>双侧海马显著偏低</b>(左 0.55、右 0.01 百分位),'
              '<b>小脑偏低</b>(8.1 百分位),<b>脑室相对增大</b>(88.4 百分位);同时<b>皮层表面积(62.7)、白质体积(54.1)、脑干(55.5)、背侧纹状体(尾状核 76.1/89.4、壳核 64.6/55.4)保有</b>。'
              '该组合构成一个成分选择性、弥散且大致对称的偏差模式。') % (n_th05, n_th10, mean_th)
    p61_en = ('The structural phenotype of the subject (24 y, female) can be summarised as: <b>diffuse cortical thinning</b> (%d of 68 parcels with thickness centile &lt; 0.05, %d &lt; 0.10, mean %.3f), '
              '<b>low whole-brain and subcortical grey-matter volume</b> (GMV 9.5, sGMV 14.5 centile), <b>markedly low bilateral hippocampus</b> (L 0.55, R 0.01 centile), '
              '<b>low cerebellum</b> (8.1 centile), <b>relatively enlarged ventricles</b> (88.4 centile); with <b>preserved cortical surface area (62.7), white-matter volume (54.1), brainstem (55.5) and dorsal striatum (caudate 76.1/89.4, putamen 64.6/55.4)</b>. '
              'This combination forms a component-selective, diffuse and roughly symmetric deviation pattern.') % (n_th05, n_th10, mean_th)
    discuss = ('<section id="discuss"><h2><span class="no">6</span>%s</h2>' % bi('医学与神经科学意义综合解读', 'Medical and neuroscience integrated interpretation')
               + '<h3>%s</h3>' % bi('6.1 表型归纳', '6.1 Phenotype summary')
               + '<p>%s</p>' % it('discuss61', p61_zh, p61_en)
               + '<h3>%s</h3>' % bi('6.2 神经发育视角', '6.2 Neurodevelopmental perspective')
               + '<p>%s</p>' % it('discuss62',
                   '厚度-面积解离符合放射单元假说(radial unit hypothesis):切向扩张(柱数量、折叠)相对保有,而柱内放射结构或修剪/髓鞘平衡异常。'
                   '海马、小脑与丘脑的同时受累提示影响多系统的发育程序(产前-围产期神经元生成、迁移或修剪节律)可能早在围产期前后即已启动;'
                   '弥散对称的模式不支持血管性、肿瘤性或炎症后等局灶病因。年轻年龄、白质保有与缺乏区域选择性萎缩梯度,'
                   '使典型神经退行性病程(AD 型内嗅-内侧颞叶梯度、FTD 型额颞梯度)的可能性较低,'
                   '但海马的极低百分位要求以基线-随访对比排除极早发型退行过程。',
                   'The thickness-area dissociation fits the radial unit hypothesis: tangential expansion (column number, folding) relatively preserved, while intra-column radial structure or pruning/myelination balance is abnormal. '
                   'Simultaneous involvement of hippocampus, cerebellum and thalamus suggests that a multi-system developmental programme (prenatal-perinatal neurogenesis, migration or pruning rhythm) may have started around the perinatal period; '
                   'the diffuse symmetric pattern argues against vascular, neoplastic or post-inflammatory focal causes. Young age, preserved white matter and absence of a regionally selective atrophy gradient '
                   'make typical neurodegenerative courses (AD-type entorhinal-medial temporal gradient, FTD-type frontotemporal gradient) unlikely, '
                   'but the extremely low hippocampal centile requires baseline-follow-up comparison to exclude very-early-onset degeneration.')
               + '<h3>%s</h3>' % bi('6.3 技术与测量视角(优先排除)', '6.3 Technical and measurement perspective (exclude first)')
               + '<p>%s</p>' % it('discuss63',
                   '皮层厚度估计对序列、SNR、场强与算法版本高度敏感;Euler 数(-16/-24)提示表面重建存在轻度拓扑缺陷。'
                   '在给出任何生物学结论前,应以第二管线(如 CAT12/其他分割策略)与重复扫描交叉验证;'
                   '海马-杏仁核边界是自动分割的已知误差源,极低海马百分位需人工核查或多算法共识确认。',
                   'Cortical thickness estimation is highly sensitive to sequence, SNR, field strength and algorithm version; Euler numbers (-16/-24) indicate mild topological defects in surface reconstruction. '
                   'Before any biological conclusion, cross-validation with a second pipeline (e.g. CAT12/other segmentation strategy) and repeat scanning is required; '
                   'the hippocampus-amygdala boundary is a known automatic-segmentation error source, and the extremely low hippocampal centile needs manual check or multi-algorithm consensus.')
               + '<h3>%s</h3>' % bi('6.4 临床关联与建议', '6.4 Clinical associations and recommendations')
               + '<ul class="tight">'
               + '<li>%s</li>' % it('discuss64[0]', '认知筛查得分偏低(见第 1 节)与海马-皮层低体积构成可检验的关联;建议按记忆、执行、语言、视空间分域的成套神经心理评估。',
                                   'Low cognitive-screening performance (see section 1) and the low hippocampal-cortical volumes form a testable association; a domain-wise neuropsychological battery (memory, executive, language, visuospatial) is recommended.')
               + '<li>%s</li>' % it('discuss64[1]', '以脑室百分位与平均厚度为纵向锚点,12–24 个月随访复扫:稳定则支持发育性/静止性表型,上移/下移则提示进行性过程。',
                                   'Use ventricular centile and mean thickness as longitudinal anchors with repeat MRI at 12-24 months: stability supports a developmental/static phenotype, upward/downward drift indicates a progressive process.')
               + '<li>%s</li></ul>' % it('discuss64[2]', '登记诊断为 HC(健康对照),但影像与认知筛查均提示偏离;研究层面应保留其作为"影像-认知偏离"候选个案的随访价值,临床层面建议神经内科门诊评估。',
                                         'The registered diagnosis is HC (healthy control), yet both imaging and cognitive screening indicate deviation; at research level the case retains follow-up value as an "imaging-cognitive deviation" candidate, and at clinical level a neurology outpatient evaluation is suggested.')
               + '<div class="callout">%s</div>' % bi('本报告全部结论基于单次横断面结构 MRI 的形态计量,仅供科学研究使用;不构成诊断,不替代临床面诊与纵向评估。',
                                                      'All conclusions of this report are based on a single cross-sectional structural MRI morphometry, for scientific research only; they do not constitute a diagnosis and do not replace clinical examination or longitudinal assessment.')
               + '</section>')

    body = toc + '<div class="wrap">' + s0 + SUBJ_BLOCK + ''.join(g) + q + ''.join(m) + ''.join(t) + discuss + conclusion_section(cen_pop, _pop_keys) + '</div>'
    return page('个体化脑健康结构影像报告(I)— 全人群通用常模参照', '', hero + body)


# ================================================================ report 2
def build_subtype():
    hero = ('<header class="hero"><div class="wrap"><h1>%s</h1>'
            '<p class="sub">%s</p>'
            '<div class="meta"><span>%s</span><span>%s</span>'
            '<span>%s</span><span>%s</span>'
            '<span>%s</span></div></div></header>'
            % (bi('个体化脑健康结构影像报告(II)', 'Individualized Brain-Health Structural Imaging Report (II)'),
               bi('亚型特异性常模参照 · subtype-specific normative trajectories · 个体特异性偏差识别',
                  'Subtype-specific normative reference · subtype-specific normative trajectories · individual-specific deviation detection'),
               bi('受检者:' + META['id'], 'Subject: ' + META['id']),
               bi('%s 岁 / %s' % (META['age'], META['sex_zh']),
                  '%s y / %s' % (META['age'], META['sex_en'])),
               SCANNER,
               bi('报告日期:' + _REPORT_DATE, 'Report date: ' + _REPORT_DATE),
               bi('常模类型:亚型特异性常模', 'Norm type: subtype-specific')))
    toc = ('<div class="wrap"><nav class="toc"><b>%s</b>'
           '<a href="#legend">%s</a><a href="#subj">%s</a>'
           '<a href="#global">%s</a><a href="#quant">%s</a>'
           '<a href="#maps">%s</a><a href="#tables">%s</a>'
           '<a href="#discuss">%s</a><a href="#conclusion">%s</a></nav></div>'
           % (bi('目录', 'Contents'), bi('阅读说明', 'Reading guide'), bi('受检者概览', 'Subject overview'),
              bi('全局指标曲线', 'Global metric curves'), bi('定量百分位汇总', 'Quantitative centile summary'),
              bi('分区百分位映射图', 'Regional centile maps'), bi('汇总表', 'Summary tables'),
              bi('综合解读', 'Integrated interpretation'), bi('结论总结', 'Conclusions summary')))

    s0 = """
<section id="legend"><h2><span class="no">0</span>%s</h2>
<ul class="tight">
<li>%s</li>
<li>%s</li>
<li>%s</li>
<li>%s</li>
<li>%s</li>
</ul>
</section>""" % (
        bi('阅读说明与图例', 'Reading guide and legend'),
        bi('<b>曲线:</b>浅蓝 = 男性、红 = 女性<b>亚型特异性</b>参考轨迹;实线 = 中位数,点线 = 95% 参考区间。',
           '<b>Curves:</b> light blue = male, red = female <b>subtype-specific</b> reference trajectories; solid = median, dotted = 95% reference interval.'),
        bi('<b>散点:</b>暗红点 = 受检者实测值;绿点 = 受检者个体伪正常态数值。',
           '<b>Scatter:</b> dark-red dot = subject measured value; green dot = subject individual pseudo-normal value.'),
        bi('<b>与报告(I)的关系:</b>本报告全部曲线与映射图基于<b>亚型特异性常模</b>;受检者对<b>全人群通用常模</b>的百分位定量见报告(I),可与本报告并列判读。第 3 节定量表所采用的常模口径以该节说明为准。',
           '<b>Relation to report (I):</b> all curves and maps in this report use the <b>subtype-specific norm</b>; the subject\'s quantification against the <b>population-general norm</b> is in report (I), for side-by-side reading. The norm calibre used by the section 3 tables is stated in that section.'),
        bi('<b>映射图色标:</b>蓝 = 高百分位,红 = 低百分位;皮层表面视图悬停按 DKT atlas 显示脑区名称,皮层下冠状面映射图悬停按 FreeSurfer aseg 模板显示核团名称(均为英文)。',
           '<b>Map colour scale:</b> blue = high centile, red = low centile; hovering cortical surface views shows DKT atlas region names and hovering the subcortical coronal map shows FreeSurfer aseg structure names (both in English).'),
        bi('<b>单位标注:</b>全文单位一律以括号形式标注于指标名后(如 Volume(ml)、Thickness(mm)、Surface.area(cm²))。',
           '<b>Unit notation:</b> units are uniformly parenthesised after the metric name (e.g. Volume(ml), Thickness(mm), Surface.area(cm2)).'))


    if not INTERP:
        s0 += '<div class="callout">%s</div>' % bi(
            '注意:解读引擎(Qoder CLI 或已配置大模型)不可用,以下解读为冻结演示模板文本,不适用于其他受检者。',
            'Note: no interpretation engine (Qoder CLI or a configured model) was available; the interpretations below are frozen demo template text and must not be applied to other subjects.')
    g = ['<section id="global"><h2><span class="no">2</span>%s</h2>' % bi('全局脑结构指标(含小脑与脑干):亚型特异性常模对照', 'Global brain structural metrics (incl. cerebellum and brainstem): subtype-specific normative comparison'),
         '<p class="mut">%s</p>' % bi('曲线为亚型特异性参考轨迹(含小脑与脑干);暗红点 = 受检者实测值,绿点 = 受检者个体伪正常态数值。解读侧重"亚型参考是否改变结论"。',
                                      'Curves are subtype-specific reference trajectories (including cerebellum and brainstem); dark-red dot = subject measured value, green dot = subject individual pseudo-normal value. Interpretation focuses on whether the subtype reference changes conclusions.')]
    g.append(fig('GMV_Subtype_specific.png') + interp(it('GMV_sub',
        '亚型女性中位数在 24 岁处约 525 ml、下界约 495 ml;受检者 374.7 ml 深陷下界之下,与中位数的差距约 150 ml,'
        '明显大于通用常模下的差距(约 105 ml)。结论方向不变(灰质不足)但严重度上调:即使以"结构指标更高的亚型人群"为参照,'
        '该灰质缺口仍远超亚型变异范围,属个体特异性偏差,而非亚型背景可解释的偏移。',
        'The subtype female median at 24 y is about 525 ml with a lower bound of about 495 ml; the subject\'s 374.7 ml falls far below the bound, with a gap of about 150 ml to the median, '
        'clearly larger than under the general norm (about 105 ml). The direction is unchanged (grey-matter insufficiency) but severity is upgraded: even against a "subtype population with higher structural metrics", '
        'the grey-matter deficit far exceeds subtype variation - an individual-specific deviation not explainable by subtype background.'), sub=True))
    g.append(fig('sGMV_Subtype_specific.png') + interp(it('sGMV_sub',
        '皮层下灰质在亚型参考下同样位于区间下缘或以下。皮层下核团是皮层环路的枢纽,其低值与海马、丘脑单核团结果互证;'
        '在亚型常模下该偏差亦不能被亚型背景吸收,提示皮层-皮层下两级灰质不足是个体层面的稳定特征。',
        'Subcortical grey matter also lies at or below the interval edge under the subtype reference. Subcortical nuclei are hubs of cortical circuits; this low value corroborates the single-nucleus results of hippocampus and thalamus. '
        'Under the subtype norm the deviation is likewise not absorbed by subtype background, indicating that two-level cortical-subcortical grey-matter insufficiency is a stable individual-level feature.'), sub=True))
    g.append(fig('WMV_Subtype_specific.png') + interp(it('WMV_sub',
        '白质体积在通用常模下居中(54.1 百分位);在位置更高的亚型参考下,同一实测值的相对位置下移,'
        '但多仍落于亚型区间之内。该指标是"常模选择改变百分位、但不改变定性结论"的典型例子:'
        '白质保有这一特征在两套常模下均成立,支持偏差的灰质选择性。',
        'White-matter volume is mid-range under the general norm (54.1 centile); under the higher-positioned subtype reference the relative position of the same measured value shifts down '
        'but mostly stays within the subtype interval. This metric is a typical example of "norm choice changes centiles but not the qualitative conclusion": '
        'preserved white matter holds under both norms, supporting grey-matter selectivity of the deviation.'), sub=True))
    g.append(fig('Ventricles_Subtype_specific.png') + interp(it('Ventricles_sub',
        '脑室体积的亚型轨迹与通用轨迹形态相近(随龄加速增长)。受检者 27.0 ml 在两套参考下均处区间上缘附近,'
        '作为实质体积偏低的代偿性镜像(ex vacuo)与纵向监测锚点的解读在亚型框架下同样成立。',
        'The subtype ventricular trajectory resembles the general one (accelerating growth with age). The subject\'s 27.0 ml lies near the upper edge under both references; '
        'the reading as a compensatory mirror of low parenchymal volume (ex vacuo) and as a longitudinal monitoring anchor holds equally in the subtype frame.'), sub=True))
    g.append(fig('mean_thickness_Subtype_specific.png') + interp(it('mean_thickness_sub',
        '平均皮层厚度在亚型参考下进一步远离中位数曲线。厚度是皮层完整性最敏感的宏观指标;'
        '双常模下一致的极低位置(通用常模第 0.03 百分位)使"弥散性皮层变薄"成为本例最稳健的表型结论,'
        '技术因素排除后应视为核心生物学特征。',
        'Mean cortical thickness moves further from the median curve under the subtype reference. Thickness is the most sensitive macroscopic marker of cortical integrity; '
        'the consistently extreme position under both norms (0.03rd centile under the general norm) makes "diffuse cortical thinning" the most robust phenotypic conclusion of this case, '
        'to be regarded as a core biological feature once technical factors are excluded.'), sub=True))
    g.append(fig('total_surface_area_Subtype_specific.png') + interp(it('total_surface_area_sub',
        '皮层总面积在亚型参考下的相对位置较通用常模下移,但仍多处于区间内。'
        '面积-厚度解离(切向扩张保有、放射维度不足)在亚型框架下依旧可见,支持发育机制定位的稳定性。',
        'The relative position of total cortical area under the subtype reference shifts down compared with the general norm but mostly remains within the interval. '
        'The area-thickness dissociation (preserved tangential expansion, insufficient radial dimension) persists in the subtype frame, supporting stability of the developmental-mechanism localisation.'), sub=True))
    g.append(fig('Cerebellum_Subtype_specific.png') + interp(it('Cerebellum_sub',
        '小脑亚型曲线峰值更高(女约 143 ml)、早期上升更陡;受检者 109.7 ml 落于亚型女性 95% 下界(24 岁约 118 ml)之下,'
        '而通用常模下仅贴近下界。即:小脑偏低不能由亚型背景解释,属个体特异性偏差。'
        '其经小脑-丘脑-皮层环路参与工作记忆、语言与情绪调节,与皮层变薄的共现提示共享的发育性背景。',
        'The subtype cerebellar curve peaks higher (female about 143 ml) with a steeper early rise; the subject\'s 109.7 ml falls below the subtype female 95% lower bound (about 118 ml at 24 y), '
        'whereas under the general norm it merely approaches the bound. Thus the low cerebellar volume cannot be explained by subtype background and is an individual-specific deviation. '
        'Via the cerebello-thalamo-cortical loop it participates in working memory, language and emotional regulation; its co-occurrence with cortical thinning suggests a shared developmental background.'), sub=True))
    g.append(fig('BrainStem_Subtype_specific.png') + interp(it('BrainStem_sub',
        '脑干体积在亚型参考下仍处区间中部附近,与通用常模结论一致(55.5 百分位):'
        '生命中枢与传导束形态保有,是整体模式中的一枚稳定阴性证据。',
        'Brainstem volume remains near the middle of the interval under the subtype reference, consistent with the general norm (55.5 centile): '
        'preserved vital centres and tracts - a stable piece of negative evidence in the overall pattern.'), sub=True))
    g.append('</section>')

    q = quant_section(3)
    _sub_native = (cen_sub_src == 'individual_centile_score_subtype_specific.csv')
    if _sub_native:
        _loll = ('皮层分区厚度百分位分布(受检者个体,亚型特异性常模口径)',
                 'Cortical parcel thickness centile distribution (subject, subtype-specific calibre)')
        _ttl = ('受检者定量百分位汇总(亚型特异性常模口径)',
                'Subject quantitative centile summary (subtype-specific calibre)')
        _note = bi(
            '下表百分位以<b>亚型特异性常模</b>计算,与本报告第 2 节曲线、第 4 节映射图同一常模口径;全人群通用常模下的百分位见报告(I),可并列判读。'
            '双常模并列的意义:若某指标在两套常模下一致偏离,则为稳健的个体特异性偏差;若仅在通用常模下偏离,则可能由亚型/队列背景解释。',
            'The centiles in the tables below are computed with the <b>subtype-specific norm</b>, the same calibre as the section 2 curves and section 4 maps of this report; centiles under the population-general norm are given in report (I) for side-by-side reading. '
            'The point of dual-norm comparison: if a metric deviates consistently under both norms it is a robust individual-specific deviation; if it deviates only under the general norm it may be explained by subtype/cohort background.')
    else:
        _loll = ('皮层分区厚度百分位分布(受检者个体,通用常模回退)',
                 'Cortical parcel thickness centile distribution (subject, general-norm fallback)')
        _ttl = ('受检者定量百分位汇总(通用常模回退)',
                'Subject quantitative centile summary (general-norm fallback)')
        _note = bi(
            '注意:亚型特异性百分位文件(individual_centile_score_subtype_specific.csv)暂缺,下表百分位暂以<b>全人群通用常模</b>回退显示(与报告(I)一致);待管线输出该文件后将自动切换为亚型特异性口径。',
            'Note: the subtype-specific centile export (individual_centile_score_subtype_specific.csv) is not present yet, so the centiles below temporarily fall back to the <b>population-general norm</b> (consistent with report (I)); they switch to the subtype-specific calibre automatically once the pipeline writes that file.')
    q = q.replace('皮层分区厚度百分位分布(受检者个体)</span>', _loll[0] + '</span>')
    q = q.replace('Cortical parcel thickness centile distribution (subject)</span>', _loll[1] + '</span>')
    q = q.replace('受检者定量百分位汇总</span>', _ttl[0] + '</span>')
    q = q.replace('Subject quantitative centile summary</span>', _ttl[1] + '</span>')
    q = q.replace('</h2>', '</h2><p class="mut">%s</p>' % _note, 1)

    m = ['<section id="maps"><h2><span class="no">4</span>%s</h2>' % bi('分区百分位映射图:亚型常模下的系统性红移', 'Regional centile mapping: systematic red shift under the subtype norm'),
         '<p class="mut">%s</p>' % bi('以下映射图为<b>亚型特异性常模</b>下的分区百分位映射。与报告(I)在通用常模下的映射图相比,'
         '红色(低百分位)占比系统性上升:皮层体积 17% → 67%、皮层厚度 54% → 64%、皮层面积 6% → 69%、皮层下核团 17% → 98%。'
         '该红移与"亚型曲线整体位置更高"的形态学观察方向一致,是<b>常模选择会系统性改变百分位结论</b>的直接可视化证据,'
         '也是个体化评估必须同时报告双常模的理由;鼠标悬停时,皮层表面四视图按 DKT atlas 显示脑区名称(英文),'
         '皮层下冠状面映射图按 FreeSurfer aseg 模板显示核团名称(英文)。',
         'The maps below show regional centile mapping under the <b>subtype-specific normative model</b>. Compared with the maps under the general norm in report (I), '
         'the red (low-centile) fraction rises systematically: cortical volume 17% -> 67%, cortical thickness 54% -> 64%, cortical area 6% -> 69%, subcortical nuclei 17% -> 98%. '
         'This red shift matches the morphological observation that subtype curves sit higher overall, and is direct visual evidence that <b>norm choice systematically changes centile conclusions</b> - '
         'the reason individualized assessment must report both norms; on hover, the four cortical surface views show DKT atlas region names (English) '
         'and the subcortical coronal map shows FreeSurfer aseg structure names (English).')]
    for name in MAPS:
        m.append(fig(name.replace('.png', '_Subtype_specific.png'), extra_class='map', max_dim=2000,
                     hover='aseg' if name.startswith('subcortical') else 'dkt'))
    m.append('</section>')

    t = ['<section id="tables"><h2><span class="no">5</span>%s</h2>' % bi('局部指标汇总表(亚型常模)', 'Regional summary tables (subtype norm)'),
         '<p class="mut">%s</p>' % bi('常模输出的分区汇总表(亚型常模口径),表头单位统一以括号标注(Volume(ml)、Volume(centile)、Thickness(mm)、Thickness(centile)、Area(cm²)、Area(centile))。',
         'Normative-model regional summary tables under the subtype norm; header units uniformly parenthesised (Volume(ml), Volume(centile), Thickness(mm), Thickness(centile), Area(cm2), Area(centile)).')]
    t.append(fig('global_summary_Subtype_specific.png', white_bg=False))
    t.append(fig('local_summary_Subtype_specific.png', extra_class='map', white_bg=False))
    t.append('</section>')

    discuss = ('<section id="discuss"><h2><span class="no">6</span>%s</h2>' % bi('医学与神经科学意义综合解读(亚型视角)', 'Medical and neuroscience integrated interpretation (subtype perspective)')
               + '<h3>%s</h3>' % bi('6.1 亚型对照的核心增益', '6.1 Subtype comparison: the core gain')
               + '<p>%s</p>' % it('sub61',
                   '亚型特异性常模把问题从"个体是否偏离人群"推进到"个体是否偏离其亚型背景的典型范围"。'
                   '本例中,总灰质、平均皮层厚度、小脑等核心指标在亚型参考下不仅未回归区间,反而偏差幅度扩大'
                   '(如小脑由"贴近通用下界"变为"低于亚型下界"),说明这些偏差<b>不能</b>由亚型/队列背景吸收,属个体特异性偏差;'
                   '而白质体积、皮层总面积等指标在亚型参考下仍居区间内,其"保有"结论对常模选择稳健。'
                   '这一分离正是亚型常模的临床价值所在:它把干预与随访的优先级集中到真正个体特异的轴上。',
                   'The subtype-specific norm advances the question from "does the individual deviate from the population" to "does the individual deviate from the typical range of his/her subtype background". '
                   'Here, core metrics such as total grey matter, mean cortical thickness and cerebellum do not return into the interval under the subtype reference; their deviation amplitudes even enlarge '
                   '(e.g. cerebellum changes from "near the general lower bound" to "below the subtype lower bound"), showing these deviations <b>cannot</b> be absorbed by subtype/cohort background and are individual-specific; '
                   'whereas white-matter volume and total cortical area remain within the interval, so their "preserved" conclusion is robust to norm choice. '
                   'This separation is the clinical value of subtype norms: it concentrates intervention and follow-up priority on the truly individual-specific axes.')
               + '<h3>%s</h3>' % bi('6.2 机制与临床含义', '6.2 Mechanisms and clinical implications')
               + '<ul class="tight">'
               + '<li>%s</li>' % it('sub62[0]', '<b>发育机制:</b>厚度-面积解离在双常模下稳定,支持放射单元假说(radial unit hypothesis)框架下柱内(放射)维度异常、切向扩张保有;海马(左 0.55、右 0.01 百分位)、小脑、丘脑的同时受累提示多系统发育程序异常早在围产期前后即已启动。',
                                   '<b>Developmental mechanism:</b> the thickness-area dissociation is stable under both norms, supporting intra-column (radial) dimension abnormality with preserved tangential expansion under the radial unit hypothesis; simultaneous involvement of hippocampus (L 0.55, R 0.01 centile), cerebellum and thalamus suggests a multi-system developmental programme anomaly starting around the perinatal period.')
               + '<li>%s</li>' % it('sub62[1]', '<b>退行可能性:</b>年轻年龄、白质保有、无区域选择性萎缩梯度,使典型退行性病程可能性低;但海马极低值需基线-随访对比排除极早发型过程。',
                                   '<b>Degenerative likelihood:</b> young age, preserved white matter and absence of a regionally selective atrophy gradient make a typical degenerative course unlikely; the extremely low hippocampal value still requires baseline-follow-up comparison to exclude a very-early-onset process.')
               + '<li>%s</li>' % it('sub62[2]', '<b>亚型分层:</b>若该受检者被归入某影像/疾病亚型,其"低于亚型下界"的轴(灰质、厚度、小脑、海马)应作为该亚型内的极端表型标记,用于预后富集与机制亚组分析。',
                                   '<b>Subtype stratification:</b> if this subject is assigned to an imaging/disease subtype, the axes "below the subtype lower bound" (grey matter, thickness, cerebellum, hippocampus) should serve as extreme phenotype markers within that subtype for prognostic enrichment and mechanistic subgroup analysis.')
               + '<li>%s</li></ul>' % it('sub62[3]', '<b>临床关联:</b>认知筛查偏低与海马-皮层低体积构成可检验关联;建议分域神经心理评估与 12–24 个月纵向复扫,以脑室百分位与平均厚度为锚点。',
                                         '<b>Clinical association:</b> low cognitive-screening performance and low hippocampal-cortical volumes form a testable association; domain-wise neuropsychological assessment and 12-24 month longitudinal repeat MRI are suggested, anchored on ventricular centile and mean thickness.')
               + '<h3>%s</h3>' % bi('6.3 判读注意事项', '6.3 Interpretation caveats')
               + '<div class="callout">%s</div>' % it('sub63',
                   '亚型常模的参考人群定义、样本量与分层稳定性直接决定区间宽度与位置;在亚型队列样本有限时,'
                   '区间估计的不确定性大于通用常模,本报告的红移/偏差结论应结合报告(I)并列判读。全部输出仅供科学研究,不构成诊断依据。',
                   'The reference population definition, sample size and stratification stability of a subtype norm directly determine interval width and position; with limited subtype cohort samples, '
                   'interval uncertainty exceeds that of the general norm, and the red-shift/deviation conclusions of this report should be read side by side with report (I). All outputs are for scientific research only and do not constitute diagnostic evidence.')
               + '</section>')

    body = toc + '<div class="wrap">' + s0 + SUBJ_BLOCK + ''.join(g) + q + ''.join(m) + ''.join(t) + discuss + conclusion_section(cen_sub, _sub_keys) + '</div>'
    return page('个体化脑健康结构影像报告(II)— 亚型特异性常模参照', '', hero + body)


set_centile(cen_pop)
_pop_html = build_population()
set_centile(cen_sub)
_sub_html = build_subtype()
_built = {'population-normative': _pop_html, 'subtype-specific': _sub_html}
for _kind, _html in _built.items():
    with open(os.path.join(BASE, 'Individualized-brain-health-report_%s.html' % _kind),
              'w', encoding='utf-8') as fh:
        fh.write(_html)
    for _lang in ('zh', 'en'):
        with open(os.path.join(BASE, 'Individualized-brain-health-report_%s_%s.html'
                               % (_kind, _lang)), 'w', encoding='utf-8') as fh:
            fh.write(mono(_html, _lang))
print('OK')
print('centile sources: population=%s | subtype=%s' % (cen_pop_src, cen_sub_src))
set_centile(cen_pop)
print('[pop] thickness <0.05:', n_th05, '<0.10:', n_th10, 'mean %.3f' % mean_th)
print('[pop] area mean %.3f >0.9:%d' % (mean_ar, n_ar90))
print('[pop] vol <0.10:', n_vol10)
print('[pop] lowest thickness:', low_th_names)
print('[pop] lowest volume:', low_vol_names)
set_centile(cen_sub)
print('[sub] thickness <0.05:', n_th05, '<0.10:', n_th10, 'mean %.3f' % mean_th)
print('[sub] area mean %.3f >0.9:%d' % (mean_ar, n_ar90))
print('[sub] vol <0.10:', n_vol10)
print('[sub] lowest thickness:', low_th_names)
print('[sub] lowest volume:', low_vol_names)
print('[pop] abnormal indicators:', abnormal_indicators(cen_pop, _pop_keys))
print('[sub] abnormal indicators:', abnormal_indicators(cen_sub, _sub_keys))
