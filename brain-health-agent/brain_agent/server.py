# -*- coding: utf-8 -*-
import json
import os
import time

from flask import Flask, jsonify, request, send_from_directory

from . import llm
from .i18n import msg
from .runner import Runner, TERMINAL
from .webui import CONSOLE_HTML

TOOLS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools.json')


def _lang(req):
    lang = req.args.get('lang')
    if not lang and req.is_json:
        lang = (req.get_json(silent=True) or {}).get('lang')
    return lang if lang in ('zh', 'en') else 'zh'


def _clean(part):
    return part and '/' not in part and '\\' not in part and part not in ('.', '..')


def create_app(cfg):
    app = Flask(__name__)
    runner = Runner(cfg)

    @app.get('/')
    def console():
        return CONSOLE_HTML, 200, {'Content-Type': 'text/html; charset=utf-8'}

    @app.get('/health')
    def health():
        lang = _lang(request)
        return jsonify({'status': 'ok', 'lang': lang, 'message': msg('ok', lang),
                        'version': '0.1.0', 'time': time.time()})

    @app.post('/api/runs')
    def start():
        lang = _lang(request)
        body = request.get_json(silent=True) or {}
        root = body.get('root') or ''
        if not root or not os.path.isdir(root):
            return jsonify({'error': msg('root_missing', lang)}), 400
        run_id = runner.start_run(root, body.get('subjects'),
                                  bool(body.get('dry_run')), body.get('run_id'))
        key = 'run_planned' if body.get('dry_run') else 'run_started'
        return jsonify({'run_id': run_id, 'message': msg(key, lang),
                        'status_url': '/api/runs/%s' % run_id}), 202

    @app.get('/api/runs')
    def runs():
        lang = _lang(request)
        out = []
        for rid in runner.list_runs():
            st = runner.load_state(rid)
            out.append({'run_id': rid, 'root': st.get('root'),
                        'created': st.get('created'), 'dry_run': st.get('dry_run'),
                        'subjects': len(st.get('subjects', [])),
                        'done': sum(1 for s in st.get('subjects', []) if s['status'] == 'done'),
                        'failed': sum(1 for s in st.get('subjects', []) if s['status'] == 'failed')})
        return jsonify({'lang': lang, 'runs': out})

    @app.get('/api/runs/<run_id>')
    def run_detail(run_id):
        lang = _lang(request)
        if run_id not in runner.list_runs():
            return jsonify({'error': msg('run_not_found', lang)}), 404
        st = runner.load_state(run_id)
        st['lang'] = lang
        st['running'] = (not st.get('dry_run')) and any(
            s['status'] not in TERMINAL for s in st.get('subjects', []))
        return jsonify(st)

    @app.get('/api/runs/<run_id>/files/<subject>/<name>')
    def run_file(run_id, subject, name):
        lang = _lang(request)
        if not (_clean(subject) and _clean(name)) or run_id not in runner.list_runs():
            return jsonify({'error': msg('path_escape', lang)}), 400
        directory = os.path.join(runner.run_dir(run_id), subject)
        if not os.path.exists(os.path.join(directory, name)):
            return jsonify({'error': msg('file_not_found', lang)}), 404
        return send_from_directory(directory, name)

    @app.get('/api/tools')
    def tools():
        with open(TOOLS_PATH, encoding='utf-8') as fh:
            return jsonify(json.load(fh))

    @app.get('/api/llm')
    def llm_get():
        return jsonify(llm.masked(llm.load_config()))

    @app.post('/api/llm')
    def llm_set():
        lang = _lang(request)
        body = request.get_json(silent=True) or {}
        cfg = llm.load_config() or {}
        for key in ('base_url', 'model', 'timeout'):
            if body.get(key) is not None:
                cfg[key] = body[key]
        if body.get('api_key'):
            cfg['api_key'] = body['api_key']
        if body.get('clear'):
            cfg.pop('api_key', None)
        llm.save_config(cfg)
        out = llm.masked(cfg)
        out['message'] = msg('llm_saved', lang)
        return jsonify(out)

    @app.post('/api/llm/test')
    def llm_test():
        lang = _lang(request)
        cfg = llm.load_config()
        if not cfg:
            return jsonify({'error': msg('llm_not_configured', lang)}), 400
        try:
            return jsonify({'reply': llm.ping(cfg)})
        except Exception as exc:
            return jsonify({'error': str(exc)}), 502

    return app
