# -*- coding: utf-8 -*-
CATALOG = {
    'ok':                {'zh': '服务正常', 'en': 'service ok'},
    'run_started':       {'zh': '运行已启动', 'en': 'run started'},
    'run_planned':       {'zh': '运行已规划(dry-run)', 'en': 'run planned (dry-run)'},
    'run_not_found':     {'zh': '运行不存在', 'en': 'run not found'},
    'root_missing':      {'zh': '根目录不存在或不是目录', 'en': 'root missing or not a directory'},
    'path_escape':       {'zh': '非法路径参数', 'en': 'illegal path parameter'},
    'file_not_found':    {'zh': '文件不存在', 'en': 'file not found'},
    'bad_lang':          {'zh': '不支持的语言,已回退中文', 'en': 'unsupported language, fell back to zh'},
    'llm_not_configured': {'zh': '未配置大模型(base_url/api_key/model)', 'en': 'LLM not configured (base_url/api_key/model)'},
    'llm_saved':         {'zh': '大模型配置已保存', 'en': 'LLM configuration saved'},
}


def msg(key, lang='zh'):
    entry = CATALOG.get(key, {})
    return entry.get(lang, entry.get('zh', key))
