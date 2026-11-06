# brain-health-agent

个体化脑健康报告 agent + Web 控制台 (中英双语，默认中文)。
输入一个文件夹路径 → 自动扫描 DICOM 系列 → 串行调用 FreeSurfer 管线 → 快照共享产物 →
解读引擎 (Qoder CLI 或已配置大模型) 生成解读 → 为 population-normative 与 subtype-specific
各生成**中、英两份独立报告**(每受检者 4 个 HTML,页顶提供另一语言版本的切换链接)。

## 目录结构

- `brain_agent/` agent 包 (discovery / pipeline / snapshot / llm / reports / runner / server / webui / cli)
- `report_builder/build_reports.py` 报告构建器 (agent 运行时复制到每个受检者快照目录执行)
- `agent_out/` 运行状态与产物 (`agent_out/<run_id>/<subject>/` = 快照 + 4 份 HTML(中/英 × 两常模)+ 日志;不入库)
- `run.sh` 一键启动 web 服务
- `brain_agent/config.json` 本机路径配置 (不入库;模板见 `config.example.json`)
- `brain_agent/llm_config.json` 大模型配置 (base_url / api_key / model，权限 600，不入库;模板见 `llm_config.example.json`)

## 启动

```bash
./run.sh 8000          # 或：conda activate py38 && python -m brain_agent.cli serve --port 8000
```

浏览器打开 http://127.0.0.1:8000/ (页面右上「中文 / EN」切换，默认中文)。

## CLI

```bash
python -m brain_agent.cli discover <root>            # 只读预检
python -m brain_agent.cli run <root> [--subjects A,B] [--dry-run] [--wait]
python -m brain_agent.cli status <run_id>
python -m brain_agent.cli report-only <snapshot_dir> # 仅重建两份 HTML
python -m brain_agent.cli serve [--port 8000]
```

## HTTP API

- `GET /` 双语 web 控制台
- `POST /api/runs {root, subjects?, dry_run?}` / `GET /api/runs` / `GET /api/runs/<id>`
- `GET /api/runs/<id>/files/<subject>/<name>` 报告下载/打开
- `GET /api/tools` OpenAI function-calling schema(供大模型/CLI 工具调用本 agent)
- `GET|POST /api/llm`、`POST /api/llm/test` 大模型配置 (掩码回显)/连通测试
- 所有 JSON 接口支持 `?lang=zh|en`(默认 zh)

## 解读引擎

- 默认：本机 **Qoder CLI** 无头模式 (`qoder -p`,禁用工具) 按受检者定量摘要实时生成中英双语解读，无需任何 key;
- 在页面右上「⚙ 设置」弹窗中保存 OpenAI 兼容配置 (base_url / api_key / model) 后，改用该大模型;
- 两者均不可用时，报告回退到冻结的演示模板文本并在页首显式标注;引擎可用但调用失败时该受检者判 failed。

## 隐私说明

- 受检者身份信息 (ID/性别/年龄/检查日期/采集中心) 在运行时由 agent 从 DICOM 头表
  (`<series>_sub/T1.csv`) 提取为 `subject_meta.json`,报告构建器据此渲染;仓库内不含任何
  真实受检者信息，仅保留匿名演示回退值。
- 本机路径与 API key 均在 git 忽略的配置文件中;运行产物 `agent_out/` 不入库。

## 依赖

- 服务/agent:conda 环境 `py38`(flask);
- 报告构建:系统 `python3`(numpy + Pillow);
- 管线 (外部，不在本项目内):`freesurf_luhao.py`、dcm2niix、FreeSurfer、Rscript(见 config.json)。

## 并发保护

管线共享 scratch 目录 (`tem_data`) 为覆盖式:agent 启动每个受检者前检测外部
`freesurf_luhao.py` 进程，存在则标记 `blocked`;同一 agent 并发运行由 `agent_out/run.lock` 互斥。
`blocked/failed` 的受检者可用相同 run_id 重新发起以续跑。

---

# brain-health-agent (English)

An individualized brain-health report agent with a bilingual web console (Simplified Chinese by
default, English toggle). Given a folder path it scans DICOM series, runs the FreeSurfer pipeline
serially per subject, snapshots the shared pipeline outputs, asks the interpretation engine
(Qoder CLI by default, or a configured model) for the medical/neuroscience interpretation, and
emits **separate Chinese and English reports** for both population-normative and subtype-specific
norms (four HTML files per subject, each with a top-right link to the other language version).

## Layout

- `brain_agent/` the agent package (discovery / pipeline / snapshot / llm / reports / runner / server / webui / cli)
- `report_builder/build_reports.py` report builder (copied into each subject snapshot dir at run time)
- `agent_out/` run state and outputs (`agent_out/<run_id>/<subject>/` = snapshot + four HTML reports (zh/en x two norms) + logs; not in git)
- `run.sh` one-command web service launcher
- `brain_agent/config.json` machine-local paths (not in git; template: `config.example.json`)
- `brain_agent/llm_config.json` LLM settings (base_url / api_key / model, mode 600, not in git; template: `llm_config.example.json`)

## Start

```bash
./run.sh 8000          # or: conda activate py38 && python -m brain_agent.cli serve --port 8000
```

Open http://127.0.0.1:8000/ (top-right 中文 / EN toggle, Chinese by default).

## CLI

```bash
python -m brain_agent.cli discover <root>            # read-only pre-check
python -m brain_agent.cli run <root> [--subjects A,B] [--dry-run] [--wait]
python -m brain_agent.cli status <run_id>
python -m brain_agent.cli report-only <snapshot_dir> # rebuild the two HTML reports only
python -m brain_agent.cli serve [--port 8000]
```

## HTTP API

- `GET /` bilingual web console
- `POST /api/runs {root, subjects?, dry_run?}` / `GET /api/runs` / `GET /api/runs/<id>`
- `GET /api/runs/<id>/files/<subject>/<name>` open/download a report
- `GET /api/tools` OpenAI function-calling schema so LLMs/CLI tools can drive this agent
- `GET|POST /api/llm`, `POST /api/llm/test` LLM settings (masked key echo) / connectivity test
- All JSON endpoints accept `?lang=zh|en` (default zh)

## Interpretation engines

- Default: the local **Qoder CLI** in headless mode (`qoder -p`, tools disabled) generates the
  bilingual interpretation live from each subject's quantitative summary - no key required;
- Saving an OpenAI-compatible configuration (base_url / api_key / model) in the gear settings
  modal (top right) switches interpretation to that model;
- If neither is available the report falls back to frozen demo template text with a banner; if an
  engine is available but fails, the subject is marked failed.

## Privacy

- Subject identifiers (id/sex/age/study date/site) are extracted at run time from the DICOM header
  table (`<series>_sub/T1.csv`) into `subject_meta.json` and rendered by the report builder; the
  repository contains no real subject data, only anonymised demo fallbacks.
- Machine paths and API keys live in git-ignored config files; run outputs under `agent_out/` are
  not committed.

## Dependencies

- service/agent: conda env `py38` (flask);
- report build: system `python3` (numpy + Pillow);
- pipeline (external, not part of this project): `freesurf_luhao.py`, dcm2niix, FreeSurfer, Rscript (see config.json).

## Concurrency guard

The pipeline's shared scratch dir (`tem_data`) is overwrite-style: before each subject the agent
checks for an external `freesurf_luhao.py` process and marks the subject `blocked` if one is
running; concurrent agent runs are mutually excluded via `agent_out/run.lock`. Blocked or failed
subjects can be retried by re-launching the same run_id.
