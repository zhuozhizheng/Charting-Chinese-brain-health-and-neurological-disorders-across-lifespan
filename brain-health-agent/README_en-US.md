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

Open http://127.0.0.1:8000/ (top-right CH / EN toggle, Chinese by default).

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

---

# brain-health-agent

Individualized Brain Health Report Agent + Web Console (Bilingual Chinese-English, Default Chinese).
Input a folder path → automatically scan DICOM series → serially invoke the FreeSurfer pipeline → snapshot shared artifacts →
the interpretation engine (Qoder CLI or a configured LLM) generates interpretations → for population-normative and subtype-specific, respectively, generate two separate reports, one in Chinese and one in English (4 HTML files per subject, with a link at the top of each page to switch to the other language version).

## Directory Structure

- brain_agent/ agent package (discovery / pipeline / snapshot / llm / reports / runner / server / webui / cli)

- report_builder/build_reports.py report builder (copied by the agent into each subject's snapshot directory at runtime for execution)

- agent_out/ run status and artifacts (agent_out/<run_id>/<subject>/ = snapshot + 4 HTML files (Chinese/English × two norms) + logs; not committed to the repository)

- run.sh one-command startup for the web service

- brain_agent/config.json local path configuration (not committed; see config.example.json for the template)

- brain_agent/llm_config.json LLM configuration (base_url / api_key / model, permissions 600, not committed; see llm_config.example.json for the template)

## 启动

```bash
./run.sh 8000          # 或：conda activate py38 && python -m brain_agent.cli serve --port 8000
```

Open http://127.0.0.1:8000/ in your browser (use the "中文 / EN" toggle in the upper-right corner of the page; Chinese is the default).

## CLI

```bash
python -m brain_agent.cli discover <root>            
python -m brain_agent.cli run <root> [--subjects A,B] [--dry-run] [--wait]
python -m brain_agent.cli status <run_id>
python -m brain_agent.cli report-only <snapshot_dir> 
python -m brain_agent.cli serve [--port 8000]
```

## HTTP API

- GET / bilingual web console

- POST /api/runs {root, subjects?, dry_run?} / GET /api/runs / GET /api/runs/<id>

- GET /api/runs/<id>/files/<subject>/<name> download/open reports

- GET /api/tools OpenAI function-calling schema (for LLMs/CLI tools to invoke this agent)

- GET|POST /api/llm, POST /api/llm/test LLM configuration (masked echo) / connectivity test

- All JSON APIs support ?lang=zh|en (default zh)

## Interpretation Engine

- Default: local Qoder CLI in headless mode (qoder -p, tools disabled) generates bilingual Chinese-English interpretations in real time based on each subject's quantitative summary, requiring no key;

- After an OpenAI-compatible configuration (base_url / api_key / model) is saved in the "⚙ Settings" dialog in the upper-right corner of the page, that LLM is used instead;

- When neither is available, reports fall back to frozen demo template text with an explicit notice at the top of the page; if the engine is available but the call fails, that subject is marked as failed.

## Privacy Notes

- Subject identity information (ID / sex / age / examination date / acquisition center) is extracted at runtime by the agent from the DICOM header table (<series>_sub/T1.csv) into subject_meta.json, based on which the report builder renders the reports; the repository contains no real subject information, only anonymous demo fallback values.

- Local paths and API keys are stored in git-ignored configuration files; run artifacts in agent_out/ are not committed to the repository.

## Dependencies

- Service/agent: conda environment py38 (flask);

- Report building: system python3 (numpy + Pillow);

- Pipeline (external, not part of this project): freesurf_luhao.py, dcm2niix, FreeSurfer, Rscript (see config.json).



## Concurrency Protection

The pipeline's shared scratch directory (tem_data) is overwrite-based: before starting each subject, the agent checks for an external freesurf_luhao.py process and marks the subject as blocked if one exists; concurrent runs of the same agent are mutually excluded via agent_out/run.lock. Subjects in blocked/failed state can be re-launched with the same run_id to resume.
