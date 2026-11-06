# -*- coding: utf-8 -*-
"""Bilingual (zh default, EN toggle) web console for the brain-health agent."""

CONSOLE_HTML = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>脑健康报告 Agent 控制台</title>
<style>
:root{--ink:#1c2733;--mut:#5b6b7c;--bg:#f2f5f9;--card:#fff;--line:#e3e9f0;--acc:#2563eb;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.7 -apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
header{background:linear-gradient(135deg,#0f2b46,#155e94 55%,#0e7490);color:#fff;padding:34px 22px 26px}
header h1{margin:0 0 6px;font-size:24px}
header p{margin:0;opacity:.9;font-size:14px}
main{max-width:1000px;margin:22px auto 60px;padding:0 18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px 22px;margin-bottom:18px}
.card h2{font-size:17px;margin:0 0 12px;color:#123c63}
label{display:block;font-size:13.5px;color:var(--mut);margin:8px 0 4px}
input[type=text]{width:100%;padding:9px 12px;border:1px solid var(--line);border-radius:9px;font-size:14px}
button{cursor:pointer;border:0;background:var(--acc);color:#fff;border-radius:9px;padding:9px 18px;font-size:14px;margin-top:12px}
button.ghost{background:#eef4fb;color:#123c63}
.langbar{position:fixed;top:12px;right:14px;z-index:50;display:flex;gap:6px}
.langbar button{margin:0;padding:4px 13px;border-radius:999px;font-size:12.5px;
background:rgba(16,28,40,.55);border:1px solid rgba(255,255,255,.55);color:#fff}
.langbar button.on{background:#fff;color:#123c63;font-weight:600}
table{border-collapse:collapse;width:100%;font-size:13.5px;background:#fff}
th,td{border:1px solid var(--line);padding:6px 9px;text-align:left}
th{background:#eef4fb}
.mut{color:var(--mut);font-size:12.5px}
.chip{border-radius:999px;padding:1px 10px;font-size:12px;color:#fff}
.st-done{background:#27ae60}.st-failed{background:#c0392b}.st-running{background:#e67e22}
.st-planned{background:#8fa3b8}.st-blocked{background:#7f8c8d}
a{color:var(--acc)}
#msg{margin-left:12px;font-size:13.5px}
#settingsModal{position:fixed;inset:0;background:rgba(10,20,32,.55);z-index:100;
display:flex;align-items:center;justify-content:center;padding:18px}
#settingsModal .modal{background:#fff;border-radius:14px;padding:22px 26px;max-width:520px;width:100%;
max-height:88vh;overflow:auto;box-shadow:0 8px 40px rgba(0,0,0,.35)}
#settingsModal h2{font-size:17px;margin:0 0 8px;color:#123c63}
html[lang="en"] [data-lang="zh"]{display:none}
html:not([lang="en"]) [data-lang="en"]{display:none}
</style>
</head>
<body>
<div class="langbar"><button id="langZh" type="button">中文</button><button id="langEn" type="button">EN</button><button id="settingsBtn" type="button" title="Settings">⚙</button></div>
<div id="settingsModal" style="display:none">
<div class="modal">
<h2><span data-lang="zh">设置 · 大模型解读</span><span data-lang="en">Settings · LLM interpretation</span></h2>
<p class="mut"><span data-lang="zh">默认使用本机 Qoder CLI 生成解读;填写并保存 API 配置后改用该大模型。</span><span data-lang="en">Interpretations use the local Qoder CLI by default; once an API configuration is saved, this model is used instead.</span></p>
<label>Base URL (OpenAI 兼容 / compatible)</label>
<input id="llmUrl" type="text" placeholder="https://api.openai.com/v1"/>
<label>Model</label>
<input id="llmModel" type="text" placeholder="gpt-4o-mini"/>
<label>API Key <span id="llmMask" class="mut"></span></label>
<input id="llmKey" type="password" placeholder="sk-..."/>
<button id="llmSave"><span data-lang="zh">保存</span><span data-lang="en">Save</span></button>
<button id="llmTest" class="ghost"><span data-lang="zh">测试连接</span><span data-lang="en">Test</span></button>
<button id="settingsClose" class="ghost"><span data-lang="zh">关闭</span><span data-lang="en">Close</span></button>
<span id="llmMsg"></span>
</div>
</div>
<header>
<h1><span data-lang="zh">脑健康报告 Agent 控制台</span><span data-lang="en">Brain-Health Report Agent Console</span></h1>
<p><span data-lang="zh">输入文件夹路径 → 自动串行分析 → 每位受检者生成 population-normative 与 subtype-specific 两份报告</span><span data-lang="en">Folder path in → serial pipeline per subject → population-normative and subtype-specific reports out</span></p>
</header>
<main>
<div class="card">
<h2><span data-lang="zh">启动运行</span><span data-lang="en">Start a run</span></h2>
<label><span data-lang="zh">文件夹路径(含 DICOM 系列目录的根目录)</span><span data-lang="en">Folder path (root containing DICOM series directories)</span></label>
<input id="root" type="text" placeholder="/data/cohort/root"/>
<label><input id="dry" type="checkbox"/> <span data-lang="zh">仅预检(dry-run,只扫描不跑管线)</span><span data-lang="en">Discovery only (dry-run, no pipeline)</span></label>
<button id="go"><span data-lang="zh">启动</span><span data-lang="en">Run</span></button>
<button id="refresh" class="ghost"><span data-lang="zh">刷新</span><span data-lang="en">Refresh</span></button>
<span id="msg"></span>
</div>
<div class="card">
<h2><span data-lang="zh">运行列表</span><span data-lang="en">Runs</span></h2>
<div id="runs" class="mut"><span data-lang="zh">加载中…</span><span data-lang="en">loading…</span></div>
</div>
<div class="card" id="detailCard" style="display:none">
<h2><span data-lang="zh">受检者明细</span><span data-lang="en">Subject detail</span> <span id="detailId" class="mut"></span></h2>
<div id="detail"></div>
</div>
</main>
<script>
var LANG='zh';
try{LANG=localStorage.getItem('bh_lang')||'zh';}catch(e){}
function setLang(l){LANG=l;document.documentElement.lang=l;
try{localStorage.setItem('bh_lang',l);}catch(e){}
document.getElementById('langZh').className=(l==='zh'?'on':'');
document.getElementById('langEn').className=(l==='en'?'on':'');
loadRuns();if(openId)loadDetail(openId);}
document.getElementById('langZh').onclick=function(){setLang('zh');};
document.getElementById('langEn').onclick=function(){setLang('en');};
var ST={done:['完成','done'],failed:['失败','failed'],running:['运行中','running'],
planned:['排队中','planned'],blocked:['被阻塞','blocked']};
function stLabel(s){var p=ST[s]||[s,s];return '<span class="chip st-'+s+'">'+(LANG==='zh'?p[0]:p[1])+'</span>';}
function esc(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');}
var openId=null;
function loadRuns(){
fetch('/api/runs?lang='+LANG).then(function(r){return r.json();}).then(function(d){
var el=document.getElementById('runs');
if(!d.runs.length){el.innerHTML='<span data-lang="zh">暂无运行记录</span><span data-lang="en">no runs yet</span>';return;}
var h='<table><tr><th>run_id</th><th>'+ (LANG==='zh'?'根目录':'root') +'</th><th>'+ (LANG==='zh'?'受检者':'subjects') +'</th><th>'+ (LANG==='zh'?'完成':'done') +'</th><th>'+ (LANG==='zh'?'失败':'failed') +'</th><th>'+ (LANG==='zh'?'状态':'state') +'</th></tr>';
d.runs.forEach(function(r){
var state=r.dry_run?(LANG==='zh'?'预检':'dry-run'):(r.failed?(LANG==='zh'?'有失败':'has failures'):(r.done===r.subjects?(LANG==='zh'?'已完成':'finished'):(LANG==='zh'?'进行中':'running')));
h+='<tr style="cursor:pointer" data-id="'+esc(r.run_id)+'"><td>'+esc(r.run_id)+'</td><td class="mut">'+esc(r.root)+'</td><td>'+r.subjects+'</td><td>'+r.done+'</td><td>'+r.failed+'</td><td>'+esc(state)+'</td></tr>';});
el.innerHTML=h+'</table>';
el.querySelectorAll('tr[data-id]').forEach(function(tr){tr.onclick=function(){openId=tr.getAttribute('data-id');loadDetail(openId);};});
});}
function loadDetail(id){
fetch('/api/runs/'+encodeURIComponent(id)+'?lang='+LANG).then(function(r){return r.json();}).then(function(d){
document.getElementById('detailCard').style.display='';
document.getElementById('detailId').textContent=id;
var h='<table><tr><th>'+ (LANG==='zh'?'受检者':'subject') +'</th><th>'+ (LANG==='zh'?'状态':'status') +'</th><th>'+ (LANG==='zh'?'报告':'reports') +'</th><th>'+ (LANG==='zh'?'备注':'note') +'</th></tr>';
(d.subjects||[]).forEach(function(s){
var links='';
if(s.status==='done'){
links='<a target="_blank" href="/api/runs/'+encodeURIComponent(id)+'/files/'+encodeURIComponent(s.subject)+'/Individualized-brain-health-report_population-normative_zh.html">pop 中</a> · '
+'<a target="_blank" href="/api/runs/'+encodeURIComponent(id)+'/files/'+encodeURIComponent(s.subject)+'/Individualized-brain-health-report_population-normative_en.html">pop EN</a> · '
+'<a target="_blank" href="/api/runs/'+encodeURIComponent(id)+'/files/'+encodeURIComponent(s.subject)+'/Individualized-brain-health-report_subtype-specific_zh.html">sub 中</a> · '
+'<a target="_blank" href="/api/runs/'+encodeURIComponent(id)+'/files/'+encodeURIComponent(s.subject)+'/Individualized-brain-health-report_subtype-specific_en.html">sub EN</a>';}
h+='<tr><td>'+esc(s.subject)+'</td><td>'+stLabel(s.status)+'</td><td>'+links+'</td><td class="mut">'+esc(s.error||'')+'</td></tr>';});
document.getElementById('detail').innerHTML=h+'</table>';
});}
document.getElementById('go').onclick=function(){
var root=document.getElementById('root').value.trim();
var msg=document.getElementById('msg');
if(!root){msg.textContent=LANG==='zh'?'请填写文件夹路径':'please enter a folder path';return;}
fetch('/api/runs?lang='+LANG,{method:'POST',headers:{'Content-Type':'application/json'},
body:JSON.stringify({root:root,dry_run:document.getElementById('dry').checked})})
.then(function(r){return r.json();}).then(function(d){
msg.textContent=(d.message||d.error||'')+(d.run_id?'  ['+d.run_id+']':'');
loadRuns();if(d.run_id){openId=d.run_id;loadDetail(d.run_id);}});};
document.getElementById('refresh').onclick=function(){loadRuns();if(openId)loadDetail(openId);};
function loadLLM(){fetch('/api/llm?lang='+LANG).then(function(r){return r.json();}).then(function(d){
document.getElementById('llmUrl').value=d.base_url||'';
document.getElementById('llmModel').value=d.model||'';
document.getElementById('llmMask').textContent=d.api_key_masked||'';});}
var modal=document.getElementById('settingsModal');
document.getElementById('settingsBtn').onclick=function(){modal.style.display='flex';loadLLM();};
document.getElementById('settingsClose').onclick=function(){modal.style.display='none';};
modal.onclick=function(e){if(e.target===modal)modal.style.display='none';};
document.getElementById('llmSave').onclick=function(){
var body={base_url:document.getElementById('llmUrl').value.trim(),
model:document.getElementById('llmModel').value.trim()};
var k=document.getElementById('llmKey').value.trim();if(k)body.api_key=k;
fetch('/api/llm?lang='+LANG,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)})
.then(function(r){return r.json();}).then(function(d){
document.getElementById('llmMsg').textContent=d.message||d.error||'';
document.getElementById('llmMask').textContent=d.api_key_masked||'';
document.getElementById('llmKey').value='';loadLLM();});};
document.getElementById('llmTest').onclick=function(){
var m=document.getElementById('llmMsg');m.textContent=LANG==='zh'?'测试中…':'testing…';
fetch('/api/llm/test?lang='+LANG,{method:'POST'}).then(function(r){return r.json();}).then(function(d){
m.textContent=d.reply?((LANG==='zh'?'连接成功: ':'connected: ')+d.reply):((LANG==='zh'?'失败: ':'failed: ')+d.error);});};
setLang(LANG==='en'?'en':'zh');
loadLLM();
setInterval(function(){loadRuns();if(openId)loadDetail(openId);},4000);
</script>
</body>
</html>
"""
