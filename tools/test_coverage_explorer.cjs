#!/usr/bin/env node
'use strict';
// Usage: node tools/test_coverage_explorer.cjs [report.html]
//        [--screenshot-dir directory] [--plotly-js local/plotly.min.js]
// Requires Playwright + Chromium (installed Chrome/Edge are also supported).
// Always tests independent synthetic fixtures, then the optional real report.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const os=require('node:os');
const path=require('node:path');
const {pathToFileURL}=require('node:url');
let chromium;
try{({chromium}=require('playwright'));}
catch{({chromium}=require(path.resolve(path.dirname(process.execPath),'../node_modules/playwright')));}

const template=fs.readFileSync(path.join(__dirname,'../web/coverage_explorer.html'),'utf8');
const names=['perigee_altitude_km','inclination_deg','eccentricity'];
const axisNames=['perigee_altitudes_km','inclinations_deg','eccentricities'];
const near=(actual,expected,message)=>assert(Math.abs(actual-expected)<=1e-11*Math.max(1,Math.abs(expected)),message+': '+actual+' != '+expected);
function syntheticData(extended=false){
 const metadata={samples:20000,seed:20260919,perigee_altitudes_km:extended?[1000,2000,3000,7000,12500,20000,30000]:[1000,2000,3000],inclinations_deg:[0,45,90],eccentricities:[0,.05,.1],position_sigmas_km:[1,10,100],velocity_sigmas_m_s:[.01,.1,1],coverage_definition:'Synthetic test criterion, not simulation results.',model_scope:'Synthetic browser test values only.',synthetic_test_data:true};
 if(extended){metadata.altitude_interpolation='log_geocentric_perigee_radius';metadata.earth_radius_km=6378.137;}
 const runs=[];
 for(const a of metadata.perigee_altitudes_km)for(const i of metadata.inclinations_deg)for(const e of metadata.eccentricities)for(const p of metadata.position_sigmas_km)for(const v of metadata.velocity_sigmas_m_s){const t=Math.exp(a/5000+i/180+e*2)*70/Math.hypot(p,v*3);runs.push({run_id:String(runs.length),samples:20000,perigee_altitude_km:a,inclination_deg:i,eccentricity:e,position_sigma_km:p,velocity_sigma_m_s:v,status:'ok',coverage_time_days:t,coverage_lower_days:.99*t,coverage_upper_days:t,min_initial_perigee_km:a-8*p,particles_below_300km:a-8*p<300?5:0,particles_below_500km:a-8*p<500?15:0});}
 const data={metadata,runs};
 if(extended){
  const referenceMeta={...metadata,perigee_altitudes_km:[19650.663],inclinations_deg:[55],eccentricities:[.02]};
  const referenceRuns=runs.slice(0,9).map(r=>({...r,run_id:'reference_'+r.run_id,perigee_altitude_km:19650.663,inclination_deg:55,eccentricity:.02,coverage_time_days:643.539/Math.hypot(r.position_sigma_km,r.velocity_sigma_m_s),coverage_lower_days:637.992/Math.hypot(r.position_sigma_km,r.velocity_sigma_m_s),coverage_upper_days:643.539/Math.hypot(r.position_sigma_km,r.velocity_sigma_m_s)}));
  data.reference_surface={label:'Original MEO / PDF',description:'Synthetic reference values preserved exactly.',metadata:referenceMeta,runs:referenceRuns};
  metadata.initial_view='reference';
 }
 return data;
}
function embeddedData(html){const match=html.match(/<script id="coverage-data" type="application\/json">([\s\S]*?)<\/script>/);assert(match,'Report must embed the result grid.');return JSON.parse(match[1]);}
function localLibrary(explicit,report){
 if(explicit)return fs.readFileSync(path.resolve(explicit),'utf8');
 if(report){const html=fs.readFileSync(report,'utf8'),match=html.match(/<script>([\s\S]*?)<\/script>/);assert(match&&match[1].includes('Plotly'),'Report must embed Plotly.');return match[1];}
 const candidates=[process.env.PLOTLY_JS,path.resolve(__dirname,'../../outputs/coverage_surface_20000/.deps/plotly/package_data/plotly.min.js'),path.resolve(__dirname,'../external/plotly.min.js')].filter(Boolean);
 const file=candidates.find(candidate=>fs.existsSync(candidate));
 assert(file,'Supply --plotly-js or a generated report containing Plotly.');return fs.readFileSync(file,'utf8');
}
function writeFixture(directory,name,data,library){
 const payload=JSON.stringify(data).replaceAll('&','\\u0026').replaceAll('<','\\u003c').replaceAll('>','\\u003e');
 const file=path.join(directory,name+'.html');
 fs.writeFileSync(file,template.replace('__COVERAGE_EXPLORER_PLOTLY__',()=>library.replaceAll('</script','<\\/script')).replace('__COVERAGE_EXPLORER_DATA__',()=>payload));return file;
}
const key=values=>values.map(v=>Number(v).toPrecision(12)).join(',');
function lookup(data){return new Map(data.runs.map(r=>[key([...names.map(n=>r[n]),r.position_sigma_km,r.velocity_sigma_m_s]),r]));}
function assertMeasured(g,data){
 assert.equal(g.mode,'measured');const table=lookup(data),meta=data.metadata;
 for(let v=0;v<meta.velocity_sigmas_m_s.length;v++)for(let p=0;p<meta.position_sigmas_km.length;p++){
  const r=table.get(key([...names.map(n=>g.values[n]),meta.position_sigmas_km[p],meta.velocity_sigmas_m_s[v]])),c=g.cells[v][p];
  if(!r||r.status!=='ok'){assert.equal(c.time,null);assert(c.reasons.includes(r?.status||'missing'));}
  else{assert.equal(c.status,'measured');assert.equal(c.time,r.coverage_time_days);assert.equal(c.lower,r.coverage_lower_days??null);assert.equal(c.upper,r.coverage_upper_days??null);}
 }
}
function assertReference(g,data){
 assert.equal(g.mode,'reference');assert.equal(g.corners.length,0);
 const reference=data.reference_surface;
 for(const [j,name] of names.entries())assert.equal(g.values[name],reference.metadata[axisNames[j]][0]);
 for(const row of g.cells)for(const c of row)assert.equal(c.status,'reference');
 assertMeasured({...g,mode:'measured',cells:g.cells.map(row=>row.map(c=>({...c,status:'measured'})))},reference);
}
function firstCellMidpoint(meta){
 return Object.fromEntries(names.map((name,j)=>{const axis=meta[axisNames[j]],radius=meta.earth_radius_km??6378.137;return [name,j===0&&meta.altitude_interpolation==='log_geocentric_perigee_radius'?Math.sqrt((radius+axis[0])*(radius+axis[1]))-radius:(axis[0]+axis[1])/2];}));
}
function assertMidpoint(g,data){
 const meta=data.metadata,table=lookup(data),axes=axisNames.map(n=>meta[n]);
 assert.equal(g.mode,'interpolated');assert.equal(g.corners.length,8);
 for(let v=0;v<meta.velocity_sigmas_m_s.length;v++)for(let p=0;p<meta.position_sigmas_km.length;p++){
  const rs=[];for(let a=0;a<2;a++)for(let i=0;i<2;i++)for(let e=0;e<2;e++)rs.push(table.get(key([axes[0][a],axes[1][i],axes[2][e],meta.position_sigmas_km[p],meta.velocity_sigmas_m_s[v]])));
  const c=g.cells[v][p],bad=rs.filter(r=>!r||r.status!=='ok');
  if(bad.length){assert.equal(c.time,null);assert.deepEqual([...c.reasons].sort(),[...new Set(bad.map(r=>r?.status||'missing'))].sort());}
  else{assert.equal(c.status,'interpolated');for(const [output,input] of [['time','coverage_time_days'],['lower','coverage_lower_days'],['upper','coverage_upper_days']]){
   const values=rs.map(r=>r[input]);if(values.some(x=>x==null)||values.some(x=>x===0)&&!values.every(x=>x===0))assert.equal(c[output],null);else if(values.every(x=>x===0))assert.equal(c[output],0);else near(c[output],Math.exp(values.reduce((s,x)=>s+Math.log(x),0)/8),'Independent geometric-mean oracle');
  }}
 }
}
async function openReport(browser,file){
 const page=await browser.newPage({viewport:{width:1440,height:1080},deviceScaleFactor:1}),errors=[],network=[];
 page.on('pageerror',error=>errors.push(error.message));page.on('request',request=>{if(!['file:','data:','about:','blob:'].includes(new URL(request.url()).protocol))network.push(request.url());});
 await page.goto(pathToFileURL(file).href);await page.waitForFunction(()=>window.coverageExplorer?.ready===true);assert.equal(await page.locator('#error').isVisible(),false);
 return {page,errors,network};
}
async function waitRender(page){await page.waitForFunction(()=>window.coverageExplorer.ready===true);}
async function assertResponsive(page,shots,prefix){
 let desktopEyeDistance=null;
 for(const [label,width,height] of [['desktop',1440,1080],['tablet',768,1024],['mobile',390,844]]){
  await page.setViewportSize({width,height});
  // Wait for Plotly's asynchronous resize and a painted frame, rather than a
  // timing guess that captures stale WebGL geometry on a busy simulation host.
  await page.evaluate(async()=>{await Plotly.Plots.resize(document.getElementById('surface'));await new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));});
  assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),label+' must not overflow horizontally.');
  assert(await page.locator('#surface canvas').count()>0,'Plotly must render a canvas.');
  const eyeDistance=await page.evaluate(()=>{const c=document.getElementById('surface')._fullLayout.scene._scene.getCamera();return Math.hypot(c.eye.x,c.eye.y,c.eye.z);});
  if(label==='desktop')desktopEyeDistance=eyeDistance;
  if(label==='mobile')assert(eyeDistance>desktopEyeDistance*1.2,'Mobile camera must pull back to keep 3D labels in view.');
  if(shots)await page.screenshot({path:path.join(shots,prefix+'-'+label+'.png'),fullPage:true});
 }
}
async function checkReport(browser,file,data,shots,prefix){
 const {page,errors,network}=await openReport(browser,file);
 try{
  const axes=axisNames.map(n=>data.metadata[n]);
  let g=await page.evaluate(()=>coverageExplorer.lastGrid);
  if(data.metadata.initial_view==='reference')assertReference(g,data);else assertMeasured(g,data);
  assert.equal(Number(await page.locator('#altitude').getAttribute('max')),axes[0].at(-1));
  // A fitted scale must follow the selected surface; a shared scale includes
  // all grid and reference cases without changing any measured values.
  const plottedMax=()=>page.evaluate(()=>document.getElementById('surface')._fullLayout.scene.zaxis.range[1]);
  near(await plottedMax(),g.max*1.04,'Initial fitted height scale');
  await page.locator('#shared-scale').check();await waitRender(page);
  const allRuns=[...data.runs,...(data.reference_surface?.runs??[])];
  const maximum=Math.max(1,...allRuns.filter(r=>r.status==='ok').map(r=>r.coverage_time_days));
  near(await plottedMax(),maximum*1.04,'Shared height scale');
  near(await page.evaluate(()=>document.getElementById('surface').data[0].cmax),maximum,'Shared color scale');
  await page.locator('#shared-scale').uncheck();await waitRender(page);
  if(data.reference_surface){
   await page.locator('#snap').check();await waitRender(page);
   await page.click('#reference-orbit');await waitRender(page);g=await page.evaluate(()=>coverageExplorer.lastGrid);assertReference(g,data);
   assert((await page.locator('#plot-state').textContent()).includes('no interpolation'));
   const [referenceDownload]=await Promise.all([page.waitForEvent('download'),page.click('#download')]);assert.equal(referenceDownload.suggestedFilename(),'coverage_surface_reference.csv');
   const csv=fs.readFileSync(await referenceDownload.path(),'utf8');assert(csv.includes(',reference,'));assert(csv.includes(String(data.reference_surface.runs[0].coverage_time_days)));
   if(shots)await page.screenshot({path:path.join(shots,prefix+'-reference.png'),fullPage:true});
   await page.locator('#altitude').fill(String(axes[0].at(-1)));await waitRender(page);assertMeasured(await page.evaluate(()=>coverageExplorer.lastGrid),data);
   await page.locator('#snap').uncheck();await waitRender(page);
  }else assert.equal(await page.locator('#reference-orbit').isVisible(),false);
  await page.click('#reset-orbit');await waitRender(page);
  // Exercise the actual input controls, including both ends of every range.
  for(const [j,id] of ['altitude','inclination','eccentricity'].entries())for(const value of [axes[j][0],axes[j].at(-1)]){
   await page.locator('#'+id).fill(String(value));await waitRender(page);g=await page.evaluate(()=>coverageExplorer.lastGrid);near(g.values[names[j]],value,'Slider value');assertMeasured(g,data);
  }
  if(axes.every(a=>a.length>1)){
   const midpoint=firstCellMidpoint(data.metadata);await page.evaluate(p=>coverageExplorer.setParameters(p),midpoint);g=await page.evaluate(()=>coverageExplorer.lastGrid);assertMidpoint(g,data);
   assert((await page.locator('#node-note').textContent()).includes('No new simulation'));
   // User camera changes must survive parameter updates (uirevision contract).
   const chosen={eye:{x:1.9,y:1.2,z:1.5},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}};
   await page.evaluate(camera=>Plotly.relayout('surface',{'scene.camera':camera}),chosen);
   await page.evaluate(p=>coverageExplorer.setParameters(p),midpoint);
   const actual=await page.evaluate(()=>document.getElementById('surface')._fullLayout.scene.camera);
   for(const n of ['x','y','z'])near(actual.eye[n],chosen.eye[n],'Preserve user camera');
   await page.click('#reset-camera');assert.notDeepEqual(await page.evaluate(()=>document.getElementById('surface')._fullLayout.scene.camera.eye),chosen.eye);
   await page.locator('#snap').check();await waitRender(page);g=await page.evaluate(()=>coverageExplorer.lastGrid);assertMeasured(g,data);await page.locator('#snap').uncheck();await waitRender(page);
  }
  await page.click('#reset-orbit');await waitRender(page);g=await page.evaluate(()=>coverageExplorer.lastGrid);assertMeasured(g,data);near(await plottedMax(),g.max*1.04,'Updated fitted height scale');
  const [download]=await Promise.all([page.waitForEvent('download'),page.click('#download')]);assert.equal(download.suggestedFilename(),'coverage_surface_measured.csv');
  await assertResponsive(page,shots,prefix);
  assert.deepEqual(errors,[],'No browser errors.');assert.deepEqual(network,[],'No network requests from the offline report.');
 }finally{await page.close();}
}
async function checkHoles(browser,file,data){
 const {page,errors}=await openReport(browser,file);
 try{
  await page.evaluate(()=>coverageExplorer.setParameters({perigee_altitude_km:1500,inclination_deg:22.5,eccentricity:.025}));
  let g=await page.evaluate(()=>coverageExplorer.lastGrid);assertMidpoint(g,data);assert.equal(g.validCount,4);assert.equal(g.totalCount,9);
  const warning=await page.locator('#warning').textContent();assert(warning.includes('required corner is unavailable'));assert(warning.includes('does not classify the selected cloud'));
  await page.evaluate(()=>coverageExplorer.setParameters({perigee_altitude_km:1000,inclination_deg:0,eccentricity:0}));g=await page.evaluate(()=>coverageExplorer.lastGrid);assertMeasured(g,data);assert((await page.locator('#warning').textContent()).includes('measured ensembles unavailable'));
  // Drag-like bursts must converge on the final parameter rather than race.
  await page.evaluate(()=>Promise.all(Array.from({length:40},(_,j)=>coverageExplorer.setParameters({perigee_altitude_km:2000+j*5}))));
  assert.equal(await page.evaluate(()=>coverageExplorer.lastGrid.values.perigee_altitude_km),2195);
  assert.equal(await page.evaluate(()=>coverageExplorer.values.perigee_altitude_km),2195);assert.deepEqual(errors,[]);
 }finally{await page.close();}
}
async function main(){
 const args=process.argv.slice(2);function option(name){const i=args.indexOf(name);if(i<0)return null;assert(args[i+1],name+' needs a value.');return args.splice(i,2)[1];}
 const shotsValue=option('--screenshot-dir'),libraryValue=option('--plotly-js');assert(args.length<=1,'Expected at most one report HTML path.');
 const report=args[0]?path.resolve(args[0]):null,shots=shotsValue?path.resolve(shotsValue):null;
 if(shots)fs.mkdirSync(shots,{recursive:true});const library=localLibrary(libraryValue,report),directory=fs.mkdtempSync(path.join(os.tmpdir(),'coverage-explorer-tests-'));
 let browser;
 try{
  const candidates=[process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE,chromium.executablePath(),'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe','C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe'].filter(Boolean),executablePath=candidates.find(p=>fs.existsSync(p));
  browser=await chromium.launch({headless:true,...(executablePath?{executablePath}:{})});
  const clean=syntheticData(),cleanFile=writeFixture(directory,'synthetic',clean,library);await checkReport(browser,cleanFile,clean,shots,'synthetic');
  const extended=syntheticData(true);await checkReport(browser,writeFixture(directory,'synthetic-extended',extended,library),extended,shots,'synthetic-extended');
  const holes=syntheticData(),statuses=['earth_intersection','mean_earth_intersection','no_phase_shear','unobserved'];
  for(let j=0;j<statuses.length;j++)Object.assign(holes.runs[j],{status:statuses[j],coverage_time_days:null,coverage_lower_days:null,coverage_upper_days:null});holes.runs.splice(4,1);
  await checkHoles(browser,writeFixture(directory,'synthetic-holes',holes,library),holes);
  if(report)await checkReport(browser,report,embeddedData(fs.readFileSync(report,'utf8')),shots,'actual');
  console.log('PASS: renderer UI; all sliders through 30000 km; measured values; independent linear/log-radius midpoint interpolation; original reference surface and CSV; fitted/shared height and color scales; all invalid/missing corner masks; camera preservation/reset; rapid updates; CSV download; desktop/tablet/mobile; offline operation'+(report?'; actual report':'')+'.');
  if(shots)console.log('Screenshots: '+shots);
 }finally{if(browser)await browser.close();const target=fs.realpathSync(directory),temporaryRoot=fs.realpathSync(os.tmpdir());assert.equal(path.dirname(target),temporaryRoot,'Cleanup target must stay inside the temporary directory.');assert(path.basename(target).startsWith('coverage-explorer-tests-'));fs.rmSync(target,{recursive:true,force:true});}
}
main().catch(error=>{console.error(error);process.exitCode=1;});
