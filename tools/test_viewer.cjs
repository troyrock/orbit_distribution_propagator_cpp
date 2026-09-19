#!/usr/bin/env node
'use strict';
// Usage: node tools/test_viewer.cjs [report.html] [--screenshot-dir directory]
// Requires Playwright + Chromium. With no report, verifies a deterministic synthetic run.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
let chromium;
try { ({ chromium } = require('playwright')); }
catch { ({ chromium } = require(path.resolve(path.dirname(process.execPath), '../node_modules/playwright'))); }

function syntheticData() {
  const count=600, epochs=16, a=2.656e7, inc=55*Math.PI/180;
  const point=(angle,r,z=0)=>[r*Math.cos(angle),r*Math.sin(angle)*Math.cos(inc)-z*Math.sin(inc),r*Math.sin(angle)*Math.sin(inc)+z*Math.cos(inc)];
  const reference=Array.from({length:181},(_,i)=>point(i/180*2*Math.PI,a));
  const frames=Array.from({length:epochs},(_,j)=>{
    const width=.00001+j*j/18;
    const angles=Array.from({length:count},(_,i)=>(i/(count-1)-.5)*width);
    const histogram=Array(72).fill(0);
    for (const p of angles) histogram[Math.min(71,Math.floor(((p%(2*Math.PI)+2*Math.PI)%(2*Math.PI))/(2*Math.PI)*72))]++;
    const resultants=[1,2,3,4].map(k=>Math.hypot(angles.reduce((s,p)=>s+Math.cos(k*p),0),angles.reduce((s,p)=>s+Math.sin(k*p),0))/count);
    return {time_s:j*30*86400,reference_elements:[a,0,0,Math.tan(inc/2),0,j*.3],positions_m:angles.map((p,i)=>point(p+j*.3,a+150*Math.sin(i*7.31),100*Math.cos(i*5.17))),reference_orbit_m:reference,
      metrics:{phase_sigma_rad:width/Math.sqrt(12),phase_q95_width_rad:width*.95,resultants,max_gap_deg:Math.max(.6,360-width*180/Math.PI),occupied_fraction:histogram.filter(x=>x>0).length/72,histogram,rtn_sigma_m:[120,Math.max(150,width*a/3),70],coverage:j>=11,mixed:j>=14}};
  });
  return {schema_version:1,metadata:{samples:count,visual_samples:count,force_model:'Synthetic viewer test (not a propagation result)',input_type:'mean',output_type:'osculating',seed:Number('9007199254740997'),seed_string:'9007199254740997',mu:3.986004418e14,earth_radius_m:6378137,elapsed_seconds:1.23,phase_definition:'Relative mean longitude modulo 360°. Counts use all samples.',coverage_definition:'Synthetic coverage flag for UI testing.',mixing_definition:'Synthetic mixing flag for UI testing.'},frames,summary:{coverage_time_s:frames[11].time_s,mixing_time_s:frames[14].time_s,analytic_mixing_time_s:350*86400}};
}

// Independent classical-element oracle: do not reuse the viewer's equinoctial
// conversion or its section-selection functions when checking the geometry.
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const norm=a=>Math.hypot(...a);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const unit=a=>a.map(x=>x/norm(a));
const radians=degrees=>degrees*Math.PI/180;
const wrappedDifference=(a,b)=>((a-b+540)%360)-180;
function classicalState(angleDeg,eccentricity=.6) {
  const a=2.656e7,mu=3.986004418e14,inc=.8,node=eccentricity?.45:0,arg=eccentricity?.7:0,f=radians(angleDeg);
  const rotate=v=>{
    const x=Math.cos(arg)*v[0]-Math.sin(arg)*v[1],y=Math.sin(arg)*v[0]+Math.cos(arg)*v[1];
    return [Math.cos(node)*x-Math.sin(node)*Math.cos(inc)*y,
      Math.sin(node)*x+Math.cos(node)*Math.cos(inc)*y,Math.sin(inc)*y];
  };
  const p=a*(1-eccentricity**2),r=p/(1+eccentricity*Math.cos(f)),speed=Math.sqrt(mu/p);
  const anchor=rotate([r*Math.cos(f),r*Math.sin(f)]),velocity=rotate([-speed*Math.sin(f),speed*(eccentricity+Math.cos(f))]);
  const normal=unit(cross(anchor,velocity)),xAxis=unit(cross(velocity,normal)),travel=unit(velocity);
  const E=2*Math.atan2(Math.sqrt(1-eccentricity)*Math.sin(f/2),Math.sqrt(1+eccentricity)*Math.cos(f/2));
  const M=E-eccentricity*Math.sin(E);
  return {anchor,velocity,normal,xAxis,travel,elements:[a,eccentricity*Math.cos(node+arg),eccentricity*Math.sin(node+arg),Math.tan(inc/2)*Math.cos(node),Math.tan(inc/2)*Math.sin(node),M+node+arg]};
}
function sectionFixture(eccentricity=.6) {
  const states=[0,30,60].map(angle=>classicalState(angle,eccentricity));
  const sixty=classicalState(60,eccentricity),zero=classicalState(0,eccentricity);
  const offset=(state,x,y,depth)=>state.anchor.map((v,j)=>v+x*state.xAxis[j]+y*state.normal[j]+depth*state.travel[j]);
  const positions=[offset(sixty,120,-80,23),offset(sixty,-240,160,-31),offset(sixty,0,0,0),
    classicalState(240,eccentricity).anchor,classicalState(180,eccentricity).anchor,
    offset(zero,50,-90,0),classicalState(359,eccentricity).anchor,classicalState(1,eccentricity).anchor,
    classicalState(4,eccentricity).anchor];
  const reference=Array.from({length:181},(_,i)=>classicalState(i*2,eccentricity).anchor);
  const data=syntheticData();
  data.metadata.samples=positions.length;data.metadata.visual_samples=positions.length;
  data.frames=states.map((state,index)=>({time_s:index*86400,reference_elements:state.elements,
    positions_m:positions,reference_orbit_m:reference,metrics:data.frames[0].metrics}));
  return {data,states,sixty};
}
function near(actual,expected,tolerance,label) {
  assert(Number.isFinite(actual),label+' must be finite.');
  assert(Math.abs(actual-expected)<=tolerance,`${label}: ${actual} versus ${expected}, tolerance ${tolerance}`);
}
function nearVector(actual,expected,tolerance,label) {
  assert.equal(actual?.length,3,label+' must have three components.');
  actual.forEach((value,index)=>near(value,expected[index],tolerance,`${label}[${index}]`));
}
function checkSectionFrame(section,expected) {
  assert.equal(section.available,true,'A valid orbital reference must enable the cross-section.');
  assert.equal(section.source,'elements','Exported elements must provide the exact reference geometry.');
  nearVector(section.anchor,expected.anchor,2e-7,'Anchor position [m]');
  nearVector(section.velocity,expected.velocity,2e-9,'Reference velocity [m/s]');
  nearVector(section.normal,expected.normal,2e-14,'Orbit normal');
  nearVector(section.xAxis,expected.xAxis,2e-14,'In-plane perpendicular axis');
  near(dot(section.xAxis,unit(section.velocity)),0,2e-14,'Cross-section perpendicular to true velocity');
  near(dot(section.normal,unit(section.velocity)),0,2e-14,'Normal perpendicular to true velocity');
}
async function setRange(page,id,value) {
  await page.locator('#'+id).fill(String(value));
  // Wait for the same animation-frame rendering path used by real interaction.
  await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
}
async function checkSectionGeometry(browser,template,screenshotDir) {
  const directory=fs.mkdtempSync(path.join(os.tmpdir(),'distribution-section-'));
  const page=await browser.newPage({viewport:{width:1440,height:1120},deviceScaleFactor:1});
  const errors=[],network=[];
  page.on('pageerror',e=>errors.push(e.message));
  page.on('request',req=>{if(!['file:','data:','about:'].includes(new URL(req.url()).protocol))network.push(req.url());});
  const load=async(data,name)=>{
    const file=path.join(directory,name+'.html');
    fs.writeFileSync(file,template.replace('__DISTRIBUTION_DATA__',JSON.stringify(data)));
    await page.goto(pathToFileURL(file).href);
    await page.waitForFunction(()=>window.distributionViewer?.ready===true&&window.distributionViewer.section);
  };
  const section=async()=>{
    await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
    return page.evaluate(()=>window.distributionViewer.section);
  };
  try {
    for(const eccentricity of [0,.6]) {
      const fixture=sectionFixture(eccentricity);
      await load(fixture.data,'elements-'+eccentricity);
      assert.equal(await page.locator('#section-follow').isChecked(),true,'Follow mode is the default.');
      assert.equal(await page.locator('#section-curvature').isChecked(),true,'Curvature correction is the default.');
      near(wrappedDifference((await section()).angleDeg,0),0,1e-10,'Initial nominal true anomaly');
      await page.click('#next');
      near(wrappedDifference((await section()).angleDeg,30),0,1e-10,'Following advances with the saved nominal state');
      await setRange(page,'section-angle',60);
      await page.locator('#section-curvature').uncheck();
      assert.equal(await page.locator('#section-follow').isChecked(),false,'Choosing an angle switches to manual placement.');
      let state=await section();
      assert.equal(state.curvatureCorrected,false,'Raw projection must be explicitly selectable.');
      checkSectionFrame(state,fixture.sixty);
      assert.equal(state.following,false);assert.equal(state.displayedCount,9);
      assert.equal(state.selectedCount,3,'A local phase band must exclude the opposite arc.');
      assert.deepEqual(state.points.map(p=>p.id).sort((a,b)=>a-b),[0,1,2]);
      const expected=[[120,-80,23],[-240,160,-31],[0,0,0]];
      state.points.forEach(p=>{near(p.x,expected[p.id][0],3e-7,'Projected in-plane offset [m]');near(p.y,expected[p.id][1],3e-7,'Projected normal offset [m]');near(p.depth,expected[p.id][2],3e-7,'Travel-direction depth [m]');});
      if(eccentricity)assert(Math.abs(dot(unit(state.anchor),unit(state.velocity)))>.3,'Fixture must distinguish velocity-normal from radial/RTN geometry.');
      await page.click('#next');
      state=await section();near(state.angleDeg,60,1e-10,'Manual angle survives epoch changes');
      assert.equal(state.following,false);
      await setRange(page,'section-angle',90);
      state=await section();assert.equal(state.available,true);assert.equal(state.selectedCount,0);assert.deepEqual(state.points,[]);
      assert([...state.anchor,...state.velocity,...state.normal,...state.xAxis].every(Number.isFinite),'Empty selections retain finite reference axes.');
      assert(!/NaN|Infinity/.test(await page.locator('body').innerText()),'Empty sections must not display invalid numeric labels.');
      const pixels=await page.locator('#section-scene').evaluate(c=>{const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0;for(let i=3;i<d.length;i+=4)if(d[i]>0)n++;return n;});
      assert(pixels>100,'An empty local band still renders axes and its empty-state message.');
      await setRange(page,'section-angle',0);const zero=await section();
      assert.deepEqual(zero.points.map(p=>p.id).sort((a,b)=>a-b),[5,6,7,8],'The phase window must wrap across zero.');
      await setRange(page,'section-angle',360);const full=await section();
      nearVector(full.anchor,zero.anchor,2e-7,'0 and 360 degree anchors');
      assert.deepEqual(full.points.map(p=>p.id),zero.points.map(p=>p.id),'0 and 360 select the same particles.');
      await setRange(page,'section-width',2);state=await section();
      assert(!state.points.some(p=>p.id===8),'Narrowing the angular band excludes the four-degree point.');
      assert(state.points.some(p=>p.id===5),'Narrowing preserves the central point.');
      await setRange(page,'section-width',10);
      await page.locator('#section-follow').check();
      near(wrappedDifference((await section()).angleDeg,60),0,1e-10,'Re-enabling follow selects the current nominal anomaly');
      await setRange(page,'time',0);await page.selectOption('#rate','8');await page.click('#play');
      await page.waitForFunction(()=>window.distributionViewer.index===2&&!window.distributionViewer.playing);
      near(wrappedDifference((await section()).angleDeg,60),0,1e-10,'Follow remains synchronized throughout saved-epoch playback');
      if(screenshotDir&&eccentricity)await page.screenshot({path:path.join(screenshotDir,'desktop-cross-section.png'),fullPage:true});
      await page.setViewportSize({width:390,height:844});
      await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
      assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),'Cross-section controls fit mobile width.');
      assert(await page.locator('#section-scene').isVisible(),'Cross-section remains visible on mobile.');
      if(screenshotDir&&eccentricity)await page.screenshot({path:path.join(screenshotDir,'mobile-cross-section.png'),fullPage:true});
      await page.setViewportSize({width:1440,height:1120});
    }
    // Below the circular-label threshold, the angle origin changes to the
    // equinoctial axis. The small physical eccentricity and its periapsis must
    // still be retained in both the anchor radius and nominal follow position.
    const nearCircular=sectionFixture(5e-13);
    await load(nearCircular.data,'near-circular');
    checkSectionFrame(await section(),nearCircular.states[0]);
    await page.click('#next');
    checkSectionFrame(await section(),nearCircular.states[1]);
    await setRange(page,'section-angle',120);
    const e=nearCircular.states[0].elements;
    const trueAnomaly=120-Math.atan2(e[2],e[1])*180/Math.PI;
    checkSectionFrame(await section(),classicalState(trueAnomaly,5e-13));
    const ellipse=sectionFixture().data;
    const arc=Array.from({length:41},(_,i)=>classicalState(60-18+i*.9).anchor);
    ellipse.frames.forEach(f=>{f.positions_m=arc;});
    ellipse.metadata.samples=arc.length;ellipse.metadata.visual_samples=arc.length;
    await load(ellipse,'perfect-ellipse');
    await setRange(page,'section-angle',60);await setRange(page,'section-width',40);
    let perfect=await section();assert.equal(perfect.curvatureCorrected,true);assert.equal(perfect.selectedCount,arc.length);
    perfect.points.forEach(p=>{near(p.x,0,3e-7,'Curvature-corrected ellipse in-plane residual [m]');near(p.y,0,3e-7,'Curvature-corrected ellipse normal residual [m]');});
    await page.locator('#section-curvature').uncheck();perfect=await section();
    assert(Math.max(...perfect.points.map(p=>Math.abs(p.x)))>10000,'Raw ellipse projection must retain finite-aperture orbital curvature.');
    const legacy=sectionFixture().data;legacy.frames.forEach(f=>delete f.reference_elements);
    await load(legacy,'legacy-polyline');
    const fallback=await section();assert.equal(fallback.available,false);assert.equal(fallback.source,'unavailable');assert.deepEqual(fallback.points,[]);
    assert(/regenerat/i.test(await page.locator('body').innerText()),'Older reports must explain that reference states require regeneration.');
    assert.deepEqual(network,[],'Cross-section must remain completely offline.');
    assert.deepEqual(errors,[],'Cross-section must not produce browser errors.');
  } finally {await page.close();fs.rmSync(directory,{recursive:true,force:true});}
}

async function main() {
  const args=process.argv.slice(2), shotIndex=args.indexOf('--screenshot-dir');
  const screenshotDir=shotIndex>=0?path.resolve(args[shotIndex+1]):null;
  if(shotIndex>=0)args.splice(shotIndex,2);
  assert(args.length<=1,'Expected at most one report HTML path.');
  const temporary=args.length===0?fs.mkdtempSync(path.join(os.tmpdir(),'distribution-viewer-')):null;
  let html=args.length?path.resolve(args[0]):path.join(temporary,'synthetic.html');
  const template=fs.readFileSync(path.join(__dirname,'../web/viewer.html'),'utf8');
  assert.equal(template.split('__DISTRIBUTION_DATA__').length,2,'Template must have one data insertion token.');
  if(temporary)fs.writeFileSync(html,template.replace('__DISTRIBUTION_DATA__',JSON.stringify(syntheticData())));
  if(screenshotDir)fs.mkdirSync(screenshotDir,{recursive:true});
  let browser;
  try {
    const candidates=[process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE,chromium.executablePath(),
      'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
      'C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe'].filter(Boolean);
    const executablePath=candidates.find(candidate=>fs.existsSync(candidate));
    browser=await chromium.launch({headless:true,...(executablePath?{executablePath}:{})});
    const page=await browser.newPage({viewport:{width:1440,height:1120},deviceScaleFactor:1});
    const errors=[],network=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('request',req=>{if(!['file:','data:','about:'].includes(new URL(req.url()).protocol))network.push(req.url());});
    await page.goto(pathToFileURL(html).href);
    await page.waitForFunction(()=>window.distributionViewer?.ready===true);
    await page.waitForFunction(()=>Number(document.getElementById('scene').dataset.renderedPoints)>0);
    assert.equal(await page.locator('#error').isVisible(),false);
    if(temporary)assert((await page.locator('#normalization').textContent()).includes('9007199254740997'),
      'The exact uint64 seed must survive JavaScript number rounding.');
    const total=await page.evaluate(()=>JSON.parse(document.getElementById('simulation-data').textContent).frames.length);
    const capture=async name=>{if(screenshotDir)await page.screenshot({path:path.join(screenshotDir,name+'.png'),fullPage:true});};
    const noOverflow=async()=>assert(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth+1),'Layout must fit viewport width.');
    await noOverflow();
    await capture('desktop-initial');
    await page.selectOption('#framing','cloud');
    await page.waitForTimeout(80);
    const cloudDrawn=await page.locator('#scene').getAttribute('data-rendered-points');
    assert(Number(cloudDrawn)>0,'Focus cloud must render sample points.');
    await capture('desktop-cloud');
    const bounds=await page.locator('#scene').boundingBox(),mid={x:bounds.x+bounds.width/2,y:bounds.y+bounds.height/2};
    await page.mouse.move(mid.x,mid.y);await page.mouse.wheel(0,-450);
    await page.waitForFunction(()=>window.distributionViewer.zoom>1.1);
    const anglesBefore=await page.evaluate(()=>window.distributionViewer.angles);
    await page.mouse.move(mid.x,mid.y);await page.mouse.down();await page.mouse.move(mid.x+70,mid.y+35,{steps:5});await page.mouse.up();
    const anglesAfter=await page.evaluate(()=>window.distributionViewer.angles);
    assert.notDeepEqual(anglesBefore,anglesAfter,'Drag must rotate the camera.');
    await page.keyboard.down('Shift');await page.mouse.down();await page.mouse.move(mid.x+110,mid.y+70,{steps:5});await page.mouse.up();await page.keyboard.up('Shift');
    assert((await page.evaluate(()=>window.distributionViewer.pan)).some(Math.abs),'Shift-drag must pan.');
    await page.click('#reset');assert.equal(await page.evaluate(()=>window.distributionViewer.zoom),1);
    await page.selectOption('#framing','orbit');await page.selectOption('#view','plane');
    // Prefer an epoch near a quarter-orbit phase width to show the actual early arc.
    // This only selects a saved snapshot; it does not manufacture intermediate positions.
    const earlyIndex=await page.evaluate(()=>{
      const frames=JSON.parse(document.getElementById('simulation-data').textContent).frames;
      let best=0,distance=Infinity;
      frames.forEach((f,i)=>{const d=Math.abs(f.metrics.phase_q95_width_rad-Math.PI/2);if(d<distance){best=i;distance=d;}});
      return best;
    });
    await page.locator('#time').fill(String(earlyIndex));
    await page.waitForFunction(i=>Number(document.getElementById('scene').dataset.epochIndex)===i,earlyIndex);
    await capture('desktop-early-banana');
    await page.locator('#time').fill(String(total-1));
    await page.waitForFunction(i=>window.distributionViewer.index===i,total-1);
    assert.equal(await page.locator('#next').isDisabled(),true);
    await capture('desktop-final-plane');
    if(total>1){
      await page.click('#previous');assert.equal(await page.evaluate(()=>window.distributionViewer.index),total-2);
      await page.locator('#time').fill('0');await page.selectOption('#rate','8');await page.click('#play');
      await page.waitForFunction(()=>window.distributionViewer.index>0);
      await page.evaluate(()=>window.distributionViewer.stop());
      assert.equal(await page.evaluate(()=>window.distributionViewer.playing),false);
      await page.locator('#time').fill(String(total-2));await page.click('#play');
      await page.waitForFunction(i=>window.distributionViewer.index===i&&!window.distributionViewer.playing,total-1);
    }
    await page.selectOption('#view','edge');
    await page.locator('#reference').uncheck();await page.locator('#hidden-particles').uncheck();
    await page.selectOption('#view','plane');await page.locator('#reference').check();await page.locator('#hidden-particles').check();
    await page.setViewportSize({width:390,height:844});await page.waitForTimeout(100);await noOverflow();await capture('mobile-final');
    const pixels=await page.locator('#scene').evaluate(c=>{const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0;for(let i=3;i<d.length;i+=4)if(d[i]>0)n++;return n;});
    assert(pixels>100,'Scene canvas must contain rendered geometry.');
    assert.deepEqual(network,[],'Offline artifact must make no network requests.');
    assert.deepEqual(errors,[],'Viewer must not produce browser errors.');
    await checkSectionGeometry(browser,template,screenshotDir);
    console.log('PASS: '+total+' epochs; desktop/mobile layout; cloud/orbit framing; rotation/pan/zoom; saved-epoch playback; circular/eccentric cross-section geometry, curvature correction, angular wrap, local selection, empty bins and follow controls; legacy-report handling; offline canvas rendering.');
    if(screenshotDir)console.log('Screenshots: '+screenshotDir);
  } finally {if(browser)await browser.close();if(temporary)fs.rmSync(temporary,{recursive:true,force:true});}
}
main().catch(error=>{console.error(error);process.exitCode=1;});
