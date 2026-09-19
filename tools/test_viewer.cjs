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
    return {time_s:j*30*86400,positions_m:angles.map((p,i)=>point(p+j*.3,a+150*Math.sin(i*7.31),100*Math.cos(i*5.17))),reference_orbit_m:reference,
      metrics:{phase_sigma_rad:width/Math.sqrt(12),phase_q95_width_rad:width*.95,resultants,max_gap_deg:Math.max(.6,360-width*180/Math.PI),occupied_fraction:histogram.filter(x=>x>0).length/72,histogram,rtn_sigma_m:[120,Math.max(150,width*a/3),70],coverage:j>=11,mixed:j>=14}};
  });
  return {schema_version:1,metadata:{samples:count,visual_samples:count,force_model:'Synthetic viewer test (not a propagation result)',input_type:'mean',output_type:'osculating',seed:42,mu:3.986004418e14,earth_radius_m:6378137,elapsed_seconds:1.23,phase_definition:'Relative mean longitude modulo 360°. Counts use all samples.',coverage_definition:'Synthetic coverage flag for UI testing.',mixing_definition:'Synthetic mixing flag for UI testing.'},frames,summary:{coverage_time_s:frames[11].time_s,mixing_time_s:frames[14].time_s,analytic_mixing_time_s:350*86400}};
}

async function main() {
  const args=process.argv.slice(2), shotIndex=args.indexOf('--screenshot-dir');
  const screenshotDir=shotIndex>=0?path.resolve(args[shotIndex+1]):null;
  if(shotIndex>=0)args.splice(shotIndex,2);
  assert(args.length<=1,'Expected at most one report HTML path.');
  const temporary=args.length===0?fs.mkdtempSync(path.join(os.tmpdir(),'distribution-viewer-')):null;
  let html=args.length?path.resolve(args[0]):path.join(temporary,'synthetic.html');
  if(temporary){const template=fs.readFileSync(path.join(__dirname,'../web/viewer.html'),'utf8');assert.equal(template.split('__DISTRIBUTION_DATA__').length,2,'Template must have one data insertion token.');fs.writeFileSync(html,template.replace('__DISTRIBUTION_DATA__',JSON.stringify(syntheticData())));}
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
    console.log('PASS: '+total+' epochs; desktop/mobile layout; cloud/orbit framing; rotation/pan/zoom; saved-epoch playback; offline canvas rendering.');
    if(screenshotDir)console.log('Screenshots: '+screenshotDir);
  } finally {if(browser)await browser.close();if(temporary)fs.rmSync(temporary,{recursive:true,force:true});}
}
main().catch(error=>{console.error(error);process.exitCode=1;});
