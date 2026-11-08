const { chromium } = require('playwright');
const path = require('path');
const assert = require('assert');

(async () => {
  const browser = await chromium.launch({headless:true, executablePath:process.env.LEAF_DOCS_CHROMIUM});
  try {
    const page = await browser.newPage({viewport:{width:1280,height:960}});
    const errors=[], requests=[];
    page.on('pageerror', e=>errors.push(e.message));
    page.on('request', r=>{if(/^https?:/.test(r.url()))requests.push(r.url());});
    await page.goto('file://'+path.join(__dirname,'examples.html'));
    await page.waitForFunction(()=>document.querySelectorAll('canvas').length===9);
    for(const precision of ['Float64','Float32']) {
      await page.selectOption('#precision',precision);
      for(const index of ['0','1','2']) {
        await page.selectOption('#oxygen-temperature',index);
        const checked=await page.evaluate(({precision,index})=>{
          const all=window.leafExamples.data[precision].oxygen_response;
          const temperatures=[...new Set(all.map(r=>r.temperature))];
          const rows=all.filter(r=>r.temperature===temperatures[Number(index)]);
          return ['oxygen-gamma','oxygen-carbon'].every((id,plotIndex)=>{
            const c=Chart.getChart(id), keys=plotIndex===0?['scaled_gamma','prescribed_gamma']:['scaled_assimilation','prescribed_assimilation'];
            return c.scales.x.type==='linear' && c.data.datasets.every((s,j)=>s.data.length===rows.length &&
              s.data.every((p,i)=>p.x===rows[i].oxygen&&p.y===rows[i][keys[j]])) &&
              c.getDatasetMeta(0).data.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y));
          }) && Array.from(document.querySelectorAll('canvas')).every(canvas=>{
            const c=Chart.getChart(canvas);return c.data.datasets.every(s=>s.data.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)));
          });
        },{precision,index});
        assert(checked);
      }
      const card=page.locator('article').filter({has:page.locator('#oxygen-carbon')});
      await card.locator('details summary').click();
      await card.locator('tbody tr').first().waitFor();
      assert.strictEqual(await card.locator('tbody tr').count(),5);
      await card.locator('input[type=checkbox]').first().uncheck();
      assert(await page.evaluate(()=>!Chart.getChart('oxygen-carbon').isDatasetVisible(0)));
      await card.locator('input[type=checkbox]').first().check();
    }
    const card=page.locator('article').filter({has:page.locator('#oxygen-carbon')});
    const download=page.waitForEvent('download');
    await card.getByRole('button',{name:'Save PNG'}).click();
    assert((await download).suggestedFilename().includes('oxygen-carbon-Float32'));
    await page.setViewportSize({width:390,height:844});
    await page.waitForTimeout(300);
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));
    assert.deepStrictEqual(errors,[]);assert.deepStrictEqual(requests,[]);
    console.log('PASS: nine plots, FP32/64, three oxygen temperatures, exact data, rendered points, tables, series toggles, PNG download, mobile and offline');
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exitCode=1;});
