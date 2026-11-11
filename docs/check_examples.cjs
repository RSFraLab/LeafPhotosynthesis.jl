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
    await page.waitForFunction(()=>document.querySelectorAll('canvas').length===17);
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
            const c=Chart.getChart(canvas);return c.data.datasets.every(s=>s.data.every(p=>Number.isFinite(p.x)&&(Number.isFinite(p.y)||p.x===0&&p.y===null)));
          });
        },{precision,index});
        assert(checked);
      }
      assert(await page.evaluate(precision=>{
        const rows=window.leafExamples.data[precision].fluorescence;
        const specs=[['fluorescence-carbon',['net_assimilation']],
          ['fluorescence-photochemistry',['photochemical_yield']],
          ['fluorescence-yield',['unquenched_yield','regulated_yield']],
          ['fluorescence-photons',['unquenched_photons','regulated_photons']]];
        return rows[0].dark && rows[0].regulated_photons===0 && rows[0].regulated_yield===null &&
          specs.every(([id,keys])=>{
            const chart=Chart.getChart(id);
            return chart.data.datasets.every((s,k)=>s.data.length===rows.length &&
              s.data.every((p,i)=>p.x===rows[i].absorbed_par&&p.y===rows[i][keys[k]]) &&
              chart.getDatasetMeta(k).data.every((p,i)=>s.data[i].y===null?p.skip:Number.isFinite(p.x)&&Number.isFinite(p.y)));
          });
      },precision));
      assert(await page.evaluate(precision=>{
        const cs=window.leafExamples.data[precision].cytochrome_stress;
        const light=cs.cytochrome_light, drought=cs.cytochrome_drought, cap=cs.cytochrome_capacity;
        const specs=[
          ['cytochrome-light',light,'absorbed_par',['jb_psii_yield','jb_total_yield','vdt_cotton_yield','vdt_drought_yield']],
          ['cytochrome-drought-high',drought,'biochemical_co2',['high_jb_psii_yield','high_jb_total_yield','high_vdt_cotton_yield','high_vdt_drought_yield']],
          ['cytochrome-drought-moderate',drought,'biochemical_co2',['mod_jb_psii_yield','mod_jb_total_yield','mod_vdt_cotton_yield','mod_vdt_drought_yield']],
          ['cytochrome-capacity',cap,'capacity_fraction',['jb_rubisco_psii_yield','jb_b6f_psii_yield','vdt_rubisco_yield','vdt_jmax_yield']]
        ];
        const droughtSignCheck = drought[0].high_jb_psii_yield < drought[drought.length-1].high_jb_psii_yield &&
          drought[0].high_vdt_cotton_yield > drought[drought.length-1].high_vdt_cotton_yield;
        const capacityForkCheck = cap[0].jb_rubisco_psii_yield < cap[cap.length-1].jb_rubisco_psii_yield &&
          cap[0].jb_b6f_psii_yield > cap[cap.length-1].jb_b6f_psii_yield;
        return light[0].dark && light[0].jb_total_photons===0 && light[0].jb_psii_yield===null &&
          droughtSignCheck && capacityForkCheck &&
          specs.every(([id,rows,xKey,keys])=>{
            const chart=Chart.getChart(id);
            return chart.data.datasets.length===4 && chart.data.datasets.every((s,k)=>s.data.length===rows.length &&
              s.data.every((p,i)=>p.x===rows[i][xKey]&&p.y===rows[i][keys[k]]) &&
              chart.getDatasetMeta(k).data.every((p,i)=>s.data[i].y===null?p.skip:Number.isFinite(p.x)&&Number.isFinite(p.y)));
          });
      },precision));
      const fluorescenceCard=page.locator('article').filter({has:page.locator('#fluorescence-yield')});
      await fluorescenceCard.locator('details summary').click();
      await fluorescenceCard.locator('tbody tr').first().waitFor();
      assert((await fluorescenceCard.locator('tbody tr').first().innerText()).includes('Undefined (dark)'));
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
    console.log('PASS: seventeen plots, FP32/64, oxygen temperatures, exact fluorescence and Cytochrome stress data, dark gaps/zero source, rendered points, tables, series toggles, PNG download, mobile and offline');
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exitCode=1;});
