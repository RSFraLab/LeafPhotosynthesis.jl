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
    await page.waitForFunction(()=>document.querySelectorAll('canvas').length===25);
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
      assert(await page.evaluate(precision=>{
        const jb=window.leafExamples.data[precision].johnson_berry;
        const c9=jb.jb_connectivity, c10=jb.jb_state_transitions, c12=jb.jb_capacity, pam=jb.jb_pam;
        const specs=[
          ['jb-connectivity',c9,'incident_par',['ku2_psii_yield','ku10_psii_yield','puddle_psii_yield']],
          ['jb-state-transitions',c10,'incident_par',['qt_ku2_yield','static_ku2_yield','qt_ku10_yield','static_ku10_yield']],
          ['jb-capacity',c12,'incident_par',['cap_10x_yield','cap_05x_yield','cap_15x_yield','cap_25x_yield']],
          ['jb-pam',pam,'incident_par',['true_psii_photochemical_yield','pam_apparent_photochemical_yield','ideal_closed_lake_1_minus_ql','pam_closed_lake_1_minus_ql']]
        ];
        // Verify Case 12 low-light vs high-light ΦF2 capacity inversion and PAM PSI-leakage invariance of qP:
        const lowIdx=10, highIdx=100;
        const capInversion = c12[lowIdx].cap_05x_yield > c12[lowIdx].cap_25x_yield &&
          c12[highIdx].cap_05x_yield < c12[highIdx].cap_25x_yield;
        const pamLeakCheck = pam[lowIdx].pam_apparent_photochemical_yield < pam[lowIdx].true_psii_photochemical_yield &&
          pam[lowIdx].pam_closed_lake_1_minus_ql < pam[lowIdx].ideal_closed_lake_1_minus_ql &&
          Math.abs(pam[lowIdx].pam_closed_puddle_1_minus_qp - pam[lowIdx].ideal_closed_puddle_1_minus_qp) < 1e-5;
        return c9[0].dark && c9[0].ku2_psii_yield===null && capInversion && pamLeakCheck &&
          specs.every(([id,rows,xKey,keys])=>{
            const chart=Chart.getChart(id);
            return chart.data.datasets.length===keys.length && chart.data.datasets.every((s,k)=>s.data.length===rows.length &&
              s.data.every((p,i)=>p.x===rows[i][xKey]&&p.y===rows[i][keys[k]]) &&
              chart.getDatasetMeta(k).data.every((p,i)=>s.data[i].y===null?p.skip:Number.isFinite(p.x)&&Number.isFinite(p.y)));
          });
      },precision));
      assert(await page.evaluate(precision=>{
        const fx=window.leafExamples.data[precision].flexas_2002;
        const fig2=fx.flexas_mechanistic_fig2, emp=fx.flexas_empirical_fig2, inset=fx.flexas_inset_conductance, ctrl=fx.flexas_controls, dig=fx.flexas_digitized_fig2;
        const cFig2=Chart.getChart('flexas-mechanistic-fig2');
        const cEmp=Chart.getChart('flexas-empirical-fig2');
        const cIns=Chart.getChart('flexas-inset-conductance');
        const cCtrl=Chart.getChart('flexas-controls');
        const peak150=Math.max(...fig2.map(r=>r.g150_fs_over_fo));
        const peak10=Math.max(...fig2.map(r=>r.g10_fs_over_fo));
        const empPeakDiff=Math.abs(Math.max(...emp.slice(1).map(r=>r.cotton_g150))-Math.max(...emp.slice(1).map(r=>r.cotton_g10)));
        const obsPointCount=cFig2.data.datasets.slice(6).reduce((acc,ds)=>acc+ds.data.length,0);
        return cFig2.data.datasets.length===12 && obsPointCount===dig.length && dig.length===55 &&
          cEmp.data.datasets.length===4 &&
          cIns.data.datasets.length===6 && cIns.data.datasets[5].data.length===6 &&
          cCtrl.data.datasets.length===4 &&
          fig2[0].g150_fs_over_fo===1 && peak150>peak10 && fig2[fig2.length-1].g10_fs_over_fo<1 &&
          empPeakDiff<0.02 && emp[emp.length-1].cotton_g10>emp[emp.length-1].cotton_g150 &&
          inset[0].pam2000_combined<inset[inset.length-1].pam2000_combined &&
          inset[0].empirical_cotton>inset[inset.length-1].empirical_cotton &&
          ctrl[ctrl.length-1].high_co2_reversion_pam>1.3 && ctrl[ctrl.length-1].dtt_unquenched>2.3;
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

    // Verify standalone Flexas et al. (2002) interactive explorer (docs/flexas_2002.html)
    const fxPage=await browser.newPage({viewport:{width:1280,height:960}});
    fxPage.on('pageerror', e=>errors.push(e.message));
    fxPage.on('request', r=>{if(/^https?:/.test(r.url()))requests.push(r.url());});
    await fxPage.goto('file://'+path.join(__dirname,'flexas_2002.html'));
    await fxPage.waitForFunction(()=>document.querySelectorAll('canvas').length===7 && document.getElementById('status').textContent.includes('Julia output'));
    assert(await fxPage.evaluate(()=>document.querySelectorAll('.katex').length>=6 && document.querySelectorAll('.katex-error').length===0));
    for(const precision of ['Float64','Float32']) {
      await fxPage.selectOption('#precision',precision);
      for(const form of ['pam2000_combined','pam2000_stomatal_only','pure_psii_combined','qt_pam2000_combined','cotton_combined','drought_combined']) {
        await fxPage.locator(`button[data-formulation="${form}"]`).click();
        assert(await fxPage.evaluate(form=>{
          const c=Chart.getChart('hero-canvas');
          const obsCount=c.data.datasets.slice(6).reduce((a,ds)=>a+ds.data.length,0);
          return c.data.datasets.length===12 && obsCount===55 &&
            document.querySelector(`button[data-formulation="${form}"]`).getAttribute('aria-pressed')==='true';
        },form));
      }
      await fxPage.locator('button[data-formulation="pam2000_combined"]').click();
      for(const g of ['150','80','62','55','39','10']) {
        await fxPage.locator(`button[data-g="${g}"]`).click();
        assert(await fxPage.evaluate(g=>{
          const m=Chart.getChart('spotlight-models-canvas');
          const intern=Chart.getChart('spotlight-internal-canvas');
          return m.data.datasets.length===7 && intern.data.datasets.length===3 &&
            document.querySelector(`button[data-g="${g}"]`).getAttribute('aria-pressed')==='true';
        },g));
      }
    }
    await fxPage.locator('#internal-mode-quenching').click();
    assert(await fxPage.evaluate(()=>Boolean(Chart.getChart('spotlight-internal-canvas').scales.y1)));
    await fxPage.locator('#domain-extended').click();
    assert(await fxPage.evaluate(()=>Chart.getChart('hero-canvas').scales.x.max===2100));
    await fxPage.locator('#domain-paper').click();
    assert(await fxPage.evaluate(()=>Chart.getChart('hero-canvas').scales.x.max===1400));
    await fxPage.locator('#hero-table-details summary').click();
    await fxPage.locator('#hero-table-mount tbody tr').first().waitFor();
    assert.strictEqual(await fxPage.locator('#hero-table-mount tbody tr').count(),55);
    const fxDownload=fxPage.waitForEvent('download');
    await fxPage.locator('#hero-actions').getByRole('button',{name:'Save PNG'}).click();
    assert((await fxDownload).suggestedFilename().includes('flexas-fig2-pam2000_combined-Float32'));
    await fxPage.setViewportSize({width:390,height:844});
    await fxPage.waitForTimeout(300);
    assert(await fxPage.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));

    assert.deepStrictEqual(errors,[]);assert.deepStrictEqual(requests,[]);
    console.log('PASS: examples.html (25 plots + 55 digitized Fig. 2 points) & flexas_2002.html (7 interactive canvases, 6 SIF formulations, 6 spotlight conductances, offline KaTeX, FP32/64, PNG export, mobile & offline)');
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exitCode=1;});
