(function () {
  const css = () => getComputedStyle(document.documentElement);
  const tok = (n) => css().getPropertyValue(n).trim();
  const T = () => ({ ink: tok('--ink'), ink2: tok('--ink-2'), muted: tok('--muted'), grid: tok('--grid'), axis: tok('--axis'), surface: tok('--surface'),
    rule: tok('--rule'), ref: tok('--ref'), heat: tok('--heat'), accent: tok('--accent'),
    s: [tok('--s1'), tok('--s2'), tok('--s3'), tok('--s4'), tok('--s5'), tok('--s6'), tok('--s7')], q: [tok('--q1'), tok('--q2'), tok('--q3'), tok('--q4'), tok('--q5')],
    h: [tok('--h0'), tok('--h1'), tok('--h2'), tok('--h3'), tok('--h4'), tok('--h5')], mono: tok('--f-mono'), body: tok('--f-body') });
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = (v, d) => (v === null || v === undefined || !isFinite(v)) ? '–' : Number(v).toLocaleString('en-GB', { minimumFractionDigits: d, maximumFractionDigits: d }).replace('-', '−');
  const dec = (v) => { const a = Math.abs(v); return a >= 100 ? 0 : a >= 10 ? 1 : a >= 1 ? 2 : 3; };
  const LINE = 'path://M0,4 L16,4 L16,6 L0,6 Z', DASH = 'path://M0,4 L6,4 L6,6 L0,6 Z M10,4 L16,4 L16,6 L10,6 Z';
  const MARK = { diamond: 'path://M5,0 L10,5 L5,10 L0,5 Z', triangle: 'path://M5,0 L10,10 L0,10 Z', circle: 'path://M5,0 A5,5 0 1,1 4.99,0 Z' };
  const icon = (s) => s.pts ? MARK[s.symbol || 'circle'] : (s.dash ? DASH : LINE);
  const meas = document.createElement('canvas').getContext('2d');
  // rows the legend will take at this width (items wrap as in the chart library: icon 16 + gap 5 + text, 12 between items)
  function legendRows(names, width, t) { meas.font = `11.5px ${t.body}`; let rows = 1, x = 0;
    for (const n of names) { const w = 16 + 5 + meas.measureText(n).width * 1.04 + 2; if (x > 0 && x + w > width) { rows++; x = 0; } x += w + 12; }
    return rows; }
  const caption = (t, text, top) => ({ type: 'text', left: 8, top, silent: true, style: { text, fill: t.ink2, font: `12px ${t.body}` } });
  const charts = [], tables = {};
  const colour = (t, c, i) => c === 'ref' ? t.ref : c === 'ink' ? t.ink : c === 'heat' ? t.heat : c === 'muted' ? t.muted : (typeof c === 'string' && c[0] === 'q') ? t.q[+c[1]] : t.s[(c === undefined || c === null ? i : c) % 7];
  function axis(t, a, o) {
    a = a || {};
    const base = { type: a.log ? 'log' : (a.cat ? 'category' : 'value'), name: a.name || '', nameLocation: 'middle', nameGap: a.gap || 30, logBase: 10,
      nameTextStyle: { color: t.ink2, fontFamily: t.body, fontSize: 12 }, axisLine: { show: true, lineStyle: { color: t.axis } }, axisTick: { lineStyle: { color: t.axis } },
      axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: a.log ? (v) => String(v) : undefined }, splitLine: { show: !a.cat, lineStyle: { color: t.grid } }, scale: true };
    if (a.min !== undefined) base.min = a.min; if (a.max !== undefined) base.max = a.max; if (a.cat) base.data = a.cat; if (a.inverse) base.inverse = true;
    if (a.hideMax) base.axisLabel.showMaxLabel = false; if (a.interval) base.interval = a.interval;
    return Object.assign(base, o || {});
  }
  const tip = (t, o) => Object.assign({ backgroundColor: t.surface, borderColor: t.rule, borderWidth: 1, padding: [8, 10], textStyle: { color: t.ink, fontFamily: t.body, fontSize: 12 },
    extraCssText: 'box-shadow:0 2px 10px rgba(0,0,0,.12);border-radius:6px;' }, o || {});
  const legend = (t, o) => Object.assign({ top: 0, left: 0, itemWidth: 16, itemHeight: 8, itemGap: 12, textStyle: { color: t.ink2, fontFamily: t.body, fontSize: 11.5 }, icon: 'roundRect' }, o || {});

  function build(spec, t, W) {
    if (spec.kind === 'field') {
      const nx = spec.r.length, ny = spec.z.length;
      const data = []; for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data.push([i, j, spec.T[j][i]]);
      const small = innerWidth < 600;
      const lab = (pos) => ({ color: t.ink, fontFamily: t.body, fontSize: 10, formatter: (p) => p.name, position: pos, backgroundColor: t.surface, padding: [1, 4], borderRadius: 2 });
      const marks = (spec.marks || []).map(m => m.x !== undefined ? { xAxis: m.x, name: m.name, label: lab('insideEndTop') } : { yAxis: m.y, name: m.name, label: lab('insideStartTop') });
      return { grid: { left: 36, right: 14, top: 30, bottom: 100, containLabel: false }, graphic: [caption(t, 'Height above the base (cm)', 6)],
        tooltip: tip(t, { trigger: 'item', formatter: (p) => `<div style="color:${t.muted};font-size:11px">${esc(Math.abs(spec.r[p.value[0]]))} cm from the axis, ${esc(spec.z[p.value[1]])} cm above the base</div><b style="font-family:${t.mono}">${fmt(p.value[2], 1)} °C</b>` }),
        xAxis: axis(t, { cat: spec.r.map(String), name: 'Distance from the well axis (cm)', gap: 26 }, { splitLine: { show: false }, axisTick: { show: false },
          axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, interval: 0, formatter: (v) => (Number(v) % 20 === 0 ? String(Math.abs(Number(v))) : '') } }),
        yAxis: axis(t, { cat: spec.z.map(String) }, { splitLine: { show: false }, axisTick: { show: false },
          axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, interval: 0, formatter: (v) => (Number(v) % 20 === 0 ? v : '') } }),
        visualMap: { type: 'continuous', min: spec.min, max: spec.max, calculable: false, orient: 'horizontal', left: 'center', bottom: 0, itemWidth: 12, itemHeight: small ? 150 : 220,
          text: [fmt(spec.max, 0) + ' °C', fmt(spec.min, 0) + ' °C'], textStyle: { color: t.ink2, fontFamily: t.mono, fontSize: 11 }, inRange: { color: t.h.slice().reverse() } },
        series: [{ type: 'heatmap', data, progressive: 0, emphasis: { itemStyle: { borderColor: t.ink, borderWidth: 1 } },
          markLine: marks.length ? { silent: true, symbol: 'none', lineStyle: { color: t.ink, width: 1, type: 'solid', opacity: 0.55 },
            data: marks } : undefined }] };
    }
    if (spec.kind === 'tornado') {
      const rows = spec.rows, cats = rows.map(r => r.label);
      return { grid: { left: 8, right: 24, top: 30, bottom: 44, containLabel: true },
        tooltip: tip(t, { trigger: 'item', formatter: (p) => { const r = rows[p.value[0]]; return `<b>${esc(r.label)}</b><br>${esc(r.lowLabel || 'low')}: <b style="font-family:${t.mono}">${fmt(r.low, dec(r.low))}</b><br>${esc(r.highLabel || 'high')}: <b style="font-family:${t.mono}">${fmt(r.high, dec(r.high))}</b><br>base: ${fmt(spec.base, dec(spec.base))} ${esc(spec.unit || '')}`; } }),
        yAxis: axis(t, { cat: cats, inverse: true }, { axisLabel: { color: t.ink2, fontFamily: t.body, fontSize: 11, width: innerWidth < 600 ? 120 : 210, overflow: 'break' } }),
        xAxis: axis(t, { name: spec.x, min: spec.min, max: spec.max }, innerWidth < 600 ? { splitNumber: 3 } : {}),
        series: [{ type: 'custom', renderItem: (params, api) => { const i = api.value(0), lo = api.value(1), hi = api.value(2); const a = api.coord([Math.min(lo, hi), i]), b = api.coord([Math.max(lo, hi), i]), c = api.coord([spec.base, i]);
            const h = 14; const kids = [];
            if (c[0] > a[0]) kids.push({ type: 'rect', shape: { x: a[0], y: a[1] - h / 2, width: Math.min(c[0], b[0]) - a[0], height: h, r: 2 }, style: { fill: t.q[1] } });
            if (b[0] > c[0]) kids.push({ type: 'rect', shape: { x: Math.max(c[0], a[0]), y: a[1] - h / 2, width: b[0] - Math.max(c[0], a[0]), height: h, r: 2 }, style: { fill: t.heat } });
            return { type: 'group', children: kids }; },
          data: rows.map((r, i) => [i, r.low, r.high]), encode: { x: [1, 2], y: 0 } },
          { type: 'line', data: [], markLine: { silent: true, symbol: 'none', lineStyle: { color: t.ink2, width: 1, type: 'dashed' }, label: { position: 'start', formatter: 'base ' + fmt(spec.base, dec(spec.base)), color: t.ink2, fontFamily: t.mono, fontSize: 10 }, data: [{ xAxis: spec.base }] } }] };
    }
    if (spec.kind === 'bar') {
      return { grid: { left: 8, right: 16, top: spec.series.length > 1 ? 40 : 20, bottom: 44, containLabel: true }, legend: spec.series.length > 1 ? legend(t) : undefined,
        tooltip: tip(t, { trigger: 'axis', axisPointer: { type: 'shadow' }, formatter: (ps) => `<b>${esc(ps[0].name)}</b><br>` + ps.map(p => `${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value, dec(p.value))}</b>`).join('<br>') }),
        yAxis: spec.horizontal ? axis(t, { cat: spec.cats, inverse: true }, { axisLabel: { color: t.ink2, fontFamily: t.body, fontSize: 11, width: 220, overflow: 'break' } }) : axis(t, { name: spec.y, gap: 40, min: 0 }),
        xAxis: spec.horizontal ? axis(t, { name: spec.y, min: 0 }) : axis(t, { cat: spec.cats, name: spec.x || '', gap: 26 }),
        series: spec.series.map((s, i) => ({ name: s.name, type: 'bar', data: s.data, barMaxWidth: 22, itemStyle: { color: colour(t, s.c, i), borderRadius: spec.horizontal ? [0, 3, 3, 0] : [3, 3, 0, 0] } })) };
    }
    // xy lines / points
    const nleg = spec.series.filter(s => !s.hide).length;
    const ser = spec.series.map((s, i) => { const col = colour(t, s.c, i);
      const o = { name: s.name, type: s.pts ? 'scatter' : 'line', data: s.data, showSymbol: !!s.sym, symbolSize: s.pts ? 11 : 6, symbol: s.symbol || 'circle',
        itemStyle: s.pts ? { color: col, borderColor: t.surface, borderWidth: 1.5 } : { color: col }, z: s.pts ? 5 : 3, connectNulls: false };
      if (!s.pts) o.lineStyle = { width: s.w || 2, color: col, type: s.dash ? 'dashed' : 'solid' };
      if (s.area) o.areaStyle = { color: col, opacity: 0.12 };
      if (s.step) o.step = 'end';
      return o; });
    if (spec.bands) ser.push({ type: 'line', data: [], markArea: { silent: true, itemStyle: { color: t.ink, opacity: 0.05 },
      label: { color: t.ink2, fontFamily: t.body, fontSize: 11, position: 'insideBottom', distance: 8 }, data: spec.bands.map(b => [{ xAxis: b.x0, name: b.name }, { xAxis: b.x1 }]) } });
    if (spec.marks) ser.push({ type: 'line', data: [], markLine: { silent: true, symbol: 'none', lineStyle: { color: t.ink2, width: 1, type: 'dashed' },
      label: { color: t.ink2, fontFamily: t.body, fontSize: 10, formatter: (p) => p.name, position: 'insideEndTop' }, data: spec.marks.map(m => m.y !== undefined ? { yAxis: m.y, name: m.name } : { xAxis: m.x, name: m.name }) } });
    const shown = spec.series.filter(s => !s.hide), rows = nleg > 1 ? legendRows(shown.map(s => s.name), (W || 600) - 4, t) : 0, capTop = 25 * rows + 4;
    return { grid: { left: 8, right: 18, top: capTop + 24, bottom: 44, containLabel: true }, graphic: [caption(t, spec.y.name || '', capTop)],
      legend: nleg > 1 ? legend(t, { data: shown.map(s => ({ name: s.name, icon: icon(s) })) }) : { show: false },
      tooltip: tip(t, { trigger: 'axis', axisPointer: { type: 'line', lineStyle: { color: t.muted, width: 1 } }, formatter: (ps) => { if (!ps.length) return ''; const x = ps[0].value[0];
        return `<div style="color:${t.muted};font-size:11px">${esc(spec.xl || spec.x.name || '')} ${fmt(x, dec(x))}</div>` + ps.filter(p => p.value[1] != null).map(p => `<div style="display:flex;gap:8px;align-items:center;margin-top:3px"><span style="width:14px;height:2px;background:${p.color};display:inline-block"></span><b style="font-family:${t.mono}">${fmt(p.value[1], dec(p.value[1]))}</b><span style="color:${t.ink2}">${esc(p.seriesName)}</span></div>`).join(''); } }),
      xAxis: axis(t, spec.x), yAxis: axis(t, Object.assign({}, spec.y, { name: '' })), series: ser };
  }
  function tableOf(spec) {
    if (spec.kind === 'field') { const step = Math.max(1, Math.round(spec.z.length / 15));
      let h = '<table><thead><tr><th>Height (cm) \\ radius (cm)</th>' + spec.r.filter((_, i) => i % 3 === 0).map(r => `<th>${esc(r)}</th>`).join('') + '</tr></thead><tbody>';
      for (let j = spec.z.length - 1; j >= 0; j -= step) h += `<tr><td>${esc(spec.z[j])}</td>` + spec.r.map((_, i) => i).filter(i => i % 3 === 0).map(i => `<td>${fmt(spec.T[j][i], 1)}</td>`).join('') + '</tr>';
      return h + '</tbody></table>'; }
    if (spec.kind === 'tornado') return '<table><thead><tr><th>Variable</th><th>Low case</th><th>High case</th></tr></thead><tbody>' + spec.rows.map(r => `<tr><td>${esc(r.label)}</td><td>${fmt(r.low, dec(r.low))}</td><td>${fmt(r.high, dec(r.high))}</td></tr>`).join('') + `</tbody></table>`;
    if (spec.kind === 'bar') return '<table><thead><tr><th></th>' + spec.series.map(s => `<th>${esc(s.name)}</th>`).join('') + '</tr></thead><tbody>' + spec.cats.map((c, i) => `<tr><td>${esc(c)}</td>` + spec.series.map(s => `<td>${fmt(s.data[i], dec(s.data[i] || 0))}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
    let h = `<table><thead><tr><th>Series</th><th>${esc(spec.x.name || 'x')}</th><th>${esc(spec.y.name || 'y')}</th></tr></thead><tbody>`;
    for (const s of spec.series) for (const p of s.data) if (p[1] != null) h += `<tr><td>${esc(s.name)}</td><td>${fmt(p[0], dec(p[0]))}</td><td>${fmt(p[1], dec(p[1]))}</td></tr>`;
    return h + '</tbody></table>';
  }
  for (const [id, spec] of Object.entries(DATA.charts)) {
    const el = document.getElementById(id);
    if (!el || typeof echarts === 'undefined') continue;
    const c = echarts.init(el, null, { renderer: 'canvas' });
    const render = () => c.setOption(build(spec, T(), el.clientWidth), true);
    render(); charts.push({ c, render }); tables[id] = () => tableOf(spec);
  }
  document.querySelectorAll('[data-table]').forEach(b => { const id = b.dataset.table; b.setAttribute('aria-expanded', 'false');
    b.addEventListener('click', () => { let box = document.getElementById(id + '-tbl');
      if (!box) { box = document.createElement('div'); box.id = id + '-tbl'; box.className = 'dtable'; box.hidden = true; b.closest('figure').appendChild(box); }
      box.hidden = !box.hidden; b.setAttribute('aria-expanded', String(!box.hidden)); b.textContent = box.hidden ? 'Table' : 'Hide table';
      if (!box.hidden && tables[id]) box.innerHTML = tables[id](); }); });
  const rerender = () => charts.forEach(x => x.render());
  let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => charts.forEach(x => { x.c.resize(); x.render(); }), 120); });
  try { document.fonts.ready.then(rerender); } catch (e) {}
  try { window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', rerender); } catch (e) {}
  new MutationObserver(rerender).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
})();
