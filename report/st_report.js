(function () {
  const css = () => getComputedStyle(document.documentElement);
  const tok = (n) => css().getPropertyValue(n).trim();
  const T = () => ({ ink: tok('--ink'), ink2: tok('--ink-2'), muted: tok('--muted'), grid: tok('--grid'), axis: tok('--axis'), surface: tok('--surface'),
    rule: tok('--rule'), ref: tok('--ref'), heat: tok('--heat'), accent: tok('--accent'),
    s: [tok('--s1'), tok('--s2'), tok('--s3'), tok('--s4'), tok('--s5'), tok('--s6'), tok('--s7')], q: [tok('--q1'), tok('--q2'), tok('--q3'), tok('--q4'), tok('--q5')],
    mono: tok('--f-mono'), body: tok('--f-body') });
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = (v, d) => (v === null || v === undefined || !isFinite(v)) ? '–' : Number(v).toLocaleString('en-GB', { minimumFractionDigits: d, maximumFractionDigits: d }).replace('-', '−');
  const dec = (v) => { const a = Math.abs(v); return a >= 100 ? 0 : a >= 10 ? 1 : a >= 1 ? 2 : 3; };
  const LINE = 'path://M0,4 L16,4 L16,6 L0,6 Z';
  const charts = [], tables = {};
  const colour = (t, c, i) => c === 'ref' ? t.ref : c === 'heat' ? t.heat : c === 'muted' ? t.muted : (typeof c === 'string' && c[0] === 'q') ? t.q[+c[1]] : t.s[(c === undefined || c === null ? i : c) % 7];
  function axis(t, a, o) {
    a = a || {};
    const base = { type: a.log ? 'log' : (a.cat ? 'category' : 'value'), name: a.name || '', nameLocation: 'middle', nameGap: a.gap || 30, logBase: 10,
      nameTextStyle: { color: t.ink2, fontFamily: t.body, fontSize: 12 }, axisLine: { show: true, lineStyle: { color: t.axis } }, axisTick: { lineStyle: { color: t.axis } },
      axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: a.log ? (v) => String(v) : undefined }, splitLine: { show: !a.cat, lineStyle: { color: t.grid } }, scale: true };
    if (a.min !== undefined) base.min = a.min; if (a.max !== undefined) base.max = a.max; if (a.cat) base.data = a.cat; if (a.inverse) base.inverse = true;
    return Object.assign(base, o || {});
  }
  const tip = (t, o) => Object.assign({ backgroundColor: t.surface, borderColor: t.rule, borderWidth: 1, padding: [8, 10], textStyle: { color: t.ink, fontFamily: t.body, fontSize: 12 },
    extraCssText: 'box-shadow:0 2px 10px rgba(0,0,0,.12);border-radius:6px;' }, o || {});
  const legend = (t, o) => Object.assign({ top: 0, left: 0, itemWidth: 16, itemHeight: 8, itemGap: 12, textStyle: { color: t.ink2, fontFamily: t.body, fontSize: 11.5 }, icon: 'roundRect' }, o || {});

  function build(spec, t) {
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
      const o = { name: s.name, type: s.pts ? 'scatter' : 'line', data: s.data, showSymbol: !!s.sym, symbolSize: s.pts ? 8 : 6, symbol: s.symbol || 'circle', itemStyle: { color: col }, z: 3, connectNulls: false };
      if (!s.pts) o.lineStyle = { width: s.w || 2, color: col, type: s.dash ? 'dashed' : 'solid' };
      if (s.area) o.areaStyle = { color: col, opacity: 0.12 };
      if (s.step) o.step = 'end';
      return o; });
    if (spec.marks) ser.push({ type: 'line', data: [], markLine: { silent: true, symbol: 'none', lineStyle: { color: t.ink2, width: 1, type: 'dashed' },
      label: { color: t.ink2, fontFamily: t.body, fontSize: 10, formatter: (p) => p.name, position: 'insideEndTop' }, data: spec.marks.map(m => m.y !== undefined ? { yAxis: m.y, name: m.name } : { xAxis: m.x, name: m.name }) } });
    return { grid: { left: 8, right: 18, top: 18 + 25 * Math.ceil(nleg / (spec.legcols || 2)), bottom: 44, containLabel: true }, legend: legend(t, { icon: LINE, data: spec.series.filter(s => !s.hide).map(s => s.name) }),
      tooltip: tip(t, { trigger: 'axis', axisPointer: { type: 'line', lineStyle: { color: t.muted, width: 1 } }, formatter: (ps) => { if (!ps.length) return ''; const x = ps[0].value[0];
        return `<div style="color:${t.muted};font-size:11px">${esc(spec.xl || spec.x.name || '')} ${fmt(x, dec(x))}</div>` + ps.filter(p => p.value[1] != null).map(p => `<div style="display:flex;gap:8px;align-items:center;margin-top:3px"><span style="width:14px;height:2px;background:${p.color};display:inline-block"></span><b style="font-family:${t.mono}">${fmt(p.value[1], dec(p.value[1]))}</b><span style="color:${t.ink2}">${esc(p.seriesName)}</span></div>`).join(''); } }),
      xAxis: axis(t, spec.x), yAxis: axis(t, Object.assign({ gap: 42 }, spec.y)), series: ser };
  }
  function tableOf(spec) {
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
    const render = () => c.setOption(build(spec, T()), true);
    render(); charts.push({ c, render }); tables[id] = () => tableOf(spec);
  }
  document.querySelectorAll('[data-table]').forEach(b => { const id = b.dataset.table; b.setAttribute('aria-expanded', 'false');
    b.addEventListener('click', () => { let box = document.getElementById(id + '-tbl');
      if (!box) { box = document.createElement('div'); box.id = id + '-tbl'; box.className = 'dtable'; box.hidden = true; b.closest('figure').appendChild(box); }
      box.hidden = !box.hidden; b.setAttribute('aria-expanded', String(!box.hidden)); b.textContent = box.hidden ? 'Table' : 'Hide table';
      if (!box.hidden && tables[id]) box.innerHTML = tables[id](); }); });
  let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => charts.forEach(x => x.c.resize()), 120); });
  const rerender = () => charts.forEach(x => x.render());
  try { window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', rerender); } catch (e) {}
  new MutationObserver(rerender).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
})();
