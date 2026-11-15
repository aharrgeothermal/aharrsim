(function () {
  const css = () => getComputedStyle(document.documentElement);
  const tok = (n) => css().getPropertyValue(n).trim();
  const T = () => ({
    ink: tok('--ink'), ink2: tok('--ink-2'), muted: tok('--muted'), grid: tok('--grid'), axis: tok('--axis'),
    surface: tok('--surface'), rule: tok('--rule'), ref: tok('--ref'), heat: tok('--heat'), accent: tok('--accent'),
    s: [tok('--s1'), tok('--s2'), tok('--s3'), tok('--s4'), tok('--s5'), tok('--s6'), tok('--s7')],
    q: [tok('--q1'), tok('--q2'), tok('--q3'), tok('--q4'), tok('--q5')],
    mono: tok('--f-mono'), body: tok('--f-body')
  });
  const charts = [];
  const tables = {};
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = (v, d) => (v === null || v === undefined || !isFinite(v)) ? '–' :
    Number(v).toLocaleString('en-GB', { minimumFractionDigits: d, maximumFractionDigits: d }).replace('-', '−');
  const LINE_ICON = 'path://M0,4 L16,4 L16,6 L0,6 Z';

  function axis(t, o) {
    return Object.assign({
      nameLocation: 'middle', nameGap: 30, nameTextStyle: { color: t.ink2, fontFamily: t.body, fontSize: 12 },
      axisLine: { show: true, lineStyle: { color: t.axis } }, axisTick: { lineStyle: { color: t.axis } },
      axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11 },
      splitLine: { show: true, lineStyle: { color: t.grid, width: 1 } }, scale: true
    }, o || {});
  }
  function tip(t, o) {
    return Object.assign({
      backgroundColor: t.surface, borderColor: t.rule, borderWidth: 1, padding: [8, 10],
      textStyle: { color: t.ink, fontFamily: t.body, fontSize: 12 }, extraCssText: 'box-shadow:0 2px 10px rgba(0,0,0,.12);border-radius:6px;'
    }, o || {});
  }
  function legend(t, o) {
    return Object.assign({ top: 0, left: 0, itemWidth: 16, itemHeight: 8, itemGap: 12, textStyle: { color: t.ink2, fontFamily: t.body, fontSize: 11.5 },
      icon: 'roundRect' }, o || {});
  }
  function lineS(name, data, color, extra) {
    return Object.assign({ name, type: 'line', data, showSymbol: false, symbolSize: 8, lineStyle: { width: 2, color }, itemStyle: { color },
      emphasis: { focus: 'series' }, z: 3 }, extra || {});
  }
  function ptsS(name, data, color, extra) {
    return Object.assign({ name, type: 'scatter', data, symbolSize: 9, itemStyle: { color, borderColor: T().surface, borderWidth: 2 }, z: 4 }, extra || {});
  }
  function axisTip(t, xl, yl, xd, yd, swap) {
    return tip(t, { trigger: 'axis', axisPointer: { type: 'line', axis: swap ? 'y' : 'x', lineStyle: { color: t.muted, width: 1 } },
      formatter: (ps) => {
        if (!ps.length) return '';
        const v0 = ps[0].value; const head = swap ? `${yl} ${fmt(v0[1], yd)}` : `${xl} ${fmt(v0[0], xd)}`;
        const rows = ps.filter(p => p.value && p.value[swap ? 0 : 1] != null).map(p => `<div style="display:flex;align-items:center;gap:8px;margin-top:3px"><span style="display:inline-block;width:14px;height:2px;background:${p.color}"></span><b style="font-family:${t.mono}">${fmt(swap ? p.value[0] : p.value[1], swap ? xd : yd)}</b><span style="color:${t.ink2}">${esc(p.seriesName)}</span></div>`);
        return `<div style="color:${t.muted};font-size:11px">${esc(head)}</div>` + rows.join('');
      } });
  }
  function mount(id, build, tableFn) {
    const el = document.getElementById(id);
    if (!el || typeof echarts === 'undefined') return;
    const c = echarts.init(el, null, { renderer: 'canvas' });
    const render = () => { c.setOption(build(T()), true); };
    render();
    charts.push({ c, render });
    if (tableFn) tables[id] = tableFn;
  }
  function seriesTable(obj, xl, yl, xd, yd) {
    return () => {
      let h = `<table><thead><tr><th>Series</th><th>${esc(xl)}</th><th>${esc(yl)}</th></tr></thead><tbody>`;
      for (const [k, arr] of Object.entries(obj)) for (const p of arr) if (p[1] != null) h += `<tr><td>${esc(k)}</td><td>${fmt(p[0], xd)}</td><td>${fmt(p[1], yd)}</td></tr>`;
      return h + '</tbody></table>';
    };
  }
  const zip = (x, y) => x.map((v, i) => [v, y[i]]);

  // =====================================================================================================
  // Schematics: a small declarative renderer (depths in km, x in 0..1)
  // =====================================================================================================
  const HOT = 'var(--heat)', COLD = 'var(--accent)';
  function drawScene(el, spec, title) {
    const W = 360, H = spec.surface_only ? 110 : 300, gy = 86, gb = 286, x0 = 46, x1 = 334;
    const X = (f) => x0 + f * (x1 - x0);
    const Z = (z) => spec.depth > 0 ? gy + (z / spec.depth) * (gb - gy) : gy;
    const parts = [];
    const text = (x, y, s, o) => { o = o || {}; parts.push(`<text x="${x.toFixed(1)}" y="${y.toFixed(1)}" font-size="${o.fs || 9.5}" fill="${o.fill || 'currentColor'}" text-anchor="${o.anchor || 'start'}" ${o.weight ? `font-weight="${o.weight}"` : ''} ${o.extra || ''}>${esc(s)}</text>`); };
    parts.push(`<defs><marker id="${title}-ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="${HOT}"/></marker>
      <marker id="${title}-ac" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="${COLD}"/></marker>
      <marker id="${title}-ag" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>`);
    if (!spec.surface_only) {
      parts.push(`<rect x="${x0 - 10}" y="${gy}" width="${x1 - x0 + 20}" height="${gb - gy}" fill="var(--rock)"/>`);
      // bands
      for (const b of (spec.bands || [])) {
        const y = Z(b.z), h = Math.max(6, Z(b.z + b.h) - y);
        const fill = { aquifer: 'var(--water)', store: 'var(--store)', hot: 'var(--steam)', steam: 'var(--steam)' }[b.kind] || 'var(--rock-2)';
        parts.push(`<rect x="${x0 - 10}" y="${y.toFixed(1)}" width="${x1 - x0 + 20}" height="${h.toFixed(1)}" fill="${fill}" ${b.kind === 'steam' ? 'fill-opacity="0.8"' : ''}/>`);
        if (b.label) text(X(0.02), y + h - 4 > y + 11 ? y + 11 : y - 3, b.label, { fill: 'var(--muted)', fs: 9 });
      }
      // depth ticks
      for (const z of (spec.ticks || [])) {
        const y = Z(z);
        parts.push(`<line x1="${x0 - 10}" y1="${y}" x2="${x0 - 4}" y2="${y}" stroke="currentColor" stroke-width="1"/>`);
        text(x0 - 13, y + 3.5, spec.depth < 1 ? `${Math.round(z * 1000)} m` : `${z} km`, { anchor: 'end', fs: 9, fill: 'var(--muted)' });
      }
      parts.push(`<line x1="${x0 - 10}" y1="${gy}" x2="${x0 - 10}" y2="${gb}" stroke="currentColor" stroke-width="1"/>`);
      // isotherms
      for (const [z, lab] of (spec.isotherms || [])) {
        const y = Z(z);
        parts.push(`<line x1="${x1 - 26}" y1="${y}" x2="${x1 + 10}" y2="${y}" stroke="${HOT}" stroke-width="1" stroke-dasharray="3 2"/>`);
        text(x1 + 8, y - 3, lab, { fill: HOT, fs: 9, anchor: 'end' });
      }
    }
    const wellhead = {}, unitpos = {};
    // fractures
    const fr = spec.fractures;
    if (fr && !spec.surface_only) {
      const n = fr.n, xa = X(fr.x1) + 10, xb = X(fr.x2) - 10, hpx = Math.max(10, (fr.h / spec.depth) * (gb - gy));
      if (fr.at) {
        for (const fx of fr.at) parts.push(`<ellipse cx="${X(fx).toFixed(1)}" cy="${Z(fr.z).toFixed(1)}" rx="3" ry="${(hpx / 2).toFixed(1)}" fill="${HOT}" fill-opacity="0.45" stroke="${HOT}" stroke-width="0.8"/>`);
      } else if (fr.horizontal_pair) {
        const ya = Z(fr.z) - 11, yb = Z(fr.z) + 11;
        for (let i = 0; i < n; i++) { const x = xa + (xb - xa) * (i + 0.5) / n; parts.push(`<line x1="${x.toFixed(1)}" y1="${ya}" x2="${x.toFixed(1)}" y2="${yb}" stroke="${HOT}" stroke-width="1.6" stroke-opacity="0.7"/>`); }
      } else {
        for (let i = 0; i < n; i++) {
          const x = n === 1 ? (xa + xb) / 2 : xa + (xb - xa) * i / (n - 1);
          parts.push(`<ellipse cx="${x.toFixed(1)}" cy="${Z(fr.z).toFixed(1)}" rx="${n === 1 ? 3 : 1.6}" ry="${(hpx / 2).toFixed(1)}" fill="${HOT}" fill-opacity="0.45" stroke="${HOT}" stroke-width="0.8"/>`);
        }
      }
    }
    // wells
    (spec.wells || []).forEach((w, i) => {
      const x = X(w.x), yb = Z(w.z);
      wellhead['w' + i] = x;
      const hot = w.dir === 'up' || w.dir === 'loop' || w.dir === 'cycle';
      const col = w.dir === 'down' ? COLD : HOT;
      const mk = (c) => `url(#${title}-${c === HOT ? 'ah' : 'ac'})`;
      if (w.kind === 'single') {
        const wd = w.tubing ? 5 : 7;
        parts.push(`<rect x="${(x - wd / 2).toFixed(1)}" y="${gy - 2}" width="${wd}" height="${(yb - gy + 2).toFixed(1)}" fill="var(--surface)" stroke="currentColor" stroke-width="1"/>`);
        if (w.steam) parts.push(`<rect x="${(x - wd / 2 + 1).toFixed(1)}" y="${gy}" width="${wd - 2}" height="${(yb - gy).toFixed(1)}" fill="var(--steam)"/>`);
        if (w.lateral) {
          const x2 = X(w.lateral.x2), n = w.lateral.n || 1, dy = (w.lateral.dy || 0) * 100;
          for (let k = 0; k < n; k++) {
            const off = n > 1 ? (k - (n - 1) / 2) * (22 / (n - 1)) : 0;
            parts.push(`<path d="M${x.toFixed(1)},${(yb + dy).toFixed(1)} Q${x.toFixed(1)},${(yb + dy + off).toFixed(1)} ${(x + 14).toFixed(1)},${(yb + dy + off).toFixed(1)} H${x2.toFixed(1)}" fill="none" stroke="currentColor" stroke-width="${n > 1 ? 1.4 : 3}"/>`);
          }
          if (dy !== 0) parts.push(`<line x1="${x}" y1="${yb}" x2="${x}" y2="${(yb + dy).toFixed(1)}" stroke="currentColor" stroke-width="3"/>`);
        }
        // flow chevrons
        const L = yb - gy, nA = Math.max(2, Math.min(4, Math.round(L / 45)));
        for (let k = 1; k <= nA; k++) {
          const y = gy + L * k / (nA + 1);
          if (w.dir === 'cycle') parts.push(`<path d="M${(x - 2.5).toFixed(1)},${(y + 3).toFixed(1)} L${x},${(y - 1).toFixed(1)} L${(x + 2.5).toFixed(1)},${(y + 3).toFixed(1)}" fill="none" stroke="${HOT}" stroke-width="1.4"/><path d="M${(x - 2.5).toFixed(1)},${(y + 6).toFixed(1)} L${x},${(y + 10).toFixed(1)} L${(x + 2.5).toFixed(1)},${(y + 6).toFixed(1)}" fill="none" stroke="${COLD}" stroke-width="1.4"/>`);
          else if (w.dir === 'up') parts.push(`<path d="M${(x - 2.5).toFixed(1)},${(y + 3).toFixed(1)} L${x},${(y - 2).toFixed(1)} L${(x + 2.5).toFixed(1)},${(y + 3).toFixed(1)}" fill="none" stroke="${HOT}" stroke-width="1.5"/>`);
          else parts.push(`<path d="M${(x - 2.5).toFixed(1)},${(y - 3).toFixed(1)} L${x},${(y + 2).toFixed(1)} L${(x + 2.5).toFixed(1)},${(y - 3).toFixed(1)}" fill="none" stroke="${COLD}" stroke-width="1.5"/>`);
        }
        if (w.flash) { const yf = Z(w.flash); for (let k = 0; k < 3; k++) parts.push(`<circle cx="${(x + (k % 2 ? 1.5 : -1.5)).toFixed(1)}" cy="${(yf - 8 * k).toFixed(1)}" r="1.4" fill="none" stroke="${HOT}" stroke-width="0.8"/>`); text(x + 7, yf + 3, 'flash', { fill: HOT, fs: 8.5 }); }
      } else if (w.kind === 'coaxial') {
        parts.push(`<rect x="${(x - 6).toFixed(1)}" y="${gy - 2}" width="12" height="${(yb - gy + 2).toFixed(1)}" fill="var(--surface)" stroke="currentColor" stroke-width="1"/>`);
        parts.push(`<rect x="${(x - 2).toFixed(1)}" y="${gy - 2}" width="4" height="${(yb - gy - (w.open_bottom ? 6 : 6)).toFixed(1)}" fill="${HOT}" fill-opacity="0.85" stroke="currentColor" stroke-width="0.6"/>`);
        if (w.open_bottom) { parts.push(`<line x1="${(x - 6).toFixed(1)}" y1="${yb}" x2="${(x + 6).toFixed(1)}" y2="${yb}" stroke="var(--rock)" stroke-width="2"/>`); }
        const L = yb - gy, nA = Math.max(2, Math.min(4, Math.round(L / 45)));
        for (let k = 1; k <= nA; k++) {
          const y = gy + L * k / (nA + 1);
          for (const dx of [-4, 4]) parts.push(`<path d="M${(x + dx - 1.4).toFixed(1)},${(y - 2.5).toFixed(1)} L${(x + dx).toFixed(1)},${(y + 1.5).toFixed(1)} L${(x + dx + 1.4).toFixed(1)},${(y - 2.5).toFixed(1)}" fill="none" stroke="${COLD}" stroke-width="1.2"/>`);
          parts.push(`<path d="M${(x - 1.4).toFixed(1)},${(y + 2.5).toFixed(1)} L${x},${(y - 1.5).toFixed(1)} L${(x + 1.4).toFixed(1)},${(y + 2.5).toFixed(1)}" fill="none" stroke="var(--surface)" stroke-width="1.2"/>`);
        }
      } else if (w.kind === 'utube') {
        const d = 5;
        parts.push(`<path d="M${(x - d).toFixed(1)},${gy - 2} V${(yb - 6).toFixed(1)} A${d},6 0 0 0 ${(x + d).toFixed(1)},${(yb - 6).toFixed(1)} V${gy - 2}" fill="none" stroke="currentColor" stroke-width="3"/>`);
        parts.push(`<path d="M${(x - d).toFixed(1)},${gy - 2} V${(yb - 6).toFixed(1)} A${d},6 0 0 0 ${(x + d).toFixed(1)},${(yb - 6).toFixed(1)} V${gy - 2}" fill="none" stroke="var(--surface)" stroke-width="1.4"/>`);
        const L = yb - gy, nA = 3;
        for (let k = 1; k <= nA; k++) {
          const y = gy + L * k / (nA + 1);
          parts.push(`<path d="M${(x - d - 1.6).toFixed(1)},${(y - 2.5).toFixed(1)} L${(x - d).toFixed(1)},${(y + 1.5).toFixed(1)} L${(x - d + 1.6).toFixed(1)},${(y - 2.5).toFixed(1)}" fill="none" stroke="${COLD}" stroke-width="1.3"/>`);
          parts.push(`<path d="M${(x + d - 1.6).toFixed(1)},${(y + 2.5).toFixed(1)} L${(x + d).toFixed(1)},${(y - 1.5).toFixed(1)} L${(x + d + 1.6).toFixed(1)},${(y + 2.5).toFixed(1)}" fill="none" stroke="${HOT}" stroke-width="1.3"/>`);
        }
      }
      const others = (spec.wells || []).filter((o, j) => j !== i);
      const nearRight = others.some(o => o.x > w.x && o.x - w.x < 0.2), nearLeft = others.some(o => o.x < w.x && w.x - o.x < 0.2);
      const labelLeft = w.label_side ? w.label_side === 'left' : (w.x < 0.5 ? !nearLeft : nearRight);
      if (w.pump) { const yp = Z(w.pump); parts.push(`<circle cx="${x.toFixed(1)}" cy="${yp.toFixed(1)}" r="5" fill="var(--surface)" stroke="currentColor" stroke-width="1.2"/><path d="M${(x - 2.2).toFixed(1)},${(yp + 2).toFixed(1)} L${x},${(yp - 2.5).toFixed(1)} L${(x + 2.2).toFixed(1)},${(yp + 2).toFixed(1)} Z" fill="currentColor"/>`);
        text(labelLeft ? x + 9 : x - 9, yp + 3.5, 'pump', { fs: 8.5, fill: 'var(--muted)', anchor: labelLeft ? 'start' : 'end' }); }
      if (w.label) text(labelLeft ? x - 9 : x + 9, gy + 12, w.label, { fs: 8.8, fill: 'var(--ink-2)', anchor: labelLeft ? 'end' : 'start' });
    });
    // notes
    for (const n of (spec.notes || [])) text(X(n.x), Z(n.z), n.text, { fs: 8.8, fill: 'var(--ink-2)', anchor: 'middle' });
    // ground line
    if (!spec.surface_only) parts.push(`<line x1="${x0 - 10}" y1="${gy}" x2="${x1 + 10}" y2="${gy}" stroke="currentColor" stroke-width="1.4"/>`);
    // units
    const uy = spec.surface_only ? 30 : 16, uh = 24, uw = 72;
    (spec.units || []).forEach((u, j) => {
      const x = X(u.x); unitpos['u' + j] = x;
      const bar = { orc: COLD, flash: HOT, steam: HOT, hp: 'var(--s3)', heat: 'var(--q2)', hx: 'var(--ref)', sep: 'var(--ref)', heater: 'var(--s4)' }[u.kind] || 'currentColor';
      const w = Math.max(uw, 6.2 * u.label.length + 10);
      parts.push(`<rect x="${(x - w / 2).toFixed(1)}" y="${uy}" width="${w.toFixed(1)}" height="${uh}" rx="3" fill="var(--surface)" stroke="currentColor" stroke-width="1"/>`);
      parts.push(`<rect x="${(x - w / 2).toFixed(1)}" y="${uy}" width="${w.toFixed(1)}" height="3" rx="1" fill="${bar}"/>`);
      text(x, uy + 16.5, u.label, { fs: 9.2, anchor: 'middle', weight: 500, fill: 'var(--ink)' });
    });
    // flows
    const ub = uy + uh;
    (spec.flows || []).forEach((f, k) => {
      const col = f.hot ? HOT : COLD, m = `url(#${title}-${f.hot ? 'ah' : 'ac'})`;
      const isU = (s) => s[0] === 'u';
      if (isU(f.a) && isU(f.b)) {
        const xa = unitpos[f.a], xb = unitpos[f.b], y = uy + uh / 2;
        const la = Math.max(uw, 6.2 * spec.units[+f.a[1]].label.length + 10) / 2, lb = Math.max(uw, 6.2 * spec.units[+f.b[1]].label.length + 10) / 2;
        parts.push(`<line x1="${(xa + la).toFixed(1)}" y1="${y}" x2="${(xb - lb - 1).toFixed(1)}" y2="${y}" stroke="${col}" stroke-width="1.6" marker-end="${m}"/>`);
        text((xa + la + xb - lb) / 2, y - 5, f.label, { fs: 8.8, anchor: 'middle', fill: 'var(--ink-2)' });
        return;
      }
      const fromWell = !isU(f.a);
      const wx = fromWell ? wellhead[f.a] : wellhead[f.b];
      const ux = fromWell ? unitpos[f.b] : unitpos[f.a];
      const off = f.back ? (f.hot ? -4 : 4) : 0;
      const yr = f.hot ? 58 : 72;  // route height
      const xw = wx + off;
      let d;
      if (fromWell) d = `M${xw.toFixed(1)},${gy - 1} V${yr} H${(ux + (f.back ? -4 : 0)).toFixed(1)} V${ub + 1}`;
      else d = `M${(ux + (f.back ? 4 : 0)).toFixed(1)},${ub} V${yr} H${xw.toFixed(1)} V${gy - 2}`;
      parts.push(`<path d="${d}" fill="none" stroke="${col}" stroke-width="1.6" marker-end="${m}"/>`);
      if (f.label) { const xm = (xw + ux) / 2; text(xm, f.hot ? yr - 4 : yr + 10, f.label, { fs: 8.8, anchor: 'middle', fill: 'var(--ink-2)' }); }
    });
    el.innerHTML = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(title)}" xmlns="http://www.w3.org/2000/svg">${parts.join('')}</svg>`;
  }
  for (const c of DATA.cases) {
    const el = document.getElementById('sc-' + c.key);
    if (el && c.spec) drawScene(el, c.spec, 'Schematic of ' + c.short);
  }

  // =====================================================================================================
  // Benchmarks
  // =====================================================================================================
  const B = DATA.bench;
  mount('b-steady', (t) => {
    const ser = B.steady.map((s, i) => lineS(s.label, s.rows.map(r => [r.dz, r.err_max]), t.s[i], { showSymbol: true, symbol: 'circle', symbolSize: 7 }));
    const r0 = B.steady[0].rows; const ref = [[400, r0[0].err_max * 1.6], [25, r0[0].err_max * 1.6 * Math.pow(25 / 400, 2)]];
    ser.push(lineS('slope 2 (reference)', ref, t.ref, { lineStyle: { width: 1.2, color: t.ref, type: 'dashed' } }));
    return { grid: { left: 60, right: 14, top: 56, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: tip(t, { trigger: 'item', formatter: (p) => `<b style="font-family:${t.mono}">${Number(p.value[1]).toExponential(2)} K</b> at dz ${p.value[0]} m<br><span style="color:${t.ink2}">${esc(p.seriesName)}</span>` }),
      xAxis: axis(t, { type: 'log', name: 'Cell size dz (m)', min: 20, max: 500, logBase: 10 }),
      yAxis: axis(t, { type: 'log', name: 'Max error (K)', nameGap: 44, logBase: 10, axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: (v) => Number(v).toExponential(0) } }),
      series: ser };
  }, () => { let h = '<table><thead><tr><th>Case</th><th>dz (m)</th><th>Max error (K)</th><th>Order</th></tr></thead><tbody>';
    for (const s of B.steady) s.rows.forEach((r, i) => { h += `<tr><td>${esc(s.label)}</td><td>${r.dz}</td><td>${Number(r.err_max).toExponential(2)}</td><td>${i ? fmt(s.order[i - 1], 3) : '–'}</td></tr>`; }); return h + '</tbody></table>'; });

  mount('b-ramey', (t) => {
    const R = B.ramey;
    return { grid: { left: 54, right: 14, top: 56, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: axisTip(t, 'Day', '', 0, 2, false),
      xAxis: axis(t, { type: 'log', name: 'Time (days)', min: 0.8, max: 4000, logBase: 10 }),
      yAxis: axis(t, { type: 'value', name: 'Wellhead temperature (°C)', nameGap: 40, min: 148, max: 160 }),
      series: [lineS('aharrsim, IAPWS water', R.map(r => [r.t_days, r.T_wh_model]), t.s[0], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('aharrsim, constant properties', R.map(r => [r.t_days, r.T_wh_model_constant_fluid]), t.s[2], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('Ramey / Hasan–Kabir', R.map(r => [r.t_days, r.T_wh_ramey]), t.ref, { showSymbol: true, symbol: 'diamond', symbolSize: 8, lineStyle: { width: 1.5, color: t.ref, type: 'dashed' } })] };
  }, () => { let h = '<table><thead><tr><th>Day</th><th>aharrsim IAPWS (°C)</th><th>aharrsim constant (°C)</th><th>Ramey (°C)</th></tr></thead><tbody>';
    for (const r of B.ramey) h += `<tr><td>${r.t_days}</td><td>${fmt(r.T_wh_model, 2)}</td><td>${fmt(r.T_wh_model_constant_fluid, 2)}</td><td>${fmt(r.T_wh_ramey, 2)}</td></tr>`; return h + '</tbody></table>'; });

  mount('b-utube', (t) => {
    const R = B.utube;
    return { grid: { left: 50, right: 14, top: 56, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: axisTip(t, 'Day', '', 0, 2, false),
      xAxis: axis(t, { type: 'log', name: 'Time (days)', min: 0.8, max: 2500, logBase: 10 }),
      yAxis: axis(t, { type: 'value', name: 'Temperature (°C)', nameGap: 36, min: 3, max: 9 }),
      series: [lineS('borehole wall, aharrsim', R.map(r => [r.t_days, r.T_b_aharrsim]), t.s[1], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('borehole wall, pygfunction', R.map(r => [r.t_days, r.T_b_pyg]), t.ref, { showSymbol: true, symbol: 'diamond', symbolSize: 8, lineStyle: { width: 1.5, color: t.ref, type: 'dashed' } }),
               lineS('outlet, aharrsim', R.map(r => [r.t_days, r.T_out_aharrsim]), t.s[0], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('outlet, pygfunction', R.map(r => [r.t_days, r.T_out_pyg]), t.q[0], { showSymbol: true, symbol: 'diamond', symbolSize: 8, lineStyle: { width: 1.5, color: t.q[0], type: 'dashed' } })] };
  }, () => { let h = '<table><thead><tr><th>Day</th><th>Q (W)</th><th>T wall aharrsim</th><th>T wall pyg</th><th>T out aharrsim</th><th>T out pyg</th></tr></thead><tbody>';
    for (const r of B.utube) h += `<tr><td>${r.t_days}</td><td>${fmt(r.Q_W, 0)}</td><td>${fmt(r.T_b_aharrsim, 2)}</td><td>${fmt(r.T_b_pyg, 2)}</td><td>${fmt(r.T_out_aharrsim, 2)}</td><td>${fmt(r.T_out_pyg, 2)}</td></tr>`; return h + '</tbody></table>'; });

  for (const w of B.twophase) {
    mount('b-tp-' + w.key, (t) => ({ grid: { left: 54, right: 14, top: 50, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: tip(t, { trigger: 'item', formatter: (p) => `<b style="font-family:${t.mono}">${fmt(p.value[0], 1)} bar</b> at ${fmt(p.value[1], 0)} m<br><span style="color:${t.ink2}">${esc(p.seriesName)}</span>` }),
      xAxis: axis(t, { type: 'value', name: 'Pressure (bar g)', min: 0 }),
      yAxis: axis(t, { type: 'value', name: 'Depth (m)', nameGap: 40, inverse: true, min: 0 }),
      series: [lineS('aharrsim', zip(w.profile.p_bar, w.profile.md), t.s[0]),
               ptsS('measured', zip(w.measured.p_bar, w.measured.md), t.ref, { symbol: 'diamond', symbolSize: 9 })] }),
      () => { let h = '<table><thead><tr><th>Depth (m)</th><th>Measured (bar)</th><th>aharrsim (bar)</th></tr></thead><tbody>';
        w.measured.md.forEach((z, i) => { const pm = interp(z, w.profile.md, w.profile.p_bar); h += `<tr><td>${fmt(z, 0)}</td><td>${fmt(w.measured.p_bar[i], 1)}</td><td>${fmt(pm, 1)}</td></tr>`; }); return h + '</tbody></table>'; });
  }
  function interp(x, xs, ys) { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) { const f = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + f * (ys[i] - ys[i - 1]); } return ys[ys.length - 1]; }

  mount('b-uloop', (t) => {
    const ser = [];
    B.uloop.forEach((c, i) => {
      const col = [t.s[0], t.s[1], t.s[2], t.s[3], t.s[4], t.s[5]][i];
      ser.push(lineS(c.label, c.comparison.map(r => [r[0], r[1]]), col, { showSymbol: false }));
      ser.push(ptsS(c.label + ' (GeoCLUSTER)', c.comparison.map(r => [r[0], r[2]]), col, { symbol: 'diamond', symbolSize: 9, legendHoverLink: false }));
    });
    return { grid: { left: 54, right: 14, top: 86, bottom: 44 }, legend: legend(t, { icon: LINE_ICON, data: B.uloop.map(c => c.label) }),
      tooltip: axisTip(t, 'Year', '', 0, 2, false),
      xAxis: axis(t, { type: 'log', name: 'Time (years)', min: 0.8, max: 50, logBase: 10 }),
      yAxis: axis(t, { type: 'value', name: 'Outlet temperature (°C)', nameGap: 38 }),
      series: ser };
  }, () => { let h = '<table><thead><tr><th>Case</th><th>Year</th><th>aharrsim (°C)</th><th>GeoCLUSTER (°C)</th><th>Δ (K)</th></tr></thead><tbody>';
    for (const c of B.uloop) for (const r of c.comparison) h += `<tr><td>${esc(c.label)}</td><td>${r[0]}</td><td>${fmt(r[1], 2)}</td><td>${fmt(r[2], 2)}</td><td>${fmt(r[1] - r[2], 2)}</td></tr>`; return h + '</tbody></table>'; });
  mount('b-uloop-par', (t) => {
    const pts = []; B.uloop.forEach(c => c.comparison.forEach(r => pts.push([r[2], r[1], c.label + ', year ' + r[0]])));
    const lo = Math.floor(Math.min(...pts.map(p => p[0])) / 20) * 20, hi = Math.ceil(Math.max(...pts.map(p => p[0])) / 20) * 20;
    return { grid: { left: 54, right: 14, top: 30, bottom: 44 },
      tooltip: tip(t, { trigger: 'item', formatter: (p) => `<b style="font-family:${t.mono}">${fmt(p.value[1], 2)} °C</b> aharrsim<br><b style="font-family:${t.mono}">${fmt(p.value[0], 2)} °C</b> GeoCLUSTER<br><span style="color:${t.ink2}">${esc(p.value[2])}</span>` }),
      xAxis: axis(t, { type: 'value', name: 'GeoCLUSTER outlet (°C)', min: lo, max: hi }),
      yAxis: axis(t, { type: 'value', name: 'aharrsim outlet (°C)', nameGap: 40, min: lo, max: hi }),
      series: [lineS('1:1', [[lo, lo], [hi, hi]], t.ref, { lineStyle: { width: 1, color: t.ref, type: 'dashed' } }), ptsS('points', pts, t.s[0], { symbolSize: 8 })] };
  }, () => { let h = '<table><thead><tr><th>Case, year</th><th>GeoCLUSTER</th><th>aharrsim</th></tr></thead><tbody>';
    for (const c of B.uloop) for (const r of c.comparison) h += `<tr><td>${esc(c.label)}, ${r[0]}</td><td>${fmt(r[2], 2)}</td><td>${fmt(r[1], 2)}</td></tr>`; return h + '</tbody></table>'; });

  mount('b-doublet', (t) => {
    const R = B.doublet.rows;
    return { grid: { left: 54, right: 14, top: 56, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: axisTip(t, 'Year', '', 0, 3, false),
      xAxis: axis(t, { type: 'value', name: 'Time (years)', min: 0, max: 100 }),
      yAxis: axis(t, { type: 'value', name: 'Dimensionless cooling of the producer', nameGap: 38, min: 0, max: 0.7 }),
      series: [lineS('model, no loss', R.map(r => [r.years, r.TD_no_loss]), t.s[0], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('arrived streamline fraction (exact)', R.map(r => [r.years, r.arrived_fraction]), t.ref, { showSymbol: true, symbol: 'diamond', symbolSize: 8, lineStyle: { width: 1.5, color: t.ref, type: 'dashed' } }),
               lineS('model, with caprock loss', R.map(r => [r.years, r.TD_with_caprock]), t.s[1], { showSymbol: true, symbol: 'circle', symbolSize: 7 })] };
  }, () => { let h = '<table><thead><tr><th>Year</th><th>T_D no loss</th><th>Arrived fraction</th><th>T_D with caprock</th><th>T produced (°C)</th></tr></thead><tbody>';
    for (const r of B.doublet.rows) h += `<tr><td>${r.years}</td><td>${fmt(r.TD_no_loss, 4)}</td><td>${fmt(r.arrived_fraction, 4)}</td><td>${fmt(r.TD_with_caprock, 4)}</td><td>${fmt(r.T_prod_C, 2)}</td></tr>`; return h + '</tbody></table>'; });

  mount('b-lauw', (t) => {
    const L = B.lauwerier;
    return { grid: { left: 54, right: 14, top: 56, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: axisTip(t, 't/τ', '', 0, 2, false),
      xAxis: axis(t, { type: 'log', name: 'Time / fracture transit time', min: 1.5, max: 100, logBase: 10 }),
      yAxis: axis(t, { type: 'value', name: 'Produced temperature (°C)', nameGap: 40, min: 190, max: 201 }),
      series: [lineS('cyclic-fracture model', zip(L.t_over_tau, L.model_C), t.s[0], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('Lauwerier', zip(L.t_over_tau, L.lauwerier_C), t.ref, { showSymbol: true, symbol: 'diamond', symbolSize: 8, lineStyle: { width: 1.5, color: t.ref, type: 'dashed' } })] };
  }, () => { let h = '<table><thead><tr><th>t/τ</th><th>Model (°C)</th><th>Lauwerier (°C)</th></tr></thead><tbody>';
    B.lauwerier.t_over_tau.forEach((x, i) => { h += `<tr><td>${x}</td><td>${fmt(B.lauwerier.model_C[i], 3)}</td><td>${fmt(B.lauwerier.lauwerier_C[i], 3)}</td></tr>`; }); return h + '</tbody></table>'; });

  mount('b-plant', (t) => {
    const P = B.plant;
    return { grid: { left: 50, right: 14, top: 84, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
      tooltip: axisTip(t, '°C', '', 0, 3, false),
      xAxis: axis(t, { type: 'value', name: 'Inlet temperature (°C)', min: 90, max: 310 }),
      yAxis: axis(t, { type: 'value', name: 'Utilisation efficiency', nameGap: 36, min: 0.1, max: 0.65 }),
      series: [lineS('ORC cycle model', P.binary.map(r => [r.T, r.cycle]), t.s[0], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('GEOPHIRES binary, subcritical', P.binary.map(r => [r.T, r.sub]), t.ref, { showSymbol: true, symbol: 'diamond', symbolSize: 7, lineStyle: { width: 1.5, color: t.ref, type: 'dashed' } }),
               lineS('GEOPHIRES binary, supercritical', P.binary.map(r => [r.T, r.sup]), t.q[0], { showSymbol: true, symbol: 'diamond', symbolSize: 7, lineStyle: { width: 1.5, color: t.q[0], type: 'dashed' } }),
               lineS('double flash, 35 °C condenser', P.flash35.map(r => [r.T, r.double]), t.s[1], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('double flash, 40 °C condenser', P.flash40.map(r => [r.T, r.double]), t.s[3], { showSymbol: true, symbol: 'circle', symbolSize: 7 }),
               lineS('GEOPHIRES double flash', P.flash35.map(r => [r.T, r.g_double]), t.s[4], { showSymbol: true, symbol: 'diamond', symbolSize: 7, lineStyle: { width: 1.5, color: t.s[4], type: 'dashed' } })] };
  }, () => { let h = '<table><thead><tr><th>T in (°C)</th><th>ORC cycle</th><th>GEOPHIRES sub</th><th>GEOPHIRES super</th></tr></thead><tbody>';
    for (const r of B.plant.binary) h += `<tr><td>${r.T}</td><td>${fmt(r.cycle, 3)}</td><td>${fmt(r.sub, 3)}</td><td>${fmt(r.sup, 3)}</td></tr>`;
    h += '</tbody></table><table style="margin-top:8px"><thead><tr><th>T in (°C)</th><th>single 35</th><th>double 35</th><th>single 40</th><th>double 40</th><th>GEOPHIRES single</th><th>GEOPHIRES double</th></tr></thead><tbody>';
    B.plant.flash35.forEach((r, i) => { const q = B.plant.flash40[i]; h += `<tr><td>${r.T}</td><td>${fmt(r.single, 3)}</td><td>${fmt(r.double, 3)}</td><td>${fmt(q.single, 3)}</td><td>${fmt(q.double, 3)}</td><td>${fmt(r.g_single, 3)}</td><td>${fmt(r.g_double, 3)}</td></tr>`; });
    return h + '</tbody></table>'; });

  // =====================================================================================================
  // Cases
  // =====================================================================================================
  for (const c of DATA.cases) {
    const s = c.series;
    if (!s) continue;
    const key = c.key;
    const xName = c.x_kind === 'day' ? 'Day of the first year' : 'Time (years)';
    const xAxisOpt = (t) => c.x_kind === 'day' ? axis(t, { type: 'value', name: xName, min: 0, max: 365 })
      : c.x_kind === 'lin' ? axis(t, { type: 'value', name: xName, min: 0, max: c.horizon })
      : axis(t, { type: 'log', name: xName, min: c.t_min || 0.005, max: c.horizon, logBase: 10 });
    const X = c.x_kind === 'day' ? s.day : s.t_yr;
    // temperature chart
    mount('ts-' + key, (t) => {
      const ser = [];
      if (c.id === 12 || c.id === 13) {
        ser.push(lineS('water at the wellhead, daily mean', zip(X, s.T_wellhead_C), t.s[0]));
        ser.push(lineS('leaving the fracture, daily mean', zip(X, s.T_fracture_C), t.s[1]));
        ser.push(lineS('first of the day', zip(X, s.T_first_C), t.q[0], { lineStyle: { width: 1.2, color: t.q[0], type: 'dashed' } }));
        ser.push(lineS('last of the day', zip(X, s.T_last_C), t.q[3], { lineStyle: { width: 1.2, color: t.q[3], type: 'dashed' } }));
      } else if (c.id === 5) {
        ser.push(lineS('loop outlet (heating seasons)', zip(X, s.T_C.map((v, i) => s.on[i] ? v : null)), t.s[0], { connectNulls: false }));
        ser.push(lineS('loop outlet (off season, no flow)', zip(X, s.T_C.map((v, i) => s.on[i] ? null : v)), t.ref, { lineStyle: { width: 1, color: t.ref, type: 'dotted' } }));
      } else {
        ser.push(lineS(c.id === 17 ? 'water leaving the chamber' : 'produced water at the wellhead', zip(X, s.T_C), t.s[0]));
        if (s.T_res_C) ser.push(lineS('reservoir outlet (fracture face)', zip(X, s.T_res_C), t.ref, { lineStyle: { width: 1.4, color: t.ref, type: 'dashed' } }));
        if (s.T_toe_C) ser.push(lineS('leaving the laterals', zip(X, s.T_toe_C), t.ref, { lineStyle: { width: 1.4, color: t.ref, type: 'dashed' } }));
        for (const v of (c.variants || [])) ser.push(lineS(v.label, zip(v.series.t_yr, v.series.T_C), t.s[1], { lineStyle: { width: 1.6, color: t.s[1] } }));
      }
      return { grid: { left: 50, right: 14, top: 36 + 22 * Math.ceil(ser.length / 2), bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
        tooltip: axisTip(t, c.x_kind === 'day' ? 'Day' : 'Year', '', c.x_kind === 'day' ? 0 : 2, 1, false),
        xAxis: xAxisOpt(t), yAxis: axis(t, { type: 'value', name: '°C', nameGap: 36 }), series: ser };
    }, () => { let h = `<table><thead><tr><th>${esc(xName)}</th><th>Water (°C)</th>${s.T_res_C ? '<th>Reservoir (°C)</th>' : ''}${s.T_toe_C ? '<th>Lateral outlet (°C)</th>' : ''}</tr></thead><tbody>`;
      const TT = (c.id === 12 || c.id === 13) ? s.T_wellhead_C : s.T_C;
      X.forEach((x, i) => { h += `<tr><td>${fmt(x, c.x_kind === 'day' ? 0 : 3)}</td><td>${fmt(TT[i], 1)}</td>${s.T_res_C ? `<td>${fmt(s.T_res_C[i], 1)}</td>` : ''}${s.T_toe_C ? `<td>${fmt(s.T_toe_C[i], 1)}</td>` : ''}</tr>`; }); return h + '</tbody></table>'; });
    // power chart
    mount('pw-' + key, (t) => {
      const ser = []; let yname = 'MW';
      if (c.id === 5) { ser.push(lineS('heat from the ground (kW)', zip(X, s.Q_kW.map((v, i) => s.on[i] ? v : null)), t.heat)); yname = 'kW'; }
      else if (c.id === 6) { ser.push(lineS('heat from the ground (kW)', zip(X, s.Q_wh_MW.map(v => v * 1e3)), t.heat)); yname = 'kW'; }
      else if (c.id === 12) { ser.push(lineS('heat, daily average (MWth)', zip(X, s.P_daily_avg_MWth), t.heat)); ser.push(lineS('gross electricity, binary (MWe)', zip(X, s.MWe), t.s[0])); }
      else if (c.id === 13) { ser.push(lineS('field heat, 20 wells (MWth)', zip(X, s.field_MWth), t.heat)); ser.push(lineS('field gross electricity (MWe)', zip(X, s.field_MWe), t.s[0])); }
      else {
        if (s.heat_MWth && s.heat_MWth.some(v => v != null)) ser.push(lineS('heat above return temperature (MWth)', zip(X, s.heat_MWth), t.heat));
        if (s.loss_MW) ser.push(lineS('conductive loss from the chamber (MW)', zip(X, s.loss_MW), t.ref, { lineStyle: { width: 1.2, color: t.ref, type: 'dashed' } }));
        if (s.Q_res_MW) ser.push(lineS('heat from the fractures (MW)', zip(X, s.Q_res_MW), t.q[3], { lineStyle: { width: 1.4, color: t.q[3], type: 'dashed' } }));
        ser.push(lineS(`gross electricity, ${c.plant} (MWe)`, zip(X, s.MWe), t.s[0]));
        if (s.pump_kW && s.pump_kW.some(v => v > 0)) ser.push(lineS('downhole pump (MWe)', zip(X, s.pump_kW.map(v => v / 1e3)), t.muted, { lineStyle: { width: 1.4, color: t.muted, type: 'dotted' } }));
        for (const v of (c.variants || [])) ser.push(lineS(v.label + ', heat (MWth)', zip(v.series.t_yr, v.series.heat_MWth), t.s[1], { lineStyle: { width: 1.6, color: t.s[1], type: 'dashed' } }), lineS(v.label + ', gross MWe', zip(v.series.t_yr, v.series.MWe), t.s[1], { lineStyle: { width: 1.6, color: t.s[1] } }));
      }
      return { grid: { left: 50, right: 14, top: 36 + 22 * Math.ceil(ser.length / 2), bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
        tooltip: axisTip(t, c.x_kind === 'day' ? 'Day' : 'Year', '', c.x_kind === 'day' ? 0 : 2, 2, false),
        xAxis: xAxisOpt(t), yAxis: axis(t, { type: 'value', name: yname, nameGap: 36, min: 0 }), series: ser };
    }, () => { let h = `<table><thead><tr><th>${esc(xName)}</th><th>Heat</th><th>Gross MWe</th><th>Pump kW</th></tr></thead><tbody>`;
      const Hh = s.heat_MWth || s.P_daily_avg_MWth || s.Q_kW || s.Q_wh_MW || [], E = s.MWe || [], Pk = s.pump_kW || [];
      X.forEach((x, i) => { h += `<tr><td>${fmt(x, c.x_kind === 'day' ? 0 : 3)}</td><td>${fmt(Hh[i], 3)}</td><td>${fmt(E[i], 3)}</td><td>${fmt(Pk[i], 0)}</td></tr>`; }); return h + '</tbody></table>'; });
    // steam fraction chart (two-phase producers)
    if (document.getElementById('xs-' + key)) {
      mount('xs-' + key, (t) => ({ grid: { left: 50, right: 14, top: 40, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
        tooltip: axisTip(t, 'Year', '', 2, 3, false), xAxis: xAxisOpt(t), yAxis: axis(t, { type: 'value', name: 'steam mass fraction', nameGap: 40, min: 0 }),
        series: [lineS('steam fraction at the wellhead', zip(X, s.x), t.heat)] }),
        () => { let h = '<table><thead><tr><th>Year</th><th>Steam fraction</th></tr></thead><tbody>'; X.forEach((x, i) => { h += `<tr><td>${fmt(x, 3)}</td><td>${fmt(s.x[i], 3)}</td></tr>`; }); return h + '</tbody></table>'; });
    }
    // profile chart
    const pr = c.profile;
    if (pr && document.getElementById('pf-' + key)) {
      mount('pf-' + key, (t) => {
        const ser = [lineS('undisturbed rock', zip(pr.T0, pr.md), t.ref, { lineStyle: { width: 1.4, color: t.ref, type: 'dashed' } })];
        if (c.well_kind === 'coaxial') { ser.push(lineS('annulus (down)', zip(pr.T_a, pr.md), t.s[0])); ser.push(lineS('tubing (up)', zip(pr.T_t, pr.md), t.heat)); }
        else if (c.well_kind === 'utube') { ser.push(lineS('down leg', zip(pr.T_a, pr.md), t.s[0])); ser.push(lineS('up leg', zip(pr.T_t, pr.md), t.heat)); }
        else ser.push(lineS(c.well_kind === 'injection' ? 'injected water' : 'produced stream', zip(pr.T_a, pr.md), c.well_kind === 'injection' ? t.s[0] : t.heat));
        if (c.profile_inj) ser.push(lineS('injector', zip(c.profile_inj.T_a, c.profile_inj.md), t.s[0]));
        return { grid: { left: 54, right: 14, top: 60, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
          tooltip: axisTip(t, '°C', 'MD (m)', 1, 0, true),
          xAxis: axis(t, { type: 'value', name: 'Temperature (°C)' }),
          yAxis: axis(t, { type: 'value', name: 'Measured depth (m)', nameGap: 44, inverse: true, min: 0 }), series: ser };
      }, () => { let h = '<table><thead><tr><th>MD (m)</th><th>Rock (°C)</th><th>Stream A (°C)</th><th>Stream B (°C)</th></tr></thead><tbody>';
        pr.md.forEach((z, i) => { h += `<tr><td>${fmt(z, 0)}</td><td>${fmt(pr.T0[i], 1)}</td><td>${fmt(pr.T_a[i], 1)}</td><td>${fmt(pr.T_t[i], 1)}</td></tr>`; }); return h + '</tbody></table>'; });
    }
    if (c.id === 13 && document.getElementById('st-' + key)) {
      mount('st-' + key, (t) => ({ grid: { left: 50, right: 14, top: 40, bottom: 44 }, legend: legend(t, { icon: LINE_ICON }),
        tooltip: axisTip(t, 'fraction produced', '', 2, 1, false),
        xAxis: axis(t, { type: 'value', name: 'Fraction of the 200 m³ produced', min: 0, max: 1 }),
        yAxis: axis(t, { type: 'value', name: '°C', nameGap: 36 }),
        series: [lineS('discharge after a 150 °C charge and 24 h soak', zip(s.storage_discharge_frac, s.storage_discharge_T_C), t.heat)] }),
        () => { let h = '<table><thead><tr><th>Fraction</th><th>T (°C)</th></tr></thead><tbody>'; s.storage_discharge_frac.forEach((f, i) => { h += `<tr><td>${fmt(f, 3)}</td><td>${fmt(s.storage_discharge_T_C[i], 1)}</td></tr>`; }); return h + '</tbody></table>'; });
    }
  }

  // =====================================================================================================
  // Overview
  // =====================================================================================================
  const O = DATA.overview;
  mount('o-scatter', (t) => {
    const groups = [...new Set(O.points.map(p => p.group))];
    return { grid: { left: 56, right: 14, top: 64, bottom: 44 }, legend: legend(t, { icon: 'circle' }),
      tooltip: tip(t, { trigger: 'item', formatter: (p) => `<b style="font-family:${t.mono}">${fmt(p.value[1] > 0.011 ? p.value[1] : 0, 2)} MWe</b> net at ${fmt(p.value[0], 0)} °C<br><span style="color:${t.ink2}">${esc(p.value[2])}</span>` }),
      xAxis: axis(t, { type: 'value', name: 'Water temperature at the horizon (°C)', min: 0, max: 350 }),
      yAxis: axis(t, { type: 'log', name: 'Net electricity (MWe)', nameGap: 44, logBase: 10, min: 0.01, max: 100, axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: (v) => String(v) } }),
      series: groups.map((g, i) => { const col = [t.s[0], t.s[1], t.s[2], t.s[3], t.s[4], t.s[5], t.s[6]][i % 7];
        return ptsS(g, O.points.filter(p => p.group === g).map(p => ({ value: [p.T, Math.max(p.net, 0.011), p.label], symbol: p.net > 0 ? 'circle' : 'emptyCircle' })), col, { symbolSize: 10,
          label: { show: true, position: 'right', formatter: (q) => q.value[2].split(' ')[0], color: t.muted, fontFamily: t.mono, fontSize: 10 } }); }) };
  }, () => { let h = '<table><thead><tr><th>Example</th><th>Water T (°C)</th><th>Net MWe</th></tr></thead><tbody>';
    for (const p of O.points) h += `<tr><td>${esc(p.label)}</td><td>${fmt(p.T, 0)}</td><td>${fmt(p.net, 2)}</td></tr>`; return h + '</tbody></table>'; });

  const Bd = O.dots;
  mount('o-dots', (t) => ({ grid: { left: 56, right: 14, top: 44, bottom: 40 }, legend: legend(t, { data: ['Heat (MWth)', 'Gross electricity (MWe)', 'Net electricity (MWe)'] }),
    tooltip: tip(t, { trigger: 'axis', formatter: (ps) => { const i = ps[0].dataIndex; return `<b>${esc(Bd.labels[i])}</b> ${esc(Bd.names[i])}<br>` + ps.filter(p => p.seriesName !== 'range' && p.value[1] != null).map(p => `${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value[1], 2)}</b>`).join('<br>'); } }),
    xAxis: axis(t, { type: 'category', data: Bd.labels, name: 'Type', nameGap: 24, splitLine: { show: false } }),
    yAxis: axis(t, { type: 'log', name: 'MW', nameGap: 44, logBase: 10, min: 0.005, axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: (v) => String(v) } }),
    series: [{ name: 'Heat (MWth)', type: 'scatter', symbol: 'diamond', symbolSize: 12, data: Bd.heat.map((v, i) => [i, v > 0.001 ? v : null]), itemStyle: { color: t.heat } },
             { name: 'Gross electricity (MWe)', type: 'scatter', symbol: 'emptyCircle', symbolSize: 11, data: Bd.gross.map((v, i) => [i, v > 0.005 ? v : null]), itemStyle: { color: t.s[0], borderWidth: 2 } },
             { name: 'Net electricity (MWe)', type: 'scatter', symbol: 'circle', symbolSize: 8, data: Bd.net.map((v, i) => [i, v > 0.005 ? v : null]), itemStyle: { color: t.s[0] } },
             { name: 'range', type: 'custom', silent: true, tooltip: { show: false },
               renderItem: (params, api) => { const i = api.value(0), lo = api.value(1), hi = api.value(2); if (!(lo > 0) || !(hi > 0)) return null;
                 const a = api.coord([i, lo]), b = api.coord([i, hi]); return { type: 'line', shape: { x1: a[0], y1: a[1], x2: b[0], y2: b[1] }, style: { stroke: t.axis, lineWidth: 1 } }; },
               data: Bd.labels.map((l, i) => [i, Math.min(...[Bd.net[i], Bd.gross[i], Bd.heat[i]].filter(v => v > 0.005)), Math.max(Bd.heat[i], Bd.gross[i])]) }] }),
    () => { let h = '<table><thead><tr><th>#</th><th>Type</th><th>MWth</th><th>Gross MWe</th><th>Net MWe</th></tr></thead><tbody>';
      Bd.labels.forEach((l, i) => { h += `<tr><td>${esc(l)}</td><td>${esc(Bd.names[i])}</td><td>${fmt(Bd.heat[i], 2)}</td><td>${fmt(Bd.gross[i], 2)}</td><td>${fmt(Bd.net[i], 2)}</td></tr>`; }); return h + '</tbody></table>'; });

  const Cl = O.claims;
  mount('o-claims', (t) => ({ grid: { left: 140, right: 20, top: 44, bottom: 44 }, legend: legend(t, { data: ['quoted range', 'model, year 1', 'model, horizon'] }),
    tooltip: tip(t, { trigger: 'item', formatter: (p) => p.seriesName === 'quoted range' ? `<b>${esc(Cl[p.dataIndex].label)}</b><br>quoted ${esc(Cl[p.dataIndex].quote)}` : `<b>${esc(Cl[p.dataIndex].label)}</b><br>${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value[0], 2)} MWe</b>` }),
    yAxis: axis(t, { type: 'category', data: Cl.map(c => c.label), splitLine: { show: false }, axisLabel: { color: t.ink2, fontFamily: t.body, fontSize: 11 } }),
    xAxis: axis(t, { type: 'log', name: 'Gross electricity (MWe), log scale', logBase: 10, min: 0.01, max: 200, axisLabel: { color: t.muted, fontFamily: t.mono, fontSize: 11, formatter: (v) => String(v) } }),
    series: [{ name: 'quoted range', type: 'custom', renderItem: (params, api) => { const i = api.value(0); const lo = Math.max(Cl[i].lo, 0.01), hi = Cl[i].hi;
        const a = api.coord([lo, i]), b = api.coord([hi, i]); const w = Math.max(b[0] - a[0], 6); return { type: 'rect', shape: { x: a[0] - (b[0] - a[0] < 6 ? 3 : 0), y: a[1] - 7, width: w, height: 14, r: 3 }, style: { fill: t.axis, opacity: 0.7 } }; },
        data: Cl.map((c, i) => [i, c.lo, c.hi]), itemStyle: { color: t.axis }, z: 1 },
      { name: 'model, year 1', type: 'scatter', symbol: 'emptyCircle', symbolSize: 11, data: Cl.map((c, i) => [Math.max(c.first, 0.01), i]), itemStyle: { color: t.s[0], borderWidth: 2 }, z: 3 },
      { name: 'model, horizon', type: 'scatter', symbol: 'circle', symbolSize: 9, data: Cl.map((c, i) => [Math.max(c.last, 0.01), i]), itemStyle: { color: t.heat, borderColor: t.surface, borderWidth: 2 }, z: 4 }] }),
    () => { let h = '<table><thead><tr><th>Type</th><th>Quoted</th><th>Model year 1 (MWe)</th><th>Model horizon (MWe)</th></tr></thead><tbody>';
      for (const c of Cl) h += `<tr><td>${esc(c.label)}</td><td>${esc(c.quote)}</td><td>${fmt(c.first, 2)}</td><td>${fmt(c.last, 2)}</td></tr>`; return h + '</tbody></table>'; });

  const Dc = O.decline;
  mount('o-decline', (t) => ({ grid: { left: 140, right: 20, top: 44, bottom: 44 }, legend: legend(t, { data: ['year 1', 'horizon'] }),
    tooltip: tip(t, { trigger: 'axis', formatter: (ps) => `<b>${esc(Dc[ps[0].dataIndex].label)}</b><br>` + ps.filter(p => p.seriesName !== 'link').map(p => `${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value[0], 0)} °C</b>`).join('<br>') }),
    yAxis: axis(t, { type: 'category', data: Dc.map(c => c.label), splitLine: { show: false }, axisLabel: { color: t.ink2, fontFamily: t.body, fontSize: 11 } }),
    xAxis: axis(t, { type: 'value', name: 'Produced water temperature (°C)', min: 0, max: 350 }),
    series: [{ name: 'link', type: 'custom', silent: true, tooltip: { show: false }, renderItem: (params, api) => { const i = api.value(0); const a = api.coord([Dc[i].first, i]), b = api.coord([Dc[i].last, i]); return { type: 'line', shape: { x1: a[0], y1: a[1], x2: b[0], y2: b[1] }, style: { stroke: t.axis, lineWidth: 2 } }; }, data: Dc.map((c, i) => [i]) },
      { name: 'year 1', type: 'scatter', symbol: 'emptyCircle', symbolSize: 10, data: Dc.map((c, i) => [c.first, i]), itemStyle: { color: t.s[0], borderWidth: 2 }, z: 3 },
      { name: 'horizon', type: 'scatter', symbol: 'circle', symbolSize: 9, data: Dc.map((c, i) => [c.last, i]), itemStyle: { color: t.heat, borderColor: t.surface, borderWidth: 2 }, z: 4 }] }),
    () => { let h = '<table><thead><tr><th>Type</th><th>Year 1 (°C)</th><th>Horizon (°C)</th></tr></thead><tbody>'; for (const c of Dc) h += `<tr><td>${esc(c.label)}</td><td>${fmt(c.first, 1)}</td><td>${fmt(c.last, 1)}</td></tr>`; return h + '</tbody></table>'; });

  const Am = O.ambient;
  mount('o-ambient', (t) => ({ grid: { left: 56, right: 14, top: 64, bottom: 60 }, legend: legend(t, {}),
    tooltip: tip(t, { trigger: 'axis', formatter: (ps) => `<b>${esc(Am.names[ps[0].dataIndex])}</b><br>` + ps.map(p => `${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value, 0)} %</b>`).join('<br>') }),
    xAxis: axis(t, { type: 'category', data: Am.labels, name: 'Type', nameGap: 24, splitLine: { show: false } }),
    yAxis: axis(t, { type: 'value', name: 'Output relative to the reference (%)', nameGap: 40, min: 0, max: 110 }),
    series: Am.series.map((s, i) => ({ name: s.name, type: 'bar', data: s.values, itemStyle: { color: [t.q[1], t.q[3], t.q[4], t.s[1]][i], borderRadius: [3, 3, 0, 0] }, barGap: '10%', barCategoryGap: '30%' })) }),
    () => { let h = `<table><thead><tr><th>Type</th>${Am.series.map(s => `<th>${esc(s.name)}</th>`).join('')}</tr></thead><tbody>`; Am.labels.forEach((l, i) => { h += `<tr><td>${esc(l)} ${esc(Am.names[i])}</td>${Am.series.map(s => `<td>${fmt(s.values[i], 0)}</td>`).join('')}</tr>`; }); return h + '</tbody></table>'; });

  const Pa = O.parasitic;
  mount('o-para', (t) => ({ grid: { left: 56, right: 14, top: 44, bottom: 60 }, legend: legend(t, {}),
    tooltip: tip(t, { trigger: 'axis', formatter: (ps) => `<b>${esc(Pa.names[ps[0].dataIndex])}</b><br>` + ps.map(p => `${esc(p.seriesName)}: <b style="font-family:${t.mono}">${fmt(p.value, 2)} MWe</b>`).join('<br>') }),
    xAxis: axis(t, { type: 'category', data: Pa.labels, name: 'Type', nameGap: 24, splitLine: { show: false } }),
    yAxis: axis(t, { type: 'value', name: 'MWe', nameGap: 40, min: 0 }),
    series: [{ name: 'net', type: 'bar', stack: 'a', data: Pa.net, itemStyle: { color: t.s[0] } }, { name: 'downhole pump', type: 'bar', stack: 'a', data: Pa.pump, itemStyle: { color: t.q[0] } }, { name: 'circulation against impedance', type: 'bar', stack: 'a', data: Pa.circ, itemStyle: { color: t.axis, borderRadius: [3, 3, 0, 0] } }] }),
    () => { let h = '<table><thead><tr><th>Type</th><th>Net</th><th>Downhole pump</th><th>Circulation</th><th>Share of gross</th></tr></thead><tbody>'; Pa.labels.forEach((l, i) => { h += `<tr><td>${esc(l)} ${esc(Pa.names[i])}</td><td>${fmt(Pa.net[i], 2)}</td><td>${fmt(Pa.pump[i], 2)}</td><td>${fmt(Pa.circ[i], 2)}</td><td>${fmt(Pa.share[i], 0)} %</td></tr>`; }); return h + '</tbody></table>'; });

  // ---------------- tables, resize, theme ----------------
  document.querySelectorAll('[data-table]').forEach(b => {
    const id = b.dataset.table;
    b.setAttribute('aria-expanded', 'false');
    b.addEventListener('click', () => {
      let box = document.getElementById(id + '-tbl');
      if (!box) { box = document.createElement('div'); box.id = id + '-tbl'; box.className = 'dtable'; box.hidden = true; b.closest('figure').appendChild(box); }
      box.hidden = !box.hidden;
      b.setAttribute('aria-expanded', String(!box.hidden));
      b.textContent = box.hidden ? 'Table' : 'Hide table';
      if (!box.hidden && tables[id]) box.innerHTML = tables[id]();
    });
  });
  let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => charts.forEach(x => x.c.resize()), 120); });
  const rerender = () => charts.forEach(x => x.render());
  try { window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', rerender); } catch (e) {}
  new MutationObserver(rerender).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
})();
