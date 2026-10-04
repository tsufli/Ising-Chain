// UI: sliders -> worker -> arrow chain + three plots. No dependencies.

const $ = (id) => document.getElementById(id);
const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });

const H_GRID = Array.from({ length: 41 }, (_, i) => +(i * 0.05).toFixed(2)); // 0 … 2
const state = { N: 12, h: 0.5, periodic: true };

// Sweep order: coarse first (ends, middle, quarters, …) so a rough curve appears
// almost immediately and then fills in.
const SWEEP_ORDER = (() => {
  const seen = new Set();
  const order = [];
  for (const step of [40, 20, 10, 5, 1]) {
    for (let i = 0; i < H_GRID.length; i += step) if (!seen.has(i)) { seen.add(i); order.push(i); }
  }
  return order;
})();

// ---------- worker plumbing ----------------------------------------------

const sweeps = new Map();       // "12p" -> { data: [sweepPoint | null …], count, done }
const sweepOwner = new Map();   // request id -> cache key
let sweepSeq = 0;
let solveSeq = 0;
let solveBusy = false;
let solvePending = false;
let live = null;
let liveMs = 0;
let sweepTimer = null;

const keyOf = (N, periodic) => `${N}${periodic ? 'p' : 'o'}`;
const currentKey = () => keyOf(state.N, state.periodic);

function requestSolve() {
  if (solveBusy) { solvePending = true; return; }
  solveBusy = true;
  solvePending = false;
  worker.postMessage({ type: 'solve', id: ++solveSeq, N: state.N, h: state.h, periodic: state.periodic });
}

function requestSweep() {
  const key = currentKey();
  const hit = sweeps.get(key);
  if (hit && hit.done) return;
  sweeps.set(key, { data: new Array(H_GRID.length).fill(null), count: 0, done: false });
  const id = ++sweepSeq;
  sweepOwner.set(id, key);
  worker.postMessage({
    type: 'sweep', id, N: state.N, periodic: state.periodic,
    hs: SWEEP_ORDER.map((i) => H_GRID[i]), idx: SWEEP_ORDER,
  });
}

worker.onmessage = (e) => {
  const msg = e.data;
  if (msg.type === 'solve') {
    solveBusy = false;
    if (msg.res.N === state.N && msg.res.periodic === state.periodic) {
      live = msg.res;
      liveMs = msg.ms;
      render();
    }
    if (solvePending) requestSolve();
  } else if (msg.type === 'sweep') {
    const entry = sweeps.get(sweepOwner.get(msg.id));
    if (!entry) return;
    entry.data[msg.i] = msg.res;
    entry.count++;
    if (sweepOwner.get(msg.id) === currentKey()) drawCharts();
    updateStatus();
  } else if (msg.type === 'sweepDone') {
    const entry = sweeps.get(sweepOwner.get(msg.id));
    if (entry) entry.done = true;
    updateStatus();
  } else if (msg.type === 'error') {
    $('status').textContent = `error: ${msg.message}`;
  }
};
worker.onerror = (e) => { $('status').textContent = `worker error: ${e.message}`; };

// ---------- canvas helpers -------------------------------------------------

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

function prepare(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  const pw = Math.round(w * dpr), ph = Math.round(h * dpr);
  if (canvas.width !== pw || canvas.height !== ph) { canvas.width = pw; canvas.height = ph; }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function hexToRgb(hex) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex);
  if (!m) return [128, 128, 128];
  const n = parseInt(m[1], 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}
const mix = (a, b, t) => `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(',')})`;

function arrow(ctx, x0, y0, x1, y1, color, lw) {
  const dx = x1 - x0, dy = y1 - y0;
  const len = Math.hypot(dx, dy);
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = lw; ctx.lineCap = 'round';
  if (len < 2) {
    ctx.beginPath(); ctx.arc((x0 + x1) / 2, (y0 + y1) / 2, lw * 0.9, 0, 2 * Math.PI); ctx.fill();
    return;
  }
  const ux = dx / len, uy = dy / len;
  const head = Math.min(11, len * 0.45);
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1 - ux * head * 0.6, y1 - uy * head * 0.6); ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(x1 - ux * head - uy * head * 0.5, y1 - uy * head + ux * head * 0.5);
  ctx.lineTo(x1 - ux * head + uy * head * 0.5, y1 - uy * head - ux * head * 0.5);
  ctx.closePath(); ctx.fill();
}

// ---------- chain of arrows ------------------------------------------------

function drawChain() {
  const canvas = $('chain');
  const { ctx, w, h } = prepare(canvas);
  const text = css('--muted');
  const up = hexToRgb(css('--up'));
  const right = hexToRgb(css('--right'));
  const N = state.N;

  const padX = 28;
  const spacing = (w - 2 * padX) / N;
  const cy = h / 2;
  const maxLen = Math.min(spacing * 0.95, h * 0.8);

  // chain backbone
  ctx.strokeStyle = css('--line'); ctx.lineWidth = 2; ctx.setLineDash([]);
  ctx.beginPath(); ctx.moveTo(padX, cy); ctx.lineTo(w - padX, cy); ctx.stroke();
  if (state.periodic) {
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(padX, cy); ctx.lineTo(padX - 14, cy);
    ctx.moveTo(w - padX, cy); ctx.lineTo(w - padX + 14, cy);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  const ok = live && live.N === N && live.periodic === state.periodic;
  for (let i = 0; i < N; i++) {
    const cx = padX + spacing * (i + 0.5);
    ctx.fillStyle = css('--line');
    ctx.beginPath(); ctx.arc(cx, cy, 3, 0, 2 * Math.PI); ctx.fill();
    if (!ok) continue;
    const mz = live.sz[i], mx = live.sxSites[i];
    const t = Math.min(1, Math.max(0, Math.atan2(mx, mz) / (Math.PI / 2)));
    const color = mix(up, right, t);
    const L = maxLen;
    const dx = mx * L, dy = -mz * L;
    arrow(ctx, cx - dx / 2, cy - dy / 2, cx + dx / 2, cy + dy / 2, color, Math.max(2, Math.min(4, spacing * 0.08)));
  }

  ctx.fillStyle = text; ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('↑ ⟨σᶻ⟩    → ⟨σˣ⟩', 8, 16);
  ctx.textAlign = 'right';
  ctx.fillText(state.periodic ? 'periodic ring' : 'open chain', w - 8, 16);
}

// ---------- shared plot frame ----------------------------------------------

const PAD = { L: 46, R: 14, T: 30, B: 32 };

/** Draws title, grid, tick labels; returns coordinate maps X(), Y(). */
function frame(canvas, { title, y, x, xLabel }) {
  const { ctx, w, h } = prepare(canvas);
  const { L, R, T, B } = PAD;
  const X = (v) => L + ((v - x.lo) / (x.hi - x.lo)) * (w - L - R);
  const Y = (v) => T + (1 - (v - y.lo) / (y.hi - y.lo)) * (h - T - B);

  ctx.fillStyle = css('--text'); ctx.textAlign = 'left';
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillText(title, 12, 16);

  ctx.font = '11px system-ui, sans-serif';
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.fillStyle = css('--muted');
  ctx.textAlign = 'right';
  for (const v of y.ticks) {
    const py = Math.round(Y(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(L, py); ctx.lineTo(w - R, py); ctx.stroke();
    ctx.fillText(y.fmt(v), L - 6, py + 4);
  }
  ctx.textAlign = 'center';
  for (const v of x.ticks) {
    const px = Math.round(X(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(px, T); ctx.lineTo(px, h - B); ctx.stroke();
    ctx.fillText(String(v), px, h - B + 14);
  }
  ctx.fillText(xLabel, (L + w - R) / 2, h - 4);
  return { ctx, w, h, X, Y };
}

const H_AXIS = { lo: 0, hi: 2, ticks: [0, 0.5, 1, 1.5, 2] };

/** Dashed critical line at h = J and the live slider marker. */
function hMarkers({ ctx, X, h }) {
  const { T, B } = PAD;
  ctx.strokeStyle = css('--crit'); ctx.setLineDash([4, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(X(1), T); ctx.lineTo(X(1), h - B); ctx.stroke();
  ctx.setLineDash([]);
  ctx.strokeStyle = css('--accent'); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(X(state.h), T); ctx.lineTo(X(state.h), h - B); ctx.stroke();
}

function emptyNote({ ctx, w, h }, text) {
  ctx.fillStyle = css('--muted'); ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'center';
  ctx.fillText(text, (PAD.L + w - PAD.R) / 2, h / 2);
}

const liveOk = () => live && live.N === state.N && live.periodic === state.periodic;
const sweepData = () => (sweeps.get(currentKey()) || { data: [] }).data;

/** Polyline through (h, value) points, breaking the line where a value is missing. */
function polyline(ctx, X, Y, pts) {
  ctx.beginPath();
  let pen = false;
  for (const p of pts) {
    if (p == null) { pen = false; continue; }
    if (pen) ctx.lineTo(X(p[0]), Y(p[1])); else ctx.moveTo(X(p[0]), Y(p[1]));
    pen = true;
  }
  ctx.stroke();
}

/** Linear interpolation of sorted (x, y) points at x; null if x is outside the data. */
function interpAt(pts, x) {
  const p = pts.filter((q) => q != null);
  if (p.length && x < p[0][0] && p[0][0] <= 0.05 + 1e-9) return p[0][1]; // curve starts at h = 0.05
  for (let i = 0; i + 1 < p.length; i++) {
    if (x >= p[i][0] - 1e-9 && x <= p[i + 1][0] + 1e-9) {
      const t = p[i + 1][0] === p[i][0] ? 0 : (x - p[i][0]) / (p[i + 1][0] - p[i][0]);
      return p[i][1] + t * (p[i + 1][1] - p[i][1]);
    }
  }
  return null;
}

function niceTicks(max) {
  const rough = Math.max(max, 1e-9) / 4;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= rough);
  const ticks = [];
  for (let v = 0; v < max + step * 0.999; v += step) ticks.push(+v.toFixed(10));
  return ticks;
}

// ---------- plot 1: low-lying spectrum ---------------------------------------

const SPEC_SERIES = [
  { sector: 'even', k: 1, color: '--even', dash: [] },
  { sector: 'even', k: 2, color: '--even', dash: [] },
  { sector: 'odd', k: 0, color: '--odd', dash: [6, 4] },
  { sector: 'odd', k: 1, color: '--odd', dash: [6, 4] },
  { sector: 'odd', k: 2, color: '--odd', dash: [6, 4] },
];

function drawSpectrum() {
  const data = sweepData();
  const series = SPEC_SERIES.map((s) =>
    H_GRID.map((hv, i) => {
      // h = 0 is the degenerate classical limit where each band of near-degenerate
      // levels collapses to one value, so the spectrum curves start at h = 0.05.
      if (i === 0) return null;
      const v = data[i] && data[i][s.sector][s.k];
      return v == null ? null : [hv, Math.max(0, v)];
    }));

  let max = 0;
  for (const pts of series) for (const p of pts) if (p) max = Math.max(max, p[1]);
  const ticks = niceTicks(max || 4);
  const top = ticks[ticks.length - 1];

  const f = frame($('chartSpec'), {
    title: 'Low-lying spectrum  E − E₀',
    y: { lo: -0.05 * top, hi: top, ticks, fmt: (v) => (Number.isInteger(v) ? String(v) : v.toFixed(1)) },
    x: H_AXIS,
    xLabel: 'h / J',
  });
  const { ctx, X, Y, w } = f;

  // ground state = zero line
  ctx.strokeStyle = css('--muted'); ctx.lineWidth = 1.5; ctx.setLineDash([]);
  ctx.beginPath(); ctx.moveTo(X(0), Y(0)); ctx.lineTo(X(2), Y(0)); ctx.stroke();

  hMarkers(f);

  ctx.lineWidth = 2; ctx.lineJoin = 'round';
  SPEC_SERIES.forEach((s, i) => {
    ctx.strokeStyle = css(s.color);
    ctx.setLineDash(s.dash);
    polyline(ctx, X, Y, series[i]);
  });
  ctx.setLineDash([]);

  // dots at the slider position
  SPEC_SERIES.forEach((s, i) => {
    const v = interpAt(series[i], state.h);
    if (v == null) return;
    ctx.fillStyle = css(s.color);
    ctx.beginPath(); ctx.arc(X(state.h), Y(v), 3.5, 0, 2 * Math.PI); ctx.fill();
  });

  // legend
  ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'right'; ctx.fillStyle = css('--muted');
  const lx = w - PAD.R;
  ctx.fillText('odd', lx, 16);
  ctx.strokeStyle = css('--odd'); ctx.lineWidth = 2; ctx.setLineDash([5, 3]);
  ctx.beginPath(); ctx.moveTo(lx - 56, 12); ctx.lineTo(lx - 30, 12); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillText('even', lx - 66, 16);
  ctx.strokeStyle = css('--even');
  ctx.beginPath(); ctx.moveTo(lx - 124, 12); ctx.lineTo(lx - 98, 12); ctx.stroke();

  if (!data.some(Boolean)) emptyNote(f, 'computing…');
}

// ---------- plot 2: order parameter m² ---------------------------------------

function drawOrder() {
  const data = sweepData();
  const pts = H_GRID.map((hv, i) => (data[i] ? [hv, data[i].m2] : null));
  const f = frame($('chartM'), {
    title: 'Order parameter  m² = ⟨(Σσᶻ)²⟩ / N²',
    y: { lo: 0, hi: 1, ticks: [0, 0.25, 0.5, 0.75, 1], fmt: (v) => v.toFixed(2) },
    x: H_AXIS,
    xLabel: 'h / J',
  });
  const { ctx, X, Y } = f;
  hMarkers(f);
  ctx.strokeStyle = css('--curve'); ctx.lineWidth = 2; ctx.lineJoin = 'round'; ctx.setLineDash([]);
  polyline(ctx, X, Y, pts);
  if (liveOk()) {
    ctx.fillStyle = css('--accent');
    ctx.beginPath(); ctx.arc(X(state.h), Y(live.m2), 4.5, 0, 2 * Math.PI); ctx.fill();
  }
  if (!data.some(Boolean)) emptyNote(f, 'computing…');
}

// ---------- plot 3: correlation function C(r) --------------------------------

function drawCorr() {
  const ok = liveOk() && live.corr;
  const rmax = state.N >> 1;
  const C = ok ? Array.from(live.corr).slice(1) : [];

  const minC = C.length ? Math.min(...C) : 0;
  const lo = Math.min(0, Math.floor(minC * 4) / 4);
  const ticks = [];
  for (let v = lo; v <= 1 + 1e-9; v += 0.25) ticks.push(+v.toFixed(2));
  const f = frame($('chartC'), {
    title: 'Correlation  C(r) = ⟨σᶻᵢ σᶻᵢ₊ᵣ⟩',
    y: { lo: lo - 0.06, hi: 1.14, ticks, fmt: (v) => v.toFixed(2) },
    x: { lo: 0.5, hi: rmax + 0.5, ticks: Array.from({ length: rmax }, (_, i) => i + 1) },
    xLabel: 'distance r (sites)',
  });
  const { ctx, X, Y } = f;

  ctx.strokeStyle = css('--muted'); ctx.lineWidth = 1; ctx.setLineDash([]);
  ctx.beginPath(); ctx.moveTo(X(0.5), Y(0)); ctx.lineTo(X(rmax + 0.5), Y(0)); ctx.stroke();

  if (!ok) { emptyNote(f, 'computing…'); return; }

  const accent = css('--accent');
  ctx.strokeStyle = accent; ctx.lineWidth = 1.5; ctx.globalAlpha = 0.35; ctx.lineJoin = 'round';
  ctx.beginPath();
  C.forEach((c, i) => (i ? ctx.lineTo(X(i + 1), Y(c)) : ctx.moveTo(X(i + 1), Y(c))));
  ctx.stroke();
  ctx.globalAlpha = 1;

  ctx.font = '10px system-ui, sans-serif'; ctx.textAlign = 'center';
  C.forEach((c, i) => {
    const px = X(i + 1), py = Y(c);
    ctx.strokeStyle = accent; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(px, Y(0)); ctx.lineTo(px, py); ctx.stroke();
    ctx.fillStyle = accent;
    ctx.beginPath(); ctx.arc(px, py, 4.5, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = css('--muted');
    ctx.fillText(c.toFixed(2), px, c >= 0 ? py - 9 : py + 16);
  });
}

function drawCharts() { drawSpectrum(); drawOrder(); drawCorr(); }

// ---------- status + render --------------------------------------------------

function phaseLabel(h) {
  if (h < 0.9) return 'Ordered phase (h < J): ferromagnet';
  if (h <= 1.1) return 'Near the critical point (h ≈ J)';
  return 'Disordered phase (h > J): paramagnet';
}

function updateStatus() {
  const entry = sweeps.get(currentKey());
  const sweep = !entry ? '' : entry.done ? ' · curves complete' : ` · computing curves ${entry.count}/${H_GRID.length}`;
  const t = live ? `solved in ${liveMs.toFixed(0)} ms` : 'solving…';
  $('status').textContent = `dim = 2^${state.N} = ${(1 << state.N).toLocaleString()} · ${t}${sweep}`;
}

function render() {
  $('phase').textContent = phaseLabel(state.h);
  drawChain();
  drawCharts();
  updateStatus();
}

// ---------- events ----------------------------------------------------------

$('hSlider').addEventListener('input', (e) => {
  state.h = parseFloat(e.target.value);
  $('hOut').textContent = state.h.toFixed(2);
  $('phase').textContent = phaseLabel(state.h);
  requestSolve();
  drawCharts();
});

$('nSlider').addEventListener('input', (e) => {
  state.N = parseInt(e.target.value, 10);
  $('nOut').textContent = String(state.N);
  live = null;
  requestSolve();
  clearTimeout(sweepTimer);
  sweepTimer = setTimeout(requestSweep, 150);
  render();
});

$('bcSelect').addEventListener('change', (e) => {
  state.periodic = e.target.value === 'periodic';
  live = null;
  requestSolve();
  requestSweep();
  render();
});

new ResizeObserver(() => render()).observe(document.querySelector('main'));
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);

requestSolve();
requestSweep();
render();
