// UI: sliders -> worker -> arrows + charts. No dependencies.

const $ = (id) => document.getElementById(id);
const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });

const H_GRID = Array.from({ length: 41 }, (_, i) => +(i * 0.05).toFixed(2)); // 0 … 2
const state = { N: 12, h: 0.5, periodic: true };

// ---------- worker plumbing ----------------------------------------------

const sweeps = new Map();       // "12p" -> { E0perSite:[], m2:[], gap:[], sx:[], count, done }
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
  const entry = {
    E0perSite: new Array(H_GRID.length).fill(null),
    m2: new Array(H_GRID.length).fill(null),
    gap: new Array(H_GRID.length).fill(null),
    sx: new Array(H_GRID.length).fill(null),
    count: 0,
    done: false,
  };
  sweeps.set(key, entry);
  const id = ++sweepSeq;
  sweepOwner.set(id, key);
  worker.postMessage({ type: 'sweep', id, N: state.N, periodic: state.periodic, hs: H_GRID });
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
    const r = msg.res;
    entry.E0perSite[msg.i] = r.E0perSite;
    entry.m2[msg.i] = r.m2;
    entry.gap[msg.i] = r.gap;
    entry.sx[msg.i] = r.sx;
    entry.count = entry.m2.filter((v) => v != null).length;
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

// ---------- charts ---------------------------------------------------------

const CHARTS = [
  { id: 'chartM',   key: 'm2',        title: 'Order parameter  m² = ⟨(Σσᶻ)²⟩ / N²', y: [0, 1],    fmt: (v) => v.toFixed(2) },
  { id: 'chartGap', key: 'gap',       title: 'Energy gap  E₁ − E₀',                  y: [0, null], fmt: (v) => v.toFixed(2) },
  { id: 'chartE',   key: 'E0perSite', title: 'Ground-state energy per site  E₀ / N', y: [null, null], fmt: (v) => v.toFixed(2) },
  { id: 'chartX',   key: 'sx',        title: 'Transverse magnetization  ⟨σˣ⟩',       y: [0, 1],    fmt: (v) => v.toFixed(2) },
];

function drawChart(cfg) {
  const { ctx, w, h } = prepare($(cfg.id));
  const L = 44, R = 14, T = 28, B = 30;
  const entry = sweeps.get(currentKey());
  const ys = entry ? entry[cfg.key] : [];
  const pts = [];
  H_GRID.forEach((x, i) => { if (ys[i] != null && Number.isFinite(ys[i])) pts.push([x, ys[i]]); });

  let lo = cfg.y[0], hi = cfg.y[1];
  const vals = pts.map((p) => p[1]);
  if (live && live.N === state.N) vals.push(live[cfg.key]);
  if (lo == null) lo = vals.length ? Math.min(...vals) : 0;
  if (hi == null) hi = vals.length ? Math.max(...vals) : 1;
  if (cfg.y[0] == null || cfg.y[1] == null) {
    const pad = (hi - lo) * 0.08 || 0.1;
    if (cfg.y[0] == null) lo -= pad;
    if (cfg.y[1] == null) hi += pad;
  }
  if (hi - lo < 1e-9) hi = lo + 1;

  const X = (x) => L + (x / 2) * (w - L - R);
  const Y = (y) => T + (1 - (y - lo) / (hi - lo)) * (h - T - B);

  ctx.font = '11px system-ui, sans-serif';
  ctx.fillStyle = css('--text'); ctx.textAlign = 'left';
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillText(cfg.title, L - 30, 16);
  ctx.font = '11px system-ui, sans-serif';

  // grid + y labels
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.fillStyle = css('--muted'); ctx.textAlign = 'right';
  for (let k = 0; k <= 4; k++) {
    const v = lo + ((hi - lo) * k) / 4;
    const y = Math.round(Y(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(w - R, y); ctx.stroke();
    ctx.fillText(cfg.fmt(v), L - 6, y + 4);
  }
  ctx.textAlign = 'center';
  for (const x of [0, 0.5, 1, 1.5, 2]) {
    const px = Math.round(X(x)) + 0.5;
    ctx.beginPath(); ctx.moveTo(px, T); ctx.lineTo(px, h - B); ctx.stroke();
    ctx.fillText(String(x), px, h - B + 14);
  }
  ctx.fillText('h / J', (L + w - R) / 2, h - 4);

  // critical line h = J
  ctx.strokeStyle = css('--crit'); ctx.setLineDash([4, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(X(1), T); ctx.lineTo(X(1), h - B); ctx.stroke();
  ctx.setLineDash([]);

  // curve
  if (pts.length > 1) {
    ctx.strokeStyle = css('--curve'); ctx.lineWidth = 2; ctx.lineJoin = 'round';
    ctx.beginPath();
    pts.forEach(([x, y], i) => (i ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y))));
    ctx.stroke();
  }

  // slider marker + live dot
  ctx.strokeStyle = css('--accent'); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(X(state.h), T); ctx.lineTo(X(state.h), h - B); ctx.stroke();
  if (live && live.N === state.N && live.periodic === state.periodic) {
    ctx.fillStyle = css('--accent');
    ctx.beginPath(); ctx.arc(X(state.h), Y(live[cfg.key]), 4.5, 0, 2 * Math.PI); ctx.fill();
  }
}

function drawCharts() { CHARTS.forEach(drawChart); }

// ---------- readout + render -----------------------------------------------

function phaseLabel(h) {
  if (h < 0.9) return 'Ordered phase (h < J): ferromagnet';
  if (h <= 1.1) return 'Near the critical point (h ≈ J)';
  return 'Disordered phase (h > J): paramagnet';
}

function updateStatus() {
  const entry = sweeps.get(currentKey());
  const n = H_GRID.length;
  const sweep = !entry ? '' : entry.done ? ' · curves complete' : ` · computing curves ${entry.count}/${n}`;
  const t = live ? `solved in ${liveMs.toFixed(0)} ms` : 'solving…';
  $('status').textContent = `dim = 2^${state.N} = ${(1 << state.N).toLocaleString()} · ${t}${sweep}`;
}

function render() {
  $('phase').textContent = phaseLabel(state.h);
  if (live) {
    $('rE0').textContent = live.E0.toFixed(4);
    $('rEN').textContent = live.E0perSite.toFixed(4);
    $('rGap').textContent = live.gap < 1e-4 ? live.gap.toExponential(1) : live.gap.toFixed(4);
    $('rM2').textContent = live.m2.toFixed(3);
    $('rSx').textContent = live.sx.toFixed(3);
  }
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
