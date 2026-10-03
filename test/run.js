// Run with:  npm test   (plain Node, no dependencies)
import assert from 'node:assert/strict';
import { solveChain, applyH, bondDiagonal, buildDenseH } from '../src/physics.js';
import { jacobiEigen } from '../src/linalg.js';

let passed = 0;
function test(name, fn) {
  const t0 = performance.now();
  try {
    fn();
    passed++;
    console.log(`  ok   ${name}  (${(performance.now() - t0).toFixed(0)} ms)`);
  } catch (err) {
    console.error(`  FAIL ${name}\n${err.stack}`);
    process.exitCode = 1;
  }
}
const close = (a, b, tol, msg) =>
  assert.ok(Math.abs(a - b) < tol, `${msg ?? ''} expected ${b}, got ${a} (tol ${tol})`);

function dense(N, J, h, periodic) {
  const dim = 1 << N;
  const { values, vectors } = jacobiEigen(buildDenseH(N, J, h, periodic), dim);
  const mask = dim - 1;
  const parity = vectors.map((v) => {
    let p = 0;
    for (let s = 0; s < dim; s++) p += v[s] * v[mask ^ s];
    return p;
  });
  return { values, vectors, parity };
}

console.log('Transverse-field Ising chain tests');

test('analytic limit h = 0 (classical Ising)', () => {
  for (const N of [4, 8, 12]) {
    close(solveChain({ N, h: 0, periodic: true, sites: false }).E0, -N, 1e-9, `periodic N=${N}`);
    close(solveChain({ N, h: 0, periodic: false, sites: false }).E0, -(N - 1), 1e-9, `open N=${N}`);
  }
});

test('analytic limit J = 0 (free spins in a transverse field)', () => {
  const N = 10, h = 0.7;
  const r = solveChain({ N, J: 0, h, sites: false });
  close(r.E0, -h * N, 1e-9, 'E0');
  close(r.gap, 2 * h, 1e-8, 'gap');
  close(r.sx, 1, 1e-8, '<sx>');
});

test('Lanczos matches dense diagonalization (E0, odd-sector energy)', () => {
  const cases = [
    [4, true, [0.3, 0.8, 1.0, 1.4, 2.0]],
    [4, false, [0.3, 1.0, 2.0]],
    [6, true, [0.5, 1.0, 1.5]],
    [6, false, [0.5, 1.0, 1.5]],
    [7, true, [0.7, 1.2]],
  ];
  for (const [N, periodic, hs] of cases) {
    for (const h of hs) {
      const d = dense(N, 1, h, periodic);
      const r = solveChain({ N, h, periodic, sites: false });
      close(r.E0, d.values[0], 1e-9, `E0 N=${N} h=${h} periodic=${periodic}`);
      const iOdd = d.parity.findIndex((p) => p < 0);
      close(r.Eodd, d.values[iOdd], 1e-7, `E_odd N=${N} h=${h} periodic=${periodic}`);
    }
  }
});

test('observables match dense ground state (m², <σˣ>, per-site <σᶻ>)', () => {
  const N = 6;
  for (const periodic of [true, false]) {
    for (const h of [0.6, 1.0, 1.6]) {
      const d = dense(N, 1, h, periodic);
      const dim = 1 << N;
      const g = d.vectors[0];
      const iOdd = d.parity.findIndex((p) => p < 0);
      const o = d.vectors[iOdd];
      let m2 = 0, sx = 0;
      const sz = new Array(N).fill(0);
      for (let s = 0; s < dim; s++) {
        let M = 0;
        for (let i = 0; i < N; i++) M += 1 - 2 * ((s >> i) & 1);
        m2 += g[s] * g[s] * M * M;
        for (let i = 0; i < N; i++) {
          sx += g[s] * g[s ^ (1 << i)];
          sz[i] += g[s] * o[s] * (1 - 2 * ((s >> i) & 1));
        }
      }
      const sgn = sz.reduce((a, b) => a + b, 0) < 0 ? -1 : 1;
      const r = solveChain({ N, h, periodic, sites: true });
      close(r.m2, m2 / (N * N), 1e-6, `m2 h=${h}`);
      close(r.sx, sx / N, 1e-6, `sx h=${h}`);
      for (let i = 0; i < N; i++) close(r.sz[i], sgn * sz[i], 1e-5, `sz[${i}] h=${h} periodic=${periodic}`);
    }
  }
});

test('eigen-residual ||Hψ − Eψ|| is small at N = 12', () => {
  const N = 12;
  const diag = bondDiagonal(N, true);
  for (const h of [0.4, 1.0, 1.7]) {
    const r = solveChain({ N, h, periodic: true, sites: false, keepVectors: true });
    for (const [name, psi, E] of [['even', r.vectors.even, r.E0], ['odd', r.vectors.odd, r.Eodd]]) {
      const out = new Float64Array(psi.length);
      applyH(psi, out, N, 1, h, diag);
      let res = 0;
      for (let s = 0; s < psi.length; s++) res += (out[s] - E * psi[s]) ** 2;
      assert.ok(Math.sqrt(res) < 1e-5, `${name} sector residual ${Math.sqrt(res)} at h=${h}`);
    }
  }
});

test('periodic chain is translation invariant', () => {
  const r = solveChain({ N: 10, h: 0.5, periodic: true, sites: true });
  for (let i = 1; i < 10; i++) {
    close(r.sz[i], r.sz[0], 1e-6, 'sz');
    close(r.sxSites[i], r.sxSites[0], 1e-6, 'sx');
  }
});

test('phase-transition signatures at N = 12', () => {
  const lo = solveChain({ N: 12, h: 0.2, sites: false });
  const mid = solveChain({ N: 12, h: 1.0, sites: false });
  const hi = solveChain({ N: 12, h: 2.0, sites: false });
  assert.ok(lo.m2 > 0.9, `ordered phase m² = ${lo.m2}`);
  assert.ok(hi.m2 < 0.2, `disordered phase m² = ${hi.m2}`);
  assert.ok(lo.m2 > mid.m2 && mid.m2 > hi.m2, 'm² decreases with h');
  assert.ok(lo.gap < 1e-3, `gap closes in ordered phase: ${lo.gap}`);
  assert.ok(hi.gap > mid.gap && mid.gap > lo.gap, 'gap grows with h');
  assert.ok(lo.sx < mid.sx && mid.sx < hi.sx, '<σˣ> grows with h');
});

console.log(`${passed} test(s) passed`);
