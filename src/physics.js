// Transverse-field Ising chain:  H = -J Σ σᶻᵢσᶻᵢ₊₁ - h Σ σˣᵢ
//
// Basis: computational (σᶻ) basis, bit i of the state index s is spin i
// (bit = 0 -> σᶻ = +1, bit = 1 -> σᶻ = -1).
//
// H is never stored as a matrix. Its action on a vector is
//   (Hv)[s] = -J·D[s]·v[s] - h·Σᵢ v[s xor 2^i],      D[s] = Σ zᵢzᵢ₊₁
// and the lowest eigenpairs are found with Lanczos. H commutes with the
// global spin flip P = Πσˣ, so each parity sector is solved separately:
// the ground state lives in the even sector, and the lowest odd state gives
// the gap (it becomes degenerate with the ground state in the ordered phase).

import { jacobiEigen } from './linalg.js';

const diagCache = new Map();

/** D[s] = Σ zᵢzᵢ₊₁ for every basis state. Cached per (N, boundary). */
export function bondDiagonal(N, periodic) {
  const key = `${N}:${periodic}`;
  const hit = diagCache.get(key);
  if (hit) return hit;
  const dim = 1 << N;
  const d = new Float64Array(dim);
  const nb = periodic ? N : N - 1;
  for (let s = 0; s < dim; s++) {
    let sum = 0;
    for (let i = 0; i < nb; i++) {
      const j = (i + 1) % N;
      sum += (1 - 2 * ((s >> i) & 1)) * (1 - 2 * ((s >> j) & 1));
    }
    d[s] = sum;
  }
  if (diagCache.size > 6) diagCache.clear();
  diagCache.set(key, d);
  return d;
}

/** out = H v, matrix-free. */
export function applyH(v, out, N, J, h, diag) {
  const dim = 1 << N;
  for (let s = 0; s < dim; s++) {
    let acc = 0;
    for (let b = 1; b < dim; b <<= 1) acc += v[s ^ b];
    out[s] = -J * diag[s] * v[s] - h * acc;
  }
}

/** Dense H, only for tests / tiny N. */
export function buildDenseH(N, J, h, periodic) {
  const dim = 1 << N;
  const H = new Float64Array(dim * dim);
  const diag = bondDiagonal(N, periodic);
  for (let s = 0; s < dim; s++) {
    H[s * dim + s] = -J * diag[s];
    for (let i = 0; i < N; i++) H[s * dim + (s ^ (1 << i))] += -h;
  }
  return H;
}

function dot(a, b) {
  let t = 0;
  for (let i = 0; i < a.length; i++) t += a[i] * b[i];
  return t;
}

function normalize(v) {
  const n = Math.sqrt(dot(v, v));
  for (let i = 0; i < v.length; i++) v[i] /= n;
  return n;
}

/** Project onto the sector with P v = parity·v, where P flips every spin. */
function project(v, parity) {
  const mask = v.length - 1;
  for (let s = 0; s < v.length; s++) {
    const t = mask ^ s;
    if (s < t) {
      const m = 0.5 * (v[s] + parity * v[t]);
      v[s] = m;
      v[t] = parity * m;
    }
  }
}

/** Deterministic pseudo-random start vector inside the chosen sector. */
function startVector(q, parity) {
  let x = parity > 0 ? 0x9e3779b9 : 0x7f4a7c15;
  for (let s = 0; s < q.length; s++) {
    x ^= x << 13; x >>>= 0;
    x ^= x >>> 17;
    x ^= x << 5; x >>>= 0;
    q[s] = x / 4294967296 - 0.5;
  }
  project(q, parity);
  normalize(q);
}

/** Lowest Ritz pair of the k×k tridiagonal matrix (alpha, beta). */
function ritzLowest(alpha, beta, k) {
  const T = new Float64Array(k * k);
  for (let i = 0; i < k; i++) {
    T[i * k + i] = alpha[i];
    if (i + 1 < k) {
      T[i * k + i + 1] = beta[i];
      T[(i + 1) * k + i] = beta[i];
    }
  }
  const { values, vectors } = jacobiEigen(T, k);
  return { value: values[0], vector: vectors[0] };
}

/**
 * Lowest eigenpair inside one parity sector via Lanczos.
 * Two passes: pass 1 builds the tridiagonal matrix (no basis stored), pass 2
 * replays the same recurrence to assemble the eigenvector. Memory stays at
 * four state vectors, so N = 16–18 fits comfortably in a browser.
 */
export function lowestInSector({
  N, J = 1, h, periodic = true, parity = 1,
  wantVector = false, tol = 1e-7, maxIter = 300,
}) {
  const dim = 1 << N;
  const diag = bondDiagonal(N, periodic);
  const maxK = Math.min(maxIter, dim / 2);

  let q = new Float64Array(dim);
  let qprev = new Float64Array(dim);
  const w = new Float64Array(dim);
  startVector(q, parity);

  const alpha = [];
  const beta = [];
  let ritz = null;
  let k = 0;

  for (let j = 0; j < maxK; j++) {
    applyH(q, w, N, J, h, diag);
    const a = dot(q, w);
    alpha.push(a);
    const bp = j > 0 ? beta[j - 1] : 0;
    for (let s = 0; s < dim; s++) w[s] -= a * q[s] + bp * qprev[s];
    project(w, parity);
    const b = Math.sqrt(dot(w, w));

    k = j + 1;
    const last = k === maxK || b < 1e-12;
    if (last || (k >= 6 && k % 3 === 0)) {
      ritz = ritzLowest(alpha, beta, k);
      if (last || b * Math.abs(ritz.vector[k - 1]) < tol) break;
    }
    beta.push(b);
    const tmp = qprev; qprev = q; q = tmp;      // qprev <- q
    for (let s = 0; s < dim; s++) q[s] = w[s] / b;
  }

  const result = { energy: ritz.value, iterations: k, vector: null };
  if (!wantVector) return result;

  // Pass 2: replay and accumulate the Ritz vector.
  const y = ritz.vector;
  const x = new Float64Array(dim);
  q = new Float64Array(dim);
  qprev = new Float64Array(dim);
  startVector(q, parity);
  for (let j = 0; j < k; j++) {
    const c = y[j];
    for (let s = 0; s < dim; s++) x[s] += c * q[s];
    if (j === k - 1) break;
    applyH(q, w, N, J, h, diag);
    const bp = j > 0 ? beta[j - 1] : 0;
    for (let s = 0; s < dim; s++) w[s] -= alpha[j] * q[s] + bp * qprev[s];
    project(w, parity);
    const tmp = qprev; qprev = q; q = tmp;
    for (let s = 0; s < dim; s++) q[s] = w[s] / beta[j];
  }
  normalize(x);
  result.vector = x;
  return result;
}

function popcount(x) {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  return Math.imul((x + (x >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24;
}

/**
 * Ground state, gap and observables for one (N, h).
 *
 * Observables:
 *  - E0, gap = E_odd - E_even
 *  - m2 = ⟨(Σσᶻ)²⟩ / N²   (order parameter, well defined for the symmetric GS)
 *  - sx = ⟨σˣ⟩ averaged over sites
 *  - sites=true also returns per-site ⟨σᶻᵢ⟩ and ⟨σˣᵢ⟩ in the
 *    symmetry-broken state (ψ_even + ψ_odd)/√2, with the overall sign chosen
 *    so the magnetization points up. (In the exact finite-size ground state
 *    ⟨σᶻᵢ⟩ is identically zero by the Z₂ symmetry.)
 */
export function solveChain({ N, J = 1, h, periodic = true, sites = true, tol = 1e-7, keepVectors = false }) {
  const dim = 1 << N;
  const even = lowestInSector({ N, J, h, periodic, parity: +1, wantVector: true, tol });
  const odd = lowestInSector({ N, J, h, periodic, parity: -1, wantVector: sites || keepVectors, tol });
  const psi = even.vector;

  let m2 = 0;
  let sxSum = 0;
  for (let s = 0; s < dim; s++) {
    const p = psi[s];
    const M = N - 2 * popcount(s);
    m2 += p * p * M * M;
    let nb = 0;
    for (let b = 1; b < dim; b <<= 1) nb += psi[s ^ b];
    sxSum += p * nb;
  }

  const res = {
    N, J, h, periodic,
    E0: even.energy,
    E0perSite: even.energy / N,
    Eodd: odd.energy,
    gap: odd.energy - even.energy,
    m2: m2 / (N * N),
    sx: sxSum / N,
    iterations: even.iterations + odd.iterations,
  };

  if (sites) {
    const o = odd.vector;
    const sz = new Array(N).fill(0);
    const sx = new Array(N).fill(0);
    for (let s = 0; s < dim; s++) {
      const es = psi[s], os = o[s];
      for (let i = 0; i < N; i++) {
        const bit = 1 << i;
        sz[i] += es * os * (1 - 2 * ((s >> i) & 1));
        sx[i] += 0.5 * (es * psi[s ^ bit] + os * o[s ^ bit]);
      }
    }
    if (sz.reduce((a, b) => a + b, 0) < 0) for (let i = 0; i < N; i++) sz[i] = -sz[i];
    res.sz = sz;
    res.sxSites = sx;
  }
  if (keepVectors) res.vectors = { even: psi, odd: odd.vector };
  return res;
}
