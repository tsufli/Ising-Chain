// Small dense symmetric eigensolver (cyclic Jacobi).
// Used for (a) the tiny tridiagonal Lanczos matrices and (b) the dense
// reference diagonalization in the tests.

/**
 * @param {Float64Array} A  symmetric n×n matrix, row-major. Destroyed.
 * @param {number} n
 * @returns {{values: Float64Array, vectors: Float64Array[]}} eigenvalues in
 *          ascending order; vectors[k] is the eigenvector of values[k].
 */
export function jacobiEigen(A, n, maxSweeps = 60) {
  const V = new Float64Array(n * n);
  for (let i = 0; i < n; i++) V[i * n + i] = 1;

  let total = 0;
  for (let i = 0; i < A.length; i++) total += A[i] * A[i];
  const stop = 1e-26 * Math.max(total, 1e-300);

  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++)
      for (let q = p + 1; q < n; q++) off += A[p * n + q] * A[p * n + q];
    if (off <= stop) break;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p * n + q];
        if (Math.abs(apq) < 1e-300) continue;
        const theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
        const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = A[k * n + p], akq = A[k * n + q];
          A[k * n + p] = c * akp - s * akq;
          A[k * n + q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = A[p * n + k], aqk = A[q * n + k];
          A[p * n + k] = c * apk - s * aqk;
          A[q * n + k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }

  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => A[a * n + a] - A[b * n + b]);
  const values = new Float64Array(n);
  const vectors = [];
  order.forEach((col, k) => {
    values[k] = A[col * n + col];
    const v = new Float64Array(n);
    for (let r = 0; r < n; r++) v[r] = V[r * n + col];
    vectors.push(v);
  });
  return { values, vectors };
}

/**
 * Eigenvalues (ascending) of the symmetric tridiagonal matrix with diagonal
 * alpha[0..k-1] and off-diagonal beta[0..k-2]. Implicit QL, O(k²), no vectors.
 */
export function tridiagEigenvalues(alpha, beta, k) {
  const d = new Float64Array(k);
  const e = new Float64Array(k);
  for (let i = 0; i < k; i++) d[i] = alpha[i];
  for (let i = 0; i < k - 1; i++) e[i] = beta[i];

  for (let l = 0; l < k; l++) {
    let iter = 0;
    let m;
    do {
      for (m = l; m < k - 1; m++) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m]) <= 1e-16 * dd) break;
      }
      if (m !== l) {
        if (iter++ === 60) break;
        let g = (d[l + 1] - d[l]) / (2 * e[l]);
        let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? Math.abs(r) : -Math.abs(r)));
        let s = 1, c = 1, p = 0;
        let i;
        for (i = m - 1; i >= l; i--) {
          const f = s * e[i];
          const b = c * e[i];
          r = Math.hypot(f, g);
          e[i + 1] = r;
          if (r === 0) { d[i + 1] -= p; e[m] = 0; break; }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
        }
        if (r === 0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0;
      }
    } while (m !== l);
  }
  d.sort();
  return d;
}

/**
 * Unit eigenvector of the same tridiagonal matrix for a known eigenvalue
 * `lambda`, by inverse iteration (O(k) per step).
 */
export function tridiagEigenvector(alpha, beta, k, lambda) {
  const mu = lambda - 1e-10 * Math.max(1, Math.abs(lambda));
  const dd = new Float64Array(k);
  let y = new Float64Array(k).fill(1 / Math.sqrt(k));
  for (let it = 0; it < 3; it++) {
    const r = Float64Array.from(y);
    dd[0] = alpha[0] - mu || 1e-30;
    for (let i = 1; i < k; i++) {
      const m = beta[i - 1] / dd[i - 1];
      dd[i] = alpha[i] - mu - m * beta[i - 1] || 1e-30;
      r[i] -= m * r[i - 1];
    }
    const x = new Float64Array(k);
    x[k - 1] = r[k - 1] / dd[k - 1];
    for (let i = k - 2; i >= 0; i--) x[i] = (r[i] - beta[i] * x[i + 1]) / dd[i];
    let n = 0;
    for (let i = 0; i < k; i++) n += x[i] * x[i];
    n = Math.sqrt(n);
    for (let i = 0; i < k; i++) x[i] /= n;
    y = x;
  }
  return y;
}
