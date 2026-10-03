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
