# Transverse-field Ising chain: exact diagonalization in the browser

Interactive simulation of the 1D quantum Ising chain

```
H = −J Σᵢ σᶻᵢ σᶻᵢ₊₁ − h Σᵢ σˣᵢ
```

Drag the transverse field **h/J** and watch the ground state go from ferromagnetic (h < J) to
paramagnetic (h > J). Everything is computed live in the browser: exact ground state, gap, order
parameter and per-site spin expectation values, for chains of up to 16 spins (65,536 basis states).

**Live demo:** `https://<your-username>.github.io/<repo-name>/` (see [Deploy](#deploy))

## What you see

- **Chain of arrows:** the local spin expectation (⟨σˣᵢ⟩, ⟨σᶻᵢ⟩). Up = ordered, tilted toward x = driven by the field.
- **Low-lying spectrum:** the lowest levels of the even (solid) and odd (dashed) parity sectors, measured from
  the ground state, versus h/J. The lowest odd level merges with the ground state in the ordered phase.
- **Order parameter** m² = ⟨(Σσᶻ)²⟩/N² versus h/J.
- **Two-spin correlation** C(r) = ⟨σᶻᵢσᶻᵢ₊ᵣ⟩ versus distance r, recomputed live for the current h, N and boundary.
- **Sliders:** h/J moves the transition; N (4–16 spins) shows finite-size rounding; a toggle switches between a
  periodic ring and an open chain. The curves have a dashed line at h = J and a marker at the slider position.

## Method

- **Matrix-free Lanczos.** H is never stored. Its action on a vector uses bit operations:
  `(Hv)[s] = −J·D[s]·v[s] − h·Σᵢ v[s ⊕ 2ⁱ]`, where `D[s] = Σ zᵢzᵢ₊₁` is precomputed once per (N, boundary).
- **Parity sectors.** H commutes with the global spin flip Πσˣ. The even and odd sectors are solved
  separately: the ground state is the lowest even state, and the gap is `E_odd − E_even`.
- **Two-pass Lanczos.** Pass 1 builds the tridiagonal matrix and stops when the Ritz residual is below 1e‑7;
  pass 2 replays the recurrence to assemble the eigenvector. Only four state vectors are in memory at once.
- **Web Worker.** All numerics run off the main thread. The sliders stay responsive, and the curves are
  filled in progressively with a cancellable sweep.

### Spectrum

Each parity sector is run through Lanczos until the lowest three Ritz pairs have residual < 1e‑6. Plain Lanczos
without re-orthogonalization produces "ghost" copies of converged eigenvalues, so Ritz values closer than the sum
of their residual bounds are merged into one level. A single start vector sees each distinct eigenvalue once, so
degenerate levels (±k momentum pairs on a ring) appear once. The sweep over h runs coarse-first so a rough curve
shows up immediately. The spectrum curves start at h = 0.05, because h = 0 is the degenerate classical limit.

### Correlations

C(r) is evaluated directly from the ground-state amplitudes, `Σₛ |ψ(s)|² (N − 2·popcount(s ⊕ rot_r(s)))/N`
(periodic) or the analogous open-chain average over the N − r pairs that fit. In the symmetric ground state
⟨σᶻ⟩ = 0, so this is also the connected correlator.

### Observables

- `m² = ⟨(Σᵢ σᶻᵢ)²⟩ / N²`: order parameter, well defined for the exact (symmetric) ground state.
- `⟨σˣ⟩ = (1/N) Σᵢ ⟨σˣᵢ⟩`.
- **Why the arrows use a superposition.** On a finite chain the exact ground state is symmetric under the
  spin flip, so ⟨σᶻᵢ⟩ = 0 identically. The arrows show `(ψ_even + ψ_odd)/√2`, which is the symmetry-broken
  ground state when the two lowest states are nearly degenerate (h < J). In the disordered phase this
  state is not physically meaningful as a "magnetized" state and a small residual ⟨σᶻ⟩ remains. It is a
  finite-size effect that shrinks as N grows (compare N = 8 and N = 16 at h = 1.7).

## Validation

`npm test` (plain Node, no dependencies) checks:

- the analytic limits h = 0 (`E₀ = −JN` periodic, `−J(N−1)` open) and J = 0 (`E₀ = −hN`, gap `2h`)
- Lanczos vs. a dense Jacobi diagonalization for N = 4–7, both boundary conditions: ground-state energy,
  lowest odd-sector energy, m², ⟨σˣ⟩ and per-site ⟨σᶻ⟩
- eigen-residuals ‖Hψ − Eψ‖ at N = 12
- the tridiagonal QL eigen-solver against Jacobi
- the lowest three distinct levels of each parity sector against dense diagonalization (N = 6, 7; both boundaries)
- C(r) against dense diagonalization, plus the limits h = 0 (C = 1 for all r) and J = 0 (C = 0 for r ≥ 1)
- translation invariance of the periodic chain
- qualitative phase-transition signatures (m² and gap trends)

## Run locally

ES modules and workers need an HTTP server (opening `index.html` as a file won't work):

```bash
npm start              # python3 -m http.server 8000
# or: npx serve .
```

Then open http://localhost:8000.

## Deploy

The repo has no build step. The included workflow runs the tests and publishes the site to GitHub Pages
on every push to `main`.

1. Create a new repository on GitHub and push this folder to `main`.
2. In the repository: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
3. The site appears at `https://<your-username>.github.io/<repo-name>/` after the workflow finishes.

## Project layout

```
index.html, style.css      page
src/main.js                UI: sliders, arrow chain, charts
src/worker.js              Web Worker wrapper
src/physics.js             H·v, parity-sector Lanczos, spectrum, observables, C(r)
src/linalg.js              Jacobi and tridiagonal eigen-solvers
test/run.js                tests
.github/workflows/pages.yml  test + deploy
```

## Possible extensions

- Correlation length ξ(h) extracted from C(r), and a log-scale view of C(r)
- Larger N (18+) with translation-symmetry sectors
- Quench dynamics: time evolution after a sudden change of h
