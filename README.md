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
- **Four curves** with a marker at the slider position and a dashed line at the critical point h = J:
  order parameter m², energy gap, ground-state energy per site, transverse magnetization ⟨σˣ⟩.
- **N slider** to see finite-size rounding: the transition sharpens as N grows.
- **Boundary toggle:** periodic ring or open chain (edge effects show up in the arrows).

## Method

- **Matrix-free Lanczos.** H is never stored. Its action on a vector uses bit operations:
  `(Hv)[s] = −J·D[s]·v[s] − h·Σᵢ v[s ⊕ 2ⁱ]`, where `D[s] = Σ zᵢzᵢ₊₁` is precomputed once per (N, boundary).
- **Parity sectors.** H commutes with the global spin flip Πσˣ. The even and odd sectors are solved
  separately: the ground state is the lowest even state, and the gap is `E_odd − E_even`.
- **Two-pass Lanczos.** Pass 1 builds the tridiagonal matrix and stops when the Ritz residual is below 1e‑7;
  pass 2 replays the recurrence to assemble the eigenvector. Only four state vectors are in memory at once.
- **Web Worker.** All numerics run off the main thread. The sliders stay responsive, and the curves are
  filled in progressively with a cancellable sweep.

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
src/physics.js             H·v, parity-sector Lanczos, observables
src/linalg.js              small Jacobi eigensolver
test/run.js                tests
.github/workflows/pages.yml  test + deploy
```

## Possible extensions

- Correlation function ⟨σᶻᵢσᶻᵢ₊ᵣ⟩ and correlation length near h = J
- Larger N (18+) with translation-symmetry sectors
- Quench dynamics: time evolution after a sudden change of h
