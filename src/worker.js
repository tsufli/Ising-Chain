// Web Worker: all heavy numerics live here so the sliders never block.
import { solveChain, sweepPoint } from './physics.js';

let latestSweep = 0;

self.onmessage = async (e) => {
  const msg = e.data;
  try {
    if (msg.type === 'solve') {
      const t0 = performance.now();
      const res = solveChain({ N: msg.N, h: msg.h, periodic: msg.periodic, sites: true, corr: true });
      self.postMessage({ type: 'solve', id: msg.id, res, ms: performance.now() - t0 });
    } else if (msg.type === 'sweep') {
      latestSweep = msg.id;
      for (let k = 0; k < msg.hs.length; k++) {
        if (latestSweep !== msg.id) return; // superseded by a newer sweep
        const res = sweepPoint({ N: msg.N, h: msg.hs[k], periodic: msg.periodic });
        self.postMessage({ type: 'sweep', id: msg.id, i: msg.idx[k], res });
        // yield so queued 'solve' messages (slider moves) get handled between points
        await new Promise((r) => setTimeout(r, 0));
      }
      if (latestSweep === msg.id) self.postMessage({ type: 'sweepDone', id: msg.id });
    }
  } catch (err) {
    self.postMessage({ type: 'error', message: String(err && err.message ? err.message : err) });
  }
};
