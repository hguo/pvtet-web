/**
 * PV curve sampling: sample the PV curve λ → μ(λ) → 3D inside a tet.
 * Ported from pv_tet_visualize.py.
 *
 * This module provides ONLY the visualization sampling functions.
 * Combinatorial detection (punctures, pairing) is done by triangle_solver.js
 * and the stitching logic in state.js.
 */
import { polyEval } from './polynomials.js';

// Regular tetrahedron vertices
export const TET_VERTS = [
  [0, 0, 0],
  [1, 0, 0],
  [0.5, Math.sqrt(3) / 2, 0],
  [0.5, Math.sqrt(3) / 6, Math.sqrt(6) / 3],
];

export const TET_EDGES = [[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]];

export const FACE_VERTS = [
  [1, 3, 2],  // face 0 (opposite vertex 0)
  [0, 2, 3],  // face 1
  [0, 3, 1],  // face 2
  [0, 1, 2],  // face 3
];

export const SEGMENT_COLORS = [
  '#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3', '#a65628',
];

/** Convert tet barycentric coords (4 values) to 3D. */
export function baryTetTo3D(mu) {
  const p = [0, 0, 0];
  for (let i = 0; i < 4; i++)
    for (let d = 0; d < 3; d++)
      p[d] += mu[i] * TET_VERTS[i][d];
  return p;
}

/** λ to tet barycentric coords: μ_i = P_i(λ) / Q(λ). Returns null if Q≈0. */
export function lambdaToBary(lam, Q, P) {
  const qVal = polyEval(Q, lam);
  if (Math.abs(qVal) < 1e-30) return null;
  return [
    polyEval(P[0], lam) / qVal,
    polyEval(P[1], lam) / qVal,
    polyEval(P[2], lam) / qVal,
    polyEval(P[3], lam) / qVal,
  ];
}

/** Check if barycentric coords are inside tet (with tolerance). */
function isInside(mu) {
  if (!mu) return false;
  for (let i = 0; i < 4; i++) {
    if (mu[i] < -1e-6 || mu[i] > 1 + 1e-6) return false;
  }
  return true;
}

/**
 * Sample PV curve in a λ-interval for visualization.
 * Returns array of { pts: number[][], lamEntry, lamExit } sub-segments.
 * @param lamHints - optional array of known puncture λ values for dense sampling
 */
export function samplePVCurve(Q, P, lamLo, lamHi, isInfinity, nSamples = 200, lamHints = []) {
  let lamVals;

  if (isInfinity) {
    const loInf = lamLo === null;
    const hiInf = lamHi === null;
    if (loInf && hiInf) {
      const t = linspace(-0.499 * Math.PI, 0.499 * Math.PI, nSamples);
      lamVals = t.map(v => Math.tan(v));
    } else if (loInf) {
      const t = linspace(-0.499 * Math.PI, 0, nSamples);
      const s = Math.max(Math.abs(lamHi), 1) * 10;
      lamVals = t.map(v => lamHi + s * Math.tan(v));
    } else {
      const t = linspace(0, 0.499 * Math.PI, nSamples);
      const s = Math.max(Math.abs(lamLo), 1) * 10;
      lamVals = t.map(v => lamLo + s * Math.tan(v));
    }
  } else {
    const mid = linspace(lamLo, lamHi, nSamples);
    const eps = (lamHi - lamLo) * 0.01;
    const nearLo = linspace(lamLo, lamLo + eps, 10);
    const nearHi = linspace(lamHi - eps, lamHi, 10);
    lamVals = sortedUnique([...mid, ...nearLo, ...nearHi]);
  }

  // Add dense samples around hint λ values (known puncture locations)
  if (lamHints.length > 0) {
    let extra = [];
    for (const h of lamHints) {
      if (!isFinite(h)) continue;
      // Dense samples in neighborhood of each hint
      const r = Math.max(Math.abs(h) * 0.1, 0.1);
      extra.push(...linspace(h - r, h + r, 40));
      // Very fine near the hint
      const r2 = Math.max(Math.abs(h) * 0.01, 0.01);
      extra.push(...linspace(h - r2, h + r2, 20));
    }
    // Also densely sample between hint pairs
    const finiteHints = lamHints.filter(h => isFinite(h)).sort((a, b) => a - b);
    for (let i = 0; i + 1 < finiteHints.length; i++) {
      extra.push(...linspace(finiteHints[i], finiteHints[i + 1], 60));
    }
    // Clip hints to the valid sampling range to avoid spurious sub-segments
    const loInf2 = lamLo === null;
    const hiInf2 = lamHi === null;
    if (!loInf2 && !hiInf2) {
      extra = extra.filter(v => v >= lamLo && v <= lamHi);
    } else if (!loInf2) {
      extra = extra.filter(v => v >= lamLo);
    } else if (!hiInf2) {
      extra = extra.filter(v => v <= lamHi);
    }
    lamVals = sortedUnique([...lamVals, ...extra]);
  }

  function bisectBoundary(lamIn, lamOut) {
    for (let i = 0; i < 30; i++) {
      const lamMid = 0.5 * (lamIn + lamOut);
      const muMid = lambdaToBary(lamMid, Q, P);
      if (isInside(muMid)) lamIn = lamMid;
      else lamOut = lamMid;
    }
    const mu = lambdaToBary(lamIn, Q, P);
    return { pt: baryTetTo3D(clipMu(mu)), lam: lamIn };
  }

  const segments = [];
  let current = [];
  let currentLams = [];
  let segLamEntry = null;
  let segLamExit = null;
  let prevLam = null;
  let prevInside = false;

  for (const lam of lamVals) {
    const mu = lambdaToBary(lam, Q, P);
    const inside = isInside(mu);

    if (inside) {
      if (!prevInside && prevLam !== null) {
        const { pt, lam: bl } = bisectBoundary(lam, prevLam);
        current.push(pt);
        currentLams.push(bl);
        segLamEntry = bl;
      } else if (!prevInside) {
        segLamEntry = lam;
      }
      current.push(baryTetTo3D(mu));
      currentLams.push(lam);
      segLamExit = lam;
    } else {
      if (prevInside && prevLam !== null) {
        const { pt, lam: bl } = bisectBoundary(prevLam, lam);
        current.push(pt);
        currentLams.push(bl);
        segLamExit = bl;
      }
      if (current.length > 1) {
        segments.push({ pts: [...current], lams: [...currentLams], lamEntry: segLamEntry, lamExit: segLamExit });
      }
      current = [];
      currentLams = [];
    }
    prevLam = lam;
    prevInside = inside;
  }

  if (current.length > 1) {
    segments.push({ pts: [...current], lams: [...currentLams], lamEntry: segLamEntry, lamExit: segLamExit });
  }

  return segments;
}

/** Clip barycentric coords to [0,1]. */
function clipMu(mu) {
  return mu.map(v => Math.max(0, Math.min(1, v)));
}

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const arr = [];
  for (let i = 0; i < n; i++) arr.push(a + (b - a) * i / (n - 1));
  return arr;
}

function sortedUnique(arr) {
  const sorted = arr.slice().sort((a, b) => a - b);
  const out = [sorted[0]];
  for (let i = 1; i < sorted.length; i++) {
    if (sorted[i] - sorted[i - 1] > 1e-15) out.push(sorted[i]);
  }
  return out;
}

/**
 * Check for bubble: closed PV curve inside tet with 0 punctures.
 * Returns { pts, lams } or null.
 */
export function sampleBubble(Q, P) {
  const nSamples = 400;
  const t = linspace(-0.499 * Math.PI, 0.499 * Math.PI, nSamples);
  const pts = [], lams = [];

  for (const ti of t) {
    const lam = Math.tan(ti);
    const mu = lambdaToBary(lam, Q, P);
    if (mu && mu.every(m => m >= -0.01)) {
      const muClip = clipMu(mu);
      const s = muClip.reduce((a, b) => a + b, 0);
      if (s > 1e-10) {
        const muNorm = muClip.map(v => v / s);
        pts.push(baryTetTo3D(muNorm));
        lams.push(lam);
      }
    }
  }

  if (pts.length > 2) {
    pts.push(pts[0]);
    lams.push(lams[0]);
    return { pts, lams };
  }
  return null;
}
