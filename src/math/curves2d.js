/**
 * 2D Triangle geometry, curve sampling, and coordinate conversions.
 */
// Equilateral triangle in 2D
export const TRI_VERTS_2D = [
  [0, 0],
  [1, 0],
  [0.5, Math.sqrt(3) / 2],
];

export const TRI_EDGES = [[0, 1], [1, 2], [2, 0]];

// Face k is opposite vertex k
// FACE_EDGES_2D[k] = [va, vb] — the two vertices of the edge opposite vertex k
export const FACE_EDGES_2D = [
  [1, 2],  // face 0: edge between v₁ and v₂ (opposite v₀)
  [2, 0],  // face 1: edge between v₂ and v₀ (opposite v₁)
  [0, 1],  // face 2: edge between v₀ and v₁ (opposite v₂)
];

export const SEGMENT_COLORS_2D = [
  '#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3', '#a65628',
];

/** Convert triangle barycentric coords (3 values) to 2D position. */
export function baryTriTo2D(mu) {
  const p = [0, 0];
  for (let i = 0; i < 3; i++)
    for (let d = 0; d < 2; d++)
      p[d] += mu[i] * TRI_VERTS_2D[i][d];
  return p;
}

/** λ to triangle barycentric coords: μ_i = P_i(λ) / Q(λ). */
export function lambdaToBary2D(lam, Q, P) {
  const qVal = Q[0] + Q[1] * lam + Q[2] * lam * lam;
  if (Math.abs(qVal) < 1e-30) return null;
  return [
    (P[0][0] + P[0][1] * lam + P[0][2] * lam * lam) / qVal,
    (P[1][0] + P[1][1] * lam + P[1][2] * lam * lam) / qVal,
    (P[2][0] + P[2][1] * lam + P[2][2] * lam * lam) / qVal,
  ];
}

/** Check if 2D barycentric coords are inside triangle. */
function isInside2D(mu) {
  if (!mu) return false;
  for (let i = 0; i < 3; i++)
    if (mu[i] < -1e-6 || mu[i] > 1 + 1e-6) return false;
  return true;
}

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
  for (let i = 1; i < sorted.length; i++)
    if (sorted[i] - sorted[i - 1] > 1e-15) out.push(sorted[i]);
  return out;
}

/**
 * Sample PV curve in 2D for visualization.
 * Returns array of { pts: number[][], lams: number[], lamEntry, lamExit } sub-segments.
 */
export function samplePVCurve2D(Q, P, lamLo, lamHi, isInfinity, nSamples = 200, lamHints = []) {
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

  if (lamHints.length > 0) {
    let extra = [];
    for (const h of lamHints) {
      if (!isFinite(h)) continue;
      const r = Math.max(Math.abs(h) * 0.1, 0.1);
      extra.push(...linspace(h - r, h + r, 40));
      const r2 = Math.max(Math.abs(h) * 0.01, 0.01);
      extra.push(...linspace(h - r2, h + r2, 20));
    }
    const finiteHints = lamHints.filter(h => isFinite(h)).sort((a, b) => a - b);
    for (let i = 0; i + 1 < finiteHints.length; i++)
      extra.push(...linspace(finiteHints[i], finiteHints[i + 1], 60));
    const loInf2 = lamLo === null;
    const hiInf2 = lamHi === null;
    if (!loInf2 && !hiInf2) extra = extra.filter(v => v >= lamLo && v <= lamHi);
    else if (!loInf2) extra = extra.filter(v => v >= lamLo);
    else if (!hiInf2) extra = extra.filter(v => v <= lamHi);
    lamVals = sortedUnique([...lamVals, ...extra]);
  }

  function bisectBoundary(lamIn, lamOut) {
    for (let i = 0; i < 30; i++) {
      const mid = 0.5 * (lamIn + lamOut);
      const mu = lambdaToBary2D(mid, Q, P);
      if (isInside2D(mu)) lamIn = mid; else lamOut = mid;
    }
    const mu = lambdaToBary2D(lamIn, Q, P);
    return { pt: baryTriTo2D(clipMu(mu)), lam: lamIn };
  }

  const segments = [];
  let current = [], currentLams = [];
  let segLamEntry = null, segLamExit = null;
  let prevLam = null, prevInside = false;

  for (const lam of lamVals) {
    const mu = lambdaToBary2D(lam, Q, P);
    const inside = isInside2D(mu);

    if (inside) {
      if (!prevInside && prevLam !== null) {
        const { pt, lam: bl } = bisectBoundary(lam, prevLam);
        current.push(pt); currentLams.push(bl);
        segLamEntry = bl;
      } else if (!prevInside) {
        segLamEntry = lam;
      }
      current.push(baryTriTo2D(mu));
      currentLams.push(lam);
      segLamExit = lam;
    } else {
      if (prevInside && prevLam !== null) {
        const { pt, lam: bl } = bisectBoundary(prevLam, lam);
        current.push(pt); currentLams.push(bl);
        segLamExit = bl;
      }
      if (current.length > 1) {
        segments.push({ pts: [...current], lams: [...currentLams],
                        lamEntry: segLamEntry, lamExit: segLamExit });
      }
      current = []; currentLams = [];
    }
    prevLam = lam;
    prevInside = inside;
  }
  if (current.length > 1) {
    segments.push({ pts: [...current], lams: [...currentLams],
                    lamEntry: segLamEntry, lamExit: segLamExit });
  }
  return segments;
}

/**
 * Sample bubble (closed PV curve) in 2D triangle.
 * Returns { pts, lams } or null.
 */
export function sampleBubble2D(Q, P) {
  const nSamples = 400;
  const t = linspace(-0.499 * Math.PI, 0.499 * Math.PI, nSamples);
  const pts = [], lams = [];

  for (const ti of t) {
    const lam = Math.tan(ti);
    const mu = lambdaToBary2D(lam, Q, P);
    if (mu && mu.every(m => m >= -0.01)) {
      const muClip = clipMu(mu);
      const s = muClip.reduce((a, b) => a + b, 0);
      if (s > 1e-10) {
        const muNorm = muClip.map(v => v / s);
        pts.push(baryTriTo2D(muNorm));
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
