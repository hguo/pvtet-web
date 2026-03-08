/**
 * Full single-tet PV classification and stitching.
 * Uses ExactPV2 (PRS + resultants + GCD) for all topological decisions.
 *
 * Given integer V[4][3], W[4][3], produces:
 *   - category string (e.g. "T4_(2,2)_Q3+_SR_Cv2")
 *   - Q/P polynomials, Q roots, discriminant
 *   - punctures with interval assignment, edge/vertex flags
 *   - paired segments via ExactPV2 pure-integer stitching
 *   - degeneracy tags, Cv/Cw positions, bubble detection
 */

import { characteristicPolynomials, polyEval } from './polynomials.js';
import { solveCubic } from './roots.js';
import {
  characteristicPolynomials_bigint, discriminantSign_bigint,
  hasSharedRoot_bigint, allParallel,
} from './bigint_poly.js';
import { solvePVTetV2 } from './exactpv2.js';
import { FACE_VERTS, TET_VERTS, baryTetTo3D } from './curves.js';

// ═══════════════════════════════════════════════════════════════════════
// check_field_zero_in_tet — exact integer Cramer's rule
// ═══════════════════════════════════════════════════════════════════════

/**
 * Check if a 3D vector field F[4][3] vanishes inside the tet interior.
 * Solves sum_i mu_i * F_i = 0, sum mu_i = 1 via exact integer arithmetic.
 * Returns null if no zero inside, or mu[4] (tet barycentric) if found.
 */
function checkFieldZeroInTet(F) {
  // (F[0]-F[3])*mu0 + (F[1]-F[3])*mu1 + (F[2]-F[3])*mu2 = -F[3]
  const A = [
    [F[0][0] - F[3][0], F[1][0] - F[3][0], F[2][0] - F[3][0]],
    [F[0][1] - F[3][1], F[1][1] - F[3][1], F[2][1] - F[3][1]],
    [F[0][2] - F[3][2], F[1][2] - F[3][2], F[2][2] - F[3][2]],
  ];
  const b = [-F[3][0], -F[3][1], -F[3][2]];

  // Exact integer determinant
  const det =
    A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
    A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
    A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  if (det === 0) return null;

  // Cramer numerators
  const n0 =
    b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
    A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) +
    A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2]);
  const n1 =
    A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) -
    b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
    A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0]);
  const n2 =
    A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) -
    A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) +
    b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  const n3 = det - n0 - n1 - n2;

  // Inside test: all mu_k = n_k/det in [0,1]
  if (det > 0) {
    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0) return null;
    if (n0 > det || n1 > det || n2 > det || n3 > det) return null;
  } else {
    if (n0 > 0 || n1 > 0 || n2 > 0 || n3 > 0) return null;
    if (n0 < det || n1 < det || n2 < det || n3 < det) return null;
  }

  return [n0 / det, n1 / det, n2 / det, n3 / det];
}

// ═══════════════════════════════════════════════════════════════════════
// Main classify function
// ═══════════════════════════════════════════════════════════════════════

/**
 * Classify a single tet PV case.
 *
 * @param {number[][]} V - V[4][3] integer field values
 * @param {number[][]} W - W[4][3] integer field values
 * @returns {ClassifiedCase}
 */
export function classifyTetCase(V, W) {
  const cc = {
    V, W,
    category: '',
    Q_coeffs: [0, 0, 0, 0],
    P_coeffs: [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
    Q_bigint: [0n, 0n, 0n, 0n],
    P_bigint: [[0n,0n,0n,0n],[0n,0n,0n,0n],[0n,0n,0n,0n],[0n,0n,0n,0n]],
    Q_disc_sign: 0,
    Q_roots: [],
    n_Q_roots: 0,
    qDegree: 0,
    punctures: [],
    intervals: [],
    pairs: [],
    has_shared_root: false,
    has_B: false,
    has_Cv_pos: false,
    has_Cw_pos: false,
    Cv_mu: null,
    Cw_mu: null,
    SR_lambda: null,
    SR_pos3d: null,
    SR_face: -1,
  };

  // ── Phase 1: Q, P polynomials (float + BigInt) ──
  const { Q, P } = characteristicPolynomials(V, W);
  cc.Q_coeffs = Q;
  cc.P_coeffs = P;

  const { Q_bigint, P_bigint } = characteristicPolynomials_bigint(V, W);
  cc.Q_bigint = Q_bigint;
  cc.P_bigint = P_bigint;

  // ── Phase 2: Q structure ──
  let degQ = 3;
  while (degQ > 0 && Math.abs(Q[degQ]) < 1e-30) degQ--;
  cc.qDegree = degQ;

  const qRoots = solveCubic(Q);
  qRoots.sort((a, b) => a - b);
  cc.Q_roots = qRoots;
  cc.n_Q_roots = qRoots.length;

  // Discriminant sign (exact for degree-3)
  if (degQ === 3) {
    cc.Q_disc_sign = discriminantSign_bigint(Q_bigint);
  } else if (degQ === 2) {
    const qdisc = Q[1] * Q[1] - 4 * Q[0] * Q[2];
    cc.Q_disc_sign = qdisc > 0 ? 1 : qdisc < 0 ? -1 : 0;
  }

  // Q type string
  let qType;
  if (degQ === 0 && Q[0] === 0) qType = 'Qz';
  else if (degQ === 0) qType = 'Q0';
  else if (degQ === 1) qType = 'Q1';
  else if (degQ === 2) {
    if (cc.Q_disc_sign > 0) qType = 'Q2';
    else if (cc.Q_disc_sign < 0) qType = 'Q2-';
    else qType = 'Q2o';
  } else {
    if (cc.Q_disc_sign > 0) qType = 'Q3+';
    else if (cc.Q_disc_sign < 0) qType = 'Q3-';
    else qType = 'Q3o';
  }

  // ── Phase 3: ExactPV2 solver ──
  const pv2 = solvePVTetV2(Q_bigint, P_bigint);

  // Map ExactPV2 punctures to display format
  const punctures = [];
  for (const ep of pv2.punctures) {
    const fv = FACE_VERTS[ep.face];
    let lambda, mu;

    if (ep.isInfinity) {
      lambda = Infinity;
      if (Q[3] !== 0) {
        mu = P.map(pk => pk[3] / Q[3]);
      } else {
        mu = [0.25, 0.25, 0.25, 0.25];
      }
    } else {
      // Solve float P[face] polynomial for roots
      const roots = solveCubic(P[ep.face]).sort((a, b) => a - b);
      lambda = (ep.rootIdx < roots.length) ? roots[ep.rootIdx] : (roots[roots.length - 1] || 0);
      const qVal = polyEval(Q, lambda);
      if (Math.abs(qVal) > 1e-30) {
        mu = P.map(pk => polyEval(pk, lambda) / qVal);
      } else {
        mu = [0.25, 0.25, 0.25, 0.25];
      }
    }

    // Clip and normalize mu
    const muClip = mu.map(m => Math.max(0, m));
    const s = muClip.reduce((a, b) => a + b, 0);
    const muNorm = s > 1e-10 ? muClip.map(m => m / s) : [0.25, 0.25, 0.25, 0.25];

    const pos3d = baryTetTo3D(muNorm);
    const bary = [muNorm[fv[0]], muNorm[fv[1]], muNorm[fv[2]]];

    // Determine tetEdge/tetVertex
    let tetEdge = null, tetVertex = -1;
    if (ep.isVertex) {
      let maxK = 0;
      for (let k = 1; k < 4; k++) if (muNorm[k] > muNorm[maxK]) maxK = k;
      tetVertex = maxK;
    } else if (ep.isEdge && ep.edgeFaces[0] >= 0) {
      tetEdge = [0, 1, 2, 3]
        .filter(v => v !== ep.edgeFaces[0] && v !== ep.edgeFaces[1])
        .sort((a, b) => a - b);
    }

    punctures.push({
      face: ep.face,
      lambda,
      bary,
      mu: muNorm,
      pos3d,
      isEdge: ep.isEdge,
      isVertex: ep.isVertex,
      tetEdge,
      tetVertex,
      faceVertices: [...fv],
      intervalIdx: ep.qInterval,
    });
  }

  cc.punctures = punctures;

  // ── Phase 4: Build intervals from Q roots (for display) ──
  cc.intervals = buildIntervals(qRoots);

  // Update interval occupancy for display
  for (const iv of cc.intervals) iv.n_pv = 0;
  for (const pi of cc.punctures) {
    if (pi.intervalIdx >= 0 && pi.intervalIdx < cc.intervals.length)
      cc.intervals[pi.intervalIdx].n_pv++;
  }

  // ── Phase 5: T-category ──
  const tags = [];

  const nFace = punctures.length;

  // Occupancy tuple from ExactPV2 qInterval
  const occMap = new Map();
  for (const pi of punctures) {
    if (pi.intervalIdx >= 0)
      occMap.set(pi.intervalIdx, (occMap.get(pi.intervalIdx) || 0) + 1);
  }
  const occTuple = [...occMap.values()].sort((a, b) => a - b);

  let tCat = `T${nFace}`;
  if (occTuple.length > 1) tCat += `_(${occTuple.join(',')})`;
  cc.category = tCat + '_' + qType;

  // ── Phase 6: Shared-root detection ──
  cc.has_shared_root = hasSharedRoot_bigint(Q_bigint, P_bigint);
  if (cc.has_shared_root) {
    tags.push('SR');
    findSRLocation(cc);
  }

  // ── Phase 7: Critical-point degeneracies ──
  detectCriticalPoints(cc, tags);
  detectDmd(cc, tags);
  detectBubble(cc, tags);

  if (tags.length > 0) cc.category += '_' + tags.join('_');

  // ── Phase 8: Pairing from ExactPV2 ──
  cc.pairs = [];
  for (const pair of pv2.pairs) {
    const pi_a = pair.a, pi_b = pair.b;
    const a_iv = punctures[pi_a] ? punctures[pi_a].intervalIdx : -1;
    const b_iv = punctures[pi_b] ? punctures[pi_b].intervalIdx : -1;
    const is_cross = a_iv !== b_iv;
    cc.pairs.push({ pi_a, pi_b, is_cross, interval_idx: Math.min(a_iv, b_iv) });
  }

  return cc;
}

// ═══════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════

function buildIntervals(qRoots) {
  const roots = [...qRoots].sort((a, b) => a - b);
  const intervals = [];
  if (roots.length === 0) {
    intervals.push({ lb: null, ub: null, isInfinity: true, n_pv: 0 });
  } else {
    intervals.push({ lb: null, ub: roots[0], isInfinity: true, n_pv: 0 });
    for (let i = 0; i < roots.length - 1; i++)
      intervals.push({ lb: roots[i], ub: roots[i + 1], isInfinity: false, n_pv: 0 });
    intervals.push({ lb: roots[roots.length - 1], ub: null, isInfinity: true, n_pv: 0 });
  }
  return intervals;
}

// ═══════════════════════════════════════════════════════════════════════
// Degeneracy detection
// ═══════════════════════════════════════════════════════════════════════

/** Detect Cv/Cw critical points (interior + boundary). */
function detectCriticalPoints(cc, tags) {
  const { V, W, Q_bigint, P_bigint } = cc;

  // Interior: check if v=0 or w=0 inside tet
  const cvMu = checkFieldZeroInTet(V);
  const cwMu = checkFieldZeroInTet(W);

  let hasCv = cvMu !== null;
  let hasCw = cwMu !== null;
  let hasC2v = false, hasC1v = false, hasC0v = false;
  let hasC2w = false, hasC1w = false, hasC0w = false;

  if (hasCv) { cc.has_Cv_pos = true; cc.Cv_mu = cvMu; }
  if (hasCw) { cc.has_Cw_pos = true; cc.Cw_mu = cwMu; }

  // Boundary Cv: λ=0 on face/edge/vertex (check BigInt P[k][0])
  if (Q_bigint[0] !== 0n) {
    for (let k = 0; k < 4; k++) {
      if (P_bigint[k][0] !== 0n) continue;
      let inside = true;
      let nzero = 0;
      for (let j = 0; j < 4; j++) {
        if (j === k) continue;
        if (P_bigint[j][0] === 0n) { nzero++; continue; }
        if (P_bigint[j][0] * Q_bigint[0] < 0n) { inside = false; break; }
      }
      if (!inside) continue;
      if (nzero >= 2) hasC0v = true;
      else if (nzero >= 1) hasC1v = true;
      else hasC2v = true;
      break;
    }
  }

  // Boundary Cw: λ=∞ on face/edge/vertex (check BigInt P[k][3])
  if (Q_bigint[3] !== 0n) {
    for (let k = 0; k < 4; k++) {
      if (P_bigint[k][3] !== 0n) continue;
      let inside = true;
      let nzero = 0;
      for (let j = 0; j < 4; j++) {
        if (j === k) continue;
        if (P_bigint[j][3] === 0n) { nzero++; continue; }
        if (P_bigint[j][3] * Q_bigint[3] < 0n) { inside = false; break; }
      }
      if (!inside) continue;
      if (nzero >= 2) hasC0w = true;
      else if (nzero >= 1) hasC1w = true;
      else hasC2w = true;
      break;
    }
  }

  // Emit single most-specific Cv{d} tag
  if (hasCv || hasC2v || hasC1v || hasC0v) {
    if (hasC0v) tags.push('Cv0');
    else if (hasC1v) tags.push('Cv1');
    else if (hasC2v) tags.push('Cv2');
    else tags.push('Cv');
  }

  // Emit single most-specific Cw{d} tag
  if (hasCw || hasC2w || hasC1w || hasC0w) {
    if (hasC0w) tags.push('Cw0');
    else if (hasC1w) tags.push('Cw1');
    else if (hasC2w) tags.push('Cw2');
    else tags.push('Cw');
  }
}

/** Detect D00/D01 vertex/edge degeneracies. */
function detectDmd(cc, tags) {
  const { V, W, punctures } = cc;

  // D00: vertex where V×W = 0
  let hasD00 = false;
  for (let i = 0; i < 4; i++) {
    const cx = V[i][1] * W[i][2] - V[i][2] * W[i][1];
    const cy = V[i][2] * W[i][0] - V[i][0] * W[i][2];
    const cz = V[i][0] * W[i][1] - V[i][1] * W[i][0];
    if (cx === 0 && cy === 0 && cz === 0) hasD00 = true;
  }

  // D01: edge puncture at generic λ (ExactPV2 already excludes waypoints)
  let hasD01 = false;
  for (const pi of punctures) {
    if (pi.isEdge) hasD01 = true;
  }

  if (hasD00) tags.push('D00');
  if (hasD01) tags.push('D01');
}

/** Detect bubble (closed PV curve with no punctures). */
function detectBubble(cc, tags) {
  const { Q_coeffs: Q, P_coeffs: P, n_Q_roots } = cc;

  // ExactPV2 already excludes waypoints, so all punctures are valid
  if (cc.punctures.length > 0) return;
  if (n_Q_roots > 0) return;
  if (Q[0] === 0) return;

  // Check all μ_k(0) = P_k[0]/Q[0] ≥ 0
  const inside = P.every(pk => pk[0] / Q[0] >= -1e-10);
  if (inside) {
    cc.has_B = true;
    tags.push('B');
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Shared root location
// ═══════════════════════════════════════════════════════════════════════

/**
 * Find the SR (shared root) location: which Q-root is shared with some P_k.
 * Sets cc.SR_lambda, cc.SR_pos3d, cc.SR_face.
 */
function findSRLocation(cc) {
  const { Q_coeffs: Q, P_coeffs: P, Q_roots: roots, punctures } = cc;

  function polyDerivEval(coeffs, x) {
    return coeffs[1] + 2 * coeffs[2] * x + 3 * coeffs[3] * x * x;
  }

  for (const root of roots) {
    for (let k = 0; k < 4; k++) {
      const pk_val = Math.abs(polyEval(P[k], root));
      const scale = Math.max(1, Math.abs(P[k][3]) + Math.abs(P[k][2]) + Math.abs(P[k][1]) + Math.abs(P[k][0]));
      if (pk_val < 1e-4 * scale) {
        cc.SR_lambda = root;
        cc.SR_face = k;

        // Position: find nearest puncture to this root
        let bestDist = Infinity, bestPi = -1;
        for (let i = 0; i < punctures.length; i++) {
          if (!isFinite(punctures[i].lambda)) continue;
          const d = Math.abs(punctures[i].lambda - root);
          if (d < bestDist) { bestDist = d; bestPi = i; }
        }
        if (bestPi >= 0 && bestDist < 1e-3) {
          cc.SR_pos3d = punctures[bestPi].pos3d;
          return;
        }

        // L'Hôpital: μ_j(root) = P_j'(root)/Q'(root)
        const qDeriv = polyDerivEval(Q, root);
        if (Math.abs(qDeriv) > 1e-10) {
          const mu = P.map(pk => polyDerivEval(pk, root) / qDeriv);
          if (mu.every(m => m >= -0.05) && mu.every(m => m <= 1.05)) {
            const muClip = mu.map(m => Math.max(0, Math.min(1, m)));
            const s = muClip.reduce((a, b) => a + b, 0);
            if (s > 1e-10) {
              cc.SR_pos3d = baryTetTo3D(muClip.map(m => m / s));
              return;
            }
          }
        }

        // Fallback: try nearby evaluation on both sides
        for (const eps of [1e-4, -1e-4, 1e-2, -1e-2]) {
          const lam = root + eps;
          const qv = polyEval(Q, lam);
          if (Math.abs(qv) < 1e-20) continue;
          const mu = P.map(pk => polyEval(pk, lam) / qv);
          if (mu.every(m => m >= -0.05) && mu.every(m => m <= 1.05)) {
            const muClip = mu.map(m => Math.max(0, Math.min(1, m)));
            const s = muClip.reduce((a, b) => a + b, 0);
            if (s > 1e-10) {
              cc.SR_pos3d = baryTetTo3D(muClip.map(m => m / s));
              return;
            }
          }
        }

        cc.SR_pos3d = null;
        return;
      }
    }
  }
}
