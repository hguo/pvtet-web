/**
 * Full single-tet PV classification and stitching.
 * Mirrors ftk2 pv_tet_case_finder.cu classify_case().
 *
 * Given integer V[4][3], W[4][3], produces:
 *   - category string (e.g. "T4_(2,2)_Q3+_SR_Cv2")
 *   - Q/P polynomials, Q roots, discriminant
 *   - punctures with interval assignment, edge/vertex flags
 *   - paired segments via Sturm-based stitching
 *   - degeneracy tags, Cv/Cw positions, bubble detection
 */

import { characteristicPolynomials, polyEval } from './polynomials.js';
import { solveCubic } from './roots.js';
import {
  characteristicPolynomials_bigint, discriminantSign_bigint,
  hasSharedRoot_bigint, gcdNormalizePoly, allParallel,
} from './bigint_poly.js';
import {
  buildSturmCubic, sturmCountCubic, sturmCountCubicCertified,
} from './sturm.js';
import { solvePVTriangle } from './triangle_solver.js';
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
    // Quadratic discriminant: Q[1]²-4Q[0]Q[2]
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

  // ── Phase 3: Puncture extraction (triangle solver on each face) ──
  const rawPunctures = [];
  for (let f = 0; f < 4; f++) {
    const fv = FACE_VERTS[f];
    const V_face = [V[fv[0]], V[fv[1]], V[fv[2]]];
    const W_face = [W[fv[0]], W[fv[1]], W[fv[2]]];

    const result = solvePVTriangle(V_face, W_face, fv);
    if (result === null) continue; // entire face is PV (degenerate)

    for (const p of result) {
      // Face bary → 3D position
      const pos3d = [0, 0, 0];
      for (let i = 0; i < 3; i++)
        for (let d = 0; d < 3; d++)
          pos3d[d] += p.bary[i] * TET_VERTS[fv[i]][d];

      // Face bary → tet bary
      const mu = [0, 0, 0, 0];
      for (let i = 0; i < 3; i++) mu[fv[i]] += p.bary[i];

      // Determine tet edge/vertex for edge/vertex punctures
      let tetEdge = null, tetVertex = -1;
      if (p.isVertex) {
        // Find which tet vertex (bary ≈ 1)
        let maxK = 0;
        for (let k = 1; k < 3; k++) if (p.bary[k] > p.bary[maxK]) maxK = k;
        tetVertex = fv[maxK];
      } else if (p.isEdge) {
        // Find which tet edge (two nonzero bary coords)
        const nonzero = [];
        for (let k = 0; k < 3; k++) if (Math.abs(p.bary[k]) > 1e-10) nonzero.push(fv[k]);
        if (nonzero.length === 2) tetEdge = nonzero.sort((a, b) => a - b);
      }

      rawPunctures.push({
        face: f,
        lambda: p.lambda,
        bary: p.bary,
        mu,
        pos3d,
        isEdge: p.isEdge,
        isVertex: p.isVertex,
        tetEdge,
        tetVertex,
        faceVertices: [...fv],
        intervalIdx: -1,
      });
    }
  }

  // Deduplicate edge punctures (same tet edge from two faces)
  const deduped = deduplicatePunctures(rawPunctures);

  // ── Phase 4: Infinity punctures (λ=∞ and λ=0 from leading/constant coeffs) ──
  addInfinityPunctures(deduped, cc);

  cc.punctures = deduped;

  // ── Phase 5: Build intervals from Q roots ──
  cc.intervals = buildIntervals(qRoots);

  // ── Phase 6: Assign punctures to intervals ──
  assignPuncturesToIntervals(cc);

  // ── Phase 7: Waypoint filtering and T-category ──
  const tags = [];

  // Recount interval occupancy excluding waypoints
  for (const iv of cc.intervals) iv.n_pv = 0;
  for (const pi of cc.punctures) {
    if (!isWaypoint(pi) && pi.intervalIdx >= 0)
      cc.intervals[pi.intervalIdx].n_pv++;
  }

  let nFace = 0;
  for (const pi of cc.punctures) {
    if (!isWaypoint(pi)) nFace++;
  }

  // Occupancy tuple (sorted non-zero n_pv counts)
  const occTuple = cc.intervals
    .map(iv => iv.n_pv)
    .filter(x => x > 0)
    .sort((a, b) => a - b);

  let tCat = `T${nFace}`;
  if (occTuple.length > 1) tCat += `_(${occTuple.join(',')})`;
  cc.category = tCat + '_' + qType;

  // ── Phase 8: Shared-root detection ──
  cc.has_shared_root = hasSharedRoot_bigint(Q_bigint, P_bigint);
  if (cc.has_shared_root) {
    tags.push('SR');
    findSRLocation(cc);
  }

  // ── Phase 9: Critical-point degeneracies ──
  detectCriticalPoints(cc, tags);
  detectDmd(cc, tags);
  detectBubble(cc, tags);

  if (tags.length > 0) cc.category += '_' + tags.join('_');

  // ── Phase 10: Stitching (pairing) ──
  cc.pairs = stitchPunctures(cc);

  return cc;
}

// ═══════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════

function isWaypoint(p) {
  return (p.isEdge || p.isVertex) &&
         (p.lambda === 0.0 || !isFinite(p.lambda));
}

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

/** Deduplicate edge and vertex punctures. */
function deduplicatePunctures(punctures) {
  // Vertex dedup: same tet vertex from multiple faces → keep lex-smallest face
  const vtxGroups = new Map();
  for (let i = 0; i < punctures.length; i++) {
    const p = punctures[i];
    if (!p.isVertex || p.tetVertex < 0) continue;
    if (!vtxGroups.has(p.tetVertex)) vtxGroups.set(p.tetVertex, []);
    vtxGroups.get(p.tetVertex).push(i);
  }
  const toRemove = new Set();
  for (const [, group] of vtxGroups) {
    if (group.length <= 1) continue;
    let canonical = group[0];
    let canFace = [...punctures[canonical].faceVertices].sort((a, b) => a - b);
    for (const pi of group) {
      const f = [...punctures[pi].faceVertices].sort((a, b) => a - b);
      if (lexLess(f, canFace)) { canonical = pi; canFace = f; }
    }
    for (const pi of group) if (pi !== canonical) toRemove.add(pi);
  }

  // Edge dedup: same tet edge from multiple faces → keep first
  const edgeSeen = new Set();
  for (let i = 0; i < punctures.length; i++) {
    if (toRemove.has(i)) continue;
    const p = punctures[i];
    if (!p.isEdge || !p.tetEdge) continue;
    const key = `${p.tetEdge[0]},${p.tetEdge[1]}`;
    if (edgeSeen.has(key)) toRemove.add(i);
    else edgeSeen.add(key);
  }

  return punctures.filter((_, i) => !toRemove.has(i));
}

function lexLess(a, b) {
  for (let i = 0; i < Math.min(a.length, b.length); i++) {
    if (a[i] < b[i]) return true;
    if (a[i] > b[i]) return false;
  }
  return a.length < b.length;
}

/**
 * Add punctures at λ=∞ and λ=0 from polynomial leading/constant coefficients.
 * Only adds when P_k[deg]==0 for some k (point is on tet face/edge/vertex).
 * Mirrors ftk2: punctures exist only where the curve crosses the tet boundary.
 */
function addInfinityPunctures(punctures, cc) {
  const { Q_coeffs: Q, P_coeffs: P, Q_bigint, P_bigint } = cc;

  // λ=∞: check each face k where P_k[3]==0 (exact integer check)
  if (Q[3] !== 0) {
    const q3 = Q_bigint[3];
    for (let k = 0; k < 4; k++) {
      if (P_bigint[k][3] !== 0n) continue; // μ_k(∞) ≠ 0 → not on face k
      // Check remaining coords: P_j[3]*Q[3] >= 0 for j≠k
      let inside = true;
      for (let j = 0; j < 4; j++) {
        if (j === k) continue;
        if (P_bigint[j][3] * q3 < 0n) { inside = false; break; }
      }
      if (!inside) continue;
      // Already found by face solver?
      if (punctures.some(p => !isFinite(p.lambda))) continue;
      // Determine edge/vertex
      const nZero = P_bigint.filter(pk => pk[3] === 0n).length;
      const isVertex = nZero >= 3;
      const isEdge = nZero >= 2 && !isVertex;
      const mu = P.map(pk => pk[3] / Q[3]);
      const muClip = mu.map(m => Math.max(0, m));
      const s = muClip.reduce((a, b) => a + b, 0);
      if (s < 1e-10) continue;
      const muNorm = muClip.map(m => m / s);
      punctures.push({
        face: -1, lambda: Infinity, bary: [0, 0, 0],
        mu: muNorm, pos3d: baryTetTo3D(muNorm),
        isEdge, isVertex, tetEdge: null, tetVertex: -1,
        faceVertices: [], intervalIdx: -1,
      });
      break; // one puncture per λ=∞
    }
  }

  // λ=0: check each face k where P_k[0]==0 (exact integer check)
  if (Q[0] !== 0) {
    const q0 = Q_bigint[0];
    for (let k = 0; k < 4; k++) {
      if (P_bigint[k][0] !== 0n) continue; // μ_k(0) ≠ 0 → not on face k
      let inside = true;
      for (let j = 0; j < 4; j++) {
        if (j === k) continue;
        if (P_bigint[j][0] * q0 < 0n) { inside = false; break; }
      }
      if (!inside) continue;
      if (punctures.some(p => p.lambda === 0.0)) continue;
      const nZero = P_bigint.filter(pk => pk[0] === 0n).length;
      const isVertex = nZero >= 3;
      const isEdge = nZero >= 2 && !isVertex;
      const mu = P.map(pk => pk[0] / Q[0]);
      const muClip = mu.map(m => Math.max(0, m));
      const s = muClip.reduce((a, b) => a + b, 0);
      if (s < 1e-10) continue;
      const muNorm = muClip.map(m => m / s);
      punctures.push({
        face: -1, lambda: 0.0, bary: [0, 0, 0],
        mu: muNorm, pos3d: baryTetTo3D(muNorm),
        isEdge, isVertex, tetEdge: null, tetVertex: -1,
        faceVertices: [], intervalIdx: -1,
      });
      break;
    }
  }
}

/** Assign each puncture to its Q-interval via Sturm counts or sign-of-Q. */
function assignPuncturesToIntervals(cc) {
  const { Q_coeffs: Q, Q_roots: roots, intervals, punctures, qDegree } = cc;

  if (intervals.length <= 1 || qDegree === 0) {
    // Single interval: all punctures in interval 0
    for (const p of punctures) p.intervalIdx = 0;
    return;
  }

  if (cc.n_Q_roots === 1) {
    // Q3- or Q1: sign-of-Q method
    const signQ3 = Q[qDegree] > 0 ? 1 : -1;
    for (const p of punctures) {
      if (!isFinite(p.lambda)) {
        // λ=∞ → rightmost interval
        p.intervalIdx = intervals.length - 1;
        continue;
      }
      const qVal = polyEval(Q, p.lambda);
      const signQ = qVal > 0 ? 1 : qVal < 0 ? -1 : 0;

      if (signQ === 0) {
        // At a Q-root boundary — perturb
        const delta = 4 * Number.EPSILON * Math.max(1, Math.abs(p.lambda));
        const qPerturbed = polyEval(Q, p.lambda + delta);
        p.intervalIdx = (qPerturbed > 0) === (signQ3 > 0) ? 1 : 0;
      } else {
        p.intervalIdx = (signQ > 0) === (signQ3 > 0) ? 1 : 0;
      }
    }
    return;
  }

  // Q3+ or Q2: Sturm sequence method
  const Q_d = gcdNormalizePoly(cc.Q_bigint);
  let degQ_d = 3;
  while (degQ_d > 0 && Q_d[degQ_d] === 0) degQ_d--;
  if (degQ_d === 0) {
    for (const p of punctures) p.intervalIdx = 0;
    return;
  }

  const seq = buildSturmCubic(Q_d);

  for (const p of punctures) {
    if (!isFinite(p.lambda)) {
      p.intervalIdx = intervals.length - 1;
      continue;
    }
    const { count, certified } = sturmCountCubicCertified(seq, p.lambda);
    let cnt;
    if (certified) {
      cnt = count;
    } else {
      const delta = 4 * Number.EPSILON * Math.max(1, Math.abs(p.lambda));
      cnt = sturmCountCubic(seq, p.lambda + delta);
    }
    p.intervalIdx = Math.max(0, Math.min(intervals.length - 1, cc.n_Q_roots - cnt));
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Degeneracy detection
// ═══════════════════════════════════════════════════════════════════════

/** Detect Cv/Cw critical points (interior + face/edge/vertex). */
function detectCriticalPoints(cc, tags) {
  const { V, W, punctures } = cc;

  // Interior: check if v=0 or w=0 inside tet
  const cvMu = checkFieldZeroInTet(V);
  const cwMu = checkFieldZeroInTet(W);

  let hasCv = cvMu !== null;
  let hasCw = cwMu !== null;
  let hasC2v = false, hasC1v = false, hasC0v = false;
  let hasC2w = false, hasC1w = false, hasC0w = false;

  if (hasCv) { cc.has_Cv_pos = true; cc.Cv_mu = cvMu; }
  if (hasCw) { cc.has_Cw_pos = true; cc.Cw_mu = cwMu; }

  // Scan punctures for face/edge/vertex at critical λ
  for (const pi of punctures) {
    const isLam0 = pi.lambda === 0.0;
    const isLamInf = !isFinite(pi.lambda);

    if (isLam0) {
      if (pi.isVertex) hasC0v = true;
      else if (pi.isEdge) hasC1v = true;
      else hasC2v = true;
    }
    if (isLamInf) {
      if (pi.isVertex) hasC0w = true;
      else if (pi.isEdge) hasC1w = true;
      else hasC2w = true;
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

  // D01: edge puncture at generic λ (not critical)
  let hasD01 = false;
  for (const pi of punctures) {
    if (pi.isEdge && pi.lambda !== 0.0 && isFinite(pi.lambda)) hasD01 = true;
  }

  if (hasD00) tags.push('D00');
  if (hasD01) tags.push('D01');
}

/** Detect bubble (closed PV curve with no punctures). */
function detectBubble(cc, tags) {
  const { Q_coeffs: Q, P_coeffs: P, n_Q_roots } = cc;

  const nFace = cc.punctures.filter(p => !isWaypoint(p)).length;
  if (nFace > 0) return;
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

  // Derivative of polynomial: d/dλ (c0 + c1λ + c2λ² + c3λ³) = c1 + 2c2λ + 3c3λ²
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

        // L'Hôpital: μ_j(root) = P_j'(root)/Q'(root) for shared component,
        // and nearby evaluation for others. Only valid if result is inside tet.
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

        // SR root exists but curve is not inside tet here — no marker
        cc.SR_pos3d = null;
        return;
      }
    }
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Stitching: pair punctures into curve segments
// ═══════════════════════════════════════════════════════════════════════

/**
 * Full stitching algorithm matching ftk2 pv_tet_case_finder.cu.
 *
 * 1. Group non-waypoint punctures by interval
 * 2. Sort within each interval by λ
 * 3. Infinity merging (check asymptotic point inside tet)
 * 4. Pair (0,1),(2,3),... within each group
 * 5. SR pass-through for unpaired remainders
 */
function stitchPunctures(cc) {
  const { punctures, Q_coeffs: Q, P_coeffs: P, Q_bigint, P_bigint } = cc;
  const pairs = [];

  // Step 1: Group non-waypoint punctures by interval
  const ivPuncs = new Map(); // intervalIdx → [puncture indices]
  for (let i = 0; i < punctures.length; i++) {
    if (isWaypoint(punctures[i])) continue;
    const iv = punctures[i].intervalIdx;
    if (iv < 0) continue;
    if (!ivPuncs.has(iv)) ivPuncs.set(iv, []);
    ivPuncs.get(iv).push(i);
  }

  // Step 2: Sort within each interval by λ (∞ sorts to end)
  for (const [, pis] of ivPuncs) {
    pis.sort((a, b) => {
      const la = punctures[a].lambda, lb = punctures[b].lambda;
      if (!isFinite(la) && !isFinite(lb)) return 0;
      if (!isFinite(la)) return 1;
      if (!isFinite(lb)) return -1;
      return la - lb;
    });
  }

  // Step 3: Infinity merging
  const nIntervals = cc.intervals.length;
  const leftIdx = 0;
  const rightIdx = nIntervals - 1;
  let mergedInfinity = false;
  // Track original interval membership for is_cross detection (matching ftk2)
  const rightOriginSet = new Set();
  const leftOriginSet = new Set();

  if (nIntervals >= 2 && (ivPuncs.has(leftIdx) || ivPuncs.has(rightIdx))) {
    // Check if asymptotic point is inside tet: sign(P_k[3]) matches sign(Q[3])
    let mergeOk;
    if (Q[3] === 0) {
      mergeOk = true;
    } else {
      mergeOk = true;
      const q3 = Q_bigint[3];
      for (let k = 0; k < 4; k++) {
        if (P_bigint[k][3] * q3 < 0n) { mergeOk = false; break; }
      }
    }

    if (mergeOk) {
      const rightPis = ivPuncs.get(rightIdx) || [];
      const leftPis = ivPuncs.get(leftIdx) || [];
      for (const pi of rightPis) rightOriginSet.add(pi);
      for (const pi of leftPis) leftOriginSet.add(pi);
      // Merge: right then left (projective order through infinity)
      const merged = [...rightPis, ...leftPis];
      ivPuncs.set(rightIdx, merged);
      ivPuncs.delete(leftIdx);
      mergedInfinity = true;
    }
  }

  // Step 4: Pair within each interval group
  const unpaired = [];
  for (const [ivIdx, pis] of ivPuncs) {
    for (let j = 0; j + 1 < pis.length; j += 2) {
      const pi_a = pis[j], pi_b = pis[j + 1];
      // is_cross: one puncture from original right interval, other from left
      let isCross = false;
      if (mergedInfinity && ivIdx === rightIdx) {
        isCross = (rightOriginSet.has(pi_a) && leftOriginSet.has(pi_b)) ||
                  (leftOriginSet.has(pi_a) && rightOriginSet.has(pi_b));
      }
      pairs.push({ pi_a, pi_b, is_cross: isCross, interval_idx: ivIdx });
    }
    if (pis.length % 2 === 1) {
      unpaired.push(pis[pis.length - 1]);
    }
  }

  // Step 5: SR pass-through pairing
  if (unpaired.length >= 2) {
    unpaired.sort((a, b) => {
      const la = punctures[a].lambda, lb = punctures[b].lambda;
      if (!isFinite(la)) return 1;
      if (!isFinite(lb)) return -1;
      return la - lb;
    });
    for (let j = 0; j + 1 < unpaired.length; j += 2) {
      pairs.push({
        pi_a: unpaired[j], pi_b: unpaired[j + 1],
        is_cross: false, interval_idx: -1,
      });
    }
  }

  return pairs;
}
