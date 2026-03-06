/**
 * Exact PV triangle solver.
 * Ported from ftk2 solve_pv_triangle.
 *
 * Given V[3][3], W[3][3] (integer field values at 3 triangle vertices)
 * and SoS vertex indices, finds puncture points where v × w = 0.
 */
import {
  triangleCharPoly, discriminantSign_bigint, allParallel,
} from './bigint_poly.js';
import {
  buildSturmCubic, tightenRootInterval,
  buildSturmDeg4, sturmCountDeg4, computeBaryNumerators,
} from './sturm.js';

const EPS = Number.EPSILON;

/**
 * Solve cubic with exact discriminant sign (SoS-aware).
 * P_float[4] = float polynomial, P_bigint[4] = BigInt polynomial.
 * indices[3] = SoS vertex indices (or null).
 * Returns array of float roots.
 */
function solveCubicSoS(P_float, P_bigint, indices) {
  const P = P_float;

  if (P[3] === 0) {
    if (P[2] === 0) {
      if (P[1] === 0) return [];
      return [-P[0] / P[1]];
    }
    const disc = P[1] * P[1] - 4 * P[2] * P[0];
    if (disc < 0) return [];
    if (disc === 0) return [-P[1] / (2 * P[2])];
    const sq = Math.sqrt(disc);
    return [(-P[1] + sq) / (2 * P[2]), (-P[1] - sq) / (2 * P[2])];
  }

  const b = P[2] / P[3], c = P[1] / P[3], d = P[0] / P[3];

  // Special case: P[0]==P[1]==0 → λ²(P[2]+P[3]λ)
  if (c === 0 && d === 0) {
    const minIdx = indices ? Math.min(indices[0], indices[1], indices[2]) : 0;
    if (minIdx % 2 === 0) return [-b]; // tangent excluded
    if (-b === 0) return [0];
    return [-b, 0, 0];
  }

  const q = (3 * c - b * b) / 9;
  const r = (-27 * d + b * (9 * c - 2 * b * b)) / 54;
  const disc = q * q * q + r * r;
  const term1 = b / 3;

  // Use exact discriminant
  let exactSign = P_bigint ? discriminantSign_bigint(P_bigint) : 0;
  if (exactSign === 0) {
    if (disc > 0) exactSign = -1;
    else if (disc < 0) exactSign = 1;
  }

  if (exactSign < 0) {
    // 1 real root (Cardano)
    const safe = disc > 0 ? disc : 0;
    const sq = Math.sqrt(safe);
    let s = r + sq; s = s < 0 ? -Math.cbrt(-s) : Math.cbrt(s);
    let t = r - sq; t = t < 0 ? -Math.cbrt(-t) : Math.cbrt(t);
    return [-term1 + s + t];
  } else if (exactSign > 0) {
    // 3 real roots (trigonometric)
    let mq = -q > 0 ? -q : 0;
    if (mq === 0) {
      const safe = disc > 0 ? disc : 0;
      const sq = Math.sqrt(safe);
      let s = r + sq; s = s < 0 ? -Math.cbrt(-s) : Math.cbrt(s);
      let t = r - sq; t = t < 0 ? -Math.cbrt(-t) : Math.cbrt(t);
      return [-term1 + s + t];
    }
    let arg = r / Math.sqrt(mq * mq * mq);
    arg = Math.max(-1, Math.min(1, arg));
    const dum1 = Math.acos(arg);
    const r13 = 2 * Math.sqrt(mq);
    return [
      -term1 + r13 * Math.cos(dum1 / 3),
      -term1 + r13 * Math.cos((dum1 + 2 * Math.PI) / 3),
      -term1 + r13 * Math.cos((dum1 + 4 * Math.PI) / 3),
    ];
  } else {
    // Repeated root → SoS
    const minIdx = indices ? Math.min(indices[0], indices[1], indices[2]) : 0;
    const r13 = r < 0 ? -Math.cbrt(-r) : Math.cbrt(r);
    const root0 = -term1 + 2 * r13;
    const root1 = -(r13 + term1);
    if (minIdx % 2 === 0) return [root0]; // tangent excluded
    if (r === 0) return [root0];
    return [root0, root1, root1];
  }
}

/**
 * Full triangle PV solver.
 *
 * @param V_face - V[3][3] integer field values at face vertices
 * @param W_face - W[3][3] integer field values at face vertices
 * @param indices - [3] global vertex indices for SoS (integers)
 * @returns Array of puncture objects { lambda, bary[3], faceVertices }
 */
export function solvePVTriangle(V_face, W_face, indices) {
  // 1. All-parallel check (exact)
  if (allParallel(V_face, W_face)) return null; // entire face is PV

  // 2. Exact integer characteristic polynomial
  const { P_bigint, P_float } = triangleCharPoly(V_face, W_face);

  // 3. Solve cubic with exact discriminant
  const floatRoots = solveCubicSoS(P_float, P_bigint, indices);

  // 4. Sturm root isolation
  const intervals = [];
  const seq = buildSturmCubic(P_float);
  for (const rf of floatRoots) {
    const result = tightenRootInterval(seq, rf);
    intervals.push(result);
  }

  // 5. Compute barycentric numerator polynomials N[3][5], D[5]
  // Build Mlin_q and blin_q from integer field values
  // VT[component][vertex], WT[component][vertex]
  const Mlin_q = Array.from({ length: 3 }, () =>
    Array.from({ length: 2 }, () => [0, 0])
  );
  const blin_q = Array.from({ length: 3 }, () => [0, 0]);

  for (let r = 0; r < 3; r++) {
    for (let c = 0; c < 2; c++) {
      Mlin_q[r][c][0] = V_face[c][r] - V_face[2][r]; // VT[r][c] - VT[r][2]
      Mlin_q[r][c][1] = -(W_face[c][r] - W_face[2][r]);
    }
    blin_q[r][0] = -V_face[2][r]; // -VT[r][2]
    blin_q[r][1] = W_face[2][r];
  }

  const { N: N_poly, D: D_poly } = computeBaryNumerators(Mlin_q, blin_q);

  // Pre-build Sturm sequence for D
  let degD = 4;
  while (degD > 0 && D_poly[degD] === 0.0) degD--;
  const seq_D = buildSturmDeg4(D_poly, degD);

  // Exact λ=0 exclusion
  const zeroIsExactRoot = P_bigint[0] === 0n;

  // 6. For each isolated root, certify barycentric signs
  const punctures = [];

  for (let i = 0; i < intervals.length; i++) {
    const { lo: lambdaLo, hi: lambdaHi, success: haveInterval } = intervals[i];
    const lambda = (lambdaLo + lambdaHi) * 0.5;

    // Skip exact λ=0 root
    if (zeroIsExactRoot && lambdaLo <= 0 && 0 <= lambdaHi) continue;

    // Certification window
    let winLo, winHi;
    if (haveInterval && lambdaLo < lambdaHi) {
      winLo = lambdaLo;
      winHi = lambdaHi;
    } else {
      const delta = Math.max(Math.abs(lambda) * 4 * EPS, Number.MIN_VALUE);
      winLo = lambda - delta;
      winHi = lambda + delta;
    }

    // Certify sign of each N_k
    const bsign = [0, 0, 0];
    let nBnd = 0;
    for (let k = 0; k < 3; k++) {
      bsign[k] = tryCertifyNkSign(k, winLo, winHi, N_poly, seq_D);
      if (bsign[k] === 0) nBnd++;
    }

    // Any certified negative → outside simplex
    if (bsign.some(s => s < 0)) continue;

    // SoS ownership rules
    if (nBnd === 2) {
      // G2: vertex puncture
      if (!indices) continue;
      const m = bsign.findIndex(s => s !== 0);
      if (m < 0) continue;
      const a = (m + 1) % 3, b = (m + 2) % 3;
      if (!(indices[m] < Math.min(indices[a], indices[b]))) continue;
    } else {
      let reject = false;
      for (let k = 0; k < 3; k++) {
        if (bsign[k] > 0) continue;
        // bsign[k] === 0: boundary → apply SoS
        if (!indices) { reject = true; break; }
        const ii = (k + 1) % 3, jj = (k + 2) % 3;
        if (!(indices[k] < Math.min(indices[ii], indices[jj]))) {
          reject = true; break;
        }
      }
      if (reject) continue;
    }

    // Compute barycentric coordinates
    const lamD = lambda;
    function evalPolySturm(c, deg, x) {
      let r = c[deg];
      for (let d = deg - 1; d >= 0; d--) r = r * x + c[d];
      return r;
    }
    const dVal = evalPolySturm(D_poly, 4, lamD);
    const nu = [0, 0, 0];
    for (let k = 0; k < 3; k++) {
      const nkVal = evalPolySturm(N_poly[k], 4, lamD);
      nu[k] = dVal > 0 ? nkVal / dVal : 0;
    }

    // Detect edge/vertex puncture
    const isEdge = nBnd >= 1;
    const isVertex = nBnd >= 2;

    punctures.push({ lambda, bary: nu, isEdge, isVertex });
  }

  return punctures;
}

/**
 * Try to certify the sign of N_k(λ*) in [lo, hi].
 * Returns +1/-1 if certified, 0 if uncertain (SoS needed).
 */
function tryCertifyNkSign(k, lo, hi, N_poly, seq_D) {
  let degNk = 4;
  while (degNk > 0 && N_poly[k][degNk] === 0.0) degNk--;

  const seq_nk = buildSturmDeg4(N_poly[k], degNk);
  if (sturmCountDeg4(seq_nk, lo) - sturmCountDeg4(seq_nk, hi) !== 0)
    return 0; // N_k has a root in window

  // Check D
  if (sturmCountDeg4(seq_D, lo) - sturmCountDeg4(seq_D, hi) !== 0)
    return 0; // D has a root

  // Certified Horner sign at lo
  function evalPolySturm(c, deg, x) {
    let r = c[deg];
    for (let d = deg - 1; d >= 0; d--) r = r * x + c[d];
    return r;
  }

  const nkLo = evalPolySturm(N_poly[k], degNk, lo);
  const ax = Math.abs(lo);
  let condNk = Math.abs(N_poly[k][degNk]);
  for (let d = degNk - 1; d >= 0; d--)
    condNk = condNk * ax + Math.abs(N_poly[k][d]);

  const evalGamma = (2 * degNk + 2) * EPS;
  if (Math.abs(nkLo) > evalGamma * condNk)
    return nkLo > 0 ? 1 : -1;

  return 0;
}
