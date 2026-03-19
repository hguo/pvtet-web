/**
 * 2D Triangle PV Classification.
 * Ported from ftk2 pv_tri_classify_2d.hpp classify_case_v2_2d().
 * All topological decisions use exact BigInt arithmetic.
 */
import { solvePVTri2D } from './exactpv2_2d.js';
import {
  effectiveDegree, resultantSign, polyGcdFull, compareRoots,
} from './exactpv2.js';
import { characteristicPolynomials2D_bigint } from './bigint_poly.js';
import { computeTriQP2D_float } from './polynomials.js';
import { solveQuadratic } from './roots.js';

/**
 * Check if 2D field F[3][2] vanishes inside triangle.
 * Returns 0=not found, 1=interior, 2=edge, 3=vertex.
 */
export function checkFieldZeroInTri2D(F) {
  for (let i = 0; i < 3; i++)
    if (F[i][0] === 0 && F[i][1] === 0) return 3;

  const a00 = F[0][0] - F[2][0], a01 = F[1][0] - F[2][0];
  const a10 = F[0][1] - F[2][1], a11 = F[1][1] - F[2][1];
  const b0 = -F[2][0], b1 = -F[2][1];

  const det = a00 * a11 - a01 * a10;
  if (det === 0) {
    let ui = -1;
    for (let i = 0; i < 3; i++)
      if (F[i][0] !== 0 || F[i][1] !== 0) { ui = i; break; }
    if (ui < 0) return 1;
    let hasPos = false, hasNeg = false;
    for (let i = 0; i < 3; i++) {
      const dot = F[i][0] * F[ui][0] + F[i][1] * F[ui][1];
      if (dot > 0) hasPos = true;
      if (dot < 0) hasNeg = true;
    }
    if (hasPos && hasNeg) return 2;
    return 0;
  }

  const n0 = b0 * a11 - a01 * b1;
  const n1 = a00 * b1 - b0 * a10;
  const n2 = det - n0 - n1;

  if (det > 0) {
    if (n0 > 0 && n1 > 0 && n2 > 0) return 1;
    if (n0 >= 0 && n1 >= 0 && n2 >= 0) {
      const nz = (n0 === 0 ? 1 : 0) + (n1 === 0 ? 1 : 0) + (n2 === 0 ? 1 : 0);
      if (nz >= 2) return 3;
      if (nz === 1) return 2;
    }
  } else {
    if (n0 < 0 && n1 < 0 && n2 < 0) return 1;
    if (n0 <= 0 && n1 <= 0 && n2 <= 0) {
      const nz = (n0 === 0 ? 1 : 0) + (n1 === 0 ? 1 : 0) + (n2 === 0 ? 1 : 0);
      if (nz >= 2) return 3;
      if (nz === 1) return 2;
    }
  }
  return 0;
}

/**
 * Classify a 2D triangle PV case.
 * @param {number[][]} V - V[3][2]
 * @param {number[][]} W - W[3][2]
 */
export function classifyTriCase2D(V, W) {
  const cc = {
    V, W,
    dim: '2d',
    category: '',
    Q_coeffs: [0, 0, 0],
    P_coeffs: [[0,0,0],[0,0,0],[0,0,0]],
    Q_bigint: [0n, 0n, 0n],
    P_bigint: [[0n,0n,0n],[0n,0n,0n],[0n,0n,0n]],
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

  // Phase 1: Q, P polynomials (BigInt + float)
  const { Q_bigint, P_bigint } = characteristicPolynomials2D_bigint(V, W);
  cc.Q_bigint = Q_bigint;
  cc.P_bigint = P_bigint;

  const { Q: Qf, P: Pf } = computeTriQP2D_float(V, W);
  cc.Q_coeffs = Qf;
  cc.P_coeffs = Pf;

  // Phase 2: Q structure (exact BigInt degree, matching C++)
  let degQ = effectiveDegree(Q_bigint, 2);
  cc.qDegree = degQ;

  if (degQ === 2) {
    const disc2 = Q_bigint[1] * Q_bigint[1] - 4n * Q_bigint[0] * Q_bigint[2];
    cc.Q_disc_sign = disc2 > 0n ? 1 : disc2 < 0n ? -1 : 0;
    cc.n_Q_roots = cc.Q_disc_sign > 0 ? 2 : cc.Q_disc_sign === 0 ? 1 : 0;
  } else if (degQ === 1) {
    cc.Q_disc_sign = 0;
    cc.n_Q_roots = 1;
  } else {
    cc.Q_disc_sign = 0;
    cc.n_Q_roots = 0;
  }

  const qRoots = solveQuadratic(Qf);
  qRoots.sort((a, b) => a - b);
  cc.Q_roots = qRoots;

  // Q type string
  let qType;
  if (degQ === 0 && Q_bigint[0] === 0n) qType = 'Qz';
  else if (degQ === 0) qType = 'Q0';
  else if (degQ === 1) qType = 'Q1';
  else {
    if (cc.Q_disc_sign > 0) qType = 'Q2+';
    else if (cc.Q_disc_sign < 0) qType = 'Q2-';
    else qType = 'Q2o';
  }

  // Phase 3: ExactPV2 2D solver
  const pv2 = solvePVTri2D(Q_bigint, P_bigint);

  // Map punctures to display format
  const { polyEval } = await_polyEval();
  const punctures = [];
  for (const ep of pv2.punctures) {
    let lambda, mu;
    if (ep.isInfinity) {
      lambda = Infinity;
      if (Qf[degQ] !== 0) {
        mu = Pf.map(pk => {
          const dPk = pk.length - 1;
          return dPk >= degQ ? pk[degQ] / Qf[degQ] : 0;
        });
      } else {
        mu = [1/3, 1/3, 1/3];
      }
    } else {
      const roots = solveQuadratic(Pf[ep.face]).sort((a, b) => a - b);
      lambda = (ep.rootIdx < roots.length) ? roots[ep.rootIdx] : (roots[roots.length - 1] || 0);
      const qVal = Qf[0] + Qf[1] * lambda + Qf[2] * lambda * lambda;
      if (Math.abs(qVal) > 1e-30) {
        mu = Pf.map(pk => (pk[0] + pk[1] * lambda + pk[2] * lambda * lambda) / qVal);
      } else {
        mu = [1/3, 1/3, 1/3];
      }
    }

    const muClip = mu.map(m => Math.max(0, m));
    const s = muClip.reduce((a, b) => a + b, 0);
    const muNorm = s > 1e-10 ? muClip.map(m => m / s) : [1/3, 1/3, 1/3];

    let triVertex = -1;
    if (ep.isEdge || ep.isVertex) {
      const ef0 = ep.edgeFaces[0], ef1 = ep.edgeFaces[1];
      for (let v = 0; v < 3; v++)
        if (v !== ef0 && v !== ef1) triVertex = v;
    }

    // Compute h_below: count passthrough roots below this puncture (matching C++)
    let hBelow = 0;
    if (pv2.hNRoots > 0) {
      if (ep.rootIdx < 0) {
        hBelow = pv2.hNRoots; // infinity puncture: all h roots below
      } else {
        for (let hi = 0; hi < pv2.hNRoots; hi++) {
          const cmp = compareRoots(
            pv2.P_red[ep.face], pv2.degP_red[ep.face],
            pv2.nDistinctRed[ep.face], ep.rootIdx,
            pv2.h, pv2.hDeg, pv2.hNRoots, hi);
          if (cmp > 0) hBelow++;
        }
      }
    }
    const intervalIdx = ep.qInterval + hBelow;

    punctures.push({
      face: ep.face,
      lambda,
      bary: [...muNorm],
      mu: muNorm,
      isEdge: ep.isEdge,
      isVertex: ep.isVertex,
      triVertex,
      intervalIdx,
      edgeFaces: ep.edgeFaces,
    });
  }
  cc.punctures = punctures;

  // Phase 4: Build intervals (from exact BigInt root count, matching C++)
  // C++ builds intervals from n_Q_roots, not from float root values.
  if (cc.n_Q_roots > 0) {
    cc.intervals = [{ lb: null, ub: null, isInfinity: true, n_pv: 0 }];
    for (let i = 0; i + 1 < cc.n_Q_roots; i++)
      cc.intervals.push({ lb: null, ub: null, isInfinity: false, n_pv: 0 });
    cc.intervals.push({ lb: null, ub: null, isInfinity: true, n_pv: 0 });
  } else {
    cc.intervals = [{ lb: null, ub: null, isInfinity: true, n_pv: 0 }];
  }

  // Fill interval bounds from float roots for display
  if (qRoots.length > 0 && cc.intervals.length > 1) {
    cc.intervals[0].ub = qRoots[0];
    for (let i = 0; i < qRoots.length - 1 && i + 1 < cc.intervals.length; i++) {
      cc.intervals[i + 1].lb = qRoots[i];
      cc.intervals[i + 1].ub = qRoots[i + 1];
    }
    cc.intervals[cc.intervals.length - 1].lb = qRoots[qRoots.length - 1];
  }

  // Count punctures per interval
  for (const pi of punctures) {
    if (pi.intervalIdx >= 0 && pi.intervalIdx < cc.intervals.length)
      cc.intervals[pi.intervalIdx].n_pv++;
  }

  // Phase 5: T-category
  const tags = [];
  const nPunc = punctures.length;
  const occMap = new Map();
  for (const pi of punctures) {
    if (pi.intervalIdx >= 0)
      occMap.set(pi.intervalIdx, (occMap.get(pi.intervalIdx) || 0) + 1);
  }
  const occTuple = [...occMap.values()].sort((a, b) => a - b);
  let tCat = `T${nPunc}`;
  if (occTuple.length > 1) tCat += `_(${occTuple.join(',')})`;
  cc.category = tCat + '_' + qType;

  // Phase 6: SR detection
  if (degQ > 0) {
    const degQi = effectiveDegree(Q_bigint, 2);
    for (let k = 0; k < 3; k++) {
      const degPk = effectiveDegree(P_bigint[k], 2);
      if (degPk <= 0) continue;
      if (resultantSign(Q_bigint, degQi, P_bigint[k], degPk) === 0) {
        cc.has_shared_root = true;
        break;
      }
    }
  }

  // SR L'Hôpital inside-triangle check
  if (cc.has_shared_root) {
    const degQi = effectiveDegree(Q_bigint, 2);
    let anyInside = false;
    for (let k = 0; k < 3; k++) {
      const degPk = effectiveDegree(P_bigint[k], 2);
      if (degPk <= 0) continue;
      if (resultantSign(Q_bigint, degQi, P_bigint[k], degPk) !== 0) continue;
      const gRes = polyGcdFull(Q_bigint, degQi, P_bigint[k], degPk);
      if (gRes.dh === 1) {
        const g0 = gRes.h[0], g1 = gRes.h[1];
        const denom = g1 * Q_bigint[1] - 2n * Q_bigint[2] * g0;
        if (denom !== 0n) {
          let inside = true;
          for (let j = 0; j < 3; j++) {
            const muJ = g1 * P_bigint[j][1] - 2n * P_bigint[j][2] * g0;
            if (muJ * denom < 0n) { inside = false; break; }
          }
          if (inside) anyInside = true;
        }
      } else if (gRes.dh >= 2) {
        anyInside = true;
      }
    }
    if (pv2.hDeg === 1) {
      const h0 = pv2.h[0], h1 = pv2.h[1];
      const denom = h1 * Q_bigint[1] - 2n * Q_bigint[2] * h0;
      if (denom !== 0n) {
        let inside = true;
        for (let j = 0; j < 3; j++) {
          const muJ = h1 * P_bigint[j][1] - 2n * P_bigint[j][2] * h0;
          if (muJ * denom < 0n) { inside = false; break; }
        }
        if (inside) anyInside = true;
      }
    } else if (pv2.hDeg >= 2) {
      anyInside = true;
    }
    if (!anyInside) cc.has_shared_root = false;
  }

  if (cc.has_shared_root) {
    const isISR = pv2.passthroughDeg >= 2;
    if (isISR) tags.push('ISR');
    else tags.push('SR');
  }

  // Phase 7: Cv/Cw detection
  {
    let denom = Q_bigint[0];
    let mu = [0n, 0n, 0n];
    if (denom === 0n) {
      denom = Q_bigint[1];
      for (let k = 0; k < 3; k++) mu[k] = P_bigint[k][1];
    } else {
      for (let k = 0; k < 3; k++) mu[k] = P_bigint[k][0];
    }
    if (denom !== 0n) {
      let inside = true, nz = 0;
      for (let k = 0; k < 3; k++) {
        if (mu[k] * denom < 0n) { inside = false; break; }
        if (mu[k] === 0n) nz++;
      }
      if (inside) {
        cc.has_Cv_pos = true;
        if (nz >= 2) tags.push('Cv0');
        else if (nz === 1) tags.push('Cv1');
        else tags.push('Cv');
      }
    } else {
      const cvRes = checkFieldZeroInTri2D(V);
      if (cvRes > 0) {
        cc.has_Cv_pos = true;
        if (cvRes === 3) tags.push('Cv0');
        else if (cvRes === 2) tags.push('Cv1');
        else tags.push('Cv');
      }
    }
  }
  {
    const dQ = effectiveDegree(Q_bigint, 2);
    let cwValid = dQ > 0;
    if (cwValid) {
      for (let k = 0; k < 3; k++) {
        if (effectiveDegree(P_bigint[k], 2) > dQ) { cwValid = false; break; }
      }
    }
    if (cwValid) {
      const denom = Q_bigint[dQ];
      let inside = true, nz = 0;
      for (let k = 0; k < 3; k++) {
        const muK = (effectiveDegree(P_bigint[k], 2) === dQ) ? P_bigint[k][dQ] : 0n;
        if (muK * denom < 0n) { inside = false; break; }
        if (muK === 0n) nz++;
      }
      if (inside) {
        cc.has_Cw_pos = true;
        if (nz >= 2) tags.push('Cw0');
        else if (nz === 1) tags.push('Cw1');
        else tags.push('Cw');
      }
    }
  }

  // D00 detection (exact integer, matching C++ int64_t)
  let hasD00 = false;
  for (let i = 0; i < 3; i++) {
    const det = BigInt(V[i][0]) * BigInt(W[i][1]) - BigInt(V[i][1]) * BigInt(W[i][0]);
    if (det === 0n) hasD00 = true;
  }
  for (const p of punctures)
    if (p.isEdge || p.isVertex) hasD00 = true;
  if (hasD00) tags.push('D00');

  // D11 detection
  let d11Edge = -1;
  for (let k = 0; k < 3; k++) {
    if (P_bigint[k][0] === 0n && P_bigint[k][1] === 0n && P_bigint[k][2] === 0n) {
      d11Edge = k; break;
    }
  }
  if (d11Edge < 0) {
    const te2d = [[1,2],[0,2],[0,1]];
    for (let e = 0; e < 3; e++) {
      const a = te2d[e][0], b = te2d[e][1];
      const detA = V[a][0] * W[a][1] - V[a][1] * W[a][0];
      const detB = V[b][0] * W[b][1] - V[b][1] * W[b][0];
      if (detA !== 0 || detB !== 0) continue;
      // Both PV, check lambda compat
      let lamCompatible = false;
      const vAzero = V[a][0] === 0 && V[a][1] === 0;
      const wAzero = W[a][0] === 0 && W[a][1] === 0;
      const vBzero = V[b][0] === 0 && V[b][1] === 0;
      const wBzero = W[b][0] === 0 && W[b][1] === 0;
      if ((vAzero && wAzero) || (vBzero && wBzero)) {
        lamCompatible = true;
      } else {
        let lamNumA = 0, lamDenA = 1, lamNumB = 0, lamDenB = 1;
        if (vAzero) { lamNumA = 0; lamDenA = 1; }
        else if (wAzero) { lamNumA = 1; lamDenA = 0; }
        else { for (let c = 0; c < 2; c++) if (W[a][c] !== 0) { lamNumA = -V[a][c]; lamDenA = W[a][c]; break; } }
        if (vBzero) { lamNumB = 0; lamDenB = 1; }
        else if (wBzero) { lamNumB = 1; lamDenB = 0; }
        else { for (let c = 0; c < 2; c++) if (W[b][c] !== 0) { lamNumB = -V[b][c]; lamDenB = W[b][c]; break; } }
        lamCompatible = lamNumA * lamDenB === lamNumB * lamDenA;
      }
      if (lamCompatible) { d11Edge = e; break; }
    }
  }
  if (d11Edge >= 0) tags.push('D11');

  // D22 detection
  {
    const allCompat = (() => {
      for (let i = 0; i < 3; i++) {
        const det = BigInt(V[i][0]) * BigInt(W[i][1]) - BigInt(V[i][1]) * BigInt(W[i][0]);
        if (det !== 0n) return false;
      }
      // Check all pairs have same lambda
      function getLam(i) {
        const vz = V[i][0] === 0 && V[i][1] === 0;
        const wz = W[i][0] === 0 && W[i][1] === 0;
        if (vz && wz) return { num: 0, den: 0, any: true };
        if (vz) return { num: 0, den: 1, any: false };
        if (wz) return { num: 1, den: 0, any: false };
        for (let c = 0; c < 2; c++)
          if (W[i][c] !== 0) return { num: -V[i][c], den: W[i][c], any: false };
        return { num: 0, den: 1, any: false };
      }
      const l0 = getLam(0), l1 = getLam(1), l2 = getLam(2);
      function compat(a, b) {
        if (a.any || b.any) return true;
        return a.num * b.den === b.num * a.den;
      }
      return compat(l0, l1) && compat(l1, l2);
    })();
    let d22 = allCompat;
    if (!d22 && qType === 'Qz') {
      let allPZero = true;
      for (let k = 0; k < 3; k++)
        if (P_bigint[k][0] !== 0n || P_bigint[k][1] !== 0n || P_bigint[k][2] !== 0n)
          allPZero = false;
      d22 = allPZero;
    }
    if (d22) tags.push('D22');
  }

  // Bubble detection
  if (pv2.punctures.length === 0 && cc.n_Q_roots === 0 && degQ > 0) {
    let bounded = true;
    for (let k = 0; k < 3; k++)
      if (effectiveDegree(P_bigint[k], 2) > degQ) { bounded = false; break; }
    if (bounded) {
      let atZero = true, atInf = true;
      for (let k = 0; k < 3; k++) {
        if (P_bigint[k][0] * Q_bigint[0] <= 0n) atZero = false;
        const pkLead = (effectiveDegree(P_bigint[k], 2) === degQ) ? P_bigint[k][degQ] : 0n;
        if (pkLead * Q_bigint[degQ] <= 0n) atInf = false;
      }
      if (atZero && atInf) {
        cc.has_B = true;
        tags.push('B');
      }
    }
  }

  // TN detection
  let nTN = 0;
  for (let k = 0; k < 3; k++) {
    const degPk = effectiveDegree(P_bigint[k], 2);
    if (degPk < 2) continue;
    const a = P_bigint[k][2], b = P_bigint[k][1];
    const discPk = b * b - 4n * P_bigint[k][0] * a;
    if (discPk !== 0n) continue;
    const a2_4 = 4n * a * a, ab_2 = 2n * a * b, b2 = b * b;
    const q4a2 = a2_4 * Q_bigint[0] - ab_2 * Q_bigint[1] + b2 * Q_bigint[2];
    if (q4a2 === 0n) continue;
    let inside = true, nz = 0;
    for (let j = 0; j < 3; j++) {
      if (j === k) continue;
      const muJ = a2_4 * P_bigint[j][0] - ab_2 * P_bigint[j][1] + b2 * P_bigint[j][2];
      if (muJ * q4a2 < 0n) { inside = false; break; }
      if (muJ === 0n) nz++;
    }
    if (!inside) continue;
    if (nz + 1 >= 2) continue;
    nTN++;
  }
  if (nTN > 0) tags.push('TN');

  if (tags.length > 0) cc.category += '_' + tags.join('_');

  // Phase 8: Pairing
  cc.pairs = [];
  for (const pair of pv2.pairs) {
    const pi_a = pair.a, pi_b = pair.b;
    // C++: is_cross = merge_infinity && (qa != qb) — uses solver's Q_red intervals
    const qa = pv2.punctures[pi_a] ? pv2.punctures[pi_a].qInterval : -1;
    const qb = pv2.punctures[pi_b] ? pv2.punctures[pi_b].qInterval : -1;
    const is_cross = pv2.mergeInfinity && (qa !== qb);
    // C++: interval_idx = punctures[a].interval_idx (h_below-adjusted)
    const iv_idx = punctures[pi_a] ? punctures[pi_a].intervalIdx : -1;
    cc.pairs.push({ pi_a, pi_b, is_cross, contains_infinity: pair.contains_infinity,
                    interval_idx: iv_idx });
  }

  return cc;
}

// Inline polyEval helper to avoid circular dependency
function await_polyEval() {
  return {
    polyEval(coeffs, x) {
      let val = 0;
      for (let i = coeffs.length - 1; i >= 0; i--)
        val = val * x + coeffs[i];
      return val;
    }
  };
}
