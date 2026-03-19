/**
 * ExactPV2 2D Triangle Solver.
 * Ported from ftk2 parallel_vector_solver.hpp solve_pv_tri_2d().
 * All topological decisions use exact BigInt arithmetic.
 */
import {
  bigSign,
  effectiveDegree,
  polyGcdFull, polyExactDiv,
  resultantSign, signsAtRoots, compareRoots,
  pairPuncturesRP1,
} from './exactpv2.js';

// ── 2D-specific QP computation ──

/**
 * Compute Q[3] and P[3][3] from V[3][2], W[3][2] (all BigInt).
 * Uses det(A+λB) convention so P[0]+P[1]+P[2] = Q.
 */
export function computeTriQP2D(V, W) {
  const a00 = V[0][0] - V[2][0], a01 = V[1][0] - V[2][0];
  const a10 = V[0][1] - V[2][1], a11 = V[1][1] - V[2][1];
  const b00 = W[0][0] - W[2][0], b01 = W[1][0] - W[2][0];
  const b10 = W[0][1] - W[2][1], b11 = W[1][1] - W[2][1];

  const Q = [0n, 0n, 0n];
  Q[0] = a00 * a11 - a10 * a01;
  Q[2] = b00 * b11 - b10 * b01;
  Q[1] = a00 * b11 + b00 * a11 - a01 * b10 - b01 * a10;

  const P = [[0n, 0n, 0n], [0n, 0n, 0n], [0n, 0n, 0n]];
  // P[0] = det(U₁, U₂)
  P[0][0] = V[1][0] * V[2][1] - V[1][1] * V[2][0];
  P[0][2] = W[1][0] * W[2][1] - W[1][1] * W[2][0];
  P[0][1] = V[1][0] * W[2][1] + W[1][0] * V[2][1]
          - V[1][1] * W[2][0] - W[1][1] * V[2][0];
  // P[1] = det(U₂, U₀)
  P[1][0] = V[2][0] * V[0][1] - V[2][1] * V[0][0];
  P[1][2] = W[2][0] * W[0][1] - W[2][1] * W[0][0];
  P[1][1] = V[2][0] * W[0][1] + W[2][0] * V[0][1]
          - V[2][1] * W[0][0] - W[2][1] * V[0][0];
  // P[2] = det(U₀, U₁)
  P[2][0] = V[0][0] * V[1][1] - V[0][1] * V[1][0];
  P[2][2] = W[0][0] * W[1][1] - W[0][1] * W[1][0];
  P[2][1] = V[0][0] * W[1][1] + W[0][0] * V[1][1]
          - V[0][1] * W[1][0] - W[0][1] * V[1][0];

  return { Q, P };
}

// ── 2D-specific sign evaluation ──

/** Exact sign of A + Bs*sqrt(D), all BigInt, D >= 0. */
export function signABsSqrtD(A, Bs, D) {
  if (Bs === 0n) return bigSign(A);
  if (D <= 0n) return bigSign(A);
  if (Bs > 0n) {
    if (A >= 0n) return 1;
    return (Bs * Bs * D > A * A) ? 1 : (Bs * Bs * D < A * A) ? -1 : 0;
  } else {
    if (A <= 0n) return -1;
    return (A * A > Bs * Bs * D) ? 1 : (A * A < Bs * Bs * D) ? -1 : 0;
  }
}

/**
 * Sign of polynomial pj (degree dj) at root ri of polynomial pk (degree dk).
 * Both degree <= 2. Pure integer arithmetic.
 */
export function evalSignAtRoot2D(pj, dj, pk, dk, ri) {
  if (dj === 0) return bigSign(pj[0]);
  if (dk === 1) {
    const a = pk[0], b = pk[1];
    if (dj === 1) {
      const val = pj[0] * b - pj[1] * a;
      const prod = val * b;
      return bigSign(prod);
    }
    // dj === 2
    const val = pj[0] * b * b - pj[1] * a * b + pj[2] * a * a;
    return bigSign(val);
  }
  if (dk === 2) {
    const a = pk[0], b = pk[1], c = pk[2];
    const D = b * b - 4n * a * c;
    const s = (ri === 0) ? -1n : 1n;
    if (dj === 1) {
      const A = 2n * c * pj[0] - b * pj[1];
      const Bs = pj[1] * s;
      const sc = (c > 0n) ? 1 : -1;
      return signABsSqrtD(A, Bs, D) * sc;
    }
    // dj === 2
    const d = pj[0], e = pj[1], f = pj[2];
    const A = 4n * c * c * d - 2n * b * c * e + 2n * f * b * b - 4n * a * c * f;
    const B = 2n * (c * e - b * f);
    return signABsSqrtD(A, B * s, D);
  }
  return 0;
}

// ── Main 2D Solver ──

export function solvePVTri2D(Q_raw, P_raw) {
  const result = {
    punctures: [],
    pairs: [],
    hasPassthrough: false,
    passthroughDeg: 0,
    mergeInfinity: false,
    nQrRoots: 0,
    h: [0n, 0n, 0n],
    hDeg: 0,
    hNRoots: 0,
    P_red: [[0n,0n,0n],[0n,0n,0n],[0n,0n,0n]],
    degP_red: [0, 0, 0],
    nDistinctRed: [0, 0, 0],
  };

  // Step 0: Copy and determine effective degrees
  const P = P_raw.map(pk => pk.slice());
  const Q = Q_raw.slice();
  const degP = P.map(pk => effectiveDegree(pk, 2));
  let degQ = effectiveDegree(Q, 2);

  // Step 1: Pass-through factoring
  let { h, dh } = polyGcdFull(P[0], degP[0], P[1], degP[1]);
  {
    const res = polyGcdFull(h, dh, P[2], degP[2]);
    h = res.h; dh = res.dh;
  }

  let Q_red = [0n, 0n, 0n];
  let degQ_red = degQ;
  if (dh >= 1) {
    result.hasPassthrough = true;
    result.passthroughDeg = dh;
    const { q: qDiv, dq } = polyExactDiv(Q, degQ, h, dh);
    for (let i = 0; i <= dq && i < 3; i++) Q_red[i] = qDiv[i];
    degQ_red = dq;
    for (let k = 0; k < 3; k++) {
      const { q: pDiv, dq: dp } = polyExactDiv(P[k], degP[k], h, dh);
      for (let i = 0; i < 3; i++) P[k][i] = (i <= dp) ? pDiv[i] : 0n;
      degP[k] = dp;
    }
  } else {
    for (let i = 0; i <= degQ && i < 3; i++) Q_red[i] = Q[i];
  }

  // Step 2: Per-face root validity
  const allRoots = [];
  const nDistinct = [0, 0, 0];

  for (let k = 0; k < 3; k++) {
    if (degP[k] <= 0) continue;
    if (degP[k] === 1) {
      nDistinct[k] = 1;
    } else {
      const d2 = P[k][1] * P[k][1] - 4n * P[k][2] * P[k][0];
      if (d2 > 0n) nDistinct[k] = 2;
      else if (d2 === 0n) nDistinct[k] = 1;
      else nDistinct[k] = 0;
    }

    for (let ri = 0; ri < nDistinct[k]; ri++) {
      let valid = true;
      let firstSign = 0;
      for (let j = 0; j < 3; j++) {
        if (j === k) continue;
        const sJ = signsAtRoots(P[k], degP[k], P[j], degP[j], 2, 0);
        const s = (ri < sJ.length) ? sJ[ri] : 0;
        if (s === 0) { valid = false; break; }
        if (firstSign === 0) firstSign = s;
        else if (s !== firstSign) { valid = false; break; }
      }

      if (!valid) {
        let nZero = 0, redoValid = true, redoFirst = 0;
        for (let j = 0; j < 3; j++) {
          if (j === k) continue;
          const sJ = signsAtRoots(P[k], degP[k], P[j], degP[j], 2, 0);
          const s = (ri < sJ.length) ? sJ[ri] : 0;
          if (s === 0) nZero++;
          else {
            if (redoFirst === 0) redoFirst = s;
            else if (s !== redoFirst) redoValid = false;
          }
        }
        if (nZero >= 1 && redoValid) valid = true;
      }

      if (!valid) continue;
      allRoots.push({
        face: k, rootIdx: ri, qInterval: -1, valid: true,
        isEdge: false, isVertex: false, isInfinity: false,
        edgeFaces: [-1, -1],
      });
    }
  }

  // Step 2b: Infinity punctures
  {
    const dQ = degQ_red;
    const Q_lead = (dQ >= 0 && dQ <= 2) ? Q_red[dQ] : 0n;
    if (dQ >= 1 && Q_lead !== 0n) {
      let allPBounded = true;
      for (let kb = 0; kb < 3; kb++)
        if (degP[kb] > dQ) { allPBounded = false; break; }

      if (allPBounded) {
        for (let k = 0; k < 3; k++) {
          if (degP[k] >= dQ) continue;

          let firstSign3 = 0, validInf = true, nZero3 = 0;
          for (let j = 0; j < 3; j++) {
            if (j === k) continue;
            const pjLead = P[j][dQ];
            if (pjLead === 0n) { nZero3++; }
            else {
              const s = pjLead > 0n ? 1 : -1;
              if (firstSign3 === 0) firstSign3 = s;
              else if (s !== firstSign3) { validInf = false; break; }
            }
          }
          if (!validInf || firstSign3 === 0) continue;

          if (nZero3 >= 1) {
            // Cw1 edge at infinity
            let jZero = -1;
            for (let j = 0; j < 3; j++) {
              if (j === k) continue;
              if (P[j][dQ] === 0n) { jZero = j; break; }
            }
            if (jZero < 0) continue;
            if (k > jZero) continue; // dedup

            const pkSub = (dQ >= 1) ? P[k][dQ - 1] : 0n;
            const pjSub = (dQ >= 1) ? P[jZero][dQ - 1] : 0n;
            if (pkSub * pjSub < 0n) continue;
            if (pkSub === 0n && pjSub === 0n) continue;
            if (pkSub === 0n && P[k][0] * Q_lead <= 0n) continue;
            if (pjSub === 0n && P[jZero][0] * Q_lead <= 0n) continue;

            allRoots.push({
              face: k, rootIdx: -1, qInterval: -1, valid: true,
              isEdge: true, isVertex: false, isInfinity: true,
              edgeFaces: [Math.min(k, jZero), Math.max(k, jZero)],
            });
            continue;
          }

          // Face-interior at infinity
          let dK = degP[k];
          while (dK > 0 && P[k][dK] === 0n) dK--;
          const gap = dQ - dK;
          let nInfPunc = 1;
          if (gap >= 2 && gap % 2 === 0) {
            const signProd = P[k][dK] * Q_lead;
            if (signProd < 0n) nInfPunc = 0;
            else if (signProd > 0n) nInfPunc = 2;
          }

          for (let ip = 0; ip < nInfPunc; ip++) {
            allRoots.push({
              face: k, rootIdx: -1, qInterval: -1, valid: true,
              isEdge: false, isVertex: false, isInfinity: true,
              edgeFaces: [-1, -1],
            });
          }
        }
      }
    }
  }

  // Step 3: Edge detection
  for (let i = 0; i < 3; i++) {
    for (let j = i + 1; j < 3; j++) {
      if (degP[i] === 0 || degP[j] === 0) continue;
      const res = resultantSign(P[i], degP[i], P[j], degP[j]);
      if (res === 0) {
        const owner = i;
        for (const ri of allRoots) {
          if (ri.isInfinity) continue;
          if (ri.face === i || ri.face === j) {
            if (ri.face === i) {
              const sJ = signsAtRoots(P[i], degP[i], P[j], degP[j], 2, 0);
              if (ri.rootIdx < sJ.length && sJ[ri.rootIdx] === 0) {
                ri.isEdge = true;
                ri.edgeFaces = [i, j];
                if (ri.face !== owner) ri.valid = false;
              }
            } else {
              const sI = signsAtRoots(P[j], degP[j], P[i], degP[i], 2, 0);
              if (ri.rootIdx < sI.length && sI[ri.rootIdx] === 0) {
                ri.isEdge = true;
                ri.edgeFaces = [i, j];
                if (ri.face !== owner) ri.valid = false;
              }
            }
          }
        }
      }
    }
  }

  // Step 4: Vertex detection (3-way shared root)
  if (degP[0] > 0 && degP[1] > 0 && degP[2] > 0) {
    const r01 = resultantSign(P[0], degP[0], P[1], degP[1]);
    if (r01 === 0) {
      const r12 = resultantSign(P[1], degP[1], P[2], degP[2]);
      if (r12 === 0) {
        const r02 = resultantSign(P[0], degP[0], P[2], degP[2]);
        if (r02 === 0) {
          for (const ri of allRoots) {
            if (ri.isInfinity) continue;
            let isShared = true;
            for (let j = 0; j < 3; j++) {
              if (j === ri.face) continue;
              const sO = signsAtRoots(P[ri.face], degP[ri.face], P[j], degP[j], 2, 0);
              if (ri.rootIdx < sO.length && sO[ri.rootIdx] !== 0)
                isShared = false;
            }
            if (isShared) {
              ri.isVertex = true;
              ri.isEdge = false;
              if (ri.face !== 0) ri.valid = false;
            }
          }
        }
      }
    }
  }

  // Step 4c: Edge/vertex pass-through exclusion
  const dP = [];
  for (let k = 0; k < 3; k++) {
    dP.push([P[k][1], 2n * P[k][2]]);
  }

  for (const ri of allRoots) {
    if (!ri.valid || ri.isInfinity) continue;

    if (ri.isEdge) {
      const k1 = ri.edgeFaces[0], k2 = ri.edgeFaces[1];
      if (k1 < 0 || k2 < 0) continue;
      const prod = [0n, 0n, 0n];
      for (let i = 0; i <= 1; i++)
        for (let j = 0; j <= 1; j++)
          prod[i + j] += dP[k1][i] * dP[k2][j];
      let dp = 2;
      while (dp > 0 && prod[dp] === 0n) dp--;

      const sProd = signsAtRoots(P[ri.face], degP[ri.face], prod, dp, 2, 0);
      if (ri.rootIdx < sProd.length && sProd[ri.rootIdx] < 0) {
        ri.valid = false;
      } else if (ri.rootIdx < sProd.length && sProd[ri.rootIdx] === 0) {
        for (let ki = 0; ki < 2; ki++) {
          const kf = ri.edgeFaces[ki];
          const degDk = effectiveDegree(dP[kf], 1);
          if (degDk < 0) continue;
          const sDk = signsAtRoots(P[ri.face], degP[ri.face], dP[kf], degDk, 2, 0);
          if (sDk[ri.rootIdx] !== 0) continue;
          // P''[kf] = 2*P[kf][2] (constant)
          const pp = 2n * P[kf][2];
          const pqProd = new Array(3).fill(0n);
          for (let jj = 0; jj <= degQ_red; jj++)
            pqProd[jj] += pp * Q_red[jj];
          let dpq = degQ_red;
          while (dpq > 0 && pqProd[dpq] === 0n) dpq--;
          const sPQ = signsAtRoots(P[ri.face], degP[ri.face], pqProd, dpq, 2, 0);
          if (sPQ[ri.rootIdx] < 0) { ri.valid = false; break; }
        }
      }
    } else if (ri.isVertex) {
      const f = ri.face;
      const otherFaces = [];
      for (let j = 0; j < 3; j++) {
        if (j === f || otherFaces.length >= 2) continue;
        if (degP[j] === 0) continue;
        const sJ = signsAtRoots(P[f], degP[f], P[j], degP[j], 2, 0);
        if (ri.rootIdx < sJ.length && sJ[ri.rootIdx] === 0)
          otherFaces.push(j);
      }
      if (otherFaces.length === 2) {
        const k1 = otherFaces[0], k2 = otherFaces[1];
        const prod1 = [0n, 0n, 0n];
        for (let i = 0; i <= 1; i++)
          for (let j = 0; j <= 1; j++)
            prod1[i + j] += dP[f][i] * dP[k1][j];
        let dp1 = 2;
        while (dp1 > 0 && prod1[dp1] === 0n) dp1--;
        const prod2 = [0n, 0n, 0n];
        for (let i = 0; i <= 1; i++)
          for (let j = 0; j <= 1; j++)
            prod2[i + j] += dP[f][i] * dP[k2][j];
        let dp2 = 2;
        while (dp2 > 0 && prod2[dp2] === 0n) dp2--;
        const s1 = signsAtRoots(P[f], degP[f], prod1, dp1, 2, 0);
        const s2 = signsAtRoots(P[f], degP[f], prod2, dp2, 2, 0);
        if ((ri.rootIdx < s1.length && s1[ri.rootIdx] < 0) ||
            (ri.rootIdx < s2.length && s2[ri.rootIdx] < 0))
          ri.valid = false;
      }
    }
  }

  // Step 4d: TN face-interior handling
  {
    const nAllOrig = allRoots.length;
    for (let r = 0; r < nAllOrig; r++) {
      const ri = allRoots[r];
      if (!ri.valid || ri.isInfinity || ri.isEdge || ri.isVertex) continue;
      if (ri.rootIdx >= 2 || degP[ri.face] < 2) continue;

      const face = ri.face;
      const degDp = effectiveDegree(dP[face], 1);
      if (degDp === 0 && dP[face][0] === 0n) continue;
      const sDp = signsAtRoots(P[face], degP[face], dP[face], degDp, 2, 0);
      if (sDp[ri.rootIdx] !== 0) continue;

      const Ppp = 2n * P[face][2];
      if (Ppp === 0n) {
        ri.valid = false;
        result.hasPassthrough = true;
        continue;
      }

      const sQ = signsAtRoots(P[face], degP[face], Q, degQ, 2, 0);
      const ppSign = Ppp > 0n ? 1 : -1;
      const qSign = sQ[ri.rootIdx] || 0;

      if (ppSign * qSign < 0) {
        ri.valid = false;
        result.hasPassthrough = true;
      } else if (ppSign * qSign > 0) {
        allRoots.push({ ...ri });
      }
    }
  }

  // Step 5: Collect valid punctures
  const validRoots = allRoots.filter(ri => ri.valid);
  if (validRoots.length === 0) {
    result.h = h; result.hDeg = dh;
    for (let k = 0; k < 3; k++) {
      result.P_red[k] = P[k].slice();
      result.degP_red[k] = degP[k];
      result.nDistinctRed[k] = nDistinct[k];
    }
    return result;
  }

  // Step 6: Sort by λ
  validRoots.sort((a, b) => {
    if (a.isInfinity && b.isInfinity) return 0;
    if (a.isInfinity) return 1;
    if (b.isInfinity) return -1;
    if (a.face === b.face)
      return a.rootIdx - b.rootIdx;
    return compareRoots(P[a.face], degP[a.face], nDistinct[a.face], a.rootIdx,
                        P[b.face], degP[b.face], nDistinct[b.face], b.rootIdx);
  });

  // Step 7: Q_red-interval assignment
  degQ_red = effectiveDegree(Q_red, degQ_red);
  let nQrRoots = 0;
  if (degQ_red >= 1) {
    if (degQ_red === 1) nQrRoots = 1;
    else {
      const d2 = Q_red[1] * Q_red[1] - 4n * Q_red[2] * Q_red[0];
      nQrRoots = (d2 > 0n) ? 2 : (d2 === 0n) ? 1 : 0;
    }
  }

  for (const ri of validRoots) {
    if (ri.isInfinity) { ri.qInterval = nQrRoots; continue; }
    if (nQrRoots === 0 || degQ_red === 0) { ri.qInterval = 0; continue; }
    let countBelow = 0;
    for (let qi = 0; qi < nQrRoots; qi++) {
      const cmp = compareRoots(P[ri.face], degP[ri.face],
        nDistinct[ri.face], ri.rootIdx,
        Q_red, degQ_red, nQrRoots, qi);
      if (cmp > 0) countBelow++;
    }
    ri.qInterval = countBelow;
  }

  // Step 8: Fill result
  for (const ri of validRoots) {
    result.punctures.push({
      face: ri.face, rootIdx: ri.rootIdx, qInterval: ri.qInterval,
      isEdge: ri.isEdge, isVertex: ri.isVertex, isInfinity: ri.isInfinity,
      edgeFaces: ri.edgeFaces,
    });
  }

  // Step 9: RP1 pairing
  let mergeInfinity = false;
  if (Q[2] === 0n) {
    mergeInfinity = true;
  } else {
    mergeInfinity = true;
    for (let k = 0; k < 3; k++) {
      if ((P[k][2] > 0n && Q[2] < 0n) || (P[k][2] < 0n && Q[2] > 0n)) {
        mergeInfinity = false;
        break;
      }
    }
  }

  const nPunctures = result.punctures.length;
  if (nPunctures >= 2) {
    // Pad P[3][3] → P_pad[3][4] for pair_punctures_rp1
    const P_pad = P.map(pk => [pk[0], pk[1], pk[2], 0n]);

    const sorted = [], pFaceArr = [], pRidxArr = [], pQiArr = [];
    let nFin = 0;
    for (let i = 0; i < nPunctures; i++) {
      sorted.push(i);
      pFaceArr.push(result.punctures[i].face);
      pRidxArr.push(result.punctures[i].rootIdx);
      pQiArr.push(result.punctures[i].qInterval);
      if (result.punctures[i].rootIdx >= 0) nFin++;
    }

    const Q_red4 = [Q_red[0], Q_red[1], Q_red[2], 0n];
    const degP_pad = degP.slice();
    const rp1Pairs = pairPuncturesRP1(3, sorted, nPunctures, nFin,
      pFaceArr, pRidxArr, pQiArr, nQrRoots, mergeInfinity,
      P_pad, degP_pad, nDistinct, Q_red4, degQ_red);

    for (const pair of rp1Pairs) {
      result.pairs.push({ a: pair.a, b: pair.b, contains_infinity: pair.contains_infinity });
    }
  }

  // Fill infrastructure
  result.mergeInfinity = mergeInfinity;
  result.nQrRoots = nQrRoots;
  result.h = h.slice(0, 3);
  result.hDeg = dh;
  if (dh === 0) result.hNRoots = 0;
  else if (dh === 1) result.hNRoots = 1;
  else {
    const d2h = h[1] * h[1] - 4n * h[2] * h[0];
    result.hNRoots = (d2h > 0n) ? 2 : (d2h === 0n) ? 1 : 0;
  }
  for (let k = 0; k < 3; k++) {
    result.P_red[k] = P[k].slice();
    result.degP_red[k] = degP[k];
    result.nDistinctRed[k] = nDistinct[k];
  }

  return result;
}
