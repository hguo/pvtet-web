/**
 * ExactPV2: Pure-Integer PV Solver using BigInt.
 * Ported from ftk2 parallel_vector_solver.hpp solve_pv_tet_v2().
 *
 * All topological decisions use exact BigInt arithmetic (PRS, Sylvester
 * resultants, polynomial GCD). No floats for topology. JavaScript BigInt
 * = arbitrary precision, so no overflow concerns unlike C++ __int128.
 */

// ── Helpers ──

export function bigAbs(x) { return x < 0n ? -x : x; }

export function bigGcd(a, b) {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b !== 0n) { const t = b; b = a % b; a = t; }
  return a;
}

export function bigSign(x) { return x > 0n ? 1 : x < 0n ? -1 : 0; }

// ── 1.1 Content Reduction ──

export function contentReduce(poly, deg) {
  if (deg < 0) return;
  let g = 0n;
  for (let i = 0; i <= deg; i++) {
    const v = bigAbs(poly[i]);
    if (g === 0n) g = v;
    else g = bigGcd(g, v);
  }
  if (g > 1n) for (let i = 0; i <= deg; i++) poly[i] /= g;
}

export function effectiveDegree(poly, maxDeg) {
  let d = maxDeg;
  while (d > 0 && poly[d] === 0n) d--;
  return d;
}

// ── 1.2 Pseudo-Remainder ──

export function prem(f, df, g, dg) {
  const r = new Array(8).fill(0n);
  for (let i = 0; i <= df; i++) r[i] = f[i];

  let dr = df;
  while (dr >= dg && dr >= 0) {
    if (r[dr] === 0n) { dr--; continue; }
    const rdr = r[dr];
    const gdq = g[dg];
    const shift = dr - dg;
    for (let i = 0; i <= dr; i++) r[i] *= gdq;
    for (let i = 0; i <= dg; i++) r[shift + i] -= rdr * g[i];
    dr--;
  }
  if (dr < 0) { dr = 0; r[0] = 0n; }
  while (dr > 0 && r[dr] === 0n) dr--;

  contentReduce(r, dr);
  return { r, dr };
}

// ── 1.3 Polynomial GCD ──

export function polyGcdFull(f_in, df, g_in, dg) {
  const p = new Array(8).fill(0n);
  const q = new Array(8).fill(0n);
  for (let i = 0; i <= df; i++) p[i] = f_in[i];
  for (let i = 0; i <= dg; i++) q[i] = g_in[i];

  df = effectiveDegree(p, df);
  dg = effectiveDegree(q, dg);
  contentReduce(p, df);
  contentReduce(q, dg);

  if (df === 0 && p[0] === 0n) return { h: q.slice(), dh: dg };
  if (dg === 0 && q[0] === 0n) return { h: p.slice(), dh: df };

  let dp = df, dq = dg;
  while (!(dq === 0 && q[0] === 0n)) {
    if (dp < dq) {
      for (let i = 0; i < 8; i++) { const t = p[i]; p[i] = q[i]; q[i] = t; }
      const t = dp; dp = dq; dq = t;
    }
    const { r, dr } = prem(p, dp, q, dq);
    dp = dq;
    for (let i = 0; i < 8; i++) p[i] = q[i];
    dq = dr;
    for (let i = 0; i < 8; i++) q[i] = r[i];
  }

  if (dp > 0 && p[dp] < 0n)
    for (let i = 0; i <= dp; i++) p[i] = -p[i];

  return { h: p.slice(), dh: dp };
}

// ── 1.4 Exact Polynomial Division ──

export function polyExactDiv(f, df, g, dg) {
  const rem = new Array(8).fill(0n);
  for (let i = 0; i <= df; i++) rem[i] = f[i];

  const dq = df - dg;
  const q = new Array(dq + 1).fill(0n);

  for (let i = dq; i >= 0; i--) {
    q[i] = rem[i + dg] / g[dg];
    for (let j = 0; j <= dg; j++)
      rem[i + j] -= q[i] * g[j];
  }
  return { q, dq };
}

// ── 1.5 Resultant Sign via Bareiss ──

export function resultantSign(f_in, df, g_in, dg) {
  if (df <= 0 || dg <= 0) {
    if (df === 0 && dg === 0) return 1;
    if (df === 0) {
      if (f_in[0] === 0n) return 0;
      const s = f_in[0] > 0n ? 1 : -1;
      return (dg % 2 === 0) ? 1 : s;
    }
    if (dg === 0) {
      if (g_in[0] === 0n) return 0;
      const s = g_in[0] > 0n ? 1 : -1;
      return (df % 2 === 0) ? 1 : s;
    }
    return 0;
  }

  const fc = new Array(8).fill(0n);
  const gc = new Array(8).fill(0n);
  for (let i = 0; i <= df; i++) fc[i] = f_in[i];
  for (let i = 0; i <= dg; i++) gc[i] = g_in[i];
  contentReduce(fc, df);
  contentReduce(gc, dg);

  const N = df + dg;
  // BigInt has arbitrary precision, no overflow guard needed

  // Build Sylvester matrix
  const M = Array.from({ length: N }, () => new Array(N).fill(0n));
  for (let i = 0; i < dg; i++)
    for (let j = 0; j <= df; j++)
      M[i][i + df - j] = fc[j];
  for (let i = 0; i < df; i++)
    for (let j = 0; j <= dg; j++)
      M[dg + i][i + dg - j] = gc[j];

  // Bareiss fraction-free elimination
  let prevPivot = 1n;
  let signSwaps = 0;
  for (let col = 0; col < N; col++) {
    let pivot = -1;
    for (let row = col; row < N; row++)
      if (M[row][col] !== 0n) { pivot = row; break; }
    if (pivot < 0) return 0;
    if (pivot !== col) {
      for (let j = 0; j < N; j++)
        [M[col][j], M[pivot][j]] = [M[pivot][j], M[col][j]];
      signSwaps++;
    }
    for (let row = col + 1; row < N; row++) {
      for (let j = col + 1; j < N; j++)
        M[row][j] = (M[col][col] * M[row][j] - M[row][col] * M[col][j]) / prevPivot;
      M[row][col] = 0n;
    }
    prevPivot = M[col][col];
  }

  const det = M[N - 1][N - 1];
  if (det === 0n) return 0;
  let s = det > 0n ? 1 : -1;
  if (signSwaps % 2) s = -s;
  return s;
}

// ── 1.6 Sign at Unique Root ──

export function signAtUniqueRoot(f, df, g, dg) {
  if (dg === 0) return bigSign(g[0]);
  const resSign = resultantSign(f, df, g, dg);
  if (resSign === 0) return 0;
  const lcSign = f[df] > 0n ? 1 : -1;
  const lcPowerSign = (dg % 2 === 0) ? 1 : lcSign;
  return resSign * lcPowerSign;
}

// ── Sign of f(p/q) via BigInt Horner ──

export function signPolyAtRational(f, df, p, q) {
  if (df <= 0) return bigSign(f[0]);
  if (q === 0n) {
    const lcS = bigSign(f[df]);
    const pS = bigSign(p);
    return lcS * ((df % 2 === 0) ? 1 : pS);
  }
  // val = f[0]*q^df + f[1]*p*q^(df-1) + ... + f[df]*p^df
  let val = f[0];
  let pPow = p;
  for (let i = 1; i <= df; i++) {
    val = val * q + pPow * f[i];
    if (i < df) pPow = pPow * p;
  }
  const valSign = bigSign(val);
  if (valSign === 0) return 0;
  const qSign = q > 0n ? 1 : -1;
  const qdfSign = (df % 2 === 0) ? 1 : qSign;
  return valSign * qdfSign;
}

// ── 1.7a Derivative ──

export function polyDerivative(f, df) {
  if (df <= 0) return { fp: [0n], dfp: 0 };
  const fp = new Array(df).fill(0n);
  for (let i = 0; i < df; i++) fp[i] = BigInt(i + 1) * f[i + 1];
  return { fp, dfp: df - 1 };
}

// ── 1.7b Square-free ──

export function polySqfree(f, df) {
  if (df <= 1) {
    const sf = f.slice(0, df + 1);
    return { sf, dsf: df };
  }
  const { fp, dfp } = polyDerivative(f, df);
  const { h: g, dh: dg } = polyGcdFull(f, df, fp, dfp);

  if (dg === 0) {
    const sf = f.slice(0, df + 1);
    contentReduce(sf, df);
    return { sf, dsf: df };
  }

  const { q: sf, dq: dsf } = polyExactDiv(f, df, g, dg);
  contentReduce(sf, dsf);
  return { sf, dsf };
}

// ── 1.7c Count roots below rational ──

export function countRootsBelowRational(f, df, p, q) {
  const signFt = signPolyAtRational(f, df, p, q);
  if (signFt === 0) return -1;

  const { fp, dfp } = polyDerivative(f, df);
  const signFpt = signPolyAtRational(fp, dfp, p, q);

  const signLc = f[df] > 0n ? 1 : -1;
  const eft = signFt * signLc;
  const efpt = signFpt * signLc;

  if (eft > 0) {
    if (efpt > 0) {
      const { fp: fpp, dfp: dfpp } = polyDerivative(fp, dfp);
      const efppt = signPolyAtRational(fpp, dfpp, p, q) * signLc;
      return (efppt <= 0) ? 1 : 3;
    }
    return 1;
  } else {
    if (efpt > 0) {
      const { fp: fpp, dfp: dfpp } = polyDerivative(fp, dfp);
      const efppt = signPolyAtRational(fpp, dfpp, p, q) * signLc;
      return (efppt <= 0) ? 0 : 2;
    }
    return 2;
  }
}

// ── Discriminant sign (BigInt, no overflow) ──

export function discriminantSign(P) {
  if (P[3] === 0n) return 0;
  let g = bigAbs(P[0]);
  for (let i = 1; i < 4; i++) g = bigGcd(g, bigAbs(P[i]));
  if (g === 0n) return 0;
  const a = P[0] / g, b = P[1] / g, c = P[2] / g, d = P[3] / g;
  const disc = 18n*a*b*c*d - 4n*b*b*b*d + b*b*c*c - 4n*a*c*c*c - 27n*a*a*d*d;
  return disc > 0n ? 1 : disc < 0n ? -1 : 0;
}

// ── 1.7 Signs at Roots ──

export function signsAtRoots(f_in, dfMax, g_in, dgMax, maxSigns, _depth) {
  if (_depth === undefined) _depth = 0;
  if (_depth > 20) return [];

  const f = new Array(8).fill(0n);
  const g = new Array(8).fill(0n);
  for (let i = 0; i <= dfMax && i < 8; i++) f[i] = f_in[i];
  for (let i = 0; i <= dgMax && i < 8; i++) g[i] = g_in[i];
  const df = effectiveDegree(f, Math.min(dfMax, 7));
  const dg = effectiveDegree(g, Math.min(dgMax, 7));

  contentReduce(f, df);
  contentReduce(g, dg);

  if (df === 0) return [];

  // Linear f: one root at -f[0]/f[1]
  if (df === 1) {
    if (maxSigns < 1) return [];
    return [signPolyAtRational(g, dg, -f[0], f[1])];
  }

  // Quadratic f
  if (df === 2) {
    const disc = f[1] * f[1] - 4n * f[2] * f[0];
    if (disc < 0n) return [];
    if (disc === 0n) {
      if (maxSigns < 1) return [];
      return [signPolyAtRational(g, dg, -f[1], 2n * f[2])];
    }
    // Two distinct roots
    if (dg >= df) {
      const { r, dr } = prem(g, dg, f, df);
      const exp = dg - df + 1;
      const lcFSign = f[df] > 0n ? 1 : -1;
      const scaleSign = (exp % 2 === 0) ? 1 : lcFSign;
      const subSigns = signsAtRoots(f, df, r, dr, 2, _depth + 1);
      return subSigns.map(s => s * scaleSign);
    }
    if (dg === 0) {
      const gs = bigSign(g[0]);
      return [gs, gs];
    }
    // dg === 1
    const signFt = signPolyAtRational(f, df, -g[0], g[1]);
    const signLcF = f[df] > 0n ? 1 : -1;
    const signGLc = g[dg] > 0n ? 1 : -1;

    if (signFt === 0) {
      const { fp } = polyDerivative(f, df);
      const signFpt = signPolyAtRational(fp, df - 1, -g[0], g[1]);
      if (signFpt * signLcF < 0) {
        return [0, signGLc];
      } else {
        return [-signGLc, 0];
      }
    }

    if (signFt * signLcF > 0) {
      const { fp } = polyDerivative(f, df);
      const signFpt = signPolyAtRational(fp, df - 1, -g[0], g[1]);
      if (signFpt * signLcF < 0) {
        return [signGLc, signGLc];
      } else {
        return [-signGLc, -signGLc];
      }
    } else {
      return [-signGLc, signGLc];
    }
  }

  // Cubic f
  const discSign = discriminantSign(f);

  if (discSign < 0) {
    // One real root
    if (maxSigns < 1) return [];
    let scale1 = 1;
    let r1, dr1;
    if (dg >= df) {
      const res = prem(g, dg, f, df);
      r1 = res.r; dr1 = res.dr;
      const exp1 = dg - df + 1;
      const lcSign1 = f[df] > 0n ? 1 : -1;
      scale1 = (exp1 % 2 === 0) ? 1 : lcSign1;
    } else {
      r1 = new Array(8).fill(0n);
      for (let j = 0; j <= dg; j++) r1[j] = g[j];
      dr1 = dg;
      contentReduce(r1, dr1);
    }
    dr1 = effectiveDegree(r1, dr1);

    if (dr1 === 0) {
      return [bigSign(r1[0]) * scale1];
    }
    if (dr1 === 1) {
      const sft = signPolyAtRational(f, df, -r1[0], r1[1]);
      const signLcF = f[df] > 0n ? 1 : -1;
      const r1LcSign = r1[dr1] > 0n ? 1 : -1;
      if (sft === 0) return [0];
      if (sft * signLcF > 0)
        return [-r1LcSign * scale1];
      else
        return [r1LcSign * scale1];
    }
    // dr1 === 2
    {
      const discR1 = r1[1] * r1[1] - 4n * r1[2] * r1[0];
      const signR12 = r1[2] > 0n ? 1 : -1;
      if (discR1 < 0n) return [signR12 * scale1];
      if (discR1 === 0n) {
        const sft = signPolyAtRational(f, df, -r1[1], 2n * r1[2]);
        return [(sft === 0) ? 0 : signR12 * scale1];
      }
      const fAtRho = signsAtRoots(r1, dr1, f, df, 2, _depth + 1);
      const signLcF = f[df] > 0n ? 1 : -1;
      const pos1 = (fAtRho[0] || 0) * signLcF;
      const pos2 = (fAtRho[1] || 0) * signLcF;
      if (fAtRho[0] === 0 || fAtRho[1] === 0) {
        return [0];
      } else if (pos1 > 0) {
        return [signR12 * scale1];
      } else if (pos2 > 0) {
        return [-signR12 * scale1];
      } else {
        return [signR12 * scale1];
      }
    }
  }

  if (discSign === 0) {
    const { sf, dsf } = polySqfree(f, df);
    return signsAtRoots(sf, dsf, g, dg, maxSigns, _depth + 1);
  }

  // disc > 0: three distinct real roots
  if (maxSigns < 3) return [];

  let r, dr;
  let scaleSign = 1;
  if (dg >= df) {
    const res = prem(g, dg, f, df);
    r = res.r; dr = res.dr;
    const exp = dg - df + 1;
    const lcFSign = f[df] > 0n ? 1 : -1;
    scaleSign = (exp % 2 === 0) ? 1 : lcFSign;
  } else {
    r = new Array(8).fill(0n);
    for (let i = 0; i <= dg; i++) r[i] = g[i];
    dr = dg;
    contentReduce(r, dr);
  }

  // Case A: constant
  if (dr === 0) {
    const rs = bigSign(r[0]) * scaleSign;
    return [rs, rs, rs];
  }

  // Case B: linear
  if (dr === 1) {
    const nBelow = countRootsBelowRational(f, df, -r[0], r[1]);
    const signRLc = r[1] > 0n ? 1 : -1;

    if (nBelow < 0) {
      // Shared root
      const { fp } = polyDerivative(f, df);
      const signFpt = signPolyAtRational(fp, df - 1, -r[0], r[1]);
      const { fp: fpp } = polyDerivative(fp, df - 1);
      const signFppt = signPolyAtRational(fpp, df - 2, -r[0], r[1]);
      const signLc = f[df] > 0n ? 1 : -1;
      const eft = signFpt * signLc;
      const efppt = signFppt * signLc;

      let which;
      if (eft < 0) which = 1;
      else if (efppt > 0) which = 2;
      else which = 0;

      const signs = [];
      for (let i = 0; i < 3; i++) {
        if (i === which) signs.push(0);
        else if (i < which) signs.push(-signRLc * scaleSign);
        else signs.push(signRLc * scaleSign);
      }
      return signs;
    }

    const signs = [];
    for (let i = 0; i < 3; i++) {
      if (i < nBelow) signs.push(-signRLc * scaleSign);
      else signs.push(signRLc * scaleSign);
    }
    return signs;
  }

  // Case C: quadratic
  {
    const discR = r[1] * r[1] - 4n * r[2] * r[0];
    const signR2 = r[2] > 0n ? 1 : -1;

    if (discR < 0n) {
      const s = signR2 * scaleSign;
      return [s, s, s];
    }

    if (discR === 0n) {
      const nBelow = countRootsBelowRational(f, df, -r[1], 2n * r[2]);
      const ftSign = signPolyAtRational(f, df, -r[1], 2n * r[2]);
      const signs = [signR2 * scaleSign, signR2 * scaleSign, signR2 * scaleSign];
      if (ftSign === 0 && nBelow >= 0 && nBelow < 3) signs[nBelow] = 0;
      return signs;
    }

    // Two distinct roots of r
    const fAtRho = signsAtRoots(r, dr, f, df, 2, _depth + 1);

    const { fp, dfp } = polyDerivative(f, df);
    const fpAtRho = signsAtRoots(r, dr, fp, dfp, 2, _depth + 1);

    const { fp: fpp, dfp: dfpp } = polyDerivative(fp, dfp);
    const edeg = effectiveDegree(fpp, dfpp);
    const fppAtRho = signsAtRoots(r, dr, fpp, edeg, 2, _depth + 1);

    const signLc = f[df] > 0n ? 1 : -1;

    let nb1 = 0, nb2 = 0;
    for (let j = 0; j < 2; j++) {
      const sf = (fAtRho.length > j) ? fAtRho[j] : 0;
      const sfp = (fpAtRho.length > j) ? fpAtRho[j] : 0;
      const sfpp = (fppAtRho.length > j) ? fppAtRho[j] : 0;
      const eft = sf * signLc;
      const efpt = sfp * signLc;
      const efppt = sfpp * signLc;

      let nb;
      if (sf === 0) {
        if (efpt < 0) nb = 1;
        else if (efppt > 0) nb = 2;
        else nb = 0;
      } else if (eft > 0) {
        if (efpt > 0) nb = (efppt <= 0) ? 1 : 3;
        else nb = 1;
      } else {
        if (efpt > 0) nb = (efppt <= 0) ? 0 : 2;
        else nb = 2;
      }
      if (j === 0) nb1 = nb; else nb2 = nb;
    }

    const signs = [];
    for (let i = 0; i < 3; i++) {
      if (i >= nb1 && i < nb2)
        signs.push(-signR2 * scaleSign);
      else
        signs.push(signR2 * scaleSign);
    }
    if (fAtRho[0] === 0 && nb1 < 3) signs[nb1] = 0;
    if (fAtRho.length >= 2 && fAtRho[1] === 0 && nb2 < 3) signs[nb2] = 0;

    return signs;
  }
}

// ── 1.8 Root Comparison ──

export function compareRoots(f, df, fNroots, fRootIdx, g, dg, gNroots, gRootIdx) {
  df = effectiveDegree(f, df);
  dg = effectiveDegree(g, dg);
  if (df === 0 || dg === 0) return 0;

  const gAtF = signsAtRoots(f, df, g, dg, 3, 0);
  if (fRootIdx >= gAtF.length) return 0;

  const sg = gAtF[fRootIdx];

  if (sg === 0) {
    const fAtG = signsAtRoots(g, dg, f, df, 3, 0);
    for (let j = 0; j < fAtG.length; j++) {
      if (fAtG[j] === 0) {
        if (j === gRootIdx) return 0;
        return (j < gRootIdx) ? -1 : 1;
      }
    }
    return 0;
  }

  const lcGSign = g[dg] > 0n ? 1 : -1;

  if (gNroots === 1) {
    const countBelow = (sg * lcGSign > 0) ? 1 : 0;
    return (countBelow > gRootIdx) ? 1 : -1;
  }

  // Need g' at roots of f
  const { fp: gp, dfp: dgp } = polyDerivative(g, dg);
  const gpAtF = signsAtRoots(f, df, gp, dgp, 3, 0);
  const sgp = (fRootIdx < gpAtF.length) ? gpAtF[fRootIdx] : 0;

  const egt = sg * lcGSign;
  const egpt = sgp * lcGSign;

  let countBelow = 0;

  if (gNroots === 2) {
    if (egt > 0) {
      countBelow = (egpt > 0) ? 2 : 0;
    } else {
      countBelow = 1;
    }
  } else if (gNroots === 3) {
    if (egt > 0) {
      if (egpt > 0) {
        const { fp: gpp, dfp: dgpp } = polyDerivative(gp, dgp);
        const edeg = effectiveDegree(gpp, dgpp);
        const gppAtF = signsAtRoots(f, df, gpp, edeg, 3, 0);
        const egppt = ((fRootIdx < gppAtF.length) ? gppAtF[fRootIdx] : 0) * lcGSign;
        countBelow = (egppt <= 0) ? 1 : 3;
      } else {
        countBelow = 1;
      }
    } else {
      if (egpt > 0) {
        const { fp: gpp, dfp: dgpp } = polyDerivative(gp, dgp);
        const edeg = effectiveDegree(gpp, dgpp);
        const gppAtF = signsAtRoots(f, df, gpp, edeg, 3, 0);
        const egppt = ((fRootIdx < gppAtF.length) ? gppAtF[fRootIdx] : 0) * lcGSign;
        countBelow = (egppt <= 0) ? 0 : 2;
      } else {
        countBelow = 2;
      }
    }
  }

  if (countBelow > gRootIdx) return 1;
  if (countBelow < gRootIdx) return -1;
  return -1; // countBelow === gRootIdx and no shared root → α < β
}

// ── RP1 Sign Helpers ──

/** Sign of polynomial p(λ) as λ→+∞: sign of leading coefficient. */
export function signAtPlusInf(p, deg) {
  let d = deg;
  while (d > 0 && p[d] === 0n) d--;
  if (d === 0) return bigSign(p[0]);
  return p[d] > 0n ? 1 : -1;
}

/** Sign of polynomial p(λ) at the RP1 infinity point: sign(lc)*(-1)^deg. */
export function signAtInf(p, deg) {
  let d = deg;
  while (d > 0 && p[d] === 0n) d--;
  if (d === 0) return bigSign(p[0]);
  const signLc = p[d] > 0n ? 1 : -1;
  return (d % 2 === 0) ? signLc : -signLc;
}

/** Sign of square-free polynomial just after its rootIdx-th root. */
export function signJustAfterRoot(p, deg, nDistinct, rootIdx) {
  let d = deg;
  while (d > 0 && p[d] === 0n) d--;
  if (d === 0) return bigSign(p[0]);
  const signLc = p[d] > 0n ? 1 : -1;
  const exp = nDistinct - 1 - rootIdx;
  return (exp % 2 === 0) ? signLc : -signLc;
}

/** Sign of g on interior of RP1 arc whose left endpoint is (leftFace, leftRootIdx). */
export function rp1SignOnArc(g, degG, nDistG, gFace, leftFace, leftRootIdx, precomputedSign) {
  if (leftRootIdx < 0)
    return signAtInf(g, degG);
  if (gFace === leftFace)
    return signJustAfterRoot(g, degG, nDistG, leftRootIdx);
  return precomputedSign;
}

/** Check if RP1 arc is inside simplex: all P_red[k]/Q_red >= 0. */
export function rp1ArcIsInside(nFaces, sp, sq) {
  if (sq === 0) return false;
  for (let k = 0; k < nFaces; k++)
    if (sp[k] * sq < 0) return false;
  return true;
}

/** Q-interval compatibility: same interval or merged through infinity. */
export function rp1ArcsSameInterval(qiA, qiB, nQrRoots, mergeInfinity) {
  if (qiA === qiB) return true;
  if (mergeInfinity)
    return (qiA === 0 && qiB === nQrRoots) || (qiA === nQrRoots && qiB === 0);
  return false;
}

/** Build sign tables for RP1 interval arithmetic. */
function rp1BuildSignTable(nFaces, PRed, degPRed, nDistinctRed, QRed, degQRed) {
  const signPk = Array.from({ length: nFaces }, () =>
    Array.from({ length: 4 }, () => new Array(nFaces).fill(0))
  );
  const signQr = Array.from({ length: nFaces }, () => new Array(4).fill(0));

  for (let fi = 0; fi < nFaces; fi++) {
    if (nDistinctRed[fi] === 0) continue;
    for (let k = 0; k < nFaces; k++) {
      if (k === fi) continue;
      const signs = signsAtRoots(PRed[fi], degPRed[fi], PRed[k], degPRed[k], 4, 0);
      for (let ri = 0; ri < signs.length && ri < 4; ri++)
        signPk[fi][ri][k] = signs[ri];
    }
    {
      const signs = signsAtRoots(PRed[fi], degPRed[fi], QRed, degQRed, 4, 0);
      for (let ri = 0; ri < signs.length && ri < 4; ri++)
        signQr[fi][ri] = signs[ri];
    }
  }
  return { signPk, signQr };
}

/** Group consecutive inside arcs into connected regions and pair consecutively. */
function rp1GroupInsideArcs(arcInside, nSorted, nFinite, sorted) {
  const pairs = [];
  const visited = new Array(nSorted).fill(false);

  for (let start = 0; start < nSorted; start++) {
    if (!arcInside[start] || visited[start]) continue;
    const group = [sorted[start]];
    let ci = false;
    let ai = start;
    while (arcInside[ai] && !visited[ai]) {
      visited[ai] = true;
      if (ai >= nFinite - 1) ci = true;
      const next = (ai + 1) % nSorted;
      group.push(sorted[next]);
      ai = next;
      if (ai === start) break;
    }
    for (let i = 0; i + 1 < group.length; i += 2) {
      pairs.push({ a: group[i], b: group[i + 1], contains_infinity: ci });
    }
  }
  return pairs;
}

/** Full RP1 pairing of punctures. */
export function pairPuncturesRP1(nFaces, sorted, nSorted, nFinite,
  pFace, pRootIdx, pQInterval, nQrRoots, mergeInfinity,
  PRed, degPRed, nDistinctRed, QRed, degQRed) {
  if (nSorted < 2) return [];

  const { signPk, signQr } = rp1BuildSignTable(
    nFaces, PRed, degPRed, nDistinctRed, QRed, degQRed);

  const arcInside = new Array(nSorted).fill(false);
  for (let ai = 0; ai < nSorted; ai++) {
    const pi = sorted[ai];
    const fi = pFace[pi];
    const ri = pRootIdx[pi];
    const next = (ai + 1) % nSorted;
    const pj = sorted[next];

    const sp = new Array(nFaces);
    for (let k = 0; k < nFaces; k++)
      sp[k] = rp1SignOnArc(PRed[k], degPRed[k], nDistinctRed[k], k,
                            fi, ri, signPk[fi][ri < 0 ? 0 : ri][k]);
    const sq = rp1SignOnArc(QRed, degQRed, 0, -1,
                             fi, ri, signQr[fi][ri < 0 ? 0 : ri]);

    let inside = rp1ArcIsInside(nFaces, sp, sq);

    if (inside)
      inside = rp1ArcsSameInterval(pQInterval[pi], pQInterval[pj],
                                    nQrRoots, mergeInfinity);

    if (inside && ai >= nFinite - 1 && !mergeInfinity) {
      const hasInfEndpoint = pRootIdx[pi] < 0 || pRootIdx[pj] < 0;
      if (!hasInfEndpoint) inside = false;
    }

    arcInside[ai] = inside;
  }

  return rp1GroupInsideArcs(arcInside, nSorted, nFinite, sorted);
}

// ── Main Solver ──

export function solvePVTetV2(Q_raw, P_raw) {
  const result = {
    punctures: [],
    pairs: [],
    hasPassthrough: false,
    passthroughDeg: 0,
  };

  // Step 0: Copy and determine effective degrees
  const P = P_raw.map(pk => pk.slice());
  const Q = Q_raw.slice();

  const degP = P.map(pk => effectiveDegree(pk, 3));
  let degQ = effectiveDegree(Q, 3);

  // Step 1: Pass-through factoring
  let { h, dh } = polyGcdFull(P[0], degP[0], P[1], degP[1]);
  for (let k = 2; k < 4; k++) {
    const res = polyGcdFull(h, dh, P[k], degP[k]);
    h = res.h; dh = res.dh;
  }

  let Q_red, degQ_red;
  if (dh >= 1) {
    result.hasPassthrough = true;
    result.passthroughDeg = dh;
    const { q: qDiv, dq } = polyExactDiv(Q, degQ, h, dh);
    Q_red = new Array(4).fill(0n);
    for (let i = 0; i <= dq; i++) Q_red[i] = qDiv[i];
    degQ_red = dq;
    // Replace P[k] with P[k]/h to exclude pass-through roots (matching C++)
    for (let k = 0; k < 4; k++) {
      const { q: pDiv, dq: dp } = polyExactDiv(P[k], degP[k], h, dh);
      for (let i = 0; i < 4; i++) P[k][i] = (i <= dp) ? pDiv[i] : 0n;
      degP[k] = dp;
    }
  } else {
    Q_red = Q.slice();
    degQ_red = degQ;
  }

  // Step 2: For each face, determine root count and validity
  const allRoots = [];
  const nDistinct = [0, 0, 0, 0];

  for (let k = 0; k < 4; k++) {
    if (degP[k] === 0) continue;

    if (degP[k] === 1) {
      nDistinct[k] = 1;
    } else if (degP[k] === 2) {
      const d2 = P[k][1] * P[k][1] - 4n * P[k][2] * P[k][0];
      if (d2 > 0n) nDistinct[k] = 2;
      else if (d2 === 0n) nDistinct[k] = 1;
      else nDistinct[k] = 0;
    } else {
      const disc = discriminantSign(P[k]);
      if (disc > 0) nDistinct[k] = 3;
      else if (disc < 0) nDistinct[k] = 1;
      else {
        const { sf, dsf } = polySqfree(P[k], degP[k]);
        nDistinct[k] = dsf;
      }
    }

    for (let ri = 0; ri < nDistinct[k]; ri++) {
      let valid = true;
      let firstSign = 0;

      for (let j = 0; j < 4; j++) {
        if (j === k) continue;
        const signsJ = signsAtRoots(P[k], degP[k], P[j], degP[j], 3, 0);
        const s = (ri < signsJ.length) ? signsJ[ri] : 0;
        if (s === 0) { valid = false; break; }
        if (firstSign === 0) firstSign = s;
        else if (s !== firstSign) { valid = false; break; }
      }

      if (!valid) {
        // Re-check with edge/vertex awareness
        let nZero = 0;
        let redoValid = true;
        let redoFirst = 0;
        for (let j = 0; j < 4; j++) {
          if (j === k) continue;
          const signsJ = signsAtRoots(P[k], degP[k], P[j], degP[j], 3, 0);
          const s = (ri < signsJ.length) ? signsJ[ri] : 0;
          if (s === 0) {
            nZero++;
          } else {
            if (redoFirst === 0) redoFirst = s;
            else if (s !== redoFirst) redoValid = false;
          }
        }
        if (nZero >= 1 && redoValid) valid = true;
      }

      if (!valid) continue;

      allRoots.push({
        face: k,
        rootIdx: ri,
        qInterval: -1,
        valid: true,
        isEdge: false,
        isVertex: false,
        isInfinity: false,
        edgeFaces: [-1, -1],
      });
    }
  }

  // Step 2b: Infinity punctures
  for (let k = 0; k < 4; k++) {
    if (degP[k] >= 3 || P[k][3] !== 0n) continue;
    if (Q[3] === 0n) continue;

    let firstSign3 = 0;
    let validInf = true;
    let nZero3 = 0;
    for (let j = 0; j < 4; j++) {
      if (j === k) continue;
      const pj3 = P[j][3];
      if (pj3 === 0n) {
        nZero3++;
      } else {
        const s = pj3 > 0n ? 1 : -1;
        if (firstSign3 === 0) firstSign3 = s;
        else if (s !== firstSign3) { validInf = false; break; }
      }
    }
    if (!validInf || firstSign3 === 0) continue;
    if (nZero3 >= 1) continue;

    allRoots.push({
      face: k,
      rootIdx: -1,
      qInterval: -1,
      valid: true,
      isEdge: false,
      isVertex: false,
      isInfinity: true,
      edgeFaces: [-1, -1],
    });
  }

  // Step 3: Edge detection
  for (let i = 0; i < 4; i++) {
    for (let j = i + 1; j < 4; j++) {
      if (degP[i] === 0 || degP[j] === 0) continue;
      const res = resultantSign(P[i], degP[i], P[j], degP[j]);
      if (res === 0) {
        const owner = i; // smaller face index
        for (const ri of allRoots) {
          if (ri.face === i || ri.face === j) {
            if (ri.face === i) {
              const signsJ = signsAtRoots(P[i], degP[i], P[j], degP[j], 3, 0);
              if (ri.rootIdx < signsJ.length && signsJ[ri.rootIdx] === 0) {
                ri.isEdge = true;
                ri.edgeFaces = [i, j];
                if (ri.face !== owner) ri.valid = false;
              }
            } else {
              const signsI = signsAtRoots(P[j], degP[j], P[i], degP[i], 3, 0);
              if (ri.rootIdx < signsI.length && signsI[ri.rootIdx] === 0) {
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

  // Step 4: Vertex detection
  for (let i = 0; i < 4; i++) {
    for (let j = i + 1; j < 4; j++) {
      for (let l = j + 1; l < 4; l++) {
        if (degP[i] === 0 || degP[j] === 0 || degP[l] === 0) continue;
        if (resultantSign(P[i], degP[i], P[j], degP[j]) !== 0) continue;
        if (resultantSign(P[j], degP[j], P[l], degP[l]) !== 0) continue;
        if (resultantSign(P[i], degP[i], P[l], degP[l]) !== 0) continue;

        const others = [i, j, l];
        for (const ri of allRoots) {
          if (ri.face === i || ri.face === j || ri.face === l) {
            let isShared = true;
            for (const other of others) {
              if (other === ri.face) continue;
              const signsO = signsAtRoots(P[ri.face], degP[ri.face],
                P[other], degP[other], 3, 0);
              if (ri.rootIdx < signsO.length && signsO[ri.rootIdx] !== 0)
                isShared = false;
            }
            if (isShared) {
              ri.isVertex = true;
              ri.isEdge = false;
              if (ri.face !== i) ri.valid = false;
            }
          }
        }
      }
    }
  }

  // Step 4b: Cv waypoint exclusion
  const lambdaPoly = [0n, 1n];
  for (const ri of allRoots) {
    if (!ri.valid) continue;
    if (!ri.isEdge && !ri.isVertex) continue;
    if (ri.isInfinity) continue;
    if (P[ri.face][0] !== 0n) continue;
    const signs = signsAtRoots(P[ri.face], degP[ri.face], lambdaPoly, 1, 3, 0);
    if (ri.rootIdx < signs.length && signs[ri.rootIdx] === 0)
      ri.valid = false;
  }

  // Step 4c: Edge/vertex pass-through exclusion
  {
    const dP = [];
    for (let k = 0; k < 4; k++) {
      const { fp, dfp } = polyDerivative(P[k], degP[k]);
      dP.push({ fp, dfp });
    }

    for (const ri of allRoots) {
      if (!ri.valid || ri.isInfinity) continue;

      if (ri.isEdge) {
        const k1 = ri.edgeFaces[0], k2 = ri.edgeFaces[1];
        if (k1 < 0 || k2 < 0) continue;
        // Product P'[k1] * P'[k2]
        const prod = new Array(8).fill(0n);
        for (let i = 0; i <= dP[k1].dfp; i++)
          for (let j = 0; j <= dP[k2].dfp; j++)
            prod[i + j] += dP[k1].fp[i] * dP[k2].fp[j];
        let dp = dP[k1].dfp + dP[k2].dfp;
        while (dp > 0 && prod[dp] === 0n) dp--;

        const sProd = signsAtRoots(P[ri.face], degP[ri.face], prod, dp, 3, 0);
        if (ri.rootIdx < sProd.length && sProd[ri.rootIdx] < 0) {
          ri.valid = false;
        } else if (ri.rootIdx < sProd.length && sProd[ri.rootIdx] === 0) {
          // Tangency at edge: check P''*Q
          for (let ki = 0; ki < 2; ki++) {
            const kf = ri.edgeFaces[ki];
            const { fp: fpp, dfp: dfpp } = polyDerivative(dP[kf].fp, dP[kf].dfp);
            const edeg = effectiveDegree(fpp, dfpp);
            const sDk = signsAtRoots(P[ri.face], degP[ri.face], dP[kf].fp, dP[kf].dfp, 3, 0);
            if (ri.rootIdx < sDk.length && sDk[ri.rootIdx] !== 0) continue;
            // P'' * Q_red product
            const pqProd = new Array(8).fill(0n);
            for (let a = 0; a <= edeg; a++)
              for (let b = 0; b <= degQ_red; b++)
                pqProd[a + b] += fpp[a] * Q_red[b];
            let dpq = edeg + degQ_red;
            while (dpq > 0 && pqProd[dpq] === 0n) dpq--;
            const sPQ = signsAtRoots(P[ri.face], degP[ri.face], pqProd, dpq, 3, 0);
            if (ri.rootIdx < sPQ.length && sPQ[ri.rootIdx] < 0) {
              ri.valid = false;
              break;
            }
          }
        }
      } else if (ri.isVertex) {
        const f = ri.face;
        const otherFaces = [];
        for (let j = 0; j < 4; j++) {
          if (j === f || otherFaces.length >= 2) continue;
          if (degP[j] === 0) continue;
          const sJ = signsAtRoots(P[f], degP[f], P[j], degP[j], 3, 0);
          if (ri.rootIdx < sJ.length && sJ[ri.rootIdx] === 0)
            otherFaces.push(j);
        }
        if (otherFaces.length === 2) {
          const k1 = otherFaces[0], k2 = otherFaces[1];
          const prod1 = new Array(8).fill(0n);
          for (let i = 0; i <= dP[f].dfp; i++)
            for (let j = 0; j <= dP[k1].dfp; j++)
              prod1[i + j] += dP[f].fp[i] * dP[k1].fp[j];
          let dp1 = dP[f].dfp + dP[k1].dfp;
          while (dp1 > 0 && prod1[dp1] === 0n) dp1--;
          const prod2 = new Array(8).fill(0n);
          for (let i = 0; i <= dP[f].dfp; i++)
            for (let j = 0; j <= dP[k2].dfp; j++)
              prod2[i + j] += dP[f].fp[i] * dP[k2].fp[j];
          let dp2 = dP[f].dfp + dP[k2].dfp;
          while (dp2 > 0 && prod2[dp2] === 0n) dp2--;
          const s1 = signsAtRoots(P[f], degP[f], prod1, dp1, 3, 0);
          const s2 = signsAtRoots(P[f], degP[f], prod2, dp2, 3, 0);
          if ((ri.rootIdx < s1.length && s1[ri.rootIdx] < 0) ||
              (ri.rootIdx < s2.length && s2[ri.rootIdx] < 0))
            ri.valid = false;
        }
      }
    }
  }

  // Step 4d: TN face-interior handling
  {
    const nAllOrig = allRoots.length;
    for (let r = 0; r < nAllOrig; r++) {
      const ri = allRoots[r];
      if (!ri.valid || ri.isInfinity || ri.isEdge || ri.isVertex) continue;

      const face = ri.face;
      const ridx = ri.rootIdx;
      if (degP[face] < 2) continue;

      const { fp: dPf, dfp: ddPf } = polyDerivative(P[face], degP[face]);
      const edeg = effectiveDegree(dPf, ddPf);
      if (edeg === 0 && dPf[0] === 0n) continue;
      const sDp = signsAtRoots(P[face], degP[face], dPf, edeg, 3, 0);
      if (ridx >= sDp.length || sDp[ridx] !== 0) continue;

      // TN candidate: P'(α)=0
      // P''[face] = [2*P[face][2], 6*P[face][3]] (matching C++ hardcoded)
      const Ppp = [2n * P[face][2], 6n * P[face][3]];
      const degPpp = effectiveDegree(Ppp, 1);
      if (degPpp === 0 && Ppp[0] === 0n) {
        ri.valid = false;
        result.hasPassthrough = true;
        continue;
      }

      // Evaluate P'' and Q separately at root, then multiply signs (matching C++)
      const sPP = signsAtRoots(P[face], degP[face], Ppp, degPpp, 3, 0);
      const sQ = signsAtRoots(P[face], degP[face], Q, degQ, 3, 0);
      const ppSign = (ridx < sPP.length) ? sPP[ridx] : 0;
      const qSign = (ridx < sQ.length) ? sQ[ridx] : 0;

      if (ppSign * qSign < 0) {
        ri.valid = false;
        result.hasPassthrough = true;
      } else if (ppSign * qSign > 0) {
        // Non-isolated TN → duplicate puncture
        allRoots.push({ ...ri });
      }
    }
  }

  // Step 5: Collect valid punctures
  const validRoots = allRoots.filter(ri => ri.valid);
  if (validRoots.length === 0) return result;

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
    else if (degQ_red === 2) {
      const d2 = Q_red[1] * Q_red[1] - 4n * Q_red[2] * Q_red[0];
      nQrRoots = (d2 > 0n) ? 2 : (d2 === 0n) ? 1 : 0;
    } else {
      const discQr = discriminantSign(Q_red);
      if (discQr > 0) nQrRoots = 3;
      else if (discQr < 0) nQrRoots = 1;
      else {
        const { sf, dsf } = polySqfree(Q_red, degQ_red);
        nQrRoots = dsf;
      }
    }
  }

  for (const ri of validRoots) {
    if (ri.isInfinity) {
      ri.qInterval = nQrRoots;
      continue;
    }
    if (nQrRoots === 0 || degQ_red === 0) {
      ri.qInterval = 0;
      continue;
    }

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
  // P is already reduced by h (Step 1), so P = P_red, degP = degP_red, nDistinct = nDistinctRed.
  for (const ri of validRoots) {
    result.punctures.push({
      face: ri.face,
      rootIdx: ri.rootIdx,
      qInterval: ri.qInterval,
      isEdge: ri.isEdge,
      isVertex: ri.isVertex,
      isInfinity: ri.isInfinity,
      edgeFaces: ri.edgeFaces,
    });
  }

  // Store infrastructure for classifier
  result.mergeInfinity = false;
  result.nQrRoots = nQrRoots;
  result.h = h;
  result.hDeg = dh;
  result.P_red = P.map(pk => pk.slice());
  result.degP_red = degP.slice();
  result.nDistinctRed = nDistinct.slice();

  // Step 9: RP1 Pairing
  let mergeInfinity = false;
  if (Q[3] === 0n) {
    mergeInfinity = true;
  } else {
    mergeInfinity = true;
    for (let k = 0; k < 4; k++) {
      if ((P[k][3] > 0n && Q[3] < 0n) || (P[k][3] < 0n && Q[3] > 0n)) {
        mergeInfinity = false;
        break;
      }
    }
  }
  result.mergeInfinity = mergeInfinity;

  const nPunctures = result.punctures.length;
  if (nPunctures >= 2) {
    // Build arrays for pairPuncturesRP1
    const sorted = [];
    const pFaceArr = [], pRidxArr = [], pQiArr = [];
    let nFin = 0;
    for (let i = 0; i < nPunctures; i++) {
      sorted.push(i);
      pFaceArr.push(result.punctures[i].face);
      pRidxArr.push(result.punctures[i].rootIdx);
      pQiArr.push(result.punctures[i].qInterval);
      if (result.punctures[i].rootIdx >= 0) nFin++;
    }

    const rp1Pairs = pairPuncturesRP1(4, sorted, nPunctures, nFin,
      pFaceArr, pRidxArr, pQiArr, nQrRoots, mergeInfinity,
      P, degP, nDistinct, Q_red, degQ_red);

    for (const pair of rp1Pairs) {
      result.pairs.push({ a: pair.a, b: pair.b, contains_infinity: pair.contains_infinity });
    }
  }

  return result;
}
