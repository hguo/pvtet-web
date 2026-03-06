/**
 * Sturm sequences for root isolation and sign certification.
 * Ported from ftk2 parallel_vector_solver.hpp.
 *
 * Two types:
 *   SturmSeqCubic  - for the cubic Q polynomial (degree 3)
 *   SturmSeqDeg4   - for the degree-4 barycentric N/D polynomials
 */

const EPS = Number.EPSILON;
const SQRT_EPS = Math.sqrt(EPS);

/** Evaluate polynomial with ascending-degree coefficients at x (Horner). */
function evalPoly(c, deg, x) {
  let r = c[deg];
  for (let i = deg - 1; i >= 0; i--) r = r * x + c[i];
  return r;
}

// ========================================================================
// Cubic Sturm Sequence (degree 3)
// ========================================================================

/**
 * Build Sturm sequence for cubic P (ascending-degree: P[0] + P[1]λ + P[2]λ² + P[3]λ³).
 * Returns { c: number[4][4], deg: number[4], n: number }
 */
export function buildSturmCubic(P) {
  const seq = {
    c: [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]],
    deg: [0, 0, 0, 0],
    n: 0,
  };

  const [p0, p1, p2, p3] = P;

  // S₀ = P
  seq.c[0] = [p0, p1, p2, p3];
  seq.deg[0] = 3;

  // S₁ = P'
  seq.c[1] = [p1, 2 * p2, 3 * p3, 0];
  seq.deg[1] = 2;

  // S₂ = -prem(S₀, S₁)
  const s20 = p3 * (p1 * p2 - 9 * p0 * p3);
  const s21 = 2 * p3 * (p2 * p2 - 3 * p1 * p3);
  seq.c[2] = [s20, s21, 0, 0];
  seq.deg[2] = (s21 !== 0.0) ? 1 : 0;

  // S₃ = -prem(S₁, S₂)
  const s30 = -(p1 * s21 * s21 - 2 * p2 * s21 * s20 + 3 * p3 * s20 * s20);
  seq.c[3] = [s30, 0, 0, 0];
  seq.deg[3] = 0;

  // Effective length
  if (s21 === 0.0 && s20 === 0.0) seq.n = 2;
  else if (s30 === 0.0) seq.n = 3;
  else seq.n = 4;

  return seq;
}

/** Count sign changes at x in the cubic Sturm sequence. */
export function sturmCountCubic(seq, x) {
  let changes = 0, prev = 0;
  for (let i = 0; i < seq.n; i++) {
    const v = evalPoly(seq.c[i], seq.deg[i], x);
    if (v !== 0.0) {
      if (prev !== 0.0 && ((prev > 0) !== (v > 0))) changes++;
      prev = v;
    }
  }
  return changes;
}

/**
 * Certified Sturm count at x with Higham error bound.
 * Returns { count, certified }.
 */
export function sturmCountCubicCertified(seq, x) {
  const ax = Math.abs(x);
  let changes = 0, certified = true, prev = 0;

  for (let i = 0; i < seq.n; i++) {
    const deg = seq.deg[i];
    const c = seq.c[i];

    // Horner evaluation
    let val = c[deg];
    for (let d = deg - 1; d >= 0; d--) val = val * x + c[d];

    // Higham condition number
    let cond = Math.abs(c[deg]);
    for (let d = deg - 1; d >= 0; d--) cond = cond * ax + Math.abs(c[d]);

    const gammaDeg = (2 * deg + 2) * EPS;

    if (val !== 0.0) {
      if (Math.abs(val) <= gammaDeg * cond) certified = false;
      if (prev !== 0.0 && ((prev > 0) !== (val > 0))) changes++;
      prev = val;
    } else {
      certified = false;
    }
  }

  return { count: changes, certified };
}

/**
 * Tighten a root into an isolating interval [lo, hi] via Sturm bisection.
 * Returns { lo, hi, success }.
 */
export function tightenRootInterval(seq, rf) {
  const scale = Math.max(Math.abs(rf), 1.0);
  let delta = scale * SQRT_EPS;
  let lo = rf - delta, hi = rf + delta;
  let cnt = sturmCountCubic(seq, lo) - sturmCountCubic(seq, hi);

  // Phase 1: expand or shrink until exactly 1 root in bracket
  for (let iter = 0; iter < 120 && cnt !== 1; iter++) {
    if (cnt === 0) delta *= 2;
    else delta *= 0.5;
    lo = rf - delta;
    hi = rf + delta;
    if (!isFinite(lo) || !isFinite(hi) || delta === 0) break;
    cnt = sturmCountCubic(seq, lo) - sturmCountCubic(seq, hi);
  }
  if (cnt !== 1) return { lo: rf, hi: rf, success: false };

  // Special: rf === 0.0 exactly
  if (rf === 0.0) return { lo: 0, hi: 0, success: false };

  // Phase 2: bisect to ULP convergence
  for (let iter = 0; iter < 200; iter++) {
    const mid = lo + (hi - lo) * 0.5;
    if (mid <= lo || mid >= hi) break;
    const halfCnt = sturmCountCubic(seq, lo) - sturmCountCubic(seq, mid);
    if (halfCnt === 1) hi = mid;
    else lo = mid;
  }

  return { lo, hi, success: true };
}

/**
 * Isolate cubic roots using Sturm sequences.
 * P[4] = float polynomial, floatRoots = initial root estimates.
 * Returns array of { lo, hi } intervals for isolated roots.
 */
export function isolateCubicRoots(P, floatRoots) {
  if (floatRoots.length === 0) return [];
  const seq = buildSturmCubic(P);
  const intervals = [];
  for (const rf of floatRoots) {
    const result = tightenRootInterval(seq, rf);
    intervals.push(result);
  }
  return intervals;
}

// ========================================================================
// Degree-4 Sturm Sequence (for barycentric N/D polynomials)
// ========================================================================

/** Float polynomial remainder: R = A mod B. Returns degree of R, or -1 if zero. */
function polyRemD(A, dA, B, dB) {
  const R = new Array(5).fill(0);
  for (let k = 0; k <= dA; k++) R[k] = A[k];

  for (let d = dA; d >= dB; d--) {
    if (R[d] === 0.0) continue;
    const coeff = R[d] / B[dB];
    const shift = d - dB;
    for (let i = 0; i <= dB; i++) R[i + shift] -= coeff * B[i];
    R[d] = 0.0;
  }

  let dR = dB - 1;
  while (dR > 0 && R[dR] === 0.0) dR--;
  return { R, deg: (R[dR] === 0.0 && dR === 0) ? -1 : dR };
}

/**
 * Build Sturm sequence for polynomial P of degree ≤ 4.
 * Returns { c: number[5][5], deg: number[5], n: number }
 */
export function buildSturmDeg4(P, degP) {
  const seq = {
    c: Array.from({ length: 5 }, () => new Array(5).fill(0)),
    deg: new Array(5).fill(0),
    n: 0,
  };

  // S₀ = P
  for (let k = 0; k <= degP; k++) seq.c[0][k] = P[k];
  seq.deg[0] = degP;
  seq.n = 1;
  if (degP === 0) return seq;

  // S₁ = P'
  for (let k = 1; k <= degP; k++) seq.c[1][k - 1] = k * seq.c[0][k];
  seq.deg[1] = degP - 1;
  seq.n = 2;
  if (degP === 1) return seq;

  // Sᵢ₊₁ = -rem(Sᵢ₋₁, Sᵢ)
  for (let i = 2; i <= degP && i <= 4; i++) {
    if (seq.deg[i - 1] === 0) break;
    const { R: rem, deg: dRem } = polyRemD(seq.c[i - 2], seq.deg[i - 2],
                                             seq.c[i - 1], seq.deg[i - 1]);
    if (dRem < 0) break;
    for (let k = 0; k <= dRem; k++) seq.c[i][k] = -rem[k];
    seq.deg[i] = dRem;
    seq.n = i + 1;
    if (dRem === 0) break;
  }

  return seq;
}

/** Sturm sign-change count at x for degree-4 Sturm sequence. */
export function sturmCountDeg4(seq, x) {
  let changes = 0, prev = 0;
  for (let i = 0; i < seq.n; i++) {
    const v = evalPoly(seq.c[i], seq.deg[i], x);
    if (v !== 0.0) {
      if (prev !== 0.0 && ((prev > 0) !== (v > 0))) changes++;
      prev = v;
    }
  }
  return changes;
}

// ========================================================================
// Barycentric numerator polynomials N[3][5] and D[5]
// ========================================================================

/** Multiply two degree-2 polynomials → degree-4. */
function mul2Poly(P, Q) {
  const R = [0, 0, 0, 0, 0];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      R[i + j] += P[i] * Q[j];
  return R;
}

/** Same but with BigInt. */
function mul2PolyBigInt(P, Q) {
  const R = [0n, 0n, 0n, 0n, 0n];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      R[i + j] += P[i] * Q[j];
  return R;
}

/**
 * Compute degree-4 barycentric numerator polynomials from integer Mlin/blin.
 * Uses BigInt for exact Gram matrix computation, then convert to float for
 * the degree-4 products (matching C++ compute_bary_numerators_from_integers).
 *
 * Mlin_q[3][2][2]: integer. M(λ)[r][c] = Mlin_q[r][c][0] + λ·Mlin_q[r][c][1]
 * blin_q[3][2]: integer.    b(λ)[r]    = blin_q[r][0]     + λ·blin_q[r][1]
 *
 * Returns { N: number[3][5], D: number[5] }
 */
export function computeBaryNumerators(Mlin_q, blin_q) {
  // Stage A: Gram matrix A[p][q][k] and rhs g[p][k] exactly in BigInt
  const A_i = Array.from({ length: 2 }, () =>
    Array.from({ length: 2 }, () => [0n, 0n, 0n])
  );
  for (let r = 0; r < 3; r++)
    for (let p = 0; p < 2; p++)
      for (let q = 0; q < 2; q++) {
        const m0p = BigInt(Mlin_q[r][p][0]), m1p = BigInt(Mlin_q[r][p][1]);
        const m0q = BigInt(Mlin_q[r][q][0]), m1q = BigInt(Mlin_q[r][q][1]);
        A_i[p][q][0] += m0p * m0q;
        A_i[p][q][1] += m0p * m1q + m1p * m0q;
        A_i[p][q][2] += m1p * m1q;
      }

  const g_i = Array.from({ length: 2 }, () => [0n, 0n, 0n]);
  for (let r = 0; r < 3; r++)
    for (let p = 0; p < 2; p++) {
      const m0 = BigInt(Mlin_q[r][p][0]), m1 = BigInt(Mlin_q[r][p][1]);
      const b0 = BigInt(blin_q[r][0]), b1 = BigInt(blin_q[r][1]);
      g_i[p][0] += m0 * b0;
      g_i[p][1] += m0 * b1 + m1 * b0;
      g_i[p][2] += m1 * b1;
    }

  // Stage B: degree-4 products.
  // For small integer inputs, everything fits in BigInt. Do exact path.
  const t0 = mul2PolyBigInt(A_i[0][0], A_i[1][1]);
  const t1 = mul2PolyBigInt(A_i[0][1], A_i[0][1]);
  const D_i = [0n, 0n, 0n, 0n, 0n];
  for (let k = 0; k < 5; k++) D_i[k] = t0[k] - t1[k];

  const a11g0 = mul2PolyBigInt(A_i[1][1], g_i[0]);
  const a01g1 = mul2PolyBigInt(A_i[0][1], g_i[1]);
  const N_i = Array.from({ length: 3 }, () => [0n, 0n, 0n, 0n, 0n]);
  for (let k = 0; k < 5; k++) N_i[0][k] = a11g0[k] - a01g1[k];

  const a00g1 = mul2PolyBigInt(A_i[0][0], g_i[1]);
  const a01g0 = mul2PolyBigInt(A_i[0][1], g_i[0]);
  for (let k = 0; k < 5; k++) N_i[1][k] = a00g1[k] - a01g0[k];

  for (let k = 0; k < 5; k++) N_i[2][k] = D_i[k] - N_i[0][k] - N_i[1][k];

  // Convert to double
  const D = D_i.map(Number);
  const N = N_i.map(row => row.map(Number));

  return { N, D };
}
