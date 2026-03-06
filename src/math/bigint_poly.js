/**
 * Exact integer polynomial computation using BigInt (replacing __int128).
 * Ported from ftk2 parallel_vector_solver.hpp.
 *
 * Since web demo inputs are already small integers ([-50, 50]),
 * no quantization is needed — we use them directly as BigInts.
 */

/** 3x3 integer determinant using BigInt. */
export function det3_bigint(a) {
  return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
       - a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2])
       + a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
}

/**
 * Exact characteristic polynomial det(A - λB) with BigInt coefficients.
 * A[3][3], B[3][3] are BigInt matrices.
 * Returns P[4] (BigInt): P[0] + P[1]λ + P[2]λ² + P[3]λ³
 */
export function charPoly3x3_bigint(A, B) {
  const P0 = det3_bigint(A);

  const P1 =
    A[1][2]*A[2][1]*B[0][0] - A[1][1]*A[2][2]*B[0][0] - A[1][2]*A[2][0]*B[0][1]
   +A[1][0]*A[2][2]*B[0][1] + A[1][1]*A[2][0]*B[0][2] - A[1][0]*A[2][1]*B[0][2]
   -A[0][2]*A[2][1]*B[1][0] + A[0][1]*A[2][2]*B[1][0] + A[0][2]*A[2][0]*B[1][1]
   -A[0][0]*A[2][2]*B[1][1] - A[0][1]*A[2][0]*B[1][2] + A[0][0]*A[2][1]*B[1][2]
   +A[0][2]*A[1][1]*B[2][0] - A[0][1]*A[1][2]*B[2][0] - A[0][2]*A[1][0]*B[2][1]
   +A[0][0]*A[1][2]*B[2][1] + A[0][1]*A[1][0]*B[2][2] - A[0][0]*A[1][1]*B[2][2];

  const P2 =
   -A[2][2]*B[0][1]*B[1][0] + A[2][1]*B[0][2]*B[1][0] + A[2][2]*B[0][0]*B[1][1]
   -A[2][0]*B[0][2]*B[1][1] - A[2][1]*B[0][0]*B[1][2] + A[2][0]*B[0][1]*B[1][2]
   +A[1][2]*B[0][1]*B[2][0] - A[1][1]*B[0][2]*B[2][0] - A[0][2]*B[1][1]*B[2][0]
   +A[0][1]*B[1][2]*B[2][0] - A[1][2]*B[0][0]*B[2][1] + A[1][0]*B[0][2]*B[2][1]
   +A[0][2]*B[1][0]*B[2][1] - A[0][0]*B[1][2]*B[2][1] + A[1][1]*B[0][0]*B[2][2]
   -A[1][0]*B[0][1]*B[2][2] - A[0][1]*B[1][0]*B[2][2] + A[0][0]*B[1][1]*B[2][2];

  const P3 = -det3_bigint(B);

  return [P0, P1, P2, P3];
}

/** Convert int arrays to BigInt 3x3 matrix. */
export function toBigInt3x3(M) {
  return M.map(row => row.map(v => BigInt(v)));
}

/**
 * Exact discriminant sign of cubic P[0]+P[1]λ+P[2]λ²+P[3]λ³ using BigInt.
 * Returns +1 (3 real roots), -1 (1 real root), or 0 (repeated/degenerate).
 *
 * Standard discriminant: Δ = 18abcd - 4b³d + b²c² - 4ac³ - 27a²d²
 * where a=P[0], b=P[1], c=P[2], d=P[3]
 */
export function discriminantSign_bigint(P) {
  if (P[3] === 0n) return 0;

  // GCD-normalize to reduce coefficient magnitudes
  let g = bigAbs(P[0]);
  for (let i = 1; i < 4; i++) g = bigGcd(g, bigAbs(P[i]));
  if (g === 0n) return 0;

  const a = P[0] / g;
  const b = P[1] / g;
  const c = P[2] / g;
  const d = P[3] / g;

  // Exact discriminant in BigInt (no overflow possible)
  const disc = 18n*a*b*c*d - 4n*b*b*b*d + b*b*c*c - 4n*a*c*c*c - 27n*a*a*d*d;

  if (disc > 0n) return 1;   // 3 distinct real roots
  if (disc < 0n) return -1;  // 1 real root
  return 0;                   // repeated root
}

/** GCD of two BigInt values (both non-negative). */
function bigGcd(a, b) {
  while (b !== 0n) { const t = b; b = a % b; a = t; }
  return a;
}

/** Absolute value of BigInt. */
function bigAbs(x) {
  return x < 0n ? -x : x;
}

/**
 * Compute tet-level Q and P polynomials in BigInt.
 * V[4][3], W[4][3] are integer arrays.
 * Returns { Q_bigint: BigInt[4], P_bigint: BigInt[4][4] }
 */
export function characteristicPolynomials_bigint(V, W) {
  // Edge-difference matrices (BigInt)
  const A = [
    [BigInt(V[0][0] - V[3][0]), BigInt(V[1][0] - V[3][0]), BigInt(V[2][0] - V[3][0])],
    [BigInt(V[0][1] - V[3][1]), BigInt(V[1][1] - V[3][1]), BigInt(V[2][1] - V[3][1])],
    [BigInt(V[0][2] - V[3][2]), BigInt(V[1][2] - V[3][2]), BigInt(V[2][2] - V[3][2])],
  ];
  const B = [
    [BigInt(W[0][0] - W[3][0]), BigInt(W[1][0] - W[3][0]), BigInt(W[2][0] - W[3][0])],
    [BigInt(W[0][1] - W[3][1]), BigInt(W[1][1] - W[3][1]), BigInt(W[2][1] - W[3][1])],
    [BigInt(W[0][2] - W[3][2]), BigInt(W[1][2] - W[3][2]), BigInt(W[2][2] - W[3][2])],
  ];

  const a = [BigInt(V[3][0]), BigInt(V[3][1]), BigInt(V[3][2])];
  const b = [BigInt(W[3][0]), BigInt(W[3][1]), BigInt(W[3][2])];

  // Q = det(A - λB)
  const Q = charPoly3x3_bigint(A, B);

  // rhs[i] = -a[i] + b[i]*λ (linear polynomials in BigInt)
  const rhs = [[-a[0], b[0]], [-a[1], b[1]], [-a[2], b[2]]];

  // Adjugate of (A - λB): each entry is degree-2 polynomial in BigInt
  function cp2x2(a00, a01, a10, a11, b00, b01, b10, b11) {
    return [
      a00 * a11 - a10 * a01,
      -(a00 * b11 - a10 * b01 + b00 * a11 - b10 * a01),
      b00 * b11 - b10 * b01,
    ];
  }

  const adj = Array.from({ length: 3 }, () =>
    Array.from({ length: 3 }, () => [0n, 0n, 0n])
  );

  adj[0][0] = cp2x2(A[1][1], A[1][2], A[2][1], A[2][2], B[1][1], B[1][2], B[2][1], B[2][2]);
  adj[1][0] = cp2x2(A[1][0], A[1][2], A[2][0], A[2][2], B[1][0], B[1][2], B[2][0], B[2][2]);
  adj[2][0] = cp2x2(A[1][0], A[1][1], A[2][0], A[2][1], B[1][0], B[1][1], B[2][0], B[2][1]);
  adj[0][1] = cp2x2(A[0][1], A[0][2], A[2][1], A[2][2], B[0][1], B[0][2], B[2][1], B[2][2]);
  adj[1][1] = cp2x2(A[0][0], A[0][2], A[2][0], A[2][2], B[0][0], B[0][2], B[2][0], B[2][2]);
  adj[2][1] = cp2x2(A[0][0], A[0][1], A[2][0], A[2][1], B[0][0], B[0][1], B[2][0], B[2][1]);
  adj[0][2] = cp2x2(A[0][1], A[0][2], A[1][1], A[1][2], B[0][1], B[0][2], B[1][1], B[1][2]);
  adj[1][2] = cp2x2(A[0][0], A[0][2], A[1][0], A[1][2], B[0][0], B[0][2], B[1][0], B[1][2]);
  adj[2][2] = cp2x2(A[0][0], A[0][1], A[1][0], A[1][1], B[0][0], B[0][1], B[1][0], B[1][1]);

  // Fix signs for adjugate
  adj[0][1] = adj[0][1].map(x => -x);
  adj[1][0] = adj[1][0].map(x => -x);
  adj[1][2] = adj[1][2].map(x => -x);
  adj[2][1] = adj[2][1].map(x => -x);

  // P[i] = sum_j adj[i][j] * rhs[j] (polynomial multiplication + addition)
  const P = Array.from({ length: 4 }, () => [0n, 0n, 0n, 0n]);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      // Multiply degree-2 adj[i][j] by degree-1 rhs[j] → degree-3
      const prod = [0n, 0n, 0n, 0n];
      for (let a = 0; a < 3; a++)
        for (let b = 0; b < 2; b++)
          prod[a + b] += adj[i][j][a] * rhs[j][b];
      for (let k = 0; k < 4; k++) P[i][k] += prod[k];
    }
  }

  // P[3] = Q - P[0] - P[1] - P[2]
  for (let k = 0; k < 4; k++) {
    P[3][k] = Q[k] - P[0][k] - P[1][k] - P[2][k];
  }

  return { Q_bigint: Q, P_bigint: P };
}

/**
 * GCD-normalize a BigInt[4] polynomial to fit safely in float64.
 * Returns a float64[4] array.
 */
export function gcdNormalizePoly(P_bigint) {
  let g = bigAbs(P_bigint[0]);
  for (let i = 1; i < P_bigint.length; i++) g = bigGcd(g, bigAbs(P_bigint[i]));
  if (g === 0n) g = 1n;
  return P_bigint.map(c => Number(c / g));
}

/**
 * Shared-root detection via Sylvester resultant (exact BigInt).
 * Returns true iff Q and some P[k] share a common root.
 * Q_bigint[4], P_bigint[4][4] are BigInt polynomial coefficients.
 */
export function hasSharedRoot_bigint(Q_bigint, P_bigint) {
  let degQ = 3;
  while (degQ > 0 && Q_bigint[degQ] === 0n) degQ--;
  if (degQ === 0) return false;

  for (let k = 0; k < 4; k++) {
    let degP = 3;
    while (degP > 0 && P_bigint[k][degP] === 0n) degP--;
    if (degP === 0) continue;

    const N = degQ + degP;
    // Build Sylvester matrix (N x N)
    const M = Array.from({ length: N }, () => new Array(N).fill(0n));
    for (let i = 0; i < degP; i++)
      for (let j = 0; j <= degQ; j++)
        M[i][i + degQ - j] = Q_bigint[j];
    for (let i = 0; i < degQ; i++)
      for (let j = 0; j <= degP; j++)
        M[degP + i][i + degP - j] = P_bigint[k][j];

    // Bareiss fraction-free elimination (exact BigInt determinant)
    let prevPivot = 1n;
    let zeroDet = false;
    for (let col = 0; col < N; col++) {
      let pivot = -1;
      for (let row = col; row < N; row++)
        if (M[row][col] !== 0n) { pivot = row; break; }
      if (pivot < 0) { zeroDet = true; break; }
      if (pivot !== col)
        for (let j = 0; j < N; j++)
          [M[col][j], M[pivot][j]] = [M[pivot][j], M[col][j]];
      for (let row = col + 1; row < N; row++) {
        for (let j = col + 1; j < N; j++)
          M[row][j] = (M[col][col] * M[row][j] - M[row][col] * M[col][j]) / prevPivot;
        M[row][col] = 0n;
      }
      prevPivot = M[col][col];
    }
    if (zeroDet || M[N - 1][N - 1] === 0n) return true;
  }
  return false;
}

/**
 * Compute characteristic polynomial for a triangle face.
 * V_face[3][3], W_face[3][3] are integer arrays (vertex × component).
 * Returns { P_bigint: BigInt[4], P_float: number[4] }
 *
 * The polynomial is det(VT - λ·WT) where VT = transpose of V_face.
 */
export function triangleCharPoly(V_face, W_face) {
  // Transpose: VT[component][vertex]
  const VT = Array.from({ length: 3 }, (_, i) =>
    V_face.map(row => BigInt(row[i]))
  );
  const WT = Array.from({ length: 3 }, (_, i) =>
    W_face.map(row => BigInt(row[i]))
  );

  const P_bigint = charPoly3x3_bigint(VT, WT);

  // Float version (for root solving)
  const VTf = Array.from({ length: 3 }, (_, i) =>
    V_face.map(row => row[i])
  );
  const WTf = Array.from({ length: 3 }, (_, i) =>
    W_face.map(row => row[i])
  );

  // Float char poly
  function det3f(a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
  }
  const A = VTf, B = WTf;
  const P0 = det3f(A);
  const P1 =
    A[1][2]*A[2][1]*B[0][0] - A[1][1]*A[2][2]*B[0][0] - A[1][2]*A[2][0]*B[0][1]
   +A[1][0]*A[2][2]*B[0][1] + A[1][1]*A[2][0]*B[0][2] - A[1][0]*A[2][1]*B[0][2]
   -A[0][2]*A[2][1]*B[1][0] + A[0][1]*A[2][2]*B[1][0] + A[0][2]*A[2][0]*B[1][1]
   -A[0][0]*A[2][2]*B[1][1] - A[0][1]*A[2][0]*B[1][2] + A[0][0]*A[2][1]*B[1][2]
   +A[0][2]*A[1][1]*B[2][0] - A[0][1]*A[1][2]*B[2][0] - A[0][2]*A[1][0]*B[2][1]
   +A[0][0]*A[1][2]*B[2][1] + A[0][1]*A[1][0]*B[2][2] - A[0][0]*A[1][1]*B[2][2];
  const P2 =
   -A[2][2]*B[0][1]*B[1][0] + A[2][1]*B[0][2]*B[1][0] + A[2][2]*B[0][0]*B[1][1]
   -A[2][0]*B[0][2]*B[1][1] - A[2][1]*B[0][0]*B[1][2] + A[2][0]*B[0][1]*B[1][2]
   +A[1][2]*B[0][1]*B[2][0] - A[1][1]*B[0][2]*B[2][0] - A[0][2]*B[1][1]*B[2][0]
   +A[0][1]*B[1][2]*B[2][0] - A[1][2]*B[0][0]*B[2][1] + A[1][0]*B[0][2]*B[2][1]
   +A[0][2]*B[1][0]*B[2][1] - A[0][0]*B[1][2]*B[2][1] + A[1][1]*B[0][0]*B[2][2]
   -A[1][0]*B[0][1]*B[2][2] - A[0][1]*B[1][0]*B[2][2] + A[0][0]*B[1][1]*B[2][2];
  const P3 = -det3f(B);

  return { P_bigint, P_float: [P0, P1, P2, P3] };
}

/**
 * Check if V[i] × W[i] = 0 at all vertices (exact BigInt cross product).
 * V[n][3], W[n][3] are integer arrays.
 */
export function allParallel(V, W) {
  for (let i = 0; i < V.length; i++) {
    const cx = BigInt(V[i][1]) * BigInt(W[i][2]) - BigInt(V[i][2]) * BigInt(W[i][1]);
    const cy = BigInt(V[i][2]) * BigInt(W[i][0]) - BigInt(V[i][0]) * BigInt(W[i][2]);
    const cz = BigInt(V[i][0]) * BigInt(W[i][1]) - BigInt(V[i][1]) * BigInt(W[i][0]);
    if (cx !== 0n || cy !== 0n || cz !== 0n) return false;
  }
  return true;
}
