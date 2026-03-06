/**
 * Polynomial computation for PV tetrahedron: Q(λ) and P_k(λ).
 * Ported from ftk2 parallel_vector_solver.hpp.
 *
 * Conventions: polynomials are arrays [c0, c1, c2, c3]
 * representing c0 + c1*λ + c2*λ² + c3*λ³.
 */

/** Evaluate polynomial at x using Horner's method. */
export function polyEval(coeffs, x) {
  let val = 0;
  for (let i = coeffs.length - 1; i >= 0; i--) {
    val = val * x + coeffs[i];
  }
  return val;
}

/** 3x3 determinant. */
function det3(A) {
  return (
    A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
    A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
    A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
  );
}

/** Characteristic polynomial of 2x2: det([[a00-λb00, a01-λb01],[a10-λb10, a11-λb11]]) */
function charPoly2x2(a00, a01, a10, a11, b00, b01, b10, b11) {
  return [
    a00 * a11 - a10 * a01,                                           // P[0]
    -(a00 * b11 - a10 * b01 + b00 * a11 - b10 * a01),               // P[1]
    b00 * b11 - b10 * b01,                                           // P[2]
  ];
}

/** Characteristic polynomial of 3x3: det(A - λB) = P[0] + P[1]λ + P[2]λ² + P[3]λ³ */
function charPoly3x3(A, B) {
  const P0 = det3(A);

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

  const P3 = -det3(B);

  return [P0, P1, P2, P3];
}

/** Multiply two polynomials (degree m and n). */
function polyMul(P, m, Q, n) {
  const R = new Array(m + n + 1).fill(0);
  for (let i = 0; i <= m; i++)
    for (let j = 0; j <= n; j++)
      R[i + j] += P[i] * Q[j];
  return R;
}

/**
 * Compute Q(λ) and P[4](λ) from V[4][3], W[4][3].
 * Returns { Q: number[4], P: number[4][4] }
 */
export function characteristicPolynomials(V, W) {
  // Edge-difference matrices
  const A = [
    [V[0][0] - V[3][0], V[1][0] - V[3][0], V[2][0] - V[3][0]],
    [V[0][1] - V[3][1], V[1][1] - V[3][1], V[2][1] - V[3][1]],
    [V[0][2] - V[3][2], V[1][2] - V[3][2], V[2][2] - V[3][2]],
  ];
  const B = [
    [W[0][0] - W[3][0], W[1][0] - W[3][0], W[2][0] - W[3][0]],
    [W[0][1] - W[3][1], W[1][1] - W[3][1], W[2][1] - W[3][1]],
    [W[0][2] - W[3][2], W[1][2] - W[3][2], W[2][2] - W[3][2]],
  ];

  const a = [V[3][0], V[3][1], V[3][2]];
  const b = [W[3][0], W[3][1], W[3][2]];

  // rhs[i] = -a[i] + b[i]*λ  (linear polynomial)
  const rhs = [
    [-a[0], b[0]],
    [-a[1], b[1]],
    [-a[2], b[2]],
  ];

  // Q = det(A - λB)
  const Q = charPoly3x3(A, B);

  // Adjugate matrix of (A - λB), each entry is degree-2 polynomial
  const adj = Array.from({ length: 3 }, () =>
    Array.from({ length: 3 }, () => [0, 0, 0])
  );

  adj[0][0] = charPoly2x2(A[1][1], A[1][2], A[2][1], A[2][2], B[1][1], B[1][2], B[2][1], B[2][2]);
  adj[1][0] = charPoly2x2(A[1][0], A[1][2], A[2][0], A[2][2], B[1][0], B[1][2], B[2][0], B[2][2]);
  adj[2][0] = charPoly2x2(A[1][0], A[1][1], A[2][0], A[2][1], B[1][0], B[1][1], B[2][0], B[2][1]);
  adj[0][1] = charPoly2x2(A[0][1], A[0][2], A[2][1], A[2][2], B[0][1], B[0][2], B[2][1], B[2][2]);
  adj[1][1] = charPoly2x2(A[0][0], A[0][2], A[2][0], A[2][2], B[0][0], B[0][2], B[2][0], B[2][2]);
  adj[2][1] = charPoly2x2(A[0][0], A[0][1], A[2][0], A[2][1], B[0][0], B[0][1], B[2][0], B[2][1]);
  adj[0][2] = charPoly2x2(A[0][1], A[0][2], A[1][1], A[1][2], B[0][1], B[0][2], B[1][1], B[1][2]);
  adj[1][2] = charPoly2x2(A[0][0], A[0][2], A[1][0], A[1][2], B[0][0], B[0][2], B[1][0], B[1][2]);
  adj[2][2] = charPoly2x2(A[0][0], A[0][1], A[1][0], A[1][1], B[0][0], B[0][1], B[1][0], B[1][1]);

  // Fix signs for adjugate (cofactor sign pattern)
  for (const c of adj[0][1]) { /* negate inline below */ }
  adj[0][1] = adj[0][1].map(x => -x);
  adj[1][0] = adj[1][0].map(x => -x);
  adj[1][2] = adj[1][2].map(x => -x);
  adj[2][1] = adj[2][1].map(x => -x);

  // P[i] = sum_j adj[i][j] * rhs[j]  (transpose: adj[i][j]^T * rhs)
  const P = Array.from({ length: 4 }, () => [0, 0, 0, 0]);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const prod = polyMul(adj[i][j], 2, rhs[j], 1);
      for (let k = 0; k < 4; k++) {
        P[i][k] += (prod[k] || 0);
      }
    }
  }

  // P[3] = Q - P[0] - P[1] - P[2]
  for (let k = 0; k < 4; k++) {
    P[3][k] = Q[k] - P[0][k] - P[1][k] - P[2][k];
  }

  return { Q, P };
}

/**
 * Format polynomial coefficients as LaTeX string.
 * coeffs = [c0, c1, c2, c3], name = "Q"
 */
export function polyToLatex(coeffs, name = 'Q') {
  const terms = [];
  for (let i = coeffs.length - 1; i >= 0; i--) {
    const c = Math.round(coeffs[i]);
    if (c === 0) continue;

    let varPart = '';
    if (i === 0) varPart = '';
    else if (i === 1) varPart = '\\lambda';
    else varPart = `\\lambda^{${i}}`;

    let term;
    if (terms.length === 0) {
      // First term
      if (Math.abs(c) === 1 && i > 0) {
        term = c < 0 ? `-${varPart}` : varPart;
      } else {
        term = `${c}${varPart}`;
      }
    } else {
      if (Math.abs(c) === 1 && i > 0) {
        term = c < 0 ? ` - ${varPart}` : ` + ${varPart}`;
      } else {
        term = c < 0 ? ` - ${Math.abs(c)}${varPart}` : ` + ${c}${varPart}`;
      }
    }
    terms.push(term);
  }

  const expr = terms.length > 0 ? terms.join('') : '0';
  return `${name}(\\lambda) = ${expr}`;
}
