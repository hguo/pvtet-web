/**
 * Cubic solver and root utilities.
 * Ported from ftk2 parallel_vector_solver.hpp solve_cubic_real.
 */

/**
 * Solve cubic P[0] + P[1]x + P[2]x² + P[3]x³ = 0.
 * Returns sorted array of real roots.
 */
export function solveCubic(P) {
  const eps = 1e-12;

  if (Math.abs(P[3]) < eps) {
    // Degenerate: quadratic or lower
    if (Math.abs(P[2]) < eps) {
      if (Math.abs(P[1]) < eps) return [];
      return [-P[0] / P[1]];
    }
    const disc = P[1] * P[1] - 4 * P[2] * P[0];
    if (disc < -eps) return [];
    if (Math.abs(disc) < eps) return [-P[1] / (2 * P[2])];
    const sd = Math.sqrt(disc);
    const roots = [(-P[1] + sd) / (2 * P[2]), (-P[1] - sd) / (2 * P[2])];
    roots.sort((a, b) => a - b);
    return roots;
  }

  // Normalize: x³ + bx² + cx + d = 0
  const b = P[2] / P[3];
  const c = P[1] / P[3];
  const d = P[0] / P[3];

  const q = (3 * c - b * b) / 9;
  const r = (-(27 * d) + b * (9 * c - 2 * b * b)) / 54;
  const disc = q * q * q + r * r;
  const term1 = b / 3;

  if (disc > eps) {
    // One real root (Cardano)
    let s = r + Math.sqrt(disc);
    s = s < 0 ? -Math.pow(-s, 1 / 3) : Math.pow(s, 1 / 3);
    let t = r - Math.sqrt(disc);
    t = t < 0 ? -Math.pow(-t, 1 / 3) : Math.pow(t, 1 / 3);
    return [-term1 + s + t];
  } else if (Math.abs(disc) < eps) {
    // Two or three equal roots
    const r13 = r < 0 ? -Math.pow(-r, 1 / 3) : Math.pow(r, 1 / 3);
    const roots = [-term1 + 2 * r13, -(r13 + term1)];
    if (Math.abs(roots[0] - roots[1]) < eps) return [roots[0]];
    roots.sort((a, b) => a - b);
    return roots;
  } else {
    // Three distinct real roots (trigonometric)
    const qNeg = -q;
    const dum1 = Math.acos(r / Math.sqrt(qNeg * qNeg * qNeg));
    const r13 = 2 * Math.sqrt(qNeg);
    const roots = [
      -term1 + r13 * Math.cos(dum1 / 3),
      -term1 + r13 * Math.cos((dum1 + 2 * Math.PI) / 3),
      -term1 + r13 * Math.cos((dum1 + 4 * Math.PI) / 3),
    ];
    roots.sort((a, b) => a - b);
    return roots;
  }
}

/**
 * Compute discriminant sign of cubic: +1 (3 roots), -1 (1 root), 0 (repeated).
 */
export function cubicDiscriminantSign(P) {
  if (Math.abs(P[3]) < 1e-12) return 0;
  const b = P[2] / P[3];
  const c = P[1] / P[3];
  const d = P[0] / P[3];
  const q = (3 * c - b * b) / 9;
  const r = (-(27 * d) + b * (9 * c - 2 * b * b)) / 54;
  const disc = q * q * q + r * r;
  if (disc > 1e-12) return -1; // 1 root
  if (disc < -1e-12) return 1; // 3 roots
  return 0;
}

/**
 * Effective degree of polynomial.
 */
export function polyDegree(coeffs) {
  let deg = coeffs.length - 1;
  while (deg > 0 && Math.abs(coeffs[deg]) < 1e-12) deg--;
  return deg;
}
