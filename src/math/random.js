/**
 * Seeded PRNG matching ftk2 pv_tet_case_finder.cu LCG.
 *
 * LCG: state = state * 1664525 + 1013904223 (mod 2^32)
 * Seed: global_id ^ (base_seed * 2654435761) with 4 warmup steps
 * Default base_seed: 42
 */

const BASE_SEED = 42;

export class SeededRNG {
  constructor(seed = 0) {
    // Initialize: global_id XOR (base_seed * 2654435761ULL)
    // In JS, use BigInt for the 64-bit multiply, then truncate to 32 bits
    const init = BigInt(seed) ^ (BigInt(BASE_SEED) * 2654435761n);
    this.state = Number(init & 0xFFFFFFFFn) >>> 0;

    // Warmup: 4 LCG iterations
    for (let i = 0; i < 4; i++) this.next();
  }

  next() {
    // LCG with Numerical Recipes constants
    this.state = (Math.imul(this.state, 1664525) + 1013904223) >>> 0;
    return this.state;
  }

  /** Return integer in [-R, R] inclusive, matching ftk2 rand_int_dev. */
  randInt(lo, hi) {
    const R = hi; // assumes symmetric range [-R, R]
    const r = this.next();
    return (r % (2 * R + 1)) - R;
  }
}

/**
 * Generate random integer V[4][3] and W[4][3] with values in [-range, range].
 * Matches ftk2 pv_tet_case_finder.cu exactly.
 *
 * Seed 0: analytically constructed bubble case (T0_Q2-_Cv_B)
 * Seed 1: non-isolated SR case (Q3+ with gcd(Q,P_2) deg 2)
 */
export function randomVW(seed, range = 9) {
  if (seed === 0) {
    return {
      V: [[2,3,-1],[-1,-2,-1],[0,-1,2],[-2,-2,0]],
      W: [[3,0,3],[-1,0,-1],[-3,-2,-1],[0,3,-3]],
    };
  }
  if (seed === 1) {
    return {
      V: [[5,0,1],[3,2,4],[1,2,-2],[2,4,4]],
      W: [[0,1,3],[3,-5,3],[0,4,5],[1,1,0]],
    };
  }
  const rng = new SeededRNG(seed);
  const V = Array.from({ length: 4 }, () =>
    Array.from({ length: 3 }, () => rng.randInt(-range, range))
  );
  const W = Array.from({ length: 4 }, () =>
    Array.from({ length: 3 }, () => rng.randInt(-range, range))
  );
  return { V, W };
}
