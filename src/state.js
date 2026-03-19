/**
 * Reactive state management with event bus.
 * Central store for V, W, computed results, and hover state.
 *
 * Uses classifyTetCase() for full combinatorial classification and stitching,
 * then builds visual segments via curve sampling.
 */
import { classifyTetCase } from './math/classify.js';
import { classifyTriCase2D } from './math/classify2d.js';
import { polyEval } from './math/polynomials.js';
import {
  SEGMENT_COLORS,
  baryTetTo3D, lambdaToBary, samplePVCurve, sampleBubble,
} from './math/curves.js';
import {
  SEGMENT_COLORS_2D,
  baryTriTo2D, lambdaToBary2D, samplePVCurve2D, sampleBubble2D,
} from './math/curves2d.js';
import { randomVW } from './math/random.js';

// Default V, W: seed 0 (bubble case, matching ftk2)
const { V: DEFAULT_V, W: DEFAULT_W } = randomVW(0, 20);

class State {
  constructor() {
    this._listeners = {};
    this.dim = '3d'; // '2d' or '3d'
    this.V = DEFAULT_V.map(r => [...r]);
    this.W = DEFAULT_W.map(r => [...r]);
    // Separate 2D field storage (seed 0 defaults)
    this.V2d = [[2,3],[-1,-2],[0,-1]];
    this.W2d = [[3,0],[-1,0],[-3,-2]];

    // Computed fields (populated by recompute)
    this.Q = [0, 0, 0, 0];
    this.P = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]];
    this.qRoots = [];
    this.qDiscSign = 0;
    this.qDegree = 3;
    this.segments = [];
    this.punctures = [];
    this.intervals = [];
    this.bubble = null;

    // Classification info
    this.category = '';
    this.hasSR = false;
    this.hasB = false;
    this.hasCvPos = false;
    this.hasCwPos = false;
    this.CvMu = null;
    this.CwMu = null;
    this.srLambda = null;
    this.srPos3d = null;

    // Hover state
    this.hoverPunctureIdx = -1;
    this.hoverSegmentIdx = -1;
    this.hoverLambda = null;
    this.hoverRingPos3d = null;

    this.recompute();
  }

  on(event, fn) {
    if (!this._listeners[event]) this._listeners[event] = [];
    this._listeners[event].push(fn);
  }

  off(event, fn) {
    if (!this._listeners[event]) return;
    this._listeners[event] = this._listeners[event].filter(f => f !== fn);
  }

  emit(event, data) {
    if (this._listeners[event]) {
      for (const fn of this._listeners[event]) fn(data);
    }
  }

  setDim(dim) {
    if (this.dim === dim) return;
    this.dim = dim;
    this.recompute();
    this.emit('modeChanged');
    this.emit('dataChanged');
  }

  setVW(V, W) {
    if (this.dim === '2d') {
      this.V2d = V.map(r => [...r]);
      this.W2d = W.map(r => [...r]);
    } else {
      this.V = V.map(r => [...r]);
      this.W = W.map(r => [...r]);
    }
    this.recompute();
    this.emit('dataChanged');
  }

  setHoverPuncture(idx) {
    if (this.hoverPunctureIdx === idx) return;
    this.hoverPunctureIdx = idx;
    this.hoverSegmentIdx = -1;
    this.hoverLambda = null;
    this.hoverRingPos3d = null;
    this.emit('hoverChanged');
  }

  setHoverSegment(idx) {
    if (this.hoverSegmentIdx === idx) return;
    this.hoverSegmentIdx = idx;
    this.hoverPunctureIdx = -1;
    this.hoverLambda = null;
    this.hoverRingPos3d = null;
    this.emit('hoverChanged');
  }

  setHoverLambda(lam) {
    this.hoverLambda = lam;
    this.hoverPunctureIdx = -1;
    this.hoverSegmentIdx = -1;
    if (lam !== null) {
      if (this.dim === '2d') {
        let mu;
        if (!isFinite(lam)) {
          const degQ = this.Q.length - 1;
          if (this.Q[degQ] !== 0) {
            mu = this.P.map(pk => (pk.length > degQ ? pk[degQ] : 0) / this.Q[degQ]);
          }
        } else {
          mu = lambdaToBary2D(lam, this.Q, this.P);
        }
        if (mu && mu.every(m => m >= -1e-6) && mu.every(m => m <= 1 + 1e-6)) {
          this.hoverRingPos3d = baryTriTo2D(mu);
        } else {
          this.hoverRingPos3d = null;
        }
      } else {
        let mu;
        if (!isFinite(lam)) {
          if (this.Q[3] !== 0) {
            mu = this.P.map(pk => pk[3] / this.Q[3]);
          }
        } else {
          mu = lambdaToBary(lam, this.Q, this.P);
        }
        if (mu && mu.every(m => m >= -1e-6) && mu.every(m => m <= 1 + 1e-6)) {
          this.hoverRingPos3d = baryTetTo3D(mu);
        } else {
          this.hoverRingPos3d = null;
        }
      }
    } else {
      this.hoverRingPos3d = null;
    }
    this.emit('hoverChanged');
  }

  clearHover() {
    if (this.hoverPunctureIdx === -1 && this.hoverSegmentIdx === -1 &&
        this.hoverLambda === null) return;
    this.hoverPunctureIdx = -1;
    this.hoverSegmentIdx = -1;
    this.hoverLambda = null;
    this.hoverRingPos3d = null;
    this.emit('hoverChanged');
  }

  recompute() {
    if (this.dim === '2d') {
      this.recompute2D();
    } else {
      this.recompute3D();
    }
  }

  recompute3D() {
    const cc = classifyTetCase(this.V, this.W);
    this.Q = cc.Q_coeffs;
    this.P = cc.P_coeffs;
    this.qRoots = cc.Q_roots;
    this.qDegree = cc.qDegree;
    this.qDiscSign = cc.Q_disc_sign;
    this.punctures = cc.punctures;
    this.intervals = cc.intervals;
    this.category = cc.category;
    this.hasSR = cc.has_shared_root;
    this.hasB = cc.has_B;
    this.hasCvPos = cc.has_Cv_pos;
    this.hasCwPos = cc.has_Cw_pos;
    this.CvMu = cc.Cv_mu;
    this.CwMu = cc.Cw_mu;
    this.srLambda = cc.SR_lambda;
    this.srPos3d = cc.SR_pos3d;

    this.segments = buildSegments(cc, this.Q, this.P);

    if (cc.has_B) {
      const bub = sampleBubble(this.Q, this.P);
      this.bubble = bub ? bub.pts : null;
      if (bub && this.segments.length === 0) {
        this.segments = [{
          ptsList: [bub.pts],
          lamsList: [bub.lams],
          color: SEGMENT_COLORS[0],
          piEntry: -1, piExit: -1,
          lamEntry: null, lamExit: null,
          infinitySpanning: true,
          infPos3d: null, zeroPos3d: null,
        }];
      }
    } else {
      this.bubble = null;
    }
  }

  recompute2D() {
    let cc;
    try {
      cc = classifyTriCase2D(this.V2d, this.W2d);
    } catch (e) {
      console.error('classifyTriCase2D error:', e);
      return;
    }
    this.Q = cc.Q_coeffs;
    this.P = cc.P_coeffs;
    this.qRoots = cc.Q_roots;
    this.qDegree = cc.qDegree;
    this.qDiscSign = cc.Q_disc_sign;
    this.punctures = cc.punctures;
    this.intervals = cc.intervals;
    this.category = cc.category;
    this.hasSR = cc.has_shared_root;
    this.hasB = cc.has_B;
    this.hasCvPos = cc.has_Cv_pos;
    this.hasCwPos = cc.has_Cw_pos;
    this.CvMu = null;
    this.CwMu = null;
    this.srLambda = null;
    this.srPos3d = null;

    try {
      this.segments = buildSegments2D(cc, this.Q, this.P);
    } catch (e) {
      console.error('buildSegments2D error:', e);
      this.segments = [];
    }

    if (cc.has_B) {
      const bub = sampleBubble2D(this.Q, this.P);
      this.bubble = bub ? bub.pts : null;
      if (bub && this.segments.length === 0) {
        this.segments = [{
          ptsList: [bub.pts],
          lamsList: [bub.lams],
          color: SEGMENT_COLORS_2D[0],
          piEntry: -1, piExit: -1,
          lamEntry: null, lamExit: null,
          infinitySpanning: true,
          infPos3d: null, zeroPos3d: null,
        }];
      }
    } else {
      this.bubble = null;
    }
  }

  getDefaultV() { return DEFAULT_V.map(r => [...r]); }
  getDefaultW() { return DEFAULT_W.map(r => [...r]); }
}

/**
 * Compute the 3D position of the asymptotic point μ(∞) = P_k[3]/Q[3].
 * Returns null if Q[3]==0 or point is outside tet.
 */
function computeInfinityPos(Q, P) {
  if (Q[3] === 0) return null;
  const mu = P.map(pk => pk[3] / Q[3]);
  const muClip = mu.map(m => Math.max(0, m));
  const s = muClip.reduce((a, b) => a + b, 0);
  if (s < 1e-10) return null;
  const muNorm = muClip.map(m => m / s);
  return baryTetTo3D(muNorm);
}

/**
 * Compute the 3D position at λ=0: μ_k(0) = P_k[0]/Q[0].
 * Returns null if Q[0]==0 or point is outside tet.
 */
function computeZeroPos(Q, P) {
  if (Q[0] === 0) return null;
  const mu = P.map(pk => pk[0] / Q[0]);
  if (!mu.every(m => m >= -1e-6)) return null;
  const muClip = mu.map(m => Math.max(0, m));
  const s = muClip.reduce((a, b) => a + b, 0);
  if (s < 1e-10) return null;
  const muNorm = muClip.map(m => m / s);
  return baryTetTo3D(muNorm);
}

/**
 * Build visual segments from classify pairs, sampling curves for 3D rendering.
 * Infinity-spanning segments are directly sampled in two halves and joined
 * through the λ=∞ point.
 */
function buildSegments(cc, Q, P) {
  const { pairs, punctures, intervals } = cc;
  if (pairs.length === 0) return [];

  // Collect finite puncture λ values as sampling hints
  const lamHints = punctures
    .map(p => p.lambda)
    .filter(l => isFinite(l));

  // Compute critical positions
  const infPos3d = computeInfinityPos(Q, P);
  const zeroPos3d = computeZeroPos(Q, P);

  // Sample all intervals for non-cross segments
  const allSubsegs = [];
  for (const iv of intervals) {
    const segs = samplePVCurve(Q, P, iv.lb, iv.ub, iv.isInfinity, 200, lamHints);
    allSubsegs.push(...segs);
  }

  const segments = [];
  for (let idx = 0; idx < pairs.length; idx++) {
    const { pi_a, pi_b, contains_infinity } = pairs[idx];
    const p1 = punctures[pi_a];
    const p2 = punctures[pi_b];
    const color = SEGMENT_COLORS[idx % SEGMENT_COLORS.length];

    if (contains_infinity && infPos3d) {
      // ── Infinity-spanning segment ──
      const l1 = p1.lambda, l2 = p2.lambda;
      const l1Inf = !isFinite(l1), l2Inf = !isFinite(l2);

      if (l1Inf || l2Inf) {
        // One puncture at λ=∞: sample from finite puncture toward ∞
        const finiteLam = l1Inf ? l2 : l1;

        // Sample both directions from finiteLam and pick the half
        // that reaches toward infinity (has more inside points)
        const rSegs = samplePVCurve(Q, P, finiteLam, null, true, 200, [finiteLam]);
        const lSegs = samplePVCurve(Q, P, null, finiteLam, true, 200, [finiteLam]);
        const rSub = pickClosest(rSegs, finiteLam);
        const lSub = pickClosest(lSegs, finiteLam);

        // Pick the half with more points (the valid direction)
        const sub = (rSub && lSub)
          ? (rSub.pts.length >= lSub.pts.length ? rSub : lSub)
          : (rSub || lSub);

        // Infinity-spanning curves go through ∞, not through 0
        const segZeroPos = null;

        if (sub) {
          const pts = [...sub.pts];
          const lams = [...sub.lams];
          // Order: finiteLam at start, ∞ at end
          if (Math.abs(lams[lams.length - 1] - finiteLam) < Math.abs(lams[0] - finiteLam)) {
            pts.reverse();
            lams.reverse();
          }
          pts.push(infPos3d);
          lams.push(Infinity);
          segments.push({
            ptsList: [pts], lamsList: [lams], color,
            piEntry: pi_a, piExit: pi_b,
            lamEntry: p1.lambda, lamExit: p2.lambda,
            infinitySpanning: true, infPos3d, zeroPos3d: segZeroPos,
          });
        } else {
          // Fallback: just mark the infinity point
          segments.push({
            ptsList: [[infPos3d]], lamsList: [[Infinity]], color,
            piEntry: pi_a, piExit: pi_b,
            lamEntry: p1.lambda, lamExit: p2.lambda,
            infinitySpanning: true, infPos3d, zeroPos3d: null,
          });
        }
      } else {
        // Both punctures finite: sample each half and join through ∞
        const hiLam = Math.max(l1, l2);
        const loLam = Math.min(l1, l2);

        // Right half: from hiLam toward +∞
        const rightSegs = samplePVCurve(Q, P, hiLam, null, true, 200, [hiLam]);
        // Left half: from -∞ toward loLam
        const leftSegs = samplePVCurve(Q, P, null, loLam, true, 200, [loLam]);

        // Pick the sub-seg closest to each puncture, preferring correct direction
        const rightSub = pickForHalf(rightSegs, hiLam, 'right');
        const leftSub = pickForHalf(leftSegs, loLam, 'left');

        // Infinity-spanning curves go through ∞, not through 0
        const segZeroPos = null;

        if (rightSub && leftSub) {
          const rightAsc = rightSub.lamEntry <= rightSub.lamExit;
          const rightPts = rightAsc ? [...rightSub.pts] : [...rightSub.pts].reverse();
          const rightLams = rightAsc ? [...rightSub.lams] : [...rightSub.lams].reverse();
          const leftAsc = leftSub.lamEntry <= leftSub.lamExit;
          const leftPts = leftAsc ? [...leftSub.pts] : [...leftSub.pts].reverse();
          const leftLams = leftAsc ? [...leftSub.lams] : [...leftSub.lams].reverse();
          const joined = [...rightPts, infPos3d, ...leftPts];
          const joinedLams = [...rightLams, Infinity, ...leftLams];
          segments.push({
            ptsList: [joined], lamsList: [joinedLams], color,
            piEntry: pi_a, piExit: pi_b,
            lamEntry: p1.lambda, lamExit: p2.lambda,
            infinitySpanning: true, infPos3d, zeroPos3d: segZeroPos,
          });
        } else {
          const pts = [], lams = [];
          if (rightSub) { pts.push(rightSub.pts); lams.push(rightSub.lams); }
          if (leftSub) { pts.push(leftSub.pts); lams.push(leftSub.lams); }
          segments.push({
            ptsList: pts, lamsList: lams, color,
            piEntry: pi_a, piExit: pi_b,
            lamEntry: p1.lambda, lamExit: p2.lambda,
            infinitySpanning: true, infPos3d, zeroPos3d: segZeroPos,
          });
        }
      }
    } else {
      // ── Normal segment: match sub-segments by λ proximity ──
      const ptsList = [];
      const lamsList = [];
      for (const sub of allSubsegs) {
        const d = Math.min(
          lamDist(sub.lamEntry, p1.lambda), lamDist(sub.lamEntry, p2.lambda),
          lamDist(sub.lamExit, p1.lambda), lamDist(sub.lamExit, p2.lambda)
        );
        let bestPairDist = Infinity;
        for (let otherIdx = 0; otherIdx < pairs.length; otherIdx++) {
          if (otherIdx === idx) continue;
          const op1 = punctures[pairs[otherIdx].pi_a];
          const op2 = punctures[pairs[otherIdx].pi_b];
          bestPairDist = Math.min(bestPairDist,
            lamDist(sub.lamEntry, op1.lambda), lamDist(sub.lamEntry, op2.lambda),
            lamDist(sub.lamExit, op1.lambda), lamDist(sub.lamExit, op2.lambda)
          );
        }
        if (d <= bestPairDist) {
          ptsList.push(sub.pts);
          lamsList.push(sub.lams);
        }
      }

      const l1 = p1.lambda, l2 = p2.lambda;
      const spansZero = isFinite(l1) && isFinite(l2) &&
        ((l1 < 0 && l2 > 0) || (l1 > 0 && l2 < 0));
      const segZeroPos = spansZero ? zeroPos3d : null;

      segments.push({
        ptsList, lamsList, color,
        piEntry: pi_a, piExit: pi_b,
        lamEntry: p1.lambda, lamExit: p2.lambda,
        infinitySpanning: false, infPos3d: null, zeroPos3d: segZeroPos,
      });
    }
  }

  return segments;
}

/** Pick the sub-segment for one half of an infinity-spanning segment.
 *  Prefers sub-segments extending in the correct direction from the puncture.
 */
function pickForHalf(subsegs, punctureLam, direction) {
  let best = null, bestDist = Infinity;
  for (const sub of subsegs) {
    const midLam = (sub.lamEntry + sub.lamExit) / 2;
    const correctDir = direction === 'right'
      ? midLam > punctureLam
      : midLam < punctureLam;
    if (!correctDir) continue;
    const d = Math.min(
      Math.abs(sub.lamEntry - punctureLam),
      Math.abs(sub.lamExit - punctureLam)
    );
    if (d < bestDist) { best = sub; bestDist = d; }
  }
  // Fallback to closest regardless of direction
  return best || pickClosest(subsegs, punctureLam);
}

/** Pick the sub-segment closest to a given λ value. */
function pickClosest(subsegs, lam) {
  let best = null, bestDist = Infinity;
  for (const sub of subsegs) {
    const d = Math.min(
      Math.abs(sub.lamEntry - lam),
      Math.abs(sub.lamExit - lam)
    );
    if (d < bestDist) { best = sub; bestDist = d; }
  }
  return best;
}

/**
 * Build visual segments from 2D classify pairs.
 */
function buildSegments2D(cc, Q, P) {
  const { pairs, punctures, intervals } = cc;
  if (pairs.length === 0) return [];

  const lamHints = punctures
    .map(p => p.lambda)
    .filter(l => isFinite(l));

  const allSubsegs = [];
  for (const iv of intervals) {
    const segs = samplePVCurve2D(Q, P, iv.lb, iv.ub, iv.isInfinity, 200, lamHints);
    allSubsegs.push(...segs);
  }

  const segments = [];
  for (let idx = 0; idx < pairs.length; idx++) {
    const { pi_a, pi_b, contains_infinity } = pairs[idx];
    const p1 = punctures[pi_a];
    const p2 = punctures[pi_b];
    const color = SEGMENT_COLORS_2D[idx % SEGMENT_COLORS_2D.length];
    const isInfSpanning = contains_infinity || false;

    if (isInfSpanning) {
      // Infinity-spanning: sample both halves beyond each puncture
      const l1 = p1.lambda, l2 = p2.lambda;
      const hiLam = Math.max(l1, l2);
      const loLam = Math.min(l1, l2);
      const rightSegs = samplePVCurve2D(Q, P, hiLam, null, true, 200, [hiLam]);
      const leftSegs = samplePVCurve2D(Q, P, null, loLam, true, 200, [loLam]);
      const ptsList = [], lamsList = [];
      for (const sub of rightSegs) { ptsList.push(sub.pts); lamsList.push(sub.lams); }
      for (const sub of leftSegs) { ptsList.push(sub.pts); lamsList.push(sub.lams); }
      segments.push({
        ptsList, lamsList, color,
        piEntry: pi_a, piExit: pi_b,
        lamEntry: p1.lambda, lamExit: p2.lambda,
        infinitySpanning: true,
        infPos3d: null, zeroPos3d: null,
      });
    } else {
      // Normal segment: match sub-segments by λ proximity
      const ptsList = [];
      const lamsList = [];
      for (const sub of allSubsegs) {
        const d = Math.min(
          lamDist(sub.lamEntry, p1.lambda), lamDist(sub.lamEntry, p2.lambda),
          lamDist(sub.lamExit, p1.lambda), lamDist(sub.lamExit, p2.lambda)
        );
        let bestPairDist = Infinity;
        for (let otherIdx = 0; otherIdx < pairs.length; otherIdx++) {
          if (otherIdx === idx) continue;
          const op1 = punctures[pairs[otherIdx].pi_a];
          const op2 = punctures[pairs[otherIdx].pi_b];
          bestPairDist = Math.min(bestPairDist,
            lamDist(sub.lamEntry, op1.lambda), lamDist(sub.lamEntry, op2.lambda),
            lamDist(sub.lamExit, op1.lambda), lamDist(sub.lamExit, op2.lambda)
          );
        }
        if (d <= bestPairDist) {
          ptsList.push(sub.pts);
          lamsList.push(sub.lams);
        }
      }
      segments.push({
        ptsList, lamsList, color,
        piEntry: pi_a, piExit: pi_b,
        lamEntry: p1.lambda, lamExit: p2.lambda,
        infinitySpanning: false,
        infPos3d: null, zeroPos3d: null,
      });
    }
  }
  return segments;
}

function lamDist(a, b) {
  if (a === null && b === null) return 0;
  if (a === null || b === null) return 1e15;
  if (!isFinite(a) && !isFinite(b)) return 0;
  if (!isFinite(a) || !isFinite(b)) return 1e15;
  return Math.abs(a - b);
}

export const state = new State();
