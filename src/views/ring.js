/**
 * D3.js SVG λ-ring view.
 * Ring circle, segment bands, Q-root markers, puncture ticks, λ labels.
 */
import * as d3 from 'd3';
import { state } from '../state.js';
import { SEGMENT_COLORS } from '../math/curves.js';

let svg, g; // D3 selections
let width, height;
const R_RING = 1.0;
const R_INNER = 0.85;
const R_OUTER = 0.98;

// User-adjustable scale multiplier (mouse scroll)
let scaleMultiplier = 1.0;
let showQRoots = true;

export function initRing() {
  const container = document.getElementById('panel-ring');
  svg = d3.select('#ring-svg');

  // Responsive sizing
  const ro = new ResizeObserver(() => updateSize());
  ro.observe(container);

  g = svg.append('g');

  state.on('dataChanged', () => { scaleMultiplier = 1.0; redraw(); });
  state.on('hoverChanged', updateHoverHighlights);

  // Scroll to adjust λ scale
  svg.on('wheel', (event) => {
    event.preventDefault();
    const factor = event.deltaY > 0 ? 1.15 : 1 / 1.15;
    scaleMultiplier *= factor;
    scaleMultiplier = Math.max(0.01, Math.min(100, scaleMultiplier));
    redraw();
    updateHoverHighlights();
  });

  // Q-roots visibility checkbox
  const chk = document.getElementById('chk-show-qroots');
  showQRoots = chk.checked;
  chk.addEventListener('change', () => {
    showQRoots = chk.checked;
    redraw();
  });

  updateSize();
  redraw();
}

function updateSize() {
  const container = document.getElementById('panel-ring');
  width = container.clientWidth;
  height = container.clientHeight;
  svg.attr('viewBox', '-2.2 -1.8 4.4 4.0');
  svg.attr('width', width).attr('height', height);
}

function computeScale() {
  const allLams = [];
  for (const p of state.punctures) {
    if (p.lambda !== null && isFinite(p.lambda)) allLams.push(Math.abs(p.lambda));
  }
  for (const r of state.qRoots) {
    allLams.push(Math.abs(r));
  }
  const absVals = allLams.filter(v => v > 1e-15);
  let scale = absVals.length > 0 ? Math.max(...absVals) : 1.0;
  return Math.max(scale, 0.5) * scaleMultiplier;
}

function lambdaToAngle(lam, scale) {
  // Maps λ=0 → bottom (π/2), λ→±∞ → top (-π/2)
  return 2 * Math.atan(lam / scale) + Math.PI / 2;
}

function angleToXY(angle, radius) {
  return [radius * Math.cos(angle), radius * Math.sin(angle)];
}

function redraw() {
  g.selectAll('*').remove();

  const scale = computeScale();

  // Ring circle
  g.append('circle')
    .attr('cx', 0).attr('cy', 0).attr('r', R_RING)
    .attr('fill', 'none').attr('stroke', '#333').attr('stroke-width', 0.02);

  // Infinity label at top
  g.append('text')
    .attr('x', 0).attr('y', -(R_RING + 0.10))
    .attr('text-anchor', 'middle').attr('dominant-baseline', 'auto')
    .attr('font-size', '0.14px').attr('font-weight', 'bold').attr('fill', '#333')
    .text('\u221E');

  // Category label above the ring
  g.append('text')
    .attr('x', 0).attr('y', -(R_RING + 0.28))
    .attr('text-anchor', 'middle').attr('dominant-baseline', 'auto')
    .attr('font-size', '0.10px').attr('fill', '#555').attr('font-family', 'monospace')
    .text(state.category);

  // Zero tick at bottom
  const a0 = lambdaToAngle(0, scale);
  const [x0, y0] = angleToXY(a0, R_RING);
  g.append('line')
    .attr('x1', x0 * 0.94).attr('y1', y0 * 0.94)
    .attr('x2', x0 * 1.06).attr('y2', y0 * 1.06)
    .attr('stroke', '#666').attr('stroke-width', 0.02);
  g.append('text')
    .attr('x', x0).attr('y', y0 + 0.14)
    .attr('text-anchor', 'middle').attr('dominant-baseline', 'hanging')
    .attr('font-size', '0.10px').attr('fill', '#666')
    .text('0');

  // Segment bands
  drawSegmentBands(g, scale);

  // Q-root markers
  if (showQRoots) drawQRoots(g, scale);

  // SR marker
  if (state.hasSR && state.srLambda !== null) {
    drawSRMarker(g, scale);
  }

  // Puncture ticks and labels
  drawPunctureTicks(g, scale);

  // Hover highlight layer (empty, filled by updateHoverHighlights)
  g.append('g').attr('class', 'hover-layer');

  // Interactive overlay: invisible circle to capture mouse events for ring hover
  // Large transparent overlay for hover interaction across the full ring area
  const overlay = g.append('circle')
    .attr('cx', 0).attr('cy', 0).attr('r', 1.8)
    .attr('fill', 'transparent').attr('stroke', 'none')
    .style('cursor', 'crosshair');

  overlay.on('mousemove', (event) => {
    const [mx, my] = d3.pointer(event, g.node());
    const dist = Math.sqrt(mx * mx + my * my);
    if (dist < 0.05 || dist > 2.0) {
      state.clearHover();
      return;
    }
    // Convert mouse position to angle, then to lambda
    let angle = Math.atan2(my, mx);
    // Unwrap: forward mapping produces [-π/2, 3π/2), atan2 returns (-π, π]
    // Angles < -π/2 (left-top region) need +2π to match
    if (angle < -Math.PI / 2) angle += 2 * Math.PI;
    // Invert: angle = 2*atan(lam/scale) + pi/2
    // lam = scale * tan((angle - pi/2) / 2)
    const halfAngle = (angle - Math.PI / 2) / 2;
    if (Math.abs(halfAngle) > 1.55) {
      // Very close to infinity — use Infinity
      state.setHoverLambda(Infinity);
      return;
    }
    const lam = scale * Math.tan(halfAngle);
    state.setHoverLambda(lam);
  });

  overlay.on('mouseleave', () => {
    state.clearHover();
  });
}

function drawSegmentBands(g, scale) {
  for (const seg of state.segments) {
    let a1 = seg.lamEntry !== null ? lambdaToAngle(seg.lamEntry, scale) : Math.PI / 2;
    let a2 = seg.lamExit !== null ? lambdaToAngle(seg.lamExit, scale) : Math.PI / 2;

    if (seg.infinitySpanning) {
      const aHi = Math.max(a1, a2);
      const aLo = Math.min(a1, a2);
      drawBand(g, aHi, aLo + 2 * Math.PI, seg.color);
    } else {
      if (a1 > a2) [a1, a2] = [a2, a1];
      drawBand(g, a1, a2, seg.color);
    }
  }
}

function drawBand(g, aStart, aEnd, color) {
  const nPts = 80;
  const arcTh = [];
  for (let i = 0; i < nPts; i++) {
    arcTh.push(aStart + (aEnd - aStart) * i / (nPts - 1));
  }

  const outerPts = arcTh.map(a => angleToXY(a, R_OUTER));
  const innerPts = arcTh.map(a => angleToXY(a, R_INNER)).reverse();
  const allPts = [...outerPts, ...innerPts];

  const pathData = 'M' + allPts.map(p => `${p[0]},${p[1]}`).join('L') + 'Z';

  g.append('path')
    .attr('d', pathData)
    .attr('fill', color).attr('fill-opacity', 0.45)
    .attr('stroke', color).attr('stroke-width', 0.005)
    .attr('class', 'seg-band')
    .on('mouseenter', () => {
      const idx = state.segments.indexOf(
        state.segments.find(s => s.color === color)
      );
      state.setHoverSegment(idx);
    })
    .on('mouseleave', () => state.clearHover());
}

function drawQRoots(g, scale) {
  const qAngles = state.qRoots.map(r => lambdaToAngle(r, scale));
  const labelRadii = state.qRoots.map(() => 0.60);

  // Stagger close labels
  const sorted = state.qRoots.map((_, i) => i).sort((a, b) => qAngles[a] - qAngles[b]);
  for (let j = 1; j < sorted.length; j++) {
    const ip = sorted[j - 1], ic = sorted[j];
    if (Math.abs(qAngles[ic] - qAngles[ip]) < 0.30) {
      labelRadii[ic] = labelRadii[ip] >= 0.55 ? 0.42 : 0.60;
    }
  }

  for (let i = 0; i < state.qRoots.length; i++) {
    const r = state.qRoots[i];
    const a = qAngles[i];
    const [rx, ry] = angleToXY(a, R_RING);

    // Hollow circle on ring
    g.append('circle')
      .attr('cx', rx).attr('cy', ry).attr('r', 0.05)
      .attr('fill', 'white').attr('stroke', 'black').attr('stroke-width', 0.015);

    // Leader line + label inside
    const lr = labelRadii[i];
    const [lx, ly] = angleToXY(a, lr);
    const [mx, my] = angleToXY(a, R_RING - 0.06);
    g.append('line')
      .attr('x1', mx).attr('y1', my).attr('x2', lx).attr('y2', ly)
      .attr('stroke', '#888').attr('stroke-width', 0.006);
    g.append('text')
      .attr('x', lx).attr('y', ly)
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
      .attr('font-size', '0.08px').attr('fill', '#333')
      .text(r.toFixed(2));
  }
}

function drawSRMarker(g, scale) {
  const a = lambdaToAngle(state.srLambda, scale);
  const [rx, ry] = angleToXY(a, R_RING);

  // Diamond marker on ring
  const size = 0.06;
  const diamond = `M${rx},${ry - size} L${rx + size * 0.6},${ry} L${rx},${ry + size} L${rx - size * 0.6},${ry} Z`;
  g.append('path')
    .attr('d', diamond)
    .attr('fill', '#ff00ff').attr('fill-opacity', 0.7)
    .attr('stroke', 'black').attr('stroke-width', 0.01);

  // Label inside ring
  const [lx, ly] = angleToXY(a, 0.65);
  const [mx, my] = angleToXY(a, R_RING - 0.06);
  g.append('line')
    .attr('x1', mx).attr('y1', my).attr('x2', lx).attr('y2', ly)
    .attr('stroke', '#ff00ff').attr('stroke-width', 0.008).attr('opacity', 0.6);
  g.append('text')
    .attr('x', lx).attr('y', ly)
    .attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
    .attr('font-size', '0.09px').attr('fill', '#ff00ff').attr('font-weight', 'bold')
    .text('SR');
}

function drawPunctureTicks(g, scale) {
  // Build puncture → color map
  const puncColor = {};
  for (const seg of state.segments) {
    if (!(seg.piEntry in puncColor)) puncColor[seg.piEntry] = seg.color;
    if (!(seg.piExit in puncColor)) puncColor[seg.piExit] = seg.color;
  }

  // Compute angles
  const puncAngles = state.punctures.map(p =>
    p.lambda !== null ? lambdaToAngle(p.lambda, scale) : Math.PI / 2
  );

  // Stagger label radii
  const labelRadii = state.punctures.map(() => R_RING + 0.35);
  const sortedIdx = state.punctures.map((_, i) => i).sort((a, b) => puncAngles[a] - puncAngles[b]);
  for (let j = 1; j < sortedIdx.length; j++) {
    const iPrev = sortedIdx[j - 1], iCurr = sortedIdx[j];
    if (Math.abs(puncAngles[iCurr] - puncAngles[iPrev]) < 0.25) {
      labelRadii[iCurr] = labelRadii[iPrev] < R_RING + 0.50
        ? R_RING + 0.55
        : R_RING + 0.35;
    }
  }

  for (let i = 0; i < state.punctures.length; i++) {
    const p = state.punctures[i];
    const color = puncColor[i] || 'black';
    const a = puncAngles[i];
    const lamStr = p.lambda !== null ? p.lambda.toFixed(2) : '\u221E';

    // Tick mark on ring
    const [x1, y1] = angleToXY(a, 0.88);
    const [x2, y2] = angleToXY(a, 1.07);
    g.append('line')
      .attr('x1', x1).attr('y1', y1).attr('x2', x2).attr('y2', y2)
      .attr('stroke', color).attr('stroke-width', 0.025)
      .attr('class', `punc-tick punc-tick-${i}`)
      .on('mouseenter', () => state.setHoverPuncture(i))
      .on('mouseleave', () => state.clearHover());

    // Lambda label outside
    const lr = labelRadii[i];
    const [lx, ly] = angleToXY(a, lr);
    const [mx, my] = angleToXY(a, R_RING + 0.08);
    g.append('line')
      .attr('x1', mx).attr('y1', my).attr('x2', lx).attr('y2', ly)
      .attr('stroke', color).attr('stroke-width', 0.006).attr('opacity', 0.5);
    g.append('text')
      .attr('x', lx).attr('y', ly)
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
      .attr('font-size', '0.08px').attr('fill', color)
      .attr('class', `punc-label punc-label-${i}`)
      .text(lamStr);
  }
}

function updateHoverHighlights() {
  const hover = g.select('.hover-layer');
  hover.selectAll('*').remove();

  const scale = computeScale();

  // Highlight hovered puncture tick
  if (state.hoverPunctureIdx >= 0 && state.hoverPunctureIdx < state.punctures.length) {
    const p = state.punctures[state.hoverPunctureIdx];
    const a = p.lambda !== null ? lambdaToAngle(p.lambda, scale) : Math.PI / 2;
    const [x1, y1] = angleToXY(a, 0.82);
    const [x2, y2] = angleToXY(a, 1.13);
    hover.append('line')
      .attr('x1', x1).attr('y1', y1).attr('x2', x2).attr('y2', y2)
      .attr('stroke', '#ffcc00').attr('stroke-width', 0.04).attr('opacity', 0.8);
  }

  // Highlight hovered segment band
  if (state.hoverSegmentIdx >= 0 && state.hoverSegmentIdx < state.segments.length) {
    const seg = state.segments[state.hoverSegmentIdx];
    let a1 = seg.lamEntry !== null ? lambdaToAngle(seg.lamEntry, scale) : Math.PI / 2;
    let a2 = seg.lamExit !== null ? lambdaToAngle(seg.lamExit, scale) : Math.PI / 2;

    if (seg.infinitySpanning) {
      drawHighlightBand(hover, Math.max(a1, a2), Math.min(a1, a2) + 2 * Math.PI, seg.color);
    } else {
      if (a1 > a2) [a1, a2] = [a2, a1];
      drawHighlightBand(hover, a1, a2, seg.color);
    }
  }

  // Highlight ring hover lambda position
  if (state.hoverLambda !== null) {
    const a = lambdaToAngle(state.hoverLambda, scale);
    const [rx, ry] = angleToXY(a, R_RING);
    hover.append('circle')
      .attr('cx', rx).attr('cy', ry).attr('r', 0.04)
      .attr('fill', '#ff8800').attr('fill-opacity', 0.8)
      .attr('stroke', 'black').attr('stroke-width', 0.01);

    // Show lambda value
    const [lx, ly] = angleToXY(a, R_RING + 0.18);
    hover.append('text')
      .attr('x', lx).attr('y', ly)
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
      .attr('font-size', '0.07px').attr('fill', '#ff8800').attr('font-weight', 'bold')
      .text(`\u03BB=${state.hoverLambda.toFixed(2)}`);
  }
}

function drawHighlightBand(parent, aStart, aEnd, color) {
  const nPts = 40;
  const arcTh = [];
  for (let i = 0; i < nPts; i++) {
    arcTh.push(aStart + (aEnd - aStart) * i / (nPts - 1));
  }
  const outerPts = arcTh.map(a => angleToXY(a, R_OUTER + 0.03));
  const innerPts = arcTh.map(a => angleToXY(a, R_INNER - 0.03)).reverse();
  const allPts = [...outerPts, ...innerPts];
  const pathData = 'M' + allPts.map(p => `${p[0]},${p[1]}`).join('L') + 'Z';

  parent.append('path')
    .attr('d', pathData)
    .attr('fill', 'none')
    .attr('stroke', color).attr('stroke-width', 0.03)
    .attr('stroke-opacity', 0.9);
}
