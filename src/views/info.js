/**
 * Info panel: KaTeX polynomial display + classification info + segment info.
 */
import katex from 'katex';
import { state } from '../state.js';
import { polyToLatex } from '../math/polynomials.js';

export function initInfo() {
  state.on('dataChanged', updateInfo);
  updateInfo();
}

function updateInfo() {
  const panel = document.getElementById('panel-info');

  // Build polynomial display
  const qLatex = polyToLatex(state.Q, 'Q');
  const pLatex = state.P.map((p, i) => polyToLatex(p, `P_{${i}}`));

  // Discriminant info
  let discStr;
  if (state.qDegree < 3) {
    discStr = `\\text{deg}(Q)=${state.qDegree}\\;(${state.qRoots.length}\\;\\text{root${state.qRoots.length !== 1 ? 's' : ''}})`;
  } else {
    const labels = {
      1: '\\Delta_Q > 0\\;(3\\;\\text{roots})',
      '-1': '\\Delta_Q < 0\\;(1\\;\\text{root})',
      0: '\\Delta_Q = 0',
    };
    discStr = labels[state.qDiscSign] || `\\Delta_Q\\;\\text{sign}=${state.qDiscSign}`;
  }

  // Segment info
  const segLines = [];
  for (let i = 0; i < state.segments.length; i++) {
    const seg = state.segments[i];
    const l1 = seg.lamEntry;
    const l2 = seg.lamExit;
    const l1s = (l1 !== null && isFinite(l1)) ? l1.toFixed(3) : '\u221E';
    const l2s = (l2 !== null && isFinite(l2)) ? l2.toFixed(3) : '\u221E';

    let desc;
    if (l1 === null && l2 === null) {
      desc = `S${i + 1}: (bubble)`;
    } else if (seg.infinitySpanning) {
      desc = `S${i + 1}: (${l1s}, +\u221E) \u222A (-\u221E, ${l2s})`;
    } else {
      const lo = Math.min(l1, l2).toFixed(3);
      const hi = Math.max(l1, l2).toFixed(3);
      desc = `S${i + 1}: (${lo}, ${hi})`;
    }

    const dot = `<span style="color:${seg.color}; font-size:16px">\u25CF</span>`;
    segLines.push(`${dot} ${desc}`);
  }

  // Render
  let html = '<div style="display:flex; gap:16px; align-items:flex-start; flex-wrap:wrap">';

  // Polynomials
  html += '<div>';
  html += `<div class="poly-line">${renderKatex(qLatex)}  &nbsp; ${renderKatex(discStr)}</div>`;
  for (const l of pLatex) {
    html += `<div class="poly-line" style="font-size:12px">${renderKatex(l)}</div>`;
  }
  html += '</div>';

  // Classification + segments
  html += '<div class="seg-info">';
  html += `<strong>Category:</strong> <code>${state.category}</code><br>`;
  html += `#Punctures: ${state.punctures.length} &nbsp; #Segments: ${state.segments.length}`;
  if (state.hasB) html += ' &nbsp; (Bubble)';
  if (state.hasSR) html += ' &nbsp; (SR)';
  if (segLines.length > 0) html += '<br>' + segLines.join(' &nbsp; ');
  html += '</div>';

  html += '</div>';

  panel.innerHTML = html;
}

function renderKatex(latex) {
  try {
    return katex.renderToString(latex, { throwOnError: false, displayMode: false });
  } catch {
    return `<code>${latex}</code>`;
  }
}
