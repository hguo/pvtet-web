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

  // Build aligned polynomial table: Q and P₀–P₃ with coefficients lined up
  const polyTableLatex = buildAlignedPolyLatex(
    state.Q, state.P, 'Q', ['P_0','P_1','P_2','P_3']
  );

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

  // Polynomials (aligned) + discriminant
  html += '<div>';
  html += `<div class="poly-line">${renderKatex(polyTableLatex, true)}</div>`;
  html += `<div class="poly-line" style="margin-top:4px">${renderKatex(discStr)}</div>`;
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

/**
 * Build a KaTeX aligned array with Q and P₀–P₃ polynomials,
 * coefficients lined up by power of λ.
 */
function buildAlignedPolyLatex(Q, P, qName, pNames) {
  const allPolys = [{ name: qName, coeffs: Q }, ...P.map((c, i) => ({ name: pNames[i], coeffs: c }))];

  // Find max degree with nonzero coefficient across all polys
  let maxDeg = 0;
  for (const { coeffs } of allPolys) {
    for (let d = coeffs.length - 1; d >= 0; d--) {
      if (Math.abs(coeffs[d]) > 1e-30) { maxDeg = Math.max(maxDeg, d); break; }
    }
  }

  const rows = [];
  for (const { name, coeffs } of allPolys) {
    const cols = [];
    let firstTerm = true;
    for (let d = maxDeg; d >= 0; d--) {
      const c = Math.round(coeffs[d]);
      const varPart = d === 0 ? '' : d === 1 ? '\\lambda' : `\\lambda^{${d}}`;

      if (c === 0) {
        // Empty column to maintain alignment
        cols.push('');
        cols.push('');
      } else if (firstTerm) {
        cols.push(''); // no sign column for first term
        if (Math.abs(c) === 1 && d > 0) {
          cols.push(c < 0 ? `-${varPart}` : `${varPart}`);
        } else {
          cols.push(`${c}${varPart}`);
        }
        firstTerm = false;
      } else {
        const sign = c > 0 ? '+' : '-';
        const absC = Math.abs(c);
        cols.push(sign);
        if (absC === 1 && d > 0) {
          cols.push(varPart);
        } else {
          cols.push(`${absC}${varPart}`);
        }
      }
    }
    rows.push(`${name}(\\lambda) &=& ${cols.join(' & ')}`);
  }

  // Column spec: name = col, then for each degree: sign col (c) + value col (r)
  const nCols = (maxDeg + 1) * 2;
  const colSpec = 'rl' + 'rl'.repeat(maxDeg + 1);

  return `\\begin{array}{${colSpec}} ${rows.join(' \\\\ ')} \\end{array}`;
}

function renderKatex(latex, display = false) {
  try {
    return katex.renderToString(latex, { throwOnError: false, displayMode: display });
  } catch {
    return `<code>${latex}</code>`;
  }
}
