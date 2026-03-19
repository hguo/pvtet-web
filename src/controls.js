/**
 * V/W matrix input grids, Randomize/Reset buttons, Seed/Range controls,
 * Gallery, Animation, Search, Export, Keyboard shortcuts, URL hash state.
 */
import { state } from './state.js';
import { randomVW, randomVW2D } from './math/random.js';
import { classifyTetCase } from './math/classify.js';
import { classifyTriCase2D } from './math/classify2d.js';

let vInputs = []; // vInputs[i][j] = <input> for V[i][j]
let wInputs = [];

// Animation state
let animTimer = null;
let animSeed = 0;

export function initControls() {
  rebuildMatrixGrids();

  document.getElementById('btn-randomize').addEventListener('click', doRandomize);
  document.getElementById('btn-reset').addEventListener('click', doReset);
  const seedInput = document.getElementById('input-seed');
  seedInput.addEventListener('change', doLoadSeed);
  seedInput.addEventListener('input', doLoadSeed);
  seedInput.addEventListener('keyup', doLoadSeed);
  document.getElementById('select-range').addEventListener('change', doLoadSeed);

  // Gallery
  document.getElementById('select-gallery').addEventListener('change', (e) => {
    const val = e.target.value;
    if (val === '') return;
    const seed = parseInt(val);
    document.getElementById('input-seed').value = seed;
    doLoadSeed();
  });

  // Animation
  document.getElementById('btn-animate').addEventListener('click', toggleAnimation);

  // Search
  document.getElementById('input-search').addEventListener('input', doSearch);

  // Export
  document.getElementById('btn-export').addEventListener('click', doExport);

  // Dark mode
  document.getElementById('btn-dark-mode').addEventListener('click', toggleDarkMode);

  // Mode toggle
  const btnMode = document.getElementById('btn-mode-toggle');
  if (btnMode) btnMode.addEventListener('click', toggleMode);

  // Keyboard shortcuts
  document.addEventListener('keydown', onKeyDown);

  // Sync inputs from state
  state.on('dataChanged', syncInputsFromState);
  state.on('modeChanged', onModeChanged);
  syncInputsFromState();

  // Load from URL hash if present
  loadFromHash();
  window.addEventListener('hashchange', loadFromHash);

  // Update hash when data changes
  state.on('dataChanged', updateHash);

  updateModeUI();
  updateGalleryOptions();
}

function toggleMode() {
  const newDim = state.dim === '2d' ? '3d' : '2d';
  state.setDim(newDim);
}

function onModeChanged() {
  updateModeUI();
  rebuildMatrixGrids();
  updateGalleryOptions();
  // Load seed 0 for the new mode to get valid data
  doLoadSeed();
  syncInputsFromState();
}

function updateModeUI() {
  const btn = document.getElementById('btn-mode-toggle');
  if (btn) btn.textContent = state.dim === '2d' ? 'Switch to 3D' : 'Switch to 2D';

  document.body.classList.toggle('mode-2d', state.dim === '2d');
  document.body.classList.toggle('mode-3d', state.dim === '3d');

  // Update field labels
  const vLabel = document.querySelector('#panel-controls h3 .field-label');
  if (vLabel) vLabel.textContent = state.dim === '2d' ? '(3\u00d72)' : '(4\u00d73)';
  const wParent = document.querySelectorAll('#panel-controls h3');
  for (const h of wParent) {
    const fl = h.querySelector('.field-label');
    if (fl && h.textContent.startsWith('W'))
      fl.textContent = state.dim === '2d' ? '(3\u00d72)' : '(4\u00d73)';
  }
}

function updateGalleryOptions() {
  const sel = document.getElementById('select-gallery');
  sel.innerHTML = '<option value="">-- Select a case --</option>';
  if (state.dim === '2d') {
    // Curated from ftk2 pv_tri_cases_2d/curated_v3.jsonl
    const opts = [
      [78187, 'T0_Q0'],
      [8751, 'T0_Q2-'],
      [7000, 'T0_Q2-_Cv_Cw_B'],
      [8744, 'T0_Q2+'],
      [545, 'T0_Q2+_SR'],
      [25644, 'T0_Q2+_ISR_Cv_Cw'],
      [31432, 'T0_Q2+_TN'],
      [49064, 'T0_Qz_Cv1_D11'],
      [1892, 'T2_Q1'],
      [7266, 'T2_Q2+'],
      [8764, 'T2_Q2-'],
      [11620, 'T2_(1,1)_Q2+_SR'],
      [11899, 'T2_(1,1)_Q2+_Cv_Cw'],
      [14219, 'T2_(1,1)_Q2+_Cv_Cw1'],
      [53518, 'T2_(1,1)_Q1_Cv_Cw'],
      [32116, 'T2_(1,1)_Q2+_SR_Cv_Cw'],
      [13395, 'T1_Q2+_Cv0'],
      [12388, 'T2_Q1'],
      [303, 'T2_(1,1)_Q2+_Cw0_D00'],
      [73915, 'T6_(2,4)_Q2+'],
    ];
    for (const [seed, cat] of opts)
      sel.innerHTML += `<option value="${seed}">${seed}: ${cat}</option>`;
  } else {
    // Curated from ftk2 paper_selection_v21.jsonl
    const opts = [
      [0, 'T0_Q2-_Cv_B'],
      [1, 'T2_Q3+_SR (NISR)'],
      [3618, 'T0_Q3+'],
      [991, 'T0_Q3-'],
      [993, 'T2_Q3+'],
      [990, 'T2_Q3-'],
      [2360, 'T2_Q2'],
      [5, 'T2_(1,1)_Q3+_Cw'],
      [17657, 'T2_(1,1)_Q3+_Cv_Cw'],
      [4, 'T4_(2,2)_Q3+_Cv'],
      [9579, 'T4_(2,2)_Q3+'],
      [13109, 'T4_(2,2)_Q3-'],
      [33, 'T4_(1,1,2)_Q3+_Cw'],
      [101980, 'T4_(1,1,2)_Q3+_SR'],
      [5976, 'T4_(1,3)_Q3+_Cw'],
      [1611, 'T4_Q3+'],
      [3617, 'T4_Q3-'],
      [1704, 'T4_Q3-_Cw0_D00'],
      [11021, 'T6_(2,2,2)_Q3+'],
      [4988, 'T6_(2,4)_Q3+'],
      [2397, 'T6_Q3+'],
      [10322, 'T6_Q3-'],
      [292, 'T6_(1,2,3)_Q3+_Cw'],
      [10553, 'T8_(2,2,4)_Q3+'],
      [25710, 'T8_(2,6)_Q3-'],
      [30810, 'T8_Q3-'],
      [4325, 'T0_Q3-_D00'],
      [65893, 'T2_Q3+_D00'],
      [23330, 'T6_(2,2,2)_Q3+_D00'],
    ];
    for (const [seed, cat] of opts)
      sel.innerHTML += `<option value="${seed}">${seed}: ${cat}</option>`;
  }
}

function rebuildMatrixGrids() {
  const nRows = state.dim === '2d' ? 3 : 4;
  const nCols = state.dim === '2d' ? 2 : 3;
  buildMatrixGrid('v-grid', 'V', true, nRows, nCols);
  buildMatrixGrid('w-grid', 'W', false, nRows, nCols);
}

function buildMatrixGrid(containerId, fieldName, isV, nRows = 4, nCols = 3) {
  const container = document.getElementById(containerId);
  container.innerHTML = '';
  container.style.gridTemplateColumns = `auto repeat(${nCols}, 1fr)`;
  const inputs = [];

  for (let i = 0; i < nRows; i++) {
    const label = document.createElement('span');
    label.className = 'row-label';
    label.textContent = `${fieldName.toLowerCase()}${String.fromCharCode(0x2080 + i)}`;
    container.appendChild(label);

    const row = [];
    for (let j = 0; j < nCols; j++) {
      const inp = document.createElement('input');
      inp.type = 'number';
      inp.step = 'any';
      inp.dataset.row = i;
      inp.dataset.col = j;
      inp.addEventListener('change', () => onMatrixInput(isV));
      container.appendChild(inp);
      row.push(inp);
    }
    inputs.push(row);
  }

  if (isV) vInputs = inputs;
  else wInputs = inputs;
}

function onMatrixInput(isV) {
  if (state.dim === '2d') {
    const V = state.V2d.map(r => [...r]);
    const W = state.W2d.map(r => [...r]);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 2; j++) {
        if (isV) V[i][j] = parseFloat(vInputs[i][j].value) || 0;
        else W[i][j] = parseFloat(wInputs[i][j].value) || 0;
      }
    }
    state.setVW(V, W);
  } else {
    const V = state.V.map(r => [...r]);
    const W = state.W.map(r => [...r]);
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 3; j++) {
        if (isV) V[i][j] = parseFloat(vInputs[i][j].value) || 0;
        else W[i][j] = parseFloat(wInputs[i][j].value) || 0;
      }
    }
    state.setVW(V, W);
  }
}

function syncInputsFromState() {
  if (state.dim === '2d') {
    if (vInputs.length < 3 || vInputs[0].length < 2) return;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 2; j++) {
        vInputs[i][j].value = state.V2d[i][j];
        wInputs[i][j].value = state.W2d[i][j];
      }
    }
  } else {
    if (vInputs.length < 4 || vInputs[0].length < 3) return;
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 3; j++) {
        vInputs[i][j].value = state.V[i][j];
        wInputs[i][j].value = state.W[i][j];
      }
    }
  }
}

function doLoadSeed() {
  const seedStr = document.getElementById('input-seed').value;
  if (seedStr === '') return;
  const seed = parseInt(seedStr);
  if (isNaN(seed) || seed < 0) return;
  const range = parseInt(document.getElementById('select-range').value) || 20;
  if (state.dim === '2d') {
    const { V, W } = randomVW2D(seed, range);
    state.setVW(V, W);
  } else {
    const { V, W } = randomVW(seed, range);
    state.setVW(V, W);
  }
}

function doRandomize() {
  doLoadSeed();
  const seedInput = document.getElementById('input-seed');
  seedInput.value = (parseInt(seedInput.value) || 0) + 1;
}

function doReset() {
  document.getElementById('input-seed').value = '0';
  doLoadSeed();
}

// ── Animation ──

function toggleAnimation() {
  const btn = document.getElementById('btn-animate');
  if (animTimer) {
    clearInterval(animTimer);
    animTimer = null;
    btn.textContent = 'Play';
  } else {
    animSeed = parseInt(document.getElementById('input-seed').value) || 0;
    btn.textContent = 'Stop';
    stepAnimation();
    animTimer = setInterval(stepAnimation, getAnimDelay());
  }
}

function stepAnimation() {
  animSeed++;
  document.getElementById('input-seed').value = animSeed;
  doLoadSeed();
}

function getAnimDelay() {
  const speed = parseInt(document.getElementById('input-anim-speed').value) || 30;
  // speed 1 → 2000ms, speed 100 → 50ms
  return Math.round(2000 / (1 + speed * 0.39));
}

// ── Search ──

function doSearch() {
  const query = document.getElementById('input-search').value.trim().toUpperCase();
  const results = document.getElementById('search-results');
  if (!query) {
    results.innerHTML = '';
    return;
  }

  const range = parseInt(document.getElementById('select-range').value) || 20;
  const matches = [];
  const maxScan = 5000;

  for (let seed = 0; seed < maxScan && matches.length < 20; seed++) {
    try {
      let cc;
      if (state.dim === '2d') {
        const { V, W } = randomVW2D(seed, range);
        cc = classifyTriCase2D(V, W);
      } else {
        const { V, W } = randomVW(seed, range);
        cc = classifyTetCase(V, W);
      }
      if (cc.category.toUpperCase().includes(query)) {
        matches.push({ seed, category: cc.category });
      }
    } catch {
      // skip errors
    }
  }

  if (matches.length === 0) {
    results.innerHTML = `<em>No matches in seeds 0-${maxScan - 1}</em>`;
  } else {
    results.innerHTML = matches.map(m =>
      `<a href="#" data-seed="${m.seed}" style="color:inherit">${m.seed}</a>: ${m.category}`
    ).join('<br>');

    // Click to load
    results.querySelectorAll('a').forEach(a => {
      a.addEventListener('click', (e) => {
        e.preventDefault();
        document.getElementById('input-seed').value = a.dataset.seed;
        doLoadSeed();
      });
    });
  }
}

// ── Export ──

function doExport() {
  const data = {
    seed: parseInt(document.getElementById('input-seed').value) || 0,
    range: parseInt(document.getElementById('select-range').value) || 20,
    V: state.V,
    W: state.W,
    Q: state.Q,
    P: state.P,
    qRoots: state.qRoots,
    qDegree: state.qDegree,
    qDiscSign: state.qDiscSign,
    category: state.category,
    hasSR: state.hasSR,
    hasB: state.hasB,
    punctures: state.punctures.map(p => ({
      lambda: p.lambda,
      face: p.face,
      pos3d: p.pos3d,
    })),
    segments: state.segments.map(s => ({
      lamEntry: s.lamEntry,
      lamExit: s.lamExit,
      color: s.color,
      infinitySpanning: s.infinitySpanning,
    })),
  };

  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `pvtet_seed${data.seed}.json`;
  a.click();
  URL.revokeObjectURL(url);
}

// ── Dark mode ──

function toggleDarkMode() {
  document.body.classList.toggle('dark');
  const btn = document.getElementById('btn-dark-mode');
  btn.textContent = document.body.classList.contains('dark') ? 'Light' : 'Dark';
  localStorage.setItem('pvtet-dark', document.body.classList.contains('dark') ? '1' : '0');
}

// Restore dark mode from localStorage
if (localStorage.getItem('pvtet-dark') === '1') {
  document.body.classList.add('dark');
  // Button text updated after DOM ready
  setTimeout(() => {
    const btn = document.getElementById('btn-dark-mode');
    if (btn) btn.textContent = 'Light';
  }, 0);
}

// ── Keyboard shortcuts ──

function onKeyDown(e) {
  // Don't intercept when typing in inputs
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.tagName === 'TEXTAREA') return;

  const seedInput = document.getElementById('input-seed');
  const seed = parseInt(seedInput.value) || 0;

  if (e.key === 'ArrowUp' || e.key === 'ArrowRight') {
    e.preventDefault();
    seedInput.value = seed + 1;
    doLoadSeed();
  } else if (e.key === 'ArrowDown' || e.key === 'ArrowLeft') {
    e.preventDefault();
    seedInput.value = Math.max(0, seed - 1);
    doLoadSeed();
  }
}

// ── URL hash state ──

function updateHash() {
  const seed = document.getElementById('input-seed').value || '0';
  const range = document.getElementById('select-range').value || '20';
  let hash = `s=${seed}&r=${range}`;
  if (state.dim === '2d') hash += '&d=2d';
  if (location.hash !== '#' + hash) {
    history.replaceState(null, '', '#' + hash);
  }
}

function loadFromHash() {
  const hash = location.hash.replace('#', '');
  if (!hash) return;
  const params = new URLSearchParams(hash);
  const seed = params.get('s');
  const range = params.get('r');
  const dim = params.get('d');

  if (dim === '2d' && state.dim !== '2d') {
    state.setDim('2d');
  }

  if (seed !== null) {
    document.getElementById('input-seed').value = seed;
  }
  if (range !== null) {
    const sel = document.getElementById('select-range');
    for (const opt of sel.options) {
      if (opt.value === range) { sel.value = range; break; }
    }
  }
  if (seed !== null) {
    doLoadSeed();
  }
}
