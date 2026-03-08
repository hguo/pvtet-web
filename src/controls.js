/**
 * V/W matrix input grids, Randomize/Reset buttons, Seed/Range controls,
 * Gallery, Animation, Search, Export, Keyboard shortcuts, URL hash state.
 */
import { state } from './state.js';
import { randomVW } from './math/random.js';
import { classifyTetCase } from './math/classify.js';

let vInputs = []; // vInputs[i][j] = <input> for V[i][j]
let wInputs = [];

// Animation state
let animTimer = null;
let animSeed = 0;

export function initControls() {
  buildMatrixGrid('v-grid', 'V', true);
  buildMatrixGrid('w-grid', 'W', false);

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

  // Keyboard shortcuts
  document.addEventListener('keydown', onKeyDown);

  // Sync inputs from state
  state.on('dataChanged', syncInputsFromState);
  syncInputsFromState();

  // Load from URL hash if present
  loadFromHash();
  window.addEventListener('hashchange', loadFromHash);

  // Update hash when data changes
  state.on('dataChanged', updateHash);
}

function buildMatrixGrid(containerId, fieldName, isV) {
  const container = document.getElementById(containerId);
  container.innerHTML = '';
  const inputs = [];

  for (let i = 0; i < 4; i++) {
    // Row label
    const label = document.createElement('span');
    label.className = 'row-label';
    label.textContent = `${fieldName.toLowerCase()}${String.fromCharCode(0x2080 + i)}`;
    container.appendChild(label);

    const row = [];
    for (let j = 0; j < 3; j++) {
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
  const V = state.V.map(r => [...r]);
  const W = state.W.map(r => [...r]);

  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 3; j++) {
      if (isV) {
        V[i][j] = parseFloat(vInputs[i][j].value) || 0;
      } else {
        W[i][j] = parseFloat(wInputs[i][j].value) || 0;
      }
    }
  }

  state.setVW(V, W);
}

function syncInputsFromState() {
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 3; j++) {
      vInputs[i][j].value = state.V[i][j];
      wInputs[i][j].value = state.W[i][j];
    }
  }
}

function doLoadSeed() {
  const seedStr = document.getElementById('input-seed').value;
  if (seedStr === '') return;
  const seed = parseInt(seedStr);
  if (isNaN(seed) || seed < 0) return;
  const range = parseInt(document.getElementById('select-range').value) || 20;
  const { V, W } = randomVW(seed, range);
  state.setVW(V, W);
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
    const { V, W } = randomVW(seed, range);
    try {
      const cc = classifyTetCase(V, W);
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
  const hash = `s=${seed}&r=${range}`;
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
