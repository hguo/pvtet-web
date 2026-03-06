/**
 * V/W matrix input grids, Randomize/Reset buttons, Seed/Range controls.
 */
import { state } from './state.js';
import { randomVW } from './math/random.js';

let vInputs = []; // vInputs[i][j] = <input> for V[i][j]
let wInputs = [];

export function initControls() {
  buildMatrixGrid('v-grid', 'V', true);
  buildMatrixGrid('w-grid', 'W', false);

  document.getElementById('btn-randomize').addEventListener('click', doRandomize);
  document.getElementById('btn-reset').addEventListener('click', doReset);
  document.getElementById('input-seed').addEventListener('change', doLoadSeed);

  // Sync inputs from state
  state.on('dataChanged', syncInputsFromState);
  syncInputsFromState();
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
  const seed = parseInt(document.getElementById('input-seed').value) || 0;
  const range = parseInt(document.getElementById('select-range').value) || 20;
  const { V, W } = randomVW(seed, range);
  state.setVW(V, W);
}

function doRandomize() {
  doLoadSeed();
  // Advance seed for next click
  const seedInput = document.getElementById('input-seed');
  seedInput.value = (parseInt(seedInput.value) || 0) + 1;
}

function doReset() {
  document.getElementById('input-seed').value = '0';
  state.setVW(state.getDefaultV(), state.getDefaultW());
}
