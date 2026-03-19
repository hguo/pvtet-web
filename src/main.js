/**
 * Entry point: initialize all views and wire events.
 */
import { initTet3D } from './views/tet3d.js';
import { initTri2D } from './views/tri2d.js';
import { initRing } from './views/ring.js';
import { initInfo } from './views/info.js';
import { initControls } from './controls.js';

function init() {
  initControls();
  initTet3D();
  initTri2D();
  initRing();
  initInfo();
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', init);
} else {
  init();
}
