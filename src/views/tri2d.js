/**
 * 2D Triangle visualization using Three.js orthographic camera.
 * Feature-parity with tet3d.js: cursor↔ring linking, Cv/Cw/SR/∞/λ=0 markers,
 * puncture coloring, segment hover highlighting.
 */
import * as THREE from 'three';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { state } from '../state.js';
import {
  TRI_VERTS_2D, baryTriTo2D, lambdaToBary2D, SEGMENT_COLORS_2D,
} from '../math/curves2d.js';

let renderer, labelRenderer, scene, camera;
let container;
let triGroup, arrowGroup, curveGroup, punctureGroup, hoverGroup, markerGroup;

// Store curve Line2 meshes for raycasting
const curveMeshInfo = [];

export function initTri2D() {
  container = document.getElementById('panel-2d');
  if (!container) return;

  const canvas = document.getElementById('tri-canvas');
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(window.devicePixelRatio);

  labelRenderer = new CSS2DRenderer();
  labelRenderer.domElement.style.position = 'absolute';
  labelRenderer.domElement.style.top = '0';
  labelRenderer.domElement.style.left = '0';
  labelRenderer.domElement.style.pointerEvents = 'none';
  container.appendChild(labelRenderer.domElement);

  scene = new THREE.Scene();
  camera = new THREE.OrthographicCamera(-0.3, 1.3, 1.2, -0.3, -1, 1);
  camera.position.set(0, 0, 1);
  camera.lookAt(0, 0, 0);

  triGroup = new THREE.Group();
  arrowGroup = new THREE.Group();
  curveGroup = new THREE.Group();
  punctureGroup = new THREE.Group();
  hoverGroup = new THREE.Group();
  markerGroup = new THREE.Group();
  scene.add(triGroup, arrowGroup, curveGroup, punctureGroup, markerGroup, hoverGroup);

  buildTriangle();

  state.on('dataChanged', () => { if (state.dim === '2d') updateAll(); });
  state.on('hoverChanged', () => { if (state.dim === '2d') updateHover(); });
  state.on('modeChanged', onModeChanged);

  // Mouse interaction for cursor↔ring linking
  renderer.domElement.addEventListener('mousemove', onMouseMove);
  renderer.domElement.addEventListener('mouseleave', () => {
    if (state.dim === '2d') state.clearHover();
  });

  window.addEventListener('resize', onResize);
  onResize();
  if (state.dim === '2d') updateAll();
  animate();
}

function buildTriangle() {
  while (triGroup.children.length > 0) triGroup.remove(triGroup.children[0]);

  // Edges
  const triGeo = new THREE.BufferGeometry();
  const verts = [];
  for (const [a, b] of [[0,1],[1,2],[2,0]]) {
    verts.push(TRI_VERTS_2D[a][0], TRI_VERTS_2D[a][1], 0);
    verts.push(TRI_VERTS_2D[b][0], TRI_VERTS_2D[b][1], 0);
  }
  triGeo.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
  triGroup.add(new THREE.LineSegments(triGeo, new THREE.LineBasicMaterial({ color: 0x888888 })));

  // Semi-transparent fill
  const fillGeo = new THREE.BufferGeometry();
  fillGeo.setAttribute('position', new THREE.Float32BufferAttribute([
    TRI_VERTS_2D[0][0], TRI_VERTS_2D[0][1], -0.01,
    TRI_VERTS_2D[1][0], TRI_VERTS_2D[1][1], -0.01,
    TRI_VERTS_2D[2][0], TRI_VERTS_2D[2][1], -0.01,
  ], 3));
  triGroup.add(new THREE.Mesh(fillGeo, new THREE.MeshBasicMaterial({
    color: 0xdddddd, transparent: true, opacity: 0.15, side: THREE.DoubleSide,
  })));

  // Vertex labels
  const labels = ['v\u2080', 'v\u2081', 'v\u2082'];
  const offsets = [[-0.04, -0.04], [1.04, -0.04], [0.5, 0.93]];
  for (let i = 0; i < 3; i++) {
    const div = document.createElement('div');
    div.textContent = labels[i];
    div.style.cssText = 'font-size:14px; color:#555; font-style:italic; font-family:serif;';
    const obj = new CSS2DObject(div);
    obj.position.set(offsets[i][0], offsets[i][1], 0);
    triGroup.add(obj);
  }
}

function updateAll() {
  updateArrows();
  updateCurves();
  updatePunctures();
  updateMarkers();
  updateHover();
}

// ── Arrows ──

function updateArrows() {
  while (arrowGroup.children.length > 0) arrowGroup.remove(arrowGroup.children[0]);
  if (state.dim !== '2d') return;

  const V = state.V2d, W = state.W2d;
  const scale = 0.04;
  const headLen = 0.04, headWid = 0.03;

  for (let i = 0; i < 3; i++) {
    const origin = new THREE.Vector3(TRI_VERTS_2D[i][0], TRI_VERTS_2D[i][1], 0);

    const vLen = Math.sqrt(V[i][0] ** 2 + V[i][1] ** 2);
    if (vLen > 1e-10) {
      const len = Math.min(vLen * scale, 0.18);
      const arrow = new THREE.ArrowHelper(
        new THREE.Vector3(V[i][0] / vLen, V[i][1] / vLen, 0),
        origin, len, 0xcc0000, headLen, headWid);
      arrowGroup.add(arrow);
    }

    const wLen = Math.sqrt(W[i][0] ** 2 + W[i][1] ** 2);
    if (wLen > 1e-10) {
      const len = Math.min(wLen * scale, 0.18);
      const arrow = new THREE.ArrowHelper(
        new THREE.Vector3(W[i][0] / wLen, W[i][1] / wLen, 0),
        origin, len, 0x0066cc, headLen, headWid);
      arrowGroup.add(arrow);
    }
  }
}

// ── Curves + feature markers ──

function updateCurves() {
  while (curveGroup.children.length > 0) curveGroup.remove(curveGroup.children[0]);
  curveMeshInfo.length = 0;
  if (state.dim !== '2d') return;

  const w = container ? container.clientWidth : 800;
  const h = container ? container.clientHeight : 600;

  for (let si = 0; si < state.segments.length; si++) {
    const seg = state.segments[si];
    const color = new THREE.Color(seg.color);
    for (let pi = 0; pi < seg.ptsList.length; pi++) {
      const pts = seg.ptsList[pi];
      if (pts.length < 2) continue;
      const positions = [];
      for (const pt of pts) positions.push(pt[0], pt[1], 0.01);
      const geo = new LineGeometry();
      geo.setPositions(positions);
      const mat = new LineMaterial({ color, linewidth: 3, resolution: new THREE.Vector2(w, h) });
      const line = new Line2(geo, mat);
      curveGroup.add(line);
      curveMeshInfo.push({ mesh: line, segIdx: si, subIdx: pi });
    }
  }

  // Bubble is now drawn as a regular segment (has ptsList+lamsList)
}

// ── Feature markers: Cv, Cw, SR, λ=0, ∞ ──

function updateMarkers() {
  while (markerGroup.children.length > 0) markerGroup.remove(markerGroup.children[0]);
  if (state.dim !== '2d') return;

  const Q = state.Q, P = state.P;

  // λ=0 marker (Cv) — only if classifier confirmed inside
  if (state.hasCvPos && Q[0] !== 0) {
    const mu0 = P.map(pk => pk[0] / Q[0]);
    const muClip = mu0.map(m => Math.max(0, m));
    const s = muClip.reduce((a, b) => a + b, 0);
    if (s > 1e-10) {
      const pos = baryTriTo2D(muClip.map(m => m / s));
      addMarker(pos, '\u03BB=0', '#228B22', 0.012);
    }
  }

  // λ=∞ marker (Cw) — only if classifier confirmed inside
  if (state.hasCwPos) {
    let dQ = Q.length - 1;
    while (dQ > 0 && Math.abs(Q[dQ]) < 1e-30) dQ--;
    if (dQ > 0 && Q[dQ] !== 0) {
      const muInf = P.map(pk => (pk.length > dQ ? pk[dQ] : 0) / Q[dQ]);
      const muClip = muInf.map(m => Math.max(0, m));
      const s = muClip.reduce((a, b) => a + b, 0);
      if (s > 1e-10) {
        const pos = baryTriTo2D(muClip.map(m => m / s));
        addMarker(pos, '\u221E', '#8B4513', 0.012);
      }
    }
  }

  // SR marker — use exact BigInt resultant to find which Q-root is shared
  if (state.hasSR && state.qRoots) {
    // For 3D: state.srPos3d is set by classifier
    // For 2D: find SR via BigInt GCD of Q and P[k]
    if (state.srLambda !== null && state.srLambda !== undefined) {
      // 3D path — srLambda is populated
      const root = state.srLambda;
      const qDeriv = Q[1] + 2 * (Q[2] || 0) * root + 3 * (Q[3] || 0) * root * root;
      if (Math.abs(qDeriv) > 0) {
        const mu = P.map(pk => {
          let d = 0;
          for (let i = pk.length - 1; i >= 0; i--) d = d * root + (i >= 1 ? i * pk[i] : 0);
          return d / qDeriv;
        });
        const muClip = mu.map(m => Math.max(0, m));
        const s = muClip.reduce((a, b) => a + b, 0);
        if (s > 0) addMarker(baryTriTo2D(muClip.map(m => m / s)), 'SR', '#ff00ff', 0.015);
      }
    } else {
      // 2D path — compute SR position from BigInt polynomials via L'Hôpital
      // Use float root as approximation, position via derivative ratio
      for (const root of state.qRoots) {
        const qDeriv = Q[1] + 2 * Q[2] * root;
        if (qDeriv === 0) continue;
        const mu = P.map(pk => (pk[1] + 2 * pk[2] * root) / qDeriv);
        const muClip = mu.map(m => Math.max(0, m));
        const s = muClip.reduce((a, b) => a + b, 0);
        if (s > 0 && mu.every(m => m >= -0.2)) {
          addMarker(baryTriTo2D(muClip.map(m => m / s)), 'SR', '#ff00ff', 0.015);
          break;
        }
      }
    }
  }
}

function addMarker(pos, label, colorStr, size) {
  // Filled circle
  const geo = new THREE.CircleGeometry(size, 20);
  const mat = new THREE.MeshBasicMaterial({ color: new THREE.Color(colorStr) });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.position.set(pos[0], pos[1], 0.025);
  markerGroup.add(mesh);

  // Outline ring
  const ringGeo = new THREE.RingGeometry(size, size * 1.3, 20);
  const ringMat = new THREE.MeshBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.5 });
  const ring = new THREE.Mesh(ringGeo, ringMat);
  ring.position.set(pos[0], pos[1], 0.026);
  markerGroup.add(ring);

  // CSS2D label
  const div = document.createElement('div');
  div.textContent = label;
  div.style.color = colorStr;
  div.style.fontSize = '14px';
  div.style.fontWeight = 'bold';
  div.style.fontFamily = 'sans-serif';
  div.style.pointerEvents = 'none';
  div.style.textShadow = '0 0 3px white, 0 0 3px white';
  const label2d = new CSS2DObject(div);
  label2d.position.set(pos[0], pos[1] + 0.035, 0);
  markerGroup.add(label2d);
}

// ── Punctures (colored by segment) ──

function updatePunctures() {
  while (punctureGroup.children.length > 0) punctureGroup.remove(punctureGroup.children[0]);
  if (state.dim !== '2d') return;

  // Build puncture → color map from segments
  const puncColor = {};
  for (const seg of state.segments) {
    if (!(seg.piEntry in puncColor)) puncColor[seg.piEntry] = seg.color;
    if (!(seg.piExit in puncColor)) puncColor[seg.piExit] = seg.color;
  }

  for (let i = 0; i < state.punctures.length; i++) {
    const p = state.punctures[i];
    const pos2d = baryTriTo2D(p.bary || p.mu);
    const color = new THREE.Color(puncColor[i] || '#666666');
    const size = p.isEdge ? 0.02 : 0.015;

    const geo = new THREE.CircleGeometry(size, 16);
    const mat = new THREE.MeshBasicMaterial({ color });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(pos2d[0], pos2d[1], 0.02);
    mesh.userData.punctureIdx = i;
    punctureGroup.add(mesh);

    // Black outline
    const ringGeo = new THREE.RingGeometry(size * 0.9, size * 1.15, 16);
    const ringMat = new THREE.MeshBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.5 });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.set(pos2d[0], pos2d[1], 0.021);
    punctureGroup.add(ring);
  }
}

// ── Hover (linked to ring) ──

function updateHover() {
  while (hoverGroup.children.length > 0) hoverGroup.remove(hoverGroup.children[0]);
  if (state.dim !== '2d') return;

  // Highlight hovered puncture
  if (state.hoverPunctureIdx >= 0 && state.hoverPunctureIdx < state.punctures.length) {
    const p = state.punctures[state.hoverPunctureIdx];
    const pos2d = baryTriTo2D(p.bary || p.mu);
    const geo = new THREE.RingGeometry(0.025, 0.04, 24);
    const mat = new THREE.MeshBasicMaterial({ color: 0xffff00, transparent: true, opacity: 0.8 });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(pos2d[0], pos2d[1], 0.03);
    hoverGroup.add(mesh);
  }

  // Highlight hovered segment (thicker re-draw)
  if (state.hoverSegmentIdx >= 0 && state.hoverSegmentIdx < state.segments.length) {
    const seg = state.segments[state.hoverSegmentIdx];
    const color = new THREE.Color(seg.color);
    const w = container ? container.clientWidth : 800;
    const h = container ? container.clientHeight : 600;
    for (const pts of seg.ptsList) {
      if (pts.length < 2) continue;
      const positions = [];
      for (const pt of pts) positions.push(pt[0], pt[1], 0.015);
      const geo = new LineGeometry();
      geo.setPositions(positions);
      const mat = new LineMaterial({ color, linewidth: 6, resolution: new THREE.Vector2(w, h) });
      hoverGroup.add(new Line2(geo, mat));
    }
  }

  // Ring hover → 2D position + λ label
  if (state.hoverRingPos3d && state.dim === '2d') {
    const pos = state.hoverRingPos3d;

    // Cursor dot
    const geo = new THREE.CircleGeometry(0.015, 16);
    const mat = new THREE.MeshBasicMaterial({ color: 0xff8800, transparent: true, opacity: 0.9 });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(pos[0], pos[1], 0.03);
    hoverGroup.add(mesh);

    // Crosshair ring
    const ringGeo = new THREE.RingGeometry(0.02, 0.03, 20);
    const ringMat = new THREE.MeshBasicMaterial({ color: 0xff8800, transparent: true, opacity: 0.5 });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.set(pos[0], pos[1], 0.031);
    hoverGroup.add(ring);

    // Lambda value label
    if (state.hoverLambda !== null) {
      const lamStr = isFinite(state.hoverLambda)
        ? `\u03BB=${state.hoverLambda.toFixed(3)}`
        : '\u03BB=\u221E';
      const div = document.createElement('div');
      div.textContent = lamStr;
      div.style.color = '#ff8800';
      div.style.fontSize = '13px';
      div.style.fontWeight = 'bold';
      div.style.fontFamily = 'sans-serif';
      div.style.pointerEvents = 'none';
      div.style.textShadow = '0 0 3px white, 0 0 3px white, 0 0 5px white';
      const label2d = new CSS2DObject(div);
      label2d.position.set(pos[0], pos[1] + 0.04, 0);
      hoverGroup.add(label2d);
    }
  }
}

// ── Mouse → cursor↔ring linking ──

function onMouseMove(event) {
  if (state.dim !== '2d') return;
  const rect = renderer.domElement.getBoundingClientRect();
  const mx = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  const my = -((event.clientY - rect.top) / rect.height) * 2 + 1;

  // Unproject mouse to world coords (orthographic)
  const worldX = (mx - (-1)) / 2 * (camera.right - camera.left) + camera.left;
  const worldY = (my - (-1)) / 2 * (camera.top - camera.bottom) + camera.bottom;

  // Find closest curve point
  let bestDist = Infinity, bestLam = null, bestSegIdx = -1;
  for (const { segIdx, subIdx } of curveMeshInfo) {
    const seg = state.segments[segIdx];
    if (!seg) continue;
    const pts = seg.ptsList[subIdx];
    const lams = seg.lamsList && seg.lamsList[subIdx];
    if (!pts || !lams) continue;
    for (let i = 0; i < pts.length; i++) {
      const dx = pts[i][0] - worldX, dy = pts[i][1] - worldY;
      const d = dx * dx + dy * dy;
      if (d < bestDist) { bestDist = d; bestLam = lams[i]; bestSegIdx = segIdx; }
    }
  }

  // Threshold: ~15px in world coords
  const pxSize = (camera.right - camera.left) / (container.clientWidth || 1);
  const threshold = (15 * pxSize) ** 2;

  if (bestDist < threshold && bestLam !== null) {
    state.setHoverLambda(bestLam);
    return;
  }

  // No hit
  if (state.hoverLambda !== null && state.hoverRingPos3d === null) {
    state.clearHover();
  }
}

// ── Mode / resize / animate ──

function onModeChanged() {
  if (state.dim === '2d') {
    setTimeout(() => { onResize(); updateAll(); }, 50);
  }
}

function onResize() {
  if (!container || !renderer) return;
  const w = container.clientWidth, h = container.clientHeight;
  if (w === 0 || h === 0) return;
  renderer.setSize(w, h);
  labelRenderer.setSize(w, h);

  const aspect = w / h;
  const pad = 0.3;
  if (aspect > 1) {
    camera.left = -pad;
    camera.right = 1 + pad;
    const vSize = (1 + 2 * pad) / aspect;
    camera.bottom = -pad;
    camera.top = -pad + vSize;
  } else {
    camera.bottom = -pad;
    camera.top = Math.sqrt(3) / 2 + pad;
    const hSize = (Math.sqrt(3) / 2 + 2 * pad) * aspect;
    camera.left = 0.5 - hSize / 2;
    camera.right = 0.5 + hSize / 2;
  }
  camera.updateProjectionMatrix();
}

function animate() {
  requestAnimationFrame(animate);
  if (!renderer || !container) return;
  if (state.dim !== '2d') return;
  if (container.offsetParent === null) return;

  const isDark = document.body.classList.contains('dark');
  scene.background = new THREE.Color(isDark ? 0x1e1e3e : 0xffffff);

  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
}
