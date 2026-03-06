/**
 * Three.js 3D tetrahedron view: wireframe, arrows, PV curves (tubes),
 * puncture markers, infinity markers.
 * Uses TrackballControls for gimbal-free quaternion rotation.
 */
import * as THREE from 'three';
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { state } from '../state.js';
import { TET_VERTS, TET_EDGES, SEGMENT_COLORS, baryTetTo3D } from '../math/curves.js';

let renderer, labelRenderer, scene, camera, controls;
let container;
let initialCamPos, initialCamUp, initialTarget;

// Groups for dynamic content
let arrowGroup, curveGroup, punctureGroup, hoverGroup;
// Static tet wireframe + faces
let tetGroup;

const ARROW_SCALE = 0.10;
const TUBE_RADIUS = 0.006;
const TUBE_SEGMENTS = 6;

// Raycasting for curve hover
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();
// Map tube mesh → { segIdx, subIdx } for λ lookup
const tubeMeshMap = new Map();

export function initTet3D() {
  container = document.getElementById('panel-3d');
  const canvas = document.getElementById('tet-canvas');

  // Renderer
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setClearColor(0xffffff, 1);

  // CSS2D label renderer
  labelRenderer = new CSS2DRenderer();
  labelRenderer.domElement.style.position = 'absolute';
  labelRenderer.domElement.style.top = '0';
  labelRenderer.domElement.style.left = '0';
  labelRenderer.domElement.style.pointerEvents = 'none';
  container.appendChild(labelRenderer.domElement);

  // Scene
  scene = new THREE.Scene();

  // Camera
  camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);
  const cx = 0.5, cy = 0.3, cz = 0.2; // tet centroid approx
  camera.position.set(cx - 1.0, cy + 0.5, cz + 1.5);
  camera.lookAt(cx, cy, cz);

  // Save initial camera state
  initialCamPos = camera.position.clone();
  initialCamUp = camera.up.clone();
  initialTarget = new THREE.Vector3(cx, cy, cz);

  // TrackballControls: quaternion-based, no gimbal lock
  controls = new TrackballControls(camera, renderer.domElement);
  controls.target.copy(initialTarget);
  controls.rotateSpeed = 3.0;
  controls.zoomSpeed = 1.5;
  controls.panSpeed = 0.8;
  controls.dynamicDampingFactor = 0.15;
  controls.update();

  // Light — balanced for visible specular on Phong materials
  scene.add(new THREE.AmbientLight(0xffffff, 0.4));
  const dlight = new THREE.DirectionalLight(0xffffff, 0.7);
  dlight.position.set(2, 3, 2);
  scene.add(dlight);
  const dlight2 = new THREE.DirectionalLight(0xffffff, 0.3);
  dlight2.position.set(-2, -1, 3);
  scene.add(dlight2);

  // Static tet wireframe
  tetGroup = new THREE.Group();
  buildTetWireframe(tetGroup);
  scene.add(tetGroup);

  // Dynamic groups
  arrowGroup = new THREE.Group();
  scene.add(arrowGroup);
  curveGroup = new THREE.Group();
  scene.add(curveGroup);
  punctureGroup = new THREE.Group();
  scene.add(punctureGroup);
  hoverGroup = new THREE.Group();
  scene.add(hoverGroup);

  // Event listeners
  state.on('dataChanged', updateAll);
  state.on('hoverChanged', updateHover);
  window.addEventListener('resize', onResize);

  // Reset view button — restore initial camera state
  document.getElementById('btn-reset-view').addEventListener('click', () => {
    camera.position.copy(initialCamPos);
    camera.up.copy(initialCamUp);
    controls.target.copy(initialTarget);
    controls.update();
  });

  // Curve hover via raycasting
  renderer.domElement.addEventListener('mousemove', onMouseMove);
  renderer.domElement.addEventListener('mouseleave', () => state.clearHover());

  onResize();
  updateAll();
  animate();
}

function onMouseMove(event) {
  const rect = renderer.domElement.getBoundingClientRect();
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera(mouse, camera);
  const intersects = raycaster.intersectObjects([...curveGroup.children], false);

  if (intersects.length > 0) {
    const hit = intersects[0];
    const info = tubeMeshMap.get(hit.object);
    if (info) {
      const seg = state.segments[info.segIdx];
      const pts = seg.ptsList[info.subIdx];
      const lams = seg.lamsList && seg.lamsList[info.subIdx];
      if (pts && lams) {
        // Find closest stored point to hit position
        const hp = hit.point;
        let bestDist = Infinity, bestLam = null;
        for (let i = 0; i < pts.length; i++) {
          const dx = pts[i][0] - hp.x, dy = pts[i][1] - hp.y, dz = pts[i][2] - hp.z;
          const d = dx * dx + dy * dy + dz * dz;
          if (d < bestDist) { bestDist = d; bestLam = lams[i]; }
        }
        if (bestLam !== null) {
          state.setHoverLambda(bestLam);
          return;
        }
      }
    }
  }
  // No hit — only clear if we were hovering from 3D (not from ring)
  if (state.hoverLambda !== null && state.hoverRingPos3d !== null) {
    // Don't clear if ring is driving the hover
  } else if (state.hoverLambda !== null) {
    state.clearHover();
  }
}

function buildTetWireframe(group) {
  // Edges
  const edgeGeo = new THREE.BufferGeometry();
  const edgeVerts = [];
  for (const [i, j] of TET_EDGES) {
    edgeVerts.push(...TET_VERTS[i], ...TET_VERTS[j]);
  }
  edgeGeo.setAttribute('position', new THREE.Float32BufferAttribute(edgeVerts, 3));
  const edgeMat = new THREE.LineBasicMaterial({ color: 0x555555, linewidth: 1 });
  group.add(new THREE.LineSegments(edgeGeo, edgeMat));

  // Transparent faces
  const faceIndices = [
    [1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2],
  ];
  for (const tri of faceIndices) {
    const geo = new THREE.BufferGeometry();
    const verts = new Float32Array(9);
    for (let k = 0; k < 3; k++) {
      verts[k * 3] = TET_VERTS[tri[k]][0];
      verts[k * 3 + 1] = TET_VERTS[tri[k]][1];
      verts[k * 3 + 2] = TET_VERTS[tri[k]][2];
    }
    geo.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
    const mat = new THREE.MeshBasicMaterial({
      color: 0xe0e0e0, transparent: true, opacity: 0.08,
      side: THREE.DoubleSide, depthWrite: false,
    });
    group.add(new THREE.Mesh(geo, mat));
  }

  // Vertex labels
  for (let i = 0; i < 4; i++) {
    const div = document.createElement('div');
    div.textContent = `v${String.fromCharCode(0x2080 + i)}`;
    div.style.color = '#333';
    div.style.fontSize = '18px';
    div.style.fontWeight = '500';
    div.style.fontFamily = 'serif';
    div.style.pointerEvents = 'none';
    const label = new CSS2DObject(div);
    const cx = TET_VERTS.reduce((s, v) => s + v[0], 0) / 4;
    const cy = TET_VERTS.reduce((s, v) => s + v[1], 0) / 4;
    const cz = TET_VERTS.reduce((s, v) => s + v[2], 0) / 4;
    const dx = TET_VERTS[i][0] - cx;
    const dy = TET_VERTS[i][1] - cy;
    const dz = TET_VERTS[i][2] - cz;
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1;
    label.position.set(
      TET_VERTS[i][0] + dx / len * 0.08,
      TET_VERTS[i][1] + dy / len * 0.08,
      TET_VERTS[i][2] + dz / len * 0.08,
    );
    group.add(label);
  }
}

function updateAll() {
  updateArrows();
  updateCurves();
  updatePunctures();
  updateHover();
}

function updateArrows() {
  arrowGroup.clear();
  const V = state.V;
  const W = state.W;
  let vMax = 0, wMax = 0;
  for (let i = 0; i < 4; i++) {
    for (let d = 0; d < 3; d++) {
      vMax = Math.max(vMax, Math.abs(V[i][d]));
      wMax = Math.max(wMax, Math.abs(W[i][d]));
    }
  }
  vMax = Math.max(vMax, 1e-10);
  wMax = Math.max(wMax, 1e-10);

  for (let i = 0; i < 4; i++) {
    const origin = new THREE.Vector3(...TET_VERTS[i]);

    // V arrow (red)
    const vDir = new THREE.Vector3(V[i][0] / vMax, V[i][1] / vMax, V[i][2] / vMax);
    const vLen = vDir.length() * ARROW_SCALE;
    if (vLen > 1e-6) {
      vDir.normalize();
      const arrow = new THREE.ArrowHelper(vDir, origin, vLen, 0xcc0000, vLen * 0.3, vLen * 0.15);
      arrow.line.material.transparent = true;
      arrow.line.material.opacity = 0.8;
      arrow.cone.material.transparent = true;
      arrow.cone.material.opacity = 0.8;
      arrowGroup.add(arrow);
    }

    // W arrow (blue)
    const wDir = new THREE.Vector3(W[i][0] / wMax, W[i][1] / wMax, W[i][2] / wMax);
    const wLen = wDir.length() * ARROW_SCALE;
    if (wLen > 1e-6) {
      wDir.normalize();
      const arrow = new THREE.ArrowHelper(wDir, origin, wLen, 0x0044cc, wLen * 0.3, wLen * 0.15);
      arrow.line.material.transparent = true;
      arrow.line.material.opacity = 0.8;
      arrow.cone.material.transparent = true;
      arrow.cone.material.opacity = 0.8;
      arrowGroup.add(arrow);
    }
  }
}

/** Build a TubeGeometry from an array of [x,y,z] points. */
function buildTube(pts, color, radius = TUBE_RADIUS) {
  const vectors = pts.map(p => new THREE.Vector3(p[0], p[1], p[2]));
  const curve = new THREE.CatmullRomCurve3(vectors, false, 'centripetal', 0.3);
  const geo = new THREE.TubeGeometry(curve, Math.max(pts.length, 20), radius, TUBE_SEGMENTS, false);
  const mat = new THREE.MeshPhongMaterial({
    color, shininess: 60, transparent: true, opacity: 0.9,
  });
  return new THREE.Mesh(geo, mat);
}

function updateCurves() {
  curveGroup.clear();
  tubeMeshMap.clear();

  // Draw PV curve segments as illuminated tubes
  for (let si = 0; si < state.segments.length; si++) {
    const seg = state.segments[si];
    const color = new THREE.Color(seg.color);
    for (let pi = 0; pi < seg.ptsList.length; pi++) {
      const pts = seg.ptsList[pi];
      if (pts.length < 2) continue;
      const tube = buildTube(pts, color);
      tubeMeshMap.set(tube, { segIdx: si, subIdx: pi });
      curveGroup.add(tube);
    }

    // Infinity marker: sphere + label at λ=∞ position
    if (seg.infPos3d) {
      const pos = seg.infPos3d;
      const segColor = new THREE.Color(seg.color);

      // Sphere marker
      const geo = new THREE.SphereGeometry(0.012, 12, 8);
      const mat = new THREE.MeshPhongMaterial({
        color: segColor, shininess: 80, transparent: true, opacity: 0.9,
      });
      const sphere = new THREE.Mesh(geo, mat);
      sphere.position.set(pos[0], pos[1], pos[2]);
      curveGroup.add(sphere);

      // Black outline ring
      const ringGeo = new THREE.RingGeometry(0.013, 0.017, 16);
      const ringMat = new THREE.MeshBasicMaterial({
        color: 0x000000, side: THREE.DoubleSide, transparent: true, opacity: 0.6,
      });
      const ring = new THREE.Mesh(ringGeo, ringMat);
      ring.position.copy(sphere.position);
      ring.lookAt(camera.position);
      curveGroup.add(ring);

      // CSS2D label "\u221E"
      const div = document.createElement('div');
      div.textContent = '\u221E';
      div.style.color = seg.color;
      div.style.fontSize = '18px';
      div.style.fontWeight = 'bold';
      div.style.fontFamily = 'serif';
      div.style.pointerEvents = 'none';
      div.style.textShadow = '0 0 3px white, 0 0 3px white';
      const label2d = new CSS2DObject(div);
      label2d.position.set(pos[0], pos[1] + 0.03, pos[2]);
      curveGroup.add(label2d);
    }

    // Zero marker: sphere + label at λ=0 position
    if (seg.zeroPos3d) {
      const pos = seg.zeroPos3d;
      const segColor = new THREE.Color(seg.color);

      const geo = new THREE.SphereGeometry(0.012, 12, 8);
      const mat = new THREE.MeshPhongMaterial({
        color: segColor, shininess: 80, transparent: true, opacity: 0.9,
      });
      const sphere = new THREE.Mesh(geo, mat);
      sphere.position.set(pos[0], pos[1], pos[2]);
      curveGroup.add(sphere);

      const ringGeo = new THREE.RingGeometry(0.013, 0.017, 16);
      const ringMat = new THREE.MeshBasicMaterial({
        color: 0x000000, side: THREE.DoubleSide, transparent: true, opacity: 0.6,
      });
      const ring = new THREE.Mesh(ringGeo, ringMat);
      ring.position.copy(sphere.position);
      ring.lookAt(camera.position);
      curveGroup.add(ring);

      const div = document.createElement('div');
      div.textContent = '\u03BB=0';
      div.style.color = seg.color;
      div.style.fontSize = '18px';
      div.style.fontWeight = 'bold';
      div.style.fontFamily = 'serif';
      div.style.pointerEvents = 'none';
      div.style.textShadow = '0 0 3px white, 0 0 3px white';
      const label2d = new CSS2DObject(div);
      label2d.position.set(pos[0], pos[1] + 0.03, pos[2]);
      curveGroup.add(label2d);
    }
  }

  // SR marker: sphere + label at shared-root position
  if (state.srPos3d) {
    const pos = state.srPos3d;
    const geo = new THREE.SphereGeometry(0.015, 12, 8);
    const mat = new THREE.MeshPhongMaterial({
      color: 0xff00ff, shininess: 80, transparent: true, opacity: 0.9,
    });
    const sphere = new THREE.Mesh(geo, mat);
    sphere.position.set(pos[0], pos[1], pos[2]);
    curveGroup.add(sphere);

    const ringGeo = new THREE.RingGeometry(0.016, 0.022, 16);
    const ringMat = new THREE.MeshBasicMaterial({
      color: 0x000000, side: THREE.DoubleSide, transparent: true, opacity: 0.6,
    });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.copy(sphere.position);
    ring.lookAt(camera.position);
    curveGroup.add(ring);

    const div = document.createElement('div');
    div.textContent = 'SR';
    div.style.color = '#ff00ff';
    div.style.fontSize = '18px';
    div.style.fontWeight = 'bold';
    div.style.fontFamily = 'sans-serif';
    div.style.pointerEvents = 'none';
    div.style.textShadow = '0 0 3px white, 0 0 3px white';
    const label2d = new CSS2DObject(div);
    label2d.position.set(pos[0], pos[1] + 0.04, pos[2]);
    curveGroup.add(label2d);
  }

  // Bubble
  if (state.bubble && state.bubble.length > 2) {
    const color = new THREE.Color(SEGMENT_COLORS[0]);
    curveGroup.add(buildTube(state.bubble, color));
  }
}

function updatePunctures() {
  punctureGroup.clear();

  // Build puncture → color map
  const puncColor = {};
  for (const seg of state.segments) {
    if (!(seg.piEntry in puncColor)) puncColor[seg.piEntry] = seg.color;
    if (!(seg.piExit in puncColor)) puncColor[seg.piExit] = seg.color;
  }

  for (let i = 0; i < state.punctures.length; i++) {
    const p = state.punctures[i];
    const color = new THREE.Color(puncColor[i] || '#666666');
    const size = p.isVertex ? 0.025 : p.isEdge ? 0.02 : 0.015;

    const geo = new THREE.SphereGeometry(size, 12, 8);
    const mat = new THREE.MeshPhongMaterial({ color, shininess: 80 });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(p.pos3d[0], p.pos3d[1], p.pos3d[2]);
    mesh.userData.punctureIdx = i;
    punctureGroup.add(mesh);

    // Black edge ring
    const ringGeo = new THREE.RingGeometry(size * 0.9, size * 1.1, 16);
    const ringMat = new THREE.MeshBasicMaterial({
      color: 0x000000, side: THREE.DoubleSide, transparent: true, opacity: 0.5,
    });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.copy(mesh.position);
    ring.lookAt(camera.position);
    punctureGroup.add(ring);
  }
}

function updateHover() {
  hoverGroup.clear();

  // Highlight hovered puncture
  if (state.hoverPunctureIdx >= 0 && state.hoverPunctureIdx < state.punctures.length) {
    const p = state.punctures[state.hoverPunctureIdx];
    const geo = new THREE.SphereGeometry(0.03, 16, 12);
    const mat = new THREE.MeshBasicMaterial({
      color: 0xffff00, transparent: true, opacity: 0.7,
    });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(p.pos3d[0], p.pos3d[1], p.pos3d[2]);
    hoverGroup.add(mesh);
  }

  // Highlight hovered segment (thicker tube)
  if (state.hoverSegmentIdx >= 0 && state.hoverSegmentIdx < state.segments.length) {
    const seg = state.segments[state.hoverSegmentIdx];
    const color = new THREE.Color(seg.color);
    for (const pts of seg.ptsList) {
      if (pts.length < 2) continue;
      hoverGroup.add(buildTube(pts, color, TUBE_RADIUS * 2.0));
    }
  }

  // Highlight ring hover position in 3D
  if (state.hoverRingPos3d) {
    const pos = state.hoverRingPos3d;
    const geo = new THREE.SphereGeometry(0.02, 12, 8);
    const mat = new THREE.MeshBasicMaterial({
      color: 0xff8800, transparent: true, opacity: 0.8,
    });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(pos[0], pos[1], pos[2]);
    hoverGroup.add(mesh);

    // Small cross-hair ring
    const ringGeo = new THREE.RingGeometry(0.025, 0.035, 16);
    const ringMat = new THREE.MeshBasicMaterial({
      color: 0xff8800, side: THREE.DoubleSide, transparent: true, opacity: 0.5,
    });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.set(pos[0], pos[1], pos[2]);
    ring.lookAt(camera.position);
    hoverGroup.add(ring);

    // Lambda value label
    if (state.hoverLambda !== null) {
      const lamStr = isFinite(state.hoverLambda)
        ? `\u03BB=${state.hoverLambda.toFixed(3)}`
        : '\u03BB=\u221E';
      const div = document.createElement('div');
      div.textContent = lamStr;
      div.style.color = '#ff8800';
      div.style.fontSize = '16px';
      div.style.fontWeight = 'bold';
      div.style.fontFamily = 'sans-serif';
      div.style.pointerEvents = 'none';
      div.style.textShadow = '0 0 3px white, 0 0 3px white, 0 0 5px white';
      const label2d = new CSS2DObject(div);
      label2d.position.set(pos[0], pos[1] + 0.05, pos[2]);
      hoverGroup.add(label2d);
    }
  }
}

function onResize() {
  const w = container.clientWidth;
  const h = container.clientHeight;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
  labelRenderer.setSize(w, h);
  controls.handleResize();
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
}
