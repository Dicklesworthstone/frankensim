/**
 * CHAPTER 1 — "THE HARDWARE": Interactive 3D Commercial Equipment Turntable & Articulated Cutaway (#equipment).
 *
 * Loads real Blender `.glb` assets from the verified 67-asset catalog (`catalog_index.json`),
 * supports category tabs, instant search across all 67 assets, LIVE Articulated Joint
 * animation (hinged oven/combi/reach-in doors, twin fryer basket lift, tilting braising pan,
 * rotating slicer blade), live 3D SI caliper dimension cage, Thermal Surface Preview,
 * and Wireframe X-Ray topology mode.
 */

import * as THREE from "three";
import {
  FLAGSHIP_EQUIPMENT,
  categorizeCatalogEntry,
  prettyEquipmentName,
  type CatalogEntry,
  type CatalogIndex,
  type EquipmentTabId,
} from "../sim/exemplars";
import {
  ASSET_BASE,
  bindCanvasResize,
  createScrollSafeOrbit,
  createStudioRenderer,
  createStudioScene,
  disposeObjectTree,
  loadEquipmentModel,
  observeInView,
  thermalColor,
} from "./kit";

interface AnimatedJointHandle {
  readonly obj: THREE.Object3D;
  readonly baseRot: THREE.Euler;
  readonly basePos: THREE.Vector3;
  readonly kind: "door_z" | "door_x" | "lid_x" | "lift_y" | "spin_x" | "spin_y" | "slide_z";
  readonly sign: number;
}

export function mountEquipmentTurntable(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-eq-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-eq-canvas]")!;
  const maybeRenderer = createStudioRenderer(canvas, true);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;

  const scene = createStudioScene(renderer, 0.72);
  const camera = new THREE.PerspectiveCamera(36, 16 / 10, 0.04, 60);
  camera.position.set(1.9, 1.45, 2.2);

  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.65, 0], {
    minDistance: 0.25,
    maxDistance: 12,
    maxPolarAngle: Math.PI * 0.495,
  });
  bindCanvasResize(stageEl, renderer, camera);

  // Sculpted forged-iron turntable plinth + engraved copper measurement rings + tick marks
  const plinthGroup = new THREE.Group();
  const baseDisc = new THREE.Mesh(
    new THREE.CylinderGeometry(1.12, 1.22, 0.055, 64),
    new THREE.MeshStandardMaterial({ color: 0x161210, roughness: 0.44, metalness: 0.68 }),
  );
  baseDisc.position.y = -0.0275;
  baseDisc.receiveShadow = true;
  plinthGroup.add(baseDisc);

  const bevelLip = new THREE.Mesh(
    new THREE.TorusGeometry(1.12, 0.008, 12, 64),
    new THREE.MeshStandardMaterial({ color: 0xb45309, roughness: 0.34, metalness: 0.85 }),
  );
  bevelLip.rotation.x = Math.PI / 2;
  bevelLip.position.y = 0.001;
  plinthGroup.add(bevelLip);

  for (const [radius, col, op] of [
    [0.48, 0xd97706, 0.42],
    [0.82, 0x94a3b8, 0.24],
    [1.06, 0xf59e0b, 0.48],
  ] as const) {
    const ring = new THREE.Mesh(
      new THREE.RingGeometry(radius - 0.005, radius, 64),
      new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: op, side: THREE.DoubleSide }),
    );
    ring.rotation.x = -Math.PI / 2;
    ring.position.y = 0.002;
    plinthGroup.add(ring);
  }
  scene.add(plinthGroup);

  const turntable = new THREE.Group();
  scene.add(turntable);
  const modelAnchor = new THREE.Group();
  const cageAnchor = new THREE.Group();
  turntable.add(modelAnchor, cageAnchor);

  // DOM bindings
  const tabBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-eq-tab]"));
  const pillGrid = root.querySelector<HTMLElement>("[data-eq-pills]")!;
  const searchInput = root.querySelector<HTMLInputElement>("[data-eq-search]");
  const catalogSelect = root.querySelector<HTMLSelectElement>("[data-eq-select]");
  const toggleThermal = root.querySelector<HTMLButtonElement>("[data-eq-thermal]");
  const toggleWire = root.querySelector<HTMLButtonElement>("[data-eq-wire]");
  const toggleCage = root.querySelector<HTMLButtonElement>("[data-eq-cage]");
  const toggleSpin = root.querySelector<HTMLButtonElement>("[data-eq-spin]");
  const toggleJoints = root.querySelector<HTMLButtonElement>("[data-eq-joints]");

  const metaName = root.querySelector<HTMLElement>("[data-eq-name]");
  const metaId = root.querySelector<HTMLElement>("[data-eq-id]");
  const metaDims = root.querySelector<HTMLElement>("[data-eq-dims]");
  const metaTris = root.querySelector<HTMLElement>("[data-eq-tris]");
  const metaCard = root.querySelector<HTMLElement>("[data-eq-card]");
  const metaSockets = root.querySelector<HTMLElement>("[data-eq-sockets]");
  const metaNote = root.querySelector<HTMLElement>("[data-eq-note]");

  let catalog: CatalogEntry[] = [];
  let activeTab: EquipmentTabId = "hot_line";
  let activeId = "NSF-RANGE-GAS-6B";
  let thermalPreview = false;
  let wireframePreview = false;
  let showCage = true;
  let autoSpin = true;
  let animateJoints = true;
  let inView = true;
  let userDragging = false;
  let loadSeq = 0;
  let activeModel: THREE.Group | null = null;
  let activeJoints: AnimatedJointHandle[] = [];

  controls.addEventListener("start", () => {
    userDragging = true;
  });
  controls.addEventListener("end", () => {
    userDragging = false;
  });

  function discoverArticulatedJoints(model: THREE.Group): AnimatedJointHandle[] {
    const found: AnimatedJointHandle[] = [];
    model.traverse((obj) => {
      const nm = (obj.name || "").toLowerCase();
      if (!nm || !obj.parent || obj === model) return;
      if (nm.includes("oven_door") || nm.includes("bin_door")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "door_x", sign: 1 });
      } else if (nm.includes("door_right")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "door_z", sign: -1 });
      } else if (nm.includes("door")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "door_z", sign: 1 });
      } else if (nm.includes("lid")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "lid_x", sign: 1 });
      } else if (nm.includes("basket_left") || nm.includes("basket_right") || nm.includes("lift_hood")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "lift_y", sign: nm.includes("right") ? -1 : 1 });
      } else if (nm.includes("drawer") || nm.includes("carriage")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "slide_z", sign: 1 });
      } else if (nm.includes("rotary_blade")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "spin_x", sign: 1 });
      } else if (nm.includes("s_blade") || nm.includes("caster_")) {
        found.push({ obj, baseRot: obj.rotation.clone(), basePos: obj.position.clone(), kind: "spin_y", sign: 1 });
      }
    });
    return found;
  }

  function applySurfaceMode(): void {
    if (!activeModel) return;
    const flagship = FLAGSHIP_EQUIPMENT.find((f) => f.id === activeId);
    const nomTempC = flagship?.nominalTempC ?? 140;
    const thermCol = thermalColor(nomTempC, 0, 340);

    activeModel.traverse((obj) => {
      if (!(obj instanceof THREE.Mesh)) return;
      const mats = Array.isArray(obj.material) ? obj.material : [obj.material];
      for (const m of mats) {
        if (!(m instanceof THREE.MeshStandardMaterial)) continue;
        m.wireframe = wireframePreview;
        m.transparent = wireframePreview || Boolean(m.userData.wasTransparent);
        m.opacity = wireframePreview ? 0.65 : (m.userData.baseOpacity as number | undefined) ?? 1.0;

        if (thermalPreview) {
          m.map = null;
          m.color.copy(thermCol);
          m.emissive.copy(thermCol);
          m.emissiveIntensity = nomTempC > 100 ? 0.85 : 0.35;
        } else {
          if (m.userData.baseColor instanceof THREE.Color) m.color.copy(m.userData.baseColor);
          if ("baseMap" in m.userData) m.map = m.userData.baseMap as THREE.Texture | null;
          if (m.userData.baseEmissive instanceof THREE.Color) m.emissive.copy(m.userData.baseEmissive);
          m.emissiveIntensity = (m.userData.baseEmissiveIntensity as number | undefined) ?? 0;
        }
        m.needsUpdate = true;
      }
    });
  }

  function rebuildCage(box: THREE.Box3, entry: CatalogEntry | null): void {
    while (cageAnchor.children.length > 0) {
      const c = cageAnchor.children[0];
      cageAnchor.remove(c);
      disposeObjectTree(c);
    }
    const size = box.getSize(new THREE.Vector3());
    const center = box.getCenter(new THREE.Vector3());

    // Crisp SI envelope wireframe
    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(new THREE.BoxGeometry(size.x, size.y, size.z)),
      new THREE.LineBasicMaterial({ color: 0xf59e0b, transparent: true, opacity: 0.52 }),
    );
    edges.position.copy(center);
    cageAnchor.add(edges);

    // Subtle socket beacons around perimeter indicating registered interaction sockets
    const socketCount = entry?.socket_ids.length ?? 4;
    const pinR = Math.max(0.011, Math.min(0.026, size.length() * 0.015));
    const pinGeo = new THREE.SphereGeometry(pinR, 14, 10);
    const pinMat = new THREE.MeshBasicMaterial({ color: 0x22d3ee });
    const ringGeo = new THREE.RingGeometry(pinR * 1.4, pinR * 1.9, 18);
    const ringMat = new THREE.MeshBasicMaterial({ color: 0x22d3ee, transparent: true, opacity: 0.45, side: THREE.DoubleSide });
    for (let i = 0; i < Math.min(socketCount, 8); i += 1) {
      const angle = (i / Math.max(1, Math.min(socketCount, 8))) * Math.PI * 2;
      const px = center.x + Math.cos(angle) * (size.x * 0.53);
      const py = Math.max(0.04, box.min.y + size.y * (0.28 + 0.46 * ((i % 3) / 2)));
      const pz = center.z + Math.sin(angle) * (size.z * 0.53);
      const pin = new THREE.Mesh(pinGeo, pinMat);
      pin.position.set(px, py, pz);
      const halo = new THREE.Mesh(ringGeo, ringMat);
      halo.position.set(px, py, pz);
      halo.lookAt(px * 2, py, pz * 2);
      cageAnchor.add(pin, halo);
    }
    cageAnchor.visible = showCage;
  }

  async function selectAsset(id: string): Promise<void> {
    const seq = ++loadSeq;
    activeId = id;
    const entry = catalog.find((c) => c.archetype_id === id) ?? null;
    const flagship = FLAGSHIP_EQUIPMENT.find((f) => f.id === id);

    for (const btn of Array.from(pillGrid.querySelectorAll<HTMLButtonElement>("button"))) {
      btn.setAttribute("aria-pressed", String(btn.dataset.assetId === id));
    }
    if (catalogSelect && catalogSelect.value !== id) {
      catalogSelect.value = id;
    }

    const glbFile = entry?.glb_path ?? `${id}.glb`;
    const model = await loadEquipmentModel(glbFile);
    if (seq !== loadSeq || !model) return;

    while (modelAnchor.children.length > 0) {
      const c = modelAnchor.children[0];
      modelAnchor.remove(c);
      disposeObjectTree(c);
    }
    // Stash initial transparency flags so wireframe toggle restores glass cleanly
    model.traverse((obj) => {
      if (obj instanceof THREE.Mesh) {
        const mats = Array.isArray(obj.material) ? obj.material : [obj.material];
        for (const m of mats) {
          if (m instanceof THREE.MeshStandardMaterial) {
            m.userData.wasTransparent = m.transparent;
            m.userData.baseOpacity = m.opacity;
          }
        }
      }
    });

    activeModel = model;
    activeJoints = discoverArticulatedJoints(model);
    modelAnchor.add(model);
    turntable.rotation.y = 0.52;

    model.updateMatrixWorld(true);
    const box = new THREE.Box3().setFromObject(model);
    const size = box.getSize(new THREE.Vector3());
    const center = box.getCenter(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z, 0.16);

    const plinthScale = Math.max(0.28, Math.min(2.4, maxDim * 0.85));
    plinthGroup.scale.set(plinthScale, 1, plinthScale);

    const dist = maxDim * 2.05;
    camera.position.set(dist * 0.72, center.y + maxDim * 0.52, dist * 0.88);
    controls.target.set(center.x, center.y, center.z);
    controls.minDistance = maxDim * 0.65;
    controls.maxDistance = maxDim * 6.5;
    controls.update();

    rebuildCage(box, entry);
    applySurfaceMode();

    const ext = entry?.measured_extents_m ?? [size.x, size.z, size.y];
    if (metaName) metaName.textContent = prettyEquipmentName(id);
    if (metaId) metaId.textContent = id;
    if (metaDims) {
      metaDims.textContent = `${Math.round(ext[0] * 1000)} × ${Math.round(ext[1] * 1000)} × ${Math.round(ext[2] * 1000)} mm`;
    }
    if (metaTris) {
      metaTris.textContent = entry ? `${entry.triangle_count.toLocaleString()} tris` : "verified .glb";
    }
    if (metaCard) {
      metaCard.textContent = entry?.matdb_card ?? "surfaces/ss304-brushed";
    }
    if (metaSockets) {
      const sc = entry?.socket_ids.length ?? 0;
      const jc = entry?.joint_count ?? activeJoints.length;
      metaSockets.textContent = `${sc} sockets · ${jc} joints`;
    }
    if (metaNote) {
      metaNote.textContent =
        flagship?.thermalNote ??
        `${entry?.asset_class ?? "equipment"} · ${((entry?.glb_bytes ?? 0) / 1024).toFixed(1)} KB binary verified`;
    }
    // Paint immediately so even quick scroll captures show the loaded model on frame 1!
    renderer.render(scene, camera);
  }

  function renderTabPills(query = ""): void {
    const q = query.trim().toLowerCase();
    pillGrid.innerHTML = "";

    let items: { id: string; label: string }[];
    if (q.length > 0 && catalog.length > 0) {
      items = catalog
        .filter(
          (c) =>
            c.archetype_id.toLowerCase().includes(q) ||
            prettyEquipmentName(c.archetype_id).toLowerCase().includes(q) ||
            c.matdb_card.toLowerCase().includes(q),
        )
        .slice(0, 12)
        .map((c) => ({ id: c.archetype_id, label: prettyEquipmentName(c.archetype_id) }));
    } else {
      items = FLAGSHIP_EQUIPMENT.filter((f) => f.tab === activeTab).map((f) => ({ id: f.id, label: f.label }));
    }

    if (items.length === 0) {
      const empty = document.createElement("span");
      empty.className = "micro";
      empty.textContent = "No matching asset in 67-item catalog";
      pillGrid.appendChild(empty);
      return;
    }

    for (const item of items) {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "chip-btn";
      b.dataset.assetId = item.id;
      b.textContent = item.label;
      b.setAttribute("aria-pressed", String(item.id === activeId));
      b.addEventListener("click", () => {
        void selectAsset(item.id);
      });
      pillGrid.appendChild(b);
    }
  }

  for (const tabBtn of tabBtns) {
    tabBtn.addEventListener("click", () => {
      activeTab = tabBtn.dataset.eqTab as EquipmentTabId;
      for (const tb of tabBtns) tb.setAttribute("aria-pressed", String(tb.dataset.eqTab === activeTab));
      if (searchInput) searchInput.value = "";
      renderTabPills("");
      const first = FLAGSHIP_EQUIPMENT.find((f) => f.tab === activeTab);
      if (first) void selectAsset(first.id);
    });
  }
  searchInput?.addEventListener("input", () => {
    renderTabPills(searchInput.value);
  });
  catalogSelect?.addEventListener("change", () => {
    const chosen = catalogSelect.value;
    if (chosen) {
      const entry = catalog.find((c) => c.archetype_id === chosen);
      if (entry) {
        activeTab = categorizeCatalogEntry(entry);
        for (const tb of tabBtns) tb.setAttribute("aria-pressed", String(tb.dataset.eqTab === activeTab));
      }
      void selectAsset(chosen);
    }
  });
  toggleThermal?.addEventListener("click", () => {
    thermalPreview = !thermalPreview;
    toggleThermal.setAttribute("aria-pressed", String(thermalPreview));
    applySurfaceMode();
  });
  toggleWire?.addEventListener("click", () => {
    wireframePreview = !wireframePreview;
    toggleWire.setAttribute("aria-pressed", String(wireframePreview));
    applySurfaceMode();
  });
  toggleCage?.addEventListener("click", () => {
    showCage = !showCage;
    toggleCage.setAttribute("aria-pressed", String(showCage));
    cageAnchor.visible = showCage;
  });
  toggleSpin?.addEventListener("click", () => {
    autoSpin = !autoSpin;
    toggleSpin.setAttribute("aria-pressed", String(autoSpin));
  });
  toggleJoints?.addEventListener("click", () => {
    animateJoints = !animateJoints;
    toggleJoints.setAttribute("aria-pressed", String(animateJoints));
    if (!animateJoints) {
      for (const j of activeJoints) {
        j.obj.rotation.copy(j.baseRot);
        j.obj.position.copy(j.basePos);
      }
    }
  });

  // Load default asset immediately in parallel with catalog_index.json so turntable paints fast!
  renderTabPills("");
  void selectAsset(activeId);
  void (async () => {
    const resp = await fetch(`${ASSET_BASE}/equipment/catalog_index.json`);
    if (resp.ok) {
      const idx = (await resp.json()) as CatalogIndex;
      catalog = [...idx.assets];
      if (catalogSelect) {
        catalogSelect.innerHTML = "";
        for (const item of catalog) {
          const opt = document.createElement("option");
          opt.value = item.archetype_id;
          opt.textContent = `${prettyEquipmentName(item.archetype_id)} (${item.archetype_id})`;
          catalogSelect.appendChild(opt);
        }
        catalogSelect.value = activeId;
      }
      // Refresh spec metadata now that full catalog is loaded
      const cur = catalog.find((c) => c.archetype_id === activeId);
      if (cur && metaTris) metaTris.textContent = `${cur.triangle_count.toLocaleString()} tris`;
      if (cur && metaSockets) metaSockets.textContent = `${cur.socket_ids.length} sockets · ${cur.joint_count} joints`;
    }
  })();

  observeInView(stageEl, (visible) => {
    inView = visible;
  });

  let lastTs = performance.now();
  let tAnim = 0;
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    lastTs = now;
    if (!inView || document.hidden) return;
    tAnim += dt;

    if (autoSpin && !userDragging && activeModel) {
      turntable.rotation.y += dt * 0.28;
    }

    // Live articulation of real glTF door hinges, fry baskets, lids, and spinning blades!
    if (animateJoints && activeJoints.length > 0) {
      const cycle = Math.pow(0.5 - 0.5 * Math.cos(tAnim * 1.15), 1.4); // Smooth dwell at closed & open
      for (const j of activeJoints) {
        if (j.kind === "door_x") {
          j.obj.rotation.set(j.baseRot.x + cycle * 0.72, j.baseRot.y, j.baseRot.z);
        } else if (j.kind === "door_z") {
          j.obj.rotation.set(j.baseRot.x, j.baseRot.y + j.sign * cycle * 0.92, j.baseRot.z);
        } else if (j.kind === "lid_x") {
          j.obj.rotation.set(j.baseRot.x - cycle * 0.65, j.baseRot.y, j.baseRot.z);
        } else if (j.kind === "lift_y") {
          const offset = j.sign < 0 ? Math.pow(0.5 - 0.5 * Math.cos(tAnim * 1.15 + 1.5), 1.4) : cycle;
          j.obj.position.set(j.basePos.x, j.basePos.y + offset * 0.14, j.basePos.z);
        } else if (j.kind === "slide_z") {
          j.obj.position.set(j.basePos.x, j.basePos.y, j.basePos.z + cycle * 0.14);
        } else if (j.kind === "spin_x") {
          j.obj.rotation.x += dt * 5.0;
        } else if (j.kind === "spin_y") {
          j.obj.rotation.y += dt * 3.2;
        }
      }
    }

    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
}
