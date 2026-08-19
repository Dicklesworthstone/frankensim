// Imported from the owner's classic-patents.com project (2026-08-18) as the
// E2.2 articulated-airframe base: parametric build at dossier dimensions
// (FLYER_DIM matches flyer-reference.json values), canard/rudder groups
// articulated. Physics NEVER reads this mesh (role boundary,
// flyer-model-provenance-v1).
import * as THREE from "three";

/**
 * 1903 Wright Flyer Historical Engineering Specifications (US Patent 821,393)
 * Reference: Smithsonian National Air and Space Museum & 1903 Wright Drawings
 * 1 world unit = 1 metre.
 */
export const FLYER_DIM = {
  span: 12.29, // Smithsonian NASM A19610048000 — 40 ft 4 in
  chord: 1.981, // 6 ft 6 in rib chord
  gap: 1.89, // 6.2 ft sourced interplane gap (wright-brothers.org Flyer I)
  camberRatio: 0.05, // 1 in 20 as flown 1903
  thicknessRatio: 0.045, // 0.089m max wing thickness at 30% chord
  anhedralDeg: 2.37, // atan(10 in droop / half-span)
  length: 6.43, // 21 ft 1 in
  canardSpan: 3.66, // 12 ft 0 in
  canardChord: 0.76, // 2 ft 6 in
  canardGap: 0.55, // 1 ft 9.5 in
  canardArm: 2.23, // wing-to-canard, AIAA 2004-0211 / Culick
  rudderHeight: 1.35, // 4 ft 5 in
  rudderChord: 0.72, // 2 ft 4.5 in
  rudderSep: 0.72, // Distance between twin rudder fins
  rudderArm: 2.05, // Aft distance from rear wing spar
  propDiameter: 2.59, // 8 ft 6 in diameter laminated spruce propellers
  propX: 2.15, // Propeller shaft offset from centerline (port & starboard)
  ribSpacing: 0.381, // 15-inch rib pitch along span
} as const;

export interface FlyerAirframe {
  group: THREE.Group;
  upperWing: THREE.Group;
  lowerWing: THREE.Group;
  canardGroup: THREE.Group;
  rudderGroup: THREE.Group;
  cradleGroup: THREE.Group;
  leftPropBlades: THREE.Group;
  rightPropBlades: THREE.Group;
  leftBayWireMat: THREE.MeshStandardMaterial;
  rightBayWireMat: THREE.MeshStandardMaterial;
  muslinMat: THREE.MeshStandardMaterial;
  spruceMat: THREE.MeshStandardMaterial;
  textures: THREE.Texture[];
}

/**
 * Stable texture variation. Museum rendering must not depend on ambient
 * randomness: the same airframe source should produce the same material
 * pattern when a visitor replays the same control sequence.
 */
function deterministicUnit(index: number, channel: number): number {
  const sample = Math.sin((index + 1) * 12.9898 + (channel + 1) * 78.233) * 43758.5453;
  return sample - Math.floor(sample);
}

/**
 * Procedural Quarter-Sawn White Spruce Grain Texture
 */
function spruceTexture(): THREE.CanvasTexture {
  if (typeof document === "undefined") {
    return new THREE.Texture() as unknown as THREE.CanvasTexture;
  }
  const canvas = document.createElement("canvas");
  canvas.width = 512;
  canvas.height = 512;
  const ctx = canvas.getContext("2d");
  if (!ctx) return new THREE.CanvasTexture(canvas);

  // Warm honey-spruce base
  ctx.fillStyle = "#c89f68";
  ctx.fillRect(0, 0, 512, 512);

  // Fine longitudinal growth rings
  for (let i = 0; i < 96; i++) {
    const x = i * 5.3 + (deterministicUnit(i, 0) - 0.5) * 2;
    const alpha = 0.08 + (i % 6 === 0 ? 0.12 : 0.04);
    ctx.strokeStyle = `rgba(88, 48, 16, ${alpha})`;
    ctx.lineWidth = 1 + (i % 4) * 0.4;
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.bezierCurveTo(x + 12, 140, x - 10, 360, x + 8, 512);
    ctx.stroke();
  }

  // Subtle wood pores and medullary rays
  for (let j = 0; j < 300; j++) {
    const px = deterministicUnit(j, 1) * 512;
    const py = deterministicUnit(j, 2) * 512;
    ctx.fillStyle = "rgba(60, 30, 10, 0.15)";
    ctx.fillRect(px, py, 1.5, 6 + deterministicUnit(j, 3) * 12);
  }

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/**
 * Procedural Pride of the West Unbleached Muslin Linen Texture
 */
function muslinTexture(): THREE.CanvasTexture {
  if (typeof document === "undefined") {
    return new THREE.Texture() as unknown as THREE.CanvasTexture;
  }
  const canvas = document.createElement("canvas");
  canvas.width = 256;
  canvas.height = 256;
  const ctx = canvas.getContext("2d");
  if (!ctx) return new THREE.CanvasTexture(canvas);

  // Historical raw unbleached cotton-linen tone
  ctx.fillStyle = "#f3ebd7";
  ctx.fillRect(0, 0, 256, 256);

  // Micro-weave grid: warp and weft threads
  ctx.strokeStyle = "rgba(145, 120, 80, 0.16)";
  ctx.lineWidth = 0.65;
  for (let i = 0; i < 256; i += 3) {
    ctx.beginPath();
    ctx.moveTo(i, 0);
    ctx.lineTo(i, 256);
    ctx.stroke();

    ctx.beginPath();
    ctx.moveTo(0, i);
    ctx.lineTo(256, i);
    ctx.stroke();
  }

  // Irregular slub yarns in unbleached fabric
  for (let s = 0; s < 40; s++) {
    ctx.strokeStyle = "rgba(110, 85, 50, 0.22)";
    ctx.lineWidth = 1.2;
    const pos = deterministicUnit(s, 4) * 256;
    ctx.beginPath();
    if (deterministicUnit(s, 5) > 0.5) {
      ctx.moveTo(pos, deterministicUnit(s, 6) * 100);
      ctx.lineTo(pos, pos + 20 + deterministicUnit(s, 7) * 40);
    } else {
      ctx.moveTo(deterministicUnit(s, 8) * 100, pos);
      ctx.lineTo(pos + 20 + deterministicUnit(s, 9) * 40, pos);
    }
    ctx.stroke();
  }

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  tex.repeat.set(12, 4);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/**
 * Laminated Dark Walnut / Hard Ash Texture (for Propeller Blades & Skids)
 */
function walnutTexture(): THREE.CanvasTexture {
  if (typeof document === "undefined") {
    return new THREE.Texture() as unknown as THREE.CanvasTexture;
  }
  const canvas = document.createElement("canvas");
  canvas.width = 256;
  canvas.height = 256;
  const ctx = canvas.getContext("2d");
  if (!ctx) return new THREE.CanvasTexture(canvas);

  ctx.fillStyle = "#5c3317";
  ctx.fillRect(0, 0, 256, 256);

  for (let i = 0; i < 32; i++) {
    const x = i * 8;
    ctx.strokeStyle = `rgba(32, 14, 4, ${0.15 + (i % 3) * 0.1})`;
    ctx.lineWidth = 2 + (i % 2) * 1.5;
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.bezierCurveTo(x + 18, 90, x - 14, 180, x + 10, 256);
    ctx.stroke();
  }

  const tex = new THREE.CanvasTexture(canvas);
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/**
 * Wright 1903 thin parabolic camber airfoil calculation
 * s: 0.0 (leading edge) to 1.0 (trailing edge)
 */
function wrightAirfoilPoint(
  s: number,
  chord: number,
  camberRatio: number,
  halfThick: number,
  upper: boolean,
) {
  const z = chord * 0.5 - s * chord; // +z is leading edge, -z is trailing edge
  // Parabolic camber distribution: max camber at 0.38 chord
  const camber = 4 * camberRatio * s * (1 - s) * chord;
  // Authentic teardrop nose with razor thin trailing edge
  const thick = halfThick * Math.sqrt(Math.max(s, 0.008)) * (1 - s * 0.94) * 2.1;
  return {
    y: camber + (upper ? thick : -thick),
    z,
  };
}

/**
 * Lofts a realistic 3D fabric wing panel with subtle rib-ridge scalloping
 */
function loftWrightWingPanel(opts: {
  x0: number;
  x1: number;
  chord: number;
  camberRatio: number;
  thicknessRatio: number;
  anhedralRad: number;
  worldX0?: number;
  stations?: number;
  airfoilPts?: number;
  ribSpacing?: number;
}): THREE.BufferGeometry {
  const stations = opts.stations ?? 36;
  const n = opts.airfoilPts ?? 24;
  const halfT = opts.thicknessRatio * opts.chord * 0.5;
  const ribPitch = opts.ribSpacing ?? 0.381;
  const positions: number[] = [];
  const uvs: number[] = [];
  const indices: number[] = [];

  for (let i = 0; i <= stations; i++) {
    const t = i / stations;
    const x = opts.x0 + (opts.x1 - opts.x0) * t;
    const worldX = (opts.worldX0 ?? 0) + x;
    const droop = -Math.tan(opts.anhedralRad) * Math.abs(worldX);

    // Subtle transverse fabric sag between ribs (0.38m pitch)
    const ribPhase = (worldX % ribPitch) / ribPitch;
    const ribBulge = Math.sin(ribPhase * Math.PI) * 0.0038;

    // Upper airfoil contour (LE to TE)
    for (let k = 0; k <= n; k++) {
      const s = k / n;
      const u = wrightAirfoilPoint(s, opts.chord, opts.camberRatio, halfT, true);
      const sag = s > 0.05 && s < 0.95 ? -ribBulge : 0;
      positions.push(x, u.y + droop + sag, u.z);
      uvs.push(t, s * 0.5);
    }
    // Lower airfoil contour (TE to LE)
    for (let k = n; k >= 0; k--) {
      const s = k / n;
      const l = wrightAirfoilPoint(s, opts.chord, opts.camberRatio, halfT, false);
      const sag = s > 0.05 && s < 0.95 ? ribBulge * 0.6 : 0;
      positions.push(x, l.y + droop + sag, l.z);
      uvs.push(t, 1 - s * 0.5);
    }
  }

  const ringVerts = (n + 1) * 2;
  const flipWinding = opts.x1 < opts.x0;

  for (let i = 0; i < stations; i++) {
    const a = i * ringVerts;
    const b = (i + 1) * ringVerts;
    for (let k = 0; k < ringVerts - 1; k++) {
      if (flipWinding) {
        indices.push(a + k, a + k + 1, b + k);
        indices.push(b + k, a + k + 1, b + k + 1);
      } else {
        indices.push(a + k, b + k, a + k + 1);
        indices.push(b + k, b + k + 1, a + k + 1);
      }
    }
  }

  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
  geo.setAttribute("uv", new THREE.Float32BufferAttribute(uvs, 2));
  geo.setIndex(indices);
  geo.computeVertexNormals();
  return geo;
}

/**
 * Creates an authentic spruce interplane strut with streamlined cross section
 * and brass socket hardware at both ends.
 */
function createStreamlineStrut(
  length: number,
  spruceMat: THREE.Material,
  brassMat: THREE.Material,
): THREE.Group {
  const g = new THREE.Group();

  // Aerodynamic teardrop strut body
  const bodyGeo = new THREE.CylinderGeometry(0.024, 0.024, length - 0.06, 12);
  bodyGeo.scale(0.65, 1.0, 1.45); // 1.25" x 2.5" streamlined aspect ratio
  const body = new THREE.Mesh(bodyGeo, spruceMat);
  body.castShadow = true;
  g.add(body);

  // Top and bottom brass mounting sockets
  for (const y of [(length - 0.04) / 2, -(length - 0.04) / 2]) {
    const socket = new THREE.Mesh(new THREE.CylinderGeometry(0.026, 0.026, 0.04, 10), brassMat);
    socket.position.y = y;
    g.add(socket);

    // Socket retaining bolt pin
    const pin = new THREE.Mesh(new THREE.CylinderGeometry(0.005, 0.005, 0.06, 6), brassMat);
    pin.rotation.z = Math.PI / 2;
    pin.position.y = y;
    g.add(pin);
  }

  return g;
}

/**
 * Creates high-strength hard steel bracing wire with authentic brass turnbuckle
 */
function addRiggingWire(
  parent: THREE.Group,
  p1: [number, number, number],
  p2: [number, number, number],
  steelMat: THREE.Material,
  brassMat: THREE.Material,
) {
  const v1 = new THREE.Vector3(...p1);
  const v2 = new THREE.Vector3(...p2);
  const dist = v1.distanceTo(v2);
  if (dist < 0.01) return;

  const mid = new THREE.Vector3().addVectors(v1, v2).multiplyScalar(0.5);
  const dir = new THREE.Vector3().subVectors(v2, v1).normalize();

  // Steel wire
  const wireGeo = new THREE.CylinderGeometry(0.0045, 0.0045, dist, 6);
  const wire = new THREE.Mesh(wireGeo, steelMat);
  wire.position.copy(mid);
  wire.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
  parent.add(wire);

  // Miniature brass turnbuckle near lower anchor
  const turnbucklePos = new THREE.Vector3().lerpVectors(v1, v2, 0.18);
  const tbGeo = new THREE.CylinderGeometry(0.011, 0.011, 0.065, 8);
  const turnbuckle = new THREE.Mesh(tbGeo, brassMat);
  turnbuckle.position.copy(turnbucklePos);
  turnbuckle.quaternion.copy(wire.quaternion);
  parent.add(turnbuckle);
}

/**
 * Creates an authentic Wright 1903 laminated spruce 8ft 6in propeller
 */
function createWrightPropeller(
  diameter: number,
  woodMat: THREE.Material,
  ironMat: THREE.Material,
  brassMat: THREE.Material,
): THREE.Group {
  const propGroup = new THREE.Group();
  const radius = diameter / 2;

  // Center hub flange & bolt circle
  const hub = new THREE.Mesh(new THREE.CylinderGeometry(0.075, 0.075, 0.12, 16), ironMat);
  hub.rotation.x = Math.PI / 2;
  propGroup.add(hub);

  // 6 Hub clamping bolts
  for (let b = 0; b < 6; b++) {
    const angle = (b * Math.PI) / 3;
    const bolt = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.14, 6), brassMat);
    bolt.rotation.x = Math.PI / 2;
    bolt.position.set(Math.cos(angle) * 0.052, Math.sin(angle) * 0.052, 0);
    propGroup.add(bolt);
  }

  // Driven roller sprocket on shaft
  const sprocket = new THREE.Mesh(new THREE.CylinderGeometry(0.105, 0.105, 0.024, 20), ironMat);
  sprocket.rotation.x = Math.PI / 2;
  sprocket.position.z = 0.075;
  propGroup.add(sprocket);

  // 2 Laminated spruce scimitar blades
  for (let bladeIdx = 0; bladeIdx < 2; bladeIdx++) {
    const bladeAngle = bladeIdx * Math.PI;
    const stations = 18;
    const ringPts = 10;
    const positions: number[] = [];
    const indices: number[] = [];

    for (let i = 0; i <= stations; i++) {
      const r = 0.08 + (i / stations) * (radius - 0.08);
      const rn = (r - 0.08) / (radius - 0.08);

      // Authentic Wright pitch distribution (twist from 30° at root to 14° at 75% to 8° at tip)
      const twist = (1 - rn * 0.8) * 0.58 + 0.14;
      const chord = 0.22 * Math.sin(Math.PI * Math.max(0.08, rn)) * (1.1 - rn * 0.35);
      const thick = 0.022 * (1 - rn * 0.75);
      const sweep = rn * rn * 0.095; // Wright scimitar trailing edge sweep

      for (let k = 0; k <= ringPts; k++) {
        const theta = (k / ringPts) * Math.PI * 2;
        const localX = Math.cos(theta) * (chord * 0.5);
        const localZ = Math.sin(theta) * thick;

        // Apply pitch twist rotation
        const rotX = localX * Math.cos(twist) - localZ * Math.sin(twist);
        const rotZ = localX * Math.sin(twist) + localZ * Math.cos(twist) + sweep;

        positions.push(rotX, r, rotZ);
      }
    }

    for (let i = 0; i < stations; i++) {
      for (let k = 0; k < ringPts; k++) {
        const a = i * (ringPts + 1) + k;
        const b = a + ringPts + 1;
        indices.push(a, b, a + 1);
        indices.push(b, b + 1, a + 1);
      }
    }

    const bladeGeo = new THREE.BufferGeometry();
    bladeGeo.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    bladeGeo.setIndex(indices);
    bladeGeo.computeVertexNormals();

    const bladeMesh = new THREE.Mesh(bladeGeo, woodMat);
    bladeMesh.rotation.z = bladeAngle;
    bladeMesh.castShadow = true;
    propGroup.add(bladeMesh);
  }

  return propGroup;
}

/**
 * Builds the museum-grade 1903 Wright Flyer 3D Airframe
 */
export function buildWrightFlyerAirframe(): FlyerAirframe {
  const d = FLYER_DIM;
  const textures: THREE.Texture[] = [];

  const spruceMap = spruceTexture();
  const muslinMap = muslinTexture();
  const walnutMap = walnutTexture();
  textures.push(spruceMap, muslinMap, walnutMap);

  // Materials with authentic 1903 archival rendering properties
  const muslinMat = new THREE.MeshPhysicalMaterial({
    map: muslinMap,
    color: 0xf3ead5,
    roughness: 0.84,
    metalness: 0.0,
    transmission: 0.08,
    thickness: 0.025,
    transparent: true,
    opacity: 0.96,
    side: THREE.DoubleSide,
  });

  const spruceMat = new THREE.MeshStandardMaterial({
    map: spruceMap,
    color: 0xd2a76e,
    roughness: 0.44,
    metalness: 0.02,
  });

  const ashMat = new THREE.MeshStandardMaterial({
    map: walnutMap,
    color: 0x8a582b,
    roughness: 0.55,
    metalness: 0.03,
  });

  const walnutMat = new THREE.MeshStandardMaterial({
    map: walnutMap,
    color: 0x5a3118,
    roughness: 0.38,
    metalness: 0.06,
  });

  const steelMat = new THREE.MeshStandardMaterial({
    color: 0x94a3b8,
    metalness: 0.95,
    roughness: 0.18,
  });
  const leftBayWireMat = steelMat.clone();
  const rightBayWireMat = steelMat.clone();

  const brassMat = new THREE.MeshStandardMaterial({
    color: 0xd4af37,
    metalness: 0.88,
    roughness: 0.28,
  });

  const alumMat = new THREE.MeshStandardMaterial({
    color: 0xa1a1aa,
    metalness: 0.86,
    roughness: 0.35,
  });

  const ironMat = new THREE.MeshStandardMaterial({
    color: 0x334155,
    metalness: 0.78,
    roughness: 0.48,
  });

  const copperMat = new THREE.MeshStandardMaterial({
    color: 0xc25e1a,
    metalness: 0.85,
    roughness: 0.25,
  });

  const woolMat = new THREE.MeshStandardMaterial({
    color: 0x1e293b,
    roughness: 0.92,
  });

  const skinMat = new THREE.MeshStandardMaterial({
    color: 0xd4a373,
    roughness: 0.65,
  });

  const group = new THREE.Group();
  const anhedralRad = (d.anhedralDeg * Math.PI) / 180;
  const yUpper = d.gap / 2;
  const yLower = -d.gap / 2;
  const zFront = d.chord * 0.42;
  const zRear = -d.chord * 0.38;
  const centerHalf = d.span * 0.22;
  const tipInboard = centerHalf;

  // --- 1. UPPER & LOWER BIPLANE WINGS ---
  const makeWingAssembly = (yPos: number) => {
    const wingGroup = new THREE.Group();
    wingGroup.position.y = yPos;

    // Center Fixed Section (Spanning between engine/pilot cradle)
    const centerGeo = loftWrightWingPanel({
      x0: -centerHalf,
      x1: centerHalf,
      chord: d.chord,
      camberRatio: d.camberRatio,
      thicknessRatio: d.thicknessRatio,
      anhedralRad: anhedralRad,
      stations: 24,
      airfoilPts: 20,
    });
    const centerMesh = new THREE.Mesh(centerGeo, muslinMat);
    centerMesh.castShadow = true;
    centerMesh.receiveShadow = true;
    wingGroup.add(centerMesh);

    // Add Internal Spars & Cambered Rib Trusses
    const addRibsAndSpars = (parent: THREE.Group, x0: number, x1: number, worldX0: number) => {
      const spanLen = Math.abs(x1 - x0);
      const ribCount = Math.max(4, Math.round(spanLen / d.ribSpacing));

      for (let i = 0; i <= ribCount; i++) {
        const localX = x0 + (i / ribCount) * (x1 - x0);
        const droop = -Math.tan(anhedralRad) * Math.abs(worldX0 + localX);

        // Cambered Rib Truss Capstrip
        const ribCurvePts: THREE.Vector3[] = [];
        for (let r = 0; r <= 16; r++) {
          const s = r / 16;
          const pt = wrightAirfoilPoint(s, d.chord, d.camberRatio, 0.015, true);
          ribCurvePts.push(new THREE.Vector3(localX, pt.y + droop, pt.z));
        }
        const ribCurve = new THREE.CatmullRomCurve3(ribCurvePts);
        const ribGeo = new THREE.TubeGeometry(ribCurve, 14, 0.009, 6, false);
        const ribMesh = new THREE.Mesh(ribGeo, spruceMat);
        ribMesh.castShadow = true;
        parent.add(ribMesh);
      }

      // Front Spar (Spruce Beam: 1.25" x 1.75")
      const sparMidX = (x0 + x1) / 2;
      const frontSpar = new THREE.Mesh(
        new THREE.CylinderGeometry(0.022, 0.022, spanLen, 8),
        spruceMat,
      );
      frontSpar.rotation.z = Math.PI / 2;
      frontSpar.position.set(sparMidX, 0.016, zFront);
      parent.add(frontSpar);

      // Rear Spar (Spruce Beam: 1" x 1.5")
      const rearSpar = new THREE.Mesh(
        new THREE.CylinderGeometry(0.019, 0.019, spanLen, 8),
        spruceMat,
      );
      rearSpar.rotation.z = Math.PI / 2;
      rearSpar.position.set(sparMidX, 0.012, zRear);
      parent.add(rearSpar);
    };

    addRibsAndSpars(wingGroup, -centerHalf, centerHalf, 0);

    // Flexible Wing Tips (Hinged for Wing Warping Deflection)
    const addFlexibleTip = (sign: -1 | 1) => {
      const tipGroup = new THREE.Group();
      tipGroup.position.x = sign * tipInboard;
      const tipSpan = sign * (d.span / 2 - tipInboard);

      const tipGeo = loftWrightWingPanel({
        x0: 0,
        x1: tipSpan,
        chord: d.chord,
        camberRatio: d.camberRatio,
        thicknessRatio: d.thicknessRatio,
        anhedralRad: anhedralRad,
        worldX0: sign * tipInboard,
        stations: 20,
        airfoilPts: 20,
      });
      const tipMesh = new THREE.Mesh(tipGeo, muslinMat);
      tipMesh.castShadow = true;
      tipMesh.receiveShadow = true;
      tipGroup.add(tipMesh);

      addRibsAndSpars(tipGroup, 0, tipSpan, sign * tipInboard);

      // Bent Wooden Wingtip Bow
      const bow = new THREE.Mesh(
        new THREE.TorusGeometry(d.chord * 0.46, 0.018, 8, 20, Math.PI),
        ashMat,
      );
      bow.rotation.x = Math.PI / 2;
      bow.rotation.z = sign > 0 ? -Math.PI / 2 : Math.PI / 2;
      bow.position.set(tipSpan, -Math.tan(anhedralRad) * (d.span / 2), 0);
      tipGroup.add(bow);

      tipGroup.name = sign < 0 ? "leftTip" : "rightTip";
      wingGroup.add(tipGroup);
    };

    addFlexibleTip(-1);
    addFlexibleTip(1);
    return wingGroup;
  };

  const upperWing = makeWingAssembly(yUpper);
  const lowerWing = makeWingAssembly(yLower);
  group.add(upperWing, lowerWing);

  // --- 2. INTERPLANE STRUTS & TRUSS BRACING RIGGING ---
  const strutBayFractions = [-0.48, -0.3, -0.12, 0.12, 0.3, 0.48];
  const strutXs = strutBayFractions.map((f) => f * d.span);

  for (const x of strutXs) {
    for (const z of [zFront, zRear]) {
      const strut = createStreamlineStrut(d.gap, spruceMat, brassMat);
      strut.position.set(x, 0, z);
      group.add(strut);
    }
  }

  // Cross-bracing piano wire network with turnbuckles
  for (let i = 0; i < strutXs.length - 1; i++) {
    const xA = strutXs[i]!;
    const xB = strutXs[i + 1]!;

    const bayMat = (xA + xB) / 2 < 0 ? leftBayWireMat : rightBayWireMat;
    // Front bay X-wires
    addRiggingWire(group, [xA, yUpper, zFront], [xB, yLower, zFront], bayMat, brassMat);
    addRiggingWire(group, [xA, yLower, zFront], [xB, yUpper, zFront], bayMat, brassMat);

    // Rear bay X-wires
    addRiggingWire(group, [xA, yUpper, zRear], [xB, yLower, zRear], bayMat, brassMat);
    addRiggingWire(group, [xA, yLower, zRear], [xB, yUpper, zRear], bayMat, brassMat);

    // Fore-and-aft bay cross-wires
    addRiggingWire(group, [xA, yUpper, zFront], [xA, yLower, zRear], steelMat, brassMat);
    addRiggingWire(group, [xA, yLower, zFront], [xA, yUpper, zRear], steelMat, brassMat);
  }

  // --- 3. LANDING SKIDS & ASH RUNNERS ---
  const createLandingSkid = (xPos: number) => {
    const skidGroup = new THREE.Group();
    // Curved ash runner extending forward from under the wings to support the canard
    const runnerCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(xPos, yLower - 0.06, zFront + 1.45),
      new THREE.Vector3(xPos, yLower - 0.38, zFront + 0.45),
      new THREE.Vector3(xPos, yLower - 0.48, 0),
      new THREE.Vector3(xPos, yLower - 0.46, zRear - 0.25),
      new THREE.Vector3(xPos, yLower - 0.35, zRear - 0.65),
    ]);
    const runner = new THREE.Mesh(
      new THREE.TubeGeometry(runnerCurve, 32, 0.032, 10, false),
      ashMat,
    );
    runner.castShadow = true;
    skidGroup.add(runner);

    // Vertical skid uprights from lower wing front and rear spars
    for (const z of [zFront, zRear]) {
      const post = new THREE.Mesh(new THREE.CylinderGeometry(0.022, 0.022, 0.48, 8), spruceMat);
      post.position.set(xPos, yLower - 0.24, z);
      skidGroup.add(post);
    }

    // Skid tip forward diagonal stay wire
    addRiggingWire(
      skidGroup,
      [xPos, yLower, zFront],
      [xPos, yLower - 0.06, zFront + 1.45],
      steelMat,
      brassMat,
    );
    return skidGroup;
  };
  group.add(createLandingSkid(-0.95), createLandingSkid(0.95));

  // --- 4. FORWARD BIPLANE CANARD ELEVATOR (PITCH CONTROL) ---
  const canardGroup = new THREE.Group();
  canardGroup.position.set(0, -0.05, d.chord / 2 + d.canardArm);

  const createCanardPlane = (yOffset: number) => {
    const geo = loftWrightWingPanel({
      x0: -d.canardSpan / 2,
      x1: d.canardSpan / 2,
      chord: d.canardChord,
      camberRatio: 0.045,
      thicknessRatio: 0.05,
      anhedralRad: 0,
      stations: 16,
      airfoilPts: 16,
    });
    const mesh = new THREE.Mesh(geo, muslinMat);
    mesh.position.y = yOffset;
    mesh.castShadow = true;
    return mesh;
  };

  canardGroup.add(createCanardPlane(d.canardGap / 2), createCanardPlane(-d.canardGap / 2));

  // Canard vertical struts
  for (const x of [-d.canardSpan * 0.4, 0, d.canardSpan * 0.4]) {
    const s = createStreamlineStrut(d.canardGap, spruceMat, brassMat);
    s.position.set(x, 0, 0);
    canardGroup.add(s);
  }

  // Forward canard outrigger boom framework
  for (const x of [-0.95, 0.95]) {
    for (const y of [yLower + 0.08, 0.16]) {
      const boom = new THREE.Mesh(
        new THREE.CylinderGeometry(0.02, 0.02, d.canardArm + 0.4, 8),
        spruceMat,
      );
      boom.rotation.x = Math.PI / 2;
      boom.position.set(x, y, d.chord / 2 + d.canardArm * 0.48);
      group.add(boom);
    }
  }
  group.add(canardGroup);

  // --- 5. AFT TWIN VERTICAL RUDDERS (YAW CONTROL & CLAIM 1 COUPLING) ---
  const rudderGroup = new THREE.Group();
  rudderGroup.position.set(0, 0.05, -d.chord / 2 - d.rudderArm);

  const createRudderFin = (xOffset: number) => {
    const finGroup = new THREE.Group();
    finGroup.position.x = xOffset;

    // Muslin covered double fin
    const fabric = new THREE.Mesh(
      new THREE.PlaneGeometry(d.rudderChord, d.rudderHeight),
      muslinMat,
    );
    fabric.rotation.y = Math.PI / 2;
    fabric.castShadow = true;
    finGroup.add(fabric);

    // Spruce framing around rudder edges
    const frame = new THREE.Mesh(
      new THREE.BoxGeometry(0.026, d.rudderHeight + 0.02, d.rudderChord + 0.02),
      spruceMat,
    );
    finGroup.add(frame);

    // Brass hinge pintles
    for (const hY of [-0.4, 0.4]) {
      const hinge = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.04, 6), brassMat);
      hinge.position.set(0, hY, -d.rudderChord / 2);
      finGroup.add(hinge);
    }
    return finGroup;
  };

  rudderGroup.add(createRudderFin(-d.rudderSep / 2), createRudderFin(d.rudderSep / 2));

  // Aft rudder outrigger booms
  for (const x of [-d.rudderSep / 2, d.rudderSep / 2]) {
    const boom = new THREE.Mesh(
      new THREE.CylinderGeometry(0.02, 0.02, d.rudderArm + 0.35, 8),
      spruceMat,
    );
    boom.rotation.x = Math.PI / 2;
    boom.position.set(x, 0, -d.chord / 2 - d.rudderArm * 0.48);
    group.add(boom);
  }
  group.add(rudderGroup);

  // --- 6. CHARLIE TAYLOR 12-HP 4-CYLINDER INLINE ENGINE ---
  const engine = new THREE.Group();
  engine.position.set(0.48, yLower + 0.24, 0.06);

  // Cast 8%-copper aluminum alloy crankcase
  const crankcase = new THREE.Mesh(new THREE.BoxGeometry(0.42, 0.26, 0.94), alumMat);
  crankcase.castShadow = true;
  engine.add(crankcase);

  // 4 Cast-iron horizontal cylinder barrels
  for (let i = 0; i < 4; i++) {
    const cylZ = -0.34 + i * 0.23;
    const cyl = new THREE.Mesh(new THREE.CylinderGeometry(0.076, 0.076, 0.32, 16), ironMat);
    cyl.rotation.z = Math.PI / 2;
    cyl.position.set(0.26, 0.02, cylZ);
    cyl.castShadow = true;
    engine.add(cyl);

    // Valve rocker arms and overhead springs
    const valve = new THREE.Mesh(new THREE.CylinderGeometry(0.022, 0.022, 0.09, 8), brassMat);
    valve.position.set(0.4, 0.09, cylZ);
    engine.add(valve);

    // Spark igniter trip rods
    const igniter = new THREE.Mesh(new THREE.CylinderGeometry(0.006, 0.006, 0.08, 6), steelMat);
    igniter.position.set(0.26, 0.12, cylZ);
    engine.add(igniter);
  }

  // Heavy 28-lb cast iron flywheel
  const flywheel = new THREE.Mesh(new THREE.CylinderGeometry(0.21, 0.21, 0.05, 28), ironMat);
  flywheel.rotation.x = Math.PI / 2;
  flywheel.position.set(-0.18, 0, -0.42);
  flywheel.castShadow = true;
  engine.add(flywheel);

  // Vertical 66-tube copper radiator mounted on front strut
  const radiatorHeight = d.gap * 0.62;
  const radiatorLocalY = yUpper - (yLower + 0.24) - radiatorHeight * 0.2;
  for (let i = 0; i < 10; i++) {
    const tube = new THREE.Mesh(
      new THREE.CylinderGeometry(0.009, 0.009, radiatorHeight, 6),
      copperMat,
    );
    tube.position.set(-0.12, radiatorLocalY, zFront - 0.02 + (i - 4.5) * 0.024);
    engine.add(tube);
  }

  // Brass gravity fuel tank mounted high on upper center strut
  const fuelTank = new THREE.Mesh(new THREE.CylinderGeometry(0.07, 0.07, 0.38, 14), brassMat);
  fuelTank.rotation.z = Math.PI / 2;
  fuelTank.position.set(0.05, yUpper - (yLower + 0.24) - 0.15, 0.06);
  engine.add(fuelTank);

  // Copper fuel line running down to carburetor
  const fuelLine = new THREE.Mesh(new THREE.CylinderGeometry(0.004, 0.004, 0.85, 6), copperMat);
  fuelLine.position.set(0.05, yUpper - (yLower + 0.24) - 0.6, 0.06);
  engine.add(fuelLine);

  group.add(engine);

  // --- 7. PRONE PILOT HIP CRADLE & ORVILLE WRIGHT FIGURE ---
  const cradleGroup = new THREE.Group();
  cradleGroup.position.set(-0.35, yLower + 0.08, 0.06);

  // Ash guide rails on lower wing surface
  const railL = new THREE.Mesh(new THREE.CylinderGeometry(0.012, 0.012, 1.05, 8), ashMat);
  railL.rotation.x = Math.PI / 2;
  railL.position.set(-0.26, 0.01, 0);
  cradleGroup.add(railL);

  const railR = new THREE.Mesh(new THREE.CylinderGeometry(0.012, 0.012, 1.05, 8), ashMat);
  railR.rotation.x = Math.PI / 2;
  railR.position.set(0.26, 0.01, 0);
  cradleGroup.add(railR);

  // Sliding wooden hip cradle
  const cradleBox = new THREE.Mesh(new THREE.BoxGeometry(0.48, 0.09, 0.35), ashMat);
  cradleBox.position.set(0, 0.06, 0.05);
  cradleGroup.add(cradleBox);

  // Warping cable attachment eyes on cradle flanks
  for (const xEye of [-0.25, 0.25]) {
    const eye = new THREE.Mesh(new THREE.TorusGeometry(0.022, 0.006, 6, 8), brassMat);
    eye.position.set(xEye, 0.09, 0.05);
    cradleGroup.add(eye);
  }

  // Prone pilot figure (Orville Wright)
  const pilotTorso = new THREE.Mesh(new THREE.CapsuleGeometry(0.13, 0.58, 8, 10), woolMat);
  pilotTorso.rotation.x = Math.PI / 2;
  pilotTorso.position.set(0, 0.16, -0.06);
  pilotTorso.castShadow = true;
  cradleGroup.add(pilotTorso);

  const pilotHead = new THREE.Mesh(new THREE.SphereGeometry(0.095, 12, 10), skinMat);
  pilotHead.position.set(0, 0.22, 0.44);
  cradleGroup.add(pilotHead);

  const pilotCap = new THREE.Mesh(new THREE.CylinderGeometry(0.105, 0.105, 0.045, 12), woolMat);
  pilotCap.position.set(0, 0.28, 0.44);
  cradleGroup.add(pilotCap);

  // Left hand elevator control lever
  const elevatorLever = new THREE.Mesh(new THREE.CylinderGeometry(0.014, 0.014, 0.52, 8), ashMat);
  elevatorLever.position.set(-0.28, 0.26, 0.26);
  elevatorLever.rotation.z = -0.22;
  cradleGroup.add(elevatorLever);

  group.add(cradleGroup);

  // --- 8. ANEMOMETER & FLIGHT INSTRUMENTS ---
  const instCluster = new THREE.Group();
  instCluster.position.set(-0.12, yLower + 0.58, zFront);

  const instMount = new THREE.Mesh(new THREE.BoxGeometry(0.11, 0.15, 0.035), ashMat);
  instCluster.add(instMount);

  // Richard anemometer cups wheel
  const anemometerHub = new THREE.Mesh(
    new THREE.CylinderGeometry(0.016, 0.016, 0.045, 8),
    brassMat,
  );
  anemometerHub.rotation.x = Math.PI / 2;
  anemometerHub.position.set(0, 0.045, 0.035);
  instCluster.add(anemometerHub);

  for (let c = 0; c < 4; c++) {
    const cupArm = new THREE.Mesh(new THREE.CylinderGeometry(0.0035, 0.0035, 0.07, 6), steelMat);
    cupArm.rotation.z = (c * Math.PI) / 2;
    cupArm.position.set(0, 0.045, 0.05);
    instCluster.add(cupArm);
  }

  // Stopwatch dial
  const stopwatch = new THREE.Mesh(new THREE.CylinderGeometry(0.026, 0.026, 0.016, 14), brassMat);
  stopwatch.rotation.x = Math.PI / 2;
  stopwatch.position.set(0, -0.04, 0.025);
  instCluster.add(stopwatch);

  group.add(instCluster);

  // --- 9. TWIN COUNTER-ROTATING PUSHER PROPELLERS & CHAIN CASINGS ---
  const makePropAssembly = (xPos: number, isPortCrossed: boolean) => {
    const propMount = new THREE.Group();
    propMount.position.set(xPos, 0.02, zRear - 0.22);

    const prop = createWrightPropeller(d.propDiameter, walnutMat, ironMat, brassMat);
    propMount.add(prop);

    // Tubular steel chain drive casing running from engine to propeller shaft
    const driveDist = Math.abs(xPos - 0.48);
    const chainCasing = new THREE.Mesh(
      new THREE.CylinderGeometry(0.014, 0.014, driveDist, 8),
      steelMat,
    );
    chainCasing.rotation.z = Math.PI / 2;
    chainCasing.position.set(-(xPos - 0.48) / 2, -0.1, 0.07);
    if (isPortCrossed) {
      chainCasing.rotation.y = 0.045; // Port chain crossed in figure-8 guide tube
    }
    propMount.add(chainCasing);

    group.add(propMount);
    return prop;
  };

  const leftPropBlades = makePropAssembly(-d.propX, true);
  const rightPropBlades = makePropAssembly(d.propX, false);

  // --- 10. CLAIM 1 WING-WARPING TO RUDDER CROSS-CABLES ---
  const outerLeftX = strutXs[0]!;
  const outerRightX = strutXs[strutXs.length - 1]!;

  addRiggingWire(
    group,
    [outerLeftX, yLower, zRear],
    [-d.rudderSep / 2, 0.05, -d.chord / 2 - d.rudderArm],
    steelMat,
    brassMat,
  );
  addRiggingWire(
    group,
    [outerRightX, yLower, zRear],
    [d.rudderSep / 2, 0.05, -d.chord / 2 - d.rudderArm],
    steelMat,
    brassMat,
  );

  return {
    group,
    upperWing,
    lowerWing,
    canardGroup,
    rudderGroup,
    cradleGroup,
    leftPropBlades,
    rightPropBlades,
    leftBayWireMat,
    rightBayWireMat,
    muslinMat,
    spruceMat,
    textures,
  };
}

/**
 * Updates Wright Flyer propeller rotation, wing warp deflection, elevator canard pitch, rudder yaw, wire tension colors, and fabric cutaway.
 */
export function updateWrightFlyerKinematics(
  airframe: FlyerAirframe,
  delta: number,
  wingWarpDeg: number,
  rudderYawDeg: number,
  elevatorPitchDeg: number,
  propDisplayOmegaRadPerS: number,
  cradleStudioX: number,
  leftBayTension: number,
  rightBayTension: number,
  isCutaway = false,
): void {
  // Propellers Rotation (Counter-Rotating to eliminate gyroscopic torque)
  airframe.leftPropBlades.rotation.z += propDisplayOmegaRadPerS * delta;
  airframe.rightPropBlades.rotation.z -= propDisplayOmegaRadPerS * delta;

  // Animate Wing Warping Deflection on Mesh Tips
  const warpRad = (wingWarpDeg * Math.PI) / 180;
  const leftTipUpper = airframe.upperWing.getObjectByName("leftTip");
  const rightTipUpper = airframe.upperWing.getObjectByName("rightTip");
  const leftTipLower = airframe.lowerWing.getObjectByName("leftTip");
  const rightTipLower = airframe.lowerWing.getObjectByName("rightTip");

  if (leftTipUpper && rightTipUpper && leftTipLower && rightTipLower) {
    leftTipUpper.rotation.x = warpRad * 0.6;
    leftTipLower.rotation.x = warpRad * 0.6;
    rightTipUpper.rotation.x = -warpRad * 0.6;
    rightTipLower.rotation.x = -warpRad * 0.6;
  }

  // Pilot hip cradle sliding sideways during wing warping
  airframe.cradleGroup.position.x = cradleStudioX;

  // Animate Elevator & Rudder
  airframe.canardGroup.rotation.x = (-elevatorPitchDeg * Math.PI) / 180;
  airframe.rudderGroup.rotation.y = (-rudderYawDeg * Math.PI) / 180;

  // Interplane X-wires: the high-AoA tip carries extra lift, so that bay's
  // piano wire goes amber, then red. Slack bay stays steel-grey.
  const leftTension = leftBayTension;
  const rightTension = rightBayTension;
  const paintBay = (mat: THREE.MeshStandardMaterial, tension: number) => {
    if (tension > 1.15) mat.color.setHex(0xef4444);
    else if (tension > 0.55) mat.color.setHex(0xf59e0b);
    else mat.color.setHex(0x94a3b8);
  };
  paintBay(airframe.leftBayWireMat, leftTension);
  paintBay(airframe.rightBayWireMat, rightTension);

  // Cutaway muslin transparency for inspecting internal rib truss and control cables
  airframe.muslinMat.opacity = isCutaway ? 0.35 : 1.0;
  airframe.muslinMat.transparent = isCutaway;
}
