// Wright-brother figure meshes (bead guzez.14). THIN consumer of the
// tested joint math in figure.ts: lofted tapered limbs with real
// elbow/knee hinges, 1903 sack suits, flat caps, and PUBLIC-DOMAIN
// Library-of-Congress portrait faces (src/assets, ~12 KB each) front-
// projected onto the head — Orville's mustache ships in the photograph.
// The photos are the only binary assets in the app; provenance:
// Wikimedia Commons "Wilbur Wright-crop.jpg" / "Orville Wright
// 1905-crop.jpg", both tagged Public domain (pre-1928 publication).

import * as THREE from "three";
import {
  BINOCULAR_POSE,
  PRONE_POSE,
  armAimAngles,
  figureSpec,
  gaitPose,
  type Brother,
  type FigureSpec,
} from "./figure.ts";
import wilburFaceUrl from "./assets/wilbur-face.jpg";
import orvilleFaceUrl from "./assets/orville-face.jpg";

const SUIT = new THREE.MeshStandardMaterial({ color: 0x2b2a31, roughness: 0.92 });
const SUIT_DARK = new THREE.MeshStandardMaterial({ color: 0x1f1e24, roughness: 0.95 });
const SHIRT = new THREE.MeshStandardMaterial({ color: 0xe9e4d6, roughness: 0.85 });
const SKIN = new THREE.MeshStandardMaterial({ color: 0xc9a284, roughness: 0.75 });
const IRON = new THREE.MeshStandardMaterial({ color: 0x3a3f45, roughness: 0.6, metalness: 0.5 });

/** Head material: painted hair/skin immediately, the PD portrait
 * drawn into the front band when the image decodes (async; headless
 * and offline runs keep the painted fallback — never a blank face). */
function headMaterial(brother: Brother): THREE.MeshStandardMaterial {
  const canvas = document.createElement("canvas");
  canvas.width = 256;
  canvas.height = 128;
  const ctx = canvas.getContext("2d");
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  const mat = new THREE.MeshStandardMaterial({ map: tex, roughness: 0.8 });
  if (ctx === null) {
    return mat;
  }
  const paint = (img: HTMLImageElement | null): void => {
    const w = canvas.width;
    const h = canvas.height;
    // Hair crown + back, skin below the hat line.
    ctx.fillStyle = brother === "wilbur" ? "#4a3c30" : "#33291f";
    ctx.fillRect(0, 0, w, h);
    ctx.fillStyle = "#c9a284";
    ctx.fillRect(0, h * 0.3, w, h * 0.7);
    if (img !== null) {
      // Front band (sphere u=0.5 faces +x): the portrait, upright.
      ctx.drawImage(img, w * 0.3, h * 0.06, w * 0.4, h * 0.94);
    }
    tex.needsUpdate = true;
  };
  paint(null);
  const img = new Image();
  img.onload = (): void => paint(img);
  img.src = brother === "wilbur" ? wilburFaceUrl : orvilleFaceUrl;
  return mat;
}


/** Every figure mesh casts AND receives in the scene's shadow map
 * (the sun's shadow pass is presentation-owned by flyerScene): a
 * figure that casts nothing reads as floating, and one that receives
 * nothing stays lit under the wing that shades him. */
function castAll(root: THREE.Object3D): void {
  root.traverse((obj) => {
    if ((obj as THREE.Mesh).isMesh) {
      obj.castShadow = true;
      obj.receiveShadow = true;
    }
  });
}

function taperedLimb(
  topR: number,
  botR: number,
  len: number,
  mat: THREE.Material,
): THREE.Mesh {
  const geo = new THREE.CylinderGeometry(topR, botR, len, 10);
  geo.translate(0, -len / 2, 0); // pivot at the TOP joint
  return new THREE.Mesh(geo, mat);
}

interface Limb {
  root: THREE.Group; // hip/shoulder hinge
  mid: THREE.Group; // knee/elbow hinge
}

function buildLimb(
  rootPos: readonly [number, number, number],
  upper: THREE.Mesh,
  upperLen: number,
  lower: THREE.Mesh,
  tip: THREE.Mesh | null,
  lowerLen: number,
): Limb {
  const root = new THREE.Group();
  root.position.set(rootPos[0], rootPos[1], rootPos[2]);
  root.add(upper);
  const mid = new THREE.Group();
  mid.position.y = -upperLen;
  mid.add(lower);
  if (tip !== null) {
    tip.position.y = -lowerLen;
    mid.add(tip);
  }
  root.add(mid);
  return { root, mid };
}

/** One posable brother. Joint conventions match figure.ts (facing +x:
 * +rotation.z on a root swings the limb forward). */
export interface BrotherFigure {
  readonly group: THREE.Group;
  readonly spec: FigureSpec;
  /** Drive the run cycle (phase in rad, ground speed m/s). */
  setGait(phaseRad: number, speedMps: number): void;
  /** Raise/lower the field glasses (overrides arm gait while up). */
  setGlasses(up: boolean): void;
  /** Aim the LEFT arm at a target in figure-local metres (x forward,
   * y up, z right), e.g. the wingtip he steadies; null releases. */
  aimLeftArm(target: readonly [number, number, number] | null): void;
}

const LEATHER = new THREE.MeshStandardMaterial({ color: 0x1c1712, roughness: 0.75, metalness: 0.1 });
const BRASS_MAT = new THREE.MeshStandardMaterial({ color: 0xbf9b40, roughness: 0.35, metalness: 0.8 });

export function createBrotherFigure(brother: Brother): BrotherFigure {
  const s = figureSpec(brother);
  const group = new THREE.Group();
  const body = new THREE.Group();
  group.add(body);
  const g = s.build;

  // Torso: tailored 1903 sack suit jacket with lapels and waist flare
  const torsoLen = s.torsoLenM;
  const torso = new THREE.Mesh(
    new THREE.CapsuleGeometry(0.165 * g, torsoLen * 0.72, 6, 12),
    SUIT,
  );
  torso.scale.set(1.3, 1, 0.82);
  torso.position.y = s.hipHeightM + torsoLen * 0.52;
  body.add(torso);

  // Jacket skirt / bottom flare
  const jacketSkirt = new THREE.Mesh(
    new THREE.CylinderGeometry(0.18 * g, 0.21 * g, 0.22, 12, 1, true),
    SUIT,
  );
  jacketSkirt.scale.set(1.25, 1, 0.85);
  jacketSkirt.position.set(0, s.hipHeightM + 0.08, 0);
  body.add(jacketSkirt);

  // Jacket lapels & buttons
  for (const side of [-1, 1] as const) {
    const lapel = new THREE.Mesh(new THREE.BoxGeometry(0.02, 0.24, 0.06), SUIT_DARK);
    lapel.position.set(0.14 * g, s.shoulderHeightM - 0.12, side * 0.07);
    lapel.rotation.set(0, side * 0.35, -0.2);
    body.add(lapel);
  }
  for (let b = 0; b < 3; b += 1) {
    const button = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.006, 6), BRASS_MAT);
    button.rotation.z = Math.PI / 2;
    button.position.set(0.145 * g, s.hipHeightM + 0.14 + b * 0.08, 0);
    body.add(button);
  }

  // Shirt collar + necktie
  const collar = new THREE.Mesh(new THREE.CylinderGeometry(0.055, 0.06, 0.055, 12), SHIRT);
  collar.position.set(0.01, s.shoulderHeightM + 0.02, 0);
  body.add(collar);
  const tie = new THREE.Mesh(new THREE.BoxGeometry(0.012, 0.18, 0.04), SUIT_DARK);
  tie.position.set(0.12 * g, s.shoulderHeightM - 0.1, 0);
  body.add(tie);

  // Head + sculpted neck, chin, nose, ears, and flat cap
  const headG = new THREE.Group();
  headG.position.y = s.shoulderHeightM + 0.065 + s.headRadiusM;
  const neck = new THREE.Mesh(new THREE.CylinderGeometry(0.048, 0.054, 0.07, 8), SKIN);
  neck.position.y = -s.headRadiusM * 0.65;
  headG.add(neck);

  const head = new THREE.Mesh(
    new THREE.SphereGeometry(s.headRadiusM, 18, 14),
    headMaterial(brother),
  );
  head.scale.set(0.92, 1.06, 0.88);
  headG.add(head);

  // Ears
  for (const side of [-1, 1] as const) {
    const ear = new THREE.Mesh(new THREE.CapsuleGeometry(0.015, 0.022, 2, 4), SKIN);
    ear.position.set(-0.01, 0, side * s.headRadiusM * 0.92);
    ear.rotation.set(0, 0, 0.2);
    headG.add(ear);
  }

  // Flat Tweed Cap
  const cap = new THREE.Mesh(
    new THREE.CylinderGeometry(s.headRadiusM * 1.02, s.headRadiusM * 1.06, 0.048, 16),
    SUIT_DARK,
  );
  cap.position.y = s.headRadiusM * 0.82;
  headG.add(cap);
  const capCrown = new THREE.Mesh(new THREE.SphereGeometry(s.headRadiusM * 1.04, 12, 8), SUIT_DARK);
  capCrown.scale.set(1.05, 0.35, 1.05);
  capCrown.position.set(0.01, s.headRadiusM * 0.88, 0);
  headG.add(capCrown);
  const brim = new THREE.Mesh(new THREE.BoxGeometry(0.11, 0.014, 0.12), SUIT_DARK);
  brim.position.set(s.headRadiusM * 0.88, s.headRadiusM * 0.62, 0);
  brim.rotation.z = -0.15;
  headG.add(brim);
  body.add(headG);

  // Legs with trouser cuffs and leather Oxford boots
  const hipY = s.hipHeightM;
  const legs: Limb[] = [];
  for (const side of [-1, 1] as const) {
    const thigh = taperedLimb(0.082 * g, 0.065 * g, s.thighLenM, SUIT);
    const shin = taperedLimb(0.064 * g, 0.052 * g, s.shinLenM, SUIT);
    // Boot: leather upper + sole + heel
    const bootG = new THREE.Mesh(new THREE.BoxGeometry(s.footLenM * 1.05, 0.075, 0.082), LEATHER);
    bootG.position.set(s.footLenM * 0.28, 0, 0);
    const limb = buildLimb(
      [0, hipY, (side * s.hipWidthM) / 3.2],
      thigh,
      s.thighLenM,
      shin,
      bootG,
      s.shinLenM,
    );
    body.add(limb.root);
    legs.push(limb);
  }

  // Arms with white cuffs and hands
  const arms: Limb[] = [];
  for (const side of [-1, 1] as const) {
    const hand = new THREE.Mesh(new THREE.BoxGeometry(0.05, 0.08, 0.04), SKIN);
    hand.position.set(0, -0.04, 0);
    const limb = buildLimb(
      [0, s.shoulderHeightM, (side * s.shoulderWidthM) / 2.1],
      taperedLimb(0.062 * g, 0.048 * g, s.upperArmLenM, SUIT),
      s.upperArmLenM,
      taperedLimb(0.046 * g, 0.038 * g, s.forearmLenM, SUIT),
      hand,
      s.forearmLenM,
    );
    body.add(limb.root);
    arms.push(limb);
  }
  const [legL, legR] = [legs[0]!, legs[1]!];
  const [armL, armR] = [arms[0]!, arms[1]!];

  // Authentic 1903 Brass Field Glasses (Binoculars)
  const glasses = new THREE.Group();
  for (const dz of [-0.035, 0.035]) {
    const barrel = new THREE.Mesh(new THREE.CylinderGeometry(0.024, 0.018, 0.13, 10), LEATHER);
    barrel.rotation.z = Math.PI / 2;
    barrel.position.set(0, 0, dz);
    glasses.add(barrel);
    const brassRim = new THREE.Mesh(new THREE.CylinderGeometry(0.026, 0.026, 0.015, 10), BRASS_MAT);
    brassRim.rotation.z = Math.PI / 2;
    brassRim.position.set(0.065, 0, dz);
    glasses.add(brassRim);
    const eyePiece = new THREE.Mesh(new THREE.CylinderGeometry(0.016, 0.016, 0.018, 10), BRASS_MAT);
    eyePiece.rotation.z = Math.PI / 2;
    eyePiece.position.set(-0.065, 0, dz);
    glasses.add(eyePiece);
  }
  const centerBridge = new THREE.Mesh(new THREE.CylinderGeometry(0.008, 0.008, 0.08, 6), BRASS_MAT);
  centerBridge.rotation.x = Math.PI / 2;
  glasses.add(centerBridge);
  const focusKnob = new THREE.Mesh(new THREE.CylinderGeometry(0.014, 0.014, 0.02, 10), BRASS_MAT);
  focusKnob.position.set(0, 0.016, 0);
  glasses.add(focusKnob);

  glasses.position.set(s.headRadiusM + 0.08, 0.02, 0);
  glasses.visible = false;
  headG.add(glasses);

  let glassesUp = false;
  let leftAim: readonly [number, number, number] | null = null;
  const armLShoulderY = s.shoulderHeightM;
  castAll(group);

  return {
    group,
    spec: s,
    setGait(phaseRad: number, speedMps: number): void {
      const p = gaitPose(phaseRad, speedMps);
      legL.root.rotation.z = p.hipL;
      legL.mid.rotation.z = -p.kneeL;
      legR.root.rotation.z = p.hipR;
      legR.mid.rotation.z = -p.kneeR;
      torso.rotation.z = p.leanRad;
      body.position.y = p.bobM;
      if (glassesUp) {
        for (const arm of [armL, armR]) {
          arm.root.rotation.set(0, 0, BINOCULAR_POSE.shoulderForwardRad);
          arm.mid.rotation.z = BINOCULAR_POSE.elbowFlexRad;
        }
      } else {
        if (leftAim !== null) {
          const [tx, ty, tz] = leftAim;
          const aim = armAimAngles(
            tx,
            ty - armLShoulderY,
            tz - (-s.shoulderWidthM / 2.15),
          );
          armL.root.rotation.order = "YZX";
          armL.root.rotation.set(0, aim.yawRad, aim.pitchRad);
          armL.mid.rotation.z = 0.1;
        } else {
          armL.root.rotation.set(0, 0, p.shoulderL);
          armL.mid.rotation.z = p.elbowL;
        }
        armR.root.rotation.set(0, 0, p.shoulderR);
        armR.mid.rotation.z = p.elbowR;
      }
      glasses.visible = glassesUp;
      headG.rotation.z = glassesUp ? 0.14 : 0; // eyes looking up at aircraft
    },
    setGlasses(up: boolean): void {
      glassesUp = up;
    },
    aimLeftArm(target: readonly [number, number, number] | null): void {
      leftAim = target;
    },
  };
}

/** Prone pilot on the cradle: lying along +x, head fore and UP,
 * hands forward on the canard lever. `setLever` couples the arms to
 * the REAL control deflection so the player sees Wilbur fly. */
export interface ProneFigure {
  readonly group: THREE.Group;
  setLever(dcRad: number): void;
  setWarp?(warpRad: number): void;
  flutter?(timeS: number, airspeedMps: number): void;
}

export function createProneBrother(brother: Brother): ProneFigure {
  const s = figureSpec(brother);
  const g = s.build;
  const group = new THREE.Group();
  // Lower body / hip pivot for cradle shifting (roll-yaw wing warping)
  const hipPivot = new THREE.Group();
  group.add(hipPivot);

  // Torso lying +x with a slight back arch; scale keeps the cradle fit.
  const torso = new THREE.Mesh(new THREE.CapsuleGeometry(0.14 * g, s.torsoLenM * 0.7, 4, 10), SUIT);
  // Capsule +y end fore-and-up = chest raised on the cradle.
  torso.rotation.z = -Math.PI / 2 + PRONE_POSE.backArchRad;
  torso.scale.set(1.25, 1, 0.7);
  torso.position.set(0.05, 0.11, 0);
  group.add(torso);

  // Fluttering coat tails on the lower back/hips
  const tailGeo = new THREE.PlaneGeometry(0.22, 0.18, 2, 2);
  const coatTail = new THREE.Mesh(tailGeo, new THREE.MeshStandardMaterial({
    color: 0x2b2a31,
    roughness: 0.95,
    side: THREE.DoubleSide,
  }));
  coatTail.rotation.set(-Math.PI / 2 + 0.2, 0, Math.PI / 2);
  coatTail.position.set(-0.24, 0.15, 0);
  hipPivot.add(coatTail);

  // Necktie fluttering under collar
  const tie = new THREE.Mesh(new THREE.BoxGeometry(0.14, 0.012, 0.035), SUIT_DARK);
  tie.position.set(0.24, 0.08, 0);
  group.add(tie);

  const headG = new THREE.Group();
  headG.position.set(s.torsoLenM * 0.5 + 0.1, 0.16, 0);
  headG.rotation.z = PRONE_POSE.headPitchRad;
  const head = new THREE.Mesh(new THREE.SphereGeometry(s.headRadiusM, 18, 14), headMaterial(brother));
  head.scale.set(0.92, 1.06, 0.88);
  headG.add(head);
  const cap = new THREE.Mesh(
    new THREE.CylinderGeometry(s.headRadiusM * 0.98, s.headRadiusM * 1.04, 0.04, 14),
    SUIT_DARK,
  );
  cap.position.y = s.headRadiusM * 0.8;
  headG.add(cap);
  const brim = new THREE.Mesh(new THREE.BoxGeometry(0.09, 0.01, 0.08), SUIT_DARK);
  brim.position.set(s.headRadiusM * 0.82, s.headRadiusM * 0.6, 0);
  headG.add(brim);
  group.add(headG);

  // Legs trail aft with a soft knee, parented to hipPivot for live warp motion.
  for (const side of [-1, 1] as const) {
    const thigh = taperedLimb(0.065 * g, 0.052 * g, s.thighLenM, SUIT);
    thigh.rotation.z = -Math.PI / 2 - PRONE_POSE.hipFlexRad;
    thigh.position.set(-0.18, 0.1, side * 0.075);
    hipPivot.add(thigh);
    const shin = taperedLimb(0.05 * g, 0.038 * g, s.shinLenM, SUIT);
    shin.rotation.z = -Math.PI / 2 - PRONE_POSE.hipFlexRad - PRONE_POSE.kneeFlexRad;
    shin.position.set(-0.18 - s.thighLenM * 0.98, 0.12, side * 0.075);
    hipPivot.add(shin);
    const boot = new THREE.Mesh(new THREE.BoxGeometry(0.16, 0.06, 0.08), SUIT_DARK);
    boot.position.set(-0.18 - s.thighLenM * 0.98 - s.shinLenM * 0.95, 0.14, side * 0.075);
    hipPivot.add(boot);
  }

  // Arms reach forward-down to the lever; hinged for setLever.
  const armRoots: THREE.Group[] = [];
  for (const side of [-1, 1] as const) {
    const hand = new THREE.Mesh(new THREE.SphereGeometry(0.04, 8, 6), SKIN);
    const limb = buildLimb(
      [s.torsoLenM * 0.42, 0.14, side * 0.19],
      taperedLimb(0.05 * g, 0.04 * g, s.upperArmLenM * 0.9, SUIT),
      s.upperArmLenM * 0.9,
      taperedLimb(0.038 * g, 0.03 * g, s.forearmLenM * 0.9, SUIT),
      hand,
      s.forearmLenM * 0.9,
    );
    // Rest pose: shoulders rolled forward, elbows soft.
    limb.root.rotation.z = PRONE_POSE.shoulderForwardRad;
    limb.mid.rotation.z = PRONE_POSE.elbowFlexRad;
    group.add(limb.root);
    armRoots.push(limb.root);
  }
  castAll(group);
  return {
    group,
    setLever(dcRad: number): void {
      // The lever pulls the hands with the REAL canard deflection —
      // bounded so the pose stays human at the mechanical stop.
      // Pull (+dcRad, nose-up) draws the hands AFT; after the mount yaw
      // (-pi/2 at the flyerScene mount site) that is rotation.z DECREASING.
      const pull = Math.max(-0.35, Math.min(0.35, dcRad * 0.6));
      for (const root of armRoots) {
        root.rotation.z = PRONE_POSE.shoulderForwardRad - pull;
      }
      headG.rotation.z = PRONE_POSE.headPitchRad + pull * 0.25;
    },
    setWarp(warpRad: number): void {
      // Shifting hips in the cradle drives wing warping.
      const shiftZ = Math.max(-0.06, Math.min(0.06, warpRad * 0.4));
      hipPivot.position.z = shiftZ;
      torso.position.z = shiftZ * 0.4;
      torso.rotation.y = -shiftZ * 1.5;
    },
    flutter(timeS: number, airspeedMps: number): void {
      const spd = Math.max(4, airspeedMps);
      const wave = Math.sin(timeS * (spd * 1.2)) * 0.12;
      coatTail.rotation.x = -Math.PI / 2 + 0.2 + wave;
      tie.rotation.y = Math.sin(timeS * 14.5) * 0.15;
    },
  };
}
