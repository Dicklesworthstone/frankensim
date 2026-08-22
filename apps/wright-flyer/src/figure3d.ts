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

export function createBrotherFigure(brother: Brother): BrotherFigure {
  const s = figureSpec(brother);
  const group = new THREE.Group();
  // The caller owns group.position (terrain seat); the gait bob rides
  // this inner pivot so the two can never fight.
  const body = new THREE.Group();
  group.add(body);
  const g = s.build;
  // Torso: capsule hip->shoulder, jacket-colored, squashed front-back.
  const torsoLen = s.torsoLenM;
  const torso = new THREE.Mesh(
    new THREE.CapsuleGeometry(0.155 * g, torsoLen * 0.72, 4, 10),
    SUIT,
  );
  torso.scale.set(1.35, 1, 0.78);
  torso.position.y = s.hipHeightM + torsoLen * 0.52;
  body.add(torso);
  // Shirt collar + tie.
  const collar = new THREE.Mesh(new THREE.CylinderGeometry(0.052, 0.056, 0.05, 12), SHIRT);
  collar.position.set(0.01, s.shoulderHeightM + 0.02, 0);
  body.add(collar);
  const tie = new THREE.Mesh(new THREE.BoxGeometry(0.012, 0.16, 0.035), SUIT_DARK);
  tie.position.set(0.115 * g, s.shoulderHeightM - 0.1, 0);
  body.add(tie);
  // Head + flat cap (u=0.5 of the sphere faces +x, where the photo is).
  const headG = new THREE.Group();
  headG.position.y = s.shoulderHeightM + 0.065 + s.headRadiusM;
  const head = new THREE.Mesh(
    new THREE.SphereGeometry(s.headRadiusM, 18, 14),
    headMaterial(brother),
  );
  head.scale.set(0.92, 1.06, 0.88);
  headG.add(head);
  const cap = new THREE.Mesh(
    new THREE.CylinderGeometry(s.headRadiusM * 0.98, s.headRadiusM * 1.04, 0.045, 14),
    SUIT_DARK,
  );
  cap.position.y = s.headRadiusM * 0.82;
  headG.add(cap);
  const brim = new THREE.Mesh(new THREE.BoxGeometry(0.1, 0.012, 0.09), SUIT_DARK);
  brim.position.set(s.headRadiusM * 0.85, s.headRadiusM * 0.62, 0);
  headG.add(brim);
  body.add(headG);
  // Legs (trouser taper), feet.
  const hipY = s.hipHeightM;
  const legs: Limb[] = [];
  for (const side of [-1, 1] as const) {
    const limb = buildLimb(
      [0, hipY, (side * s.hipWidthM) / 3.4],
      taperedLimb(0.072 * g, 0.058 * g, s.thighLenM, SUIT),
      s.thighLenM,
      taperedLimb(0.056 * g, 0.04 * g, s.shinLenM, SUIT),
      new THREE.Mesh(new THREE.BoxGeometry(s.footLenM, 0.05, 0.075), SUIT_DARK),
      s.shinLenM,
    );
    const foot = limb.mid.children[1] as THREE.Mesh;
    foot.position.x = s.footLenM * 0.28; // toes forward
    body.add(limb.root);
    legs.push(limb);
  }
  // Arms with hands.
  const arms: Limb[] = [];
  for (const side of [-1, 1] as const) {
    const hand = new THREE.Mesh(new THREE.SphereGeometry(0.042, 8, 6), SKIN);
    const limb = buildLimb(
      [0, s.shoulderHeightM, (side * s.shoulderWidthM) / 2.15],
      taperedLimb(0.055 * g, 0.042 * g, s.upperArmLenM, SUIT),
      s.upperArmLenM,
      taperedLimb(0.04 * g, 0.032 * g, s.forearmLenM, SUIT),
      hand,
      s.forearmLenM,
    );
    body.add(limb.root);
    arms.push(limb);
  }
  const [legL, legR] = [legs[0]!, legs[1]!];
  const [armL, armR] = [arms[0]!, arms[1]!];
  // Field glasses (hidden until raised).
  const glasses = new THREE.Mesh(new THREE.BoxGeometry(0.13, 0.045, 0.065), IRON);
  glasses.position.set(s.headRadiusM + 0.05, 0.01, 0);
  glasses.visible = false;
  headG.add(glasses);
  let glassesUp = false;
  let leftAim: readonly [number, number, number] | null = null;
  const armLShoulderY = s.shoulderHeightM;
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
      // Arms: glasses pose wins, then a live aim, then the gait swing.
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
      headG.rotation.z = glassesUp ? -0.12 : 0;
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
}

export function createProneBrother(brother: Brother): ProneFigure {
  const s = figureSpec(brother);
  const g = s.build;
  const group = new THREE.Group();
  // Torso lying +x with a slight back arch; scale keeps the cradle fit.
  const torso = new THREE.Mesh(new THREE.CapsuleGeometry(0.14 * g, s.torsoLenM * 0.7, 4, 10), SUIT);
  torso.rotation.z = Math.PI / 2 - PRONE_POSE.backArchRad;
  torso.scale.set(1.25, 1, 0.7);
  torso.position.set(0.05, 0.11, 0);
  group.add(torso);
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
  group.add(headG);
  // Legs trail aft with a soft knee.
  for (const side of [-1, 1] as const) {
    const thigh = taperedLimb(0.065 * g, 0.052 * g, s.thighLenM, SUIT);
    thigh.rotation.z = -Math.PI / 2 - PRONE_POSE.hipFlexRad;
    thigh.position.set(-0.18, 0.1, side * 0.075);
    group.add(thigh);
    const shin = taperedLimb(0.05 * g, 0.038 * g, s.shinLenM, SUIT);
    shin.rotation.z = -Math.PI / 2 - PRONE_POSE.hipFlexRad - PRONE_POSE.kneeFlexRad;
    shin.position.set(-0.18 - s.thighLenM * 0.98, 0.12, side * 0.075);
    group.add(shin);
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
  return {
    group,
    setLever(dcRad: number): void {
      // The lever pulls the hands with the REAL canard deflection —
      // bounded so the pose stays human at the mechanical stop.
      const pull = Math.max(-0.35, Math.min(0.35, dcRad * 0.6));
      for (const root of armRoots) {
        root.rotation.z = PRONE_POSE.shoulderForwardRad + pull;
      }
      headG.rotation.z = PRONE_POSE.headPitchRad - pull * 0.25;
    },
  };
}
