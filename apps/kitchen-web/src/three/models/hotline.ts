/**
 * The hot line: ranges, burners and grates. Dimensions are the archetypes' (assets/kitchens/specs), the layouts those
 * of the equipment class: a 36 in restaurant range is six open star burners under cast-iron pot grates, a valve panel
 * of knobs over a standard oven, a flue riser and high shelf behind, all on 6 in legs.
 */
import * as THREE from "three";
import { ringPorts, FlameSet, type Port } from "./flames";
import { register, type Builder, type Hooks } from "./index";
import { palette } from "./materials";
import { bar, cyl, knob, lathe, plate, rbox, screw, taper, turnJoint, type Shop, type V3 } from "./shop";

/* ------------------------------------------------------------------------------------------------------------ */
/* Parts                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/** The outline of a six-lobed star burner head of mean radius `r`: points and outward normals, counter-clockwise. */
function starOutline(r: number, lobes: number, depth: number, n: number): { p: [number, number][]; nrm: [number, number][] } {
  const p: [number, number][] = [];
  const nrm: [number, number][] = [];
  for (let i = 0; i < n; i++) {
    const a = (i / n) * Math.PI * 2;
    const rr = r + depth * Math.cos(lobes * a);
    const dr = -depth * lobes * Math.sin(lobes * a);
    const [x, y] = [rr * Math.cos(a), rr * Math.sin(a)];
    // Tangent d/da (r cos a, r sin a); the outward normal is it turned clockwise.
    const [tx, ty] = [dr * Math.cos(a) - rr * Math.sin(a), dr * Math.sin(a) + rr * Math.cos(a)];
    const l = Math.hypot(tx, ty);
    p.push([x, y]);
    nrm.push([ty / l, -tx / l]);
  }
  return { p, nrm };
}

/**
 * A lift-off cast-iron star burner centred at (x, z) with its port line at height y, on its venturi throat. Returns
 * the ports, numbered for `burner`. `r` is the head's mean radius (a 30 000 BTU/h head is about 55 mm).
 */
export function starBurner(s: Shop, x: number, y: number, z: number, burner: number, r = 0.052): Port[] {
  const p = palette();
  const lobes = 6;
  const depth = r * 0.2;
  const { p: outline, nrm } = starOutline(r, lobes, depth, s.seg(96, 36));
  // The head: a star-shaped ring 14 mm tall (ports drilled through its wall), a domed cap inside it.
  const shape = new THREE.Shape(outline.map(([a, b]) => new THREE.Vector2(a, b)));
  const inner = starOutline(r - 0.008, lobes, depth * 0.9, s.seg(96, 36)).p.reverse();
  shape.holes.push(new THREE.Path(inner.map(([a, b]) => new THREE.Vector2(a, b))));
  const wall = plate(s, shape, 0.014, 0.0012);
  wall.rotateX(-Math.PI / 2);
  s.add(p.castIron, wall, [x, y - 0.008, z]);
  s.add(p.castIron, lathe(s, [
    [0, 0.009],
    [r * 0.55, 0.008],
    [r * 0.82, 0.004],
    [r * 0.86, 0],
    [r * 0.86, -0.004],
    [0, -0.004],
  ], 40), [x, y, z]);
  // Throat and bowl beneath, where the venturi delivers the premixed gas.
  s.add(p.castIron, lathe(s, [
    [r * 0.45, -0.008],
    [r * 0.7, -0.012],
    [r * 0.72, -0.024],
    [0.019, -0.034],
    [0.017, -0.07],
  ], 32, 0, Math.PI * 2), [x, y, z]);
  // The ports: a hole every few millimetres round the outline, each a dark dot close up.
  const n = s.fine(1) ? 36 : 24;
  const ports: Port[] = [];
  for (let i = 0; i < n; i++) {
    const k = Math.floor((i / n) * outline.length);
    const [px, pz] = outline[k]!;
    const [nx, nz] = nrm[k]!;
    // Shape y maps to -z after the plate's turn.
    ports.push({ at: [x + px, y - 0.001, z - pz], dir: [nx * 0.85, 0.52, -nz * 0.85], burner });
    // The drilled hole: a dark dot on the wall, its axis along the outward normal (nx, 0, -nz).
    if (s.fine(2)) s.add(p.enamel, cyl(s, 0.0011, 0.003, "x", 6), [x + px, y - 0.001, z - pz], [0, Math.atan2(nz, nx), 0]);
  }
  return ports;
}

/**
 * A cast-iron pot grate, `w` by `d`, its flat top at `top` and its feet `h` below, centred at (x, z): a square frame
 * with eight fingers reaching in toward the burner, their tips the rests a pan stands on over the flame.
 */
export function potGrate(s: Shop, x: number, z: number, w: number, d: number, top: number, h = 0.034): void {
  const p = palette();
  const t = 0.012;
  const fh = 0.017;
  const yc = top - fh / 2;
  // The frame.
  for (const sz of [-1, 1]) s.add(p.castIron, rbox(s, w, fh, t, 0.0035), [x, yc, z + sz * (d / 2 - t / 2)]);
  for (const sx of [-1, 1]) s.add(p.castIron, rbox(s, t, fh, d - 2 * t + 0.002, 0.0035), [x + sx * (w / 2 - t / 2), yc, z]);
  // Feet at the corners, down to the burner box rim.
  const footH = h - fh + 0.002;
  for (const [fx, fz] of [
    [-1, -1],
    [1, -1],
    [-1, 1],
    [1, 1],
  ] as const)
    s.add(p.castIron, rbox(s, 0.02, footH, 0.02, 0.003), [x + fx * (w / 2 - 0.01), top - fh - footH / 2 + 0.001, z + fz * (d / 2 - 0.01)]);
  // Fingers: from the corners and mid-sides in toward the burner, tapering to their tips.
  const rest = Math.min(w, d) * 0.21;
  for (let i = 0; i < 8; i++) {
    const a = (i / 8) * Math.PI * 2 + Math.PI / 4;
    const [c, sn] = [Math.cos(a), Math.sin(a)];
    // Where the ray from the centre meets the frame's inner edge, overlapping it a little.
    const k = Math.min((w / 2 - t * 0.5) / Math.max(1e-6, Math.abs(c)), (d / 2 - t * 0.5) / Math.max(1e-6, Math.abs(sn)));
    const [ox, oz] = [c * k, sn * k];
    const [ix, iz] = [c * rest, sn * rest];
    const len = Math.hypot(ox - ix, oz - iz);
    // Root at the frame, pointing in: +x of the taper along (ix - ox, iz - oz).
    const yaw = -Math.atan2(iz - oz, ix - ox);
    s.add(p.castIron, taper(s, len, t * 1.15, t * 0.7, fh, 0.0032), [x + ox, top - fh, z + oz], [0, yaw, 0]);
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Ranges                                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

interface RangeLayout {
  /** Burner centres across, and their depths behind the front edge. */
  xs: number[];
  zs: number[];
}

/** Knob angle for a valve opening: off points up; turned 90° counter-clockwise is high, on round to low at 180°. */
export function valveAngle(flow: number): number {
  return flow <= 0.01 ? 0 : Math.PI / 2 + (1 - flow) * (Math.PI / 2);
}

/** A restaurant range: open burners over a standard oven. The 36 in archetype has six burners. */
function openRange(L: RangeLayout): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const legH = 0.152;
    const zf = -0.056; // the body's front face; knobs and the door handle stand forward of it
    const deck = 0.88; // the top of the burner box, under the grates
    const top = 0.914; // the grates' top: the cooking height
    const back = -d;
    const riserD = 0.095;
    const bodyBack = back + 0.012;

    // Legs: 6 in stainless tube legs with adjustable flanged feet.
    for (const lx of [0.055, w - 0.055])
      for (const lz of [zf - 0.06, bodyBack + 0.06]) {
        s.add(p.ss, cyl(s, 0.019, legH - 0.03, "y", 20), [lx, 0.03 + (legH - 0.03) / 2, lz]);
        s.add(p.ssBright, lathe(s, [
          [0, 0],
          [0.028, 0],
          [0.03, 0.006],
          [0.022, 0.012],
          [0.019, 0.03],
          [0, 0.03],
        ], 24), [lx, 0, lz]);
      }

    // The body: sides, back, bottom. The sides are single sheets with a folded front flange.
    const bodyH = deck - legH;
    for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, rbox(s, 0.0016, bodyH, zf - bodyBack, 0.0006), [sx, legH + bodyH / 2, (zf + bodyBack) / 2]);
    s.add(p.ssDull, rbox(s, w - 0.004, bodyH, 0.0016, 0.0006), [w / 2, legH + bodyH / 2, bodyBack + 0.0008]);
    s.add(p.ssDull, rbox(s, w - 0.004, 0.0016, zf - bodyBack, 0.0006), [w / 2, legH + 0.0008, (zf + bodyBack) / 2]);
    for (const sx of [0.012, w - 0.012]) s.add(p.ss, rbox(s, 0.024, bodyH - 0.003, 0.0016, 0.0006), [sx, legH + bodyH / 2, zf]);

    // Kick plate, set 8 mm back.
    s.add(p.ss, rbox(s, w - 0.05, 0.06, 0.0016, 0.0006), [w / 2, legH + 0.03, zf - 0.008]);
    if (s.fine(2)) for (const sx of [0.04, w - 0.04]) screw(s, p.ssBright, [sx, legH + 0.03, zf - 0.0072], 0.0028);

    // The oven door on its bottom hinge: a double-skinned panel, its bar handle on two standoffs.
    const doorY0 = legH + 0.064;
    const doorY1 = 0.735;
    const doorW = w - 0.05;
    const door = s.joint("oven_door", [w / 2, doorY0, zf]);
    const dh = doorY1 - doorY0;
    door.add(p.ss, rbox(door, doorW, dh, 0.03, 0.004), [0, dh / 2, 0.015 - 0.008]);
    // A pressed panel in the door face and the brand-free badge plate.
    if (s.fine(1)) door.add(p.ss, rbox(door, doorW - 0.07, dh - 0.11, 0.002, 0.004), [0, dh / 2 - 0.02, 0.0075]);
    door.add(p.ssBright, bar(door, [[-doorW / 2 + 0.07, dh - 0.045, 0.006], [-doorW / 2 + 0.07, dh - 0.045, 0.044], [doorW / 2 - 0.07, dh - 0.045, 0.044], [doorW / 2 - 0.07, dh - 0.045, 0.006]], 0.0095, 0.02, 14));
    for (const sx of [-1, 1]) door.add(p.ssBright, cyl(door, 0.015, 0.006, "z", 20), [sx * (doorW / 2 - 0.07), dh - 0.045, 0.009]);
    // The door's inner liner and its gasket, seen when it opens.
    door.add(p.enamel, rbox(door, doorW - 0.06, dh - 0.06, 0.004, 0.003), [0, dh / 2, -0.009]);
    // Hinge knuckles at the bottom corners.
    for (const sx of [-1, 1]) door.add(p.ssBright, cyl(door, 0.007, 0.03, "x", 14), [sx * (doorW / 2 - 0.04), 0.004, -0.002]);

    // The oven cavity behind the door: porcelain-enamelled, two wire racks, a lamp.
    if (s.fine(1)) {
      const cw = doorW - 0.08;
      const cd = zf - bodyBack - 0.1;
      const cy0 = doorY0 + 0.03;
      const ch = dh - 0.06;
      const cz = zf - 0.02 - cd / 2;
      s.add(p.enamel, rbox(s, cw, 0.003, cd, 0.001), [w / 2, cy0, cz]);
      s.add(p.enamel, rbox(s, cw, 0.003, cd, 0.001), [w / 2, cy0 + ch, cz]);
      for (const sx of [-1, 1]) s.add(p.enamel, rbox(s, 0.003, ch, cd, 0.001), [w / 2 + (sx * cw) / 2, cy0 + ch / 2, cz]);
      s.add(p.enamel, rbox(s, cw, ch, 0.003, 0.001), [w / 2, cy0 + ch / 2, cz - cd / 2]);
      for (const ry of [0.33, 0.62]) {
        const y = cy0 + ch * ry;
        const rods = s.fine(2) ? 14 : 7;
        for (let i = 0; i < rods; i++) s.add(p.chrome, cyl(s, 0.0022, cd - 0.02, "z", 8), [w / 2 - cw / 2 + 0.02 + ((cw - 0.04) * i) / (rods - 1), y, cz]);
        for (const sz of [-1, 1]) s.add(p.chrome, cyl(s, 0.003, cw - 0.01, "x", 8), [w / 2, y + 0.002, cz + sz * (cd / 2 - 0.015)]);
      }
      s.add(p.ledAmber, cyl(s, 0.012, 0.004, "z", 16), [w / 2 + cw / 2 - 0.05, cy0 + ch - 0.04, cz - cd / 2 + 0.003]);
    }

    // The valve panel: stainless, a little proud of the door line, its knobs on chrome bezels.
    const panelY0 = 0.742;
    const panelY1 = 0.86;
    s.add(p.ss, rbox(s, w, panelY1 - panelY0, 0.02, 0.003), [w / 2, (panelY0 + panelY1) / 2, zf - 0.006]);
    // The bull-nose rail above it.
    s.add(p.ss, cyl(s, 0.012, w, "x", 24), [w / 2, deck - 0.006, zf - 0.006]);
    const valves = L.xs.length * L.zs.length;
    const knobs = valves + 1;
    const knobY = (panelY0 + panelY1) / 2 - 0.004;
    const knobJoints: string[] = [];
    for (let i = 0; i < knobs; i++) {
      const kx = 0.07 + ((w - 0.14) * i) / (knobs - 1);
      const name = `knob_${i}`;
      const j = s.joint(name, [kx, knobY, zf + 0.004]);
      knob(j, p.bakelite, p.ssBright, [0, 0, 0], 0.021);
      knobJoints.push(name);
      // A printed dial ring round each: a thin dark annulus.
      if (s.fine(2)) s.add(p.enamel, cyl(s, 0.031, 0.0006, "z", 32), [kx, knobY, zf + 0.0045]);
    }
    if (s.fine(2)) for (const sx of [0.025, w - 0.025]) for (const sy of [panelY0 + 0.015, panelY1 - 0.015]) screw(s, p.ssBright, [sx, sy, zf + 0.004], 0.0028);

    // The cooktop: a burner box with its stainless rim, the drip tray under the burners, and the grates.
    s.add(p.ss, rbox(s, w, 0.004, 0.03, 0.0015), [w / 2, deck, zf - 0.015]);
    for (const sx of [0.012, w - 0.012]) s.add(p.ss, rbox(s, 0.024, 0.004, zf - bodyBack - 0.03, 0.0015), [sx, deck, (zf - 0.03 + bodyBack) / 2]);
    s.add(p.ssHeat, rbox(s, w - 0.04, 0.003, zf - bodyBack - 0.05, 0.001), [w / 2, deck - 0.075, (zf - 0.03 + bodyBack) / 2]);
    s.add(p.ss, rbox(s, w - 0.03, 0.05, 0.0016, 0.001), [w / 2, deck - 0.05, zf - 0.03]);
    // Grates abut, one over each burner: as wide and deep as the burners' pitch allows, 4 mm apart.
    const pitch = (v: number[]) => (v.length > 1 ? Math.min(...v.slice(1).map((x, i) => x - v[i]!)) : 0.3);
    const gw = pitch(L.xs) - 0.004;
    const gd = pitch(L.zs) - 0.004;
    const ports: Port[] = [];
    let burner = 0;
    for (const bz of L.zs)
      for (const bx of L.xs) {
        ports.push(...starBurner(s, bx, deck + 0.008, -bz, burner++));
        potGrate(s, bx, -bz, gw, gd, top);
      }
    // Standing pilots, one between each front and back pair, a pinpoint flame each; the flame set numbers them last.
    for (const bx of L.xs) {
      const pz = -(L.zs[0]! + L.zs[L.zs.length - 1]!) / 2;
      s.add(p.ssBright, cyl(s, 0.0035, 0.03, "y", 10), [bx, deck - 0.01, pz]);
      s.add(p.ssBright, bar(s, [[bx, deck - 0.025, pz], [bx, deck - 0.04, pz + 0.03], [bx, deck - 0.06, zf - 0.04]], 0.0022, 0.01, 8));
      ports.push(...ringPorts(valves, [bx, deck + 0.006, pz], 0.0035, 3, 1.25));
    }

    // Flue riser and high shelf: the oven's flue vents through a slotted grille at the back of the cooktop.
    const riserFront = back + riserD;
    s.add(p.ss, rbox(s, w, h - deck, 0.0016, 0.0006), [w / 2, (deck + h) / 2, riserFront]);
    s.add(p.ssDull, rbox(s, w, h - deck, 0.0016, 0.0006), [w / 2, (deck + h) / 2, back + 0.001]);
    for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, rbox(s, 0.0016, h - deck, riserD, 0.0006), [sx, (deck + h) / 2, back + riserD / 2]);
    // The oven flue: a dark slot across the riser's foot, behind three angled louvre slats.
    s.add(p.enamel, rbox(s, w - 0.1, 0.064, 0.004, 0.002), [w / 2, deck + 0.05, riserFront + 0.0015]);
    for (let i = 0; i < 3; i++) s.add(p.ss, rbox(s, w - 0.1, 0.014, 0.0016, 0.0006), [w / 2, deck + 0.03 + i * 0.02, riserFront + 0.008], [-0.6, 0, 0]);
    for (const sx of [0.05, w - 0.05]) s.add(p.ss, rbox(s, 0.006, 0.07, 0.014, 0.001), [sx, deck + 0.05, riserFront + 0.006]);
    // The high shelf: a folded sheet, its front edge turned down.
    const shelfY = h - 0.018;
    const shelfD = 0.27;
    s.add(p.ss, rbox(s, w, 0.0016, shelfD, 0.0006), [w / 2, shelfY, back + shelfD / 2]);
    s.add(p.ss, rbox(s, w, 0.036, 0.0016, 0.0008), [w / 2, shelfY - 0.017, back + shelfD]);
    s.add(p.ss, rbox(s, w, 0.012, 0.0016, 0.0006), [w / 2, shelfY + 0.006, back + 0.0008]);
    // Shelf brackets.
    for (const sx of [0.03, w - 0.03]) {
      const tri = new THREE.Shape([new THREE.Vector2(0, 0), new THREE.Vector2(shelfD - riserD - 0.01, 0), new THREE.Vector2(0, -0.12)]);
      const g = plate(s, tri, 0.0016);
      g.rotateY(-Math.PI / 2);
      s.add(p.ss, g, [sx, shelfY - 0.001, riserFront]);
    }

    // The gas supply behind: a 3/4 in black-iron stub with its brass shut-off.
    if (s.fine(1)) {
      s.add(p.enamel, cyl(s, 0.0135, 0.12, "z", 16), [w * 0.8, 0.5, back - 0.03]);
      s.add(p.brass, rbox(s, 0.05, 0.035, 0.035, 0.004), [w * 0.8, 0.5, back - 0.09]);
      s.add(p.enamel, rbox(s, 0.09, 0.012, 0.012, 0.003), [w * 0.8 + 0.045, 0.522, back - 0.09]);
    }

    return (b) => {
      const flames = new FlameSet(ports, 0.024);
      b.group.add(flames.mesh);
      flames.setAll(0);
      flames.setFlow(valves, 0.12);
      const hooks: Hooks = {
        flames,
        burners: valves,
        open: (f) => turnJoint(b.joints.get("oven_door"), "x", f * (Math.PI / 2) * 0.98),
        tick: (t) => flames.tick(t),
      };
      // Knobs follow their valves.
      const set = flames.setFlow.bind(flames);
      flames.setFlow = (i: number, f: number) => {
        set(i, f);
        if (i < valves) turnJoint(b.joints.get(knobJoints[i]!), "z", valveAngle(f));
      };
      return hooks;
    };
  };
}

register({
  "NSF-RANGE-GAS-6B": openRange({ xs: [0.18, 0.457, 0.734], zs: [0.22, 0.58] }),
});

export type { V3 };
