/**
 * Refrigeration and storage: walk-in coolers and freezers, reach-ins, an undercounter worktop, a sandwich makeline,
 * an ice cuber on its bin, a blast chiller, wire shelving and an aluminium speed rack. Dimensions are the archetypes'
 * (assets/kitchens/specs); the construction that of the equipment class: 4 in foamed-in-place panels cam-locked
 * together, cam-lift doors that rise as they open and fall shut, finned coils, axial fans behind wire guards,
 * magnetic gaskets, louvred condensing units with the coil and compressor showing through.
 *
 * Joints:
 *   walk-ins: cam_lift_door (rises 12 mm on its cam as it swings out), closer_body and closer_rod (the door closer
 *     follows the door), pvc_strip_<i>, evap_fan_<u>_<i>, cond_fan; and the removable parts part_roof,
 *     part_wall_right, part_front_right, part_frost_<u> (not moving: shown or hidden by set());
 *   reach-ins: front_door, or door_left and door_right; evap_fan_<i> (one per section), cond_fan;
 *   undercounter: door_left, door_right, evap_fan_0, cond_fan; makeline: flip_lid (about x, at the back of the
 *     rail), door_left, door_right, evap_fan_0, cond_fan; blast chiller: front_door, evap_fan, cond_fan;
 *   ice machine: bin_flip_door (about x, along the top of the sloped opening), cond_fan, and part_ice_0..2 (the
 *     bin's three heaps, not moving: one shown by set());
 *   wire shelving: stock (shown or hidden); speed rack: caster_swivel_0..3, pan_<i>.
 * set() names:
 *   walk-ins: door 0..1 (as open), frost 0..1 (white frost on the coils), cutaway 0|1 (hide the right wall, the front
 *     wall right of the door frame, and the roof with its condensing unit), roof 0|1 (hide only the roof), temp °C
 *     (the door-side thermometer's reading; default 3 for a cooler, -20 for a freezer);
 *   reach-ins, undercounter, makeline, blast chiller: door 0..1 (every door, as open), or one door by its joint name;
 *     temp °C (the thermometer's live reading; default 3 for a cooler, 2 for the makeline, -18 for a freezer and
 *     the chiller); makeline also lid (or flip_lid) 0..1, and its open() lifts the lid as it opens the doors;
 *   ice machine: door (or bin_flip_door) 0..1, ice 0..1 (how full the bin is: under 0.15 empty; default 0.7);
 *   wire shelving: stock 0..1 (cases, tubs and cans on the shelves; 0 hides them);
 *   speed rack: pans 0..20 (tiers loaded from the bottom, each a sheet pan carrying a hotel pan of food or a lidded
 *     one), swivel 0..2π (all four casters' heading).
 * tick(t, dt, load) spins fans with load; in a walk-in it also sways the strip curtain while the door is open.
 */
import * as THREE from "three";
import { register, type Builder, type Hooks } from "./index";
import { palette } from "./materials";
import { bar, cyl, lathe, plate, rbox, rectPath, roundLeg, roundRect, slideJoint, torus, turnJoint, type Built, type Shop, type V3 } from "./shop";

/* ------------------------------------------------------------------------------------------------------------ */
/* Atlases. Food, packaging and printed plates come from canvases drawn once, each one material for every model,  */
/* the geometry pointing its UVs at a tile.                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/**
 * Food and packaging tiles, 8 across and 6 down. The last two rows are the bake shop's: crusts and doughs mapped
 * onto each roll from above, and whole sheet pans seen from across the room.
 */
export const T = {
  tomato: 0, lettuce: 1, onion: 2, cheese: 3, pepper: 4, olive: 5, jalapeno: 6, redOnion: 7,
  cucumber: 8, ham: 9, chicken: 10, mushroom: 11, corn: 12, pickle: 13, banana: 14, spinach: 15,
  chili: 16, rice: 17, braise: 18, mac: 19, roastVeg: 20, mash: 21, sauce: 22, soup: 23,
  kraft: 24, kraftLabel: 25, canRed: 26, canGreen: 27, lidRed: 28, lidBlue: 29, lidGreen: 30, white: 31,
  crust: 32, cookie: 33, skin: 34, dough: 35, frosting: 36, parchment: 37, flourDust: 38, glaze: 39,
  panRolls: 40, panCookies: 41, panChicken: 42, panDough: 43, panCake: 44, tubWhite: 45, canLid: 46, spare: 47,
} as const;
const AT_COLS = 8;
const AT_ROWS = 6;
const TILE = 128;

type Piece = { c: string[]; n: number; w: [number, number]; h: [number, number]; shape: "rect" | "ellipse" | "ring" | "strip" | "disc" };

function drawTile(x: CanvasRenderingContext2D, t: number, base: string, pieces: Piece[], seed: number): void {
  const r = rng(seed);
  const ox = (t % AT_COLS) * TILE;
  const oy = Math.floor(t / AT_COLS) * TILE;
  x.save();
  x.beginPath();
  x.rect(ox, oy, TILE, TILE);
  x.clip();
  x.fillStyle = base;
  x.fillRect(ox, oy, TILE, TILE);
  for (const p of pieces)
    for (let i = 0; i < p.n; i++) {
      const cx = ox + r() * TILE;
      const cy = oy + r() * TILE;
      const w = p.w[0] + r() * (p.w[1] - p.w[0]);
      const h = p.h[0] + r() * (p.h[1] - p.h[0]);
      x.fillStyle = p.c[Math.floor(r() * p.c.length)]!;
      x.save();
      x.translate(cx, cy);
      x.rotate(r() * Math.PI * 2);
      x.beginPath();
      if (p.shape === "rect" || p.shape === "strip") x.rect(-w / 2, -h / 2, w, h);
      else x.ellipse(0, 0, w / 2, h / 2, 0, 0, Math.PI * 2);
      // Each piece casts a soft shadow on those under it, so the heap reads as pieces lying on pieces.
      x.shadowColor = "rgba(0,0,0,0.42)";
      x.shadowBlur = 2.5;
      x.shadowOffsetX = 0.8;
      x.shadowOffsetY = 1.3;
      x.fill();
      x.shadowColor = "rgba(0,0,0,0)";
      if (p.shape === "ring") {
        x.fillStyle = "rgba(0,0,0,0.55)";
        x.beginPath();
        x.ellipse(0, 0, w * 0.22, h * 0.22, 0, 0, Math.PI * 2);
        x.fill();
      }
      if (p.shape === "disc" || p.shape === "ring") {
        // A shading rim so each piece reads as a solid lying on the heap.
        x.strokeStyle = "rgba(0,0,0,0.22)";
        x.lineWidth = 1;
        x.beginPath();
        x.ellipse(0, 0, w / 2, h / 2, 0, 0, Math.PI * 2);
        x.stroke();
      }
      x.restore();
    }
  x.restore();
}

/** The tile's top-left corner on the atlas canvas. */
function tileAt(t: number): [number, number] {
  return [(t % AT_COLS) * TILE, Math.floor(t / AT_COLS) * TILE];
}

/** A soft radial blob: `inner` at the centre fading to `outer` at radius `r` (canvas px). */
function blob(x: CanvasRenderingContext2D, cx: number, cy: number, rx: number, ry: number, stops: [number, string][]): void {
  x.save();
  x.translate(cx, cy);
  x.scale(1, ry / rx);
  const g = x.createRadialGradient(0, 0, 0, 0, 0, rx);
  for (const [o, c] of stops) g.addColorStop(o, c);
  x.fillStyle = g;
  x.beginPath();
  x.arc(0, 0, rx, 0, Math.PI * 2);
  x.fill();
  x.restore();
}

/**
 * The bake shop's tiles. The single-piece tiles (crust, cookie, skin, dough) are seen from straight above, centred,
 * so a dome mapped from above takes the colour of a roll: deep brown on the crown paling to the cream of its
 * sides. The pan tiles are whole 18 × 26 in sheet pans (u across, v along), for the rack seen from the aisle.
 */
function drawBakery(x: CanvasRenderingContext2D): void {
  const r = rng(77);
  const speck = (ox: number, oy: number, n: number, cols: string[], s0: number, s1: number, rad = TILE / 2) => {
    for (let i = 0; i < n; i++) {
      const a = r() * Math.PI * 2;
      const d = Math.sqrt(r()) * rad;
      x.fillStyle = cols[Math.floor(r() * cols.length)]!;
      const sz = s0 + r() * (s1 - s0);
      x.beginPath();
      x.ellipse(ox + TILE / 2 + Math.cos(a) * d, oy + TILE / 2 + Math.sin(a) * d, sz, sz * (0.5 + r() * 0.5), r() * 3, 0, Math.PI * 2);
      x.fill();
    }
  };
  const fill = (t: number, col: string) => {
    const [ox, oy] = tileAt(t);
    x.fillStyle = col;
    x.fillRect(ox, oy, TILE, TILE);
    return [ox, oy] as const;
  };
  {
    // A baked roll's crust: glossy mahogany crown, golden shoulders, pale where it rose against its neighbours.
    const [ox, oy] = fill(T.crust, "#e9cf9a");
    blob(x, ox + 64, oy + 64, 64, 64, [[0, "#7c3c12"], [0.3, "#94521c"], [0.55, "#b4702c"], [0.75, "#cf9547"], [0.9, "#e0b36e"], [1, "#ebcc96"]]);
    speck(ox, oy, 60, ["rgba(255,230,180,0.18)", "rgba(60,25,5,0.16)"], 1, 3, 40);
  }
  {
    // A chocolate-chip cookie: golden, darker at the rim, cracked, studded with chips.
    const [ox, oy] = fill(T.cookie, "#b87a35");
    blob(x, ox + 64, oy + 64, 64, 64, [[0, "#e0b66e"], [0.55, "#d5a256"], [0.85, "#b97c36"], [1, "#8f5822"]]);
    x.strokeStyle = "rgba(120,70,25,0.45)";
    x.lineWidth = 1.2;
    for (let i = 0; i < 9; i++) {
      x.beginPath();
      let px = ox + 64 + (r() - 0.5) * 60;
      let py = oy + 64 + (r() - 0.5) * 60;
      x.moveTo(px, py);
      for (let k = 0; k < 3; k++) {
        px += (r() - 0.5) * 22;
        py += (r() - 0.5) * 22;
        x.lineTo(px, py);
      }
      x.stroke();
    }
    speck(ox, oy, 16, ["#3a1e10", "#2a140a", "#4a2814"], 3.5, 6, 46);
  }
  {
    // Roast chicken skin: mottled amber and brown, crisp patches, a few herb flecks.
    const [ox, oy] = fill(T.skin, "#b8732e");
    for (let i = 0; i < 70; i++) blob(x, ox + r() * TILE, oy + r() * TILE, 8 + r() * 18, 6 + r() * 14, [[0, r() < 0.5 ? "rgba(122,58,18,0.55)" : "rgba(230,170,90,0.5)"], [1, "rgba(0,0,0,0)"]]);
    speck(ox, oy, 40, ["#2a1a0c", "#3a2412", "#1c120a"], 0.7, 1.4, 64);
  }
  {
    // Proofing dough: pale and matt, a dusting of flour.
    const [ox, oy] = fill(T.dough, "#e6d2aa");
    blob(x, ox + 64, oy + 64, 64, 64, [[0, "#efdfbf"], [0.7, "#e8d4ae"], [1, "#dcc49a"]]);
    speck(ox, oy, 90, ["rgba(255,255,255,0.5)"], 0.6, 1.4, 60);
  }
  {
    // Chocolate frosting swirled on with a palette knife: overlapping arcs, glossy ridges and soft hollows.
    const [ox, oy] = fill(T.frosting, "#4a2716");
    x.save();
    x.beginPath();
    x.rect(ox, oy, TILE, TILE);
    x.clip();
    for (let i = 0; i < 46; i++) {
      const cx = ox + r() * TILE;
      const cy = oy + r() * TILE;
      const rad = 8 + r() * 18;
      const a0 = r() * Math.PI * 2;
      x.lineWidth = 2 + r() * 3;
      x.strokeStyle = r() < 0.5 ? "rgba(122,70,40,0.5)" : "rgba(34,15,6,0.45)";
      x.beginPath();
      x.arc(cx, cy, rad, a0, a0 + 1.2 + r() * 2);
      x.stroke();
    }
    x.restore();
  }
  {
    // Parchment: off-white, faint creases.
    const [ox, oy] = fill(T.parchment, "#ece4d0");
    x.strokeStyle = "rgba(150,130,100,0.18)";
    for (let i = 0; i < 7; i++) {
      x.beginPath();
      x.moveTo(ox + r() * TILE, oy);
      x.lineTo(ox + r() * TILE, oy + TILE);
      x.stroke();
    }
  }
  {
    const [ox, oy] = fill(T.flourDust, "#f3ead8");
    speck(ox, oy, 160, ["rgba(255,255,255,0.7)", "rgba(220,200,160,0.3)"], 0.6, 1.6, 64);
  }
  {
    const [ox, oy] = fill(T.glaze, "#c88a3a");
    blob(x, ox + 64, oy + 64, 64, 64, [[0, "#f1d49a"], [0.7, "#d9a45a"], [1, "#b8782e"]]);
  }
  // Whole pans for the far view: 4 across, 6 along, on parchment.
  const pan = (t: number, piece: (cx: number, cy: number, rx: number, ry: number) => void, nu = 4, nv = 6) => {
    const [ox, oy] = fill(t, "#e8dfc9");
    x.fillStyle = "#9da3a8";
    x.fillRect(ox, oy, TILE, 3);
    x.fillRect(ox, oy + TILE - 3, TILE, 3);
    x.fillRect(ox, oy, 3, TILE);
    x.fillRect(ox + TILE - 3, oy, 3, TILE);
    for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) piece(ox + ((i + 0.5) / nu) * TILE, oy + ((j + 0.5) / nv) * TILE, (TILE / nu) * 0.42, (TILE / nv) * 0.42);
  };
  pan(T.panRolls, (cx, cy, rx, ry) => blob(x, cx, cy, rx, ry, [[0, "#7a3f16"], [0.5, "#a8662a"], [0.85, "#d9a75e"], [1, "rgba(230,200,150,0.6)"]]));
  pan(T.panCookies, (cx, cy, rx, ry) => {
    blob(x, cx, cy, rx, ry, [[0, "#e0b66e"], [0.8, "#c58c44"], [1, "#94602a"]]);
    x.fillStyle = "#3a1e10";
    for (let k = 0; k < 3; k++) x.fillRect(cx + (r() - 0.5) * rx, cy + (r() - 0.5) * ry, 2, 2);
  });
  pan(T.panChicken, (cx, cy, rx, ry) => blob(x, cx + (r() - 0.5) * 3, cy, rx * 1.05, ry * 0.95, [[0, "#9a5220"], [0.6, "#b8732e"], [1, "#d9a060"]]), 3, 5);
  pan(T.panDough, (cx, cy, rx, ry) => blob(x, cx, cy, rx, ry, [[0, "#f7efdc"], [0.8, "#ecdcbc"], [1, "#d8c49c"]]));
  {
    const [ox, oy] = tileAt(T.panCake);
    x.drawImage(x.canvas, ...tileAt(T.frosting), TILE, TILE, ox, oy, TILE, TILE);
  }
  {
    // A white HDPE tub's wall: faint mould lines and the graduations moulded in.
    const [ox, oy] = fill(T.tubWhite, "#eeede7");
    x.fillStyle = "rgba(0,0,0,0.06)";
    for (let i = 0; i < 8; i++) x.fillRect(ox + 10, oy + 16 + i * 12, 14, 1.5);
    x.fillStyle = "rgba(0,0,0,0.035)";
    x.fillRect(ox + 60, oy, 1, TILE);
  }
  {
    // A can's lid: tinplate rings pressed in.
    const [ox, oy] = fill(T.canLid, "#b9bcbf");
    for (let i = 0; i < 6; i++) {
      x.strokeStyle = i % 2 ? "rgba(255,255,255,0.35)" : "rgba(0,0,0,0.22)";
      x.lineWidth = 2;
      x.beginPath();
      x.arc(ox + 64, oy + 64, 60 - i * 7, 0, Math.PI * 2);
      x.stroke();
    }
  }
  fill(T.spare, "#808080");
}

let goodsCache: THREE.MeshStandardMaterial | null = null;
/** The food and packaging atlas: diced and shredded ingredients, cooked pans, kraft cases, can labels, lid colours. */
export function goodsMaterial(): THREE.MeshStandardMaterial {
  if (goodsCache) return goodsCache;
  const c = document.createElement("canvas");
  c.width = AT_COLS * TILE;
  c.height = AT_ROWS * TILE;
  const x = c.getContext("2d")!;
  const dice = (c0: string[], n = 260, s: [number, number] = [5, 9]): Piece => ({ c: c0, n, w: s, h: s, shape: "rect" });
  const shred = (c0: string[], n = 420, l: [number, number] = [10, 26], t: [number, number] = [1.5, 3]): Piece => ({ c: c0, n, w: l, h: t, shape: "strip" });
  drawTile(x, T.tomato, "#7a130c", [dice(["#c8261a", "#b01d12", "#e0402a", "#9c160d"], 300, [6, 10]), dice(["#f3c3a0", "#e8b18a"], 40, [2, 3])], 11);
  drawTile(x, T.lettuce, "#3f6a1f", [shred(["#8fc54a", "#a9d66a", "#6ea534", "#c8e88f", "#5d8f2a"], 520, [12, 30], [2, 4])], 12);
  drawTile(x, T.onion, "#cfc8b6", [dice(["#f3efe2", "#e9e3cf", "#fbf8ef", "#ddd5bd"], 320, [5, 8])], 13);
  drawTile(x, T.cheese, "#c98a12", [shred(["#f2b630", "#f7c84a", "#e8a21c", "#fbd76a"], 520, [10, 22], [2, 3.5])], 14);
  drawTile(x, T.pepper, "#24501a", [shred(["#3f8a2a", "#2f7320", "#5aa33a", "#7cbf55"], 240, [16, 30], [4, 6])], 15);
  drawTile(x, T.olive, "#121212", [{ c: ["#1d1b1a", "#2a2624", "#151413"], n: 150, w: [9, 12], h: [9, 12], shape: "ring" }], 16);
  drawTile(x, T.jalapeno, "#2c5c1d", [{ c: ["#4e9a2e", "#3d8524", "#62ad3c"], n: 120, w: [11, 15], h: [11, 15], shape: "ring" }, dice(["#e9e2b8"], 80, [1.5, 2.5])], 17);
  drawTile(x, T.redOnion, "#6d2a4a", [shred(["#a8487a", "#c56b98", "#e8d2de", "#8c3a66"], 300, [14, 26], [2, 3.5])], 18);
  drawTile(x, T.cucumber, "#2e5a24", [{ c: ["#cfe3a8", "#bcd690", "#d8eab8"], n: 90, w: [14, 18], h: [14, 18], shape: "disc" }], 19);
  drawTile(x, T.ham, "#c2707a", [dice(["#e9a3a8", "#f0b6b8", "#d98a92", "#f7d0cf"], 260, [6, 10])], 20);
  drawTile(x, T.chicken, "#b28a5a", [dice(["#e6d2b0", "#dcc29a", "#f0e2c6", "#c9a777"], 260, [6, 11])], 21);
  drawTile(x, T.mushroom, "#6b5a46", [{ c: ["#cdbba0", "#bba684", "#e0d2bb"], n: 110, w: [12, 18], h: [8, 12], shape: "disc" }], 22);
  drawTile(x, T.corn, "#b68a10", [{ c: ["#f6cf3a", "#f2c22a", "#fadb63"], n: 380, w: [4, 6], h: [4, 6], shape: "ellipse" }], 23);
  drawTile(x, T.pickle, "#4a5a18", [{ c: ["#8ea43c", "#a1b64e", "#7a9030"], n: 100, w: [13, 17], h: [13, 17], shape: "disc" }], 24);
  drawTile(x, T.banana, "#a8901c", [{ c: ["#e8d04a", "#f1dc6a", "#d8bc32"], n: 120, w: [11, 15], h: [11, 15], shape: "ring" }], 25);
  drawTile(x, T.spinach, "#1f3d14", [{ c: ["#2f5e1e", "#3d7228", "#28521a", "#4b8233"], n: 160, w: [16, 26], h: [10, 16], shape: "ellipse" }], 26);
  drawTile(x, T.chili, "#4a160c", [dice(["#7a2a14", "#8c3418", "#5e1e0e", "#a1461f"], 300, [4, 8]), { c: ["#6a3a1c", "#3a1a0c"], n: 80, w: [5, 7], h: [3, 5], shape: "ellipse" }], 27);
  drawTile(x, T.rice, "#e7e1d2", [{ c: ["#fbf9f2", "#f1ecdf", "#ffffff", "#e2dbc8"], n: 900, w: [5, 8], h: [2, 3], shape: "ellipse" }], 28);
  drawTile(x, T.braise, "#3a1c10", [dice(["#5c2e18", "#6e3a1e", "#4a2412", "#7e4a2a"], 120, [12, 20]), { c: ["#e08a2a", "#c76f1c"], n: 30, w: [6, 9], h: [6, 9], shape: "disc" }], 29);
  drawTile(x, T.mac, "#d29a22", [{ c: ["#f2c048", "#e8ad2e", "#f7d26a"], n: 260, w: [8, 12], h: [4, 6], shape: "ellipse" }, dice(["#9c5a14"], 40, [3, 6])], 30);
  drawTile(x, T.roastVeg, "#5a3a1a", [dice(["#d9822a", "#8a3a18", "#e2b24a", "#5e7a2a", "#b0402a"], 160, [10, 16])], 31);
  drawTile(x, T.mash, "#e6dcc0", [{ c: ["#f1ead4", "#ddd1b0", "#f7f1e0"], n: 140, w: [16, 30], h: [10, 18], shape: "ellipse" }], 32);
  drawTile(x, T.sauce, "#8c1c10", [{ c: ["#a32616", "#b8341e", "#7a160c"], n: 90, w: [14, 30], h: [10, 20], shape: "ellipse" }], 33);
  drawTile(x, T.soup, "#a8661e", [{ c: ["#e08a2a", "#3f7a28", "#f2d26a", "#c24a1e"], n: 90, w: [5, 9], h: [5, 9], shape: "rect" }], 34);
  const fibre: Piece = { c: ["rgba(90,60,30,0.10)", "rgba(255,236,200,0.10)"], n: 500, w: [2, 5], h: [1, 2], shape: "ellipse" };
  drawTile(x, T.kraft, "#ad8659", [fibre], 35);
  drawTile(x, T.kraftLabel, "#ad8659", [fibre], 36);
  {
    // A printed panel on the case: a white label with bars of plain type, no names.
    const ox = (T.kraftLabel % AT_COLS) * TILE;
    const oy = Math.floor(T.kraftLabel / AT_COLS) * TILE;
    x.fillStyle = "#efeae0";
    x.fillRect(ox + 22, oy + 30, 84, 56);
    x.fillStyle = "#2a2a2a";
    for (let i = 0; i < 6; i++) x.fillRect(ox + 28, oy + 36 + i * 8, i === 0 ? 60 : 34 + ((i * 17) % 36), i === 0 ? 5 : 3);
    x.fillStyle = "#1a5c9c";
    x.fillRect(ox + 22, oy + 92, 84, 6);
  }
  const can = (t: number, col: string) => {
    const ox = (t % AT_COLS) * TILE;
    const oy = Math.floor(t / AT_COLS) * TILE;
    x.fillStyle = col;
    x.fillRect(ox, oy, TILE, TILE);
    x.fillStyle = "#f4efe2";
    x.fillRect(ox, oy + 44, TILE, 40);
    x.fillStyle = "#222";
    for (let i = 0; i < 3; i++) x.fillRect(ox + 12, oy + 52 + i * 9, 50 + i * 12, 4);
  };
  can(T.canRed, "#a51e1a");
  can(T.canGreen, "#2c6a2a");
  for (const [t, col] of [[T.lidRed, "#b8261c"], [T.lidBlue, "#1e4f9c"], [T.lidGreen, "#2c7a3a"], [T.white, "#f2f1ec"]] as const) {
    x.fillStyle = col;
    x.fillRect((t % AT_COLS) * TILE, Math.floor(t / AT_COLS) * TILE, TILE, TILE);
  }
  drawBakery(x);
  const map = new THREE.CanvasTexture(c);
  map.colorSpace = THREE.SRGBColorSpace;
  map.anisotropy = 8;
  const m = new THREE.MeshStandardMaterial({ map, bumpMap: map, bumpScale: 0.6, metalness: 0, roughness: 0.62 });
  m.userData.uv = "own";
  goodsCache = m;
  return m;
}

/** The UV rectangle of tile `t`, inset half a texel per side so neighbours never bleed in. */
function tileUV(t: number, cols = AT_COLS, rows = AT_ROWS, px = TILE): [number, number, number, number] {
  const e = 2 / (cols * px);
  const u0 = (t % cols) / cols + e;
  const u1 = ((t % cols) + 1) / cols - e;
  const v1 = 1 - Math.floor(t / cols) / rows - e * (cols / rows);
  const v0 = 1 - (Math.floor(t / cols) + 1) / rows + e * (cols / rows);
  return [u0, v0, u1, v1];
}
/** Remap a geometry's 0..1 UVs into a tile rectangle. */
function intoTile(g: THREE.BufferGeometry, r: [number, number, number, number], su = 1, sv = 1): THREE.BufferGeometry {
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) {
    const u = (uv.getX(i) * su) % 1.0001;
    const v = (uv.getY(i) * sv) % 1.0001;
    uv.setXY(i, r[0] + (r[2] - r[0]) * Math.min(1, u), r[1] + (r[3] - r[1]) * Math.min(1, v));
  }
  return g;
}

/**
 * Printed plates: the engraved data plate every NSF unit carries (ratings and refrigerant, no maker's name), and the
 * graphics of a control panel. Lit readouts are a second atlas, unlit, so they glow.
 */
let plateCache: { plate: THREE.MeshStandardMaterial; led: THREE.MeshBasicMaterial } | null = null;
export const P = { data: 0, panel: 1, ice: 2, chill: 3 } as const;
export const L = { cool: 0, frz: 1, probe: 2, ice: 3, cool2: 4, chill: 5, set: 6, blank: 7 } as const;
function plateMaterials(): { plate: THREE.MeshStandardMaterial; led: THREE.MeshBasicMaterial } {
  if (plateCache) return plateCache;
  const c = document.createElement("canvas");
  c.width = 1024;
  c.height = 256;
  const x = c.getContext("2d")!;
  // Tile 0: the data plate, brushed with black-filled engraving.
  x.fillStyle = "#c9ccd0";
  x.fillRect(0, 0, 256, 256);
  x.strokeStyle = "#222";
  x.lineWidth = 3;
  x.strokeRect(8, 8, 240, 240);
  x.fillStyle = "#1c1c1c";
  x.font = "bold 20px monospace";
  x.fillText("COMMERCIAL REFRIGERATION", 18, 38);
  x.font = "15px monospace";
  const lines = ["MODEL  CR-2D-49S", "SERIAL 26-114870", "115 V  60 Hz  1 PH", "FULL LOAD  6.8 A", "REFRIGERANT  R-290", "CHARGE  150 g", "DESIGN P  HI 2.6 MPa", "          LO 1.0 MPa", "FOR USE BY TRAINED", "PERSONNEL ONLY"];
  lines.forEach((l, i) => x.fillText(l, 20, 66 + i * 18));
  // Tile 1: a control panel's printed graphics: keys, arrows, a bezel (no model shows it yet). The arrows and the power
  // mark are drawn, not typed: the page's sans-serif lacks them, and typing the first sent the browser through the
  // system's fonts for a fallback (about 17 ms on a desktop, in the middle of the hero kitchen's build).
  x.fillStyle = "#1a1c20";
  x.fillRect(256, 0, 256, 256);
  x.strokeStyle = "#8a9099";
  x.lineWidth = 2;
  for (let i = 0; i < 4; i++) x.strokeRect(276 + i * 58, 170, 44, 44);
  x.strokeRect(270, 24, 228, 120);
  {
    x.save();
    const kc = (i: number) => 298 + i * 58;
    x.fillStyle = x.strokeStyle = "#cfd4da";
    for (const [i, dir] of [[0, -1], [1, 1]] as const) {
      x.beginPath();
      x.moveTo(kc(i) - 10, 192 - dir * 8.5);
      x.lineTo(kc(i) + 10, 192 - dir * 8.5);
      x.lineTo(kc(i), 192 + dir * 8.5);
      x.closePath();
      x.fill();
    }
    x.font = "bold 22px sans-serif";
    x.textAlign = "center";
    x.fillText("SET", kc(2), 200, 38);
    x.lineWidth = 2.4;
    x.lineCap = "round";
    x.beginPath();
    x.arc(kc(3), 193, 7.5, -Math.PI / 2 + 0.7, -Math.PI / 2 - 0.7 + 2 * Math.PI);
    x.moveTo(kc(3), 183);
    x.lineTo(kc(3), 192);
    x.stroke();
    x.restore();
  }
  // Tile 2: the ice machine's front panel graphics.
  x.fillStyle = "#e7e9ec";
  x.fillRect(512, 0, 256, 256);
  x.fillStyle = "#3a3e44";
  x.font = "bold 18px sans-serif";
  x.fillText("ICE  ·  WASH  ·  OFF", 536, 120);
  x.fillRect(536, 140, 200, 3);
  // Tile 3: the blast chiller's mode legends.
  x.fillStyle = "#1a1c20";
  x.fillRect(768, 0, 256, 256);
  x.fillStyle = "#cfd4da";
  x.font = "bold 17px sans-serif";
  ["SOFT CHILL", "HARD CHILL", "FREEZE", "PROBE", "HOLD"].forEach((l, i) => x.fillText(l, 790, 40 + i * 40));
  const map = new THREE.CanvasTexture(c);
  map.colorSpace = THREE.SRGBColorSpace;
  map.anisotropy = 8;
  const pm = new THREE.MeshStandardMaterial({ map, metalness: 0.55, roughness: 0.4 });
  pm.userData.uv = "own";
  // The lit readouts: seven-segment style digits on a dark glass, eight tiles in a row.
  const d = document.createElement("canvas");
  d.width = 1024;
  d.height = 128;
  const y = d.getContext("2d")!;
  y.fillStyle = "#04070a";
  y.fillRect(0, 0, 1024, 128);
  const texts: [string, string][] = [["38°F", "#ff3b2e"], ["-4°F", "#ff3b2e"], ["68°C", "#3cff7a"], ["ICE", "#3cff7a"], ["36", "#3cff7a"], ["-18°", "#5cc8ff"], ["SET 3", "#ffb347"], ["", "#000"]];
  texts.forEach(([t, col], i) => {
    y.fillStyle = col;
    y.font = "bold 74px monospace";
    y.textAlign = "center";
    y.textBaseline = "middle";
    y.fillText(t, i * 128 + 64, 66, 118);
  });
  const lm = new THREE.CanvasTexture(d);
  lm.colorSpace = THREE.SRGBColorSpace;
  const led = new THREE.MeshBasicMaterial({ map: lm });
  led.userData.uv = "own";
  plateCache = { plate: pm, led };
  return plateCache;
}
/** A flat printed or lit face, `w` by `h`, facing +z at z = 0, showing a tile of the plate (4 across) or LED (8) atlas. */
function decal(w: number, h: number, tile: number, cols: number): THREE.BufferGeometry {
  const g = new THREE.PlaneGeometry(w, h);
  const r = tileUV(tile, cols, 1, cols === 4 ? 256 : 128);
  return intoTile(g, r);
}

/**
 * A cutting board's HDPE: off-white, satin, crossed by fine knife scores that roughen it where they run and catch
 * a little grey, a faint stain or two. UVs in metres (a tile is 0.3 m).
 */
let boardCache: THREE.MeshStandardMaterial | null = null;
function boardMaterial(): THREE.MeshStandardMaterial {
  if (boardCache) return boardCache;
  const N = 512;
  const c = document.createElement("canvas");
  c.width = c.height = N;
  const x = c.getContext("2d")!;
  const r = rng(91);
  x.fillStyle = "#d6d2c6";
  x.fillRect(0, 0, N, N);
  for (let i = 0; i < 6; i++) {
    const g = x.createRadialGradient(r() * N, r() * N, 0, r() * N, r() * N, 60 + r() * 120);
    g.addColorStop(0, "rgba(196,188,170,0.10)");
    g.addColorStop(1, "rgba(196,188,170,0)");
    x.fillStyle = g;
    x.fillRect(0, 0, N, N);
  }
  const rough = document.createElement("canvas");
  rough.width = rough.height = N;
  const y = rough.getContext("2d")!;
  y.fillStyle = "#8a8a8a";
  y.fillRect(0, 0, N, N);
  for (let i = 0; i < 900; i++) {
    const a = (r() < 0.6 ? 0 : Math.PI / 2) + (r() - 0.5) * 0.7;
    const len = 10 + r() * 70;
    const px = r() * N;
    const py = r() * N;
    for (const [ctx, col, lw] of [[x, `rgba(112,106,94,${0.08 + r() * 0.16})`, 0.8], [y, "rgba(235,235,235,0.5)", 1.2]] as const) {
      ctx.strokeStyle = col;
      ctx.lineWidth = lw;
      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(px + Math.cos(a) * len, py + Math.sin(a) * len);
      ctx.stroke();
    }
  }
  const map = new THREE.CanvasTexture(c);
  map.colorSpace = THREE.SRGBColorSpace;
  const rmap = new THREE.CanvasTexture(rough);
  for (const t of [map, rmap]) {
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
    t.repeat.set(1 / 0.3, 1 / 0.3);
    t.anisotropy = 8;
  }
  const m = new THREE.MeshStandardMaterial({ map, roughnessMap: rmap, roughness: 0.72, metalness: 0, color: 0xffffff });
  m.userData.uv = "own";
  return (boardCache = m);
}

/** The photoluminescent plastic of a walk-in's inside safety release: pale green, faintly aglow. */
let glowCache: THREE.MeshStandardMaterial | null = null;
function glowMaterial(): THREE.MeshStandardMaterial {
  return (glowCache ??= new THREE.MeshStandardMaterial({ color: 0xd4ecbd, emissive: 0x5d9a46, emissiveIntensity: 0.3, roughness: 0.42 }));
}

let pvcCache: THREE.MeshPhysicalMaterial | null = null;
/** Translucent PVC strip, the clear blue-tinged 2 mm film of a doorway curtain. */
function pvcMaterial(): THREE.MeshPhysicalMaterial {
  if (pvcCache) return pvcCache;
  const m = new THREE.MeshPhysicalMaterial({ color: 0xb8d4dc, metalness: 0, roughness: 0.18, transparent: true, opacity: 0.38, side: THREE.DoubleSide, depthWrite: false, clearcoat: 0.6 });
  m.userData.uv = "own";
  pvcCache = m;
  return m;
}

/** The wire mat of a shelf at room distance: chrome wires on a transparent ground, alpha-tested. One tile is 8 wires. */
let wireCache: THREE.MeshStandardMaterial | null = null;
function wireMatMaterial(): THREE.MeshStandardMaterial {
  if (wireCache) return wireCache;
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const x = c.getContext("2d")!;
  x.fillStyle = "#000";
  x.fillRect(0, 0, 128, 128);
  x.fillStyle = "#fff";
  for (let i = 0; i < 8; i++) x.fillRect(i * 16 + 6, 0, 4, 128);
  x.fillRect(0, 60, 128, 5);
  const a = new THREE.CanvasTexture(c);
  a.wrapS = a.wrapT = THREE.RepeatWrapping;
  const p = palette().chrome;
  const m = new THREE.MeshStandardMaterial({ color: p.color, metalness: p.metalness, roughness: 0.2, alphaMap: a, alphaTest: 0.5, side: THREE.DoubleSide });
  m.userData.uv = "own";
  wireCache = m;
  return m;
}

/** Fins and frost: the coil's aluminium, whitening and dulling as frost builds. One per model: its frost is live. */
function frostMaterial(): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: 0xb8bdc2, metalness: 0.7, roughness: 0.38, side: THREE.DoubleSide });
  m.userData.uv = "box";
  return m;
}
function setFrost(m: THREE.MeshStandardMaterial, f: number): void {
  const k = Math.max(0, Math.min(1, f));
  m.color.setRGB(0.72 + 0.24 * k, 0.74 + 0.23 * k, 0.76 + 0.22 * k);
  m.metalness = 0.7 * (1 - k) * (1 - k);
  m.roughness = 0.38 + 0.55 * k;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Refrigeration parts                                                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A finned-tube coil filling x 0..len, y 0..hgt, z -dep..0 of the current frame: aluminium plate fins threaded on
 * staggered 3/8 in copper tubes, galvanised tube sheets at the ends with the return bends outside them. A cooler's
 * evaporator is about 6 fins per inch, a freezer's 4 (frost needs the room); a condenser 12. Close up every fin is
 * there; from the aisle one in four.
 */
function finCoil(s: Shop, fin: THREE.Material, len: number, hgt: number, dep: number, fpi = 6, bends = true): void {
  const p = palette();
  if (s.detail === 0) {
    s.add(fin, new THREE.BoxGeometry(len, hgt, dep), [len / 2, hgt / 2, -dep / 2]);
    return;
  }
  const pitch = (0.0254 / fpi) * (s.fine(2) ? 1 : 5);
  const n = Math.max(2, Math.floor((len - 0.004) / pitch));
  const t = s.fine(2) ? 0.0003 : 0.0008;
  {
    const pos: number[] = [];
    const idx: number[] = [];
    const quad = (a: V3, b: V3, c: V3, d: V3) => {
      const k = pos.length / 3;
      pos.push(...a, ...b, ...c, ...d);
      idx.push(k, k + 1, k + 2, k, k + 2, k + 3);
    };
    for (let i = 0; i <= n; i++) {
      const x = 0.002 + (i * (len - 0.004)) / n;
      const [x0, x1] = [x - t / 2, x + t / 2];
      quad([x1, 0, 0], [x1, 0, -dep], [x1, hgt, -dep], [x1, hgt, 0]);
      quad([x0, 0, -dep], [x0, 0, 0], [x0, hgt, 0], [x0, hgt, -dep]);
      quad([x0, 0, 0], [x1, 0, 0], [x1, hgt, 0], [x0, hgt, 0]);
      quad([x1, 0, -dep], [x0, 0, -dep], [x0, hgt, -dep], [x1, hgt, -dep]);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setIndex(idx);
    g.computeVertexNormals();
    s.add(fin, g);
  }
  // Tubes: rows 25 mm apart up the face, ranks 22 mm deep, every other rank staggered half a row.
  const rows = Math.max(1, Math.floor(hgt / 0.0254));
  const ranks = Math.max(1, Math.floor(dep / 0.022));
  const tr = 0.0048;
  const tubeAt = (i: number, k: number): [number, number] => [(i + 0.5 + (k % 2 ? 0.25 : -0.25)) * (hgt / rows), -(k + 0.5) * (dep / ranks)];
  for (let k = 0; k < ranks; k++) {
    // The inner ranks hide deep in the pack; the face ranks are the ones seen between the fins.
    if (k > 0 && k < ranks - 1) continue;
    for (let i = 0; i < rows; i++) {
      const [y, z] = tubeAt(i, k);
      s.add(p.copper, new THREE.CylinderGeometry(tr, tr, len + 0.024, s.fine(2) ? 7 : 5, 1, true).rotateZ(Math.PI / 2), [len / 2, y, z]);
    }
  }
  // Tube sheets, and the U-bends outside them that join each tube to the next in its circuit.
  for (const x of [-0.006, len + 0.006]) s.add(p.ssDull, bx(s, 0.0016, hgt + 0.003, dep + 0.003, 0.0006), [x, hgt / 2, -dep / 2]);
  if (bends && s.fine(1))
    for (const [x, side] of [[-0.012, -1], [len + 0.012, 1]] as const)
      for (let k = 0; k < ranks; k++) {
        if (!s.fine(2) && k > 0) continue;
        for (let i = side < 0 ? 0 : 1; i + 1 < rows; i += 2) {
          const [y0, z] = tubeAt(i, k);
          const [y1] = tubeAt(i + 1, k);
          const R = (y1 - y0) / 2;
          const g = torus(s, R, tr, "z", Math.PI, 12, 6);
          // The half ring opens toward the tube sheet, bulging out of it.
          g.rotateZ(-Math.PI / 2);
          s.add(p.copper, g, [x, (y0 + y1) / 2, z], [0, 0, 0], [side, 1, 1]);
        }
      }
}

/**
 * An axial fan's blades and hub facing +z, about the joint's z axis: `blades` sheet-steel blades of constant pitch,
 * each twisted along its span as a real propeller fan's are, riveted to a pressed spider on the motor shaft.
 */
function axialFan(s: Shop, blade: THREE.Material, hub: THREE.Material, r: number, blades = 4, rh = r * 0.28): void {
  const nb = s.fine(2) ? (r < 0.07 ? 9 : 16) : 4;
  const sweep = 0.55;
  const P = r * 0.7; // pitch: the axial advance per turn
  for (let b = 0; b < blades; b++) {
    const a0 = (b / blades) * Math.PI * 2;
    const pts: THREE.Vector2[] = [];
    const lead = (k: number) => a0 + sweep * k * k;
    const span = (k: number) => 0.9 + 0.55 * k; // the blade widens toward its tip
    for (let i = 0; i <= nb; i++) {
      const k = i / nb;
      const rr = rh + (r - rh) * k;
      pts.push(new THREE.Vector2(rr * Math.cos(lead(k)), rr * Math.sin(lead(k))));
    }
    for (let i = nb; i >= 0; i--) {
      const k = i / nb;
      const rr = rh + (r - rh) * k;
      const a = lead(k) - span(k) * (0.4 + 0.05 * Math.sin(k * Math.PI));
      pts.push(new THREE.Vector2(rr * Math.cos(a), rr * Math.sin(a)));
    }
    const g = plate(s, new THREE.Shape(pts), 0.0014);
    const pos = g.attributes.position as THREE.BufferAttribute;
    const mid = a0 - 0.35;
    for (let i = 0; i < pos.count; i++) {
      let a = Math.atan2(pos.getY(i), pos.getX(i)) - mid;
      a = Math.atan2(Math.sin(a), Math.cos(a));
      pos.setZ(i, pos.getZ(i) - a * (P / (Math.PI * 2)));
    }
    g.deleteAttribute("normal");
    g.computeVertexNormals();
    s.add(blade, g);
  }
  // The profile runs bottom to top, so the lathe winds outward and reads from outside.
  const hg = lathe(s, [
    [0, -0.01],
    [rh * 0.9, -0.01],
    [rh, -0.006],
    [rh * 0.95, 0.008],
    [rh * 0.5, 0.017],
    [0, 0.018],
  ], s.fine(2) ? (r < 0.07 ? 18 : 28) : 10);
  hg.rotateX(Math.PI / 2);
  s.add(hub, hg);
  if (s.fine(2))
    for (let b = 0; b < blades; b++) {
      const a = (b / blades) * Math.PI * 2 - 0.1;
      screwLite(s, hub, [Math.cos(a) * rh * 0.8, Math.sin(a) * rh * 0.8, 0.006], 0.0025);
    }
}

/**
 * A wire fan guard facing +z, dished by `dish`: concentric rings welded to radial spokes, the outer ring on four
 * mounting tabs. Zinc-plated steel wire about 2.5 mm.
 */
function fanGuard(s: Shop, mat: THREE.Material, r: number, dish = 0.02, wr = 0.0013): void {
  if (s.detail === 0) return;
  const rings = s.fine(2) ? Math.max(3, Math.round(r / 0.027)) : r > 0.12 ? 2 : 1;
  const dz = (rho: number) => dish * (1 - (rho / r) * (rho / r));
  for (let i = 1; i <= rings; i++) {
    const rho = (r * i) / rings;
    s.add(mat, torus(s, rho, i === rings ? wr * 1.8 : wr, "z", Math.PI * 2, s.fine(2) ? Math.max(12, Math.round((32 * rho) / Math.max(r, 0.1))) : 12, 4), [0, 0, dz(rho)]);
  }
  const spokes = s.fine(2) ? (r < 0.07 ? 8 : 12) : 4;
  for (let i = 0; i < spokes; i++) {
    const a = (i / spokes) * Math.PI * 2;
    let prev: V3 | null = null;
    for (let k = 0; k <= 3; k++) {
      const rho = r * (0.12 + 0.88 * (k / 3));
      const pt: V3 = [Math.cos(a) * rho, Math.sin(a) * rho, dz(rho) + wr * 1.6];
      if (prev) wire(s, mat, prev, pt, wr * 1.1, 5);
      prev = pt;
    }
  }
  s.add(mat, cyl(s, r * 0.14, 0.003, "z", 20), [0, 0, dish + wr]);
  for (let i = 0; i < 4; i++) {
    const a = (i / 4) * Math.PI * 2 + Math.PI / 4;
    s.add(mat, new THREE.BoxGeometry(0.016, 0.03, 0.0015), [Math.cos(a) * (r + 0.008), Math.sin(a) * (r + 0.008), 0.0008], [0, 0, a + Math.PI / 2]);
  }
}

/** A louvred grille in the current frame's xy, facing +z: `n` slats angled down, in a folded frame, `w` by `h`. */
function louvres(s: Shop, mat: THREE.Material, w: number, h: number, n: number, depth = 0.012, back: THREE.Material | null = null): void {
  const p = palette();
  if (back) s.add(back, new THREE.PlaneGeometry(w - 0.004, h - 0.004), [0, 0, -depth - 0.004]);
  const pitch = h / n;
  for (let i = 0; i < n; i++) {
    const y = -h / 2 + pitch * (i + 0.5);
    s.add(mat, bx(s, w - 0.006, pitch * 1.15, 0.0012, 0.0005), [0, y, -depth / 2], [-0.75, 0, 0]);
  }
  for (const sx of [-1, 1]) s.add(mat, bx(s, 0.008, h, depth + 0.004, 0.001), [(sx * (w - 0.008)) / 2, 0, -depth / 2]);
  for (const sy of [-1, 1]) s.add(mat, bx(s, w, 0.008, depth + 0.004, 0.001), [0, (sy * (h - 0.008)) / 2, -depth / 2]);
  void p;
}

/**
 * A hermetic compressor standing at the current origin: a welded steel dome on three rubber-grommeted feet, the
 * terminal box on its side, the suction and discharge stubs in copper. `r` is the shell's radius.
 */
function compressor(s: Shop, r = 0.08, h = 0.24, n = 32): void {
  const p = palette();
  s.add(p.enamel, lathe(s, [
    [0, 0.03],
    [r * 0.9, 0.03],
    [r, 0.05],
    [r, h * 0.55],
    [r * 1.02, h * 0.57],
    [r, h * 0.6],
    [r * 0.95, h * 0.8],
    [r * 0.7, h * 0.95],
    [0, h],
  ], n));
  // The base plate and its feet on grommets.
  for (let i = 0; i < 3; i++) {
    const a = (i / 3) * Math.PI * 2 + 0.5;
    s.add(p.enamel, bx(s, r * 0.7, 0.006, 0.04, 0.002), [Math.cos(a) * r * 0.75, 0.033, Math.sin(a) * r * 0.75], [0, -a, 0]);
    s.add(p.rubber, cyl(s, 0.013, 0.022, "y", 14), [Math.cos(a) * r * 1.1, 0.02, Math.sin(a) * r * 1.1]);
    if (s.fine(2)) s.add(p.ssBright, cyl(s, 0.006, 0.006, "y", 8), [Math.cos(a) * r * 1.1, 0.034, Math.sin(a) * r * 1.1]);
  }
  // Terminal box with its cover and the leads out.
  s.add(p.plastic, bx(s, r * 0.7, r * 0.6, r * 0.35, 0.006), [0, h * 0.5, r * 0.95]);
  if (s.fine(2)) s.add(p.plastic, cyl(s, 0.006, 0.05, "y", 8), [r * 0.18, h * 0.5 + r * 0.45, r * 0.95]);
  // Suction (larger) and discharge stubs.
  pipe(s, p.copper, [[-r * 0.6, h * 0.45, r * 0.6], [-r * 1.1, h * 0.45, r * 1.0], [-r * 1.1, h * 0.95, r * 1.0]], 0.0079, 0.02, 10);
  pipe(s, p.copper, [[r * 0.3, h * 0.9, -r * 0.3], [r * 0.3, h * 1.12, -r * 0.3], [r * 0.9, h * 1.12, -r * 0.3]], 0.0048, 0.015, 8);
}

/** A data plate, riveted, `w` by `h`, facing +z at the current origin. */
function dataPlate(s: Shop, w = 0.09, h = 0.09): void {
  const p = palette();
  if (!s.fine(1)) return;
  // The printed face stands half a millimetre proud of its backing plate, so the two never share a depth.
  s.add(plateMaterials().plate, decal(w, h, P.data, 4), [0, 0, 0.0011]);
  s.add(p.ssDull, new THREE.BoxGeometry(w + 0.004, h + 0.004, 0.0006), [0, 0, 0.0003]);
  if (s.fine(2)) for (const [sx, sy] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) screwLite(s, p.ssBright, [(sx * (w - 0.006)) / 2, (sy * (h - 0.006)) / 2, 0.0012], 0.0018);
}

/**
 * A digital thermometer facing +z at the current origin: a black ABS bezel round a lit three-digit readout, and its
 * two keys below.
 */
function thermometer(s: Shop, tile: number, w = 0.07, h = 0.04, live: THREE.Material | null = null): void {
  const p = palette();
  s.add(p.plastic, bx(s, w, h, 0.012, 0.003), [0, 0, 0.006]);
  if (live) s.add(live, decalFull(w * 0.72, h * 0.5), [0, h * 0.1, 0.0125]);
  else if (s.fine(1)) s.add(plateMaterials().led, decal(w * 0.72, h * 0.5, tile, 8), [0, h * 0.1, 0.0125]);
  else s.add(p.screen, new THREE.BoxGeometry(w * 0.7, h * 0.5, 0.001), [0, h * 0.1, 0.0125]);
  if (s.fine(2)) for (const sx of [-1, 1]) s.add(p.charcoal, bx(s, 0.008, 0.005, 0.002, 0.001), [sx * 0.012, -h * 0.32, 0.0125]);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Pans, tubs and wires                                                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/** A part's box: edges rounded close up, where they catch the light; a plain box from the aisle, where they cannot. */
function bx(s: Shop, w: number, h: number, d: number, r = 0.0015): THREE.BufferGeometry {
  return s.fine(2) ? rbox(s, w, h, d, r) : new THREE.BoxGeometry(w, h, d);
}

/**
 * A box without the faces named (+x, -x, +y, -y, +z, -z): a panel's foam core shows only at its cut edges, and a
 * face left lying a millimetre behind a skin would fight it for depth from across a room.
 */
function openBox(w: number, h: number, d: number, drop: ("px" | "nx" | "py" | "ny" | "pz" | "nz")[]): THREE.BufferGeometry {
  const g = new THREE.BoxGeometry(w, h, d);
  const order = ["px", "nx", "py", "ny", "pz", "nz"];
  const src = g.index!.array;
  const keep: number[] = [];
  order.forEach((f, i) => {
    if (drop.includes(f as never)) return;
    for (let k = i * 6; k < i * 6 + 6; k++) keep.push(src[k]!);
  });
  g.setIndex(keep);
  g.clearGroups();
  return g;
}

/** A shape through corners [x, y, r], each corner rounded by its own r (0 sharp), closed. */
function roundedPoly(pts: readonly (readonly [number, number, number])[]): THREE.Shape {
  const sh = new THREE.Shape();
  const n = pts.length;
  const along = (i: number, j: number, d: number): [number, number] => {
    const a = pts[i]!;
    const b = pts[j]!;
    const L = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
    return [a[0] + ((b[0] - a[0]) * d) / L, a[1] + ((b[1] - a[1]) * d) / L];
  };
  for (let i = 0; i < n; i++) {
    const [x, y, r] = pts[i]!;
    const prev = (i + n - 1) % n;
    const next = (i + 1) % n;
    if (r > 0) {
      const a = along(i, prev, r);
      const b = along(i, next, r);
      if (i === 0) sh.moveTo(a[0], a[1]);
      else sh.lineTo(a[0], a[1]);
      sh.quadraticCurveTo(x, y, b[0], b[1]);
    } else if (i === 0) sh.moveTo(x, y);
    else sh.lineTo(x, y);
  }
  sh.closePath();
  return sh;
}

/**
 * A plate cut from `shape` (its own xy), `t` thick along +z from 0, its edges eased by `bev` without growing past the
 * outline: a pressed or cast part whose edges catch the light.
 */
function bevelPlate(s: Shop, shape: THREE.Shape, t: number, bev: number, curve = 8): THREE.BufferGeometry {
  const b = s.detail > 0 ? Math.min(bev, t / 2 - 1e-5) : 0;
  const g = new THREE.ExtrudeGeometry(shape, {
    depth: Math.max(1e-5, t - 2 * b),
    bevelEnabled: b > 0,
    bevelThickness: b,
    bevelSize: b,
    bevelOffset: -b,
    bevelSegments: s.fine(2) ? 2 : 1,
    curveSegments: s.seg(curve, 2),
  });
  if (b > 0) g.translate(0, 0, b);
  return g;
}

/**
 * A formed angle trimming an outside corner. In its own frame the corner runs along z from 0 to `len`; one face of
 * the corner is the plane y = 0 reaching toward +x, the other the plane x = 0 reaching toward -y. The sheet, `t`
 * thick, stands `gap` off both faces, `legA` along the first and `legB` down the second, its heel rounded to `R`.
 */
function angleBar(s: Shop, legA: number, legB: number, len: number, t = 0.0012, R = 0.005, gap = 0.0006): THREE.BufferGeometry {
  const sh = new THREE.Shape();
  const cx = -gap - t + R;
  const cy = gap + t - R;
  sh.moveTo(legA, gap);
  sh.lineTo(legA, gap + t);
  sh.lineTo(cx, gap + t);
  sh.absarc(cx, cy, R, Math.PI / 2, Math.PI, false);
  sh.lineTo(-gap - t, -legB);
  sh.lineTo(-gap, -legB);
  sh.lineTo(-gap, cy);
  sh.absarc(cx, cy, R - t, Math.PI, Math.PI / 2, true);
  sh.lineTo(legA, gap);
  return new THREE.ExtrudeGeometry(sh, { depth: len, bevelEnabled: false, curveSegments: s.seg(6, 2) });
}

/** A polyline in a plane with each interior corner rounded to `r` by `k` segments. */
function roundPath(pts: readonly (readonly [number, number])[], r: number, k: number): [number, number][] {
  const out: [number, number][] = [[pts[0]![0], pts[0]![1]]];
  for (let i = 1; i + 1 < pts.length; i++) {
    const [ax, ay] = pts[i - 1]!;
    const [px, py] = pts[i]!;
    const [bx2, by2] = pts[i + 1]!;
    const la = Math.hypot(px - ax, py - ay);
    const lb = Math.hypot(bx2 - px, by2 - py);
    const rr = Math.min(r, la * 0.45, lb * 0.45);
    const a: [number, number] = [px + ((ax - px) * rr) / la, py + ((ay - py) * rr) / la];
    const b: [number, number] = [px + ((bx2 - px) * rr) / lb, py + ((by2 - py) * rr) / lb];
    for (let j = 0; j <= k; j++) {
      const t = j / k;
      const u = 1 - t;
      out.push([u * u * a[0] + 2 * u * t * px + t * t * b[0], u * u * a[1] + 2 * u * t * py + t * t * b[1]]);
    }
  }
  const last = pts[pts.length - 1]!;
  out.push([last[0], last[1]]);
  return out;
}

/**
 * A closed section swept along a polyline lying in the plane z = z0 of the current frame: `prof` points are (u, v),
 * u across the path in its plane, v out of the plane along `out`; mitred at each vertex, its ends left open (they
 * always die into another part). Wound to face away from the section's centre. A gasket round a door, a bead.
 */
function sweepPath(path: readonly (readonly [number, number])[], prof: readonly (readonly [number, number])[], z0: number, out: 1 | -1): THREE.BufferGeometry {
  const n = path.length;
  const m = prof.length;
  const pos: number[] = [];
  const unit = (x: number, y: number): [number, number] => {
    const L = Math.hypot(x, y) || 1;
    return [x / L, y / L];
  };
  for (let i = 0; i < n; i++) {
    const P = path[i]!;
    const t0 = i > 0 ? unit(P[0] - path[i - 1]![0], P[1] - path[i - 1]![1]) : null;
    const t1 = i + 1 < n ? unit(path[i + 1]![0] - P[0], path[i + 1]![1] - P[1]) : null;
    const [tx, ty] = unit((t0?.[0] ?? 0) + (t1?.[0] ?? 0), (t0?.[1] ?? 0) + (t1?.[1] ?? 0));
    const nx = -ty;
    const ny = tx;
    const sg = t1 ?? t0!;
    const c = Math.max(0.35, nx * -sg[1] + ny * sg[0]);
    for (const [u, v] of prof) pos.push(P[0] + (nx * u) / c, P[1] + (ny * u) / c, z0 + out * v);
  }
  const idx: number[] = [];
  for (let i = 0; i + 1 < n; i++)
    for (let k = 0; k < m; k++) {
      const k2 = (k + 1) % m;
      const a = i * m + k;
      const b = i * m + k2;
      const c = (i + 1) * m + k;
      const d = (i + 1) * m + k2;
      idx.push(a, c, b, b, c, d);
    }
  // Face outward: test the first quad against the direction from the section's centre to it.
  const V = (q: number) => new THREE.Vector3(pos[q * 3]!, pos[q * 3 + 1]!, pos[q * 3 + 2]!);
  const cu = prof.reduce((a, q) => a + q[0], 0) / m;
  const cv = prof.reduce((a, q) => a + q[1], 0) / m;
  void cu;
  const A = V(0);
  const C = V(m);
  const B = V(1);
  const nrm = new THREE.Vector3().subVectors(C, A).cross(new THREE.Vector3().subVectors(B, A));
  const mid = A.clone().add(B).multiplyScalar(0.5);
  const ctr = new THREE.Vector3(path[0]![0], path[0]![1], z0 + out * cv);
  if (nrm.dot(mid.sub(ctr)) < 0)
    for (let i = 0; i < idx.length; i += 3) {
      const tmp = idx[i + 1]!;
      idx[i + 1] = idx[i + 2]!;
      idx[i + 2] = tmp;
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** A pan-head screw's head facing +z from z = 0, eight-sided: the kit's screw is too dear for a walk-in's hundreds. */
function screwLite(s: Shop, mat: THREE.Material, at: V3, r = 0.003, rot?: V3): void {
  if (!s.fine(2)) return;
  const g = new THREE.LatheGeometry([new THREE.Vector2(0, 0), new THREE.Vector2(r, 0), new THREE.Vector2(r * 0.8, r * 0.45), new THREE.Vector2(0, r * 0.55)], 8);
  g.rotateX(Math.PI / 2);
  s.add(mat, g, at, rot);
}

/**
 * A pipe along a polyline: straight cylinders joined by torus bends of radius `bend`, as copper is bent with a
 * lever bender. Far fewer triangles than a swept tube for long runs.
 */
function pipe(s: Shop, mat: THREE.Material, pts: readonly V3[], r: number, bend = r * 4, radial = 12): void {
  const v = pts.map((q) => new THREE.Vector3(q[0], q[1], q[2]));
  const rad = s.seg(radial, 5);
  let from = v[0]!.clone();
  const straight = (a: THREE.Vector3, b: THREE.Vector3) => {
    const dd = b.clone().sub(a);
    const len = dd.length();
    if (len < 1e-4) return;
    const g = new THREE.CylinderGeometry(r, r, len, rad, 1, true);
    g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), dd.normalize()));
    g.translate((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
    s.add(mat, g);
  };
  for (let i = 1; i < v.length; i++) {
    const p1 = v[i]!;
    if (i === v.length - 1) {
      straight(from, p1);
      break;
    }
    const u = p1.clone().sub(v[i - 1]!).normalize();
    const w = v[i + 1]!.clone().sub(p1).normalize();
    const phi = Math.acos(Math.max(-1, Math.min(1, u.dot(w))));
    if (phi < 1e-3) continue;
    const cut0 = Math.min(bend * Math.tan(phi / 2), p1.distanceTo(v[i - 1]!) * 0.49, p1.distanceTo(v[i + 1]!) * 0.49);
    const R = cut0 / Math.tan(phi / 2);
    const a = p1.clone().addScaledVector(u, -cut0);
    const b = p1.clone().addScaledVector(w, cut0);
    straight(from, a);
    const m = w.clone().addScaledVector(u, -w.dot(u)).normalize();
    const C = a.clone().addScaledVector(m, R);
    const xA = m.clone().negate();
    const zA = new THREE.Vector3().crossVectors(xA, u);
    const g = new THREE.TorusGeometry(R, r, rad, Math.max(3, Math.round((phi / (Math.PI / 2)) * s.seg(8, 3))), phi);
    g.applyMatrix4(new THREE.Matrix4().makeBasis(xA, u, zA).setPosition(C));
    s.add(mat, g);
    from = b;
  }
}

/** A straight round wire from a to b, open-ended (its ends always land in another part), few-sided. */
function wire(s: Shop, mat: THREE.Material, a: V3, b: V3, r: number, radial?: number): void {
  const d = new THREE.Vector3(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const len = d.length();
  if (len < 1e-5) return;
  const g = new THREE.CylinderGeometry(r, r, len, radial ?? (s.fine(2) ? 6 : 4), 1, true);
  const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.normalize());
  const e = new THREE.Euler().setFromQuaternion(q, "XYZ");
  s.add(mat, g, [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2], [e.x, e.y, e.z]);
}

/** A rounded rectangle's outline in plan (x, z), `k` points a corner, counter-clockwise seen from above. */
function outline(w: number, d: number, r: number, k: number): [number, number][] {
  const rr = Math.max(0, Math.min(r, w / 2 - 1e-5, d / 2 - 1e-5));
  const out: [number, number][] = [];
  // Each corner's quarter arc about its own centre, swept the same way round.
  const corners: [number, number, number][] = [
    [w / 2 - rr, -d / 2 + rr, 0],
    [-w / 2 + rr, -d / 2 + rr, -Math.PI / 2],
    [-w / 2 + rr, d / 2 - rr, -Math.PI],
    [w / 2 - rr, d / 2 - rr, Math.PI / 2],
  ];
  for (const [cx, cz, a0] of corners)
    for (let i = 0; i <= k; i++) {
      const a = a0 - (i / k) * (Math.PI / 2);
      out.push([cx + rr * Math.cos(a), cz + rr * Math.sin(a)]);
    }
  return out;
}

export interface Ring {
  w: number;
  d: number;
  r: number;
  y: number;
}
/**
 * A surface lofted through rounded-rectangle rings in order (a pan's outer wall up, its rim, its inner wall down),
 * with flat caps on the first and last rings facing up (+1) or down (-1) when asked. A ring listed twice makes a
 * sharp crease.
 */
function loft(s: Shop, rings: Ring[], capFirst: -1 | 0 | 1, capLast: -1 | 0 | 1, k?: number): THREE.BufferGeometry {
  const kk = k ?? s.seg(6, 2);
  const loops = rings.map((g) => outline(g.w, g.d, g.r, kk));
  const n = loops[0]!.length;
  const pos: number[] = [];
  const idx: number[] = [];
  loops.forEach((l, j) => l.forEach(([x, z]) => pos.push(x, rings[j]!.y, z)));
  for (let j = 0; j + 1 < loops.length; j++)
    for (let i = 0; i < n; i++) {
      const a = j * n + i;
      const b = j * n + ((i + 1) % n);
      const c = (j + 1) * n + i;
      const d = (j + 1) * n + ((i + 1) % n);
      idx.push(a, b, c, b, d, c);
    }
  const cap = (j: number, up: number) => {
    const base = pos.length / 3;
    pos.push(0, rings[j]!.y, 0);
    const l = loops[j]!;
    l.forEach(([x, z]) => pos.push(x, rings[j]!.y, z));
    for (let i = 0; i < n; i++) {
      const a = base + 1 + i;
      const b = base + 1 + ((i + 1) % n);
      if (up > 0) idx.push(base, a, b);
      else idx.push(base, b, a);
    }
  };
  // Caps face up (+1) or down (-1).
  if (capFirst) cap(0, capFirst);
  if (capLast) cap(rings.length - 1, capLast);
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** Smooth value noise in the plane, about -1..1, a lattice of unit cells. */
function vnoise(x: number, z: number, seed: number): number {
  const h = (i: number, j: number) => {
    let n = (i * 374761393 + j * 668265263 + seed * 982451653) | 0;
    n = Math.imul(n ^ (n >>> 13), 1274126177);
    return ((n ^ (n >>> 16)) >>> 0) / 2147483648 - 1;
  };
  const i = Math.floor(x);
  const j = Math.floor(z);
  const fx = x - i;
  const fz = z - j;
  const u = fx * fx * (3 - 2 * fx);
  const v = fz * fz * (3 - 2 * fz);
  const a = h(i, j) + (h(i + 1, j) - h(i, j)) * u;
  const b = h(i, j + 1) + (h(i + 1, j + 1) - h(i, j + 1)) * u;
  return a + (b - a) * v;
}

/** How each food lies in a pan: the relief's height and the size of its lumps, in metres. */
function foodRelief(tile: number): [number, number] {
  switch (tile) {
    case T.soup:
    case T.sauce:
      return [0.0008, 0.06];
    case T.mash:
      return [0.0055, 0.055];
    case T.rice:
    case T.corn:
      return [0.003, 0.035];
    case T.lettuce:
    case T.spinach:
    case T.redOnion:
      return [0.007, 0.04];
    case T.cheese:
    case T.pepper:
      return [0.005, 0.034];
    case T.roastVeg:
    case T.braise:
      return [0.0065, 0.034];
    case T.chili:
    case T.mac:
      return [0.004, 0.03];
    case T.frosting:
      return [0.0015, 0.05];
    default:
      return [0.0042, 0.03];
  }
}

/**
 * A food surface filling a pan's inside at height y: heaped toward the middle as a ladle leaves it and lumped by the
 * pieces it is made of (diced, shredded, spooned), trimmed to the pan's rounded corners and settling to the fill line
 * at the wall; its UVs on a food tile. `heap` is the crown in metres.
 */
function foodSurface(s: Shop, w: number, d: number, r: number, heap: number, tile: number, seed: number): THREE.BufferGeometry {
  const [amp, lam] = foodRelief(tile);
  // Fine enough to carry its lumps, no finer.
  const step = Math.max(0.011, Math.min(0.02, lam * 0.4));
  const nx = s.fine(2) ? Math.max(6, Math.min(44, Math.round(w / step))) : s.fine(1) ? 4 : 1;
  const nz = s.fine(2) ? Math.max(6, Math.min(56, Math.round(d / step))) : s.fine(1) ? 3 : 1;
  const g = new THREE.PlaneGeometry(w, d, nx, nz);
  g.rotateX(-Math.PI / 2);
  const pos = g.attributes.position as THREE.BufferAttribute;
  const rnd = rng(seed);
  const ox = rnd() * 50;
  const oz = rnd() * 50;
  // Where the ladle last dropped: the heap's peak sits off centre.
  const px = (rnd() - 0.5) * w * 0.3;
  const pz = (rnd() - 0.5) * d * 0.3;
  for (let i = 0; i < pos.count; i++) {
    let x = pos.getX(i);
    let z = pos.getZ(i);
    const cx = Math.max(-(w / 2 - r), Math.min(w / 2 - r, x));
    const cz = Math.max(-(d / 2 - r), Math.min(d / 2 - r, z));
    const dl = Math.hypot(x - cx, z - cz);
    if (dl > r) {
      x = cx + ((x - cx) * r) / dl;
      z = cz + ((z - cz) * r) / dl;
    }
    const edge = Math.abs(Math.abs(x) - w / 2) < 1e-4 || Math.abs(Math.abs(z) - d / 2) < 1e-4 || dl > r * 0.98;
    // Distance in from the wall, 0 at it: the food climbs the wall a little and is lumpiest well in from it.
    const inX = w / 2 - Math.abs(x);
    const inZ = d / 2 - Math.abs(z);
    const wall = Math.min(1, Math.min(inX, inZ) / 0.025);
    const qx = (x - px) / (w / 2 + Math.abs(px));
    const qz = (z - pz) / (d / 2 + Math.abs(pz));
    const crown = Math.max(0, 1 - (qx * qx + qz * qz) * 0.85);
    let y = 0;
    if (!edge && s.fine(1)) {
      const n1 = vnoise((x + ox) / lam, (z + oz) / lam, seed);
      const n2 = vnoise((x + ox) / (lam * 0.45), (z + oz) / (lam * 0.45), seed + 7);
      const lumps = s.fine(2) ? amp * (0.65 * n1 + 0.35 * n2) : 0;
      y = heap * crown * (0.85 + 0.3 * rnd()) * wall + lumps * wall + 0.0015 * (1 - wall);
    }
    pos.setXYZ(i, x, y, z);
  }
  g.computeVertexNormals();
  // A tile is about 0.15 m of food; repeat across a big pan, mirrored so seams match.
  const reps = Math.max(1, Math.round(Math.max(w, d) / 0.16));
  const uv = g.attributes.uv as THREE.BufferAttribute;
  const [u0, v0, u1, v1] = tileUV(tile);
  for (let i = 0; i < uv.count; i++) {
    const fu = (uv.getX(i) * reps * (w / Math.max(w, d))) % 2;
    const fv = (uv.getY(i) * reps * (d / Math.max(w, d))) % 2;
    const mu = fu > 1 ? 2 - fu : fu;
    const mv = fv > 1 ? 2 - fv : fv;
    uv.setXY(i, u0 + (u1 - u0) * mu, v0 + (v1 - v0) * mv);
  }
  return g;
}

/**
 * A stainless steam-table pan (hotel pan), top `w` by `d`, `h` deep, standing on its base at the current origin: a
 * drawn shell tapering 4° to its floor, the rim turned out into a flange that hangs it in a well, `fill` 0..1 of it
 * food from `tile` (or empty when tile < 0), or covered by a flat lid.
 */
export function hotelPan(s: Shop, w: number, d: number, h: number, tile: number, fill = 0.85, lid = false, seed = 1): void {
  const p = palette();
  const t = 0.0008;
  const tap = h * 0.07;
  const r1 = Math.min(0.03, Math.min(w, d) * 0.12);
  const r0 = r1 * 0.75;
  const f = 0.011; // flange
  if (s.detail === 0) {
    s.add(p.ss, new THREE.BoxGeometry(w, h, d), [0, h / 2, 0]);
    return;
  }
  const far = !s.fine(2);
  const rings: Ring[] = far
    ? [
        { w: w - 2 * tap, d: d - 2 * tap, r: r0, y: 0 },
        { w: w + 2 * f, d: d + 2 * f, r: r1 + f, y: h },
        { w: w - 2 * t, d: d - 2 * t, r: r1 - t, y: h },
        { w: w - 2 * tap - 2 * t, d: d - 2 * tap - 2 * t, r: r0 - t, y: t },
      ]
    : [
    { w: w - 2 * tap, d: d - 2 * tap, r: r0, y: 0 },
    { w: w - 2 * tap, d: d - 2 * tap, r: r0, y: 0 },
    { w, d, r: r1, y: h - 0.004 },
    { w: w + 2 * f, d: d + 2 * f, r: r1 + f, y: h - 0.001 },
    { w: w + 2 * f + 0.002, d: d + 2 * f + 0.002, r: r1 + f, y: h - 0.006 },
    { w: w + 2 * f - 0.001, d: d + 2 * f - 0.001, r: r1 + f, y: h - 0.0068 },
    { w: w + 2 * f - 0.002, d: d + 2 * f - 0.002, r: r1 + f, y: h + 0.0008 },
    { w: w - 2 * t, d: d - 2 * t, r: r1 - t, y: h - 0.003 },
    { w: w - 2 * tap - 2 * t, d: d - 2 * tap - 2 * t, r: r0 - t, y: t },
  ];
  s.add(p.ss, loft(s, rings, -1, 1, s.fine(2) ? 4 : 1));
  if (lid) {
    // A flat lid with a stiffening panel pressed up in it and a bent wire handle... steam-table lids are solid.
    const lw = w + 2 * f;
    const ld = d + 2 * f;
    s.add(p.ss, loft(s, far ? [
      { w: lw, d: ld, r: r1 + f, y: h + 0.0012 },
      { w: lw - 0.04, d: ld - 0.04, r: r1, y: h + 0.011 },
    ] : [
      { w: lw, d: ld, r: r1 + f, y: h + 0.0012 },
      { w: lw, d: ld, r: r1 + f, y: h + 0.004 },
      { w: lw - 0.03, d: ld - 0.03, r: r1, y: h + 0.006 },
      { w: lw - 0.04, d: ld - 0.04, r: r1, y: h + 0.011 },
    ], -1, 1, far ? 1 : 3));
    if (!far) s.add(p.ssBright, bar(s, [[-0.05, h + 0.011, 0], [-0.04, h + 0.03, 0], [0.04, h + 0.03, 0], [0.05, h + 0.011, 0]], 0.003, 0.008, 8));
  } else if (tile >= 0 && fill > 0 && s.detail > 0) {
    const fy = t + (h - 0.004 - t) * fill;
    const k = fy / h;
    const fw = w - 2 * tap * (1 - k) - 2 * t - 0.0015;
    const fd = d - 2 * tap * (1 - k) - 2 * t - 0.0015;
    s.add(goodsMaterial(), foodSurface(s, fw, fd, Math.max(0.004, r1 * 0.9 - t), Math.min(0.016, h * 0.18), tile, seed), [0, fy, 0]);
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Wire shelving                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

export interface ShelfOpts {
  /** Shelf heights (the top of the wire mat) above the floor. Default four, 0.2 m to 0.13 m under the post tops. */
  tiers?: number[];
  /** Cut the posts' inch grooves. On by default close up; a walk-in turns them off to share its budget. */
  grooves?: boolean;
  /** The shops (joints, so they can be hidden) the stock goes in, item by item in turn; none when absent. */
  stock?: Shop[] | null;
  /** How full the shelves are, 0..1. */
  fill?: number;
  /** Cases among the stock (default true); without, tubs and cans only. */
  cartons?: boolean;
  /** Fewer-sided wires and collars, for the many units in a walk-in. */
  lite?: boolean;
  seed?: number;
}

/**
 * NSF wire shelving, `L` long (x 0..L), `D` deep (z -D..0), on posts `H` tall, in the current frame: four 1 in
 * posts grooved every inch, each shelf clamped on tapered split sleeves snapped into a groove; the mats chrome
 * wires front to back on longitudinal supports, their edges zig-zag trusses that carry the load.
 */
export function wireShelving(s: Shop, L: number, D: number, H: number, o: ShelfOpts = {}): void {
  const p = palette();
  const pr = 0.0127;
  const c = 0.017;
  const top = H - 0.13;
  const tiers = o.tiers ?? [0.2, 0.2 + (top - 0.2) / 3, 0.2 + ((top - 0.2) * 2) / 3, top];
  const posts: [number, number][] = [[c, -c], [L - c, -c], [c, -D + c], [L - c, -D + c]];
  const lite = o.lite ?? false;
  const grooves = (o.grooves ?? true) && s.fine(2);
  const cs = lite ? 10 : 18;
  for (const [x, z] of posts) {
    if (grooves) {
      const prof: [number, number][] = [[0, 0.028], [pr, 0.028]];
      for (let gy = 0.05; gy < H - 0.03; gy += 0.0254) prof.push([pr, gy - 0.0011], [pr - 0.0008, gy], [pr, gy + 0.0011]);
      prof.push([pr, H - 0.006], [0, H - 0.006]);
      s.add(p.chrome, lathe(s, prof, 8), [x, 0, z]);
    } else if (s.detail === 0) {
      s.add(p.chrome, new THREE.BoxGeometry(pr * 1.8, H, pr * 1.8), [x, H / 2, z]);
      continue;
    } else s.add(p.chrome, new THREE.CylinderGeometry(pr, pr, H - 0.034, s.fine(2) ? 12 : 6, 1, true), [x, 0.028 + (H - 0.034) / 2, z]);
    // Levelling foot, and the cap on top.
    if (s.fine(2) && !lite) {
      s.add(p.plastic, lathe(s, [[0, 0], [0.016, 0], [0.017, 0.004], [0.014, 0.012], [0, 0.012]], 16), [x, 0, z]);
      s.add(p.chrome, cyl(s, 0.006, 0.02, "y", 8), [x, 0.02, z]);
      s.add(p.plastic, lathe(s, [[0, 0], [pr + 0.0006, 0], [pr + 0.0006, 0.008], [pr * 0.7, 0.012], [0, 0.013]], 16), [x, H - 0.013, z]);
    } else {
      s.add(p.plastic, new THREE.CylinderGeometry(0.014, 0.016, 0.028, 6), [x, 0.014, z]);
      s.add(p.plastic, new THREE.CylinderGeometry(pr * 0.8, pr + 0.0006, 0.012, 6), [x, H - 0.006, z]);
    }
  }
  const wr = 0.0022;
  for (const y of tiers) {
    if (s.detail === 0) {
      s.add(p.chrome, new THREE.BoxGeometry(L - 0.01, 0.01, D - 0.01), [L / 2, y - 0.005, -D / 2]);
      continue;
    }
    // Sleeves and collars at the posts.
    for (const [x, z] of posts) {
      if (!s.fine(2) || lite) {
        s.add(p.chrome, new THREE.CylinderGeometry(0.0205, 0.0186, 0.045, lite && s.fine(2) ? 10 : 6, 1, true), [x, y - 0.0235, z]);
        continue;
      }
      s.add(p.chrome, lathe(s, [[0.0142, y - 0.046], [0.0186, y - 0.046], [0.0208, y - 0.004], [0.0198, y - 0.001], [0.0142, y - 0.001]], cs), [x, 0, z]);
      s.add(p.poly, lathe(s, [[0.0128, y - 0.054], [0.0139, y - 0.054], [0.0151, y - 0.045], [0.0128, y - 0.045]], cs - 4), [x, 0, z]);
      if (s.fine(2) && !lite) s.add(p.enamel, new THREE.BoxGeometry(0.0006, 0.007, 0.0012), [x + 0.0142, y - 0.0496, z]);
    }
    const xa = c + 0.019;
    const xb = L - c - 0.019;
    const za = -c - 0.019;
    const zb = -D + c + 0.019;
    // Edge trusses along the front and back: top and bottom chords, zig-zag between.
    for (const z of [-0.006, -D + 0.006]) {
      wire(s, p.chrome, [xa, y - 0.0028, z], [xb, y - 0.0028, z], 0.0028);
      wire(s, p.chrome, [xa, y - 0.03, z], [xb, y - 0.03, z], 0.0024);
      const n = Math.max(2, Math.round((xb - xa) / (s.fine(2) ? (lite ? 0.08 : 0.05) : 0.3)));
      for (let i = 0; i < n; i++) {
        const x0 = xa + ((xb - xa) * i) / n;
        const x1 = xa + ((xb - xa) * (i + 1)) / n;
        wire(s, p.chrome, i % 2 ? [x0, y - 0.03, z] : [x0, y - 0.003, z], i % 2 ? [x1, y - 0.003, z] : [x1, y - 0.03, z], 0.0019);
      }
    }
    // Side frames: top and lower wires with three struts.
    for (const x of [0.006, L - 0.006]) {
      wire(s, p.chrome, [x, y - 0.0028, za], [x, y - 0.0028, zb], 0.0028);
      wire(s, p.chrome, [x, y - 0.022, za], [x, y - 0.022, zb], 0.0022);
      if (s.fine(2)) for (const f of [0.25, 0.5, 0.75]) wire(s, p.chrome, [x, y - 0.022, za + (zb - za) * f], [x, y - 0.003, za + (zb - za) * f], 0.0019);
    }
    // Longitudinal supports under the mat.
    if (s.fine(2)) for (const f of [0.2, 0.5, 0.8]) wire(s, p.chrome, [0.008, y - 0.0068, -D * f], [L - 0.008, y - 0.0068, -D * f], 0.0025);
    // The mat: wires front to back an inch apart, close up; an alpha-tested sheet from the aisle.
    if (s.fine(2)) {
      for (let x = 0.022; x < L - 0.02; x += 0.0254) wire(s, p.chrome, [x, y - wr, -0.007], [x, y - wr, -D + 0.007], wr, lite ? 3 : 5);
    } else {
      const g = new THREE.PlaneGeometry(L - 0.012, D - 0.012);
      g.rotateX(-Math.PI / 2);
      const uv = g.attributes.uv as THREE.BufferAttribute;
      for (let i = 0; i < uv.count; i++) uv.setXY(i, uv.getX(i) * ((L - 0.012) / (8 * 0.0254)), uv.getY(i) * ((D - 0.012) / 0.3));
      s.add(wireMatMaterial(), g, [L / 2, y - wr, -D / 2]);
    }
  }
  if (o.stock?.length && s.detail > 0) shelfStock(o.stock, L, D, tiers, o.fill ?? 0.7, o.seed ?? 7, o.cartons ?? true);
}

/**
 * Dry and cold stock on a shelving unit: kraft cases with printed panels, white tubs with coloured lids, #10 cans,
 * a lidded steam-table pan or two, placed across each shelf to `fill`.
 */
function shelfStock(shops: Shop[], L: number, D: number, tiers: number[], fill: number, seed: number, cartons: boolean): void {
  const p = palette();
  const gm = goodsMaterial();
  const r = rng(seed);
  let item = 0;
  let s = shops[0]!;
  const cs = s.fine(2) ? 28 : s.fine(1) ? 8 : 6;
  const wallUV = (g: THREE.BufferGeometry, t: number) => intoTile(g, tileUV(t));
  const caseBox = (w: number, h: number, d: number, at: V3) => {
    const g = new THREE.BoxGeometry(w, h, d);
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let f = 0; f < 6; f++) {
      const t = tileUV(f === 4 ? T.kraftLabel : T.kraft);
      for (let i = f * 4; i < f * 4 + 4; i++) uv.setXY(i, t[0] + (t[2] - t[0]) * uv.getX(i), t[1] + (t[3] - t[1]) * uv.getY(i));
    }
    s.add(gm, g, at);
    // The taped seam over the top flaps.
    if (s.fine(2)) s.add(gm, intoTile(new THREE.BoxGeometry(0.05, 0.0008, d + 0.001), tileUV(T.white)), [at[0], at[1] + h / 2 + 0.0004, at[2]]);
  };
  for (let ti = 0; ti < tiers.length; ti++) {
    const y = tiers[ti]!;
    const clear = ti + 1 < tiers.length ? tiers[ti + 1]! - y - 0.06 : 0.36;
    let x = 0.045;
    while (x < L - 0.08) {
      s = shops[item++ % shops.length]!;
      const kind = cartons ? r() : 0.45 + r() * 0.55;
      if (r() > fill) {
        x += 0.12;
        continue;
      }
      if (kind < 0.45) {
        const w = 0.25 + r() * 0.15;
        const h = Math.min(clear, 0.18 + r() * 0.16);
        const d = Math.min(D - 0.06, 0.3 + r() * 0.15);
        if (x + w > L - 0.04) break;
        caseBox(w, h, d, [x + w / 2, y + h / 2 + 0.0005, -0.03 - d / 2]);
        if (r() < 0.4 && h * 2 < clear) caseBox(w * 0.96, h, d * 0.96, [x + w / 2, y + h * 1.5 + 0.001, -0.03 - d / 2]);
        x += w + 0.02;
      } else if (kind < 0.7) {
        // White HDPE tubs, two deep.
        const rr = 0.09 + r() * 0.03;
        const h = Math.min(clear, 0.16 + r() * 0.1);
        if (x + 2 * rr > L - 0.04) break;
        const lid = [T.lidRed, T.lidBlue, T.lidGreen][Math.floor(r() * 3)]!;
        for (const z of [-0.03 - rr, -0.03 - 3 * rr - 0.01]) {
          if (-z + rr > D - 0.02) continue;
          if (!s.fine(2)) {
            s.add(p.hdpe, new THREE.CylinderGeometry(rr, rr * 0.92, h, cs, 1, true), [x + rr, y + h / 2 + 0.0005, z]);
            s.add(gm, intoTile(new THREE.CylinderGeometry(rr + 0.004, rr + 0.004, 0.012, cs), tileUV(lid)), [x + rr, y + h + 0.006, z]);
            continue;
          }
          // A moulded tub: tapered wall, a stacking ring low down, the lip rolled out at the top; its lid's skirt
          // snapped over the lip, the top recessed in a ring.
          const b0 = rr * 0.9;
          s.add(p.hdpe, lathe(s, [[0, 0.0005], [b0 - 0.006, 0.0005], [b0, 0.004], [b0 + 0.0012, 0.016], [b0 + 0.0028, 0.019], [b0 + 0.003, 0.024], [rr - 0.002, h - 0.016], [rr + 0.002, h - 0.012], [rr + 0.003, h - 0.008]], cs), [x + rr, y, z]);
          s.add(gm, wallUV(lathe(s, [[rr + 0.0045, h - 0.0115], [rr + 0.005, h - 0.004], [rr + 0.004, h + 0.002], [rr - 0.004, h + 0.0028], [rr - 0.008, h + 0.0012], [rr - 0.014, h + 0.0012], [rr - 0.018, h + 0.003], [0, h + 0.003]], cs), lid), [x + rr, y, z]);
        }
        x += 2 * rr + 0.02;
      } else if (kind < 0.88) {
        // A row of #10 cans.
        const rr = 0.078;
        const h = 0.178;
        if (h > clear || x + 2 * rr > L - 0.04) {
          x += 0.1;
          continue;
        }
        const label = r() < 0.5 ? T.canRed : T.canGreen;
        for (let k = 0; k < 2; k++) {
          const z = -0.03 - rr - k * (2 * rr + 0.006);
          if (-z + rr > D - 0.02) continue;
          if (!s.fine(2)) {
            s.add(gm, intoTile(new THREE.CylinderGeometry(rr, rr, h, cs, 1, true), tileUV(label)), [x + rr, y + h / 2, z]);
            s.add(p.alu, new THREE.CircleGeometry(rr, cs).rotateX(-Math.PI / 2), [x + rr, y + h, z]);
            continue;
          }
          // A #10 can: the labelled body between beaded shoulders, double-seamed rims top and bottom, the ends
          // dished inside the seams.
          s.add(gm, intoTile(new THREE.CylinderGeometry(rr, rr, h - 0.03, cs, 1, true), tileUV(label)), [x + rr, y + h / 2, z]);
          const seam: [number, number][] = [[rr - 0.0035, 0.0005], [rr + 0.0012, 0.0005], [rr + 0.0016, 0.002], [rr + 0.0016, 0.0045], [rr + 0.0006, 0.006], [rr, 0.0085], [rr - 0.0003, 0.011], [rr, 0.0135], [rr - 0.0003, 0.015]];
          s.add(p.alu, lathe(s, seam, cs), [x + rr, y, z]);
          s.add(p.alu, lathe(s, seam.map(([a, b]) => [a, h - b] as [number, number]).reverse(), cs), [x + rr, y, z]);
          s.add(p.alu, lathe(s, [[rr - 0.0035, h - 0.0005], [rr - 0.0055, h - 0.003], [0, h - 0.0026]], cs), [x + rr, y, z]);
        }
        x += 2 * rr + 0.012;
      } else {
        const pw = 0.325;
        if (x + pw > L - 0.04 || clear < 0.12 || D < 0.58) {
          x += 0.1;
          continue;
        }
        s.at([x + pw / 2, y + 0.0005, -0.03 - 0.265], [0, Math.PI / 2, 0], () => hotelPan(s, 0.53, 0.325, 0.1, T.chili, 0.8, true, seed + ti));
        x += pw + 0.04;
      }
    }
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Casters and the speed rack                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A 5 in plate-mount swivel caster in a joint whose origin is on the floor under its swivel axis: the top plate
 * (underside at `yPlate`) on a double ball raceway round the kingpin, a pressed-steel yoke trailing the wheel
 * behind the axis so it castors, a polyurethane tread on a hub, the axle bolt and nut; a toe brake when asked.
 */
function swivelCaster(s: Shop, yPlate: number, r = 0.0635, brake = false, trail = 0.036): void {
  const p = palette();
  const ww = 0.032;
  if (s.detail === 0) {
    s.add(p.ssDull, new THREE.BoxGeometry(0.05, yPlate - r, 0.04), [0, (yPlate + r) / 2, -trail / 2]);
    s.add(p.rubber, new THREE.CylinderGeometry(r, r, ww, 6).rotateZ(Math.PI / 2), [0, r, -trail]);
    return;
  }
  if (!s.fine(2)) {
    const ybr = yPlate - 0.02;
    s.add(p.ssDull, new THREE.BoxGeometry(0.1, 0.006, 0.085), [0, yPlate - 0.003, 0]);
    s.add(p.ssDull, new THREE.BoxGeometry(ww + 0.02, 0.014, 0.07), [0, ybr + 0.004, -trail * 0.35]);
    for (const sx of [-1, 1]) s.add(p.ssDull, new THREE.BoxGeometry(0.004, ybr - r, 0.03), [sx * (ww / 2 + 0.006), (ybr + r) / 2, -trail]);
    s.add(p.rubber, new THREE.CylinderGeometry(r, r, ww, 10).rotateZ(Math.PI / 2), [0, r, -trail]);
    return;
  }
  s.add(p.ssDull, bx(s, 0.1, 0.006, 0.085, 0.002), [0, yPlate - 0.003, 0]);
  if (s.fine(2))
    for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) {
      s.add(p.ssBright, cyl(s, 0.0055, 0.004, "y", 6), [sx * 0.038, yPlate - 0.008, sz * 0.03]);
    }
  // The raceway: two pressed cups about the kingpin.
  s.add(p.ssDull, lathe(s, [[0, -0.008], [0.037, -0.008], [0.04, -0.004], [0.036, 0], [0, 0]], s.fine(2) ? 16 : 8), [0, yPlate - 0.006, 0]);
  s.add(p.ssBright, cyl(s, 0.009, 0.01, "y", 6), [0, yPlate - 0.017, 0]);
  // The yoke: a bridge under the raceway and two legs down to the axle behind it.
  const ya = r;
  const yb = yPlate - 0.02;
  s.add(p.ssDull, new THREE.BoxGeometry(ww + 0.02, 0.005, 0.07), [0, yb, -trail * 0.35]);
  // Each leg a pressed plate, broad under the bridge, tapering down and back to a round boss about the axle.
  const leg = roundedPoly([
    [0.022, yb + 0.0025, 0.003],
    [-0.046, yb + 0.0025, 0.003],
    [-trail - 0.016, ya + 0.006, 0.014],
    [-trail - 0.004, ya - 0.016, 0.012],
    [-trail + 0.016, ya - 0.002, 0.012],
    [0.012, yb - 0.026, 0.012],
  ]);
  const lg = new THREE.ExtrudeGeometry(leg, { depth: 0.004, bevelEnabled: false, curveSegments: 4 });
  lg.rotateY(-Math.PI / 2);
  for (const sx of [-1, 1]) s.add(p.ssDull, lg.clone(), [sx * (ww / 2 + 0.006) + 0.002, 0, 0]);
  lg.dispose();
  // Wheel: tread on a hub, the axle bolt through the yoke.
  const tread = lathe(s, [
    [r * 0.62, -ww / 2],
    [r * 0.93, -ww / 2],
    [r, -ww * 0.3],
    [r, ww * 0.3],
    [r * 0.93, ww / 2],
    [r * 0.62, ww / 2],
  ], s.fine(2) ? 24 : 10);
  tread.rotateZ(Math.PI / 2);
  s.add(p.rubber, tread, [0, ya, -trail]);
  s.add(p.poly, cyl(s, r * 0.64, ww * 0.9, "x", s.fine(2) ? 20 : 8), [0, ya, -trail]);
  s.add(p.ssBright, cyl(s, 0.005, ww + 0.022, "x", 8), [0, ya, -trail]);
  if (s.fine(2)) for (const sx of [-1, 1]) s.add(p.ssBright, cyl(s, 0.008, 0.005, "x", 6), [sx * (ww / 2 + 0.011), ya, -trail]);
  if (brake && s.fine(1)) {
    s.add(p.plastic, bx(s, 0.03, 0.006, 0.04, 0.002), [0, ya + r * 0.55, -trail - r * 0.95], [0.5, 0, 0]);
    s.add(p.ssDull, bx(s, ww + 0.004, 0.004, 0.03, 0.0015), [0, ya + r + 0.004, -trail - r * 0.5], [0.25, 0, 0]);
  }
}

/**
 * The 20 slide heights of a 20-tier speed rack `h` tall, bottom first: the height of the rim of a sheet pan resting
 * on each slide, fixture frame y (the slide's tread runs 4.5 mm under it, where the rolled bead bears).
 */
export function speedRackSlides(h = 1.753, tiers = 20): number[] {
  const pitch = 0.076;
  const top = h - 0.09;
  return Array.from({ length: tiers }, (_, i) => top - (tiers - 1 - i) * pitch);
}

/** How far a slide's tread runs under its pan-rim height: the depth of the sheet pan's bead below its rim's top. */
const BEAD_DROP = 0.0045;

/**
 * The aluminium end-load speed rack: two welded side frames, each a pair of formed channel uprights (their open
 * side out, the slides welded across their webs) capped by a tube rail top and bottom, tied front and back by cross
 * tubes, on four 5 in swivel casters (the front pair braked). The 20 slides are formed runners: a web welded to the
 * uprights, a narrow tread turned in on top and a stiffening lip turned down at its edge, at 3 in. Sheet pans slide
 * in from the end and ride the treads by their rolled beads; each loaded tier carries a pan from the bake shop.
 */
const speedRack: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const yb = 0.147; // caster plates' underside: the frame's base
  const tx = 0.0222; // upright (and rail) width across the rack
  const tz = 0.0318; // upright depth along it
  const th = 0.0024; // the channel's gauge
  const yBase = yb + 0.006;
  const rail = 0.0254;
  const y0 = yBase + rail; // uprights from the base rail
  const y1 = h - rail; // to the top rail
  // Uprights: formed channel, open to the outside.
  const upright = (left: boolean, back: boolean) => {
    const R = s.fine(2) ? 0.003 : 0;
    const ri = s.fine(2) ? 0.0008 : 0;
    const zz0 = back ? d - tz : 0;
    const X = (x: number) => (left ? x : w - x);
    const sh = roundedPoly([
      [X(0), zz0, 0],
      [X(tx), zz0, R],
      [X(tx), zz0 + tz, R],
      [X(0), zz0 + tz, 0],
      [X(0), zz0 + tz - th, 0],
      [X(tx - th), zz0 + tz - th, ri],
      [X(tx - th), zz0 + th, ri],
      [X(0), zz0 + th, 0],
    ]);
    if (s.detail === 0) {
      s.add(p.alu, new THREE.BoxGeometry(tx, y1 - y0, tz), [X(tx / 2), (y0 + y1) / 2, -zz0 - tz / 2]);
      return;
    }
    const g = new THREE.ExtrudeGeometry(sh, { depth: y1 - y0, bevelEnabled: false, curveSegments: s.seg(3, 1) });
    g.rotateX(-Math.PI / 2);
    s.add(p.alu, g, [0, y0, 0]);
  };
  for (const left of [true, false]) for (const back of [false, true]) upright(left, back);
  // Rails: a tube along each side frame's top and foot, cross tubes front and back between them.
  for (const [y, hh] of [[h - rail / 2, rail], [yBase + rail / 2, rail]] as const) {
    for (const x of [tx / 2, w - tx / 2]) s.add(p.alu, bx(s, tx, hh, d, 0.002), [x, y, -d / 2]);
    for (const z of [-0.0115, -d + 0.0115]) s.add(p.alu, bx(s, w - 2 * tx - 0.0008, 0.019, 0.019, 0.0018), [w / 2, y, z]);
  }
  // Caster plates welded under the frame's corners.
  for (const x of [0.065, w - 0.065])
    for (const z of [-0.075, -d + 0.075]) s.add(p.alu, bx(s, 0.13, 0.006, 0.15, 0.0015), [x, yb + 0.003, z]);
  // The slides.
  const ys = speedRackSlides(h);
  {
    const tw = 0.0096; // tread
    const st = 0.0016;
    const lip = 0.005;
    const web = 0.019;
    const len = d - 0.001;
    const prof = (sx: number) => {
      const b = s.fine(2) ? 1 : 0;
      return roundedPoly([
        [sx * 0, 0, 0.0016 * b],
        [sx * tw, 0, 0.0012 * b],
        [sx * tw, -lip, 0],
        [sx * (tw - st), -lip, 0],
        [sx * (tw - st), -st, 0.0004 * b],
        [sx * st, -st, 0.0004 * b],
        [sx * st, -web, 0],
        [sx * 0, -web, 0],
      ]);
    };
    const runner = [1, -1].map((sx) => {
      if (s.detail === 0) return null;
      const g = new THREE.ExtrudeGeometry(prof(sx), { depth: len, bevelEnabled: false, curveSegments: s.seg(3, 1) });
      g.translate(0, 0, -len - 0.0005);
      return g;
    });
    for (const y of ys)
      for (const [k, sx] of [[0, 1], [1, -1]] as const) {
        const x0 = sx > 0 ? tx : w - tx;
        const yt = y - BEAD_DROP;
        const g = runner[k];
        if (g) s.add(p.alu, g.clone(), [x0, yt, 0]);
        else s.add(p.alu, new THREE.BoxGeometry(tw, 0.0016, d - 0.002), [x0 + (sx * tw) / 2, yt - 0.0008, -d / 2]);
        // Weld beads where the web meets each upright.
        if (s.fine(2)) for (const z of [-tz - 0.0012, -d + tz + 0.0012]) s.add(p.alu, new THREE.BoxGeometry(0.0026, web - 0.002, 0.0024), [x0 + sx * 0.0024, yt - web / 2, z]);
      }
    runner.forEach((g) => g?.dispose());
  }
  // Casters, each a swivel joint.
  const casters: [number, number][] = [[0.065, -0.075], [w - 0.065, -0.075], [0.065, -d + 0.075], [w - 0.065, -d + 0.075]];
  casters.forEach(([x, z], i) => {
    const j = s.joint(`caster_swivel_${i}`, [x, 0, z]);
    swivelCaster(j, yb, 0.0635, i < 2);
  });
  // Loads: a sheet pan on each tier from the bake shop, rolls, cookies, roast chicken, dough proofing, a sheet cake.
  ys.forEach((y, i) => {
    const j = s.joint(`pan_${i}`, [w / 2, y, -d / 2]);
    sheetPan(j, 0.457, 0.66);
    bakeLoad(j, i % 5, 0.457, 0.66, 41 + i);
  });
  return (b) => {
    const show = (n: number) => ys.forEach((_, i) => {
      const g = b.joints.get(`pan_${i}`);
      if (g) g.visible = i < n;
    });
    show(0);
    return {
      set: (name, v) => {
        if (name === "pans") show(Math.round(Math.max(0, Math.min(20, v))));
        else if (name === "swivel") for (let i = 0; i < 4; i++) turnJoint(b.joints.get(`caster_swivel_${i}`), "y", v);
      },
    };
  };
};

/** The sheet pan's section, [inset from the outline, y] from the outer floor up over the bead and down to the inner. */
const SHEET_PROFILE: readonly (readonly [number, number])[] = [
  [0.006, -0.025],
  [0.003, -0.0236],
  [0.0012, -0.005],
  [-0.0002, -0.0046],
  [-0.0014, -0.0037],
  [-0.0024, -0.0023],
  [-0.0026, -0.0012],
  [-0.0019, 0.0001],
  [-0.0006, 0.0004],
  [0.0008, -0.0002],
  [0.0016, -0.0012],
  [0.0026, -0.0022],
  [0.0036, -0.0236],
  [0.0066, -0.0246],
];
const SHEET_FAR: readonly (readonly [number, number])[] = [
  [0.006, -0.025],
  [0.0012, -0.005],
  [-0.0026, -0.0018],
  [0.0012, 0.0004],
  [0.0036, -0.0236],
  [0.0066, -0.0246],
];
/** The inner floor of a sheet pan, below its rim. */
const SHEET_FLOOR = -0.0246;

/**
 * An 18 × 26 in aluminium sheet pan, its rolled bead's top at the current origin's height (the bead's underside
 * 4.5 mm below, where it rides a rack's slides), its floor 25 mm below: drawn in one piece, the wall nearly upright.
 */
function sheetPan(s: Shop, w: number, d: number): void {
  const p = palette();
  if (s.detail === 0) {
    s.add(p.alu, new THREE.BoxGeometry(w - 0.004, 0.022, d - 0.004), [0, -0.0136, 0]);
    return;
  }
  const prof = s.fine(2) ? SHEET_PROFILE : SHEET_FAR;
  const R = 0.012;
  s.add(p.alu, loft(s, prof.map(([i, y]) => ({ w: w - 2 * i, d: d - 2 * i, r: R - i, y })), -1, 1, s.fine(2) ? 4 : 1));
}

/**
 * What a sheet pan carries off the bake shop's ovens: 0 dinner rolls, 1 cookies, 2 roast chicken thighs, 3 dough
 * balls proofing, 4 a sheet cake. Close up each piece is shaped and baked; from across the room the pan is one
 * picture of its load.
 */
function bakeLoad(s: Shop, kind: number, w: number, d: number, seed: number): void {
  if (s.detail === 0) return;
  const gm = goodsMaterial();
  const fy = SHEET_FLOOR + 0.0008;
  const iw = w - 0.016;
  const id = d - 0.016;
  if (!s.fine(2)) {
    const t = [T.panRolls, T.panCookies, T.panChicken, T.panDough, T.panCake][kind]!;
    // The picture lies just under the rim, where the tops of the pieces are.
    const hgt = [0.021, 0.011, 0.021, 0.021, 0.02][kind]!;
    const g = new THREE.PlaneGeometry(iw, id).rotateX(-Math.PI / 2);
    s.add(gm, intoTile(g, tileUV(t)), [0, fy + hgt, 0]);
    return;
  }
  const rnd = rng(seed);
  // Parchment under everything but the cake.
  if (kind !== 4) s.add(gm, intoTile(new THREE.PlaneGeometry(iw, id).rotateX(-Math.PI / 2), tileUV(T.parchment)), [0, fy, 0]);
  /** A piece seen from above: a lathe dome `r` across and `hgt` tall, squashed and lumped, its UVs projected down. */
  const piece = (tile: number, r: number, hgt: number, at: [number, number], sx: number, sz: number, lump: number, prof: [number, number][]) => {
    const g = lathe(s, prof.map(([a, b]) => [a * r, b * hgt] as [number, number]), 12);
    const pos = g.attributes.position as THREE.BufferAttribute;
    const k = rnd() * 9;
    for (let i = 0; i < pos.count; i++) {
      const x = pos.getX(i);
      const y = pos.getY(i);
      const z = pos.getZ(i);
      const a = Math.atan2(z, x);
      const f = 1 + lump * (Math.sin(a * 3 + k) * 0.6 + Math.sin(a * 5 + k * 2) * 0.4);
      pos.setXYZ(i, x * sx * f, y * (1 + lump * 0.5 * Math.sin(a * 2 + k)), z * sz * f);
    }
    g.computeVertexNormals();
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < uv.count; i++) uv.setXY(i, 0.5 + pos.getX(i) / (2 * r * sx * 1.1), 0.5 - pos.getZ(i) / (2 * r * sz * 1.1));
    s.add(gm, intoTile(g, tileUV(tile)), [at[0], fy, at[1]], [0, rnd() * Math.PI * 2, 0]);
  };
  const grid = (nu: number, nv: number, fn: (x: number, z: number) => void) => {
    for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) fn(-iw / 2 + ((i + 0.5) * iw) / nu, -id / 2 + ((j + 0.5) * id) / nv);
  };
  const dome: [number, number][] = [[0, 0], [0.93, 0], [1, 0.14], [0.96, 0.42], [0.82, 0.68], [0.58, 0.89], [0.28, 0.99], [0, 1]];
  const disc: [number, number][] = [[0, 0], [0.97, 0], [1, 0.35], [0.9, 0.8], [0.55, 0.97], [0, 1]];
  switch (kind) {
    case 0:
      grid(4, 6, (x, z) => piece(T.crust, 0.04, 0.042, [x + (rnd() - 0.5) * 0.006, z + (rnd() - 0.5) * 0.006], 1, 1.04, 0.035, dome));
      break;
    case 1:
      grid(4, 6, (x, z) => piece(T.cookie, 0.036, 0.011, [x + (rnd() - 0.5) * 0.01, z + (rnd() - 0.5) * 0.01], 1, 1, 0.06, disc));
      break;
    case 2:
      grid(3, 5, (x, z) => piece(T.skin, 0.052, 0.034, [x + (rnd() - 0.5) * 0.012, z + (rnd() - 0.5) * 0.012], 1.18, 0.86, 0.12, dome));
      break;
    case 3:
      grid(4, 6, (x, z) => piece(T.dough, 0.034, 0.036, [x, z], 1, 1, 0.02, dome));
      break;
    default: {
      // A sheet cake: baked to the pan's walls, its frosting combed.
      const g = foodSurface(s, w - 0.012, d - 0.012, 0.009, 0.002, T.frosting, seed);
      s.add(gm, g, [0, SHEET_FLOOR + 0.028, 0]);
      break;
    }
  }
}

const shelving: Builder = (s, e) => {
  const stock = [0, 1, 2].map((i) => s.joint(`stock_${i}`, [0, 0, 0]));
  wireShelving(s, e.w, e.d, e.h, { stock, fill: 0.8, seed: 11 });
  return (b) => ({
    set: (name, v) => {
      if (name === "stock") for (let i = 0; i < 3; i++) {
        const g = b.joints.get(`stock_${i}`);
        if (g) g.visible = v > i / 3 + 0.01;
      }
    },
  });
};

const cooler8x10 = walkIn({ freezer: false, door0: 0.76, doorW: 0.914, units: 1, fans: 2 });
const freezer8x8 = walkIn({ freezer: true, door0: 0.76, doorW: 0.914, units: 1, fans: 2 });
register({
  "NSF-WALKIN-COOLER": cooler8x10,
  "NSF-WALKIN-COOLER-8X10": cooler8x10,
  "NSF-WALKIN-COOLER-14X16": walkIn({ freezer: false, door0: 1.6, doorW: 1.067, units: 2, fans: 2 }),
  "NSF-WALKIN-FREEZER": freezer8x8,
  "NSF-WALKIN-FREEZER-8X8": freezer8x8,
  "NSF-WALKIN-FREEZER-8X10": walkIn({ freezer: true, door0: 0.76, doorW: 0.914, units: 1, fans: 2 }),
  "NSF-WALKIN-FREEZER-10X14": walkIn({ freezer: true, door0: 0.99, doorW: 1.067, units: 1, fans: 3 }),
  // The same 8 x 10 boxes with empty shelving along the back, for meat and fish placed case by case.
  "NSF-WALKIN-COOLER-8X10-BACK": walkIn({ freezer: false, door0: 0.76, doorW: 0.914, units: 1, fans: 2, back: true }),
  "NSF-WALKIN-FREEZER-8X10-BACK": walkIn({ freezer: true, door0: 0.76, doorW: 0.914, units: 1, fans: 2, back: true }),
  "NSF-WIRE-SHELVING-48": shelving,
  "NSF-WIRE-SHELVING-48IN": shelving,
  "NSF-SPEEDRACK-20PAN": speedRack,
});

/* ------------------------------------------------------------------------------------------------------------ */
/* Walk-ins                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/** 4 in foamed-in-place polyurethane panels: R-28 in a cooler, R-32 in a freezer. */
const PANEL = 0.1016;
/** Panels come 46 in wide, cam-locked edge to edge. */
const MODULE = 1.168;
/** The front wall's outer face stands back of the box front so the door's hardware, handle foremost, stays in it. */
const WALL_FRONT = -0.087;

interface WalkInSpec {
  freezer: boolean;
  /** Clear door opening: its left edge and width (36 in or 42 in), 78 in high. */
  door0: number;
  doorW: number;
  /** Unit coolers on the back wall, and fans on each. */
  units: number;
  fans: number;
  /** Shelving along the back wall as well, under the unit cooler: BACK_SHELF, left empty for cases put on it one by
   *  one. The shelving down the sides stops short of it, one unit a side, leaving room to stand at it. */
  back?: boolean;
}

/** A walk-in's back-wall shelving (WalkInSpec.back), in the box's frame: two 42 in units side by side from `x0`, 18 in
 *  deep, their backs `gap` off the back wall's inside, on three shelves (the posts clear of the unit cooler's drain,
 *  which goes into the wall 1.29 m up). */
export const BACK_SHELF = { L: 1.067, D: 0.457, H: 1.27, tiers: [0.2, 0.67, 1.14] as readonly number[], x0: [PANEL + 0.012, PANEL + 0.012 + 1.067 + 0.01] as readonly number[], gap: 0.012 } as const;

/** A panel joint's groove: the width of the dark line it draws, its skins turned in either side. */
const GROOVE = 0.004;

interface WallOpts {
  /** Run the outer skin on past the wall's left or right end, over a meeting wall's foam (negative: stop short). */
  extL?: number;
  extR?: number;
  /** Draw the foam's end faces at the left and right ends (off where a meeting wall's skin covers them). */
  ends?: [boolean, boolean];
}

/**
 * A wall of panels between (x0, y0) and (x1, y1) in a frame whose +z is the wall's outer face (at z = 0) and whose
 * x runs along it: a foam core between an outer and an inner skin. At every panel joint (`seams`, along x) each skin's
 * edge turns in to make a shallow V-groove with the joint's dark sealant at its root, the line a walk-in's panels
 * read by from across a room; cam-lock plug buttons inside along the joints at the heights `plugs`.
 */
function panelWall(s: Shop, x0: number, x1: number, y0: number, y1: number, outer: THREE.Material, inner: THREE.Material, seams: number[], plugs: number[], o: WallOpts = {}): void {
  const p = palette();
  const t = PANEL;
  const sk = 0.0008;
  const h = y1 - y0;
  const [endL, endR] = o.ends ?? [true, true];
  const drop: ("px" | "nx" | "py" | "ny" | "pz" | "nz")[] = ["pz", "nz"];
  if (!endL) drop.push("nx");
  if (!endR) drop.push("px");
  // Where a meeting wall's skin covers an end, the foam stops clear of it.
  const f0 = x0 + (endL ? 0 : 0.001);
  const f1 = x1 - (endR ? 0 : 0.001);
  s.add(p.ceramic, openBox(f1 - f0, h, t - 0.003, drop), [(f0 + f1) / 2, y0 + h / 2, -t / 2]);
  // From across the room the skins run unbroken: a gap there would show through the wall.
  const joints = s.fine(1) ? seams.filter((x) => x > x0 + 0.02 && x < x1 - 0.02) : [];
  const g = GROOVE / 2;
  const cuts = [x0 - (o.extL ?? 0), ...joints, x1 + (o.extR ?? 0)];
  for (let i = 0; i + 1 < cuts.length; i++) {
    const a = cuts[i]! + (i > 0 ? g : 0);
    const b = cuts[i + 1]! - (i + 2 < cuts.length ? g : 0);
    s.add(outer, new THREE.BoxGeometry(b - a, h, sk), [(a + b) / 2, y0 + h / 2, -sk / 2]);
    const ai = Math.max(a, x0);
    const bi = Math.min(b, x1);
    if (bi > ai) s.add(inner, new THREE.BoxGeometry(bi - ai, h, sk), [(ai + bi) / 2, y0 + h / 2, -t + sk / 2]);
  }
  for (const x of joints) {
    if (s.fine(2))
      for (const sx of [-1, 1]) {
        s.add(outer, new THREE.BoxGeometry(sk, h - 0.002, 0.0022), [x + sx * (g + sk / 2), y0 + h / 2, -sk - 0.0011]);
        s.add(inner, new THREE.BoxGeometry(sk, h - 0.002, 0.0022), [x + sx * (g + sk / 2), y0 + h / 2, -t + sk + 0.0011]);
      }
    const dz = s.fine(2) ? 0.003 : 0.0016;
    s.add(p.rubber, new THREE.PlaneGeometry(GROOVE + 0.0012, h - 0.002), [x, y0 + h / 2, -dz]);
    s.add(p.rubber, new THREE.PlaneGeometry(GROOVE + 0.0012, h - 0.002), [x, y0 + h / 2, -t + dz], [0, Math.PI, 0]);
  }
  if (s.fine(2))
    for (const x of joints)
      for (const y of plugs) if (y > y0 + 0.05 && y < y1 - 0.05) camPlug(s, [x - 0.045, y, -t]);
}

/** A cam-lock access hole's snap-in plug button on an inner skin, facing -z (into the room). */
function camPlug(s: Shop, at: V3, rot: V3 = [0, 0, 0]): void {
  const g = lathe(s, [[0, 0.0004], [0.0115, 0.0004], [0.0105, 0.0016], [0.007, 0.0032], [0, 0.0036]], 10);
  g.rotateX(-Math.PI / 2);
  s.at(at, rot, () => s.add(palette().hdpe, g, [0, 0, 0]));
}

/**
 * A heavy cam-lift hinge's barrel, split along a two-start helix: the lower half fixed to the frame, the upper half
 * to the door. As the door swings out the upper half rides up the helix and lifts the door (12 mm at 90°); let go,
 * the door's weight runs it back down and shut. `lower` picks the half; `y` is the split's mid height.
 */
function camBarrel(s: Shop, mat: THREE.Material, r: number, y: number, len: number, lower: boolean, lift = 0.012): void {
  const n = s.fine(2) ? 14 : 8;
  const L = lift;
  const f = (phi: number) => y + L - ((2 * L) / Math.PI) * phi; // phi in [0, π)
  const pos: number[] = [];
  const idx: number[] = [];
  const yEnd = lower ? y - len / 2 : y + len / 2;
  const ch = 0.0008; // the chamfer either side of the split that shows it as a line
  // Columns: two helix sweeps of π each, the second a copy turned half round; each column is end, chamfer, split.
  for (let half = 0; half < 2; half++)
    for (let i = 0; i <= n; i++) {
      const phi = (i / n) * Math.PI * 0.999;
      const a = phi + half * Math.PI;
      const ys = f(phi);
      const c = Math.cos(a);
      const sn = -Math.sin(a);
      const yc = lower ? ys - 0.0012 : ys + 0.0012;
      pos.push(c * r, yEnd, sn * r, c * r, yc, sn * r, c * (r - ch), ys, sn * (r - ch));
    }
  const cols = 2 * (n + 1);
  for (let k = 0; k < cols; k++) {
    const k2 = (k + 1) % cols;
    for (let row = 0; row < 2; row++) {
      const a = k * 3 + row;
      const b = k2 * 3 + row;
      if (lower) idx.push(a, b, a + 1, b, b + 1, a + 1);
      else idx.push(a, a + 1, b, b, a + 1, b + 1);
    }
  }
  // The end cap, and a pin head on the top.
  const base = pos.length / 3;
  pos.push(0, yEnd, 0);
  for (let k = 0; k < cols; k++) {
    const k2 = (k + 1) % cols;
    if (lower) idx.push(base, k2 * 3, k * 3);
    else idx.push(base, k * 3, k2 * 3);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  s.add(mat, g);
  if (!lower) s.add(mat, lathe(s, [[0, 0], [r, 0], [r * 0.95, 0.002], [r * 0.6, 0.006], [0, 0.007]], 18), [0, yEnd, 0]);
}

/** Diamond tread plate: raised lozenges in a herringbone on a sheet `w` by `h`, facing +z at z = 0. */
function diamondPlate(s: Shop, w: number, h: number): void {
  const p = palette();
  s.add(p.alu, bx(s, w, h, 0.0016, 0.0006), [0, 0, -0.0008]);
  if (!s.fine(2)) return;
  const pitch = 0.034;
  const loz = new THREE.OctahedronGeometry(1, 0);
  for (let j = 0; j * pitch < h - 0.012; j++)
    for (let i = 0; i * pitch < w - 0.012; i++) {
      const x = -w / 2 + 0.012 + i * pitch + (j % 2 ? pitch / 2 : 0);
      const y = -h / 2 + 0.012 + j * pitch;
      if (x > w / 2 - 0.01) continue;
      const g = loz.clone();
      s.add(p.alu, g, [x, y, 0.0008], [0, 0, (i + j) % 2 ? Math.PI / 4 : -Math.PI / 4], [0.0024, 0.0095, 0.0012]);
    }
}

/** A thermometer's live readout: a seven-segment face (red unless `col` says) redrawn when the temperature moves. */
function liveReadout(col = "#ff3b2e"): { mat: THREE.MeshBasicMaterial; set: (c: number) => void } {
  const c = document.createElement("canvas");
  c.width = 256;
  c.height = 96;
  const x = c.getContext("2d")!;
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  const mat = new THREE.MeshBasicMaterial({ map: tex });
  mat.userData.uv = "own";
  let last = NaN;
  const set = (deg: number) => {
    const f = Math.round((deg * 1.8 + 32) * 10) / 10;
    if (f === last) return;
    last = f;
    x.fillStyle = "#06090b";
    x.fillRect(0, 0, 256, 96);
    x.globalAlpha = 0.08;
    x.fillStyle = col;
    x.font = "bold 70px monospace";
    x.textAlign = "right";
    x.textBaseline = "middle";
    x.fillText("88.8", 196, 52);
    x.globalAlpha = 1;
    x.fillText(f.toFixed(1), 196, 52);
    x.font = "bold 34px monospace";
    x.fillText("°F", 250, 34);
    tex.needsUpdate = true;
  };
  return { mat, set };
}

interface DoorRig {
  /** Hinge pin axis (x, z), and the closer's frame pivot and door bracket (door-local, from the pin). */
  hx: number;
  hz: number;
  ax: number;
  az: number;
  bx: number;
  bz: number;
  rod0: number;
}

/**
 * The walk-in's door and its frame. The frame: a stainless casing round the opening, standing 3 mm proud of the
 * panel with its edges eased, the anti-sweat heater wire under it (its junction box and pilot over the head), the
 * jambs and head lined in stainless through the wall, a lighter casing inside, a bevelled threshold across the foot
 * (heated in a freezer). The leaf: a stainless slab with rolled edges lapping the casing, the liner's plug filling the
 * opening behind it with a black breaker strip round its root, a magnetic gasket on three sides of the slab's back
 * sealing on the casing, a neoprene sweep on an aluminium retainer down to the threshold, diamond tread kick plates
 * outside and in; three cam-lift hinges; a cast latch with its key cylinder and a pull grip, its bolt into a roller
 * strike on the frame; the glow-in-the-dark safety release through the leaf inside; a hydraulic closer across the
 * head. Everything that swings is under `cam_lift_door`.
 */
function walkinDoor(s: Shop, k: WalkInSpec, floorY: number, doorH: number, skinOut: THREE.Material, skinIn: THREE.Material): DoorRig {
  const p = palette();
  const zf = WALL_FRONT;
  const zi = zf - PANEL;
  const lap = 0.045;
  const x0 = k.door0;
  const x1 = k.door0 + k.doorW;
  const dTop = floorY + doorH;
  const sl = x0 - lap;
  const sr = x1 + lap;
  const sy0 = floorY + 0.02;
  const sy1 = dTop + lap;
  const st = 0.0272;
  const zs = zf + 0.035; // the slab's outer face
  const zB = zs - st; // and its back, 7.8 mm off the wall: the casing and the gasket between
  const hx = sl - 0.004;
  const hz = zs + 0.013;
  const lift = 0.012;
  const band = 0.058;
  const cz = 0.003;
  // The frame: casing outside and in, jamb and head liners through the wall.
  const casing = (b: number, r: number) =>
    roundedPoly([
      [x0 - b, floorY, 0],
      [x0, floorY, 0],
      [x0, dTop, r * 0.4],
      [x1, dTop, r * 0.4],
      [x1, floorY, 0],
      [x1 + b, floorY, 0],
      [x1 + b, dTop + b, r],
      [x0 - b, dTop + b, r],
    ]);
  s.add(p.ss, bevelPlate(s, casing(band, 0.01), cz, 0.0012), [0, 0, zf + 0.0003]);
  if (s.detail > 0) s.add(p.ss, bevelPlate(s, casing(0.04, 0.006), 0.0025, 0.0008), [0, 0, zi - 0.0028]);
  for (const x of [x0 + 0.0006, x1 - 0.0006]) s.add(p.ss, new THREE.BoxGeometry(0.0012, doorH, PANEL), [x, floorY + doorH / 2, zf - PANEL / 2]);
  s.add(p.ss, new THREE.BoxGeometry(k.doorW, 0.0012, PANEL), [x0 + k.doorW / 2, dTop - 0.0006, zf - PANEL / 2]);
  // The threshold: a stainless plate bevelled both edges, across the foot of the opening and over the wall's depth.
  const thTop = floorY + (k.freezer ? 0.006 : 0.005);
  {
    const za = zi - 0.012;
    const zb = zf + 0.03;
    const y0 = floorY + 0.0004;
    const prof: [number, number][] = [[za, y0], [za + 0.014, thTop], [zb - 0.014, thTop], [zb, y0]];
    const len = k.doorW + 0.03;
    s.add(p.ss, sweepX(prof, len), [x0 + k.doorW / 2, 0, 0]);
    if (s.fine(1)) {
      s.add(p.ss, capX(prof, -1), [x0 + k.doorW / 2 - len / 2, 0, 0]);
      s.add(p.ss, capX(prof, 1), [x0 + k.doorW / 2 + len / 2, 0, 0]);
    }
    if (s.fine(2)) for (const f of [0.08, 0.5, 0.92]) screwLite(s, p.ssBright, [x0 - 0.015 + len * f, thTop, (za + zb) / 2], 0.0028, [-Math.PI / 2, 0, 0]);
  }
  // The roller strike on the frame beside the latch, the heater's junction box and pilot lamp over the head.
  const ly = floorY + 1.08;
  {
    const zb = zf + cz + 0.0003;
    const dz = zs + 0.024 - zb;
    s.add(p.ssBright, bx(s, 0.03, 0.12, dz, 0.004), [sr + 0.021, ly, zb + dz / 2]);
    if (s.detail > 0) s.add(p.chrome, cyl(s, 0.008, 0.05, "y", 16), [sr + 0.009, ly, zs + 0.012]);
    if (s.fine(2)) for (const y of [-0.045, 0.045]) screwLite(s, p.ssBright, [sr + 0.021, ly + y, zb + dz], 0.0025);
  }
  s.add(p.ssDull, bx(s, 0.1, 0.1, 0.045, 0.003), [sl + 0.16, sy1 + 0.13, zf + 0.0225]);
  s.add(p.ledAmber, cyl(s, 0.006, 0.004, "z", 12), [sl + 0.16, sy1 + 0.13, zf + 0.0468]);
  s.add(p.ssBright, cyl(s, 0.0095, 0.003, "z", 16), [sl + 0.16, sy1 + 0.13, zf + 0.0459]);
  if (s.fine(1)) s.at([sl + 0.33, sy1 + 0.13, zf + 0.0002], [0, 0, 0], () => dataPlate(s, 0.1, 0.07));

  // Hinges: the frame leaf, offset to the slab's face by its block, and the lower barrel half.
  const hys = [sy0 + 0.36, (sy0 + sy1) / 2 + 0.05, sy1 - 0.24];
  for (const y of hys) {
    s.add(p.chrome, bx(s, 0.05, 0.17, zs + 0.004 - zf, 0.005), [hx - 0.03, y, (zf + zs + 0.004) / 2]);
    if (s.detail === 0) s.add(p.chrome, new THREE.CylinderGeometry(0.016, 0.016, 0.17, 6), [hx, y, hz]);
    else s.at([hx, 0, hz], [0, 0, 0], () => camBarrel(s, p.chrome, 0.016, y, 0.17, true, lift));
    if (s.fine(2)) for (const dy of [-0.06, -0.02, 0.02, 0.06]) screwLite(s, p.ssBright, [hx - 0.036, y + dy, zs + 0.004], 0.003);
  }

  const j = s.joint("cam_lift_door", [hx, 0, hz]);
  const X = (x: number) => x - hx;
  const Z = (z: number) => z - hz;
  const sw = sr - sl;
  const cx = X(sl + sw / 2);
  const kick = 0.31;
  // The slab: stainless wrapped round its foam, the edges rolled.
  j.add(skinOut, bx(j, sw, sy1 - sy0, st, 0.006), [cx, (sy0 + sy1) / 2, Z(zs - st / 2)]);
  // The plug: the liner's pan filling the opening, its face just shy of the wall's inner face.
  const px0 = x0 + 0.007;
  const px1 = x1 - 0.007;
  const py1 = dTop - 0.007;
  const zP = zi + 0.012;
  const pcx = X((px0 + px1) / 2);
  const pw = px1 - px0;
  j.add(skinIn, bx(j, pw, py1 - sy0, zB + 0.002 - zP, 0.004), [pcx, (sy0 + py1) / 2, Z((zB + 0.002 + zP) / 2)]);
  if (j.detail > 0) {
    // The black ABS breaker round the plug's root, where the cold liner meets the warm slab.
    const bw = 0.009;
    const bz = Z(zB - 0.0006);
    j.add(p.plastic, bx(j, bw, py1 - sy0, 0.0022, 0.0008), [X(px0 - bw / 2), (sy0 + py1) / 2, bz]);
    j.add(p.plastic, bx(j, bw, py1 - sy0, 0.0022, 0.0008), [X(px1 + bw / 2), (sy0 + py1) / 2, bz]);
    j.add(p.plastic, bx(j, pw + 2 * bw, bw, 0.0022, 0.0008), [pcx, py1 + bw / 2, bz]);
    // The magnetic gasket: a grey PVC bulb over its magnet, up the hinge side, across the head, down the latch side.
    const gi = 0.022;
    const path = roundPath([[X(sl + gi), sy0 + 0.004], [X(sl + gi), sy1 - gi], [X(sr - gi), sy1 - gi], [X(sr - gi), sy0 + 0.004]], 0.014, j.fine(2) ? 4 : 1);
    const prof: [number, number][] = j.fine(2)
      ? [[-0.006, 0], [0.006, 0], [0.0075, 0.0012], [0.0058, 0.0021], [0.0082, 0.003], [0.0068, 0.004], [-0.0068, 0.004], [-0.0082, 0.003], [-0.0058, 0.0021], [-0.0075, 0.0012]]
      : [[-0.007, 0], [0.007, 0], [0.0075, 0.004], [-0.0075, 0.004]];
    j.add(p.gasket, sweepPath(path, prof, Z(zB) + 0.0003, -1));
    // The sweep: an aluminium retainer across the foot of the plug's face, its neoprene blade down to the threshold.
    j.add(p.ssDull, bx(j, pw - 0.012, 0.05, 0.005, 0.0015), [pcx, sy0 + 0.019, Z(zP - 0.0024)]);
    j.add(p.rubber, bx(j, pw - 0.016, sy0 - 0.003 - (thTop + 0.0008), 0.0026, 0.0008), [pcx, (sy0 - 0.003 + thTop + 0.0008) / 2, Z(zP - 0.0024)]);
    if (j.fine(2)) for (let i = 0; i < 6; i++) screwLite(j, p.ssBright, [pcx - pw / 2 + 0.04 + (i * (pw - 0.08)) / 5, sy0 + 0.03, Z(zP - 0.0049)], 0.0022, [0, Math.PI, 0]);
    // Kick plate inside.
    j.at([pcx, sy0 + 0.05 + 0.125, Z(zP - 0.0019)], [0, Math.PI, 0], () => diamondPlate(j, pw - 0.03, 0.25));
  }
  // Kick plate outside, set in from the slab's rolled edges, screwed at its corners.
  j.at([cx, sy0 + 0.008 + kick / 2, Z(zs + 0.0019)], [0, 0, 0], () => diamondPlate(j, sw - 0.016, kick));
  if (j.fine(2)) for (const x of [-1, 1]) for (const y of [-0.13, 0.13]) screwLite(j, p.ssBright, [cx + x * (sw / 2 - 0.022), sy0 + 0.008 + kick / 2 + y, Z(zs + 0.0019)], 0.0025);
  // Hinge straps on the door, and the upper barrel halves riding the cam.
  for (const y of hys) {
    j.add(p.chrome, bx(j, 0.15, 0.09, 0.009, 0.004), [X(sl + 0.07), y + 0.02, Z(zs + 0.0047)]);
    j.add(p.chrome, bx(j, 0.03, 0.09, 0.02, 0.006), [0.012, y + 0.02, -0.002]);
    if (j.detail > 0) camBarrel(j, p.chrome, 0.016, y, 0.17, false, lift);
    if (j.fine(2)) for (const x of [0.04, 0.08, 0.12]) screwLite(j, p.ssBright, [X(sl + x), y + 0.02, Z(zs + 0.0092)], 0.003);
  }
  // The latch: a die-cast body near the free edge with its key cylinder, the pull grip held off it on two arms
  // where the hand closes, the bolt out of the edge into the strike.
  const lx = sr - 0.06;
  const gx = lx - 0.075;
  const gz = zs + 0.047;
  j.add(p.chrome, bevelPlate(j, roundRect(0.058, 0.25, 0.022), 0.024, j.fine(2) ? 0.006 : 0.003, 10), [X(lx), ly, Z(zs + 0.0002)]);
  j.add(p.chrome, cyl(j, 0.012, 0.012, "z", 20), [X(lx), ly + 0.085, Z(zs + 0.029)]);
  if (j.fine(2)) {
    j.add(p.enamel, bx(j, 0.0022, 0.012, 0.002, 0.0006), [X(lx), ly + 0.085, Z(zs + 0.0352)]);
    j.add(p.chrome, torus(j, 0.009, 0.0026, "x", Math.PI * 2, 18, 8), [X(lx), ly - 0.112, Z(zs + 0.018)]);
  }
  j.add(p.chrome, bar(j, [[X(lx - 0.016), ly + 0.07, Z(zs + 0.02)], [X(gx), ly + 0.07, Z(gz)], [X(gx), ly - 0.07, Z(gz)], [X(lx - 0.016), ly - 0.07, Z(zs + 0.02)]], 0.011, 0.022, 14));
  j.add(p.ssBright, bx(j, 0.05, 0.03, 0.012, 0.003), [X(sr + 0.004), ly, Z(zs - 0.012)]);
  // Inside: the safety release, a glow-in-the-dark knob on a stainless rose that anyone shut in can push.
  j.add(p.ssBright, cyl(j, 0.03, 0.006, "z", 24), [X(lx), ly, Z(zP - 0.0028)]);
  j.add(glowMaterial(), lathe(j, [[0, 0], [0.011, 0], [0.011, 0.014], [0.022, 0.017], [0.026, 0.024], [0.024, 0.031], [0.016, 0.037], [0, 0.039]], j.fine(2) ? 28 : 12), [X(lx), ly, Z(zP - 0.0058)], [-Math.PI / 2, 0, 0]);
  if (j.fine(2)) for (const dy of [-0.09, 0.09]) screwLite(j, p.ssBright, [X(lx), ly + dy, Z(zs + 0.0242)], 0.0026);

  // The closer: body on a bracket on the head, its rod to a bracket on the door's top rail.
  const yc = sy1 + 0.035;
  const zc = zf + 0.03;
  const ax = hx + 0.48;
  const cbx = 0.13;
  const bz = zc - hz;
  s.add(p.ssDull, bx(s, 0.05, 0.07, zc - zf + 0.01, 0.003), [ax + 0.02, yc, (zf + zc) / 2]);
  j.add(p.ssDull, bx(j, 0.035, 0.08, 0.0045, 0.002), [cbx, yc - 0.03, Z(zs + 0.00225)]);
  j.add(p.ssDull, bx(j, 0.035, 0.02, zc - zs + 0.01, 0.002), [cbx, yc, (bz + Z(zs)) / 2]);
  const body = s.joint("closer_body", [ax, yc, zc]);
  body.add(p.charcoal, cyl(body, 0.019, 0.23, "x", 24), [-0.125, 0, 0]);
  body.add(p.ssBright, cyl(body, 0.021, 0.012, "x", 24), [-0.236, 0, 0]);
  body.add(p.ssBright, cyl(body, 0.008, 0.03, "y", 12), [0, 0, 0]);
  const rod0 = ax - hx - cbx;
  const rod = body.joint("closer_rod", [0, 0, 0]);
  rod.add(p.chrome, cyl(rod, 0.0055, 0.3, "x", 12), [-rod0 + 0.15, 0, 0]);
  rod.add(p.ssBright, bx(rod, 0.022, 0.022, 0.022, 0.005), [-rod0, 0, 0]);
  return { hx, hz, ax, az: zc, bx: cbx, bz, rod0 };
}

/**
 * A low-profile unit cooler, `fans` across, hung from the ceiling on threaded rod with its intake to the wall behind:
 * aluminium casing, the finned coil across the back, a fan panel with spun venturi rings, a fan and wire guard in
 * each, the sloped drain pan under it all and its drain. Local frame: x 0..Lh, y 0..Hh, z -Dh..0 (fan face at 0).
 * Returns its size.
 */
function unitCooler(s: Shop, u: number, fans: number, freezer: boolean, fin: THREE.MeshStandardMaterial): { L: number; H: number; D: number } {
  const p = palette();
  const fr = fans > 2 ? 0.18 : 0.165;
  const pitch = fr * 2 + 0.1;
  const L = fans * pitch + 0.24;
  const yP = 0.075; // the drain pan's depth below the casing
  const H = 2 * (fr + 0.036) + yP;
  const D = 0.5;
  // Casing: top, ends, and the fan panel with its rings, folded edges, the drain pan under.
  s.add(p.alu, bx(s, L, 0.0016, D, 0.0006), [L / 2, H - 0.0008, -D / 2]);
  s.add(p.alu, bx(s, L, 0.03, 0.0016, 0.0006), [L / 2, H - 0.015, -D + 0.0008]);
  for (const x of [0.0008, L - 0.0008]) s.add(p.alu, bx(s, 0.0016, H - yP, D, 0.0006), [x, yP + (H - yP) / 2, -D / 2]);
  if (s.detail > 0) {
    const shape = roundRect(L - 0.004, H - yP - 0.004, 0.003);
    for (let i = 0; i < fans; i++) shape.holes.push(circleHoleCW(0.12 + pitch * (i + 0.5) - L / 2, 0, fr + 0.012));
    const g = plate(s, shape, 0.0016);
    s.add(p.alu, g, [L / 2, yP + (H - yP) / 2, -0.0016]);
  } else s.add(p.alu, new THREE.BoxGeometry(L, H - yP, 0.002), [L / 2, yP + (H - yP) / 2, -0.001]);
  // The coil across the back, intake face open to the wall.
  const cl = L - 0.24;
  const ch = H - yP - 0.05;
  const cd = 0.13;
  s.at([0.12, yP + 0.02, -D + cd + 0.004], [0, 0, 0], () => finCoil(s, fin, cl, ch, cd, freezer ? 4 : 6, false));
  // Frost on the intake face: a white crust that thickens (the joint scales in z) once the fins are furred over.
  const ff = s.joint(`part_frost_${u}`, [0, 0, -D + 0.004]);
  ff.add(fin, bx(ff, cl + 0.01, ch + 0.01, 0.01, 0.003), [0.12 + cl / 2, yP + 0.02 + ch / 2, -0.005]);
  // The drain pan: a sloped tray with turned-up lips, its outlet at the right end.
  s.add(p.alu, bx(s, L - 0.004, 0.0016, D - 0.004, 0.0006), [L / 2, yP - 0.045, -D / 2], [0.0, 0, -0.012]);
  for (const z of [-0.002, -D + 0.002]) s.add(p.alu, bx(s, L - 0.004, 0.05, 0.0016, 0.0006), [L / 2, yP - 0.022, z]);
  for (const x of [0.002, L - 0.002]) s.add(p.alu, bx(s, 0.0016, 0.05, D - 0.004, 0.0006), [x, yP - 0.022, -D / 2]);
  s.add(p.ssDull, cyl(s, 0.017, 0.035, "y", 16), [L - 0.06, yP - 0.065, -D / 2]);
  // Fans: a spun venturi ring, the motor on its spider behind, the blades (a joint), the guard.
  for (let i = 0; i < fans; i++) {
    const cx = 0.12 + pitch * (i + 0.5);
    const cy = yP + (H - yP) / 2;
    if (s.detail === 0) {
      s.add(p.charcoal, new THREE.CircleGeometry(fr, 8), [cx, cy, 0.001]);
      continue;
    }
    s.add(p.alu, new THREE.LatheGeometry([[fr + 0.013, -0.05], [fr + 0.008, -0.03], [fr + 0.012, -0.004], [fr + 0.03, 0.0], [fr + 0.032, 0.002], [fr + 0.011, 0.002], [fr + 0.004, -0.03], [fr + 0.009, -0.05]].map(([a, b]) => new THREE.Vector2(a, b)), s.fine(2) ? 32 : 12), [cx, cy, 0], [Math.PI / 2, 0, 0]);
    if (s.fine(2)) {
      s.add(p.charcoal, cyl(s, 0.055, 0.09, "z", 24), [cx, cy, -0.13]);
      for (let k = 0; k < 3; k++) {
        const a = (k / 3) * Math.PI * 2 + Math.PI / 2;
        wire(s, p.ssDull, [cx + Math.cos(a) * 0.05, cy + Math.sin(a) * 0.05, -0.12], [cx + Math.cos(a) * (fr + 0.01), cy + Math.sin(a) * (fr + 0.01), -0.035], 0.004, 6);
      }
    }
    const fj = s.joint(`evap_fan_${u}_${i}`, [cx, cy, -0.03]);
    axialFan(fj, p.alu, p.ssDull, fr - 0.008, 5);
    s.at([cx, cy, 0.003], [0, 0, 0], () => fanGuard(s, p.chrome, fr + 0.02, 0.015));
  }
  // The end compartment's access cover screws, the nameplate, an electrical box with its conduit stub.
  if (s.fine(2)) for (const x of [0.0, L]) for (const y of [yP + 0.04, H - 0.04]) for (const z of [-0.06, -D + 0.06]) screwLite(s, p.ssBright, [x + (x > 0 ? 0.0016 : -0.0016), y, z], 0.0025, [0, x > 0 ? Math.PI / 2 : -Math.PI / 2, 0]);
  s.at([L + 0.0016, H * 0.6, -D * 0.7], [0, Math.PI / 2, 0], () => dataPlate(s, 0.07, 0.05));
  s.add(p.ssDull, bx(s, 0.05, 0.1, 0.1, 0.003), [L + 0.025, H * 0.55, -0.12]);
  return { L, H, D };
}

/** A hole for a shape traced clockwise (the kit's `circleHole` with its own segment count). */
function circleHoleCW(cx: number, cy: number, r: number): THREE.Path {
  return new THREE.Path().absarc(cx, cy, r, 0, Math.PI * 2, true);
}

/**
 * A top-mount air-cooled condensing unit on its base rails, standing at the current origin: x 0..0.95, y up, z
 * -0.55..0. The hermetic compressor and its receiver in the open bay, the condenser coil across the back with the
 * fan drawing through it and out past its guard, the copper between them, a sight glass and filter drier, the
 * service valves, and the electrical box.
 */
function condensingUnit(s: Shop, fin: THREE.Material): void {
  const p = palette();
  const W = 0.95;
  const D = 0.55;
  const H = 0.52;
  for (const z of [-0.04, -D + 0.04]) s.add(p.ssDull, bx(s, W, 0.05, 0.05, 0.002), [W / 2, 0.025, z]);
  s.add(p.ssDull, bx(s, W, 0.002, D, 0.001), [W / 2, 0.051, -D / 2]);
  // The condenser section on the right: a sheet-metal shroud, the coil on its back and its open side, the fan out front.
  const cx0 = 0.42;
  const cw = W - cx0;
  s.add(p.charcoal, bx(s, cw, 0.0016, D, 0.0006), [cx0 + cw / 2, H, -D / 2]);
  s.add(p.charcoal, bx(s, 0.0016, H - 0.052, D, 0.0006), [cx0, (H + 0.052) / 2, -D / 2]);
  if (s.detail > 0) {
    const shape = roundRect(cw - 0.004, H - 0.056, 0.004);
    shape.holes.push(circleHoleCW(0, 0, 0.2));
    s.add(p.charcoal, plate(s, shape, 0.0016), [cx0 + cw / 2, (H + 0.052) / 2, -0.0016]);
  }
  s.at([cx0 + 0.02, 0.07, -D + 0.05], [0, 0, 0], () => finCoil(s, fin, cw - 0.04, H - 0.09, 0.044, 12, false));
  const fj = s.joint("cond_fan", [cx0 + cw / 2, (H + 0.052) / 2, -0.04]);
  if (s.detail === 0) {
    s.add(p.charcoal, bx(s, cw, H - 0.052, 0.004), [cx0 + cw / 2, (H + 0.052) / 2, -0.002]);
    s.add(p.enamel, new THREE.CylinderGeometry(0.085, 0.085, 0.26, 8), [0.17, 0.18, -0.3]);
    return;
  }
  axialFan(fj, p.enamel, p.enamel, 0.185, 4);
  s.add(p.charcoal, cyl(s, 0.06, 0.1, "z", 24), [cx0 + cw / 2, (H + 0.052) / 2, -0.12]);
  s.at([cx0 + cw / 2, (H + 0.052) / 2, 0.002], [0, 0, 0], () => fanGuard(s, p.enamel, 0.21, 0.012, 0.0015));
  // The compressor bay: compressor, receiver, drier, sight glass, valves, the electrical box on the bulkhead.
  s.at([0.17, 0.052, -0.3], [0, 0, 0], () => compressor(s, 0.085, 0.26));
  s.add(p.enamel, lathe(s, [[0, 0], [0.04, 0.004], [0.05, 0.03], [0.05, 0.27], [0.04, 0.296], [0, 0.3]], 24), [0.34, 0.11, -0.46], [0, 0, Math.PI / 2]);
  pipe(s, p.copper, [[0.34, 0.11, -0.46], [0.38, 0.13, -0.46], [0.38, 0.13, -0.2], [0.33, 0.2, -0.1]], 0.0048, 0.025, 8);
  s.add(p.copper, cyl(s, 0.02, 0.12, "y", 18), [0.33, 0.26, -0.1]);
  s.add(p.glass, cyl(s, 0.011, 0.004, "z", 16), [0.33, 0.35, -0.083]);
  s.add(p.brass, cyl(s, 0.014, 0.03, "z", 16), [0.33, 0.35, -0.1]);
  pipe(s, p.copper, [[0.33, 0.32, -0.1], [0.33, 0.42, -0.1], [0.22, 0.42, -0.1], [0.22, 0.42, -0.02]], 0.0048, 0.02, 8);
  pipe(s, p.copper, [[0.087, 0.3, -0.22], [0.087, 0.42, -0.22], [0.1, 0.42, -0.08], [0.1, 0.42, -0.02]], 0.0079, 0.025, 10);
  for (const x of [0.1, 0.22]) s.add(p.brass, bx(s, 0.03, 0.035, 0.03, 0.004), [x, 0.42, -0.015]);
  pipe(s, p.copper, [[0.24, 0.3, -0.33], [0.3, 0.42, -0.33], [cx0 - 0.01, 0.42, -0.33]], 0.0064, 0.02, 8);
  s.add(p.charcoal, bx(s, 0.16, 0.2, 0.08, 0.004), [cx0 - 0.09, 0.4, -D + 0.05]);
  s.at([cx0 - 0.0008, 0.3, -0.13], [0, -Math.PI / 2, 0], () => dataPlate(s, 0.08, 0.06));
}

/** The vinyl coved base NSF asks for where wall meets floor: a 25 mm fillet and a 100 mm skirt, along local x. */
function wallBase(s: Shop, x0: number, x1: number, y: number): void {
  if (x1 - x0 < 0.02) return;
  const p = palette();
  const r = 0.025;
  s.add(p.gasket, new THREE.BoxGeometry(x1 - x0, 0.1, 0.002), [(x0 + x1) / 2, y + 0.05, -PANEL - 0.001]);
  if (!s.fine(1)) return;
  const sh = new THREE.Shape();
  sh.moveTo(0, 0);
  sh.lineTo(r, 0);
  sh.absarc(r, r, r, -Math.PI / 2, Math.PI, true);
  sh.lineTo(0, 0);
  const g = plate(s, sh, x1 - x0, 0, 6);
  g.rotateY(Math.PI / 2);
  s.add(p.gasket, g, [x0, y, -PANEL - 0.0018]);
}

/**
 * A modular walk-in of 4 in panels. Cooler: floorless on a coved screed. Freezer: an insulated floor with an
 * aluminium wear skin and a heated door. Front wall stainless, the rest of the skins embossed cream; the door on the
 * left of the front wall swinging out, its strip curtain inside; the unit coolers on the back wall, shelving down the
 * sides, the condensing unit on the roof above.
 */
function walkIn(k: WalkInSpec): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const T = PANEL;
    const zf = WALL_FRONT;
    const floorY = k.freezer ? T : 0;
    const doorH = 1.981;
    const top = h - T;
    const zb = -d + T;
    const zi = zf - T;
    const fin = frostMaterial();
    const plugs = [floorY + 0.35, floorY + 1.2, top - 0.3];
    const seamsFrom = (a: number, b: number) => {
      const out: number[] = [];
      for (let x = a + MODULE; x < b - 0.3; x += MODULE) out.push(x);
      return out;
    };
    const roof = s.joint("part_roof", [0, 0, 0]);
    const wallR = s.joint("part_wall_right", [0, 0, 0]);
    const frontR = s.joint("part_front_right", [0, 0, 0]);

    // Side and back walls. The side walls run the box's full depth; the front and back walls' outer skins wrap over
    // their foam ends, and a formed angle trims each outside corner.
    const sideLen = d + zf;
    const sk = 0.0008;
    s.at([0, 0, -d], [0, -Math.PI / 2, 0], () => {
      panelWall(s, 0, sideLen, 0, top, p.panel, p.panel, seamsFrom(0, sideLen), plugs, { extL: -sk, extR: -sk, ends: [false, false] });
      wallBase(s, T, sideLen - T, floorY);
    });
    wallR.at([w, 0, zf], [0, Math.PI / 2, 0], () => {
      panelWall(wallR, 0, sideLen, 0, top, p.panel, p.panel, seamsFrom(0, sideLen), plugs, { extL: -sk, extR: -sk, ends: [false, false] });
      wallBase(wallR, T, sideLen - T, floorY);
    });
    // The back wall's frame runs from x = w (its 0) to x = 0: its skin's run over the left wall's end is its right.
    s.at([w, 0, -d], [0, Math.PI, 0], () => {
      panelWall(s, T, w - T, 0, top, p.panel, p.panel, seamsFrom(T, w - T), plugs, { extR: T });
      wallBase(s, T, w - T, floorY);
    });
    wallR.at([w, 0, -d], [0, Math.PI, 0], () => wallR.add(p.panel, new THREE.BoxGeometry(T, top, sk), [T / 2, top / 2, -sk / 2]));
    // The front wall: the door section (jambs and head), the panel left of it, and right of the latch-side jamb.
    const xr = k.door0 + k.doorW + 0.34;
    const dTop = floorY + doorH;
    const frontSeams = [k.door0 - 0.16, xr];
    s.at([0, 0, zf], [0, 0, 0], () => {
      panelWall(s, T, k.door0, 0, top, p.ss, p.panel, frontSeams, plugs, { extL: T });
      panelWall(s, k.door0, k.door0 + k.doorW, dTop, top, p.ss, p.panel, [], plugs);
      if (floorY > 0) panelWall(s, k.door0, k.door0 + k.doorW, 0, floorY, p.ss, p.panel, [], []);
      panelWall(s, k.door0 + k.doorW, xr, 0, top, p.ss, p.panel, [], plugs);
      wallBase(s, T, k.door0 - 0.045, floorY);
      wallBase(s, k.door0 + k.doorW + 0.045, xr, floorY);
    });
    frontR.at([0, 0, zf], [0, 0, 0], () => {
      panelWall(frontR, xr, w - T, 0, top, p.ss, p.panel, seamsFrom(xr, w - T), plugs, { extR: T });
      wallBase(frontR, xr, w - T, floorY);
    });
    // Corner trims: a 2 in formed angle up each outside corner, stainless on the front, painted at the back; they
    // stop under the roof's flashing, which laps over their heads.
    {
      const trim = (sh: Shop, x: number, z: number, rot: number, mat: THREE.Material) =>
        sh.at([x, 0, z], [0, rot, 0], () => {
          const g = angleBar(sh, 0.05, 0.05, top - 0.005, 0.0012, 0.005, 0.0006);
          g.rotateX(Math.PI / 2);
          g.translate(0, top - 0.002, 0);
          sh.add(mat, g);
          if (sh.fine(2)) for (const y of [0.15, top / 2, top - 0.15]) for (const [a, b] of [[0.032, 0.0018], [-0.0018, -0.032]] as const) screwLite(sh, p.ssBright, a > 0 ? [a, y, b] : [a, y, b], 0.0024, a > 0 ? [0, 0, 0] : [0, -Math.PI / 2, 0]);
        });
      trim(s, 0, zf, 0, p.ss);
      trim(s, 0, -d, -Math.PI / 2, p.panel);
      trim(wallR, w, -d, Math.PI, p.panel);
      trim(wallR, w, zf, Math.PI / 2, p.ss);
    }
    // The ceiling: panels spanning side to side, cam-locked and grooved at their joints like the walls, cream both
    // faces; a formed flashing angle round the top edge covering the panels' edges and lapping over the walls' heads,
    // stainless on the front and painted elsewhere, the side runs lapped over the front and back ones.
    {
      const cz = (zf - d) / 2;
      const cd = zf + d;
      roof.add(p.ceramic, openBox(w, T - 0.003, cd, ["py", "ny"]), [w / 2, top + T / 2, cz]);
      const seams = roof.fine(1) ? seamsFrom(0, cd).map((v) => zf - v) : [];
      const cuts = [zf, ...seams, -d];
      const g = GROOVE / 2;
      for (let i = 0; i + 1 < cuts.length; i++) {
        const a = cuts[i]! - (i > 0 ? g : 0);
        const b = cuts[i + 1]! + (i + 2 < cuts.length ? g : 0);
        roof.add(p.panel, new THREE.BoxGeometry(w, sk, a - b), [w / 2, h - sk / 2, (a + b) / 2]);
        roof.add(p.panel, new THREE.BoxGeometry(w, sk, a - b), [w / 2, top + sk / 2, (a + b) / 2]);
      }
      for (const z of seams) {
        if (roof.fine(2))
          for (const sz of [-1, 1]) {
            roof.add(p.panel, new THREE.BoxGeometry(w - 0.002, 0.0022, sk), [w / 2, h - sk - 0.0011, z + sz * (g + sk / 2)]);
            roof.add(p.panel, new THREE.BoxGeometry(w - 0.002, 0.0022, sk), [w / 2, top + sk + 0.0011, z + sz * (g + sk / 2)]);
          }
        const dy = roof.fine(2) ? 0.003 : 0.0016;
        roof.add(p.rubber, new THREE.PlaneGeometry(w - 0.002, GROOVE + 0.0012).rotateX(-Math.PI / 2), [w / 2, h - dy, z]);
        roof.add(p.rubber, new THREE.PlaneGeometry(w - 0.002, GROOVE + 0.0012).rotateX(Math.PI / 2), [w / 2, top + dy, z]);
        if (roof.fine(2)) for (let x = T + 0.3; x < w - T - 0.1; x += 0.6) camPlug(roof, [x, top, z + 0.045], [-Math.PI / 2, 0, 0]);
      }
      const leg = T + 0.03;
      const ft = 0.0012;
      const gF = 0.0024;
      const gS = 0.004;
      const eF = gF + ft;
      const eS = eF + 0.001;
      // The front and back runs lap over the sides' ends, so the stainless runs unbroken round the front corners.
      roof.at([-eS, h, zf], [0, Math.PI / 2, 0], () => roof.add(p.ss, angleBar(roof, 0.05, leg, w + 2 * eS, ft, 0.004, gS)));
      roof.at([w + eS, h, -d], [0, -Math.PI / 2, 0], () => roof.add(p.panel, angleBar(roof, 0.05, leg, w + 2 * eS, ft, 0.004, gS)));
      roof.at([0, h, -d - eF], [0, 0, 0], () => roof.add(p.panel, angleBar(roof, 0.05, leg, cd + 2 * eF, ft, 0.004, gF)));
      roof.at([w, h, zf + eF], [0, Math.PI, 0], () => roof.add(p.panel, angleBar(roof, 0.05, leg, cd + 2 * eF, ft, 0.004, gF)));
    }
    if (k.freezer) {
      // The insulated floor: 4 in panel, its wear skin of aluminium.
      s.add(p.ceramic, openBox(w - 2 * T, T - 0.003, zi - zb, ["py", "ny"]), [w / 2, T / 2, (zi + zb) / 2]);
      s.add(p.alu, new THREE.BoxGeometry(w - 2 * T, 0.0012, zi - zb), [w / 2, T - 0.0006, (zi + zb) / 2]);
    }

    // The door and its frame.
    const rig = walkinDoor(s, k, floorY, doorH, p.ss, p.panel);

    // Beside the latch: thermometer, light switch with its pilot, the pressure-relief port above.
    const sideX = k.door0 + k.doorW + 0.2;
    const readout = s.fine(1) ? liveReadout() : null;
    s.at([sideX, floorY + 1.55, zf], [0, 0, 0], () => {
      s.add(p.plastic, bx(s, 0.11, 0.07, 0.022, 0.005), [0, 0, 0.011]);
      if (readout) s.add(readout.mat, decalFull(0.084, 0.034), [0, 0.006, 0.0225]);
      else s.add(p.screen, new THREE.BoxGeometry(0.08, 0.034, 0.002), [0, 0.006, 0.022]);
      if (s.fine(2)) for (const x of [-0.03, 0, 0.03]) s.add(p.charcoal, bx(s, 0.016, 0.008, 0.003, 0.001), [x, -0.024, 0.0225]);
    });
    s.at([sideX, floorY + 1.25, zf], [0, 0, 0], () => {
      s.add(p.ssDull, bx(s, 0.07, 0.115, 0.04, 0.004), [0, 0, 0.02]);
      s.add(p.ss, bx(s, 0.076, 0.12, 0.004, 0.002), [0, 0, 0.042]);
      s.add(p.bakelite, bx(s, 0.012, 0.03, 0.014, 0.003), [0, 0.008, 0.05], [0.35, 0, 0]);
      s.add(p.ledRed, cyl(s, 0.006, 0.004, "z", 12), [0, -0.038, 0.045]);
      if (s.fine(2)) for (const y of [-0.05, 0.05]) screwLite(s, p.ssBright, [0, y, 0.044], 0.0022);
    });
    const vent = (zz: number, out: boolean) =>
      s.at([sideX, floorY + 1.92, zz], [0, out ? 0 : Math.PI, 0], () => {
        s.add(p.hdpe, lathe(s, [[0, 0.0005], [0.06, 0.0005], [0.058, 0.012], [0.05, 0.028], [0, 0.03]], 32), [0, 0, 0], [Math.PI / 2, 0, 0]);
        if (s.fine(1)) for (let i = -2; i <= 2; i++) s.add(p.charcoal, bx(s, 0.07 - Math.abs(i) * 0.012, 0.005, 0.004, 0.0015), [0, i * 0.011, 0.029]);
        if (k.freezer && out) s.add(p.ledAmber, cyl(s, 0.004, 0.003, "z", 10), [0.04, -0.035, 0.016]);
      });
    vent(zf, true);
    vent(zi, false);

    // The strip curtain on its rail across the inside of the head.
    const strips: string[] = [];
    {
      const rail = k.doorW + 0.08;
      const ry = dTop + 0.07;
      const rz = zi - 0.022;
      s.add(p.ssDull, bx(s, rail, 0.04, 0.012, 0.002), [k.door0 + k.doorW / 2, ry, zi - 0.0064]);
      s.add(p.ssBright, bx(s, rail - 0.012, 0.012, 0.03, 0.003), [k.door0 + k.doorW / 2, ry - 0.012, rz + 0.004]);
      if (s.fine(2)) for (const f of [0.08, 0.5, 0.92]) screwLite(s, p.ssBright, [k.door0 - 0.04 + rail * f, ry + 0.008, zi - 0.012], 0.003, [Math.PI, 0, 0]);
      const sw = 0.203;
      const n = Math.max(5, Math.round((k.doorW + 0.06) / 0.16));
      const pitch = (k.doorW + 0.06 - sw) / (n - 1);
      const len = ry - 0.02 - floorY - 0.006;
      const pvc = pvcMaterial();
      for (let i = 0; i < n; i++) {
        const x = k.door0 - 0.03 + sw / 2 + i * pitch;
        const z = rz - (i % 2) * 0.004;
        const name = `pvc_strip_${i}`;
        const j = s.joint(name, [x, ry - 0.02, z]);
        j.add(pvc, bx(j, sw, len, 0.002, 0.0008), [0, -len / 2, 0]);
        // The strip folds over a hanger plate: a stainless clamp across its top with two bolts.
        j.add(p.ssBright, bx(j, sw - 0.01, 0.03, 0.005, 0.0015), [0, -0.008, i % 2 ? -0.0042 : 0.0042]);
        // Its bolts run from just behind the strip out through the clamp; a symmetric bolt would end flush with the
        // face of the strip hung behind it.
        if (j.fine(2)) for (const sx of [-0.06, 0.06]) j.add(p.ssBright, cyl(j, 0.004, 0.0095, "z", 8), [sx, -0.008, i % 2 ? -0.00225 : 0.00225]);
        strips.push(name);
      }
    }

    // Unit coolers on the back wall, hung on rod from the ceiling, intake to the wall, with their lines up through
    // the ceiling and drain down the wall.
    const units: { x: number; L: number; H: number; D: number }[] = [];
    for (let u = 0; u < k.units; u++) {
      const cx = (w * (u + 1)) / (k.units + 1);
      const fans = k.fans;
      const pr = fans > 2 ? 0.18 : 0.165;
      const L = fans * (pr * 2 + 0.1) + 0.24;
      const H = 2 * (pr + 0.036) + 0.075;
      const D = 0.5;
      const ox = cx - L / 2;
      const oy = top - 0.07 - H;
      const oz = zb + 0.3 + D;
      s.at([ox, oy, oz], [0, 0, 0], () => {
        unitCooler(s, u, fans, k.freezer, fin);
        // Hanger rails and rods with their nuts, and the ceiling plates.
        for (const z of [-0.08, -D + 0.08]) s.add(p.ssDull, new THREE.BoxGeometry(L + 0.02, 0.03, 0.04), [L / 2, H + 0.015, z]);
        for (const x of [0.06, L - 0.06])
          for (const z of [-0.08, -D + 0.08]) {
            s.add(p.ssDull, new THREE.CylinderGeometry(0.0048, 0.0048, 0.07, 6, 1, true), [x, H + 0.035, z]);
            s.add(p.ssBright, new THREE.CylinderGeometry(0.009, 0.009, 0.01, 6), [x, H + 0.035, z]);
            s.add(p.ssDull, new THREE.BoxGeometry(0.06, 0.005, 0.06), [x, H + 0.0675, z]);
          }
        // Suction line in black closed-cell insulation, liquid line bare, conduit: up from the right end.
        const ex = L + 0.05;
        pipe(s, p.rubber, [[L - 0.005, H * 0.5, -0.3], [ex, H * 0.5, -0.3], [ex, H + 0.07, -0.3]], 0.019, 0.04, 14);
        pipe(s, p.copper, [[L - 0.005, H * 0.38, -0.36], [ex + 0.01, H * 0.38, -0.36], [ex + 0.01, H + 0.07, -0.36]], 0.0048, 0.02, 8);
        pipe(s, p.ssDull, [[L + 0.05, H * 0.55, -0.12], [L + 0.07, H * 0.55, -0.12], [L + 0.07, H + 0.07, -0.12]], 0.009, 0.03, 8);
        // Drain: from the pan's outlet back to the wall and down it, in copper wrapped with heat tape (a freezer)
        // or white PVC (a cooler), out through the wall in a sealed sleeve.
        const dm = k.freezer ? p.copper : p.hdpe;
        pipe(s, dm, [[L - 0.06, 0.0, -D / 2], [L - 0.06, -0.06, -D / 2], [L - 0.06, -0.06, -D - 0.27], [L - 0.06, -0.5, -D - 0.27], [L - 0.06, -0.5, -D - 0.3]], 0.0127, 0.04, 12);
        if (k.freezer && s.fine(2)) {
          const pts: V3[] = [];
          for (let i = 0; i <= 60; i++) {
            const a = i * 0.9;
            pts.push([L - 0.06 + Math.cos(a) * 0.0135, -0.08 - i * 0.0068, -D - 0.27 + Math.sin(a) * 0.0135]);
          }
          s.add(p.charcoal, bar(s, pts, 0.0018, 0.004, 5));
        }
        s.add(p.ssBright, cyl(s, 0.025, 0.006, "z", 16), [L - 0.06, -0.5, -D - 0.297]);
      });
      units.push({ x: ox, L, H, D });
      // Through the ceiling: a pitch-pocket boot on the roof.
      roof.add(p.rubber, cyl(roof, 0.045, 0.05, "y", 20, 0.03), [ox + L + 0.06, h + 0.025, zb + 0.3 + D - 0.3]);
    }

    // The condensing unit on the roof, its lines across to the boots.
    const cu = { x: Math.min(w - 1.05, w / 2 + 0.15), z: -d + 0.72 };
    roof.at([cu.x, h, cu.z], [0, 0, 0], () => condensingUnit(roof, fin === null ? p.alu : p.alu));
    for (const un of units) {
      const bx0 = un.x + un.L + 0.06;
      const bz0 = zb + 0.3 + un.D - 0.3;
      const yr = h + 0.07;
      pipe(roof, p.rubber, [[bx0, h + 0.03, bz0], [bx0, yr, bz0], [cu.x + 0.1, yr, bz0], [cu.x + 0.1, yr, cu.z + 0.03], [cu.x + 0.1, h + 0.4, cu.z + 0.03]], 0.019, 0.05, 14);
      pipe(roof, p.copper, [[bx0 + 0.01, h + 0.03, bz0 - 0.06], [bx0 + 0.01, yr - 0.02, bz0 - 0.06], [cu.x + 0.22, yr - 0.02, bz0 - 0.06], [cu.x + 0.22, yr - 0.02, cu.z + 0.03], [cu.x + 0.22, h + 0.4, cu.z + 0.03]], 0.0048, 0.03, 8);
    }

    // The vapour-proof light: a sealed polycarbonate tube fixture on the ceiling, stainless clips.
    {
      const lz = (zi + zb) / 2 + 0.25;
      const ll = Math.min(1.22, (zi - zb) * 0.45);
      roof.add(p.poly, bx(roof, 0.1, 0.05, ll, 0.008), [w / 2, top - 0.025, lz]);
      if (roof.detail > 0) {
        const lens = new THREE.CylinderGeometry(0.045, 0.045, ll - 0.04, roof.seg(20, 6), 1, false, Math.PI / 2, Math.PI);
        lens.rotateX(Math.PI / 2);
        roof.add(p.ledWhite, lens, [w / 2, top - 0.05, lz]);
        if (roof.fine(2)) for (let i = 0; i < 6; i++) roof.add(p.ssBright, bx(roof, 0.1, 0.012, 0.016, 0.002), [w / 2, top - 0.052, lz - ll / 2 + 0.07 + (i * (ll - 0.14)) / 5], [0, 0, 0]);
        for (const sz of [-1, 1]) roof.add(p.poly, bx(roof, 0.096, 0.09, 0.02, 0.006), [w / 2, top - 0.05, lz + (sz * (ll - 0.02)) / 2]);
      }
    }

    // Shelving down both sides, and along the back if asked for.
    const Ds = w < 2.7 ? 0.457 : 0.61;
    const run = zi - 0.06 - (k.back ? zb + BACK_SHELF.gap + BACK_SHELF.D + 0.52 : zb + 0.03);
    if (k.back) for (const x of BACK_SHELF.x0) s.at([x, floorY, zb + BACK_SHELF.gap + BACK_SHELF.D], [0, 0, 0], () => wireShelving(s, BACK_SHELF.L, BACK_SHELF.D, BACK_SHELF.H, { grooves: false, lite: true, tiers: [...BACK_SHELF.tiers] }));
    const lens: number[] = [];
    {
      let left = run;
      for (const L of [1.219, 1.219, 1.219, 0.914, 0.762, 0.61]) if (L <= left + 0.001 && lens.length < (w > 2.7 ? 2 : 3)) {
        lens.push(L);
        left -= L + 0.01;
      }
    }
    let z0 = zi - 0.06;
    lens.forEach((L, i) => {
      s.at([T + 0.012 + Ds, floorY, z0], [0, Math.PI / 2, 0], () => wireShelving(s, L, Ds, 1.88, { grooves: false, lite: true, stock: [s], fill: w > 2.7 ? 0.4 : 0.55, seed: 101 + i, cartons: !k.back }));
      wallR.at([w - T - 0.012 - Ds, floorY, z0 - L], [0, -Math.PI / 2, 0], () => wireShelving(wallR, L, Ds, 1.88, { grooves: false, lite: true, stock: [wallR], fill: w > 2.7 ? 0.35 : 0.5, seed: 201 + i, cartons: !k.back }));
      z0 -= L + 0.01;
    });

    return (b: Built) => {
      const J = (n: string) => b.joints.get(n);
      const door = J("cam_lift_door");
      const door0 = door ? door.position.clone() : new THREE.Vector3();
      const body = J("closer_body");
      const rod = J("closer_rod");
      const evapFans: THREE.Object3D[] = [];
      for (let u = 0; u < k.units; u++) for (let i = 0; i < k.fans; i++) {
        const f = J(`evap_fan_${u}_${i}`);
        if (f) evapFans.push(f);
      }
      const condFan = J("cond_fan");
      const stripJ = strips.map((n) => J(n)).filter((x): x is THREE.Group => !!x);
      const frostJ = Array.from({ length: k.units }, (_, u) => J(`part_frost_${u}`));
      let open = 0;
      let sway = 0;
      let spin = 0;
      let cut = false;
      let roofOff = false;
      const eu = new THREE.Euler();
      const q = new THREE.Quaternion();
      const setDoor = (f: number) => {
        open = Math.max(0, Math.min(1, f));
        const th = open * 1.92;
        turnJoint(door, "y", -th);
        if (door) door.position.set(door0.x, door0.y + 0.024 * (th / Math.PI), door0.z);
        // The closer follows: its body swings to aim at the door's bracket and the rod draws out.
        const Bx = rig.hx + rig.bx * Math.cos(th) - rig.bz * Math.sin(th);
        const Bz = rig.hz + rig.bx * Math.sin(th) + rig.bz * Math.cos(th);
        const dx = Bx - rig.ax;
        const dz = Bz - rig.az;
        turnJoint(body, "y", Math.atan2(dz, -dx));
        slideJoint(rod, "x", -(Math.hypot(dx, dz) - rig.rod0), [0, 0, 0]);
      };
      const setStrips = (t: number) => {
        stripJ.forEach((j, i) => {
          // Cold air spills out under the strips and warm air comes in over them: the feet blow out, flutter, part.
          const a = -sway * (0.07 + 0.035 * Math.sin(1.9 * t + i * 1.7) + 0.02 * Math.sin(3.3 * t + i * 0.6));
          const part = sway * 0.035 * Math.sin(1.1 * t + i * 2.3) * (i % 2 ? 1 : -1);
          const twist = sway * 0.12 * Math.sin(1.4 * t + i);
          q.setFromEuler(eu.set(a, twist, part));
          const rest = j.userData.rest as THREE.Quaternion | undefined;
          if (rest) j.quaternion.copy(rest).multiply(q);
          else j.quaternion.copy(q);
        });
      };
      const vis = () => {
        const r = J("part_roof");
        if (r) r.visible = !(cut || roofOff);
        for (const n of ["part_wall_right", "part_front_right"]) {
          const g = J(n);
          if (g) g.visible = !cut;
        }
      };
      const setFrostAll = (f: number) => {
        setFrost(fin, f);
        const crust = Math.max(0, (f - 0.5) / 0.5);
        for (const g of frostJ) if (g) {
          g.visible = crust > 0.01;
          g.scale.set(1, 1, Math.max(0.01, crust));
        }
      };
      setFrostAll(0);
      readout?.set(k.freezer ? -20 : 3);
      const hooks: Hooks = {
        open: setDoor,
        set: (name, v) => {
          if (name === "door") setDoor(v);
          else if (name === "frost") setFrostAll(v);
          else if (name === "cutaway") {
            cut = v > 0.5;
            vis();
          } else if (name === "roof") {
            roofOff = v > 0.5;
            vis();
          } else if (name === "temp") readout?.set(v);
        },
        tick: (t, dt, load) => {
          spin += dt * Math.PI * 2 * 5 * Math.max(0, Math.min(1, load));
          for (const f of evapFans) turnJoint(f, "z", -spin);
          turnJoint(condFan, "z", -spin * 0.8);
          sway += (open - sway) * (1 - Math.exp(-dt / 0.5));
          setStrips(t);
        },
      };
      return hooks;
    };
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Cabinets: the kit. Reach-ins, the undercounter worktop, the makeline and the blast chiller share one door,     */
/* one body and one way of breathing; each builder places them.                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** A cabinet door's liner plug lies this far in from the slab's edges; the opening it fills is 3 mm bigger all round. */
const PLUG = 0.034;
/** A door gasket's depth, slab back to sealing face: the face frame lies this far behind a closed door. */
const SEAL = 0.0125;

/**
 * A magnetic door gasket round a rounded rectangle `w` by `h` (its inner edge) in the current frame's xy: the PVC
 * extrusion's dart in the door, two bellows and the bulb round the magnet strip, its sealing face toward -z, the
 * door's back at z = 0.
 */
function doorGasket(s: Shop, w: number, h: number, r = 0.01, gw = 0.019): void {
  if (s.detail === 0) return;
  const dp = SEAL;
  const prof: [number, number][] = s.fine(2)
    ? [[0.003, 0], [gw - 0.006, 0], [gw - 0.003, 0.004], [gw + 0.001, dp - 0.004], [gw - 0.001, dp], [0.001, dp], [-0.001, dp - 0.004], [0.002, 0.005], [0.003, 0]]
    : [[0.002, 0], [gw - 0.004, 0], [gw, dp], [0, dp], [0.002, 0]];
  const rings: Ring[] = prof.map(([o, y]) => ({ w: w + 2 * o, d: h + 2 * o, r: r + o, y }));
  s.add(palette().gasket, loft(s, rings, 0, 0, s.fine(2) ? 3 : 1), [0, 0, 0.0008], [-Math.PI / 2, 0, 0]);
}

interface CabDoorOpts {
  /** The slab's width, height and thickness. */
  w: number;
  h: number;
  t: number;
  /** +1: hinged on its left, the free edge toward +x; -1: hinged on its right. */
  side: 1 | -1;
  /** How far the slab's hinge edge lies beyond the pin. */
  inset: number;
  /** A bar pull standing off the face near the free edge, or a pull extruded along the top edge. */
  handle: "bar" | "lip";
  /** Handle length; a bar's centre height above the pin and its distance in from the free edge. */
  hl: number;
  hy?: number;
  hx?: number;
  lock?: boolean;
}

/**
 * A foamed cabinet door in its joint, the origin on the hinge pin at the slab's mid-thickness and mid-height: a
 * stainless slab with rolled edges, the liner's plug stepped to fill the opening, the gasket round the plug, the
 * pull and a cylinder lock, the hinge leaves on the top and bottom edges with their pins.
 */
function cabDoor(j: Shop, o: CabDoorOpts): void {
  const p = palette();
  const { w, h, t, side } = o;
  const X = (x: number) => side * x;
  const cx = X(w / 2 - o.inset);
  const hx = X(w - o.inset - (o.hx ?? 0.035));
  const hy = o.hy ?? 0;
  j.add(p.ss, bx(j, w, h, t, 0.005), [cx, 0, 0]);
  if (j.detail === 0) {
    if (o.handle === "bar") j.add(p.ssBright, new THREE.BoxGeometry(0.02, o.hl, 0.026), [hx, hy, t / 2 + 0.013]);
    return;
  }
  const pd = 0.022;
  j.add(p.ssDull, bx(j, w - 2 * PLUG, h - 2 * PLUG, pd + 0.004, 0.006), [cx, 0, -t / 2 - pd / 2 + 0.002]);
  j.at([cx, 0, -t / 2], [0, 0, 0], () => doorGasket(j, w - 2 * PLUG + 0.004, h - 2 * PLUG + 0.004));
  if (o.handle === "bar") {
    const off = 0.028;
    const r = 0.0095;
    const [y0, y1] = [hy - o.hl / 2, hy + o.hl / 2];
    pipe(j, p.ssBright, [[hx, y0, t / 2], [hx, y0, t / 2 + off], [hx, y1, t / 2 + off], [hx, y1, t / 2]], r, 0.016, 14);
    for (const y of [y0, y1]) j.add(p.ssBright, cyl(j, r * 1.6, 0.005, "z", 20), [hx, y, t / 2 + 0.0025]);
    if (o.lock && j.fine(2)) {
      const ly = y0 - 0.055;
      j.add(p.ssBright, cyl(j, 0.0105, 0.006, "z", 20), [hx, ly, t / 2 + 0.003]);
      j.add(p.charcoal, new THREE.BoxGeometry(0.0018, 0.009, 0.0008), [hx, ly, t / 2 + 0.0064]);
    }
  } else {
    // A flange screwed to the face, a lip out over the fingers, a drip edge down: the fingers hook behind it.
    const sh = new THREE.Shape([[0, -0.042], [0.004, -0.042], [0.004, -0.005], [0.02, -0.005], [0.02, -0.018], [0.024, -0.018], [0.024, 0], [0, 0]].map(([a, b]) => new THREE.Vector2(a, b)));
    j.add(p.ssBright, plate(j, sh, o.hl), [cx + o.hl / 2, h / 2 - 0.004, t / 2], [0, -Math.PI / 2, 0]);
  }
  for (const sy of [-1, 1]) {
    j.add(p.ssBright, bx(j, 0.05, 0.004, 0.03, 0.0012), [X(0.015), sy * (h / 2 + 0.002), 0]);
    j.add(p.ssBright, cyl(j, 0.0055, 0.012, "y", 12), [0, sy * (h / 2 + 0.006), 0]);
  }
}

/**
 * The cabinet's half of a door's hinges: a bracket off the face frame out to each pin, its leaf screwed to the
 * frame; under a worktop (its underside at `roof`) the top bracket is a block screwed up into the top instead.
 */
function hingeBrackets(s: Shop, px: number, side: 1 | -1, y0: number, y1: number, pz: number, zf: number, roof?: number): void {
  const p = palette();
  const L = pz - zf + 0.014;
  const x = px + side * 0.008;
  for (const [y, sy] of [[y1, 1], [y0, -1]] as const) {
    if (sy > 0 && roof !== undefined) {
      const hb = roof - (y + 0.008);
      s.add(p.ssBright, bx(s, 0.046, hb, L, 0.0015), [x, y + 0.008 + hb / 2, zf + L / 2]);
      continue;
    }
    s.add(p.ssBright, bx(s, 0.046, 0.007, L, 0.0015), [x, y + sy * 0.0115, zf + L / 2]);
    s.add(p.ssBright, bx(s, 0.046, 0.016, 0.003, 0.001), [x, y + sy * 0.016, zf + 0.0015]);
    if (s.fine(2)) for (const dx of [-0.013, 0.013]) screwLite(s, p.ss, [x + dx, y + sy * 0.016, zf + 0.003], 0.0022);
  }
}

interface DoorPlan {
  name: string;
  side: 1 | -1;
  /** The pin's x, and the slab's left edge and width. */
  px: number;
  xa: number;
  w: number;
}
/** One door hinged left, or a pair hinged at the outer edges meeting in the middle `gap` apart. */
function planDoors(w: number, names: string[], pinX: number, inset: number, gap = 0.006): DoorPlan[] {
  if (names.length === 1) return [{ name: names[0]!, side: 1, px: pinX, xa: pinX - inset, w: w - 2 * (pinX - inset) }];
  const dw = w / 2 - gap / 2 - (pinX - inset);
  return [
    { name: names[0]!, side: 1, px: pinX, xa: pinX - inset, w: dw },
    { name: names[1]!, side: -1, px: w - pinX, xa: w - pinX + inset - dw, w: dw },
  ];
}
/** The face-frame opening a door's plug fills: [x0, x1, y0, y1]. */
function openingOf(o: DoorPlan, y0: number, y1: number): [number, number, number, number] {
  const k = PLUG - 0.003;
  return [o.xa + k, o.xa + o.w - k, y0 + k, y1 - k];
}

interface CabBody {
  w: number;
  d: number;
  y0: number;
  y1: number;
  /** The face frame's front plane. */
  zf: number;
  openings: [number, number, number, number][];
  /** Skin the top (a makeline carries its rail there instead). */
  top?: boolean;
}
interface Cavity {
  x0: number;
  x1: number;
  y0: number;
  y1: number;
  /** The liner's back, and its front where the breakers end. */
  zb: number;
  zl: number;
}

/**
 * A foamed cabinet body, x 0..w, y y0..y1, from its back at z = -d to the face frame's front at z = zf: stainless
 * ends and top, a galvanised back and bottom, the face frame cut for the door openings, black ABS breaker strips
 * lining each opening back to the stainless liner, the mullion between a pair. Returns the cavity in the liner.
 */
function cabBody(s: Shop, o: CabBody): Cavity {
  const p = palette();
  const { w, d, y0, y1, zf } = o;
  const t = 0.0012;
  const H = y1 - y0;
  const zs = zf - 0.0015;
  const sd = d + zs;
  const zc = -d + sd / 2;
  for (const x of [t / 2, w - t / 2]) s.add(p.ss, bx(s, t, H, sd, 0.0005), [x, y0 + H / 2, zc]);
  if (o.top ?? true) s.add(p.ss, bx(s, w - 2 * t, t, sd, 0.0005), [w / 2, y1 - t / 2, zc]);
  s.add(p.ssDull, new THREE.BoxGeometry(w - 2 * t, t, sd), [w / 2, y0 + t / 2, zc]);
  s.add(p.alu, new THREE.BoxGeometry(w - 2 * t, H - 2 * t, t), [w / 2, y0 + H / 2, -d + t / 2]);
  // The face frame: one sheet across the front, the openings cut from it.
  const ff = roundRect(w, H, 0.0015, w / 2, y0 + H / 2);
  for (const [a, b, c, e] of o.openings) ff.holes.push(rectPath(new THREE.Path(), b - a, e - c, 0.002, (a + b) / 2, (c + e) / 2, true));
  s.add(p.ss, plate(s, ff, 0.0015), [0, 0, zs]);
  // Breakers: the cold liner never touches the warm skin; ABS strips bridge the throat between them.
  const bd = 0.045;
  const zl = zs - bd;
  const zm = zs - bd / 2;
  for (const [a, b, c, e] of o.openings) {
    for (const y of [c + 0.001, e - 0.001]) s.add(p.plastic, new THREE.BoxGeometry(b - a, 0.002, bd), [(a + b) / 2, y, zm]);
    for (const x of [a + 0.001, b - 0.001]) s.add(p.plastic, new THREE.BoxGeometry(0.002, e - c - 0.004, bd), [x, (c + e) / 2, zm]);
  }
  const cav: Cavity = {
    x0: Math.min(...o.openings.map((q) => q[0])),
    x1: Math.max(...o.openings.map((q) => q[1])),
    y0: Math.min(...o.openings.map((q) => q[2])),
    y1: Math.max(...o.openings.map((q) => q[3])),
    zb: -d + 0.055,
    zl,
  };
  const ld = cav.zl - cav.zb;
  const lz = (cav.zl + cav.zb) / 2;
  const lw = cav.x1 - cav.x0;
  const lh = cav.y1 - cav.y0;
  for (const x of [cav.x0 - t / 2, cav.x1 + t / 2]) s.add(p.ssDull, new THREE.BoxGeometry(t, lh + 2 * t, ld), [x, cav.y0 + lh / 2, lz]);
  for (const y of [cav.y0 - t / 2, cav.y1 + t / 2]) s.add(p.ssDull, new THREE.BoxGeometry(lw, t, ld), [cav.x0 + lw / 2, y, lz]);
  s.add(p.ssDull, new THREE.BoxGeometry(lw + 2 * t, lh + 2 * t, t), [cav.x0 + lw / 2, cav.y0 + lh / 2, cav.zb - t / 2]);
  const ops = [...o.openings].sort((a, b) => a[0] - b[0]);
  for (let i = 0; i + 1 < ops.length; i++) {
    const a = ops[i]![1];
    const b = ops[i + 1]![0];
    if (b - a > 0.01) s.add(p.ssDull, new THREE.BoxGeometry(b - a, lh, t), [(a + b) / 2, cav.y0 + lh / 2, zl - t / 2]);
  }
  return cav;
}

/** A slotted pilaster strip from y0, `h` tall, its back at z = 0 of the current frame and its face toward +z. */
function pilaster(s: Shop, y0: number, h: number): void {
  const p = palette();
  s.add(p.ss, new THREE.BoxGeometry(0.022, h, 0.005), [0, y0 + h / 2, 0.0025]);
  if (!s.fine(2)) return;
  for (let y = y0 + 0.02; y < y0 + h - 0.02; y += 0.0254) s.add(p.charcoal, new THREE.PlaneGeometry(0.0045, 0.013), [0, y, 0.0056]);
}

/**
 * An epoxy-coated wire shelf `w` by `d`, its top at y = 0, centred in the current frame: a 6 mm frame round deck
 * wires front to back at an inch, two cross supports welded under them.
 */
function epoxyShelf(s: Shop, w: number, d: number): void {
  const p = palette();
  const m = p.poly;
  if (s.detail === 0) {
    s.add(m, new THREE.BoxGeometry(w, 0.006, d), [0, -0.003, 0]);
    return;
  }
  const fine = s.fine(2);
  const fr = 0.003;
  const wr = 0.0021;
  const rad = fine ? 8 : 4;
  for (const z of [d / 2 - fr, -d / 2 + fr]) wire(s, m, [-w / 2, -fr, z], [w / 2, -fr, z], fr, rad);
  for (const x of [w / 2 - fr, -w / 2 + fr]) wire(s, m, [x, -fr, -d / 2 + 2 * fr], [x, -fr, d / 2 - 2 * fr], fr, rad);
  const n = Math.max(3, Math.round((w - 0.02) / (fine ? 0.0254 : 0.076)));
  for (let i = 1; i < n; i++) {
    const x = -w / 2 + (i * w) / n;
    wire(s, m, [x, -wr, -d / 2 + fr], [x, -wr, d / 2 - fr], wr, fine ? 6 : 3);
  }
  for (const z of [-d / 4, d / 4]) wire(s, m, [-w / 2 + fr, -2 * wr - 0.0025, z], [w / 2 - fr, -2 * wr - 0.0025, z], 0.0025, fine ? 6 : 3);
}

/** Shelves on pilasters down the cavity's walls (and up the middle of a two-section one), `tiers` the shelf tops. */
function cabShelves(s: Shop, c: Cavity, sections: [number, number][], tiers: number[], top: number): void {
  const zf = c.zl - 0.035;
  const zb = c.zb + 0.03;
  const ph = top - c.y0 - 0.01;
  for (const z of [zf, zb]) {
    s.at([c.x0, 0, z], [0, Math.PI / 2, 0], () => pilaster(s, c.y0 + 0.005, ph));
    s.at([c.x1, 0, z], [0, -Math.PI / 2, 0], () => pilaster(s, c.y0 + 0.005, ph));
    for (let i = 0; i + 1 < sections.length; i++) {
      const xm = (sections[i]![1] + sections[i + 1]![0]) / 2;
      s.add(palette().ss, new THREE.BoxGeometry(0.012, ph, 0.024), [xm, c.y0 + 0.005 + ph / 2, z]);
    }
  }
  for (const [a, b] of sections)
    for (const y of tiers) {
      s.at([(a + b) / 2, y, (zf + zb) / 2], [0, 0, 0], () => epoxyShelf(s, b - a - 0.014, zf - zb + 0.03));
      if (s.fine(2))
        for (const x of [a + 0.008, b - 0.008]) for (const z of [zf, zb]) s.add(palette().ssBright, new THREE.BoxGeometry(0.008, 0.016, 0.02), [x, y - 0.012, z + (z === zf ? -0.012 : 0.012)]);
    }
}

/**
 * A cabinet's evaporator housing hung under its ceiling at the back, x 0..L of the current frame, its top at y = 0,
 * its face at z = 0 and its back D behind: the coil across the back drawing air up through slots in the floor,
 * `fans` fans blowing out over the shelves through wire guards in the face.
 */
function cabEvaporator(s: Shop, L: number, H: number, D: number, fans: number, name: (i: number) => string): void {
  const p = palette();
  const fr = Math.min(0.085, H / 2 - 0.018);
  const fx = (i: number) => (L * (i + 0.5)) / fans;
  if (s.detail > 0) {
    const sh = roundRect(L, H, 0.003);
    for (let i = 0; i < fans; i++) sh.holes.push(circleHoleCW(fx(i) - L / 2, 0, fr + 0.005));
    s.add(p.ss, plate(s, sh, 0.0012), [L / 2, -H / 2, -0.0012]);
  } else s.add(p.ss, new THREE.BoxGeometry(L, H, 0.002), [L / 2, -H / 2, -0.001]);
  s.add(p.ss, new THREE.BoxGeometry(L, 0.0012, D), [L / 2, -H + 0.0006, -D / 2]);
  for (const x of [0.0006, L - 0.0006]) s.add(p.ss, new THREE.BoxGeometry(0.0012, H - 0.0012, D), [x, -H / 2 + 0.0006, -D / 2]);
  if (s.detail === 0) return;
  if (s.fine(2)) {
    // The coil runs the housing's length; only the lengths behind the fans can be seen, so only they are drawn.
    for (let i = 0; i < fans; i++) s.at([fx(i) - fr - 0.03, -H + 0.006, -D + 0.07], [0, 0, 0], () => finCoil(s, p.alu, 2 * fr + 0.06, H - 0.02, 0.05, 6, false));
    for (let x = 0.03; x < L - 0.04; x += 0.03) s.add(p.charcoal, new THREE.PlaneGeometry(0.016, 0.05).rotateX(Math.PI / 2), [x + 0.008, -H - 0.0006, -D + 0.045]);
  }
  for (let i = 0; i < fans; i++) {
    const x = fx(i);
    const fj = s.joint(name(i), [x, -H / 2, -0.022]);
    if (s.fine(2)) {
      s.add(p.charcoal, cyl(s, 0.03, 0.05, "z", 20), [x, -H / 2, -0.07]);
      axialFan(fj, p.alu, p.ssDull, fr - 0.006, 4);
    } else s.add(p.charcoal, new THREE.CircleGeometry(fr, 10), [x, -H / 2, -0.03]);
    s.at([x, -H / 2, 0.002], [0, 0, 0], () => fanGuard(s, p.chrome, fr + 0.004, 0.008, 0.0011));
  }
}

/**
 * The condenser behind a louvred grille facing +z at z = zg, the grille from x = gx to x1 and y0 to y1: the coil just
 * behind the louvres, the shroud and plate the fan (a joint, discharging back) turns in, its motor and guard.
 */
function grilleCondenser(s: Shop, gx: number, x1: number, y0: number, y1: number, zg: number, fanName: string, cwMax: number): { cx0: number; zc0: number } {
  const p = palette();
  const gw = x1 - gx;
  const H = y1 - y0;
  s.at([gx + gw / 2, y0 + H / 2, zg], [0, 0, 0], () => louvres(s, p.ss, gw, H, s.fine(1) ? Math.max(3, Math.round(H / 0.018)) : 3, 0.014));
  const cw = Math.min(cwMax, gw - 0.06);
  const ch = H - 0.024;
  const cx0 = x1 - 0.03 - cw;
  const zc0 = zg - 0.024;
  if (s.detail === 0) {
    s.add(p.charcoal, new THREE.BoxGeometry(gw - 0.01, H - 0.01, 0.002), [gx + gw / 2, y0 + H / 2, zg - 0.03]);
    return { cx0, zc0 };
  }
  s.at([cx0, y0 + 0.012, zc0], [0, 0, 0], () => finCoil(s, p.alu, cw, ch, 0.032, 12, true));
  const zs = zc0 - 0.032;
  const zp = zs - 0.026;
  const fr = Math.min(0.08, ch / 2 - 0.006);
  for (const y of [y0 + 0.012, y0 + 0.012 + ch]) s.add(p.charcoal, new THREE.BoxGeometry(cw, 0.001, zs - zp), [cx0 + cw / 2, y, (zs + zp) / 2]);
  for (const x of [cx0, cx0 + cw]) s.add(p.charcoal, new THREE.BoxGeometry(0.001, ch, zs - zp), [x, y0 + 0.012 + ch / 2, (zs + zp) / 2]);
  const fx = cx0 + cw / 2;
  const fy = y0 + 0.012 + ch / 2;
  const sh = roundRect(cw, ch, 0.002);
  sh.holes.push(circleHoleCW(0, 0, fr + 0.005));
  s.add(p.charcoal, plate(s, sh, 0.0012), [fx, fy, zp - 0.0012]);
  const fj = s.joint(fanName, [fx, fy, zp - 0.014], [0, Math.PI, 0]);
  axialFan(fj, p.enamel, p.enamel, fr - 0.004, 4);
  if (s.fine(2)) {
    s.add(p.charcoal, cyl(s, 0.028, 0.045, "z", 20), [fx, fy, zp - 0.05]);
    s.at([fx, fy, zp - 0.076], [0, Math.PI, 0], () => fanGuard(s, p.enamel, fr + 0.004, 0.005, 0.0012));
  }
  return { cx0, zc0 };
}

/**
 * A top-mount condensing unit in its compartment on a cabinet's roof, x 0..w, y y0..y1, from the back at z = -d to
 * its face at z = zg, overhanging the face frame at zf: stainless ends and cover, the control bezel left of x = gx
 * (`bezel` draws on it, centred, facing +z), the louvred grille across the rest with the condenser behind it, the
 * compressor, drier and copper in the bay behind the bezel, the back open for the fan's air.
 */
function topCondenser(s: Shop, w: number, d: number, y0: number, y1: number, zg: number, zf: number, gx: number, fanName: string, bezel: () => void): void {
  const p = palette();
  const H = y1 - y0;
  const t = 0.0012;
  const L = d + zg;
  const zc = -d + L / 2;
  s.add(p.ss, bx(s, w, t, L, 0.0005), [w / 2, y1 - t / 2, zc]);
  for (const x of [t / 2, w - t / 2]) s.add(p.ss, bx(s, t, H - t, L, 0.0005), [x, y0 + (H - t) / 2, zc]);
  s.add(p.ss, new THREE.BoxGeometry(w - 2 * t, t, zg - zf), [w / 2, y0 + t / 2, (zg + zf) / 2]);
  s.add(p.ss, bx(s, gx - t, H - t, 0.0012, 0.0005), [t + (gx - t) / 2, y0 + (H - t) / 2, zg - 0.0006]);
  s.at([gx / 2, y0 + H / 2, zg], [0, 0, 0], bezel);
  const { cx0, zc0 } = grilleCondenser(s, gx, w - t, y0 + t, y1 - t, zg, fanName, 0.42);
  if (s.detail === 0) return;
  // The compressor in the bay behind the bezel, its discharge to the coil's top header, the liquid line back from
  // the coil through the drier and down into the cabinet, the insulated suction line up out of it.
  const cr = 0.058;
  const ch = (H - 0.016) / 1.12;
  const cpx = Math.max(gx * 0.5, (t + cx0) / 2);
  const cpz = -d * 0.55;
  s.at([cpx, y0, cpz], [0, 0, 0], () => compressor(s, cr, ch));
  const dy = y0 + ch * 1.12;
  pipe(s, p.copper, [[cpx + cr * 0.9, dy, cpz - cr * 0.3], [cx0 - 0.03, dy, cpz - cr * 0.3], [cx0 - 0.03, dy, zc0 - 0.016], [cx0 - 0.012, dy, zc0 - 0.016]], 0.0048, 0.02, 8);
  const lx = cx0 - 0.055;
  const ly = y0 + 0.028;
  pipe(s, p.copper, [[cx0 - 0.012, ly, zc0 - 0.016], [lx, ly, zc0 - 0.016], [lx, ly, zc0 - 0.07]], 0.0032, 0.015, 8);
  s.add(p.copper, cyl(s, 0.0135, 0.1, "z", 16), [lx, ly, zc0 - 0.12]);
  pipe(s, p.copper, [[lx, ly, zc0 - 0.17], [lx, ly, -d * 0.85], [lx, y0 - 0.02, -d * 0.85]], 0.0032, 0.015, 8);
  const sx = cpx - cr * 1.1;
  const sy = y0 + ch * 0.95;
  const sz = cpz + cr;
  pipe(s, p.copper, [[sx, sy, sz], [sx, sy, sz + 0.06]], 0.0064, 0.02, 10);
  pipe(s, p.rubber, [[sx, sy, sz + 0.03], [sx, sy, sz + 0.09], [sx, y0 - 0.02, sz + 0.09]], 0.0115, 0.025, 12);
  s.add(p.charcoal, bx(s, 0.09, 0.08, 0.06, 0.003), [t + 0.046, y0 + 0.05, -d + 0.1]);
  if (s.fine(2)) pipe(s, p.plastic, [[t + 0.09, y0 + 0.05, -d + 0.1], [cpx - 0.01, y0 + 0.05, -d + 0.1], [cpx, y0 + ch * 0.5, cpz + cr * 0.95 + 0.02]], 0.004, 0.03, 6);
}

/**
 * A front-breathing condensing unit in the bay under a counter-height cabinet, x 0..w, y y0..y1, from z = -d to the
 * louvred kick panel at z = zg: the control bezel left of x = gx, the low condenser coil behind the louvres with its
 * fan, the compressor at the back right (x = cx) rising into the doghouse the caller lines in the cabinet's floor.
 */
function kickCondenser(s: Shop, w: number, d: number, y0: number, y1: number, zg: number, gx: number, cx: number, fanName: string, bezel: () => void): void {
  const p = palette();
  const H = y1 - y0;
  const t = 0.0012;
  const L = d + zg;
  const zc = -d + L / 2;
  for (const x of [t / 2, w - t / 2]) s.add(p.ss, bx(s, t, H, L, 0.0005), [x, y0 + H / 2, zc]);
  // The floor pan stops behind the bezel and grille.
  s.add(p.ssDull, new THREE.BoxGeometry(w - 2 * t, t, L - 0.0014), [w / 2, y0 + t / 2, zc - 0.0007]);
  s.add(p.alu, new THREE.BoxGeometry(w - 2 * t, H - t, t), [w / 2, y0 + (H + t) / 2, -d + t / 2]);
  s.add(p.ss, bx(s, gx - t, H, 0.0012, 0.0005), [t + (gx - t) / 2, y0 + H / 2, zg - 0.0006]);
  s.at([gx / 2, y0 + H / 2, zg], [0, 0, 0], bezel);
  const { cx0, zc0 } = grilleCondenser(s, gx, w - t, y0 + t, y1, zg, fanName, 0.42);
  if (!s.fine(2)) return;
  const cr = 0.055;
  const cz = -d + 0.16;
  s.at([cx, y0 + t, cz], [0, 0, 0], () => compressor(s, cr, 0.15, 18));
  pipe(s, p.copper, [[cx + cr * 0.9, y0 + 0.168, cz - cr * 0.3], [cx + cr * 0.9, y0 + 0.168, cz + 0.1], [cx0 + 0.3, y0 + 0.07, zc0 - 0.05], [cx0 + 0.02, y0 + 0.07, zc0 - 0.016], [cx0 - 0.012, y0 + 0.07, zc0 - 0.016]], 0.0048, 0.02, 8);
}

/** A 6 in stainless leg on its adjustable bullet foot under a cabinet whose base is at y = yb. */
function leg6(s: Shop, x: number, z: number, yb: number): void {
  const p = palette();
  if (s.detail === 0) {
    s.add(p.ss, new THREE.CylinderGeometry(0.02, 0.02, yb, 6), [x, yb / 2, z]);
    return;
  }
  if (s.fine(2)) roundLeg(s, p.ss, x, z, 0, yb, 0.0195);
  else {
    s.add(p.ss, new THREE.CylinderGeometry(0.0195, 0.0195, yb - 0.06, 8, 1, true), [x, (yb + 0.06) / 2, z]);
    s.add(p.ss, new THREE.CylinderGeometry(0.015, 0.02, 0.06, 8), [x, 0.03, z]);
  }
  s.add(p.ssDull, bx(s, 0.07, 0.006, 0.07, 0.002), [x, yb - 0.0035, z]);
}

/** A rocker switch facing +z at the current origin, its lamp lit. */
function rocker(s: Shop, at: V3): void {
  const p = palette();
  s.add(p.plastic, bx(s, 0.024, 0.034, 0.006, 0.002), [at[0], at[1], at[2] + 0.003]);
  s.add(p.bakelite, bx(s, 0.016, 0.026, 0.008, 0.003), [at[0], at[1], at[2] + 0.006], [0.18, 0, 0]);
  if (s.fine(1)) s.add(p.ledGreen, new THREE.PlaneGeometry(0.008, 0.004), [at[0], at[1] + 0.008, at[2] + 0.0111]);
}

/**
 * Hooks for a cabinet's doors (swung outward about their pins by up to `max`), its fans and its readout. set():
 * door (all), or a door's joint name; temp °C.
 */
function cabHooks(b: Built, doors: DoorPlan[], max: number, fans: [string, number][], readout: { set: (c: number) => void } | null, temp: number): Hooks & { setDoor: (i: number, f: number) => void } {
  const dj = doors.map((o) => ({ j: b.joints.get(o.name), side: o.side, name: o.name }));
  const fj = fans.map(([n, k]) => [b.joints.get(n), k] as const);
  let spin = 0;
  const setDoor = (i: number, f: number) => {
    const o = dj[i];
    if (o) turnJoint(o.j, "y", -o.side * Math.max(0, Math.min(1, f)) * max);
  };
  readout?.set(temp);
  return {
    setDoor,
    open: (f) => dj.forEach((_, i) => setDoor(i, f)),
    set: (name, v) => {
      if (name === "door") dj.forEach((_, i) => setDoor(i, v));
      else if (name === "temp") readout?.set(v);
      else {
        const i = dj.findIndex((o) => o.name === name);
        if (i >= 0) setDoor(i, v);
      }
    },
    tick: (_t, dt, load) => {
      spin += dt * Math.PI * 2 * 5 * Math.max(0, Math.min(1, load));
      for (const [f, k] of fj) turnJoint(f, "z", -spin * k);
    },
  };
}

/**
 * A sheet swept along x (from -len/2 to len/2) through a profile of (z, y) points, smooth along the profile: its
 * face is up for a profile running toward +z, toward +z for one running down.
 */
function sweepX(pts: readonly (readonly [number, number])[], len: number): THREE.BufferGeometry {
  const pos: number[] = [];
  const idx: number[] = [];
  for (const [z, y] of pts) pos.push(-len / 2, y, z, len / 2, y, z);
  for (let i = 0; i + 1 < pts.length; i++) {
    const a = 2 * i;
    idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}
/** The points of a quadratic curve in (z, y) from a to b about control c, `n` segments, both ends included. */
function quad2(a: readonly [number, number], c: readonly [number, number], b: readonly [number, number], n: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    const u = 1 - t;
    out.push([u * u * a[0] + 2 * u * t * c[0] + t * t * b[0], u * u * a[1] + 2 * u * t * c[1] + t * t * b[1]]);
  }
  return out;
}
/** A flat cap in the plane x = 0 through a closed (z, y) outline, facing +x or -x. */
function capX(pts: readonly (readonly [number, number])[], face: 1 | -1): THREE.BufferGeometry {
  const g = new THREE.ShapeGeometry(new THREE.Shape(pts.map(([z, y]) => new THREE.Vector2(z, y))));
  // The shape lies in xy facing +z; turned so its x runs along +z it faces -x.
  g.rotateY(-Math.PI / 2);
  if (face > 0) {
    const ix = g.index!.array;
    for (let i = 0; i < ix.length; i += 3) {
      const tmp = ix[i + 1]!;
      ix[i + 1] = ix[i + 2]!;
      ix[i + 2] = tmp;
    }
    const n = g.attributes.normal as THREE.BufferAttribute;
    for (let i = 0; i < n.count; i++) n.setX(i, -n.getX(i));
  }
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Reach-ins                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A top-mount reach-in on 5 in casters with one or two solid doors: the foamed stainless body, the doors on their pins
 * with bar pulls and locks, epoxy shelves on pilasters, the evaporator housing under the ceiling at the back, an LED
 * light, the data plate; on the roof the condensing unit behind its louvred grille and the controller.
 */
function reachIn(doorsN: 1 | 2, freezer: boolean): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const yb = 0.152;
    const yt = h - 0.195;
    const dy0 = 0.192;
    const dy1 = yt - 0.025;
    const pinZ = -0.065;
    const t = 0.052;
    const zf = pinZ - t / 2 - SEAL;
    const inset = 0.012;
    const doors = planDoors(w, doorsN === 1 ? ["front_door"] : ["door_left", "door_right"], 0.035, inset);
    const cav = cabBody(s, { w, d, y0: yb, y1: yt, zf, openings: doors.map((o) => openingOf(o, dy0, dy1)) });
    // Plate casters, 5 in, the front pair with toe brakes.
    for (const x of [0.065, w - 0.065])
      for (const [z, brake] of [[zf - 0.065, true], [-d + 0.065, false]] as const) s.at([x, 0, z], [0, 0, 0], () => swivelCaster(s, yb, 0.0635, brake, 0.036));
    for (const o of doors) {
      const j = s.joint(o.name, [o.px, (dy0 + dy1) / 2, pinZ]);
      cabDoor(j, { w: o.w, h: dy1 - dy0, t, side: o.side, inset, handle: "bar", hl: 0.34, hx: 0.032, lock: true });
      hingeBrackets(s, o.px, o.side, dy0, dy1, pinZ, zf);
    }
    // Inside: three shelves a section, the evaporator housing, a light a section, the data plate.
    const mid = w / 2;
    const sections: [number, number][] = doorsN === 1 ? [[cav.x0, cav.x1]] : [[cav.x0, mid - 0.012], [mid + 0.012, cav.x1]];
    const eh = 0.19;
    const ed = 0.3;
    const top = cav.y1 - eh - 0.03;
    const tiers = [0.27, 0.53, 0.79].map((k) => cav.y0 + (top - cav.y0) * k);
    if (s.detail > 0) cabShelves(s, cav, sections, tiers, top);
    else for (const [a, b] of sections) for (const y of tiers.slice(0, 2)) s.add(p.poly, new THREE.BoxGeometry(b - a - 0.02, 0.006, cav.zl - cav.zb - 0.06), [(a + b) / 2, y, (cav.zl + cav.zb) / 2]);
    s.at([cav.x0 + 0.025, cav.y1, cav.zb + ed], [0, 0, 0], () => cabEvaporator(s, cav.x1 - cav.x0 - 0.05, eh, ed, doorsN, (i) => `evap_fan_${i}`));
    for (const [a, b] of sections) {
      const x = (a + b) / 2;
      s.add(p.ss, bx(s, 0.3, 0.018, 0.034, 0.003), [x, cav.y1 - 0.009, cav.zl - 0.07]);
      if (s.fine(1)) s.add(p.ledWhite, new THREE.PlaneGeometry(0.28, 0.018).rotateX(Math.PI / 2), [x, cav.y1 - 0.0186, cav.zl - 0.07]);
    }
    s.at([cav.x0 + 0.0006, cav.y1 - 0.32, cav.zl - 0.12], [0, Math.PI / 2, 0], () => dataPlate(s, 0.08, 0.08));
    const readout = s.fine(1) ? liveReadout("#3cff7a") : null;
    topCondenser(s, w, d, yt, h, pinZ + t / 2 - 0.001, zf, 0.15, "cond_fan", () => {
      s.at([0, 0.025, 0], [0, 0, 0], () => thermometer(s, freezer ? L.frz : L.cool, 0.084, 0.046, readout?.mat ?? null));
      rocker(s, [0, -0.05, 0]);
    });
    const fans: [string, number][] = doors.map((_, i) => [`evap_fan_${i}`, 1]);
    fans.push(["cond_fan", 0.8]);
    return (b) => cabHooks(b, doors, 2.094, fans, readout, freezer ? -18 : 3);
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The undercounter worktop and the makeline                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/** A concave fillet `r` along x from 0 to -len, in the corner of a floor at y = 0 and a wall at z = 0 behind it. */
function cove(s: Shop, mat: THREE.Material, len: number, r: number): void {
  const sh = new THREE.Shape();
  sh.moveTo(0, 0);
  sh.lineTo(r, 0);
  sh.absarc(r, r, r, -Math.PI / 2, Math.PI, true);
  sh.lineTo(0, 0);
  s.add(mat, plate(s, sh, len, 0, 6), [0, 0, 0], [0, -Math.PI / 2, 0]);
}

/**
 * The two-door undercounter worktop on 4 in casters: the stainless top with its nosing and coved rear upturn, the
 * foamed body and its doors with bar pulls, a shelf a section, the evaporator under the ceiling, the doghouse over
 * the compressor, and the louvred kick panel the condensing unit breathes through, the controller beside it.
 */
const undercounter: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const yc = 0.13;
  const yk = 0.232;
  const yw = h - 0.068;
  const pinZ = -0.065;
  const t = 0.046;
  const zf = pinZ - t / 2 - SEAL;
  const inset = 0.012;
  const dy0 = yk + 0.003;
  const dy1 = yw - 0.014;
  const doors = planDoors(w, ["door_left", "door_right"], 0.03, inset);
  const cav = cabBody(s, { w, d, y0: yk, y1: yw, zf, openings: doors.map((o) => openingOf(o, dy0, dy1)) });
  for (const o of doors) {
    const j = s.joint(o.name, [o.px, (dy0 + dy1) / 2, pinZ]);
    cabDoor(j, { w: o.w, h: dy1 - dy0, t, side: o.side, inset, handle: "bar", hl: 0.24, hx: 0.037, hy: 0.507 - (dy0 + dy1) / 2 });
    hingeBrackets(s, o.px, o.side, dy0, dy1, pinZ, zf, yw);
  }
  // The worktop: a 38 mm nosing round a sound-deadened deck, the rear upturn coved into it.
  s.add(p.ss, bx(s, w, 0.038, d - 0.003, 0.006), [w / 2, yw + 0.019, -d / 2 - 0.0015]);
  s.add(p.ss, bx(s, w, h - yw - 0.036, 0.012, 0.004), [w / 2, (yw + 0.036 + h) / 2, -d + 0.006]);
  if (s.fine(1)) s.at([w, yw + 0.038, -d + 0.012], [0, 0, 0], () => cove(s, p.ss, w, 0.008));
  for (const x of [0.07, w - 0.07])
    for (const [z, brake] of [[zf - 0.06, true], [-d + 0.07, false]] as const) s.at([x, 0, z], [0, 0, 0], () => swivelCaster(s, yc, 0.05, brake, 0.03));
  const readout = s.fine(1) ? liveReadout("#3cff7a") : null;
  const cx = w - 0.22;
  kickCondenser(s, w, d, yc, yk, zf, 0.15, cx, "cond_fan", () => {
    s.at([-0.02, 0, 0], [0, 0, 0], () => thermometer(s, L.cool, 0.07, 0.04, readout?.mat ?? null));
    rocker(s, [0.05, 0, 0]);
  });
  // Inside: the doghouse over the compressor, a shelf a section, the evaporator under the ceiling.
  const dw = cav.x1 - (cx - 0.13);
  s.add(p.ssDull, openBox(dw, 0.075, 0.24, ["ny", "nz", "px"]), [cav.x1 - dw / 2, cav.y0 + 0.0375, cav.zb + 0.12]);
  const mid = w / 2;
  const sections: [number, number][] = [[cav.x0, mid - 0.012], [mid + 0.012, cav.x1]];
  const eh = 0.12;
  const ed = 0.2;
  const top = cav.y1 - eh - 0.02;
  if (s.detail > 0) cabShelves(s, cav, sections, [cav.y0 + (top - cav.y0) * 0.52], top);
  s.at([cav.x0 + 0.025, cav.y1, cav.zb + ed], [0, 0, 0], () => cabEvaporator(s, cav.x1 - cav.x0 - 0.05, eh, ed, 1, (i) => `evap_fan_${i}`));
  s.at([cav.x0 + 0.0006, cav.y1 - 0.12, cav.zl - 0.1], [0, Math.PI / 2, 0], () => dataPlate(s, 0.07, 0.07));
  return (b) => cabHooks(b, doors, 2.007, [["evap_fan_0", 1], ["cond_fan", 0.8]], readout, 3);
};

/**
 * The 72 in sandwich makeline on 5 in casters: the refrigerated base with two doors, a flush pan rail of eighteen
 * pans of toppings (third-size in front, sixth-size behind) hanging in its well, the insulated flip lid hinged on
 * the rear riser and resting on the front dam, the white HDPE cutting board along the front on its nosing, and the
 * condensing unit breathing through the kick panel.
 */
const makeline: Builder = (s, e) => {
  const p = palette();
  const { w, d } = e;
  const t1 = 0.0012;
  const yc = 0.15;
  const yk = 0.232;
  const yt = 0.911; // the deck under the board
  const yr = 0.934; // the pan rail's top
  const yw = 0.828; // the well's floor, over the cabinet's ceiling
  const yd = 0.975; // the dam and riser tops, where the lid rests
  const yn = 0.876; // the nosing's underside
  const pinZ = -0.06;
  const t = 0.046;
  const zf = pinZ - t / 2 - SEAL;
  const inset = 0.012;
  const dy0 = yk + 0.003;
  const dy1 = yw - 0.003 + (PLUG - 0.003);
  const doors = planDoors(w, ["door_left", "door_right"], 0.03, inset);
  const cav = cabBody(s, { w, d, y0: yk, y1: yt, zf, openings: doors.map((o) => openingOf(o, dy0, dy1)), top: false });
  for (const o of doors) {
    const j = s.joint(o.name, [o.px, (dy0 + dy1) / 2, pinZ]);
    cabDoor(j, { w: o.w, h: dy1 - dy0, t, side: o.side, inset, handle: "lip", hl: 0.5 });
    hingeBrackets(s, o.px, o.side, dy0, dy1, pinZ, zf, yn);
  }
  for (const x of [0.08, w - 0.08])
    for (const [z, brake] of [[zf - 0.07, true], [-d + 0.08, false]] as const) s.at([x, 0, z], [0, 0, 0], () => swivelCaster(s, yc, 0.0635, brake, 0.036));
  // The deck and nosing the board lies on, the board.
  const zr = -0.33;
  s.add(p.ss, new THREE.BoxGeometry(w, t1, -0.012 - zr), [w / 2, yt - t1 / 2, (zr - 0.012) / 2]);
  s.add(p.ss, bx(s, w, yt - yn, t1, 0.0005), [w / 2, (yt + yn) / 2, -0.0126]);
  s.add(p.ss, new THREE.BoxGeometry(w, t1, -0.012 - zf), [w / 2, yn + t1 / 2, (zf - 0.012) / 2]);
  for (const x of [t1 / 2, w - t1 / 2]) s.add(p.ss, new THREE.BoxGeometry(t1, yt - yn, -0.012 - zf), [x, (yt + yn) / 2, (zf - 0.012) / 2]);
  {
    // The cutting board: 3/4 in HDPE, its edges chamfered, lying on the deck.
    const bw = w - 0.024;
    const bd = 0.3;
    const ch = s.detail > 0 ? 0.003 : 0;
    const g = new THREE.ExtrudeGeometry(roundRect(bw - 2 * ch, bd - 2 * ch, 0.006), {
      depth: 0.019 - 2 * ch,
      bevelEnabled: ch > 0,
      bevelThickness: ch,
      bevelSize: ch,
      bevelSegments: 1,
      curveSegments: s.seg(3, 1),
    });
    g.translate(0, 0, ch);
    g.rotateX(-Math.PI / 2);
    s.add(boardMaterial(), g, [w / 2, yt + 0.0005, -0.016 - bd / 2]);
  }
  // The rail: the dam in front, the riser behind with the hinge, the ends, the strips round the well, the well.
  const n = 9;
  const pw = 0.176;
  const xw0 = (w - n * pw) / 2 - 0.003;
  const xw1 = w - xw0;
  const zw0 = -0.3372;
  const zw1 = -0.829;
  const zR = -0.8305;
  s.add(p.ss, bx(s, w, yd - yt, 0.006, 0.0025), [w / 2, (yt + yd) / 2, zr - 0.003]);
  s.add(p.ss, bx(s, w, yd - yt, d + zR, 0.004), [w / 2, (yt + yd) / 2, (zR - d) / 2]);
  for (const x of [t1 / 2, w - t1 / 2]) s.add(p.ss, new THREE.BoxGeometry(t1, yd - yt, zr - 0.006 - zR), [x, (yt + yd) / 2, (zr - 0.006 + zR) / 2]);
  for (const [a, b] of [[t1, xw0], [xw1, w - t1]] as const) s.add(p.ss, new THREE.BoxGeometry(b - a, 0.0015, zr - 0.006 - zR), [(a + b) / 2, yr - 0.00075, (zr - 0.006 + zR) / 2]);
  const wh = yr - 0.0008 - yw;
  const wd = zw0 - zw1;
  for (const z of [zw0 + t1 / 2, zw1 - t1 / 2]) s.add(p.ssDull, new THREE.BoxGeometry(xw1 - xw0 + 2 * t1, wh, t1), [w / 2, yw + wh / 2, z]);
  for (const x of [xw0 - t1 / 2, xw1 + t1 / 2]) s.add(p.ssDull, new THREE.BoxGeometry(t1, wh, wd), [x, yw + wh / 2, (zw0 + zw1) / 2]);
  s.add(p.ssDull, new THREE.BoxGeometry(xw1 - xw0, t1, wd), [w / 2, yw + t1 / 2, (zw0 + zw1) / 2]);
  // The pans: a third in front, a sixth behind, nine across, each on its flange.
  const thirds = [T.tomato, T.lettuce, T.onion, T.cheese, T.pepper, T.olive, T.jalapeno, T.cucumber, T.ham];
  const sixths = [T.chicken, T.mushroom, T.corn, T.pickle, T.banana, T.spinach, T.redOnion, T.olive, T.jalapeno];
  const ph = 0.1;
  if (s.detail === 0) s.add(p.ss, new THREE.BoxGeometry(xw1 - xw0, ph, wd), [w / 2, yr - ph / 2, (zw0 + zw1) / 2]);
  else for (let i = 0; i < n; i++) {
    const x = xw0 + 0.003 + pw / 2 + i * pw;
    s.at([x, yr - ph + 0.0015, -0.34 - 0.162], [0, 0, 0], () => hotelPan(s, pw - 0.022, 0.324 - 0.022, ph, thirds[i]!, 0.93, false, 31 + i));
    s.at([x, yr - ph + 0.0015, -0.664 - 0.081], [0, 0, 0], () => hotelPan(s, pw - 0.022, 0.162 - 0.022, ph, sixths[i]!, 0.9, false, 51 + i));
  }
  // The lid: an insulated shell swept across the rail, domed over the pans, its pull along the front edge, on three
  // knuckles of the piano hinge along the riser.
  const lid = s.joint("flip_lid", [w / 2, 0.974, -0.824]);
  {
    const lw = w - 0.024;
    const k = lid.fine(2) ? 16 : lid.fine(1) ? 6 : 2;
    const outer = quad2([-0.006, 0.03], [0.24, 0.257], [0.492, 0.032], k);
    const inner = quad2([0.468, 0.002], [0.24, 0.222], [0.014, 0.002], k);
    lid.add(p.ss, sweepX(outer, lw));
    lid.add(p.ssDull, sweepX(inner, lw));
    lid.add(p.ss, sweepX([[0.492, 0.032], [0.492, 0.002], [0.468, 0.002]], lw));
    lid.add(p.ss, sweepX([[0.014, 0.002], [-0.006, 0.002], [-0.006, 0.03]], lw));
    // Solid end plates close the dome, faced both ways so they also read from under the open lid.
    const end: [number, number][] = [[-0.006, 0.002], ...outer, [0.492, 0.002]];
    lid.add(p.ss, capX(end, -1), [-lw / 2, 0, 0]);
    lid.add(p.ss, capX(end, 1), [lw / 2, 0, 0]);
    lid.add(p.ssDull, capX(end, 1), [-lw / 2 + 0.0015, 0, 0]);
    lid.add(p.ssDull, capX(end, -1), [lw / 2 - 0.0015, 0, 0]);
    lid.add(p.ssBright, bx(lid, lw * 0.45, 0.012, 0.018, 0.004), [0, 0.027, 0.5]);
    for (const x of [-0.6, 0, 0.6]) {
      lid.add(p.ssBright, cyl(lid, 0.006, 0.09, "x", 14), [x, 0, 0]);
      lid.add(p.ssBright, new THREE.BoxGeometry(0.09, 0.0025, 0.03), [x, 0.0045, 0.012]);
    }
  }
  for (const x of [-0.65, -0.55, 0.55, 0.65, -0.05, 0.05]) s.add(p.ssBright, bx(s, 0.01, 0.006, 0.02, 0.001), [w / 2 + x, yd - 0.004, -0.824 - 0.004]);
  const readout = s.fine(1) ? liveReadout("#3cff7a") : null;
  const cx = w - 0.25;
  kickCondenser(s, w, d, yc, yk, zf, 0.15, cx, "cond_fan", () => {
    s.at([-0.02, -0.004, 0], [0, 0, 0], () => thermometer(s, L.cool, 0.066, 0.038, readout?.mat ?? null));
    rocker(s, [0.05, -0.004, 0]);
  });
  // Inside: the doghouse, a shelf a section, the evaporator at the back under the well.
  const dw = cav.x1 - (cx - 0.13);
  s.add(p.ssDull, openBox(dw, 0.07, 0.24, ["ny", "nz", "px"]), [cav.x1 - dw / 2, cav.y0 + 0.035, cav.zb + 0.12]);
  const mid = w / 2;
  const sections: [number, number][] = [[cav.x0, mid - 0.012], [mid + 0.012, cav.x1]];
  const eh = 0.12;
  const ed = 0.2;
  const ceil = cav.y1;
  const top = ceil - eh - 0.02;
  if (s.detail > 0) cabShelves(s, cav, sections, [cav.y0 + (top - cav.y0) * 0.52], top);
  s.at([cav.x0 + 0.025, ceil, cav.zb + ed], [0, 0, 0], () => cabEvaporator(s, cav.x1 - cav.x0 - 0.05, eh, ed, 1, (i) => `evap_fan_${i}`));
  s.at([cav.x0 + 0.0006, ceil - 0.12, cav.zl - 0.1], [0, Math.PI / 2, 0], () => dataPlate(s, 0.07, 0.07));
  return (b) => {
    const hk = cabHooks(b, doors, 1.92, [["evap_fan_0", 1], ["cond_fan", 0.8]], readout, 2);
    const lj = b.joints.get("flip_lid");
    const setLid = (f: number) => turnJoint(lj, "x", -Math.max(0, Math.min(1, f)) * 1.833);
    return {
      tick: hk.tick,
      open: (f) => {
        setLid(f);
        hk.open?.(f);
      },
      set: (name, v) => (name === "lid" || name === "flip_lid" ? setLid(v) : hk.set?.(name, v)),
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* The blast chiller                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A reach-in blast chiller for `levels` pans on 6 in legs: the foamed body and its heavy door, rod pan racks down
 * the sides, the evaporator behind its shroud across the back with the big fan blowing over the pans, four hotel
 * pans loaded; on the roof the control panel (readout, mode legends, keys), the core probe in its clips on a
 * coiled lead, and the condensing unit behind its louvres.
 */
function blastChiller(levels: number): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const yb = 0.152;
    const yt = h - 0.195;
    const dy0 = 0.192;
    const dy1 = yt - 0.025;
    const pinZ = -0.065;
    const t = 0.06;
    const zf = pinZ - t / 2 - SEAL;
    const inset = 0.012;
    const doors = planDoors(w, ["front_door"], 0.04, inset);
    const cav = cabBody(s, { w, d, y0: yb, y1: yt, zf, openings: doors.map((o) => openingOf(o, dy0, dy1)) });
    for (const x of [0.05, w - 0.05]) for (const z of [zf - 0.05, -d + 0.05]) leg6(s, x, z, yb);
    const o = doors[0]!;
    const j = s.joint(o.name, [o.px, (dy0 + dy1) / 2, pinZ]);
    cabDoor(j, { w: o.w, h: dy1 - dy0, t, side: 1, inset, handle: "bar", hl: 0.5, hx: 0.03, lock: true });
    hingeBrackets(s, o.px, 1, dy0, dy1, pinZ, zf);
    // The evaporator: a shroud across the back, the fan in its ring blowing forward, the coil behind.
    const ze = cav.zb + 0.17;
    const ew = cav.x1 - cav.x0;
    const eh = cav.y1 - cav.y0;
    const fr = 0.15;
    const fx = (cav.x0 + cav.x1) / 2;
    const fy = cav.y0 + eh * 0.55;
    if (s.detail > 0) {
      const sh = roundRect(ew, eh, 0.002, fx, cav.y0 + eh / 2);
      sh.holes.push(circleHoleCW(fx, fy, fr + 0.006));
      s.add(p.ss, plate(s, sh, 0.0012), [0, 0, ze - 0.0012]);
    } else s.add(p.ss, new THREE.BoxGeometry(ew, eh, 0.002), [fx, cav.y0 + eh / 2, ze - 0.001]);
    // Behind the fan its motor on three struts from the shroud, and behind the motor the coil, clear of the back.
    if (s.fine(2)) {
      s.at([fx - 0.2, fy - 0.2, ze - 0.105], [0, 0, 0], () => finCoil(s, p.alu, 0.4, 0.4, 0.06, 4, false));
      s.add(p.charcoal, cyl(s, 0.05, 0.06, "z", 24), [fx, fy, ze - 0.065]);
      s.add(p.ssBright, cyl(s, 0.008, 0.008, "z", 12), [fx, fy, ze - 0.031]);
      for (let k = 0; k < 3; k++) {
        const a = Math.PI / 2 + (k * Math.PI * 2) / 3;
        const [c, sn] = [Math.cos(a), Math.sin(a)];
        s.add(p.ssDull, bar(s, [[fx + c * 0.049, fy + sn * 0.049, ze - 0.065], [fx + c * (fr + 0.03), fy + sn * (fr + 0.03), ze - 0.065], [fx + c * (fr + 0.03), fy + sn * (fr + 0.03), ze - 0.0016]], 0.004, 0.012, 8));
      }
    } else if (s.fine(1)) s.add(p.alu, new THREE.BoxGeometry(0.4, 0.4, 0.06), [fx, fy, ze - 0.135]);
    const ej = s.joint("evap_fan", [fx, fy, ze - 0.02]);
    axialFan(ej, p.alu, p.ssDull, fr - 0.008, 5);
    s.at([fx, fy, ze + 0.002], [0, 0, 0], () => fanGuard(s, p.chrome, fr + 0.006, 0.014, 0.0013));
    // Rod racks down each side: front and back uprights, a runner each level the pans' flanges ride on.
    const pitch = Math.min(0.09, (eh - 0.12) / (levels - 1 + 0.6));
    const ys = Array.from({ length: levels }, (_, i) => cav.y0 + 0.07 + i * pitch);
    const zr0 = cav.zl - 0.03;
    const zr1 = ze + 0.03;
    const rx = 0.258;
    if (s.detail > 0)
      for (const sx of [-1, 1]) {
        const x = fx + sx * rx;
        for (const z of [zr0, zr1]) wire(s, p.ss, [x + sx * 0.008, cav.y0, z], [x + sx * 0.008, ys[ys.length - 1]! + 0.04, z], 0.005, s.fine(2) ? 10 : 5);
        for (const y of ys) {
          wire(s, p.ss, [x, y - 0.004, zr0 + 0.005], [x, y - 0.004, zr1 - 0.005], 0.004, s.fine(2) ? 8 : 4);
          wire(s, p.ss, [x, y - 0.004, zr0], [x + sx * 0.008, y - 0.004, zr0], 0.0035, 5);
          if (s.fine(2)) wire(s, p.ss, [x, y - 0.004, zr1], [x + sx * 0.008, y - 0.004, zr1], 0.0035, 5);
        }
      }
    else for (const sx of [-1, 1]) s.add(p.ss, new THREE.BoxGeometry(0.01, eh - 0.1, zr0 - zr1), [fx + sx * rx, cav.y0 + eh / 2, (zr0 + zr1) / 2]);
    const loaded: [number, number, boolean][] = [[1, T.braise, false], [4, T.mac, false], [7, T.roastVeg, true], [10, T.rice, false]];
    const ph = 0.065;
    for (const [lv, tile, lidded] of loaded) {
      const y = ys[Math.min(lv, ys.length - 1)]!;
      s.at([fx, y - ph + 0.001, (zr0 + zr1) / 2], [0, 0, 0], () => hotelPan(s, 0.53 - 0.022, 0.325 - 0.022, ph, tile, 0.75, lidded, 70 + lv));
    }
    // The roof: control panel, probe, condensing unit.
    const readout = s.fine(1) ? liveReadout("#3cff7a") : null;
    const pm = plateMaterials();
    topCondenser(s, w, d, yt, h, pinZ + t / 2 - 0.001, zf, 0.36, "cond_fan", () => {
      s.add(p.bakelite, bx(s, 0.32, 0.105, 0.004, 0.003), [0, 0.035, 0.002]);
      if (readout) s.add(readout.mat, decalFull(0.11, 0.041), [-0.085, 0.058, 0.0047]);
      else s.add(p.screen, new THREE.BoxGeometry(0.11, 0.041, 0.001), [-0.085, 0.058, 0.0045]);
      if (s.fine(1)) {
        s.add(pm.plate, decal(0.08, 0.08, P.chill, 4), [0.055, 0.035, 0.0047]);
        s.add(pm.led, decal(0.05, 0.02, L.chill, 8), [-0.085, 0.012, 0.0047]);
      }
      for (const [kx, ky] of [[0.115, 0.06], [0.145, 0.06], [0.115, 0.02], [0.145, 0.02]] as const) s.add(p.charcoal, cyl(s, 0.009, 0.004, "z", 16), [kx, ky, 0.006]);
      if (s.fine(2)) for (let i = 0; i < 5; i++) s.add(p.ledAmber, new THREE.CircleGeometry(0.0016, 8), [0.017, 0.0644 - i * 0.0125, 0.0053]);
      // The core probe in its two clips, its lead coiled back to the gland.
      const py = -0.06;
      const pz = 0.012;
      for (const cxp of [0.035, 0.135]) s.add(p.plastic, bx(s, 0.012, 0.022, 0.012, 0.003), [cxp, py, 0.006]);
      s.add(p.plastic, cyl(s, 0.0075, 0.055, "x", 16, 0.0065), [0.0425, py, pz]);
      s.add(p.ssBright, cyl(s, 0.0019, 0.1, "x", 10, 0.0007), [0.12, py, pz]);
      s.add(p.plastic, cyl(s, 0.009, 0.014, "z", 16), [-0.14, py, 0.007]);
      const x0 = -0.14;
      const x1 = 0.015;
      if (s.fine(2)) {
        const turns = 24;
        const per = 8;
        const pts: THREE.Vector3[] = [];
        for (let i = 0; i <= turns * per; i++) {
          const u = i / (turns * per);
          const a = (i / per) * Math.PI * 2;
          const cxh = x0 + 0.012 + (x1 - x0 - 0.016) * u;
          const cyh = py - 0.012 * Math.sin(Math.PI * u);
          pts.push(new THREE.Vector3(cxh, cyh + 0.0055 * Math.cos(a), pz + 0.0055 * Math.sin(a)));
        }
        s.add(p.plastic, new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), turns * per, 0.0014, 4, false));
        pipe(s, p.plastic, [[x0, py, 0.014], [x0 + 0.012, py, pz]], 0.0018, 0.005, 6);
      } else if (s.fine(1)) pipe(s, p.plastic, [[x0, py, pz], [(x0 + x1) / 2, py - 0.01, pz], [x1, py, pz]], 0.005, 0.03, 6);
    });
    return (b) => cabHooks(b, doors, 2.094, [["evap_fan", 1], ["cond_fan", 0.8]], readout, -18);
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The ice machine                                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * An air-cooled cube ice machine: the cuber head on its storage bin, the bin on 6 in legs. The bin's stainless
 * shell and white polyethylene liner round a heap of half-dice cubes; the sloped flip-up door set back under the
 * lip, its moulded pull along the bottom edge. The head's louvred front with the condenser coil and fan behind it,
 * its control panel and readout, the compressor and the exhaust grille at the back.
 */
const iceMachine: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const t = 0.0012;
  const yb = 0.152;
  const yh = 1.17;
  const zb = -0.002; // the bin's front, a step behind the head's
  for (const x of [0.045, w - 0.045]) for (const z of [-0.045, -d + 0.045]) leg6(s, x, z, yb);
  // The bin: ends, back, bottom, the front below the door, the sill, the slope round the door, the lip above.
  const yL = 0.78;
  const zL = -0.03;
  const yU = 1.11;
  const zU = -0.085;
  const alpha = Math.atan2(zL - zU, yU - yL);
  const Ls = Math.hypot(yU - yL, zU - zL);
  const zAt = (y: number) => zL + ((zU - zL) * (y - yL)) / (yU - yL);
  for (const x of [t / 2, w - t / 2]) s.add(p.ss, bx(s, t, yh - yb, d + zb, 0.0005), [x, (yb + yh) / 2, (zb - d) / 2]);
  s.add(p.alu, new THREE.BoxGeometry(w - 2 * t, yh - yb - 2 * t, t), [w / 2, (yb + yh) / 2, -d + t / 2]);
  s.add(p.ssDull, new THREE.BoxGeometry(w - 2 * t, t, d + zb - t), [w / 2, yb + t / 2, (zb - d - t) / 2]);
  s.add(p.ss, bx(s, w - 2 * t, yL - yb, t, 0.0005), [w / 2, (yb + yL) / 2, zb - t / 2]);
  s.add(p.ss, new THREE.BoxGeometry(w - 2 * t, t, zb - zL), [w / 2, yL - t / 2, (zb + zL) / 2]);
  s.add(p.ss, bx(s, w - 2 * t, yh - yU, t, 0.0005), [w / 2, (yU + yh) / 2, zb - t / 2]);
  s.add(p.ss, new THREE.BoxGeometry(w - 2 * t, t, zb - zU), [w / 2, yU + t / 2, (zb + zU) / 2]);
  const ow = 0.56;
  const o0 = 0.03;
  const o1 = Ls - 0.045;
  s.at([0, yL, zL], [-alpha, 0, 0], () => {
    if (s.detail > 0) {
      const sh = roundRect(w - 2 * t, Ls, 0.001, w / 2, Ls / 2);
      sh.holes.push(rectPath(new THREE.Path(), ow, o1 - o0, 0.012, w / 2, (o0 + o1) / 2, true));
      s.add(p.ss, plate(s, sh, t), [0, 0, -t]);
    } else s.add(p.ss, new THREE.BoxGeometry(w - 2 * t, Ls, t), [w / 2, Ls / 2, -t / 2]);
  });
  s.at([w / 2 - 0.12, (yU + yh) / 2, zb], [0, 0, 0], () => dataPlate(s, 0.07, 0.04));
  // The liner: white polyethylene, the throat from the door's opening back to it.
  const lx0 = 0.045;
  const lx1 = w - 0.045;
  const ly0 = 0.24;
  const ly1 = 1.13;
  const lz0 = -d + 0.045;
  const lz1 = -0.1;
  const yo0 = yL + o0 * Math.cos(alpha);
  const yo1 = yL + o1 * Math.cos(alpha);
  const lw = lx1 - lx0;
  const lh = ly1 - ly0;
  const ld = lz1 - lz0;
  for (const x of [lx0 - t / 2, lx1 + t / 2]) s.add(p.hdpe, new THREE.BoxGeometry(t, lh, ld), [x, ly0 + lh / 2, lz0 + ld / 2]);
  for (const y of [ly0 - t / 2, ly1 + t / 2]) s.add(p.hdpe, new THREE.BoxGeometry(lw, t, ld), [w / 2, y, lz0 + ld / 2]);
  s.add(p.hdpe, new THREE.BoxGeometry(lw, lh, t), [w / 2, ly0 + lh / 2, lz0 - t / 2]);
  if (s.detail > 0) {
    const sh = roundRect(lw, lh, 0.001, w / 2, ly0 + lh / 2);
    sh.holes.push(rectPath(new THREE.Path(), ow, yo1 - yo0, 0.01, w / 2, (yo0 + yo1) / 2, true));
    s.add(p.hdpe, plate(s, sh, t), [0, 0, lz1]);
    s.add(p.hdpe, new THREE.BoxGeometry(ow, t, zAt(yo0) - lz1), [w / 2, yo0 - t / 2, (zAt(yo0) + lz1) / 2]);
    s.add(p.hdpe, new THREE.BoxGeometry(ow, t, zAt(yo1) - lz1), [w / 2, yo1 + t / 2, (zAt(yo1) + lz1) / 2]);
    // The throat's sides stop a millimetre behind the slope's face, under the opening's edge.
    const side = new THREE.Shape([new THREE.Vector2(lz1, yo0), new THREE.Vector2(zAt(yo0) - 0.001, yo0), new THREE.Vector2(zAt(yo1) - 0.001, yo1), new THREE.Vector2(lz1, yo1)]);
    for (const x of [w / 2 - ow / 2 + t, w / 2 + ow / 2]) s.add(p.hdpe, plate(s, side, t), [x, 0, 0], [0, -Math.PI / 2, 0]);
  } else s.add(p.hdpe, new THREE.BoxGeometry(lw, lh, t), [w / 2, ly0 + lh / 2, lz1 + t / 2]);
  // The ice: three heaps (one shown, by how full), each piled under the drop at the back and slumping forward at
  // its angle of repose, half-dice cubes loose on top.
  const tops = [0.6, 0.84, 1.05];
  tops.forEach((top, k) => {
    const ij = s.joint(`part_ice_${k}`, [0, 0, 0]);
    const heap = (x: number, z: number) => Math.max(ly0 + 0.012, Math.min(ly1 - 0.08, top - 0.55 * Math.hypot((x - w / 2) * 1.25, z + d * 0.58)));
    const nx = ij.fine(2) ? 14 : ij.fine(1) ? 6 : 1;
    const g = new THREE.PlaneGeometry(lw - 0.004, ld - 0.004, nx, nx);
    g.rotateX(-Math.PI / 2);
    const pa = g.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < pa.count; i++) {
      const x = pa.getX(i) + w / 2;
      const z = pa.getZ(i) + lz0 + ld / 2;
      pa.setXYZ(i, x, heap(x, z), z);
    }
    g.computeVertexNormals();
    ij.add(p.frp, g);
    if (!ij.fine(1)) return;
    const rnd = rng(90 + k);
    const cube = new THREE.BoxGeometry(0.022, 0.012, 0.022);
    const count = ij.fine(2) ? 90 : 28;
    for (let c = 0; c < count; c++) {
      const x = lx0 + 0.02 + rnd() * (lw - 0.04);
      const z = lz0 + 0.02 + rnd() * (ld - 0.04);
      ij.add(p.glass, cube.clone(), [x, heap(x, z) + 0.004, z], [rnd() * 0.8 - 0.4, rnd() * Math.PI, rnd() * 0.8 - 0.4]);
    }
  });
  // The door: a foamed stainless panel on the slope hung from a piano hinge along its top edge, its knuckles proud
  // of its face; the liner behind; a black moulded pocket let into its face along the bottom edge for the fingers.
  const pinY = 1.08;
  const zH = 0.0225; // the hinge's knuckles stand this far off the slope, on the door's face
  const pv: V3 = [w / 2, pinY + zH * Math.sin(alpha), zAt(pinY) + zH * Math.cos(alpha)];
  const dj = s.joint("bin_flip_door", pv, [-alpha, 0, 0]);
  {
    const Ld = (pinY - (yL + 0.012)) / Math.cos(alpha);
    const dw = ow + 0.04;
    const Z = (z: number) => z - zH; // the door's depths measured off the slope
    const yTop = -0.007;
    const zf = 0.0242; // its face
    const zc = 0.0058; // its core's back
    const sk = 0.0012;
    // The pull pocket, centred where a hand (or a claw) finds it.
    const pkW = 0.3;
    const pkH = 0.032;
    const yP = (0.82 - pv[1] - Z(zf) * Math.sin(alpha)) / Math.cos(alpha);
    const bot = -Ld;
    if (dj.detail === 0) dj.add(p.ss, new THREE.BoxGeometry(dw, yTop - bot, zf - zc), [0, (yTop + bot) / 2, Z((zf + zc) / 2)]);
    else {
      // Its face: a skin with the pocket's mouth cut in it, edges eased.
      const face = roundRect(dw, yTop - bot, 0.004, 0, (yTop + bot) / 2);
      face.holes.push(rectPath(new THREE.Path(), pkW, pkH, 0.012, 0, yP, true));
      dj.add(p.ss, bevelPlate(dj, face, sk, 0.0005), [0, 0, Z(zf - sk)]);
      // The core, wrapped at its edges, built round the pocket.
      const cz = zf - sk - zc;
      const czm = Z(zc + cz / 2);
      const sideW = (dw - pkW) / 2 - 0.001;
      for (const sx of [-1, 1]) dj.add(p.ss, new THREE.BoxGeometry(sideW, yTop - bot, cz), [sx * (dw / 2 - sideW / 2), (yTop + bot) / 2, czm]);
      const iw = dw - 2 * sideW;
      const pTop = yP + pkH / 2 + 0.001;
      const pBot = yP - pkH / 2 - 0.001;
      dj.add(p.ss, new THREE.BoxGeometry(iw, yTop - pTop, cz), [0, (yTop + pTop) / 2, czm]);
      if (pBot - bot > 0.0005) dj.add(p.ss, new THREE.BoxGeometry(iw, pBot - bot, cz), [0, (pBot + bot) / 2, czm]);
      const pd = 0.0135; // the pocket's depth
      dj.add(p.ss, new THREE.BoxGeometry(iw, pTop - pBot, zf - sk - pd - 0.0004 - zc), [0, yP, Z((zc + zf - sk - pd - 0.0004) / 2)]);
      // The pocket: a moulded cup, narrowing to its floor.
      const cup = loft(dj, [
        { w: pkW + 0.002, d: pkH + 0.002, r: 0.013, y: zf - sk * 0.5 },
        { w: pkW - 0.004, d: pkH - 0.004, r: 0.011, y: zf - sk - 0.003 },
        { w: pkW - 0.012, d: pkH - 0.01, r: 0.008, y: zf - sk - pd },
      ], 0, 1, dj.fine(2) ? 3 : 1);
      cup.rotateX(Math.PI / 2);
      dj.add(p.plastic, cup, [0, yP, -zH]);
    }
    if (dj.detail > 0) dj.add(p.hdpe, bx(dj, dw - 0.06, Ld - 0.05, 0.004, 0.0015), [0, -Ld / 2, Z(0.0036)]);
    // The piano hinge: knuckles alternating between the door's leaf and the bin's.
    const kr = 0.0045;
    const nk = dj.fine(2) ? 15 : 3;
    const kl = (dw - 0.02) / nk;
    const leafT = 0.001;
    const zl = Z(zf);
    const bin = (fn: () => void) => s.at(pv, [-alpha, 0, 0], fn);
    for (let i = 0; i < nk; i++) {
      const x = -dw / 2 + 0.01 + (i + 0.5) * kl;
      if (i % 2 === 0) dj.add(p.ss, cyl(dj, kr, kl - 0.0008, "x", dj.fine(2) ? 14 : 8), [x, 0, 0]);
      else bin(() => s.add(p.ss, cyl(s, kr, kl - 0.0008, "x", s.fine(2) ? 14 : 8), [x, 0, 0]));
    }
    dj.add(p.ss, new THREE.BoxGeometry(dw - 0.02, 0.018, leafT), [0, -0.004 - 0.009, zl + leafT / 2]);
    bin(() => {
      // The bin's leaf on a stand-off rail across the slope above the door.
      s.add(p.ss, new THREE.BoxGeometry(dw - 0.02, 0.018, leafT), [0, 0.004 + 0.009, zl + leafT / 2]);
      s.add(p.ss, bx(s, dw + 0.02, 0.02, zl + zH - 0.0004, 0.002), [0, 0.005 + 0.01, (zl - zH) / 2 + 0.0002]);
      if (s.fine(2)) for (let i = 0; i < 6; i++) screwLite(s, p.ssBright, [-dw / 2 + 0.04 + (i * (dw - 0.08)) / 5, 0.014, zl + leafT], 0.0022);
    });
    if (dj.fine(2)) for (let i = 0; i < 6; i++) screwLite(dj, p.ssBright, [-dw / 2 + 0.04 + (i * (dw - 0.08)) / 5, -0.014, zl + leafT], 0.0022);
  }
  // The head: bottom pan, ends, top, the louvred intake and the condenser behind it, the panel above.
  const H = h - yh;
  s.add(p.gasket, new THREE.BoxGeometry(w - 0.004, 0.0018, d - 0.008), [w / 2, yh + 0.001, -d / 2]);
  // (The pan stops behind the kick panel, whose face it would otherwise share.)
  s.add(p.ssDull, new THREE.BoxGeometry(w - 2 * t, t, d - t), [w / 2, yh + 0.002 + t / 2, -(d + t) / 2]);
  for (const x of [t / 2, w - t / 2]) s.add(p.ss, bx(s, t, H - 0.002 - t, d, 0.0005), [x, yh + 0.002 + (H - 0.002 - t) / 2, -d / 2]);
  s.add(p.ss, bx(s, w, t, d, 0.0005), [w / 2, h - t / 2, -d / 2]);
  const gy0 = yh + 0.018;
  const gy1 = yh + 0.29;
  for (const [a, b] of [[t, 0.02], [w - 0.02, w - t]] as const) s.add(p.ss, bx(s, b - a, gy1 - yh - 0.002, t, 0.0005), [(a + b) / 2, (yh + 0.002 + gy1) / 2, -t / 2]);
  s.add(p.ss, bx(s, w - 0.04, gy0 - yh - 0.002, t, 0.0005), [w / 2, (yh + 0.002 + gy0) / 2, -t / 2]);
  grilleCondenser(s, 0.02, w - 0.02, gy0, gy1, 0, "cond_fan", 0.56);
  s.add(p.ss, bx(s, w - 2 * t, h - t - gy1, t, 0.0005), [w / 2, (gy1 + h - t) / 2, -t / 2]);
  // The control panel: a dark glass strip with the readout, the mode legend, three keys, the bin-full lamp.
  const cy = h - 0.12;
  s.add(p.bakelite, bx(s, 0.3, 0.07, 0.005, 0.003), [w / 2, cy, -0.0019]);
  const pm = plateMaterials();
  if (s.fine(1)) {
    s.add(pm.led, decal(0.07, 0.03, L.ice, 8), [w / 2 - 0.09, cy + 0.008, 0.0012]);
    // The legend band of the panel tile, cropped to its printing.
    const [u0, v0, u1, v1] = tileUV(P.ice, 4, 1, 256);
    s.add(pm.plate, intoTile(new THREE.PlaneGeometry(0.12, 0.028), [u0 + (u1 - u0) * 0.07, v0 + (v1 - v0) * 0.42, u0 + (u1 - u0) * 0.83, v0 + (v1 - v0) * 0.62]), [w / 2 + 0.05, cy + 0.012, 0.0012]);
    s.add(p.ledGreen, new THREE.CircleGeometry(0.003, 10), [w / 2 + 0.13, cy + 0.02, 0.0012]);
  }
  for (const k of [-0.03, 0, 0.03]) s.add(p.charcoal, cyl(s, 0.007, 0.004, "z", 14), [w / 2 + 0.05 + k, cy - 0.022, 0.0018]);
  // Behind: the compressor, and the exhaust louvres in the back panel the fan blows out through.
  if (s.detail > 0) {
    s.at([0.2, yh + 0.0032, -0.5], [0, 0, 0], () => compressor(s, 0.07, 0.2));
    s.at([w / 2, (gy0 + gy1) / 2, -d + t], [0, Math.PI, 0], () => louvres(s, p.ss, w - 0.08, gy1 - gy0, s.fine(2) ? 14 : 5, 0.014));
    const sh = roundRect(w - 2 * t, H - 0.004 - t, 0.001, w / 2, yh + 0.002 + (H - t) / 2);
    sh.holes.push(rectPath(new THREE.Path(), w - 0.08, gy1 - gy0, 0.002, w / 2, (gy0 + gy1) / 2, true));
    s.add(p.alu, plate(s, sh, t), [0, 0, -d]);
  } else s.add(p.alu, new THREE.BoxGeometry(w - 2 * t, H - t, t), [w / 2, yh + H / 2, -d + t / 2]);
  return (b) => {
    const door = b.joints.get("bin_flip_door");
    const heaps = tops.map((_, k) => b.joints.get(`part_ice_${k}`));
    const fan = b.joints.get("cond_fan");
    let spin = 0;
    const setDoor = (f: number) => turnJoint(door, "x", -Math.max(0, Math.min(1, f)) * 1.222);
    const setIce = (v: number) => {
      const k = v < 0.15 ? -1 : v < 0.5 ? 0 : v < 0.83 ? 1 : 2;
      heaps.forEach((g, i) => {
        if (g) g.visible = i === k;
      });
    };
    setIce(0.7);
    return {
      open: setDoor,
      set: (name, v) => {
        if (name === "door" || name === "bin_flip_door") setDoor(v);
        else if (name === "ice") setIce(v);
      },
      tick: (_t, dt, load) => {
        spin += dt * Math.PI * 2 * 5 * Math.max(0, Math.min(1, load));
        turnJoint(fan, "z", -spin);
      },
    };
  };
};

/** A plane `w` by `h` facing +z with plain 0..1 UVs (for a live canvas face). */
function decalFull(w: number, h: number): THREE.BufferGeometry {
  return new THREE.PlaneGeometry(w, h);
}

const reachIn1 = reachIn(1, false);
const reachIn2 = reachIn(2, false);
const blast15 = blastChiller(15);
register({
  "NSF-REACHIN-1D": reachIn1,
  "NSF-REACHIN-FREEZER-1DOOR": reachIn(1, true),
  "NSF-REACHIN-2D": reachIn2,
  "NSF-REACHIN-COOLER-2DOOR": reachIn2,
  "NSF-UNDERCOUNTER-2D": undercounter,
  "NSF-UNDERCOUNTER-COOLER-48IN": undercounter,
  "NSF-MAKELINE-COLD-72": makeline,
  "NSF-ICEMACHINE-AIRCOOLED": iceMachine,
  "NSF-ICE-MACHINE-CUBER-30IN": iceMachine,
  "NSF-BLAST-CHILLER-15PAN": blast15,
  "NSF-BLAST-CHILLER-10PAN": blastChiller(10),
});
