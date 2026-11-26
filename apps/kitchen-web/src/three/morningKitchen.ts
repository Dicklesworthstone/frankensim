/**
 * The morning's kitchen (chapter 07): a catering commissary's hot side, two combis and a six-burner range across the
 * aisle from a worktable, insulated carriers by the door, and the cooks who work it. `MorningKitchen` builds it and
 * acts a planned morning out on it (planMorning, sim/cateringMorning.ts):
 *   - a pot or a skillet stands on each burner for as long as the plan has a batch cooking there, its flame lit;
 *   - each combi's door is open while a cook loads or unloads it, and its fan, heat and steam run with its program;
 *   - a hotel pan stands at a table place while a cook packs it, and a cutting board when they prep;
 *   - the insulated carriers by the door fill with the hot pans as they are packed, six to a carrier, and (with
 *     `rollOut`) roll out of the door when the van leaves;
 *   - the cooks walk to each station the plan sends them to and work there for its minutes, merged into stretches
 *     of at least a second of animation: a kinematic illustration (people.ts Routines).
 * The morning figure shows one kitchen; the two-kitchen figure two, the second across town.
 */
import * as THREE from "three";
import { BURNERS, OVENS, SLOTS, TABLE_SPOTS, menuDishes, type Job, type Morning } from "../sim/cateringMorning";
import { makeTool, type ToolKind } from "./heroKitchen";
import { buildFixture, type Model, type Piece } from "./models/all";
import { Routine, makeStaff, type Cohort, type RoutineStep, type Staff, type WorkKind } from "./models/people";
import { place } from "./weekLab";

/** Animation seconds per minute of the morning; the pause on the finished morning before it starts again. */
export const SEC_PER_MIN = 0.2;
export const END_HOLD_S = 3;
/** The shortest stretch of one activity the cooks act out, in animation seconds. */
const MIN_ACT_S = 1.6;

/** The line's fronts face +z at LINE_Z: two combis and the range. The worktable faces it across the aisle. */
export const LINE_Z = -0.62;
export const OVEN = { id: "NSF-COMBI-10PAN", w: 0.847, d: 0.776, h: 1.79 };
export const OVEN_X0 = [-2.08, -1.2] as const;
export const RANGE = { id: "NSF-RANGE-GAS-6B", w: 0.914, d: 0.864, h: 1.473, x0: -0.32 };
const GRATE_Y = 0.914;
const BURNER_X = [0.18, 0.457, 0.734] as const;
const BURNER_Z = [0.22, 0.58] as const;
export const TABLE = { id: "NSF-PREPTABLE-96", w: 2.438, d: 0.762, h: 0.914, x0: -1.98 };
export const TABLE_Z = 0.8;
const POT = { id: "NSF-STOCKPOT-8QT", w: 0.24, d: 0.3, h: 0.215 };
const SKILLET = { id: "NSF-SKILLET-10IN", w: 0.254, d: 0.42, h: 0.05 };
const BOARD = { id: "NSF-CUTTING-BOARD-HDPE-24X18", w: 0.61, d: 0.457, h: 0.019 };
const HPAN = { id: "NSF-HOTEL-PAN-FULL-4IN", w: 0.53, d: 0.325, h: 0.1 };
export const CARRIER = { id: "NSF-CARRIER-INSULATED", w: 0.46, d: 0.635, h: 0.625 };
export const CARRIER_X0 = 0.7;
export const CARRIER_COLS = 3;
const CARRIER_HIGH = 3;
/** Hot pans to a carrier, as the van's carriers are loaded. */
export const PER_CARRIER = 6;
/** The pieces a kitchen builds, for main.ts to build ahead. */
export const MORNING_PIECES: readonly Piece[] = [
  [OVEN.id, OVEN, 2],
  [OVEN.id, OVEN, 2],
  [RANGE.id, RANGE, 2],
  [TABLE.id, TABLE, 2],
  ...Array.from({ length: BURNERS }, (): Piece => [POT.id, POT, 1]),
  ...Array.from({ length: BURNERS }, (): Piece => [SKILLET.id, SKILLET, 1]),
  ...Array.from({ length: TABLE_SPOTS }, (): Piece => [BOARD.id, BOARD, 1]),
  ...Array.from({ length: TABLE_SPOTS }, (): Piece => [HPAN.id, HPAN, 2]),
  [CARRIER.id, CARRIER, 1],
];

/** A cook on shift: their build and their look. */
export interface CookSpec {
  cohort: Cohort;
  seed: number;
  apron: number;
}

/** Each menu place's colour (the chicken, the starch, the third dish, the sides), on the charts and their keys. */
export const SLOT_COLOR = ["#f2b134", "#5ec8f2", "#e0905f", "#86d39b"] as const;

const burnerAt = (k: number) => ({ x: RANGE.x0 + BURNER_X[k % 3]!, z: LINE_Z - BURNER_Z[k >= 3 ? 1 : 0]! });
const tableX = (k: number) => TABLE.x0 + (TABLE.w * (k + 0.5)) / TABLE_SPOTS;
const ovenX = (o: number) => OVEN_X0[o]! + OVEN.w / 2;
export const fmtLate = (m: number) => (m < 0.5 ? "On time" : m < 59.5 ? `${Math.round(m)} min late` : `${Math.floor(Math.round(m) / 60)} h ${String(Math.round(m) % 60).padStart(2, "0")} late`);
export const shortName = (name: string) => name.replace(/^(Chicken|Mushroom|Baked|Carne asada|Roast chicken|Short rib) /, "").replace(/^./, (c) => c.toUpperCase());

/* ------------------------------------------------------------------------------------------------------------ */
/* The cooks' morning: what the plan has each doing, where, minute by minute, as a walk-and-work routine.         */
/* ------------------------------------------------------------------------------------------------------------ */

interface Spot {
  key: string;
  x: number;
  z: number;
  yaw: number;
}
interface Act {
  t0: number;
  t1: number;
  spot: Spot;
  work: WorkKind | null;
  at?: THREE.Vector3;
  surf?: number;
}

const SPOT_Z = { line: LINE_Z + 0.3, oven: LINE_Z + 0.5, table: TABLE_Z - 0.36 };

/** What cook `c` is doing over the morning, in animation seconds: one act per stretch of the same work at one place. */
function actsFor(m: Morning, c: number, endMin: number): Act[] {
  const jobs = m.jobs.filter((j) => j.cook === c);
  const cuts = new Set<number>([0, endMin]);
  for (const j of jobs) {
    cuts.add(j.t0);
    cuts.add(Math.min(endMin, j.t1));
  }
  const ts = [...cuts].filter((t) => t <= endMin).sort((a, b) => a - b);
  const acts: Act[] = [];
  let spot: Spot | null = null;
  for (let i = 0; i + 1 < ts.length; i++) {
    const a = ts[i]!;
    const b = ts[i + 1]!;
    if (b - a < 1e-6) continue;
    const mid = (a + b) / 2;
    const on = jobs.filter((j) => j.t0 <= mid && j.t1 > mid);
    const whole = on.find((j) => j.hands >= 1 - 1e-9 && j.at !== "range");
    let act: Act;
    if (whole) act = wholeAct(whole, a, b);
    else if (on.length) {
      // Minding pots and pans on the range: stood at their half of it, working the front one; with only a pot on the
      // back row (minded, not stirred, and out of easy reach), watching it.
      const ks = on.map((j) => j.spot);
      const side = on[0]!.side ?? 0;
      const x = RANGE.x0 + (RANGE.w * (side + 0.5)) / 2;
      const front = ks.find((k) => k < 3) ?? ks[0]!;
      const skillet = on.some((j) => m.burns.some((u) => u.burner === j.spot && u.t0 <= mid && u.t1 > mid && u.vessel === "skillet"));
      const v = burnerAt(front);
      act = {
        t0: a,
        t1: b,
        spot: { key: `range:${side}`, x, z: SPOT_Z.line, yaw: Math.PI },
        work: front >= 3 ? null : skillet ? "saute" : "stir",
        at: new THREE.Vector3(v.x, GRATE_Y + (skillet ? 0.03 : POT.h), v.z + (skillet ? 0 : 0.05)),
        surf: GRATE_Y,
      };
    } else act = { t0: a, t1: b, spot: spot ?? { key: `table:${c % TABLE_SPOTS}`, x: tableX(c % TABLE_SPOTS), z: SPOT_Z.table, yaw: 0 }, work: null };
    spot = act.spot;
    const prev = acts[acts.length - 1];
    if (prev && prev.spot.key === act.spot.key && prev.work === act.work) prev.t1 = b;
    else acts.push(act);
  }
  // To animation seconds, and slivers merged into the longer stretch beside them, so no motion flickers.
  for (const x of acts) {
    x.t0 *= SEC_PER_MIN;
    x.t1 *= SEC_PER_MIN;
  }
  for (let pass = 0; pass < 4; pass++)
    for (let k = 0; k < acts.length; k++) {
      const x = acts[k]!;
      if (x.t1 - x.t0 >= MIN_ACT_S || acts.length < 2) continue;
      const p = acts[k - 1];
      const n = acts[k + 1];
      const into = !p ? n! : !n ? p : p.t1 - p.t0 >= n.t1 - n.t0 ? p : n;
      into.t0 = Math.min(into.t0, x.t0);
      into.t1 = Math.max(into.t1, x.t1);
      acts.splice(k, 1);
      k--;
    }
  for (let k = acts.length - 1; k > 0; k--) {
    const p = acts[k - 1]!;
    const x = acts[k]!;
    if (p.spot.key === x.spot.key && p.work === x.work) {
      p.t1 = x.t1;
      acts.splice(k, 1);
    }
  }
  return acts;
}

function wholeAct(j: Job, a: number, b: number): Act {
  if (j.at === "table") {
    const x = tableX(j.spot);
    const pack = j.kind === "pack";
    return {
      t0: a,
      t1: b,
      spot: { key: `table:${j.spot}`, x, z: SPOT_Z.table, yaw: 0 },
      work: pack ? "plate" : "chop",
      at: new THREE.Vector3(x, TABLE.h + (pack ? 0.07 : BOARD.h), TABLE_Z + 0.27),
      surf: TABLE.h + (pack ? 0 : BOARD.h),
    };
  }
  const x = ovenX(j.spot) + 0.12;
  return {
    t0: a,
    t1: b,
    spot: { key: `oven:${j.spot}`, x, z: SPOT_Z.oven, yaw: Math.PI },
    work: "place",
    at: new THREE.Vector3(ovenX(j.spot), 1.02, LINE_Z - 0.18),
    surf: 1.0,
  };
}

export interface CookRig {
  st: Staff;
  routine: Routine | null;
  walks: Map<string, number>;
  /** The cook's own tools, made the first time they are wanted, and the one in the right hand. */
  tools: Map<ToolKind, THREE.Object3D>;
  held: ToolKind | null;
}

/** The tool each kind of work is done with. */
const TOOL_FOR: Partial<Record<WorkKind, ToolKind>> = { stir: "ladle", saute: "skillet", chop: "knife", plate: "spoon" };

/** How a cook gets from one place to another: a few short steps sideways along the same counter, else a walk that
 * comes in square to the counter (through a point a stride in front of it), so it ends facing the work. */
function moveSteps(a: Spot, b: Spot): RoutineStep[] {
  const along = Math.abs(a.yaw - b.yaw) < 1e-6 && Math.abs(a.z - b.z) < 0.15 && Math.abs(a.x - b.x) < 0.62;
  if (along) return [{ kind: "shuffle", to: [b.x, b.z] }];
  const back = b.yaw === 0 ? -0.38 : 0.38;
  return [
    { kind: "walk", path: [[b.x, b.z + back], [b.x, b.z]], speed: 1.0 },
    { kind: "turn", yaw: b.yaw },
  ];
}

/** The cook's routine for the acts, square with the morning's clock: each move starts early enough to arrive on time
 * when there is time for it, and eats into the stretch it moves to when there is not. Built twice: the first build
 * measures each move from where it really starts (its stepTimes), the second lays the stretches out with those.
 * `need`: how long the morning and the walk back to the first place take, for the loop to cover. */
function routineFor(rig: CookRig, acts: Act[], loopS: number): { routine: Routine; need: number } {
  const first = acts[0]!.spot;
  const guess = (a: Spot, b: Spot) => {
    const k = `${a.key}>${b.key}`;
    let w = rig.walks.get(k);
    if (w === undefined) {
      w = new Routine(rig.st, { x: a.x, z: a.z, yaw: a.yaw }, moveSteps(a, b), false).duration;
      rig.walks.set(k, w);
    }
    return w;
  };
  const lay = (moveS: (i: number, a: Spot, b: Spot) => number) => {
    const steps: RoutineStep[] = [];
    /** Each move's first step index and its step count. */
    const moves: [number, number][] = [];
    let cur = first;
    let t = 0;
    const go = (to: Spot, notBefore: number) => {
      const w = moveS(moves.length, cur, to);
      const start = Math.max(t, notBefore - w);
      if (start > t + 0.02) steps.push({ kind: "stand", duration: start - t });
      const ms = moveSteps(cur, to);
      moves.push([steps.length, ms.length]);
      steps.push(...ms);
      t = start + w;
      cur = to;
    };
    for (const a of acts) {
      if (a.spot.key !== cur.key) go(a.spot, a.t0);
      const d = a.t1 - t;
      if (d > 0.05) {
        steps.push(a.work ? { kind: "work", work: a.work, duration: d, opts: { at: a.at, surfaceY: a.surf } } : { kind: "stand", duration: d });
        t = a.t1;
      }
    }
    if (cur.key !== first.key) go(first, loopS - 0.3);
    const need = t + 0.3;
    steps.push({ kind: "stand", duration: Math.max(0.1, loopS - t) });
    return { steps, moves, need };
  };
  const start = { x: first.x, z: first.z, yaw: first.yaw };
  const one = lay((_, a, b) => guess(a, b));
  const r1 = new Routine(rig.st, start, one.steps, true);
  const real = one.moves.map(([i, n]) => r1.stepTimes[i + n]! - r1.stepTimes[i]!);
  const two = lay((i, a, b) => real[i] ?? guess(a, b));
  return { routine: new Routine(rig.st, start, two.steps, true), need: two.need };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The kitchen.                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** How far the carriers roll out of the door when the van leaves, metres, and in how many animation seconds. */
const ROLL_M = 3;
const ROLL_S = 2.4;

/** What the kitchen shows at a moment, for the page's readouts. */
export interface KitchenMoment {
  /** Hot pans packed so far, of all the morning's. */
  packed: number;
  hot: number;
  /** Burners with a pot or a pan cooking. */
  burning: number;
  /** A new minute of the morning since the last call (the readouts change only then). */
  minute: boolean;
}

export class MorningKitchen {
  /** Everything built, in the kitchen's own frame: the line's fronts at LINE_Z, the table's at TABLE_Z. */
  readonly group = new THREE.Group();
  readonly cooks: CookRig[];
  /** The carriers roll out of the door when the van leaves. */
  rollOut = false;
  private morning: Morning | null = null;
  /** The hot pans in the order they are packed (the sides go cold, in boxes, not in the carriers). */
  private hotPans: number[] = [];
  private lastMinute = -1;
  private packed = 0;
  private readonly ovens: Model[];
  private readonly range: Model;
  private readonly pots: Model[] = [];
  private readonly skillets: THREE.Group[] = [];
  private readonly boards: THREE.Group[] = [];
  private readonly packPans: Model[] = [];
  private readonly carrier: Model;
  private readonly carriers: THREE.Object3D[] = [];
  /** The carriers' stack, which rolls out of the door when the van leaves. */
  private readonly stack = new THREE.Group();
  private readonly door = [0, 0];
  private readonly heat = [0, 0];
  /** The burners whose skillet a cook holds: each cook holding one, the front skillet among their range jobs now. */
  private readonly inHand = new Set<number>();

  constructor(specs: readonly CookSpec[]) {
    const add = (o: THREE.Object3D) => this.group.add(o);

    /* The line: two combis and the range. */
    this.ovens = OVEN_X0.map((x0) => {
      const m = buildFixture(OVEN.id, OVEN, 2);
      add(place(m, x0, LINE_Z, true));
      return m;
    });
    this.range = buildFixture(RANGE.id, RANGE, 2);
    add(place(this.range, RANGE.x0, LINE_Z, true));
    this.range.hooks.flames?.setAll?.(0);

    /* A pot and a skillet for each burner, shown while the plan has one there. Skillets on the back burners turn their
     * handles back and out, clear of the pans in front. */
    for (let k = 0; k < BURNERS; k++) {
      const b = burnerAt(k);
      const pot = buildFixture(POT.id, POT, 1);
      pot.group.position.set(b.x - POT.w / 2, GRATE_Y + 0.0006, b.z + POT.d / 2);
      pot.group.visible = false;
      pot.hooks.set?.("level", 0.7);
      add(pot.group);
      this.pots.push(pot);
      const sk = buildFixture(SKILLET.id, SKILLET, 1);
      const pivot = new THREE.Group();
      // The fixture frame puts the pan's centre at (w/2, 0, -d + w/2), its handle towards +z.
      sk.group.position.set(-SKILLET.w / 2, 0, SKILLET.d - SKILLET.w / 2);
      pivot.add(sk.group);
      pivot.position.set(b.x, GRATE_Y + 0.0006, b.z);
      pivot.rotation.y = k < 3 ? 0 : k % 3 === 2 ? (3 * Math.PI) / 4 : (-3 * Math.PI) / 4;
      pivot.visible = false;
      add(pivot);
      this.skillets.push(pivot);
    }

    /* The worktable across the aisle, a board and a hotel pan at each of its three places. */
    const table = buildFixture(TABLE.id, TABLE, 2);
    add(place(table, TABLE.x0, TABLE_Z, false));
    for (let k = 0; k < TABLE_SPOTS; k++) {
      const x = tableX(k);
      const bd = buildFixture(BOARD.id, BOARD, 1);
      bd.group.position.set(x - BOARD.w / 2, TABLE.h + 0.0005, TABLE_Z + 0.06 + BOARD.d);
      add(bd.group);
      this.boards.push(bd.group);
      const pan = buildFixture(HPAN.id, HPAN, 2);
      pan.group.position.set(x - HPAN.w / 2, TABLE.h + 0.0005, TABLE_Z + 0.1 + HPAN.d);
      pan.hooks.set?.("fill", 0.8);
      pan.group.visible = false;
      add(pan.group);
      this.packPans.push(pan);
    }

    /* The van's carriers, stacked by the door: full ones shut, the one being loaded open. */
    add(this.stack);
    this.carrier = buildFixture(CARRIER.id, CARRIER, 1);
    this.carrier.hooks.open?.(0);
    this.carrier.hooks.set?.("pans", 0);
    for (let k = 0; k < CARRIER_COLS * CARRIER_HIGH; k++) {
      const g = this.carrier.group.clone();
      g.rotation.y = Math.PI;
      g.position.copy(this.carrierSlot(k));
      g.visible = false;
      this.stack.add(g);
      this.carriers.push(g);
    }
    this.carrier.group.rotation.y = Math.PI;
    this.carrier.group.visible = false;
    this.stack.add(this.carrier.group);

    /* The cooks. */
    this.cooks = specs.map((c, i) => {
      const st = makeStaff(c.cohort, 2, { seed: c.seed, apron: c.apron, sleeves: i === 1 ? "rolled" : "cuffed" });
      st.place(tableX(i % TABLE_SPOTS), SPOT_Z.table, 0);
      st.group.visible = false;
      add(st.group);
      return { st, routine: null, walks: new Map(), tools: new Map(), held: null };
    });
  }

  private carrierSlot(k: number): THREE.Vector3 {
    const col = Math.floor(k / CARRIER_HIGH);
    return new THREE.Vector3(CARRIER_X0 + col * (CARRIER.w + 0.03) + CARRIER.w, (k % CARRIER_HIGH) * (CARRIER.h + 0.002), TABLE_Z);
  }

  /**
   * Act out this morning with its first `crew` cooks (null: nobody on shift), over a loop of `loopS` animation
   * seconds that covers minutes 0 to `endMin`. The loop covers the morning, its closing pause and every cook's walk
   * back to where they started; if a walk back runs past it, the loop grows and the routines are laid out again (the
   * walks themselves do not change). Returns the loop's length.
   */
  plan(m: Morning | null, crew: number, endMin: number, loopS: number): number {
    this.morning = m;
    this.hotPans = [];
    if (m) for (const b of m.batches) if (b.slot < SLOTS.length) for (let p = 0; p < b.pans; p++) this.hotPans.push(b.packed);
    this.hotPans.sort((a, b) => a - b);
    this.lastMinute = -1;
    const acts = this.cooks.map((_, i) => (m && i < crew ? actsFor(m, i, endMin) : []));
    for (let pass = 0; pass < 2; pass++) {
      let need = 0;
      this.cooks.forEach((k, i) => {
        k.st.group.visible = !!m && i < crew;
        const a = acts[i]!;
        if (!a.length) {
          k.routine = null;
          return;
        }
        const r = routineFor(k, a, loopS);
        k.routine = r.routine;
        need = Math.max(need, r.need);
      });
      if (need <= loopS + 1e-6) break;
      loopS = need;
    }
    return loopS;
  }

  /** Pose each cook at loop time `at`, the tool for their work in the right hand (none walking or standing). */
  pose(at: number): void {
    for (const k of this.cooks) {
      const r = k.routine?.apply(at);
      const want = r && r.activity === "work" && r.work ? (TOOL_FOR[r.work] ?? null) : null;
      if (want === k.held) continue;
      if (k.held) k.tools.get(k.held)?.removeFromParent();
      k.held = want;
      if (!want) continue;
      let tool = k.tools.get(want);
      if (!tool) k.tools.set(want, (tool = makeTool(want)));
      k.st.hand("R").add(tool);
      tool.position.set(0, 0, 0);
      tool.quaternion.identity();
    }
  }

  /** The kitchen at minute `t` of the morning, `dt` animation seconds after the last call. */
  show(t: number, dt: number): KitchenMoment {
    const morning = this.morning;
    const { range, ovens, pots, skillets, inHand, door, heat } = this;
    inHand.clear();
    if (morning)
      this.cooks.forEach((k, c) => {
        if (k.held !== "skillet") return;
        const j = morning.jobs.find((x) => x.cook === c && x.at === "range" && x.spot < 3 && x.t0 <= t && x.t1 > t && morning.burns.some((u) => u.burner === x.spot && u.t0 <= t && u.t1 > t && u.vessel === "skillet"));
        if (j) inHand.add(j.spot);
      });
    // Burners: the vessel on each, its flame lit while it cooks.
    let burning = 0;
    for (let k = 0; k < BURNERS; k++) {
      const b = morning?.burns.find((u) => u.burner === k && u.t0 <= t && u.t1 > t);
      if (b) burning++;
      const pot = pots[k]!;
      const sk = skillets[k]!;
      pot.group.visible = !!b && b.vessel === "pot";
      // Risotto and cheese sauce are stirred: no lid. A skillet in a cook's hand is not on the grate as well.
      const lid = pot.joints.get("pot_lid_lift");
      if (lid) lid.visible = false;
      sk.visible = !!b && b.vessel === "skillet" && !inHand.has(k);
      range.hooks.flames?.setFlow(k, b ? (b.vessel === "pot" ? 0.55 : 0.8) : 0);
      if (pot.group.visible) {
        pot.hooks.set?.("boil", 0.6);
        pot.hooks.tick?.(t, dt, 0.6);
      }
    }
    range.hooks.tick?.(t, dt, 0.5);
    // Ovens: the door open to load and unload, the program running between.
    for (let o = 0; o < OVENS; o++) {
      const l = morning?.loads.find((x) => x.oven === o && x.t0 <= t && x.t3 > t);
      const open = l && (t < l.t1 || t >= l.tu) ? 1 : 0;
      const run = l && t >= l.t1 && t < l.t2 ? 1 : 0;
      door[o]! += (open - door[o]!) * Math.min(1, dt * 4);
      heat[o]! += (run - heat[o]!) * Math.min(1, dt * 2);
      const ov = ovens[o]!;
      ov.hooks.open?.(door[o]!);
      ov.hooks.set?.("heat", heat[o]!);
      ov.hooks.set?.("fan", run ? 1 : 0);
      ov.hooks.set?.("steam", run && l && /Braise|Reheat/.test(l.verb) ? 0.6 : 0.15 * run);
      ov.hooks.tick?.(t, dt, heat[o]!);
    }
    // The carriers out of the door once the van has gone.
    const gone = this.rollOut && morning ? Math.max(0, t - morning.depart) * SEC_PER_MIN : 0;
    const roll = Math.min(1, gone / ROLL_S);
    this.stack.position.x = ROLL_M * roll * roll * (3 - 2 * roll);
    this.stack.visible = roll < 1;
    const m = Math.floor(t);
    const minute = m !== this.lastMinute;
    if (minute) {
      this.lastMinute = m;
      // The table: a pan where a place is packing (its dish's food), else the board.
      const dishes = morning ? menuDishes(morning.settings.menu) : [];
      for (let k = 0; k < TABLE_SPOTS; k++) {
        const j = morning?.jobs.find((x) => x.at === "table" && x.spot === k && x.kind === "pack" && x.t0 <= t && x.t1 > t);
        this.packPans[k]!.group.visible = !!j;
        this.boards[k]!.visible = !j;
        if (j) this.packPans[k]!.hooks.set?.("food", dishes[j.slot]!.food);
      }
      // The carriers: full ones shut, the one being loaded open with its pans.
      const hot = this.hotPans;
      let packed = 0;
      while (packed < hot.length && hot[packed]! <= t) packed++;
      this.packed = packed;
      const full = Math.min(this.carriers.length, Math.floor(packed / PER_CARRIER));
      const part = packed - full * PER_CARRIER;
      this.carriers.forEach((g, k) => (g.visible = k < full));
      if (full < this.carriers.length && part > 0) {
        this.carrier.group.visible = true;
        this.carrier.group.position.copy(this.carrierSlot(full));
        this.carrier.hooks.set?.("pans", part);
        this.carrier.hooks.open?.(1);
      } else this.carrier.group.visible = false;
    }
    return { packed: this.packed, hot: this.hotPans.length, burning, minute };
  }
}
