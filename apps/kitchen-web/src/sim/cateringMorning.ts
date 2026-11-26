/**
 * A catering commissary's morning: every meal leaves on one van at 10:30, and lunch is served at the client at 11:30.
 *
 * The kitchen is a six-burner range, two 10-pan combi ovens, a worktable with room for three cooks, and two to four
 * cooks. The menu is three dishes, each a chain of stages on those stations: prep at the table, a pot or a pan on the
 * range, a program in an oven, and packing at the table into pans for the van. A stage at the table takes a whole
 * cook for its minutes; a pot or a pan on the range takes half a cook while it cooks (one cook minds two); an oven
 * takes a whole cook to load and to unload, and nobody while its program runs.
 *
 * The plan is a list schedule, stepped from event to event:
 *   - a stage starts once the stage before it is done, its station has room and a cook has the hands for it;
 *   - an oven runs one program at a time: a load is one dish, every batch of it ready for that stage, up to 10 pans;
 *     whoever is free unloads it when the program ends;
 *   - each dish's cooking starts as early as leaves it at its best at lunch (its holding time, counted back from
 *     11:30), but no more than 90 minutes before it would have to start to be packed just in time for 10:30; a dish
 *     packed late is started earlier, by as much as it was late, as far as its holding time allows. That is the plan,
 *     on time or not. When it is late, the same again with no holding limit says what making the van would cost:
 *     how many portions would be past their best;
 *   - free hands take the ready stage with the least time to spare before 10:30, at the station they are already
 *     at when they can.
 *
 * "Cooked the day before" moves the prep that keeps a day, and the whole cooking of the dishes that keep, to the
 * quiet day before: those dishes are only reheated in an oven this morning.
 *
 * `planMorning` runs one morning; `morningLimit` finds the most meals a plan gets onto the van at 10:30 with every
 * portion still at its best at lunch. Everything is deterministic.
 */

/** Minutes from 6:00: the van leaves at 10:30, lunch is served at 11:30, and the model gives up at 16:00. */
export const DOOR = 270;
export const SERVE = 330;
export const HORIZON = 600;
export const BURNERS = 6;
export const OVENS = 2;
export const OVEN_PANS = 10;
export const TABLE_SPOTS = 3;
/** Cooks who can work the 36-inch range at once, shoulder to shoulder: one at its left half, one at its right. */
export const RANGE_PLACES = 2;
/** The most minutes before just in time a dish's cooking is started, however well it holds. */
export const EARLY = 90;

export type Station = "table" | "range" | "oven";
export type Vessel = "pot" | "skillet";

export interface Stage {
  /** What the cooks are doing, for the chart. */
  verb: string;
  at: Station;
  /** Table: hands-on minutes per portion. Range: minutes a batch cooks. Oven: minutes the program runs. */
  min: number;
  /** Range: what the batch cooks in, and the share of a cook's hands it needs the whole time (0.5: one minds two). */
  vessel?: Vessel;
  hands?: number;
  /** Oven: hands-on minutes per pan, to load it, and again to unload it. */
  perPan?: number;
  /** Table prep that keeps a day in the walk-in, so it can be done the day before. */
  keeps?: boolean;
}

export interface MorningDish {
  key: string;
  name: string;
  /** Portions in a batch: a pot or a skillet on the range, or the batch the prep makes up for an oven. */
  batch: number;
  /** Portions in a 4-inch full hotel pan, packed for the van. */
  pan: number;
  /** Minutes it stays at its best once cooked, held hot in an insulated carrier. */
  holdMin: number;
  stages: readonly Stage[];
  /** A dish that keeps: its stages this morning when it was cooked the day before and chilled. */
  reheat?: readonly Stage[];
  /** The hotel pan's food (smallwares food index) and the 3D look's key. */
  food: number;
}

const pack = (min: number, verb = "Pack"): Stage => ({ verb, at: "table", min });

export const DISHES = {
  piccata: {
    key: "piccata",
    name: "Chicken piccata",
    batch: 8,
    pan: 24,
    holdMin: 90,
    stages: [{ verb: "Pound and dredge", at: "table", min: 0.15, keeps: true }, { verb: "Sear", at: "range", min: 8, vessel: "skillet", hands: 0.5 }, pack(0.04, "Shingle and sauce")],
    food: 0,
  },
  thighs: {
    key: "thighs",
    name: "Roast chicken thighs",
    batch: 80,
    pan: 20,
    holdMin: 150,
    stages: [{ verb: "Trim and season", at: "table", min: 0.06, keeps: true }, { verb: "Roast", at: "oven", min: 35, perPan: 0.4 }, pack(0.01, "Lid and label")],
    food: 0,
  },
  risotto: {
    key: "risotto",
    name: "Mushroom risotto",
    batch: 30,
    pan: 30,
    holdMin: 75,
    stages: [{ verb: "Dice and slice", at: "table", min: 0.05, keeps: true }, { verb: "Stir", at: "range", min: 25, vessel: "pot", hands: 1 }, pack(0.02, "Pan up")],
    food: 1,
  },
  mac: {
    key: "mac",
    name: "Baked mac and cheese",
    batch: 60,
    pan: 30,
    holdMin: 180,
    stages: [
      { verb: "Grate the cheese", at: "table", min: 0.04, keeps: true },
      { verb: "Cheese sauce", at: "range", min: 20, vessel: "pot", hands: 0.5 },
      { verb: "Fold and pan", at: "table", min: 0.04 },
      { verb: "Bake", at: "oven", min: 25, perPan: 0.3 },
      pack(0.01, "Lid and label"),
    ],
    reheat: [{ verb: "Reheat", at: "oven", min: 30, perPan: 0.3 }, pack(0.01, "Lid and label")],
    food: 1,
  },
  tacos: {
    key: "tacos",
    name: "Carne asada tacos",
    batch: 20,
    pan: 20,
    holdMin: 80,
    stages: [
      { verb: "Slice, make salsa", at: "table", min: 0.12, keeps: true },
      { verb: "Sear the steak", at: "range", min: 8, vessel: "skillet", hands: 0.5 },
      { verb: "Warm tortillas", at: "range", min: 6, vessel: "skillet", hands: 1 },
      pack(0.08, "Build and wrap"),
    ],
    food: 0,
  },
  braise: {
    key: "braise",
    name: "Short rib braise",
    batch: 100,
    pan: 20,
    holdMin: 240,
    stages: [{ verb: "Trim and season", at: "table", min: 0.06, keeps: true }, { verb: "Braise", at: "oven", min: 150, perPan: 0.5 }, pack(0.03, "Portion and sauce")],
    reheat: [{ verb: "Reheat", at: "oven", min: 30, perPan: 0.3 }, pack(0.03, "Portion and sauce")],
    food: 0,
  },
} as const satisfies Record<string, MorningDish>;

/** What every meal gets whatever the menu: a salad and bread, cut and boxed at the table, kept cold. */
export const SIDES: MorningDish = {
  key: "sides",
  name: "Salad and bread",
  batch: 50,
  pan: 25,
  holdMin: Infinity,
  stages: [{ verb: "Wash and cut", at: "table", min: 0.1, keeps: true }, pack(0.12, "Box the sides")],
  food: 2,
};

/** The menu's three places, the share of the meals each takes, the dish as written and the one that does its job. */
export const SLOTS: readonly { share: number; written: MorningDish; swap: MorningDish }[] = [
  { share: 0.4, written: DISHES.piccata, swap: DISHES.thighs },
  { share: 0.35, written: DISHES.risotto, swap: DISHES.mac },
  { share: 0.25, written: DISHES.tacos, swap: DISHES.braise },
];
/** Menus by number: bit k swaps place k. */
export const MENUS = 1 << SLOTS.length;
/** The morning's dishes: the menu's three, then the sides (place 3). */
export const menuDishes = (menu: number): MorningDish[] => [...SLOTS.map((s, k) => (menu & (1 << k) ? s.swap : s.written)), SIDES];

export type Crew = 2 | 3 | 4;
export interface MorningSettings {
  meals: number;
  cooks: Crew;
  menu: number;
  /** Prep that keeps, and the dishes that keep, done the day before. */
  ahead: boolean;
}

/** One stretch of a cook's work. */
export interface Job {
  cook: number;
  /** The batch's number in the morning. */
  batch: number;
  at: Station;
  /** The burner, the oven or the table place it is at; on the range, the half of it the cook stands at. */
  spot: number;
  side?: number;
  /** The menu place (0..2) whose dish it is. */
  slot: number;
  verb: string;
  kind: "prep" | "range" | "load" | "unload" | "pack";
  t0: number;
  t1: number;
  hands: number;
}
/** A pot or a pan on a burner. */
export interface Burn {
  burner: number;
  batch: number;
  slot: number;
  vessel: Vessel;
  t0: number;
  t1: number;
}
/** An oven's load: loading from t0 to t1, its program to t2, unloading to t3. */
export interface Load {
  oven: number;
  batches: number[];
  slot: number;
  verb: string;
  pans: number;
  t0: number;
  t1: number;
  t2: number;
  t3: number;
}
export interface BatchOut {
  slot: number;
  portions: number;
  pans: number;
  /** When its cooking finished, and when it was packed. `pans`: the hotel pans its packing fills (a pan part filled
   * by one batch and topped up by the next counts at the batch that starts it). */
  cooked: number;
  packed: number;
}

export interface Morning {
  settings: MorningSettings;
  /** Per place, how much earlier than its holding time allows its cooking may start. */
  shift: number[];
  jobs: Job[];
  burns: Burn[];
  loads: Load[];
  batches: BatchOut[];
  /** The van leaves when the last pan is packed, at 10:30 at the soonest. */
  depart: number;
  late: number;
  /** Everything packed: false only when the horizon ran out first. */
  finished: boolean;
  /** Portions held longer than they stay at their best, by lunch: per menu place, and in all. */
  pastBest: number;
  pastBestBySlot: number[];
  /** When each place's last pan was packed. */
  lastPacked: number[];
  /** Cooked earlier to make the van at 10:30: the portions past their best by lunch, or null if no start makes it. */
  onTimePastBest: number | null;
  pans: number;
  /** Crew-minutes of hands-on work this morning, and moved to the day before. */
  handsMin: number;
  aheadMin: number;
  /** The station that was fullest over the last 90 minutes of work: range burners, oven time, hands, table places. */
  bottleneck: Station | "hands";
  use: Record<Station | "hands", number>;
}

/* ---------------------------------------------------------------------------------------------------------------- */

interface Batch {
  id: number;
  slot: number;
  i: number;
  portions: number;
  pans: number;
  stages: readonly Stage[];
  next: number;
  busy: boolean;
  ready: number;
  /** The earliest the first stage off the table may start (its holding time counted back from lunch). */
  release: number;
  firstCook: number;
  cooked: number;
  packed: number;
}

const EPS = 1e-6;
const stageMin = (s: Stage, b: { portions: number; pans: number }) =>
  s.at === "table" ? s.min * b.portions : s.at === "range" ? s.min : s.min + 2 * b.pans * (s.perPan ?? 0);
/** A stage's hands-on crew-minutes for a batch. */
const handsOf = (s: Stage, b: { portions: number; pans: number }) =>
  s.at === "table" ? s.min * b.portions : s.at === "range" ? s.min * (s.hands ?? 1) : 2 * b.pans * (s.perPan ?? 0);

/** The batches a dish is cooked in for these portions: as few as its batch size allows, the portions shared evenly. */
function batchesFor(d: MorningDish, slot: number, portions: number, ahead: boolean): Batch[] {
  if (portions <= 0) return [];
  const stages = ahead && d.reheat ? d.reheat : ahead ? d.stages.filter((s) => !s.keeps) : d.stages;
  const n = Math.ceil(portions / d.batch - EPS);
  const out: Batch[] = [];
  let left = portions;
  for (let i = 0; i < n; i++) {
    const p = Math.round(left / (n - i));
    left -= p;
    const firstCook = stages.findIndex((s) => s.at !== "table");
    out.push({ id: 0, slot, i, portions: p, pans: Math.ceil(p / d.pan - EPS), stages, next: 0, busy: false, ready: 0, release: 0, firstCook, cooked: -1, packed: -1 });
  }
  return out;
}

/** What each place's dish takes of the meals (shares rounded so they sum to the meals), then the sides: every meal. */
export function portionsFor(meals: number): number[] {
  const out = SLOTS.map((s) => Math.floor(meals * s.share));
  let left = meals - out.reduce((a, b) => a + b, 0);
  for (let k = 0; left > 0; k = (k + 1) % out.length, left--) out[k]!++;
  return [...out, meals];
}

function schedule(s: MorningSettings, shift: readonly number[], fresh: boolean): Morning {
  const dishes = menuDishes(s.menu);
  const portions = portionsFor(s.meals);
  const batches = dishes.flatMap((d, k) => batchesFor(d, k, portions[k]!, s.ahead));
  batches.forEach((b, i) => (b.id = i));
  // The work moved to the day before: the prep that keeps, or all but the packing of a dish that keeps.
  let aheadMin = 0;
  if (s.ahead)
    dishes.forEach((d, k) => {
      const moved = d.reheat ? d.stages.slice(0, -1) : d.stages.filter((st) => st.keeps);
      for (const b of batchesFor(d, k, portions[k]!, false)) for (const st of moved) aheadMin += handsOf(st, b);
    });
  // Each batch's cooking starts as early as lets it be cooked (out of its pot, its oven program done) just as its
  // holding time runs back from lunch, but no more than EARLY before just in time for the van (its stages from the
  // first off the table to the end); its shift moves that earlier, only as far as its holding time allows if `fresh`.
  for (const b of batches) {
    if (b.firstCook < 0) continue;
    const d = dishes[b.slot]!;
    let last = b.firstCook;
    let tailMin = 0;
    for (let j = b.firstCook; j < b.stages.length; j++) {
      if (b.stages[j]!.at !== "table") last = j;
      tailMin += stageMin(b.stages[j]!, b);
    }
    let chain = 0;
    for (let j = b.firstCook; j <= last; j++) {
      const st = b.stages[j]!;
      chain += st.at === "oven" ? st.min + b.pans * (st.perPan ?? 0) : stageMin(st, b);
    }
    const floor = SERVE - d.holdMin - chain;
    const soon = DOOR - tailMin - EARLY;
    const sh = shift[b.slot] ?? 0;
    b.release = Math.max(0, fresh ? Math.max(floor, soon - sh) : Math.max(floor, soon) - sh);
  }
  const tail = (b: Batch) => {
    let m = 0;
    for (let j = b.next; j < b.stages.length; j++) m += stageMin(b.stages[j]!, b);
    return m;
  };

  const spare = new Array<number>(s.cooks).fill(1);
  /** Where each cook last worked: the station, and on the range the burner column (0..2). */
  const where: { at: Station | null; col: number; spot: number }[] = Array.from({ length: s.cooks }, (_, c) => ({ at: null, col: c % 3, spot: c % TABLE_SPOTS }));
  const burner: (Batch | null)[] = new Array(BURNERS).fill(null);
  /** The cook at each half of the range (-1: free), and the pans each cook has on it. */
  const place: number[] = new Array(RANGE_PLACES).fill(-1);
  const onRange = new Array<number>(s.cooks).fill(0);
  const table: boolean[] = new Array(TABLE_SPOTS).fill(false);
  type Oven = { phase: "idle" | "load" | "run" | "wait" | "unload"; load: Load | null; batches: Batch[]; stage: Stage | null; until: number };
  const ovens: Oven[] = Array.from({ length: OVENS }, () => ({ phase: "idle", load: null, batches: [], stage: null, until: 0 }));
  type Run = { end: number; done: () => void };
  const running: Run[] = [];
  const jobs: Job[] = [];
  const burns: Burn[] = [];
  const loads: Load[] = [];
  let handsMin = 0;

  const pickCook = (need: number, at: Station, prefer: (c: number) => number): number => {
    let best = -1;
    let bestScore = Infinity;
    for (let c = 0; c < s.cooks; c++) {
      if (spare[c]! < need - EPS) continue;
      // At the station already first, then (on the range) a cook minding a pot already, then the lowest number.
      const score = (where[c]!.at === at ? 0 : 10) + prefer(c) + c * 0.01;
      if (score < bestScore) {
        bestScore = score;
        best = c;
      }
    }
    return best;
  };
  const advance = (b: Batch, t: number) => {
    const st = b.stages[b.next]!;
    if (st.at !== "table" && b.stages.slice(b.next + 1).every((x) => x.at === "table")) b.cooked = t;
    b.next++;
    b.busy = false;
    b.ready = t;
    if (b.next >= b.stages.length) b.packed = t;
  };

  const tryStart = (b: Batch, t: number): boolean => {
    const st = b.stages[b.next]!;
    const d = stageMin(st, b);
    if (st.at === "table") {
      const spot = table.indexOf(false);
      if (spot < 0) return false;
      const c = pickCook(1, "table", () => 0);
      if (c < 0) return false;
      // The place this cook last had, when it is free.
      const mine = !table[where[c]!.spot] ? where[c]!.spot : spot;
      table[mine] = true;
      spare[c] = 0;
      where[c] = { at: "table", col: where[c]!.col, spot: mine };
      b.busy = true;
      const kind = b.next === b.stages.length - 1 ? "pack" : "prep";
      jobs.push({ cook: c, batch: b.id, at: "table", spot: mine, slot: b.slot, verb: st.verb, kind, t0: t, t1: t + d, hands: 1 });
      handsMin += d;
      running.push({
        end: t + d,
        done: () => {
          table[mine] = false;
          spare[c] = 1;
          advance(b, t + d);
        },
      });
      return true;
    }
    if (st.at === "range") {
      const hands = st.hands ?? 1;
      if (!burner.includes(null)) return false;
      // Only two cooks fit at the range: one there already, minding a pan, takes the second beside it; a new one
      // comes only to a free half.
      const free = place.indexOf(-1);
      const c = pickCook(hands, "range", (k) => (place.includes(k) ? -5 : free < 0 ? Infinity : 0));
      if (c < 0 || (!place.includes(c) && free < 0)) return false;
      if (!place.includes(c)) place[free] = c;
      const side = place.indexOf(c);
      let bi = -1;
      let bd = Infinity;
      for (let k = 0; k < BURNERS; k++) {
        if (burner[k]) continue;
        // On the cook's half (the middle column between them); what is worked all the time (a skillet, a pot stirred
        // throughout) to the front row, in reach, and a pot only minded to the back.
        const col = k % 3;
        const back = k >= 3;
        const front = (st.vessel ?? "pot") === "skillet" || hands >= 1;
        const dist = (side === 0 ? col : 2 - col) + (front ? (back ? 0.6 : 0) : back ? 0 : 0.3);
        if (dist < bd) {
          bd = dist;
          bi = k;
        }
      }
      burner[bi] = b;
      spare[c] -= hands;
      onRange[c]!++;
      where[c] = { at: "range", col: side, spot: where[c]!.spot };
      b.busy = true;
      jobs.push({ cook: c, batch: b.id, at: "range", spot: bi, side, slot: b.slot, verb: st.verb, kind: "range", t0: t, t1: t + d, hands });
      burns.push({ burner: bi, batch: b.id, slot: b.slot, vessel: st.vessel ?? "pot", t0: t, t1: t + d });
      handsMin += d * hands;
      running.push({
        end: t + d,
        done: () => {
          burner[bi] = null;
          spare[c] += hands;
          if (--onRange[c]! === 0) place[place.indexOf(c)] = -1;
          advance(b, t + d);
        },
      });
      return true;
    }
    // An oven: the first idle one; the load is every batch of this dish ready for this stage, up to its pans.
    const o = ovens.findIndex((x) => x.phase === "idle");
    if (o < 0) return false;
    const c = pickCook(1, "oven", (k) => (where[k]!.at === "oven" && where[k]!.spot === o ? -1 : 0));
    if (c < 0) return false;
    const group: Batch[] = [b];
    let pans = b.pans;
    for (const x of batches) {
      if (x === b || x.busy || x.slot !== b.slot || x.next !== b.next || x.ready > t + EPS || x.next >= x.stages.length) continue;
      if (x.next === x.firstCook && x.release > t + EPS) continue;
      if (pans + x.pans > OVEN_PANS) continue;
      group.push(x);
      pans += x.pans;
    }
    const per = st.perPan ?? 0;
    const tl = pans * per;
    const ov = ovens[o]!;
    const load: Load = { oven: o, batches: group.map((x) => x.id), slot: b.slot, verb: st.verb, pans, t0: t, t1: t + tl, t2: t + tl + st.min, t3: t + tl + st.min };
    loads.push(load);
    ov.phase = "load";
    ov.load = load;
    ov.batches = group;
    ov.stage = st;
    for (const x of group) x.busy = true;
    spare[c] = 0;
    where[c] = { at: "oven", col: where[c]!.col, spot: o };
    jobs.push({ cook: c, batch: b.id, at: "oven", spot: o, slot: b.slot, verb: `Load: ${st.verb.toLowerCase()}`, kind: "load", t0: t, t1: t + tl, hands: 1 });
    handsMin += tl;
    running.push({
      end: t + tl,
      done: () => {
        spare[c] = 1;
        ov.phase = "run";
        running.push({
          end: load.t2,
          done: () => {
            ov.phase = "wait";
            for (const x of group) if (x.stages.slice(x.next + 1).every((y) => y.at === "table")) x.cooked = load.t2;
          },
        });
      },
    });
    return true;
  };

  const tryUnload = (o: number, t: number): boolean => {
    const ov = ovens[o]!;
    const c = pickCook(1, "oven", (k) => (where[k]!.at === "oven" && where[k]!.spot === o ? -1 : 0));
    if (c < 0) return false;
    const load = ov.load!;
    const tl = load.pans * (ov.stage!.perPan ?? 0);
    ov.phase = "unload";
    spare[c] = 0;
    where[c] = { at: "oven", col: where[c]!.col, spot: o };
    jobs.push({ cook: c, batch: load.batches[0]!, at: "oven", spot: o, slot: load.slot, verb: "Unload", kind: "unload", t0: t, t1: t + tl, hands: 1 });
    handsMin += tl;
    running.push({
      end: t + tl,
      done: () => {
        spare[c] = 1;
        load.t3 = t + tl;
        const group = ov.batches;
        ov.phase = "idle";
        ov.load = null;
        ov.batches = [];
        for (const x of group) {
          const cooked = x.cooked;
          advance(x, t + tl);
          x.cooked = cooked >= 0 ? cooked : x.cooked;
        }
      },
    });
    return true;
  };

  let t = 0;
  let finished = false;
  for (let guard = 0; guard < 200000; guard++) {
    // Finish what ends now, in the order it ends.
    running.sort((a, b) => a.end - b.end);
    while (running.length && running[0]!.end <= t + EPS) running.shift()!.done();
    // Start what can: ovens waiting to be unloaded first, then the ready stages, least time to spare first.
    for (let again = true; again; ) {
      again = false;
      for (let o = 0; o < OVENS; o++) if (ovens[o]!.phase === "wait" && tryUnload(o, t)) again = true;
      const ready = batches.filter((b) => !b.busy && b.next < b.stages.length && b.ready <= t + EPS && !(b.next === b.firstCook && b.release > t + EPS));
      // Food just off the heat is panned up first, then the stage with the most work still behind it.
      const hot = (b: Batch) => (b.next > 0 && b.stages[b.next - 1]!.at !== "table" ? 1 : 0);
      ready.sort((a, b) => hot(b) - hot(a) || tail(b) - tail(a) || a.slot - b.slot || a.i - b.i);
      for (const b of ready)
        if (tryStart(b, t)) {
          again = true;
          break;
        }
    }
    if (batches.every((b) => b.next >= b.stages.length)) {
      finished = true;
      break;
    }
    // The next moment anything can change: a stage ending, or a release coming due.
    let next = Infinity;
    for (const r of running) next = Math.min(next, r.end);
    for (const b of batches) if (!b.busy && b.next === b.firstCook && b.release > t + EPS) next = Math.min(next, b.release);
    if (!Number.isFinite(next) || next > HORIZON) break;
    t = Math.max(t, next);
  }

  const last = batches.reduce((m, b) => Math.max(m, b.packed), 0);
  // The hotel pans each batch fills when packed: a dish's portions go into whole pans in the order they are packed,
  // so three skillets of cutlets fill one pan between them.
  const packPans = new Array<number>(batches.length).fill(0);
  dishes.forEach((d, k) => {
    let cum = 0;
    const mine = batches.filter((b) => b.slot === k).sort((a, b) => (a.packed < 0 ? Infinity : a.packed) - (b.packed < 0 ? Infinity : b.packed) || a.i - b.i);
    for (const b of mine) {
      const before = Math.ceil(cum / d.pan - EPS);
      cum += b.portions;
      packPans[b.id] = Math.ceil(cum / d.pan - EPS) - before;
    }
  });
  const depart = finished ? Math.max(DOOR, last) : HORIZON;
  const late = depart - DOOR;
  const lunch = SERVE + late;
  const pastBestBySlot = dishes.map(() => 0);
  const lastPacked = dishes.map(() => 0);
  for (const b of batches) lastPacked[b.slot] = Math.max(lastPacked[b.slot]!, b.packed < 0 ? HORIZON : b.packed);
  for (const b of batches) if (b.cooked >= 0 && lunch - b.cooked > dishes[b.slot]!.holdMin + EPS) pastBestBySlot[b.slot]! += b.portions;

  // How full each station was over the last 90 minutes of work.
  const w1 = Math.max(60, last);
  const w0 = w1 - 90;
  const over = (a: number, b: number) => Math.max(0, Math.min(b, w1) - Math.max(a, w0));
  let burnerMin = 0;
  for (const x of burns) burnerMin += over(x.t0, x.t1);
  let ovenMin = 0;
  for (const x of loads) ovenMin += over(x.t0, x.t3);
  let handMin = 0;
  let tableMin = 0;
  for (const j of jobs) {
    handMin += over(j.t0, j.t1) * j.hands;
    if (j.at === "table") tableMin += over(j.t0, j.t1);
  }
  const use = { range: burnerMin / (BURNERS * 90), oven: ovenMin / (OVENS * 90), hands: handMin / (s.cooks * 90), table: tableMin / (TABLE_SPOTS * 90) };
  const bottleneck = (Object.keys(use) as (keyof typeof use)[]).reduce((a, b) => (use[b] > use[a] + EPS ? b : a), "range" as keyof typeof use);

  return {
    settings: s,
    shift: [...shift],
    jobs,
    burns,
    loads,
    batches: batches.map((b, i) => ({ slot: b.slot, portions: b.portions, pans: packPans[i]!, cooked: b.cooked, packed: b.packed })),
    depart,
    late,
    finished,
    pastBest: pastBestBySlot.reduce((a, b) => a + b, 0),
    pastBestBySlot,
    lastPacked,
    onTimePastBest: late <= 0.5 ? pastBestBySlot.reduce((a, b) => a + b, 0) : null,
    pans: packPans.reduce((a, b) => a + b, 0),
    handsMin,
    aheadMin,
    bottleneck,
    use,
  };
}

/** Run the morning, moving each dish that is packed late earlier by its lateness, until the van is made or nothing moves. */
function settle(s: MorningSettings, fresh: boolean, from: readonly number[] | null = null): Morning {
  let shift = from ? [...from] : new Array<number>(SLOTS.length + 1).fill(0);
  let cur = schedule(s, shift, fresh);
  let best = cur;
  for (let k = 0; k < 12 && cur.late > 0.5; k++) {
    const next = shift.map((x, d) => Math.min(SERVE, x + Math.max(0, cur.lastPacked[d]! - DOOR)));
    if (next.every((x, d) => Math.abs(x - shift[d]!) < 0.5)) break;
    shift = next;
    cur = schedule(s, shift, fresh);
    if (cur.late < best.late - 0.5 || (Math.abs(cur.late - best.late) <= 0.5 && cur.pastBest < best.pastBest)) best = cur;
  }
  return best;
}

/**
 * The morning for these settings: each dish cooked just in time for the van, moved earlier where that is late, but
 * never so early that it is past its best at lunch. When that misses the van, the same without the limit says what
 * making the van would cost.
 */
export function planMorning(s: MorningSettings): Morning {
  const plan = settle(s, true);
  if (plan.late <= 0.5) return plan;
  const rush = settle(s, false, plan.shift);
  return { ...plan, onTimePastBest: rush.late <= 0.5 ? rush.pastBest : null };
}

/** Whether a morning makes the van with everything still at its best at lunch. */
export const madeIt = (m: Morning) => m.finished && m.late <= 0.5 && m.pastBest === 0;

/** The most meals, in tens, that these settings get onto the van at 10:30 with every portion at its best at lunch. */
export function morningLimit(s: Omit<MorningSettings, "meals">, lo = 10, hi = 2400): number {
  if (!madeIt(planMorning({ ...s, meals: lo }))) return 0;
  // Bisection on tens of meals: lo makes it, hi does not.
  let a = lo / 10;
  let b = hi / 10;
  if (madeIt(planMorning({ ...s, meals: hi }))) return hi;
  while (b - a > 1) {
    const m = Math.floor((a + b) / 2);
    if (madeIt(planMorning({ ...s, meals: m * 10 }))) a = m;
    else b = m;
  }
  return a * 10;
}

/** Clock time for minutes from 6:00. */
export const clockOf = (m: number): string => {
  const t = Math.max(0, Math.round(m));
  const h = 6 + Math.floor(t / 60);
  return `${h > 12 ? h - 12 : h}:${String(t % 60).padStart(2, "0")}`;
};
