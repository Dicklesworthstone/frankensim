/**
 * Every fixture, smallware and food a kitchen holds, built procedurally from its archetype's dimensions.
 *
 * `buildFixture(id, extents, detail)` returns the model in its fixture frame (see shop.ts), sized to the manifest's
 * extents, with live hooks for whatever on it moves or glows. Archetypes without a builder of their own fall back to
 * their family's generic shape, so a manifest never shows a hole.
 */
import * as THREE from "three";
import type { FlameSet } from "./flames";
import { palette } from "./materials";
import { Shop, type Built, type Detail } from "./shop";

export interface Ext {
  /** Width across the front, depth front to back, height: metres. */
  w: number;
  d: number;
  h: number;
}

/** What a built model can be told to do while it runs. */
export interface Hooks {
  /** Gas flames, one burner per valve. */
  flames?: FlameSet;
  /** Burner valve count (excluding pilots), in the order the flame set numbers them. */
  burners?: number;
  /** Materials whose glow tracks a heat (elements, radiants, embers), 0 cold to 1 bright. */
  glows?: THREE.MeshStandardMaterial[];
  /** Advance anything continuous (fans, conveyor belts, bubbles). `load` 0..1 is how hard the fixture works. */
  tick?(t: number, dt: number, load: number): void;
  /** Open the fixture's doors, lids or drawers by a fraction 0..1, where it has them. */
  open?(frac: number): void;
  /**
   * Set a named live quantity the model shows (each family file documents its names): a steak's `browning` and
   * `shrink`, a loaf's `spring` and `crust`, an oven's `fan` and `steam`, a walk-in's `door` and `frost`, a fryer
   * basket's `lift`. Values are whatever the kernel computed; the model only draws them.
   */
  set?(name: string, value: number): void;
}

export interface Model extends Built {
  id: string;
  ext: Ext;
  hooks: Hooks;
}

/** A builder adds parts to the shop for an archetype of extents `e`, and may return a function that, given the
 *  built result (its joints), returns the hooks. */
export type Builder = (s: Shop, e: Ext) => ((b: Built) => Hooks) | void;

const registry = new Map<string, Builder>();
const fallbacks: [RegExp, Builder][] = [];

/** Register builders by exact archetype id. */
export function register(builders: Record<string, Builder>): void {
  for (const [id, b] of Object.entries(builders)) registry.set(id, b);
}
/** Register a builder for every archetype id matching `re` that has none of its own. */
export function registerFallback(re: RegExp, b: Builder): void {
  fallbacks.push([re, b]);
}

export function hasBuilder(id: string): boolean {
  return registry.has(id) || fallbacks.some(([re]) => re.test(id));
}

/** Build one archetype in its fixture frame: one built ahead with the same arguments (see buildAhead), or now. */
export function buildFixture(id: string, e: Ext, detail: Detail, shadows = true): Model {
  const made = ahead.get(pieceKey(id, e, detail, shadows))?.pop();
  return made ? { ...made, ext: e } : buildNow(id, e, detail, shadows);
}

/** What buildFixture takes: one piece for a chapter to build. */
export type Piece = readonly [id: string, e: Ext, detail: Detail, shadows?: boolean];

/** Pieces built ahead, by their arguments, each for one buildFixture call to take. */
const ahead = new Map<string, Model[]>();
const pieceKey = (id: string, e: Ext, detail: Detail, shadows: boolean) => `${id}|${e.w}|${e.d}|${e.h}|${detail}|${shadows}`;

/**
 * Build a piece now, for the next buildFixture call with the same arguments to take. A chapter's pieces built ahead,
 * one to an idle moment, come off the long task that mounts it. Each is built once and taken once, so its taker owns
 * it exactly as if it had built it itself.
 */
export function buildAhead([id, e, detail, shadows = true]: Piece): void {
  const key = pieceKey(id, e, detail, shadows);
  const made = buildNow(id, e, detail, shadows);
  const list = ahead.get(key);
  if (list) list.push(made);
  else ahead.set(key, [made]);
}

function buildNow(id: string, e: Ext, detail: Detail, shadows: boolean): Model {
  const b = registry.get(id) ?? fallbacks.find(([re]) => re.test(id))?.[1] ?? plainCase;
  const shop = new Shop(detail, id);
  const after = b(shop, e);
  const built = shop.build(shadows);
  built.group.name = id;
  const hooks = after ? after(built) : {};
  return { ...built, id, ext: e, hooks };
}

/** The last resort: a stainless case of the right size, so that nothing is ever missing. */
const plainCase: Builder = (s, e) => {
  s.add(palette().ssDull, new THREE.BoxGeometry(e.w, e.h, e.d), [e.w / 2, e.h / 2, -e.d / 2]);
};
