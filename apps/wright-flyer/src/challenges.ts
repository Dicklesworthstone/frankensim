// E9.2 challenge rail. Presets are typed menu selections and therefore
// travel through the same query-to-scenario path as every free launch;
// no challenge-only tuning or hidden scenario object exists.

import { menuQuery, scenarioFromQuery, type MenuSelection } from "./menu.ts";
import type { ScenarioInit } from "./sim/protocol.ts";

export type ChallengeId =
  | "assisted-first-flight"
  | "authentic-fourth-flight"
  | "modeled-fourth-flight"
  | "huffman-catapult";

export interface ChallengePreset {
  readonly id: ChallengeId;
  readonly title: string;
  readonly description: string;
  readonly selection: MenuSelection;
}

export const CHALLENGE_PRESETS: readonly ChallengePreset[] = [
  {
    id: "assisted-first-flight",
    title: "ASSISTED FIRST FLIGHT",
    description: "Fly the first Dec-17 ensemble with the bounded 30% canard assist.",
    selection: { mode: "human", site: "kdh", assist: true, flight: 1 },
  },
  {
    id: "authentic-fourth-flight",
    title: "AUTHENTIC FOURTH FLIGHT",
    description: "Raw human controls in the fourth-flight wind ensemble; no outcome is promised.",
    selection: { mode: "human", site: "kdh", assist: false, flight: 4 },
  },
  {
    id: "modeled-fourth-flight",
    title: "MODELED PILOT · FLIGHT 4",
    description: "Watch the registered modeled-pilot member in the fourth-flight ensemble.",
    selection: { mode: "historical", site: "kdh", assist: false, flight: 4 },
  },
  {
    id: "huffman-catapult",
    title: "HUFFMAN CATAPULT",
    description: "Load the 1904-05 catapult preset; its current model envelope may refuse honestly.",
    selection: { mode: "fixed", site: "huffman", assist: false },
  },
];

export function challengeById(id: string): ChallengePreset | null {
  return CHALLENGE_PRESETS.find((preset) => preset.id === id) ?? null;
}

export function challengeQuery(preset: ChallengePreset): string {
  return menuQuery(preset.selection);
}

export function challengeScenario(preset: ChallengePreset): ScenarioInit {
  return scenarioFromQuery(new URLSearchParams(challengeQuery(preset).slice(1)));
}

/**
 * Lossless configuration identity, not a released PhysicalScenarioId.
 * It names every ScenarioInit field explicitly so query/config drift is
 * observable without pretending this local descriptor is an engine hash.
 */
export function scenarioConfigIdentity(scenario: ScenarioInit): string {
  return JSON.stringify({
    schema: "wf-scenario-init-v1",
    seed: scenario.seed.toString(),
    rhoKgM3: scenario.rhoKgM3,
    headwindMps: scenario.headwindMps,
    mode: scenario.mode,
    member: scenario.member,
    railLengthM: scenario.railLengthM,
    maxTicks: scenario.maxTicks.toString(),
    assist: scenario.assist,
    catapult: scenario.catapult,
  });
}

export function challengeScenarioIdentity(preset: ChallengePreset): string {
  return scenarioConfigIdentity(challengeScenario(preset));
}
