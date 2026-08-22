// The four powered flights of December 17, 1903 as MISSION presets
// (plan §2.2 historian's loop, §3.2 validation anchors). PURE data +
// comparison math: no DOM, no sim imports, headless-tested.
//
// Honesty law (plan §2.1 copy rule): the historical flights are
// DISTRIBUTIONS, never "beat this" targets, and every verdict line
// names the comparison for what it is — a modeled reconstruction
// against period accounts, not a scoreboard or a record claim.
//
// Band units are METRES. Precision classes follow the source lineages:
//   - "measured": one authoritative surveyed figure (flight 4 —
//     Orville's diary distance reconciled with the 1904 track survey).
//   - "accounts": period accounts disagree; the band spans them.

export interface HistoricalFlight {
  readonly id: 1 | 2 | 3 | 4;
  readonly pilot: "Orville" | "Wilbur";
  /** Representative duration [s] (accounts round coarsely). */
  readonly durationS: number;
  /** Conservative band edge [m]. */
  readonly lowM: number;
  /** Generous band edge [m]. */
  readonly highM: number;
  readonly precision: "measured" | "accounts";
  /** One-line provenance phrasing (card-ready). */
  readonly note: string;
}

/** 120 ft = 36.58 m; 100 ft = 30.48 m; 175 ft = 53.34 m;
 * 200 ft = 60.96 m; 852 ft = 259.69 m. */
export const DEC17_FLIGHTS: readonly HistoricalFlight[] = [
  {
    id: 1,
    pilot: "Orville",
    durationS: 12,
    lowM: 30.5,
    highM: 36.6,
    precision: "accounts",
    note: "~12 s; 120 ft (~37 m) canonical, diary lineage ~100 ft (~30 m)",
  },
  {
    id: 2,
    pilot: "Wilbur",
    durationS: 11.5,
    lowM: 48.8,
    highM: 57.9,
    precision: "accounts",
    note: "~11-12 s; ~175 ft (~53 m) class (accounts vary)",
  },
  {
    id: 3,
    pilot: "Orville",
    durationS: 15,
    lowM: 55.8,
    highM: 66.8,
    precision: "accounts",
    note: "~15 s; ~200 ft (~61 m) class (accounts vary)",
  },
  {
    id: 4,
    pilot: "Wilbur",
    /** 852 ft surveyed to the nearest foot (~±0.3 m); the band carries
     * a ±3 ft (0.9 m) reconciliation tolerance around 259.7 m. */
    durationS: 59,
    lowM: 258.8,
    highM: 260.6,
    precision: "measured",
    note: "59 s; 852 ft (~260 m) over ground — the only precisely measured flight",
  },
];

/** Flight by 1-based id, or null outside [1, 4]. */
export function flightByIndex(id: number): HistoricalFlight | null {
  if (!Number.isInteger(id)) {
    return null;
  }
  return DEC17_FLIGHTS.find((f) => f.id === id) ?? null;
}

/** Deterministic per-mission scenario seed: each flight draws a
 * DIFFERENT (but perfectly reproducible) wind ensemble from the same
 * atmosphere machinery. Base 1903 + a prime stride keeps ids apart. */
export function flightSeed(id: number): bigint {
  return 1903n + 7919n * BigInt(id);
}

export interface MissionOutcome {
  /** Short verdict phrase (honest, non-gamey). */
  readonly verdict: string;
  /** Card-ready lines (mission header, run numbers, verdict). */
  readonly lines: readonly string[];
}

/** Compare a finished run against the flight's historical band. The
 * wording never implies the historical value was a target to beat —
 * it is context the run landed near, short of, or beyond, and the
 * trailing line always names the comparison as modeled-vs-accounts. */
export function missionOutcome(
  f: HistoricalFlight,
  downrangeM: number,
  airborneS: number,
): MissionOutcome {
  const run = `your run: ${downrangeM.toFixed(1)} m downrange, ${airborneS.toFixed(1)} s airborne`;
  const band =
    f.precision === "measured"
      ? `historical band: ${f.lowM.toFixed(1)}–${f.highM.toFixed(1)} m (surveyed)`
      : `historical band: ${f.lowM.toFixed(1)}–${f.highM.toFixed(1)} m (accounts vary)`;
  let verdict: string;
  if (downrangeM >= f.lowM && downrangeM <= f.highM) {
    verdict = "WITHIN the historical band";
  } else if (downrangeM < f.lowM) {
    verdict = "SHORT of the historical band — the drawn conditions and your inputs differ from the day";
  } else {
    verdict = "BEYOND the historical band — a modeled outcome, not a record claim";
  }
  return {
    verdict,
    lines: [
      `— MISSION: flight ${f.id} · ${f.pilot} — ${f.note}`,
      run,
      band,
      verdict,
      "(modeled reconstruction vs historical accounts — not a scoreboard)",
    ],
  };
}
