// Final validated-receipt population (bead wf-root-guzez.9.4,
// E8.3b): wire the vv-scorecard's Wright Flyer rows into the badge
// system, with the interpolation-vs-extrapolation status LIVE.
//
// Laws:
//   - every subsystem badge reflects its ACTUAL frozen-registry
//     receipt: a row with a receipt digest surfaces Evidenced with
//     the link; NO receipt surfaces "Estimated" — never blank,
//     never green;
//   - the operating point's standing against the row's registered
//     context is labeled live: inside = interpolated, outside =
//     EXTRAPOLATED (and the word appears in the sentence — an
//     extrapolated read must say so where the user is looking).

export interface BridgeRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type BridgeResult<T> = { ok: true; value: T } | { ok: false; refusal: BridgeRefusal };

function refuse<T>(code: string, message: string, repair: string): BridgeResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** One WF scorecard row (the corpus registration projection). */
export interface WfScorecardRow {
  readonly datasetId: string;
  readonly metric: string;
  /** the binding receipt digest, when the registry carries one. */
  readonly receiptDigest?: string;
  /** registered context axis. */
  readonly contextName: string;
  readonly contextLo: number;
  readonly contextHi: number;
}

/** Row cap. */
export const MAX_ROWS = 64;

export type ReceiptStanding = "evidenced" | "estimated";
export type QueryStanding = "interpolated" | "extrapolated";

export interface SubsystemBadgeRow {
  readonly datasetId: string;
  readonly metric: string;
  readonly standing: ReceiptStanding;
  readonly query: QueryStanding;
  readonly sentence: string;
  readonly receiptLink: string | null;
}

/**
 * Bridge the scorecard rows to display rows at an operating point.
 * The point supplies each row's context coordinate by axis name;
 * missing coordinates refuse (a badge computed at nowhere is a lie).
 */
export function bridgeScorecard(
  rows: readonly WfScorecardRow[],
  point: Readonly<Record<string, number>>,
): BridgeResult<SubsystemBadgeRow[]> {
  if (rows.length === 0 || rows.length > MAX_ROWS) {
    return refuse(
      "bridge-rows-invalid",
      `${rows.length} rows outside [1, ${MAX_ROWS}]`,
      "the registry projection is bounded",
    );
  }
  const out: SubsystemBadgeRow[] = [];
  for (const row of rows) {
    if (row.datasetId.trim() === "" || !(row.contextLo <= row.contextHi)) {
      return refuse(
        "bridge-row-malformed",
        `row '${row.datasetId}'`,
        "named dataset, ordered context",
      );
    }
    const digest = row.receiptDigest;
    const digestValid = digest !== undefined && /^[0-9a-f]{64}$/.test(digest);
    if (digest !== undefined && !digestValid) {
      return refuse(
        "bridge-digest-malformed",
        `${row.datasetId}: '${digest}'`,
        "receipt digests are 64 lowercase hex chars",
      );
    }
    const x = point[row.contextName];
    if (x === undefined || !Number.isFinite(x)) {
      return refuse(
        "bridge-point-invalid",
        `missing coordinate '${row.contextName}' for ${row.datasetId}`,
        "supply the operating point for every registered axis",
      );
    }
    const query: QueryStanding =
      x >= row.contextLo && x <= row.contextHi ? "interpolated" : "extrapolated";
    const standing: ReceiptStanding = digestValid ? "evidenced" : "estimated";
    const sentence =
      standing === "evidenced"
        ? `${row.metric}: evidenced by frozen-registry receipt (${query})`
        : `${row.metric}: Estimated — no receipt in the frozen registry (${query})`;
    out.push({
      datasetId: row.datasetId,
      metric: row.metric,
      standing,
      query,
      sentence,
      receiptLink: digestValid ? `receipt:${digest}` : null,
    });
  }
  return { ok: true, value: out };
}
