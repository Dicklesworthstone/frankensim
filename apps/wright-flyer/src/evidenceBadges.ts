// Evidence badges (bead wf-root-guzez.11.6.2, E10.2-ii): the beta
// build's live V/H badges, built VERBATIM from real receipt payloads.
//
// Laws:
//   - a badge NEVER renders a pass state without a receipt digest
//     (green without evidence is the exact failure mode the evidence
//     program exists to prevent);
//   - unknown verdicts render NO-DATA, never blank and never green;
//   - badges are passthrough: the app re-scores NOTHING.

export interface BadgeRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type BadgeResult<T> = { ok: true; value: T } | { ok: false; refusal: BadgeRefusal };

function refuse<T>(code: string, message: string, repair: string): BadgeResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** A receipt row as the evidence plane publishes it. */
export interface ReceiptRow {
  /** V/H case label (e.g. "V-08b1", "H-07"). */
  readonly caseId: string;
  /** Verdict string from the receipt ("pass", "reported-only", …). */
  readonly verdict: string;
  /** The receipt's content digest (64 hex chars) — or absent. */
  readonly receiptDigest?: string;
  /** Declared comparison class. */
  readonly comparisonClass: string;
}

export type BadgeState = "evidenced-pass" | "reported" | "no-data";

export interface EvidenceBadge {
  readonly caseId: string;
  readonly state: BadgeState;
  readonly receiptDigest: string | null;
  readonly comparisonClass: string;
  /** the exact verdict string, verbatim. */
  readonly verdict: string;
}

/** Badge budget per view. */
export const MAX_BADGES = 64;

/**
 * Build badges from receipt rows. The digest-required law lives
 * here: a "pass" verdict WITHOUT a digest refuses — it never demotes
 * silently, because a silently demoted forgery would still teach the
 * user something false about the evidence plane.
 */
export function buildEvidenceBadges(rows: readonly ReceiptRow[]): BadgeResult<EvidenceBadge[]> {
  if (rows.length > MAX_BADGES) {
    return refuse(
      "badge-count-exceeded",
      `${rows.length} rows > ${MAX_BADGES}`,
      "page the evidence registry",
    );
  }
  const badges: EvidenceBadge[] = [];
  for (const row of rows) {
    if (row.caseId.trim() === "") {
      return refuse("badge-case-missing", "a row without a case id", "name the V/H case");
    }
    const digest = row.receiptDigest ?? null;
    const digestValid = digest !== null && /^[0-9a-f]{64}$/.test(digest);
    if (digest !== null && !digestValid) {
      return refuse(
        "badge-digest-malformed",
        `${row.caseId}: digest '${digest}'`,
        "receipt digests are 64 lowercase hex chars",
      );
    }
    const verdict = row.verdict.trim().toLowerCase();
    let state: BadgeState;
    if (verdict === "pass") {
      if (!digestValid) {
        return refuse(
          "badge-pass-without-receipt",
          `${row.caseId} claims pass with no receipt digest`,
          "a pass state requires its receipt; re-emit from the evidence plane",
        );
      }
      state = "evidenced-pass";
    } else if (verdict === "reported-only" || verdict === "reported") {
      state = "reported";
    } else {
      // Unknown/absent verdicts are NO-DATA — never blank, never green.
      state = "no-data";
    }
    badges.push({
      caseId: row.caseId,
      state,
      receiptDigest: digestValid ? digest : null,
      comparisonClass: row.comparisonClass,
      verdict: row.verdict,
    });
  }
  return { ok: true, value: badges };
}
