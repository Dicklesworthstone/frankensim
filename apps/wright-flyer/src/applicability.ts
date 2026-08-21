// Evidence/applicability plumbing + empty-receipt UX (bead
// wf-root-guzez.9.3, E8.3a). ApplicabilityDomain intersection
// machinery with LIMITING-SUBSYSTEM attribution per bound, badge
// composition (color + icon + sentence + subsystem breakdown +
// receipt link), and honest empty states before receipts exist.
//
// Laws:
//   - the composed badge is GREEN only for an evidenced pass INSIDE
//     the intersected applicability domain; outside is amber and the
//     sentence NAMES the limiting subsystem and bound;
//   - NO-DATA composes gray with the honest empty-state sentence —
//     "no receipt yet" claims nothing, and is never blank;
//   - a receipt link exists only when a receipt digest exists.

import type { EvidenceBadge } from "./evidenceBadges.ts";

export interface ApplicRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type ApplicResult<T> = { ok: true; value: T } | { ok: false; refusal: ApplicRefusal };

function refuse<T>(code: string, message: string, repair: string): ApplicResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

export interface DomainAxis {
  readonly name: string;
  readonly lo: number;
  readonly hi: number;
}

export interface ApplicabilityDomain {
  readonly subsystem: string;
  readonly axes: readonly DomainAxis[];
}

/** Axis cap per domain. */
export const MAX_AXES = 8;
/** Subsystem cap per intersection. */
export const MAX_SUBSYSTEMS = 16;

export interface IntersectedAxis {
  readonly name: string;
  readonly lo: number;
  readonly hi: number;
  /** the subsystem whose lower bound binds. */
  readonly loLimitedBy: string;
  /** the subsystem whose upper bound binds. */
  readonly hiLimitedBy: string;
}

export type Intersection =
  | { readonly kind: "domain"; readonly axes: readonly IntersectedAxis[] }
  | {
      readonly kind: "empty";
      readonly axis: string;
      /** the two subsystems whose bounds cannot both hold. */
      readonly conflict: readonly [string, string];
    };

/**
 * Intersect the subsystems' applicability boxes axis by axis,
 * recording WHICH subsystem binds each bound. An empty intersection
 * is an honest result naming the conflicting pair — never an error
 * to hide.
 */
export function intersectDomains(
  domains: readonly ApplicabilityDomain[],
): ApplicResult<Intersection> {
  if (domains.length === 0 || domains.length > MAX_SUBSYSTEMS) {
    return refuse(
      "applicability-subsystems-invalid",
      `${domains.length} subsystems outside [1, ${MAX_SUBSYSTEMS}]`,
      "intersect a bounded declared set",
    );
  }
  const first = domains[0];
  if (first === undefined) {
    return refuse("applicability-subsystems-invalid", "empty", "unreachable");
  }
  for (const d of domains) {
    if (d.axes.length === 0 || d.axes.length > MAX_AXES) {
      return refuse(
        "applicability-axes-invalid",
        `${d.subsystem}: ${d.axes.length} axes outside [1, ${MAX_AXES}]`,
        "declared axes are few and named",
      );
    }
    if (d.axes.some((a) => !(Number.isFinite(a.lo) && Number.isFinite(a.hi) && a.lo <= a.hi))) {
      return refuse(
        "applicability-axes-invalid",
        `${d.subsystem}: malformed axis`,
        "finite lo <= hi per axis",
      );
    }
  }
  // Every domain must declare the same axis names (same order).
  const names = first.axes.map((a) => a.name);
  for (const d of domains) {
    const dn = d.axes.map((a) => a.name);
    if (dn.length !== names.length || dn.some((n, i) => n !== names[i])) {
      return refuse(
        "applicability-axes-mismatched",
        `${d.subsystem} declares [${dn.join(",")}] vs [${names.join(",")}]`,
        "subsystems declare the shared axis set",
      );
    }
  }
  const axes: IntersectedAxis[] = [];
  for (let i = 0; i < names.length; i += 1) {
    let lo = Number.NEGATIVE_INFINITY;
    let hi = Number.POSITIVE_INFINITY;
    let loBy = "";
    let hiBy = "";
    for (const d of domains) {
      const a = d.axes[i];
      if (a === undefined) continue;
      if (a.lo > lo) {
        lo = a.lo;
        loBy = d.subsystem;
      }
      if (a.hi < hi) {
        hi = a.hi;
        hiBy = d.subsystem;
      }
    }
    if (lo > hi) {
      return {
        ok: true,
        value: {
          kind: "empty",
          axis: names[i] ?? "",
          conflict: [loBy, hiBy],
        },
      };
    }
    axes.push({ name: names[i] ?? "", lo, hi, loLimitedBy: loBy, hiLimitedBy: hiBy });
  }
  return { ok: true, value: { kind: "domain", axes } };
}

export type Standing =
  | { readonly inside: true }
  | {
      readonly inside: false;
      readonly axis: string;
      readonly bound: "lo" | "hi";
      readonly limitedBy: string;
    };

/**
 * Where does an operating point stand? AT a bound is INSIDE (the
 * declared intervals are closed); the first violated axis names its
 * limiting subsystem.
 */
export function standingAt(
  inter: Intersection,
  point: Readonly<Record<string, number>>,
): ApplicResult<Standing> {
  if (inter.kind === "empty") {
    return refuse(
      "applicability-empty-domain",
      `no operating point exists (${inter.axis}: ${inter.conflict[0]} vs ${inter.conflict[1]})`,
      "the conflict is the finding; widen one subsystem's evidence",
    );
  }
  for (const a of inter.axes) {
    const x = point[a.name];
    if (x === undefined || !Number.isFinite(x)) {
      return refuse(
        "applicability-point-invalid",
        `missing/non-finite coordinate '${a.name}'`,
        "supply every declared axis",
      );
    }
    if (x < a.lo) {
      return { ok: true, value: { inside: false, axis: a.name, bound: "lo", limitedBy: a.loLimitedBy } };
    }
    if (x > a.hi) {
      return { ok: true, value: { inside: false, axis: a.name, bound: "hi", limitedBy: a.hiLimitedBy } };
    }
  }
  return { ok: true, value: { inside: true } };
}

export interface ComposedBadge {
  readonly caseId: string;
  readonly color: "green" | "amber" | "gray";
  readonly icon: "check" | "boundary" | "empty";
  readonly sentence: string;
  /** per-subsystem breakdown (verbatim names). */
  readonly subsystems: readonly string[];
  /** present ONLY with a receipt digest. */
  readonly receiptLink: string | null;
}

/**
 * Compose the display badge: evidence state × applicability standing.
 * Green is EARNED (evidenced pass AND inside); amber names the
 * limiting subsystem; gray is the honest empty state.
 */
export function composeBadge(
  badge: EvidenceBadge,
  standing: Standing,
  subsystems: readonly string[],
): ComposedBadge {
  const link = badge.receiptDigest === null ? null : `receipt:${badge.receiptDigest}`;
  if (badge.state === "no-data") {
    return {
      caseId: badge.caseId,
      color: "gray",
      icon: "empty",
      sentence: `no receipt yet for ${badge.caseId} — nothing is claimed`,
      subsystems,
      receiptLink: link,
    };
  }
  if (!standing.inside) {
    return {
      caseId: badge.caseId,
      color: "amber",
      icon: "boundary",
      sentence: `outside declared applicability: ${standing.axis} ${standing.bound} bound, limited by ${standing.limitedBy}`,
      subsystems,
      receiptLink: link,
    };
  }
  if (badge.state === "evidenced-pass") {
    return {
      caseId: badge.caseId,
      color: "green",
      icon: "check",
      sentence: `${badge.caseId} evidenced (${badge.comparisonClass}) within declared applicability`,
      subsystems,
      receiptLink: link,
    };
  }
  return {
    caseId: badge.caseId,
    color: "gray",
    icon: "empty",
    sentence: `${badge.caseId} ${badge.verdict} — reported, not a pass claim`,
    subsystems,
    receiptLink: link,
  };
}
