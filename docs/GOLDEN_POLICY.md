# Golden-evidence policy: couplings and justified bumps (bead y4pt)

A golden hash is a SENTINEL over bit semantics, not a checksum to be
refreshed when it gets in the way. Two incidents on 2026-07-09
motivate the mechanics here: the fs-la `rand_nla` golden was re-pinned
with a plausible-sounding but WRONG justification (an fs-rand commit
that only changes behavior at the u64::MAX counter boundary), while
the actual mover was optimization-level-dependent `f64::powi` — the
re-pin froze release-mode bits while debug still produced the old
value, so the sentinel's verdict depended on the build profile.

## The coupling registry

`golden-couplings.json` (workspace root) declares, for every golden:
the file and const that carry it, the upstream SEMANTIC SURFACES it
depends on, and the surface version it was frozen against. Surfaces
declare a `pub const <NAME>_VERSION: u32` in source and bump it on ANY
change that can move downstream bits.

`cargo run -p xtask -- check-goldens` (part of `check-all`) fails
when a surface's source const drifts from the registry, or when a
golden's pin lags a surface row — each failure names every dependent
golden that must be deliberately re-frozen. An upstream semantic
change therefore POINTS AT its downstream goldens instead of
surprising them.

Identity-owning surfaces additionally pin their exact domain constant and
literal in the same surface row. `check-goldens` rejects a domain/value drift,
and `check-identities` requires every generated identity declaration to name a
matching surface/version. Changing an identity schema or domain therefore
requires one deliberate version/coupling update; regenerating
`identity-schemas.json` alone cannot bless a silent re-key.

## The justified-bump protocol

A golden re-pin is valid only when ALL of the following hold:

1. **Committed tree.** The new value is reproduced at a clean, committed
   tree — never a dirty working tree (sweeper commits make "what code
   produced this?" unanswerable otherwise).
2. **Both modes.** The value is reproduced in BOTH debug and release
   (and on both reference ISAs where the golden claims cross-ISA
   validity). A value that differs by profile means the OBSERVED
   QUANTITY is build-mode-dependent — that is a bug to fix (see the
   powi incident, bead 4xnt), not a value to freeze.
3. **Plausible root cause.** The named cause must be able to move the
   observed bits. "A dependency changed" is not a cause; the specific
   semantic change is. If the cause cannot plausibly reach the observed
   value, STOP — the golden is telling you about a different bug.
4. **Registry updated.** The golden's `depends_on` pins and
   `justification` in `golden-couplings.json` are updated in the same
   commit, including the evidence trail (where the two-mode logs live).

## Reviewer checklist for any commit touching a golden const

- [ ] Does the commit/bead record two-mode, committed-tree reproduction?
- [ ] Does the stated root cause plausibly move these bits?
- [ ] Is `golden-couplings.json` updated in the same commit?
- [ ] If an upstream surface changed semantics, was its version const
      bumped (so `check-goldens` catches OTHER dependents)?

A re-pin that skips any row above is the incident pattern this policy
exists to prevent.
