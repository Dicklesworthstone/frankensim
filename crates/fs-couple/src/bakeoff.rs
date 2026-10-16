//! Bake-off protocol and receipt schema (music bead
//! `frankensim-music-v8-root-3ez8g.1.4`).
//!
//! Adding an image to an instrument menu REQUIRES a bake-off: same cards,
//! same fixture, same QoIs; report residual, budget, failure mode. The
//! outcomes are exactly three — KEEP BOTH (the default when images are
//! orthogonal: time vs frequency domain, phrase vs hold, 1-DOF vs
//! two-mass), keep one FOR A SUBSET of claims, or REFUSE THE NEWCOMER.
//! Doctrine D21 (menus, not winners) forbids deleting a passing orthogonal
//! image, so [`BakeoffOutcome`] structurally CANNOT express "delete the
//! loser" — there is no such variant to construct.
//!
//! Division of honesty: the HARNESS measures (per-QoI values and residuals
//! against a caller-supplied reference, plus deterministic cost counters);
//! the OUTCOME is a reviewed judgment the caller passes in. A bake-off
//! verdict weighs failure modes and claim scopes, which no residual table
//! can adjudicate alone — the receipt exists so the judgment is inspectable
//! and re-derivable, not so it is automated away.
//!
//! Receipts are deterministic byte artifacts under a canonical
//! line-oriented encoding, content-addressed through `fs-blake3` with a
//! domain-separated hash. Wall-clock time never enters a receipt: the
//! budget fields are logical (state count, steps, iteration counts) so the
//! bytes are bit-stable across machines. Real measured samples/sec belongs
//! to the budget lane (bead 3ez8g.2.2), which the registry's `budget_row`
//! field references separately.
//!
//! The registry (`instrument-claims.json`, `xtask check-instrument-claims`)
//! carries `{"kind":"bakeoff","ref":...}` evidence entries pointing at
//! committed receipt files; the `bakeoff` kind is recorded-only there, so
//! the receipt file itself (golden-pinned by the executing test) is the
//! durable evidence object.

use fs_blake3::{ContentHash, hash_domain};
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Domain separator for receipt content addressing.
pub const BAKEOFF_HASH_DOMAIN: &str = "org.frankensim.fs-couple.bakeoff-receipt.v1";
/// First line of every canonical receipt.
pub const BAKEOFF_SCHEMA_LINE: &str = "frankensim-bakeoff-receipt-v1";

/// Typed refusal from receipt construction or decoding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BakeoffError {
    /// A field violates its structural invariant.
    Invalid {
        /// Failed invariant, named.
        what: &'static str,
    },
    /// Canonical bytes failed to decode at a named line.
    Decode {
        /// 1-based line number of the first offense.
        line: usize,
        /// What was expected there.
        what: &'static str,
    },
}

/// One contender's measured half of the comparison.
#[derive(Clone, Debug, PartialEq)]
pub struct ContenderResult {
    /// Image identity as the claims registry spells it (e.g. `modal-zoh`).
    pub image: String,
    /// Owner crates, D23-style.
    pub owner_crates: Vec<String>,
    /// Measured value per QoI. Keys must exactly match the fixture's
    /// reference map.
    pub measured: BTreeMap<String, f64>,
    /// Logical state count (budget-shaped, deterministic).
    pub states: usize,
    /// Steps executed on the fixture.
    pub steps: usize,
    /// Total implicit-solver iterations (0 for explicit/exact methods).
    pub solver_iterations: usize,
    /// Observed failure modes, free text, empty when none observed.
    pub failure_modes: Vec<String>,
}

/// The reviewed verdict. There is deliberately no variant that removes an
/// image: D21 as a type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BakeoffOutcome {
    /// Both images stay on the menu, each with its claim scope.
    KeepBoth {
        /// Claim scope kept by the first contender.
        scope_a: String,
        /// Claim scope kept by the second contender.
        scope_b: String,
    },
    /// One image keeps only a named subset of claims; the other keeps the
    /// rest. Nobody is deleted.
    KeepForSubset {
        /// Image whose claims narrow.
        narrowed: String,
        /// The subset it keeps.
        subset: String,
    },
    /// The newcomer is refused admission to the menu (it never had rows to
    /// delete). The incumbent is untouched.
    RefuseNewcomer {
        /// The refused image id.
        newcomer: String,
        /// Why, with the counter-argument where one exists.
        reason: String,
    },
}

/// A complete, content-addressable bake-off record.
#[derive(Clone, Debug, PartialEq)]
pub struct BakeoffReceipt {
    /// Filling as the registry spells it (e.g. `string`).
    pub filling: String,
    /// Fixture identity: what was run, spelled so a reader can find the
    /// executing test. Never a wall-clock or commit stamp — the receipt's
    /// git history carries those.
    pub fixture: String,
    /// Digest of the shared cards/parameters both contenders consumed.
    pub shared_cards: ContentHash,
    /// Reference (analytic where possible) value per QoI.
    pub reference: BTreeMap<String, f64>,
    /// Exactly two contenders; a bake-off is pairwise by construction.
    pub contenders: [ContenderResult; 2],
    /// The reviewed verdict.
    pub outcome: BakeoffOutcome,
    /// Why the verdict, in reviewable prose.
    pub rationale: String,
    /// Listening-receipt digests where perceptual QoIs were judged.
    pub listening_receipts: Vec<ContentHash>,
}

fn finite(value: f64, what: &'static str) -> Result<f64, BakeoffError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BakeoffError::Invalid { what })
    }
}

fn clean_text<'a>(text: &'a str, what: &'static str) -> Result<&'a str, BakeoffError> {
    if text.trim().is_empty() {
        return Err(BakeoffError::Invalid { what });
    }
    if text.contains('\n') || text.contains('\t') {
        return Err(BakeoffError::Invalid { what });
    }
    Ok(text)
}

impl BakeoffReceipt {
    /// Relative residual of a contender's measurement against the
    /// reference, per QoI: `|measured - reference| / max(|reference|, 1)`.
    /// The `max(.., 1)` floor keeps near-zero references from manufacturing
    /// infinite residuals; references are expected in natural units where
    /// O(1) is meaningful scale.
    #[must_use]
    pub fn residuals(&self, contender: usize) -> BTreeMap<String, f64> {
        let mut out = BTreeMap::new();
        for (qoi, reference) in &self.reference {
            if let Some(measured) = self.contenders[contender].measured.get(qoi) {
                let denom = reference.abs().max(1.0);
                out.insert(qoi.clone(), (measured - reference).abs() / denom);
            }
        }
        out
    }

    /// Validate structural invariants: non-empty identities, matching QoI
    /// key sets, finite values, clean single-line text fields.
    ///
    /// # Errors
    /// The first violated invariant, named.
    pub fn validate(&self) -> Result<(), BakeoffError> {
        clean_text(&self.filling, "filling must be non-empty single-line")?;
        clean_text(&self.fixture, "fixture must be non-empty single-line")?;
        clean_text(&self.rationale, "rationale must be non-empty single-line")?;
        if self.reference.is_empty() {
            return Err(BakeoffError::Invalid {
                what: "reference QoI map must be non-empty",
            });
        }
        for (qoi, value) in &self.reference {
            clean_text(qoi, "QoI names must be non-empty single-line")?;
            finite(*value, "reference values must be finite")?;
        }
        for contender in &self.contenders {
            clean_text(&contender.image, "image id must be non-empty single-line")?;
            if contender.owner_crates.is_empty() {
                return Err(BakeoffError::Invalid {
                    what: "owner_crates must be non-empty (D23)",
                });
            }
            for owner in &contender.owner_crates {
                clean_text(owner, "owner crate names must be non-empty single-line")?;
            }
            if contender.measured.len() != self.reference.len()
                || contender
                    .measured
                    .keys()
                    .zip(self.reference.keys())
                    .any(|(a, b)| a != b)
            {
                return Err(BakeoffError::Invalid {
                    what: "each contender must measure exactly the reference QoI set",
                });
            }
            for value in contender.measured.values() {
                finite(*value, "measured values must be finite")?;
            }
            for mode in &contender.failure_modes {
                clean_text(mode, "failure modes must be non-empty single-line")?;
            }
        }
        if self.contenders[0].image == self.contenders[1].image {
            return Err(BakeoffError::Invalid {
                what: "contenders must be distinct images",
            });
        }
        match &self.outcome {
            BakeoffOutcome::KeepBoth { scope_a, scope_b } => {
                clean_text(scope_a, "scope_a must be non-empty single-line")?;
                clean_text(scope_b, "scope_b must be non-empty single-line")?;
            }
            BakeoffOutcome::KeepForSubset { narrowed, subset } => {
                let named = clean_text(narrowed, "narrowed must be non-empty single-line")?;
                clean_text(subset, "subset must be non-empty single-line")?;
                if named != self.contenders[0].image && named != self.contenders[1].image {
                    return Err(BakeoffError::Invalid {
                        what: "narrowed must name one of the contenders",
                    });
                }
            }
            BakeoffOutcome::RefuseNewcomer { newcomer, reason } => {
                let named = clean_text(newcomer, "newcomer must be non-empty single-line")?;
                clean_text(reason, "refusal reason must be non-empty single-line")?;
                if named != self.contenders[0].image && named != self.contenders[1].image {
                    return Err(BakeoffError::Invalid {
                        what: "newcomer must name one of the contenders",
                    });
                }
            }
        }
        Ok(())
    }

    /// Canonical byte encoding: line-oriented `key\tvalue`, fixed field
    /// order, floats as `{:e}` (round-trip exact for f64), no timestamps.
    /// Identical receipts produce identical bytes on every host.
    ///
    /// # Errors
    /// Structural invariants via [`Self::validate`].
    pub fn to_canonical_bytes(&self) -> Result<Vec<u8>, BakeoffError> {
        self.validate()?;
        let mut out = String::new();
        let _ = writeln!(out, "{BAKEOFF_SCHEMA_LINE}");
        let _ = writeln!(out, "filling\t{}", self.filling);
        let _ = writeln!(out, "fixture\t{}", self.fixture);
        let _ = writeln!(out, "shared-cards\t{}", hex(&self.shared_cards));
        for (qoi, value) in &self.reference {
            let _ = writeln!(out, "reference\t{qoi}\t{value:e}");
        }
        for contender in &self.contenders {
            let _ = writeln!(out, "contender\t{}", contender.image);
            let _ = writeln!(out, "owners\t{}", contender.owner_crates.join(","));
            for (qoi, value) in &contender.measured {
                let _ = writeln!(out, "measured\t{qoi}\t{value:e}");
            }
            let _ = writeln!(
                out,
                "budget\tstates={}\tsteps={}\tsolver-iterations={}",
                contender.states, contender.steps, contender.solver_iterations
            );
            for mode in &contender.failure_modes {
                let _ = writeln!(out, "failure-mode\t{mode}");
            }
        }
        match &self.outcome {
            BakeoffOutcome::KeepBoth { scope_a, scope_b } => {
                let _ = writeln!(out, "outcome\tkeep-both\t{scope_a}\t{scope_b}");
            }
            BakeoffOutcome::KeepForSubset { narrowed, subset } => {
                let _ = writeln!(out, "outcome\tkeep-for-subset\t{narrowed}\t{subset}");
            }
            BakeoffOutcome::RefuseNewcomer { newcomer, reason } => {
                let _ = writeln!(out, "outcome\trefuse-newcomer\t{newcomer}\t{reason}");
            }
        }
        let _ = writeln!(out, "rationale\t{}", self.rationale);
        for receipt in &self.listening_receipts {
            let _ = writeln!(out, "listening\t{}", hex(receipt));
        }
        Ok(out.into_bytes())
    }

    /// Content identity of the canonical bytes under the receipt domain.
    ///
    /// # Errors
    /// Structural invariants via [`Self::validate`].
    pub fn content_hash(&self) -> Result<ContentHash, BakeoffError> {
        Ok(hash_domain(
            BAKEOFF_HASH_DOMAIN,
            &self.to_canonical_bytes()?,
        ))
    }

    /// Decode canonical bytes back into a receipt (the round-trip half).
    /// Strict: unknown lines, wrong order classes, or malformed fields
    /// refuse with the offending line number.
    ///
    /// # Errors
    /// [`BakeoffError::Decode`] naming the first offending line.
    #[allow(clippy::too_many_lines)] // one coherent canonical decoder
    #[allow(clippy::items_after_statements)] // local decode helpers live next to their use
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, BakeoffError> {
        let text = std::str::from_utf8(bytes).map_err(|_| BakeoffError::Decode {
            line: 0,
            what: "utf-8",
        })?;
        let mut lines = text.lines().enumerate();
        let mut next = || lines.next().map(|(index, line)| (index + 1, line));

        let (line, schema) = next().ok_or(BakeoffError::Decode {
            line: 1,
            what: "schema line",
        })?;
        if schema != BAKEOFF_SCHEMA_LINE {
            return Err(BakeoffError::Decode {
                line,
                what: "schema line",
            });
        }
        fn tagged<'a>(
            entry: Option<(usize, &'a str)>,
            tag: &'static str,
        ) -> Result<(usize, &'a str), BakeoffError> {
            let (line, text) = entry.ok_or(BakeoffError::Decode {
                line: usize::MAX,
                what: tag,
            })?;
            text.strip_prefix(tag)
                .and_then(|rest| rest.strip_prefix('\t'))
                .map(|rest| (line, rest))
                .ok_or(BakeoffError::Decode { line, what: tag })
        }
        fn parse_f64(line: usize, text: &str) -> Result<f64, BakeoffError> {
            text.parse().map_err(|_| BakeoffError::Decode {
                line,
                what: "float",
            })
        }
        fn parse_hash(line: usize, text: &str) -> Result<ContentHash, BakeoffError> {
            let bytes = text.as_bytes();
            if bytes.len() != 64 {
                return Err(BakeoffError::Decode {
                    line,
                    what: "64-hex hash",
                });
            }
            let mut out = [0u8; 32];
            for (index, chunk) in bytes.as_chunks::<2>().0.iter().enumerate() {
                let hi = hex_val(chunk[0]).ok_or(BakeoffError::Decode {
                    line,
                    what: "hex digit",
                })?;
                let lo = hex_val(chunk[1]).ok_or(BakeoffError::Decode {
                    line,
                    what: "hex digit",
                })?;
                out[index] = (hi << 4) | lo;
            }
            Ok(ContentHash(out))
        }

        let (_, filling) = tagged(next(), "filling")?;
        let (_, fixture) = tagged(next(), "fixture")?;
        let (line, cards) = tagged(next(), "shared-cards")?;
        let shared_cards = parse_hash(line, cards)?;

        let mut reference = BTreeMap::new();
        let mut pending = next();
        while let Some((line, text)) = pending {
            let Some(rest) = text.strip_prefix("reference\t") else {
                pending = Some((line, text));
                break;
            };
            let (qoi, value) = rest.split_once('\t').ok_or(BakeoffError::Decode {
                line,
                what: "reference qoi\\tvalue",
            })?;
            reference.insert(qoi.to_string(), parse_f64(line, value)?);
            pending = next();
        }

        let mut contenders = Vec::new();
        while contenders.len() < 2 {
            let (line, text) = pending.ok_or(BakeoffError::Decode {
                line: usize::MAX,
                what: "contender",
            })?;
            let image = text
                .strip_prefix("contender\t")
                .ok_or(BakeoffError::Decode {
                    line,
                    what: "contender",
                })?
                .to_string();
            let (_, owners) = tagged(next(), "owners")?;
            let owner_crates: Vec<String> = owners.split(',').map(str::to_string).collect();
            let mut measured = BTreeMap::new();
            pending = next();
            while let Some((line, text)) = pending {
                let Some(rest) = text.strip_prefix("measured\t") else {
                    pending = Some((line, text));
                    break;
                };
                let (qoi, value) = rest.split_once('\t').ok_or(BakeoffError::Decode {
                    line,
                    what: "measured qoi\\tvalue",
                })?;
                measured.insert(qoi.to_string(), parse_f64(line, value)?);
                pending = next();
            }
            let (line, budget) = tagged(pending, "budget")?;
            let mut states = None;
            let mut steps = None;
            let mut solver_iterations = None;
            for part in budget.split('\t') {
                let (key, value) = part.split_once('=').ok_or(BakeoffError::Decode {
                    line,
                    what: "budget key=value",
                })?;
                let parsed: usize = value.parse().map_err(|_| BakeoffError::Decode {
                    line,
                    what: "budget integer",
                })?;
                match key {
                    "states" => states = Some(parsed),
                    "steps" => steps = Some(parsed),
                    "solver-iterations" => solver_iterations = Some(parsed),
                    _ => {
                        return Err(BakeoffError::Decode {
                            line,
                            what: "budget key",
                        });
                    }
                }
            }
            let mut failure_modes = Vec::new();
            pending = next();
            while let Some((line, text)) = pending {
                let Some(rest) = text.strip_prefix("failure-mode\t") else {
                    pending = Some((line, text));
                    break;
                };
                failure_modes.push(rest.to_string());
                pending = next();
            }
            contenders.push(ContenderResult {
                image,
                owner_crates,
                measured,
                states: states.ok_or(BakeoffError::Decode {
                    line,
                    what: "budget states",
                })?,
                steps: steps.ok_or(BakeoffError::Decode {
                    line,
                    what: "budget steps",
                })?,
                solver_iterations: solver_iterations.ok_or(BakeoffError::Decode {
                    line,
                    what: "budget solver-iterations",
                })?,
                failure_modes,
            });
        }

        let (line, outcome_text) = tagged(pending, "outcome")?;
        let mut parts = outcome_text.splitn(3, '\t');
        let kind = parts.next().unwrap_or_default();
        let first = parts.next().ok_or(BakeoffError::Decode {
            line,
            what: "outcome fields",
        })?;
        let second = parts.next().ok_or(BakeoffError::Decode {
            line,
            what: "outcome fields",
        })?;
        let outcome = match kind {
            "keep-both" => BakeoffOutcome::KeepBoth {
                scope_a: first.to_string(),
                scope_b: second.to_string(),
            },
            "keep-for-subset" => BakeoffOutcome::KeepForSubset {
                narrowed: first.to_string(),
                subset: second.to_string(),
            },
            "refuse-newcomer" => BakeoffOutcome::RefuseNewcomer {
                newcomer: first.to_string(),
                reason: second.to_string(),
            },
            _ => {
                return Err(BakeoffError::Decode {
                    line,
                    what: "outcome kind",
                });
            }
        };

        let (_, rationale) = tagged(next(), "rationale")?;
        let mut listening_receipts = Vec::new();
        while let Some((line, text)) = next() {
            let rest = text
                .strip_prefix("listening\t")
                .ok_or(BakeoffError::Decode {
                    line,
                    what: "listening",
                })?;
            listening_receipts.push(parse_hash(line, rest)?);
        }

        let [a, b]: [ContenderResult; 2] =
            contenders.try_into().map_err(|_| BakeoffError::Decode {
                line: usize::MAX,
                what: "two contenders",
            })?;
        let receipt = Self {
            filling: filling.to_string(),
            fixture: fixture.to_string(),
            shared_cards,
            reference,
            contenders: [a, b],
            outcome,
            rationale: rationale.to_string(),
            listening_receipts,
        };
        receipt.validate().map_err(|_| BakeoffError::Decode {
            line: usize::MAX,
            what: "validated receipt",
        })?;
        Ok(receipt)
    }
}

fn hex(hash: &ContentHash) -> String {
    let mut out = String::with_capacity(64);
    for byte in hash.0 {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

const fn hex_val(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contender(image: &str, offset: f64) -> ContenderResult {
        let mut measured = BTreeMap::new();
        measured.insert("q1".to_string(), 1.0 + offset);
        measured.insert("v1".to_string(), -0.5 + offset);
        ContenderResult {
            image: image.to_string(),
            owner_crates: vec!["fs-couple".to_string()],
            measured,
            states: 2,
            steps: 100,
            solver_iterations: 0,
            failure_modes: vec![],
        }
    }

    fn receipt() -> BakeoffReceipt {
        let mut reference = BTreeMap::new();
        reference.insert("q1".to_string(), 1.0);
        reference.insert("v1".to_string(), -0.5);
        BakeoffReceipt {
            filling: "string".to_string(),
            fixture: "unit fixture".to_string(),
            shared_cards: hash_domain(BAKEOFF_HASH_DOMAIN, b"cards"),
            reference,
            contenders: [
                contender("modal-zoh", 0.0),
                contender("phs-modal-bank", 1e-3),
            ],
            outcome: BakeoffOutcome::KeepBoth {
                scope_a: "performance image".to_string(),
                scope_b: "authority cross-check".to_string(),
            },
            rationale: "orthogonal roles; both pass".to_string(),
            listening_receipts: vec![],
        }
    }

    #[test]
    fn canonical_round_trip_is_exact() {
        let original = receipt();
        let bytes = original.to_canonical_bytes().expect("encode");
        let decoded = BakeoffReceipt::from_canonical_bytes(&bytes).expect("decode");
        assert_eq!(original, decoded);
        let re_encoded = decoded.to_canonical_bytes().expect("re-encode");
        assert_eq!(bytes, re_encoded, "canonical bytes must be stable");
        assert_eq!(
            original.content_hash().expect("hash"),
            decoded.content_hash().expect("hash")
        );
    }

    #[test]
    fn residuals_are_relative_with_unit_floor() {
        let receipt = receipt();
        let residuals = receipt.residuals(1);
        assert!((residuals["q1"] - 1e-3).abs() < 1e-12);
        assert!((residuals["v1"] - 1e-3).abs() < 1e-12);
        let exact = receipt.residuals(0);
        #[allow(clippy::float_cmp)] // EXACT zero-residual pin
        {
            assert_eq!(exact["q1"], 0.0);
        }
    }

    #[test]
    fn outcome_has_no_deletion_arm() {
        // D21 as a type: enumerate every variant; each names surviving
        // scopes or a refused NEWCOMER (who never had rows), never a
        // deletion of an admitted image. This test is the structural
        // assertion the bead demands — adding a `Delete` variant breaks
        // this match and forces the reviewer to read this comment.
        let outcomes = [
            BakeoffOutcome::KeepBoth {
                scope_a: "a".into(),
                scope_b: "b".into(),
            },
            BakeoffOutcome::KeepForSubset {
                narrowed: "modal-zoh".into(),
                subset: "polyphony".into(),
            },
            BakeoffOutcome::RefuseNewcomer {
                newcomer: "modal-zoh".into(),
                reason: "budget strictly worse on every claim".into(),
            },
        ];
        for outcome in outcomes {
            match outcome {
                BakeoffOutcome::KeepBoth { .. }
                | BakeoffOutcome::KeepForSubset { .. }
                | BakeoffOutcome::RefuseNewcomer { .. } => {}
            }
        }
    }

    #[test]
    fn validation_refuses_structural_defects() {
        let mut bad = receipt();
        bad.contenders[1].image = "modal-zoh".to_string();
        assert!(matches!(bad.validate(), Err(BakeoffError::Invalid { .. })));

        let mut bad = receipt();
        bad.contenders[0].measured.remove("q1");
        assert!(matches!(bad.validate(), Err(BakeoffError::Invalid { .. })));

        let mut bad = receipt();
        bad.reference.insert("q1".to_string(), f64::NAN);
        assert!(matches!(bad.validate(), Err(BakeoffError::Invalid { .. })));

        let mut bad = receipt();
        bad.rationale = "two\nlines".to_string();
        assert!(matches!(bad.validate(), Err(BakeoffError::Invalid { .. })));

        let mut bad = receipt();
        bad.outcome = BakeoffOutcome::RefuseNewcomer {
            newcomer: "someone-else".to_string(),
            reason: "not a contender".to_string(),
        };
        assert!(matches!(bad.validate(), Err(BakeoffError::Invalid { .. })));
    }

    #[test]
    fn decode_refuses_malformed_lines() {
        let bytes = receipt().to_canonical_bytes().expect("encode");
        let text = String::from_utf8(bytes).expect("utf8");
        let tampered = text.replace(BAKEOFF_SCHEMA_LINE, "frankensim-bakeoff-receipt-v9");
        assert!(matches!(
            BakeoffReceipt::from_canonical_bytes(tampered.as_bytes()),
            Err(BakeoffError::Decode { line: 1, .. })
        ));
        let truncated = &text[..text.len() / 2];
        assert!(BakeoffReceipt::from_canonical_bytes(truncated.as_bytes()).is_err());
    }
}
