//! Human listening receipts (music bead `frankensim-music-v8-root-3ez8g.1.2`).
//!
//! The crate's pinned law — [`crate::LISTENING_LAW`]: psychoacoustic
//! metrics are never a substitute for human listening — needs the human
//! half to be a durable ARTIFACT, not a vibe in a chat log. A listening
//! receipt records who listened, to exactly which rendered artifact
//! (content digest), judging exactly which question, with what verdict,
//! with psychoacoustic metrics attached as EVIDENCE FIELDS — supporting
//! context, never the verdict.
//!
//! Three verdicts, and the third is load-bearing:
//! [`ListeningVerdict::Unadjudicated`] is a first-class state — a rendered
//! fixture awaiting an ear is a real, recordable situation — but the type
//! system makes it USELESS as pass evidence: only an adjudicated receipt
//! answers [`ListeningReceipt::supports_pass`], so a claims-registry row
//! can never cite an unadjudicated receipt as if a human approved it.
//!
//! Calibration honesty is inherited structurally: the metrics block's
//! absolute-SPL field is an `Option` that stays `None` without a
//! [`crate::Calibration`]-derived value — empty, never fabricated (the
//! same refusal `spl_from_pcm_rms` enforces). Canonical bytes are
//! line-oriented, deterministic, and content-addressed by the CALLER's
//! hash (this crate depends only on fs-fft/fs-math; the digest of the
//! receipt bytes is computed where fs-blake3 is in scope, exactly like
//! WAV artifacts). Wall-clock never enters the bytes: `session` is a
//! caller-declared label (a date string is fine — it is data the human
//! wrote, not a machine clock read).

use std::fmt::Write as _;

/// Schema line of every canonical receipt.
pub const LISTENING_RECEIPT_SCHEMA: &str = "frankensim-listening-receipt-v1";

/// The human verdict. `Unadjudicated` is first-class and useless for
/// gates by construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ListeningVerdict {
    /// A named human judged the question and passed it.
    Pass,
    /// A named human judged the question and failed it.
    Fail,
    /// Rendered and metric-annotated, awaiting a human ear.
    Unadjudicated,
}

impl ListeningVerdict {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Fail => "fail",
            Self::Unadjudicated => "unadjudicated",
        }
    }
}

/// Psychoacoustic metrics attached as evidence. Every field is optional:
/// absent means "not computed", never zero. `spl_db` additionally means
/// "not computable without a calibration" when `None` — the calibration
/// refusal made structural.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AttachedMetrics {
    /// Stationary loudness [sone], if computed.
    pub loudness_sone: Option<f64>,
    /// DIN 45692 sharpness [acum], if computed.
    pub sharpness_acum: Option<f64>,
    /// Log-attack-time [log10 s], if computed.
    pub log_attack_time: Option<f64>,
    /// Absolute SPL [dB]; `None` without a calibration — never fabricated.
    pub spl_db: Option<f64>,
}

/// Typed refusal from receipt construction or decoding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReceiptError {
    /// A field violates its structural invariant.
    Invalid {
        /// The named invariant.
        what: &'static str,
    },
    /// Canonical bytes failed to decode.
    Decode {
        /// 1-based line of the first offense (0 = not line-addressable).
        line: usize,
        /// What was expected.
        what: &'static str,
    },
}

/// One durable human-listening record.
#[derive(Clone, Debug, PartialEq)]
pub struct ListeningReceipt {
    /// Who listened (a name/handle, human-owned).
    pub listener: String,
    /// Caller-declared session label (e.g. an ISO date the human wrote).
    pub session: String,
    /// What was listened to: the artifact's content digest as lowercase
    /// hex (the music renderer's provenance `wav_blake3`).
    pub artifact_hex: String,
    /// Where the artifact's provenance lives (repo-relative path or URI).
    pub artifact_ref: String,
    /// The exact question judged ("does the attack read as reed?").
    pub question: String,
    /// The verdict.
    pub verdict: ListeningVerdict,
    /// Free-text observations (single line; may be empty for pass).
    pub observations: String,
    /// Evidence metrics (context, never the verdict).
    pub metrics: AttachedMetrics,
}

fn clean_line<'a>(text: &'a str, what: &'static str) -> Result<&'a str, ReceiptError> {
    if text.trim().is_empty() {
        return Err(ReceiptError::Invalid { what });
    }
    if text.contains('\n') || text.contains('\t') {
        return Err(ReceiptError::Invalid { what });
    }
    Ok(text)
}

fn opt_field(out: &mut String, key: &str, value: Option<f64>) {
    match value {
        Some(v) => {
            let _ = writeln!(out, "{key}\t{v:e}");
        }
        None => {
            let _ = writeln!(out, "{key}\tabsent");
        }
    }
}

impl ListeningReceipt {
    /// Whether this receipt can serve as PASS evidence for a gate. Only
    /// an adjudicated pass qualifies; unadjudicated receipts exist to be
    /// listened to, not cited.
    #[must_use]
    pub fn supports_pass(&self) -> bool {
        self.verdict == ListeningVerdict::Pass
    }

    /// Validate structural invariants.
    ///
    /// # Errors
    /// The first violated invariant, named. An adjudicated FAIL without
    /// observations refuses: a human who failed a fixture owes the next
    /// agent the reason.
    pub fn validate(&self) -> Result<(), ReceiptError> {
        clean_line(&self.listener, "listener must be a non-empty single line")?;
        if self.verdict != ListeningVerdict::Unadjudicated {
            // An adjudicated verdict names a real human, not a placeholder.
            let listener = self.listener.to_ascii_lowercase();
            if listener == "pending" || listener == "unadjudicated" || listener == "nobody" {
                return Err(ReceiptError::Invalid {
                    what: "adjudicated verdicts need a real listener name",
                });
            }
        }
        clean_line(&self.session, "session must be a non-empty single line")?;
        clean_line(&self.question, "question must be a non-empty single line")?;
        clean_line(
            &self.artifact_ref,
            "artifact_ref must be a non-empty single line",
        )?;
        if self.artifact_hex.len() != 64
            || !self
                .artifact_hex
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
        {
            return Err(ReceiptError::Invalid {
                what: "artifact_hex must be 64 lowercase hex chars",
            });
        }
        if self.verdict == ListeningVerdict::Fail && self.observations.trim().is_empty() {
            return Err(ReceiptError::Invalid {
                what: "a FAIL verdict owes observations (why it failed)",
            });
        }
        if self.observations.contains('\n') || self.observations.contains('\t') {
            return Err(ReceiptError::Invalid {
                what: "observations must be a single line",
            });
        }
        for metric in [
            self.metrics.loudness_sone,
            self.metrics.sharpness_acum,
            self.metrics.log_attack_time,
            self.metrics.spl_db,
        ] {
            if let Some(value) = metric
                && !value.is_finite()
            {
                return Err(ReceiptError::Invalid {
                    what: "attached metrics must be finite when present",
                });
            }
        }
        Ok(())
    }

    /// Canonical deterministic bytes: fixed field order, `{:e}` floats,
    /// `absent` for uncomputed metrics. The caller content-addresses these
    /// bytes where a hash is in scope.
    ///
    /// # Errors
    /// Structural invariants via [`Self::validate`].
    pub fn to_canonical_bytes(&self) -> Result<Vec<u8>, ReceiptError> {
        self.validate()?;
        let mut out = String::new();
        let _ = writeln!(out, "{LISTENING_RECEIPT_SCHEMA}");
        let _ = writeln!(out, "listener\t{}", self.listener);
        let _ = writeln!(out, "session\t{}", self.session);
        let _ = writeln!(out, "artifact\t{}", self.artifact_hex);
        let _ = writeln!(out, "artifact-ref\t{}", self.artifact_ref);
        let _ = writeln!(out, "question\t{}", self.question);
        let _ = writeln!(out, "verdict\t{}", self.verdict.name());
        let _ = writeln!(out, "observations\t{}", self.observations);
        opt_field(&mut out, "loudness-sone", self.metrics.loudness_sone);
        opt_field(&mut out, "sharpness-acum", self.metrics.sharpness_acum);
        opt_field(&mut out, "log-attack-time", self.metrics.log_attack_time);
        opt_field(&mut out, "spl-db", self.metrics.spl_db);
        let _ = writeln!(out, "law\t{}", crate::LISTENING_LAW);
        Ok(out.into_bytes())
    }

    /// Strict round-trip decoder.
    ///
    /// # Errors
    /// [`ReceiptError::Decode`] at the first offending line.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, ReceiptError> {
        let text = std::str::from_utf8(bytes).map_err(|_| ReceiptError::Decode {
            line: 0,
            what: "utf-8",
        })?;
        let lines: Vec<&str> = text.lines().collect();
        if lines.first().copied() != Some(LISTENING_RECEIPT_SCHEMA) {
            return Err(ReceiptError::Decode {
                line: 1,
                what: "schema line",
            });
        }
        let mut cursor = 1usize;
        let mut expect = |tag: &'static str| -> Result<(usize, String), ReceiptError> {
            let line = cursor + 1;
            let content = lines
                .get(cursor)
                .ok_or(ReceiptError::Decode { line, what: tag })?;
            cursor += 1;
            content
                .strip_prefix(tag)
                .and_then(|rest| rest.strip_prefix('\t'))
                .map(|rest| (line, rest.to_string()))
                .ok_or(ReceiptError::Decode { line, what: tag })
        };
        let (_, listener) = expect("listener")?;
        let (_, session) = expect("session")?;
        let (_, artifact_hex) = expect("artifact")?;
        let (_, artifact_ref) = expect("artifact-ref")?;
        let (_, question) = expect("question")?;
        let (line, verdict_text) = expect("verdict")?;
        let verdict = match verdict_text.as_str() {
            "pass" => ListeningVerdict::Pass,
            "fail" => ListeningVerdict::Fail,
            "unadjudicated" => ListeningVerdict::Unadjudicated,
            _ => {
                return Err(ReceiptError::Decode {
                    line,
                    what: "verdict",
                });
            }
        };
        let (_, observations) = expect("observations")?;
        let mut metric = |tag: &'static str| -> Result<Option<f64>, ReceiptError> {
            let (line, value) = expect(tag)?;
            if value == "absent" {
                Ok(None)
            } else {
                value
                    .parse()
                    .map(Some)
                    .map_err(|_| ReceiptError::Decode { line, what: tag })
            }
        };
        let metrics = AttachedMetrics {
            loudness_sone: metric("loudness-sone")?,
            sharpness_acum: metric("sharpness-acum")?,
            log_attack_time: metric("log-attack-time")?,
            spl_db: metric("spl-db")?,
        };
        let (line, law) = expect("law")?;
        if law != crate::LISTENING_LAW {
            return Err(ReceiptError::Decode { line, what: "law" });
        }
        let receipt = Self {
            listener,
            session,
            artifact_hex,
            artifact_ref,
            question,
            verdict,
            observations,
            metrics,
        };
        receipt.validate().map_err(|_| ReceiptError::Decode {
            line: 0,
            what: "validated receipt",
        })?;
        Ok(receipt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receipt(verdict: ListeningVerdict) -> ListeningReceipt {
        ListeningReceipt {
            listener: "jemanuel".to_string(),
            session: "2026-08-15".to_string(),
            artifact_hex: "05ca73192645e1664d5175435bc611a4469d9478f7fc8e2b0c57aae4b7f9491b"
                .to_string(),
            artifact_ref: "data/listening/reed-fixture.provenance.json".to_string(),
            question: "does the sustained tone read as a reed instrument?".to_string(),
            verdict,
            observations: "square-ish clarinet-like timbre; abrupt attack".to_string(),
            metrics: AttachedMetrics {
                loudness_sone: Some(14.2),
                sharpness_acum: Some(1.4),
                log_attack_time: Some(-2.1),
                spl_db: None, // no calibration: stays empty, never fabricated
            },
        }
    }

    #[test]
    fn round_trip_is_exact_for_every_verdict() {
        for verdict in [
            ListeningVerdict::Pass,
            ListeningVerdict::Fail,
            ListeningVerdict::Unadjudicated,
        ] {
            let original = receipt(verdict);
            let bytes = original.to_canonical_bytes().expect("encode");
            let decoded = ListeningReceipt::from_canonical_bytes(&bytes).expect("decode");
            assert_eq!(original, decoded);
            assert_eq!(bytes, decoded.to_canonical_bytes().expect("re-encode"));
        }
    }

    #[test]
    fn only_adjudicated_pass_supports_gates() {
        assert!(receipt(ListeningVerdict::Pass).supports_pass());
        assert!(!receipt(ListeningVerdict::Fail).supports_pass());
        assert!(
            !receipt(ListeningVerdict::Unadjudicated).supports_pass(),
            "an unadjudicated receipt must be structurally useless as pass evidence"
        );
    }

    #[test]
    fn uncalibrated_spl_stays_absent_in_the_bytes() {
        let bytes = receipt(ListeningVerdict::Pass)
            .to_canonical_bytes()
            .expect("encode");
        let text = String::from_utf8(bytes).expect("utf8");
        assert!(
            text.contains("spl-db\tabsent"),
            "no calibration => absent, never a fabricated number"
        );
        assert!(
            text.contains(crate::LISTENING_LAW),
            "the law travels in the bytes"
        );
    }

    #[test]
    fn refusals_fire_by_name() {
        let mut bad = receipt(ListeningVerdict::Fail);
        bad.observations = String::new();
        assert!(matches!(bad.validate(), Err(ReceiptError::Invalid { .. })));

        let mut bad = receipt(ListeningVerdict::Pass);
        bad.listener = "pending".to_string();
        assert!(
            matches!(bad.validate(), Err(ReceiptError::Invalid { .. })),
            "an adjudicated verdict cannot hide behind a placeholder listener"
        );

        let mut bad = receipt(ListeningVerdict::Unadjudicated);
        bad.artifact_hex = "shorty".to_string();
        assert!(matches!(bad.validate(), Err(ReceiptError::Invalid { .. })));

        let mut bad = receipt(ListeningVerdict::Pass);
        bad.metrics.loudness_sone = Some(f64::NAN);
        assert!(matches!(bad.validate(), Err(ReceiptError::Invalid { .. })));

        let mut bad = receipt(ListeningVerdict::Pass);
        bad.question = "two\nlines".to_string();
        assert!(matches!(bad.validate(), Err(ReceiptError::Invalid { .. })));
    }

    #[test]
    fn decode_refuses_tampered_law_and_schema() {
        let bytes = receipt(ListeningVerdict::Pass)
            .to_canonical_bytes()
            .expect("encode");
        let text = String::from_utf8(bytes).expect("utf8");
        let wrong_schema =
            text.replace(LISTENING_RECEIPT_SCHEMA, "frankensim-listening-receipt-v9");
        assert!(ListeningReceipt::from_canonical_bytes(wrong_schema.as_bytes()).is_err());
        let wrong_law = text.replace("never a substitute", "sometimes a substitute");
        assert!(
            ListeningReceipt::from_canonical_bytes(wrong_law.as_bytes()).is_err(),
            "a receipt that rewrites the listening law is refused"
        );
    }

    #[test]
    fn placeholder_listener_is_legal_while_unadjudicated() {
        let mut pending = receipt(ListeningVerdict::Unadjudicated);
        pending.listener = "pending".to_string();
        assert!(
            pending.validate().is_ok(),
            "awaiting-an-ear is a real state"
        );
        assert!(!pending.supports_pass());
    }
}
