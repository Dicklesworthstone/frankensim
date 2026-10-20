//! Contact-law `(K, alpha)` packs (music bead
//! `frankensim-music-v8-root-3ez8g.13.1`).
//!
//! THE SCHEMA'S THESIS: a contact law without its PAIR context is not
//! a datum. A `(K, alpha)` penalty pair is identified for MATERIAL A
//! against MATERIAL B in a GEOMETRY CLASS by a NAMED method, and it
//! is valid over the force/velocity regime the identification
//! exercised — so every one of those is a required schema field and
//! an incomplete row REFUSES. Lookups outside the declared validity
//! REFUSE rather than extrapolate (the falsifier below executes
//! this).
//!
//! The `graze-advisory` field is load-bearing doctrine, recorded
//! from the executed jawari fixture: `1 < alpha < 2` laws carry an
//! unbounded contact Hessian at the boundary (`d²f ~ p^{alpha-2}`),
//! which stalls FD-Jacobian Newton in tight-graze regimes — stiff
//! grazing fixtures should carry `alpha >= 2` cards.
//!
//! File format (`contact.tsv`, header `frankensim.matdb-contact.v1`),
//! one card per pack directory next to the usual `manifest.tsv`:
//!
//! ```text
//! frankensim.matdb-contact.v1
//! pair<TAB>material_a<TAB><a><TAB>material_b<TAB><b><TAB>geometry<TAB><class>
//! law<TAB>penalty-power<TAB>k_n_per_m_alpha<TAB><K><TAB>alpha<TAB><v><TAB>chi_s_per_m<TAB><v>
//! identification<TAB><method-tag><TAB><detail>
//! validity<TAB>force_n<TAB><lo><TAB><hi>
//! validity<TAB>velocity_m_s<TAB><lo><TAB><hi>
//! graze-advisory<TAB><text>
//! citation<TAB><text>
//! ```

use std::fmt;
use std::path::Path;

/// Typed contact-pack failures.
#[derive(Debug, PartialEq, Eq)]
pub enum ContactPackError {
    /// Pack directory or file missing (the named absent-pack refusal —
    /// never a default `K`).
    Missing {
        /// What was looked for.
        what: String,
    },
    /// Schema violation, by field.
    Schema {
        /// Diagnosis.
        what: String,
    },
    /// A lookup outside the card's identified validity regime.
    OutsideValidity {
        /// Which axis and the offending value.
        what: String,
    },
    /// The card's `(K, alpha)` were fitted to a FrankenSim fixture (BR-03).
    /// A fixture is not a material datum: only the named fixture path reads it.
    FixtureCalibrated {
        /// Pack id.
        pack: String,
    },
}

impl fmt::Display for ContactPackError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ContactPackError::Missing { what } => write!(f, "contact pack missing: {what}"),
            ContactPackError::Schema { what } => write!(f, "contact pack schema: {what}"),
            ContactPackError::OutsideValidity { what } => {
                write!(f, "contact lookup outside validity: {what}")
            }
            ContactPackError::FixtureCalibrated { pack } => write!(
                f,
                "contact pack {pack} is fixture-calibrated, not material data; \
                 tests read it through ContactLawCard::load_named_fixture"
            ),
        }
    }
}

/// Identification method tags whose `(K, alpha)` were fitted to a FrankenSim
/// fixture. [`ContactLawCard::load`] and [`ContactLawCard::parse`] refuse such
/// cards; tags compare trimmed and case-insensitively.
pub const FIXTURE_CALIBRATED_IDENTIFICATIONS: [&str; 2] =
    ["identified-from-fixture", "fixture-calibrated"];

fn fixture_calibrated(method: &str) -> bool {
    let method = method.trim();
    FIXTURE_CALIBRATED_IDENTIFICATIONS
        .iter()
        .any(|tag| tag.eq_ignore_ascii_case(method))
}

impl std::error::Error for ContactPackError {}

/// A loaded contact-law card: the pair + geometry context IS the
/// identity, never just the numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactLawCard {
    /// Pack directory name (receipt spelling).
    pub pack_id: String,
    /// Material doing the contacting (e.g. `steel-string`).
    pub material_a: String,
    /// Material being contacted (e.g. `shaped-bone-bridge`).
    pub material_b: String,
    /// Geometry class the identification used.
    pub geometry_class: String,
    /// Penalty stiffness K [N/m^alpha].
    pub k_n_per_m_alpha: f64,
    /// Penalty exponent alpha (>= 1).
    pub alpha: f64,
    /// Hunt–Crossley internal loss chi [s/m] (0 = elastic).
    pub chi_s_per_m: f64,
    /// Identification method tag + detail.
    pub identification: String,
    /// Contact-force validity [N], inclusive.
    pub force_validity_n: (f64, f64),
    /// Approach-velocity validity [m/s], inclusive.
    pub velocity_validity_m_s: (f64, f64),
    /// The graze-regime advisory (required; doctrine, not decoration).
    pub graze_advisory: String,
    /// Source citation / authorship statement.
    pub citation: String,
}

/// A validated lookup: what [`crate::contact_pack`] hands a consumer
/// for one operating regime. The receipt is the typed provenance the
/// fs-dcontact CONTRACT's wiring note reserved.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactReceipt {
    /// Pack id.
    pub pack_id: String,
    /// `<a>-on-<b> [<geometry>]`, display-ready.
    pub pair_label: String,
    /// Penalty stiffness K [N/m^alpha].
    pub k_n_per_m_alpha: f64,
    /// Penalty exponent.
    pub alpha: f64,
    /// Internal loss chi [s/m].
    pub chi_s_per_m: f64,
    /// The identification statement carried through.
    pub identification: String,
}

/// The pack id (directory name) and `contact.tsv` text.
fn read_pack(dir: &Path) -> Result<(String, String), ContactPackError> {
    let pack_id = dir
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| ContactPackError::Missing {
            what: "pack directory name".to_string(),
        })?
        .to_string();
    let path = dir.join("contact.tsv");
    let text = std::fs::read_to_string(&path).map_err(|_| ContactPackError::Missing {
        what: format!("{} (no default K exists; wire a pack)", path.display()),
    })?;
    Ok((pack_id, text))
}

fn parse_f64(s: &str, what: &str) -> Result<f64, ContactPackError> {
    s.parse::<f64>().map_err(|_| ContactPackError::Schema {
        what: format!("{what}: unparseable number {s:?}"),
    })
}

impl ContactLawCard {
    /// Load `<dir>/contact.tsv`.
    ///
    /// # Errors
    /// [`ContactPackError::Missing`] when the pack or file is absent;
    /// [`ContactPackError::Schema`] on any incomplete or non-physical
    /// row; [`ContactPackError::FixtureCalibrated`] for a card fitted to a
    /// fixture.
    pub fn load(dir: &Path) -> Result<ContactLawCard, ContactPackError> {
        let (pack_id, text) = read_pack(dir)?;
        Self::parse(&pack_id, &text)
    }

    /// Parse the `contact.tsv` grammar.
    ///
    /// # Errors
    /// [`ContactPackError::Schema`] with the offending field;
    /// [`ContactPackError::FixtureCalibrated`] for a card fitted to a
    /// fixture.
    pub fn parse(pack_id: &str, text: &str) -> Result<ContactLawCard, ContactPackError> {
        Self::parse_with(pack_id, text, false)
    }

    /// The named fixture path: load a fixture-calibrated card for a test
    /// that reproduces its fixture. Only under `cfg(test)` or the
    /// non-default `fixture-packs` feature, which only dev-dependencies
    /// enable; a card that is not fixture-calibrated refuses here.
    ///
    /// # Errors
    /// As [`ContactLawCard::load`], with [`ContactPackError::Schema`] for a
    /// card that is not fixture-calibrated.
    #[cfg(any(test, feature = "fixture-packs"))]
    pub fn load_named_fixture(dir: &Path) -> Result<ContactLawCard, ContactPackError> {
        let (pack_id, text) = read_pack(dir)?;
        Self::parse_named_fixture(&pack_id, &text)
    }

    /// [`ContactLawCard::load_named_fixture`] on text.
    ///
    /// # Errors
    /// As [`ContactLawCard::load_named_fixture`].
    #[cfg(any(test, feature = "fixture-packs"))]
    pub fn parse_named_fixture(
        pack_id: &str,
        text: &str,
    ) -> Result<ContactLawCard, ContactPackError> {
        Self::parse_with(pack_id, text, true)
    }

    fn parse_with(
        pack_id: &str,
        text: &str,
        named_fixture: bool,
    ) -> Result<ContactLawCard, ContactPackError> {
        let mut lines = text.lines();
        if lines.next() != Some("frankensim.matdb-contact.v1") {
            return Err(ContactPackError::Schema {
                what: "missing frankensim.matdb-contact.v1 header".to_string(),
            });
        }
        let mut state = ContactLawRowState::default();
        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            parse_contact_row(&mut state, &cols)?;
        }
        // PAIR-CONTEXT COMPLETENESS: a row without both materials AND
        // the geometry class refuses — that is the schema's thesis.
        let ContactLawRowState {
            material_a,
            material_b,
            geometry,
            law,
            method,
            identification,
            force_v,
            vel_v,
            graze,
            citation,
        } = state;
        let material_a = material_a.ok_or_else(|| ContactPackError::Schema {
            what: "pair context incomplete: material_a missing".to_string(),
        })?;
        let material_b = material_b.ok_or_else(|| ContactPackError::Schema {
            what: "pair context incomplete: material_b missing".to_string(),
        })?;
        let geometry_class = geometry.ok_or_else(|| ContactPackError::Schema {
            what: "pair context incomplete: geometry class missing".to_string(),
        })?;
        let (k, alpha, chi) = law.ok_or_else(|| ContactPackError::Schema {
            what: "law row missing".to_string(),
        })?;
        if !(k.is_finite() && k > 0.0) {
            return Err(ContactPackError::Schema {
                what: format!("K must be positive, got {k}"),
            });
        }
        if !(alpha.is_finite() && alpha >= 1.0) {
            return Err(ContactPackError::Schema {
                what: format!("alpha must be >= 1, got {alpha}"),
            });
        }
        if !(chi.is_finite() && chi >= 0.0) {
            return Err(ContactPackError::Schema {
                what: format!("chi must be non-negative, got {chi}"),
            });
        }
        // BR-03: a law fitted to a fixture is the fixture's, not a material's.
        let fixture = method.as_deref().is_some_and(fixture_calibrated);
        if fixture && !named_fixture {
            return Err(ContactPackError::FixtureCalibrated {
                pack: pack_id.to_string(),
            });
        }
        if named_fixture && !fixture {
            return Err(ContactPackError::Schema {
                what: format!("{pack_id} is not fixture-calibrated; load it with load"),
            });
        }
        Ok(ContactLawCard {
            pack_id: pack_id.to_string(),
            material_a,
            material_b,
            geometry_class,
            k_n_per_m_alpha: k,
            alpha,
            chi_s_per_m: chi,
            identification: identification.ok_or_else(|| ContactPackError::Schema {
                what: "identification row missing".to_string(),
            })?,
            force_validity_n: force_v.ok_or_else(|| ContactPackError::Schema {
                what: "force validity missing".to_string(),
            })?,
            velocity_validity_m_s: vel_v.ok_or_else(|| ContactPackError::Schema {
                what: "velocity validity missing".to_string(),
            })?,
            graze_advisory: graze.ok_or_else(|| ContactPackError::Schema {
                what: "graze advisory missing (doctrine, not decoration)".to_string(),
            })?,
            citation: citation.ok_or_else(|| ContactPackError::Schema {
                what: "citation missing".to_string(),
            })?,
        })
    }

    /// Display label `<a>-on-<b> [<geometry>]`.
    #[must_use]
    pub fn pair_label(&self) -> String {
        format!(
            "{}-on-{} [{}]",
            self.material_a, self.material_b, self.geometry_class
        )
    }

    /// Validate an operating regime against the card and mint the
    /// typed receipt. THE FALSIFIER LIVES HERE: outside the declared
    /// force/velocity validity the lookup REFUSES, never extrapolates.
    ///
    /// # Errors
    /// [`ContactPackError::OutsideValidity`] by axis.
    pub fn lookup(
        &self,
        peak_force_n: f64,
        peak_velocity_m_s: f64,
    ) -> Result<ContactReceipt, ContactPackError> {
        let (flo, fhi) = self.force_validity_n;
        if !(peak_force_n >= flo && peak_force_n <= fhi) {
            return Err(ContactPackError::OutsideValidity {
                what: format!(
                    "force {peak_force_n} N outside [{flo}, {fhi}] N ({})",
                    self.pack_id
                ),
            });
        }
        let (vlo, vhi) = self.velocity_validity_m_s;
        if !(peak_velocity_m_s >= vlo && peak_velocity_m_s <= vhi) {
            return Err(ContactPackError::OutsideValidity {
                what: format!(
                    "velocity {peak_velocity_m_s} m/s outside [{vlo}, {vhi}] m/s ({})",
                    self.pack_id
                ),
            });
        }
        Ok(ContactReceipt {
            pack_id: self.pack_id.clone(),
            pair_label: self.pair_label(),
            k_n_per_m_alpha: self.k_n_per_m_alpha,
            alpha: self.alpha,
            chi_s_per_m: self.chi_s_per_m,
            identification: self.identification.clone(),
        })
    }
}

/// Accumulated row state for [`ContactLawCard::parse`].
#[derive(Default)]
struct ContactLawRowState {
    material_a: Option<String>,
    material_b: Option<String>,
    geometry: Option<String>,
    law: Option<(f64, f64, f64)>,
    method: Option<String>,
    identification: Option<String>,
    force_v: Option<(f64, f64)>,
    vel_v: Option<(f64, f64)>,
    graze: Option<String>,
    citation: Option<String>,
}
/// Apply one tab-separated `contact.tsv` row to the accumulated state.
///
/// # Errors
/// [`ContactPackError::Schema`] for unknown row kinds, drifted field
/// names, malformed numbers, or unordered validity ranges — byte-for-byte
/// the same refusals the monolithic parser produced.
fn parse_contact_row(
    state: &mut ContactLawRowState,
    cols: &[&str],
) -> Result<(), ContactPackError> {
    match cols[0] {
        "pair" => parse_pair_row(state, cols),
        "law" => parse_law_row(state, cols),
        "identification" => {
            if cols.len() < 3 {
                return Err(ContactPackError::Schema {
                    what: "identification needs a method tag and detail".to_string(),
                });
            }
            state.method = Some(cols[1].to_string());
            state.identification = Some(format!("{}: {}", cols[1], cols[2]));
            Ok(())
        }
        "validity" => parse_validity_row(state, cols),
        "graze-advisory" | "citation" => parse_advisory_row(state, cols),
        other => Err(ContactPackError::Schema {
            what: format!("unknown row kind {other:?}"),
        }),
    }
}

/// `pair material_a <a> material_b <b> geometry <g>`
fn parse_pair_row(state: &mut ContactLawRowState, cols: &[&str]) -> Result<(), ContactPackError> {
    let mut i = 1;
    while i + 1 < cols.len() {
        match cols[i] {
            "material_a" => state.material_a = Some(cols[i + 1].to_string()),
            "material_b" => state.material_b = Some(cols[i + 1].to_string()),
            "geometry" => state.geometry = Some(cols[i + 1].to_string()),
            other => {
                return Err(ContactPackError::Schema {
                    what: format!("unknown pair field {other:?}"),
                });
            }
        }
        i += 2;
    }
    Ok(())
}

/// `law penalty-power k_n_per_m_alpha <k> alpha <a> chi_s_per_m <chi>`
fn parse_law_row(state: &mut ContactLawRowState, cols: &[&str]) -> Result<(), ContactPackError> {
    if cols.len() < 8 || cols[1] != "penalty-power" {
        return Err(ContactPackError::Schema {
            what: "law row must be penalty-power with k/alpha/chi".to_string(),
        });
    }
    if cols[2] != "k_n_per_m_alpha" || cols[4] != "alpha" || cols[6] != "chi_s_per_m" {
        return Err(ContactPackError::Schema {
            what: "law row field names drifted".to_string(),
        });
    }
    state.law = Some((
        parse_f64(cols[3], "K")?,
        parse_f64(cols[5], "alpha")?,
        parse_f64(cols[7], "chi")?,
    ));
    Ok(())
}

/// `validity force_n|velocity_m_s <lo> <hi>`
fn parse_validity_row(
    state: &mut ContactLawRowState,
    cols: &[&str],
) -> Result<(), ContactPackError> {
    if cols.len() < 4 {
        return Err(ContactPackError::Schema {
            what: "validity row needs axis lo hi".to_string(),
        });
    }
    let lo = parse_f64(cols[2], "validity lo")?;
    let hi = parse_f64(cols[3], "validity hi")?;
    if !(lo.is_finite() && hi.is_finite() && hi > lo) {
        return Err(ContactPackError::Schema {
            what: format!("validity range not ordered: {lo} .. {hi}"),
        });
    }
    match cols[1] {
        "force_n" => state.force_v = Some((lo, hi)),
        "velocity_m_s" => state.vel_v = Some((lo, hi)),
        other => {
            return Err(ContactPackError::Schema {
                what: format!("unknown validity axis {other:?}"),
            });
        }
    }
    Ok(())
}

/// Single-cell advisory rows (`graze-advisory`, `citation`).
fn parse_advisory_row(
    state: &mut ContactLawRowState,
    cols: &[&str],
) -> Result<(), ContactPackError> {
    if cols.len() < 2 || cols[1].trim().is_empty() {
        return Err(ContactPackError::Schema {
            what: if cols[0] == "graze-advisory" {
                "empty graze advisory".to_string()
            } else {
                "empty citation".to_string()
            },
        });
    }
    let cell = Some(cols[1].to_string());
    if cols[0] == "graze-advisory" {
        state.graze = cell;
    } else {
        state.citation = cell;
    }
    Ok(())
}

#[cfg(test)]
mod contact_pack_tests {
    use super::*;

    fn repo_root() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("workspace root")
            .to_path_buf()
    }

    const GOOD: &str = "frankensim.matdb-contact.v1\n\
        pair\tmaterial_a\tsteel-string\tmaterial_b\tshaped-bone-bridge\tgeometry\tjawari-graded-parabolic\n\
        law\tpenalty-power\tk_n_per_m_alpha\t1.0e8\talpha\t2.0\tchi_s_per_m\t0.0\n\
        identification\tidentified-from-fixture\tdetail\n\
        validity\tforce_n\t0.0\t50.0\n\
        validity\tvelocity_m_s\t0.0\t5.0\n\
        graze-advisory\talpha >= 2 for tight-graze fixtures\n\
        citation\tauthored\n";

    #[test]
    fn cp_001_committed_packs_load_and_round_trip() {
        let seed = repo_root().join("data/matdb/seed-v1");
        let mut loaded = 0usize;
        for pack in [
            "contact-jawari-bone-bridge",
            "contact-cane-reed-on-mouthpiece-lay",
            "contact-string-on-fret-nickel",
        ] {
            // The jawari card is fixture-calibrated (BR-03): only the named
            // fixture path reads it (cp-005 executes the refusal).
            let card = if pack == "contact-jawari-bone-bridge" {
                ContactLawCard::load_named_fixture(&seed.join(pack))
            } else {
                ContactLawCard::load(&seed.join(pack))
            }
            .expect(pack);
            assert!(!card.material_a.is_empty() && !card.material_b.is_empty());
            assert!(card.k_n_per_m_alpha > 0.0 && card.alpha >= 1.0);
            assert!(!card.graze_advisory.is_empty());
            // Round-trip identity through the receipt at an in-validity
            // operating point.
            let mid_f = f64::midpoint(card.force_validity_n.0, card.force_validity_n.1);
            let mid_v = f64::midpoint(card.velocity_validity_m_s.0, card.velocity_validity_m_s.1);
            let receipt = card.lookup(mid_f, mid_v).expect("in-validity lookup");
            assert_eq!(
                receipt.k_n_per_m_alpha.to_bits(),
                card.k_n_per_m_alpha.to_bits()
            );
            assert_eq!(receipt.alpha.to_bits(), card.alpha.to_bits());
            assert!(receipt.pair_label.contains("-on-"));
            loaded += 1;
            println!(
                "{{\"suite\":\"fs-matdb\",\"case\":\"cp-001-pack\",\"pack\":\"{}\",\
                 \"pair\":\"{}\",\"k\":{:.3e},\"alpha\":{},\"chi\":{}}}",
                card.pack_id,
                card.pair_label(),
                card.k_n_per_m_alpha,
                card.alpha,
                card.chi_s_per_m
            );
        }
        assert_eq!(loaded, 3);
    }

    #[test]
    fn cp_002_pair_context_completeness_refuses() {
        // The schema's thesis, executed: strip each pair field in turn.
        for (strip, expect) in [
            ("material_a\tsteel-string\t", "material_a missing"),
            ("material_b\tshaped-bone-bridge\t", "material_b missing"),
            (
                "\tgeometry\tjawari-graded-parabolic",
                "geometry class missing",
            ),
        ] {
            let text = GOOD.replace(strip, "");
            let err = ContactLawCard::parse("x", &text).expect_err("must refuse");
            match err {
                ContactPackError::Schema { what } => {
                    assert!(what.contains(expect), "{what} vs {expect}");
                }
                other => panic!("wrong refusal class: {other:?}"),
            }
        }
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-002-pair-completeness\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn cp_003_out_of_validity_lookup_refuses() {
        // GOOD is the jawari fixture card, so it reads through the named path.
        let card = ContactLawCard::parse_named_fixture("x", GOOD).expect("card");
        card.lookup(10.0, 1.0).expect("inside is fine");
        // The falsifier: epsilon outside each axis refuses.
        assert!(matches!(
            card.lookup(50.0 + 1e-9, 1.0),
            Err(ContactPackError::OutsideValidity { .. })
        ));
        assert!(matches!(
            card.lookup(10.0, 5.0 + 1e-9),
            Err(ContactPackError::OutsideValidity { .. })
        ));
        assert!(matches!(
            card.lookup(-1.0, 1.0),
            Err(ContactPackError::OutsideValidity { .. })
        ));
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-003-validity-falsifier\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn cp_004_law_and_absence_refusals() {
        assert!(matches!(
            ContactLawCard::load(&repo_root().join("data/matdb/seed-v1/contact-nonexistent")),
            Err(ContactPackError::Missing { .. })
        ));
        for (from, to) in [
            ("alpha\t2.0", "alpha\t0.5"),
            ("k_n_per_m_alpha\t1.0e8", "k_n_per_m_alpha\t-1.0"),
            ("chi_s_per_m\t0.0", "chi_s_per_m\t-0.1"),
            (
                "validity\tforce_n\t0.0\t50.0",
                "validity\tforce_n\t50.0\t0.0",
            ),
        ] {
            let text = GOOD.replace(from, to);
            assert!(
                matches!(
                    ContactLawCard::parse("x", &text),
                    Err(ContactPackError::Schema { .. })
                ),
                "must refuse {to:?}"
            );
        }
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-004-law-refusals\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn cp_005_fixture_calibrated_card_refuses_through_load_and_parse() {
        // BR-03: the committed jawari card can no longer reach a consumer
        // as a material datum, by file or by text.
        let pack = "contact-jawari-bone-bridge";
        let err = ContactLawCard::load(&repo_root().join("data/matdb/seed-v1").join(pack))
            .expect_err("fixture-calibrated card must refuse");
        assert_eq!(
            err,
            ContactPackError::FixtureCalibrated {
                pack: pack.to_string()
            }
        );
        assert!(err.to_string().contains("fixture-calibrated"), "{err}");
        assert_eq!(
            ContactLawCard::parse("x", GOOD),
            Err(ContactPackError::FixtureCalibrated {
                pack: "x".to_string()
            })
        );
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-005-fixture-refusal\",\"pack\":\"{pack}\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn cp_006_named_fixture_path_reads_the_card_bitwise() {
        let seed = repo_root().join("data/matdb/seed-v1");
        let card = ContactLawCard::load_named_fixture(&seed.join("contact-jawari-bone-bridge"))
            .expect("named fixture path");
        assert_eq!(card.k_n_per_m_alpha.to_bits(), 1.0e8_f64.to_bits());
        assert_eq!(card.alpha.to_bits(), 2.0_f64.to_bits());
        assert!(card.identification.starts_with("identified-from-fixture: "));
        // The named path is not a general loader: a card that is not
        // fixture-calibrated refuses there.
        assert!(matches!(
            ContactLawCard::load_named_fixture(&seed.join("contact-string-on-fret-nickel")),
            Err(ContactPackError::Schema { .. })
        ));
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-006-named-fixture\",\"k\":{:e},\"alpha\":{},\"verdict\":\"pass\"}}",
            card.k_n_per_m_alpha, card.alpha
        );
    }

    #[test]
    fn cp_007_refusal_follows_the_identification_not_the_name() {
        // Renaming the pack or recasing the tag does not admit a fixture;
        // a source-backed identification of the same law parses.
        for tag in [
            "identified-from-fixture",
            "Identified-From-Fixture",
            " fixture-calibrated ",
        ] {
            let text = GOOD.replace("identified-from-fixture", tag);
            assert_eq!(
                ContactLawCard::parse("renamed-pack", &text),
                Err(ContactPackError::FixtureCalibrated {
                    pack: "renamed-pack".to_string()
                }),
                "{tag:?}"
            );
        }
        let literature = GOOD.replace("identified-from-fixture", "literature");
        let card = ContactLawCard::parse("x", &literature).expect("source-backed card");
        assert!(card.identification.starts_with("literature: "));
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-007-identification-not-name\",\"verdict\":\"pass\"}}"
        );
    }

    fn rust_sources(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                rust_sources(&path, out);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }

    /// Every `Cargo.toml` under `dir`: the workspace root, crates, apps and
    /// nested example or fuzz manifests. Build output and data are skipped.
    fn manifests(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if path.is_dir() {
                if !(name.starts_with('.')
                    || name == "target"
                    || name == "node_modules"
                    || name == "data")
                {
                    manifests(&path, out);
                }
            } else if name == "Cargo.toml" {
                out.push(path);
            }
        }
    }

    #[test]
    fn cp_008_named_fixture_path_is_not_reachable_from_runtime_code() {
        // "Only tests" made checkable. No manifest anywhere in the tree may
        // enable fixture-packs outside a dev-dependencies table; that includes
        // [workspace.dependencies], through which a member could inherit it.
        // fs-matdb itself may declare it under [features]. No crate's src/
        // may name the fixture path.
        let root = repo_root();
        let own_manifest = root.join("crates/fs-matdb/Cargo.toml");
        let mut found = Vec::new();
        manifests(&root, &mut found);
        found.sort();
        let mut checked = 0usize;
        for path in &found {
            let Ok(manifest) = std::fs::read_to_string(path) else {
                continue;
            };
            let own = *path == own_manifest;
            let mut section = String::new();
            for line in manifest.lines().map(str::trim) {
                if line.starts_with('[') {
                    section = line.to_owned();
                } else if line.contains("fixture-packs") && !line.starts_with('#') {
                    let dev = section == "[dev-dependencies]"
                        || section.starts_with("[dev-dependencies.")
                        || (section.starts_with("[target.")
                            && section.contains(".dev-dependencies"));
                    assert!(
                        dev || (own && section == "[features]"),
                        "{} enables fixture-packs under {section}",
                        path.display()
                    );
                }
            }
            if !own {
                let mut sources = Vec::new();
                rust_sources(&path.parent().unwrap().join("src"), &mut sources);
                for source in &sources {
                    let text = std::fs::read_to_string(source).unwrap_or_default();
                    assert!(
                        !text.contains("load_named_fixture")
                            && !text.contains("parse_named_fixture"),
                        "{} names the fixture path",
                        source.display()
                    );
                }
            }
            checked += 1;
        }
        println!(
            "{{\"suite\":\"fs-matdb\",\"case\":\"cp-008-fixture-path-unreachable\",\"manifests_checked\":{checked},\"verdict\":\"pass\"}}"
        );
        assert!(checked > 100, "the scan must see the workspace ({checked})");
        assert!(
            found.contains(&root.join("Cargo.toml")),
            "the workspace root manifest must be scanned"
        );
    }
}
