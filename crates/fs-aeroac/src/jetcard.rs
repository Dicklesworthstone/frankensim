//! Jet cards — the minted boundary object between the fs-aeroac lab
//! lane and the audio-rate performance lane (bead
//! frankensim-music-v8-root-3ez8g.10.2).
//!
//! The menu law downstream is "no card -> refuse": a jet island may
//! only speak from a card, and a card is MINTED from lab receipts,
//! never typed from vibes. This module owns the card schema, the
//! minting rules, and the authority-class law:
//!
//! - `X-Struct` ONLY when the card carries a residual issued against
//!   lab data; `X-Est` otherwise. Minting enforces this; a caller
//!   cannot construct a structural card without its residual.
//! - The fs-aeroac scope law is inherited verbatim: ABSOLUTE SPL IS
//!   NEVER CLAIMED — cards carry shape/scaling authority only, and
//!   every card embeds [`crate::SCOPE_STATEMENT`] (the
//!   marketing-mutation guard asserts its presence).
//!
//! Three claim kinds exist:
//! - the TONAL INTERIM card (explicitly permitted by the bead):
//!   minted from the validated 2-D stage-I edge-tone point
//!   (St 0.03662 vs Brown's 0.03554, +3.0% inside the record's ±6%
//!   bin quantization; hysteresis and multi-stability recorded) with
//!   a TYPED NARROW CLAIM — edge-tone class, NOT flute broadband;
//! - the BROADBAND card: mintable only from a demonstrated 3-D
//!   broadband regime (classified sweep rungs);
//! - the REFUSAL-BOUNDARY card: the honest artifact when the 3-D
//!   sweep produces no broadband regime — "no broadband below Re X
//!   on this rig family" is a useful, typed truth.

use crate::jetlab::JetLabiumConfig;
use crate::noisetable::N_BANDS;
use crate::slot_jet_3d::SlotJet3dRung;
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_math::det;

/// Schema tag embedded in every serialized card; any other value
/// refuses on parse (no silent reinterpretation of foreign bytes).
pub const JET_CARD_SCHEMA: &str = "fs-aeroac.jet-card/v1";

/// Authority class of a jet card.
///
/// The bead law: `X-Struct` only when a residual is issued against
/// lab data; `X-Est` otherwise. The vocabulary matches the
/// instrument-claims registry exactness classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CardAuthority {
    /// Structural authority: the card's numbers carry a measured
    /// residual against lab/published data.
    XStruct,
    /// Estimate authority: no residual is attached.
    XEst,
}

impl CardAuthority {
    /// Registry-vocabulary label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            CardAuthority::XStruct => "X-Struct",
            CardAuthority::XEst => "X-Est",
        }
    }
}

/// A residual issued against lab or published reference data — the
/// admission ticket for `X-Struct` authority.
#[derive(Debug, Clone, PartialEq)]
pub struct CardResidual {
    /// Which quantity the residual is measured on.
    pub quantity: String,
    /// The lab-measured value.
    pub measured: f64,
    /// The reference value the residual is issued against.
    pub reference: f64,
    /// Where the reference comes from (citation text).
    pub reference_source: String,
    /// Relative half-width of the measurement's own quantization or
    /// stated uncertainty band (e.g. FFT bin quantization).
    pub bin_halfwidth_rel: f64,
}

impl CardResidual {
    /// Signed relative deviation `(measured - reference)/reference`.
    #[must_use]
    pub fn relative_deviation(&self) -> f64 {
        (self.measured - self.reference) / self.reference
    }
}

/// The validity region a card claims. Queries outside it refuse —
/// a card is a narrow truth, not an extrapolation license.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidityRegion {
    /// Lowest jet Reynolds number the claim was demonstrated at.
    pub reynolds_lo: f64,
    /// Highest jet Reynolds number the claim was demonstrated at.
    pub reynolds_hi: f64,
    /// Edge distance in slot heights (`h/delta`) of the rig family.
    pub h_over_delta: f64,
    /// Whether the claim requires the jet-root receptivity edge
    /// (nozzle wall): the executed structural finding is that
    /// without it the rig locks to the free jet's own mode instead
    /// of the Brown stage ladder.
    pub requires_receptivity_edge: bool,
    /// Whether amplitude-bearing consumption additionally requires
    /// the deterministic symmetry-breaking seed (the executed
    /// vacuous-oscillation trap: mirror-symmetric rigs oscillate in
    /// amplified roundoff).
    pub seeded_amplitude_claims_only: bool,
}

/// Mean-jet profile parameters at the flue exit (the quantities the
/// island's jet model consumes).
#[derive(Debug, Clone, PartialEq)]
pub struct MeanJetProfile {
    /// Centerline (peak) velocity [lu/step].
    pub u_centerline: f64,
    /// Slot half-height parameter of the smoothed top-hat [lu].
    pub slot_half: f64,
    /// Momentum thickness of the exit profile [lu], computed from
    /// the SAME discrete profile the rig imposes (see
    /// [`momentum_thickness_smoothed_tophat`]).
    pub momentum_thickness: f64,
}

/// Receptivity/gain characterization: the MEASURED edge-tone
/// feedback quantities the island's gain/delay model consumes.
/// Everything here is a recorded lab observation, not a fit.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeToneFeedback {
    /// Brown stage index of the locked state (stage I = 1).
    pub stage: u8,
    /// Locked Strouhal number `f * delta / u` of the recorded state.
    pub locked_strouhal: f64,
    /// The gated ladder band as multiples of Brown's stage
    /// prediction (lo, hi) — the fine-lattice regression gates
    /// membership of this band, not digit equality.
    pub ladder_band: (f64, f64),
    /// Recorded co-existing locked states (multi-stability): the
    /// attractor reached varies with seed and resolution.
    pub multi_stable_strouhal: Vec<f64>,
    /// Whether bitwise-reproducible hysteresis was recorded on the
    /// adiabatic ramp protocol.
    pub hysteresis_recorded: bool,
    /// Saturated transverse-force RMS of the locked limit cycle
    /// [lattice units, shape/scaling authority only], when the lab
    /// recorded one under the seeded-amplitude law. `None` when no
    /// recorded amplitude is citable — the island's saturation
    /// model then has no card authority to lean on.
    pub saturated_force_rms: Option<f64>,
}

/// The typed narrow claim a card makes.
#[derive(Debug, Clone, PartialEq)]
pub enum JetCardClaim {
    /// Edge-tone class tonal claim (the interim card): a locked
    /// stage-ladder oscillation. NOT flute broadband; consuming this
    /// card as a noise source is a type error by construction.
    EdgeToneTonal {
        /// The measured feedback characterization.
        feedback: EdgeToneFeedback,
    },
    /// A demonstrated 3-D broadband regime (spectral flatness above
    /// the tonal ceiling on an amplitude-qualified, in-regime rung).
    Broadband {
        /// Measured spectral flatness of the demonstrating rung.
        flatness: f64,
        /// Reynolds number of the demonstrating rung.
        reynolds: f64,
    },
    /// The honest refusal artifact: the executed sweep found NO
    /// broadband regime below the probed boundary.
    BroadbandRefusalBoundary {
        /// Highest Reynolds number probed with a tonal outcome.
        max_reynolds_probed: f64,
        /// Number of executed rungs backing the boundary.
        rungs_probed: usize,
    },
}

impl JetCardClaim {
    /// Serialization tag.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        match self {
            JetCardClaim::EdgeToneTonal { .. } => "edge-tone-tonal",
            JetCardClaim::Broadband { .. } => "broadband",
            JetCardClaim::BroadbandRefusalBoundary { .. } => "broadband-refusal-boundary",
        }
    }
}

/// Provenance: which rig, which receipts, which classification.
#[derive(Debug, Clone, PartialEq)]
pub struct CardProvenance {
    /// FNV-1a fingerprint of the exact rig configuration the claim
    /// was minted from ([`jet_labium_fingerprint`] for the 2-D rig;
    /// [`crate::slot_jet_3d::config_fingerprint`] for the 3-D rig).
    pub rig_fingerprint: u64,
    /// Receipt references (test names, receipt file paths, bead
    /// threads) a reviewer can chase.
    pub receipts: Vec<String>,
    /// Classification outcome the lab recorded for this claim.
    pub classification: String,
}

/// A minted jet card.
///
/// Constructed only through the minting functions; the invariants
/// (`X-Struct` implies residual, broadband implies band content,
/// scope statement present) are enforced there and re-checked by
/// [`JetCard::validate`].
#[derive(Debug, Clone, PartialEq)]
pub struct JetCard {
    /// Schema tag ([`JET_CARD_SCHEMA`]).
    pub schema: String,
    /// The typed narrow claim.
    pub claim: JetCardClaim,
    /// Authority class under the bead law.
    pub authority: CardAuthority,
    /// Mean-jet exit profile parameters.
    pub profile: MeanJetProfile,
    /// Relative band power-density levels [dB] over the 16
    /// log-spaced Strouhal bands of [`crate::noisetable`], when the
    /// claim carries band noise content. `None` for a tonal claim
    /// (a locked line's location lives in the feedback record, and
    /// minting a noise spectrum for it would launder a limit cycle).
    pub band_db: Option<[f64; N_BANDS]>,
    /// Residual against lab data (present iff authority is
    /// `X-Struct`).
    pub residual: Option<CardResidual>,
    /// Where the claim is valid; queries outside refuse.
    pub validity: ValidityRegion,
    /// Rig fingerprint, receipts, classification.
    pub provenance: CardProvenance,
    /// The scope law, embedded verbatim ([`crate::SCOPE_STATEMENT`]).
    pub scope: String,
}

/// Momentum thickness of the rig's smoothed top-hat exit profile,
/// computed by the SAME discrete row sum the rig imposes
/// (`u(y)/U = (1 + tanh((b - |y - yc|)/w))/2` on integer rows):
/// `theta = sum_y (u/U)(1 - u/U)`. For `b >> w` this approaches the
/// analytic `w` (two tanh edges at `w/2` each).
///
/// # Errors
/// [`AeroacError::InvalidParameter`] on non-positive `slot_half`,
/// `slot_smoothing`, or a domain too small to contain the jet.
pub fn momentum_thickness_smoothed_tophat(
    slot_half: f64,
    slot_smoothing: f64,
    ny: usize,
) -> Result<f64, AeroacError> {
    if !(slot_half > 0.0 && slot_half.is_finite()) {
        return Err(AeroacError::InvalidParameter { what: "slot_half" });
    }
    if !(slot_smoothing > 0.0 && slot_smoothing.is_finite()) {
        return Err(AeroacError::InvalidParameter {
            what: "slot_smoothing",
        });
    }
    #[allow(clippy::cast_precision_loss)]
    let half_domain = ny as f64 / 2.0;
    if half_domain < slot_half + 6.0 * slot_smoothing {
        return Err(AeroacError::InvalidParameter {
            what: "ny too small: the jet profile must decay inside the domain",
        });
    }
    #[allow(clippy::cast_precision_loss)]
    let yc = ny as f64 / 2.0 - 0.5;
    let mut theta = 0.0;
    for y in 0..ny {
        #[allow(clippy::cast_precision_loss)]
        let dy = (y as f64 - yc).abs();
        let frac = 0.5 * (1.0 + det::tanh((slot_half - dy) / slot_smoothing));
        theta += frac * (1.0 - frac);
    }
    Ok(theta)
}

/// Canonical FNV-1a fingerprint of a 2-D jet-labium configuration
/// (the same discipline as
/// [`crate::slot_jet_3d::config_fingerprint`] for the 3-D rig),
/// binding a card to exactly one rig setup.
#[must_use]
pub fn jet_labium_fingerprint(cfg: &JetLabiumConfig) -> u64 {
    let tag = format!(
        "jetcard-v1|{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:?}",
        cfg.nx,
        cfg.ny,
        cfg.slot_half,
        cfg.slot_smoothing,
        cfg.u_jet,
        cfg.tau,
        cfg.edge_distance,
        cfg.plate_length,
        cfg.fringe_width,
        cfg.fringe_sigma,
        cfg.steps_settle,
        cfg.steps_record,
        cfg.seed_amplitude,
        cfg.nozzle_thickness,
        cfg.collision,
    );
    fnv1a(tag.as_bytes())
}

/// FNV-1a over raw bytes (deterministic correlation identity, not
/// cryptographic — the workspace-wide `ProvenanceHash` caveat
/// applies verbatim).
#[must_use]
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

impl JetCard {
    /// Re-check the minting invariants on an existing card (used
    /// after deserialization and by consumers before trusting a
    /// card handed across a boundary).
    ///
    /// # Errors
    /// [`AeroacError::InvalidParameter`] naming the violated law.
    pub fn validate(&self) -> Result<(), AeroacError> {
        if self.schema != JET_CARD_SCHEMA {
            return Err(AeroacError::InvalidParameter {
                what: "unknown jet-card schema (refused by name)",
            });
        }
        if self.scope != SCOPE_STATEMENT {
            return Err(AeroacError::InvalidParameter {
                what: "scope statement missing or mutated (marketing-mutation guard)",
            });
        }
        match self.authority {
            CardAuthority::XStruct => {
                if self.residual.is_none() {
                    return Err(AeroacError::InvalidParameter {
                        what: "X-Struct authority requires a residual against lab data",
                    });
                }
            }
            CardAuthority::XEst => {}
        }
        if matches!(self.claim, JetCardClaim::Broadband { .. }) && self.band_db.is_none() {
            return Err(AeroacError::InvalidParameter {
                what: "a broadband claim requires band noise content",
            });
        }
        if !(self.validity.reynolds_lo > 0.0
            && self.validity.reynolds_hi >= self.validity.reynolds_lo)
        {
            return Err(AeroacError::InvalidParameter {
                what: "validity region must satisfy 0 < reynolds_lo <= reynolds_hi",
            });
        }
        if !(self.profile.u_centerline > 0.0
            && self.profile.slot_half > 0.0
            && self.profile.momentum_thickness > 0.0)
        {
            return Err(AeroacError::InvalidParameter {
                what: "mean-jet profile parameters must be positive",
            });
        }
        Ok(())
    }

    /// Refuse consumption outside the card's validity region.
    ///
    /// # Errors
    /// [`AeroacError::InvalidParameter`] when `reynolds` is outside
    /// the demonstrated range or non-finite.
    pub fn admit_query(&self, reynolds: f64) -> Result<(), AeroacError> {
        if !reynolds.is_finite() {
            return Err(AeroacError::NonFinite {
                what: "query Reynolds number",
            });
        }
        if !(self.validity.reynolds_lo..=self.validity.reynolds_hi).contains(&reynolds) {
            return Err(AeroacError::InvalidParameter {
                what: "query outside the card's demonstrated Reynolds range (no extrapolation)",
            });
        }
        Ok(())
    }

    /// Deterministic content hash: FNV-1a over the canonical JSON
    /// serialization (field order pinned by [`JetCard::to_json`]).
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        fnv1a(self.to_json().as_bytes())
    }
}

/// Mint the TONAL INTERIM card from the recorded, validated 2-D
/// stage-I edge-tone point (bead clause 3, explicitly permitted).
///
/// Every number here is a recorded lab observation from the executed
/// edge-tone staging battery (fs-aeroac CONTRACT invariants 10/13
/// and `tests/edgetone_staging.rs`): stage-I lock St 0.03662 at
/// Re 144, `h/delta = 10`, vs Brown's (1937) 0.03554 (+3.0%, inside
/// the ±6% bin quantization and the published spread); the 0.02-seed
/// neighboring state St 0.0458; bitwise-reproducible ramp hysteresis
/// over Re 144→264→144. The claim is edge-tone class ONLY.
///
/// # Errors
/// Momentum-thickness refusals (structurally impossible for the
/// pinned staging geometry, but the arithmetic is not duplicated).
pub fn mint_tonal_interim_card() -> Result<JetCard, AeroacError> {
    let cfg = staging_rig_config();
    let theta = momentum_thickness_smoothed_tophat(cfg.slot_half, cfg.slot_smoothing, cfg.ny)?;
    let card = JetCard {
        schema: JET_CARD_SCHEMA.to_owned(),
        claim: JetCardClaim::EdgeToneTonal {
            feedback: EdgeToneFeedback {
                stage: 1,
                locked_strouhal: 0.036_62,
                ladder_band: (0.7, 1.4),
                multi_stable_strouhal: vec![0.036_62, 0.045_8],
                hysteresis_recorded: true,
                saturated_force_rms: None,
            },
        },
        authority: CardAuthority::XStruct,
        profile: MeanJetProfile {
            u_centerline: cfg.u_jet,
            slot_half: cfg.slot_half,
            momentum_thickness: theta,
        },
        band_db: None,
        residual: Some(CardResidual {
            quantity: "stage-I locked Strouhal".to_owned(),
            measured: 0.036_62,
            reference: 0.035_54,
            reference_source: "Brown (1937) stage-I edge-tone law at h/delta = 10 \
                               (Vaik/Paal Part I Table 1 coefficients)"
                .to_owned(),
            bin_halfwidth_rel: 0.06,
        }),
        validity: ValidityRegion {
            reynolds_lo: 144.0,
            reynolds_hi: 264.0,
            h_over_delta: 10.0,
            requires_receptivity_edge: true,
            seeded_amplitude_claims_only: true,
        },
        provenance: CardProvenance {
            rig_fingerprint: jet_labium_fingerprint(&cfg),
            receipts: vec![
                "fs-aeroac/tests/edgetone_staging.rs::edge_tone_stage_one_strouhal_matches_published"
                    .to_owned(),
                "fs-aeroac/CONTRACT.md invariant 10 (edge-tone staging record)".to_owned(),
                "fs-aeroac/CONTRACT.md invariant 13 (adiabatic ramp hysteresis record)".to_owned(),
                "bead frankensim-music-v8-root-3ez8g.10.2".to_owned(),
            ],
            classification: "tonal stage-I edge-tone lock (2-D CentralMoment/BGK landscape; \
                             the 2-D broadband refusal stands)"
                .to_owned(),
        },
        scope: SCOPE_STATEMENT.to_owned(),
    };
    card.validate()?;
    Ok(card)
}

/// The exact staging rig configuration the tonal interim card is
/// minted from (`tests/edgetone_staging.rs`, Re 144, h/delta = 10).
#[must_use]
pub fn staging_rig_config() -> JetLabiumConfig {
    JetLabiumConfig {
        nx: 192,
        ny: 64,
        slot_half: 3.0,
        slot_smoothing: 1.2,
        u_jet: 0.08,
        tau: 0.51,
        edge_distance: 60,
        plate_length: 50,
        fringe_width: 32,
        fringe_sigma: 0.3,
        steps_settle: 4000,
        steps_record: 16_384,
        seed_amplitude: 0.005,
        nozzle_thickness: 2,
        collision: fs_lbm::core2::CollisionModel2::Bgk,
    }
}

/// Regime gates a 3-D rung must pass before it can back any card
/// (an out-of-regime rung would mint a fabricated authority).
/// Bins from the record's temporal-Nyquist edge within which a peak is
/// the parity artifact, not a tone (the 2-D staging battery's law:
/// `bin < n/2 - 8`).
pub const NYQUIST_EDGE_BINS: usize = 8;

/// Recover the FFT record length from a rung's bin-width disclosure and
/// the profile it was measured on: `bin_width = (1/n) 2 slot_half / u`.
fn record_len_of(rung: &SlotJet3dRung, profile: &MeanJetProfile) -> Option<usize> {
    let n = 2.0 * profile.slot_half / (profile.u_centerline * rung.strouhal_bin_width);
    if !n.is_finite() || n < 64.0 {
        return None;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rounded = n.round() as usize;
    ((n - n.round()).abs() < 1e-6 && rounded.is_power_of_two()).then_some(rounded)
}

/// Whether the rung's peak sits at the temporal-Nyquist edge of its
/// record (the LBM bounce-back parity artifact class): such a "tone"
/// is numerical and no card may cite it.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] when the rung's bin-width
/// disclosure and the profile do not reproduce a power-of-two record.
pub fn peak_at_nyquist_edge(
    rung: &SlotJet3dRung,
    profile: &MeanJetProfile,
) -> Result<bool, AeroacError> {
    let n = record_len_of(rung, profile).ok_or(AeroacError::InvalidParameter {
        what: "rung bin-width disclosure does not reproduce a power-of-two record on this profile",
    })?;
    Ok(rung.peak_bin + NYQUIST_EDGE_BINS >= n / 2)
}

fn admit_rung(rung: &SlotJet3dRung, profile: &MeanJetProfile) -> Result<(), AeroacError> {
    if peak_at_nyquist_edge(rung, profile)? {
        return Err(AeroacError::InvalidParameter {
            what: "rung peak sits at the temporal-Nyquist edge (parity artifact, not a tone) and cannot back a card",
        });
    }
    if !rung.amplitude_qualified {
        return Err(AeroacError::InvalidParameter {
            what: "rung below the force-RMS amplitude floor cannot back a card",
        });
    }
    if rung.mach_max_lattice > 0.25 {
        return Err(AeroacError::InvalidParameter {
            what: "rung left the low-Mach regime",
        });
    }
    if !(rung.flux_imbalance.is_finite() && rung.flux_imbalance < 0.05) {
        return Err(AeroacError::InvalidParameter {
            what: "rung flux imbalance too large (outlet-reflection pathology class)",
        });
    }
    Ok(())
}

/// Mint a BROADBAND card from executed, classified 3-D sweep rungs.
///
/// Requires at least one admitted rung classified broadband (the
/// demonstrating rung; the lowest-Re broadband rung is chosen so
/// the validity floor is honest). Band content and the residual are
/// the caller's lab artifacts; without a residual the card is
/// capped at `X-Est` (the bead's authority law).
///
/// # Errors
/// [`AeroacError::InvalidParameter`] when no admitted broadband
/// rung exists, or rung admission fails on the demonstrating rung.
pub fn mint_broadband_card(
    rungs: &[SlotJet3dRung],
    profile: MeanJetProfile,
    band_db: [f64; N_BANDS],
    residual: Option<CardResidual>,
    rig_fingerprint: u64,
    receipts: Vec<String>,
) -> Result<JetCard, AeroacError> {
    let mut demo: Option<&SlotJet3dRung> = None;
    for rung in rungs {
        if !rung.tonal && admit_rung(rung, &profile).is_ok() {
            let better = demo.is_none_or(|d| rung.reynolds < d.reynolds);
            if better {
                demo = Some(rung);
            }
        }
    }
    let Some(rung) = demo else {
        return Err(AeroacError::InvalidParameter {
            what: "no admitted broadband rung: a broadband card cannot be minted \
                   (the refusal-boundary card is the honest artifact)",
        });
    };
    let authority = if residual.is_some() {
        CardAuthority::XStruct
    } else {
        CardAuthority::XEst
    };
    let hi = rungs
        .iter()
        .filter(|r| !r.tonal && admit_rung(r, &profile).is_ok())
        .map(|r| r.reynolds)
        .fold(rung.reynolds, f64::max);
    let card = JetCard {
        schema: JET_CARD_SCHEMA.to_owned(),
        claim: JetCardClaim::Broadband {
            flatness: rung.flatness,
            reynolds: rung.reynolds,
        },
        authority,
        profile,
        band_db: Some(band_db),
        residual,
        validity: ValidityRegion {
            reynolds_lo: rung.reynolds,
            reynolds_hi: hi,
            h_over_delta: 0.0,
            requires_receptivity_edge: false,
            seeded_amplitude_claims_only: true,
        },
        provenance: CardProvenance {
            rig_fingerprint,
            receipts,
            classification: "broadband (3-D slot-jet sweep classification)".to_owned(),
        },
        scope: SCOPE_STATEMENT.to_owned(),
    };
    card.validate()?;
    Ok(card)
}

/// Mint the REFUSAL-BOUNDARY card from an executed sweep in which
/// every admitted rung classified tonal — the honest "no broadband
/// below Re X on this rig family" artifact (bead clause 4).
///
/// # Errors
/// [`AeroacError::InvalidParameter`] when fewer than two admitted
/// rungs exist (a single point is not a boundary) or when any
/// admitted rung is broadband (the broadband card is then the
/// truthful mint, not this one).
pub fn mint_refusal_boundary_card(
    rungs: &[SlotJet3dRung],
    profile: MeanJetProfile,
    rig_fingerprint: u64,
    receipts: Vec<String>,
) -> Result<JetCard, AeroacError> {
    let admitted: Vec<&SlotJet3dRung> = rungs
        .iter()
        .filter(|r| admit_rung(r, &profile).is_ok())
        .collect();
    if admitted.len() < 2 {
        return Err(AeroacError::InvalidParameter {
            what: "a refusal boundary needs at least two admitted rungs",
        });
    }
    if admitted.iter().any(|r| !r.tonal) {
        return Err(AeroacError::InvalidParameter {
            what: "an admitted broadband rung exists: mint the broadband card instead",
        });
    }
    let lo = admitted
        .iter()
        .map(|r| r.reynolds)
        .fold(f64::INFINITY, f64::min);
    let hi = admitted
        .iter()
        .map(|r| r.reynolds)
        .fold(f64::NEG_INFINITY, f64::max);
    let card = JetCard {
        schema: JET_CARD_SCHEMA.to_owned(),
        claim: JetCardClaim::BroadbandRefusalBoundary {
            max_reynolds_probed: hi,
            rungs_probed: admitted.len(),
        },
        authority: CardAuthority::XEst,
        profile,
        band_db: None,
        residual: None,
        validity: ValidityRegion {
            reynolds_lo: lo,
            reynolds_hi: hi,
            h_over_delta: 0.0,
            requires_receptivity_edge: false,
            seeded_amplitude_claims_only: true,
        },
        provenance: CardProvenance {
            rig_fingerprint,
            receipts,
            classification: "tonal at every admitted rung (broadband refused below the boundary)"
                .to_owned(),
        },
        scope: SCOPE_STATEMENT.to_owned(),
    };
    card.validate()?;
    Ok(card)
}

// ------------------------------------------------------------------
// Serialization: pinned-field-order JSON writer + strict fail-closed
// parser. The parser accepts EXACTLY the writer's shape (schema v1);
// any other bytes refuse by name — no silent reinterpretation.
// ------------------------------------------------------------------

fn esc(s: &str, out: &mut String) {
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
}

#[allow(clippy::format_push_string)] // export builder, clarity over micro-alloc
impl JetCard {
    /// Canonical JSON export (field order pinned; the content hash
    /// is defined over these exact bytes).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\"schema\":\"");
        esc(&self.schema, &mut s);
        s.push_str("\",\"claim_kind\":\"");
        s.push_str(self.claim.kind());
        s.push_str("\",\"claim\":{");
        match &self.claim {
            JetCardClaim::EdgeToneTonal { feedback } => {
                s.push_str(&format!(
                    "\"stage\":{},\"locked_strouhal\":{},\"ladder_band_lo\":{},\
                     \"ladder_band_hi\":{},\"multi_stable_strouhal\":[",
                    feedback.stage,
                    feedback.locked_strouhal,
                    feedback.ladder_band.0,
                    feedback.ladder_band.1
                ));
                for (i, v) in feedback.multi_stable_strouhal.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    s.push_str(&format!("{v}"));
                }
                s.push_str(&format!(
                    "],\"hysteresis_recorded\":{}",
                    feedback.hysteresis_recorded
                ));
                s.push_str(",\"saturated_force_rms\":");
                match feedback.saturated_force_rms {
                    None => s.push_str("null"),
                    Some(v) => s.push_str(&format!("{v}")),
                }
            }
            JetCardClaim::Broadband { flatness, reynolds } => {
                s.push_str(&format!("\"flatness\":{flatness},\"reynolds\":{reynolds}"));
            }
            JetCardClaim::BroadbandRefusalBoundary {
                max_reynolds_probed,
                rungs_probed,
            } => {
                s.push_str(&format!(
                    "\"max_reynolds_probed\":{max_reynolds_probed},\"rungs_probed\":{rungs_probed}"
                ));
            }
        }
        s.push_str("},\"authority\":\"");
        s.push_str(self.authority.label());
        s.push_str("\",\"profile\":{");
        s.push_str(&format!(
            "\"u_centerline\":{},\"slot_half\":{},\"momentum_thickness\":{}",
            self.profile.u_centerline, self.profile.slot_half, self.profile.momentum_thickness
        ));
        s.push_str("},\"band_db\":");
        match &self.band_db {
            None => s.push_str("null"),
            Some(bands) => {
                s.push('[');
                for (i, v) in bands.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    s.push_str(&format!("{v}"));
                }
                s.push(']');
            }
        }
        s.push_str(",\"residual\":");
        match &self.residual {
            None => s.push_str("null"),
            Some(r) => {
                s.push_str("{\"quantity\":\"");
                esc(&r.quantity, &mut s);
                s.push_str(&format!(
                    "\",\"measured\":{},\"reference\":{},\"reference_source\":\"",
                    r.measured, r.reference
                ));
                esc(&r.reference_source, &mut s);
                s.push_str(&format!(
                    "\",\"bin_halfwidth_rel\":{}}}",
                    r.bin_halfwidth_rel
                ));
            }
        }
        s.push_str(&format!(
            ",\"validity\":{{\"reynolds_lo\":{},\"reynolds_hi\":{},\"h_over_delta\":{},\
             \"requires_receptivity_edge\":{},\"seeded_amplitude_claims_only\":{}}}",
            self.validity.reynolds_lo,
            self.validity.reynolds_hi,
            self.validity.h_over_delta,
            self.validity.requires_receptivity_edge,
            self.validity.seeded_amplitude_claims_only
        ));
        s.push_str(&format!(
            ",\"provenance\":{{\"rig_fingerprint\":{},\"receipts\":[",
            self.provenance.rig_fingerprint
        ));
        for (i, r) in self.provenance.receipts.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push('"');
            esc(r, &mut s);
            s.push('"');
        }
        s.push_str("],\"classification\":\"");
        esc(&self.provenance.classification, &mut s);
        s.push_str("\"},\"scope\":\"");
        esc(&self.scope, &mut s);
        s.push_str("\"}");
        s
    }

    /// Parse a card serialized by [`JetCard::to_json`]. Strict:
    /// exactly the pinned v1 shape, refusing anything else by name.
    /// The parsed card is re-validated before it is returned, so a
    /// byte stream cannot smuggle an invariant violation through
    /// deserialization.
    ///
    /// # Errors
    /// [`AeroacError::InvalidParameter`] naming the first violated
    /// expectation.
    pub fn from_json(text: &str) -> Result<Self, AeroacError> {
        let mut p = Parser { s: text, pos: 0 };
        p.lit("{\"schema\":\"")?;
        let schema = p.string()?;
        if schema != JET_CARD_SCHEMA {
            return Err(AeroacError::InvalidParameter {
                what: "unknown jet-card schema (refused by name)",
            });
        }
        p.lit("\"claim_kind\":\"")?;
        let kind = p.string()?;
        p.lit("\"claim\":{")?;
        let claim = match kind.as_str() {
            "edge-tone-tonal" => {
                p.lit("\"stage\":")?;
                let stage = p.number()?;
                p.lit("\"locked_strouhal\":")?;
                let locked = p.number()?;
                p.lit("\"ladder_band_lo\":")?;
                let lo = p.number()?;
                p.lit("\"ladder_band_hi\":")?;
                let hi = p.number()?;
                p.lit("\"multi_stable_strouhal\":[")?;
                let states = p.number_list(']')?;
                p.lit(",\"hysteresis_recorded\":")?;
                let hyst = p.boolean()?;
                p.lit(",\"saturated_force_rms\":")?;
                let sat = if p.peek_null() {
                    p.lit("null")?;
                    None
                } else {
                    Some(p.number()?)
                };
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let stage_u8 = stage as u8;
                JetCardClaim::EdgeToneTonal {
                    feedback: EdgeToneFeedback {
                        stage: stage_u8,
                        locked_strouhal: locked,
                        ladder_band: (lo, hi),
                        multi_stable_strouhal: states,
                        hysteresis_recorded: hyst,
                        saturated_force_rms: sat,
                    },
                }
            }
            "broadband" => {
                p.lit("\"flatness\":")?;
                let flatness = p.number()?;
                p.lit("\"reynolds\":")?;
                let reynolds = p.number()?;
                JetCardClaim::Broadband { flatness, reynolds }
            }
            "broadband-refusal-boundary" => {
                p.lit("\"max_reynolds_probed\":")?;
                let max_re = p.number()?;
                p.lit("\"rungs_probed\":")?;
                let n = p.number()?;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let n_rungs = n as usize;
                JetCardClaim::BroadbandRefusalBoundary {
                    max_reynolds_probed: max_re,
                    rungs_probed: n_rungs,
                }
            }
            _ => {
                return Err(AeroacError::InvalidParameter {
                    what: "unknown jet-card claim kind (refused by name)",
                });
            }
        };
        p.lit("},\"authority\":\"")?;
        let auth = p.string()?;
        let authority = match auth.as_str() {
            "X-Struct" => CardAuthority::XStruct,
            "X-Est" => CardAuthority::XEst,
            _ => {
                return Err(AeroacError::InvalidParameter {
                    what: "unknown jet-card authority class",
                });
            }
        };
        p.lit("\"profile\":{\"u_centerline\":")?;
        let u_centerline = p.number()?;
        p.lit("\"slot_half\":")?;
        let slot_half = p.number()?;
        p.lit("\"momentum_thickness\":")?;
        let theta = p.number()?;
        p.lit("},\"band_db\":")?;
        let band_db = if p.peek_null() {
            p.lit("null")?;
            None
        } else {
            p.lit("[")?;
            let list = p.number_list(']')?;
            let arr: [f64; N_BANDS] =
                list.try_into().map_err(|_| AeroacError::InvalidParameter {
                    what: "band_db must have exactly N_BANDS entries",
                })?;
            Some(arr)
        };
        p.lit(",\"residual\":")?;
        let residual = if p.peek_null() {
            p.lit("null")?;
            None
        } else {
            p.lit("{\"quantity\":\"")?;
            let quantity = p.string()?;
            p.lit("\"measured\":")?;
            let measured = p.number()?;
            p.lit("\"reference\":")?;
            let reference = p.number()?;
            p.lit("\"reference_source\":\"")?;
            let source = p.string()?;
            p.lit("\"bin_halfwidth_rel\":")?;
            let bin = p.number()?;
            p.lit("}")?;
            Some(CardResidual {
                quantity,
                measured,
                reference,
                reference_source: source,
                bin_halfwidth_rel: bin,
            })
        };
        p.lit(",\"validity\":{\"reynolds_lo\":")?;
        let re_lo = p.number()?;
        p.lit("\"reynolds_hi\":")?;
        let re_hi = p.number()?;
        p.lit("\"h_over_delta\":")?;
        let h_over_delta = p.number()?;
        p.lit("\"requires_receptivity_edge\":")?;
        let edge = p.boolean()?;
        p.lit(",\"seeded_amplitude_claims_only\":")?;
        let seeded = p.boolean()?;
        p.lit("},\"provenance\":{\"rig_fingerprint\":")?;
        let fingerprint = p.integer()?;
        p.lit("\"receipts\":[")?;
        let receipts = p.string_list(']')?;
        p.lit(",\"classification\":\"")?;
        let classification = p.string()?;
        p.lit("},\"scope\":\"")?;
        let scope = p.string()?;
        p.lit("}")?;
        if p.pos != p.s.len() {
            return Err(AeroacError::InvalidParameter {
                what: "trailing bytes after the jet-card document",
            });
        }
        let card = JetCard {
            schema,
            claim,
            authority,
            profile: MeanJetProfile {
                u_centerline,
                slot_half,
                momentum_thickness: theta,
            },
            band_db,
            residual,
            validity: ValidityRegion {
                reynolds_lo: re_lo,
                reynolds_hi: re_hi,
                h_over_delta,
                requires_receptivity_edge: edge,
                seeded_amplitude_claims_only: seeded,
            },
            provenance: CardProvenance {
                rig_fingerprint: fingerprint,
                receipts,
                classification,
            },
            scope,
        };
        card.validate()?;
        Ok(card)
    }
}

/// Minimal strict cursor over a pinned document shape. Shared with
/// the sweep-receipt reader in [`crate::slot_jet_3d`]: both formats are
/// writer-pinned field orders, and a strict cursor is exactly the
/// fail-closed reader they want.
pub(crate) struct Parser<'a> {
    s: &'a str,
    pos: usize,
}

impl<'a> Parser<'a> {
    pub(crate) fn new(s: &'a str) -> Self {
        Self { s, pos: 0 }
    }

    fn refuse(what: &'static str) -> AeroacError {
        AeroacError::InvalidParameter { what }
    }

    /// Whether the cursor has consumed the whole input (trailing
    /// bytes after a document are a refusal for line-framed receipts).
    pub(crate) fn at_end(&self) -> bool {
        self.pos == self.s.len()
    }

    /// A number token that may be non-finite (`inf`/`NaN` as the
    /// writer's `Display` emits them). Only for fields whose
    /// non-finite value is itself a recorded fact, never for
    /// quantities a downstream gate consumes without re-checking.
    pub(crate) fn number_allow_nonfinite(&mut self) -> Result<f64, AeroacError> {
        let tok = self.number_token()?;
        tok.parse()
            .map_err(|_| Self::refuse("jet-card parse: malformed number"))
    }

    /// Expect the literal `lit` at the cursor. Number/bool fields
    /// are separated by `,`/`}` consumed by the value readers, so
    /// literals starting with a field key also swallow a single
    /// leading separator if present.
    pub(crate) fn lit(&mut self, lit: &str) -> Result<(), AeroacError> {
        let rest = &self.s[self.pos..];
        let rest = match rest.strip_prefix(',') {
            Some(stripped)
                if !(lit.starts_with(',') || lit.starts_with('}') || lit.starts_with(']')) =>
            {
                stripped
            }
            _ => rest,
        };
        if let Some(after) = rest.strip_prefix(lit) {
            self.pos = self.s.len() - after.len();
            Ok(())
        } else {
            Err(Self::refuse("jet-card parse: pinned field order violated"))
        }
    }

    fn peek_null(&self) -> bool {
        self.s[self.pos..].starts_with("null") || self.s[self.pos..].starts_with(",null")
    }

    /// A quoted string body up to the closing quote, unescaping the
    /// writer's escapes; consumes the closing `"` and a following
    /// `,` when present.
    pub(crate) fn string(&mut self) -> Result<String, AeroacError> {
        let bytes = self.s.as_bytes();
        let mut out = String::new();
        let mut i = self.pos;
        while i < bytes.len() {
            match bytes[i] {
                b'\\' => {
                    let next = bytes
                        .get(i + 1)
                        .ok_or_else(|| Self::refuse("jet-card parse: dangling escape"))?;
                    match next {
                        b'\\' => out.push('\\'),
                        b'"' => out.push('"'),
                        b'n' => out.push('\n'),
                        _ => return Err(Self::refuse("jet-card parse: unknown escape")),
                    }
                    i += 2;
                }
                b'"' => {
                    self.pos = i + 1;
                    if self.s[self.pos..].starts_with(',') {
                        self.pos += 1;
                    }
                    return Ok(out);
                }
                _ => {
                    let ch_end = self.s[i..]
                        .char_indices()
                        .nth(1)
                        .map_or(self.s.len(), |(o, _)| i + o);
                    out.push_str(&self.s[i..ch_end]);
                    i = ch_end;
                }
            }
        }
        Err(Self::refuse("jet-card parse: unterminated string"))
    }

    /// A number token ending at `,`, `}`, or `]` (delimiter kept for
    /// `lit` to consume, except a plain `,` which is swallowed).
    fn number_token(&mut self) -> Result<&str, AeroacError> {
        let rest = &self.s[self.pos..];
        let end = rest
            .find([',', '}', ']'])
            .ok_or_else(|| Self::refuse("jet-card parse: unterminated number"))?;
        let tok = &rest[..end];
        self.pos += end;
        Ok(tok)
    }

    pub(crate) fn number(&mut self) -> Result<f64, AeroacError> {
        let tok = self.number_token()?;
        let v: f64 = tok
            .parse()
            .map_err(|_| Self::refuse("jet-card parse: malformed number"))?;
        if !v.is_finite() {
            return Err(Self::refuse("jet-card parse: non-finite number"));
        }
        Ok(v)
    }

    pub(crate) fn integer(&mut self) -> Result<u64, AeroacError> {
        let tok = self.number_token()?;
        tok.parse()
            .map_err(|_| Self::refuse("jet-card parse: malformed integer"))
    }

    pub(crate) fn boolean(&mut self) -> Result<bool, AeroacError> {
        if self.lit("true").is_ok() {
            return Ok(true);
        }
        if self.lit("false").is_ok() {
            return Ok(false);
        }
        Err(Self::refuse("jet-card parse: malformed boolean"))
    }

    /// A `[..]` list of numbers whose opening bracket was already
    /// consumed; consumes the closing bracket.
    fn number_list(&mut self, close: char) -> Result<Vec<f64>, AeroacError> {
        let mut out = Vec::new();
        loop {
            if self.s[self.pos..].starts_with(close) {
                self.pos += 1;
                return Ok(out);
            }
            out.push(self.number()?);
            if self.s[self.pos..].starts_with(',') {
                self.pos += 1;
            }
        }
    }

    /// A `[..]` list of quoted strings whose opening bracket was
    /// already consumed; consumes the closing bracket.
    fn string_list(&mut self, close: char) -> Result<Vec<String>, AeroacError> {
        let mut out = Vec::new();
        loop {
            if self.s[self.pos..].starts_with(close) {
                self.pos += 1;
                return Ok(out);
            }
            if !self.s[self.pos..].starts_with('"') {
                return Err(Self::refuse("jet-card parse: expected quoted receipt"));
            }
            self.pos += 1;
            out.push(self.string()?);
        }
    }
}
