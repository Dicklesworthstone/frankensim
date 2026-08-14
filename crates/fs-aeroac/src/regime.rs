//! Broadband-regime admission and the recorded 2D tonal-lock refusal
//! (bead frankensim-l011o).
//!
//! The executed CentralMoment 2D jet-labium sweep (2026-08-10) stays an
//! essentially pure tone at every reachable Reynolds number. This module
//! pins that landscape as data, refuses to catalog those spectra as
//! broadband flute noise, and names the 3D slot-jet follow-up that would
//! be required to reopen the claim.

use crate::AeroacError;
use crate::SCOPE_STATEMENT;
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};

/// Mechanism recorded for the 2D tonal lock.
pub const TWO_D_INVERSE_CASCADE: &str = "2D Navier-Stokes inverse-cascades kinetic energy \
to large scales, so a planar slot-jet/labium pair organizes into a coherent \
limit-cycle oscillation rather than a broadband cascade; this is a 2D physics \
boundary, not an under-resolved 3D jet.";

/// Spectral-flatness ceiling below which a spectrum is classified tonal.
///
/// Executed CentralMoment rungs reported flatness ~1e-18. A decade of
/// margin keeps the pin from flickering on rounding while still refusing
/// anything that is not a near-line spectrum.
pub const TONAL_FLATNESS_CEILING: f64 = 1.0e-6;

/// One pinned 2D probe row from the executed CentralMoment sweep.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TonalLockRow {
    /// Jet Reynolds number of the executed run.
    pub reynolds: f64,
    /// Slot half-width in lattice units.
    pub slot_half: f64,
    /// Locked Strouhal number (f * delta / u).
    pub strouhal: f64,
    /// Upper bound on measured spectral flatness (geometric/arithmetic
    /// power-mean ratio; 1 = white).
    pub flatness_upper: f64,
    /// Whether the run completed inside the low-Mach contract.
    pub ran_in_regime: bool,
}

/// Executed CentralMoment 2D landscape (bead l011o comment, 2026-08-10).
///
/// These are recorded observations, not a claim that every future operator
/// change will reprint the same digits. A later probe that finds flatness
/// above [`TONAL_FLATNESS_CEILING`] at these (Re, delta) must reopen the
/// bead rather than silently overwrite this table.
pub const PINNED_2D_CENTRAL_MOMENT_TONAL: &[TonalLockRow] = &[
    TonalLockRow {
        reynolds: 432.0,
        slot_half: 6.0,
        strouhal: 0.092,
        flatness_upper: 1.0e-15,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 576.0,
        slot_half: 6.0,
        strouhal: 0.101,
        flatness_upper: 1.0e-15,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 1_152.0,
        slot_half: 6.0,
        strouhal: 0.467,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 2_304.0,
        slot_half: 6.0,
        strouhal: 0.467,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 4_608.0,
        slot_half: 6.0,
        strouhal: 0.467,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 9_216.0,
        slot_half: 6.0,
        strouhal: 0.467,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 2_304.0,
        slot_half: 12.0,
        strouhal: 0.458,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
    TonalLockRow {
        reynolds: 9_216.0,
        slot_half: 12.0,
        strouhal: 0.476,
        flatness_upper: 1.0e-18,
        ran_in_regime: true,
    },
];

/// Typed 2D broadband refusal: the sweep exists, every row is tonal,
/// and the mechanism is named.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoDBroadbandRefusal {
    /// Executed rows that close the 2D hunt.
    pub rows: &'static [TonalLockRow],
    /// Physics mechanism, not a numerical excuse.
    pub mechanism: &'static str,
    /// Scope law copied onto the refusal so a consumer cannot drop it.
    pub scope: &'static str,
}

/// The recorded 2D refusal.
#[must_use]
pub const fn two_d_broadband_refusal() -> TwoDBroadbandRefusal {
    TwoDBroadbandRefusal {
        rows: PINNED_2D_CENTRAL_MOMENT_TONAL,
        mechanism: TWO_D_INVERSE_CASCADE,
        scope: SCOPE_STATEMENT,
    }
}

/// Geometric-to-arithmetic mean of a one-sided power spectrum.
///
/// White noise → 1. A pure tone → 0. This is the measurement the
/// 3-D jet hunt must pass; it does not mint a broadband table.
///
/// # Errors
/// [`AeroacError`] if fewer than two finite positive bins remain.
pub fn measure_spectral_flatness(power: &[f64]) -> Result<f64, AeroacError> {
    let mut n = 0.0;
    let mut sum = 0.0;
    let mut log_sum = 0.0;
    for &p in power {
        if !p.is_finite() {
            return Err(AeroacError::NonFinite {
                what: "spectral power",
            });
        }
        if p > 0.0 {
            n += 1.0;
            sum += p;
            log_sum += p.ln();
        }
    }
    if n < 2.0 || !(sum > 0.0) {
        return Err(AeroacError::InvalidParameter {
            what: "spectral flatness needs at least two positive power bins",
        });
    }
    let arith = sum / n;
    Ok((log_sum / n).exp() / arith)
}

/// A measured spectrum classified against [`TONAL_FLATNESS_CEILING`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectrumClass {
    /// Flatness at or below the tonal ceiling.
    Tonal {
        /// Measured geometric/arithmetic ratio.
        flatness: f64,
    },
    /// Flatness above the tonal ceiling.
    Broadband {
        /// Measured geometric/arithmetic ratio.
        flatness: f64,
    },
}

/// Classify a power spectrum without cataloging it as flute noise.
///
/// # Errors
/// Measurement refusals from [`measure_spectral_flatness`].
pub fn classify_spectrum(power: &[f64]) -> Result<SpectrumClass, AeroacError> {
    let flatness = measure_spectral_flatness(power)?;
    if flatness <= TONAL_FLATNESS_CEILING {
        Ok(SpectrumClass::Tonal { flatness })
    } else {
        Ok(SpectrumClass::Broadband { flatness })
    }
}

/// Refuse to treat a spectrum as broadband flute noise.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] when `flatness` is non-finite or
/// at or below [`TONAL_FLATNESS_CEILING`].
pub fn admit_broadband_spectrum(flatness: f64) -> Result<(), AeroacError> {
    if !flatness.is_finite() {
        return Err(AeroacError::NonFinite {
            what: "spectral flatness",
        });
    }
    if flatness <= TONAL_FLATNESS_CEILING {
        return Err(AeroacError::InvalidParameter {
            what: "spectrum is tonal (flatness at or below the broadband floor); \
                   cataloging it as flute-noise would launder a limit cycle",
        });
    }
    Ok(())
}

/// Geometry and operator contract for a future 3D slot-jet broadband hunt.
///
/// This is the concrete follow-up named by l011o: 3D is where a jet
/// cascade can exist. It does not claim a demonstrated 3D broadband
/// regime; constructing and validating the spec is the evaluation of
/// the D3Q19 path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlotJet3dFollowUp {
    /// Lattice sites in x (must be a multiple of the D3Q19 tile size 4).
    pub nx: u32,
    /// Lattice sites in y.
    pub ny: u32,
    /// Lattice sites in z (span; broadband lives in the third direction).
    pub nz: u32,
    /// Slot half-width [lu].
    pub slot_half: f64,
    /// Collision law. Central-moment is the unlocked 3D operator.
    pub collision: CollisionModel3,
    /// Flatness that would reopen a broadband claim.
    pub broadband_flatness_floor: f64,
}

impl SlotJet3dFollowUp {
    /// Minimal tile-aligned slot-jet specification using D3Q19
    /// [`CollisionModel3::CentralMoment`].
    #[must_use]
    pub const fn minimal_central_moment() -> Self {
        Self {
            nx: 32,
            ny: 16,
            nz: 16,
            slot_half: 2.0,
            collision: CollisionModel3::CentralMoment {
                second_order_rate: 1.5,
                higher_order_rate: 1.5,
            },
            broadband_flatness_floor: TONAL_FLATNESS_CEILING,
        }
    }

    /// Check tile alignment, physical slot, and collision-parameter window.
    ///
    /// # Errors
    /// [`AeroacError`] when the spec cannot be executed even as a smoke.
    pub fn validate(self) -> Result<(), AeroacError> {
        if self.nx == 0
            || self.ny == 0
            || self.nz == 0
            || self.nx % 4 != 0
            || self.ny % 4 != 0
            || self.nz % 4 != 0
        {
            return Err(AeroacError::InvalidParameter {
                what: "3D slot-jet extents must be positive multiples of the D3Q19 tile size 4",
            });
        }
        if !(self.slot_half > 0.0 && self.slot_half.is_finite()) {
            return Err(AeroacError::InvalidParameter {
                what: "3D slot half-width must be positive and finite",
            });
        }
        if !(self.broadband_flatness_floor > 0.0 && self.broadband_flatness_floor.is_finite()) {
            return Err(AeroacError::InvalidParameter {
                what: "broadband flatness floor must be positive and finite",
            });
        }
        self.collision
            .validate()
            .map_err(|_| AeroacError::InvalidParameter {
                what: "3D collision relaxation is outside the D3Q19 physical window",
            })?;
        Ok(())
    }
}

/// Result of evaluating the 3D operator path without claiming a regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SlotJet3dOperatorSmoke {
    /// The D3Q19 central-moment law admits the follow-up spec.
    pub collision_admits: bool,
    /// A rest-frame equilibrium is a positive, finite probability set.
    pub equilibrium_live: bool,
    /// Whether a broadband 3D regime was demonstrated (always false here).
    pub broadband_demonstrated: bool,
}

/// Smoke-evaluate the D3Q19 central-moment path named by the follow-up.
///
/// This is the l011o (b) evaluation: the operator exists and is
/// executable. It is not a Reynolds sweep and it does not mint a
/// broadband table.
///
/// # Errors
/// Spec validation failures.
pub fn evaluate_slot_jet_3d_operator(
    spec: SlotJet3dFollowUp,
) -> Result<SlotJet3dOperatorSmoke, AeroacError> {
    spec.validate()?;
    let f = equilibrium3(1.0, [0.05, 0.0, 0.0]);
    let mut mass = 0.0;
    let mut all_finite_positive = true;
    for &fi in &f {
        if !(fi.is_finite() && fi > 0.0) {
            all_finite_positive = false;
        }
        mass += fi;
    }
    if !((mass - 1.0).abs() < 1.0e-12) {
        all_finite_positive = false;
    }
    Ok(SlotJet3dOperatorSmoke {
        collision_admits: true,
        equilibrium_live: all_finite_positive,
        broadband_demonstrated: false,
    })
}

impl SlotJet3dOperatorSmoke {
    /// Attach a *measured* power spectrum. `broadband_demonstrated`
    /// becomes true only if the spectrum is admitted above both the
    /// tonal ceiling and the follow-up floor. A missing 3-D run
    /// cannot call this; a white-noise fixture can, and that is the
    /// honest measurement primitive, not a minted jet table.
    ///
    /// # Errors
    /// [`measure_spectral_flatness`] refusals.
    pub fn incorporate_measured_spectrum(
        self,
        power: &[f64],
        floor: f64,
    ) -> Result<Self, AeroacError> {
        let flatness = measure_spectral_flatness(power)?;
        let demonstrated = admit_broadband_spectrum(flatness).is_ok() && flatness > floor;
        Ok(Self {
            broadband_demonstrated: demonstrated,
            ..self
        })
    }
}
