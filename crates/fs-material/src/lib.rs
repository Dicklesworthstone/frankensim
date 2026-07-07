//! fs-material (plan patch Rev E): the constitutive-law kernel.
//! Materials are NOT constants — they are mathematical objects with
//! calibration domains, CONSISTENT tangent operators (algorithmic,
//! matching the stress update exactly, or Newton dies), thermodynamic
//! guardrails, hysteresis, and uncertainty. Owning them in one crate is
//! what makes structural claims credible.
//!
//! Layer: L3 (FLUX support). Runtime deps: `std`, fs-ad (duals for
//! energy-derived stresses/tangents), fs-evidence (model cards +
//! Evidence), fs-qty, fs-math.

pub mod calibrate;
pub mod elastic;
pub mod fiber;
pub mod hyper;
pub mod plastic;
pub mod tensor;

pub use calibrate::{CalibrationFit, calibrate_bilinear};
pub use elastic::{IsotropicElastic, OrthotropicElastic};
pub use fiber::{ManderConcrete, MenegottoPintoSteel, Uniaxial};
pub use hyper::{Hyperelastic, HyperelasticModel};
pub use plastic::{J2Plasticity, J2State};

use core::fmt;
use fs_evidence::ModelCard;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Voigt-notation symmetric tensor: [xx, yy, zz, xy, yz, zx] with TENSOR
/// (not engineering) shear components.
pub type Voigt = [f64; 6];

/// A 6×6 tangent operator in Voigt layout.
pub type Tangent6 = [[f64; 6]; 6];

/// Thermodynamic admissibility declarations — each law states what it
/// guarantees so checks can be machine-run where possible.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialAdmissibility {
    /// The law derives from a stored-energy function.
    pub has_stored_energy: bool,
    /// Dissipation is non-negative along every admissible path.
    pub dissipation_nonnegative: bool,
    /// The energy is polyconvex (hyperelastic laws; None = unknown).
    pub polyconvex: Option<bool>,
    /// The consistent tangent is symmetric (associative flow/elastic).
    pub tangent_symmetric: bool,
    /// Human-readable failure-envelope statement.
    pub failure_envelope: &'static str,
}

/// A small-strain constitutive law (Voigt 6-space) with internal state.
/// The TANGENT CONTRACT: `tangent` must be the exact derivative of the
/// stress returned by `stress` at the SAME converged state update — the
/// conformance suite enforces this against finite differences for every
/// law (merge-gate discipline).
pub trait SmallStrainLaw {
    /// Internal-variable state (empty for elasticity).
    type State: Clone + PartialEq + fmt::Debug;

    /// The virgin state.
    fn initial_state(&self) -> Self::State;

    /// Stress at `strain`, given the state at the LAST converged step
    /// (the update is algorithmic: predictor/corrector inside).
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt;

    /// The algorithmically consistent tangent dσ/dε at the same update.
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6;

    /// Commit the state update for `strain`.
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State;

    /// The law's thermodynamic declarations.
    fn admissibility(&self) -> MaterialAdmissibility;

    /// The law's model card (assumptions, calibration domain, failures).
    fn card(&self) -> ModelCard;
}

/// Evaluate a law and wrap the stress in `Evidence`: exact numerics (the
/// update is deterministic arithmetic), the law's model card attached,
/// and `in_domain` reflecting whether the strain sits inside the card's
/// CALIBRATION DOMAIN — outside-domain use is flagged, not refused.
#[must_use]
pub fn evidence_stress<L: SmallStrainLaw>(
    law: &L,
    strain: &Voigt,
    state: &L::State,
) -> fs_evidence::Evidence<Voigt> {
    use std::collections::BTreeMap;
    let stress = law.stress(strain, state);
    let qoi = tensor::von_mises(&stress);
    let card = law.card();
    let magnitude = strain.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let mut point = BTreeMap::new();
    point.insert("strain-magnitude".to_string(), magnitude);
    let mut canon = format!("stress:{}", card.name);
    for v in strain {
        let _ = core::fmt::Write::write_fmt(&mut canon, format_args!(";{v}"));
    }
    fs_evidence::Evidence {
        value: stress,
        qoi,
        numerical: fs_evidence::NumericalCertificate::exact(qoi),
        statistical: fs_evidence::StatisticalCertificate::None,
        model: fs_evidence::ModelEvidence::from_card(&card, &point),
        sensitivity: fs_evidence::SensitivitySummary::default(),
        provenance: fs_evidence::ProvenanceHash::of_bytes(canon.as_bytes()),
        adjoint_ref: None,
    }
}

/// Structured material failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialError {
    /// Parameters outside physical bounds (e.g. ν ≥ 0.5, negative modulus).
    Parameters {
        /// Diagnosis.
        what: String,
    },
    /// A deformation state the law cannot evaluate (e.g. det F ≤ 0).
    State {
        /// Diagnosis.
        what: String,
    },
    /// Calibration failure (insufficient/degenerate data).
    Calibration {
        /// Diagnosis.
        what: String,
    },
}

impl fmt::Display for MaterialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaterialError::Parameters { what } => write!(f, "bad material parameters: {what}"),
            MaterialError::State { what } => write!(f, "inadmissible state: {what}"),
            MaterialError::Calibration { what } => write!(f, "calibration failed: {what}"),
        }
    }
}

impl std::error::Error for MaterialError {}
