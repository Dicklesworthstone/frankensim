//! fs-material: constitutive laws and authority-separated identifiability.
//! Materials are NOT constants — they are mathematical objects with
//! calibration domains, CONSISTENT tangent operators (algorithmic,
//! matching the stress update exactly, or Newton dies), thermodynamic
//! guardrails, hysteresis, and uncertainty. Owning them in one crate is
//! what makes structural claims credible.
//!
//! The [`identifiability`] module separately closes multi-case physical
//! questions, source admission, execution plans, and evidence-bearing
//! assessments without conflating a content digest with a scientific theorem.
//!
//! Layer: L3 (FLUX support). Runtime deps: `std`, fs-ad (duals for
//! energy-derived stresses/tangents), fs-blake3 (typed content identities),
//! fs-evidence (model cards, Evidence, and V&V artifacts), fs-matdb (immutable
//! material/model cards), fs-qty, fs-math, fs-thermal (shared transport),
//! and fs-thermochem (species thermodynamics).

/// Prescribed active fibre stress on top of the passive hyperelastic laws (M0.3-W3).
pub mod active_tissue;
/// Finite gas/solid adsorption with sourced equilibrium, kinetics and calorics.
pub mod adsorption;
/// Source-bound monotonic axial comparisons from reported tensile/compression moduli.
pub mod axial_reference;
/// Bounded compression references from original nominal stress/strain observations.
pub mod compression_reference;
pub mod calibrate;
/// Source-bound liquid/vapor properties shared by cavitation consumers.
pub mod cavitation;
/// Irreversible mixed-mode cohesive interfaces with work and fracture accounting.
pub mod cohesive;
pub mod conductor;
/// Permanent sub-yield Norton/Arrhenius creep and consistent physical-step updates.
pub mod creep;
/// Ordered source-backed Ohmic contact resistance, distinct from bulk resistivity.
pub mod electrical_interface;
pub mod liquid;
pub mod liquid_temperature;
pub mod elastic;
/// Bounded streaming rainflow and explicitly scoped empirical stress-life damage.
pub mod fatigue;
pub mod fiber;
/// Source-derived dielectric and magnetic storage for shared circuit consumers.
pub mod field_storage;
/// Shared gas/liquid transport for material-driven flow and heat consumers.
pub mod fluid_transport;
/// Shared checked adiabatic acoustic response for gases and liquids.
pub mod fluid_acoustics;
pub mod gas;
/// Prescribed axial temperature/composition fields on the shared gas laws.
pub mod gas_field;
pub mod gas_mixture;
/// Molecular Rayleigh scattering from explicit refractivity and the shared gas state.
pub mod gas_scattering;
pub mod gas_species;
mod gas_state;
pub mod moist_air;
/// Shared dilute-solution activities, finite inventories and selective membrane transport.
pub mod osmosis;
pub mod graph;
/// Contact-temperature history and explicit strength handoff for polymer welding.
pub mod healing;
mod heat_capacity_state;
pub mod hyper;
pub mod identifiability;
pub mod iso9613;
/// Equilibrium charge carriers and sourced collisional DC transport.
pub mod ionized_gas;
pub mod phase;
/// Source-bound reciprocal longitudinal piezoelectric stress, field and storage.
pub mod piezoelectric;
/// Source-bound reciprocal incremental magnetostriction about a declared bias.
pub mod magnetostrictive;
/// Source-bound anisotropic reduced viscosity and work-conjugate plane transforms.
pub mod plane_stress_viscosity;
pub mod plastic;
pub mod porous;
/// Molecular relaxation absorption and dispersion of any gas mixture (M0.1-W2).
pub mod relaxation;
pub mod section;
/// Laminate (layered) plate/shell sections (IRP M1.2-W7).
pub mod laminate;
/// Source-bound viscoelastic surface films and explicit isotropic tensor lifts.
pub mod coating;
/// Temperature refresh of complete source-bound layered mechanical sections.
pub mod section_temperature;
/// Source-bound effective sheet mass, membrane and independent bending response.
pub mod sheet;
pub mod state_point;
/// Bounded solid-vapor pressure and latent-heat laws for shared phase transfer.
pub mod sublimation;
/// Original-water calorics and shared continuum evaporation/condensation.
pub mod evaporation;
pub mod tensor;
pub mod thermal;
/// Source-bound reciprocal Seebeck and Peltier cross-domain transport.
pub mod thermoelectric;
/// Shared species thermodynamics, acoustic medium and source-backed gas transport.
pub mod thermochemical_gas;
/// Bounded temperature-dependent Maxwell history, thermal strain and work balance.
pub mod thermovisco;
pub mod visco;
/// Finite-strain Maxwell branches on a hyperelastic equilibrium (multiplicative split, incremental potential).
pub mod visco_hyper;
/// Passive full-tensor Prony memory with work-conjugate stress and heat.
pub mod visco_tensor;
/// Wood equilibrium moisture content and the FPL moisture-content adjustments of stiffness and density.
pub mod wood_moisture;

pub use calibrate::{CalibrationFit, calibrate_bilinear};
pub use elastic::{IsotropicElastic, OrthotropicElastic};
pub use fiber::{ManderConcrete, MenegottoPintoSteel, Uniaxial, WoolFelt, WoolFeltState};
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

/// Borrowed constitutive laws retain their original history and source owner.
/// Reductions and element assembly can reuse an immutable law without cloning
/// its coefficient tables or source receipts at every integration point.
impl<L: SmallStrainLaw + ?Sized> SmallStrainLaw for &L {
    type State = L::State;

    fn initial_state(&self) -> Self::State {
        (**self).initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        (**self).stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        (**self).tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        (**self).update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        (**self).admissibility()
    }
    fn card(&self) -> ModelCard {
        (**self).card()
    }
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
