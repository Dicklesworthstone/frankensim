//! Isotropic vector dielectric switching through the original vector play solver.
//!
//! One branch owns one three-dimensional polarization, with the Euclidean
//! pinning norm and spherical saturation constraint. Three independent scalar
//! memories would introduce preferred axes and are not this model. The native
//! constrained return map, consistent tangent and bounded displacement inverse
//! are reused without duplicating their numerical algorithms.
//!
//! `D = epsilon_b E + sum(P)`,
//! `W = epsilon_b |E|²/2 + sum(|P|²/(2chi))`, and
//! `Q = sum(Ec |P_new-P_old|)` retain electric SI quantities. Endpoint work
//! `E_new.delta D = delta W + Q + numerical_dissipation` is distinguished from
//! an exact continuous rotating-path integral. Refine turns and reversals;
//! numerical dissipation is never physical heat.
//!
//! The E-controlled convex incremental coenergy is
//! `Phi = epsilon_b |E|²/2 + sum(E.P - |P|²/(2chi) - Ec|P-P_old|)`.
//! Its gradient is the returned D and its Hessian is dD/dE. The D-controlled
//! conjugate potential is `W + Q`. Neither is additional physical storage.
//! Every Newton trial must keep the SAME accepted `P_old`; a trial never
//! becomes the next trial's hysteresis history.
//!
//! Isotropy is an explicit model choice, not inferred rotating-field evidence
//! from a uniaxial source. Fits remain fixed within their admitted source
//! domain. This supplies no domain kinetics, texture anisotropy, breakdown,
//! switching strain, temperature-dependent coefficients or empirical validation.

use super::{
    DielectricHysteresisBranch, FieldStorageError, ScalarDielectricHysteresis,
    source::{DielectricHysteresisSources, ResolvedDielectricHysteresis},
};
use crate::field_storage::hysteresis::vector::{
    VectorHysteresis, VectorHysteresisControls, VectorHysteresisError, VectorHysteresisState,
    VectorHysteresisStep, inverse_tangent,
};
use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::QueryPoint;

type Vector = [f64; 3];
type Matrix = [[f64; 3]; 3];

/// Physical electric residual allowances and bounded local work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorDielectricHysteresisControls {
    pub relative_tolerance: f64,
    pub absolute_field_tolerance_v_m: f64,
    pub absolute_displacement_tolerance_c_m2: f64,
    pub work_tolerance_j_m3: f64,
    pub maximum_branch_iterations: usize,
    pub maximum_inverse_iterations: usize,
    pub maximum_backtracks: usize,
}

impl VectorDielectricHysteresisControls {
    fn play(self) -> VectorHysteresisControls {
        VectorHysteresisControls {
            relative_tolerance: self.relative_tolerance,
            absolute_field_tolerance_a_m: self.absolute_field_tolerance_v_m,
            absolute_flux_tolerance_t: self.absolute_displacement_tolerance_c_m2,
            work_tolerance_j_m3: self.work_tolerance_j_m3,
            maximum_branch_iterations: self.maximum_branch_iterations,
            maximum_inverse_iterations: self.maximum_inverse_iterations,
            maximum_backtracks: self.maximum_backtracks,
        }
    }
}

/// Source, original constitutive, electric residual or cancellation refusal.
#[derive(Debug)]
pub enum VectorDielectricHysteresisError {
    Material(FieldStorageError),
    Invalid(&'static str),
    Cancelled,
    DisplacementRange {
        observed_c_m2: f64,
        maximum_c_m2: f64,
    },
    BranchNonconvergence {
        branch: usize,
        iterations: usize,
    },
    DisplacementNonconvergence {
        iterations: usize,
        residual_c_m2: f64,
    },
}

// Preserve native diagnoses while expressing the few quantity-specific labels
// in the facade's units. This translation runs only on the cold refusal path.
fn electric_reason(reason: &'static str) -> &'static str {
    match reason {
        "hysteresis needs positive finite background permeability and field domain, with at most 1024 branches" => {
            "dielectric hysteresis needs positive finite background permittivity and field domain, with at most 1024 branches"
        }
        "magnetic history is not representable on its declared field domain" => {
            "dielectric history is not representable on its declared field domain"
        }
        "vector magnetic arithmetic is unrepresentable" => {
            "vector dielectric arithmetic is unrepresentable"
        }
        "vector magnetic field exceeds its original magnitude domain" => {
            "electric field exceeds its original magnitude domain"
        }
        "vector incremental permeability is not representably SPD" => {
            "vector incremental permittivity is not representably SPD"
        }
        "finite vector flux density required" => "finite vector electric displacement required",
        "retained vector inverse missed its physical flux residual" => {
            "retained vector inverse missed its physical electric-displacement residual"
        }
        other => other,
    }
}

impl From<FieldStorageError> for VectorDielectricHysteresisError {
    fn from(value: FieldStorageError) -> Self {
        match value {
            FieldStorageError::Invalid(reason) => Self::Invalid(electric_reason(reason)),
            other => Self::Material(other),
        }
    }
}

impl From<VectorHysteresisError> for VectorDielectricHysteresisError {
    fn from(value: VectorHysteresisError) -> Self {
        match value {
            VectorHysteresisError::Material(error) => error.into(),
            VectorHysteresisError::Invalid(reason) => Self::Invalid(electric_reason(reason)),
            VectorHysteresisError::Cancelled => Self::Cancelled,
            VectorHysteresisError::FluxRange {
                observed_t,
                maximum_t,
            } => Self::DisplacementRange {
                observed_c_m2: observed_t,
                maximum_c_m2: maximum_t,
            },
            VectorHysteresisError::BranchNonconvergence { branch, iterations } => {
                Self::BranchNonconvergence { branch, iterations }
            }
            VectorHysteresisError::FluxNonconvergence {
                iterations,
                residual_t,
            } => Self::DisplacementNonconvergence {
                iterations,
                residual_c_m2: residual_t,
            },
        }
    }
}

impl core::fmt::Display for VectorDielectricHysteresisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "vector dielectric history: {self:?}")
    }
}
impl std::error::Error for VectorDielectricHysteresisError {}
type Result<T> = core::result::Result<T, VectorDielectricHysteresisError>;

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(VectorDielectricHysteresisError::Invalid(
            "vector dielectric arithmetic is unrepresentable",
        ))
    }
}
fn dot(a: Vector, b: Vector) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

/// Immutable isotropic fit sharing the original scalar owner's parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorDielectricHysteresis {
    play: VectorHysteresis,
}

/// Private accepted three-dimensional electric field and polarization memory.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorDielectricHysteresisState {
    play: VectorHysteresisState,
}
impl VectorDielectricHysteresisState {
    pub const fn field_v_m(&self) -> Vector {
        self.play.field_a_m()
    }
    pub const fn displacement_c_m2(&self) -> Vector {
        self.play.flux_density_t()
    }
    pub fn polarizations_c_m2(&self) -> &[Vector] {
        self.play.polarizations_t()
    }
    pub const fn stored_energy_j_m3(&self) -> f64 {
        self.play.stored_energy_j_m3()
    }
}

/// Pure variational increment, with physical and numerical losses kept apart.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorDielectricHysteresisStep {
    pub state: VectorDielectricHysteresisState,
    pub irreversible_heat_j_m3: f64,
    pub energy_change_j_m3: f64,
    /// E_new dot the actually retained delta D, not an exact path integral.
    pub endpoint_work_j_m3: f64,
    pub numerical_dissipation_j_m3: f64,
    pub energy_residual_j_m3: f64,
    /// Frozen-history dD/dE [F/m]; selected active-set derivative at kinks.
    pub incremental_permittivity_f_m: Matrix,
    pub maximum_constitutive_residual_v_m: f64,
    pub inverse_iterations: usize,
    pub branch_iterations: usize,
    pub inverse_residual_c_m2: f64,
    coenergy_change_j_m3: f64,
}

impl VectorDielectricHysteresisStep {
    fn from_play(
        previous: &VectorDielectricHysteresisState,
        step: VectorHysteresisStep,
    ) -> Result<Self> {
        let de = std::array::from_fn(|i| step.state.field_a_m()[i] - previous.field_v_m()[i]);
        // Phi(E)-Phi(E_old) = delta E.D_old + endpoint_work - delta W - Q.
        // Use the native nonnegative remainder plus its admitted residual to
        // avoid subtracting complete stored energies about a polarized bias.
        let coenergy_change_j_m3 = finite(
            dot(de, previous.displacement_c_m2())
                + step.numerical_dissipation_j_m3
                + step.energy_residual_j_m3,
        )?;
        Ok(Self {
            state: VectorDielectricHysteresisState { play: step.state },
            irreversible_heat_j_m3: step.irreversible_heat_j_m3,
            energy_change_j_m3: step.energy_change_j_m3,
            endpoint_work_j_m3: step.endpoint_work_j_m3,
            numerical_dissipation_j_m3: step.numerical_dissipation_j_m3,
            energy_residual_j_m3: step.energy_residual_j_m3,
            incremental_permittivity_f_m: step.incremental_permeability_h_m,
            maximum_constitutive_residual_v_m: step.maximum_constitutive_residual_a_m,
            inverse_iterations: step.inverse_iterations,
            branch_iterations: step.branch_iterations,
            inverse_residual_c_m2: step.inverse_residual_t,
            coenergy_change_j_m3,
        })
    }

    /// D-controlled potential W_new+Q; its D-gradient is the endpoint E.
    pub fn incremental_potential_j_m3(&self) -> Result<f64> {
        finite(self.state.stored_energy_j_m3() + self.irreversible_heat_j_m3)
    }

    /// Same D-controlled potential relative to the accepted stored energy.
    pub fn incremental_potential_change_j_m3(&self) -> Result<f64> {
        finite(self.energy_change_j_m3 + self.irreversible_heat_j_m3)
    }

    /// E-controlled coenergy Phi=E.D-W-Q, with the SAME frozen old memory.
    /// This is a Newton objective, not reversible energy or physical heating.
    pub fn incremental_coenergy_j_m3(&self) -> Result<f64> {
        finite(
            dot(self.state.field_v_m(), self.state.displacement_c_m2())
                - self.state.stored_energy_j_m3()
                - self.irreversible_heat_j_m3,
        )
    }

    /// Stable Phi(E)-Phi(E_old), anchored to this increment's accepted history.
    /// Use the law's pairwise method for nearby trials far from that anchor.
    pub const fn incremental_coenergy_change_j_m3(&self) -> f64 {
        self.coenergy_change_j_m3
    }

    /// Consistent dE/dD [m/F] from the same native SPD inverse, no second fit.
    pub fn incremental_inverse_permittivity_m_f(&self) -> Result<Matrix> {
        inverse_tangent(self.incremental_permittivity_f_m).map_err(Into::into)
    }
}

impl VectorDielectricHysteresis {
    /// Explicit isotropic quadratic branches; zero branches is linear dielectric.
    pub fn new(
        background_permittivity_f_m: f64,
        field_limit_v_m: f64,
        branches: Vec<DielectricHysteresisBranch>,
    ) -> Result<Self> {
        Self::from_isotropic_scalar(&ScalarDielectricHysteresis::new(
            background_permittivity_f_m,
            field_limit_v_m,
            branches,
        )?)
    }

    /// Explicit isotropic interpretation of the SAME original scalar fit.
    /// Shares immutable parameters; it does not infer vector calibration.
    pub fn from_isotropic_scalar(law: &ScalarDielectricHysteresis) -> Result<Self> {
        Ok(Self {
            play: VectorHysteresis::from_isotropic_scalar(&law.play)?,
        })
    }

    pub fn background_permittivity_f_m(&self) -> f64 {
        self.play.background_permeability_h_m()
    }
    pub fn field_limit_v_m(&self) -> f64 {
        self.play.field_limit_a_m()
    }
    pub fn initial_state(&self) -> VectorDielectricHysteresisState {
        VectorDielectricHysteresisState {
            play: self.play.initial_state(),
        }
    }

    /// One backward variational increment from unchanged accepted memory.
    pub fn candidate_field(
        &self,
        previous: &VectorDielectricHysteresisState,
        field_v_m: Vector,
        controls: VectorDielectricHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorDielectricHysteresisStep> {
        let step =
            self.play
                .candidate_field(&previous.play, field_v_m, controls.play(), cancelled)?;
        VectorDielectricHysteresisStep::from_play(previous, step)
    }

    /// Bounded inverse of the SAME frozen-history D(E) response.
    /// No inverse trial commits; returned D is the constitutive result.
    pub fn candidate_displacement(
        &self,
        previous: &VectorDielectricHysteresisState,
        displacement_c_m2: Vector,
        controls: VectorDielectricHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorDielectricHysteresisStep> {
        let step = self.play.candidate_flux_density(
            &previous.play,
            displacement_c_m2,
            controls.play(),
            cancelled,
        )?;
        VectorDielectricHysteresisStep::from_play(previous, step)
    }

    /// Stable Phi(to)-Phi(from), with BOTH trials evaluated from `previous`.
    /// Factored quadratic and pinning-norm differences preserve small changes
    /// beneath unrelated field energy. Work is bounded by the admitted branch
    /// count, allocates nothing, and never replaces the frozen old memory.
    /// The caller owns candidate stationarity; this is not a new history step.
    pub fn incremental_coenergy_difference_j_m3(
        &self,
        previous: &VectorDielectricHysteresisState,
        from: &VectorDielectricHysteresisState,
        to: &VectorDielectricHysteresisState,
    ) -> Result<f64> {
        self.play
            .coenergy_difference_j_m3(&previous.play, &from.play, &to.play)
            .map_err(Into::into)
    }
}

/// Original scalar source fit with an explicitly selected isotropic vector law.
/// Uniaxial source evidence is retained, never promoted to rotational evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedVectorDielectricHysteresis {
    source: ResolvedDielectricHysteresis,
    law: VectorDielectricHysteresis,
}

impl ResolvedVectorDielectricHysteresis {
    pub fn from_isotropic_source(source: ResolvedDielectricHysteresis) -> Result<Self> {
        let law = VectorDielectricHysteresis::from_isotropic_scalar(source.law())?;
        Ok(Self { source, law })
    }
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source.material()
    }
    pub const fn sources(&self) -> &DielectricHysteresisSources {
        self.source.sources()
    }
    pub const fn scalar_source(&self) -> &ResolvedDielectricHysteresis {
        &self.source
    }
    pub const fn law(&self) -> &VectorDielectricHysteresis {
        &self.law
    }

    /// Recheck the same selected constant fit without altering vector memory.
    pub fn at_query_point(&self, point: &QueryPoint) -> Result<Self> {
        Ok(Self {
            source: self.source.at_query_point(point)?,
            law: self.law.clone(),
        })
    }
}
