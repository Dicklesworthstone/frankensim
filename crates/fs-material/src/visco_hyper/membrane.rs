//! Finite-strain membrane relaxation with retained three-dimensional memory.
//!
//! Engineering Green strain `[E11,E22,2E12]` supplies the same objective
//! Cholesky representative and exact kinematic chain rule as an elastic sheet.
//! The original finite-viscoelastic potential relaxes the transverse traction,
//! with either finite compliance or an explicit elastic-volume constraint.
//! Viscous metrics stay in the fixed natural material frame: removing spatial
//! rotation from the representative does not rotate or reset material memory.
//!
//! One frozen adapter borrows one accepted history and one actual time interval.
//! Trial `response` queries cannot advance either. `step_response` proposes the
//! native update only after the transverse solve succeeds; its physical storage,
//! loading work and relaxation loss remain distinct from the incremental sheet
//! potential. Spatial loading quadrature error is not part of this local loss.
//! No source qualification, bending, wrinkling, damage or thermal evolution is
//! inferred from the supplied explicit law and reference thickness.

use super::{
    FrozenViscoHyperelastic, FrozenViscoRemainder, ViscoHyperelastic, ViscoRemainderState,
    ViscoState,
};
use crate::hyper::incompressible_plane_stress::incompressible_plane_stress_with_seed;
use crate::hyper::mixed::VolumetricTreatment;
use crate::hyper::plane_stress::{
    PlaneStressResponse, PlaneStressSettings, plane_stress_with_seed,
};
use crate::sheet::hyperelastic::{admitted_representative, response_from_plane_stress};
use crate::sheet::{SheetError, SheetLaw, SheetResponse};

/// Original material history and its explicitly selected volume treatment.
/// The exact-volume variant retains a native remainder ledger; it never
/// creates a full-energy state or computes an unused bulk penalty.
#[derive(Clone, Debug, PartialEq)]
pub enum ViscoMembraneState {
    Full(ViscoState),
    Incompressible(ViscoRemainderState),
}

impl ViscoMembraneState {
    /// Native stress-free history at the law's natural reference deformation.
    /// A membrane with a different initial metric must first propose a zero-time
    /// elastic jump at that metric. The frozen adapter admits the selected law
    /// and history before any spatial trial; unsupported remainder splits refuse.
    pub fn initial(law: &ViscoHyperelastic, treatment: VolumetricTreatment) -> Self {
        match treatment {
            VolumetricTreatment::OriginalMaterial => Self::Full(law.initial_state()),
            VolumetricTreatment::ExactIncompressible => {
                Self::Incompressible(law.initial_remainder_state())
            }
        }
    }

    pub const fn volumetric_treatment(&self) -> VolumetricTreatment {
        match self {
            Self::Full(_) => VolumetricTreatment::OriginalMaterial,
            Self::Incompressible(_) => VolumetricTreatment::ExactIncompressible,
        }
    }

    /// Actual canonical deformation at the accepted material point.
    pub const fn deformation_gradient(&self) -> &[f64; 9] {
        match self {
            Self::Full(state) => &state.deformation_gradient,
            Self::Incompressible(state) => &state.deformation_gradient,
        }
    }

    /// Physical storage per original reference volume, not incremental potential.
    pub const fn stored_pa(&self) -> f64 {
        match self {
            Self::Full(state) => state.stored,
            Self::Incompressible(state) => state.stored,
        }
    }

    /// Cumulative native frozen-metric loading work per original reference volume.
    pub const fn work_pa(&self) -> f64 {
        match self {
            Self::Full(state) => state.work,
            Self::Incompressible(state) => state.work,
        }
    }

    /// Cumulative signed natural-reference work, separate from mechanical work.
    pub const fn reference_work_pa(&self) -> f64 {
        match self {
            Self::Full(state) => state.reference_work,
            Self::Incompressible(state) => state.reference_work,
        }
    }

    /// Cumulative physical branch relaxation loss per original reference volume.
    pub const fn dissipated_pa(&self) -> f64 {
        match self {
            Self::Full(state) => state.dissipated,
            Self::Incompressible(state) => state.dissipated,
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum FrozenPotential<'a> {
    Full(FrozenViscoHyperelastic<'a>),
    Incompressible(FrozenViscoRemainder<'a>),
}

/// Immutable original-history membrane potential for one spatial solve.
#[derive(Clone, Copy, Debug)]
pub struct FrozenViscoMembrane<'a> {
    potential: FrozenPotential<'a>,
    reference_thickness_m: f64,
    settings: PlaneStressSettings,
}

/// Complete proposed membrane endpoint. Dropping it leaves the accepted history
/// unchanged. The spatial owner publishes all points and its configuration only
/// after its own equilibrium, work and cancellation checks pass.
#[derive(Clone, Debug)]
pub struct ViscoMembraneStepResponse {
    pub state: ViscoMembraneState,
    /// Incremental potential, second-Piola resultants and exact Green-strain
    /// tangent at this endpoint. `sheet.energy_j_m2` is NOT stored energy.
    pub sheet: SheetResponse,
    /// Physical full Piola stress in the canonical spatial frame, in Pa.
    /// The exact-volume variant includes its pressure multiplier. A spatial
    /// consumer rotates this stress if it needs Cartesian three-dimensional
    /// components; the objective sheet resultants need no such rotation.
    pub first_piola_pa: [f64; 9],
    /// Compression-positive pressure multiplier for exact incompressibility.
    /// Finite-compliance stress needs no independent pressure and returns `None`.
    pub pressure_pa: Option<f64>,
    pub stored_energy_j_m2: f64,
    /// Original frozen-metric loading work for this interval, not cumulative.
    pub loading_work_j_m2: f64,
    /// Nonnegative native relaxation loss for this interval, not cumulative.
    pub dissipated_j_m2: f64,
    pub relative_transverse_residual: f64,
}

struct ReducedResponse {
    plane_stress: PlaneStressResponse,
    constrained_piola: Option<[f64; 9]>,
    pressure_pa: Option<f64>,
}

impl<'a> FrozenViscoMembrane<'a> {
    /// Freeze the original material history and one nonnegative physical time
    /// interval. Reference thickness is positive, finite, and belongs to the
    /// ORIGINAL mesh reference; it never follows current area or relaxed normal
    /// stretch. The native plane-stress owner checks its solve settings/domain
    /// on each query. Local iteration/branch caps bound work between spatial
    /// cancellation checkpoints.
    pub fn new(
        law: &'a ViscoHyperelastic,
        accepted: &'a ViscoMembraneState,
        dt_s: f64,
        reference_thickness_m: f64,
        settings: PlaneStressSettings,
    ) -> Result<Self, SheetError> {
        if !reference_thickness_m.is_finite() || reference_thickness_m <= 0.0 {
            return Err(SheetError::Invalid(
                "finite positive original reference thickness required",
            ));
        }
        let potential = match accepted {
            ViscoMembraneState::Full(state) => FrozenPotential::Full(law.freeze(state, dt_s)?),
            ViscoMembraneState::Incompressible(state) => {
                FrozenPotential::Incompressible(law.freeze_remainder(state, dt_s)?)
            }
        };
        Ok(Self {
            potential,
            reference_thickness_m,
            settings,
        })
    }

    fn reduced_response(&self, strain: [f64; 3]) -> Result<ReducedResponse, SheetError> {
        let in_plane = admitted_representative(strain)?;
        match self.potential {
            FrozenPotential::Full(frozen) => {
                let accepted = frozen.accepted.deformation_gradient;
                let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
                let seed = [
                    accepted[2],
                    accepted[5],
                    frozen.law.reference_deformation().volume_ratio() / determinant,
                ];
                Ok(ReducedResponse {
                    plane_stress: plane_stress_with_seed(&frozen, &in_plane, seed, self.settings)?,
                    constrained_piola: None,
                    pressure_pa: None,
                })
            }
            FrozenPotential::Incompressible(frozen) => {
                let accepted = frozen.accepted.deformation_gradient;
                let response = incompressible_plane_stress_with_seed(
                    &frozen,
                    &in_plane,
                    [accepted[2], accepted[5]],
                    self.settings,
                )?;
                Ok(ReducedResponse {
                    plane_stress: PlaneStressResponse {
                        deformation_gradient: response.deformation_gradient,
                        energy_density_pa: response.energy_density_pa,
                        piola: response.piola,
                        tangent: response.tangent,
                        relative_transverse_residual: response.relative_transverse_residual,
                        iterations: response.iterations,
                    },
                    constrained_piola: Some(response.full_piola),
                    pressure_pa: Some(response.pressure_pa),
                })
            }
        }
    }

    /// Propose one original material update at the converged transverse state.
    /// The same accepted metrics, actual interval, native reference mapping and
    /// exact chain rule as `response` are retained. This method does not mutate
    /// the borrowed accepted state or promote local work to spatial nodal work.
    pub fn step_response(&self, strain: [f64; 3]) -> Result<ViscoMembraneStepResponse, SheetError> {
        let reduced = self.reduced_response(strain)?;
        let sheet =
            response_from_plane_stress(strain, reduced.plane_stress, self.reference_thickness_m)?;
        let f = &reduced.plane_stress.deformation_gradient;
        let (state, first_piola_pa, loading_work_pa, dissipated_pa) = match self.potential {
            FrozenPotential::Full(frozen) => {
                let material = frozen.law.step_response(frozen.accepted, f, frozen.dt_s)?;
                (
                    ViscoMembraneState::Full(material.state),
                    material.first_piola_pa,
                    material.loading_work_pa,
                    material.dissipated_pa,
                )
            }
            FrozenPotential::Incompressible(frozen) => {
                let material =
                    frozen
                        .law
                        .step_remainder_response(frozen.accepted, f, frozen.dt_s)?;
                (
                    ViscoMembraneState::Incompressible(material.state),
                    reduced
                        .constrained_piola
                        .expect("exact-volume response carries pressure stress"),
                    material.loading_work_pa,
                    material.dissipated_pa,
                )
            }
        };
        let stored_energy_j_m2 = self.reference_thickness_m * state.stored_pa();
        let loading_work_j_m2 = self.reference_thickness_m * loading_work_pa;
        let dissipated_j_m2 = self.reference_thickness_m * dissipated_pa;
        if ![stored_energy_j_m2, loading_work_j_m2, dissipated_j_m2]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(SheetError::Invalid(
                "unrepresentable physical membrane work, storage or relaxation loss",
            ));
        }
        Ok(ViscoMembraneStepResponse {
            state,
            sheet,
            first_piola_pa,
            pressure_pa: reduced.pressure_pa,
            stored_energy_j_m2,
            loading_work_j_m2,
            dissipated_j_m2,
            relative_transverse_residual: reduced.plane_stress.relative_transverse_residual,
        })
    }
}

impl SheetLaw for FrozenViscoMembrane<'_> {
    fn response(&self, strain: [f64; 3]) -> Result<SheetResponse, SheetError> {
        let reduced = self.reduced_response(strain)?;
        response_from_plane_stress(strain, reduced.plane_stress, self.reference_thickness_m)
    }
}
