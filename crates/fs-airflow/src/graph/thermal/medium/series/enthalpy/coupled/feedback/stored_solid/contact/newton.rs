//! Spatial contact/radiation feedback in the original solid/tank/pump transaction.
//!
//! Each solid callback performs one joint primary-enthalpy Newton solve. Source
//! contact resistance, local P1 emission and reflected enclosure irradiation are
//! evaluated at the same endpoint. The shared cooling transaction still owns
//! fluid properties, hydraulic flow, return enthalpy and atomic publication.

use std::collections::BTreeSet;

use fs_conduction::interface::temperature::{
    TemperatureContacts, enthalpy::ContactEnthalpyNewtonConfig,
};
use fs_conduction::radiation::{GrayDiffuseEnclosure, OpenDiffuseEnclosure};
use fs_conduction::{
    AmbientRadiationPatch, ConductionError, MaterialAmbientRadiationPatch,
    MaterialPointwiseAmbientRadiation, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    PointwiseEnclosureTolerance, ThermalBoundary,
};
use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;

use super::super::radiation::{add_work, admit_heat_ownership, partition_convection};
use super::super::{
    HeterogeneousStoredCoolingState, LiquidResistanceModel, ParallelLiquidNetwork, PumpHeadCurve,
    StoredCoolingControls, StoredCoolingError, StoredSolidProblem, poll,
};
use super::{ContactStoredCoolingStep, ContactStoredSolidResponse};

/// Original spatial radiation model composed with the contact Newton response.
/// Emissivities, view factors and source support stay on their existing owners.
#[derive(Clone, Copy, Debug)]
pub enum StoredContactRadiation<'a> {
    /// Conduction, source-backed contacts and coolant convection only.
    None,
    /// Local P1 emission to independently prescribed ambient reservoirs.
    Ambient(&'a [AmbientRadiationPatch]),
    /// Original material emissivity refreshed at the actual nodal enthalpy.
    /// Each patch retains the exact caloric object assigned to its storage nodes.
    MaterialAmbient(&'a [MaterialAmbientRadiationPatch<'a>]),
    /// Reciprocal internal surface exchange; its net power is zero to roundoff.
    ClosedEnclosure {
        enclosure: &'a GrayDiffuseEnclosure,
        tolerance: PointwiseEnclosureTolerance,
    },
    /// Internal reflection plus one explicit external black reservoir.
    OpenEnclosure {
        enclosure: &'a OpenDiffuseEnclosure,
        ambient_temperature_k: f64,
        tolerance: PointwiseEnclosureTolerance,
    },
}

impl StoredContactRadiation<'_> {
    fn patch_names(&self) -> Result<BTreeSet<&str>, StoredCoolingError> {
        match self {
            Self::None => Ok(BTreeSet::new()),
            Self::Ambient(patches) => {
                if patches.is_empty() {
                    return Err(StoredCoolingError::Invalid(
                        "use None when stored contact radiation is absent",
                    ));
                }
                let mut names = BTreeSet::new();
                for patch in *patches {
                    if !names.insert(patch.region()) {
                        return Err(StoredCoolingError::Invalid(
                            "duplicate stored-contact radiation patch",
                        ));
                    }
                }
                Ok(names)
            }
            Self::MaterialAmbient(patches) => {
                if patches.is_empty() {
                    return Err(StoredCoolingError::Invalid(
                        "use None when stored contact radiation is absent",
                    ));
                }
                let mut names = BTreeSet::new();
                for patch in *patches {
                    if !names.insert(patch.region) {
                        return Err(StoredCoolingError::Invalid(
                            "duplicate stored-contact radiation patch",
                        ));
                    }
                }
                Ok(names)
            }
            Self::ClosedEnclosure { enclosure, .. } => Ok(enclosure
                .surfaces()
                .iter()
                .map(|surface| surface.name())
                .collect()),
            Self::OpenEnclosure { enclosure, .. } => Ok(enclosure
                .surfaces()
                .iter()
                .map(|surface| surface.name())
                .collect()),
        }
    }
}

/// One thermal integration policy, plus the original contact Newton work gates.
/// `thermal.solid` is the sole step configuration; no second Newton or energy
/// policy is selected implicitly inside the cooling adapter.
#[derive(Clone, Copy, Debug)]
pub struct ContactNewtonStoredCoolingControls {
    pub thermal: StoredCoolingControls,
    /// Residual callbacks per solid solve, including rejected source/line trials.
    pub max_residual_evaluations: usize,
    /// Original matching/overlap contact-tangent pair allowance per solve.
    pub max_face_pairs: usize,
    /// Independent complete nodal energy residual norm limit, J.
    pub residual_tolerance_j: f64,
    /// Actual joint Newton solves across all inlet, pump and wall trials.
    pub max_total_enthalpy_solves: usize,
}

impl<'n, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'n, 'a, 'card, H> {
    /// Advance a spatial contact/radiation assembly and finite coolant together.
    ///
    /// The original contact Newton owner combines its full tangent with the
    /// storage, nonlinear conductivity and selected local radiation tangent.
    /// Every candidate starts from the same old nodal enthalpy and tank state;
    /// latent fractions are carried by h, never reconstructed from T. All
    /// source, residual, energy, hydraulic and regional coolant gates remain.
    ///
    /// `contacts` owns every interface, so `problem.interfaces` must be absent.
    /// A pointwise radiation owner is prepared against EACH actual mapped
    /// boundary, retaining its exact mesh/partition binding. Additional external
    /// Robin surfaces must be named by the radiation model. Only the original
    /// coolant convection heats the return stream. Enclosure exchange between
    /// solid patches is internal; its signed net surface power, including any
    /// explicit black bath, enters the external account exactly once.
    /// Material ambient patches also refresh their pinned emissivity at every
    /// trial h, retaining temperature and latent-fraction feedback in Newton.
    ///
    /// One joint enthalpy solve is counted per solid trial. Work exhaustion or
    /// any physics/cancellation refusal publishes neither history. Fixed
    /// geometry/pressure and the existing quasi-steady pipe approximation remain;
    /// this adapter adds no contact mechanics or fluid phase-change model.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stored_solid_with_pump_and_contact_newton<'p, F>(
        &self,
        cx: &Cx<'_>,
        state: &mut HeterogeneousStoredCoolingState<'a, '_, '_, '_>,
        end_time_s: f64,
        reservoir_external_heat_j: f64,
        initial_walls_k: &[f64],
        pump: &'p PumpHeadCurve,
        problem: StoredSolidProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        radiation: StoredContactRadiation<'_>,
        controls: ContactNewtonStoredCoolingControls,
        boundary: F,
    ) -> Result<ContactStoredCoolingStep<'p, 'a>, StoredCoolingError>
    where
        F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<ThermalBoundary, ConductionError>,
    {
        poll(cx)?;
        if problem.interfaces.is_some() {
            return Err(StoredCoolingError::Invalid(
                "temperature contacts own the complete interface set; include fixed interfaces in that owner",
            ));
        }
        if controls.max_total_enthalpy_solves == 0 {
            return Err(StoredCoolingError::Invalid(
                "positive total contact Newton solve allowance required",
            ));
        }
        let config = ContactEnthalpyNewtonConfig {
            step: controls.thermal.solid,
            max_residual_evaluations: controls.max_residual_evaluations,
            max_face_pairs: controls.max_face_pairs,
            residual_tolerance_j: controls.residual_tolerance_j,
        };
        let patch_names = radiation.patch_names()?;
        let coolant_names: Vec<&str> = self
            .branches()
            .iter()
            .flat_map(|branch| {
                branch
                    .path
                    .sections
                    .iter()
                    .map(|section| section.name.as_str())
            })
            .collect();
        let storage = state.storage();
        let mut enthalpy_solves = 0usize;
        let mut contact_iterations = 0usize;
        let mut newton_iterations = 0usize;
        let mut krylov_iterations = 0usize;
        let step = self.advance_stored_using(
            cx,
            state,
            end_time_s,
            reservoir_external_heat_j,
            initial_walls_k,
            pump,
            problem,
            controls.thermal,
            boundary,
            |cx, mapped, old, dt| {
                if enthalpy_solves >= controls.max_total_enthalpy_solves {
                    return Err(StoredCoolingError::ThermalWorkBudget {
                        performed: enthalpy_solves,
                        maximum: controls.max_total_enthalpy_solves,
                    });
                }
                admit_heat_ownership(mapped.boundary, &coolant_names, &patch_names)?;
                let thermal = match radiation {
                    StoredContactRadiation::None => {
                        contacts.advance_enthalpy_newton(cx, storage, mapped, old, dt, config)?
                    }
                    StoredContactRadiation::Ambient(patches) => {
                        let pointwise = PointwiseAmbientRadiation::new(
                            cx,
                            mapped.mesh,
                            mapped.boundary,
                            patches,
                        )?;
                        contacts.advance_enthalpy_newton_with_pointwise_ambient_radiation(
                            cx, storage, mapped, old, dt, config, &pointwise,
                        )?
                    }
                    StoredContactRadiation::MaterialAmbient(patches) => {
                        let pointwise = MaterialPointwiseAmbientRadiation::new(
                            cx,
                            mapped.mesh,
                            mapped.boundary,
                            patches,
                        )?;
                        contacts.advance_enthalpy_newton_with_material_ambient_radiation(
                            cx, storage, mapped, old, dt, config, &pointwise,
                        )?
                    }
                    StoredContactRadiation::ClosedEnclosure {
                        enclosure,
                        tolerance,
                    } => {
                        let pointwise = PointwiseEnclosureRadiation::new(
                            cx,
                            mapped.mesh,
                            mapped.boundary,
                            enclosure,
                            tolerance,
                        )?;
                        contacts.advance_enthalpy_newton_with_pointwise_enclosure_radiation(
                            cx, storage, mapped, old, dt, config, &pointwise,
                        )?
                    }
                    StoredContactRadiation::OpenEnclosure {
                        enclosure,
                        ambient_temperature_k,
                        tolerance,
                    } => {
                        let pointwise = PointwiseEnclosureRadiation::new_open(
                            cx,
                            mapped.mesh,
                            mapped.boundary,
                            enclosure,
                            ambient_temperature_k,
                            tolerance,
                        )?;
                        contacts.advance_enthalpy_newton_with_pointwise_enclosure_radiation(
                            cx, storage, mapped, old, dt, config, &pointwise,
                        )?
                    }
                };
                enthalpy_solves = add_work(enthalpy_solves, 1)?;
                contact_iterations = add_work(contact_iterations, thermal.iterations)?;
                newton_iterations = add_work(newton_iterations, thermal.nonlinear_iterations)?;
                krylov_iterations = add_work(krylov_iterations, thermal.krylov_iterations)?;
                let (coolant_robin_fluxes, coolant_out_w, external_convection_out_w) =
                    partition_convection(
                        &thermal.physical.convective_robin_fluxes,
                        &coolant_names,
                        &patch_names,
                        thermal.physical.radiation_out_w,
                    )?;
                Ok(ContactStoredSolidResponse {
                    thermal,
                    coolant_robin_fluxes,
                    coolant_out_w,
                    external_convection_out_w,
                })
            },
        )?;
        // The common owner has already committed both histories; no fallible
        // work or uncounted final re-solve follows this publication.
        Ok(ContactStoredCoolingStep {
            step,
            total_enthalpy_solves: enthalpy_solves,
            total_contact_iterations: contact_iterations,
            total_newton_iterations: newton_iterations,
            total_krylov_iterations: krylov_iterations,
        })
    }
}
