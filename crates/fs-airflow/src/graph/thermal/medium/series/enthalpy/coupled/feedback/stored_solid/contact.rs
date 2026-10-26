//! Temperature-dependent material contacts in the existing stored cooling loop.
//!
//! The conduction owner resolves each interface's pinned R(T) source law and
//! checks the full physical enthalpy residual. The original cooling transaction
//! supplies pump/flow feedback and atomically publishes solid and tank history.
//! Contact exchange is internal; only original coolant convection enters return
//! enthalpy. Optional radiation remains external even on a radiation-only trace.

use std::collections::BTreeSet;

use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionError,
    RobinFlux, ThermalBoundary};
use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts,
    enthalpy::ContactEnthalpyStep};
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;

use super::{HeterogeneousStoredCoolingState, ParallelLiquidNetwork, PumpHeadCurve,
    StoredCoolingControls, StoredCoolingError, StoredCoolingStep, StoredSolidProblem,
    StoredThermalResponse, poll};
use super::radiation::{add_work, admit_heat_ownership, partition_convection};

/// The original coupled controls plus physical contact closure and one shared
/// allowance for actual enthalpy solves across ALL inlet/pressure/wall trials.
#[derive(Clone, Copy, Debug)]
pub struct ContactStoredCoolingControls {
    pub thermal: StoredCoolingControls,
    pub contact: ContactIterationConfig,
    pub max_total_enthalpy_solves: usize,
}

/// Original FEM/contact solution and its independently re-evaluated heat split.
/// Applied and endpoint contact receipts remain available on `thermal`.
#[derive(Debug)]
pub struct ContactStoredSolidResponse {
    pub thermal: ContactEnthalpyStep,
    /// Physical original coolant convection, in network branch-major order.
    pub coolant_robin_fluxes: Vec<RobinFlux>,
    /// Sum of coolant convection only, W outward from the solid.
    pub coolant_out_w: f64,
    /// Original convection on explicitly radiating non-coolant patches, W outward.
    pub external_convection_out_w: f64,
}

impl StoredThermalResponse for ContactStoredSolidResponse {
    fn conduction(&self) -> &EnthalpyStepSolution { &self.thermal.conduction }
    fn coolant_fluxes(&self) -> &[RobinFlux] { &self.coolant_robin_fluxes }
    fn coolant_out_w(&self) -> f64 { self.coolant_out_w }
    fn external_boundary_out_w(&self) -> f64 {
        self.external_convection_out_w + self.thermal.physical.radiation_out_w
    }
}

/// Accepted histories and cumulative physical-solver work, not just the final
/// contact response's counters. No additional solve follows publication.
#[derive(Debug)]
pub struct ContactStoredCoolingStep<'p, 'a> {
    pub step: StoredCoolingStep<'p, 'a, ContactStoredSolidResponse>,
    pub total_enthalpy_solves: usize,
    pub total_contact_iterations: usize,
    pub total_newton_iterations: usize,
    pub total_krylov_iterations: usize,
}

impl<'n, 'a, 'card> ParallelLiquidNetwork<'n, 'a, 'card> {
    /// Advance source-dependent solid interfaces, coolant, pump and reservoir
    /// on ONE immutable-history interval, optionally with ambient radiation.
    ///
    /// `contacts` is the complete contact owner, including any fixed interfaces;
    /// `problem.interfaces` must be absent instead of being silently ignored or
    /// counted twice. Use the original heterogeneous enthalpy owner also for a
    /// one-material contact assembly. Distinct material traces keep their own
    /// charts, masses and phase states, including latent plateaus.
    ///
    /// Every thermal response delegates to TemperatureContacts::advance_enthalpy
    /// or its radiation variant. The original R(T) source selection, full nodal
    /// joule residual and energy gates stay active. Returned coolant rows come
    /// from that response's fresh physical endpoint, never its combined
    /// convection/radiation rows or an internal contact heat-flow sum.
    ///
    /// Supply `None` when radiation is absent, otherwise a nonempty patch set
    /// and its bounded original controls. Additional external Robin surfaces
    /// must be named by that patch set, as in the existing radiation adapter.
    ///
    /// The shared enthalpy-work remainder caps the product of contact and
    /// radiation iterations BEFORE a response starts. Near exhaustion, this
    /// conservative cap may refuse before the caller's per-response maxima;
    /// it never enlarges tolerances or resets work per pressure/inlet trial.
    /// Any refusal or cancellation leaves both accepted histories unchanged.
    /// This is fixed-geometry thermal contact, not contact mechanics, evolving
    /// normal pressure, fluid boiling, or deposited hydraulic/motor heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stored_solid_with_pump_and_contacts<'p, F>(
        &self, cx: &Cx<'_>, state: &mut HeterogeneousStoredCoolingState<'a, '_, '_, '_>,
        end_time_s: f64, reservoir_external_heat_j: f64, initial_walls_k: &[f64],
        pump: &'p PumpHeadCurve, problem: StoredSolidProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        radiation: Option<(&[AmbientRadiationPatch], AmbientRadiationConfig)>,
        controls: ContactStoredCoolingControls, boundary: F,
    ) -> Result<ContactStoredCoolingStep<'p, 'a>, StoredCoolingError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<ThermalBoundary, ConductionError>,
    {
        poll(cx)?;
        if problem.interfaces.is_some() {
            return Err(StoredCoolingError::Invalid(
                "temperature contacts own the complete interface set; include fixed interfaces in that owner",
            ));
        }
        if controls.max_total_enthalpy_solves == 0 || controls.contact.max_iterations == 0
            || !(controls.contact.relaxation.is_finite()
                && controls.contact.relaxation > 0.0 && controls.contact.relaxation <= 1.0)
            || !(controls.contact.residual_tolerance_j.is_finite()
                && controls.contact.residual_tolerance_j > 0.0)
        {
            return Err(StoredCoolingError::Invalid("positive contact work/residual limits and relaxation in (0,1] required"));
        }
        let mut patch_names = BTreeSet::new();
        if let Some((patches, config)) = radiation {
            config.validate()?;
            if patches.is_empty() {
                return Err(StoredCoolingError::Invalid("use None when contact radiation is absent"));
            }
            for patch in patches {
                if !patch_names.insert(patch.region()) {
                    return Err(StoredCoolingError::Invalid("duplicate stored-contact radiation patch"));
                }
            }
        }
        let coolant_names: Vec<&str> = self.branches().iter()
            .flat_map(|b| b.path.sections.iter().map(|s| s.name.as_str())).collect();
        let storage = state.storage();
        let mut enthalpy_solves = 0usize;
        let mut contact_iterations = 0usize;
        let mut newton_iterations = 0usize;
        let mut krylov_iterations = 0usize;
        let step = self.advance_stored_using(cx, state, end_time_s, reservoir_external_heat_j,
            initial_walls_k, pump, problem, controls.thermal, boundary, |cx, mapped, old, dt| {
                let remaining = controls.max_total_enthalpy_solves.saturating_sub(enthalpy_solves);
                if remaining == 0 {
                    return Err(StoredCoolingError::ThermalWorkBudget {
                        performed: enthalpy_solves, maximum: controls.max_total_enthalpy_solves,
                    });
                }
                admit_heat_ownership(mapped.boundary, &coolant_names, &patch_names)?;
                let radiation_limit = radiation.map_or(1, |(_, config)| config.max_iterations).min(remaining);
                let contact = ContactIterationConfig {
                    max_iterations: controls.contact.max_iterations.min(remaining / radiation_limit),
                    ..controls.contact
                };
                let thermal = match radiation {
                    Some((patches, config)) => contacts.advance_enthalpy_with_ambient_radiation(
                        cx, storage, mapped, old, dt, controls.thermal.solid, contact, patches,
                        AmbientRadiationConfig { max_iterations: radiation_limit, ..config },
                    )?,
                    None => contacts.advance_enthalpy(
                        cx, storage, mapped, old, dt, controls.thermal.solid, contact,
                    )?,
                };
                let responses = if radiation.is_some() { thermal.radiation_iterations } else { thermal.iterations };
                enthalpy_solves = add_work(enthalpy_solves, responses)?;
                if enthalpy_solves > controls.max_total_enthalpy_solves {
                    return Err(StoredCoolingError::ThermalWorkBudget {
                        performed: enthalpy_solves, maximum: controls.max_total_enthalpy_solves,
                    });
                }
                contact_iterations = add_work(contact_iterations, thermal.iterations)?;
                newton_iterations = add_work(newton_iterations, thermal.nonlinear_iterations)?;
                krylov_iterations = add_work(krylov_iterations, thermal.krylov_iterations)?;
                let (coolant_robin_fluxes, coolant_out_w, external_convection_out_w) = partition_convection(
                    &thermal.physical.convective_robin_fluxes, &coolant_names, &patch_names,
                    thermal.physical.radiation_out_w,
                )?;
                Ok(ContactStoredSolidResponse { thermal, coolant_robin_fluxes,
                    coolant_out_w, external_convection_out_w })
            })?;
        Ok(ContactStoredCoolingStep { step, total_enthalpy_solves: enthalpy_solves,
            total_contact_iterations: contact_iterations, total_newton_iterations: newton_iterations,
            total_krylov_iterations: krylov_iterations })
    }
}
