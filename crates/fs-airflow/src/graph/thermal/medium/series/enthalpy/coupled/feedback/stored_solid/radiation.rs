//! Implicit material-backed radiation in the stored solid/coolant transaction.
//!
//! Radiation is external heat, never coolant enthalpy. Original convection on
//! network regions alone feeds the liquid solver; additional Robin regions are
//! allowed only when explicitly named by a radiation patch. Their original
//! convection belongs to an external reservoir as well. Named zero-flux traces
//! use the same radiation-only boundary supported by the conduction owner.
//! Uniform and heterogeneous solid charts share the complete transaction.
//! The existing area-mean gray-patch model and emissivity snapshot are retained;
//! this adapter adds no radiation, source-material or time-integration law.

use std::collections::BTreeSet;

use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionError, RobinFlux,
    ThermalBoundary,
};
use fs_conduction::transient::enthalpy::{EnthalpyStepSolution,
    radiation::EnthalpyRadiationStepSolution};
use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;

use super::{ParallelLiquidNetwork, PumpHeadCurve, StoredCoolingControls,
    StoredCoolingError, StoredMaterialCoolingState, StoredCoolingStep, StoredSolidProblem,
    SolidEnthalpyStorage,
    StoredThermalResponse, poll};

/// All previous gates remain active. Actual inner enthalpy solves have their
/// own whole-step cap because one radiative response may need many of them.
#[derive(Clone, Copy, Debug)]
pub struct RadiatingStoredCoolingControls {
    pub thermal: StoredCoolingControls,
    pub radiation: AmbientRadiationConfig,
    pub max_total_enthalpy_solves: usize,
}

/// Original radiative FEM solution plus an explicit, non-overlapping heat split.
#[derive(Debug)]
pub struct RadiatingStoredSolidResponse {
    /// Complete original field, physical nodal/joule checks, applied boundary,
    /// nonlinear radiative heat and all material/patch records.
    pub thermal: EnthalpyRadiationStepSolution,
    /// Original convection on coolant regions, in network branch-major order.
    pub coolant_robin_fluxes: Vec<RobinFlux>,
    /// Sum of the above, W outward; excludes every radiative contribution.
    pub coolant_out_w: f64,
    /// Original convection on explicitly declared non-coolant patches, W outward.
    pub external_convection_out_w: f64,
}

impl StoredThermalResponse for RadiatingStoredSolidResponse {
    fn conduction(&self) -> &EnthalpyStepSolution { &self.thermal.conduction }
    fn coolant_fluxes(&self) -> &[RobinFlux] { &self.coolant_robin_fluxes }
    fn coolant_out_w(&self) -> f64 { self.coolant_out_w }
    fn external_boundary_out_w(&self) -> f64 {
        self.external_convection_out_w + self.thermal.radiation.nonlinear_radiation_out_w
    }
}

/// Accepted joint history and exact cumulative physical-solver work. Counts
/// span all inlet, pump-pressure, wall and radiation iterations in this step.
#[derive(Debug)]
pub struct RadiatingStoredCoolingStep<'p, 'a> {
    pub step: StoredCoolingStep<'p, 'a, RadiatingStoredSolidResponse>,
    pub total_enthalpy_solves: usize,
    pub total_newton_iterations: usize,
    pub total_krylov_iterations: usize,
}

impl<'n, 'a, 'card> ParallelLiquidNetwork<'n, 'a, 'card> {
    /// Advance both original histories with implicit ambient surface radiation.
    ///
    /// `boundary` maps current coolant rows and any explicitly radiating
    /// external Robin regions to the existing mesh. Coolant trace names and
    /// coefficients remain checked by the original network. Additional Robin
    /// regions MUST have a named radiation patch; they cannot silently become
    /// coolant or disappear from the external energy account. Radiation and
    /// convection may share a region only when that physical model is intended.
    ///
    /// All thermal calls use the same old nodal enthalpy and duration. The
    /// existing radiative owner additionally reassembles the physical endpoint
    /// nodal residual and nonlinear radiative energy; the joint transaction
    /// counts only separated coolant convection as reservoir return energy.
    /// Failed source, boundary, convergence, work or cancellation gates publish
    /// neither history. No enclosure, participating medium, evolving emissivity
    /// law, fluid residence time or hydraulic-heat deposition is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stored_solid_with_pump_and_radiation<'p, F, M: SolidEnthalpyStorage + ?Sized>(
        &self, cx: &Cx<'_>, state: &mut StoredMaterialCoolingState<'a, '_, M>,
        end_time_s: f64, reservoir_external_heat_j: f64, initial_walls_k: &[f64],
        pump: &'p PumpHeadCurve, problem: StoredSolidProblem<'_>,
        patches: &[AmbientRadiationPatch], controls: RadiatingStoredCoolingControls,
        boundary: F,
    ) -> Result<RadiatingStoredCoolingStep<'p, 'a>, StoredCoolingError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<ThermalBoundary, ConductionError>,
    {
        poll(cx)?;
        controls.radiation.validate()?;
        if controls.max_total_enthalpy_solves == 0 || patches.is_empty() {
            return Err(StoredCoolingError::Invalid("radiative storage needs patches and positive total enthalpy work"));
        }
        let mut patch_names = BTreeSet::new();
        for patch in patches {
            if !patch_names.insert(patch.region()) {
                return Err(StoredCoolingError::Invalid("duplicate stored-solid radiation patch"));
            }
        }
        let coolant_names: Vec<&str> = self.branches().iter()
            .flat_map(|b| b.path.sections.iter().map(|s| s.name.as_str())).collect();
        let storage = state.storage();
        let mut enthalpy_solves = 0usize;
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
                let radiation = AmbientRadiationConfig {
                    max_iterations: controls.radiation.max_iterations.min(remaining),
                    ..controls.radiation
                };
                let thermal = storage.advance_with_ambient_radiation(cx, mapped,
                    problem.interfaces, old, dt, controls.thermal.solid, patches, radiation)?;
                enthalpy_solves = add_work(enthalpy_solves, thermal.radiation.iterations)?;
                newton_iterations = add_work(newton_iterations, thermal.radiation.solid_iterations)?;
                krylov_iterations = add_work(krylov_iterations, thermal.radiation.krylov_iterations)?;
                partition(thermal, &coolant_names, &patch_names)
            })?;
        // No fallible operation after the shared owner publishes both histories.
        Ok(RadiatingStoredCoolingStep { step, total_enthalpy_solves: enthalpy_solves,
            total_newton_iterations: newton_iterations, total_krylov_iterations: krylov_iterations })
    }
}

pub(super) fn add_work(total: usize, increment: usize) -> Result<usize, StoredCoolingError> {
    total.checked_add(increment).ok_or(StoredCoolingError::Invalid("radiative thermal work count overflow"))
}

// Shared by the fixed-contact and source-dependent-contact response owners.
// A new constitutive response must not invent another destination for heat.
pub(super) fn admit_heat_ownership(
    boundary: &ThermalBoundary, coolant_names: &[&str], patch_names: &BTreeSet<&str>,
) -> Result<(), StoredCoolingError> {
    for (name, condition) in boundary.region_names().iter().zip(boundary.conditions()) {
        if matches!(condition, fs_conduction::ThermalBc::Robin { .. })
            && !coolant_names.contains(&name.as_str())
            && !patch_names.contains(name.as_str())
        {
            return Err(StoredCoolingError::Invalid("non-coolant Robin heat requires an explicit radiation patch"));
        }
    }
    Ok(())
}

fn partition(
    thermal: EnthalpyRadiationStepSolution, coolant_names: &[&str], patch_names: &BTreeSet<&str>,
) -> Result<RadiatingStoredSolidResponse, StoredCoolingError> {
    let (coolant_robin_fluxes, coolant_out_w, external_convection_out_w) = partition_convection(
        &thermal.convective_robin_fluxes, coolant_names, patch_names,
        thermal.radiation.nonlinear_radiation_out_w,
    )?;
    Ok(RadiatingStoredSolidResponse { thermal, coolant_robin_fluxes,
        coolant_out_w, external_convection_out_w })
}

pub(super) fn partition_convection(
    fluxes: &[RobinFlux], coolant_names: &[&str], patch_names: &BTreeSet<&str>,
    radiation_out_w: f64,
) -> Result<(Vec<RobinFlux>, f64, f64), StoredCoolingError> {
    let mut seen = BTreeSet::new();
    let mut external = 0.0;
    for flux in fluxes {
        if !seen.insert(flux.region.as_str()) || !flux.heat_rate_w.is_finite() {
            return Err(StoredCoolingError::Invalid("duplicate or nonfinite radiative convection record"));
        }
        if !coolant_names.contains(&flux.region.as_str()) {
            if !patch_names.contains(flux.region.as_str()) {
                return Err(StoredCoolingError::Invalid("unowned non-coolant convection record"));
            }
            external += flux.heat_rate_w;
        }
    }
    let coolant_robin_fluxes: Vec<_> = coolant_names.iter().map(|name| {
        fluxes.iter().find(|f| f.region == *name).cloned()
            .ok_or(StoredCoolingError::Invalid("radiative response lost a coolant trace"))
    }).collect::<Result<_, _>>()?;
    let coolant_out_w = coolant_robin_fluxes.iter().map(|f| f.heat_rate_w).sum::<f64>();
    if !external.is_finite() || !coolant_out_w.is_finite()
        || !(external + radiation_out_w).is_finite()
    {
        return Err(StoredCoolingError::Invalid("radiative heat partition is unrepresentable"));
    }
    Ok((coolant_robin_fluxes, coolant_out_w, external))
}
