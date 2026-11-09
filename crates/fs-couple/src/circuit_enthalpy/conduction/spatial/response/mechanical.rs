//! One nonlinear current, rigid momentum and finite-enthalpy endpoint.
//!
//! The electrical response sees the actual mechanical midpoint. Its original
//! local Joule heat enters the same endpoint thermal residual, and the native
//! impulse owner accounts for kinetic storage separately from dissipation.

use super::*;
use crate::circuit_enthalpy::{
    conduction::{
        heat::ConductionHeatError,
        spatial::mechanical::{
            SpatialConductionKick,
            response::{
                SpatialConductionResponseKick, SpatialConductionResponseKickControls,
                SpatialConductionResponseKickError, kick_response_with_materials_candidate,
            },
        },
    },
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
};
use fs_em::spatial::conduction::response::motion::SpatialMovingConductionResponseDrive;
use fs_exec::Cx;
use fs_ga::Wrench;
use fs_mbd::articulated::{FreeFloatingBaseState, FreeFloatingImpulseResponse};
use fs_mbd_material::deformation::TetrahedralMass;
use fs_solver::NewtonReport;

/// Native current/force and mechanical endpoint from the same thermal trial.
/// The original source observations remain in the enclosing field report.
#[derive(Debug)]
pub struct SpatialConductionResponseMechanicalField {
    /// Original impulse result, current snapshot and independently closed work.
    pub kick: SpatialConductionKick,
    /// Mechanical midpoint solve; distinct from electrical source iterations
    /// and the enclosing finite-enthalpy Newton solve.
    pub nonlinear: NewtonReport,
}

impl ResponseHeatField for SpatialConductionResponseMechanicalField {
    fn joule_power_w(&self) -> f64 {
        self.kick.field.joule_power_w()
    }
    fn cell_joule_power_w(&self) -> &[f64] {
        self.kick.field.cell_joule_power_w()
    }
    fn contact_response(&self) -> &[SpatialConductionContactResponse] {
        self.kick.field.contact_response()
    }
    fn mechanical_power_w(&self) -> f64 {
        self.kick.field.mechanical_power_w()
    }
}

/// Accepted original source owner, finite thermal history and native mechanical
/// impulse endpoint. The report's total energy residual includes kinetic
/// storage and independently authored mechanical work.
pub type SpatialConductionResponseMechanicalEnthalpyCandidate<'t, S> =
    SpatialConductionResponseEnthalpyCandidate<'t, S, SpatialConductionResponseMechanicalField>;

impl<'t, S>
    SpatialConductionResponseEnthalpyCandidate<'t, S, SpatialConductionResponseMechanicalField>
{
    /// Native proposed mechanical endpoint; consume it together with the field
    /// source owner and thermal history from this same candidate.
    pub fn mechanical(&self) -> &FreeFloatingImpulseResponse {
        &self.report.electrical.field.field.kick.mechanical
    }
}

fn kick_error<E: From<ConductionThermalError>>(
    error: SpatialConductionResponseKickError,
    bounded: usize,
    requested: usize,
) -> CircuitEnthalpyError<E> {
    match error {
        SpatialConductionResponseKickError::Cancelled => CircuitEnthalpyError::Cancelled,
        SpatialConductionResponseKickError::FieldSolveBudget { .. } if bounded < requested => {
            CircuitEnthalpyError::CircuitBudget
        }
        other => CircuitEnthalpyError::Source(
            ConductionThermalError::SpatialResponseMechanical(other).into(),
        ),
    }
}

fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

// Move native arrays and original source observations into the common thermal
// field seam. No current, force, material or mechanical result is reconstructed.
fn thermal_field(
    report: SpatialConductionResponseKick,
) -> SpatialConductionResponseReport<SpatialConductionResponseMechanicalField> {
    let SpatialConductionResponseReport {
        field,
        terminals,
        evaluated,
        evaluated_contacts,
        iterations,
        contact_quadrature_refinement,
        ..
    } = report.field;
    SpatialConductionResponseReport {
        field: SpatialConductionResponseMechanicalField {
            kick: SpatialConductionKick {
                mechanical: report.mechanical,
                field,
                midpoint_twist_body: report.midpoint_twist_body,
                lorentz_wrench_body: report.lorentz_wrench_body,
                terminal_work_j: report.terminal_work_j,
                joule_heat_j: report.joule_heat_j,
                mechanical_work_j: report.mechanical_work_j,
                external_work_j: report.external_work_j,
                angular_velocity_residual_rad_s: report.angular_velocity_residual_rad_s,
                linear_velocity_residual_m_s: report.linear_velocity_residual_m_s,
                mechanical_work_residual_j: report.mechanical_work_residual_j,
                total_energy_residual_j: report.total_energy_residual_j,
                field_solves: report.field_solves,
            },
            nonlinear: report.nonlinear,
        },
        terminals,
        evaluated,
        evaluated_contacts,
        iterations,
        // Charge all mechanical residual/Jacobian trials and final admission,
        // not merely the final converged material-response field.
        field_solves: report.field_solves,
        source_evaluations: report.source_evaluations,
        contact_quadrature_refinement,
    }
}

/// Advance nonlinear moving current, a sourced rigid momentum kick and any
/// complete heat assembly on one common finite interval.
///
/// At every proposed T/phase endpoint, the original electrical sources are
/// refreshed, the original nonlinear field and actual mechanical midpoint are
/// jointly solved, and only their local conductor-frame Joule heat enters the
/// original finite cells. Every trial starts from the same supplied mechanical
/// state and thermal history. Lorentz work changes native kinetic storage;
/// independently authored external mechanical work stays separate. The final
/// global account includes both thermal and kinetic storage.
///
/// The original tetrahedral mass must cover the same contact-free electrical
/// mesh and material cards. Its original pinned density is re-evaluated at the
/// proposed material state and must remain unchanged and independent of owned
/// electrical coordinates. Existing caloric admission additionally retains its
/// constant-density, fixed-volume chart requirements. Thermal contacts,
/// radiation, reservoirs and finite links retain their original heat owners.
///
/// This is a fixed-configuration momentum substep: pose, geometry, inertia and
/// imposed B remain held. Position/rotation transport, induction and advective
/// inventory require their own owners and a new configuration admission.
/// Every native electrical solve, including discarded midpoint Newton probes,
/// consumes the shared thermal circuit-trial allowance. Failure, cancellation,
/// unsupported source states or failed work gates accepts none of the owners.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_conduction_response_mechanical_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConductionResponse,
    mass: &TetrahedralMass<'_>,
    drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
    flux_density_t: &[[f64; 3]],
    initial: FreeFloatingBaseState,
    external_wrench_body: Wrench,
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    electromechanical: SpatialConductionResponseKickControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialConductionResponseMechanicalEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let drive = drive.into();
    let dt = target_s - thermal.time_s();
    // The common thermal owner bounds each whole coupled midpoint solve.
    // Its inner material iterations retain their separate authored limits.
    let response_budget = SpatialConductionResponseControls {
        maximum_field_solves: electromechanical.physical.maximum_field_solves,
        ..electromechanical.response
    };
    let mut next = heat::step_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        problem,
        regions,
        volume_tolerance_m3,
        assembly,
        controls,
        &cancel,
        |layout, resolve| {
            layout.step_response_driver_candidate(
                thermal,
                target_s,
                problem,
                &[],
                assembly,
                response_budget,
                controls,
                initial_energy_increments_j,
                resolve,
                |sample: &A::Sample| &sample.as_ref().material,
                |sources, contacts, bounded| {
                    if !contacts.is_empty() {
                        return Err(map_conduction_error(invalid(
                            "rigid nonlinear current requires a contact-free electrical mesh",
                        )));
                    }
                    let bounded_kick = SpatialConductionResponseKickControls {
                        physical: super::super::mechanical::SpatialConductionKickControls {
                            maximum_field_solves: bounded.maximum_field_solves,
                            ..electromechanical.physical
                        },
                        ..electromechanical
                    };
                    let candidate = kick_response_with_materials_candidate(
                        problem,
                        sources,
                        mass,
                        drive,
                        flux_density_t,
                        dt,
                        initial,
                        external_wrench_body,
                        assembly,
                        bounded_kick,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| {
                        kick_error(
                            error,
                            bounded.maximum_field_solves,
                            electromechanical.physical.maximum_field_solves,
                        )
                    })?;
                    let (problem, report) = candidate.into_parts();
                    let supplied_power_w = match &report.field.terminals {
                        Some(terminals) => terminals.supplied_power_w,
                        None => report.field.field.terminal_power_w(),
                    };
                    Ok((problem, thermal_field(report), supplied_power_w))
                },
                &cancel,
            )
        },
    )?;
    // The shared thermal residual used terminal - Lorentz work as its supply.
    // Add the independent native delta_K - external_work - Lorentz_work
    // residual to obtain delta_H + delta_K - terminal - external - other_heat.
    let mechanical_residual = next
        .report
        .electrical
        .field
        .field
        .kick
        .mechanical_work_residual_j;
    next.report.total_energy_residual_j =
        finite(next.report.total_energy_residual_j + mechanical_residual)?;
    if next.report.total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
        return Err(map_conduction_error(invalid(
            "combined nonlinear current, kinetic and thermal energy do not close",
        )));
    }
    closed_balance(
        dt,
        next.report.total_energy_residual_j,
        next.report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}
