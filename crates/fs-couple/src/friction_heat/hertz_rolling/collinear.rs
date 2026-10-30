//! Event-resolved native Hertz/rolling trajectories with original finite heat.
//!
//! Every thermal trial resolves the original surface source at its own face
//! temperature and the rolling card at the sphere's bulk temperature. The
//! original collinear native owner integrates sliding, stops and subsequent
//! rolling/sticking phases with those held trial coefficients. Integrated
//! surface and bulk work enter their distinct original caloric destinations.

use super::RigidHertzRollingHeatContact;
use crate::friction_heat::hertz::{RigidHertzHeatContact, RigidHertzHeatSource, validate_geometry};
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use crate::friction_heat::rolling::RigidRollingHeatBinding;
use crate::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances, finite};
use fs_contact::interface_binding::thermal_friction::hertz::{
    SourceHertzWrench, SourceHertzWrenchError,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::hertz_rolling::RigidHertzRollingContact;
use fs_mbd::tangential_friction::hertz_rolling::collinear::{
    PreparedRigidHertzRollingCollinear, RigidHertzRollingCollinearCandidate,
    RigidHertzRollingCollinearControls, RigidHertzRollingCollinearError,
    RigidHertzRollingCollinearHistory,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// Original event/mechanical gates and the existing finite thermal root controls.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingCollinearHeatControls {
    pub mechanical: RigidHertzRollingCollinearControls,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical root scaling [J/K], never another physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidHertzRollingCollinearHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidHertzRollingCollinearError),
    FrictionSource(SourceHertzWrenchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidHertzRollingCollinearError> for RigidHertzRollingCollinearHeatError<E> {
    fn from(error: RigidHertzRollingCollinearError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceHertzWrenchError> for RigidHertzRollingCollinearHeatError<E> {
    fn from(error: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidHertzRollingCollinearHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidHertzRollingCollinearHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidHertzRollingCollinearHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "collinear Hertz/rolling finite heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidHertzRollingCollinearHeatError<E> {}

/// One complete trajectory, including its opaque stop continuation and original
/// rolling history, together with the original finite thermal endpoint.
pub struct HeatedRigidHertzRollingCollinearCandidate<'m, P> {
    pub mechanical: RigidHertzRollingCollinearCandidate,
    pub thermal: DepositionSimulation<'m>,
    /// Surface-network audit. Its injected source includes internal bulk work;
    /// use external_heat_j below for the physical system's outside supply.
    /// This adapter does not classify the bulk source's entropy.
    pub heat: FrictionHeatNetworkReport,
    pub friction_source: SourceHertzWrench,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    pub external_heat_j: f64,
    pub sphere_temperature_before_k: f64,
    pub sphere_temperature_after_k: f64,
    /// Actual caloric change minus outside heat and integrated physical losses.
    pub heat_transfer_residual_j: f64,
    /// Actual native kinetic plus original caloric change minus outside heat.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: RigidHertzRollingCollinearCandidate,
    friction_source: SourceHertzWrench,
    external_heat_j: f64,
    sphere_temperature_k: f64,
    source: P,
}

/// Propose event-resolved original shear/rolling work and finite body heating.
///
/// The native owner admits its invariant longitudinal-shear/lateral-rolling
/// subspace and retains the full finite-body recoil. Genuine stop continuation
/// comes from its opaque history, never a thermal wrapper's velocity threshold.
/// Every root trial restarts the complete trajectory from the same accepted
/// bodies and history, using one prepared native geometry and actual duration.
/// No intermediate stop, source receipt, body or caloric endpoint publishes.
///
/// The original Hertz source is queried at the solved common face temperature
/// and explicit nominal pressure. The rolling card is admitted at the sphere's
/// actual bulk temperature. Their coefficients are held throughout that trial's
/// mechanical window; its exact phase integrals supply surface and bulk heat.
/// Surface work enters the original massless face. Bulk rolling work enters
/// exactly one original whole-sphere cell with matching native mass, COM and
/// declared body/material identities. The callback supplies other heat only.
///
/// Thermal storage/transport uses the original backward-Euler endpoint model
/// over the whole window, including the portion after a mechanical stop. Thus
/// exact stopping for held trial coefficients does not establish the exact
/// stop for a continuously evolving temperature-dependent material. Refine
/// thermal windows independently. Source/topology events must already be cut
/// by the caller. Poses, normal load and contact radius remain fixed; no moving
/// thermal geometry, normal work, partial-slip memory, variable viscosity,
/// entropy production or root uniqueness is supplied here.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_rolling_collinear_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    history: &RigidHertzRollingCollinearHistory,
    contact: RigidHertzRollingHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    bulk_binding: RigidRollingHeatBinding,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzRollingCollinearHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRollingCollinearCandidate<'m, P>, RigidHertzRollingCollinearHeatError<E>>
{
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    controls.mechanical.validate()?;
    let target = accepted_time_s + dt_s;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(RigidHertzRollingCollinearHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || cells == 0
        || cells > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || !temperature_guess_k.is_finite()
        || temperature_guess_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || contact.interface != contact.source.interface().interface()
    {
        return Err(RigidHertzRollingCollinearHeatError::Invalid(
            "synchronized bounded thermal coordinates and the same original interface required",
        ));
    }
    let duration_s = target - accepted_time_s;
    let radius = det::sqrt(contact.patch.contact_area_m2() / core::f64::consts::PI);
    let initial_source = contact.source.hertz_at(
        temperature_guess_k,
        contact.nominal_pressure_pa,
        contact.patch.normal_load_n(),
        radius,
    )?;
    let temperature_before = thermal.temperature_k(bulk_binding.cell_id)?;
    let prepared = PreparedRigidHertzRollingCollinear::new(
        bodies,
        properties,
        history,
        RigidHertzRollingContact {
            binding: contact.binding,
            frame: contact.frame,
            hertz: initial_source.law(),
            patch: contact.patch,
            interface: contact.interface,
            rolling: contact.rolling,
            temperature_k: temperature_before,
            excitation_frequency_hz: contact.excitation_frequency_hz,
            work: contact.work,
        },
        duration_s,
        controls.mechanical,
        gate,
    )?;
    let sphere = prepared.deforming_body();
    let side = usize::from(sphere == contact.binding.body_b);
    let expected_cell = [route.patch.link.a, route.patch.link.b][side];
    let cell = thermal.cells().get(&bulk_binding.cell_id).ok_or(
        RigidHertzRollingCollinearHeatError::Invalid("inactive sphere bulk cell"),
    )?;
    if bulk_binding.body_index != sphere
        || bulk_binding.cell_id != expected_cell
        || cell.body_id != bulk_binding.thermal_body_id
        || cell.material_id != bulk_binding.material_id
        || cell.mass_kg != properties[sphere].mass()
    {
        return Err(RigidHertzRollingCollinearHeatError::Invalid(
            "whole-sphere mass and declared body/material must occupy the original route side",
        ));
    }
    for (&id, other) in thermal.cells() {
        checkpoint()?;
        if id != bulk_binding.cell_id && other.body_id == bulk_binding.thermal_body_id {
            return Err(RigidHertzRollingCollinearHeatError::Invalid(
                "whole-sphere thermal body must have exactly one active caloric cell",
            ));
        }
    }
    let point = prepared.contact_points_world()[0];
    let point = [point.x, point.y, point.z];
    let tolerance = controls
        .mechanical
        .mechanical
        .mechanical
        .geometry_tolerance_m;
    validate_geometry(
        point,
        RigidHertzHeatContact {
            binding: contact.binding,
            frame: contact.frame,
            source: contact.source,
            normal_force_n: contact.patch.normal_load_n(),
            contact_radius_m: radius,
            nominal_pressure_pa: contact.nominal_pressure_pa,
        },
        route,
        tolerance,
    )
    .map_err(RigidHertzRollingCollinearHeatError::Invalid)?;
    let sign = if side == 0 { 1.0 } else { -1.0 };
    let normal = contact.frame.normal();
    let delta = std::array::from_fn::<_, 3, _>(|i| {
        0.5 * cell.geometry.start_m[i] + 0.5 * cell.geometry.end_m[i]
            - (point[i] + sign * prepared.sphere_radius_m() * normal[i])
    });
    if delta.iter().any(|v| !v.is_finite())
        || det::hypot(det::hypot(delta[0], delta[1]), delta[2]) > tolerance
    {
        return Err(RigidHertzRollingCollinearHeatError::Invalid(
            "whole-sphere caloric centre must coincide with its original native COM",
        ));
    }
    let step = FrictionHeatNetworkStep::new(thermal, target, &[route], 1, checkpoint)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidHertzRollingCollinearHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidHertzRollingCollinearHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
        step.evaluate(
            &x[..cells],
            &[temperature],
            |endpoint, temperatures| {
                let sphere_temperature = endpoint.temperature_k(bulk_binding.cell_id)?;
                let friction_source = contact.source.hertz_at(
                    temperatures[0],
                    contact.nominal_pressure_pa,
                    contact.patch.normal_load_n(),
                    radius,
                )?;
                let mechanical = prepared.advance_at_temperature(
                    friction_source.law(),
                    sphere_temperature,
                    gate,
                )?;
                let RigidHertzHeatSource {
                    conductivity_w_m_k,
                    transport,
                } = resolve(endpoint, temperatures[0])
                    .map_err(RigidHertzRollingCollinearHeatError::Source)?;
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                checkpoint()?;
                let external_heat_j = finite(
                    duration_s
                        * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                )?;
                // The trajectory supplies integral c*omega(t)^2 dt. Squaring a
                // mean angular speed here would lose variance and stopping phases.
                let bulk_j = mechanical.report.bulk_heat_j;
                let bulk_w = finite(bulk_j / duration_s)?;
                if bulk_j > 0.0 && bulk_w == 0.0 {
                    return Err(DepositionError::Arithmetic.into());
                }
                if bulk_w != 0.0 {
                    let power = sources_w.entry(bulk_binding.cell_id).or_default();
                    let total = finite(*power + bulk_w)?;
                    if total == *power {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    *power = total;
                }
                Ok(FrictionHeatNetworkInput {
                    interfaces: vec![InterfaceHeatInput {
                        dissipated_j: mechanical.report.surface_heat_j,
                        conductivity_w_m_k,
                    }],
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            mechanical,
                            friction_source,
                            external_heat_j,
                            sphere_temperature_k: sphere_temperature,
                            source: payload,
                        },
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidHertzRollingCollinearHeatError::Interrupted(
            result.progress,
        ));
    }
    if !result.report.converged {
        return Err(RigidHertzRollingCollinearHeatError::NotConverged(
            result.report,
        ));
    }
    let trial = evaluate(&solver.x)?;
    let staged = trial.payload();
    let native = &staged.mechanical.report;
    let caloric_change = trial.report().thermal.stored_energy_change_j;
    let transfer = finite(
        caloric_change - staged.external_heat_j - native.surface_heat_j - native.bulk_heat_j,
    )?;
    let combined = finite(
        native.kinetic_after_j - native.kinetic_before_j + caloric_change - staged.external_heat_j,
    )?;
    let summed = finite(native.total_energy_residual_j + trial.report().total_energy_residual_j)?;
    if transfer.abs() > controls.heat.total_j
        || combined.abs() > controls.combined_energy_tolerance_j
        || summed.abs() > controls.combined_energy_tolerance_j
    {
        return Err(RigidHertzRollingCollinearHeatError::Invalid(
            "integrated surface/bulk work and original caloric energy do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidHertzRollingCollinearCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        friction_source: staged.friction_source,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        external_heat_j: staged.external_heat_j,
        sphere_temperature_before_k: temperature_before,
        sphere_temperature_after_k: staged.sphere_temperature_k,
        heat_transfer_residual_j: transfer,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<RigidHertzRollingCollinearHeatError<E>>,
) -> RigidHertzRollingCollinearHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => {
            RigidHertzRollingCollinearHeatError::Invalid(field)
        }
    }
}
