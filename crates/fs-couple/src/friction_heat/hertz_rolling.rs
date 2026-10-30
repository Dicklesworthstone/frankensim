//! Joint Hertz surface friction, sphere bulk rolling loss and finite heat.
//!
//! One original five-coordinate native response supplies both physical losses.
//! Hertz heat enters its massless common-temperature interface; deformation
//! heat enters the original whole-sphere caloric cell. Original coefficients
//! are resolved at interface temperature, while the original rolling card is
//! admitted at the sphere's actual bulk temperature. Nothing publishes early.

pub mod collinear;
pub mod multiple;

use super::hertz::{RigidHertzHeatContact, RigidHertzHeatSource, validate_geometry};
use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::rolling::RigidRollingHeatBinding;
use super::{FrictionHeatRoute, FrictionHeatTolerances, finite};
use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_contact::interface_binding::thermal_friction::hertz::{
    SourceHertzWrench, SourceHertzWrenchError,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz_rolling::{
    PreparedRigidHertzRollingPair, RigidHertzRollingCandidate, RigidHertzRollingContact,
    RigidHertzRollingControls, RigidHertzRollingError, RigidHertzRollingReport,
};
use fs_mbd::tangential_friction::rolling::RigidRollingHistory;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::InterfaceSystemRef;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::{RollingLossLaw, RollingPatchReceipt, RollingWorkOwnership};
use std::cell::Cell;

/// Original ordered interface, normal patch and distinct rolling work owner.
/// The caller binds these sources to the same physical interface/materials.
#[derive(Clone, Copy)]
pub struct RigidHertzRollingHeatContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceThermalCoulomb,
    pub patch: &'a RollingPatchReceipt,
    pub interface: &'a InterfaceSystemRef,
    pub rolling: &'a RollingLossLaw,
    pub excitation_frequency_hz: f64,
    pub work: &'a RollingWorkOwnership,
    /// An explicit original source coordinate, never inferred from load/area.
    pub nominal_pressure_pa: Option<f64>,
}

/// Bounded original mechanical/thermal work and independent joule gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingHeatControls {
    pub mechanical: RigidHertzRollingControls,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical scaling only, not another physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidHertzRollingHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidHertzRollingError),
    FrictionSource(SourceHertzWrenchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidHertzRollingError> for RigidHertzRollingHeatError<E> {
    fn from(error: RigidHertzRollingError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceHertzWrenchError> for RigidHertzRollingHeatError<E> {
    fn from(error: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidHertzRollingHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidHertzRollingHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidHertzRollingHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint Hertz/rolling finite heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidHertzRollingHeatError<E> {}

/// Publish both indexed bodies, original rolling history and thermal state together.
pub struct HeatedRigidHertzRollingCandidate<'m, P> {
    pub binding: RigidContactBinding,
    pub bodies_after: [RigidBodyState; 2],
    pub rolling_history: RigidRollingHistory,
    pub mechanical: RigidHertzRollingReport,
    pub thermal: DepositionSimulation<'m>,
    /// Inner surface-interface audit. Its external heat includes injected bulk
    /// rolling work; the actual outside heat for the combined system is below.
    /// Entropy classification is not supplied by this combined energy owner.
    pub heat: FrictionHeatNetworkReport,
    pub friction_source: SourceHertzWrench,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Authored and caller outside heat, measured before bulk injection.
    pub external_heat_j: f64,
    pub sphere_temperature_before_k: f64,
    pub sphere_temperature_after_k: f64,
    /// Actual caloric change minus outside heat and both physical losses.
    pub heat_transfer_residual_j: f64,
    /// Actual contacted-pair kinetic plus caloric change minus outside heat.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: RigidHertzRollingCandidate,
    friction_source: SourceHertzWrench,
    external_heat_j: f64,
    sphere_temperature_k: f64,
    source: P,
}

/// Propose simultaneous surface/bulk friction work and original finite heat.
///
/// The normal patch must retain a positive area matching the original circular
/// thermal interface. Its actual load and radius feed the pinned friction
/// source. The whole-sphere cell must occupy the corresponding A/B side, carry
/// the native sphere mass and declared material/body identities, and have its
/// centre at the sphere COM. A second active cell for that sphere refuses.
/// Other body/material correspondence remains the caller's responsibility.
///
/// One prepared five-coordinate response is reused at every thermal trial.
/// Surface friction resolves at common interface temperature; the original
/// constant-coefficient rolling card is admitted at actual sphere temperature.
/// Their two heat channels enter the original cells exactly once. The callback
/// supplies endpoint conductivities and other heat only. All accepted owners
/// stay borrowed through trial failure, exhaustion and cancellation.
///
/// Native normal/shear decoupling, normal-motion and unresolved-event limits
/// remain unchanged. This is a fixed-geometry phase, uniform sphere bulk cell
/// and common-temperature surface approximation, not a full mixed-contact
/// continuum or spatial heat solution. Thermal transport is backward Euler;
/// nonlinear convergence does not replace time/spatial refinement or establish
/// root uniqueness. No variable-viscosity, elastic-history or entropy owner is
/// invented by this adapter.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_rolling_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    history: &RigidRollingHistory,
    contact: RigidHertzRollingHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    bulk_binding: RigidRollingHeatBinding,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzRollingHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRollingCandidate<'m, P>, RigidHertzRollingHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let target = accepted_time_s + dt_s;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(RigidHertzRollingHeatError::Invalid(
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
    {
        return Err(RigidHertzRollingHeatError::Invalid(
            "synchronized finite bounded mechanical/thermal inputs required",
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
    let prepared = PreparedRigidHertzRollingPair::new(
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
    let cell =
        thermal
            .cells()
            .get(&bulk_binding.cell_id)
            .ok_or(RigidHertzRollingHeatError::Invalid(
                "inactive sphere bulk cell",
            ))?;
    if bulk_binding.body_index != sphere
        || bulk_binding.cell_id != expected_cell
        || cell.body_id != bulk_binding.thermal_body_id
        || cell.material_id != bulk_binding.material_id
        || cell.mass_kg != properties[sphere].mass()
    {
        return Err(RigidHertzRollingHeatError::Invalid(
            "whole-sphere mass and declared body/material must occupy the original route side",
        ));
    }
    for (&id, other) in thermal.cells() {
        checkpoint()?;
        if id != bulk_binding.cell_id && other.body_id == bulk_binding.thermal_body_id {
            return Err(RigidHertzRollingHeatError::Invalid(
                "whole-sphere thermal body must have exactly one active caloric cell",
            ));
        }
    }
    let point = prepared.contact_points_world()[0];
    let point = [point.x, point.y, point.z];
    let tolerance = controls.mechanical.mechanical.geometry_tolerance_m;
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
    .map_err(RigidHertzRollingHeatError::Invalid)?;
    let sign = if side == 0 { 1.0 } else { -1.0 };
    let normal = contact.frame.normal();
    let delta = std::array::from_fn::<_, 3, _>(|i| {
        0.5 * cell.geometry.start_m[i] + 0.5 * cell.geometry.end_m[i]
            - (point[i] + sign * prepared.sphere_radius_m() * normal[i])
    });
    if delta.iter().any(|v| !v.is_finite())
        || det::hypot(det::hypot(delta[0], delta[1]), delta[2]) > tolerance
    {
        return Err(RigidHertzRollingHeatError::Invalid(
            "whole-sphere caloric centre must coincide with its original native COM",
        ));
    }
    let step = FrictionHeatNetworkStep::new(thermal, target, &[route], 1, checkpoint)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidHertzRollingHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidHertzRollingHeatError::Invalid(
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
                    .map_err(RigidHertzRollingHeatError::Source)?;
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                checkpoint()?;
                // Preserve small outside heat before adding the distinct bulk
                // work source; subtracting it from a larger total can erase it.
                let external_heat_j = finite(
                    duration_s
                        * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                )?;
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
        return Err(RigidHertzRollingHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidHertzRollingHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x)?;
    let staged = trial.payload();
    let mechanical = &staged.mechanical.report;
    let caloric_change = trial.report().thermal.stored_energy_change_j;
    let heat_transfer_residual_j = finite(
        caloric_change
            - staged.external_heat_j
            - mechanical.surface_heat_j
            - mechanical.bulk_heat_j,
    )?;
    let combined = finite(
        (mechanical.kinetic_after_j - mechanical.kinetic_before_j) + caloric_change
            - staged.external_heat_j,
    )?;
    let owner_residual =
        finite(mechanical.total_energy_residual_j + trial.report().total_energy_residual_j)?;
    if heat_transfer_residual_j.abs() > controls.heat.total_j
        || combined.abs() > controls.combined_energy_tolerance_j
        || owner_residual.abs() > controls.combined_energy_tolerance_j
    {
        return Err(RigidHertzRollingHeatError::Invalid(
            "original surface/bulk work and finite caloric energy do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidHertzRollingCandidate {
        binding: contact.binding,
        bodies_after: staged.mechanical.bodies_after,
        rolling_history: staged.mechanical.rolling_history,
        mechanical: staged.mechanical.report,
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
        heat_transfer_residual_j,
        combined_energy_residual_j: combined,
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<RigidHertzRollingHeatError<E>>,
) -> RigidHertzRollingHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => RigidHertzRollingHeatError::Invalid(s),
    }
}
