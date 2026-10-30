//! Sourced curved sliding and finite arrest coupled to the original heat cells.
//!
//! The native trajectory integrates held coefficients through the actual
//! stopping event, retaining its curved displacement, unsigned path and physical
//! traction-slip work separately. Every thermal trial resolves the same original
//! source claims at its endpoint temperatures and at the trajectory's mean
//! speed/path, then closes the source-motion residual. The zero-capacity common
//! interface distributes only that physical work to the existing finite cells.
//!
//! This is an endpoint thermal/source discretization: coefficients are held
//! within each mechanical phase, including its arrest event. Subdivision is
//! required for changing material coefficients or the timing of heat release.
//! The native quadrature estimates and nonlinear residuals do not certify a
//! continuously varying source law, a unique source root or spatial convergence.
//! Poses, contact frame, normal load and body inertia remain fixed.

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::trajectory::SourceRelaxationError;
use fs_mbd::tangential_friction::directional::trajectory::{
    RigidDirectionalTrajectoryError, RigidDirectionalTrajectoryReport,
    apply_rigid_directional_trajectory_with,
};
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

#[derive(Debug)]
pub enum DirectionalTrajectoryHeatError<E> {
    Coupling(DirectionalHeatError<E>),
    Mechanical(RigidDirectionalTrajectoryError<SourceRelaxationError>),
}
impl<E> From<DirectionalHeatError<E>> for DirectionalTrajectoryHeatError<E> {
    fn from(error: DirectionalHeatError<E>) -> Self {
        Self::Coupling(error)
    }
}
impl<E> From<DepositionError> for DirectionalTrajectoryHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Coupling(DirectionalHeatError::Thermal(error))
    }
}
impl<E> From<NewtonError> for DirectionalTrajectoryHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Coupling(DirectionalHeatError::Newton(error))
    }
}
impl<E> From<RigidDirectionalTrajectoryError<SourceRelaxationError>>
    for DirectionalTrajectoryHeatError<E>
{
    fn from(error: RigidDirectionalTrajectoryError<SourceRelaxationError>) -> Self {
        Self::Mechanical(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for DirectionalTrajectoryHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sourced directional trajectory/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DirectionalTrajectoryHeatError<E> {}

/// One jointly accepted continuous mechanical phase and finite thermal endpoint.
pub struct HeatedDirectionalTrajectoryReport<P> {
    pub mechanical: RigidDirectionalTrajectoryReport,
    pub heat: FrictionHeatReport,
    pub phases: [EquilibriumPhaseState; 2],
    pub pressure_pa: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub source_evaluations: usize,
    pub queried_distance_m: f64,
    /// Unsigned travel of the native curve, not its net displacement norm.
    pub actual_distance_m: f64,
    pub source_motion_residual_m: f64,
    /// Delta(K + H_cells) - Q_external; friction has one internal debit/credit.
    pub combined_energy_residual_j: f64,
}

struct StagedTrajectory<P> {
    bodies: [RigidBodyState; 2],
    history: SourceThermalEllipticState,
    mechanical: RigidDirectionalTrajectoryReport,
    phases: [EquilibriumPhaseState; 2],
    source: P,
    motion: [f64; 3],
}

/// Advance a sourced finite-body slide through arrest while closing finite heat.
///
/// The surface/pressure ownership, original caloric objects and thermal face
/// geometry are exactly those of `apply_source_directional_heat`. Mechanical
/// travel follows the native full-matrix curve. The source speed is its mean
/// unsigned speed over this phase, including any time after arrest; source path
/// increments by that curve's travel. Endpoint temperature and mean-speed
/// coefficient sampling require time refinement. Convergence depends on the
/// source law's regularity and a consistent physical branch; this phase does
/// not integrate a pointwise varying source law.
///
/// All source and quadrature trials start from the original accepted bodies and
/// history. Exactly two bodies, their source history and the original thermal
/// state publish together after heat/impulse/energy checks and cancellation.
/// Thermal source callbacks must stage their payload without external mutation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_source_directional_trajectory_heat<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut SourceThermalEllipticState,
    law: &SourceThermalEllipticCoulomb,
    contact: RigidDirectionalContact,
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    material: DirectionalHeatMaterial<'_>,
    dt_s: f64,
    controls: DirectionalHeatSolve,
    relaxation: EllipticRelaxationControls,
    energy_guess_j: &[f64],
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedDirectionalTrajectoryReport<P>, DirectionalTrajectoryHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let binding = contact.binding;
    let target_s = thermal.time_s() + dt_s;
    if bodies.len() != properties.len()
        || bodies.len() > controls.mechanical.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target_s.is_finite()
        || target_s <= thermal.time_s()
        || contact.contact_temperature_k.is_some()
        || contact.nominal_pressure_pa.is_some()
        || !contact.normal_load_n.is_finite()
        || contact.normal_load_n < 0.0
        || controls.maximum_cells == 0
        || thermal.cells().len() > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || [
            controls.difference_step_j,
            controls.temperature_scale_j_per_k,
            interface_temperature_guess_k,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !controls.combined_energy_tolerance_j.is_finite()
        || controls.combined_energy_tolerance_j < 0.0
        || energy_guess_j.len() != thermal.cells().len()
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(DirectionalHeatError::Invalid(
            "bounded original bodies, thermal coordinates and exclusive temperature/pressure ownership required",
        )
        .into());
    }
    material.validate(law)?;
    validate_geometry(
        bodies,
        properties,
        contact,
        route,
        controls.mechanical.geometry_tolerance_m,
    )?;
    let pressure = finite(contact.normal_load_n / route.contact_area_m2())?;
    if contact.normal_load_n > 0.0 && pressure == 0.0 {
        return Err(DirectionalHeatError::Invalid("normal pressure underflow").into());
    }
    let initial_bodies = [bodies[binding.body_a], bodies[binding.body_b]];
    let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
    let mut pair_contact = contact;
    pair_contact.binding.body_a = 0;
    pair_contact.binding.body_b = 1;
    pair_contact.nominal_pressure_pa = Some(pressure);
    let step = FrictionHeatStep::new(thermal, target_s, route, checkpoint)?;
    // The represented clock difference owns both mechanical and thermal time.
    let dt_s = step.duration_s();
    let cells = step.cell_ids().len();
    let trials = Cell::new(0usize);
    let source_evaluations = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(DirectionalHeatError::TrialBudget.into());
        }
        trials.set(trials.get() + 1);
        if x.len() != cells + 1 {
            return Err(DirectionalHeatError::Invalid(
                "incomplete directional trajectory thermal iterate",
            )
            .into());
        }
        let temperature =
            finite(interface_temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k)?;
        step.evaluate(
            &x[..cells],
            temperature,
            |endpoint, interface_temperature| {
                let (phases, query) =
                    material.sample(endpoint, route, law, interface_temperature, pressure)?;
                let primary = &law.principal_sources()[0].axes().temperature;
                let contact_temperature =
                    *query
                        .axes()
                        .get(primary)
                        .ok_or(DirectionalHeatError::Invalid(
                            "missing original contact temperature",
                        ))?;
                let mut resolved_contact = pair_contact;
                resolved_contact.contact_temperature_k = Some(contact_temperature);
                let mut next_bodies = initial_bodies;
                let mut next_history = history.clone();
                let mut motion = None;
                let mechanical = apply_rigid_directional_trajectory_with(
                    &mut next_bodies,
                    &pair_properties,
                    resolved_contact,
                    &mut next_history,
                    dt_s,
                    controls.mechanical,
                    |accepted, input| {
                        let resolved = law.relax_motion_consistent(
                            accepted,
                            input,
                            &query,
                            controls.motion,
                            relaxation,
                            &mut || gate.is_requested(),
                        )?;
                        source_evaluations.set(
                            source_evaluations
                                .get()
                                .checked_add(resolved.evaluations)
                                .ok_or(SourceRelaxationError::Query(MotionQueryError::Invalid(
                                    "source evaluation count overflow",
                                )))?,
                        );
                        motion = Some([
                            resolved.queried_distance_m,
                            resolved.actual_distance_m,
                            resolved.distance_residual_m,
                        ]);
                        Ok::<_, SourceRelaxationError>((resolved.relaxation, resolved.state))
                    },
                    gate,
                )?;
                let DirectionalHeatSource {
                    conductivity_w_m_k,
                    transport,
                } = resolve(endpoint, phases).map_err(DirectionalHeatError::Source)?;
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                Ok::<_, DirectionalTrajectoryHeatError<E>>(FrictionHeatInput {
                    dissipated_j: mechanical.relaxation.dissipated_j,
                    conductivity_w_m_k,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: StagedTrajectory {
                            bodies: next_bodies,
                            history: next_history,
                            mechanical,
                            phases,
                            source: payload,
                            motion: motion.ok_or(DirectionalHeatError::Invalid(
                                "missing source trajectory receipt",
                            ))?,
                        },
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(cells + 1, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(trajectory_numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(trajectory_numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(trajectory_numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(DirectionalHeatError::Interrupted(result.progress).into());
    }
    if !result.report.converged {
        return Err(DirectionalHeatError::NotConverged(result.report).into());
    }
    let trial = evaluate(&solver.x)?;
    let mechanical = &trial.payload().mechanical;
    let heat = trial.report();
    let combined = finite(
        mechanical.kinetic_energy_change_j + heat.thermal.stored_energy_change_j
            - heat.external_heat_j,
    )?;
    if combined.abs() > controls.combined_energy_tolerance_j {
        return Err(DirectionalHeatError::Invalid(
            "native trajectory kinetic and finite thermal energies do not balance",
        )
        .into());
    }
    let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    drop(problem);
    // All fallible work and callbacks precede this joint publication boundary.
    bodies[binding.body_a] = staged.bodies[0];
    bodies[binding.body_b] = staged.bodies[1];
    *history = staged.history;
    *thermal = next_thermal;
    Ok(HeatedDirectionalTrajectoryReport {
        mechanical: staged.mechanical,
        heat,
        phases: staged.phases,
        pressure_pa: pressure,
        source: staged.source,
        nonlinear: result.report,
        thermal_trials: trials.get(),
        source_evaluations: source_evaluations.get(),
        queried_distance_m: staged.motion[0],
        actual_distance_m: staged.motion[1],
        source_motion_residual_m: staged.motion[2],
        combined_energy_residual_j: combined,
    })
}

fn trajectory_numerical<E>(
    error: FiniteDifferenceError<DirectionalTrajectoryHeatError<E>>,
) -> DirectionalTrajectoryHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(message) => DirectionalHeatError::Invalid(message).into(),
    }
}
