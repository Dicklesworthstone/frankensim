//! Original directional material friction and finite rigid/thermal endpoints.
//!
//! The existing rigid owner supplies full tensor mobility and physical work.
//! Its source callback closes every declared speed/path coordinate against the
//! actual matrix-resolved slip. The existing zero-capacity interface then
//! distributes that work by the original cells' geometry and conductivities.
//! Both surface temperatures and an explicitly declared interface temperature
//! can feed back into the SAME selected coefficient claims at each thermal trial.
//!
//! This is one held-load, frozen-geometry midpoint phase. Thermal cells can be
//! local contact volumes, rather than complete rigid-body caloric inventories.
//! Their prismatic discretization and perfect-contact route retain the limits
//! of `FrictionHeatStep`. Pose/normal evolution and unresolved slip events stay
//! with the trajectory owner. Endpoint thermal coupling needs time refinement;
//! neither nonlinear convergence nor a source-motion root proves uniqueness.

/// Continuous directional trajectories, including finite arrest, in the same
/// original material/finite-enthalpy transaction.
pub mod trajectory;

/// Static hold and driven directional friction with the original heat network.
pub mod driven;

/// Original source claims through moving-sphere impacts and sustained contact.
pub mod sphere_plane;

/// Original interface sources through simultaneous normal/directional impacts.
pub mod impact;

/// Original interface sources through joint-constrained simultaneous impacts.
pub mod articulated;

use super::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};
use fs_contact::interface_binding::thermal_friction::FrictionQueryAxisKind;
use fs_contact::interface_binding::thermal_friction::elliptic::motion::{
    MotionQueryError, MotionQuerySolve,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_exec::CancelGate;
use fs_matdb::QueryPoint;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialError;
use fs_math::det;
use fs_mbd::tangential_friction::directional::{
    RigidDirectionalContact, RigidDirectionalControls, RigidDirectionalError,
    RigidDirectionalReport, apply_rigid_directional_with,
};
use fs_mbd::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// Original ordered surface material/enthalpy objects and temperature ownership.
#[derive(Clone, Copy, Debug)]
pub struct DirectionalHeatMaterial<'a> {
    pub point: ThermalInterfacePoint<'a>,
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    /// Optional existing absolute-temperature axis owned by the solved common
    /// interface. It must be distinct from all surface/mean/phase axes.
    pub interface_temperature_axis: Option<&'a str>,
}

/// Material transport at the actual two finite-cell endpoint states.
pub struct DirectionalHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    /// Other heat and link conductances only; the bound route owns its link and
    /// the mechanical owner supplies friction heat exactly once.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct DirectionalHeatSolve {
    pub mechanical: RigidDirectionalControls,
    pub motion: MotionQuerySolve,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Newton coordinate scale [J/K], never a physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum DirectionalHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidDirectionalError<MotionQueryError>),
    Dynamics(DynamicsError),
    Material(ThermalMaterialError),
    Phase(PhaseStateError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for DirectionalHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<DynamicsError> for DirectionalHeatError<E> {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl<E> From<ThermalMaterialError> for DirectionalHeatError<E> {
    fn from(e: ThermalMaterialError) -> Self {
        Self::Material(e)
    }
}
impl<E> From<PhaseStateError> for DirectionalHeatError<E> {
    fn from(e: PhaseStateError) -> Self {
        Self::Phase(e)
    }
}
impl<E> From<NewtonError> for DirectionalHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E> From<RigidDirectionalError<MotionQueryError>> for DirectionalHeatError<E> {
    fn from(e: RigidDirectionalError<MotionQueryError>) -> Self {
        Self::Mechanical(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for DirectionalHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional rigid/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DirectionalHeatError<E> {}

/// Source, body and thermal results from the same accepted endpoint.
pub struct HeatedDirectionalReport<P> {
    pub mechanical: RigidDirectionalReport<SourceThermalEllipticState>,
    pub heat: FrictionHeatReport,
    pub phases: [EquilibriumPhaseState; 2],
    pub pressure_pa: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub source_evaluations: usize,
    pub queried_distance_m: f64,
    pub actual_distance_m: f64,
    pub source_motion_residual_m: f64,
    /// Delta(K + U_contact + H_cells) - Q_external. Physical friction cancels
    /// internally; it is not an additional external energy supply.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    bodies: [RigidBodyState; 2],
    mechanical: RigidDirectionalReport<SourceThermalEllipticState>,
    phases: [EquilibriumPhaseState; 2],
    source: P,
    motion: [f64; 3],
}

/// Advance source-owned directional friction, original bodies and finite heat.
///
/// `contact` A/B must follow `route` and `material` A/B. Its temperature and
/// pressure fields must be `None`: this solve owns temperature and derives
/// pressure from held normal load / actual route area. Both body contact points
/// must meet the shared thermal face. Only the selected two bodies are staged.
///
/// The original source checkpoint, body momenta and thermal state publish
/// together after all gates. Source callbacks stage receipts without changing
/// external state. Motion-dependent source brackets are explicit physical
/// inputs; full-matrix slip need not lie below released Euclidean displacement.
/// Work is bounded by thermal trial count times the source evaluation budget.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_source_directional_heat<'m, P, E>(
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
    energy_guess_j: &[f64],
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedDirectionalReport<P>, DirectionalHeatError<E>> {
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
        ));
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
        return Err(DirectionalHeatError::Invalid("normal pressure underflow"));
    }
    let initial_bodies = [bodies[binding.body_a], bodies[binding.body_b]];
    let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
    let mut pair_contact = contact;
    pair_contact.binding.body_a = 0;
    pair_contact.binding.body_b = 1;
    pair_contact.nominal_pressure_pa = Some(pressure);
    let step = FrictionHeatStep::new(thermal, target_s, route, checkpoint)?;
    let dt_s = step.duration_s();
    let cells = step.cell_ids().len();
    let trials = Cell::new(0usize);
    let source_evaluations = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(DirectionalHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != cells + 1 {
            return Err(DirectionalHeatError::Invalid(
                "incomplete directional thermal iterate",
            ));
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
                let mechanical = apply_rigid_directional_with(
                    &mut next_bodies,
                    &pair_properties,
                    resolved_contact,
                    &mut next_history,
                    dt_s,
                    controls.mechanical,
                    |accepted, input, compliance, mode| {
                        let resolved = law.advance_motion_consistent(
                            accepted,
                            input,
                            &query,
                            mode,
                            compliance,
                            controls.motion,
                            || gate.is_requested(),
                        )?;
                        source_evaluations.set(
                            source_evaluations
                                .get()
                                .checked_add(resolved.evaluations)
                                .ok_or(MotionQueryError::Invalid(
                                    "source evaluation count overflow",
                                ))?,
                        );
                        motion = Some([
                            resolved.queried_distance_m,
                            resolved.actual_distance_m,
                            resolved.distance_residual_m,
                        ]);
                        Ok::<_, MotionQueryError>(resolved.step)
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
                Ok::<_, DirectionalHeatError<E>>(FrictionHeatInput {
                    dissipated_j: mechanical.contact.dissipated_j,
                    conductivity_w_m_k,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            bodies: next_bodies,
                            mechanical,
                            phases,
                            source: payload,
                            motion: motion.ok_or(DirectionalHeatError::Invalid(
                                "missing source motion receipt",
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
        return Err(DirectionalHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(DirectionalHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x)?;
    let mechanical = &trial.payload().mechanical;
    let heat = trial.report();
    let combined = finite(
        mechanical.kinetic_energy_change_j
            + mechanical.contact.stored_change_j
            + heat.thermal.stored_energy_change_j
            - heat.external_heat_j,
    )?;
    if combined.abs() > controls.combined_energy_tolerance_j {
        return Err(DirectionalHeatError::Invalid(
            "native kinetic, contact and thermal energies do not balance",
        ));
    }
    let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    let next_history = staged.mechanical.contact.state.clone();
    checkpoint()?;
    drop(problem);
    // All fallible work and callbacks precede this joint publication boundary.
    bodies[binding.body_a] = staged.bodies[0];
    bodies[binding.body_b] = staged.bodies[1];
    *history = next_history;
    *thermal = next_thermal;
    Ok(HeatedDirectionalReport {
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

impl DirectionalHeatMaterial<'_> {
    fn validate<E>(
        &self,
        law: &SourceThermalEllipticCoulomb,
    ) -> Result<(), DirectionalHeatError<E>> {
        if let Some(axis) = self.interface_temperature_axis {
            let absolute = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            ));
            if axis.trim().is_empty()
                || !self.point.context.axes().contains_key(axis)
                || self.point.temperature_axes.contains(&axis)
                || self.point.mean_temperature_axis == Some(axis)
                || self.point.liquid_fraction_axes.contains(&Some(axis))
                || self
                    .point
                    .context
                    .axis_quantities()
                    .get(axis)
                    .is_some_and(|q| {
                        *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)
                    })
            {
                return Err(DirectionalHeatError::Invalid(
                    "solved interface temperature needs a distinct existing Kelvin axis",
                ));
            }
        }
        let owns_temperature = |axis: &str| {
            self.point.is_thermal_axis(axis, Temperature::DIMS)
                || self.interface_temperature_axis == Some(axis)
        };
        for source in law.principal_sources() {
            if source.interface().card_identity() != self.point.interface.content_hash()
                || !owns_temperature(&source.axes().temperature)
                || source.query_axes().iter().any(|axis| {
                    axis.kind == FrictionQueryAxisKind::Temperature && !owns_temperature(&axis.name)
                })
            {
                return Err(DirectionalHeatError::Invalid(
                    "original interface and all evolving temperature owners must agree",
                ));
            }
        }
        Ok(())
    }

    fn sample<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        route: &FrictionHeatRoute,
        law: &SourceThermalEllipticCoulomb,
        interface_temperature: f64,
        pressure: f64,
    ) -> Result<([EquilibriumPhaseState; 2], QueryPoint), DirectionalHeatError<E>> {
        let phase = |i: usize| {
            let cell =
                thermal
                    .cells()
                    .get(&route.cells()[i])
                    .ok_or(DirectionalHeatError::Invalid(
                        "inactive original thermal surface",
                    ))?;
            if !thermal.uses_material_law(cell.material_id, self.caloric[i]) {
                return Err(DirectionalHeatError::Invalid(
                    "directional heat must retain original caloric objects",
                ));
            }
            let state = self.caloric[i].state_at_specific_enthalpy(cell.specific_h_j_kg)?;
            if state.liquid_mass_fraction() != 0.0 {
                return Err(DirectionalHeatError::Invalid(
                    "dry directional contact refuses liquid surfaces",
                ));
            }
            Ok::<_, DirectionalHeatError<E>>(state)
        };
        let phases = [phase(0)?, phase(1)?];
        let mut point = self.point.at(phases)?;
        let coordinate = |q: QueryPoint, axis: &str, value: f64| {
            q.with(axis, value).map_err(|e| {
                DirectionalHeatError::Material(ThermalMaterialError::Coordinate(e.to_string()))
            })
        };
        if let Some(axis) = self.interface_temperature_axis {
            point = coordinate(point, axis, interface_temperature)?;
        }
        let source = &law.principal_sources()[0];
        if let Some(axis) = &source.axes().nominal_pressure {
            point = coordinate(point, axis, pressure)?;
        }
        for axis in source.query_axes() {
            if axis.kind == FrictionQueryAxisKind::NominalPressure {
                point = coordinate(point, &axis.name, pressure)?;
            }
        }
        Ok((phases, point))
    }
}

fn finite<E>(x: f64) -> Result<f64, DirectionalHeatError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DirectionalHeatError::Invalid(
            "unrepresentable directional work or coordinate",
        ))
    }
}
fn magnitude(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
fn validate_geometry<E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    route: &FrictionHeatRoute,
    tolerance_m: f64,
) -> Result<(), DirectionalHeatError<E>> {
    if !tolerance_m.is_finite() || tolerance_m < 0.0 {
        return Err(DirectionalHeatError::Invalid(
            "finite nonnegative geometry tolerance required",
        ));
    }
    let geometry = route.geometry;
    let center = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or(DirectionalHeatError::Invalid(
            "thermal face has no shared centre",
        ))?;
    let vector = |v: [f64; 3]| Vec3::new(v[0], v[1], v[2]);
    let normal = vector(contact.frame.normal());
    let axis = vector(geometry[0].end_m).sub(vector(geometry[0].start_m));
    if magnitude(axis.cross(normal)) > tolerance_m {
        return Err(DirectionalHeatError::Invalid(
            "thermal face and contact tangent frame disagree",
        ));
    }
    for (i, arm) in [
        (contact.binding.body_a, contact.binding.arm_a_body),
        (contact.binding.body_b, contact.binding.arm_b_body),
    ] {
        let point = bodies[i].point_kinematics(properties[i], arm)?.point_world;
        if magnitude(point.sub(vector(center))) > tolerance_m {
            return Err(DirectionalHeatError::Invalid(
                "rigid contact and thermal face centres disagree",
            ));
        }
    }
    Ok(())
}
fn numerical<E>(error: FiniteDifferenceError<DirectionalHeatError<E>>) -> DirectionalHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(message) => DirectionalHeatError::Invalid(message),
    }
}
