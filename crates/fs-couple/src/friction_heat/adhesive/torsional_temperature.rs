//! Finite caloric temperatures feed the original loaded I/III JKR source.
//!
//! A native mechanical endpoint is staged first. This endpoint thermal phase
//! reaches that same time and changes coefficients at held approach and twist.
//! Each trial uses the two actual receiving-cell temperatures and one solved
//! massless interface temperature. Its original source supplies relaxed area
//! and stiffness; the existing variable-area heat owner supplies transport.
//!
//! JKR supplies mechanical parameter work, not a contact entropy/internal-energy
//! law. That signed work therefore has an explicitly named external driver port
//! and is NEVER debited from enthalpy or called dissipation. This is an open,
//! endpoint operator split, not closed thermoelastic feedback or a guarantee
//! of stability along an unspecified temperature/geometry path.

use super::super::{FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep};
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::motion::temperature::JkrTorsionalTemperatureReport;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeError, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatureError, AdhesiveTemperatures,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::{dd::Dd, det};
use fs_mbd::tangential_friction::adhesive::jkr_torsional::{
    PreparedRigidJkrTorsionalPair, RigidJkrTorsionalError, RigidJkrTorsionalState,
};
use fs_mbd::{RigidBodyState, Vec3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::{StepCoefficients, resolvable_increment_j};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

pub use super::super::cohesive::CohesiveHeatControls as JkrTorsionalThermalControls;

/// Original native A/B materials and actual adjoining finite receivers.
/// The route describes the staged mechanical endpoint before thermal refresh.
#[derive(Clone, Copy)]
pub struct JkrTorsionalTemperatureHeatBinding<'a> {
    pub body_ids: [u64; 2],
    pub material_cards: [&'a MaterialCard; 2],
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    pub route: &'a FrictionHeatRoute,
    /// Named physical driver accounting for the source's signed coefficient work.
    pub parameter_work_source: &'a str,
}

#[derive(Clone, Copy, Debug)]
pub struct JkrTorsionalTemperatureHeatControls {
    pub thermal: JkrTorsionalThermalControls,
    pub geometry_tolerance_m: f64,
    /// Bound on the endpoint area change relative to the incoming route.
    /// This is not a bound on an unobserved path between thermal endpoints.
    pub maximum_fractional_area_change: f64,
    pub maximum_absolute_parameter_work_j: f64,
}

/// Endpoint bulk conductivity and OTHER thermal supply only.
pub struct JkrTorsionalTemperatureHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    pub transport: StepCoefficients<P>,
}

#[derive(Debug)]
pub enum JkrTorsionalTemperatureHeatError<E> {
    Invalid(&'static str),
    Native(RigidJkrTorsionalError),
    Temperature(AdhesiveTemperatureError),
    MixedMode(JkrMixedModeError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    EnergyBalance {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl<E> From<RigidJkrTorsionalError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: RigidJkrTorsionalError) -> Self {
        Self::Native(e)
    }
}
impl<E> From<AdhesiveTemperatureError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: AdhesiveTemperatureError) -> Self {
        Self::Temperature(e)
    }
}
impl<E> From<JkrMixedModeError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: JkrMixedModeError) -> Self {
        Self::MixedMode(e)
    }
}
impl<E> From<DepositionError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<NewtonError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for JkrTorsionalTemperatureHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR torsional temperature heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for JkrTorsionalTemperatureHeatError<E> {}

pub struct JkrTorsionalTemperatureHeatReport<P> {
    pub temperature: JkrTorsionalTemperatureReport,
    pub temperatures: AdhesiveTemperatures,
    pub parameter_work_source: String,
    pub heat: FrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    /// Delta U_contact + Delta H - parameter work - other thermal supply.
    /// Native kinetic energy is unchanged by this phase.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_residual_j: f64,
}

/// Complete source, sealed native history and finite heat endpoint, unpublished.
/// Owns the successor source; a subsequent native capture borrows `source()`
/// and continues `state()` instead of initializing a new twist or clock.
pub struct JkrTorsionalTemperatureHeatCandidate<'t, P> {
    source: SourceJkrMixedMode,
    bodies: [RigidBodyState; 2],
    state: RigidJkrTorsionalState,
    thermal: DepositionSimulation<'t>,
    route: FrictionHeatRoute,
    report: JkrTorsionalTemperatureHeatReport<P>,
}
impl<'t, P> JkrTorsionalTemperatureHeatCandidate<'t, P> {
    pub const fn source(&self) -> &SourceJkrMixedMode {
        &self.source
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn state(&self) -> RigidJkrTorsionalState {
        self.state
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn route(&self) -> &FrictionHeatRoute {
        &self.route
    }
    pub const fn report(&self) -> &JkrTorsionalTemperatureHeatReport<P> {
        &self.report
    }
    /// Adopt all parts together. No mechanical impulse or geometry update occurs.
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrMixedMode,
        [RigidBodyState; 2],
        RigidJkrTorsionalState,
        DepositionSimulation<'t>,
        FrictionHeatRoute,
        JkrTorsionalTemperatureHeatReport<P>,
    ) {
        (
            self.source,
            self.bodies,
            self.state,
            self.thermal,
            self.route,
            self.report,
        )
    }
}

struct RefreshedNative {
    source: SourceJkrMixedMode,
    bodies: [RigidBodyState; 2],
    state: RigidJkrTorsionalState,
    temperatures: AdhesiveTemperatures,
    report: JkrTorsionalTemperatureReport,
}
struct Payload<P> {
    native: RefreshedNative,
    source: P,
}

/// Solve original finite temperatures, source coefficients and relaxed area.
///
/// `previous` is an already staged mechanical endpoint. Its time must be later
/// than the accepted thermal time, which this phase advances to exactly that
/// endpoint. Native geometry must already admit the state; pending compatible
/// travel is retained, never consumed here. The incoming route must describe
/// this same mechanical endpoint before changing source temperature. A moving
/// face, stale normal gap or different area needs its original geometry owner.
///
/// Thermal transport uses the endpoint area and actual ambient coverage, not an
/// invented mean-area path. Both finite-body temperatures and the solved common
/// interface Kelvin refresh their own original selected claims independently.
/// Incoming bulk source T coordinates must already agree with the accepted
/// receiving-cell temperatures, including for temperature-independent claims.
/// A source without an original T coordinate cannot use this feedback owner.
/// The callback receives those refreshed sources and endpoint covered geometry.
/// It stages original conductivity/other-supply data; it must not publish state,
/// reinsert parameter work as heat or override this route's conduction link.
///
/// Each trial internally authenticates native/source continuation. The signed
/// full normal/torsional/surface parameter work retains its named driver port.
/// No frictional heat or mechanical entropy is inferred. Losing an attached
/// equilibrium, source support, receiver identity or any acceptance gate refuses
/// the entire candidate, including on final-trial cancellation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_jkr_torsional_temperature_heat<'t, P, E>(
    native: &PreparedRigidJkrTorsionalPair<'_>,
    previous: RigidJkrTorsionalState,
    thermal: &DepositionSimulation<'t>,
    binding: JkrTorsionalTemperatureHeatBinding<'_>,
    controls: JkrTorsionalTemperatureHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        AdhesiveTemperatures,
        &SourceJkrMixedMode,
    ) -> Result<JkrTorsionalTemperatureHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<JkrTorsionalTemperatureHeatCandidate<'t, P>, JkrTorsionalTemperatureHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let c = controls.thermal;
    c.newton.validate()?;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(JkrTorsionalTemperatureHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    if !previous.time_s().is_finite()
        || previous.time_s() <= thermal.time_s()
        || cells == 0
        || cells > c.maximum_cells
        || c.maximum_interfaces == 0
        || dimension > c.maximum_unknowns
        || c.maximum_trials == 0
        || c.maximum_newton_iterations == 0
        || !c.difference_step_j.is_finite()
        || c.difference_step_j <= 0.0
        || !c.temperature_scale_j_per_k.is_finite()
        || c.temperature_scale_j_per_k <= 0.0
        || !temperature_guess_k.is_finite()
        || temperature_guess_k <= 0.0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || binding.parameter_work_source.trim().is_empty()
        || [
            controls.geometry_tolerance_m,
            controls.maximum_fractional_area_change,
            controls.maximum_absolute_parameter_work_j,
            c.heat.cell_j,
            c.heat.interface_j,
            c.heat.thermal_j,
            c.heat.total_j,
            c.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(JkrTorsionalTemperatureHeatError::Invalid(
            "bounded thermal endpoint and explicit parameter-work port required",
        ));
    }
    validate_receivers(native, previous, thermal, binding)?;
    let contact = native.contact();
    let normal = contact.source.normal_source();
    let normal_state = normal
        .evaluate_attached(previous.coordinates()[0])
        .map_err(AdhesiveTemperatureError::from)?
        .state();
    let step = FrictionHeatStep::new(thermal, previous.time_s(), binding.route, check)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= c.maximum_trials {
            return Err(JkrTorsionalTemperatureHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(JkrTorsionalTemperatureHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let interface_k = temperature_guess_k + x[cells] / c.temperature_scale_j_per_k;
        let increments: Vec<_> = thermal
            .cells()
            .values()
            .zip(&x[..cells])
            .map(|(cell, &q)| resolvable_increment_j(cell.specific_h_j_kg, q, cell.mass_kg))
            .collect();
        let prepare = |endpoint: &DepositionSimulation<'t>, interface_k| {
            let ids = binding.route.cells();
            let temperatures = AdhesiveTemperatures {
                surface_a_k: endpoint.temperature_k(ids[0])?,
                surface_b_k: endpoint.temperature_k(ids[1])?,
                interface_k,
            };
            let refresh = normal.refresh_at_held_indentation(
                normal_state,
                binding.material_cards,
                temperatures,
            )?;
            let source =
                SourceJkrMixedMode::new(refresh.source().clone(), contact.source.envelope())?;
            let continuation = native.continue_at_temperature(previous, &refresh, &source, gate)?;
            let report = *continuation.report();
            let bodies = *continuation.bodies();
            validate_geometry(native, bodies, binding.route, controls.geometry_tolerance_m)?;
            let area = report.after.contact_area_m2;
            let fractional =
                (area - binding.route.contact_area_m2()).abs() / binding.route.contact_area_m2();
            if !fractional.is_finite()
                || fractional > controls.maximum_fractional_area_change
                || report.parameter_work_j.abs() > controls.maximum_absolute_parameter_work_j
            {
                return Err(JkrTorsionalTemperatureHeatError::Invalid(
                    "source area or parameter work exceeds the explicit endpoint allowance",
                ));
            }
            let state = continuation.state();
            drop(continuation);
            Ok((
                area,
                RefreshedNative {
                    source,
                    bodies,
                    state,
                    temperatures,
                    report,
                },
            ))
        };
        let transport =
            |endpoint: &DepositionSimulation<'t>, _temperature, staged: RefreshedNative| {
                let input = resolve(endpoint, staged.temperatures, &staged.source)
                    .map_err(JkrTorsionalTemperatureHeatError::Source)?;
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = input.transport;
                Ok(FrictionHeatInput {
                    dissipated_j: 0.0,
                    conductivity_w_m_k: input.conductivity_w_m_k,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Payload {
                            native: staged,
                            source: payload,
                        },
                    },
                })
            };
        if entropy {
            step.evaluate_with_area_and_entropy(&increments, interface_k, prepare, transport, check)
        } else {
            step.evaluate_with_area(&increments, interface_k, prepare, transport, check)
        }
    };
    let problem = FallibleFiniteDifference::new(dimension, c.difference_step_j, |x: &[f64]| {
        evaluate(x, false).map(|trial| trial.residual_j().to_vec())
    })
    .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
    let construction = NewtonKrylovState::new(&problem, initial, c.newton);
    if let Some(e) = problem.take_failure() {
        return Err(numerical(e));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, c.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(e) = problem.take_failure() {
        return Err(numerical(e));
    }
    check()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(JkrTorsionalTemperatureHeatError::Interrupted(
            result.progress,
        ));
    }
    if !result.report.converged {
        return Err(JkrTorsionalTemperatureHeatError::NotConverged(
            result.report,
        ));
    }
    let trial = evaluate(&solver.x, c.audit_entropy)?;
    let temperature = &trial.payload().native.report;
    let heat = trial.report();
    let combined = (Dd::from_f64(temperature.storage_change_j)
        + Dd::from_f64(heat.thermal.stored_energy_change_j)
        - Dd::from_f64(temperature.parameter_work_j)
        - Dd::from_f64(heat.external_heat_j))
    .to_f64();
    let absolute_stages = (Dd::from_f64(temperature.energy_residual_j.abs())
        + Dd::from_f64(heat.total_energy_residual_j.abs()))
    .to_f64();
    if heat.dissipated_j != 0.0
        || [combined, absolute_stages]
            .iter()
            .any(|r| !r.is_finite() || r.abs() > c.combined_energy_tolerance_j)
    {
        return Err(JkrTorsionalTemperatureHeatError::EnergyBalance {
            residual_j: combined,
            absolute_stage_residual_j: absolute_stages,
        });
    }
    let (thermal, payload, heat, route) = trial.accept(c.heat, check)?;
    if thermal.time_s() != previous.time_s() {
        return Err(JkrTorsionalTemperatureHeatError::Invalid(
            "native and thermal endpoint clocks differ",
        ));
    }
    let staged = payload.native;
    check()?;
    Ok(JkrTorsionalTemperatureHeatCandidate {
        source: staged.source,
        bodies: staged.bodies,
        state: staged.state,
        thermal,
        route,
        report: JkrTorsionalTemperatureHeatReport {
            temperature: staged.report,
            temperatures: staged.temperatures,
            parameter_work_source: binding.parameter_work_source.to_owned(),
            heat,
            source: payload.source,
            nonlinear: result.report,
            thermal_trials: trials.get(),
            combined_energy_residual_j: combined,
            absolute_stage_residual_j: absolute_stages,
        },
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<JkrTorsionalTemperatureHeatError<E>>,
) -> JkrTorsionalTemperatureHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(m) => JkrTorsionalTemperatureHeatError::Invalid(m),
    }
}

fn validate_receivers<E>(
    native: &PreparedRigidJkrTorsionalPair<'_>,
    previous: RigidJkrTorsionalState,
    thermal: &DepositionSimulation<'_>,
    binding: JkrTorsionalTemperatureHeatBinding<'_>,
) -> Result<(), JkrTorsionalTemperatureHeatError<E>> {
    binding.route.validate(thermal)?;
    let elastic = native
        .contact()
        .source
        .normal_source()
        .normal_state()
        .elastic();
    let identities = [
        elastic.surface_a_state().resolved().card_identity(),
        elastic.surface_b_state().resolved().card_identity(),
    ];
    let source_temperatures = [
        elastic.surface_a_state().resolved().state_coordinate("T"),
        elastic.surface_b_state().resolved().state_coordinate("T"),
    ];
    let area = previous.motion_state().contact().contact_area_m2;
    if binding.body_ids[0] == binding.body_ids[1]
        || (area - binding.route.contact_area_m2()).abs()
            > 64.0 * f64::EPSILON * area.max(binding.route.contact_area_m2())
    {
        return Err(JkrTorsionalTemperatureHeatError::Invalid(
            "original ordered native bodies and endpoint area required",
        ));
    }
    for (i, id) in binding.route.cells().into_iter().enumerate() {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(JkrTorsionalTemperatureHeatError::Invalid(
                "missing original thermal receiver",
            ))?;
        if cell.body_id != binding.body_ids[i]
            || binding.caloric[i].material_card_identity() != identities[i]
            || !thermal.uses_material_law(cell.material_id, binding.caloric[i])
        {
            return Err(JkrTorsionalTemperatureHeatError::Invalid(
                "original ordered registered caloric materials required",
            ));
        }
        let caloric_temperature = thermal.temperature_k(id)?;
        let Some(source_temperature) = source_temperatures[i] else {
            return Err(JkrTorsionalTemperatureHeatError::Invalid(
                "original bulk source must own its absolute-temperature T coordinate",
            ));
        };
        // Only representation roundoff is admitted, not an independently stale
        // source temperature or a caller-selected thermomechanical tolerance.
        if !source_temperature.is_finite()
            || source_temperature <= 0.0
            || (source_temperature - caloric_temperature).abs()
                > 64.0 * f64::EPSILON * source_temperature.abs().max(caloric_temperature.abs())
        {
            return Err(JkrTorsionalTemperatureHeatError::Invalid(
                "incoming bulk source temperatures differ from their accepted caloric cells",
            ));
        }
    }
    Ok(())
}

fn validate_geometry<E>(
    native: &PreparedRigidJkrTorsionalPair<'_>,
    bodies: [RigidBodyState; 2],
    route: &FrictionHeatRoute,
    tolerance: f64,
) -> Result<(), JkrTorsionalTemperatureHeatError<E>> {
    let g = route.geometry;
    let centre = [g[0].start_m, g[0].end_m]
        .into_iter()
        .find(|p| *p == g[1].start_m || *p == g[1].end_m)
        .ok_or(JkrTorsionalTemperatureHeatError::Invalid(
            "original thermal face does not adjoin",
        ))?;
    let n = native.contact().frame.normal();
    let normal = Vec3::new(n[0], n[1], n[2]);
    let arms = native.reaction_arms_body();
    for i in 0..2 {
        let point = bodies[i]
            .pose()
            .point_world_from_body(arms[i])
            .map_err(RigidJkrTorsionalError::from)?;
        let axis = Vec3::new(
            g[i].end_m[0] - g[i].start_m[0],
            g[i].end_m[1] - g[i].start_m[1],
            g[i].end_m[2] - g[i].start_m[2],
        );
        let interior = if g[i].start_m == centre {
            g[i].end_m
        } else {
            g[i].start_m
        };
        let inward = Vec3::new(
            0.5 * (interior[0] - centre[0]),
            0.5 * (interior[1] - centre[1]),
            0.5 * (interior[2] - centre[2]),
        );
        let depth = inward.dot(normal) * if i == 0 { -1.0 } else { 1.0 };
        if magnitude(point.sub(Vec3::new(centre[0], centre[1], centre[2]))) > tolerance
            || magnitude(axis.cross(normal)) > tolerance
            || !depth.is_finite()
            || depth <= tolerance
        {
            return Err(JkrTorsionalTemperatureHeatError::Invalid(
                "native reaction plane and ordered thermal face must coincide",
            ));
        }
    }
    Ok(())
}
fn magnitude(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
