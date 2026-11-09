//! Shared finite-temperature/area solve for original loaded adhesive models.
//!
//! Native mode adapters contribute only their sealed source continuation. The
//! original variable-area heat route, finite receivers, massless interface
//! temperature, nonlinear solve and independent work/heat acceptance remain
//! one implementation. Mechanical coefficient work has an explicit external
//! supply; it is never inferred caloric heat or contact entropy.

use super::super::{FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep};
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeError, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatureError, AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::{dd::Dd, det};
use fs_mbd::{DynamicsError, RigidBodyState, Vec3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::{StepCoefficients, resolvable_increment_j};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

pub use super::super::cohesive::CohesiveHeatControls as JkrAdhesiveThermalControls;

/// Original native A/B materials and actual adjoining finite receivers.
/// The route describes the staged mechanical endpoint before thermal refresh.
#[derive(Clone, Copy)]
pub struct JkrAdhesiveTemperatureHeatBinding<'a> {
    pub body_ids: [u64; 2],
    pub material_cards: [&'a MaterialCard; 2],
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    pub route: &'a FrictionHeatRoute,
    /// Named physical driver accounting for the source's signed coefficient work.
    pub parameter_work_source: &'a str,
}

#[derive(Clone, Copy, Debug)]
pub struct JkrAdhesiveTemperatureHeatControls {
    pub thermal: JkrAdhesiveThermalControls,
    pub geometry_tolerance_m: f64,
    /// Bound on the endpoint area change relative to the incoming route.
    /// This is not a bound on an unobserved path between thermal endpoints.
    pub maximum_fractional_area_change: f64,
    pub maximum_absolute_parameter_work_j: f64,
}

/// Endpoint bulk conductivity and OTHER thermal supply only.
pub struct JkrAdhesiveTemperatureHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    pub transport: StepCoefficients<P>,
}

#[derive(Debug)]
pub enum JkrAdhesiveTemperatureHeatError<N, E> {
    Invalid(&'static str),
    Native(N),
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
impl<N, E> From<AdhesiveTemperatureError> for JkrAdhesiveTemperatureHeatError<N, E> {
    fn from(e: AdhesiveTemperatureError) -> Self {
        Self::Temperature(e)
    }
}
impl<N, E> From<JkrMixedModeError> for JkrAdhesiveTemperatureHeatError<N, E> {
    fn from(e: JkrMixedModeError) -> Self {
        Self::MixedMode(e)
    }
}
impl<N, E> From<DepositionError> for JkrAdhesiveTemperatureHeatError<N, E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<N, E> From<NewtonError> for JkrAdhesiveTemperatureHeatError<N, E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<N: core::fmt::Debug, E: core::fmt::Debug> core::fmt::Display
    for JkrAdhesiveTemperatureHeatError<N, E>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR temperature heat refused: {self:?}")
    }
}
impl<N: std::error::Error + 'static, E: std::error::Error + 'static> std::error::Error
    for JkrAdhesiveTemperatureHeatError<N, E>
{
}

pub struct JkrAdhesiveTemperatureHeatReport<R, P> {
    pub temperature: R,
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
pub struct JkrAdhesiveTemperatureHeatCandidate<'t, S, R, P> {
    source: SourceJkrMixedMode,
    bodies: [RigidBodyState; 2],
    state: S,
    thermal: DepositionSimulation<'t>,
    route: FrictionHeatRoute,
    report: JkrAdhesiveTemperatureHeatReport<R, P>,
}
impl<'t, S: Copy, R, P> JkrAdhesiveTemperatureHeatCandidate<'t, S, R, P> {
    pub const fn source(&self) -> &SourceJkrMixedMode {
        &self.source
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn state(&self) -> S {
        self.state
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn route(&self) -> &FrictionHeatRoute {
        &self.route
    }
    pub const fn report(&self) -> &JkrAdhesiveTemperatureHeatReport<R, P> {
        &self.report
    }
    /// Adopt all parts together. No mechanical impulse or geometry update occurs.
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrMixedMode,
        [RigidBodyState; 2],
        S,
        DepositionSimulation<'t>,
        FrictionHeatRoute,
        JkrAdhesiveTemperatureHeatReport<R, P>,
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

/// Captured native projection. No public callback can replace source history.
pub(super) struct NativeView<'a> {
    pub source: &'a SourceJkrMixedMode,
    pub normal: [f64; 3],
    pub reaction_arms: [Vec3; 2],
}
#[derive(Clone, Copy)]
pub(super) struct NativeEndpoint {
    pub time_s: f64,
    pub approach_m: f64,
    pub contact_area_m2: f64,
}
#[derive(Clone, Copy)]
pub(super) struct ParameterWork {
    pub storage_change_j: f64,
    pub parameter_work_j: f64,
    pub energy_residual_j: f64,
}
/// Supplied only by an original sealed native temperature continuation.
pub(super) struct NativeContinuation<S, R> {
    pub bodies: [RigidBodyState; 2],
    pub state: S,
    pub report: R,
    pub area_m2: f64,
    pub work: ParameterWork,
}
pub(super) trait NativeTemperature {
    type State: Copy;
    type Report: Copy;
    type Error: From<DynamicsError>;

    fn view(&self) -> NativeView<'_>;
    fn endpoint(previous: Self::State) -> NativeEndpoint;
    fn refresh(
        &self,
        previous: Self::State,
        refresh: &SourceJkrThermalRefresh,
        successor: &SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<NativeContinuation<Self::State, Self::Report>, Self::Error>;
}
struct RefreshedNative<S, R> {
    source: SourceJkrMixedMode,
    bodies: [RigidBodyState; 2],
    state: S,
    temperatures: AdhesiveTemperatures,
    report: R,
    work: ParameterWork,
}
struct Payload<S, R, P> {
    native: RefreshedNative<S, R>,
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
/// complete original contact parameter work retains its named driver port.
/// No frictional heat or mechanical entropy is inferred. Losing an attached
/// equilibrium, source support, receiver identity or any acceptance gate refuses
/// the entire candidate, including on final-trial cancellation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn propose<'t, T: NativeTemperature, P, E>(
    native: &T,
    previous: T::State,
    thermal: &DepositionSimulation<'t>,
    binding: JkrAdhesiveTemperatureHeatBinding<'_>,
    controls: JkrAdhesiveTemperatureHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        AdhesiveTemperatures,
        &SourceJkrMixedMode,
    ) -> Result<JkrAdhesiveTemperatureHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<
    JkrAdhesiveTemperatureHeatCandidate<'t, T::State, T::Report, P>,
    JkrAdhesiveTemperatureHeatError<T::Error, E>,
> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let endpoint = T::endpoint(previous);
    let c = controls.thermal;
    c.newton.validate()?;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(JkrAdhesiveTemperatureHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    if !endpoint.time_s.is_finite()
        || endpoint.time_s <= thermal.time_s()
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
        return Err(JkrAdhesiveTemperatureHeatError::Invalid(
            "bounded thermal endpoint and explicit parameter-work port required",
        ));
    }
    validate_receivers(native, previous, thermal, binding)?;
    let contact = native.view();
    let normal = contact.source.normal_source();
    let normal_state = normal
        .evaluate_attached(endpoint.approach_m)
        .map_err(AdhesiveTemperatureError::from)?
        .state();
    let step = FrictionHeatStep::new(thermal, endpoint.time_s, binding.route, check)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= c.maximum_trials {
            return Err(JkrAdhesiveTemperatureHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(JkrAdhesiveTemperatureHeatError::Invalid(
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
            let continuation = native
                .refresh(previous, &refresh, &source, gate)
                .map_err(JkrAdhesiveTemperatureHeatError::Native)?;
            let report = continuation.report;
            let bodies = continuation.bodies;
            validate_geometry(native, bodies, binding.route, controls.geometry_tolerance_m)?;
            let area = continuation.area_m2;
            let fractional =
                (area - binding.route.contact_area_m2()).abs() / binding.route.contact_area_m2();
            if !fractional.is_finite()
                || fractional > controls.maximum_fractional_area_change
                || continuation.work.parameter_work_j.abs()
                    > controls.maximum_absolute_parameter_work_j
            {
                return Err(JkrAdhesiveTemperatureHeatError::Invalid(
                    "source area or parameter work exceeds the explicit endpoint allowance",
                ));
            }
            let state = continuation.state;
            let work = continuation.work;
            Ok((
                area,
                RefreshedNative {
                    source,
                    bodies,
                    state,
                    temperatures,
                    report,
                    work,
                },
            ))
        };
        let transport = |endpoint: &DepositionSimulation<'t>,
                         _temperature,
                         staged: RefreshedNative<T::State, T::Report>| {
            let input = resolve(endpoint, staged.temperatures, &staged.source)
                .map_err(JkrAdhesiveTemperatureHeatError::Source)?;
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
        return Err(JkrAdhesiveTemperatureHeatError::Interrupted(
            result.progress,
        ));
    }
    if !result.report.converged {
        return Err(JkrAdhesiveTemperatureHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, c.audit_entropy)?;
    let temperature = &trial.payload().native.work;
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
        return Err(JkrAdhesiveTemperatureHeatError::EnergyBalance {
            residual_j: combined,
            absolute_stage_residual_j: absolute_stages,
        });
    }
    let (thermal, payload, heat, route) = trial.accept(c.heat, check)?;
    if thermal.time_s() != endpoint.time_s {
        return Err(JkrAdhesiveTemperatureHeatError::Invalid(
            "native and thermal endpoint clocks differ",
        ));
    }
    let staged = payload.native;
    check()?;
    Ok(JkrAdhesiveTemperatureHeatCandidate {
        source: staged.source,
        bodies: staged.bodies,
        state: staged.state,
        thermal,
        route,
        report: JkrAdhesiveTemperatureHeatReport {
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

fn numerical<N, E>(
    error: FiniteDifferenceError<JkrAdhesiveTemperatureHeatError<N, E>>,
) -> JkrAdhesiveTemperatureHeatError<N, E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(m) => JkrAdhesiveTemperatureHeatError::Invalid(m),
    }
}

fn validate_receivers<T: NativeTemperature, E>(
    native: &T,
    previous: T::State,
    thermal: &DepositionSimulation<'_>,
    binding: JkrAdhesiveTemperatureHeatBinding<'_>,
) -> Result<(), JkrAdhesiveTemperatureHeatError<T::Error, E>> {
    binding.route.validate(thermal)?;
    let elastic = native
        .view()
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
    let area = T::endpoint(previous).contact_area_m2;
    if binding.body_ids[0] == binding.body_ids[1]
        || (area - binding.route.contact_area_m2()).abs()
            > 64.0 * f64::EPSILON * area.max(binding.route.contact_area_m2())
    {
        return Err(JkrAdhesiveTemperatureHeatError::Invalid(
            "original ordered native bodies and endpoint area required",
        ));
    }
    for (i, id) in binding.route.cells().into_iter().enumerate() {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(JkrAdhesiveTemperatureHeatError::Invalid(
                "missing original thermal receiver",
            ))?;
        if cell.body_id != binding.body_ids[i]
            || binding.caloric[i].material_card_identity() != identities[i]
            || !thermal.uses_material_law(cell.material_id, binding.caloric[i])
        {
            return Err(JkrAdhesiveTemperatureHeatError::Invalid(
                "original ordered registered caloric materials required",
            ));
        }
        let caloric_temperature = thermal.temperature_k(id)?;
        let Some(source_temperature) = source_temperatures[i] else {
            return Err(JkrAdhesiveTemperatureHeatError::Invalid(
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
            return Err(JkrAdhesiveTemperatureHeatError::Invalid(
                "incoming bulk source temperatures differ from their accepted caloric cells",
            ));
        }
    }
    Ok(())
}

fn validate_geometry<T: NativeTemperature, E>(
    native: &T,
    bodies: [RigidBodyState; 2],
    route: &FrictionHeatRoute,
    tolerance: f64,
) -> Result<(), JkrAdhesiveTemperatureHeatError<T::Error, E>> {
    let g = route.geometry;
    let centre = [g[0].start_m, g[0].end_m]
        .into_iter()
        .find(|p| *p == g[1].start_m || *p == g[1].end_m)
        .ok_or(JkrAdhesiveTemperatureHeatError::Invalid(
            "original thermal face does not adjoin",
        ))?;
    let n = native.view().normal;
    let normal = Vec3::new(n[0], n[1], n[2]);
    let arms = native.view().reaction_arms;
    for i in 0..2 {
        let point = bodies[i]
            .pose()
            .point_world_from_body(arms[i])
            .map_err(|error| JkrAdhesiveTemperatureHeatError::Native(T::Error::from(error)))?;
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
            return Err(JkrAdhesiveTemperatureHeatError::Invalid(
                "native reaction plane and ordered thermal face must coincide",
            ));
        }
    }
    Ok(())
}
fn magnitude(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
