//! Material/flow-dependent closure around the native steady buoyancy solver.
//!
//! The callback re-evaluates original hydraulic and heat-transfer laws from
//! actual branch/segment mean fluid temperatures and signed flow. This driver
//! preserves topology, exchange area and reservoir temperatures. Its inner
//! owner closes density/gravity/flow/heat; the outer owner closes drag and hA.
//! Neither relaxation nor IQN-ILS supplies an acceptance criterion.
//!
//! Reference mass/Cp remain Boussinesq quantities. Material/correlation support
//! stays with the callback; this driver supplies no missing source properties,
//! film convention, turbulence model, coupled adjoint or stability claim.

use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransport;
use fs_qty::{Temperature, VolumetricFlowRate};

use super::{NaturalCirculationBranch, NaturalCirculationConfig, NaturalCirculationError,
    NaturalCirculationProblem, NaturalCirculationSolution, poll, positive, relative,
    solve_with_density_seed};
use super::super::transport::BranchThermalModel;
use crate::graph::LossGraph;

/// Actual fluid means and signed flow in the ORIGINAL branch/segment order.
/// Segment means are axial means of the exponential profile, not wall values.
/// A correlation's film or bulk property convention is an explicit callback
/// decision; the driver never silently evaluates every property at (Tw+Tb)/2.
#[derive(Debug, Clone, PartialEq)]
pub struct CirculationClosureState {
    pub branch_mean_temperatures_k: Vec<f64>,
    pub segment_mean_temperatures_k: Vec<Vec<f64>>,
    pub volume_flows_m3_s: Vec<f64>,
}

/// Rebuilt native coefficients and the caller's original material/correlation
/// owners. Payload is returned from the exact response used by the final solve.
#[derive(Debug)]
pub struct CirculationResponse<P> {
    pub graph: LossGraph,
    pub models: Vec<BranchThermalModel>,
    pub payload: P,
}

/// Independent physical residuals plus explicit work bounds. Inner budgets
/// are further capped by the remaining total budgets before each solve.
#[derive(Debug, Clone, Copy)]
pub struct CirculationClosureConfig {
    pub inner: NaturalCirculationConfig,
    pub max_iterations: usize,
    pub max_response_evaluations: usize,
    pub max_hydraulic_solves: usize,
    /// Only density-profile material evaluations; response callbacks retain
    /// their own property-query budgets and are counted separately above.
    pub max_buoyancy_material_evaluations: usize,
    pub temperature_tolerance_k: f64,
    pub absolute_flow_tolerance: VolumetricFlowRate,
    pub relative_flow_tolerance: f64,
    /// Symmetric relative change in EACH linear/quadratic drag and hA term.
    /// A mechanism appearing/disappearing has relative change one.
    pub relative_coefficient_tolerance: f64,
    pub relaxation: f64,
    pub temperature_scale_k: f64,
    pub flow_scale: VolumetricFlowRate,
    pub acceleration: Option<IqnIlsConfig>,
}

/// One accepted coupled field with its applied source response and residuals.
#[derive(Debug)]
pub struct CoupledCirculationSolution<P> {
    pub solution: NaturalCirculationSolution,
    pub response: CirculationResponse<P>,
    /// Exact input state at which the applied response was evaluated.
    pub response_state: CirculationClosureState,
    pub iterations: usize,
    pub response_evaluations: usize,
    pub hydraulic_solves: usize,
    pub buoyancy_material_evaluations: usize,
    pub max_temperature_residual_k: f64,
    pub max_flow_residual_m3_s: f64,
    pub max_relative_coefficient_residual: f64,
}

#[derive(Debug)]
pub enum CirculationClosureError<E> {
    Invalid(&'static str),
    Source(E),
    Natural(NaturalCirculationError),
    Budget(&'static str),
    NotConverged { iterations: usize, temperature_residual_k: f64,
        flow_residual_m3_s: f64, coefficient_residual: f64 },
}
impl<E: core::fmt::Debug> core::fmt::Display for CirculationClosureError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "circulation closure: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CirculationClosureError<E> {}
impl<E> From<NaturalCirculationError> for CirculationClosureError<E> {
    fn from(error: NaturalCirculationError) -> Self { Self::Natural(error) }
}

/// Close temperature-dependent drag and flow-dependent wall exchange with the
/// original native natural-circulation solve. Each response callback must return
/// the same topology/physical regions. Supplied seed flows are property guesses,
/// not imposed boundary flows. No partially converged field is published.
///
/// The output passes unrelaxed temperature/flow tests AND a fresh coefficient
/// evaluation at the actual solved field. The inner density solver also retains
/// its independent hydraulic confirmation. Warm-starting its integrated density
/// avoids throwing away the last physical density field between outer updates.
#[allow(clippy::too_many_arguments)]
pub fn solve_circulation_with_response<M, P, E>(
    cx: &Cx<'_>,
    problem: &NaturalCirculationProblem<'_>,
    config: CirculationClosureConfig,
    initial: &CirculationClosureState,
    mut material: impl FnMut(&Cx<'_>, usize, Temperature) -> Result<M, MaterialError>,
    mut response: impl FnMut(&Cx<'_>, &CirculationClosureState) -> Result<CirculationResponse<P>, E>,
) -> Result<CoupledCirculationSolution<P>, CirculationClosureError<E>>
where M: FluidTransport,
{
    poll(cx)?;
    admit(problem, config, initial)?;
    super::admit(cx, problem, config.inner, &initial.branch_mean_temperatures_k)?;
    let mut state = initial.clone();
    let dimension = encode(&state, config).len();
    let mut acceleration = config.acceleration
        .map(|policy| IqnIls::new(dimension, policy))
        .transpose().map_err(NaturalCirculationError::from)?;
    let mut warm_density: Option<Vec<f64>> = None;
    let mut responses = 0;
    let mut hydraulic_solves = 0;
    let mut material_evaluations = 0;
    let mut last_temperature = f64::INFINITY;
    let mut last_flow = f64::INFINITY;
    let mut last_coefficients = f64::INFINITY;
    for iteration in 0..config.max_iterations {
        poll(cx)?;
        let applied = evaluate_response(cx, problem, &state, &mut response,
            &mut responses, config.max_response_evaluations)?;
        let branches: Vec<_> = applied.models.iter().zip(problem.branches)
            .map(|(model, original)| NaturalCirculationBranch {
                model: model.clone(), axial_fractions: original.axial_fractions.clone(),
                wall_temperatures_k: original.wall_temperatures_k.clone(),
            }).collect();
        let current = NaturalCirculationProblem {
            graph: &applied.graph, reference: problem.reference,
            node_positions_m: problem.node_positions_m,
            acceleration_m_s2: problem.acceleration_m_s2,
            maximum_relative_density_departure: problem.maximum_relative_density_departure,
            pressure_datums: problem.pressure_datums, branches: &branches,
        };
        let mut inner = config.inner;
        inner.max_iterations = inner.max_iterations.min(
            (config.max_hydraulic_solves-hydraulic_solves)/2);
        inner.max_material_evaluations = inner.max_material_evaluations.min(
            config.max_buoyancy_material_evaluations-material_evaluations);
        if inner.max_iterations == 0 { return Err(CirculationClosureError::Budget("hydraulic solves")); }
        if inner.max_material_evaluations == 0 {
            return Err(CirculationClosureError::Budget("buoyancy material evaluations"));
        }
        let solved = solve_with_density_seed(cx, &current, inner,
            &state.branch_mean_temperatures_k, warm_density.as_deref(), &mut material)?;
        hydraulic_solves += solved.hydraulic_solves;
        material_evaluations += solved.material_evaluations;
        let next = solved_state(&solved);
        let (temperature_residual, flow_residual, flow_closed) = residual(&state, &next, config);
        last_temperature = temperature_residual;
        last_flow = flow_residual;
        // Re-evaluate the ORIGINAL source/correlation at the solved field.
        // A small state change cannot hide a discontinuous constitutive jump.
        let confirmed = evaluate_response(cx, problem, &next, &mut response,
            &mut responses, config.max_response_evaluations)?;
        last_coefficients = coefficient_residual(&applied, &confirmed);
        if !last_coefficients.is_finite() {
            return Err(CirculationClosureError::Invalid("nonfinite closure coefficients"));
        }
        if last_temperature <= config.temperature_tolerance_k && flow_closed
            && last_coefficients <= config.relative_coefficient_tolerance
        {
            poll(cx)?;
            return Ok(CoupledCirculationSolution {
                solution: solved, response: applied, response_state: state,
                iterations: iteration+1, response_evaluations: responses,
                hydraulic_solves, buoyancy_material_evaluations: material_evaluations,
                max_temperature_residual_k: last_temperature,
                max_flow_residual_m3_s: last_flow,
                max_relative_coefficient_residual: last_coefficients,
            });
        }
        warm_density = Some(solved.thermal_mean_densities_kg_m3);
        let reversed = state.volume_flows_m3_s.iter().zip(&next.volume_flows_m3_s)
            .any(|(a,b)| a.is_sign_negative() != b.is_sign_negative());
        if reversed {
            // Flow guesses choose source evaluation points, not a physical
            // interpolation path. Averaging opposite flows can invent zero
            // flow and reject a perfectly nonzero reversed circulation.
            if let Some(iqn) = acceleration.as_mut() { iqn.reset(); }
            state = next;
            continue;
        }
        let x = encode(&state, config);
        let image = encode(&next, config);
        let mut candidate: Vec<_> = x.iter().zip(&image)
            .map(|(a,b)| a + config.relaxation*(b-a)).collect();
        if let Some(iqn) = acceleration.as_mut() {
            let step = iqn.step(&x, &image, config.relaxation)
                .map_err(NaturalCirculationError::from)?;
            let extrapolated = decode(&state, &step.values, config);
            if valid_state(problem, &extrapolated)
                && extrapolated.volume_flows_m3_s.iter().zip(&next.volume_flows_m3_s)
                    .all(|(a,b)| a.is_sign_negative() == b.is_sign_negative())
            { candidate = step.values; } else { iqn.reset(); }
        }
        let candidate = decode(&state, &candidate, config);
        if !valid_state(problem, &candidate) {
            return Err(CirculationClosureError::Invalid("nonrepresentable relaxed closure state"));
        }
        state = candidate;
    }
    Err(CirculationClosureError::NotConverged {
        iterations: config.max_iterations, temperature_residual_k: last_temperature,
        flow_residual_m3_s: last_flow, coefficient_residual: last_coefficients,
    })
}

fn evaluate_response<P,E>(
    cx: &Cx<'_>, problem: &NaturalCirculationProblem<'_>, state: &CirculationClosureState,
    callback: &mut impl FnMut(&Cx<'_>, &CirculationClosureState) -> Result<CirculationResponse<P>,E>,
    count: &mut usize, budget: usize,
) -> Result<CirculationResponse<P>, CirculationClosureError<E>> {
    poll(cx)?;
    if *count >= budget { return Err(CirculationClosureError::Budget("material responses")); }
    *count += 1;
    let result = callback(cx, state).map_err(CirculationClosureError::Source)?;
    poll(cx)?;
    if result.graph.node_count() != problem.graph.node_count()
        || result.graph.branches().len() != problem.branches.len()
        || result.models.len() != problem.branches.len()
        || result.graph.branches().iter().zip(problem.graph.branches())
            .any(|(a,b)| a.from != b.from || a.to != b.to || a.loss.name != b.loss.name)
    { return Err(CirculationClosureError::Invalid("response changed hydraulic topology or branch identity")); }
    for (model, original) in result.models.iter().zip(problem.branches) {
        poll(cx)?;
        match (model, &original.model) {
            (BranchThermalModel::Adiabatic, BranchThermalModel::Adiabatic) => {}
            (BranchThermalModel::Exchange(a), BranchThermalModel::Exchange(b))
                if a.len() == b.len() && a.iter().zip(b).all(|(a,b)|
                    a.region() == b.region() && a.area_m2() == b.area_m2()) => {}
            _ => return Err(CirculationClosureError::Invalid("response changed exchange region, area or ordering")),
        }
    }
    Ok(result)
}

fn solved_state(solution: &NaturalCirculationSolution) -> CirculationClosureState {
    let segment_mean_temperatures_k = solution.transport.branches.iter()
        .zip(&solution.hydraulics.branches).map(|(heat, flow)| {
            let mut means: Vec<_> = heat.march.as_ref().map(|m|
                m.segments.iter().map(|s| s.reference_temperature_k).collect()).unwrap_or_default();
            if flow.flow.value().is_sign_negative() { means.reverse(); }
            means
        }).collect();
    CirculationClosureState {
        branch_mean_temperatures_k: solution.branch_mean_temperatures_k.clone(),
        segment_mean_temperatures_k,
        volume_flows_m3_s: solution.hydraulics.branches.iter().map(|b| b.flow.value()).collect(),
    }
}

fn relative_change(a: f64, b: f64) -> f64 {
    if !a.is_finite() || !b.is_finite() || a < 0.0 || b < 0.0 { return f64::NAN; }
    let scale = a.max(b);
    if scale == 0.0 { 0.0 } else { (a/scale-b/scale).abs() }
}

fn coefficient_residual<P>(a: &CirculationResponse<P>, b: &CirculationResponse<P>) -> f64 {
    let mut maximum = 0.0_f64;
    for i in 0..a.graph.branches().len() {
        let pairs = [
            (a.graph.branches()[i].loss.resistance.value(), b.graph.branches()[i].loss.resistance.value()),
            (a.graph.linear_losses()[i].as_ref().map_or(0.0, |l| l.resistance.value()),
             b.graph.linear_losses()[i].as_ref().map_or(0.0, |l| l.resistance.value())),
        ];
        for (x,y) in pairs {
            let value = relative_change(x,y);
            if !value.is_finite() { return f64::NAN; }
            maximum = maximum.max(value);
        }
        if let (BranchThermalModel::Exchange(x), BranchThermalModel::Exchange(y)) =
            (&a.models[i], &b.models[i])
        {
            for (x,y) in x.iter().zip(y) {
                // Areas have already been checked identical, so comparing h
                // is exactly the relative hA change without product overflow.
                let value = relative_change(x.htc_w_per_m2_k(), y.htc_w_per_m2_k());
                if !value.is_finite() { return f64::NAN; }
                maximum = maximum.max(value);
            }
        }
    }
    maximum
}

fn residual(a: &CirculationClosureState, b: &CirculationClosureState,
    config: CirculationClosureConfig) -> (f64, f64, bool) {
    let mut temperature = 0.0_f64;
    for (x,y) in a.branch_mean_temperatures_k.iter().chain(a.segment_mean_temperatures_k.iter().flatten())
        .zip(b.branch_mean_temperatures_k.iter().chain(b.segment_mean_temperatures_k.iter().flatten()))
    { temperature = temperature.max((x-y).abs()); }
    let mut flow = 0.0_f64;
    let mut closed = true;
    for (x,y) in a.volume_flows_m3_s.iter().zip(&b.volume_flows_m3_s) {
        let delta = (x-y).abs();
        flow = flow.max(delta);
        let threshold = config.absolute_flow_tolerance.value()
            + config.relative_flow_tolerance*x.abs().max(y.abs());
        closed &= threshold.is_finite() && delta <= threshold
            && x.is_sign_negative() == y.is_sign_negative();
    }
    (temperature, flow, closed)
}

fn encode(state: &CirculationClosureState, config: CirculationClosureConfig) -> Vec<f64> {
    state.branch_mean_temperatures_k.iter().chain(state.segment_mean_temperatures_k.iter().flatten())
        .map(|value| value/config.temperature_scale_k)
        .chain(state.volume_flows_m3_s.iter().map(|value| value/config.flow_scale.value())).collect()
}
fn decode(shape: &CirculationClosureState, values: &[f64],
    config: CirculationClosureConfig) -> CirculationClosureState {
    let mut values = values.iter();
    CirculationClosureState {
        branch_mean_temperatures_k: shape.branch_mean_temperatures_k.iter().map(|_|
            values.next().copied().unwrap_or(f64::NAN)*config.temperature_scale_k).collect(),
        segment_mean_temperatures_k: shape.segment_mean_temperatures_k.iter().map(|segments|
            segments.iter().map(|_| values.next().copied().unwrap_or(f64::NAN)
                *config.temperature_scale_k).collect()).collect(),
        volume_flows_m3_s: shape.volume_flows_m3_s.iter().map(|_|
            values.next().copied().unwrap_or(f64::NAN)*config.flow_scale.value()).collect(),
    }
}

fn valid_state(problem: &NaturalCirculationProblem<'_>, state: &CirculationClosureState) -> bool {
    let n = problem.branches.len();
    state.branch_mean_temperatures_k.len() == n && state.segment_mean_temperatures_k.len() == n
        && state.volume_flows_m3_s.len() == n
        && state.branch_mean_temperatures_k.iter().all(|t| positive(*t))
        && state.segment_mean_temperatures_k.iter().zip(problem.branches).all(|(t,b)| {
            let count = match &b.model { BranchThermalModel::Adiabatic => 0,
                BranchThermalModel::Exchange(s) => s.len() };
            t.len() == count && t.iter().all(|t| positive(*t))
        })
        && state.volume_flows_m3_s.iter().all(|q| q.is_finite() && *q != 0.0)
}
fn admit<E>(problem: &NaturalCirculationProblem<'_>, config: CirculationClosureConfig,
    state: &CirculationClosureState) -> Result<(), CirculationClosureError<E>> {
    if !valid_state(problem, state) || problem.branches.is_empty()
        || config.max_iterations == 0 || config.max_response_evaluations < 2
        || config.max_hydraulic_solves < 2 || config.max_buoyancy_material_evaluations == 0
        || !positive(config.temperature_tolerance_k) || !positive(config.absolute_flow_tolerance.value())
        || !relative(config.relative_flow_tolerance) || !relative(config.relative_coefficient_tolerance)
        || !positive(config.relaxation) || config.relaxation > 1.0
        || !positive(config.temperature_scale_k) || !positive(config.flow_scale.value())
        || encode(state, config).iter().any(|v| !v.is_finite())
    { return Err(CirculationClosureError::Invalid("invalid coupled closure state, residual scales or budgets")); }
    Ok(())
}

#[cfg(test)]
mod tests;
