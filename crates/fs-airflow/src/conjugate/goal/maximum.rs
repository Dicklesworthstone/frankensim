//! Cached maximum-error bounds with the actual ordered air-reference feedback.
//! Flow, material coefficients and geometry are fixed, but wall temperatures
//! and every downstream air reference are variables of the coupled system.
//!
//! Lowering uses the production segment's admitted effectiveness and NTU ratio.
//! The affine coefficient assembly is rounded binary64: ensuing enclosures
//! concern this STORED model, not its assembly error or exact transcendental
//! coefficients, nonlinear conductivity, radiation, hydraulics or validation.

use fs_conduction::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalAnalyzer, LinearRobinFeedbackAnalyzer,
    RobinFeedbackAnalysisConfig,
};
use fs_conduction::{ConductionProblem, LinearConfig, ThermalInterfaces};
use fs_exec::Cx;

use super::{AirPath, BTreeSet, Result, admitted_exchange_terms, bad, finite, poll};

/// Immutable affine reference law in branch-major, stream-wise region order.
/// Created only from admitted physical AirPath values, never from unit probes
/// or a finite-difference Jacobian at one temperature field.
#[derive(Debug, Clone, PartialEq)]
pub struct AirReferenceLaw {
    regions: Vec<String>,
    offset_k: Vec<f64>,
    wall_matrix: Vec<f64>,
}
impl AirReferenceLaw {
    /// Exact coordinate identities in order.
    #[must_use]
    pub fn regions(&self) -> &[String] { &self.regions }
    /// Affine offsets in kelvin; separate inlets are never mixed.
    #[must_use]
    pub fn offset_k(&self) -> &[f64] { &self.offset_k }
    /// Row-major square map from wall means to Robin references.
    #[must_use]
    pub fn wall_matrix(&self) -> &[f64] { &self.wall_matrix }
    /// Evaluate the rounded affine model for diagnostics or independent replay.
    /// This is not an interval enclosure of the original march's arithmetic.
    ///
    /// # Errors
    /// Wrong wall count, nonfinite arithmetic or caller cancellation.
    pub fn evaluate(&self, cx: &Cx<'_>, walls_k: &[f64]) -> Result<Vec<f64>> {
        poll(cx)?;
        let n = self.regions.len();
        if walls_k.len() != n { return Err(bad("one wall temperature per air-reference coordinate required")); }
        for &value in walls_k { poll(cx)?; finite(value)?; }
        let mut out = zeros(cx, n)?;
        for i in 0..n {
            let mut value = self.offset_k[i];
            for (j, &wall) in walls_k.iter().enumerate() {
                poll(cx)?;
                value = finite(self.wall_matrix[i*n+j].mul_add(wall, value))?;
            }
            out[i] = value;
        }
        poll(cx)?;
        Ok(out)
    }
}

fn zeros(cx: &Cx<'_>, count: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| bad("air-reference coefficient allocation refused"))?;
    for i in 0..count { if i % 512 == 0 { poll(cx)?; } values.push(0.0); }
    Ok(values)
}

/// Lower all branches using the production exponential segment law. A segment
/// depends on itself and all upstream walls on its OWN branch, never another
/// branch or a downstream wall. Input order, not sorting, defines coordinates.
///
/// At segment j, r_j = g_j T_in,j + (1-g_j) T_wall,j and
/// T_out,j = (1-eps_j) T_in,j + eps_j T_wall,j. Carry the complete affine
/// inlet to the next segment. No primal run or differentiated iteration trace
/// is needed. Geometry/flow and their coefficients remain fixed.
///
/// # Errors
/// Empty paths, repeated region ownership across branches, count/work overflow,
/// explicit port/coefficient cap, nonfinite coefficients and cancellation.
pub fn affine_reference_law(
    cx: &Cx<'_>, paths: &[AirPath], max_ports: usize, max_coefficients: usize,
) -> Result<AirReferenceLaw> {
    poll(cx)?;
    if paths.is_empty() { return Err(bad("at least one air path is required")); }
    let mut n = 0_usize;
    for path in paths {
        poll(cx)?;
        n = n.checked_add(path.segments().len()).ok_or_else(|| bad("air port count overflow"))?;
        if n > max_ports { return Err(bad("air-reference port limit exceeded")); }
    }
    let square = n.checked_mul(n).ok_or_else(|| bad("air-reference coefficient count overflow"))?;
    if square > max_coefficients { return Err(bad("air-reference coefficient limit exceeded")); }
    let mut seen = BTreeSet::new();
    let mut regions = Vec::new();
    regions.try_reserve_exact(n).map_err(|_| bad("air-reference name allocation refused"))?;
    for path in paths {
        for segment in path.segments() {
            poll(cx)?;
            if !seen.insert(segment.region()) { return Err(bad("an air region may belong to only one branch")); }
            regions.push(segment.region().to_string());
        }
    }
    let mut offset_k = zeros(cx, n)?;
    let mut wall_matrix = zeros(cx, square)?;
    let mut start = 0;
    for path in paths {
        let mut constant = path.inlet_temperature_k();
        let mut upstream = zeros(cx, path.segments().len())?;
        for (j, segment) in path.segments().iter().enumerate() {
            poll(cx)?;
            let (_, eps, g) = admitted_exchange_terms(segment, path.capacity_rate_w_per_k())?;
            let row = start+j;
            offset_k[row] = finite(g * constant)?;
            for k in 0..j {
                poll(cx)?;
                wall_matrix[row*n+start+k] = finite(g * upstream[k])?;
            }
            wall_matrix[row*n+row] = finite(1.0-g)?;
            let carry = finite(1.0-eps)?;
            constant = finite(carry * constant)?;
            for value in &mut upstream[..j] { poll(cx)?; *value = finite(carry * *value)?; }
            upstream[j] = eps;
        }
        start += path.segments().len();
    }
    poll(cx)?;
    Ok(AirReferenceLaw { regions, offset_k, wall_matrix })
}

/// Prepare a real linear solid/air maximum analyzer, reusable across iterates.
/// The full air-reference map is retained, not the references of the supplied
/// field. Production material, fixed contact, Robin faces and prescribed-node
/// lifts are consumed by the existing conduction analyzer. Per-port h and area
/// must match the actual AirPath, as in the existing coupled-goal comparison.
///
/// Structural admission and optional response work are explicit; response
/// iterations are shared across all branches, not a new budget for each port.
/// Failure to establish contraction leaves a typed no-bound analysis. It does
/// not trigger a frozen-solid fallback. No nonlinear or radiation admission is
/// implied by this API.
///
/// # Errors
/// Existing linear-solid and field admission, duplicate/missing/mismatched air
/// ports, exceeded preparation budgets, nonfinite arithmetic and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn prepare_linear_maximum<'m>(
    cx: &Cx<'_>, problem: ConductionProblem<'m>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], linear: LinearConfig, temperature: &[f64],
    solid_config: LinearGoalAnalysisConfig, feedback_config: RobinFeedbackAnalysisConfig,
) -> Result<LinearRobinFeedbackAnalyzer<'m>> {
    let law = affine_reference_law(cx, paths, feedback_config.residual.max_ports,
        feedback_config.max_lowering_entries)?;
    let regions: Vec<&str> = law.regions.iter().map(String::as_str).collect();
    let solid = LinearGoalAnalyzer::new_for_maximum(
        cx, problem, interfaces, linear, temperature, solid_config,
    )?;
    let analyzer = solid.with_robin_feedback(cx, &regions, &law.offset_k, &law.wall_matrix, feedback_config)?;
    for (port, segment) in analyzer.ports().iter().zip(paths.iter().flat_map(|path| path.segments())) {
        poll(cx)?;
        if port.name != segment.region() || port.htc_w_m2_k != segment.htc_w_per_m2_k()
            || (port.area_m2-segment.area_m2()).abs()
                > 128.0 * f64::EPSILON * port.area_m2.max(segment.area_m2())
        { return Err(bad("air maximum analysis requires matching solid/air names, h and wetted area")); }
    }
    poll(cx)?;
    Ok(analyzer)
}
