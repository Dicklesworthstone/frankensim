//! Temperature-dependent buoyancy on a closed native hydraulic graph.
//!
//! Each branch is geometrically straight; its thermal segments occupy explicit
//! axial fractions. The existing exponential channel law supplies T(s), and a
//! caller's material owner supplies rho(T). Integrating rho along each branch
//! supplies its gravitational head. Mixing, hydraulic loss and heat transport
//! remain owned by the existing graph and closed-circulation solvers.
//!
//! This is a steady Boussinesq model: one reference rho and Cp conserve volume
//! and sensible heat, while only gravity uses the temperature-dependent density.
//! Drag and hA remain frozen. No compressibility, stagnant-fluid conduction,
//! inertial stability, coupled adjoint or experimental validation is implied.

use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig, IqnIlsError};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::det;
use fs_qty::{Density, Temperature, VolumetricFlowRate};

use crate::graph::{FixedPressure, GraphError, GraphSolution, GraphSolveConfig, LossGraph, Sum};
use super::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportError,
    TransportMarch, TransportNetwork, MAX_CIRCULATION_NODES};

/// Thermal geometry in the branch's DECLARED from-to order. Axial fractions
/// describe length, not wetted-area weights; each segment has uniform h and
/// wetted perimeter along that length. Split a bent branch into graph edges.
#[derive(Debug, Clone, PartialEq)]
pub struct NaturalCirculationBranch {
    pub model: BranchThermalModel,
    /// One positive fraction per exchange segment, summing to one. Empty for
    /// an adiabatic branch, whose flowing fluid has uniform temperature.
    pub axial_fractions: Vec<f64>,
    /// Prescribed wall reservoirs in the SAME declared segment order.
    pub wall_temperatures_k: Vec<f64>,
}

/// Immutable source/topology inputs. The material callback retains ownership
/// of its source data; this checked reference snapshot adds no data authority.
#[derive(Debug)]
pub struct NaturalCirculationProblem<'a> {
    pub graph: &'a LossGraph,
    pub reference: FluidTransportState,
    pub node_positions_m: &'a [[f64; 3]],
    pub acceleration_m_s2: [f64; 3],
    /// Bound in [0,1) on every sampled abs(rho-rho0)/rho0.
    pub maximum_relative_density_departure: f64,
    /// Pressure datums only. The closed thermal owner independently checks
    /// that no datum supplies or withdraws fluid.
    pub pressure_datums: &'a [FixedPressure],
    pub branches: &'a [NaturalCirculationBranch],
}

/// All nonlinear work and stopping controls are explicit. IQN-ILS is optional
/// acceleration; only unrelaxed material/head closure and a fresh hydraulic
/// solve can accept a result.
#[derive(Debug, Clone, Copy)]
pub struct NaturalCirculationConfig {
    pub hydraulics: GraphSolveConfig,
    pub transport: TransportConfig,
    pub max_iterations: usize,
    pub max_material_evaluations: usize,
    /// Maximum composite density quadrature panels per segment, in 2..=1024.
    /// Panels must also resolve the thermal decay: NTU/panels <= 1.
    pub max_density_quadrature_panels: usize,
    /// Positive absolute and nonnegative relative density residual allowances.
    pub absolute_density_tolerance_kg_m3: f64,
    pub relative_density_tolerance: f64,
    pub absolute_flow_change_tolerance: VolumetricFlowRate,
    pub relative_flow_change_tolerance: f64,
    /// Fallback density relaxation in (0,1]. Never used as a stopping test.
    pub relaxation: f64,
    pub acceleration: Option<IqnIlsConfig>,
    /// The sensible-heat model omits viscous heating only within this explicit
    /// total watt allowance. No hydraulic dissipation is silently added twice.
    pub maximum_neglected_dissipation_w: f64,
}

/// A mutually consistent nominal buoyant-flow/thermal state within its stated
/// residuals. Physical branch/source labels retain their declared ordering.
#[derive(Debug, Clone, PartialEq)]
pub struct NaturalCirculationSolution {
    pub hydraulics: GraphSolution,
    pub transport: TransportMarch,
    /// Exact axial means of the exponential temperature profiles, K.
    pub branch_mean_temperatures_k: Vec<f64>,
    /// Numerically integrated source density at the accepted temperature field.
    /// The gravity binding retains the applied densities, which can differ by
    /// at most the reported nonlinear residual.
    pub thermal_mean_densities_kg_m3: Vec<f64>,
    pub iterations: usize,
    pub hydraulic_solves: usize,
    pub material_evaluations: usize,
    pub max_density_residual_kg_m3: f64,
    pub max_gravity_head_residual_pa: f64,
    /// Difference from re-solving hydraulics at the unrelaxed thermal density.
    pub max_confirmed_flow_change_m3_s: f64,
    pub neglected_dissipation_w: f64,
    /// Sum of segment-weighted paired quadrature differences over all branches;
    /// this is a numerical estimator, not a rigorous source-error enclosure.
    pub density_quadrature_error_kg_m3: f64,
}

#[derive(Debug)]
pub enum NaturalCirculationError {
    InvalidInput(&'static str),
    Branch { branch: usize, reason: &'static str },
    Graph(GraphError),
    Transport(TransportError),
    Material(MaterialError),
    Acceleration(IqnIlsError),
    MaterialBudget { evaluations: usize },
    Quadrature { branch: usize, segment: usize },
    Dissipation { watts: f64, allowance_w: f64 },
    NotConverged { iterations: usize, max_density_residual_kg_m3: f64,
        max_confirmed_flow_change_m3_s: f64 },
    Interrupted,
}
impl core::fmt::Display for NaturalCirculationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "natural circulation: {self:?}")
    }
}
impl std::error::Error for NaturalCirculationError {}
impl From<GraphError> for NaturalCirculationError {
    fn from(error: GraphError) -> Self { Self::Graph(error) }
}
impl From<TransportError> for NaturalCirculationError {
    fn from(error: TransportError) -> Self { Self::Transport(error) }
}
impl From<MaterialError> for NaturalCirculationError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<IqnIlsError> for NaturalCirculationError {
    fn from(error: IqnIlsError) -> Self { Self::Acceleration(error) }
}

/// Solve closed-network buoyancy feedback from explicit initial branch
/// temperatures. Seeds initialize gravity; they are never thermal reservoirs.
/// The source is queried at actual fluid temperatures and must retain its own
/// pressure, composition, phase and validity rules. Source errors are terminal.
///
/// Density uses resolved composite four-point Gauss quadrature with a paired
/// two-point estimate. Thermal-boundary-layer resolution, panel and total source
/// budgets prevent silent under-resolution. Sampling does not certify arbitrary
/// material functions; source uncertainty remains the caller's authority.
///
/// Zero-flow heated or thermally undetermined adiabatic branches refuse. This
/// selects an admitted steady circulation from the seed, not its dynamical
/// stability or uniqueness among multiple nonlinear steady states.
pub fn solve_natural_circulation<M: FluidTransport>(
    cx: &Cx<'_>, problem: &NaturalCirculationProblem<'_>,
    config: NaturalCirculationConfig, initial_branch_temperatures_k: &[f64],
    material: impl FnMut(&Cx<'_>, usize, Temperature) -> Result<M, MaterialError>,
) -> Result<NaturalCirculationSolution, NaturalCirculationError> {
    poll(cx)?;
    admit(cx, problem, config, initial_branch_temperatures_k)?;
    let rho0 = problem.reference.density_kg_m3();
    let density_tolerance = finite(config.absolute_density_tolerance_kg_m3
        + config.relative_density_tolerance * rho0)?;
    let mut accelerator = config.acceleration.map(|policy|
        IqnIls::new(problem.branches.len(), policy)).transpose()?;
    let mut sampler = Sampler { material, reference: problem.reference,
        maximum_departure: problem.maximum_relative_density_departure,
        evaluations: 0, limit: config.max_material_evaluations };
    let mut densities = Vec::with_capacity(problem.branches.len());
    for (branch, &temperature) in initial_branch_temperatures_k.iter().enumerate() {
        densities.push(sampler.at(cx, branch, temperature)?);
    }
    let air = TransportAir { density: Density::new(rho0),
        specific_heat_j_kg_k: problem.reference.specific_heat_j_kg_k() };
    let mut previous_directions = Vec::new();
    let mut hydraulic_solves = 0;
    let mut last_density_residual = f64::INFINITY;
    let mut last_flow_change = f64::INFINITY;
    for iteration in 0..config.max_iterations {
        poll(cx)?;
        let hydraulics = solve_hydraulics(cx, problem, config.hydraulics, &densities)?;
        hydraulic_solves += 1;
        let directions: Vec<_> = hydraulics.branches.iter()
            .map(|branch| branch.flow.value().is_sign_negative()).collect();
        if directions != previous_directions {
            if let Some(accelerator) = &mut accelerator { accelerator.reset(); }
            previous_directions.clone_from(&directions);
        }
        let mut walls = Vec::new();
        for (branch, &reverse) in problem.branches.iter().zip(&directions) {
            if reverse { walls.extend(branch.wall_temperatures_k.iter().rev().copied()); }
            else { walls.extend_from_slice(&branch.wall_temperatures_k); }
        }
        let network = TransportNetwork::new_circulating(cx, &hydraulics, air,
            problem.branches.iter().map(|branch| branch.model.clone()).collect(),
            Temperature::new(initial_branch_temperatures_k[0]), config.transport)?;
        let transport = network.march(cx, &walls)?;
        let dissipation = dissipated_power(&hydraulics)?;
        if dissipation > config.maximum_neglected_dissipation_w {
            return Err(NaturalCirculationError::Dissipation {
                watts: dissipation, allowance_w: config.maximum_neglected_dissipation_w });
        }
        let mut thermal_density = Vec::with_capacity(densities.len());
        let mut mean_temperatures = Vec::with_capacity(densities.len());
        let mut quadrature_error = Sum::default();
        let mut max_head_residual = 0.0_f64;
        last_density_residual = 0.0;
        for (edge, branch) in problem.branches.iter().enumerate() {
            poll(cx)?;
            let result = &transport.branches[edge];
            let mut mean = Sum::default();
            let mut departure = Sum::default();
            if let Some(march) = &result.march {
                for (slot, segment) in march.segments.iter().enumerate() {
                    let authored = if directions[edge] { march.segments.len()-1-slot } else { slot };
                    let fraction = branch.axial_fractions[authored];
                    let wall = branch.wall_temperatures_k[authored];
                    let (rho, error) = sampler.integrate(cx, edge, authored,
                        segment.inlet_temperature_k, wall, segment.ntu,
                        config.max_density_quadrature_panels, density_tolerance * 0.125)?;
                    mean.add(fraction * segment.reference_temperature_k)?;
                    departure.add(fraction * (rho-rho0))?;
                    quadrature_error.add(fraction * error)?;
                }
            } else {
                let temperature = result.inlet_temperature_k.ok_or(
                    NaturalCirculationError::Branch { branch: edge,
                        reason: "stagnant adiabatic branch has no determined thermal density" })?;
                mean.add(temperature)?;
                departure.add(sampler.at(cx, edge, temperature)? - rho0)?;
            }
            mean_temperatures.push(mean.value()?);
            let rho = finite(rho0 + departure.value()?)?;
            let residual = (rho-densities[edge]).abs();
            last_density_residual = last_density_residual.max(residual);
            let topology = &problem.graph.branches()[edge];
            let mut acceleration_distance = Sum::default();
            for axis in 0..3 {
                acceleration_distance.add(problem.acceleration_m_s2[axis]
                    * (problem.node_positions_m[topology.to][axis]
                        - problem.node_positions_m[topology.from][axis]))?;
            }
            max_head_residual = max_head_residual.max(finite(
                residual * acceleration_distance.value()?.abs())?);
            thermal_density.push(rho);
        }
        if last_density_residual <= density_tolerance {
            // Confirm the UNRELAXED map with a fresh physical hydraulic solve.
            // A tiny relaxation step cannot conceal a large flow discrepancy.
            let confirmed = solve_hydraulics(cx, problem, config.hydraulics, &thermal_density)?;
            hydraulic_solves += 1;
            let mut admitted = true;
            last_flow_change = 0.0;
            for (old, new) in hydraulics.branches.iter().zip(&confirmed.branches) {
                poll(cx)?;
                let (q, next) = (old.flow.value(), new.flow.value());
                let change = finite((next-q).abs())?;
                last_flow_change = last_flow_change.max(change);
                let tolerance = finite(config.absolute_flow_change_tolerance.value()
                    + config.relative_flow_change_tolerance * q.abs().max(next.abs()))?;
                admitted &= change <= tolerance && next != 0.0
                    && q.is_sign_negative() == next.is_sign_negative();
            }
            if admitted {
                poll(cx)?;
                return Ok(NaturalCirculationSolution { hydraulics, transport,
                    branch_mean_temperatures_k: mean_temperatures,
                    thermal_mean_densities_kg_m3: thermal_density,
                    iterations: iteration+1, hydraulic_solves,
                    material_evaluations: sampler.evaluations,
                    max_density_residual_kg_m3: last_density_residual,
                    max_gravity_head_residual_pa: max_head_residual,
                    max_confirmed_flow_change_m3_s: last_flow_change,
                    neglected_dissipation_w: dissipation,
                    density_quadrature_error_kg_m3: quadrature_error.value()? });
            }
        }
        let fallback: Vec<_> = densities.iter().zip(&thermal_density)
            .map(|(&old, &new)| old + config.relaxation * (new-old)).collect();
        densities = if let Some(accelerator) = &mut accelerator {
            // Scale every coordinate by the same reference density.
            let x: Vec<_> = densities.iter().map(|rho| (rho-rho0)/rho0).collect();
            let image: Vec<_> = thermal_density.iter().map(|rho| (rho-rho0)/rho0).collect();
            let proposal = accelerator.step(&x, &image, config.relaxation)?;
            if proposal.values.iter().all(|value| value.is_finite()
                && value.abs() <= problem.maximum_relative_density_departure) {
                proposal.values.into_iter().map(|value| rho0 * (1.0+value)).collect()
            } else {
                accelerator.reset();
                fallback
            }
        } else { fallback };
    }
    Err(NaturalCirculationError::NotConverged { iterations: config.max_iterations,
        max_density_residual_kg_m3: last_density_residual,
        max_confirmed_flow_change_m3_s: last_flow_change })
}

fn admit(cx: &Cx<'_>, problem: &NaturalCirculationProblem<'_>, config: NaturalCirculationConfig,
    temperatures: &[f64]) -> Result<(), NaturalCirculationError> {
    if problem.graph.node_count() > MAX_CIRCULATION_NODES
        || problem.branches.is_empty() || problem.branches.len() != problem.graph.branches().len()
        || temperatures.len() != problem.branches.len()
        || temperatures.iter().any(|value| !value.is_finite() || *value <= 0.0)
        || problem.node_positions_m.len() != problem.graph.node_count()
        || problem.node_positions_m.iter().flatten().any(|value| !value.is_finite())
        || problem.acceleration_m_s2.iter().any(|value| !value.is_finite())
        || !relative(problem.maximum_relative_density_departure)
        || config.max_iterations == 0 || config.max_material_evaluations == 0
        || !(2..=1024).contains(&config.max_density_quadrature_panels)
        || !positive(config.absolute_density_tolerance_kg_m3)
        || !relative(config.relative_density_tolerance)
        || !positive(config.absolute_flow_change_tolerance.value())
        || !relative(config.relative_flow_change_tolerance)
        || !positive(config.relaxation) || config.relaxation > 1.0
        || !config.maximum_neglected_dissipation_w.is_finite()
        || config.maximum_neglected_dissipation_w < 0.0
    { return Err(NaturalCirculationError::InvalidInput("invalid circulation geometry, seed or work/residual budget")); }
    // Admit the original hydraulic controls before invoking a source callback.
    config.hydraulics.validate()?;
    for (edge, branch) in problem.branches.iter().enumerate() {
        poll(cx)?;
        let n = match &branch.model { BranchThermalModel::Adiabatic => 0,
            BranchThermalModel::Exchange(rows) => {
                if rows.is_empty() { return Err(NaturalCirculationError::Branch {
                    branch: edge, reason: "empty exchange branch" }); }
                rows.len()
            }
        };
        if branch.axial_fractions.len() != n || branch.wall_temperatures_k.len() != n
            || branch.axial_fractions.iter().any(|value| !positive(*value))
            || branch.wall_temperatures_k.iter().any(|value| !positive(*value))
            || (n > 0 && (branch.axial_fractions.iter().sum::<f64>()-1.0).abs() > 64.0*f64::EPSILON)
        { return Err(NaturalCirculationError::Branch { branch: edge,
            reason: "exchange segments need positive axial fractions summing to one and matching positive wall temperatures" }); }
        let topology = &problem.graph.branches()[edge];
        if problem.node_positions_m[topology.from] == problem.node_positions_m[topology.to] {
            return Err(NaturalCirculationError::Branch { branch: edge,
                reason: "straight branch endpoints must have distinct physical positions" });
        }
    }
    Ok(())
}

fn solve_hydraulics(cx: &Cx<'_>, problem: &NaturalCirculationProblem<'_>,
    config: GraphSolveConfig, densities: &[f64]) -> Result<GraphSolution, NaturalCirculationError> {
    let reference = problem.reference;
    let states: Result<Vec<_>, _> = densities.iter().map(|&rho| FluidTransportState::try_new(
        rho, reference.specific_heat_j_kg_k(), reference.thermal_conductivity_w_m_k(),
        reference.dynamic_viscosity_pa_s())).collect();
    Ok(problem.graph.clone().with_boussinesq_gravity(&reference, &states?,
        problem.node_positions_m, problem.acceleration_m_s2,
        problem.maximum_relative_density_departure, cx)?
        .solve(problem.pressure_datums, config, cx)?)
}

fn dissipated_power(flow: &GraphSolution) -> Result<f64, NaturalCirculationError> {
    let mut sum = Sum::default();
    for branch in &flow.branches {
        sum.add(branch.flow.value() * branch.dissipative_pressure_drop.value())?;
    }
    Ok(sum.value()?)
}

struct Sampler<F> {
    material: F, reference: FluidTransportState, maximum_departure: f64,
    evaluations: usize, limit: usize,
}
impl<M: FluidTransport, F: FnMut(&Cx<'_>, usize, Temperature) -> Result<M, MaterialError>> Sampler<F> {
    fn at(&mut self, cx: &Cx<'_>, branch: usize, temperature: f64)
        -> Result<f64, NaturalCirculationError> {
        poll(cx)?;
        if self.evaluations == self.limit {
            return Err(NaturalCirculationError::MaterialBudget { evaluations: self.evaluations });
        }
        if !positive(temperature) { return Err(NaturalCirculationError::InvalidInput("nonpositive sampled fluid temperature")); }
        self.evaluations += 1;
        let rho = (self.material)(cx, branch, Temperature::new(temperature))?
            .transport_state()?.density_kg_m3();
        poll(cx)?;
        if ((rho-self.reference.density_kg_m3())/self.reference.density_kg_m3()).abs()
            > self.maximum_departure {
            return Err(NaturalCirculationError::Branch { branch,
                reason: "sampled density exceeds the declared Boussinesq allowance" });
        }
        Ok(rho)
    }

    #[allow(clippy::too_many_arguments)]
    fn integrate(&mut self, cx: &Cx<'_>, branch: usize, segment: usize,
        inlet: f64, wall: f64, ntu: f64, limit: usize, tolerance: f64)
        -> Result<(f64, f64), NaturalCirculationError> {
        // Gauss sites are interior. The physically present endpoint states
        // must independently satisfy the source domain and density allowance.
        let inlet_density = self.at(cx, branch, inlet)?;
        let outlet = finite(inlet + (wall-inlet)*(-det::expm1(-ntu)))?;
        if outlet == inlet { return Ok((inlet_density, 0.0)); }
        self.at(cx, branch, outlet)?;
        let mut panels = 1;
        while ntu/(panels as f64) > 1.0 {
            panels *= 2;
            if panels > limit { return Err(NaturalCirculationError::Quadrature { branch, segment }); }
        }
        loop {
            let mut fourth = Sum::default();
            let mut error = Sum::default();
            let rho0 = self.reference.density_kg_m3();
            for panel in 0..panels {
                let width = 1.0/panels as f64;
                let midpoint = (panel as f64+0.5)*width;
                let mut coarse = Sum::default();
                let mut fine = Sum::default();
                for (point, weight, high) in [
                    (-0.5773502691896257, 1.0, false), (0.5773502691896257, 1.0, false),
                    (-0.8611363115940526, 0.34785484513745385, true),
                    (-0.3399810435848563, 0.6521451548625461, true),
                    (0.3399810435848563, 0.6521451548625461, true),
                    (0.8611363115940526, 0.34785484513745385, true),
                ] {
                    let x = midpoint + 0.5*width*point;
                    let temperature = inlet + (wall-inlet)*(-det::expm1(-ntu*x));
                    let value = weight*(self.at(cx, branch, temperature)?-rho0);
                    if high { fine.add(value)?; } else { coarse.add(value)?; }
                }
                let integral = 0.5*width*fine.value()?;
                fourth.add(integral)?;
                error.add((integral-0.5*width*coarse.value()?).abs())?;
            }
            let error = error.value()?;
            if error <= tolerance { return Ok((finite(rho0+fourth.value()?)?, error)); }
            if panels > limit/2 { return Err(NaturalCirculationError::Quadrature { branch, segment }); }
            panels *= 2;
        }
    }
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn relative(value: f64) -> bool { value.is_finite() && (0.0..1.0).contains(&value) }
fn finite(value: f64) -> Result<f64, NaturalCirculationError> {
    if value.is_finite() { Ok(value) }
    else { Err(NaturalCirculationError::InvalidInput("unrepresentable circulation arithmetic")) }
}
fn poll(cx: &Cx<'_>) -> Result<(), NaturalCirculationError> {
    cx.checkpoint().map_err(|_| NaturalCirculationError::Interrupted)
}

#[cfg(test)]
mod tests;
