//! Implicit shared-solid/air-network derivatives, not frozen-wall gradients.
//!
//! With solid response T=S(r,s,f), wall means W=C*T, and air references
//! r=R(W,b,s), solve (I-R_W*C*S_r) dr = R_W*C*(S_s ds+S_f df)+R_b db+R_s ds.
//! Reverse mode solves the transpose interface equation before accumulating
//! every inlet, conductance and assembled-load gradient. The two existing
//! producer linearizations do all FEM and transport differentiation.
//!
//! Both stationary and explicit vector IQN-ILS methods check the UNRELAXED
//! equation residual using fresh producer evaluations. Acceleration solves the
//! same affine interface map; it does not differentiate the primal iteration
//! history or assemble a dense interface Jacobian. Either method may exhaust
//! its budget; neither claims interval accuracy or passivity.
//! Hydraulics, mesh, material laws and prescribed solid temperatures stay
//! fixed. A log(h) control changes BOTH the Robin matrix and air-side hA.

use std::fmt;
use fs_conduction::ConductionError;
use fs_conduction::adjoint::robin::{RobinDifferential, RobinLinearization};
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig, IqnIlsError};
use fs_exec::Cx;

use crate::conjugate::{ConjugateConfig, Relaxation, SolidRegionState, solve_conjugate_from};
use super::super::transport::{TransportError, TransportNetwork};
use super::super::transport::sensitivity::{TransportDifferential, TransportLinearization, TransportObjective};

mod flow_scale;
pub use flow_scale::CoupledFlowScaleGradient;

/// Explicit budget for the reduced interface linear equation.
#[derive(Debug, Clone, Copy)]
pub struct InterfaceSolveConfig {
    /// Maximum solid/air derivative sweeps, including the successful sweep.
    pub max_iterations: usize,
    /// Positive absolute equation-residual floor, in the unknown's units.
    /// Tangents use kelvin; adjoints use objective units per kelvin.
    pub absolute_tolerance: f64,
    /// Relative residual tolerance in [0,1), against max(|current|,|proposal|).
    pub relative_tolerance: f64,
    /// Fixed relaxation in (0,1], also the IQN startup/rank-zero fallback.
    /// Convergence is tested before relaxation or vector acceleration.
    pub relaxation: f64,
}

/// External controls; wall temperatures and Robin references are NOT controls.
#[derive(Debug, Clone)]
pub struct CoupledDirection {
    /// Supply-temperature differences per hydraulic node, zero off supplies.
    pub inlets_k: Vec<f64>,
    /// ln(h) differences in network region order, at fixed mesh/area.
    pub log_htc: Vec<f64>,
    /// Assembled solid nodal-load differences, W, in full mesh order.
    pub nodal_load_w: Vec<f64>,
}

/// Weights on solid and air outputs. All explicit heat-rate signs are outward
/// from the solid / absorbed by the air, as in the underlying producer APIs.
#[derive(Debug, Clone)]
pub struct CoupledObjective {
    /// Weights on full nodal solid temperatures.
    pub nodal_temperatures: Vec<f64>,
    /// Weights on area-mean solid temperatures in network region order.
    pub wall_temperatures: Vec<f64>,
    /// Weights on selected outward solid heat rates in network region order.
    pub solid_heat_rates: Vec<f64>,
    /// Weights on the transported outputs, including effective references.
    pub air: TransportObjective,
}

/// Tangent of the mutually coupled model, with the interface residual retained.
#[derive(Debug, Clone)]
pub struct CoupledDifferential {
    /// Solid derivative evaluated at solved_reference_tangent_k.
    pub solid: RobinDifferential,
    /// Air derivative evaluated from that solid's wall-temperature derivative.
    pub air: TransportDifferential,
    /// Reference tangent actually used by the solid; the air proposal is separate.
    pub solved_reference_tangent_k: Vec<f64>,
    /// Maximum unrelaxed interface-equation residual.
    pub interface_residual: f64,
    /// Completed derivative sweeps.
    pub iterations: usize,
}

/// Total gradient after closing the transpose interface equation.
#[derive(Debug, Clone)]
pub struct CoupledGradient {
    /// Derivatives with respect to external supply temperatures; zero off supplies.
    pub inlets: Vec<f64>,
    /// Derivatives with respect to ln(h), including solid and air contributions.
    pub log_htc: Vec<f64>,
    /// Derivatives with respect to assembled solid nodal loads.
    pub nodal_load: Vec<f64>,
    /// Interface multiplier actually used for the returned pullbacks.
    pub interface_adjoint: Vec<f64>,
    /// Maximum unrelaxed transpose-equation residual.
    pub interface_residual: f64,
    /// Completed derivative sweeps.
    pub iterations: usize,
}

/// Typed producer, binding, interruption or implicit-solve refusal.
#[derive(Debug)]
pub enum CoupledSensitivityError {
    /// A selected producer, vector or numerical budget is incompatible.
    InvalidInput(&'static str),
    /// Underlying solid refusal, preserved without coercing it to convergence.
    Solid(ConductionError),
    /// Underlying transport or conjugate admission refusal.
    Air(TransportError),
    /// Invalid acceleration policy or unrepresentable vector-update arithmetic.
    Acceleration(IqnIlsError),
    /// A context checkpoint refused more work; no partial derivative is returned.
    Interrupted,
    /// The interface equation did not meet its budget.
    DidNotConverge { iterations: usize, residual: f64, tolerance: f64 },
}
impl fmt::Display for CoupledSensitivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "coupled sensitivity: {self:?}") }
}
impl std::error::Error for CoupledSensitivityError {}
impl From<ConductionError> for CoupledSensitivityError {
    fn from(error: ConductionError) -> Self { Self::Solid(error) }
}
impl From<TransportError> for CoupledSensitivityError {
    fn from(error: TransportError) -> Self { Self::Air(error) }
}
type Result<T> = std::result::Result<T, CoupledSensitivityError>;

/// Concrete FEM/air binding at a checked coupled fixed point. No user-supplied
/// derivative callback or publicly assembled primal result can replace a side.
pub struct CoupledLinearization<'a, 'flow> {
    solid: &'a RobinLinearization,
    air: TransportLinearization<'a, 'flow>,
}
impl<'a, 'flow> CoupledLinearization<'a, 'flow> {
    /// Bind the same names, order, area and h on both sides, and recheck primal
    /// temperature and per-branch watt balance with the existing conjugate gate.
    /// The solid must have been solved against its actual final references,
    /// e.g. those returned by solve_coupled_transport, not the next proposal.
    pub fn new(cx: &Cx<'_>, network: &'a TransportNetwork<'flow>, solid: &'a RobinLinearization,
        primal_gate: &ConjugateConfig) -> Result<Self> {
        poll(cx)?;
        let names = network.regions();
        if names.is_empty() || names.len() != solid.ports().len()
            || names.iter().zip(solid.ports()).any(|(a,b)| *a != b.name.as_str())
        { return Err(bad("solid ports must equal all network regions in exact order")); }
        let walls = solid.wall_means(cx, &solid.primal().temperature)?;
        let air = network.linearize(cx, &walls)?;
        let probe = ConjugateConfig { max_iterations: 1, relaxation: Relaxation::Fixed { omega: 1.0 }, ..*primal_gate };
        for branch in 0..network.hydraulics().branches.len() {
            poll(cx)?;
            let range = network.row_range(branch);
            if range.is_empty() { continue; }
            let inlet = air.primal().branches[branch].inlet_temperature_k
                .ok_or_else(|| bad("missing active branch inlet"))?;
            let path = network.path(branch, inlet)?;
            let mut states = Vec::new();
            let mut refs = Vec::new();
            for (index, segment) in range.clone().zip(path.segments()) {
                poll(cx)?;
                let port = &solid.ports()[index];
                if port.htc_w_m2_k != segment.htc_w_per_m2_k()
                    || (port.area_m2 - segment.area_m2()).abs() > 128.0*f64::EPSILON*port.area_m2.max(segment.area_m2())
                { return Err(bad("solid and air must share h and area, not merely a region name")); }
                let flux = solid.primal().report.robin_fluxes.iter().find(|flux| flux.region == port.name)
                    .ok_or_else(|| bad("solid report is missing the bound Robin region"))?;
                states.push(SolidRegionState::from_robin_flux(flux));
                refs.push(port.reference_k);
                if (air.primal().reference_temperatures_k[index]-port.reference_k).abs() > primal_gate.temperature_tolerance_k {
                    return Err(bad("solid references are not a converged coupled fixed point"));
                }
            }
            solve_conjugate_from(cx, &path, &probe, &refs, |_, _| Ok(states.clone()))
                .map_err(TransportError::from)?;
        }
        poll(cx)?;
        Ok(Self { solid, air })
    }

    /// Correctly sized zero external perturbation.
    #[must_use]
    pub fn zero_direction(&self) -> CoupledDirection {
        let d = self.air.zero_direction();
        CoupledDirection { inlets_k: d.inlets_k, log_htc: d.log_conductances,
            nodal_load_w: vec![0.0; self.solid.primal().temperature.len()] }
    }

    /// Correctly sized zero solid/air objective.
    #[must_use]
    pub fn zero_objective(&self) -> CoupledObjective {
        CoupledObjective { nodal_temperatures: vec![0.0; self.solid.primal().temperature.len()],
            wall_temperatures: vec![0.0; self.solid.ports().len()], solid_heat_rates: vec![0.0; self.solid.ports().len()],
            air: self.air.zero_objective() }
    }

    /// Solve the forward implicit equation by stationary relaxation with true,
    /// unrelaxed residual checks. No partially converged field is returned.
    pub fn apply(&self, cx: &Cx<'_>, direction: &CoupledDirection, config: InterfaceSolveConfig) -> Result<CoupledDifferential> {
        self.apply_driver(cx, direction, config, None)
    }

    /// Solve the same forward equation using bounded vector IQN-ILS. All
    /// secants belong to this tangent solve at the fixed admitted primal;
    /// they are not recycled from primal or previous derivative iterations.
    /// The true interface equation, not update size, decides convergence.
    ///
    /// # Errors
    /// Producer, input, acceleration and interruption refusals propagate;
    /// exhausting the declared sweep budget returns `DidNotConverge`.
    pub fn apply_iqn(&self, cx: &Cx<'_>, direction: &CoupledDirection,
        config: InterfaceSolveConfig, acceleration: IqnIlsConfig) -> Result<CoupledDifferential> {
        self.apply_driver(cx, direction, config, Some(acceleration))
    }

    fn apply_driver(&self, cx: &Cx<'_>, direction: &CoupledDirection,
        config: InterfaceSolveConfig, acceleration: Option<IqnIlsConfig>) -> Result<CoupledDifferential> {
        validate(config)?;
        poll(cx)?;
        let mut accelerator = acceleration.map(|policy| IqnIls::new(self.solid.ports().len(), policy))
            .transpose().map_err(CoupledSensitivityError::Acceleration)?;
        let mut current = vec![0.0; self.solid.ports().len()];
        for iteration in 1..=config.max_iterations {
            poll(cx)?;
            let mut ds = self.solid.zero_direction();
            ds.references_k.clone_from(&current);
            ds.log_htc.clone_from(&direction.log_htc);
            ds.nodal_load_w.clone_from(&direction.nodal_load_w);
            let solid = self.solid.apply(cx, &ds)?;
            let mut da = self.air.zero_direction();
            da.walls_k.clone_from(&solid.mean_wall_temperatures_k);
            da.inlets_k.clone_from(&direction.inlets_k);
            da.log_conductances.clone_from(&direction.log_htc);
            let air = self.air.apply(cx, &da)?;
            let (residual, tolerance) = equation(&current, &air.reference_temperatures_k, config)?;
            poll(cx)?;
            if residual <= tolerance {
                return Ok(CoupledDifferential { solid, air, solved_reference_tangent_k: current,
                    interface_residual: residual, iterations: iteration });
            }
            if iteration == config.max_iterations {
                return Err(CoupledSensitivityError::DidNotConverge { iterations: iteration, residual, tolerance });
            }
            update(cx, &mut current, &air.reference_temperatures_k, config.relaxation, &mut accelerator)?;
        }
        unreachable!("positive bounded iteration returns")
    }

    /// Solve the transpose interface equation by stationary relaxation, then
    /// accumulate all controls, including direct heat-functional terms and
    /// solid feedback. Never differentiate the primal iteration history.
    pub fn pullback(&self, cx: &Cx<'_>, objective: &CoupledObjective, config: InterfaceSolveConfig) -> Result<CoupledGradient> {
        self.pullback_driver(cx, objective, config, None)
    }

    /// Solve the same transpose equation using its own bounded vector IQN-ILS
    /// history. Direct solid/air objective terms and all coupled controls are
    /// accumulated at the accepted adjoint, not at an unchecked extrapolation.
    ///
    /// # Errors
    /// As for [`Self::apply_iqn`]; no partial gradient is published on failure.
    pub fn pullback_iqn(&self, cx: &Cx<'_>, objective: &CoupledObjective,
        config: InterfaceSolveConfig, acceleration: IqnIlsConfig) -> Result<CoupledGradient> {
        self.pullback_driver(cx, objective, config, Some(acceleration))
    }

    fn pullback_driver(&self, cx: &Cx<'_>, objective: &CoupledObjective,
        config: InterfaceSolveConfig, acceleration: Option<IqnIlsConfig>) -> Result<CoupledGradient> {
        validate(config)?;
        poll(cx)?;
        let mut accelerator = acceleration.map(|policy| IqnIls::new(self.solid.ports().len(), policy))
            .transpose().map_err(CoupledSensitivityError::Acceleration)?;
        // Validate every air objective slot before augmenting its references.
        self.air.pullback(cx, &objective.air)?;
        if objective.wall_temperatures.len() != self.solid.ports().len() {
            return Err(bad("one wall-objective weight per port required"));
        }
        let mut current = vec![0.0; self.solid.ports().len()];
        for iteration in 1..=config.max_iterations {
            poll(cx)?;
            let mut weights = objective.air.clone();
            for (weight, value) in weights.references.iter_mut().zip(&current) { *weight = finite(*weight + value)?; }
            let air = self.air.pullback(cx, &weights)?;
            let walls = objective.wall_temperatures.iter().zip(&air.walls)
                .map(|(a,b)| finite(a+b)).collect::<Result<Vec<_>>>()?;
            let solid = self.solid.pullback(cx, &objective.nodal_temperatures, &walls, &objective.solid_heat_rates)?;
            let (residual, tolerance) = equation(&current, &solid.references, config)?;
            poll(cx)?;
            if residual <= tolerance {
                let log_htc = solid.log_htc.iter().zip(&air.log_conductances)
                    .map(|(a,b)| finite(a+b)).collect::<Result<Vec<_>>>()?;
                return Ok(CoupledGradient { inlets: air.inlets, log_htc, nodal_load: solid.nodal_load,
                    interface_adjoint: current, interface_residual: residual, iterations: iteration });
            }
            if iteration == config.max_iterations {
                return Err(CoupledSensitivityError::DidNotConverge { iterations: iteration, residual, tolerance });
            }
            update(cx, &mut current, &solid.references, config.relaxation, &mut accelerator)?;
        }
        unreachable!("positive bounded iteration returns")
    }
}
fn update(cx: &Cx<'_>, current: &mut Vec<f64>, next: &[f64], omega: f64,
    accelerator: &mut Option<IqnIls>) -> Result<()> {
    poll(cx)?;
    if let Some(accelerator) = accelerator.as_mut() {
        let proposal = accelerator.step(current, next, omega)
            .map_err(CoupledSensitivityError::Acceleration)?;
        poll(cx)?;
        *current = proposal.values;
    } else {
        relax(current, next, omega)?;
    }
    Ok(())
}
fn validate(config: InterfaceSolveConfig) -> Result<()> {
    if config.max_iterations == 0 || !config.absolute_tolerance.is_finite() || config.absolute_tolerance <= 0.0
        || !(0.0..1.0).contains(&config.relative_tolerance)
        || !config.relaxation.is_finite() || config.relaxation <= 0.0 || config.relaxation > 1.0
    { return Err(bad("invalid implicit interface budget, tolerance or relaxation")); }
    Ok(())
}
fn equation(current: &[f64], next: &[f64], config: InterfaceSolveConfig) -> Result<(f64, f64)> {
    if current.len() != next.len() { return Err(bad("interface vector arity mismatch")); }
    let mut residual = 0.0_f64;
    let mut scale = 0.0_f64;
    for (&a,&b) in current.iter().zip(next) {
        residual = residual.max(finite(b-a)?.abs());
        scale = scale.max(finite(a)?.abs()).max(finite(b)?.abs());
    }
    Ok((residual, finite(config.absolute_tolerance + config.relative_tolerance*scale)?))
}
fn relax(current: &mut [f64], next: &[f64], omega: f64) -> Result<()> {
    for (a,b) in current.iter_mut().zip(next) { *a = finite(*a + omega*(b-*a))?; }
    Ok(())
}
fn bad(reason: &'static str) -> CoupledSensitivityError { CoupledSensitivityError::InvalidInput(reason) }
fn finite(value: f64) -> Result<f64> { if value.is_finite() { Ok(value) } else { Err(bad("nonfinite implicit derivative arithmetic")) } }
fn poll(cx: &Cx<'_>) -> Result<()> { cx.checkpoint().map_err(|_| CoupledSensitivityError::Interrupted) }

#[cfg(test)]
mod tests;