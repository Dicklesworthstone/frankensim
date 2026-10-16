//! Implicit mixed-air feedback for a checked total-enthalpy endpoint.
//! The direct h cotangent is independent of temperature seeds, including where
//! dT/dh=0. The solid's complete radiative tangent stays inside every sweep;
//! only original convection references and wall temperatures connect to air.

use super::*;
use fs_conduction::transient::enthalpy::adjoint::{EnthalpyRobinResponse, EnthalpyStepGradient};

/// Total endpoint gradient after the true transpose interface equation closes.
#[derive(Debug, Clone)]
pub struct CoupledEnthalpyGradient {
    /// Full physical multiplier and previous-h / nodal-source pullbacks.
    /// `iterations` inside this field counts only the final solid solve.
    pub solid: EnthalpyStepGradient,
    /// External supply-temperature derivatives; zero off supply nodes.
    pub inlets: Vec<f64>,
    /// Original convective ln(HTC) derivatives, including air-side hA effects.
    pub log_htc: Vec<f64>,
    /// Assembled nodal watt-load derivatives, equal to dt times solid.adjoint.
    pub nodal_load: Vec<f64>,
    /// Interface multiplier actually used by the returned solid and air VJPs.
    pub interface_adjoint: Vec<f64>,
    /// Maximum unrelaxed transpose-interface residual.
    pub interface_residual: f64,
    /// Completed solid/air derivative sweeps.
    pub iterations: usize,
    /// Sum of actual solid Krylov columns over every interface sweep.
    pub solid_krylov_iterations: usize,
}

/// Total enthalpy endpoint pullback extended by a common signed-flow scale.
#[derive(Debug, Clone)]
pub struct CoupledEnthalpyFlowScaleGradient {
    /// Original controls, including the direct latent-history cotangent.
    pub thermal: CoupledEnthalpyGradient,
    /// Objective units per ln(s) when every branch and external heat-capacity
    /// rate scales by s at fixed HTC, areas, topology and inlet temperatures.
    /// This is not an individual branch-flow or fan operating-point derivative.
    pub log_flow_scale: f64,
}

/// Concrete checked enthalpy/transport binding, with no user callback escape.
/// Geometry, masses, material laws, time step and old h are frozen. Flows stay
/// fixed except for the explicit common-flow-scale pullbacks.
/// Its source/history gradients include mixed-air and any radiative feedback.
pub struct CoupledEnthalpyLinearization<'a, 'step, 'm, 'flow> {
    solid: &'a EnthalpyRobinResponse<'step, 'm>,
    air: TransportLinearization<'a, 'flow>,
}

impl<'a, 'step, 'm, 'flow> CoupledEnthalpyLinearization<'a, 'step, 'm, 'flow> {
    /// Recheck names/order, exact convection HTC/area, effective references and
    /// every branch's convective heat balance at the accepted fixed point.
    /// Radiation is neither added to air heat nor used as air conductance.
    pub fn new(
        cx: &Cx<'_>,
        network: &'a TransportNetwork<'flow>,
        solid: &'a EnthalpyRobinResponse<'step, 'm>,
        primal_gate: &ConjugateConfig,
    ) -> Result<Self> {
        poll(cx)?;
        let walls = solid.wall_means(cx, solid.temperature())?;
        let air = bind_solid_transport(
            cx,
            network,
            solid.ports(),
            solid.robin_fluxes(),
            &walls,
            primal_gate,
        )?;
        Ok(Self { solid, air })
    }

    /// Correctly sized temperature/wall/convection/air objective. Pass direct
    /// h carry separately to pullback; this deliberately does not hide a T->h inverse.
    #[must_use]
    pub fn zero_objective(&self) -> CoupledObjective {
        CoupledObjective {
            nodal_temperatures: vec![0.0; self.solid.temperature().len()],
            wall_temperatures: vec![0.0; self.solid.ports().len()],
            solid_heat_rates: vec![0.0; self.solid.ports().len()],
            air: self.air.zero_objective(),
        }
    }

    /// Close the endpoint transpose by stationary relaxation. h_carry has one
    /// entry per mesh vertex and is added after constitutive T'(h) scaling.
    /// The declared linear budget applies per solid sweep; the interface budget
    /// bounds sweeps. Failed or cancelled work returns no partial gradient.
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        h_carry: &[f64],
        config: InterfaceSolveConfig,
    ) -> Result<CoupledEnthalpyGradient> {
        self.pullback_driver(cx, objective, h_carry, config, None)
    }

    /// Same physical transpose using fresh bounded IQN-ILS history. Acceptance
    /// uses the unrelaxed equation residual, never update size or acceleration history.
    pub fn pullback_iqn(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        h_carry: &[f64],
        config: InterfaceSolveConfig,
        acceleration: IqnIlsConfig,
    ) -> Result<CoupledEnthalpyGradient> {
        self.pullback_driver(cx, objective, h_carry, config, Some(acceleration))
    }

    /// Extend the stationary coupled pullback by uniform scaling of all signed
    /// flow/capacity rates, with their proportions, directions and HTC fixed.
    /// The accepted interface multiplier includes the direct h carry and the
    /// full solid/radiation response. One additional air reverse sweep computes
    /// the scale derivative; no further solid Krylov work is required.
    ///
    /// # Errors
    /// Preserves all endpoint, transpose-residual, budget and cancellation
    /// refusals. A failed or nonfinite contraction publishes no gradient.
    pub fn pullback_flow_scale(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        h_carry: &[f64],
        config: InterfaceSolveConfig,
    ) -> Result<CoupledEnthalpyFlowScaleGradient> {
        let thermal = self.pullback(cx, objective, h_carry, config)?;
        self.with_flow_scale(cx, objective, thermal)
    }

    /// Same uniform-flow control using fresh bounded IQN-ILS for the coupled
    /// transpose. Acceptance checks the true interface equation before the
    /// air-only contraction, never a frozen-wall or frozen-history derivative.
    ///
    /// # Errors
    /// The refusals of [`Self::pullback_iqn`] and nonfinite/cancelled contraction.
    pub fn pullback_flow_scale_iqn(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        h_carry: &[f64],
        config: InterfaceSolveConfig,
        acceleration: IqnIlsConfig,
    ) -> Result<CoupledEnthalpyFlowScaleGradient> {
        let thermal = self.pullback_iqn(cx, objective, h_carry, config, acceleration)?;
        self.with_flow_scale(cx, objective, thermal)
    }

    fn with_flow_scale(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        thermal: CoupledEnthalpyGradient,
    ) -> Result<CoupledEnthalpyFlowScaleGradient> {
        let log_flow_scale = super::flow_scale::contract_flow_scale(
            cx,
            &self.air,
            &objective.air,
            &thermal.interface_adjoint,
        )?;
        Ok(CoupledEnthalpyFlowScaleGradient {
            thermal,
            log_flow_scale,
        })
    }

    fn pullback_driver(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        h_carry: &[f64],
        config: InterfaceSolveConfig,
        acceleration: Option<IqnIlsConfig>,
    ) -> Result<CoupledEnthalpyGradient> {
        validate(config)?;
        poll(cx)?;
        let mut accelerator = acceleration
            .map(|policy| IqnIls::new(self.solid.ports().len(), policy))
            .transpose()
            .map_err(CoupledSensitivityError::Acceleration)?;
        // Validate every air objective slot before augmenting its references.
        self.air.pullback(cx, &objective.air)?;
        if objective.wall_temperatures.len() != self.solid.ports().len() {
            return Err(bad("one wall-objective weight per port required"));
        }
        let mut current = vec![0.0; self.solid.ports().len()];
        let mut solid_krylov_iterations = 0_usize;
        for iteration in 1..=config.max_iterations {
            poll(cx)?;
            let mut weights = objective.air.clone();
            for (weight, value) in weights.references.iter_mut().zip(&current) {
                *weight = finite(*weight + value)?;
            }
            let air = self.air.pullback(cx, &weights)?;
            let walls = objective
                .wall_temperatures
                .iter()
                .zip(&air.walls)
                .map(|(a, b)| finite(a + b))
                .collect::<Result<Vec<_>>>()?;
            let solid = self.solid.pullback(
                cx,
                h_carry,
                &objective.nodal_temperatures,
                &walls,
                &objective.solid_heat_rates,
            )?;
            solid_krylov_iterations = solid_krylov_iterations
                .checked_add(solid.transport.iterations)
                .ok_or_else(|| bad("enthalpy derivative work counter overflow"))?;
            let (residual, tolerance) = equation(&current, &solid.references, config)?;
            poll(cx)?;
            if residual <= tolerance {
                let log_htc = solid
                    .log_htc
                    .iter()
                    .zip(&air.log_conductances)
                    .map(|(a, b)| finite(a + b))
                    .collect::<Result<Vec<_>>>()?;
                return Ok(CoupledEnthalpyGradient {
                    solid: solid.transport,
                    inlets: air.inlets,
                    log_htc,
                    nodal_load: solid.nodal_load,
                    interface_adjoint: current,
                    interface_residual: residual,
                    iterations: iteration,
                    solid_krylov_iterations,
                });
            }
            if iteration == config.max_iterations {
                return Err(CoupledSensitivityError::DidNotConverge {
                    iterations: iteration,
                    residual,
                    tolerance,
                });
            }
            update(
                cx,
                &mut current,
                &solid.references,
                config.relaxation,
                &mut accelerator,
            )?;
        }
        unreachable!("positive bounded iteration returns")
    }
}
