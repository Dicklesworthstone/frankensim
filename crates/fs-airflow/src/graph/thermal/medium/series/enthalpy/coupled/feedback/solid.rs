//! Joint solid / enthalpy / temperature-dependent convection fixed point.
//!
//! The iteration unknown is the wall-temperature vector, not independently
//! tunable h and coolant references. Each trial material march produces both
//! boundary coefficients and references; the solid response then determines
//! the next wall state. Original material/correlation rows remain retained.

use fs_convection::CorrelationRobinBoundary;
use fs_couple::AitkenRelaxation;
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig};
use fs_exec::Cx;
use fs_qty::{Temperature, VolumetricFlowRate};

use crate::AirflowError;
use crate::conjugate::{ConjugateConfig, Relaxation};
use crate::graph::thermal::{admit, checked, iqn_error};
use crate::graph::thermal::medium::{FlowHeatSection, FrozenFlowMedium};
use super::{LiquidFeedbackBudget, LiquidFeedbackError, LiquidHeatMarch, LiquidHeatPath, checkpoint};
use super::super::LiquidSolidResponse;

impl From<AirflowError> for LiquidFeedbackError {
    fn from(error: AirflowError) -> Self { Self::Enthalpy(error.into()) }
}

/// Joint result after every local and global coupling gate passed.
#[derive(Debug)]
pub struct CoupledLiquidFeedbackSolution<'a, R> {
    /// Original solid result; its accepted state is never advanced in trial calls.
    pub solid: LiquidSolidResponse<R>,
    /// EXACT material march whose boundary rows were sent to the solid callback.
    pub boundary_state: LiquidHeatMarch<'a>,
    /// Re-evaluated material march at the returned solid's actual wall temperatures.
    pub fluid: LiquidHeatMarch<'a>,
    /// Number of solid callbacks performed.
    pub iterations: usize,
    /// Maximum unrelaxed wall-temperature residual, K.
    pub maximum_wall_residual_k: f64,
    /// Maximum effective coolant-reference residual, K.
    pub maximum_reference_residual_k: f64,
    /// Maximum relative coefficient change across the last solid response.
    pub maximum_coefficient_residual_relative: f64,
    /// Maximum absolute local solid/coolant heat imbalance, W.
    pub maximum_heat_imbalance_w: f64,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Solve the solid and temperature-dependent coolant laws together.
    ///
    /// Uses existing IQN-ILS and configured fixed/Aitken fallback. Every trial
    /// invokes the same material, caloric and correlation owners. The callback
    /// MUST solve from the same accepted solid history and return the actual
    /// sent reference for each region. Both wall and coolant-reference changes,
    /// relative h change, and EACH region's watt imbalance must pass. Opposing
    /// local errors cannot cancel. A small temperature step alone never accepts
    /// an unclosed power balance, and iterations continue while progress remains.
    ///
    /// An IQN extrapolation outside the caloric wall span discards its secants
    /// and retries the configured relaxation; actual source/solver refusals are
    /// terminal. No constitutive fallback, pressure update or flow redistribution
    /// is implied. This is a regional heat-transfer model, not a boundary-layer
    /// or boiling solver. Inner and outer numerical tolerances remain explicit.
    ///
    /// # Errors
    /// Preserves owner refusals; rejects malformed responses, unsupported states,
    /// unclosed balances, cancellation and iteration exhaustion without a result.
    pub fn solve_coupled_iqn<R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], config: &ConjugateConfig,
        local: LiquidFeedbackBudget, acceleration: IqnIlsConfig, solid: F,
    ) -> Result<CoupledLiquidFeedbackSolution<'a, R>, LiquidFeedbackError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        let (thermal, _, _) = self.solve_coupled_iqn_march(
            cx, initial_walls_k, config, local, acceleration, solid,
            |cx, walls| self.march(cx, walls, local).map(|thermal| (thermal, ())),
        )?;
        Ok(thermal)
    }

    // One wall iteration, response admission and regional power gate for both
    // ordinary and dissipative liquid trajectories. Retain each march's native
    // payload beside its thermal state, including the exact rows sent to the
    // solid and the independent re-evaluation at the returned wall temperatures.
    pub(super) fn solve_coupled_iqn_march<R, P, E, F, M>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], config: &ConjugateConfig,
        local: LiquidFeedbackBudget, acceleration: IqnIlsConfig, mut solid: F,
        mut march: M,
    ) -> Result<(CoupledLiquidFeedbackSolution<'a, R>, P, P), E>
    where
        E: From<LiquidFeedbackError>,
        F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
        M: FnMut(&Cx<'_>, &[f64]) -> Result<(LiquidHeatMarch<'a>, P), E>,
    {
        checkpoint(cx).map_err(LiquidFeedbackError::from)?;
        local.admit()?;
        let descriptions: Vec<_> = self.sections.iter().map(|s| FlowHeatSection {
            region: &s.name, geometry: s.geometry, correlation: s.correlation,
            direction: fs_convection::ThermalDirection::HeatingFluid,
        }).collect();
        let material = self.law.reference_state();
        let template = FrozenFlowMedium::Liquid(material).channel_series(
            VolumetricFlowRate::new(self.mass_flow / material.density_kg_m3()),
            Temperature::new(self.inlet.temperature_k()), &descriptions,
        ).map_err(LiquidFeedbackError::from)?;
        admit(core::slice::from_ref(template.path()), config, initial_walls_k)
            .map_err(LiquidFeedbackError::from)?;
        let mut accelerator = IqnIls::new(self.sections.len(), acceleration)
            .map_err(iqn_error).map_err(LiquidFeedbackError::from)?;
        let mut relaxer = match config.relaxation {
            Relaxation::Fixed { .. } => None,
            Relaxation::Aitken { omega_init, omega_max } => Some(AitkenRelaxation::new(omega_init, omega_max)),
        };
        let total_area: f64 = self.sections.iter().map(|s| s.geometry.wetted_area.value()).sum();
        let mut walls = initial_walls_k.to_vec();
        let mut last = [0.0_f64; 4];
        for iteration in 0..config.max_iterations {
            let (boundary_state, boundary_payload) = march(cx, &walls)?;
            let rows: Vec<_> = boundary_state.regions.iter().map(|r| r.boundary.clone()).collect();
            let references: Vec<_> = rows.iter().map(|r| r.reference_temperature().value()).collect();
            checkpoint(cx).map_err(LiquidFeedbackError::from)?;
            let response = solid(cx, &rows).map_err(LiquidFeedbackError::from)?;
            checkpoint(cx).map_err(LiquidFeedbackError::from)?;
            // Reuse the existing source-bound response ownership/area/reference gate.
            template.check_enthalpy_solid_response(&references, &response.regions)
                .map_err(LiquidFeedbackError::from)?;
            let next: Vec<_> = response.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
            let (fluid, fluid_payload) = march(cx, &next)?;
            let mut weighted = 0.0;
            let mut failed_balance = None;
            last = [0.0; 4];
            for (index, ((before, after), body)) in boundary_state.regions.iter()
                .zip(&fluid.regions).zip(&response.regions).enumerate()
            {
                let wall_delta = checked("material wall residual", next[index] - walls[index])
                    .map_err(LiquidFeedbackError::from)?;
                weighted = checked("material weighted wall residual", weighted
                    + self.sections[index].geometry.wetted_area.value() / total_area * wall_delta)
                    .map_err(LiquidFeedbackError::from)?;
                last[0] = last[0].max(wall_delta.abs());
                let reference_delta = checked("material coolant reference residual",
                    after.boundary.reference_temperature().value() - references[index])
                    .map_err(LiquidFeedbackError::from)?;
                last[1] = last[1].max(reference_delta.abs());
                let old_h = before.boundary.coefficient().value.value();
                let new_h = after.boundary.coefficient().value.value();
                last[2] = last[2].max((new_h - old_h).abs() / new_h.max(old_h));
                let imbalance = checked("material local heat imbalance", body.heat_rate_w - after.heat_rate_w)
                    .map_err(LiquidFeedbackError::from)?.abs();
                last[3] = last[3].max(imbalance);
                let threshold = config.balance_tolerance_w.max(config.balance_relative_tolerance
                    * body.heat_rate_w.abs().max(after.heat_rate_w.abs()));
                if imbalance > threshold { failed_balance = Some((imbalance, threshold)); }
            }
            let temperature_closed = last[0] <= config.temperature_tolerance_k
                && last[1] <= config.temperature_tolerance_k;
            let coefficients_closed = last[2] <= local.coefficient_relative_tolerance;
            if temperature_closed && coefficients_closed && failed_balance.is_none() {
                checkpoint(cx).map_err(LiquidFeedbackError::from)?;
                return Ok((CoupledLiquidFeedbackSolution {
                    solid: response, boundary_state, fluid, iterations: iteration + 1,
                    maximum_wall_residual_k: last[0], maximum_reference_residual_k: last[1],
                    maximum_coefficient_residual_relative: last[2], maximum_heat_imbalance_w: last[3],
                }, boundary_payload, fluid_payload));
            }
            if (last[0] == 0.0 || iteration + 1 == config.max_iterations)
                && temperature_closed && coefficients_closed
                && let Some((imbalance, threshold)) = failed_balance
            {
                return Err(LiquidFeedbackError::from(AirflowError::ConjugateBalanceUnclosed {
                    iterations: iteration + 1, max_region_imbalance_bits: imbalance.to_bits(), tolerance_bits: threshold.to_bits(),
                }).into());
            }
            let omega = match (&config.relaxation, relaxer.as_mut()) {
                (Relaxation::Fixed { omega }, _) => *omega,
                (Relaxation::Aitken { .. }, Some(r)) => r.next_omega(weighted),
                (Relaxation::Aitken { omega_init, .. }, None) => *omega_init,
            };
            let mut proposed = accelerator.step(&walls, &next, omega)
                .map_err(iqn_error).map_err(LiquidFeedbackError::from)?.values;
            let knots = self.source.curve().knots();
            let inside = |t: &f64| t.is_finite() && *t >= knots[0].temperature_k
                && *t <= knots[knots.len()-1].temperature_k;
            if !proposed.iter().all(inside) {
                accelerator.reset();
                proposed = walls.iter().zip(&next).map(|(&old, &new)| old + omega * (new-old)).collect();
                if !proposed.iter().all(inside) {
                    return Err(LiquidFeedbackError::Invalid("configured relaxation leaves the supported wall-temperature span").into());
                }
            }
            checkpoint(cx).map_err(LiquidFeedbackError::from)?;
            walls = proposed;
        }
        Err(LiquidFeedbackError::CouplingNotConverged {
            iterations: config.max_iterations, maximum_wall_residual_k: last[0],
            maximum_reference_residual_k: last[1], maximum_coefficient_residual_relative: last[2],
            maximum_heat_imbalance_w: last[3],
        }.into())
    }
}
