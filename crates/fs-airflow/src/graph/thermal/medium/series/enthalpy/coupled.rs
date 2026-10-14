//! Close source-backed coolant enthalpy against an actual solid solve.
//!
//! Reuses the existing IQN-ILS accelerator, conjugate configuration admission,
//! typed convection-to-conduction boundary lowering and enthalpy march. No
//! alternate caloric law or solid discretization is introduced. Each callback
//! evaluates a candidate from the same accepted solid state; advancing that
//! accepted state during trial iterations would describe a different problem.

use fs_convection::CorrelationRobinBoundary;
use fs_couple::AitkenRelaxation;
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig};
use fs_exec::Cx;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_qty::Temperature;

use crate::AirflowError;
use crate::conjugate::{ConjugateConfig, Relaxation, SolidRegionState};
use crate::graph::thermal::{admit, checked, checkpoint, iqn_error};
use super::{LiquidEnthalpyError, LiquidEnthalpyMarch, MaterialFlowPath};

/// One candidate solid solve plus its original, unerased solver state/report.
#[derive(Debug)]
pub struct LiquidSolidResponse<R> {
    /// One region in the supplied path order, including the solved reference.
    pub regions: Vec<SolidRegionState>,
    /// The actual solid field, transient-step proposal or solver report.
    pub state: R,
}

/// Accepted common solid/coolant fixed point. No partial candidate is returned.
#[derive(Debug)]
pub struct CoupledLiquidEnthalpySolution<'a, 'c, R> {
    /// Original solid result evaluated at `boundary_rows`.
    pub solid: LiquidSolidResponse<R>,
    /// Source-backed coolant trajectory at the returned solid wall temperatures.
    pub fluid: LiquidEnthalpyMarch<'a, 'c>,
    /// Exact rows sent to the accepted solid solve, not the next proposed rows.
    pub boundary_rows: Vec<CorrelationRobinBoundary>,
    /// Number of actual solid solves performed.
    pub iterations: usize,
    /// Maximum fixed-point temperature residual before relaxation, K.
    pub maximum_reference_residual_k: f64,
    /// Maximum absolute regional solid/fluid heat-rate difference, W.
    pub maximum_heat_imbalance_w: f64,
}

impl<'a> MaterialFlowPath<'a> {
    /// Couple a solid solve to the SAME material-bound sensible-enthalpy march.
    ///
    /// The callback receives concrete, evidence-bearing Robin rows in path
    /// order. It returns the original solid payload and region summaries with
    /// `mean_reference_temperature_k = Some(sent_reference)`. Every candidate
    /// must start from the same accepted solid history. The payload is published
    /// only when both the temperature residual and EACH region's watt balance
    /// pass. Opposing regional heat errors cannot cancel into acceptance.
    ///
    /// Geometry, hydraulic flow and h remain frozen as in `march_liquid_enthalpy`;
    /// IQN acceleration changes convergence, not source/phase applicability.
    /// `config.relaxation` supplies the existing fixed/Aitken fallback. No
    /// guarantee of convergence, calibrated physical accuracy or time accuracy
    /// is added by closing this discrete interface.
    ///
    /// # Errors
    /// Propagates source, transport, solid, cancellation and accelerator errors.
    /// Refuses mismatched region names/areas/references, nonphysical results,
    /// unclosed regional heat or exhaustion of the explicit iteration limit.
    pub fn solve_liquid_enthalpy_iqn<'c, R, F>(
        &self, cx: &Cx<'_>, source: &'c ResolvedHeatCapacityCurve, temperature_axis: &str,
        config: &ConjugateConfig, acceleration: IqnIlsConfig, mut solid: F,
    ) -> Result<CoupledLiquidEnthalpySolution<'a, 'c, R>, LiquidEnthalpyError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        let mut reference = vec![self.path.inlet_temperature_k(); self.exchanges.len()];
        checkpoint(cx, 0, &reference)?;
        // This module is a descendant of the existing thermal driver: use its
        // actual configuration admission and accelerator-error mapping.
        admit(core::slice::from_ref(&self.path), config, &reference)?;
        self.admit_caloric_source(cx, source, temperature_axis)?;
        super::state_at_temperature(source.curve(), self.path.inlet_temperature_k())?;
        let mut accelerator = IqnIls::new(reference.len(), acceleration).map_err(iqn_error)?;
        let mut relaxer = match config.relaxation {
            Relaxation::Fixed { .. } => None,
            Relaxation::Aitken { omega_init, omega_max } => Some(AitkenRelaxation::new(omega_init, omega_max)),
        };
        let total_area: f64 = self.path.segments().iter().map(|s| s.area_m2()).sum();
        let mut last_change = 0.0_f64;
        for iteration in 0..config.max_iterations {
            checkpoint(cx, iteration, &reference)?;
            let rows = self.exchanges.iter().zip(&reference)
                .map(|(exchange, &t)| exchange.robin_boundary(Temperature::new(t)))
                .collect::<Result<Vec<_>, _>>().map_err(LiquidEnthalpyError::Flow)?;
            let response = solid(cx, &rows)?;
            checkpoint(cx, iteration, &reference)?;
            self.check_enthalpy_solid_response(&reference, &response.regions)?;
            let walls: Vec<_> = response.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
            let fluid = self.march_liquid_enthalpy(cx, source, temperature_axis, &walls)?;
            let mut updated = Vec::with_capacity(reference.len());
            let mut weighted = 0.0;
            last_change = 0.0;
            let mut maximum_heat_imbalance = 0.0_f64;
            let mut failed_balance = None;
            for (index, ((exchange, body), &now)) in fluid.regions().iter()
                .zip(&response.regions).zip(&reference).enumerate()
            {
                let next = exchange.reference_temperature_k();
                let delta = checked("liquid reference residual", next - now)?;
                last_change = last_change.max(delta.abs());
                weighted = checked("liquid weighted reference residual", weighted
                    + (self.path.segments()[index].area_m2() / total_area) * delta)?;
                let imbalance = checked("liquid regional heat imbalance", body.heat_rate_w - exchange.heat_rate_w())?.abs();
                maximum_heat_imbalance = maximum_heat_imbalance.max(imbalance);
                let scale = body.heat_rate_w.abs().max(exchange.heat_rate_w().abs());
                let threshold = config.balance_tolerance_w.max(config.balance_relative_tolerance * scale);
                if imbalance > threshold && failed_balance.is_none_or(|(worst, _)| imbalance > worst) {
                    failed_balance = Some((imbalance, threshold));
                }
                updated.push(next);
            }
            if last_change <= config.temperature_tolerance_k {
                checkpoint(cx, iteration, &reference)?;
                if let Some((imbalance, threshold)) = failed_balance {
                    return Err(AirflowError::ConjugateBalanceUnclosed {
                        iterations: iteration + 1,
                        max_region_imbalance_bits: imbalance.to_bits(), tolerance_bits: threshold.to_bits(),
                    }.into());
                }
                return Ok(CoupledLiquidEnthalpySolution {
                    solid: response, fluid, boundary_rows: rows, iterations: iteration + 1,
                    maximum_reference_residual_k: last_change, maximum_heat_imbalance_w: maximum_heat_imbalance,
                });
            }
            let omega = match (&config.relaxation, relaxer.as_mut()) {
                (Relaxation::Fixed { omega }, _) => *omega,
                (Relaxation::Aitken { .. }, Some(relaxer)) => relaxer.next_omega(weighted),
                (Relaxation::Aitken { omega_init, .. }, None) => *omega_init,
            };
            reference = accelerator.step(&reference, &updated, omega).map_err(iqn_error)?.values;
        }
        Err(AirflowError::ConjugateNotConverged { iterations: config.max_iterations,
            max_change_bits: last_change.to_bits(), tolerance_bits: config.temperature_tolerance_k.to_bits() }.into())
    }

    fn check_enthalpy_solid_response(
        &self, references: &[f64], response: &[SolidRegionState],
    ) -> Result<(), LiquidEnthalpyError> {
        if response.len() != self.path.segments().len() {
            return Err(AirflowError::SolidResponseArity { expected: self.path.segments().len(), found: response.len() }.into());
        }
        for (index, (segment, body)) in self.path.segments().iter().zip(response).enumerate() {
            if segment.region() != body.region {
                return Err(AirflowError::SegmentRegionMismatch { index, expected: segment.region().to_owned(), found: body.region.clone() }.into());
            }
            if !super::positive(body.mean_wall_temperature_k) || !super::positive(body.area_m2)
                || !body.heat_rate_w.is_finite()
            {
                return Err(LiquidEnthalpyError::Invalid("solid reply must have positive finite temperature/area and finite heat"));
            }
            let area = segment.area_m2();
            if (body.area_m2 - area).abs() > 1e-9 * area {
                return Err(AirflowError::SegmentAreaMismatch { index, region: body.region.clone(),
                    declared_bits: area.to_bits(), reported_bits: body.area_m2.to_bits() }.into());
            }
            let reported = body.mean_reference_temperature_k.ok_or(LiquidEnthalpyError::Invalid(
                "enthalpy coupling requires the solid's solved-against Robin reference"))?;
            if !reported.is_finite() || (reported - references[index]).abs() > 1e-9 * references[index].abs().max(1.0) {
                return Err(AirflowError::ReferenceTemperatureMismatch { index, region: body.region.clone(),
                    sent_bits: references[index].to_bits(), reported_bits: reported.to_bits() }.into());
            }
        }
        Ok(())
    }
}
