//! Finite-step thermal entropy admission on already evaluated passive stages.
use super::*;

impl ImplicitConductionChain<'_> {
    pub(super) fn admit_passive_sdirk2_entropy(
        &self,
        thermal: &DepositionSimulation<'_>,
        first: &Endpoint,
        second: &BackwardEulerTrial<'_, Endpoint>,
        dt_s: f64,
        cell_tolerance_j: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        let audit = second
            .entropy_audit()
            .ok_or(ImplicitConductionError::Invalid(
                "missing final entropy account",
            ))?;
        // Never treat a stream or an unclassified source/work term as a heat
        // bath. The caller also excludes generalized boundary-exchange laws,
        // whose signed powers do not assert passive constitutive behavior.
        if audit.has_unclassified_sources || audit.advection.is_some() {
            return Ok(());
        }
        let entropy_change_j_k =
            finite(audit.thermal_storage_change_j_k - audit.reservoir_entropy_input_j_k)?;
        if entropy_change_j_k >= 0.0 {
            return Ok(());
        }

        // Only a negative candidate needs this O(cells + links) coefficient
        // reduction. Reuse original coverage and actual stage conductances;
        // do not integrate calorics again or re-query a nonlinear surface law.
        let exposed = thermal.exposed_areas_m2(&self.ids, &mut *checkpoint)?;
        let stages = [
            (first, dt_s * (1.0 - GAMMA)),
            (second.payload(), dt_s * GAMMA),
        ];
        let mut minimum_temperature_k = f64::INFINITY;
        let mut absolute_heat_j = Sum::default();
        for (i, &id) in self.ids.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let cell = &thermal.cells[&id];
            minimum_temperature_k = minimum_temperature_k.min(thermal.temperature_k(id)?);
            for (stage, duration_s) in stages {
                let temperature_k = stage.states[i].temperature_k;
                minimum_temperature_k = minimum_temperature_k.min(temperature_k);
                let ambient = cell.boundary.ambient_k;
                let conductance = stage
                    .boundary_values
                    .ambient_w_k
                    .get(&id)
                    .copied()
                    .unwrap_or(exposed[&id] * cell.boundary.htc_w_m2_k);
                let heat = finite(duration_s * finite(conductance * (ambient - temperature_k))?)?;
                absolute_heat_j.add(heat.abs())?;
                if heat != 0.0 {
                    minimum_temperature_k = minimum_temperature_k.min(ambient);
                }
                if let Some(reservoir) = cell.boundary.reservoir {
                    let conductance = stage
                        .boundary_values
                        .reservoir_w_k
                        .get(&id)
                        .copied()
                        .unwrap_or(reservoir.conductance_w_k);
                    let heat = finite(
                        duration_s
                            * finite(conductance * (reservoir.temperature_k - temperature_k))?,
                    )?;
                    absolute_heat_j.add(heat.abs())?;
                    if heat != 0.0 {
                        minimum_temperature_k = minimum_temperature_k.min(reservoir.temperature_k);
                    }
                }
            }
        }
        for (i, (&id, link)) in self.conduction.links.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let a = self.index[&link.a];
            let b = self.index[&link.b];
            for (stage, duration_s) in stages {
                let power = finite(
                    stage.link_values[&id]
                        * (stage.states[a].temperature_k - stage.states[b].temperature_k),
                )?;
                absolute_heat_j.add(finite(duration_s * power)?.abs())?;
            }
        }

        // Both native stages already passed absolute LOCAL joule gates. The
        // retained first-stage weight is (1-gamma)/gamma, hence their combined
        // local acceptance envelope is n*tolerance/gamma. Convert that finite
        // energy-increment budget through dS=dH/T, never through absolute h or
        // an entropy datum. A total-joule tolerance alone could hide cancelling
        // cell errors with different reciprocal temperatures.
        let residual_allowance_j_k = finite(
            finite(self.ids.len() as f64 * cell_tolerance_j / GAMMA)? / minimum_temperature_k,
        )?;
        // Each passive heat exchange has two entropy endpoints (two cells, or
        // one cell and its bath). Twice its absolute heat / Tmin bounds the
        // magnitude of their entropy increments, including bath cancellation.
        // Compensated reductions and a 64-epsilon operation allowance keep this
        // scale tied to transported heat, not large reference-energy offsets.
        let roundoff_j_k = finite(
            128.0 * f64::EPSILON * finite(absolute_heat_j.total()? / minimum_temperature_k)?
                + 64.0 * f64::EPSILON * residual_allowance_j_k,
        )?;
        let allowance_j_k = finite(residual_allowance_j_k + roundoff_j_k)?;
        checkpoint()?;
        if entropy_change_j_k < -allowance_j_k {
            return Err(ImplicitConductionError::EntropyDecrease {
                entropy_change_j_k,
                allowance_j_k,
            });
        }
        Ok(())
    }
}
