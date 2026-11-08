//! Steady adiabatic water/steam mixing at one common pressure.
//!
//! Mass and original-source enthalpy flux determine one equilibrium outlet.
//! Steam condensation, vaporization and sensible heating remain in the original
//! p,h inverse. No temperature average or separate latent-heat source is used.
//! Flow work is already in h; an independent u+p*v balance checks that account.
//!
//! This ideal junction has no storage, external heat or shaft work. Inlet and
//! outlet kinetic/potential energy and momentum are explicitly neglected. It
//! supplies no pressure-loss, entrainment, finite-rate condensation or hydraulic
//! flow law. Pressure-changing devices must prepare streams on the same chart
//! before mixing; a common nominal pressure alone cannot replace that binding.

use super::{
    AcrossWaterFlowError, WaterFlowChart, WaterFlowError, WaterFlowSolve, WaterFlowState,
    check_closure, poll,
};

/// One accepted pure-water source state and its prescribed positive flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterMixingInlet {
    pub state: WaterFlowState,
    pub mass_flow_kg_per_s: f64,
}

/// Completed equilibrium junction; all rates use the supplied flow scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterStreamMixing {
    pub outlet: WaterFlowState,
    pub mass_flow_kg_per_s: f64,
    pub inlet_enthalpy_rate_w: f64,
    pub outlet_enthalpy_rate_w: f64,
    /// Outgoing minus incoming original-source enthalpy rate, evaluated from
    /// compensated source differences before the gross rates are rounded.
    pub enthalpy_residual_w: f64,
    pub internal_energy_change_rate_w: f64,
    /// Outgoing p*v flux minus the incoming p*v flux, already included in h.
    pub pressure_flow_work_change_rate_w: f64,
    /// Independent Delta U_dot + Delta(p*v)_dot for this adiabatic junction.
    pub first_law_residual_w: f64,
    /// Outgoing minus incoming RAW source entropy rate. Never clipped or
    /// corrected by a chart's pressure-specific thermal entropy primitive.
    pub entropy_generation_w_per_k: f64,
    /// Largest absolute endpoint h-u-p_requested*v defect, in J/kg.
    pub maximum_pressure_closure_defect_j_kg: f64,
    pub enthalpy_iterations: usize,
    pub dense_source_evaluations: usize,
}

impl WaterFlowChart {
    /// Join prescribed water/steam streams using their original enthalpy flux.
    ///
    /// Every inlet must retain this exact prepared chart, including pressure,
    /// source window and numerical/source-matching controls. Supply at least
    /// one positive finite flow; omit inactive ports rather than inventing a
    /// trace flow. A single inlet is an exact source-state bypass.
    ///
    /// The existing p,h inverse owns all phase and source-region selection.
    /// Its gap/ambiguity/domain refusal propagates without a replacement fit.
    /// Physical energy and raw-entropy allowances are per kg of total outlet
    /// throughput, as for the original pressure devices. The chart's inverse
    /// accuracy must support the requested energy tolerance.
    ///
    /// Work is linear in inlet count plus one bounded caloric inverse. Success
    /// allocates no memory. Cancellation and any late balance failure return
    /// no partial junction and never mutate the supplied source states.
    pub fn mix_streams(
        &self,
        inlets: &[WaterMixingInlet],
        controls: WaterFlowSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterStreamMixing, AcrossWaterFlowError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if inlets.is_empty() {
            return Err(WaterFlowError::Invalid("mixing requires at least one inlet").into());
        }
        // Anchor the four caloric accounts before multiplying by throughput.
        // The gross source reference never obscures a small physical change.
        let reference = properties(inlets[0].state);
        // Mass, then reference-relative h, u, p*v and raw s fluxes.
        let mut sums = [Sum::default(); 5];
        let mut maximum_closure = 0.0_f64;
        for inlet in inlets {
            poll(&mut cancelled)?;
            let mass = inlet.mass_flow_kg_per_s;
            if !mass.is_finite() || mass <= 0.0 {
                return Err(
                    WaterFlowError::Invalid("mixing flows must be finite positive kg/s").into(),
                );
            }
            if inlet.state.law() != *self {
                return Err(WaterFlowError::Invalid(
                    "mixing requires the same original pressure chart",
                )
                .into());
            }
            let state = inlet.state;
            let defect = state.pressure_closure_defect_j_kg();
            check_closure([defect, defect], controls)?;
            maximum_closure = maximum_closure.max(defect.abs());
            sums[0].add(mass);
            for ((sum, value), anchor) in sums[1..].iter_mut().zip(properties(state)).zip(reference)
            {
                let value = mass * (value - anchor);
                sum.add(finite(value, "mixing inlet flux")?);
            }
        }
        for sum in sums {
            finite(sum.value(), "mixing total flux")?;
        }
        let mass = sums[0].value();
        let target_h = finite(
            reference[0] + (sums[1].sum / mass + sums[1].correction / mass),
            "mixing specific enthalpy",
        )?;
        let (outlet, enthalpy_iterations, dense_source_evaluations) = if inlets.len() == 1 {
            (inlets[0].state, 0, 0)
        } else {
            self.enthalpy(target_h, &mut cancelled)?
        };
        check_closure(
            [maximum_closure, outlet.pressure_closure_defect_j_kg()],
            controls,
        )?;
        maximum_closure = maximum_closure.max(outlet.pressure_closure_defect_j_kg().abs());
        let outlet_h = finite(
            mass * outlet.specific_enthalpy_j_kg(),
            "mixing outlet enthalpy rate",
        )?;
        let mut inlet_h = sums[1];
        inlet_h.add(finite(
            mass * reference[0],
            "mixing reference enthalpy rate",
        )?);
        let inlet_h = finite(inlet_h.value(), "mixing inlet enthalpy rate")?;
        let output = properties(outlet);
        let residuals: [Sum; 4] = core::array::from_fn(|i| {
            let mut residual = Sum::default();
            residual.add(mass * (output[i] - reference[i]));
            residual.add(-sums[i + 1].sum);
            residual.add(-sums[i + 1].correction);
            residual
        });
        let [dh, du, dpv, entropy] = residuals.map(Sum::value);
        for value in [dh, du, dpv, entropy] {
            finite(value, "mixing original-source balance")?;
        }
        let mut first_law = residuals[1];
        first_law.add(residuals[2].sum);
        first_law.add(residuals[2].correction);
        let first_law = finite(first_law.value(), "mixing first law")?;
        let enthalpy_residual = finite(dh / mass, "mixing specific enthalpy residual")?;
        let first_law_residual = finite(first_law / mass, "mixing specific first law")?;
        if enthalpy_residual.abs() > controls.energy_tolerance_j_kg
            || first_law_residual.abs() > controls.energy_tolerance_j_kg
        {
            return Err(WaterFlowError::EnergyBalance {
                enthalpy_residual_j_kg: enthalpy_residual,
                first_law_residual_j_kg: first_law_residual,
            }
            .into());
        }
        let specific_entropy = finite(entropy / mass, "mixing specific entropy production")?;
        if specific_entropy < -controls.entropy_tolerance_j_kg_k {
            return Err(WaterFlowError::EntropyDecrease {
                change_j_kg_k: specific_entropy,
            }
            .into());
        }
        poll(&mut cancelled)?;
        Ok(WaterStreamMixing {
            outlet,
            mass_flow_kg_per_s: mass,
            inlet_enthalpy_rate_w: inlet_h,
            outlet_enthalpy_rate_w: outlet_h,
            enthalpy_residual_w: dh,
            internal_energy_change_rate_w: du,
            pressure_flow_work_change_rate_w: dpv,
            first_law_residual_w: first_law,
            entropy_generation_w_per_k: entropy,
            maximum_pressure_closure_defect_j_kg: maximum_closure,
            enthalpy_iterations,
            dense_source_evaluations,
        })
    }
}

fn properties(state: WaterFlowState) -> [f64; 4] {
    [
        state.specific_enthalpy_j_kg(),
        state.specific_internal_energy_j_kg(),
        state.pressure_pa() * state.specific_volume_m3_kg(),
        state.specific_entropy_j_kg_k(),
    ]
}

fn finite(value: f64, stage: &'static str) -> Result<f64, AcrossWaterFlowError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WaterFlowError::Unrepresentable(stage).into())
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests;
