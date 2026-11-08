//! Constant applied-source and linear-bath flow on the original phase chart.

use super::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseHeatExchangeError,
    PhaseStateError, ResponseAccumulator, add, invalid,
};
use crate::thermal::transport::isobaric::{
    source_driven_heat_exchange, source_driven_heat_exchange_duration,
};

/// One immutable source-driven proposal, including exact latent residence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseDrivenHeatExchange {
    pub initial: EquilibriumPhaseState,
    pub state: EquilibriumPhaseState,
    /// Initial signed driving power integrated over the full window [J/kg].
    pub initial_drive_j_kg: f64,
    /// Nonnegative conductance-time per mass [J/(kg K)].
    pub exposure_j_kg_k: f64,
    pub specific_heat_j_kg: f64,
    pub represented_enthalpy_change_j_kg: f64,
    pub heat_residual_j_kg: f64,
    /// Endpoint minus initial T, before absolute-temperature rounding [K].
    pub temperature_change_k: f64,
    pub mean_temperature_k: f64,
    /// Time mean minus initial T, before absolute-temperature rounding [K].
    pub mean_temperature_change_k: f64,
}

/// Derivatives of continuous enthalpy flow, before native-state quantization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseDrivenHeatExchangeJacobian {
    /// Physical fixed-source/fixed-bath d(h_out)/d(h_in).
    pub inlet_enthalpy_gain: f64,
    pub inlet_enthalpy_loss: f64,
    /// d(h_out)/d(Q0) at fixed initial h and conductance exposure.
    pub source_heat_gain: f64,
    /// d(T_mean)/d(Q0) [kg K/J].
    pub mean_source_derivative_kg_k_j: f64,
    /// Q0-X*(T_out-T_initial), retained through the physical attenuation [J/kg].
    pub remaining_drive_j_kg: f64,
    /// d(T_mean)/d(log duration), scaling Q0 and X together [K].
    pub mean_log_duration_derivative_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseDrivenHeatExchangeResponse {
    pub exchange: PhaseDrivenHeatExchange,
    pub jacobian: PhaseDrivenHeatExchangeJacobian,
}

impl EquilibriumEnthalpyPhaseCurve {
    /// Exact dh/ds = Q0-X*(T(h)-T_initial), 0 <= s <= 1.
    ///
    /// The caller supplies Q0=dt*(P+sum G_i*(T_i-T_initial))/m and
    /// X=dt*sum G_i/m. Actual baths, applied power, fixed mass, mechanical
    /// constraint and the lumped transfer model remain caller declarations.
    /// No absolute virtual bath is constructed. Negative virtual equilibria
    /// and zero conductance are supported while the actual chart path is valid.
    ///
    /// Original affine Cp intervals reuse the shared driven-capacity kernel.
    /// Authored mushy spans retain their original dh/dT; sharp latent plateaus
    /// use constant-temperature enthalpy flow exactly once. Traversal allocates
    /// nothing and checks cancellation at each original interval. Completely
    /// lost native-h increments refuse; partial rounding remains in the residual.
    pub fn driven_heat_exchange(
        &self,
        initial: EquilibriumPhaseState,
        initial_drive_j_kg: f64,
        exposure_j_kg_k: f64,
    ) -> Result<PhaseDrivenHeatExchange, PhaseHeatExchangeError> {
        self.driven_heat_exchange_with_checkpoint(initial, initial_drive_j_kg, exposure_j_kg_k, || false)
    }

    pub fn driven_heat_exchange_with_checkpoint(
        &self,
        initial: EquilibriumPhaseState,
        initial_drive_j_kg: f64,
        exposure_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseDrivenHeatExchange, PhaseHeatExchangeError> {
        self.driven_exchange_impl(initial, initial_drive_j_kg, exposure_j_kg_k, None, cancelled)
    }

    /// The same proposal with analytic source, inlet and time-mean response.
    ///
    /// T(h) is continuous across phase events, so the variational flow has
    /// identity event saltation. Latent residence is included without an
    /// apparent Cp, finite differences or time subdivision. At zero drive an
    /// interior corner with unequal one-sided caloric slopes has no unique
    /// response, even when X=0; only this response API refuses that corner.
    pub fn driven_heat_exchange_response(
        &self,
        initial: EquilibriumPhaseState,
        initial_drive_j_kg: f64,
        exposure_j_kg_k: f64,
    ) -> Result<PhaseDrivenHeatExchangeResponse, PhaseHeatExchangeError> {
        self.driven_heat_exchange_response_with_checkpoint(
            initial, initial_drive_j_kg, exposure_j_kg_k, || false,
        )
    }

    pub fn driven_heat_exchange_response_with_checkpoint(
        &self,
        initial: EquilibriumPhaseState,
        initial_drive_j_kg: f64,
        exposure_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseDrivenHeatExchangeResponse, PhaseHeatExchangeError> {
        let mut response = ResponseAccumulator::new();
        let exchange = self.driven_exchange_impl(
            initial, initial_drive_j_kg, exposure_j_kg_k, Some(&mut response), cancelled,
        )?;
        let remaining = initial_drive_j_kg * response.gain;
        if !remaining.is_finite() {
            return Err(invalid("remaining phase driving power is not representable"));
        }
        Ok(PhaseDrivenHeatExchangeResponse {
            exchange,
            jacobian: PhaseDrivenHeatExchangeJacobian {
                inlet_enthalpy_gain: response.gain,
                inlet_enthalpy_loss: response.loss,
                source_heat_gain: if exposure_j_kg_k == 0.0 { 1.0 } else { response.reservoir },
                mean_source_derivative_kg_k_j: response.mean,
                remaining_drive_j_kg: remaining,
                mean_log_duration_derivative_k: response.mean_log_exposure,
            },
        })
    }

    fn driven_exchange_impl(
        &self,
        initial: EquilibriumPhaseState,
        drive: f64,
        feedback: f64,
        mut response: Option<&mut ResponseAccumulator>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseDrivenHeatExchange, PhaseHeatExchangeError> {
        if cancelled() { return Err(PhaseHeatExchangeError::Cancelled); }
        self.advance_specific_energy(initial, 0.0)?;
        if !drive.is_finite() || !feedback.is_finite() || feedback < 0.0 {
            return Err(invalid("phase driven exchange needs finite drive and nonnegative exposure"));
        }
        let (initial_h, initial_t) = (initial.specific_enthalpy_j_kg(), initial.temperature_k());
        if drive == 0.0 {
            if let Some(response) = response.as_deref_mut() {
                self.driven_equilibrium_response(initial, feedback, response)?;
            }
            if cancelled() { return Err(PhaseHeatExchangeError::Cancelled); }
            return Ok(PhaseDrivenHeatExchange {
                initial, state: initial, initial_drive_j_kg: drive, exposure_j_kg_k: feedback,
                specific_heat_j_kg: 0.0, represented_enthalpy_change_j_kg: 0.0,
                heat_residual_j_kg: 0.0, temperature_change_k: 0.0,
                mean_temperature_k: initial_t, mean_temperature_change_k: 0.0,
            });
        }
        // Source-only storage has one authoritative update. Segment travel and
        // residence retain heat/temperature moments; their rounding must not
        // change the final native h or its original support admission.
        let source_state = if feedback == 0.0 {
            Some(self.advance_specific_energy(initial, drive)?)
        } else { None };
        let heating = drive > 0.0;
        let knots = self.knots();
        let upper = if heating {
            knots.partition_point(|k| k.specific_enthalpy_j_kg <= initial_h)
        } else {
            knots.partition_point(|k| k.specific_enthalpy_j_kg < initial_h)
        };
        if upper == 0 || upper == knots.len() {
            return Err(invalid("driven exchange would leave original enthalpy support"));
        }
        let mut index = upper - 1;
        let (mut h, mut temperature) = (initial_h, initial_t);
        let (mut used, mut used_correction) = (0.0, 0.0);
        let (mut heat, mut heat_correction) = (0.0, 0.0);
        let (mut mean_change, mut mean_correction) = (0.0, 0.0);
        for _ in 0..knots.len() - 1 {
            if cancelled() { return Err(PhaseHeatExchangeError::Cancelled); }
            let (lower, upper) = (knots[index], knots[index + 1]);
            let boundary = if heating { upper } else { lower };
            let plateau = lower.temperature_k == upper.temperature_k;
            let driving = (-feedback).mul_add(temperature - initial_t, drive);
            let remaining = 1.0 - used;
            if !driving.is_finite() || driving == 0.0 || (driving > 0.0) != heating
                || !remaining.is_finite() || remaining <= 0.0
            {
                return Err(invalid("remaining phase drive or residence is not representable"));
            }
            let cp = self.driven_sensible_capacities(index);
            let sensible = [
                [lower.temperature_k, cp.map_or(0.0, |v| v[0])],
                [upper.temperature_k, cp.map_or(0.0, |v| v[1])],
            ];
            let to_boundary = if plateau {
                (boundary.specific_enthalpy_j_kg - h) / driving
            } else {
                source_driven_heat_exchange_duration(
                    &sensible, temperature, boundary.temperature_k, driving, feedback,
                )?
            };
            if to_boundary.is_nan() || to_boundary <= 0.0 {
                return Err(invalid("phase boundary residence is not representable"));
            }
            // Exact source endpoints can differ by one residence-time ULP from
            // 1-used. Rounded h alone is insufficient: partial storage
            // quantization must not lengthen the physical source trajectory.
            let source_endpoint = source_state.is_some_and(|state|
                state.specific_enthalpy_j_kg() == boundary.specific_enthalpy_j_kg
                    && drive == boundary.specific_enthalpy_j_kg - initial_h);
            let finished = source_endpoint || remaining <= to_boundary;
            let duration = if finished && !source_endpoint { remaining } else { to_boundary };
            let (local_heat, delta_t, local_mean) = if plateau {
                if let Some(response) = response.as_deref_mut() {
                    response.push_driven(1.0, 0.0, duration, 0.0, 0.0, driving, duration)?;
                }
                (
                    if !finished || remaining == to_boundary || source_endpoint {
                        boundary.specific_enthalpy_j_kg - h
                    } else {
                        duration * driving
                    },
                    0.0, 0.0,
                )
            } else {
                let local = source_driven_heat_exchange(
                    &sensible, temperature, driving, feedback, duration,
                )?;
                if let Some(response) = response.as_deref_mut() {
                    response.push_driven(
                        local.inlet_energy_gain, local.inlet_energy_loss,
                        local.source_gain, local.mean_source_derivative,
                        local.temperature_change, driving, duration,
                    )?;
                    response.push_temperature_moment(
                        used, duration, 1.0, local.temperature_change,
                        local.mean_log_duration_derivative,
                    )?;
                }
                (local.energy_change, local.temperature_change, local.mean_temperature_change)
            };
            if !local_heat.is_finite() || local_heat == 0.0 || (local_heat > 0.0) != heating {
                return Err(invalid("nonzero phase-driven heat is not representable"));
            }
            let next_h = if finished && let Some(state) = source_state {
                state.specific_enthalpy_j_kg()
            } else if !finished || remaining == to_boundary {
                boundary.specific_enthalpy_j_kg
            } else {
                let next = h + local_heat;
                if next == h { return Err(PhaseStateError::LostSpecificEnergyIncrement.into()); }
                next
            };
            if !next_h.is_finite() || next_h < lower.specific_enthalpy_j_kg
                || next_h > upper.specific_enthalpy_j_kg
            {
                return Err(invalid("rounded driven enthalpy leaves its original interval"));
            }
            add(&mut heat, &mut heat_correction, local_heat)?;
            add(&mut mean_change, &mut mean_correction,
                duration * ((temperature - initial_t) + local_mean))?;
            if finished {
                if feedback == 0.0 { heat = drive; }
                let represented = next_h - initial_h;
                if represented == 0.0 { return Err(PhaseStateError::LostSpecificEnergyIncrement.into()); }
                let state = if let Some(state) = source_state { state }
                    else { self.state_at_specific_enthalpy(next_h)? };
                let mean = initial_t + mean_change;
                let delta = (temperature - initial_t) + delta_t;
                let residual = heat - represented;
                let roundoff = 8.0 * f64::EPSILON * initial_t.max(state.temperature_k());
                if !represented.is_finite() || !mean.is_finite() || !residual.is_finite()
                    || mean < initial_t.min(state.temperature_k()) - roundoff
                    || mean > initial_t.max(state.temperature_k()) + roundoff
                {
                    return Err(invalid("phase-driven state or mean violates its heat account"));
                }
                if cancelled() { return Err(PhaseHeatExchangeError::Cancelled); }
                return Ok(PhaseDrivenHeatExchange {
                    initial, state, initial_drive_j_kg: drive, exposure_j_kg_k: feedback,
                    specific_heat_j_kg: heat, represented_enthalpy_change_j_kg: represented,
                    heat_residual_j_kg: residual, temperature_change_k: delta,
                    mean_temperature_k: mean, mean_temperature_change_k: mean_change,
                });
            }
            add(&mut used, &mut used_correction, duration)?;
            h = boundary.specific_enthalpy_j_kg;
            temperature = boundary.temperature_k;
            if heating {
                if index + 2 == knots.len() { break; }
                index += 1;
            } else {
                if index == 0 { break; }
                index -= 1;
            }
        }
        Err(invalid("driven exchange exceeds original enthalpy support"))
    }

    fn driven_sensible_capacities(&self, index: usize) -> Option<[f64; 2]> {
        let (a, b) = (self.knots()[index], self.knots()[index + 1]);
        (a.temperature_k != b.temperature_k).then(|| {
            self.sensible_heat_capacities_on_interval(index).unwrap_or(
                [(b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg)
                    / (b.temperature_k - a.temperature_k); 2],
            )
        })
    }

    fn driven_equilibrium_response(
        &self,
        initial: EquilibriumPhaseState,
        feedback: f64,
        response: &mut ResponseAccumulator,
    ) -> Result<(), PhaseHeatExchangeError> {
        let knots = self.knots();
        let h = initial.specific_enthalpy_j_kg();
        let index = knots.partition_point(|k| k.specific_enthalpy_j_kg <= h)
            .saturating_sub(1).min(knots.len() - 2);
        let current = self.driven_sensible_capacities(index);
        if index > 0 && h == knots[index].specific_enthalpy_j_kg
            && self.driven_sensible_capacities(index - 1).map(|cp| cp[1])
                != current.map(|cp| cp[0])
        {
            return Err(invalid("zero-drive phase corner has unequal one-sided caloric responses"));
        }
        if let Some(cp) = current {
            let local = source_driven_heat_exchange(
                &[[knots[index].temperature_k, cp[0]], [knots[index + 1].temperature_k, cp[1]]],
                initial.temperature_k(), 0.0, feedback, 1.0,
            )?;
            response.gain = local.inlet_energy_gain;
            response.loss = local.inlet_energy_loss;
            response.reservoir = local.source_gain;
            response.mean = local.mean_source_derivative;
        } else {
            response.reservoir = 1.0;
        }
        Ok(())
    }
}

impl ResponseAccumulator {
    #[allow(clippy::too_many_arguments)]
    fn push_driven(
        &mut self, gain: f64, loss: f64, source: f64, mean: f64,
        delta_t: f64, driving: f64, duration: f64,
    ) -> Result<(), PhaseHeatExchangeError> {
        let previous = if self.reservoir == 0.0 {
            0.0
        } else {
            let ratio = delta_t / driving;
            if ratio.is_finite() { ratio * self.reservoir } else { (delta_t * self.reservoir) / driving }
        };
        add(&mut self.mean, &mut self.mean_correction, previous + duration * mean)?;
        self.reservoir = gain.mul_add(self.reservoir, source);
        self.loss = gain.mul_add(self.loss, loss);
        self.gain *= gain;
        if !self.reservoir.is_finite() || self.reservoir < 0.0
            || !self.gain.is_finite() || self.gain < 0.0
            || !self.loss.is_finite() || self.loss < 0.0 || self.mean < 0.0
        {
            return Err(invalid("phase-driven variational response is not representable"));
        }
        Ok(())
    }
}

