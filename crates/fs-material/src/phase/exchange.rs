//! Fixed-reservoir exchange on the original sensible and latent enthalpy path.
//!
//! The caller supplies exposure G*dt/m [J/(kg K)]. Each sensible interval reuses
//! the original affine-Cp exchange kernel; an isothermal interval consumes heat
//! at its actual constant driving temperature. No apparent heat capacity,
//! temporal subdivision, extra latent source or extrapolated material is used.

use super::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use crate::thermal::transport::{
    ThermalTransportError,
    isobaric::{source_linear_heat_exchange, source_linear_heat_exchange_exposure},
};

/// Refusal leaves the source chart and supplied accepted state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseHeatExchangeError {
    Invalid(&'static str),
    Phase(PhaseStateError),
    Caloric(ThermalTransportError),
    Cancelled,
}

impl core::fmt::Display for PhaseHeatExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "phase heat exchange: {what}"),
            Self::Phase(error) => core::fmt::Display::fmt(error, f),
            Self::Caloric(error) => core::fmt::Display::fmt(error, f),
            Self::Cancelled => write!(f, "phase heat exchange cancelled"),
        }
    }
}
impl std::error::Error for PhaseHeatExchangeError {}
impl From<PhaseStateError> for PhaseHeatExchangeError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<ThermalTransportError> for PhaseHeatExchangeError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Caloric(error)
    }
}

/// One proposed exchange, retaining the original material and phase identities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseHeatExchange {
    pub initial: EquilibriumPhaseState,
    pub state: EquilibriumPhaseState,
    pub reservoir_temperature_k: f64,
    /// Total conductance-time per fixed material mass [J/(kg K)].
    pub exposure_j_kg_k: f64,
    /// Signed integrated heat per material mass, positive into the body [J/kg].
    /// This is computed from relative source integrals, without subtracting
    /// large absolute enthalpy primitives. Latent heat appears exactly once.
    pub specific_heat_j_kg: f64,
    /// Actual represented endpoint difference in the original h coordinate.
    pub represented_enthalpy_change_j_kg: f64,
    /// Integrated heat minus the represented endpoint difference [J/kg].
    /// The caller must admit this against its physical energy budget.
    pub heat_residual_j_kg: f64,
    /// Exposure-weighted material temperature [K]. With fixed G and mass this
    /// is also the time mean. Zero exposure returns the initial temperature.
    pub mean_temperature_k: f64,
}

/// Jacobian of the continuous original caloric flow, before absolute h is rounded.
///
/// All unmentioned chart data and inputs remain fixed. These are enthalpy and
/// integrated-heat derivatives; T_out need not be differentiable at a phase edge.
/// They do not differentiate floating-point state quantization. The exchange's
/// represented h change and actual energy residual remain separately binding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseHeatExchangeJacobian {
    /// d h_out / d h_in. Nonnegative attenuation; may underflow to zero.
    pub inlet_enthalpy_gain: f64,
    /// Stable 1 - gain, so d q / d h_in = -loss for heat per fixed mass.
    pub inlet_enthalpy_loss: f64,
    /// d h_out / d T_reservoir = d q / d T_reservoir [J/(kg K)].
    pub reservoir_enthalpy_derivative_j_kg_k: f64,
    /// d h_out / d X = d q / d X [K]. Retains the final driving temperature
    /// through its attenuation even when the absolute T_out rounds to the bath.
    pub exposure_enthalpy_derivative_k: f64,
    /// d T_mean / d T_reservoir. A positive moment preserves weak exchange.
    pub mean_reservoir_derivative: f64,
}

/// One exact exchange and its analytic response, from the same chart traversal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseHeatExchangeResponse {
    pub exchange: PhaseHeatExchange,
    pub jacobian: PhaseHeatExchangeJacobian,
}

impl EquilibriumEnthalpyPhaseCurve {
    /// Exchange against a fixed-temperature reservoir through the full chart.
    ///
    /// Solves dh/dX = T_reservoir - T(h), where X = G*dt/m. The original
    /// affine-Cp sensible intervals, authored affine T(h) intervals (including
    /// mushy spans), and sharp latent plateaus retain their own exact laws.
    /// Starting on a plateau therefore requires an enthalpy state, not just T.
    ///
    /// The caller owns fixed pressure/composition, fixed mass, spatial lumping
    /// and the constant conductance reduction. No radiation, other heat source,
    /// pressure evolution or independent pressure-work account is inferred.
    /// The reservoir may be outside source support; the material may not.
    ///
    /// Visits each crossed chart interval once, with a bounded final sensible
    /// inverse and no allocation. Compact exact-source ingress avoids optional
    /// chord subdivisions. Returned relative heat does not make a rounded
    /// absolute h state independent of an arbitrarily large datum: completely
    /// lost updates refuse and partial quantization remains in the residual.
    pub fn linear_heat_exchange(
        &self,
        initial: EquilibriumPhaseState,
        reservoir_k: f64,
        exposure_j_kg_k: f64,
    ) -> Result<PhaseHeatExchange, PhaseHeatExchangeError> {
        self.linear_heat_exchange_with_checkpoint(initial, reservoir_k, exposure_j_kg_k, || false)
    }

    /// Same exchange with cancellation before work, at each crossed interval,
    /// and before publication. Returning true from the callback cancels.
    pub fn linear_heat_exchange_with_checkpoint(
        &self,
        initial: EquilibriumPhaseState,
        reservoir_k: f64,
        exposure_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseHeatExchange, PhaseHeatExchangeError> {
        self.linear_heat_exchange_impl(initial, reservoir_k, exposure_j_kg_k, None, cancelled)
    }

    /// Exchange and analytic enthalpy/heat response for implicit coupling.
    ///
    /// Sensible response reuses the original affine-Cp kernel. An isothermal
    /// segment has unit inlet gain and bath derivative X, so changes in latent
    /// residence are retained across melting or freezing. Continuous T(h)
    /// makes the switching-event saltation identity: no perturbed marches,
    /// finite differences, apparent Cp or extra transition source are used.
    ///
    /// At positive-time equilibrium, unequal one-sided caloric slopes have no
    /// unique two-sided Jacobian. Only this response API refuses that corner;
    /// the value-only exchange remains admissible. At a chart support endpoint
    /// the derivative is one-sided within the admitted source domain.
    pub fn linear_heat_exchange_response(
        &self,
        initial: EquilibriumPhaseState,
        reservoir_k: f64,
        exposure_j_kg_k: f64,
    ) -> Result<PhaseHeatExchangeResponse, PhaseHeatExchangeError> {
        self.linear_heat_exchange_response_with_checkpoint(
            initial,
            reservoir_k,
            exposure_j_kg_k,
            || false,
        )
    }

    /// Response with the same bounded work and cancellation points as exchange.
    pub fn linear_heat_exchange_response_with_checkpoint(
        &self,
        initial: EquilibriumPhaseState,
        reservoir_k: f64,
        exposure_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseHeatExchangeResponse, PhaseHeatExchangeError> {
        let mut response = ResponseAccumulator::new();
        let exchange = self.linear_heat_exchange_impl(
            initial,
            reservoir_k,
            exposure_j_kg_k,
            Some(&mut response),
            cancelled,
        )?;
        Ok(PhaseHeatExchangeResponse {
            exchange,
            jacobian: response.finish(exchange)?,
        })
    }

    fn linear_heat_exchange_impl(
        &self,
        initial: EquilibriumPhaseState,
        reservoir_k: f64,
        exposure_j_kg_k: f64,
        mut response: Option<&mut ResponseAccumulator>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseHeatExchange, PhaseHeatExchangeError> {
        if cancelled() {
            return Err(PhaseHeatExchangeError::Cancelled);
        }
        self.advance_specific_energy(initial, 0.0)?;
        if !reservoir_k.is_finite()
            || reservoir_k <= 0.0
            || !exposure_j_kg_k.is_finite()
            || exposure_j_kg_k < 0.0
        {
            return Err(invalid(
                "positive finite reservoir temperature and finite nonnegative exposure required",
            ));
        }
        let initial_t = initial.temperature_k();
        let initial_h = initial.specific_enthalpy_j_kg();
        if exposure_j_kg_k == 0.0 || initial_t == reservoir_k {
            if exposure_j_kg_k > 0.0
                && let Some(response) = response.as_deref_mut()
            {
                self.equilibrium_exchange_response(initial, exposure_j_kg_k, response)?;
            }
            if cancelled() {
                return Err(PhaseHeatExchangeError::Cancelled);
            }
            return Ok(PhaseHeatExchange {
                initial,
                state: initial,
                reservoir_temperature_k: reservoir_k,
                exposure_j_kg_k,
                specific_heat_j_kg: 0.0,
                represented_enthalpy_change_j_kg: 0.0,
                heat_residual_j_kg: 0.0,
                mean_temperature_k: initial_t,
            });
        }
        let heating = reservoir_k > initial_t;
        let knots = self.knots();
        let upper = if heating {
            knots.partition_point(|k| k.specific_enthalpy_j_kg <= initial_h)
        } else {
            knots.partition_point(|k| k.specific_enthalpy_j_kg < initial_h)
        };
        if upper == 0 || upper == knots.len() {
            return Err(invalid(
                "heat exchange would leave the original enthalpy support",
            ));
        }
        let mut index = upper - 1;
        let (mut h, mut temperature) = (initial_h, initial_t);
        let (mut used, mut used_correction) = (0.0, 0.0);
        let (mut heat, mut heat_correction) = (0.0, 0.0);
        let (mut mean_change, mut mean_correction) = (0.0, 0.0);
        for _ in 0..knots.len() - 1 {
            if cancelled() {
                return Err(PhaseHeatExchangeError::Cancelled);
            }
            let (lower, upper) = (knots[index], knots[index + 1]);
            let boundary = if heating { upper } else { lower };
            let plateau = lower.temperature_k == upper.temperature_k;
            let driving = reservoir_k - temperature;
            let remaining = exposure_j_kg_k - used;
            if !remaining.is_finite() || remaining <= 0.0 {
                return Err(invalid("remaining exchange exposure is not representable"));
            }
            let capacities = if plateau {
                [0.0; 2]
            } else {
                self.sensible_heat_capacities_on_interval(index)
                    .unwrap_or_else(|| {
                        [(upper.specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg)
                            / (upper.temperature_k - lower.temperature_k);
                            2]
                    })
            };
            let sensible = [
                [lower.temperature_k, capacities[0]],
                [upper.temperature_k, capacities[1]],
            ];
            let to_boundary = if plateau {
                (boundary.specific_enthalpy_j_kg - h) / driving
            } else if (heating && boundary.temperature_k >= reservoir_k)
                || (!heating && boundary.temperature_k <= reservoir_k)
            {
                f64::INFINITY
            } else {
                source_linear_heat_exchange_exposure(
                    &sensible,
                    temperature,
                    boundary.temperature_k,
                    reservoir_k,
                )?
            };
            if to_boundary.is_nan() || to_boundary <= 0.0 {
                return Err(invalid(
                    "exchange to the next chart boundary is not representable",
                ));
            }
            let finished = remaining <= to_boundary;
            let local_exposure = if finished { remaining } else { to_boundary };
            let (local_heat, local_mean) = if plateau {
                if let Some(response) = response.as_deref_mut() {
                    response.push(
                        1.0,
                        0.0,
                        local_exposure,
                        0.0,
                        local_exposure,
                        exposure_j_kg_k,
                    )?;
                }
                (
                    if !finished || remaining == to_boundary {
                        boundary.specific_enthalpy_j_kg - h
                    } else {
                        local_exposure * driving
                    },
                    temperature,
                )
            } else {
                let exchange = source_linear_heat_exchange(
                    &sensible,
                    temperature,
                    reservoir_k,
                    local_exposure,
                )?;
                if let Some(response) = response.as_deref_mut() {
                    response.push(
                        exchange.inlet_enthalpy_gain,
                        exchange.inlet_enthalpy_loss,
                        exchange.reservoir_enthalpy_derivative_j_kg_k,
                        exchange.mean_reservoir_derivative,
                        local_exposure,
                        exposure_j_kg_k,
                    )?;
                }
                (
                    exchange.specific_enthalpy_change_j_kg,
                    exchange.mean_temperature_k,
                )
            };
            if !local_heat.is_finite()
                || local_heat == 0.0
                || (local_heat > 0.0) != heating
                || !local_mean.is_finite()
            {
                return Err(invalid(
                    "nonzero exchange heat or its mean is not representable",
                ));
            }
            let next_h = if !finished || remaining == to_boundary {
                boundary.specific_enthalpy_j_kg
            } else {
                let next = h + local_heat;
                if next == h {
                    return Err(PhaseStateError::LostSpecificEnergyIncrement.into());
                }
                next
            };
            if !next_h.is_finite()
                || next_h < lower.specific_enthalpy_j_kg
                || next_h > upper.specific_enthalpy_j_kg
            {
                return Err(invalid(
                    "rounded exchange enthalpy leaves its original interval",
                ));
            }
            add(&mut heat, &mut heat_correction, local_heat)?;
            add(
                &mut mean_change,
                &mut mean_correction,
                (local_exposure / exposure_j_kg_k) * (local_mean - initial_t),
            )?;
            if finished {
                let represented = next_h - initial_h;
                if represented == 0.0 {
                    return Err(PhaseStateError::LostSpecificEnergyIncrement.into());
                }
                let state = self.state_at_specific_enthalpy(next_h)?;
                let mean = initial_t + mean_change;
                let residual = heat - represented;
                let roundoff = 8.0 * f64::EPSILON * initial_t.max(reservoir_k);
                if !represented.is_finite()
                    || !mean.is_finite()
                    || !residual.is_finite()
                    || state.temperature_k() < initial_t.min(reservoir_k) - roundoff
                    || state.temperature_k() > initial_t.max(reservoir_k) + roundoff
                    || mean < initial_t.min(state.temperature_k()) - roundoff
                    || mean > initial_t.max(state.temperature_k()) + roundoff
                {
                    return Err(invalid(
                        "rounded exchange violates its heat or temperature account",
                    ));
                }
                if cancelled() {
                    return Err(PhaseHeatExchangeError::Cancelled);
                }
                return Ok(PhaseHeatExchange {
                    initial,
                    state,
                    reservoir_temperature_k: reservoir_k,
                    exposure_j_kg_k,
                    specific_heat_j_kg: heat,
                    represented_enthalpy_change_j_kg: represented,
                    heat_residual_j_kg: residual,
                    mean_temperature_k: mean,
                });
            }
            add(&mut used, &mut used_correction, local_exposure)?;
            h = boundary.specific_enthalpy_j_kg;
            temperature = boundary.temperature_k;
            if heating {
                if index + 2 == knots.len() {
                    break;
                }
                index += 1;
            } else {
                if index == 0 {
                    break;
                }
                index -= 1;
            }
        }
        Err(invalid(
            "heat exchange exceeds the original enthalpy support",
        ))
    }
    fn equilibrium_exchange_response(
        &self,
        initial: EquilibriumPhaseState,
        exposure: f64,
        response: &mut ResponseAccumulator,
    ) -> Result<(), PhaseHeatExchangeError> {
        let knots = self.knots();
        let h = initial.specific_enthalpy_j_kg();
        let t = initial.temperature_k();
        let index = knots
            .partition_point(|k| k.specific_enthalpy_j_kg <= h)
            .saturating_sub(1)
            .min(knots.len() - 2);
        let capacities = |i: usize| {
            let (a, b) = (knots[i], knots[i + 1]);
            if a.temperature_k == b.temperature_k {
                None
            } else {
                Some(self.sensible_heat_capacities_on_interval(i).unwrap_or(
                    [(b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg)
                        / (b.temperature_k - a.temperature_k); 2],
                ))
            }
        };
        let current = capacities(index);
        if index > 0
            && h == knots[index].specific_enthalpy_j_kg
            && capacities(index - 1).map(|cp| cp[1]) != current.map(|cp| cp[0])
        {
            return Err(invalid(
                "equilibrium at unequal one-sided caloric slopes has no unique response",
            ));
        }
        if let Some(cp) = current {
            let local = source_linear_heat_exchange(
                &[
                    [knots[index].temperature_k, cp[0]],
                    [knots[index + 1].temperature_k, cp[1]],
                ],
                t,
                t,
                exposure,
            )?;
            response.push(
                local.inlet_enthalpy_gain,
                local.inlet_enthalpy_loss,
                local.reservoir_enthalpy_derivative_j_kg_k,
                local.mean_reservoir_derivative,
                exposure,
                exposure,
            )
        } else {
            response.push(1.0, 0.0, exposure, 0.0, exposure, exposure)
        }
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), PhaseHeatExchangeError> {
    let increment = value - *correction;
    let next = *sum + increment;
    if !value.is_finite() || !next.is_finite() {
        return Err(invalid("integrated exchange quantity is not representable"));
    }
    *correction = (next - *sum) - increment;
    *sum = next;
    Ok(())
}

// Compose the exact local variational flow. T(h) is continuous at every chart
// boundary, so event-time terms cancel; the saltation multiplier is one.
struct ResponseAccumulator {
    gain: f64,
    loss: f64,
    reservoir: f64,
    mean: f64,
    mean_correction: f64,
}

impl ResponseAccumulator {
    fn new() -> Self {
        Self {
            gain: 1.0,
            loss: 0.0,
            reservoir: 0.0,
            mean: 0.0,
            mean_correction: 0.0,
        }
    }

    fn push(
        &mut self,
        gain: f64,
        loss: f64,
        reservoir: f64,
        mean: f64,
        local_exposure: f64,
        total_exposure: f64,
    ) -> Result<(), PhaseHeatExchangeError> {
        // X-B is a positive moment. Normalize before multiplying to retain a
        // representable weak response even if the unnormalized moment underflows.
        add(
            &mut self.mean,
            &mut self.mean_correction,
            loss * (self.reservoir / total_exposure) + (local_exposure / total_exposure) * mean,
        )?;
        self.reservoir = gain.mul_add(self.reservoir, reservoir);
        self.loss = gain.mul_add(self.loss, loss);
        self.gain *= gain;
        if !self.reservoir.is_finite()
            || self.reservoir < 0.0
            || !self.loss.is_finite()
            || self.loss < 0.0
            || !self.gain.is_finite()
            || self.gain < 0.0
            || self.mean < 0.0
        {
            return Err(invalid("exchange Jacobian is not representable"));
        }
        Ok(())
    }

    fn finish(
        self,
        exchange: PhaseHeatExchange,
    ) -> Result<PhaseHeatExchangeJacobian, PhaseHeatExchangeError> {
        let exposure_derivative =
            (exchange.reservoir_temperature_k - exchange.initial.temperature_k()) * self.gain;
        if !exposure_derivative.is_finite() {
            return Err(invalid("exchange exposure derivative is not representable"));
        }
        Ok(PhaseHeatExchangeJacobian {
            inlet_enthalpy_gain: self.gain,
            inlet_enthalpy_loss: self.loss,
            reservoir_enthalpy_derivative_j_kg_k: self.reservoir,
            exposure_enthalpy_derivative_k: exposure_derivative,
            mean_reservoir_derivative: self.mean,
        })
    }
}

fn invalid(what: &'static str) -> PhaseHeatExchangeError {
    PhaseHeatExchangeError::Invalid(what)
}

#[cfg(test)]
mod tests;
