//! Heat-source evaluation for the original stage and adaptive thermal owner.
//!
//! A prescribed reservoir contributes G dt (T_res - T_end) to each implicit
//! stage. Its derivative adds G dt to the positive heat-balance derivative.
//! This module supplies boundary heat and receipts, not another integrator.

use super::super::Sum;
use super::finite_reservoir::FiniteReservoir;
use super::{
    SolubilityCaloricErrorV1, SolubilityThermalKineticControlsV1,
    SolubilityThermalKineticOutcomeV1, SolubilityThermalKineticStateV1,
    SolubilityThermalKineticsV1, WorkingState, checked, poll, positive,
};
use fs_qty::{Temperature, Time};

impl SolubilityThermalKineticsV1 {
    /// Exchange heat with a prescribed-temperature reservoir through a finite
    /// nonnegative conductance [W/K], using the SAME phase state, scalar heat
    /// solve and first-order adaptive controller as [`Self::advance`].
    ///
    /// Each implicit stage supplies `G*dt*(T_res-T_end)` joules. Only the
    /// accepted two-half-stage powers contribute to the compensated window
    /// receipt; rejected and coarse trials supply no physical heat. The final
    /// energy residual independently compares original-to-final enthalpy with
    /// that integrated receipt. The reservoir is externally prescribed; no
    /// finite reservoir heat capacity or additional thermal state is inferred.
    ///
    /// The reservoir temperature must be positive and finite but may lie
    /// outside the material's source support. Every solved bath temperature
    /// must remain inside the original root bracket and caloric/interaction
    /// domains. An unreachable supported endpoint refuses without clipping
    /// the heat. Zero conductance is the exact existing adiabatic path.
    pub fn advance_with_heat_transfer(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        duration: Time,
        reservoir_temperature: Temperature,
        conductance_w_k: f64,
        controls: SolubilityThermalKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityThermalKineticOutcomeV1, SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        if !positive(reservoir_temperature.value())
            || !conductance_w_k.is_finite()
            || conductance_w_k < 0.0
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "positive finite reservoir temperature and nonnegative finite thermal conductance",
            ));
        }
        if conductance_w_k == 0.0 {
            return self.advance(initial, duration, 0.0, controls, cancelled);
        }
        self.advance_with_source(
            initial,
            duration,
            WindowHeat::Reservoir(Reservoir {
                temperature_k: reservoir_temperature.value(),
                conductance_w_k,
            }),
            controls,
            cancelled,
        )
        .map(|(balance, _)| balance)
    }
}

#[derive(Clone, Copy)]
pub(super) enum WindowHeat {
    Uniform(f64),
    Reservoir(Reservoir),
    Finite(FiniteReservoir),
}

#[derive(Clone, Copy)]
pub(super) struct Reservoir {
    temperature_k: f64,
    conductance_w_k: f64,
}
impl Reservoir {
    fn heat(self, duration: Time, endpoint: Temperature) -> Result<f64, SolubilityCaloricErrorV1> {
        let difference = self.temperature_k - endpoint.value();
        if difference == 0.0 || duration.value() == 0.0 {
            return Ok(0.0);
        }
        // Do not require G*dt or G*deltaT to be representable separately when
        // their complete physical heat product is representable.
        let g = self.conductance_w_k;
        let dt = duration.value();
        for value in [
            (g * dt) * difference,
            (g * difference) * dt,
            (dt * difference) * g,
        ] {
            if value.is_finite() && value != 0.0 {
                return Ok(value);
            }
        }
        Err(SolubilityCaloricErrorV1::Arithmetic(
            "nonzero reservoir stage heat is not representable",
        ))
    }
}

/// Independent sum of actual accepted endpoint-power contributions. Keeping
/// its compensation in trials prevents many small heat increments disappearing
/// into a large accumulated receipt before the final energy comparison.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct HeatReceipt {
    main: f64,
    correction: f64,
}
impl HeatReceipt {
    pub(super) fn from_value(value: f64) -> Self {
        Self {
            main: value,
            correction: 0.0,
        }
    }
    pub(super) fn add(self, value: f64) -> Result<Self, SolubilityCaloricErrorV1> {
        let mut sum = Sum {
            main: self.main,
            correction: self.correction,
        };
        sum.add(value)?;
        Ok(Self {
            main: sum.main,
            correction: sum.correction,
        })
    }
    pub(super) fn value(self) -> Result<f64, SolubilityCaloricErrorV1> {
        Sum {
            main: self.main,
            correction: self.correction,
        }
        .value()
    }
    pub(super) fn subtract(self, other: Self) -> Result<Self, SolubilityCaloricErrorV1> {
        self.add(-other.main)?.add(-other.correction)
    }
    pub(super) fn divided_by(self, divisor: f64) -> Result<f64, SolubilityCaloricErrorV1> {
        self.shifted_quotient(divisor, 0.0)
    }
    pub(super) fn shifted_quotient(
        self,
        divisor: f64,
        offset: f64,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        let mut sum = Sum::default();
        sum.add(offset)?;
        sum.add(checked(
            self.main / divisor,
            "reservoir temperature increment",
        )?)?;
        sum.add(checked(
            self.correction / divisor,
            "reservoir temperature correction",
        )?)?;
        sum.value()
    }
    /// Required prefix heat relative to independently evaluated bath enthalpy.
    fn offset(self, baseline: f64) -> Result<f64, SolubilityCaloricErrorV1> {
        self.add(-baseline)?.value()
    }
    pub(super) fn residual(self, enthalpy: f64) -> Result<f64, SolubilityCaloricErrorV1> {
        checked(
            -self.offset(enthalpy)?,
            "reservoir whole-window heat residual",
        )
    }
}

#[derive(Clone, Copy)]
pub(super) enum StageHeat {
    Prescribed(f64),
    Reservoir {
        source: Reservoir,
        prefix: HeatReceipt,
    },
    Finite {
        source: FiniteReservoir,
        prefix: HeatReceipt,
    },
}

pub(super) struct StageResult {
    pub(super) state: WorkingState,
    pub(super) heat: HeatReceipt,
}
impl StageHeat {
    pub(super) fn required(
        self,
        baseline: f64,
        duration: Time,
        endpoint: Temperature,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        match self {
            Self::Prescribed(target) => checked(target - baseline, "thermal kinetic stage heat"),
            Self::Reservoir { source, prefix } => {
                let mut sum = Sum::default();
                sum.add(prefix.offset(baseline)?)?;
                sum.add(source.heat(duration, endpoint)?)?;
                sum.value()
            }
            Self::Finite { source, prefix } => {
                let mut sum = Sum::default();
                sum.add(prefix.offset(baseline)?)?;
                sum.add(source.heat(prefix, duration, endpoint)?)?;
                sum.value()
            }
        }
    }

    pub(super) fn finish(
        self,
        state: WorkingState,
        duration: Time,
    ) -> Result<StageResult, SolubilityCaloricErrorV1> {
        let heat = match self {
            Self::Prescribed(target) => HeatReceipt::from_value(target),
            Self::Reservoir { source, prefix } => {
                prefix.add(source.heat(duration, state.temperature)?)?
            }
            Self::Finite { source, prefix } => {
                prefix.add(source.heat(prefix, duration, state.temperature)?)?
            }
        };
        Ok(StageResult { state, heat })
    }
}

#[cfg(test)]
mod tests;
