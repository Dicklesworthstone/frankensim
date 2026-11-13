//! Original physical supplied/resistive power and complete state covectors.

use super::{
    ForcedMeanObjective, ForcedMeanObjectiveControls, ForcedSensitivityError,
    PhsForcedSensitivityError, PhsOrbitError, PhsStorageHessianAction,
    PreparedPhsForcedSensitivity, bad,
};
use crate::phs_orbit::tangent::StorageTangent;
use fs_math::dd::Dd;
use fs_phs::PortHamiltonian;

mod adjoint;
pub use adjoint::PreparedPhsForcedPowerAdjoint;

/// Original physical power, averaged over a cycle or integrated over its
/// actual period. Supplied work and resistance loss remain separate signed
/// observations; neither is reconstructed from the other's balance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhsPowerObjective {
    MeanSuppliedPower,
    MeanResistancePower,
    SuppliedEnergyPerCycle,
    ResistanceEnergyPerCycle,
}

impl PhsPowerObjective {
    pub(in crate::phs_orbit) fn supplied(self) -> bool {
        matches!(self, Self::MeanSuppliedPower | Self::SuppliedEnergyPerCycle)
    }

    pub(in crate::phs_orbit) fn per_cycle(self) -> bool {
        matches!(
            self,
            Self::SuppliedEnergyPerCycle | Self::ResistanceEnergyPerCycle
        )
    }

    // Sample the normalized objective itself, so value AND every covector
    // component are accepted in the requested power/energy units on both grids.
    pub(in crate::phs_orbit) fn normalize(
        self,
        omega: f64,
        value: &mut f64,
        gradient: &mut [f64],
        frequency: Option<&mut f64>,
    ) -> bool {
        if self.per_cycle() {
            let period = core::f64::consts::TAU / omega;
            if !period.is_finite() || period <= 0.0 {
                return false;
            }
            for entry in core::iter::once(&mut *value).chain(gradient.iter_mut()) {
                let scaled = (Dd::from_f64(period) * Dd::from_f64(*entry)).to_f64();
                if !scaled.is_finite() || (scaled == 0.0 && *entry != 0.0) {
                    return false;
                }
                *entry = scaled;
            }
        }
        if let Some(frequency) = frequency {
            *frequency = if self.per_cycle() {
                -*value / omega
            } else {
                0.0
            };
            if !frequency.is_finite() || (*frequency == 0.0 && self.per_cycle() && *value != 0.0) {
                return false;
            }
        }
        value.is_finite() && gradient.iter().all(|v| v.is_finite())
    }

    pub(in crate::phs_orbit) fn direct_partial(
        self,
        omega: f64,
        omega_derivative: f64,
        direct_power: f64,
        original_power: f64,
    ) -> Option<f64> {
        if !self.per_cycle() {
            return direct_power.is_finite().then_some(direct_power);
        }
        let period = core::f64::consts::TAU / omega;
        if !period.is_finite() || period <= 0.0 {
            return None;
        }
        let mut derivative = Dd::from_f64(direct_power);
        if omega_derivative != 0.0 {
            derivative = derivative
                - Dd::from_f64(omega_derivative) * Dd::from_f64(original_power)
                    / Dd::from_f64(omega);
        }
        let value = (Dd::from_f64(period) * derivative).to_f64();
        (value.is_finite() && (value != 0.0 || derivative.to_f64() == 0.0)).then_some(value)
    }
}

impl PreparedPhsForcedSensitivity<'_> {
    /// Cycle mean of the ORIGINAL `grad(H)^T R grad(H)` and its complete state
    /// covector, including nonlinear storage and explicit memory. The stored
    /// source's exact continuous Hessian is required; no quadratic substitute
    /// or finite-difference fallback is introduced.
    ///
    /// Direct parameter dependence of R or H remains an explicit objective
    /// partial for the adjoint. Signed raw-parts resistance stays signed; this
    /// observation neither certifies passivity nor assigns loss to a heat sink.
    pub fn mean_resistance_power(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        self.mean_physical_power(
            PhsPowerObjective::MeanResistancePower,
            controls,
            None,
            cancelled,
        )
    }

    /// The same original resistance-power observation with an explicitly
    /// authored COMPLETE original Hamiltonian Hessian action.
    pub fn mean_resistance_power_analytic(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        self.mean_physical_power(
            PhsPowerObjective::MeanResistancePower,
            controls,
            Some(hessian),
            cancelled,
        )
    }

    /// Cycle mean of the original supplied power `grad(H)^T G u(t)` and its
    /// complete state covector. Uses the EXACT drive retained by preparation,
    /// with its DC, phases and every admitted harmonic. Positive and negative
    /// port contributions remain signed and are contracted by the pHS owner.
    ///
    /// The state gradient is `Hxx G u`; source/port/material dependence is a
    /// separate explicit adjoint partial. In the original fixed-phase Fourier
    /// chart this observation has no direct clock dependence at fixed source
    /// coefficients. It is not reconstructed from resistance power or energy
    /// balance, so a finite-order orbit's accounting defect remains observable.
    pub fn mean_supplied_power(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        self.mean_physical_power(
            PhsPowerObjective::MeanSuppliedPower,
            controls,
            None,
            cancelled,
        )
    }

    /// Supplied power under the same retained drive with an explicit complete
    /// original Hessian. No source or tangent is inferred from an accepted orbit.
    pub fn mean_supplied_power_analytic(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        self.mean_physical_power(
            PhsPowerObjective::MeanSuppliedPower,
            controls,
            Some(hessian),
            cancelled,
        )
    }

    fn mean_physical_power(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective, PhsForcedSensitivityError> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let mut observation = PowerObservation::new(self.system, hessian)?;
        let mut input = vec![f64::NAN; self.system.port_dim()];
        let drive = &self.input;
        let omega = self.prepared.angular_frequency();
        self.prepared
            .mean_objective(
                controls,
                |time, state, value, gradient| {
                    if !observation.load(state) {
                        return false;
                    }
                    let complete = if objective.supplied() {
                        input.fill(f64::NAN);
                        drive.sample_into(time, &mut input)
                            && observation.supplied(state, &input, value, gradient)
                    } else {
                        observation.resistance(state, value, gradient)
                    };
                    complete && objective.normalize(omega, value, gradient, None)
                },
                cancelled,
            )
            .map_err(Into::into)
    }
}

// Shared by forced and autonomous observations. Original pHS kernels own all
// R/G contractions; the same storage tangent owner used by HB supplies Hxx.
// Scratch is O(n), plus any source-owned quadratic metadata. As with generic
// mean_objective, source callbacks/metadata are outside fs-orbit's numerical
// work/payload budget; evaluations occur only on its two admitted AFT grids.
pub(in crate::phs_orbit) struct PowerObservation<'a> {
    system: &'a PortHamiltonian,
    tangent: StorageTangent<'a>,
    effort: Vec<f64>,
    effort_gradient: Vec<f64>,
}

impl<'a> PowerObservation<'a> {
    pub(in crate::phs_orbit) fn new(
        system: &'a PortHamiltonian,
        hessian: Option<&'a PhsStorageHessianAction<'a>>,
    ) -> Result<Self, PhsOrbitError> {
        let n = system.state_dim();
        let tangent = match hessian {
            Some(hessian) => StorageTangent::hessian(n, hessian),
            None => StorageTangent::from_storage(system),
        };
        if !tangent.available() {
            return Err(bad(
                "physical power requires the complete original continuous storage tangent",
            ));
        }
        Ok(Self {
            system,
            tangent,
            effort: vec![f64::NAN; n],
            effort_gradient: vec![f64::NAN; n],
        })
    }

    pub(in crate::phs_orbit) fn load(&mut self, state: &[f64]) -> bool {
        self.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut self.effort);
        self.effort.iter().all(|v| v.is_finite())
    }

    pub(in crate::phs_orbit) fn flows(&self, out: &mut [f64]) -> bool {
        self.system
            .output_from_effort_into(&self.effort, out)
            .is_ok()
    }

    pub(in crate::phs_orbit) fn resistance(
        &mut self,
        state: &[f64],
        value: &mut f64,
        gradient: &mut [f64],
    ) -> bool {
        let Ok(power) = self.system.resistance_power_from_effort(&self.effort) else {
            return false;
        };
        *value = power;
        self.system
            .resistance_power_effort_gradient_into(&self.effort, &mut self.effort_gradient)
            .is_ok()
            && self
                .tangent
                .apply_hessian(state, &self.effort_gradient, gradient)
    }

    pub(in crate::phs_orbit) fn supplied(
        &mut self,
        state: &[f64],
        input: &[f64],
        value: &mut f64,
        gradient: &mut [f64],
    ) -> bool {
        let Ok(power) = self.system.supplied_power_from_effort(&self.effort, input) else {
            return false;
        };
        *value = power;
        self.system
            .input_flow_into(input, &mut self.effort_gradient)
            .is_ok()
            && self
                .tangent
                .apply_hessian(state, &self.effort_gradient, gradient)
    }
}

#[cfg(test)]
mod tests;
