//! Homogeneous caloric changes when a finite-volume cell exchanges species.

use super::super::{
    BinaryCaloricErrorV1, BinaryCaloricStateV1, BinaryLiquidCaloricModelV1, Sum, checked, poll,
    positive, temperature,
};
use fs_qty::Temperature;
use std::sync::Arc;

/// Actual original-source changes, including the enthalpy and entropy carried
/// by the changed species inventory. A single open cell's entropy may decrease;
/// a closed spatial owner must audit the sum over all cells and boundary fluxes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidTransportChangeV1 {
    pub enthalpy_change_j: f64,
    pub entropy_change_j_per_k: f64,
}

/// Independent temperature/energy accuracy and bounded homogeneous inverse work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidTransportSolveV1 {
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
}

/// One completed original-source homogeneous inventory/enthalpy update.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryLiquidTransportBalanceV1 {
    final_state: BinaryCaloricStateV1,
    change: BinaryLiquidTransportChangeV1,
    energy_residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
}

impl BinaryLiquidTransportBalanceV1 {
    pub const fn final_state(&self) -> &BinaryCaloricStateV1 {
        &self.final_state
    }
    pub const fn change(&self) -> &BinaryLiquidTransportChangeV1 {
        &self.change
    }
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    pub const fn temperature_bracket(&self) -> [Temperature; 2] {
        self.temperature_bracket
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

type Result<T> = core::result::Result<T, BinaryCaloricErrorV1>;

impl BinaryLiquidCaloricModelV1 {
    /// Original homogeneous H/S changes with differing local inventories.
    /// Reuses direct thermal/mixing increments, carrying each changed inventory
    /// in the original source-relative datum. No phase flash is substituted.
    pub fn homogeneous_transport_change(
        &self,
        initial: &BinaryCaloricStateV1,
        final_state: &BinaryCaloricStateV1,
    ) -> Result<BinaryLiquidTransportChangeV1> {
        self.admit_transport_state(initial)?;
        self.admit_transport_state(final_state)?;
        Ok(BinaryLiquidTransportChangeV1 {
            enthalpy_change_j: self.transport_enthalpy_change(initial, final_state)?,
            entropy_change_j_per_k: self.transport_entropy_change(initial, final_state)?,
        })
    }

    /// Recover homogeneous temperature from actual new inventories and signed
    /// enthalpy transfer using original Cp and the shared bounded scalar inverse.
    /// Failure publishes no state; closed conservation is the caller's account.
    pub fn homogeneous_transport_candidate(
        &self,
        initial: &BinaryCaloricStateV1,
        new_component_moles: [f64; 2],
        enthalpy_change_j: f64,
        controls: BinaryLiquidTransportSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidTransportBalanceV1> {
        poll(&mut cancelled)?;
        self.admit_transport_state(initial)?;
        let bracket = controls.temperature_bracket.map(|t| t.value());
        let domain = self.temperature_domain().map(|t| t.value());
        if !positive(bracket[0])
            || !positive(bracket[1])
            || bracket[0] >= bracket[1]
            || bracket[0] < domain[0]
            || bracket[1] > domain[1]
            || !positive(controls.temperature_tolerance_k)
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
            || !(1..=256).contains(&controls.maximum_iterations)
            || !enthalpy_change_j.is_finite()
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "original caloric bracket, finite transport energy and bounded homogeneous inverse controls",
            ));
        }
        let root = temperature::solve(
            temperature::TemperatureControl {
                bracket,
                temperature_tolerance_k: controls.temperature_tolerance_k,
                energy_tolerance_j: controls.energy_tolerance_j,
                maximum_iterations: controls.maximum_iterations,
            },
            |t| {
                poll(&mut cancelled)?;
                let state = self.homogeneous(Temperature::new(t), new_component_moles)?;
                let change_j = self.transport_enthalpy_change(initial, &state)?;
                Ok::<_, BinaryCaloricErrorV1>(temperature::TemperatureSample {
                    residual_j: checked(
                        change_j - enthalpy_change_j,
                        "homogeneous transport energy residual",
                    )?,
                    state,
                    change_j,
                })
            },
            |state| state.isobaric_heat_capacity_j_per_k().ok(),
        )?;
        let change = BinaryLiquidTransportChangeV1 {
            enthalpy_change_j: root.sample.change_j,
            entropy_change_j_per_k: self.transport_entropy_change(initial, &root.sample.state)?,
        };
        poll(&mut cancelled)?;
        Ok(BinaryLiquidTransportBalanceV1 {
            final_state: root.sample.state,
            change,
            energy_residual_j: root.sample.residual_j,
            temperature_bracket: root.bracket.map(Temperature::new),
            iterations: root.iterations,
        })
    }

    fn admit_transport_state(&self, state: &BinaryCaloricStateV1) -> Result<()> {
        if (!Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner)
            || state.equilibrium.is_some()
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "transport requires homogeneous states of the same original caloric owner",
            ));
        }
        Ok(())
    }

    fn transport_enthalpy_change(
        &self,
        initial: &BinaryCaloricStateV1,
        final_state: &BinaryCaloricStateV1,
    ) -> Result<f64> {
        let mut change = Sum::default();
        change.add(self.enthalpy_increment(&[initial], final_state)?)?;
        for (i, source) in self.heat_capacities().iter().enumerate() {
            let delta = final_state.component_moles[i] - initial.component_moles[i];
            if delta != 0.0 {
                change.add(checked(
                    delta
                        * source.enthalpy_change_j_per_mol(
                            source.reference_temperature(),
                            final_state.temperature,
                        )?,
                    "exchanged original pure enthalpy",
                )?)?;
            }
        }
        change.add(checked(
            (final_state.total_moles - initial.total_moles)
                * final_state.mixing.mixing_enthalpy_j_per_mol(),
            "exchanged mixing enthalpy",
        )?)?;
        change.value()
    }

    fn transport_entropy_change(
        &self,
        initial: &BinaryCaloricStateV1,
        final_state: &BinaryCaloricStateV1,
    ) -> Result<f64> {
        let mut change = Sum::default();
        change.add(self.entropy_increment(&[initial], final_state)?)?;
        for (i, source) in self.heat_capacities().iter().enumerate() {
            let delta = final_state.component_moles[i] - initial.component_moles[i];
            if delta != 0.0 {
                change.add(checked(
                    delta
                        * source.entropy_change_j_per_mol_k(
                            source.reference_temperature(),
                            final_state.temperature,
                        )?,
                    "exchanged original pure entropy",
                )?)?;
            }
        }
        change.add(checked(
            (final_state.total_moles - initial.total_moles)
                * final_state.mixing.mixing_entropy_j_per_mol_k(),
            "exchanged mixing entropy",
        )?)?;
        change.value()
    }
}
