//! Thermal cyclic branches on the existing physical-time frame solver.
//!
//! The material owner transports each branch's strain-like coordinate through
//! temperature changes. This adapter integrates the native account, retaining
//! separate yield, viscous, recovery and numerical terms.

use super::{Frame2d, FrameCreepError, FrameInelasticEnergy, SectionModel, require};
use crate::fiber::{FiberLaw, FiberState};
use crate::fiber::chaboche::thermal::{ThermalChabocheFiberLaw, ThermalChabocheFiberState};
use fs_material::creep::CreepTolerances;
use fs_material::plastic::energy::{EnergyError, EnergyTolerances};
use fs_material::plastic::thermal::chaboche::axial::{
    UniaxialThermalChaboche, UniaxialThermalChabocheState,
};

impl Frame2d {
    /// Read all accepted thermal backstresses without another constitutive step.
    pub fn thermal_chaboche_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&ThermalChabocheFiberLaw, &ThermalChabocheFiberState)> {
        let SectionModel::Fiber(section) = self.elements.get(element)?.sections.get(section)?
        else {
            return None;
        };
        match (&section.fibers.get(fiber)?.law, section.states.get(fiber)?) {
            (FiberLaw::ThermalChaboche(law), FiberState::ThermalChaboche(state)) => {
                Some((law, state))
            }
            _ => None,
        }
    }
}

pub(super) fn audit(
    law: &ThermalChabocheFiberLaw,
    old: &ThermalChabocheFiberState,
    accepted: &ThermalChabocheFiberState,
    old_strain: f64,
    new_strain: f64,
    tolerances: CreepTolerances,
) -> Result<FrameInelasticEnergy, FrameCreepError> {
    macro_rules! dispatch {
        ($($variant:ident),+ $(,)?) => {
            match (law, old, accepted) {
                $(
                    (
                        ThermalChabocheFiberLaw::$variant(l),
                        ThermalChabocheFiberState::$variant(a),
                        ThermalChabocheFiberState::$variant(b),
                    ) => audit_native(l, a, b, old_strain, new_strain, tolerances),
                )+
                _ => Err(FrameCreepError::Invalid(
                    "thermal cyclic work must retain the complete original branch history",
                )),
            }
        };
    }
    dispatch!(One, Two, Three, Four, Five, Six, Seven, Eight)
}

fn audit_native<const N: usize>(
    law: &UniaxialThermalChaboche<N>,
    old: &UniaxialThermalChabocheState<N>,
    accepted: &UniaxialThermalChabocheState<N>,
    old_strain: f64,
    new_strain: f64,
    tolerances: CreepTolerances,
) -> Result<FrameInelasticEnergy, FrameCreepError> {
    let increment = law
        .advance_with_energy(
            old_strain,
            old,
            new_strain,
            EnergyTolerances {
                yield_pa: tolerances.stress_pa,
                energy_j_m3: tolerances.energy_j_m3,
            },
        )
        .map_err(FrameCreepError::ThermalChaboche)?;
    // The scalar and full-tensor implicit solves may differ by roundoff.
    // Check every independent component in physical stress units, including
    // accumulated flow; equality of the aggregate backstress is insufficient.
    let native_state = increment.state.mechanical.tensor_history();
    let accepted_state = accepted.mechanical.tensor_history();
    let modulus = law.tensor_law().law().backbone().backbone().elastic.youngs;
    let close = |left: f64, right: f64, scale: f64| {
        let residual = (left - right) * scale;
        residual.is_finite() && residual.abs() <= tolerances.stress_pa
    };
    require(
        increment.state.temperature_k() == accepted.temperature_k()
            && close(native_state.plastic.alpha, accepted_state.plastic.alpha, modulus)
            && native_state.plastic.plastic_strain.iter()
                .zip(&accepted_state.plastic.plastic_strain)
                .all(|(&a, &b)| close(a, b, modulus))
            && native_state.backstress_pa.iter().zip(&accepted_state.backstress_pa)
                .all(|(a, b)| a.iter().zip(b).all(|(&a, &b)| close(a, b, 1.0))),
        "accepted thermal cyclic branches disagree with their native interval work account",
    )?;
    let initial_stored_j_m3 = law
        .stored_energy(old_strain, old)
        .map_err(FrameCreepError::ThermalChaboche)?
        .total_j_m3();
    let delta = increment.state.mechanical.equivalent_plastic_strain()
        - old.mechanical.equivalent_plastic_strain();
    let plastic_dissipation_j_m3 =
        law.tensor_law().law().backbone().backbone().yield_stress * delta;
    let thermal = increment.thermal;
    let native = thermal.mechanical.energy;
    let mut energy = FrameInelasticEnergy {
        initial_stored_j_m3,
        final_stored_j_m3: native.final_state.total_j_m3(),
        mechanical_work_j_m3: thermal.mechanical_work_j_m3,
        physical_dissipation_j_m3: native.physical_dissipation_j_m3,
        creep_dissipation_j_m3: 0.0,
        plastic_dissipation_j_m3,
        overstress_dissipation_j_m3: increment.overstress_dissipation_j_m3,
        recovery_dissipation_j_m3: increment.recovery_dissipation_j_m3,
        numerical_dissipation_j_m3: native.numerical_dissipation_j_m3,
        thermoelastic_storage_change_j_m3: thermal.thermoelastic_storage_change_j_m3,
        reversible_heat_j_m3: thermal.reversible_heat_j_m3,
        residual_j_m3: 0.0,
    };
    let partition_residual = energy.physical_dissipation_j_m3
        - energy.plastic_dissipation_j_m3
        - energy.overstress_dissipation_j_m3
        - energy.recovery_dissipation_j_m3;
    energy.residual_j_m3 = energy.mechanical_work_j_m3
        - (energy.final_stored_j_m3 - energy.initial_stored_j_m3)
        - energy.thermoelastic_storage_change_j_m3
        - energy.physical_dissipation_j_m3
        - energy.reversible_heat_j_m3
        - energy.numerical_dissipation_j_m3;
    for residual_j_m3 in [partition_residual, energy.residual_j_m3] {
        if !residual_j_m3.is_finite() || residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(FrameCreepError::Plastic(EnergyError::EnergyResidual {
                residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            }));
        }
    }
    Ok(energy)
}
