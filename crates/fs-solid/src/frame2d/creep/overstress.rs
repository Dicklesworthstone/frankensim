//! Exact native tabulated Perzyna work on the existing physical-time frame owner.
//!
//! Original table storage, yield heat, physical overstress work and numerical
//! integration loss remain separate. This module adapts the native account; it
//! owns no additional equilibrium solve, return map or internal history.

use super::{FrameCreepError, require};
use fs_material::creep::{CreepEnergy, CreepTolerances};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::UniaxialJ2State;
use fs_material::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;

/// Per-reference-volume native work for timed frame materials [J/m3].
/// Existing creep and yield-offset fields retain their original meanings;
/// physical Perzyna work has its own channel and contributes to heat once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameInelasticEnergy {
    pub initial_stored_j_m3: f64,
    pub final_stored_j_m3: f64,
    pub mechanical_work_j_m3: f64,
    /// Sum of creep, yield-offset and physical overstress dissipation.
    pub physical_dissipation_j_m3: f64,
    pub creep_dissipation_j_m3: f64,
    /// Yield-offset work only; hardening is retained as storage.
    pub plastic_dissipation_j_m3: f64,
    /// Endpoint viscous overstress times accumulated plastic increment.
    pub overstress_dissipation_j_m3: f64,
    /// Native endpoint integration loss, never a caloric heat source.
    pub numerical_dissipation_j_m3: f64,
    pub thermoelastic_storage_change_j_m3: f64,
    pub reversible_heat_j_m3: f64,
    pub residual_j_m3: f64,
}

impl FrameInelasticEnergy {
    pub fn caloric_heat_j_m3(self) -> f64 {
        self.physical_dissipation_j_m3 + self.reversible_heat_j_m3
    }
}

impl From<CreepEnergy> for FrameInelasticEnergy {
    fn from(value: CreepEnergy) -> Self {
        Self {
            initial_stored_j_m3: value.initial_stored_j_m3,
            final_stored_j_m3: value.final_stored_j_m3,
            mechanical_work_j_m3: value.mechanical_work_j_m3,
            physical_dissipation_j_m3: value.physical_dissipation_j_m3,
            creep_dissipation_j_m3: value.creep_dissipation_j_m3,
            plastic_dissipation_j_m3: value.plastic_dissipation_j_m3,
            overstress_dissipation_j_m3: 0.0,
            numerical_dissipation_j_m3: value.numerical_dissipation_j_m3,
            thermoelastic_storage_change_j_m3: value.thermoelastic_storage_change_j_m3,
            reversible_heat_j_m3: value.reversible_heat_j_m3,
            residual_j_m3: value.residual_j_m3,
        }
    }
}

/// Changing the actual interval preserves every authored coefficient or selected
/// original source. A source-backed curve cannot be replaced by an equal-valued
/// authored curve or by a different card, claim, kinetic selection or context.
pub(super) fn admit_rebind(
    old: &UniaxialTabulatedOverstress,
    next: &UniaxialTabulatedOverstress,
) -> Result<(), FrameCreepError> {
    let a = old.backbone();
    let b = next.backbone();
    let ta = a.backbone().hardening();
    let tb = b.backbone().hardening();
    require(
        a.backbone().backbone() == b.backbone().backbone()
            && (std::ptr::eq(ta.points(), tb.points()) || ta == tb)
            && a.viscosity_pa_s() == b.viscosity_pa_s()
            && a.exponent() == b.exponent()
            && a.stress_scale_pa() == b.stress_scale_pa(),
        "tabulated overstress history cannot change its original elasticity, table, Prager storage or rate parameters",
    )?;
    match (old.source(), next.source()) {
        (None, None) => Ok(()),
        (Some(a), Some(b)) => require(
            a.kinetics() == b.kinetics()
                && (std::ptr::eq(a.source(), b.source()) || a.source() == b.source()),
            "tabulated overstress refresh must retain the complete original source and fixed context",
        ),
        _ => Err(FrameCreepError::Invalid(
            "authored and sourced tabulated overstress histories cannot be interchanged",
        )),
    }
}

pub(super) fn audit(
    law: &UniaxialTabulatedOverstress,
    state: &UniaxialJ2State,
    accepted: &UniaxialJ2State,
    old_strain: f64,
    new_strain: f64,
    tolerances: CreepTolerances,
) -> Result<FrameInelasticEnergy, FrameCreepError> {
    let increment = law
        .advance_with_energy(
            old_strain,
            state,
            new_strain,
            EnergyTolerances {
                yield_pa: tolerances.stress_pa,
                energy_j_m3: tolerances.energy_j_m3,
            },
        )
        .map_err(FrameCreepError::Plastic)?;
    let history = accepted.tensor_history();
    let plasticity = law.backbone().backbone().backbone();
    let modulus = plasticity.elastic.youngs;
    require(
        history
            .plastic_strain
            .iter()
            .zip(increment.state.plastic_strain)
            .all(|(a, b)| (a - b).abs() * modulus <= tolerances.stress_pa)
            && (history.alpha - increment.state.alpha).abs() * modulus <= tolerances.stress_pa,
        "published axial plastic history disagrees with its native physical-time work account",
    )?;
    let delta = increment.state.alpha - state.equivalent_plastic_strain();
    let viscous = law
        .backbone()
        .endpoint_overstress_pa(delta)
        .map_err(|error| {
            FrameCreepError::Mechanical(super::SolidError::MaterialRefused {
                what: error.to_string(),
            })
        })?;
    let energy = increment.energy;
    Ok(FrameInelasticEnergy {
        initial_stored_j_m3: energy.initial.total_j_m3(),
        final_stored_j_m3: energy.final_state.total_j_m3(),
        mechanical_work_j_m3: energy.mechanical_work_j_m3,
        physical_dissipation_j_m3: energy.physical_dissipation_j_m3,
        creep_dissipation_j_m3: 0.0,
        plastic_dissipation_j_m3: plasticity.yield_stress * delta,
        overstress_dissipation_j_m3: viscous * delta,
        numerical_dissipation_j_m3: energy.numerical_dissipation_j_m3,
        thermoelastic_storage_change_j_m3: 0.0,
        reversible_heat_j_m3: 0.0,
        residual_j_m3: energy.residual_j_m3,
    })
}
