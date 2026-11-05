//! Temperature-dependent brittle strength on the original retained solid.
//!
//! The ordered physical split changes sourced coefficients at the accepted
//! strain, then re-equilibrates using the original P1 Newton owner at the new
//! temperature. Parameter work and induced damage remain separate from the
//! subsequent mechanical work, fracture loss and endpoint integration defect.
//! No intermediate source/state change is published if either stage refuses.
//!
//! Temperature changes coefficients only: no thermal eigenstrain, caloric
//! equation, heat conversion or entropy-production claim is added. The caller
//! resolves this split with physical temperature/load increments.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::{BrittleSourceError, ResolvedBrittleCrackBand};

use super::{
    BrittleTetControls, BrittleTetError, BrittleTetSolid, BrittleTetStep, checkpoint, energy_gate,
};
use crate::plastic3::{IncrementSettings, SmallStrainTetSolution};

/// Original source-transfer accounts per reference volume, J/m3.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetTemperatureElementEnergy {
    /// Signed work from changing material coefficients at fixed strain.
    pub parameter_work_j_m3: f64,
    /// Additional irreversible fracture loss during the parameter change.
    pub damage_dissipated_j_m3: f64,
    pub stored_change_j_m3: f64,
    pub defect_j_m3: f64,
}

/// Whole-solid parameter-transfer accounts, before mechanical equilibration.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BrittleTetTemperatureEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed coefficient-change work; this is not an inferred heat input.
    pub parameter_work_j: f64,
    pub fracture_dissipation_j: f64,
    /// Parameter work minus storage change and induced fracture loss.
    pub balance_error_j: f64,
}

/// A complete accepted source transition and mechanical endpoint.
#[derive(Debug)]
pub struct BrittleTetTemperatureStep<'a> {
    /// Final solid and the separate original mechanical endpoint accounts.
    pub mechanical: BrittleTetStep<'a>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub transfer: BrittleTetTemperatureEnergy,
    /// Combined parameter plus exact mechanical work minus total storage
    /// change and fracture loss. External endpoint quadrature is separate.
    pub total_balance_error_j: f64,
}

/// Original source-transfer refusal or original solid/energy/cancellation gate.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetTemperatureError {
    Source {
        element: usize,
        source: BrittleSourceError,
    },
    Solid(BrittleTetError),
}

impl core::fmt::Display for BrittleTetTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle solid temperature endpoint refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetTemperatureError {}
impl From<BrittleTetError> for BrittleTetTemperatureError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}

impl<'a> BrittleTetSolid<'a> {
    /// Transfer the original histories to checked successor source states,
    /// then solve one actual mechanical endpoint at those new coefficients.
    ///
    /// Successors must preserve the original card/claims, nonthermal state,
    /// Poisson ratio, fixed plane, band width and constitutive bounds under the
    /// material owner's rules. They must outlive the returned solid. Loads are
    /// absolute targets; supports may be held fixed or allowed to move through
    /// the same native equilibrium path. Every refusal leaves `self` unchanged.
    ///
    /// Parameter work and cumulative loss obey the original material-energy
    /// tolerance; mechanical gates retain their original meanings. A loss of
    /// static load capacity is a refusal, not an invented postfailure equilibrium.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_at_temperature<'next>(
        &self,
        next_materials: &'next [ResolvedBrittleCrackBand],
        temperature_axis: &str,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetTemperatureStep<'next>, BrittleTetTemperatureError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if next_materials.len() != self.materials.len() {
            return Err(BrittleTetError::MaterialCount {
                expected: self.materials.len(),
                got: next_materials.len(),
            }
            .into());
        }
        if [
            controls.material_energy_tolerance_j,
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative work tolerances required").into(),
            );
        }
        let mut states = Vec::with_capacity(next_materials.len());
        let mut stresses = Vec::with_capacity(next_materials.len());
        let mut element_transfer = Vec::with_capacity(next_materials.len());
        let mut transfer = BrittleTetTemperatureEnergy {
            stored_before_j: self.stored_energy_j,
            ..BrittleTetTemperatureEnergy::default()
        };
        let mut cumulative_loss = 0.0;
        for (element, (old_source, next_source)) in
            self.materials.iter().zip(next_materials).enumerate()
        {
            checkpoint(cx)?;
            let old_state = &self.solution.element_states[element];
            let step = old_source
                .transfer_temperature(
                    old_state,
                    next_source,
                    temperature_axis,
                    next_source.law().settings().energy_tolerance_j_m3,
                )
                .map_err(|source| BrittleTetTemperatureError::Source { element, source })?;
            let held = next_source
                .law()
                .trial(step.state.strain(), &step.state)
                .map_err(|source| BrittleTetTemperatureError::Source {
                    element,
                    source: source.into(),
                })?;
            if held.state != step.state || held.energy.work != 0.0 || held.energy.dissipated != 0.0
            {
                return Err(BrittleTetError::Invalid(
                    "transferred native state is not stationary at its held strain",
                )
                .into());
            }
            if !step.state.is_valid()
                || step.state.strain() != old_state.strain()
                || step.state.damage() < old_state.damage()
                || step.state.dissipated_j_m3() < old_state.dissipated_j_m3()
                || step.damage_dissipated_j_m3 < 0.0
            {
                return Err(BrittleTetError::Invalid(
                    "source transfer changed strain or healed history",
                )
                .into());
            }
            let volume = self.bands[element].reference_volume_m3();
            transfer.stored_after_j =
                volume.mul_add(step.state.stored_energy_j_m3(), transfer.stored_after_j);
            transfer.parameter_work_j =
                volume.mul_add(step.parameter_work_j_m3, transfer.parameter_work_j);
            transfer.fracture_dissipation_j =
                volume.mul_add(step.damage_dissipated_j_m3, transfer.fracture_dissipation_j);
            cumulative_loss = volume.mul_add(step.state.dissipated_j_m3(), cumulative_loss);
            element_transfer.push(BrittleTetTemperatureElementEnergy {
                parameter_work_j_m3: step.parameter_work_j_m3,
                damage_dissipated_j_m3: step.damage_dissipated_j_m3,
                stored_change_j_m3: step.stored_change_j_m3,
                defect_j_m3: step.defect_j_m3,
            });
            states.push(step.state);
            stresses.push(held.stress_pa);
        }
        transfer.balance_error_j = (transfer.parameter_work_j
            - (transfer.stored_after_j - transfer.stored_before_j))
            - transfer.fracture_dissipation_j;
        energy_gate(
            "temperature parameter work/storage/fracture",
            transfer.balance_error_j,
            controls.material_energy_tolerance_j,
        )?;
        energy_gate(
            "temperature cumulative fracture loss",
            (cumulative_loss - self.fracture_dissipation_j) - transfer.fracture_dissipation_j,
            controls.material_energy_tolerance_j,
        )?;
        if [
            transfer.stored_after_j,
            transfer.fracture_dissipation_j,
            cumulative_loss,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                BrittleTetError::Invalid("finite nonnegative transferred energy required").into(),
            );
        }
        // This possibly unequilibrated state is private. Only a successful
        // native mechanical endpoint can publish its successor material/history.
        let candidate = BrittleTetSolid {
            geometry: self.geometry,
            materials: next_materials,
            bands: Arc::clone(&self.bands),
            budget: self.budget,
            solution: SmallStrainTetSolution {
                displacement_m: self.solution.displacement_m.clone(),
                element_strain: self.solution.element_strain.clone(),
                element_stress: stresses,
                element_states: states,
                increments: Vec::new(),
            },
            stored_energy_j: transfer.stored_after_j,
            fracture_dissipation_j: cumulative_loss,
        };
        let mechanical = candidate.advance(prescribed_m, nodal_forces_n, settings, controls, cx)?;
        let total_balance_error_j = ((transfer.parameter_work_j
            + mechanical.energy.material_work_j)
            - (mechanical.solid.stored_energy_j() - transfer.stored_before_j))
            - (transfer.fracture_dissipation_j + mechanical.energy.fracture_dissipation_j);
        energy_gate(
            "combined temperature/mechanical work",
            total_balance_error_j,
            controls.material_energy_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(BrittleTetTemperatureStep {
            mechanical,
            element_transfer,
            transfer,
            total_balance_error_j,
        })
    }
}
