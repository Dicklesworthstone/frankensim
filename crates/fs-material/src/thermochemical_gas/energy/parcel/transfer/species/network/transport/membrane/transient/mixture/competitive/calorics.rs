//! Complete isothermal membrane calorics in the original gas reference.
//!
//! A new, explicit property on each ORIGINAL species/layer card supplies
//! `u_mobile - h_gas(T)` in J/mol. Neither Henry solubility nor gas Cp supplies
//! this missing condensed-phase caloric datum. Original Langmuir calorics own
//! `u_bound - h_gas(T) = -adsorption_heat` in their rigid-solid convention.
//! At held temperature the dry skeleton has no sensible storage change.
//!
//! Conserved species enthalpy datums cancel symbolically before arithmetic.
//! The remaining membrane storage is the heat its thermostat must supply;
//! the original finite gas thermostats separately remove pressure-flow work.
//! Actual original-to-final gas inventories close the combined first law,
//! while the original mass/species budgets remain independently mandatory.
//! No large formation datum is multiplied by numerical inventory residuals.
//!
//! Layer rows are storage decompositions, not local heat-flux predictions:
//! interlayer thermal conduction and interfacial heat partition are unmodeled.
//! This is a fixed-temperature, rigid-geometry closure, not a temperature
//! evolution, entropy-production, excess-mixture caloric or heat-capacity model.

use super::super::super::{admit_henry, slab::adaptive::MembraneDiffusionAdvanceWorkspace};
use super::super::adaptive::MembraneMixtureAdvanceControls;
use super::adaptive::{
    IsothermalCompetitiveGasAdvanceError, IsothermalCompetitiveGasAdvanceReport,
};
use super::*;
use crate::state_point::{MaterialStatePointError, ScalarAdmissibility, ScalarPropertyRequirement};
use fs_qty::Dims;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Signed condensed mobile internal energy minus the ORIGINAL species gas
/// molar enthalpy at the same temperature. Negative values release heat on
/// uptake. The per-species transport binding fixes the species identity.
pub const MOBILE_RELATIVE_INTERNAL_ENERGY_PROPERTY: &str =
    "gas_membrane_mobile_internal_energy_relative_to_gas_enthalpy";
pub const MOLAR_INTERNAL_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, -1]);

#[derive(Clone, Debug, PartialEq)]
pub enum MembraneCaloricError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Native(IsothermalCompetitiveGasError),
    Adaptive(IsothermalCompetitiveGasAdvanceError),
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for MembraneCaloricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound isothermal membrane calorics: {self:?}")
    }
}
impl std::error::Error for MembraneCaloricError {}
impl From<MaterialStatePointError> for MembraneCaloricError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Source(e)
    }
}
impl From<IsothermalCompetitiveGasError> for MembraneCaloricError {
    fn from(e: IsothermalCompetitiveGasError) -> Self {
        if e == IsothermalCompetitiveGasError::Cancelled {
            Self::Cancelled
        } else {
            Self::Native(e)
        }
    }
}
impl From<IsothermalGasMembraneError> for MembraneCaloricError {
    fn from(e: IsothermalGasMembraneError) -> Self {
        IsothermalCompetitiveGasError::from(e).into()
    }
}
impl From<GasEnergyError> for MembraneCaloricError {
    fn from(e: GasEnergyError) -> Self {
        IsothermalCompetitiveGasError::from(e).into()
    }
}
impl From<IsothermalCompetitiveGasAdvanceError> for MembraneCaloricError {
    fn from(e: IsothermalCompetitiveGasAdvanceError) -> Self {
        if e == IsothermalCompetitiveGasAdvanceError::Cancelled {
            Self::Cancelled
        } else {
            Self::Adaptive(e)
        }
    }
}
type Result<T> = core::result::Result<T, MembraneCaloricError>;

#[derive(Clone, Debug, PartialEq)]
pub struct MembraneLayerCaloricAccount {
    pub layer: usize,
    pub mobile_amount_change_mol: f64,
    pub trapped_amount_change_mol: f64,
    pub mobile_relative_internal_energy_j_mol: f64,
    /// None is an absent trapping law, never an assumed zero bound energy.
    pub trapped_relative_internal_energy_j_mol: Option<f64>,
    pub mobile_relative_storage_change_j: f64,
    pub trapped_relative_storage_change_j: f64,
    /// q_mobile * (delta n_mobile + delta n_bound). Includes redistribution
    /// between layers with different mobile calorics; it is not a local flux.
    pub dissolution_contribution_j: f64,
    /// (q_bound - q_mobile) * delta n_bound, so binding energy is counted once.
    pub trapping_contribution_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneSpeciesCaloricAccount {
    pub species: SpeciesId,
    pub layers: Vec<MembraneLayerCaloricAccount>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneCaloricAccount {
    pub interval_s: [f64; 2],
    pub temperature_k: f64,
    pub species: Vec<MembraneSpeciesCaloricAccount>,
    /// Original gas h_i(T) is removed once per conserved species inventory.
    pub membrane_relative_storage_change_j: f64,
    pub gas_relative_storage_changes_j: [f64; 2],
    pub dissolution_contribution_j: f64,
    pub trapping_contribution_j: f64,
    /// Positive heat enters the membrane; adsorption commonly makes it negative.
    pub thermostat_heat_into_membrane_j: f64,
    pub thermostat_heat_into_gas_j: [f64; 2],
    pub total_thermostat_heat_into_system_j: f64,
    pub total_relative_storage_change_j: f64,
    /// Computed in small reference-relative accounts, not by subtracting
    /// displayed absolute gas enthalpy or nearly equal large aggregate totals.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct IsothermalCaloricCompetitiveGasMembrane<'a> {
    state: IsothermalCompetitiveGasMembrane<'a>,
    mobile_relative: Arc<[Vec<f64>]>,
    energy_tolerance_j: f64,
}
#[derive(Clone, Debug)]
pub struct CaloricCompetitiveGasCandidate<'a, T> {
    state: IsothermalCaloricCompetitiveGasMembrane<'a>,
    transport: T,
    calorics: MembraneCaloricAccount,
}
impl<'a, T> CaloricCompetitiveGasCandidate<'a, T> {
    pub const fn state(&self) -> &IsothermalCaloricCompetitiveGasMembrane<'a> {
        &self.state
    }
    pub const fn transport(&self) -> &T {
        &self.transport
    }
    pub const fn calorics(&self) -> &MembraneCaloricAccount {
        &self.calorics
    }
    pub fn into_parts(
        self,
    ) -> (
        IsothermalCaloricCompetitiveGasMembrane<'a>,
        T,
        MembraneCaloricAccount,
    ) {
        (self.state, self.transport, self.calorics)
    }
}

impl<'a> IsothermalCompetitiveGasMembrane<'a> {
    /// Resolve the explicit mobile caloric property from every original layer
    /// card and selection. The added requirement joins the original binding,
    /// so EVERY fixed or adaptive native probe rechecks its actual pressure
    /// envelope with the same pinned source claim. No detached endpoint-only
    /// source query or assumed zero energy is used.
    ///
    /// `energy_tolerance_j` is an absolute first-law closure bound, distinct
    /// from the original controller's temporal amount/pressure tolerances.
    pub fn with_caloric_closure(
        mut self,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalCaloricCompetitiveGasMembrane<'a>> {
        poll(&mut cancelled)?;
        if !energy_tolerance_j.is_finite() || energy_tolerance_j <= 0.0 {
            return Err(MembraneCaloricError::Invalid(
                "positive finite whole-system energy tolerance required",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            MOBILE_RELATIVE_INTERNAL_ENERGY_PROPERTY,
            MOLAR_INTERNAL_ENERGY_DIMS,
            ScalarAdmissibility::Finite,
        )?;
        let mut channels = Vec::with_capacity(self.channels.len());
        let mut values = Vec::with_capacity(self.channels.len());
        let mut support = Vec::with_capacity(self.support.len());
        for (s, channel) in self.channels.iter().enumerate() {
            poll(&mut cancelled)?;
            let mut channel = channel.clone();
            let mut bindings = (*channel.bindings).clone();
            let mut row = Vec::with_capacity(bindings.layers.len());
            let mut receipts = self.support[s].clone();
            for (l, layer) in bindings.layers.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                layer.original = layer.original.with_additional_requirements(
                    layer.binding.material,
                    core::slice::from_ref(&requirement),
                    layer.binding.selection.clone(),
                )?;
                let property = layer
                    .original
                    .property(MOBILE_RELATIVE_INTERNAL_ENERGY_PROPERTY)
                    .ok_or(MembraneCaloricError::Invalid(
                        "original mobile caloric source missing",
                    ))?;
                admit_henry(property)?;
                admit_independent(property)?;
                row.push(property.value_si());
                for endpoint in &mut receipts.layers[l] {
                    *endpoint = layer
                        .original
                        .at_query_point(layer.binding.material, &endpoint.query())?;
                }
            }
            channel.bindings = Arc::new(bindings);
            channels.push(channel);
            values.push(row);
            support.push(receipts);
        }
        self.channels = channels.into();
        self.support = support;
        poll(&mut cancelled)?;
        Ok(IsothermalCaloricCompetitiveGasMembrane {
            state: self,
            mobile_relative: values.into(),
            energy_tolerance_j,
        })
    }
}

impl<'a> IsothermalCaloricCompetitiveGasMembrane<'a> {
    pub const fn state(&self) -> &IsothermalCompetitiveGasMembrane<'a> {
        &self.state
    }
    pub const fn energy_tolerance_j(&self) -> f64 {
        self.energy_tolerance_j
    }
    pub fn workspace(&self, maximum_bytes: usize) -> Result<CompetitiveMembraneWorkspace> {
        Ok(self.state.workspace(maximum_bytes)?)
    }
    pub fn adaptive_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>> {
        Ok(self.state.adaptive_workspace(maximum_bytes)?)
    }
    pub fn step_candidate(
        &self,
        duration_s: f64,
        work: &mut CompetitiveMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CaloricCompetitiveGasCandidate<'a, IsothermalCompetitiveGasStep>> {
        let (next, transport) = self
            .state
            .step_candidate(duration_s, work, &mut cancelled)?
            .into_parts();
        let calorics = self.measure(&next, transport.thermostat_heat_into_gas_j, &mut cancelled)?;
        Ok(self.finish(next, transport, calorics))
    }
    /// The original competitive controller advances the full source-qualified
    /// gas and site system. Its accepted-only gas heat is combined with direct
    /// original-to-final cell storage, which telescopes exactly at fixed T.
    /// Temporal accuracy remains controlled by the original per-species state
    /// and face estimates, not an independently claimed heat-error certificate.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: MembraneMixtureAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CaloricCompetitiveGasCandidate<'a, IsothermalCompetitiveGasAdvanceReport>> {
        let (next, transport) = self
            .state
            .advance_candidate(duration_s, controls, work, &mut cancelled)?
            .into_parts();
        let calorics = self.measure(&next, transport.thermostat_heat_into_gas_j, &mut cancelled)?;
        Ok(self.finish(next, transport, calorics))
    }
    fn finish<T>(
        &self,
        next: IsothermalCompetitiveGasMembrane<'a>,
        transport: T,
        calorics: MembraneCaloricAccount,
    ) -> CaloricCompetitiveGasCandidate<'a, T> {
        CaloricCompetitiveGasCandidate {
            state: Self {
                state: next,
                mobile_relative: self.mobile_relative.clone(),
                energy_tolerance_j: self.energy_tolerance_j,
            },
            transport,
            calorics,
        }
    }
    fn measure(
        &self,
        next: &IsothermalCompetitiveGasMembrane<'_>,
        gas_heat: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MembraneCaloricAccount> {
        poll(cancelled)?;
        let old = &self.state;
        let temperature = old.temperature_k();
        let before_gas = [
            inventory(old.gases[0].model(), old.gases[0].mass_kg())?,
            inventory(old.gases[1].model(), old.gases[1].mass_kg())?,
        ];
        let after_gas = [
            inventory(next.gases[0].model(), next.gases[0].mass_kg())?,
            inventory(next.gases[1].model(), next.gases[1].mass_kg())?,
        ];
        let mut species = Vec::with_capacity(old.channels.len());
        let mut membrane_storage = Sum::default();
        let mut dissolution = Sum::default();
        let mut trapping = Sum::default();
        let mut gas_storage = [Sum::default(), Sum::default()];
        for (s, channel) in old.channels.iter().enumerate() {
            poll(cancelled)?;
            let id = &channel.bindings.species;
            for side in 0..2 {
                let dm = finite(
                    after_gas[side].get(id).copied().unwrap_or(0.0)
                        - before_gas[side].get(id).copied().unwrap_or(0.0),
                )?;
                gas_storage[side].add(finite(-dm * (R * temperature / channel.molar_mass))?)?;
            }
            let a = old.membrane().mobile(s).expect("admitted species");
            let b = next.membrane().mobile(s).expect("admitted species");
            let ta = old.membrane().trapped_moles(s).expect("admitted species");
            let tb = next.membrane().trapped_moles(s).expect("admitted species");
            let mut layers = Vec::with_capacity(channel.bindings.layers.len());
            for (l, range) in a.model().layer_cell_ranges().iter().enumerate() {
                poll(cancelled)?;
                let source = &a.model().layers()[l];
                let capacity = finite(
                    (source.solubility.value_si() * (source.thickness_m / source.cells as f64))
                        * a.model().area_m2(),
                )?;
                let mut mobile = Sum::default();
                let mut bound = Sum::default();
                for cell in range.clone() {
                    if cell % 128 == 0 {
                        poll(cancelled)?;
                    }
                    mobile.add(finite(
                        capacity
                            * (b.partial_pressures_pa()[cell] - a.partial_pressures_pa()[cell]),
                    )?)?;
                    bound.add(finite(tb[cell] - ta[cell])?)?;
                }
                let qm = self.mobile_relative[s][l];
                let qb = old.trap_support[s][l]
                    .as_ref()
                    .map(|receipt| -receipt[0].adsorption_heat_j_mol);
                if qb.is_none() && bound.value != 0.0 {
                    return Err(MembraneCaloricError::Invalid(
                        "bound storage changed without its original caloric source",
                    ));
                }
                let mobile_energy = finite(qm * mobile.value)?;
                let bound_energy = finite(qb.unwrap_or(0.0) * bound.value)?;
                let dissolved = finite(qm * finite(mobile.value + bound.value)?)?;
                let trapped = finite((qb.unwrap_or(qm) - qm) * bound.value)?;
                membrane_storage.add(mobile_energy)?;
                membrane_storage.add(bound_energy)?;
                dissolution.add(dissolved)?;
                trapping.add(trapped)?;
                layers.push(MembraneLayerCaloricAccount {
                    layer: l,
                    mobile_amount_change_mol: mobile.value,
                    trapped_amount_change_mol: bound.value,
                    mobile_relative_internal_energy_j_mol: qm,
                    trapped_relative_internal_energy_j_mol: qb,
                    mobile_relative_storage_change_j: mobile_energy,
                    trapped_relative_storage_change_j: bound_energy,
                    dissolution_contribution_j: dissolved,
                    trapping_contribution_j: trapped,
                });
            }
            species.push(MembraneSpeciesCaloricAccount {
                species: id.clone(),
                layers,
            });
        }
        let membrane_heat = finite(dissolution.value + trapping.value)?;
        let mut residual = Sum::default();
        residual.add(membrane_storage.value)?;
        residual.add(-dissolution.value)?;
        residual.add(-trapping.value)?;
        let mut total_storage = Sum::default();
        total_storage.add(membrane_storage.value)?;
        let mut total_heat = Sum::default();
        total_heat.add(membrane_heat)?;
        for side in 0..2 {
            residual.add(finite(gas_storage[side].value - gas_heat[side])?)?;
            total_storage.add(gas_storage[side].value)?;
            total_heat.add(gas_heat[side])?;
        }
        if residual.value.abs() > self.energy_tolerance_j {
            return Err(MembraneCaloricError::EnergyBalance {
                residual_j: residual.value,
                tolerance_j: self.energy_tolerance_j,
            });
        }
        poll(cancelled)?;
        Ok(MembraneCaloricAccount {
            interval_s: [old.elapsed_s(), next.elapsed_s()],
            temperature_k: temperature,
            species,
            membrane_relative_storage_change_j: membrane_storage.value,
            gas_relative_storage_changes_j: gas_storage.map(|s| s.value),
            dissolution_contribution_j: dissolution.value,
            trapping_contribution_j: trapping.value,
            thermostat_heat_into_membrane_j: membrane_heat,
            thermostat_heat_into_gas_j: gas_heat,
            total_thermostat_heat_into_system_j: total_heat.value,
            total_relative_storage_change_j: total_storage.value,
            energy_residual_j: residual.value,
        })
    }
}
