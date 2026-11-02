//! Source-qualified competitive membrane storage between actual finite mixtures.
//!
//! The original independent-mixture binding supplies mobile Henry/Fick laws,
//! gas catalogs and surface identities. Each layer explicitly supplies ONE
//! original competitive lattice. Species partial pressures drive the coupled
//! stationary-matrix solve; summed polymer pressures and original gas EOS
//! totals qualify the material sources. One original reconstruction per gas
//! updates all species together, followed by independent species mass audits.
//!
//! Temperature and geometry remain fixed. Gas stream and thermostat accounts
//! do not include a mobile/trapped caloric store or heat of sorption. This is
//! neither Maxwell-Stefan gas diffusion nor a source of mixed-gas evidence:
//! supplied claims must already apply to the intended mixture, and unowned
//! composition, humidity or concentration-dependent applicability refuses.

/// Complete competitive breakthrough windows with original finite-gas accounts.
pub mod adaptive;

use super::super::{
    ResolvedMembraneDiffusionError,
    reservoir::IsothermalGasMembraneError,
    slab::{
        MembraneDiffusionError,
        competitive::{
            CompetitiveMembraneChannel, CompetitiveMembraneDiffusionState,
            CompetitiveMembraneError, CompetitiveMembraneSolve, CompetitiveMembraneStep,
            CompetitiveMembraneWorkspace, FiniteCompetitiveMembraneDiffusionState,
        },
        trapping::LangmuirMembraneTrap,
    },
    trapping::{LangmuirMembraneTrapSource, ResolvedTrappingMembraneError, admit_traps},
};
use super::{
    Bindings, ClosedGasParcel, GasSpeciesBoundaryStep, IsothermalGasMembraneMixture,
    IsothermalGasMembraneSolve, MembraneDiffusionSourceSupport, PressureSample,
    SpeciesBoundaryRelativeTransfer, SpeciesId, Sum, ThermochemicalGasModel, admit_independent,
    finite, gas_error, inventory, isothermal_stream, poll,
};
use crate::adsorption::{
    LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_DESORPTION_RATE_PROPERTY,
    LANGMUIR_REFERENCE_AFFINITY_PROPERTY, LangmuirProperties, competitive::CompetitiveAdsorbent,
};
use crate::thermochemical_gas::energy::GasEnergyError;
use std::{collections::BTreeMap, sync::Arc};

/// One explicitly shared source lattice per layer. No capacity is inferred by
/// adding independent species capacities. Nonadsorbing catalog species remain
/// present in the original mobile membrane and both actual gas mixtures.
#[derive(Clone, Debug)]
pub struct CompetitiveMembraneTrapSource {
    pub adsorbent: Arc<CompetitiveAdsorbent>,
    /// kg dry common adsorbent skeleton per m3 of the complete membrane layer.
    pub dry_skeleton_density_kg_m3: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsothermalCompetitiveGasError {
    Gas(IsothermalGasMembraneError),
    Source(ResolvedTrappingMembraneError),
    Competitive(CompetitiveMembraneError),
    Cancelled,
}
impl core::fmt::Display for IsothermalCompetitiveGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-qualified competitive gas membrane: {self:?}")
    }
}
impl std::error::Error for IsothermalCompetitiveGasError {}
impl From<IsothermalGasMembraneError> for IsothermalCompetitiveGasError {
    fn from(e: IsothermalGasMembraneError) -> Self {
        if e == IsothermalGasMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
impl From<ResolvedTrappingMembraneError> for IsothermalCompetitiveGasError {
    fn from(e: ResolvedTrappingMembraneError) -> Self {
        if e == ResolvedTrappingMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Source(e)
        }
    }
}
impl From<CompetitiveMembraneError> for IsothermalCompetitiveGasError {
    fn from(e: CompetitiveMembraneError) -> Self {
        if e == CompetitiveMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Competitive(e)
        }
    }
}
impl From<ResolvedMembraneDiffusionError> for IsothermalCompetitiveGasError {
    fn from(e: ResolvedMembraneDiffusionError) -> Self {
        ResolvedTrappingMembraneError::from(e).into()
    }
}
impl From<MembraneDiffusionError> for IsothermalCompetitiveGasError {
    fn from(e: MembraneDiffusionError) -> Self {
        CompetitiveMembraneError::from(e).into()
    }
}
impl From<GasEnergyError> for IsothermalCompetitiveGasError {
    fn from(e: GasEnergyError) -> Self {
        IsothermalGasMembraneError::from(e).into()
    }
}
type Result<T> = core::result::Result<T, IsothermalCompetitiveGasError>;
fn invalid(message: &'static str) -> IsothermalCompetitiveGasError {
    IsothermalGasMembraneError::Invalid(message).into()
}

#[derive(Clone, Debug)]
struct SourceChannel<'a> {
    bindings: Arc<Bindings<'a>>,
    traps: Arc<[Option<LangmuirMembraneTrapSource>]>,
    coefficients: Arc<[Option<LangmuirMembraneTrap>]>,
    molar_mass: f64,
    molar_enthalpy: f64,
}

#[derive(Clone, Debug)]
pub struct IsothermalCompetitiveGasMembrane<'a> {
    channels: Arc<[SourceChannel<'a>]>,
    finite: FiniteCompetitiveMembraneDiffusionState,
    gases: [ClosedGasParcel; 2],
    catalog: Arc<ThermochemicalGasModel>,
    support: Vec<MembraneDiffusionSourceSupport>,
    trap_support: Vec<Vec<Option<[LangmuirProperties; 2]>>>,
    gas_controls: IsothermalGasMembraneSolve,
    controls: CompetitiveMembraneSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveGasSpeciesAccount {
    pub species: SpeciesId,
    pub actual_partial_pressures_pa: [f64; 2],
    pub pressure_mismatch_pa: [f64; 2],
    pub pressure_mismatch_mass_kg: [f64; 2],
    pub stream_enthalpy_into_gas_j: [f64; 2],
    pub thermostat_heat_into_gas_j: [f64; 2],
    pub actual_gas_mass_changes_kg: [f64; 2],
    pub mobile_mass_change_kg: f64,
    pub trapped_mass_change_kg: f64,
    pub total_mass_before_kg: f64,
    pub total_mass_after_kg: f64,
    pub mass_residual_kg: f64,
    pub total_mass_residual_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalCompetitiveGasStep {
    pub diffusion: CompetitiveMembraneStep,
    /// Same canonical SpeciesId order as the native diffusion report.
    pub species: Vec<CompetitiveGasSpeciesAccount>,
    pub gas: [GasSpeciesBoundaryStep; 2],
    pub actual_reservoir_pressures_pa: [f64; 2],
    pub stream_enthalpy_into_gas_j: [f64; 2],
    pub thermostat_heat_into_gas_j: [f64; 2],
    pub gas_energy_residual_j: f64,
    pub actual_gas_mass_changes_kg: [f64; 2],
    pub membrane_mass_change_kg: f64,
    pub total_mass_before_kg: f64,
    pub total_mass_after_kg: f64,
    pub mass_residual_kg: f64,
    pub total_mass_residual_kg: f64,
}

#[derive(Clone, Debug)]
pub struct IsothermalCompetitiveGasCandidate<'a> {
    state: IsothermalCompetitiveGasMembrane<'a>,
    step: IsothermalCompetitiveGasStep,
}
impl<'a> IsothermalCompetitiveGasCandidate<'a> {
    pub const fn state(&self) -> &IsothermalCompetitiveGasMembrane<'a> {
        &self.state
    }
    pub const fn step(&self) -> &IsothermalCompetitiveGasStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        IsothermalCompetitiveGasMembrane<'a>,
        IsothermalCompetitiveGasStep,
    ) {
        (self.state, self.step)
    }
}

fn qualify_traps(
    bindings: &Bindings<'_>,
    traps: &[Option<LangmuirMembraneTrapSource>],
    support: &MembraneDiffusionSourceSupport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(
    Vec<Option<LangmuirMembraneTrap>>,
    Vec<Option<[LangmuirProperties; 2]>>,
)> {
    let (coefficients, receipts) = admit_traps(
        &bindings.species,
        bindings.temperature_k,
        traps,
        support,
        cancelled,
    )?;
    for receipt in receipts.iter().flatten().flatten() {
        poll(cancelled)?;
        for key in [
            LANGMUIR_CAPACITY_PROPERTY,
            LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
            LANGMUIR_DESORPTION_RATE_PROPERTY,
        ] {
            admit_independent(
                receipt
                    .source
                    .property(key)
                    .ok_or(invalid("competitive trap lost its original source claim"))?,
            )?;
        }
    }
    Ok((coefficients, receipts))
}

fn sample(
    state: &FiniteCompetitiveMembraneDiffusionState,
    totals: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PressureSample> {
    let fields: Vec<_> = state
        .reservoir_partial_pressures_pa()
        .iter()
        .enumerate()
        .map(|(s, &p)| (state.membrane().mobile(s).expect("admitted species"), p))
        .collect();
    Ok(PressureSample::collect_mobile(&fields, totals, cancelled)?)
}

impl<'a> IsothermalGasMembraneMixture<'a> {
    /// Attach one shared lattice per physical layer to the already resolved
    /// original mixture. Every catalog species needs an explicit initial
    /// trapped profile, including zeros for species without a site law.
    /// Existing gas budgets are retained; `controls.balance` owns the new
    /// simultaneous mobile/site/finite-reservoir solve.
    pub fn with_competitive_sites(
        self,
        sites: &[Option<CompetitiveMembraneTrapSource>],
        initial_trapped_moles: &[(SpeciesId, Vec<f64>)],
        controls: CompetitiveMembraneSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalCompetitiveGasMembrane<'a>> {
        poll(&mut cancelled)?;
        if sites.len() != self.channels[0].bindings.layers.len()
            || initial_trapped_moles.len() != self.channels.len()
        {
            return Err(invalid(
                "one shared-site slot per layer and one initial profile per catalog species required",
            ));
        }
        let mut initial = BTreeMap::new();
        for (id, loading) in initial_trapped_moles {
            if self.index(id).is_none() || initial.insert(id, loading).is_some() {
                return Err(invalid(
                    "initial shared-site profiles must match the unique complete gas roster",
                ));
            }
        }
        for site in sites.iter().flatten() {
            poll(&mut cancelled)?;
            if !site.dry_skeleton_density_kg_m3.is_finite()
                || site.dry_skeleton_density_kg_m3 <= 0.0
            {
                return Err(invalid(
                    "shared-site dry skeleton density must be finite positive",
                ));
            }
            for law in site.adsorbent.channels() {
                if self.index(law.species()).is_none() {
                    return Err(invalid(
                        "every competitive lattice species needs an original gas/mobile channel",
                    ));
                }
                law.admit_catalog(&self.catalog)
                    .map_err(ResolvedTrappingMembraneError::from)?;
            }
        }
        let capacities = self.channels[0].finite.capacities();
        let pressures = self
            .channels
            .iter()
            .map(|c| c.finite.reservoir_pressures_pa())
            .collect();
        let mut channels = Vec::with_capacity(self.channels.len());
        let mut numerical = Vec::with_capacity(self.channels.len());
        let mut trap_support = Vec::with_capacity(self.channels.len());
        for (s, channel) in self.channels.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let id = &channel.bindings.species;
            let traps: Vec<_> = sites
                .iter()
                .map(|site| {
                    site.as_ref().and_then(|site| {
                        site.adsorbent
                            .channels()
                            .iter()
                            .find(|law| law.species() == id)
                            .map(|law| LangmuirMembraneTrapSource {
                                adsorbent: Arc::new(law.clone()),
                                dry_skeleton_density_kg_m3: site.dry_skeleton_density_kg_m3,
                            })
                    })
                })
                .collect();
            let (coefficients, receipts) =
                qualify_traps(&channel.bindings, &traps, &self.support[s], &mut cancelled)?;
            numerical.push(CompetitiveMembraneChannel {
                species: id.clone(),
                mobile: channel.finite.into_membrane(),
                traps: coefficients.clone(),
                initial_trapped_moles: initial[id].clone(),
            });
            trap_support.push(receipts);
            channels.push(SourceChannel {
                bindings: channel.bindings,
                traps: traps.into(),
                coefficients: coefficients.into(),
                molar_mass: channel.molar_mass,
                molar_enthalpy: channel.molar_enthalpy,
            });
        }
        let membrane = CompetitiveMembraneDiffusionState::new(numerical, &mut cancelled)?;
        let finite = FiniteCompetitiveMembraneDiffusionState::new(membrane, capacities, pressures)?;
        let state = IsothermalCompetitiveGasMembrane {
            channels: channels.into(),
            finite,
            gases: self.gases,
            catalog: self.catalog,
            support: self.support,
            trap_support,
            gas_controls: self.controls,
            controls,
        };
        state.admit_capacity_mass()?;
        state.total_mass_kg()?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}

impl<'a> IsothermalCompetitiveGasMembrane<'a> {
    pub fn temperature_k(&self) -> f64 {
        self.channels[0].bindings.temperature_k
    }
    pub fn elapsed_s(&self) -> f64 {
        self.membrane().elapsed_s()
    }
    pub const fn gases(&self) -> &[ClosedGasParcel; 2] {
        &self.gases
    }
    pub fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }
    pub const fn controls(&self) -> CompetitiveMembraneSolve {
        self.controls
    }
    pub const fn gas_controls(&self) -> IsothermalGasMembraneSolve {
        self.gas_controls
    }
    pub const fn finite_state(&self) -> &FiniteCompetitiveMembraneDiffusionState {
        &self.finite
    }
    pub fn membrane(&self) -> &CompetitiveMembraneDiffusionState {
        self.finite.membrane()
    }
    pub fn species(&self) -> impl ExactSizeIterator<Item = &SpeciesId> {
        self.membrane().species()
    }
    fn index(&self, id: &SpeciesId) -> Option<usize> {
        self.channels
            .binary_search_by(|c| c.bindings.species.cmp(id))
            .ok()
    }
    pub fn source_support(&self, id: &SpeciesId) -> Option<&MembraneDiffusionSourceSupport> {
        self.index(id).map(|s| &self.support[s])
    }
    pub fn trap_source_support(
        &self,
        id: &SpeciesId,
    ) -> Option<&[Option<[LangmuirProperties; 2]>]> {
        self.index(id).map(|s| self.trap_support[s].as_slice())
    }
    pub fn workspace(&self, maximum_bytes: usize) -> Result<CompetitiveMembraneWorkspace> {
        Ok(self.finite.workspace(maximum_bytes)?)
    }
    pub fn total_mass_kg(&self) -> Result<f64> {
        let mut sum = Sum::default();
        sum.add(self.gases[0].mass_kg())?;
        sum.add(self.gases[1].mass_kg())?;
        for (s, channel) in self.channels.iter().enumerate() {
            sum.add(finite(
                self.membrane().species_amount_mol(s)? * channel.molar_mass,
            )?)?;
        }
        Ok(sum.value)
    }
    fn admit_mass(&self, residual_kg: f64) -> Result<()> {
        if !residual_kg.is_finite() || residual_kg.abs() > self.gas_controls.mass_tolerance_kg {
            Err(IsothermalGasMembraneError::MassBalance {
                residual_kg,
                tolerance_kg: self.gas_controls.mass_tolerance_kg,
            }
            .into())
        } else {
            Ok(())
        }
    }
    fn pressure_budget(&self, side: usize) -> f64 {
        self.gas_controls
            .mass_tolerance_kg
            .min(self.gas_controls.gas[side].mass_tolerance_kg)
    }
    fn admit_capacity_mass(&self) -> Result<()> {
        let inventories = [
            inventory(self.gases[0].model(), self.gases[0].mass_kg())?,
            inventory(self.gases[1].model(), self.gases[1].mass_kg())?,
        ];
        let mut sum = Sum::default();
        for (s, channel) in self.channels.iter().enumerate() {
            for side in 0..2 {
                let numerical = finite(
                    self.finite.capacities()[side].value_si()
                        * self.finite.reservoir_partial_pressures_pa()[s][side],
                )?;
                let actual = inventories[side]
                    .get(&channel.bindings.species)
                    .copied()
                    .unwrap_or(0.0);
                let error = finite(actual - finite(numerical * channel.molar_mass)?)?;
                if error.abs() > self.pressure_budget(side) {
                    return Err(IsothermalGasMembraneError::PressureClosure {
                        reservoir: side,
                        mass_error_kg: error,
                        tolerance_kg: self.pressure_budget(side),
                    }
                    .into());
                }
                sum.add(error.abs())?;
            }
        }
        self.admit_mass(sum.value)
    }

    /// Solve all species, shared vacancies and both finite baths together,
    /// reconstruct the two original gases once, then qualify every actual old,
    /// solved and EOS-reseeded source envelope. No failure publishes state.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        work: &mut CompetitiveMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalCompetitiveGasCandidate<'a>> {
        poll(&mut cancelled)?;
        let (numerical, diffusion) = self
            .finite
            .step_candidate(duration_s, self.controls, work, &mut cancelled)?
            .into_parts();
        let mut transfers: [Vec<SpeciesBoundaryRelativeTransfer>; 2] = [Vec::new(), Vec::new()];
        let mut streams = Vec::with_capacity(self.channels.len());
        let mut heats = Vec::with_capacity(self.channels.len());
        for (channel, (_, report)) in self.channels.iter().zip(&diffusion.species) {
            poll(&mut cancelled)?;
            let amounts = [
                -report.boundary_transfers_mol[0],
                report.boundary_transfers_mol[1],
            ];
            let mut stream = [0.0; 2];
            let mut heat = [0.0; 2];
            for side in 0..2 {
                let (transfer, enthalpy) = isothermal_stream(
                    &channel.bindings.species,
                    amounts[side],
                    channel.molar_mass,
                    channel.molar_enthalpy,
                    self.temperature_k(),
                )?;
                stream[side] = enthalpy;
                heat[side] = transfer.relative_energy_into_gas_j;
                transfers[side].push(transfer);
            }
            streams.push(stream);
            heats.push(heat);
        }
        let mut gases = Vec::with_capacity(2);
        let mut gas_reports = Vec::with_capacity(2);
        for side in 0..2 {
            poll(&mut cancelled)?;
            let (gas, _, report) = self.gases[side]
                .exchange_species_boundary_relative(
                    &self.catalog,
                    &transfers[side],
                    self.gas_controls.gas[side],
                    &mut cancelled,
                )
                .map_err(|e| gas_error(side, e))?
                .into_parts();
            if gas.state().state().temperature != self.temperature_k() {
                return Err(invalid(
                    "original mixture reconstruction changed held temperature",
                ));
            }
            gases.push(gas);
            gas_reports.push(report);
        }
        let gases: [ClosedGasParcel; 2] = gases
            .try_into()
            .map_err(|_| invalid("two gases required"))?;
        let gas: [GasSpeciesBoundaryStep; 2] = gas_reports
            .try_into()
            .map_err(|_| invalid("two gas reports required"))?;
        let old = sample(
            &self.finite,
            self.gases.each_ref().map(|g| g.state().state().pressure),
            &mut cancelled,
        )?;
        let mut totals = [Sum::default(), Sum::default()];
        for p in numerical.reservoir_partial_pressures_pa() {
            for side in 0..2 {
                totals[side].add(p[side])?;
            }
        }
        let solved = sample(&numerical, totals.map(|s| s.value), &mut cancelled)?;
        let before_inventory = [
            inventory(self.gases[0].model(), self.gases[0].mass_kg())?,
            inventory(self.gases[1].model(), self.gases[1].mass_kg())?,
        ];
        let after_inventory = [
            inventory(gases[0].model(), gases[0].mass_kg())?,
            inventory(gases[1].model(), gases[1].mass_kg())?,
        ];
        let mut actual_partial = Vec::with_capacity(self.channels.len());
        let mut rows = Vec::with_capacity(self.channels.len());
        let mut storage = Sum::default();
        let mut stream_total = [Sum::default(), Sum::default()];
        let mut heat_total = [Sum::default(), Sum::default()];
        for (s, channel) in self.channels.iter().enumerate() {
            poll(&mut cancelled)?;
            let id = &channel.bindings.species;
            let partial = [
                gases[0].species_partial_pressure_pa(id)?,
                gases[1].species_partial_pressure_pa(id)?,
            ];
            let mut pressure_error = [0.0; 2];
            let mut mass_error = [0.0; 2];
            let mut changes = [0.0; 2];
            let mut before = Sum::default();
            let mut after = Sum::default();
            let mut balance = Sum::default();
            for side in 0..2 {
                pressure_error[side] =
                    finite(partial[side] - numerical.reservoir_partial_pressures_pa()[s][side])?;
                mass_error[side] = finite(
                    numerical.capacities()[side].value_si()
                        * pressure_error[side]
                        * channel.molar_mass,
                )?;
                if mass_error[side].abs() > self.pressure_budget(side) {
                    return Err(IsothermalGasMembraneError::PressureClosure {
                        reservoir: side,
                        mass_error_kg: mass_error[side],
                        tolerance_kg: self.pressure_budget(side),
                    }
                    .into());
                }
                let a = before_inventory[side].get(id).copied().unwrap_or(0.0);
                let b = after_inventory[side].get(id).copied().unwrap_or(0.0);
                changes[side] = finite(b - a)?;
                balance.add(changes[side])?;
                before.add(a)?;
                after.add(b)?;
                stream_total[side].add(streams[s][side])?;
                heat_total[side].add(heats[s][side])?;
            }
            let report = &diffusion.species[s].1;
            let mobile = finite(report.mobile_change_mol * channel.molar_mass)?;
            let trapped = finite(report.trapped_change_mol * channel.molar_mass)?;
            balance.add(mobile)?;
            balance.add(trapped)?;
            storage.add(mobile)?;
            storage.add(trapped)?;
            self.admit_mass(balance.value)?;
            before.add(finite(
                self.membrane().species_amount_mol(s)? * channel.molar_mass,
            )?)?;
            after.add(finite(
                numerical.membrane().species_amount_mol(s)? * channel.molar_mass,
            )?)?;
            self.admit_mass(after.value - before.value)?;
            rows.push(CompetitiveGasSpeciesAccount {
                species: id.clone(),
                actual_partial_pressures_pa: partial,
                pressure_mismatch_pa: pressure_error,
                pressure_mismatch_mass_kg: mass_error,
                stream_enthalpy_into_gas_j: streams[s],
                thermostat_heat_into_gas_j: heats[s],
                actual_gas_mass_changes_kg: changes,
                mobile_mass_change_kg: mobile,
                trapped_mass_change_kg: trapped,
                total_mass_before_kg: before.value,
                total_mass_after_kg: after.value,
                mass_residual_kg: balance.value,
                total_mass_residual_kg: after.value - before.value,
            });
            actual_partial.push(partial);
        }
        // EOS rounding never feeds back into a fictional gas inventory: actual
        // original partial pressures become the next finite-bath coordinates.
        let capacities = numerical.capacities();
        let finite_state = FiniteCompetitiveMembraneDiffusionState::new(
            numerical.into_membrane(),
            capacities,
            actual_partial,
        )?;
        let pressures = gases.each_ref().map(|g| g.state().state().pressure);
        let actual = sample(&finite_state, pressures, &mut cancelled)?;
        let mut support = Vec::with_capacity(self.channels.len());
        let mut trap_support = Vec::with_capacity(self.channels.len());
        for (s, channel) in self.channels.iter().enumerate() {
            poll(&mut cancelled)?;
            let mobile = finite_state.membrane().mobile(s).expect("admitted species");
            let receipt = channel.bindings.admit_pressures(
                [
                    (mobile, old.cells.as_slice(), old.faces.as_slice()),
                    (mobile, solved.cells.as_slice(), solved.faces.as_slice()),
                    (mobile, actual.cells.as_slice(), actual.faces.as_slice()),
                ],
                &mut cancelled,
            )?;
            let (coefficients, traps) =
                qualify_traps(&channel.bindings, &channel.traps, &receipt, &mut cancelled)?;
            if coefficients.as_slice() != channel.coefficients.as_ref() {
                return Err(invalid(
                    "competitive candidate changed original site coefficients or density",
                ));
            }
            support.push(receipt);
            trap_support.push(traps);
        }
        let state = Self {
            channels: self.channels.clone(),
            finite: finite_state,
            gases,
            catalog: self.catalog.clone(),
            support,
            trap_support,
            gas_controls: self.gas_controls,
            controls: self.controls,
        };
        state.admit_capacity_mass()?;
        let changes = [
            state.gases[0].mass_kg() - self.gases[0].mass_kg(),
            state.gases[1].mass_kg() - self.gases[1].mass_kg(),
        ];
        let mut balance = Sum::default();
        balance.add(changes[0])?;
        balance.add(changes[1])?;
        balance.add(storage.value)?;
        state.admit_mass(balance.value)?;
        let before = self.total_mass_kg()?;
        let after = state.total_mass_kg()?;
        state.admit_mass(after - before)?;
        let gas_energy_residual_j = finite(gas[0].energy_residual_j + gas[1].energy_residual_j)?;
        let step = IsothermalCompetitiveGasStep {
            diffusion,
            species: rows,
            gas,
            actual_reservoir_pressures_pa: pressures,
            stream_enthalpy_into_gas_j: stream_total.map(|s| s.value),
            thermostat_heat_into_gas_j: heat_total.map(|s| s.value),
            gas_energy_residual_j,
            actual_gas_mass_changes_kg: changes,
            membrane_mass_change_kg: storage.value,
            total_mass_before_kg: before,
            total_mass_after_kg: after,
            mass_residual_kg: balance.value,
            total_mass_residual_kg: after - before,
        };
        poll(&mut cancelled)?;
        Ok(IsothermalCompetitiveGasCandidate { state, step })
    }
}
