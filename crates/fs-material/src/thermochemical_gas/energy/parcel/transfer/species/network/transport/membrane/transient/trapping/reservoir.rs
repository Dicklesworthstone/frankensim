//! Original pure-gas reservoirs around source-qualified kinetic membrane traps.
//!
//! The existing gas boundary owner reconstructs every proposed parcel. Its
//! source-relative isothermal stream cancels flow work before formation
//! references enter. Mobile diffusion, immobile kinetics and both gas amounts
//! share the existing implicit solve; actual EOS pressures are then reseeded.
//! These are gas-side energy accounts. No mobile/trap caloric storage, heat of
//! sorption, nonisothermal film or whole-system energy closure is inferred.

use super::super::reservoir::{
    self as gas_owner, GasBinding, IsothermalGasMembraneError, IsothermalGasMembraneSolve, finite,
    gas_error, isothermal_stream, poll,
};
use super::super::slab::trapping::{
    FiniteTrappingMembraneDiffusionState, TrappingMembraneDiffusionState, TrappingMembraneError,
    TrappingMembraneSolve, TrappingMembraneStep, TrappingMembraneWorkspace,
};
use super::super::{MembraneDiffusionError, ResolvedMembraneDiffusionError};
use super::{
    ResolvedFiniteTrappingMembraneDiffusion, ResolvedTrappingMembraneDiffusion,
    ResolvedTrappingMembraneError, admit_traps,
};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::{
            ClosedGasParcel,
            transfer::species::{
                Sum,
                boundary::{GasSpeciesBoundarySolve, GasSpeciesBoundaryStep},
            },
        },
    },
};
use fs_thermochem::SpeciesId;
use std::sync::Arc;

mod adaptive;
pub use adaptive::IsothermalTrappingGasAdvanceReport;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalTrappingGasSolve {
    pub membrane: TrappingMembraneSolve,
    pub gas: [GasSpeciesBoundarySolve; 2],
    /// Actual gas plus mobile and trapped mass, after EOS pressure reseeding.
    pub mass_tolerance_kg: f64,
}
impl IsothermalTrappingGasSolve {
    fn gas_controls(self) -> IsothermalGasMembraneSolve {
        IsothermalGasMembraneSolve {
            membrane: self.membrane.balance,
            gas: self.gas,
            mass_tolerance_kg: self.mass_tolerance_kg,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsothermalTrappingGasError {
    Source(ResolvedTrappingMembraneError),
    Gas(IsothermalGasMembraneError),
    Cancelled,
}
impl core::fmt::Display for IsothermalTrappingGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isothermal trapping membrane: {self:?}")
    }
}
impl std::error::Error for IsothermalTrappingGasError {}
impl From<ResolvedTrappingMembraneError> for IsothermalTrappingGasError {
    fn from(e: ResolvedTrappingMembraneError) -> Self {
        if e == ResolvedTrappingMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Source(e)
        }
    }
}
impl From<IsothermalGasMembraneError> for IsothermalTrappingGasError {
    fn from(e: IsothermalGasMembraneError) -> Self {
        if e == IsothermalGasMembraneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
impl From<TrappingMembraneError> for IsothermalTrappingGasError {
    fn from(e: TrappingMembraneError) -> Self {
        ResolvedTrappingMembraneError::from(e).into()
    }
}
impl From<ResolvedMembraneDiffusionError> for IsothermalTrappingGasError {
    fn from(e: ResolvedMembraneDiffusionError) -> Self {
        ResolvedTrappingMembraneError::from(e).into()
    }
}
impl From<MembraneDiffusionError> for IsothermalTrappingGasError {
    fn from(e: MembraneDiffusionError) -> Self {
        ResolvedTrappingMembraneError::from(e).into()
    }
}
impl From<GasEnergyError> for IsothermalTrappingGasError {
    fn from(e: GasEnergyError) -> Self {
        IsothermalGasMembraneError::from(e).into()
    }
}
type Result<T> = core::result::Result<T, IsothermalTrappingGasError>;
use IsothermalGasMembraneError as GasError;

/// Net fixed-temperature gas streams and independently audited actual masses.
#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalTrappingGasAccount {
    pub amount_into_gas_mol: [f64; 2],
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
pub struct IsothermalTrappingGasStep {
    pub diffusion: TrappingMembraneStep,
    pub gas: [GasSpeciesBoundaryStep; 2],
    pub accounting: IsothermalTrappingGasAccount,
    /// Actual original-EOS pressure minus the simultaneous numerical endpoint.
    pub pressure_mismatch_pa: [f64; 2],
}

#[derive(Clone, Debug)]
pub struct IsothermalTrappingGasMembrane<'a> {
    finite: ResolvedFiniteTrappingMembraneDiffusion<'a>,
    gases: [ClosedGasParcel; 2],
    catalog: Arc<ThermochemicalGasModel>,
    controls: IsothermalTrappingGasSolve,
    molar_mass_kg_mol: f64,
    molar_enthalpy_j_mol: f64,
}
#[derive(Clone, Debug)]
pub struct IsothermalTrappingGasCandidate<'a, R = IsothermalTrappingGasStep> {
    state: IsothermalTrappingGasMembrane<'a>,
    report: R,
}
impl<'a, R> IsothermalTrappingGasCandidate<'a, R> {
    pub const fn state(&self) -> &IsothermalTrappingGasMembrane<'a> {
        &self.state
    }
    pub const fn report(&self) -> &R {
        &self.report
    }
    pub fn into_parts(self) -> (IsothermalTrappingGasMembrane<'a>, R) {
        (self.state, self.report)
    }
}
impl<'a> ResolvedTrappingMembraneDiffusion<'a> {
    pub fn with_finite_gas_reservoirs(
        self,
        gases: [ClosedGasParcel; 2],
        controls: IsothermalTrappingGasSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalTrappingGasMembrane<'a>> {
        IsothermalTrappingGasMembrane::new(self, gases, controls, cancelled)
    }
}
impl<'a> IsothermalTrappingGasMembrane<'a> {
    pub fn new(
        membrane: ResolvedTrappingMembraneDiffusion<'a>,
        gases: [ClosedGasParcel; 2],
        controls: IsothermalTrappingGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        controls.gas_controls().validate(membrane.temperature_k())?;
        let GasBinding {
            catalog,
            capacities,
            pressures,
            molar_mass_kg_mol,
            molar_enthalpy_j_mol,
        } = gas_owner::bind_gases(membrane.species(), membrane.temperature_k(), &gases)?;
        membrane.bindings.admit_gas_surfaces(&catalog)?;
        for trap in membrane.trap_bindings.iter().flatten() {
            poll(&mut cancelled)?;
            trap.adsorbent
                .admit_catalog(&catalog)
                .map_err(ResolvedTrappingMembraneError::from)?;
        }
        let finite = ResolvedFiniteTrappingMembraneDiffusion::new(
            membrane,
            capacities,
            pressures,
            &mut cancelled,
        )?;
        let state = Self {
            finite,
            gases,
            catalog,
            controls,
            molar_mass_kg_mol,
            molar_enthalpy_j_mol,
        };
        state.admit_capacity_mass()?;
        state.total_mass_kg()?;
        poll(&mut cancelled)?;
        Ok(state)
    }
    pub fn species(&self) -> &SpeciesId {
        self.finite.species()
    }
    pub fn temperature_k(&self) -> f64 {
        self.finite.temperature_k()
    }
    pub fn membrane(&self) -> &TrappingMembraneDiffusionState {
        self.finite.membrane()
    }
    pub const fn finite_state(&self) -> &ResolvedFiniteTrappingMembraneDiffusion<'a> {
        &self.finite
    }
    pub const fn gases(&self) -> &[ClosedGasParcel; 2] {
        &self.gases
    }
    pub fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }
    pub const fn controls(&self) -> IsothermalTrappingGasSolve {
        self.controls
    }
    pub fn workspace(&self) -> Result<TrappingMembraneWorkspace> {
        Ok(self.finite.workspace()?)
    }
    pub fn total_mass_kg(&self) -> Result<f64> {
        let mut sum = Sum::default();
        sum.add(self.gases[0].mass_kg())?;
        sum.add(finite(
            self.membrane().stored_amount_mol()? * self.molar_mass_kg_mol,
        )?)?;
        sum.add(self.gases[1].mass_kg())?;
        Ok(sum.value)
    }
    fn admit_mass(&self, residual_kg: f64) -> Result<()> {
        if !residual_kg.is_finite() || residual_kg.abs() > self.controls.mass_tolerance_kg {
            return Err(GasError::MassBalance {
                residual_kg,
                tolerance_kg: self.controls.mass_tolerance_kg,
            }
            .into());
        }
        Ok(())
    }
    fn admit_pressure_mass(&self, reservoir: usize, mass_error_kg: f64) -> Result<()> {
        let tolerance_kg = self
            .controls
            .mass_tolerance_kg
            .min(self.controls.gas[reservoir].mass_tolerance_kg);
        if !mass_error_kg.is_finite() || mass_error_kg.abs() > tolerance_kg {
            return Err(GasError::PressureClosure {
                reservoir,
                mass_error_kg,
                tolerance_kg,
            }
            .into());
        }
        Ok(())
    }
    fn admit_capacity_mass(&self) -> Result<()> {
        let mut sum = Sum::default();
        for i in 0..2 {
            let mass = finite(
                (self.finite.capacities()[i].value_si() * self.finite.reservoir_pressures_pa()[i])
                    * self.molar_mass_kg_mol,
            )?;
            let residual = finite(self.gases[i].mass_kg() - mass)?;
            self.admit_pressure_mass(i, residual)?;
            sum.add(residual)?;
        }
        self.admit_mass(sum.value)
    }
    fn account(
        &self,
        after: &Self,
        amounts: [f64; 2],
        mobile_change: f64,
        trapped_change: f64,
    ) -> Result<IsothermalTrappingGasAccount> {
        let mut stream = [0.0; 2];
        let mut heat = [0.0; 2];
        let mut actual = [0.0; 2];
        let mut sum = Sum::default();
        for i in 0..2 {
            let (transfer, enthalpy) = isothermal_stream(
                self.species(),
                amounts[i],
                self.molar_mass_kg_mol,
                self.molar_enthalpy_j_mol,
                self.temperature_k(),
            )?;
            stream[i] = enthalpy;
            heat[i] = transfer.relative_energy_into_gas_j;
            actual[i] = finite(after.gases[i].mass_kg() - self.gases[i].mass_kg())?;
            self.admit_pressure_mass(i, finite(actual[i] - transfer.mass_into_gas_kg)?)?;
            sum.add(actual[i])?;
        }
        let mobile_mass_change_kg = finite(mobile_change * self.molar_mass_kg_mol)?;
        let trapped_mass_change_kg = finite(trapped_change * self.molar_mass_kg_mol)?;
        sum.add(mobile_mass_change_kg)?;
        sum.add(trapped_mass_change_kg)?;
        self.admit_mass(sum.value)?;
        let before = self.total_mass_kg()?;
        let after = after.total_mass_kg()?;
        self.admit_mass(after - before)?;
        Ok(IsothermalTrappingGasAccount {
            amount_into_gas_mol: amounts,
            stream_enthalpy_into_gas_j: stream,
            thermostat_heat_into_gas_j: heat,
            actual_gas_mass_changes_kg: actual,
            mobile_mass_change_kg,
            trapped_mass_change_kg,
            total_mass_before_kg: before,
            total_mass_after_kg: after,
            mass_residual_kg: sum.value,
            total_mass_residual_kg: after - before,
        })
    }
    /// One source-qualified implicit step followed by original gas reconstruction.
    /// Both gas caloric reports remain explicit; no accepted state is mutated.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        work: &mut TrappingMembraneWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalTrappingGasCandidate<'a>> {
        self.step_with_controls(duration_s, self.controls.membrane, work, cancelled)
    }
    fn step_with_controls(
        &self,
        duration_s: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalTrappingGasCandidate<'a>> {
        poll(&mut cancelled)?;
        let (mut numerical, diffusion) = self
            .finite
            .step_candidate(duration_s, controls, work, &mut cancelled)?
            .into_parts();
        let amounts = [
            -diffusion.boundary_transfers_mol[0],
            diffusion.boundary_transfers_mol[1],
        ];
        let mut gases = Vec::with_capacity(2);
        let mut rows = Vec::with_capacity(2);
        let mut mismatch = [0.0; 2];
        for i in 0..2 {
            poll(&mut cancelled)?;
            let (transfer, _) = isothermal_stream(
                self.species(),
                amounts[i],
                self.molar_mass_kg_mol,
                self.molar_enthalpy_j_mol,
                self.temperature_k(),
            )?;
            let (gas, _, report) = self.gases[i]
                .exchange_species_boundary_relative(
                    &self.catalog,
                    &[transfer],
                    self.controls.gas[i],
                    &mut cancelled,
                )
                .map_err(|e| gas_error(i, e))?
                .into_parts();
            if gas.state().state().temperature != self.temperature_k() {
                return Err(GasError::Invalid(
                    "original gas reconstruction changed the held trap-film temperature",
                )
                .into());
            }
            mismatch[i] =
                finite(gas.state().state().pressure - numerical.reservoir_pressures_pa()[i])?;
            self.admit_pressure_mass(
                i,
                finite(
                    (numerical.capacities()[i].value_si() * mismatch[i]) * self.molar_mass_kg_mol,
                )?,
            )?;
            gases.push(gas);
            rows.push(report);
        }
        let gases: [ClosedGasParcel; 2] = gases
            .try_into()
            .map_err(|_| GasError::Invalid("two original gas endpoints required"))?;
        let gas: [GasSpeciesBoundaryStep; 2] = rows
            .try_into()
            .map_err(|_| GasError::Invalid("two original gas reports required"))?;
        let pressures = gases.each_ref().map(|g| g.state().state().pressure);
        let before = self
            .membrane()
            .mobile()
            .sample_fluxes(self.finite.reservoir_pressures_pa(), &mut cancelled)?;
        let actual = numerical
            .membrane()
            .mobile()
            .sample_fluxes(pressures, &mut cancelled)?;
        let support = numerical.bindings.admit(
            [
                (
                    self.membrane().mobile(),
                    before.face_partial_pressures_pa.as_slice(),
                ),
                (
                    numerical.membrane().mobile(),
                    diffusion.transport.face_partial_pressures_pa.as_slice(),
                ),
                (
                    numerical.membrane().mobile(),
                    actual.face_partial_pressures_pa.as_slice(),
                ),
            ],
            &mut cancelled,
        )?;
        let (coefficients, trap_support) = admit_traps(
            numerical.species(),
            numerical.temperature_k(),
            &numerical.trap_bindings,
            &support,
            &mut cancelled,
        )?;
        if coefficients.as_slice() != numerical.membrane().traps() {
            return Err(ResolvedTrappingMembraneError::Invalid(
                "actual gas pressures changed original trap coefficients",
            )
            .into());
        }
        let capacities = numerical.capacities();
        numerical.state = FiniteTrappingMembraneDiffusionState::new(
            numerical.state.into_membrane(),
            capacities,
            pressures,
        )?;
        numerical.support = support;
        numerical.trap_support = trap_support;
        let state = Self {
            finite: numerical,
            gases,
            catalog: self.catalog.clone(),
            controls: self.controls,
            molar_mass_kg_mol: self.molar_mass_kg_mol,
            molar_enthalpy_j_mol: self.molar_enthalpy_j_mol,
        };
        state.admit_capacity_mass()?;
        let accounting = self.account(
            &state,
            amounts,
            diffusion.mobile_change_mol,
            diffusion.trapped_change_mol,
        )?;
        poll(&mut cancelled)?;
        Ok(IsothermalTrappingGasCandidate {
            state,
            report: IsothermalTrappingGasStep {
                diffusion,
                gas,
                accounting,
                pressure_mismatch_pa: mismatch,
            },
        })
    }
}
