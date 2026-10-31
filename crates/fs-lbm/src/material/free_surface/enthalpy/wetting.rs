//! Thermal and geometric state feed the original planar contact-angle closure.
//!
//! Wet/interface cells use their actual source-bound caloric temperature.
//! Genuinely dry wall stencil cells retain the declared reference query; they
//! are numerical ghosts, not liquid thermal stores. A zero-mass wet ghost uses
//! the temperature extension already owned by the enthalpy transport.
//!
//! The thermal stage holds the last accepted contact speed while re-evaluating
//! both pinned angle claims. Its angles drive the same interval's native flow.
//! The accepted advected temperature and original backward geometric speed then
//! admit the endpoint history. Sigma remains the admitted native constant.
//! This is a first-order local-cell wetting closure, not resolved interface
//! temperature, Marangoni physics, calibrated pinning or surface-energy closure.

use super::*;
use crate::d3q19::freesurface3::PlanarWall3;
use crate::material::free_surface::wetting::{
    CoxVoinovControls, CoxVoinovLaw, FreeSurfaceWettingControls, FreeSurfaceWettingSample,
    WettingThermalPoint,
};
use fs_matdb::InterfaceSystemCard;
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Couple the original Cox boundary to an explicitly sourced outer gas.
    /// Its actual mu_gas divided by the local mu_liquid(T) supplies the two-phase
    /// wedge ratio. The two phases retain separate temperatures and receipts;
    /// the bulk gas is not evaluated at the liquid temperature.
    ///
    /// Ordinary stepping holds this explicit gas snapshot. Existing finite-gas
    /// conduction and solved-interface entry points instead consume their SAME
    /// evolving parcel at each thermal endpoint and native flow preparation.
    /// Original gas transport claims stay pinned while admitted composition,
    /// temperature and pressure follow that owner. No second gas momentum or
    /// heat model is created. Microscopic interface-angle claims, constant sigma,
    /// Newtonian creeping-wedge assumptions and matching scales remain explicit.
    pub fn with_cox_voinov_gas_wetting(
        mut self,
        controls: CoxVoinovControls,
        gas: &ResolvedThermochemicalGasState,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        if self.rheology.is_some() {
            return Err(Error::Invalid("Cox wetting requires a Newtonian liquid"));
        }
        let binding = self.wetting.take()
            .ok_or(Error::Invalid("Cox wetting requires a sourced static hysteresis binding"))?
            .with_cox_voinov_gas(controls, gas, &self.flow.scaling, gate)?;
        let temperatures = temperatures(&self.states, self.law.temperature_axis(), gate)?;
        let binding = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        Ok(self)
    }

    /// The outer gas used by accepted contact-line samples. Its source T/p and
    /// composition remain distinct from each local liquid/interface query.
    #[must_use]
    pub fn wetting_outer_gas(&self) -> Option<&ResolvedThermochemicalGasState> {
        self.wetting.as_ref().and_then(WettingBinding::outer_gas)
    }

    /// Bind one ordered substrate/liquid interface to a stationary planar wall.
    /// The same native contact geometry, speed estimate and hysteresis owner is
    /// shared with isothermal flow. Source selection occurs once at the original
    /// complete liquid query; subsequent T/speed evaluations pin both claims.
    /// Other substrate/environment conditions remain explicitly fixed context.
    ///
    /// All ordinary heat, molecular heat, conduction, finite-wall and finite-gas
    /// entry points retain this binding. Thermal substeps and the final native
    /// candidate must remain in source support before any participant commits.
    /// No capillary or contact-angle change is converted into invented heat.
    pub fn with_dynamic_planar_wetting(
        mut self,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        selection: MaterialPropertySelection,
        controls: FreeSurfaceWettingControls,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        if self.thermocapillary.is_some() {
            return Err(Error::Invalid("flat-interface thermocapillarity cannot impose a separate contact-angle geometry"));
        }
        if self.wetting.is_some() {
            return Err(Error::Invalid("thermal liquid already owns a planar wetting history"));
        }
        let binding = WettingBinding::try_new(
            &mut self.flow, wall, card, selection, controls, gate,
        )?;
        let temperatures = temperatures(&self.states, self.law.temperature_axis(), gate)?;
        let binding = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        Ok(self)
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Add a hydrodynamic Cox law to an existing STATIC hysteresis binding.
    /// The actual local source viscosity follows the enthalpy trajectory; sigma
    /// remains the native constant. Shear-dependent liquids and already dynamic
    /// angle curves refuse because the Cox wedges require Newtonian fluids.
    /// A caller-declared finite outer/liquid viscosity ratio remains fixed as
    /// the liquid heats: that explicitly models proportional outer viscosity,
    /// not a resolved or independently sourced gas-temperature trajectory.
    /// Apparent angles alter momentum, without manufacturing capillary heat.
    pub fn with_cox_voinov_wetting(
        mut self,
        law: CoxVoinovLaw,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        if self.rheology.is_some() {
            return Err(Error::Invalid("Cox wetting requires a Newtonian liquid"));
        }
        let binding = self
            .wetting
            .take()
            .ok_or(Error::Invalid(
                "Cox wetting requires a sourced static hysteresis binding",
            ))?
            .with_cox_voinov(law, &self.flow.scaling, gate)?;
        let temperatures = temperatures(&self.states, self.law.temperature_axis(), gate)?;
        let binding = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        Ok(self)
    }

    /// Explicit physical matching scales and limits of the active Cox model.
    #[must_use]
    pub fn cox_voinov_wetting_law(&self) -> Option<&CoxVoinovLaw> {
        self.wetting
            .as_ref()
            .and_then(WettingBinding::cox_voinov_law)
    }

    /// Original selected interface bundle, including both angle claims.
    #[must_use]
    pub fn wetting_source(&self) -> Option<&ResolvedInterfaceStatePoint> {
        self.wetting.as_ref().map(WettingBinding::source)
    }

    /// Accepted local temperature/speed receipts and next-step angle field.
    /// Dry numerical wall ghosts retain the reference temperature context.
    #[must_use]
    pub fn wetting_samples(&self) -> Option<&BTreeMap<[usize; 3], FreeSurfaceWettingSample>> {
        self.wetting.as_ref().map(WettingBinding::samples)
    }

    pub(super) fn prepare_wetting(
        &self,
        lattice: &mut FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        outer_gas: Option<&ResolvedThermochemicalGasState>,
        gate: &CancelGate,
    ) -> Result<Option<WettingBinding<'m>>, Error> {
        let Some(binding) = &self.wetting else { return Ok(None); };
        let prepared = outer_gas.map(|gas| binding.at_prepared_outer_gas(gas));
        let binding = prepared.as_ref().unwrap_or(binding);
        let temperatures = temperatures(states, self.law.temperature_axis(), gate)?;
        let candidate = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        candidate.install(lattice, gate)?;
        Ok(Some(candidate))
    }
}

pub(super) fn temperature(state: &ResolvedLiquidState, axis: &str) -> Result<f64, Error> {
    state.material().query_point().iter()
        .find(|(name, _)| name == axis)
        .map(|(_, value)| *value)
        .filter(|value| value.is_finite() && *value > 0.0)
        .ok_or(Error::Invalid("wetting requires the actual liquid absolute-temperature coordinate"))
}

pub(super) fn temperatures(
    states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
    axis: &str,
    gate: &CancelGate,
) -> Result<BTreeMap<[usize; 3], WettingThermalPoint>, Error> {
    let mut field = BTreeMap::new();
    for (&location, state) in states {
        poll(gate)?;
        if let Some(position) = FreeSurface3::transport_coordinates(location) {
            field.insert(position, WettingThermalPoint {
                temperature_k: temperature(state, axis)?,
                dynamic_viscosity_pa_s: state.dynamic_viscosity_pa_s(),
            });
        }
    }
    Ok(field)
}
