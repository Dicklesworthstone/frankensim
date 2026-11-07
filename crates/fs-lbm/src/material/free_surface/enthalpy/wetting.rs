//! Thermal and geometric state feed the original planar contact-angle closure.
//!
//! Wet/interface cells use their actual source-bound caloric temperature.
//! Genuinely dry wall stencil cells retain the declared reference query; they
//! are numerical ghosts, not liquid thermal stores. A zero-mass wet ghost uses
//! the temperature extension already owned by the enthalpy transport.
//!
//! The thermal stage holds the last accepted contact speed while re-evaluating
//! the pinned equilibrium source or hysteresis claims. Its angles drive the
//! same interval's native flow.
//! The accepted advected temperature and original backward geometric speed then
//! admit the endpoint history. Equilibrium and dynamic adhesion compose with
//! uniform sigma(T); purely equilibrium wetting also admits the finite-slope
//! height graph's existing Marangoni owner. The graph boundary angle is installed
//! before its normal, area and stress are reconstructed. This local-cell closure
//! does not resolve interface
//! temperature, calibrated pinning or a surface-energy balance.

use super::*;
use crate::d3q19::freesurface3::PlanarWall3;
use crate::material::free_surface::wetting::{
    CoxVoinovControls, CoxVoinovLaw, EquilibriumAngleOrigin, EquilibriumWettingRateLaw,
    FreeSurfaceWettingControls, FreeSurfaceWettingSample, MolecularKineticLaw, WettingThermalPoint,
};
use fs_matdb::InterfaceSystemCard;
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Evolve an explicitly equilibrium wall using its local liquid temperature.
    ///
    /// Choose a sourced equilibrium angle or Young--Dupre from reversible
    /// adhesion work and the SAME liquid/gas tension retained by the
    /// native thermal law. The original ordered interface and selected claim
    /// stay pinned. Every heat, conduction and flow candidate refreshes this
    /// response before publication; a source-domain exit or loss of partial
    /// wetting rolls back the complete thermal/flow step.
    ///
    /// This smooth, homogeneous stationary-wall closure supplies no pinning,
    /// contact-line friction or dynamic-angle law until explicitly extended by
    /// [`Self::with_dynamic_equilibrium_wetting`]. It cannot be combined with
    /// a separate hysteresis, planar or local thermocapillary owner.
    /// The uniform thermal-capillarity owner can evolve sigma(T), provided all
    /// actual wet cells agree exactly at every admitted stage. A height graph
    /// instead admits one planar SIDE wall: every occupied interface cell in
    /// a wall column must give the same angle. Other columns may differ; the
    /// two end columns meeting another side wall must remain exactly neutral
    /// because a nonneutral two-wall corner law is not resolved. Its
    /// angle supplies the shared graph normal/area and tangential-stress metric;
    /// native pressure retains its original same-angle VOF ghost curvature.
    /// This first-cell boundary approximation needs spatial refinement. Changing
    /// adhesion alters native capillary momentum, but does not invent surface
    /// heat or establish a closed total surface/mechanical energy balance.
    /// The local liquid-cell temperature approximates the interface temperature;
    /// no separately solved solid/liquid interface temperature is inferred.
    pub fn with_equilibrium_planar_wetting(
        mut self,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        origin: EquilibriumAngleOrigin,
        selection: MaterialPropertySelection,
        maximum_wall_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let uniform_tension = self.uniform_thermal_capillarity();
        let graph_tension = self.height_graph_thermocapillary_surface().is_some();
        if self.thermocapillary.is_some() && !(uniform_tension || graph_tension)
            || self.wetting.is_some() || graph_tension && self.surface_temperature_forces
        {
            return Err(Error::Invalid(
                "equilibrium wetting needs constant, uniform or bulk-temperature height-graph tension",
            ));
        }
        let mut graph_materials = if graph_tension {
            let masses = self.flow.lattice.mass_snapshot();
            let old_h = self.h.clone();
            // Storage-free ghosts use the SAME thermal extension at initial
            // attachment and later heated endpoints; they are not extra input
            // temperatures that may independently select a wall angle.
            extend_empty_ghosts(&masses, &mut self.h, &old_h, true)?;
            Some(self.resolve_liquid_states(&self.flow.lattice, &self.h, &masses, gate)?)
        } else {
            None
        };
        let binding = WettingBinding::try_equilibrium(
            &mut self.flow, wall, card, origin, selection, self.law.temperature_axis(),
            uniform_tension || graph_tension,
            maximum_wall_cells, gate,
        )?;
        let states = graph_materials.as_ref().map_or(&self.states, |resolved| &resolved.states);
        let temperatures = temperatures(states, self.law.temperature_axis(), gate)?;
        let binding = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        if let Some(mut resolved) = graph_materials.take() {
            // Attaching an angle changes the boundary patch normal and area.
            // Reconstruct its stress now, without advancing heat or momentum.
            self.resolve_surface_fields(&self.flow.lattice, &mut resolved, &mut None, gate)?;
            if let Some(tensions) = resolved.surface_tensions {
                self.flow.lattice.rebind_surface_tensions(tensions, gate)?;
            }
            if let Some(forces) = resolved.forces {
                self.flow.lattice.rebind_material_force_field(forces, gate)?;
            }
            self.states = resolved.states;
            self.rheology_samples = resolved.rheology;
        }
        Ok(self)
    }

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
    /// heat model is created. Microscopic interface-angle claims, native sigma,
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
    /// The uniform thermal-capillarity owner may evolve sigma(T); all actual
    /// wet-cell tensions must still agree exactly at every admitted stage.
    /// Other thermocapillary geometry owners cannot share this moving boundary.
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
        if self.thermocapillary.is_some() && !self.uniform_thermal_capillarity() {
            return Err(Error::Invalid(
                "dynamic planar wetting needs constant or spatially uniform surface tension",
            ));
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
    /// Explicitly drive moving contact lines from the original equilibrium
    /// angle or reversible Young--Dupre work using one chosen Cox or molecular
    /// rate law. The same local T and sigma feed equilibrium and rate response.
    ///
    /// This homogeneous smooth-wall model has no inferred hysteresis interval
    /// or pinning threshold. Original equilibrium claims retain their selected
    /// identities and do not acquire a fictitious speed coordinate. The native
    /// accepted VOF history supplies speed; unresolved geometry, source exits,
    /// calibration exits and cancellation refuse the complete candidate.
    /// Kinetic parameters and Cox cutoffs remain explicit model inputs. Local
    /// molecular loss is diagnostic and is not silently injected as heat.
    /// Constant or uniform thermal surface tension is required, as for the
    /// existing dynamic planar closure; height-graph coupling is not admitted.
    pub fn with_dynamic_equilibrium_wetting(
        mut self,
        law: EquilibriumWettingRateLaw,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        if self.thermocapillary.is_some() && !self.uniform_thermal_capillarity() {
            return Err(Error::Invalid(
                "dynamic equilibrium wetting needs constant or spatially uniform surface tension",
            ));
        }
        if matches!(&law, EquilibriumWettingRateLaw::CoxVoinov(_)) && self.rheology.is_some() {
            return Err(Error::Invalid("Cox wetting requires a Newtonian liquid"));
        }
        let binding = self
            .wetting
            .take()
            .ok_or(Error::Invalid(
                "dynamic equilibrium wetting requires an original equilibrium source",
            ))?
            .with_equilibrium_rate_law(
                law,
                self.law.temperature_axis(),
                &self.flow.scaling,
                gate,
            )?;
        let temperatures = temperatures(&self.states, self.law.temperature_axis(), gate)?;
        let binding = binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        Ok(self)
    }

    /// Add a hydrodynamic Cox law to an existing STATIC hysteresis binding.
    /// The actual local source viscosity and native surface tension follow the
    /// enthalpy trajectory. Surface tension is constant unless the uniform
    /// thermal-capillarity owner is active. Shear-dependent liquids and already
    /// dynamic angle curves refuse because the Cox wedges require Newtonian fluids.
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

    /// Add molecular jump kinetics to an existing STATIC hysteresis binding.
    ///
    /// Every thermal preparation and endpoint uses its actual local liquid
    /// temperature and the same tension as native pressure. The caller's jump
    /// frequency and length stay constant only inside their declared calibration
    /// interval. Dry numerical ghosts retain the original reference context.
    /// Cox, another molecular model and measured rate curves cannot be combined
    /// with this alternative closure. Model dissipation is reported per contact-
    /// line length; this entry does not create a new heat source.
    pub fn with_molecular_kinetic_wetting(
        mut self,
        law: MolecularKineticLaw,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let binding = self
            .wetting
            .take()
            .ok_or(Error::Invalid(
                "molecular wetting requires a sourced static hysteresis binding",
            ))?
            .with_molecular_kinetic(
                law,
                self.law.temperature_axis(),
                &self.flow.scaling,
                gate,
            )?;
        let temperatures = temperatures(&self.states, self.law.temperature_axis(), gate)?;
        let binding =
            binding.at_temperatures(self.law.temperature_axis(), &temperatures, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        poll(gate)?;
        self.wetting = Some(binding);
        Ok(self)
    }

    /// Calibration and admitted domain of the optional molecular rate law.
    #[must_use]
    pub fn molecular_kinetic_wetting_law(&self) -> Option<&MolecularKineticLaw> {
        self.wetting
            .as_ref()
            .and_then(WettingBinding::molecular_kinetic_law)
    }

    /// Explicit physical matching scales and limits of the active Cox model.
    #[must_use]
    pub fn cox_voinov_wetting_law(&self) -> Option<&CoxVoinovLaw> {
        self.wetting
            .as_ref()
            .and_then(WettingBinding::cox_voinov_law)
    }

    /// Original selected equilibrium or hysteresis interface source bundle.
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

    pub(super) fn endpoint_wetting(
        &self,
        binding: Option<&WettingBinding<'m>>,
        lattice: &mut FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        gate: &CancelGate,
    ) -> Result<Option<WettingBinding<'m>>, Error> {
        let Some(binding) = binding else { return Ok(None); };
        let temperatures = temperatures(states, self.law.temperature_axis(), gate)?;
        Ok(Some(binding.at_endpoint(
            lattice, &self.flow.scaling,
            Some((self.law.temperature_axis(), &temperatures)), gate,
        )?))
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
                surface_tension_n_m: state.material().property(SURFACE_TENSION_PROPERTY)
                    .map(|source| source.value_si()),
            });
        }
    }
    Ok(field)
}
