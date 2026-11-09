//! Source molecular/turbulent heat and constituent enthalpy share one gas endpoint.

use super::adaptive::{TurbulentGasSpeciesAdvance, advance_prepared};
use super::{
    CancelGate, ChtError, ClosedGasParcel, FvFlow, GasEnergySolve, MaterialAdvectionStep,
    PreparedTurbulentGasMixing, SpeciesNetworkError, SpeciesNetworkSolve, TurbulentGasSpeciesStep,
    VoxelDomain, check, material_error, poll,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkAdvanceControls, SpeciesNetworkHeatTransfer, transport::GasTransportSolve,
};
use std::collections::BTreeMap;

struct PreparedTurbulentHeat<'a> {
    domain: &'a VoxelDomain,
    heat_eddy_viscosity_m2_s: &'a [f64],
    turbulent_prandtl: f64,
    include_molecular: bool,
}

fn invalid(reason: &'static str) -> SpeciesNetworkError {
    SpeciesNetworkError::Invalid(reason)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl<'a> PreparedTurbulentHeat<'a> {
    fn new(
        flow: &'a FvFlow,
        domain: &'a VoxelDomain,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        include_molecular: bool,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if !positive(flow.turbulent_prandtl)
            || flow.heat_eddy_viscosity_m2_s.len() != domain.cell_count()
        {
            return Err(material_error(invalid(
                "turbulent gas heat requires the complete field and a positive Prandtl number",
            )));
        }
        for (cell, &eddy) in flow.heat_eddy_viscosity_m2_s.iter().enumerate() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            if !eddy.is_finite() || eddy < 0.0 || (!domain.is_fluid(cell) && eddy != 0.0) {
                return Err(material_error(invalid(
                    "thermal eddy viscosity must be finite, nonnegative and zero in solids",
                )));
            }
            let diffusivity = eddy / flow.turbulent_prandtl;
            if !diffusivity.is_finite() || (eddy > 0.0 && diffusivity == 0.0) {
                return Err(material_error(invalid(
                    "positive turbulent thermal diffusivity is unrepresentable",
                )));
            }
        }
        for (id, parcel) in parcels {
            poll(gate)?;
            let endpoint = endpoints.get(id).ok_or_else(|| {
                material_error(invalid("every original gas needs its caloric controls"))
            })?;
            let temperature = parcel.state().state().temperature;
            if !(endpoint.minimum_temperature_k..=endpoint.maximum_temperature_k)
                .contains(&temperature)
            {
                return Err(material_error(invalid(
                    "original gas temperature lies outside its endpoint bracket",
                )));
            }
        }
        Ok(Self {
            domain,
            heat_eddy_viscosity_m2_s: &flow.heat_eddy_viscosity_m2_s,
            turbulent_prandtl: flow.turbulent_prandtl,
            include_molecular,
        })
    }

    fn sample(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        duration_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<SpeciesNetworkHeatTransfer>, SpeciesNetworkError> {
        check(cancelled)?;
        if !positive(duration_s) {
            return Err(invalid(
                "turbulent heat interval must be finite and positive",
            ));
        }
        let mut conductivity = BTreeMap::new();
        for (&id, parcel) in parcels {
            check(cancelled)?;
            let state = parcel.state().state();
            let eddy = self.heat_eddy_viscosity_m2_s[id as usize];
            let turbulent = if eddy == 0.0 {
                0.0
            } else {
                let k = state.density * state.specific_heat_cp * (eddy / self.turbulent_prandtl);
                if !positive(k) {
                    return Err(invalid(
                        "source-derived turbulent conductivity is unrepresentable",
                    ));
                }
                k
            };
            // Molecular and unresolved transport act in parallel INSIDE each
            // material half-cell. Put their sum through the original series
            // face resistance; summing two separately harmonic face laws would
            // impose different fictitious interface temperatures.
            let k = if self.include_molecular {
                let molecular = state.thermal_conductivity;
                let total = molecular + turbulent;
                if !positive(molecular) || !positive(total) {
                    return Err(invalid(
                        "source molecular plus turbulent conductivity is unrepresentable",
                    ));
                }
                total
            } else {
                turbulent
            };
            conductivity.insert(id, k);
        }
        // The shared source owner already admits at most 256 fluid parcels.
        // Only positive directions are visited, so every face occurs once.
        let mut heat = Vec::new();
        heat.try_reserve_exact(3 * parcels.len())
            .map_err(|_| invalid("unable to allocate bounded turbulent gas heat streams"))?;
        for (&left, &k_left) in &conductivity {
            check(cancelled)?;
            if k_left == 0.0 {
                continue;
            }
            for side in [1, 3, 5] {
                let Some(right) = self
                    .domain
                    .neighbor(left as usize, side)
                    .filter(|&cell| self.domain.is_fluid(cell))
                else {
                    continue;
                };
                let right_id = right as u64;
                let k_right = conductivity[&right_id];
                if k_right == 0.0 {
                    continue;
                }
                let difference = parcels[&left].state().state().temperature
                    - parcels[&right_id].state().state().temperature;
                if difference == 0.0 {
                    continue;
                }
                let axis = side / 2;
                let left_resistance = (0.5 * self.domain.widths(left as usize)[axis]) / k_left;
                let right_resistance = (0.5 * self.domain.widths(right)[axis]) / k_right;
                if !positive(left_resistance) || !positive(right_resistance) {
                    return Err(invalid(
                        "turbulent heat half-cell resistance is unrepresentable",
                    ));
                }
                let largest = left_resistance.max(right_resistance);
                let conductance = (self.domain.face_area(left as usize, axis) / largest)
                    / (1.0 + left_resistance.min(right_resistance) / largest);
                let energy = (conductance * duration_s) * difference;
                if !positive(conductance) || !energy.is_finite() || energy == 0.0 {
                    return Err(invalid("nonzero turbulent heat stream is unrepresentable"));
                }
                heat.push(SpeciesNetworkHeatTransfer {
                    parcels: [left, right_id],
                    heat_left_to_right_j: energy,
                });
            }
        }
        check(cancelled)?;
        Ok(heat)
    }
}

impl FvFlow {
    /// Advance turbulent heat and species from one unchanged material snapshot.
    ///
    /// Species use the original momentum `nu_t / Sc_t` exactly as in
    /// [`Self::diffuse_turbulent_gas_species`]. Heat independently uses this
    /// flow's thermal eddy field and retained positive turbulent Prandtl number:
    /// `k_t = rho cp nu_t,heat / Pr_t`. Density, heat capacity and temperature
    /// come from each SAME original gas state. The original LVEL/SGS thermal
    /// policy therefore remains explicit; it never selects species mixing.
    ///
    /// Graded half-cell thermal resistances give one antisymmetric face heat
    /// `Q = dt K (T_left - T_right)`. This reduced heat is passed alongside all
    /// species streams to one original material reconstruction. Partial species
    /// enthalpy and formation references are carried exactly once. A uniform
    /// composition can still exchange heat; uniform temperature gives zero
    /// reduced heat even while composition mixes.
    ///
    /// Solid, exterior and zero-conductivity faces are adiabatic for this
    /// component. Molecular conduction and wall heat exchange are not supplied.
    /// [`Self::diffuse_turbulent_gas_transport_with_molecular_heat`] explicitly
    /// includes the original gas conductivity in the same interior face law.
    /// This explicit snapshot is first order in time and has no unconditional
    /// stability or finite-step entropy claim. Use a resolved interval, or
    /// [`Self::advance_turbulent_gas_transport`] for source refresh, local error
    /// control and the original closed-window physical audits. A failed source
    /// solve publishes neither temperature nor composition.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_turbulent_gas_transport(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesStep, ChtError> {
        self.gas_transport_snapshot(
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            duration_s,
            endpoints,
            controls,
            false,
            gate,
        )
    }

    /// Combine original molecular conduction, turbulent heat and species.
    ///
    /// Every cell supplies `k = k_source(T,Y) + rho cp(T,Y) nu_t,heat / Pr_t`
    /// from the SAME original material snapshot. The total cell conductivities
    /// enter the graded half-cell series resistance before heat is sampled.
    /// Conducting fluid faces therefore retain molecular heat transfer when
    /// either thermal eddy viscosity vanishes. At zero SGS this recovers the
    /// original source's molecular conduction, with variable caloric inversion.
    ///
    /// All heat and turbulent species streams enter one original gas endpoint;
    /// partial species enthalpy is carried once. The heat ledger records the
    /// total paired heat and must not be combined with another molecular solve
    /// on the same interval. No molecular species diffusion, Soret/Dufour law,
    /// wall heat transfer or external reservoir is supplied. This retains the
    /// existing explicit first-order snapshot and its fixed source/work bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_turbulent_gas_transport_with_molecular_heat(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesStep, ChtError> {
        self.gas_transport_snapshot(
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            duration_s,
            endpoints,
            controls,
            true,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn gas_transport_snapshot(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        include_molecular: bool,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesStep, ChtError> {
        let prepared = PreparedTurbulentGasMixing::new(
            self,
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            endpoints,
            controls,
            gate,
        )?;
        let heat =
            PreparedTurbulentHeat::new(self, domain, parcels, endpoints, include_molecular, gate)?;
        let mut cancelled = || gate.is_requested();
        let heat = heat
            .sample(parcels, duration_s, &mut cancelled)
            .map_err(material_error)?;
        prepared
            .step_candidate(
                parcels,
                duration_s,
                endpoints,
                controls,
                &heat,
                &mut cancelled,
            )
            .map_err(material_error)
    }

    /// Complete a turbulent heat/species interval with source-aware local error control.
    ///
    /// This reuses the material network's existing coarse-versus-two-half-step
    /// controller, accepted-only mass/enthalpy/heat ledger, and independent
    /// original-to-final mass, energy and entropy checks. Rejected probes spend
    /// work but cannot contribute a heat stream or alter an accepted inventory.
    /// Every probe recomputes conductivity from its current source density and
    /// `cp(T,Y)`, then samples heat and donor partial enthalpy at that same state.
    ///
    /// The original momentum and thermal eddy fields, grid and their physical
    /// epoch are held for this interval; recomputing the turbulent flow remains
    /// the momentum owner's responsibility. The retained two-half-step path is
    /// still first order, with explicit local species/temperature tolerances and
    /// finite work/step limits. This is not a global-error or turbulence-model
    /// validation certificate. All original parcel/stream/source bounds apply.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_turbulent_gas_transport(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: GasTransportSolve,
        temporal: SpeciesNetworkAdvanceControls,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesAdvance, ChtError> {
        self.gas_transport_window(
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            duration_s,
            endpoints,
            physical,
            temporal,
            false,
            gate,
        )
    }

    /// Advance combined molecular/turbulent heat and turbulent species.
    ///
    /// This uses [`Self::diffuse_turbulent_gas_transport_with_molecular_heat`]'s
    /// total half-cell conductivity in the existing adaptive material owner.
    /// Every trial refreshes original source conductivity, density, heat
    /// capacity, composition and donor partial enthalpy together; source
    /// conductivity is never frozen at the beginning of the whole interval.
    /// Only accepted fine steps enter the signed total-heat ledger. Original
    /// source failures, thermal-bracket recovery, physical budgets and
    /// cancellation keep [`Self::advance_turbulent_gas_transport`]'s semantics.
    ///
    /// Molecular conduction persists as the thermal SGS field vanishes.
    /// Molecular species diffusion and external wall heat require their own
    /// declared models. The held flow epoch and first-order locally controlled
    /// integration have the same limits as the turbulent-only entry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_turbulent_gas_transport_with_molecular_heat(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: GasTransportSolve,
        temporal: SpeciesNetworkAdvanceControls,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesAdvance, ChtError> {
        self.gas_transport_window(
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            duration_s,
            endpoints,
            physical,
            temporal,
            true,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn gas_transport_window(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: GasTransportSolve,
        temporal: SpeciesNetworkAdvanceControls,
        include_molecular: bool,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesAdvance, ChtError> {
        let prepared = PreparedTurbulentGasMixing::new(
            self,
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            endpoints,
            physical.exchange,
            gate,
        )?;
        let heat =
            PreparedTurbulentHeat::new(self, domain, parcels, endpoints, include_molecular, gate)?;
        advance_prepared(
            &prepared,
            parcels,
            duration_s,
            endpoints,
            physical,
            temporal,
            |state, dt, cancelled| heat.sample(state, dt, cancelled),
            gate,
        )
    }
}
