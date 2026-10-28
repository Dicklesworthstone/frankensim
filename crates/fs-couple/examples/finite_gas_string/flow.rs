//! Accepted two-parcel observations; thermodynamics stays in the native owners.
use super::*;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use std::collections::BTreeSet;

pub(super) const COLUMNS: [&str; 12] = [
    "surrounding_mass_kg",
    "surrounding_internal_energy_j",
    "reservoir_mass_kg",
    "reservoir_internal_energy_j",
    "reservoir_temperature_k",
    "reservoir_pressure_pa",
    "pair_internal_energy_j",
    "mass_flow_into_surroundings_kg_s",
    "cumulative_mass_into_surroundings_kg",
    "cumulative_enthalpy_into_surroundings_j",
    "surrounding_storage_change_j",
    "reservoir_storage_change_j",
];

pub(super) struct History {
    pub species: Vec<String>,
    pub records: Vec<Vec<f64>>,
    initial: [(f64, f64); 2],
    mass_into_kg: f64,
    enthalpy_into_j: f64,
}

impl History {
    pub fn new(gas: &IsochoricGasParcel, reservoir: &ClosedGasParcel, samples: usize) -> Self {
        let species = gas
            .state()
            .thermodynamics()
            .receipt()
            .species()
            .iter()
            .chain(reservoir.state().thermodynamics().receipt().species())
            .map(|id| id.as_str().to_owned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Self {
            species,
            records: Vec::with_capacity(samples),
            initial: [
                (gas.mass_kg(), gas.state().specific_internal_energy_j_kg()),
                (
                    reservoir.mass_kg(),
                    reservoir.state().specific_internal_energy_j_kg(),
                ),
            ],
            mass_into_kg: 0.0,
            enthalpy_into_j: 0.0,
        }
    }

    /// Measure extensive endpoints independently of the flow step's report.
    /// Difference specific energies first to retain small changes against a
    /// valid large thermochemical formation-energy reference.
    pub fn energy_changes(
        &self,
        gas: &IsochoricGasParcel,
        reservoir: &ClosedGasParcel,
    ) -> [f64; 2] {
        let change = |mass: f64, u: f64, (m0, u0): (f64, f64)| mass * (u - u0) + (mass - m0) * u0;
        [
            change(
                gas.mass_kg(),
                gas.state().specific_internal_energy_j_kg(),
                self.initial[0],
            ),
            change(
                reservoir.mass_kg(),
                reservoir.state().specific_internal_energy_j_kg(),
                self.initial[1],
            ),
        ]
    }

    pub fn reference_energy_scale(
        &self,
        gas: &IsochoricGasParcel,
        reservoir: &ClosedGasParcel,
    ) -> f64 {
        self.initial.iter().map(|(m, u)| m * u.abs()).sum::<f64>()
            + gas.mass_kg() * gas.state().specific_internal_energy_j_kg().abs()
            + reservoir.mass_kg() * reservoir.state().specific_internal_energy_j_kg().abs()
    }

    pub fn record(
        &mut self,
        gas: &IsochoricGasParcel,
        reservoir: &ClosedGasParcel,
        mass_flow_into_kg_s: f64,
        mass_into_kg: f64,
        enthalpy_into_j: f64,
    ) -> Result<()> {
        self.mass_into_kg += mass_into_kg;
        self.enthalpy_into_j += enthalpy_into_j;
        let gas_energy = gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
        let reservoir_energy =
            reservoir.mass_kg() * reservoir.state().specific_internal_energy_j_kg();
        let changes = self.energy_changes(gas, reservoir);
        let mut row = vec![
            gas.mass_kg(),
            gas_energy,
            reservoir.mass_kg(),
            reservoir_energy,
            reservoir.state().state().temperature,
            reservoir.state().state().pressure,
            gas_energy + reservoir_energy,
            mass_flow_into_kg_s,
            self.mass_into_kg,
            self.enthalpy_into_j,
            changes[0],
            changes[1],
        ];
        for species in &self.species {
            row.push(mole_fraction(gas.state(), species));
            row.push(mole_fraction(reservoir.state(), species));
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err("nonfinite accepted two-parcel observation".into());
        }
        self.records.push(row);
        Ok(())
    }
}

fn mole_fraction(gas: &ResolvedThermochemicalGasState, species: &str) -> f64 {
    let receipt = gas.thermodynamics().receipt();
    receipt
        .species()
        .iter()
        .position(|id| id.as_str() == species)
        .map_or(0.0, |i| f64::from_bits(receipt.mole_fraction_bits()[i]))
}
