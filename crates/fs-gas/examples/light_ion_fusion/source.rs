//! Six original material sources using the existing example's monatomic law.

// Share the actual statistical-mechanical source builder; only its four-ion
// example entry point is unused in this six-species program.
#[allow(dead_code)]
#[path = "../dt_fusion/source.rs"]
mod calorics;

use fs_qty::{Amount, Temperature, Volume};
use fs_thermochem::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use fs_thermochem::fusion::plasma::LightIonFusionPlasmaV1;
use fs_thermochem::fusion::{ELECTRON_MASS_KG, KELVIN_PER_KEV, NuclearSpeciesV1};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, SpeciesId};

const NA: f64 = 6.022_140_76e23;

/// Authored number densities [p,D,T,He3,He4] in m^-3, V in m^3, and kT in keV.
pub fn plasma(
    density: [f64; 5],
    volume: f64,
    temperature_kev: [f64; 2],
) -> (LightIonFusionPlasmaV1, FrozenIonizedGasStateV1) {
    let names = ["p+", "D+", "T+", "He3++", "He4++", "e-"];
    let ids = names.map(|name| SpeciesId::new(name).unwrap());
    let masses = [
        NuclearSpeciesV1::Proton.mass_kg(),
        NuclearSpeciesV1::Deuterium.mass_kg(),
        NuclearSpeciesV1::Tritium.mass_kg(),
        NuclearSpeciesV1::Helium3.mass_kg(),
        NuclearSpeciesV1::Helium4.mass_kg(),
        ELECTRON_MASS_KG,
    ];
    let spin_weights = [2.0, 3.0, 2.0, 2.0, 1.0, 2.0];
    let models = (0..6)
        .map(|i| calorics::monatomic(names[i], masses[i], spin_weights[i]))
        .collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.to_vec(),
            vec![vec![1, 1, 1, 0, 0, 0], vec![0, 0, 0, 1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids.to_vec(), vec![1, 1, 1, 2, 2, -1]).unwrap(),
        models,
    )
    .unwrap();
    let model = LightIonFusionPlasmaV1::new(
        bank,
        ids[0].clone(),
        ids[1].clone(),
        ids[2].clone(),
        ids[3].clone(),
        ids[4].clone(),
    )
    .unwrap();
    let mut amounts = [Amount::new(0.0); 6];
    for (i, n) in model.indices().into_iter().zip(density) {
        amounts[i] = Amount::new(n * volume / NA);
    }
    amounts[model.electron_index()] = Amount::new(
        (density[0] + density[1] + density[2] + 2.0 * (density[3] + density[4])) * volume / NA,
    );
    let state = FrozenIonizedGasStateV1::new(
        model.thermodynamics(),
        &amounts,
        Volume::new(volume),
        Temperature::new(temperature_kev[0] * KELVIN_PER_KEV),
        Temperature::new(temperature_kev[1] * KELVIN_PER_KEV),
        64.0 * f64::EPSILON,
        || false,
    )
    .unwrap();
    (model, state)
}
