//! Explicit ideal-monatomic material sources for network evolution checks.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Volume};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::fusion::ELECTRON_MASS_KG;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn monatomic(name: &str, mass: f64, maximum_temperature: f64) -> Nasa9StandardStateModelV1 {
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip([0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1e5,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 1e4, maximum_temperature),
        sources: vec![],
        provenance: Provenance {
            source: "Manufactured ideal-monatomic sensible calorics with CODATA 2022 bare masses; no measured plasma claim".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * AVOGADRO),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1e5),
            ElementalReferenceIdV1::new("light-ion-sensible-test").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

pub(super) fn plasma(
    density: [f64; 5],
    volume: f64,
    temperature_kev: [f64; 2],
) -> (LightIonFusionPlasmaV1, FrozenIonizedGasStateV1) {
    plasma_with_helium3_limit(density, volume, temperature_kev, 2e9)
}

pub(super) fn plasma_with_helium3_limit(
    density: [f64; 5],
    volume: f64,
    temperature_kev: [f64; 2],
    helium3_maximum_temperature: f64,
) -> (LightIonFusionPlasmaV1, FrozenIonizedGasStateV1) {
    let names = ["p+", "D+", "T+", "He3++", "He4++", "e-"];
    let ids: Vec<_> = names.map(id).into();
    let masses = [
        NuclearSpeciesV1::Proton.mass_kg(),
        NuclearSpeciesV1::Deuterium.mass_kg(),
        NuclearSpeciesV1::Tritium.mass_kg(),
        NuclearSpeciesV1::Helium3.mass_kg(),
        NuclearSpeciesV1::Helium4.mass_kg(),
        ELECTRON_MASS_KG,
    ];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 1, 0, 0, 0], vec![0, 0, 0, 1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![1, 1, 1, 2, 2, -1]).unwrap(),
        names
            .into_iter()
            .zip(masses)
            .enumerate()
            .map(|(i, (name, mass))| {
                monatomic(
                    name,
                    mass,
                    if i == 3 {
                        helium3_maximum_temperature
                    } else {
                        2e9
                    },
                )
            })
            .collect(),
    )
    .unwrap();
    let model = LightIonFusionPlasmaV1::new(
        bank,
        id(names[0]),
        id(names[1]),
        id(names[2]),
        id(names[3]),
        id(names[4]),
    )
    .unwrap();
    let mut amounts = [Amount::new(0.0); 6];
    for (index, number_density) in model.indices().into_iter().zip(density) {
        amounts[index] = Amount::new(number_density * volume / AVOGADRO);
    }
    amounts[model.electron_index()] = Amount::new(
        (density[0] + density[1] + density[2] + 2.0 * density[3] + 2.0 * density[4]) * volume
            / AVOGADRO,
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
