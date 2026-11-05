//! Analytic monatomic source cards for the constructed D-T example and tests.
//! This is an explicit ideal-gas model, not a measured plasma calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Pressure, Temperature, Volume};
use fs_thermochem::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use fs_thermochem::fusion::plasma::DtFusionPlasmaV1;
use fs_thermochem::fusion::{
    BOLTZMANN_J_PER_K, ELECTRON_MASS_KG, KELVIN_PER_KEV, NuclearSpeciesV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

const NA: f64 = 6.022_140_76e23;
const PLANCK: f64 = 6.626_070_15e-34;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).expect("literal species")
}

pub(super) fn monatomic(name: &str, mass: f64, spin_weight: f64) -> Nasa9StandardStateModelV1 {
    let p0 = 100_000.0;
    // Sackur-Tetrode translational entropy with an explicit ground-state
    // spin weight; zero formation energy keeps nuclear rest energy outside U.
    let entropy = 2.5
        + fs_math::det::ln(spin_weight * BOLTZMANN_J_PER_K / p0)
        + 1.5
            * fs_math::det::ln(
                2.0 * core::f64::consts::PI * mass * BOLTZMANN_J_PER_K / (PLANCK * PLANCK),
            );
    let coefficients = [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, entropy];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: p0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
        state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 1.0e4, 2.0e9),
        sources: vec![],
        provenance: Provenance {
            source: "Ideal monatomic statistical mechanics, with CODATA 2022 bare masses: https://physics.nist.gov/cuu/Constants/Table/allascii.txt; constructed example, no measured caloric calibration".into(),
            license: "analytic model and public physical constants".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * NA),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(p0),
            ElementalReferenceIdV1::new("dt-ideal-kinetic-energy").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

/// D/T/He4 number densities [1/m3], volume [m3], and ion/electron kT [keV].
pub fn plasma(
    density: [f64; 3],
    volume: f64,
    temperature_kev: [f64; 2],
) -> (DtFusionPlasmaV1, FrozenIonizedGasStateV1) {
    let ids = vec![id("D+"), id("T+"), id("He4++"), id("e-")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![1, 1, 2, -1]).unwrap(),
        vec![
            monatomic("D+", NuclearSpeciesV1::Deuterium.mass_kg(), 3.0),
            monatomic("T+", NuclearSpeciesV1::Tritium.mass_kg(), 2.0),
            monatomic("He4++", NuclearSpeciesV1::Helium4.mass_kg(), 1.0),
            monatomic("e-", ELECTRON_MASS_KG, 2.0),
        ],
    )
    .unwrap();
    let model = DtFusionPlasmaV1::new(bank, id("D+"), id("T+"), id("He4++")).unwrap();
    let mut amounts = [Amount::new(0.0); 4];
    for (index, number_density) in model.indices().into_iter().zip(density) {
        amounts[index] = Amount::new(number_density * volume / NA);
    }
    amounts[model.electron_index()] =
        Amount::new((density[0] + density[1] + 2.0 * density[2]) * volume / NA);
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
