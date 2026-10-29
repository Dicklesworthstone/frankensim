//! Five original NASA CEA gases for the reacting-compartment example.
//!
//! NASA TP-2002-211556, ideal-gas standard pressure 1 bar.
//! Source: https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp
//! Git blob: 5ddeb5b213701d0e8b79dbda7b48ede14e31404a.
//! Original 200--1000 and 1000--6000 K regions and source molar masses.
//! The original Cp coefficients are unchanged. The existing explicit caloric
//! continuity owner anchors the low-temperature region at 300 K and applies
//! bounded constant h/s offsets to the high-temperature region. The example
//! reports those actual offsets; they are not empirical uncertainty estimates.

use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    MolarThermalQuantityV1, Nasa9RegionV1, Nasa9StandardStateModelV1, ReactionId,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};
use std::error::Error;

pub const NAMES: [&str; 5] = ["CH4", "CO2", "H2O", "N2", "O2"];
pub const MASSES: [f64; 5] = [0.01604246, 0.0440095, 0.01801528, 0.0280134, 0.0319988];
pub const STOICHIOMETRY: [f64; 5] = [-1.0, 1.0, 2.0, 0.0, -2.0];
pub const REACTION: &str = "methane-complete-oxidation-mixing-closure";
pub const SOURCE_COMMIT: &str = "2f79a647e85737742fb08a127da00d1ea1d069da";
pub const SOURCE_BLOB: &str = "5ddeb5b213701d0e8b79dbda7b48ede14e31404a";

const COEFFICIENTS: [[[f64; 9]; 2]; 5] = [
    [
        [
            -1.766850998e5,
            2.786181020e3,
            -1.202577850e1,
            3.917619290e-2,
            -3.619054430e-5,
            2.026853043e-8,
            -4.976705490e-12,
            -2.331314360e4,
            8.904322750e1,
        ],
        [
            3.730042760e6,
            -1.383501485e4,
            2.049107091e1,
            -1.961974759e-3,
            4.727313040e-7,
            -3.728814690e-11,
            1.623737207e-15,
            7.532066910e4,
            -1.219124889e2,
        ],
    ],
    [
        [
            4.943650540e4,
            -6.264116010e2,
            5.301725240,
            2.503813816e-3,
            -2.127308728e-7,
            -7.689988780e-10,
            2.849677801e-13,
            -4.528198460e4,
            -7.048279440,
        ],
        [
            1.176962419e5,
            -1.788791477e3,
            8.291523190,
            -9.223156780e-5,
            4.863676880e-9,
            -1.891053312e-12,
            6.330036590e-16,
            -3.908350590e4,
            -2.652669281e1,
        ],
    ],
    [
        [
            -3.947960830e4,
            5.755731020e2,
            9.317826530e-1,
            7.222712860e-3,
            -7.342557370e-6,
            4.955043490e-9,
            -1.336933246e-12,
            -3.303974310e4,
            1.724205775e1,
        ],
        [
            1.034972096e6,
            -2.412698562e3,
            4.646110780,
            2.291998307e-3,
            -6.836830480e-7,
            9.426468930e-11,
            -4.822380530e-15,
            -1.384286509e4,
            -7.978148510,
        ],
    ],
    [
        [
            2.210371497e4,
            -3.818461820e2,
            6.082738360,
            -8.530914410e-3,
            1.384646189e-5,
            -9.625793620e-9,
            2.519705809e-12,
            7.108460860e2,
            -1.076003744e1,
        ],
        [
            5.877124060e5,
            -2.239249073e3,
            6.066949220,
            -6.139685500e-4,
            1.491806679e-7,
            -1.923105485e-11,
            1.061954386e-15,
            1.283210415e4,
            -1.586640027e1,
        ],
    ],
    [
        [
            -3.425563420e4,
            4.847000970e2,
            1.119010961,
            4.293889240e-3,
            -6.836300520e-7,
            -2.023372700e-9,
            1.039040018e-12,
            -3.391454870e3,
            1.849699470e1,
        ],
        [
            -1.037939022e6,
            2.344830282e3,
            1.819732036,
            1.267847582e-3,
            -2.188067988e-7,
            2.053719572e-11,
            -8.193467050e-16,
            -1.689010929e4,
            1.738716506e1,
        ],
    ],
];

pub fn source() -> Result<(MixtureGas, IdealGasKineticsV1, Vec<f64>), Box<dyn Error>> {
    let models = (0..NAMES.len()).map(model).collect::<Result<Vec<_>, _>>()?;
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            ["C", "H", "N", "O"]
                .iter()
                .map(|s| ElementId::new(*s))
                .collect::<Result<Vec<_>, _>>()?,
            ids.clone(),
            vec![
                vec![1, 1, 0, 0, 0],
                vec![4, 0, 2, 0, 0],
                vec![0, 0, 0, 2, 0],
                vec![0, 2, 1, 0, 2],
            ],
        )?,
        ChargeVector::new(ids, vec![0; NAMES.len()])?,
        models,
    )?;
    // This adapter uses stoichiometry, never Arrhenius kinetics. The required
    // rate-law carrier is explicitly a unit normalization, not methane kinetics.
    let unused_law = ArrheniusRateV1::new(
        Temperature::new(300.0),
        ReactionProgressRateQuantityV1::new(1.0),
        0.0,
        MolarEnergyQuantityV1::new(0.0),
    )?;
    let reaction = ElementaryReactionV1::new(
        ReactionId::new(REACTION)?,
        vec![(SpeciesId::new("CH4")?, 1), (SpeciesId::new("O2")?, 2)],
        vec![(SpeciesId::new("CO2")?, 1), (SpeciesId::new("H2O")?, 2)],
        unused_law,
    )?;
    let mechanism = IdealGasKineticsV1::new(
        bank.clone(),
        "nasa-cea-methane-inventory-closure",
        1,
        vec![reaction],
        32.0 * f64::EPSILON,
    )?;
    let gas = MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(3500.0)),
        || false,
    )?;
    let air_mass = 0.79 * MASSES[3] + 0.21 * MASSES[4];
    let fractions = vec![
        0.0,
        0.0,
        0.0,
        0.79 * MASSES[3] / air_mass,
        0.21 * MASSES[4] / air_mass,
    ];
    Ok((gas, mechanism, fractions))
}

fn model(index: usize) -> Result<Nasa9StandardStateModelV1, Box<dyn Error>> {
    let mut regions = Vec::with_capacity(2);
    for (region, (lo, hi)) in [(200.0, 1000.0), (1000.0, 6000.0)].into_iter().enumerate() {
        let mut parameters = std::collections::BTreeMap::new();
        for ((key, dims), value) in required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(COEFFICIENTS[index][region])
        {
            parameters.insert(key.into(), LawParameter { value, dims });
        }
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: 100000.0,
                dims: Pressure::DIMS,
            },
        );
        regions.push(Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId("nasa9-standard-state".into()),
            law_version: 1, parameters, state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", lo, hi),
            sources: vec![],
            provenance: Provenance {
                source: format!("NASA TP-2002-211556; nasa/cea {SOURCE_COMMIT}; thermo.inp blob {SOURCE_BLOB}; original {} {lo}--{hi} K", NAMES[index]),
                license: "NASA thermodynamic data; original source attribution retained".into(),
                artifact: None,
            },
        })?);
    }
    Ok(Nasa9StandardStateModelV1::new(
        SpeciesId::new(NAMES[index])?,
        MolarMass::new(MASSES[index]),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("nasa-glenn-reference-elements")?,
        )?,
        regions,
    )?
    .with_caloric_continuity(
        Temperature::new(300.0),
        MolarEnergyQuantityV1::new(0.01),
        MolarThermalQuantityV1::new(1e-5),
    )?)
}
