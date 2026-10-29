//! Original low-temperature NASA9 N2/O2 calorics, with an explicit dry-air mix.
//!
//! McBride, Zehe and Gordon, NASA TP-2002-211556; original coefficients from
//! https://github.com/Cantera/cantera/blob/v3.2.0/data/airNASA9.yaml
//! Git blob 5975d1303819332c1bd77a783d9873b66f65c952. Only the original
//! 200--1000 K region is used. There is no fit, seam adjustment or extrapolation.

use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use std::error::Error;

pub fn source() -> Result<(MixtureGas, Vec<f64>), Box<dyn Error>> {
    let models = vec![
        species(
            "N2",
            0.028014,
            [
                2.210371497e4,
                -381.846182,
                6.08273836,
                -8.53091441e-3,
                1.384646189e-5,
                -9.62579362e-9,
                2.519705809e-12,
                710.846086,
                -10.76003744,
            ],
        )?,
        species(
            "O2",
            0.031998,
            [
                -3.42556342e4,
                484.700097,
                1.119010961,
                4.29388924e-3,
                -6.83630052e-7,
                -2.0233727e-9,
                1.039040018e-12,
                -3391.45487,
                18.4969947,
            ],
        )?,
    ];
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("N")?, ElementId::new("O")?],
            ids.clone(),
            vec![vec![2, 0], vec![0, 2]],
        )?,
        ChargeVector::new(ids, vec![0, 0])?,
        models,
    )?;
    let total_mass = 0.79 * 0.028014 + 0.21 * 0.031998;
    let fractions = bank
        .models()
        .iter()
        .map(|m| {
            if m.species().as_str() == "N2" {
                0.79 * 0.028014 / total_mass
            } else {
                0.21 * 0.031998 / total_mass
            }
        })
        .collect();
    Ok((
        MixtureGas::new(
            bank,
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
            || false,
        )?,
        fractions,
    ))
}

fn species(
    name: &str,
    mass: f64,
    coefficients: [f64; 9],
) -> Result<Nasa9StandardStateModelV1, Box<dyn Error>> {
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 101325.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()), law_version: 1,
        parameters, state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: vec![],
        provenance: Provenance {
            source: "NASA TP-2002-211556; Cantera v3.2.0 data/airNASA9.yaml, original 200--1000 K coefficients; blob 5975d1303819332c1bd77a783d9873b66f65c952".into(),
            license: "NASA source data distributed with Cantera; original attribution retained".into(),
            artifact: None,
        },
    };
    Ok(Nasa9StandardStateModelV1::new(
        SpeciesId::new(name)?,
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(101325.0),
            ElementalReferenceIdV1::new("nasa-glenn-reference-elements")?,
        )?,
        vec![Nasa9RegionV1::from_card(card)?],
    )?)
}
