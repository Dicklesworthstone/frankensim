//! Run `cargo run -p fs-gas --features experimental-euler --example reacting_oxygen -- 80`.
//! Emits a one-metre periodic O/O2 reacting-flow field as CSV.
//!
//! Species calorics and the isolated `2 O + M <=> O2 + M` reaction are taken
//! from GRI-Mech 3.0 as distributed in Cantera v3.2.0's `data/gri30.yaml`:
//! <https://github.com/Cantera/cantera/blob/v3.2.0/data/gri30.yaml>.
//! This is a source-data demonstration of an estimated inviscid solver, not a
//! complete air/combustion mechanism or experimental validation.

use fs_evidence::ValidityDomain;
use fs_gas::chemistry::IsochoricKineticsSolve;
use fs_gas::multicomponent::reactive::{ReactingEuler1d, ReactiveStepControl};
use fs_gas::multicomponent::{Boundaries, CaloricDomain, MixtureEuler1d, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, PressureDependenceV1,
    ReactionProgressRateQuantityV1, ThirdBodyEfficienciesV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    Nasa9RegionV1, Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1, SpeciesId,
    StandardStateConventionV1, StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    required_nasa9_parameter_dims,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let count = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(80);
    if !(2..=4096).contains(&count) {
        return Err("example cell count must be 2..=4096".into());
    }
    let mechanism = oxygen_mechanism()?;
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(1000.0), Temperature::new(3500.0)),
        || false,
    )?;
    let states = (0..count)
        .map(|i| {
            let x = (i as f64 + 0.5) / count as f64;
            let atomic_fraction = 0.02 + 0.01 * (2.0 * core::f64::consts::PI * x).sin();
            let r_mix = R * (atomic_fraction / 0.015999 + (1.0 - atomic_fraction) / 0.031998);
            let density = 101325.0 / (r_mix * 1800.0);
            gas.conserved(
                &[density * atomic_fraction, density * (1.0 - atomic_fraction)],
                100.0,
                101325.0,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let flow = MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.01,
        states,
        Boundaries::Periodic,
        || false,
    )?;
    let mut simulation = ReactingEuler1d::new(
        flow,
        mechanism,
        IsochoricKineticsSolve {
            minimum_temperature_k: 1000.0,
            maximum_temperature_k: 3500.0,
            relative_tolerance: 1e-6,
            energy_tolerance_j: 1e-9,
            absolute_amount_tolerance_mol: 1e-15,
            initial_step_s: 1e-6,
            maximum_step_s: 1e-3,
            ..IsochoricKineticsSolve::default()
        },
        || false,
    )?;
    let initial = simulation.flow().totals();
    let mut steps = 0;
    let mut chemistry_steps = 0;
    let end = 1e-3;
    while simulation.flow().time().value() < end {
        if steps >= 100000 {
            return Err("example step budget exhausted".into());
        }
        let report = simulation.step(
            Time::new(end - simulation.flow().time().value()),
            ReactiveStepControl::default(),
            || false,
        )?;
        steps += 1;
        chemistry_steps += report.chemistry.accepted_steps;
    }
    let result = simulation.flow().totals();
    let references = simulation.flow().gas().energy_references();
    let energy_residual = (result.energy - initial.energy)
        + (result.species_masses[0] - initial.species_masses[0]) * references[0]
        + (result.species_masses[1] - initial.species_masses[1]) * references[1];
    eprintln!(
        "GRI3 isolated O/O2 reacting flow: cells={count},steps={steps},chemistry_substeps={chemistry_steps},time_s={:.12e},mass_residual_kg={:.6e},chemical_total_energy_residual_j={energy_residual:.6e},atomic_oxygen_initial_kg={:.6e},atomic_oxygen_final_kg={:.6e}",
        simulation.flow().time().value(),
        result.mass() - initial.mass(),
        initial.species_masses[0],
        result.species_masses[0]
    );
    println!(
        "x_m,density_kg_m3,velocity_m_s,pressure_pa,temperature_k,mass_fraction_O,mass_fraction_O2"
    );
    for i in 0..count {
        let q = simulation.flow().cell(i).expect("valid index");
        let p = simulation.flow().gas().primitive(q, || false)?;
        println!(
            "{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            (i as f64 + 0.5) / count as f64,
            p.density,
            p.velocity,
            p.pressure,
            p.temperature,
            q.partial_densities[0] / p.density,
            q.partial_densities[1] / p.density
        );
    }
    Ok(())
}

fn oxygen_mechanism() -> Result<IdealGasKineticsV1, Box<dyn std::error::Error>> {
    let models = vec![
        oxygen_source(
            "O",
            0.015999,
            [
                2.56942078,
                -8.59741137e-5,
                4.19484589e-8,
                -1.00177799e-11,
                1.22833691e-15,
                2.92175791e4,
                4.78433864,
            ],
        )?,
        oxygen_source(
            "O2",
            0.031998,
            [
                3.28253784,
                1.48308754e-3,
                -7.57966669e-7,
                2.09470555e-10,
                -2.16717794e-14,
                -1088.45772,
                5.45323129,
            ],
        )?,
    ];
    let species = models
        .iter()
        .map(|s| s.species().clone())
        .collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("O")?],
            species.clone(),
            vec![vec![1, 2]],
        )?,
        ChargeVector::new(species, vec![0, 0])?,
        models,
    )?;
    // GRI: A=1.2e17 cm6/(mol2 s), b=-1, Ea=0. Convert cm6 to m6:
    // k(1000K)=(1.2e17 * 1e-12)/1000=120 m6/(mol2 s).
    // At c*=1 mol/m3, the API's normalized progress scale is 120 mol/(m3 s).
    // GRI's listed efficiency overrides concern absent species; O and O2 use1.
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("gri30-reaction-1")?,
        vec![(SpeciesId::new("O")?, 2)],
        vec![(SpeciesId::new("O2")?, 1)],
        ArrheniusRateV1::new(
            Temperature::new(1000.0),
            ReactionProgressRateQuantityV1::new(120.0),
            -1.0,
            MolarEnergyQuantityV1::new(0.0),
        )?,
    )?
    .with_pressure_dependence(PressureDependenceV1::ThirdBody(
        ThirdBodyEfficienciesV1::default(),
    ));
    Ok(IdealGasKineticsV1::new(
        bank,
        "gri30-cantera-3.2.0-oxygen-reaction-1",
        1,
        vec![reaction],
        1e-12,
    )?)
}

fn oxygen_source(
    name: &str,
    mass: f64,
    nasa7: [f64; 7],
) -> Result<Nasa9StandardStateModelV1, Box<dyn std::error::Error>> {
    // Exact NASA7 -> NASA9 embedding: the T^-2 and T^-1 Cp coefficients are0.
    // Retain only the original 1000..3500K region, avoiding any invented
    // reconciliation of the source's independently fitted region endpoints.
    let a = [
        0.0, 0.0, nasa7[0], nasa7[1], nasa7[2], nasa7[3], nasa7[4], nasa7[5], nasa7[6],
    ];
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    // Cantera's NASA standard-state reference pressure defaults to 1 atm.
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 101325.0,
            dims: Pressure::DIMS,
        },
    );
    let card=ConstitutiveModelCard {
        law:LawId("nasa9-standard-state".into()),law_version:1,parameters,state_schema_version:0,
        initial_state:InitialStatePolicy::ZeroInternalState,
        validity:ValidityDomain::unconstrained().with("T",1000.0,3500.0),sources:vec![],
        provenance:Provenance {source:"GRI-Mech 3.0, Cantera v3.2.0 data/gri30.yaml, original high-temperature NASA7 coefficients".into(),
            license:"GRI-Mech source data distributed with Cantera; original attribution retained".into(),artifact:None},
    };
    Ok(Nasa9StandardStateModelV1::new(
        SpeciesId::new(name)?,
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(101325.0),
            ElementalReferenceIdV1::new("gri30-original-reference-elements")?,
        )?,
        vec![Nasa9RegionV1::from_card(card)?],
    )?)
}
