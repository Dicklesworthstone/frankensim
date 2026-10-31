use super::*;
mod background;
mod transient;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Volume};
use fs_thermochem::electron_collision::{
    AVOGADRO_PER_MOL as NA, BOLTZMANN_J_PER_K as KB, ELECTRON_MASS_KG as ME, ELECTRONVOLT_J as EV,
    MaxwellianRateControlsV1, nist,
};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::AtomicStatisticalWeightsV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

const MP: f64 = 1.672_621_925_95e-27;
const CHI: f64 = 13.598_434_599_702 * EV * NA;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn species(name: &str, mass: f64, g: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    let h = 6.626_070_15e-34;
    let a8 = 2.5
        + (g * (2.0 * core::f64::consts::PI * mass * KB / (h * h)).powf(1.5) * KB / 100000.0).ln();
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip([
        0.0,
        0.0,
        2.5,
        0.0,
        0.0,
        0.0,
        0.0,
        formation / R,
        a8,
    ]) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass * NA),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("atomic-translational-test").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
            state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 1000.0, 50000.0), sources: vec![],
            provenance: Provenance { source: "analytic Sackur-Tetrode ground-state H/H+/e with NIST ionization gap; no fitted discharge claim".into(),
                license: "analytic test fixture".into(), artifact: None },
        }).unwrap()]).unwrap()
}
fn initial(
    x: f64,
    te: f64,
    th: f64,
    elemental_shift: f64,
    charge_shift: f64,
) -> (FrozenIonizedGasStateV1, AtomicIonizationV1) {
    let ids = vec![id("H"), id("H+"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        vec![
            species("H", MP + ME, 2.0, elemental_shift),
            species("H+", MP, 1.0, elemental_shift + CHI + charge_shift),
            species("e", ME, 2.0, -charge_shift),
        ],
    )
    .unwrap();
    let channel = AtomicIonizationV1::new(
        bank.clone(),
        ReactionId::new("H-impact").unwrap(),
        id("H+"),
        ElectronImpactRateV1::new(
            id("e"),
            nist::hydrogen(id("H")),
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-9,
                maximum_subdivisions: 4096,
            },
            0.1,
        )
        .unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let state = FrozenIonizedGasStateV1::new(
        &bank,
        &[Amount::new(1.0 - x), Amount::new(x), Amount::new(x)],
        Volume::new(1.0),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap();
    (state, channel)
}
fn controls() -> AtomicIonizationSolve {
    let mut controls =
        AtomicIonizationSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    controls.relative_tolerance = 1e-6;
    controls.energy_tolerance_j = 1e-7;
    controls
}
fn run(
    state: &FrozenIonizedGasStateV1,
    channel: &AtomicIonizationV1,
    duration: f64,
) -> AtomicIonizationResult {
    advance_atomic_ionization(
        state,
        channel,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(duration),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn atomic_ionization_cools_electrons_to_independent_energy_and_partial_balance_endpoint() {
    let (state, channel) = initial(0.05, 30000.0, 4000.0, 0.0, 0.0);
    let result = run(&state, &channel, 1e-3);
    // Independent 70-digit Decimal solution of the statistical Saha factor
    // and the electron energy equation; no collision-rate/time-step code used.
    let expected_x = 0.058_078_371_038_461_755;
    let expected_te = 11194.144_850_406_924;
    let end = result.state();
    assert!((end.amounts()[2].value() - expected_x).abs() < 2e-8);
    assert!((end.electron_temperature().value() - expected_te).abs() < 0.05);
    assert_eq!(end.heavy_temperature(), state.heavy_temperature());
    assert!((result.balance().entropy_change_j_k - 0.442_967_127_757_084_75).abs() < 1e-6);
    assert!((end.amounts()[0].value() + end.amounts()[1].value() - 1.0).abs() < 2e-15);
    assert_eq!(end.amounts()[1], end.amounts()[2]);
    assert!((end.mass().value() / state.mass().value() - 1.0).abs() < 2e-15);
    assert!(result.balance().energy_residual_j.abs() < 1e-7);
    assert!(result.balance().minimum_step_entropy_excess_j_k >= -controls().entropy_tolerance_j_k);
    assert!(result.balance().maximum_normalized_local_error <= 1.0);
}

#[test]
fn atomic_three_body_recombination_returns_chemical_energy_to_electrons() {
    let (state, channel) = initial(0.2, 6000.0, 8000.0, 0.0, 0.0);
    let result = run(&state, &channel, 1e-3);
    let end = result.state();
    assert!(result.reaction_extent().value() < 0.0);
    assert!(result.balance().chemical_storage_change_j < 0.0);
    assert!((end.amounts()[2].value() - 0.187_674_426_037_161_4).abs() < 2e-8);
    assert!((end.electron_temperature().value() - 13303.236_842_581_762).abs() < 0.05);
    assert_eq!(end.heavy_temperature(), state.heavy_temperature());
    assert!((result.balance().entropy_change_j_k - 0.624_134_275_063_291_4).abs() < 1e-6);
    assert!((end.internal_energy().value() - state.internal_energy().value()).abs() < 1e-7);
}

#[test]
fn held_pool_powers_close_source_energy_with_elemental_and_charge_reference_shifts() {
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        for (pe, ph) in [(2e7, 3e6), (-1e6, -2e6)] {
            let (state, channel) = initial(0.05, 20000.0, 5000.0, elemental, charge);
            let duration = 2e-7;
            let result = advance_atomic_ionization(
                &state,
                &channel,
                Power::new(pe),
                Power::new(ph),
                Time::new(duration),
                controls(),
                || false,
            )
            .unwrap();
            let end = result.state();
            let extent = result.reaction_extent().value();
            let electron_sensible = 1.5
                * R
                * (end.amounts()[2].value() * end.electron_temperature().value()
                    - state.amounts()[2].value() * state.electron_temperature().value());
            let heavy_sensible =
                1.5 * R * (end.heavy_temperature().value() - state.heavy_temperature().value());
            assert!(
                (electron_sensible + channel.ionization_energy_j_mol() * extent - pe * duration)
                    .abs()
                    < 2e-7
            );
            assert!((heavy_sensible - ph * duration).abs() < 1e-7);
            assert!(result.balance().energy_residual_j.abs() < 1e-7);
            assert!((result.balance().electron_deposited_energy_j - pe * duration).abs() < 1e-15);
            assert!(
                result.balance().minimum_step_entropy_excess_j_k
                    >= -controls().entropy_tolerance_j_k
            );
        }
    }
}

#[test]
fn zero_electrons_stay_zero_and_cannot_receive_electron_work() {
    let (state, channel) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let result = advance_atomic_ionization(
        &state,
        &channel,
        Power::new(0.0),
        Power::new(1000.0),
        Time::new(1e-3),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.state().amounts(), state.amounts());
    assert_eq!(
        result.state().electron_temperature(),
        state.electron_temperature()
    );
    assert!(result.state().heavy_temperature().value() > 4000.0);
    assert!(matches!(
        advance_atomic_ionization(
            &state,
            &channel,
            Power::new(1.0),
            Power::new(0.0),
            Time::new(1e-6),
            controls(),
            || false
        ),
        Err(AtomicIonizationError::Invalid(_))
    ));
}

#[test]
fn cancellation_source_limits_and_step_budget_refuse_without_mutating_input() {
    let (state, channel) = initial(0.05, 30000.0, 4000.0, 0.0, 0.0);
    let original = state.clone();
    assert_eq!(
        advance_atomic_ionization(
            &state,
            &channel,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-3),
            controls(),
            || true
        ),
        Err(AtomicIonizationError::Cancelled)
    );
    let mut calls = 0;
    assert_eq!(
        advance_atomic_ionization(
            &state,
            &channel,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-3),
            controls(),
            || {
                calls += 1;
                calls > 200
            }
        ),
        Err(AtomicIonizationError::Cancelled)
    );
    let mut bounded = controls();
    bounded.maximum_steps = 1;
    assert!(matches!(
        advance_atomic_ionization(
            &state,
            &channel,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-3),
            bounded,
            || false
        ),
        Err(AtomicIonizationError::StepBudget)
    ));
    let mut outside = controls();
    outside.minimum_electron_temperature_k = 31000.0;
    assert!(matches!(
        advance_atomic_ionization(
            &state,
            &channel,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-3),
            outside,
            || false
        ),
        Err(AtomicIonizationError::OutsideDomain)
    ));
    assert_eq!(state, original);
    let zero = run(&state, &channel, 0.0);
    assert_eq!(zero.state(), &state);
    assert_eq!(zero.balance().accepted_steps, 0);
}

#[test]
fn recombination_at_the_closed_electron_source_minimum_heats_into_the_domain() {
    let (state, channel) = initial(0.05, 1000.0, 4000.0, 0.0, 0.0);
    assert_eq!(
        state.electron_temperature().value(),
        controls().minimum_electron_temperature_k
    );
    // The physical reverse reaction heats this boundary state. Reconstructing
    // the unchanged state as Ue/Cve previously rounded below 1000 K, so every
    // retry refused before even evaluating the supported interior solution.
    let result = run(&state, &channel, 1e-10);
    assert!(result.reaction_extent().value() < 0.0);
    assert!(result.state().electron_temperature().value() > 1000.0);
    assert_eq!(
        result.state().heavy_temperature(),
        state.heavy_temperature()
    );
    assert_eq!(result.state().amounts()[1], result.state().amounts()[2]);
    assert!(result.balance().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(result.balance().entropy_change_j_k >= -controls().entropy_tolerance_j_k);
}
