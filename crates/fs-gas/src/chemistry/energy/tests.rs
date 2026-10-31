use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, Area, MolarMass, Pressure, Volume};
use fs_thermochem::electron_collision::{BOLTZMANN_J_PER_K, ElectronCrossSectionPointV1};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn species(name: &str, mass: f64, formation: f64) -> Nasa9StandardStateModelV1 {
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
        0.0,
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
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 1000.0, 30000.0),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured monatomic calorics and collision tests; no measured plasma claim"
                .into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * AVOGADRO_PER_MOL),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("energy-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn state(te: f64, th: f64, fraction: f64, formation_shift: f64) -> FrozenIonizedGasStateV1 {
    let proton = 1.672_621_925_95e-27;
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
            species("H", proton + ELECTRON_MASS_KG, formation_shift),
            species("H+", proton, formation_shift + 1.3e6),
            species("e", ELECTRON_MASS_KG, 0.0),
        ],
    )
    .unwrap();
    FrozenIonizedGasStateV1::new(
        &bank,
        &[
            Amount::new(1.0 - fraction),
            Amount::new(fraction),
            Amount::new(fraction),
        ],
        Volume::new(1.0),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap()
}
fn collision(sigma: f64, maximum_ev: f64) -> ElectronCrossSectionV1 {
    ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "manufactured constant momentum cross section".into(),
            license: "fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(sigma),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(maximum_ev * 1.602_176_634e-19),
                cross_section: Area::new(sigma),
            },
        ],
        Area::new(sigma),
    )
    .unwrap()
}
fn controls() -> ElectronHeavyEnergySolve {
    let mut c = ElectronHeavyEnergySolve::new(
        Temperature::new(1000.0),
        Temperature::new(30000.0),
        0.01,
        0.1,
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    );
    c.initial_step_s = 1e-4;
    c.maximum_step_s = 1e-3;
    c.energy_tolerance_j = 1e-8;
    c
}

// For constant sigma, G(Te)=g*sqrt(Te); this exact nonlinear two-pool solution
// is stronger than freezing a collision frequency to manufacture an exponential.
fn analytic_constant_sigma(
    initial: &FrozenIonizedGasStateV1,
    sigma: f64,
    duration: f64,
) -> (f64, f64) {
    let ch = initial.heavy_heat_capacity_j_k();
    let ce = initial.electron_heat_capacity_j_k();
    let te = initial.electron_temperature().value();
    let th = initial.heavy_temperature().value();
    let eq = (ch * th + ce * te) / (ch + ce);
    let mean_speed = (8.0 * BOLTZMANN_J_PER_K / (core::f64::consts::PI * ELECTRON_MASS_KG)).sqrt();
    let mass = initial.thermodynamics().models()[0].molar_mass().value() / AVOGADRO_PER_MOL;
    let electrons = initial.amounts()[2].value() * AVOGADRO_PER_MOL;
    let g = 3.0
        * BOLTZMANN_J_PER_K
        * electrons
        * initial.number_density_per_m3(0).unwrap()
        * (ELECTRON_MASS_KG / mass)
        * (4.0 / 3.0)
        * sigma
        * mean_speed;
    let y = (te / eq).sqrt();
    let z = (y - 1.0) / (y + 1.0) * (-g * (1.0 / ce + 1.0 / ch) * eq.sqrt() * duration).exp();
    let end_e = eq * ((1.0 + z) / (1.0 - z)).powi(2);
    let end_h = th + (ce / ch) * (te - end_e);
    (end_e, end_h)
}

#[test]
fn sourced_two_temperature_relaxation_converges_to_independent_nonlinear_solution() {
    let initial = state(12000.0, 4000.0, 1e-4, 0.0);
    let sigma = 1e-24;
    let duration = 0.002;
    let expected = analytic_constant_sigma(&initial, sigma, duration);
    let mut errors = Vec::new();
    for tolerance in [1e-3, 1e-5, 1e-7] {
        let mut c = controls();
        c.relative_tolerance = tolerance;
        let result = advance_frozen_energy(
            &initial,
            &[collision(sigma, 200.0)],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(duration),
            c,
            || false,
        )
        .unwrap();
        let te = result.state().electron_temperature().value();
        let th = result.state().heavy_temperature().value();
        errors.push((te - expected.0).abs());
        assert!(te < 12000.0 && te > th && th > 4000.0);
        assert_eq!(result.state().amounts(), initial.amounts());
        assert!(result.balance().electron_to_heavy_energy_j > 0.0);
        assert!(result.balance().electron_energy_residual_j.abs() <= c.energy_tolerance_j);
        assert!(result.balance().heavy_energy_residual_j.abs() <= c.energy_tolerance_j);
        assert!(result.balance().maximum_normalized_local_error <= 1.0);
    }
    assert!(
        errors[1] < errors[0] * 0.3 && errors[2] < errors[1] * 0.3,
        "{errors:?}"
    );
    assert!(errors[2] / expected.0 < 3e-4, "{errors:?}");
}

#[test]
fn independent_external_powers_close_each_source_pool_without_formation_cancellation() {
    for shift in [0.0, 1e12] {
        let initial = state(12000.0, 4000.0, 1e-4, shift);
        let c = controls();
        let result = advance_frozen_energy(
            &initial,
            &[collision(1e-24, 200.0)],
            Power::new(2.0),
            Power::new(4.0),
            Time::new(0.001),
            c,
            || false,
        )
        .unwrap();
        let b = result.balance();
        assert_eq!(b.electron_deposited_energy_j, 0.002);
        assert_eq!(b.heavy_deposited_energy_j, 0.004);
        let mut changes = [0.0; 2];
        for (i, model) in initial.thermodynamics().models().iter().enumerate() {
            let pool = usize::from(Some(i) == initial.electron_index());
            changes[pool] += initial.amounts()[i].value()
                * model
                    .molar_internal_energy_change(
                        initial.species_temperature(i).unwrap(),
                        result.state().species_temperature(i).unwrap(),
                    )
                    .unwrap()
                    .value();
        }
        assert!((changes[0] - (0.004 + b.electron_to_heavy_energy_j)).abs() < c.energy_tolerance_j);
        assert!((changes[1] - (0.002 - b.electron_to_heavy_energy_j)).abs() < c.energy_tolerance_j);
        assert!((changes[0] + changes[1] - 0.006).abs() < 2.0 * c.energy_tolerance_j);
    }
}

#[test]
fn zero_exchange_follows_original_calorics_and_absent_electrons_refuse_heat() {
    let initial = state(6000.0, 4000.0, 1e-4, 0.0);
    let result = advance_frozen_energy(
        &initial,
        &[collision(0.0, 200.0)],
        Power::new(1.0),
        Power::new(2.0),
        Time::new(0.01),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.exchanged_energy().value(), 0.0);
    assert!(
        (result.state().electron_temperature().value()
            - (6000.0 + 0.01 / initial.electron_heat_capacity_j_k()))
        .abs()
            < 1e-6
    );
    assert!(
        (result.state().heavy_temperature().value()
            - (4000.0 + 0.02 / initial.heavy_heat_capacity_j_k()))
        .abs()
            < 1e-6
    );
    let neutral = state(4000.0, 4000.0, 0.0, 0.0);
    assert!(matches!(
        advance_frozen_energy(
            &neutral,
            &[],
            Power::new(1.0),
            Power::new(0.0),
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(ElectronHeavyEnergyError::Invalid(_))
    ));
    let heated = advance_frozen_energy(
        &neutral,
        &[],
        Power::new(0.0),
        Power::new(2.0),
        Time::new(0.01),
        controls(),
        || false,
    )
    .unwrap();
    assert!(heated.state().heavy_temperature().value() > 4000.0);
    assert_eq!(heated.state().amounts()[2].value(), 0.0);
}

#[test]
fn trace_electron_temperature_is_resolved_below_absolute_parcel_energy_tolerance() {
    let initial = state(6000.0, 4000.0, 1e-24, 0.0);
    let duration = 0.01;
    let electron_power = initial.electron_heat_capacity_j_k() * 1000.0 / duration;
    let c = controls();
    assert!(electron_power * duration < c.energy_tolerance_j * 1e-10);
    let result = advance_frozen_energy(
        &initial,
        &[collision(0.0, 200.0)],
        Power::new(electron_power),
        Power::new(0.0),
        Time::new(duration),
        c,
        || false,
    )
    .unwrap();
    assert!((result.state().electron_temperature().value() - 7000.0).abs() < 1e-6);
    assert_eq!(result.state().heavy_temperature().value(), 4000.0);
    // Relaxation time depends on the neutral density, not on the tiny electron
    // inventory. The same source-driven temperature decay must remain active.
    let trace = state(12000.0, 4000.0, 1e-20, 0.0);
    let expected = analytic_constant_sigma(&trace, 1e-24, 0.002).0;
    let mut c = controls();
    c.relative_tolerance = 1e-7;
    let cooled = advance_frozen_energy(
        &trace,
        &[collision(1e-24, 200.0)],
        Power::new(0.0),
        Power::new(0.0),
        Time::new(0.002),
        c,
        || false,
    )
    .unwrap();
    assert!((cooled.state().electron_temperature().value() - expected).abs() / expected < 3e-4);
}

#[test]
fn missing_or_inadequate_collision_source_and_cancellation_leave_original_state_intact() {
    let initial = state(12000.0, 4000.0, 1e-4, 0.0);
    let original = initial.clone();
    assert!(matches!(
        advance_frozen_energy(
            &initial,
            &[],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(ElectronHeavyEnergyError::MissingSource(_))
    ));
    assert!(matches!(
        advance_frozen_energy(
            &initial,
            &[collision(1e-24, 1.0)],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(ElectronHeavyEnergyError::Collision { .. })
    ));
    let mut calls = 0;
    assert_eq!(
        advance_frozen_energy(
            &initial,
            &[collision(1e-24, 200.0)],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(0.01),
            controls(),
            || {
                calls += 1;
                calls > 40
            }
        ),
        Err(ElectronHeavyEnergyError::Cancelled)
    );
    assert_eq!(initial, original);
    let mut bounded = controls();
    bounded.maximum_steps = 1;
    assert_eq!(
        advance_frozen_energy(
            &initial,
            &[collision(1e-24, 200.0)],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1.0),
            bounded,
            || false
        ),
        Err(ElectronHeavyEnergyError::StepBudget)
    );
}

#[test]
fn slow_but_light_target_does_not_acquire_massive_target_physics() {
    let ids = vec![id("H"), id("H+"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![vec![1, 1, 0]]).unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        vec![species("H", ELECTRON_MASS_KG, 0.0),
             species("H+", 1.672_621_925_95e-27, 0.0),
             species("e", ELECTRON_MASS_KG, 0.0)],
    ).unwrap();
    let initial = FrozenIonizedGasStateV1::new(&bank,
        &[Amount::new(0.9999), Amount::new(1e-4), Amount::new(1e-4)],
        Volume::new(1.0), Temperature::new(1000.0), Temperature::new(30000.0), 1e-12, || false).unwrap();
    let mut c = controls();
    c.maximum_target_thermal_speed_ratio = 0.2;
    assert!((1000.0_f64 / 30000.0).sqrt() < c.maximum_target_thermal_speed_ratio);
    assert_eq!(advance_frozen_energy(&initial, &[collision(1e-24,200.0)], Power::new(0.0),
        Power::new(0.0), Time::new(0.001), c, || false), Err(ElectronHeavyEnergyError::OutsideDomain));
}
