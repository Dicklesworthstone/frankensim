//! Coupled chemistry and thermal exchange against independent source equations.

use super::*;
use crate::chemistry::energy::{ElectronHeavyCollisionSources, ElectronHeavyEnergyError};
use fs_qty::Area;
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
};

fn parcel(
    x: f64,
    te: f64,
    th: f64,
    elemental_shift: f64,
    charge_shift: f64,
) -> (FrozenIonizedGasStateV1, AtomicIonizationV1) {
    let (source, channel) = initial(x, te, th, elemental_shift, charge_shift);
    let state = FrozenIonizedGasStateV1::new(
        source.thermodynamics(),
        source.amounts(),
        Volume::new(1000.0),
        source.heavy_temperature(),
        source.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    (state, channel)
}

fn sources(state: &FrozenIonizedGasStateV1) -> ElectronHeavyCollisionSources {
    // A synthetic constant momentum cross section supplies an analytic neutral
    // Maxwellian moment. It does not claim measured hydrogen transport data.
    let momentum = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "synthetic constant momentum cross section for coupled source oracle".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(1e-20),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(500.0 * EV),
                cross_section: Area::new(1e-20),
            },
        ],
        Area::new(1e-20),
    )
    .unwrap();
    ElectronHeavyCollisionSources {
        momentum_sources: vec![momentum],
        coulomb: Some(
            CoulombModelV1::new(
                state,
                CoulombScreeningV1::StaticAllChargedDebye,
                CoulombLimitsV1 {
                    minimum_coulomb_logarithm: 4.0,
                    minimum_particles_per_debye_sphere: 10.0,
                    maximum_coulomb_coupling: 0.1,
                    maximum_thermal_phase_space_density: 0.1,
                    maximum_thermal_energy_rest_energy_ratio: 0.001,
                    maximum_quantum_to_classical_ratio: 0.5,
                    maximum_relative_drift_thermal_ratio: 0.1,
                },
                || false,
            )
            .unwrap(),
        ),
        maximum_charged_mole_fraction: 1.0,
        maximum_target_thermal_speed_ratio: 0.1,
        rate_controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    }
}

fn coupled_controls() -> AtomicIonizationSolve {
    let mut c = controls();
    c.initial_step_s = 1e-9;
    c.maximum_step_s = 1e-6;
    c
}

fn assert_budgets(result: &AtomicIonizationResult) {
    let initial = result.initial_state();
    let end = result.state();
    let b = result.balance();
    let q = b.electron_to_heavy_energy_j;
    let electron_sensible = 1.5
        * R
        * (end.amounts()[2].value() * end.electron_temperature().value()
            - initial.amounts()[2].value() * initial.electron_temperature().value());
    let heavy_sensible =
        1.5 * R * (end.heavy_temperature().value() - initial.heavy_temperature().value());
    let chemical = result.channel().ionization_energy_j_mol() * result.reaction_extent().value();
    assert!(
        (electron_sensible + chemical + q - b.electron_deposited_energy_j).abs() < 3e-7,
        "electron sensible/chemical/exchange balance: {b:?}"
    );
    assert!((heavy_sensible - q - b.heavy_deposited_energy_j).abs() < 3e-7);
    assert!(
        (b.electron_to_heavy_coulomb_energy_j + b.electron_to_heavy_neutral_energy_j - q).abs()
            < 3e-11
    );
    assert!(b.energy_residual_j.abs() <= result.controls().energy_tolerance_j);
    assert!(b.electron_sensible_energy_residual_j.abs() <= result.controls().energy_tolerance_j);
    assert!(b.heavy_sensible_energy_residual_j.abs() <= result.controls().energy_tolerance_j);
    assert!(b.minimum_step_entropy_excess_j_k >= -result.controls().entropy_tolerance_j_k);
    assert!(b.maximum_normalized_local_error <= 1.0);
    assert!((end.amounts()[0].value() + end.amounts()[1].value() - 1.0).abs() < 2e-15);
    assert_eq!(end.amounts()[1], end.amounts()[2]);
    assert!((end.mass().value() / initial.mass().value() - 1.0).abs() < 3e-15);
    assert!(end.relative_charge_residual() <= result.controls().relative_charge_tolerance);
    assert_eq!(end.volume(), initial.volume());
    assert_eq!(end.thermodynamics(), initial.thermodynamics());
}

#[test]
fn coupled_transient_refines_to_independent_beb_landau_and_neutral_oracle() {
    // Independent Python/SciPy development oracle; no Rust calls or rate
    // surrogate. QUADPACK integrates NIST generalized BEB Eq.1 for H 1s:
    // B=U=13.6057 eV, N=1, Q=.5668, a0=.52918e-10 m. Sources:
    // https://physics.nist.gov/PhysRefData/Ionization/Eqs/latex.html
    // https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I
    // Kc is the Sackur-Tetrode ground-state formula with gH=2,gI=1,ge=2,
    // actual me,mp,mH=me+mp, and chi=13.598434599702 eV. The independent
    // finite-mass Landau moment uses static all-carrier Debye screening and
    // b90=|qe*qi|/(4*pi*epsilon0*mu*S), S=kB*(Te/me+Th/mp):
    // https://farside.ph.utexas.edu/teaching/plasma/lectures/node42.html
    // Neutral km=(4/3)*sigma_mt*sqrt(8*kB*Te/(pi*me)) with synthetic constant
    // sigma_mt=1e-20 m2 and the Banks stationary-heavy energy-transfer factor:
    // https://ntrs.nasa.gov/api/citations/19660020172/downloads/19660020172.pdf
    //
    // DOP853 (rtol=3e-12) and Radau (rtol=1e-12) independently evolve
    // [x,Qc,Qn], using exact changing-population electron/heavy energy budgets.
    // At t=1e-6 s their Te values agree within4e-8 K and x within2e-15.
    // This is a numerical oracle for the declared model, not plasma data.
    let expected = [
        0.050_732_428_248_191_45,
        23_470.789_371_955_572,
        4232.216_808_612_829,
        118.434_841_588_653_46,
        2777.702_120_188_658_7,
        0.696_340_020_483_489_8,
    ];
    let (state, channel) = parcel(0.05, 30000.0, 4000.0, 0.0, 0.0);
    let collision_sources = sources(&state);
    let mut errors = [[0.0; 6]; 3];
    let mut accepted = [0; 3];
    for (i, relative) in [1e-3, 1e-5, 1e-8].into_iter().enumerate() {
        let mut c = coupled_controls();
        c.relative_tolerance = relative;
        let result = advance_atomic_ionization_with_collisions(
            &state,
            &channel,
            &collision_sources,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-6),
            c,
            || false,
        )
        .unwrap();
        let end = result.state();
        let b = result.balance();
        let actual = [
            end.amounts()[2].value(),
            end.electron_temperature().value(),
            end.heavy_temperature().value(),
            b.electron_to_heavy_neutral_energy_j,
            b.electron_to_heavy_coulomb_energy_j,
            b.entropy_change_j_k,
        ];
        errors[i] = core::array::from_fn(|j| (actual[j] - expected[j]).abs());
        accepted[i] = b.accepted_steps;
        assert_budgets(&result);
        assert_eq!(result.elapsed_time(), Time::new(1e-6));
        assert_eq!(result.collision_sources(), Some(&collision_sources));
        assert!(result.reaction_extent().value() > 0.0);
        assert!(b.electron_to_heavy_neutral_energy_j > 0.0);
        assert!(b.electron_to_heavy_coulomb_energy_j > 0.0);
        assert!(b.entropy_change_j_k > 0.0);
    }
    eprintln!("coupled absolute [x,Te,Th,Qn,Qc,S] errors={errors:?}, steps={accepted:?}");
    // First-order backward Euler with local step doubling has a global error
    // roughly proportional to sqrt(local tolerance), not to that tolerance.
    for column in 0..6 {
        assert!(errors[1][column] < 0.25 * errors[0][column], "{errors:?}");
        assert!(errors[2][column] < 0.25 * errors[1][column], "{errors:?}");
    }
    assert!(errors[2][0] < 1e-7, "{errors:?}");
    assert!(errors[2][1] < 0.5 && errors[2][2] < 0.5, "{errors:?}");
    assert!(errors[2][3] < 0.02 && errors[2][4] < 0.5, "{errors:?}");
    assert!(errors[2][5] < 5e-5, "{errors:?}");
}

#[test]
fn coupled_reverse_exchange_closes_held_work_independently_of_formation_gauge() {
    let mut endpoints = Vec::new();
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        let (state, channel) = parcel(0.05, 6000.0, 12000.0, elemental, charge);
        let collision_sources = sources(&state);
        let result = advance_atomic_ionization_with_collisions(
            &state,
            &channel,
            &collision_sources,
            Power::new(2e7),
            Power::new(-3e6),
            Time::new(2e-7),
            coupled_controls(),
            || false,
        )
        .unwrap();
        assert_budgets(&result);
        let b = result.balance();
        assert_eq!(result.collision_sources(), Some(&collision_sources));
        // Independent DOP853 gives Qc=-826.6745637647523 J and
        // Qn=-2.7568032495634767 J. Its reaction extent is only
        // -2.6684363052531523e-11 mol, below this fixture's requested amount
        // resolution. Resolved recombination is tested separately below.
        assert!(b.electron_to_heavy_neutral_energy_j < 0.0);
        assert!(b.electron_to_heavy_coulomb_energy_j < 0.0);
        assert!(result.state().electron_temperature().value() > 6000.0);
        assert!(result.state().heavy_temperature().value() < 12000.0);
        assert!((b.electron_deposited_energy_j - 4.0).abs() < 1e-14);
        assert!((b.heavy_deposited_energy_j + 0.6).abs() < 1e-14);
        endpoints.push([
            result.state().amounts()[2].value(),
            result.state().electron_temperature().value(),
            result.state().heavy_temperature().value(),
        ]);
    }
    assert!((endpoints[0][0] - endpoints[1][0]).abs() < 1e-12);
    for j in 1..3 {
        assert!((endpoints[0][j] - endpoints[1][j]).abs() < 1e-4);
    }
}

#[test]
fn recombination_requires_the_initially_absent_neutral_momentum_source() {
    let (state, channel) = parcel(1.0, 10000.0, 12000.0, 0.0, 0.0);
    assert_eq!(state.amounts()[0].value(), 0.0);
    let complete = sources(&state);
    let mut missing = complete.clone();
    missing.momentum_sources.clear();
    assert!(matches!(
        advance_atomic_ionization_with_collisions(
            &state,
            &channel,
            &missing,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-7),
            coupled_controls(),
            || false,
        ),
        Err(AtomicIonizationError::Caloric(ElectronHeavyEnergyError::MissingSource(target)))
            if target == id("H")
    ));
    let result = advance_atomic_ionization_with_collisions(
        &state,
        &channel,
        &complete,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(1e-7),
        coupled_controls(),
        || false,
    )
    .unwrap();
    assert!(result.state().amounts()[0].value() > 1e-9);
    assert!(result.reaction_extent().value() < 0.0);
    assert!(result.balance().electron_to_heavy_neutral_energy_j < 0.0);
    assert_budgets(&result);
}

#[test]
fn coupled_sources_refuse_stale_population_receipts_and_cancel_atomically() {
    let (state, channel) = parcel(0.05, 30000.0, 4000.0, 0.0, 0.0);
    let original = state.clone();
    let collision_sources = sources(&state);
    let changed = FrozenIonizedGasStateV1::new(
        state.thermodynamics(),
        state.amounts(),
        Volume::new(1001.0),
        state.heavy_temperature(),
        state.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(matches!(
        advance_atomic_ionization_with_collisions(
            &changed,
            &channel,
            &collision_sources,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-6),
            coupled_controls(),
            || false,
        ),
        Err(AtomicIonizationError::Caloric(
            ElectronHeavyEnergyError::Invalid(_)
        ))
    ));
    for cancel_after in [0, 600] {
        let mut calls = 0;
        assert_eq!(
            advance_atomic_ionization_with_collisions(
                &state,
                &channel,
                &collision_sources,
                Power::new(0.0),
                Power::new(0.0),
                Time::new(1e-6),
                coupled_controls(),
                || {
                    calls += 1;
                    calls > cancel_after
                },
            ),
            Err(AtomicIonizationError::Cancelled)
        );
    }
    assert_eq!(state, original);
}

#[test]
fn coupled_zero_electrons_remain_absent_and_cannot_absorb_electron_work() {
    let (state, channel) = parcel(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let collision_sources = sources(&state);
    let result = advance_atomic_ionization_with_collisions(
        &state,
        &channel,
        &collision_sources,
        Power::new(0.0),
        Power::new(1e6),
        Time::new(1e-6),
        coupled_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.state().amounts(), state.amounts());
    assert_eq!(result.reaction_extent().value(), 0.0);
    assert_eq!(
        result.state().electron_temperature(),
        state.electron_temperature()
    );
    assert!((result.state().heavy_temperature().value() - 4000.0 - 1.0 / (1.5 * R)).abs() < 2e-9);
    let b = result.balance();
    assert_eq!(b.electron_to_heavy_energy_j, 0.0);
    assert_eq!(b.electron_to_heavy_coulomb_energy_j, 0.0);
    assert_eq!(b.electron_to_heavy_neutral_energy_j, 0.0);
    assert_eq!(result.collision_sources(), Some(&collision_sources));
    assert_budgets(&result);
    assert!(matches!(
        advance_atomic_ionization_with_collisions(
            &state,
            &channel,
            &collision_sources,
            Power::new(1.0),
            Power::new(0.0),
            Time::new(1e-6),
            coupled_controls(),
            || false,
        ),
        Err(AtomicIonizationError::Invalid(_))
    ));
}
