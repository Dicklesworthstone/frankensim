use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

#[path = "tests/fixture.rs"]
mod fixture;
use fixture::*;

const T: f64 = 298.15;
const VM: f64 = 18e-6;
const VOLUME: f64 = 1e-3;

fn basis() -> SolventActivityBasis {
    SolventActivityBasis::new(VM, AmountConcentration::new(1000.0)).unwrap()
}
fn options() -> PrecipitationOptions {
    let mut options = PrecipitationOptions::default();
    options.aqueous.absolute_charge_tolerance = Amount::new(1e-18);
    options.aqueous.relative_charge_tolerance = fs_qty::Dimensionless::new(1e-12);
    options.aqueous.ph_tolerance = fs_qty::Dimensionless::new(1e-13);
    options.aqueous.max_iterations = 256;
    options
}
fn solution() -> ElectrolyteSolubilitySolution {
    ElectrolyteSolubilitySolution::from_salt(
        salt(T, 4e-6),
        basis(),
        Volume::new(VOLUME),
        Amount::new(0.11),
        vec![(FAMILY.into(), Amount::new(0.24))],
        vec![],
        options(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "actual {actual:.17e}, expected {expected:.17e}",
    );
}
fn phase(solution: &ElectrolyteSolubilitySolution) -> &SaltSpeciation {
    let ElectrolytePhaseEquilibrium::Salt(phase) = solution.equilibrium() else {
        panic!("expected original salt owner")
    };
    phase
}
fn amount(solution: &ElectrolyteSolubilitySolution, species: &str) -> f64 {
    let index = solution
        .species()
        .iter()
        .position(|id| id.as_str() == species)
        .unwrap();
    solution.solution().solute_moles()[index]
}
fn snapshot(solution: &ElectrolyteSolubilitySolution) -> Vec<u64> {
    let mut bits = vec![
        solution.solution().temperature_k().to_bits(),
        solution.solution().solvent_volume_m3().to_bits(),
        solution.component_total().value().to_bits(),
        solution.equilibrium().solid_amount().value().to_bits(),
        solution.equilibrium().aqueous().ph.value().to_bits(),
        solution
            .equilibrium()
            .aqueous()
            .amount_at_unit_activity
            .value()
            .to_bits(),
    ];
    bits.extend(
        solution
            .solution()
            .solute_moles()
            .iter()
            .map(|n| n.to_bits()),
    );
    bits.extend(
        solution
            .family_totals()
            .iter()
            .map(|(_, n)| n.value().to_bits()),
    );
    bits.extend(
        solution
            .spectators()
            .iter()
            .map(|ion| ion.amount.value().to_bits()),
    );
    bits
}

#[test]
fn original_non_unit_salt_and_ph_supply_actual_dissolved_osmotic_particles_once() {
    let solution = solution();
    let phase = phase(&solution);
    // Manufactured independent equilibrium: h=Ka=1e-7, [HCN]=[CN]=0.02
    // and [Ca]=0.01 on n_ref=1. Charge closes at2*Ca-CN=0 and
    // Kdiss=Ca*CN^2=4e-6. The remaining0.1 mol Ca(CN)2 is not dissolved.
    assert_eq!(
        solution.equilibrium().status(),
        PrecipitationStatus::SolidPresent
    );
    near(phase.dissolved_cation.value(), 0.01, 2e-8, 0.0);
    near(phase.dissolved_family.value(), 0.04, 2e-8, 0.0);
    near(phase.solid_amount.value(), 0.1, 2e-9, 0.0);
    near(phase.aqueous.ph.value(), 7.0, 0.0, 2e-8);
    near(amount(&solution, H), 1e-7, 3e-8, 0.0);
    near(amount(&solution, OH), 1e-7, 3e-8, 0.0);
    near(amount(&solution, HC), 0.02, 3e-8, 0.0);
    near(amount(&solution, CN), 0.02, 3e-8, 0.0);
    near(amount(&solution, CA), 0.01, 3e-8, 0.0);
    assert_eq!(solution.species(), [id(H), id(OH), id(HC), id(CN), id(CA)]);
    assert!(!solution.species().iter().any(|id| id.as_str() == SALT));
    assert!(!solution.species().iter().any(|id| id.as_str() == WATER));
    assert_eq!(solution.component_total().value(), 0.11);
    assert_eq!(solution.family_totals()[0].1.value(), 0.24);
    assert!(solution.spectators().is_empty());
    assert_eq!(phase.aqueous.spectators.len(), 1);
    assert_eq!(solution.solvent(), &id(WATER));
    assert_eq!(solution.solution().temperature_k(), T);
    assert_eq!(
        solution
            .equilibrium()
            .standard()
            .reference_pressure()
            .value(),
        100_000.0
    );
    near(
        solution.solution().osmotic_pressure_pa().unwrap(),
        R * T * (0.05 + 2e-7) / VOLUME,
        3e-8,
        0.0,
    );
    near(
        solution.solution().water_activity().unwrap(),
        (-VM * (0.05 + 2e-7) / VOLUME).exp(),
        1e-10,
        0.0,
    );
    near(
        phase.dissolved_cation.value() + phase.solid_amount.value(),
        0.11,
        1e-14,
        0.0,
    );
    near(
        phase.dissolved_family.value() + 2.0 * phase.solid_amount.value(),
        0.24,
        1e-14,
        0.0,
    );
}

#[test]
fn solvent_exchange_and_updated_analytical_stocks_preserve_sources_and_reversible_solids() {
    let initial = solution();
    let before = snapshot(&initial);
    let concentrated = initial.with_solvent_volume_m3(5e-4).unwrap();
    assert_eq!(
        concentrated
            .equilibrium()
            .aqueous()
            .amount_at_unit_activity
            .value(),
        0.5
    );
    assert!(concentrated.equilibrium().solid_amount().value() > 0.1);
    assert!(concentrated.equilibrium().aqueous().ph.value() < 7.0);
    assert_eq!(concentrated.component_total(), initial.component_total());
    assert_eq!(concentrated.family_totals(), initial.family_totals());
    assert!(std::ptr::eq(concentrated.species(), initial.species()));
    assert!(std::ptr::eq(
        concentrated.family_totals(),
        initial.family_totals()
    ));
    for (current, original) in concentrated
        .equilibrium()
        .standard()
        .terms()
        .iter()
        .zip(initial.equilibrium().standard().terms())
    {
        assert_eq!(current.term.species, original.term.species);
        assert_eq!(current.term.coefficient, original.term.coefficient);
        assert_eq!(current.pack_id, original.pack_id);
        assert_eq!(current.answer, original.answer);
    }
    let direct = phase(&initial)
        .model()
        .speciate(
            initial.component_total(),
            initial.family_totals(),
            initial.spectators(),
            Amount::new(0.5),
            options(),
            || false,
        )
        .unwrap();
    assert_eq!(
        phase(&concentrated).dissolved_cation,
        direct.dissolved_cation
    );
    assert_eq!(
        phase(&concentrated).dissolved_family,
        direct.dissolved_family
    );
    assert_eq!(phase(&concentrated).solid_amount, direct.solid_amount);
    assert_eq!(phase(&concentrated).aqueous.ph, direct.aqueous.ph);

    let diluted = concentrated.with_solvent_volume_m3(0.02).unwrap();
    assert_eq!(
        diluted.equilibrium().status(),
        PrecipitationStatus::SolidAbsent
    );
    assert_eq!(phase(&diluted).dissolved_cation.value(), 0.11);
    assert_eq!(phase(&diluted).dissolved_family.value(), 0.24);
    assert_eq!(phase(&diluted).solid_amount.value(), 0.0);
    let returned = diluted.with_solvent_volume_m3(VOLUME).unwrap();
    assert_eq!(snapshot(&returned), before);

    // Adding0.1 mol of the same formula stock adds0.1 mol Ca and0.2 mol
    // analytical cyanide. At this saturated equilibrium it all stays solid.
    let updated = initial
        .with_analytical_totals(Amount::new(0.21), &[Amount::new(0.44)], &[], || false)
        .unwrap();
    near(phase(&updated).solid_amount.value(), 0.2, 2e-9, 0.0);
    near(amount(&updated, CA), 0.01, 3e-8, 0.0);
    assert_eq!(updated.component_total().value(), 0.21);
    assert_eq!(updated.family_totals()[0].1.value(), 0.44);
    let depleted = updated.with_solvent_volume_m3(0.03).unwrap();
    assert_eq!(
        depleted.equilibrium().status(),
        PrecipitationStatus::SolidAbsent
    );
    assert_eq!(phase(&depleted).dissolved_cation.value(), 0.21);
    assert_eq!(phase(&depleted).dissolved_family.value(), 0.44);
    assert_eq!(snapshot(&initial), before);
}

#[test]
fn original_hydrolysis_result_counts_its_solved_cation_and_fixed_counterion_once() {
    let model = hydrolysis(T, 1e-9);
    let chloride = SpectatorIon {
        species: id("Cl[-](aq)"),
        charge: -1,
        amount: Amount::new(2e-5),
    };
    let solution = ElectrolyteSolubilitySolution::from_hydrolysis(
        model,
        basis(),
        Volume::new(VOLUME),
        Amount::new(1.1e-4),
        vec![],
        vec![chloride.clone()],
        options(),
        || false,
    )
    .unwrap();
    let ElectrolytePhaseEquilibrium::Hydrolysis(result) = solution.equilibrium() else {
        panic!("expected original hydrolysis model")
    };
    // Kprec=h^2/aCa=1e-9 at h=1e-7,aCa=1e-5; Cl=2e-5 closes charge.
    near(result.dissolved_amount.value(), 1e-5, 2e-8, 0.0);
    near(result.solid_amount.value(), 1e-4, 2e-9, 0.0);
    near(result.aqueous.ph.value(), 7.0, 0.0, 2e-8);
    assert_eq!(solution.spectators(), &[chloride]);
    assert_eq!(solution.species(), [id(H), id(OH), id(CA), id("Cl[-](aq)")]);
    near(
        solution.solution().osmotic_pressure_pa().unwrap(),
        R * T * 3.02e-5 / VOLUME,
        3e-8,
        0.0,
    );
    assert_eq!(solution.solution().solvent_volume_m3(), VOLUME);
    assert!(
        solution
            .equilibrium()
            .reaction()
            .terms()
            .iter()
            .any(|term| { term.species.as_str() == WATER && term.coefficient < 0 })
    );
}

#[test]
fn activity_basis_layout_duplicates_and_cancelled_trials_leave_accepted_state_unchanged() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(SolventActivityBasis::new(value, AmountConcentration::new(1000.0)).is_err());
        assert!(SolventActivityBasis::new(VM, AmountConcentration::new(value)).is_err());
    }
    assert!(
        SolventActivityBasis::new(VM, AmountConcentration::new(1e-300))
            .unwrap()
            .amount_at_unit_activity(Volume::new(1e-300))
            .is_err()
    );
    let initial = solution();
    let before = snapshot(&initial);
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(initial.with_solvent_volume_m3(volume).is_err());
        assert_eq!(snapshot(&initial), before);
    }
    assert!(
        initial
            .with_analytical_totals(Amount::new(0.1), &[], &[], || false)
            .is_err()
    );
    assert!(
        initial
            .with_analytical_totals(Amount::new(-0.1), &[Amount::new(0.24)], &[], || false)
            .is_err()
    );
    for duplicate in [CA, CN, H, OH, WATER, SALT] {
        let invalid = ElectrolyteSolubilitySolution::from_salt(
            salt(T, 4e-6),
            basis(),
            Volume::new(VOLUME),
            Amount::new(0.11),
            vec![(FAMILY.into(), Amount::new(0.24))],
            vec![SpectatorIon {
                species: id(duplicate),
                charge: 1,
                amount: Amount::new(0.01),
            }],
            options(),
            || false,
        );
        assert!(invalid.is_err(), "duplicate {duplicate} must refuse");
    }
    let too_many = (0..60)
        .map(|i| SpectatorIon {
            species: id(&format!("spectator-{i}")),
            charge: if i % 2 == 0 { 1 } else { -1 },
            amount: Amount::new(0.0),
        })
        .collect();
    assert!(
        ElectrolyteSolubilitySolution::from_salt(
            salt(T, 4e-6),
            basis(),
            Volume::new(VOLUME),
            Amount::new(0.11),
            vec![(FAMILY.into(), Amount::new(0.24))],
            too_many,
            options(),
            || false,
        )
        .is_err()
    );

    let mut calls = 0;
    let accepted = initial
        .with_solvent_volume_m3_cancelled(2e-3, || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 2);
    for stop in [1, calls] {
        let mut observed = 0;
        let trial = initial.with_solvent_volume_m3_cancelled(2e-3, || {
            observed += 1;
            observed >= stop
        });
        assert!(matches!(
            trial,
            Err(ElectrolyteSolubilityError::Speciation(
                AcidBaseError::Cancelled
            ))
        ));
        assert_eq!(observed, stop);
        assert_eq!(snapshot(&initial), before);
    }
    let replay = initial.with_solvent_volume_m3(2e-3).unwrap();
    assert_eq!(snapshot(&accepted), snapshot(&replay));
}
