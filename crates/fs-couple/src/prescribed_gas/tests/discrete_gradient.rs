//! G1 source isentropes and cycle work; G0 conservative composition and atomicity.
use super::*;
use crate::prescribed_gas::reactive::{self, ReactivePrescribedGasSolve};

fn single_discrete(
    gas: &ClosedGasParcel,
    drive: PrescribedGasDrive,
) -> PrescribedGasDiscreteGradientCandidate {
    step_discrete_gradient_candidate(
        &BTreeMap::from([(1, gas.clone())]),
        &[],
        0.001,
        &BTreeMap::from([(1, energy_solve())]),
        &BTreeMap::from([(1, drive)]),
        solve(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_one_step_variable_cv_compression_and_expansion_follow_source_isentropes() {
    let initial = gas(0.75, 0.002, 400.0, 100000.0, 0.001);
    for volume in [0.0005, 0.0015] {
        // Independent analytic integral: Cv/R=2.25+0.0015*T.
        let (mut lo, mut hi) = (200.0_f64, 1000.0_f64);
        for _ in 0..64 {
            let t = 0.5 * (lo + hi);
            let entropy = 2.25 * (t / 400.0).ln()
                + 0.0015 * (t - 400.0)
                + (volume / initial.volume_m3()).ln();
            if entropy > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let expected_temperature = 0.5 * (lo + hi);
        let actual = single_discrete(&initial, drive(volume, 0.0));
        let report = &actual.step().parcels[&1];
        let end = &actual.parcels()[&1];
        let r = initial.state().state().specific_gas_constant;
        let t = end.state().state().temperature;
        let expected_energy =
            initial.mass_kg() * r * (2.25 * (t - 400.0) + 0.00075 * (t * t - 400.0 * 400.0));
        near(t, expected_temperature, 2e-8);
        near(
            end.state().state().pressure,
            initial.mass_kg() * r * t / volume,
            1e-7,
        );
        near(expected_energy + actual.step().gross_gas_work_j, 0.0, 2e-8);
        near(report.material.entropy_change_j_k, 0.0, 2e-10);
        assert_eq!(report.material.entropy_transfer_j_k, 0.0);
        near(
            report.material.work_pressure_pa * (volume - initial.volume_m3()),
            actual.step().gross_gas_work_j,
            1e-12,
        );
        // The work sample is a thermodynamic mean, not the old midpoint effort.
        let midpoint = single(&initial, drive(volume, 0.0));
        assert!((midpoint.parcels()[&1].state().state().temperature - t).abs() > 0.05);
        assert_eq!(end.mass_kg(), initial.mass_kg());
        assert_eq!(end.point().axes().get("specimen-lot"), Some(&42.0));

        let returned = single_discrete(end, drive(initial.volume_m3(), 0.0));
        near(
            returned.parcels()[&1].state().state().temperature,
            400.0,
            2e-8,
        );
        near(
            returned.step().gross_gas_work_j + actual.step().gross_gas_work_j,
            0.0,
            2e-8,
        );
    }
}

#[test]
fn g1_one_step_otto_legs_match_analytic_efficiency_without_volume_subdivision() {
    let initial = gas(0.75, 0.0, 300.0, 100000.0, 0.001);
    let ratio = 3.0_f64;
    let cv = initial.state().state().specific_gas_constant * 2.25;
    let compressed = single_discrete(&initial, drive(initial.volume_m3() / ratio, 0.0));
    let gas_b = &compressed.parcels()[&1];
    near(
        gas_b.state().state().temperature,
        300.0 * ratio.powf(1.0 / 2.25),
        2e-8,
    );
    let supplied_heat = initial.mass_kg() * cv * (850.0 - gas_b.state().state().temperature);
    let heated = single_discrete(gas_b, drive(gas_b.volume_m3(), supplied_heat));
    let expanded = single_discrete(&heated.parcels()[&1], drive(initial.volume_m3(), 0.0));
    let gas_d = &expanded.parcels()[&1];
    let rejected_heat = initial.mass_kg() * cv * (300.0 - gas_d.state().state().temperature);
    let cooled = single_discrete(gas_d, drive(initial.volume_m3(), rejected_heat));
    let work = compressed.step().gross_gas_work_j + expanded.step().gross_gas_work_j;
    assert!(compressed.step().gross_gas_work_j < 0.0);
    assert!(expanded.step().gross_gas_work_j > 0.0);
    assert_eq!(heated.step().gross_gas_work_j, 0.0);
    assert_eq!(cooled.step().gross_gas_work_j, 0.0);
    near(work / supplied_heat, 1.0 - ratio.powf(-1.0 / 2.25), 2e-10);
    near(work, supplied_heat + rejected_heat, 3e-8);
    near(
        cooled.parcels()[&1].state().state().temperature,
        300.0,
        2e-8,
    );
    near(
        cooled.parcels()[&1].state().state().pressure,
        100000.0,
        1e-5,
    );
}

#[test]
fn g0_open_network_keeps_flow_enthalpy_work_and_post_flow_entropy_distinct() {
    let (mut gases, links, endpoints, drives) = open_network();
    // Incoming gas changes the composition as well as mass and energy.
    gases.insert(1, gas(0.25, 0.0, 400.0, 180000.0, 0.002));
    let inventory = species_inventory(&gases);
    let flow =
        snapshot::step_candidate(&gases, &links, 1e-5, &endpoints, solve().exchange, || false)
            .unwrap();
    let actual = step_discrete_gradient_candidate(
        &gases,
        &links,
        1e-5,
        &endpoints,
        &drives,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(&actual.step().flow, flow.step());
    for (species, mass) in inventory {
        near(species_inventory(actual.parcels())[&species], mass, 1e-14);
    }
    for (&id, end) in actual.parcels() {
        let old = &gases[&id];
        let moved = &flow.parcels()[&id];
        let report = &actual.step().parcels[&id];
        let independent_energy = end.mass_kg() * end.state().specific_internal_energy_j_kg()
            - old.mass_kg() * old.state().specific_internal_energy_j_kg();
        near(
            independent_energy + report.material.work_by_gas_j
                - drives[&id].heat_j
                - report.net_enthalpy_into_gas_j,
            0.0,
            2e-8,
        );
        let ds = moved
            .model()
            .thermodynamics()
            .specific_entropy_change(
                fs_qty::Temperature::new(moved.state().state().temperature),
                fs_qty::Pressure::new(moved.state().state().pressure),
                fs_qty::Temperature::new(end.state().state().temperature),
                fs_qty::Pressure::new(end.state().state().pressure),
            )
            .unwrap()
            .value()
            * moved.mass_kg();
        let theta = report.thermodynamic_mean.state().temperature;
        near(ds, drives[&id].heat_j / theta, 2e-10);
        near(ds, report.material.entropy_change_j_k, 2e-10);
        near(
            report.material.work_pressure_pa * (end.volume_m3() - old.volume_m3()),
            report.material.work_by_gas_j,
            1e-12,
        );
    }
    near(actual.step().total_mass_residual_kg, 0.0, 1e-14);
    near(actual.step().total_energy_residual_j, 0.0, 2e-8);
    near(
        actual.step().net_drive_work_j + actual.step().exterior_pressure_work_j,
        actual.step().gross_gas_work_j,
        1e-12,
    );
    let mut vacuum = drives.clone();
    for drive in vacuum.values_mut() {
        drive.exterior_pressure_pa = 0.0;
    }
    let alternate = step_discrete_gradient_candidate(
        &gases,
        &links,
        1e-5,
        &endpoints,
        &vacuum,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(actual.parcels(), alternate.parcels());
    near(
        alternate.step().net_drive_work_j,
        actual.step().gross_gas_work_j,
        1e-12,
    );

    let checked = reactive::step_discrete_gradient_candidate(
        &gases,
        &links,
        1e-5,
        &endpoints,
        &drives,
        &[],
        ReactivePrescribedGasSolve {
            network: solve(),
            species_tolerance_mol: 1e-12,
        },
        || false,
    )
    .unwrap();
    assert_eq!(checked.parcels(), actual.parcels());
    assert_eq!(checked.step(), actual.step());
    assert!(checked.reactions().is_empty());
    assert!(checked.maximum_species_residual_mol() < 1e-12);
}

#[test]
fn g0_discrete_updates_retain_noop_and_publish_nothing_on_late_refusal_or_cancel() {
    let initial = gas(0.75, 0.002, 400.0, 100000.0, 0.001);
    let idle = single_discrete(&initial, drive(initial.volume_m3(), 0.0));
    assert_eq!(idle.parcels()[&1], initial);
    assert_eq!(idle.step().gross_gas_work_j, 0.0);
    let (gases, links, endpoints, drives) = open_network();
    let saved = gases.clone();
    let mut bad = drives.clone();
    bad.get_mut(&3).unwrap().heat_j = 1e12;
    assert!(matches!(
        step_discrete_gradient_candidate(&gases, &links, 1e-5, &endpoints, &bad, solve(), || false,),
        Err(Error::Parcel { id: 3, .. })
    ));
    bad = drives.clone();
    bad.get_mut(&3).unwrap().final_volume_m3 *= 0.01;
    assert!(matches!(
        step_discrete_gradient_candidate(&gases, &links, 1e-5, &endpoints, &bad, solve(), || false,),
        Err(Error::Parcel { id: 3, .. })
    ));
    let mut total = 0;
    step_discrete_gradient_candidate(&gases, &links, 1e-5, &endpoints, &drives, solve(), || {
        total += 1;
        false
    })
    .unwrap();
    for stop in [1, total / 2, total] {
        let mut count = 0;
        assert_eq!(
            step_discrete_gradient_candidate(
                &gases,
                &links,
                1e-5,
                &endpoints,
                &drives,
                solve(),
                || {
                    count += 1;
                    count == stop
                },
            ),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(gases, saved);
}

#[test]
fn g1_reacting_compression_uses_post_reaction_capacity_and_no_second_heat_supply() {
    use crate::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
    use fs_thermochem::equilibrium::{
        DeclaredStoichiometricMassActionV1, IdealGasReactionEquilibriumV1,
        MassActionKineticOrderConventionV1, ReactionProgressRateQuantityV1,
        ReactionProgressRateScaleV1, ThermodynamicReverseRateClosureV1,
    };
    use fs_thermochem::{
        ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
        UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, verify_conservation,
    };

    // Manufactured A -> 7 B conserves the source masses 28 g/mol -> 7*4 g/mol.
    // Both original constant-Cv species retain zero formation-energy references.
    let initial = gas(0.75, 0.0, 400.0, 100000.0, 0.001);
    let species = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let element = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![7, 1]],
    )
    .unwrap();
    let id = ReactionId::new("synthetic-fragmentation").unwrap();
    let matrix =
        StoichiometricMatrix::new(species.clone(), vec![id.clone()], vec![vec![-1], vec![7]])
            .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&element, &matrix, &charge).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        matrix,
        certificate,
        id,
        initial.model().thermodynamics().components().to_vec(),
    )
    .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(1000.0)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    let law = GasReaction::new(initial.model().clone(), kinetics).unwrap();
    let gases = BTreeMap::from([(1, initial.clone())]);
    let endpoints = BTreeMap::from([(1, energy_solve())]);
    let drives = BTreeMap::from([(1, drive(0.0009, 0.0))]);
    let reactions = [reactive::DrivenGasReaction {
        chamber: 1,
        reaction: &law,
        controls: GasReactionSolve {
            gas: energy_solve(),
            maximum_consumed_fraction: 0.1,
            species_tolerance_mol: 1e-12,
            entropy_tolerance_j_k: 1e-10,
        },
    }];
    let controls = ReactivePrescribedGasSolve {
        network: solve(),
        species_tolerance_mol: 1e-12,
    };
    let mut polls = 0;
    let actual = reactive::step_discrete_gradient_candidate(
        &gases,
        &[],
        5e-5,
        &endpoints,
        &drives,
        &reactions,
        controls,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let initial_moles = initial.mass_kg() / (0.75 * 0.028 + 0.25 * 0.004);
    let original_moles = [0.75 * initial_moles, 0.25 * initial_moles];
    let extent = actual.reactions()[&1].extent_mol;
    assert!(extent > 0.0);
    let amounts = [original_moles[0] - extent, original_moles[1] + 7.0 * extent];
    let cv = R * (2.5 * amounts[0] + 1.5 * amounts[1]);
    let initial_energy = R * 400.0 * (2.5 * original_moles[0] + 1.5 * original_moles[1]);
    let reacted_temperature = initial_energy / cv;
    let expected_temperature = reacted_temperature
        * (initial.volume_m3() / drives[&1].final_volume_m3)
            .powf(R * (amounts[0] + amounts[1]) / cv);
    let end = &actual.parcels()[&1];
    near(end.state().state().temperature, expected_temperature, 2e-8);
    let final_moles = end.mass_kg()
        / end
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
    for (&fraction, expected) in end
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()
        .iter()
        .zip(amounts)
    {
        near(final_moles * fraction, expected, 1e-12);
    }
    near(
        actual.step().gross_gas_work_j,
        -cv * (expected_temperature - reacted_temperature),
        2e-8,
    );
    near(
        actual.step().parcels[&1].material.entropy_change_j_k,
        0.0,
        2e-10,
    );
    assert_eq!(actual.step().supplied_heat_j, 0.0);
    assert!(actual.reactions()[&1].entropy_change_j_k > 0.0);
    assert!(actual.maximum_species_residual_mol() < 1e-12);
    near(actual.step().total_energy_residual_j, 0.0, 2e-8);

    let mut count = 0;
    assert_eq!(
        reactive::step_discrete_gradient_candidate(
            &gases,
            &[],
            5e-5,
            &endpoints,
            &drives,
            &reactions,
            controls,
            || {
                count += 1;
                count == polls
            },
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(gases[&1], initial);
}
