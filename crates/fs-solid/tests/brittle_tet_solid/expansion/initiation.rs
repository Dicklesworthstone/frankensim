//! Source-integrated thermal tensile onset through the original native solver.
//! Manufactured coefficients provide analytic checks, not material validation.

use super::*;
use fs_material::brittle::initiation::{TensileInitiationError, TensileInitiationSettings};
use fs_material::brittle::source::initiation::ResolvedTensileInitiationError;
use fs_solid::brittle::solid::expansion::initiation::{
    ThermalTensileElasticTetSolid, ThermalTensileElasticTetStep, ThermalTensileTetInitiationError,
    ThermalTensileTetOnset, ThermalTensileTetProbe,
};

fn initiation_settings() -> TensileInitiationSettings {
    TensileInitiationSettings {
        maximum_strain: 0.03,
        relative_stress_tolerance: 1e-8,
        fraction_tolerance: 1e-10,
        maximum_bisections: 64,
    }
}

fn elastic<'a>(probe: ThermalTensileTetProbe<'a>) -> ThermalTensileElasticTetStep<'a> {
    match probe {
        ThermalTensileTetProbe::Elastic(step) => step,
        ThermalTensileTetProbe::Onset(_) => panic!("expected an elastic thermal endpoint"),
    }
}

fn onset<'a>(probe: ThermalTensileTetProbe<'a>) -> ThermalTensileTetOnset<'a> {
    match probe {
        ThermalTensileTetProbe::Onset(onset) => onset,
        ThermalTensileTetProbe::Elastic(_) => panic!("expected a thermal tensile crossing"),
    }
}

fn y_restraint() -> Vec<(usize, f64)> {
    // The x and z normal coordinates expand freely; y is restrained. With
    // nu=0 this gives one tensile principal stress E*(-free), without a hint.
    (0..12)
        .filter(|&dof| dof != 3 && dof != 11)
        .map(|dof| (dof, 0.0))
        .collect()
}

#[test]
fn thermal_onset_free_expansion_above_total_strain_threshold_stays_elastic() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = expansion_card("thermal onset free expansion");
        let reference = [resolved_expansion(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let cold = [reference[0].at_temperature(&card, 200.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut solid = virgin;
        for next in [&hot, &cold, &reference] {
            let step = elastic(
                solid
                    .probe(
                        next,
                        &RIGID_SUPPORTS,
                        &[],
                        IncrementSettings::default(),
                        initiation_settings(),
                        controls(),
                        None,
                        cx,
                    )
                    .unwrap(),
            );
            let free = next[0].free_linear_strain();
            for (node, xyz) in nodes.iter().enumerate() {
                for axis in 0..3 {
                    assert!(
                        (step.solid.solution().displacement_m[node][axis] - free * xyz[axis]).abs()
                            < 1e-13
                    );
                }
            }
            let mechanical = next[0]
                .mechanical_strain(&step.solid.solution().element_strain[0])
                .unwrap();
            assert!(mechanical.iter().all(|e| e.abs() < 5e-14));
            assert!(
                step.solid.solution().element_stress[0]
                    .iter()
                    .all(|s| s.abs() < 1e-4)
            );
            assert!(step.solid.stored_energy_j() < 1e-12);
            assert_eq!(step.energy.fracture_dissipation_j, 0.0);
            assert_eq!(step.energy.parameter_fracture_dissipation_j, 0.0);
            check_balances(&step.energy);
            solid = step.solid;
        }
        assert!(hot[0].free_linear_strain() > 10.0 * STRENGTH / YOUNG);
        assert!(cold[0].free_linear_strain().abs() > 5.0 * STRENGTH / YOUNG);
        assert_eq!(solid.materials()[0].expansion(), reference[0].expansion());
    });
}

#[test]
fn restrained_cooling_discovers_a_plane_and_replays_from_conserved_old_thermal_state() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = expansion_card("thermal onset cooling");
        let reference = [resolved_expansion(&card, &geometry, 400.0, 400.0)];
        let warm = [reference[0].at_temperature(&card, 398.0).unwrap()];
        let cold = [reference[0].at_temperature(&card, 360.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let accepted = elastic(
            virgin
                .probe(
                    &warm,
                    &y_restraint(),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let strain0 = 2e-5 * 2.0 - 0.5 * 1e-7 * 2.0_f64.powi(2);
        let strain1 = 2e-5 * 40.0 - 0.5 * 1e-7 * 40.0_f64.powi(2);
        let volume = geometry.element_volumes_m3()[0];
        close(
            accepted.solid.solution().element_stress[0][1],
            YOUNG * strain0,
            2e-10,
        );
        close(
            accepted.solid.stored_energy_j(),
            0.5 * YOUNG * strain0 * strain0 * volume,
            2e-10,
        );
        close(
            accepted.energy.eigenstrain_work_j,
            accepted.solid.stored_energy_j(),
            2e-10,
        );
        check_balances(&accepted.energy);
        let snapshot = accepted.solid.solution().clone();
        let prediction = onset(
            accepted
                .solid
                .probe(
                    &cold,
                    &y_restraint(),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let expected_fraction = (STRENGTH / YOUNG - strain0) / (strain1 - strain0);
        let bracket = prediction.first_fraction_bracket();
        assert!(bracket[0] <= expected_fraction + 1e-9 && bracket[1] >= expected_fraction - 1e-9);
        assert!(bracket[1] - bracket[0] <= initiation_settings().fraction_tolerance);
        let event = prediction.events()[0].as_ref().unwrap();
        assert!(event.native().normal()[0].abs() < 1e-12);
        close(event.native().normal()[1].abs(), 1.0, 2e-12);
        assert!(event.native().normal()[2].abs() < 1e-12);
        assert_eq!(event.source().material(), cold[0].band().material());
        assert_eq!(prediction.initial_solution(), &snapshot);
        close(
            prediction.elastic_predictor().element_stress[0][1],
            YOUNG * strain1,
            2e-10,
        );

        let bands = prediction.resolve_bands(cx).unwrap();
        assert_eq!(bands.initial[0].expansion(), warm[0].expansion());
        assert_eq!(bands.next[0].expansion(), cold[0].expansion());
        let width = 0.02 / 3.0;
        close(bands.initial[0].band().law().band_width_m(), width, 2e-12);
        let handoff = prediction.handoff(&bands.initial, 1e-9, cx).unwrap();
        assert_eq!(
            handoff.solid.solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            handoff.solid.solution().element_strain,
            snapshot.element_strain
        );
        close(
            handoff.solid.solution().element_states[0].strain()[1],
            strain0,
            2e-10,
        );
        close(
            handoff.solid.stored_energy_j(),
            accepted.solid.stored_energy_j(),
            2e-10,
        );
        assert!(handoff.storage_defect_j.abs() < 1e-12);
        assert!(handoff.maximum_stress_defect_pa < 1e-6 * STRENGTH);
        assert_eq!(handoff.solid.fracture_dissipation_j(), 0.0);
        assert_eq!(handoff.solid.solution().element_states[0].damage(), 0.0);
        assert!(handoff.solid.solution().increments.is_empty());
        assert!(prediction.handoff(&warm, 1e-9, cx).is_err());
        assert!(prediction.handoff(&bands.next, 1e-9, cx).is_err());

        let replayed = prediction.replay(&bands, controls(), cx).unwrap();
        let limit = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * width);
        let stress = STRENGTH * (failure - strain1) / (failure - limit);
        close(replayed.solid.solution().element_stress[0][1], stress, 2e-8);
        assert!(replayed.solid.solution().element_stress[0][0].abs() < 1e-4);
        assert!(replayed.solid.solution().element_stress[0][2].abs() < 1e-4);
        let full_opening_work =
            0.5 * STRENGTH * limit + 0.5 * (STRENGTH + stress) * (strain1 - limit);
        let expected_loss = volume * (full_opening_work - 0.5 * stress * strain1);
        close(replayed.energy.fracture_dissipation_j, expected_loss, 2e-8);
        close(
            replayed.energy.eigenstrain_work_j,
            volume * full_opening_work - accepted.solid.stored_energy_j(),
            2e-8,
        );
        assert_eq!(replayed.energy.parameter_work_j, 0.0);
        assert!(replayed.energy.mechanical_work_j.abs() < 1e-9);
        check_balances(&replayed.energy);
        assert_eq!(accepted.solid.solution(), &snapshot);
        with_cx(true, |cancelled| {
            assert!(matches!(
                prediction.resolve_bands(cancelled),
                Err(ThermalTensileTetInitiationError::Thermal(
                    ThermalBrittleTetError::Solid(BrittleTetError::Cancelled)
                ))
            ));
            assert!(prediction.replay(&bands, controls(), cancelled).is_err());
        });
    });
}

fn weakening_card() -> MaterialCard {
    let original = expansion_card("thermal onset coefficient crossing");
    let mut claims = ClaimSet::new();
    for (_, original_claim) in original.claims().claims_ordered() {
        let mut claim = original_claim.clone();
        let variation = match claim.key.name() {
            YOUNG_MODULUS_PROPERTY => Some((YOUNG, 0.5 * YOUNG)),
            BRITTLE_TENSILE_STRENGTH_PROPERTY => Some((STRENGTH, 0.25 * STRENGTH)),
            _ => None,
        };
        if let Some((reference, hot)) = variation {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, reference), (400.0, hot)],
                dims: Pressure::DIMS,
            };
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                400.0,
            );
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

#[test]
fn thermal_coefficient_work_is_separate_and_already_critical_start_is_typed_refusal() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = weakening_card();
        let reference = [resolved_expansion(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        for fraction in [0.1, 0.75] {
            let strain = fraction * STRENGTH / YOUNG;
            let target = affine_targets(&nodes, [0.0, strain, 0.0]);
            let accepted = elastic(
                virgin
                    .probe(
                        &reference,
                        &target,
                        &[],
                        IncrementSettings::default(),
                        initiation_settings(),
                        controls(),
                        None,
                        cx,
                    )
                    .unwrap(),
            );
            let snapshot = accepted.solid.solution().clone();
            let result = accepted.solid.probe(
                &hot,
                &target,
                &[],
                IncrementSettings::default(),
                initiation_settings(),
                controls(),
                None,
                cx,
            );
            if fraction < 0.5 {
                let step = elastic(result.unwrap());
                close(
                    step.energy.parameter_work_j,
                    -0.5 * accepted.solid.stored_energy_j(),
                    2e-10,
                );
                assert_eq!(step.energy.mechanical_work_j, 0.0);
                assert_eq!(step.energy.fracture_dissipation_j, 0.0);
                check_balances(&step.energy);
            } else {
                assert!(matches!(
                    result,
                    Err(ThermalTensileTetInitiationError::Element {
                        element: 0,
                        source: ResolvedTensileInitiationError::Native(
                            TensileInitiationError::StartNotElastic { .. }
                        ),
                    })
                ));
            }
            assert_eq!(accepted.solid.solution(), &snapshot);
        }
    });
}

#[test]
fn fully_restrained_isotropic_cooling_requires_a_physical_principal_direction_hint() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = expansion_card("thermal onset repeated principal tension");
        let reference = [resolved_expansion(&card, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&card, 360.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            virgin.probe(
                &cold,
                &fixed(),
                &[],
                IncrementSettings::default(),
                initiation_settings(),
                controls(),
                None,
                cx
            ),
            Err(ThermalTensileTetInitiationError::Element {
                element: 0,
                source: ResolvedTensileInitiationError::Native(
                    TensileInitiationError::AmbiguousPlane { multiplicity: 3 }
                ),
            })
        ));
        let prediction = onset(
            virgin
                .probe(
                    &cold,
                    &fixed(),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    Some(&[Some([0.0, 1.0, 0.0])]),
                    cx,
                )
                .unwrap(),
        );
        assert_eq!(
            prediction.events()[0].as_ref().unwrap().native().normal(),
            [0.0, 1.0, 0.0]
        );
        assert_eq!(virgin.stored_energy_j(), 0.0);
    });
}
