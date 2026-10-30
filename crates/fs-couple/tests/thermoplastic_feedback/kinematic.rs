//! Cyclic backstress retains its own storage while native transport receives heat.
use super::*;
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;

fn fixture(variable: bool) -> Fixture {
    let mut data = properties(500.0, false);
    data[2] = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6);
    let mut kinematic = scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 2e8);
    if variable {
        // Equal endpoints cannot hide an interior derivative jump: the
        // admitted affine free energy has no entropy discontinuity.
        curve(
            &mut kinematic,
            "T",
            vec![(280.0, 2e8), (290.0, 2.2e8), (300.0, 2e8), (320.0, 2e8)],
        );
    }
    data.push(kinematic);
    Fixture::from_properties(
        data,
        500.0,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    )
}

#[test]
fn cyclic_mixed_hardening_heats_the_original_caloric_store_without_heating_backstress_storage() {
    context(|gate, cx| {
        let fixture = fixture(false);
        let mut query = fixture.query();
        query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
        let mut point = CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).unwrap();
        let ambient = environment(300.0, 0.0);
        for strain in [shear(0.01), shear(-0.012), shear(0.008)] {
            let old_temperature = point.caloric_state().temperature_k();
            let old_alpha = point.history().alpha;
            let step = point
                .advance(
                    strain,
                    1.0,
                    config(),
                    |input| ambient.advance(cx, input),
                    || gate.is_requested(),
                )
                .unwrap();
            let alpha_increment = step.next.history().alpha - old_alpha;
            let expected_heat = VOLUME * 3e6 * alpha_increment;
            close(step.physical_heat_j, expected_heat, 1e-10);
            close(
                step.next.caloric_state().temperature_k() - old_temperature,
                expected_heat / (point.mass_kg() * 500.0),
                2e-8,
            );
            assert!(step.mechanical.energy.final_state.kinematic_hardening_j_m3 > 0.0);
            assert!(step.mechanical.energy.numerical_dissipation_j_m3 > 0.0);
            assert_eq!(step.transport.external_heat_j, 0.0);
            assert!(step.combined_residual_j.abs() <= config().energy_tolerance_j);
            let source = step.next.resolved().kinematic_hardening_source().unwrap();
            assert_eq!(source.value_si(), 2e8);
            assert_eq!(
                source.answer().receipt.query_point,
                vec![("T".to_owned(), step.next.caloric_state().temperature_k())]
            );
            assert_eq!(step.next.resolved().material().properties().len(), 7);
            point = step.next;
        }
        assert_eq!(point.time_s(), 3.0);
        assert!(point.history().alpha > 0.02);
    });
}

#[test]
fn thermal_source_refuses_interior_kinematic_modulus_kinks_before_creating_history() {
    let fixture = fixture(true);
    let mut query = fixture.query();
    query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    assert!(matches!(
        CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE),
        Err(FeedbackError::Invalid(_))
    ));
}

fn affine_kinematic(slope: f64) -> Fixture {
    let mut data = properties(500.0, false);
    data[2] = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6);
    let mut kinematic = scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 2e8);
    curve(
        &mut kinematic,
        "T",
        [280.0, 290.0, 300.0, 320.0]
            .map(|t| (t, 2e8 + slope * (t - 300.0)))
            .to_vec(),
    );
    data.push(kinematic);
    Fixture::from_properties(
        data,
        500.0,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    )
}

#[test]
fn affine_backstress_storage_and_entropy_cancel_at_frozen_history_even_at_zero_modulus() {
    use fs_material::SmallStrainLaw;
    use fs_material::plastic::J2State;
    use fs_material::plastic::energy::EnergyTolerances;
    use fs_material::plastic::thermal::ThermoplasticSource;

    let slope = -1e7;
    let fixture = affine_kinematic(slope);
    let mut query = fixture.query();
    query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    let source = ThermoplasticSource::try_new(query, 0.0).unwrap();
    let history = J2State {
        plastic_strain: [0.002, -0.001, -0.001, 0.0003, -0.0002, 0.0001],
        alpha: 0.01,
    };
    let p = history.plastic_strain;
    // Tensor-Voigt shear weights belong to the original Prager potential.
    let shape =
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + 2.0 * (p[3] * p[3] + p[4] * p[4] + p[5] * p[5]))
            / 3.0;
    for (old_t, new_t) in [(290.0, 320.0), (320.0, 290.0)] {
        let step = source
            .advance_with_energy(
                500.0 * (old_t - 300.0),
                500.0 * (new_t - 300.0),
                &p,
                &history,
                &p,
                EnergyTolerances {
                    energy_j_m3: 1e-6,
                    yield_pa: 1e-3,
                },
            )
            .unwrap();
        let coefficient = slope * (new_t - old_t) * shape;
        assert_eq!(step.mechanical.state, history);
        close(step.coefficient_storage_change_j_m3, coefficient, 1e-9);
        close(step.stored_change_j_m3(), coefficient, 1e-9);
        close(step.thermoelastic_storage_change_j_m3, -coefficient, 1e-9);
        close(step.reversible_heat_j_m3, 0.0, 1e-9);
        close(step.caloric_heat_j_m3(), 0.0, 1e-9);
        close(step.mechanical_work_j_m3, 0.0, 1e-9);
        let (_, resolved) = source.resolve_enthalpy(500.0 * (new_t - 300.0)).unwrap();
        let claim = resolved.kinematic_hardening_source().unwrap();
        close(claim.value_si(), 2e8 + slope * (new_t - 300.0), 1e-6);
        fixture
            .card
            .claims()
            .verify_receipt(&claim.answer().receipt)
            .unwrap();
    }
    // The zero-H_k endpoint must also admit an active native plastic return
    // and its consistent mechanical tangent, rather than only an elastic rest.
    let strain = [0.01, -0.005, -0.005, 0.0, 0.0, 0.0];
    let yielded = source
        .advance_with_energy(
            10000.0,
            10000.0,
            &p,
            &history,
            &strain,
            EnergyTolerances {
                energy_j_m3: 1e-6,
                yield_pa: 1e-3,
            },
        )
        .unwrap();
    assert!(yielded.mechanical.state.alpha > history.alpha);
    let (_, zero) = source.resolve_enthalpy(10000.0).unwrap();
    assert_eq!(zero.law().kinematic_hardening, 0.0);
    assert!(
        zero.law()
            .tangent(&strain, &history)
            .iter()
            .flatten()
            .all(|x| x.is_finite())
    );
}

#[test]
fn thermal_backstress_curve_cannot_cross_into_negative_stored_energy() {
    let fixture = affine_kinematic(-1.1e7);
    let mut query = fixture.query();
    query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    assert!(CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).is_err());
}
