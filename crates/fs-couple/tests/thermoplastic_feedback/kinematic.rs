//! Cyclic backstress retains its own storage while native transport receives heat.
use super::*;
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;

fn fixture(variable: bool) -> Fixture {
    let mut data = properties(500.0, false);
    data[2] = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6);
    let mut kinematic = scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 2e8);
    if variable {
        // Equal endpoint and initial values cannot hide an interior change in
        // free-energy coefficient; this slice has no d(H_k)/dT entropy term.
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
fn thermal_source_refuses_interior_kinematic_modulus_changes_before_creating_history() {
    let fixture = fixture(true);
    let mut query = fixture.query();
    query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    assert!(matches!(
        CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE),
        Err(FeedbackError::Invalid(_))
    ));
}
