//! Source-exact separable nonlinear thermodynamics, using manufactured data.
#[path = "nonlinear_magnetic_thermal/common.rs"]
mod common;
use common::*;
use fs_material::field_storage::nonlinear::thermal::ThermalReversibleMagneticMaterial;

#[test]
fn nonlinear_source_scale_preserves_energy_entropy_inverse_and_all_source_roles() {
    for slope in [-0.001, 0.0, 0.001] {
        let card = card(slope);
        let original = reference(&card);
        let law = law(&card);
        let reference = law.at_temperature(300.0).unwrap();
        assert_eq!(
            reference.magnetic().evaluate(2.0).unwrap(),
            original.evaluate(2.0).unwrap()
        );
        for t in [290.0, 330.0, 360.0] {
            let at = law.at_temperature(t).unwrap();
            let scale = 1.0 + slope * (t - 300.0);
            let volume = 0.3;
            for b in [0.0, 0.5, 1.0, 2.5, 4.0] {
                let response = at.magnetic().evaluate(b).unwrap();
                let state = at.state(b, volume).unwrap();
                near(response.field_strength_a_m, scale * h(b));
                near(response.energy_density_j_m3, scale * psi(b));
                near(
                    at.magnetic()
                        .flux_density_at_field_strength(response.field_strength_a_m)
                        .unwrap(),
                    b,
                );
                near(state.free_energy_j, scale * psi(b) * volume);
                near(state.entropy_j_k, -slope * psi(b) * volume);
                near(
                    state.internal_energy_j,
                    (1.0 - 300.0 * slope) * psi(b) * volume,
                );
                assert_eq!(state.heat_capacity_j_k, 0.0);
                assert!(response.differential_reluctivity_m_h > 0.0);
            }
            let before = 2.5;
            let after = 0.5;
            near(
                at.magnetic()
                    .bregman_density_j_m3([before, 0.0], [after, 0.0])
                    .unwrap(),
                scale * (h(after) * (after - before) - psi(after) + psi(before)),
            );
            assert_eq!(at.material().properties().len(), 7);
            for (a, b) in law
                .material()
                .properties()
                .iter()
                .zip(at.material().properties())
            {
                assert_eq!(a.requirement(), b.requirement());
                assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
            }
            reference
                .magnetic()
                .admit_temperature_refresh(at.magnetic())
                .unwrap();
            if t != 300.0 && slope != 0.0 {
                assert!(
                    reference
                        .magnetic()
                        .admit_fixed_material_refresh(at.material())
                        .is_err()
                );
            }
        }
        let lo = law.at_temperature(329.99).unwrap().state(2.5, 0.3).unwrap();
        let hi = law.at_temperature(330.01).unwrap().state(2.5, 0.3).unwrap();
        near(
            -(hi.free_energy_j - lo.free_energy_j) / 0.02,
            law.at_temperature(330.0)
                .unwrap()
                .state(2.5, 0.3)
                .unwrap()
                .entropy_j_k,
        );
        assert_eq!(law.minimum_heat_capacity_j_k(4.0, 0.3).unwrap(), 0.0);
        assert!(law.at_temperature(360.01).is_err());
        assert!(reference.magnetic().evaluate(4.01).is_err());
    }
}

#[test]
fn source_admission_refuses_missing_scale_nonaffine_path_and_changed_original_law() {
    let source = card(0.001);
    let reference = reference(&source);
    assert!(
        ThermalReversibleMagneticMaterial::new(
            &reference,
            reference.material(),
            &scale_key(),
            "T",
            [290.0, 360.0]
        )
        .is_err()
    );
    for knots in [
        vec![(200.0, 0.9), (300.0, 1.0), (320.0, 1.2), (500.0, 1.2)],
        vec![(200.0, 0.9), (300.0, 1.0), (320.0, -0.1), (500.0, 1.2)],
    ] {
        let other = card_with_scale(knots);
        assert!(
            ThermalReversibleMagneticMaterial::new(
                &common::reference(&other),
                &bundle(&other, 300.0),
                &scale_key(),
                "T",
                [290.0, 360.0]
            )
            .is_err()
        );
    }
    let a = law(&source).at_temperature(330.0).unwrap();
    let other = card(-0.001);
    let b = law(&other).at_temperature(330.0).unwrap();
    assert!(
        a.magnetic()
            .admit_temperature_refresh(b.magnetic())
            .is_err()
    );
    assert!(
        ThermalReversibleMagneticMaterial::new(
            &reference,
            &bundle(&source, 300.0),
            &scale_key(),
            "T",
            [190.0, 360.0]
        )
        .is_err()
    );
}
