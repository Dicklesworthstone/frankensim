//! The alternative line-friction model drives the original dense capillary owner.

use super::*;
use fs_lbm::material::free_surface::wetting::{
    FreeSurfaceWettingError, MolecularKineticControls, MolecularKineticError,
    MolecularKineticLaw,
};

fn molecular_law(maximum_speed_m_s: f64) -> MolecularKineticLaw {
    MolecularKineticLaw::try_new(MolecularKineticControls {
        // Authored synthetic calibration, not measurements of this fixture.
        jump_length_m: 5e-10,
        equilibrium_jump_frequency_per_s: 1e9,
        minimum_temperature_k: 290.0,
        maximum_temperature_k: 500.0,
        maximum_speed_m_s,
        minimum_angle_radians: 0.01,
        maximum_angle_radians: 3.0,
    })
    .unwrap()
}

#[test]
fn g1_molecular_predictions_keep_static_receipts_and_drive_native_wetting() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let source = static_card(&bulk, 1.2);
    let constant = dynamic_card(&bulk, [1.2; 2], [1.2; 2], 100.0, false);
    let gate = CancelGate::new_clock_free();
    let mut molecular = flow2(&bulk, &source, PlanarWall2::YMin, 0.00002)
        .with_molecular_kinetic(molecular_law(100.0), "T", &gate)
        .unwrap();
    let mut reference = flow2(&bulk, &constant, PlanarWall2::YMin, 0.00002);
    let mass = molecular.flow().ledger_mass_kg_per_m().unwrap();

    molecular.step_with_gate(&gate).unwrap();
    reference.step_with_gate(&gate).unwrap();
    same_state(molecular.flow().lattice(), reference.flow().lattice());
    assert!(molecular.samples().values().any(|sample| sample.signed_speed_m_s != 0.0));
    for sample in molecular.samples().values() {
        assert!(sample.hydrodynamic.is_none());
        assert_eq!(sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS), Some(0.0));
        let prediction = sample.molecular_kinetic.unwrap();
        assert_eq!(prediction.temperature_k, 300.0);
        assert_eq!(prediction.surface_tension_n_m, 0.05);
        let speed = if sample.regime == ContactLineRegime2::Pinned {
            0.0
        } else {
            sample.signed_speed_m_s
        };
        assert_eq!(prediction.signed_speed_m_s, speed);
        let expected_force =
            (2.0 * 1.380_649e-23 * 300.0 / 5e-10_f64.powi(2)) * speed.asinh();
        assert!(
            (prediction.driving_force_per_length_n_m - expected_force).abs()
                < 2e-12 * expected_force.abs().max(1e-20)
        );
        assert_eq!(prediction.applied_angle_radians, sample.applied_angle_radians);
        assert_eq!(
            prediction.dissipation_per_length_w_m,
            prediction.driving_force_per_length_n_m * speed,
        );
        assert!(prediction.dissipation_per_length_w_m >= 0.0);
        for name in [RECEDING_CONTACT_ANGLE_PROPERTY, ADVANCING_CONTACT_ANGLE_PROPERTY] {
            assert_eq!(
                sample.material.property(name).unwrap().answer().receipt.selected,
                molecular.source().property(name).unwrap().answer().receipt.selected,
            );
        }
    }
    molecular.step_with_gate(&gate).unwrap();
    reference.step_with_gate(&gate).unwrap();
    assert_ne!(molecular.flow().lattice().grid.f, reference.flow().lattice().grid.f);
    assert!((molecular.flow().ledger_mass_kg_per_m().unwrap() - mass).abs() < mass * 1e-10);
}

#[test]
fn g3_activated_molecular_wetting_changes_native_capillary_momentum() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let source = static_card(&bulk, 1.2);
    let gate = CancelGate::new_clock_free();
    let fixed_law = molecular_law(100.0);
    let activated_law = MolecularKineticLaw::try_new_arrhenius(
        *fixed_law.controls(), 400.0, 2000.0 * 1.380_649e-23,
    ).unwrap();
    let make = || flow2(&bulk, &source, PlanarWall2::YMin, 0.00002);
    let mut fixed = make().with_molecular_kinetic(fixed_law, "T", &gate).unwrap();
    let mut activated = make().with_molecular_kinetic(activated_law, "T", &gate).unwrap();
    let mass = activated.flow().ledger_mass_kg_per_m().unwrap();
    fixed.step_with_gate(&gate).unwrap();
    activated.step_with_gate(&gate).unwrap();
    same_state(fixed.flow().lattice(), activated.flow().lattice());
    assert!(activated.samples().iter().any(|(position, sample)| {
        sample.applied_angle_radians != fixed.samples()[position].applied_angle_radians
    }));
    let expected_frequency = 1e9 * (2000.0_f64 * (1.0 / 400.0 - 1.0 / 300.0)).exp();
    for (position, sample) in activated.samples() {
        assert_eq!(sample.material, fixed.samples()[position].material);
        let actual = sample.molecular_kinetic.unwrap().equilibrium_jump_frequency_per_s;
        assert!((actual / expected_frequency - 1.0).abs() < 2e-14);
    }
    fixed.step_with_gate(&gate).unwrap();
    activated.step_with_gate(&gate).unwrap();
    assert_ne!(fixed.flow().lattice().grid.f, activated.flow().lattice().grid.f);
    assert!((activated.flow().ledger_mass_kg_per_m().unwrap() - mass).abs() < mass * 1e-10);
}

#[test]
fn g4_molecular_model_rejects_double_rates_and_rolls_back_domain_exit() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let source = static_card(&bulk, 1.2);
    let dynamic = dynamic_card(&bulk, [1.2, 1.0], [1.2, 1.4], 100.0, false);
    let gate = CancelGate::new_clock_free();
    let make = || flow2(&bulk, &source, PlanarWall2::YMin, 0.00002);
    let cox = || CoxVoinovLaw::try_new(cox_controls(0.01), &gate).unwrap();
    assert!(
        make().with_cox_voinov(cox(), &gate).unwrap()
            .with_molecular_kinetic(molecular_law(100.0), "T", &gate).is_err()
    );
    assert!(
        make().with_molecular_kinetic(molecular_law(100.0), "T", &gate).unwrap()
            .with_cox_voinov(cox(), &gate).is_err()
    );
    assert!(
        flow2(&bulk, &dynamic, PlanarWall2::YMin, 0.00002)
            .with_molecular_kinetic(molecular_law(100.0), "T", &gate).is_err()
    );
    assert!(
        make().with_molecular_kinetic(
            molecular_law(100.0), "not_the_liquid_temperature", &gate
        ).is_err()
    );

    let mut limited = make()
        .with_molecular_kinetic(molecular_law(1e-18), "T", &gate)
        .unwrap();
    let before = limited.flow().lattice().clone();
    let samples = limited.samples().clone();
    assert!(matches!(
        limited.step_with_gate(&gate),
        Err(FreeSurfaceWettingError::MolecularKinetic(
            MolecularKineticError::OutsideDomain(_)
        ))
    ));
    same_state(limited.flow().lattice(), &before);
    assert_eq!(limited.samples(), &samples);
}
