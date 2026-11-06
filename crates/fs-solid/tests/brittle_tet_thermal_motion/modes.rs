//! G1/G3: thermal damage changes physical vibration without changing matter.

use super::*;
use fs_solid::brittle::solid::expansion::dynamics::modes::{
    FrozenBrittleModesControls, FrozenBrittlePencil,
};

fn mode_controls() -> FrozenBrittleModesControls {
    FrozenBrittleModesControls {
        maximum_dofs: 12,
        maximum_staged_entries: 144,
    }
}

fn axial_dofs() -> Vec<usize> {
    (0..12).filter(|&dof| dof != 3).collect()
}

fn check_axial_mode(pencil: &FrozenBrittlePencil, modulus: f64) {
    assert_eq!(pencil.free_dofs(), &[3]);
    let stiffness = modulus * VOLUME;
    let mass = DENSITY * VOLUME / 10.0;
    near(pencil.stiffness_n_m().get(0, 0), stiffness, 1e-10);
    near(pencil.mass_kg().get(0, 0), mass, 1e-15);
    near(pencil.reference_mass_kg(), DENSITY * VOLUME, 1e-15);
    // The existing eigenproblem owner consumes the actual extracted pencil;
    // the P1 one-free-DOF oracle distinguishes consistent mass from lumping.
    let modes = fs_modal::eigh_gen_dense(
        &pencil.stiffness_n_m().to_dense(),
        &pencil.mass_kg().to_dense(),
        1,
    )
    .unwrap();
    assert_eq!(modes.len(), 1);
    near(modes[0].lambda, stiffness / mass, 1e-9);
    near(modes[0].phi[0].powi(2) * mass, 1.0, 1e-14);
}

#[test]
fn cooling_damage_and_thermal_reclosure_use_mechanical_strain_for_vibration() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("frozen-thermal-modes", 4.5e-5, false);
        let reference = [source(&card, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let warm = [reference[0].at_temperature(&card, 450.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 3, cx).unwrap();
        let fixed = axial_dofs();
        check_axial_mode(
            &initial
                .frozen_vibration_pencil(&fixed, mode_controls(), cx)
                .unwrap(),
            YOUNG,
        );
        let motion = initial
            .advance(
                &cold,
                0.01,
                &axial_supports(),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .motion;
        let state = motion.solid().solution().element_states[0];
        let snapshot = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let fracture = motion.solid().fracture_dissipation_j();
        let kinetic = motion.kinetic_energy_j();
        assert!(state.damage() > 0.9 && state.damage() < 1.0);
        // Cooling creates tensile mechanical strain even while the actual
        // nodal displacement gives a compressive TOTAL axial strain.
        assert!(state.strain()[0] > 0.0);
        assert!(snapshot.element_strain[0][0] < 0.0);
        let pencil = motion
            .frozen_vibration_pencil(&fixed, mode_controls(), cx)
            .unwrap();
        check_axial_mode(&pencil, YOUNG * (1.0 - state.damage()));

        // Reverse the distinction: positive TOTAL expansion can leave the
        // damaged plane compressed mechanically and recover its full normal
        // stiffness. Damage itself remains irreversible.
        let compressed = motion
            .advance(
                &warm,
                0.02,
                &affine_targets(0.00025),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .motion;
        let compressed_state = compressed.solid().solution().element_states[0];
        assert!(compressed_state.strain()[0] < 0.0);
        assert!(compressed.solid().solution().element_strain[0][0] > 0.0);
        near(compressed_state.damage(), state.damage(), 1e-14);
        check_axial_mode(
            &compressed
                .frozen_vibration_pencil(&fixed, mode_controls(), cx)
                .unwrap(),
            YOUNG,
        );

        // At zero MECHANICAL opening, one two-sided stiffness would suppress
        // the real opening/reclosure switch. Positive thermal strain does not
        // make that nondifferentiable point an admitted tensile branch.
        let closure = compressed
            .advance(
                &warm,
                0.02,
                &affine_targets(warm[0].free_linear_strain()),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .motion;
        assert_eq!(
            closure.solid().solution().element_states[0].strain()[0],
            0.0
        );
        assert!(matches!(
            closure.frozen_vibration_pencil(&fixed, mode_controls(), cx),
            Err(BrittleTetError::Material { element: 0, .. })
        ));

        let repeat = motion
            .frozen_vibration_pencil(&fixed, mode_controls(), cx)
            .unwrap();
        assert_eq!(pencil.stiffness_n_m(), repeat.stiffness_n_m());
        assert_eq!(pencil.mass_kg(), repeat.mass_kg());
        assert!(matches!(
            motion.frozen_vibration_pencil(
                &fixed,
                FrozenBrittleModesControls {
                    maximum_dofs: 11,
                    ..mode_controls()
                },
                cx,
            ),
            Err(BrittleTetError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            motion.frozen_vibration_pencil(&fixed, mode_controls(), cx),
            Err(BrittleTetError::Cancelled)
        ));
        assert_eq!(motion.solid().solution(), &snapshot);
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.solid().fracture_dissipation_j(), fracture);
        assert_eq!(motion.kinetic_energy_j(), kinetic);
        assert_eq!(motion.time_s(), 0.01);
        assert_eq!(motion.accepted_steps(), 1);
        assert!(std::ptr::eq(motion.solid().materials(), cold.as_slice()));
    });
}

fn changing_modulus_card() -> MaterialCard {
    let original = card("thermal-modal-stiffness-and-mass", 0.0, true);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == YOUNG_MODULUS_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, YOUNG), (400.0, 0.5 * YOUNG)],
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
fn hot_modal_stiffness_uses_current_coefficients_and_original_reference_mass() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = changing_modulus_card();
        let reference = [source(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1, cx).unwrap();
        let fixed = axial_dofs();
        let cold_pencil = initial
            .frozen_vibration_pencil(&fixed, mode_controls(), cx)
            .unwrap();
        let motion = initial
            .advance(
                &hot,
                0.1,
                &affine_targets(0.00025),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .motion;
        assert_eq!(motion.remaining_steps(), 0);
        assert_eq!(motion.solid().solution().element_states[0].damage(), 0.0);
        near(hot[0].band().law().elastic().youngs, 0.5 * YOUNG, 1e-14);
        near(
            hot[0].band().reference_density_kg_m3().unwrap(),
            0.5 * DENSITY,
            1e-14,
        );
        let hot_pencil = motion
            .frozen_vibration_pencil(&fixed, mode_controls(), cx)
            .unwrap();
        check_axial_mode(&cold_pencil, YOUNG);
        check_axial_mode(&hot_pencil, 0.5 * YOUNG);
        assert_eq!(cold_pencil.mass_kg(), hot_pencil.mass_kg());
        assert_eq!(motion.reference_density_kg_m3(), &[DENSITY]);
        // Extracting a spectrum is not another physical timestep or source
        // query; exhausted motion allowance does not erase the accepted state.
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 0);
        assert!(std::ptr::eq(motion.solid().materials(), hot.as_slice()));
    });
}
