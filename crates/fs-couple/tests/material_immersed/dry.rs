//! Original source density changes dry inertia and reduced gravity together.
use super::*;
use fs_couple::material_immersed::dry::ResolvedImmersedBody;
use fs_lbm::d3q19::DryBodyCouplingConfig3;

fn resting() -> SphereStateSi {
    SphereStateSi {
        center_m: [0.04; 3],
        orientation: So3::identity(),
        velocity_m_s: [0.0; 3],
        angular_momentum_body_kg_m2_s: [0.0; 3],
    }
}
fn dry_snapshot(model: &ResolvedImmersedBody<'_, ResolvedLiquidState>) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                result.push(model.grid().populations(x, y, z));
            }
        }
    }
    result
}

#[test]
fn source_light_neutral_heavy_solids_keep_dry_mass_and_apply_reduced_gravity_once() {
    let source = fluid(1000.0, 0.11);
    for density in [500.0, 1000.0, 2000.0] {
        let solid = body(density);
        let mut model = ResolvedImmersedBody::sphere(
            &source,
            domain(),
            &solid,
            RADIUS,
            resting(),
            [0.0, 0.0, -1e-4],
            &CancelGate::new(),
        )
        .unwrap();
        assert!(core::ptr::eq(model.material(), &solid));
        assert!(core::ptr::eq(model.scaling().material(), &source));
        let dry_mass = solid.inertia().mass();
        close(
            model.mass_kg(model.body_lattice().mass).unwrap(),
            dry_mass,
            1e-16,
        );
        let displaced = 1000.0 * VolumeMoments::sphere(RADIUS).unwrap().volume_m3();
        close(model.displaced_mass_kg(), displaced, 1e-16);
        close(
            model.reduced_gravity_force_n()[2],
            -(dry_mass - displaced) * 1e-4,
            1e-20,
        );
        let report = model
            .step_forces_si(
                [0.0; 6],
                [0.011, -0.013, 0.009],
                DryBodyCouplingConfig3::default(),
                &CancelGate::new(),
            )
            .unwrap();
        let state = model.state_si().unwrap();
        if density == 1000.0 {
            close(state.velocity_m_s[2], 0.0, 1e-13);
            assert_eq!(model.reduced_gravity_force_n(), [0.0; 3]);
        } else {
            assert_eq!(state.velocity_m_s[2].is_sign_positive(), density < 1000.0);
            assert!(model.gravity_potential_change_j(&report).unwrap() < 0.0);
        }
        for x in report.combined_momentum_residual {
            close(x, 0.0, 2e-10);
        }
        close(
            model.mass_kg(model.body_lattice().mass).unwrap(),
            dry_mass,
            1e-16,
        );
        assert_eq!(model.accepted_steps(), 1);
    }
}

#[test]
fn neutral_body_keeps_existing_motion_and_real_reciprocal_drag() {
    let source = fluid(1000.0, 0.11);
    let solid = body(1000.0);
    let mut initial = resting();
    initial.velocity_m_s[0] = 1e-4;
    let mut model = ResolvedImmersedBody::sphere(
        &source,
        domain(),
        &solid,
        RADIUS,
        initial,
        [0.0, 0.0, -1e-4],
        &CancelGate::new(),
    )
    .unwrap();
    let mass = model.body_lattice().mass;
    let inertia = model.body_lattice().inertia;
    let report = model
        .step_forces_si(
            [0.0; 6],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    let after = model.state_si().unwrap();
    assert!(after.velocity_m_s[0] > 0.0);
    assert!(after.velocity_m_s[0] < initial.velocity_m_s[0]);
    assert_eq!(model.body_lattice().mass, mass);
    assert_eq!(model.body_lattice().inertia, inertia);
    assert!(report.fluid[1].momentum[0] > report.fluid[0].momentum[0]);
    assert!(report.interior_correction_impulse[0] > 0.0);
    assert!(report.interior_slip_rms[1] > 0.0);
    close(report.combined_momentum_residual[0], 0.0, 2e-10);
}

#[test]
fn late_implicit_refusal_preserves_source_body_populations_and_clock() {
    let source = fluid(1000.0, 0.11);
    let solid = body(500.0);
    let mut model = ResolvedImmersedBody::sphere(
        &source,
        domain(),
        &solid,
        RADIUS,
        resting(),
        [0.0, 0.0, -1e-4],
        &CancelGate::new(),
    )
    .unwrap();
    let before = *model.body_lattice();
    let fields = dry_snapshot(&model);
    let controls = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..DryBodyCouplingConfig3::default()
    };
    let error = model
        .step_forces_si([0.0; 6], [0.0; 3], controls, &CancelGate::new())
        .unwrap_err();
    assert!(matches!(error, MaterialImmersedError::Native(
        PartialSaturationError3::Invalid { reason }) if reason.contains("budget")));
    assert_eq!(*model.body_lattice(), before);
    assert_eq!(dry_snapshot(&model), fields);
    assert_eq!(model.accepted_steps(), 0);
    assert!(core::ptr::eq(model.scaling().material(), &source));
    let accepted = model
        .step_forces_si(
            [0.0; 6],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    assert!(accepted.candidate_evaluations > 1);
    assert_eq!(model.accepted_steps(), 1);
}
