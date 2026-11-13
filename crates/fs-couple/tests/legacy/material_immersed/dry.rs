//! Original source density changes dry inertia and reduced gravity together.
use super::*;
use fs_couple::material_immersed::dry::ResolvedImmersedBody;
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, DryBodyCouplingConfig3,
};

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

#[test]
fn source_scaled_endpoint_history_replays_native_dry_motion_and_commits_one_si_clock() {
    let density = 1000.0;
    let source = fluid(density, 0.11);
    let solid = body(500.0);
    let gate = CancelGate::new();
    let mut initial = resting();
    initial.velocity_m_s[0] = 1e-4;
    let mut model =
        ResolvedImmersedBody::sphere(&source, domain(), &solid, RADIUS, initial, [0.0; 3], &gate)
            .unwrap();
    let uniform = [1e-4, -2e-4, 1.5e-4];
    let field: Vec<[f64; 3]> = (0..512)
        .map(|index| {
            [
                1e-4 * (1 + index % 8) as f64,
                -2e-4 * (1 + (index / 8) % 8) as f64,
                3e-4 * (1 + index / 64) as f64,
            ]
        })
        .collect();
    model.rebind_uniform_fluid_force_si(uniform, &gate).unwrap();
    model.rebind_fluid_force_field_si(&field, &gate).unwrap();
    // Reconstruct the ordinary native grid from explicit SI rho, dx and dt,
    // then use the ORIGINAL source-derived dry mass/tensor and geometry.
    let factor = (DT / DX) * (DT / density);
    let mut native = BoundaryGrid3::new(
        8,
        8,
        8,
        model.scaling().flow().tau,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    native
        .rebind_uniform_force(uniform.map(|f| f * factor), &gate)
        .unwrap();
    native.rebind_force_field(|[x, y, z]| field[(z * 8 + y) * 8 + x].map(|f| f * factor));
    let mut native_body = *model.body_lattice();
    let mut sampler = BodyFrameSampler3::new(&native, 4).unwrap();
    let shape = BodyFrameEllipsoid3::sphere(RADIUS / DX).unwrap();
    let before = *model.body_lattice();
    let populations = dry_snapshot(&model);
    let forces: Vec<_> = (0..512)
        .map(|i| model.grid().force_density(i % 8, (i / 8) % 8, i / 64))
        .collect();
    for (scale, controls, cancel) in [
        (f64::NAN, DryBodyCouplingConfig3::default(), false),
        (
            2.0,
            DryBodyCouplingConfig3 {
                maximum_evaluations: 1,
                ..DryBodyCouplingConfig3::default()
            },
            false,
        ),
        (2.0, DryBodyCouplingConfig3::default(), true),
    ] {
        let refused = CancelGate::new();
        if cancel {
            refused.request();
        }
        assert!(
            model
                .step_forces_si_with_endpoint_fluid_force_scale(
                    [0.0; 6], scale, [0.0; 3], controls, &refused,
                )
                .is_err()
        );
        assert_eq!(*model.body_lattice(), before);
        assert_eq!(dry_snapshot(&model), populations);
        assert_eq!(model.accepted_steps(), 0);
        for (i, force) in forces.iter().enumerate() {
            assert_eq!(
                model.grid().force_density(i % 8, (i / 8) % 8, i / 64),
                *force
            );
        }
    }
    let mut previous_scale = 1.0;
    for (step, scale) in [2.0, 0.0, -1.0, 0.5].into_iter().enumerate() {
        let report = model
            .step_forces_si_with_endpoint_fluid_force_scale(
                [0.0; 6],
                scale,
                [0.0; 3],
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let original = native
            .step_with_dry_partial_saturation_body_and_endpoint_force_scale(
                &mut sampler,
                &shape,
                &mut native_body,
                [0.0; 6],
                scale,
                [0.0; 3],
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        assert_eq!(report, original);
        assert_eq!(*model.body_lattice(), native_body);
        assert_eq!(model.accepted_steps(), (step + 1) as u64);
        assert!(core::ptr::eq(model.material(), &solid));
        assert!(core::ptr::eq(model.scaling().material(), &source));
        for i in 0..512 {
            let cell = [i % 8, (i / 8) % 8, i / 64];
            assert_eq!(
                model.grid().populations(cell[0], cell[1], cell[2]),
                native.populations(cell[0], cell[1], cell[2])
            );
            assert_eq!(
                model.grid().force_density(cell[0], cell[1], cell[2]),
                native.force_density(cell[0], cell[1], cell[2])
            );
        }
        for axis in 0..3 {
            let physical_impulse = field
                .iter()
                .map(|f| {
                    (uniform[axis] + 0.5 * (previous_scale + scale) * f[axis]) * DX.powi(3) * DT
                })
                .sum::<f64>();
            close(
                model
                    .impulse_ns(report.native.coupling.forcing.impulse[axis])
                    .unwrap(),
                physical_impulse,
                2e-18,
            );
        }
        close(
            model.mass_kg(model.body_lattice().mass).unwrap(),
            solid.inertia().mass(),
            1e-16,
        );
        previous_scale = scale;
    }
    assert_ne!(model.state_si().unwrap().center_m, initial.center_m);
}
