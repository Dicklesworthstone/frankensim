//! G1/G3: same finite equilibrium, original mass and conservative field Hessian.
use super::*;
use fs_solid::hyper3::gravity::GravityField3;
use fs_solid::hyper3::vibration::{VibrationError, VibrationSettings};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const GRADED: [f64; 4] = [3.0, 3.0, 3.0, 15.0];

fn controls() -> VibrationSettings {
    VibrationSettings {
        maximum_free_dofs: 12,
        relative_residual_tolerance: 1e-9,
    }
}

fn supports() -> Vec<(usize, f64)> {
    (0..12).filter(|&d| d != 11).map(|d| (d, 0.0)).collect()
}

fn plan(lambda: f64) -> fs_modal::WindowPlan {
    let half = 0.4 * lambda.abs();
    fs_modal::WindowPlan {
        windows: vec![(lambda - half, lambda + half)],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 4 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}

fn options() -> fs_modal::SliceOptions {
    fs_modal::SliceOptions {
        max_lanczos: 32,
        ritz_tol: 1e-12,
        ldlt: fs_sparse::LdltOptions {
            max_factor_bytes: Some(2 * 1024 * 1024),
            ..Default::default()
        },
        ..Default::default()
    }
}

// Differentiate V0 W(diag(1,1,s)) directly, with the parent's mu=2, lambda=3.
fn material_force_stiffness(s: f64) -> (f64, f64) {
    (
        (2.0 * s + (3.0 * s.ln() - 2.0) / s) / 6.0,
        (2.0 + (2.0 + 3.0 * (1.0 - s.ln())) / (s * s)) / 6.0,
    )
}

fn field(x: [f64; 3], constant: f64, gradient: f64) -> Result<GravitySample3, HyperTetError> {
    Ok(GravitySample3 {
        potential_j_kg: -constant * x[2] - 0.5 * gradient * x[2] * x[2],
        acceleration_m_s2: [0.0, 0.0, constant + gradient * x[2]],
        acceleration_gradient_s2: [[0.0; 3], [0.0; 3], [0.0, 0.0, gradient]],
    })
}

#[test]
fn same_loaded_equilibrium_retains_gravity_pressure_and_graded_modal_stiffness() {
    let law = material();
    let fixed = supports();
    let mut p = problem(&NODES, &TETS, &law, &fixed);
    p.settings.load_steps = 3;
    let stretch = 1.25;
    let (internal, material_k) = material_force_stiffness(stretch);
    with_cx(|cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, Default::default(), 4, cx).unwrap();
        for graded in [false, true] {
            // Independent simplex moments of rho0*N3 and rho0*N3^2. Both
            // inventories weigh 1 kg; a top-heavy field has different moments.
            let (share, mass) = if graded { (0.325, 0.15) } else { (0.25, 0.1) };
            for gravity_gradient in [0.0, -3.0] {
                for pressure in [false, true] {
                    let (p0, dp) = if pressure { (0.3, -0.4) } else { (0.0, 0.0) };
                    // Gravity force = g0*share + dg*mass*s. Pressure potential
                    // is p0*s/6 + dp*s^2/24. No production force enters this oracle.
                    let constant = (internal + p0 / 6.0 + dp * stretch / 12.0
                        - gravity_gradient * mass * stretch)
                        / share;
                    let gravity = |x| field(x, constant, gravity_gradient);
                    let gravity: &GravityField3<'_> = &gravity;
                    let held_pressure = pressure.then_some((
                        &surface,
                        TetPressure::Affine {
                            origin_m: [0.0; 3],
                            pressure_pa: p0,
                            gradient_pa_m: [0.0, 0.0, dp],
                        },
                    ));
                    let dense = if graded {
                        p.solve_vibration_with_nodal_density_gravity(
                            &GRADED,
                            gravity,
                            held_pressure,
                            controls(),
                            cx,
                        )
                    } else {
                        p.solve_vibration_with_gravity(
                            &[6.0],
                            gravity,
                            held_pressure,
                            controls(),
                            cx,
                        )
                    }
                    .unwrap();
                    let expected_k = material_k + dp / 12.0 - gravity_gradient * mass;
                    let expected_lambda = expected_k / mass;
                    near(
                        dense.equilibrium().displacement_m[3][2],
                        stretch - 1.0,
                        2e-9,
                    );
                    assert_eq!(dense.free_dofs(), &[11]);
                    near(dense.stiffness_n_m()[0], expected_k, 2e-8);
                    near(dense.mass_kg()[0], mass, 2e-14);
                    near(dense.modes()[0].lambda, expected_lambda, 3e-7);
                    assert!(core::ptr::eq(dense.held_gravity_field().unwrap(), gravity));
                    assert_eq!(dense.held_pressure_field().is_some(), pressure);
                    assert_eq!(dense.reference_density_kg_m3(), &[6.0]);
                    assert_eq!(
                        dense.reference_nodal_density_kg_m3(),
                        graded.then_some(&GRADED[..])
                    );
                    let sparse = if graded {
                        p.solve_vibration_windows_with_nodal_density_gravity(
                            &GRADED,
                            gravity,
                            held_pressure,
                            &plan(expected_lambda),
                            &options(),
                            cx,
                        )
                    } else {
                        p.solve_vibration_windows_with_gravity(
                            &[6.0],
                            gravity,
                            held_pressure,
                            &plan(expected_lambda),
                            &options(),
                            cx,
                        )
                    }
                    .unwrap();
                    assert_eq!(sparse.equilibrium(), dense.equilibrium());
                    assert_eq!(sparse.stiffness_n_m().to_dense(), dense.stiffness_n_m());
                    assert_eq!(sparse.mass_kg().to_dense(), dense.mass_kg());
                    assert_eq!(sparse.spectrum().reports[0].expected, 1);
                    near(
                        sparse.spectrum().reports[0].modes[0].lambda,
                        expected_lambda,
                        3e-7,
                    );
                    assert!(core::ptr::eq(sparse.held_gravity_field().unwrap(), gravity));
                }
            }
        }
    });
}

#[test]
fn uniform_nodal_gravity_keeps_original_bits_and_zero_field_keeps_legacy_spectrum() {
    let law = material();
    let fixed = supports();
    let p = problem(&NODES, &TETS, &law, &fixed);
    with_cx(|cx| {
        let gravity = |x| field(x, 0.4, -1.0);
        let scalar = p
            .solve_vibration_with_gravity(&[6.0], &gravity, None, controls(), cx)
            .unwrap();
        let nodal = p
            .solve_vibration_with_nodal_density_gravity(&[6.0; 4], &gravity, None, controls(), cx)
            .unwrap();
        assert_eq!(scalar.equilibrium(), nodal.equilibrium());
        assert_eq!(scalar.stiffness_n_m(), nodal.stiffness_n_m());
        assert_eq!(scalar.mass_kg(), nodal.mass_kg());
        assert_eq!(
            scalar.modes()[0].lambda.to_bits(),
            nodal.modes()[0].lambda.to_bits()
        );
        let zero = |x| field(x, 0.0, 0.0);
        let retained = p
            .solve_vibration_with_gravity(&[6.0], &zero, None, controls(), cx)
            .unwrap();
        let legacy = p.solve_vibration(&[6.0], controls(), cx).unwrap();
        assert_eq!(retained.equilibrium(), legacy.equilibrium());
        assert_eq!(retained.stiffness_n_m(), legacy.stiffness_n_m());
        assert_eq!(retained.mass_kg(), legacy.mass_kg());
        assert_eq!(
            retained.modes()[0].lambda.to_bits(),
            legacy.modes()[0].lambda.to_bits()
        );
        assert!(retained.held_gravity_field().is_some());
        assert!(legacy.held_gravity_field().is_none());
    });
}

#[test]
fn destabilizing_gravity_retains_negative_squared_frequency_in_both_spectral_owners() {
    let law = material();
    let fixed = supports();
    let mut p = problem(&NODES, &TETS, &law, &fixed);
    p.settings.load_steps = 1;
    with_cx(|cx| {
        // A finite equilibrium gives the relative static residual a nonzero
        // force scale while gravity overwhelms material restoring stiffness.
        let stretch = 1.1;
        let gradient = 20.0;
        let (internal, material_stiffness) = material_force_stiffness(stretch);
        let constant = (internal - gradient * 0.1 * stretch) / 0.25;
        let gravity = |x| field(x, constant, gradient);
        let expected = (material_stiffness - gradient * 0.1) / 0.1;
        assert!(expected < 0.0);
        let dense = p
            .solve_vibration_with_gravity(&[6.0], &gravity, None, controls(), cx)
            .unwrap();
        near(
            dense.equilibrium().displacement_m[3][2],
            stretch - 1.0,
            2e-10,
        );
        near(dense.modes()[0].lambda, expected, 3e-8);
        assert!(dense.modes()[0].interval.1 < 0.0);
        let sparse = p
            .solve_vibration_windows_with_gravity(
                &[6.0],
                &gravity,
                None,
                &plan(expected),
                &options(),
                cx,
            )
            .unwrap();
        near(sparse.spectrum().reports[0].modes[0].lambda, expected, 3e-8);
        assert!(sparse.spectrum().reports[0].modes[0].interval.1 < 0.0);
    });
}

#[test]
fn gravity_vibration_preserves_original_field_density_pressure_and_resource_refusals() {
    let law = material();
    let fixed = supports();
    let p = problem(&NODES, &TETS, &law, &fixed);
    with_cx(|cx| {
        let outside = |_| {
            Err(HyperTetError::InvalidInput {
                what: "outside original gravity domain".into(),
            })
        };
        let asymmetric = |x| {
            let mut sample = field(x, 0.0, 0.0)?;
            sample.acceleration_gradient_s2[0][1] = 1.0;
            Ok(sample)
        };
        let nonfinite = |x| {
            let mut sample = field(x, 0.0, 0.0)?;
            sample.potential_j_kg = f64::NAN;
            Ok(sample)
        };
        for gravity in [&outside as &GravityField3<'_>, &asymmetric, &nonfinite] {
            assert!(
                p.solve_vibration_with_gravity(&[6.0], gravity, None, controls(), cx)
                    .is_err()
            );
            assert!(
                p.solve_vibration_windows_with_nodal_density_gravity(
                    &GRADED,
                    gravity,
                    None,
                    &plan(10.0),
                    &options(),
                    cx,
                )
                .is_err()
            );
        }
        let cancelled = |_| Err(HyperTetError::Cancelled);
        assert!(matches!(
            p.solve_vibration_with_gravity(&[6.0], &cancelled, None, controls(), cx),
            Err(VibrationError::Solid(HyperTetError::Cancelled))
        ));
        let gravity = |x| field(x, 0.0, 0.0);
        for rho in [&[][..], &[0.0], &[f64::NAN]] {
            assert!(
                p.solve_vibration_with_gravity(rho, &gravity, None, controls(), cx)
                    .is_err()
            );
        }
        assert!(
            p.solve_vibration_with_nodal_density_gravity(&[6.0; 3], &gravity, None, controls(), cx)
                .is_err()
        );
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, Default::default(), 4, cx).unwrap();
        let per_face = [0.0; 4];
        assert!(
            p.solve_vibration_with_gravity(
                &[6.0],
                &gravity,
                Some((&surface, TetPressure::PerFace(&per_face))),
                controls(),
                cx
            )
            .is_err()
        );
        assert!(
            p.solve_vibration_with_gravity(
                &[6.0],
                &gravity,
                None,
                VibrationSettings {
                    maximum_free_dofs: 0,
                    ..controls()
                },
                cx
            )
            .is_err()
        );
        let mut invalid_plan = plan(10.0);
        invalid_plan.maximum_bytes = 0;
        assert!(matches!(
            p.solve_vibration_windows_with_gravity(
                &[6.0],
                &gravity,
                None,
                &invalid_plan,
                &options(),
                cx,
            ),
            Err(VibrationError::Modal(_))
        ));
    });
}
