//! G1/G3: exact soft traces, standing waves and original heterogeneous work.
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::fixtures::masked_cube_grid;
use fs_rep_mesh::TetComplex;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};

fn brick(nx: usize) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] /= nx as f64;
        point[1] *= 0.2;
        point[2] *= 0.3;
    }
    (complex, positions)
}

fn end(geometry: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| geometry.positions()[v as usize][0] == x)
        })
        .collect()
}

fn velocity(faces: &[[u32; 3]], value: C64) -> Vec<VelocityTrace> {
    faces.iter().map(|&face| (face, [value; 3])).collect()
}

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn error(
    model: &HeterogeneousHelmholtzFe<'_, '_>,
    pressure: &[C64],
    expected: impl Fn(f64) -> C64,
) -> f64 {
    let geometry = model.geometry();
    let mut worst = 0.0_f64;
    for (i, point) in geometry.positions().iter().enumerate() {
        worst = worst.max((pressure[i] - expected(point[0])).abs());
    }
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let x = tet
            .iter()
            .map(|&v| geometry.positions()[v as usize][0])
            .sum::<f64>()
            / 4.0;
        let value = model.pressure_at(pressure, element, [0.25; 4]).unwrap();
        worst = worst.max((value - expected(x)).abs());
    }
    worst
}

#[test]
fn pressure_release_driven_duct_converges_to_standing_wave_and_stored_energy() {
    // With exp(+i wt), p(L)=0 and prescribed inward v(0), momentum gives
    // p(x)=i rho c v(0) sin(k(L-x))/cos(kL). No radiation end correction is
    // implied. MIT 6.013 chapter 13, section 13.2 and Fig. 13.2.3 give the
    // same closed/open standing-wave boundary conditions.
    let omega = 2.2;
    let k: f64 = omega / MEDIUM.sound_speed_m_s;
    let impedance = MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s;
    let drive = C64::new(0.6, -0.4);
    let amplitude = C64::new(-drive.im, drive.re).scale(impedance / k.cos());
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (complex, positions) = brick(nx);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let prepared = model
            .prepare_with_pressure_release(
                omega / core::f64::consts::TAU,
                &[],
                &end(&geometry, 1.0),
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let field = prepared
            .solve(&velocity(&end(&geometry, 0.0), drive), &mut || false)
            .unwrap();
        let relative = error(&model, &field.pressure_pa, |x| {
            amplitude.scale((k * (1.0 - x)).sin())
        }) / amplitude.abs();
        errors.push(relative);
        assert_eq!(field.pressure_pa.len(), geometry.ndof());
        assert!(field.backward_error < 5e-13);
        assert!(field.input_power_w.abs() < 2e-12);
        assert_eq!(field.boundary_dissipation_w, 0.0);
        for (i, &constrained) in prepared.pressure_release.constrained.iter().enumerate() {
            if constrained {
                assert_eq!(field.pressure_pa[i], C64::ZERO);
            }
        }
        if nx == 8 {
            let scale = MEDIUM.density_kg_m3 * 0.06 * drive.abs().powi(2) / (4.0 * k.cos().powi(2));
            let kinetic = scale * (0.5 + (2.0 * k).sin() / (4.0 * k));
            let compression = scale * (0.5 - (2.0 * k).sin() / (4.0 * k));
            assert!((field.kinetic_energy_j / kinetic - 1.0).abs() < 3e-4);
            assert!((field.compressional_energy_j / compression - 1.0).abs() < 3e-4);
        }
    }
    assert!(
        errors[1] < errors[0] / 5.0 && errors[2] < errors[1] / 5.0,
        "{errors:?}"
    );
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn pressure_release_preserves_material_interface_reflection_and_passive_boundary_work() {
    let right = Medium {
        density_kg_m3: 4.1,
        sound_speed_m_s: 1.6,
    };
    let (complex, positions) = brick(8);
    let geometry = HelmholtzFe::new(&complex, &positions, 3, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        5,
        &mut |_, point| Ok::<_, AcousticFeError>(if point[0] < 0.5 { MEDIUM } else { right }),
        &mut || false,
    )
    .unwrap();
    let omega = 1.2;
    let (k1, k2) = (
        omega / MEDIUM.sound_speed_m_s,
        omega / right.sound_speed_m_s,
    );
    let (z1, z2) = (
        MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s,
        right.density_kg_m3 * right.sound_speed_m_s,
    );
    let y = C64::new(0.17, -0.06);
    let drive = C64::new(0.3, 0.7);
    // Propagate pressure/velocity analytically from a unit-velocity release
    // end through two physical layers, then close v_imposed=v(0)+Y*p(0).
    let (s1, c1) = (0.5 * k1).sin_cos();
    let (s2, c2) = (0.5 * k2).sin_cos();
    let p0 = C64::new(0.0, z2 * s2 * c1 + z1 * c2 * s1);
    let v0 = c2 * c1 - z2 / z1 * s2 * s1;
    let outlet_velocity = drive / (C64::from_re(v0) + y * p0);
    let expected = |x: f64| {
        let value = if x <= 0.5 {
            z2 * s2 * (k1 * (0.5 - x)).cos() + z1 * c2 * (k1 * (0.5 - x)).sin()
        } else {
            z2 * (k2 * (1.0 - x)).sin()
        };
        C64::new(0.0, value) * outlet_velocity
    };
    let inlet = end(&geometry, 0.0);
    let prepared = model
        .prepare_with_pressure_release(
            omega / core::f64::consts::TAU,
            &[BoundaryAdmittance {
                faces: inlet.clone(),
                admittance_m_pa_s: y,
            }],
            &end(&geometry, 1.0),
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let field = prepared
        .solve(&velocity(&inlet, drive), &mut || false)
        .unwrap();
    assert!(error(&model, &field.pressure_pa, expected) < 3e-6);
    let boundary_loss = 0.5 * 0.06 * y.re * expected(0.0).abs().powi(2);
    assert!(boundary_loss > 0.0);
    assert!((field.boundary_dissipation_w / boundary_loss - 1.0).abs() < 3e-6);
    assert!(field.power_balance_residual_w().abs() < 2e-12 * boundary_loss);
}

#[test]
fn pressure_release_constrains_high_order_traces_and_keeps_original_volume_sources() {
    let (complex, positions) = brick(2);
    let geometry = HelmholtzFe::new(&complex, &positions, 4, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = end(&geometry, 1.0);
    let original = model
        .volume_velocity_source(
            6,
            &mut |_, p: [f64; 3]| Ok::<_, HeterogeneousError>(1.0 + p[0] + p[1] * p[2]),
            &mut || false,
        )
        .unwrap();
    let prepare = |faces: &[[u32; 3]], limits| {
        model
            .prepare_with_pressure_release(0.25, &[], faces, limits, &mut || false)
            .unwrap()
    };
    let first = prepare(&release, HarmonicLimits::default());
    let free = first.pressure_release.free_dofs.len();
    assert!(
        first
            .pressure_release
            .constrained
            .iter()
            .enumerate()
            .any(|(i, &fixed)| fixed && i >= positions.len()),
        "edge/face traces must be constrained"
    );
    // Budget the matrix that is actually factored, not a full zero-padded one.
    let limited = HarmonicLimits {
        max_dofs: free,
        max_dense_bytes: 2 * size_of::<C64>() * free * free,
        ..HarmonicLimits::default()
    };
    let mut permuted = release.clone();
    permuted.reverse();
    for face in &mut permuted {
        face.swap(0, 2);
    }
    let second = prepare(&permuted, limited);
    let drive = C64::new(0.7, -0.3);
    let a = first
        .solve_with_volume_sources(&[], &[(&original, drive)], &mut || false)
        .unwrap();
    let b = second
        .solve_with_volume_sources(&[], &[(&original, drive)], &mut || false)
        .unwrap();
    assert_eq!(a.pressure_pa, b.pressure_pa);
    assert!(a.pressure_pa.iter().any(|p| p.abs() > 0.0));
    assert!(a.input_power_w.abs() < 1e-12);
    assert!((a.volume_source_power_w - a.input_power_w).abs() < 1e-12);
    for face in release {
        let rule = geometry.face_rule(face, &tri_quadrature(6)).unwrap();
        for &dof in &rule.dofs {
            assert_eq!(a.pressure_pa[dof], C64::ZERO);
        }
    }
    let too_small = HarmonicLimits {
        max_dense_bytes: limited.max_dense_bytes - 1,
        ..limited
    };
    assert!(
        model
            .prepare_with_pressure_release(0.25, &[], &permuted, too_small, &mut || false)
            .is_err()
    );
}

#[test]
fn pressure_release_rejects_conflicting_laws_and_preserves_reusable_factor_after_cancellation() {
    let (complex, positions) = brick(2);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = end(&geometry, 1.0);
    let conflict = BoundaryAdmittance {
        faces: release.clone(),
        admittance_m_pa_s: C64::ZERO,
    };
    assert!(
        model
            .prepare_with_pressure_release(
                0.25,
                &[conflict],
                &release,
                HarmonicLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let prepared = model
        .prepare_with_pressure_release(0.25, &[], &release, HarmonicLimits::default(), &mut || {
            false
        })
        .unwrap();
    let conflict = velocity(&release, C64::ZERO);
    assert!(prepared.solve(&conflict, &mut || false).is_err());
    assert!(
        prepared
            .solve_with_volume_sources(&conflict, &[], &mut || false)
            .is_err()
    );
    let drive = velocity(&end(&geometry, 0.0), C64::ONE);
    let expected = prepared.solve(&drive, &mut || false).unwrap();
    let mut calls = 0;
    assert!(matches!(
        prepared.solve(&drive, &mut || {
            calls += 1;
            calls == 20
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(
        prepared.solve(&drive, &mut || false).unwrap().pressure_pa,
        expected.pressure_pa
    );
    assert!(matches!(
        model.prepare_with_pressure_release(
            0.25,
            &[],
            &release,
            HarmonicLimits::default(),
            &mut || true
        ),
        Err(HeterogeneousError::Cancelled)
    ));
    // Existing rigid/impedance entry point is exactly the empty constraint map.
    let a = model
        .prepare(0.25, &[], HarmonicLimits::default(), &mut || false)
        .unwrap()
        .solve(&drive, &mut || false)
        .unwrap();
    let b = model
        .prepare_with_pressure_release(0.25, &[], &[], HarmonicLimits::default(), &mut || false)
        .unwrap()
        .solve(&drive, &mut || false)
        .unwrap();
    assert_eq!(a.pressure_pa, b.pressure_pa);
    assert_eq!(a.backward_error.to_bits(), b.backward_error.to_bits());
}
