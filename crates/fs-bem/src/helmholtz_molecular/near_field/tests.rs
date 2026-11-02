//! G0/G1/G3 source integration. Synthetic caloric data isolate the propagation
//! law; these are numerical fixtures, not measured atmospheric validation.
use super::MolecularNearFieldOptions;
use crate::{
    helmholtz::{
        Formulation, HelmholtzError, Medium,
        laplace::complex::ComplexOptions,
        molecular::{MolecularError, MolecularExterior, MolecularRadiation},
        solve_radiation,
    },
    near_field::{Geometry, Options, VectorObservation},
    panel3d::SpherePanels,
};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::c64::C64;
use std::sync::Arc;

fn material(relaxing: bool) -> MolecularExterior {
    let text = if relaxing {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

fn surface() -> SpherePanels {
    let points = [
        [0.0, 0.0, -0.002],
        [0.1, 0.0, -0.002],
        [0.1, 0.1, -0.002],
        [0.0, 0.1, -0.002],
        [0.0, 0.0, 0.002],
        [0.1, 0.0, 0.002],
        [0.1, 0.1, 0.002],
        [0.0, 0.1, 0.002],
    ];
    let triangles = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    SpherePanels::from_triangles(triangles.map(|t| t.map(|i| points[i])).to_vec()).unwrap()
}

fn controls() -> MolecularNearFieldOptions {
    MolecularNearFieldOptions {
        quadrature: Options {
            relative_tolerance: 1e-8,
            maximum_depth: 14,
            maximum_kernel_evaluations: 5_000_000,
        },
        maximum_row_bytes: 1_000_000,
    }
}

fn solve(
    material: &MolecularExterior,
    surface: &SpherePanels,
    omega: f64,
    velocities: &[&[C64]],
) -> Vec<MolecularRadiation> {
    let n = surface.areas().len();
    material
        .solve_batch(
            surface,
            omega,
            velocities,
            Formulation::PlainCbie,
            &ComplexOptions {
                max_panels: n,
                max_bytes: 32 * n * n,
                max_work: 512 * (n * n) as u64,
            },
            || false,
        )
        .unwrap()
}

fn close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-30),
        "{a:?} != {b:?}"
    );
}

fn unchanged(a: &VectorObservation, b: &VectorObservation) {
    assert_eq!(a.scalar.pressure, b.scalar.pressure);
    assert_eq!(a.particle_velocity_m_s, b.particle_velocity_m_s);
    assert_eq!(
        a.scalar.quadrature_error_estimate_pa,
        b.scalar.quadrature_error_estimate_pa
    );
    assert_eq!(
        a.quadrature_error_estimate_m_s,
        b.quadrature_error_estimate_m_s
    );
}

#[test]
fn empty_reservoir_adaptive_receiver_recovers_original_real_pressure_and_velocity() {
    let surface = surface();
    let material = material(false);
    let omega = 650.0;
    let drive: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::new(n[2], 0.3 * n[0]))
        .collect();
    let solved = solve(&material, &surface, omega, &[&drive]).remove(0);
    let k = solved.wavenumber_per_m();
    assert_eq!(k.im, 0.0);
    let medium = Medium {
        density: material.source().gas().density,
        sound_speed: omega / k.re,
    };
    let native = solve_radiation(&surface, k.re, medium, &drive, Formulation::PlainCbie).unwrap();
    let points = [[0.031, 0.043, 0.027], [0.059, 0.024, -0.026]];
    let original = Geometry::new(&surface, &points, 0.01)
        .unwrap()
        .prepare_velocity(k.re, medium, controls().quadrature)
        .unwrap()
        .evaluate_velocity(&native)
        .unwrap();
    let prepared = solved
        .prepare_near_field(&surface, &points, 0.01, controls(), || false)
        .unwrap();
    let actual = prepared.evaluate(&solved, || false).unwrap();
    assert!(prepared.kernel_evaluations() > 80 * surface.areas().len() * points.len());
    for i in 0..points.len() {
        close(
            actual.scalar.pressure[i],
            original.scalar.pressure[i],
            original.scalar.pressure[i].abs(),
            2e-10,
        );
        let velocity_scale = original.particle_velocity_m_s[i]
            .iter()
            .map(|u| u.norm_sq())
            .sum::<f64>()
            .sqrt();
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[i][axis],
                original.particle_velocity_m_s[i][axis],
                velocity_scale,
                2e-10,
            );
        }
        assert!(actual.scalar.quadrature_error_estimate_pa[i].is_finite());
        assert!(actual.scalar.quadrature_error_estimate_pa[i] >= 0.0);
        assert!(
            actual.quadrature_error_estimate_m_s[i]
                .iter()
                .all(|e| e.is_finite() && *e >= 0.0)
        );
    }
    let repeated = solved
        .observe_near_field(&surface, &points, 0.01, controls(), || false)
        .unwrap();
    unchanged(&actual, &repeated);
}

#[test]
fn relaxing_source_rows_reuse_coherent_fields_and_preserve_physical_euler_momentum() {
    let surface = surface();
    let material = material(true);
    let omega = 650.0;
    let first: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::new(n[2], 0.3 * n[0]))
        .collect();
    let second: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::new(0.5 + 0.2 * n[0], -0.1 * n[1]))
        .collect();
    let (a, b) = (C64::new(1.0, -0.3), C64::new(-0.2, 0.7));
    let combined: Vec<_> = first
        .iter()
        .zip(&second)
        .map(|(&x, &y)| a * x + b * y)
        .collect();
    let solved = solve(&material, &surface, omega, &[&first, &second, &combined]);
    assert!(solved[0].wavenumber_per_m().im > 0.0);
    let point = [0.031, 0.043, 0.027];
    let prepared = solved[0]
        .prepare_near_field(&surface, &[point], 0.01, controls(), || false)
        .unwrap();
    let fields: Vec<_> = solved
        .iter()
        .map(|source| prepared.evaluate(source, || false).unwrap())
        .collect();
    let pressure_scale =
        a.abs() * fields[0].scalar.pressure[0].abs() + b.abs() * fields[1].scalar.pressure[0].abs();
    close(
        fields[2].scalar.pressure[0],
        a * fields[0].scalar.pressure[0] + b * fields[1].scalar.pressure[0],
        pressure_scale,
        2e-11,
    );
    for axis in 0..3 {
        let scale = a.abs() * fields[0].particle_velocity_m_s[0][axis].abs()
            + b.abs() * fields[1].particle_velocity_m_s[0][axis].abs();
        close(
            fields[2].particle_velocity_m_s[0][axis],
            a * fields[0].particle_velocity_m_s[0][axis]
                + b * fields[1].particle_velocity_m_s[0][axis],
            scale,
            2e-11,
        );
    }
    unchanged(
        &fields[0],
        &prepared.evaluate(&solved[0].clone(), || false).unwrap(),
    );

    // Independently differentiate scalar receiver pressure. This specifically
    // tests conversion with the ORIGINAL real rho*omega in a dispersive gas.
    // A complex k*c_eq replacement produces a different physical velocity.
    let h = 4e-4;
    let mut points = Vec::new();
    for axis in 0..3 {
        for step in [-1.0, -0.5, 0.5, 1.0] {
            let mut x = point;
            x[axis] += step * h;
            points.push(x);
        }
    }
    let scalar = solved[2]
        .observe_near_field(&surface, &points, 0.01, controls(), || false)
        .unwrap();
    let momentum = material.source().gas().density * omega;
    let scale = fields[2].particle_velocity_m_s[0]
        .iter()
        .map(|u| u.norm_sq())
        .sum::<f64>()
        .sqrt();
    for axis in 0..3 {
        let p = &scalar.scalar.pressure[4 * axis..4 * axis + 4];
        let derivative = (p[0] - p[1].scale(8.0) + p[2].scale(8.0) - p[3]).scale(1.0 / (6.0 * h));
        let expected = derivative * C64::new(0.0, -1.0 / momentum);
        close(
            fields[2].particle_velocity_m_s[0][axis],
            expected,
            scale,
            5e-6,
        );
    }
}

#[test]
fn adaptive_source_binding_budgets_and_cancellation_preserve_reusable_rows() {
    let surface = surface();
    let material = material(true);
    let drive: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::from_re(n[2]))
        .collect();
    let solved = solve(&material, &surface, 650.0, &[&drive]).remove(0);
    let points = [[0.031, 0.043, 0.027], [0.059, 0.024, -0.026]];
    let mut preparation_polls = 0;
    let prepared = solved
        .prepare_near_field(&surface, &points, 0.01, controls(), || {
            preparation_polls += 1;
            false
        })
        .unwrap();
    let expected = prepared.evaluate(&solved, || false).unwrap();
    for changed in 0..7 {
        let mut wrong = solved.clone();
        match changed {
            0 => wrong.source = Arc::new(wrong.source().clone()),
            1 => wrong.omega *= 1.01,
            2 => wrong.k.re *= 1.01,
            3 => wrong.fingerprint ^= 1,
            4 => {
                wrong.pressure.pop();
            }
            5 => wrong.pressure[0] = C64::new(f64::NAN, 0.0),
            6 => wrong.velocity[0] = C64::new(0.0, f64::INFINITY),
            _ => unreachable!(),
        }
        assert!(
            prepared.evaluate(&wrong, || false).is_err(),
            "accepted mismatched case {changed}"
        );
    }
    for options in [
        MolecularNearFieldOptions {
            maximum_row_bytes: 0,
            ..controls()
        },
        MolecularNearFieldOptions {
            quadrature: Options {
                maximum_kernel_evaluations: 80,
                ..controls().quadrature
            },
            ..controls()
        },
        MolecularNearFieldOptions {
            quadrature: Options {
                maximum_depth: 0,
                ..controls().quadrature
            },
            ..controls()
        },
    ] {
        assert!(
            solved
                .prepare_near_field(&surface, &points, 0.01, options, || false)
                .is_err()
        );
    }
    assert!(
        solved
            .prepare_near_field(&surface, &[[0.05, 0.05, 0.0]], 0.01, controls(), || false)
            .is_err()
    );
    let other = SpherePanels::icosphere(0.1, 0).unwrap();
    assert!(
        solved
            .prepare_near_field(&other, &points, 0.01, controls(), || false)
            .is_err()
    );
    assert!(preparation_polls > 2);
    let mut polls = 0;
    assert!(matches!(
        solved.prepare_near_field(&surface, &points, 0.01, controls(), || {
            polls += 1;
            polls >= preparation_polls / 2
        }),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let mut evaluation_polls = 0;
    prepared
        .evaluate(&solved, || {
            evaluation_polls += 1;
            false
        })
        .unwrap();
    assert!(evaluation_polls > 1);
    let mut polls = 0;
    assert!(matches!(
        prepared.evaluate(&solved, || {
            polls += 1;
            polls >= evaluation_polls
        }),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    unchanged(&expected, &prepared.evaluate(&solved, || false).unwrap());
}
