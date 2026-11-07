//! Observe actual pressure-jump BEM solutions with the shared adaptive P1
//! Green derivatives. These checks concern receiver accuracy for the solved
//! mesh, not convergence of the boundary discretization to a continuum screen.
use fs_bem::helmholtz::bc::{Admittance, SurfaceBc};
use fs_bem::helmholtz::{HelmholtzError, Medium, molecular::MolecularExterior};
use fs_bem::near_field::{FirstOrder, Options, VectorObservation};
use fs_bem::reciprocal_radiation::open_screen::{
    MolecularScreenPortLimits, OpenOptions, OpenPanels, PreparedScreen, ScreenComplexValueLimits,
    ScreenPortLimits,
};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::c64::C64;

fn surface(curved: bool) -> OpenPanels {
    let mut vertices = Vec::new();
    for j in 0..4 {
        for i in 0..4 {
            let x = (i as f64 - 1.5) * 0.25;
            let y = (j as f64 - 1.5) * 0.25;
            vertices.push([x, y, if curved { 0.125 * x * y } else { 0.0 }]);
        }
    }
    let mut triangles = Vec::new();
    for j in 0..3 {
        for i in 0..3 {
            let a = 4 * j + i;
            triangles.push([a, a + 1, a + 5]);
            triangles.push([a, a + 5, a + 4]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}

fn options() -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 4,
        angular_refinement: 1,
        maximum_unknowns: 16,
        maximum_cache_bytes: 1_000_000,
        maximum_work: u64::MAX,
        workers: 1,
    }
}

fn quadrature() -> Options {
    Options {
        relative_tolerance: 1e-10,
        maximum_depth: 14,
        maximum_kernel_evaluations: 20_000_000,
    }
}

fn limits(surface: &OpenPanels) -> ScreenPortLimits {
    ScreenPortLimits {
        maximum_ports: 2,
        maximum_matrix_entries: 4,
        maximum_jump_entries: 2 * surface.vertices().len(),
        maximum_work: u64::MAX,
    }
}

fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= tolerance * scale.max(1e-300),
        "error {error}, scale {scale}"
    );
}

fn close_fields(actual: &VectorObservation, expected: &VectorObservation, tolerance: f64) {
    close(
        &actual.scalar.pressure,
        &expected.scalar.pressure,
        tolerance,
    );
    for (a, e) in actual
        .particle_velocity_m_s
        .iter()
        .zip(&expected.particle_velocity_m_s)
    {
        close(a, e, tolerance);
    }
}

#[test]
fn screen_receivers_retain_two_sided_euler_velocity_and_pressure_null_pickup() {
    let surface = surface(false);
    let medium = Medium::air();
    let k = 1.2;
    let owner = PreparedScreen::new(&surface, 2.0, options(), || false).unwrap();
    let drive: Vec<_> = surface
        .vertices()
        .iter()
        .map(|x| C64::new(1.0 + 0.3 * x[0], 0.2 * x[1]))
        .collect();
    let source = owner
        .solve(
            k,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    let points = [[0.18, -0.13, 0.46], [0.18, -0.13, -0.46], [0.7, 0.12, 0.0]];
    let rows = source
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap();
    let actual = rows.evaluate(|| false).unwrap();
    assert_eq!(rows.points_m(), points);
    close(
        &actual.scalar.pressure,
        &source
            .pressure_at(&points, 1e-5, quadrature(), || false)
            .unwrap(),
        2e-9,
    );
    let p = &actual.scalar.pressure;
    let v = &actual.particle_velocity_m_s;
    close(&[p[1]], &[-p[0]], 2e-12);
    close(&v[1], &[-v[0][0], -v[0][1], v[0][2]], 2e-12);

    // Differentiate the independently integrated scalar pressure, on BOTH
    // sides, through three second-order finite-difference rungs.
    let rho_omega = medium.density * medium.sound_speed * k;
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let mut displaced = Vec::new();
        for &point in &points[..2] {
            for axis in 0..3 {
                for sign in [-1.0, 1.0] {
                    let mut x = point;
                    x[axis] += sign * h;
                    displaced.push(x);
                }
            }
        }
        let pressure = source
            .pressure_at(&displaced, 1e-5, quadrature(), || false)
            .unwrap();
        let mut error2 = 0.0;
        let mut scale2 = 0.0;
        for receiver in 0..2 {
            for axis in 0..3 {
                let index = 6 * receiver + 2 * axis;
                let fd = (pressure[index + 1] - pressure[index]).scale(0.5 / h);
                let analytic = v[receiver][axis] * C64::new(0.0, rho_omega);
                error2 += (fd - analytic).abs().powi(2);
                scale2 += analytic.abs().powi(2);
            }
        }
        errors.push((error2 / scale2).sqrt());
    }
    assert!(
        errors[1] < 0.28 * errors[0] && errors[2] < 0.28 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-3, "{errors:?}");
    println!("open-screen Euler derivative relative errors: {errors:?}");

    // A planar double layer has a pressure node in its plane outside the
    // screen. The normal velocity and a directional microphone remain live.
    assert!(p[2].abs() < 1e-12 * p[0].abs());
    assert!(v[2][2].abs() > 1e-4);
    let figure_eight = FirstOrder::new(0.0, [0.0, 0.0, 1.0])
        .unwrap()
        .observe(p[2], v[2], medium)
        .unwrap();
    let cardioid = FirstOrder::new(0.5, [0.0, 0.0, 1.0])
        .unwrap()
        .observe(p[2], v[2], medium)
        .unwrap();
    close(
        &[figure_eight],
        &[v[2][2].scale(-medium.density * medium.sound_speed)],
        1e-14,
    );
    close(&[cardioid], &[figure_eight.scale(0.5)], 1e-14);
    assert!(
        actual
            .scalar
            .quadrature_error_estimate_pa
            .iter()
            .all(|x| x.is_finite() && *x >= 0.0)
    );
    assert!(
        actual
            .quadrature_error_estimate_m_s
            .iter()
            .flatten()
            .all(|x| x.is_finite() && *x >= 0.0)
    );
}

#[test]
fn screen_port_receivers_keep_coherent_fields_and_reuse_quadrature() {
    let surface = surface(true);
    let medium = Medium::air();
    let k = 1.2;
    let owner = PreparedScreen::new(&surface, 2.0, options(), || false).unwrap();
    let nodal: Vec<Vec<f64>> = vec![
        vec![1.0; surface.vertices().len()],
        surface
            .vertices()
            .iter()
            .map(|x| 0.3 + x[0] - 0.4 * x[1])
            .collect(),
    ];
    let traces: Vec<Vec<_>> = nodal
        .iter()
        .map(|row| {
            surface
                .triangles()
                .iter()
                .map(|ids| ids.map(|v| row[v]))
                .collect()
        })
        .collect();
    let refs: Vec<_> = traces.iter().map(Vec::as_slice).collect();
    let source = owner
        .solve_velocity_ports(k, medium, &refs, limits(&surface), || false)
        .unwrap();
    let points = [[0.2, -0.1, 0.45], [0.1, 0.2, -0.6]];
    let rows = source
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap();
    let work = rows.kernel_evaluations();
    assert!(work >= 80 * points.len() * surface.triangles().len());
    assert_eq!(rows.points_m(), points);
    let rates = [C64::new(0.6, -0.25), C64::new(-0.4, 0.7)];
    let actual = rows.evaluate(&rates, || false).unwrap();
    let drive: Vec<_> = (0..surface.vertices().len())
        .map(|v| rates[0].scale(nodal[0][v]) + rates[1].scale(nodal[1][v]))
        .collect();
    let combined = owner
        .solve(
            k,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    let expected = combined
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap()
        .evaluate(|| false)
        .unwrap();
    close_fields(&actual, &expected, 3e-12);
    let one = rows.evaluate(&[rates[0], C64::ZERO], || false).unwrap();
    let two = rows.evaluate(&[C64::ZERO, rates[1]], || false).unwrap();
    let coherent = actual.scalar.pressure[0].abs().powi(2);
    let incoherent = one.scalar.pressure[0].abs().powi(2) + two.scalar.pressure[0].abs().powi(2);
    assert!((coherent - incoherent).abs() > 0.01 * incoherent);
    let zero = rows.evaluate(&[C64::ZERO; 2], || false).unwrap();
    assert!(
        zero.scalar
            .pressure
            .iter()
            .chain(zero.particle_velocity_m_s.iter().flatten())
            .all(|&value| value == C64::ZERO)
    );
    let replay = rows.evaluate(&rates, || false).unwrap();
    assert_eq!(replay.scalar.pressure, actual.scalar.pressure);
    assert_eq!(replay.particle_velocity_m_s, actual.particle_velocity_m_s);
    assert_eq!(rows.kernel_evaluations(), work);
    assert!(rows.evaluate(&[], || false).is_err());
    assert!(
        rows.evaluate(&[C64::new(f64::NAN, 0.0), C64::ZERO], || false)
            .is_err()
    );
    assert!(matches!(
        rows.evaluate(&rates, || true),
        Err(HelmholtzError::Cancelled)
    ));
}

#[test]
fn screen_receivers_refuse_unresolved_geometry_resources_and_cancelled_work() {
    let surface = surface(false);
    let medium = Medium::air();
    let owner = PreparedScreen::new(&surface, 2.0, options(), || false).unwrap();
    let drive = vec![C64::from_re(1.0); surface.vertices().len()];
    let bc = SurfaceBc {
        drive_velocity: &drive,
        admittance: Admittance::Rigid,
    };
    let source = owner.solve(1.2, medium, bc, || false).unwrap();
    let points = [[0.1, 0.2, 0.5], [-0.2, 0.1, -0.5]];
    for refused in [
        vec![[0.0; 3]],
        vec![surface.vertices()[0]],
        vec![[f64::NAN, 0.0, 1.0]],
        vec![points[0]; fs_bem::near_field::MAX_RECEIVERS + 1],
    ] {
        assert!(
            source
                .prepare_receivers(&refused, 1e-5, quadrature(), 1_000_000, || false)
                .is_err()
        );
    }
    assert!(
        source
            .prepare_receivers(&points, 1e-5, quadrature(), 0, || false)
            .is_err()
    );
    assert!(
        source
            .prepare_receivers(&points, 0.0, quadrature(), 1_000_000, || false)
            .is_err()
    );
    let short = Options {
        maximum_kernel_evaluations: 80 * surface.triangles().len() * points.len() - 1,
        ..quadrature()
    };
    assert!(
        source
            .prepare_receivers(&points, 1e-5, short, 1_000_000, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        source.prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || {
            polls += 1;
            polls >= 20
        }),
        Err(HelmholtzError::Cancelled)
    ));
    let rows = source
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap();
    let original = rows.evaluate(|| false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        rows.evaluate(|| {
            polls += 1;
            polls >= 4
        }),
        Err(HelmholtzError::Cancelled)
    ));
    let replay = rows.evaluate(|| false).unwrap();
    assert_eq!(replay.scalar.pressure, original.scalar.pressure);
    assert_eq!(replay.particle_velocity_m_s, original.particle_velocity_m_s);
    let damped = owner
        .solve_laplace(C64::new(0.1, -1.2 * medium.sound_speed), medium, bc, || {
            false
        })
        .unwrap();
    assert!(
        damped
            .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
            .is_err()
    );
}

#[test]
fn screen_receivers_match_the_zero_relaxation_material_owner() {
    let law = read_relaxation(&mut "declared 1 synthetic 0".split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    let material = MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap();
    let omega = 450.0;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    assert_eq!(k.im, 0.0);
    let surface = surface(true);
    let owner = PreparedScreen::new(&surface, 2.0, options(), || false).unwrap();
    let traces: Vec<Vec<_>> = [0.0, 1.0]
        .iter()
        .map(|&slope| {
            surface
                .triangles()
                .iter()
                .map(|ids| ids.map(|v| 1.0 + slope * surface.vertices()[v][0]))
                .collect()
        })
        .collect();
    let refs: Vec<_> = traces.iter().map(Vec::as_slice).collect();
    let plain = owner
        .solve_velocity_ports(
            k.re,
            Medium {
                density: material.source().gas().density,
                sound_speed: omega / k.re,
            },
            &refs,
            limits(&surface),
            || false,
        )
        .unwrap();
    let molecular = owner
        .solve_molecular_velocity_ports(
            &material,
            omega,
            &refs,
            MolecularScreenPortLimits {
                ports: limits(&surface),
                operator: ScreenComplexValueLimits {
                    maximum_matrix_entries: 1024,
                    maximum_workspace_bytes: 4_000_000,
                    maximum_work: u64::MAX,
                    maximum_growth_exponent: 8.0,
                },
                maximum_workspace_bytes: 4_000_000,
            },
            || false,
        )
        .unwrap();
    let points = [[0.2, 0.1, 0.5], [-0.1, 0.15, -0.4]];
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    let actual = plain
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap()
        .evaluate(&rates, || false)
        .unwrap();
    let reference = molecular
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap()
        .evaluate(&rates, || false)
        .unwrap();
    close_fields(&actual, &reference, 3e-12);
}
