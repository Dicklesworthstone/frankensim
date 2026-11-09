//! G0/G1/G2/G3/G5: physical soft-screen boundaries, a static-disc oracle,
//! optical theorem, finite-source reciprocity and original Euler gradients.
use super::*;
use crate::{
    helmholtz::bc::{IncidentSource, Monopole, PlaneWave},
    near_field,
};

fn medium() -> Medium {
    Medium {
        density: 1.7,
        sound_speed: 287.0,
    }
}
fn options() -> OpenOptions {
    OpenOptions {
        edge_level: 1,
        maximum_edge_level: 4,
        triangle_order: 4,
        angular_refinement: 1,
        maximum_unknowns: 1024,
        maximum_cache_bytes: 64 * 1024 * 1024,
        maximum_work: u64::MAX,
        workers: 1,
    }
}
fn limits() -> SoftScreenScatteringLimits {
    SoftScreenScatteringLimits {
        maximum_cases: 8,
        maximum_density_entries: 8192,
        maximum_source_point_products: 10_000_000,
        minimum_source_clearance_m: 1e-5,
    }
}
fn plane(amplitude: C64) -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: amplitude,
    })
}
// Unit disc with quadratic edge grading. Refining changes the physical
// polygon and edge-density representation independently of the Green kernel.
fn disc(radial: usize, azimuth: usize) -> OpenPanels {
    let mut vertices = vec![[0.0; 3]];
    for r in 1..=radial {
        let radius = det::sin(PI * r as f64 / (2.0 * radial as f64));
        for a in 0..azimuth {
            let phi = TAU * a as f64 / azimuth as f64;
            vertices.push([radius * det::cos(phi), radius * det::sin(phi), 0.0]);
        }
    }
    let at = |r: usize, a: usize| 1 + (r - 1) * azimuth + a % azimuth;
    let mut triangles = Vec::new();
    for a in 0..azimuth {
        triangles.push([0, at(1, a), at(1, a + 1)]);
    }
    for r in 1..radial {
        for a in 0..azimuth {
            triangles.push([at(r, a), at(r + 1, a), at(r + 1, a + 1)]);
            triangles.push([at(r, a), at(r + 1, a + 1), at(r, a + 1)]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}
fn receiver_options() -> near_field::Options {
    near_field::Options {
        relative_tolerance: 1e-8,
        maximum_depth: 12,
        maximum_kernel_evaluations: 2_000_000,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-20),
        "{a:?} versus {b:?}"
    );
}

#[test]
fn static_disc_capacity_and_axial_field_converge_without_clamping_edge_density() {
    // Independent Laplace disc: sigma=-4/(pi*sqrt(a²-r²)), total charge=-8a,
    // p_scattered(0,0,z)=-(2/pi)*atan(a/|z|). These follow by integrating the
    // equilibrium 1/sqrt(a²-r²) density against 1/(4*pi*R), not from BEM.
    let k = 1e-4;
    let exact = C64::from_re(-2.0 / PI * det::atan2(1.0, 0.7));
    let mut charge_errors = Vec::new();
    let mut field_errors = Vec::new();
    for (radial, azimuth) in [(2, 12), (3, 18), (4, 24)] {
        let screen = disc(radial, azimuth);
        let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
        let frequency = prepared
            .prepare_frequency(k, medium(), 64 * 1024 * 1024, || false)
            .unwrap();
        let solution = frequency
            .scatter(plane(C64::ONE), limits(), || false)
            .unwrap();
        let charge = solution
            .density
            .iter()
            .zip(&prepared.facets)
            .fold(C64::ZERO, |sum, (&v, f)| sum + v.scale(f.area));
        charge_errors.push((charge + C64::from_re(8.0)).abs() / 8.0);
        let receivers = solution
            .prepare_receivers(
                &[[0.0, 0.0, 0.7], [0.0, 0.0, -0.7]],
                1e-4,
                receiver_options(),
                2_000_000,
                || false,
            )
            .unwrap();
        let field = receivers.evaluate_scattered(|| false).unwrap();
        field_errors.push((field.scalar.pressure[0] - exact).abs() / exact.abs());
        close(field.scalar.pressure[0], field.scalar.pressure[1], 1e-12);
        close(
            field.particle_velocity_m_s[0][2],
            -field.particle_velocity_m_s[1][2],
            1e-12,
        );
        assert!(solution.density.last().unwrap().abs() > solution.density[0].abs());
        assert!(solution.backward_residual < 1e-12);
    }
    eprintln!("soft-disc charge errors {charge_errors:?}; axial-field errors {field_errors:?}");
    assert!(charge_errors.windows(2).all(|e| e[1] < e[0]));
    assert!(field_errors.windows(2).all(|e| e[1] < e[0]));
    assert!(charge_errors[2] < 0.04 && field_errors[2] < 0.04);
}

#[test]
fn optical_theorem_and_both_sides_retain_coherent_scattered_power() {
    let screen = disc(2, 12);
    let k = 0.7;
    let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 2_000_000, || false)
        .unwrap();
    let amplitude = C64::new(0.8, -0.6);
    let solution = frequency
        .scatter(plane(amplitude), limits(), || false)
        .unwrap();
    let forward = solution
        .scattered_far_field([0.0, 0.0, 1.0], || false)
        .unwrap();
    let backward = solution
        .scattered_far_field([0.0, 0.0, -1.0], || false)
        .unwrap();
    close(forward, backward, 1e-12);
    let extinction =
        2.0 * PI / (medium().density * medium().sound_speed * k) * (amplitude.conj() * forward).im;
    assert!(solution.scattered_power_w > 0.0);
    assert!((solution.scattered_power_w - extinction).abs() < 2e-12 * extinction);
    // Independent, finer outgoing angular rule, with both antipodal partners.
    let mut power = 0.0;
    for (direction, weight) in directions(2) {
        for d in [direction, direction.map(|x| -x)] {
            let f = solution.scattered_far_field(d, || false).unwrap();
            power += weight * f.abs().powi(2) / (2.0 * medium().density * medium().sound_speed);
        }
    }
    assert!((power - extinction).abs() < 1e-9 * extinction);
}

#[test]
fn finite_source_reciprocity_superposition_and_euler_velocity() {
    let screen = disc(2, 12);
    let k = 0.5;
    let prepared = PreparedSoftScreen::new(
        &screen,
        k,
        OpenOptions {
            triangle_order: 8,
            ..options()
        },
        || false,
    )
    .unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 2_000_000, || false)
        .unwrap();
    let a = [Monopole {
        position_m: [0.3, -0.1, 1.2],
        volume_velocity_m3_s: C64::new(0.001, -0.002),
    }];
    let b = [Monopole {
        position_m: [-0.2, 0.4, -1.1],
        volume_velocity_m3_s: a[0].volume_velocity_m3_s,
    }];
    let cases = [IncidentField::Monopoles(&a), IncidentField::Monopoles(&b)];
    let solved = frequency.scatter_many(&cases, limits(), || false).unwrap();
    let evaluate = |source: &SoftScreenScattering<'_, '_, '_>, point| {
        source
            .prepare_receivers(&[point], 1e-4, receiver_options(), 2_000_000, || false)
            .unwrap()
            .evaluate_scattered(|| false)
            .unwrap()
            .scalar
            .pressure[0]
    };
    close(
        evaluate(&solved[0], b[0].position_m),
        evaluate(&solved[1], a[0].position_m),
        2e-7,
    );
    let terms = [
        IncidentSource::Monopole(a[0]),
        IncidentSource::Monopole(b[0]),
    ];
    let combined = frequency
        .scatter(IncidentField::Sources(&terms), limits(), || false)
        .unwrap();
    for i in 0..combined.density.len() {
        close(
            combined.density[i],
            solved[0].density[i] + solved[1].density[i],
            1e-11,
        );
    }
    let center = [0.27, 0.31, 0.8];
    let h = 2e-4;
    let mut points = vec![center];
    for axis in 0..3 {
        let mut minus = center;
        let mut plus = center;
        minus[axis] -= h;
        plus[axis] += h;
        points.extend([minus, plus]);
    }
    let rows = combined
        .prepare_receivers(&points, 1e-4, receiver_options(), 2_000_000, || false)
        .unwrap();
    let total = rows.evaluate(|| false).unwrap();
    let inverse = C64::new(0.0, -1.0 / (k * medium().sound_speed * medium().density));
    for axis in 0..3 {
        let derivative = (total.scalar.pressure[2 + 2 * axis]
            - total.scalar.pressure[1 + 2 * axis])
            .scale(0.5 / h);
        close(
            total.particle_velocity_m_s[0][axis],
            derivative * inverse,
            2e-7,
        );
    }
}

#[test]
fn free_edge_only_triangles_are_unknowns_and_budgets_refuse_before_work() {
    let screen = OpenPanels::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    assert_eq!(screen.unknowns(), 0); // No hard-screen P1 degrees of freedom.
    let prepared = PreparedSoftScreen::new(&screen, 0.1, options(), || false).unwrap();
    assert_eq!(prepared.unknowns(), 1); // Soft-screen P0 density remains live.
    assert!(matches!(
        PreparedSoftScreen::new(&screen, 0.1, options(), || true),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(
        PreparedSoftScreen::new(
            &screen,
            0.1,
            OpenOptions {
                maximum_unknowns: 0,
                ..options()
            },
            || false
        )
        .is_err()
    );
    assert!(
        prepared
            .prepare_frequency(0.1, medium(), 1, || false)
            .is_err()
    );
    assert!(
        prepared
            .prepare_frequency(0.2, medium(), 100_000, || false)
            .is_err()
    );
    let frequency = prepared
        .prepare_frequency(0.1, medium(), 100_000, || false)
        .unwrap();
    let solution = frequency
        .scatter(plane(C64::ONE), limits(), || false)
        .unwrap();
    assert!(solution.density[0].abs() > 0.0);
    assert!(
        frequency
            .scatter(
                plane(C64::ONE),
                SoftScreenScatteringLimits {
                    maximum_density_entries: 0,
                    ..limits()
                },
                || false
            )
            .is_err()
    );
    assert!(
        frequency
            .scatter(
                plane(C64::ONE),
                SoftScreenScatteringLimits {
                    maximum_source_point_products: 1,
                    ..limits()
                },
                || false
            )
            .is_err()
    );
    let source = [Monopole {
        position_m: [0.2, 0.2, 0.0],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        frequency
            .scatter(IncidentField::Monopoles(&source), limits(), || false)
            .is_err()
    );
    assert!(
        solution
            .prepare_receivers(
                &[[0.2, 0.2, 0.0]],
                1e-4,
                receiver_options(),
                100_000,
                || false
            )
            .is_err()
    );
    assert!(
        solution
            .prepare_receivers(&[[0.2, 0.2, 1.0]], 1e-4, receiver_options(), 1, || false)
            .is_err()
    );
    assert!(matches!(
        frequency.scatter(plane(C64::ONE), limits(), || true),
        Err(HelmholtzError::Cancelled)
    ));
}

#[test]
fn prepared_replay_parallel_rows_and_orientation_reversal_preserve_scattering() {
    let screen = disc(1, 8);
    let reverse = OpenPanels::new(
        screen.vertices.clone(),
        screen
            .triangles
            .iter()
            .map(|t| [t[0], t[2], t[1]])
            .collect(),
    )
    .unwrap();
    let k = 0.25;
    let original = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let parallel = PreparedSoftScreen::new(
        &screen,
        k,
        OpenOptions {
            workers: 3,
            ..options()
        },
        || false,
    )
    .unwrap();
    let reversed = PreparedSoftScreen::new(&reverse, k, options(), || false).unwrap();
    let frequency = original
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let second = parallel
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let reversed = reversed
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let a = frequency
        .scatter(plane(C64::ONE), limits(), || false)
        .unwrap();
    let b = second.scatter(plane(C64::ONE), limits(), || false).unwrap();
    let c = reversed
        .scatter(plane(C64::ONE), limits(), || false)
        .unwrap();
    let replay = frequency
        .scatter_many(&[plane(C64::ONE), plane(C64::ONE)], limits(), || false)
        .unwrap();
    for i in 0..a.density.len() {
        assert_eq!(a.density[i], b.density[i]);
        assert_eq!(a.density[i], replay[0].density[i]);
        assert_eq!(a.density[i], replay[1].density[i]);
        close(a.density[i], c.density[i], 2e-10);
    }
}
