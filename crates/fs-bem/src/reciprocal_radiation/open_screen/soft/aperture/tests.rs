//! G1/G3: complementary fields, finite sources, flux, geometry and refinement.
use super::*;
use crate::helmholtz::bc::{ForceDipole, Monopole, StressQuadrupole};

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
        maximum_source_point_products: 1_000_000,
        minimum_source_clearance_m: 1e-5,
    }
}
fn receiver_options() -> near_field::Options {
    near_field::Options {
        relative_tolerance: 1e-7,
        maximum_depth: 12,
        maximum_kernel_evaluations: 5_000_000,
    }
}
fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.6, 0.0, 0.8],
        amplitude_pa: C64::new(0.8, -0.6),
    })
}
fn mesh(n: usize, centers: &[f64]) -> OpenPanels {
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for &center in centers {
        let start = vertices.len();
        for j in 0..=n {
            for i in 0..=n {
                vertices.push([
                    center - 0.5 + i as f64 / n as f64,
                    -0.5 + j as f64 / n as f64,
                    0.0,
                ]);
            }
        }
        let at = |i, j| start + j * (n + 1) + i;
        for j in 0..n {
            for i in 0..n {
                triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
                triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
            }
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}
fn close(actual: C64, expected: C64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "{actual:?} versus {expected:?}, error {:e}",
        (actual - expected).abs()
    );
}

#[test]
fn babinet_fields_and_transmitted_power_reuse_the_coupled_density() {
    let screen = mesh(2, &[-0.8, 0.8]);
    let k = 0.7;
    let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 64 * 1024 * 1024, || false)
        .unwrap();
    let source = frequency.scatter(plane(), limits(), || false).unwrap();
    let aperture = source.as_rigid_aperture([0.0, 0.0, 1.0], || false).unwrap();
    let points = [
        [0.2, 0.1, 0.4],
        [0.2, 0.1, -0.4],
        [-0.7, 0.2, 0.6],
        [-0.7, 0.2, -0.6],
    ];
    let rows = source
        .prepare_receivers(&points, 1e-5, receiver_options(), 1_000_000, || false)
        .unwrap();
    let scattered = rows.evaluate_scattered(|| false).unwrap();
    let receiver = aperture
        .prepare_receivers(&points, 1e-5, receiver_options(), 1_000_000, || false)
        .unwrap();
    let field = receiver.evaluate(|| false).unwrap();
    assert_eq!(receiver.points_m(), &points);
    assert_eq!(receiver.kernel_evaluations(), rows.kernel_evaluations());
    for (i, &point) in points.iter().enumerate() {
        let mut p = -scattered.scalar.pressure[i];
        let mut u = scattered.particle_velocity_m_s[i].map(|v| -v);
        if point[2] < 0.0 {
            let direct = plane().pressure_velocity_at(k, medium(), point).unwrap();
            let reflected = plane()
                .pressure_velocity_at(k, medium(), [point[0], point[1], -point[2]])
                .unwrap();
            p = scattered.scalar.pressure[i] + direct.pressure_pa + reflected.pressure_pa;
            u = core::array::from_fn(|a| {
                scattered.particle_velocity_m_s[i][a]
                    + direct.velocity_m_s[a]
                    + reflected.velocity_m_s[a].scale(if a == 2 { -1.0 } else { 1.0 })
            });
        }
        close(field.scalar.pressure[i], p, 3e-15);
        for a in 0..3 {
            close(field.particle_velocity_m_s[i][a], u[a], 3e-17);
        }
        assert_eq!(
            field.scalar.quadrature_error_estimate_pa[i],
            scattered.scalar.quadrature_error_estimate_pa[i]
        );
    }
    let mut power = 0.0;
    for (d, w) in directions(2) {
        let d = if d[2] > 0.0 { d } else { d.map(|v| -v) };
        let f = aperture.transmitted_far_field(d, || false).unwrap();
        close(f, -source.scattered_far_field(d, || false).unwrap(), 0.0);
        power += w * f.abs().powi(2) / (2.0 * medium().density * medium().sound_speed);
    }
    assert!(power > 0.0);
    assert!((power - aperture.transmitted_power_w()).abs() < 2e-10 * power);
}

#[test]
fn finite_source_mixture_preserves_reflected_force_and_stress_components() {
    let screen = mesh(1, &[-0.8, 0.8]);
    let k = 0.5;
    let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 2_000_000, || false)
        .unwrap();
    let sources = [
        IncidentSource::Monopole(Monopole {
            position_m: [0.1, 0.2, -0.9],
            volume_velocity_m3_s: C64::new(0.003, -0.001),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [-0.3, 0.1, -0.7],
            force_on_fluid_n: [C64::new(0.2, 0.3), C64::from_re(0.4), C64::new(-0.1, 0.2)],
        }),
        IncidentSource::StressQuadrupole(StressQuadrupole {
            position_m: [0.2, -0.1, -0.8],
            stress_moment_n_m: [
                C64::from_re(0.2),
                C64::new(0.1, 0.3),
                C64::from_re(-0.1),
                C64::from_re(0.04),
                C64::new(0.07, -0.03),
                C64::from_re(-0.02),
            ],
        }),
    ];
    let mut images = sources;
    for image in &mut images {
        match image {
            IncidentSource::Monopole(s) => s.position_m[2] = -s.position_m[2],
            IncidentSource::ForceDipole(s) => {
                s.position_m[2] = -s.position_m[2];
                s.force_on_fluid_n[2] = -s.force_on_fluid_n[2];
            }
            IncidentSource::StressQuadrupole(s) => {
                s.position_m[2] = -s.position_m[2];
                s.stress_moment_n_m[4] = -s.stress_moment_n_m[4];
                s.stress_moment_n_m[5] = -s.stress_moment_n_m[5];
            }
            _ => unreachable!(),
        }
    }
    let incident = IncidentField::Sources(&sources);
    let source = frequency.scatter(incident, limits(), || false).unwrap();
    let aperture = source.as_rigid_aperture([0.0, 0.0, 1.0], || false).unwrap();
    let points = [[0.4, 0.3, -0.4], [0.4, 0.3, 0.4]];
    let field = aperture
        .prepare_receivers(&points, 1e-5, receiver_options(), 1_000_000, || false)
        .unwrap()
        .evaluate(|| false)
        .unwrap();
    let scattered = source
        .prepare_receivers(&points, 1e-5, receiver_options(), 1_000_000, || false)
        .unwrap()
        .evaluate_scattered(|| false)
        .unwrap();
    let direct = incident
        .pressure_velocity_at(k, medium(), points[0])
        .unwrap();
    let image = IncidentField::Sources(&images)
        .pressure_velocity_at(k, medium(), points[0])
        .unwrap();
    close(
        field.scalar.pressure[0],
        scattered.scalar.pressure[0] + direct.pressure_pa + image.pressure_pa,
        3e-15,
    );
    for a in 0..3 {
        close(
            field.particle_velocity_m_s[0][a],
            scattered.particle_velocity_m_s[0][a] + direct.velocity_m_s[a] + image.velocity_m_s[a],
            3e-17,
        );
    }
    close(field.scalar.pressure[1], -scattered.scalar.pressure[1], 0.0);
    let mut sum = [ComplexSum::default(); 4];
    for term in &sources {
        let solved = frequency
            .scatter(
                IncidentField::Sources(core::slice::from_ref(term)),
                limits(),
                || false,
            )
            .unwrap();
        let adapter = solved.as_rigid_aperture([0.0, 0.0, 1.0], || false).unwrap();
        let one = adapter
            .prepare_receivers(&points[1..], 1e-5, receiver_options(), 1_000_000, || false)
            .unwrap()
            .evaluate(|| false)
            .unwrap();
        sum[0].add(one.scalar.pressure[0]);
        for a in 0..3 {
            sum[a + 1].add(one.particle_velocity_m_s[0][a]);
        }
    }
    close(field.scalar.pressure[1], sum[0].finish(), 4e-14);
    for a in 0..3 {
        close(
            field.particle_velocity_m_s[1][a],
            sum[a + 1].finish(),
            4e-16,
        );
    }
}

#[test]
fn rational_rotated_plane_retains_translation_and_vector_reflection() {
    // Exactly authored x+z=7; this is not an epsilon-flat rotated quad.
    let screen = OpenPanels::new(
        vec![
            [2.0, 3.0, 5.0],
            [3.0, 3.0, 4.0],
            [3.0, 4.0, 4.0],
            [2.0, 4.0, 5.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    let k = 0.5;
    let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let source = frequency.scatter(plane(), limits(), || false).unwrap();
    let n = core::f64::consts::FRAC_1_SQRT_2;
    let aperture = source.as_rigid_aperture([n, 0.0, n], || false).unwrap();
    let point = [2.5, 3.5, 4.0];
    let field = aperture
        .prepare_receivers(&[point], 1e-5, receiver_options(), 100_000, || false)
        .unwrap()
        .evaluate(|| false)
        .unwrap();
    let scattered = source
        .prepare_receivers(&[point], 1e-5, receiver_options(), 100_000, || false)
        .unwrap()
        .evaluate_scattered(|| false)
        .unwrap();
    let direct = plane().pressure_velocity_at(k, medium(), point).unwrap();
    let reflected = plane()
        .pressure_velocity_at(k, medium(), [7.0 - point[2], point[1], 7.0 - point[0]])
        .unwrap();
    close(
        field.scalar.pressure[0],
        scattered.scalar.pressure[0] + direct.pressure_pa + reflected.pressure_pa,
        4e-15,
    );
    let reflected_velocity = [
        -reflected.velocity_m_s[2],
        reflected.velocity_m_s[1],
        -reflected.velocity_m_s[0],
    ];
    for a in 0..3 {
        close(
            field.particle_velocity_m_s[0][a],
            scattered.particle_velocity_m_s[0][a] + direct.velocity_m_s[a] + reflected_velocity[a],
            4e-17,
        );
    }
    assert!(source.as_rigid_aperture([-n, 0.0, -n], || false).is_err());
    assert!(source.as_rigid_aperture([0.0, 0.0, 1.0], || false).is_err());
}

#[test]
fn aperture_pressure_continuity_refines_and_rigid_complement_has_zero_normal_limit() {
    let k = 0.2;
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let screen = mesh(n, &[0.0]);
        let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
        let frequency = prepared
            .prepare_frequency(k, medium(), 64 * 1024 * 1024, || false)
            .unwrap();
        let source = frequency.scatter(plane(), limits(), || false).unwrap();
        let aperture = source.as_rigid_aperture([0.0, 0.0, 1.0], || false).unwrap();
        let points = [
            [0.13, 0.09, 0.025],
            [0.13, 0.09, -0.025],
            [2.0, 0.1, 0.001],
            [2.0, 0.1, -0.001],
        ];
        let field = aperture
            .prepare_receivers(&points, 1e-5, receiver_options(), 4_000_000, || false)
            .unwrap()
            .evaluate(|| false)
            .unwrap();
        errors.push((field.scalar.pressure[0] - field.scalar.pressure[1]).abs());
        for i in [2, 3] {
            assert!(field.particle_velocity_m_s[i][2].abs() < 3e-6);
        }
        assert!(
            (field.particle_velocity_m_s[0][2] - field.particle_velocity_m_s[1][2]).abs() < 3e-5
        );
    }
    eprintln!("rigid-aperture pressure continuity errors {errors:?}");
    assert!(errors.windows(2).all(|e| e[1] < e[0]));
    assert!(errors[2] < 0.15);
}

#[test]
fn nonplanar_wrong_side_unresolved_and_excess_work_inputs_refuse() {
    // Exact dyadic c=a+3*(b-a), although the float cross product is nonzero.
    // Without rank admission orient3d would call EVERY fourth point coplanar.
    let collinear = [
        [
            1.1492890483454067e-08,
            3.563152819184494e-07,
            0.08935167410790057,
        ],
        [
            1.887444353038774e-08,
            7.532563712378997e-07,
            -0.8506507290805212,
        ],
        [
            3.363754962425508e-08,
            1.5471385498768002e-06,
            -2.7306555354573647,
        ],
    ];
    assert!(
        norm(cross(
            sub(collinear[1], collinear[0]),
            sub(collinear[2], collinear[0])
        )) > 0.0
    );
    assert!(admit_plane(collinear).is_err());
    let triangle = [[0.0; 3], [1e-200, 0.0, 0.0], [0.0, 1e-200, 0.0]];
    admit_plane(triangle).unwrap();
    assert_eq!(
        plane_side(triangle, [1e-201, 1e-201, 1e-201]).unwrap(),
        Sign::Negative
    );
    assert!(
        plane_side(
            [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [0.2, 0.2, 1e-200]
        )
        .is_err()
    );
    let screen = mesh(1, &[0.0]);
    let k = 0.4;
    let prepared = PreparedSoftScreen::new(&screen, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let mut bound = limits();
    bound.maximum_source_point_products = 34;
    let source = frequency.scatter(plane(), bound, || false).unwrap();
    let aperture = source.as_rigid_aperture([0.0, 0.0, 1.0], || false).unwrap();
    for point in [[2.0, 0.0, 0.0], [f64::NAN, 0.0, 0.2]] {
        assert!(
            aperture
                .prepare_receivers(&[point], 1e-5, receiver_options(), 1_000_000, || false)
                .is_err()
        );
    }
    assert!(
        aperture
            .prepare_receivers(
                &[[0.2, 0.1, -0.5]; 18],
                1e-5,
                receiver_options(),
                1_000_000,
                || false
            )
            .is_err()
    );
    assert!(
        aperture
            .prepare_receivers(&[[0.2, 0.1, 0.5]], 1e-5, receiver_options(), 0, || false)
            .is_err()
    );
    assert!(matches!(
        aperture.prepare_receivers(
            &[[0.2, 0.1, 0.5]],
            1e-5,
            receiver_options(),
            1_000_000,
            || true
        ),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(
        aperture
            .transmitted_far_field([0.0, 0.0, -1.0], || false)
            .is_err()
    );
    for position in [[0.0, 0.0, 1.0], [2.0, 0.0, 0.0]] {
        let input = [Monopole {
            position_m: position,
            volume_velocity_m3_s: C64::ONE,
        }];
        let source = frequency
            .scatter(IncidentField::Monopoles(&input), limits(), || false)
            .unwrap();
        assert!(source.as_rigid_aperture([0.0, 0.0, 1.0], || false).is_err());
    }
    let bent = OpenPanels::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1e-12],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    let prepared = PreparedSoftScreen::new(&bent, k, options(), || false).unwrap();
    let frequency = prepared
        .prepare_frequency(k, medium(), 1_000_000, || false)
        .unwrap();
    let source = frequency.scatter(plane(), limits(), || false).unwrap();
    assert!(source.as_rigid_aperture([0.0, 0.0, 1.0], || false).is_err());
}
