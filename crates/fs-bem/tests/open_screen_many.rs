//! PreparedScreen::solve_many (.372): one assembly and LU per frequency for
//! many drives, bitwise N separate solve calls at every declared worker count,
//! with each drive's own checks and whole-call refusal.
use fs_bem::helmholtz::bc::{Admittance, SurfaceBc};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::reciprocal_radiation::open_screen::{
    OpenOptions, OpenPanels, PreparedScreen, ScreenPortLimits,
};
use fs_math::{c64::C64, det};
use std::f64::consts::{PI, TAU};
use std::time::Instant;

fn options(workers: usize) -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 4,
        triangle_order: 4,
        angular_refinement: 2,
        maximum_unknowns: 8192,
        maximum_cache_bytes: 1_000_000_000,
        maximum_work: u64::MAX,
        workers,
    }
}

/// The open-disc fixture of tests/open_and_impedance.rs: r = sin(pi t / 2)
/// rings, quadratic grading toward the free edge.
fn disc(radial: usize, azimuth: usize) -> OpenPanels {
    let mut vertices = vec![[0.0; 3]];
    for j in 1..=radial {
        let r = det::sin(PI * j as f64 / (2.0 * radial as f64));
        for m in 0..azimuth {
            let phi = TAU * m as f64 / azimuth as f64;
            vertices.push([r * det::cos(phi), r * det::sin(phi), 0.0]);
        }
    }
    let at = |j: usize, m: usize| 1 + (j - 1) * azimuth + m % azimuth;
    let mut triangles = Vec::new();
    for m in 0..azimuth {
        triangles.push([0, at(1, m), at(1, m + 1)]);
    }
    for j in 1..radial {
        for m in 0..azimuth {
            triangles.push([at(j, m), at(j + 1, m), at(j + 1, m + 1)]);
            triangles.push([at(j, m), at(j + 1, m + 1), at(j, m + 1)]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}

fn bits(a: &[C64], b: &[C64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert_eq!(
            (x.re.to_bits(), x.im.to_bits()),
            (y.re.to_bits(), y.im.to_bits())
        );
    }
}

/// Four distinct drives: piston, a dipole-like tilt, an m = 2 pattern and a
/// complex mix, so the back-substitutions differ.
fn drives(surface: &OpenPanels) -> Vec<Vec<C64>> {
    let v = surface.vertices();
    vec![
        v.iter().map(|_| C64::new(1.0, 0.0)).collect(),
        v.iter().map(|x| C64::new(x[0], 0.0)).collect(),
        v.iter()
            .map(|x| C64::new(x[0] * x[0] - x[1] * x[1], 0.0))
            .collect(),
        v.iter().map(|x| C64::new(1.0 + x[0], 0.5 * x[1])).collect(),
    ]
}

#[test]
fn solve_many_is_bitwise_separate_solves_at_every_worker_count() {
    let start = Instant::now();
    let surface = disc(5, 24);
    let medium = Medium::air();
    let fields = drives(&surface);
    let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
    let beta = vec![C64::new(0.3, 0.1); surface.vertices().len()];
    let admittances = [
        Admittance::Rigid,
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
    ];
    for workers in [1, 2, 3, 7] {
        let prepared = PreparedScreen::new(&surface, 2.0, options(workers), || false).unwrap();
        for (law, admittance) in admittances.iter().enumerate() {
            let many = prepared
                .solve_many(2.0, medium, *admittance, &refs, 8, || false)
                .unwrap();
            assert_eq!(many.len(), refs.len());
            for (drive, batch) in refs.iter().zip(&many) {
                let single = prepared
                    .solve(
                        2.0,
                        medium,
                        SurfaceBc {
                            drive_velocity: drive,
                            admittance: *admittance,
                        },
                        || false,
                    )
                    .unwrap();
                bits(single.pressure_jump_pa(), batch.pressure_jump_pa());
                for (s, b) in [
                    (single.absorbed_power_w, batch.absorbed_power_w),
                    (single.driving_power_w, batch.driving_power_w),
                    (single.backward_residual, batch.backward_residual),
                    (
                        single.radiated_power_w().unwrap(),
                        batch.radiated_power_w().unwrap(),
                    ),
                ] {
                    assert_eq!(s.to_bits(), b.to_bits(), "workers {workers}, law {law}");
                }
                assert_eq!(single.force_n().re.to_bits(), batch.force_n().re.to_bits());
                assert_eq!(single.force_n().im.to_bits(), batch.force_n().im.to_bits());
            }
        }
    }
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"open_screen_solve_many_bitwise\",\"drives\":{},\"workers\":[1,2,3,7],\"admittances\":[\"rigid\",\"local\"],\"elapsed_s\":{:.3}}}",
        refs.len(),
        start.elapsed().as_secs_f64()
    );
}

#[test]
fn solve_many_refuses_whole_calls() {
    let surface = disc(3, 12);
    let medium = Medium::air();
    let fields = drives(&surface);
    let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
    let prepared = PreparedScreen::new(&surface, 2.0, options(1), || false).unwrap();
    let refused = |r: Result<Vec<_>, HelmholtzError>| r.err().expect("must refuse");
    // The declared ceiling, cap and cap + 1.
    assert!(
        prepared
            .solve_many(2.0, medium, Admittance::Rigid, &refs, refs.len(), || false)
            .is_ok()
    );
    assert!(matches!(
        refused(prepared.solve_many(
            2.0,
            medium,
            Admittance::Rigid,
            &refs,
            refs.len() - 1,
            || false
        )),
        HelmholtzError::BadParameter { .. }
    ));
    assert!(matches!(
        refused(prepared.solve_many(2.0, medium, Admittance::Rigid, &[], 8, || false)),
        HelmholtzError::BadParameter { .. }
    ));
    // A frequency outside the prepared band.
    assert!(matches!(
        refused(prepared.solve_many(2.5, medium, Admittance::Rigid, &refs, 8, || false)),
        HelmholtzError::BadParameter { .. }
    ));
    // One malformed drive refuses the call before any factor is formed.
    let short = vec![C64::new(1.0, 0.0); surface.vertices().len() - 1];
    let mixed: Vec<&[C64]> = vec![refs[0], &short, refs[1]];
    assert!(matches!(
        refused(prepared.solve_many(2.0, medium, Admittance::Rigid, &mixed, 8, || false)),
        HelmholtzError::BadParameter { .. }
    ));
    // Cancellation is polled between back-substitutions too: a gate that trips
    // late in the call still refuses it whole.
    let mut polls = 0usize;
    let total = {
        let mut count = 0usize;
        prepared
            .solve_many(2.0, medium, Admittance::Rigid, &refs, 8, || {
                count += 1;
                false
            })
            .unwrap();
        count
    };
    let late = prepared.solve_many(2.0, medium, Admittance::Rigid, &refs, 8, || {
        polls += 1;
        polls >= total - 1
    });
    assert!(matches!(late.err(), Some(HelmholtzError::Cancelled)));
}

fn port_limits(surface: &OpenPanels, ports: usize) -> ScreenPortLimits {
    ScreenPortLimits {
        maximum_ports: ports,
        maximum_matrix_entries: ports * ports,
        maximum_jump_entries: ports * surface.vertices().len(),
        maximum_work: u64::MAX,
    }
}

fn area_normal(surface: &OpenPanels, facet: usize) -> (f64, [f64; 3]) {
    let t = surface.triangles()[facet].map(|i| surface.vertices()[i]);
    let a: [f64; 3] = std::array::from_fn(|i| t[1][i] - t[0][i]);
    let b: [f64; 3] = std::array::from_fn(|i| t[2][i] - t[0][i]);
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let length = n.iter().map(|v| v * v).sum::<f64>().sqrt();
    (0.5 * length, n.map(|v| v / length))
}

// Independent degree-two triangle quadrature, rather than the production
// triangle mass coefficients or its retained weak-drive dot product.
const WORK_POINTS: [[f64; 3]; 3] = [
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];
fn quadrature_force(surface: &OpenPanels, trace: &[[f64; 3]], jump: &[C64]) -> C64 {
    let mut force = C64::ZERO;
    for (f, &ids) in surface.triangles().iter().enumerate() {
        let (area, _) = area_normal(surface, f);
        for q in WORK_POINTS {
            let velocity = (0..3).map(|a| q[a] * trace[f][a]).sum::<f64>();
            let pressure = (0..3).fold(C64::ZERO, |v, a| v + jump[ids[a]].scale(q[a]));
            force = force + pressure.scale(area / 3.0 * velocity);
        }
    }
    force
}
fn close_fields(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= tolerance * scale.max(1e-300),
        "field error {error}, scale {scale}"
    );
}

#[test]
fn velocity_ports_retain_full_p1_work_and_coherent_two_sided_fields() {
    let surface = disc(3, 12);
    let medium = Medium::air();
    let prepared = PreparedScreen::new(&surface, 0.7, options(1), || false).unwrap();
    let nodal: Vec<Vec<f64>> = vec![
        surface.vertices().iter().map(|_| 1.0).collect(),
        surface
            .vertices()
            .iter()
            .map(|x| 0.3 + x[0] - 0.4 * x[1])
            .collect(),
        vec![0.0; surface.vertices().len()],
    ];
    let rows: Vec<Vec<[f64; 3]>> = nodal
        .iter()
        .map(|row| {
            surface
                .triangles()
                .iter()
                .map(|ids| ids.map(|v| row[v]))
                .collect()
        })
        .collect();
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let source = prepared
        .solve_velocity_ports(0.7, medium, &refs, port_limits(&surface, 3), || false)
        .unwrap();
    assert_eq!(source.ports(), 3);
    assert_eq!(
        source.frequency_work(),
        prepared.velocity_port_work(3).unwrap()
    );
    assert!(source.maximum_backward_residual() < 1e-10);
    let scale = source
        .impedance()
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max);
    for (j, solution) in source.solutions().iter().enumerate() {
        let drive: Vec<_> = nodal[j].iter().copied().map(C64::from_re).collect();
        let original = prepared
            .solve(
                0.7,
                medium,
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: Admittance::Rigid,
                },
                || false,
            )
            .unwrap();
        close_fields(
            solution.pressure_jump_pa(),
            original.pressure_jump_pa(),
            2e-13,
        );
        for (i, row) in rows.iter().enumerate() {
            let expected = quadrature_force(&surface, row, solution.pressure_jump_pa());
            assert!((source.impedance()[3 * i + j] - expected).abs() <= 2e-13 * scale);
            assert!(
                (source.impedance()[3 * i + j] - source.impedance()[3 * j + i]).abs()
                    <= 1e-12 * scale
            );
        }
        for (v, jump) in solution.pressure_jump_pa().iter().enumerate() {
            if surface.is_free_edge_vertex(v) {
                assert_eq!(*jump, C64::ZERO);
            }
        }
    }
    assert!(source.impedance()[1].abs() > 1e-3 * source.impedance()[0].abs());
    assert!(
        source.impedance()[0].im < 0.0,
        "exp(-iwt) added mass is negative imaginary"
    );
    assert!((source.impedance()[0] - source.solutions()[0].force_n()).abs() < 2e-13 * scale);
    let rates = [
        C64::new(0.6, -0.2),
        C64::new(-0.4, 0.7),
        C64::new(3.0, -5.0),
    ];
    let combined: Vec<_> = (0..surface.vertices().len())
        .map(|v| (0..3).fold(C64::ZERO, |sum, j| sum + rates[j].scale(nodal[j][v])))
        .collect();
    let original = prepared
        .solve(
            0.7,
            medium,
            SurfaceBc {
                drive_velocity: &combined,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    close_fields(
        &source.pressure_jump(&rates, || false).unwrap(),
        original.pressure_jump_pa(),
        3e-13,
    );
    let work = source.input_power_w(&rates).unwrap();
    let radiation = source.radiated_power_w(&rates, || false).unwrap();
    assert!(work > 0.0);
    assert!((work - original.driving_power_w).abs() < 2e-12 * work);
    assert!((work - radiation).abs() < 1e-9 * work);
    let points = [[0.2, -0.1, 2.0], [0.2, -0.1, -2.0]];
    let receivers = fs_bem::near_field::Options::default();
    let pressure = source
        .pressure_at(&rates, &points, 1e-5, receivers, || false)
        .unwrap();
    close_fields(
        &pressure,
        &original
            .pressure_at(&points, 1e-5, receivers, || false)
            .unwrap(),
        4e-13,
    );
    assert!((pressure[0] + pressure[1]).abs() <= 1e-12 * pressure[0].abs());
    close_fields(
        &[source.far_field(&rates, [0.0, 0.0, 1.0], || false).unwrap()],
        &[original.far_field([0.0, 0.0, 1.0]).unwrap()],
        4e-13,
    );
}

#[test]
fn folded_facet_normals_keep_the_exact_pressure_test_pairing() {
    let flat = disc(3, 12);
    let surface = OpenPanels::new(
        flat.vertices()
            .iter()
            .map(|p| [p[0], p[1], 0.35 * (p[0] + 0.2 * p[1]).abs()])
            .collect(),
        flat.triangles().to_vec(),
    )
    .unwrap();
    let prepared = PreparedScreen::new(&surface, 0.7, options(1), || false).unwrap();
    let medium = Medium::air();
    let directions = [[0.6, -0.2, 1.0], [-0.4, 0.8, 0.1]];
    let rows: Vec<Vec<[f64; 3]>> = directions
        .iter()
        .map(|direction| {
            (0..surface.triangles().len())
                .map(|f| {
                    let (_, normal) = area_normal(&surface, f);
                    [(0..3).map(|a| direction[a] * normal[a]).sum::<f64>(); 3]
                })
                .collect()
        })
        .collect();
    let range = rows[0]
        .iter()
        .map(|r| r[0])
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        });
    assert!(
        range.1 - range.0 > 0.1,
        "physical nodal translation has different facet normal traces"
    );
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let source = prepared
        .solve_velocity_ports(0.7, medium, &refs, port_limits(&surface, 2), || false)
        .unwrap();
    let active: Vec<_> = (0..surface.vertices().len())
        .filter(|&v| !surface.is_free_edge_vertex(v))
        .collect();
    let mut dof = vec![None; surface.vertices().len()];
    for (i, &v) in active.iter().enumerate() {
        dof[v] = Some(i);
    }
    let n = active.len();
    let mut mass = vec![C64::ZERO; n * n];
    let mut loads = vec![vec![C64::ZERO; n]; 2];
    for (f, &ids) in surface.triangles().iter().enumerate() {
        let (area, _) = area_normal(&surface, f);
        for q in WORK_POINTS {
            for a in 0..3 {
                let Some(i) = dof[ids[a]] else {
                    continue;
                };
                for b in 0..3 {
                    if let Some(j) = dof[ids[b]] {
                        mass[i * n + j] = mass[i * n + j] + C64::from_re(area / 3.0 * q[a] * q[b]);
                    }
                }
                for r in 0..2 {
                    let velocity = (0..3).map(|b| q[b] * rows[r][f][b]).sum::<f64>();
                    loads[r][i] = loads[r][i] + C64::from_re(area / 3.0 * q[a] * velocity);
                }
            }
        }
    }
    // A continuous nodal drive with the same weak load must have the SAME
    // original BEM response, even though the facet-local fields are different.
    let mass_factor = fs_la::eigen_complex::lu_complex(&mass, n).unwrap();
    for (j, load) in loads.iter_mut().enumerate() {
        mass_factor.solve(load);
        let mut drive = vec![C64::ZERO; surface.vertices().len()];
        for (i, &v) in active.iter().enumerate() {
            drive[v] = load[i];
        }
        let original = prepared
            .solve(
                0.7,
                medium,
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: Admittance::Rigid,
                },
                || false,
            )
            .unwrap();
        close_fields(
            source.solutions()[j].pressure_jump_pa(),
            original.pressure_jump_pa(),
            2e-12,
        );
        let scale = source
            .impedance()
            .iter()
            .map(|z| z.abs())
            .fold(0.0, f64::max);
        for i in 0..2 {
            let expected = quadrature_force(&surface, &rows[i], original.pressure_jump_pa());
            assert!((source.impedance()[2 * i + j] - expected).abs() <= 3e-12 * scale);
        }
    }
}

#[test]
fn velocity_port_limits_and_late_cancellation_refuse_the_whole_batch() {
    let surface = disc(2, 8);
    let prepared = PreparedScreen::new(&surface, 0.7, options(1), || false).unwrap();
    let rows = [
        vec![[1.0; 3]; surface.triangles().len()],
        vec![[-0.4; 3]; surface.triangles().len()],
    ];
    let refs = [rows[0].as_slice(), rows[1].as_slice()];
    let limits = port_limits(&surface, 2);
    for missing in 0..4 {
        let mut short = limits;
        match missing {
            0 => short.maximum_ports -= 1,
            1 => short.maximum_matrix_entries -= 1,
            2 => short.maximum_jump_entries -= 1,
            _ => short.maximum_work = prepared.velocity_port_work(2).unwrap() - 1,
        }
        assert!(
            prepared
                .solve_velocity_ports(0.7, Medium::air(), &refs, short, || false)
                .is_err()
        );
    }
    let malformed = [refs[0], &refs[1][1..]];
    assert!(
        prepared
            .solve_velocity_ports(0.7, Medium::air(), &malformed, limits, || false)
            .is_err()
    );
    let mut total = 0;
    let reference = prepared
        .solve_velocity_ports(0.7, Medium::air(), &refs, limits, || {
            total += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    let result = prepared.solve_velocity_ports(0.7, Medium::air(), &refs, limits, || {
        count += 1;
        count >= total - 1
    });
    assert!(matches!(result.err(), Some(HelmholtzError::Cancelled)));
    let replay = prepared
        .solve_velocity_ports(0.7, Medium::air(), &refs, limits, || false)
        .unwrap();
    bits(reference.impedance(), replay.impedance());
    let rates = [C64::ONE, C64::ONE];
    assert!(matches!(
        reference.pressure_jump(&rates, || true),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(reference.pressure_jump(&rates[..1], || false).is_err());
}
