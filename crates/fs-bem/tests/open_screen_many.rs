//! PreparedScreen::solve_many (.372): one assembly and LU per frequency for
//! many drives, bitwise N separate solve calls at every declared worker count,
//! with each drive's own checks and whole-call refusal.
use fs_bem::helmholtz::bc::{Admittance, SurfaceBc};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::reciprocal_radiation::open_screen::{OpenOptions, OpenPanels, PreparedScreen};
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
