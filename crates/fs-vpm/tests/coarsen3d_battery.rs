//! E4.7-ii battery (bead wf-root-guzez.5.18.2): every coarsening
//! invariant EXECUTED per item — Kelvin closure re-checked per cell on
//! the coarsened lattice; per-STATION vector impulse to roundoff (never a
//! totals-only sum); connectivity; geometric core smoothing;
//! symmetry preserved bitwise-per-pair;
//! the FORBIDDEN naive-decimation falsifier violates the impulse
//! invariant and is CAUGHT; the mixed-norm metric is REPORTED and
//! nonzero (never forced to vanish); caps at cap AND cap+1;
//! determinism golden.
//! Repro: cargo test -p fs-vpm --test coarsen3d_battery

use fs_vpm::coarsen3d::{coarsen_oldest, naive_decimate, station_impulse};
use fs_vpm::filament3d::{FilamentWake, ShedRow, WakeRateCertificate};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-vpm-coarsen3d\",\"case\":\"{case}\",{payload}}}");
}

/// Bitwise f64 equality; both sides are bit-identical by construction.
fn bits_equal(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

fn line(n: usize) -> Vec<[f64; 3]> {
    (0..=n)
        .map(|i| [0.0, i as f64 - n as f64 / 2.0, 0.0])
        .collect()
}

fn built_wake(rows: usize) -> FilamentWake {
    let cert = WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: 8,
        max_rows: rows + 4,
    };
    let mut wake = FilamentWake::new(cert, line(8)).unwrap();
    for t in 0..rows as u64 {
        let g: Vec<f64> = (0..8)
            .map(|s| {
                let y = (f64::from(s) - 3.5) / 4.0;
                (1.0 - y * y).max(0.0) * (1.0 + 0.15 * (t as f64 * 0.37).sin())
            })
            .collect();
        wake.shed(&g, [-0.11, 0.0, -0.004]).unwrap();
    }
    wake
}

#[test]
fn per_station_impulse_exact_and_kelvin_survives() {
    let mut wake = built_wake(32);
    let before: Vec<f64> = (0..8).map(|s| station_impulse(&wake, s)).collect();
    let metric = coarsen_oldest(&mut wake, 8).unwrap();
    assert_eq!(metric.rows_before, 32);
    assert_eq!(metric.rows_after, 24);
    // Per-STATION impulse oracle (a totals-only sum is blind to
    // permutation — workspace law).
    for s in 0..8 {
        let after = station_impulse(&wake, s);
        let scale = before[s].abs().max(1e-12);
        assert!(
            ((after - before[s]) / scale).abs() < 1e-12,
            "station {s}: {} -> {after}",
            before[s]
        );
    }
    // Signed incident-edge balances still hold on the COARSENED lattice.
    for r in 0..=wake.rows.len() {
        for s in 0..=8 {
            assert!(
                wake.node_circulation_residual(r, s).unwrap().abs() < 2e-15,
                "closure at ({r},{s}) after coarsening"
            );
        }
    }
    // Connectivity + core-moment growth on merged rows.
    for (i, row) in wake.rows.iter().enumerate() {
        assert_eq!(row.nodes.len(), 9);
        if i < 8 {
            assert!(row.core2_m2 > 0.0, "merged row {i} retains spread");
        } else {
            assert!(bits_equal(row.core2_m2, 0.0), "fresh row {i} unspread");
        }
    }
    // The mixed-norm metric is REPORTED and honest (nonzero).
    assert!(metric.near_rms_mps > 0.0);
    assert!(metric.far_delta_mps >= 0.0);
    jlog(
        "invariants",
        &format!(
            "\"near_rms\":{},\"far_delta\":{}",
            metric.near_rms_mps, metric.far_delta_mps
        ),
    );
}

#[test]
fn symmetry_is_preserved_pairwise() {
    let mut wake = built_wake(16);
    coarsen_oldest(&mut wake, 4).unwrap();
    for row in &wake.rows {
        for s in 0..4 {
            assert_eq!(
                row.gamma[s].to_bits(),
                row.gamma[7 - s].to_bits(),
                "gamma mirror"
            );
            assert_eq!(
                row.nodes[s][1].to_bits(),
                (-row.nodes[9 - 1 - s][1]).to_bits(),
                "node y mirror"
            );
        }
    }
    jlog("symmetry", "\"pairwise_bitwise\":true");
}

/// Pulsed shedding (strong tick-to-tick alternation — the 120 Hz
/// content naive decimation ALIASES away while lawful pair-merging
/// integrates it exactly).
fn pulsed_wake(rows: usize) -> FilamentWake {
    let cert = WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: 8,
        max_rows: rows + 4,
    };
    let mut wake = FilamentWake::new(cert, line(8)).unwrap();
    for t in 0..rows as u64 {
        let pulse = if t % 2 == 0 { 2.0 } else { 0.2 };
        let g: Vec<f64> = (0..8)
            .map(|s| {
                let y = (f64::from(s) - 3.5) / 4.0;
                (1.0 - y * y).max(0.0) * pulse
            })
            .collect();
        wake.shed(&g, [-0.11, 0.0, -0.004]).unwrap();
    }
    wake
}

#[test]
fn naive_decimation_falsifier_violates_impulse_and_is_caught() {
    let wake0 = pulsed_wake(32);
    let before: Vec<f64> = (0..8).map(|s| station_impulse(&wake0, s)).collect();
    // The FORBIDDEN scheme.
    let mut naive = wake0.clone();
    naive_decimate(&mut naive);
    let mut worst = 0.0f64;
    for s in 0..8 {
        let after = station_impulse(&naive, s);
        worst = worst.max(((after - before[s]) / before[s].abs().max(1e-12)).abs());
    }
    assert!(
        worst > 0.2,
        "naive decimation must visibly violate impulse: {worst}"
    );
    // And the LAWFUL coarsener on the same wake does not.
    let mut lawful = wake0.clone();
    coarsen_oldest(&mut lawful, 16).unwrap();
    let mut worst_lawful = 0.0f64;
    for s in 0..8 {
        let after = station_impulse(&lawful, s);
        worst_lawful = worst_lawful.max(((after - before[s]) / before[s].abs().max(1e-12)).abs());
    }
    assert!(worst_lawful < 1e-12, "lawful coarsen exact: {worst_lawful}");
    jlog(
        "falsifier",
        &format!("\"naive_worst\":{worst},\"lawful_worst\":{worst_lawful}"),
    );
}

#[test]
fn caps_and_determinism() {
    // Exactly 2k rows admits; 2k-1 refuses (cap AND cap+1 class).
    let mut wake = built_wake(8);
    assert!(coarsen_oldest(&mut wake, 4).is_ok());
    let mut wake = built_wake(7);
    assert_eq!(
        coarsen_oldest(&mut wake, 4).unwrap_err().code,
        "coarsen-invalid"
    );
    let mut wake = built_wake(8);
    assert_eq!(
        coarsen_oldest(&mut wake, 0).unwrap_err().code,
        "coarsen-invalid"
    );
    // Determinism golden over the coarsened digest.
    let run = || {
        let mut w = built_wake(32);
        coarsen_oldest(&mut w, 8).unwrap();
        coarsen_oldest(&mut w, 4).unwrap();
        w.digest()
    };
    let a = run();
    assert_eq!(a, run(), "bit-identical twice");
    jlog("golden", &format!("\"digest\":\"{a}\""));
    // Separate area/circulation normalization changes final rounding in this pin.
    // Independent impulse/frame/refusal oracles validate the merge physics.
    assert_eq!(
        a, "26d8bfdb64df343cc9e3a7a4940dc528471c96dae8e1a49d70f5154f1ba15f82",
        "coarsen golden moved — determinism regression or an \
         intentional merge-rule change requiring the golden-bump protocol"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Independent polygon edge integral, without the production vector-area helper.
fn station_vector_impulse(wake: &FilamentWake, station: usize) -> [f64; 3] {
    let origin = wake.line_nodes[station];
    let mut impulse = [0.0; 3];
    for (r, row) in wake.rows.iter().enumerate() {
        let next = wake.rows.get(r + 1).map_or(&wake.line_nodes, |r| &r.nodes);
        let vertices = [
            next[station],
            next[station + 1],
            row.nodes[station + 1],
            row.nodes[station],
        ];
        for edge in 0..4 {
            let a = std::array::from_fn(|c| vertices[edge][c] - origin[c]);
            let b = std::array::from_fn(|c| vertices[(edge + 1) % 4][c] - origin[c]);
            let moment = cross(a, b);
            for c in 0..3 {
                impulse[c] += 0.5 * row.gamma[station] * moment[c];
            }
        }
    }
    impulse
}

fn close_vector(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        assert!((a[c] - b[c]).abs() <= tolerance, "{a:?} != {b:?}");
    }
}

fn map_geometry(wake: &mut FilamentWake, map: impl Fn([f64; 3]) -> [f64; 3]) {
    for node in wake
        .line_nodes
        .iter_mut()
        .chain(wake.rows.iter_mut().flat_map(|r| &mut r.nodes))
    {
        *node = map(*node);
    }
}

#[test]
fn full_vector_impulse_and_core_smoothing_are_rotation_covariant() {
    let original = built_wake(16);
    let mut reference = original.clone();
    coarsen_oldest(&mut reference, 4).unwrap();
    for frame in 0..2 {
        let rotate = |p: [f64; 3]| {
            if frame == 0 {
                [p[2], p[0], p[1]]
            } else {
                let q = [0.6 * p[0] - 0.8 * p[1], 0.8 * p[0] + 0.6 * p[1], p[2]];
                [q[0], 0.8 * q[1] - 0.6 * q[2], 0.6 * q[1] + 0.8 * q[2]]
            }
        };
        let transform = |p| {
            let r = rotate(p);
            [r[0] + 2.0, r[1] - 4.0, r[2] + 1.0]
        };
        let mut moved = original.clone();
        map_geometry(&mut moved, transform);
        let before: Vec<_> = (0..8).map(|s| station_vector_impulse(&moved, s)).collect();
        let total = moved.hydrodynamic_impulse_per_density().unwrap();
        coarsen_oldest(&mut moved, 4).unwrap();
        for (s, &expected) in before.iter().enumerate() {
            close_vector(station_vector_impulse(&moved, s), expected, 5e-13);
        }
        close_vector(
            moved.hydrodynamic_impulse_per_density().unwrap(),
            total,
            5e-13,
        );
        for (actual, expected) in moved.rows.iter().zip(&reference.rows) {
            assert!((actual.core2_m2 - expected.core2_m2).abs() < 1e-13);
            assert!((actual.age_ticks - expected.age_ticks).abs() < 1e-13);
            for (&a, &b) in actual.gamma.iter().zip(&expected.gamma) {
                assert!((a - b).abs() < 1e-12, "circulation {a} != {b}");
            }
            for (&a, &b) in actual.nodes.iter().zip(&expected.nodes) {
                close_vector(a, transform(b), 1e-13);
            }
        }
        let probe = [0.4, 0.3, 1.2];
        close_vector(
            moved.induced_velocity(transform(probe)),
            rotate(reference.induced_velocity(probe)),
            1e-12,
        );
    }
}

#[test]
fn transverse_convection_and_small_length_units_remain_admissible() {
    for length_scale in [1.0, 1.0 / 1_099_511_627_776.0] {
        let mut wake = built_wake(8);
        // A proper cyclic rotation of a flat sheet: every x is exactly zero.
        map_geometry(&mut wake, |p| {
            [0.0, length_scale * p[0], length_scale * p[1]]
        });
        let before = wake.hydrodynamic_impulse_per_density().unwrap();
        coarsen_oldest(&mut wake, 2).unwrap();
        close_vector(
            wake.hydrodynamic_impulse_per_density().unwrap(),
            before,
            1e-13 * length_scale * length_scale,
        );
        assert!(wake.rows[0].core2_m2 > 0.0);
    }
}

#[test]
fn unrepresentable_warped_pair_refuses_without_changing_the_wake() {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 1,
            max_rows: 4,
        },
        vec![[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
    )
    .unwrap();
    wake.rows = vec![
        ShedRow {
            nodes: vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            gamma: vec![1.0],
            core2_m2: 0.01,
            age_ticks: 2.0,
        },
        ShedRow {
            nodes: vec![[1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
            gamma: vec![2.0],
            core2_m2: 0.02,
            age_ticks: 1.0,
        },
    ];
    let before = wake.clone();
    assert_eq!(
        coarsen_oldest(&mut wake, 1).unwrap_err().code,
        "coarsen-impulse-unrepresentable"
    );
    assert_eq!(wake, before);
}

#[test]
fn overflowing_pair_count_and_invalid_geometry_refuse_atomically() {
    let mut wake = built_wake(8);
    let before = wake.clone();
    assert_eq!(
        coarsen_oldest(&mut wake, usize::MAX).unwrap_err().code,
        "coarsen-invalid"
    );
    assert_eq!(wake, before);
    wake.rows[0].nodes.pop();
    let before = wake.clone();
    assert_eq!(
        coarsen_oldest(&mut wake, 2).unwrap_err().code,
        "filament-geometry-invalid"
    );
    assert_eq!(wake, before);
}

fn physical_impulse(wake: &FilamentWake) -> [f64; 3] {
    let mut result = [0.0; 3];
    wake.visit_segments_with_core_growth(0.0, |edge| {
        let a = edge.start;
        let b = edge.end;
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        for i in 0..3 {
            result[i] += 0.5 * edge.circulation * cross[i];
        }
    });
    result
}

fn strip_impulse(wake: &FilamentWake, station: usize) -> [f64; 3] {
    let mut isolated = wake.clone();
    for row in &mut isolated.rows {
        for (s, gamma) in row.gamma.iter_mut().enumerate() {
            if s != station {
                *gamma = 0.0;
            }
        }
    }
    physical_impulse(&isolated)
}

fn vector_near(actual: [f64; 3], expected: [f64; 3]) {
    let scale = actual
        .iter()
        .chain(&expected)
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        (0..3).all(|i| (actual[i] - expected[i]).abs() <= 3e-12 * scale),
        "{actual:?} != {expected:?}"
    );
}

fn sheared_planar_wake(length_scale: f64) -> FilamentWake {
    let nodes = |row: usize| {
        (0..=3)
            .map(|s| {
                let x = row as f64 + 0.07 * (s * s) as f64;
                let y = (1.0 + 0.06 * row as f64) * s as f64 + 0.2 * row as f64;
                [
                    x * length_scale,
                    y * length_scale,
                    (0.3 * x - 0.2 * y) * length_scale,
                ]
            })
            .collect::<Vec<_>>()
    };
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 120.0,
            n_stations: 3,
            max_rows: 12,
        },
        nodes(6),
    )
    .unwrap();
    wake.rows = (0..6)
        .map(|i| ShedRow {
            nodes: nodes(i),
            gamma: vec![
                1.0 + 0.2 * i as f64,
                -(2.0 + 0.1 * i as f64),
                0.5 - 0.03 * i as f64,
            ],
            core2_m2: 0.0,
            age_ticks: (6 - i) as f64,
        })
        .collect();
    wake
}

#[test]
fn g0_g3_nonuniform_width_multiple_pairs_preserve_each_physical_strip() {
    for scale in [1e-6, 1.0, 1e6] {
        let mut wake = sheared_planar_wake(scale);
        let before: Vec<_> = (0..3).map(|s| strip_impulse(&wake, s)).collect();
        coarsen_oldest(&mut wake, 2).unwrap();
        for s in 0..3 {
            vector_near(strip_impulse(&wake, s), before[s]);
        }
        coarsen_oldest(&mut wake, 1).unwrap();
        for s in 0..3 {
            vector_near(strip_impulse(&wake, s), before[s]);
        }
    }
}

#[test]
fn g0_nonrepresentable_vector_impulse_refuses_without_partial_coarsening() {
    for late in [false, true] {
        let count = if late { 4 } else { 2 };
        let mut wake = FilamentWake::new(
            WakeRateCertificate {
                shed_hz: 120.0,
                n_stations: 1,
                max_rows: 8,
            },
            vec![[count as f64, 0.0, 0.0], [count as f64, 1.0, 0.0]],
        )
        .unwrap();
        let heights = if late {
            vec![0.0, 0.0, 1.0, -1.0]
        } else {
            vec![0.0, 1.0]
        };
        let gamma = if late {
            vec![1.0, 0.0, 1.0, 1.0]
        } else {
            vec![1.0, 2.0]
        };
        wake.rows = (0..count)
            .map(|i| ShedRow {
                nodes: vec![[i as f64, 0.0, 0.0], [i as f64, 1.0, heights[i]]],
                gamma: vec![gamma[i]],
                core2_m2: 0.0,
                age_ticks: (count - i) as f64,
            })
            .collect();
        if !late {
            vector_near(physical_impulse(&wake), [0.5, -1.5, 3.0]);
        }
        let before = wake.digest();
        assert_eq!(
            coarsen_oldest(&mut wake, count / 2).unwrap_err().code,
            "coarsen-impulse-unrepresentable"
        );
        assert_eq!(wake.digest(), before);
    }
}

#[test]
fn circulation_survives_overflowing_and_underflowing_dimensional_impulse_products() {
    // The circulation and every geometric coordinate are representable even
    // when Gamma * area is not. The analytical one-span merge has strength
    // (2/3) * (Gamma_a + Gamma_b), independent of either physical scale.
    for (length, strength) in [(1e-60, 1e-240), (1e-60, 1e-200), (1e60, 1e200)] {
        for sign in [1.0, -1.0] {
            let nodes = |x| vec![[x * length, 0.0, 0.0], [x * length, length, 0.0]];
            let mut wake = FilamentWake::new(
                WakeRateCertificate {
                    shed_hz: 1.0,
                    n_stations: 1,
                    max_rows: 2,
                },
                nodes(2.0),
            )
            .unwrap();
            wake.rows = (0..2)
                .map(|r| ShedRow {
                    nodes: nodes(r as f64),
                    gamma: vec![if r == 0 { strength } else { sign * strength }],
                    core2_m2: 0.0,
                    age_ticks: (2 - r) as f64,
                })
                .collect();
            coarsen_oldest(&mut wake, 1).unwrap();
            let expected = (2.0 / 3.0) * (1.0 + sign);
            let actual = wake.rows[0].gamma[0] / strength;
            assert!(
                (actual - expected).abs() < 2e-14,
                "length {length:e}, strength {strength:e}, sign {sign}: {actual} != {expected}"
            );
            assert!((wake.rows[0].core2_m2 / (length * length) - 0.25).abs() < 2e-14);
        }
    }
}
