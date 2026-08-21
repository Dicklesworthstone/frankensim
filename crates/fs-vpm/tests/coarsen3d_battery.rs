//! E4.7-ii battery (bead wf-root-guzez.5.18.2): every coarsening
//! invariant EXECUTED per item — Kelvin closure re-checked per cell on
//! the coarsened lattice; per-STATION impulse exact (never a
//! totals-only sum); connectivity; core second-moment growth
//! (parallel-axis, retained); symmetry preserved bitwise-per-pair;
//! the FORBIDDEN naive-decimation falsifier violates the impulse
//! invariant and is CAUGHT; the mixed-norm metric is REPORTED and
//! nonzero (never forced to vanish); caps at cap AND cap+1;
//! determinism golden.
//! Repro: cargo test -p fs-vpm --test coarsen3d_battery

use fs_vpm::coarsen3d::{coarsen_oldest, naive_decimate, station_impulse};
use fs_vpm::filament3d::{FilamentWake, WakeRateCertificate};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-vpm-coarsen3d\",\"case\":\"{case}\",{payload}}}");
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
                let y = (s as f64 - 3.5) / 4.0;
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
    // Kelvin closure still holds per cell on the COARSENED lattice.
    for r in 0..wake.rows.len().saturating_sub(1) {
        for s in 0..8 {
            assert_eq!(
                wake.cell_net_circulation(r, s).unwrap(),
                0.0,
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
            assert_eq!(row.core2_m2, 0.0, "fresh row {i} unspread");
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
                let y = (s as f64 - 3.5) / 4.0;
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
    assert_eq!(
        a, "8f1e9a12690e2cd76b4062c2f7979f06616ed122d801b8cc2811f8d5207c2022",
        "coarsen golden moved — determinism regression or an \
         intentional merge-rule change requiring the golden-bump protocol"
    );
}
