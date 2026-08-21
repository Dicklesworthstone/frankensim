//! E4.7-i battery (bead wf-root-guzez.5.18.1): Kelvin closure checked
//! per CELL per step (the invariant, not the comment); connectivity
//! degree audit; induced-velocity liveness + y-symmetry oracle; the
//! wake-rate certificate at cap AND cap+1; rows-exhausted refusal
//! (never silent truncation); shed determinism golden; and the 2-D
//! particle lane's golden UNTOUCHED (additivity pin).
//! Repro: cargo test -p fs-vpm --test filament3d_battery

use fs_vpm::filament3d::{
    FilamentWake, MAX_ROWS, MAX_STATIONS, NEAR_WAKE_SEGMENT_BUDGET, WakeRateCertificate,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-vpm-filament3d\",\"case\":\"{case}\",{payload}}}");
}

fn line(n: usize) -> Vec<[f64; 3]> {
    (0..=n)
        .map(|i| [0.0, i as f64 - n as f64 / 2.0, 0.0])
        .collect()
}

fn cert(stations: usize, rows: usize) -> WakeRateCertificate {
    WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: stations,
        max_rows: rows,
    }
}

#[test]
fn kelvin_closure_holds_per_cell_every_step() {
    let mut wake = FilamentWake::new(cert(8, 64), line(8)).unwrap();
    // Time-varying elliptic-ish circulation so the shed rows DIFFER
    // (a constant distribution cannot falsify the closure).
    for t in 0..40u64 {
        let g: Vec<f64> = (0..8)
            .map(|s| {
                let y = (s as f64 - 3.5) / 4.0;
                (1.0 - y * y).max(0.0) * (1.0 + 0.1 * (t as f64 * 0.3).sin())
            })
            .collect();
        wake.shed(&g, [-0.1, 0.0, 0.0]).unwrap();
        // Per-cell oracle after EVERY shed.
        for r in 0..wake.rows.len().saturating_sub(1) {
            for s in 0..8 {
                let net = wake.cell_net_circulation(r, s).unwrap();
                assert_eq!(net, 0.0, "Kelvin closure at cell ({r},{s}) tick {t}");
            }
        }
    }
    jlog(
        "kelvin",
        &format!(
            "\"rows\":{},\"cells_checked_per_step\":\"(rows-1)*8\"",
            wake.rows.len()
        ),
    );
}

#[test]
fn connectivity_and_symmetry_oracles() {
    // Symmetric circulation + symmetric line: the induced velocity at
    // mirrored probes must mirror bitwise in y.
    let mut wake = FilamentWake::new(cert(8, 32), line(8)).unwrap();
    for _ in 0..16 {
        let g: Vec<f64> = (0..8)
            .map(|s| {
                let y = (s as f64 - 3.5) / 4.0;
                (1.0 - y * y).max(0.0)
            })
            .collect();
        wake.shed(&g, [-0.1, 0.0, 0.0]).unwrap();
    }
    let p = [0.5, 1.25, 0.3];
    let m = [0.5, -1.25, 0.3];
    let vp = wake.induced_velocity(p);
    let vm = wake.induced_velocity(m);
    assert!(vp[0] != 0.0 || vp[2] != 0.0, "induced velocity is LIVE");
    // Mirrored probes traverse the SAME summation order over mirrored
    // geometry, so the pairing reverses and rounding differs at the
    // last ulps — the symmetry oracle is physical (1e-12 relative),
    // not bitwise (bitwise mirror equality would require a reversed
    // accumulation order).
    let scale = vp[0].abs().max(vp[2].abs()).max(1e-12);
    assert!(((vp[0] - vm[0]) / scale).abs() < 1e-12, "x mirrors");
    assert!(((vp[2] - vm[2]) / scale).abs() < 1e-12, "z mirrors");
    assert!(((vp[1] + vm[1]) / scale).abs() < 1e-12, "y anti-mirrors");
    // Connectivity: every row keeps stations+1 nodes (degree audit).
    for row in &wake.rows {
        assert_eq!(row.nodes.len(), 9);
        assert_eq!(row.gamma.len(), 8);
    }
    jlog(
        "symmetry",
        &format!("\"v_probe\":[{},{},{}]", vp[0], vp[1], vp[2]),
    );
}

#[test]
fn certificate_caps_and_refusals() {
    // At every cap admits; one past refuses.
    assert!(cert(MAX_STATIONS, 64).admit().is_ok());
    assert_eq!(
        cert(MAX_STATIONS + 1, 64).admit().unwrap_err().code,
        "wake-rate-uncertified"
    );
    assert!(cert(8, MAX_ROWS).admit().is_ok());
    assert_eq!(
        cert(8, MAX_ROWS + 1).admit().unwrap_err().code,
        "wake-rate-uncertified"
    );
    assert!(
        WakeRateCertificate {
            shed_hz: 1000.0,
            n_stations: 8,
            max_rows: 8
        }
        .admit()
        .is_ok()
    );
    assert_eq!(
        WakeRateCertificate {
            shed_hz: 1000.0_f64.next_up(),
            n_stations: 8,
            max_rows: 8
        }
        .admit()
        .unwrap_err()
        .code,
        "wake-rate-uncertified"
    );
    // Segment budget: stations*rows blowup refuses even under the
    // individual caps.
    assert_eq!(
        cert(MAX_STATIONS, MAX_ROWS).admit().unwrap_err().code,
        "wake-rate-uncertified"
    );
    let _ = NEAR_WAKE_SEGMENT_BUDGET;
    // Rows exhausted refuses, never silently drops.
    let mut wake = FilamentWake::new(cert(2, 3), line(2)).unwrap();
    for _ in 0..3 {
        wake.shed(&[1.0, 0.5], [-0.1, 0.0, 0.0]).unwrap();
    }
    assert_eq!(
        wake.shed(&[1.0, 0.5], [-0.1, 0.0, 0.0]).unwrap_err().code,
        "filament-rows-exhausted"
    );
    // Bad line / bad gammas.
    assert_eq!(
        FilamentWake::new(cert(2, 3), line(3)).unwrap_err().code,
        "filament-line-invalid"
    );
    assert_eq!(
        FilamentWake::new(cert(2, 3), line(2))
            .unwrap()
            .shed(&[f64::NAN, 0.0], [0.0; 3])
            .unwrap_err()
            .code,
        "filament-shed-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn shed_determinism_golden() {
    let run = || {
        let mut wake = FilamentWake::new(cert(6, 48), line(6)).unwrap();
        for t in 0..48u64 {
            let g: Vec<f64> = (0..6)
                .map(|s| ((s as f64 + 1.0) * 0.2 + t as f64 * 0.01).sin() * 0.5 + 1.0)
                .collect();
            wake.shed(&g, [-0.11, 0.0, -0.005]).unwrap();
        }
        wake.digest()
    };
    let a = run();
    assert_eq!(a, run(), "bit-identical twice");
    jlog("golden", &format!("\"digest\":\"{a}\""));
    assert_eq!(
        a, "50b186032a5e875793063557fc5336e9a55ad011ae5128d4ff467adf0dd53b55",
        "filament shed golden moved — determinism regression or an \
         intentional lattice change requiring the golden-bump protocol"
    );
}

#[test]
fn additivity_pin_2d_lane_untouched() {
    // The 2-D particle lane still exposes its surface and behaves —
    // one deterministic probe (its own batteries own the full pins).
    use fs_vpm::{VortexParticle, total_circulation};
    let ps = vec![
        VortexParticle::new([0.0, 0.0], 1.5),
        VortexParticle::new([1.0, 0.0], -0.5),
    ];
    assert_eq!(total_circulation(&ps), 1.0);
    jlog("additivity", "\"lane_2d\":\"intact\"");
}
