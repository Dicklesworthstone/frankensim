//! E4.3b2-i battery (bead wf-root-guzez.5.8.1): the A1 FOM state-space
//! at frozen grid points. EXECUTED: dimension oracles; simulated
//! stability (bounded + settled step response) at EVERY frozen-grid
//! point; DC gain vs an INDEPENDENT steady recompute per output;
//! ground-vs-free DC discriminator (images are in the operator, so
//! they must be in the LTI); caps at cap AND cap+1; determinism golden.
//! Repro: cargo test -p fs-wing --test rom_battery --release

use fs_wing::images::CertifiedGround;
use fs_wing::prescribedwake::frozen_grid_v1;
use fs_wing::rom::{A1Lti, assemble_a1_lti, wright_a1_layout_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-a1lti\",\"case\":\"{case}\",{payload}}}");
}

fn ground() -> CertifiedGround {
    CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.000606,
        certificate_rms_m: 0.801,
    }
}

const V: f64 = 13.0;
const ROWS: usize = 120;

#[test]
fn dimensions_and_determinism_golden() {
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let lti = assemble_a1_lti(&layout, &grid.points[0], &ground(), V, ROWS).unwrap();
    assert_eq!(lti.n_stations, 12);
    assert_eq!(lti.order, 24);
    assert_eq!(lti.a.len(), 24 * 24);
    assert_eq!(lti.b.len(), 24 * 2);
    assert_eq!(lti.c.len(), 3 * 24);
    let again = assemble_a1_lti(&layout, &grid.points[0], &ground(), V, ROWS).unwrap();
    assert_eq!(lti.digest, again.digest, "bit-identical twice");
    jlog("golden", &format!("\"digest\":\"{}\"", lti.digest));
    assert_eq!(
        lti.digest, "4a4d55b8d585f60e46900e2fb5e2c7546f880aa08f5339e849f69aa85543ca81",
        "A1 LTI golden moved — determinism regression or an intentional \
         assembly change requiring the golden-bump protocol"
    );
}

#[test]
fn stability_simulated_at_every_frozen_grid_point() {
    // The reduction (5.8.2) presumes a stable FOM: at EVERY registered
    // operating point the step response must stay bounded and settle.
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let mut checked = 0usize;
    for point in &grid.points {
        let lti = assemble_a1_lti(&layout, point, &ground(), V, ROWS).unwrap();
        let y = lti.simulate(&|_| [0.1, 0.02], 1.0 / 240.0, 4_000).unwrap();
        let last = y[y.len() - 1];
        for (o, v) in last.iter().enumerate() {
            assert!(
                v.is_finite() && v.abs() < 1.0e6,
                "point {point:?} out {o}: {v}"
            );
        }
        // Settled: final 5% window spread under 1% of the final value.
        let tail = &y[y.len() - 200..];
        for o in 0..3 {
            let (lo, hi) = tail
                .iter()
                .fold((f64::MAX, f64::MIN), |(l, h), s| (l.min(s[o]), h.max(s[o])));
            let scale = last[o].abs().max(1.0);
            assert!((hi - lo) / scale < 0.01, "unsettled at {point:?} out {o}");
        }
        checked += 1;
    }
    assert_eq!(checked, grid.points.len());
    jlog("stability", &format!("\"points\":{checked}"));
}

#[test]
fn dc_gain_matches_the_independent_steady_recompute() {
    // Per-output oracle at a spread of grid points: the LTI's settled
    // step output must equal the DIRECT steady solve (which never
    // touches A/B/C/D).
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let u = [0.1, 0.02];
    for idx in [0usize, 17, 55, 100, 143] {
        let point = &grid.points[idx];
        let lti = assemble_a1_lti(&layout, point, &ground(), V, ROWS).unwrap();
        let y = lti.simulate(&|_| u, 1.0 / 240.0, 6_000).unwrap();
        let settled = y[y.len() - 1];
        let direct = A1Lti::dc_direct(&layout, point, &ground(), V, ROWS, u).unwrap();
        for o in 0..3 {
            let scale = direct[o].abs().max(1.0);
            assert!(
                (settled[o] - direct[o]).abs() / scale < 5.0e-3,
                "point {idx} out {o}: lti {} vs direct {}",
                settled[o],
                direct[o]
            );
        }
    }
    jlog("dc", "\"points\":5,\"outputs\":3");
}

#[test]
fn ground_point_dc_differs_from_free_air() {
    // Images live in the operator; the LTI inherits them — a ground
    // point (h/b 0.1) and the free-air twin (h/b 10) must give
    // different wing-lift DC (the discriminator that the image path
    // reached the state-space).
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let low = grid
        .points
        .iter()
        .find(|p| p.h_over_b == 0.1 && p.pitch_rad == 0.05 && p.roll_rad == 0.0)
        .unwrap();
    let high = grid
        .points
        .iter()
        .find(|p| {
            p.h_over_b == 10.0
                && p.pitch_rad == 0.05
                && p.roll_rad == 0.0
                && p.canard_rad == low.canard_rad
                && p.warp_rad == low.warp_rad
                && p.convection == low.convection
        })
        .unwrap();
    let u = [0.0, 0.05];
    let dl = A1Lti::dc_direct(&wright_a1_layout_v1(), low, &ground(), V, ROWS, u).unwrap();
    let dh = A1Lti::dc_direct(&layout, high, &ground(), V, ROWS, u).unwrap();
    let rel = (dl[0] - dh[0]).abs() / dh[0].abs().max(1.0);
    assert!(
        rel > 0.01,
        "ground images must move the DC: {dl:?} vs {dh:?}"
    );
    jlog(
        "images",
        &format!("\"low_n\":{},\"high_n\":{}", dl[0], dh[0]),
    );
}

#[test]
fn caps_at_cap_and_cap_plus_one() {
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let p = &grid.points[0];
    assert!(assemble_a1_lti(&layout, p, &ground(), 40.0, ROWS).is_ok());
    for v in [40.0_f64.next_up(), 5.0_f64.next_down(), f64::NAN] {
        assert!(matches!(
            assemble_a1_lti(&layout, p, &ground(), v, ROWS),
            Err(e) if e.code == "a1-speed-invalid"
        ));
    }
    // Wake-row caps pass through from the operator (cap AND cap+1).
    assert!(assemble_a1_lti(&layout, p, &ground(), V, 512).is_ok());
    assert!(matches!(
        assemble_a1_lti(&layout, p, &ground(), V, 513),
        Err(e) if e.code == "wake-rows-invalid"
    ));
    // Simulation domain caps.
    let lti = assemble_a1_lti(&layout, p, &ground(), V, ROWS).unwrap();
    assert!(lti.simulate(&|_| [0.0, 0.0], 0.01, 1).is_ok());
    assert!(matches!(
        lti.simulate(&|_| [0.0, 0.0], 0.01_f64.next_up(), 1),
        Err(e) if e.code == "a1-sim-invalid"
    ));
    assert!(matches!(
        lti.simulate(&|_| [0.0, 0.0], 0.001, 200_001),
        Err(e) if e.code == "a1-sim-invalid"
    ));
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}
