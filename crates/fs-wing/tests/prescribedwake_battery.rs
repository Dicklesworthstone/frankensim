//! E4.3b1 battery (bead wf-root-guzez.5.7): deterministic assembly over
//! the FROZEN grid with the registered identity, image symmetry at the
//! ground plane (with the image-free liveness twin), the far-height
//! limit (non-vacuous via a low-height twin), the straight-wake
//! cross-kernel oracle vs the production horseshoe legs, caps at cap
//! AND cap+1, golden. Repro: cargo test -p fs-wing --test prescribedwake_battery

use fs_wing::images::CertifiedGround;
use fs_wing::prescribedwake::{
    MAX_WAKE_ROWS, WakeOperatingPoint, assemble_operator, assemble_operator_no_images,
    frozen_grid_v1, probe_velocity,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-e43b1\",\"case\":\"{case}\",{payload}}}");
}

const SPAN_B: f64 = 12.29;

fn ground() -> CertifiedGround {
    CertifiedGround {
        z_m: 0.0,
        certificate_slope: 0.001,
        certificate_rms_m: 0.3,
    }
}

fn stations() -> Vec<[f64; 3]> {
    // Trailing-edge shed points across the span (5 stations).
    (0..5)
        .map(|i| [-1.98, -4.9 + 2.45 * f64::from(i), 0.0])
        .collect()
}

fn probes() -> Vec<[f64; 3]> {
    vec![[0.0, -3.0, 0.0], [0.0, 0.0, 0.0], [0.0, 3.0, -0.5]]
}

fn op_point(h_over_b: f64) -> WakeOperatingPoint {
    WakeOperatingPoint {
        h_over_b,
        pitch_rad: 0.05,
        roll_rad: 0.0,
        canard_rad: 0.0,
        warp_rad: 0.0,
        convection: 1.0,
    }
}

#[test]
fn frozen_grid_identity_is_registered_and_stable() {
    let g1 = frozen_grid_v1();
    let g2 = frozen_grid_v1();
    assert_eq!(g1.digest, g2.digest, "identity stable");
    assert_eq!(g1.points.len(), 4 * 3 * 2 * 3 * 2 * 2, "full product");
    // Canonical nesting: first point is the first entry of every axis.
    let p0 = g1.points[0];
    assert_eq!(p0.h_over_b, 0.1);
    assert_eq!(p0.convection, 0.85);
    jlog(
        "grid",
        &format!(
            "\"n\":{},\"digest\":\"{}\"",
            g1.points.len(),
            &g1.digest[..16]
        ),
    );
}

#[test]
fn assembles_deterministically_over_the_grid() {
    // DONE-WHEN: the operator assembles at EVERY frozen point (all 288),
    // deterministically (bitwise repeat spot-checked on a stride).
    let g = frozen_grid_v1();
    let (st, pr, gr) = (stations(), probes(), ground());
    let mut assembled = 0usize;
    for (i, pt) in g.points.iter().enumerate() {
        let op = assemble_operator(&st, &pr, &gr, pt, 24, 0.5, SPAN_B).unwrap();
        assert_eq!(op.n_probes * op.n_stations, op.w_normal.len());
        assert!(op.w_normal.iter().all(|v| v.is_finite()));
        assert!(op.image_aware);
        if i % 37 == 0 {
            let op2 = assemble_operator(&st, &pr, &gr, pt, 24, 0.5, SPAN_B).unwrap();
            assert_eq!(op.w_normal, op2.w_normal, "bitwise repeat at point {i}");
        }
        assembled += 1;
    }
    assert_eq!(assembled, g.points.len());
    jlog("assembly", &format!("\"points\":{assembled}"));
}

#[test]
fn image_symmetry_zero_normal_flow_at_the_plane() {
    // The DONE-WHEN image-symmetry clause: real + mirror gives ZERO
    // z-velocity ON the ground plane. The image-free twin violates it
    // (liveness — the check would be vacuous otherwise).
    let pt = op_point(0.15);
    let dz = -(0.15 * SPAN_B);
    let shed = [-1.98, -2.0, dz];
    for probe_xy in [[-3.0, 0.5], [2.0, -1.5], [-8.0, 3.0]] {
        let on_plane = [probe_xy[0], probe_xy[1], 0.0];
        let v = probe_velocity(shed, on_plane, 0.0, &pt, 48, 0.5);
        let scale = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-12);
        assert!(
            v[2].abs() < 1e-10 * scale.max(1.0),
            "normal flow at the plane must vanish: {v:?}"
        );
        // Liveness: the raw (un-mirrored) line does NOT vanish there.
        let line_only = {
            // probe_velocity with ground pushed to -infinity emulates
            // no-image: mirror at z = -1e6 contributes ~0.
            probe_velocity(shed, on_plane, -1.0e6, &pt, 48, 0.5)
        };
        assert!(
            line_only[2].abs() > 1e-6,
            "the un-mirrored wake must have normal flow at z=0: {line_only:?}"
        );
    }
    jlog("image-symmetry", "\"zero_normal_flow\":true");
}

#[test]
fn far_height_limit_and_low_height_twin() {
    let (st, pr, gr) = (stations(), probes(), ground());
    let hi = assemble_operator(&st, &pr, &gr, &op_point(10.0), 24, 0.5, SPAN_B).unwrap();
    let hi_free =
        assemble_operator_no_images(&st, &pr, &gr, &op_point(10.0), 24, 0.5, SPAN_B).unwrap();
    let lo = assemble_operator(&st, &pr, &gr, &op_point(0.12), 24, 0.5, SPAN_B).unwrap();
    let lo_free =
        assemble_operator_no_images(&st, &pr, &gr, &op_point(0.12), 24, 0.5, SPAN_B).unwrap();
    let rel = |a: &[f64], b: &[f64]| -> f64 {
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (x, y) in a.iter().zip(b) {
            num += (x - y) * (x - y);
            den += y * y;
        }
        (num / den.max(1e-300)).sqrt()
    };
    let hi_dev = rel(&hi.w_normal, &hi_free.w_normal);
    let lo_dev = rel(&lo.w_normal, &lo_free.w_normal);
    assert!(
        hi_dev < 1e-3,
        "at h/b = 10 the image must be negligible: {hi_dev}"
    );
    assert!(
        lo_dev > 0.05,
        "at h/b = 0.12 the image must matter (non-vacuity): {lo_dev}"
    );
    jlog(
        "far-limit",
        &format!("\"hi_dev\":{hi_dev},\"lo_dev\":{lo_dev}"),
    );
}

#[test]
fn straight_wake_matches_the_production_horseshoe_kernel() {
    // Cross-kernel oracle: one station, no pitch, convection 1, many
    // rows + far cap = a semi-infinite straight trailing line — the
    // SAME object the production horseshoe legs integrate. Compare the
    // z-velocity at probes against fs-wing's trailing-leg kernel via
    // induced_velocity_free on a degenerate one-panel system.
    let pt = WakeOperatingPoint {
        h_over_b: 10.0, // effectively free air
        pitch_rad: 0.0,
        roll_rad: 0.0,
        canard_rad: 0.0,
        warp_rad: 0.0,
        convection: 1.0,
    };
    let dz = -(10.0 * SPAN_B);
    let shed = [0.0, 0.0, dz];
    // Production reference: a horseshoe with endpoints far apart in y so
    // only ONE leg is nearby; unit circulation; stream along +x means
    // legs extend toward −x (the frozen convention).
    let panel = fs_wing::Panel {
        surface: fs_wing::SurfaceId::WingLower,
        a: [0.0, 0.0, dz],
        b: [0.0, 5000.0, dz],
        ctrl: [0.0, 0.5, dz],
        normal: [0.0, 0.0, -1.0],
        width_m: 5000.0,
    };
    // Independent inline Biot-Savart for the BOUND segment (subtracted
    // from the production horseshoe to leave legs only — the oracle's
    // own arithmetic, not the library's).
    let bound_only = |p: [f64; 3], a: [f64; 3], b: [f64; 3]| -> f64 {
        let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
        let n1 = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
        let n2 = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
        let c = [
            r1[1] * r2[2] - r1[2] * r2[1],
            r1[2] * r2[0] - r1[0] * r2[2],
            r1[0] * r2[1] - r1[1] * r2[0],
        ];
        let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
        let r0 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let k = ((r0[0] * r1[0] + r0[1] * r1[1] + r0[2] * r1[2]) / n1
            - (r0[0] * r2[0] + r0[1] * r2[1] + r0[2] * r2[2]) / n2)
            / (4.0 * core::f64::consts::PI * c2);
        c[2] * k
    };
    for probe_off in [[-3.0f64, 1.5, -0.8], [-6.0, -2.0, 1.2], [-1.5, 0.7, 0.4]] {
        let probe = [
            shed[0] + probe_off[0],
            shed[1] + probe_off[1],
            shed[2] + probe_off[2],
        ];
        let wake_v = probe_velocity(shed, probe, -1.0e9, &pt, 400, 0.5);
        let full = fs_wing::induced_velocity_free(
            probe,
            core::slice::from_ref(&panel),
            &[1.0],
            [13.86, 0.0, 0.0],
        );
        let legs_z = full[2] - bound_only(probe, panel.a, panel.b);
        // Traversal sense: the horseshoe's near leg runs TOWARD the
        // bound vortex (+x); the wake line runs downstream (-x) — the
        // same physical filament with opposite traversal, so the oracle
        // compares against MINUS the leg field.
        assert!(
            (wake_v[2] + legs_z).abs() < 0.02 * legs_z.abs().max(1e-9),
            "cross-kernel mismatch at {probe_off:?}: wake {} vs -legs {}",
            wake_v[2],
            -legs_z
        );
    }
    jlog("cross-kernel", "\"straight_wake_matches\":true");
}

#[test]
fn caps_at_cap_and_cap_plus_one() {
    let (st, pr, gr) = (stations(), probes(), ground());
    let pt = op_point(0.2);
    assert!(assemble_operator(&st, &pr, &gr, &pt, MAX_WAKE_ROWS, 0.5, SPAN_B).is_ok());
    assert_eq!(
        assemble_operator(&st, &pr, &gr, &pt, MAX_WAKE_ROWS + 1, 0.5, SPAN_B)
            .unwrap_err()
            .code,
        "wake-rows-invalid",
        "cap+1 refuses"
    );
    assert_eq!(
        assemble_operator(&st, &pr, &gr, &pt, 0, 0.5, SPAN_B)
            .unwrap_err()
            .code,
        "wake-rows-invalid"
    );
    assert_eq!(
        assemble_operator(&[], &pr, &gr, &pt, 24, 0.5, SPAN_B)
            .unwrap_err()
            .code,
        "wake-grid-empty"
    );
    // Below-ground refusal: h/b that puts shed points under the plane.
    let below = WakeOperatingPoint {
        h_over_b: -0.1,
        ..pt
    };
    assert_eq!(
        assemble_operator(&st, &pr, &gr, &below, 24, 0.5, SPAN_B)
            .unwrap_err()
            .code,
        "wake-geometry-invalid"
    );
    // Uncertified ground refuses (pass-through of the E4.4a law).
    let bad_ground = CertifiedGround {
        z_m: 0.0,
        certificate_slope: 1.0,
        certificate_rms_m: 50.0,
    };
    assert!(assemble_operator(&st, &pr, &bad_ground, &pt, 24, 0.5, SPAN_B).is_err());
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn golden_digest() {
    let (st, pr, gr) = (stations(), probes(), ground());
    let g = frozen_grid_v1();
    let mut payload = Vec::new();
    payload.extend_from_slice(g.digest.as_bytes());
    for pt in g.points.iter().step_by(29) {
        let op = assemble_operator(&st, &pr, &gr, pt, 24, 0.5, SPAN_B).unwrap();
        for v in &op.w_normal {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-wing.e43b1-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "2d9eaa0a810eb3317f12a55cb60b60538846faa5e21ac8f0239d69a0c9d2b196",
        "prescribed-wake golden moved — determinism regression or an \
         intentional operator/grid change requiring the golden-bump protocol"
    );
}
