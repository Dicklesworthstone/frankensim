//! V-06a battery (bead wf-root-guzez.5.10, E4.4a): image EXACTNESS —
//! the solved system's z-velocity residual ON the plane at machine
//! precision; the ground-effect lift trend vs h/b (monotone toward the
//! free-air limit); the free-air limit itself; the certificate gate at
//! the band AND the next float; below-ground refusal; golden.
//! Repro: cargo test -p fs-wing --test images_battery

use fs_wing::images::{
    CERT_MAX_RMS_M, CERT_MAX_SLOPE, CertifiedGround, induced_velocity_with_images,
    solve_weissinger_ground,
};
use fs_wing::{SurfaceId, flat_surface, solve_weissinger_linear};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-v06a\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const V: f64 = 13.86;

fn freestream(alpha: f64) -> [f64; 3] {
    [V * alpha.cos(), 0.0, V * alpha.sin()]
}

fn kdh_ground() -> CertifiedGround {
    // The E1.3/E3.4-iii launch-flat certificate numbers.
    CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.000606,
        certificate_rms_m: 0.801,
    }
}

#[test]
fn image_exactness_zero_normal_flow_on_the_plane() {
    // Solve the 1903 biplane 1.2 m above the plane; then probe the TOTAL
    // z-velocity (induced + freestream-normal component is zero for a
    // level plane: freestream w is uniform, the plane is a streamsurface
    // of the IMAGE construction for the induced field) at points ON the
    // plane: the induced w must cancel between real and image systems to
    // machine precision RELATIVE to the induced speed scale.
    let ground = kdh_ground();
    let mut p = flat_surface(
        SurfaceId::WingLower,
        12.29,
        1.981,
        0.0,
        ground.z_m - 1.2,
        8,
        2,
    )
    .unwrap();
    p.extend(
        flat_surface(
            SurfaceId::WingUpper,
            12.29,
            1.981,
            0.0,
            ground.z_m - 3.09,
            8,
            2,
        )
        .unwrap(),
    );
    let fs_v = freestream(0.07);
    let r = solve_weissinger_ground(&p, fs_v, RHO, &ground).unwrap();
    let mut worst_rel = 0.0f64;
    for i in 0..7 {
        for j in 0..5 {
            let probe = [
                -30.0 + 12.0 * f64::from(i),
                -20.0 + 9.5 * f64::from(j),
                ground.z_m,
            ];
            let v = induced_velocity_with_images(probe, &p, &r.gamma, fs_v, ground.z_m);
            let scale = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-12);
            worst_rel = worst_rel.max(v[2].abs() / scale.max(0.01));
        }
    }
    assert!(
        worst_rel < 1e-9,
        "plane residual {worst_rel:e} not machine-exact"
    );
    jlog(
        "exactness",
        &format!("\"worst_rel_w_on_plane\":{worst_rel:e}"),
    );
}

#[test]
fn ground_effect_trend_and_free_air_limit() {
    // Fixed alpha: lift GROWS as the wing nears the certified plane, and
    // approaches the free-air value from above as h/b grows.
    let ground = kdh_ground();
    let fs_v = freestream(0.05);
    let wing_at = |h: f64| {
        flat_surface(
            SurfaceId::WingLower,
            12.29,
            1.981,
            0.0,
            ground.z_m - h,
            8,
            2,
        )
        .unwrap()
    };
    let free = solve_weissinger_linear(&wing_at(1000.0), fs_v, RHO)
        .unwrap()
        .total_lift_n;
    let l_low = solve_weissinger_ground(&wing_at(1.2), fs_v, RHO, &ground)
        .unwrap()
        .total_lift_n;
    let l_mid = solve_weissinger_ground(&wing_at(6.0), fs_v, RHO, &ground)
        .unwrap()
        .total_lift_n;
    let l_high = solve_weissinger_ground(&wing_at(60.0), fs_v, RHO, &ground)
        .unwrap()
        .total_lift_n;
    assert!(
        l_low > l_mid && l_mid > l_high,
        "ground effect must decay monotonically ({l_low}, {l_mid}, {l_high})"
    );
    assert!(
        l_low / free > 1.05,
        "near-ground lift gain {} must be real",
        l_low / free
    );
    assert!(
        (l_high / free - 1.0).abs() < 0.01,
        "far-field must recover free air ({l_high} vs {free})"
    );
    jlog(
        "trend",
        &format!(
            "\"k_low\":{},\"k_mid\":{},\"k_high\":{}",
            l_low / free,
            l_mid / free,
            l_high / free
        ),
    );
}

#[test]
fn certificate_gate_at_band_and_next_float() {
    let fs_v = freestream(0.05);
    let p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 1.5, 8, 2).unwrap();
    let at_band = CertifiedGround {
        z_m: 3.0,
        certificate_slope: CERT_MAX_SLOPE,
        certificate_rms_m: CERT_MAX_RMS_M,
    };
    assert!(
        solve_weissinger_ground(&p, fs_v, RHO, &at_band).is_ok(),
        "AT the band admits"
    );
    let over_slope = CertifiedGround {
        certificate_slope: f64::from_bits(CERT_MAX_SLOPE.to_bits() + 1),
        ..at_band
    };
    assert_eq!(
        solve_weissinger_ground(&p, fs_v, RHO, &over_slope)
            .unwrap_err()
            .code,
        "ground-uncertified"
    );
    let over_rms = CertifiedGround {
        certificate_rms_m: f64::from_bits(CERT_MAX_RMS_M.to_bits() + 1),
        ..at_band
    };
    assert_eq!(
        solve_weissinger_ground(&p, fs_v, RHO, &over_rms)
            .unwrap_err()
            .code,
        "ground-uncertified"
    );
    // Bands must MATCH the fs-flyer prelaunch issuance bands (cross-link).
    assert_eq!(CERT_MAX_SLOPE, 0.005);
    assert_eq!(CERT_MAX_RMS_M, 1.2);
    // Below-ground refusal (frd z-down: panel z beyond the plane).
    let sunk = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 3.5, 8, 2).unwrap();
    assert_eq!(
        solve_weissinger_ground(&sunk, fs_v, RHO, &kdh_ground())
            .unwrap_err()
            .code,
        "aircraft-below-ground"
    );
    jlog(
        "certificate",
        "\"bands\":\"admit at band, refuse next float; cross-linked\"",
    );
}

#[test]
fn images_golden_digest() {
    // The 1903 biplane at the flight height over the certified flat.
    let ground = kdh_ground();
    let mut p = flat_surface(
        SurfaceId::WingLower,
        12.29,
        1.981,
        0.0,
        ground.z_m - 1.2,
        8,
        2,
    )
    .unwrap();
    p.extend(
        flat_surface(
            SurfaceId::WingUpper,
            12.29,
            1.981,
            0.0,
            ground.z_m - 3.09,
            8,
            2,
        )
        .unwrap(),
    );
    let r = solve_weissinger_ground(&p, freestream(0.07), RHO, &ground).unwrap();
    let mut payload = Vec::new();
    for g in &r.gamma {
        payload.extend_from_slice(&g.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&r.total_lift_n.to_bits().to_le_bytes());
    let digest = fs_blake3::hash_domain("org.frankensim.fs-wing.v06a-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"lift\":{}", r.total_lift_n),
    );
    assert_eq!(
        digest, "9e7686445dac182f563826d8631baafb123b01179b3389e9d4da83bb5d37f58b",
        "image golden moved — determinism regression or an intentional \
         kernel change requiring the golden-bump protocol"
    );
}
