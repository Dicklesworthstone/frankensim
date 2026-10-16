//! V-04b1 analytic-half battery (bead wf-root-guzez.4.6.2.1,
//! E3.3b-ii-a): the Γ→0 isotropic limit against an INDEPENDENTLY coded
//! von Kármán tensor (with a Γ>0 liveness twin), incompressibility of
//! the SHEARED tensor, Hermitian symmetry + positive semidefiniteness,
//! the k₁→0 branch continuity, hypergeometric cross-path oracle,
//! Reynolds-stress anisotropy + the uw<0 shear stress (the quantity a
//! diagonal model cannot produce), caps at cap AND cap+1, golden.
//! Repro: cargo test -p fs-atmo --test mann_battery

use fs_atmo::mann::{
    MANN_TARGET_V1, MAX_GAMMA, MIN_KL, MannParams, hypergeom_mann, isotropic_tensor, mann_tensor,
    stress_integrals,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-atmo-v04b1a\",\"case\":\"{case}\",{payload}}}");
}

const SAMPLE_K: [[f64; 3]; 6] = [
    [0.7, 0.3, 0.5],
    [1.5, -0.8, 0.2],
    [0.1, 0.0, 1.1],
    [-0.6, 0.9, -0.4],
    [2.5, 1.5, 3.0],
    [0.05, 0.4, 0.02],
];

#[test]
fn gamma_zero_reduces_exactly_to_von_karman() {
    let iso = MannParams {
        gamma: 0.0,
        ..MANN_TARGET_V1
    };
    for k in SAMPLE_K {
        let a = mann_tensor(&iso, k).unwrap();
        let b = isotropic_tensor(&iso, k);
        for i in 0..3 {
            for j in 0..3 {
                let scale = b[i][j].abs().max(1e-12);
                assert!(
                    (a[i][j] - b[i][j]).abs() < 1e-10 * scale,
                    "iso limit mismatch at {k:?}[{i}][{j}]: {} vs {}",
                    a[i][j],
                    b[i][j]
                );
            }
        }
    }
    // Liveness twin: Γ = 3.9 must DIFFER from isotropic (the limit test
    // would be vacuous if the shear never did anything).
    let k = SAMPLE_K[0];
    let sheared = mann_tensor(&MANN_TARGET_V1, k).unwrap();
    let iso_t = isotropic_tensor(&MANN_TARGET_V1, k);
    assert!(
        (sheared[0][2] - iso_t[0][2]).abs() > 1e-6 * iso_t[0][0].abs(),
        "shear must move the 13 component"
    );
    jlog("iso-limit", "\"exact_with_liveness_twin\":true");
}

#[test]
fn sheared_tensor_is_incompressible_symmetric_and_psd() {
    for gamma in [0.0, 3.9] {
        let p = MannParams {
            gamma,
            ..MANN_TARGET_V1
        };
        for k in SAMPLE_K {
            let phi = mann_tensor(&p, k).unwrap();
            let scale = (0..3).map(|i| phi[i][i].abs()).fold(0.0f64, f64::max)
                * (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
            // Continuity: Φ·k = 0 (RDT preserves incompressibility with
            // the CURRENT wavevector) — the strongest ζ-algebra oracle.
            for i in 0..3 {
                let div: f64 = (0..3).map(|j| phi[i][j] * k[j]).sum();
                assert!(
                    div.abs() < 1e-9 * scale.max(1e-12),
                    "continuity violated (gamma {gamma}) at {k:?} row {i}: {div}"
                );
            }
            // Symmetry (the tensor is real: Hermitian = symmetric here).
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (phi[i][j] - phi[j][i]).abs() < 1e-14 * scale.max(1e-12),
                        "symmetry at {k:?}"
                    );
                }
            }
            // PSD: diagonal nonnegative + Cauchy-Schwarz on off-diagonals.
            for i in 0..3 {
                assert!(phi[i][i] >= -1e-15, "diagonal must be >= 0");
                for j in 0..3 {
                    if i != j {
                        assert!(
                            phi[i][j] * phi[i][j] <= phi[i][i] * phi[j][j] * (1.0 + 1e-9) + 1e-30,
                            "Cauchy-Schwarz at {k:?} ({i},{j})"
                        );
                    }
                }
            }
        }
    }
    jlog("structure", "\"continuity_symmetry_psd\":true");
}

#[test]
fn k1_limit_branch_is_continuous() {
    // The k₁ → 0 branch (ζ1 → −β) must join the regular path smoothly.
    let k_reg = [1e-7, 0.9, 0.7];
    let k_lim = [0.0, 0.9, 0.7];
    let a = mann_tensor(&MANN_TARGET_V1, k_reg).unwrap();
    let b = mann_tensor(&MANN_TARGET_V1, k_lim).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let scale = b[i][j].abs().max(1e-9);
            assert!(
                (a[i][j] - b[i][j]).abs() < 1e-4 * scale,
                "k1 branch discontinuity [{i}][{j}]: {} vs {}",
                a[i][j],
                b[i][j]
            );
        }
    }
    jlog("k1-branch", "\"continuous\":true");
}

#[test]
fn hypergeom_cross_path_oracle() {
    // Independent path: for |z| < 1 the DIRECT Gauss series converges —
    // compare against the Pfaff-transformed evaluation.
    let direct = |z: f64| -> f64 {
        let (a, b, c) = (1.0 / 3.0, 17.0 / 6.0, 4.0 / 3.0);
        let mut term = 1.0f64;
        let mut sum = 1.0f64;
        for n in 0..500 {
            let nf = n as f64;
            term *= (a + nf) * (b + nf) / (c + nf) * z / (nf + 1.0);
            sum += term;
        }
        sum
    };
    for z in [-0.05, -0.2, -0.5, -0.8] {
        let p = hypergeom_mann(z).unwrap();
        let d = direct(z);
        assert!(
            (p - d).abs() < 1e-10 * d.abs(),
            "2F1({z}): pfaff {p} vs direct {d}"
        );
    }
    assert!(
        (hypergeom_mann(0.0).unwrap() - 1.0).abs() < 1e-14,
        "2F1(0) = 1"
    );
    jlog("hypergeom", "\"cross_path\":true");
}

#[test]
fn stresses_show_surface_layer_anisotropy_and_negative_uw() {
    let s = stress_integrals(&MANN_TARGET_V1).unwrap();
    let (u2, v2, w2, uw) = (s[0][0], s[1][1], s[2][2], s[0][2]);
    assert!(u2 > 0.0 && v2 > 0.0 && w2 > 0.0);
    assert!(
        u2 > v2 && v2 > w2,
        "surface-layer ordering u2>v2>w2: {u2} {v2} {w2}"
    );
    assert!(
        (0.4..0.95).contains(&(v2 / u2)),
        "v2/u2 {} outside the Mann-class band",
        v2 / u2
    );
    assert!(
        (0.2..0.8).contains(&(w2 / u2)),
        "w2/u2 {} outside the Mann-class band",
        w2 / u2
    );
    assert!(uw < 0.0, "the shear stress must be NEGATIVE: {uw}");
    assert!(
        (0.05..0.5).contains(&(-uw / u2)),
        "-uw/u2 {} outside the band",
        -uw / u2
    );
    // The diagonal-only impossibility twin: the ISOTROPIC tensor's uw
    // integral is zero — the shear distortion is what buys the stress.
    let iso = MannParams {
        gamma: 0.0,
        ..MANN_TARGET_V1
    };
    let s0 = stress_integrals(&iso).unwrap();
    assert!(
        s0[0][2].abs() < 1e-3 * s0[0][0],
        "isotropic uw must vanish: {}",
        s0[0][2]
    );
    let iso_ratio = s0[1][1] / s0[0][0];
    assert!(
        (iso_ratio - 1.0).abs() < 0.05,
        "isotropic must be near-isotropic in the ratios: {iso_ratio}"
    );
    jlog(
        "stresses",
        &format!(
            "\"u2\":{u2},\"v2\":{v2},\"w2\":{w2},\"uw\":{uw},\"v_ratio\":{},\"w_ratio\":{},\"uw_ratio\":{}",
            v2 / u2,
            w2 / u2,
            -uw / u2
        ),
    );
}

#[test]
fn caps_at_cap_and_cap_plus_one() {
    let mk = |gamma: f64| MannParams {
        gamma,
        ..MANN_TARGET_V1
    };
    assert!(mk(MAX_GAMMA).admit().is_ok(), "gamma cap admits");
    assert_eq!(
        mk(MAX_GAMMA.next_up()).admit().unwrap_err().code,
        "mann-params-invalid",
        "gamma cap+1 refuses"
    );
    assert_eq!(mk(-1e-300).admit().unwrap_err().code, "mann-params-invalid");
    // kL floor: at the floor admits, one ulp under refuses.
    let k_at = MIN_KL / MANN_TARGET_V1.length_m;
    assert!(mann_tensor(&MANN_TARGET_V1, [k_at, 0.0, 0.0]).is_ok());
    let k_under = (MIN_KL / MANN_TARGET_V1.length_m) * (1.0 - 1e-16);
    let under = mann_tensor(&MANN_TARGET_V1, [k_under.next_down(), 0.0, 0.0]);
    assert_eq!(under.unwrap_err().code, "mann-wavevector-invalid");
    let bad_len = MannParams {
        length_m: 0.0,
        ..MANN_TARGET_V1
    };
    assert_eq!(bad_len.admit().unwrap_err().code, "mann-params-invalid");
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn determinism_and_golden() {
    let mut payload = Vec::new();
    for k in SAMPLE_K {
        let a = mann_tensor(&MANN_TARGET_V1, k).unwrap();
        let b = mann_tensor(&MANN_TARGET_V1, k).unwrap();
        assert_eq!(a, b, "bitwise repeat");
        for row in a {
            for v in row {
                payload.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
    }
    let s = stress_integrals(&MANN_TARGET_V1).unwrap();
    for row in s {
        for v in row {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-atmo.v04b1a-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "7f5cbd11e9d3493df8b4cfde0859b2e98ca28c3e8f078bdbc6bcdf0f8e5a2097",
        "Mann-tensor golden moved — determinism regression or an \
         intentional target change requiring the golden-bump protocol"
    );
}
