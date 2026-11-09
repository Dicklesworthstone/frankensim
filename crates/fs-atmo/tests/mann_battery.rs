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
fn hypergeom_converges_across_the_energy_scale_and_branch_join() {
    // Independent Euler integral (DLMF 15.6.1, https://dlmf.nist.gov/15.6.E1):
    // interchange a,b and put t=s^3 to obtain integral_0^1 (1-z*s^3)^(-17/6) ds.
    // Composite Simpson integrates a smooth positive function without either
    // production transformation, its series, or its gamma-function constants.
    let integral = |z: f64| {
        const BINS: usize = 4096;
        let mut sum = 1.0 + (1.0 - z).powf(-17.0 / 6.0);
        for i in 1..BINS {
            let s = i as f64 / BINS as f64;
            let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
            sum += weight * (1.0 - z * s * s * s).powf(-17.0 / 6.0);
        }
        sum / (3 * BINS) as f64
    };
    for z in [
        -50.0,
        -2.001,
        (-2.0_f64).next_down(),
        -2.0,
        (-2.0_f64).next_up(),
        -1.9,
        -1.1,
        -1.02,
        -1.01,
        -1.005,
        -1.001,
        -1.000_001,
        (-1.0_f64).next_down(),
        -1.0,
        (-1.0_f64).next_up(),
        -0.5,
        0.0,
    ] {
        let actual = hypergeom_mann(z).unwrap();
        let expected = integral(z);
        assert!(
            (actual - expected).abs() < 2e-12 * expected,
            "2F1({z}): transformed series {actual} versus Euler integral {expected}"
        );
    }
    for boundary in [-2.0_f64, -1.0] {
        let left = hypergeom_mann(boundary.next_down()).unwrap();
        let right = hypergeom_mann(boundary.next_up()).unwrap();
        assert!((left - right).abs() < 1e-13 * right);
    }
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
fn stress_quadrature_resolves_shear_through_the_admission_cap() {
    // Independent NumPy/SciPy reference: Gauss-Legendre integration of
    // Mann's tensor using scipy.special.hyp2f1 and ordinary libm, with
    // 192 log-radius × 192 polar-angle × 96 quadrant-azimuth nodes.
    // Refinement to 256×256×128 changes these values by <2e-12 relative.
    // Both retain the declared band 0.01 <= |k|L <= 1000. These fixtures
    // test the integral's physical accuracy, independently of its hash.
    let cases = [
        (
            0.0,
            [
                0.317_625_560_911_890_37,
                0.317_625_560_911_890_37,
                0.317_625_560_911_891_04,
                0.0,
            ],
        ),
        (
            1.0,
            [
                0.367_994_061_030_218_14,
                0.335_403_756_677_867_26,
                0.330_752_124_964_028_84,
                -0.075_196_526_288_425_81,
            ],
        ),
        (
            3.9,
            [
                1.027_819_833_051_831_6,
                0.522_071_024_212_877_9,
                0.277_147_452_086_774_2,
                -0.249_150_708_943_220_9,
            ],
        ),
        (
            7.0,
            [
                2.159_201_107_275_951_6,
                0.692_647_913_447_887_3,
                0.216_598_473_670_035_35,
                -0.348_515_664_680_424_16,
            ],
        ),
        (
            MAX_GAMMA,
            [
                3.444_100_052_670_84,
                0.818_790_914_962_963_4,
                0.179_016_908_145_662_65,
                -0.406_421_363_317_368_7,
            ],
        ),
    ];
    for (gamma, expected) in cases {
        let stress = stress_integrals(&MannParams {
            gamma,
            ..MANN_TARGET_V1
        })
        .unwrap();
        let actual = [stress[0][0], stress[1][1], stress[2][2], stress[0][2]];
        for i in 0..4 {
            assert!(
                (actual[i] - expected[i]).abs() <= 2e-6 * expected[i].abs().max(1e-12),
                "Gamma {gamma}, stress {i}: {} versus converged {}",
                actual[i],
                expected[i]
            );
        }
        assert_eq!(stress[0][1], 0.0);
        assert_eq!(stress[1][0], 0.0);
        assert_eq!(stress[1][2], 0.0);
        assert_eq!(stress[2][1], 0.0);
        assert_eq!(stress[0][2], stress[2][0]);
        assert!(stress[0][2] * stress[0][2] <= stress[0][0] * stress[2][2]);
    }
}

#[test]
fn isotropic_stresses_match_independent_radial_energy_integral() {
    // det-ok: the oracle deliberately uses platform libm and composite
    // Simpson integration, sharing neither the production Gauss nodes nor
    // its dimensionless spectral-tensor reduction.
    let p = MannParams {
        gamma: 0.0,
        ..MANN_TARGET_V1
    };
    let lo = (1e-2 / p.length_m).ln();
    let hi = (1e3 / p.length_m).ln();
    let integral = |bins: usize| {
        let step = (hi - lo) / bins as f64;
        let mut sum = 0.0;
        for i in 0..=bins {
            let k = (lo + step * i as f64).exp();
            let kl = k * p.length_m;
            let energy = p.alpha_eps23 * p.length_m.powf(5.0 / 3.0) * kl.powi(4)
                / (1.0 + kl * kl).powf(17.0 / 6.0);
            let weight = if i == 0 || i == bins {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * energy * k;
        }
        (2.0 / 3.0) * sum * step / 3.0
    };
    let reference = integral(8192);
    assert!((reference - integral(4096)).abs() < 1e-12 * reference);
    let stress = stress_integrals(&p).unwrap();
    for (i, row) in stress.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if i == j {
                assert!((value - reference).abs() < 2e-10 * reference);
                assert_eq!(value, stress[0][0]);
            } else {
                assert_eq!(value, 0.0);
            }
        }
    }
}

#[test]
fn stress_scaling_is_independent_of_dimensional_wavevector_conditioning() {
    for gamma in [0.0, 3.9, MAX_GAMMA] {
        let p = MannParams {
            gamma,
            ..MANN_TARGET_V1
        };
        let original = stress_integrals(&p).unwrap();
        let rescaled = stress_integrals(&MannParams {
            alpha_eps23: 2.5 * p.alpha_eps23,
            length_m: 8.0 * p.length_m,
            ..p
        })
        .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = 10.0 * original[i][j];
                assert!((rescaled[i][j] - expected).abs() <= 2e-14 * expected.abs());
            }
        }
    }
    let original = stress_integrals(&MANN_TARGET_V1).unwrap();
    for length_m in [1e-180, 1e180] {
        // det-ok: independent unit-rescaling oracle. Dimensional k² and
        // k0^4 overflow/underflow here even though every stress is ordinary.
        let alpha_eps23 =
            MANN_TARGET_V1.alpha_eps23 * (MANN_TARGET_V1.length_m / length_m).powf(2.0 / 3.0);
        let rescaled = stress_integrals(&MannParams {
            alpha_eps23,
            length_m,
            ..MANN_TARGET_V1
        })
        .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = original[i][j];
                assert!((rescaled[i][j] - expected).abs() <= 3e-13 * expected.abs());
            }
        }
    }
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
        digest, "488b5eda04b8c2fd7321bdee3458d388bbb75c1cfe48617cd798f284880a52f1",
        "Mann-tensor golden moved — determinism regression or an \
         intentional target change requiring the golden-bump protocol"
    );
}
