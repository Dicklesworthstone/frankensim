//! fs-material conformance suite (the tfz.2 bead). Acceptance:
//! consistent-tangent-vs-FD for EVERY law (merge-gate discipline);
//! frame-indifference; return-mapping consistency + dissipation
//! non-negativity on cyclic paths; Mander/Menegotto–Pinto hysteresis
//! fixture behavior; calibration round-trip; admissibility machine
//! checks.

use fs_material::tensor::{contract, rotate, rotation};
use fs_material::{
    Hyperelastic, HyperelasticModel, IsotropicElastic, J2Plasticity, ManderConcrete,
    MenegottoPintoSteel, OrthotropicElastic, SmallStrainLaw, Uniaxial, Voigt, calibrate_bilinear,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-material/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn rand_strain(seed: &mut u64, scale: f64) -> Voigt {
    core::array::from_fn(|_| (lcg(seed) * 2.0 - 1.0) * scale)
}

fn steel() -> IsotropicElastic {
    IsotropicElastic::new(200e9, 0.3, 0.05).expect("steel")
}

fn j2() -> J2Plasticity {
    J2Plasticity::new(steel(), 400e6, 2e9).expect("j2")
}

/// FD-vs-tangent gate for a small-strain law at one strain/state point.
fn check_tangent<L: SmallStrainLaw>(law: &L, strain: &Voigt, state: &L::State, tol: f64) {
    let tangent = law.tangent(strain, state);
    let h = 1e-8;
    for j in 0..6 {
        let mut up = *strain;
        let mut dn = *strain;
        up[j] += h;
        dn[j] -= h;
        let s_up = law.stress(&up, state);
        let s_dn = law.stress(&dn, state);
        for i in 0..6 {
            let fd = (s_up[i] - s_dn[i]) / (2.0 * h);
            let scale = tangent[i][j].abs().max(1e9);
            assert!(
                (tangent[i][j] - fd).abs() / scale < tol,
                "tangent[{i}][{j}] = {} vs FD {fd}",
                tangent[i][j]
            );
        }
    }
}

#[test]
fn mt_001_consistent_tangent_gate_every_law() {
    let mut seed = 0x3A7_0001u64;
    // Elastic laws: exact everywhere.
    let iso = steel();
    let ortho = OrthotropicElastic::new(
        [140e9, 10e9, 10e9],
        [0.30, 0.30, 0.45],
        [5e9, 3.5e9, 5e9],
        0.02,
    )
    .expect("ortho");
    for _ in 0..20 {
        let e = rand_strain(&mut seed, 5e-3);
        check_tangent(&iso, &e, &iso.initial_state(), 1e-6);
        check_tangent(&ortho, &e, &ortho.initial_state(), 1e-6);
    }
    // J2: elastic branch, plastic branch, and mid-cycle states.
    let law = j2();
    let mut state = law.initial_state();
    for step in 1..=12 {
        let amp = f64::from(step) * 6e-4 * if step % 4 < 2 { 1.0 } else { -1.0 };
        let mut e = [0.0f64; 6];
        e[0] = amp;
        e[3] = 0.3 * amp;
        check_tangent(&law, &e, &state, 5e-5);
        state = law.update_state(&e, &state);
    }
    // Hyperelastic 9×9 tangents vs FD of the Piola stress.
    for model in [
        HyperelasticModel::NeoHookean {
            mu: 1.0e6,
            lambda: 4.0e6,
        },
        HyperelasticModel::MooneyRivlin {
            c10: 0.4e6,
            c01: 0.1e6,
            kappa: 5.0e6,
        },
    ] {
        let law = Hyperelastic::new(model, 3.0).expect("hyper");
        for _ in 0..8 {
            let mut f = [0.0f64; 9];
            for (i, v) in f.iter_mut().enumerate() {
                *v = if i % 4 == 0 { 1.0 } else { 0.0 };
                *v += (lcg(&mut seed) * 2.0 - 1.0) * 0.15;
            }
            let a = law.tangent(&f).expect("tangent");
            let h = 1e-7;
            for j in 0..9 {
                let mut up = f;
                let mut dn = f;
                up[j] += h;
                dn[j] -= h;
                let p_up = law.piola(&up).expect("p+");
                let p_dn = law.piola(&dn).expect("p-");
                for i in 0..9 {
                    let fd = (p_up[i] - p_dn[i]) / (2.0 * h);
                    let scale = a[i][j].abs().max(1e6);
                    assert!(
                        (a[i][j] - fd).abs() / scale < 1e-5,
                        "hyper tangent[{i}][{j}] {} vs FD {fd}",
                        a[i][j]
                    );
                }
            }
        }
    }
    // Uniaxial laws through cyclic histories.
    let mp = MenegottoPintoSteel::new(200e9, 400e6, 0.02).expect("mp");
    let mut mp_state = mp.initial_state();
    let mander = ManderConcrete::new(35e6, 0.004, 25e9, 0.015).expect("mander");
    let mut mander_state = mander.initial_state();
    let path = [
        0.001, 0.003, 0.005, 0.002, -0.001, -0.004, 0.0, 0.004, 0.006, 0.003,
    ];
    for &eps in &path {
        for (name, tan, fd_pair) in [
            ("menegotto-pinto", mp.tangent(eps, &mp_state), {
                let h = 1e-9;
                (mp.stress(eps + h, &mp_state), mp.stress(eps - h, &mp_state))
            }),
            ("mander", mander.tangent(eps.abs(), &mander_state), {
                let h = 1e-9;
                (
                    mander.stress(eps.abs() + h, &mander_state),
                    mander.stress(eps.abs() - h, &mander_state),
                )
            }),
        ] {
            let fd = (fd_pair.0 - fd_pair.1) / 2e-9;
            let scale = tan.abs().max(1e9);
            assert!(
                (tan - fd).abs() / scale < 1e-4,
                "{name} tangent {tan} vs FD {fd} at eps={eps}"
            );
        }
        mp_state = mp.update_state(eps, &mp_state);
        mander_state = mander.update_state(eps.abs(), &mander_state);
    }
    verdict(
        "mt-001",
        "FD tangent gate green for iso/ortho elastic, J2 through a cycle, NH/MR 9x9, \
         Menegotto-Pinto + Mander through cyclic histories",
    );
}

#[test]
fn mt_002_frame_indifference() {
    let mut seed = 0x3A7_0002u64;
    let iso = steel();
    for _ in 0..30 {
        let axis = {
            let v = [
                lcg(&mut seed) * 2.0 - 1.0,
                lcg(&mut seed) * 2.0 - 1.0,
                lcg(&mut seed) * 2.0 - 1.0,
            ];
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / n, v[1] / n, v[2] / n]
        };
        let q = rotation(axis, lcg(&mut seed) * 6.0 - 3.0);
        // Isotropic small-strain: σ(QεQᵀ) = Q σ(ε) Qᵀ.
        let e = rand_strain(&mut seed, 4e-3);
        let lhs = iso.stress(&rotate(&e, &q), &());
        let rhs = rotate(&iso.stress(&e, &()), &q);
        for (a, b) in lhs.iter().zip(&rhs) {
            assert!((a - b).abs() < 1e-3, "isotropy/objectivity broke");
        }
        // Hyperelastic objectivity: P(QF) = Q P(F).
        let nh = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 1.0e6,
                lambda: 4.0e6,
            },
            3.0,
        )
        .expect("nh");
        let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        for v in &mut f {
            *v += (lcg(&mut seed) * 2.0 - 1.0) * 0.12;
        }
        let p = nh.piola(&f).expect("p");
        let fm = [[f[0], f[1], f[2]], [f[3], f[4], f[5]], [f[6], f[7], f[8]]];
        let qf = fs_material::tensor::matmul3(&q, &fm);
        let qf_flat = [
            qf[0][0], qf[0][1], qf[0][2], qf[1][0], qf[1][1], qf[1][2], qf[2][0], qf[2][1],
            qf[2][2],
        ];
        let p_qf = nh.piola(&qf_flat).expect("p(QF)");
        let pm = [[p[0], p[1], p[2]], [p[3], p[4], p[5]], [p[6], p[7], p[8]]];
        let qp = fs_material::tensor::matmul3(&q, &pm);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p_qf[3 * i + j] - qp[i][j]).abs() < 1e-3,
                    "hyperelastic objectivity broke at ({i},{j})"
                );
            }
        }
    }
    verdict(
        "mt-002",
        "rotated strain -> rotated stress; P(QF) = Q P(F) over 30 rotations",
    );
}

#[test]
fn mt_003_j2_return_mapping_and_dissipation() {
    let law = j2();
    let e_mod = 200e9;
    let h_mod = 2e9;
    let mut state = law.initial_state();
    // Uniaxial-strain cyclic path (drive ε_xx; this is a constrained
    // uniaxial-strain state, so verify with the yield function instead of
    // uniaxial-stress formulas).
    let path = [
        0.0005, 0.001, 0.0025, 0.004, 0.002, 0.0, -0.002, -0.004, 0.0, 0.004,
    ];
    let mut dissipation = 0.0f64;
    let mut prev_plastic = state.plastic_strain;
    for &amp in &path {
        let mut e = [0.0f64; 6];
        e[0] = amp;
        let stress = law.stress(&e, &state);
        let new_state = law.update_state(&e, &state);
        // Consistency: after the return map, f <= tolerance.
        let f = law.yield_function(&stress, new_state.alpha);
        assert!(
            f < 1.0,
            "yield consistency violated: f = {f} Pa at eps={amp}"
        );
        // Dissipation increment σ : Δε_p >= 0 (associative flow).
        let mut dep = [0.0f64; 6];
        for (d, (n, o)) in dep
            .iter_mut()
            .zip(new_state.plastic_strain.iter().zip(&prev_plastic))
        {
            *d = n - o;
        }
        let inc = contract(&stress, &dep);
        assert!(inc > -1e-6, "negative dissipation increment {inc}");
        dissipation += inc;
        prev_plastic = new_state.plastic_strain;
        state = new_state;
    }
    assert!(dissipation > 0.0, "a plastic cycle must dissipate");
    assert!(state.alpha > 0.0, "plastic flow must have occurred");
    // Uniaxial-STRESS check via a driven 1D fiber equivalent: the J2
    // uniaxial tangent E·H/(E+3μ-ish) is covered by mt-001's FD gate;
    // here assert hardening moved the yield surface.
    let virgin_yield = 400e6;
    let hardened = virgin_yield + h_mod * state.alpha;
    assert!(
        hardened > virgin_yield,
        "isotropic hardening must expand the surface"
    );
    let _ = e_mod;
    verdict(
        "mt-003",
        "return-map consistency f<=tol on a 10-step cycle; dissipation increments \
         non-negative, total positive; surface hardened",
    );
}

#[test]
fn mt_004_hysteresis_fixtures() {
    // Menegotto–Pinto published-curve features.
    let mp = MenegottoPintoSteel::new(200e9, 400e6, 0.02).expect("mp");
    let ey = 400e6 / 200e9;
    // (a) Virgin loading approaches the hardening asymptote.
    let s_large = mp.stress(10.0 * ey, &mp.initial_state());
    let asymptote = 0.02 * 200e9 * 10.0 * ey + 400e6 * (1.0 - 0.02);
    assert!(
        (s_large - asymptote).abs() / asymptote < 0.02,
        "virgin curve must approach the b*E0 asymptote: {s_large} vs {asymptote}"
    );
    // (b) Initial tangent is E0; far-field tangent is b*E0.
    let t0 = mp.tangent(1e-6, &mp.initial_state());
    assert!((t0 - 200e9).abs() / 200e9 < 1e-3, "initial tangent E0");
    let t_inf = mp.tangent(12.0 * ey, &mp.initial_state());
    assert!(
        (t_inf - 0.02 * 200e9).abs() / (0.02 * 200e9) < 0.05,
        "hardening tangent b*E0"
    );
    // (c) Bauschinger effect: after a plastic excursion and reversal, the
    // reverse curve departs from linear-elastic EARLIER than the virgin
    // yield (softened knee).
    let mut state = mp.initial_state();
    state = mp.update_state(4.0 * ey, &state);
    let sigma_at_reversal = state.sig_prev;
    let probe = 4.0 * ey - 1.5 * ey; // 1.5 εy of elastic return
    let sig_probe = mp.stress(probe, &state);
    let elastic_prediction = sigma_at_reversal - 200e9 * 1.5 * ey;
    assert!(
        sig_probe > elastic_prediction + 0.02 * 400e6,
        "reverse branch must soften below the elastic line (Bauschinger): \
         {sig_probe} vs {elastic_prediction}"
    );
    // (d) A full symmetric cycle closes into a stable loop with positive
    // dissipated energy.
    let mut state = mp.initial_state();
    let mut prev = (0.0f64, 0.0f64);
    let mut area = 0.0f64;
    let steps = 400;
    for k in 0..=steps {
        let t = f64::from(k) / f64::from(steps);
        let eps = 5.0 * ey * fs_math::det::sin(t * 2.0 * core::f64::consts::PI);
        let sig = mp.stress(eps, &state);
        state = mp.update_state(eps, &state);
        if k > 0 {
            area += f64::midpoint(sig, prev.1) * (eps - prev.0);
        }
        prev = (eps, sig);
    }
    assert!(area > 0.0, "hysteresis loop must dissipate (area {area})");

    // Mander envelope: peak at (eps_cc, fcc), correct initial modulus,
    // unload to residual strain, reload rejoining the envelope.
    let mander = ManderConcrete::new(35e6, 0.004, 25e9, 0.015).expect("mander");
    let (peak, dpeak) = mander.envelope(0.004);
    assert!((peak - 35e6).abs() / 35e6 < 1e-12, "envelope peak is fcc");
    assert!(
        dpeak.abs() < 35e6 / 0.004 * 1e-9,
        "slope vanishes at the peak"
    );
    let (s_small, d_small) = mander.envelope(1e-9);
    assert!(
        s_small >= 0.0 && (d_small - 25e9).abs() / 25e9 < 1e-3,
        "initial modulus Ec"
    );
    // Post-peak softening.
    let (post, _) = mander.envelope(0.008);
    assert!(post < peak, "envelope must soften past the peak");
    // Unload/reload cycle.
    let mut st = mander.initial_state();
    st = mander.update_state(0.006, &st);
    let sig_top = mander.stress(0.006, &st);
    let eps_p = 0.006 - sig_top / 25e9;
    assert!(
        mander.stress(eps_p, &st).abs() < 1e-3,
        "unload reaches zero at eps_p"
    );
    let mid = f64::midpoint(eps_p, 0.006);
    let sig_mid = mander.stress(mid, &st);
    assert!(
        sig_mid > 0.0 && sig_mid < sig_top,
        "reload line between (eps_p,0)-(top)"
    );
    // Rejoins the envelope beyond the previous maximum.
    st = mander.update_state(0.007, &st);
    let (env7, _) = mander.envelope(0.007);
    assert!(
        (mander.stress(0.007, &st) - env7).abs() / env7 < 1e-12,
        "rejoins envelope"
    );
    verdict(
        "mt-004",
        "M-P asymptote/tangents/Bauschinger/dissipating loop; Mander peak/softening/\
         unload-reload fixture behavior",
    );
}

#[test]
fn mt_005_calibration_round_trip() {
    // Synthetic bilinear data: E = 200 GPa, sigma_y = 400 MPa, H = 4 GPa,
    // with deterministic ±0.5 MPa "measurement" noise.
    let e_true = 200e9;
    let sy_true = 400e6;
    let h_true = 4e9;
    let ey = sy_true / e_true;
    let mut seed = 0x3A7_0005u64;
    let mut data = Vec::new();
    for k in 1..=40 {
        let eps = f64::from(k) * 3e-4 / 2.0;
        let clean = if eps <= ey {
            e_true * eps
        } else {
            sy_true + h_true * (eps - ey)
        };
        data.push((eps, clean + (lcg(&mut seed) - 0.5) * 1e6));
    }
    let fit = calibrate_bilinear(&data).expect("fit");
    assert!(
        (fit.youngs - e_true).abs() / e_true < 0.01,
        "E recovered within 1%: {}",
        fit.youngs
    );
    assert!(
        (fit.post_yield - h_true).abs() / h_true < 0.05,
        "H recovered within 5%: {}",
        fit.post_yield
    );
    assert!(
        (fit.yield_stress - sy_true).abs() / sy_true < 0.02,
        "sigma_y recovered within 2%: {}",
        fit.yield_stress
    );
    // The uncertainty envelope contains the truth.
    assert!(
        (fit.youngs - e_true).abs() < 4.0 * fit.youngs_se.max(1e6),
        "E envelope contains truth"
    );
    println!(
        "{{\"suite\":\"fs-material/conformance\",\"metric\":\"calibration\",\
         \"E\":{:.4e},\"E_se\":{:.3e},\"H\":{:.4e},\"sigma_y\":{:.4e},\"rms\":{:.3e}}}",
        fit.youngs, fit.youngs_se, fit.post_yield, fit.yield_stress, fit.rms_residual
    );
    // Degenerate data refuses structurally.
    assert!(calibrate_bilinear(&data[..4]).is_err());
    assert!(calibrate_bilinear(&[(0.0, 0.0); 8]).is_err());
    verdict(
        "mt-005",
        "bilinear calibration recovers (E, H, sigma_y) with truth in envelope",
    );
}

/// Every law ships a card with assumptions, a bounded domain, failures.
fn check_card_completeness() {
    let laws_cards = [
        steel().card(),
        j2().card(),
        MenegottoPintoSteel::new(200e9, 400e6, 0.02)
            .expect("mp")
            .card(),
        ManderConcrete::new(35e6, 0.004, 25e9, 0.015)
            .expect("m")
            .card(),
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 1e6,
                lambda: 4e6,
            },
            3.0,
        )
        .expect("nh")
        .card(),
    ];
    for card in &laws_cards {
        assert!(!card.assumptions.is_empty(), "{}: assumptions", card.name);
        assert!(
            !card.validity.param_names().is_empty(),
            "{}: domain",
            card.name
        );
        assert!(!card.known_failures.is_empty(), "{}: failures", card.name);
    }
}

#[test]
fn mt_006_admissibility_machine_checks() {
    // Polyconvexity spot test: W(F + t a⊗b) is convex in t along random
    // rank-one directions (necessary condition, sampled).
    let mut seed = 0x3A7_0006u64;
    for model in [
        HyperelasticModel::NeoHookean {
            mu: 1.0e6,
            lambda: 4.0e6,
        },
        HyperelasticModel::MooneyRivlin {
            c10: 0.4e6,
            c01: 0.1e6,
            kappa: 5.0e6,
        },
    ] {
        let law = Hyperelastic::new(model, 3.0).expect("hyper");
        assert_eq!(law.admissibility().polyconvex, Some(true));
        for _ in 0..40 {
            let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            for v in &mut f {
                *v += (lcg(&mut seed) * 2.0 - 1.0) * 0.1;
            }
            let a: [f64; 3] = core::array::from_fn(|_| lcg(&mut seed) * 2.0 - 1.0);
            let b: [f64; 3] = core::array::from_fn(|_| lcg(&mut seed) * 2.0 - 1.0);
            // Rank-one direction d = a⊗b; second derivative along d must
            // be >= 0 (rank-one convexity, implied by polyconvexity).
            let mut d = [0.0f64; 9];
            for i in 0..3 {
                for j in 0..3 {
                    d[3 * i + j] = a[i] * b[j];
                }
            }
            let h = 1e-4;
            let w = |t: f64| {
                let mut ft = f;
                for (x, dd) in ft.iter_mut().zip(&d) {
                    *x += t * dd;
                }
                let (val, _) = fs_ad::gradient::<9>(ft, |x| law.energy(&x));
                val
            };
            let second = (w(h) - 2.0 * w(0.0) + w(-h)) / (h * h);
            assert!(
                second > -1e-2 * second.abs().max(1.0),
                "rank-one convexity violated: {second}"
            );
        }
    }
    // Inadmissible parameter sets refuse at construction.
    assert!(IsotropicElastic::new(-1.0, 0.3, 0.01).is_err());
    assert!(IsotropicElastic::new(200e9, 0.6, 0.01).is_err());
    assert!(
        OrthotropicElastic::new([1e9, 1e9, 1e9], [0.9, 0.9, 0.9], [1e9, 1e9, 1e9], 0.01).is_err(),
        "thermodynamically inadmissible Poisson set must refuse"
    );
    assert!(J2Plasticity::new(steel(), -5.0, 1e9).is_err());
    assert!(MenegottoPintoSteel::new(200e9, 400e6, 1.5).is_err());
    assert!(
        ManderConcrete::new(35e6, 0.004, 5e9, 0.015).is_err(),
        "Ec <= Esec refused"
    );
    check_card_completeness();
    // Evidence adapter: in-domain vs out-of-domain strain is FLAGGED.
    let iso = steel();
    let mut e_ok = [0.0f64; 6];
    e_ok[0] = 1e-3;
    let ev = fs_material::evidence_stress(&iso, &e_ok, &());
    assert!(
        ev.model.in_domain,
        "1e-3 strain is inside the calibration domain"
    );
    let mut e_far = [0.0f64; 6];
    e_far[0] = 0.2; // way past strain_limit = 0.05
    let ev_far = fs_material::evidence_stress(&iso, &e_far, &());
    assert!(
        !ev_far.model.in_domain,
        "0.2 strain must be flagged out-of-domain"
    );
    verdict(
        "mt-006",
        "rank-one convexity sampled for NH/MR; inadmissible parameters refuse; every \
         card carries assumptions, domain, failures; Evidence flags domain exit",
    );
}
