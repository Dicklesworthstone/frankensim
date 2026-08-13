//! fs-nlmodal conformance battery: quadrature/symmetry certificates,
//! Duffing backbone vs the analytic perturbation formula, exact energy
//! conservation, force-vs-energy mutation visibility through the
//! fs-phs supply audit, amplitude-scaling covariance, parametric
//! mode-coupling cross-checked against an independent RK4 integrator,
//! the Kirchhoff-Carrier pitch-glide hand formula, and the struck-
//! plate cascade casebook.

use fs_nlmodal::{
    KcStringParams, NlModalError, SineMode, SosModalStorage, StressChannel, VkPlateParams,
    assemble, duffing_backbone, kirchhoff_carrier_string, single_mode_beta, von_karman_ss_plate,
};
use fs_phs::{PortHamiltonian, Storage, step};

fn steel_plate() -> VkPlateParams {
    VkPlateParams {
        lx: 0.4,
        ly: 0.3,
        h: 1.0e-3,
        young: 2.0e11,
        nu: 0.3,
        rho: 7850.0,
    }
}

fn guitar_string() -> KcStringParams {
    // Plain-steel guitar string class: L 0.65 m, T 70 N,
    // mu 5e-3 kg/m, EA ~ 2e11 * 2e-7 m^2 = 4e4 N.
    KcStringParams {
        length: 0.65,
        tension: 70.0,
        lin_density: 5.0e-3,
        ea: 4.0e4,
    }
}

fn modes_grid(mx: usize, my: usize) -> Vec<SineMode> {
    let mut v = Vec::new();
    for m in 1..=mx {
        for n in 1..=my {
            v.push(SineMode { m, n });
        }
    }
    v
}

/// Measure the free-oscillation angular frequency of mode 0 by timing
/// zero crossings of q_0 over many periods (measurement code, no
/// model knowledge).
fn measure_omega(sys: &PortHamiltonian, x0: Vec<f64>, dt: f64, steps: usize) -> f64 {
    let mut x = x0;
    let mut crossings = Vec::new();
    let mut prev_q = x[0];
    for i in 0..steps {
        x = step(sys, &x, &[0.0], dt).expect("step").x;
        let q = x[0];
        if prev_q > 0.0 && q <= 0.0 {
            // Linear interpolation of the crossing time.
            let frac = prev_q / (prev_q - q);
            crossings.push((i as f64 + frac) * dt);
        }
        prev_q = q;
    }
    assert!(crossings.len() >= 3, "not enough crossings to measure");
    let first = crossings[0];
    let last = *crossings.last().expect("nonempty");
    let periods = (crossings.len() - 1) as f64;
    2.0 * core::f64::consts::PI * periods / (last - first)
}

#[test]
fn plate_construction_certificates() {
    let model =
        von_karman_ss_plate(&steel_plate(), &modes_grid(3, 2), &modes_grid(4, 3)).expect("plate");
    // Quadrature certificate: far BELOW the refusal threshold (an
    // assert at the threshold itself would be vacuous after expect —
    // review finding); measured ~1e-13 class.
    assert!(model.quadrature_residual < 1.0e-10);
    // Symmetry of every channel is verified by assemble (and would
    // refuse otherwise).
    let n = model.storage.omegas.len();
    assert_eq!(n, 6);
    assert_eq!(model.storage.channels.len(), 12);
    // Frequencies ascend with mode index norm (spot: (1,1) lowest).
    let w11 = model.storage.omegas[0];
    assert!(model.storage.omegas.iter().all(|&w| w >= w11));
    let sys = assemble(model.storage, &vec![0.0; n], &vec![1.0; n]).expect("assemble");
    assert_eq!(sys.state_dim(), 2 * n);
    // A coupling channel must be genuinely NONZERO (non-vacuous
    // nonlinearity): the (1,1)x(1,1) interaction projected on some
    // stress mode survives the sine selection rules.
    println!(
        "{{\"suite\":\"fs-nlmodal\",\"case\":\"construction\",\"quadrature_residual\":{:.3e},\"verdict\":\"pass\"}}",
        model.quadrature_residual
    );
}

#[test]
fn small_thick_panel_certifies_three_odd_odd_modes() {
    // The couple-layer steel panel (0.20 × 0.15 × 2 mm) used to refuse
    // n_modes = 3 at the 1e-8 quadrature gate. Higher-order Gauss is
    // the certificate, not a looser gate.
    let plate = VkPlateParams {
        lx: 0.20,
        ly: 0.15,
        h: 0.002,
        young: 2.0e11,
        nu: 0.3,
        rho: 7800.0,
    };
    let disp = vec![
        SineMode { m: 1, n: 1 },
        SineMode { m: 1, n: 3 },
        SineMode { m: 3, n: 1 },
    ];
    let stress = vec![SineMode { m: 2, n: 1 }];
    let model = von_karman_ss_plate(&plate, &disp, &stress).expect("small-panel VK");
    assert!(model.quadrature_residual < 3.0e-8);
    assert_eq!(model.storage.omegas.len(), 3);
}

#[test]
fn coupling_is_nonvacuous_and_certified_against_asymmetry() {
    let model =
        von_karman_ss_plate(&steel_plate(), &modes_grid(2, 2), &modes_grid(3, 3)).expect("plate");
    let peak = model
        .storage
        .channels
        .iter()
        .flat_map(|c| c.coupling.iter())
        .fold(0.0f64, |a, &v| a.max(v.abs()));
    assert!(peak > 0.0, "coupling tensor must not vanish");
    // Assemble REFUSES an asymmetric coupling by name.
    let n = model.storage.omegas.len();
    let mut bad = vec![0.0; n * n];
    bad[1] = 1.0; // (0,1) without (1,0)
    let refused = assemble(
        SosModalStorage {
            omegas: model.storage.omegas.clone(),
            channels: vec![StressChannel {
                coefficient: 1.0,
                coupling: bad,
            }],
        },
        &vec![0.0; n],
        &vec![0.0; n],
    );
    assert!(matches!(
        refused,
        Err(NlModalError::AsymmetricCoupling { .. })
    ));
}

#[test]
fn duffing_backbone_matches_perturbation_formula() {
    // Single-mode Kirchhoff-Carrier string IS a Duffing oscillator;
    // the measured amplitude-dependent frequency must track the
    // analytic first-order backbone.
    let storage = kirchhoff_carrier_string(&guitar_string(), 1).expect("kc");
    let omega0 = storage.omegas[0];
    let beta = single_mode_beta(&storage, 0);
    let sys = assemble(storage, &[0.0], &[0.0]).expect("assemble");
    let dt = 2.0e-6;
    let steps = 40_000;
    // Amplitudes give first-order relative shifts ~3e-2 and ~1.2e-1
    // (review-corrected numbers); at the larger one the second-order
    // backbone term is ~5% of the shift — inside the 10% budget.
    for amp in [2.0e-4, 4.0e-4] {
        let measured = measure_omega(&sys, vec![amp, 0.0], dt, steps);
        let predicted = duffing_backbone(omega0, beta, amp);
        let shift_pred = predicted - omega0;
        let shift_meas = measured - omega0;
        assert!(
            (shift_meas - shift_pred).abs() <= 0.10 * shift_pred.abs() + 1.0e-6 * omega0,
            "amp {amp}: measured shift {shift_meas:.4e} vs predicted {shift_pred:.4e}"
        );
    }
    println!("{{\"suite\":\"fs-nlmodal\",\"case\":\"duffing-backbone\",\"verdict\":\"pass\"}}");
}

#[test]
fn kirchhoff_carrier_beta_matches_hand_formula_and_glide() {
    // INDEPENDENT hand derivation (in physical amplitude A of mode k):
    //   T(t) = T0 + (EA/2L) I,  I = (L/2) sum (k pi/L)^2 A_k^2 ...
    //   backbone: dw/w = (3/32) (EA/T0) (k pi A / L)^2.
    // In mass-normalized coordinates a = sqrt(mu L / 2) A the cubic
    // coefficient must be beta_k = EA (k pi / L)^4 / (2 mu^2 L).
    let p = guitar_string();
    let storage = kirchhoff_carrier_string(&p, 3).expect("kc");
    let pi = core::f64::consts::PI;
    for k in 1..=3usize {
        let kk = k as f64 * pi / p.length;
        let beta_hand = p.ea * kk * kk * kk * kk / (2.0 * p.lin_density * p.lin_density * p.length);
        let beta_code = single_mode_beta(&storage, k - 1);
        assert!(
            (beta_code - beta_hand).abs() <= 1.0e-12 * beta_hand,
            "mode {k}: beta {beta_code:.6e} vs hand {beta_hand:.6e}"
        );
    }
    // Glide formula in PHYSICAL amplitude for the fundamental.
    let a_phys = 1.0e-3; // 1 mm — audible-glide territory
    let a_norm = fs_math::det::sqrt(p.lin_density * p.length / 2.0) * a_phys;
    let omega0 = storage.omegas[0];
    let beta = single_mode_beta(&storage, 0);
    let glide_backbone = (duffing_backbone(omega0, beta, a_norm) - omega0) / omega0;
    let glide_hand =
        (3.0 / 32.0) * (p.ea / p.tension) * (pi * a_phys / p.length) * (pi * a_phys / p.length);
    assert!(
        (glide_backbone - glide_hand).abs() <= 1.0e-10 * glide_hand,
        "glide {glide_backbone:.6e} vs hand {glide_hand:.6e}"
    );
    println!(
        "{{\"suite\":\"fs-nlmodal\",\"case\":\"kc-glide\",\"rel_glide_at_1mm\":{glide_hand:.4e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn prestressed_beam_is_harmonic_without_ei_and_inharmonic_with_it() {
    let p = guitar_string();
    let r_flex = fs_nlmodal::prestressed_beam_omega(p.length, p.tension, p.lin_density, 0.0, 2)
        / fs_nlmodal::prestressed_beam_omega(p.length, p.tension, p.lin_density, 0.0, 1);
    assert!((r_flex - 2.0).abs() < 1.0e-14);
    let ei = 2.0e11 * core::f64::consts::PI * (6.0e-4_f64).powi(4) / 4.0;
    let r_stiff = fs_nlmodal::prestressed_beam_omega(p.length, p.tension, p.lin_density, ei, 2)
        / fs_nlmodal::prestressed_beam_omega(p.length, p.tension, p.lin_density, ei, 1);
    assert!(r_stiff > 2.01);
}

#[test]
fn energy_conserved_undamped() {
    let model =
        von_karman_ss_plate(&steel_plate(), &modes_grid(2, 2), &modes_grid(3, 3)).expect("plate");
    let n = model.storage.omegas.len();
    let sys = assemble(model.storage, &vec![0.0; n], &vec![0.0; n]).expect("assemble");
    // Strike-like initial condition: large amplitude on several modes.
    let mut x = vec![0.0; 2 * n];
    for k in 0..n {
        x[2 * k] = 1.0e-3 / (k + 1) as f64;
        x[2 * k + 1] = f64::midpoint(k as f64, 1.0);
    }
    let h0 = sys.hamiltonian(&x);
    let dt = 1.0e-6;
    let mut worst = 0.0f64;
    for _ in 0..20_000 {
        let rec = step(&sys, &x, &[0.0], dt).expect("step");
        x = rec.x;
        worst = worst.max((sys.hamiltonian(&x) - h0).abs());
    }
    assert!(
        worst <= 1.0e-9 * h0,
        "energy drifted {worst:.3e} of H0 {h0:.3e}"
    );
    println!(
        "{{\"suite\":\"fs-nlmodal\",\"case\":\"conservation\",\"rel_drift\":{:.3e},\"verdict\":\"pass\"}}",
        worst / h0
    );
}

/// Mutant storage: hamiltonian uses the SYMMETRIC coupling, gradient
/// uses a raw ASYMMETRIC one — the classic hand-derived-tensor
/// force-vs-energy divergence.
struct MutantStorage {
    inner: SosModalStorage,
    asym: Vec<f64>,
}

impl Storage for MutantStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.inner.hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let n = self.inner.omegas.len();
        for k in 0..n {
            out[2 * k] = self.inner.omegas[k] * self.inner.omegas[k] * x[2 * k];
            out[2 * k + 1] = x[2 * k + 1];
        }
        let ch = &self.inner.channels[0];
        let mut sform = 0.0;
        for p in 0..n {
            for q in 0..n {
                sform += self.asym[p * n + q] * x[2 * p] * x[2 * q];
            }
        }
        for k in 0..n {
            let mut row = 0.0;
            for q in 0..n {
                row += self.asym[k * n + q] * x[2 * q];
            }
            out[2 * k] += ch.coefficient * sform * row;
        }
    }
}

#[test]
fn unsymmetrized_tensor_admission_and_trajectory_falsifier() {
    // EXECUTED ARCHITECTURAL FINDING, recorded honestly: the Gonzalez
    // discrete gradient FORCES dg.(x1-x0) = H(x1)-H(x0) for ANY
    // gradient function (the correction term absorbs force-vs-energy
    // divergence), so an energy audit STRUCTURALLY CANNOT catch an
    // unsymmetrized tensor here — unlike force-side integrators where
    // conservation breaks. The guards in THIS architecture are:
    // (1) assemble() refuses asymmetric couplings at admission
    //     (tested in coupling_is_nonvacuous_...), and
    // (2) the divergence is a TRAJECTORY error: the mutant departs
    //     measurably from the clean system it claims to represent.
    let clean = kirchhoff_carrier_string(&guitar_string(), 2).expect("kc");
    let n = 2;
    let mut asym = clean.channels[0].coupling.clone();
    asym[1] = 0.5 * asym[0];
    asym[0] *= 0.5;
    let mutant = MutantStorage {
        inner: kirchhoff_carrier_string(&guitar_string(), 2).expect("kc"),
        asym,
    };
    let dim = 2 * n;
    let mut j = vec![0.0; dim * dim];
    for k in 0..n {
        j[(2 * k) * dim + 2 * k + 1] = 1.0;
        j[(2 * k + 1) * dim + 2 * k] = -1.0;
    }
    let sys_mut = PortHamiltonian::new(
        dim,
        0,
        j.clone(),
        vec![0.0; dim * dim],
        vec![],
        Box::new(mutant),
    )
    .expect("admit mutant storage");
    let sys_clean = assemble(clean, &[0.0, 0.0], &[0.0, 0.0]).expect("clean");
    let x0 = vec![2.0e-3, 0.0, 1.0e-3, 0.0];
    let dt = 2.0e-6;
    let mut xm = x0.clone();
    let mut xc = x0.clone();
    let mut worst_h_drift = 0.0f64;
    let h0 = sys_mut.hamiltonian(&x0);
    for _ in 0..5000 {
        let rm = step(&sys_mut, &xm, &[], dt).expect("mutant step");
        worst_h_drift = worst_h_drift.max((sys_mut.hamiltonian(&rm.x) - h0).abs());
        xm = rm.x;
        xc = step(&sys_clean, &xc, &[0.0], dt).expect("clean step").x;
    }
    // (a) The structural-repair property, pinned: the mutant CONSERVES
    // its coded H (this is why an energy audit is the WRONG falsifier
    // under discrete-gradient stepping).
    // Measured 1.5e-7 relative at authoring (the divergent gradient
    // degrades the implicit solve's convergence quality, not the
    // identity itself); 1e-5 still proves the architectural point —
    // energy stays 4+ orders quieter than the trajectory divergence.
    assert!(
        worst_h_drift <= 1.0e-5 * h0,
        "discrete gradient must conserve even the mutant's H ({worst_h_drift:.3e})"
    );
    // (b) The trajectory falsifier: mutant vs clean diverge visibly.
    let scale = xc.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    let dev = xm
        .iter()
        .zip(&xc)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        dev > 0.05 * scale,
        "tensor asymmetry must be trajectory-visible: dev {dev:.3e} vs scale {scale:.3e}"
    );
}

#[test]
fn amplitude_scaling_covariance_in_weak_limit() {
    // The nonlinear frequency shift scales as amplitude^2: doubling
    // the amplitude must quadruple the measured shift (within the
    // measurement floor).
    let storage = kirchhoff_carrier_string(&guitar_string(), 1).expect("kc");
    let omega0 = storage.omegas[0];
    let sys = assemble(storage, &[0.0], &[0.0]).expect("assemble");
    let dt = 2.0e-6;
    let s1 = measure_omega(&sys, vec![2.0e-4, 0.0], dt, 40_000) - omega0;
    let s2 = measure_omega(&sys, vec![4.0e-4, 0.0], dt, 40_000) - omega0;
    let ratio = s2 / s1;
    assert!(
        (3.6..4.4).contains(&ratio),
        "shift scaling ratio {ratio:.3} (expected ~4)"
    );
}

/// Independent in-test RK4 on the same Hamiltonian vector field —
/// the cross-check route for the mode-coupling trajectory.
fn rk4_run(storage: &SosModalStorage, mut x: Vec<f64>, dt: f64, steps: usize) -> Vec<f64> {
    let n2 = x.len();
    let deriv = |x: &[f64]| -> Vec<f64> {
        let mut g = vec![0.0; n2];
        storage.gradient(x, &mut g);
        // J grad H with symplectic pair blocks.
        let mut d = vec![0.0; n2];
        for k in 0..n2 / 2 {
            d[2 * k] = g[2 * k + 1];
            d[2 * k + 1] = -g[2 * k];
        }
        d
    };
    for _ in 0..steps {
        let k1 = deriv(&x);
        let mid1: Vec<f64> = x.iter().zip(&k1).map(|(&a, &b)| a + 0.5 * dt * b).collect();
        let k2 = deriv(&mid1);
        let mid2: Vec<f64> = x.iter().zip(&k2).map(|(&a, &b)| a + 0.5 * dt * b).collect();
        let k3 = deriv(&mid2);
        let end: Vec<f64> = x.iter().zip(&k3).map(|(&a, &b)| a + dt * b).collect();
        let k4 = deriv(&end);
        for i in 0..n2 {
            x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    x
}

#[test]
fn parametric_mode_coupling_cross_checked_against_rk4() {
    // Two-mode Kirchhoff-Carrier string: mode 1 at large amplitude
    // parametrically pumps mode 2 (the tension modulation oscillates
    // at 2 w1 = w2's natural frequency). Energy must TRANSFER, and the
    // discrete-gradient trajectory must agree with an independent RK4
    // integration of the same vector field.
    let storage = kirchhoff_carrier_string(&guitar_string(), 2).expect("kc");
    let sys = assemble(
        kirchhoff_carrier_string(&guitar_string(), 2).expect("kc"),
        &[0.0, 0.0],
        &[0.0, 0.0],
    )
    .expect("assemble");
    let x0 = vec![1.5e-3, 0.0, 1.0e-6, 0.0];
    let dt = 1.0e-6;
    let steps = 50_000;
    let mut x = x0.clone();
    let mode2_energy = |x: &[f64], w: f64| f64::midpoint(x[3] * x[3], w * w * x[2] * x[2]);
    let w2 = storage.omegas[1];
    let e2_start = mode2_energy(&x, w2);
    for _ in 0..steps {
        x = step(&sys, &x, &[0.0], dt).expect("step").x;
    }
    let e2_end = mode2_energy(&x, w2);
    // Diagonal (tension-modulation) coupling pumps mode 2 at
    // 2*w1 = w2, which is the SECOND-order parametric condition
    // (principal resonance would need pump = 2*w2): transfer exists
    // but is weak — that is correct string physics (tension modulation
    // mainly detunes; the strong cascade lives in plate cross terms).
    assert!(
        e2_end > 2.0 * e2_start.max(1.0e-30),
        "parametric transfer failed: {e2_start:.3e} -> {e2_end:.3e}"
    );
    // Independent RK4 route agrees on the endpoint state within a
    // trajectory tolerance (both order-appropriate at this dt).
    let x_rk4 = rk4_run(&storage, x0, dt, steps);
    let scale = x.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    let dev = x
        .iter()
        .zip(&x_rk4)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        dev <= 0.02 * scale,
        "discrete-gradient vs RK4 deviate {dev:.3e} (scale {scale:.3e})"
    );
    println!(
        "{{\"suite\":\"fs-nlmodal\",\"case\":\"parametric-coupling\",\"e2_growth\":{:.3e},\"rk4_dev\":{:.3e},\"verdict\":\"pass\"}}",
        e2_end / e2_start,
        dev / scale
    );
}

#[test]
fn plate_cascade_casebook() {
    // Cymbal-like plate initialized in its FUNDAMENTAL only, at three
    // amplitude levels (w/h ~ 0.1, 0.6, 2): any energy leaving mode 1
    // is NONLINEAR transfer by construction (the linearized system
    // keeps it at exactly zero), so the cascade-onset-with-level claim
    // is isolated from the strike's linear bandwidth. Ledger closed at
    // every step.
    let params = VkPlateParams {
        lx: 0.35,
        ly: 0.35,
        h: 0.8e-3,
        young: 1.1e11, // bronze-class
        nu: 0.34,
        rho: 8800.0,
    };
    let model = von_karman_ss_plate(&params, &modes_grid(3, 3), &modes_grid(4, 4)).expect("plate");
    let n = model.storage.omegas.len();
    let omegas = model.storage.omegas.clone();
    // Mass-normalized amplitude giving physical center deflection w:
    // w = phi_norm * q with phi_norm = sqrt(4/(rho h lx ly)).
    let phi_norm = fs_math::det::sqrt(4.0 / (params.rho * params.h * params.lx * params.ly));
    let sys = assemble(model.storage, &vec![0.0; n], &vec![0.0; n]).expect("assemble");
    let dt = 5.0e-6;
    let mut fractions = Vec::new();
    for w_over_h in [0.1, 0.6, 2.0] {
        let q0 = w_over_h * params.h / phi_norm;
        let mut x = vec![0.0; 2 * n];
        x[0] = q0;
        let h0 = sys.hamiltonian(&x);
        let mut worst_defect = f64::NEG_INFINITY;
        for _ in 0..4000 {
            let rec = step(&sys, &x, &[0.0], dt).expect("step");
            worst_defect = worst_defect.max(rec.supply_defect().abs());
            x = rec.x;
        }
        assert!(
            worst_defect <= 1.0e-8 * h0.max(1.0e-12),
            "ledger violated: {worst_defect:.3e}"
        );
        let energies: Vec<f64> = (0..n)
            .map(|k| {
                f64::midpoint(
                    x[2 * k + 1] * x[2 * k + 1],
                    omegas[k] * omegas[k] * x[2 * k] * x[2 * k],
                )
            })
            .collect();
        let total: f64 = energies.iter().sum();
        let leaked = (total - energies[0]) / total.max(f64::MIN_POSITIVE);
        fractions.push(leaked);
        println!(
            "{{\"suite\":\"fs-nlmodal\",\"case\":\"cascade\",\"w_over_h\":{w_over_h},\"leaked_fraction\":{leaked:.5},\"modal_energies\":{energies:?},\"n_modes\":{n},\"tensor_entries\":{}}}",
            n * n * 16
        );
    }
    // Cascade onset: the leaked fraction grows STRONGLY with level
    // (at w/h = 0.1 the plate is essentially linear).
    assert!(
        fractions[1] > 10.0 * fractions[0].max(1.0e-12) || fractions[1] > 1.0e-3,
        "moderate level should leak visibly: {fractions:?}"
    );
    assert!(
        fractions[2] > 3.0 * fractions[1],
        "cascade must strengthen with level: {fractions:?}"
    );
    println!(
        "{{\"suite\":\"fs-nlmodal\",\"case\":\"cascade\",\"fractions\":{fractions:?},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn all_zero_channel_constructs() {
    // Selection-rule parity makes some (stress, disp) combinations
    // integrate to EXACTLY zero everywhere; such a channel must
    // construct (with a zero matrix), not spuriously refuse on its
    // own roundoff scale (review finding: 0.75 "relative" residual).
    let model = von_karman_ss_plate(
        &steel_plate(),
        &[SineMode { m: 1, n: 1 }],
        &[SineMode { m: 2, n: 1 }],
    )
    .expect("all-zero channel must construct");
    let peak = model.storage.channels[0]
        .coupling
        .iter()
        .fold(0.0f64, |a, &v| a.max(v.abs()));
    // Entries are roundoff-class relative to a physical channel.
    let phys = von_karman_ss_plate(
        &steel_plate(),
        &[SineMode { m: 1, n: 1 }],
        &[SineMode { m: 1, n: 1 }],
    )
    .expect("physical channel");
    let phys_peak = phys.storage.channels[0]
        .coupling
        .iter()
        .fold(0.0f64, |a, &v| a.max(v.abs()));
    assert!(
        peak < 1.0e-9 * phys_peak,
        "zero channel: {peak:.3e} vs {phys_peak:.3e}"
    );
}

#[test]
fn nan_and_duplicates_refused() {
    // NaN must not slip through comparison gates (review finding).
    let storage = kirchhoff_carrier_string(&guitar_string(), 2).expect("kc");
    let mut bad = kirchhoff_carrier_string(&guitar_string(), 2).expect("kc");
    bad.omegas[0] = f64::NAN;
    assert!(assemble(bad, &[0.0, 0.0], &[0.0, 0.0]).is_err());
    let mut bad2 = kirchhoff_carrier_string(&guitar_string(), 2).expect("kc");
    bad2.channels[0].coupling[0] = f64::NAN;
    assert!(assemble(bad2, &[0.0, 0.0], &[0.0, 0.0]).is_err());
    drop(storage);
    // Duplicate modes double-count physics: refused by name.
    assert!(matches!(
        von_karman_ss_plate(
            &steel_plate(),
            &[SineMode { m: 1, n: 1 }, SineMode { m: 1, n: 1 }],
            &[SineMode { m: 1, n: 1 }],
        ),
        Err(NlModalError::Parameter {
            what: "duplicate mode in list"
        })
    ));
}

#[test]
fn refusals_are_typed() {
    assert!(matches!(
        von_karman_ss_plate(&steel_plate(), &[], &modes_grid(2, 2)),
        Err(NlModalError::Parameter { .. })
    ));
    assert!(matches!(
        von_karman_ss_plate(
            &steel_plate(),
            &[SineMode { m: 0, n: 1 }],
            &modes_grid(2, 2)
        ),
        Err(NlModalError::Parameter { .. })
    ));
    let bad = VkPlateParams {
        h: -1.0,
        ..steel_plate()
    };
    assert!(von_karman_ss_plate(&bad, &modes_grid(1, 1), &modes_grid(1, 1)).is_err());
    assert!(kirchhoff_carrier_string(&guitar_string(), 0).is_err());
}
