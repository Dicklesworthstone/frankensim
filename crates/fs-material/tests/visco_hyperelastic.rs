//! M0.3-W1b: finite-strain Maxwell branches on a hyperelastic equilibrium law. Every oracle is written out here and
//! shares no code with the law under test: the generalized-Maxwell time series of `fs_material::visco`, a scalar RK4
//! of the same evolution law on the axisymmetric principal stretch, the Ogden `isochoric_energy` gradient of the
//! equilibrium module (with the viscous map applied by a Cholesky factor), and central finite differences. The
//! integrator is backward Euler in the logarithmic elastic stretches, formal order 1; where a test needs more than
//! that it extrapolates two step sizes and says so. No parameter is fitted; the branch sets are fixtures.

use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::tensor::rotation;
use fs_material::visco::GeneralizedMaxwell;
use fs_material::visco_hyper::{
    BranchEnergy, FiniteStrainBranch, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic, ViscoState,
};
use fs_material::MaterialError;
use fs_math::det;
use std::time::Instant;

const I: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];

fn nh_branch(mu: f64, tau: f64) -> FiniteStrainBranch {
    FiniteStrainBranch {
        energy: BranchEnergy::NeoHookean { mu_pa: mu },
        tau_s: tau,
    }
}

fn nh_equilibrium(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 6.0).unwrap()
}

fn nh_law(mu_inf: f64, lambda: f64, branches: &[(f64, f64)]) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        nh_equilibrium(mu_inf, lambda),
        branches.iter().map(|&(m, t)| nh_branch(m, t)).collect(),
    )
    .unwrap()
}

fn uniaxial(l: f64) -> [f64; 9] {
    let s = 1.0 / det::sqrt(l);
    [l, 0., 0., 0., s, 0., 0., 0., s]
}

fn flat(m: [[f64; 3]; 3]) -> [f64; 9] {
    core::array::from_fn(|i| m[i / 3][i % 3])
}

fn matmul(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| (0..3).map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3]).sum())
}

fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64; 9]) -> f64 {
    det::sqrt(dot(a, a))
}

fn det3(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn stream(name: &str) -> fs_rand::Stream {
    let seed = name.bytes().fold(0xcbf29ce484222325u64, |s, b| {
        (s ^ u64::from(b)).wrapping_mul(0x100000001b3)
    });
    fs_rand::StreamKey {
        seed,
        kernel: 27,
        tile: 0,
    }
    .stream()
}

fn random_f(rng: &mut fs_rand::Stream, amplitude: f64) -> [f64; 9] {
    loop {
        let f: [f64; 9] = core::array::from_fn(|i| I[i] + amplitude * (2. * rng.next_f64() - 1.));
        if det3(&f) > 0.3 {
            return f;
        }
    }
}

fn close(a: f64, b: f64, scale: f64, tol: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tol * scale,
        "{a:e} vs {b:e}, scale {scale:e}, tolerance {tol:e}"
    );
}

fn close_all(a: &[f64; 9], b: &[f64; 9], tol: f64, floor: f64) {
    let scale = norm(a).max(norm(b)).max(floor);
    for (x, y) in a.iter().zip(b) {
        close(*x, *y, scale, tol);
    }
}

// Axial minus lateral Cauchy stress of an isochoric uniaxial state, per unit Hencky strain amplitude.
fn deviatoric_axial(p: &[f64; 9], f: &[f64; 9], e0: f64) -> f64 {
    (p[0] * f[0] - p[4] * f[4]) / e0
}

// ---------------------------------------------------------------------------------------------------------------
// G1: the small-strain limit is the generalized Maxwell law.
// ---------------------------------------------------------------------------------------------------------------

fn hencky_uniaxial(e: f64) -> [f64; 9] {
    [
        det::exp(e),
        0.,
        0.,
        0.,
        det::exp(-0.5 * e),
        0.,
        0.,
        0.,
        det::exp(-0.5 * e),
    ]
}

#[test]
fn g1_small_strain_limit_is_generalized_maxwell() {
    let t0 = Instant::now();
    let g_inf = 1.0e5;
    let branches = [(2.0e4, 1.0e-3), (5.0e4, 1.0e-2), (1.0e5, 1.0e-1)];
    let law = nh_law(g_inf, 1.0e9, &branches);
    // Incompressible uniaxial: E = 3 G, branch by branch.
    let gm = GeneralizedMaxwell::new(
        3.0 * g_inf,
        branches.iter().map(|&(g, t)| (3.0 * g, t)).collect(),
    )
    .unwrap();
    let e0 = 1.0e-6;

    // Relaxation modulus after a strain step: backward Euler at dt and dt/2, extrapolated to second order.
    let dt = 1.0e-6;
    let samples = [1_000usize, 10_000, 100_000];
    let run_relaxation = |h: f64, factor: usize| -> Vec<f64> {
        let f = hencky_uniaxial(e0);
        let (mut s, _) = law.step(&law.initial_state(), &I, 0.0).unwrap();
        let (jump, p) = law.step(&s, &f, 0.0).unwrap();
        s = jump;
        let mut out = vec![deviatoric_axial(&p, &f, e0)];
        let mut next = 0;
        for k in 1..=samples[2] * factor {
            let (ns, p) = law.step(&s, &f, h).unwrap();
            s = ns;
            if next < samples.len() && k == samples[next] * factor {
                out.push(deviatoric_axial(&p, &f, e0));
                next += 1;
            }
        }
        out
    };
    let coarse = run_relaxation(dt, 1);
    let fine = run_relaxation(0.5 * dt, 2);
    let mut gm_state = gm.state();
    let mut reference = vec![gm.step(&mut gm_state, e0, 0.0) / e0];
    let mut next = 0;
    for k in 1..=samples[2] {
        let stress = gm.step(&mut gm_state, e0, dt);
        if next < samples.len() && k == samples[next] {
            reference.push(stress / e0);
            next += 1;
        }
    }
    let mut worst: f64 = 0.0;
    let mut first_order: f64 = 0.0;
    for i in 0..reference.len() {
        let extrapolated = 2.0 * fine[i] - coarse[i];
        let err = (extrapolated - reference[i]).abs() / reference[i].abs();
        worst = worst.max(err);
        first_order = first_order.max((coarse[i] - reference[i]).abs() / reference[i].abs());
        assert!(
            err <= 1e-5,
            "relaxation modulus sample {i}: {extrapolated:e} vs {:e}, relative {err:e}",
            reference[i]
        );
    }

    // Complex modulus at five frequencies across the branch band: the same sinusoidal history through both laws from
    // rest; the in-phase and quadrature stress over the second period, with identical quadrature.
    let omegas = [10.0, 30.0, 100.0, 300.0, 1000.0];
    let mut worst_complex: f64 = 0.0;
    for &omega in &omegas {
        let period = 2.0 * core::f64::consts::PI / omega;
        let moduli = |h_steps: usize| -> (f64, f64, f64, f64) {
            let h = period / h_steps as f64;
            let mut s = law.initial_state();
            let mut g = gm.state();
            let (mut a_ip, mut a_q, mut g_ip, mut g_q) = (0.0, 0.0, 0.0, 0.0);
            for k in 1..=2 * h_steps {
                let t = k as f64 * h;
                let e = e0 * det::sin(omega * t);
                let f = hencky_uniaxial(e);
                let (ns, p) = law.step(&s, &f, h).unwrap();
                s = ns;
                let ours = deviatoric_axial(&p, &f, e0) * e0;
                let theirs = gm.step(&mut g, e, h);
                if k > h_steps {
                    a_ip += ours * det::sin(omega * t) * h;
                    a_q += ours * det::cos(omega * t) * h;
                    g_ip += theirs * det::sin(omega * t) * h;
                    g_q += theirs * det::cos(omega * t) * h;
                }
            }
            let c = 2.0 / (e0 * period);
            (c * a_ip, c * a_q, c * g_ip, c * g_q)
        };
        let (a1, q1, g1, gq1) = moduli(2_000);
        let (a2, q2, g2, gq2) = moduli(4_000);
        // Richardson: the law's backward Euler is first order; the generalized Maxwell step is exact for the sampled
        // piecewise-linear strain, so its own two resolutions differ only by the quadrature.
        let ours = (2.0 * a2 - a1, 2.0 * q2 - q1);
        let theirs = (2.0 * g2 - g1, 2.0 * gq2 - gq1);
        let amplitude = det::sqrt(theirs.0 * theirs.0 + theirs.1 * theirs.1);
        let err = det::sqrt((ours.0 - theirs.0).powi(2) + (ours.1 - theirs.1).powi(2)) / amplitude;
        worst_complex = worst_complex.max(err);
        assert!(
            err <= 1e-5,
            "complex modulus at omega {omega}: ours {ours:?}, generalized Maxwell {theirs:?}, relative {err:e}"
        );
    }
    println!(
        "{{\"bead\":\"M0.3-W1b\",\"case\":\"g1_small_strain_limit_is_generalized_maxwell\",\"inputs\":{{\"branches\":3,\"equilibrium\":\"neo-hookean\",\"strain_amplitude\":{e0:e},\"omegas\":{omegas:?}}},\"provenance\":\"test fixture, no source\",\"rungs\":[{dt:e},{:e}],\"observed\":{{\"relaxation_extrapolated_worst\":{worst:e},\"relaxation_first_order_worst\":{first_order:e},\"complex_worst\":{worst_complex:e}}},\"tolerance\":1e-5,\"verdict\":\"pass\",\"elapsed_s\":{:.2}}}",
        0.5 * dt,
        t0.elapsed().as_secs_f64()
    );
}

// ---------------------------------------------------------------------------------------------------------------
// G1: step-and-hold relaxation against an independently integrated reference of the same evolution law.
// ---------------------------------------------------------------------------------------------------------------

// Axisymmetric principal log stretch a of a neo-Hookean branch (elastic principal stretches e^a, e^-a/2, e^-a/2):
// da/dt = -(e^{2a} - e^{-a}) / (3 tau), from the principal Kirchhoff stress (2 mu/3)(e^{2a} - e^{-a}) and
// D_v = tau_dev / (2 G tau). RK4, written without the law.
fn axial_log_stretch(a0: f64, tau: f64, t_end: f64, h: f64) -> f64 {
    let rate = |a: f64| -(det::exp(2. * a) - det::exp(-a)) / (3. * tau);
    let mut a = a0;
    let steps = (t_end / h).round() as usize;
    for _ in 0..steps {
        let k1 = rate(a);
        let k2 = rate(a + 0.5 * h * k1);
        let k3 = rate(a + 0.5 * h * k2);
        let k4 = rate(a + h * k3);
        a += h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
    }
    a
}

fn reference_axial_piola(mu_inf: f64, branches: &[(f64, f64)], lambda: f64, t: f64, h: f64) -> f64 {
    let a0 = det::ln(lambda);
    let mut p = mu_inf * (lambda - 1.0 / lambda);
    for &(mu, tau) in branches {
        let a = axial_log_stretch(a0, tau, t, h);
        p += (2. * mu / 3.) * (det::exp(2. * a) - det::exp(-a)) / lambda;
    }
    p
}

#[test]
fn g1_step_stretch_relaxation_reference() {
    let t0 = Instant::now();
    let mu_inf = 1.0e5;
    let branches = [(3.0e4, 0.02), (8.0e4, 0.2)];
    let law = nh_law(mu_inf, 1.0e9, &branches);
    let tau_min = 0.02;
    let times = [0.1 * tau_min, tau_min, 5.0 * tau_min];
    let mut cases = Vec::new();
    let mut worst_order_gap: f64 = 0.0;
    for &lambda in &[1.5, 2.0, 3.0] {
        let f = uniaxial(lambda);
        let refs: Vec<f64> = times
            .iter()
            .map(|&t| reference_axial_piola(mu_inf, &branches, lambda, t, tau_min / 1.0e4))
            .collect();
        // errors[rung][sample]
        let mut errors = [[0.0f64; 3]; 3];
        for rung in 0..3 {
            let h = tau_min / (10.0 * f64::from(1u32 << rung));
            let (mut s, _) = law.step(&law.initial_state(), &I, 0.0).unwrap();
            s = law.step(&s, &f, 0.0).unwrap().0;
            let mut t = 0.0;
            for (i, &target) in times.iter().enumerate() {
                let steps = ((target - t) / h).round() as usize;
                let mut p = [0.0; 9];
                for _ in 0..steps {
                    let (ns, pk) = law.step(&s, &f, h).unwrap();
                    s = ns;
                    p = pk;
                }
                t = target;
                errors[rung][i] = (p[0] - refs[i]).abs() / refs[i].abs();
            }
        }
        for i in 0..3 {
            let p12 = (errors[0][i] / errors[1][i]).log2();
            let p23 = (errors[1][i] / errors[2][i]).log2();
            worst_order_gap = worst_order_gap.max((p12 - 1.0).abs()).max((p23 - 1.0).abs());
            assert!(
                (p12 - 1.0).abs() <= 0.2 && (p23 - 1.0).abs() <= 0.2,
                "lambda {lambda}, t = {:e}: observed orders {p12:.3}, {p23:.3} (errors {:e}, {:e}, {:e})",
                times[i],
                errors[0][i],
                errors[1][i],
                errors[2][i]
            );
        }
        // The long-time stress is the equilibrium stress.
        let (mut s, _) = law.step(&law.initial_state(), &I, 0.0).unwrap();
        s = law.step(&s, &f, 0.0).unwrap().0;
        for _ in 0..4 {
            s = law.step(&s, &f, 0.1).unwrap().0;
        }
        let (_, p_long) = law.step(&s, &f, 1.0e12).unwrap();
        let p_eq = law.equilibrium.piola(&f).unwrap();
        close_all(&p_long, &p_eq, 1e-9, 1.0);
        cases.push((lambda, errors));
    }
    println!(
        "{{\"bead\":\"M0.3-W1b\",\"case\":\"g1_step_stretch_relaxation_reference\",\"inputs\":{{\"branches\":{},\"equilibrium\":\"neo-hookean\",\"stretch\":[1.5,2.0,3.0],\"times_s\":{times:?}}},\"provenance\":\"test fixture, no source\",\"rungs\":[{:e},{:e},{:e}],\"observed\":{{\"errors\":{cases:?}}},\"order\":{{\"formal\":1,\"worst_gap\":{worst_order_gap:.3}}},\"tolerance\":0.2,\"verdict\":\"pass\",\"elapsed_s\":{:.2}}}",
        branches.len(),
        tau_min / 10.0,
        tau_min / 20.0,
        tau_min / 40.0,
        t0.elapsed().as_secs_f64()
    );
}

// ---------------------------------------------------------------------------------------------------------------
// G1: a cyclic history. The ledger closes at roundoff; the stress-power quadrature closes at the integrator's order.
// ---------------------------------------------------------------------------------------------------------------

#[test]
fn g1_cyclic_loop_area_matches_ledger_at_integrator_order() {
    let law = nh_law(1.0e5, 1.0e8, &[(5.0e4, 0.1), (8.0e4, 1.0)]);
    let omega = 10.0;
    let period = 2.0 * core::f64::consts::PI / omega;
    let cycles = 6usize;
    let mut defects = Vec::new();
    for &per in &[100usize, 200, 400] {
        let h = period / per as f64;
        let mut s = law.initial_state();
        let mut p_prev = [0.0; 9];
        let mut f_prev = I;
        let (mut quad, mut quad_last) = (0.0, 0.0);
        let (mut stored_last0, mut diss_last0) = (0.0, 0.0);
        for k in 1..=cycles * per {
            if k == (cycles - 1) * per + 1 {
                stored_last0 = s.stored;
                diss_last0 = s.dissipated;
            }
            let f = uniaxial(1.0 + 0.4 * det::sin(omega * k as f64 * h));
            let (ns, p) = law.step(&s, &f, h).unwrap();
            let df: [f64; 9] = core::array::from_fn(|i| f[i] - f_prev[i]);
            let pm: [f64; 9] = core::array::from_fn(|i| 0.5 * (p[i] + p_prev[i]));
            let w = dot(&pm, &df);
            quad += w;
            if k > (cycles - 1) * per {
                quad_last += w;
            }
            // The ledger closes per step: loading work = change of stored energy + dissipation.
            let dw = ns.work - s.work;
            let ds = ns.stored - s.stored;
            let dd = ns.dissipated - s.dissipated;
            assert!(dd >= 0.0, "dissipation {dd:e}");
            close(dw, ds + dd, dw.abs().max(ds.abs()).max(dd).max(1e-6), 1e-10);
            s = ns;
            p_prev = p;
            f_prev = f;
        }
        assert!(s.dissipated > 0.0);
        close(s.work, s.stored + s.dissipated, s.dissipated, 1e-10);
        let total = (quad - (s.stored + s.dissipated)).abs() / s.dissipated;
        let last = (quad_last - ((s.stored - stored_last0) + (s.dissipated - diss_last0))).abs()
            / (s.dissipated - diss_last0);
        defects.push((total, last));
    }
    // Backward Euler and the trapezoid rule differ at first order in the step: the defect halves with it.
    for i in 0..2 {
        let order = (defects[i].0 / defects[i + 1].0).log2();
        assert!(
            (0.75..=1.25).contains(&order),
            "cumulative quadrature defects {defects:?}, order {order:.3}"
        );
    }
    assert!(defects[2].1 < 2.0e-2, "last-cycle defects {defects:?}");
    println!(
        "{{\"bead\":\"M0.3-W1b\",\"case\":\"g1_cyclic_loop_area_matches_ledger_at_integrator_order\",\"observed\":{{\"quadrature_defects\":{defects:?}}},\"verdict\":\"pass\"}}"
    );
}

// ---------------------------------------------------------------------------------------------------------------
// G0: dissipation and the ledger over random multiaxial histories.
// ---------------------------------------------------------------------------------------------------------------

fn random_law() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        nh_equilibrium(1.0e5, 3.0e5),
        vec![
            nh_branch(4.0e4, 0.05),
            FiniteStrainBranch {
                energy: BranchEnergy::Ogden {
                    terms: vec![
                        OgdenTerm {
                            mu_pa: 3.0e4,
                            alpha: 3.0,
                        },
                        OgdenTerm {
                            mu_pa: -1.0e4,
                            alpha: -2.0,
                        },
                    ],
                },
                tau_s: 0.5,
            },
        ],
    )
    .unwrap()
}

#[test]
fn g0_dissipation_nonnegative_random_multiaxial() {
    let t0 = Instant::now();
    let law = random_law();
    let mut rng = stream("g0_dissipation");
    let (mut min_dissipation, mut worst_defect): (f64, f64) = (f64::INFINITY, 0.0);
    let histories = 10_000;
    for _ in 0..histories {
        let mut s = law.initial_state();
        for _ in 0..4 {
            let f = random_f(&mut rng, 0.25);
            let dt = 0.05 * det::exp(det::ln(10.0) * (5.0 * rng.next_f64() - 3.0));
            let (ns, _) = law.step(&s, &f, dt).unwrap();
            let dw = ns.work - s.work;
            let ds = ns.stored - s.stored;
            let dd = ns.dissipated - s.dissipated;
            assert!(dd >= 0.0, "negative dissipation {dd:e} at dt {dt:e}");
            min_dissipation = min_dissipation.min(dd);
            let scale = dw.abs().max(ds.abs()).max(dd).max(1.0);
            worst_defect = worst_defect.max((dw - ds - dd).abs() / scale);
            assert!(
                (dw - ds - dd).abs() <= 1e-10 * scale,
                "ledger defect {:e} (work {dw:e}, stored {ds:e}, dissipated {dd:e})",
                dw - ds - dd
            );
            s = ns;
        }
    }
    println!(
        "{{\"bead\":\"M0.3-W1b\",\"case\":\"g0_dissipation_nonnegative_random_multiaxial\",\"inputs\":{{\"histories\":{histories},\"steps\":4,\"branches\":\"neo-hookean + two-term ogden\",\"seed\":\"fnv(g0_dissipation)\"}},\"provenance\":\"test fixture, no source\",\"observed\":{{\"min_step_dissipation\":{min_dissipation:e},\"ledger_defect_worst\":{worst_defect:e}}},\"tolerance\":1e-10,\"verdict\":\"pass\",\"elapsed_s\":{:.2}}}",
        t0.elapsed().as_secs_f64()
    );
}

// ---------------------------------------------------------------------------------------------------------------
// G0: objectivity, rest state, instantaneous and equilibrium limits.
// ---------------------------------------------------------------------------------------------------------------

#[test]
fn g0_objectivity_and_rest_state() {
    let law = random_law();
    // Rest: P = 0, no flow, empty ledger.
    let s0 = law.initial_state();
    let (s1, p) = law.step(&s0, &I, 0.3).unwrap();
    assert!(norm(&p) <= 1e-9, "rest stress {:e}", norm(&p));
    for (a, b) in s1
        .inverse_viscous_cauchy_green
        .iter()
        .zip(&s0.inverse_viscous_cauchy_green)
    {
        for (x, y) in a.iter().zip(b) {
            close(*x, *y, 1.0, 1e-14);
        }
    }
    assert!(s1.stored.abs() <= 1e-12 && s1.dissipated <= 1e-18);

    let mut rng = stream("g0_objectivity");
    for _ in 0..100 {
        let fs: Vec<[f64; 9]> = (0..6).map(|_| random_f(&mut rng, 0.2)).collect();
        let dts: Vec<f64> = (0..6).map(|_| 0.02 + 0.2 * rng.next_f64()).collect();
        let mut axis: [f64; 3] = core::array::from_fn(|_| 2. * rng.next_f64() - 1.);
        let n = det::sqrt(axis.iter().map(|x| x * x).sum());
        for x in &mut axis {
            *x /= n;
        }
        let q = flat(rotation(axis, 6. * rng.next_f64() - 3.));
        let (mut a, mut b) = (law.initial_state(), law.initial_state());
        for k in 0..6 {
            let fb = if k >= 3 { matmul(&q, &fs[k]) } else { fs[k] };
            let (na, pa) = law.step(&a, &fs[k], dts[k]).unwrap();
            let (nb, pb) = law.step(&b, &fb, dts[k]).unwrap();
            a = na;
            b = nb;
            // The referential branch state is unchanged by a superposed rotation; the stress rotates with it.
            for (ca, cb) in a
                .inverse_viscous_cauchy_green
                .iter()
                .zip(&b.inverse_viscous_cauchy_green)
            {
                let scale = ca.iter().map(|v| v.abs()).fold(1.0, f64::max);
                for (x, y) in ca.iter().zip(cb) {
                    close(*x, *y, scale, 1e-12);
                }
            }
            let expected = if k >= 3 { matmul(&q, &pa) } else { pa };
            close_all(&pb, &expected, 1e-12, 1.0e3);
        }
    }
}

// Oracle branch stress for a relaxed or evolved state: the Ogden isochoric-energy gradient of the equilibrium module
// at F_e = F L, with C_v^{-1} = L L^T (Cholesky), mapped back by P = P_e L^T.
fn oracle_branch_piola(branch: &FiniteStrainBranch, cv_inv: &[f64; 6], f: &[f64; 9]) -> [f64; 9] {
    let terms = match &branch.energy {
        BranchEnergy::NeoHookean { mu_pa } => vec![OgdenTerm {
            mu_pa: *mu_pa,
            alpha: 2.0,
        }],
        BranchEnergy::Ogden { terms } => terms.clone(),
    };
    let ogden = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms,
            volumetric: Volumetric::Quadratic { kappa: 1.0e6 },
        },
        1.0e6,
    )
    .unwrap();
    let c = [
        [cv_inv[0], cv_inv[3], cv_inv[4]],
        [cv_inv[3], cv_inv[1], cv_inv[5]],
        [cv_inv[4], cv_inv[5], cv_inv[2]],
    ];
    let mut l = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let mut s = c[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            l[i][j] = if i == j { det::sqrt(s) } else { s / l[j][j] };
        }
    }
    let lf = flat(l);
    let fe = matmul(f, &lf);
    let (_, pe) = fs_ad::gradient::<9>(fe, |x| ogden.isochoric_energy(&x));
    let lt: [f64; 9] = core::array::from_fn(|ij| lf[3 * (ij % 3) + ij / 3]);
    matmul(&pe, &lt)
}

#[test]
fn g0_instantaneous_and_equilibrium_limits() {
    let law = random_law();
    let mut rng = stream("g0_limits");
    for _ in 0..50 {
        let f = random_f(&mut rng, 0.25);
        let s = law.initial_state();
        // dt = 0: the equilibrium response plus the instantaneous response of every branch.
        let (_, p0) = law.step(&s, &f, 0.0).unwrap();
        let mut expected = law.equilibrium.piola(&f).unwrap();
        for (b, cv) in law.branches.iter().zip(&s.inverse_viscous_cauchy_green) {
            let pb = oracle_branch_piola(b, cv, &f);
            for (e, v) in expected.iter_mut().zip(pb) {
                *e += v;
            }
        }
        close_all(&p0, &expected, 1e-10, 1.0);
        // dt -> infinity: the equilibrium response only.
        let (_, pinf) = law.step(&s, &f, 1.0e15).unwrap();
        close_all(&pinf, &law.equilibrium.piola(&f).unwrap(), 1e-10, 1.0);
        // The stress at the end of a finite step is the frozen stress at the relaxed state.
        let (relaxed, p_step) = law.step(&s, &f, 0.07).unwrap();
        let (_, p_frozen) = law.step(&relaxed, &f, 0.0).unwrap();
        close_all(&p_step, &p_frozen, 1e-10, 1.0);
        let mut oracle = law.equilibrium.piola(&f).unwrap();
        for (b, cv) in law.branches.iter().zip(&relaxed.inverse_viscous_cauchy_green) {
            let pb = oracle_branch_piola(b, cv, &f);
            for (e, v) in oracle.iter_mut().zip(pb) {
                *e += v;
            }
        }
        close_all(&p_step, &oracle, 1e-10, 1.0);
    }
    // The frozen stored energy has the frozen stress as its gradient (integrated along a straight path).
    let s = law.step(&law.initial_state(), &random_f(&mut rng, 0.1), 0.05).unwrap().0;
    let f0 = s.deformation_gradient;
    let f1 = random_f(&mut rng, 0.1);
    let gl_x = [-0.906_179_845_938_664, -0.538_469_310_105_683_1, 0.0, 0.538_469_310_105_683_1, 0.906_179_845_938_664];
    let gl_w = [0.236_926_885_056_189_1, 0.478_628_670_499_366_5, 0.568_888_888_888_888_9, 0.478_628_670_499_366_5, 0.236_926_885_056_189_1];
    let df: [f64; 9] = core::array::from_fn(|i| f1[i] - f0[i]);
    let mut work = 0.0;
    for (x, w) in gl_x.iter().zip(gl_w.iter()) {
        let t = 0.5 * (1.0 + x);
        let ft: [f64; 9] = core::array::from_fn(|i| f0[i] + t * df[i]);
        let (_, p) = law.step(&s, &ft, 0.0).unwrap();
        work += 0.5 * w * dot(&p, &df);
    }
    let de = law.free_energy(&s, &f1).unwrap() - law.free_energy(&s, &f0).unwrap();
    close(work, de, de.abs().max(1.0), 1e-9);
}

// ---------------------------------------------------------------------------------------------------------------
// G0: the algorithmic tangent is the exact derivative of the step map.
// ---------------------------------------------------------------------------------------------------------------

fn fd_tangent(law: &ViscoHyperelastic, s: &ViscoState, f: &[f64; 9], dt: f64, h: f64) -> [[f64; 9]; 9] {
    let mut a = [[0.0; 9]; 9];
    for col in 0..9 {
        let (mut fp, mut fm) = (*f, *f);
        fp[col] += h;
        fm[col] -= h;
        let pp = law.step(s, &fp, dt).unwrap().1;
        let pm = law.step(s, &fm, dt).unwrap().1;
        for row in 0..9 {
            a[row][col] = (pp[row] - pm[row]) / (2.0 * h);
        }
    }
    a
}

#[test]
fn g0_algorithmic_tangent_matches_fd() {
    let law = random_law();
    let mut rng = stream("g0_tangent");
    let mut worst: f64 = 0.0;
    let mut cases: Vec<(String, [f64; 9], ViscoState, f64)> = Vec::new();
    for i in 0..8 {
        let mut s = law.initial_state();
        for _ in 0..3 {
            s = law.step(&s, &random_f(&mut rng, 0.2), 0.04).unwrap().0;
        }
        cases.push((format!("generic {i}"), random_f(&mut rng, 0.2), s, 0.05));
    }
    // Isotropic trial stretches: exactly equal roots, and roots a hair apart (the swap-path integral).
    let rest = law.initial_state();
    cases.push(("rest".to_string(), I, rest.clone(), 0.05));
    cases.push(("uniform stretch".to_string(), core::array::from_fn(|i| 1.2 * I[i]), rest.clone(), 0.05));
    let near: [f64; 9] = core::array::from_fn(|i| 1.1 * I[i] + 1e-7 * (2. * rng.next_f64() - 1.));
    cases.push(("near-isotropic".to_string(), near, rest.clone(), 0.05));
    cases.push(("frozen (dt = 0)".to_string(), random_f(&mut rng, 0.2), rest, 0.0));
    for (name, f, s, dt) in &cases {
        let a = law.algorithmic_tangent(s, f, *dt).unwrap();
        let fd = fd_tangent(&law, s, f, *dt, 1e-6);
        let scale = a.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        let mut err: f64 = 0.0;
        let mut asym: f64 = 0.0;
        for r in 0..9 {
            for c in 0..9 {
                err = err.max((a[r][c] - fd[r][c]).abs() / scale);
                asym = asym.max((a[r][c] - a[c][r]).abs() / scale);
            }
        }
        worst = worst.max(err);
        assert!(err <= 1e-6, "{name}: tangent versus central differences, relative {err:e}");
        assert!(asym <= 1e-9, "{name}: tangent asymmetry {asym:e}");
    }
    assert!(law.admissibility().tangent_symmetric);
    println!("{{\"bead\":\"M0.3-W1b\",\"case\":\"g0_algorithmic_tangent_matches_fd\",\"observed\":{{\"worst_relative\":{worst:e}}},\"tolerance\":1e-6,\"verdict\":\"pass\"}}");
}

// ---------------------------------------------------------------------------------------------------------------
// Refusals.
// ---------------------------------------------------------------------------------------------------------------

#[test]
fn branch_cap_and_cap_plus_one() {
    let eq = nh_equilibrium(1.0e5, 1.0e8);
    let make = |n: usize| {
        ViscoHyperelastic::new(eq.clone(), (0..n).map(|i| nh_branch(1.0e4, 0.01 * (1 + i) as f64)).collect())
    };
    assert!(make(MAX_FINITE_STRAIN_BRANCHES).is_ok());
    assert!(matches!(make(MAX_FINITE_STRAIN_BRANCHES + 1), Err(MaterialError::Parameters { .. })));
}

#[test]
fn refuses_nonpositive_tau_or_negative_branch_modulus() {
    let eq = nh_equilibrium(1.0e5, 1.0e8);
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ViscoHyperelastic::new(eq.clone(), vec![nh_branch(1.0e4, tau)]),
            Err(MaterialError::Parameters { .. })
        ));
    }
    for mu in [0.0, -1.0e4, f64::NAN] {
        assert!(matches!(
            ViscoHyperelastic::new(eq.clone(), vec![nh_branch(mu, 0.1)]),
            Err(MaterialError::Parameters { .. })
        ));
    }
    let bad_ogden = FiniteStrainBranch {
        energy: BranchEnergy::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 1.0e4,
                alpha: -2.0,
            }],
        },
        tau_s: 0.1,
    };
    assert!(matches!(
        ViscoHyperelastic::new(eq, vec![bad_ogden]),
        Err(MaterialError::Parameters { .. })
    ));
}

#[test]
fn refuses_det_f_nonpositive_and_bad_steps() {
    let law = nh_law(1.0e5, 1.0e8, &[(1.0e4, 0.1)]);
    let s = law.initial_state();
    let flipped = [-1., 0., 0., 0., 1., 0., 0., 0., 1.];
    let flat = [1., 0., 0., 0., 1., 0., 0., 0., 0.];
    for f in [flipped, flat] {
        assert!(matches!(law.step(&s, &f, 0.1), Err(MaterialError::State { .. })));
        assert!(matches!(law.algorithmic_tangent(&s, &f, 0.1), Err(MaterialError::State { .. })));
    }
    // A negative or non-finite step is a typed error, not a panic.
    for dt in [-1.0e-9, f64::NAN, f64::INFINITY] {
        assert!(matches!(law.step(&s, &I, dt), Err(MaterialError::State { .. })));
    }
    // A state that does not belong to the law.
    let other = nh_law(1.0e5, 1.0e8, &[(1.0e4, 0.1), (2.0e4, 1.0)]).initial_state();
    assert!(matches!(law.step(&other, &I, 0.1), Err(MaterialError::State { .. })));
}

#[test]
fn local_newton_nonconvergence_is_typed_refusal() {
    let law = nh_law(1.0e5, 1.0e8, &[(1.0e4, 0.1)]);
    let f = uniaxial(2.0);
    assert!(law.step(&law.initial_state(), &f, 0.05).is_ok());
    // A tolerance no double-precision residual can reach: the solve refuses; it never returns its best iterate.
    let strict = law.clone().with_local_newton(1.0e-300, 64).unwrap();
    assert!(matches!(
        strict.step(&strict.initial_state(), &f, 0.05),
        Err(MaterialError::State { .. })
    ));
    // Too few iterations for a large first step.
    let capped = law.clone().with_local_newton(1.0e-13, 1).unwrap();
    assert!(matches!(
        capped.step(&capped.initial_state(), &uniaxial(3.0), 10.0),
        Err(MaterialError::State { .. })
    ));
    assert!(law.clone().with_local_newton(0.0, 8).is_err());
    assert!(law.clone().with_local_newton(1.0e-12, 0).is_err());
}
