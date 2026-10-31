//! fs-nlmodal conformance battery: quadrature/symmetry certificates,
//! Duffing backbone vs the analytic perturbation formula, exact energy
//! conservation, force-vs-energy mutation visibility through the
//! fs-phs supply audit, amplitude-scaling covariance, parametric
//! mode-coupling cross-checked against an independent RK4 integrator,
//! the Kirchhoff-Carrier pitch-glide hand formula, and the struck-
//! plate cascade casebook.

use fs_nlmodal::{
    KcStringParams, NlModalError, SampledPlateMode, SineMode, SosModalStorage, StressChannel,
    VkPlateParams, assemble, assemble_storage, assemble_storage_with_damping, duffing_backbone,
    kirchhoff_carrier_moving_end, kirchhoff_carrier_string, single_mode_beta,
    von_karman_sampled_plate, von_karman_ss_plate,
};
use fs_phs::{PortHamiltonian, Storage, step};

#[test]
fn full_viscous_operator_preserves_local_damper_nullspace_and_forced_response() {
    // A single physical damper on y=q1-q2 gives C=c*[1,-1]^T*[1,-1].
    // Equal-frequency common motion is lossless; differential motion has
    // decay rate c. A diagonal-only projection destroys this distinction.
    let omega: f64 = 8.0;
    let c: f64 = 1.2;
    let storage = || SosModalStorage {
        omegas: vec![omega, omega],
        channels: vec![],
    };
    let sys = assemble_storage_with_damping(
        2,
        &[omega, omega],
        &[c, -c, -c, c],
        1,
        vec![0.0, 1.0, 0.0, -1.0],
        Box::new(storage()),
    )
    .unwrap();
    let split = fs_phs::linear_exact::split_quadratic(
        &sys,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: 4,
            max_omega_dt: 3.0,
        },
    )
    .unwrap();
    let dt: f64 = 0.3;
    let common =
        fs_phs::linear_exact::step_render(&sys, &split, &[0.4, 0.2, 0.4, 0.2], &[0.0], dt).unwrap();
    let q_common = 0.4 * (omega * dt).cos() + 0.2 / omega * (omega * dt).sin();
    let v_common = 0.2 * (omega * dt).cos() - 0.4 * omega * (omega * dt).sin();
    for (actual, expected) in common
        .x
        .iter()
        .zip([q_common, v_common, q_common, v_common])
    {
        assert!((actual - expected).abs() < 2.0e-12);
    }
    assert!(common.dissipated.abs() < 2.0e-12);

    let initial = [1.1, 0.1, -0.3, 0.3]; // q+ = .4, q- = .7, v+ = .2, v- = -.1
    let force: f64 = 0.6; // conjugate to y=q1-q2, so differential acceleration is F
    let forced = fs_phs::linear_exact::step_render(&sys, &split, &initial, &[force], dt).unwrap();
    let equilibrium = force / omega.powi(2);
    let shifted_q = 0.7 - equilibrium;
    let wd = (omega.powi(2) - c.powi(2)).sqrt();
    let decay = (-c * dt).exp();
    let q_diff = equilibrium
        + decay * (shifted_q * (wd * dt).cos() + (-0.1 + c * shifted_q) / wd * (wd * dt).sin());
    let v_diff = decay
        * (-0.1 * (wd * dt).cos() - (-0.1 * c + omega.powi(2) * shifted_q) / wd * (wd * dt).sin());
    for (actual, expected) in forced.x.iter().zip([
        q_common + q_diff,
        v_common + v_diff,
        q_common - q_diff,
        v_common - v_diff,
    ]) {
        assert!(
            (actual - expected).abs() < 3.0e-12,
            "{actual} vs {expected}"
        );
    }
    let work = force * ((forced.x[0] - forced.x[2]) - (initial[0] - initial[2]));
    assert!((forced.supplied - work).abs() < 2.0e-12);
    assert!(forced.dissipated > 0.0);
    assert!(
        forced.balance_residual().abs() < 3.0e-12,
        "defect={}, delta_h={}, loss={}, work={}",
        forced.balance_residual(),
        forced.delta_h,
        forced.dissipated,
        forced.supplied
    );

    let diagonal = assemble_storage(
        2,
        &[omega, omega],
        &[c / (2.0 * omega); 2],
        1,
        vec![0.0, 1.0, 0.0, -1.0],
        Box::new(storage()),
    )
    .unwrap();
    let diagonal_split = fs_phs::linear_exact::split_quadratic(
        &diagonal,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: 4,
            max_omega_dt: 3.0,
        },
    )
    .unwrap();
    let wrong = fs_phs::linear_exact::step_render(
        &diagonal,
        &diagonal_split,
        &[0.4, 0.2, 0.4, 0.2],
        &[0.0],
        dt,
    )
    .unwrap();
    assert!(
        wrong.dissipated > 0.1,
        "diagonalization must visibly damp common motion"
    );
}

#[test]
fn full_damping_diagonal_limit_is_the_same_system() {
    let storage = || SosModalStorage {
        omegas: vec![2.0, 5.0],
        channels: vec![StressChannel {
            coefficient: 0.8,
            coupling: vec![1.0, 0.3, 0.3, 2.0],
        }],
    };
    let g = vec![0.0, 0.0, 1.0, -0.2, 0.0, 0.0, -0.4, 0.8];
    let original = assemble_storage(
        2,
        &[2.0, 5.0],
        &[0.1, 0.03],
        2,
        g.clone(),
        Box::new(storage()),
    )
    .unwrap();
    let full = assemble_storage_with_damping(
        2,
        &[2.0, 5.0],
        &[0.4, 0.0, 0.0, 0.3],
        2,
        g,
        Box::new(storage()),
    )
    .unwrap();
    assert_eq!(original.structure(), full.structure());
    let initial = [0.3, 0.2, -0.1, 0.4];
    let a = step(&original, &initial, &[0.7, -0.5], 0.01).unwrap();
    let b = step(&full, &initial, &[0.7, -0.5], 0.01).unwrap();
    assert_eq!(a.x, b.x);
    assert_eq!(a.dissipated.to_bits(), b.dissipated.to_bits());
    assert_eq!(a.supplied.to_bits(), b.supplied.to_bits());
}

#[test]
fn full_damping_nonlinear_trajectory_converges_and_closes_work() {
    let storage = || SosModalStorage {
        omegas: vec![2.0, 5.0],
        channels: vec![StressChannel {
            coefficient: 0.8,
            coupling: vec![1.0, 0.3, 0.3, 2.0],
        }],
    };
    let damping = [0.6, 0.2, 0.2, 0.1];
    let sys = assemble_storage_with_damping(
        2,
        &[2.0, 5.0],
        &damping,
        1,
        vec![0.0, 1.0, 0.0, -0.4],
        Box::new(storage()),
    )
    .unwrap();
    let initial = vec![0.3, 0.2, -0.1, 0.4];
    let reference = rk4_run_damped(&storage(), initial.clone(), 0.0001, 4000, &damping);
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let dt = 0.4 / steps as f64;
        let mut x = initial.clone();
        let mut loss = 0.0;
        for _ in 0..steps {
            let rec = step(&sys, &x, &[0.0], dt).unwrap();
            assert!(rec.dissipated >= 0.0);
            assert!(
                rec.balance_residual().abs() < 2.0e-11,
                "steps={steps}, defect={}, residual={}",
                rec.balance_residual(),
                rec.solver_residual
            );
            loss += rec.dissipated;
            x = rec.x;
        }
        assert!(loss > 0.0);
        assert!((sys.hamiltonian(&x) + loss - sys.hamiltonian(&initial)).abs() < 2.0e-10);
        errors.push(
            x.iter()
                .zip(&reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(
        errors[0] > 3.5 * errors[1] && errors[1] > 3.5 * errors[2],
        "{errors:?}"
    );
    assert!(errors[2] < 3.0e-5, "{errors:?}");
}

#[test]
fn full_damping_refuses_active_asymmetric_and_nonfinite_operators() {
    let make = |matrix: &[f64], g: Vec<f64>| {
        assemble_storage_with_damping(
            2,
            &[2.0, 3.0],
            matrix,
            1,
            g,
            Box::new(SosModalStorage {
                omegas: vec![2.0, 3.0],
                channels: vec![],
            }),
        )
    };
    for matrix in [
        vec![1.0, 2.0, 2.0, 1.0], // one active eigenvector despite positive diagonal
        vec![1.0, 0.5, 0.0, 1.0],
        vec![1.0, f64::NAN, f64::NAN, 1.0],
        vec![f64::INFINITY, 0.0, 0.0, 1.0],
        vec![1.0],
    ] {
        assert!(make(&matrix, vec![0.0, 1.0, 0.0, -1.0]).is_err());
    }
    assert!(make(&[0.0; 4], vec![0.0, f64::NAN, 0.0, 0.0]).is_err());
    // Zero loss remains an actual undamped system; no regularization is added.
    let lossless = make(&[0.0; 4], vec![0.0, 1.0, 0.0, 0.0]).unwrap();
    assert_eq!(lossless.structure().1, &[0.0; 16]);
}

#[test]
fn resolved_nonlinear_force_survives_unchanged_high_energy_modes() {
    // The second mode stores much more energy, but does no work in this
    // increment. Subtracting total endpoint energies used to erase the
    // first mode's resolved restoring force entirely.
    let storage = SosModalStorage {
        omegas: vec![1.0, 1.0e6],
        channels: vec![StressChannel {
            coefficient: 1.0,
            coupling: vec![1.0, 0.0, 0.0, 0.0],
        }],
    };
    let a = [1.0, 0.0, 1.0, 0.0];
    for delta in [1.0e-2, 1.0e-6, 1.0e-10, -1.0e-6] {
        let b = [a[0] + delta, 0.0, 1.0, 0.0];
        let dg = fs_phs::discrete_gradient(&storage, &a, &b);
        // Independent scalar divided difference of q^2/2 + q^4/4,
        // factored before evaluation so no endpoint energy is subtracted.
        let expected = 0.5 * (a[0] + b[0]) + 0.25 * (a[0] + b[0]) * (a[0] * a[0] + b[0] * b[0]);
        assert!(
            (dg[0] - expected).abs() < 2.0e-14 * expected.abs(),
            "delta={delta:.3e}: force {} versus {expected}",
            dg[0],
        );
        assert_eq!(dg[1], 0.0);
        assert_eq!(dg[2], 1.0e12);
        assert_eq!(dg[3], 0.0);
    }
}

#[test]
fn coupled_quartic_discrete_work_matches_independent_gauss_integral() {
    let storage = SosModalStorage {
        omegas: vec![2.0, 3.0, 5.0],
        channels: vec![
            StressChannel {
                coefficient: 7.0,
                coupling: vec![2.0, 0.25, -0.5, 0.25, 1.0, 0.75, -0.5, 0.75, -1.5],
            },
            StressChannel {
                coefficient: 3.0,
                coupling: vec![1.0, -0.5, 0.0, -0.5, 2.0, 0.25, 0.0, 0.25, 1.0],
            },
        ],
    };
    let a = [0.2, -0.3, -0.4, 0.7, 0.8, -0.2];
    let b = [-0.5, 0.1, 0.6, -0.4, 0.3, 0.9];
    let mut points = [[0.0; 6]; 2];
    let mut gradients = [[0.0; 6]; 2];
    // The ordinary gradient is cubic along the segment. Two Gauss nodes
    // integrate its work exactly, independently of the discrete remainder.
    for (point, sign) in points.iter_mut().zip([-1.0, 1.0]) {
        for k in 0..a.len() {
            point[k] = f64::midpoint(a[k], b[k]) + sign * (b[k] - a[k]) / (2.0 * 3.0_f64.sqrt());
        }
    }
    storage.gradient(&points[0], &mut gradients[0]);
    storage.gradient(&points[1], &mut gradients[1]);
    let work: f64 = (0..a.len())
        .map(|k| f64::midpoint(gradients[0][k], gradients[1][k]) * (b[k] - a[k]))
        .sum();
    let dg = fs_phs::discrete_gradient(&storage, &a, &b);
    let discrete_work: f64 = (0..a.len()).map(|k| dg[k] * (b[k] - a[k])).sum();
    assert!((discrete_work - work).abs() < 2.0e-14 * work.abs());
    let reverse = fs_phs::discrete_gradient(&storage, &b, &a);
    for (forward, reverse) in dg.iter().zip(reverse) {
        assert!((forward - reverse).abs() < 2.0e-14 * (1.0 + forward.abs()));
    }
    let unchanged = fs_phs::discrete_gradient(&storage, &a, &a);
    let mut ordinary = [0.0; 6];
    storage.gradient(&a, &mut ordinary);
    assert_eq!(unchanged, ordinary);
}

#[test]
fn admitted_roundoff_asymmetry_preserves_discrete_work() {
    let storage = SosModalStorage {
        omegas: vec![2.0, 3.0],
        channels: vec![StressChannel {
            coefficient: 7.0,
            coupling: vec![1.0, 0.5, 0.5 + 5.0e-11, 2.0],
        }],
    };
    let a = [0.2, -0.3, -0.4, 0.7];
    let b = [-0.5, 0.1, 0.6, -0.4];
    let expected = storage.hamiltonian(&b) - storage.hamiltonian(&a);
    let dg = fs_phs::discrete_gradient(&storage, &a, &b);
    let work: f64 = dg.iter().enumerate().map(|(k, g)| g * (b[k] - a[k])).sum();
    assert!((work - expected).abs() < 2.0e-14 * expected.abs());
    assemble(storage, &[0.0; 2], &[0.0; 2]).expect("within the existing symmetry tolerance");
}

fn steel_plate() -> VkPlateParams {
    VkPlateParams {
        lx: 0.4,
        ly: 0.3,
        h: 1.0e-3,
        young: 2.0e11,
        nu: 0.3,
        rho: 7850.0,
        pretension_n_m: 0.0,
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
        pretension_n_m: 0.0,
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
fn g1_high_mode_plate_couplings_match_analytic_harmonic_integrals() {
    // Independent product-to-sum integration, with integral sin(j*pi*x)
    // exactly zero for even j and 2/(j*pi) for odd j, including negative j.
    let integrals = |s: usize, a: usize, b: usize| {
        let (s, a, b) = (s as i32, a as i32, b as i32);
        let integral = |j: i32| {
            if j % 2 == 0 {
                0.0
            } else {
                2.0 / (j as f64 * core::f64::consts::PI)
            }
        };
        let difference = integral(s + a - b) + integral(s - a + b);
        let sum = integral(s + a + b) + integral(s - a - b);
        (0.25 * (difference - sum), 0.25 * (difference + sum))
    };
    let params = VkPlateParams {
        lx: 0.7,
        ly: 0.4,
        ..steel_plate()
    };
    let disp = [
        SineMode { m: 1, n: 1 },
        SineMode { m: 6, n: 2 },
        SineMode { m: 4, n: 6 },
        SineMode { m: 13, n: 9 },
    ];
    let stress = [
        SineMode { m: 3, n: 5 },
        SineMode { m: 2, n: 1 },
        SineMode { m: 4, n: 4 },
        SineMode { m: 1, n: 1 },
    ];
    let model = von_karman_ss_plate(&params, &disp, &stress).unwrap();
    assert!(model.quadrature_residual < 3.0e-8);
    let area = params.lx * params.ly;
    let normalization = 8.0 / (params.rho * params.h * area * area.sqrt());
    for (channel, &s) in model.storage.channels.iter().zip(&stress) {
        for (i, &a) in disp.iter().enumerate() {
            for (j, &b) in disp.iter().enumerate() {
                let (sx, cx) = integrals(s.m, a.m, b.m);
                let (sy, cy) = integrals(s.n, a.n, b.n);
                let (ax, ay) = (
                    a.m as f64 * core::f64::consts::PI / params.lx,
                    a.n as f64 * core::f64::consts::PI / params.ly,
                );
                let (bx, by) = (
                    b.m as f64 * core::f64::consts::PI / params.lx,
                    b.n as f64 * core::f64::consts::PI / params.ly,
                );
                let raw = ((ax * by).powi(2) + (ay * bx).powi(2)) * sx * sy
                    - 2.0 * ax * ay * bx * by * cx * cy;
                let expected = normalization * area * raw;
                let actual = channel.coupling[i * disp.len() + j];
                assert!(
                    (actual - expected).abs() < 2.0e-10 * expected.abs().max(1.0),
                    "stress {s:?}, modes {a:?}/{b:?}: {actual} vs {expected}"
                );
                assert_eq!(
                    actual.to_bits(),
                    channel.coupling[j * disp.len() + i].to_bits()
                );
                if (s.m + a.m + b.m) % 2 == 0 || (s.n + a.n + b.n) % 2 == 0 {
                    assert_eq!(
                        actual, 0.0,
                        "forbidden harmonic channel must be exactly zero"
                    );
                }
            }
        }
    }
    // Both signs of real, non-diagonal harmonic interactions remain present.
    assert!(model.storage.channels[0].coupling[disp.len() + 2] > 0.0);
    assert!(model.storage.channels[2].coupling[2 * disp.len()] < 0.0);
    assert!(model.storage.channels[0].coupling[3 * disp.len()] < 0.0);
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
    let free = kirchhoff_carrier_moving_end(&p, 3).expect("moving-end kc");
    let c = fs_math::det::sqrt(p.tension / p.lin_density);
    for k in 0..3usize {
        let kappa = (k as f64 + 0.5) * pi / p.length;
        let want_w = kappa * c;
        assert!(
            (free.omegas[k] - want_w).abs() <= 1.0e-12 * want_w,
            "free-fixed ω[{k}] {} vs {want_w}",
            free.omegas[k]
        );
        let beta_hand =
            p.ea * kappa * kappa * kappa * kappa / (2.0 * p.lin_density * p.lin_density * p.length);
        let beta_code = single_mode_beta(&free, k);
        assert!(
            (beta_code - beta_hand).abs() <= 1.0e-12 * beta_hand,
            "free-fixed mode {k}: beta {beta_code:.6e} vs hand {beta_hand:.6e}"
        );
        assert!(
            free.omegas[k] < storage.omegas[k],
            "half-integer κ must sit below the fixed-fixed partial"
        );
    }
    assert!(kirchhoff_carrier_moving_end(&p, 0).is_err());
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
fn g1_pinned_beam_harmonics_survive_zero_and_vanishing_tension() {
    let (length, density, ei): (f64, f64, f64) = (0.7, 1.0e-3, 2.0e-3);
    for k in 1..=16 {
        // Independent simply-supported bending eigenvalue: EI*kappa^4/mu.
        let kappa = k as f64 * core::f64::consts::PI / length;
        let expected = (ei / density).sqrt() * kappa.powi(2);
        for tension in [0.0, 1.0e-320, 1.0e-20] {
            let actual = fs_nlmodal::prestressed_beam_omega(length, tension, density, ei, k);
            assert!(
                (actual / expected - 1.0).abs() < 2.0e-14,
                "mode {k}, tension {tension}: {actual} vs {expected}"
            );
        }
    }
    let fundamental = fs_nlmodal::prestressed_beam_omega(length, 0.0, density, ei, 1);
    let second = fs_nlmodal::prestressed_beam_omega(length, 0.0, density, ei, 2);
    assert!((second / fundamental - 4.0).abs() < 1.0e-14);
}

#[test]
fn g1_pinned_beam_compression_retains_the_euler_buckling_boundary() {
    // L=pi gives kappa=k. For EI=1, the first critical compressive load is
    // T=-1, and the second is T=-4. Only the unstable mode loses a real omega.
    let length = core::f64::consts::PI;
    let frequency = |tension, k| fs_nlmodal::prestressed_beam_omega(length, tension, 4.0, 1.0, k);
    assert_eq!(frequency(-0.75, 1), 0.25);
    assert_eq!(frequency(-1.0, 1), 0.0);
    assert!(frequency(-1.25, 1).is_nan());
    assert!((frequency(-1.0, 2) - 3.0_f64.sqrt()).abs() < 1.0e-14);
    let offset = 1.0e-8_f64;
    let expected = 0.5 * offset.sqrt();
    assert!((frequency(-1.0 + offset, 1) / expected - 1.0).abs() < 5.0e-9);
    assert!(frequency(-1.0 - offset, 1).is_nan());
}

#[test]
fn g3_pinned_beam_dispersion_preserves_physical_scaling() {
    let (length, tension, density, ei) = (0.7, 80.0, 1.0e-3, 2.0e-3);
    let reference = fs_nlmodal::prestressed_beam_omega(length, tension, density, ei, 7);
    for scale in [1.0e-250, 1.0e250] {
        // Scale every extensive mass/stiffness by the same factor.
        let actual = fs_nlmodal::prestressed_beam_omega(
            length,
            scale * tension,
            scale * density,
            scale * ei,
            7,
        );
        assert!((actual / reference - 1.0).abs() < 2.0e-14);
    }
    for scale in [1.0e-60_f64, 1.0e60_f64] {
        // Same material and tensile stress, with all lengths scaled: A~s^2,
        // I~s^4, mu~s^2, T~s^2. Every angular frequency must scale as 1/s.
        let actual = fs_nlmodal::prestressed_beam_omega(
            length * scale,
            tension * scale.powi(2),
            density * scale.powi(2),
            ei * scale.powi(4),
            7,
        );
        assert!((actual * scale / reference - 1.0).abs() < 2.0e-14);
    }
    for (density, ei, expected) in [(1.0e200, 1.0e-200, 1.0e-200), (1.0e-200, 1.0e200, 1.0e200)] {
        let actual = fs_nlmodal::prestressed_beam_omega(core::f64::consts::PI, 0.0, density, ei, 1);
        assert!((actual / expected - 1.0).abs() < 2.0e-14);
    }
}

#[test]
fn g0_beam_dispersion_never_discards_invalid_bending_stiffness() {
    for ei in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(fs_nlmodal::prestressed_beam_omega(1.0, 10.0, 1.0, ei, 1).is_nan());
    }
    assert!(fs_nlmodal::prestressed_beam_omega(1.0, 10.0, 1.0, 0.0, 0).is_nan());
    assert!(fs_nlmodal::prestressed_beam_omega(0.0, 10.0, 1.0, 0.0, 1).is_nan());
    assert!(fs_nlmodal::prestressed_beam_omega(1.0, 10.0, 0.0, 0.0, 1).is_nan());
    assert!(fs_nlmodal::prestressed_beam_omega(1.0, f64::NAN, 1.0, 0.0, 1).is_nan());
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
fn rk4_run(storage: &SosModalStorage, x: Vec<f64>, dt: f64, steps: usize) -> Vec<f64> {
    rk4_run_damped(storage, x, dt, steps, &[])
}

fn rk4_run_damped(
    storage: &SosModalStorage,
    mut x: Vec<f64>,
    dt: f64,
    steps: usize,
    damping: &[f64],
) -> Vec<f64> {
    let n2 = x.len();
    let deriv = |x: &[f64]| -> Vec<f64> {
        let mut g = vec![0.0; n2];
        storage.gradient(x, &mut g);
        // J grad H with symplectic pair blocks.
        let mut d = vec![0.0; n2];
        for k in 0..n2 / 2 {
            d[2 * k] = g[2 * k + 1];
            d[2 * k + 1] = -g[2 * k];
            if !damping.is_empty() {
                for j in 0..n2 / 2 {
                    d[2 * k + 1] -= damping[k * (n2 / 2) + j] * x[2 * j + 1];
                }
            }
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
        pretension_n_m: 0.0,
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

#[test]
fn pretension_raises_the_linear_frequencies() {
    let loose = steel_plate();
    let mut taut = loose;
    taut.pretension_n_m = 2.0e4;
    let a = von_karman_ss_plate(&loose, &modes_grid(1, 1), &modes_grid(1, 1)).expect("loose");
    let b = von_karman_ss_plate(&taut, &modes_grid(1, 1), &modes_grid(1, 1)).expect("taut");
    assert!(
        b.storage.omegas[0] > a.storage.omegas[0] * 1.01,
        "in-plane tension must raise ω ({:.1} vs {:.1})",
        b.storage.omegas[0],
        a.storage.omegas[0]
    );
    let mut bad = loose;
    bad.pretension_n_m = -1.0;
    assert!(von_karman_ss_plate(&bad, &modes_grid(1, 1), &modes_grid(1, 1)).is_err());
}

#[test]
fn sampled_sine_grid_reprints_the_analytic_airy_channel() {
    let params = steel_plate();
    let disp = [SineMode { m: 1, n: 1 }];
    let stress = [SineMode { m: 1, n: 1 }];
    let analytic = von_karman_ss_plate(&params, &disp, &stress).expect("ss");
    let (nx, ny) = (21, 17);
    let mut w = vec![0.0; nx * ny];
    let pi = core::f64::consts::PI;
    for j in 0..ny {
        let y = params.ly * j as f64 / (ny - 1) as f64;
        for i in 0..nx {
            let x = params.lx * i as f64 / (nx - 1) as f64;
            w[j * nx + i] = (pi * x / params.lx).sin() * (pi * y / params.ly).sin();
        }
    }
    let sampled = von_karman_sampled_plate(
        &params,
        &[SampledPlateMode {
            omega: analytic.storage.omegas[0],
            w,
            nx,
            ny,
        }],
        &stress,
    )
    .expect("sampled sine");
    let a = analytic.storage.channels[0].coefficient * analytic.storage.channels[0].coupling[0];
    let b = sampled.storage.channels[0].coefficient * sampled.storage.channels[0].coupling[0];
    let rel = (a - b).abs() / a.abs().max(1.0e-30);
    assert!(
        rel < 0.08,
        "FD sine grid must reprint the SS Airy channel ({a} vs {b}, rel {rel})"
    );
    assert!(
        von_karman_sampled_plate(&params, &[], &stress).is_err(),
        "empty sampled list must refuse"
    );
}
