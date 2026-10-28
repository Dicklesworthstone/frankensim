//! Generalized-α convergence regression across the ρ∞ range (kept from
//! the tfz.12 probe that diagnosed the battery's period-point metric
//! blindness: at t = 2π, cos′ = 0, so a q-only error measures phase
//! error QUADRATICALLY and fakes order ≈ 4; max(q, v) error is the
//! honest metric and shows clean order 2 at every ρ∞).

use fs_time::{GeneralizedAlpha, galpha_step};

fn err_at(rho: f64, h_nom: f64, t_end: f64) -> f64 {
    let steps = (t_end / h_nom).round() as usize;
    let hh = t_end / steps as f64;
    let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, hh, rho);
    let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0f64]);
    let f = vec![0.0f64];
    for _ in 0..steps {
        galpha_step(&ga, &mut q, &mut v, &mut a, &f);
    }
    let (qe, ve) = (t_end.cos(), -t_end.sin());
    (q[0] - qe).abs().max((v[0] - ve).abs())
}

#[test]
fn galpha_order_two_across_rho_range() {
    for &rho in &[0.0f64, 0.3, 0.5, 0.8, 1.0] {
        for &t in &[1.7f64, std::f64::consts::TAU] {
            let (e1, e2, e3) = (
                err_at(rho, 0.04, t),
                err_at(rho, 0.02, t),
                err_at(rho, 0.01, t),
            );
            let (o1, o2) = ((e1 / e2).log2(), (e2 / e3).log2());
            println!(
                "{{\"suite\":\"fs-time\",\"case\":\"galpha-order-sweep\",\"verdict\":\"info\",\
                 \"detail\":\"rho={rho} T={t:.3} errs={e1:.3e}/{e2:.3e}/{e3:.3e} \
                 orders={o1:.2},{o2:.2}\"}}"
            );
            assert!(
                (o1 - 2.0).abs() < 0.15 && (o2 - 2.0).abs() < 0.15,
                "generalized-alpha order off at rho={rho}, T={t}: {o1:.2}, {o2:.2}"
            );
        }
    }
}

#[test]
fn galpha_reaches_the_static_response_under_a_constant_load() {
    // Regression: the external load enters the equilibrium RHS with COEFFICIENT
    // 1 (f_next is already the load at t+(1-af)h). A stray (1-af) scaling drove
    // the steady state to (1-af)*f/K instead of f/K — wrong for any nonzero load
    // with rho_inf > 0 (load effectively halved at rho_inf=1), and invisible to
    // the free-vibration (f=0) tests. 1-DOF damped SDOF settles to the static
    // K*q = f, i.e. q -> f/K for EVERY rho_inf.
    let (mm, cc, kk, load) = (1.0, 0.5, 2.0, 3.0);
    let target = load / kk;
    for &rho in &[0.0, 0.3, 0.5, 0.8, 1.0] {
        let ga = GeneralizedAlpha::new(&[mm], &[cc], &[kk], 1, 0.05, rho);
        let (mut q, mut v, mut a) = (vec![0.0], vec![0.0], vec![0.0]);
        let f = vec![load];
        for _ in 0..4000 {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
        }
        assert!(
            (q[0] - target).abs() < 1e-6,
            "rho_inf={rho}: static response {} != f/K = {target}",
            q[0]
        );
    }
}

#[test]
fn dense_galpha_retains_vibration_below_position_resolution() {
    // G0: cos(h) rounds to one, while sin(h) and acceleration are resolved.
    // Recovering a_next from q_next-q_old either loses the vibration or
    // amplifies position roundoff by 1/h². The tiny-step case also exercises
    // the acceleration matrix without an overflowing reciprocal h².
    for h in [1e-8_f64, 1e-160] {
        for rho in [0.0, 0.3, 0.5, 0.8, 1.0] {
            let method = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, h, rho);
            let (mut q, mut v, mut a) = ([1.0], [0.0], [-1.0]);
            galpha_step(&method, &mut q, &mut v, &mut a, &[0.0]);
            assert_eq!(q[0], 1.0, "rho={rho}, h={h}: unresolved position moved");
            assert!(
                (v[0] + h.sin()).abs() <= 8.0 * f64::EPSILON * h,
                "rho={rho}, h={h}: lost resolved velocity {}",
                v[0]
            );
            assert!(
                (a[0] + h.cos()).abs() <= 8.0 * f64::EPSILON,
                "rho={rho}, h={h}: lost resolved acceleration {}",
                a[0]
            );
        }
    }
}

#[test]
fn dense_galpha_rigid_acceleration_is_independent_of_coordinate_origin() {
    // G0/G3: a coupled mass matrix under a constant load has exactly constant
    // acceleration. Translating the coordinates cannot change velocity or
    // acceleration, even when q+h*v rounds at a very different scale.
    let mass = [2.0, 0.25, 0.25, 3.0];
    let h = 1e-7;
    let acceleration = [2.0, -0.5];
    let force = [3.875, -1.0];
    for rho in [0.0, 0.3, 0.5, 0.8, 1.0] {
        let method = GeneralizedAlpha::new(&mass, &[0.0; 4], &[0.0; 4], 2, h, rho);
        let mut reference: Option<([f64; 2], [f64; 2])> = None;
        for origin in [0.0, 1.0, 1e6] {
            let (mut q, mut v, mut a) = ([origin; 2], [0.25, -0.125], acceleration);
            for step in 1..=16 {
                galpha_step(&method, &mut q, &mut v, &mut a, &force);
                let t = step as f64 * h;
                for axis in 0..2 {
                    let initial_velocity = [0.25, -0.125][axis];
                    let exact_q = origin + initial_velocity * t + 0.5 * acceleration[axis] * t * t;
                    let exact_v = initial_velocity + acceleration[axis] * t;
                    assert!(
                        (q[axis] - exact_q).abs() <= 32.0 * f64::EPSILON * exact_q.abs().max(1.0)
                    );
                    assert!((v[axis] - exact_v).abs() <= 16.0 * f64::EPSILON);
                    assert!((a[axis] - acceleration[axis]).abs() <= 8.0 * f64::EPSILON);
                }
            }
            if let Some((expected_v, expected_a)) = reference {
                assert_eq!(v.map(f64::to_bits), expected_v.map(f64::to_bits));
                assert_eq!(a.map(f64::to_bits), expected_a.map(f64::to_bits));
            } else {
                reference = Some((v, a));
            }
        }
    }
}

#[test]
fn dense_galpha_rejects_malformed_systems_before_assembly() {
    let rejects = [
        std::panic::catch_unwind(|| GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 0, 0.1, 0.8))
            .is_err(),
        std::panic::catch_unwind(|| {
            GeneralizedAlpha::new(&[1.0], &[0.0; 4], &[1.0; 4], 2, 0.1, 0.8)
        })
        .is_err(),
        std::panic::catch_unwind(|| {
            GeneralizedAlpha::new(&[1.0; 4], &[0.0], &[1.0; 4], 2, 0.1, 0.8)
        })
        .is_err(),
        std::panic::catch_unwind(|| {
            GeneralizedAlpha::new(&[1.0; 4], &[0.0; 4], &[1.0], 2, 0.1, 0.8)
        })
        .is_err(),
        std::panic::catch_unwind(|| GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, 0.0, 0.8))
            .is_err(),
        std::panic::catch_unwind(|| {
            GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, f64::INFINITY, 0.8)
        })
        .is_err(),
    ];
    assert!(rejects.into_iter().all(|rejects| rejects));
}

#[test]
#[should_panic(expected = "generalized-alpha mass shape")]
fn dense_galpha_rejects_oversized_mass_matrix() {
    let _ = GeneralizedAlpha::new(&[1.0; 5], &[0.0; 4], &[1.0; 4], 2, 0.1, 0.8);
}

#[test]
#[should_panic(expected = "generalized-alpha damping shape")]
fn dense_galpha_rejects_oversized_damping_matrix() {
    let _ = GeneralizedAlpha::new(&[1.0; 4], &[0.0; 5], &[1.0; 4], 2, 0.1, 0.8);
}

#[test]
#[should_panic(expected = "generalized-alpha stiffness shape")]
fn dense_galpha_rejects_oversized_stiffness_matrix() {
    let _ = GeneralizedAlpha::new(&[1.0; 4], &[0.0; 4], &[1.0; 5], 2, 0.1, 0.8);
}

#[test]
fn galpha_step_rejects_oversized_state_without_mutation() {
    let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, 0.1, 0.8);
    let mut q = vec![1.0, 99.0];
    let mut v = vec![0.0];
    let mut a = vec![-1.0];
    let f = vec![0.0];
    let before = (q.clone(), v.clone(), a.clone());

    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
        }))
        .is_err()
    );
    assert_eq!(q, before.0);
    assert_eq!(v, before.1);
    assert_eq!(a, before.2);
}

#[test]
#[should_panic(expected = "generalized-alpha matrix dimension overflows usize")]
fn dense_galpha_rejects_overflowing_dimension_before_allocation() {
    let _ = GeneralizedAlpha::new(&[], &[], &[], usize::MAX, 0.1, 0.8);
}

#[test]
fn dense_galpha_rejects_nan_timestep() {
    assert!(
        std::panic::catch_unwind(|| {
            GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, f64::NAN, 0.8)
        })
        .is_err()
    );
}
