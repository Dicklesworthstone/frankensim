//! fs-dcontact conformance battery: polyline geometry, the bouncing
//! analytic fixture (restitution + closed-form max penetration),
//! string-fret rattle energy exactness with an explicit-integrator
//! mutation, one-sidedness mutation, iteration budgets across a
//! velocity sweep, collocation refinement convergence, the jawari
//! spectral-enrichment casebook, determinism, and typed refusals.

use fs_dcontact::{ContactStorage, DContactError, Obstacle, polyline_heights, string_collocation};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage, step};

/// G0: SI line measures, not probability weights, determine total reaction.
#[test]
fn g0_weights_are_physical_measure_refinement_invariant() {
    let (length, density_k, penetration) = (0.73, 2e8, 1e-4_f64);
    let expected = length * density_k * penetration * penetration;
    for cells in [8, 16, 32] {
        let obstacle = Obstacle::new(
            vec![1.0; cells],
            cells,
            1,
            vec![0.0; cells],
            vec![length / cells as f64; cells],
            density_k,
            2.0,
            "synthetic uniform line; K per metre".into(),
        )
        .unwrap();
        let storage =
            ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![obstacle]).unwrap();
        let mut gradient = [0.0; 2];
        storage.gradient(&[penetration, 0.0], &mut gradient);
        assert!((gradient[0] - expected).abs() <= 1e-12 * expected);
        let energy = length * density_k * penetration.powi(3) / 3.0;
        assert!((storage.hamiltonian(&[penetration, 0.0]) - energy).abs() <= 1e-12 * energy);
    }
}

#[test]
fn g0_dissipative_forces_refuse_mismatched_velocity_length() {
    for chi in [0.0, 0.5] {
        let law = Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.0],
            vec![1.0],
            100.0,
            1.5,
            "synthetic shape refusal".into(),
        )
        .unwrap()
        .with_internal_loss(chi)
        .unwrap();
        for wrong in [vec![], vec![0.0, 0.0]] {
            assert!(matches!(
                law.dissipative_modal_forces(1, &[0.01, 0.0], &wrong),
                Err(DContactError::Shape { .. })
            ));
        }
        assert!(
            law.dissipative_modal_forces(1, &[0.01, 0.0], &[0.0])
                .is_ok()
        );
    }
}

// The production potential and production implicit loss port are integrated by
// the existing Gonzalez owner. No explicit force lag or restitution fit enters.
fn isolated_impact(mass: f64, k: f64, chi: f64, incoming: f64, dt: f64) -> (f64, f64) {
    let floor = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        k,
        1.5,
        "synthetic Hertz/HC impact benchmark".into(),
    )
    .unwrap()
    .with_internal_loss(chi)
    .unwrap();
    let loss = ContactStorage::new(Box::new(FreeMass { m: mass }), 1, vec![floor.clone()]).unwrap();
    let storage = ContactStorage::new(Box::new(FreeMass { m: mass }), 1, vec![floor]).unwrap();
    let sys = PortHamiltonian::new(
        2,
        1,
        symplectic_j(1),
        vec![0.0; 4],
        vec![0.0, 1.0],
        Box::new(storage),
    )
    .unwrap();
    let mut workspace = fs_phs::StepWorkspace::new(&sys).unwrap();
    let mut x = [0.0, mass * incoming];
    let mut time = 0.0;
    let port = |state: &[f64], effort: &[f64], out: &mut [f64]| {
        loss.dissipative_flow_into(state, effort, out)
    };
    for _ in 0..10_000 {
        let mut next = [0.0; 2];
        let mut output = [0.0];
        workspace
            .step_into_dissipative_controlled(
                &sys,
                &x,
                &[0.0],
                dt,
                &mut next,
                &mut output,
                &port,
                None,
                None,
                || false,
            )
            .unwrap();
        if x[0] > 0.0 && next[0] <= 0.0 && next[1] < 0.0 {
            // Contact starts exactly at t=0, delta=0, v=incoming. Locate its
            // END by solving the original implicit evolution from this same
            // pre-event state, with the trial duration as the root variable.
            // Interpolating endpoint q/p does not solve that evolution and can
            // pollute the contact-time order and restitution near separation.
            let (mut before, mut after) = (0.0, dt);
            let mut event = next;
            let mut event_dt = dt;
            for _ in 0..56 {
                let trial_dt = f64::midpoint(before, after);
                if trial_dt == before || trial_dt == after {
                    break;
                }
                let mut trial = [0.0; 2];
                workspace
                    .step_into_dissipative_controlled(
                        &sys,
                        &x,
                        &[0.0],
                        trial_dt,
                        &mut trial,
                        &mut output,
                        &port,
                        None,
                        None,
                        || false,
                    )
                    .unwrap();
                if trial[0].abs() < event[0].abs() {
                    event = trial;
                    event_dt = trial_dt;
                }
                if trial[0] > 0.0 {
                    before = trial_dt;
                } else {
                    after = trial_dt;
                }
            }
            // One ppm of a step's travel is negligible against the 0.5%/1e-3
            // physical gates and remains above the Newton roundoff floor.
            assert!(event[0].abs() <= 1e-6 * incoming * dt);
            return (time + event_dt, -event[1] / (mass * incoming));
        }
        x = next;
        time += dt;
    }
    panic!("impact never separated within the declared 10,000-step budget");
}

/// G1: Johnson 1985 section 11.4; the time formula's 2.943 is rounded.
/// Order is measured from rung differences, so that rounding cannot fake a floor.
#[test]
fn g1_hertz_sphere_impact_contact_time_matches_johnson() {
    let (mass, modulus, radius, incoming) = (0.05, 2e9, 0.02_f64, 0.4);
    let k = (4.0 / 3.0) * modulus * radius.sqrt();
    let maximum = fs_math::det::pow(
        15.0 * mass * incoming * incoming / (16.0 * modulus * radius.sqrt()),
        0.4,
    );
    let reference = 2.943 * maximum / incoming;
    // Retain the original three-rung gate while exposing finer rungs to
    // distinguish an event error from a pre-asymptotic discretization band.
    // Native results, not a presumed method order, must adjudicate that claim.
    let mut times = [0.0; 7];
    let rungs = [32, 64, 128, 256, 512, 1024, 2048];
    for (i, cells) in rungs.into_iter().enumerate() {
        let dt = reference / f64::from(cells);
        let (time, restitution) = isolated_impact(mass, k, 0.0, incoming, dt);
        times[i] = time;
        assert!((restitution - 1.0).abs() < 1e-3);
        println!(
            "{{\"suite\":\"normal-contact\",\"case\":\"hertz-impact\",\"dt_s\":{dt},\"contact_time_s\":{time},\"reference_s\":{reference},\"relative_error\":{}}}",
            (time / reference - 1.0).abs()
        );
    }
    for (i, window) in times.windows(3).enumerate() {
        let observed = ((window[0] - window[1]).abs() / (window[1] - window[2]).abs()).log2();
        println!(
            "{{\"suite\":\"normal-contact\",\"case\":\"hertz-impact-order\",\"coarse_cells\":{},\"middle_cells\":{},\"fine_cells\":{},\"observed_order\":{observed},\"order_gate\":1.8}}",
            rungs[i],
            rungs[i + 1],
            rungs[i + 2]
        );
    }
    let order = ((times[0] - times[1]).abs() / (times[1] - times[2]).abs()).log2();
    assert!(
        order.is_finite() && order >= 1.8,
        "contact time order {order}: {times:?}"
    );
    assert!((times[2] / reference - 1.0).abs() <= 0.005);
    assert!((times[6] / reference - 1.0).abs() <= 0.005);
}

#[test]
fn g1_hunt_crossley_restitution_matches_exact_relation() {
    let (mass, modulus, radius, incoming) = (0.05, 2e9, 0.02_f64, 0.4);
    let k = (4.0 / 3.0) * modulus * radius.sqrt();
    let maximum = fs_math::det::pow(
        15.0 * mass * incoming * incoming / (16.0 * modulus * radius.sqrt()),
        0.4,
    );
    let time_scale = 2.943 * maximum / incoming;
    for chi_v in [0.01, 0.1, 0.5] {
        let left = chi_v - fs_math::det::ln(1.0 + chi_v);
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..80 {
            let e = f64::midpoint(lo, hi);
            let right = -chi_v * e - fs_math::det::ln(1.0 - chi_v * e);
            if right < left {
                lo = e;
            } else {
                hi = e;
            }
        }
        let reference = f64::midpoint(lo, hi);
        let mut errors = [0.0; 3];
        for (i, cells) in [64, 128, 256].into_iter().enumerate() {
            let dt = time_scale / f64::from(cells);
            let (_, actual) = isolated_impact(mass, k, chi_v / incoming, incoming, dt);
            errors[i] = (actual - reference).abs();
            println!(
                "{{\"suite\":\"normal-contact\",\"case\":\"hunt-crossley\",\"chi_v\":{chi_v},\"dt_s\":{dt},\"restitution\":{actual},\"reference\":{reference},\"absolute_error\":{}}}",
                errors[i]
            );
        }
        assert!(errors[2] <= 1e-3, "chi*v={chi_v}: {errors:?}");
        assert!(
            errors[2] < errors[0],
            "the HC timestep axis must be live: {errors:?}"
        );
    }
}

/// Free 1-DOF mass (ball): H = p^2/(2 m) — no spring.
struct FreeMass {
    m: f64,
}

impl Storage for FreeMass {
    // The free mass's mathematical midpoint remainder is exactly zero. Supply
    // it explicitly; no nonlinear contact is classified as quadratic and no
    // Newton rounding allowance or physical parameter changes.
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        (a.len() == 2 && b.len() == 2 && mid.len() == 2).then_some(0.0)
    }

    fn hamiltonian(&self, x: &[f64]) -> f64 {
        x[1] * x[1] / (2.0 * self.m)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 0.0;
        out[1] = x[1] / self.m;
    }
}

/// G0: the exact central HC refusal's pre-step state. An opaque wrapper keeps
/// the old subtractive residual as a live mutation; the physical law is shared.
#[test]
fn g0_hunt_crossley_cancellation_converges_under_the_unchanged_request() {
    struct OpaqueMass(FreeMass);
    impl Storage for OpaqueMass {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            self.0.hamiltonian(x)
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            self.0.gradient(x, out);
        }
    }
    let (mass, modulus, radius, incoming) = (0.05, 2e9, 0.02_f64, 0.4);
    let k = (4.0 / 3.0) * modulus * radius.sqrt();
    let maximum = fs_math::det::pow(
        15.0 * mass * incoming * incoming / (16.0 * modulus * radius.sqrt()),
        0.4,
    );
    let dt = 2.943 * maximum / incoming / 256.0;
    let obstacle = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        k,
        1.5,
        "synthetic HC central-refusal regression; SI free mass".into(),
    )
    .unwrap()
    .with_internal_loss(0.5 / incoming)
    .unwrap();
    let storage = |declared| {
        let inner: Box<dyn Storage> = if declared {
            Box::new(FreeMass { m: mass })
        } else {
            Box::new(OpaqueMass(FreeMass { m: mass }))
        };
        ContactStorage::new(inner, 1, vec![obstacle.clone()]).unwrap()
    };
    let loss = storage(true);
    let port = |state: &[f64], effort: &[f64], out: &mut [f64]| {
        loss.dissipative_flow_into(state, effort, out)
    };
    // Recorded at the failing 256-cell trajectory's step120, not an oracle.
    let initial = [5.258_730_595_626_957e-5, -5.472_353_310_761_396e-5];
    let budget = fs_phs::SolveBudget {
        relative: 1e-13,
        absolute: 0.0,
        max_iterations: 50,
    };
    let scale = initial.iter().copied().map(f64::abs).fold(0.0, f64::max);
    let request = budget.threshold(scale);
    for declared in [false, true] {
        let contact = storage(declared);
        assert!(contact.exact_quadratic_matrix().is_none());
        assert!(contact.quadratic_part().is_none_or(|part| !part.exact));
        let sys = PortHamiltonian::new(
            2,
            1,
            symplectic_j(1),
            vec![0.0; 4],
            vec![0.0, 1.0],
            Box::new(contact),
        )
        .unwrap();
        let mut work = fs_phs::StepWorkspace::new(&sys).unwrap();
        work.set_solve_budget(budget).unwrap();
        let (mut next, mut output) = ([123.0; 2], [456.0]);
        let result = work.step_into_dissipative_controlled(
            &sys,
            &initial,
            &[0.0],
            dt,
            &mut next,
            &mut output,
            &port,
            None,
            None,
            || false,
        );
        if !declared {
            assert!(matches!(result,Err(fs_phs::PreparedStepError::Solver(
                fs_phs::PhsError::NewtonRefused {certificate})) if certificate.residual>certificate.requested_budget));
            assert_eq!(next, [123.0; 2]);
            assert_eq!(output, [456.0]);
        } else {
            let record = result.unwrap();
            let certificate = record.solver_certificate.unwrap();
            assert_eq!(certificate.requested_budget.to_bits(), request.to_bits());
            assert_eq!(certificate.effective_budget.to_bits(), request.to_bits());
            assert!(!certificate.floor_raised);
            assert!(record.solver_residual <= request);
            assert!(record.newton_iters < budget.max_iterations);
            let energy_scale = sys.hamiltonian(&initial).abs()
                + sys.hamiltonian(&next).abs()
                + record.dissipated.abs();
            assert!(record.balance_residual().abs() <= 128.0 * f64::EPSILON * energy_scale);
            let hessian = |x: &[f64], d: &[f64], out: &mut [f64]| {
                loss.hessian_vector_with(x, d, out, |_, direction, value| {
                    value[0] = 0.0;
                    value[1] = direction[1] / mass;
                    true
                })
            };
            let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
                loss.dissipative_flow_tangent_into(x, e, dx, de, out)
            };
            let mut analytic = fs_phs::StepWorkspace::new(&sys).unwrap();
            analytic.set_solve_budget(budget).unwrap();
            let (mut other, mut other_output) = ([0.0; 2], [0.0]);
            let other_record = analytic
                .step_into_dissipative_controlled(
                    &sys,
                    &initial,
                    &[0.0],
                    dt,
                    &mut other,
                    &mut other_output,
                    &port,
                    Some(&hessian),
                    Some(&tangent),
                    || false,
                )
                .unwrap();
            assert!(other_record.solver_residual <= request);
            assert_eq!(
                other_record.solver_certificate.unwrap().effective_budget,
                request
            );
            for (a, b) in next.into_iter().zip(other) {
                assert!((a - b).abs() <= 32.0 * request);
            }
            println!(
                "{{\"suite\":\"normal-contact\",\"case\":\"hc-cancellation-source-fix\",\"dt_s\":{dt},\"residual\":{},\"requested_budget\":{request},\"effective_budget\":{},\"floor_raised\":false,\"iterations\":{},\"balance_residual_j\":{},\"verdict\":\"pass\"}}",
                record.solver_residual,
                certificate.effective_budget,
                record.newton_iters,
                record.balance_residual()
            );
        }
    }
}

fn symplectic_j(n_modes: usize) -> Vec<f64> {
    let dim = 2 * n_modes;
    let mut j = vec![0.0; dim * dim];
    for k in 0..n_modes {
        j[(2 * k) * dim + 2 * k + 1] = 1.0;
        j[(2 * k + 1) * dim + 2 * k] = -1.0;
    }
    j
}

/// Modal string as a QuadraticStorage in the interleaved layout
/// (mass-normalized: H = 1/2 sum p^2 + w^2 q^2).
fn string_storage(
    length: f64,
    tension: f64,
    mu: f64,
    n_modes: usize,
) -> (QuadraticStorage, Vec<f64>) {
    let c = fs_math::det::sqrt(tension / mu);
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (1..=n_modes).map(|k| k as f64 * pi * c / length).collect();
    let dim = 2 * n_modes;
    let mut q = vec![0.0; dim * dim];
    for k in 0..n_modes {
        q[(2 * k) * dim + 2 * k] = omegas[k] * omegas[k];
        q[(2 * k + 1) * dim + 2 * k + 1] = 1.0;
    }
    (QuadraticStorage::new(q, dim).expect("q"), omegas)
}

#[test]
fn polyline_heights_match_analytic_geometry() {
    // A tent profile: heights interpolate linearly on each segment.
    let verts = [(0.0, 0.0), (1.0, 2.0), (3.0, 0.0)];
    let h = polyline_heights(&verts, &[0.0, 0.5, 1.0, 2.0, 3.0]).expect("heights");
    let expect = [0.0, 1.0, 2.0, 1.0, 0.0];
    for (a, b) in h.iter().zip(&expect) {
        assert!((a - b).abs() < 1.0e-14);
    }
    // Typed refusals: unsorted, short, out of span.
    assert!(matches!(
        polyline_heights(&[(1.0, 0.0), (0.5, 0.0)], &[0.7]),
        Err(DContactError::Parameter { .. })
    ));
    assert!(polyline_heights(&[(0.0, 0.0)], &[0.0]).is_err());
    assert!(polyline_heights(&verts, &[3.5]).is_err());
}

#[test]
fn bouncing_mass_restitution_and_max_penetration() {
    // Ball (1 kg) falling onto a power-law floor: pure potential
    // contact means restitution EXACTLY 1 (to solver tolerance), and
    // the max penetration follows the closed form
    //   E_impact = w K/(alpha+1) p_max^(alpha+1)
    // up to the gravity work during penetration (bounded in-test).
    let m = 1.0;
    let g = 9.81;
    let (k_c, alpha) = (1.0e6, 1.5);
    let floor = Obstacle::new(
        vec![1.0], // collocation: q IS the displacement toward the floor
        1,
        1,
        vec![0.5], // gap: floor 0.5 m below the start
        vec![1.0],
        k_c,
        alpha,
        "test-fixture: authored (K, alpha), no material claim".to_string(),
    )
    .expect("floor");
    let storage = ContactStorage::new(Box::new(FreeMass { m }), 1, vec![floor]).expect("storage");
    let sys = PortHamiltonian::new(
        2,
        1,
        symplectic_j(1),
        vec![0.0; 4],
        vec![0.0, 1.0],
        Box::new(storage),
    )
    .expect("sys");
    let dt = 1.0e-5;
    let mut x = vec![0.0, 0.0];
    let mut v_impact = 0.0f64;
    let mut v_out = 0.0f64;
    let mut max_pen = 0.0f64;
    let mut in_contact_prev = false;
    for _ in 0..40_000 {
        let rec = step(&sys, &x, &[m * g], dt).expect("step");
        x = rec.x;
        let pen = (x[0] - 0.5).max(0.0);
        max_pen = max_pen.max(pen);
        let in_contact = pen > 0.0;
        if in_contact && !in_contact_prev {
            v_impact = x[1] / m; // velocity entering contact
        }
        if !in_contact && in_contact_prev {
            v_out = x[1] / m; // velocity leaving contact
            break;
        }
        in_contact_prev = in_contact;
    }
    assert!(v_impact > 0.0, "ball never hit the floor");
    let restitution = (v_out / v_impact).abs(); // v_out is negative (rebound)
    assert!(
        (restitution - 1.0).abs() < 1.0e-3,
        "elastic potential restitution {restitution:.6} != 1"
    );
    // Closed-form max penetration from the impact kinetic energy.
    let e_impact = 0.5 * m * v_impact * v_impact;
    let p_closed = fs_math::det::pow(e_impact * (alpha + 1.0) / k_c, 1.0 / (alpha + 1.0));
    // Gravity adds m g p of work during penetration; bound the
    // correction and require agreement within it plus 2%.
    let gravity_correction = m * g * p_closed / e_impact;
    assert!(
        (max_pen - p_closed).abs() <= (0.02 + gravity_correction) * p_closed,
        "max penetration {max_pen:.5e} vs closed form {p_closed:.5e}"
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"bouncing\",\"restitution\":{restitution:.6},\"p_max\":{max_pen:.4e},\"p_closed\":{p_closed:.4e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn hunt_crossley_lowers_restitution() {
    let m = 1.0;
    let g = 9.81;
    let (k_c, alpha) = (1.0e6, 2.0);
    let floor = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.5],
        vec![1.0],
        k_c,
        alpha,
        "test-fixture: authored Hunt-Crossley".to_string(),
    )
    .expect("floor")
    .with_internal_loss(0.8)
    .expect("chi");
    let storage =
        ContactStorage::new(Box::new(FreeMass { m }), 1, vec![floor.clone()]).expect("st");
    let sys = PortHamiltonian::new(
        2,
        1,
        symplectic_j(1),
        vec![0.0; 4],
        vec![0.0, 1.0],
        Box::new(storage),
    )
    .expect("sys");
    let dt = 1.0e-5;
    let mut x = vec![0.0, 0.0];
    let mut v_impact = 0.0f64;
    let mut v_out = 0.0f64;
    let mut in_contact_prev = false;
    for _ in 0..40_000 {
        let v = x[1] / m;
        let extra = floor.dissipative_modal_forces(1, &x, &[v]).unwrap();
        let rec = step(&sys, &x, &[m * g + extra[0]], dt).expect("step");
        x = rec.x;
        let pen = (x[0] - 0.5).max(0.0);
        let in_contact = pen > 0.0;
        if in_contact && !in_contact_prev {
            v_impact = x[1] / m;
        }
        if !in_contact && in_contact_prev {
            v_out = x[1] / m;
            break;
        }
        in_contact_prev = in_contact;
    }
    assert!(v_impact > 0.0, "ball never hit the floor");
    let restitution = (v_out / v_impact).abs();
    assert!(
        restitution < 0.97,
        "Hunt-Crossley must dissipate, restitution {restitution:.4}"
    );
    assert!(
        restitution > 0.2,
        "Hunt-Crossley must still rebound, restitution {restitution:.4}"
    );
    assert!(
        Obstacle::new(vec![1.0], 1, 1, vec![0.0], vec![1.0], 1.0, 2.0, "x".into(),)
            .unwrap()
            .with_internal_loss(-1.0)
            .is_err()
    );
}

#[test]
fn hunt_crossley_unloading_saturates_at_nonattractive_reaction() {
    let obstacle = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.5],
        vec![1.0],
        1.0e6,
        2.0,
        "test-fixture: unilateral Hunt-Crossley saturation".to_string(),
    )
    .expect("obstacle")
    .with_internal_loss(0.8)
    .expect("chi");
    let storage = ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![obstacle.clone()])
        .expect("storage");
    // p = 0.5 m and chi * p_dot = -1.6: the raw 1 + chi * p_dot
    // multiplier would reverse the stored unilateral reaction.
    let x = [1.0, -2.0];
    let mut gradient = [0.0; 2];
    storage.gradient(&x, &mut gradient);
    let port = obstacle.dissipative_modal_forces(1, &x, &[-2.0]).unwrap();
    let net_local_reaction = -gradient[0] + port[0];
    assert_eq!(
        net_local_reaction, 0.0,
        "unloading contact must not attract"
    );

    let extreme = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.5],
        vec![1.0],
        1.0e200,
        2.0,
        "test-fixture: finite extreme-scale Hunt-Crossley".to_string(),
    )
    .expect("extreme obstacle")
    .with_internal_loss(1.0e200)
    .expect("extreme chi");
    let extreme_storage =
        ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![extreme.clone()])
            .expect("extreme storage");
    let extreme_x = [1.0, -5.0e-201];
    let mut extreme_gradient = [0.0; 2];
    extreme_storage.gradient(&extreme_x, &mut extreme_gradient);
    let extreme_port = extreme
        .dissipative_modal_forces(1, &extreme_x, &[-5.0e-201])
        .unwrap();
    let extreme_net = -extreme_gradient[0] + extreme_port[0];
    assert!(extreme_net.is_finite());
    assert_eq!(extreme_net, -0.5 * extreme_gradient[0]);

    // Two independently active collocation points exercise the ordering that
    // a single-point test cannot: clamp each local loss factor before modal
    // projection. The sign-varying dense projection makes the two local
    // factors observable: point 0 fully unloads while point 1 retains half its
    // reaction.
    let multimode = Obstacle::new(
        vec![1.0, 1.0, 1.0, -1.0],
        2,
        2,
        vec![0.5, 0.5],
        vec![1.0, 2.0],
        64.0,
        2.0,
        "test-fixture: independent local Hunt-Crossley saturation".to_string(),
    )
    .expect("multimode obstacle")
    .with_internal_loss(1.0)
    .expect("multimode chi");
    let multimode_port = multimode
        .dissipative_modal_forces(2, &[1.0, 0.0, 0.0, 0.0], &[-1.25, -0.75])
        .unwrap();
    assert_eq!(
        multimode_port,
        [32.0, 0.0],
        "local unloading factors must be clamped before modal projection"
    );
    let nonfinite_port = obstacle.dissipative_modal_forces(1, &[f64::NAN, 0.0], &[0.0]);
    assert!(
        matches!(nonfinite_port, Err(DContactError::Parameter { .. })),
        "DEF-07: a fallible port must refuse invalid state instead of erasing contact"
    );
}

/// The string + flat-fret fixture used by several tests.
fn fret_system(n_modes: usize, gap: f64, k_c: f64) -> (PortHamiltonian, ContactStorage, Vec<f64>) {
    let (length, tension, mu) = (0.65, 70.0, 5.0e-3);
    let (storage, omegas) = string_storage(length, tension, mu, n_modes);
    // 24 fret-line collocation points along the first two thirds.
    let points: Vec<f64> = (1i32..=24).map(|i| length * f64::from(i) / 36.0).collect();
    let coll = string_collocation(length, mu, &points, n_modes).expect("collocation");
    let seg = length / 36.0;
    let fret = Obstacle::new(
        coll.clone(),
        points.len(),
        n_modes,
        vec![gap; points.len()],
        vec![seg; points.len()],
        k_c,
        1.5,
        "test-fixture: authored fret line".to_string(),
    )
    .expect("fret");
    let contact = ContactStorage::new(
        Box::new(string_storage(0.65, 70.0, 5.0e-3, n_modes).0),
        n_modes,
        vec![fret.clone()],
    )
    .expect("contact");
    let sys = PortHamiltonian::new(
        2 * n_modes,
        0,
        symplectic_j(n_modes),
        vec![0.0; 4 * n_modes * n_modes],
        vec![],
        Box::new(ContactStorage::new(Box::new(storage), n_modes, vec![fret]).expect("cs")),
    )
    .expect("sys");
    (sys, contact, omegas)
}

/// Pluck-like initial state: displacement toward the fret in the low
/// modes.
fn pluck(n_modes: usize, amp: f64) -> Vec<f64> {
    let mut x = vec![0.0; 2 * n_modes];
    for k in 0..4.min(n_modes) {
        x[2 * k] = amp / (k + 1) as f64;
    }
    x
}

#[test]
fn fret_rattle_conserves_energy_and_explicit_mutation_grows_it() {
    let n_modes = 8;
    let (sys, contact, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let x0 = pluck(n_modes, 6.0e-3);
    let h0 = sys.hamiltonian(&x0);
    let dt = 2.0e-6;
    let steps = 20_000;
    let mut x = x0.clone();
    let mut worst = 0.0f64;
    let mut contact_events = 0usize;
    let mut was_active = false;
    let mut max_pen = 0.0f64;
    for _ in 0..steps {
        let rec = step(&sys, &x, &[], dt).expect("step");
        x = rec.x;
        let probe = contact.probe(&x);
        max_pen = max_pen.max(probe.max_penetration);
        let active = probe.active_points > 0;
        if active && !was_active {
            contact_events += 1;
        }
        was_active = active;
        worst = worst.max((sys.hamiltonian(&x) - h0).abs());
    }
    assert!(
        contact_events >= 3,
        "rattle must actually rattle: {contact_events}"
    );
    assert!(
        worst <= 1.0e-8 * h0,
        "discrete-gradient contact drifted {worst:.3e} of {h0:.3e}"
    );
    // MUTATION: explicit symplectic-Euler on the same vector field at
    // the same dt visibly pumps energy through the stiff contact.
    let (_, contact2, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let mut xe = x0.clone();
    let mut grad = vec![0.0; 2 * n_modes];
    let mut h_explicit_max = h0;
    for _ in 0..steps {
        contact2.gradient(&xe, &mut grad);
        for k in 0..n_modes {
            xe[2 * k + 1] -= dt * grad[2 * k];
        }
        contact2.gradient(&xe, &mut grad);
        for k in 0..n_modes {
            xe[2 * k] += dt * grad[2 * k + 1];
        }
        h_explicit_max = h_explicit_max.max(contact2.hamiltonian(&xe));
        if !h_explicit_max.is_finite() {
            break;
        }
    }
    let explicit_growth = (h_explicit_max - h0) / h0;
    assert!(
        !h_explicit_max.is_finite() || explicit_growth > 1.0e3 * (worst / h0).max(1.0e-12),
        "explicit integration must visibly grow energy: growth {explicit_growth:.3e} vs DG drift {:.3e}",
        worst / h0
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"fret-rattle\",\"events\":{contact_events},\"max_penetration\":{max_pen:.3e},\"dg_drift\":{:.3e},\"explicit_growth\":{explicit_growth:.3e},\"verdict\":\"pass\"}}",
        worst / h0
    );
}

/// Mutant potential WITHOUT the one-sided clamp: separation becomes
/// attraction — the exact bug class the [.]_+ exists to prevent.
struct TwoSidedMutant {
    inner: ContactStorage,
}

impl Storage for TwoSidedMutant {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        // Same structure energy; contact term evaluated WITHOUT clamp
        // (the |p|-weighted odd form the bug class produces).
        let ob = &self.inner.obstacles()[0];
        let nm = self.inner.n_modes();
        let mut h = self.inner.inner_storage().hamiltonian(x);
        for i in 0..ob.n_points() {
            let mut disp = 0.0;
            for k in 0..nm {
                disp += ob.collocation()[i * nm + k] * x[2 * k];
            }
            let p = disp - ob.gaps()[i]; // NO clamp
            h += ob.weights()[i] * ob.stiffness() / 3.0 * p * p * p.abs();
        }
        h
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.inner.inner_storage().gradient(x, out);
        let ob = &self.inner.obstacles()[0];
        let nm = self.inner.n_modes();
        for i in 0..ob.n_points() {
            let mut disp = 0.0;
            for k in 0..nm {
                disp += ob.collocation()[i * nm + k] * x[2 * k];
            }
            let p = disp - ob.gaps()[i];
            let f = ob.weights()[i] * ob.stiffness() * p * p.abs();
            for k in 0..nm {
                out[2 * k] += f * ob.collocation()[i * nm + k];
            }
        }
    }
}

#[test]
fn dropped_one_sidedness_detected_as_attraction() {
    // Probe: a state SEPARATED from the obstacle must feel zero
    // contact force from the true storage, nonzero (attractive) from
    // the mutant.
    let floor = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.5],
        vec![1.0],
        1.0e6,
        2.0,
        "test".to_string(),
    )
    .expect("floor");
    let real =
        ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![floor.clone()]).expect("real");
    let mutant = TwoSidedMutant {
        inner: ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![floor]).expect("m"),
    };
    // Sweep separated depths: true force identically zero everywhere
    // on the free side, mutant force present and pulling TOWARD the
    // obstacle (positive q-gradient = restoring toward +q contact).
    for sep in [0.1, 0.3, 0.6] {
        let x_separated = vec![0.5 - sep - 0.5, 0.0]; // sep below the gap
        let mut g_real = vec![0.0; 2];
        let mut g_mut = vec![0.0; 2];
        real.gradient(&x_separated, &mut g_real);
        mutant.gradient(&x_separated, &mut g_mut);
        assert!(
            g_real[0].abs() < f64::MIN_POSITIVE,
            "separated true contact force must vanish at sep {sep}"
        );
        assert!(
            g_mut[0].abs() > 1.0e2,
            "mutant must exert spurious force at sep {sep}: {:.3e}",
            g_mut[0]
        );
    }
}

#[test]
fn iteration_budget_held_across_velocity_sweep_and_stall_is_typed() {
    let n_modes = 8;
    let (sys, _, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let dt = 2.0e-6;
    let mut histogram = [0usize; 8]; // buckets of 4 iterations
    let mut max_iters = 0usize;
    for amp in [2.0e-3, 6.0e-3, 1.5e-2] {
        let mut x = pluck(n_modes, amp);
        for _ in 0..4000 {
            let rec = step(&sys, &x, &[], dt).expect("step within budget");
            max_iters = max_iters.max(rec.newton_iters);
            histogram[(rec.newton_iters / 4).min(7)] += 1;
            x = rec.x;
        }
    }
    // Interior bound (review finding: <= 50 is vacuous — a successful
    // step CANNOT report more than the fs-phs cap; a stall is already
    // a panic above). Histogram-measured headroom at authoring: all
    // steps <= 8 iterations.
    assert!(max_iters <= 24, "iteration budget regressed: {max_iters}");
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"iteration-budget\",\"max_iters\":{max_iters},\"histogram\":{histogram:?},\"verdict\":\"pass\"}}"
    );
    // The stall path is TYPED: absurd stiffness with a huge step
    // refuses by name instead of looping or lying.
    let (sys_stiff, _, _) = fret_system(4, 1.0e-6, 1.0e18);
    let x = pluck(4, 5.0e-2);
    let out = step(&sys_stiff, &x, &[], 1.0e-2);
    assert!(
        matches!(out, Err(fs_phs::PhsError::NewtonStalled { .. })),
        "expected NewtonStalled, got {out:?}"
    );
}

#[test]
fn collocation_refinement_converges() {
    // Contact energy of a FIXED penetrating string shape under a flat
    // obstacle: refining the collocation quadrature converges (the
    // 8->32 difference dominates the 32->128 difference).
    let (length, mu) = (0.65, 5.0e-3);
    let n_modes = 4;
    let shape = |x: &mut Vec<f64>| {
        x[0] = 8.0e-3; // fundamental pushed into the obstacle
        x[4] = 2.0e-3;
    };
    let energy_at = |n_pts: usize| -> f64 {
        let pts: Vec<f64> = (1..=n_pts)
            .map(|i| length * i as f64 / (n_pts as f64 + 1.0))
            .collect();
        let coll = string_collocation(length, mu, &pts, n_modes).expect("collocation");
        let seg = length / (n_pts as f64 + 1.0);
        let ob = Obstacle::new(
            coll,
            n_pts,
            n_modes,
            vec![1.0e-4; n_pts],
            vec![seg; n_pts],
            1.0e7,
            1.5,
            "test".to_string(),
        )
        .expect("ob");
        let cs = ContactStorage::new(
            Box::new(string_storage(length, 70.0, mu, n_modes).0),
            n_modes,
            vec![ob],
        )
        .expect("cs");
        let mut x = vec![0.0; 2 * n_modes];
        shape(&mut x);
        cs.probe(&x).contact_energy
    };
    let (e8, e32, e128) = (energy_at(8), energy_at(32), energy_at(128));
    let d_coarse = (e32 - e8).abs();
    let d_fine = (e128 - e32).abs();
    assert!(e128 > 0.0, "fixture must penetrate");
    // Measured ratio 0.0099 at authoring (near O(h^2)); 0.1 keeps 10x
    // headroom while catching an O(h)-degradation (ratio ~0.25).
    assert!(
        d_fine < 0.1 * d_coarse.max(1.0e-30),
        "refinement not converging: {e8:.6e} {e32:.6e} {e128:.6e}"
    );
}

#[test]
fn jawari_bridge_enriches_high_band_vs_clean_termination() {
    // The sitar/tanpura jawari: a GRADED bridge surface near the
    // termination repeatedly grazes the string and pumps energy
    // upward. The CONTROL is the clean knife-edge termination — the
    // plain string with NO collision surface (a first attempt used a
    // tiny-gap rattling point as "hard bridge" and measured MORE high
    // band than the jawari, executed: a rattle point is itself a
    // collision exciter, not a termination). Quantified as the
    // high-band energy fraction after the transient.
    let (length, tension, mu) = (0.65, 70.0, 5.0e-3);
    // n = 12 modes / 16k steps at 2.5 us: ~4 fundamental periods of
    // grazing — enough for the enrichment signature while keeping the
    // debug-build implicit solve inside the suite budget (n = 20 with
    // 60k steps measured ~25 minutes).
    let n_modes = 12;
    let dt = 2.5e-6;
    let steps = 16_000;
    let run = |obstacle: Obstacle| -> Vec<f64> {
        let (storage, omegas) = string_storage(length, tension, mu, n_modes);
        let cs = ContactStorage::new(Box::new(storage), n_modes, vec![obstacle]).expect("cs");
        let sys = PortHamiltonian::new(
            2 * n_modes,
            0,
            symplectic_j(n_modes),
            vec![0.0; 4 * n_modes * n_modes],
            vec![],
            Box::new(cs),
        )
        .expect("sys");
        // Pluck AWAY from the obstacle side: the string swings back
        // through neutral and GRAZES the bridge each period (a
        // positive pluck started 20 gap-depths INSIDE the profile —
        // executed unphysical preload, energies 100x the pluck's).
        let mut x = pluck(n_modes, -4.0e-3);
        for _ in 0..steps {
            x = step(&sys, &x, &[], dt).expect("step").x;
        }
        (0..n_modes)
            .map(|k| {
                f64::midpoint(
                    x[2 * k + 1] * x[2 * k + 1],
                    omegas[k] * omegas[k] * x[2 * k] * x[2 * k],
                )
            })
            .collect()
    };
    // Jawari: parabolic graded profile over x in [0.02, 0.08] m with
    // gaps from 0.15 mm (near the end) to 1.5 mm, via the polyline
    // helper.
    let jaw_pts: Vec<f64> = (0i32..12)
        .map(|i| 0.02 + 0.06 * f64::from(i) / 11.0)
        .collect();
    let profile: Vec<(f64, f64)> = (0i32..12)
        .map(|i| {
            let t = f64::from(i) / 11.0;
            (0.02 + 0.06 * t, 1.5e-4 + 1.35e-3 * t * t)
        })
        .collect();
    let gaps = polyline_heights(&profile, &jaw_pts).expect("gaps");
    let seg = 0.06 / 12.0;
    let jawari = Obstacle::new(
        string_collocation(length, mu, &jaw_pts, n_modes).expect("collocation"),
        jaw_pts.len(),
        n_modes,
        gaps,
        vec![seg; jaw_pts.len()],
        1.0e8,
        // alpha = 2 (C^2 force): alpha = 1.5 has an UNBOUNDED contact
        // Hessian at the boundary (d2f ~ p^-0.5) and the tight-graze
        // jawari regime stalled the FD-Jacobian Newton on it
        // (executed); exponent 2 is inside the published musical-
        // collision range and removes the singularity.
        2.0,
        "test-fixture: jawari-class graded profile (authored)".to_string(),
    )
    .expect("jawari");
    // Clean termination: zero contact stiffness = the plain linear
    // string (the model's fixed ends ARE the knife-edge bridge).
    let clean = Obstacle::new(
        string_collocation(length, mu, &[0.02], n_modes).expect("collocation"),
        1,
        n_modes,
        vec![1.5e-4],
        vec![1.0],
        0.0,
        2.0,
        "test-fixture: zero-stiffness control (clean termination)".to_string(),
    )
    .expect("clean");
    let e_jaw = run(jawari);
    let e_hard = run(clean);
    let frac = |e: &[f64]| -> f64 {
        let hi: f64 = e[6..].iter().sum();
        let total: f64 = e.iter().sum();
        hi / total.max(f64::MIN_POSITIVE)
    };
    let (f_jaw, f_clean) = (frac(&e_jaw), frac(&e_hard));
    // The linear control PRESERVES the pluck's band split exactly
    // (nothing couples modes), so its high fraction is the initial
    // one; the jawari must sit far above it AND carry an absolute
    // high-band share.
    assert!(
        f_jaw > 10.0 * f_clean.max(1.0e-9) && f_jaw > 0.05,
        "jawari must enrich the high band: {f_jaw:.4} vs clean {f_clean:.4}"
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"jawari\",\"high_band_fraction_jawari\":{f_jaw:.4},\"high_band_fraction_clean\":{f_clean:.4},\"band_energies_jawari\":{e_jaw:?},\"band_energies_clean\":{e_hard:?},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn contact_gradient_matches_finite_difference_of_h() {
    // THE oracle the conservation test cannot be (review finding: the
    // Gonzalez correction forces dg.dx = dH for ANY gradient, so
    // energy conservation is structurally blind to contact-gradient
    // errors): central finite differences of the coded H must match
    // gradient() on a MULTI-POINT, NON-UNIFORM-WEIGHT obstacle at a
    // state with a mixed active/inactive contact set.
    let (length, mu) = (0.65, 5.0e-3);
    let n_modes = 5;
    let pts: Vec<f64> = (1i32..=9).map(|i| length * f64::from(i) / 10.0).collect();
    let coll = string_collocation(length, mu, &pts, n_modes).expect("collocation");
    let weights: Vec<f64> = (1i32..=9).map(|i| 0.01 * f64::from(i)).collect(); // non-uniform
    let gaps: Vec<f64> = (0..9)
        .map(|i| if i % 2 == 0 { 1.0e-4 } else { 5.0e-2 }) // half active, half far
        .collect();
    let ob = Obstacle::new(
        coll,
        9,
        n_modes,
        gaps,
        weights,
        3.0e6,
        1.7,
        "test".to_string(),
    )
    .expect("ob");
    let cs = ContactStorage::new(
        Box::new(string_storage(length, 70.0, mu, n_modes).0),
        n_modes,
        vec![ob],
    )
    .expect("cs");
    let mut x = vec![0.0; 2 * n_modes];
    for k in 0..n_modes {
        x[2 * k] = 3.0e-3 / (k + 1) as f64;
        x[2 * k + 1] = 0.1 * (k as f64 - 2.0);
    }
    let probe = cs.probe(&x);
    assert!(
        probe.active_points >= 2 && probe.active_points <= 8,
        "fixture needs a MIXED active set, got {}",
        probe.active_points
    );
    let mut g = vec![0.0; 2 * n_modes];
    cs.gradient(&x, &mut g);
    let scale = g.iter().fold(1.0e-30f64, |a, &v| a.max(v.abs()));
    for i in 0..2 * n_modes {
        let h_step = 1.0e-7 * (1.0 + x[i].abs());
        let mut xp = x.clone();
        let mut xm = x.clone();
        xp[i] += h_step;
        xm[i] -= h_step;
        let fd = (cs.hamiltonian(&xp) - cs.hamiltonian(&xm)) / (2.0 * h_step);
        assert!(
            (fd - g[i]).abs() <= 1.0e-5 * scale,
            "gradient[{i}] {:.6e} vs FD {fd:.6e}",
            g[i]
        );
    }
    // Probe-vs-Hamiltonian consistency (review finding: duplicated
    // potential formula): contact energy equals the storage split.
    let inner_h = string_storage(length, 70.0, mu, n_modes).0.hamiltonian(&x);
    assert!(
        (probe.contact_energy - (cs.hamiltonian(&x) - inner_h)).abs()
            <= 1.0e-12 * probe.contact_energy.max(1.0),
        "probe energy diverged from the Hamiltonian split"
    );
}

#[test]
fn determinism_bitwise() {
    let n_modes = 6;
    let (sys, _, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let x0 = pluck(n_modes, 6.0e-3);
    let runs: Vec<Vec<f64>> = (0..2)
        .map(|_| {
            let mut x = x0.clone();
            for _ in 0..2000 {
                x = step(&sys, &x, &[], 2.0e-6).expect("step").x;
            }
            x
        })
        .collect();
    for (a, b) in runs[0].iter().zip(&runs[1]) {
        assert_eq!(a.to_bits(), b.to_bits(), "nondeterministic step");
    }
}

#[test]
fn refusals_are_typed() {
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            -1.0,
            1.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            1.0,
            0.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            f64::NAN,
            1.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(
        Obstacle::new(
            vec![1.0, 2.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            1.0,
            1.5,
            String::new()
        )
        .is_err()
    );
}

#[test]
fn contact_pack_migration_is_provenance_only() {
    // The 3ez8g.13.1 wiring e2e: the jawari fixture, rebuilt from the
    // RECEIPTED contact card (matdb pack contact-jawari-bone-bridge)
    // instead of authored (K, alpha), must produce BITWISE-identical
    // physics — the migration touches provenance and nothing else.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf();
    let card =
        fs_matdb::ContactLawCard::load(&root.join("data/matdb/seed-v1/contact-jawari-bone-bridge"))
            .expect("committed jawari pack");
    // In-validity operating point (the fixture's regime); the
    // out-of-validity refusal falsifier executes in fs-matdb cp-003.
    let receipt = card.lookup(1.0, 1.0).expect("in-validity lookup");
    let maximum = card.force_validity_n.1;
    assert!(card.lookup(maximum, 1.0).is_ok());
    assert!(matches!(
        card.lookup(maximum.next_up(), 1.0),
        Err(fs_matdb::ContactPackError::OutsideValidity { .. })
    ));
    let (length, tension, mu) = (0.65, 70.0, 5.0e-3);
    let n_modes = 12;
    let dt = 2.5e-6;
    let steps = 2_000;
    let jaw_pts: Vec<f64> = (0i32..12)
        .map(|i| 0.02 + 0.06 * f64::from(i) / 11.0)
        .collect();
    let profile: Vec<(f64, f64)> = (0i32..12)
        .map(|i| {
            let t = f64::from(i) / 11.0;
            (0.02 + 0.06 * t, 1.5e-4 + 1.35e-3 * t * t)
        })
        .collect();
    let gaps = polyline_heights(&profile, &jaw_pts).expect("gaps");
    let seg = 0.06 / 12.0;
    let coll = string_collocation(length, mu, &jaw_pts, n_modes).expect("collocation");
    let authored = Obstacle::new(
        coll.clone(),
        jaw_pts.len(),
        n_modes,
        gaps.clone(),
        vec![seg; jaw_pts.len()],
        1.0e8,
        2.0,
        "test-fixture: jawari-class graded profile (authored)".to_string(),
    )
    .expect("authored");
    let receipted = Obstacle::from_receipt(
        coll,
        jaw_pts.len(),
        n_modes,
        gaps,
        vec![seg; jaw_pts.len()],
        &receipt,
    )
    .expect("receipted");
    assert!(
        receipted
            .provenance()
            .starts_with("matdb:contact-jawari-bone-bridge"),
        "typed receipt provenance: {}",
        receipted.provenance()
    );
    let run = |obstacle: Obstacle| -> Vec<f64> {
        let (storage, _) = string_storage(length, tension, mu, n_modes);
        let cs = ContactStorage::new(Box::new(storage), n_modes, vec![obstacle]).expect("cs");
        let sys = PortHamiltonian::new(
            2 * n_modes,
            0,
            symplectic_j(n_modes),
            vec![0.0; 4 * n_modes * n_modes],
            vec![],
            Box::new(cs),
        )
        .expect("sys");
        let mut x = pluck(n_modes, -4.0e-3);
        for _ in 0..steps {
            x = step(&sys, &x, &[], dt).expect("step").x;
        }
        x
    };
    let (xa, xr) = (run(authored), run(receipted));
    for (k, (a, r)) in xa.iter().zip(&xr).enumerate() {
        assert!(
            a.to_bits() == r.to_bits(),
            "state {k} diverged: {a:e} vs {r:e} — the migration is not provenance-only"
        );
    }
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"pack-migration\",\"pack\":\"{}\",\
         \"pair\":\"{}\",\"steps\":{steps},\"verdict\":\"bitwise-identical\"}}",
        receipt.pack_id, receipt.pair_label
    );
}
