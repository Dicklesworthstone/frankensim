//! DEF-05, G0/G4. The old f9c6f1fcb fallbacks admitted residuals below
//! 1e-6*scale even after the declared target had failed. A consistent staircase
//! Hamiltonian manufactures a missing smooth Newton root at a gradient jump.
//! Native comparison with that historical revision remains batch verification;
//! the tests below prove the current refusal and its quantitative certificate.

use fs_phs::{
    NewtonSite, PhsError, PortHamiltonian, QuadraticStorage, SolveBudget, SolveCertificate,
    StepWorkspace, Storage, admit_default_newton_iterate, admit_newton_iterate, discrete_gradient,
    join_port, kirchhoff_parallel_step, mass_spring_damper, step, step_descriptor,
};

const DELTA: f64 = 1e-9;
const Q: f64 = 1e-7;
const INITIAL: [f64; 2] = [0.1, 0.3];
const DT: f64 = 0.01;

struct Staircase;
impl Storage for Staircase {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let k = (x[0] / Q).floor();
        0.5 * (x[0] * x[0] + x[1] * x[1]) + DELTA * (Q * k * (k - 1.0) / 2.0 + k * (x[0] - k * Q))
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + DELTA * (x[0] / Q).floor();
        out[1] = x[1];
    }
}

fn staircase(ports: usize) -> PortHamiltonian {
    let mut g = vec![0.0; 2 * ports];
    g[0] = 1.0;
    g[ports + 1] = 1.0;
    PortHamiltonian::new(
        2,
        ports,
        vec![0.0; 4],
        vec![1.0, 0.0, 0.0, 1.0],
        g,
        Box::new(Staircase),
    )
    .unwrap()
}

fn input() -> [f64; 2] {
    // Manufacture forcing halfway between the two one-sided DG efforts at the
    // same physical target, never from a solver's accepted/error return value.
    let target = [0.2, 0.4];
    let mut left = target;
    let mut right = target;
    left[0] -= Q * 1e-4;
    right[0] += Q * 1e-4;
    let a = discrete_gradient(&Staircase, &INITIAL, &left);
    let b = discrete_gradient(&Staircase, &INITIAL, &right);
    [
        (target[0] - INITIAL[0]) / DT + f64::midpoint(a[0], b[0]),
        (target[1] - INITIAL[1]) / DT + f64::midpoint(a[1], b[1]),
    ]
}

fn refused(error: PhsError, site: NewtonSite) -> SolveCertificate {
    let PhsError::NewtonRefused { certificate } = error else {
        panic!("expected certificate refusal: {error:?}")
    };
    assert_eq!(certificate.site, site);
    assert!(certificate.residual > certificate.effective_budget);
    assert!(
        certificate.residual < 1e-6 * certificate.scale,
        "fixture must remain inside the historical relaxed band: {certificate:?}"
    );
    assert!(certificate.iterations <= 50);
    emit("stalled-above-budget", certificate, "refuse");
    certificate
}

fn emit(case: &str, c: SolveCertificate, verdict: &str) {
    let git = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let git_sha = String::from_utf8(git.stdout).unwrap();
    println!(
        "{{\"suite\":\"newton-refusal\",\"case\":\"{case}\",\"site\":\"{:?}\",\"residual\":{},\"budget\":{},\"requested_budget\":{},\"floor\":{},\"floor_raised\":{},\"iterations\":{},\"verdict\":\"{verdict}\",\"git_sha\":\"{}\"}}",
        c.site,
        c.residual,
        c.effective_budget,
        c.requested_budget,
        c.rounding_floor,
        c.floor_raised,
        c.iterations,
        git_sha.trim()
    );
}

#[test]
fn step_refuses_stalled_iterate_above_budget() {
    refused(
        step(&staircase(2), &INITIAL, &input(), DT).unwrap_err(),
        NewtonSite::Step,
    );
}

#[test]
fn step_descriptor_refuses_stalled_iterate_above_budget() {
    // A zero-coupling extra port attaches an independent constrained scalar.
    // It makes a nonsingular descriptor block without changing the two-state
    // staircase equation or its two external forcing columns.
    let b = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![1.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let sys = join_port(staircase(3), b, 2, 0).unwrap();
    refused(
        step_descriptor(&sys, &[INITIAL[0], INITIAL[1], 0.0, 0.0], &input(), DT).unwrap_err(),
        NewtonSite::Descriptor,
    );
}

#[test]
fn prepared_step_refuses_stalled_iterate_above_budget() {
    for mode in 0..3 {
        let sys = staircase(2);
        let mut work = StepWorkspace::new(&sys).unwrap();
        if mode == 2 {
            work.set_condensed_pairs(&[[0, 1]]).unwrap();
        }
        let (mut next, mut y) = ([123.0, 456.0], [789.0, 1011.0]);
        let error = if mode == 0 {
            work.step_into(&sys, &INITIAL, &input(), DT, &mut next, &mut y)
        } else {
            work.step_into_analytic(
                &sys,
                &INITIAL,
                &input(),
                DT,
                &mut next,
                &mut y,
                &|_, direction, out| {
                    out.copy_from_slice(direction);
                    true
                },
            )
        }
        .unwrap_err();
        refused(error, NewtonSite::Prepared);
        assert_eq!(next, [123.0, 456.0]);
        assert_eq!(y, [789.0, 1011.0]);
    }
}

#[test]
fn kirchhoff_parallel_step_refuses_relaxed_split() {
    let system = || {
        PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![0.0],
            vec![1.0],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap()
    };
    // At dt=0 no flow split can change the two unequal existing efforts. The
    // old final 1e-6 gate accepted their 2e-7 mismatch after a zero derivative.
    let c = refused(
        kirchhoff_parallel_step(&system(), &[2e-7], &system(), &[0.0], 0.0, 0.0).unwrap_err(),
        NewtonSite::KirchhoffParallel,
    );
    assert_eq!(c.residual, 2e-7);
    assert_eq!(c.iterations, 0);
}

fn certificate(residual: f64, budget: f64, floor: f64) -> SolveCertificate {
    SolveBudget {
        relative: 0.0,
        absolute: budget,
        max_iterations: 50,
    }
    .certificate(NewtonSite::Step, residual, 1.0, floor, 1, 0)
}

#[test]
fn budget_boundary_accepts_at_cap_refuses_at_next_up() {
    for (budget, floor) in [(1e-12, 0.0), (1e-12, 1e-12), (2e-12, 1e-12)] {
        let cap = certificate(budget, budget, floor);
        assert_eq!(admit_newton_iterate(cap).unwrap(), cap);
        assert!(matches!(
            admit_newton_iterate(certificate(budget.next_up(), budget, floor)),
            Err(PhsError::NewtonRefused { .. })
        ));
        emit("at-explicit-cap", cap, "accept");
    }
}

#[test]
fn budget_below_rounding_floor_refuses_typed() {
    let c = certificate(0.0, 1e-12, 1e-12_f64.next_up());
    assert_eq!(
        admit_newton_iterate(c),
        Err(PhsError::BudgetBelowRoundingFloor { certificate: c })
    );
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(admit_newton_iterate(certificate(bad, 1e-12, 0.0)).is_err());
    }
    assert!(
        SolveBudget {
            relative: f64::NAN,
            ..SolveBudget::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        SolveBudget {
            max_iterations: 51,
            ..SolveBudget::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        SolveBudget {
            max_iterations: 50,
            ..SolveBudget::default()
        }
        .validate()
        .is_ok()
    );
    emit("explicit-below-floor", c, "refuse");
    let sys = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
    let mut work = StepWorkspace::new(&sys).unwrap();
    work.set_solve_budget(SolveBudget {
        relative: 0.0,
        absolute: 0.0,
        max_iterations: 50,
    })
    .unwrap();
    let (mut next, mut output) = ([123.0, 456.0], [789.0]);
    let error = work
        .step_into(&sys, &[1.0, 0.0], &[0.0], 0.0, &mut next, &mut output)
        .unwrap_err();
    let PhsError::BudgetBelowRoundingFloor { certificate: c } = error else {
        panic!("{error:?}")
    };
    assert_eq!(c.requested_budget, 0.0);
    assert!(c.rounding_floor > 0.0);
    assert_eq!(next, [123.0, 456.0]);
    assert_eq!(output, [789.0]);
    emit("prepared-explicit-below-floor", c, "refuse");
}

#[test]
fn iteration_limit_n_accepts_n_plus_one_refuses() {
    let sys = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([123.0], [456.0]);
    let hessian = |_: &[f64], d: &[f64], out: &mut [f64]| {
        out.copy_from_slice(d);
        true
    };
    work.set_iteration_limit(0).unwrap();
    let error = work
        .step_into_analytic(&sys, &[0.0], &[1.0], 1.0, &mut next, &mut y, &hessian)
        .unwrap_err();
    let PhsError::NewtonRefused { certificate: c } = error else {
        panic!("{error:?}")
    };
    assert_eq!(c.iterations, 0);
    assert_eq!(next, [123.0]);
    assert_eq!(y, [456.0]);
    work.set_iteration_limit(1).unwrap();
    let record = work
        .step_into_analytic(&sys, &[0.0], &[1.0], 1.0, &mut next, &mut y, &hessian)
        .unwrap();
    assert_eq!(record.newton_iters, 1);
    assert_eq!(next, [1.0]);
    assert_eq!(y, [0.5]);
    emit(
        "at-update-cap",
        record.solver_certificate.unwrap(),
        "accept",
    );
}

#[test]
fn noise_floor_converged_iterate_still_accepted() {
    // The reed casebook's actual physical lamella and load; no fabricated
    // residual or hand-typed allowance enters the stepping path.
    let sys = mass_spring_damper(3e-6, 120.0, 0.02).unwrap();
    let mut x = vec![0.0, 0.0];
    let mut raised = None;
    // Run through the stationary tail. Stopping at the first raised certificate
    // hid the later best-at-entry refusal in the coupled reed casebook.
    for _ in 0..60_000 {
        let record = step(&sys, &x, &[0.03], 2e-5).unwrap();
        let c = record.solver_certificate.unwrap();
        assert_eq!(c.requested_budget, 1e-13 * c.scale);
        assert_eq!(c.effective_budget, c.requested_budget.max(c.rounding_floor));
        assert!(c.residual <= c.effective_budget);
        if c.floor_raised
            && c.residual > c.requested_budget
            && (record.x[0] - 2.5e-4).abs() <= 1e-6 * 2.5e-4
        {
            raised = Some(c);
        }
        x = record.x;
    }
    assert!((x[0] - 2.5e-4).abs() <= 1e-6 * 2.5e-4);
    assert!(x[1].abs() <= 1e-8);
    let c = raised.expect("quasi-static reed must expose its computed roundoff allowance");
    assert!(c.rounding_floor > c.requested_budget);
    assert_eq!(admit_default_newton_iterate(c).unwrap(), c);
    let mut explicit = c;
    explicit.floor_raised = false;
    explicit.effective_budget = explicit.requested_budget;
    assert_eq!(
        admit_newton_iterate(explicit),
        Err(PhsError::BudgetBelowRoundingFloor {
            certificate: explicit
        })
    );
    let mut above_cap = c;
    above_cap.residual = above_cap.effective_budget.next_up();
    assert!(matches!(
        admit_default_newton_iterate(above_cap),
        Err(PhsError::NewtonRefused { .. })
    ));
    let mut mutant = c;
    mutant.rounding_floor = 0.0;
    assert!(
        matches!(
            admit_default_newton_iterate(mutant),
            Err(PhsError::NewtonRefused { .. })
        ),
        "dropping the floor must not leave a self-consistent raised certificate"
    );
    emit("quasi-static-default-floor", c, "accept");
}

#[test]
fn refusal_leaves_state_and_outputs_untouched() {
    let sys = staircase(2);
    let initial = INITIAL;
    let u = input();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = (
        [f64::from_bits(0x7ff8_0000_0000_1234), -0.0],
        [789.0, 1011.0],
    );
    let before = (next.map(f64::to_bits), y.map(f64::to_bits));
    let error = work
        .step_into(&sys, &initial, &u, DT, &mut next, &mut y)
        .unwrap_err();
    let c = refused(error, NewtonSite::Prepared);
    assert_eq!(next.map(f64::to_bits), before.0);
    assert_eq!(y.map(f64::to_bits), before.1);
    assert_eq!(initial, INITIAL);
    assert_eq!(u, input());
    emit("refusal-keeps-state-and-buffers", c, "refuse");
}
