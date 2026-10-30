//! G1/G3: exact common-flow evolution keeps physical histories and accounts.
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedCommonFlowZoh, admit_midpoint_band, descriptor_pole_rate_bound,
};
use fs_phs::{
    AcousticSection, ChainOutlet, DescriptorPortHamiltonian, PhsError, PortHamiltonian,
    QuadraticStorage, RelaxationBranch, ViscothermalPin, acoustic_chain_outlet_wall,
    duffing_oscillator, join_port, mass_spring_damper,
};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

fn mechanical_join() -> DescriptorPortHamiltonian {
    let left = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.1],
        vec![0.0, 0.0, 1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![3.0, 0.0, 0.0, 0.5], 2).unwrap()),
    )
    .unwrap();
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    // Positive definite material C, whose first row is not diagonally dominant.
    r[5] = 0.4;
    r[7] = 0.5;
    r[13] = 0.5;
    r[15] = 0.8;
    let right = PortHamiltonian::new(
        4,
        3,
        j,
        r,
        vec![
            0.0, 0.0, 0.0, 1.0, 0.7, -0.2, 0.0, 0.0, 0.0, -0.3, -0.4, 1.1,
        ],
        Box::new(
            QuadraticStorage::new(
                vec![
                    4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                ],
                4,
            )
            .unwrap(),
        ),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 4.0,
        relaxation_time_s: 0.15,
    }])
    .unwrap();
    join_port(left, right, 0, 0).unwrap()
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    time: f64,
    steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = time / steps as f64;
    for _ in 0..steps {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..N {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

// Independent Newton equations in original displacement/momentum/Maxwell
// coordinates. Differentiating va = vb - 0.3 vc gives the scalar reaction;
// this oracle never constructs the production projection or exponential.
fn mechanical_rhs(x: &[f64; 12], u: &[f64; 3]) -> ([f64; 12], f64) {
    let va = x[1] / 2.0;
    let (vb, vc) = (x[3], x[5]);
    let arm = 2.0 * x[2] - x[4] - x[6];
    let fa = -3.0 * x[0] - 0.1 * va + u[0];
    let fb = -4.0 * x[2] - 0.4 * vb - 0.5 * vc - 2.0 * arm + 0.7 * u[1] - 0.2 * u[2];
    let fc = -9.0 * x[4] - 0.5 * vb - 0.8 * vc + arm - 0.4 * u[1] + 1.1 * u[2];
    let reaction = -(fa / 2.0 - fb + 0.3 * fc) / 1.59;
    let y = [va, 0.7 * vb - 0.4 * vc, -0.2 * vb + 1.1 * vc];
    let heat = 0.1 * va * va + 0.4 * vb * vb + vb * vc + 0.8 * vc * vc + arm * arm / 0.15;
    (
        [
            va,
            fa + reaction,
            vb,
            fb - reaction,
            vc,
            fc + 0.3 * reaction,
            arm / 0.15,
            y[0],
            y[1],
            y[2],
            heat,
            u[0] * y[0] + u[1] * y[1] + u[2] * y[2],
        ],
        reaction,
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn full_material_constraint_matches_independent_motion_heat_work_and_reaction() {
    let initial = [0.04, 0.2, 0.02, 0.112, -0.03, 0.04, 0.03, -123.0];
    let input = [0.3, 0.07, -0.2];
    let dt = 0.037;
    let mut runtime = PreparedCommonFlowZoh::new_owned(mechanical_join(), dt, budget(8)).unwrap();
    assert_eq!(
        (
            runtime.state_dim(),
            runtime.differential_dim(),
            runtime.port_dim()
        ),
        (8, 7, 3)
    );
    assert_eq!(runtime.duration(), dt);
    let mut oracle = [0.0; 12];
    oracle[..7].copy_from_slice(&initial[..7]);
    let initial_reaction = mechanical_rhs(&oracle, &input).1;
    close(
        runtime.instantaneous_reaction(&initial, &input).unwrap(),
        initial_reaction,
        2e-14,
    );
    let reference = rk4(oracle, dt, 4096, |x| mechanical_rhs(x, &input).0);
    let step = runtime.step(&initial, &input).unwrap();
    for (&actual, &expected) in step.x[..7].iter().zip(&reference[..7]) {
        close(actual, expected, 2e-12);
    }
    for i in 0..3 {
        close(step.y[i], reference[7 + i] / dt, 2e-12);
    }
    close(step.dissipated, reference[10], 2e-13);
    close(step.supplied, reference[11], 2e-13);
    close(step.x[7], mechanical_rhs(&reference, &input).1, 2e-12);
    close(step.x[1] / 2.0, step.x[3] - 0.3 * step.x[5], 2e-14);
    assert!(step.dissipated > 0.0);
    assert!(step.balance_residual().abs() < 2e-14);
    assert!(step.solver_certificate.is_none());
    close(
        runtime.hamiltonian(&step.x) - runtime.hamiltonian(&initial),
        step.delta_h,
        2e-14,
    );

    // The algebraic seed owns no energy and cannot alter a single physical bit.
    let mut alternate = initial;
    alternate[7] = 1e80;
    assert_eq!(runtime.step(&alternate, &input).unwrap().x, step.x);
    let mut next = [17.0; 8];
    let mut average = [19.0; 3];
    let record = runtime
        .step_into(&initial, &input, &mut next, &mut average)
        .unwrap();
    assert_eq!(next.as_slice(), step.x);
    assert_eq!(average.as_slice(), step.y);
    assert_eq!(record.dissipated.to_bits(), step.dissipated.to_bits());

    let half = PreparedCommonFlowZoh::new_owned(mechanical_join(), dt / 2.0, budget(8)).unwrap();
    let first = half.step(&initial, &input).unwrap();
    let second = half.step(&first.x, &input).unwrap();
    for (&actual, &expected) in second.x.iter().zip(&step.x) {
        close(actual, expected, 3e-13);
    }
    for i in 0..3 {
        close((first.y[i] + second.y[i]) / 2.0, step.y[i], 2e-13);
    }
    close(first.dissipated + second.dissipated, step.dissipated, 2e-14);
    close(first.supplied + second.supplied, step.supplied, 2e-14);
}

fn acoustic_join() -> DescriptorPortHamiltonian {
    let sections = [AcousticSection {
        length: 0.3,
        radius: 0.01,
        outlet_radius: 0.01,
        cells: 2,
    }];
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8134058821488234e-5,
        gamma: 1.4,
        prandtl: 0.7084164102556271,
        foster_branches: 3,
    };
    let duct = acoustic_chain_outlet_wall(
        &sections,
        1.204,
        343.2,
        ChainOutlet::PressureRelease,
        2,
        &[],
        Some(&pin),
        None,
    )
    .unwrap();
    let chamber = acoustic_chain_outlet_wall(
        &sections,
        1.204,
        343.2,
        ChainOutlet::Rigid,
        2,
        &[],
        None,
        None,
    )
    .unwrap();
    // In acoustic units the shared output is pressure and the split inputs
    // are volume flow. The untouched far-face ports remain external.
    join_port(duct, chamber, 0, 0).unwrap()
}

#[test]
fn retained_foster_history_evolves_at_audio_clock_with_original_constraint() {
    let system = acoustic_join();
    let n = system.state_dim();
    assert_eq!(n, 21); // 16 original duct states + four chamber states + lambda.
    let dt = 1.0 / 48_000.0;
    let rate = descriptor_pole_rate_bound(&system).unwrap();
    assert!(admit_midpoint_band(rate, dt).is_err());
    let q = system.storage().quadratic_part().unwrap().matrix;
    let constraint = system.dirac_j()[20 * n..21 * n].to_vec();
    let mut initial = vec![0.0; n];
    initial[4] = 1e-3 / q[4 * n + 4].sqrt(); // Warm series Foster history.
    initial[7] = -0.4e-3 / q[7 * n + 7].sqrt(); // Warm thermal Foster history.
    let input = [1e-5, -2e-5];
    let energy = system.hamiltonian(&initial);
    let runtime = PreparedCommonFlowZoh::new_owned(system, dt, budget(n)).unwrap();
    let step = runtime.step(&initial, &input).unwrap();
    assert_eq!(step.x.len(), n);
    assert!(step.dissipated > 0.0);
    assert!((step.x[4] - initial[4]).abs() > 1e-6 * initial[4].abs());
    let mut effort = vec![0.0; n];
    for i in 0..n {
        effort[i] = (0..n).map(|j| q[i * n + j] * step.x[j]).sum();
    }
    let mismatch: f64 = constraint.iter().zip(&effort).map(|(a, b)| a * b).sum();
    let scale: f64 = constraint
        .iter()
        .zip(&effort)
        .map(|(a, b)| (a * b).abs())
        .sum();
    assert!(mismatch.abs() <= 2e-12 * scale);
    assert!(step.balance_residual().abs() < 2e-12 * energy.max(step.supplied.abs()));

    let half = PreparedCommonFlowZoh::new_owned(acoustic_join(), dt / 2.0, budget(n)).unwrap();
    let a = half.step(&initial, &input).unwrap();
    let b = half.step(&a.x, &input).unwrap();
    let weighted_error: f64 = (0..n)
        .map(|i| q[i * n + i] * (b.x[i] - step.x[i]).powi(2))
        .sum();
    assert!(weighted_error < 1e-20 * energy);
    close(a.dissipated + b.dissipated, step.dissipated, 2e-11 * energy);
    close(a.supplied + b.supplied, step.supplied, 2e-11 * energy);
}

#[test]
fn unsupported_or_incompatible_states_refuse_without_publishing_or_projecting() {
    let nonlinear = join_port(
        duffing_oscillator(1.0, 4.0, 0.1, 0.02).unwrap(),
        mass_spring_damper(1.0, 4.0, 0.03).unwrap(),
        0,
        0,
    )
    .unwrap();
    assert!(matches!(
        PreparedCommonFlowZoh::new_owned(nonlinear, 0.01, budget(5)),
        Err(PhsError::LinearSplitBudget { .. })
    ));
    assert!(PreparedCommonFlowZoh::new_owned(mechanical_join(), 0.01, budget(7)).is_err());
    assert!(PreparedCommonFlowZoh::new_owned(mechanical_join(), -0.01, budget(8)).is_err());
    let mut runtime = PreparedCommonFlowZoh::new_owned(mechanical_join(), 0.01, budget(8)).unwrap();
    let initial = [0.04, 0.2, 0.02, 0.112, -0.03, 0.04, 0.03, 0.0];
    let mut bad = initial;
    bad[1] += 0.01;
    let mut next = [17.0; 8];
    let mut output = [19.0; 3];
    assert!(
        runtime
            .step_into(&bad, &[0.0; 3], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, [17.0; 8]);
    assert_eq!(output, [19.0; 3]);
    assert_eq!(bad[1], initial[1] + 0.01);
    assert!(
        runtime
            .step_into(&initial, &[f64::NAN; 3], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, [17.0; 8]);
    assert_eq!(output, [19.0; 3]);
    assert!(
        runtime
            .step_into(&initial[..7], &[0.0; 3], &mut next, &mut output)
            .is_err()
    );
    assert!(
        runtime
            .step_into(&initial, &[0.0; 3], &mut next[..7], &mut output)
            .is_err()
    );
    assert_eq!(next, [17.0; 8]);
    assert_eq!(output, [19.0; 3]);
    let zero = PreparedCommonFlowZoh::new_owned(mechanical_join(), 0.0, budget(8)).unwrap();
    let record = zero.step(&initial, &[0.0; 3]).unwrap();
    assert_eq!(&record.x[..7], &initial[..7]);
    assert_eq!(
        (record.delta_h, record.dissipated, record.supplied),
        (0.0, 0.0, 0.0)
    );
}
