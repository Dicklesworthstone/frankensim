//! G1/G3: positive, independently integrated physical heat at very small motion.
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedCommonFlowZoh, PreparedExactZoh, integrated_resistance_loss,
    split_quadratic, step_render,
};
use fs_phs::{
    AcousticSection, ChainOutlet, PortHamiltonian, QuadraticStorage, RelaxationBranch,
    ViscothermalPin, acoustic_chain_outlet_wall, join_port, mass_spring_damper, transformer,
};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && actual >= 0.0,
        "physical heat must be finite and nonnegative: {actual:e}"
    );
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "physical heat {actual:.17e} vs {expected:.17e}"
    );
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = duration / steps as f64;
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

#[test]
fn held_equilibrium_is_silent_and_tiny_vibration_heat_remains_resolved() {
    for damping in [1e-200, 1e-100, 1e-8, 0.02, 2.0, 20.0] {
        let sys = mass_spring_damper(1.0, 1.0, damping).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let prepared = PreparedExactZoh::new(&sys, 1.0, budget(2)).unwrap();
        assert_eq!(prepared.step(&[1.0, 0.0], &[1.0]).unwrap().dissipated, 0.0);
        assert_eq!(
            step_render(&sys, &split, &[1.0, 0.0], &[1.0], 1.0)
                .unwrap()
                .dissipated,
            0.0
        );
        // Integrate the independently written physical perturbation equation.
        // Its true q perturbation is the represented subtraction, not 1e-10.
        let q = 1.0 + 1e-10;
        let velocity = 0.0;
        let reference = rk4([q - 1.0, velocity, 0.0], 1.0, 4096, |x| {
            [x[1], -x[0] - damping * x[1], damping * x[1] * x[1]]
        });
        let direct = step_render(&sys, &split, &[q, velocity], &[1.0], 1.0).unwrap();
        let cached = prepared.step(&[q, velocity], &[1.0]).unwrap();
        relative(direct.dissipated, reference[2], 3e-11);
        assert_eq!(direct.dissipated.to_bits(), cached.dissipated.to_bits());
    }
}

#[test]
fn high_quality_factor_and_stiff_relaxation_keep_independent_positive_heat() {
    // One full lossless-period limit: integral sin(t)^2 dt = pi. The energy
    // subtraction is zero at this tiny damping, but the physical heat is not.
    for damping in [1e-24, 1e-100, 1e-200] {
        let sys = mass_spring_damper(1.0, 1.0, damping).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let step = step_render(&sys, &split, &[0.1, 0.0], &[0.0], core::f64::consts::TAU).unwrap();
        relative(
            step.dissipated,
            damping * 0.01 * core::f64::consts::PI,
            2e-12,
        );
    }
    for rate in [1e-12, 0.25, 1.0, 100.0, 1e6] {
        let sys = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![rate],
            vec![0.0],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap();
        let split = split_quadratic(&sys, budget(1)).unwrap();
        let prepared = PreparedExactZoh::new(&sys, 1.0, budget(1)).unwrap();
        let expected = -0.5 * (-2.0 * rate).exp_m1();
        relative(
            prepared.step(&[1.0], &[0.0]).unwrap().dissipated,
            expected,
            8e-12,
        );
        relative(
            step_render(&sys, &split, &[1.0], &[0.0], 1.0)
                .unwrap()
                .dissipated,
            expected,
            8e-12,
        );
        relative(
            integrated_resistance_loss(&sys, &split, &[1.0], &[0.0], 1.0, &[0.25 * rate]).unwrap(),
            0.25 * expected,
            8e-12,
        );
    }
}

fn relative_damper() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    r[5] = 0.3;
    r[7] = -0.3;
    r[13] = -0.3;
    r[15] = 0.3;
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 1.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
}

#[test]
fn a_damping_nullspace_and_its_small_relative_motion_are_not_replaced_by_heat_noise() {
    let sys = relative_damper();
    let prepared = PreparedExactZoh::new(&sys, 0.3, budget(4)).unwrap();
    let initial = [0.4, -0.3, 0.4, -0.3];
    let common = prepared.step(&initial, &[0.7]).unwrap();
    assert!(common.dissipated >= 0.0 && common.dissipated < 1e-28);
    let mut moving = initial;
    moving[3] += 1e-10;
    let reference = rk4([0.0, moving[1] - moving[3], 0.0], 0.3, 2048, |x| {
        [x[1], -4.0 * x[0] - 0.6 * x[1], 0.3 * x[1] * x[1]]
    });
    relative(
        prepared.step(&moving, &[0.7]).unwrap().dissipated,
        reference[2],
        2e-6,
    );
    let lossless = mass_spring_damper(1.0, 4.0, 0.0).unwrap();
    assert_eq!(
        PreparedExactZoh::new(&lossless, 0.3, budget(2))
            .unwrap()
            .step(&[0.4, -0.3], &[0.7])
            .unwrap()
            .dissipated,
        0.0
    );
}

fn material_body() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    r[5] = 0.4;
    r[7] = -0.1;
    r[13] = -0.1;
    r[15] = 0.2;
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 1.0, 0.0, -0.3],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 4.0,
        relaxation_time_s: 0.15,
    }])
    .unwrap()
}

#[test]
fn complete_material_heat_matches_physical_rhs_and_clock_subdivision_with_atomic_outputs() {
    let sys = material_body();
    let initial = [0.02, -0.03, 0.01, 0.04, 0.003];
    let force = 0.7;
    let duration = 0.3;
    let reference = rk4(
        [
            initial[0], initial[1], initial[2], initial[3], initial[4], 0.0, 0.0,
        ],
        duration,
        4096,
        |x| {
            let arm = 2.0 * x[0] - x[2] - x[4];
            [
                x[1],
                -4.0 * x[0] - 0.4 * x[1] + 0.1 * x[3] - 2.0 * arm + force,
                x[3],
                -9.0 * x[2] + 0.1 * x[1] - 0.2 * x[3] + arm - 0.3 * force,
                arm / 0.15,
                0.4 * x[1] * x[1] - 0.2 * x[1] * x[3] + 0.2 * x[3] * x[3] + arm * arm / 0.15,
                force * (x[1] - 0.3 * x[3]),
            ]
        },
    );
    let mut full = PreparedExactZoh::new(&sys, duration, budget(5)).unwrap();
    let record = full.step(&initial, &[force]).unwrap();
    for (actual, expected) in record.x.iter().zip(reference) {
        assert!((actual - expected).abs() < 2e-12, "{actual} vs {expected}");
    }
    relative(record.dissipated, reference[5], 3e-11);
    assert!((record.supplied - reference[6]).abs() < 2e-12);
    for count in [2, 4] {
        let sub = PreparedExactZoh::new(&sys, duration / count as f64, budget(5)).unwrap();
        let mut x = initial.to_vec();
        let mut heat = 0.0;
        for _ in 0..count {
            let r = sub.step(&x, &[force]).unwrap();
            x = r.x;
            heat += r.dissipated;
        }
        relative(heat, record.dissipated, 3e-12);
        for (a, b) in x.iter().zip(&record.x) {
            assert!((a - b).abs() < 2e-12);
        }
    }
    let mut next = vec![123.0; 5];
    let mut output = vec![456.0];
    assert!(
        full.step_into(&initial, &[f64::NAN], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, vec![123.0; 5]);
    assert_eq!(output, vec![456.0]);
    let accepted = full
        .step_into(&initial, &[force], &mut next, &mut output)
        .unwrap();
    assert_eq!(next, record.x);
    assert_eq!(accepted.dissipated, record.dissipated);
}

#[test]
fn separate_general_graph_and_tiny_loss_components_keep_their_physical_scales() {
    let n = 7;
    let mut r = vec![0.0; n * n];
    // Positive definite, but not diagonally dominant: needs the general factor.
    r[0] = 1.0;
    r[1] = 2.0;
    r[n] = 2.0;
    r[n + 1] = 5.0;
    // Six edge rows in four coordinates exercise the existing QR compression.
    for i in 2..6 {
        for j in i + 1..6 {
            r[i * n + i] += 0.25;
            r[j * n + j] += 0.25;
            r[i * n + j] = -0.25;
            r[j * n + i] = -0.25;
        }
    }
    r[6 * n + 6] = 1e-200;
    let mut q = vec![0.0; n * n];
    for i in 0..n {
        q[i * n + i] = 1.0;
    }
    let sys = PortHamiltonian::new(
        n,
        1,
        vec![0.0; n * n],
        r.clone(),
        vec![0.0; n],
        Box::new(QuadraticStorage::new(q, n).unwrap()),
    )
    .unwrap();
    let dt = 0.25;
    let prepared = PreparedExactZoh::new(&sys, dt, budget(n)).unwrap();
    assert_eq!(sys.structure().1, r);
    let mut initial = [0.0; 7];
    initial[6] = 1.0;
    relative(
        prepared.step(&initial, &[0.0]).unwrap().dissipated,
        dt * 1e-200,
        5e-12,
    );
    initial[6] = 0.0;
    initial[2] = 1.0;
    initial[3] = -1.0;
    relative(
        prepared.step(&initial, &[0.0]).unwrap().dissipated,
        -(-2.0 * dt).exp_m1(),
        5e-12,
    );
    initial[2..6].fill(1.0);
    assert!(prepared.step(&initial, &[0.0]).unwrap().dissipated < 1e-28);
    // Also check the admitted non-diagonally-dominant component's heat against
    // independently integrated coupled relaxation, not an energy difference.
    initial.fill(0.0);
    initial[0] = 0.3;
    initial[1] = -0.1;
    let reference = rk4([0.3, -0.1, 0.0], dt, 2048, |x| {
        [
            -x[0] - 2.0 * x[1],
            -2.0 * x[0] - 5.0 * x[1],
            x[0] * x[0] + 4.0 * x[0] * x[1] + 5.0 * x[1] * x[1],
        ]
    });
    relative(
        prepared.step(&initial, &[0.0]).unwrap().dissipated,
        reference[2],
        5e-12,
    );
}

#[test]
fn a_negative_ground_hidden_below_the_diagonal_ulp_is_not_admitted_as_passive() {
    let tiny = 2.0_f64.powi(-55);
    // The rounded row-0 sum is 1, but its exact ground is -tiny. At [1,1,1]
    // this represented matrix has negative physical power. An eigensolver's
    // stopping tolerance must not erase the weak off-diagonal and admit it.
    for (n, r) in [
        (3, vec![1.0, -1.0, -tiny, -1.0, 1.0, 0.0, -tiny, 0.0, tiny]),
        (2, vec![1.0, -1.0, -1.0, 1.0 - 1e-13]),
    ] {
        let mut q = vec![0.0; n * n];
        for i in 0..n {
            q[i * n + i] = 1.0;
        }
        let sys = PortHamiltonian::new(
            n,
            1,
            vec![0.0; n * n],
            r.clone(),
            vec![0.0; n],
            Box::new(QuadraticStorage::new(q, n).unwrap()),
        )
        .unwrap();
        assert!(PreparedExactZoh::new(&sys, 0.1, budget(n)).is_err());
        assert_eq!(sys.structure().1, r);
    }
}

#[test]
fn actual_foster_duct_heat_matches_edge_power_quadrature_and_subdivision() {
    // Evaluated dry air at 293.15 K and 101325 Pa. Keep the actual constructor:
    // its two thermal graph nullspaces previously produced rounded negative
    // eigenvalues beside the much stronger viscous series conductances.
    let density = 1.204_105_429_879_052_6;
    let sound_speed = 343.233_814_329_002;
    let viscosity = 1.813_405_882_148_823_4e-5;
    let pin = ViscothermalPin {
        dynamic_viscosity: viscosity,
        gamma: 1.4,
        prandtl: 0.708_416_410_255_627_1,
        foster_branches: 3,
    };
    let sys = acoustic_chain_outlet_wall(
        &[AcousticSection {
            length: 0.3,
            radius: 0.01,
            outlet_radius: 0.01,
            cells: 2,
        }],
        density,
        sound_speed,
        ChainOutlet::PressureRelease,
        1,
        &[],
        Some(&pin),
        None,
    )
    .unwrap();
    let n = sys.state_dim();
    assert_eq!(n, 16);
    let q = sys.storage().exact_quadratic_matrix().unwrap();
    let (j, r, g) = sys.structure();
    let original_r = r.to_vec();
    let initial_effort = [
        30.0, 1e-4, -10.0, -0.5e-4, 0.4e-4, 0.1e-4, -0.3e-4, 10.0, -2.0, 3.0, -0.1e-4, 0.5e-4,
        -0.6e-4, -1.0, 4.0, 2.0,
    ];
    let initial: [f64; 16] = core::array::from_fn(|i| initial_effort[i] / q[i * n + i]);
    let input = 1e-5;
    let duration = 1.0 / 48_000.0;
    let prepared = PreparedExactZoh::new(&sys, duration, budget(n)).unwrap();
    let record = prepared.step(&initial, &[input]).unwrap();
    assert_eq!(sys.structure().1, original_r);
    let mut extended = [0.0; 18];
    extended[..n].copy_from_slice(&initial);
    let inertance = density * 0.15 / (core::f64::consts::PI * 0.01_f64.powi(2));
    let steady_resistance = 8.0 * viscosity / (density * 0.01_f64.powi(2)) * inertance;
    let reference = rk4(extended, duration, 2048, |x| {
        let effort: [f64; 16] = core::array::from_fn(|i| q[i * n + i] * x[i]);
        let mut result = [0.0; 18];
        for i in 0..n {
            result[i] = g[i] * input
                + (0..n)
                    .map(|k| (j[i * n + k] - r[i * n + k]) * effort[k])
                    .sum::<f64>();
            result[17] += input * g[i] * effort[i];
            for k in i + 1..n {
                assert!(r[i * n + k] <= 0.0);
                result[16] += -r[i * n + k] * (effort[i] - effort[k]).powi(2);
            }
        }
        // Independently integrate the actual Foster edge powers and two
        // Poiseuille grounds; no endpoint energy subtraction or prepared root.
        result[16] += steady_resistance * (effort[1].powi(2) + effort[3].powi(2));
        result
    });
    relative(record.dissipated, reference[16], 3e-10);
    let energy = sys.hamiltonian(&initial);
    let state_error: f64 = (0..n)
        .map(|i| q[i * n + i] * (record.x[i] - reference[i]).powi(2))
        .sum();
    assert!(state_error < 1e-20 * energy);
    assert!((record.supplied - reference[17]).abs() < 2e-11 * energy);
    let closure = sys.hamiltonian(&record.x) - energy + record.dissipated - record.supplied;
    assert!(closure.abs() < 2e-11 * energy);
    let half = PreparedExactZoh::new(&sys, duration / 2.0, budget(n)).unwrap();
    let a = half.step(&initial, &[input]).unwrap();
    let b = half.step(&a.x, &[input]).unwrap();
    relative(a.dissipated + b.dissipated, record.dissipated, 3e-11);
    let subdivision_error: f64 = (0..n)
        .map(|i| q[i * n + i] * (b.x[i] - record.x[i]).powi(2))
        .sum();
    assert!(subdivision_error < 1e-20 * energy);
}

fn sea_level_foster(inlets: usize) -> PortHamiltonian {
    // Original dry-air source at 288.15 K and 101325 Pa, including its transport.
    acoustic_chain_outlet_wall(
        &[AcousticSection {
            length: 0.5,
            radius: 0.01,
            outlet_radius: 0.0,
            cells: 2,
        }],
        1.2249991558877122,
        340.2941077869353,
        ChainOutlet::PressureRelease,
        inlets,
        &[],
        Some(&ViscothermalPin {
            dynamic_viscosity: 1.7893802780775831e-5,
            gamma: 1.4,
            prandtl: 0.7098527545298142,
            foster_branches: 3,
        }),
        None,
    )
    .unwrap()
}

#[test]
fn source_foster_squares_survive_rounded_assembly_sharing_and_physical_composition() {
    let mut source = sea_level_foster(1);
    let n = source.state_dim();
    assert_eq!(n, 16);
    let q = source.storage().exact_quadratic_matrix().unwrap().to_vec();
    let (j, r, g) = source.structure();
    // The physical thermal star has no ground. Summing its positive resistor
    // edges rounded its diagonal down by less than an ulp. A raw matrix cannot
    // assert that physical provenance: its represented ground is negative.
    let (mut ground, mut tail) = (r[0], 0.0);
    for &edge in &r[1..n] {
        let (sum, error) = fs_math::eft::two_sum(ground, edge);
        ground = sum;
        tail += error;
    }
    assert_eq!(ground + tail, -4.1359030627651384e-25);
    let original_r = r.to_vec();
    let raw = PortHamiltonian::new(
        n,
        1,
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(QuadraticStorage::new(q.clone(), n).unwrap()),
    )
    .unwrap();
    let dt = 1.0 / 48_000.0;
    assert!(
        PreparedExactZoh::new(&raw, dt, budget(n)).is_err(),
        "source witnesses must never relax admission of a caller's raw negative ground"
    );
    let shared = source.share_storage();
    let initial_effort = [
        30.0, 1e-4, -10.0, -0.5e-4, 0.4e-4, 0.1e-4, -0.3e-4, 10.0, -2.0, 3.0, -0.1e-4, 0.5e-4,
        -0.6e-4, -1.0, 4.0, 2.0,
    ];
    let initial: [f64; 16] = core::array::from_fn(|i| initial_effort[i] / q[i * n + i]);
    let input = 1e-5;
    let record = PreparedExactZoh::new(&source, dt, budget(n))
        .unwrap()
        .step(&initial, &[input])
        .unwrap();
    let copy = PreparedExactZoh::new(&shared, dt, budget(n))
        .unwrap()
        .step(&initial, &[input])
        .unwrap();
    assert_eq!(record.x, copy.x);
    assert_eq!(record.dissipated.to_bits(), copy.dissipated.to_bits());
    assert_eq!(source.structure().1, original_r);
    let (j, r, g) = source.structure();
    let mut extended = [0.0; 18];
    extended[..n].copy_from_slice(&initial);
    let inertance = 1.2249991558877122 * 0.25 / (core::f64::consts::PI * 0.01_f64.powi(2));
    let ground_resistance =
        8.0 * 1.7893802780775831e-5 / (1.2249991558877122 * 0.01_f64.powi(2)) * inertance;
    let oracle = rk4(extended, dt, 2048, |x| {
        let e: [f64; 16] = core::array::from_fn(|i| q[i * n + i] * x[i]);
        let mut dx = [0.0; 18];
        for i in 0..n {
            dx[i] = g[i] * input
                + (0..n)
                    .map(|k| (j[i * n + k] - r[i * n + k]) * e[k])
                    .sum::<f64>();
            dx[17] += input * g[i] * e[i];
            for k in i + 1..n {
                dx[16] += -r[i * n + k] * (e[i] - e[k]).powi(2);
            }
        }
        dx[16] += ground_resistance * (e[1].powi(2) + e[3].powi(2));
        dx
    });
    relative(record.dissipated, oracle[16], 3e-10);
    let energy = source.hamiltonian(&initial);
    assert!((record.supplied - oracle[17]).abs() < 2e-11 * energy);
    let state_error: f64 = (0..n)
        .map(|i| q[i * n + i] * (record.x[i] - oracle[i]).powi(2))
        .sum();
    assert!(state_error < 1e-20 * energy);

    // Keep original full-C/Maxwell body, separate mechanical and acoustic
    // source ports, then add the common-velocity junction. Neither composition
    // may erase the source's physical resistor squares or any material state.
    let body = material_body()
        .with_additional_ports(1, &[0.0, 0.03, 0.0, -0.02, 0.0])
        .unwrap();
    let loaded = transformer(body, sea_level_foster(2), 1, 1, 1.0).unwrap();
    let joined = join_port(mass_spring_damper(2.0, 3.0, 0.1).unwrap(), loaded, 0, 0).unwrap();
    let n = joined.state_dim();
    assert_eq!(n, 24);
    let mut initial = vec![0.0; n];
    initial[2] = 0.01;
    initial[6] = 0.02;
    let prepared = PreparedCommonFlowZoh::new_owned(joined, dt, budget(n)).unwrap();
    let step = prepared.step(&initial, &[input]).unwrap();
    assert!(step.dissipated > 0.0);
    assert_eq!(step.x.len(), n);
    assert!(step.balance_residual().abs() < 1e-12 * prepared.hamiltonian(&initial));
}
