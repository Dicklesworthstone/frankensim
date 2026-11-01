//! G0/G2/G3: physical direct loads share the complete continuous trajectory.
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage, ResistanceWitness};

fn body(g: Vec<f64>, ports: usize) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        ports,
        vec![0.0; 4],
        vec![0.0; 4],
        g,
        Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap()
}
fn exact(system: PortHamiltonian, dt: f64) -> PreparedExactZoh<'static> {
    let n = system.state_dim();
    PreparedExactZoh::new_owned(
        system,
        dt,
        LinearSplitBudget {
            max_modes: n,
            max_omega_dt: 0.04,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-12 * scale, "{a:.17e} != {b:.17e}");
}

#[test]
fn mutual_load_retains_both_decay_rates_and_independent_heat() {
    let system = body(vec![1.0, 0.0, 0.0, 1.0], 2)
        .with_port_conductance(&[0, 1], &[2.0, -1.0, -1.0, 2.0])
        .unwrap();
    assert_eq!(system.state_dim(), 2);
    assert_eq!(system.port_dim(), 2);
    assert_eq!(system.structure().1, &[2.0, -1.0, -1.0, 2.0]);
    let dt: f64 = 0.37;
    let step = exact(system, dt).step(&[1.0, 0.0], &[0.0, 0.0]).unwrap();
    close(step.x[0], ((-dt).exp() + (-3.0 * dt).exp()) / 2.0, 1.0);
    close(step.x[1], ((-dt).exp() - (-3.0 * dt).exp()) / 2.0, 1.0);
    let heat = (1.0 - (-2.0 * dt).exp()) / 4.0 + (1.0 - (-6.0 * dt).exp()) / 4.0;
    close(step.dissipated, heat, 1.0);
    close(step.delta_h + step.dissipated, 0.0, 1.0);
}

#[test]
fn direct_and_original_source_squares_share_forced_dynamics() {
    let source = [0.3, -0.6, 0.15];
    let r: Vec<_> = source
        .iter()
        .flat_map(|a| source.iter().map(move |b| a * b))
        .collect();
    let witness = ResistanceWitness::try_from_source_rows(3, r.clone(), source.to_vec()).unwrap();
    // Original outputs: y0=x0+x2, y1=x1, y2=2*x0-x1.
    let system = PortHamiltonian::new(
        3,
        3,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.5, 0.0, -0.5, 0.0],
        r,
        vec![1.0, 0.0, 2.0, 0.0, 1.0, -1.0, 1.0, 0.0, 0.0],
        Box::new(
            QuadraticStorage::new(vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 3).unwrap(),
        ),
    )
    .unwrap()
    .with_resistance_witness(witness)
    .unwrap();
    let row = [0.4, -0.7];
    let d: Vec<_> = row
        .iter()
        .flat_map(|a| row.iter().map(move |b| a * b))
        .collect();
    let load = ResistanceWitness::try_from_source_rows(2, d, row.to_vec()).unwrap();
    let original_h = system.hamiltonian(&[0.2, -0.3, 0.1]);
    let system = system.with_port_conductance_witness(&[2, 0], load).unwrap();
    assert_eq!(system.hamiltonian(&[0.2, -0.3, 0.1]), original_h);
    let dt = 0.4;
    let step = exact(system, dt)
        .step(&[0.2, -0.3, 0.1], &[0.1, 0.2, -0.1])
        .unwrap();
    // Independent continuous equations: original gyroscopic coupling,
    // source resistor, and the physical rank-one load in output coordinates.
    let rhs = |x: &[f64; 5]| {
        let y = [x[0] + x[2], x[1], 2.0 * x[0] - x[1]];
        let source_effort = 0.3 * x[0] - 0.6 * x[1] + 0.15 * x[2];
        let load_effort = 0.4 * y[2] - 0.7 * y[0];
        let u = [0.1 + 0.7 * load_effort, 0.2, -0.1 - 0.4 * load_effort];
        [
            x[1] - 0.3 * source_effort + u[0] + 2.0 * u[2],
            -x[0] + 0.5 * x[2] + 0.6 * source_effort + u[1] - u[2],
            -0.5 * x[1] - 0.15 * source_effort + u[0],
            source_effort * source_effort + load_effort * load_effort,
            0.1 * y[0] + 0.2 * y[1] - 0.1 * y[2],
        ]
    };
    let mut x = [0.2, -0.3, 0.1, 0.0, 0.0];
    let h = dt / 2048.0;
    for _ in 0..2048 {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..5 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    for i in 0..3 {
        close(step.x[i], x[i], 1.0);
    }
    close(step.dissipated, x[3], 1.0);
    close(step.supplied, x[4], 1.0);
}

#[test]
fn physical_port_scaling_preserves_small_and_large_conductance() {
    let system = body(vec![1e10, 0.0, 0.0, 1e-10], 2)
        .with_port_conductance(&[0, 1], &[1e-20, 0.0, 0.0, 1e20])
        .unwrap();
    for (&a, &b) in system.structure().1.iter().zip(&[1.0, 0.0, 0.0, 1.0]) {
        close(a, b, 1.0);
    }
    let y = system.output(&[0.4, -0.7]);
    // Compare the physical port law directly across forty decades in D.
    // Temporal admission remains the exact propagator's separate obligation.
    close(1e-20 * y[0] * y[0] + 1e20 * y[1] * y[1], 0.65, 1.0);
}

#[test]
fn selected_order_and_inactive_ports_preserve_the_physical_load() {
    let g = vec![1.0, 0.0, 2.0, 0.0, 1.0, -1.0];
    let a = body(g.clone(), 3)
        .with_port_conductance(&[2, 0], &[2.0, 0.5, 0.5, 1.0])
        .unwrap();
    let b = body(g.clone(), 3)
        .with_port_conductance(&[0, 2], &[1.0, 0.5, 0.5, 2.0])
        .unwrap();
    assert_eq!(a.structure(), b.structure());
    let step_a = exact(a, 0.1).step(&[0.3, -0.1], &[0.1, 0.2, 0.3]).unwrap();
    let step_b = exact(b, 0.1).step(&[0.3, -0.1], &[0.1, 0.2, 0.3]).unwrap();
    assert_eq!(step_a.x, step_b.x);
    let identity = body(g.clone(), 3)
        .with_port_conductance(&[], &[])
        .unwrap()
        .with_port_conductance(&[2, 1, 0], &[0.0; 9])
        .unwrap();
    assert_eq!(identity.structure().2, g);
    assert_eq!(identity.structure().1, &[0.0; 4]);
}

#[test]
fn invalid_or_nonpassive_loads_refuse_without_repair() {
    let make = || body(vec![1.0, 0.0, 0.0, 1.0], 2);
    for (ports, d) in [
        (vec![0, 0], vec![1.0, 0.0, 0.0, 1.0]),
        (vec![2], vec![1.0]),
        (vec![0, 1], vec![1.0]),
        (vec![0], vec![f64::NAN]),
        (vec![0], vec![f64::INFINITY]),
        (vec![0], vec![-1e-200]),
        (vec![0, 1], vec![1.0, 2.0, 2.0, 1.0]),
        (vec![0, 1], vec![1.0, 0.0, 0.1, 1.0]),
    ] {
        assert!(make().with_port_conductance(&ports, &d).is_err());
    }
}
