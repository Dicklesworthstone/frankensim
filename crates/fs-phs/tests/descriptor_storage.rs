//! Coupled material storage uses the existing descriptor structure and solver.

use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::{
    DescriptorPortHamiltonian, PhsError, QuadraticStorage, Storage, discrete_gradient,
    step_descriptor,
};
use std::{cell::Cell, rc::Rc};

struct CountedQuadratic {
    inner: QuadraticStorage,
    declare_exact: bool,
    gradients: Rc<Cell<usize>>,
}

impl Storage for CountedQuadratic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.inner.hamiltonian(x)
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.gradients.set(self.gradients.get() + 1);
        self.inner.gradient(x, out);
    }

    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.declare_exact.then(|| self.inner.matrix())
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        self.inner.midpoint_energy_remainder(a, b, mid)
    }
}

fn coupled_rlc(scale: f64, declare_exact: bool) -> (DescriptorPortHamiltonian, Rc<Cell<usize>>) {
    use fs_phs::circuit::{MutualInductance, assemble_circuit_with_mutual_inductance};
    let circuit = assemble_circuit_with_mutual_inductance(
        &CircuitGraph {
            node_count: 4,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Inductor {
                        henries: 2.0 * scale,
                    },
                ),
                (1, 0, Branch::Resistor { ohms: 1.0 }),
                (
                    2,
                    0,
                    Branch::Inductor {
                        henries: 0.5 * scale,
                    },
                ),
                (2, 0, Branch::Resistor { ohms: 2.0 }),
                (
                    3,
                    0,
                    Branch::Capacitor {
                        farads: 4.0 * scale,
                    },
                ),
                (3, 0, Branch::Resistor { ohms: 3.0 }),
                (1, 0, Branch::CurrentSource { port: 0 }),
            ],
            transformers: vec![],
        },
        &[MutualInductance {
            branches: (0, 2),
            henries: 0.5 * scale,
        }],
    )
    .unwrap();
    let native = &circuit.system;
    let n = native.state_dim();
    let nd = native.differential_dim();
    let q = native.storage().exact_quadratic_matrix().unwrap();
    let mut differential = vec![0.0; nd * nd];
    for row in 0..nd {
        differential[row * nd..(row + 1) * nd].copy_from_slice(&q[row * n..row * n + nd]);
    }
    let gradients = Rc::new(Cell::new(0));
    let (j, r, g) = native.structure();
    let system = DescriptorPortHamiltonian::new(
        n,
        nd,
        native.port_dim(),
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(CountedQuadratic {
            inner: QuadraticStorage::new(differential, nd).unwrap(),
            declare_exact,
            gradients: Rc::clone(&gradients),
        }),
    )
    .unwrap();
    (system, gradients)
}

#[test]
fn exact_quadratic_coupled_rlc_uses_native_ledger_without_coordinate_probes() {
    // The same physical trajectory under fifteen orders of magnitude of L/C/time
    // scaling. An independent two-winding KVL solve includes both M terms.
    for scale in [1e-9, 1.0, 1e6] {
        let (system, gradients) = coupled_rlc(scale, true);
        let initial = [1.8 * scale, 0.3 * scale, 0.8 * scale, 0.0, 0.0, 0.0, 0.0];
        let record = step_descriptor(&system, &initial, &[0.2], 0.1 * scale).unwrap();
        let determinant = 2.05 * 0.6 - 0.25;
        let currents = [
            (0.6 * 1.77 - 0.5 * 0.34) / determinant,
            (2.05 * 0.34 - 0.5 * 1.77) / determinant,
        ];
        let charge = 0.8 * (1.0 - 0.1 / 24.0) / (1.0 + 0.1 / 24.0);
        let expected = [
            2.0 * currents[0] + 0.5 * currents[1],
            0.5 * currents[0] + 0.5 * currents[1],
            charge,
        ];
        for (actual, value) in record.x[..3].iter().zip(expected) {
            assert!((actual / scale - value).abs() < 2e-12);
        }
        let mid_i = [
            f64::midpoint(1.0, currents[0]),
            f64::midpoint(-0.4, currents[1]),
        ];
        let voltages = [
            0.2 - mid_i[0],
            -2.0 * mid_i[1],
            f64::midpoint(0.8, charge) / 4.0,
        ];
        let heat =
            0.1 * (voltages[0].powi(2) + voltages[1].powi(2) / 2.0 + voltages[2].powi(2) / 3.0);
        assert!((record.dissipated / scale - heat).abs() < 2e-12);
        assert!((record.supplied / scale - 0.1 * 0.2 * voltages[0]).abs() < 2e-12);
        assert!((record.balance_residual() / scale).abs() < 3e-12);
        assert!(record.newton_iters <= 3);
        assert!(
            gradients.get() < 2 * system.state_dim(),
            "no centered gradient column sweep"
        );
    }
}

#[test]
fn exact_quadratic_descriptor_matches_the_undeclared_storage_fallback() {
    let (exact, exact_calls) = coupled_rlc(1.0, true);
    let (probed, probed_calls) = coupled_rlc(1.0, false);
    let mut a = vec![1.8, 0.3, 0.8, 0.0, 0.0, 0.0, 0.0];
    let mut b = a.clone();
    for input in [0.2, -0.7, 1.1, 0.0] {
        let fast = step_descriptor(&exact, &a, &[input], 0.1).unwrap();
        let reference = step_descriptor(&probed, &b, &[input], 0.1).unwrap();
        for (x, y) in fast.x.iter().zip(&reference.x) {
            assert!((x - y).abs() < 2e-11);
        }
        for (x, y) in [fast.dissipated, fast.supplied, fast.delta_h, fast.y[0]]
            .into_iter()
            .zip([
                reference.dissipated,
                reference.supplied,
                reference.delta_h,
                reference.y[0],
            ])
        {
            assert!((x - y).abs() < 2e-11);
        }
        a = fast.x;
        b = reference.x;
    }
    assert!(probed_calls.get() >= exact_calls.get() + 8 * exact.state_dim());
}

struct QuarticCapacitor;

impl Storage for QuarticCapacitor {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        assert_eq!(x.len(), 1, "the law must never receive a multiplier");
        x[0].powi(4) / 4.0
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        assert_eq!((x.len(), out.len()), (1, 1));
        out[0] = x[0].powi(3);
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        assert_eq!((a.len(), b.len(), mid.len()), (1, 1, 1));
        Some(mid[0] * (b[0] - a[0]).powi(3) / 4.0)
    }
}

fn nonlinear_rc(current_scale: f64) -> DescriptorPortHamiltonian {
    // q_dot = i, 0 = voltage_source - q^3 - i. The multiplier is
    // current_scale*i; changing its coordinates must not change physics.
    DescriptorPortHamiltonian::new(
        2,
        1,
        1,
        vec![0.0, 1.0 / current_scale, -1.0 / current_scale, 0.0],
        vec![0.0, 0.0, 0.0, 1.0 / current_scale.powi(2)],
        vec![0.0, 1.0 / current_scale],
        Box::new(QuarticCapacitor),
    )
    .unwrap()
}

#[test]
fn nonlinear_descriptor_keeps_the_complete_energy_remainder_on_charge() {
    let q0 = 0.5_f64;
    let voltage = 2.0;
    let dt = 0.3;
    // Independent scalar reference: the quartic energy divided difference
    // is (q1^3 + q1^2*q0 + q1*q0^2 + q0^3)/4, not midpoint q^3.
    let mut lower = q0;
    let mut upper = 2.0;
    for _ in 0..60 {
        let q1 = f64::midpoint(lower, upper);
        let effort = (q1.powi(3) + q1.powi(2) * q0 + q1 * q0.powi(2) + q0.powi(3)) / 4.0;
        if q1 - q0 - dt * (voltage - effort) > 0.0 {
            upper = q1;
        } else {
            lower = q1;
        }
    }
    let expected_charge = f64::midpoint(lower, upper);
    let expected_current = (expected_charge - q0) / dt;
    for current_scale in [1.0, 8.0] {
        let system = nonlinear_rc(current_scale);
        let before = [q0, current_scale * (voltage - q0.powi(3))];
        let record = step_descriptor(&system, &before, &[voltage], dt).unwrap();
        assert!((record.x[0] - expected_charge).abs() < 2e-12);
        assert!((record.y[0] - expected_current).abs() < 2e-12);
        assert!((record.dissipated - dt * expected_current.powi(2)).abs() < 5e-12);
        assert!((record.supplied - dt * voltage * expected_current).abs() < 5e-12);
        assert!((record.delta_h + record.dissipated - record.supplied).abs() < 5e-12);

        // Reproduce the old defect independently: forming Gonzalez over charge
        // AND the changing current multiplier assigns energy to the latter.
        // Replacing that algebraic effort then discards this nonzero remainder.
        let wrong = discrete_gradient(system.storage(), &before, &record.x);
        let lost = record.delta_h - wrong[0] * (record.x[0] - q0);
        assert!(lost.abs() > 1e-4);
    }
}

#[test]
fn storage_is_differential_only_and_quadratic_metadata_has_zero_algebraic_block() {
    let nonlinear = nonlinear_rc(1.0);
    assert_eq!(nonlinear.hamiltonian(&[2.0, 1e150]), 4.0);
    let mut gradient = [7.0, 7.0];
    nonlinear.storage().gradient(&[2.0, -1e150], &mut gradient);
    assert_eq!(gradient, [8.0, 0.0]);
    let quadratic = DescriptorPortHamiltonian::new(
        2,
        1,
        0,
        vec![0.0; 4],
        vec![0.0; 4],
        vec![],
        Box::new(QuadraticStorage::new(vec![3.0], 1).unwrap()),
    )
    .unwrap();
    let expected = [3.0, 0.0, 0.0, 0.0];
    assert_eq!(
        quadratic.storage().exact_quadratic_matrix().unwrap(),
        &expected
    );
    let part = quadratic.storage().quadratic_part().unwrap();
    assert!(part.exact);
    assert_eq!(part.matrix, expected);
    assert_eq!(quadratic.hamiltonian(&[2.0, 1e150]), 6.0);
}

#[test]
fn reusing_native_rlc_stamps_and_storage_preserves_ports_and_trajectory() {
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (2, 0, Branch::Capacitor { farads: 3.0 }),
            (1, 2, Branch::Resistor { ohms: 4.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let native = &circuit.system;
    let n = native.state_dim();
    let nd = native.differential_dim();
    let q = native.storage().exact_quadratic_matrix().unwrap();
    let mut differential_q = vec![0.0; nd * nd];
    for row in 0..nd {
        differential_q[row * nd..(row + 1) * nd].copy_from_slice(&q[row * n..row * n + nd]);
    }
    let (j, r, g) = native.structure();
    let rebuilt = DescriptorPortHamiltonian::new(
        n,
        nd,
        native.port_dim(),
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(QuadraticStorage::new(differential_q, nd).unwrap()),
    )
    .unwrap();
    assert_eq!(native.structure(), rebuilt.structure());
    let ports = [0.15];
    let mut state = vec![0.0; n];
    state[circuit.flux_index[0]] = 0.25;
    state[circuit.charge_index[0]] = 0.4;
    let initial = circuit.consistent_initial_state(&state, &ports).unwrap();
    let mut a = initial.clone();
    let mut b = initial;
    for _ in 0..8 {
        let original = circuit.step_audited(&a, &ports, 0.1).unwrap().0;
        let coupled = step_descriptor(&rebuilt, &b, &ports, 0.1).unwrap();
        for (left, right) in original.x.iter().zip(&coupled.x) {
            assert!((left - right).abs() < 2e-12);
        }
        assert!((original.y[0] - coupled.y[0]).abs() < 2e-12);
        assert!((original.dissipated - coupled.dissipated).abs() < 2e-12);
        assert!((coupled.delta_h + coupled.dissipated - coupled.supplied).abs() < 2e-12);
        a = original.x;
        b = coupled.x;
    }
}

#[test]
fn descriptor_constructor_rejects_invalid_partition_and_structure() {
    let make =
        |nd, j, r, g| DescriptorPortHamiltonian::new(2, nd, 1, j, r, g, Box::new(QuarticCapacitor));
    assert!(make(3, vec![0.0; 4], vec![0.0; 4], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0, 1.0, 1.0, 0.0], vec![0.0; 4], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0; 4], vec![0.0, 0.0, 0.0, -1.0], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0; 4], vec![0.0; 4], vec![f64::NAN, 0.0]).is_err());
    assert!(
        DescriptorPortHamiltonian::new(
            usize::MAX,
            1,
            1,
            vec![],
            vec![],
            vec![],
            Box::new(QuarticCapacitor),
        )
        .is_err()
    );
}

struct BoundedPlateStorage {
    epsilon_area: f64,
    mass: f64,
    gap_limits: [f64; 2],
}

impl BoundedPlateStorage {
    fn supported(&self, x: &[f64]) -> bool {
        self.gap_limits[0] < x[1] && x[1] < self.gap_limits[1]
    }
}

impl Storage for BoundedPlateStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !self.supported(x) {
            return f64::NAN;
        }
        0.5 * x[0] * (x[0] / self.epsilon_area) * x[1] + 0.5 * x[2] * (x[2] / self.mass)
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !self.supported(x) {
            out.fill(f64::NAN);
            return;
        }
        out[0] = (x[0] / self.epsilon_area) * x[1];
        out[1] = 0.5 * x[0] * (x[0] / self.epsilon_area);
        out[2] = x[2] / self.mass;
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], _: &[f64]) -> Option<f64> {
        if !self.supported(a) || !self.supported(b) {
            return Some(f64::NAN);
        }
        let dq = b[0] - a[0];
        Some(0.125 * dq * (dq / self.epsilon_area) * (b[1] - a[1]))
    }
}

const PLATE_EPSILON_AREA: f64 = 8.854_187_812_8e-12 * 0.01;
const PLATE_MASS: f64 = 0.1;

fn bounded_plate() -> DescriptorPortHamiltonian {
    // [isolated charge, gap, momentum, source voltage, source current].
    // The last pair imposes voltage=u and current=0. This unloaded source
    // changes the global Newton scale without adding charge or energy.
    let mut j = vec![0.0; 25];
    j[7] = 1.0;
    j[11] = -1.0;
    j[19] = 1.0;
    j[23] = -1.0;
    DescriptorPortHamiltonian::new(
        5,
        3,
        1,
        j,
        vec![0.0; 25],
        vec![0.0, 0.0, 0.0, 0.0, 1.0],
        Box::new(BoundedPlateStorage {
            epsilon_area: PLATE_EPSILON_AREA,
            mass: PLATE_MASS,
            gap_limits: [5e-5, 2e-4],
        }),
    )
    .unwrap()
}

#[test]
fn narrow_physical_gap_with_high_voltage_uses_finite_centered_probes() {
    let gap = 1e-4;
    let voltage = 100.0;
    let charge = PLATE_EPSILON_AREA * voltage / gap;
    let dt = 1e-3;
    let system = bounded_plate();
    let before = [charge, gap, 0.0, 0.0, 0.0];

    // At fixed isolated charge, plate attraction is constant, so midpoint
    // exactly gives the elementary constant-acceleration trajectory.
    let attraction = charge.powi(2) / (2.0 * PLATE_EPSILON_AREA);
    let expected_momentum = -attraction * dt;
    let expected_gap = gap - 0.5 * attraction / PLATE_MASS * dt * dt;
    let record = step_descriptor(&system, &before, &[voltage], dt).unwrap();
    assert!((record.x[0] - charge).abs() < 1e-18);
    assert!((record.x[1] - expected_gap).abs() < 2e-12);
    assert!((record.x[2] - expected_momentum).abs() < 5e-11);
    assert!(record.y[0].abs() < 1e-12);
    assert!(record.delta_h.abs() < 1e-14);
    assert!(record.balance_residual().abs() < 1e-14);
    assert_eq!(record.dissipated, 0.0);
    assert_eq!(record.supplied, 0.0);
    assert!(record.x[1] > 5e-5 && record.x[1] < 2e-4);
}

struct InvalidEffort;
impl Storage for InvalidEffort {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
    }
}

#[test]
fn descriptor_refuses_nonfinite_state_storage_effort_and_flow() {
    let system = bounded_plate();
    for (before, input) in [
        ([1e-8, 5e-5, 0.0, 0.0, 0.0], 100.0),
        ([f64::NAN, 1e-4, 0.0, 0.0, 0.0], 100.0),
        ([1e-8, 1e-4, 0.0, 0.0, 0.0], f64::NAN),
    ] {
        assert!(matches!(
            step_descriptor(&system, &before, &[input], 1e-3),
            Err(PhsError::DescriptorDomain { .. })
        ));
    }

    let invalid_effort = DescriptorPortHamiltonian::new(
        1,
        1,
        0,
        vec![0.0],
        vec![0.0],
        vec![],
        Box::new(InvalidEffort),
    )
    .unwrap();
    assert!(matches!(
        step_descriptor(&invalid_effort, &[0.0], &[], 1.0),
        Err(PhsError::DescriptorDomain { .. })
    ));
    let overflowing_flow = DescriptorPortHamiltonian::new(
        2,
        2,
        0,
        vec![0.0, f64::MAX, -f64::MAX, 0.0],
        vec![0.0; 4],
        vec![],
        Box::new(QuadraticStorage::new(vec![2.0, 0.0, 0.0, 2.0], 2).unwrap()),
    )
    .unwrap();
    assert!(matches!(
        step_descriptor(&overflowing_flow, &[1.0, 1.0], &[], 1.0),
        Err(PhsError::DescriptorDomain { .. })
    ));
}

struct BoundedQuartic;

impl BoundedQuartic {
    fn supported(x: &[f64]) -> bool {
        0.01 < x[0] && x[0] < 1.5
    }
}

impl Storage for BoundedQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if Self::supported(x) {
            QuarticCapacitor.hamiltonian(x)
        } else {
            f64::NAN
        }
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if Self::supported(x) {
            QuarticCapacitor.gradient(x, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        if Self::supported(a) && Self::supported(b) {
            QuarticCapacitor.midpoint_energy_remainder(a, b, mid)
        } else {
            Some(f64::NAN)
        }
    }
}

#[test]
fn descriptor_backtracks_only_to_regain_supported_nonlinear_storage() {
    let system = DescriptorPortHamiltonian::new(
        2,
        1,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0],
        vec![0.0, 1.0],
        Box::new(BoundedQuartic),
    )
    .unwrap();
    let q0 = 0.1_f64;
    let voltage = 2.0;
    let dt = 1.0;
    // The first undamped Newton update exceeds the explicit upper bound,
    // but the independent implicit scalar root lies inside support.
    let mut lower = q0;
    let mut upper = 1.5;
    for _ in 0..60 {
        let q1 = f64::midpoint(lower, upper);
        let effort = (q1.powi(3) + q1.powi(2) * q0 + q1 * q0.powi(2) + q0.powi(3)) / 4.0;
        if q1 - q0 - dt * (voltage - effort) > 0.0 {
            upper = q1;
        } else {
            lower = q1;
        }
    }
    let record = step_descriptor(&system, &[q0, voltage - q0.powi(3)], &[voltage], dt).unwrap();
    assert!((record.x[0] - f64::midpoint(lower, upper)).abs() < 2e-12);
    assert!(record.balance_residual().abs() < 5e-12);
}
