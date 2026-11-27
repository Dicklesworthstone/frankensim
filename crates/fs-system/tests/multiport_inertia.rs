//! fs-system's multi-port inertia, `elements::MultiportInertia`: conformance, tangents and exact
//! linearity over a sweep, its step and receipt against sums by hand, one row and one port against
//! `Mass` to the bit, its refusals, and two coupled ports on springs through the kernel against
//! the modal closed form.
//! Run with `cargo test -p fs-system --test multiport_inertia`.

// Exact float comparisons here assert bitwise agreement with `Mass`.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{
    Case, Observed, check_component, check_exact_linear, check_port_tangent,
};
use fs_ports::component::{Component, Proposal, Receipt};
use fs_system::SystemBuilder;
use fs_system::elements::{ElementError, Mass, Motion, MultiportInertia, Spring};

const KIND: PortKind = PortKind::MechanicalForceVelocity;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 12,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Deterministic numbers in `[lo, hi)` (splitmix64), so a sweep is the same on every run.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        lo + (hi - lo) * ((z >> 11) as f64 / (1_u64 << 53) as f64)
    }

    fn many(&mut self, count: usize, lo: f64, hi: f64) -> Vec<f64> {
        (0..count).map(|_| self.next(lo, hi)).collect()
    }
}

/// As in tests/resistance.rs: `A A^T` plus a quarter on the diagonal, each entry and its mirror
/// summed the same way, so the two agree bit for bit.
fn definite(numbers: &mut Numbers, n: usize) -> Vec<f64> {
    let a = numbers.many(n * n, -2.0, 2.0);
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let s: f64 = (0..n).map(|k| a[i * n + k] * a[j * n + k]).sum();
            m[i * n + j] = s + if i == j { 0.25 } else { 0.0 };
        }
    }
    m
}

/// One step proposed directly: the new state, the outputs and the receipt.
fn step<C: Component>(
    element: &C,
    old: &C::State,
    dt_s: f64,
    inputs: &[f64],
) -> (C::State, Vec<f64>, Receipt)
where
    C::State: Clone,
{
    let (mut new, mut outputs) = (old.clone(), vec![0.0; inputs.len()]);
    let proposal = with_cx(|cx| element.propose(old, dt_s, inputs, &mut new, &mut outputs, cx));
    let Ok(Proposal::Step(receipt)) = proposal else {
        panic!("no step at {inputs:?}");
    };
    (new, outputs, receipt)
}

#[test]
fn conforms_and_steps_by_the_midpoint_rule_with_a_lossless_receipt() {
    let mut numbers = Numbers(20_261_127);
    for (n, m) in [(1_usize, 1_usize), (1, 3), (2, 1), (2, 2), (3, 3), (5, 2)] {
        let mass = definite(&mut numbers, n);
        let map = numbers.many(n * m, -2.0, 2.0);
        let element = MultiportInertia::new(KIND, mass.clone(), map.clone(), m).expect("admitted");
        assert_eq!(element.ports().len(), m);
        for _ in 0..30 {
            let dt = numbers.next(1.0e-5, 2.0e-2);
            let old = numbers.many(n, -3.0, 3.0);
            let inputs = numbers.many(m, -40.0, 40.0);
            let stale = vec![9.0; n];
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s: dt,
                inputs: &inputs,
            };
            let (other_state, other_inputs) =
                (numbers.many(n, -3.0, 3.0), numbers.many(m, -4.0, 4.0));
            let other = Case {
                state: &other_state,
                stale: &stale,
                dt_s: dt,
                inputs: &other_inputs,
            };
            let (new, outputs, receipt) = step(&element, &old, dt, &inputs);
            with_cx(|cx| {
                let observed = check_component(&element, &case, cx)
                    .unwrap_or_else(|failure| panic!("conformance ({n}, {m}): {failure}"));
                assert!(matches!(observed, Observed::Stepped { .. }));
                // The kit's central difference steps by h = max(1, |inputs|) eps^(1/3); an output
                // with a large offset loses about 2 eps |output| / (2 h) of its slope there.
                let largest = |values: &[f64]| values.iter().fold(1.0_f64, |a, v| a.max(v.abs()));
                let h = largest(&inputs) * f64::EPSILON.cbrt();
                let absolute = 1.0e-12_f64.max(8.0 * f64::EPSILON * largest(&outputs) / h);
                check_port_tangent(&element, &case, 1.0e-6, absolute, cx)
                    .unwrap_or_else(|failure| panic!("tangent ({n}, {m}): {failure}"));
                check_exact_linear(&element, &case, &other, cx)
                    .unwrap_or_else(|failure| panic!("exact linear ({n}, {m}): {failure}"));
            });
            // By hand: M (v1 - v0) = dt G e row by row (the solve errs by a few n units of
            // |L| |D| |L^T| |v1 - v0|), the flows G^T (v0 + v1) / 2, and the kinetic energy's
            // change and the port work against the receipt.
            let (mut kinetic, mut size, mut work) = (0.0, 0.0, 0.0);
            for i in 0..n {
                let (mut momentum, mut force, mut row_size) = (0.0, 0.0, 0.0);
                for j in 0..n {
                    momentum += mass[i * n + j] * (new[j] - old[j]);
                    // The factor's |L| |D| |L^T| is bounded entrywise by sqrt(M_ii M_jj).
                    let root = (mass[i * n + i] * mass[j * n + j]).sqrt();
                    row_size += root * (new[j].abs() + old[j].abs());
                    kinetic += 0.5 * mass[i * n + j] * (new[i] * new[j] - old[i] * old[j]);
                    size += 0.5 * mass[i * n + j].abs() * (new[i] * new[j]).abs()
                        + 0.5 * mass[i * n + j].abs() * (old[i] * old[j]).abs();
                }
                for j in 0..m {
                    force += dt * map[i * m + j] * inputs[j];
                    row_size += (dt * map[i * m + j] * inputs[j]).abs();
                }
                let rows = 8.0 * (n + m) as f64 * f64::EPSILON;
                assert!(
                    (momentum - force).abs() <= rows * row_size,
                    "({n}, {m}) row {i}"
                );
            }
            for j in 0..m {
                let terms = (0..n).map(|i| map[i * m + j] * 0.5 * (old[i] + new[i]));
                let (flow, flow_size) = terms.fold((0.0, 0.0), |(f, s), t| (f + t, s + t.abs()));
                let bound = 4.0 * n as f64 * f64::EPSILON * flow_size;
                assert!((outputs[j] - flow).abs() <= bound, "({n}, {m}) flow {j}");
                work += dt * inputs[j] * outputs[j];
                size += (dt * inputs[j] * outputs[j]).abs();
            }
            let tolerance = receipt.energy_tolerance_j;
            assert_eq!((receipt.dissipated_j, receipt.numerical_loss_j), (0.0, 0.0));
            assert!(tolerance <= 64.0 * (n * n * m) as f64 * f64::EPSILON * size * 1.000_001);
            assert!(
                (receipt.stored_change_j - kinetic).abs() <= tolerance,
                "({n}, {m}) kinetic"
            );
            assert!(
                (receipt.stored_change_j - work).abs() <= tolerance,
                "({n}, {m}) work"
            );
        }
    }
}

#[test]
fn one_row_and_one_port_is_the_mass_step_to_the_bit() {
    let mut numbers = Numbers(1_234);
    for _ in 0..400 {
        let (m, v0) = (numbers.next(0.01, 100.0), numbers.next(-5.0, 5.0));
        let (dt, force) = (numbers.next(1.0e-6, 0.1), numbers.next(-50.0, 50.0));
        let inertia = MultiportInertia::new(KIND, vec![m], vec![1.0], 1).expect("one row");
        let mass = Mass::new(m).expect("mass");
        let start = Motion {
            position_m: 0.0,
            velocity_m_s: v0,
        };
        let (motion, by_mass, receipt_mass) = step(&mass, &start, dt, &[force]);
        let (velocity, by_inertia, receipt) = step(&inertia, &vec![v0], dt, &[force]);
        assert_eq!(velocity[0].to_bits(), motion.velocity_m_s.to_bits());
        assert_eq!(by_inertia[0].to_bits(), by_mass[0].to_bits());
        let bits = |r: Receipt| [r.stored_change_j.to_bits(), r.energy_tolerance_j.to_bits()];
        assert_eq!(bits(receipt), bits(receipt_mass), "receipt");
        assert_eq!(receipt, receipt_mass);
        let energy = inertia.stored_energy_j(&velocity).expect("energy");
        let by_mass = mass.stored_energy_j(&motion).expect("energy");
        assert_eq!(energy.to_bits(), by_mass.to_bits(), "stored energy");
    }
}

#[test]
fn refuses_a_matrix_that_is_not_symmetric_positive_definite_and_a_map_of_the_wrong_shape() {
    let refused = |mass: Vec<f64>, map: Vec<f64>, ports| {
        MultiportInertia::new(KIND, mass, map, ports).unwrap_err()
    };
    let count = |name, value| ElementError::Parameter { name, value };
    let square = "mass matrix entries (n x n, n >= 1)";
    let shape = "port map entries (n x ports, ports >= 1)";
    assert_eq!(
        refused(vec![1.0, 0.0, 1.0], vec![1.0], 1),
        count(square, 3.0)
    );
    assert_eq!(refused(Vec::new(), Vec::new(), 1), count(square, 0.0));
    assert_eq!(refused(vec![1.0], vec![1.0, 2.0], 1), count(shape, 2.0));
    assert_eq!(
        refused(vec![2.0, 0.0, 0.0, 2.0], vec![1.0; 3], 2),
        count(shape, 3.0)
    );
    assert_eq!(refused(vec![1.0], Vec::new(), 0), count(shape, 0.0));
    let nan = refused(vec![1.0, f64::NAN, f64::NAN, 1.0], vec![1.0; 2], 1);
    assert!(
        matches!(nan, ElementError::Parameter { name: "mass matrix entry", value } if value.is_nan())
    );
    assert_eq!(
        refused(vec![1.0], vec![f64::INFINITY], 1),
        count("port map entry", f64::INFINITY)
    );
    assert_eq!(
        refused(
            vec![1.0, 0.5, 0.5 + f64::EPSILON / 2.0, 1.0],
            vec![1.0; 2],
            1
        ),
        ElementError::Asymmetric { row: 1, column: 0 }
    );
    // A signed zero equals its mirror.
    MultiportInertia::new(KIND, vec![1.0, 0.0, -0.0, 1.0], vec![1.0; 2], 1).expect("symmetric");
    // Pivots by hand, as for Resistance: 1 - 2 * 2 = -3; 1 - 1 * 1 = 0; -1.
    for (mass, pivot, value) in [
        (vec![1.0, 2.0, 2.0, 1.0], 1, -3.0),
        (vec![1.0, 1.0, 1.0, 1.0], 1, 0.0),
        (vec![-1.0], 0, -1.0),
    ] {
        let rows = mass.len().isqrt();
        assert_eq!(
            refused(mass, vec![1.0; rows], 1),
            ElementError::NotPositiveDefinite { pivot, value }
        );
    }
    // A state with a velocity for each row only, finite: the step and the codec refuse others.
    let element =
        MultiportInertia::new(KIND, vec![2.0, 1.0, 1.0, 2.0], vec![1.0; 2], 1).expect("ok");
    let (mut new, mut out) = (Vec::new(), vec![0.0]);
    let short = with_cx(|cx| element.propose(&vec![1.0], 1.0e-3, &[1.0], &mut new, &mut out, cx));
    assert!(matches!(short, Err(ElementError::Parameter { .. })));
    let mut words = Vec::new();
    assert!(!element.encode_state(&vec![1.0, f64::NAN], &mut words));
    assert!(element.encode_state(&vec![1.0, -2.0], &mut words));
    assert_eq!(element.decode_state(&words), Some(vec![1.0, -2.0]));
    assert_eq!(element.decode_state(&words[..1]), None);
}

/// Two coordinates with the coupled mass matrix `M = [[2, 1], [1, 2]]` kg, each joined to a spring
/// of `k = 3` N/m anchored at the origin, started at rest length with velocities `[1, 0]` m/s. The
/// modes are `[1, 1]` at `omega^2 = k / 3` and `[1, -1]` at `omega^2 = k`, half the start in
/// each. The kernel's step is the implicit midpoint rule, which turns each mode's `(x, v / omega)`
/// by `2 atan(omega dt / 2)` per step, so after `n` steps the velocities are
/// `(cos n theta_1 +- cos n theta_2) / 2`: the second coordinate beats.
#[test]
fn two_coupled_ports_on_springs_beat_as_the_midpoint_rule_turns_each_mode() {
    let (k, dt, steps) = (3.0_f64, 1.0e-2_f64, 2_000_u32);
    let inertia =
        MultiportInertia::new(KIND, vec![2.0, 1.0, 1.0, 2.0], vec![1.0, 0.0, 0.0, 1.0], 2)
            .expect("inertia");
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", inertia, vec![1.0, 0.0])
        .expect("add body");
    let mut springs = Vec::new();
    for (name, port) in [("spring0", "p0"), ("spring1", "p1")] {
        let spring = builder
            .add_with_tangent(name, Spring::new(k, 0.0).expect("spring"), 0.0)
            .expect("add spring");
        let (p, a) = (
            builder.port(body, port).expect("p"),
            builder.port(spring, "a").expect("a"),
        );
        builder.common_flow(p, &[a]).expect("junction");
        let b = builder.port(spring, "b").expect("b");
        builder
            .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
            .expect("anchor");
        springs.push(spring);
    }
    let mut system = builder.seal().expect("seal");
    let turn = |omega: f64| 2.0 * (0.5 * omega * dt).atan();
    let (slow, fast) = (turn((k / 3.0).sqrt()), turn(k.sqrt()));
    let mut worst = 0.0_f64;
    with_cx(|cx| {
        for n in 1..=steps {
            system.step(dt, cx).expect("step");
            let v = system.state::<Vec<f64>>(body).expect("velocities");
            let (a, b) = ((f64::from(n) * slow).cos(), (f64::from(n) * fast).cos());
            worst = worst
                .max((v[0] - 0.5 * (a + b)).abs())
                .max((v[1] - 0.5 * (a - b)).abs());
        }
    });
    // Each step rounds the velocities, of size at most 1, and the closed form's angles by a few
    // ulps; the start's energy is 1 J.
    let rounding = 8.0 * f64::from(steps) * f64::EPSILON;
    assert!(worst <= rounding, "worst {worst:e} against {rounding:e}");
    let v = system.state::<Vec<f64>>(body).expect("velocities").clone();
    let x: Vec<f64> = springs
        .iter()
        .map(|&s| -*system.state::<f64>(s).expect("extension"))
        .collect();
    let energy = 0.5 * (2.0 * v[0] * v[0] + 2.0 * v[0] * v[1] + 2.0 * v[1] * v[1])
        + 0.5 * k * (x[0] * x[0] + x[1] * x[1]);
    assert!((energy - 1.0).abs() <= rounding, "energy {energy}");
    let ledger = system.ledger();
    assert!(ledger.stored_change_j().abs() <= rounding);
    assert!(ledger.energy_residual_j().abs() <= rounding);
}
