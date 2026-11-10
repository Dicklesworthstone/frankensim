//! fs-system's resistor N-port, `elements::Resistance`: conformance and tangents over a sweep, its
//! efforts and heat against sums by hand, the diagonal case against `Damper`, its refusals, and one
//! mass whose motion it damps through the kernel.
//! Run with `cargo test -p fs-system --test resistance`.

// Exact float comparisons here assert bitwise agreement with a sum by hand or with `Damper`.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent};
use fs_system::SystemBuilder;
use fs_system::elements::{Damper, ElementError, HeatCapacity, Mass, Motion, Resistance};

const KIND: PortKind = PortKind::MechanicalForceVelocity;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 11,
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
}

/// The mandatory and tangent checks at one operating point; the outputs of the step.
fn conforms<C>(
    element: &C,
    state: &C::State,
    stale: &C::State,
    dt_s: f64,
    inputs: &[f64],
) -> Vec<f64>
where
    C: PortTangent,
    C::State: PartialEq + core::fmt::Debug,
{
    let case = Case {
        state,
        stale,
        dt_s,
        inputs,
    };
    with_cx(|cx| {
        let observed = check_component(element, &case, cx)
            .unwrap_or_else(|failure| panic!("conformance at {inputs:?}: {failure}"));
        let outputs = match observed {
            Observed::Stepped { outputs, .. } => outputs,
            Observed::Shortened { dt_s } => panic!("shortened to {dt_s} at {inputs:?}"),
        };
        // As in tests/elements.rs: the kit's central difference steps by
        // h = max(1, |inputs|) eps^(1/3), and an output with a large offset loses about
        // 2 eps |output| / (2 h) of its slope to roundoff there.
        let largest = |values: &[f64]| values.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        let h = largest(inputs) * f64::EPSILON.cbrt();
        let absolute = 1.0e-12_f64.max(8.0 * f64::EPSILON * largest(&outputs) / h);
        check_port_tangent(element, &case, 1.0e-6, absolute, cx)
            .unwrap_or_else(|failure| panic!("tangent at {inputs:?}: {failure}"));
        outputs
    })
}

/// A symmetric positive definite matrix: `A A^T` plus a quarter on the diagonal, each entry and
/// its mirror summed the same way, so the two agree bit for bit.
fn definite(numbers: &mut Numbers, n: usize) -> Vec<f64> {
    let a: Vec<f64> = (0..n * n).map(|_| numbers.next(-2.0, 2.0)).collect();
    let mut r = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let s: f64 = (0..n).map(|k| a[i * n + k] * a[j * n + k]).sum();
            r[i * n + j] = s + if i == j { 0.25 } else { 0.0 };
        }
    }
    r
}

#[test]
fn conforms_and_returns_r_u_with_all_its_power_as_heat() {
    let mut numbers = Numbers(20_261_110);
    for n in [1_usize, 2, 3, 5] {
        let matrix = definite(&mut numbers, n);
        let element = Resistance::new(KIND, matrix.clone()).expect("definite");
        assert_eq!((element.rows(), element.ports().len()), (n, 2));
        for _ in 0..40 {
            let dt = numbers.next(1.0e-5, 2.0e-2);
            let mut inputs: Vec<f64> = (0..n).map(|_| numbers.next(-3.0, 3.0)).collect();
            inputs.push(numbers.next(250.0, 900.0));
            let out = conforms(&element, &(), &(), dt, &inputs);
            let (mut power, mut scale, mut diagonal) = (0.0, 0.0, 0.0);
            for i in 0..n {
                let mut effort = 0.0;
                for j in 0..n {
                    effort += matrix[i * n + j] * inputs[j];
                    scale += (matrix[i * n + j] * inputs[i] * inputs[j]).abs();
                }
                assert_eq!(out[i], effort, "n {n} row {i}");
                power += effort * inputs[i];
                diagonal += matrix[i * n + i].sqrt() * inputs[i].abs();
            }
            let bound = 64.0 * (n as f64 + 2.0) * f64::EPSILON * (scale + diagonal * diagonal);
            assert!(
                (out[n] + power).abs() <= bound,
                "n {n}: heat {} against u.Ru {power}",
                out[n]
            );
        }
    }
}

#[test]
fn a_diagonal_matrix_returns_each_row_the_damper_effort_exactly() {
    let element = Resistance::new(KIND, vec![2.5, 0.0, 0.0, 0.7]).expect("diagonal");
    let (a, b) = (Damper::new(2.5).expect("a"), Damper::new(0.7).expect("b"));
    let mut numbers = Numbers(7);
    for _ in 0..50 {
        let (u, v) = (numbers.next(-3.0, 3.0), numbers.next(-3.0, 3.0));
        let out = conforms(&element, &(), &(), 1.0e-3, &[u, v, 300.0]);
        let along_a = conforms(&a, &(), &(), 1.0e-3, &[u, 0.0, 300.0]);
        let along_b = conforms(&b, &(), &(), 1.0e-3, &[v, 0.0, 300.0]);
        assert_eq!(out[0], along_a[0]);
        assert_eq!(out[1], along_b[0]);
    }
}

#[test]
fn refuses_what_is_not_a_symmetric_positive_definite_square() {
    let refused = |matrix: Vec<f64>| Resistance::new(KIND, matrix).unwrap_err();
    assert!(matches!(
        refused(vec![1.0, 0.0, 1.0]),
        ElementError::Parameter { .. }
    ));
    assert!(matches!(
        refused(Vec::new()),
        ElementError::Parameter { .. }
    ));
    assert!(matches!(
        refused(vec![1.0, f64::NAN, f64::NAN, 1.0]),
        ElementError::Parameter { .. }
    ));
    assert_eq!(
        refused(vec![1.0, 0.5, 0.5 + f64::EPSILON / 2.0, 1.0]),
        ElementError::Asymmetric { row: 1, column: 0 }
    );
    // A signed zero equals its mirror.
    Resistance::new(KIND, vec![1.0, 0.0, -0.0, 1.0]).expect("symmetric");
    // Pivots by hand: 1 - 2 * 2 = -3; 1 - 1 * 1 = 0; -1.
    assert_eq!(
        refused(vec![1.0, 2.0, 2.0, 1.0]),
        ElementError::NotPositiveDefinite {
            pivot: 1,
            value: -3.0
        }
    );
    assert_eq!(
        refused(vec![1.0, 1.0, 1.0, 1.0]),
        ElementError::NotPositiveDefinite {
            pivot: 1,
            value: 0.0
        }
    );
    assert_eq!(
        refused(vec![-1.0]),
        ElementError::NotPositiveDefinite {
            pivot: 0,
            value: -1.0
        }
    );
}

/// A mass on a common-flow junction with a one-row resistance: the kernel's step is the
/// trapezoidal rule, `v' = v (1 - r dt / 2m) / (1 + r dt / 2m)`, worked here apart from the
/// kernel; the heat the resistance books is the kinetic energy the mass lost, and the body that
/// takes the heat gains it. Tolerances as in tests/elements.rs's damped oscillator: the interface
/// is solved to 1e-9 of the scale each step, and the body's temperature near 300 K is rounded to
/// an ulp of 300 K at each step.
#[test]
fn a_mass_loses_to_it_exactly_what_the_trapezoidal_rule_says() {
    let (m, r, v0) = (0.25_f64, 1.5_f64, -0.3_f64);
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(m).expect("mass"),
            Motion {
                position_m: 0.0,
                velocity_m_s: v0,
            },
        )
        .expect("add mass");
    let loss = builder
        .add_with_tangent("loss", Resistance::new(KIND, vec![r]).expect("r"), ())
        .expect("add resistance");
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(2.0).expect("body"), 300.0)
        .expect("add body");
    let p = builder.port(mass, "p").expect("p");
    let u = builder.port(loss, "u").expect("u");
    builder.common_flow(p, &[u]).expect("junction");
    let (body_heat, loss_heat) = (
        builder.port(body, "heat").expect("body heat"),
        builder.port(loss, "heat").expect("loss heat"),
    );
    builder
        .common_effort(body_heat, &[loss_heat])
        .expect("heat route");
    let mut system = builder.seal().expect("seal");
    assert!(system.guess(loss_heat, &[300.0]));
    let (dt, steps) = (2.0e-3_f64, 500_u32);
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).expect("step");
        }
    });
    let factor = (1.0 - r * dt / (2.0 * m)) / (1.0 + r * dt / (2.0 * m));
    let expected = v0 * factor.powi(steps as i32);
    let v = system.state::<Motion>(mass).expect("motion").velocity_m_s;
    assert!(
        (v - expected).abs() <= 1.0e-9 * v0.abs(),
        "{v} against {expected}"
    );
    let lost = 0.5 * m * (v0 * v0 - v * v);
    let ledger = system.ledger();
    let dissipated = ledger.nodes[loss.index()].dissipated_j.value();
    assert!(
        (dissipated - lost).abs() <= 1.0e-9 * lost,
        "{dissipated} against {lost}"
    );
    let account = ledger.nodes[body.index()].stored_change_j.value();
    assert!(
        (account - dissipated).abs() <= 1.0e-12 * dissipated,
        "{account} against {dissipated}"
    );
    let warmed = 2.0 * (*system.state::<f64>(body).expect("temperature") - 300.0);
    let rounding = 2.0 * f64::from(steps) * f64::EPSILON * 300.0;
    assert!(
        (warmed - dissipated).abs() <= rounding,
        "{warmed} against {dissipated}"
    );
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * 0.5 * m * v0 * v0);
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}
