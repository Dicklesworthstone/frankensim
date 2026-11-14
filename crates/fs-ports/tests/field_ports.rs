//! G0/G3: field work uses the declared duality, including off-diagonal terms.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Class, Component, ContractViolation, DualityMetric, PortLayout, PortRole, PortSpec, Proposal,
    Receipt, audit_step,
};
use fs_ports::{PortKind, StableId};

fn width(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).expect("positive width")
}

fn metric(entries: [f64; 4]) -> DualityMetric {
    DualityMetric::new(width(2), entries.to_vec()).expect("positive-definite metric")
}

fn port(role: PortRole, entries: [f64; 4]) -> PortSpec {
    PortSpec::field_power(
        StableId::new("surface").expect("id"),
        PortKind::MechanicalForceVelocity,
        role,
        metric(entries),
    )
    .expect("field port")
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} within {tolerance}"
    );
}

#[test]
fn off_diagonal_pairing_is_independent_of_which_side_is_the_input() {
    let entries = [5.0, -2.0, -2.0, 4.0];
    let effort = [3.0, -4.0];
    let flow = [2.0, -3.0];
    // M f = [16, -16], so e^T M f = 112 W, not e.f = 18 W.
    for role in [PortRole::EffortIn, PortRole::FlowIn] {
        let field = port(role, entries);
        let (input, output) = if role == PortRole::EffortIn {
            (&effort, &flow)
        } else {
            (&flow, &effort)
        };
        assert_eq!(field.energy_in_j(0.125, input, output), 14.0);
        assert_eq!(field.width(), 2);
        assert_eq!(field.field_metric().expect("metric").entries(), &entries);
        assert!(field.energy_in_j(0.125, &input[..1], output).is_nan());
    }
}

#[test]
fn identity_metric_keeps_the_ordinary_pairing() {
    for role in [PortRole::EffortIn, PortRole::FlowIn] {
        let ordinary = PortSpec::power(
            StableId::new("ordinary").expect("id"),
            PortKind::MechanicalForceVelocity,
            role,
            width(2),
        )
        .expect("power port");
        let field = port(role, [1.0, 0.0, 0.0, 1.0]);
        assert!(ordinary.field_metric().is_none());
        for scale in [1.0e-100, 1.0, 1.0e100] {
            let input = [3.25 * scale, -0.125 * scale];
            let output = [-1.75 / scale, 5.5 / scale];
            assert_eq!(
                field.energy_in_j(0.025, &input, &output).to_bits(),
                ordinary.energy_in_j(0.025, &input, &output).to_bits()
            );
        }
    }
}

#[test]
fn metric_solve_recovers_an_independently_known_field() {
    let field_metric = metric([4.0, 2.0, 2.0, 3.0]);
    let mut product = [f64::NAN; 2];
    field_metric
        .multiply(&[1.25, -2.5], &mut product)
        .expect("multiply");
    assert_eq!(product, [0.0, -5.0]);
    let mut recovered = [f64::NAN; 2];
    field_metric
        .solve(&[0.0, -5.0], &mut recovered)
        .expect("solve");
    close(recovered[0], 1.25, 1.0e-15);
    close(recovered[1], -2.5, 2.0e-15);
    assert_eq!(metric([1.0, -0.0, 0.0, 2.0]), metric([1.0, 0.0, -0.0, 2.0]));
}

#[test]
fn invalid_metrics_and_unrepresentable_factorizations_are_refused() {
    assert!(matches!(
        DualityMetric::new(width(2), vec![1.0, 0.0, 1.0]),
        Err(ContractViolation::MetricShape { .. })
    ));
    assert!(matches!(
        DualityMetric::new(width(usize::MAX), vec![]),
        Err(ContractViolation::MetricShape { .. })
    ));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            DualityMetric::new(width(2), vec![1.0, bad, bad, 1.0]),
            Err(ContractViolation::MetricCoefficient { .. })
        ));
    }
    assert!(matches!(
        DualityMetric::new(width(2), vec![2.0, 0.5, 0.25, 1.0]),
        Err(ContractViolation::MetricAsymmetric { .. })
    ));
    for entries in [[1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 2.0, 1.0]] {
        assert!(matches!(
            DualityMetric::new(width(2), entries.to_vec()),
            Err(ContractViolation::MetricPivot { index: 1, .. })
        ));
    }
    assert!(matches!(
        DualityMetric::new(width(1), vec![-1.0]),
        Err(ContractViolation::MetricPivot { index: 0, .. })
    ));
    assert!(matches!(
        DualityMetric::new(width(2), vec![1.0e-300, 1.0e300, 1.0e300, 1.0]),
        Err(ContractViolation::MetricFactorization { .. })
    ));
    assert!(matches!(
        PortSpec::field_power(
            StableId::new("bad_role").expect("id"),
            PortKind::MechanicalForceVelocity,
            PortRole::SignalIn,
            metric([1.0, 0.0, 0.0, 1.0]),
        ),
        Err(ContractViolation::PortRole { .. })
    ));
}

#[test]
fn metric_actions_refuse_wrong_shapes_and_nonfinite_arithmetic() {
    let field_metric = metric([2.0, 0.0, 0.0, 1.0]);
    let mut output = [31.0, 37.0];
    assert!(matches!(
        field_metric.multiply(&[1.0], &mut output),
        Err(ContractViolation::VectorLength { .. })
    ));
    assert_eq!(output, [31.0, 37.0]);
    assert!(matches!(
        field_metric.solve(&[1.0, 2.0], &mut output[..1]),
        Err(ContractViolation::VectorLength { .. })
    ));
    assert_eq!(output, [31.0, 37.0]);
    assert!(matches!(
        field_metric.solve(&[f64::NAN, 1.0], &mut output),
        Err(ContractViolation::NonFinite { .. })
    ));
    assert!(matches!(
        field_metric.multiply(&[f64::MAX, 1.0], &mut output),
        Err(ContractViolation::NonFinite { .. })
    ));
    let tiny = DualityMetric::new(width(1), vec![1.0e-300]).expect("positive scalar");
    assert!(matches!(
        tiny.solve(&[1.0e300], &mut [0.0]),
        Err(ContractViolation::NonFinite { .. })
    ));
}

/// A surface with areal mass sigma and the trace mass matrix [[2,1],[1,2]].
/// Its acceleration coefficients are traction / sigma, and the kinetic
/// energy is sigma/2 * v^T M v. The law does not call the metric helpers.
struct InertialSurface {
    sigma: f64,
    ports: [PortSpec; 1],
}

impl Component for InertialSurface {
    type State = [f64; 2];
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(self.sigma * (state[0] * state[0] + state[0] * state[1] + state[1] * state[1]))
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint()?;
        let dv = [dt_s * inputs[0] / self.sigma, dt_s * inputs[1] / self.sigma];
        *new = [old[0] + dv[0], old[1] + dv[1]];
        let mean = [0.5 * (old[0] + new[0]), 0.5 * (old[1] + new[1])];
        outputs.copy_from_slice(&mean);
        let stored_change =
            self.sigma * (mean[0] * (2.0 * dv[0] + dv[1]) + mean[1] * (dv[0] + 2.0 * dv[1]));
        let tolerance = 64.0
            * f64::EPSILON
            * (self.stored_energy_j(old)? + self.stored_energy_j(new)? + stored_change.abs());
        Ok(Proposal::Step(Receipt::lossless(stored_change, tolerance)))
    }
}

#[test]
fn an_inertial_field_obeys_traction_work_and_rejects_a_euclidean_account() {
    let component = InertialSurface {
        sigma: 3.0,
        ports: [port(PortRole::EffortIn, [2.0, 1.0, 1.0, 2.0])],
    };
    let old = [0.2, -0.4];
    let traction = [4.0, -1.0];
    let dt_s = 0.125;
    let case = Case {
        state: &old,
        stale: &[7.0, -8.0],
        dt_s,
        inputs: &traction,
    };
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        let Observed::Stepped {
            receipt,
            audit,
            outputs,
        } = check_component(&component, &case, &cx).expect("field component conforms")
        else {
            panic!("unexpected shortening");
        };
        close(outputs[0], old[0] + 0.5 * dt_s * traction[0] / 3.0, 1.0e-16);
        close(outputs[1], old[1] + 0.5 * dt_s * traction[1] / 3.0, 1.0e-16);
        close(audit.balance_defect_j, 0.0, receipt.energy_tolerance_j);

        let mut new = [0.0; 2];
        let mut output = [0.0; 2];
        component
            .propose(&old, dt_s, &traction, &mut new, &mut output, &cx)
            .expect("proposal");
        let wrong_receipt = Receipt::lossless(
            dt_s * (traction[0] * output[0] + traction[1] * output[1]),
            receipt.energy_tolerance_j,
        );
        let layout = PortLayout::new(component.ports()).expect("layout");
        assert!(matches!(
            audit_step(
                component.ports(),
                &layout,
                dt_s,
                &traction,
                &output,
                &wrong_receipt,
                component.stored_energy_j(&old).expect("old energy"),
                component.stored_energy_j(&new).expect("new energy"),
            ),
            Err(ContractViolation::EnergyBalance { .. })
        ));
    });
}

/// A small generator of the same numbers on every build (splitmix64), as
/// fs-system's load-line tests draw theirs.
struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn between(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * (self.next() >> 11) as f64 / 9_007_199_254_740_992.0
    }

    /// A number in `[-1, 1)` times a power of ten drawn in `[-3, 3)`.
    fn spread(&mut self) -> f64 {
        self.between(-1.0, 1.0) * 10.0_f64.powf(self.between(-3.0, 3.0))
    }
}

#[test]
fn a_field_port_pairs_drawn_vectors_by_its_metric() {
    let mut draw = Draw(20_261_114);
    let mut worst = 0.0_f64;
    for case in 0..400 {
        let n = 1 + (draw.next() % 6) as usize;
        // A^T A + n I, scaled: positive definite, and symmetric to the bit,
        // since a mirrored entry sums the same products in the same order.
        let a: Vec<f64> = (0..n * n).map(|_| draw.between(-1.0, 1.0)).collect();
        let scale = 10.0_f64.powf(draw.between(-6.0, 6.0));
        let mut entries = vec![0.0; n * n];
        for row in 0..n {
            for column in 0..n {
                let mut sum = if row == column { n as f64 } else { 0.0 };
                for k in 0..n {
                    sum += a[k * n + row] * a[k * n + column];
                }
                entries[row * n + column] = scale * sum;
            }
        }
        let field_metric = DualityMetric::new(width(n), entries.clone()).expect("metric");
        let effort: Vec<f64> = (0..n).map(|_| draw.spread()).collect();
        let flow: Vec<f64> = (0..n).map(|_| draw.spread()).collect();
        let dt = draw.between(1.0e-6, 0.1);
        // Independent: every term e_i M_ij f_j in one running sum, columns
        // outermost, the other way round from the port's row sums.
        let (mut want, mut size) = (0.0_f64, 0.0_f64);
        for column in 0..n {
            for row in 0..n {
                let term = effort[row] * entries[row * n + column] * flow[column];
                want += term;
                size += term.abs();
            }
        }
        want *= dt;
        size *= dt;
        // With u = eps / 2 and S the sum of |e_i M_ij f_j|: the port's row
        // sums, its products with e and its sum of them are off by at most
        // (2n + 1) u S; the judge's n^2 terms of two products each by
        // (n^2 + 2) u S; dt rounds each once more. That is (n^2 + 2n + 5) u S,
        // and (n^2 + 2n + 4) eps S holds it about twice over.
        let allowed = (n * n + 2 * n + 4) as f64 * f64::EPSILON * size;
        let mut paired = Vec::new();
        for role in [PortRole::EffortIn, PortRole::FlowIn] {
            let field = PortSpec::field_power(
                StableId::new("surface").expect("id"),
                PortKind::MechanicalForceVelocity,
                role,
                field_metric.clone(),
            )
            .expect("field port");
            let (input, output) = if role == PortRole::EffortIn {
                (&effort, &flow)
            } else {
                (&flow, &effort)
            };
            let got = field.energy_in_j(dt, input, output);
            assert!(
                (got - want).abs() <= allowed,
                "case {case}, width {n}: {got:e} J against {want:e} within {allowed:e}"
            );
            worst = worst.max((got - want).abs() / allowed);
            paired.push(got.to_bits());
        }
        // Whichever side the port receives, the same bits.
        assert_eq!(paired[0], paired[1], "case {case}");
    }
    eprintln!("drawn field pairings: at most {worst:.4} of the allowance from e^T M f");
}
