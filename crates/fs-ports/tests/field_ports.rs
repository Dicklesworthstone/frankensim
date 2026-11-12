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
