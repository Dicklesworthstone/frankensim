//! G0: the Component contract v0 on a reference implementation.
//!
//! The reference is a linear spring-mass-damper with a mechanical port and a
//! heat port, stepped by the implicit midpoint rule. For a linear system that
//! rule balances energy exactly: the work of the held force on the mean
//! velocity equals the stored change plus the damper loss, and the loss
//! leaves through the heat port. Three deliberately broken variants each fail
//! the law they break.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{
    Case, Failure, Observed, check_component, check_port_tangent,
};
use fs_ports::component::{
    Carrier, Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, PortTangent,
    Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Oscillator {
    position_m: f64,
    velocity_m_s: f64,
}

/// Which law, if any, the component breaks on purpose.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fault {
    None,
    /// Reports only 90 % of the damper loss as heat.
    HeatLeak,
    /// Adds the stale candidate's position to the result.
    ReadsCandidate,
    /// Returns the VJP with the wrong sign.
    VjpSign,
}

#[derive(Debug)]
struct SpringMassDamper {
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    ports: [PortSpec; 2],
    fault: Fault,
}

impl SpringMassDamper {
    fn new(fault: Fault) -> Self {
        let drive = PortSpec::power(
            StableId::new("drive").expect("id"),
            PortKind::MechanicalForceVelocity,
            PortRole::EffortIn,
            NonZeroUsize::MIN,
        )
        .expect("power port");
        let heat =
            PortSpec::heat(StableId::new("heat").expect("id"), PortRole::EffortIn).expect("heat");
        Self {
            mass_kg: 0.25,
            stiffness_n_m: 400.0,
            damping_n_s_m: 1.5,
            ports: [drive, heat],
            fault,
        }
    }

    /// Mean velocity over the step and its derivative with respect to force.
    fn mean_velocity(&self, old: &Oscillator, dt_s: f64, force_n: f64) -> (f64, f64) {
        let half = 0.5 * dt_s / self.mass_kg;
        let denominator = 1.0 + half * self.damping_n_s_m + 0.5 * dt_s * half * self.stiffness_n_m;
        let mean = (old.velocity_m_s + half * (force_n - self.stiffness_n_m * old.position_m))
            / denominator;
        (mean, half / denominator)
    }
}

impl Component for SpringMassDamper {
    type State = Oscillator;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Oscillator) -> Result<f64, Cancelled> {
        Ok(0.5 * self.mass_kg * state.velocity_m_s * state.velocity_m_s
            + 0.5 * self.stiffness_n_m * state.position_m * state.position_m)
    }

    fn propose(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Oscillator,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let (force_n, temperature_k) = (inputs[0], inputs[1]);
        let stale_position = new.position_m;
        let (mean, _) = self.mean_velocity(old, dt_s, force_n);
        let mid_position = old.position_m + 0.5 * dt_s * mean;
        new.position_m = old.position_m + dt_s * mean;
        new.velocity_m_s = 2.0 * mean - old.velocity_m_s;
        if self.fault == Fault::ReadsCandidate {
            new.position_m += 1.0e-3 * stale_position;
        }
        let loss_w = self.damping_n_s_m * mean * mean;
        let heat_out_w = if self.fault == Fault::HeatLeak {
            0.9 * loss_w
        } else {
            loss_w
        };
        outputs[0] = mean;
        outputs[1] = -heat_out_w;
        // Stored change without cancellation: m v_mean dv + k x_mid dx.
        let stored_change_j = self.mass_kg * mean * (new.velocity_m_s - old.velocity_m_s)
            + self.stiffness_n_m * mid_position * (dt_s * mean);
        let scale_j = (dt_s * force_n * mean).abs()
            + dt_s * loss_w
            + self.stored_energy_j(old)?
            + self.stored_energy_j(new)?;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: dt_s * loss_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * loss_w / temperature_k,
            energy_tolerance_j: 16.0 * f64::EPSILON * scale_j,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for SpringMassDamper {
    fn jvp(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        let (mean, slope) = self.mean_velocity(old, dt_s, inputs[0]);
        d_outputs[0] = slope * d_inputs[0];
        d_outputs[1] = -2.0 * self.damping_n_s_m * mean * slope * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        let (mean, slope) = self.mean_velocity(old, dt_s, inputs[0]);
        let sign = if self.fault == Fault::VjpSign {
            -1.0
        } else {
            1.0
        };
        bar_inputs[0] = sign
            * (slope * bar_outputs[0] - 2.0 * self.damping_n_s_m * mean * slope * bar_outputs[1]);
        bar_inputs[1] = 0.0;
        Ok(())
    }
}

const START: Oscillator = Oscillator {
    position_m: 2.0e-3,
    velocity_m_s: -0.3,
};
const STALE: Oscillator = Oscillator {
    position_m: 7.0,
    velocity_m_s: 11.0,
};
const INPUTS: [f64; 2] = [1.25, 293.15];

fn case() -> Case<'static, Oscillator> {
    Case {
        state: &START,
        stale: &STALE,
        dt_s: 1.0e-3,
        inputs: &INPUTS,
    }
}

#[test]
fn reference_component_passes_conformance() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::None);
    let observed = with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect("conforms");
    let Observed::Stepped {
        receipt,
        audit,
        outputs,
    } = observed
    else {
        panic!("the reference takes the step");
    };
    assert!(receipt.dissipated_j > 0.0);
    // The damper loss leaves as heat: heat-port energy in is minus the loss.
    let heat_in_j = component.ports()[1].energy_in_j(1.0e-3, &INPUTS[1..], &outputs[1..]);
    assert!((heat_in_j + receipt.dissipated_j).abs() <= receipt.energy_tolerance_j);
    assert!(audit.balance_defect_j.abs() <= receipt.energy_tolerance_j);
    // Entropy leaves through the heat port at the port temperature.
    let entropy_out = -audit.port_entropy_in_j_per_k;
    assert!((entropy_out - receipt.entropy_generated_j_per_k).abs() <= 1.0e-18);
    with_cx(&gate, |cx| {
        check_port_tangent(&component, &case(), 1.0e-7, 1.0e-12, cx)
    })
    .expect("tangent conforms");
}

#[test]
fn energy_closes_over_many_steps_without_drift() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::None);
    let layout = PortLayout::new(component.ports()).expect("layout");
    let dt_s = 2.0e-4;
    let (mut state, mut next) = (START, START);
    let mut outputs = vec![0.0; layout.output_len()];
    let start_j = component.stored_energy_j(&state).expect("energy");
    let (mut work_j, mut heat_j) = (0.0, 0.0);
    with_cx(&gate, |cx| {
        for step in 0..20_000_u32 {
            let force_n = 0.8 * (f64::from(step) * dt_s * 180.0).sin();
            let inputs = [force_n, 300.0];
            let Proposal::Step(_) = component
                .propose(&state, dt_s, &inputs, &mut next, &mut outputs, cx)
                .expect("step")
            else {
                panic!("the reference never shortens");
            };
            work_j += component.ports()[0].energy_in_j(dt_s, &inputs[..1], &outputs[..1]);
            heat_j -= component.ports()[1].energy_in_j(dt_s, &inputs[1..], &outputs[1..]);
            core::mem::swap(&mut state, &mut next);
        }
    });
    let end_j = component.stored_energy_j(&state).expect("energy");
    let defect_j = work_j - heat_j - (end_j - start_j);
    let scale_j = work_j.abs() + heat_j.abs() + start_j + end_j;
    assert!(heat_j > 0.0);
    assert!(
        defect_j.abs() <= 1.0e-11 * scale_j,
        "defect {defect_j} J of {scale_j} J"
    );
}

#[test]
fn a_heat_leak_fails_the_first_law() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::HeatLeak);
    let failure = with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("leak");
    assert!(
        matches!(
            failure,
            Failure::Contract(ContractViolation::EnergyBalance { .. })
        ),
        "{failure}"
    );
}

#[test]
fn reading_the_candidate_slot_fails_purity() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::ReadsCandidate);
    let failure =
        with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("impure");
    let Failure::Law { law, .. } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    assert_eq!(law, "propose reads nothing from the candidate slot");
}

#[test]
fn a_wrong_vjp_sign_fails_the_transpose_test() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::VjpSign);
    let failure = with_cx(&gate, |cx| {
        check_port_tangent(&component, &case(), 1.0e-7, 1.0e-12, cx)
    })
    .expect_err("wrong sign");
    let Failure::Law { law, .. } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    assert_eq!(law, "vjp is the transpose of jvp");
}

#[test]
fn cancellation_reaches_the_caller_as_the_owner_error() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    let component = SpringMassDamper::new(Fault::None);
    let failure =
        with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("cancelled");
    assert!(matches!(failure, Failure::Owner(Cancelled)), "{failure}");
}

#[test]
fn layout_orders_ports_and_refuses_duplicates() {
    let component = SpringMassDamper::new(Fault::None);
    let layout = PortLayout::new(component.ports()).expect("layout");
    assert_eq!(
        (layout.ports(), layout.input_len(), layout.output_len()),
        (2, 2, 2)
    );
    assert_eq!(layout.input_range(1), 1..2);
    assert_eq!(component.ports()[1].carrier(), Carrier::Heat);
    let twice = [component.ports()[0].clone(), component.ports()[0].clone()];
    assert!(matches!(
        PortLayout::new(&twice),
        Err(ContractViolation::DuplicatePort { .. })
    ));
    let signal = PortSpec::signal_out(StableId::new("probe").expect("id"), NonZeroUsize::MIN);
    let mixed = PortLayout::new(&[component.ports()[0].clone(), signal]).expect("layout");
    assert_eq!((mixed.input_len(), mixed.output_len()), (1, 2));
    assert!(mixed.input_range(1).is_empty());
}
