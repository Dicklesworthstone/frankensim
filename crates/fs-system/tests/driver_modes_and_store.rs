//! Direct verification of beads rc-1109 .5.4 (state store, accelerator &
//! step-size rollback, cache non-observability, snapshots) and .5.5 (all solve
//! modes by tier, monolithic StageResidual through fallible fs-solver adapter,
//! typed budget partial outcomes, cancellation at every phase, per-step JSONL).

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;
use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::identity::CanonicalLimits;
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotAlgorithmIdV2,
    SnapshotBudgetStateIdV2, SnapshotDeterminismV2, SnapshotExecutionFingerprintIdV2,
    SnapshotLimitsV2, SnapshotPauseRequestIdV2, SnapshotProvenanceIdV2, SnapshotRngCounterIdV2,
};
use fs_exec::{Budget, CancelGate, Cx, DrainTracker, ExecMode, RunId, StreamKey};
use fs_ports::component::conformance::{Case, Failure, check_stage_residual};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{Damper, HeatCapacity, Mass, Motion, Spring};
use fs_system::{
    AcceptedStateCache, AdaptiveConfig, AitkenConfig, BlockSolveMode, BlockSolverPolicy,
    BudgetKind, ErrorNorm, ErrorTolerance, FixedPointMethod, IqnIlsConfig, MonolithicConfig,
    NodeId, Phase, RouteTemperature, StateCacheKey, StepBudgets, StepVerdict, System,
    SystemBuilder, SystemError, SystemImage,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn id(s: &str) -> StableId {
    StableId::new(s).expect("valid id")
}

fn mech_port(name: &str, role: PortRole) -> PortSpec {
    PortSpec::power(
        id(name),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::MIN,
    )
    .expect("mech port")
}

fn heat_port(name: &str, role: PortRole) -> PortSpec {
    PortSpec::heat(id(name), role).expect("heat port")
}

fn roundoff(scale: f64) -> f64 {
    64.0 * f64::EPSILON * scale.max(1.0)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StageOwnerRefusal(u32);

/// Point mass implementing Component + PortTangent + StageResidual (exact midpoint).
/// State = velocity v [m/s]. Stage unknown z = [v_new].
struct StageMass {
    mass_kg: f64,
    ports: [PortSpec; 1],
    refuse_residual: Rc<Cell<Option<u32>>>,
    refuse_jvp: Rc<Cell<Option<u32>>>,
    cancel_in_propose: Rc<Cell<bool>>,
    refuse_admit: Rc<Cell<bool>>,
}

impl StageMass {
    fn new(mass_kg: f64) -> Self {
        Self {
            mass_kg,
            ports: [mech_port("p", PortRole::EffortIn)],
            refuse_residual: Rc::new(Cell::new(None)),
            refuse_jvp: Rc::new(Cell::new(None)),
            cancel_in_propose: Rc::new(Cell::new(false)),
            refuse_admit: Rc::new(Cell::new(false)),
        }
    }
}

impl Component for StageMass {
    type State = f64;
    type Error = StageOwnerRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, v: &f64) -> Result<f64, StageOwnerRefusal> {
        Ok(0.5 * self.mass_kg * v * v)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, StageOwnerRefusal> {
        if self.cancel_in_propose.get() {
            cx.cancel_gate().request();
        }
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        let f = inputs[0];
        *new = old + dt_s * f / self.mass_kg;
        let v_mid = 0.5 * (old + *new);
        outputs[0] = v_mid;
        let d_e = 0.5 * self.mass_kg * (*new * *new - old * old);
        let scale = 0.5 * self.mass_kg * (old * old + *new * *new) + (dt_s * f * v_mid).abs();
        Ok(Proposal::Step(Receipt::lossless(d_e, roundoff(scale))))
    }

    fn admit(&self, _old: &f64, _new: &f64, _receipt: &Receipt) -> Result<(), StageOwnerRefusal> {
        if self.refuse_admit.get() {
            Err(StageOwnerRefusal(77))
        } else {
            Ok(())
        }
    }
}

impl PortTangent for StageMass {
    fn jvp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        d_outputs[0] = 0.5 * dt_s / self.mass_kg * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        bar_inputs[0] = 0.5 * dt_s / self.mass_kg * bar_outputs[0];
        Ok(())
    }
}

impl StageResidual for StageMass {
    fn unknowns(&self) -> usize {
        1
    }

    fn initial_guess(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), StageOwnerRefusal> {
        z[0] = old + dt_s * inputs[0] / self.mass_kg;
        Ok(())
    }

    fn residual(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        if self.cancel_in_propose.get() {
            cx.cancel_gate().request();
        }
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        if let Some(code) = self.refuse_residual.get() {
            return Err(StageOwnerRefusal(code));
        }
        r[0] = self.mass_kg * (z[0] - old) - dt_s * inputs[0];
        outputs[0] = 0.5 * (old + z[0]);
        Ok(())
    }

    fn jvp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        if let Some(code) = self.refuse_jvp.get() {
            return Err(StageOwnerRefusal(code));
        }
        dr[0] = self.mass_kg * dz[0];
        d_outputs[0] = 0.5 * dz[0];
        Ok(())
    }

    fn vjp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        bar_z[0] = self.mass_kg * bar_r[0] + 0.5 * bar_outputs[0];
        Ok(())
    }

    fn accept(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut f64,
    ) -> Result<Receipt, StageOwnerRefusal> {
        *new = z[0];
        let v_mid = 0.5 * (old + *new);
        let d_e = 0.5 * self.mass_kg * (*new * *new - old * old);
        let scale =
            0.5 * self.mass_kg * (old * old + *new * *new) + (dt_s * inputs[0] * v_mid).abs();
        Ok(Receipt::lossless(d_e, roundoff(scale)))
    }
}

/// Nonlinear Duffing spring: V(x) = 0.5 k x^2 + 0.25 beta x^4.
/// Energy-exact discrete-gradient average force:
/// F_avg(x0, x1) = 0.5 k (x0 + x1) + 0.25 beta (x0 + x1) (x0^2 + x1^2).
/// Two ports: `a` and `b` (FlowIn -> EffortOut). Extension rate = v_a - v_b.
struct StageDuffingSpring {
    k: f64,
    beta: f64,
    ports: [PortSpec; 2],
}

impl StageDuffingSpring {
    fn new(k: f64, beta: f64) -> Self {
        Self {
            k,
            beta,
            ports: [
                mech_port("a", PortRole::FlowIn),
                mech_port("b", PortRole::FlowIn),
            ],
        }
    }

    fn force(&self, x0: f64, x1: f64) -> f64 {
        let s = x0 + x1;
        0.5 * self.k * s + 0.25 * self.beta * s * (x0 * x0 + x1 * x1)
    }

    fn dforce_dx1(&self, x0: f64, x1: f64) -> f64 {
        0.5 * self.k + 0.25 * self.beta * (x0 * x0 + 2.0 * x0 * x1 + 3.0 * x1 * x1)
    }
}

impl Component for StageDuffingSpring {
    type State = f64;
    type Error = StageOwnerRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, x: &f64) -> Result<f64, StageOwnerRefusal> {
        let x2 = x * x;
        Ok(0.5 * self.k * x2 + 0.25 * self.beta * x2 * x2)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        let v_rel = inputs[0] - inputs[1];
        *new = old + dt_s * v_rel;
        let f = self.force(*old, *new);
        outputs[0] = f;
        outputs[1] = -f;
        let d_e = self.stored_energy_j(new)? - self.stored_energy_j(old)?;
        let scale =
            self.stored_energy_j(old)? + self.stored_energy_j(new)? + (dt_s * f * v_rel).abs();
        Ok(Proposal::Step(Receipt::lossless(d_e, roundoff(scale))))
    }
}

impl PortTangent for StageDuffingSpring {
    fn jvp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        let x1 = old + dt_s * (inputs[0] - inputs[1]);
        let df = self.dforce_dx1(*old, x1) * dt_s * (d_inputs[0] - d_inputs[1]);
        d_outputs[0] = df;
        d_outputs[1] = -df;
        Ok(())
    }

    fn vjp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        let x1 = old + dt_s * (inputs[0] - inputs[1]);
        let g = self.dforce_dx1(*old, x1) * dt_s * (bar_outputs[0] - bar_outputs[1]);
        bar_inputs[0] = g;
        bar_inputs[1] = -g;
        Ok(())
    }
}

impl StageResidual for StageDuffingSpring {
    fn unknowns(&self) -> usize {
        1
    }

    fn initial_guess(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), StageOwnerRefusal> {
        z[0] = old + dt_s * (inputs[0] - inputs[1]);
        Ok(())
    }

    fn residual(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        r[0] = z[0] - old - dt_s * (inputs[0] - inputs[1]);
        let f = self.force(*old, z[0]);
        outputs[0] = f;
        outputs[1] = -f;
        Ok(())
    }

    fn jvp_z(
        &self,
        old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        dr[0] = dz[0];
        let df = self.dforce_dx1(*old, z[0]) * dz[0];
        d_outputs[0] = df;
        d_outputs[1] = -df;
        Ok(())
    }

    fn vjp_z(
        &self,
        old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        bar_z[0] = bar_r[0] + self.dforce_dx1(*old, z[0]) * (bar_outputs[0] - bar_outputs[1]);
        Ok(())
    }

    fn accept(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut f64,
    ) -> Result<Receipt, StageOwnerRefusal> {
        *new = z[0];
        let f = self.force(*old, *new);
        let d_e = self.stored_energy_j(new)? - self.stored_energy_j(old)?;
        let scale = self.stored_energy_j(old)?
            + self.stored_energy_j(new)?
            + (dt_s * f * (inputs[0] - inputs[1])).abs();
        Ok(Receipt::lossless(d_e, roundoff(scale)))
    }
}

/// Viscous damper with heat outlet and no internal stage unknowns. Its
/// algebraic force and heat still close the coupled interface equations.
struct StageDamper {
    c: f64,
    ports: [PortSpec; 3],
}

impl StageDamper {
    fn new(c: f64) -> Self {
        Self {
            c,
            ports: [
                mech_port("a", PortRole::FlowIn),
                mech_port("b", PortRole::FlowIn),
                heat_port("heat", PortRole::EffortIn),
            ],
        }
    }
}

impl Component for StageDamper {
    type State = ();
    type Error = StageOwnerRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, StageOwnerRefusal> {
        Ok(0.0)
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        let dv = inputs[0] - inputs[1];
        let t_k = inputs[2];
        if !(t_k.is_finite() && t_k > 0.0) {
            return Err(StageOwnerRefusal(3));
        }
        let f = self.c * dv;
        let p_diss = self.c * dv * dv;
        outputs[0] = f;
        outputs[1] = -f;
        outputs[2] = -p_diss;
        let diss_j = dt_s * p_diss;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: diss_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: diss_j / t_k,
            residual: 0.0,
            tolerance: 0.0,
            energy_tolerance_j: roundoff(diss_j),
        }))
    }
}

impl PortTangent for StageDamper {
    fn jvp(
        &self,
        (): &(),
        _dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        let dv = inputs[0] - inputs[1];
        let ddv = d_inputs[0] - d_inputs[1];
        d_outputs[0] = self.c * ddv;
        d_outputs[1] = -self.c * ddv;
        d_outputs[2] = -2.0 * self.c * dv * ddv;
        Ok(())
    }

    fn vjp(
        &self,
        (): &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        let dv = inputs[0] - inputs[1];
        let g = self.c * (bar_outputs[0] - bar_outputs[1]) - 2.0 * self.c * dv * bar_outputs[2];
        bar_inputs[0] = g;
        bar_inputs[1] = -g;
        bar_inputs[2] = 0.0;
        Ok(())
    }
}

impl StageResidual for StageDamper {
    fn unknowns(&self) -> usize {
        0
    }

    fn initial_guess(
        &self,
        (): &(),
        _dt_s: f64,
        _inputs: &[f64],
        _z: &mut [f64],
    ) -> Result<(), StageOwnerRefusal> {
        Ok(())
    }

    fn residual(
        &self,
        (): &(),
        _dt_s: f64,
        inputs: &[f64],
        _z: &[f64],
        _r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        let dv = inputs[0] - inputs[1];
        let t_k = inputs[2];
        if !(t_k.is_finite() && t_k > 0.0) {
            return Err(StageOwnerRefusal(3));
        }
        let force = self.c * dv;
        outputs[0] = force;
        outputs[1] = -force;
        outputs[2] = -force * dv;
        Ok(())
    }

    fn jvp_z(
        &self,
        (): &(),
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _dz: &[f64],
        _dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        d_outputs.fill(0.0);
        Ok(())
    }

    fn vjp_z(
        &self,
        (): &(),
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _bar_r: &[f64],
        _bar_outputs: &[f64],
        _bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        Ok(())
    }

    fn accept(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        _z: &[f64],
        (): &mut (),
    ) -> Result<Receipt, StageOwnerRefusal> {
        let dv = inputs[0] - inputs[1];
        let t_k = inputs[2];
        let diss_j = dt_s * self.c * dv * dv;
        Ok(Receipt {
            stored_change_j: 0.0,
            dissipated_j: diss_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: diss_j / t_k,
            residual: 0.0,
            tolerance: 0.0,
            energy_tolerance_j: roundoff(diss_j),
        })
    }
}

/// Thermal capacity implementing Component + PortTangent + StageResidual.
struct StageThermalBody {
    c_j_k: f64,
    ports: [PortSpec; 1],
}

impl StageThermalBody {
    fn new(c_j_k: f64) -> Self {
        Self {
            c_j_k,
            ports: [heat_port("heat", PortRole::FlowIn)],
        }
    }
}

impl Component for StageThermalBody {
    type State = f64;
    type Error = StageOwnerRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, t_k: &f64) -> Result<f64, StageOwnerRefusal> {
        Ok(self.c_j_k * *t_k)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        let q = inputs[0];
        *new = old + dt_s * q / self.c_j_k;
        if !(*new > 0.0) {
            return Err(StageOwnerRefusal(4));
        }
        let t_mid = 0.5 * (old + *new);
        outputs[0] = t_mid;
        let d_e = self.c_j_k * (*new - old);
        Ok(Proposal::Step(Receipt::lossless(
            d_e,
            roundoff(self.c_j_k * (*old + *new) + (dt_s * q).abs()),
        )))
    }
}

impl PortTangent for StageThermalBody {
    fn jvp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        d_outputs[0] = 0.5 * dt_s / self.c_j_k * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        bar_inputs[0] = 0.5 * dt_s / self.c_j_k * bar_outputs[0];
        Ok(())
    }
}

impl StageResidual for StageThermalBody {
    fn unknowns(&self) -> usize {
        1
    }

    fn initial_guess(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), StageOwnerRefusal> {
        z[0] = old + dt_s * inputs[0] / self.c_j_k;
        Ok(())
    }

    fn residual(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        if !(z[0] > 0.0) {
            return Err(StageOwnerRefusal(4));
        }
        r[0] = self.c_j_k * (z[0] - old) - dt_s * inputs[0];
        outputs[0] = 0.5 * (old + z[0]);
        Ok(())
    }

    fn jvp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        dr[0] = self.c_j_k * dz[0];
        d_outputs[0] = 0.5 * dz[0];
        Ok(())
    }

    fn vjp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageOwnerRefusal> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        bar_z[0] = self.c_j_k * bar_r[0] + 0.5 * bar_outputs[0];
        Ok(())
    }

    fn accept(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut f64,
    ) -> Result<Receipt, StageOwnerRefusal> {
        *new = z[0];
        let d_e = self.c_j_k * (*new - old);
        Ok(Receipt::lossless(
            d_e,
            roundoff(self.c_j_k * (*old + *new) + (dt_s * inputs[0]).abs()),
        ))
    }
}

fn build_two_component_oscillator() -> (System, NodeId, NodeId) {
    let mut b = SystemBuilder::new();
    let m = b
        .add_with_tangent_and_stage_residual("mass", StageMass::new(1.5), 0.9)
        .unwrap();
    let s = b
        .add_with_tangent_and_stage_residual("spring", StageDuffingSpring::new(18.0, 6.0), 0.25)
        .unwrap();
    b.bond(b.port(m, "p").unwrap(), b.port(s, "a").unwrap())
        .unwrap();
    b.boundary(b.port(s, "b").unwrap(), Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    (b.seal().unwrap(), m, s)
}

fn build_five_component_chain() -> (System, [NodeId; 5]) {
    let mut b = SystemBuilder::new();
    let m1 = b
        .add_with_tangent_and_stage_residual("m1", StageMass::new(1.2), 1.1)
        .unwrap();
    let s1 = b
        .add_with_tangent_and_stage_residual("s1", StageDuffingSpring::new(14.0, 4.0), 0.15)
        .unwrap();
    let m2 = b
        .add_with_tangent_and_stage_residual("m2", StageMass::new(2.0), -0.4)
        .unwrap();
    let d = b
        .add_with_tangent_and_stage_residual("damper", StageDamper::new(1.8), ())
        .unwrap();
    let body = b
        .add_with_tangent_and_stage_residual("body", StageThermalBody::new(25.0), 300.0)
        .unwrap();

    b.bond(b.port(m1, "p").unwrap(), b.port(s1, "a").unwrap())
        .unwrap();
    b.common_flow(
        b.port(m2, "p").unwrap(),
        &[b.port(s1, "b").unwrap(), b.port(d, "a").unwrap()],
    )
    .unwrap();
    b.boundary(b.port(d, "b").unwrap(), Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    let dh = b.port(d, "heat").unwrap();
    b.heat_route(
        dh,
        &[(b.port(body, "heat").unwrap(), 1.0)],
        None,
        RouteTemperature::SinkEquivalent,
    )
    .unwrap();
    let mut sys = b.seal().unwrap();
    assert!(sys.guess(dh, &[300.0]));
    (sys, [m1, s1, m2, d, body])
}

/// An algebraic signal relation. The accepted counter lets the test detect
/// publication of an algebraic owner's candidate after a later refusal.
struct AlgebraicSignal {
    gain: f64,
    offset: f64,
    ports: [PortSpec; 2],
    fault: Rc<Cell<u8>>,
}

impl AlgebraicSignal {
    fn new(gain: f64, offset: f64) -> Self {
        Self {
            gain,
            offset,
            ports: [
                PortSpec::signal_in(id("in"), NonZeroUsize::MIN),
                PortSpec::signal_out(id("out"), NonZeroUsize::MIN),
            ],
            fault: Rc::new(Cell::new(0)),
        }
    }
}

impl Component for AlgebraicSignal {
    type State = usize;
    type Error = StageOwnerRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &usize) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &usize,
        _dt_s: f64,
        inputs: &[f64],
        new: &mut usize,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        *new = old + 1;
        outputs[0] = self.offset + self.gain * inputs[0];
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

impl StageResidual for AlgebraicSignal {
    fn unknowns(&self) -> usize {
        0
    }

    fn initial_guess(
        &self,
        _old: &usize,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &mut [f64],
    ) -> Result<(), Self::Error> {
        Ok(())
    }

    fn residual(
        &self,
        old: &usize,
        dt_s: f64,
        inputs: &[f64],
        _z: &[f64],
        _r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        if self.fault.get() == 1 {
            return Err(StageOwnerRefusal(42));
        }
        if self.fault.get() == 3 {
            cx.cancel_gate().request();
        }
        let mut candidate = *old;
        self.propose(old, dt_s, inputs, &mut candidate, outputs, cx)?;
        Ok(())
    }

    fn jvp_z(
        &self,
        _old: &usize,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _dz: &[f64],
        _dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))?;
        d_outputs.fill(if self.fault.get() == 4 { 1.0 } else { 0.0 });
        Ok(())
    }

    fn vjp_z(
        &self,
        _old: &usize,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _bar_r: &[f64],
        _bar_outputs: &[f64],
        _bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        cx.checkpoint().map_err(|_| StageOwnerRefusal(999))
    }

    fn accept(
        &self,
        old: &usize,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        new: &mut usize,
    ) -> Result<Receipt, Self::Error> {
        if self.fault.get() == 2 {
            return Err(StageOwnerRefusal(84));
        }
        *new = old + 1;
        Ok(Receipt::lossless(0.0, 0.0))
    }
}

#[test]
fn zero_unknown_stage_conformance_checks_outputs_and_acceptance() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        check_stage_residual(
            &Damper::new(3.0).unwrap(),
            &Case {
                state: &(),
                stale: &(),
                dt_s: 0.01,
                inputs: &[0.8, -0.4, 300.0],
            },
            &[],
            1.0e-7,
            1.0e-12,
            cx,
        )
        .expect("native zero-unknown damper conforms");

        let signal = AlgebraicSignal::new(0.25, 1.0);
        let case = Case {
            state: &10,
            stale: &99,
            dt_s: 0.01,
            inputs: &[2.0],
        };
        check_stage_residual(&signal, &case, &[], 1.0e-7, 1.0e-12, cx)
            .expect("algebraic owner with an accepted counter conforms");
        signal.fault.set(4);
        assert!(matches!(
            check_stage_residual(&signal, &case, &[], 1.0e-7, 1.0e-12, cx),
            Err(Failure::Law {
                law: "jvp_z output derivative matches finite differences",
                ..
            })
        ));
    });
}

#[test]
fn monolithic_all_algebraic_loop_solves_and_refuses_atomically() {
    let policies = [
        BlockSolverPolicy::Auto,
        BlockSolverPolicy::ByTier {
            fixed_point: FixedPointMethod::IqnIls(IqnIlsConfig::default()),
            monolithic: MonolithicConfig::default(),
        },
        BlockSolverPolicy::MonolithicNewtonKrylov(MonolithicConfig::default()),
    ];
    for policy in policies {
        let mut b = SystemBuilder::new();
        let a = b
            .add_with_stage_residual("a", AlgebraicSignal::new(0.25, 1.0), 10)
            .unwrap();
        let signal_b = AlgebraicSignal::new(0.5, 2.0);
        let fault = Rc::clone(&signal_b.fault);
        let other = b.add_with_stage_residual("b", signal_b, 20).unwrap();
        let ao = b.port(a, "out").unwrap();
        let bo = b.port(other, "out").unwrap();
        b.signal(ao, b.port(other, "in").unwrap()).unwrap();
        b.signal(bo, b.port(a, "in").unwrap()).unwrap();
        let mut sys = b.seal().unwrap();
        sys.set_solver_policy(policy);
        assert_eq!(
            sys.block_modes(),
            vec![BlockSolveMode::MonolithicNewtonKrylov]
        );
        let plan = sys.execution_plan();
        assert_eq!(plan.blocks[0].stage_unknowns, 0);
        assert_eq!(plan.blocks[0].kept_unknowns, 2);
        assert_eq!(plan.blocks[0].tear_variables.len(), 2);
        assert!(!plan.explain().contains("missing StageResidual"));

        for code in 1..=3 {
            fault.set(code);
            let gate = CancelGate::new_clock_free();
            let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("injected refusal");
            match (code, err) {
                (1 | 2, SystemError::Component { node, error, .. }) => {
                    assert_eq!(node, "b");
                    let expected = if code == 1 { 42 } else { 84 };
                    assert_eq!(
                        error.downcast_ref::<StageOwnerRefusal>(),
                        Some(&StageOwnerRefusal(expected))
                    );
                }
                (3, SystemError::Cancelled) => {}
                (_, other) => panic!("incorrect refusal: {other:?}"),
            }
            assert_eq!(*sys.state::<usize>(a).unwrap(), 10);
            assert_eq!(*sys.state::<usize>(other).unwrap(), 20);
            assert_eq!(sys.steps(), 0);
            assert_eq!(sys.time_s(), 0.0);
            assert_eq!(sys.last_output(ao), &[0.0]);
            assert_eq!(sys.last_output(bo), &[0.0]);
            assert_eq!(sys.ledger().energy_residual_j(), 0.0);
        }

        fault.set(0);
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| sys.step(0.01, cx)).expect("clean retry");
        // Independent solution of y_a = 1 + y_b/4, y_b = 2 + y_a/2.
        assert!((sys.last_output(ao)[0] - 12.0 / 7.0).abs() < 1.0e-13);
        assert!((sys.last_output(bo)[0] - 20.0 / 7.0).abs() < 1.0e-13);
        assert_eq!(*sys.state::<usize>(a).unwrap(), 11);
        assert_eq!(*sys.state::<usize>(other).unwrap(), 21);
        assert_eq!(sys.steps(), 1);
        assert_eq!(sys.last_receipts().len(), 2);
        assert_eq!(sys.ledger().energy_residual_j(), 0.0);
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn two_component_and_five_component_graphs_converge_identically_across_all_modes() {
    let gate = CancelGate::new_clock_free();
    let policies = [
        (
            "interface_newton",
            BlockSolverPolicy::InterfaceNewton,
            BlockSolveMode::InterfaceNewton,
        ),
        (
            "fixed_point_aitken",
            BlockSolverPolicy::FixedPointAitken(AitkenConfig::default()),
            BlockSolveMode::FixedPointAitken,
        ),
        (
            "fixed_point_iqn_ils",
            BlockSolverPolicy::FixedPointIqnIls(IqnIlsConfig::default()),
            BlockSolveMode::FixedPointIqnIls,
        ),
        (
            "monolithic_nk",
            BlockSolverPolicy::MonolithicNewtonKrylov(MonolithicConfig::default()),
            BlockSolveMode::MonolithicNewtonKrylov,
        ),
    ];

    // 1. Two-component nonlinear Duffing oscillator across all 4 modes.
    let mut ref_two: Option<(f64, f64)> = None;
    for (name, policy, expected_mode) in policies {
        let (mut sys, m, s) = build_two_component_oscillator();
        sys.set_solver_policy(policy);
        assert_eq!(sys.block_modes(), vec![expected_mode]);
        with_cx(&gate, |cx| {
            for _ in 0..20 {
                sys.step(0.01, cx)
                    .unwrap_or_else(|e| panic!("{name} 2-comp step failed: {e}"));
            }
        });
        let v = *sys.state::<f64>(m).unwrap();
        let x = *sys.state::<f64>(s).unwrap();
        assert!(
            sys.ledger().energy_residual_j().abs() <= 1.0e-11,
            "{name} 2-comp energy defect {:e}",
            sys.ledger().energy_residual_j()
        );
        if let Some((rv, rx)) = ref_two {
            assert!(
                (v - rv).abs() < 1.0e-10 && (x - rx).abs() < 1.0e-10,
                "{name} diverged from reference: ({v}, {x}) vs ({rv}, {rx})"
            );
        } else {
            ref_two = Some((v, x));
        }
    }

    // 2. Five-component coupled mechanical-dissipative-thermal chain across all 4 modes.
    let mut ref_five: Option<[f64; 4]> = None;
    for (name, policy, expected_mode) in policies {
        let (mut sys, [m1, s1, m2, _d, body]) = build_five_component_chain();
        sys.set_solver_policy(policy);
        assert!(sys.block_modes().contains(&expected_mode));
        with_cx(&gate, |cx| {
            for _ in 0..20 {
                sys.step(0.01, cx)
                    .unwrap_or_else(|e| panic!("{name} 5-comp step failed: {e}"));
            }
        });
        let vals = [
            *sys.state::<f64>(m1).unwrap(),
            *sys.state::<f64>(s1).unwrap(),
            *sys.state::<f64>(m2).unwrap(),
            *sys.state::<f64>(body).unwrap(),
        ];
        assert!(vals[3] > 300.0, "{name} damper must heat thermal body");
        assert!(
            sys.ledger().energy_residual_j().abs() <= 1.0e-10,
            "{name} 5-comp energy defect {:e}",
            sys.ledger().energy_residual_j()
        );
        assert!(sys.ledger().entropy_generated_j_per_k() > 0.0);
        if let Some(rv) = ref_five {
            for (i, (&a, &b)) in vals.iter().zip(&rv).enumerate() {
                assert!(
                    (a - b).abs() < 1.0e-9,
                    "{name} 5-comp coord {i}: {a} vs {b}"
                );
            }
        } else {
            ref_five = Some(vals);
        }
    }

    // 3. Verify ByTier selects OrderedSweep on acyclic, FixedPointIqnIls on plain Component-only
    //    loops, InterfaceNewton on PortTangent loops, and MonolithicNewtonKrylov on StageResidual.
    let mut b = SystemBuilder::new();
    let m = b
        .add(
            "m",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.0,
            },
        )
        .unwrap();
    let sp = b.add("s", Spring::new(12.0, 0.0).unwrap(), 0.2).unwrap();
    let h = b.add("h", HeatCapacity::new(5.0).unwrap(), 300.0).unwrap();
    b.bond(b.port(m, "p").unwrap(), b.port(sp, "a").unwrap())
        .unwrap();
    b.boundary(b.port(sp, "b").unwrap(), Box::new(|_, _, o| o[0] = 0.0))
        .unwrap();
    b.boundary(b.port(h, "heat").unwrap(), Box::new(|_, _, o| o[0] = 2.0))
        .unwrap();
    let mut tiered = b.seal().unwrap();
    tiered.set_solver_policy(BlockSolverPolicy::ByTier {
        fixed_point: FixedPointMethod::IqnIls(IqnIlsConfig::default()),
        monolithic: MonolithicConfig::default(),
    });
    let modes = tiered.block_modes();
    assert!(modes.contains(&BlockSolveMode::OrderedSweep));
    assert!(modes.contains(&BlockSolveMode::FixedPointIqnIls));
    with_cx(&gate, |cx| {
        for _ in 0..10 {
            tiered.step(0.01, cx).expect("ByTier mixed sweep + IQN-ILS");
        }
    });
    assert!(tiered.ledger().energy_residual_j().abs() <= 1.0e-12);
}

#[test]
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn accelerator_and_step_size_memory_roll_back_bitwise_and_round_trip_snapshots() {
    let gate = CancelGate::new_clock_free();

    fn build_native_oscillator() -> (System, NodeId, NodeId) {
        let mut b = SystemBuilder::new();
        let m = b
            .add_with_tangent(
                "m",
                Mass::new(1.2).unwrap(),
                Motion {
                    position_m: 0.2,
                    velocity_m_s: 0.8,
                },
            )
            .unwrap();
        let s = b
            .add_with_tangent("s", Spring::new(20.0, 0.0).unwrap(), -0.2)
            .unwrap();
        let d = b
            .add_with_tangent("d", Damper::new(0.6).unwrap(), ())
            .unwrap();
        b.common_flow(
            b.port(m, "p").unwrap(),
            &[b.port(s, "a").unwrap(), b.port(d, "a").unwrap()],
        )
        .unwrap();
        b.boundary(b.port(s, "b").unwrap(), Box::new(|_, _, o| o[0] = 0.0))
            .unwrap();
        b.boundary(b.port(d, "b").unwrap(), Box::new(|_, _, o| o[0] = 0.0))
            .unwrap();
        b.boundary(b.port(d, "heat").unwrap(), Box::new(|_, _, o| o[0] = 300.0))
            .unwrap();
        (b.seal().unwrap(), m, s)
    }

    // 1. Run with IQN-ILS + adaptive stepping so both accelerator history and suggested_dt_s are populated.
    let (mut sys, m, s) = build_native_oscillator();
    sys.set_solver_policy(BlockSolverPolicy::FixedPointIqnIls(IqnIlsConfig::default()));
    with_cx(&gate, |cx| {
        for _ in 0..5 {
            sys.step(0.01, cx).unwrap();
        }
        let adaptive_rep = sys
            .step_adaptive(
                0.02,
                AdaptiveConfig {
                    order: 2,
                    min_step_s: 1.0e-12,
                    max_attempts: 16,
                },
                |old, coarse, fine| {
                    let mut n = ErrorNorm::new();
                    let a = old.state::<Motion>(m).unwrap();
                    let b = coarse.state::<Motion>(m).unwrap();
                    let c = fine.state::<Motion>(m).unwrap();
                    n.include(
                        a.velocity_m_s,
                        b.velocity_m_s,
                        c.velocity_m_s,
                        ErrorTolerance {
                            absolute: 1.0e-4,
                            relative: 0.0,
                        },
                    );
                    n.value()
                },
                cx,
            )
            .unwrap();
        assert_eq!(
            sys.suggested_dt_s().map(f64::to_bits),
            Some(adaptive_rep.suggested_dt_s.to_bits())
        );
    });
    assert!(sys.block_accelerator_history_len(0).unwrap() > 0);
    let saved_suggested_bits = sys.suggested_dt_s().unwrap().to_bits();
    let saved_omega_bits = sys.block_relaxation_omega(0).unwrap().to_bits();
    let saved_hist_len = sys.block_accelerator_history_len(0).unwrap();
    let saved_gen = sys.accepted_generation();

    // Capture both in-memory Checkpoint and durable SolverStateV2 SystemImage!
    let cp = sys.checkpoint().unwrap();
    assert_eq!(cp.suggested_dt_s().unwrap().to_bits(), saved_suggested_bits);
    assert_eq!(
        cp.block_relaxation_omega(0).unwrap().to_bits(),
        saved_omega_bits
    );

    let drain_gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(RunId(99), &drain_gate);
    let worker = tracker.register_worker().unwrap();
    drain_gate.request();
    drop(worker);
    let pause = PausedSnapshotBoundaryV2::from_drain_report(
        tracker.finalize().unwrap(),
        SnapshotPauseRequestIdV2::from_bytes([0x42; 32]),
        1,
    );
    let ctx = ExpectedResumeContextV2::for_paused_state::<SystemImage>(
        SnapshotAlgorithmIdV2::from_bytes([0x11; 32]),
        1,
        sys.problem_id(),
        SnapshotRngCounterIdV2::from_bytes([0; 32]),
        SnapshotDeterminismV2::Deterministic,
        SnapshotExecutionFingerprintIdV2::from_bytes([0x22; 32]),
        SnapshotBudgetStateIdV2::from_bytes([0x33; 32]),
        SnapshotProvenanceIdV2::from_bytes([0x44; 32]),
        pause,
    );
    let limits = SnapshotLimitsV2::new(
        1 << 20,
        4096,
        CanonicalLimits::new(16384, 4096, 32, 32, 64),
        1 << 16,
        1 << 20,
        4096,
    );
    let sealed = sys.snapshot(&ctx, limits, || false).unwrap();

    // 2. Budget-exhausted trial must leave accepted controller, history, and generation bitwise untouched!
    sys.set_budgets(StepBudgets {
        max_evaluations: Some(1),
        ..StepBudgets::default()
    });
    with_cx(&gate, |cx| {
        assert!(matches!(
            sys.step(0.01, cx),
            Err(SystemError::BudgetExhausted(_))
        ));
    });
    assert_eq!(
        sys.suggested_dt_s().unwrap().to_bits(),
        saved_suggested_bits
    );
    assert_eq!(
        sys.block_relaxation_omega(0).unwrap().to_bits(),
        saved_omega_bits
    );
    assert_eq!(
        sys.block_accelerator_history_len(0).unwrap(),
        saved_hist_len
    );
    assert_eq!(sys.accepted_generation(), saved_gen);

    // 3. Continue uninterrupted for 10 steps after clearing budget.
    sys.set_budgets(StepBudgets::default());
    with_cx(&gate, |cx| {
        for _ in 0..10 {
            sys.step(0.01, cx).unwrap();
        }
    });
    let expected_v = sys.state::<Motion>(m).unwrap().velocity_m_s.to_bits();
    let expected_x = sys.state::<f64>(s).unwrap().to_bits();
    let expected_omega = sys.block_relaxation_omega(0).unwrap().to_bits();

    // 4. Restore Checkpoint -> replay 10 steps -> BITWISE IDENTICAL!
    sys.restore(&cp).unwrap();
    assert!(sys.accepted_generation() > saved_gen);
    with_cx(&gate, |cx| {
        for _ in 0..10 {
            sys.step(0.01, cx).unwrap();
        }
    });
    assert_eq!(
        sys.state::<Motion>(m).unwrap().velocity_m_s.to_bits(),
        expected_v
    );
    assert_eq!(sys.state::<f64>(s).unwrap().to_bits(), expected_x);
    assert_eq!(
        sys.block_relaxation_omega(0).unwrap().to_bits(),
        expected_omega
    );

    // 5. Resume durable snapshot into brand-new System -> replay 10 steps -> BITWISE IDENTICAL!
    let opened =
        SystemImage::unseal_v2_expected(sealed.bytes(), &sealed.expectation(), limits, || false)
            .unwrap();
    let (mut resumed, rm, rs) = build_native_oscillator();
    resumed.set_solver_policy(BlockSolverPolicy::FixedPointIqnIls(IqnIlsConfig::default()));
    resumed.resume(&opened).unwrap();
    assert_eq!(
        resumed.suggested_dt_s().unwrap().to_bits(),
        saved_suggested_bits
    );
    assert_eq!(
        resumed.block_relaxation_omega(0).unwrap().to_bits(),
        saved_omega_bits
    );
    assert_eq!(
        resumed.block_accelerator_history_len(0).unwrap(),
        saved_hist_len
    );
    with_cx(&gate, |cx| {
        for _ in 0..10 {
            resumed.step(0.01, cx).unwrap();
        }
    });
    assert_eq!(
        resumed.state::<Motion>(rm).unwrap().velocity_m_s.to_bits(),
        expected_v
    );
    assert_eq!(resumed.state::<f64>(rs).unwrap().to_bits(), expected_x);
    assert_eq!(
        resumed.block_relaxation_omega(0).unwrap().to_bits(),
        expected_omega
    );

    // 6. Cache non-observability + AcceptedStateCache generation keying!
    let (mut warm, wm, ws) = build_native_oscillator();
    let (mut cold, cm, cs) = build_native_oscillator();
    let mut cache: AcceptedStateCache<f64> = AcceptedStateCache::new();
    with_cx(&gate, |cx| {
        for k in 0..15 {
            let key = StateCacheKey::new(warm.accepted_generation(), 0.005);
            assert!(
                cache.get(key).is_none(),
                "generation advanced so prior entry must miss"
            );
            cache.insert(key, f64::from(k));
            assert_eq!(cache.get(key), Some(&f64::from(k)));

            warm.step(0.005, cx).unwrap();
            cold.invalidate_caches();
            cold.step(0.005, cx).unwrap();
            assert_eq!(
                warm.state::<Motion>(wm).unwrap().velocity_m_s.to_bits(),
                cold.state::<Motion>(cm).unwrap().velocity_m_s.to_bits()
            );
            assert_eq!(
                warm.state::<f64>(ws).unwrap().to_bits(),
                cold.state::<f64>(cs).unwrap().to_bits()
            );
        }
    });
}

#[test]
fn monolithic_stage_residual_propagates_owner_refusal_without_nan_and_cancellation_at_every_phase()
{
    let mut b = SystemBuilder::new();
    let mass = StageMass::new(1.5);
    let refuse_residual = Rc::clone(&mass.refuse_residual);
    let refuse_jvp = Rc::clone(&mass.refuse_jvp);
    let cancel_in_propose = Rc::clone(&mass.cancel_in_propose);
    let refuse_admit = Rc::clone(&mass.refuse_admit);
    let m = b.add_with_stage_residual("mass", mass, 0.75).unwrap();
    let s = b
        .add_with_stage_residual("spring", StageDuffingSpring::new(15.0, 3.0), 0.3)
        .unwrap();
    b.bond(b.port(m, "p").unwrap(), b.port(s, "a").unwrap())
        .unwrap();
    b.boundary(b.port(s, "b").unwrap(), Box::new(|_, _, o| o[0] = 0.0))
        .unwrap();
    let mut sys = b.seal().unwrap();
    assert_eq!(
        sys.block_modes(),
        vec![BlockSolveMode::MonolithicNewtonKrylov]
    );

    let orig_v = sys.state::<f64>(m).unwrap().to_bits();
    let orig_x = sys.state::<f64>(s).unwrap().to_bits();

    // 1. Residual refusal propagates exact StageOwnerRefusal(42), never NaN!
    let gate = CancelGate::new_clock_free();
    refuse_residual.set(Some(42));
    let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("residual refusal");
    match err {
        SystemError::Component { node, phase, error } => {
            assert_eq!(node, "mass");
            assert_eq!(phase, Phase::Propose);
            assert_eq!(
                error.downcast_ref::<StageOwnerRefusal>(),
                Some(&StageOwnerRefusal(42))
            );
        }
        other => panic!("expected typed owner refusal, got {other:?}"),
    }
    assert_eq!(sys.state::<f64>(m).unwrap().to_bits(), orig_v);
    assert_eq!(sys.state::<f64>(s).unwrap().to_bits(), orig_x);
    assert_eq!(sys.steps(), 0);
    refuse_residual.set(None);

    // 2. JVP refusal propagates exact StageOwnerRefusal(84), never NaN!
    refuse_jvp.set(Some(84));
    let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("jvp refusal");
    match err {
        SystemError::Component { node, phase, error } => {
            assert_eq!(node, "mass");
            assert_eq!(phase, Phase::Propose);
            assert_eq!(
                error.downcast_ref::<StageOwnerRefusal>(),
                Some(&StageOwnerRefusal(84))
            );
        }
        other => panic!("expected typed jvp owner refusal, got {other:?}"),
    }
    assert_eq!(sys.state::<f64>(m).unwrap().to_bits(), orig_v);
    assert_eq!(sys.steps(), 0);
    refuse_jvp.set(None);

    // 3. Late admit refusal leaves all states untouched!
    refuse_admit.set(true);
    let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("admit refusal");
    assert!(matches!(
        err,
        SystemError::Component {
            phase: Phase::Admit,
            ..
        }
    ));
    assert_eq!(sys.state::<f64>(m).unwrap().to_bits(), orig_v);
    assert_eq!(sys.state::<f64>(s).unwrap().to_bits(), orig_x);
    refuse_admit.set(false);

    // 4. Mid-propose cancellation returns Cancelled and leaves state untouched!
    cancel_in_propose.set(true);
    let cancel_gate = CancelGate::new_clock_free();
    let err = with_cx(&cancel_gate, |cx| sys.step(0.01, cx)).expect_err("mid-propose cancel");
    assert!(matches!(err, SystemError::Cancelled));
    assert_eq!(sys.state::<f64>(m).unwrap().to_bits(), orig_v);
    assert_eq!(sys.steps(), 0);
    cancel_in_propose.set(false);

    // 5. Pre-cancelled context aborts immediately at Negotiate!
    let pre_gate = CancelGate::new_clock_free();
    pre_gate.request();
    let err = with_cx(&pre_gate, |cx| sys.step(0.01, cx)).expect_err("pre-cancel");
    assert!(matches!(err, SystemError::Cancelled));
    assert_eq!(sys.steps(), 0);

    // 6. Clean retry succeeds identically to a first attempt!
    with_cx(&gate, |cx| sys.step(0.01, cx)).expect("clean retry");
    assert_eq!(sys.steps(), 1);
    assert_eq!(sys.last_receipts().len(), 2);
}

#[test]
#[allow(clippy::result_large_err)]
fn budgets_return_typed_partial_outcomes_and_step_diagnosed_emits_valid_jsonl() {
    let gate = CancelGate::new_clock_free();
    let (mut sys, _nodes) = build_five_component_chain();

    // 1. Evaluations budget exhaustion -> typed PartialStepOutcome
    sys.set_budgets(StepBudgets {
        max_evaluations: Some(2),
        ..StepBudgets::default()
    });
    let (err, diag) = with_cx(&gate, |cx| sys.step_diagnosed(0.01, cx)).expect_err("eval budget");
    match err {
        SystemError::BudgetExhausted(partial) => {
            assert_eq!(partial.budget, BudgetKind::Evaluations);
            assert!(partial.spent_evaluations > 2);
        }
        other => panic!("expected BudgetExhausted(Evaluations), got {other:?}"),
    }
    assert_eq!(
        diag.verdict,
        StepVerdict::BudgetExhausted {
            budget: BudgetKind::Evaluations
        }
    );
    let jsonl = diag.to_jsonl();
    assert!(!jsonl.contains('\n'));
    assert!(jsonl.contains("\"verdict\":{\"budget_exhausted\":\"evaluations\"}"));
    assert_eq!(sys.steps(), 0);

    // 2. Work budget exhaustion -> typed PartialStepOutcome
    sys.set_budgets(StepBudgets {
        max_work_units: Some(1),
        ..StepBudgets::default()
    });
    let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("work budget");
    assert!(matches!(
        err,
        SystemError::BudgetExhausted(p) if p.budget == BudgetKind::Work
    ));
    assert_eq!(sys.steps(), 0);

    // 3. Wall-time budget exhaustion (0 ns) -> typed PartialStepOutcome
    sys.set_budgets(StepBudgets {
        max_wall_ns: Some(0),
        ..StepBudgets::default()
    });
    let err = with_cx(&gate, |cx| sys.step(0.01, cx)).expect_err("wall budget");
    assert!(matches!(
        err,
        SystemError::BudgetExhausted(p) if p.budget == BudgetKind::WallTime
    ));
    assert_eq!(sys.steps(), 0);

    // 4. Accepted diagnosed step emits complete JSONL with dt, mode, iterations,
    //    residual_norms, receipts (in declared component order), and verdict!
    sys.set_budgets(StepBudgets::default());
    sys.set_solver_policy(BlockSolverPolicy::FixedPointIqnIls(IqnIlsConfig::default()));
    let (report, diag) = with_cx(&gate, |cx| sys.step_diagnosed(0.01, cx)).expect("accepted step");
    assert_eq!(diag.step, 1);
    assert_eq!(diag.dt_s, report.dt_s);
    assert_eq!(diag.verdict, StepVerdict::Accepted);
    assert_eq!(diag.receipts.len(), 5);
    assert_eq!(
        diag.receipts
            .iter()
            .map(|r| r.node.as_str())
            .collect::<Vec<_>>(),
        vec!["m1", "s1", "m2", "damper", "body"]
    );
    assert_eq!(sys.last_receipts().len(), 5);
    let line = diag.to_jsonl();
    assert!(!line.contains('\n'), "JSONL must be single-line");
    for required in [
        "\"dt\":",
        "\"mode\":\"fixed_point_iqn_ils\"",
        "\"iterations\":",
        "\"residual_norms\":{\"l2\":",
        "\"linf\":",
        "\"receipts\":[",
        "\"verdict\":\"accepted\"",
    ] {
        assert!(line.contains(required), "JSONL missing {required}: {line}");
    }
}
