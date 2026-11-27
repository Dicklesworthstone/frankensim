//! Lock-free `LiveSource` component (`Component + PortTangent + StageResidual`)
//! and power-port lowering hook for interactive runtime inputs.

use core::num::NonZeroUsize;
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_exec::Cx;
use fs_ports::component::{
    Carrier, Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};
use fs_ports::{PortKind, StableId};
use fs_qty::Dims;
use fs_scenario::EntityId;
use fs_system::lowering::{
    InterfaceKindTag, InterfaceLoweringKey, LoweredInterfaceOutput, ProvenanceRecorder,
};
use fs_system::{LoweringError, LoweringRegistry};

use crate::descriptor::RuntimeRefusal;

/// Power-port card family registered by [`crate::AssemblyRuntime`] for interactive inputs.
pub const RUNTIME_INPUT_FAMILY: &str = "runtime-input";

/// Declarative specification of one live interactive input attached to a draft port.
#[derive(Debug, Clone, PartialEq)]
pub struct LiveInputSpec {
    /// Unique input name exposed through `inputs()`.
    pub name: String,
    /// Target entity owning `port`.
    pub target: EntityId,
    /// Target entity display name.
    pub target_name: String,
    /// Exposed port name on `target` (e.g. `"mechanics"`).
    pub port: String,
    /// Power domain for `LiveSource` attachment.
    pub domain: PortKind,
    /// Physical dimensions of the commanded value.
    pub dims: Dims,
    /// Initial value held before any `set_input` call.
    pub default_value: f64,
}

pub(crate) struct RegisteredInput {
    pub(crate) port: String,
    pub(crate) dims: Dims,
    pub(crate) default_value: f64,
    pub(crate) cell: Arc<AtomicU64>,
}

/// Lock-free interactive power source implementing `StageResidual` so untouched
/// constant sources stay on `BlockSolveMode::MonolithicNewtonKrylov` and match static CLI
/// lowering bit-for-bit (`to_bits()`).
///
/// `AssemblyRuntime` latches `pending -> live` once before each `System::step`;
/// callers constructing `LiveSource` directly must not mutate `cell` mid-step.
#[derive(Debug, Clone)]
pub struct LiveSource {
    ports: [PortSpec; 1],
    cell: Arc<AtomicU64>,
}

fn err_nf(msg: &'static str) -> RuntimeRefusal {
    RuntimeRefusal::new("non-finite-input", "step", msg)
}

impl LiveSource {
    /// Create a lock-free power-port source backed by `cell`.
    #[must_use]
    pub fn new(kind: PortKind, role: PortRole, cell: Arc<AtomicU64>) -> Self {
        let id = StableId::new("p").expect("static port id");
        let spec = PortSpec::power(id, kind, role, NonZeroUsize::MIN).expect("static power port");
        Self {
            ports: [spec],
            cell,
        }
    }

    fn value(&self) -> Result<f64, RuntimeRefusal> {
        let v = f64::from_bits(self.cell.load(Ordering::Relaxed));
        v.is_finite()
            .then_some(v)
            .ok_or_else(|| err_nf("live input is non-finite"))
    }

    fn book(&self, old: f64, dt_s: f64, input: f64) -> Result<(f64, f64, Receipt), RuntimeRefusal> {
        if !input.is_finite() {
            return Err(err_nf("port input non-finite"));
        }
        let val = self.value()?;
        let change_j = dt_s * val * input;
        let next = old + change_j;
        let tol = 8.0 * f64::EPSILON * (old.abs() + next.abs() + change_j.abs());
        Ok((val, next, Receipt::lossless(change_j, tol)))
    }
}

fn ck(cx: &Cx<'_>) -> Result<(), RuntimeRefusal> {
    cx.checkpoint()
        .map_err(|_| RuntimeRefusal::new("cancelled", "step", "cancelled"))
}

impl Component for LiveSource {
    type State = f64;
    type Error = RuntimeRefusal;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, state: &f64) -> Result<f64, RuntimeRefusal> {
        Ok(*state)
    }
    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        out: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, RuntimeRefusal> {
        ck(cx)?;
        let (val, next, r) = self.book(*old, dt_s, inputs[0])?;
        *new = next;
        out[0] = val;
        Ok(Proposal::Step(r))
    }
}

impl PortTangent for LiveSource {
    fn jvp(
        &self,
        _: &f64,
        _: f64,
        _: &[f64],
        _: &[f64],
        d_out: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), RuntimeRefusal> {
        ck(cx)?;
        d_out[0] = 0.0;
        Ok(())
    }
    fn vjp(
        &self,
        _: &f64,
        _: f64,
        _: &[f64],
        _: &[f64],
        bar_in: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), RuntimeRefusal> {
        ck(cx)?;
        bar_in[0] = 0.0;
        Ok(())
    }
}

impl StageResidual for LiveSource {
    fn unknowns(&self) -> usize {
        0
    }
    fn initial_guess(
        &self,
        _: &f64,
        _: f64,
        inputs: &[f64],
        _: &mut [f64],
    ) -> Result<(), RuntimeRefusal> {
        inputs
            .iter()
            .all(|x| x.is_finite())
            .then_some(())
            .ok_or_else(|| err_nf("port input non-finite"))
    }
    fn residual(
        &self,
        _: &f64,
        _: f64,
        _: &[f64],
        _: &[f64],
        _: &mut [f64],
        out: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), RuntimeRefusal> {
        ck(cx)?;
        out[0] = self.value()?;
        Ok(())
    }
    fn jvp_z(
        &self,
        _: &f64,
        _: f64,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: &mut [f64],
        d_out: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), RuntimeRefusal> {
        ck(cx)?;
        d_out[0] = 0.0;
        Ok(())
    }
    fn vjp_z(
        &self,
        _: &f64,
        _: f64,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), RuntimeRefusal> {
        ck(cx)
    }
    fn accept(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        _: &[f64],
        new: &mut f64,
        cx: &Cx<'_>,
    ) -> Result<Receipt, RuntimeRefusal> {
        ck(cx)?;
        let (_, next, r) = self.book(*old, dt_s, inputs[0])?;
        *new = next;
        Ok(r)
    }
}

pub(crate) fn dims_match_power_port(kind: PortKind, role: PortRole, dims: Dims) -> bool {
    let want = match (kind, role) {
        (PortKind::MechanicalForceVelocity, PortRole::EffortIn) => [1, 1, -2, 0, 0, 0],
        (PortKind::MechanicalForceVelocity, PortRole::FlowIn) => [1, 0, -1, 0, 0, 0],
        (PortKind::RotationalTorqueAngularVelocity, PortRole::EffortIn) => [2, 1, -2, 0, 0, 0],
        (PortKind::RotationalTorqueAngularVelocity, PortRole::FlowIn) => [0, 0, -1, 0, 0, 0],
        (PortKind::ElectricalVoltageCurrent, PortRole::EffortIn) => [2, 1, -3, 0, -1, 0],
        (PortKind::ElectricalVoltageCurrent, PortRole::FlowIn) => [0, 0, 0, 0, 1, 0],
        (PortKind::FluidPressureFlux, PortRole::EffortIn) => [-1, 1, -2, 0, 0, 0],
        (PortKind::FluidPressureFlux, PortRole::FlowIn) => [3, 0, -1, 0, 0, 0],
        _ => return false,
    };
    dims.0 == want
}

pub(crate) fn register_runtime_inputs(
    registry: &mut LoweringRegistry,
    cells_shared: Arc<BTreeMap<String, RegisteredInput>>,
) {
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::PowerPort, RUNTIME_INPUT_FAMILY),
        move |ctx, graph| {
            let reg =
                cells_shared
                    .get(&ctx.interface.name)
                    .ok_or_else(|| LoweringError::Invalid {
                        item: ctx.interface.name.clone(),
                        detail: "missing runtime input binding".to_owned(),
                    })?;
            let tp = graph
                .resolve_entity_port(
                    &ctx.draft.input().catalog,
                    ctx.interface.feature_b,
                    &reg.port,
                )
                .ok_or_else(|| LoweringError::MissingPort {
                    item: ctx.interface.name.clone(),
                    entity: ctx.interface.feature_b,
                    port: reg.port.clone(),
                })?;
            let (carrier, role) = graph.port_meta(tp).ok_or_else(|| LoweringError::Invalid {
                item: ctx.interface.name.clone(),
                detail: format!("port `{}` has no carrier metadata", reg.port),
            })?;
            let cell = Arc::clone(&reg.cell);
            let mut rec = ProvenanceRecorder::default();
            let mut created = Vec::new();
            match (carrier, role) {
                (Carrier::Power(kind), PortRole::EffortIn | PortRole::FlowIn) => {
                    if !dims_match_power_port(kind, role, reg.dims) {
                        return Err(LoweringError::Invalid {
                            item: ctx.interface.name.clone(),
                            detail: format!(
                                "input dimensions do not match `{}` port role",
                                reg.port
                            ),
                        });
                    }
                    let nname = format!("runtime.{}", ctx.interface.name);
                    let src_role = if role == PortRole::EffortIn {
                        PortRole::FlowIn
                    } else {
                        PortRole::EffortIn
                    };
                    let node = graph.add_with_tangent_and_residual(
                        &nname,
                        None,
                        LiveSource::new(kind, src_role, cell),
                        0.0,
                    )?;
                    let sp = graph.port(node, "p")?;
                    if role == PortRole::EffortIn {
                        graph.attach_common_flow(tp, sp, 1.0);
                    } else {
                        graph.attach_common_effort(tp, sp, 1.0);
                    }
                    created.push(nname);
                }
                (Carrier::Signal, PortRole::SignalIn) => {
                    graph.builder_mut().boundary(
                        tp,
                        Box::new(move |_, _, out| {
                            out[0] = f64::from_bits(cell.load(Ordering::Relaxed));
                        }),
                    )?;
                }
                _ => {
                    return Err(LoweringError::Invalid {
                        item: ctx.interface.name.clone(),
                        detail: format!("unsupported interactive carrier/role on `{}`", reg.port),
                    });
                }
            }
            rec.record_declared_field(
                "interactive_default",
                reg.default_value,
                reg.dims,
                &ctx.interface.name,
                "runtime.input",
            );
            let (property_receipts, parameter_sources) = rec.finish();
            Ok(LoweredInterfaceOutput {
                created_nodes: created,
                property_receipts,
                parameter_sources,
            })
        },
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_system::elements::{Force, Mass, Motion, SourceAccount, Spring, Waveform};
    use fs_system::{BlockSolveMode, SystemBuilder};

    #[test]
    fn power_port_role_dimension_guard_rejects_velocity_into_effort_port() {
        use PortKind::{MechanicalForceVelocity as M, RotationalTorqueAngularVelocity as R};
        use PortRole::{EffortIn as E, FlowIn as F};
        assert!(dims_match_power_port(M, E, Dims([1, 1, -2, 0, 0, 0])));
        assert!(!dims_match_power_port(M, E, Dims([1, 0, -1, 0, 0, 0])));
        assert!(dims_match_power_port(M, F, Dims([1, 0, -1, 0, 0, 0])));
        assert!(dims_match_power_port(R, E, Dims([2, 1, -2, 0, 0, 0])));
    }

    #[test]
    fn live_source_matches_static_force_bit_for_bit_and_books_work() {
        let build = |mode: u8| {
            let mut b = SystemBuilder::new();
            let m0 = Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            };
            let m = b
                .add_with_tangent_and_stage_residual("m", Mass::new(2.0).unwrap(), m0)
                .unwrap();
            let k = b
                .add_with_tangent_and_stage_residual("k", Spring::new(100.0, 0.0).unwrap(), 0.0)
                .unwrap();
            let (mp, kp, kb) = (
                b.port(m, "p").unwrap(),
                b.port(k, "a").unwrap(),
                b.port(k, "b").unwrap(),
            );
            if mode == 2 {
                b.common_flow(mp, &[kp]).unwrap();
            } else {
                let f = if mode == 1 {
                    let cell = Arc::new(AtomicU64::new(10.0_f64.to_bits()));
                    let src =
                        LiveSource::new(PortKind::MechanicalForceVelocity, PortRole::FlowIn, cell);
                    b.add_with_tangent_and_stage_residual("f", src, 0.0)
                        .unwrap()
                } else {
                    let src = Force::new(Waveform::constant(10.0).unwrap());
                    b.add_with_tangent_and_stage_residual("f", src, SourceAccount::default())
                        .unwrap()
                };
                b.common_flow(mp, &[kp, b.port(f, "p").unwrap()]).unwrap();
            }
            b.boundary(kb, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
            (b.seal().unwrap(), m)
        };
        let (mut s_live, m_live) = build(1);
        let (mut s_ref, m_ref) = build(0);
        let (s_free, _) = build(2);
        assert_eq!(
            s_free.block_modes(),
            vec![BlockSolveMode::MonolithicNewtonKrylov]
        );
        assert_eq!(s_live.block_modes(), s_ref.block_modes());
        assert_eq!(
            s_live.block_modes(),
            vec![BlockSolveMode::MonolithicNewtonKrylov]
        );
        let (pool, gate) = (
            ArenaPool::new(ArenaConfig::default()),
            CancelGate::new_clock_free(),
        );
        let key = StreamKey {
            seed: 1,
            kernel_id: 119,
            tile: 0,
            iteration: 0,
        };
        for _ in 0..20 {
            pool.scope(|a| {
                let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
                s_live.step(1.0e-3, &cx).unwrap();
                s_ref.step(1.0e-3, &cx).unwrap();
            });
            assert_eq!(
                s_live.state::<Motion>(m_live).unwrap().position_m.to_bits(),
                s_ref.state::<Motion>(m_ref).unwrap().position_m.to_bits()
            );
        }
        assert_eq!(
            s_live.ledger().energy_residual_j().to_bits(),
            s_ref.ledger().energy_residual_j().to_bits()
        );
        assert!(s_live.ledger().energy_residual_j().abs() < 1.0e-12);
    }
}
