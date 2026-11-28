//! Core zero-allocation runtime over a sealed [`System`].

use std::collections::BTreeSet;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::acoustic::AcousticMicrophone;
use fs_system::elements::Motion;
use fs_system::lowering::{Observed, Observer};
use fs_system::{BlockSolveMode, Ledger, NodeId, StepReport, System};

use crate::descriptor::{
    InputDescriptor, ProbeDescriptor, RuntimeRefusal, SNAPSHOT_DESCRIPTOR_SCHEMA,
    SnapshotDescriptor, inputs_to_json,
};

/// Pre-bound interactive input channel for [`AssemblyRuntime::from_system`].
pub struct InputBinding {
    /// Descriptor exported by `inputs()`.
    pub descriptor: InputDescriptor,
    /// Lock-free bit cell read by the kernel source or boundary closure.
    pub cell: Arc<AtomicU64>,
}

/// Reader closure writing one probe's slice directly into the snapshot buffer.
pub type CustomProbeFn = Box<dyn Fn(&System, &mut [f64])>;

/// Pre-bound probe slot for [`AssemblyRuntime::from_system`].
pub struct ProbeBinding {
    /// Descriptor exported by `probes()`.
    pub descriptor: ProbeDescriptor,
    /// Zero-allocation reader populating `out[..descriptor.len]`.
    pub reader: CustomProbeFn,
}

struct InputSlot {
    pending: AtomicU64,
    live: Arc<AtomicU64>,
}

pub(crate) enum ProbeTarget {
    MotionPos(NodeId),
    MotionVel(NodeId),
    ScalarState(NodeId),
    ObserverCell(usize),
    ObserverHeatRate(usize),
    ObserverBody(usize, usize),
    AcousticMic(usize),
    Custom(CustomProbeFn),
}

pub(crate) struct ResolvedProbe {
    pub(crate) offset: usize,
    pub(crate) len: usize,
    pub(crate) target: ProbeTarget,
}

fn err(code: &'static str, stage: &'static str, msg: impl Into<String>) -> RuntimeRefusal {
    RuntimeRefusal::new(code, stage, msg)
}

fn next_offset(d: &ProbeDescriptor, off: usize) -> Option<usize> {
    let sh = &d.shape;
    let ok = d.len > 0
        && (sh.is_empty() || sh.iter().try_fold(1usize, |a, &x| a.checked_mul(x)) == Some(d.len));
    ok.then(|| off.checked_add(d.len)).flatten()
}

/// Generic deterministic runtime over a sealed [`System`].
pub struct AssemblyRuntime {
    system: System,
    observers: Vec<Observer>,
    microphones: Vec<AcousticMicrophone>,
    pool: ArenaPool,
    gate: CancelGate,
    stream_key: StreamKey,
    input_descs: Vec<InputDescriptor>,
    input_slots: Vec<InputSlot>,
    descriptor: SnapshotDescriptor,
    probes: Vec<ResolvedProbe>,
    halted: bool,
    last_refusal: Option<RuntimeRefusal>,
}

impl AssemblyRuntime {
    /// Wrap an already-sealed [`System`] with explicit inputs and probe bindings.
    pub fn from_system(
        assembly: impl Into<String>,
        system: System,
        inputs: Vec<InputBinding>,
        probes: Vec<ProbeBinding>,
        seed: u64,
    ) -> Result<Self, RuntimeRefusal> {
        let (mut input_descs, mut input_cells, mut seen) =
            (Vec::new(), Vec::new(), BTreeSet::new());
        for (i, mut b) in inputs.into_iter().enumerate() {
            let d = &mut b.descriptor;
            let ok = !d.name.trim().is_empty()
                && seen.insert(d.name.clone())
                && d.default_value.is_finite()
                && d.width > 0;
            if !ok {
                return Err(err("invalid-input", "load", format!("bad `{}`", d.name)));
            }
            d.index = i;
            b.cell.store(d.default_value.to_bits(), Ordering::Relaxed);
            input_descs.push(b.descriptor);
            input_cells.push(b.cell);
        }
        let (mut probe_descs, mut resolved, mut off) = (Vec::new(), Vec::new(), 0usize);
        seen.clear();
        for (i, mut b) in probes.into_iter().enumerate() {
            let d = &mut b.descriptor;
            let ok = !d.name.trim().is_empty() && seen.insert(d.name.clone());
            let Some(next) = ok.then(|| next_offset(d, off)).flatten() else {
                return Err(err("invalid-probe", "load", format!("bad `{}`", d.name)));
            };
            (d.index, d.offset) = (i, off);
            let len = d.len;
            probe_descs.push(b.descriptor);
            resolved.push(ResolvedProbe {
                offset: off,
                len,
                target: ProbeTarget::Custom(b.reader),
            });
            off = next;
        }
        Ok(Self::build(
            assembly.into(),
            system,
            Vec::new(),
            Vec::new(),
            input_descs,
            input_cells,
            probe_descs,
            resolved,
            off,
            seed,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        assembly: String,
        system: System,
        observers: Vec<Observer>,
        microphones: Vec<AcousticMicrophone>,
        input_descs: Vec<InputDescriptor>,
        input_cells: Vec<Arc<AtomicU64>>,
        probe_descs: Vec<ProbeDescriptor>,
        probes: Vec<ResolvedProbe>,
        total_len: usize,
        seed: u64,
    ) -> Self {
        let input_slots = input_descs
            .iter()
            .zip(input_cells)
            .map(|(d, live)| InputSlot {
                pending: AtomicU64::new(d.default_value.to_bits()),
                live,
            })
            .collect();
        Self {
            system,
            observers,
            microphones,
            pool: ArenaPool::new(ArenaConfig::default()),
            gate: CancelGate::new_clock_free(),
            stream_key: StreamKey {
                seed,
                kernel_id: 119,
                tile: 0,
                iteration: 0,
            },
            input_descs,
            input_slots,
            descriptor: SnapshotDescriptor {
                schema: SNAPSHOT_DESCRIPTOR_SCHEMA.to_owned(),
                assembly,
                total_len,
                probes: probe_descs,
            },
            probes,
            halted: false,
            last_refusal: None,
        }
    }

    /// Declared interactive input descriptors.
    #[must_use]
    pub fn inputs(&self) -> &[InputDescriptor] {
        &self.input_descs
    }

    /// Canonical JSON array of interactive input descriptors.
    #[must_use]
    pub fn inputs_json(&self) -> String {
        inputs_to_json(&self.input_descs)
    }

    /// Stage one interactive input value without allocating; latched at the next [`Self::step`].
    pub fn set_input(&self, index: usize, value: f64) -> Result<(), RuntimeRefusal> {
        let slot = self
            .input_slots
            .get(index)
            .ok_or_else(|| err("unknown-input", "input", format!("idx {index}")))?;
        if !value.is_finite() {
            return Err(err("non-finite-input", "input", format!("idx {index}")));
        }
        slot.pending.store(value.to_bits(), Ordering::Relaxed);
        Ok(())
    }

    /// Stage one interactive input by declared name without allocating; latched at the next [`Self::step`].
    pub fn set_input_by_name(&self, name: &str, value: f64) -> Result<(), RuntimeRefusal> {
        let idx = self
            .input_descs
            .iter()
            .position(|d| d.name == name)
            .ok_or_else(|| err("unknown-input", "input", format!("`{name}`")))?;
        self.set_input(idx, value)
    }

    /// Assembly-generated snapshot descriptor.
    #[must_use]
    pub const fn probes(&self) -> &SnapshotDescriptor {
        &self.descriptor
    }

    /// Canonical JSON encoding of [`Self::probes`].
    #[must_use]
    pub fn probes_json(&self) -> String {
        self.descriptor.to_json()
    }

    /// Required length of the flat `f64` buffer passed to [`Self::read_probes`].
    #[must_use]
    pub const fn snapshot_len(&self) -> usize {
        self.descriptor.total_len
    }

    /// Latch staged inputs and advance the kernel by `dt_s` seconds in deterministic mode.
    pub fn step(&mut self, dt_s: f64) -> Result<StepReport, RuntimeRefusal> {
        if self.halted {
            return Err(err("runtime-halted", "step", "runtime halted"));
        }
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(self.halt("invalid-dt", format!("dt={dt_s}")));
        }
        for mic in &self.microphones {
            let want = mic.mechanics_dt_s();
            if (dt_s - want).abs() > 1.0e-9 * want {
                return Err(self.halt("microphone-dt-mismatch", format!("dt={dt_s} != {want}")));
            }
        }
        for slot in &self.input_slots {
            slot.live
                .store(slot.pending.load(Ordering::Relaxed), Ordering::Relaxed);
        }
        let res = self.pool.scope(|arena| {
            let cx = Cx::new(
                &self.gate,
                arena,
                self.stream_key,
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            self.system.step(dt_s, &cx)
        });
        let rep = res.map_err(|e| self.halt("kernel-step-refused", e.to_string()))?;
        for mic in &mut self.microphones {
            if let Err(why) = mic.advance(&self.system) {
                return Err(self.halt("microphone-step-refused", why.to_owned()));
            }
        }
        Ok(rep)
    }

    /// Fill `buffer[..snapshot_len()]` with current probe values without allocating.
    pub fn read_probes(&self, buffer: &mut [f64]) -> Result<usize, RuntimeRefusal> {
        let need = self.descriptor.total_len;
        if buffer.len() < need {
            return Err(err("short-buffer", "probe", format!("need {need}")));
        }
        for p in &self.probes {
            self.eval_probe(p, &mut buffer[p.offset..p.offset + p.len]);
        }
        Ok(need)
    }

    /// AudioWorklet-friendly block render into `out` (`f32`) with zero allocations.
    ///
    /// On any refusal, halts cleanly, fills remaining frames with `0.0`, and queues
    /// the refusal for [`Self::take_refusal`]. Returns rendered frame count.
    pub fn render_f32(
        &mut self,
        probe_index: usize,
        dt_s: f64,
        substeps: usize,
        out: &mut [f32],
    ) -> usize {
        let bad_dt = !(dt_s.is_finite() && dt_s > 0.0);
        let bad_pr = self.probes.get(probe_index).is_none_or(|p| p.len != 1);
        if self.halted || bad_dt || substeps == 0 || bad_pr {
            if !self.halted {
                self.halt(
                    if bad_dt {
                        "invalid-dt"
                    } else {
                        "invalid-render"
                    },
                    format!("dt={dt_s} sub={substeps} p={probe_index}"),
                );
            }
            out.fill(0.0);
            return 0;
        }
        let mut slot = [0.0_f64; 1];
        for i in 0..out.len() {
            for _ in 0..substeps {
                if self.step(dt_s).is_err() {
                    out[i..].fill(0.0);
                    return i;
                }
            }
            self.eval_probe(&self.probes[probe_index], &mut slot);
            let sample = slot[0] as f32;
            if !sample.is_finite() {
                self.halt("non-finite-probe", format!("probe {probe_index}"));
                out[i..].fill(0.0);
                return i;
            }
            out[i] = sample;
        }
        out.len()
    }

    /// Borrow the kernel's energy/conservation ledger.
    #[must_use]
    pub fn ledger(&self) -> &Ledger {
        self.system.ledger()
    }

    /// Canonical JSON export of the kernel window ledger.
    #[must_use]
    pub fn ledger_json(&self) -> String {
        self.system.ledger().export().canonical_json
    }

    /// Poll and clear the refusal that halted the runtime, reported at most once.
    pub fn take_refusal(&mut self) -> Option<RuntimeRefusal> {
        self.last_refusal.take()
    }

    /// Whether the runtime has halted after a refusal.
    #[must_use]
    pub const fn is_halted(&self) -> bool {
        self.halted
    }

    /// Current accepted simulation time [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.system.time_s()
    }

    /// Total accepted kernel steps.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.system.steps()
    }

    /// Per-block solve modes selected by the sealed kernel.
    #[must_use]
    pub fn block_modes(&self) -> Vec<BlockSolveMode> {
        self.system.block_modes()
    }

    fn halt(&mut self, code: &'static str, msg: String) -> RuntimeRefusal {
        let r = err(code, "step", msg);
        self.halted = true;
        self.last_refusal = Some(r.clone());
        r
    }

    fn eval_probe(&self, p: &ResolvedProbe, dst: &mut [f64]) {
        match &p.target {
            ProbeTarget::MotionPos(n) => {
                dst[0] = self
                    .system
                    .state::<Motion>(*n)
                    .map_or(f64::NAN, |m| m.position_m);
            }
            ProbeTarget::MotionVel(n) => {
                dst[0] = self
                    .system
                    .state::<Motion>(*n)
                    .map_or(f64::NAN, |m| m.velocity_m_s);
            }
            ProbeTarget::ScalarState(n) => {
                dst[0] = self.system.state::<f64>(*n).copied().unwrap_or(f64::NAN);
            }
            ProbeTarget::ObserverCell(i) => {
                dst[0] = match (self.observers[*i])(&self.system) {
                    Some(Observed::Cell(t)) => t,
                    _ => f64::NAN,
                };
            }
            ProbeTarget::ObserverHeatRate(i) => {
                dst[0] = match (self.observers[*i])(&self.system) {
                    Some(Observed::HeatRate(q)) => q,
                    _ => f64::NAN,
                };
            }
            ProbeTarget::ObserverBody(i, c) => {
                dst[0] = match (self.observers[*i])(&self.system) {
                    Some(Observed::Body {
                        velocity_m_s: v,
                        angular_velocity_rad_s: w,
                    }) => [v[0], v[1], v[2], w[0], w[1], w[2]][*c],
                    _ => f64::NAN,
                };
            }
            ProbeTarget::AcousticMic(i) => dst[0] = self.microphones[*i].last_pressure_pa(),
            ProbeTarget::Custom(f) => f(&self.system, dst),
        }
    }
}
