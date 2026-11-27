//! Pre-resolved zero-allocation probe readers over a lowered [`System`].

use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use fs_fsim::dims_to_unit;
use fs_qty::Dims;
use fs_scenario::EntityId;
use fs_scenario::assembly::draft::{DraftAssembly, PhysicsDomain, ProbeDecl, ProbeKind};
use fs_system::elements::Motion;
use fs_system::lowering::{Observed, Observer};
use fs_system::{LoweredSystem, NodeId, System};

use crate::descriptor::{InputDescriptor, ProbeDescriptor, RuntimeRefusal};
use crate::runtime::{AssemblyRuntime, ProbeTarget, ResolvedProbe};

const LEN_D: [i8; 6] = [1, 0, 0, 0, 0, 0];
const VEL_D: [i8; 6] = [1, 0, -1, 0, 0, 0];
const OMG_D: [i8; 6] = [0, 0, -1, 0, 0, 0];
const TMP_D: [i8; 6] = [0, 0, 0, 1, 0, 0];
const PWR_D: [i8; 6] = [2, 1, -3, 0, 0, 0];

pub(crate) fn load_err(code: &'static str, msg: impl Into<String>) -> RuntimeRefusal {
    RuntimeRefusal::new(code, "load", msg)
}

pub(crate) fn from_lowered(
    draft: &DraftAssembly,
    lowered: LoweredSystem,
    input_descs: Vec<InputDescriptor>,
    input_cells: Vec<Arc<AtomicU64>>,
    seed: u64,
) -> Result<AssemblyRuntime, RuntimeRefusal> {
    let (mut obs, mut obs_idx) = (Vec::new(), BTreeMap::new());
    for (eid, o) in lowered.observers {
        obs_idx.insert(eid, obs.len());
        obs.push(o);
    }
    let ctx = ProbeCtx {
        draft,
        sys: &lowered.system,
        nodes: &lowered.entity_nodes,
        obs: &obs,
        obs_idx: &obs_idx,
    };
    let (mut descs, mut probes, mut offset) = (Vec::new(), Vec::new(), 0usize);
    for (i, decl) in draft.input().probes.iter().enumerate() {
        let (desc, target) = ctx.resolve_probe(i, offset, decl)?;
        offset += desc.len;
        probes.push(ResolvedProbe {
            offset: desc.offset,
            len: desc.len,
            target,
        });
        descs.push(desc);
    }
    Ok(AssemblyRuntime::build(
        draft.input().name.clone(),
        lowered.system,
        obs,
        input_descs,
        input_cells,
        descs,
        probes,
        offset,
        seed,
    ))
}

struct ProbeCtx<'a> {
    draft: &'a DraftAssembly,
    sys: &'a System,
    nodes: &'a BTreeMap<EntityId, Vec<NodeId>>,
    obs: &'a [Observer],
    obs_idx: &'a BTreeMap<EntityId, usize>,
}

fn motion_axis(name: &str, d: [i8; 6]) -> Option<usize> {
    let (qty, axis) = name.rsplit_once('-')?;
    let k = ["x", "y", "z"].iter().position(|a| *a == axis)?;
    match qty {
        "velocity" if d == VEL_D => Some(k),
        "angular-velocity" if d == OMG_D => Some(k + 3),
        _ => None,
    }
}

impl ProbeCtx<'_> {
    fn is_thermal(&self, eid: EntityId) -> bool {
        self.draft
            .input()
            .parts
            .iter()
            .any(|p| p.entity == eid && p.physics.contains(&PhysicsDomain::Thermal))
    }

    fn scalar(
        &self,
        index: usize,
        offset: usize,
        name: &str,
        m: (&str, EntityId, &str, Dims),
    ) -> ProbeDescriptor {
        let p = self.draft.input().parts.iter().find(|p| p.entity == m.1);
        let (target, unit) = (
            p.map_or_else(|| m.1.token(), |p| p.name.clone()),
            dims_to_unit(m.3).unwrap_or("SI").to_owned(),
        );
        ProbeDescriptor {
            index,
            name: name.to_owned(),
            kind: m.0.to_owned(),
            target,
            quantity: m.2.to_owned(),
            unit,
            dims: m.3.0,
            offset,
            len: 1,
            shape: Vec::new(),
        }
    }

    fn resolve_probe(
        &self,
        idx: usize,
        off: usize,
        decl: &ProbeDecl,
    ) -> Result<(ProbeDescriptor, ProbeTarget), RuntimeRefusal> {
        let miss = |msg: String| load_err("probe-unlowered", msg);
        match &decl.kind {
            ProbeKind::Thermocouple {
                time_constant: Some(_),
                ..
            } => Err(load_err(
                "probe-unsupported",
                format!("`{}` lag", decl.name),
            )),
            ProbeKind::Thermocouple {
                location,
                time_constant: None,
                ..
            } => {
                let t = self
                    .find_temp(*location)
                    .ok_or_else(|| miss(format!("`{}`", decl.name)))?;
                let m = ("thermocouple", *location, "temperature", Dims(TMP_D));
                Ok((self.scalar(idx, off, &decl.name, m), t))
            }
            ProbeKind::FieldQuantity {
                target, quantity, ..
            } => {
                let qn = quantity.name.as_str();
                let t = self
                    .resolve_field(*target, qn, quantity.dims.0)
                    .ok_or_else(|| miss(format!("`{}` `{qn}`", decl.name)))?;
                let m = ("field-quantity", *target, qn, quantity.dims);
                Ok((self.scalar(idx, off, &decl.name, m), t))
            }
            _ => Err(load_err(
                "probe-unsupported",
                format!("`{}` kind", decl.name),
            )),
        }
    }

    fn resolve_field(&self, eid: EntityId, name: &str, d: [i8; 6]) -> Option<ProbeTarget> {
        let nodes = self.nodes.get(&eid).map_or(&[][..], Vec::as_slice);
        let obs = |pred: fn(&Observed) -> bool| {
            self.obs_idx
                .get(&eid)
                .copied()
                .filter(|&i| (self.obs[i])(self.sys).as_ref().is_some_and(pred))
        };
        let motion = || {
            nodes
                .iter()
                .find(|&&n| self.sys.state::<Motion>(n).is_some())
                .copied()
        };
        match (name, d) {
            ("position" | "displacement", LEN_D) => {
                motion().map(ProbeTarget::MotionPos).or_else(|| {
                    (!self.is_thermal(eid))
                        .then(|| nodes.iter().find(|&&n| self.sys.state::<f64>(n).is_some()))
                        .flatten()
                        .map(|&n| ProbeTarget::ScalarState(n))
                })
            }
            ("velocity", VEL_D) => motion().map(ProbeTarget::MotionVel),
            ("temperature", TMP_D) => self.find_temp(eid),
            ("heat_rate", PWR_D) => {
                obs(|o| matches!(o, Observed::HeatRate(_))).map(ProbeTarget::ObserverHeatRate)
            }
            (other, d) => {
                let c = motion_axis(other, d)?;
                obs(|o| matches!(o, Observed::Body { .. })).map(|i| ProbeTarget::ObserverBody(i, c))
            }
        }
    }

    fn find_temp(&self, eid: EntityId) -> Option<ProbeTarget> {
        if let Some(&i) = self.obs_idx.get(&eid)
            && matches!((self.obs[i])(self.sys), Some(Observed::Cell(_)))
        {
            return Some(ProbeTarget::ObserverCell(i));
        }
        let nodes = self
            .is_thermal(eid)
            .then(|| self.nodes.get(&eid))
            .flatten()?;
        nodes
            .iter()
            .find(|&&n| self.sys.state::<f64>(n).is_some())
            .map(|&n| ProbeTarget::ScalarState(n))
    }
}
