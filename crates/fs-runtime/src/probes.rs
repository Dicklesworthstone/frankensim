//! Pre-resolved zero-allocation probe readers over a lowered [`System`].

use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use fs_fsim::dims_to_unit;
use fs_qty::Dims;
use fs_scenario::EntityId;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftGeometrySource, PhysicsDomain, ProbeDecl, ProbeKind,
};
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

    fn is_translational_mechanical(&self, eid: EntityId) -> bool {
        self.draft.input().parts.iter().any(|p| {
            p.entity == eid
                && p.physics.contains(&PhysicsDomain::Mechanics)
                && !p.physics.contains(&PhysicsDomain::Thermal)
                && !matches!(p.geometry, DraftGeometrySource::SpiralSpring(_))
        })
    }

    fn unique_motion(&self, eid: EntityId, decl: &str) -> Result<Option<NodeId>, RuntimeRefusal> {
        let mut found = None;
        for &n in self.nodes.get(&eid).map_or(&[][..], Vec::as_slice) {
            if self.sys.state::<Motion>(n).is_some() {
                if found.is_some() {
                    return Err(load_err(
                        "probe-ambiguous",
                        format!("`{decl}` multiple Motion nodes"),
                    ));
                }
                found = Some(n);
            }
        }
        Ok(found)
    }

    fn unique_scalar(&self, eid: EntityId, decl: &str) -> Result<Option<NodeId>, RuntimeRefusal> {
        let mut found = None;
        for &n in self.nodes.get(&eid).map_or(&[][..], Vec::as_slice) {
            if self.sys.state::<f64>(n).is_some() {
                if found.is_some() {
                    return Err(load_err(
                        "probe-ambiguous",
                        format!("`{decl}` multiple scalar nodes"),
                    ));
                }
                found = Some(n);
            }
        }
        Ok(found)
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
                    .find_temp(*location, &decl.name)?
                    .ok_or_else(|| miss(format!("`{}`", decl.name)))?;
                let m = ("thermocouple", *location, "temperature", Dims(TMP_D));
                Ok((self.scalar(idx, off, &decl.name, m), t))
            }
            ProbeKind::FieldQuantity {
                target, quantity, ..
            } => {
                let qn = quantity.name.as_str();
                let t = self
                    .resolve_field(*target, &decl.name, qn, quantity.dims.0)?
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

    fn resolve_field(
        &self,
        eid: EntityId,
        decl: &str,
        name: &str,
        d: [i8; 6],
    ) -> Result<Option<ProbeTarget>, RuntimeRefusal> {
        let obs = |pred: fn(&Observed) -> bool| {
            self.obs_idx
                .get(&eid)
                .copied()
                .filter(|&i| (self.obs[i])(self.sys).as_ref().is_some_and(pred))
        };
        match (name, d) {
            ("position" | "displacement", LEN_D) => {
                if let Some(n) = self.unique_motion(eid, decl)? {
                    Ok(Some(ProbeTarget::MotionPos(n)))
                } else if self.is_translational_mechanical(eid) {
                    Ok(self.unique_scalar(eid, decl)?.map(ProbeTarget::ScalarState))
                } else {
                    Ok(None)
                }
            }
            ("velocity", VEL_D) => Ok(self.unique_motion(eid, decl)?.map(ProbeTarget::MotionVel)),
            ("temperature", TMP_D) => self.find_temp(eid, decl),
            ("heat_rate", PWR_D) => {
                Ok(obs(|o| matches!(o, Observed::HeatRate(_))).map(ProbeTarget::ObserverHeatRate))
            }
            (other, d) => Ok(motion_axis(other, d)
                .and_then(|c| obs(|o| matches!(o, Observed::Body { .. })).map(|i| (i, c)))
                .map(|(i, c)| ProbeTarget::ObserverBody(i, c))),
        }
    }

    fn find_temp(&self, eid: EntityId, decl: &str) -> Result<Option<ProbeTarget>, RuntimeRefusal> {
        if let Some(&i) = self.obs_idx.get(&eid)
            && matches!((self.obs[i])(self.sys), Some(Observed::Cell(_)))
        {
            return Ok(Some(ProbeTarget::ObserverCell(i)));
        }
        if !self.is_thermal(eid) {
            return Ok(None);
        }
        Ok(self.unique_scalar(eid, decl)?.map(ProbeTarget::ScalarState))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use fs_fsim::{PackCatalog, admit_assembly_source};
    use fs_ports::PortKind;
    use fs_qty::QtyAny;
    use fs_scenario::Scenario;
    use fs_scenario::assembly::draft::SpiralSpringGenerator;
    use fs_system::SystemBuilder;
    use fs_system::elements::{Mass, Storage};

    use super::*;

    fn base_draft() -> (DraftAssembly, EntityId) {
        let fsim = include_bytes!("../../../examples/assemblies/forced_oscillator.fsim");
        let pack = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let adm =
            admit_assembly_source(fsim, &PackCatalog::from_raw_bytes(&[pack]).unwrap()).unwrap();
        let eid = adm.draft.input().parts[0].entity;
        (adm.draft, eid)
    }

    #[test]
    fn refuses_non_mechanical_and_spiral_scalar_and_multi_motion() {
        let (draft, eid) = base_draft();
        let (mut b, m0) = (
            SystemBuilder::new(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        );
        let s_node = b
            .add(
                "s",
                Storage::in_domain(PortKind::ElectricalVoltageCurrent, 1.0).unwrap(),
                1.0,
            )
            .unwrap();
        let m1 = b.add("m1", Mass::new(1.0).unwrap(), m0).unwrap();
        let m2 = b.add("m2", Mass::new(1.0).unwrap(), m0).unwrap();
        for n in [s_node, m1, m2] {
            let p = b.port(n, "p").unwrap();
            b.boundary(p, Box::new(|_, _, o| o[0] = 0.0)).unwrap();
        }
        let sys = b.seal().unwrap();
        let empty = BTreeMap::new();

        // 1. Multi-Motion part refuses `position` and `velocity` as `probe-ambiguous`.
        let multi = BTreeMap::from([(eid, vec![m1, m2])]);
        let ctx_m = ProbeCtx {
            draft: &draft,
            sys: &sys,
            nodes: &multi,
            obs: &[],
            obs_idx: &empty,
        };
        assert_eq!(
            ctx_m
                .resolve_field(eid, "p", "position", LEN_D)
                .err()
                .unwrap()
                .code,
            "probe-ambiguous"
        );
        assert_eq!(
            ctx_m
                .resolve_field(eid, "v", "velocity", VEL_D)
                .err()
                .unwrap()
                .code,
            "probe-ambiguous"
        );

        // 2. Non-mechanical, non-thermal part with scalar f64 refuses `position`.
        let mut inp_el = draft.input().clone();
        inp_el.parts[0].physics = BTreeSet::from([PhysicsDomain::Electromagnetics]);
        let scen = Scenario::new(&inp_el.name, 1, inp_el.environment.base.clone());
        let d_el = DraftAssembly::admit(inp_el, &scen).unwrap();
        let single = BTreeMap::from([(eid, vec![s_node])]);
        let ctx_el = ProbeCtx {
            draft: &d_el,
            sys: &sys,
            nodes: &single,
            obs: &[],
            obs_idx: &empty,
        };
        assert!(
            ctx_el
                .resolve_field(eid, "p", "position", LEN_D)
                .unwrap()
                .is_none()
        );

        // 3. Rotational SpiralSpring part with scalar f64 refuses `position`.
        let mut inp_sp = draft.input().clone();
        inp_sp.parts[0].geometry = DraftGeometrySource::SpiralSpring(SpiralSpringGenerator {
            inner_radius: QtyAny::new(0.01, Dims(LEN_D)),
            outer_radius: QtyAny::new(0.05, Dims(LEN_D)),
            turns: 4.0,
            blade_thickness: QtyAny::new(0.001, Dims(LEN_D)),
            blade_height: QtyAny::new(0.005, Dims(LEN_D)),
        });
        let d_sp = DraftAssembly::admit(inp_sp, &scen).unwrap();
        let ctx_sp = ProbeCtx {
            draft: &d_sp,
            sys: &sys,
            nodes: &single,
            obs: &[],
            obs_idx: &empty,
        };
        assert!(
            ctx_sp
                .resolve_field(eid, "p", "position", LEN_D)
                .unwrap()
                .is_none()
        );
    }
}
