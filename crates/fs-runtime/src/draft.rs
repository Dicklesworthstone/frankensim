//! Lowering `.fsim` sources and admitted [`DraftAssembly`] graphs with
//! lock-free interactive inputs into [`AssemblyRuntime`].

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use fs_fsim::{PackCatalog, admit_assembly_source, dims_to_unit};
use fs_matdb::Provenance;
use fs_ports::PortKind;
use fs_qty::Dims;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftInterfaceDecl, DraftInterfaceKind, DraftPartKind,
    SourceTarget, TimeHistory,
};
use fs_scenario::{EntityDeclaration, InterfacePair, Scenario};
use fs_system::lowering::PowerPortSchemaCard;
use fs_system::{CardStore, LoweringRegistry};

use crate::descriptor::{InputDescriptor, RuntimeRefusal};
use crate::live_source::{
    LiveInputSpec, RUNTIME_INPUT_FAMILY, RegisteredInput, register_runtime_inputs,
};
use crate::probes::{from_lowered, load_err};
use crate::runtime::AssemblyRuntime;

impl AssemblyRuntime {
    /// Parse and admit a `.fsim` source with canonical `fs-matdb` card packs
    /// (`FSMCDPK\0` / `FSINTPK\0`), then lower into a runtime.
    pub fn load(
        fsim_source: &[u8],
        packs: &[&[u8]],
        seed_override: Option<u64>,
    ) -> Result<Self, RuntimeRefusal> {
        let cat = PackCatalog::from_raw_bytes(packs).map_err(|r| load_err(r.code, r.what))?;
        let adm = admit_assembly_source(fsim_source, &cat).map_err(|r| load_err(r.code, r.what))?;
        let seed = seed_override.unwrap_or(adm.seed);
        Self::from_draft(&adm.draft, &adm.card_store, seed)
    }

    /// Lower an admitted [`DraftAssembly`], promoting constant power/signal
    /// `SourceDecl` entries on non-driven parts into live interactive inputs.
    pub fn from_draft(
        draft: &DraftAssembly,
        cards: &CardStore,
        seed: u64,
    ) -> Result<Self, RuntimeRefusal> {
        Self::from_draft_inputs(draft, cards, &[], seed)
    }

    /// Lower an admitted [`DraftAssembly`], promoting constant power/signal
    /// `SourceDecl` entries on non-driven parts alongside explicit `live_inputs`.
    pub fn from_draft_inputs(
        draft: &DraftAssembly,
        cards: &CardStore,
        live_inputs: &[LiveInputSpec],
        seed: u64,
    ) -> Result<Self, RuntimeRefusal> {
        let mut reg = LoweringRegistry::with_kernel_defaults();
        let (draft, cards, descs, cells) =
            attach_live_inputs(draft, cards, live_inputs, &mut reg, seed)?;
        let lowered = reg
            .lower(&draft, &cards)
            .map_err(|e| load_err("lowering-refused", e.to_string()))?;
        from_lowered(&draft, lowered, descs, cells, seed)
    }
}

fn infer_port_kind(dims: Dims) -> Option<PortKind> {
    use PortKind as P;
    match dims.0 {
        [1, 1, -2, 0, 0, 0] | [1, 0, -1, 0, 0, 0] | [0, 0, 0, 0, 0, 0] => {
            Some(P::MechanicalForceVelocity)
        }
        [2, 1, -2, 0, 0, 0] | [0, 0, -1, 0, 0, 0] => Some(P::RotationalTorqueAngularVelocity),
        [2, 1, -3, 0, -1, 0] | [0, 0, 0, 0, 1, 0] => Some(P::ElectricalVoltageCurrent),
        [-1, 1, -2, 0, 0, 0] | [3, 0, -1, 0, 0, 0] => Some(P::FluidPressureFlux),
        _ => None,
    }
}

fn collect_live_specs(inp: &mut DraftAssemblyInput, extra: &[LiveInputSpec]) -> Vec<LiveInputSpec> {
    let (mut specs, mut kept) = (Vec::new(), Vec::with_capacity(inp.sources.len()));
    for src in inp.sources.drain(..) {
        let promoted = match (&src.target, &src.history) {
            (SourceTarget::Port { entity, port }, TimeHistory::Constant(q)) => {
                let part = inp.parts.iter().find(|p| p.entity == *entity);
                part.is_some_and(|p| p.kind != DraftPartKind::DrivenBody)
                    .then(|| infer_port_kind(q.dims))
                    .flatten()
                    .map(|domain| LiveInputSpec {
                        name: src.name.clone(),
                        target: *entity,
                        target_name: part.map_or_else(|| entity.token(), |p| p.name.clone()),
                        port: port.clone(),
                        domain,
                        dims: q.dims,
                        default_value: q.value,
                    })
            }
            _ => None,
        };
        if let Some(s) = promoted {
            specs.push(s);
        } else {
            kept.push(src);
        }
    }
    inp.sources = kept;
    specs.extend_from_slice(extra);
    specs
}

fn declare_input_interface(
    inp: &mut DraftAssemblyInput,
    cards: &mut CardStore,
    spec: &LiveInputSpec,
) -> Result<String, RuntimeRefusal> {
    let bad = |e: &dyn core::fmt::Display| load_err("invalid-input", e.to_string());
    let card = PowerPortSchemaCard {
        family: RUNTIME_INPUT_FAMILY.to_owned(),
        schema_name: spec.name.clone(),
        domain: spec.domain,
        description: spec.port.clone(),
        provenance: Provenance {
            source: "fs-runtime interactive boundary".to_owned(),
            license: "workspace-internal".to_owned(),
            artifact: None,
        },
    };
    let schema = cards.insert_power_port_schema(card).map_err(|e| bad(&e))?;
    let dname = format!("runtime_input:{}", spec.name);
    let mut decl = |d| inp.catalog.declare(d).map_err(|e| bad(&e));
    let sa = decl(EntityDeclaration::surface(
        spec.target,
        &format!("{dname}:ref"),
    ))?;
    let sb = decl(EntityDeclaration::surface(
        spec.target,
        &format!("{dname}:port"),
    ))?;
    let eid = decl(EntityDeclaration::interface(
        inp.root,
        &dname,
        InterfacePair::ordered(sa, sb),
    ))?;
    inp.interfaces.push(DraftInterfaceDecl {
        name: dname.clone(),
        entity: eid,
        feature_a: sa,
        feature_b: sb,
        kind: DraftInterfaceKind::PowerPort { schema },
    });
    Ok(dname)
}

type Attached = (
    DraftAssembly,
    CardStore,
    Vec<InputDescriptor>,
    Vec<Arc<AtomicU64>>,
);

fn attach_live_inputs(
    draft: &DraftAssembly,
    cards: &CardStore,
    extra: &[LiveInputSpec],
    reg: &mut LoweringRegistry,
    seed: u64,
) -> Result<Attached, RuntimeRefusal> {
    let mut input = draft.input().clone();
    let specs = collect_live_specs(&mut input, extra);
    if specs.is_empty() {
        return Ok((draft.clone(), cards.clone(), Vec::new(), Vec::new()));
    }
    let (mut cards_out, mut seen, mut by_name) = (cards.clone(), BTreeSet::new(), BTreeMap::new());
    let (mut descs, mut cells) = (Vec::new(), Vec::new());
    for (index, spec) in specs.iter().enumerate() {
        if spec.name.trim().is_empty()
            || !seen.insert(spec.name.clone())
            || !spec.default_value.is_finite()
        {
            return Err(load_err("invalid-input", format!("bad `{}`", spec.name)));
        }
        let dname = declare_input_interface(&mut input, &mut cards_out, spec)?;
        let cell = Arc::new(AtomicU64::new(spec.default_value.to_bits()));
        descs.push(InputDescriptor {
            index,
            name: spec.name.clone(),
            target: spec.target_name.clone(),
            port: spec.port.clone(),
            unit: dims_to_unit(spec.dims).unwrap_or("SI").to_owned(),
            dims: spec.dims.0,
            width: 1,
            default_value: spec.default_value,
        });
        by_name.insert(
            dname,
            RegisteredInput {
                port: spec.port.clone(),
                dims: spec.dims,
                default_value: spec.default_value,
                cell: Arc::clone(&cell),
            },
        );
        cells.push(cell);
    }
    let scen = Scenario::new(&input.name, seed, input.environment.base.clone());
    let adm =
        DraftAssembly::admit(input, &scen).map_err(|e| load_err("admit-refused", e.to_string()))?;
    register_runtime_inputs(reg, Arc::new(by_name));
    Ok((adm, cards_out, descs, cells))
}
