//! Data-only multi-part bodies on the original material shell and mortar owners.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::shell::SharedMaterialShellAssemblyBody;
use fs_couple::shear_shell::assembly::{
    MaterialShellAssembly, MaterialShellAssemblyImage, ShellAssemblyForcePort, ShellAssemblyReductionOptions,
};
use fs_plate::shell::junction::{
    EdgeTie, FaceTie, RotationCoupling, TieBudget, TieKind, TieOffset, TieSide,
};

/// A shell-assembly table written as this v4 profile.
#[path = "assembly_table.rs"]
pub(crate) mod table;
#[path = "assembly/glue.rs"]
mod glue;

// The il Cannone violin body at 6 mm
// (data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h6.tsv) has 64 components and
// 157 ties, the largest table on file; both caps are exactly what it needs.
const MAX_COMPONENTS: usize = 64;
// Totals across every component of one assembled body: four components at the per-component ceiling.
const MAX_TOTAL_NODES: usize = 4 * MAX_NODES;
const MAX_TOTAL_TRIANGLES: usize = 4 * MAX_TRIANGLES;
const MAX_TIES: usize = 157;

#[derive(Clone, Debug)]
struct Component {
    name: String,
    shell: super::Spec,
}
#[derive(Clone, Debug)]
struct AssemblyBridge {
    component: usize,
    footprint: Bridge,
}

#[derive(Clone, Debug)]
pub(crate) struct Spec {
    components: Vec<Component>,
    edges: Vec<EdgeTie>,
    faces: Vec<FaceTie>,
    glue: glue::Inputs,
    ties: TieBudget,
    bridges: Vec<AssemblyBridge>,
    window: (f64, f64),
    projection_work: usize,
    coating_memory: Option<coating::MemoryInput>,
}

impl Spec {
    pub(crate) fn parse(bytes: &[u8]) -> Result<Self> {
        let mut r = Records {
            lines: std::str::from_utf8(bytes)?.lines(),
            line: 0,
        };
        r.row("frankensim-string-body-v4", 0)?;
        let (window, projection_work) = read_window(&mut r)?;
        let coating_memory = coating::MemoryInput::read(&mut r)?;
        let n = count(r.row("components", 1)?[0], 1, MAX_COMPONENTS)?;
        let mut components = Vec::<Component>::with_capacity(n);
        let (mut nodes, mut triangles) = (0, 0);
        for _ in 0..n {
            let name = r.row("component", 1)?[0];
            if name.len() > 128
                || name.chars().any(char::is_control)
                || components.iter().any(|c| c.name == name)
            {
                return Err("assembly component names must be unique bounded tokens".into());
            }
            let mut shell = super::Spec::read(&mut r, Some((window, projection_work)))?;
            shell.coating_memory = coating_memory.clone();
            nodes += shell.nodes.len();
            triangles += shell.triangles.len();
            if nodes > MAX_TOTAL_NODES || triangles > MAX_TOTAL_TRIANGLES {
                return Err(format!(
                    "assembled body exceeds {MAX_TOTAL_NODES} total nodes or {MAX_TOTAL_TRIANGLES} total triangles"
                )
                .into());
            }
            r.row("end-component", 0)?;
            components.push(Component {
                name: name.into(),
                shell,
            });
        }
        if !components
            .iter()
            .any(|c| c.shell.regions.iter().any(|r| r.pack == "@body"))
        {
            return Err("assembled body must retain the performance body material pack".into());
        }
        if coating_memory.is_some() && components.iter().all(|c| c.shell.coatings.is_empty()) {
            return Err("coating-memory requires at least one physical film in the assembly".into());
        }
        let v = r.row("tie-budget", 4)?;
        let ties = TieBudget {
            max_segments: count(v[0], 1, 65536)?,
            max_slave_dofs: count(v[1], 1, 6 * MAX_TOTAL_NODES)?,
            gap_tolerance_m: scalar(v[2])?,
            angle_tolerance_rad: scalar(v[3])?,
        };
        if ties.gap_tolerance_m < 0.0
            || !(0.0..=core::f64::consts::PI).contains(&ties.angle_tolerance_rad)
        {
            return Err("tie tolerances need nonnegative gap and angle in [0,pi]".into());
        }
        let mut glue = glue::Inputs::default();
        let n = count(r.row("edge-ties", 1)?[0], 0, MAX_TIES)?;
        let mut edges = Vec::with_capacity(n);
        for _ in 0..n {
            let v = r.row("edge-tie", 6)?;
            let a = component(&components, v[0])?;
            let b = component(&components, v[1])?;
            if a == b {
                return Err("assembly ties must address distinct components".into());
            }
            let offset = offset(v[2])?;
            let rotations = rotations(v[3])?;
            let na = count(v[4], 2, components[a].shell.nodes.len())?;
            let nb = count(v[5], 2, components[b].shell.nodes.len())?;
            let master = indices(&mut r, "master-node", na, components[a].shell.nodes.len())?;
            let slave = indices(&mut r, "slave-node", nb, components[b].shell.nodes.len())?;
            edges.push(EdgeTie {
                a: TieSide {
                    component: a,
                    nodes: master,
                },
                b: TieSide {
                    component: b,
                    nodes: slave,
                },
                offset,
                rotations,
                kind: TieKind::Rigid,
            });
            glue.read_edge(&mut r, &components, &edges, edges.len() - 1)?;
        }
        let n = count(r.row("face-ties", 1)?[0], 0, MAX_TIES - edges.len())?;
        let mut faces = Vec::with_capacity(n);
        for _ in 0..n {
            let v = r.row("face-tie", 6)?;
            let a = component(&components, v[0])?;
            let b = component(&components, v[1])?;
            if a == b {
                return Err("assembly ties must address distinct components".into());
            }
            let offset = offset(v[2])?;
            let rotations = rotations(v[3])?;
            let na = count(v[4], 1, components[a].shell.triangles.len())?;
            let nb = count(v[5], 1, components[b].shell.triangles.len())?;
            let master_tris = indices(
                &mut r,
                "master-triangle",
                na,
                components[a].shell.triangles.len(),
            )?;
            let slave_tris = indices(
                &mut r,
                "slave-triangle",
                nb,
                components[b].shell.triangles.len(),
            )?;
            faces.push(FaceTie {
                master_component: a,
                master_tris,
                slave_component: b,
                slave_tris,
                offset,
                rotations,
                kind: TieKind::Rigid,
            });
            glue.read_face(&mut r, &components, &faces, faces.len() - 1)?;
        }
        let n = count(r.row("bridges", 1)?[0], 1, 8)?;
        let mut bridges = Vec::with_capacity(n);
        for _ in 0..n {
            let (record, fields) = Bridge::record(&r);
            let v = r.row(record, fields + 1)?;
            let component = component(&components, v[0])?;
            let shell = &components[component].shell;
            bridges.push(AssemblyBridge {
                component,
                footprint: Bridge::read_fields(&v[1..], shell.nodes.len(), &shell.triangles, &shell.supports)?,
            });
        }
        if r.lines.next().is_some() {
            return Err("unconsumed assembly body profile record".into());
        }
        Ok(Self {
            components,
            edges,
            faces,
            glue,
            ties,
            bridges,
            window,
            projection_work,
            coating_memory,
        })
    }

    pub(crate) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if self.bridges.len() != spec.strings.len() {
            return Err("assembled body requires one physical bridge per actual string".into());
        }
        if !(1..=MAX_MODES).contains(&spec.body.modes) {
            return Err(format!("assembled body admits 1..={MAX_MODES} retained global body modes").into());
        }
        spec.network_dimensions_with_body_connections(true, self.additional_body_connections())?;
        if let Some(memory) = &self.coating_memory { memory.preflight(spec)?; }
        for component in &self.components {
            component.shell.geometry(spec)?;
            component
                .shell
                .interfaces
                .preflight(&component.shell, spec)?;
        }
        for bridge in &self.bridges {
            bridge.footprint.position(&self.components[bridge.component].shell, spec)?;
        }
        Ok(())
    }

    pub(crate) fn additional_body_connections(&self) -> usize {
        self.components
            .iter()
            .map(|c| c.shell.additional_body_connections())
            .sum()
    }

    pub(crate) fn preflight_rotation(&self, spec: &input::Spec, geometry: &ShellRotationGeometry) -> Result<()> {
        let component = self.components.get(geometry.shell_component)
            .ok_or("shell rotation names an absent original assembly component")?;
        component.shell.preflight_rotation(spec,geometry.triangle,geometry.barycentric,geometry.axis)
    }

    pub(crate) fn material_memory_states(&self) -> usize {
        self.coating_memory.as_ref().map_or(0, |memory| memory.options.max_memory_coordinates)
    }

    pub(crate) fn bind_material_memory(&self, bank: ReciprocalStringBank, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        Ok(if self.coating_memory.is_some() {
            bank.with_material_shell_coating_relaxation(None, gate)?
        } else { bank })
    }

    pub(crate) fn load_interfaces(
        &mut self,
        spec: &input::Spec,
        gate: &CancelGate,
        loader: &mut impl FnMut(&str) -> Result<fs_matdb::NormalizedInterfacePack>,
    ) -> Result<()> {
        // Commit a complete source set only after cumulative source admission.
        // Cache caller reads across all components; each part retains its own
        // exact ordered interface identity and application geometry.
        let mut candidate = self.components.clone();
        let mut glue = self.glue.clone();
        let mut packs = BTreeMap::<String, fs_matdb::NormalizedInterfacePack>::new();
        let mut cached = |name: &str| {
            if let Some(pack) = packs.get(name) {
                return Ok(pack.clone());
            }
            let pack = loader(name)?;
            packs.insert(name.into(), pack.clone());
            Ok(pack)
        };
        let mut bytes = 0usize;
        for c in &mut candidate {
            c.shell.load_interfaces(gate, &mut cached)?;
            bytes = bytes
                .checked_add(c.shell.interface_pack_bytes()?)
                .ok_or("assembly interface source bytes overflow")?;
            if bytes > input::MAX_TOTAL_PACK_BYTES {
                return Err(
                    "assembled interface packs exceed the original total source byte budget".into(),
                );
            }
        }
        glue.load(spec, gate, &mut cached)?;
        bytes = bytes.checked_add(glue.loaded_pack_bytes()?)
            .ok_or("assembly interface source bytes overflow")?;
        if bytes > input::MAX_TOTAL_PACK_BYTES {
            return Err("assembled interface packs exceed the original total source byte budget".into());
        }
        poll(gate)?;
        self.components = candidate;
        self.glue = glue;
        Ok(())
    }

    pub(crate) fn has_interface_source(&self, name: &str) -> bool {
        self.glue.has_source(name)
            || self.components.iter().any(|component| component.shell.has_interface_source(name))
    }

    pub(crate) fn interface_pack_bytes(&self) -> Result<usize> {
        self.components.iter().try_fold(self.glue.loaded_pack_bytes()?, |n, c| {
            n.checked_add(c.shell.interface_pack_bytes()?)
                .ok_or_else(|| "assembly interface source bytes overflow".into())
        })
    }

    pub(crate) fn bridge_separation_m(
        &self,
        spec: &input::Spec,
        positive: usize,
        negative: usize,
    ) -> Result<f64> {
        if self.bridges.len() != spec.strings.len() {
            return Err("assembly rotation requires one physical bridge per string".into());
        }
        let p = self
            .bridges
            .get(positive)
            .ok_or("unknown positive assembly bridge")?;
        let n = self
            .bridges
            .get(negative)
            .ok_or("unknown negative assembly bridge")?;
        if p.footprint.direction != n.footprint.direction {
            return Err(
                "assembly chord rotation requires a common physical translation direction".into(),
            );
        }
        bridge_lever(
            p.footprint.position(&self.components[p.component].shell, spec)?,
            n.footprint.position(&self.components[n.component].shell, spec)?,
            p.footprint.direction,
        )
    }

    pub(crate) fn compile(
        &self,
        spec: &input::Spec,
        sources: Vec<MaterialStringSource>,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        settings: StringBankOptions,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<ReciprocalStringBank> {
        self.compile_with_preparation(spec, sources, catalog, loader, settings,
            spec.prepared_linear.into(), cx, gate)
    }

    pub(crate) fn compile_with_preparation(
        &self, spec: &input::Spec, sources: Vec<MaterialStringSource>,
        catalog: &mut Catalog, loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        settings: StringBankOptions, preparation: ModalPreparation, cx: &Cx<'_>, gate: &CancelGate,
    ) -> Result<ReciprocalStringBank> {
        let body = self.resolve_body(spec, catalog, loader, gate)?;
        let bank = ReciprocalStringBank::from_material_shell_assembly_prepared(
            cx,
            sources,
            spec.gas,
            body,
            spec.cavity.clone(),
            preparation,
            spec.bow_law(0)?.clone(),
            settings,
            gate,
        )?;
        let image = bank.source().material_shell_assembly_image()
            .ok_or("assembled body lost its original modal image")?;
        self.glue.require_modes(image.eigenvalue_intervals(), gate)?;
        Ok(bank)
    }

    /// Project the original assembled body independently of playback preparation.
    pub(crate) fn project(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<MaterialShellAssemblyImage>
    {
        let body = self.resolve_body(spec, catalog, loader, gate)?;
        let image = body.specimen.reduce(&body.reduction, &body.force_ports, gate)?;
        self.glue.require_modes(image.eigenvalue_intervals(), gate)?;
        Ok(image)
    }

    fn resolve_body(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<SharedMaterialShellAssemblyBody>
    {
        poll(gate)?;
        self.preflight(spec)?;
        let components = self
            .components
            .iter()
            .map(|c| {
                c.shell
                    .resolve_component(&c.name, spec, catalog, loader, gate)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut edges = self.edges.clone();
        let mut faces = self.faces.clone();
        self.glue.bind(spec, &components, &mut edges, &mut faces, gate)?;
        let force_ports = self
            .bridges
            .iter()
            .map(|bridge| {
                let mut rows: Vec<_> = self
                    .components
                    .iter()
                    .map(|c| vec![0.0; 6 * c.shell.nodes.len()])
                    .collect();
                rows[bridge.component] = bridge.footprint.force_port(
                    self.components[bridge.component].shell.nodes.len())?;
                Ok(ShellAssemblyForcePort { components: rows })
            })
            .collect::<Result<Vec<_>>>()?;
        let body = SharedMaterialShellAssemblyBody {
            specimen: MaterialShellAssembly {
                components,
                edge_ties: edges,
                face_ties: faces,
                tie_budget: self.ties,
            },
            reduction: ShellAssemblyReductionOptions {
                eigenvalue_window: self.window,
                n_modes: spec.body.modes,
                projection_budget: ShellProjectionBudget {
                    max_modes: MAX_MODES,
                    max_projection_terms: self.projection_work,
                    reference_resistance_n_s_m: 1.0,
                },
                max_total_nodes: MAX_TOTAL_NODES,
                max_total_triangles: MAX_TOTAL_TRIANGLES,
            },
            force_ports,
            bridges: spec.strings.iter().map(|s| s.bridge).collect(),
            modal_budget: modal_budget(),
            listener_m: spec.body.listener,
        };
        Ok(body)
    }
}

fn component(components: &[Component], name: &str) -> Result<usize> {
    components
        .iter()
        .position(|c| c.name == name)
        .ok_or_else(|| "assembly reference names an unknown component".into())
}
fn offset(value: &str) -> Result<TieOffset> {
    match value {
        "midsurface" => Ok(TieOffset::Midsurface),
        "positive-skin" => Ok(TieOffset::Skin { side: 1 }),
        "negative-skin" => Ok(TieOffset::Skin { side: -1 }),
        _ => Err("tie offset must be midsurface, positive-skin or negative-skin".into()),
    }
}
fn rotations(value: &str) -> Result<RotationCoupling> {
    match value {
        "full" => Ok(RotationCoupling::Full),
        "translations" => Ok(RotationCoupling::TranslationsOnly),
        _ => Err("tie rotation policy must be full or translations".into()),
    }
}
fn indices(r: &mut Records<'_>, record: &str, n: usize, bound: usize) -> Result<Vec<usize>> {
    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        let i = index(r.row(record, 1)?[0], bound)?;
        if result.contains(&i) {
            return Err("tie node or triangle list repeats an index".into());
        }
        result.push(i);
    }
    Ok(result)
}
