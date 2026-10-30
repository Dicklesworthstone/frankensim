//! Explicit curved midsurface data feeding the shared material shear-shell owner.
//! Coordinates, directors, homogeneous supports and physical bridge directions
//! are data. No modal shape, damping coefficient or acoustic gain is accepted.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::shell::SharedMaterialShellBody;
use fs_couple::shear_shell::{
    compile_material_shear_shell, MaterialShellImage, ShellProjectionBudget, ShellReductionOptions, ShellRegionMaterial,
    CONSTANT_LOSS_REFERENCE_OMEGA, CONSTANT_LOSS_REFERENCE_OMEGA_DIMS, TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES,
};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
};
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::shell::survey::validate_midsurface;

#[path = "shell/interfaces.rs"]
mod interfaces;
#[path = "shell/pressure.rs"]
mod pressure;
#[path = "shell/coating.rs"]
mod coating;
#[path = "shell/assembly.rs"]
pub(super) mod assembly;

// Per-component ceilings of one shell profile. The largest component on file is the il Cannone rib set at
// 6 mm (data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h6.tsv): 2904
// nodes and 5082 triangles.
pub(super) const MAX_NODES: usize = 3072;
pub(super) const MAX_TRIANGLES: usize = 6144;
pub(super) const MAX_MODES: usize = input::MAX_BODY_MODES;
const _: () = assert!(MAX_MODES <= fs_couple::shear_shell::MAX_SHELL_MODES);
// Interface bindings keep their own allowance.
const MAX_PROJECTION_WORK: usize = 268_435_456;
// The shell projection work is `triangles * (10^6 + 324 modes^2) + 6 nodes modes` plus force and pressure
// projections (about 1.4e11 for MAX_TRIANGLES triangles and MAX_MODES modes). An assembly adds the owner's
// conservative dense mortar bound `(6 total nodes)^2 (modes + bridges)`: about 2.0e12 in all at the v4 caps
// (12288 nodes, 24576 triangles, 256 modes, 8 bridges). At 128 modes and four bridges it is 7.7e11 for the
// 11769-node il Cannone at 6 mm and 2.3e11 for the 5805-node Esteso guitar. 2^41 on 64-bit hosts.
pub(super) const MAX_SHELL_PROJECTION_WORK: usize = if (1u64 << 41) > usize::MAX as u64 { usize::MAX } else { (1u64 << 41) as usize };

#[derive(Clone, Debug)]
struct Node {
    // x/y scale with the original performance dimensions; z is in metres.
    position: [f64; 3],
    director: [f64; 3],
}
#[derive(Clone, Debug)]
struct RegionInput {
    name: String,
    pack: String,
    model: PlateMaterialModel,
    thickness_ratio: f64,
    grain_rad: f64,
    shear_correction: f64,
    axis1_global: [f64; 3],
    triangles: Vec<usize>,
}
#[derive(Clone, Debug)]
struct Bridge {
    node: usize,
    direction: [f64; 3],
}
#[derive(Clone, Debug)]
pub(super) struct Spec {
    nodes: Vec<Node>,
    triangles: Vec<[usize; 3]>,
    regions: Vec<RegionInput>,
    coatings: coating::Inputs,
    coating_memory: Option<coating::MemoryInput>,
    supports: Vec<ShellConstraint>,
    interfaces: interfaces::Inputs,
    pressure_surfaces: Option<fs_couple::shear_shell::pressure::ShellPressureSurfaces>,
    bridges: Vec<Bridge>,
    minimum_facet_dot: f64,
    mass: ShellMass,
    drilling_alpha: f64,
    window: (f64, f64),
    projection_work: usize,
}

impl Spec {
    pub(super) fn parse(bytes: &[u8]) -> Result<Self> {
        let mut r = Records { lines: std::str::from_utf8(bytes)?.lines(), line: 0 };
        r.row("frankensim-string-body-v3", 0)?;
        let result = Self::read(&mut r, None)?;
        if r.lines.next().is_some() {
            return Err("unconsumed shell body profile record".into());
        }
        Ok(result)
    }

    /// One geometry/material reader for standalone bodies and assembled parts.
    /// Components inherit the assembly's spectral budget and declare their
    /// physical bridge footprints at the assembly level.
    fn read(r: &mut Records<'_>, assembly_window: Option<((f64, f64), usize)>) -> Result<Self> {
        let minimum_facet_dot = scalar(r.row("director-min-dot", 1)?[0])?;
        if minimum_facet_dot <= 0.0 || minimum_facet_dot > 1.0 {
            return Err("shell director-min-dot must lie in (0,1]".into());
        }
        let mass = match r.row("mass", 1)?[0] {
            "consistent" => ShellMass::Consistent,
            "lumped" => ShellMass::Lumped,
            _ => return Err("shell mass must be consistent or lumped".into()),
        };
        let drilling_alpha = scalar(r.row("drilling", 1)?[0])?;
        if drilling_alpha <= 0.0 {
            return Err("shell drilling strength must be explicitly positive".into());
        }
        let (window, projection_work) = if let Some(settings) = assembly_window {
            settings
        } else {
            read_window(r)?
        };
        let coating_memory = if assembly_window.is_none() { coating::MemoryInput::read(r)? } else { None };
        let n = count(r.row("regions", 1)?[0], 1, 16)?;
        let mut regions = Vec::<RegionInput>::with_capacity(n);
        for _ in 0..n {
            let v = r.row("region", 9)?;
            if v[0].chars().any(char::is_control) || regions.iter().any(|x| x.name == v[0]) {
                return Err("shell region labels must be unique nonempty tokens".into());
            }
            let pack = v[1];
            if pack != "@body" && (pack.chars().any(char::is_control)
                || Path::new(pack).is_absolute() || pack.contains('\\')
                || Path::new(pack).components().any(|c| !matches!(c, std::path::Component::Normal(_)))) {
                return Err("shell material pack must be a relative path without traversal".into());
            }
            let thickness_ratio = scalar(v[3])?;
            let shear_correction = scalar(v[5])?;
            if thickness_ratio <= 0.0 || shear_correction <= 0.0 {
                return Err("shell region requires positive thickness and shear correction".into());
            }
            regions.push(RegionInput {
                name: v[0].into(), pack: pack.into(), model: model(v[2])?,
                thickness_ratio, grain_rad: scalar(v[4])?, shear_correction,
                axis1_global: unit_direction(&v[6..9])?,
                triangles: Vec::new(),
            });
        }
        if assembly_window.is_none() && !regions.iter().any(|x| x.pack == "@body") {
            return Err("shell regions must retain the performance body material pack".into());
        }
        let coatings = coating::Inputs::parse(r, &regions)?;
        if coating_memory.is_some() && coatings.is_empty() {
            return Err("coating-memory requires at least one declared physical film".into());
        }
        let n = count(r.row("nodes", 1)?[0], 3, MAX_NODES)?;
        let mut nodes = Vec::with_capacity(n);
        for _ in 0..n {
            let v = r.row("node", 6)?;
            let position = [scalar(v[0])?, scalar(v[1])?, scalar(v[2])?];
            if !(0.0..=1.0).contains(&position[0]) || !(0.0..=1.0).contains(&position[1]) {
                return Err("shell x/y node fractions must lie in [0,1]".into());
            }
            let director = unit_direction(&v[3..6])?;
            nodes.push(Node { position, director });
        }
        let n = count(r.row("triangles", 1)?[0], 1, MAX_TRIANGLES)?;
        let mut triangles = Vec::with_capacity(n);
        for element in 0..n {
            let v = r.row("triangle", 4)?;
            let t = [index(v[0], nodes.len())?, index(v[1], nodes.len())?, index(v[2], nodes.len())?];
            let region = regions.iter_mut().find(|x| x.name == v[3])
                .ok_or("shell triangle refers to an unknown material region")?;
            region.triangles.push(element);
            triangles.push(t);
        }
        if regions.iter().any(|x| x.triangles.is_empty()) {
            return Err("every declared shell material region must own a triangle".into());
        }
        let n = count(r.row("supports", 1)?[0], 0, nodes.len())?;
        let mut supports = Vec::<ShellConstraint>::with_capacity(n);
        for _ in 0..n {
            let v = r.row("support", 2)?;
            let node = index(v[0], nodes.len())?;
            let bytes = v[1].as_bytes();
            if supports.iter().any(|x| x.node == node) || bytes.len() != 6
                || bytes.iter().any(|x| *x != b'0' && *x != b'1') || !bytes.contains(&b'1') {
                return Err("shell support needs a unique node and six explicit 0/1 DOF flags".into());
            }
            let mask = core::array::from_fn(|i| bytes[i] == b'1');
            supports.push(ShellConstraint { node, mask, frame: ConstraintFrame::Global });
        }
        let interfaces = interfaces::Inputs::parse(r, &nodes, &triangles, &regions)?;
        if assembly_window.is_none() && supports.is_empty() && interfaces.connection_count() == 0 {
            return Err("shell body requires explicit supports or material interfaces".into());
        }
        let pressure_surfaces = pressure::parse(r, triangles.len())?;
        if assembly_window.is_some() && pressure_surfaces.is_none() {
            return Err("every assembled shell component must declare its pressure surfaces".into());
        }
        let n = if assembly_window.is_some() { 0 } else { count(r.row("bridges", 1)?[0], 1, 8)? };
        let mut bridges = Vec::with_capacity(n);
        for _ in 0..n {
            let v = r.row("bridge", 4)?;
            let node = index(v[0], nodes.len())?;
            let direction = unit_direction(&v[1..4])?;
            if let Some(support) = supports.iter().find(|s| s.node == node) {
                if !(0..3).any(|i| direction[i] != 0.0 && !support.mask[i]) {
                    return Err("shell bridge has no unconstrained translational direction".into());
                }
            }
            bridges.push(Bridge { node, direction });
        }
        Ok(Self { nodes, triangles, regions, coatings, coating_memory, supports, interfaces, pressure_surfaces, bridges, minimum_facet_dot,
            mass, drilling_alpha, window, projection_work })
    }

    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if self.bridges.len() != spec.strings.len() {
            return Err("shell profile requires one physical bridge per actual string".into());
        }
        if spec.body.modes == 0 || spec.body.modes > MAX_MODES {
            return Err(format!("shell body profile admits 1..={MAX_MODES} retained body modes").into());
        }
        self.geometry(spec)?;
        if let Some(memory) = &self.coating_memory { memory.preflight(spec)?; }
        self.interfaces.preflight(self, spec)?;
        Ok(())
    }

    pub(super) fn additional_body_connections(&self) -> usize {
        self.interfaces.connection_count()
    }

    pub(super) fn material_memory_states(&self) -> usize {
        self.coating_memory.as_ref().map_or(0, |memory| memory.options.max_memory_coordinates)
    }

    pub(super) fn bind_material_memory(&self, bank: ReciprocalStringBank, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        Ok(if self.coating_memory.is_some() {
            bank.with_material_shell_coating_relaxation(None, gate)?
        } else { bank })
    }

    pub(super) fn load_interfaces(&mut self, gate: &CancelGate,
        loader: &mut impl FnMut(&str) -> Result<fs_matdb::NormalizedInterfacePack>) -> Result<()>
    {
        self.interfaces.load(gate, loader)
    }

    pub(super) fn interface_pack_bytes(&self) -> Result<usize> {
        self.interfaces.loaded_pack_bytes()
    }

    fn position(&self, node: usize, spec: &input::Spec) -> [f64; 3] {
        let p = self.nodes[node].position;
        [p[0] * spec.body.length, p[1] * spec.body.width, p[2]]
    }

    fn geometry(&self, spec: &input::Spec) -> Result<(ShellMesh, NodalDirectors)> {
        let mesh = ShellMesh::new((0..self.nodes.len()).map(|i| self.position(i, spec)).collect(),
            self.triangles.clone())?;
        validate_midsurface(&mesh, ProfileBudget { max_nodes: MAX_NODES,
            max_triangles: MAX_TRIANGLES, max_feature_evaluations: 0 })?;
        let directors = NodalDirectors::new(&mesh,
            self.nodes.iter().map(|n| n.director).collect(), self.minimum_facet_dot)?;
        Ok((mesh, directors))
    }

    pub(super) fn bridge_separation_m(&self, spec: &input::Spec, positive: usize, negative: usize) -> Result<f64> {
        if self.bridges.len() != spec.strings.len() {
            return Err("shell rotational joint requires one physical bridge per string".into());
        }
        let p = self.bridges.get(positive).ok_or("unknown positive shell bridge")?;
        let n = self.bridges.get(negative).ok_or("unknown negative shell bridge")?;
        if p.direction != n.direction {
            return Err("shell chord rotation requires a common physical translation direction".into());
        }
        let p = self.position(p.node, spec);
        let n = self.position(n.node, spec);
        let d = self.bridges[positive].direction;
        bridge_lever(p, n, d)
    }

    pub(super) fn compile(&self, spec: &input::Spec, sources: Vec<MaterialStringSource>,
        catalog: &mut Catalog, loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        settings: StringBankOptions, cx: &Cx<'_>, gate: &CancelGate) -> Result<ReciprocalStringBank>
    {
        self.compile_with_preparation(spec, sources, catalog, loader, settings,
            spec.prepared_linear.into(), cx, gate)
    }

    pub(super) fn compile_with_preparation(&self, spec: &input::Spec, sources: Vec<MaterialStringSource>,
        catalog: &mut Catalog, loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        settings: StringBankOptions, preparation: ModalPreparation, cx: &Cx<'_>, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        let body = self.resolve_body(spec, catalog, loader, gate)?;
        // The prepared record, when the data declares one, must reach the constructor: a body above the ordinary
        // 42-coordinate block cannot be built without it.
        Ok(ReciprocalStringBank::from_material_shear_shell_prepared(cx, sources, spec.gas, body,
            spec.cavity.clone(), preparation, spec.bow_law(0)?.clone(), settings, gate)?)
    }

    /// Resolve the same physical body without constructing strings or a time propagator.
    pub(super) fn project(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<MaterialShellImage>
    {
        let body = self.resolve_body(spec, catalog, loader, gate)?;
        let reduction = ShellReductionOptions {
            constraints: body.constraints, mass: body.mass, drilling_alpha: body.drilling_alpha,
            eigenvalue_window: body.eigenvalue_window, n_modes: body.n_modes,
            assembly_budget: body.assembly_budget, projection_budget: body.projection_budget,
        };
        Ok(body.specimen.reduce(&body.directors, &reduction, &body.force_ports, gate)?)
    }

    fn resolve_body(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<SharedMaterialShellBody>
    {
        poll(gate)?;
        self.preflight(spec)?;
        let component = self.resolve_component("body", spec, catalog, loader, gate)?;
        let force_ports = self.bridges.iter().map(|bridge| {
            let mut port = vec![0.0; 6 * self.nodes.len()];
            port[6 * bridge.node..6 * bridge.node + 3].copy_from_slice(&bridge.direction);
            port
        }).collect();
        let body = SharedMaterialShellBody {
            specimen: component.specimen, directors: component.directors,
            constraints: component.constraints, mass: component.mass,
            drilling_alpha: component.drilling_alpha, eigenvalue_window: self.window,
            n_modes: spec.body.modes, force_ports,
            bridges: spec.strings.iter().map(|s| s.bridge).collect(),
            assembly_budget: component.assembly_budget,
            projection_budget: ShellProjectionBudget { max_modes: MAX_MODES,
                max_projection_terms: self.projection_work, reference_resistance_n_s_m: 1.0 },
            modal_budget: modal_budget(), listener_m: spec.body.listener,
        };
        Ok(body)
    }

    fn resolve_component(&self, name: &str, spec: &input::Spec,
        catalog: &mut Catalog, loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        gate: &CancelGate) -> Result<fs_couple::shear_shell::assembly::MaterialShellComponent>
    {
        poll(gate)?;
        let (mesh, directors) = self.geometry(spec)?;
        let assembly_budget = MitcBudget { max_nodes: MAX_NODES, max_triangles: MAX_TRIANGLES };
        let mut regions = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            poll(gate)?;
            let pack = if region.pack == "@body" { &spec.body.pack } else { &region.pack };
            let mut requested = requirements(region.model)?;
            if !matches!(region.model, PlateMaterialModel::Isotropic) {
                for key in ["shear_modulus_13", "shear_modulus_23"] {
                    requested.push(Requirement::try_new(key, Pressure::DIMS, Admit::StrictlyPositive)?);
                }
            }
            for key in TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES {
                requested.push(Requirement::try_new(key, DynViscosity::DIMS, Admit::Finite)?);
            }
            let mut material = resolve_properties(pack, catalog, loader, spec, &requested)?;
            // A pack that states its viscosities' constant-loss anchor is read with it; it is never assumed.
            if !catalog[pack.as_str()].pack.claims().claims_for(CONSTANT_LOSS_REFERENCE_OMEGA).is_empty() {
                requested.push(Requirement::try_new(CONSTANT_LOSS_REFERENCE_OMEGA,
                    CONSTANT_LOSS_REFERENCE_OMEGA_DIMS, Admit::StrictlyPositive)?);
                material = resolve_properties(pack, catalog, loader, spec, &requested)?;
            }
            let thickness = spec.body.thickness * region.thickness_ratio;
            if !thickness.is_finite() || thickness <= 0.0 {
                return Err("shell material thickness is unrepresentable".into());
            }
            regions.push(ShellRegionMaterial {
                region: PlateRegionMaterial {
                    region: PlateRegion { name: region.name.clone(), triangle_indices: region.triangles.clone() },
                    material, model: region.model, material_angle_rad: region.grain_rad,
                    thickness_constraint: PlateThicknessConstraint::FixedThickness(thickness),
                },
                shear_correction: region.shear_correction, axis1_global: region.axis1_global,
            });
        }
        let specimen = compile_material_shear_shell(mesh, regions, assembly_budget, gate)?;
        let specimen = self.coatings.bind(specimen, spec, catalog, loader, gate)?;
        let specimen = if let Some(memory) = &self.coating_memory {
            specimen.with_coating_relaxation(memory.options, gate)?
        } else { specimen };
        let specimen = self.interfaces.bind(specimen, gate)?;
        let specimen = if let Some(surfaces) = &self.pressure_surfaces {
            specimen.with_pressure_surfaces(surfaces.clone(), gate)?
        } else { specimen };
        Ok(fs_couple::shear_shell::assembly::MaterialShellComponent {
            name: name.into(), specimen, directors, constraints: self.supports.clone(),
            mass: self.mass, drilling_alpha: self.drilling_alpha, assembly_budget,
        })
    }
}

fn read_window(r: &mut Records<'_>) -> Result<((f64, f64), usize)> {
    let v = r.row("window", 2)?;
    let window = (scalar(v[0])?, scalar(v[1])?);
    if window.0 < 0.0 || window.1 <= window.0 {
        return Err("shell eigenvalue window must be nonnegative and increasing".into());
    }
    let work = count(r.row("projection-work", 1)?[0], 1, MAX_SHELL_PROJECTION_WORK)?;
    Ok((window, work))
}

fn bridge_lever(p: [f64; 3], n: [f64; 3], d: [f64; 3]) -> Result<f64> {
    let chord = [p[0] - n[0], p[1] - n[1], p[2] - n[2]];
    // d dot (theta cross chord) = theta dot (chord cross d): the
    // work-conjugate angular coordinate uses the perpendicular lever.
    let lever = norm([
        chord[1] * d[2] - chord[2] * d[1],
        chord[2] * d[0] - chord[0] * d[2],
        chord[0] * d[1] - chord[1] * d[0],
    ]);
    let chord_length = norm(chord);
    if !lever.is_finite() || !chord_length.is_finite() || lever <= 0.0
        || lever <= 64.0 * f64::EPSILON * chord_length {
        return Err("shell rotational bridges need a resolved lever perpendicular to their force direction".into());
    }
    Ok(lever)
}

fn count(value: &str, min: usize, max: usize) -> Result<usize> {
    let n = value.parse()?;
    if n < min || n > max { return Err("shell profile count exceeds its explicit range".into()); }
    Ok(n)
}
fn index(value: &str, length: usize) -> Result<usize> {
    count(value, 0, length - 1)
}
fn norm(v: [f64; 3]) -> f64 {
    fs_math::det::hypot(fs_math::det::hypot(v[0], v[1]), v[2])
}
fn unit_direction(values: &[&str]) -> Result<[f64; 3]> {
    let v = [scalar(values[0])?, scalar(values[1])?, scalar(values[2])?];
    // Match the original director owner's squared-unit tolerance; retain the
    // supplied components without normalization or inferred smoothing.
    let square: f64 = v.iter().map(|x| x * x).sum();
    if !square.is_finite() || (square - 1.0).abs() > 1e-12 {
        return Err("shell direction must be an explicitly supplied unit vector".into());
    }
    Ok(v)
}
