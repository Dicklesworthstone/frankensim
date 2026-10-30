//! Source/data adapter for directional elasticity, viscosity and material memory.
//! Every physical operator, tensor rotation, mode, bridge and pressure projection
//! remains in the existing plate/shell/material/reciprocal owners.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::chart::{SharedMaterialPlateChartBody, PlateBridgeFootprints};
use fs_couple::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
use fs_couple::thin_plate::{compile_plate_material_chart, PlateChartRadiation, PlateRegionMaterial};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as VISCOSITY;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::one_way::{
    OneWaySetupStore, ShellAcousticSkin,
};

pub const MAX_BYTES: usize = 65_536;
/// Shell (v3) and assembled (v4) profiles carry node and triangle rows, so they take a larger byte
/// ceiling; flat plate profiles keep [`MAX_BYTES`]. An assembly at its own record caps (8192 nodes, 16384
/// triangles, both pressure roles on every triangle) writes about 2.1 MB at full float precision.
pub const MAX_SHELL_BYTES: usize = 4 * 1_048_576;
const MAX_GRID: usize = 12;

#[path = "body/shell.rs"]
mod shell;
#[path = "body/relaxation.rs"]
mod relaxation;
#[path = "body/cavity.rs"]
mod cavity;
#[path = "body/skin.rs"]
mod skin;
use relaxation::PlateLoss;
/// The choices `material_strings assemble-shell` reads for a shell-assembly table.
pub(crate) use shell::assembly::table::OptionsFile as ShellTableOptions;

/// An optional complete body-model selection; the performance still owns its
/// material pack, physical dimensions, thickness, mode count and listener.
#[derive(Clone, Debug)]
pub struct Profile {
    pub bytes: Vec<u8>,
    model: BodyModel,
    interior: Option<cavity::Spec>,
    skin: Option<skin::Spec>,
    observer_store: Option<skin::SetupDir>,
}
#[derive(Clone, Debug)]
enum BodyModel {
    Plate(PlateProfile),
    Shell(shell::Spec),
    Assembly(shell::assembly::Spec),
}

/// The original retained structural image, without a string bank or propagator.
pub(super) enum ProjectedShellBody {
    Shell(fs_couple::shear_shell::MaterialShellImage),
    Assembly(fs_couple::shear_shell::assembly::MaterialShellAssemblyImage),
}
impl ProjectedShellBody {
    pub(super) fn shell(&self) -> Option<&fs_couple::shear_shell::MaterialShellImage> {
        match self { Self::Shell(image) => Some(image), Self::Assembly(_) => None }
    }
    pub(super) fn assembly(&self) -> Option<&fs_couple::shear_shell::assembly::MaterialShellAssemblyImage> {
        match self { Self::Assembly(image) => Some(image), Self::Shell(_) => None }
    }
}
impl Profile {
    pub fn has_planar_material_basis(&self) -> bool { matches!(self.model, BodyModel::Plate(_)) }
    pub fn has_shell_material_basis(&self) -> bool {
        matches!(self.model, BodyModel::Shell(_) | BodyModel::Assembly(_))
    }

    /// Use the same source resolution and eigensolve as playback, without its
    /// strings, cavity, bow, acoustic loading or prepared transition matrices.
    pub(super) fn project_shell(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<ProjectedShellBody>
    {
        poll(gate)?;
        match &self.model {
            BodyModel::Shell(shell) => Ok(ProjectedShellBody::Shell(shell.project(spec, catalog, loader, gate)?)),
            BodyModel::Assembly(assembly) => Ok(ProjectedShellBody::Assembly(assembly.project(spec, catalog, loader, gate)?)),
            BodyModel::Plate(_) => Err("body-response requires a curved shell or assembled shell body".into()),
        }
    }

    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_SHELL_BYTES { return Err("body profile exceeds 4 MiB".into()); }
        let text = std::str::from_utf8(bytes)?;
        let (body_text, interior) = cavity::split(text)?;
        let (body_text, skin) = skin::split(body_text)?;
        let body_bytes = body_text.as_bytes();
        let version = text.lines().next().and_then(|s| s.split_ascii_whitespace().next())
            .ok_or("missing body profile version")?;
        let model = match version {
            "frankensim-string-body-v1" | "frankensim-string-body-v2" =>
                BodyModel::Plate(PlateProfile::parse(body_bytes)?),
            "frankensim-string-body-v3" => BodyModel::Shell(shell::Spec::parse(body_bytes)?),
            "frankensim-string-body-v4" => BodyModel::Assembly(shell::assembly::Spec::parse(body_bytes)?),
            _ => return Err("unsupported body profile version".into()),
        };
        if interior.is_some() && matches!(model, BodyModel::Plate(_)) {
            return Err("interior FE requires a curved shell or assembled shell body profile".into());
        }
        if skin.is_some() && matches!(model, BodyModel::Plate(_)) {
            return Err("an acoustic skin requires a curved shell or assembled shell body profile".into());
        }
        Ok(Self { bytes: bytes.to_vec(), model, interior, skin, observer_store: None })
    }
    pub fn has_acoustic_skin(&self) -> bool { self.skin.is_some() }
    pub fn has_spatial_cavity(&self) -> bool { self.interior.is_some() }
    /// The declared closed acoustic skin in the performance body's frame, if any.
    pub fn acoustic_skin(&self, spec: &input::Spec) -> Option<ShellAcousticSkin> {
        self.skin.as_ref().map(|skin| skin.resolve(spec))
    }
    /// Keep the skin's one-way observer setups in the existing directory `dir` (`--observer-store`).
    pub fn keep_observer_setups(&mut self, dir: &std::path::Path) -> Result<()> {
        if self.skin.is_none() {
            return Err("--observer-store requires a body profile that declares its closed acoustic skin".into());
        }
        self.observer_store = Some(skin::SetupDir::open(dir)?);
        Ok(())
    }
    pub fn observer_store(&self) -> Option<&dyn OneWaySetupStore> {
        self.observer_store.as_ref().map(|store| store as &dyn OneWaySetupStore)
    }
    /// Where the setup under `key` is kept, with `--observer-store`.
    pub fn observer_store_entry(&self, key: &fs_blake3::ContentHash) -> Option<std::path::PathBuf> {
        self.observer_store.as_ref().map(|store| store.entry(key))
    }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        spec.validate_cavity_gas()?;
        if spec.cavity_gas_record.is_some() && spec.cavity.is_none() && self.interior.is_none() {
            return Err("cavity_gas requires an actual lumped or spatial body cavity".into());
        }
        match &self.model {
            BodyModel::Plate(_) if spec.body.modes > input::MAX_PLATE_BODY_MODES => Err(format!(
                "flat plate body profiles admit 1..={} retained body modes", input::MAX_PLATE_BODY_MODES).into()),
            BodyModel::Plate(plate) => plate.loss.preflight(spec),
            BodyModel::Shell(shell) => shell.preflight(spec),
            BodyModel::Assembly(assembly) => assembly.preflight(spec),
        }?;
        if let Some(interior) = &self.interior { interior.preflight(spec)?; }
        match &self.skin {
            Some(skin) => skin.preflight(spec)?,
            None if spec.one_way_observer.is_some() => return Err(
                "one_way_observer requires a body profile that declares its closed acoustic skin".into()),
            None if matches!(spec.body_radiation, Some(input::BodyRadiation::Skin(_))) => return Err(
                "a computed skin load requires a body profile that declares its closed acoustic skin".into()),
            None => {}
        }
        if self.observer_store.is_some() && spec.one_way_observer.is_none() {
            return Err("--observer-store keeps one_way_observer setups, and this performance has none".into());
        }
        spec.network_dimensions_with_body_modes(true, self.additional_body_connections(),
            self.additional_body_modes())?;
        Ok(())
    }
    /// Conservative interface, pressure and internal-loss connection allowance
    /// before material-pack I/O.
    pub fn additional_body_connections(&self) -> usize {
        let structural = match &self.model {
            BodyModel::Plate(_) => 0,
            BodyModel::Shell(shell) => shell.additional_body_connections(),
            BodyModel::Assembly(assembly) => assembly.additional_body_connections(),
        };
        structural + self.interior.as_ref().map_or(0, |air|
            air.pressure_allowance() + air.opening_count())
    }
    /// Every possible positive-frequency pressure mode has an air inertia;
    /// the exact uniform mode has only its compression spring. Each declared
    /// opening contributes one original compact-neck fluid inertia.
    pub fn additional_body_modes(&self) -> usize {
        self.interior.as_ref().map_or(0, |air|
            air.pressure_allowance() - 1 + air.opening_count())
    }
    pub fn has_thermal_cavity(&self) -> bool {
        self.interior.as_ref().is_some_and(cavity::Spec::has_thermal_loss)
    }
    /// Load original ordered interface packs before the ordinary bulk loader.
    /// Geometry and complete-network admission precede the first source read.
    pub fn load_interfaces(&mut self, spec: &input::Spec, gate: &CancelGate,
        mut loader: impl FnMut(&str) -> Result<fs_matdb::NormalizedInterfacePack>) -> Result<()>
    {
        poll(gate)?;
        self.preflight(spec)?;
        match &mut self.model {
            BodyModel::Plate(_) => {},
            BodyModel::Shell(shell) => shell.load_interfaces(gate, &mut loader)?,
            BodyModel::Assembly(assembly) => assembly.load_interfaces(spec, gate, &mut loader)?,
        }
        Ok(())
    }
    pub fn interface_pack_bytes(&self) -> Result<usize> {
        match &self.model {
            BodyModel::Plate(_) => Ok(0),
            BodyModel::Shell(shell) => shell.interface_pack_bytes(),
            BodyModel::Assembly(assembly) => assembly.interface_pack_bytes(),
        }
    }
    pub fn bridge_separation_m(&self, spec: &input::Spec, positive: usize, negative: usize) -> Result<f64> {
        match &self.model {
            BodyModel::Plate(plate) => plate.bridge_separation_m(spec, positive, negative),
            BodyModel::Shell(shell) => shell.bridge_separation_m(spec, positive, negative),
            BodyModel::Assembly(assembly) => assembly.bridge_separation_m(spec, positive, negative),
        }
    }
    pub fn compile(&self, spec:&input::Spec, sources:Vec<MaterialStringSource>,
        catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<MaterialPack>,
        settings:StringBankOptions, cx:&Cx<'_>, gate:&CancelGate) -> Result<ReciprocalStringBank>
    {
        let bank = match &self.model {
            BodyModel::Plate(plate) => plate.compile(spec, sources, catalog, loader, settings, cx, gate),
            BodyModel::Shell(shell) => shell.compile(spec, sources, catalog, loader, settings, cx, gate),
            BodyModel::Assembly(assembly) => assembly.compile(spec, sources, catalog, loader, settings, cx, gate),
        }?;
        // Original string/body material histories are the fixed prefix of the
        // shared memory owner. Cavity thermal diffusion appends to that prefix
        // before rotational and nonlinear coordinates are attached.
        let bank = if self.has_thermal_cavity() {
            // Explicit layered memory retains the final mass-normalized string
            // basis, so source-gas inertia must precede that selection.
            let bank = match spec.cylinder_gas_loading {
                Some(resolution) => bank.with_cylinder_gas_loading(resolution, gate)?,
                None => bank,
            };
            let bank = self.bind_material_memory(bank, gate)?;
            match &spec.bending_memory {
                Some(memory) => memory.bind(bank, gate)?,
                None => bank,
            }
        } else { bank };
        match &self.interior {
            Some(interior) => interior.bind(bank, spec, gate),
            None => Ok(bank),
        }
    }
    /// Reserve body and string histories together before material-pack I/O.
    pub fn material_memory_states(&self) -> usize {
        let solid = match &self.model {
            BodyModel::Plate(plate) => plate.loss.max_memory_states(),
            BodyModel::Shell(shell) => shell.material_memory_states(),
            BodyModel::Assembly(assembly) => assembly.material_memory_states(),
        };
        solid + self.interior.as_ref().map_or(0, cavity::Spec::thermal_memory_states)
    }
    /// All solid/gas/axial coordinates are now fixed; memory enters the SAME
    /// coupled mechanical owner before conversion to multiple bows and PCM.
    pub fn bind_material_memory(&self, bank: ReciprocalStringBank, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        match &self.model {
            BodyModel::Plate(plate) => plate.loss.bind(bank, gate),
            BodyModel::Shell(shell) => shell.bind_material_memory(bank, gate),
            BodyModel::Assembly(assembly) => assembly.bind_material_memory(bank, gate),
        }
    }
}
#[derive(Clone, Debug)]
struct PlateProfile {
    nx: usize, ny: usize,
    layout: Layout,
    support: EdgeSupport,
    clamp_x_min: bool,
    pretension_n_m: f64,
    window: (f64, f64),
    loss: PlateLoss,
}
#[derive(Clone, Debug)]
enum Layout {
    Uniform { model:PlateMaterialModel, grain_rad:f64, bridge:(usize,usize) },
    Regional { regions:Vec<Region>, bridges:Vec<(usize,usize)> },
}
#[derive(Clone, Debug)]
struct Region {
    name:String, pack:String, model:PlateMaterialModel, thickness_ratio:f64, grain_rad:f64,
    // Half-open structured CELL ranges, never independent meshes or glued parts.
    x:(usize,usize), y:(usize,usize),
}
impl PlateProfile {
    fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_BYTES { return Err("body profile exceeds 64 KiB".into()); }
        let mut input = Records { lines: std::str::from_utf8(bytes)?.lines(), line: 0 };
        let mut version=input.lines.next().ok_or("missing body profile version")?.split_ascii_whitespace();
        input.line=1;
        let regional=match version.next() {
            Some("frankensim-string-body-v1")=>false,
            Some("frankensim-string-body-v2")=>true,
            _=>return Err("unsupported body profile version".into()),
        };
        if version.next().is_some() { return Err("unconsumed body profile version field".into()); }
        let v = input.row("mesh", 2)?;
        let nx: usize = v[0].parse()?; let ny: usize = v[1].parse()?;
        if !(1..=MAX_GRID).contains(&nx) || !(1..=MAX_GRID).contains(&ny) {
            return Err("body grid requires 1..=12 intervals along each axis".into());
        }
        let uniform=if regional { None } else {
            Some((model(input.row("elastic",1)?[0])?,scalar(input.row("grain",1)?[0])?))
        };
        let (support,clamp_x_min) = match input.row("support", 1)?[0] {
            "simply-supported" => (EdgeSupport::SimplySupported,false),
            "clamped" => (EdgeSupport::Clamped,false),
            "clamped-x-min" => (EdgeSupport::Clamped,true),
            _ => return Err("body support must be simply-supported, clamped or clamped-x-min".into()),
        };
        if !clamp_x_min && (nx<2 || ny<2) {
            return Err("a perimeter-supported body requires at least 2 intervals along each axis".into());
        }
        let pretension_n_m = scalar(input.row("pretension", 1)?[0])?;
        if pretension_n_m < 0.0 { return Err("body profile requires nonnegative conservative pretension".into()); }
        let v = input.row("window", 2)?;
        let window = (scalar(v[0])?, scalar(v[1])?);
        if window.0 < 0.0 || window.1 <= window.0 { return Err("body eigenvalue window must be nonnegative and increasing".into()); }
        let layout=if let Some((model,grain_rad))=uniform {
            Layout::Uniform { model,grain_rad,bridge:read_bridge(&mut input,nx,ny,clamp_x_min)? }
        } else { read_regions(&mut input,nx,ny,clamp_x_min)? };
        let record=input.lines.next();
        let loss = PlateLoss::parse(record,&mut input)?;
        if input.lines.next().is_some() { return Err("unconsumed body profile record".into()); }
        Ok(Self { nx,ny,layout,support,clamp_x_min,pretension_n_m,window,loss })
    }

    /// Actual length of the ordered point-footprint chord [m]. Its order sets
    /// the sign of the scalar rocking coordinate; no 3D axis is inferred here.
    pub fn bridge_separation_m(&self, spec: &input::Spec, positive: usize, negative: usize) -> Result<f64> {
        let Layout::Regional { bridges, .. } = &self.layout else {
            return Err("rotational bridges require a v2 body profile with separate physical footprints".into());
        };
        if bridges.len() != spec.strings.len() {
            return Err("rotational body profile must retain one footprint per string".into());
        }
        let p = bridges.get(positive).ok_or("unknown positive body footprint")?;
        let n = bridges.get(negative).ok_or("unknown negative body footprint")?;
        let dx = (p.0 as f64 - n.0 as f64) * (spec.body.length / self.nx as f64);
        let dy = (p.1 as f64 - n.1 as f64) * (spec.body.width / self.ny as f64);
        let separation = fs_math::det::hypot(dx, dy);
        if !separation.is_finite() || separation <= 0.0 {
            return Err("rotational body footprints must have a finite nonzero physical separation".into());
        }
        Ok(separation)
    }

    /// Build the same original elastic basis for every bridge/radiation map.
    /// Tensor viscosity is projected here; selected material memory is attached
    /// after the complete mechanical coordinates are fixed. No scalar loss is added.
    pub fn compile(&self, spec:&input::Spec, sources:Vec<MaterialStringSource>,
        catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<MaterialPack>,
        settings:StringBankOptions, cx:&Cx<'_>, gate:&CancelGate) -> Result<ReciprocalStringBank>
    {
        poll(gate)?;
        let inputs=self.region_inputs();
        let mesh = PlateMesh::rectangle(spec.body.length,spec.body.width,self.nx,self.ny);
        let nodes = mesh.nodes.len(); let elements = mesh.tris.len();
        let spatial=match &self.layout {
            Layout::Uniform { .. }=>None,
            Layout::Regional { bridges,.. }=> {
                if bridges.len()!=spec.strings.len() {
                    return Err("body profile requires one bridge node per actual material string".into());
                }
                Some(PlateBridgeFootprints { weights:bridges.iter().map(|&b|self.footprint(b)).collect(),
                    max_projection_terms:input::MAX_WORK as usize })
            },
        };
        let reference=match &self.layout {
            Layout::Uniform { bridge,.. }=>*bridge,
            Layout::Regional { bridges,.. }=>bridges[0],
        };
        let mut regions=Vec::with_capacity(inputs.len());
        for region in inputs {
            poll(gate)?;
            let pack=if region.pack=="@body" { &spec.body.pack } else { &region.pack };
            let material=resolve_properties(pack,catalog,loader,spec,&self.loss.requirements(region.model)?)?;
            let thickness=spec.body.thickness*region.thickness_ratio;
            if !thickness.is_finite() || thickness<=0.0 { return Err("regional body thickness is unrepresentable".into()); }
            regions.push(PlateRegionMaterial {
                region:PlateRegion { name:region.name.clone(),triangle_indices:self.triangles(&region) },
                material,model:region.model,material_angle_rad:region.grain_rad,
                thickness_constraint:PlateThicknessConstraint::FixedThickness(thickness),
            });
        }
        let fixed=if self.clamp_x_min {
            (0..=self.ny).map(|j|j*(self.nx+1)).collect()
        } else {PlateMesh::rectangle_boundary(self.nx,self.ny)};
        let specimen=compile_plate_material_chart(mesh,fixed,regions)?;
        let viscosity = if self.loss.is_tensor() {
            Some(PlateViscosity::from_chart(&specimen,VISCOSITY,gate)?)
        } else { None };
        let body = SharedMaterialPlateChartBody {
            specimen,reduction:PlateChartRadiation {
                assembly:AssemblyOptions { pretension:self.pretension_n_m,support:self.support },
                eigenvalue_window:self.window,n_modes:spec.body.modes,damping_ratio:0.0,
                unit_force_weights:self.footprint(reference),
            },
            bridges:spec.strings.iter().map(|w|w.bridge).collect(),
            max_plate_dofs:3*nodes,max_plate_elements:elements,modal_budget:modal_budget(),listener_m:spec.body.listener,
        };
        let limits=PlateViscosityLimits { max_modes:16,max_projection_terms:input::MAX_WORK as usize,
            reference_resistance_n_s_m:1.0 };
        Ok(match (spatial, viscosity) {
            (Some(footprints), Some(viscosity)) =>
                ReciprocalStringBank::from_viscous_material_plate_chart_with_force_ports(cx,sources,spec.gas,
                    body,footprints,viscosity,limits,spec.cavity.clone(),spec.bow_law(0)?.clone(),settings,gate)?,
            (None, Some(viscosity)) =>
                ReciprocalStringBank::from_viscous_material_plate_chart(cx,sources,spec.gas,body,viscosity,
                    limits,spec.cavity.clone(),spec.bow_law(0)?.clone(),settings,gate)?,
            (Some(footprints), None) =>
                ReciprocalStringBank::from_material_plate_chart_with_force_ports(cx,sources,spec.gas,
                    body,footprints,spec.cavity.clone(),spec.bow_law(0)?.clone(),settings,gate)?,
            (None, None) =>
                ReciprocalStringBank::from_material_plate_chart(cx,sources,spec.gas,body,
                    spec.cavity.clone(),spec.bow_law(0)?.clone(),settings,gate)?,
        })
    }
    fn footprint(&self,node:(usize,usize))->Vec<f64> {
        let mut weights=vec![0.0;(self.nx+1)*(self.ny+1)];
        weights[node.1*(self.nx+1)+node.0]=1.0;weights
    }
    fn region_inputs(&self)->Vec<Region> {
        match &self.layout {
            Layout::Uniform { model,grain_rad,.. }=>vec![Region {name:"body".into(),pack:"@body".into(),
                model:*model,thickness_ratio:1.0,grain_rad:*grain_rad,x:(0,self.nx),y:(0,self.ny)}],
            Layout::Regional { regions,.. }=>regions.clone(),
        }
    }
    fn triangles(&self,region:&Region)->Vec<usize> {
        // PlateMesh::rectangle emits two triangles per cell, x varying fastest.
        (region.y.0..region.y.1).flat_map(|j|(region.x.0..region.x.1)
            .flat_map(move |i|[2*(j*self.nx+i),2*(j*self.nx+i)+1])).collect()

    }
}

fn read_bridge(input:&mut Records<'_>,nx:usize,ny:usize,clamp_x_min:bool)->Result<(usize,usize)> {
    let v=input.row("bridge",2)?;let b=(v[0].parse()?,v[1].parse()?);
    if b.0==0 || b.0>nx || b.1>ny || !clamp_x_min && (b.0==nx || b.1==0 || b.1==ny) {
        return Err("body bridge must name an actual unconstrained transverse grid node".into());
    }
    Ok(b)
}
fn read_regions(input:&mut Records<'_>,nx:usize,ny:usize,clamp_x_min:bool)->Result<Layout> {
    let n:usize=input.row("regions",1)?[0].parse()?;
    if n==0 || n>16 {return Err("body profile requires 1..=16 explicit regions".into());}
    let mut regions=Vec::<Region>::with_capacity(n);let mut occupied=vec![false;nx*ny];
    let mut uses_body=false;
    for _ in 0..n {
        let v=input.row("region",9)?;
        if v[0].is_empty() || v[0].chars().any(char::is_control) || regions.iter().any(|r|r.name==v[0]) {
            return Err("body region labels must be unique nonempty tokens".into());
        }
        if v[1]!="@body" && (v[1].chars().any(char::is_control) || Path::new(v[1]).is_absolute() || v[1].contains('\\')
            || Path::new(v[1]).components().any(|c|!matches!(c,std::path::Component::Normal(_)))) {
            return Err("regional material pack must be a relative path without traversal".into());
        }
        let region=Region {name:v[0].into(),pack:v[1].into(),model:model(v[2])?,
            thickness_ratio:scalar(v[3])?,grain_rad:scalar(v[4])?,
            x:(v[5].parse()?,v[6].parse()?),y:(v[7].parse()?,v[8].parse()?)};
        if region.thickness_ratio<=0.0 || region.x.0>=region.x.1 || region.x.1>nx
            || region.y.0>=region.y.1 || region.y.1>ny {
            return Err("body region needs positive thickness and a nonempty in-grid cell rectangle".into());
        }
        for j in region.y.0..region.y.1 { for i in region.x.0..region.x.1 {
            if std::mem::replace(&mut occupied[j*nx+i],true) {return Err("body regions overlap".into());}
        } }
        uses_body|=region.pack=="@body";regions.push(region);
    }
    if !uses_body || occupied.contains(&false) {
        return Err("body regions must cover every cell exactly once and use the performance body pack".into());
    }
    let n:usize=input.row("bridges",1)?[0].parse()?;
    if n==0 || n>8 {return Err("body profile requires 1..=8 ordered bridge nodes".into());}
    let bridges=(0..n).map(|_|read_bridge(input,nx,ny,clamp_x_min)).collect::<Result<Vec<_>>>()?;
    Ok(Layout::Regional {regions,bridges})
}

fn requirements(model:PlateMaterialModel) -> Result<Vec<Requirement>> {
    let mut r = elastic_requirements(model)?;
    for key in VISCOSITY { r.push(Requirement::try_new(key,DynViscosity::DIMS,Admit::Finite)?); }
    Ok(r)
}
fn elastic_requirements(model:PlateMaterialModel) -> Result<Vec<Requirement>> {
    let mut r = vec![Requirement::try_new("density",Density::DIMS,Admit::StrictlyPositive)?];
    let elastic: &[(&str,Dims)] = match model {
        PlateMaterialModel::Isotropic => &[("young_modulus",Pressure::DIMS),("poisson_ratio",Dims::NONE)],
        _ => &[("young_modulus_1",Pressure::DIMS),("young_modulus_2",Pressure::DIMS),
            ("poisson_ratio_12",Dims::NONE),("shear_modulus_12",Pressure::DIMS)],
    };
    for &(key,dims) in elastic { r.push(Requirement::try_new(key,dims,Admit::Finite)?); }
    Ok(r)
}
fn model(value:&str) -> Result<PlateMaterialModel> {
    match value {
        "isotropic" => Ok(PlateMaterialModel::Isotropic),
        "orthotropic12" => Ok(PlateMaterialModel::Orthotropic12),
        "orthotropic21" => Ok(PlateMaterialModel::Orthotropic21),
        _ => Err("body elastic model must be isotropic, orthotropic12 or orthotropic21".into()),
    }
}
fn scalar(value:&str) -> Result<f64> {
    let value:f64 = value.parse()?;
    if !value.is_finite() { return Err("body profile requires finite scalars".into()); }
    Ok(value)
}
struct Records<'a> { lines:std::str::Lines<'a>,line:usize }
impl<'a> Records<'a> {
    fn row(&mut self, name:&str, fields:usize) -> Result<Vec<&'a str>> {
        self.line += 1;
        let line = self.lines.next().ok_or_else(||format!("missing {name} at body line {}",self.line))?;
        // Never allocate proportional to an untrusted field count.
        let mut tokens = line.split_ascii_whitespace();
        if tokens.next()!=Some(name) { return Err(format!("expected {name} at body line {}",self.line).into()); }
        let values:Vec<_> = tokens.by_ref().take(fields).collect();
        if values.len()!=fields || tokens.next().is_some() { return Err(format!("wrong fields at body line {}",self.line).into()); }
        Ok(values)
    }
}

#[cfg(test)]
#[path="body/tests.rs"]
mod tests;


