//! Source/data adapter for directional body elasticity and full viscous coupling.
//! Every physical operator, tensor rotation, mode, bridge and pressure projection
//! remains in the existing plate/material/reciprocal owners.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::chart::SharedMaterialPlateChartBody;
use fs_couple::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
use fs_couple::thin_plate::{compile_plate_material_chart, PlateChartRadiation, PlateRegionMaterial};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as VISCOSITY;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};

pub const MAX_BYTES: usize = 65_536;
const MAX_GRID: usize = 12;

/// An optional complete body-model selection; the performance still owns its
/// material pack, physical dimensions, thickness, mode count and listener.
#[derive(Clone, Debug)]
pub struct Profile {
    pub bytes: Vec<u8>,
    nx: usize, ny: usize,
    model: PlateMaterialModel,
    grain_rad: f64,
    support: EdgeSupport,
    pretension_n_m: f64,
    window: (f64, f64),
    bridge: (usize, usize),
}
impl Profile {
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_BYTES { return Err("body profile exceeds 64 KiB".into()); }
        let mut input = Records { lines: std::str::from_utf8(bytes)?.lines(), line: 0 };
        input.row("frankensim-string-body-v1", 0)?;
        let v = input.row("mesh", 2)?;
        let nx: usize = v[0].parse()?; let ny: usize = v[1].parse()?;
        if !(2..=MAX_GRID).contains(&nx) || !(2..=MAX_GRID).contains(&ny) {
            return Err("body grid requires 2..=12 intervals along each axis".into());
        }
        let model = model(input.row("elastic", 1)?[0])?;
        let grain_rad = scalar(input.row("grain", 1)?[0])?;
        let support = match input.row("support", 1)?[0] {
            "simply-supported" => EdgeSupport::SimplySupported,
            "clamped" => EdgeSupport::Clamped,
            _ => return Err("body support must be simply-supported or clamped".into()),
        };
        let pretension_n_m = scalar(input.row("pretension", 1)?[0])?;
        if pretension_n_m < 0.0 { return Err("body profile requires nonnegative conservative pretension".into()); }
        let v = input.row("window", 2)?;
        let window = (scalar(v[0])?, scalar(v[1])?);
        if window.0 < 0.0 || window.1 <= window.0 { return Err("body eigenvalue window must be nonnegative and increasing".into()); }
        let v = input.row("bridge", 2)?;
        let bridge = (v[0].parse()?, v[1].parse()?);
        if bridge.0 == 0 || bridge.0 >= nx || bridge.1 == 0 || bridge.1 >= ny {
            return Err("body bridge must be an interior grid node".into());
        }
        if input.lines.next().is_some() { return Err("unconsumed body profile record".into()); }
        Ok(Self { bytes:bytes.to_vec(),nx,ny,model,grain_rad,support,pretension_n_m,window,bridge })
    }

    /// Reduce one sourced tensor on the same elastic basis as every physical
    /// bridge/radiation map; no scalar damping is added to these body modes.
    pub fn compile(&self, spec:&input::Spec, sources:Vec<MaterialStringSource>,
        catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<NormalizedPack>,
        settings:StringBankOptions, cx:&Cx<'_>, gate:&CancelGate) -> Result<ReciprocalStringBank>
    {
        poll(gate)?;
        let requirements = requirements(self.model)?;
        let material = resolve_properties(&spec.body.pack,catalog,loader,spec,&requirements)?;
        poll(gate)?;
        let mesh = PlateMesh::rectangle(spec.body.length,spec.body.width,self.nx,self.ny);
        let nodes = mesh.nodes.len(); let elements = mesh.tris.len();
        let region = PlateRegionMaterial {
            region:PlateRegion { name:"body".into(),triangle_indices:(0..elements).collect() },
            material,model:self.model,
            thickness_constraint:PlateThicknessConstraint::FixedThickness(spec.body.thickness),
            material_angle_rad:self.grain_rad,
        };
        let mut footprint = vec![0.0;nodes];
        footprint[self.bridge.1*(self.nx+1)+self.bridge.0] = 1.0;
        let specimen = compile_plate_material_chart(mesh,PlateMesh::rectangle_boundary(self.nx,self.ny),vec![region])?;
        let viscosity = PlateViscosity::from_chart(&specimen,VISCOSITY,gate)?;
        let body = SharedMaterialPlateChartBody {
            specimen,reduction:PlateChartRadiation {
                assembly:AssemblyOptions { pretension:self.pretension_n_m,support:self.support },
                eigenvalue_window:self.window,n_modes:spec.body.modes,damping_ratio:0.0,
                unit_force_weights:footprint,
            },
            bridges:spec.strings.iter().map(|w|w.bridge).collect(),
            max_plate_dofs:3*nodes,max_plate_elements:elements,modal_budget:modal_budget(),listener_m:spec.body.listener,
        };
        Ok(ReciprocalStringBank::from_viscous_material_plate_chart(cx,sources,spec.gas,body,viscosity,
            PlateViscosityLimits { max_modes:16,max_projection_terms:input::MAX_WORK as usize,
                reference_resistance_n_s_m:1.0 },spec.cavity.clone(),spec.law.clone(),settings,gate)?)
    }
}

fn requirements(model:PlateMaterialModel) -> Result<Vec<Requirement>> {
    let mut r = vec![Requirement::try_new("density",Density::DIMS,Admit::StrictlyPositive)?];
    let elastic: &[(&str,Dims)] = match model {
        PlateMaterialModel::Isotropic => &[("young_modulus",Pressure::DIMS),("poisson_ratio",Dims::NONE)],
        _ => &[("young_modulus_1",Pressure::DIMS),("young_modulus_2",Pressure::DIMS),
            ("poisson_ratio_12",Dims::NONE),("shear_modulus_12",Pressure::DIMS)],
    };
    for &(key,dims) in elastic { r.push(Requirement::try_new(key,dims,Admit::Finite)?); }
    for key in VISCOSITY { r.push(Requirement::try_new(key,DynViscosity::DIMS,Admit::Finite)?); }
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
