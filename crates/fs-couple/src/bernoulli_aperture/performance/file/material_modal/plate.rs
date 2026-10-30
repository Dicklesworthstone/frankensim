//! Original regional material sources and meshes for the common modal valve.
use super::{PlateValveInputError, Reader, bad, checkpoint};
use crate::bernoulli_aperture::modal::material::ModalPlateOptions;
use crate::bernoulli_aperture::modal::relaxation::{InitialModalMemory, TensorModalRelaxationLimits};
use crate::bernoulli_aperture::modal::ModalAperture;
use crate::bernoulli_aperture::paired::BladeState;
use crate::thin_plate::{PlateChartRadiation, PlateMaterialModel, PlateRegionMaterial,
    PlateThicknessConstraint, ResolvedPlateChart, compile_plate_material_chart};
use crate::plate_viscosity::PlateViscosityLimits;
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, resolve_material_state_point};
use fs_material::visco_tensor::source::{TensorBranchProperties, TensorElasticity, TensorPronyProperties};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use std::{collections::BTreeMap, sync::Arc};

/// Maximum encoded normalized artifact; the CLI bounds bytes before decoding.
pub const MAX_MODAL_VALVE_PACK_BYTES: usize = 8 * 1024 * 1024;
const MAX_TOTAL_PACK_BYTES: usize = 64 * 1024 * 1024;

/// One material assignment on an actual plate, retaining the complete artifact.
#[derive(Clone, Debug)]
pub struct ModalValveMaterialRecord {
    blade: usize,
    section: usize,
    region: usize,
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
}
impl ModalValveMaterialRecord {
    /// Original independent blade address.
    pub const fn blade(&self) -> usize { self.blade }
    /// Original authored material section address within that blade.
    pub const fn section(&self) -> usize { self.section }
    /// Resolved chart region address. A thickness map can subdivide an authored
    /// section; every subdivision retains that same section, pack and state.
    pub const fn region(&self) -> usize { self.region }
    /// Relative artifact name, interpreted only by the caller.
    pub fn path(&self) -> &str { &self.path }
    /// Original claims, material state and source provenance.
    pub fn pack(&self) -> &NormalizedMaterialCardPack { &self.pack }
}

pub(super) enum LossSelection {
    Elastic,
    Tensor { properties: [String; 6], limits: PlateViscosityLimits },
    Maxwell { properties: Vec<(String, String)>, initial: InitialModalMemory, regional: bool },
    TensorMaxwell { strain_limit: String, branches: Vec<TensorBranchProperties>,
        initial: InitialModalMemory, limits: TensorModalRelaxationLimits },
}
impl LossSelection {
    pub fn read(rows: &mut Reader<'_>) -> Result<Self, PlateValveInputError> {
        let mut row = rows.row("loss")?;
        let selected = match row.word()? {
            "elastic" => Self::Elastic,
            "tensor" => {
                let mut properties: [String; 6] = core::array::from_fn(|_| String::new());
                for key in &mut properties {
                    let value = row.word()?;
                    if value.len() > 128 { return Err(bad(row.line, "viscosity property name exceeds 128 bytes")); }
                    *key = value.into();
                }
                let reference_resistance_n_s_m = row.scalar()?;
                let max_projection_terms = row.count(16_777_216)?;
                Self::Tensor { properties, limits: PlateViscosityLimits {
                    max_modes: 32, max_projection_terms, reference_resistance_n_s_m,
                } }
            }
            selection @ ("proportional-maxwell" | "regional-maxwell") => {
                let regional = selection == "regional-maxwell";
                let initial = match row.word()? {
                    "relaxed" => InitialModalMemory::Relaxed,
                    "unrelaxed" => InitialModalMemory::Unrelaxed,
                    _ => return Err(bad(row.line, "Maxwell initial memory must be relaxed or unrelaxed")),
                };
                let count = row.count(8)?;
                if count == 0 { return Err(bad(row.line, "Maxwell material needs 1..=8 explicit source terms")); }
                let mut properties = Vec::with_capacity(count);
                for _ in 0..count {
                    let modulus = row.word()?.to_string();
                    let time = row.word()?.to_string();
                    if modulus.len() > 128 || time.len() > 128 {
                        return Err(bad(row.line, "Maxwell property name exceeds 128 bytes"));
                    }
                    properties.push((modulus, time));
                }
                Self::Maxwell { properties, initial, regional }
            }
            "tensor-maxwell" => {
                let initial = match row.word()? {
                    "relaxed" => InitialModalMemory::Relaxed,
                    "unrelaxed" => InitialModalMemory::Unrelaxed,
                    _ => return Err(bad(row.line, "tensor Maxwell initial memory must be relaxed or unrelaxed")),
                };
                let strain_limit = source_key(&mut row)?;
                let count = row.count(8)?;
                if count == 0 { return Err(bad(row.line, "tensor Maxwell requires 1..=8 complete source branches")); }
                let limits = TensorModalRelaxationLimits {
                    max_material_points: row.count(4096)?,
                    max_memory_states: row.count(256)?,
                    max_projection_terms: row.count(67_108_864)?,
                };
                row.finish()?;
                let mut branches = Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row = rows.row("tensor_branch")?;
                    let relaxation_time = source_key(&mut row)?;
                    let mut stiffness_upper = core::array::from_fn(|_| String::new());
                    for key in &mut stiffness_upper { *key = source_key(&mut row)?; }
                    row.finish()?;
                    branches.push(TensorBranchProperties { stiffness_upper, relaxation_time });
                }
                return Ok(Self::TensorMaxwell { strain_limit, branches, initial, limits });
            }
            _ => return Err(bad(row.line, "loss requires elastic, tensor viscosity, proportional-maxwell, regional-maxwell or tensor-maxwell source properties")),
        };
        row.finish()?;
        Ok(selected)
    }

    fn tensor_properties(&self, model: PlateMaterialModel) -> Option<TensorPronyProperties> {
        let Self::TensorMaxwell { strain_limit, branches, .. } = self else { return None; };
        Some(TensorPronyProperties {
            elasticity: match model {
                PlateMaterialModel::Isotropic => TensorElasticity::Isotropic,
                _ => TensorElasticity::Orthotropic,
            },
            strain_limit: strain_limit.clone(), branches: branches.clone(),
        })
    }

    /// Both native duct routes attach the same material law to the same valve.
    pub(super) fn attach(self, valve: ModalAperture, gate: &CancelGate)
        -> Result<ModalAperture, PlateValveInputError>
    {
        match &self {
            Self::Maxwell { properties, initial, regional } => {
                if *regional { valve.with_regional_material_relaxation(properties, *initial, gate) }
                else { valve.with_material_relaxation(properties, *initial, gate) }
                    .map_err(PlateValveInputError::Physics)
            }
            Self::TensorMaxwell { initial, limits, .. } => {
                let mut laws = Vec::with_capacity(valve.material_plates().len());
                for plate in valve.material_plates() {
                    let mut regional = Vec::with_capacity(plate.source().regions().len());
                    for region in plate.source().regions() {
                        checkpoint(gate)?;
                        let properties = self.tensor_properties(region.input.model)
                            .expect("tensor Maxwell selection");
                        regional.push(properties.bind(&region.input.material).map_err(|error|
                            PlateValveInputError::MaterialSource { line: 0, source: error.to_string() })?);
                    }
                    laws.push(regional);
                }
                valve.with_tensor_material_relaxation(laws, *initial, *limits, gate)
                    .map_err(PlateValveInputError::Physics)
            }
            _ => Ok(valve),
        }
    }
}

fn source_key(row: &mut super::super::Row<'_>) -> Result<String, PlateValveInputError> {
    let key = row.word()?;
    if key.len() > 128 { return Err(bad(row.line, "tensor source property name exceeds 128 bytes")); }
    Ok(key.to_string())
}

#[derive(Default)]
pub(super) struct MaterialCache {
    packs: BTreeMap<String, Arc<NormalizedMaterialCardPack>>,
    bytes: usize,
}
impl MaterialCache {
    pub(super) fn load(&mut self, path: &str, line: usize,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate) -> Result<Arc<NormalizedMaterialCardPack>, PlateValveInputError>
    {
        checkpoint(gate)?;
        if !self.packs.contains_key(path) {
            let pack = loader(path).map_err(|source| PlateValveInputError::MaterialSource { line, source })?;
            checkpoint(gate)?;
            let size = pack.to_bytes().len();
            if size > MAX_MODAL_VALVE_PACK_BYTES || size > MAX_TOTAL_PACK_BYTES.saturating_sub(self.bytes) {
                return Err(bad(line, "material artifacts exceed 8 MiB per pack or 64 MiB combined"));
            }
            self.bytes += size;
            self.packs.insert(path.into(), Arc::new(pack));
        }
        Ok(Arc::clone(&self.packs[path]))
    }
}

pub(super) fn relative_path(path: &str, line: usize) -> Result<(), PlateValveInputError> {
    if path.len() > 512 || path.contains('\\') || path.contains(':')
        || path.split('/').any(|p| p.is_empty() || p == "." || p == "..") {
        return Err(bad(line, "material pack path must stay relative to the performance directory"));
    }
    Ok(())
}

type PlateInput = (ResolvedPlateChart, ModalPlateOptions, Vec<BladeState>, Vec<ModalValveMaterialRecord>);

pub(super) fn read(
    rows: &mut Reader<'_>, blade: usize, existing_modes: usize, loss: &LossSelection,
    cache: &mut MaterialCache,
    loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    gate: &CancelGate,
) -> Result<PlateInput, PlateValveInputError> {
    let mut row = rows.row("blade")?;
    let n_modes = row.count(32usize.saturating_sub(existing_modes))?;
    let low = row.scalar()?;
    let high = row.scalar()?;
    let maximum_angular_step = row.scalar()?;
    let maximum_trace_slope = row.scalar()?;
    let maximum_slit_deviation_m = row.scalar()?;
    let max_projection_terms = row.count(16_777_216)?;
    row.finish()?;
    if n_modes == 0 || low < 0.0 || high <= low
        || !(core::f64::consts::TAU * high).powi(2).is_finite() {
        return Err(bad(rows.line, "each blade requires positive explicit modes and an ordered finite frequency window"));
    }
    let mut row = rows.row("support")?;
    let support = match row.word()? {
        "clamped" => EdgeSupport::Clamped,
        "simply-supported" => EdgeSupport::SimplySupported,
        _ => return Err(bad(row.line, "unknown physical plate support")),
    };
    let pretension = row.scalar()?;
    row.finish()?;
    let section_count = rows.count("sections", 64)?;
    if section_count == 0 { return Err(bad(rows.line, "every blade needs material sections")); }
    let mut sections = Vec::with_capacity(section_count);
    let mut records = Vec::with_capacity(section_count);
    for section in 0..section_count {
        let mut row = rows.row("section")?;
        let path = row.word()?.to_string();
        relative_path(&path, row.line)?;
        let model = match row.word()? {
            "isotropic" => PlateMaterialModel::Isotropic,
            "orthotropic12" => PlateMaterialModel::Orthotropic12,
            "orthotropic21" => PlateMaterialModel::Orthotropic21,
            _ => return Err(bad(row.line, "material section needs an explicit isotropic or ordered orthotropic plane")),
        };
        let temperature = row.scalar()?;
        let omega = row.scalar()?;
        let thickness = row.scalar()?;
        let angle = row.scalar()?;
        let line = row.line;
        let environment = crate::material_source::read_fixed_environment(&mut row.fields)
            .map_err(|source| PlateValveInputError::MaterialSource { line, source: source.to_string() })?;
        row.finish()?;
        if temperature <= 0.0 || omega < 0.0 || thickness <= 0.0 {
            return Err(bad(line, "material state needs positive kelvin/thickness and nonnegative angular frequency"));
        }
        let pack = cache.load(&path, line, loader, gate)?;
        let point = fs_material::state_point::query_at_fixed_environment(
            pack.card(), temperature, omega, &environment)
            .map_err(|e| PlateValveInputError::MaterialSource { line, source: e.to_string() })?;
        let tensor_requirements = loss.tensor_properties(model).map(|p| p.requirements())
            .transpose().map_err(|error| PlateValveInputError::MaterialSource {
                line, source: error.to_string(),
            })?;
        let material = (|| {
            let mut required = vec![Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?];
            match model {
                PlateMaterialModel::Isotropic => {
                    required.push(Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?);
                    required.push(Requirement::try_new("poisson_ratio", Dims::NONE,
                        Admit::OpenInterval { lower: -1.0, upper: 0.5 })?);
                }
                _ => {
                    for key in ["young_modulus_1", "young_modulus_2", "shear_modulus_12"] {
                        required.push(Requirement::try_new(key, Pressure::DIMS, Admit::StrictlyPositive)?);
                    }
                    required.push(Requirement::try_new("poisson_ratio_12", Dims::NONE, Admit::Finite)?);
                }
            }
            if let LossSelection::Tensor { properties, .. } = loss {
                for key in properties {
                    required.push(Requirement::try_new(key, DynViscosity::DIMS, Admit::Finite)?);
                }
            }
            if let LossSelection::Maxwell { properties, .. } = loss {
                for (modulus, time) in properties {
                    // Several arms may intentionally share a source time or
                    // modulus. Resolve each claim once; the law still checks
                    // every use's dimensions and retains each separate arm.
                    if !required.iter().any(|r| r.name() == modulus) {
                        required.push(Requirement::try_new(modulus, Pressure::DIMS, Admit::NonNegative)?);
                    }
                    if !required.iter().any(|r| r.name() == time) {
                        required.push(Requirement::try_new(time, Dims([0, 0, 1, 0, 0, 0]), Admit::StrictlyPositive)?);
                    }
                }
            }
            if let Some(requirements) = &tensor_requirements {
                for requirement in requirements {
                    if let Some(existing) = required.iter().find(|r| r.name() == requirement.name()) {
                        if existing.quantity() != requirement.quantity() {
                            return Err(fs_material::state_point::MaterialStatePointError::InvalidRequirement {
                                property: requirement.name().into(), reason: "tensor source property dimensions conflict with the structural chart",
                            });
                        }
                    } else { required.push(requirement.clone()); }
                }
            }
            resolve_material_state_point(pack.card(), &point, &required, MaterialPropertySelection::SingleClaimOnly)
        })().map_err(|source| PlateValveInputError::Material { line, source })?;
        sections.push(PlateRegionMaterial {
            region: PlateRegion { name: format!("section-{section}"), triangle_indices: Vec::new() },
            material, model, thickness_constraint: PlateThicknessConstraint::FixedThickness(thickness),
            material_angle_rad: angle,
        });
        records.push(ModalValveMaterialRecord { blade, section, region: section, path, pack });
    }
    let node_count = rows.count("nodes", 512)?;
    if node_count < 3 { return Err(bad(rows.line, "blade needs at least three original nodes")); }
    let mut nodes = Vec::with_capacity(node_count);
    let mut supports = Vec::new();
    let mut unit_force_weights = Vec::with_capacity(node_count);
    for node in 0..node_count {
        checkpoint(gate)?;
        let mut row = rows.row("node")?;
        nodes.push((row.scalar()?, row.scalar()?));
        if row.flag()? { supports.push(node); }
        unit_force_weights.push(row.scalar()?);
        row.finish()?;
    }
    let triangle_count = rows.count("triangles", 2048)?;
    if triangle_count == 0 { return Err(bad(rows.line, "blade needs original triangles")); }
    let mut triangles = Vec::with_capacity(triangle_count);
    let mut thicknesses = Vec::with_capacity(triangle_count);
    for element in 0..triangle_count {
        checkpoint(gate)?;
        let mut row = rows.row("triangle")?;
        triangles.push([row.parse()?, row.parse()?, row.parse()?]);
        let section: usize = row.parse()?;
        let thickness = if let Some(selector) = row.fields.next() {
            if selector != "thickness" {
                return Err(bad(row.line, "triangle geometry accepts only thickness THICKNESS_M"));
            }
            let value = row.scalar()?;
            if value <= 0.0 {
                return Err(bad(row.line, "triangle thickness must be finite and positive"));
            }
            Some(value)
        } else { None };
        row.finish()?;
        sections.get_mut(section).ok_or_else(|| bad(rows.line, "triangle references an absent material section"))?
            .region.triangle_indices.push(element);
        thicknesses.push(thickness);
    }
    let (sections, records) = thickness_regions(sections, records, &thicknesses, gate)?;
    let edge_count = rows.count("slit_edges", node_count)?;
    let mut slit_edges = Vec::with_capacity(edge_count);
    for _ in 0..edge_count {
        let mut row = rows.row("edge")?;
        slit_edges.push([row.parse()?, row.parse()?]);
        row.finish()?;
    }
    let mut initial = Vec::with_capacity(n_modes);
    if rows.count("initial", 32)? != n_modes {
        return Err(bad(rows.line, "initial states must cover exactly the requested blade modes"));
    }
    for _ in 0..n_modes {
        let mut row = rows.row("state")?;
        initial.push(BladeState { displacement_m: row.scalar()?, velocity_m_s: row.scalar()? });
        row.finish()?;
    }
    let mesh = PlateMesh::from_unstructured(nodes, triangles).map_err(PlateValveInputError::Geometry)?;
    let chart = compile_plate_material_chart(mesh, supports, sections).map_err(PlateValveInputError::Physics)?;
    let options = ModalPlateOptions {
        reduction: PlateChartRadiation {
            assembly: AssemblyOptions { pretension, support },
            eigenvalue_window: ((core::f64::consts::TAU * low).powi(2), (core::f64::consts::TAU * high).powi(2)),
            n_modes, damping_ratio: 0.0, unit_force_weights,
        },
        slit_edges, maximum_trace_slope, maximum_slit_deviation_m, maximum_angular_step,
        max_nodes: 512, max_triangles: 2048, max_projection_terms,
    };
    Ok((chart, options, initial, records))
}


/// Lower authored element geometry through the same regional material chart.
/// There are at most as many occupied regions as original triangles; pack
/// queries and source-byte budgets remain those of the original declarations.
/// Equal thicknesses share a region. The original uniform layout is retained,
/// including when a file explicitly repeats its section's thickness.
fn thickness_regions(
    sections: Vec<PlateRegionMaterial>,
    records: Vec<ModalValveMaterialRecord>,
    thicknesses: &[Option<f64>],
    gate: &CancelGate,
) -> Result<(Vec<PlateRegionMaterial>, Vec<ModalValveMaterialRecord>), PlateValveInputError> {
    let mut regions = Vec::new();
    let mut sources = Vec::new();
    for (section, record) in sections.into_iter().zip(records) {
        checkpoint(gate)?;
        let PlateThicknessConstraint::FixedThickness(default) = section.thickness_constraint else {
            return Err(bad(0, "blade thickness maps require an explicit section thickness"));
        };
        let mut groups = Vec::<(f64, Vec<usize>)>::new();
        let mut indices = BTreeMap::new();
        for &element in &section.region.triangle_indices {
            checkpoint(gate)?;
            let thickness = thicknesses[element].unwrap_or(default);
            let index = *indices.entry(thickness.to_bits()).or_insert_with(|| {
                groups.push((thickness, Vec::new()));
                groups.len() - 1
            });
            groups[index].1.push(element);
        }
        // Keep empty declarations too: the existing chart owner must still
        // refuse an unused physical section instead of silently dropping it.
        if groups.is_empty() || (groups.len() == 1 && groups[0].0.to_bits() == default.to_bits()) {
            let mut record = record;
            record.region = regions.len();
            regions.push(section);
            sources.push(record);
        } else {
            for (index, (thickness, triangles)) in groups.into_iter().enumerate() {
                let mut region = section.clone();
                region.region.name = format!("{}-thickness-{index}", section.region.name);
                region.region.triangle_indices = triangles;
                region.thickness_constraint = PlateThicknessConstraint::FixedThickness(thickness);
                let mut source = record.clone();
                source.region = regions.len();
                regions.push(region);
                sources.push(source);
            }
        }
    }
    Ok((regions, sources))
}
