//! Source-defined three-dimensional valves in the ordinary material DG host.
//!
//! This adapter composes the existing tensor source resolver, P1 exterior-load
//! integral, complete solid modal image and coupled aperture. It owns no
//! constitutive law, eigenproblem, damping fit or pressure/flow integrator.
//! A declaration selects complete source histories or an explicitly admitted
//! frozen-frequency image. Both retain the same geometry and conjugate ports.

use super::plate::{MaterialCache, relative_path};
use super::{
    BladeState, NormalizedMaterialCardPack, PlateValveInputError, Reader, bad, checkpoint,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::modal::{
    ModalAperture,
    relaxation::InitialModalMemory,
    solid::{HereditarySolidOptions, MaterialSolidOptions},
};
use fs_dcontact::Obstacle;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::QueryPoint;
use fs_material::gas::GasState;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint,
};
use fs_material::visco_tensor::source::{
    ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
    TensorSourceError,
};
use fs_modal::ModalDampingBudget;
use fs_solid::linear3::{
    TetAssemblyBudget,
    visco::{SolidProjectionBudget, TetViscoField, TetViscoMaterial, ViscoStateBudget},
};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::reduce::{
    OrificeSpec, ValveCardRequest, ValveElementOrder,
    visco::{
        Controls, HereditaryModalValve, HereditaryValveControls, MaterialModalValve,
        reduce_hereditary_pressure_loaded_valve, reduce_modal_pressure_loaded_valve,
    },
};
use std::sync::Arc;

mod closure;

/// One complete source selection at its original external material state.
pub struct SolidValveMaterialRecord {
    section: usize,
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
    source: ResolvedTensorMaterial,
    axes: [[f64; 3]; 3],
}
impl SolidValveMaterialRecord {
    pub const fn section(&self) -> usize {
        self.section
    }
    pub fn path(&self) -> &str {
        &self.path
    }
    pub fn pack(&self) -> &NormalizedMaterialCardPack {
        &self.pack
    }
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.source.state()
    }
    /// Same source binding used by ordinary solid and membrane dynamics.
    pub const fn source(&self) -> &ResolvedTensorMaterial {
        &self.source
    }
    /// Material-to-mesh rotation, admitted by the original tensor law.
    pub const fn material_axes(&self) -> &[[f64; 3]; 3] {
        &self.axes
    }
}

/// Original source records; the complete structural image lives in the valve.
pub struct SolidValveRecord {
    tetrahedra: usize,
    materials: Vec<SolidValveMaterialRecord>,
}
impl SolidValveRecord {
    pub const fn element_order(&self) -> &'static str {
        "tet4"
    }
    pub const fn tetrahedra(&self) -> usize {
        self.tetrahedra
    }
    pub fn materials(&self) -> &[SolidValveMaterialRecord] {
        &self.materials
    }
}

pub(super) struct PreparedSolid {
    pub image: SolidImage,
    initial: Vec<BladeState>,
    force_per_n: Vec<[f64; 3]>,
    record: SolidValveRecord,
    closure: Option<crate::bernoulli_aperture::modal::solid::SolidSurfaceClosureSpec>,
}
pub(super) enum SolidImage {
    Frozen(MaterialModalValve),
    Hereditary {
        source: HereditaryModalValve,
        memory: InitialModalMemory,
        maximum_angular_step: f64,
    },
}
impl SolidImage {
    pub fn node_count(&self) -> usize {
        match self {
            Self::Frozen(source) => source.node_count(),
            Self::Hereditary { source, .. } => source.node_count(),
        }
    }
    pub fn mode_count(&self) -> usize {
        match self {
            Self::Frozen(source) => source.modes().len(),
            Self::Hereditary { source, .. } => source.modes().len(),
        }
    }
    fn rest_opening_m(&self) -> f64 {
        match self {
            Self::Frozen(source) => source.rest_opening_m(),
            Self::Hereditary { source, .. } => source.rest_opening_m(),
        }
    }
}
impl PreparedSolid {
    pub fn into_aperture(
        self,
        gas: GasState,
        options: MaterialSolidOptions,
        contact: Obstacle,
        gate: &CancelGate,
    ) -> Result<(ModalAperture, SolidValveRecord), PlateValveInputError> {
        let force = [Some(self.force_per_n.as_slice()), None];
        let mut valve = match self.image {
            SolidImage::Frozen(source) => ModalAperture::from_material_solid(
                source,
                force,
                gas,
                options,
                self.initial,
                contact,
                gate,
            ),
            SolidImage::Hereditary {
                source,
                memory,
                maximum_angular_step,
            } => ModalAperture::from_hereditary_material_solid(
                source,
                force,
                gas,
                HereditarySolidOptions {
                    aperture: options,
                    maximum_angular_step,
                },
                self.initial,
                memory,
                contact,
                gate,
            ),
        }
        .map_err(PlateValveInputError::Physics)?;
        if let Some(closure) = self.closure {
            valve = valve.with_solid_surface_closure(closure, gate)
                .map_err(PlateValveInputError::Physics)?;
        }
        Ok((valve, self.record))
    }
}

enum Evolution {
    Frozen {
        reference_omega_rad_s: f64,
        maximum_storage_sample_error: f64,
        maximum_loss_sample_error: f64,
        maximum_modal_loss_factor: f64,
        damping: ModalDampingBudget,
    },
    Hereditary {
        memory: InitialModalMemory,
        maximum_angular_step: f64,
        projection: SolidProjectionBudget,
        maximum_factor_scalar_terms: usize,
        maximum_basis_error: f64,
    },
}
impl Evolution {
    fn read(rows: &mut Reader<'_>, maximum_modes: usize) -> Result<Self, PlateValveInputError> {
        if rows
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            == Some("solid_memory")
        {
            let mut row = rows.row("solid_memory")?;
            let memory = match row.word()? {
                "relaxed" => InitialModalMemory::Relaxed,
                "unrelaxed" => InitialModalMemory::Unrelaxed,
                _ => {
                    return Err(bad(
                        row.line,
                        "solid memory requires relaxed or unrelaxed source history",
                    ));
                }
            };
            let maximum_angular_step = row.scalar()?;
            let projection = SolidProjectionBudget {
                max_coordinates: maximum_modes,
                max_basis_scalars: row.count(1_048_576)?,
                max_material_points: row.count(2048)?,
                max_scalar_terms: row.count(67_108_864)?,
            };
            let maximum_factor_scalar_terms = row.count(67_108_864)?;
            let maximum_basis_error = row.scalar()?;
            row.finish()?;
            Ok(Self::Hereditary {
                memory,
                maximum_angular_step,
                projection,
                maximum_factor_scalar_terms,
                maximum_basis_error,
            })
        } else {
            let mut row = rows.row("solid_approximation")?;
            let selected = Self::Frozen {
                reference_omega_rad_s: row.scalar()?,
                maximum_storage_sample_error: row.scalar()?,
                maximum_loss_sample_error: row.scalar()?,
                maximum_modal_loss_factor: row.scalar()?,
                damping: ModalDampingBudget {
                    max_rel_decay_error: row.scalar()?,
                    max_freq_error_cents: row.scalar()?,
                    max_cluster_size: row.count(32)?,
                },
            };
            row.finish()?;
            Ok(selected)
        }
    }
}

fn failure(line: usize, error: impl core::fmt::Display, gate: &CancelGate) -> PlateValveInputError {
    if gate.is_requested() {
        PlateValveInputError::Physics(AcousticRealizeError::Cancelled)
    } else {
        PlateValveInputError::MaterialSource {
            line,
            source: error.to_string(),
        }
    }
}
fn vector(row: &mut super::super::Row<'_>) -> Result<[f64; 3], PlateValveInputError> {
    Ok([row.scalar()?, row.scalar()?, row.scalar()?])
}
fn property(row: &mut super::super::Row<'_>) -> Result<String, PlateValveInputError> {
    let key = row.word()?;
    if key.len() > 128 {
        return Err(bad(row.line, "solid property name exceeds 128 bytes"));
    }
    Ok(key.into())
}

/// Read a bounded P1 volume, original tensor properties, physical pressure
/// faces and force footprint. Every selected mode survives the source image.
/// Surface-pressure integration is shared with ordinary three-dimensional
/// solid mechanics; unsupported higher-order faces never become P1 loads.
pub(super) fn read(
    rows: &mut Reader<'_>,
    rest_opening_m: f64,
    cache: &mut MaterialCache,
    loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    gate: &CancelGate,
) -> Result<PreparedSolid, PlateValveInputError> {
    checkpoint(gate)?;
    let mut row = rows.row("solid")?;
    if row.word()? != "tet4" {
        return Err(bad(
            row.line,
            "native solid pressure geometry requires explicit tet4",
        ));
    }
    let source_id = row.word()?.to_string();
    if source_id.len() > 512 {
        return Err(bad(row.line, "solid source identity exceeds 512 bytes"));
    }
    row.finish()?;
    let mut row = rows.row("solid_window")?;
    let window_hz = (row.scalar()?, row.scalar()?);
    let maximum_modes = row.count(32)?;
    let retained_compliance_floor = row.scalar()?;
    row.finish()?;
    let evolution = Evolution::read(rows, maximum_modes)?;
    let mut row = rows.row("solid_limits")?;
    let assembly = TetAssemblyBudget {
        maximum_nodes: row.count(512)?,
        maximum_tetrahedra: row.count(2048)?,
        maximum_free_dofs: row.count(96)?,
        minimum_scaled_jacobian: row.scalar()?,
    };
    let internal_variables = ViscoStateBudget {
        max_internal_variables: row.count(1_048_576)?,
    };
    row.finish()?;
    let section_count = rows.count("solid_materials", 64)?;
    if section_count == 0 {
        return Err(bad(rows.line, "solid needs original material sections"));
    }
    let mut materials = Vec::with_capacity(section_count);
    let mut records = Vec::with_capacity(section_count);
    for section in 0..section_count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_material")?;
        let path = row.word()?.to_string();
        relative_path(&path, row.line)?;
        let elasticity = match row.word()? {
            "isotropic" => TensorElasticity::Isotropic,
            "orthotropic" => TensorElasticity::Orthotropic,
            _ => {
                return Err(bad(
                    row.line,
                    "solid material needs isotropic or full orthotropic elasticity",
                ));
            }
        };
        let temperature = row.scalar()?;
        let omega = row.scalar()?;
        let strain_limit = property(&mut row)?;
        let line = row.line;
        row.finish()?;
        if temperature <= 0.0 || omega < 0.0 {
            return Err(bad(
                line,
                "solid state requires positive kelvin and nonnegative angular frequency",
            ));
        }
        let mut row = rows.row("solid_axes")?;
        let axes = [vector(&mut row)?, vector(&mut row)?, vector(&mut row)?];
        row.finish()?;
        let count = rows.count("solid_branches", 8)?;
        let mut branches = Vec::with_capacity(count);
        for _ in 0..count {
            let mut row = rows.row("tensor_branch")?;
            let relaxation_time = property(&mut row)?;
            let mut stiffness_upper = core::array::from_fn(|_| String::new());
            for key in &mut stiffness_upper {
                *key = property(&mut row)?;
            }
            row.finish()?;
            branches.push(TensorBranchProperties {
                stiffness_upper,
                relaxation_time,
            });
        }
        let properties = TensorPronyProperties {
            elasticity,
            strain_limit,
            branches,
        };
        let pack = cache.load(&path, line, loader, gate)?;
        let query = QueryPoint::new()
            .with("T", temperature)
            .and_then(|p| p.with("omega", omega))
            .map_err(|e| failure(line, e, gate))?;
        let source = ResolvedTensorMaterial::resolve(
            pack.card(),
            &query,
            properties,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|e| match e {
            TensorSourceError::StatePoint(source) => PlateValveInputError::Material { line, source },
            e => failure(line, e, gate),
        })?;
        materials.push(TetViscoMaterial::from_source(&source, axes)
            .map_err(|e| failure(line, e, gate))?);
        records.push(SolidValveMaterialRecord {
            section,
            path,
            pack,
            source,
            axes,
        });
    }
    let node_count = rows.count("solid_nodes", assembly.maximum_nodes)?;
    if node_count < 4 {
        return Err(bad(
            rows.line,
            "tetrahedral solid needs at least four nodes",
        ));
    }
    let mut nodes_m = Vec::with_capacity(node_count);
    let mut fixed_dofs = Vec::new();
    let mut force_per_n = Vec::with_capacity(node_count);
    for node in 0..node_count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_node")?;
        nodes_m.push(vector(&mut row)?);
        for axis in 0..3 {
            if row.flag()? {
                fixed_dofs.push(3 * node + axis);
            }
        }
        force_per_n.push(vector(&mut row)?);
        row.finish()?;
    }
    let count = rows.count("solid_elements", assembly.maximum_tetrahedra)?;
    if count == 0 {
        return Err(bad(rows.line, "solid needs original tetrahedra"));
    }
    let mut tetrahedra = Vec::with_capacity(count);
    let mut field = Vec::with_capacity(count);
    let mut occupied = vec![false; section_count];
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_tet")?;
        tetrahedra.push([row.parse()?, row.parse()?, row.parse()?, row.parse()?]);
        let section: usize = row.parse()?;
        let material = materials
            .get(section)
            .ok_or_else(|| bad(row.line, "tetrahedron names an absent material section"))?;
        field.push(material.clone());
        occupied[section] = true;
        row.finish()?;
    }
    if occupied.iter().any(|v| !v) {
        return Err(bad(
            rows.line,
            "every declared solid material must own an original element",
        ));
    }
    let mut row = rows.row("solid_orifice")?;
    let opening_axis = vector(&mut row)?;
    let width_axis = vector(&mut row)?;
    let opposing_plane_offset_m = row.scalar()?;
    row.finish()?;
    let face_count = rows.count("solid_face_nodes", node_count)?;
    let mut face_nodes = Vec::with_capacity(face_count);
    for _ in 0..face_count {
        let mut row = rows.row("solid_face_node")?;
        let node: usize = row.parse()?;
        if node >= node_count || face_nodes.contains(&node) {
            return Err(bad(
                row.line,
                "solid opening face nodes must be unique original nodes",
            ));
        }
        face_nodes.push(node);
        row.finish()?;
    }
    let pressure_count = rows.count("solid_pressure_faces", count * 4)?;
    let mut pressure_faces = Vec::with_capacity(pressure_count);
    for _ in 0..pressure_count {
        let mut row = rows.row("solid_pressure_face")?;
        let mut nodes: [usize; 3] = [row.parse()?, row.parse()?, row.parse()?];
        nodes.sort_unstable();
        let sign = row.scalar()?;
        if sign != 1.0 && sign != -1.0 {
            return Err(bad(
                row.line,
                "solid pressure face sign must be +1 inward or -1 outward",
            ));
        }
        pressure_faces.push((nodes, sign));
        row.finish()?;
    }
    let count = rows.count("initial", 32)?;
    let mut initial = Vec::with_capacity(count);
    for _ in 0..count {
        let mut row = rows.row("state")?;
        initial.push(BladeState {
            displacement_m: row.scalar()?,
            velocity_m_s: row.scalar()?,
        });
        row.finish()?;
    }
    let line = rows.line;
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig {
        limit_bytes: Some(64 * 1024 * 1024),
        free_list_max_bytes: 0,
        ..fs_alloc::ArenaConfig::default()
    });
    let image = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x534f4c494456,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let geometry = TetLoadGeometry3::new(
            &nodes_m,
            &tetrahedra,
            assembly,
            assembly.maximum_tetrahedra * 4,
            &cx,
        )
        .map_err(|e| failure(line, e, gate))?;
        let mut ids = Vec::with_capacity(pressure_faces.len());
        for (nodes, sign) in &pressure_faces {
            checkpoint(gate)?;
            let id = geometry
                .faces()
                .iter()
                .position(|face| face.nodes == *nodes)
                .ok_or_else(|| bad(line, "pressure face must be an original exterior triangle"))?;
            ids.push(([id], *sign));
        }
        let surfaces: Vec<_> = ids
            .iter()
            .map(|(ids, sign)| SurfaceLoad3::Pressure {
                faces: ids.as_slice(),
                pascals: *sign,
            })
            .collect();
        let load = geometry
            .assemble(&surfaces, VolumeLoad3::None, &[], &cx)
            .map_err(|e| failure(line, e, gate))?;
        let mut pressure = vec![[0.0; 3]; node_count];
        for &(dof, force) in load.nodal_forces_n() {
            pressure[dof / 3][dof % 3] = force;
        }
        let request = ValveCardRequest {
            nodes_m: &nodes_m,
            tetrahedra: &tetrahedra,
            element_order: ValveElementOrder::Tet4,
            fixed_dofs: &fixed_dofs,
            material: TetViscoField::PerElement(&field),
            orifice: OrificeSpec {
                face_nodes: &face_nodes,
                opening_axis,
                width_axis,
                opposing_plane_offset_m,
            },
            window_hz,
            retained_compliance_floor,
            loss_factor: 0.0,
            source_id: &source_id,
        };
        match evolution {
            Evolution::Frozen {
                reference_omega_rad_s,
                maximum_storage_sample_error,
                maximum_loss_sample_error,
                maximum_modal_loss_factor,
                damping,
            } => reduce_modal_pressure_loaded_valve(
                &request,
                &pressure,
                &Controls {
                    reference_omega_rad_s,
                    maximum_modes,
                    assembly,
                    internal_variables,
                    maximum_storage_sample_error,
                    maximum_loss_sample_error,
                    maximum_modal_loss_factor,
                    damping,
                },
                &cx,
            )
            .map(SolidImage::Frozen)
            .map_err(|e| failure(line, e, gate)),
            Evolution::Hereditary {
                memory,
                maximum_angular_step,
                projection,
                maximum_factor_scalar_terms,
                maximum_basis_error,
            } => reduce_hereditary_pressure_loaded_valve(
                &request,
                &pressure,
                &HereditaryValveControls {
                    assembly,
                    internal_variables,
                    projection,
                    maximum_factor_scalar_terms,
                    maximum_basis_error,
                },
                &cx,
            )
            .map(|source| SolidImage::Hereditary {
                source,
                memory,
                maximum_angular_step,
            })
            .map_err(|e| failure(line, e, gate)),
        }
    })?;
    if image.mode_count() != initial.len() {
        return Err(bad(
            line,
            "initial states must cover every selected solid mode",
        ));
    }
    let gap = image.rest_opening_m();
    if !rest_opening_m.is_finite()
        || rest_opening_m <= 0.0
        || (gap - rest_opening_m).abs() > 1e-12 * gap.abs().max(rest_opening_m.abs())
    {
        return Err(bad(
            line,
            "valve rest gap must agree with the original three-dimensional opening geometry",
        ));
    }
    let closure = closure::read(rows, &image, gate)?;
    checkpoint(gate)?;
    Ok(PreparedSolid {
        image,
        initial,
        force_per_n,
        record: SolidValveRecord {
            tetrahedra: tetrahedra.len(),
            materials: records,
        },
        closure,
    })
}
