//! Material-derived cantilever slits at retained chimney endpoints.
//! Geometry and source cards enter existing DKT, exact aperture, contact and
//! Maxwell owners. This adapter neither fits a frequency nor switches a load.
use super::material::Sources;
use super::{FlueInputError, Row, Rows, bad, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use crate::bernoulli_aperture::network::{NetworkNode, boundaries::BoundaryAperture};
use crate::bernoulli_aperture::paired::BladePropagator;
use crate::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use crate::bernoulli_aperture::viscous::ViscousSlit;
use crate::flue_loop::{
    FlueError, FlueVoice,
    network::{FlueBoundaryAperture, FlueSourcePaths},
};
use crate::thin_plate::{
    PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint, compile_plate_material_chart,
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, resolve_material_state_point,
};
use fs_material::visco::GeneralizedMaxwell;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion, SliceOptions};
use fs_qty::{Density, Dims, Pressure};
use fs_scenario::assembly::{ParameterKey, ParameterPath};
use fs_scenario::entity::EntityDeclaration;
use std::sync::Arc;

/// Stable path used by the source declaration and canonical v2 force gestures.
/// Names are declaration identities; changing a display label does not retarget
/// an existing runtime. The parameter is a nonnegative closing force in newtons.
#[must_use]
pub fn boundary_closing_force_path(assembly: &str, part: &str) -> ParameterPath {
    let assembly = EntityDeclaration::assembly(assembly).identity();
    ParameterPath {
        entity: EntityDeclaration::part(assembly, part).identity(),
        parameter: ParameterKey::new("closing_force", fs_scenario::gesture::path::FORCE),
    }
}

/// Typed pressure/angle paths matching `control_source ASSEMBLY PART`.
#[must_use]
pub fn flue_source_control_paths(assembly: &str, part: &str) -> FlueSourcePaths {
    let assembly = EntityDeclaration::assembly(assembly).identity();
    FlueSourcePaths::for_entity(EntityDeclaration::part(assembly, part).identity())
}

pub(super) fn read_control_source(
    rows: &mut Rows<'_>,
) -> Result<Option<FlueSourcePaths>, FlueInputError> {
    if !rows.next_is("control_source") {
        return Ok(None);
    }
    let mut row = rows.row("control_source")?;
    let assembly = word(&mut row)?;
    let part = word(&mut row)?;
    row.end()?;
    Ok(Some(flue_source_control_paths(&assembly, &part)))
}

/// Complete normalized source and physical address of one moving boundary.
#[derive(Clone, Debug)]
pub struct FlueBoundaryRecord {
    node: usize,
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
    force: ParameterPath,
}
impl FlueBoundaryRecord {
    /// Original driven endpoint, preserved during loss-cell lowering.
    pub const fn node(&self) -> usize {
        self.node
    }
    /// Relative source artifact name.
    pub fn path(&self) -> &str {
        &self.path
    }
    /// Exact immutable original material card and its evidence.
    pub fn pack(&self) -> &NormalizedMaterialCardPack {
        &self.pack
    }
    /// Actual dimensioned consumer path for the force schedule.
    pub const fn closing_force(&self) -> &ParameterPath {
        &self.force
    }
}

struct MaxwellInput {
    properties: Vec<(String, String)>,
    band: (f64, f64),
    max_dt_over_tau: f64,
    max_angular_step: f64,
}
struct Input {
    node: usize,
    line: usize,
    path: String,
    force: ParameterPath,
    closing_force: f64,
    exterior: f64,
    temperature: f64,
    omega: f64,
    length: f64,
    width: f64,
    thickness: f64,
    nx: usize,
    ny: usize,
    rest: f64,
    initial: ApertureState,
    mode: usize,
    band: (f64, f64),
    slit_variation: f64,
    slope: f64,
    contact: (f64, f64, f64),
    slit_length: Option<f64>,
    maxwell: Option<MaxwellInput>,
}
pub(super) struct ParsedBoundaries(Vec<Input>);
fn word(row: &mut Row<'_>) -> Result<String, FlueInputError> {
    let value = row
        .fields
        .next()
        .ok_or_else(|| bad(row.line, "missing boundary declaration name"))?;
    if value.len() > 128 {
        return Err(bad(row.line, "boundary name exceeds 128 bytes"));
    }
    Ok(value.into())
}
fn physical(source: AcousticRealizeError) -> FlueInputError {
    FlueInputError::Physics(FlueError::Physics(source))
}
impl ParsedBoundaries {
    pub fn read(rows: &mut Rows<'_>) -> Result<Self, FlueInputError> {
        if !rows.next_is("boundaries") {
            return Ok(Self(Vec::new()));
        }
        let declaration = rows.record("boundaries")?;
        if declaration == "none" {
            return Ok(Self(Vec::new()));
        }
        let mut row = Row {
            fields: declaration.split_ascii_whitespace(),
            line: rows.line,
        };
        let count: usize = row.number()?;
        row.end()?;
        if !(1..=64).contains(&count) {
            return Err(bad(
                rows.line,
                "boundaries needs 1..=64 explicit moving slits",
            ));
        }
        let mut inputs = Vec::<Input>::with_capacity(count);
        for _ in 0..count {
            let mut row = rows.row("boundary")?;
            let node = row.number()?;
            let assembly = word(&mut row)?;
            let part = word(&mut row)?;
            let force = boundary_closing_force_path(&assembly, &part);
            let closing_force = row.scalar()?;
            let exterior = row.scalar()?;
            let line = row.line;
            row.end()?;
            if closing_force < 0.0
                || inputs.last().is_some_and(|previous| previous.node >= node)
                || inputs.iter().any(|previous| previous.force == force)
            {
                return Err(bad(
                    line,
                    "boundaries need distinct ordered nodes, unique physical paths and nonnegative closing force",
                ));
            }
            let mut row = rows.row("boundary_plate")?;
            let path = row
                .fields
                .next()
                .ok_or_else(|| bad(row.line, "boundary needs a material pack"))?
                .to_string();
            super::material::relative_path(&path, row.line)?;
            let temperature = row.scalar()?;
            let omega = row.scalar()?;
            let length = row.scalar()?;
            let width = row.scalar()?;
            let thickness = row.scalar()?;
            let nx: usize = row.number()?;
            let ny: usize = row.number()?;
            let rest = row.scalar()?;
            let initial = ApertureState {
                opening_m: row.scalar()?,
                opening_velocity_m_s: row.scalar()?,
            };
            row.end()?;
            if temperature <= 0.0
                || omega < 0.0
                || length <= 0.0
                || width <= 0.0
                || thickness <= 0.0
                || rest <= 0.0
                || !(1..=32).contains(&nx)
                || !(1..=32).contains(&ny)
                || (nx + 1) * (ny + 1) > 512
                || 2 * nx * ny > 2048
            {
                return Err(bad(
                    rows.line,
                    "boundary plate needs positive SI geometry/state, bounded explicit mesh and rest gap",
                ));
            }
            let mut row = rows.row("boundary_mode")?;
            let band = (row.scalar()?, row.scalar()?);
            let mode = row.number()?;
            let slit_variation = row.scalar()?;
            let slope = row.scalar()?;
            row.end()?;
            if band.0 < 0.0
                || band.1 <= band.0
                || !(core::f64::consts::TAU * band.1).powi(2).is_finite()
            {
                return Err(bad(
                    rows.line,
                    "boundary mode needs an ordered finite search band",
                ));
            }
            let mut row = rows.row("boundary_contact")?;
            let contact = (row.scalar()?, row.scalar()?, row.scalar()?);
            row.end()?;
            let slit_length = if rows.next_is("boundary_slit") {
                let record = rows.record("boundary_slit")?;
                if record == "none" {
                    None
                } else {
                    let mut row = Row {
                        fields: record.split_ascii_whitespace(),
                        line: rows.line,
                    };
                    let length = row.scalar()?;
                    row.end()?;
                    if length <= 0.0 {
                        return Err(bad(rows.line, "boundary slit length must be positive"));
                    }
                    Some(length)
                }
            } else {
                None
            };
            let mut row = rows.row("boundary_loss")?;
            let maxwell = match row.fields.next() {
                Some("elastic") => {
                    row.end()?;
                    None
                }
                Some("maxwell") => {
                    let terms: usize = row.number()?;
                    let band = (row.scalar()?, row.scalar()?);
                    let max_dt_over_tau = row.scalar()?;
                    let max_angular_step = row.scalar()?;
                    row.end()?;
                    if !(1..=64).contains(&terms) {
                        return Err(bad(rows.line, "Maxwell source needs 1..=64 named terms"));
                    }
                    let mut properties = Vec::with_capacity(terms);
                    for _ in 0..terms {
                        let mut row = rows.row("boundary_term")?;
                        properties.push((word(&mut row)?, word(&mut row)?));
                        row.end()?;
                    }
                    Some(MaxwellInput {
                        properties,
                        band,
                        max_dt_over_tau,
                        max_angular_step,
                    })
                }
                _ => {
                    return Err(bad(
                        rows.line,
                        "boundary_loss needs elastic or explicit material Maxwell terms",
                    ));
                }
            };
            inputs.push(Input {
                node,
                line,
                path,
                force,
                closing_force,
                exterior,
                temperature,
                omega,
                length,
                width,
                thickness,
                nx,
                ny,
                rest,
                initial,
                mode,
                band,
                slit_variation,
                slope,
                contact,
                slit_length,
                maxwell,
            });
        }
        Ok(Self(inputs))
    }

    pub fn bind(
        self,
        voice: FlueVoice,
        max_steps: u64,
        sources: &mut Sources,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate,
    ) -> Result<(FlueVoice, Vec<FlueBoundaryRecord>), FlueInputError> {
        let network = voice.network().expect("file construction selected a graph");
        let expected = network
            .graph()
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(node, kind)| matches!(kind, NetworkNode::Driven).then_some(node));
        if !expected.eq(self.0.iter().map(|input| input.node)) {
            return Err(bad(
                0,
                "every driven chimney requires exactly one material boundary in node order",
            ));
        }
        if self.0.is_empty() {
            return Ok((voice, Vec::new()));
        }
        let mut mechanics = Vec::with_capacity(self.0.len());
        let mut records = Vec::with_capacity(self.0.len());
        for input in self.0 {
            poll(gate)?;
            let section = network
                .graph()
                .sections
                .iter()
                .position(|section| section.nodes.contains(&input.node))
                .ok_or_else(|| bad(input.line, "boundary must meet its own represented chimney"))?;
            let gas = network.section_gases()[section];
            let impedance = network.represented_sections()[section].impedance_pa_s_m3;
            let pack = sources.load(&input.path, input.line, loader, gate)?;
            let point = QueryPoint::new()
                .with("T", input.temperature)
                .and_then(|p| p.with("omega", input.omega))
                .map_err(|source| FlueInputError::MaterialSource {
                    line: input.line,
                    source: source.to_string(),
                })?;
            let material = (|| {
                let mut requirements = vec![
                    Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
                    Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?,
                    Requirement::try_new(
                        "poisson_ratio",
                        Dims::NONE,
                        Admit::OpenInterval {
                            lower: -1.0,
                            upper: 0.5,
                        },
                    )?,
                ];
                if let Some(maxwell) = &input.maxwell {
                    for (modulus, time) in &maxwell.properties {
                        requirements.push(Requirement::try_new(
                            modulus,
                            Pressure::DIMS,
                            Admit::NonNegative,
                        )?);
                        requirements.push(Requirement::try_new(
                            time,
                            Dims([0, 0, 1, 0, 0, 0]),
                            Admit::StrictlyPositive,
                        )?);
                    }
                }
                resolve_material_state_point(
                    pack.card(),
                    &point,
                    &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                )
            })()
            .map_err(|source| FlueInputError::Material {
                line: input.line,
                source,
            })?;
            let use_band = input
                .maxwell
                .as_ref()
                .map_or(input.band, |maxwell| maxwell.band);
            if use_band.0 < 0.0 || use_band.1 <= use_band.0 {
                return Err(bad(
                    input.line,
                    "boundary material needs an ordered nonnegative use band",
                ));
            }
            for property in material.properties() {
                let band = property.angular_frequency_band_rad_s().map_err(|source| {
                    FlueInputError::Material {
                        line: input.line,
                        source,
                    }
                })?;
                if !property.is_constant_at_fixed_temperature()
                    || band.is_some_and(|(lo, hi)| {
                        lo > core::f64::consts::TAU * use_band.0
                            || hi < core::f64::consts::TAU * use_band.1
                    })
                {
                    return Err(bad(
                        input.line,
                        "boundary material parameters require constant source values applicable over the declared use band",
                    ));
                }
            }
            let mesh = PlateMesh::rectangle(input.length, input.width, input.nx, input.ny);
            // Original rectangle-grid addresses, without floating-point
            // coordinate equality or an inferred boundary on arbitrary meshes.
            let supports = (0..=input.ny).map(|j| j * (input.nx + 1)).collect();
            let slit_edges = (0..input.ny)
                .map(|j| {
                    let a = j * (input.nx + 1) + input.nx;
                    [a, a + input.nx + 1]
                })
                .collect();
            let triangles: Vec<_> = (0..mesh.tris.len()).collect();
            let relaxation = if let Some(maxwell) = input.maxwell {
                let terms = maxwell
                    .properties
                    .iter()
                    .map(|(e, tau)| {
                        (
                            material.property(e).expect("required modulus").value_si(),
                            material
                                .property(tau)
                                .expect("required relaxation time")
                                .value_si(),
                        )
                    })
                    .collect();
                let law = GeneralizedMaxwell::new(
                    material
                        .property("young_modulus")
                        .expect("required equilibrium modulus")
                        .value_si(),
                    terms,
                )
                .map_err(|source| FlueInputError::MaterialSource {
                    line: input.line,
                    source: source.to_string(),
                })?;
                Some(PlateRelaxationSpec {
                    regions: vec![PlateRelaxationRegion {
                        triangles: triangles.clone(),
                        material: law,
                        poisson_ratio: material
                            .property("poisson_ratio")
                            .expect("required Poisson ratio")
                            .value_si(),
                        band_hz: maxwell.band,
                        provenance: format!(
                            "material card {} at T={} K",
                            material.card_identity(),
                            input.temperature
                        ),
                    }],
                    max_branches: 64,
                    max_dt_over_tau: maxwell.max_dt_over_tau,
                    max_angular_step: maxwell.max_angular_step,
                })
            } else {
                None
            };
            let chart = compile_plate_material_chart(
                mesh,
                supports,
                vec![PlateRegionMaterial {
                    region: PlateRegion {
                        name: "boundary-plate".into(),
                        triangle_indices: triangles,
                    },
                    material,
                    model: PlateMaterialModel::Isotropic,
                    thickness_constraint: PlateThicknessConstraint::FixedThickness(input.thickness),
                    material_angle_rad: 0.0,
                }],
            )
            .map_err(physical)?;
            let plate = PlateApertureReduction::from_material_chart(
                chart,
                PlateApertureOptions {
                    assembly: AssemblyOptions {
                        pretension: 0.0,
                        support: EdgeSupport::Clamped,
                    },
                    eigenvalue_window: (
                        (core::f64::consts::TAU * input.band.0).powi(2),
                        (core::f64::consts::TAU * input.band.1).powi(2),
                    ),
                    mode_index: input.mode,
                    eigensolver: SliceOptions::default(),
                    max_nodes: 512,
                    max_triangles: 2048,
                    slit_edges,
                    rest_opening_m: input.rest,
                    damping_ratio: 0.0,
                    max_slit_mode_variation: input.slit_variation,
                    max_slope: input.slope,
                },
                gate,
            )
            .map_err(physical)?;
            let contact = Obstacle::new(
                vec![-1.0],
                1,
                1,
                vec![0.0],
                vec![1.0],
                input.contact.0,
                input.contact.1,
                format!("declared boundary contact {}", input.force.entity),
            )
            .and_then(|lay| lay.with_internal_loss(input.contact.2))
            .map_err(|source| FlueInputError::MaterialSource {
                line: input.line,
                source: source.to_string(),
            })?;
            let mut aperture = DynamicAperture::from_plate_with_propagator(
                plate,
                gas.density,
                impedance,
                voice.sample_period_s(),
                max_steps,
                input.initial,
                contact,
                BladePropagator::ExactZoh,
            )
            .map_err(physical)?;
            if let Some(relaxation) = relaxation {
                aperture = aperture
                    .with_plate_relaxation(relaxation, InitialApertureMemory::Relaxed)
                    .map_err(physical)?;
            }
            if let Some(length) = input.slit_length {
                aperture
                    .set_viscous_slit(ViscousSlit::from_gas(&gas, length).map_err(physical)?)
                    .map_err(physical)?;
            }
            records.push(FlueBoundaryRecord {
                node: input.node,
                path: input.path,
                pack,
                force: input.force.clone(),
            });
            mechanics.push(FlueBoundaryAperture {
                boundary: BoundaryAperture {
                    node: input.node,
                    aperture,
                },
                closing_force: input.force,
                initial_closing_force_n: input.closing_force,
                exterior_pressure_pa: input.exterior,
            });
        }
        let count = mechanics.len();
        let voice = voice
            .with_boundary_apertures(mechanics, count, gate)
            .map_err(FlueInputError::Physics)?;
        Ok((voice, records))
    }
}
