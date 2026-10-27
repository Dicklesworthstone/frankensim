//! Humidity history changes the original DKT stiffness and inertia.
//!
//! The slab coordinate is the plate's IN-PLANE x direction. Moisture is uniform
//! across width and thickness; those faces are sealed in this one-dimensional
//! reduction. Each x strip owns its local source state, with no spatial average.
//! This is small-amplitude fixed-geometry bending, not hygroscopic deformation
//! or transfer of an already vibrating state between different modal bases.

use super::*;
use crate::thin_plate::{PlateMaterialModel, plate_elasticity};
use fs_plate::{
    AssemblyOptions, PlateMesh, PlateModel, PlateSection, PlateSectionField, assemble_with_sections,
};
use fs_qty::Pressure;

/// Bounds checked before native mesh or DKT assembly allocation.
#[derive(Clone, Copy, Debug)]
pub struct MoisturePlateBudget {
    /// Maximum mesh vertices.
    pub maximum_nodes: usize,
    /// Maximum triangles and sections.
    pub maximum_triangles: usize,
}

/// Immutable geometry, supports, source selection and principal-axis assignment.
#[derive(Clone, Debug)]
pub struct MoisturePlate {
    definition: ContentHash,
    mesh: PlateMesh,
    thickness_m: f64,
    transverse_elements: usize,
    model: PlateMaterialModel,
    material_angle_rad: f64,
    elastic: ResolvedMaterialStatePoint,
    boundary: Vec<usize>,
    options: AssemblyOptions,
}

/// One complete reassembled material/mechanical snapshot.
#[derive(Clone, Debug)]
pub struct MoisturePlateState {
    time_s: f64,
    elastic: Vec<ResolvedMaterialStatePoint>,
    sections: Vec<PlateSection>,
    model: PlateModel,
    wet_mass_kg: f64,
}
impl MoisturePlateState {
    /// The native pencil consumed by ordinary modal, force and pressure owners.
    pub const fn model(&self) -> &PlateModel {
        &self.model
    }
    /// Same accepted local material state for each x strip.
    pub fn elastic_states(&self) -> &[ResolvedMaterialStatePoint] {
        &self.elastic
    }
    /// Actual sections in native triangle order.
    pub fn sections(&self) -> &[PlateSection] {
        &self.sections
    }
    /// Conserved dry mass plus transported water, exactly the modeled inertia.
    pub const fn wet_mass_kg(&self) -> f64 {
        self.wet_mass_kg
    }
    /// Accepted moisture time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
}

/// Moisture and actual mechanics are published together after all admissions.
pub struct MoisturePlateStep {
    /// New accepted moisture state; the input remains unchanged.
    pub moisture: MoistureSlab,
    /// Original water transfers and solver report.
    pub transport: MoistureStepReport,
    /// Native material-dependent DKT state.
    pub plate: MoisturePlateState,
}

impl MoisturePlate {
    /// Bind the actual slab volume to a rectangular plate on the native mesh.
    /// Width = slab cross-section / thickness. There is one x element per moisture
    /// cell and the requested number across width. No dry or wet mass is inferred
    /// from a new material density evaluation.
    ///
    /// Elasticity uses the existing isotropic/orthotropic plane-stress lowering.
    /// The reference selection is pinned, then every local accepted moisture state
    /// re-evaluates ALL selected coefficients. Missing moisture-domain data refuses.
    /// Damping is not inferred from E or from moisture; use an admitted loss owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        slab: &MoistureSlab,
        thickness_m: f64,
        transverse_elements: usize,
        model: PlateMaterialModel,
        material_angle_rad: f64,
        selection: MaterialPropertySelection,
        boundary: Vec<usize>,
        options: AssemblyOptions,
        budget: MoisturePlateBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let nx = slab.moisture.len();
        let ny = transverse_elements;
        let nodes = nx
            .checked_add(1)
            .and_then(|x| ny.checked_add(1).and_then(|y| x.checked_mul(y)))
            .ok_or(MoistureError::Invalid("plate node count overflow"))?;
        let triangles = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(2))
            .ok_or(MoistureError::Invalid("plate triangle count overflow"))?;
        if ny == 0
            || nodes > budget.maximum_nodes
            || triangles > budget.maximum_triangles
            || !positive(thickness_m)
            || !material_angle_rad.is_finite()
        {
            return Err(MoistureError::Invalid(
                "bounded nonempty plate mesh, positive thickness and finite material angle required",
            ));
        }
        let width = slab.geometry.cross_section_m2 / thickness_m;
        if !positive(width) {
            return Err(MoistureError::Invalid("plate width is unrepresentable"));
        }
        let mut requirements = Vec::new();
        match model {
            PlateMaterialModel::Isotropic => {
                requirements.push(ScalarPropertyRequirement::try_new(
                    "young_modulus",
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )?);
                requirements.push(ScalarPropertyRequirement::try_new(
                    "poisson_ratio",
                    Dims::NONE,
                    ScalarAdmissibility::OpenInterval {
                        lower: -1.0,
                        upper: 0.5,
                    },
                )?);
            }
            PlateMaterialModel::Orthotropic12 | PlateMaterialModel::Orthotropic21 => {
                for key in ["young_modulus_1", "young_modulus_2", "shear_modulus_12"] {
                    requirements.push(ScalarPropertyRequirement::try_new(
                        key,
                        Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )?);
                }
                requirements.push(ScalarPropertyRequirement::try_new(
                    "poisson_ratio_12",
                    Dims::NONE,
                    ScalarAdmissibility::Finite,
                )?);
            }
        }
        let elastic = resolve_material_state_point(
            slab.card(),
            &slab.reference.query(),
            &requirements,
            selection,
        )?;
        if elastic
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(MoistureError::Invalid(
                "moisture plate requires static elastic source laws",
            ));
        }
        let raw = PlateMesh::rectangle(slab.geometry.length_m, width, nx, ny);
        let mesh = PlateMesh::from_unstructured(raw.nodes, raw.tris)?;
        let result = Self {
            definition: slab.definition_identity(),
            mesh,
            thickness_m,
            transverse_elements,
            model,
            material_angle_rad,
            elastic,
            boundary,
            options,
        };
        result.resolve(slab, cx)?;
        poll(cx)?;
        Ok(result)
    }

    /// Actual reference mesh, shared by every returned DKT operator.
    pub const fn mesh(&self) -> &PlateMesh {
        &self.mesh
    }

    /// Rebuild the original native operators from accepted local material data.
    /// Cancellation is polled for every strip and around the bounded synchronous
    /// DKT call; it cannot preempt an individual original assembly invocation.
    pub fn resolve(
        &self,
        slab: &MoistureSlab,
        cx: &Cx<'_>,
    ) -> Result<MoisturePlateState, MoistureError> {
        poll(cx)?;
        if slab.definition_identity() != self.definition {
            return Err(MoistureError::Invalid(
                "plate requires its original moisture source and reference inventory",
            ));
        }
        let mut elastic = Vec::with_capacity(slab.moisture.len());
        let mut strip_sections = Vec::with_capacity(slab.moisture.len());
        for &m in &slab.moisture {
            poll(cx)?;
            let state = self
                .elastic
                .at_query_point(slab.card(), &slab.point_at(m)?)?;
            let (e1, e2, nu, g) = plate_elasticity(&state, self.model)?;
            let density = slab.dry_density_kg_m3 * (1.0 + m);
            strip_sections.push(PlateSection::orthotropic_plane_stress_at_angle(
                e1,
                e2,
                nu,
                g,
                self.thickness_m,
                density,
                self.material_angle_rad,
            )?);
            elastic.push(state);
        }
        let mut sections = Vec::with_capacity(self.mesh.tris.len());
        // EXACT native rectangle order: y row, x cell, then its two triangles.
        for _ in 0..self.transverse_elements {
            poll(cx)?;
            for &section in &strip_sections {
                sections.extend([section, section]);
            }
        }
        let model = assemble_with_sections(
            &self.mesh,
            PlateSectionField::PerElement(&sections),
            &self.boundary,
            &[],
            &self.options,
        )?;
        let mut modeled_mass = 0.0;
        for (triangle, section) in self.mesh.tris.iter().zip(&sections) {
            let [a, b, c] = triangle.map(|i| self.mesh.nodes[i]);
            let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0));
            modeled_mass += area * section.thickness * section.density;
        }
        let wet_mass = slab.wet_mass_kg();
        if !positive(modeled_mass)
            || !positive(wet_mass)
            || (modeled_mass - wet_mass).abs()
                > 512.0 * f64::EPSILON * wet_mass * (self.mesh.tris.len() as f64).sqrt()
        {
            return Err(MoistureError::Invalid(
                "plate inertia does not match the transported wet inventory",
            ));
        }
        poll(cx)?;
        Ok(MoisturePlateState {
            time_s: slab.time_s(),
            elastic,
            sections,
            model,
            wet_mass_kg: wet_mass,
        })
    }

    /// Propose transport AND the material-dependent mechanical endpoint.
    /// A failed mechanical source query/assembly discards the moisture candidate.
    pub fn advance(
        &self,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        config: MoistureStepConfig,
        cx: &Cx<'_>,
    ) -> Result<MoisturePlateStep, MoistureError> {
        let (moisture, transport) = slab.advance(dt_s, faces, config, cx)?;
        let plate = self.resolve(&moisture, cx)?;
        poll(cx)?;
        Ok(MoisturePlateStep {
            moisture,
            transport,
            plate,
        })
    }

    /// Accept finite-atmosphere transport and its native mechanical endpoint
    /// together. A late elastic-source/assembly refusal discards both proposed
    /// inventories. The existing fixed-geometry isothermal plate model applies.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_atmosphere(
        &self,
        atmosphere: &super::atmosphere::FiniteMoistureAtmosphere,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<super::atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        controls: super::atmosphere::AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        (super::atmosphere::AtmosphereMoistureCandidate, MoisturePlateState),
        super::atmosphere::AtmosphereMoistureError,
    > {
        let candidate = atmosphere.advance(slab, dt_s, faces, transport, controls, cx)?;
        let plate = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, plate))
    }
}
