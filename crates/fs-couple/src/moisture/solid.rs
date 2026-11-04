//! Local transported moisture drives actual three-dimensional small-strain shape.
//!
//! One reference prism or Cartesian box per moisture cell is tessellated into six native
//! tetrahedra. Each cell's selected E/nu or orthotropic tensor and three source
//! swelling derivatives are evaluated at its own moisture. The original solid
//! owner computes compatible displacement, constraint reactions, stored elastic
//! energy, deformed geometry and a signed initial-stress vibration pencil.
//!
//! This is a slow, isothermal, quasistatic response. Prescribed surface moisture
//! or the shared finite-atmosphere owner supplies its water boundary. Diffusion
//! ordinarily remains on the original reference-volume chart; deformation does
//! not silently replace its water inventory or transport metric. The layered
//! consumer's explicit spatial-diffusion path instead couples the metric through
//! compatible axial swelling. Chemical reservoir work is
//! outside this mechanical account. Sorption heat, stress-assisted diffusion,
//! hysteresis, finite deformation and dynamic mode-state remapping are separate.

use super::*;
#[path = "solid/cartesian.rs"]
mod cartesian;
pub use cartesian::MoistureCartesianSolidStep;
#[path = "solid/layered.rs"]
mod layered;
use fs_material::state_point::{
    IntegratedOrthotropicHygroscopicExpansion, integrate_orthotropic_hygroscopic_expansion,
};
use fs_solid::{
    linear3::{
        TetAssemblyBudget, TetDeformedMesh, TetEigenstrainDisplacementSolution,
        TetEigenstrainStressReport, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
        TetStaticSolveConfig,
        eigenstrain::{TetEigenstrainField, TetEigenstrainState},
        source::SolidStiffnessModel,
    },
    prestress::PrestressedTetAssembly,
};
pub use layered::{
    MoistureLayerSolidSource, MoistureLayeredSolid, MoistureLayeredSolidStep,
    MoistureSpatialAtmosphereCandidate, MoistureSpatialAtmosphereControls,
    MoistureSpatialAtmosphereReport, MoistureSpatialDiffusionReport,
};

/// Fixed source-principal and mesh frames, independent of elastic symmetry.
#[derive(Clone, Copy, Debug)]
pub struct HygroscopicAxes {
    /// Proper rotation; columns are the three source directions in mesh axes.
    pub principal_to_mesh: [[f64; 3]; 3],
    /// Nonzero identity of the authored swelling directions.
    pub principal_frame: ContentHash,
    /// Nonzero identity of the reference mesh frame.
    pub mesh_frame: ContentHash,
}

/// Native work controls plus explicit mechanical and mass admission limits.
#[derive(Clone, Copy, Debug)]
pub struct MoistureSolidControls {
    /// Allocation, free-DOF and element-quality limits.
    pub assembly: TetAssemblyBudget,
    /// Original bounded, true-residual-checked PCG solve.
    pub equilibrium: TetStaticSolveConfig,
    /// Bound on free, total and elastic mesh-frame Mandel strain.
    pub maximum_absolute_mandel_strain: f64,
    /// Maximum compatible node displacement from the ORIGINAL geometry [m].
    pub maximum_nodal_displacement_m: f64,
    /// Absolute agreement of native reference inertia and transported wet mass [kg].
    pub mass_tolerance_kg: f64,
}

/// Original geometry/source/support definition; accepted moisture is external.
#[derive(Clone, Debug)]
pub struct MoistureSolid {
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    reference_cell_volumes_m3: Vec<f64>,
    fixed_dofs: Vec<usize>,
    source: MoistureSolidSource,
    controls: MoistureSolidControls,
}

/// One pinned material/datum, shared by homogeneous and heterogeneous meshes.
#[derive(Clone, Debug)]
struct MoistureSolidSource {
    definition: ContentHash,
    stiffness: SolidStiffnessModel,
    elastic: ResolvedMaterialStatePoint,
    expansion_selection: MaterialPropertySelection,
    axes: HygroscopicAxes,
}

#[derive(Default)]
struct MoistureSolidCells {
    elastic: Vec<ResolvedMaterialStatePoint>,
    expansion: Vec<IntegratedOrthotropicHygroscopicExpansion>,
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetEigenstrainState>,
    element_mass_kg: Vec<f64>,
}

/// A complete accepted physical projection of one moisture endpoint.
#[derive(Clone, Debug)]
pub struct MoistureSolidState {
    time_s: f64,
    elastic: Vec<ResolvedMaterialStatePoint>,
    expansion: Vec<IntegratedOrthotropicHygroscopicExpansion>,
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetEigenstrainState>,
    solution: TetEigenstrainDisplacementSolution,
    stress: TetEigenstrainStressReport,
    geometry: TetDeformedMesh,
    tangent: Option<PrestressedTetAssembly>,
    current_density_kg_m3: Vec<f64>,
    wet_mass_kg: f64,
}
impl MoistureSolidState {
    /// Accepted physical moisture time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Same selected local elastic properties, in moisture-cell order.
    pub fn elastic_states(&self) -> &[ResolvedMaterialStatePoint] {
        &self.elastic
    }
    /// Actual source-integrated directional strains, in moisture-cell order.
    pub fn expansion_paths(&self) -> &[IntegratedOrthotropicHygroscopicExpansion] {
        &self.expansion
    }
    /// Native elastic tensors and conserved reference inertia, in tet order.
    pub fn materials(&self) -> &[TetElasticMaterial] {
        &self.materials
    }
    /// Actual mesh-frame initial strain, in tet order.
    pub fn eigenstrains(&self) -> &[TetEigenstrainState] {
        &self.strains
    }
    /// Original native equilibrium, including its true residual.
    pub const fn solution(&self) -> &TetEigenstrainDisplacementSolution {
        &self.solution
    }
    /// Physical local stress, reactions and stored elastic energy.
    pub const fn stress(&self) -> &TetEigenstrainStressReport {
        &self.stress
    }
    /// Compatible deformed geometry; topology and original datum remain fixed.
    pub const fn geometry(&self) -> &TetDeformedMesh {
        &self.geometry
    }
    /// Original signed reference-configuration K_e+K_g and wet M pencil.
    /// None denotes a fully prescribed body with no vibration DOFs. This is not
    /// a stability certificate; nonpositive eigenvalues must not become sound.
    pub const fn tangent(&self) -> Option<&PrestressedTetAssembly> {
        self.tangent.as_ref()
    }
    /// Current physical density = each conserved wet element mass/current volume.
    pub fn current_density_kg_m3(&self) -> &[f64] {
        &self.current_density_kg_m3
    }
    /// Original dry skeleton plus explicitly transported water [kg].
    pub const fn wet_mass_kg(&self) -> f64 {
        self.wet_mass_kg
    }
}

/// Moisture and all dependent mechanical outputs are returned together.
pub struct MoistureSolidStep {
    /// New accepted moisture; original input remains unchanged.
    pub moisture: MoistureSlab,
    /// Actual boundary water transfers and original transport solve report.
    pub transport: MoistureStepReport,
    /// Compatible source-bound mechanical state at the same endpoint.
    pub solid: MoistureSolidState,
}

impl MoistureSolidSource {
    fn bind(
        field: binding::MoistureField<'_>,
        stiffness: SolidStiffnessModel,
        elastic_selection: MaterialPropertySelection,
        expansion_selection: MaterialPropertySelection,
        axes: HygroscopicAxes,
    ) -> Result<Self, MoistureError> {
        let point = field.reference.query();
        let elastic = stiffness.resolve(field.card, &point, elastic_selection)?;
        if elastic
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(MoistureError::Invalid(
                "hygroelasticity needs static environmental elastic laws",
            ));
        }
        let expansion = integrate_orthotropic_hygroscopic_expansion(
            field.card,
            &point,
            &point,
            expansion_selection,
        )?;
        let expansion_selection = MaterialPropertySelection::PinnedByProperty(
            expansion
                .reference()
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        Ok(Self {
            definition: field.identity,
            stiffness,
            elastic,
            expansion_selection,
            axes,
        })
    }

    /// Append only original local source responses. Mechanical compatibility is
    /// solved once on the shared mesh after every participating layer succeeds.
    fn append_cells(
        &self,
        field: binding::MoistureField<'_>,
        volumes_m3: &[f64],
        strain_limit: f64,
        cells: &mut MoistureSolidCells,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        poll(cx)?;
        if field.identity != self.definition || field.moisture.len() != volumes_m3.len() {
            return Err(MoistureError::Invalid(
                "solid requires its original moisture card, datum and inventory",
            ));
        }
        cells.elastic.reserve(field.moisture.len());
        cells.expansion.reserve(field.moisture.len());
        cells.materials.reserve(6 * field.moisture.len());
        cells.strains.reserve(6 * field.moisture.len());
        cells.element_mass_kg.reserve(6 * field.moisture.len());
        for (i, (&m, &volume)) in field.moisture.iter().zip(volumes_m3).enumerate() {
            poll(cx)?;
            let point = field.point_at(m)?;
            let state = self.elastic.at_query_point(field.card, &point)?;
            let path = integrate_orthotropic_hygroscopic_expansion(
                field.card,
                &field.reference.query(),
                &point,
                self.expansion_selection.clone(),
            )?;
            let mass = field.dry_mass_per_cell_kg * (1. + m);
            let density = mass / volume;
            let mut hash =
                DomainHasher::new("org.frankensim.fs-couple.moisture-solid-local-state.v1");
            hash.update(field.identity.as_bytes());
            hash.update(state.identity().as_bytes());
            hash.update(&(i as u64).to_le_bytes());
            hash.update(&mass.to_bits().to_le_bytes());
            hash.update(&volume.to_bits().to_le_bytes());
            let material =
                self.stiffness
                    .material(&state, density, hash.finalize(), strain_limit)?;
            let strain = TetEigenstrainState::try_from_principal_strains(
                &material,
                *path.free_principal_strains(),
                self.axes.principal_to_mesh,
                self.axes.principal_frame,
                self.axes.mesh_frame,
                path.identity(),
                strain_limit,
            )?;
            for _ in 0..6 {
                cells.materials.push(material.clone());
                cells.strains.push(strain);
                cells.element_mass_kg.push(mass / 6.);
            }
            cells.elastic.push(state);
            cells.expansion.push(path);
        }
        Ok(())
    }
}

/// Common native equilibrium and inertia owner for every moisture mesh.
#[allow(clippy::too_many_arguments)]
fn solve_native_moisture_solid(
    nodes_m: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    fixed_dofs: &[usize],
    mesh_frame: ContentHash,
    controls: MoistureSolidControls,
    time_s: f64,
    wet_mass_kg: f64,
    cells: MoistureSolidCells,
    cx: &Cx<'_>,
) -> Result<MoistureSolidState, MoistureError> {
    let MoistureSolidCells {
        elastic,
        expansion,
        materials,
        strains,
        element_mass_kg,
    } = cells;
    let problem = TetLinearElasticProblem {
        nodes_m,
        tetrahedra,
        materials: TetMaterialField::PerElement(&materials),
        fixed_dofs,
        budget: controls.assembly,
    };
    let eigenstrain = TetEigenstrainField::PerElement(&strains);
    let solution = problem.solve_eigenstrain_displacement(eigenstrain, controls.equilibrium, cx)?;
    let geometry = problem.update_geometry_from_displacement(
        &solution,
        controls.maximum_nodal_displacement_m,
        cx,
    )?;
    let tangent = if fixed_dofs.len() == nodes_m.len() * 3 {
        None
    } else {
        Some(problem.assemble_eigenstrain_tangent(
            &solution,
            eigenstrain,
            mesh_frame,
            controls.maximum_absolute_mandel_strain,
            cx,
        )?)
    };
    let stress = match &tangent {
        Some(t) => t.stress().clone(),
        None => problem.recover_eigenstrain_stress(
            &solution,
            eigenstrain,
            mesh_frame,
            controls.maximum_absolute_mandel_strain,
            cx,
        )?,
    };
    let mut current_density_kg_m3 = Vec::with_capacity(element_mass_kg.len());
    for (&mass, &volume) in element_mass_kg.iter().zip(geometry.element_volumes_m3()) {
        poll(cx)?;
        let density = mass / volume;
        if !positive(density) {
            return Err(MoistureError::Invalid(
                "deformed wet density is unrepresentable",
            ));
        }
        current_density_kg_m3.push(density);
    }
    let native_mass = tangent.as_ref().map_or_else(
        || element_mass_kg.iter().sum(),
        |t| t.pencil().total_mass_kg,
    );
    if !positive(native_mass) || (native_mass - wet_mass_kg).abs() > controls.mass_tolerance_kg {
        return Err(MoistureError::Invalid(
            "native solid inertia differs from transported wet inventory",
        ));
    }
    poll(cx)?;
    Ok(MoistureSolidState {
        time_s,
        elastic,
        expansion,
        materials,
        strains,
        solution,
        stress,
        geometry,
        tangent,
        current_density_kg_m3,
        wet_mass_kg,
    })
}

impl MoistureSolid {
    /// Bind a source/geometry/support definition to the original slab inventory.
    /// Diffusion runs along x; width = slab cross-section / thickness. The other
    /// two directions are well mixed and sealed. Every tetrahedron lies wholly
    /// inside exactly one original moisture cell; no centroid averaging occurs.
    ///
    /// `elastic_selection` selects only the chosen elastic schema;
    /// `expansion_selection` selects three coefficients plus their dry-basis
    /// reference-moisture datum. Both are pinned before any evolution. Constraints
    /// are zero displacement in full node-major xyz order, with original native
    /// rigid-mode admission. No damping or chemical energy is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        slab: &MoistureSlab,
        thickness_m: f64,
        stiffness: SolidStiffnessModel,
        elastic_selection: MaterialPropertySelection,
        expansion_selection: MaterialPropertySelection,
        axes: HygroscopicAxes,
        fixed_dofs: Vec<usize>,
        controls: MoistureSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let cells = slab.moisture.len();
        let nodes = cells
            .checked_add(1)
            .and_then(|n| n.checked_mul(4))
            .ok_or(MoistureError::Invalid("moisture solid node count overflow"))?;
        let tets = cells
            .checked_mul(6)
            .ok_or(MoistureError::Invalid("moisture solid tet count overflow"))?;
        if !positive(thickness_m)
            || nodes > controls.assembly.maximum_nodes
            || tets > controls.assembly.maximum_tetrahedra
            || !positive(controls.maximum_absolute_mandel_strain)
            || !positive(controls.maximum_nodal_displacement_m)
            || !positive(controls.mass_tolerance_kg)
        {
            return Err(MoistureError::Invalid(
                "positive mechanical geometry/limits and bounded slab mesh required",
            ));
        }
        let width = slab.geometry.cross_section_m2 / thickness_m;
        if !positive(width) {
            return Err(MoistureError::Invalid("solid width is unrepresentable"));
        }
        let mut nodes_m = Vec::with_capacity(nodes);
        for i in 0..=cells {
            poll(cx)?;
            let x = slab.geometry.length_m * (i as f64 / cells as f64);
            nodes_m.extend([
                [x, 0., 0.],
                [x, width, 0.],
                [x, width, thickness_m],
                [x, 0., thickness_m],
            ]);
        }
        let mut tetrahedra = Vec::with_capacity(tets);
        let mut reference_cell_volumes_m3 = Vec::with_capacity(cells);
        for i in 0..cells {
            poll(cx)?;
            let a = 4 * i;
            let b = 4 * (i + 1);
            tetrahedra.extend([
                [a, a + 1, a + 2, b + 2],
                [a, a + 2, a + 3, b + 2],
                [a, a + 3, b + 3, b + 2],
                [a, b + 3, b, b + 2],
                [a, b, b + 1, b + 2],
                [a, b + 1, a + 1, b + 2],
            ]);
            let volume = (nodes_m[b][0] - nodes_m[a][0]) * width * thickness_m;
            if !positive(volume) {
                return Err(MoistureError::Invalid(
                    "reference solid strip volume is unrepresentable",
                ));
            }
            reference_cell_volumes_m3.push(volume);
        }
        Self::bind_mesh(
            slab.field(),
            nodes_m,
            tetrahedra,
            reference_cell_volumes_m3,
            stiffness,
            elastic_selection,
            expansion_selection,
            axes,
            fixed_dofs,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn bind_mesh(
        field: binding::MoistureField<'_>,
        nodes_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[usize; 4]>,
        reference_cell_volumes_m3: Vec<f64>,
        stiffness: SolidStiffnessModel,
        elastic_selection: MaterialPropertySelection,
        expansion_selection: MaterialPropertySelection,
        axes: HygroscopicAxes,
        fixed_dofs: Vec<usize>,
        controls: MoistureSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let source = MoistureSolidSource::bind(
            field,
            stiffness,
            elastic_selection,
            expansion_selection,
            axes,
        )?;
        let result = Self {
            nodes_m,
            tetrahedra,
            reference_cell_volumes_m3,
            fixed_dofs,
            source,
            controls,
        };
        result.resolve_field(field, cx)?;
        Ok(result)
    }
    /// Original reference vertices. Slabs have four yz corners per x plane;
    /// Cartesian fields use their complete x-fastest lattice of grid corners.
    pub fn reference_nodes_m(&self) -> &[[f64; 3]] {
        &self.nodes_m
    }
    /// Native connectivity, six consecutive tetrahedra per moisture cell.
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.tetrahedra
    }

    /// Re-evaluate local materials and solve the original native equilibrium.
    /// All outputs are candidates until source, solver, geometry and mass gates
    /// pass. The material reference and water inventory are never reinitialized.
    pub fn resolve(
        &self,
        slab: &MoistureSlab,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidState, MoistureError> {
        self.resolve_field(slab.field(), cx)
    }

    fn resolve_field(
        &self,
        field: binding::MoistureField<'_>,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidState, MoistureError> {
        let mut cells = MoistureSolidCells::default();
        self.source.append_cells(
            field,
            &self.reference_cell_volumes_m3,
            self.controls.maximum_absolute_mandel_strain,
            &mut cells,
            cx,
        )?;
        solve_native_moisture_solid(
            &self.nodes_m,
            &self.tetrahedra,
            &self.fixed_dofs,
            self.source.axes.mesh_frame,
            self.controls,
            field.time_s,
            field.wet_mass_kg(),
            cells,
            cx,
        )
    }

    /// Conservative transport followed by source-driven compatible swelling.
    /// A failed mechanical/source/cancellation gate leaves the original slab
    /// accepted; the proposed water transfer and geometry are returned together.
    pub fn advance(
        &self,
        slab: &MoistureSlab,
        target_s: f64,
        faces: &[MoistureFace; 2],
        transport: MoistureStepConfig,
        cx: &Cx<'_>,
    ) -> Result<MoistureSolidStep, MoistureError> {
        let (moisture, report) = slab.advance(target_s, faces, transport, cx)?;
        let solid = self.resolve(&moisture, cx)?;
        poll(cx)?;
        Ok(MoistureSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }

    /// Accept finite-atmosphere transport and compatible native solid mechanics
    /// together. Source-integrated swelling, displacement, stress, geometry and
    /// the prestressed wet-inertia pencil all resolve at the same endpoint.
    /// A late source, equilibrium, geometry or cancellation refusal discards
    /// both proposed water inventories. The original fixed-reference diffusion
    /// and small-strain isothermal mechanics remain the physical owners; no
    /// sorption heat, stress-driven diffusion or dynamic state remap is inferred.
    /// Gas internal-energy change is not total sorption heat, and the changed
    /// gas pressure is not automatically applied as a mechanical surface traction.
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
        (
            super::atmosphere::AtmosphereMoistureCandidate,
            MoistureSolidState,
        ),
        super::atmosphere::AtmosphereMoistureError,
    > {
        let candidate = atmosphere.advance(slab, dt_s, faces, transport, controls, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }

    /// Couple finite solution activity, original bound-water transport and
    /// compatible native swelling as one unpublished endpoint. Nonvolatile
    /// solutes remain in the equal-pressure solution; the original sorption
    /// law is consumed through its equilibrium vapor-activity convention.
    /// A late mechanical/source refusal leaves both water inventories accepted
    /// at their originals. Osmotic pressure is not a substituted surface load.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_solution(
        &self,
        reservoir: &super::solution::FiniteMoistureSolution,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<super::atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        controls: super::solution::SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        (
            super::solution::SolutionMoistureCandidate,
            MoistureSolidState,
        ),
        super::solution::SolutionMoistureError,
    > {
        let candidate = reservoir.advance(slab, dt_s, faces, transport, controls, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }
}
