//! Native nodal enthalpy transport and native tetrahedral thermal stress.
//!
//! The heterogeneous backward-Euler owner retains the only heat history,
//! original caloric charts, reference masses and mesh. Conductivity, stiffness
//! and integrated expansion come from the same material cards. This consumer
//! does not copy nodal heat into a second finite-cell simulation.
//!
//! Mechanics samples each P1 element at its **enthalpy centroid**: invert the
//! arithmetic mean of its four nodal specific enthalpies, then resolve the
//! material and integrated free strain there. This is an explicit one-point
//! quadrature rung; nonlinear material gradients require mesh refinement.
//! Every corner is admitted separately before that averaging, so a cool mean
//! cannot conceal a molten or unsupported node. The native heat owner keeps
//! its own conductivity quadrature and caloric balance unchanged.
//!
//! This path is one-way, stationary-reference heat to quasistatic mechanics.
//! Elastic storage is reported by the native stress owner, separately from
//! the native caloric audit. Mechanical work is not silently converted into
//! heat, and deformation does not silently change the transport geometry.

mod source;

use fs_blake3::ContentHash;
use fs_conduction::bc::ThermalBoundary;
use fs_conduction::field::ScalarField;
use fs_conduction::interface::ThermalInterfaces;
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use fs_conduction::transient::enthalpy::{EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution};
use fs_conduction::{ConductionError, ConductionProblem};
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, PropertyKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_solid::linear3::source::SolidMaterialError;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField, TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::loads3::{LoadError3, TetLoadCase3, TetLoadGeometry3};

use super::{
    SolidEnthalpyError, SolidEquilibrium, SolidEquilibriumConfig, ThermalExpansionModel,
    ThermalExpansionSource, ThermalStiffnessModel, same_measure,
};
use source::BoundSource;

/// Source roles for one entry in the native heterogeneous material table.
/// The explicit reference enthalpy sets the stress-free material temperature;
/// it is not another evolving heat store. Fixed pressure/composition and all
/// query quantity descriptors are inherited from `point.context`.
#[derive(Clone, Debug)]
pub struct SpatialMaterialSource<'a> {
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub point: ThermalMaterialPoint<'a>,
    pub reference_specific_enthalpy_j_kg: f64,
    /// Closed declared interval with continuous support for each pinned source.
    pub temperature_interval_k: [f64; 2],
    /// Exact isotropic conductivity property key, including quantity metadata.
    pub conductivity: PropertyKey,
    pub conductivity_selection: ClaimSelection,
    pub mechanical_selection: MaterialPropertySelection,
    pub stiffness_model: ThermalStiffnessModel,
    pub expansion_selection: MaterialPropertySelection,
    pub expansion_model: ThermalExpansionModel,
}

/// Refusal publishes neither a heat endpoint nor a mechanical endpoint.
#[derive(Debug)]
pub enum SpatialSolidEnthalpyError {
    Invalid(&'static str),
    Vertex { vertex: usize, reason: &'static str },
    Support { property: String, reason: String },
    Thermal(EnthalpyError),
    Conduction(ConductionError),
    Phase(PhaseStateError),
    Solid(SolidEnthalpyError),
    Cancelled,
}
impl core::fmt::Display for SpatialSolidEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial enthalpy/solid coupling refused: {self:?}")
    }
}
impl std::error::Error for SpatialSolidEnthalpyError {}
impl From<ConductionError> for SpatialSolidEnthalpyError {
    fn from(e: ConductionError) -> Self {
        if matches!(e, ConductionError::Cancelled { .. }) {
            Self::Cancelled
        } else {
            Self::Conduction(e)
        }
    }
}
impl From<EnthalpyError> for SpatialSolidEnthalpyError {
    fn from(e: EnthalpyError) -> Self {
        if matches!(
            e,
            EnthalpyError::Conduction(ConductionError::Cancelled { .. })
        ) {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl From<PhaseStateError> for SpatialSolidEnthalpyError {
    fn from(e: PhaseStateError) -> Self {
        Self::Phase(e)
    }
}
impl From<SolidEnthalpyError> for SpatialSolidEnthalpyError {
    fn from(e: SolidEnthalpyError) -> Self {
        match e {
            SolidEnthalpyError::Cancelled
            | SolidEnthalpyError::Solid(TetElasticError::Cancelled)
            | SolidEnthalpyError::Loads(LoadError3::Cancelled) => Self::Cancelled,
            _ => Self::Solid(e),
        }
    }
}
impl From<MaterialStatePointError> for SpatialSolidEnthalpyError {
    fn from(e: MaterialStatePointError) -> Self {
        SolidEnthalpyError::from(e).into()
    }
}
impl From<ThermalMaterialError> for SpatialSolidEnthalpyError {
    fn from(e: ThermalMaterialError) -> Self {
        SolidEnthalpyError::from(e).into()
    }
}
impl From<SolidMaterialError> for SpatialSolidEnthalpyError {
    fn from(e: SolidMaterialError) -> Self {
        SolidEnthalpyError::from(e).into()
    }
}
impl From<TetElasticError> for SpatialSolidEnthalpyError {
    fn from(e: TetElasticError) -> Self {
        SolidEnthalpyError::from(e).into()
    }
}
impl From<LoadError3> for SpatialSolidEnthalpyError {
    fn from(e: LoadError3) -> Self {
        SolidEnthalpyError::from(e).into()
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), SpatialSolidEnthalpyError> {
    cx.checkpoint()
        .map_err(|_| SpatialSolidEnthalpyError::Cancelled)
}

/// A material binding to the existing heat owner. Connectivity is converted
/// once from the native mesh's u32 indices; positions and caloric curves stay
/// borrowed from that exact owner. No new transport or elasticity kernel exists.
pub struct SpatialSolidEnthalpy<'a, 'm, 'c> {
    thermal: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    tetrahedra: Vec<[usize; 4]>,
    sources: Vec<BoundSource<'a>>,
    conductivity: ElementMaterials,
    budget: TetAssemblyBudget,
    mesh_frame: ContentHash,
    maximum_absolute_mandel_strain: f64,
}
impl<'a, 'm: 'a, 'c: 'a> SpatialSolidEnthalpy<'a, 'm, 'c> {
    /// Bind one source per ORIGINAL native material index. Reference density
    /// must match that prepared inventory and the exact reference material
    /// source; source charts cannot be substituted after heat preparation.
    ///
    /// Each selected law is pinned across the declared continuous T interval,
    /// retaining all source knots. The existing integrated-expansion law uses
    /// canonical T. A separate declared phase conductivity multiplier refuses:
    /// the conductivity card here is the sole bulk-conductivity authority.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
        sources: Vec<SpatialMaterialSource<'a>>,
        budget: TetAssemblyBudget,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SpatialSolidEnthalpyError> {
        poll(cx)?;
        if sources.len() != thermal.material_count()
            || sources.is_empty()
            || u32::try_from(sources.len() - 1).is_err()
            || mesh_frame == ContentHash([0; 32])
            || !maximum_absolute_mandel_strain.is_finite()
            || maximum_absolute_mandel_strain <= 0.0
            || maximum_absolute_mandel_strain >= 1.0
        {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "native material table, representable IDs, nonzero frame and strain bound in (0,1) required",
            ));
        }
        if thermal.has_phase_conductivity() {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "source-bound conductivity does not admit a separate declared phase multiplier",
            ));
        }
        let mesh = thermal.mesh();
        if mesh.vertex_count() > budget.maximum_nodes
            || mesh.element_count() > budget.maximum_tetrahedra
        {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "solid mesh budget exceeded",
            ));
        }
        let mut tetrahedra = Vec::with_capacity(mesh.element_count());
        for tet in &mesh.complex().tets {
            poll(cx)?;
            tetrahedra.push(tet.map(|vertex| vertex as usize));
        }
        let face_entries = tetrahedra
            .len()
            .checked_mul(4)
            .ok_or(SpatialSolidEnthalpyError::Invalid("face budget overflow"))?;
        let geometry =
            TetLoadGeometry3::new(mesh.positions(), &tetrahedra, budget, face_entries, cx)?;
        for (element, &volume) in geometry.element_volumes_m3().iter().enumerate() {
            poll(cx)?;
            if !same_measure(volume, mesh.element_volume(element)) {
                return Err(SpatialSolidEnthalpyError::Invalid(
                    "heat and solid reference volume measures disagree",
                ));
            }
        }
        let mut bound = Vec::with_capacity(sources.len());
        let mut conductivities = Vec::with_capacity(sources.len());
        for (index, source) in sources.into_iter().enumerate() {
            poll(cx)?;
            let reference =
                thermal
                    .reference_material(index)
                    .ok_or(SpatialSolidEnthalpyError::Invalid(
                        "native reference material missing",
                    ))?;
            let entry = BoundSource::try_new(
                source,
                reference,
                mesh_frame,
                maximum_absolute_mandel_strain,
                cx,
            )?;
            let id = MaterialId(
                u32::try_from(index)
                    .map_err(|_| SpatialSolidEnthalpyError::Invalid("material ID overflow"))?,
            );
            conductivities.push((id, entry.conductivity.clone()));
            bound.push(entry);
        }
        let assignments = thermal
            .element_material_ids()
            .iter()
            .map(|&index| {
                u32::try_from(index)
                    .map(MaterialId)
                    .map_err(|_| SpatialSolidEnthalpyError::Invalid("material ID overflow"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let conductivity = ElementMaterials::new(MaterialTable::new(conductivities)?, assignments)?;
        poll(cx)?;
        Ok(Self {
            thermal,
            tetrahedra,
            sources: bound,
            conductivity,
            budget,
            mesh_frame,
            maximum_absolute_mandel_strain,
        })
    }

    /// The original heat owner, including its frozen nodal mass inventory.
    pub const fn thermal(&self) -> &HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
        self.thermal
    }
    /// Exact source-bound conductivity table passed into native heat assembly.
    pub const fn conductivity(&self) -> &ElementMaterials {
        &self.conductivity
    }
    /// Native load geometry on precisely the same reference mesh.
    pub fn geometry(&self, cx: &Cx<'_>) -> Result<TetLoadGeometry3<'_>, SpatialSolidEnthalpyError> {
        poll(cx)?;
        let entries = self
            .tetrahedra
            .len()
            .checked_mul(4)
            .ok_or(SpatialSolidEnthalpyError::Invalid("face budget overflow"))?;
        Ok(TetLoadGeometry3::new(
            self.thermal.mesh().positions(),
            &self.tetrahedra,
            self.budget,
            entries,
            cx,
        )?)
    }

    /// Resolve the complete spatial heat endpoint. Corner checks precede
    /// enthalpy averaging, and never replace a current density with rho(T).
    pub fn resolve(
        &self,
        specific_h: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialSolidState<'_>, SpatialSolidEnthalpyError> {
        poll(cx)?;
        if specific_h.len() != self.thermal.mesh().vertex_count() {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "one specific enthalpy per heat vertex required",
            ));
        }
        for (vertex, (&h, &id)) in specific_h
            .iter()
            .zip(self.thermal.vertex_material_ids())
            .enumerate()
        {
            poll(cx)?;
            let source = &self.sources[id];
            let phase = source.source.curve.state_at_specific_enthalpy(h)?;
            if !phase.is_fully_solid() {
                return Err(SpatialSolidEnthalpyError::Vertex {
                    vertex,
                    reason: "each corner must remain fully solid; centroid averaging cannot hide melt",
                });
            }
            if !source.contains(phase.temperature_k()) {
                return Err(SpatialSolidEnthalpyError::Vertex {
                    vertex,
                    reason: "corner leaves the continuously supported material interval",
                });
            }
            source.resolve(phase, self.mesh_frame, self.maximum_absolute_mandel_strain)?;
        }
        let mut materials = Vec::with_capacity(self.tetrahedra.len());
        let mut strains = Vec::with_capacity(self.tetrahedra.len());
        let mut samples = Vec::with_capacity(self.tetrahedra.len());
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let id = self.thermal.element_material_ids()[element];
            let source = &self.sources[id];
            // Scaling first avoids an otherwise unnecessary overflowing sum.
            let h = tet
                .iter()
                .map(|&vertex| 0.25 * specific_h[vertex])
                .sum::<f64>();
            let phase = source.source.curve.state_at_specific_enthalpy(h)?;
            let (elastic, expansion, material) =
                source.resolve(phase, self.mesh_frame, self.maximum_absolute_mandel_strain)?;
            strains.push(*expansion.thermal_strain());
            materials.push(material);
            samples.push(SpatialTetMaterialSample {
                element,
                material_index: id,
                phase,
                elastic,
                expansion,
                reference_elastic: source.reference_elastic.clone(),
                reference_density_kg_m3: source.reference_density_kg_m3,
                reference_volume_m3: self.thermal.mesh().element_volume(element),
            });
        }
        poll(cx)?;
        Ok(SpatialSolidState {
            nodes_m: self.thermal.mesh().positions(),
            tetrahedra: &self.tetrahedra,
            materials,
            strains,
            samples,
            budget: self.budget,
            mesh_frame: self.mesh_frame,
            maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain,
        })
    }

    /// Native heat and quasistatic stress/geometry accept as one candidate.
    /// Caller-supplied old enthalpy remains unchanged on every refusal. The
    /// caller supplies actual thermal boundary/contact/source owners; bulk
    /// conductivity cannot be replaced independently of the material binding.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        old_specific_h: &[f64],
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        loads: &TetLoadCase3<'_, '_>,
        fixed_dofs: &[usize],
        solve: SolidEquilibriumConfig,
        cx: &Cx<'_>,
    ) -> Result<SpatialSolidCandidate<'_>, SpatialSolidEnthalpyError> {
        self.resolve(old_specific_h, cx)?;
        let problem = ConductionProblem {
            mesh: self.thermal.mesh(),
            boundary,
            material: self.conductivity.model_for(0)?,
            element_materials: Some(&self.conductivity),
            source,
        };
        let thermal = self.thermal.advance(
            cx,
            problem,
            interfaces,
            old_specific_h,
            dt_s,
            thermal_config,
        )?;
        let state = self.resolve(&thermal.specific_enthalpy_j_kg, cx)?;
        let equilibrium = state.solve(loads, fixed_dofs, solve, cx)?;
        poll(cx)?;
        Ok(SpatialSolidCandidate {
            thermal,
            state,
            equilibrium,
        })
    }
}

/// Exact source evidence at the declared enthalpy-centroid quadrature point.
#[derive(Clone, Debug)]
pub struct SpatialTetMaterialSample {
    pub element: usize,
    pub material_index: usize,
    pub phase: EquilibriumPhaseState,
    pub elastic: ResolvedMaterialStatePoint,
    pub expansion: ThermalExpansionSource,
    pub reference_elastic: ResolvedMaterialStatePoint,
    pub reference_density_kg_m3: f64,
    pub reference_volume_m3: f64,
}

/// Source-resolved fields for native loads, stress and other solid consumers.
#[derive(Clone, Debug)]
pub struct SpatialSolidState<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetThermalStrainState>,
    samples: Vec<SpatialTetMaterialSample>,
    budget: TetAssemblyBudget,
    mesh_frame: ContentHash,
    maximum_absolute_mandel_strain: f64,
}
impl SpatialSolidState<'_> {
    pub fn materials(&self) -> &[TetElasticMaterial] {
        &self.materials
    }
    pub fn samples(&self) -> &[SpatialTetMaterialSample] {
        &self.samples
    }
    pub fn thermal_strains(&self) -> TetThermalStrainField<'_> {
        TetThermalStrainField::PerElement(&self.strains)
    }
    pub fn problem<'s>(&'s self, fixed_dofs: &'s [usize]) -> TetLinearElasticProblem<'s> {
        TetLinearElasticProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            materials: TetMaterialField::PerElement(&self.materials),
            fixed_dofs,
            budget: self.budget,
        }
    }
    pub fn solve(
        &self,
        loads: &TetLoadCase3<'_, '_>,
        fixed_dofs: &[usize],
        config: SolidEquilibriumConfig,
        cx: &Cx<'_>,
    ) -> Result<SolidEquilibrium, SpatialSolidEnthalpyError> {
        poll(cx)?;
        if config.mesh_frame != self.mesh_frame {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "equilibrium and source mesh frames must agree",
            ));
        }
        let problem = self.problem(fixed_dofs);
        let solution = loads.solve(
            &problem,
            self.thermal_strains(),
            config.static_solve,
            config.mesh_frame,
            self.maximum_absolute_mandel_strain,
            cx,
        )?;
        let deformed_mesh = problem.update_geometry_from_displacement(
            &solution.displacement,
            config.maximum_nodal_displacement_m,
            cx,
        )?;
        poll(cx)?;
        Ok(SolidEquilibrium {
            solution,
            deformed_mesh,
        })
    }
}

/// The original caloric report and native mechanical response are both present
/// only after all gates accept. They are separate energy accounts on this rung.
pub struct SpatialSolidCandidate<'a> {
    pub thermal: EnthalpyStepSolution,
    pub state: SpatialSolidState<'a>,
    pub equilibrium: SolidEquilibrium,
}
