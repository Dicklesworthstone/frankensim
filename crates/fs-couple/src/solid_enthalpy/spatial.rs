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
//! This path is one-way, stationary-reference heat to quasistatic mechanics
//! or the original continuing tetrahedral Newmark trajectory. Timed states
//! retain the explicit thermal clock across the inertial handoff.
//! Elastic storage is reported by the native stress owner, separately from
//! the native caloric audit. Mechanical work is not silently converted into
//! heat, and deformation does not silently change the transport geometry.

mod conductivity;
/// Native nodal heat and continuing pressure-driven modal motion.
pub mod motion;
mod source;

pub use conductivity::SpatialConductivitySource;

use fs_blake3::ContentHash;
use fs_conduction::bc::ThermalBoundary;
use fs_conduction::field::ScalarField;
use fs_conduction::interface::ThermalInterfaces;
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyHeatStep, EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler,
};
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
use fs_solid::loads3::transient::checkpoint::TetTransientCheckpoint3;
use fs_solid::loads3::transient::thermal::ThermalTransientStep3;
use fs_solid::loads3::transient::{TetTransientRun3, TransientConfig3, TransientError3};
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
            | SolidEnthalpyError::Loads(LoadError3::Cancelled)
            | SolidEnthalpyError::Transient(TransientError3::Cancelled) => Self::Cancelled,
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
impl From<TransientError3> for SpatialSolidEnthalpyError {
    fn from(e: TransientError3) -> Self {
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
        let conductivity = vec![SpatialConductivitySource::Isotropic; sources.len()];
        Self::try_new_with_conductivity(
            thermal,
            sources,
            conductivity,
            budget,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// Bind explicit thermal symmetry per ORIGINAL native material index.
    /// Orthotropic entries replace only the scalar conductivity role: each
    /// principal k_i(T) is pinned on the same card and complete query context,
    /// then the native heat owner uses K = R diag(k_i) R^T and its derivative.
    /// Thermal, elastic and expansion principal frames remain independent.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_conductivity(
        thermal: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
        sources: Vec<SpatialMaterialSource<'a>>,
        conductivity_sources: Vec<SpatialConductivitySource>,
        budget: TetAssemblyBudget,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SpatialSolidEnthalpyError> {
        poll(cx)?;
        if sources.len() != thermal.material_count()
            || conductivity_sources.len() != sources.len()
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
        for (index, (source, conductivity_source)) in
            sources.into_iter().zip(conductivity_sources).enumerate()
        {
            poll(cx)?;
            let curve =
                thermal
                    .reference_curve(index)
                    .ok_or(SpatialSolidEnthalpyError::Invalid(
                        "native caloric chart missing",
                    ))?;
            // A P1 profile has no scalar reference-material descriptor.
            // Admit its original chart and constitutive source independently;
            // resolve()/begin_transient() retain the actual native inventory.
            let scalar_density = thermal
                .reference_material(index)
                .map(|reference| reference.reference_density_kg_m3);
            let entry = BoundSource::try_new(
                source,
                conductivity_source,
                curve,
                scalar_density,
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
    /// Admitted frame and pinned principal selections, by original material
    /// index. The native conductivity receipts retain the exact source queries.
    pub fn conductivity_source(&self, material_index: usize) -> Option<&SpatialConductivitySource> {
        self.sources
            .get(material_index)
            .map(|s| &s.conductivity_source)
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
                stiffness_model: source.source.stiffness_model,
                expansion_model: source.source.expansion_model,
            });
        }
        poll(cx)?;
        Ok(SpatialSolidState {
            nodes_m: self.thermal.mesh().positions(),
            tetrahedra: &self.tetrahedra,
            nodal_reference_density_kg_m3: self.thermal.nodal_reference_density_kg_m3(),
            materials,
            strains,
            samples,
            budget: self.budget,
            mesh_frame: self.mesh_frame,
            time_s: None,
            maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain,
        })
    }

    /// Resolve a native nodal heat field at its explicitly supplied physical
    /// time. The heat owner has no implicit absolute clock: the caller must use
    /// the start/endpoint times of the actual admitted heat interval. Untimed
    /// [`Self::resolve`] states remain valid for equilibrium and modal queries,
    /// but cannot start or advance a full inertial trajectory.
    pub fn resolve_at(
        &self,
        specific_h: &[f64],
        time_s: f64,
        cx: &Cx<'_>,
    ) -> Result<SpatialSolidState<'_>, SpatialSolidEnthalpyError> {
        poll(cx)?;
        if !time_s.is_finite() {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "spatial thermal state requires a finite explicit physical time",
            ));
        }
        let mut state = self.resolve(specific_h, cx)?;
        state.time_s = Some(time_s);
        Ok(state)
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
        let thermal = self.advance_thermal(
            old_specific_h,
            dt_s,
            boundary,
            interfaces,
            source,
            thermal_config,
            cx,
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

    /// Native ambient or closed-enclosure radiation heats the same solid.
    ///
    /// The selected mean/pointwise operator remains owned by conduction.
    /// Its complete physical heat report survives the mechanical handoff;
    /// radiation is never converted into an extra source or convection port.
    /// Material, phase, strain or load refusal returns neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_radiation(
        &self,
        old_specific_h: &[f64],
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
        loads: &TetLoadCase3<'_, '_>,
        fixed_dofs: &[usize],
        solve: SolidEquilibriumConfig,
        cx: &Cx<'_>,
    ) -> Result<SpatialSolidCandidate<'_, EnthalpyHeatStep>, SpatialSolidEnthalpyError> {
        let thermal = self.advance_thermal_with_radiation(
            old_specific_h,
            dt_s,
            boundary,
            interfaces,
            source,
            thermal_config,
            radiation,
            cx,
        )?;
        let state = self.resolve(&thermal.conduction().specific_enthalpy_j_kg, cx)?;
        let equilibrium = state.solve(loads, fixed_dofs, solve, cx)?;
        poll(cx)?;
        Ok(SpatialSolidCandidate {
            thermal,
            state,
            equilibrium,
        })
    }

    /// Solve native nodal heat and advance the original full inertial motion
    /// over exactly the SAME supplied timestep. The start clock comes from
    /// `run`; no endpoint time or independently refreshed material is supplied.
    /// `old_specific_h` is the caller's complete accepted nodal heat history,
    /// whose projected material and eigenstrain must match the accepted run.
    ///
    /// The explicit radiation selector uses the existing heat owner, including
    /// `EnthalpyRadiation::None` for prescribed conduction boundaries alone.
    /// Its original physical energy audit is returned without turning radiation
    /// into another heat source. Mechanics keeps its independent signed material
    /// work and damping accounts; this one-way, fixed-geometry partition does
    /// not infer reciprocal thermoelastic heat or a coupled accuracy bound.
    ///
    /// The native checkpoint moves accepted motion buffers and shares its
    /// existing matrices; no self-referential problem or second trajectory is
    /// created. Bind that checkpoint to `candidate.state.problem(fixed_dofs)`
    /// to continue. Heat, source state, motion and work are returned together
    /// only after every solve, source, displacement, energy and cancellation
    /// gate succeeds. A late refusal leaves both supplied histories unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transient<'g, 'mesh>(
        &self,
        old_specific_h: &[f64],
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
        run: &TetTransientRun3<'g, 'mesh, '_, '_>,
        loads: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<SpatialSolidTransientCandidate<'_, 'g, 'mesh>, SpatialSolidEnthalpyError> {
        poll(cx)?;
        let start_s = run.state().time_s;
        let end_s = start_s + dt_s;
        if !dt_s.is_finite() || dt_s <= 0.0 || !end_s.is_finite() || end_s <= start_s {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "coupled spatial step requires a finite positive duration and a distinct endpoint",
            ));
        }
        let previous = self.resolve_at(old_specific_h, start_s, cx)?;
        if !previous.matches_transient(run) {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "accepted motion must consume the material, thermal strain and reference density of the supplied nodal heat history",
            ));
        }
        // The source and every old corner have already been admitted above.
        // Reuse the exact native heat owner without resolving them a second time.
        let thermal = self.thermal.advance_with_radiation(
            cx,
            self.thermal_problem(boundary, source)?,
            interfaces,
            old_specific_h,
            dt_s,
            thermal_config,
            radiation,
        )?;
        let state = self.resolve_at(&thermal.conduction().specific_enthalpy_j_kg, end_s, cx)?;
        let (motion, mechanical) = {
            let problem = state.problem(run.material_problem().fixed_dofs);
            // Keep the exact heat duration. In general (start + dt) - start
            // does not reproduce dt, even for an ordinary clock like 2.0 s.
            let (motion, report) =
                state.advance_transient_interval(&previous, run, &problem, loads, dt_s, cx)?;
            (motion.into_checkpoint(cx)?, report)
        };
        poll(cx)?;
        Ok(SpatialSolidTransientCandidate {
            thermal,
            state,
            motion,
            mechanical,
        })
    }

    // All consumers share the original source-bound operator and inventory.
    fn thermal_problem<'p>(
        &'p self,
        boundary: &'p ThermalBoundary,
        source: &'p ScalarField,
    ) -> Result<ConductionProblem<'p>, SpatialSolidEnthalpyError> {
        Ok(ConductionProblem {
            mesh: self.thermal.mesh(),
            boundary,
            material: self.conductivity.model_for(0)?,
            element_materials: Some(&self.conductivity),
            source,
        })
    }

    // Shared by the quasistatic and continuing-modal consumers. Neither can
    // substitute an independent conductivity or a second thermal inventory.
    #[allow(clippy::too_many_arguments)]
    fn advance_thermal(
        &self,
        old_specific_h: &[f64],
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        cx: &Cx<'_>,
    ) -> Result<EnthalpyStepSolution, SpatialSolidEnthalpyError> {
        self.resolve(old_specific_h, cx)?;
        Ok(self.thermal.advance(
            cx,
            self.thermal_problem(boundary, source)?,
            interfaces,
            old_specific_h,
            dt_s,
            thermal_config,
        )?)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_thermal_with_radiation(
        &self,
        old_specific_h: &[f64],
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<EnthalpyHeatStep, SpatialSolidEnthalpyError> {
        self.resolve(old_specific_h, cx)?;
        Ok(self.thermal.advance_with_radiation(
            cx,
            self.thermal_problem(boundary, source)?,
            interfaces,
            old_specific_h,
            dt_s,
            thermal_config,
            radiation,
        )?)
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
    pub stiffness_model: ThermalStiffnessModel,
    pub expansion_model: ThermalExpansionModel,
}

/// Source-resolved fields for native loads, stress and other solid consumers.
#[derive(Clone, Debug)]
pub struct SpatialSolidState<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    nodal_reference_density_kg_m3: Option<&'a [f64]>,
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetThermalStrainState>,
    samples: Vec<SpatialTetMaterialSample>,
    budget: TetAssemblyBudget,
    mesh_frame: ContentHash,
    time_s: Option<f64>,
    maximum_absolute_mandel_strain: f64,
}
impl SpatialSolidState<'_> {
    /// Explicit thermal time, when resolved through `resolve_at`.
    pub const fn time_s(&self) -> Option<f64> {
        self.time_s
    }
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

    fn matches_problem(&self, problem: &TetLinearElasticProblem<'_>) -> bool {
        self.nodes_m == problem.nodes_m
            && self.tetrahedra == problem.tetrahedra
            && match problem.materials {
                TetMaterialField::Uniform(material) => self.materials.iter().all(|m| m == material),
                TetMaterialField::PerElement(materials) => materials == self.materials,
            }
    }

    fn matches_transient(&self, run: &TetTransientRun3<'_, '_, '_, '_>) -> bool {
        let thermal_matches = match run.thermal_strain_field() {
            Some(TetThermalStrainField::Uniform(state)) => self.strains.iter().all(|s| s == state),
            Some(TetThermalStrainField::PerElement(states)) => states == self.strains,
            None => false,
        };
        self.time_s
            .is_some_and(|time| time.to_bits() == run.state().time_s.to_bits())
            && self.matches_problem(run.material_problem())
            && match self.nodal_reference_density_kg_m3 {
                Some(density) => run.nodal_reference_density_kg_m3() == Some(density),
                None => super::reference_inventory_matches(run, run.material_problem()),
            }
            && thermal_matches
    }

    /// Start native inertial motion at this state's explicit thermal time.
    ///
    /// Resolve through [`SpatialSolidEnthalpy::resolve_at`] and supply this
    /// state's `problem()`; it remains borrowed by the run. The heat owner's
    /// exact reference inventory supplies inertia, including any explicit P1
    /// density field. Material cards still supply stiffness and source identity.
    /// Motion, supports and allowances retain their native meanings; no static
    /// equilibration or replacement thermal inventory is introduced.
    ///
    /// `loads` already contains assembled forces. Acceleration or gravity
    /// loads must be assembled with the same reference profile through the
    /// native geometry/history APIs; this initializer does not rebuild them.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient<'g, 'm, 'p, 's>(
        &self,
        problem: &'p TetLinearElasticProblem<'s>,
        loads: &TetLoadCase3<'g, 'm>,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        config: TransientConfig3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, SpatialSolidEnthalpyError> {
        poll(cx)?;
        let Some(time_s) = self.time_s else {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "initial inertial motion requires an explicitly timed spatial state",
            ));
        };
        if !self.matches_problem(problem) {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "initial motion must use the exact spatial material problem",
            ));
        }
        let thermal = Some(self.thermal_strains());
        let run = match self.nodal_reference_density_kg_m3 {
            Some(density) => loads.begin_transient_with_nodal_reference_density(
                problem,
                density,
                thermal,
                time_s,
                displacement_m,
                velocity_m_s,
                config,
                cx,
            ),
            None => loads.begin_transient(
                problem,
                thermal,
                time_s,
                displacement_m,
                velocity_m_s,
                config,
                cx,
            ),
        }?;
        poll(cx)?;
        Ok(run)
    }

    /// Continue the original full tetrahedral motion through this native nodal
    /// heat endpoint. Both states must be resolved with `resolve_at`; the prior
    /// state's clock, material identities and thermal strain must match the
    /// accepted run exactly. Supply this endpoint's `problem()` and retain it
    /// alongside both resolved states for the lifetime of the returned run.
    ///
    /// The native owner keeps reference mass and stationary zero supports,
    /// rebuilds changed stiffness and retains displacement/velocity/acceleration.
    /// Expansion work, stiffness work and supplied Rayleigh damping retain
    /// their original separate accounts. This is prescribed one-way thermal
    /// loading: none of those accounts is silently debited from nodal enthalpy.
    /// The caller accepts the original heat result and this motion candidate
    /// together. Every refusal leaves both supplied histories unchanged.
    pub fn advance_transient_candidate<'g, 'm, 'p, 's>(
        &self,
        previous: &Self,
        run: &TetTransientRun3<'g, 'm, '_, '_>,
        problem: &'p TetLinearElasticProblem<'s>,
        loads: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), SpatialSolidEnthalpyError>
    {
        poll(cx)?;
        let (Some(start_s), Some(end_s)) = (previous.time_s, self.time_s) else {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "inertial thermal loading requires explicitly timed spatial endpoints",
            ));
        };
        self.advance_transient_interval(previous, run, problem, loads, end_s - start_s, cx)
    }

    // Shared checked handoff. The integrated consumer supplies the native heat
    // duration directly; the manual endpoint API derives it from explicit times.
    #[allow(clippy::too_many_arguments)]
    fn advance_transient_interval<'g, 'm, 'p, 's>(
        &self,
        previous: &Self,
        run: &TetTransientRun3<'g, 'm, '_, '_>,
        problem: &'p TetLinearElasticProblem<'s>,
        loads: &TetLoadCase3<'_, '_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), SpatialSolidEnthalpyError>
    {
        poll(cx)?;
        let (Some(start_s), Some(end_s)) = (previous.time_s, self.time_s) else {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "inertial thermal loading requires explicitly timed spatial endpoints",
            ));
        };
        if !previous.matches_transient(run)
            || !self.matches_problem(problem)
            || self.nodes_m != previous.nodes_m
            || self.tetrahedra != previous.tetrahedra
            || self.mesh_frame != previous.mesh_frame
            || self.nodal_reference_density_kg_m3 != previous.nodal_reference_density_kg_m3
            || self.maximum_absolute_mandel_strain.to_bits()
                != previous.maximum_absolute_mandel_strain.to_bits()
            || self.samples.len() != previous.samples.len()
        {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "spatial transient endpoints must retain exact mesh, frame, material, thermal, reference inventory and clock bindings",
            ));
        }
        for (old, next) in previous.samples.iter().zip(&self.samples) {
            poll(cx)?;
            if old.element != next.element
                || old.material_index != next.material_index
                || old.phase.phase_curve_identity() != next.phase.phase_curve_identity()
                || old.reference_elastic != next.reference_elastic
                || old.expansion.reference() != next.expansion.reference()
                || old.stiffness_model != next.stiffness_model
                || old.expansion_model != next.expansion_model
                || old.reference_density_kg_m3.to_bits() != next.reference_density_kg_m3.to_bits()
                || old.reference_volume_m3.to_bits() != next.reference_volume_m3.to_bits()
            {
                return Err(SpatialSolidEnthalpyError::Invalid(
                    "spatial heating must retain the original material assignment, caloric chart, source laws and reference inventory",
                ));
            }
        }
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || end_s <= start_s
            || (run.state().time_s + dt_s).to_bits() != end_s.to_bits()
        {
            return Err(SpatialSolidEnthalpyError::Invalid(
                "spatial heat and motion must share one distinct representable endpoint",
            ));
        }
        let candidate = run.advance_thermal_candidate(
            problem,
            self.thermal_strains(),
            loads,
            dt_s,
            self.mesh_frame,
            self.maximum_absolute_mandel_strain,
            cx,
        )?;
        poll(cx)?;
        Ok(candidate)
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
pub struct SpatialSolidCandidate<'a, Thermal = EnthalpyStepSolution> {
    pub thermal: Thermal,
    pub state: SpatialSolidState<'a>,
    pub equilibrium: SolidEquilibrium,
}

/// Jointly accepted native heat, source state and continuing full motion.
///
/// The checkpoint retains the original numerical controls, step count, loads,
/// thermal fields, reference mass and shared operators. Its existing `bind`
/// method admits only this exact material problem and the original supports.
/// Caloric storage and mechanical material work remain separate accounts.
pub struct SpatialSolidTransientCandidate<'a, 'g, 'm> {
    pub thermal: EnthalpyHeatStep,
    pub state: SpatialSolidState<'a>,
    pub motion: TetTransientCheckpoint3<'g, 'm>,
    pub mechanical: ThermalTransientStep3,
}
