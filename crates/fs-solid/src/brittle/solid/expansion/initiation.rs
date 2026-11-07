//! Thermal elastic prehistory and source-bound tensile-plane discovery.
//!
//! The native P1 solver evaluates the original elastic backbone at TOTAL strain
//! minus source-integrated free expansion. Free expansion cannot initiate a
//! crack; restraint and nonuniform expansion can. As in the original thermal
//! brittle owner, coefficients change first at held MECHANICAL strain, followed
//! by a straight mechanical-strain increment under the next coefficients.
//! The native Rankine detector inspects that chord, not a continuous temperature
//! trajectory. Its fractions are neither Kelvin nor physical event times.
//!
//! A coefficient change that activates the held mechanical state has its own
//! result before any displacement solve. It selects a plane at the explicitly
//! supplied endpoint without inventing a strain chord or physical event time.
//! Source extrema, reversals and nonlinear temperature paths need physical
//! endpoints. This owner supplies no heat equation, caloric heat or entropy.
//!
//! A crossing leaves the accepted elastic state untouched. Original and next
//! sources are prepared with the same discovered planes and geometric widths.
//! Handoff transfers the OLD accepted thermal state, with no new work or loss;
//! the existing thermal brittle owner replays the requested endpoint. Selected
//! planes then stay fixed. Repeated principal tension requires an explicit
//! physical direction hint; later predicted events ignore fracture-induced
//! redistribution, as in the nonthermal tensile predictor.

/// Retained thermal midpoint motion through either tensile initiation stage.
pub mod dynamics;

use std::sync::Arc;

use fs_evidence::ModelCard;
use fs_exec::Cx;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::source::initiation::{
    ResolvedTensileActivation, ResolvedTensileInitiationError, ResolvedTensileOnset,
};
use fs_material::tensor::contract;
use fs_material::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};

use super::{
    ThermalBrittleTetEnergy, ThermalBrittleTetError, ThermalBrittleTetSolid, ThermalBrittleTetStep,
};
use crate::brittle::TetCrackBandGeometry;
use crate::brittle::solid::{BrittleTetControls, BrittleTetError, checkpoint, energy_gate};
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};

/// Accepted tensile-only elastic history with an immutable expansion reference.
/// Total FE strain is retained; mechanical strain is always recovered from its
/// exact original thermal source. No damaged or plastic history can be imported.
#[derive(Debug)]
pub struct ThermalTensileElasticTetSolid<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedThermalBrittle],
    budget: TetAssemblyBudget,
    solution: SmallStrainTetSolution<()>,
    stored_energy_j: f64,
    supports: Option<Vec<usize>>,
}

/// Accepted elastic endpoint; physical fracture accounts are exactly zero.
#[derive(Debug)]
pub struct ThermalTensileElasticTetStep<'a> {
    pub solid: ThermalTensileElasticTetSolid<'a>,
    pub energy: ThermalBrittleTetEnergy,
}

/// Accepted elasticity, a strain-chord crossing, or coefficient activation.
#[derive(Debug)]
pub enum ThermalTensileTetProbe<'a> {
    Elastic(ThermalTensileElasticTetStep<'a>),
    Onset(ThermalTensileTetOnset<'a>),
    /// The new coefficients activate the OLD mechanical state before solving
    /// displacement. No candidate strain chord or event fraction is implied.
    Activation(ThermalTensileTetActivation<'a>),
}

/// Original accepted state and candidate thermal endpoint, both source-bound.
#[derive(Debug)]
pub struct ThermalTensileTetOnset<'a> {
    request: ThermalTensileTetRequest<'a>,
    predictor: SmallStrainTetSolution<()>,
    events: Vec<Option<ResolvedTensileOnset<'a>>>,
    first_fraction: [f64; 2],
    first_candidate: usize,
}

/// Tensile activation caused by an explicitly supplied coefficient endpoint.
/// Each observation uses the held OLD mechanical strain and NEW coefficients.
/// There is no elastic endpoint predictor, fraction, temperature or time claim.
#[derive(Debug)]
pub struct ThermalTensileTetActivation<'a> {
    request: ThermalTensileTetRequest<'a>,
    activations: Vec<Option<ResolvedTensileActivation<'a>>>,
}

/// Original source/state transfer and replay shared by both initiation stages.
/// Only the enclosing owner assigns endpoint-force or interval-mean semantics.
#[derive(Debug)]
struct ThermalTensileTetRequest<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedThermalBrittle],
    next_materials: &'a [ResolvedThermalBrittle],
    budget: TetAssemblyBudget,
    initial: SmallStrainTetSolution<()>,
    elastic_storage_j: f64,
    normals: Vec<[f64; 3]>,
    settings: TensileInitiationSettings,
    prescribed_end_m: Vec<(usize, f64)>,
    nodal_forces_end_n: Vec<(usize, f64)>,
    solver_settings: IncrementSettings,
}

/// A solved elastic endpoint is inspected without copying its nodal arrays.
#[derive(Debug)]
enum ThermalTensileEndpoint<'a> {
    Elastic(SmallStrainTetSolution<()>),
    Onset(ThermalTensileTetOnset<'a>),
}

/// Same prepared physical planes at the original and proposed temperatures.
/// These own source preparation only. Neither contains accepted fracture history.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalTensileTetBands {
    pub initial: Vec<ResolvedThermalBrittle>,
    pub next: Vec<ResolvedThermalBrittle>,
}

/// Original thermal endpoint under the native brittle representation.
#[derive(Debug)]
pub struct ThermalBrittleTetHandoff<'a> {
    pub solid: ThermalBrittleTetSolid<'a>,
    pub elastic_storage_j: f64,
    pub storage_defect_j: f64,
    pub maximum_stress_defect_pa: f64,
}

/// Native source, tensile detection, equilibrium, work or continuity refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalTensileTetInitiationError {
    Invalid(&'static str),
    Element {
        element: usize,
        source: ResolvedTensileInitiationError,
    },
    Thermal(ThermalBrittleTetError),
}

impl core::fmt::Display for ThermalTensileTetInitiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal tensile initiation refused: {self:?}")
    }
}
impl std::error::Error for ThermalTensileTetInitiationError {}
impl From<ThermalBrittleTetError> for ThermalTensileTetInitiationError {
    fn from(error: ThermalBrittleTetError) -> Self {
        Self::Thermal(error)
    }
}
impl From<BrittleTetError> for ThermalTensileTetInitiationError {
    fn from(error: BrittleTetError) -> Self {
        Self::Thermal(error.into())
    }
}
impl From<SmallStrainTetError> for ThermalTensileTetInitiationError {
    fn from(error: SmallStrainTetError) -> Self {
        BrittleTetError::from(error).into()
    }
}

/// An assembly view, not a second thermal constitutive or damage law.
#[derive(Clone, Copy)]
struct ThermalElasticView<'a>(&'a ResolvedThermalBrittle);

impl SmallStrainLaw for ThermalElasticView<'_> {
    type State = ();
    fn initial_state(&self) {}
    fn stress(&self, total: &Voigt, (): &()) -> Voigt {
        self.0
            .mechanical_strain(total)
            .map(|strain| self.0.band().law().elastic().stress(&strain, &()))
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, total: &Voigt, (): &()) -> Tangent6 {
        self.0
            .mechanical_strain(total)
            .map(|strain| self.0.band().law().elastic().tangent(&strain, &()))
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, _total: &Voigt, (): &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        self.0.band().law().elastic().admissibility()
    }
    fn card(&self) -> ModelCard {
        self.0.band().law().elastic().card()
    }
}

impl<'a> ThermalTensileElasticTetSolid<'a> {
    /// Bind virgin reference-temperature sources and the original native mesh.
    /// Dormant source planes are explicit inputs; their actual widths are
    /// derived from the receiving geometry only when fracture is selected.
    pub fn try_new(
        geometry: &'a TetLoadGeometry3<'a>,
        materials: &'a [ResolvedThermalBrittle],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermalTensileTetInitiationError> {
        checkpoint(cx)?;
        count(materials.len(), geometry.tetrahedra().len())?;
        if materials.is_empty()
            || geometry.nodes_m().is_empty()
            || geometry.nodes_m().len() > budget.maximum_nodes
            || materials.len() > budget.maximum_tetrahedra
            || geometry.nodes_m().len().checked_mul(3).is_none()
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
            || materials.iter().any(|m| m.free_linear_strain() != 0.0)
        {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "virgin zero-expansion sources and nonempty budgeted geometry required",
            ));
        }
        let solution = SmallStrainTetSolution {
            displacement_m: vec![[0.0; 3]; geometry.nodes_m().len()],
            element_strain: vec![[0.0; 6]; materials.len()],
            element_stress: vec![[0.0; 6]; materials.len()],
            element_states: vec![(); materials.len()],
            increments: Vec::new(),
        };
        checkpoint(cx)?;
        Ok(Self {
            geometry,
            materials,
            budget,
            solution,
            stored_energy_j: 0.0,
            supports: None,
        })
    }

    #[must_use]
    pub const fn geometry(&self) -> &'a TetLoadGeometry3<'a> {
        self.geometry
    }
    #[must_use]
    pub const fn materials(&self) -> &'a [ResolvedThermalBrittle] {
        self.materials
    }
    /// FE strains are TOTAL, including source-integrated free expansion.
    #[must_use]
    pub const fn solution(&self) -> &SmallStrainTetSolution<()> {
        &self.solution
    }
    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }

    /// Solve one source-bound endpoint and inspect its mechanical strain chord.
    /// The unique prescribed DOF set stays fixed; values are absolute targets.
    /// Coefficient activation at the held mechanical start returns separately
    /// before a displacement solve; neither prediction accepts physical state.
    #[allow(clippy::too_many_arguments)]
    pub fn probe<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        solver_settings: IncrementSettings,
        initiation_settings: TensileInitiationSettings,
        controls: BrittleTetControls,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetProbe<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        admit_controls(controls)?;
        let (supports, activation) = self.prepare_endpoint(
            next_materials,
            prescribed_m,
            nodal_forces_n,
            solver_settings,
            initiation_settings,
            direction_hints,
            cx,
        )?;
        if let Some(activation) = activation {
            return Ok(ThermalTensileTetProbe::Activation(activation));
        }
        let laws: Vec<_> = next_materials.iter().map(ThermalElasticView).collect();
        let problem = SmallStrainTetProblem {
            nodes_m: self.geometry.nodes_m(),
            tetrahedra: self.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings: solver_settings,
        };
        let solution = problem.solve_from(
            self.solution.displacement_m.as_flattened(),
            &self.solution.element_states,
            cx,
        )?;
        let solution = match self.inspect_endpoint(
            next_materials,
            solution,
            prescribed_m,
            nodal_forces_n,
            solver_settings,
            initiation_settings,
            direction_hints,
            cx,
        )? {
            ThermalTensileEndpoint::Elastic(solution) => solution,
            ThermalTensileEndpoint::Onset(onset) => {
                return Ok(ThermalTensileTetProbe::Onset(onset));
            }
        };
        let energy = self.recover_elastic(
            next_materials,
            &solution,
            prescribed_m,
            nodal_forces_n,
            controls,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(ThermalTensileTetProbe::Elastic(
            ThermalTensileElasticTetStep {
                solid: ThermalTensileElasticTetSolid {
                    geometry: self.geometry,
                    materials: next_materials,
                    budget: self.budget,
                    solution,
                    stored_energy_j: energy.stored_after_j,
                    supports: Some(supports),
                },
                energy,
            },
        ))
    }

    /// Admit the exact source/load request and inspect its coefficient stage.
    /// Native static and inertial consumers share this check before any solve.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn prepare_endpoint<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        solver_settings: IncrementSettings,
        initiation_settings: TensileInitiationSettings,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<
        (Vec<usize>, Option<ThermalTensileTetActivation<'next>>),
        ThermalTensileTetInitiationError,
    >
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        count(next_materials.len(), self.materials.len())?;
        // Activation can precede the FE solve; retain its native settings
        // admission even when no matrix assembly is needed for observation.
        if solver_settings.max_newton_iterations == 0
            || !solver_settings.relative_tolerance.is_finite()
            || solver_settings.relative_tolerance <= 0.0
            || !solver_settings.linear_tolerance.is_finite()
            || solver_settings.linear_tolerance <= 0.0
        {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "native Newton settings must admit the retained thermal request",
            ));
        }
        let dofs = self.solution.displacement_m.len() * 3;
        if prescribed_m.len() > dofs
            || nodal_forces_n.len() > dofs
            || direction_hints.is_some_and(|h| h.len() != self.materials.len())
            || prescribed_m
                .iter()
                .chain(nodal_forces_n)
                .any(|&(dof, value)| dof >= dofs || !value.is_finite())
        {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "finite endpoint arrays must fit the original DOF and element counts",
            ));
        }
        let mut supports: Vec<_> = prescribed_m.iter().map(|&(dof, _)| dof).collect();
        supports.sort_unstable();
        if supports.windows(2).any(|p| p[0] == p[1])
            || self.supports.as_ref().is_some_and(|old| old != &supports)
        {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "thermal elastic history requires the same unique prescribed DOF set",
            ));
        }
        // Ordinary elastic steps retain no per-element activation array. Only
        // a real activation allocates the original-order observation slots.
        let mut activations: Vec<Option<ResolvedTensileActivation<'next>>> = Vec::new();
        for (element, (old, next)) in self.materials.iter().zip(next_materials).enumerate() {
            checkpoint(cx)?;
            old.require_temperature_successor(next)
                .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
            let start = mechanical(old, &self.solution.element_strain[element], element)?;
            let activation = next
                .band()
                .tensile_initiator(initiation_settings)
                .and_then(|detector| {
                    detector.activation_at_controlled(
                        &start,
                        direction_hints.and_then(|h| h[element]),
                        &mut || cx.checkpoint().is_err(),
                    )
                })
                .map_err(|source| ThermalTensileTetInitiationError::Element { element, source })?;
            if let Some(activation) = activation {
                if activations.is_empty() {
                    activations.resize_with(next_materials.len(), || None);
                }
                activations[element] = Some(activation);
            }
        }
        let activation = if !activations.is_empty() {
            let normals = activations
                .iter()
                .enumerate()
                .map(|(element, event)| {
                    event.as_ref().map_or_else(
                        || self.materials[element].band().law().normal(),
                        |event| event.native().normal(),
                    )
                })
                .collect();
            Some(ThermalTensileTetActivation {
                request: self.initiation_request(
                    next_materials,
                    normals,
                    prescribed_m,
                    nodal_forces_n,
                    solver_settings,
                    initiation_settings,
                ),
                activations,
            })
        } else {
            None
        };
        checkpoint(cx)?;
        Ok((supports, activation))
    }

    /// Shared native strain-chord inspection after true endpoint stress recovery.
    #[allow(clippy::too_many_arguments)]
    fn inspect_endpoint<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        solution: SmallStrainTetSolution<()>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        solver_settings: IncrementSettings,
        initiation_settings: TensileInitiationSettings,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileEndpoint<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        let mut events: Vec<Option<ResolvedTensileOnset<'next>>> = Vec::new();
        let mut first_fraction = [1.0_f64; 2];
        let mut first_candidate = None;
        let mut smallest_upper = f64::INFINITY;
        for element in 0..next_materials.len() {
            checkpoint(cx)?;
            let next: &'next ResolvedThermalBrittle = &next_materials[element];
            let start = mechanical(
                &self.materials[element],
                &self.solution.element_strain[element],
                element,
            )?;
            let end = mechanical(next, &solution.element_strain[element], element)?;
            let event = next
                .band()
                .tensile_initiator(initiation_settings)
                .and_then(|detector| {
                    detector.between_controlled(
                        &start,
                        &end,
                        direction_hints.and_then(|h| h[element]),
                        &mut || cx.checkpoint().is_err(),
                    )
                })
                .map_err(|source| ThermalTensileTetInitiationError::Element { element, source })?;
            if let Some(event) = event {
                let bracket = event.native().fraction_bracket();
                first_fraction[0] = first_fraction[0].min(bracket[0]);
                first_fraction[1] = first_fraction[1].min(bracket[1]);
                if bracket[1] < smallest_upper {
                    smallest_upper = bracket[1];
                    first_candidate = Some(element);
                }
                if events.is_empty() {
                    events.resize_with(next_materials.len(), || None);
                }
                events[element] = Some(event);
            }
        }
        checkpoint(cx)?;
        if let Some(first_candidate) = first_candidate {
            let normals = events
                .iter()
                .enumerate()
                .map(|(element, event)| {
                    event.as_ref().map_or_else(
                        || self.materials[element].band().law().normal(),
                        |event| event.native().normal(),
                    )
                })
                .collect();
            return Ok(ThermalTensileEndpoint::Onset(ThermalTensileTetOnset {
                request: self.initiation_request(
                    next_materials,
                    normals,
                    prescribed_m,
                    nodal_forces_n,
                    solver_settings,
                    initiation_settings,
                ),
                predictor: solution,
                events,
                first_fraction,
                first_candidate,
            }));
        }
        Ok(ThermalTensileEndpoint::Elastic(solution))
    }

    #[allow(clippy::too_many_arguments)]
    fn initiation_request<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        normals: Vec<[f64; 3]>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        solver_settings: IncrementSettings,
        settings: TensileInitiationSettings,
    ) -> ThermalTensileTetRequest<'next>
    where
        'a: 'next,
    {
        ThermalTensileTetRequest {
            geometry: self.geometry,
            materials: self.materials,
            next_materials,
            budget: self.budget,
            initial: self.solution.clone(),
            elastic_storage_j: self.stored_energy_j,
            normals,
            settings,
            prescribed_end_m: prescribed_m.to_vec(),
            nodal_forces_end_n: nodal_forces_n.to_vec(),
            solver_settings,
        }
    }

    /// Shared exact elastic coefficient and strain-work observations. Applied
    /// forces, support reactions and inertia remain with the receiving owner.
    fn elastic_material_energy(
        &self,
        next_materials: &[ResolvedThermalBrittle],
        solution: &SmallStrainTetSolution<()>,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetEnergy, ThermalTensileTetInitiationError> {
        let mut energy = ThermalBrittleTetEnergy {
            stored_before_j: self.stored_energy_j,
            ..ThermalBrittleTetEnergy::default()
        };
        for (element, next) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            let total0 = &self.solution.element_strain[element];
            let total1 = &solution.element_strain[element];
            let strain0 = mechanical(&self.materials[element], total0, element)?;
            let strain1 = mechanical(next, total1, element)?;
            let stress0 = next.band().law().elastic().stress(&strain0, &());
            let stress1 = &solution.element_stress[element];
            let mean = core::array::from_fn(|i| 0.5 * stress0[i] + 0.5 * stress1[i]);
            let old_density = 0.5 * contract(&self.solution.element_stress[element], &strain0);
            let prepared_density = 0.5 * contract(&stress0, &strain0);
            let density = 0.5 * contract(stress1, &strain1);
            let total_change = core::array::from_fn(|i| total1[i] - total0[i]);
            let mechanical_change = core::array::from_fn(|i| strain1[i] - strain0[i]);
            let free_change =
                next.free_linear_strain() - self.materials[element].free_linear_strain();
            let mechanical_work = contract(&mean, &total_change);
            let eigenstrain_work = -free_change * (mean[0] + mean[1] + mean[2]);
            let constitutive_work = contract(&mean, &mechanical_change);
            let internal_work = contract(stress1, &total_change);
            if [old_density, prepared_density, density]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
                || [
                    mechanical_work,
                    eigenstrain_work,
                    constitutive_work,
                    internal_work,
                ]
                .iter()
                .any(|v| !v.is_finite())
            {
                return Err(ThermalTensileTetInitiationError::Invalid(
                    "finite elastic storage and thermal work required",
                ));
            }
            let volume = self.geometry.element_volumes_m3()[element];
            energy.stored_after_j = volume.mul_add(density, energy.stored_after_j);
            energy.parameter_work_j =
                volume.mul_add(prepared_density - old_density, energy.parameter_work_j);
            energy.mechanical_work_j = volume.mul_add(mechanical_work, energy.mechanical_work_j);
            energy.eigenstrain_work_j = volume.mul_add(eigenstrain_work, energy.eigenstrain_work_j);
            energy.constitutive_strain_work_j =
                volume.mul_add(constitutive_work, energy.constitutive_strain_work_j);
            energy.internal_endpoint_work_j =
                volume.mul_add(internal_work, energy.internal_endpoint_work_j);
        }
        energy.strain_split_error_j = energy.constitutive_strain_work_j
            - (energy.mechanical_work_j + energy.eigenstrain_work_j);
        energy.material_balance_error_j = energy.parameter_work_j
            + energy.constitutive_strain_work_j
            - (energy.stored_after_j - energy.stored_before_j);
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - energy.mechanical_work_j;
        Ok(energy)
    }

    /// Quadratic native elasticity gives the exact mean stress under the next
    /// coefficients. This is the same work split as the original thermal owner,
    /// without creating a dormant fracture history to obtain elastic response.
    #[allow(clippy::too_many_arguments)]
    fn recover_elastic(
        &self,
        next_materials: &[ResolvedThermalBrittle],
        solution: &SmallStrainTetSolution<()>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetEnergy, ThermalTensileTetInitiationError> {
        let mut energy = self.elastic_material_energy(next_materials, solution, cx)?;
        let old_u = self.solution.displacement_m.as_flattened();
        let new_u = solution.displacement_m.as_flattened();
        for &(dof, force) in nodal_forces_n {
            checkpoint(cx)?;
            energy.dead_load_endpoint_work_j =
                force.mul_add(new_u[dof] - old_u[dof], energy.dead_load_endpoint_work_j);
        }
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&solution.increments[0].reactions_n) {
            checkpoint(cx)?;
            energy.support_endpoint_work_j =
                reaction.mul_add(new_u[dof] - old_u[dof], energy.support_endpoint_work_j);
        }
        energy.strain_split_error_j = energy.constitutive_strain_work_j
            - (energy.mechanical_work_j + energy.eigenstrain_work_j);
        energy.material_balance_error_j = energy.parameter_work_j
            + energy.constitutive_strain_work_j
            - (energy.stored_after_j - energy.stored_before_j);
        let external = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j;
        energy.equilibrium_work_residual_j = external - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - energy.mechanical_work_j;
        energy.total_balance_error_j =
            external + energy.parameter_work_j + energy.eigenstrain_work_j
                - (energy.stored_after_j - energy.stored_before_j)
                - energy.endpoint_work_defect_j;
        for (account, defect, tolerance) in [
            (
                "thermal elastic strain work split",
                energy.strain_split_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic material balance",
                energy.material_balance_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic complete work balance",
                energy.total_balance_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic equilibrium work",
                energy.equilibrium_work_residual_j,
                controls.equilibrium_work_tolerance_j,
            ),
            (
                "thermal elastic endpoint defect",
                energy.endpoint_work_defect_j,
                controls.maximum_endpoint_work_defect_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        Ok(energy)
    }
}

impl<'a> ThermalTensileTetOnset<'a> {
    #[must_use]
    pub fn events(&self) -> &[Option<ResolvedTensileOnset<'a>>] {
        &self.events
    }
    #[must_use]
    pub const fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        &self.request.initial
    }
    /// Disposable elastic endpoint; its over-strength stress is never accepted.
    #[must_use]
    pub const fn elastic_predictor(&self) -> &SmallStrainTetSolution<()> {
        &self.predictor
    }
    /// Numerical first crossing bracket of the MECHANICAL strain chord only.
    #[must_use]
    pub const fn first_fraction_bracket(&self) -> [f64; 2] {
        self.first_fraction
    }
    #[must_use]
    pub const fn first_candidate_element(&self) -> usize {
        self.first_candidate
    }

    /// Prepare the selected planes at both exact original source states.
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetBands, ThermalTensileTetInitiationError> {
        self.request.resolve_bands(cx)
    }

    /// Preserve the OLD accepted state without charging preload work again.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedThermalBrittle],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetHandoff<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        self.request.handoff(materials, energy_tolerance_j, cx)
    }

    /// Replay the retained source/load endpoint through native thermal fracture.
    pub fn replay<'next>(
        &self,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetStep<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        self.request.replay(bands, controls, cx)
    }
}

impl<'a> ThermalTensileTetActivation<'a> {
    /// Held mechanical-strain observations in original element order. An empty
    /// entry keeps that element's explicitly declared dormant crack plane.
    #[must_use]
    pub fn activations(&self) -> &[Option<ResolvedTensileActivation<'a>>] {
        &self.activations
    }

    /// Original accepted TOTAL strain/stress and displacement, before weakening.
    #[must_use]
    pub const fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        &self.request.initial
    }

    /// Prepare the selected planes at both exact original source states.
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetBands, ThermalTensileTetInitiationError> {
        self.request.resolve_bands(cx)
    }

    /// Preserve the OLD accepted state without charging preload work again.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedThermalBrittle],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetHandoff<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        self.request.handoff(materials, energy_tolerance_j, cx)
    }

    /// Replay the retained source/load endpoint through native thermal fracture.
    pub fn replay<'next>(
        &self,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetStep<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        self.request.replay(bands, controls, cx)
    }
}

impl<'a> ThermalTensileTetRequest<'a> {
    /// Preserve both exact thermal sources while preparing discovered planes.
    /// Untriggered elements retain their original declared planes. No source
    /// lookup, temperature interpolation or history mutation occurs here.
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetBands, ThermalTensileTetInitiationError> {
        checkpoint(cx)?;
        let mut initial = Vec::with_capacity(self.materials.len());
        let mut next = Vec::with_capacity(self.materials.len());
        for element in 0..self.materials.len() {
            checkpoint(cx)?;
            let normal = self.normal(element);
            let geometry = TetCrackBandGeometry::from_element(self.geometry, element, normal)
                .map_err(BrittleTetError::Geometry)?;
            let prepare = |source: &ResolvedThermalBrittle| {
                let mut settings = source.band().law().settings();
                settings.normal = normal;
                settings.band_width_m = geometry.band_width_m();
                let band = ResolvedBrittleCrackBand::try_from_material(
                    source.band().material().clone(),
                    settings,
                )?;
                source.with_resolved_band(band)
            };
            initial.push(
                prepare(&self.materials[element])
                    .map_err(|source| ThermalBrittleTetError::Source { element, source })?,
            );
            next.push(
                prepare(&self.next_materials[element])
                    .map_err(|source| ThermalBrittleTetError::Source { element, source })?,
            );
        }
        checkpoint(cx)?;
        Ok(ThermalTensileTetBands { initial, next })
    }

    /// Transfer the OLD accepted total/mechanical strain, stress and storage.
    /// The new-temperature predictor contributes no accepted work or history.
    /// Native initial trials must remain undamaged and stationary, including
    /// mixed-mode history; a previously inelastic response cannot be erased.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedThermalBrittle],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetHandoff<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        count(materials.len(), self.materials.len())?;
        if !energy_tolerance_j.is_finite() || energy_tolerance_j < 0.0 {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "finite nonnegative handoff energy tolerance required",
            ));
        }
        let mut bands = Vec::with_capacity(materials.len());
        let mut states = Vec::with_capacity(materials.len());
        let mut stresses = Vec::with_capacity(materials.len());
        let mut storage = 0.0;
        let mut maximum_stress_defect_pa = 0.0_f64;
        for (element, receiving) in materials.iter().enumerate() {
            checkpoint(cx)?;
            let geometry = self.admit_band(&self.materials[element], receiving, element)?;
            let strain = mechanical(
                &self.materials[element],
                &self.initial.element_strain[element],
                element,
            )?;
            let law = receiving.band().law();
            let state = law
                .trial(&strain, &law.initial_state())
                .map_err(|source| BrittleTetError::Material { element, source })?
                .state;
            if state.damage() != 0.0
                || state.dissipated_j_m3() != 0.0
                || state.cohesive_point().mixed_mode_history().is_some()
            {
                return Err(ThermalTensileTetInitiationError::Invalid(
                    "thermal handoff cannot erase prior mixed-mode initiation or damage",
                ));
            }
            let held = law
                .trial(&strain, &state)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            if held.state != state || held.energy.work != 0.0 || held.energy.dissipated != 0.0 {
                return Err(ThermalTensileTetInitiationError::Invalid(
                    "thermal handoff history must be stationary",
                ));
            }
            for (&new, &old) in held
                .stress_pa
                .iter()
                .zip(&self.initial.element_stress[element])
            {
                let defect = (new - old).abs();
                if !defect.is_finite()
                    || defect
                        > self.settings.relative_stress_tolerance
                            * self.materials[element]
                                .band()
                                .law()
                                .parameters()
                                .tensile_strength_pa
                {
                    return Err(ThermalTensileTetInitiationError::Invalid(
                        "thermal handoff changes original elastic stress",
                    ));
                }
                maximum_stress_defect_pa = maximum_stress_defect_pa.max(defect);
            }
            let density = state.stored_energy_j_m3();
            let old_density = 0.5 * contract(&self.initial.element_stress[element], &strain);
            if !density.is_finite()
                || density < 0.0
                || !old_density.is_finite()
                || (density - old_density).abs() > law.settings().energy_tolerance_j_m3
            {
                return Err(ThermalTensileTetInitiationError::Invalid(
                    "thermal handoff changes original elastic storage",
                ));
            }
            storage = geometry.reference_volume_m3().mul_add(density, storage);
            bands.push(geometry);
            states.push(state);
            stresses.push(held.stress_pa);
        }
        let storage_defect_j = storage - self.elastic_storage_j;
        energy_gate(
            "thermal elastic-to-brittle storage transfer",
            storage_defect_j,
            energy_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(ThermalBrittleTetHandoff {
            solid: ThermalBrittleTetSolid {
                geometry: self.geometry,
                materials,
                bands: bands.into(),
                released_faces: Arc::from([]),
                budget: self.budget,
                solution: SmallStrainTetSolution {
                    displacement_m: self.initial.displacement_m.clone(),
                    element_strain: self.initial.element_strain.clone(),
                    element_stress: stresses,
                    element_states: states,
                    increments: Vec::new(),
                },
                stored_energy_j: storage,
                fracture_dissipation_j: 0.0,
            },
            elastic_storage_j: self.elastic_storage_j,
            storage_defect_j,
            maximum_stress_defect_pa,
        })
    }

    /// Replay exactly the retained thermal/load request through the native
    /// brittle owner. The supplied preparations must match BOTH original sources.
    pub fn replay<'next>(
        &self,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetStep<'next>, ThermalTensileTetInitiationError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        admit_controls(controls)?;
        self.admit_next_materials(&bands.next, cx)?;
        let handoff = self.handoff(&bands.initial, controls.material_energy_tolerance_j, cx)?;
        Ok(handoff.solid.advance(
            &bands.next,
            &self.prescribed_end_m,
            &self.nodal_forces_end_n,
            self.solver_settings,
            controls,
            cx,
        )?)
    }

    fn admit_next_materials(
        &self,
        materials: &[ResolvedThermalBrittle],
        cx: &Cx<'_>,
    ) -> Result<(), ThermalTensileTetInitiationError> {
        count(materials.len(), self.next_materials.len())?;
        for (element, next) in materials.iter().enumerate() {
            checkpoint(cx)?;
            self.admit_band(&self.next_materials[element], next, element)?;
        }
        Ok(())
    }

    fn normal(&self, element: usize) -> [f64; 3] {
        self.normals[element]
    }

    fn admit_band(
        &self,
        original: &ResolvedThermalBrittle,
        receiving: &ResolvedThermalBrittle,
        element: usize,
    ) -> Result<TetCrackBandGeometry, ThermalTensileTetInitiationError> {
        let rebound = original
            .with_resolved_band(receiving.band().clone())
            .map_err(|source| ThermalBrittleTetError::Source { element, source })?;
        if rebound != *receiving {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "receiving band changed the original thermal expansion source",
            ));
        }
        let expected = self.normal(element);
        let actual = receiving.band().law().normal();
        let sign = if expected.iter().zip(actual).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
            -1.0
        } else {
            1.0
        };
        if (0..3).any(|i| (expected[i] - sign * actual[i]).abs() > 128.0 * f64::EPSILON) {
            return Err(ThermalTensileTetInitiationError::Invalid(
                "receiving thermal band changed the discovered or declared plane",
            ));
        }
        let geometry = TetCrackBandGeometry::from_element(self.geometry, element, actual)
            .map_err(BrittleTetError::Geometry)?;
        let required = geometry.band_width_m();
        let supplied = receiving.band().law().band_width_m();
        if (required - supplied).abs() > 64.0 * f64::EPSILON * required.max(supplied) {
            return Err(BrittleTetError::WidthMismatch {
                element,
                required_m: required,
                material_m: supplied,
            }
            .into());
        }
        Ok(geometry)
    }
}

fn count(actual: usize, expected: usize) -> Result<(), ThermalTensileTetInitiationError> {
    if actual != expected {
        return Err(BrittleTetError::MaterialCount {
            expected,
            got: actual,
        }
        .into());
    }
    Ok(())
}

fn mechanical(
    source: &ResolvedThermalBrittle,
    total: &Voigt,
    element: usize,
) -> Result<Voigt, ThermalTensileTetInitiationError> {
    Ok(source
        .mechanical_strain(total)
        .map_err(|source| ThermalBrittleTetError::Source { element, source })?)
}

fn admit_controls(controls: BrittleTetControls) -> Result<(), ThermalTensileTetInitiationError> {
    if [
        controls.material_energy_tolerance_j,
        controls.equilibrium_work_tolerance_j,
        controls.maximum_endpoint_work_defect_j,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(ThermalTensileTetInitiationError::Invalid(
            "finite nonnegative thermal work tolerances required",
        ));
    }
    Ok(())
}
