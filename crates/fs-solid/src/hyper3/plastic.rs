//! Retained finite-deformation inelasticity on the native P1 tetrahedra.
//!
//! Each physical endpoint uses the original hyper3 geometry, assembly, tangent,
//! Newton and globalization. Every trial sees the SAME accepted plastic metric
//! and accumulated flow; a successful endpoint returns a new complete solid.
//! Rejection or cancellation leaves the original solid untouched. Both native
//! and source-resolved laws remain borrowed, without extracting source data.
//!
//! The native models use multiplicative Hencky elasticity with Kirchhoff J2,
//! Norton creep, optional primary resistance/recovery and explicit rate flow.
//! Their reduced incremental potential drives Newton. Physical storage,
//! inelastic dissipation, integrated flow potential, endpoint work
//! quadrature error and equilibrium work residual are separate quantities.
//! Endpoint work is a first-order approximation along a physical increment;
//! its signed defect need not be positive or invariant under a finite rotation
//! of a prestressed body. It is never added to material heat or dissipation.
//!
//! This is a compressible displacement-only P1 formulation. Isochoric plastic
//! flow can cause volumetric locking even when the elastic Poisson ratio is
//! moderate. Accurate bending in strongly plastic or nearly incompressible
//! regimes requires a mixed/stabilized formulation and mesh evidence; an affine
//! patch does not establish that accuracy. [`super::mixed::plastic`] supplies
//! the native P2/P1 mixed successor for J2, Perzyna, Norton and combined or
//! recovering laws, with original bulk compliance or exact incompressibility.
//! Thermal evolution here retains its displacement-only formulation.
//! [`dynamics`] adds endpoint inertia
//! on this same solve. [`thermal`] adds changing temperature and exact native
//! thermoelastic/plastic accounts. Contact, follower loads and cyclic kinematic
//! hardening are outside these paths.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::plastic::finite::{
    FiniteInelasticMaterial, FiniteJ2, FiniteJ2Error, FiniteJ2State,
};

use super::{
    Element, HyperTetError, HyperTetProblem, HyperTetSettings, HyperTetSolution,
    HyperelasticLawField, checkpoint,
};
use crate::linear3::TetAssemblyBudget;

/// Progressive ductile softening through the shared nonpotential Newton owner.
pub mod ductile;
/// Consistent reference inertia with atomic finite-plastic motion acceptance.
pub mod dynamics;
/// Fully integrated finite Hex8 plasticity with optional mean dilatation.
pub mod hex;
/// Shared frictionless plane tooling on retained native tetrahedral responses.
pub mod plane_contact;
/// Retained pressure-sensitive finite compaction on the same tetrahedral owner.
pub mod porous;
/// Retained finite plastic history with evolving thermal coefficients and stretch.
pub mod thermal;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

type AcceptedMotion = (Vec<[f64; 3]>, crate::inertia::InertialEnergy);

/// Exact original laws, in original element order; no property averaging.
#[derive(Debug)]
pub enum FinitePlasticLawField<'a, L = FiniteJ2> {
    /// One native or source-resolved law shared by the mesh.
    Uniform(&'a L),
    /// One native or source-resolved law per element.
    PerElement(&'a [L]),
}

impl<L> Copy for FinitePlasticLawField<'_, L> {}

impl<L> Clone for FinitePlasticLawField<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<L> FinitePlasticLawField<'_, L> {
    fn validate_count(&self, expected: usize) -> Result<(), HyperTetError> {
        if let Self::PerElement(laws) = self
            && laws.len() != expected
        {
            return Err(HyperTetError::LawCountMismatch {
                expected,
                got: laws.len(),
            });
        }
        Ok(())
    }

    fn at(&self, element: usize) -> &L {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
        }
    }
}

impl<L: FiniteInelasticMaterial> FinitePlasticLawField<'_, L> {
    fn common_duration_s(
        &self,
        elements: usize,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, FinitePlasticTetError> {
        let mut common = None;
        for element in 0..elements {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            if let Some(duration) = self.at(element).duration_s() {
                if !duration.is_finite()
                    || duration <= 0.0
                    || common.is_some_and(|old| old != duration)
                {
                    return Err(invalid(
                        "every rate-dependent element must share one finite positive physical duration",
                    ));
                }
                common = Some(duration);
            }
        }
        Ok(common)
    }
}

/// Explicit absolute physical-step gates, in addition to native force gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePlasticStepControls {
    /// Maximum absolute support/dead-load versus internal endpoint work residual.
    /// Also bounds independent constitutive-versus-assembled potential checks.
    pub equilibrium_work_tolerance_j: f64,
    /// Maximum absolute signed endpoint integration defect, J.
    /// Reduce the physical increment when this gate refuses it.
    pub maximum_endpoint_work_defect_j: f64,
    /// Maximum allowed positive yield residual of any accepted material state.
    pub yield_tolerance_pa: f64,
}

/// Physical accounting for one endpoint, all energies per whole reference mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct FinitePlasticTetEnergy {
    /// Accepted elastic plus hardening storage before this increment.
    pub stored_before_j: f64,
    /// Independently evaluated elastic plus hardening storage after acceptance.
    pub stored_after_j: f64,
    /// Native plastic/viscous channel. For the native sole-Norton convention
    /// this is creep heat; combined laws report creep in its separate channel.
    pub plastic_dissipation_j: f64,
    /// Physical creep heat from the native combined creep/plastic channel.
    pub creep_dissipation_j: f64,
    /// Integrated native yield/kinetic potential, distinct from physical heat.
    pub flow_potential_j: f64,
    /// Reduced incremental potential used by Newton, not physical storage.
    pub incremental_potential_j: f64,
    /// Endpoint dead nodal/body force dotted with the actual displacement change.
    pub dead_load_endpoint_work_j: f64,
    /// Endpoint support reaction dotted with the actual prescribed motion.
    pub support_endpoint_work_j: f64,
    /// Sum(V P_after : (F_after - F_before)).
    pub internal_endpoint_work_j: f64,
    /// External endpoint work minus material internal work and, for dynamics,
    /// inertial endpoint work; signed and independently gate checked.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus storage change minus total inelastic heat.
    /// Signed quadrature defect; never heat or an assumed nonnegative quantity.
    pub endpoint_work_defect_j: f64,
}

impl FinitePlasticTetEnergy {
    /// Complete physical heat from the native inelastic mechanisms. Numerical
    /// endpoint defects and the reduced flow potential are never caloric heat.
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

/// Original mesh, original law binding and complete immutable accepted history.
///
/// Construction admits only the virgin reference state. Advance returns a new
/// solid after all gates pass; neither histories nor displacement can be edited
/// independently. Reference geometry is prepared once and shared across steps.
#[derive(Debug)]
pub struct FinitePlasticTetSolid<'a, L: FiniteInelasticMaterial = FiniteJ2> {
    material: FinitePlasticLawField<'a, L>,
    duration_s: Option<f64>,
    state: FinitePlasticTetState<'a, L::State>,
}

/// Numerical history shared by fixed-temperature, porous and thermal owners.
/// No constitutive source is stored or reinterpreted by this common kernel.
#[derive(Debug)]
struct FinitePlasticTetState<'a, S = FiniteJ2State> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    elements: Arc<Vec<Element>>,
    budget: TetAssemblyBudget,
    displacement_m: Vec<[f64; 3]>,
    element_deformation: Vec<[f64; 9]>,
    element_states: Vec<S>,
    element_first_piola_pa: Vec<[f64; 9]>,
    stored_energy_j: f64,
}

struct EndpointWork {
    internal_j: f64,
    dead_load_j: f64,
    support_j: f64,
    residual_j: f64,
}

fn reference_size(
    nodes: usize,
    elements: usize,
    budget: TetAssemblyBudget,
) -> Result<usize, FinitePlasticTetError> {
    if nodes == 0
        || elements == 0
        || nodes > budget.maximum_nodes
        || elements > budget.maximum_tetrahedra
    {
        return Err(invalid(
            "finite plastic mesh is empty or exceeds its reference budget",
        ));
    }
    nodes
        .checked_mul(3)
        .ok_or_else(|| invalid("finite plastic displacement size overflow"))
}

fn validate_controls(
    settings: HyperTetSettings,
    controls: FinitePlasticStepControls,
) -> Result<(), FinitePlasticTetError> {
    if settings.load_steps != 1
        || [
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
            controls.yield_tolerance_pa,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(invalid(
            "one physical endpoint and finite nonnegative work/yield tolerances required",
        ));
    }
    Ok(())
}

impl<'a> FinitePlasticTetState<'a> {
    fn reference_elements<L: FiniteStrainPotential>(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        material: HyperelasticLawField<'_, L>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Arc<Vec<Element>>, FinitePlasticTetError> {
        let dofs = reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        let prescribed: Vec<_> = (0..dofs).map(|dof| (dof, 0.0)).collect();
        let problem = HyperTetProblem {
            nodes_m,
            tetrahedra,
            material,
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget,
            settings: HyperTetSettings::default(),
        };
        let elements = problem.validate()?;
        checkpoint(cx)?;
        Ok(Arc::new(elements))
    }
}

impl<S> FinitePlasticTetState<'_, S> {
    #[allow(clippy::too_many_arguments)]
    fn solve_endpoint<L: FiniteStrainPotential>(
        &self,
        material: HyperelasticLawField<'_, L>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(HyperTetSolution, Vec<f64>), FinitePlasticTetError> {
        let problem = HyperTetProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            material,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            budget: self.budget,
            settings,
        };
        let equilibrium = problem.solve_loaded_from(
            Some(self.displacement_m.as_flattened()),
            Some(&self.elements),
            None,
            None,
            inertia,
            cx,
        )?;
        let force = problem.external(&self.elements, 1.0);
        Ok((equilibrium, force))
    }

    fn endpoint_work(
        &self,
        equilibrium: &HyperTetSolution,
        force: &[f64],
        deformation: &[[f64; 9]],
        stress: &[[f64; 9]],
        inertial_work_j: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<EndpointWork, FinitePlasticTetError> {
        self.endpoint_work_fields(
            equilibrium.displacement_m.as_flattened(),
            equilibrium.reactions_n.iter().copied(),
            force,
            deformation,
            stress,
            inertial_work_j,
            cx,
        )
    }

    /// Shared physical work for potential and nonpotential endpoint solvers.
    /// No scalar material potential is inferred from an accepted force field.
    #[allow(clippy::too_many_arguments)]
    fn endpoint_work_fields(
        &self,
        new_u: &[f64],
        reactions_n: impl IntoIterator<Item = (usize, f64)>,
        force: &[f64],
        deformation: &[[f64; 9]],
        stress: &[[f64; 9]],
        inertial_work_j: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<EndpointWork, FinitePlasticTetError> {
        let mut internal_j = 0.0;
        for (index, element) in self.elements.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let density = (0..9).fold(0.0, |sum, component| {
                stress[index][component].mul_add(
                    deformation[index][component] - self.element_deformation[index][component],
                    sum,
                )
            });
            internal_j = element.volume.mul_add(density, internal_j);
        }
        let old_u = self.displacement_m.as_flattened();
        let dead_load_j = force
            .iter()
            .zip(new_u.iter().zip(old_u))
            .fold(0.0, |sum, (force, (after, before))| {
                force.mul_add(after - before, sum)
            });
        let support_j = reactions_n.into_iter().fold(0.0, |sum, (dof, force)| {
            force.mul_add(new_u[dof] - old_u[dof], sum)
        });
        let mut residual_j = dead_load_j + support_j - internal_j;
        if let Some(inertial_work_j) = inertial_work_j {
            residual_j -= inertial_work_j;
        }
        Ok(EndpointWork {
            internal_j,
            dead_load_j,
            support_j,
            residual_j,
        })
    }
}

/// A fully accepted next solid and its native equilibrium/physical-step evidence.
#[derive(Debug)]
pub struct FinitePlasticTetStep<'a, L: FiniteInelasticMaterial = FiniteJ2> {
    /// Accepted full next configuration, source binding and material history.
    pub solid: FinitePlasticTetSolid<'a, L>,
    /// Native displacement, support reactions and Newton convergence history.
    /// `strain_energy_j` here is the reduced incremental potential.
    pub equilibrium: HyperTetSolution,
    /// Independently checked physical storage, dissipation and endpoint work.
    pub energy: FinitePlasticTetEnergy,
}

/// Refusal at the original numerical/material owner or explicit physical gate.
#[derive(Debug, Clone, PartialEq)]
pub enum FinitePlasticTetError {
    /// Original geometry, load, force-equilibrium, budget or cancellation refusal.
    Native(HyperTetError),
    /// Original finite plastic law refused an accepted or proposed history.
    Material {
        /// Index in the original tetrahedron order.
        element: usize,
        /// Unmodified native material diagnosis.
        source: FiniteJ2Error,
    },
    /// A finite physical-work tolerance was exceeded, or its value overflowed.
    EnergyGate {
        /// Independently checked quantity.
        quantity: &'static str,
        /// Signed value of the refused quantity.
        value_j: f64,
        /// Caller-specified absolute allowance.
        maximum_abs_j: f64,
    },
}

impl core::fmt::Display for FinitePlasticTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "finite plastic tetrahedron {element}: {source}")
            }
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "finite plastic {quantity}: {value_j} J exceeds absolute budget {maximum_abs_j} J"
            ),
        }
    }
}

impl std::error::Error for FinitePlasticTetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Native(error) => Some(error),
            Self::Material { source, .. } => Some(source),
            Self::EnergyGate { .. } => None,
        }
    }
}

impl From<HyperTetError> for FinitePlasticTetError {
    fn from(error: HyperTetError) -> Self {
        Self::Native(error)
    }
}

struct FrozenPoint<'a, L: FiniteInelasticMaterial = FiniteJ2> {
    law: &'a L,
    accepted: &'a L::State,
}

fn material_error(error: FiniteJ2Error) -> MaterialError {
    MaterialError::State {
        what: error.to_string(),
    }
}

impl<L: FiniteInelasticMaterial> FiniteStrainPotential for FrozenPoint<'_, L> {
    fn energy_density(&self, deformation: &[f64; 9]) -> f64 {
        self.law
            .stress_and_potential(deformation, self.accepted)
            .map_or(f64::NAN, |(_, potential)| potential)
    }

    /// The increment of the same frozen-history potential. Endpoint
    /// subtraction keeps its bits wherever it resolves the change. Below 2^20
    /// ulps of the larger endpoint potential, the subtraction is rounding of
    /// two large potentials, while a static line search compares the increment
    /// with a slope that falls quadratically near equilibrium. There the
    /// increment is the work of the conjugate stress along the straight path
    /// (Simpson's rule, exact for a cubic potential).
    fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
        let (Ok((stress_before, potential_before)), Ok((stress_after, potential_after))) = (
            self.law.stress_and_potential(before, self.accepted),
            self.law.stress_and_potential(after, self.accepted),
        ) else {
            return f64::NAN;
        };
        let difference = potential_after - potential_before;
        let resolution =
            1_048_576.0 * f64::EPSILON * potential_after.abs().max(potential_before.abs());
        if !(difference.abs() < resolution) {
            return difference;
        }
        let middle: [f64; 9] = core::array::from_fn(|i| 0.5 * (before[i] + after[i]));
        let Ok((stress_middle, _)) = self.law.stress_and_potential(&middle, self.accepted) else {
            return difference;
        };
        (0..9).fold(0.0, |work, i| {
            (stress_before[i] + 4.0 * stress_middle[i] + stress_after[i])
                .mul_add(after[i] - before[i], work)
        }) / 6.0
    }

    fn first_piola(&self, deformation: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.law
            .stress_and_potential(deformation, self.accepted)
            .map(|(stress, _)| stress)
            .map_err(material_error)
    }

    fn consistent_tangent(&self, deformation: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.law
            .response(deformation, self.accepted)
            .map(|response| response.tangent_pa)
            .map_err(material_error)
    }

    fn stress_and_potential(
        &self,
        deformation: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        self.law
            .stress_and_potential(deformation, self.accepted)
            .map_err(material_error)
    }
}

fn invalid(what: impl Into<String>) -> FinitePlasticTetError {
    HyperTetError::InvalidInput { what: what.into() }.into()
}

fn energy_gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FinitePlasticTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FinitePlasticTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

impl<'a, L: FiniteInelasticMaterial> FinitePlasticTetSolid<'a, L> {
    /// Prepare an undeformed, stress-free mesh with virgin native histories.
    ///
    /// # Errors
    /// Original mesh/material-count/size/quality refusals or cancellation.
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        material: FinitePlasticLawField<'a, L>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePlasticTetError> {
        checkpoint(cx)?;
        material.validate_count(tetrahedra.len())?;
        reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        let duration_s = material.common_duration_s(tetrahedra.len(), cx)?;
        let element_states: Vec<_> = (0..tetrahedra.len())
            .map(|index| material.at(index).initial_state())
            .collect();
        let frozen: Vec<_> = element_states
            .iter()
            .enumerate()
            .map(|(index, accepted)| FrozenPoint {
                law: material.at(index),
                accepted,
            })
            .collect();
        let elements = FinitePlasticTetState::reference_elements(
            nodes_m,
            tetrahedra,
            HyperelasticLawField::PerElement(&frozen),
            budget,
            cx,
        )?;
        Ok(Self {
            material,
            duration_s,
            state: FinitePlasticTetState {
                nodes_m,
                tetrahedra,
                elements,
                budget,
                displacement_m: vec![[0.0; 3]; nodes_m.len()],
                element_deformation: vec![IDENTITY; tetrahedra.len()],
                element_states,
                element_first_piola_pa: vec![[0.0; 9]; tetrahedra.len()],
                stored_energy_j: 0.0,
            },
        })
    }

    /// Original immutable reference coordinates [m], in native node order.
    pub fn reference_nodes_m(&self) -> &'a [[f64; 3]] {
        self.state.nodes_m
    }
    /// Original ordered material connectivity; deformation never rewrites it.
    pub fn reference_tetrahedra(&self) -> &'a [[usize; 4]] {
        self.state.tetrahedra
    }
    /// Original size and geometric-quality admission bounds.
    pub fn assembly_budget(&self) -> TetAssemblyBudget {
        self.state.budget
    }

    /// Accepted total displacement, in original reference-node order.
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.state.displacement_m
    }
    /// Accepted total deformation gradients, row-major, in tetrahedron order.
    pub fn element_deformation(&self) -> &[[f64; 9]] {
        &self.state.element_deformation
    }
    /// Complete native histories, including independent resistance when selected.
    pub fn element_states(&self) -> &[L::State] {
        &self.state.element_states
    }
    /// Shared physical interval for explicit kinetics; `None` for a completely
    /// rate-independent field. Each accepted call advances this same interval.
    pub const fn duration_s(&self) -> Option<f64> {
        self.duration_s
    }
    /// Accepted first Piola stress, row-major, in tetrahedron order.
    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.state.element_first_piola_pa
    }
    /// Physical elastic plus hardening storage over the whole reference mesh.
    pub fn stored_energy_j(&self) -> f64 {
        self.state.stored_energy_j
    }
    /// Original borrowed law owners, including complete model cards when supplied.
    pub fn material(&self) -> FinitePlasticLawField<'a, L> {
        self.material
    }

    /// Solve and accept one physical endpoint with absolute supports/dead loads.
    ///
    /// Every Newton and line-search evaluation borrows the old accepted history.
    /// The full metric and accumulated flow are committed together only after
    /// native equilibrium, material-domain, yield and independent work gates.
    /// Loads may be reversed or removed; this is not a ramp from virgin state.
    ///
    /// # Errors
    /// Requires one native load step and explicit finite nonnegative work/yield
    /// tolerances. Every refusal, including cancellation after equilibrium,
    /// returns no candidate and leaves `self` unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticTetStep<'a, L>, FinitePlasticTetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        clippy::type_complexity
    )]
    fn advance_loaded(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(FinitePlasticTetStep<'a, L>, Option<AcceptedMotion>), FinitePlasticTetError> {
        checkpoint(cx)?;
        validate_controls(settings, controls)?;
        for (index, state) in self.state.element_states.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            self.material
                .at(index)
                .validate_committed(
                    &self.state.element_deformation[index],
                    state,
                    controls.yield_tolerance_pa,
                )
                .map_err(|source| FinitePlasticTetError::Material {
                    element: index,
                    source,
                })?;
        }
        let frozen: Vec<_> = self
            .state
            .element_states
            .iter()
            .enumerate()
            .map(|(index, accepted)| FrozenPoint {
                law: self.material.at(index),
                accepted,
            })
            .collect();
        let (equilibrium, force) = self.state.solve_endpoint(
            HyperelasticLawField::PerElement(&frozen),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            inertia,
            cx,
        )?;
        let new_u = equilibrium.displacement_m.as_flattened();
        let motion = inertia
            .map(|inertia| inertia.motion(new_u, cx))
            .transpose()
            .map_err(super::inertia_error)?;
        let mut element_states = Vec::with_capacity(self.state.elements.len());
        let mut element_deformation = Vec::with_capacity(self.state.elements.len());
        let mut element_first_piola_pa = Vec::with_capacity(self.state.elements.len());
        let mut energy = FinitePlasticTetEnergy {
            stored_before_j: self.state.stored_energy_j,
            stored_after_j: 0.0,
            plastic_dissipation_j: 0.0,
            creep_dissipation_j: 0.0,
            flow_potential_j: 0.0,
            incremental_potential_j: 0.0,
            dead_load_endpoint_work_j: 0.0,
            support_endpoint_work_j: 0.0,
            internal_endpoint_work_j: 0.0,
            equilibrium_work_residual_j: 0.0,
            endpoint_work_defect_j: 0.0,
        };
        for (index, element) in self.state.elements.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let law = self.material.at(index);
            let f = HyperTetProblem::<FrozenPoint<'_, L>>::deformation_gradient(element, new_u);
            let refused = |source| FinitePlasticTetError::Material {
                element: index,
                source,
            };
            let response = law
                .response(&f, &self.state.element_states[index])
                .map_err(refused)?;
            law.validate_committed(&f, &response.state, controls.yield_tolerance_pa)
                .map_err(refused)?;
            // Independent stored-energy evaluation uses the proposed committed
            // metric without a second plastic update or a surrogate elastic law.
            let stored = law.stored_energy(&f, &response.state).map_err(refused)?;
            energy_gate(
                "constitutive stored-energy consistency",
                element.volume * (stored - response.stored_energy_j_m3),
                controls.equilibrium_work_tolerance_j,
            )?;
            energy.stored_after_j = element.volume.mul_add(stored, energy.stored_after_j);
            energy.plastic_dissipation_j = element.volume.mul_add(
                response.plastic_dissipation_j_m3,
                energy.plastic_dissipation_j,
            );
            energy.creep_dissipation_j = element
                .volume
                .mul_add(response.creep_dissipation_j_m3, energy.creep_dissipation_j);
            energy.flow_potential_j = element
                .volume
                .mul_add(response.flow_potential_j_m3, energy.flow_potential_j);
            energy.incremental_potential_j = element.volume.mul_add(
                response.incremental_potential_j_m3,
                energy.incremental_potential_j,
            );
            element_deformation.push(f);
            element_first_piola_pa.push(response.first_piola_pa);
            element_states.push(response.state);
        }
        let work = self.state.endpoint_work(
            &equilibrium,
            &force,
            &element_deformation,
            &element_first_piola_pa,
            motion.as_ref().map(|(_, account)| account.inertial_work_j),
            cx,
        )?;
        energy.internal_endpoint_work_j = work.internal_j;
        energy.dead_load_endpoint_work_j = work.dead_load_j;
        energy.support_endpoint_work_j = work.support_j;
        energy.equilibrium_work_residual_j = work.residual_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j
            - (energy.stored_after_j - energy.stored_before_j)
            - energy.inelastic_dissipation_j();
        energy_gate(
            "assembled incremental-potential consistency",
            equilibrium.strain_energy_j - energy.incremental_potential_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "storage/flow-potential consistency",
            energy.incremental_potential_j - energy.stored_after_j - energy.flow_potential_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "support/dead-load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "signed endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if !energy.stored_after_j.is_finite()
            || energy.stored_after_j < 0.0
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
            || !energy.creep_dissipation_j.is_finite()
            || energy.creep_dissipation_j < 0.0
            || !energy.flow_potential_j.is_finite()
            || energy.flow_potential_j < 0.0
        {
            return Err(invalid(
                "finite nonnegative storage, inelastic heat and flow potential required",
            ));
        }
        checkpoint(cx)?;
        let solid = Self {
            material: self.material,
            duration_s: self.duration_s,
            state: FinitePlasticTetState {
                nodes_m: self.state.nodes_m,
                tetrahedra: self.state.tetrahedra,
                elements: Arc::clone(&self.state.elements),
                budget: self.state.budget,
                displacement_m: equilibrium.displacement_m.clone(),
                element_deformation,
                element_states,
                element_first_piola_pa,
                stored_energy_j: energy.stored_after_j,
            },
        };
        Ok((
            FinitePlasticTetStep {
                solid,
                equilibrium,
                energy,
            },
            motion,
        ))
    }
}
