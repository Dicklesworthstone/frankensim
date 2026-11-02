//! Finite-strain viscoelastic equilibrium with retained native material memory.
//!
//! Each tetrahedron owns one `ViscoState`; its original law, accepted history
//! and actual time interval stay fixed through every Newton/backtracking trial.
//! The existing `HyperTetProblem` solves the original reduced incremental
//! potential. Only after equilibrium and work gates pass are displacement,
//! viscous metrics, stresses, physical storage and time published together.
//!
//! Native frozen-state loading work is distinct from nodal endpoint work.
//! Their signed integration discrepancy is never viscous heat. Isothermal
//! quasistatics, finite bulk modulus and P1 tetrahedra only: this path does not
//! mitigate volumetric locking or provide mixed incompressibility, inertia,
//! contact, evolving temperature or source-parameter interpolation.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::visco_hyper::{
    FrozenViscoHyperelastic, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic, ViscoState,
};

use super::{
    Element, HyperTetError, HyperTetProblem, HyperTetSettings, HyperTetSolution,
    HyperelasticLawField, checkpoint,
};
use crate::linear3::TetAssemblyBudget;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Original immutable laws, without averaging regional relaxation spectra.
#[derive(Debug, Clone, Copy)]
pub enum ViscoelasticLawField<'a> {
    Uniform(&'a ViscoHyperelastic),
    /// Exactly one law per tetrahedron in connectivity order.
    PerElement(&'a [ViscoHyperelastic]),
}

impl<'a> ViscoelasticLawField<'a> {
    fn at(self, element: usize) -> &'a ViscoHyperelastic {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
        }
    }
}

/// Absolute physical tolerances; no integration error is silently accepted as heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteViscoStepControls {
    /// Native material balance, assembled potential and support/load work [J].
    pub equilibrium_work_tolerance_j: f64,
    /// Maximum absolute signed endpoint-work quadrature defect [J]. Reduce the
    /// physical time/load increment when this gate refuses.
    pub maximum_endpoint_work_defect_j: f64,
}

/// Physical and algorithmic accounts for one accepted endpoint.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteViscoTetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Nonnegative native relaxation loss during this interval, not cumulative.
    pub dissipated_j: f64,
    /// Native loading work at the old frozen viscous state; not nodal work.
    pub native_loading_work_j: f64,
    /// Newton's reduced potential; neither physical storage nor heat.
    pub incremental_potential_j: f64,
    /// Native loading work minus stored-energy change and relaxation loss.
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Dead-load plus support work minus internal endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus physical storage change and native loss.
    /// Signed quadrature error, reported and gated separately from dissipation.
    pub endpoint_work_defect_j: f64,
}

/// Mesh, native material bindings and complete immutable accepted state.
#[derive(Debug)]
pub struct FiniteViscoTetSolid<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    material: ViscoelasticLawField<'a>,
    elements: Arc<Vec<Element>>,
    budget: TetAssemblyBudget,
    displacement_m: Vec<[f64; 3]>,
    states: Vec<ViscoState>,
    first_piola_pa: Vec<[f64; 9]>,
    stored_energy_j: f64,
    time_s: f64,
}

/// A complete accepted candidate. Dropping it leaves the previous solid intact.
#[derive(Debug)]
pub struct FiniteViscoTetStep<'a> {
    pub solid: FiniteViscoTetSolid<'a>,
    /// The native solver's `strain_energy_j` is the incremental potential here.
    pub equilibrium: HyperTetSolution,
    pub energy: FiniteViscoTetEnergy,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteViscoTetError {
    Native(HyperTetError),
    Material {
        element: usize,
        source: MaterialError,
    },
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<HyperTetError> for FiniteViscoTetError {
    fn from(error: HyperTetError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for FiniteViscoTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "viscoelastic tetrahedron {element}: {source}")
            }
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "viscoelastic {quantity}: {value_j:e} J exceeds {maximum_abs_j:e} J"
            ),
        }
    }
}
impl std::error::Error for FiniteViscoTetError {}

fn invalid(what: &str) -> FiniteViscoTetError {
    HyperTetError::InvalidInput { what: what.into() }.into()
}
fn energy_gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteViscoTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FiniteViscoTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

impl<'a> FiniteViscoTetSolid<'a> {
    /// Prepare the undeformed relaxed solid. Geometry is prepared once and
    /// shared by later candidates; every element retains its own native history.
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        material: ViscoelasticLawField<'a>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoTetError> {
        checkpoint(cx)?;
        if nodes_m.is_empty()
            || tetrahedra.is_empty()
            || nodes_m.len() > budget.maximum_nodes
            || tetrahedra.len() > budget.maximum_tetrahedra
        {
            return Err(invalid(
                "viscoelastic reference mesh is empty or exceeds its budget",
            ));
        }
        if let ViscoelasticLawField::PerElement(laws) = material
            && laws.len() != tetrahedra.len()
        {
            return Err(HyperTetError::LawCountMismatch {
                expected: tetrahedra.len(),
                got: laws.len(),
            }
            .into());
        }
        let dofs = nodes_m
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("viscoelastic displacement size overflow"))?;
        let mut states = Vec::with_capacity(tetrahedra.len());
        let mut first_piola_pa = Vec::with_capacity(tetrahedra.len());
        for element in 0..tetrahedra.len() {
            checkpoint(cx)?;
            let law = material.at(element);
            if law.branches.len() > MAX_FINITE_STRAIN_BRANCHES {
                return Err(invalid(
                    "viscoelastic branch count exceeds the native material cap",
                ));
            }
            let response = law
                .step_response(&law.initial_state(), &IDENTITY, 0.0)
                .map_err(|source| FiniteViscoTetError::Material { element, source })?;
            states.push(response.state);
            first_piola_pa.push(response.first_piola_pa);
        }
        let frozen = Self::freeze(material, &states, 0.0, cx)?;
        let prescribed: Vec<_> = (0..dofs).map(|dof| (dof, 0.0)).collect();
        let problem = HyperTetProblem {
            nodes_m,
            tetrahedra,
            material: HyperelasticLawField::PerElement(&frozen),
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget,
            settings: HyperTetSettings::default(),
        };
        let elements = problem.validate()?;
        let stored_energy_j = elements
            .iter()
            .zip(&states)
            .map(|(element, state)| element.volume * state.stored)
            .sum::<f64>();
        if !stored_energy_j.is_finite() {
            return Err(invalid("initial viscoelastic storage is unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(Self {
            nodes_m,
            tetrahedra,
            material,
            elements: Arc::new(elements),
            budget,
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            states,
            first_piola_pa,
            stored_energy_j,
            time_s: 0.0,
        })
    }

    fn freeze<'b>(
        material: ViscoelasticLawField<'b>,
        states: &'b [ViscoState],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<FrozenViscoHyperelastic<'b>>, FiniteViscoTetError> {
        states
            .iter()
            .enumerate()
            .map(|(element, state)| {
                checkpoint(cx)?;
                material
                    .at(element)
                    .freeze(state, dt_s)
                    .map_err(|source| FiniteViscoTetError::Material { element, source })
            })
            .collect()
    }

    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    pub fn element_states(&self) -> &[ViscoState] {
        &self.states
    }
    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn material(&self) -> ViscoelasticLawField<'a> {
        self.material
    }

    /// Solve one actual interval with absolute endpoint supports and dead loads.
    /// `dt_s = 0` is the original instantaneous elastic jump. Loading, unloading,
    /// stress relaxation and creep all use the same persistent material history.
    ///
    /// Exactly one native load endpoint is required. Neither the previous solid
    /// nor any history changes on failure, including a failed post-solve work
    /// gate. Numerical substepping is the caller's choice and advances time only
    /// through accepted candidates.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoTetStep<'a>, FiniteViscoTetError> {
        checkpoint(cx)?;
        let time_s = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s < 0.0
            || !time_s.is_finite()
            || settings.load_steps != 1
            || ![
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
            ]
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "one endpoint, nonnegative finite interval and finite work allowances required",
            ));
        }
        let frozen = Self::freeze(self.material, &self.states, dt_s, cx)?;
        let problem = HyperTetProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            material: HyperelasticLawField::PerElement(&frozen),
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
            None,
            cx,
        )?;
        let u = equilibrium.displacement_m.as_flattened();
        let old_u = self.displacement_m.as_flattened();
        let mut energy = FiniteViscoTetEnergy {
            stored_before_j: self.stored_energy_j,
            ..Default::default()
        };
        let mut states = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        for (index, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let deformation =
                HyperTetProblem::<FrozenViscoHyperelastic<'_>>::deformation_gradient(element, u);
            let response = self
                .material
                .at(index)
                .step_response(&self.states[index], &deformation, dt_s)
                .map_err(|source| FiniteViscoTetError::Material {
                    element: index,
                    source,
                })?;
            let w = element.volume;
            energy.stored_after_j += w * response.state.stored;
            energy.dissipated_j += w * response.dissipated_pa;
            energy.native_loading_work_j += w * response.loading_work_pa;
            energy.incremental_potential_j += w * response.incremental_potential_pa;
            let work_density = (0..9).fold(0.0, |sum, c| {
                response.first_piola_pa[c].mul_add(
                    deformation[c] - self.states[index].deformation_gradient[c],
                    sum,
                )
            });
            energy.internal_endpoint_work_j += w * work_density;
            first_piola_pa.push(response.first_piola_pa);
            states.push(response.state);
        }
        let external = problem.external(&self.elements, 1.0);
        energy.dead_load_endpoint_work_j = external
            .iter()
            .zip(u.iter().zip(old_u))
            .fold(0.0, |sum, (force, (after, before))| {
                force.mul_add(after - before, sum)
            });
        energy.support_endpoint_work_j = equilibrium
            .reactions_n
            .iter()
            .fold(0.0, |sum, (dof, force)| {
                force.mul_add(u[*dof] - old_u[*dof], sum)
            });
        let stored_change = energy.stored_after_j - energy.stored_before_j;
        energy.material_balance_error_j =
            energy.native_loading_work_j - stored_change - energy.dissipated_j;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j =
            energy.internal_endpoint_work_j - stored_change - energy.dissipated_j;
        energy_gate(
            "native material balance",
            energy.material_balance_error_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "assembled incremental potential",
            equilibrium.strain_energy_j - energy.incremental_potential_j,
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
            || !energy.dissipated_j.is_finite()
            || energy.dissipated_j < 0.0
        {
            return Err(invalid(
                "finite storage and nonnegative native relaxation loss required",
            ));
        }
        checkpoint(cx)?;
        let solid = Self {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            material: self.material,
            elements: Arc::clone(&self.elements),
            budget: self.budget,
            displacement_m: equilibrium.displacement_m.clone(),
            states,
            first_piola_pa,
            stored_energy_j: energy.stored_after_j,
            time_s,
        };
        Ok(FiniteViscoTetStep {
            solid,
            equilibrium,
            energy,
        })
    }
}
