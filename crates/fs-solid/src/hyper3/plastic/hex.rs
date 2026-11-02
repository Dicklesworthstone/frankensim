//! Retained finite thermoplasticity on eight-node hexahedra.
//!
//! Eight fully integrated material points retain independent native histories.
//! The optional mean-dilatation map uses `Fbar_q=(Jbar/J_q)^(1/3) F_q`, where
//! `Jbar` is the reference-volume-weighted element mean. The element potential
//! is the sum of the ORIGINAL native potentials evaluated at these assumed
//! gradients. Its complete chain rule supplies force and tangent, including
//! geometric and cross-quadrature terms. The compatible physical gradients
//! remain independently available and cannot hide local inversion.
//!
//! Temperature is constant within each cell during an endpoint proposal. Every
//! Newton trial uses the same accepted history and actual physical interval;
//! the existing native J2, Perzyna, Norton and simultaneous Norton/J2 selections
//! keep their source family, kinetics and thermal energy conventions intact.
//! Plastic and creep losses plus reversible exchange are caloric heat. The
//! difference between nodal endpoint work and `sum w Pbar:delta Fbar` is a
//! signed projection integration defect, added to the native endpoint defect
//! and NEVER deposited as heat.
//!
//! This is a finite-bulk, displacement-only representation of a condensed
//! mean-dilatation formulation, with fully integrated deviatoric response.
//! It supplies no fully incompressible inf-sup guarantee, general shear-locking
//! remedy, follower/contact loads, evolving conductivity or heat solve.
//! [`dynamics`] adds consistent reference inertia on the same material solve.

use std::collections::BTreeSet;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::plastic::finite::FiniteJ2State;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteThermalJ2Material, FrozenThermalFiniteJ2,
};
use fs_sparse::{Coo, Csr};

use super::FinitePlasticStepControls;
use super::thermal::FiniteThermoplasticElementEnergy;
use crate::plastic3::hex::HexAssemblyBudget;
use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// Reference-mass momentum with complete thermal and inelastic history.
pub mod dynamics;
mod element;
use element::PreparedHex;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

type AcceptedHexMotion = (Vec<[f64; 3]>, crate::inertia::InertialEnergy);

/// Compatible full integration or the finite-strain assumed volume map.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteHexFormulation {
    Standard,
    MeanDilatation,
}

/// Preserve solver and point-material diagnoses; all refusals are transactional.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteThermoplasticHexError {
    Solid(SmallStrainTetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2ThermalError,
    },
    Energy {
        quantity: &'static str,
        error_j: f64,
        maximum_j: f64,
    },
}

impl From<SmallStrainTetError> for FiniteThermoplasticHexError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error)
    }
}
impl core::fmt::Display for FiniteThermoplasticHexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material {
                element,
                point,
                source,
            } => {
                write!(f, "finite thermal hex {element}, point {point}: {source}")
            }
            Self::Energy {
                quantity,
                error_j,
                maximum_j,
            } => {
                write!(
                    f,
                    "finite thermal hex {quantity}: {error_j:e} J exceeds {maximum_j:e} J"
                )
            }
        }
    }
}
impl std::error::Error for FiniteThermoplasticHexError {}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}
fn energy_gate(
    quantity: &'static str,
    error_j: f64,
    maximum_j: f64,
) -> Result<(), FiniteThermoplasticHexError> {
    if !error_j.is_finite() || error_j.abs() > maximum_j {
        return Err(FiniteThermoplasticHexError::Energy {
            quantity,
            error_j,
            maximum_j,
        });
    }
    Ok(())
}

/// One frozen-history trial. Forces are conjugate to actual nodal displacement;
/// the optional tangent includes every derivative of the assumed deformation.
#[derive(Debug)]
pub struct FiniteThermoplasticHexEvaluation {
    pub internal_force_n: Vec<[f64; 3]>,
    pub tangent_n_m: Option<Csr>,
    pub incremental_potential_j: f64,
}

/// Whole-reference-mesh accounts; caloric enthalpy belongs to a separate owner.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteThermoplasticHexEnergy {
    pub free_energy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub plastic_dissipation_j: f64,
    pub creep_dissipation_j: f64,
    pub reversible_heat_j: f64,
    pub caloric_heat_j: f64,
    pub coefficient_storage_change_j: f64,
    pub thermoelastic_storage_change_j: f64,
    /// Native material-coordinate integration defect, excluding projection.
    pub material_endpoint_work_defect_j: f64,
    /// Nodal endpoint work minus native work in the assumed Fbar coordinates.
    pub projection_work_defect_j: f64,
    /// Sum of native and projection defects; signed numerical error, never heat.
    pub endpoint_work_defect_j: f64,
    pub incremental_potential_j: f64,
    /// Sum of native Pbar:delta Fbar over original quadrature volume.
    pub material_endpoint_work_j: f64,
    /// Actual nodal internal force dotted with the nodal displacement increment.
    pub internal_endpoint_work_j: f64,
    /// Optional momentum work, separately resolved into kinetic change and
    /// numerical loss by the dynamic owner. Zero for a static increment.
    pub inertial_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    pub balance_error_j: f64,
}

impl FiniteThermoplasticHexEnergy {
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

/// Original source families, prepared geometry and complete accepted histories.
/// Arrays of quadrature data use element-major order, eight points per element.
#[derive(Debug)]
pub struct FiniteThermoplasticHexSolid<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    nodes_m: &'a [[f64; 3]],
    hexahedra: &'a [[usize; 8]],
    elements: Arc<[PreparedHex]>,
    families: Arc<[L]>,
    formulation: FiniteHexFormulation,
    budget: HexAssemblyBudget,
    temperatures_k: Vec<f64>,
    displacement_m: Vec<[f64; 3]>,
    physical_deformation: Vec<[f64; 9]>,
    assumed_deformation: Vec<[f64; 9]>,
    states: Vec<FiniteJ2State>,
    first_piola_pa: Vec<[f64; 9]>,
    free_energy_j: f64,
    internal_energy_j: f64,
}

/// Fully checked candidate; the borrowed previous solid remains unchanged.
#[derive(Debug)]
pub struct FiniteThermoplasticHexStep<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    pub solid: FiniteThermoplasticHexSolid<'a, L>,
    pub equilibrium: IncrementRecord,
    /// Native assumed-coordinate accounts, per original reference volume.
    pub quadrature_energy: Vec<FiniteThermoplasticElementEnergy>,
    pub energy: FiniteThermoplasticHexEnergy,
}

struct FrozenCell {
    law: FrozenThermalFiniteJ2,
    minimum_mechanical_jacobian: f64,
}

struct HexOperator<'a> {
    elements: &'a [PreparedHex],
    frozen: &'a [FrozenCell],
    initial: &'a [FiniteJ2State],
    formulation: FiniteHexFormulation,
    dofs: usize,
    inertia: Option<&'a crate::inertia::BackwardEulerInertia<'a>>,
}

// Keep the solver's flat vector throughout Newton; reshape only for the public
// inspection result, avoiding two extra force allocations on every trial.
struct AssembledHex {
    internal_force_n: Vec<f64>,
    tangent_n_m: Option<Csr>,
    incremental_potential_j: f64,
}

fn local_displacement(element: &PreparedHex, u: &[f64]) -> [f64; 24] {
    core::array::from_fn(|i| u[3 * element.nodes[i / 3] + i % 3])
}

impl HexOperator<'_> {
    fn assemble(
        &self,
        states: &[FiniteJ2State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(AssembledHex, f64), SmallStrainTetError> {
        checkpoint(cx)?;
        if states.len() != 8 * self.elements.len()
            || u.len() != self.dofs
            || u.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("finite hex trial displacement/history mismatch"));
        }
        let mut force = vec![0.0; self.dofs];
        let mut potential = 0.0;
        let mut force_scale = 0.0;
        let mut coo = need_tangent.then(|| Coo::new(self.dofs, self.dofs));
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let cell = &self.frozen[e];
            let result = element.evaluate(
                &cell.law,
                &states[8 * e..8 * e + 8],
                &local_displacement(element, u),
                self.formulation == FiniteHexFormulation::MeanDilatation,
                cell.minimum_mechanical_jacobian,
                need_tangent,
                cx,
            )?;
            potential += result.incremental_potential_j;
            for (i, value) in result.local_force_n.iter().enumerate() {
                let dof = 3 * element.nodes[i / 3] + i % 3;
                force[dof] += value;
                if prescribed.is_none_or(|fixed| !fixed[dof]) {
                    force_scale += value.abs();
                }
            }
            if let Some(matrix) = &mut coo {
                let tangent = result
                    .tangent_n_m
                    .expect("full element request returns tangent");
                for (i, row) in tangent.iter().enumerate() {
                    for (j, value) in row.iter().enumerate() {
                        if *value != 0.0 {
                            matrix.push(
                                3 * element.nodes[i / 3] + i % 3,
                                3 * element.nodes[j / 3] + j % 3,
                                *value,
                            );
                        }
                    }
                }
            }
        }
        if let (Some(matrix), Some(inertia)) = (&mut coo, self.inertia) {
            inertia
                .add_tangent(matrix, cx)
                .map_err(SmallStrainTetError::from)?;
        }
        let tangent_n_m = coo.map(|matrix| matrix.assemble());
        if !potential.is_finite()
            || !force_scale.is_finite()
            || force.iter().any(|x| !x.is_finite())
            || tangent_n_m.as_ref().is_some_and(|matrix| {
                (0..matrix.nrows()).any(|row| matrix.row(row).1.iter().any(|x| !x.is_finite()))
            })
        {
            return Err(invalid(
                "finite hex assembled force, potential or tangent is not finite",
            ));
        }
        checkpoint(cx)?;
        Ok((
            AssembledHex {
                internal_force_n: force,
                tangent_n_m,
                incremental_potential_j: potential,
            },
            force_scale,
        ))
    }
}

impl IncrementOperator for HexOperator<'_> {
    type State = FiniteJ2State;
    fn initial_states(&self) -> Vec<Self::State> {
        self.initial.to_vec()
    }
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let (result, scale) = self.assemble(states, u, prescribed, false, cx)?;
        Ok((result.internal_force_n, scale))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        Ok(self
            .assemble(states, u, None, true, cx)?
            .0
            .tangent_n_m
            .expect("full tangent request"))
    }
    fn tangent_symmetric(&self) -> bool {
        // The native finite laws and the complete geometric projection are
        // derivatives of the same scalar incremental potential.
        true
    }
    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        _final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        let mut changed = 0;
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let cell = &self.frozen[e];
            let theta = cell.law.thermal_stretch();
            let kinematics = element.kinematics(
                &local_displacement(element, u),
                self.formulation == FiniteHexFormulation::MeanDilatation,
                cell.minimum_mechanical_jacobian * theta * theta * theta,
                cx,
            )?;
            let mut evolving = false;
            for q in 0..8 {
                let old = &states[8 * e + q];
                let response = cell
                    .law
                    .response(&kinematics.assumed_deformation[q], old)
                    .map_err(|error| invalid(format!("finite hex {e} point {q}: {error}")))?;
                evolving |= response.state != *old;
                states[8 * e + q] = response.state;
            }
            changed += usize::from(evolving);
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        // Initial/tangent/commit errors still propagate. Only a rejected Newton
        // proposal may recover by shortening; cancellation is never swallowed.
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticHexSolid<'a, L> {
    /// Admit a virgin reference mesh with one complete source family per cell.
    /// Initial temperature must equal each family's declared reference.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        hexahedra: &'a [[usize; 8]],
        families: Vec<L>,
        initial_temperatures_k: &[f64],
        formulation: FiniteHexFormulation,
        budget: HexAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticHexError> {
        checkpoint(cx)?;
        let count = hexahedra.len();
        let points = count
            .checked_mul(8)
            .ok_or_else(|| invalid("hex point count overflow"))?;
        let entries = count
            .checked_mul(24 * 24)
            .ok_or_else(|| invalid("hex matrix count overflow"))?;
        nodes_m
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("hex DOF count overflow"))?;
        if nodes_m.is_empty()
            || count == 0
            || families.len() != count
            || initial_temperatures_k.len() != count
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
            || nodes_m.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(invalid("nonempty finite hex mesh, one source/T per cell and positive geometry bound required").into());
        }
        for (quantity, value, limit) in [
            ("nodes", nodes_m.len(), budget.maximum_nodes),
            ("elements", count, budget.maximum_elements),
            (
                "quadrature points",
                points,
                budget.maximum_quadrature_points,
            ),
            ("tangent entries", entries, budget.maximum_tangent_entries),
        ] {
            if value > limit {
                return Err(invalid(format!(
                    "finite hex {quantity} {value} exceeds budget {limit}"
                ))
                .into());
            }
        }
        let mut elements = Vec::with_capacity(count);
        let mut states = Vec::with_capacity(points);
        let mut seen = BTreeSet::new();
        let mut used = vec![false; nodes_m.len()];
        for (e, connectivity) in hexahedra.iter().enumerate() {
            checkpoint(cx)?;
            let mut key = *connectivity;
            key.sort_unstable();
            if key[7] >= nodes_m.len()
                || key.windows(2).any(|pair| pair[0] == pair[1])
                || !seen.insert(key)
            {
                return Err(invalid(format!(
                    "finite hex {e} has invalid or duplicate connectivity"
                ))
                .into());
            }
            for &node in connectivity {
                used[node] = true;
            }
            let family = families[e].thermal_finite_j2();
            if initial_temperatures_k[e] != family.parameters().reference_temperature_k {
                return Err(invalid(
                    "virgin finite hex temperature must equal original source reference",
                )
                .into());
            }
            family
                .at_temperature_for_inspection(initial_temperatures_k[e])
                .map_err(|source| FiniteThermoplasticHexError::Material {
                    element: e,
                    point: 0,
                    source,
                })?;
            elements.push(PreparedHex::new(
                nodes_m,
                *connectivity,
                e,
                budget.minimum_scaled_jacobian,
                cx,
            )?);
            states.extend((0..8).map(|_| family.initial_state()));
        }
        if used.iter().any(|value| !value) {
            return Err(invalid("finite hex mesh contains unused vertices").into());
        }
        Ok(Self {
            nodes_m,
            hexahedra,
            elements: Arc::from(elements),
            families: Arc::from(families),
            formulation,
            budget,
            temperatures_k: initial_temperatures_k.to_vec(),
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            physical_deformation: vec![IDENTITY; points],
            assumed_deformation: vec![IDENTITY; points],
            states,
            first_piola_pa: vec![[0.0; 9]; points],
            free_energy_j: 0.0,
            internal_energy_j: 0.0,
        })
    }

    pub fn families(&self) -> &[L] {
        &self.families
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    pub fn quadrature_deformation(&self) -> &[[f64; 9]] {
        &self.physical_deformation
    }
    pub fn quadrature_assumed_deformation(&self) -> &[[f64; 9]] {
        &self.assumed_deformation
    }
    pub fn quadrature_states(&self) -> &[FiniteJ2State] {
        &self.states
    }
    /// Native Piola stress conjugate to ASSUMED Fbar, not a physical nodal stress.
    pub fn quadrature_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    pub fn free_energy_j(&self) -> f64 {
        self.free_energy_j
    }
    pub fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub fn formulation(&self) -> FiniteHexFormulation {
        self.formulation
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements.get(element).map(|e| e.reference_volume_m3)
    }
    pub fn quadrature_volume_m3(&self, element: usize, point: usize) -> Option<f64> {
        self.elements
            .get(element)?
            .points
            .get(point)
            .map(|p| p.volume_m3)
    }

    fn freeze(
        &self,
        temperatures_k: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<FrozenCell>, FiniteThermoplasticHexError> {
        checkpoint(cx)?;
        if temperatures_k.len() != self.families.len() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid(
                "one finite endpoint temperature per hex and a positive physical interval required",
            )
            .into());
        }
        self.families
            .iter()
            .enumerate()
            .map(|(e, family)| {
                checkpoint(cx)?;
                let family = family.thermal_finite_j2();
                Ok(FrozenCell {
                    law: family
                        .at_temperature_for_step(temperatures_k[e], dt_s)
                        .map_err(|source| FiniteThermoplasticHexError::Material {
                            element: e,
                            point: 0,
                            source,
                        })?,
                    minimum_mechanical_jacobian: family.limits().minimum_jacobian,
                })
            })
            .collect()
    }

    /// Inspect a proposed displacement using immutable accepted histories.
    /// Calling this repeatedly never consumes time or commits material flow.
    pub fn evaluate_trial(
        &self,
        new_temperatures_k: &[f64],
        trial_displacement_m: &[[f64; 3]],
        dt_s: f64,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexEvaluation, FiniteThermoplasticHexError> {
        let frozen = self.freeze(new_temperatures_k, dt_s, cx)?;
        let operator = HexOperator {
            elements: &self.elements,
            frozen: &frozen,
            initial: &self.states,
            formulation: self.formulation,
            dofs: 3 * self.nodes_m.len(),
            inertia: None,
        };
        let result = operator
            .assemble(
                &self.states,
                trial_displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )?
            .0;
        Ok(FiniteThermoplasticHexEvaluation {
            internal_force_n: result
                .internal_force_n
                .chunks_exact(3)
                .map(|x| [x[0], x[1], x[2]])
                .collect(),
            tangent_n_m: result.tangent_n_m,
            incremental_potential_j: result.incremental_potential_j,
        })
    }

    /// One actual physical endpoint at supplied cell temperatures and dead loads.
    /// Supports and nodal loads are absolute endpoint values. The formulation
    /// and source family cannot change after construction.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexStep<'a, L>, FiniteThermoplasticHexError> {
        Ok(self
            .advance_with_inertia(
                new_temperatures_k,
                prescribed_m,
                nodal_forces_n,
                dt_s,
                settings,
                controls,
                None,
                cx,
            )?
            .0)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_with_inertia(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<
        (FiniteThermoplasticHexStep<'a, L>, Option<AcceptedHexMotion>),
        FiniteThermoplasticHexError,
    > {
        checkpoint(cx)?;
        let dofs = 3 * self.nodes_m.len();
        if settings.max_newton_iterations == 0
            || settings.max_linear_iterations == 0
            || self.budget.maximum_load_steps == 0
            || !settings.relative_tolerance.is_finite()
            || settings.relative_tolerance <= 0.0
            || !settings.linear_tolerance.is_finite()
            || settings.linear_tolerance <= 0.0
            || [
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
                controls.yield_tolerance_pa,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "finite hex needs explicit finite solver and nonnegative work/yield controls",
            )
            .into());
        }
        let mut fixed = vec![false; dofs];
        for &(dof, value) in prescribed_m {
            if dof >= dofs || !value.is_finite() || fixed[dof] {
                return Err(invalid("invalid or repeated finite hex prescribed DOF").into());
            }
            fixed[dof] = true;
        }
        if dofs - prescribed_m.len() > self.budget.maximum_free_dofs {
            return Err(invalid("finite hex free displacement count exceeds budget").into());
        }
        let mut external = vec![0.0; dofs];
        for &(dof, value) in nodal_forces_n {
            if dof >= dofs || !value.is_finite() {
                return Err(invalid("invalid finite hex dead nodal force").into());
            }
            external[dof] += value;
        }
        if external.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite hex summed nodal force is not finite").into());
        }
        let frozen = self.freeze(new_temperatures_k, dt_s, cx)?;
        for (e, family) in self.families.iter().enumerate() {
            checkpoint(cx)?;
            let old_law = family
                .thermal_finite_j2()
                .at_temperature_for_inspection(self.temperatures_k[e])
                .map_err(|source| FiniteThermoplasticHexError::Material {
                    element: e,
                    point: 0,
                    source,
                })?;
            for q in 0..8 {
                old_law
                    .validate_committed(
                        &self.assumed_deformation[8 * e + q],
                        &self.states[8 * e + q],
                        controls.yield_tolerance_pa,
                    )
                    .map_err(|source| FiniteThermoplasticHexError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
            }
        }
        let operator = HexOperator {
            elements: &self.elements,
            frozen: &frozen,
            initial: &self.states,
            formulation: self.formulation,
            dofs,
            inertia,
        };
        let solved = increment::solve(
            &operator,
            IncrementControls {
                dofs,
                prescribed_m,
                nodal_forces_n,
                load_path: &[1.0],
                settings,
            },
            Some((self.displacement_m.as_flattened(), &self.states)),
            inertia,
            cx,
        )?;
        let new_u = &solved.displacement_m;
        let motion = inertia
            .map(|owner| owner.motion(new_u, cx))
            .transpose()
            .map_err(SmallStrainTetError::from)?;
        let assembled = operator.assemble(&self.states, new_u, None, false, cx)?.0;
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .expect("one physical endpoint");
        let mut quadrature_energy = Vec::with_capacity(self.states.len());
        let mut physical = Vec::with_capacity(self.states.len());
        let mut assumed = Vec::with_capacity(self.states.len());
        let mut stress = Vec::with_capacity(self.states.len());
        let mut energy = FiniteThermoplasticHexEnergy {
            inertial_endpoint_work_j: motion
                .as_ref()
                .map_or(0.0, |(_, energy)| energy.inertial_work_j),
            ..FiniteThermoplasticHexEnergy::default()
        };
        let (mut old_free, mut new_free, mut old_internal, mut new_internal, mut flow_potential) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        let tolerance = controls.equilibrium_work_tolerance_j;
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let theta = frozen[e].law.thermal_stretch();
            let coordinates = element.kinematics(
                &local_displacement(element, new_u),
                self.formulation == FiniteHexFormulation::MeanDilatation,
                frozen[e].minimum_mechanical_jacobian * theta * theta * theta,
                cx,
            )?;
            let family = self.families[e].thermal_finite_j2();
            for q in 0..8 {
                let index = 8 * e + q;
                let value = family
                    .checked_increment_for_step(
                        &self.assumed_deformation[index],
                        &self.states[index],
                        self.temperatures_k[e],
                        &coordinates.assumed_deformation[q],
                        new_temperatures_k[e],
                        controls.yield_tolerance_pa,
                        dt_s,
                    )
                    .map_err(|source| FiniteThermoplasticHexError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
                if value.mechanical.state != solved.point_states[index] {
                    return Err(invalid(
                        "finite hex endpoint audit differs from frozen-history accepted return",
                    )
                    .into());
                }
                let account = FiniteThermoplasticElementEnergy::from_increment(&value);
                let weight = element.points[q].volume_m3;
                energy_gate(
                    "native point balance",
                    weight * account.balance_error_j_m3,
                    tolerance,
                )?;
                energy_gate(
                    "native stored energy",
                    weight * (account.free_after_j_m3 - value.mechanical.stored_energy_j_m3),
                    tolerance,
                )?;
                old_free = weight.mul_add(account.free_before_j_m3, old_free);
                new_free = weight.mul_add(account.free_after_j_m3, new_free);
                old_internal = weight.mul_add(account.internal_before_j_m3, old_internal);
                new_internal = weight.mul_add(account.internal_after_j_m3, new_internal);
                energy.free_energy_change_j = weight.mul_add(
                    account.free_after_j_m3 - account.free_before_j_m3,
                    energy.free_energy_change_j,
                );
                energy.internal_energy_change_j = weight.mul_add(
                    account.internal_energy_change_j_m3(),
                    energy.internal_energy_change_j,
                );
                energy.plastic_dissipation_j = weight.mul_add(
                    account.plastic_dissipation_j_m3,
                    energy.plastic_dissipation_j,
                );
                energy.creep_dissipation_j =
                    weight.mul_add(account.creep_dissipation_j_m3, energy.creep_dissipation_j);
                energy.reversible_heat_j =
                    weight.mul_add(account.reversible_heat_j_m3, energy.reversible_heat_j);
                energy.caloric_heat_j =
                    weight.mul_add(account.caloric_heat_j_m3(), energy.caloric_heat_j);
                energy.coefficient_storage_change_j = weight.mul_add(
                    account.coefficient_storage_change_j_m3,
                    energy.coefficient_storage_change_j,
                );
                energy.thermoelastic_storage_change_j = weight.mul_add(
                    account.thermoelastic_storage_change_j_m3,
                    energy.thermoelastic_storage_change_j,
                );
                energy.material_endpoint_work_defect_j = weight.mul_add(
                    account.endpoint_work_defect_j_m3,
                    energy.material_endpoint_work_defect_j,
                );
                energy.material_endpoint_work_j = weight.mul_add(
                    account.internal_endpoint_work_j_m3,
                    energy.material_endpoint_work_j,
                );
                energy.incremental_potential_j = weight.mul_add(
                    account.incremental_potential_j_m3,
                    energy.incremental_potential_j,
                );
                let flow = if let Some(kinetics) = family.overstress() {
                    let p = family.parameters();
                    let yield_pa = p.slopes.yield_stress_pa_k.mul_add(
                        new_temperatures_k[e] - p.reference_temperature_k,
                        p.reference.yield_stress_pa,
                    );
                    let yield_work = yield_pa * value.mechanical.plastic_increment;
                    yield_work
                        + (1.0 / (1.0 + 1.0 / kinetics.exponent))
                            * (account.inelastic_dissipation_j_m3() - yield_work)
                } else {
                    account.inelastic_dissipation_j_m3()
                };
                flow_potential = weight.mul_add(flow, flow_potential);
                physical.push(coordinates.physical_deformation[q]);
                assumed.push(coordinates.assumed_deformation[q]);
                stress.push(value.mechanical.first_piola_pa);
                quadrature_energy.push(account);
            }
        }
        let old_u = self.displacement_m.as_flattened();
        let force = &assembled.internal_force_n;
        for i in 0..dofs {
            let delta = new_u[i] - old_u[i];
            energy.internal_endpoint_work_j =
                force[i].mul_add(delta, energy.internal_endpoint_work_j);
            energy.dead_load_endpoint_work_j =
                external[i].mul_add(delta, energy.dead_load_endpoint_work_j);
        }
        for (&(dof, _), reaction) in prescribed_m.iter().zip(&equilibrium.reactions_n) {
            energy.support_endpoint_work_j =
                reaction.mul_add(new_u[dof] - old_u[dof], energy.support_endpoint_work_j);
        }
        energy.projection_work_defect_j =
            energy.internal_endpoint_work_j - energy.material_endpoint_work_j;
        energy.endpoint_work_defect_j =
            energy.material_endpoint_work_defect_j + energy.projection_work_defect_j;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j
            - energy.inertial_endpoint_work_j;
        energy.balance_error_j = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j
            - energy.internal_energy_change_j
            - energy.caloric_heat_j
            - energy.endpoint_work_defect_j
            - energy.inertial_endpoint_work_j;
        for (quantity, value) in [
            ("accepted free energy", old_free - self.free_energy_j),
            (
                "accepted internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "assembled incremental potential",
                assembled.incremental_potential_j - energy.incremental_potential_j,
            ),
            (
                "stored energy and integrated flow potential",
                energy.incremental_potential_j - new_free - flow_potential,
            ),
            (
                "free-energy increment",
                new_free - old_free - energy.free_energy_change_j,
            ),
            (
                "internal-energy increment",
                new_internal - old_internal - energy.internal_energy_change_j,
            ),
            (
                "caloric decomposition",
                energy.caloric_heat_j - energy.inelastic_dissipation_j() - energy.reversible_heat_j,
            ),
            (
                "support and dead-load work",
                energy.equilibrium_work_residual_j,
            ),
            ("total endpoint balance", energy.balance_error_j),
        ] {
            energy_gate(quantity, value, tolerance)?;
        }
        energy_gate(
            "signed total endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if !new_free.is_finite()
            || new_free < 0.0
            || !new_internal.is_finite()
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
            || !energy.creep_dissipation_j.is_finite()
            || energy.creep_dissipation_j < 0.0
            || [
                energy.coefficient_storage_change_j,
                energy.thermoelastic_storage_change_j,
                energy.reversible_heat_j,
                energy.caloric_heat_j,
                energy.material_endpoint_work_defect_j,
                energy.projection_work_defect_j,
            ]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "finite hex storage, heat and numerical accounts must be representable",
            )
            .into());
        }
        checkpoint(cx)?;
        Ok((
            FiniteThermoplasticHexStep {
                solid: Self {
                    nodes_m: self.nodes_m,
                    hexahedra: self.hexahedra,
                    elements: Arc::clone(&self.elements),
                    families: Arc::clone(&self.families),
                    formulation: self.formulation,
                    budget: self.budget,
                    temperatures_k: new_temperatures_k.to_vec(),
                    displacement_m: new_u.chunks_exact(3).map(|x| [x[0], x[1], x[2]]).collect(),
                    physical_deformation: physical,
                    assumed_deformation: assumed,
                    states: solved.point_states,
                    first_piola_pa: stress,
                    free_energy_j: new_free,
                    internal_energy_j: new_internal,
                },
                equilibrium,
                quadrature_energy,
                energy,
            },
            motion,
        ))
    }
}
