//! Finite progressive ductility at all eight original H8 material points.
//!
//! The native H8 geometry, full integration and mean-dilatation pullback are
//! shared with the thermal plastic owner. The original finite ductile material
//! supplies first-Piola stress and its complete nonsymmetric derivative; no
//! scalar potential, symmetric tangent, hourglass stiffness or extra return
//! map is substituted. Residual probes use the material's stress-only path.
//!
//! Plastic dissipation and damage release are physical accounts. The native
//! material endpoint defect and the Fbar projection defect remain signed
//! numerical accounts. Acceptance bounds the sum of their absolute point/cell
//! contributions, so cancellation cannot hide an unresolved increment.
//!
//! Every point retains its exact immutable material and characteristic length
//! on the original mesh. Length is an explicitly authored modeling input;
//! this owner does not infer a crack normal or certify arbitrary-mesh fracture
//! objectivity. Temperature is a held initiation condition, not a thermal solve.
//! Physical cutback belongs to the caller. No erosion, fragment/contact handoff,
//! creep coupling, arclength, general shear-locking remedy or experimental
//! qualification is supplied. Physical Jacobian probes remain independent of
//! Fbar and do not prove global injectivity of a warped cell.

use core::cell::RefCell;
use core::fmt::Debug;
use std::collections::BTreeSet;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::ductile::evolution::{DifferentiableDuctileInitiationLaw, DuctileJ2Conditions};
use fs_material::plastic::finite::ductile::{
    FiniteDuctileEnergy, FiniteDuctileError, FiniteDuctileJ2, FiniteDuctileState,
};
use fs_sparse::{Coo, Csr};

use super::element::{HexKinematics, HexPointResponse, PreparedHex};
use super::{FiniteHexFormulation, IDENTITY, checkpoint, invalid, local_displacement};
use crate::ductile::{DuctileStepConditions, DuctileTemperatures};
use crate::hyper3::plastic::FinitePlasticStepControls;
use crate::plastic3::hex::HexAssemblyBudget;
use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

#[derive(Debug)]
pub enum FiniteDuctileHexError {
    Solid(SmallStrainTetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteDuctileError,
    },
    Energy {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}
impl From<SmallStrainTetError> for FiniteDuctileHexError {
    fn from(value: SmallStrainTetError) -> Self {
        Self::Solid(value)
    }
}
impl core::fmt::Display for FiniteDuctileHexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile hex refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileHexError {}

fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteDuctileHexError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        return Err(FiniteDuctileHexError::Energy {
            quantity,
            value_j,
            maximum_abs_j,
        });
    }
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let sum = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - sum) + value
        } else {
            (value - sum) + self.value
        };
        self.value = sum;
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}

/// Force and optional full nonsymmetric derivative of one frozen-history trial.
#[derive(Debug)]
pub struct FiniteDuctileHexEvaluation {
    pub internal_force_n: Vec<[f64; 3]>,
    pub tangent_n_m: Option<Csr>,
}

/// Whole-reference-mesh physical and numerical endpoint accounts, J.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteDuctileHexEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Integrated pointwise storage differences, avoiding large-total subtraction.
    pub stored_change_j: f64,
    pub plastic_dissipation_j: f64,
    pub damage_release_j: f64,
    /// Sum w_q Pbar_q : (Fbar_after_q - Fbar_before_q).
    pub material_endpoint_work_j: f64,
    /// Actual nodal internal force dotted with the nodal displacement change.
    pub internal_endpoint_work_j: f64,
    pub material_endpoint_work_defect_j: f64,
    /// Sum of cell nodal endpoint work minus native assumed-coordinate work.
    pub projection_work_defect_j: f64,
    pub endpoint_work_defect_j: f64,
    pub absolute_material_endpoint_work_defect_j: f64,
    pub absolute_projection_work_defect_j: f64,
    /// Sum of both absolute accounts; opposing errors cannot cancel this gate.
    pub absolute_endpoint_work_defect_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    pub balance_error_j: f64,
}

/// Immutable original mesh, point material bindings and complete accepted state.
#[derive(Debug)]
pub struct FiniteDuctileHexSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    nodes_m: &'g [[f64; 3]],
    hexahedra: &'g [[usize; 8]],
    elements: Arc<[PreparedHex]>,
    materials: &'m [&'m FiniteDuctileJ2<L>],
    minimum_jacobians: Arc<[f64]>,
    formulation: FiniteHexFormulation,
    budget: HexAssemblyBudget,
    displacement_m: Vec<[f64; 3]>,
    physical_deformation: Vec<[f64; 9]>,
    assumed_deformation: Vec<[f64; 9]>,
    states: Vec<FiniteDuctileState<'m, L>>,
    first_piola_pa: Vec<[f64; 9]>,
    stored_energy_j: f64,
}

/// One complete accepted interval. No input state is modified on refusal.
#[derive(Debug)]
pub struct FiniteDuctileHexStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub solid: FiniteDuctileHexSolid<'g, 'm, L>,
    pub equilibrium: IncrementRecord,
    /// Original native accounts per reference volume, element-major then q.
    pub quadrature_energy: Vec<FiniteDuctileEnergy>,
    pub energy: FiniteDuctileHexEnergy,
}

struct Operator<'s, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    elements: &'s [PreparedHex],
    materials: &'m [&'m FiniteDuctileJ2<L>],
    accepted: &'s [FiniteDuctileState<'m, L>],
    conditions: &'s [DuctileJ2Conditions],
    minimum_jacobians: &'s [f64],
    formulation: FiniteHexFormulation,
    dofs: usize,
    refusal: RefCell<Option<(usize, usize, FiniteDuctileError)>>,
}

impl<'m, L> Operator<'_, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn material_refusal(
        &self,
        element: usize,
        point: usize,
        source: FiniteDuctileError,
    ) -> SmallStrainTetError {
        let error = invalid(format!(
            "finite ductile hex {element}, point {point}: {source}"
        ));
        *self.refusal.borrow_mut() = Some((element, point, source));
        error
    }
    fn endpoint_refusal(&self, error: SmallStrainTetError) -> FiniteDuctileHexError {
        if matches!(error, SmallStrainTetError::Cancelled) {
            return error.into();
        }
        match self.refusal.borrow_mut().take() {
            Some((element, point, source)) => FiniteDuctileHexError::Material {
                element,
                point,
                source,
            },
            None => error.into(),
        }
    }
    fn kinematics(
        &self,
        element: usize,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<HexKinematics, SmallStrainTetError> {
        self.elements[element].kinematics(
            &local_displacement(&self.elements[element], u),
            self.formulation == FiniteHexFormulation::MeanDilatation,
            self.minimum_jacobians[element],
            cx,
        )
    }
    fn assemble(
        &self,
        states: &[FiniteDuctileState<'m, L>],
        u: &[f64],
        prescribed: Option<&[bool]>,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64, Option<Csr>), SmallStrainTetError> {
        *self.refusal.borrow_mut() = None;
        checkpoint(cx)?;
        if states.len() != 8 * self.elements.len()
            || u.len() != self.dofs
            || u.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "complete finite ductile H8 displacement/history fields required",
            ));
        }
        let mut force = vec![0.0; self.dofs];
        let mut scale = 0.0;
        let mut matrix = need_tangent.then(|| Coo::new(self.dofs, self.dofs));
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let kinematics = self.kinematics(e, u, cx)?;
            let response = element.evaluate_response(
                &kinematics,
                self.formulation == FiniteHexFormulation::MeanDilatation,
                need_tangent,
                |q, f, tangent| {
                    let index = 8 * e + q;
                    let refused = |source| self.material_refusal(e, q, source);
                    if tangent {
                        let value = self.materials[index]
                            .trial(f, &states[index], self.conditions[index])
                            .map_err(refused)?;
                        Ok(HexPointResponse {
                            first_piola_pa: value.first_piola_pa,
                            tangent_pa: Some(value.tangent_pa),
                        })
                    } else {
                        let stress = self.materials[index]
                            .stress(f, &states[index], self.conditions[index])
                            .map_err(refused)?;
                        Ok(HexPointResponse {
                            first_piola_pa: stress,
                            tangent_pa: None,
                        })
                    }
                },
                cx,
            )?;
            for (i, value) in response.local_force_n.iter().enumerate() {
                let row = 3 * element.nodes[i / 3] + i % 3;
                force[row] += value;
                if prescribed.is_none_or(|fixed| !fixed[row]) {
                    scale += value.abs();
                }
            }
            if let (Some(matrix), Some(tangent)) = (&mut matrix, &response.tangent_n_m) {
                for (i, row) in tangent.iter().enumerate() {
                    for (j, &value) in row.iter().enumerate() {
                        if value != 0.0 {
                            matrix.push(
                                3 * element.nodes[i / 3] + i % 3,
                                3 * element.nodes[j / 3] + j % 3,
                                value,
                            );
                        }
                    }
                }
            }
        }
        let matrix = matrix.map(|coo| coo.assemble());
        if !scale.is_finite()
            || force.iter().any(|x| !x.is_finite())
            || matrix.as_ref().is_some_and(|a| {
                (0..a.nrows()).any(|row| a.row(row).1.iter().any(|x| !x.is_finite()))
            })
        {
            return Err(invalid(
                "finite ductile H8 force or complete derivative overflow",
            ));
        }
        checkpoint(cx)?;
        Ok((force, scale, matrix))
    }
}

impl<'m, L> IncrementOperator for Operator<'_, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    type State = FiniteDuctileState<'m, L>;
    fn initial_states(&self) -> Vec<Self::State> {
        self.accepted.to_vec()
    }
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let (force, scale, _) = self.assemble(states, u, prescribed, false, cx)?;
        Ok((force, scale))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        Ok(self
            .assemble(states, u, None, true, cx)?
            .2
            .expect("requested complete H8 tangent"))
    }
    fn tangent_symmetric(&self) -> bool {
        false
    }
    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        _stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        *self.refusal.borrow_mut() = None;
        let mut changed = 0;
        for e in 0..self.elements.len() {
            checkpoint(cx)?;
            let kinematics = self.kinematics(e, u, cx)?;
            let mut evolving = false;
            for q in 0..8 {
                let i = 8 * e + q;
                let step = self.materials[i]
                    .trial(
                        &kinematics.assumed_deformation[q],
                        &states[i],
                        self.conditions[i],
                    )
                    .map_err(|source| self.material_refusal(e, q, source))?;
                evolving |= step.state.native_state() != states[i].native_state()
                    || step.state.damage() != states[i].damage();
                states[i] = step.state;
            }
            changed += usize::from(evolving);
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}

impl<'g, 'm, L> FiniteDuctileHexSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Bind exactly eight original material references per cell in the native
    /// `q=4*i+2*j+k` two-point Gauss order. Repeated references are supported;
    /// their histories remain independent. Each positive characteristic length
    /// is explicitly supplied by its immutable material, not guessed from mesh
    /// volume. This assignment and the reference mesh persist through all steps.
    pub fn try_new(
        nodes_m: &'g [[f64; 3]],
        hexahedra: &'g [[usize; 8]],
        materials: &'m [&'m FiniteDuctileJ2<L>],
        formulation: FiniteHexFormulation,
        budget: HexAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileHexError> {
        checkpoint(cx)?;
        let count = hexahedra.len();
        let points = count
            .checked_mul(8)
            .ok_or_else(|| invalid("H8 point count overflow"))?;
        let entries = count
            .checked_mul(24 * 24)
            .ok_or_else(|| invalid("H8 tangent count overflow"))?;
        nodes_m
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("H8 DOF count overflow"))?;
        if nodes_m.is_empty()
            || count == 0
            || materials.len() != points
            || nodes_m.iter().flatten().any(|x| !x.is_finite())
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
        {
            return Err(invalid("nonempty finite H8 mesh, eight original materials per cell and positive geometry bound required").into());
        }
        for (quantity, requested, maximum) in [
            ("nodes", nodes_m.len(), budget.maximum_nodes),
            ("elements", count, budget.maximum_elements),
            (
                "quadrature points",
                points,
                budget.maximum_quadrature_points,
            ),
            ("tangent entries", entries, budget.maximum_tangent_entries),
        ] {
            if requested > maximum {
                return Err(invalid(format!(
                    "finite ductile H8 {quantity} {requested} exceeds budget {maximum}"
                ))
                .into());
            }
        }
        let mut elements = Vec::with_capacity(count);
        let mut states = Vec::with_capacity(points);
        let mut minimum_jacobians = Vec::with_capacity(count);
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
                    "finite ductile H8 {e} has invalid or duplicate connectivity"
                ))
                .into());
            }
            for &node in connectivity {
                used[node] = true;
            }
            elements.push(PreparedHex::new(
                nodes_m,
                *connectivity,
                e,
                budget.minimum_scaled_jacobian,
                cx,
            )?);
            let mut minimum = 0.0_f64;
            for q in 0..8 {
                let model = materials[8 * e + q];
                let state = model.initial_state();
                model
                    .validate_committed(&state, model.settings().energy_tolerances.yield_pa)
                    .map_err(|source| FiniteDuctileHexError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
                minimum = minimum.max(model.native().limits().minimum_jacobian);
                states.push(state);
            }
            // Independent physical probes use the strictest original point
            // domain in this cell. Fbar cannot hide a local physical inversion.
            elements[e].kinematics(
                &[0.0; 24],
                formulation == FiniteHexFormulation::MeanDilatation,
                minimum,
                cx,
            )?;
            minimum_jacobians.push(minimum);
        }
        if used.iter().any(|used| !used) {
            return Err(invalid("finite ductile H8 mesh contains unused vertices").into());
        }
        Ok(Self {
            nodes_m,
            hexahedra,
            elements: Arc::from(elements),
            materials,
            minimum_jacobians: Arc::from(minimum_jacobians),
            formulation,
            budget,
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            physical_deformation: vec![IDENTITY; points],
            assumed_deformation: vec![IDENTITY; points],
            states,
            first_piola_pa: vec![[0.0; 9]; points],
            stored_energy_j: 0.0,
        })
    }

    pub fn reference_nodes_m(&self) -> &'g [[f64; 3]] {
        self.nodes_m
    }
    pub fn reference_hexahedra(&self) -> &'g [[usize; 8]] {
        self.hexahedra
    }
    pub fn materials(&self) -> &'m [&'m FiniteDuctileJ2<L>] {
        self.materials
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
    pub fn quadrature_states(&self) -> &[FiniteDuctileState<'m, L>] {
        &self.states
    }
    /// Native Piola conjugate to assumed Fbar; nodal forces use the full pullback.
    pub fn quadrature_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn formulation(&self) -> FiniteHexFormulation {
        self.formulation
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements
            .get(element)
            .map(|value| value.reference_volume_m3)
    }
    pub fn quadrature_volume_m3(&self, element: usize, point: usize) -> Option<f64> {
        self.elements
            .get(element)?
            .points
            .get(point)
            .map(|value| value.volume_m3)
    }

    fn point_conditions(
        &self,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<Vec<DuctileJ2Conditions>, FiniteDuctileHexError> {
        checkpoint(cx)?;
        let covered = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(values) => values.len() == self.states.len(),
            DuctileTemperatures::PerHexElement(values) => values.len() == self.elements.len(),
        };
        if !covered || !conditions.duration_s.is_finite() || conditions.duration_s <= 0.0 {
            return Err(
                invalid("complete H8 temperatures and positive actual duration required").into(),
            );
        }
        let mut result = Vec::with_capacity(self.states.len());
        for (index, model) in self.materials.iter().enumerate() {
            if index % 128 == 0 {
                checkpoint(cx)?;
            }
            let temperature_k = match conditions.temperatures_k {
                DuctileTemperatures::Uniform(value) => value,
                DuctileTemperatures::PerPoint(values) => values[index],
                DuctileTemperatures::PerHexElement(values) => values[index / 8],
            };
            let at = DuctileJ2Conditions {
                duration_s: conditions.duration_s,
                temperature_k,
            };
            model
                .validate_conditions(at)
                .map_err(|source| FiniteDuctileHexError::Material {
                    element: index / 8,
                    point: index % 8,
                    source,
                })?;
            result.push(at);
        }
        Ok(result)
    }
    fn operator<'s>(&'s self, conditions: &'s [DuctileJ2Conditions]) -> Operator<'s, 'm, L> {
        Operator {
            elements: &self.elements,
            materials: self.materials,
            accepted: &self.states,
            conditions,
            minimum_jacobians: &self.minimum_jacobians,
            formulation: self.formulation,
            dofs: 3 * self.nodes_m.len(),
            refusal: RefCell::new(None),
        }
    }

    /// Evaluate a physical trial from frozen accepted history, without consuming
    /// time or publishing flow. A rate-dependent law uses the supplied actual
    /// interval, even at unchanged deformation. No potential is claimed.
    pub fn evaluate_trial(
        &self,
        trial_displacement_m: &[[f64; 3]],
        conditions: DuctileStepConditions<'_>,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileHexEvaluation, FiniteDuctileHexError> {
        let conditions = self.point_conditions(conditions, cx)?;
        let operator = self.operator(&conditions);
        let (force, _, tangent_n_m) = operator
            .assemble(
                &self.states,
                trial_displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )
            .map_err(|error| operator.endpoint_refusal(error))?;
        Ok(FiniteDuctileHexEvaluation {
            internal_force_n: force.chunks_exact(3).map(|x| [x[0], x[1], x[2]]).collect(),
            tangent_n_m,
        })
    }

    /// Advance one actual physical interval under absolute endpoint supports and
    /// dead nodal loads. Full geometry and every native history remain original.
    /// The work-defect limit bounds absolute MATERIAL plus PROJECTION defects,
    /// stronger than a gate on their potentially cancelling signed total.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileHexStep<'g, 'm, L>, FiniteDuctileHexError> {
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
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid(
                "finite H8 solver and nonnegative finite work/yield controls required",
            )
            .into());
        }
        let mut fixed = vec![false; dofs];
        for &(dof, value) in prescribed_m {
            if dof >= dofs || !value.is_finite() || fixed[dof] {
                return Err(invalid("invalid or repeated finite H8 support DOF").into());
            }
            fixed[dof] = true;
        }
        if dofs - prescribed_m.len() > self.budget.maximum_free_dofs {
            return Err(invalid("finite H8 free DOF budget exceeded").into());
        }
        let mut external = vec![0.0; dofs];
        for &(dof, value) in nodal_forces_n {
            if dof >= dofs || !value.is_finite() {
                return Err(invalid("invalid finite H8 dead nodal force").into());
            }
            external[dof] += value;
        }
        if external.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite H8 summed nodal force overflow").into());
        }
        let conditions = self.point_conditions(conditions, cx)?;
        for (i, (model, state)) in self.materials.iter().zip(&self.states).enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            model
                .validate_committed(state, controls.yield_tolerance_pa)
                .map_err(|source| FiniteDuctileHexError::Material {
                    element: i / 8,
                    point: i % 8,
                    source,
                })?;
        }
        let operator = self.operator(&conditions);
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
            None,
            cx,
        )
        .map_err(|error| operator.endpoint_refusal(error))?;
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .ok_or_else(|| invalid("native H8 endpoint record is absent"))?;
        if equilibrium.reactions_n.len() != prescribed_m.len() {
            return Err(invalid("native H8 reaction count changed").into());
        }
        let new_u = &solved.displacement_m;
        let old_u = self.displacement_m.as_flattened();
        let count = self.states.len();
        let mut physical = Vec::with_capacity(count);
        let mut assumed = Vec::with_capacity(count);
        let mut stress = Vec::with_capacity(count);
        let mut quadrature_energy = Vec::with_capacity(count);
        // old storage, new storage, storage change, plastic, damage, material
        // work, signed material defect, absolute material defect, point residual.
        let mut sums = [Sum::default(); 9];
        let mut projection = Sum::default();
        let mut absolute_projection = Sum::default();
        let mut local_work_sum = Sum::default();
        let mut internal_force = vec![0.0; dofs];
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let kinematics = operator.kinematics(e, new_u, cx)?;
            let mut cell_material_work = Sum::default();
            for q in 0..8 {
                let index = 8 * e + q;
                let model = self.materials[index];
                let refused = |source| FiniteDuctileHexError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let step = model
                    .trial(
                        &kinematics.assumed_deformation[q],
                        &self.states[index],
                        conditions[index],
                    )
                    .map_err(refused)?;
                if step.state != solved.point_states[index] {
                    return Err(invalid("finite H8 accepted history differs from its original frozen-history return").into());
                }
                model
                    .validate_committed(&step.state, controls.yield_tolerance_pa)
                    .map_err(refused)?;
                let a = step.energy;
                let volume = element.points[q].volume_m3;
                for (sum, value) in sums.iter_mut().zip([
                    a.initial_stored_j_m3,
                    a.final_stored_j_m3,
                    a.final_stored_j_m3 - a.initial_stored_j_m3,
                    a.plastic_dissipation_j_m3,
                    a.damage_release_j_m3,
                    a.mechanical_work_j_m3,
                    a.endpoint_work_defect_j_m3,
                    a.endpoint_work_defect_j_m3.abs(),
                    a.residual_j_m3.abs(),
                ]) {
                    sum.add(volume * value);
                }
                cell_material_work.add(volume * a.mechanical_work_j_m3);
                physical.push(kinematics.physical_deformation[q]);
                assumed.push(kinematics.assumed_deformation[q]);
                stress.push(step.first_piola_pa);
                quadrature_energy.push(a);
            }
            // Reuse accepted point stresses without another constitutive probe.
            let response = element.evaluate_response(
                &kinematics,
                self.formulation == FiniteHexFormulation::MeanDilatation,
                false,
                |q, _, _| {
                    Ok(HexPointResponse {
                        first_piola_pa: stress[8 * e + q],
                        tangent_pa: None,
                    })
                },
                cx,
            )?;
            let mut cell_work = Sum::default();
            for (i, force) in response.local_force_n.iter().enumerate() {
                let dof = 3 * element.nodes[i / 3] + i % 3;
                internal_force[dof] += force;
                cell_work.add(force * (new_u[dof] - old_u[dof]));
            }
            let cell_work = cell_work.total();
            let defect = cell_work - cell_material_work.total();
            local_work_sum.add(cell_work);
            projection.add(defect);
            absolute_projection.add(defect.abs());
        }
        let totals = sums.map(Sum::total);
        if totals.iter().any(|v| !v.is_finite())
            || [totals[0], totals[1], totals[3], totals[4]]
                .iter()
                .any(|v| *v < 0.0)
            || internal_force.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("finite H8 storage/loss/work is not representable").into());
        }
        let mut internal = Sum::default();
        let mut dead = Sum::default();
        let mut support = Sum::default();
        for dof in 0..dofs {
            let delta = new_u[dof] - old_u[dof];
            internal.add(internal_force[dof] * delta);
            dead.add(external[dof] * delta);
        }
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&equilibrium.reactions_n) {
            support.add(reaction * (new_u[dof] - old_u[dof]));
        }
        let mut energy = FiniteDuctileHexEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: totals[1],
            stored_change_j: totals[2],
            plastic_dissipation_j: totals[3],
            damage_release_j: totals[4],
            material_endpoint_work_j: totals[5],
            internal_endpoint_work_j: internal.total(),
            material_endpoint_work_defect_j: totals[6],
            projection_work_defect_j: projection.total(),
            absolute_material_endpoint_work_defect_j: totals[7],
            absolute_projection_work_defect_j: absolute_projection.total(),
            dead_load_endpoint_work_j: dead.total(),
            support_endpoint_work_j: support.total(),
            ..FiniteDuctileHexEnergy::default()
        };
        energy.endpoint_work_defect_j =
            energy.material_endpoint_work_defect_j + energy.projection_work_defect_j;
        energy.absolute_endpoint_work_defect_j = energy.absolute_material_endpoint_work_defect_j
            + energy.absolute_projection_work_defect_j;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        energy.balance_error_j = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j
            - energy.stored_change_j
            - energy.plastic_dissipation_j
            - energy.damage_release_j
            - energy.endpoint_work_defect_j;
        for (quantity, value) in [
            ("retained stored energy", totals[0] - self.stored_energy_j),
            (
                "stored energy difference",
                (totals[1] - totals[0]) - totals[2],
            ),
            ("sum of absolute point energy residuals", totals[8]),
            (
                "cell/global internal work",
                local_work_sum.total() - energy.internal_endpoint_work_j,
            ),
            (
                "material/projection work",
                energy.internal_endpoint_work_j
                    - energy.material_endpoint_work_j
                    - energy.projection_work_defect_j,
            ),
            ("support/dead-load work", energy.equilibrium_work_residual_j),
            ("complete physical work balance", energy.balance_error_j),
        ] {
            gate(quantity, value, controls.equilibrium_work_tolerance_j)?;
        }
        gate(
            "sum of absolute material and projection defects",
            energy.absolute_endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteDuctileHexStep {
            solid: Self {
                nodes_m: self.nodes_m,
                hexahedra: self.hexahedra,
                elements: Arc::clone(&self.elements),
                materials: self.materials,
                minimum_jacobians: Arc::clone(&self.minimum_jacobians),
                formulation: self.formulation,
                budget: self.budget,
                displacement_m: new_u.chunks_exact(3).map(|x| [x[0], x[1], x[2]]).collect(),
                physical_deformation: physical,
                assumed_deformation: assumed,
                states: solved.point_states,
                first_piola_pa: stress,
                stored_energy_j: energy.stored_after_j,
            },
            equilibrium,
            quadrature_energy,
            energy,
        })
    }
}
