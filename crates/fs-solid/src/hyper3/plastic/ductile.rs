//! Finite-deformation ductile softening on the original P1 reference mesh.
//!
//! Multiplicative plastic history, finite rotations and the complete damaged
//! first-Piola tangent reach the existing Newton/FGMRES increment owner. An
//! evolving damage response is generally nonpotential: it is never presented
//! as a `FiniteStrainPotential` or symmetrized to enter an energy minimizer.
//! Reference geometry, directional band length and immutable material owners
//! survive every step. Failed proposals publish no displacement or history.
//! Explicit native overstress uses the same physical duration as initiation
//! and attached inertia. Accepted stress inspection never advances that clock.
//!
//! Physical plastic dissipation and damage release remain separate from the
//! signed endpoint work-integration defect. The work budget bounds the SUM of
//! absolute element defects, so opposing errors cannot conceal underresolution.
//! Drive unstable softening by prescribed displacement or use smaller physical
//! intervals; explicit finite-support fracture transfer and matched-trace contact
//! are separate child owners. This owner does not supply arclength, arbitrary
//! remeshing or thermal feedback. P1 volumetric locking and
//! the parent's displacement-only applicability limits remain in force.

use core::cell::RefCell;
use core::fmt::Debug;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::ductile::evolution::{DifferentiableDuctileInitiationLaw, DuctileJ2Conditions};
use fs_material::plastic::finite::ductile::{
    FiniteDuctileEnergy, FiniteDuctileError, FiniteDuctileJ2, FiniteDuctileState,
};
use fs_sparse::{Coo, Csr};

use super::{
    Element, FinitePlasticStepControls, FinitePlasticTetError, FinitePlasticTetState, FrozenPoint,
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, IDENTITY,
    TetAssemblyBudget,
};
use crate::brittle::{TetCrackBandGeometry, TetCrackBandGeometryError};
use crate::ductile::{DuctileStepConditions, DuctileTemperatures};
use crate::inertia::{BackwardEulerInertia, InertialEnergy};
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// Frictionless translating tools on the same progressive material and momentum.
pub mod contact;
/// Consistent reference inertia on the same nonpotential finite ductile solve.
pub mod dynamics;
/// Native complete-failure separation with original histories and consistent mass.
pub mod fracture;
/// Bounded physical-interval continuation through native finite softening.
pub mod path;

type Motion = (Vec<[f64; 3]>, InertialEnergy);

/// Original geometry, material, solver or physical-work refusal.
#[derive(Debug)]
pub enum FiniteDuctileTetError {
    Reference(FinitePlasticTetError),
    Band(TetCrackBandGeometryError),
    Solve(SmallStrainTetError),
    Material {
        element: usize,
        source: FiniteDuctileError,
    },
    Energy {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl core::fmt::Display for FiniteDuctileTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Reference(error) => core::fmt::Display::fmt(error, f),
            Self::Band(error) => core::fmt::Display::fmt(error, f),
            Self::Solve(error) => core::fmt::Display::fmt(error, f),
            Self::Material { element, source } => {
                write!(f, "finite ductile tetrahedron {element}: {source}")
            }
            Self::Energy {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "finite ductile {quantity}: {value_j:e} J exceeds {maximum_abs_j:e} J"
            ),
        }
    }
}
impl std::error::Error for FiniteDuctileTetError {}
impl From<FinitePlasticTetError> for FiniteDuctileTetError {
    fn from(value: FinitePlasticTetError) -> Self {
        Self::Reference(value)
    }
}
impl From<HyperTetError> for FiniteDuctileTetError {
    fn from(value: HyperTetError) -> Self {
        Self::Reference(value.into())
    }
}
impl From<SmallStrainTetError> for FiniteDuctileTetError {
    fn from(value: SmallStrainTetError) -> Self {
        Self::Solve(value)
    }
}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}
fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteDuctileTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FiniteDuctileTetError::Energy {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

/// Whole-reference-mesh joules for one accepted physical interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteDuctileTetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub plastic_dissipation_j: f64,
    pub damage_release_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    pub internal_endpoint_work_j: f64,
    /// Signed material work quadrature defect, never heat or fracture energy.
    pub endpoint_work_defect_j: f64,
    /// Sum of absolute element defects; this quantity consumes the work budget.
    pub absolute_endpoint_work_defect_j: f64,
    /// External/support work minus material and any attached inertial work.
    pub equilibrium_work_residual_j: f64,
    /// Complete material balance after subtracting any attached inertial work.
    pub balance_error_j: f64,
}

/// Immutable accepted specimen, including the original per-element band records.
#[derive(Debug)]
pub struct FiniteDuctileTetSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    materials: &'m [&'m FiniteDuctileJ2<L>],
    bands: &'g [TetCrackBandGeometry],
    state: FinitePlasticTetState<'g, FiniteDuctileState<'m, L>>,
}

/// Complete accepted configuration and native nonlinear-equilibrium record.
#[derive(Debug)]
pub struct FiniteDuctileTetStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub solid: FiniteDuctileTetSolid<'g, 'm, L>,
    /// Reactions follow the input prescribed-DOF order.
    pub equilibrium: IncrementRecord,
    pub energy: FiniteDuctileTetEnergy,
    /// Original final material proposals, per reference volume in element order.
    /// Coupled consumers need not difference cumulative loss or repeat a return.
    pub element_energy: Vec<FiniteDuctileEnergy>,
}

struct Operator<'s, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    elements: &'s [Element],
    materials: &'m [&'m FiniteDuctileJ2<L>],
    accepted: &'s [FiniteDuctileState<'m, L>],
    conditions: &'s [DuctileJ2Conditions],
    inertia: Option<&'s BackwardEulerInertia<'s>>,
    dofs: usize,
    /// Only the CURRENT failed evaluation; successful trials clear this record.
    /// Preserve typed damage/source refusals across the shared solver boundary.
    refusal: RefCell<Option<(usize, FiniteDuctileError)>>,
}

#[derive(Clone, Copy)]
enum MaterialEvaluation {
    Stress,
    Tangent,
    /// Retained physical stress, used only to inspect an initial equilibrium.
    Accepted,
}

impl<'m, L> Operator<'_, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn material_refusal(&self, element: usize, source: FiniteDuctileError) -> SmallStrainTetError {
        let error = invalid(format!("finite ductile element {element}: {source}"));
        *self.refusal.borrow_mut() = Some((element, source));
        error
    }

    fn endpoint_refusal(&self, error: SmallStrainTetError) -> FiniteDuctileTetError {
        if matches!(error, SmallStrainTetError::Cancelled) {
            return error.into();
        }
        match self.refusal.borrow_mut().take() {
            Some((element, source)) => FiniteDuctileTetError::Material { element, source },
            None => error.into(),
        }
    }

    fn assemble(
        &self,
        states: &[FiniteDuctileState<'m, L>],
        u: &[f64],
        prescribed: Option<&[bool]>,
        evaluation: MaterialEvaluation,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64, Option<Csr>), SmallStrainTetError> {
        *self.refusal.borrow_mut() = None;
        checkpoint(cx)?;
        if states.len() != self.elements.len()
            || u.len() != self.dofs
            || u.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite ductile displacement/history fields required",
            ));
        }
        let mut force = vec![0.0; self.dofs];
        let mut scale = 0.0;
        let mut matrix = matches!(evaluation, MaterialEvaluation::Tangent)
            .then(|| Coo::new(self.dofs, self.dofs));
        for (e, element) in self.elements.iter().enumerate() {
            if e % 128 == 0 {
                checkpoint(cx)?;
            }
            let f = HyperTetProblem::<FrozenPoint<'_>>::deformation_gradient(element, u);
            let refused = |error| self.material_refusal(e, error);
            let (p, derivative) = match evaluation {
                MaterialEvaluation::Tangent => {
                    let step = self.materials[e]
                        .trial(&f, &states[e], self.conditions[e])
                        .map_err(refused)?;
                    (step.first_piola_pa, Some(step.tangent_pa))
                }
                MaterialEvaluation::Stress => (
                    self.materials[e]
                        .stress(&f, &states[e], self.conditions[e])
                        .map_err(refused)?,
                    None,
                ),
                MaterialEvaluation::Accepted => {
                    if !self.materials[e].owns_state(&states[e])
                        || !states[e].is_valid()
                        || &f != states[e].deformation()
                    {
                        return Err(invalid(
                            "retained stress requires its original accepted material and deformation",
                        ));
                    }
                    (*states[e].first_piola_pa(), None)
                }
            };
            for a in 0..4 {
                for i in 0..3 {
                    let row = 3 * element.nodes[a] + i;
                    let value = element.volume
                        * (0..3).fold(0.0, |sum, j| {
                            p[3 * i + j].mul_add(element.gradients[a][j], sum)
                        });
                    force[row] += value;
                    if prescribed.is_none_or(|fixed| !fixed[row]) {
                        scale += value.abs();
                    }
                    if let (Some(coo), Some(derivative)) = (&mut matrix, &derivative) {
                        for b in 0..4 {
                            for k in 0..3 {
                                let mut value = 0.0;
                                for j in 0..3 {
                                    for l in 0..3 {
                                        value = (element.gradients[a][j]
                                            * derivative[3 * i + j][3 * k + l])
                                            .mul_add(element.gradients[b][l], value);
                                    }
                                }
                                value *= element.volume;
                                if !value.is_finite() {
                                    return Err(invalid(
                                        "finite ductile tangent cannot be represented",
                                    ));
                                }
                                if value != 0.0 {
                                    coo.push(row, 3 * element.nodes[b] + k, value);
                                }
                            }
                        }
                    }
                }
            }
        }
        if !scale.is_finite() || force.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite ductile assembled force cannot be represented",
            ));
        }
        // The shared increment owner adds inertia to residual forces. Its
        // operator must add the SAME original derivative exactly once.
        if let (Some(matrix), Some(inertia)) = (&mut matrix, self.inertia) {
            inertia.add_tangent(matrix, cx)?;
        }
        checkpoint(cx)?;
        Ok((force, scale, matrix.map(|coo| coo.assemble())))
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
        let (force, scale, _) =
            self.assemble(states, u, prescribed, MaterialEvaluation::Stress, cx)?;
        Ok((force, scale))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        Ok(self
            .assemble(states, u, None, MaterialEvaluation::Tangent, cx)?
            .2
            .expect("requested full tangent"))
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
        for (e, (element, state)) in self.elements.iter().zip(states).enumerate() {
            if e % 128 == 0 {
                checkpoint(cx)?;
            }
            let f = HyperTetProblem::<FrozenPoint<'_>>::deformation_gradient(element, u);
            let step = self.materials[e]
                .trial(&f, state, self.conditions[e])
                .map_err(|error| self.material_refusal(e, error))?;
            changed += usize::from(
                step.state.native_state() != state.native_state()
                    || step.state.damage() != state.damage(),
            );
            *state = step.state;
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}

impl<'g, 'm, L> FiniteDuctileTetSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Admit the virgin reference configuration and original directional lengths.
    /// Geometry is prepared once by the existing finite tetrahedral owner.
    pub fn try_new(
        geometry: &'g TetLoadGeometry3<'g>,
        materials: &'m [&'m FiniteDuctileJ2<L>],
        bands: &'g [TetCrackBandGeometry],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileTetError> {
        checkpoint(cx)?;
        let nodes_m = geometry.nodes_m();
        let tetrahedra = geometry.tetrahedra();
        if !budget.minimum_scaled_jacobian.is_finite() || budget.minimum_scaled_jacobian <= 0.0 {
            return Err(
                invalid("positive finite reference-geometry quality bound required").into(),
            );
        }
        super::reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        if materials.len() != tetrahedra.len() || bands.len() != tetrahedra.len() {
            return Err(invalid(
                "one original finite ductile material and band per tetrahedron required",
            )
            .into());
        }
        for (e, (model, band)) in materials.iter().zip(bands).enumerate() {
            if e % 128 == 0 {
                checkpoint(cx)?;
            }
            let original = TetCrackBandGeometry::from_element(geometry, e, band.normal())
                .map_err(FiniteDuctileTetError::Band)?;
            if original.element() != band.element()
                || original.nodes() != band.nodes()
                || original.reference_volume_m3() != band.reference_volume_m3()
                || (original.projected_area_m2() / band.projected_area_m2() - 1.0).abs()
                    > 64.0 * f64::EPSILON
                || (original.band_width_m() / band.band_width_m() - 1.0).abs() > 64.0 * f64::EPSILON
                || model.settings().characteristic_length_m != band.band_width_m()
            {
                return Err(invalid(
                    "finite ductile band and material length must match the original element",
                )
                .into());
            }
        }
        let element_states: Vec<_> = materials
            .iter()
            .map(|model| model.initial_state())
            .collect();
        let frozen: Vec<_> = materials
            .iter()
            .zip(&element_states)
            .map(|(model, state)| FrozenPoint {
                law: model.native(),
                accepted: state.native_state(),
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
            materials,
            bands,
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

    pub fn reference_nodes_m(&self) -> &'g [[f64; 3]] {
        self.state.nodes_m
    }
    pub fn reference_tetrahedra(&self) -> &'g [[usize; 4]] {
        self.state.tetrahedra
    }
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.state.displacement_m
    }
    pub fn element_deformation(&self) -> &[[f64; 9]] {
        &self.state.element_deformation
    }
    pub fn element_states(&self) -> &[FiniteDuctileState<'m, L>] {
        &self.state.element_states
    }
    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.state.element_first_piola_pa
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.state.stored_energy_j
    }
    pub fn bands(&self) -> &'g [TetCrackBandGeometry] {
        self.bands
    }
    pub fn materials(&self) -> &'m [&'m FiniteDuctileJ2<L>] {
        self.materials
    }

    /// Advance one physical endpoint with absolute supports and dead/body loads.
    /// The work-defect allowance bounds sum(V*abs(element defect)), separately
    /// from force-equilibrium work. Reduce the interval after a damage/work refusal.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileTetStep<'g, 'm, L>, FiniteDuctileTetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            conditions,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    /// Common material, equilibrium and work owner for static and moving
    /// specimens. Only the original inertia owner supplies momentum or velocity.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_loaded(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(FiniteDuctileTetStep<'g, 'm, L>, Option<Motion>), FiniteDuctileTetError> {
        checkpoint(cx)?;
        let count = self.materials.len();
        let covered = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(values) => values.len() == count,
            DuctileTemperatures::PerHexElement(_) => false,
        };
        if !covered || !conditions.duration_s.is_finite() || conditions.duration_s <= 0.0 {
            return Err(invalid(
                "complete tetrahedral temperatures and positive physical duration required",
            )
            .into());
        }
        let settings_hyper = HyperTetSettings {
            load_steps: 1,
            max_newton_iterations: settings.max_newton_iterations,
            relative_tolerance: settings.relative_tolerance,
            max_backtracks: settings.max_backtracks,
            linear_tolerance: settings.linear_tolerance,
            max_linear_iterations: settings.max_linear_iterations,
        };
        super::validate_controls(settings_hyper, controls)?;
        let mut point_conditions = Vec::with_capacity(count);
        for (e, (model, state)) in self
            .materials
            .iter()
            .zip(&self.state.element_states)
            .enumerate()
        {
            if e % 128 == 0 {
                checkpoint(cx)?;
            }
            let temperature_k = match conditions.temperatures_k {
                DuctileTemperatures::Uniform(value) => value,
                DuctileTemperatures::PerPoint(values) => values[e],
                DuctileTemperatures::PerHexElement(_) => unreachable!("coverage refused above"),
            };
            let at = DuctileJ2Conditions {
                duration_s: conditions.duration_s,
                temperature_k,
            };
            let refused = |source| FiniteDuctileTetError::Material { element: e, source };
            model
                .validate_committed(state, controls.yield_tolerance_pa)
                .map_err(refused)?;
            model.validate_conditions(at).map_err(refused)?;
            point_conditions.push(at);
        }
        let frozen: Vec<_> = self
            .materials
            .iter()
            .zip(&self.state.element_states)
            .map(|(model, state)| FrozenPoint {
                law: model.native(),
                accepted: state.native_state(),
            })
            .collect();
        // Original admission and body-force integration; damage contributes no
        // fictitious potential to this geometry/load preparation.
        let load_problem = HyperTetProblem {
            nodes_m: self.state.nodes_m,
            tetrahedra: self.state.tetrahedra,
            material: HyperelasticLawField::PerElement(&frozen),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            budget: self.state.budget,
            settings: settings_hyper,
        };
        load_problem.validate_inputs()?;
        let force = load_problem.external(&self.state.elements, 1.0);
        let nodal: Vec<_> = force
            .iter()
            .copied()
            .enumerate()
            .filter(|(_, v)| *v != 0.0)
            .collect();
        let operator = Operator {
            elements: &self.state.elements,
            materials: self.materials,
            accepted: &self.state.element_states,
            conditions: &point_conditions,
            inertia,
            dofs: self.state.nodes_m.len() * 3,
            refusal: RefCell::new(None),
        };
        let solved = increment::solve(
            &operator,
            IncrementControls {
                dofs: operator.dofs,
                prescribed_m,
                nodal_forces_n: &nodal,
                load_path: &[1.0],
                settings,
            },
            Some((
                self.state.displacement_m.as_flattened(),
                &self.state.element_states,
            )),
            inertia,
            cx,
        )
        .map_err(|error| operator.endpoint_refusal(error))?;
        let motion = inertia
            .map(|inertia| inertia.motion(&solved.displacement_m, cx))
            .transpose()
            .map_err(SmallStrainTetError::from)?;
        let record = solved
            .increments
            .into_iter()
            .last()
            .ok_or_else(|| invalid("native finite ductile endpoint record is absent"))?;
        if record.reactions_n.len() != prescribed_m.len() {
            return Err(invalid("native finite ductile support reaction count changed").into());
        }
        let mut deformation = Vec::with_capacity(count);
        let mut piola = Vec::with_capacity(count);
        let mut element_energy = Vec::with_capacity(count);
        let mut totals = [Sum::default(); 7];
        for (e, element) in self.state.elements.iter().enumerate() {
            if e % 128 == 0 {
                checkpoint(cx)?;
            }
            let model = self.materials[e];
            let f = HyperTetProblem::<FrozenPoint<'_>>::deformation_gradient(
                element,
                &solved.displacement_m,
            );
            let refused = |source| FiniteDuctileTetError::Material { element: e, source };
            let step = model
                .trial(&f, &self.state.element_states[e], point_conditions[e])
                .map_err(refused)?;
            if step.state != solved.point_states[e] {
                return Err(invalid("native finite ductile committed history differs from its frozen-history proposal").into());
            }
            model
                .validate_committed(&step.state, controls.yield_tolerance_pa)
                .map_err(refused)?;
            for (sum, value) in totals.iter_mut().zip([
                step.energy.final_stored_j_m3,
                step.energy.plastic_dissipation_j_m3,
                step.energy.damage_release_j_m3,
                step.energy.mechanical_work_j_m3,
                step.energy.endpoint_work_defect_j_m3,
                step.energy.endpoint_work_defect_j_m3.abs(),
                step.energy.residual_j_m3.abs(),
            ]) {
                sum.add(element.volume * value);
            }
            deformation.push(f);
            piola.push(step.first_piola_pa);
            element_energy.push(step.energy);
        }
        let values = totals.map(Sum::value);
        if values.iter().any(|v| !v.is_finite()) || values[..3].iter().any(|v| *v < 0.0) {
            return Err(invalid("finite nonnegative storage and physical losses required").into());
        }
        let reactions = prescribed_m
            .iter()
            .zip(&record.reactions_n)
            .map(|(&(dof, _), &reaction)| (dof, reaction));
        let work = self.state.endpoint_work_fields(
            &solved.displacement_m,
            reactions,
            &force,
            &deformation,
            &piola,
            motion.as_ref().map(|(_, account)| account.inertial_work_j),
            cx,
        )?;
        let energy = FiniteDuctileTetEnergy {
            stored_before_j: self.state.stored_energy_j,
            stored_after_j: values[0],
            plastic_dissipation_j: values[1],
            damage_release_j: values[2],
            dead_load_endpoint_work_j: work.dead_load_j,
            support_endpoint_work_j: work.support_j,
            internal_endpoint_work_j: work.internal_j,
            endpoint_work_defect_j: values[4],
            absolute_endpoint_work_defect_j: values[5],
            equilibrium_work_residual_j: work.residual_j,
            balance_error_j: ((work.dead_load_j + work.support_j)
                - (values[0] - self.state.stored_energy_j))
                - values[1]
                - values[2]
                - values[4]
                - motion
                    .as_ref()
                    .map_or(0.0, |(_, account)| account.inertial_work_j),
        };
        for (name, value) in [
            ("element work consistency", values[6]),
            (
                "constitutive/assembled endpoint work",
                values[3] - work.internal_j,
            ),
            (
                "support/dead-load equilibrium work",
                energy.equilibrium_work_residual_j,
            ),
            ("complete physical work balance", energy.balance_error_j),
        ] {
            gate(name, value, controls.equilibrium_work_tolerance_j)?;
        }
        gate(
            "sum of absolute element integration defects",
            energy.absolute_endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        let displacement_m = solved
            .displacement_m
            .chunks_exact(3)
            .map(|v| [v[0], v[1], v[2]])
            .collect();
        Ok((
            FiniteDuctileTetStep {
                solid: Self {
                    materials: self.materials,
                    bands: self.bands,
                    state: FinitePlasticTetState {
                        nodes_m: self.state.nodes_m,
                        tetrahedra: self.state.tetrahedra,
                        elements: Arc::clone(&self.state.elements),
                        budget: self.state.budget,
                        displacement_m,
                        element_deformation: deformation,
                        element_states: solved.point_states,
                        element_first_piola_pa: piola,
                        stored_energy_j: energy.stored_after_j,
                    },
                },
                equilibrium: record,
                energy,
                element_energy,
            },
            motion,
        ))
    }
}
