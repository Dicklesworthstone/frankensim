//! Progressive finite ductility on the original P1/Q1 sheet continuum.
//!
//! The original three-dimensional material owns plastic flow, stress-state/rate
//! initiation, irreversible attenuation and physical loss. At every reference
//! quadrature point its transverse Piola tractions are relaxed before the full
//! nonsymmetric tangent is condensed. The existing sheet geometry/assembly and
//! Newton/FGMRES increment owner are reused; no scalar potential or symmetric
//! stiffness is invented for evolving damage.
//!
//! A fixed directional reference band width is area / projected transverse
//! span. Original thickness supplies integration volume, not the band width.
//! This projected-footprint regularization does not insert a crack or qualify
//! mesh independence for arbitrary localization. Finite in-plane forming and
//! softening are supplied; bending, wrinkling, separation/contact and thermal
//! coefficient evolution remain outside this owner.

use core::cell::RefCell;
use core::fmt::Debug;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::ductile::evolution::{DifferentiableDuctileInitiationLaw, DuctileJ2Conditions};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::ductile::plane_stress::FiniteDuctilePlaneStressStep;
use fs_material::plastic::finite::ductile::{
    FiniteDuctileEnergy, FiniteDuctileError, FiniteDuctileJ2, FiniteDuctileState,
};
use fs_sparse::Csr;

use super::{
    Element, SheetAssemblyBudget, assemble_sheet, checkpoint, endpoint_forces, invalid,
    prepare_geometry,
};
use crate::ductile::{DuctileStepConditions, DuctileTemperatures};
use crate::hyper3::plastic::FinitePlasticStepControls;
use crate::inertia::{BackwardEulerInertia, InertialEnergy};
use crate::mesh2::{Mesh2, checked_reference_element};
use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// Reference-mass dynamics on the same progressive finite sheet solve.
pub mod dynamics;

type Motion = (Vec<[f64; 2]>, InertialEnergy);

/// Native geometry/solver, original point constitutive, or physical work refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteDuctileSheetError {
    Solid(SmallStrainTetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteDuctileError,
    },
    Energy {
        quantity: &'static str,
        error_j: f64,
        maximum_j: f64,
    },
}
impl From<SmallStrainTetError> for FiniteDuctileSheetError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error)
    }
}
impl core::fmt::Display for FiniteDuctileSheetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile sheet refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileSheetError {}

/// Immutable geometry binding for one in-plane reference localization direction.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetCrackBandGeometry {
    nodes: [usize; 4],
    arity: usize,
    reference: [[f64; 2]; 4],
    normal: [f64; 2],
    area_m2: f64,
    projected_length_m: f64,
    width_m: f64,
}

impl SheetCrackBandGeometry {
    /// Mean directional width of the actual convex P1/Q1 polygon in meters.
    /// The declared normal remains fixed in the reference material frame.
    pub fn from_element(
        mesh: &Mesh2,
        element: usize,
        normal: [f64; 2],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        let conn = mesh
            .elems
            .get(element)
            .ok_or_else(|| invalid("sheet band element is absent"))?;
        let reference = checked_reference_element(&mesh.nodes, conn)
            .ok_or_else(|| invalid("sheet band requires an oriented finite P1/Q1 polygon"))?;
        let norm = normal[0].hypot(normal[1]);
        if !norm.is_finite() || norm <= 0.0 || normal.iter().any(|x| !x.is_finite()) {
            return Err(invalid("sheet band needs a finite nonzero in-plane normal").into());
        }
        let normal = normal.map(|x| x / norm);
        let relative = reference.map(|x| [x[0] - reference[0][0], x[1] - reference[0][1]]);
        let arity = conn.len();
        let mut twice_area = 0.0;
        for a in 1..arity - 1 {
            twice_area += relative[a][0] * relative[a + 1][1] - relative[a][1] * relative[a + 1][0];
        }
        let area_m2 = 0.5 * twice_area;
        let mut low = 0.0_f64;
        let mut high = 0.0_f64;
        for x in &relative[..arity] {
            let projection = -normal[1] * x[0] + normal[0] * x[1];
            if !projection.is_finite() {
                return Err(invalid("sheet band projection is not representable").into());
            }
            low = low.min(projection);
            high = high.max(projection);
        }
        let projected_length_m = high - low;
        let width_m = area_m2 / projected_length_m;
        if [area_m2, projected_length_m, width_m]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(
                invalid("sheet band area, projected span and width must be positive").into(),
            );
        }
        let mut nodes = [0; 4];
        nodes[..arity].copy_from_slice(conn);
        Ok(Self {
            nodes,
            arity,
            reference,
            normal,
            area_m2,
            projected_length_m,
            width_m,
        })
    }
    pub fn band_width_m(&self) -> f64 {
        self.width_m
    }
    pub fn reference_normal(&self) -> [f64; 2] {
        self.normal
    }
    pub fn reference_area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn projected_length_m(&self) -> f64 {
        self.projected_length_m
    }

    fn matches(&self, element: &Element) -> bool {
        self.arity == element.arity
            && self.nodes == element.nodes
            && self.reference == element.reference
    }
}

/// Whole-sheet joules. Numerical endpoint defects are never material heat.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteDuctileSheetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub plastic_dissipation_j: f64,
    pub damage_release_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    pub internal_endpoint_work_j: f64,
    /// Residual transverse P:dF, independently bounded in absolute pointwise sum.
    pub transverse_endpoint_work_j: f64,
    pub endpoint_work_defect_j: f64,
    pub absolute_endpoint_work_defect_j: f64,
    pub equilibrium_work_residual_j: f64,
    pub balance_error_j: f64,
}

/// Force and optional complete tangent evaluated from a frozen accepted history.
#[derive(Debug)]
pub struct FiniteDuctileSheetEvaluation {
    pub internal_force_n: Vec<[f64; 2]>,
    pub tangent_n_m: Option<Csr>,
}

/// Accepted specimen. Cloning shares geometry and retains the same material owners.
#[derive(Clone, Debug)]
pub struct FiniteDuctileSheet<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    mesh: &'g Mesh2,
    materials: &'m [&'m FiniteDuctileJ2<L>],
    thickness_m: &'g [f64],
    bands: &'g [SheetCrackBandGeometry],
    elements: Arc<Vec<Element>>,
    local: PlaneStressSettings,
    budget: SheetAssemblyBudget,
    displacement_m: Vec<[f64; 2]>,
    states: Vec<FiniteDuctileState<'m, L>>,
    stored_energy_j: f64,
    physical_time_s: f64,
}

/// One complete accepted physical increment, published after every work gate.
#[derive(Debug)]
pub struct FiniteDuctileSheetStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub sheet: FiniteDuctileSheet<'g, 'm, L>,
    pub equilibrium: IncrementRecord,
    pub energy: FiniteDuctileSheetEnergy,
    /// Original per-reference-volume accounts in element/quadrature order.
    pub point_energy: Vec<FiniteDuctileEnergy>,
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
fn energy_gate(
    quantity: &'static str,
    error_j: f64,
    maximum_j: f64,
) -> Result<(), FiniteDuctileSheetError> {
    if !error_j.is_finite() || error_j.abs() > maximum_j {
        Err(FiniteDuctileSheetError::Energy {
            quantity,
            error_j,
            maximum_j,
        })
    } else {
        Ok(())
    }
}

struct Operator<'s, 'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    sheet: &'s FiniteDuctileSheet<'g, 'm, L>,
    conditions: &'s [DuctileJ2Conditions],
    inertia: Option<&'s BackwardEulerInertia<'s>>,
    refusal: RefCell<Option<(usize, usize, FiniteDuctileError)>>,
}

impl<'m, L> Operator<'_, '_, 'm, L>
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
            "finite ductile sheet element {element}, point {point}: {source}"
        ));
        *self.refusal.borrow_mut() = Some((element, point, source));
        error
    }
    fn endpoint_refusal(&self, error: SmallStrainTetError) -> FiniteDuctileSheetError {
        if matches!(error, SmallStrainTetError::Cancelled) {
            return error.into();
        }
        match self.refusal.borrow_mut().take() {
            Some((element, point, source)) => FiniteDuctileSheetError::Material {
                element,
                point,
                source,
            },
            None => error.into(),
        }
    }
    fn response(
        &self,
        e: usize,
        q: usize,
        u: &[f64],
        states: &[FiniteDuctileState<'m, L>],
    ) -> Result<FiniteDuctilePlaneStressStep<'m, L>, SmallStrainTetError> {
        let element = &self.sheet.elements[e];
        let index = element.offset + q;
        self.sheet.materials[e]
            .plane_stress_trial(
                &element.in_plane(&element.points[q], u),
                &states[index],
                [2, 5, 8].map(|i| states[index].deformation()[i]),
                self.conditions[index],
                self.sheet.local,
            )
            .map_err(|source| self.material_refusal(e, q, source))
    }
    fn assemble(
        &self,
        states: &[FiniteDuctileState<'m, L>],
        u: &[f64],
        prescribed: Option<&[bool]>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Option<Csr>, f64), SmallStrainTetError> {
        *self.refusal.borrow_mut() = None;
        if states.len() != self.sheet.states.len() {
            return Err(invalid(
                "complete original ductile sheet histories required",
            ));
        }
        assemble_sheet(
            self.sheet.mesh,
            &self.sheet.elements,
            u,
            prescribed,
            tangent,
            self.inertia,
            cx,
            |e, q| {
                let response = self.response(e, q, u, states)?;
                Ok((response.piola_pa, response.tangent_pa))
            },
        )
    }
}

impl<'m, L> IncrementOperator for Operator<'_, '_, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    type State = FiniteDuctileState<'m, L>;
    fn initial_states(&self) -> Vec<Self::State> {
        self.sheet.states.clone()
    }
    fn components_per_node(&self) -> usize {
        2
    }
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let (force, _, cap) = self.assemble(states, u, prescribed, false, cx)?;
        Ok((force, cap))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        Ok(self
            .assemble(states, u, None, true, cx)?
            .1
            .expect("requested tangent"))
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
        for (e, element) in self.sheet.elements.iter().enumerate() {
            let mut evolving = false;
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let proposal = self.response(e, q, u, states)?.step.state;
                evolving |= proposal.native_state() != states[index].native_state()
                    || proposal.damage() != states[index].damage();
                states[index] = proposal;
            }
            changed += usize::from(evolving);
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}

impl<'g, 'm, L> FiniteDuctileSheet<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Bind virgin material histories to original quadrature, thickness and band
    /// geometry. Construction does not consume a fictitious physical interval.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        mesh: &'g Mesh2,
        materials: &'m [&'m FiniteDuctileJ2<L>],
        thickness_m: &'g [f64],
        bands: &'g [SheetCrackBandGeometry],
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        if materials.len() != mesh.elems.len()
            || bands.len() != mesh.elems.len()
            || !local.relative_tolerance.is_finite()
            || local.relative_tolerance <= 0.0
            || local.relative_tolerance >= 1.0
            || local.max_iterations == 0
            || local.max_backtracks == 0
        {
            return Err(invalid(
                "one original material/band per sheet element and finite local budgets required",
            )
            .into());
        }
        let elements = prepare_geometry(mesh, thickness_m, budget, cx)?;
        let count = elements.last().map_or(0, |e| e.offset + e.points.len());
        let mut states = Vec::with_capacity(count);
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            if !bands[e].matches(element)
                || materials[e].settings().characteristic_length_m != bands[e].band_width_m()
            {
                return Err(invalid(format!(
                    "sheet element {e} damage length differs from its original directional band"
                ))
                .into());
            }
            states.extend((0..element.points.len()).map(|_| materials[e].initial_state()));
        }
        Ok(Self {
            mesh,
            materials,
            thickness_m,
            bands,
            elements,
            local,
            budget,
            displacement_m: vec![[0.0; 2]; mesh.nodes.len()],
            states,
            stored_energy_j: 0.0,
            physical_time_s: 0.0,
        })
    }

    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    pub fn quadrature_states(&self) -> &[FiniteDuctileState<'m, L>] {
        &self.states
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn physical_time_s(&self) -> f64 {
        self.physical_time_s
    }
    pub fn materials(&self) -> &'m [&'m FiniteDuctileJ2<L>] {
        self.materials
    }
    pub fn bands(&self) -> &'g [SheetCrackBandGeometry] {
        self.bands
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements
            .get(element)
            .map(|e| e.points.iter().map(|p| p.volume_m3).sum())
    }
    /// Original point volume in native element/quadrature order. This is not
    /// the current volume or an equal share inferred from an element mean.
    pub fn quadrature_volume_m3(&self, element: usize, point: usize) -> Option<f64> {
        self.elements.get(element)?.points.get(point).map(|p| p.volume_m3)
    }
    pub fn quadrature_reference_volumes_m3(&self) -> impl Iterator<Item = f64> + '_ {
        self.elements
            .iter()
            .flat_map(|e| e.points.iter().map(|p| p.volume_m3))
    }
    pub fn quadrature_thickness_m(&self) -> Vec<f64> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(e, element)| {
                (0..element.points.len()).map(move |q| {
                    self.thickness_m[e] * self.states[element.offset + q].deformation()[8]
                })
            })
            .collect()
    }

    fn point_conditions(
        &self,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<Vec<DuctileJ2Conditions>, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        let covered = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(values) => values.len() == self.states.len(),
            DuctileTemperatures::PerHexElement(_) => false,
        };
        if !covered || !conditions.duration_s.is_finite() || conditions.duration_s <= 0.0 {
            return Err(invalid(
                "ductile sheet needs actual point temperatures and one positive physical duration",
            )
            .into());
        }
        let mut points = Vec::with_capacity(self.states.len());
        for (e, element) in self.elements.iter().enumerate() {
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let temperature_k = match conditions.temperatures_k {
                    DuctileTemperatures::Uniform(t) => t,
                    DuctileTemperatures::PerPoint(values) => values[index],
                    DuctileTemperatures::PerHexElement(_) => unreachable!("coverage admitted"),
                };
                let at = DuctileJ2Conditions {
                    duration_s: conditions.duration_s,
                    temperature_k,
                };
                self.materials[e]
                    .validate_conditions(at)
                    .map_err(|source| FiniteDuctileSheetError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
                points.push(at);
            }
        }
        Ok(points)
    }

    /// Frozen-history force/derivative observation. The accepted clock is unchanged.
    pub fn evaluate_trial(
        &self,
        displacement_m: &[[f64; 2]],
        conditions: DuctileStepConditions<'_>,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileSheetEvaluation, FiniteDuctileSheetError> {
        let conditions = self.point_conditions(conditions, cx)?;
        let operator = Operator {
            sheet: self,
            conditions: &conditions,
            inertia: None,
            refusal: RefCell::new(None),
        };
        let (force, tangent, _) = operator
            .assemble(
                &self.states,
                displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )
            .map_err(|error| operator.endpoint_refusal(error))?;
        Ok(FiniteDuctileSheetEvaluation {
            internal_force_n: force.chunks_exact(2).map(|v| [v[0], v[1]]).collect(),
            tangent_n_m: tangent,
        })
    }

    /// Advance one physical interval with absolute supports and dead nodal loads.
    /// Displacement control can traverse a falling load branch; force control
    /// may refuse an unstable branch and does not claim arclength continuation.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileSheetStep<'g, 'm, L>, FiniteDuctileSheetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            conditions,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    // One nonpotential equilibrium/work owner serves both static and moving
    // specimens. The original inertia owns momentum, velocity and kinetic loss.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_loaded(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(FiniteDuctileSheetStep<'g, 'm, L>, Option<Motion>), FiniteDuctileSheetError> {
        checkpoint(cx)?;
        let external = endpoint_forces(
            self.mesh,
            self.budget,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
        )?;
        let point_conditions = self.point_conditions(conditions, cx)?;
        let time = self.physical_time_s + conditions.duration_s;
        if !time.is_finite() || time <= self.physical_time_s {
            return Err(invalid("ductile sheet physical time cannot advance representably").into());
        }
        let mut old_storage = Sum::default();
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let state = &self.states[element.offset + q];
                self.materials[e]
                    .validate_committed(state, controls.yield_tolerance_pa)
                    .map_err(|source| FiniteDuctileSheetError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
                old_storage.add(point.volume_m3 * state.stored_energy_j_m3());
            }
        }
        energy_gate(
            "retained storage",
            old_storage.value() - self.stored_energy_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        let operator = Operator {
            sheet: self,
            conditions: &point_conditions,
            inertia,
            refusal: RefCell::new(None),
        };
        let solved = increment::solve(
            &operator,
            IncrementControls {
                dofs: 2 * self.mesh.nodes.len(),
                prescribed_m,
                nodal_forces_n,
                load_path: &[1.0],
                settings,
            },
            Some((self.displacement_m.as_flattened(), &self.states)),
            inertia,
            cx,
        )
        .map_err(|error| operator.endpoint_refusal(error))?;
        let motion = inertia
            .map(|inertia| inertia.motion_flat(&solved.displacement_m, cx))
            .transpose()
            .map_err(SmallStrainTetError::from)?
            .map(|(velocity, energy)| {
                (
                    velocity.chunks_exact(2).map(|v| [v[0], v[1]]).collect(),
                    energy,
                )
            });
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .ok_or_else(|| invalid("sheet endpoint record absent"))?;
        if equilibrium.reactions_n.len() != prescribed_m.len() {
            return Err(invalid("sheet endpoint reaction coverage changed").into());
        }
        let mut totals = [Sum::default(); 9];
        let mut point_energy = Vec::with_capacity(self.states.len());
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let response = operator
                    .response(e, q, &solved.displacement_m, &self.states)
                    .map_err(|error| operator.endpoint_refusal(error))?;
                let step = response.step;
                if step.state != solved.point_states[index] {
                    return Err(invalid(
                        "sheet audit differs from original frozen-history proposal",
                    )
                    .into());
                }
                self.materials[e]
                    .validate_committed(&step.state, controls.yield_tolerance_pa)
                    .map_err(|source| FiniteDuctileSheetError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
                let transverse: f64 = [2, 5, 8]
                    .into_iter()
                    .map(|i| {
                        step.first_piola_pa[i]
                            * (response.deformation_gradient[i]
                                - self.states[index].deformation()[i])
                    })
                    .sum();
                for (sum, value) in totals.iter_mut().zip([
                    step.energy.final_stored_j_m3,
                    step.energy.plastic_dissipation_j_m3,
                    step.energy.damage_release_j_m3,
                    step.energy.mechanical_work_j_m3,
                    step.energy.endpoint_work_defect_j_m3,
                    step.energy.endpoint_work_defect_j_m3.abs(),
                    step.energy.residual_j_m3.abs(),
                    transverse,
                    transverse.abs(),
                ]) {
                    sum.add(point.volume_m3 * value);
                }
                point_energy.push(step.energy);
            }
        }
        let values = totals.map(Sum::value);
        if values.iter().any(|v| !v.is_finite()) || values[..3].iter().any(|v| *v < 0.0) {
            return Err(
                invalid("sheet storage and physical loss are not finite and nonnegative").into(),
            );
        }
        let mut dead_work = Sum::default();
        let old_u = self.displacement_m.as_flattened();
        for (dof, force) in external.iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            dead_work.add(force * (solved.displacement_m[dof] - old_u[dof]));
        }
        let mut support_work = Sum::default();
        for (index, (&(dof, _), reaction)) in prescribed_m
            .iter()
            .zip(&equilibrium.reactions_n)
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            support_work.add(reaction * (solved.displacement_m[dof] - old_u[dof]));
        }
        let external_work = dead_work.value() + support_work.value();
        let inertial_work = motion
            .as_ref()
            .map_or(0.0, |(_, energy)| energy.inertial_work_j);
        let energy = FiniteDuctileSheetEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: values[0],
            plastic_dissipation_j: values[1],
            damage_release_j: values[2],
            dead_load_endpoint_work_j: dead_work.value(),
            support_endpoint_work_j: support_work.value(),
            internal_endpoint_work_j: values[3],
            transverse_endpoint_work_j: values[7],
            endpoint_work_defect_j: values[4],
            absolute_endpoint_work_defect_j: values[5],
            equilibrium_work_residual_j: external_work - values[3] - inertial_work,
            balance_error_j: external_work
                - (values[0] - self.stored_energy_j)
                - values[1]
                - values[2]
                - values[4]
                - inertial_work,
        };
        for (name, value) in [
            ("point constitutive work identity", values[6]),
            ("sum of absolute transverse residual work", values[8]),
            (
                "support/dead-load equilibrium work",
                energy.equilibrium_work_residual_j,
            ),
            ("complete physical work balance", energy.balance_error_j),
        ] {
            energy_gate(name, value, controls.equilibrium_work_tolerance_j)?;
        }
        energy_gate(
            "sum of absolute point work defects",
            values[5],
            controls.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        Ok((
            FiniteDuctileSheetStep {
                sheet: Self {
                    mesh: self.mesh,
                    materials: self.materials,
                    thickness_m: self.thickness_m,
                    bands: self.bands,
                    elements: Arc::clone(&self.elements),
                    local: self.local,
                    budget: self.budget,
                    displacement_m: solved
                        .displacement_m
                        .chunks_exact(2)
                        .map(|v| [v[0], v[1]])
                        .collect(),
                    states: solved.point_states,
                    stored_energy_j: energy.stored_after_j,
                    physical_time_s: time,
                },
                equilibrium,
                energy,
                point_energy,
            },
            motion,
        ))
    }
}
