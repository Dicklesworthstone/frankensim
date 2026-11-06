//! Finite-strain membrane loading and relaxation with retained material memory.
//!
//! Each original triangle keeps its native three-dimensional Maxwell history.
//! Plane-stress relaxation, the objective Green-strain chain rule, reference
//! geometry and matrix-free Newton solve stay with their existing owners. Every
//! rejected Newton trial sees the same accepted metrics and physical interval.
//! Only a complete, equilibrated and work-admitted proposal publishes new
//! displacement, transverse deformation, stress, storage, loss and time.
//!
//! The equilibrium law and natural reference must exactly match the original
//! source-bound hyperelastic sheet. Branch parameters are explicit caller data;
//! this binder does not claim that the source card characterizes their spectrum.
//! Areal mass, original thickness and the selected finite/exact-volume treatment
//! remain unchanged. No inertia, bending, wrinkling, contact, damage or evolving
//! temperature is inferred by the quasistatic methods.

/// Native inertial intervals with physical relaxation and separate numerical loss.
pub mod dynamics;

/// Prescribed natural-strain events with retained material memory and inventory.
pub mod reference_change;
pub use reference_change::{
    FiniteViscoMembraneMotionReferenceChange, FiniteViscoMembraneReferenceChange,
    FiniteViscoMembraneReferenceEnergy,
};

use std::cell::RefCell;
use std::sync::Arc;

use fs_material::sheet::hyperelastic::HyperelasticMembrane;
use fs_material::sheet::{MembraneLaw, SheetError, SheetLaw, SheetResponse};
use fs_material::visco_hyper::membrane::{
    FrozenViscoMembrane, ViscoMembraneState, ViscoMembraneStepResponse,
};
use fs_material::visco_hyper::{MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic};

pub use crate::hyper3::visco::{FiniteViscoStepControls, ViscoelasticLawField};

use super::{
    MembraneError, MembraneInertialPotential, MembranePatch, MembraneSolution,
    MembraneSolveSettings, finite_vectors, poll,
};

/// Physical and numerical accounts for one accepted physical interval [J].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteViscoMembraneEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Native nonnegative branch relaxation loss, not cumulative.
    pub dissipated_j: f64,
    /// Frozen-metric constitutive loading work, not nodal quadrature.
    pub native_loading_work_j: f64,
    /// Newton's reduced material potential, not physical stored energy.
    pub incremental_potential_j: f64,
    /// Native loading work minus storage change and relaxation loss.
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    /// Current-area pressure at the accepted endpoint dotted with displacement.
    pub pressure_endpoint_work_j: f64,
    /// Exact held-pressure work on the straight accepted displacement segment.
    /// A moving perimeter makes this a path work, not a global potential.
    pub pressure_path_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Zero in quasistatics; supplied by the original inertia potential otherwise.
    pub inertial_endpoint_work_j: f64,
    /// Applied endpoint work minus internal and inertial endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus physical storage change and native loss.
    /// This signed integration defect is never added to physical heat.
    pub endpoint_work_defect_j: f64,
}

/// Original patch, explicit branch laws and complete immutable accepted state.
#[derive(Clone, Debug)]
pub struct FiniteViscoMembrane<'a> {
    patch: &'a MembranePatch,
    material: ViscoelasticLawField<'a>,
    displacement_m: Arc<[[f64; 3]]>,
    states: Arc<[ViscoMembraneState]>,
    first_piola_pa: Arc<[[f64; 9]]>,
    pressure_pa: Arc<[Option<f64>]>,
    resultants_n_m: Arc<[[f64; 3]]>,
    stored_energy_j: f64,
    time_s: f64,
}

/// Complete proposed endpoint; dropping it leaves the previous membrane intact.
#[derive(Debug)]
pub struct FiniteViscoMembraneStep<'a> {
    pub membrane: FiniteViscoMembrane<'a>,
    /// `strain_energy_j` is the frozen material's incremental potential here.
    /// Physical storage and relaxation loss belong to `energy`.
    pub equilibrium: MembraneSolution,
    pub energy: FiniteViscoMembraneEnergy,
    /// Physical loss in each original triangle during this interval [J].
    /// Retained before cumulative-ledger rounding, for conservative heat transfer.
    /// Neither endpoint-work defects nor numerical kinetic loss enter this array.
    pub element_dissipated_j: Vec<f64>,
}

/// Original material/geometry/solver refusal or an explicit interval-work gate.
#[derive(Debug)]
pub enum FiniteViscoMembraneError {
    Native(MembraneError),
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<MembraneError> for FiniteViscoMembraneError {
    fn from(error: MembraneError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for FiniteViscoMembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "viscoelastic membrane {quantity}: {value_j:e} J exceeds {maximum_abs_j:e} J",
            ),
        }
    }
}
impl std::error::Error for FiniteViscoMembraneError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Native(error) => Some(error),
            Self::EnergyGate { .. } => None,
        }
    }
}

fn invalid(reason: &'static str) -> FiniteViscoMembraneError {
    MembraneError::Invalid(reason).into()
}

fn material_error(element: usize, source: SheetError) -> FiniteViscoMembraneError {
    MembraneError::Material { element, source }.into()
}

// Both virgin initialization and retained-reference changes bind the same
// source equilibrium and fixed material frame. No constitutive state is reset
// merely to recheck that immutable binding.
fn admit_sheet_binding<'p>(
    patch: &'p MembranePatch,
    index: usize,
    law: &ViscoHyperelastic,
) -> Result<&'p HyperelasticMembrane, FiniteViscoMembraneError> {
    let MembraneLaw::Hyperelastic(sheet) = &patch.elements[index].law else {
        return Err(invalid(
            "native viscoelastic memory requires the original hyperelastic sheet",
        ));
    };
    if law.branches.len() > MAX_FINITE_STRAIN_BRANCHES {
        return Err(invalid(
            "membrane branch count exceeds the native material cap",
        ));
    }
    let g = sheet.reference_stretches();
    let reference = [g[0], 0.0, 0.0, 0.0, g[1], 0.0, 0.0, 0.0, g[2]];
    if law.equilibrium != *sheet.model().law()
        || *law.reference_deformation().deformation_gradient() != reference
    {
        return Err(invalid(
            "viscoelastic equilibrium and natural reference must match the original sheet",
        ));
    }
    Ok(sheet.as_ref())
}

fn energy_gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteViscoMembraneError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FiniteViscoMembraneError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

struct Endpoint {
    states: Vec<ViscoMembraneState>,
    first_piola_pa: Vec<[f64; 9]>,
    pressure_pa: Vec<Option<f64>>,
    responses: Vec<SheetResponse>,
    element_dissipated_j: Vec<f64>,
    energy: FiniteViscoMembraneEnergy,
}

impl Endpoint {
    fn new(count: usize, stored_before_j: f64) -> Self {
        Self {
            states: Vec::with_capacity(count),
            first_piola_pa: Vec::with_capacity(count),
            pressure_pa: Vec::with_capacity(count),
            responses: Vec::with_capacity(count),
            element_dissipated_j: Vec::with_capacity(count),
            energy: FiniteViscoMembraneEnergy {
                stored_before_j,
                ..Default::default()
            },
        }
    }

    fn push(&mut self, area_m2: f64, step: ViscoMembraneStepResponse) {
        self.energy.stored_after_j += area_m2 * step.stored_energy_j_m2;
        self.energy.dissipated_j += area_m2 * step.dissipated_j_m2;
        self.element_dissipated_j
            .push(area_m2 * step.dissipated_j_m2);
        self.energy.native_loading_work_j += area_m2 * step.loading_work_j_m2;
        self.energy.incremental_potential_j += area_m2 * step.sheet.energy_j_m2;
        self.states.push(step.state);
        self.first_piola_pa.push(step.first_piola_pa);
        self.pressure_pa.push(step.pressure_pa);
        self.responses.push(step.sheet);
    }
}

impl<'a> FiniteViscoMembrane<'a> {
    /// Place virgin branch histories at zero displacement with a native
    /// zero-duration elastic jump. A prescribed natural reference may make
    /// this initial configuration prestressed; no equilibrium is presumed.
    pub fn try_new(
        patch: &'a MembranePatch,
        material: ViscoelasticLawField<'a>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FiniteViscoMembraneError> {
        Self::try_new_from(
            patch,
            material,
            &vec![[0.0; 3]; patch.mesh.nodes.len()],
            cancelled,
        )
    }

    /// Initialize virgin metrics at an explicit installed configuration. The
    /// native initial loading work is retained without elapsed time or loss.
    /// All triangles must retain their source equilibrium law, natural reference
    /// and volume treatment; supplied Maxwell branches remain explicit inputs.
    pub fn try_new_from(
        patch: &'a MembranePatch,
        material: ViscoelasticLawField<'a>,
        displacement_m: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        if !finite_vectors(displacement_m, patch.mesh.nodes.len()) {
            return Err(invalid("finite original-node displacements required"));
        }
        if let ViscoelasticLawField::PerElement(laws) = material
            && laws.len() != patch.elements.len()
        {
            return Err(invalid(
                "one original viscoelastic law per membrane triangle required",
            ));
        }
        let mut endpoint = Endpoint::new(patch.elements.len(), 0.0);
        for (index, element) in patch.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let law = material.at(index);
            let sheet = admit_sheet_binding(patch, index, law)?;
            let strain =
                patch.admit_kinematics(index, patch.kinematics(element, displacement_m))?;
            let virgin = ViscoMembraneState::initial(law, sheet.volumetric_treatment());
            let frozen = FrozenViscoMembrane::new(
                law,
                &virgin,
                0.0,
                sheet.reference_thickness_m(),
                sheet.settings(),
            )
            .map_err(|source| material_error(index, source))?;
            let step = frozen
                .step_response(strain)
                .map_err(|source| material_error(index, source))?;
            endpoint.push(element.area, step);
        }
        if !endpoint.energy.stored_after_j.is_finite() {
            return Err(invalid("unrepresentable initial membrane storage"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            patch,
            material,
            displacement_m: displacement_m.into(),
            states: endpoint.states.into(),
            first_piola_pa: endpoint.first_piola_pa.into(),
            pressure_pa: endpoint.pressure_pa.into(),
            resultants_n_m: endpoint
                .responses
                .iter()
                .map(|r| r.resultants_n_m)
                .collect::<Vec<_>>()
                .into(),
            stored_energy_j: endpoint.energy.stored_after_j,
            time_s: 0.0,
        })
    }

    /// Original source/geometry/mass owner, borrowed without reconstruction.
    pub fn patch(&self) -> &'a MembranePatch {
        self.patch
    }
    pub fn material(&self) -> ViscoelasticLawField<'a> {
        self.material
    }
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    /// Native histories in original triangle order and fixed material axes.
    pub fn triangle_states(&self) -> &[ViscoMembraneState] {
        &self.states
    }
    /// Physical Piola stress in each canonical Cholesky representative's spatial
    /// frame, not global Cartesian axes. Exact-volume stresses include pressure.
    pub fn triangle_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    /// Exact-volume pressure multipliers; finite-compliance entries are `None`.
    pub fn triangle_pressure_pa(&self) -> &[Option<f64>] {
        &self.pressure_pa
    }
    /// Objective second-Piola resultants in original warp/weft axes [N/m].
    pub fn triangle_resultants_n_m(&self) -> &[[f64; 3]] {
        &self.resultants_n_m
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    fn sheet(&self, index: usize) -> &HyperelasticMembrane {
        match &self.patch.elements[index].law {
            MembraneLaw::Hyperelastic(sheet) => sheet,
            _ => unreachable!("constructor admitted immutable hyperelastic sheets"),
        }
    }

    fn freeze(
        &self,
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FrozenViscoMembrane<'_>>, FiniteViscoMembraneError> {
        let mut frozen = Vec::with_capacity(self.states.len());
        for (index, state) in self.states.iter().enumerate() {
            poll(cancelled)?;
            let sheet = self.sheet(index);
            frozen.push(
                FrozenViscoMembrane::new(
                    self.material.at(index),
                    state,
                    dt_s,
                    sheet.reference_thickness_m(),
                    sheet.settings(),
                )
                .map_err(|source| material_error(index, source))?,
            );
        }
        Ok(frozen)
    }

    /// Propose one actual physical interval under endpoint dead forces and
    /// prescribed displacements. A zero interval is an elastic loading jump;
    /// positive intervals relax the same accepted histories during every trial.
    /// No internal load ramp advances material memory more than once.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        prescribed: &[(usize, f64)],
        settings: MembraneSolveSettings,
        controls: FiniteViscoStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneStep<'a>, FiniteViscoMembraneError> {
        self.advance_with_pressure(
            dt_s,
            dead_load_n,
            0.0,
            prescribed,
            settings,
            controls,
            cancelled,
        )
    }

    /// The same retained-history interval with a held signed pressure jump on
    /// the CURRENT membrane surface. The original fixed-perimeter admission,
    /// exact follower tangent and pressure-work owner remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure_pa: f64,
        prescribed: &[(usize, f64)],
        settings: MembraneSolveSettings,
        controls: FiniteViscoStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneStep<'a>, FiniteViscoMembraneError> {
        self.advance_loaded(
            dt_s,
            dead_load_n,
            pressure_pa,
            prescribed,
            &self.displacement_m,
            settings,
            controls,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_loaded(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure_pa: f64,
        prescribed: &[(usize, f64)],
        initial_displacement_m: &[[f64; 3]],
        settings: MembraneSolveSettings,
        controls: FiniteViscoStepControls,
        inertia: Option<&dyn MembraneInertialPotential>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneStep<'a>, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        let time_s = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s < 0.0
            || !time_s.is_finite()
            || (dt_s > 0.0 && time_s <= self.time_s)
            || ![
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
            ]
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "representable nonnegative interval and finite nonnegative work allowances required",
            ));
        }
        let frozen = self.freeze(dt_s, &mut cancelled)?;
        // Krylov repeatedly applies B^T C B at one Newton state. Cache only an
        // exactly repeated local strain within this immutable history/interval,
        // so those products do not repeat the transverse spectral return map.
        // Failed responses are never cached; rejected global trials cannot
        // advance either the history or the physical interval.
        let cache = RefCell::new(vec![None::<([f64; 3], SheetResponse)>; frozen.len()]);
        let response = |index: usize, strain| {
            if let Some((previous, value)) = cache.borrow()[index]
                && previous == strain
            {
                return Ok(value);
            }
            let value = frozen[index].response(strain)?;
            cache.borrow_mut()[index] = Some((strain, value));
            Ok(value)
        };
        let equilibrium = self.patch.solve_with_response(
            &response,
            inertia,
            dead_load_n,
            pressure_pa,
            prescribed,
            initial_displacement_m,
            settings,
            &mut cancelled,
        )?;
        let mut endpoint = Endpoint::new(self.states.len(), self.stored_energy_j);
        for (index, element) in self.patch.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let strain = self.patch.admit_kinematics(
                index,
                self.patch.kinematics(element, &equilibrium.displacement_m),
            )?;
            let step = frozen[index]
                .step_response(strain)
                .map_err(|source| material_error(index, source))?;
            endpoint.push(element.area, step);
        }
        // Reuse the original assembly with accepted physical point responses,
        // avoiding another constitutive return map or a duplicate nodal force law.
        let physical = self.patch.evaluate_with_response(
            &equilibrium.displacement_m,
            &|index, _| Ok(endpoint.responses[index]),
            &mut cancelled,
        )?;
        let pressure = if pressure_pa == 0.0 {
            Vec::new()
        } else {
            self.patch
                .pressure_load_n(&equilibrium.displacement_m, pressure_pa, &mut cancelled)?
        };
        let inertial = if let Some(inertia) = inertia {
            let value = inertia.evaluate(&equilibrium.displacement_m, &mut cancelled)?;
            if !finite_vectors(&value.internal_gradient_n, self.displacement_m.len()) {
                return Err(invalid("finite complete endpoint inertia required"));
            }
            value.internal_gradient_n
        } else {
            Vec::new()
        };
        let mut energy = endpoint.energy;
        for (node, (u, old)) in equilibrium
            .displacement_m
            .iter()
            .zip(self.displacement_m.iter())
            .enumerate()
        {
            poll(&mut cancelled)?;
            for axis in 0..3 {
                let increment = u[axis] - old[axis];
                energy.internal_endpoint_work_j +=
                    physical.internal_gradient_n[node][axis] * increment;
                energy.dead_load_endpoint_work_j += dead_load_n[node][axis] * increment;
                if !pressure.is_empty() {
                    energy.pressure_endpoint_work_j += pressure[node][axis] * increment;
                }
                if !inertial.is_empty() {
                    energy.inertial_endpoint_work_j += inertial[node][axis] * increment;
                }
            }
        }
        for &(dof, _) in prescribed {
            let increment = equilibrium.displacement_m[dof / 3][dof % 3]
                - self.displacement_m[dof / 3][dof % 3];
            energy.support_endpoint_work_j += equilibrium.residual_n[dof / 3][dof % 3] * increment;
        }
        if pressure_pa != 0.0 {
            energy.pressure_path_work_j = self.patch.pressure_work_j(
                &self.displacement_m,
                &equilibrium.displacement_m,
                pressure_pa,
                &mut cancelled,
            )?;
        }
        let physical_increment =
            energy.stored_after_j - energy.stored_before_j + energy.dissipated_j;
        energy.material_balance_error_j = energy.native_loading_work_j - physical_increment;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.pressure_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j
            - energy.inertial_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j - physical_increment;
        energy_gate(
            "native material balance",
            energy.material_balance_error_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "assembled incremental potential",
            energy.incremental_potential_j - equilibrium.strain_energy_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "support/load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "endpoint-work integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if !energy.pressure_path_work_j.is_finite() || energy.dissipated_j < 0.0 {
            return Err(invalid(
                "unrepresentable pressure path work or negative native relaxation loss",
            ));
        }
        poll(&mut cancelled)?;
        let membrane = Self {
            patch: self.patch,
            material: self.material,
            displacement_m: equilibrium.displacement_m.clone().into(),
            states: endpoint.states.into(),
            first_piola_pa: endpoint.first_piola_pa.into(),
            pressure_pa: endpoint.pressure_pa.into(),
            resultants_n_m: endpoint
                .responses
                .iter()
                .map(|r| r.resultants_n_m)
                .collect::<Vec<_>>()
                .into(),
            stored_energy_j: energy.stored_after_j,
            time_s,
        };
        Ok(FiniteViscoMembraneStep {
            membrane,
            equilibrium,
            energy,
            element_dissipated_j: endpoint.element_dissipated_j,
        })
    }
}
