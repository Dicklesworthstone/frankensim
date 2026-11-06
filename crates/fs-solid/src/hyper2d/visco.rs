//! Persistent finite-strain viscoelastic histories for native thin solids.
//!
//! Each original integration point owns its history. One native Newton solve
//! borrows those histories at a single physical interval; only its complete,
//! admitted endpoint becomes the next solid. Original area times original
//! thickness weights physical storage, native loading work and relaxation loss.
//! The solver's incremental potential and signed endpoint-work error remain
//! separate accounts. Exact elastic incompressibility uses bulk-free histories.
//!
//! This is quasistatic in-plane P1/Q1 mechanics with local transverse traction
//! relaxation. It adds no inertia, shell bending, wrinkling or caloric law.
//! Nonuniform histories are retained at all four Q1 points; P1 has constant
//! deformation and one centroid. Time/load refinement controls the reported
//! work defect, including transverse motion during a held in-plane interval.

use super::{
    HyperProblem, HyperelasticMaterials, NewtonReport, NewtonSettings, PlaneMaterial,
    PlaneMaterialPoint, material_point_kinematics,
};
use crate::{Mesh2, Patch, SolidError};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::incompressible_plane_stress::IncompressiblePlaneStressPotential;
use fs_material::hyper::mixed::RemainderResponse;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::visco_hyper::{
    FrozenViscoHyperelastic, FrozenViscoRemainder, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic,
    ViscoRemainderState, ViscoState,
};
use std::sync::Arc;

pub use crate::hyper3::visco::{
    FiniteViscoStepControls, FiniteViscoTetEnergy as FiniteViscoPlaneEnergy, ViscoelasticLawField,
};

/// Prescribed reference-strain changes with retained quadrature histories.
pub mod reference_change;
pub use reference_change::{FiniteViscoPlaneReferenceChange, FiniteViscoPlaneReferenceEnergy};

/// The transverse constraint used throughout initialization and every interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteViscoPlaneMode {
    /// The complete finite-compliance native constitutive law.
    FiniteBulk,
    /// Exact elastic volume, with original bulk-free histories and a pressure
    /// multiplier. The native remainder owner's split admission still applies.
    Incompressible,
}

#[derive(Debug)]
enum Histories {
    Full(Vec<ViscoState>),
    Remainder(Vec<ViscoRemainderState>),
}

impl Histories {
    fn new(mode: FiniteViscoPlaneMode, capacity: usize) -> Self {
        match mode {
            FiniteViscoPlaneMode::FiniteBulk => Self::Full(Vec::with_capacity(capacity)),
            FiniteViscoPlaneMode::Incompressible => Self::Remainder(Vec::with_capacity(capacity)),
        }
    }

    fn at(&self, index: usize) -> PointHistory<'_> {
        match self {
            Self::Full(states) => PointHistory::Full(&states[index]),
            Self::Remainder(states) => PointHistory::Remainder(&states[index]),
        }
    }
}

#[derive(Clone, Copy)]
enum PointHistory<'a> {
    Full(&'a ViscoState),
    Remainder(&'a ViscoRemainderState),
}

impl PointHistory<'_> {
    fn gradient(self) -> [f64; 9] {
        match self {
            Self::Full(state) => state.deformation_gradient,
            Self::Remainder(state) => state.deformation_gradient,
        }
    }
}

enum ProposedHistory {
    Full(ViscoState),
    Remainder(ViscoRemainderState),
}

struct PointUpdate {
    history: ProposedHistory,
    first_piola_pa: [f64; 9],
    pressure_pa: Option<f64>,
    stored_pa: f64,
    loading_work_pa: f64,
    dissipated_pa: f64,
    incremental_potential_pa: f64,
}

fn update_point(
    law: &ViscoHyperelastic,
    accepted: PointHistory<'_>,
    in_plane: &[f64; 4],
    dt_s: f64,
    use_accepted_seed: bool,
    local: PlaneStressSettings,
) -> Result<PointUpdate, MaterialError> {
    let f = accepted.gradient();
    match accepted {
        PointHistory::Full(state) => {
            let response = law.step_plane_stress_response(
                state,
                in_plane,
                dt_s,
                use_accepted_seed.then_some([f[2], f[5], f[8]]),
                local,
            )?;
            let material = response.material;
            Ok(PointUpdate {
                stored_pa: material.state.stored,
                history: ProposedHistory::Full(material.state),
                first_piola_pa: material.first_piola_pa,
                pressure_pa: None,
                loading_work_pa: material.loading_work_pa,
                dissipated_pa: material.dissipated_pa,
                incremental_potential_pa: material.incremental_potential_pa,
            })
        }
        PointHistory::Remainder(state) => {
            let response = law.step_incompressible_plane_stress_response(
                state,
                in_plane,
                dt_s,
                use_accepted_seed.then_some([f[2], f[5]]),
                local,
            )?;
            let material = response.material;
            Ok(PointUpdate {
                stored_pa: material.state.stored,
                history: ProposedHistory::Remainder(material.state),
                first_piola_pa: response.plane_stress.full_piola,
                pressure_pa: Some(response.plane_stress.pressure_pa),
                loading_work_pa: material.loading_work_pa,
                dissipated_pa: material.dissipated_pa,
                incremental_potential_pa: material.incremental_potential_pa,
            })
        }
    }
}

struct Endpoint {
    histories: Histories,
    first_piola_pa: Vec<[f64; 9]>,
    pressure_pa: Vec<f64>,
    energy: FiniteViscoPlaneEnergy,
}

impl Endpoint {
    fn new(mode: FiniteViscoPlaneMode, count: usize, stored_before_j: f64) -> Self {
        Self {
            histories: Histories::new(mode, count),
            first_piola_pa: Vec::with_capacity(count),
            pressure_pa: if mode == FiniteViscoPlaneMode::Incompressible {
                Vec::with_capacity(count)
            } else {
                Vec::new()
            },
            energy: FiniteViscoPlaneEnergy {
                stored_before_j,
                ..Default::default()
            },
        }
    }

    fn push(&mut self, volume_m3: f64, update: PointUpdate) -> Result<(), FiniteViscoPlaneError> {
        if !(volume_m3.is_finite() && volume_m3 > 0.0) {
            return Err(invalid(
                "positive representable original point volume required",
            ));
        }
        self.energy.stored_after_j += volume_m3 * update.stored_pa;
        self.energy.native_loading_work_j += volume_m3 * update.loading_work_pa;
        self.energy.dissipated_j += volume_m3 * update.dissipated_pa;
        self.energy.incremental_potential_j += volume_m3 * update.incremental_potential_pa;
        if ![
            self.energy.stored_after_j,
            self.energy.native_loading_work_j,
            self.energy.dissipated_j,
            self.energy.incremental_potential_j,
        ]
        .iter()
        .all(|value| value.is_finite())
            || self.energy.dissipated_j < 0.0
        {
            return Err(invalid(
                "finite material accounts and nonnegative relaxation loss required",
            ));
        }
        self.first_piola_pa.push(update.first_piola_pa);
        if let Some(pressure) = update.pressure_pa {
            self.pressure_pa.push(pressure);
        }
        match (&mut self.histories, update.history) {
            (Histories::Full(states), ProposedHistory::Full(state)) => states.push(state),
            (Histories::Remainder(states), ProposedHistory::Remainder(state)) => states.push(state),
            _ => return Err(invalid("inconsistent point-history formulation")),
        }
        Ok(())
    }
}

// The native point-field interface includes an unconstrained potential as well
// as the exact-volume capability. Keep this adapter private to the exact owner:
// its mathematical potential is the original remainder, never a full bulk law.
// No bulk energy or pressure work is reconstructed to satisfy the interface.
struct ExactPointPotential<'a>(FrozenViscoRemainder<'a>);

impl FiniteStrainPotential for ExactPointPotential<'_> {
    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        self.0
            .incompressible_remainder_stress_and_potential(f)
            .map_or(f64::NAN, |(_, potential)| potential)
    }
    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.0
            .incompressible_remainder_stress_and_potential(f)
            .map(|(piola, _)| piola)
    }
    fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.0
            .incompressible_remainder_response(f)
            .map(|response| response.tangent_pa)
    }
}

impl IncompressiblePlaneStressPotential for ExactPointPotential<'_> {
    fn incompressible_reference_volume_ratio(&self) -> f64 {
        self.0.incompressible_reference_volume_ratio()
    }
    fn incompressible_remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        self.0.incompressible_remainder_stress_and_potential(f)
    }
    fn incompressible_remainder_response(
        &self,
        f: &[f64; 9],
    ) -> Result<RemainderResponse, MaterialError> {
        self.0.incompressible_remainder_response(f)
    }
}
impl PlaneMaterial for ExactPointPotential<'_> {}

enum FrozenViews<'a> {
    Full(Vec<FrozenViscoHyperelastic<'a>>),
    Remainder(Vec<ExactPointPotential<'a>>),
}

struct FrozenPointField<'a> {
    views: FrozenViews<'a>,
    histories: &'a Histories,
    offsets: &'a [usize],
}

impl FrozenPointField<'_> {
    fn index(&self, element: usize, local_point: usize) -> Option<usize> {
        let start = *self.offsets.get(element)?;
        let end = *self.offsets.get(element + 1)?;
        (local_point < end - start).then_some(start + local_point)
    }
}

impl HyperelasticMaterials for FrozenPointField<'_> {
    fn material_count(&self) -> Option<usize> {
        Some(self.offsets.len() - 1)
    }
    fn material_at(&self, _: usize) -> Option<&dyn PlaneMaterial> {
        None
    }
    fn material_point_count(&self, element: usize) -> Option<usize> {
        self.offsets
            .get(element + 1)
            .zip(self.offsets.get(element))
            .map(|(end, start)| end - start)
    }
    fn material_at_point(&self, element: usize, local_point: usize) -> Option<&dyn PlaneMaterial> {
        let index = self.index(element, local_point)?;
        match &self.views {
            FrozenViews::Full(views) => views.get(index).map(|view| view as &dyn PlaneMaterial),
            FrozenViews::Remainder(views) => {
                views.get(index).map(|view| view as &dyn PlaneMaterial)
            }
        }
    }
    fn transverse_seed_at(&self, element: usize, local_point: usize) -> Option<[f64; 3]> {
        let f = self
            .histories
            .at(self.index(element, local_point)?)
            .gradient();
        Some([f[2], f[5], f[8]])
    }
}

/// Original mesh/model bindings and a complete immutable accepted trajectory.
#[derive(Debug)]
pub struct FiniteViscoPlaneSolid<'a> {
    mesh: &'a Mesh2,
    material: ViscoelasticLawField<'a>,
    mode: FiniteViscoPlaneMode,
    reference_thickness_m: f64,
    local: PlaneStressSettings,
    displacement_m: Arc<[[f64; 2]]>,
    points: Arc<[PlaneMaterialPoint]>,
    offsets: Arc<[usize]>,
    histories: Histories,
    first_piola_pa: Vec<[f64; 9]>,
    pressure_pa: Vec<f64>,
    stored_energy_j: f64,
    time_s: f64,
}

/// A complete candidate; dropping it leaves the previous body and time intact.
#[derive(Debug)]
pub struct FiniteViscoPlaneStep<'a> {
    pub solid: FiniteViscoPlaneSolid<'a>,
    /// Original native nonlinear convergence report.
    pub equilibrium: NewtonReport,
    /// Shared native physical/algorithmic accounts. Inertial work is zero in
    /// this quasistatic owner. Pressure constraint work is not invented heat.
    pub energy: FiniteViscoPlaneEnergy,
    /// Physical support reactions [N], in sorted global in-plane DOF order.
    pub reactions_n: Vec<(usize, f64)>,
}

#[derive(Debug)]
pub enum FiniteViscoPlaneError {
    Native(SolidError),
    Material {
        element: usize,
        local_point: usize,
        source: MaterialError,
    },
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<SolidError> for FiniteViscoPlaneError {
    fn from(error: SolidError) -> Self {
        Self::Native(error)
    }
}

impl core::fmt::Display for FiniteViscoPlaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Material {
                element,
                local_point,
                source,
            } => write!(
                f,
                "viscoelastic element {element}, point {local_point}: {source}"
            ),
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
impl std::error::Error for FiniteViscoPlaneError {}

fn invalid(what: &str) -> FiniteViscoPlaneError {
    SolidError::InvalidInput { what: what.into() }.into()
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), FiniteViscoPlaneError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled.into())
}
fn material_error(point: &PlaneMaterialPoint, source: MaterialError) -> FiniteViscoPlaneError {
    FiniteViscoPlaneError::Material {
        element: point.element,
        local_point: point.local_point,
        source,
    }
}
fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteViscoPlaneError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FiniteViscoPlaneError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

impl<'a> FiniteViscoPlaneSolid<'a> {
    /// Initialize at zero displacement, with virgin branch metrics and locally
    /// relaxed transverse deformation. This can be a prestressed configuration.
    pub fn try_new(
        mesh: &'a Mesh2,
        material: ViscoelasticLawField<'a>,
        mode: FiniteViscoPlaneMode,
        reference_thickness_m: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoPlaneError> {
        let displacement = vec![[0.0; 2]; mesh.nodes.len()];
        Self::try_new_from(
            mesh,
            material,
            mode,
            reference_thickness_m,
            &displacement,
            local,
            cx,
        )
    }

    /// Place virgin histories at an explicit in-plane configuration through
    /// native zero-duration material steps. Initial placement work is retained
    /// in each history, without relaxation loss or a fictitious physical age.
    /// The initial body need not be in global equilibrium; its next interval
    /// uses the caller's actual supports and loads.
    pub fn try_new_from(
        mesh: &'a Mesh2,
        material: ViscoelasticLawField<'a>,
        mode: FiniteViscoPlaneMode,
        reference_thickness_m: f64,
        initial_displacement_m: &[[f64; 2]],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoPlaneError> {
        checkpoint(cx)?;
        if !(reference_thickness_m.is_finite() && reference_thickness_m > 0.0) {
            return Err(invalid("positive finite original thickness required"));
        }
        if let ViscoelasticLawField::PerElement(laws) = material
            && laws.len() != mesh.elems.len()
        {
            return Err(invalid(
                "one original viscoelastic law per element required",
            ));
        }
        let points = material_point_kinematics(mesh, initial_displacement_m, cx)?;
        let mut offsets = vec![0usize; mesh.elems.len() + 1];
        let mut endpoint = Endpoint::new(mode, points.len(), 0.0);
        for point in &points {
            checkpoint(cx)?;
            offsets[point.element + 1] += 1;
            let law = material.at(point.element);
            if law.branches.len() > MAX_FINITE_STRAIN_BRANCHES {
                return Err(invalid("viscoelastic branch count exceeds the native cap"));
            }
            let update = match mode {
                FiniteViscoPlaneMode::FiniteBulk => update_point(
                    law,
                    PointHistory::Full(&law.initial_state()),
                    &point.in_plane,
                    0.0,
                    false,
                    local,
                ),
                FiniteViscoPlaneMode::Incompressible => update_point(
                    law,
                    PointHistory::Remainder(&law.initial_remainder_state()),
                    &point.in_plane,
                    0.0,
                    false,
                    local,
                ),
            }
            .map_err(|source| material_error(point, source))?;
            endpoint.push(point.reference_area_m2 * reference_thickness_m, update)?;
        }
        for element in 0..mesh.elems.len() {
            offsets[element + 1] += offsets[element];
        }
        checkpoint(cx)?;
        Ok(Self {
            mesh,
            material,
            mode,
            reference_thickness_m,
            local,
            displacement_m: initial_displacement_m.into(),
            points: points.into(),
            offsets: offsets.into(),
            histories: endpoint.histories,
            first_piola_pa: endpoint.first_piola_pa,
            pressure_pa: endpoint.pressure_pa,
            stored_energy_j: endpoint.energy.stored_after_j,
            time_s: 0.0,
        })
    }

    pub fn reference_mesh(&self) -> &'a Mesh2 {
        self.mesh
    }
    pub fn material(&self) -> ViscoelasticLawField<'a> {
        self.material
    }
    pub fn mode(&self) -> FiniteViscoPlaneMode {
        self.mode
    }
    pub fn reference_thickness_m(&self) -> f64 {
        self.reference_thickness_m
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    /// Element-major point order shared by every native operator and history.
    pub fn material_points(&self) -> &[PlaneMaterialPoint] {
        &self.points
    }
    pub fn point_states(&self) -> Option<&[ViscoState]> {
        match &self.histories {
            Histories::Full(states) => Some(states),
            _ => None,
        }
    }
    pub fn point_remainder_states(&self) -> Option<&[ViscoRemainderState]> {
        match &self.histories {
            Histories::Remainder(states) => Some(states),
            _ => None,
        }
    }
    /// Full physical first Piola, including exact pressure where applicable.
    pub fn point_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    /// Compression-positive pressure of the exact plane-stress reduction.
    pub fn point_pressure_pa(&self) -> Option<&[f64]> {
        (self.mode == FiniteViscoPlaneMode::Incompressible).then_some(&self.pressure_pa)
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    fn freeze(
        &self,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<FrozenPointField<'_>, FiniteViscoPlaneError> {
        let views = match &self.histories {
            Histories::Full(states) => {
                let mut views = Vec::with_capacity(states.len());
                for (point, state) in self.points.iter().zip(states) {
                    checkpoint(cx)?;
                    views.push(
                        self.material
                            .at(point.element)
                            .freeze(state, dt_s)
                            .map_err(|source| material_error(point, source))?,
                    );
                }
                FrozenViews::Full(views)
            }
            Histories::Remainder(states) => {
                let mut views = Vec::with_capacity(states.len());
                for (point, state) in self.points.iter().zip(states) {
                    checkpoint(cx)?;
                    views.push(ExactPointPotential(
                        self.material
                            .at(point.element)
                            .freeze_remainder(state, dt_s)
                            .map_err(|source| material_error(point, source))?,
                    ));
                }
                FrozenViews::Remainder(views)
            }
        };
        Ok(FrozenPointField {
            views,
            histories: &self.histories,
            offsets: &self.offsets,
        })
    }

    /// Solve and propose one physical interval with absolute endpoint supports
    /// [m] and dead nominal boundary tractions [Pa]. Exactly one native load
    /// endpoint is required; time refinement uses successive accepted bodies.
    /// Zero duration admits an instantaneous loading change without relaxation.
    ///
    /// Every local/global trial reads the same accepted histories. Post-solve
    /// storage, work and loss gates must all pass before returning a successor.
    /// A held in-plane displacement can still relax thickness/shear, so native
    /// loading work need not vanish; its endpoint-work difference is numerical
    /// error and is never exported as heat.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        settings: NewtonSettings,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoPlaneStep<'a>, FiniteViscoPlaneError> {
        checkpoint(cx)?;
        let time_s = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s < 0.0
            || !time_s.is_finite()
            || (dt_s > 0.0 && time_s <= self.time_s)
            || settings.load_steps != 1
            || ![
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
            ]
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "one endpoint, representable nonnegative interval and finite work allowances required",
            ));
        }
        let frozen = self.freeze(dt_s, cx)?;
        let problem = HyperProblem {
            mesh: self.mesh,
            material: &frozen,
            dirichlet: dirichlet.to_vec(),
            traction: traction.to_vec(),
            settings,
        };
        let (displacement_m, equilibrium) = match self.mode {
            FiniteViscoPlaneMode::FiniteBulk => {
                problem.solve_plane_stress_from(&self.displacement_m, self.local, cx)
            }
            FiniteViscoPlaneMode::Incompressible => {
                problem.solve_incompressible_plane_stress_from(&self.displacement_m, self.local, cx)
            }
        }?;
        let points = material_point_kinematics(self.mesh, &displacement_m, cx)?;
        if points.len() != self.points.len() {
            return Err(invalid("accepted integration-point layout changed"));
        }
        let mut endpoint = Endpoint::new(self.mode, points.len(), self.stored_energy_j);
        for (index, (point, previous)) in points.iter().zip(self.points.iter()).enumerate() {
            checkpoint(cx)?;
            if point.element != previous.element
                || point.local_point != previous.local_point
                || point.reference_area_m2 != previous.reference_area_m2
            {
                return Err(invalid(
                    "original integration-point identity or volume changed",
                ));
            }
            let update = update_point(
                self.material.at(point.element),
                self.histories.at(index),
                &point.in_plane,
                dt_s,
                true,
                self.local,
            )
            .map_err(|source| material_error(point, source))?;
            endpoint.push(point.reference_area_m2 * self.reference_thickness_m, update)?;
        }
        let u = displacement_m.as_flattened();
        let old = self.displacement_m.as_flattened();
        // load=0 retains internal force only, without another Hessian assembly.
        let (internal, potential) = match self.mode {
            FiniteViscoPlaneMode::FiniteBulk => (
                problem.plane_stress_unconstrained_residual(u, 0.0, self.local, cx)?,
                problem.plane_stress_potential_energy(u, 0.0, self.local, cx)?,
            ),
            FiniteViscoPlaneMode::Incompressible => (
                problem
                    .incompressible_plane_stress_unconstrained_residual(u, 0.0, self.local, cx)?,
                problem.incompressible_plane_stress_potential_energy(u, 0.0, self.local, cx)?,
            ),
        };
        let mut external = vec![0.0; u.len()];
        problem.for_traction(1.0, |dof, value| external[dof] += value)?;
        let mut energy = endpoint.energy;
        for dof in 0..u.len() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            let increment = u[dof] - old[dof];
            energy.internal_endpoint_work_j +=
                self.reference_thickness_m * internal[dof] * increment;
            energy.dead_load_endpoint_work_j +=
                self.reference_thickness_m * external[dof] * increment;
        }
        let mut reactions_n = Vec::new();
        for &dof in problem.fixed_dofs(1.0)?.keys() {
            if reactions_n.len() % 256 == 0 {
                checkpoint(cx)?;
            }
            let reaction = self.reference_thickness_m * (internal[dof] - external[dof]);
            if !reaction.is_finite() {
                return Err(invalid("support reaction is unrepresentable"));
            }
            energy.support_endpoint_work_j += reaction * (u[dof] - old[dof]);
            reactions_n.push((dof, reaction));
        }
        let stored_change = energy.stored_after_j - energy.stored_before_j;
        energy.material_balance_error_j =
            energy.native_loading_work_j - stored_change - energy.dissipated_j;
        energy.endpoint_work_defect_j =
            energy.internal_endpoint_work_j - stored_change - energy.dissipated_j;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        gate(
            "native material balance",
            energy.material_balance_error_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        gate(
            "assembled incremental potential",
            self.reference_thickness_m * potential - energy.incremental_potential_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        gate(
            "support/dead-load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        gate(
            "signed endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteViscoPlaneStep {
            solid: Self {
                mesh: self.mesh,
                material: self.material,
                mode: self.mode,
                reference_thickness_m: self.reference_thickness_m,
                local: self.local,
                displacement_m: displacement_m.into(),
                points: points.into(),
                offsets: Arc::clone(&self.offsets),
                histories: endpoint.histories,
                first_piola_pa: endpoint.first_piola_pa,
                pressure_pa: endpoint.pressure_pa,
                stored_energy_j: energy.stored_after_j,
                time_s,
            },
            equilibrium,
            energy,
            reactions_n,
        })
    }
}
