//! Nonlinear reversible induction on the original planar field/conductor owners.
//!
//! Backward Euler solves H(curl A_next) + C_sigma delta A/dt = J_next in
//! the weak P1 form. The static material Newton/FEEC/checked-CG owner solves
//! this positive variational problem; the original eddy owner supplies the
//! floating conductor projection and exact spatial sigma E^2 heating.
//!
//! Energy is NOT H*B/2 or midpoint-H work. For the source-integrated convex
//! energy psi, D = H_next.(B_next-B_old) - (psi_next-psi_old) is the nonnegative
//! backward-Euler numerical dissipation. Delta psi + Joule heat + D = source
//! work for a fixed law. A source-pinned linear mu(T) refresh adds explicit
//! parameter work P=psi_next(B_old)-psi_old(B_old), so delta psi + heat + D
//! equals source work + P. D is never material heat; a thermal owner must close
//! P and the field entropy. Explicit source-scaled H(B,T) may also refresh;
//! arbitrary nonlinear/recoil source changes still refuse. This is first-order
//! fixed-geometry isotropic reversible MQS, without hysteresis or remanence.

use super::{
    Conduction, EddyCellHeat, EddyConductor, EddyControls, EddyError, EddyHeat, EddySection,
    admit_algebraic, algebraic_blocks, bind, controls, finite, heat_cells, load, poll, positive,
    regions,
};
use crate::planar::{
    MagneticError, PlanarMagnetostatics,
    nonlinear::{
        IsotropicMagneticMaterial, NonlinearMagneticControls, NonlinearMagneticError,
        NonlinearMagneticReport, PlanarNonlinearField, PlanarNonlinearMagnetostatics,
        PositiveQuadratic,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    field_storage::{
        admit_fixed_magnetic_refresh, admit_magnetic_temperature_refresh,
        linear_magnetic_permeability,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TriComplex2;
use std::{ops::Deref, sync::Arc};

/// Saturating and fixed-remanence induction with prescribed Eulerian motion.
pub mod moving;

/// Physical accounting allowances and the shared bounded nonlinear solve.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearEddyControls {
    pub nonlinear: NonlinearMagneticControls,
    pub energy_tolerance_j: f64,
    pub current_balance_tolerance_a: f64,
}
#[derive(Debug)]
pub enum NonlinearEddyError {
    Induction(EddyError),
    Magnetic(NonlinearMagneticError),
    Cancelled,
}
impl From<EddyError> for NonlinearEddyError {
    fn from(error: EddyError) -> Self {
        match error {
            EddyError::Cancelled => Self::Cancelled,
            other => Self::Induction(other),
        }
    }
}
impl From<MagneticError> for NonlinearEddyError {
    fn from(error: MagneticError) -> Self {
        NonlinearMagneticError::from(error).into()
    }
}
impl From<NonlinearMagneticError> for NonlinearEddyError {
    fn from(error: NonlinearMagneticError) -> Self {
        match error {
            NonlinearMagneticError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Magnetic(other),
        }
    }
}
impl core::fmt::Display for NonlinearEddyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear planar induction: {self:?}")
    }
}
impl std::error::Error for NonlinearEddyError {}

impl EddySection for PlanarNonlinearMagnetostatics {
    fn mesh(&self) -> &TriComplex2 {
        self.mesh()
    }
    fn thickness_m(&self) -> f64 {
        self.thickness_m()
    }
    fn admit_conductor(
        &self,
        triangle: usize,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<(), EddyError> {
        match &self.materials()[triangle] {
            IsotropicMagneticMaterial::Linear(original) => {
                admit_fixed_magnetic_refresh(original, material)?;
            }
            IsotropicMagneticMaterial::Reversible(law) => {
                law.admit_fixed_material_refresh(material)?;
            }
            IsotropicMagneticMaterial::Recoil(law) => {
                law.material().admit_fixed_material_refresh(material)?;
            }
        }
        Ok(())
    }
}

struct BackwardConduction<'a> {
    accepted: &'a [f64],
    conduction: &'a Conduction,
    dt: f64,
}
impl PositiveQuadratic for BackwardConduction<'_> {
    fn reference(&self) -> &[f64] {
        self.accepted
    }
    fn apply(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(), MagneticError> {
        if !self.conduction.apply(x, y, cancelled) {
            return Err(MagneticError::Cancelled);
        }
        for value in y {
            *value /= self.dt;
            if !value.is_finite() {
                return Err(MagneticError::Invalid(
                    "induction tangent is unrepresentable",
                ));
            }
        }
        Ok(())
    }
}

#[derive(Clone)]
enum MagneticOwner<'m> {
    Borrowed(&'m PlanarNonlinearMagnetostatics),
    Owned(Arc<PlanarNonlinearMagnetostatics>),
}
impl Deref for MagneticOwner<'_> {
    type Target = PlanarNonlinearMagnetostatics;
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Borrowed(problem) => problem,
            Self::Owned(problem) => problem,
        }
    }
}

/// Preserve every original source and law kind. An explicit separable nonlinear
/// thermal law may refresh alongside linear mu(T); its parameter work remains
/// in the native receipt and needs an external or finite thermal energy owner.
pub(crate) fn admit_material_refresh(
    before: &[IsotropicMagneticMaterial],
    after: &[IsotropicMagneticMaterial],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), NonlinearEddyError> {
    if before.len() != after.len() {
        return Err(
            EddyError::Invalid("complete original magnetic correspondence required").into(),
        );
    }
    for (before, after) in before.iter().zip(after) {
        poll(cancelled)?;
        match (before, after) {
            (IsotropicMagneticMaterial::Linear(a), IsotropicMagneticMaterial::Linear(b)) => {
                admit_magnetic_temperature_refresh(a, b).map_err(EddyError::from)?;
            }
            (
                IsotropicMagneticMaterial::Reversible(a),
                IsotropicMagneticMaterial::Reversible(b),
            ) => { a.admit_temperature_refresh(b).map_err(EddyError::from)?; }
            (IsotropicMagneticMaterial::Recoil(a), IsotropicMagneticMaterial::Recoil(b))
                if a == b => {}
            _ => return Err(EddyError::Invalid(
                "temperature refresh preserves linear/separable sources and complete fixed recoil laws",
            )
            .into()),
        }
    }
    Ok(())
}

/// Immutable accepted A/current/time history, retaining the original sources.
#[derive(Clone)]
pub struct PlanarNonlinearEddy<'m> {
    problem: MagneticOwner<'m>,
    conductors: Vec<EddyConductor>,
    insulators: Vec<usize>,
    boundary: Vec<(usize, f64)>,
    components: Vec<usize>,
    algebraic: Vec<Vec<usize>>,
    impressed_jz_a_m2: Vec<f64>,
    potential_wb_m: Vec<f64>,
    energy_j: f64,
    time_s: f64,
    report: NonlinearMagneticReport,
}
#[derive(Clone, Debug)]
pub struct NonlinearEddyMagneticCell {
    pub triangle: usize,
    pub volume_m3: f64,
    pub flux_density_t: [[f64; 2]; 2],
    /// Endpoint source-integrated energies, not secant H*B/2.
    pub free_energy_j: [f64; 2],
    /// Convex Bregman loss of backward Euler; not physical heat.
    pub numerical_dissipation_j: f64,
    /// Signed psi_next(B_old)-psi_old(B_old) at fixed original flux density.
    /// Neither Joule heat nor external impressed-current work.
    pub parameter_work_j: f64,
}
#[derive(Clone, Debug)]
pub struct NonlinearEddyStep {
    pub cells: Vec<EddyCellHeat>,
    pub magnetic_cells: Vec<NonlinearEddyMagneticCell>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Actual integrated sigma E^2; the only deposited magnetic heating here.
    pub heat_j: f64,
    /// Endpoint impressed source paired with the actual accepted delta A.
    pub supplied_work_j: f64,
    /// Nonnegative first-order time-discretization loss; never material heat.
    pub numerical_dissipation_j: f64,
    /// Sum of source-resolved material work; a thermal/external owner must close it.
    pub parameter_work_j: f64,
    /// Delta stored + heat + numerical dissipation - supplied work - parameter work.
    pub energy_residual_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
    /// Recomputed full nonlinear induction residual, including C_sigma/dt.
    pub report: NonlinearMagneticReport,
}
pub struct PlanarNonlinearEddyCandidate<'m> {
    state: PlanarNonlinearEddy<'m>,
    step: NonlinearEddyStep,
}
impl<'m> PlanarNonlinearEddyCandidate<'m> {
    pub fn state(&self) -> &PlanarNonlinearEddy<'m> {
        &self.state
    }
    pub fn step(&self) -> &NonlinearEddyStep {
        &self.step
    }
    /// Publish the accepted field together with its physical heat receiver.
    pub fn into_parts(self) -> (PlanarNonlinearEddy<'m>, NonlinearEddyStep) {
        (self.state, self.step)
    }
}

fn admit_controls(value: NonlinearEddyControls) -> Result<(), EddyError> {
    controls(EddyControls {
        linear: value.nonlinear.linear,
        energy_tolerance_j: value.energy_tolerance_j,
        current_balance_tolerance_a: value.current_balance_tolerance_a,
    })?;
    if !positive(value.nonlinear.absolute_residual_tolerance_a_m)
        || !value.nonlinear.relative_residual_tolerance.is_finite()
        || !(0.0..1.0).contains(&value.nonlinear.relative_residual_tolerance)
    {
        return Err(EddyError::Invalid(
            "nonlinear induction tolerances are invalid",
        ));
    }
    Ok(())
}

fn numerical_loss_density(
    material: &IsotropicMagneticMaterial,
    before: [f64; 2],
    after: [f64; 2],
) -> Result<f64, NonlinearEddyError> {
    let mu = match material {
        IsotropicMagneticMaterial::Linear(material) => {
            linear_magnetic_permeability(material).map_err(EddyError::from)?
        }
        IsotropicMagneticMaterial::Recoil(law) => law.material().permeability_h_m(),
        IsotropicMagneticMaterial::Reversible(law) => {
            return Ok(law
                .bregman_density_j_m3(before, after)
                .map_err(EddyError::from)?);
        }
    };
    // Fixed Br cancels from the convex increment. Do not subtract large bias
    // energies, or turn this time-discretization loss into constitutive heat.
    let difference = (after[0] - before[0]).hypot(after[1] - before[1]);
    Ok(finite((0.5 * difference) * (difference / mu))?)
}

impl<'m> PlanarNonlinearEddy<'m> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &PlanarNonlinearField<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        time_s: f64,
        limits: NonlinearEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, NonlinearEddyError> {
        Self::new_with_impressed_current(
            initial,
            conductors,
            insulating_triangles,
            boundary,
            &vec![0.0; initial.problem().mesh().faces().len()],
            time_s,
            limits,
            cancelled,
        )
    }

    /// Begin from the actual source-bound field. Insulating/floating algebraic
    /// constraints must already hold; no initialization impulse resets A.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_impressed_current(
        initial: &PlanarNonlinearField<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        impressed_jz_a_m2: &[f64],
        time_s: f64,
        limits: NonlinearEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, NonlinearEddyError> {
        poll(cancelled)?;
        admit_controls(limits)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(EddyError::Invalid("finite nonnegative initial time required").into());
        }
        let problem = initial.problem();
        let prescribed =
            PlanarMagnetostatics::admit_mesh_boundary(problem.mesh(), boundary, cancelled)?;
        if prescribed
            .iter()
            .zip(initial.potential_wb_m())
            .any(|(p, a)| p.is_some_and(|p| p.to_bits() != a.to_bits()))
        {
            return Err(
                EddyError::Invalid("fixed induction boundary must match accepted field").into(),
            );
        }
        let (conductors, _) = bind(problem, conductors, None, cancelled)?;
        let components = regions(problem, &conductors, insulating_triangles, cancelled)?;
        let algebraic = algebraic_blocks(problem, &conductors, &prescribed, cancelled)?;
        let loads = load(
            problem,
            impressed_jz_a_m2,
            &components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let magnetic = problem.magnetic_load(initial.potential_wb_m(), cancelled)?;
        admit_algebraic(
            problem,
            &magnetic,
            &loads,
            &algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            problem: MagneticOwner::Borrowed(problem),
            conductors,
            insulators: insulating_triangles.to_vec(),
            boundary: boundary.to_vec(),
            components,
            algebraic,
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            potential_wb_m: initial.potential_wb_m().to_vec(),
            energy_j: initial.energy_j(),
            time_s,
            report: initial.report().clone(),
        })
    }
    pub fn problem(&self) -> &PlanarNonlinearMagnetostatics {
        &self.problem
    }
    pub fn conductors(&self) -> &[EddyConductor] {
        &self.conductors
    }
    pub fn insulating_triangles(&self) -> &[usize] {
        &self.insulators
    }
    pub fn potential_wb_m(&self) -> &[f64] {
        &self.potential_wb_m
    }
    pub fn impressed_current_density_a_m2(&self) -> &[f64] {
        &self.impressed_jz_a_m2
    }
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Recover actual B/H/energy without re-equilibrating the evolving field.
    /// The returned field uses the existing explicit frame/position transfer.
    pub fn field(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearField<'_>, NonlinearEddyError> {
        Ok(self
            .problem
            .recover(self.potential_wb_m.clone(), self.report.clone(), cancelled)?)
    }

    /// One first-order backward-Euler candidate, retaining the accepted field.
    /// Endpoint current needs a balanced return on every mesh component.
    /// Conductivity may refresh within its pinned source/context; the magnetic
    /// law must remain unchanged over its full admitted field range. Failed
    /// Newton, accounting, source support or cancellation publishes no state.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        conductors: Vec<EddyConductor>,
        impressed_jz_a_m2: &[f64],
        limits: NonlinearEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearEddyCandidate<'m>, NonlinearEddyError> {
        self.step_with_owner(
            target_s,
            self.problem.clone(),
            conductors,
            impressed_jz_a_m2,
            limits,
            cancelled,
        )
    }

    /// Reassemble the original geometry with source-pinned linear mu(T) or an
    /// explicitly declared separable nonlinear H(B,T), retaining accepted A and
    /// the original field-energy reference. Endpoints include signed material
    /// work; this operation does not identify that work with physical heating.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        target_s: f64,
        materials: Vec<IsotropicMagneticMaterial>,
        conductors: Vec<EddyConductor>,
        impressed_jz_a_m2: &[f64],
        assembly: WeightedAssemblyLimits,
        limits: NonlinearEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearEddyCandidate<'m>, NonlinearEddyError> {
        poll(cancelled)?;
        admit_material_refresh(self.problem.materials(), &materials, cancelled)?;
        let problem = PlanarNonlinearMagnetostatics::new(
            self.problem.mesh().clone(),
            materials,
            assembly,
            cancelled,
        )?;
        self.step_with_owner(
            target_s,
            MagneticOwner::Owned(Arc::new(problem)),
            conductors,
            impressed_jz_a_m2,
            limits,
            cancelled,
        )
    }

    fn step_with_owner(
        &self,
        target_s: f64,
        next_problem: MagneticOwner<'m>,
        conductors: Vec<EddyConductor>,
        impressed_jz_a_m2: &[f64],
        limits: NonlinearEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearEddyCandidate<'m>, NonlinearEddyError> {
        poll(cancelled)?;
        admit_controls(limits)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite() || !positive(dt) {
            return Err(
                EddyError::Invalid("positive representable induction interval required").into(),
            );
        }
        let (conductors, conduction) = bind(
            &*next_problem,
            conductors,
            Some(&self.conductors),
            cancelled,
        )?;
        let loads = load(
            &*next_problem,
            impressed_jz_a_m2,
            &self.components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let quadratic = BackwardConduction {
            accepted: &self.potential_wb_m,
            conduction: &conduction,
            dt,
        };
        let field = next_problem.solve_loaded(
            &loads,
            &self.boundary,
            &self.potential_wb_m,
            Some(&quadratic),
            limits.nonlinear,
            cancelled,
        )?;
        let potential = field.potential_wb_m();
        let magnetic = next_problem.magnetic_load(potential, cancelled)?;
        admit_algebraic(
            &*next_problem,
            &magnetic,
            &loads,
            &self.algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let increments = potential
            .iter()
            .zip(&self.potential_wb_m)
            .map(|(a, b)| finite(a - b))
            .collect::<Result<Vec<_>, _>>()?;
        let EddyHeat {
            cells,
            heat_j,
            region_currents_a,
        } = heat_cells(
            &conduction,
            &conductors,
            &increments,
            dt,
            next_problem.thickness_m(),
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let old_field = self.field(cancelled)?;
        let refreshed_old =
            next_problem.recover(self.potential_wb_m.clone(), self.report.clone(), cancelled)?;
        let mut magnetic_cells = Vec::with_capacity(next_problem.materials().len());
        let mut numerical_dissipation_j = 0.0;
        let mut parameter_work_j = 0.0;
        for (triangle, material) in next_problem.materials().iter().enumerate() {
            poll(cancelled)?;
            let b0 = old_field.magnetic_flux_density_t()[triangle];
            let b1 = field.magnetic_flux_density_t()[triangle];
            let volume_m3 = self
                .problem
                .mesh()
                .face_measure(triangle)
                .ok_or(EddyError::Invalid("missing magnetic volume"))?;
            let density = numerical_loss_density(material, b0, b1)?;
            let numerical = finite(volume_m3 * density)?;
            if numerical < 0.0 {
                return Err(EddyError::Balance(
                    "convex magnetic numerical loss must be nonnegative",
                )
                .into());
            }
            numerical_dissipation_j = finite(numerical_dissipation_j + numerical)?;
            let parameter = finite(
                refreshed_old.cell_energy_j()[triangle] - old_field.cell_energy_j()[triangle],
            )?;
            parameter_work_j = finite(parameter_work_j + parameter)?;
            magnetic_cells.push(NonlinearEddyMagneticCell {
                triangle,
                volume_m3,
                flux_density_t: [b0, b1],
                free_energy_j: [
                    old_field.cell_energy_j()[triangle],
                    field.cell_energy_j()[triangle],
                ],
                numerical_dissipation_j: numerical,
                parameter_work_j: parameter,
            });
        }
        let supplied_work_j = finite(loads.iter().zip(&increments).map(|(f, d)| f * d).sum())?;
        let energy_residual_j = finite(
            field.energy_j() - self.energy_j + heat_j + numerical_dissipation_j
                - supplied_work_j
                - parameter_work_j,
        )?;
        if energy_residual_j.abs() > limits.energy_tolerance_j {
            return Err(EddyError::Balance("nonlinear magnetic storage, physical heat, numerical loss and source work do not close").into());
        }
        let step = NonlinearEddyStep {
            cells,
            magnetic_cells,
            stored_before_j: self.energy_j,
            stored_after_j: field.energy_j(),
            heat_j,
            supplied_work_j,
            numerical_dissipation_j,
            parameter_work_j,
            energy_residual_j,
            region_currents_a,
            report: field.report().clone(),
        };
        let potential_wb_m = potential.to_vec();
        let energy_j = field.energy_j();
        let report = field.report().clone();
        let state = Self {
            problem: next_problem,
            conductors,
            insulators: self.insulators.clone(),
            boundary: self.boundary.clone(),
            components: self.components.clone(),
            algebraic: self.algebraic.clone(),
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            potential_wb_m,
            energy_j,
            time_s: target_s,
            report,
        };
        poll(cancelled)?;
        Ok(PlanarNonlinearEddyCandidate { state, step })
    }
}
