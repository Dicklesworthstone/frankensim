//! Prescribed accelerating supports on the complete original elastic pencil.
//!
//! Support displacement, velocity and acceleration are retained continuously.
//! Each interval linearly interpolates support acceleration and applied force.
//! The free-coordinate source is therefore cubic, including the ORIGINAL
//! consistent cross mass: f_eff = f_f - K_fc u_c - M_fc a_c. Every free elastic
//! mode is retained by the shared complete-spectrum owner.
//!
//! Reactions use M a + K u - f in the original coordinates. Their work is
//! independently integrated from the actual support trajectory and displacement
//! moments; it is never filled from the full stored-energy difference. Both
//! subsystem and full-system energy gates precede accepted-state publication.
//! Omitted assembly coordinates remain homogeneous stationary supports; their
//! reactions cannot be recovered from this already reduced input. There is no
//! viscous law, changing topology or finite rotation. Instantaneous support
//! velocity changes and applied impulses use the same original mass partition.

use super::LinearSolidImpulseError;
use super::spectral::{
    CompleteSpectralSolid, CompleteSpectralSolidControls, CompleteSpectralSolidError,
};
use crate::linear3::TetElasticAssembly;
use fs_exec::Cx;
use fs_sparse::{Coo, Csr};

mod admission;
mod frictional_impact;
mod impact;
mod impulse;
mod sustained;
pub use frictional_impact::{
    SupportedSolidFrictionalImpactError, SupportedSolidFrictionalImpactReport,
};
pub use impact::{SupportedSolidNormalImpactError, SupportedSolidNormalImpactReport};
pub(super) use impulse::SupportedImpulseCandidate;
pub use impulse::SupportedSolidImpulse;

/// Whole-system limits in addition to the free spectral owner's limits.
#[derive(Clone, Copy, Debug)]
pub struct SupportedSolidControls {
    pub spectral: CompleteSpectralSolidControls,
    pub maximum_total_dofs: usize,
    /// Conservative simultaneous original/partition/spectral scratch envelope.
    pub maximum_preparation_bytes: usize,
    /// Whole preparation admits 12*n^3 + 128*n^2 + 16384*n product terms.
    pub maximum_preparation_multiply_terms: u64,
    /// Whole interval admits 64*n^2 + 8192*n product terms.
    pub maximum_step_multiply_terms: u64,
    /// Absolute original free-row force-balance allowance [N].
    pub maximum_free_force_residual_n: f64,
    /// Absolute interval and cumulative full-system balance allowance [J].
    pub maximum_work_error_j: f64,
}

#[derive(Debug)]
pub enum SupportedSolidError {
    Invalid(&'static str),
    Budget { required: u128, maximum: u128 },
    Operator { defect: f64, allowance: f64 },
    Energy { defect_j: f64, allowance_j: f64 },
    FreeForce { residual_n: f64, allowance_n: f64 },
    Spectral(CompleteSpectralSolidError),
    Impulse(LinearSolidImpulseError),
    Cancelled,
}
impl core::fmt::Display for SupportedSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "supported spectral solid refused: {self:?}")
    }
}
impl core::error::Error for SupportedSolidError {}
impl From<CompleteSpectralSolidError> for SupportedSolidError {
    fn from(error: CompleteSpectralSolidError) -> Self {
        match error {
            CompleteSpectralSolidError::Cancelled => Self::Cancelled,
            other => Self::Spectral(other),
        }
    }
}
impl From<LinearSolidImpulseError> for SupportedSolidError {
    fn from(error: LinearSolidImpulseError) -> Self {
        match error {
            LinearSolidImpulseError::Cancelled => Self::Cancelled,
            other => Self::Impulse(other),
        }
    }
}
type Result<T> = core::result::Result<T, SupportedSolidError>;
fn invalid(message: &'static str) -> SupportedSolidError {
    SupportedSolidError::Invalid(message)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| SupportedSolidError::Cancelled)
}
fn zeros(n: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| invalid("allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}
fn copy<T: Copy>(values: &[T]) -> Result<Vec<T>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(values.len())
        .map_err(|_| invalid("allocation refused"))?;
    result.extend_from_slice(values);
    Ok(result)
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}
fn bound(required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(SupportedSolidError::Budget { required, maximum })
    } else {
        Ok(())
    }
}

/// Full original-coordinate energy and independently integrated physical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportedSolidEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub applied_work_j: f64,
    /// Work supplied BY the prescribed supports TO the solid [J].
    pub support_work_j: f64,
    pub cumulative_applied_work_j: f64,
    pub cumulative_support_work_j: f64,
    /// W_applied + W_support - Delta H [J].
    pub balance_error_j: f64,
    pub cumulative_balance_error_j: f64,
}

#[derive(Clone, Debug)]
pub struct SupportedSolidStep {
    pub step: u64,
    pub time_s: f64,
    /// Delta u / dt in the original represented-coordinate order [m/s].
    pub mean_velocity_mps: Vec<f64>,
    /// Endpoint force exerted BY each support ON the solid [N], support order.
    pub support_reaction_n: Vec<f64>,
    pub free_force_residual_n: f64,
    pub energy: SupportedSolidEnergy,
}

/// Extract the complete free pencil used by [`SupportedSpectralSolid`].
/// `support_dofs` are sorted global labels present in assembly.free_dofs.
/// The returned free map is the original map with those labels removed.
/// This cold projection admits dimensions, finite symmetric entries and the
/// explicit whole-system capacities. Full K/M positivity is admitted when
/// preparing motion. Metadata retains original physical reference mass/volume.
pub fn partition_free_pencil(
    assembly: &TetElasticAssembly,
    support_dofs: &[usize],
    controls: SupportedSolidControls,
    cx: &Cx<'_>,
) -> Result<TetElasticAssembly> {
    let (free, _) = partition_rows(assembly, support_dofs, controls, cx)?;
    let n = assembly.free_dofs.len();
    let stiffness = symmetric(
        &assembly.stiffness,
        n,
        controls.spectral.maximum_relative_operator_error,
        cx,
    )?;
    let mass = symmetric(
        &assembly.mass,
        n,
        controls.spectral.maximum_relative_operator_error,
        cx,
    )?;
    free_pencil(assembly, &free, &stiffness, &mass, cx)
}

/// Exact conservative motion with explicitly driven support accelerations.
pub struct SupportedSpectralSolid {
    free: CompleteSpectralSolid,
    coordinates: Vec<usize>,
    support_dofs: Vec<usize>,
    free_rows: Vec<usize>,
    support_rows: Vec<usize>,
    stiffness: Vec<f64>,
    mass: Vec<f64>,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    support_acceleration: Vec<f64>,
    dt_s: f64,
    controls: SupportedSolidControls,
    stored_j: f64,
    kinetic_j: f64,
    initial_energy_j: f64,
    applied_work_j: f64,
    applied_correction_j: f64,
    support_work_j: f64,
    support_correction_j: f64,
}

// Complete unpublished original interval, including its free modal state and
// independently integrated full-coordinate support reactions/work.
struct SupportedIntervalTrial {
    free: super::spectral::polynomial::PolynomialTrial,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    support_acceleration: Vec<f64>,
    applied_correction_j: f64,
    support_correction_j: f64,
    output: SupportedSolidStep,
}

impl SupportedSpectralSolid {
    /// Supply ALL generalized eigenpairs of [`partition_free_pencil`]. Initial
    /// u and v and applied forces use the ORIGINAL assembly.free_dofs order;
    /// support accelerations use support_dofs order. Full original M must be
    /// SPD and K PSD to the explicit numerical operator allowance. Positivity
    /// checks are row/tile-polled and do not modify retained matrix entries.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_basis(
        assembly: &TetElasticAssembly,
        support_dofs: &[usize],
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        initial_support_acceleration_mps2: &[f64],
        controls: SupportedSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        let (free_rows, support_rows) = partition_rows(assembly, support_dofs, controls, cx)?;
        let n = assembly.free_dofs.len();
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || initial_displacement_m.len() != n
            || initial_velocity_mps.len() != n
            || initial_support_acceleration_mps2.len() != support_dofs.len()
            || initial_displacement_m
                .iter()
                .chain(initial_velocity_mps)
                .chain(initial_support_acceleration_mps2)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite initial motion and positive clock required",
            ));
        }
        let tolerance = controls.spectral.maximum_relative_operator_error;
        let stiffness = symmetric(&assembly.stiffness, n, tolerance, cx)?;
        let mass = symmetric(&assembly.mass, n, tolerance, cx)?;
        admission::admit_operators(&stiffness, &mass, n, tolerance, cx)?;
        let pencil = free_pencil(assembly, &free_rows, &stiffness, &mass, cx)?;
        let mut free_u = zeros(free_rows.len())?;
        let mut free_v = zeros(free_rows.len())?;
        for (i, &row) in free_rows.iter().enumerate() {
            free_u[i] = initial_displacement_m[row];
            free_v[i] = initial_velocity_mps[row];
        }
        let free = CompleteSpectralSolid::prepare_cubic_forces_with_basis(
            &pencil,
            eigenvalues_rad2_s2,
            basis_columns,
            dt_s,
            &free_u,
            &free_v,
            controls.spectral,
            cx,
        )?;
        let stored_j = quadratic(&stiffness, initial_displacement_m, cx)?;
        let kinetic_j = quadratic(&mass, initial_velocity_mps, cx)?;
        if !(stored_j + kinetic_j).is_finite()
            || stored_j < -controls.maximum_work_error_j
            || kinetic_j < 0.0
        {
            return Err(invalid("invalid original physical energy"));
        }
        let result = Self {
            free,
            coordinates: copy(&assembly.free_dofs)?,
            support_dofs: copy(support_dofs)?,
            free_rows,
            support_rows,
            stiffness,
            mass,
            displacement: copy(initial_displacement_m)?,
            velocity: copy(initial_velocity_mps)?,
            support_acceleration: copy(initial_support_acceleration_mps2)?,
            dt_s,
            controls,
            stored_j,
            kinetic_j,
            initial_energy_j: stored_j + kinetic_j,
            applied_work_j: 0.0,
            applied_correction_j: 0.0,
            support_work_j: 0.0,
            support_correction_j: 0.0,
        };
        poll(cx)?;
        Ok(result)
    }

    pub fn represented_dofs(&self) -> &[usize] {
        &self.coordinates
    }
    pub fn free_dofs(&self) -> &[usize] {
        self.free.free_dofs()
    }
    pub fn support_dofs(&self) -> &[usize] {
        &self.support_dofs
    }
    pub fn displacement_m(&self) -> &[f64] {
        &self.displacement
    }
    pub fn velocity_mps(&self) -> &[f64] {
        &self.velocity
    }
    pub fn support_acceleration_mps2(&self) -> &[f64] {
        &self.support_acceleration
    }
    pub fn time_s(&self) -> f64 {
        self.free.time_s()
    }
    pub fn accepted_steps(&self) -> u64 {
        self.free.steps()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_j
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_j
    }
    pub fn cumulative_applied_work_j(&self) -> f64 {
        self.applied_work_j
    }
    pub fn cumulative_support_work_j(&self) -> f64 {
        self.support_work_j
    }

    /// Interpolate retained support acceleration to the supplied endpoint,
    /// integrating it to continuous support velocity and displacement. Forces
    /// are affine over the interval; finite force jumps apply no impulse.
    /// All original-coordinate work/reaction gates precede either owner commit.
    pub fn advance_affine_supports(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        force_start_n: &[f64],
        force_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SupportedSolidStep> {
        let trial = self.trial_affine_supports(
            support_acceleration_end_mps2,
            force_start_n,
            force_end_n,
            cx,
        )?;
        Ok(self.accept_affine_supports(trial))
    }

    fn trial_affine_supports(
        &self,
        support_acceleration_end_mps2: &[f64],
        force_start_n: &[f64],
        force_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SupportedIntervalTrial> {
        poll(cx)?;
        let n = self.coordinates.len();
        let count = self.support_rows.len();
        if support_acceleration_end_mps2.len() != count
            || force_start_n.len() != n
            || force_end_n.len() != n
            || support_acceleration_end_mps2
                .iter()
                .chain(force_start_n)
                .chain(force_end_n)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite support acceleration and force endpoints required",
            ));
        }
        let h = self.dt_s;
        let mut support_coefficients = [zeros(count)?, zeros(count)?, zeros(count)?, zeros(count)?];
        let mut acceleration_increment = zeros(count)?;
        let mut force_increment = zeros(n)?;
        let mut displacement = copy(&self.displacement)?;
        let mut velocity = copy(&self.velocity)?;
        let mut acceleration = zeros(n)?;
        let mut moments = [zeros(n)?, zeros(n)?, zeros(n)?];
        for i in 0..n {
            force_increment[i] = force_end_n[i] - force_start_n[i];
            if !force_increment[i].is_finite() {
                return Err(invalid("force increment overflow"));
            }
        }
        for (c, &row) in self.support_rows.iter().enumerate() {
            poll(cx)?;
            let da = support_acceleration_end_mps2[c] - self.support_acceleration[c];
            acceleration_increment[c] = da;
            support_coefficients[0][c] = self.displacement[row];
            support_coefficients[1][c] = h * self.velocity[row];
            support_coefficients[2][c] = 0.5 * h * (h * self.support_acceleration[c]);
            support_coefficients[3][c] = (h / 6.0) * (h * da);
            displacement[row] = (0..4).map(|r| support_coefficients[r][c]).sum();
            velocity[row] = self.velocity[row] + h * (self.support_acceleration[c] + 0.5 * da);
            acceleration[row] = support_acceleration_end_mps2[c];
            for (k, moment) in moments.iter_mut().enumerate() {
                moment[row] = (0..4)
                    .map(|r| support_coefficients[r][c] / (k + r + 1) as f64)
                    .sum();
            }
            if !da.is_finite()
                || !displacement[row].is_finite()
                || !velocity[row].is_finite()
                || support_coefficients
                    .iter()
                    .any(|values| !values[c].is_finite())
                || moments.iter().any(|values| !values[row].is_finite())
            {
                return Err(invalid("unrepresentable prescribed support trajectory"));
            }
        }
        let nf = self.free_rows.len();
        let mut force = [zeros(nf)?, zeros(nf)?, zeros(nf)?, zeros(nf)?];
        for (f, &row) in self.free_rows.iter().enumerate() {
            poll(cx)?;
            force[0][f] = force_start_n[row];
            force[1][f] = force_increment[row];
            for (c, &column) in self.support_rows.iter().enumerate() {
                for r in 0..4 {
                    force[r][f] = (-self.stiffness[row * n + column])
                        .mul_add(support_coefficients[r][c], force[r][f]);
                }
                force[0][f] = (-self.mass[row * n + column])
                    .mul_add(self.support_acceleration[c], force[0][f]);
                force[1][f] =
                    (-self.mass[row * n + column]).mul_add(acceleration_increment[c], force[1][f]);
            }
        }
        let trial = self
            .free
            .trial_polynomial([&force[0], &force[1], &force[2], &force[3]], cx)?;
        for (f, &row) in self.free_rows.iter().enumerate() {
            displacement[row] = trial.displacement_m()[f];
            velocity[row] = trial.velocity_mps()[f];
            acceleration[row] = trial.acceleration_mps2()[f];
            for k in 0..3 {
                moments[k][row] = trial.displacement_moments_m()[k][f];
            }
        }
        let mut residual = zeros(n)?;
        for i in 0..n {
            poll(cx)?;
            residual[i] = (0..n).fold(-force_end_n[i], |sum, j| {
                self.mass[i * n + j].mul_add(
                    acceleration[j],
                    self.stiffness[i * n + j].mul_add(displacement[j], sum),
                )
            });
            if !residual[i].is_finite() {
                return Err(invalid("reaction overflow"));
            }
        }
        let free_force_residual_n = self
            .free_rows
            .iter()
            .fold(0.0_f64, |m, &i| m.max(residual[i].abs()));
        if free_force_residual_n > self.controls.maximum_free_force_residual_n {
            return Err(SupportedSolidError::FreeForce {
                residual_n: free_force_residual_n,
                allowance_n: self.controls.maximum_free_force_residual_n,
            });
        }
        let mut applied_by_row = zeros(n)?;
        let (mut applied_work_j, mut applied_roundoff) = (0.0, 0.0);
        let mut mean_velocity_mps = zeros(n)?;
        for i in 0..n {
            poll(cx)?;
            let du = displacement[i] - self.displacement[i];
            applied_by_row[i] = force_start_n[i]
                .mul_add(du, force_increment[i] * (displacement[i] - moments[0][i]));
            add(
                &mut applied_work_j,
                &mut applied_roundoff,
                applied_by_row[i],
            );
            mean_velocity_mps[i] = du / h;
        }
        let mut support_reaction_n = zeros(count)?;
        let (mut support_work_j, mut support_roundoff) = (0.0, 0.0);
        for (c, &row) in self.support_rows.iter().enumerate() {
            poll(cx)?;
            support_reaction_n[c] = residual[row];
            let v0 = self.velocity[row];
            let v1 = velocity[row];
            let a0 = self.support_acceleration[c];
            let a1 = support_acceleration_end_mps2[c];
            let da = acceleration_increment[c];
            for j in 0..n {
                // Twice integrating the ORIGINAL reaction's inertial term by
                // parts: integral vc*aj = [vc*vj]-[ac*uj]+Delta ac*I0j.
                // This rearrangement avoids subtracting large position terms.
                let mass_work = v1 * (velocity[j] - self.velocity[j])
                    + (v1 - v0) * self.velocity[j]
                    - a1 * (displacement[j] - self.displacement[j])
                    + da * (moments[0][j] - self.displacement[j]);
                let stiffness_work =
                    h * (v0 * moments[0][j] + h * (a0 * moments[1][j] + 0.5 * da * moments[2][j]));
                add(
                    &mut support_work_j,
                    &mut support_roundoff,
                    self.mass[row * n + j]
                        .mul_add(mass_work, self.stiffness[row * n + j] * stiffness_work),
                );
            }
            add(
                &mut support_work_j,
                &mut support_roundoff,
                -applied_by_row[row],
            );
        }
        let stored_after_j = quadratic(&self.stiffness, &displacement, cx)?;
        let kinetic_after_j = quadratic(&self.mass, &velocity, cx)?;
        if !(stored_after_j + kinetic_after_j).is_finite()
            || stored_after_j < -self.controls.maximum_work_error_j
            || kinetic_after_j < 0.0
        {
            return Err(invalid("invalid original endpoint physical energy"));
        }
        let change = quadratic_change(&self.stiffness, &self.displacement, &displacement, cx)?
            + quadratic_change(&self.mass, &self.velocity, &velocity, cx)?;
        let (mut cumulative_applied_work_j, mut applied_correction_j) =
            (self.applied_work_j, self.applied_correction_j);
        let (mut cumulative_support_work_j, mut support_correction_j) =
            (self.support_work_j, self.support_correction_j);
        add(
            &mut cumulative_applied_work_j,
            &mut applied_correction_j,
            applied_work_j,
        );
        add(
            &mut cumulative_support_work_j,
            &mut support_correction_j,
            support_work_j,
        );
        let energy = SupportedSolidEnergy {
            stored_before_j: self.stored_j,
            stored_after_j,
            kinetic_before_j: self.kinetic_j,
            kinetic_after_j,
            applied_work_j,
            support_work_j,
            cumulative_applied_work_j,
            cumulative_support_work_j,
            balance_error_j: applied_work_j + support_work_j - change,
            cumulative_balance_error_j: cumulative_applied_work_j + cumulative_support_work_j
                - (stored_after_j + kinetic_after_j - self.initial_energy_j),
        };
        if [
            applied_work_j,
            support_work_j,
            cumulative_applied_work_j,
            cumulative_support_work_j,
            applied_correction_j,
            support_correction_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || mean_velocity_mps.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("unrepresentable physical work or mean velocity"));
        }
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > self.controls.maximum_work_error_j {
                return Err(SupportedSolidError::Energy {
                    defect_j,
                    allowance_j: self.controls.maximum_work_error_j,
                });
            }
        }
        let step = trial.motion().step;
        let time_s = trial.motion().time_s;
        let accepted_acceleration = copy(support_acceleration_end_mps2)?;
        poll(cx)?;
        Ok(SupportedIntervalTrial {
            free: trial,
            displacement,
            velocity,
            support_acceleration: accepted_acceleration,
            applied_correction_j,
            support_correction_j,
            output: SupportedSolidStep {
                step,
                time_s,
                mean_velocity_mps,
                support_reaction_n,
                free_force_residual_n,
                energy,
            },
        })
    }

    fn accept_affine_supports(&mut self, trial: SupportedIntervalTrial) -> SupportedSolidStep {
        debug_assert_eq!(self.accepted_steps() + 1, trial.output.step);
        // All original free and full-system checks precede this shared,
        // infallible publication boundary. No prepared operator is cloned.
        self.free.commit_polynomial(trial.free);
        self.displacement = trial.displacement;
        self.velocity = trial.velocity;
        self.support_acceleration = trial.support_acceleration;
        self.stored_j = trial.output.energy.stored_after_j;
        self.kinetic_j = trial.output.energy.kinetic_after_j;
        self.applied_work_j = trial.output.energy.cumulative_applied_work_j;
        self.applied_correction_j = trial.applied_correction_j;
        self.support_work_j = trial.output.energy.cumulative_support_work_j;
        self.support_correction_j = trial.support_correction_j;
        trial.output
    }
}

fn partition_rows(
    assembly: &TetElasticAssembly,
    support_dofs: &[usize],
    controls: SupportedSolidControls,
    cx: &Cx<'_>,
) -> Result<(Vec<usize>, Vec<usize>)> {
    poll(cx)?;
    let n = assembly.free_dofs.len();
    n.checked_mul(n)
        .ok_or_else(|| invalid("matrix shape overflow"))?;
    if n < 2
        || support_dofs.is_empty()
        || support_dofs.len() >= n
        || !controls
            .spectral
            .maximum_relative_operator_error
            .is_finite()
        || !(0.0..1.0).contains(&controls.spectral.maximum_relative_operator_error)
        || !controls.maximum_free_force_residual_n.is_finite()
        || controls.maximum_free_force_residual_n < 0.0
        || !controls.maximum_work_error_j.is_finite()
        || controls.maximum_work_error_j < 0.0
    {
        return Err(invalid(
            "nonempty free/support partition and valid numerical controls required",
        ));
    }
    bound(n as u128, controls.maximum_total_dofs as u128)?;
    let nn = n as u128;
    bound(
        192 * nn * nn + 8192 * nn + 8 * assembly.element_volumes_m3.len() as u128,
        controls.maximum_preparation_bytes as u128,
    )?;
    bound(
        12 * nn * nn * nn + 128 * nn * nn + 16384 * nn,
        controls.maximum_preparation_multiply_terms as u128,
    )?;
    bound(
        64 * nn * nn + 8192 * nn,
        controls.maximum_step_multiply_terms as u128,
    )?;
    if assembly.free_dofs.windows(2).any(|w| w[0] >= w[1])
        || support_dofs.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(invalid(
            "coordinate and support labels must be strictly increasing",
        ));
    }
    let mut supports = Vec::new();
    supports
        .try_reserve_exact(support_dofs.len())
        .map_err(|_| invalid("support map allocation refused"))?;
    for &dof in support_dofs {
        poll(cx)?;
        supports.push(
            assembly
                .free_dofs
                .binary_search(&dof)
                .map_err(|_| invalid("support is absent from the original represented pencil"))?,
        );
    }
    let mut free = Vec::new();
    free.try_reserve_exact(n - supports.len())
        .map_err(|_| invalid("free map allocation refused"))?;
    let mut c = 0;
    for i in 0..n {
        poll(cx)?;
        if c < supports.len() && supports[c] == i {
            c += 1;
        } else {
            free.push(i);
        }
    }
    Ok((free, supports))
}

fn free_pencil(
    assembly: &TetElasticAssembly,
    free: &[usize],
    stiffness: &[f64],
    mass: &[f64],
    cx: &Cx<'_>,
) -> Result<TetElasticAssembly> {
    let (n, nf) = (assembly.free_dofs.len(), free.len());
    let (mut k, mut m) = (Coo::new(nf, nf), Coo::new(nf, nf));
    let mut free_dofs = Vec::new();
    free_dofs
        .try_reserve_exact(nf)
        .map_err(|_| invalid("free map allocation refused"))?;
    for (i, &row) in free.iter().enumerate() {
        poll(cx)?;
        free_dofs.push(assembly.free_dofs[row]);
        for (j, &column) in free.iter().enumerate() {
            if stiffness[row * n + column] != 0.0 {
                k.push(i, j, stiffness[row * n + column]);
            }
            if mass[row * n + column] != 0.0 {
                m.push(i, j, mass[row * n + column]);
            }
        }
    }
    let stiffness = k.assemble();
    poll(cx)?;
    let mass = m.assemble();
    poll(cx)?;
    Ok(TetElasticAssembly {
        stiffness,
        mass,
        free_dofs,
        total_mass_kg: assembly.total_mass_kg,
        minimum_scaled_jacobian: assembly.minimum_scaled_jacobian,
        element_volumes_m3: copy(&assembly.element_volumes_m3)?,
    })
}

fn symmetric(matrix: &Csr, n: usize, tolerance: f64, cx: &Cx<'_>) -> Result<Vec<f64>> {
    if matrix.nrows() != n || matrix.ncols() != n {
        return Err(invalid(
            "source matrix shape differs from original coordinate map",
        ));
    }
    let mut dense = zeros(n * n)?;
    for i in 0..n {
        poll(cx)?;
        let (columns, values) = matrix.row(i);
        for (&j, &value) in columns.iter().zip(values) {
            if !value.is_finite() {
                return Err(invalid("source matrix must be finite"));
            }
            dense[i * n + j] = value;
        }
    }
    let scale = dense.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    for i in 0..n {
        poll(cx)?;
        for j in 0..i {
            let defect = if scale == 0.0 {
                0.0
            } else {
                (dense[i * n + j] / scale - dense[j * n + i] / scale).abs()
            };
            if defect > tolerance {
                return Err(SupportedSolidError::Operator {
                    defect,
                    allowance: tolerance,
                });
            }
            let value = 0.5 * dense[i * n + j] + 0.5 * dense[j * n + i];
            dense[i * n + j] = value;
            dense[j * n + i] = value;
        }
    }
    Ok(dense)
}

fn quadratic(matrix: &[f64], x: &[f64], cx: &Cx<'_>) -> Result<f64> {
    let n = x.len();
    let (mut total, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(cx)?;
        let row = (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s));
        add(&mut total, &mut correction, 0.5 * x[i] * row);
    }
    if !total.is_finite() {
        return Err(invalid("physical energy overflow"));
    }
    Ok(total)
}
fn quadratic_change(matrix: &[f64], before: &[f64], after: &[f64], cx: &Cx<'_>) -> Result<f64> {
    let n = before.len();
    let (mut total, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(cx)?;
        let row = (0..n).fold(0.0, |s, j| {
            matrix[i * n + j].mul_add(after[j] - before[j], s)
        });
        add(
            &mut total,
            &mut correction,
            (0.5 * before[i] + 0.5 * after[i]) * row,
        );
    }
    if !total.is_finite() {
        return Err(invalid("physical energy increment overflow"));
    }
    Ok(total)
}
