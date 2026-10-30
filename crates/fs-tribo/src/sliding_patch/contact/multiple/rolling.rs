//! Simultaneous original Hertz surfaces and linear isotropic rolling loss.
//!
//! Per patch the conjugate coordinates are
//! `[du_x, du_y, a*dtheta_spin, R*dtheta_roll_x, R*dtheta_roll_y]` [m],
//! with acting efforts `[Fx, Fy, Mspin/a, Mroll_x/R, Mroll_y/R]` [N].
//! The mechanical owner supplies the complete positive-semidefinite response
//! `x = free + C*w`. The rolling owner supplies its original fixed displacement
//! coefficient `d >= 0` [N/m], so `w_roll = -d*x_roll`; for the admitted
//! viscoelastic sphere this is `d=A*N/(R*dt)`. This algebra does not infer a
//! material, normal load, radius, interval, or applicability domain.
//!
//! One symmetric, diagonally equilibrated factor of
//! `B=I+sqrt(D)*Crr*sqrt(D)` eliminates every rolling coordinate together.
//! The identity is the original viscous equation, not added body compliance.
//! The unchanged multi-contact Hertz owner receives the resulting full Schur
//! response, including dependent static load redistribution. No individual
//! contact advances time, changes a checkpoint, or thermalizes another's work.
//!
//! Numerical symmetry is admitted by the original full-response tolerance.
//! Final compatibility, rolling law and work use the unmodified supplied C.
//! A singular/unresolved numerical Schur response refuses: no negative pivot,
//! compliance, or dissipation is clamped into a favorable answer. This bounded
//! reduction inherits the original Hertz solver's nonuniqueness/convergence
//! limits and the mechanical/material owners' held-geometry approximation.

use core::fmt;

use super::super::HertzWrenchPhase;
use fs_math::det;

use crate::contact_response::{ContactResponseError, PositiveFactor};

use super::super::{HertzWrenchConfig, HertzWrenchInput, HertzWrenchMode, Problem};
use super::{
    HertzContactInput, HertzWrenchContact, MultipleHertzConfig, MultipleHertzError,
    MultipleHertzStep, PreparedHertzContacts,
};

/// Explicit full-response allocation and physical residual budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleHertzRollingConfig {
    /// Original Hertz work/effort/support controls and the reduced 3N matrix cap.
    pub hertz: MultipleHertzConfig,
    /// Maximum entries of the supplied 5N-by-5N matrix, from 25 to 409600.
    /// Preparation retains at most three full-size dense arrays at once;
    /// repeated coefficient trials use the retained factors and O(N) scratch.
    pub maximum_matrix_terms: usize,
    /// Resultant residual of `w_roll + d*x_roll` in scaled force units [N].
    pub rolling_force_tolerance_n: f64,
    /// Absolute total physical heat minus original full-response work [J].
    pub total_energy_tolerance_j: f64,
}

impl MultipleHertzRollingConfig {
    pub fn validate(self) -> Result<(), MultipleHertzRollingError> {
        self.hertz.validate()?;
        if !(25..=409600).contains(&self.maximum_matrix_terms)
            || !self.rolling_force_tolerance_n.is_finite()
            || self.rolling_force_tolerance_n < 0.0
            || !self.total_energy_tolerance_j.is_finite()
            || self.total_energy_tolerance_j < 0.0
        {
            return Err(invalid("joint Hertz/rolling budgets"));
        }
        Ok(())
    }
}

/// One original surface law and released five-coordinate motion [m].
#[derive(Clone, Copy, Debug)]
pub struct HertzRollingContactInput<'a> {
    pub law: &'a HertzWrenchContact,
    pub free_coordinates_m: [f64; 5],
}

/// One original linear rolling response. No history or material identity is made.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzRollingStep {
    pub scaled_displacement_m: [f64; 2],
    pub scaled_force_n: [f64; 2],
    pub dissipated_energy_j: f64,
    /// Positive absorbed work at the original full-response motion [J].
    pub contact_work_j: f64,
    pub work_residual_j: f64,
}

/// One indivisible constitutive candidate with separate surface/bulk losses.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipleHertzRollingStep {
    /// Original 3N Hertz solution on the eliminated response. Its released
    /// motion is a numerical Schur coordinate, never external applied work.
    pub hertz: MultipleHertzStep,
    pub rolling: Vec<HertzRollingStep>,
    /// Constitutive motion; compatibility with the original full C is audited.
    pub scaled_motion_m: Vec<[f64; 5]>,
    pub scaled_acting_wrench_n: Vec<[f64; 5]>,
    pub maximum_compatibility_residual_m: f64,
    pub maximum_rolling_force_residual_n: f64,
    pub maximum_surface_force_residual_n: f64,
    pub maximum_surface_spin_moment_residual_nm: f64,
    pub surface_heat_j: f64,
    pub bulk_heat_j: f64,
    pub surface_work_j: f64,
    pub rolling_work_j: f64,
    pub surface_work_residual_j: f64,
    pub rolling_work_residual_j: f64,
    /// Positive absorbed surface plus rolling work at original full-C motion.
    pub contact_work_j: f64,
    pub work_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MultipleHertzRollingError {
    Invalid(&'static str),
    Hertz(MultipleHertzError),
    Cancelled,
    Audit {
        maximum_compatibility_residual_m: f64,
        maximum_rolling_force_residual_n: f64,
        maximum_surface_force_residual_n: f64,
        maximum_surface_spin_moment_residual_nm: f64,
        maximum_patch_work_residual_j: f64,
        surface_work_residual_j: f64,
        rolling_work_residual_j: f64,
        total_work_residual_j: f64,
    },
}

impl From<MultipleHertzError> for MultipleHertzRollingError {
    fn from(value: MultipleHertzError) -> Self {
        match value {
            MultipleHertzError::Cancelled => Self::Cancelled,
            error => Self::Hertz(error),
        }
    }
}
impl From<ContactResponseError> for MultipleHertzRollingError {
    fn from(value: ContactResponseError) -> Self {
        match value {
            ContactResponseError::Invalid(field) => invalid(field),
            ContactResponseError::Cancelled => Self::Cancelled,
        }
    }
}
impl fmt::Display for MultipleHertzRollingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "simultaneous Hertz and linear rolling contacts: {self:?}"
        )
    }
}
impl std::error::Error for MultipleHertzRollingError {}

/// Immutable full response, original rolling coefficients and reusable factors.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedHertzRollingContacts {
    count: usize,
    compliance: Vec<f64>,
    damping: Vec<f64>,
    root_damping: Vec<f64>,
    rolling_factor: EquilibratedFactor,
    /// L^-1*S^-1*sqrt(D)*Crh, row-major 2N by 3N.
    transfer: Vec<f64>,
    hertz: PreparedHertzContacts,
    config: MultipleHertzRollingConfig,
}

impl PreparedHertzRollingContacts {
    pub fn new(
        count: usize,
        scaled_compliance_m_per_n: Vec<f64>,
        rolling_damping_n_per_m: Vec<f64>,
        config: MultipleHertzRollingConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MultipleHertzRollingError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if count == 0 || count > config.hertz.maximum_contacts {
            return Err(invalid("joint contact count"));
        }
        let n = 5 * count;
        let h = 3 * count;
        let r = 2 * count;
        if n * n > config.maximum_matrix_terms
            || scaled_compliance_m_per_n.len() != n * n
            || h * h > config.hertz.maximum_matrix_terms
            || rolling_damping_n_per_m.len() != count
            || rolling_damping_n_per_m
                .iter()
                .any(|d| !d.is_finite() || *d < 0.0)
        {
            return Err(invalid("joint matrix shape, damping or allocation budget"));
        }
        // Admission only. Dropping this factor before Schur preparation bounds
        // peak memory and prevents a numerical range projection replacing C.
        drop(PositiveFactor::new(
            &scaled_compliance_m_per_n,
            n,
            config.hertz.local.relative_matrix_tolerance,
            &mut cancelled,
        )?);
        let root_damping: Vec<_> = rolling_damping_n_per_m
            .iter()
            .map(|&d| det::sqrt(d))
            .collect();
        let symmetric = |i: usize, j: usize| {
            let a = scaled_compliance_m_per_n[i * n + j];
            let b = scaled_compliance_m_per_n[j * n + i];
            if a == b { a } else { 0.5 * a + 0.5 * b }
        };
        let mut b = vec![0.0; r * r];
        for i in 0..r {
            poll(&mut cancelled)?;
            for j in 0..=i {
                let value = number(
                    root_damping[i / 2]
                        * symmetric(rolling_index(i), rolling_index(j))
                        * root_damping[j / 2]
                        + if i == j { 1.0 } else { 0.0 },
                )?;
                b[i * r + j] = value;
                b[j * r + i] = value;
            }
        }
        let rolling_factor = EquilibratedFactor::new(
            b,
            r,
            config.hertz.local.relative_matrix_tolerance,
            &mut cancelled,
        )?;
        let mut transfer = vec![0.0; r * h];
        let mut rhs = vec![0.0; r];
        for j in 0..h {
            poll(&mut cancelled)?;
            for i in 0..r {
                rhs[i] = number(root_damping[i / 2] * symmetric(rolling_index(i), hertz_index(j)))?;
            }
            let column = rolling_factor.forward(&rhs, &mut cancelled)?;
            for i in 0..r {
                transfer[i * h + j] = column[i];
            }
        }
        let mut effective = vec![0.0; h * h];
        for i in 0..h {
            poll(&mut cancelled)?;
            for j in 0..=i {
                let original = symmetric(hertz_index(i), hertz_index(j));
                let mut value = original;
                for k in 0..r {
                    value = number(value - transfer[k * h + i] * transfer[k * h + j])?;
                }
                // Finite nonnegative rolling damping cannot remove a positive
                // Hertz diagonal exactly. Refuse subtraction that cannot
                // resolve it, rather than treating rounded zero as prescribed
                // motion and silently changing static/kinetic selection.
                if i == j
                    && original > 0.0
                    && value <= (8.0 * f64::EPSILON * (r + 1) as f64) * original
                {
                    return Err(invalid("unresolved eliminated Hertz diagonal"));
                }
                effective[i * h + j] = value;
                effective[j * h + i] = value;
            }
        }
        // The original Hertz owner admits the numerical Schur response. Tiny
        // negative pivots are not silently repaired, even at very stiff damping.
        // Reserve one original native support call per patch for final
        // full-C kinetic admission. Static candidates need no kinetic query.
        // The reserve is part of the original total budget, never extra work.
        if config.hertz.maximum_evaluations <= count {
            return Err(invalid("joint original support audit budget"));
        }
        let hertz_config = MultipleHertzConfig {
            maximum_evaluations: config.hertz.maximum_evaluations - count,
            ..config.hertz
        };
        let hertz = PreparedHertzContacts::new(count, effective, hertz_config, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            count,
            compliance: scaled_compliance_m_per_n,
            damping: rolling_damping_n_per_m,
            root_damping,
            rolling_factor,
            transfer,
            hertz,
            config,
        })
    }

    #[must_use]
    pub const fn contact_count(&self) -> usize {
        self.count
    }
    #[must_use]
    pub fn scaled_compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    #[must_use]
    pub fn effective_hertz_compliance_m_per_n(&self) -> &[f64] {
        self.hertz.scaled_compliance_m_per_n()
    }
    #[must_use]
    pub fn rolling_damping_n_per_m(&self) -> &[f64] {
        &self.damping
    }
    #[must_use]
    pub const fn config(&self) -> MultipleHertzRollingConfig {
        self.config
    }

    /// Preserve each accepted Hertz phase while eliminating all original
    /// rolling moments together. Rolling-induced free motion does not erase
    /// the mechanical owner's distinction between moving and stationary patches.
    pub fn advance_phases(
        &self,
        sites: &[HertzRollingContactInput<'_>],
        phases: &[HertzWrenchPhase],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleHertzRollingStep, MultipleHertzRollingError> {
        poll(&mut cancelled)?;
        if sites.len() != self.count || phases.len() != self.count {
            return Err(invalid("joint contact phase count"));
        }
        let mut laws = Vec::with_capacity(self.count);
        for (site, &phase) in sites.iter().zip(phases) {
            poll(&mut cancelled)?;
            laws.push(site.law.phase_law(phase));
        }
        let inputs: Vec<_> = sites
            .iter()
            .zip(&laws)
            .map(|(site, law)| HertzRollingContactInput { law, ..*site })
            .collect();
        self.advance(&inputs, cancelled)
    }

    /// Reuse the same accepted full response and rolling coefficients for each
    /// trial of the original Hertz laws. No partial response is published.
    pub fn advance(
        &self,
        sites: &[HertzRollingContactInput<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleHertzRollingStep, MultipleHertzRollingError> {
        poll(&mut cancelled)?;
        if sites.len() != self.count
            || sites
                .iter()
                .flat_map(|site| site.free_coordinates_m)
                .any(|x| !x.is_finite())
        {
            return Err(invalid("joint input count or released coordinates"));
        }
        let n = 5 * self.count;
        let h = 3 * self.count;
        let r = 2 * self.count;
        let mut rhs = vec![0.0; r];
        for i in 0..r {
            rhs[i] = number(self.root_damping[i / 2] * sites[i / 2].free_coordinates_m[3 + i % 2])?;
        }
        let mut transformed = self.rolling_factor.forward(&rhs, &mut cancelled)?;
        let mut inputs = Vec::with_capacity(self.count);
        for (index, site) in sites.iter().enumerate() {
            poll(&mut cancelled)?;
            let mut free = [0.0; 3];
            for (axis, value) in free.iter_mut().enumerate() {
                *value = site.free_coordinates_m[axis];
                for k in 0..r {
                    *value =
                        number(*value - self.transfer[k * h + 3 * index + axis] * transformed[k])?;
                }
            }
            let spin = number(free[2] / site.law.parameters().contact_radius_m)?;
            if free[2] != 0.0 && spin == 0.0 {
                return Err(invalid("unrepresentable eliminated spin coordinate"));
            }
            inputs.push(HertzContactInput {
                law: site.law,
                free_displacement_m: [free[0], free[1]],
                free_spin_angle_rad: spin,
            });
        }
        let mut hertz = self.hertz.advance(&inputs, &mut cancelled)?;
        let mut wrench = vec![[0.0; 5]; self.count];
        let mut motion = vec![[0.0; 5]; self.count];
        for i in 0..self.count {
            let step = hertz.contacts[i];
            let radius = sites[i].law.parameters().contact_radius_m;
            wrench[i][..3].copy_from_slice(&[
                step.force_n[0],
                step.force_n[1],
                number(step.spin_moment_nm / radius)?,
            ]);
            motion[i][..3].copy_from_slice(&[
                step.slip_displacement_m[0],
                step.slip_displacement_m[1],
                number(step.spin_angle_rad * radius)?,
            ]);
        }
        // The same forward factor used for the Schur free motion recovers the
        // rolling response after every surface wrench has been solved together.
        for i in 0..r {
            poll(&mut cancelled)?;
            for j in 0..h {
                transformed[i] =
                    number(transformed[i] + self.transfer[i * h + j] * wrench[j / 3][j % 3])?;
            }
        }
        let recovered = self.rolling_factor.backward(&transformed, &mut cancelled)?;
        for i in 0..r {
            let root = self.root_damping[i / 2];
            wrench[i / 2][3 + i % 2] = number(-root * recovered[i])?;
            if root > 0.0 {
                motion[i / 2][3 + i % 2] = number(recovered[i] / root)?;
            }
        }
        let mut actual = vec![[0.0; 5]; self.count];
        let mut max_compatibility = 0.0_f64;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut value = sites[i / 5].free_coordinates_m[i % 5];
            for j in 0..n {
                value = number(value + self.compliance[i * n + j] * wrench[j / 5][j % 5])?;
            }
            actual[i / 5][i % 5] = value;
            // A released rolling coordinate retains every cross-patch recoil;
            // zero damping constrains its effort, never its motion.
            if i % 5 >= 3 && self.damping[i / 5] == 0.0 {
                motion[i / 5][i % 5] = value;
            }
            max_compatibility = max_compatibility.max((motion[i / 5][i % 5] - value).abs());
        }
        let mut rolling = Vec::with_capacity(self.count);
        let mut surface_heat = 0.0;
        let mut bulk_heat = 0.0;
        let mut surface_work = 0.0;
        let mut rolling_work = 0.0;
        let mut max_rolling_residual = 0.0_f64;
        let mut max_surface_force_residual = 0.0_f64;
        let mut max_surface_moment_residual = 0.0_f64;
        let mut max_patch_work_residual = 0.0_f64;
        for i in 0..self.count {
            poll(&mut cancelled)?;
            let mut h_work = 0.0;
            let mut r_work = 0.0;
            let mut dissipation = 0.0;
            let mut residual = [0.0; 2];
            for j in 0..3 {
                h_work = number(h_work - wrench[i][j] * actual[i][j])?;
            }
            let mut h_heat = 0.0;
            if hertz.contacts[i].mode == HertzWrenchMode::Sliding {
                let actual_hertz = [actual[i][0], actual[i][1], actual[i][2]];
                let native_input = HertzWrenchInput {
                    free_displacement_m: [actual[i][0], actual[i][1]],
                    free_spin_angle_rad: actual[i][2] / sites[i].law.parameters().contact_radius_m,
                    scaled_compliance_m_per_n: [[0.0; 3]; 3],
                };
                let mut problem = Problem {
                    law: sites[i].law,
                    input: &native_input,
                    free: actual_hertz,
                    config: HertzWrenchConfig {
                        maximum_evaluations: 1,
                        ..self.config.hertz.local
                    },
                    evaluations: 0,
                    cancelled: &mut cancelled,
                };
                let support = problem.support(actual_hertz).map_err(|error| {
                    MultipleHertzRollingError::from(MultipleHertzError::Contact { index: i, error })
                })?;
                hertz.evaluations += problem.evaluations;
                hertz.contacts[i].evaluations += problem.evaluations;
                max_surface_force_residual = max_surface_force_residual.max(norm([
                    number(wrench[i][0] + support[0])?,
                    number(wrench[i][1] + support[1])?,
                ])?);
                max_surface_moment_residual = max_surface_moment_residual.max(number(
                    (wrench[i][2] + support[2]).abs() * sites[i].law.parameters().contact_radius_m,
                )?);
                for j in 0..3 {
                    h_heat = number(h_heat + support[j] * actual_hertz[j])?;
                }
            }
            // Static candidates retain the original feasible reaction at zero
            // constitutive motion. Tiny full-C arithmetic residuals do not
            // spuriously query the discontinuous kinetic law at a new direction.
            for j in 0..2 {
                r_work = number(r_work - wrench[i][j + 3] * actual[i][j + 3])?;
                let root_motion = number(self.root_damping[i] * motion[i][j + 3])?;
                dissipation = number(dissipation + root_motion * root_motion)?;
                residual[j] = number(wrench[i][j + 3] + self.damping[i] * actual[i][j + 3])?;
            }
            max_rolling_residual = max_rolling_residual.max(norm(residual)?);
            let rolling_residual = number(dissipation - r_work)?;
            max_patch_work_residual = max_patch_work_residual
                .max(rolling_residual.abs())
                .max((h_heat - h_work).abs())
                .max((h_heat - hertz.contacts[i].dissipated_energy_j).abs());
            surface_heat = number(surface_heat + h_heat)?;
            bulk_heat = number(bulk_heat + dissipation)?;
            surface_work = number(surface_work + h_work)?;
            rolling_work = number(rolling_work + r_work)?;
            rolling.push(HertzRollingStep {
                scaled_displacement_m: [motion[i][3], motion[i][4]],
                scaled_force_n: [wrench[i][3], wrench[i][4]],
                dissipated_energy_j: dissipation,
                contact_work_j: r_work,
                work_residual_j: rolling_residual,
            });
        }
        let surface_residual = number(surface_heat - surface_work)?;
        let rolling_residual = number(bulk_heat - rolling_work)?;
        let work = number(surface_work + rolling_work)?;
        let total_residual = number(surface_heat + bulk_heat - work)?;
        let tol = self.config;
        if max_compatibility > tol.hertz.local.displacement_tolerance_m
            || max_rolling_residual > tol.rolling_force_tolerance_n
            || max_surface_force_residual > tol.hertz.local.force_tolerance_n
            || max_surface_moment_residual > tol.hertz.local.moment_tolerance_nm
            || max_patch_work_residual > tol.hertz.local.energy_tolerance_j
            || surface_residual.abs() > tol.hertz.total_energy_tolerance_j
            || rolling_residual.abs() > tol.total_energy_tolerance_j
            || total_residual.abs() > tol.total_energy_tolerance_j
        {
            return Err(MultipleHertzRollingError::Audit {
                maximum_compatibility_residual_m: max_compatibility,
                maximum_rolling_force_residual_n: max_rolling_residual,
                maximum_surface_force_residual_n: max_surface_force_residual,
                maximum_surface_spin_moment_residual_nm: max_surface_moment_residual,
                maximum_patch_work_residual_j: max_patch_work_residual,
                surface_work_residual_j: surface_residual,
                rolling_work_residual_j: rolling_residual,
                total_work_residual_j: total_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(MultipleHertzRollingStep {
            surface_heat_j: surface_heat,
            hertz,
            rolling,
            scaled_motion_m: motion,
            scaled_acting_wrench_n: wrench,
            maximum_compatibility_residual_m: max_compatibility,
            maximum_rolling_force_residual_n: max_rolling_residual,
            maximum_surface_force_residual_n: max_surface_force_residual,
            maximum_surface_spin_moment_residual_nm: max_surface_moment_residual,
            bulk_heat_j: bulk_heat,
            surface_work_j: surface_work,
            rolling_work_j: rolling_work,
            surface_work_residual_j: surface_residual,
            rolling_work_residual_j: rolling_residual,
            contact_work_j: work,
            work_residual_j: total_residual,
        })
    }
}

/// B=S*L*L^T*S with S=sqrt(diag B). Equilibration keeps an inactive rolling
/// coordinate (Bii=1) well scaled beside a strongly damped coordinate.
#[derive(Clone, Debug, PartialEq)]
struct EquilibratedFactor {
    n: usize,
    scale: Vec<f64>,
    lower: Vec<f64>,
}
impl EquilibratedFactor {
    fn new(
        mut matrix: Vec<f64>,
        n: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MultipleHertzRollingError> {
        let mut scale = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            if matrix[i * n + i] < 1.0 {
                return Err(invalid("unresolved rolling elimination diagonal"));
            }
            scale.push(number(det::sqrt(matrix[i * n + i]))?);
        }
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..=i {
                let mut value = number((matrix[i * n + j] / scale[i]) / scale[j])?;
                for k in 0..j {
                    value = number(value - matrix[i * n + k] * matrix[j * n + k])?;
                }
                matrix[i * n + j] = if i == j {
                    if value <= tolerance {
                        return Err(invalid("unresolved rolling elimination pivot"));
                    }
                    number(det::sqrt(value))?
                } else {
                    number(value / matrix[j * n + j])?
                };
            }
        }
        Ok(Self {
            n,
            scale,
            lower: matrix,
        })
    }
    fn forward(
        &self,
        rhs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, MultipleHertzRollingError> {
        let mut result = vec![0.0; self.n];
        for i in 0..self.n {
            poll(cancelled)?;
            let mut value = number(rhs[i] / self.scale[i])?;
            for j in 0..i {
                value = number(value - self.lower[i * self.n + j] * result[j])?;
            }
            result[i] = number(value / self.lower[i * self.n + i])?;
        }
        Ok(result)
    }
    fn backward(
        &self,
        rhs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, MultipleHertzRollingError> {
        let mut result = rhs.to_vec();
        for i in (0..self.n).rev() {
            poll(cancelled)?;
            let mut value = result[i];
            for j in i + 1..self.n {
                value = number(value - self.lower[j * self.n + i] * result[j])?;
            }
            result[i] = number(value / self.lower[i * self.n + i])?;
        }
        for (value, &scale) in result.iter_mut().zip(&self.scale) {
            *value = number(*value / scale)?;
        }
        Ok(result)
    }
}

fn hertz_index(index: usize) -> usize {
    5 * (index / 3) + index % 3
}
fn rolling_index(index: usize) -> usize {
    5 * (index / 2) + 3 + index % 2
}
fn number(value: f64) -> Result<f64, MultipleHertzRollingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite joint Hertz/rolling arithmetic"))
    }
}
fn norm(value: [f64; 2]) -> Result<f64, MultipleHertzRollingError> {
    let scale = value[0].abs().max(value[1].abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    number(scale * det::sqrt((value[0] / scale).powi(2) + (value[1] / scale).powi(2)))
}
fn invalid(field: &'static str) -> MultipleHertzRollingError {
    MultipleHertzRollingError::Invalid(field)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultipleHertzRollingError> {
    if cancelled() {
        Err(MultipleHertzRollingError::Cancelled)
    } else {
        Ok(())
    }
}
