//! Simultaneous Hertz shear/spin on one complete mechanical response.
//!
//! Each contact keeps the original three conjugate coordinates and Hertz
//! surface. Off-diagonal blocks are physical recoil, retained throughout
//! `x = x_free + C*w_acting`. A prepared positive-definite matrix supplies the
//! simultaneous zero-slip reaction first. If that reaction is not admitted by
//! every original static surface, deterministic block Gauss--Seidel starts from
//! zero wrench and calls the original three-axis solver. No sweep advances time.
//!
//! Every accepted patch is checked again against the common final motion using
//! the original analytic support or static membership, and all physical work
//! closes against the unchanged full matrix. Static and kinetic coefficients
//! may differ; this bounded root selection does not establish global uniqueness
//! or convergence. Canonical physical ordering belongs to the mechanical owner.
//! This rung admits full positive-definite response or exact zero response.
//! Rank-deficient nonzero response is refused, never regularized.

use core::fmt;

use fs_math::det;

use super::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchStep, Matrix, Problem, Vector, admit_matrix, dot,
};

/// Explicit matrix, iteration, native support and physical closure budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleHertzConfig {
    /// Contact ceiling, from 1 through 128; every patch has three coordinates.
    pub maximum_contacts: usize,
    /// Ceiling on full row-major matrix entries, from 9 through 147456.
    /// Preparation retains two dense matrices and O(contact count) state.
    pub maximum_matrix_terms: usize,
    /// Complete block sweeps, from 1 through 4096.
    pub maximum_sweeps: usize,
    /// Total original support evaluations, including final global admission.
    /// Between 1 and 1048576. Each local solve receives only the remaining budget.
    pub maximum_evaluations: usize,
    /// Original per-patch numerical controls and accepted physical effort limits.
    pub local: HertzWrenchConfig,
    /// Absolute sum of native heat minus full-matrix physical work [J].
    pub total_energy_tolerance_j: f64,
}

impl MultipleHertzConfig {
    pub fn validate(self) -> Result<(), MultipleHertzError> {
        if !(1..=128).contains(&self.maximum_contacts)
            || !(9..=147456).contains(&self.maximum_matrix_terms)
            || !(1..=4096).contains(&self.maximum_sweeps)
            || !(1..=1048576).contains(&self.maximum_evaluations)
            || !self.total_energy_tolerance_j.is_finite()
            || self.total_energy_tolerance_j < 0.0
        {
            return Err(MultipleHertzError::Invalid("multiple Hertz budgets"));
        }
        self.local
            .validate()
            .map_err(|_| MultipleHertzError::Invalid("original Hertz controls"))
    }
}

/// One original law and released motion, in the prepared contact/axis order.
#[derive(Clone, Copy, Debug)]
pub struct HertzContactInput<'a> {
    pub law: &'a HertzWrenchContact,
    pub free_displacement_m: [f64; 2],
    pub free_spin_angle_rad: f64,
}

/// Every contact and the complete physical work belong to one candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipleHertzStep {
    pub contacts: Vec<HertzWrenchStep>,
    pub sweeps: usize,
    pub evaluations: usize,
    pub maximum_compatibility_residual_m: f64,
    pub dissipated_energy_j: f64,
    /// Work of all wrenches against the actual full-matrix common motion [J].
    pub contact_work_j: f64,
    pub work_residual_j: f64,
}

/// No failure publishes an individual patch or changes accepted body state.
#[derive(Clone, Debug, PartialEq)]
pub enum MultipleHertzError {
    Invalid(&'static str),
    Contact {
        index: usize,
        error: HertzWrenchError,
    },
    Cancelled,
    Unresolved {
        sweeps: usize,
        evaluations: usize,
        displacement_residual_m: f64,
        force_residual_n: f64,
        spin_moment_residual_nm: f64,
        energy_residual_j: f64,
    },
}

impl fmt::Display for MultipleHertzError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "simultaneous Hertz contacts: {self:?}")
    }
}
impl std::error::Error for MultipleHertzError {}

/// Immutable physical full response, factored once for repeated source trials.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedHertzContacts {
    count: usize,
    compliance: Vec<f64>,
    diagonal: Vec<Matrix>,
    factor: Option<PositiveFactor>,
    config: MultipleHertzConfig,
}

impl PreparedHertzContacts {
    /// Admit the complete row-major `(3*count)^2` radius-scaled response [m/N].
    /// The physical matrix is never replaced by its symmetric numerical factor.
    pub fn new(
        count: usize,
        scaled_compliance_m_per_n: Vec<f64>,
        config: MultipleHertzConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MultipleHertzError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if count == 0 || count > config.maximum_contacts {
            return Err(MultipleHertzError::Invalid("contact count"));
        }
        let n = 3 * count;
        if n * n > config.maximum_matrix_terms || scaled_compliance_m_per_n.len() != n * n {
            return Err(MultipleHertzError::Invalid(
                "full compliance shape/term budget",
            ));
        }
        let factor = PositiveFactor::new(
            &scaled_compliance_m_per_n,
            n,
            config.local.relative_matrix_tolerance,
            &mut cancelled,
        )?;
        let mut diagonal = Vec::with_capacity(count);
        for index in 0..count {
            poll(&mut cancelled)?;
            let block = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    scaled_compliance_m_per_n[(3 * index + i) * n + 3 * index + j]
                })
            });
            admit_matrix(block, config.local).map_err(|error| contact_error(index, error))?;
            diagonal.push(block);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            count,
            compliance: scaled_compliance_m_per_n,
            diagonal,
            factor,
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
    pub const fn config(&self) -> MultipleHertzConfig {
        self.config
    }

    /// Resolve the original laws from one common released mechanical state.
    /// An infeasible global static seed is discarded before the kinetic search.
    /// Intermediate local effort limits are unconstraining; the caller's limits
    /// apply to every final accepted physical wrench after collective closure.
    pub fn advance(
        &self,
        sites: &[HertzContactInput<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleHertzStep, MultipleHertzError> {
        poll(&mut cancelled)?;
        if sites.len() != self.count {
            return Err(MultipleHertzError::Invalid("contact input count"));
        }
        let n = 3 * self.count;
        let mut free = Vec::with_capacity(n);
        for site in sites {
            poll(&mut cancelled)?;
            let spin = number(site.free_spin_angle_rad * site.law.parameters.contact_radius_m)?;
            if !site.free_spin_angle_rad.is_finite()
                || (site.free_spin_angle_rad != 0.0 && spin == 0.0)
            {
                return Err(MultipleHertzError::Invalid("released spin coordinate"));
            }
            free.extend_from_slice(&[
                number(site.free_displacement_m[0])?,
                number(site.free_displacement_m[1])?,
                spin,
            ]);
        }
        let seed = match &self.factor {
            Some(factor) => factor.reaction(&free, &mut cancelled)?,
            None => None,
        };
        // A displacement tolerance must never turn prescribed nonzero motion
        // into stiction when the physical response is exactly zero.
        let try_static = seed.is_some() || free.iter().all(|&value| value == 0.0);
        let mut forces = seed.unwrap_or_else(|| vec![0.0; n]);
        let mut evaluations = 0;
        let mut modes: Vec<_> = sites
            .iter()
            .map(|site| {
                if site.law.static_patch.sliding_capacity_n() == 0.0 {
                    HertzWrenchMode::Frictionless
                } else {
                    HertzWrenchMode::Sticking
                }
            })
            .collect();
        let mut iterations = vec![0; self.count];
        let mut local_evaluations = vec![0; self.count];
        let mut last = Residual::default();
        if try_static
            && let Some(candidate) = self.audit(
                sites,
                &free,
                &forces,
                &modes,
                &iterations,
                &mut local_evaluations,
                0,
                &mut evaluations,
                &mut last,
                &mut cancelled,
            )?
        {
            return Ok(candidate);
        }
        forces.fill(0.0);
        for sweep in 1..=self.config.maximum_sweeps {
            for (index, site) in sites.iter().enumerate() {
                poll(&mut cancelled)?;
                let remaining = self.config.maximum_evaluations - evaluations;
                if remaining == 0 {
                    return Err(last.error(sweep - 1, evaluations));
                }
                let local_free = self.local_free(&free, &forces, index)?;
                let input = HertzWrenchInput {
                    free_displacement_m: [local_free[0], local_free[1]],
                    free_spin_angle_rad: local_free[2] / site.law.parameters.contact_radius_m,
                    scaled_compliance_m_per_n: self.diagonal[index],
                };
                let local = HertzWrenchConfig {
                    maximum_evaluations: remaining.min(self.config.local.maximum_evaluations),
                    maximum_force_n: f64::MAX,
                    maximum_spin_moment_nm: f64::MAX,
                    ..self.config.local
                };
                let step = site
                    .law
                    .advance(&input, local, &mut cancelled)
                    .map_err(|error| contact_error(index, error))?;
                evaluations += step.evaluations;
                local_evaluations[index] += step.evaluations;
                iterations[index] += step.iterations;
                modes[index] = step.mode;
                forces[3 * index..3 * index + 3].copy_from_slice(&[
                    -step.force_n[0],
                    -step.force_n[1],
                    -step.spin_moment_nm / site.law.parameters.contact_radius_m,
                ]);
            }
            if let Some(candidate) = self.audit(
                sites,
                &free,
                &forces,
                &modes,
                &iterations,
                &mut local_evaluations,
                sweep,
                &mut evaluations,
                &mut last,
                &mut cancelled,
            )? {
                return Ok(candidate);
            }
        }
        Err(last.error(self.config.maximum_sweeps, evaluations))
    }

    fn local_free(
        &self,
        free: &[f64],
        forces: &[f64],
        index: usize,
    ) -> Result<Vector, MultipleHertzError> {
        let n = free.len();
        let mut local = [0.0; 3];
        for axis in 0..3 {
            let row = 3 * index + axis;
            let mut value = free[row];
            for (column, &force) in forces.iter().enumerate() {
                if column / 3 != index {
                    value = number(value - self.compliance[row * n + column] * force)?;
                }
            }
            local[axis] = value;
        }
        Ok(local)
    }

    #[allow(clippy::too_many_arguments)]
    fn audit(
        &self,
        sites: &[HertzContactInput<'_>],
        free: &[f64],
        forces: &[f64],
        modes: &[HertzWrenchMode],
        iterations: &[usize],
        local_evaluations: &mut [usize],
        sweeps: usize,
        evaluations: &mut usize,
        last: &mut Residual,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<MultipleHertzStep>, MultipleHertzError> {
        let n = free.len();
        let mut contacts = Vec::with_capacity(self.count);
        let mut heat = 0.0;
        let mut work = 0.0;
        let mut resolved = true;
        *last = Residual::default();
        for (index, site) in sites.iter().enumerate() {
            poll(cancelled)?;
            let local_free = self.local_free(free, forces, index)?;
            let input = HertzWrenchInput {
                free_displacement_m: [local_free[0], local_free[1]],
                free_spin_angle_rad: local_free[2] / site.law.parameters.contact_radius_m,
                scaled_compliance_m_per_n: self.diagonal[index],
            };
            let wrench = [
                forces[3 * index],
                forces[3 * index + 1],
                forces[3 * index + 2],
            ];
            let mut actual = [0.0; 3];
            for axis in 0..3 {
                let row = 3 * index + axis;
                actual[axis] = number(
                    free[row]
                        - self.compliance[row * n..(row + 1) * n]
                            .iter()
                            .zip(forces)
                            .map(|(&a, &b)| a * b)
                            .sum::<f64>(),
                )?;
            }
            let remaining = self.config.maximum_evaluations - *evaluations;
            if remaining == 0 {
                return Err(last.error(sweeps, *evaluations));
            }
            let mut problem = Problem {
                law: site.law,
                input: &input,
                free: local_free,
                config: HertzWrenchConfig {
                    maximum_evaluations: remaining.min(self.config.local.maximum_evaluations),
                    ..self.config.local
                },
                evaluations: 0,
                cancelled,
            };
            let (motion, support) = match modes[index] {
                HertzWrenchMode::Sticking => {
                    if !problem
                        .static_feasible(wrench)
                        .map_err(|error| contact_error(index, error))?
                    {
                        *evaluations += problem.evaluations;
                        local_evaluations[index] += problem.evaluations;
                        return Ok(None);
                    }
                    ([0.0; 3], wrench)
                }
                HertzWrenchMode::Sliding => {
                    let support = problem
                        .support(actual)
                        .map_err(|error| contact_error(index, error))?;
                    (actual, support)
                }
                HertzWrenchMode::Frictionless => {
                    // Exact zero capacity cannot acquire a tiny reaction from
                    // the static seed merely because it is below a tolerance.
                    if wrench != [0.0; 3] {
                        return Ok(None);
                    }
                    (actual, [0.0; 3])
                }
            };
            let mut step = problem
                .candidate(wrench, motion, support, modes[index], iterations[index])
                .map_err(|error| contact_error(index, error))?;
            // Reconstruct compatibility from the ORIGINAL FULL matrix, not the
            // block's rounded released motion or its symmetric numerical factor.
            step.compatibility_residual_m =
                core::array::from_fn(|axis| motion[axis] - actual[axis]);
            step.body_work_residual_j = number(dot(wrench, step.compatibility_residual_m))?;
            *evaluations += problem.evaluations;
            local_evaluations[index] += problem.evaluations;
            step.evaluations = local_evaluations[index];
            last.displacement = last.displacement.max(
                step.compatibility_residual_m
                    .into_iter()
                    .map(f64::abs)
                    .fold(0.0, f64::max),
            );
            last.force = last.force.max(step.force_residual_n);
            last.moment = last.moment.max(step.spin_moment_residual_nm);
            last.energy = last
                .energy
                .max(step.work_residual_j.abs())
                .max(step.duality_gap_j.abs())
                .max(step.body_work_residual_j.abs());
            resolved &= problem.resolved(&step);
            heat = number(heat + step.dissipated_energy_j)?;
            work = number(work + dot(wrench, actual))?;
            contacts.push(step);
        }
        let residual = number(heat - work)?;
        last.energy = last.energy.max(residual.abs());
        if !resolved || residual.abs() > self.config.total_energy_tolerance_j {
            return Ok(None);
        }
        // Effort budgets are physical accepted-wrench budgets, never a limiter
        // on intermediate block reactions that may disappear on shared bodies.
        for (index, step) in contacts.iter().enumerate() {
            poll(cancelled)?;
            let force = det::hypot(step.force_n[0], step.force_n[1]);
            if force > self.config.local.maximum_force_n
                || step.spin_moment_nm.abs() > self.config.local.maximum_spin_moment_nm
            {
                return Err(contact_error(
                    index,
                    HertzWrenchError::EffortBudget {
                        force_n: force,
                        spin_moment_nm: step.spin_moment_nm.abs(),
                    },
                ));
            }
        }
        poll(cancelled)?;
        Ok(Some(MultipleHertzStep {
            contacts,
            sweeps,
            evaluations: *evaluations,
            maximum_compatibility_residual_m: last.displacement,
            dissipated_energy_j: heat,
            contact_work_j: work,
            work_residual_j: residual,
        }))
    }
}

/// Normalized Cholesky admission/seed only; no numerical pivot is added to C.
#[derive(Clone, Debug, PartialEq)]
struct PositiveFactor {
    n: usize,
    scale: f64,
    lower: Vec<f64>,
}
impl PositiveFactor {
    fn new(
        matrix: &[f64],
        n: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Self>, MultipleHertzError> {
        let mut scale = 0.0_f64;
        for row in matrix.chunks(n) {
            poll(cancelled)?;
            for &value in row {
                scale = scale.max(number(value)?.abs());
            }
        }
        if scale == 0.0 {
            return Ok(None);
        }
        let mut lower = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..=i {
                let a = matrix[i * n + j] / scale;
                let b = matrix[j * n + i] / scale;
                if (a - b).abs() > tolerance {
                    return Err(MultipleHertzError::Invalid("nonsymmetric full compliance"));
                }
                let mut value = 0.5 * a + 0.5 * b;
                for k in 0..j {
                    value -= lower[i * n + k] * lower[j * n + k];
                }
                lower[i * n + j] = if i == j {
                    if !value.is_finite() || value <= tolerance {
                        return Err(MultipleHertzError::Invalid(
                            "nonpositive-definite full compliance",
                        ));
                    }
                    det::sqrt(value)
                } else {
                    number(value / lower[j * n + j])?
                };
            }
        }
        Ok(Some(Self { n, scale, lower }))
    }

    fn reaction(
        &self,
        free: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Vec<f64>>, MultipleHertzError> {
        let mut force = vec![0.0; self.n];
        for i in 0..self.n {
            poll(cancelled)?;
            let mut value = free[i] / self.scale;
            for j in 0..i {
                value -= self.lower[i * self.n + j] * force[j];
            }
            force[i] = value / self.lower[i * self.n + i];
        }
        for i in (0..self.n).rev() {
            poll(cancelled)?;
            let mut value = force[i];
            for j in i + 1..self.n {
                value -= self.lower[j * self.n + i] * force[j];
            }
            force[i] = value / self.lower[i * self.n + i];
        }
        Ok(force.iter().all(|v| v.is_finite()).then_some(force))
    }
}

#[derive(Clone, Copy, Default)]
struct Residual {
    displacement: f64,
    force: f64,
    moment: f64,
    energy: f64,
}
impl Residual {
    fn error(self, sweeps: usize, evaluations: usize) -> MultipleHertzError {
        MultipleHertzError::Unresolved {
            sweeps,
            evaluations,
            displacement_residual_m: self.displacement,
            force_residual_n: self.force,
            spin_moment_residual_nm: self.moment,
            energy_residual_j: self.energy,
        }
    }
}
fn contact_error(index: usize, error: HertzWrenchError) -> MultipleHertzError {
    match error {
        HertzWrenchError::Cancelled => MultipleHertzError::Cancelled,
        error => MultipleHertzError::Contact { index, error },
    }
}
fn number(value: f64) -> Result<f64, MultipleHertzError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MultipleHertzError::Invalid(
            "nonfinite multiple-contact arithmetic",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultipleHertzError> {
    if cancelled() {
        Err(MultipleHertzError::Cancelled)
    } else {
        Ok(())
    }
}
