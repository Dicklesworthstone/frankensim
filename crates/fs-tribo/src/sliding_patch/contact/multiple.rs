//! Simultaneous Hertz shear/spin on one complete mechanical response.
//!
//! Each contact keeps the original three conjugate coordinates and Hertz
//! surface. Off-diagonal blocks are physical recoil, retained throughout
//! `x = x_free + C*w_acting`. A prepared positive-semidefinite matrix supplies a
//! supported minimum-norm zero-slip seed. Dependent contact coordinates may
//! redistribute this reaction inside the original static surfaces. A bounded
//! static-feasibility phase must find such a reaction or separate the original
//! static sets from the requested motion before kinetic alternatives are tried.
//! Deterministic block solves call the original three-axis solver throughout.
//! No sweep advances time or deposits heat.
//!
//! Every accepted patch is checked again against the common final motion using
//! the original analytic support or static membership, and all physical work
//! closes against the unchanged full matrix. Static and kinetic coefficients
//! may differ; this bounded root selection does not establish global uniqueness
//! or convergence. Canonical physical ordering belongs to the mechanical owner.
//! Full response may be rank deficient, including exactly zero. Each local
//! three-axis block must retain the original positive-definite/exact-zero
//! admission. No compliance, viscosity or stored contact energy is invented.
//! Static separation uses an explicit work margin and floating-point rounding
//! guard; it is a bounded numerical decision, not an interval proof.

use core::fmt;

pub mod rolling;

use fs_math::det;

use crate::contact_response::{ContactResponseError, PositiveFactor};

use super::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchStep, Matrix, Problem, SPIN_RADIUS, Vector, admit_matrix, dot,
};

/// Explicit matrix, iteration, native support and physical closure budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleHertzConfig {
    /// Contact ceiling, from 1 through 128; every patch has three coordinates.
    pub maximum_contacts: usize,
    /// Ceiling on full row-major matrix entries, from 9 through 147456.
    /// Preparation retains at most three dense matrices plus O(contact count).
    pub maximum_matrix_terms: usize,
    /// Complete block sweeps shared by static feasibility and the kinetic solve,
    /// from 1 through 4096.
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
    /// Dependent contacts may admit a different all-static load allocation, but
    /// neither feasibility nor a separating static-support margin was resolved.
    /// A lower kinetic branch is never substituted for this missing decision.
    StaticFeasibilityUnresolved {
        sweeps: usize,
        evaluations: usize,
        displacement_residual_m: f64,
        separation_margin_j: f64,
    },
    Unresolved {
        sweeps: usize,
        evaluations: usize,
        displacement_residual_m: f64,
        force_residual_n: f64,
        spin_moment_residual_nm: f64,
        energy_residual_j: f64,
    },
}

impl From<ContactResponseError> for MultipleHertzError {
    fn from(error: ContactResponseError) -> Self {
        match error {
            ContactResponseError::Invalid(field) => Self::Invalid(field),
            ContactResponseError::Cancelled => Self::Cancelled,
        }
    }
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
    factor: PositiveFactor,
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
    /// Dependent response resolves possible static load redistribution before
    /// admitting kinetic alternatives; unresolved feasibility explicitly refuses.
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
        let mut forces = self.factor.reaction(&free, &mut cancelled)?;
        let mut evaluations = 0;
        let static_modes: Vec<_> = sites
            .iter()
            .map(|site| {
                if site.law.static_patch.sliding_capacity_n() == 0.0 {
                    HertzWrenchMode::Frictionless
                } else {
                    HertzWrenchMode::Sticking
                }
            })
            .collect();
        let mut modes = static_modes.clone();
        let mut iterations = vec![0; self.count];
        let mut local_evaluations = vec![0; self.count];
        let mut last = Residual::default();
        let (supported, _) = self.static_supported(sites, &free, &forces)?;
        if supported
            && let Some(candidate) = self.audit(
                sites,
                &free,
                &forces,
                &static_modes,
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
        let mut spent_sweeps = 0;
        let reduced_static_constraints = sites
            .iter()
            .any(|site| site.law.static_patch.sliding_capacity_n() == 0.0);
        if self.factor.rank() > 0
            && (self.factor.rank() < n || !supported || reduced_static_constraints)
        {
            // A minimum-norm reaction outside a static set does not prove that
            // every load allocation in the nullspace is outside those sets.
            // Nor does an inaccurate full-rank solve prove infeasibility, or a
            // seed that incorrectly constrained a zero-capacity row to rest.
            // Equal static/kinetic coefficients here are a numerical projection
            // onto the ORIGINAL static sets; these temporary laws never supply
            // the returned physical mode, work or heat.
            let projection_laws: Vec<_> = sites
                .iter()
                .enumerate()
                .map(|(index, site)| {
                    let mut parameters = site.law.parameters();
                    parameters.kinetic_mu = parameters.static_mu;
                    HertzWrenchContact::new(parameters).map_err(|error| contact_error(index, error))
                })
                .collect::<Result<_, _>>()?;
            let projection_sites: Vec<_> = sites
                .iter()
                .zip(&projection_laws)
                .map(|(site, law)| HertzContactInput {
                    law,
                    free_displacement_m: site.free_displacement_m,
                    free_spin_angle_rad: site.free_spin_angle_rad,
                })
                .collect();
            let mut margin = 0.0;
            let mut separated = self.static_separated(
                sites,
                &free,
                &forces,
                &mut evaluations,
                &mut local_evaluations,
                &mut margin,
                &mut cancelled,
            )?;
            while !separated && spent_sweeps < self.config.maximum_sweeps {
                if evaluations == self.config.maximum_evaluations {
                    break;
                }
                self.block_sweep(
                    &projection_sites,
                    &free,
                    &mut forces,
                    &mut modes,
                    &mut iterations,
                    &mut local_evaluations,
                    &mut evaluations,
                    spent_sweeps,
                    &last,
                    &mut cancelled,
                )?;
                spent_sweeps += 1;
                let (supported, residual) = self.static_supported(sites, &free, &forces)?;
                last.displacement = residual;
                if supported
                    && let Some(candidate) = self.audit(
                        sites,
                        &free,
                        &forces,
                        &static_modes,
                        &iterations,
                        &mut local_evaluations,
                        spent_sweeps,
                        &mut evaluations,
                        &mut last,
                        &mut cancelled,
                    )?
                {
                    return Ok(candidate);
                }
                separated = self.static_separated(
                    sites,
                    &free,
                    &forces,
                    &mut evaluations,
                    &mut local_evaluations,
                    &mut margin,
                    &mut cancelled,
                )?;
            }
            if !separated {
                return Err(MultipleHertzError::StaticFeasibilityUnresolved {
                    sweeps: spent_sweeps,
                    evaluations,
                    displacement_residual_m: last.displacement,
                    separation_margin_j: margin,
                });
            }
            forces.fill(0.0);
            iterations.fill(0);
        }
        for sweep in spent_sweeps + 1..=self.config.maximum_sweeps {
            self.block_sweep(
                sites,
                &free,
                &mut forces,
                &mut modes,
                &mut iterations,
                &mut local_evaluations,
                &mut evaluations,
                sweep - 1,
                &last,
                &mut cancelled,
            )?;
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

    #[allow(clippy::too_many_arguments)]
    fn block_sweep(
        &self,
        sites: &[HertzContactInput<'_>],
        free: &[f64],
        forces: &mut [f64],
        modes: &mut [HertzWrenchMode],
        iterations: &mut [usize],
        local_evaluations: &mut [usize],
        evaluations: &mut usize,
        sweeps: usize,
        last: &Residual,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MultipleHertzError> {
        for (index, site) in sites.iter().enumerate() {
            poll(cancelled)?;
            let remaining = self.config.maximum_evaluations - *evaluations;
            if remaining == 0 {
                return Err(last.error(sweeps, *evaluations));
            }
            let free = self.local_free(free, forces, index)?;
            let input = HertzWrenchInput {
                free_displacement_m: [free[0], free[1]],
                free_spin_angle_rad: free[2] / site.law.parameters.contact_radius_m,
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
                .advance(&input, local, &mut *cancelled)
                .map_err(|error| contact_error(index, error))?;
            *evaluations += step.evaluations;
            local_evaluations[index] += step.evaluations;
            iterations[index] += step.iterations;
            modes[index] = step.mode;
            forces[3 * index..3 * index + 3].copy_from_slice(&[
                -step.force_n[0],
                -step.force_n[1],
                -step.spin_moment_nm / site.law.parameters.contact_radius_m,
            ]);
        }
        Ok(())
    }

    // A user displacement tolerance must not erase a small unsupported drive
    // on a nullspace. All-static admission also requires scale-aware arithmetic
    // closure of the original physical rows. Zero-capacity rows may move freely.
    fn static_supported(
        &self,
        sites: &[HertzContactInput<'_>],
        free: &[f64],
        forces: &[f64],
    ) -> Result<(bool, f64), MultipleHertzError> {
        let n = free.len();
        let mut supported = true;
        let mut maximum = 0.0_f64;
        for (row, &released) in free.iter().enumerate() {
            if sites[row / 3].law.static_patch.sliding_capacity_n() == 0.0 {
                continue;
            }
            let mut response = 0.0;
            let mut scale = released.abs();
            for (column, &force) in forces.iter().enumerate() {
                let term = number(self.compliance[row * n + column] * force)?;
                response = number(response + term)?;
                scale = number(scale + term.abs())?;
            }
            let residual = number(released - response)?.abs();
            maximum = maximum.max(residual);
            let rounding = 64.0 * f64::EPSILON * (n + 1) as f64 * scale;
            supported &= residual <= rounding.min(self.config.local.displacement_tolerance_m);
        }
        Ok((supported, maximum))
    }

    // If every active patch could stick, C*w=free on its rows and each w_i is
    // inside its original static set. Therefore y.free <= sum h_i(C^T*y) for
    // every y supported on those rows. Strict violation authorizes the kinetic
    // search; failure to resolve it is NEVER interpreted as static infeasibility.
    #[allow(clippy::too_many_arguments)]
    fn static_separated(
        &self,
        sites: &[HertzContactInput<'_>],
        free: &[f64],
        forces: &[f64],
        evaluations: &mut usize,
        local_evaluations: &mut [usize],
        margin: &mut f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, MultipleHertzError> {
        let n = free.len();
        let mut residual = vec![0.0; n];
        for row in 0..n {
            poll(cancelled)?;
            if sites[row / 3].law.static_patch.sliding_capacity_n() == 0.0 {
                continue;
            }
            residual[row] = number(
                free[row]
                    - self.compliance[row * n..(row + 1) * n]
                        .iter()
                        .zip(forces)
                        .map(|(&a, &b)| a * b)
                        .sum::<f64>(),
            )?;
        }
        let scale = residual.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if scale == 0.0 {
            return Ok(false);
        }
        let factor_direction = self.factor.reaction(&residual, cancelled)?;
        let normalized: Vec<_> = residual.iter().map(|&value| value / scale).collect();
        for mut direction in [factor_direction, normalized] {
            for row in 0..n {
                if sites[row / 3].law.static_patch.sliding_capacity_n() == 0.0 {
                    direction[row] = 0.0;
                }
            }
            let mut lhs = 0.0;
            let mut arithmetic_scale = 0.0;
            let mut dual = vec![0.0; n];
            let mut dual_absolute_sum = vec![0.0; n];
            for row in 0..n {
                poll(cancelled)?;
                let term = number(direction[row] * free[row])?;
                lhs = number(lhs + term)?;
                arithmetic_scale = number(arithmetic_scale + term.abs())?;
                for column in 0..n {
                    // Transpose of the ORIGINAL matrix, including admitted
                    // asymmetry/roundoff, never its symmetric numerical factor.
                    let term = number(self.compliance[row * n + column] * direction[row])?;
                    dual[column] = number(dual[column] + term)?;
                    dual_absolute_sum[column] = number(dual_absolute_sum[column] + term.abs())?;
                }
            }
            let mut support = 0.0;
            for (index, site) in sites.iter().enumerate() {
                poll(cancelled)?;
                // A near-null dual can cancel large matrix products. The
                // original static support is Lipschitz with force bound cap
                // and scaled spin bound cap*(3*pi/16); propagate those column
                // accumulation errors even when the rounded dual is zero.
                let capacity = number(site.law.static_patch.sliding_capacity_n())?;
                arithmetic_scale = number(
                    arithmetic_scale
                        + capacity
                            * (dual_absolute_sum[3 * index]
                                + dual_absolute_sum[3 * index + 1]
                                + SPIN_RADIUS * dual_absolute_sum[3 * index + 2]),
                )?;
                let motion = [dual[3 * index], dual[3 * index + 1], dual[3 * index + 2]];
                if motion == [0.0; 3] || site.law.static_patch.sliding_capacity_n() == 0.0 {
                    continue;
                }
                let remaining = self.config.maximum_evaluations - *evaluations;
                if remaining == 0 {
                    return Ok(false);
                }
                let projection = HertzWrenchContact {
                    parameters: site.law.parameters,
                    static_patch: site.law.static_patch,
                    kinetic_patch: site.law.static_patch,
                };
                let input = HertzWrenchInput {
                    free_displacement_m: [motion[0], motion[1]],
                    free_spin_angle_rad: motion[2] / site.law.parameters.contact_radius_m,
                    scaled_compliance_m_per_n: [[0.0; 3]; 3],
                };
                let mut problem = Problem {
                    law: &projection,
                    input: &input,
                    free: motion,
                    config: HertzWrenchConfig {
                        maximum_evaluations: remaining.min(self.config.local.maximum_evaluations),
                        ..self.config.local
                    },
                    evaluations: 0,
                    cancelled,
                };
                let native = problem
                    .support(motion)
                    .map_err(|error| contact_error(index, error))?;
                *evaluations += problem.evaluations;
                local_evaluations[index] += problem.evaluations;
                let value = number(dot(native, motion))?;
                support = number(support + value)?;
                arithmetic_scale = number(arithmetic_scale + value.abs())?;
            }
            let gap = number(lhs - support)?;
            *margin = margin.max(gap);
            let rounding = 128.0 * f64::EPSILON * (n + 1) as f64 * arithmetic_scale;
            let physical_margin = self.config.total_energy_tolerance_j
                + self.count as f64 * self.config.local.energy_tolerance_j;
            if gap > physical_margin + rounding {
                return Ok(true);
            }
        }
        Ok(false)
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
