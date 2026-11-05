//! Simultaneous maximum-dissipation impulses at scalar friction ports.
//!
//! The original physical response is `v_after = v_before + W*p`. Each port
//! has one explicitly supplied, frozen impulse set `-capacity <= p <= capacity`.
//! The unique candidate minimizes `v_before.dot(p) + p.dot(W*p)/2` over their
//! product. Coupling remains in every solve: independently clamping the
//! unconstrained stopping impulses does not solve this problem.
//!
//! Coordinates must be work-conjugate. A mechanical owner can express an
//! angular port as velocity `L*omega` and impulse `angular_impulse/L`, retaining
//! its authored length `L` throughout assembly and native application. This
//! owner supplies no normal load, coefficient, timestep, static/kinetic phase
//! selection, elastic history, prescribed-source work or heat allocation.

use fs_math::dd::Dd;

/// Hard bound on the number of simultaneous scalar friction ports.
pub const MAX_FRICTION_IMPULSE_PORTS: usize = 32;

type Vector = [f64; MAX_FRICTION_IMPULSE_PORTS];

/// Borrowed physical ports in one unchanged ordering.
#[derive(Clone, Copy, Debug)]
pub struct CoupledFrictionImpulseInput<'a> {
    /// Relative port velocities before the impulse event [m/s].
    pub velocity_before_m_s: &'a [f64],
    /// Complete row-major `n*n` mobility, mapping impulses to velocities [1/kg].
    /// Positive-definite admission includes every port, even a zero-capacity one.
    pub mobility_per_kg: &'a [f64],
    /// Nonnegative symmetric physical impulse capacities [N s]. These are not
    /// numerical effort budgets and are never reduced to satisfy one.
    pub capacity_ns: &'a [f64],
}

/// Numerical work and physical acceptance gates for one atomic candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledFrictionImpulseControls {
    /// Maximum active-set iterations, in `1..=4096`.
    pub maximum_iterations: usize,
    /// Original inverse-mobility reconstruction error for each impulse [N s].
    pub impulse_tolerance_ns: f64,
    /// Original-matrix box KKT residual in each velocity coordinate [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Absolute support-gap and independent work/kinetic closure tolerance [J].
    pub energy_tolerance_j: f64,
    /// Positive normalized symmetry and squared-pivot tolerance, at most 1e-8.
    /// The pivot gate also retains the arithmetic floor `n*EPSILON`.
    pub relative_matrix_tolerance: f64,
    /// Maximum absolute ACTUAL impulse in any port [N s]. A required solution
    /// exceeding this nonnegative budget refuses; its physical set is unchanged.
    pub maximum_impulse_ns: f64,
}

impl CoupledFrictionImpulseControls {
    /// Validate bounded work and dimensional tolerances before any solve.
    pub fn validate(self) -> Result<(), CoupledFrictionImpulseError> {
        if !(1..=4096).contains(&self.maximum_iterations)
            || !positive(self.impulse_tolerance_ns)
            || !positive(self.velocity_tolerance_m_s)
            || !positive(self.energy_tolerance_j)
            || !positive(self.relative_matrix_tolerance)
            || self.relative_matrix_tolerance > 1e-8
            || !self.maximum_impulse_ns.is_finite()
            || self.maximum_impulse_ns < 0.0
        {
            return Err(CoupledFrictionImpulseError::Invalid("impulse controls"));
        }
        Ok(())
    }
}

/// Original-system candidate; only the first `count` array entries are ports.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct CoupledFrictionImpulseReport {
    /// Number of supplied ports; unused array tails are exactly zero.
    pub count: usize,
    /// Physical impulses, each inside its original closed interval [N s].
    pub impulse_ns: Vector,
    /// Actual original `v_before + W*p`, including floating arrest residuals [m/s].
    pub velocity_after_m_s: Vector,
    /// `W^-1*(v_after-v_before)-p`, using an independent factor solve [N s].
    pub impulse_residual_ns: Vector,
    /// Signed box KKT residual of the actual outgoing velocity [m/s]. A
    /// zero-capacity coordinate is fixed and has no velocity stationarity claim.
    pub velocity_residual_m_s: Vector,
    /// Nonnegative total mechanical loss from midpoint impulse work [J].
    /// Individual port work can exchange energy between coupled relative modes.
    pub dissipated_j: f64,
    /// Independently reconstructed inverse-mobility kinetic change [J].
    pub kinetic_energy_change_j: f64,
    /// Kinetic change plus mechanical loss [J].
    pub energy_residual_j: f64,
    /// Original box support gap `sum(capacity*abs(v_after)+p*v_after)` [J].
    /// This is an optimality diagnostic, distinct from mechanical dissipation.
    pub duality_gap_j: f64,
    /// Completed bounded active-set iterations.
    pub iterations: usize,
    /// Full admission plus changed-free-set factorizations. A fixed free set
    /// reuses its factor for original-system refinement.
    pub factorizations: usize,
}

/// Refusal leaves all borrowed input data unchanged and returns no candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledFrictionImpulseError {
    /// Invalid dimensions, physical inputs, controls or matrix symmetry.
    Invalid(&'static str),
    /// Finite inputs produced an unrepresentable intermediate or output.
    NonFinite(&'static str),
    /// A nonpositive or numerically unresolved normalized pivot. No added
    /// compliance or arbitrary null-space representative repairs this mobility.
    Mobility {
        /// Original port index of the failed principal pivot.
        index: usize,
        /// Squared normalized pivot, possibly nonpositive.
        scaled_pivot: f64,
    },
    /// An actual physical solution exceeds the separate numerical effort budget.
    ImpulseBudget {
        /// Port whose returned impulse exceeded the budget.
        index: usize,
        /// Signed actual impulse [N s].
        impulse_ns: f64,
        /// Maximum permitted absolute impulse [N s].
        maximum_ns: f64,
    },
    /// The active-set budget ended before the original KKT and support gates met.
    IterationLimit {
        /// Completed bounded iterations.
        iterations: usize,
        /// Largest original velocity KKT residual [m/s].
        velocity_residual_m_s: f64,
        /// Original box support gap [J].
        duality_gap_j: f64,
    },
    /// Floating resolution could not meet an original physical residual gate.
    Unresolved {
        /// Largest original velocity KKT residual [m/s].
        velocity_residual_m_s: f64,
        /// Largest independently reconstructed impulse error [N s].
        impulse_residual_ns: f64,
        /// Work closure or unresolved support gap [J].
        energy_residual_j: f64,
    },
    /// Actual midpoint work was active or its independent kinetic check failed.
    EnergyBalance {
        /// Candidate mechanical loss [J].
        dissipated_j: f64,
        /// Kinetic change plus that loss [J].
        residual_j: f64,
    },
    /// Caller cancellation, including the final publication checkpoint.
    Cancelled,
}

impl core::fmt::Display for CoupledFrictionImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled scalar friction impulses: {self:?}")
    }
}
impl std::error::Error for CoupledFrictionImpulseError {}

/// Solve every supplied scalar friction port simultaneously.
///
/// The original positive-definite mobility admits one impulse minimizer. A
/// deterministic primal active set resolves its free ports together, advances
/// only to the first physical bound, and releases a bound only for a violating
/// outgoing velocity. Numerical diagonal scaling changes coordinates, not the
/// original physical response. The original-matrix KKT, support gap, impulse
/// reconstruction and independent kinetic-work checks govern publication.
///
/// Work is bounded by `maximum_iterations` principal solves for at most 32
/// rows. Two fixed Cholesky workspaces and fixed vectors require no allocation.
/// Cancellation is polled at matrix/factor/solve rows and immediately before
/// return. This instantaneous mechanical law does not infer finite-time local
/// friction heat, branch history, normal contact or a trajectory.
pub fn solve_coupled_friction_impulses(
    input: &CoupledFrictionImpulseInput<'_>,
    controls: CoupledFrictionImpulseControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CoupledFrictionImpulseReport, CoupledFrictionImpulseError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let count = input.velocity_before_m_s.len();
    if count > MAX_FRICTION_IMPULSE_PORTS
        || input.capacity_ns.len() != count
        || input.mobility_per_kg.len() != count * count
    {
        return Err(CoupledFrictionImpulseError::Invalid("port dimensions"));
    }
    let mut scales = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    for (i, scale) in scales.iter_mut().enumerate().take(count) {
        poll(&mut cancelled)?;
        if !input.velocity_before_m_s[i].is_finite()
            || !input.capacity_ns[i].is_finite()
            || input.capacity_ns[i] < 0.0
        {
            return Err(CoupledFrictionImpulseError::Invalid("velocity or capacity"));
        }
        let diagonal = input.mobility_per_kg[i * count + i];
        if !positive(diagonal) {
            return Err(CoupledFrictionImpulseError::Invalid(
                "positive mobility diagonal",
            ));
        }
        *scale = diagonal.sqrt();
    }
    for i in 0..count {
        poll(&mut cancelled)?;
        for j in 0..i {
            let a = normalized(input.mobility_per_kg[i * count + j], scales[i], scales[j])?;
            let b = normalized(input.mobility_per_kg[j * count + i], scales[j], scales[i])?;
            if (a - b).abs() > controls.relative_matrix_tolerance {
                return Err(CoupledFrictionImpulseError::Invalid("mobility symmetry"));
            }
        }
    }
    let all = if count == 32 {
        u32::MAX
    } else {
        (1_u32 << count) - 1
    };
    let mut full_factor = Factor::empty();
    full_factor.prepare(input, &scales, all, controls, &mut cancelled)?;
    let mut factorizations = usize::from(count != 0);
    let mut free_factor = Factor::empty();
    let mut cached_free = None;
    let mut impulse = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    // -1 lower, +1 upper, 0 free, 2 fixed zero capacity.
    let mut active = [0_i8; MAX_FRICTION_IMPULSE_PORTS];
    for (i, state) in active.iter_mut().enumerate().take(count) {
        if input.capacity_ns[i] == 0.0 {
            *state = 2;
        }
    }
    let mut iterations = 0;
    loop {
        poll(&mut cancelled)?;
        let after = response(input, &impulse, &mut cancelled)?;
        let residual = residual(input, &impulse, &after)?;
        // A velocity tolerance must not create a zero-impulse dead band:
        // initially nonzero reacting motion reaches at least one real solve.
        if (iterations != 0 || residual.maximum_velocity == 0.0)
            && residual.maximum_velocity <= controls.velocity_tolerance_m_s
            && residual.gap_j <= controls.energy_tolerance_j
        {
            return finish(
                input,
                controls,
                &scales,
                &full_factor,
                impulse,
                after,
                residual,
                iterations,
                factorizations,
                &mut cancelled,
            );
        }
        if iterations == controls.maximum_iterations {
            return Err(CoupledFrictionImpulseError::IterationLimit {
                iterations,
                velocity_residual_m_s: residual.maximum_velocity,
                duality_gap_j: residual.gap_j,
            });
        }
        iterations += 1;
        let mut free = 0_u32;
        let mut free_residual = 0.0_f64;
        let mut release = None;
        let mut worst = 0.0;
        for i in 0..count {
            if active[i] == 0 {
                free |= 1_u32 << i;
                free_residual = free_residual.max(after[i].abs());
            } else if active[i] != 2 {
                let violation = (f64::from(active[i]) * after[i]).max(0.0);
                let score = finite(violation / scales[i], "active gradient scaling")?;
                if score > worst {
                    worst = score;
                    release = Some(i);
                }
            }
        }
        if free_residual <= controls.velocity_tolerance_m_s
            && let Some(i) = release
        {
            active[i] = 0;
            continue;
        }
        if free == 0 {
            return Err(unresolved(&residual));
        }
        let factor = if free == all {
            &full_factor
        } else {
            if cached_free != Some(free) {
                free_factor.prepare(input, &scales, free, controls, &mut cancelled)?;
                factorizations += 1;
                cached_free = Some(free);
            }
            &free_factor
        };
        let rhs = after.map(|v| -v);
        let direction = factor.solve(&rhs, &scales, &mut cancelled)?;
        let mut fraction = 1.0_f64;
        let mut blocking = None;
        for i in 0..count {
            if active[i] != 0 || direction[i] == 0.0 {
                continue;
            }
            let distance = if direction[i] > 0.0 {
                input.capacity_ns[i] - impulse[i]
            } else {
                input.capacity_ns[i] + impulse[i]
            };
            if direction[i].abs() >= distance {
                let candidate = distance / direction[i].abs();
                if candidate < fraction || (candidate == fraction && blocking.is_none()) {
                    fraction = candidate;
                    blocking = Some((i, direction[i].is_sign_positive()));
                }
            }
        }
        let mut changed = false;
        for i in 0..count {
            if active[i] != 0 {
                continue;
            }
            // The feasible line search determined the step. Rounding at a
            // face is kept inside the ORIGINAL box; the coupled free equations
            // and original KKT are solved again after every changed face.
            let candidate = finite(
                fraction.mul_add(direction[i], impulse[i]),
                "active-set step",
            )?
            .clamp(-input.capacity_ns[i], input.capacity_ns[i]);
            changed |= candidate != impulse[i];
            impulse[i] = candidate;
            if candidate == -input.capacity_ns[i] {
                active[i] = -1;
            } else if candidate == input.capacity_ns[i] {
                active[i] = 1;
            }
        }
        if let Some((i, upper)) = blocking {
            impulse[i] = if upper {
                input.capacity_ns[i]
            } else {
                -input.capacity_ns[i]
            };
            active[i] = if upper { 1 } else { -1 };
        } else if !changed {
            if let Some(i) = release {
                active[i] = 0;
            } else {
                return Err(unresolved(&residual));
            }
        }
    }
}

struct Residual {
    velocity: Vector,
    maximum_velocity: f64,
    gap_j: f64,
}

fn residual(
    input: &CoupledFrictionImpulseInput<'_>,
    impulse: &Vector,
    after: &Vector,
) -> Result<Residual, CoupledFrictionImpulseError> {
    let mut velocity = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut maximum_velocity = 0.0_f64;
    let mut gap = Dd::ZERO;
    for i in 0..input.velocity_before_m_s.len() {
        let cap = input.capacity_ns[i];
        if impulse[i].abs() > cap {
            return Err(CoupledFrictionImpulseError::Invalid(
                "original impulse feasibility",
            ));
        }
        if cap == 0.0 {
            continue;
        }
        velocity[i] = if impulse[i] == -cap {
            after[i].min(0.0)
        } else if impulse[i] == cap {
            after[i].max(0.0)
        } else {
            after[i]
        };
        maximum_velocity = maximum_velocity.max(velocity[i].abs());
        // Difference from the exact support point avoids subtracting large
        // opposing support/work terms at a saturated friction port.
        let distance = if after[i] >= 0.0 {
            Dd::from_f64(cap) + Dd::from_f64(impulse[i])
        } else {
            Dd::from_f64(cap) - Dd::from_f64(impulse[i])
        };
        gap = gap + distance * Dd::from_f64(after[i].abs());
    }
    let gap_j = finite(gap.to_f64(), "support gap")?;
    Ok(Residual {
        velocity,
        maximum_velocity,
        gap_j,
    })
}

#[allow(clippy::too_many_arguments)]
fn finish(
    input: &CoupledFrictionImpulseInput<'_>,
    controls: CoupledFrictionImpulseControls,
    scales: &Vector,
    factor: &Factor,
    impulse: Vector,
    after: Vector,
    residual: Residual,
    iterations: usize,
    factorizations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CoupledFrictionImpulseReport, CoupledFrictionImpulseError> {
    let count = input.velocity_before_m_s.len();
    let mut change = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut midpoint = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    for i in 0..count {
        poll(cancelled)?;
        if impulse[i].abs() > controls.maximum_impulse_ns {
            return Err(CoupledFrictionImpulseError::ImpulseBudget {
                index: i,
                impulse_ns: impulse[i],
                maximum_ns: controls.maximum_impulse_ns,
            });
        }
        change[i] = finite(
            after[i] - input.velocity_before_m_s[i],
            "actual velocity change",
        )?;
        midpoint[i] = 0.5 * input.velocity_before_m_s[i] + 0.5 * after[i];
    }
    let reconstructed = factor.solve(&change, scales, cancelled)?;
    let mut impulse_residual_ns = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut maximum_impulse_residual = 0.0_f64;
    let mut work = Dd::ZERO;
    let mut kinetic = Dd::ZERO;
    for i in 0..count {
        poll(cancelled)?;
        impulse_residual_ns[i] = finite(reconstructed[i] - impulse[i], "impulse reconstruction")?;
        maximum_impulse_residual = maximum_impulse_residual.max(impulse_residual_ns[i].abs());
        work = work + Dd::from_f64(impulse[i]) * Dd::from_f64(midpoint[i]);
        kinetic = kinetic + Dd::from_f64(reconstructed[i]) * Dd::from_f64(midpoint[i]);
    }
    let dissipated_j = if work.to_f64() == 0.0 {
        0.0
    } else {
        finite((-work).to_f64(), "midpoint loss")?
    };
    let kinetic_energy_change_j = finite(kinetic.to_f64(), "kinetic change")?;
    let energy_residual_j = finite((kinetic - work).to_f64(), "energy residual")?;
    if dissipated_j < 0.0 || energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(CoupledFrictionImpulseError::EnergyBalance {
            dissipated_j,
            residual_j: energy_residual_j,
        });
    }
    if maximum_impulse_residual > controls.impulse_tolerance_ns {
        return Err(CoupledFrictionImpulseError::Unresolved {
            velocity_residual_m_s: residual.maximum_velocity,
            impulse_residual_ns: maximum_impulse_residual,
            energy_residual_j,
        });
    }
    poll(cancelled)?;
    Ok(CoupledFrictionImpulseReport {
        count,
        impulse_ns: impulse,
        velocity_after_m_s: after,
        impulse_residual_ns,
        velocity_residual_m_s: residual.velocity,
        dissipated_j,
        kinetic_energy_change_j,
        energy_residual_j,
        duality_gap_j: residual.gap_j,
        iterations,
        factorizations,
    })
}

/// Bounded factor storage for this small, physical friction-port owner. The
/// general dense Cholesky/orthant owners allocate; this path retains one full
/// admission factor and one reusable active-face factor on the stack.
struct Factor {
    lower: [[f64; MAX_FRICTION_IMPULSE_PORTS]; MAX_FRICTION_IMPULSE_PORTS],
    indices: [usize; MAX_FRICTION_IMPULSE_PORTS],
    count: usize,
}

impl Factor {
    fn empty() -> Self {
        Self {
            lower: [[0.0; MAX_FRICTION_IMPULSE_PORTS]; MAX_FRICTION_IMPULSE_PORTS],
            indices: [0; MAX_FRICTION_IMPULSE_PORTS],
            count: 0,
        }
    }

    fn prepare(
        &mut self,
        input: &CoupledFrictionImpulseInput<'_>,
        scales: &Vector,
        mask: u32,
        controls: CoupledFrictionImpulseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CoupledFrictionImpulseError> {
        let n = input.velocity_before_m_s.len();
        self.count = 0;
        for i in 0..n {
            if mask & (1_u32 << i) != 0 {
                self.indices[self.count] = i;
                self.count += 1;
            }
        }
        let floor = controls
            .relative_matrix_tolerance
            .max(n as f64 * f64::EPSILON);
        for a in 0..self.count {
            poll(cancelled)?;
            let i = self.indices[a];
            for b in 0..=a {
                let j = self.indices[b];
                let left = normalized(input.mobility_per_kg[i * n + j], scales[i], scales[j])?;
                let right = normalized(input.mobility_per_kg[j * n + i], scales[j], scales[i])?;
                let mut sum = Dd::from_f64(0.5 * left + 0.5 * right);
                for k in 0..b {
                    sum = sum - Dd::from_f64(self.lower[a][k]) * Dd::from_f64(self.lower[b][k]);
                }
                let value = finite(sum.to_f64(), "normalized Cholesky pivot")?;
                self.lower[a][b] = if a == b {
                    if value <= floor {
                        return Err(CoupledFrictionImpulseError::Mobility {
                            index: i,
                            scaled_pivot: value,
                        });
                    }
                    value.sqrt()
                } else {
                    finite(value / self.lower[b][b], "normalized Cholesky entry")?
                };
            }
        }
        Ok(())
    }

    fn solve(
        &self,
        rhs: &Vector,
        scales: &Vector,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vector, CoupledFrictionImpulseError> {
        let mut temporary = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        for a in 0..self.count {
            poll(cancelled)?;
            let i = self.indices[a];
            let mut sum = Dd::from_f64(finite(rhs[i] / scales[i], "scaled right-hand side")?);
            for (b, &value) in temporary.iter().enumerate().take(a) {
                sum = sum - Dd::from_f64(self.lower[a][b]) * Dd::from_f64(value);
            }
            temporary[a] = finite(sum.to_f64() / self.lower[a][a], "forward solve")?;
        }
        for a in (0..self.count).rev() {
            poll(cancelled)?;
            let mut sum = Dd::from_f64(temporary[a]);
            for (b, &value) in temporary.iter().enumerate().take(self.count).skip(a + 1) {
                sum = sum - Dd::from_f64(self.lower[b][a]) * Dd::from_f64(value);
            }
            temporary[a] = finite(sum.to_f64() / self.lower[a][a], "backward solve")?;
        }
        let mut result = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        for (a, &value) in temporary.iter().enumerate().take(self.count) {
            let i = self.indices[a];
            result[i] = finite(value / scales[i], "physical impulse solve")?;
        }
        Ok(result)
    }
}

fn response(
    input: &CoupledFrictionImpulseInput<'_>,
    impulse: &Vector,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vector, CoupledFrictionImpulseError> {
    let n = input.velocity_before_m_s.len();
    let mut after = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    for (i, value) in after.iter_mut().enumerate().take(n) {
        poll(cancelled)?;
        let mut sum = Dd::from_f64(input.velocity_before_m_s[i]);
        for (j, &p) in impulse.iter().enumerate().take(n) {
            sum = sum + Dd::from_f64(input.mobility_per_kg[i * n + j]) * Dd::from_f64(p);
        }
        *value = finite(sum.to_f64(), "original endpoint velocity")?;
    }
    Ok(after)
}

fn normalized(value: f64, first: f64, second: f64) -> Result<f64, CoupledFrictionImpulseError> {
    finite((value / first) / second, "normalized mobility")
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn finite(value: f64, label: &'static str) -> Result<f64, CoupledFrictionImpulseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CoupledFrictionImpulseError::NonFinite(label))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CoupledFrictionImpulseError> {
    if cancelled() {
        Err(CoupledFrictionImpulseError::Cancelled)
    } else {
        Ok(())
    }
}

fn unresolved(residual: &Residual) -> CoupledFrictionImpulseError {
    CoupledFrictionImpulseError::Unresolved {
        velocity_residual_m_s: residual.maximum_velocity,
        impulse_residual_ns: 0.0,
        energy_residual_j: residual.gap_j,
    }
}
