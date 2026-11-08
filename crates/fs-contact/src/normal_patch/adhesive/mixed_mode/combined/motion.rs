//! Fully coupled normal, two-shear and torsional inertia on one adhesive disk.
//!
//! q=[approach,shear1,shear2,twist] has units [m,m,m,rad]; conjugate efforts
//! have units [N,N,N,Nm]. The original full-stick I/II/III source relaxes one
//! shared contact radius, so every coordinate changes the same potential.
//! q''=a_external-W*f retains the complete reciprocal native 4x4 mobility.
//!
//! W=D*C*D with D_ii=sqrt(W_ii) gives a common energy metric without comparing
//! metres to radians. Fixed 4x4 factors solve the analytic discrete-gradient
//! equation. A bounded damped Newton continuation admits only a resolved
//! positive symmetric inertial tangent. This is local numerical admission,
//! not uniqueness or a certificate of arbitrary-step accuracy.
//!
//! The actual held-effort trajectory is quadratic. Its coordinate extrema
//! enter the original joint source box, including hidden shear/twist turns;
//! independent endpoint or mode-by-mode clipping cannot admit the phase.
//! All scratch and factors are fixed-size arrays. No allocation, restarted
//! source history, Coulomb law, irreversible loss or heat is introduced.

use super::{
    JkrCombinedResponse, SourceJkrMixedMode,
    path::{JkrCombinedPath, JkrCombinedPathError},
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Component order is [normal,shear1,shear2,twist]. Mixed-unit arrays keep
/// metre and radian tolerances, or force and torque limits, separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedMotionBudget {
    pub maximum_step_s: f64,
    /// Maximum whole-path excursion from the initial coordinate [m,m,m,rad].
    pub maximum_coordinate_increment: [f64; 4],
    pub maximum_rate: [f64; 4],
    /// Original source, discrete and equivalent external effort bounds [N,N,N,Nm].
    pub maximum_effort: [f64; 4],
    pub maximum_iterations: usize,
    /// Actual original-source queries, including rejected probes and corners.
    pub maximum_evaluations: usize,
    pub maximum_backtracks: usize,
    pub coordinate_absolute_tolerance: [f64; 4],
    pub coordinate_relative_tolerance: f64,
    pub rate_absolute_tolerance: [f64; 4],
    pub rate_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrCombinedMotionError {
    Path(JkrCombinedPathError),
    Invalid(&'static str),
    Range(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    /// Refine the physical phase; no alternative inertial root is selected.
    InertialTangent,
    NotConverged {
        residual: [f64; 4],
        tolerance: [f64; 4],
    },
    RateBalance {
        residual: [f64; 4],
        tolerance: [f64; 4],
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<JkrCombinedPathError> for JkrCombinedMotionError {
    fn from(value: JkrCombinedPathError) -> Self {
        match value {
            JkrCombinedPathError::Cancelled => Self::Cancelled,
            JkrCombinedPathError::Budget => Self::Budget,
            other => Self::Path(other),
        }
    }
}
impl core::fmt::Display for JkrCombinedMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled normal/shear/twist JKR motion: {self:?}")
    }
}
impl std::error::Error for JkrCombinedMotionError {}

/// The original endpoint is cached privately. Neither public response fields
/// nor an unrelated source, mobility or numerical policy can replace it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedMotionState {
    identity: ContentHash,
    contact: JkrCombinedResponse,
    rates: [f64; 4],
    time_s: f64,
    steps: u64,
}
impl JkrCombinedMotionState {
    pub const fn coordinates(self) -> [f64; 4] {
        self.contact.coordinates
    }
    pub const fn rates(self) -> [f64; 4] {
        self.rates
    }
    pub const fn contact(self) -> JkrCombinedResponse {
        self.contact
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn steps(self) -> u64 {
        self.steps
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedMotionStep {
    initial: JkrCombinedMotionState,
    state: JkrCombinedMotionState,
    pub initial_contact: JkrCombinedResponse,
    pub contact: JkrCombinedResponse,
    pub mean_effort: [f64; 4],
    /// Native reactions are -J^T times these impulses [Ns,Ns,Ns,Nms].
    pub effort_impulse: [f64; 4],
    /// First time moment of the held efforts [Ns²,Ns²,Ns²,Nms²].
    pub effort_time_moment: [f64; 4],
    pub contact_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub applied_work_j: f64,
    pub energy_residual_j: f64,
    pub source_work_residual_j: f64,
    pub coordinate_residual: [f64; 4],
    pub rate_residual: [f64; 4],
    pub maximum_source_effort: [f64; 4],
    pub iterations: usize,
    pub evaluations: usize,
}
impl JkrCombinedMotionStep {
    pub const fn initial_state(&self) -> JkrCombinedMotionState {
        self.initial
    }
    pub const fn state(&self) -> JkrCombinedMotionState {
        self.state
    }
}

/// Immutable source and held projection; construction and repeated proposals
/// allocate nothing. Accepted history is never used as Newton scratch.
#[derive(Clone, Debug)]
pub struct JkrCombinedMotion<'a> {
    source: &'a SourceJkrMixedMode,
    mobility: [[f64; 4]; 4],
    scale: [f64; 4],
    inverse_correlation: [[f64; 4]; 4],
    inverse_factor: [[f64; 4]; 4],
    budget: JkrCombinedMotionBudget,
    identity: ContentHash,
}

struct Trial {
    position: [f64; 4],
    path: JkrCombinedPath,
    residual: [f64; 4],
    norm: f64,
    normalized_rhs: [f64; 4],
    jacobian: [[f64; 4]; 4],
}

impl<'a> JkrCombinedMotion<'a> {
    pub fn new(
        source: &'a SourceJkrMixedMode,
        mut mobility: [[f64; 4]; 4],
        budget: JkrCombinedMotionBudget,
    ) -> Result<Self, JkrCombinedMotionError> {
        validate_budget(budget)?;
        if mobility.iter().flatten().any(|x| !x.is_finite())
            || (0..4).any(|i| mobility[i][i] <= 0.0)
        {
            return Err(JkrCombinedMotionError::Invalid(
                "finite positive mobility diagonals required",
            ));
        }
        let scale: [f64; 4] = core::array::from_fn(|i| mobility[i][i].sqrt());
        let mut correlation = [[0.0; 4]; 4];
        for i in 0..4 {
            correlation[i][i] = (mobility[i][i] / scale[i]) / scale[i];
            for j in 0..i {
                let left = (mobility[i][j] / scale[i]) / scale[j];
                let right = (mobility[j][i] / scale[j]) / scale[i];
                if !left.is_finite()
                    || !right.is_finite()
                    || (left - right).abs() > 64.0 * f64::EPSILON
                {
                    return Err(JkrCombinedMotionError::Invalid(
                        "reciprocal mixed-unit mobility required",
                    ));
                }
                // Canonicalize only admitted representation-level reciprocity
                // differences; the full physical matrix enters the identity.
                let cross = mobility[i][j].midpoint(mobility[j][i]);
                mobility[i][j] = cross;
                mobility[j][i] = cross;
                let value = (cross / scale[i]) / scale[j];
                correlation[i][j] = value;
                correlation[j][i] = value;
            }
        }
        let factor = cholesky(correlation).ok_or(JkrCombinedMotionError::Invalid(
            "resolved positive native mobility required",
        ))?;
        let inverse_factor = inverse_lower(factor).ok_or(JkrCombinedMotionError::Invalid(
            "unrepresentable native mass factor",
        ))?;
        let mut inverse_correlation = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..=i {
                let left = core::array::from_fn(|k| inverse_factor[k][i]);
                let right = core::array::from_fn(|k| inverse_factor[k][j]);
                let value = dot(left, right);
                if !value.is_finite() {
                    return Err(JkrCombinedMotionError::Invalid(
                        "unrepresentable native mass metric",
                    ));
                }
                inverse_correlation[i][j] = value;
                inverse_correlation[j][i] = value;
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-combined-motion.v1");
        hash.update(source.identity().as_bytes());
        for x in mobility
            .into_iter()
            .flatten()
            .chain(budget.maximum_coordinate_increment)
            .chain(budget.maximum_rate)
            .chain(budget.maximum_effort)
            .chain(budget.coordinate_absolute_tolerance)
            .chain(budget.rate_absolute_tolerance)
            .chain([
                budget.maximum_step_s,
                budget.coordinate_relative_tolerance,
                budget.rate_relative_tolerance,
                budget.energy_absolute_tolerance_j,
                budget.energy_relative_tolerance,
            ])
        {
            hash.update(&x.to_bits().to_le_bytes());
        }
        for x in [
            budget.maximum_iterations,
            budget.maximum_evaluations,
            budget.maximum_backtracks,
        ] {
            hash.update(&(x as u64).to_le_bytes());
        }
        Ok(Self {
            source,
            mobility,
            scale,
            inverse_correlation,
            inverse_factor,
            budget,
            identity: hash.finalize(),
        })
    }
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub const fn mobility(&self) -> [[f64; 4]; 4] {
        self.mobility
    }
    pub const fn budget(&self) -> JkrCombinedMotionBudget {
        self.budget
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Begin with explicitly authored attached coordinates; no bond nucleation
    /// or hidden twist preload is inferred. Exactly one original query is used.
    pub fn initial_state(
        &self,
        coordinates: [f64; 4],
        rates: [f64; 4],
        gate: &CancelGate,
    ) -> Result<JkrCombinedMotionState, JkrCombinedMotionError> {
        poll(gate)?;
        bounded(rates, self.budget.maximum_rate, "initial rates")?;
        let path = self
            .source
            .combined_discrete_path(coordinates, coordinates, gate)?;
        bounded(
            path.mean_effort,
            self.budget.maximum_effort,
            "initial source efforts",
        )?;
        self.kinetic_energy(rates)?;
        poll(gate)?;
        Ok(JkrCombinedMotionState {
            identity: self.identity,
            contact: path.contact,
            rates,
            time_s: 0.0,
            steps: 0,
        })
    }

    /// A held acceleration is the projection J M^-1 f_external from the actual
    /// mechanical owner. Every original normal/shear/twist cross response remains.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        previous: JkrCombinedMotionState,
        dt_s: f64,
        external_acceleration: [f64; 4],
        gate: &CancelGate,
    ) -> Result<JkrCombinedMotionStep, JkrCombinedMotionError> {
        poll(gate)?;
        self.check(previous)?;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.budget.maximum_step_s {
            return Err(JkrCombinedMotionError::Range("phase duration"));
        }
        let time_s = finite(previous.time_s + dt_s, "motion clock")?;
        if time_s <= previous.time_s {
            return Err(JkrCombinedMotionError::Range(
                "unrepresentable clock increment",
            ));
        }
        let steps = previous
            .steps
            .checked_add(1)
            .ok_or(JkrCombinedMotionError::Budget)?;
        let half_h2 = finite(0.5 * dt_s * dt_s, "squared duration")?;
        if half_h2 == 0.0 {
            return Err(JkrCombinedMotionError::Range("squared duration underflow"));
        }
        let normalized_external = checked_array(
            core::array::from_fn(|i| external_acceleration[i] / self.scale[i]),
            "external acceleration",
        )?;
        let force_scaled = mul(self.inverse_correlation, normalized_external);
        let applied_effort = checked_array(
            core::array::from_fn(|i| force_scaled[i] / self.scale[i]),
            "equivalent external effort",
        )?;
        bounded(
            applied_effort,
            self.budget.maximum_effort,
            "equivalent external effort",
        )?;
        let free = checked_array(
            core::array::from_fn(|i| {
                (Dd::from_f64(dt_s) * Dd::from_f64(previous.rates[i])
                    + Dd::from_f64(half_h2) * Dd::from_f64(external_acceleration[i]))
                .to_f64()
            }),
            "free coordinate increment",
        )?;
        let mut evaluations = 0;
        let mut trial = self.trial(
            previous,
            previous.coordinates(),
            free,
            half_h2,
            &mut evaluations,
            gate,
        )?;
        for iteration in 0..=self.budget.maximum_iterations {
            poll(gate)?;
            let tolerance = self.coordinate_tolerance(previous.coordinates(), trial.position)?;
            if within(trial.residual, tolerance) {
                return self.finish(
                    previous,
                    trial,
                    applied_effort,
                    external_acceleration,
                    dt_s,
                    time_s,
                    steps,
                    iteration,
                    evaluations,
                    gate,
                );
            }
            if iteration == self.budget.maximum_iterations {
                return Err(JkrCombinedMotionError::NotConverged {
                    residual: trial.residual,
                    tolerance,
                });
            }
            let normalized_correction = solve(trial.jacobian, trial.normalized_rhs.map(|x| -x))?;
            let correction = checked_array(
                core::array::from_fn(|i| self.scale[i] * normalized_correction[i]),
                "Newton correction",
            )?;
            let mut next = None;
            let mut fraction = 1.0;
            let mut last_error = None;
            for _ in 0..=self.budget.maximum_backtracks {
                poll(gate)?;
                let position =
                    core::array::from_fn(|i| trial.position[i] + fraction * correction[i]);
                match self.trial(previous, position, free, half_h2, &mut evaluations, gate) {
                    Ok(candidate) => {
                        if candidate.norm < trial.norm
                            || within(
                                candidate.residual,
                                self.coordinate_tolerance(previous.coordinates(), position)?,
                            )
                        {
                            next = Some(candidate);
                            break;
                        }
                    }
                    Err(JkrCombinedMotionError::Budget) => {
                        return Err(JkrCombinedMotionError::Budget);
                    }
                    Err(JkrCombinedMotionError::Cancelled) => {
                        return Err(JkrCombinedMotionError::Cancelled);
                    }
                    Err(error) => last_error = Some(error),
                }
                fraction *= 0.5;
            }
            trial = next.ok_or_else(|| {
                last_error.unwrap_or(JkrCombinedMotionError::NotConverged {
                    residual: trial.residual,
                    tolerance,
                })
            })?;
        }
        unreachable!("bounded Newton loop always returns")
    }

    fn trial(
        &self,
        previous: JkrCombinedMotionState,
        to: [f64; 4],
        free: [f64; 4],
        half_h2: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Trial, JkrCombinedMotionError> {
        poll(gate)?;
        let increment = sub(to, previous.coordinates());
        bounded(
            increment,
            self.budget.maximum_coordinate_increment,
            "coordinate increment",
        )?;
        let path = self.source.combined_path_counted(
            previous.contact,
            to,
            evaluations,
            self.budget.maximum_evaluations,
            gate,
        )?;
        bounded(
            path.mean_effort,
            self.budget.maximum_effort,
            "discrete source efforts",
        )?;
        let reaction = mul(self.mobility, path.mean_effort);
        let residual = checked_array(
            core::array::from_fn(|i| {
                (Dd::from_f64(increment[i]) - Dd::from_f64(free[i])
                    + Dd::from_f64(half_h2) * Dd::from_f64(reaction[i]))
                .to_f64()
            }),
            "coordinate residual",
        )?;
        let norm = finite(
            residual
                .into_iter()
                .zip(self.budget.maximum_coordinate_increment)
                .map(|(value, limit)| (value / limit).abs())
                .fold(0.0_f64, f64::max),
            "normalized residual",
        )?;
        let normalized_rhs = checked_array(
            mul(
                self.inverse_correlation,
                core::array::from_fn(|i| residual[i] / self.scale[i]),
            ),
            "normalized inertial residual",
        )?;
        let jacobian = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                self.inverse_correlation[i][j]
                    + half_h2 * self.scale[i] * path.end_tangent[i][j] * self.scale[j]
            })
        });
        if !coercive(jacobian) {
            return Err(JkrCombinedMotionError::InertialTangent);
        }
        Ok(Trial {
            position: to,
            path,
            residual,
            norm,
            normalized_rhs,
            jacobian,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        previous: JkrCombinedMotionState,
        trial: Trial,
        applied: [f64; 4],
        external: [f64; 4],
        dt: f64,
        time_s: f64,
        steps: u64,
        iterations: usize,
        mut evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrCombinedMotionStep, JkrCombinedMotionError> {
        let reaction = mul(self.mobility, trial.path.mean_effort);
        let expected_change = checked_array(
            core::array::from_fn(|i| dt * (external[i] - reaction[i])),
            "rate change",
        )?;
        let rates = checked_array(
            core::array::from_fn(|i| previous.rates[i] + expected_change[i]),
            "endpoint rates",
        )?;
        bounded(rates, self.budget.maximum_rate, "endpoint rates")?;
        self.kinetic_energy(rates)?;
        let midpoint: [f64; 4] = core::array::from_fn(|i| previous.rates[i].midpoint(rates[i]));
        let increment = sub(trial.position, previous.coordinates());
        let coordinate_residual = checked_array(
            core::array::from_fn(|i| increment[i] - dt * midpoint[i]),
            "midpoint coordinate balance",
        )?;
        let coordinate_tolerance =
            self.coordinate_tolerance(previous.coordinates(), trial.position)?;
        if !within(coordinate_residual, coordinate_tolerance) {
            return Err(JkrCombinedMotionError::NotConverged {
                residual: coordinate_residual,
                tolerance: coordinate_tolerance,
            });
        }
        let rate_change = sub(rates, previous.rates);
        let rate_residual = sub(rate_change, expected_change);
        let rate_tolerance = checked_array(
            core::array::from_fn(|i| {
                self.budget.rate_absolute_tolerance[i]
                    + self.budget.rate_relative_tolerance
                        * previous.rates[i].abs().max(rates[i].abs())
            }),
            "rate tolerance",
        )?;
        if !within(rate_residual, rate_tolerance) {
            return Err(JkrCombinedMotionError::RateBalance {
                residual: rate_residual,
                tolerance: rate_tolerance,
            });
        }
        let control = checked_array(
            core::array::from_fn(|i| previous.coordinates()[i] + 0.5 * dt * previous.rates[i]),
            "quadratic control point",
        )?;
        let physical_endpoint = checked_array(
            core::array::from_fn(|i| {
                (Dd::from_f64(previous.coordinates()[i])
                    + Dd::from_f64(dt) * Dd::from_f64(midpoint[i]))
                .to_f64()
            }),
            "quadratic endpoint",
        )?;
        // The native impulse/time moment produces this constant-acceleration
        // trajectory. Enclose both its endpoint and the original-source
        // Newton endpoint, whose admitted coordinate residual is explicit.
        let mut bounds = quadratic_bounds([previous.coordinates(), control, physical_endpoint])?;
        for i in 0..4 {
            bounds[0][i] = bounds[0][i].min(trial.position[i]);
            bounds[1][i] = bounds[1][i].max(trial.position[i]);
        }
        for corner in bounds {
            bounded(
                sub(corner, previous.coordinates()),
                self.budget.maximum_coordinate_increment,
                "quadratic coordinate excursion",
            )?;
        }
        let from = previous.coordinates();
        let covered = (0..4).all(|i| {
            bounds[0][i] >= from[i].min(trial.position[i])
                && bounds[1][i] <= from[i].max(trial.position[i])
        });
        // A monotone phase already inside the admitted endpoint box reuses
        // its original radius bounds; no duplicate constitutive solves.
        let radii = if covered {
            [
                trial.path.minimum_contact_radius_m,
                trial.path.maximum_contact_radius_m,
            ]
        } else {
            self.source.admit_combined_box(
                previous.contact,
                trial.path.contact,
                bounds,
                &mut evaluations,
                self.budget.maximum_evaluations,
                gate,
            )?
        };
        let maximum_source_effort = self.source_effort_bounds(bounds, radii)?;
        bounded(
            maximum_source_effort,
            self.budget.maximum_effort,
            "whole-path source efforts",
        )?;
        let scaled_midpoint = checked_array(
            core::array::from_fn(|i| midpoint[i] / self.scale[i]),
            "normalized midpoint rate",
        )?;
        let scaled_change = checked_array(
            core::array::from_fn(|i| rate_change[i] / self.scale[i]),
            "normalized rate change",
        )?;
        let kinetic_change = finite(
            dot(
                mul(self.inverse_factor, scaled_midpoint),
                mul(self.inverse_factor, scaled_change),
            ),
            "kinetic energy change",
        )?;
        let applied_work = finite(dot(applied, increment), "applied work")?;
        let energy_residual = finite(
            (Dd::from_f64(kinetic_change) + Dd::from_f64(trial.path.reversible_work_j)
                - Dd::from_f64(applied_work))
            .to_f64(),
            "energy residual",
        )?;
        let energy_tolerance = finite(
            self.budget.energy_absolute_tolerance_j
                + self.budget.energy_relative_tolerance
                    * (kinetic_change.abs()
                        + trial.path.reversible_work_j.abs()
                        + applied_work.abs()),
            "energy tolerance",
        )?;
        if energy_residual.abs().max(trial.path.work_residual_j.abs()) > energy_tolerance {
            return Err(JkrCombinedMotionError::EnergyBalance {
                residual_j: energy_residual.abs().max(trial.path.work_residual_j.abs()),
                tolerance_j: energy_tolerance,
            });
        }
        let effort_impulse =
            checked_array(trial.path.mean_effort.map(|x| x * dt), "contact impulse")?;
        let effort_time_moment = checked_array(
            effort_impulse.map(|x| x * (0.5 * dt)),
            "contact time moment",
        )?;
        poll(gate)?;
        Ok(JkrCombinedMotionStep {
            initial: previous,
            state: JkrCombinedMotionState {
                identity: self.identity,
                contact: trial.path.contact,
                rates,
                time_s,
                steps,
            },
            initial_contact: previous.contact,
            contact: trial.path.contact,
            mean_effort: trial.path.mean_effort,
            effort_impulse,
            effort_time_moment,
            contact_work_j: trial.path.reversible_work_j,
            kinetic_energy_change_j: kinetic_change,
            applied_work_j: applied_work,
            energy_residual_j: energy_residual,
            source_work_residual_j: trial.path.work_residual_j,
            coordinate_residual,
            rate_residual,
            maximum_source_effort,
            iterations,
            evaluations,
        })
    }

    fn source_effort_bounds(
        &self,
        q: [[f64; 4]; 2],
        radii: [f64; 2],
    ) -> Result<[f64; 4], JkrCombinedMotionError> {
        let normal = self.source.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let r = normal.effective_radius_m();
        let mut force = 0.0_f64;
        // N=2Ea(delta-a²/(3R)) is affine in delta. Its only positive-radius
        // stationary point at fixed delta is sqrt(R*delta).
        for delta in [q[0][0], q[1][0]] {
            for a in radii {
                let value = finite(
                    2.0 * e * a * (delta - a * (a / r) / 3.0),
                    "source normal effort bound",
                )?;
                force = force.max(value.abs());
            }
            if delta > 0.0 {
                let a = r.sqrt() * delta.sqrt();
                if a >= radii[0] && a <= radii[1] {
                    force = force.max(
                        finite((4.0 / 3.0) * e * a * delta, "source normal effort extremum")?.abs(),
                    );
                }
            }
        }
        let a = radii[1];
        let shear_stiffness = 8.0 * self.source.reduced_shear_modulus_pa() * a;
        let twist = q[0][3].abs().max(q[1][3].abs());
        let moment =
            ((16.0 / 3.0) * self.source.torsional_reduced_shear_modulus_pa() * a) * a * (a * twist);
        let bound = [
            force,
            shear_stiffness * q[0][1].abs().max(q[1][1].abs()),
            shear_stiffness * q[0][2].abs().max(q[1][2].abs()),
            moment,
        ];
        checked_array(
            bound.map(|x| if x == 0.0 { 0.0 } else { x.next_up() }),
            "original combined source effort enclosure",
        )
    }

    fn coordinate_tolerance(
        &self,
        from: [f64; 4],
        to: [f64; 4],
    ) -> Result<[f64; 4], JkrCombinedMotionError> {
        checked_array(
            core::array::from_fn(|i| {
                self.budget.coordinate_absolute_tolerance[i]
                    + self.budget.coordinate_relative_tolerance * from[i].abs().max(to[i].abs())
            }),
            "coordinate tolerance",
        )
    }
    fn check(&self, state: JkrCombinedMotionState) -> Result<(), JkrCombinedMotionError> {
        if state.identity != self.identity {
            return Err(JkrCombinedMotionError::StateMismatch);
        }
        Ok(())
    }
    fn kinetic_energy(&self, rates: [f64; 4]) -> Result<f64, JkrCombinedMotionError> {
        let scaled = checked_array(
            core::array::from_fn(|i| rates[i] / self.scale[i]),
            "normalized rates",
        )?;
        let physical = checked_array(mul(self.inverse_factor, scaled), "native kinetic factor")?;
        let energy = finite(0.5 * dot(physical, physical), "projected kinetic energy")?;
        if energy < 0.0 {
            return Err(JkrCombinedMotionError::Range(
                "negative projected kinetic energy",
            ));
        }
        Ok(energy)
    }
}

fn validate_budget(b: JkrCombinedMotionBudget) -> Result<(), JkrCombinedMotionError> {
    if [b.maximum_step_s]
        .into_iter()
        .chain(b.maximum_coordinate_increment)
        .chain(b.maximum_rate)
        .chain(b.maximum_effort)
        .any(|x| !x.is_finite() || x <= 0.0)
        || b.coordinate_absolute_tolerance
            .into_iter()
            .chain(b.rate_absolute_tolerance)
            .chain([
                b.coordinate_relative_tolerance,
                b.rate_relative_tolerance,
                b.energy_absolute_tolerance_j,
                b.energy_relative_tolerance,
            ])
            .any(|x| !x.is_finite() || x < 0.0)
        || b.maximum_iterations == 0
        || b.maximum_iterations > 256
        || b.maximum_evaluations == 0
        || b.maximum_evaluations > 1_048_576
        || b.maximum_backtracks > 64
    {
        Err(JkrCombinedMotionError::Invalid(
            "finite physical bounds, nonnegative tolerances and bounded work required",
        ))
    } else {
        Ok(())
    }
}

/// Extrema of each actual quadratic coordinate, evaluated by de Casteljau.
/// Outward ULP padding screens floating-point roundoff; no interval certificate
/// is claimed. Correlations between coordinates may make the box conservative.
fn quadratic_bounds(points: [[f64; 4]; 3]) -> Result<[[f64; 4]; 2], JkrCombinedMotionError> {
    let mut bounds = [
        core::array::from_fn(|i| points[0][i].min(points[2][i])),
        core::array::from_fn(|i| points[0][i].max(points[2][i])),
    ];
    for i in 0..4 {
        let left = finite(points[1][i] - points[0][i], "quadratic initial slope")?;
        let right = finite(points[2][i] - points[1][i], "quadratic final slope")?;
        if left != 0.0 && right != 0.0 && left.is_sign_negative() != right.is_sign_negative() {
            let scale = left.abs().max(right.abs());
            let t = (left / scale) / ((left / scale) - (right / scale));
            let a = points[0][i] + t * left;
            let b = points[1][i] + t * right;
            let extremum = finite(a + t * (b - a), "quadratic extremum")?;
            bounds[0][i] = bounds[0][i].min(extremum.next_down());
            bounds[1][i] = bounds[1][i].max(extremum.next_up());
        }
    }
    Ok(bounds)
}

/// Fixed unpivoted Cholesky for a symmetric positive metric. Each pivot must
/// resolve against its original diagonal; no unit-dependent absolute cutoff.
fn cholesky(a: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    if a.iter().flatten().any(|x| !x.is_finite()) || (0..4).any(|i| a[i][i] <= 0.0) {
        return None;
    }
    let mut lower = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            let mut residual = Dd::from_f64(a[i][j]);
            for k in 0..j {
                residual = residual - Dd::from_f64(lower[i][k]) * Dd::from_f64(lower[j][k]);
            }
            let value = residual.to_f64();
            if !value.is_finite() {
                return None;
            }
            if i == j {
                if value <= 64.0 * f64::EPSILON * a[i][i] {
                    return None;
                }
                lower[i][j] = value.sqrt();
            } else {
                lower[i][j] = value / lower[j][j];
                if !lower[i][j].is_finite() {
                    return None;
                }
            }
        }
    }
    Some(lower)
}

fn inverse_lower(lower: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut inverse = [[0.0; 4]; 4];
    for column in 0..4 {
        for i in 0..4 {
            let mut value = Dd::from_f64(if i == column { 1.0 } else { 0.0 });
            for k in 0..i {
                value = value - Dd::from_f64(lower[i][k]) * Dd::from_f64(inverse[k][column]);
            }
            inverse[i][column] = value.to_f64() / lower[i][i];
            if !inverse[i][column].is_finite() {
                return None;
            }
        }
    }
    Some(inverse)
}

/// The discrete endpoint derivative is generally nonsymmetric. Its symmetric
/// part must be resolved SPD in its own diagonal scale before Newton proceeds.
fn coercive(a: [[f64; 4]; 4]) -> bool {
    if a.iter().flatten().any(|x| !x.is_finite()) || (0..4).any(|i| a[i][i] <= 0.0) {
        return false;
    }
    let scale: [f64; 4] = core::array::from_fn(|i| a[i][i].sqrt());
    let mut symmetric = [[0.0; 4]; 4];
    for i in 0..4 {
        symmetric[i][i] = 1.0;
        for j in 0..i {
            let cross = (a[i][j].midpoint(a[j][i]) / scale[i]) / scale[j];
            symmetric[i][j] = cross;
            symmetric[j][i] = cross;
        }
    }
    cholesky(symmetric).is_some()
}

/// Diagonal equilibration and deterministic scaled partial pivoting retain the
/// actual nonsymmetric analytic tangent. All four rows live on the stack.
fn solve(a: [[f64; 4]; 4], rhs: [f64; 4]) -> Result<[f64; 4], JkrCombinedMotionError> {
    let scale: [f64; 4] = core::array::from_fn(|i| a[i][i].sqrt());
    let mut matrix: [[f64; 4]; 4] =
        core::array::from_fn(|i| core::array::from_fn(|j| (a[i][j] / scale[i]) / scale[j]));
    let mut value: [f64; 4] = core::array::from_fn(|i| rhs[i] / scale[i]);
    let mut row_scale: [f64; 4] = matrix.map(|row| {
        row.into_iter()
            .fold(0.0_f64, |maximum, x| maximum.max(x.abs()))
    });
    if matrix
        .iter()
        .flatten()
        .chain(value.iter())
        .any(|x| !x.is_finite())
        || row_scale.iter().any(|x| *x <= 0.0)
    {
        return Err(JkrCombinedMotionError::InertialTangent);
    }
    for k in 0..4 {
        let mut pivot = k;
        for i in k + 1..4 {
            if matrix[i][k].abs() / row_scale[i] > matrix[pivot][k].abs() / row_scale[pivot] {
                pivot = i;
            }
        }
        if matrix[pivot][k].abs() <= 64.0 * f64::EPSILON * row_scale[pivot] {
            return Err(JkrCombinedMotionError::InertialTangent);
        }
        matrix.swap(k, pivot);
        value.swap(k, pivot);
        row_scale.swap(k, pivot);
        for i in k + 1..4 {
            let factor = matrix[i][k] / matrix[k][k];
            matrix[i][k] = 0.0;
            for j in k + 1..4 {
                matrix[i][j] = (Dd::from_f64(matrix[i][j])
                    - Dd::from_f64(factor) * Dd::from_f64(matrix[k][j]))
                .to_f64();
            }
            value[i] =
                (Dd::from_f64(value[i]) - Dd::from_f64(factor) * Dd::from_f64(value[k])).to_f64();
        }
    }
    let mut solution = [0.0; 4];
    for i in (0..4).rev() {
        let mut residual = Dd::from_f64(value[i]);
        for j in i + 1..4 {
            residual = residual - Dd::from_f64(matrix[i][j]) * Dd::from_f64(solution[j]);
        }
        solution[i] = residual.to_f64() / matrix[i][i];
    }
    checked_array(
        core::array::from_fn(|i| solution[i] / scale[i]),
        "normalized Newton solution",
    )
}

fn dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    let mut sum = Dd::from_f64(0.0);
    for i in 0..4 {
        sum = sum + Dd::from_f64(a[i]) * Dd::from_f64(b[i]);
    }
    sum.to_f64()
}
fn mul(a: [[f64; 4]; 4], x: [f64; 4]) -> [f64; 4] {
    a.map(|row| dot(row, x))
}
fn sub(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn within(value: [f64; 4], limit: [f64; 4]) -> bool {
    (0..4).all(|i| value[i].is_finite() && value[i].abs() <= limit[i])
}
fn bounded(
    value: [f64; 4],
    limit: [f64; 4],
    field: &'static str,
) -> Result<(), JkrCombinedMotionError> {
    if within(value, limit) {
        Ok(())
    } else {
        Err(JkrCombinedMotionError::Range(field))
    }
}
fn checked_array(value: [f64; 4], field: &'static str) -> Result<[f64; 4], JkrCombinedMotionError> {
    if value.iter().all(|x| x.is_finite()) {
        Ok(value)
    } else {
        Err(JkrCombinedMotionError::Range(field))
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrCombinedMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrCombinedMotionError::Range(field))
    }
}
fn poll(gate: &CancelGate) -> Result<(), JkrCombinedMotionError> {
    if gate.is_requested() {
        Err(JkrCombinedMotionError::Cancelled)
    } else {
        Ok(())
    }
}
