//! Coupled approach/twist inertia on the original relaxed-area I/III potential.
//!
//! For q=[approach (m),twist (rad)] and conjugate f=[force (N),torque (Nm)],
//! q''=a_external-W f with the complete reciprocal native mobility. Writing
//! W=D C D, D=diag(sqrt(W_ii)), makes both normalized coordinates have the same
//! units. The midpoint/discrete-gradient equation is solved in that metric:
//! C^-1 D^-1(dq-h v0-h² a_external/2) + h² D G/2 = 0.
//! Its analytic Jacobian uses the original relaxed-area path tangent. A bounded
//! damped Newton continuation requires a positive symmetric inertial tangent;
//! this is a local numerical admission, not a global uniqueness certificate.
//!
//! The constant discrete effort defines a quadratic position path. Its actual
//! coordinate extrema, including hidden reversals, feed the original contact
//! box admission and force/torque bounds. Endpoint attachment alone is not
//! sufficient. A refusal publishes no successor. There is no inferred snap,
//! frictional slip, irreversible loss or heat. Held projection, fixed material,
//! small deformation and temporal refinement remain physical assumptions.

use super::{
    JkrTorsionalResponse, SourceJkrMixedMode,
    path::{JkrTorsionalPath, JkrTorsionalPathError},
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Component order is always [normal, torsional]; arrays retain their different
/// SI units rather than combining metres with radians or newtons with Nm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalMotionBudget {
    pub maximum_step_s: f64,
    /// Maximum whole-path displacement from the initial coordinate [m,rad].
    pub maximum_coordinate_increment: [f64; 2],
    pub maximum_rate: [f64; 2],
    /// Bounds original source, discrete and equivalent external efforts [N,Nm].
    pub maximum_effort: [f64; 2],
    pub maximum_iterations: usize,
    /// Actual original-source queries, including rejected probes and corners.
    pub maximum_evaluations: usize,
    pub maximum_backtracks: usize,
    pub coordinate_absolute_tolerance: [f64; 2],
    pub coordinate_relative_tolerance: f64,
    pub rate_absolute_tolerance: [f64; 2],
    pub rate_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrTorsionalMotionError {
    Path(JkrTorsionalPathError),
    Invalid(&'static str),
    Range(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    /// Refine the physical phase; no alternative inertial root is selected.
    InertialTangent,
    NotConverged {
        residual: [f64; 2],
        tolerance: [f64; 2],
    },
    RateBalance {
        residual: [f64; 2],
        tolerance: [f64; 2],
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<JkrTorsionalPathError> for JkrTorsionalMotionError {
    fn from(value: JkrTorsionalPathError) -> Self {
        match value {
            JkrTorsionalPathError::Cancelled => Self::Cancelled,
            JkrTorsionalPathError::Budget => Self::Budget,
            other => Self::Path(other),
        }
    }
}
impl core::fmt::Display for JkrTorsionalMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled normal/twist JKR motion: {self:?}")
    }
}
impl std::error::Error for JkrTorsionalMotionError {}

/// The original endpoint is cached privately. Neither public response fields
/// nor an unrelated source, mobility or numerical policy can replace it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalMotionState {
    identity: ContentHash,
    contact: JkrTorsionalResponse,
    rates: [f64; 2],
    time_s: f64,
    steps: u64,
}
impl JkrTorsionalMotionState {
    pub const fn coordinates(self) -> [f64; 2] {
        [self.contact.indentation_m, self.contact.twist_angle_rad]
    }
    pub const fn rates(self) -> [f64; 2] {
        self.rates
    }
    pub const fn contact(self) -> JkrTorsionalResponse {
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
pub struct JkrTorsionalMotionStep {
    initial: JkrTorsionalMotionState,
    state: JkrTorsionalMotionState,
    pub initial_contact: JkrTorsionalResponse,
    pub contact: JkrTorsionalResponse,
    pub mean_effort: [f64; 2],
    /// Positive efforts resist the two coordinates [N s,Nm s].
    pub effort_impulse: [f64; 2],
    /// First time moment over this constant-discrete-effort phase [N s²,Nm s²].
    pub effort_time_moment: [f64; 2],
    pub contact_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub applied_work_j: f64,
    pub energy_residual_j: f64,
    pub source_work_residual_j: f64,
    pub coordinate_residual: [f64; 2],
    pub rate_residual: [f64; 2],
    pub maximum_source_effort: [f64; 2],
    pub iterations: usize,
    pub evaluations: usize,
}
impl JkrTorsionalMotionStep {
    pub const fn initial_state(&self) -> JkrTorsionalMotionState {
        self.initial
    }
    pub const fn state(&self) -> JkrTorsionalMotionState {
        self.state
    }
}

/// Immutable source and held projection; construction and repeated proposals
/// allocate nothing. Accepted history is never used as Newton scratch.
#[derive(Clone, Debug)]
pub struct JkrTorsionalMotion<'a> {
    source: &'a SourceJkrMixedMode,
    mobility: [[f64; 2]; 2],
    scale: [f64; 2],
    inverse_correlation: [[f64; 2]; 2],
    budget: JkrTorsionalMotionBudget,
    identity: ContentHash,
}

struct Trial {
    position: [f64; 2],
    path: JkrTorsionalPath,
    residual: [f64; 2],
    norm: f64,
    normalized_rhs: [f64; 2],
    jacobian: [[f64; 2]; 2],
}

impl<'a> JkrTorsionalMotion<'a> {
    pub fn new(
        source: &'a SourceJkrMixedMode,
        mut mobility: [[f64; 2]; 2],
        budget: JkrTorsionalMotionBudget,
    ) -> Result<Self, JkrTorsionalMotionError> {
        validate_budget(budget)?;
        if mobility.iter().flatten().any(|x| !x.is_finite())
            || mobility[0][0] <= 0.0
            || mobility[1][1] <= 0.0
        {
            return Err(JkrTorsionalMotionError::Invalid(
                "finite positive mobility diagonals required",
            ));
        }
        let scale = [mobility[0][0].sqrt(), mobility[1][1].sqrt()];
        let cross = [
            (mobility[0][1] / scale[0]) / scale[1],
            (mobility[1][0] / scale[1]) / scale[0],
        ];
        if cross.iter().any(|x| !x.is_finite()) || (cross[0] - cross[1]).abs() > 64.0 * f64::EPSILON
        {
            return Err(JkrTorsionalMotionError::Invalid(
                "reciprocal mixed-unit mobility required",
            ));
        }
        // Average only admitted roundoff-level reciprocity differences; every
        // later reaction and identity uses this same canonical reciprocal W.
        mobility[0][1] = mobility[0][1].midpoint(mobility[1][0]);
        mobility[1][0] = mobility[0][1];
        let correlation = (mobility[0][1] / scale[0]) / scale[1];
        let determinant = (1.0 - correlation) * (1.0 + correlation);
        if !determinant.is_finite() || determinant <= 64.0 * f64::EPSILON {
            return Err(JkrTorsionalMotionError::Invalid(
                "resolved positive native mobility required",
            ));
        }
        let inverse_correlation = [
            [1.0 / determinant, -correlation / determinant],
            [-correlation / determinant, 1.0 / determinant],
        ];
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-torsional-motion.v1");
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
            budget,
            identity: hash.finalize(),
        })
    }
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub const fn mobility(&self) -> [[f64; 2]; 2] {
        self.mobility
    }
    pub const fn budget(&self) -> JkrTorsionalMotionBudget {
        self.budget
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Begin with explicitly authored attached coordinates; no bond nucleation
    /// or hidden twist preload is inferred. Exactly one original query is used.
    pub fn initial_state(
        &self,
        coordinates: [f64; 2],
        rates: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrTorsionalMotionState, JkrTorsionalMotionError> {
        poll(gate)?;
        bounded(rates, self.budget.maximum_rate, "initial rates")?;
        let path = self
            .source
            .torsional_discrete_path(coordinates, coordinates, gate)?;
        bounded(
            path.mean_effort,
            self.budget.maximum_effort,
            "initial source efforts",
        )?;
        self.kinetic_energy(rates)?;
        poll(gate)?;
        Ok(JkrTorsionalMotionState {
            identity: self.identity,
            contact: path.contact,
            rates,
            time_s: 0.0,
            steps: 0,
        })
    }

    /// A held acceleration is the projection J M^-1 f_external from the actual
    /// mechanical owner. The original normal/twist cross response is retained.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        previous: JkrTorsionalMotionState,
        dt_s: f64,
        external_acceleration: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrTorsionalMotionStep, JkrTorsionalMotionError> {
        poll(gate)?;
        self.check(previous)?;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.budget.maximum_step_s {
            return Err(JkrTorsionalMotionError::Range("phase duration"));
        }
        let time_s = finite(previous.time_s + dt_s, "motion clock")?;
        if time_s <= previous.time_s {
            return Err(JkrTorsionalMotionError::Range(
                "unrepresentable clock increment",
            ));
        }
        let steps = previous
            .steps
            .checked_add(1)
            .ok_or(JkrTorsionalMotionError::Budget)?;
        let half_h2 = finite(0.5 * dt_s * dt_s, "squared duration")?;
        if half_h2 == 0.0 {
            return Err(JkrTorsionalMotionError::Range("squared duration underflow"));
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
                return Err(JkrTorsionalMotionError::NotConverged {
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
                    Err(JkrTorsionalMotionError::Budget) => {
                        return Err(JkrTorsionalMotionError::Budget);
                    }
                    Err(JkrTorsionalMotionError::Cancelled) => {
                        return Err(JkrTorsionalMotionError::Cancelled);
                    }
                    Err(error) => last_error = Some(error),
                }
                fraction *= 0.5;
            }
            trial = next.ok_or_else(|| {
                last_error.unwrap_or(JkrTorsionalMotionError::NotConverged {
                    residual: trial.residual,
                    tolerance,
                })
            })?;
        }
        unreachable!("bounded Newton loop always returns")
    }

    fn trial(
        &self,
        previous: JkrTorsionalMotionState,
        to: [f64; 2],
        free: [f64; 2],
        half_h2: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Trial, JkrTorsionalMotionError> {
        poll(gate)?;
        let increment = sub(to, previous.coordinates());
        bounded(
            increment,
            self.budget.maximum_coordinate_increment,
            "coordinate increment",
        )?;
        let path = self.source.torsional_path_counted(
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
            (residual[0] / self.budget.maximum_coordinate_increment[0])
                .abs()
                .max((residual[1] / self.budget.maximum_coordinate_increment[1]).abs()),
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
            return Err(JkrTorsionalMotionError::InertialTangent);
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
        previous: JkrTorsionalMotionState,
        trial: Trial,
        applied: [f64; 2],
        external: [f64; 2],
        dt: f64,
        time_s: f64,
        steps: u64,
        iterations: usize,
        mut evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalMotionStep, JkrTorsionalMotionError> {
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
        let midpoint: [f64; 2] = core::array::from_fn(|i| previous.rates[i].midpoint(rates[i]));
        let increment = sub(trial.position, previous.coordinates());
        let coordinate_residual = checked_array(
            core::array::from_fn(|i| increment[i] - dt * midpoint[i]),
            "midpoint coordinate balance",
        )?;
        let coordinate_tolerance =
            self.coordinate_tolerance(previous.coordinates(), trial.position)?;
        if !within(coordinate_residual, coordinate_tolerance) {
            return Err(JkrTorsionalMotionError::NotConverged {
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
            return Err(JkrTorsionalMotionError::RateBalance {
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
        for i in 0..2 {
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
        let covered = (0..2).all(|i| {
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
            self.source.admit_torsional_box(
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
                scaled_midpoint,
                mul(self.inverse_correlation, scaled_change),
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
            return Err(JkrTorsionalMotionError::EnergyBalance {
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
        Ok(JkrTorsionalMotionStep {
            initial: previous,
            state: JkrTorsionalMotionState {
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
        q: [[f64; 2]; 2],
        radii: [f64; 2],
    ) -> Result<[f64; 2], JkrTorsionalMotionError> {
        let normal = self.source.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let r = normal.effective_radius_m();
        let mut force = 0.0_f64;
        // N=2Ea(delta-a²/(3R)) is affine in delta. At fixed delta its only
        // positive-radius stationary point is sqrt(R delta); include it when
        // it lies inside the conservative radius interval.
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
                    let value =
                        finite((4.0 / 3.0) * e * a * delta, "source normal effort extremum")?;
                    force = force.max(value.abs());
                }
            }
        }
        let a = radii[1];
        let twist = q[0][1].abs().max(q[1][1].abs());
        let moment =
            ((16.0 / 3.0) * self.source.torsional_reduced_shear_modulus_pa() * a) * a * (a * twist);
        checked_array(
            [
                force.next_up(),
                if moment == 0.0 { 0.0 } else { moment.next_up() },
            ],
            "source effort enclosure",
        )
    }

    fn coordinate_tolerance(
        &self,
        from: [f64; 2],
        to: [f64; 2],
    ) -> Result<[f64; 2], JkrTorsionalMotionError> {
        checked_array(
            core::array::from_fn(|i| {
                self.budget.coordinate_absolute_tolerance[i]
                    + self.budget.coordinate_relative_tolerance * from[i].abs().max(to[i].abs())
            }),
            "coordinate tolerance",
        )
    }
    fn check(&self, state: JkrTorsionalMotionState) -> Result<(), JkrTorsionalMotionError> {
        if state.identity != self.identity {
            return Err(JkrTorsionalMotionError::StateMismatch);
        }
        Ok(())
    }
    fn kinetic_energy(&self, rates: [f64; 2]) -> Result<f64, JkrTorsionalMotionError> {
        let scaled = checked_array(
            core::array::from_fn(|i| rates[i] / self.scale[i]),
            "normalized rates",
        )?;
        let energy = finite(
            0.5 * dot(scaled, mul(self.inverse_correlation, scaled)),
            "projected kinetic energy",
        )?;
        if energy < 0.0 {
            return Err(JkrTorsionalMotionError::Range(
                "negative projected kinetic energy",
            ));
        }
        Ok(energy)
    }

    /// Rebind an authenticated native projection without a physical step.
    /// The native owner separately verifies unchanged body momenta and energy.
    /// Changing this observable kinetic partition supplies no work or heat.
    pub fn rebind_projection(
        &self,
        previous: JkrTorsionalMotionState,
        new_mobility: [[f64; 2]; 2],
        new_rates: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrTorsionalProjectionContinuation<'a>, JkrTorsionalMotionError> {
        poll(gate)?;
        self.check(previous)?;
        let motion = Self::new(self.source, new_mobility, self.budget)?;
        bounded(new_rates, self.budget.maximum_rate, "reprojected rates")?;
        let before = self.kinetic_energy(previous.rates)?;
        let after = motion.kinetic_energy(new_rates)?;
        let report = JkrTorsionalProjectionReport {
            projected_kinetic_energy_before_j: before,
            projected_kinetic_energy_after_j: after,
            projected_kinetic_repartition_j: finite(
                (Dd::from_f64(after) - Dd::from_f64(before)).to_f64(),
                "projected kinetic repartition",
            )?,
        };
        let state = JkrTorsionalMotionState {
            identity: motion.identity,
            rates: new_rates,
            ..previous
        };
        poll(gate)?;
        Ok(JkrTorsionalProjectionContinuation {
            previous,
            motion,
            state,
            report,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalProjectionReport {
    pub projected_kinetic_energy_before_j: f64,
    pub projected_kinetic_energy_after_j: f64,
    /// Diagnostic repartition, not physical body work, loss or heat [J].
    pub projected_kinetic_repartition_j: f64,
}
#[derive(Clone, Debug)]
pub struct JkrTorsionalProjectionContinuation<'a> {
    previous: JkrTorsionalMotionState,
    motion: JkrTorsionalMotion<'a>,
    state: JkrTorsionalMotionState,
    report: JkrTorsionalProjectionReport,
}
impl<'a> JkrTorsionalProjectionContinuation<'a> {
    pub const fn previous_state(&self) -> JkrTorsionalMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrTorsionalMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrTorsionalMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrTorsionalProjectionReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrTorsionalMotion<'a>,
        JkrTorsionalMotionState,
        JkrTorsionalProjectionReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

fn validate_budget(b: JkrTorsionalMotionBudget) -> Result<(), JkrTorsionalMotionError> {
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
        Err(JkrTorsionalMotionError::Invalid(
            "finite physical bounds, nonnegative tolerances and bounded work required",
        ))
    } else {
        Ok(())
    }
}

/// Extrema of each actual quadratic coordinate, evaluated by de Casteljau.
/// Outward ULP padding screens floating-point roundoff; no interval certificate
/// is claimed. Correlations between coordinates may make the box conservative.
fn quadratic_bounds(points: [[f64; 2]; 3]) -> Result<[[f64; 2]; 2], JkrTorsionalMotionError> {
    let mut bounds = [
        core::array::from_fn(|i| points[0][i].min(points[2][i])),
        core::array::from_fn(|i| points[0][i].max(points[2][i])),
    ];
    for i in 0..2 {
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

fn coercive(a: [[f64; 2]; 2]) -> bool {
    if a.iter().flatten().any(|x| !x.is_finite()) || a[0][0] <= 0.0 || a[1][1] <= 0.0 {
        return false;
    }
    let cross = (a[0][1].midpoint(a[1][0]) / a[0][0].sqrt()) / a[1][1].sqrt();
    cross.is_finite() && (1.0 - cross) * (1.0 + cross) > 64.0 * f64::EPSILON
}
fn solve(a: [[f64; 2]; 2], rhs: [f64; 2]) -> Result<[f64; 2], JkrTorsionalMotionError> {
    let scale = a.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    let a = a.map(|row| row.map(|x| x / scale));
    let b = rhs.map(|x| x / scale);
    let determinant = (Dd::from_f64(a[0][0]) * Dd::from_f64(a[1][1])
        - Dd::from_f64(a[0][1]) * Dd::from_f64(a[1][0]))
    .to_f64();
    if !determinant.is_finite() || determinant <= 0.0 {
        return Err(JkrTorsionalMotionError::InertialTangent);
    }
    checked_array(
        [
            dot([a[1][1], -a[0][1]], b) / determinant,
            dot([-a[1][0], a[0][0]], b) / determinant,
        ],
        "normalized Newton solution",
    )
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    (Dd::from_f64(a[0]) * Dd::from_f64(b[0]) + Dd::from_f64(a[1]) * Dd::from_f64(b[1])).to_f64()
}
fn mul(a: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    a.map(|row| dot(row, x))
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn within(value: [f64; 2], limit: [f64; 2]) -> bool {
    (0..2).all(|i| value[i].is_finite() && value[i].abs() <= limit[i])
}
fn bounded(
    value: [f64; 2],
    limit: [f64; 2],
    field: &'static str,
) -> Result<(), JkrTorsionalMotionError> {
    if within(value, limit) {
        Ok(())
    } else {
        Err(JkrTorsionalMotionError::Range(field))
    }
}
fn checked_array(
    value: [f64; 2],
    field: &'static str,
) -> Result<[f64; 2], JkrTorsionalMotionError> {
    if value.iter().all(|x| x.is_finite()) {
        Ok(value)
    } else {
        Err(JkrTorsionalMotionError::Range(field))
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrTorsionalMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrTorsionalMotionError::Range(field))
    }
}
fn poll(gate: &CancelGate) -> Result<(), JkrTorsionalMotionError> {
    if gate.is_requested() {
        Err(JkrTorsionalMotionError::Cancelled)
    } else {
        Ok(())
    }
}
