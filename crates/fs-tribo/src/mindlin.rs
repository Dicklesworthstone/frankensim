//! Fixed-load, proportional Hertz--Cattaneo--Mindlin tangential contact.
//!
//! With `C = mu P`, `k0 = 8 a / sum((2-nu_i)/G_i)` and `d = 3 C/(2 k0)`,
//! virgin loading is `Q(u) = C [1-(1-u/d)^(3/2)]`, up to gross slip.
//! The Hertz continuum of elastic Coulomb elements has normalized thresholds
//! `r in [0,1]` and stiffness measure `k0 dr/(2 sqrt(1-r))`. Its retained
//! strain is piecewise affine in r; clipping each piece to `[-r,r]` implements
//! arbitrary scalar reversals and return-point memory. All force, spring-energy
//! and slip-work integrals are polynomial integrals after `t = sqrt(1-r)`.
//! There is no finite spring discretization or authored tangential stiffness.
//!
//! Each interval returns the exact path-average force, not the endpoint force.
//! Thus the existing held-effort port closes `delta(H) + D = F_average du`, and
//! its body load line is solved with that same average. Root trials allocate
//! nothing; only the accepted candidate retains new memory intervals.
//!
//! This is the classical decoupled Cattaneo--Mindlin approximation, for a
//! circular Hertz point patch at fixed normal load/radius and one fixed tangent
//! direction. It does not solve changing normal load, nonproportional shear,
//! torsion, normal/tangential elastic coupling, or spatial traction/wear fields.
//! The continuum representation follows Popov, Hess and Willert, *Handbook of
//! Contact Mechanics* (2019), sections 4.4 and 4.6.3, equations 4.21--4.25,
//! 4.57 and 4.59: <https://doi.org/10.1007/978-3-662-58709-6_4>.
//! Its relation to Hertz--Mindlin contact is also derived by Lyashenko & Popov
//! (2015), Methods: <https://doi.org/10.1038/srep08479>.

use crate::{
    FrictionRegime, TangentialLaw, TangentialStep, TangentialStepInput, TriboError, checked_add,
    checked_mul, finite, nonnegative_finite, positive_finite,
};

/// Explicit normal loading, unloading and recontact with retained shear memory.
pub mod loading;

/// Largest admitted number of retained analytic history intervals.
pub const MAX_HISTORY_INTERVALS: usize = 4096;

/// Source-neutral inputs. The caller binds the Hertz patch, two elastic bulk
/// responses and a constant local friction coefficient; none is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MindlinPointParameters {
    /// Frozen compressive patch load [N].
    pub normal_force_n: f64,
    /// Frozen circular Hertz contact radius [m], not the sphere radius.
    pub contact_radius_m: f64,
    /// Shear moduli of the two ordered bodies [Pa].
    pub shear_moduli_pa: [f64; 2],
    /// Poisson ratios of the same ordered bodies, `-1 < nu <= 1/2`.
    pub poisson_ratios: [f64; 2],
    /// Constant local Coulomb coefficient; zero is the frictionless limit.
    pub friction_coefficient: f64,
    /// Caller memory ceiling, inclusive, from 1 to [`MAX_HISTORY_INTERVALS`].
    /// The initial state uses one interval; unresolved nested reversals refuse.
    pub maximum_history_intervals: usize,
}

/// Analytically integrated Hertz--Cattaneo--Mindlin point-contact law.
#[derive(Clone, Debug, PartialEq)]
pub struct MindlinDeresiewiczPoint {
    parameters: MindlinPointParameters,
    stiffness_n_per_m: f64,
    capacity_n: f64,
    critical_slip_m: f64,
    force_scale_n: f64,
    energy_scale_j: f64,
}

// On [previous.upper, upper], z(r)/d = slope*r + offset. Slopes are -1, 0, 1.
// Adjacent identical polynomials are merged; no reversal is silently dropped.
#[derive(Clone, Copy, Debug, PartialEq)]
struct MemoryInterval {
    upper: f64,
    slope: f64,
    offset: f64,
}

/// Complete fixed-patch checkpoint. The private law binding prevents using
/// retained strain with another modulus, contact radius or friction capacity.
#[derive(Clone, Debug, PartialEq)]
pub struct MindlinState {
    owner: MindlinDeresiewiczPoint,
    memory: Vec<MemoryInterval>,
    slip_m: f64,
    endpoint_force_n: f64,
    stored_energy_j: f64,
}

impl MindlinState {
    /// Accepted endpoint force entering the interface [N]. This differs from
    /// the path-average held effort returned by [`TangentialLaw::advance`].
    #[must_use]
    pub const fn endpoint_force_n(&self) -> f64 {
        self.endpoint_force_n
    }

    /// Accumulated solved displacement along the unchanged first tangent [m].
    #[must_use]
    pub const fn slip_m(&self) -> f64 {
        self.slip_m
    }

    /// Number of analytic intervals carrying the complete reversal memory.
    #[must_use]
    pub fn history_intervals(&self) -> usize {
        self.memory.len()
    }
}

impl MindlinDeresiewiczPoint {
    /// Derives compliance from the two elastic bodies and the Hertz radius.
    pub fn new(parameters: MindlinPointParameters) -> Result<Self, TriboError> {
        positive_finite(parameters.normal_force_n, "mindlin_normal_force_n")?;
        positive_finite(parameters.contact_radius_m, "mindlin_contact_radius_m")?;
        nonnegative_finite(
            parameters.friction_coefficient,
            "mindlin_friction_coefficient",
        )?;
        if !(1..=MAX_HISTORY_INTERVALS).contains(&parameters.maximum_history_intervals) {
            return Err(TriboError::InvalidInput {
                field: "mindlin_maximum_history_intervals",
            });
        }
        let mut compliance = 0.0;
        for (shear, poisson) in parameters
            .shear_moduli_pa
            .into_iter()
            .zip(parameters.poisson_ratios)
        {
            positive_finite(shear, "mindlin_shear_modulus_pa")?;
            if !(poisson > -1.0 && poisson <= 0.5) {
                return Err(TriboError::InvalidInput {
                    field: "mindlin_poisson_ratio",
                });
            }
            compliance = checked_add(
                compliance,
                (2.0 - poisson) / shear,
                "mindlin_elastic_compliance",
            )?;
        }
        positive_finite(compliance, "mindlin_elastic_compliance")?;
        let stiffness_n_per_m = 8.0 * parameters.contact_radius_m / compliance;
        positive_finite(stiffness_n_per_m, "mindlin_initial_stiffness_n_per_m")?;
        let capacity_n = checked_mul(
            parameters.friction_coefficient,
            parameters.normal_force_n,
            "mindlin_capacity_n",
        )?;
        if parameters.friction_coefficient > 0.0 {
            positive_finite(capacity_n, "mindlin_capacity_n")?;
        }
        let force_scale_n = checked_mul(1.5, capacity_n, "mindlin_force_scale_n")?;
        let critical_slip_m = force_scale_n / stiffness_n_per_m;
        if capacity_n > 0.0 {
            positive_finite(critical_slip_m, "mindlin_gross_slip_displacement_m")?;
        }
        let energy_scale_j = checked_mul(force_scale_n, critical_slip_m, "mindlin_energy_scale_j")?;
        Ok(Self {
            parameters,
            stiffness_n_per_m,
            capacity_n,
            critical_slip_m,
            force_scale_n,
            energy_scale_j,
        })
    }

    /// Original frozen patch and ordered elastic inputs.
    #[must_use]
    pub const fn parameters(&self) -> MindlinPointParameters {
        self.parameters
    }

    /// No-slip tangent stiffness `8 a G*` [N/m].
    #[must_use]
    pub const fn initial_stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }

    /// Gross sliding capacity `mu P` [N].
    #[must_use]
    pub const fn gross_slip_capacity_n(&self) -> f64 {
        self.capacity_n
    }

    /// Virgin-loading displacement at complete slip [m]; zero if frictionless.
    #[must_use]
    pub const fn gross_slip_displacement_m(&self) -> f64 {
        self.critical_slip_m
    }

    fn response(
        &self,
        state: &MindlinState,
        increment_m: f64,
        mut memory: Option<&mut Vec<MemoryInterval>>,
    ) -> Result<Response, TriboError> {
        let increment = increment_m / self.critical_slip_m;
        finite(increment, "mindlin_normalized_increment")?;
        if increment == 0.0 {
            if let Some(out) = memory {
                out.extend_from_slice(&state.memory);
            }
            return Ok(Response {
                average_n: state.endpoint_force_n,
                endpoint_n: state.endpoint_force_n,
                energy_j: state.stored_energy_j,
                stored_change_j: 0.0,
                heat_j: 0.0,
                derivative_n_per_m: 0.5 * self.stiffness_n_per_m,
            });
        }
        let direction = increment.signum();
        let mut totals = Integrals::default();
        let mut lower = 0.0;
        for old in &state.memory {
            let shifted = checked_add(old.offset, increment, "mindlin_history_increment")?;
            let rate = 1.0 - direction * old.slope;
            let crossing = if rate == 0.0 {
                if direction * shifted > 0.0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                }
            } else {
                direction * shifted / rate
            };
            let split = crossing.clamp(lower, old.upper);
            if split > lower {
                let moments = Moments::new(lower, split);
                let old_mean = old.slope * moments.mean + old.offset;
                let elastic_mean = direction * moments.mean - old_mean;
                let elastic_slope = direction - old.slope;
                let elastic_square = moments.weight
                    * (elastic_mean * elastic_mean
                        + elastic_slope * elastic_slope * moments.variance);
                totals.average +=
                    direction * moments.weight * moments.mean - 0.5 * (elastic_square / increment);
                totals.derivative += 0.5 * ((elastic_square / increment) / increment);
                totals.stored_change += 0.5
                    * moments.weight
                    * (elastic_mean * (direction * moments.mean + old_mean)
                        + (1.0 - old.slope * old.slope) * moments.variance);
                // At an internal crossing the plastic slip is exactly zero.
                // This positive moment form avoids subtracting nearly equal
                // O(Q^2) terms to obtain the small-amplitude O(Q^3) heat.
                let plastic_at_upper = if crossing <= old.upper {
                    0.0
                } else {
                    direction * shifted - rate * split
                };
                nonnegative_finite(plastic_at_upper, "mindlin_plastic_increment")?;
                totals.heat += moments.weight
                    * (plastic_at_upper * moments.mean + rate * moments.radius_times_upper_gap);
                totals.endpoint += direction * moments.weight * moments.mean;
                totals.energy +=
                    0.5 * moments.weight * (moments.mean * moments.mean + moments.variance);
                if let Some(out) = memory.as_deref_mut() {
                    self.append_memory(
                        out,
                        MemoryInterval {
                            upper: split,
                            slope: direction,
                            offset: 0.0,
                        },
                    )?;
                }
            }
            if split < old.upper {
                let moments = Moments::new(split, old.upper);
                let old_mean = old.slope * moments.mean + old.offset;
                let average = old_mean + 0.5 * increment;
                let new_mean = old.slope * moments.mean + shifted;
                totals.average += moments.weight * average;
                totals.derivative += 0.5 * moments.weight;
                totals.stored_change += increment * moments.weight * average;
                totals.endpoint += moments.weight * new_mean;
                totals.energy += 0.5
                    * moments.weight
                    * (new_mean * new_mean + old.slope * old.slope * moments.variance);
                if let Some(out) = memory.as_deref_mut() {
                    self.append_memory(
                        out,
                        MemoryInterval {
                            upper: old.upper,
                            slope: old.slope,
                            offset: shifted,
                        },
                    )?;
                }
            }
            lower = old.upper;
        }
        let response = Response {
            average_n: self.force_scale_n * totals.average,
            endpoint_n: self.force_scale_n * totals.endpoint,
            energy_j: self.energy_scale_j * totals.energy,
            stored_change_j: self.energy_scale_j * totals.stored_change,
            heat_j: self.energy_scale_j * totals.heat,
            derivative_n_per_m: self.stiffness_n_per_m * totals.derivative,
        };
        for value in [
            response.average_n,
            response.endpoint_n,
            response.stored_change_j,
        ] {
            finite(value, "mindlin_response")?;
        }
        for value in [
            response.energy_j,
            response.heat_j,
            response.derivative_n_per_m,
        ] {
            nonnegative_finite(value, "mindlin_energy_or_tangent")?;
        }
        Ok(response)
    }

    fn append_memory(
        &self,
        out: &mut Vec<MemoryInterval>,
        interval: MemoryInterval,
    ) -> Result<(), TriboError> {
        if let Some(last) = out.last_mut()
            && last.slope == interval.slope
            && last.offset == interval.offset
        {
            last.upper = interval.upper;
            return Ok(());
        }
        if out.len() == self.parameters.maximum_history_intervals {
            return Err(TriboError::OutsideApplicability {
                field: "mindlin_history_intervals",
                value: (out.len() + 1) as f64,
                minimum: 1.0,
                maximum: self.parameters.maximum_history_intervals as f64,
            });
        }
        out.push(interval);
        Ok(())
    }

    fn solve_increment(
        &self,
        state: &MindlinState,
        x: &TangentialStepInput,
    ) -> Result<(f64, usize), TriboError> {
        let free = x.slip_increment_m[0];
        let compliance = x.compliance_m_per_n;
        if compliance == 0.0 {
            return Ok((free, 0));
        }
        let bound = self.capacity_n.min(x.solve.maximum_force_n);
        let evaluate = |force: f64| -> Result<(f64, f64, f64, f64), TriboError> {
            let increment = free - checked_mul(compliance, force, "mindlin_body_response_m")?;
            finite(increment, "mindlin_resolved_increment_m")?;
            let response = self.response(state, increment, None)?;
            Ok((
                force - response.average_n,
                x.solve.tolerance(force, response.average_n)?,
                increment,
                1.0 + compliance * response.derivative_n_per_m,
            ))
        };
        let mut lo = -bound;
        let mut hi = bound;
        for (force, lower) in [(lo, true), (hi, false)] {
            let (residual, tolerance, increment, _) = evaluate(force)?;
            if residual.abs() <= tolerance {
                return Ok((increment, 0));
            }
            if (lower && residual > 0.0) || (!lower && residual < 0.0) {
                return Err(TriboError::OutsideApplicability {
                    field: "maximum_force_n",
                    value: (force - residual).abs(),
                    minimum: 0.0,
                    maximum: x.solve.maximum_force_n,
                });
            }
        }
        let mut force = 0.0;
        let mut final_residual = 0.0;
        let mut final_tolerance = 0.0;
        for iteration in 1..=x.solve.max_iterations {
            let (residual, tolerance, increment, derivative) = evaluate(force)?;
            if residual.abs() <= tolerance {
                return Ok((increment, iteration));
            }
            final_residual = residual;
            final_tolerance = tolerance;
            if residual < 0.0 {
                lo = force;
            } else {
                hi = force;
            }
            let proposed = force - residual / derivative;
            force = if proposed.is_finite() && proposed > lo && proposed < hi {
                proposed
            } else {
                f64::midpoint(lo, hi)
            };
        }
        Err(TriboError::TangentialSolve {
            residual_n: final_residual,
            tolerance_n: final_tolerance,
            iterations: x.solve.max_iterations,
        })
    }
}

impl TangentialLaw for MindlinDeresiewiczPoint {
    type State = MindlinState;

    fn model_id(&self) -> &str {
        "cattaneo-mindlin-fixed-load-proportional-v1"
    }

    fn initial_state(&self) -> Self::State {
        MindlinState {
            owner: self.clone(),
            memory: vec![MemoryInterval {
                upper: 1.0,
                slope: 0.0,
                offset: 0.0,
            }],
            slip_m: 0.0,
            endpoint_force_n: 0.0,
            stored_energy_j: 0.0,
        }
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if state.owner != *self {
            f64::NAN
        } else {
            state.stored_energy_j
        }
    }

    fn advance(
        &self,
        state: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        x.solve.validate()?;
        if state.owner != *self {
            return Err(TriboError::InvalidInput {
                field: "mindlin_checkpoint_owner",
            });
        }
        positive_finite(x.normal_n, "normal_force")?;
        if x.normal_n != self.parameters.normal_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "mindlin_fixed_normal_force_n",
                value: x.normal_n,
                minimum: self.parameters.normal_force_n,
                maximum: self.parameters.normal_force_n,
            });
        }
        positive_finite(x.dt_s, "dt_s")?;
        nonnegative_finite(x.compliance_m_per_n, "compliance_m_per_n")?;
        if let Some(t) = x.contact_temperature_k {
            positive_finite(t, "temperature_kelvin")?;
        }
        if let Some(p) = x.nominal_pressure_pa {
            nonnegative_finite(p, "contact_pressure_pa")?;
        }
        for value in x.slip_increment_m {
            finite(value, "slip_increment_m")?;
        }
        if x.slip_increment_m[1] != 0.0 {
            return Err(TriboError::InvalidInput {
                field: "mindlin_requires_proportional_first_axis_motion",
            });
        }
        let (increment, iterations) = if self.capacity_n == 0.0 {
            (x.slip_increment_m[0], 0)
        } else {
            self.solve_increment(state, x)?
        };
        let mut memory = Vec::with_capacity(
            (state.memory.len() + 1).min(self.parameters.maximum_history_intervals),
        );
        let response = if self.capacity_n == 0.0 {
            memory.extend_from_slice(&state.memory);
            Response {
                average_n: 0.0,
                endpoint_n: 0.0,
                energy_j: 0.0,
                stored_change_j: 0.0,
                heat_j: 0.0,
                derivative_n_per_m: 0.0,
            }
        } else {
            self.response(state, increment, Some(&mut memory))?
        };
        if response.average_n.abs() > x.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "maximum_force_n",
                value: response.average_n.abs(),
                minimum: 0.0,
                maximum: x.solve.maximum_force_n,
            });
        }
        let slip_m = checked_add(state.slip_m, increment, "mindlin_accumulated_slip_m")?;
        let work = checked_mul(response.average_n, increment, "mindlin_work_j")?;
        let closure_residual_j = response.stored_change_j + response.heat_j - work;
        let compatibility_residual_m =
            increment - (x.slip_increment_m[0] - x.compliance_m_per_n * response.average_n);
        finite(closure_residual_j, "mindlin_closure_residual_j")?;
        finite(compatibility_residual_m, "mindlin_compatibility_residual_m")?;
        Ok(TangentialStep {
            traction_n: [response.average_n, 0.0],
            slip_increment_m: [increment, 0.0],
            state: MindlinState {
                owner: self.clone(),
                memory,
                slip_m,
                endpoint_force_n: response.endpoint_n,
                stored_energy_j: response.energy_j,
            },
            regime: if self.capacity_n == 0.0 || response.heat_j > 0.0 {
                FrictionRegime::Sliding
            } else {
                FrictionRegime::Sticking
            },
            stored_change_j: response.stored_change_j,
            dissipated_j: response.heat_j,
            closure_residual_j,
            compatibility_residual_m: [compatibility_residual_m, 0.0],
            iterations,
        })
    }
}

#[derive(Default)]
struct Integrals {
    average: f64,
    derivative: f64,
    endpoint: f64,
    energy: f64,
    stored_change: f64,
    heat: f64,
}

struct Response {
    average_n: f64,
    endpoint_n: f64,
    energy_j: f64,
    stored_change_j: f64,
    heat_j: f64,
    derivative_n_per_m: f64,
}

// Exact moments of r = 1-t^2 under uniform dt on this interval. Positive,
// factored forms preserve small edge-slip work and avoid subtracting powers.
struct Moments {
    weight: f64,
    mean: f64,
    variance: f64,
    radius_times_upper_gap: f64,
}

impl Moments {
    fn new(lower: f64, upper: f64) -> Self {
        let a = (1.0 - lower).sqrt();
        let b = (1.0 - upper).sqrt();
        let width = (upper - lower) / (a + b);
        let square = width * width;
        let mean = lower + width * (2.0 * a + b) / 3.0;
        let variance = square * (a * b / 3.0 + 4.0 * square / 45.0);
        let upper_gap = width * (a + 2.0 * b) / 3.0;
        let radius_times_upper_gap =
            lower * upper_gap + (2.0 / 3.0) * square * (a * b + square / 5.0);
        Self {
            weight: width,
            mean,
            variance,
            radius_times_upper_gap,
        }
    }
}
