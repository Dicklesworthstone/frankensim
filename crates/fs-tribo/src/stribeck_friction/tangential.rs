//! Retained scalar stable branches through the shared held-effort interface.
//!
//! The first loaded motion selects a tangent axis. Continued sliding and resolved
//! reversal use that same axis; objective frame transport rotates it with the
//! accepted slip vector. An oblique change while sliding refuses. This is the
//! original one-dimensional hysteresis rule, not an invented two-axis recapture
//! rule. Arrest or loss of normal load releases the occupied direction.
//!
//! Critical-point preparation and root updates share the caller's step budget.
//! Matching load/mobility reuses the prepared line. Explicit preparation is also
//! available as a bounded, zero-motion candidate for difficult changing loads.
//! Nothing mutates accepted history, and exported heat is F dot actual slip.

use super::{ContactBranch, ImplicitFrictionError, StribeckLoadLine};
use crate::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};
use crate::{
    DryFrictionState, DryInterfaceSystemCard, FrictionLaw, FrictionRegime, TangentialLaw,
    TangentialSolveConfig, TangentialStep, TangentialStepInput, TriboError, checked_add,
    checked_mul, nonnegative_finite, positive_finite,
};
use fs_math::det;

/// An original dry law with occupied-branch continuation instead of static
/// priority on every interval. Optional cards retain their original ordered
/// identity and state-domain checks; this owner admits no new material claims.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticFriction {
    law: FrictionLaw,
    card: Option<DryInterfaceSystemCard>,
}

/// Complete accepted branch/direction history and reusable numerical preparation.
/// The law binding prevents another coefficient set or interface from silently
/// inheriting this history. A prepared line changes no physical state or work.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticFrictionState {
    owner: HystereticFriction,
    branch: ContactBranch,
    tangent_axis: Option<[f64; 2]>,
    slip_m: [f64; 2],
    contact_temperature_k: Option<f64>,
    prepared: Option<StribeckLoadLine>,
    roots: u8,
    selected: u8,
}

impl HystereticFrictionState {
    /// Occupied stable branch; its sign is relative to `tangent_axis`.
    #[must_use]
    pub const fn branch(&self) -> ContactBranch {
        self.branch
    }
    /// Retained axis of a sliding branch in the checkpoint's current frame.
    #[must_use]
    pub const fn tangent_axis(&self) -> Option<[f64; 2]> {
        self.tangent_axis
    }
    /// Accumulated solved slip vector [m], in that same current tangent frame.
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }
    /// Last explicit thermal boundary; no thermal evolution is inferred.
    #[must_use]
    pub const fn contact_temperature_k(&self) -> Option<f64> {
        self.contact_temperature_k
    }
    /// Number and selected index of the last accepted discrete intersections.
    /// Released contact has `(0, 0)`; a prescribed sliding motion has `(1, 0)`.
    #[must_use]
    pub const fn intersections(&self) -> (u8, u8) {
        (self.roots, self.selected)
    }
    /// Prepared line, reusable only for its exact original load and mobility.
    #[must_use]
    pub const fn prepared_load_line(&self) -> Option<&StribeckLoadLine> {
        self.prepared.as_ref()
    }
}

impl HystereticFriction {
    /// Bind a validated constitutive law with caller-supplied coefficients.
    pub fn new(law: FrictionLaw) -> Result<Self, TriboError> {
        law.validate()?;
        Ok(Self { law, card: None })
    }

    /// Preserve an original ordered dry-interface card and its applicability.
    pub fn from_card(card: DryInterfaceSystemCard) -> Self {
        Self {
            law: card.friction_law().clone(),
            card: Some(card),
        }
    }

    /// Original immutable coefficient law.
    #[must_use]
    pub const fn friction_law(&self) -> &FrictionLaw {
        &self.law
    }
    /// Original optional interface card; its authority is never upgraded here.
    #[must_use]
    pub const fn interface_card(&self) -> Option<&DryInterfaceSystemCard> {
        self.card.as_ref()
    }

    /// Prepare a load line without advancing slip, branch, temperature or work.
    /// The returned candidate retains all accepted history. Matching preparation
    /// costs zero iterations. This is optional for ordinary steps; it permits an
    /// explicitly larger preparation budget near almost coincident folds.
    pub fn prepare_state(
        &self,
        state: &HystereticFrictionState,
        normal_n: f64,
        mobility_m_per_n_s: f64,
        maximum_preparation_iterations: usize,
    ) -> Result<(HystereticFrictionState, usize), TriboError> {
        self.require_state(state)?;
        if state
            .prepared
            .as_ref()
            .is_some_and(|line| line.matches(normal_n, mobility_m_per_n_s))
        {
            return Ok((state.clone(), 0));
        }
        let (line, iterations) = StribeckLoadLine::new_within(
            &self.law,
            normal_n,
            mobility_m_per_n_s,
            maximum_preparation_iterations,
        )
        .map_err(|error| match error {
            ImplicitFrictionError::Constitutive(error) => error,
            ImplicitFrictionError::InvalidInput(field) => TriboError::InvalidInput { field },
            _ => TriboError::InvalidInput {
                field: "hysteretic_load_line_preparation_budget",
            },
        })?;
        let mut candidate = state.clone();
        candidate.prepared = Some(line);
        Ok((candidate, iterations))
    }

    fn require_state(&self, state: &HystereticFrictionState) -> Result<(), TriboError> {
        if &state.owner != self {
            return Err(TriboError::InvalidInput {
                field: "hysteretic_friction_state_owner",
            });
        }
        Ok(())
    }

    fn admit_boundary(&self, input: &TangentialStepInput) -> Result<(), TriboError> {
        if let Some(card) = &self.card {
            let temperature = input
                .contact_temperature_k
                .ok_or(TriboError::InvalidInput {
                    field: "temperature_kelvin",
                })?;
            let pressure = input.nominal_pressure_pa.ok_or(TriboError::InvalidInput {
                field: "contact_pressure_pa",
            })?;
            card.applicability()
                .temperature_kelvin()
                .require_contains(temperature, "temperature_kelvin")?;
            card.applicability()
                .contact_pressure_pa()
                .require_contains(pressure, "contact_pressure_pa")?;
        }
        Ok(())
    }

    fn finish(
        &self,
        old: &HystereticFrictionState,
        input: &TangentialStepInput,
        mut step: TangentialStep<HystereticFrictionState>,
    ) -> Result<TangentialStep<HystereticFrictionState>, TriboError> {
        for i in 0..2 {
            step.state.slip_m[i] = checked_add(
                old.slip_m[i],
                step.slip_increment_m[i],
                "accumulated_slip_m",
            )?;
        }
        step.state.contact_temperature_k = input.contact_temperature_k;
        if let Some(card) = &self.card {
            let temperature = input
                .contact_temperature_k
                .ok_or(TriboError::InvalidInput {
                    field: "temperature_kelvin",
                })?;
            let pressure = input.nominal_pressure_pa.ok_or(TriboError::InvalidInput {
                field: "contact_pressure_pa",
            })?;
            let speed = finite(
                det::hypot(step.slip_increment_m[0], step.slip_increment_m[1]) / input.dt_s,
                "slip_speed_mps",
            )?;
            card.validate_state(DryFrictionState::new(temperature, pressure, speed)?)?;
        }
        Ok(step)
    }
}

impl TangentialLaw for HystereticFriction {
    type State = HystereticFrictionState;

    fn model_id(&self) -> &str {
        "dry-stable-branch-continuation-v1"
    }

    fn initial_state(&self) -> Self::State {
        HystereticFrictionState {
            owner: self.clone(),
            branch: ContactBranch::Fresh,
            tangent_axis: None,
            slip_m: [0.0; 2],
            contact_temperature_k: None,
            prepared: None,
            roots: 0,
            selected: 0,
        }
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if &state.owner == self { 0.0 } else { f64::NAN }
    }

    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.require_state(state)?;
        input.solve.validate()?;
        nonnegative_finite(input.normal_n, "normal_force")?;
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        positive_finite(input.dt_s, "dt_s")?;
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "temperature_kelvin")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "contact_pressure_pa")?;
        }
        for value in input.slip_increment_m {
            finite(value, "slip_increment_m")?;
        }
        self.admit_boundary(input)?;
        let distance = if input.slip_increment_m[1] == 0.0 {
            input.slip_increment_m[0].abs()
        } else {
            det::hypot(input.slip_increment_m[0], input.slip_increment_m[1])
        };
        finite(distance, "free_slip_distance_m")?;
        let mut axis = [0.0; 2];
        let mut sign = 1.0;
        if distance > 0.0 && input.normal_n > 0.0 {
            let free_axis = [
                input.slip_increment_m[0] / distance,
                input.slip_increment_m[1] / distance,
            ];
            axis = state.tangent_axis.unwrap_or(free_axis);
            let transverse = axis[0] * free_axis[1] - axis[1] * free_axis[0];
            // Only coordinate/rotation roundoff is admitted, not an authored
            // physical capture angle or a projection of an oblique motion.
            if transverse.abs() > 128.0 * f64::EPSILON {
                return Err(TriboError::InvalidInput {
                    field: "oblique_hysteretic_continuation",
                });
            }
            let alignment = axis[0] * free_axis[0] + axis[1] * free_axis[1];
            sign = if alignment < 0.0 { -1.0 } else { 1.0 };
        }

        // These limits are unique and already owned by the shared dry law.
        // In particular, release cannot leave an occupied sliding branch.
        if input.normal_n == 0.0 || input.compliance_m_per_n == 0.0 || distance == 0.0 {
            let unique = self.law.advance(&self.law.initial_state(), input)?;
            let mut candidate = state.clone();
            if input.normal_n == 0.0 || input.compliance_m_per_n == 0.0 {
                candidate.prepared = None;
            }
            candidate.branch = if input.normal_n == 0.0 {
                ContactBranch::Fresh
            } else if unique.regime == FrictionRegime::Sticking {
                ContactBranch::Stuck
            } else {
                ContactBranch::SlidingLow {
                    positive: sign > 0.0,
                }
            };
            candidate.tangent_axis = if matches!(candidate.branch, ContactBranch::SlidingLow { .. })
            {
                Some(axis)
            } else {
                None
            };
            candidate.roots = u8::from(input.normal_n > 0.0);
            candidate.selected = 0;
            return self.finish(
                state,
                input,
                TangentialStep {
                    traction_n: unique.traction_n,
                    slip_increment_m: unique.slip_increment_m,
                    state: candidate,
                    regime: unique.regime,
                    stored_change_j: unique.stored_change_j,
                    dissipated_j: unique.dissipated_j,
                    closure_residual_j: unique.closure_residual_j,
                    compatibility_residual_m: unique.compatibility_residual_m,
                    iterations: unique.iterations,
                },
            );
        }
        let mobility = finite(input.compliance_m_per_n / input.dt_s, "discrete_mobility")?;
        positive_finite(mobility, "discrete_mobility")?;
        let free_velocity = finite(sign * (distance / input.dt_s), "free_relative_velocity")?;
        let (line, preparation) = match &state.prepared {
            Some(line) if line.matches(input.normal_n, mobility) => (line.clone(), 0),
            _ => StribeckLoadLine::new_within(
                &self.law,
                input.normal_n,
                mobility,
                input.solve.max_iterations,
            )
            .map_err(|error| solve_error(error, input.solve, input.solve.max_iterations))?,
        };
        let (resolved, root_iterations) = line
            .step_within(
                free_velocity,
                state.branch,
                input.solve.max_iterations - preparation,
            )
            .map_err(|error| solve_error(error, input.solve, input.solve.max_iterations))?;
        let iterations = preparation + root_iterations;
        let magnitude = resolved.traction_n.abs();
        if magnitude > input.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "maximum_force_n",
                value: magnitude,
                minimum: 0.0,
                maximum: input.solve.maximum_force_n,
            });
        }
        let regime = if resolved.branch == ContactBranch::Stuck {
            FrictionRegime::Sticking
        } else {
            FrictionRegime::Sliding
        };
        let mut traction_n = [0.0; 2];
        let mut slip_increment_m = [0.0; 2];
        let mut compatibility_residual_m = [0.0; 2];
        for i in 0..2 {
            traction_n[i] = checked_mul(axis[i], resolved.traction_n, "traction_n")?;
            let actual = finite(
                input.slip_increment_m[i]
                    - checked_mul(input.compliance_m_per_n, traction_n[i], "body_response_m")?,
                "resolved_slip_m",
            )?;
            slip_increment_m[i] = if regime == FrictionRegime::Sticking {
                0.0
            } else {
                actual
            };
            compatibility_residual_m[i] =
                finite(slip_increment_m[i] - actual, "compatibility_residual_m")?;
        }
        if regime == FrictionRegime::Sliding {
            let speed = finite(
                det::hypot(slip_increment_m[0], slip_increment_m[1]) / input.dt_s,
                "slip_speed_mps",
            )?;
            positive_finite(speed, "zero_sliding_increment")?;
            let expected = checked_mul(
                input.normal_n,
                self.law.kinetic_coefficient(speed)?,
                "traction_n",
            )?;
            let residual_n = (magnitude - expected).abs();
            let tolerance_n = input.solve.tolerance(magnitude, expected)?;
            if residual_n > tolerance_n {
                return Err(TriboError::TangentialSolve {
                    residual_n,
                    tolerance_n,
                    iterations,
                });
            }
        }
        let work = checked_add(
            checked_mul(traction_n[0], slip_increment_m[0], "tangent_work_j")?,
            checked_mul(traction_n[1], slip_increment_m[1], "tangent_work_j")?,
            "tangent_work_j",
        )?;
        nonnegative_finite(work, "tangent_work_j")?;
        let mut candidate = state.clone();
        candidate.branch = resolved.branch;
        candidate.tangent_axis = if regime == FrictionRegime::Sliding {
            Some(axis)
        } else {
            None
        };
        candidate.prepared = Some(line);
        candidate.roots = resolved.roots;
        candidate.selected = resolved.selected;
        self.finish(
            state,
            input,
            TangentialStep {
                traction_n,
                slip_increment_m,
                state: candidate,
                regime,
                stored_change_j: 0.0,
                dissipated_j: work,
                closure_residual_j: 0.0,
                compatibility_residual_m,
                iterations,
            },
        )
    }
}

impl ObjectiveTangentialLaw for HystereticFriction {
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError> {
        self.require_state(state)?;
        let mut candidate = state.clone();
        candidate.slip_m = rotation.apply(state.slip_m)?;
        candidate.tangent_axis = state
            .tangent_axis
            .map(|axis| rotation.apply(axis))
            .transpose()?;
        Ok(candidate)
    }
}

fn solve_error(
    error: ImplicitFrictionError<TriboError>,
    config: TangentialSolveConfig,
    iterations: usize,
) -> TriboError {
    match error {
        ImplicitFrictionError::Constitutive(error) => error,
        ImplicitFrictionError::InvalidInput(field) => TriboError::InvalidInput { field },
        ImplicitFrictionError::NonDissipative => TriboError::InvalidInput {
            field: "hysteretic_friction_dissipation",
        },
        ImplicitFrictionError::NoConvergence => TriboError::TangentialSolve {
            residual_n: f64::INFINITY,
            tolerance_n: config.force_absolute_tolerance_n,
            iterations,
        },
    }
}

fn finite(value: f64, field: &'static str) -> Result<f64, TriboError> {
    crate::finite(value, field)?;
    Ok(value)
}
