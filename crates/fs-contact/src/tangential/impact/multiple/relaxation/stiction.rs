//! Coupled breakaway and sliding friction on the original affine flow owner.
//!
//! Moving ports carry the kinetic force. At an exact zero-velocity state, a
//! holding port can supply any force within its larger static interval. The
//! declared continuation preserves maximal admissible sticking sets by set
//! inclusion: a mode is discarded if another admissible mode keeps every one
//! of its held ports and at least one additional port at rest. Incomparable
//! surviving modes refuse as ambiguous; their enumeration order is not physics.
//!
//! This is an explicit hold-priority convention for a discontinuous model,
//! not a general uniqueness theorem or a velocity regularization. A nonzero
//! incoming velocity slides, however small. Compliance, microscopic junction
//! history, changing loads and material selection remain separate owners.

use super::{
    Accurate, CoupledFrictionFlowControls, CoupledFrictionFlowError, CoupledFrictionFlowInput,
    CoupledFrictionFlowReport, CoupledFrictionImpulseError, CoupledFrictionImpulseInput, Dd,
    Factor, MAX_FRICTION_IMPULSE_PORTS, Path, PhaseForces, PhaseLaw, Vector,
    advance_coupled_friction, advance_with_law, force_budget, poll, response_row,
};

/// Bounded small joint/contact islands; the existing equal-bound owner still
/// admits its original thirty-two ports.
pub const MAX_STICTION_FLOW_PORTS: usize = 6;
/// All hold/positive-release/negative-release assignments for six resting ports.
pub const MAX_STICTION_MODES_PER_PHASE: usize = 729;

/// Work-conjugate force data, held over one physical phase.
#[derive(Clone, Copy, Debug)]
pub struct CoupledStictionFlowInput<'a> {
    pub velocity_before_m_s: &'a [f64],
    /// Full row-major positive-definite mobility, never its diagonal surrogate.
    pub mobility_per_kg: &'a [f64],
    pub static_capacity_n: &'a [f64],
    pub kinetic_capacity_n: &'a [f64],
    pub drive_n: &'a [f64],
    pub duration_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledStictionFlowControls {
    /// Original event, matrix, acceleration, impulse, work and resource gates.
    pub flow: CoupledFrictionFlowControls,
    /// Total candidate balances across the complete path, not per segment.
    /// At most `729 * flow.maximum_segments` can be needed for six ports.
    pub maximum_mode_evaluations: usize,
}

/// Bounded diagnostics separate from the original continuous path receipt.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledStictionFlowDiagnostics {
    pub mode_evaluations: usize,
    /// Original principal mobility factorizations, reused across release signs.
    pub mode_factorizations: usize,
    /// Initially resting positive-capacity ports leaving rest in resolved phases.
    pub breakaways: usize,
    /// Held coordinates during the last integrated affine segment. An arrest
    /// exactly at the endpoint does not perform a speculative next-phase solve.
    pub last_phase_holding_mask: u32,
    /// Exact zero coordinates of the resolved path endpoint, including arrests.
    pub endpoint_resting_mask: u32,
    /// Largest actual force-balance residual on a selected held coordinate.
    pub maximum_holding_acceleration_residual_m_s2: f64,
    /// Equal bounds delegate to the original owner, without mode enumeration.
    /// In that case only the endpoint mask is supplied by this diagnostic.
    pub used_equal_bound_solver: bool,
}

#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct CoupledStictionFlowReport {
    /// The original endpoint, integrated impulses, path work and speed extrema.
    /// Its integrated-capacity residual uses the static bound times duration.
    pub flow: CoupledFrictionFlowReport,
    pub diagnostics: CoupledStictionFlowDiagnostics,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CoupledStictionFlowError {
    Invalid(&'static str),
    Flow(CoupledFrictionFlowError),
    ModeBudget {
        evaluations: usize,
        maximum: usize,
    },
    NoAdmissibleMode {
        resting_mask: u32,
    },
    /// Two already maximal witnesses suffice to refuse. `maximal_modes` is
    /// the witnessed count, not an exhaustive count of all possible branches.
    Ambiguous {
        resting_mask: u32,
        maximal_modes: usize,
        first_holding_mask: u32,
        second_holding_mask: u32,
    },
}
impl From<CoupledFrictionFlowError> for CoupledStictionFlowError {
    fn from(error: CoupledFrictionFlowError) -> Self {
        Self::Flow(error)
    }
}
impl From<CoupledFrictionImpulseError> for CoupledStictionFlowError {
    fn from(error: CoupledFrictionImpulseError) -> Self {
        Self::Flow(error.into())
    }
}
impl core::fmt::Display for CoupledStictionFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled static/kinetic friction: {self:?}")
    }
}
impl std::error::Error for CoupledStictionFlowError {}

/// Resolve distinct breakaway/sliding limits without an artificial capture speed.
///
/// Static bounds must be at least the nonnegative kinetic bounds. Every actual
/// moving port uses its kinetic force until an isolated zero event from the
/// original flow owner. All resting coordinates first attempt to hold together.
/// If they cannot, original principal factors test every undominated holding
/// set and release direction, subject to the explicit work budget. Strictly
/// nonzero acceleration must agree with each released port's slip direction.
///
/// The largest sticking set is not selected by cardinality: two incomparable
/// sets, even with different sizes, refuse. Static forces are not clipped or
/// admitted above their original capacities. A partial balance whose small
/// computed capacity overshoot is unresolved at the original acceleration
/// allowance refuses instead of turning numerical error into breakaway.
///
/// Endpoint momentum response, continuous work, speed extrema and all original
/// force/impulse/energy budgets remain the shared owner's. Its floating impulse
/// endpoint and exact event path endpoint remain distinct. An external owner
/// must preserve or explicitly project its own numerical resting constraints;
/// this function never infers stiction history from a small velocity.
pub fn advance_coupled_stiction(
    input: &CoupledStictionFlowInput<'_>,
    controls: CoupledStictionFlowControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CoupledStictionFlowReport, CoupledStictionFlowError> {
    poll(&mut cancelled)?;
    let n = input.velocity_before_m_s.len();
    if n > MAX_STICTION_FLOW_PORTS
        || input.static_capacity_n.len() != n
        || input.kinetic_capacity_n.len() != n
        || controls.maximum_mode_evaluations == 0
    {
        return Err(CoupledStictionFlowError::Invalid(
            "stiction dimensions or mode budget",
        ));
    }
    for (&static_n, &kinetic_n) in input.static_capacity_n.iter().zip(input.kinetic_capacity_n) {
        poll(&mut cancelled)?;
        if !static_n.is_finite()
            || !kinetic_n.is_finite()
            || kinetic_n < 0.0
            || static_n < kinetic_n
        {
            return Err(CoupledStictionFlowError::Invalid(
                "finite static >= kinetic >= zero bounds",
            ));
        }
    }
    let original = CoupledFrictionFlowInput {
        velocity_before_m_s: input.velocity_before_m_s,
        mobility_per_kg: input.mobility_per_kg,
        capacity_n: input.static_capacity_n,
        drive_n: input.drive_n,
        duration_s: input.duration_s,
    };
    let equal = input
        .static_capacity_n
        .iter()
        .zip(input.kinetic_capacity_n)
        .all(|(&a, &b)| a == b);
    let (flow, mut diagnostics) = if equal {
        (
            advance_coupled_friction(&original, controls.flow, &mut cancelled)?,
            CoupledStictionFlowDiagnostics {
                used_equal_bound_solver: true,
                ..Default::default()
            },
        )
    } else {
        let mut law = StictionPhases {
            kinetic: input.kinetic_capacity_n,
            maximum_evaluations: controls.maximum_mode_evaluations,
            diagnostics: CoupledStictionFlowDiagnostics::default(),
        };
        let flow = advance_with_law(&original, controls.flow, &mut law, &mut cancelled)?;
        (flow, law.diagnostics)
    };
    for i in 0..n {
        if flow.flow_velocity_after_m_s[i] == 0.0 {
            diagnostics.endpoint_resting_mask |= 1_u32 << i;
        }
    }
    poll(&mut cancelled)?;
    Ok(CoupledStictionFlowReport { flow, diagnostics })
}

struct StictionPhases<'a> {
    kinetic: &'a [f64],
    maximum_evaluations: usize,
    diagnostics: CoupledStictionFlowDiagnostics,
}
#[derive(Clone, Copy)]
struct Candidate {
    force: Vector,
    held: u32,
    maximum_residual: f64,
}
impl StictionPhases<'_> {
    fn available(&self) -> Result<(), CoupledStictionFlowError> {
        if self.diagnostics.mode_evaluations >= self.maximum_evaluations {
            Err(CoupledStictionFlowError::ModeBudget {
                evaluations: self.diagnostics.mode_evaluations,
                maximum: self.maximum_evaluations,
            })
        } else {
            Ok(())
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn candidate<C: FnMut() -> bool>(
        &mut self,
        input: &CoupledFrictionFlowInput<'_>,
        controls: CoupledFrictionFlowControls,
        scales: &Vector,
        fixed: &Vector,
        resting: u32,
        held: u32,
        positive_release: u32,
        factor: &Factor,
        cancelled: &mut C,
    ) -> Result<Option<Candidate>, CoupledStictionFlowError> {
        poll(cancelled)?;
        self.available()?;
        self.diagnostics.mode_evaluations += 1;
        let n = input.velocity_before_m_s.len();
        let all = (1_u32 << n) - 1;
        let mut force = *fixed;
        let closed_form = held == all;
        if closed_form {
            // W is already admitted SPD. Full rest therefore has the exact
            // force solution -drive, including the original static boundary.
            for (f, &drive) in force.iter_mut().zip(input.drive_n) {
                *f = -drive;
            }
        } else {
            for i in 0..n {
                if resting & !held & (1_u32 << i) != 0 {
                    force[i] = if positive_release & (1_u32 << i) != 0 {
                        -self.kinetic[i]
                    } else {
                        self.kinetic[i]
                    };
                }
            }
            let mut rhs = [0.0; MAX_FRICTION_IMPULSE_PORTS];
            for i in 0..n {
                if held & (1_u32 << i) != 0 {
                    rhs[i] = -response_row(input, i, &force)?.to_f64();
                }
            }
            let solved = factor.solve(&rhs, scales, cancelled)?;
            for i in 0..n {
                if held & (1_u32 << i) != 0 {
                    force[i] = solved[i];
                }
            }
        }
        for i in 0..n {
            if held & (1_u32 << i) == 0 {
                continue;
            }
            let excess = force[i].abs() - input.capacity_n[i];
            if excess > 0.0 {
                // This is a refusal test only; the actual force is never
                // replaced with a clipped force or a widened static bound.
                let acceleration_scale = (Dd::from_f64(excess)
                    * Dd::from_f64(input.mobility_per_kg[i * n + i]))
                .to_f64();
                if !closed_form && acceleration_scale <= controls.acceleration_tolerance_m_s2 {
                    return Err(CoupledFrictionFlowError::Unresolved(
                        "static capacity crossing needs a resolved force balance",
                    )
                    .into());
                }
                return Ok(None);
            }
        }
        let mut maximum_residual = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            if resting & (1_u32 << i) == 0 {
                continue;
            }
            let acceleration = response_row(input, i, &force)?;
            let value = acceleration.to_f64();
            if held & (1_u32 << i) != 0 {
                maximum_residual = maximum_residual.max(value.abs());
                if value.abs() > controls.acceleration_tolerance_m_s2 {
                    return Err(CoupledFrictionFlowError::Unresolved(
                        "original static force balance",
                    )
                    .into());
                }
            } else {
                if acceleration != Dd::ZERO && value == 0.0 {
                    return Err(CoupledFrictionFlowError::Unresolved(
                        "representable release acceleration",
                    )
                    .into());
                }
                let admissible = if positive_release & (1_u32 << i) != 0 {
                    value > 0.0
                } else {
                    value < 0.0
                };
                if !admissible {
                    return Ok(None);
                }
            }
        }
        Ok(Some(Candidate {
            force,
            held,
            maximum_residual,
        }))
    }
}

impl PhaseLaw for StictionPhases<'_> {
    type Error = CoupledStictionFlowError;

    #[allow(clippy::too_many_lines)]
    fn resolve<C: FnMut() -> bool>(
        &mut self,
        input: &CoupledFrictionFlowInput<'_>,
        controls: CoupledFrictionFlowControls,
        velocity: &Accurate,
        scales: &Vector,
        _remaining: Dd,
        sums: &mut Path,
        cancelled: &mut C,
    ) -> Result<PhaseForces, Self::Error> {
        let n = input.velocity_before_m_s.len();
        let all = (1_u32 << n) - 1;
        let mut resting = 0_u32;
        let mut fixed = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        for i in 0..n {
            poll(cancelled)?;
            if velocity[i] == Dd::ZERO {
                resting |= 1_u32 << i;
            } else {
                let speed = velocity[i].to_f64();
                if speed == 0.0 {
                    return Err(CoupledFrictionFlowError::Unresolved(
                        "representable moving direction",
                    )
                    .into());
                }
                fixed[i] = if speed > 0.0 {
                    -self.kinetic[i]
                } else {
                    self.kinetic[i]
                };
            }
        }
        let zeros = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        let original = CoupledFrictionImpulseInput {
            velocity_before_m_s: input.velocity_before_m_s,
            mobility_per_kg: input.mobility_per_kg,
            capacity_ns: &zeros[..n],
        };
        let mut factor = Factor::empty();
        let mut selected: Option<Candidate> = None;
        // Descending cardinality is only search order. Every undominated set
        // is tested; an incomparable smaller set is not discarded by its size.
        'sets: for size in (0..=resting.count_ones()).rev() {
            let mut held = resting;
            loop {
                poll(cancelled)?;
                let dominated = selected.is_some_and(|chosen| held & chosen.held == held);
                if held.count_ones() == size && !dominated {
                    self.available()?;
                    if held != 0 && held != all {
                        factor.prepare(&original, scales, held, controls.projection, cancelled)?;
                        self.diagnostics.mode_factorizations += 1;
                        sums.factorizations += 1;
                    } else if held == 0 {
                        // A previously prepared factor may belong to another
                        // set; an empty solve must not reuse its old rows.
                        factor.prepare(&original, scales, 0, controls.projection, cancelled)?;
                    }
                    let released = resting & !held;
                    let mut indices = [0usize; MAX_STICTION_FLOW_PORTS];
                    let mut count = 0;
                    for i in 0..n {
                        if released & (1_u32 << i) != 0 {
                            indices[count] = i;
                            count += 1;
                        }
                    }
                    for signs in 0..(1usize << count) {
                        let mut positive = 0_u32;
                        for (i, &port) in indices[..count].iter().enumerate() {
                            if signs & (1usize << i) != 0 {
                                positive |= 1_u32 << port;
                            }
                        }
                        if let Some(candidate) = self.candidate(
                            input, controls, scales, &fixed, resting, held, positive, &factor,
                            cancelled,
                        )? {
                            if let Some(previous) = selected {
                                return Err(CoupledStictionFlowError::Ambiguous {
                                    resting_mask: resting,
                                    maximal_modes: 2,
                                    first_holding_mask: previous.held,
                                    second_holding_mask: candidate.held,
                                });
                            }
                            selected = Some(candidate);
                            if held == resting {
                                break 'sets;
                            }
                        }
                    }
                }
                if held == 0 {
                    break;
                }
                held = (held - 1) & resting;
            }
        }
        let selected = selected.ok_or(CoupledStictionFlowError::NoAdmissibleMode {
            resting_mask: resting,
        })?;
        // Budget admission follows physical selection. A different branch is
        // never chosen merely because its forces fit a computational limit.
        for i in 0..n {
            poll(cancelled)?;
            force_budget(i, selected.force[i], controls)?;
            let total = super::checked(
                Dd::from_f64(input.drive_n[i]) + Dd::from_f64(selected.force[i]),
                "stiction resultant",
            )?;
            force_budget(i, total.to_f64(), controls)?;
            if resting & !selected.held & (1_u32 << i) != 0 && input.capacity_n[i] > 0.0 {
                self.diagnostics.breakaways += 1;
            }
        }
        self.diagnostics.last_phase_holding_mask = selected.held;
        self.diagnostics.maximum_holding_acceleration_residual_m_s2 = self
            .diagnostics
            .maximum_holding_acceleration_residual_m_s2
            .max(selected.maximum_residual);
        Ok(PhaseForces {
            friction: selected.force,
            holding_mask: Some(selected.held),
        })
    }
}
