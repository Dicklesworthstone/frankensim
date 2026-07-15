//! The shared admitted-budget accountant (bead frankensim-sj31i.6).
//!
//! Root cause this module closes: workflows hash `Cx::budget()` fields
//! into provenance identity but enforce, at best, a locally-consumed
//! poll quota — deadline and cost quota were never uniformly enforced,
//! so a `Budget::ZERO` invocation could complete successfully while its
//! provenance claimed a zero allowance. This accountant makes the
//! hashed budget and the enforced budget THE SAME CONTRACT: it is
//! admitted once from the ambient `Cx` budget plus the checked work
//! plan's declared cost, and every deterministic tile boundary then
//! passes through [`AdmittedBudget::checkpoint`] /
//! [`AdmittedBudget::charge_cost`].
//!
//! Semantics:
//! - Admission refuses an already-expired deadline (including
//!   `Budget::ZERO`'s `Some(Time::ZERO)`) and a declared cost plan that
//!   exceeds the cost quota, before any work runs.
//! - Checkpoints enforce, in order: cancellation (always takes
//!   precedence), deadline under the caller-supplied deterministic
//!   [`TimeSource`], then poll quota. `u32::MAX` polls and a `None`
//!   deadline/cost quota mean unlimited (asupersync's sentinels).
//! - Refusals latch: after the first refusal every subsequent call
//!   returns the same reason, so exhausted work cannot smuggle a
//!   partial success past its own refusal ("no-partial-success").
//! - [`AdmittedBudget::consumption`] retains the final consumption and
//!   the latched refusal reason for evidence, exactly as they were
//!   enforced.
//!
//! The accountant deliberately has no children, memory leases, or
//! receipt hashing — that is [`crate::InvocationBudget`]'s contract for
//! nested science. This type is the uniform tile-boundary enforcement
//! primitive those workflows were missing.

use asupersync::time::TimeSource;
use asupersync::types::{Budget, Time};

use crate::cx::Cx;

/// Why the accountant refused admission or further work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BudgetRefusal {
    /// The deadline had already passed at admission time.
    DeadlineExpiredAtAdmission {
        /// Admitted absolute deadline, nanoseconds.
        deadline_ns: u64,
        /// Deterministic clock observation at admission, nanoseconds.
        observed_ns: u64,
    },
    /// The budget carries a deadline but the context has no ambient
    /// time source; running un-timed would make the deadline
    /// unenforceable, so admission fails closed.
    DeadlineWithoutClock {
        /// Admitted absolute deadline, nanoseconds.
        deadline_ns: u64,
    },
    /// The checked work plan's declared cost exceeds the cost quota.
    CostPlanExceedsQuota {
        /// Cost the plan declared it needs.
        planned: u64,
        /// The admitted cost quota.
        quota: u64,
    },
    /// Cancellation observed at a checkpoint (always takes precedence).
    Cancelled {
        /// Stable phase name at the observing checkpoint.
        phase: &'static str,
    },
    /// The deadline passed mid-run.
    DeadlineExpired {
        /// Stable phase name at the observing checkpoint.
        phase: &'static str,
        /// Admitted absolute deadline, nanoseconds.
        deadline_ns: u64,
        /// Deterministic clock observation, nanoseconds.
        observed_ns: u64,
    },
    /// The poll quota is exhausted.
    PollsExhausted {
        /// Stable phase name at the observing checkpoint.
        phase: &'static str,
        /// The admitted poll quota.
        quota: u32,
    },
    /// A cost charge exceeds the remaining cost quota.
    CostExhausted {
        /// Stable phase name at the charging boundary.
        phase: &'static str,
        /// Units this charge requested.
        requested: u64,
        /// Units that remained before this charge.
        remaining: u64,
        /// The admitted cost quota.
        quota: u64,
    },
}

impl core::fmt::Display for BudgetRefusal {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DeadlineExpiredAtAdmission {
                deadline_ns,
                observed_ns,
            } => write!(
                formatter,
                "budget refused at admission: deadline {deadline_ns}ns already \
                 passed at {observed_ns}ns"
            ),
            Self::DeadlineWithoutClock { deadline_ns } => write!(
                formatter,
                "budget refused at admission: deadline {deadline_ns}ns cannot be \
                 enforced without an ambient time source"
            ),
            Self::CostPlanExceedsQuota { planned, quota } => write!(
                formatter,
                "budget refused at admission: checked work plan declares \
                 {planned} cost units but the quota admits {quota}"
            ),
            Self::Cancelled { phase } => {
                write!(
                    formatter,
                    "budget stop in '{phase}': cancellation requested"
                )
            }
            Self::DeadlineExpired {
                phase,
                deadline_ns,
                observed_ns,
            } => write!(
                formatter,
                "budget stop in '{phase}': deadline {deadline_ns}ns passed at \
                 {observed_ns}ns"
            ),
            Self::PollsExhausted { phase, quota } => write!(
                formatter,
                "budget stop in '{phase}': poll quota {quota} exhausted"
            ),
            Self::CostExhausted {
                phase,
                requested,
                remaining,
                quota,
            } => write!(
                formatter,
                "budget stop in '{phase}': cost charge of {requested} exceeds \
                 the {remaining} remaining of quota {quota}"
            ),
        }
    }
}

impl core::error::Error for BudgetRefusal {}

/// Final accounting of an admitted budget, retained in evidence.
///
/// Mirrors exactly what was enforced: the admitted contract, what was
/// consumed, and the first (latched) refusal if any. `refusal: None`
/// with work complete is the only state that may report success.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BudgetConsumption {
    /// The admitted budget contract (identical to the hashed one).
    pub admitted: Budget,
    /// Cost units the checked work plan declared at admission.
    pub planned_cost: u64,
    /// Checkpoints consumed against the poll quota.
    pub polls_used: u32,
    /// Cost units actually charged.
    pub cost_charged: u64,
    /// The first refusal, if the accountant ever refused.
    pub refusal: Option<BudgetRefusal>,
}

/// A budget admitted once against the checked work plan and enforced at
/// every deterministic tile boundary thereafter.
pub struct AdmittedBudget<'clock> {
    admitted: Budget,
    clock: Option<&'clock dyn TimeSource>,
    planned_cost: u64,
    polls_used: u32,
    cost_charged: u64,
    refusal: Option<BudgetRefusal>,
}

impl core::fmt::Debug for AdmittedBudget<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("AdmittedBudget")
            .field("planned_cost", &self.planned_cost)
            .field("polls_used", &self.polls_used)
            .field("cost_charged", &self.cost_charged)
            .field("refusal", &self.refusal)
            .finish_non_exhaustive()
    }
}

impl<'clock> AdmittedBudget<'clock> {
    /// Admit the ambient `Cx` budget against the checked work plan's
    /// declared cost, under a deterministic clock.
    ///
    /// # Errors
    /// Refuses an already-expired deadline (a `Budget::ZERO` deadline
    /// can never admit) and a cost plan exceeding the cost quota,
    /// before any work spends authority.
    pub fn admit(
        cx: &Cx<'_>,
        planned_cost: u64,
        clock: &'clock dyn TimeSource,
    ) -> Result<Self, BudgetRefusal> {
        Self::admit_with(cx.budget(), planned_cost, Some(clock))
    }

    /// Admit against the context's ambient time source (bead sj31i.6).
    ///
    /// Keeps existing entry-point signatures intact: the executor
    /// attaches the clock via `Cx::with_time_source`, and workflows
    /// admit from the context alone.
    ///
    /// # Errors
    /// Everything [`AdmittedBudget::admit`] refuses, plus a budget that
    /// carries a deadline on a context with no ambient clock (fail
    /// closed: an unenforceable deadline never admits).
    pub fn admit_ambient(cx: &Cx<'clock>, planned_cost: u64) -> Result<Self, BudgetRefusal> {
        let admitted = cx.budget();
        let clock = cx.time_source();
        if clock.is_none()
            && let Some(deadline) = admitted.deadline
        {
            return Err(BudgetRefusal::DeadlineWithoutClock {
                deadline_ns: deadline.as_nanos(),
            });
        }
        Self::admit_with(admitted, planned_cost, clock)
    }

    fn admit_with(
        admitted: Budget,
        planned_cost: u64,
        clock: Option<&'clock dyn TimeSource>,
    ) -> Result<Self, BudgetRefusal> {
        if let (Some(deadline), Some(clock)) = (admitted.deadline, clock) {
            let now = clock.now();
            if now >= deadline {
                return Err(BudgetRefusal::DeadlineExpiredAtAdmission {
                    deadline_ns: deadline.as_nanos(),
                    observed_ns: now.as_nanos(),
                });
            }
        }
        if let Some(quota) = admitted.cost_quota
            && planned_cost > quota
        {
            return Err(BudgetRefusal::CostPlanExceedsQuota {
                planned: planned_cost,
                quota,
            });
        }
        Ok(Self {
            admitted,
            clock,
            planned_cost,
            polls_used: 0,
            cost_charged: 0,
            refusal: None,
        })
    }

    /// Enforce the admitted contract at a deterministic tile boundary.
    ///
    /// Order: latched refusal, cancellation (precedence), deadline
    /// under the admitted clock, then poll quota. A quota of
    /// `u32::MAX` is unlimited and consumes nothing.
    ///
    /// # Errors
    /// Returns the first refusal and latches it; every later call
    /// returns the same reason.
    pub fn checkpoint(&mut self, phase: &'static str, cx: &Cx<'_>) -> Result<(), BudgetRefusal> {
        if let Some(refusal) = self.refusal {
            return Err(refusal);
        }
        if cx.is_cancel_requested() {
            return Err(self.latch(BudgetRefusal::Cancelled { phase }));
        }
        if let (Some(deadline), Some(clock)) = (self.admitted.deadline, self.clock) {
            let now = clock.now();
            if now >= deadline {
                return Err(self.latch(BudgetRefusal::DeadlineExpired {
                    phase,
                    deadline_ns: deadline.as_nanos(),
                    observed_ns: now.as_nanos(),
                }));
            }
        }
        if self.admitted.poll_quota != u32::MAX {
            if self.polls_used >= self.admitted.poll_quota {
                return Err(self.latch(BudgetRefusal::PollsExhausted {
                    phase,
                    quota: self.admitted.poll_quota,
                }));
            }
            self.polls_used += 1;
        }
        Ok(())
    }

    /// Enforce only cancellation and the deadline, without consuming
    /// poll quota (bead sj31i.6). For workflows whose poll ledger is
    /// shared with nested science and enforced elsewhere.
    ///
    /// # Errors
    /// Returns the first refusal and latches it, like
    /// [`AdmittedBudget::checkpoint`].
    pub fn observe_deadline(
        &mut self,
        phase: &'static str,
        cx: &Cx<'_>,
    ) -> Result<(), BudgetRefusal> {
        if let Some(refusal) = self.refusal {
            return Err(refusal);
        }
        if cx.is_cancel_requested() {
            return Err(self.latch(BudgetRefusal::Cancelled { phase }));
        }
        if let (Some(deadline), Some(clock)) = (self.admitted.deadline, self.clock) {
            let now = clock.now();
            if now >= deadline {
                return Err(self.latch(BudgetRefusal::DeadlineExpired {
                    phase,
                    deadline_ns: deadline.as_nanos(),
                    observed_ns: now.as_nanos(),
                }));
            }
        }
        Ok(())
    }

    /// Convenience: checkpoint then charge one completed tile's cost.
    ///
    /// # Errors
    /// Propagates the first refusal from either step (latched).
    pub fn tile_boundary(
        &mut self,
        phase: &'static str,
        cx: &Cx<'_>,
        tile_cost: u64,
    ) -> Result<(), BudgetRefusal> {
        self.checkpoint(phase, cx)?;
        self.charge_cost(phase, tile_cost)
    }

    /// Charge cost consumed by a completed tile against the quota.
    ///
    /// A `None` cost quota is unlimited but still accumulates
    /// `cost_charged` for evidence.
    ///
    /// # Errors
    /// Refuses (and latches) a charge exceeding the remaining quota;
    /// the failed charge is not recorded as consumed.
    pub fn charge_cost(&mut self, phase: &'static str, units: u64) -> Result<(), BudgetRefusal> {
        if let Some(refusal) = self.refusal {
            return Err(refusal);
        }
        let Some(charged) = self.cost_charged.checked_add(units) else {
            let quota = self.admitted.cost_quota.unwrap_or(u64::MAX);
            return Err(self.latch(BudgetRefusal::CostExhausted {
                phase,
                requested: units,
                remaining: quota.saturating_sub(self.cost_charged),
                quota,
            }));
        };
        if let Some(quota) = self.admitted.cost_quota
            && charged > quota
        {
            return Err(self.latch(BudgetRefusal::CostExhausted {
                phase,
                requested: units,
                remaining: quota - self.cost_charged,
                quota,
            }));
        }
        self.cost_charged = charged;
        Ok(())
    }

    /// Whether a refusal has latched (no further work may succeed).
    #[must_use]
    pub fn is_refused(&self) -> bool {
        self.refusal.is_some()
    }

    /// The deterministic clock reading right now, for evidence stamps.
    #[must_use]
    pub fn now(&self) -> Option<Time> {
        self.clock.map(TimeSource::now)
    }

    /// Final consumption and latched refusal, exactly as enforced.
    #[must_use]
    pub fn consumption(&self) -> BudgetConsumption {
        BudgetConsumption {
            admitted: self.admitted,
            planned_cost: self.planned_cost,
            polls_used: self.polls_used,
            cost_charged: self.cost_charged,
            refusal: self.refusal,
        }
    }

    fn latch(&mut self, refusal: BudgetRefusal) -> BudgetRefusal {
        self.refusal = Some(refusal);
        refusal
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cx::{CancelGate, ExecMode, StreamKey};
    use asupersync::time::VirtualClock;

    const STREAM: StreamKey = StreamKey {
        seed: 0x53_4A31_4936,
        kernel_id: 6,
        tile: 0,
        iteration: 0,
    };

    fn with_cx<R>(budget: Budget, gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(gate, arena, STREAM, budget, ExecMode::Deterministic);
            f(&cx)
        })
    }

    fn budget(deadline: Option<u64>, polls: u32, cost: Option<u64>) -> Budget {
        Budget {
            deadline: deadline.map(Time::from_secs),
            poll_quota: polls,
            cost_quota: cost,
            priority: 0,
        }
    }

    #[test]
    fn zero_budget_refuses_at_admission_and_cannot_complete() {
        let gate = CancelGate::new();
        with_cx(Budget::ZERO, &gate, |cx| {
            let clock = VirtualClock::new();
            let refusal =
                AdmittedBudget::admit(cx, 0, &clock).expect_err("Budget::ZERO must never admit");
            assert!(matches!(
                refusal,
                BudgetRefusal::DeadlineExpiredAtAdmission { deadline_ns: 0, .. }
            ));
        });
    }

    #[test]
    fn already_expired_and_mid_run_deadlines_refuse_deterministically() {
        let gate = CancelGate::new();
        with_cx(budget(Some(10), u32::MAX, None), &gate, |cx| {
            let clock = VirtualClock::new();
            clock.advance(10_000_000_000);
            assert!(matches!(
                AdmittedBudget::admit(cx, 0, &clock),
                Err(BudgetRefusal::DeadlineExpiredAtAdmission { .. })
            ));
        });
        with_cx(budget(Some(10), 8, None), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted =
                AdmittedBudget::admit(cx, 0, &clock).expect("future deadline admits");
            admitted
                .checkpoint("warm", cx)
                .expect("before the deadline the boundary passes");
            clock.advance(10_000_000_000);
            let refusal = admitted
                .checkpoint("late", cx)
                .expect_err("mid-run expiry refuses");
            assert!(matches!(
                refusal,
                BudgetRefusal::DeadlineExpired { phase: "late", .. }
            ));
            let consumption = admitted.consumption();
            assert_eq!(consumption.polls_used, 1);
            assert_eq!(consumption.refusal, Some(refusal));
        });
    }

    #[test]
    fn poll_quota_enforces_exact_boundary_and_one_below() {
        let gate = CancelGate::new();
        with_cx(budget(None, 3, None), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted = AdmittedBudget::admit(cx, 0, &clock).expect("admits");
            for _ in 0..3 {
                admitted
                    .checkpoint("tile", cx)
                    .expect("quota admits exactly its declared boundaries");
            }
            let refusal = admitted
                .checkpoint("tile", cx)
                .expect_err("one past the quota refuses");
            assert_eq!(
                refusal,
                BudgetRefusal::PollsExhausted {
                    phase: "tile",
                    quota: 3
                }
            );
            assert_eq!(admitted.consumption().polls_used, 3);
        });
        with_cx(budget(None, 0, None), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted =
                AdmittedBudget::admit(cx, 0, &clock).expect("zero polls still admits");
            assert!(matches!(
                admitted.checkpoint("first", cx),
                Err(BudgetRefusal::PollsExhausted { quota: 0, .. })
            ));
        });
    }

    #[test]
    fn unlimited_polls_consume_nothing() {
        let gate = CancelGate::new();
        with_cx(budget(None, u32::MAX, None), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted = AdmittedBudget::admit(cx, 0, &clock).expect("admits");
            for _ in 0..1000 {
                admitted.checkpoint("tile", cx).expect("unlimited");
            }
            assert_eq!(admitted.consumption().polls_used, 0);
        });
    }

    #[test]
    fn cost_plan_admission_and_exhaustion_are_exact() {
        let gate = CancelGate::new();
        with_cx(budget(None, u32::MAX, Some(10)), &gate, |cx| {
            let clock = VirtualClock::new();
            assert_eq!(
                AdmittedBudget::admit(cx, 11, &clock).expect_err("over-plan refuses"),
                BudgetRefusal::CostPlanExceedsQuota {
                    planned: 11,
                    quota: 10
                }
            );
            let mut admitted = AdmittedBudget::admit(cx, 10, &clock).expect("exact plan admits");
            admitted.charge_cost("tile-a", 9).expect("within quota");
            admitted.charge_cost("tile-b", 1).expect("exact boundary");
            let refusal = admitted
                .charge_cost("tile-c", 1)
                .expect_err("exhausted quota refuses");
            assert_eq!(
                refusal,
                BudgetRefusal::CostExhausted {
                    phase: "tile-c",
                    requested: 1,
                    remaining: 0,
                    quota: 10
                }
            );
            let consumption = admitted.consumption();
            assert_eq!(consumption.cost_charged, 10);
            assert_eq!(consumption.refusal, Some(refusal));
        });
    }

    #[test]
    fn cancellation_takes_precedence_over_every_other_refusal() {
        let gate = CancelGate::new();
        with_cx(budget(Some(1), 0, Some(0)), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted = AdmittedBudget::admit(cx, 0, &clock).expect("admits");
            gate.request();
            clock.advance(5_000_000_000);
            assert_eq!(
                admitted
                    .checkpoint("cancelled", cx)
                    .expect_err("cancellation refuses"),
                BudgetRefusal::Cancelled { phase: "cancelled" }
            );
        });
    }

    #[test]
    fn refusals_latch_and_forbid_partial_success() {
        let gate = CancelGate::new();
        with_cx(budget(None, 1, Some(5)), &gate, |cx| {
            let clock = VirtualClock::new();
            let mut admitted = AdmittedBudget::admit(cx, 5, &clock).expect("admits");
            admitted.checkpoint("tile", cx).expect("first boundary");
            let first = admitted
                .checkpoint("tile", cx)
                .expect_err("quota exhausted");
            assert!(admitted.is_refused());
            assert_eq!(admitted.charge_cost("tile", 1).expect_err("latched"), first);
            assert_eq!(
                admitted.checkpoint("later", cx).expect_err("latched"),
                first
            );
            assert_eq!(admitted.consumption().cost_charged, 0);
            assert_eq!(admitted.consumption().refusal, Some(first));
        });
    }

    #[test]
    fn ambient_admission_fails_closed_on_deadline_without_clock() {
        let gate = CancelGate::new();
        with_cx(budget(Some(10), u32::MAX, None), &gate, |cx| {
            assert!(matches!(
                AdmittedBudget::admit_ambient(cx, 0),
                Err(BudgetRefusal::DeadlineWithoutClock { .. })
            ));
        });
        with_cx(budget(None, 4, Some(7)), &gate, |cx| {
            let mut admitted = AdmittedBudget::admit_ambient(cx, 7)
                .expect("deadline-free budgets admit without a clock");
            admitted.checkpoint("tile", cx).expect("boundary");
            admitted.charge_cost("tile", 7).expect("charge");
            assert_eq!(admitted.consumption().cost_charged, 7);
            assert_eq!(admitted.now(), None);
        });
    }

    #[test]
    fn ambient_clock_enforces_the_deadline() {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let clock = VirtualClock::new();
            let cx = Cx::new(
                &gate,
                arena,
                STREAM,
                budget(Some(10), 8, None),
                ExecMode::Deterministic,
            )
            .with_time_source(&clock);
            let mut admitted =
                AdmittedBudget::admit_ambient(&cx, 0).expect("future deadline admits");
            admitted.checkpoint("warm", &cx).expect("in time");
            clock.advance(10_000_000_000);
            assert!(matches!(
                admitted.checkpoint("late", &cx),
                Err(BudgetRefusal::DeadlineExpired { phase: "late", .. })
            ));
        });
    }

    #[test]
    fn replaying_the_same_plan_yields_identical_consumption() {
        let gate = CancelGate::new();
        let run = || {
            with_cx(budget(Some(100), 4, Some(7)), &gate, |cx| {
                let clock = VirtualClock::new();
                let mut admitted = AdmittedBudget::admit(cx, 7, &clock).expect("admits");
                for _ in 0..3 {
                    admitted.checkpoint("tile", cx).expect("boundary");
                }
                admitted.charge_cost("tile", 7).expect("charge");
                admitted.consumption()
            })
        };
        assert_eq!(run(), run());
    }
}
