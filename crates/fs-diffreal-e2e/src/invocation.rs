//! Fixed, typed invocation plan for the as-built scientific transaction.

use super::*;

pub(super) const AS_BUILT_INVOCATION_POLICY: &str = "fs-diffreal-e2e/as-built-invocation-budget/v1";

const SETUP_WORK: u128 = 128;
const SETUP_POLLS: u32 = 2;
const SETUP_MEMORY: u64 = 4 * 1024;
const ASSIMILATION_POLLS: u32 = 292;
// Fixed fixture: initial/final (2) + 682 record units / 16 (42) +
// 3,051 deterministic hash bytes / 1,024 (2). Fast mode hashes nine fewer
// bytes and therefore has the same exact poll count.
const ASSIMILATION_EXPECTED_POLLS: u32 = 46;
const ASSIMILATION_MEMORY: u64 = 64 * 1024;
const ASSIMILATION_OUTPUT: u64 = 16 * 1024;
const PUBLICATION_WORK: u128 = 8 * 1024;
const PUBLICATION_POLLS: u32 = 2;
const PUBLICATION_MEMORY: u64 = 32 * 1024;
const PUBLICATION_OUTPUT: u64 = 32 * 1024;

#[derive(Debug, Clone)]
pub(super) struct AsBuiltInvocationPlan {
    pub(super) setup: fs_exec::InvocationResources,
    pub(super) registration: fs_exec::InvocationResources,
    pub(super) difference: fs_exec::InvocationResources,
    pub(super) belief: fs_exec::InvocationResources,
    pub(super) assimilation: fs_exec::InvocationResources,
    pub(super) assimilation_shape: fs_exec::InvocationResources,
    pub(super) publication: fs_exec::InvocationResources,
    pub(super) total: fs_exec::InvocationResources,
}

impl AsBuiltInvocationPlan {
    pub(super) fn preflight(execution: &DiffRealExecutionIdentity) -> Result<Self, DiffRealError> {
        let setup = resources(SETUP_WORK, SETUP_POLLS, 0, SETUP_MEMORY, 0)?;
        let registration =
            fs_asbuilt::registration_invocation_resources(AS_BUILT_DESIGN_POINTS.len())?;
        let difference = fs_asbuilt::as_built_diff_invocation_resources(
            AS_BUILT_DESIGN_POINTS.len(),
            AS_BUILT_DESIGN_POINTS.len(),
            AS_BUILT_DESIGN_TOLERANCE,
            AS_BUILT_MEASUREMENT_NOISE,
            AS_BUILT_CALIBRATION_CANDIDATE,
        )?;
        let belief =
            fs_assimilate::diagonal_belief_invocation_resources(AS_BUILT_PRIOR_MEAN.len())?;
        let observations = as_built_observations()?;
        let assimilation_shape =
            fs_assimilate::colored_assimilation_invocation_resources_for_shape(
                AS_BUILT_PRIOR_MEAN.len(),
                &observations,
                AS_BUILT_ASSIMILATION_PARAMETER,
                AS_BUILT_ASSIMILATION_BOUNDS.0,
                AS_BUILT_ASSIMILATION_BOUNDS.1,
                execution.mode(),
            )?;
        let assimilation = resources(
            assimilation_shape.work().get(),
            ASSIMILATION_POLLS,
            assimilation_shape.evaluations().get(),
            ASSIMILATION_MEMORY,
            ASSIMILATION_OUTPUT,
        )?;
        assimilation.checked_sub(assimilation_shape)?;
        let publication = resources(
            PUBLICATION_WORK,
            PUBLICATION_POLLS,
            0,
            PUBLICATION_MEMORY,
            PUBLICATION_OUTPUT,
        )?;
        let phases = [
            setup,
            registration,
            difference,
            belief,
            assimilation,
            publication,
        ];
        let total = transaction_total(setup, &phases)?;
        Ok(Self {
            setup,
            registration,
            difference,
            belief,
            assimilation,
            assimilation_shape,
            publication,
            total,
        })
    }

    pub(super) fn limits(
        &self,
        execution: &DiffRealExecutionIdentity,
    ) -> fs_exec::InvocationLimits {
        let ambient = execution.budget();
        let poll_limit = ambient.poll_quota;
        let cost_limit = ambient.cost_quota.unwrap_or(self.total.cost().get());
        let memory_limit = execution
            .operation_memory_limit_bytes()
            .unwrap_or(self.total.memory().get());
        let limits = fs_exec::InvocationResources::new(
            self.total.work(),
            fs_exec::PollUnits::new(poll_limit),
            fs_exec::CostUnits::new(cost_limit),
            self.total.evaluations(),
            fs_exec::MemoryBytes::new(memory_limit),
            self.total.output(),
        );
        fs_exec::InvocationLimits::new(
            limits,
            ambient.deadline,
            accuracy_obligation(),
            capability_scope(),
        )
    }

    pub(super) fn verifies_receipt(
        &self,
        receipt: &fs_exec::InvocationReceipt,
        execution: &DiffRealExecutionIdentity,
    ) -> bool {
        if receipt.version() != fs_exec::INVOCATION_RECEIPT_VERSION
            || receipt.invocation_id() != invocation_id(execution)
            || receipt.limits() != &self.limits(execution)
            || receipt.required() != self.total
            || receipt.disposition() != fs_exec::InvocationDisposition::Completed
            || !receipt.verifies_integrity()
        {
            return false;
        }
        let children = receipt.children();
        if children.len() != 6 {
            return false;
        }
        let transaction = children[0].id();
        if children[0].ordinal() != 0
            || children[0].phase() != "as-built.transaction"
            || children[0].granted() != self.total
            || children[0].parent().is_some()
            || children[0].disposition() != fs_exec::InvocationDisposition::Completed
        {
            return false;
        }
        let leaves = [
            (
                "as-built.registration",
                self.registration,
                expected_leaf(self.registration, self.registration.polls().get()),
            ),
            (
                "as-built.comparison",
                self.difference,
                expected_leaf(self.difference, self.difference.polls().get()),
            ),
            (
                "as-built.prior",
                self.belief,
                expected_leaf(self.belief, self.belief.polls().get()),
            ),
            (
                "as-built.assimilation",
                self.assimilation,
                expected_leaf(self.assimilation_shape, ASSIMILATION_EXPECTED_POLLS),
            ),
            (
                "as-built.publication",
                self.publication,
                expected_leaf(self.publication, self.publication.polls().get()),
            ),
        ];
        if !children[1..].iter().zip(leaves).enumerate().all(
            |(index, (child, (phase, grant, evidence)))| {
                child.ordinal() == (index + 1) as u64
                    && child.phase() == phase
                    && child.granted() == grant
                    && child.parent() == Some(transaction)
                    && child.disposition() == fs_exec::InvocationDisposition::Completed
                    && leaf_evidence_matches(child, evidence)
            },
        ) {
            return false;
        }

        let setup = expected_leaf(self.setup, self.setup.polls().get());
        let mut transaction_consumed = setup.consumed;
        let mut requested_memory = setup.memory;
        let mut retained_output = 0_u64;
        let mut nested_memory_peak = 0_u64;
        for (_, _, evidence) in leaves {
            let Ok(consumed) = transaction_consumed.checked_add(evidence.consumed) else {
                return false;
            };
            transaction_consumed = consumed;
            let Some(requested) = requested_memory.checked_add(evidence.memory) else {
                return false;
            };
            requested_memory = requested;
            let Some(output) = retained_output.checked_add(evidence.consumed.output().get()) else {
                return false;
            };
            retained_output = output;
            nested_memory_peak = nested_memory_peak.max(evidence.memory);
        }
        let Some(transaction_memory_peak) = setup.memory.checked_add(nested_memory_peak) else {
            return false;
        };
        let Ok(expected_remaining) = self.total.checked_sub(transaction_consumed) else {
            return false;
        };
        children[0].consumed() == transaction_consumed
            && children[0].direct_consumed() == setup.consumed
            && children[0].direct_memory_peak_bytes() == setup.memory
            && children[0].memory_peak_bytes() == transaction_memory_peak
            && children[0].memory_requested_bytes() == setup.memory as u128
            && children[0].memory_released_bytes() == setup.memory as u128
            && children[0].output_retained_bytes() == 0
            && receipt.remaining() == expected_remaining
            && receipt.memory_peak_bytes() == transaction_memory_peak
            && receipt.memory_requested_bytes() == requested_memory as u128
            && receipt.memory_released_bytes() == requested_memory as u128
            && receipt.output_retained_bytes() == retained_output
    }
}

#[derive(Clone, Copy)]
struct ExpectedLeafEvidence {
    consumed: fs_exec::InvocationResources,
    memory: u64,
}

fn expected_leaf(shape: fs_exec::InvocationResources, polls: u32) -> ExpectedLeafEvidence {
    ExpectedLeafEvidence {
        consumed: fs_exec::InvocationResources::new(
            shape.work(),
            fs_exec::PollUnits::new(polls),
            shape.cost(),
            shape.evaluations(),
            fs_exec::MemoryBytes::new(0),
            shape.output(),
        ),
        memory: shape.memory().get(),
    }
}

fn leaf_evidence_matches(child: &fs_exec::ChildReceipt, expected: ExpectedLeafEvidence) -> bool {
    child.consumed() == expected.consumed
        && child.direct_consumed() == expected.consumed
        && child.direct_memory_peak_bytes() == expected.memory
        && child.memory_peak_bytes() == expected.memory
        && child.memory_requested_bytes() == expected.memory as u128
        && child.memory_released_bytes() == expected.memory as u128
        && child.output_retained_bytes() == expected.consumed.output().get()
}

pub(super) fn invocation_id(execution: &DiffRealExecutionIdentity) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_str(&mut canonical, "policy", AS_BUILT_INVOCATION_POLICY);
    push_identity_str(&mut canonical, "stage", AS_BUILT_STAGE);
    push_identity_str(
        &mut canonical,
        "evidence-identity",
        AS_BUILT_EVIDENCE_IDENTITY,
    );
    encode_execution_identity(&mut canonical, execution);
    hash_domain(
        "frankensim.fs-diffreal-e2e.as-built-invocation.v1",
        &canonical,
    )
}

pub(super) fn domain_refusal(stage: &'static str, detail: &str) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_str(&mut canonical, "stage", stage);
    push_identity_str(&mut canonical, "detail", detail);
    hash_domain(
        "frankensim.fs-diffreal-e2e.invocation-domain-refusal.v1",
        &canonical,
    )
}

fn resources(
    work: u128,
    polls: u32,
    evaluations: u64,
    memory: u64,
    output: u64,
) -> Result<fs_exec::InvocationResources, DiffRealError> {
    let cost = u64::try_from(work)
        .map_err(|_| fs_exec::InvocationError::ArithmeticOverflow { resource: "cost" })?;
    Ok(fs_exec::InvocationResources::new(
        fs_exec::WorkUnits::new(work),
        fs_exec::PollUnits::new(polls),
        fs_exec::CostUnits::new(cost),
        fs_exec::EvaluationUnits::new(evaluations),
        fs_exec::MemoryBytes::new(memory),
        fs_exec::OutputBytes::new(output),
    ))
}

fn transaction_total(
    setup: fs_exec::InvocationResources,
    phases: &[fs_exec::InvocationResources],
) -> Result<fs_exec::InvocationResources, DiffRealError> {
    let mut work = 0_u128;
    let mut polls = 0_u32;
    let mut cost = 0_u64;
    let mut evaluations = 0_u64;
    let mut output = 0_u64;
    let mut nested_memory = 0_u64;
    for (index, phase) in phases.iter().enumerate() {
        work = work
            .checked_add(phase.work().get())
            .ok_or(fs_exec::InvocationError::ArithmeticOverflow { resource: "work" })?;
        polls = polls
            .checked_add(phase.polls().get())
            .ok_or(fs_exec::InvocationError::ArithmeticOverflow { resource: "polls" })?;
        cost = cost
            .checked_add(phase.cost().get())
            .ok_or(fs_exec::InvocationError::ArithmeticOverflow { resource: "cost" })?;
        evaluations = evaluations.checked_add(phase.evaluations().get()).ok_or(
            fs_exec::InvocationError::ArithmeticOverflow {
                resource: "evaluations",
            },
        )?;
        output = output.checked_add(phase.output().get()).ok_or(
            fs_exec::InvocationError::ArithmeticOverflow {
                resource: "output-bytes",
            },
        )?;
        if index != 0 {
            nested_memory = nested_memory.max(phase.memory().get());
        }
    }
    let memory = setup.memory().get().checked_add(nested_memory).ok_or(
        fs_exec::InvocationError::ArithmeticOverflow {
            resource: "memory-bytes",
        },
    )?;
    Ok(fs_exec::InvocationResources::new(
        fs_exec::WorkUnits::new(work),
        fs_exec::PollUnits::new(polls),
        fs_exec::CostUnits::new(cost),
        fs_exec::EvaluationUnits::new(evaluations),
        fs_exec::MemoryBytes::new(memory),
        fs_exec::OutputBytes::new(output),
    ))
}

fn accuracy_obligation() -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_f64(
        &mut canonical,
        "design-tolerance",
        AS_BUILT_DESIGN_TOLERANCE,
    );
    push_identity_f64(
        &mut canonical,
        "measurement-noise",
        AS_BUILT_MEASUREMENT_NOISE,
    );
    push_identity_f64(
        &mut canonical,
        "assimilation-lower",
        AS_BUILT_ASSIMILATION_BOUNDS.0,
    );
    push_identity_f64(
        &mut canonical,
        "assimilation-upper",
        AS_BUILT_ASSIMILATION_BOUNDS.1,
    );
    hash_domain(
        "frankensim.fs-diffreal-e2e.as-built-accuracy-obligation.v1",
        &canonical,
    )
}

fn capability_scope() -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_str(&mut canonical, "stage", AS_BUILT_STAGE);
    push_identity_str(&mut canonical, "policy", AS_BUILT_INVOCATION_POLICY);
    push_identity_str(
        &mut canonical,
        "as-built-calibration",
        AS_BUILT_CALIBRATION_CANDIDATE,
    );
    hash_domain(
        "frankensim.fs-diffreal-e2e.as-built-capability-scope.v1",
        &canonical,
    )
}
