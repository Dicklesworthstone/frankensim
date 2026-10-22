//! Recover learned proposals and the current physical phase without refitting.
use fs_blake3::{ContentHash, hash_domain};
use super::{ImportancePilotRound, ImportanceProposalFit, PilotImportanceConfig,
    PilotImportanceExecution, phase_seed};
use crate::product_checkpoint::execution_identity;
use crate::product_plan::admit_plan;
use crate::{ImportanceConfig, ImportanceExecution, UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSIPL001";
const HEADER: usize = 49;
const FIT_BYTES: usize = 42; // Two counts, three f64 values and two flags.
const MIN_BYTES: usize = HEADER + 8 + 81 + 32;
const DOMAIN: &str = "org.frankensim.uq.pilot-importance.checkpoint.v1";

fn invalid(reason: &'static str) -> UqCheckpointError {
    UqCheckpointError::InvalidEncoding(reason)
}
fn push_u64(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u64).to_le_bytes());
}
fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}
fn identity(plan: &UqPlan, config: &PilotImportanceConfig, model: ContentHash) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(execution_identity(plan, model).as_bytes());
    push_u64(&mut bytes, config.pilot_rounds);
    push_u64(&mut bytes, config.samples_per_pilot);
    push_f64(&mut bytes, config.initial_proposal.defensive_fraction);
    push_u64(&mut bytes, config.initial_proposal.germ_shift.len());
    for &x in &config.initial_proposal.germ_shift { push_f64(&mut bytes, x); }
    for x in [config.fit.elite_fraction, config.fit.relaxation, config.fit.maximum_shift_norm] {
        push_f64(&mut bytes, x);
    }
    // This version binds phase_seed, weighted elite selection and the frozen
    // production transition, in addition to all original plan/model fields.
    hash_domain("org.frankensim.uq.pilot-importance.execution.v1", &bytes)
}
fn wire_limit(plan: &UqPlan, config: &PilotImportanceConfig) -> Result<usize, UqCheckpointError> {
    if !(1..=16).contains(&config.pilot_rounds) {
        return Err(UqCheckpointError::InvalidPlan("importance learning requires 1..=16 pilots"));
    }
    plan.parameters.len().checked_mul(8).and_then(|n| n.checked_add(FIT_BYTES))
        .and_then(|n| n.checked_mul(config.pilot_rounds))
        .and_then(|n| plan.budget_max_samples.checked_mul(16).and_then(|m| n.checked_add(m)))
        .and_then(|n| n.checked_add(MIN_BYTES))
        .ok_or_else(|| invalid("pilot checkpoint size overflows"))
}

impl PilotImportanceExecution {
    /// Save exact learned shifts, fitting diagnostics and the current phase's
    /// original QoIs/weights. Completed pilot responses have already been
    /// reduced and need not be retained or evaluated again. A paid pilot paused
    /// before fitting, and paid production paused before publication, are valid.
    /// Failed studies cannot become resumable prefixes.
    ///
    /// Model identity must cover the physical evaluator and all fixed inputs.
    /// The original plan, fitting policy, phase-key rule and proposal semantics
    /// are also bound. Checksums detect corruption, NOT authenticity or correct
    /// physics/fitting. Restore only trusted bytes on the same runtime profile.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let status = match self.status {
            UqStatus::BudgetTruncated => 0,
            UqStatus::Cancelled => 1,
            UqStatus::Complete => 2,
            UqStatus::Refused => return Err(UqCheckpointError::NotResumable),
        };
        if self.failure.is_some() { return Err(UqCheckpointError::NotResumable); }
        let current = self.current.checkpoint(model)?;
        let size = HEADER + self.completed.len() * (FIT_BYTES + 8 * self.plan.parameters.len())
            + 8 + current.len() + 32;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(size)
            .map_err(|_| invalid("pilot checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(identity(&self.plan, &self.config, model).as_bytes());
        bytes.push(status);
        push_u64(&mut bytes, self.completed.len());
        for round in &self.completed {
            let fit = &round.fit;
            push_u64(&mut bytes, fit.selected);
            push_u64(&mut bytes, fit.observed_failures);
            push_f64(&mut bytes, fit.cutoff);
            bytes.push(u8::from(fit.targets_failure_event));
            push_f64(&mut bytes, fit.effective_samples);
            push_f64(&mut bytes, fit.norm_before_limit);
            bytes.push(u8::from(fit.norm_limited));
            for &x in &fit.proposal.germ_shift { push_f64(&mut bytes, x); }
        }
        push_u64(&mut bytes, current.len());
        bytes.extend_from_slice(&current);
        let checksum = hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Restore the SAME learning experiment without physical calls or rerunning
    /// any completed fit. The existing fixed-proposal codec owns observation
    /// and likelihood decoding. Phase seeds, sampled proposals and paid counts
    /// are reconstructed from the immutable configuration and ordered history.
    ///
    /// Framing and checksums are bounded/checked before allocating history or
    /// observation storage. An uncompleted fit can run on the next advance;
    /// its already paid observations remain available. Positive allowance can
    /// also finalize paid production without another physical call.
    pub fn restore(
        plan: &UqPlan, config: PilotImportanceConfig, model: ContentHash, bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        let maximum = wire_limit(plan, &config)?;
        if bytes.len() < MIN_BYTES || bytes.len() > maximum {
            return Err(invalid("pilot checkpoint length exceeds admitted work"));
        }
        if &bytes[..8] != MAGIC { return Err(invalid("unknown pilot importance checkpoint version")); }
        let end = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..end]).as_bytes().as_slice() != &bytes[end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, &config, model).as_bytes().as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(invalid("unknown or refused pilot status")),
        };
        let mut reader = Reader { bytes: &bytes[41..end], position: 0 };
        let count = reader.count()?;
        if count > config.pilot_rounds {
            return Err(invalid("too many completed importance pilots"));
        }
        let record_bytes = FIT_BYTES + 8 * plan.parameters.len();
        let nested_offset = HEADER + count * record_bytes;
        if nested_offset + 8 + 81 > end {
            return Err(invalid("truncated pilot history or current producer"));
        }
        // Canonical admission still validates the full initial proposal and all
        // work limits; no checkpoint can silently relax the expected policy.
        let mut execution = Self::new(plan, config).map_err(UqCheckpointError::InvalidPlan)?;
        let mut proposal = execution.config.initial_proposal.clone();
        for round in 0..count {
            let fit = read_fit(&mut reader, plan, &execution.config)?;
            execution.completed.push(ImportancePilotRound {
                seed: phase_seed(plan.seed, round), sampled_proposal: proposal,
                fit: fit.clone(),
            });
            proposal = fit.proposal;
        }
        let length = reader.count()?;
        let payload = reader.take(length)?;
        if reader.position != reader.bytes.len() {
            return Err(invalid("trailing pilot checkpoint bytes"));
        }
        let mut current_plan = plan.clone();
        current_plan.seed = phase_seed(plan.seed, count);
        current_plan.budget_max_samples = if count == execution.config.pilot_rounds {
            execution.production_samples
        } else { execution.config.samples_per_pilot };
        let current = ImportanceExecution::restore(&current_plan, proposal, model, payload)?;
        if (status == UqStatus::Complete
                && (count != execution.config.pilot_rounds || current.status != UqStatus::Complete))
            || (status == UqStatus::BudgetTruncated && current.status != UqStatus::BudgetTruncated)
        {
            return Err(invalid("pilot status disagrees with retained phase"));
        }
        execution.current = current;
        execution.status = status;
        Ok(execution)
    }
}

fn read_fit(
    reader: &mut Reader<'_>, plan: &UqPlan, config: &PilotImportanceConfig,
) -> Result<ImportanceProposalFit, UqCheckpointError> {
    let selected = reader.count()?;
    let failures = reader.count()?;
    let cutoff = reader.number()?;
    let targets_failure_event = reader.flag()?;
    let effective_samples = reader.number()?;
    let norm_before_limit = reader.number()?;
    let norm_limited = reader.flag()?;
    let n = config.samples_per_pilot;
    let wanted = ((n as f64 * config.fit.elite_fraction).ceil() as usize).max(2).min(n);
    let threshold = plan.compliance_threshold.ok_or(UqCheckpointError::InvalidPlan("missing failure ceiling"))?;
    if selected < wanted || selected > n || failures > n
        || effective_samples < 1.0 - 1e-12 || effective_samples > selected as f64 * (1.0 + 1e-12)
        || norm_before_limit < 0.0 || norm_limited != (norm_before_limit > config.fit.maximum_shift_norm)
        || targets_failure_event != (failures >= wanted)
        || (targets_failure_event && (selected != failures || cutoff.to_bits() != threshold.to_bits()))
        || (!targets_failure_event && cutoff > threshold)
    {
        return Err(invalid("inconsistent completed pilot fit"));
    }
    let mut shift = Vec::with_capacity(plan.parameters.len());
    for _ in &plan.parameters { shift.push(reader.number()?); }
    // Protect future sampler arithmetic without claiming that saved diagnostics
    // prove a real weighted fit. Original fit observations are intentionally gone.
    let scale = shift.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let norm = if scale == 0.0 { 0.0 } else {
        scale * fs_math::det::sqrt(shift.iter().map(|x| (x / scale).powi(2)).sum())
    };
    let tolerance = 64.0 * f64::EPSILON * plan.parameters.len() as f64;
    if !norm.is_finite() || norm > config.fit.maximum_shift_norm * (1.0 + tolerance) {
        return Err(invalid("saved pilot shift exceeds its numerical norm limit"));
    }
    Ok(ImportanceProposalFit {
        samples: n, selected, observed_failures: failures, cutoff, targets_failure_event,
        effective_samples, norm_before_limit, norm_limited,
        proposal: ImportanceConfig { germ_shift: shift,
            defensive_fraction: config.initial_proposal.defensive_fraction },
    })
}

struct Reader<'a> { bytes: &'a [u8], position: usize }
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], UqCheckpointError> {
        let end = self.position.checked_add(n).ok_or_else(|| invalid("pilot field length overflows"))?;
        let bytes = self.bytes.get(self.position..end).ok_or_else(|| invalid("truncated pilot field"))?;
        self.position = end;
        Ok(bytes)
    }
    fn bits(&mut self) -> Result<u64, UqCheckpointError> {
        let bytes = self.take(8)?.try_into().map_err(|_| invalid("incomplete pilot number"))?;
        Ok(u64::from_le_bytes(bytes))
    }
    fn count(&mut self) -> Result<usize, UqCheckpointError> {
        usize::try_from(self.bits()?).map_err(|_| invalid("pilot count overflows usize"))
    }
    fn number(&mut self) -> Result<f64, UqCheckpointError> {
        let value = f64::from_bits(self.bits()?);
        if value.is_finite() { Ok(value) } else { Err(invalid("nonfinite pilot fit value")) }
    }
    fn flag(&mut self) -> Result<bool, UqCheckpointError> {
        match self.take(1)?[0] {
            0 => Ok(false), 1 => Ok(true), _ => Err(invalid("invalid pilot fit flag")),
        }
    }
}
