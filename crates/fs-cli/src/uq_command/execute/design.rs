//! Uncertainty-aware selection over a PREDECLARED finite design family.
//! Each candidate retains its own ordinary UqExecution and confidence sequence.
//! The union bound spends alpha once across that fixed family; candidate
//! dependence from common random numbers needs no independence assumption.
use super::*;
use compliance::{Assessment, Decision};
use model::DesignGrid;
use fs_blake3::{ContentHash, DomainHasher};
use std::io::Read;
use std::path::Path;

const MAGIC: &[u8; 8] = b"FSUQDG1\0";

struct Verdict {
    status: &'static str,
    selected: Option<usize>,
    qualified: Option<usize>,
}

fn verdict(grid: &DesignGrid, assessments: &[Assessment]) -> Verdict {
    let mut unresolved = false;
    for index in grid.preference_order() {
        match assessments[index].decision {
            Decision::BelowTarget => {}
            Decision::Indeterminate => unresolved = true,
            Decision::MeetsTarget => return Verdict {
                status: if unresolved { "inconclusive" } else { "selected" },
                selected: (!unresolved).then_some(index),
                qualified: Some(index),
            },
        }
    }
    Verdict { status: if unresolved { "inconclusive" } else { "no-qualified-candidate" },
        selected: None, qualified: None }
}

fn candidate_policy(family: Policy, count: usize) -> Result<Policy> {
    // Downward one ulp ensures the declared allocations do not overspend the
    // family level through rounded division. The underlying CS keeps its
    // existing floating-point/no-outward-rounding authority boundary.
    let alpha = (family.alpha / count as f64).next_down();
    Policy::new(family.required_probability, alpha, family.min_samples)
}

fn binding(grid: &DesignGrid, family: Policy, config: &Config) -> Result<String> {
    let multipliers = grid.multipliers.iter().map(|value| format!("\"{:016x}\"",value.to_bits()))
        .collect::<Vec<_>>().join(",");
    Ok(format!("{{\"design_version\":1,\"allocation\":\"equal-family-alpha-next-down\",\"control\":{},\"candidate_bits\":[{}],\"policy\":{}}}",
        quote(grid.label()),multipliers,family.checkpoint_binding(&config.render_parameters()?)))
}

pub(super) fn execute(base_text: &str, base: &J, config: &Config, options: &Options,
    grid: &DesignGrid) -> Result<ExecutionOutput> {
    grid.validate(base, config.samples)?;
    let family = options.compliance.ok_or_else(|| bad("candidate selection requires the complete compliance policy"))?;
    let policy = candidate_policy(family,grid.multipliers.len())?;
    let plan = config.plan();
    policy.validate_plan(&plan)?;
    let identity = if options.resume.is_some() || options.checkpoint.is_some() {
        Some(checkpoint::model_identity(base_text,&binding(grid,family,config)?)?)
    } else { None };
    let mut executions = match &options.resume {
        Some(path) => restore(path,&plan,identity.expect("resume identity"),grid.multipliers.len())?,
        None => (0..grid.multipliers.len()).map(|_| UqExecution::new(&plan).map_err(bad))
            .collect::<Result<Vec<_>>>()?,
    };
    let mut assessments = executions.iter().map(|execution| policy.assess(execution))
        .collect::<Result<Vec<_>>>()?;
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    if let Some(output) = &output { save(output,&executions,identity.expect("output identity"))?; }
    let allowance = options.max_new_samples.unwrap_or(MAX_PRODUCT_SAMPLES);
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut accepted = 0;
    let mut termination = "sample-budget";
    'candidates: for index in grid.preference_order() {
        if assessments[index].decision == Decision::MeetsTarget { break; }
        while !assessments[index].reached() && executions[index].observations().len() < config.samples {
            if accepted == allowance { termination = "sample-chunk"; break 'candidates; }
            let before = executions[index].observations().len();
            let report = executions[index].advance_interruptible(1,
                || Instant::now() >= deadline, |values| {
                    let request = config.sample_candidate_request(base,values,grid.control,grid.multipliers[index])?;
                    match evaluate_sample_for(&request,deadline,config.qoi) {
                        Ok(value) => Ok(Some(value)),
                        Err(EvaluationError::Budget) => Ok(None),
                        Err(EvaluationError::Child(message)) => Err(model_failure(message)),
                    }
                });
            if report.status == UqStatus::Refused {
                let error = Failure { code: "cooling-network-uq-refused", message: format!(
                    "candidate {index} (multiplier {}): {}",grid.multipliers[index],
                    report.rejection_reason.as_deref().unwrap_or("model refused")) };
                if let Some(output) = &output {
                    output.invalidate(&error.to_string()).map_err(|mut failed| {
                        failed.message.push_str(&format!("; original failure: {error}")); failed
                    })?;
                }
                return Err(error);
            }
            accepted += executions[index].observations().len() - before;
            if let Some(output) = &output { save(output,&executions,identity.expect("output identity"))?; }
            assessments[index] = policy.assess(&executions[index])?;
            if report.status == UqStatus::Cancelled {
                termination = "wall-time-budget";
                break 'candidates;
            }
        }
        // Never infer a candidate verdict from its neighbours or from the
        // empirical success fraction. An unresolved preferred candidate stays
        // unresolved even if a less preferred candidate qualifies.
        if assessments[index].decision == Decision::MeetsTarget { break; }
    }
    let result = verdict(grid,&assessments);
    if result.status != "inconclusive" { termination = "candidate-family-resolved"; }
    let rows = grid.multipliers.iter().enumerate().map(|(index,&multiplier)| {
        let assessment = &assessments[index];
        let interval = match &assessment.estimate {
            Some(estimate) => format!("[{},{}]",number_json(estimate.lo)?,number_json(estimate.hi)?),
            None => "null".into(),
        };
        Ok(format!("{{\"index\":{index},\"multiplier\":{},\"samples_evaluated\":{},\"decision\":{},\"empirical_probability_of_compliance\":{},\"probability_confidence_sequence\":{interval}}}",
            number_json(multiplier)?,executions[index].observations().len(),
            quote(if executions[index].observations().is_empty() { "not-evaluated" } else { assessment.decision.label() }),
            optional_number(assessment.estimate.as_ref().map(|estimate| estimate.mean))?))
    }).collect::<Result<Vec<_>>>()?.join(",");
    // No invocation-specific field on resolved outcomes: terminal resume is
    // byte-identical and launches no child, even with no remaining time.
    let recovery = if result.status == "inconclusive" { format!(",\"checkpoint\":{}",
        options.checkpoint.as_deref().map_or_else(||"null".into(),|path|quote(&path.to_string_lossy()))) }
        else { String::new() };
    let total: usize = executions.iter().map(|execution|execution.observations().len()).sum();
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.design.v1\",\"authority\":\"estimated-model-sampling-confidence\",\"status\":{},\"termination\":{},\"control\":{},\"selected_multiplier\":{},\"qualified_multiplier\":{},\"selection_scope\":\"best predeclared candidate only when every more-preferred candidate is below target; no interpolation or continuous/global optimum\",\"qoi\":{},\"seed\":{},\"temperature_limit_k\":{},\"required_probability\":{},\"family_alpha\":{},\"alpha_per_candidate\":{},\"min_decision_samples\":{},\"samples_per_candidate\":{},\"total_sample_cap\":{},\"samples_evaluated\":{},\"candidates\":[{}],\"correlation\":{},\"parameters\":[{}]{},\"sampling\":\"same predeclared seed and parameter draw ordinal across candidates; independent observations within each fixed model under the sampler assumptions\",\"confidence_scope\":\"union bound over this fixed candidate family and all assessed sample prefixes; no independence required across candidates; no error control across separate invocations with changed families, seeds or policies\",\"no_claim\":\"qualified_multiplier is not a selected optimum when better candidates remain unresolved; no feasibility inference from neighbouring candidates or failed model solves; no physical/model-form uncertainty, outward-rounding, mesh or continuous-time-peak certificate; selected multiplier must be applied to the same sampled-load/speed semantics\"}}\n",
        quote(result.status),quote(termination),quote(grid.label()),
        optional_number(result.selected.map(|i|grid.multipliers[i]))?,
        optional_number(result.qualified.map(|i|grid.multipliers[i]))?,
        config.qoi.render(objective_kind(base)?),quote(&config.seed.to_string()),
        optional_number(config.threshold_k)?,number_json(family.required_probability)?,
        number_json(family.alpha)?,number_json(policy.alpha)?,policy.min_samples,config.samples,
        config.samples*grid.multipliers.len(),total,rows,quote(config.correlation_label),
        config.render_parameters()?,recovery,
    );
    Ok(ExecutionOutput { stdout, exit_code: if result.status == "inconclusive" {exit::BUDGET} else {exit::SUCCESS} })
}

fn candidate_identity(identity: ContentHash, index: usize) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.cooling-uq.design-candidate.v1");
    hash.update(identity.as_bytes());
    hash.update(&(index as u64).to_le_bytes());
    hash.finalize()
}

fn encode(executions: &[UqExecution], identity: ContentHash) -> Result<Vec<u8>> {
    let mut bytes = MAGIC.to_vec();
    bytes.extend_from_slice(&(executions.len() as u64).to_le_bytes());
    for (index,execution) in executions.iter().enumerate() {
        let entry = execution.checkpoint(candidate_identity(identity,index)).map_err(bad)?;
        bytes.extend_from_slice(&(entry.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&entry);
    }
    Ok(bytes)
}

fn save(output: &checkpoint::Output, executions: &[UqExecution], identity: ContentHash) -> Result<()> {
    output.publish(&encode(executions,identity)?)
}

fn restore(path: &Path, plan: &UqPlan, identity: ContentHash, count: usize) -> Result<Vec<UqExecution>> {
    let cap = 16 + count * (1032 + 8 * plan.budget_max_samples);
    let mut bytes = Vec::new();
    File::open(path).map_err(|error|bad(format!("cannot open design checkpoint: {error}")))?
        .take(cap as u64 + 1).read_to_end(&mut bytes)
        .map_err(|error|bad(format!("cannot read design checkpoint: {error}")))?;
    if bytes.len() > cap { return Err(bad("design checkpoint exceeds the admitted sample cap")); }
    decode(&bytes,plan,identity,count)
}

fn decode(mut bytes: &[u8], plan: &UqPlan, identity: ContentHash, count: usize) -> Result<Vec<UqExecution>> {
    fn take<'a>(bytes: &mut &'a [u8], n: usize) -> Result<&'a [u8]> {
        if bytes.len() < n { return Err(bad("truncated design checkpoint")); }
        let (head,tail) = bytes.split_at(n); *bytes = tail; Ok(head)
    }
    fn integer(bytes: &mut &[u8]) -> Result<usize> {
        let mut word = [0;8]; word.copy_from_slice(take(bytes,8)?);
        usize::try_from(u64::from_le_bytes(word)).map_err(|_|bad("checkpoint length overflow"))
    }
    if take(&mut bytes,8)? != MAGIC || integer(&mut bytes)? != count {
        return Err(bad("checkpoint is not for this design candidate family"));
    }
    let mut executions = Vec::with_capacity(count);
    for index in 0..count {
        let n = integer(&mut bytes)?;
        let entry = take(&mut bytes,n)?;
        executions.push(UqExecution::restore(plan,candidate_identity(identity,index),entry).map_err(bad)?);
    }
    if !bytes.is_empty() { return Err(bad("trailing bytes in design checkpoint")); }
    Ok(executions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::DesignControl;

    fn grid() -> DesignGrid { DesignGrid::parse(DesignControl::FanSpeed,"0.5,1,1.5").unwrap() }
    fn plan() -> UqPlan {
        UqPlan::new("temperature",PropagationMethod::MonteCarlo,64)
            .with_parameter(ParameterUncertainty::uniform("ambient",290.0,310.0,"K"))
            .with_correlation(CorrelationModel::Independent).with_compliance_threshold(300.0)
    }
    fn assessment(decision: Decision) -> Assessment { Assessment { estimate: None, decision } }

    #[test]
    fn unresolved_preferred_candidate_cannot_be_skipped_into_an_optimum() {
        let rows = [assessment(Decision::Indeterminate),assessment(Decision::MeetsTarget),assessment(Decision::BelowTarget)];
        let result = verdict(&grid(),&rows);
        assert_eq!(result.status,"inconclusive"); assert_eq!(result.selected,None); assert_eq!(result.qualified,Some(1));
        let rows = [assessment(Decision::BelowTarget),assessment(Decision::MeetsTarget),assessment(Decision::Indeterminate)];
        assert_eq!(verdict(&grid(),&rows).selected,Some(1));
        let rows = [assessment(Decision::BelowTarget),assessment(Decision::BelowTarget),assessment(Decision::BelowTarget)];
        assert_eq!(verdict(&grid(),&rows).status,"no-qualified-candidate");
    }

    #[test]
    fn family_allocation_never_reuses_the_whole_alpha_for_every_candidate() {
        let policy = Policy::new(0.5,0.05,32).unwrap();
        for count in 2..=64 {
            let individual = candidate_policy(policy,count).unwrap();
            assert!(individual.alpha * count as f64 <= policy.alpha);
            assert!(individual.alpha < policy.alpha);
        }
    }

    #[test]
    fn bundle_replays_each_candidate_stream_and_rejects_wrong_identity_or_framing() {
        let plan = plan(); let identity = ContentHash([17;32]);
        let mut states: Vec<_> = (0..3).map(|_|UqExecution::new(&plan).unwrap()).collect();
        states[0].advance(7,||false,|values|Ok::<_,&str>(values[0]));
        states[1].advance(3,||false,|values|Ok::<_,&str>(values[0]+1.0));
        let bytes = encode(&states,identity).unwrap();
        let mut restored = decode(&bytes,&plan,identity,3).unwrap();
        for (index,(original,resumed)) in states.iter_mut().zip(&mut restored).enumerate() {
            original.advance(9,||false,|values|Ok::<_,&str>(values[0]+index as f64));
            resumed.advance(9,||false,|values|Ok::<_,&str>(values[0]+index as f64));
        }
        assert_eq!(encode(&states,identity).unwrap(),encode(&restored,identity).unwrap());
        assert!(decode(&bytes,&plan,ContentHash([18;32]),3).is_err());
        assert!(decode(&bytes,&plan,identity,2).is_err());
        assert!(decode(&bytes[..bytes.len()-1],&plan,identity,3).is_err());
        let mut trailing = bytes.clone(); trailing.push(0);
        assert!(decode(&trailing,&plan,identity,3).is_err());
    }
}
