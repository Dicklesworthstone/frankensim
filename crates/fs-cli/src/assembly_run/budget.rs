//! Eight-term QoI uncertainty budgets and tri-state requirement verdicts for
//! assembly runs (`frankensim-rc-1109-root-3gc12.11.8` Stage 1).

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::uncertainty::{
    ComplianceVerdict, EngineeringUncertaintyBudget, EngineeringUncertaintyKind as Uk,
    EngineeringUncertaintyTerm, FlippingUnknown, RequirementRelation, ScalarRequirement, TermValue,
    UncertaintyArtifactRef,
};
use fs_matdb::{
    ClaimId, ClaimSet, InterfaceSystemCard, MaterialCard, PropertyValue, UncertaintyModel,
};
use fs_qty::Dims;
use fs_report::{BudgetTermItem, RequirementReportItem};
use fs_scenario::Scenario;
use fs_scenario::assembly::draft::{DraftAssembly, DraftInterfaceKind, ProbeKind, StepPolicy};
use fs_system::{
    BlockSolveMode, BlockSolverPolicy, CardStore, LoweredSystem, LoweringReceipt, MonolithicConfig,
    ParameterSource, SolveConfig,
};

use super::codec::{AdmittedAssembly, AssemblyRequirement, ProbeStat};
use super::report::{EvaluatedProbe, RunTrace, evaluate_probes};
use crate::import::json_string as js;
use crate::json_read::JsonValue;
use crate::solve::SolveRefusal;
use crate::stage_driver::{RetainedSideArtifact, staged_refusal};

pub(crate) const QOI_BUDGET_KIND: &str = "assembly-qoi-budget";
pub(crate) const QOI_BUDGET_SCHEMA: &str = "frankensim.cli.assembly-qoi-budget.v1";

#[derive(Debug, Clone)]
pub(crate) struct EvaluatedQoiBudget {
    pub(crate) requirements: Vec<RequirementReportItem>,
    pub(crate) terms: Vec<BudgetTermItem>,
}

fn req_refusal(detail: impl Into<String>) -> SolveRefusal {
    staged_refusal(
        "cli-assembly-requirement",
        "probe",
        detail,
        "ensure every requirement targets a non-empty probe with matching SI dimensions",
    )
}

fn unit_dims(u: &str) -> Option<Dims> {
    match u {
        "K" => Some(Dims([0, 0, 0, 1, 0, 0])),
        "m" => Some(Dims([1, 0, 0, 0, 0, 0])),
        "m/s" => Some(Dims([1, 0, -1, 0, 0, 0])),
        "s" => Some(Dims([0, 0, 1, 0, 0, 0])),
        "Hz" | "rad/s" => Some(Dims([0, 0, -1, 0, 0, 0])),
        "1" => Some(Dims::NONE),
        _ => None,
    }
}

fn probe_stat(vals: &[f64], stat: ProbeStat) -> Option<f64> {
    if vals.is_empty() || vals.iter().any(|v| !v.is_finite()) {
        return None;
    }
    match stat {
        ProbeStat::Max => vals.iter().copied().reduce(f64::max),
        ProbeStat::Min => vals.iter().copied().reduce(f64::min),
        ProbeStat::Final => vals.last().copied(),
    }
}

fn nodata(reason: impl Into<String>) -> (TermValue, String) {
    let r = reason.into();
    (TermValue::unknown(&r).unwrap(), r)
}

fn roundoff_term(
    probe: &str,
    unit: &str,
    adm: &AdmittedAssembly,
    rec: Option<&LoweringReceipt>,
    de: Option<f64>,
) -> (TermValue, String) {
    let Some(rec) = rec else {
        return nodata("no lowering receipt retained");
    };
    let Some(de) = de.filter(|d| d.is_finite() && *d >= 0.0) else {
        return nodata("no transient window-ledger energy defect retained");
    };
    let inp = adm.draft.input();
    let ent = inp
        .probes
        .iter()
        .find(|p| p.name == probe)
        .and_then(|p| match &p.kind {
            ProbeKind::Thermocouple { location, .. } => Some(*location),
            ProbeKind::Timing { target, .. }
            | ProbeKind::FieldQuantity { target, .. }
            | ProbeKind::StrainGauge { target, .. } => Some(*target),
            ProbeKind::Microphone { .. } => None,
        });
    let Some(ent) = ent else {
        return nodata("probe has no single target part");
    };
    let par = |n: &str| {
        let c = rec.components.iter().find(|c| c.entity == ent)?;
        c.parameter_sources
            .iter()
            .find(|p| p.parameter == n && p.value > 0.0)
            .map(|p| p.value)
    };
    let hw = match unit {
        "K" => par("heat_capacity").map(|c| de / c),
        "m/s" => par("mass").map(|m| (2.0 * de / m).sqrt()),
        "m" => inp
            .joints
            .iter()
            .filter(|j| j.child == ent || j.parent == ent)
            .filter_map(|j| j.stiffness.filter(|q| q.value > 0.0).map(|q| q.value))
            .reduce(f64::min)
            .map(|k| (2.0 * de / k).sqrt()),
        _ => None,
    };
    let Some(v) = hw.filter(|v| v.is_finite() && *v >= 0.0) else {
        return nodata("no dimensional energy-to-QoI conversion for this probe");
    };
    let r = format!(
        "ledger defect {de:.3e} J -> {unit}; includes solver residual (conservative overlap with solver-algebraic)"
    );
    (TermValue::interval(0.0, v).unwrap(), r)
}

type UClaim = (ContentHash, bool, ClaimId, String);

fn audit_cards(rec: &LoweringReceipt, cards: &CardStore) -> Result<Vec<UClaim>, String> {
    let (mut seen, mut unstated, mut out, mut total) =
        (BTreeSet::new(), BTreeSet::new(), Vec::new(), 0usize);
    let mut srcs: Vec<(bool, &ParameterSource)> = Vec::new();
    for c in &rec.components {
        srcs.extend(c.parameter_sources.iter().map(|p| (false, &p.source)));
    }
    for i in &rec.interfaces {
        srcs.extend(i.parameter_sources.iter().map(|p| (true, &p.source)));
    }
    for (iface, src) in srcs {
        let ParameterSource::MatDbClaim {
            card_hash: h,
            property: prop,
            selected_claim: s,
            ..
        } = src
        else {
            continue;
        };
        let cid = ClaimId(*s);
        if !seen.insert((*h, cid)) {
            continue;
        }
        total += 1;
        let cl = if iface {
            cards.interface(h).and_then(|e| e.card.claims().claim(cid))
        } else {
            cards.material(h).and_then(|c| c.claims().claim(cid))
        }
        .ok_or_else(|| format!("missing claim `{prop}` in card store"))?;
        match cl.uncertainty {
            UncertaintyModel::Unstated => {
                unstated.insert(prop.clone());
            }
            UncertaintyModel::HalfWidth { half_width: x, .. }
            | UncertaintyModel::RelativeHalfWidth { fraction: x, .. }
                if x > 0.0 =>
            {
                out.push((*h, iface, cid, prop.clone()));
            }
            _ => {}
        }
    }
    if total == 0 {
        return Err("no material or interface card claims consumed".into());
    }
    if !unstated.is_empty() {
        let ks: Vec<_> = unstated.into_iter().collect();
        return Err(format!("unstated card uncertainty on: {}", ks.join(", ")));
    }
    if out.len() > 4 {
        return Err(format!(
            "{} uncertain claims exceed 16-corner cap",
            out.len()
        ));
    }
    Ok(out)
}

fn perturb_claims(src: &ClaimSet, edits: &BTreeMap<ClaimId, f64>) -> Result<ClaimSet, String> {
    let shift = |v: f64, u: &UncertaintyModel, s: f64| -> Result<f64, String> {
        let hw = match *u {
            UncertaintyModel::HalfWidth { half_width, .. } => half_width,
            UncertaintyModel::RelativeHalfWidth { fraction, .. } => v.abs() * fraction,
            UncertaintyModel::Unstated => return Err("unstated uncertainty".into()),
        };
        let y = v + s * hw;
        y.is_finite()
            .then_some(y)
            .ok_or_else(|| "non-finite".into())
    };
    let mut cs = ClaimSet::new();
    for oid in src.observation_ids() {
        if let Some(obs) = src.observation(oid) {
            cs.register_observation(obs.clone())
                .map_err(|e| format!("{e}"))?;
        }
    }
    for (cid, cl) in src.claims_ordered() {
        let mut next = cl.clone();
        if let Some(&s) = edits.get(&cid) {
            next.value = match &cl.value {
                PropertyValue::Scalar { value, dims } => PropertyValue::Scalar {
                    value: shift(*value, &cl.uncertainty, s)?,
                    dims: *dims,
                },
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    knots,
                    dims,
                } => {
                    let mut nk = Vec::with_capacity(knots.len());
                    for &(x, y) in knots {
                        nk.push((x, shift(y, &cl.uncertainty, s)?));
                    }
                    PropertyValue::Curve {
                        abscissa: abscissa.clone(),
                        abscissa_dims: *abscissa_dims,
                        knots: nk,
                        dims: *dims,
                    }
                }
                PropertyValue::Surface2D { .. } => return Err("Surface2D unsupported".into()),
            };
        }
        cs.insert_claim(next).map_err(|e| format!("{e}"))?;
    }
    Ok(cs)
}

fn lower_corner(
    adm: &AdmittedAssembly,
    uc: &[UClaim],
    mask: usize,
) -> Result<LoweredSystem, String> {
    let mut edits: BTreeMap<(ContentHash, bool), BTreeMap<ClaimId, f64>> = BTreeMap::new();
    for (b, (h, iface, cid, _)) in uc.iter().enumerate() {
        edits
            .entry((*h, *iface))
            .or_default()
            .insert(*cid, if (mask >> b) & 1 == 0 { -1.0 } else { 1.0 });
    }
    let (mut remap, mut store, mut inp) =
        (BTreeMap::new(), CardStore::new(), adm.draft.input().clone());
    for part in &mut inp.parts {
        for mr in &mut part.materials {
            let old = mr.card;
            if let Some(&h) = remap.get(&old) {
                mr.card = h;
                continue;
            }
            let c = adm
                .card_store
                .material(&old)
                .ok_or("missing material card")?;
            let nc = match edits.get(&(old, false)) {
                Some(e) => MaterialCard::assemble(
                    c.id().clone(),
                    perturb_claims(c.claims(), e)?,
                    c.models().to_vec(),
                )
                .map_err(|e| format!("{e}"))?,
                None => c.clone(),
            };
            let h = store.insert_material(nc);
            remap.insert(old, h);
            mr.card = h;
        }
    }
    for iface in &mut inp.interfaces {
        let lc = match &mut iface.kind {
            DraftInterfaceKind::ThermalContact { law_card, .. }
            | DraftInterfaceKind::MechanicalContact { law_card, .. }
            | DraftInterfaceKind::Glue { law_card }
            | DraftInterfaceKind::ElectricalContact { law_card, .. } => law_card,
            DraftInterfaceKind::PowerPort { .. }
            | DraftInterfaceKind::HeatRoute { .. }
            | DraftInterfaceKind::ShearFilm { .. } => continue,
        };
        let old = *lc;
        if let Some(&h) = remap.get(&old) {
            *lc = h;
            continue;
        }
        let en = adm
            .card_store
            .interface(&old)
            .ok_or("missing interface card")?;
        let nc = match edits.get(&(old, true)) {
            Some(e) => InterfaceSystemCard::assemble(
                en.card.surface_a().clone(),
                en.card.surface_b().clone(),
                en.card.context().clone(),
                perturb_claims(en.card.claims(), e)?,
                en.card.models().to_vec(),
            )
            .map_err(|e| format!("{e}"))?,
            None => en.card.clone(),
        };
        let h = store
            .insert_interface(en.family.clone(), nc)
            .map_err(|e| format!("{e}"))?;
        remap.insert(old, h);
        *lc = h;
    }
    let scen = Scenario::new(&adm.name, adm.seed, inp.environment.base.clone());
    let draft = DraftAssembly::admit(inp, &scen).map_err(|e| format!("{e}"))?;
    super::lowering_registry()
        .lower(&draft, &store)
        .map_err(|e| format!("{e}"))
}

fn ev_id(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        let lo = c.to_ascii_lowercase();
        if lo == '/' {
            out.push_str("-per-");
        } else if lo.is_ascii_lowercase()
            || lo.is_ascii_digit()
            || matches!(lo, '.' | '_' | ':' | '-')
        {
            out.push(lo);
        } else {
            out.push('-');
        }
    }
    if out.is_empty() { "id".into() } else { out }
}

fn assess_one(
    req: &AssemblyRequirement,
    unit: &str,
    nom: f64,
    specs: [(Uk, TermValue, String, &'static str); 8],
    src_h: ContentHash,
) -> Result<(RequirementReportItem, Vec<BudgetTermItem>, String), SolveRefusal> {
    let qoi = format!("{}.{}", req.probe, req.stat.as_str());
    let (eq, eu, ei) = (ev_id(&qoi), ev_id(unit), ev_id(&req.id));
    let (mut ev, mut terms) = (Vec::with_capacity(8), Vec::with_capacity(8));
    for (kind, tv, reason, owner) in specs {
        let (st, val) = match &tv {
            TermValue::IntervalBound { upper, .. } => ("measured", Some(*upper)),
            TermValue::Negligible { .. } => ("negligible", Some(0.0)),
            _ => ("no-data", None),
        };
        let p =
            UncertaintyArtifactRef::new(owner, src_h).map_err(|e| req_refusal(format!("{e}")))?;
        ev.push(
            EngineeringUncertaintyTerm::try_new(kind, tv, p)
                .map_err(|e| req_refusal(format!("{e}")))?,
        );
        terms.push(BudgetTermItem {
            qoi: qoi.clone(),
            kind: kind.name().into(),
            state: st.into(),
            value: val,
            reason,
            owner: owner.into(),
        });
    }
    let b = EngineeringUncertaintyBudget::try_new(&eq, &eu, ev)
        .map_err(|e| req_refusal(format!("{e}")))?;
    let rp = UncertaintyArtifactRef::new("assembly-requirement", src_h)
        .map_err(|e| req_refusal(format!("{e}")))?;
    let sr = ScalarRequirement::try_new(&ei, &eq, &eu, req.relation, req.limit.value, rp)
        .map_err(|e| req_refusal(format!("{e}")))?;
    let v = b
        .assess_requirement(nom, &sr, &[])
        .map_err(|e| req_refusal(format!("{e}")))?;
    let margin = if req.relation == RequirementRelation::AtLeast {
        nom - req.limit.value
    } else {
        req.limit.value - nom
    };
    let (outcome, unm, hr, bh) = match v {
        ComplianceVerdict::Compliant { budget, .. } => ("pass", Vec::new(), None, budget),
        ComplianceVerdict::NonCompliant { budget, .. } => ("fail", Vec::new(), None, budget),
        ComplianceVerdict::Indeterminate {
            flipping_unknowns: fu,
            budget,
            ..
        } => (
            "indeterminate",
            terms
                .iter()
                .filter(|t| t.state == "no-data")
                .map(|t| t.kind.clone())
                .collect(),
            fu.iter()
                .map(FlippingUnknown::required_magnitude)
                .reduce(f64::min),
            budget,
        ),
    };
    let item = RequirementReportItem {
        id: req.id.clone(),
        qoi,
        region: req.probe.clone(),
        effective_limit: req.limit.value,
        required_margin: 0.0,
        nominal_margin: margin,
        unit: unit.into(),
        outcome: outcome.into(),
        unmeasured_terms: unm,
        decision_headroom: hr,
    };
    Ok((item, terms, bh.to_hex()))
}

struct ExtraAcc<'a> {
    adm: &'a AdmittedAssembly,
    t0: Instant,
    consumed: f64,
    json: Vec<String>,
    metrics: Vec<(i64, String, f64)>,
}

fn run_extra<B, F>(
    acc: &mut ExtraAcc<'_>,
    meta: (&str, &str, f64, f64, f64, Option<f64>),
    build: B,
    run_transient: &mut F,
) -> Result<Vec<f64>, String>
where
    B: FnOnce() -> Result<LoweredSystem, String>,
    F: FnMut(&mut LoweredSystem, Option<f64>) -> Result<RunTrace, SolveRefusal>,
{
    let (kind, label, dt_s, rel, mono, dt_ov) = meta;
    let (cap, mkey) = (
        acc.adm.wall_budget_s,
        format!("driver.budget.{label}.wall_s"),
    );
    let (res, steps, status, detail) =
        if cap - (acc.consumed + acc.t0.elapsed().as_secs_f64()) <= 0.0 {
            let d = format!("wall budget {cap:.3}s exhausted before {label}");
            (Err(d.clone()), 0, "budget_exceeded", d)
        } else {
            let st0 = Instant::now();
            let r = build().and_then(|mut ls| {
                let tr = run_transient(&mut ls, dt_ov).map_err(|e| e.what)?;
                let (_, _, ep) =
                    evaluate_probes(&acc.adm.draft.input().probes, &acc.adm.entity_names, &tr)
                        .map_err(|e| e.what)?;
                let mut vals = Vec::with_capacity(acc.adm.requirements.len());
                for rq in &acc.adm.requirements {
                    let v = ep
                        .iter()
                        .find(|p| p.name == rq.probe)
                        .and_then(|p| probe_stat(&p.values, rq.stat));
                    vals.push(v.ok_or_else(|| "extra run missing probe".to_string())?);
                }
                Ok((tr.times_s.len().saturating_sub(1), vals))
            });
            acc.metrics
                .push((9, mkey.clone(), st0.elapsed().as_secs_f64()));
            match r {
                Ok((n, vals)) => (Ok(vals), n, "ok", format!("{label} completed")),
                Err(e) => {
                    let d = format!("{label} refused: {e}");
                    (Err(d.clone()), 0, "refused", d)
                }
            }
        };
    let (jk, jl, js_s, jd, jm) = (js(kind), js(label), js(status), js(&detail), js(&mkey));
    acc.json.push(format!("{{\"kind\":{jk},\"label\":{jl},\"dt_s\":{dt_s:.17e},\"rel_tol\":{rel:.17e},\"mono_tol\":{mono:.17e},\"steps\":{steps},\"status\":{js_s},\"detail\":{jd},\"wall_s_metric\":{jm}}}"));
    res
}

fn tight_system(adm: &AdmittedAssembly, t_rel: f64, t_mono: f64) -> Result<LoweredSystem, String> {
    let mut ls = super::lowering_registry()
        .lower(&adm.draft, &adm.card_store)
        .map_err(|e| format!("{e}"))?;
    let def_sc = SolveConfig::default();
    ls.system.set_config(SolveConfig {
        relative_tolerance: t_rel,
        interface_energy_tolerance: def_sc.interface_energy_tolerance * 1.0e-2,
        ..def_sc
    });
    let ms = ls.system.block_modes();
    let ok = |m: &BlockSolveMode| {
        matches!(
            m,
            BlockSolveMode::OrderedSweep
                | BlockSolveMode::LoadLineNewton
                | BlockSolveMode::MonolithicNewtonKrylov
        )
    };
    if ms.contains(&BlockSolveMode::MonolithicNewtonKrylov) && ms.iter().all(ok) {
        let mc = MonolithicConfig {
            absolute_tolerance: t_mono,
            ..MonolithicConfig::default()
        };
        ls.system
            .set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(mc));
    }
    Ok(ls)
}

pub(crate) type BudgetEvalOut = (
    EvaluatedQoiBudget,
    RetainedSideArtifact,
    Vec<(i64, String, f64)>,
);

#[allow(clippy::too_many_lines)]
pub(crate) fn evaluate_qoi_budgets<F>(
    adm: &AdmittedAssembly,
    probes: &[EvaluatedProbe],
    rec: Option<&LoweringReceipt>,
    energy_defect_j: Option<f64>,
    consumed_wall_s: f64,
    mut run_transient: F,
) -> Result<BudgetEvalOut, SolveRefusal>
where
    F: FnMut(&mut LoweredSystem, Option<f64>) -> Result<RunTrace, SolveRefusal>,
{
    let reqs = &adm.requirements;
    let mut noms = Vec::with_capacity(reqs.len());
    for r in reqs {
        let p = probes
            .iter()
            .find(|p| p.name == r.probe)
            .ok_or_else(|| req_refusal(format!("probe `{}` missing", r.probe)))?;
        let d = unit_dims(&p.unit).ok_or_else(|| {
            req_refusal(format!("probe `{}` unit `{}` unsupported", p.name, p.unit))
        })?;
        if r.limit.dims != d {
            return Err(req_refusal(format!(
                "requirement `{}` limit dimensions do not match `{}` ({})",
                r.id, p.name, p.unit
            )));
        }
        let nom = probe_stat(&p.values, r.stat)
            .ok_or_else(|| req_refusal(format!("probe `{}` empty/non-finite", p.name)))?;
        noms.push((p.unit.clone(), nom));
    }
    let dt = match adm.draft.input().time_plan.step_policy {
        StepPolicy::Fixed { dt } => dt.value,
        StepPolicy::Adaptive { dt_initial, .. } => dt_initial.value,
    };
    let (d_rel, d_mono) = (
        SolveConfig::default().relative_tolerance,
        MonolithicConfig::default().absolute_tolerance,
    );
    let (t_rel, t_mono) = (d_rel * 1.0e-2, d_mono * 1.0e-1);
    let mut acc = ExtraAcc {
        adm,
        t0: Instant::now(),
        consumed: consumed_wall_s,
        json: Vec::new(),
        metrics: Vec::new(),
    };
    let relower = || {
        super::lowering_registry()
            .lower(&adm.draft, &adm.card_store)
            .map_err(|e| format!("{e}"))
    };
    let half = run_extra(
        &mut acc,
        (
            "step_halving",
            "step_halving",
            dt * 0.5,
            d_rel,
            d_mono,
            Some(dt * 0.5),
        ),
        relower,
        &mut run_transient,
    );
    let tight = run_extra(
        &mut acc,
        ("solver_tight", "solver_tight", dt, t_rel, t_mono, None),
        || tight_system(adm, t_rel, t_mono),
        &mut run_transient,
    );
    let param = match rec {
        None => Err("no lowering receipt for card uncertainty audit".into()),
        Some(lr) => match audit_cards(lr, &adm.card_store) {
            Err(w) => Err(w),
            Ok(uc) if uc.is_empty() => Ok((
                vec![0.0; reqs.len()],
                "all claims declare zero half-width".into(),
            )),
            Ok(uc) => {
                let (corners, mut md, mut fail) =
                    (1usize << uc.len(), vec![0.0f64; reqs.len()], None);
                for mask in 0..corners {
                    let lbl = format!("param_corner_{mask}");
                    match run_extra(
                        &mut acc,
                        ("param_corner", &lbl, dt, d_rel, d_mono, None),
                        || lower_corner(adm, &uc, mask),
                        &mut run_transient,
                    ) {
                        Ok(vs) => {
                            for (s, (&(_, n), v)) in md.iter_mut().zip(noms.iter().zip(vs)) {
                                *s = s.max((v - n).abs());
                            }
                        }
                        Err(e) => {
                            fail = Some(e);
                            break;
                        }
                    }
                }
                let ns: Vec<_> = uc.iter().map(|c| c.3.as_str()).collect();
                fail.map_or_else(
                    || {
                        Ok((
                            md,
                            format!("max over {corners} corners ({}, Estimated)", ns.join(", ")),
                        ))
                    },
                    Err,
                )
            }
        },
    };
    let (mut req_items, mut term_items) = (
        Vec::with_capacity(reqs.len()),
        Vec::with_capacity(reqs.len() * 8),
    );
    let mut json = format!(
        "{{\"schema\":{},\"wall_budget_s\":{:.17e},\"requirements\":[",
        js(QOI_BUDGET_SCHEMA),
        adm.wall_budget_s
    );
    let opt = |o: Option<f64>| o.map_or_else(|| "null".into(), |v| format!("{v:.17e}"));
    for (i, (r, (unit, nom))) in reqs.iter().zip(&noms).enumerate() {
        let ro = roundoff_term(&r.probe, unit, adm, rec, energy_defect_j);
        let sol = match &tight {
            Ok(v) => (
                TermValue::interval(0.0, (nom - v[i]).abs()).unwrap(),
                format!("tightened solver (rel->{t_rel:.0e}, mono->{t_mono:.0e}, Estimated)"),
            ),
            Err(w) => nodata(w),
        };
        let disc = match &half {
            Ok(v) => (
                TermValue::interval(0.0, 2.0 * (nom - v[i]).abs()).unwrap(),
                format!(
                    "2*|stat(dt)-stat(dt/2)| for nominal dt={dt}s assuming order p>=1 (Estimated)"
                ),
            ),
            Err(w) => nodata(w),
        };
        let par = match &param {
            Ok((v, w)) => (TermValue::interval(0.0, v[i]).unwrap(), w.clone()),
            Err(w) => nodata(w),
        };
        let (g, b, m, s) = (
            nodata("no as-built geometry tolerance declared"),
            nodata("no boundary-condition envelope declared"),
            nodata("no model-form discrepancy bound declared"),
            nodata("no sensor or measurement uncertainty declared"),
        );
        let specs = [
            (Uk::Roundoff, ro.0, ro.1, "window-ledger"),
            (Uk::SolverAlgebraic, sol.0, sol.1, "solver-tight"),
            (Uk::Discretization, disc.0, disc.1, "step-halving"),
            (Uk::Geometry, g.0, g.1, "geometry"),
            (Uk::Parameters, par.0, par.1, "card-store"),
            (Uk::BoundaryConditions, b.0, b.1, "boundary-conditions"),
            (Uk::ModelForm, m.0, m.1, "model-form"),
            (Uk::Measurement, s.0, s.1, "measurement"),
        ];
        let (item, terms, bh) = assess_one(r, unit, *nom, specs, adm.source_hash)?;
        if i > 0 {
            json.push(',');
        }
        let rel = if r.relation == RequirementRelation::AtLeast {
            "at-least"
        } else {
            "at-most"
        };
        let unm: Vec<_> = item.unmeasured_terms.iter().map(|u| js(u)).collect();
        let tj: Vec<_> = terms.iter().map(|t| {
            let (jk, jst, jv, jr, jo) = (js(&t.kind), js(&t.state), opt(t.value), js(&t.reason), js(&t.owner));
            format!("{{\"kind\":{jk},\"state\":{jst},\"value\":{jv},\"reason\":{jr},\"owner\":{jo}}}")
        }).collect();
        let (jid, jp, jst, jq, ju, jrel, jout, jbh, jhr) = (
            js(&item.id),
            js(&r.probe),
            js(r.stat.as_str()),
            js(&item.qoi),
            js(&item.unit),
            js(rel),
            js(&item.outcome),
            js(&bh),
            opt(item.decision_headroom),
        );
        let (lim, nm, us, ts) = (
            item.effective_limit,
            item.nominal_margin,
            unm.join(","),
            tj.join(","),
        );
        let _ = write!(
            json,
            "{{\"id\":{jid},\"probe\":{jp},\"stat\":{jst},\"qoi\":{jq},\"unit\":{ju},\"relation\":{jrel},\"limit\":{lim:.17e},\"nominal\":{nom:.17e},\"nominal_margin\":{nm:.17e},\"outcome\":{jout},\"budget_hash\":{jbh},\"decision_headroom\":{jhr},\"unmeasured_terms\":[{us}],\"terms\":[{ts}]}}"
        );
        req_items.push(item);
        term_items.extend(terms);
    }
    let _ = write!(json, "],\"extra_runs\":[{}]}}", acc.json.join(","));
    let bytes = json.into_bytes();
    let side = RetainedSideArtifact {
        kind: QOI_BUDGET_KIND,
        artifact: hash_bytes(&bytes),
        bytes,
    };
    Ok((
        EvaluatedQoiBudget {
            requirements: req_items,
            terms: term_items,
        },
        side,
        acc.metrics,
    ))
}

impl EvaluatedQoiBudget {
    pub(crate) fn from_json(text: &str) -> Result<Self, SolveRefusal> {
        let root =
            JsonValue::parse(text).map_err(|e| req_refusal(format!("bad budget json: {e}")))?;
        let (mut requirements, mut terms) = (Vec::new(), Vec::new());
        for r in root
            .get("requirements")
            .and_then(JsonValue::as_array)
            .into_iter()
            .flatten()
        {
            let qoi = r.str_field("qoi").unwrap_or("").to_string();
            for t in r
                .get("terms")
                .and_then(JsonValue::as_array)
                .into_iter()
                .flatten()
            {
                terms.push(BudgetTermItem {
                    qoi: qoi.clone(),
                    kind: t.str_field("kind").unwrap_or("").into(),
                    state: t.str_field("state").unwrap_or("no-data").into(),
                    value: t.f64_field("value"),
                    reason: t.str_field("reason").unwrap_or("").into(),
                    owner: t.str_field("owner").unwrap_or("").into(),
                });
            }
            let unm = r
                .get("unmeasured_terms")
                .and_then(JsonValue::as_array)
                .map(|a| {
                    a.iter()
                        .filter_map(|v| v.as_str().map(str::to_string))
                        .collect()
                })
                .unwrap_or_default();
            requirements.push(RequirementReportItem {
                id: r.str_field("id").unwrap_or("").into(),
                qoi,
                region: r.str_field("probe").unwrap_or("").into(),
                effective_limit: r.f64_field("limit").unwrap_or(0.0),
                required_margin: 0.0,
                nominal_margin: r.f64_field("nominal_margin").unwrap_or(0.0),
                unit: r.str_field("unit").unwrap_or("SI").into(),
                outcome: r.str_field("outcome").unwrap_or("indeterminate").into(),
                unmeasured_terms: unm,
                decision_headroom: r.f64_field("decision_headroom"),
            });
        }
        Ok(Self {
            requirements,
            terms,
        })
    }
}

pub(crate) fn attach_wall_metrics(raw: &str, metrics: &[(i64, String, f64)]) -> String {
    let mut out = raw.to_string();
    for (_, key, val) in metrics {
        if key.starts_with("driver.budget.") {
            let pat = format!("\"wall_s_metric\":{}", js(key));
            out = out.replace(&pat, &format!("{pat},\"wall_s\":{val:.17e}"));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::QtyAny;

    #[test]
    fn assess_one_emits_pass_fail_and_indeterminate_via_evidence_budget() {
        let mk = |lim| AssemblyRequirement {
            id: "t_max".into(),
            probe: "hot_temp".into(),
            stat: ProbeStat::Max,
            relation: RequirementRelation::AtMost,
            limit: QtyAny::new(lim, Dims([0, 0, 0, 1, 0, 0])),
        };
        let bnd = || Uk::ALL.map(|k| (k, TermValue::interval(0.0, 0.5).unwrap(), "b".into(), "t"));
        let h = hash_bytes(b"unit-test");
        assert_eq!(
            assess_one(&mk(310.0), "K", 300.0, bnd(), h)
                .unwrap()
                .0
                .outcome,
            "pass"
        );
        assert_eq!(
            assess_one(&mk(310.0), "K", 320.0, bnd(), h)
                .unwrap()
                .0
                .outcome,
            "fail"
        );
        assert_eq!(
            assess_one(&mk(310.0), "K", 308.0, bnd(), h)
                .unwrap()
                .0
                .outcome,
            "indeterminate"
        );
        let mut nd = bnd();
        nd[4].1 = TermValue::unknown("unstated").unwrap();
        let ind = assess_one(&mk(310.0), "K", 300.0, nd, h).unwrap().0;
        assert_eq!(
            (ind.outcome.as_str(), ind.unmeasured_terms.as_slice()),
            ("indeterminate", &["parameters".to_string()][..])
        );
        assert!(ind.decision_headroom.is_some_and(|x| x > 0.0));
    }
}
