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
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftGeometrySource, DraftInterfaceKind, PrimitiveShape,
    ProbeKind, StepPolicy, TimeHistory,
};
use fs_system::{
    BlockSolveMode, BlockSolverPolicy, CardStore, LoweredSystem, LoweringReceipt, MonolithicConfig,
    ParameterSource, SolveConfig,
};

use super::codec::{AdmittedAssembly, AssemblyRequirement, BOX_AXES, ProbeStat, TargetEnvelope};
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
        rec.components
            .iter()
            .filter(|c| c.entity == ent)
            .flat_map(|c| &c.parameter_sources)
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
        let (cl, is_iface) =
            if iface && let Some(c) = cards.interface(h).and_then(|e| e.card.claims().claim(cid)) {
                (c, true)
            } else {
                let c = cards
                    .material(h)
                    .and_then(|c| c.claims().claim(cid))
                    .ok_or_else(|| format!("missing claim `{prop}` in card store"))?;
                (c, false)
            };
        match cl.uncertainty {
            UncertaintyModel::Unstated => {
                unstated.insert(prop.clone());
            }
            UncertaintyModel::HalfWidth { half_width: x, .. }
            | UncertaintyModel::RelativeHalfWidth { fraction: x, .. }
                if x > 0.0 =>
            {
                out.push((*h, is_iface, cid, prop.clone()));
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
                PropertyValue::Surface2D {
                    axes,
                    axis_dims,
                    samples,
                    triangles,
                    dims,
                } => {
                    let mut ns = Vec::with_capacity(samples.len());
                    for &(pt, z) in samples {
                        ns.push((pt, shift(z, &cl.uncertainty, s)?));
                    }
                    PropertyValue::Surface2D {
                        axes: axes.clone(),
                        axis_dims: *axis_dims,
                        samples: ns,
                        triangles: triangles.clone(),
                        dims: *dims,
                    }
                }
            };
        }
        cs.insert_claim(next).map_err(|e| format!("{e}"))?;
    }
    Ok(cs)
}

fn relower_input(
    adm: &AdmittedAssembly,
    inp: DraftAssemblyInput,
    store: &CardStore,
) -> Result<LoweredSystem, String> {
    let scen = Scenario::new(&adm.name, adm.seed, inp.environment.base.clone());
    let draft = DraftAssembly::admit(inp, &scen).map_err(|e| format!("{e}"))?;
    super::lowering_registry()
        .lower(&draft, store)
        .map_err(|e| format!("{e}"))
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
    let (mut remap, mut store, mut inp) = (
        BTreeMap::<(ContentHash, bool), ContentHash>::new(),
        CardStore::new(),
        adm.draft.input().clone(),
    );
    {
        let mut remap_mat = |slot: &mut ContentHash| -> Result<(), String> {
            let old = *slot;
            if let Some(&h) = remap.get(&(old, false)) {
                *slot = h;
                return Ok(());
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
            remap.insert((old, false), h);
            *slot = h;
            Ok(())
        };
        for part in &mut inp.parts {
            for mr in &mut part.materials {
                remap_mat(&mut mr.card)?;
            }
        }
        for iface in &mut inp.interfaces {
            let tb = match &mut iface.kind {
                DraftInterfaceKind::ShearFilm { third_body, .. } => Some(third_body),
                DraftInterfaceKind::ThermalContact { third_body, .. }
                | DraftInterfaceKind::MechanicalContact { third_body, .. }
                | DraftInterfaceKind::ElectricalContact { third_body, .. } => third_body.as_mut(),
                _ => None,
            };
            if let Some(tb) = tb {
                remap_mat(&mut tb.material_card)?;
            }
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
        if let Some(&h) = remap.get(&(old, true)) {
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
        remap.insert((old, true), h);
        *lc = h;
    }
    relower_input(adm, inp, &store)
}

fn lower_geom_corner(adm: &AdmittedAssembly, mask: usize) -> Result<LoweredSystem, String> {
    let mut inp = adm.draft.input().clone();
    for (b, g) in adm.envelopes.geometry.iter().enumerate() {
        let s = if (mask >> b) & 1 == 0 { -1.0 } else { 1.0 };
        let part = inp
            .parts
            .iter_mut()
            .find(|p| p.name == g.part)
            .ok_or("missing geometry part")?;
        let DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) = &mut part.geometry
        else {
            return Err("geometry corner requires box primitive".into());
        };
        [lx, ly, lz][g.axis].value += s * g.half_width_m;
    }
    relower_input(adm, inp, &adm.card_store)
}

fn lower_bc_corner(adm: &AdmittedAssembly, mask: usize) -> Result<LoweredSystem, String> {
    let mut inp = adm.draft.input().clone();
    for (b, bc) in adm.envelopes.boundary_conditions.iter().enumerate() {
        let s = if (mask >> b) & 1 == 0 { -1.0 } else { 1.0 };
        let src = inp
            .sources
            .iter_mut()
            .find(|x| x.name == bc.target)
            .ok_or("missing bc source")?;
        let (TimeHistory::Constant(q) | TimeHistory::Harmonic { offset: q, .. }) = &mut src.history
        else {
            return Err("bc corner requires constant/harmonic source".into());
        };
        q.value += s * bc.half_width.value;
    }
    relower_input(adm, inp, &adm.card_store)
}

fn sweep_corners<B, F>(
    acc: &mut ExtraAcc<'_>,
    kind: &str,
    empty: &str,
    labels: &[String],
    noms: &[(String, Dims, f64)],
    mut build: B,
    run_transient: &mut F,
) -> Result<(Vec<f64>, String), String>
where
    B: FnMut(usize) -> Result<LoweredSystem, String>,
    F: FnMut(&mut LoweredSystem, Option<f64>) -> Result<RunTrace, SolveRefusal>,
{
    if labels.is_empty() {
        return Err(empty.into());
    }
    if labels.len() > 4 {
        return Err(format!("{} inputs exceed 16-corner cap", labels.len()));
    }
    let (corners, mut md) = (1usize << labels.len(), vec![0.0f64; noms.len()]);
    for mask in 0..corners {
        let lbl = format!("{kind}_{mask}");
        let meta = (kind, lbl.as_str(), acc.dt, acc.rel, acc.mono, None);
        let vs = run_extra(acc, meta, || build(mask), run_transient)?;
        for (slot, (&(_, _, nom), v)) in md.iter_mut().zip(noms.iter().zip(vs)) {
            *slot = slot.max((v - nom).abs());
        }
    }
    Ok((
        md,
        format!(
            "max over {corners} corners ({}; omits O(eps^2) cross-category terms, Estimated)",
            labels.join(", ")
        ),
    ))
}

fn direct_term(
    list: &[TargetEnvelope],
    probe: &str,
    dims: Dims,
    empty_msg: &str,
) -> Result<(TermValue, String), SolveRefusal> {
    let Some(e) = list.iter().find(|e| e.target == probe) else {
        return Ok(nodata(empty_msg));
    };
    if e.half_width.dims != dims {
        return Err(req_refusal(format!(
            "envelope on probe `{probe}` dimensions do not match probe"
        )));
    }
    Ok((
        TermValue::interval(0.0, e.half_width.value).unwrap(),
        format!("Declared: {}", e.source),
    ))
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
    dt: f64,
    rel: f64,
    mono: f64,
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
        noms.push((p.unit.clone(), d, nom));
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
        dt,
        rel: d_rel,
        mono: d_mono,
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
                let ns: Vec<_> = uc.iter().map(|c| c.3.clone()).collect();
                sweep_corners(
                    &mut acc,
                    "param_corner",
                    "",
                    &ns,
                    &noms,
                    |m| lower_corner(adm, &uc, m),
                    &mut run_transient,
                )
            }
        },
    };
    let gl: Vec<_> = adm
        .envelopes
        .geometry
        .iter()
        .map(|g| format!("{}.{} ({})", g.part, BOX_AXES[g.axis], g.source))
        .collect();
    let geom = sweep_corners(
        &mut acc,
        "geom_corner",
        "no as-built geometry tolerance declared",
        &gl,
        &noms,
        |m| lower_geom_corner(adm, m),
        &mut run_transient,
    );
    let bl: Vec<_> = adm
        .envelopes
        .boundary_conditions
        .iter()
        .map(|b| format!("{} ({})", b.target, b.source))
        .collect();
    let bc = sweep_corners(
        &mut acc,
        "bc_corner",
        "no boundary-condition envelope declared",
        &bl,
        &noms,
        |m| lower_bc_corner(adm, m),
        &mut run_transient,
    );
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
    for (i, (r, (unit, dims, nom))) in reqs.iter().zip(&noms).enumerate() {
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
        let pick = |res: &Result<(Vec<f64>, String), String>| match res {
            Ok((v, w)) => (TermValue::interval(0.0, v[i]).unwrap(), w.clone()),
            Err(w) => nodata(w),
        };
        let (par, g, b) = (pick(&param), pick(&geom), pick(&bc));
        let m = direct_term(
            &adm.envelopes.model_form,
            &r.probe,
            *dims,
            "no model-form discrepancy bound declared",
        )?;
        let s = direct_term(
            &adm.envelopes.measurement,
            &r.probe,
            *dims,
            "no sensor or measurement uncertainty declared",
        )?;
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
    let mut by_key = BTreeMap::new();
    for (_, key, val) in metrics {
        if key.starts_with("driver.budget.") {
            by_key.insert(key.as_str(), *val);
        }
    }
    let mut out = raw.to_string();
    for (key, val) in by_key {
        let pat = format!("\"wall_s_metric\":{}", js(key));
        out = out.replace(&pat, &format!("{pat},\"wall_s\":{val:.17e}"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyKey, Provenance,
    };
    use fs_qty::QtyAny;

    #[test]
    fn surface2d_perturbation_and_wall_metric_dedup() {
        let pr = Provenance {
            source: "NIST".into(),
            license: "CC0".into(),
            artifact: Some(hash_bytes(b"nist")),
        };
        let mut cs = ClaimSet::new();
        let obs = cs
            .register_observation(ObservationDataset {
                specimen: "srm".into(),
                method: "ghp".into(),
                artifact: hash_bytes(b"srm"),
                caveats: "none".into(),
                provenance: pr.clone(),
            })
            .unwrap();
        let k = Dims([1, 1, -3, -1, 0, 0]);
        let cid = cs
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("thermal_conductivity", k),
                value: PropertyValue::Surface2D {
                    axes: ["T".into(), "rho".into()],
                    axis_dims: [Dims([0, 0, 0, 1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])],
                    samples: vec![
                        ([280.0, 37.0], 0.030),
                        ([320.0, 37.0], 0.036),
                        ([320.0, 46.0], 0.038),
                    ],
                    triangles: vec![[0, 1, 2]],
                    dims: k,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("rho", 37.0, 46.0),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.015,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![obs],
                provenance: pr,
            })
            .unwrap();
        for (s, mult) in [(1.0, 1.015), (-1.0, 0.985)] {
            let p = perturb_claims(&cs, &BTreeMap::from([(cid, s)])).unwrap();
            let PropertyValue::Surface2D {
                samples, triangles, ..
            } = &p.claims_ordered().next().unwrap().1.value
            else {
                panic!()
            };
            assert_eq!(triangles, &vec![[0, 1, 2]]);
            assert!(
                (samples[0].1 - 0.030 * mult).abs() < 1.0e-15
                    && (samples[2].1 - 0.038 * mult).abs() < 1.0e-15
            );
        }
        let glued = attach_wall_metrics(
            r#"{"wall_s_metric":"driver.budget.step_halving.wall_s"}"#,
            &[
                (9, "driver.budget.step_halving.wall_s".into(), 0.1),
                (9, "driver.budget.step_halving.wall_s".into(), 0.2),
            ],
        );
        assert_eq!(glued.matches("\"wall_s\":").count(), 1);
        assert!(glued.contains("2.000000000000000"));
    }

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
