//! fs-evidence conformance suite (CONTRACT.md: any reimplementation must
//! pass). G0: conservativeness laws of the composition algebra (composed
//! enclosures contain true propagation; composed validity ⊆ intersection).
//! Plus the acceptance battery: the registration lint, the worked
//! model-discrepancy-dominates example, the out-of-distribution refusal,
//! and deterministic ledger rows. Aggregate verdicts use the canonical fs-obs
//! schema; seeded cases carry their actual input seed.

use fs_evidence::{
    Ambition, Color, ColorRank, DecisionStatus, DiscrepancyModel, EscalationAdvice, Evidence,
    FidelityPair, IntervalOp, ModelBracket, ModelCard, ModelEvidence, ModelRegistry,
    NumericalCertificate, Op, ProvenanceHash, StatisticalCertificate, UncertaintySource,
    ValidityDomain, compose, validate_color_payload,
};
use fs_propcheck::Stream;
use std::collections::{BTreeMap, BTreeSet};

const SUITE: &str = "fs-evidence/conformance";
const FIXED_INPUT_SEED: u64 = 0;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("evidence verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("evidence verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

/// In-house LCG (high bits — the low-bit period trap is documented in
/// fs-exec's suite).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }
}

fn pt(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
    pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
}

#[test]
fn evd_001_g0_composition_conservativeness_battery() {
    const SEED: u64 = 0x0E01_2026_0706_C0DE;
    let mut rng = Lcg(SEED);
    let p = ProvenanceHash::of_bytes(b"g0");
    let mut checked = 0u32;
    for _ in 0..20_000 {
        // Random enclosures around random true values.
        let ta = (rng.unit() - 0.5) * 200.0;
        let tb = (rng.unit() - 0.5) * 200.0;
        let (wa, wb) = (rng.unit() * 5.0, rng.unit() * 5.0);
        let a = Evidence::enclosed(ta, ta - wa, ta + wa, p);
        let b = Evidence::enclosed(tb, tb - wb, tb + wb, p);
        for op in [Op::Add, Op::Sub, Op::Mul, Op::Min, Op::Max] {
            let c = Evidence::combine(op, &a, &b, ());
            // Sample true values inside the operand enclosures; the
            // composed enclosure must contain every propagation.
            for _ in 0..3 {
                let va = ta - wa + 2.0 * wa * rng.unit();
                let vb = tb - wb + 2.0 * wb * rng.unit();
                let true_out = match op {
                    Op::Add => va + vb,
                    Op::Sub => va - vb,
                    Op::Mul => va * vb,
                    Op::Min => va.min(vb),
                    Op::Max => va.max(vb),
                };
                assert!(
                    c.numerical.lo <= true_out && true_out <= c.numerical.hi,
                    "containment broke: {op:?} {true_out} outside \
                     [{}, {}] (seed {SEED:#x})",
                    c.numerical.lo,
                    c.numerical.hi
                );
                checked += 1;
            }
        }
    }
    // Validity conservativeness: composed ⊆ intersection, on a fixed case.
    let ma = ModelEvidence {
        cards: vec!["a".into()],
        assumptions: vec!["assume-a".into()],
        validity: ValidityDomain::unconstrained().with("Re", 1e4, 1e5),
        discrepancy_rel: 0.02,
        in_domain: true,
    };
    let mb = ModelEvidence {
        cards: vec!["b".into()],
        assumptions: vec!["assume-b".into()],
        validity: ValidityDomain::unconstrained().with("Re", 5e4, 2e5),
        discrepancy_rel: 0.03,
        in_domain: true,
    };
    let pa = ProvenanceHash::of_bytes(b"a");
    let ea = Evidence::exact(1.0, pa).with_model(ma);
    let eb = Evidence::exact(2.0, pa).with_model(mb);
    let ec = Evidence::combine(Op::Add, &ea, &eb, ());
    let inside = pt(&[("Re", 7e4)]);
    let outside = pt(&[("Re", 2e4)]);
    let validity_law = ec.model.validity.contains(&inside)
        && !ec.model.validity.contains(&outside)
        && (ec.model.discrepancy_rel - 0.05).abs() < 1e-12
        && ec.model.assumptions == vec!["assume-a".to_string(), "assume-b".to_string()];
    verdict(
        "evd-001",
        checked == 300_000 && validity_law,
        &format!(
            "composed enclosures contain true propagation over {checked} seeded samples \
             (seed {SEED:#x}); validity intersects, assumptions union, bands add"
        ),
        SEED,
    );
}

#[test]
fn evd_002_model_card_registration_lint() {
    let mut reg = ModelRegistry::new();
    let refused = reg
        .register_solver("flux.lbm-les", "les-smagorinsky")
        .is_err();
    reg.register_card(ModelCard::new(
        "les-smagorinsky",
        "1.0.0",
        Ambition::Frontier,
        vec!["resolved-eddy regime".into()],
        ValidityDomain::unconstrained().with("Re", 1e4, 1e6),
        vec!["under-predicts separation".into()],
        0.10,
    ));
    let accepted = reg
        .register_solver("flux.lbm-les", "les-smagorinsky")
        .is_ok();
    // Rows ride the observability schema toward the ledger.
    let mut em = fs_obs::Emitter::new("fs-evidence/conformance", "evd-002/cards");
    let mut rows_valid = true;
    for row in reg.to_ledger_rows_json() {
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "fs-evidence-model-card".to_string(),
                    json: row,
                },
                None,
            )
            .to_jsonl();
        rows_valid &= fs_obs::validate_line(&line).is_ok();
        println!("{line}");
    }
    verdict(
        "evd-002",
        refused && accepted && rows_valid,
        "a solver without a card cannot register (teaching refusal); with the card it can; \
         model_cards rows validate through fs-obs",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_003_worked_example_model_discrepancy_dominates() {
    // THE core-argument scenario: mesh error 0.7%, statistical 0.5%,
    // LES closure discrepancy 10% — the drag number is NOT decision-grade
    // at a 5% threshold, and the report must say model-form dominates.
    let p = ProvenanceHash::of_bytes(b"drag-v3");
    let card = ModelCard::new(
        "les-smagorinsky",
        "1.0.0",
        Ambition::Frontier,
        vec!["resolved-eddy regime".into()],
        ValidityDomain::unconstrained().with("Re", 1e4, 1e6),
        vec![],
        0.10,
    );
    let at = pt(&[("Re", 2e5)]);
    let drag = Evidence::enclosed(0.0312, 0.0312 * (1.0 - 0.007), 0.0312 * (1.0 + 0.007), p)
        .with_statistical(StatisticalCertificate::HalfWidth {
            half_width: 0.0312 * 0.005,
            confidence: 0.95,
        })
        .with_model(ModelEvidence::from_card(&card, &at));
    let status = drag.assess(0.05);
    let (dominated, detail) = match &status {
        DecisionStatus::NotDecisionGrade { dominant, detail } => {
            (*dominant == UncertaintySource::ModelForm, detail.clone())
        }
        DecisionStatus::DecisionGrade => (false, "unexpectedly decision-grade".to_string()),
    };
    let advice_ok = drag.escalation_advice(0.05) == EscalationAdvice::EscalateModelFidelity;
    // The same numbers with a calibrated 2% closure ARE decision-grade —
    // proof the verdict tracks the model band, not the numerics.
    let better_card = ModelCard::new(
        "les-calibrated",
        "1.0.0",
        Ambition::Frontier,
        vec![],
        ValidityDomain::unconstrained().with("Re", 1e4, 1e6),
        vec![],
        0.02,
    );
    let better = drag
        .clone()
        .with_model(ModelEvidence::from_card(&better_card, &at));
    let flips = matches!(better.assess(0.05), DecisionStatus::DecisionGrade);
    let mut emitter = fs_obs::Emitter::new(SUITE, "evd-003/drag");
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: "fs-evidence-row".to_string(),
            json: drag.to_ledger_row_json(),
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("drag evidence row must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("drag evidence row must use the fs-obs wire schema");
    println!("{line}");
    verdict(
        "evd-003",
        dominated && advice_ok && flips && detail.contains("model-form"),
        &format!("beautifully-certified-wrong caught: {detail}"),
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_004_discrepancy_model_flags_out_of_distribution_queries() {
    const SEED: u64 = 0x0E04_D15C_2026_0706;
    let mut rng = Lcg(SEED);
    // Synthetic two-fidelity ledger corpus: panel-method vs LES drag over
    // Re ∈ [1e4, 1e5], lo-fi biased ~6-9%.
    let pairs: Vec<FidelityPair> = (0..50)
        .map(|_| {
            let re = 1e4 + rng.unit() * 9e4;
            let hi = 1.0 + 0.1 * rng.unit();
            FidelityPair {
                params: pt(&[("Re", re)]),
                lo_fi: hi * (1.06 + 0.03 * rng.unit()),
                hi_fi: hi,
            }
        })
        .collect();
    let model = DiscrepancyModel::fit(&pairs).expect("fit");
    let band = model.query(&pt(&[("Re", 5e4)])).expect("in-domain");
    let band_sane = band.mean_rel > 0.05 && band.max_rel < 0.12 && band.max_rel >= band.mean_rel;
    let err = model
        .query(&pt(&[("Re", 1e6)]))
        .expect_err("out-of-distribution must refuse");
    let teaching = err.to_string().contains("extrapolation") && err.param == "Re";
    // The refusal also blocks certification through evidence_at.
    let ood_blocked = model
        .evidence_at("panel-vs-les", &pt(&[("Re", 1e6)]))
        .is_err();
    let unexpected_dimension_blocked = model
        .evidence_at("panel-vs-les", &pt(&[("Re", 5e4), ("Mach", 0.08)]))
        .is_err();
    let simulated_model = model
        .evidence_at("panel-vs-les", &pt(&[("Re", 5e4)]))
        .expect("in-domain discrepancy evidence");
    let simulated_color =
        fs_evidence::color_of(&NumericalCertificate::enclosure(0.9, 1.1), &simulated_model);
    let simulation_stays_estimated = matches!(
        simulated_color,
        fs_evidence::Color::Estimated { dispersion, .. }
            if dispersion.to_bits() == band.max_rel.to_bits()
    );
    verdict(
        "evd-004",
        band_sane
            && teaching
            && ood_blocked
            && unexpected_dimension_blocked
            && simulation_stays_estimated,
        &format!(
            "trained on 50 pairs (seed {SEED:#x}): in-domain band mean {:.3}/max {:.3}; \
             out-of-domain and unexpected-dimension queries refused; paired simulations remain \
             Estimated rather than impersonating an experimental anchor",
            band.mean_rel, band.max_rel
        ),
        SEED,
    );
}

#[test]
fn evd_005_bracketing_reports_spread_and_rows_are_deterministic() {
    // The vessel flagship's stated mitigation: bracket contact-angle
    // models, report the objective's sensitivity band.
    let bracket = ModelBracket::new()
        .try_with_member("contact-angle-60", 0.90)
        .expect("valid first bracket member")
        .try_with_member("contact-angle-90", 1.00)
        .expect("valid second bracket member")
        .try_with_member("contact-angle-120", 1.16)
        .expect("valid third bracket member");
    let ev = bracket
        .evidence(ProvenanceHash::of_bytes(b"vessel-lip-v1"))
        .expect("bracket evidence");
    let enclosing = ev.numerical.lo <= 0.90 && ev.numerical.hi >= 1.16;
    let spread_reported = ev.model.discrepancy_rel > 0.2
        && ev
            .sensitivity
            .d_qoi
            .contains_key("model-choice(bracket-spread)");
    // A bracket is honest evidence but NOT certifiable past a threshold
    // that its spread violates.
    let status = ev.assess(0.05);
    let not_decision_grade = matches!(
        &status,
        DecisionStatus::NotDecisionGrade { dominant, .. }
            if *dominant == UncertaintySource::ModelForm
    );
    // Ledger rows: canonical and repeatable (G5-class determinism).
    let row1 = ev.to_ledger_row_json();
    let row2 = ev.to_ledger_row_json();
    let mut em = fs_obs::Emitter::new("fs-evidence/conformance", "evd-005/bracket");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "fs-evidence-row".to_string(),
                json: row1.clone(),
            },
            None,
        )
        .to_jsonl();
    let valid = fs_obs::validate_line(&line).is_ok();
    println!("{line}");
    verdict(
        "evd-005",
        enclosing && spread_reported && not_decision_grade && row1 == row2 && valid,
        "bracketed contact-angle models: enclosure spans members, spread carried as the \
         model band, correctly not decision-grade at 5%, rows deterministic + schema-valid",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_006_certified_discipline_composes() {
    let p = ProvenanceHash::of_bytes(b"cert");
    // Rigorous chain: exact times enclosure stays certifiable.
    let a = Evidence::exact(2.0, p);
    let b = Evidence::enclosed(3.0, 2.9, 3.1, p);
    let c = Evidence::combine(Op::Mul, &a, &b, ());
    let chain_ok = c.certified().is_ok();
    // Injecting an estimate anywhere poisons certification downstream.
    let mut est = Evidence::exact(1.0, p);
    est.numerical = NumericalCertificate::estimate(0.9, 1.1);
    let d = Evidence::combine(Op::Add, &est, &b, ());
    let poisoned = d.certified().is_err();
    verdict(
        "evd-006",
        chain_ok && poisoned,
        "Certified discipline survives rigorous composition and refuses estimate-tainted \
         chains (kind severity is monotone)",
        FIXED_INPUT_SEED,
    );
}

fn scalar_value_mismatch_refuses(p: ProvenanceHash) -> bool {
    let mut evidence = Evidence::exact(1.0, p);
    evidence.value = 2.0;
    matches!(
        evidence.certified(),
        Err(fs_evidence::CertifyError::ScalarValueMismatch {
            value: 2.0,
            qoi: 1.0
        })
    )
}

fn impossible_model_validity_refuses(p: ProvenanceHash) -> bool {
    [
        ValidityDomain::unconstrained()
            .with("Re", 1.0, 2.0)
            .intersect(&ValidityDomain::unconstrained().with("Re", 3.0, 4.0)),
        ValidityDomain::unconstrained().with("Re", f64::NAN, 2.0),
        ValidityDomain::unconstrained().with("Re", 1.0, f64::INFINITY),
    ]
    .into_iter()
    .all(|validity| {
        let evidence = Evidence::exact(1.0, p).with_model(ModelEvidence {
            cards: vec!["forged-regime".to_string()],
            assumptions: Vec::new(),
            validity,
            discrepancy_rel: 0.01,
            in_domain: true,
        });
        matches!(
            evidence.clone().certified(),
            Err(fs_evidence::CertifyError::InvalidModelValidity)
        ) && evidence.breakdown().model_rel.is_infinite()
            && matches!(
                evidence.assess(0.05),
                DecisionStatus::NotDecisionGrade { .. }
            )
    })
}

fn certification_round_trip_revalidates(p: ProvenanceHash) -> (bool, bool, bool) {
    let a = Evidence::exact(2.0, p);
    let b = Evidence::enclosed(3.0, 2.9, 3.1, p);
    let certified = Evidence::combine(Op::Mul, &a, &b, 6.0)
        .certified()
        .expect("rigorous chain certifies");
    let readable =
        certified.qoi.to_bits() == 6.0f64.to_bits() && certified.evidence().numerical.lo <= 6.0;
    let mut reopened = certified.into_evidence();
    reopened.numerical = NumericalCertificate::estimate(5.9, 6.1);
    let weakened_refused = matches!(
        reopened.clone().certified(),
        Err(fs_evidence::CertifyError::NotRigorous { .. })
    );
    reopened.numerical = NumericalCertificate::enclosure(5.9, 6.1);
    reopened.value = 7.0;
    reopened.qoi = 7.0;
    let drifted_refused = matches!(
        reopened.certified(),
        Err(fs_evidence::CertifyError::QoiOutsideEnclosure { .. })
    );
    (readable, weakened_refused, drifted_refused)
}

/// gp3.2.1 — the Certified<T> trust boundary: every public forge route
/// is refused with the structured reason; the opaque newtype has no
/// mutable access, so weakening means an EXPLICIT downgrade whose
/// re-certification re-validates (the reconstruction round trip).
#[test]
fn evd_012_certified_is_unforgeable_and_validated() {
    use fs_evidence::{CertifyError, NumericalKind};
    let p = ProvenanceHash::of_bytes(b"adversarial");
    let forge = |numerical: NumericalCertificate, qoi: f64| {
        let mut e = Evidence::exact(0.0, p);
        e.qoi = qoi;
        e.value = qoi;
        e.numerical = numerical;
        e.certified()
    };
    // (a) NaN/infinite "exact" values refuse (Evidence::exact accepts
    // them as EVIDENCE; the certification boundary is where they stop).
    let nan_exact = matches!(
        Evidence::exact(f64::NAN, p).certified(),
        Err(CertifyError::ExactNotFinite { .. })
    );
    let inf_exact = matches!(
        Evidence::exact(f64::INFINITY, p).certified(),
        Err(CertifyError::ExactNotFinite { .. })
    );
    // (b) A hand-built Exact whose bounds contradict its QoI refuses.
    let inconsistent = matches!(
        forge(
            NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: 1.0,
                hi: 2.0,
            },
            1.0,
        ),
        Err(CertifyError::ExactInconsistent { .. })
    );
    // (c) Hand-built inverted / NaN / infinite enclosures refuse: the
    // validator checks the NUMBERS, not the constructor that claimed them.
    let inverted = matches!(
        forge(
            NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: 2.0,
                hi: 1.0,
            },
            1.5,
        ),
        Err(CertifyError::InvertedBounds { .. })
    );
    let nan_bound = matches!(
        forge(
            NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: f64::NAN,
                hi: 1.0,
            },
            0.5,
        ),
        Err(CertifyError::NonFiniteBounds { .. })
    );
    let inf_bound = matches!(
        forge(
            NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: 0.0,
                hi: f64::INFINITY,
            },
            0.5,
        ),
        Err(CertifyError::NonFiniteBounds { .. })
    );
    // (d) A QoI outside its own claimed enclosure refuses.
    let escaped = matches!(
        forge(NumericalCertificate::enclosure(0.0, 1.0), 2.0),
        Err(CertifyError::QoiOutsideEnclosure { .. })
    );
    // (e) Scalar evidence cannot carry one value while certifying another.
    let scalar_mismatch = scalar_value_mismatch_refuses(p);
    // (f) A public ModelEvidence literal cannot override an impossible
    // validity box by merely asserting `in_domain: true`.
    let invalid_model_validity = impossible_model_validity_refuses(p);
    // (g-h) Legitimate composition certifies, while an explicit downgrade
    // loses the mark and every reconstruction re-enters validation.
    let (readable, weakened_refused, drifted_refused) = certification_round_trip_revalidates(p);
    verdict(
        "evd-012",
        nan_exact
            && inf_exact
            && inconsistent
            && inverted
            && nan_bound
            && inf_bound
            && escaped
            && scalar_mismatch
            && invalid_model_validity
            && readable
            && weakened_refused
            && drifted_refused,
        "Certified<T> is opaque: numerical, scalar-value, and model-validity forge routes are \
         refused with structured reasons; rigorous composition still certifies; \
         downgrade-mutate-recertify re-validates",
        FIXED_INPUT_SEED,
    );
}

fn indeterminate_arithmetic_fails_closed(p: ProvenanceHash) -> bool {
    // Public Evidence composition must preserve an honest whole-line
    // enclosure when IEEE arithmetic is indeterminate. In particular,
    // min/max folds must not silently discard the NaN corners of 0 * ±∞.
    let zero = Evidence::exact(0.0, p);
    let whole = Evidence::enclosed(1.0, f64::NEG_INFINITY, f64::INFINITY, p);
    let positive_infinity = Evidence::enclosed(f64::INFINITY, f64::INFINITY, f64::INFINITY, p);
    let negative_infinity =
        Evidence::enclosed(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, p);
    [
        Evidence::combine(Op::Mul, &zero, &whole, ()),
        Evidence::combine(Op::Add, &positive_infinity, &negative_infinity, ()),
        Evidence::combine(Op::Sub, &positive_infinity, &positive_infinity, ()),
    ]
    .iter()
    .all(|evidence| {
        evidence.numerical.lo == f64::NEG_INFINITY
            && evidence.numerical.hi == f64::INFINITY
            && evidence.breakdown().numerical_rel.is_infinite()
    })
}

fn malformed_statistics_fail_closed(p: ProvenanceHash) -> bool {
    use fs_evidence::CertifyError;

    let malformed_statistics = [
        StatisticalCertificate::EValue {
            e: -1.0,
            alpha: 0.05,
        },
        StatisticalCertificate::EValue {
            e: 1.0,
            alpha: f64::NAN,
        },
        StatisticalCertificate::HalfWidth {
            half_width: -0.1,
            confidence: 0.95,
        },
        StatisticalCertificate::HalfWidth {
            half_width: 0.1,
            confidence: 1.0,
        },
    ];
    let direct_inputs_fail_closed = malformed_statistics.into_iter().all(|statistical| {
        let evidence = Evidence::exact(1.0, p).with_statistical(statistical);
        matches!(
            evidence.clone().certified(),
            Err(CertifyError::InvalidStatistical { .. })
        ) && evidence.breakdown().statistical_rel.is_infinite()
            && matches!(
                evidence.assess(0.05),
                DecisionStatus::NotDecisionGrade { .. }
            )
    });
    let invalid_statistical_operand =
        Evidence::exact(1.0, p).with_statistical(StatisticalCertificate::EValue {
            e: f64::NAN,
            alpha: 0.05,
        });
    let valid_statistical_operand =
        Evidence::exact(1.0, p).with_statistical(StatisticalCertificate::HalfWidth {
            half_width: 0.1,
            confidence: 0.95,
        });
    let composed_statistical = Evidence::combine(
        Op::Add,
        &invalid_statistical_operand,
        &valid_statistical_operand,
        (),
    );
    let statistical_composition_fails_closed = matches!(
        composed_statistical.clone().certified(),
        Err(CertifyError::InvalidStatistical { .. })
    ) && composed_statistical
        .breakdown()
        .statistical_rel
        .is_infinite();
    direct_inputs_fail_closed && statistical_composition_fails_closed
}

fn malformed_model_discrepancies_fail_closed(p: ProvenanceHash) -> bool {
    use fs_evidence::CertifyError;

    let direct_inputs_fail_closed =
        [f64::NAN, -0.1, f64::NEG_INFINITY]
            .into_iter()
            .all(|discrepancy_rel| {
                let mut model = ModelEvidence::none();
                model.discrepancy_rel = discrepancy_rel;
                let evidence = Evidence::exact(1.0, p).with_model(model);
                matches!(
                    evidence.clone().certified(),
                    Err(CertifyError::InvalidModelDiscrepancy { .. })
                ) && evidence.breakdown().model_rel.is_infinite()
                    && matches!(
                        evidence.assess(0.05),
                        DecisionStatus::NotDecisionGrade { .. }
                    )
            });

    // Invalid negative uncertainty must not cancel a valid positive band
    // during composition and emerge as apparently exact model evidence.
    let mut invalid_model = ModelEvidence::none();
    invalid_model.discrepancy_rel = -0.2;
    let invalid_operand = Evidence::exact(1.0, p).with_model(invalid_model);
    let mut valid_model = ModelEvidence::none();
    valid_model.discrepancy_rel = 0.2;
    let valid_operand = Evidence::exact(1.0, p).with_model(valid_model);
    let composed = Evidence::combine(Op::Add, &invalid_operand, &valid_operand, ());
    direct_inputs_fail_closed
        && composed.model.discrepancy_rel == f64::INFINITY
        && matches!(
            composed.assess(0.05),
            DecisionStatus::NotDecisionGrade { .. }
        )
}

#[test]
fn evd_013_malformed_uncertainty_and_indeterminate_arithmetic_fail_closed() {
    let p = ProvenanceHash::of_bytes(b"fail-closed-evidence");

    // Positive infinity is a valid, explicit unbounded model claim. It may
    // retain the Certified integrity mark, but can never be decision-grade.
    let mut unbounded_model = ModelEvidence::none();
    unbounded_model.discrepancy_rel = f64::INFINITY;
    let unbounded = Evidence::exact(1.0, p).with_model(unbounded_model);
    let unbounded_is_honest = unbounded.clone().certified().is_ok()
        && unbounded.breakdown().model_rel == f64::INFINITY
        && matches!(
            unbounded.assess(0.05),
            DecisionStatus::NotDecisionGrade { .. }
        );

    let invalid_threshold_refuses = [f64::NAN, f64::INFINITY, -0.1]
        .into_iter()
        .all(|threshold| {
            matches!(
                Evidence::exact(1.0, p).assess(threshold),
                DecisionStatus::NotDecisionGrade { .. }
            )
        });

    verdict(
        "evd-013",
        indeterminate_arithmetic_fails_closed(p)
            && malformed_statistics_fail_closed(p)
            && malformed_model_discrepancies_fail_closed(p)
            && unbounded_is_honest
            && invalid_threshold_refuses,
        "indeterminate interval arithmetic widens to the whole line; malformed statistical \
         and model uncertainty cannot certify or become decision-grade; an explicit infinite \
         model band remains honest but non-decision-grade",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_007_disjoint_validated_regimes_demote_not_launder() {
    use fs_evidence::{Color, IntervalOp, compose};
    let a = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 1e5, 2e5),
        dataset: "A".into(),
    };
    // Overlapping regimes: stays Validated with the (narrower) intersection.
    let b_overlap = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 1.5e5, 3e5),
        dataset: "B".into(),
    };
    let overlap = compose(&a, &b_overlap, IntervalOp::Add);
    let overlap_ok = if let Color::Validated { regime, dataset } = &overlap {
        let (lo, hi) = regime.bound("Re").unwrap_or((0.0, 0.0));
        (lo - 1.5e5).abs() < 1.0
            && (hi - 2e5).abs() < 1.0
            && dataset == "derived:v2:datasets:1:A+1:B"
    } else {
        false
    };
    // DISJOINT regimes: no state satisfies both anchors, so the composition must
    // NOT stay Validated (previously it laundered a phantom [3e5,3e5] regime).
    let b_disjoint = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 3e5, 4e5),
        dataset: "B".into(),
    };
    let disjoint = compose(&a, &b_disjoint, IntervalOp::Add);
    let demoted = matches!(
        disjoint,
        Color::Estimated { estimator, dispersion }
            if estimator == "derived:v2:disjoint-regimes:1:A+1:B"
                && dispersion.is_infinite()
    );
    verdict(
        "evd-007",
        overlap_ok && demoted,
        "validated⊕validated intersects when overlapping; disjoint regimes demote to \
         Estimated with infinite/no-claim dispersion instead of laundering a phantom point \
         regime or asserting a zero spread",
        FIXED_INPUT_SEED,
    );
}

#[test]
#[allow(clippy::float_cmp)] // exact next_down/next_up bit checks
fn evd_008_verified_compose_stays_a_true_enclosure() {
    use fs_evidence::{Color, IntervalOp, compose};
    // Round-to-nearest a+b can EXCLUDE the exact real sum — 0.1+0.2 is the
    // textbook case (the exact sum lies below the nearest f64). A Verified
    // composition must OUTWARD-round so the interval still brackets the exact
    // result: lo < a+b < hi, strictly.
    let a = Color::Verified { lo: 0.1, hi: 0.1 };
    let b = Color::Verified { lo: 0.2, hi: 0.2 };
    let Color::Verified { lo, hi } = compose(&a, &b, IntervalOp::Add) else {
        panic!("Verified ⊕ Verified must stay Verified");
    };
    let mid = 0.1_f64 + 0.2_f64;
    assert!(
        lo < mid && mid < hi,
        "Add not outward-rounded: [{lo}, {hi}] vs {mid}"
    );
    assert_eq!(lo, mid.next_down());
    assert_eq!(hi, mid.next_up());

    // A 0×∞ product corner is indeterminate and the min/max fold would SILENTLY
    // DROP the NaN, reporting a bogus tight interval — it must degrade to the
    // whole real line instead.
    let zero = Color::Verified { lo: 0.0, hi: 0.0 };
    let whole = Color::Verified {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };
    let Color::Verified { lo, hi } = compose(&zero, &whole, IntervalOp::Mul) else {
        panic!("stays Verified");
    };
    assert!(
        lo == f64::NEG_INFINITY && hi == f64::INFINITY,
        "0×∞ Mul must be the whole real line, got [{lo}, {hi}]"
    );

    verdict(
        "evd-008",
        true,
        "Verified compose is outward-rounded (true enclosure) and NaN-safe",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_009_non_finite_regimes_fail_closed_and_payloads_escape() {
    use fs_evidence::{Color, check_regime, regime_demotion};

    let validated = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 1.0, 10.0),
        dataset: "anchors".to_string(),
    };
    let mut non_finite_states_demote = true;
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let state = pt(&[("Re", value)]);
        let borrowed = regime_demotion(&validated, &state);
        let (checked, flag) = check_regime(&validated, &state);
        non_finite_states_demote &= matches!(
            checked,
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ) && flag
            .as_ref()
            .is_some_and(|d| d.axis == "Re" && !d.value.is_finite())
            && matches!(
                (&borrowed, &flag),
                (Some(left), Some(right))
                    if left.dataset == right.dataset
                        && left.axis == right.axis
                        && left.value.to_bits() == right.value.to_bits()
            );
    }

    let inverted = ValidityDomain::unconstrained()
        .with("Re", 1.0, 2.0)
        .intersect(&ValidityDomain::unconstrained().with("Re", 3.0, 4.0));
    let invalid_regimes = [
        ValidityDomain::unconstrained(),
        ValidityDomain::unconstrained().with("Re", f64::NEG_INFINITY, 10.0),
        ValidityDomain::unconstrained().with("Re", f64::NAN, f64::NAN),
        ValidityDomain::unconstrained().with("Re", f64::NAN, 10.0),
        ValidityDomain::unconstrained().with("Re", 1.0, f64::NAN),
        inverted,
    ];
    let mut invalid_regimes_demote = true;
    for regime in invalid_regimes {
        let invalid = Color::Validated {
            regime,
            dataset: "anchors".to_string(),
        };
        let state = pt(&[("Re", 5.0)]);
        let borrowed = regime_demotion(&invalid, &state);
        let (checked, flag) = check_regime(&invalid, &state);
        invalid_regimes_demote &= matches!(
            checked,
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ) && flag.is_some()
            && matches!(
                (&borrowed, &flag),
                (Some(left), Some(right))
                    if left.dataset == right.dataset
                        && left.axis == right.axis
                        && left.value.to_bits() == right.value.to_bits()
            );
    }

    let hostile = Color::Estimated {
        estimator: "est\"\\\n\r\t\u{0007}".to_string(),
        dispersion: f64::INFINITY,
    };
    let payload = hostile.payload_json();
    let payload_is_escaped = !payload.chars().any(char::is_control)
        && payload.contains(r#"est\"\\\n\r\t\u0007"#)
        && payload.contains(r#""dispersion":"non-finite:inf""#)
        && payload == hostile.payload_json();

    verdict(
        "evd-009",
        non_finite_states_demote && invalid_regimes_demote && payload_is_escaped,
        "NaN and infinite state values, empty/non-finite/inverted regimes all demote \
         with infinite dispersion; hostile payload metadata is escaped deterministically \
         and non-finite floats are tagged strings",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_010_verified_gate_refuses_nan_and_inverted_bounds() {
    use fs_evidence::{Color, color_of, verified_from};
    // exact(NaN) is not an enclosure — the color gate must refuse, not mint
    // Verified{NaN,NaN} (bead wa8i E4).
    assert!(verified_from(&NumericalCertificate::exact(f64::NAN)).is_err());
    // enclosure(NaN, hi) now PROPAGATES the NaN (min/max used to drop it into a
    // razor-thin false interval), so the gate refuses it too.
    assert!(verified_from(&NumericalCertificate::enclosure(f64::NAN, 1.0)).is_err());
    // A valid finite enclosure passes; an out-of-order INPUT is normalized.
    assert!(matches!(
        verified_from(&NumericalCertificate::enclosure(1.0, 5.0)),
        Ok(Color::Verified { lo, hi }) if lo <= 1.0 && hi >= 5.0
    ));
    assert!(matches!(
        verified_from(&NumericalCertificate::enclosure(5.0, 1.0)), // normalized to [1,5]
        Ok(Color::Verified { .. })
    ));
    // Infinite bounds are a valid (loose) enclosure and pass.
    assert!(matches!(
        verified_from(&NumericalCertificate::enclosure(
            f64::NEG_INFINITY,
            f64::INFINITY
        )),
        Ok(Color::Verified { .. })
    ));
    // color_of routes uncarded numerics through the same gate: valid → Verified,
    // NaN → falls through to Estimated (never a false Verified).
    let uncarded = ModelEvidence {
        cards: vec![],
        assumptions: vec![],
        validity: ValidityDomain::unconstrained(),
        discrepancy_rel: 0.0,
        in_domain: true,
    };
    assert_eq!(
        color_of(&NumericalCertificate::enclosure(1.0, 2.0), &uncarded).rank(),
        fs_evidence::ColorRank::Verified
    );
    assert_ne!(
        color_of(&NumericalCertificate::exact(f64::NAN), &uncarded).rank(),
        fs_evidence::ColorRank::Verified
    );
    verdict(
        "evd-010",
        true,
        "verified_from / color_of refuse NaN and inverted bounds (fail closed)",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_011_color_canonical_identity_is_versioned_and_bit_exact() {
    use core::fmt::Write as _;
    use fs_evidence::Color;

    let signed_zero = Color::Verified { lo: 0.0, hi: -0.0 };
    let encoded = signed_zero.canonical_bytes();
    let mut encoded_hex = String::with_capacity(encoded.len() * 2);
    for byte in encoded {
        write!(&mut encoded_hex, "{byte:02x}").expect("writing to String cannot fail");
    }
    assert_eq!(
        encoded_hex, "02000800000000000000000000000000000008000000000000000000000000000080",
        "this vector freezes Color canonical encoding v2"
    );

    let first = Color::Verified { lo: 1.0, hi: 2.0 };
    let next = Color::Verified {
        lo: 1.0f64.next_up(),
        hi: 2.0,
    };
    assert_eq!(first.payload_json(), next.payload_json());
    assert_ne!(first.canonical_bytes(), next.canonical_bytes());

    let forward = Color::Validated {
        regime: ValidityDomain::unconstrained()
            .with("alpha", -0.0, 1.0)
            .with("beta", 2.0, 3.0),
        dataset: "anchors".to_string(),
    };
    let reverse = Color::Validated {
        regime: ValidityDomain::unconstrained()
            .with("beta", 2.0, 3.0)
            .with("alpha", -0.0, 1.0),
        dataset: "anchors".to_string(),
    };
    assert_eq!(forward.canonical_bytes(), reverse.canonical_bytes());

    verdict(
        "evd-011",
        true,
        "Color canonical encoding v2 is frozen, bit-exact for f64 payloads, and \
         deterministic across validity-domain insertion order",
        FIXED_INPUT_SEED,
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one bounded-identity adversarial composition matrix
fn evd_014_color_provenance_composition_is_bounded_without_laundering() {
    use fs_evidence::{
        Color, ColorRank, IntervalOp, MAX_COLOR_IDENTITY_BYTES, color_leaf_identity_reason,
        color_of, compose,
    };

    let source = Color::Estimated {
        estimator: "wedge-proposer-v1".to_string(),
        dispersion: 0.1,
    };
    let repeated = compose(&source, &source, IntervalOp::Hull);
    assert!(matches!(
        &repeated,
        Color::Estimated {
            estimator,
            dispersion,
        } if estimator
            == "derived:v2:composed:17:wedge-proposer-v1+17:wedge-proposer-v1"
            && dispersion.to_bits() == 0.2f64.to_bits()
    ));
    let verified = Color::Verified { lo: 1.0, hi: 1.0 };
    let validated_pass_through = compose(
        &Color::Validated {
            regime: ValidityDomain::unconstrained().with("x", 0.0, 1.0),
            dataset: "dataset-0".to_string(),
        },
        &verified,
        IntervalOp::Mul,
    );
    let Color::Validated { dataset, .. } = validated_pass_through else {
        panic!("a verified transform cannot strengthen or weaken a validated parent")
    };
    assert_eq!(
        dataset,
        "derived:v2:datasets-verified:9:dataset-0+8:verified"
    );
    assert_eq!(
        color_leaf_identity_reason(&dataset),
        Some("derived-identity-requires-lineage")
    );

    let same_estimator_twice = compose(
        &Color::Estimated {
            estimator: "shared".to_string(),
            dispersion: 0.0,
        },
        &Color::Estimated {
            estimator: "shared".to_string(),
            dispersion: 0.0,
        },
        IntervalOp::Hull,
    );
    let estimator_with_verified = compose(
        &Color::Estimated {
            estimator: "shared".to_string(),
            dispersion: 0.0,
        },
        &verified,
        IntervalOp::Hull,
    );
    assert_ne!(
        same_estimator_twice.canonical_bytes(),
        estimator_with_verified.canonical_bytes(),
        "two Estimated operands cannot alias an Estimated-plus-Verified derivation"
    );

    let shared_regime = ValidityDomain::unconstrained().with("x", 0.0, 1.0);
    let shared_validated = Color::Validated {
        regime: shared_regime,
        dataset: "shared".to_string(),
    };
    let validated_twice = compose(&shared_validated, &shared_validated, IntervalOp::Hull);
    let validated_with_verified = compose(&shared_validated, &verified, IntervalOp::Hull);
    assert_ne!(
        validated_twice.canonical_bytes(),
        validated_with_verified.canonical_bytes(),
        "two Validated operands cannot alias a Validated-plus-Verified derivation"
    );

    let estimator_named_verified = compose(
        &Color::Estimated {
            estimator: "shared".to_string(),
            dispersion: 0.0,
        },
        &Color::Estimated {
            estimator: "verified".to_string(),
            dispersion: 0.0,
        },
        IntervalOp::Hull,
    );
    assert_ne!(
        estimator_named_verified.canonical_bytes(),
        estimator_with_verified.canonical_bytes(),
        "an estimator literally named verified cannot impersonate the Verified operand class"
    );

    let dataset_named_verified = compose(
        &shared_validated,
        &Color::Validated {
            regime: ValidityDomain::unconstrained().with("x", 0.0, 1.0),
            dataset: "verified".to_string(),
        },
        IntervalOp::Hull,
    );
    assert_ne!(
        dataset_named_verified.canonical_bytes(),
        validated_with_verified.canonical_bytes(),
        "a dataset literally named verified cannot impersonate the Verified operand class"
    );

    let mut folded = source;
    for index in 0..1_024 {
        folded = compose(
            &folded,
            &Color::Estimated {
                estimator: format!("independent-estimator-{index}"),
                dispersion: 0.0,
            },
            IntervalOp::Hull,
        );
        let Color::Estimated { estimator, .. } = &folded else {
            panic!("estimated composition must remain estimated")
        };
        assert!(estimator.len() <= MAX_COLOR_IDENTITY_BYTES);
    }
    assert_eq!(folded.rank(), ColorRank::Estimated);
    assert_eq!(
        folded,
        (0..1_024).fold(
            Color::Estimated {
                estimator: "wedge-proposer-v1".to_string(),
                dispersion: 0.1,
            },
            |color, index| compose(
                &color,
                &Color::Estimated {
                    estimator: format!("independent-estimator-{index}"),
                    dispersion: 0.0,
                },
                IntervalOp::Hull,
            )
        ),
        "bounded identities replay deterministically"
    );

    let overlong = "x".repeat(MAX_COLOR_IDENTITY_BYTES + 1);
    let overlong_estimate = Color::Estimated {
        estimator: overlong.clone(),
        dispersion: 0.0,
    };
    let compact_repeated = compose(&overlong_estimate, &overlong_estimate, IntervalOp::Hull);
    let compact_with_verified = compose(
        &overlong_estimate,
        &Color::Verified { lo: 0.0, hi: 1.0 },
        IntervalOp::Hull,
    );
    assert_ne!(
        compact_repeated.canonical_bytes(),
        compact_with_verified.canonical_bytes(),
        "domain labels must also separate operand classes after compact hashing"
    );
    for bounded in [compact_repeated, compact_with_verified] {
        let Color::Estimated { estimator, .. } = bounded else {
            panic!("an estimated operand must cap the result at Estimated")
        };
        assert!(estimator.len() <= MAX_COLOR_IDENTITY_BYTES);
    }

    let left_grouping = compose(
        &Color::Estimated {
            estimator: "A+B".to_string(),
            dispersion: 0.0,
        },
        &Color::Estimated {
            estimator: "C".to_string(),
            dispersion: 0.0,
        },
        IntervalOp::Hull,
    );
    let right_grouping = compose(
        &Color::Estimated {
            estimator: "A".to_string(),
            dispersion: 0.0,
        },
        &Color::Estimated {
            estimator: "B+C".to_string(),
            dispersion: 0.0,
        },
        IntervalOp::Hull,
    );
    assert_ne!(
        left_grouping, right_grouping,
        "length framing must distinguish identities containing the separator"
    );

    let regime = ValidityDomain::unconstrained().with("x", 0.0, 1.0);
    let mut validated = Color::Validated {
        regime: regime.clone(),
        dataset: "dataset-0".to_string(),
    };
    for index in 1..1_024 {
        validated = compose(
            &validated,
            &Color::Validated {
                regime: regime.clone(),
                dataset: format!("dataset-{index}"),
            },
            IntervalOp::Hull,
        );
        let Color::Validated { dataset, .. } = &validated else {
            panic!("overlapping validated regimes must remain validated")
        };
        assert!(dataset.len() <= MAX_COLOR_IDENTITY_BYTES);
    }

    let cards = (0..1_024)
        .map(|index| format!("model-card-{index}"))
        .collect::<Vec<_>>();
    let modeled = ModelEvidence {
        cards,
        assumptions: Vec::new(),
        validity: regime,
        discrepancy_rel: 0.0,
        in_domain: true,
    };
    let Color::Estimated {
        estimator: dataset,
        dispersion,
    } = color_of(&NumericalCertificate::enclosure(0.0, 1.0), &modeled)
    else {
        panic!("bounded model cards remain Estimated without an authenticated anchor")
    };
    assert!(dataset.len() <= MAX_COLOR_IDENTITY_BYTES);
    assert_eq!(dispersion.to_bits(), 0.0_f64.to_bits());

    let one_long_card = ModelEvidence {
        cards: vec![overlong],
        assumptions: Vec::new(),
        validity: ValidityDomain::unconstrained(),
        discrepancy_rel: 0.0,
        in_domain: true,
    };
    let Color::Estimated { estimator, .. } =
        color_of(&NumericalCertificate::estimate(0.0, 1.0), &one_long_card)
    else {
        panic!("one unbounded model-card identity remains Estimated")
    };
    assert!(estimator.len() <= MAX_COLOR_IDENTITY_BYTES);
    verdict(
        "evd-014",
        true,
        "repeated and heterogeneous estimator composition stays bounded, deterministic, and Estimated",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_015_malformed_bridge_inputs_fail_closed() {
    use fs_evidence::{Color, ModelEvidence, NumericalCertificate, ValidityDomain, color_of};

    let no_claim = color_of(&NumericalCertificate::no_claim(), &ModelEvidence::none());
    assert!(matches!(
        no_claim,
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
    let malformed_exact = color_of(
        &NumericalCertificate::exact(f64::NAN),
        &ModelEvidence::none(),
    );
    assert!(matches!(
        malformed_exact,
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
    let uncarded_model_spread = ModelEvidence {
        discrepancy_rel: 0.2,
        ..ModelEvidence::none()
    };
    assert!(matches!(
        color_of(
            &NumericalCertificate::enclosure(0.0, 1.0),
            &uncarded_model_spread,
        ),
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));

    let regime = ValidityDomain::unconstrained().with("x", 0.0, 1.0);
    let out_of_domain = ModelEvidence {
        cards: vec!["card-a".to_string()],
        assumptions: Vec::new(),
        validity: regime.clone(),
        discrepancy_rel: 0.1,
        in_domain: false,
    };
    assert!(matches!(
        color_of(
            &NumericalCertificate::enclosure(0.0, 1.0),
            &out_of_domain,
        ),
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
    let malformed_carded = ModelEvidence {
        in_domain: true,
        ..out_of_domain.clone()
    };
    assert!(matches!(
        color_of(
            &NumericalCertificate::exact(f64::NAN),
            &malformed_carded,
        ),
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
    let canonical_cards = ModelEvidence {
        cards: vec!["zeta".to_string(), "alpha".to_string(), "zeta".to_string()],
        assumptions: Vec::new(),
        validity: regime.clone(),
        discrepancy_rel: 0.0,
        in_domain: true,
    };
    let sorted_cards = ModelEvidence {
        cards: vec!["alpha".to_string(), "zeta".to_string()],
        ..canonical_cards.clone()
    };
    assert_eq!(
        color_of(&NumericalCertificate::enclosure(0.0, 1.0), &canonical_cards,),
        color_of(&NumericalCertificate::enclosure(0.0, 1.0), &sorted_cards,),
        "model-card identity is a canonical sorted set"
    );

    let model_spread = ModelEvidence {
        cards: vec!["identified-model".to_string()],
        discrepancy_rel: 0.2,
        ..ModelEvidence::none()
    };
    let Color::Estimated { dispersion, .. } =
        color_of(&NumericalCertificate::estimate(0.9, 1.1), &model_spread)
    else {
        panic!("an uncarded numerical estimate remains Estimated")
    };
    assert!(
        (dispersion - 0.3).abs() < 1e-12,
        "numerical and model relative spreads add conservatively: {dispersion}"
    );
    verdict(
        "evd-015",
        true,
        "malformed and out-of-domain bridge inputs demote to infinite-dispersion Estimated evidence",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_015b_malformed_color_payloads_fail_closed() {
    use fs_evidence::{
        Color, ColorRank, IntervalOp, ValidityDomain, compose, validate_color_payload,
    };

    let valid = Color::Verified { lo: 1.0, hi: 2.0 };
    let malformed = [
        Color::Verified {
            lo: f64::NAN,
            hi: f64::NAN,
        },
        Color::Verified { lo: 2.0, hi: 1.0 },
        Color::Estimated {
            estimator: "bad-dispersion".to_string(),
            dispersion: -1.0,
        },
        Color::Estimated {
            estimator: "bad-dispersion".to_string(),
            dispersion: f64::NAN,
        },
        Color::Validated {
            regime: ValidityDomain::unconstrained(),
            dataset: "dataset-a".to_string(),
        },
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("x", 0.0, f64::INFINITY),
            dataset: "dataset-a".to_string(),
        },
    ];
    for bad in malformed {
        assert!(validate_color_payload(&bad).is_err());
        for op in [IntervalOp::Add, IntervalOp::Mul, IntervalOp::Hull] {
            let result = compose(&bad, &valid, op);
            assert_eq!(result.rank(), ColorRank::Estimated);
            assert!(matches!(
                result,
                Color::Estimated { dispersion, .. } if dispersion.is_infinite()
            ));
        }
    }

    let whole_line = Color::Verified {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };
    validate_color_payload(&whole_line).expect("whole-line enclosure is sound but vacuous");
    verdict(
        "evd-015b",
        true,
        "malformed color payloads fail closed while the sound whole-line enclosure remains valid",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn evd_015c_empty_regime_demotion_has_a_valid_bounded_identity() {
    use fs_evidence::{
        Color, ValidityDomain, check_regime, color_identity_reason, demotion_estimator_identity,
        validate_color_payload,
    };
    use std::collections::BTreeMap;

    let malformed = Color::Validated {
        regime: ValidityDomain::unconstrained(),
        dataset: "fixture-dataset".to_string(),
    };
    let (demoted, record) = check_regime(&malformed, &BTreeMap::new());
    let Color::Estimated {
        estimator,
        dispersion,
    } = &demoted
    else {
        panic!("an undeclared regime must demote")
    };
    assert!(dispersion.is_infinite());
    assert!(record.is_some());
    assert_eq!(
        estimator,
        &demotion_estimator_identity("fixture-dataset", "<undeclared-regime>")
    );
    assert_eq!(color_identity_reason(estimator), None);
    validate_color_payload(&demoted).expect("demotion must always emit a valid color payload");
}

#[test]
fn evd_015c_malformed_card_diagnostics_remain_derived() {
    use fs_evidence::{
        Color, ModelEvidence, NumericalCertificate, color_leaf_identity_reason, color_of,
    };

    let reserved_card = ModelEvidence {
        cards: vec!["derived:v2:forged-model-card".to_string()],
        ..ModelEvidence::none()
    };
    let Color::Estimated {
        estimator,
        dispersion,
    } = color_of(&NumericalCertificate::enclosure(0.0, 1.0), &reserved_card)
    else {
        panic!("a malformed model-card identity must fail closed as Estimated")
    };
    assert!(dispersion.is_infinite());
    assert!(estimator.starts_with("derived:v2:invalid-card-"));
    assert_eq!(
        color_leaf_identity_reason(&estimator),
        Some("derived-identity-requires-lineage"),
        "generated malformed-card diagnostics stay in the reserved derived namespace"
    );

    let diagnostic = |cards: &[&str]| {
        let modeled = ModelEvidence {
            cards: cards.iter().map(|card| (*card).to_string()).collect(),
            ..ModelEvidence::none()
        };
        let Color::Estimated {
            estimator,
            dispersion,
        } = color_of(&NumericalCertificate::enclosure(0.0, 1.0), &modeled)
        else {
            panic!("a malformed card set must fail closed as Estimated")
        };
        assert!(dispersion.is_infinite());
        estimator
    };
    let suffix_a = diagnostic(&["pending", "A"]);
    let suffix_b = diagnostic(&["pending", "B"]);
    assert_ne!(
        suffix_a, suffix_b,
        "the complete malformed card set must participate in diagnostic identity"
    );
    assert_eq!(
        suffix_a,
        diagnostic(&["A", "pending", "A"]),
        "card identities have deterministic set semantics across order and duplicates"
    );
    verdict(
        "evd-015c",
        true,
        "malformed card diagnostics bind the complete canonical card set, remain derived, and \n         cannot be re-rooted as evidence leaves",
        FIXED_INPUT_SEED,
    );
}

// ---------------------------------------------------------------------------
// G0 property adoption (bead frankensim-4nh8): the color-rank lattice law,
// generated + shrunk via fs-propcheck. Fixed color cases above remain pins.
// ---------------------------------------------------------------------------

type ColorRecipe = (u64, (i64, i64));
type ColorCompositionCase = (ColorRecipe, ColorRecipe, u64);

fn generate_color_composition_case(stream: &mut Stream) -> ColorCompositionCase {
    (
        (
            stream.next_u64() % 3,
            (stream.int_in(-8, 8), stream.int_in(-8, 8)),
        ),
        (
            stream.next_u64() % 3,
            (stream.int_in(-8, 8), stream.int_in(-8, 8)),
        ),
        stream.next_u64() % 3,
    )
}

fn color_from_recipe(recipe: &ColorRecipe, side: &str) -> Color {
    let &(kind, (a, b)) = recipe;
    let lo = a.min(b) as f64;
    let hi = a.max(b) as f64;
    match kind % 3 {
        0 => Color::Verified { lo, hi },
        1 => Color::Validated {
            regime: ValidityDomain::unconstrained().with("x", lo, hi),
            dataset: format!("property-dataset-{side}"),
        },
        _ => Color::Estimated {
            estimator: format!("property-estimator-{side}"),
            dispersion: a.unsigned_abs() as f64 + b.unsigned_abs() as f64,
        },
    }
}

fn interval_op_from_recipe(recipe: u64) -> IntervalOp {
    match recipe % 3 {
        0 => IntervalOp::Add,
        1 => IntervalOp::Mul,
        _ => IntervalOp::Hull,
    }
}

fn validated_pair_is_disjoint(left: &Color, right: &Color) -> bool {
    match (left, right) {
        (Color::Validated { regime: a, .. }, Color::Validated { regime: b, .. }) => {
            a.intersect(b).is_empty()
        }
        _ => false,
    }
}

fn color_composition_rank_law(case: &ColorCompositionCase) -> bool {
    let left = color_from_recipe(&case.0, "left");
    let right = color_from_recipe(&case.1, "right");
    let op = interval_op_from_recipe(case.2);
    let weakest = left.rank().min(right.rank());
    let expected = if validated_pair_is_disjoint(&left, &right) {
        ColorRank::Estimated
    } else {
        weakest
    };
    let forward = compose(&left, &right, op);
    let reverse = compose(&right, &left, op);

    validate_color_payload(&left).is_ok()
        && validate_color_payload(&right).is_ok()
        && validate_color_payload(&forward).is_ok()
        && validate_color_payload(&reverse).is_ok()
        && forward.rank() <= weakest
        && reverse.rank() <= weakest
        && forward.rank() == expected
        && reverse.rank() == expected
}

#[test]
fn evd_016_g0_color_composition_rank_never_launders() {
    const SEED: u64 = 0xE71D_4A48_0001;
    const CASES: u64 = 512;

    // Keep the fixed seed honest: it covers the full 3x3x3 kind/op matrix and
    // both branches of Validated + Validated composition.
    let mut combinations = BTreeSet::new();
    let mut saw_validated_overlap = false;
    let mut saw_validated_disjoint = false;
    for case_index in 0..CASES {
        let mut stream = Stream::for_case(SEED, case_index);
        let case = generate_color_composition_case(&mut stream);
        combinations.insert((case.0.0 % 3, case.1.0 % 3, case.2 % 3));

        let left = color_from_recipe(&case.0, "left");
        let right = color_from_recipe(&case.1, "right");
        if let (Color::Validated { .. }, Color::Validated { .. }) = (&left, &right) {
            if validated_pair_is_disjoint(&left, &right) {
                saw_validated_disjoint = true;
            } else {
                saw_validated_overlap = true;
            }
        }
    }
    assert_eq!(
        combinations.len(),
        27,
        "fixed seed must cover all operand-kind/operation combinations"
    );
    assert!(
        saw_validated_overlap && saw_validated_disjoint,
        "fixed seed must cover both validated-regime branches"
    );

    fs_propcheck::check(
        "color-composition-rank-never-launders",
        SEED,
        CASES,
        generate_color_composition_case,
        color_composition_rank_law,
    );
    verdict(
        "evd-016",
        true,
        "512 generated shrinkable compositions cover all color-kind/operation combinations, both operand orders, and both validated-regime branches",
        SEED,
    );
}
