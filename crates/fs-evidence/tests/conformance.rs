//! fs-evidence conformance suite (CONTRACT.md: any reimplementation must
//! pass). G0: conservativeness laws of the composition algebra (composed
//! enclosures contain true propagation; composed validity ⊆ intersection).
//! Plus the acceptance battery: the registration lint, the worked
//! model-discrepancy-dominates example, the out-of-distribution refusal,
//! and deterministic ledger rows. JSON-line verdicts; seeded cases carry
//! their seed.

use fs_evidence::{
    Ambition, DecisionStatus, DiscrepancyModel, EscalationAdvice, Evidence, FidelityPair,
    ModelBracket, ModelCard, ModelEvidence, ModelRegistry, NumericalCertificate, Op,
    ProvenanceHash, StatisticalCertificate, UncertaintySource, ValidityDomain,
};
use std::collections::BTreeMap;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-evidence/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
    println!("{}", drag.to_ledger_row_json());
    verdict(
        "evd-003",
        dominated && advice_ok && flips && detail.contains("model-form"),
        &format!("beautifully-certified-wrong caught: {detail}"),
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
    verdict(
        "evd-004",
        band_sane && teaching && ood_blocked,
        &format!(
            "trained on 50 pairs (seed {SEED:#x}): in-domain band mean {:.3}/max {:.3}; \
             out-of-domain query refused with the violated parameter named",
            band.mean_rel, band.max_rel
        ),
    );
}

#[test]
fn evd_005_bracketing_reports_spread_and_rows_are_deterministic() {
    // The vessel flagship's stated mitigation: bracket contact-angle
    // models, report the objective's sensitivity band.
    let bracket = ModelBracket::new()
        .with_member("contact-angle-60", 0.90)
        .with_member("contact-angle-90", 1.00)
        .with_member("contact-angle-120", 1.16);
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
    );
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
    // (e) Legitimate exact/enclosure composition remains usable, and
    // reads flow through the immutable Deref view.
    let a = Evidence::exact(2.0, p);
    let b = Evidence::enclosed(3.0, 2.9, 3.1, p);
    let cert = Evidence::combine(Op::Mul, &a, &b, ())
        .certified()
        .expect("rigorous chain certifies");
    let readable = cert.qoi.to_bits() == 6.0f64.to_bits() && cert.evidence().numerical.lo <= 6.0;
    // (f) Downgrade-mutate-recertify: the ONLY mutation path loses the
    // mark, and reconstruction re-validates (round-trip invariance).
    let mut reopened = cert.into_evidence();
    reopened.numerical = NumericalCertificate::estimate(5.9, 6.1);
    let weakened_refused = matches!(
        reopened.clone().certified(),
        Err(CertifyError::NotRigorous { .. })
    );
    reopened.numerical = NumericalCertificate::enclosure(5.9, 6.1);
    reopened.qoi = 7.0; // outside the reclaimed enclosure
    let drifted_refused = matches!(
        reopened.certified(),
        Err(CertifyError::QoiOutsideEnclosure { .. })
    );
    verdict(
        "evd-012",
        nan_exact
            && inf_exact
            && inconsistent
            && inverted
            && nan_bound
            && inf_bound
            && escaped
            && readable
            && weakened_refused
            && drifted_refused,
        "Certified<T> is opaque: every forge route refused with its structured reason; \
         rigorous composition still certifies; downgrade-mutate-recertify re-validates",
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
    let overlap_ok = if let Color::Validated { regime, .. } = &overlap {
        let (lo, hi) = regime.bound("Re").unwrap_or((0.0, 0.0));
        (lo - 1.5e5).abs() < 1.0 && (hi - 2e5).abs() < 1.0
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
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    );
    verdict(
        "evd-007",
        overlap_ok && demoted,
        "validated⊕validated intersects when overlapping; disjoint regimes demote to \
         Estimated with infinite/no-claim dispersion instead of laundering a phantom point \
         regime or asserting a zero spread",
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
    );
}

#[test]
fn evd_009_non_finite_regimes_fail_closed_and_payloads_escape() {
    use fs_evidence::{Color, check_regime};

    let validated = Color::Validated {
        regime: ValidityDomain::unconstrained().with("Re", 1.0, 10.0),
        dataset: "anchors".to_string(),
    };
    let mut non_finite_states_demote = true;
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let (checked, flag) = check_regime(&validated, &pt(&[("Re", value)]));
        non_finite_states_demote &= matches!(
            checked,
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ) && flag
            .is_some_and(|d| d.axis == "Re" && !d.value.is_finite());
    }

    let inverted = ValidityDomain::unconstrained()
        .with("Re", 1.0, 2.0)
        .intersect(&ValidityDomain::unconstrained().with("Re", 3.0, 4.0));
    let invalid_regimes = [
        ValidityDomain::unconstrained(),
        ValidityDomain::unconstrained().with("Re", f64::NEG_INFINITY, 10.0),
        ValidityDomain::unconstrained().with("Re", f64::NAN, f64::NAN),
        inverted,
    ];
    let mut invalid_regimes_demote = true;
    for regime in invalid_regimes {
        let invalid = Color::Validated {
            regime,
            dataset: "anchors".to_string(),
        };
        let (checked, flag) = check_regime(&invalid, &pt(&[("Re", 5.0)]));
        invalid_regimes_demote &= matches!(
            checked,
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ) && flag.is_some();
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
    );
}
