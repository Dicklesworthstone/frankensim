//! G1 existing finite-window reed onset; G0 type/provenance boundary.
#[allow(dead_code)]
#[path = "../examples/midi_reed.rs"]
mod fixture;
use fs_couple::description::per_note::*;
use fs_qty::QtyAny;
use fs_scenario::performance::*;

#[test]
fn pn_001_threshold_objective_resolves_on_model() {
    let (d, p, t) = fixture::fixture();
    let perf = import_smf(&fixture::authored_midi(), &d, &p, &t, fixture::limits()).unwrap();
    let c = &perf.per_note()[0];
    let budget = fixture::budget();
    let a = resolve_threshold(&mut fixture::ReedOnset, c, budget).unwrap();
    let b = resolve_threshold(&mut fixture::ReedOnset, c, budget).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.content_hash(), b.content_hash());
    assert!(a.bracket_pa.1 - a.bracket_pa.0 <= budget.pressure_tolerance_pa);
    assert!(!fixture::ReedOnset.oscillates(a.bracket_pa.0).unwrap());
    assert!(fixture::ReedOnset.oscillates(a.bracket_pa.1).unwrap());
    assert_eq!(a.resolved(), QtyAny::new(a.bracket_pa.1, fixture::PA));
    assert!(a.probes.len() <= budget.max_evaluations);
    // Independent hand bisection using exactly the same declared physical
    // observable; no synthetic threshold or recording-derived oracle.
    let (mut lo, mut hi) = (0.0, 4800.0);
    for _ in 0..8 {
        let mid = (lo + hi) * 0.5;
        if fixture::ReedOnset.oscillates(mid).unwrap() {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    assert_eq!(a.bracket_pa, (lo, hi));
}
#[test]
fn pn_002_resolver_signature_has_no_reference_input() {
    type ModelOnly = fn(
        &mut fixture::ReedOnset,
        &PerNoteControl,
        ThresholdBudget,
    ) -> Result<ThresholdReceipt, ResolutionError>;
    let _api: ModelOnly = resolve_threshold::<fixture::ReedOnset>;
    let (d, p, t) = fixture::fixture();
    let perf = import_smf(&fixture::authored_midi(), &d, &p, &t, fixture::limits()).unwrap();
    let mut control = perf.per_note()[0].clone();
    control.objective = DeclaredObjective::SlotCentre;
    assert_eq!(
        resolve_threshold(&mut fixture::ReedOnset, &control, fixture::budget()),
        Err(ResolutionError::UnsupportedObjective)
    );
}
#[test]
fn pn_003_receipt_recorded_in_provenance() {
    let (d, p, t) = fixture::fixture();
    let mut perf = import_smf(&fixture::authored_midi(), &d, &p, &t, fixture::limits()).unwrap();
    assert!(matches!(
        lower(&perf, &d, &p, &t),
        Err(PerformanceError::UnresolvedControl { .. })
    ));
    for i in 0..perf.per_note().len() {
        let r = resolve_threshold(
            &mut fixture::ReedOnset,
            &perf.per_note()[i],
            fixture::budget(),
        )
        .unwrap();
        record_resolution(&mut perf, i, &r).unwrap();
        assert_eq!(perf.per_note()[i].receipt, Some(r.content_hash()));
    }
    let replay = Performance::from_bytes(&perf.to_bytes(), &d, &p, &t, fixture::limits()).unwrap();
    let a = lower(&perf, &d, &p, &t).unwrap();
    let b = lower(&replay, &d, &p, &t).unwrap();
    assert_eq!(a.provenance.per_note, perf.per_note());
    assert_eq!(a.provenance.per_note, b.provenance.per_note);
    assert_eq!(a.provenance.performance, b.provenance.performance);
}
#[test]
fn pn_004_search_budget_units_and_bracket_refuse() {
    struct Model;
    impl ThresholdModel for Model {
        fn identity(&self) -> fs_blake3::ContentHash {
            fs_blake3::hash_bytes(b"G0 linear monotone oracle pressure=2000Pa")
        }
        fn oscillates(&mut self, p: f64) -> Result<bool, String> {
            Ok(p >= 2000.0)
        }
    }
    let (d, p, t) = fixture::fixture();
    let perf = import_smf(&fixture::authored_midi(), &d, &p, &t, fixture::limits()).unwrap();
    let mut c = perf.per_note()[0].clone();
    assert_eq!(
        resolve_threshold(
            &mut Model,
            &c,
            ThresholdBudget {
                max_evaluations: 2,
                ..fixture::budget()
            }
        ),
        Err(ResolutionError::BudgetExceeded)
    );
    c.domain.0 = QtyAny::new(2100.0, fixture::PA);
    assert_eq!(
        resolve_threshold(&mut Model, &c, fixture::budget()),
        Err(ResolutionError::Unbracketed)
    );
    c.domain.0 = QtyAny::dimensionless(0.0);
    assert!(matches!(
        resolve_threshold(&mut Model, &c, fixture::budget()),
        Err(ResolutionError::Invalid(_))
    ));
}
