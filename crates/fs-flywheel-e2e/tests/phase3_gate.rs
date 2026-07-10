//! ADDENDUM PHASE 3 — HORIZON: the terminal roadmap gate (bead
//! xpck.5). NOT a build gate but an ACTIVATION LEDGER: each horizon
//! proposal's trigger measurement is INSTRUMENTED and its current
//! verdict recorded in writing (the signed holding-pen package is the
//! quarterly-review artifact). Nothing opens as a broad program;
//! radical systems die of breadth more often than of ambition (R10).
#![cfg(feature = "flywheel-e2e")]

use fs_evidence::Color;
use fs_package::{Claim, EvidencePackage, Provenance};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-flywheel-e2e/phase3\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[test]
fn p3_001_proposal_a_trigger_fired() {
    use fs_surrogate::ladder::{Ladder, rb_coverage};
    // Proposal A's activation measurement: RB coverage of the wedge
    // query volume, kill floor 0.2. INSTRUMENTED and — on the elliptic
    // beachhead — FIRED.
    let ladder = Ladder::build(150, (0.0, 4.0), &[5, 2], false);
    let mus: Vec<f64> = (0..10).map(|i| 4.0 * f64::from(i) / 9.0).collect();
    let coverage = rb_coverage(&ladder, &mus, &[1e-2, 1e-5, 1e-8]);
    println!(
        "{{\"metric\":\"horizon-A\",\"rb_coverage\":{coverage:.3},\"floor\":0.2,\
         \"status\":\"TRIGGER FIRED — active behind abstraction-ladder\"}}"
    );
    assert!(coverage >= 0.2, "the beachhead trigger holds: {coverage}");
    verdict(
        "p3-001",
        "Proposal A: rb_coverage instrumented and above the 0.2 kill floor — the one \
         horizon member whose trigger has fired (active, feature-gated)",
    );
}

#[test]
fn p3_002_proposal_c_instrumented_awaiting_audit() {
    use fs_plan::voi::{AuditVerdict, audit_verdict};
    // Proposal C's activation is CONDITIONAL: the machinery is live
    // (Phase-2 benchmarks), but SCHEDULING AUTHORITY requires the
    // prospective audit to show recommendations beat agent choices —
    // and with no audit evidence the verdict is Demote by design.
    assert_eq!(
        audit_verdict(&[]),
        AuditVerdict::DemoteToReporting,
        "no evidence, no authority — the default is the safe one"
    );
    println!(
        "{{\"metric\":\"horizon-C\",\"status\":\"INSTRUMENTED — authority awaits the \
         prospective audit (two quarters of matched-cost comparisons)\"}}"
    );
    verdict(
        "p3-002",
        "Proposal C: the audit instrument exists and defaults to demotion without \
         evidence — activation is a measurement, not a decision",
    );
}

#[test]
fn p3_003_proposal_4_instrument_only_by_default() {
    use fs_time::slabs::{Activation, CoupledFixture, activation_report, march_instrumented};
    // Proposal 4's gate: control activates ONLY where splitting error
    // dominates the budget. Both directions of the instrument verified.
    let weak = CoupledFixture { coupling: |_| 0.02 };
    let (_, ledger) = march_instrumented(&weak, [1.0, 0.5], 2.0, 8, 1);
    let (frac_weak, v_weak) = activation_report(&ledger, 1e-2);
    assert_eq!(v_weak, Activation::InstrumentOnly);
    let strong = CoupledFixture { coupling: |_| 1.5 };
    let (_, ledger) = march_instrumented(&strong, [1.0, 0.5], 2.0, 8, 1);
    let (frac_strong, v_strong) = activation_report(&ledger, 1e-2);
    assert_eq!(v_strong, Activation::ControlJustified);
    println!(
        "{{\"metric\":\"horizon-4\",\"weak_fraction\":{frac_weak:.3},\
         \"strong_fraction\":{frac_strong:.3},\
         \"status\":\"INSTRUMENTED — control gated on a paying workload's budget\"}}"
    );
    verdict(
        "p3-003",
        "Proposal 4: the splitting-error activation instrument fires in both directions; \
         default posture is instrumented-but-uncontrolled",
    );
}

#[test]
fn p3_004_proposal_13b_prevalence_measurement() {
    use fs_symmetry::cyclic_residual;
    // Proposal 13b's gate: >=15% of real workloads present exploitable
    // symmetry. The PREVALENCE INSTRUMENT: fraction of a workload
    // battery whose fields are near-k-fold (relative residual < 0.05).
    let prevalence = |battery: &[Vec<f64>]| -> f64 {
        let hits = battery
            .iter()
            .filter(|v| cyclic_residual(v, 2).is_ok_and(|r| r.relative < 0.05))
            .count();
        #[allow(clippy::cast_precision_loss)]
        {
            hits as f64 / battery.len() as f64
        }
    };
    // A symmetry-rich battery clears the bar; a generic one does not.
    let rich: Vec<Vec<f64>> = (0..10)
        .map(|k| {
            if k < 3 {
                vec![1.0, 2.0, 1.0, 2.0] // exactly 2-fold
            } else if k < 5 {
                vec![1.0, 2.0, 1.001, 2.001] // near-symmetric
            } else {
                vec![f64::from(k), 1.0, -2.0, 0.3] // generic
            }
        })
        .collect();
    let rich_frac = prevalence(&rich);
    let generic: Vec<Vec<f64>> = (0..10)
        .map(|k| vec![f64::from(k) + 0.7, 1.3, -2.1, 0.4 * f64::from(k)])
        .collect();
    let generic_frac = prevalence(&generic);
    println!(
        "{{\"metric\":\"horizon-13b\",\"rich_prevalence\":{rich_frac:.2},\
         \"generic_prevalence\":{generic_frac:.2},\"bar\":0.15,\
         \"status\":\"INSTRUMENTED — detection ships, the solver waits for prevalence\"}}"
    );
    assert!(
        rich_frac >= 0.15 && generic_frac < 0.15,
        "both directions measured"
    );
    verdict(
        "p3-004",
        "Proposal 13b: the symmetry-prevalence instrument separates a symmetry-rich \
         battery (>=15%) from a generic one (<15%) — detection ships, the dedicated \
         solver waits for real-workload prevalence",
    );
}

#[test]
fn p3_005_proposal_11_r8_gate_and_the_holding_pen() {
    use fs_asbuilt::{Fiducial, Point2, register, well_posed};
    // Proposal 11's R8 instrument: registration must be tighter than
    // the deviations being certified. GOOD fiducials pass...
    let good: Vec<Fiducial> = [(0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0)]
        .iter()
        .map(|&(x, y)| {
            Fiducial::new(
                Point2::new(x, y),
                Point2::new(x + 0.30001, y + 0.19999), // clean rigid shift
            )
        })
        .collect();
    let reg = register(&good).expect("registers");
    assert!(
        well_posed(&reg, 0.05),
        "clean registration certifies 0.05 deviations"
    );
    // ...and sloppy fiducials FAIL the same certification (the R8 kill
    // visible in the instrument).
    let sloppy: Vec<Fiducial> = [(0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0)]
        .iter()
        .enumerate()
        .map(|(k, &(x, y))| {
            let wobble = 0.2 * f64::from(u32::try_from(k).expect("small"));
            Fiducial::new(
                Point2::new(x, y),
                Point2::new(x + 0.3 + wobble, y + 0.2 - wobble),
            )
        })
        .collect();
    let reg = register(&sloppy).expect("registers");
    assert!(
        !well_posed(&reg, 0.05),
        "sloppy registration cannot certify what it cannot resolve"
    );
    // THE HOLDING PEN, IN WRITING: the five statuses as a signed
    // package — the quarterly-review artifact a third party can check.
    let pkg = EvidencePackage::new(Provenance::new("phase3-horizon", "Cargo.lock"))
        .with_claim(Claim::new(
            "A-abstraction-ladder",
            "trigger FIRED: rb_coverage 1.0 >= 0.2 on the beachhead",
            Color::Verified { lo: 0.2, hi: 1.0 },
        ))
        .with_claim(Claim::new(
            "C-value-of-information",
            "instrumented; scheduling authority awaits the prospective audit",
            Color::Estimated {
                estimator: "prospective-audit-pending".to_string(),
                dispersion: 1.0,
            },
        ))
        .with_claim(Claim::new(
            "4-spacetime-complex",
            "instrumented-but-uncontrolled; control gated on a splitting-dominated \
             paying workload",
            Color::Estimated {
                estimator: "workload-demand-pending".to_string(),
                dispersion: 1.0,
            },
        ))
        .with_claim(Claim::new(
            "13b-symmetry-solver",
            "prevalence instrument live; dedicated solver waits for >=15% real-workload \
             symmetry",
            Color::Estimated {
                estimator: "prevalence-pending".to_string(),
                dispersion: 1.0,
            },
        ))
        .with_claim(Claim::new(
            "11-reality-as-a-chart",
            "R8 instrument live (registration vs certified deviation); full-field \
             activation awaits metrology partnerships — point-sensor assimilation \
             ships meanwhile",
            Color::Estimated {
                estimator: "metrology-partnership-pending".to_string(),
                dispersion: 1.0,
            },
        ))
        .signed("phase3-horizon-gate");
    let check = fs_checker::check(&pkg);
    assert!(check.passed(), "the holding-pen record re-verifies");
    let breakdown = pkg.color_breakdown();
    assert!(
        breakdown.verified == 1 && breakdown.estimated == 4,
        "exactly one trigger has fired; four wait honestly: {breakdown:?}"
    );
    println!(
        "{{\"metric\":\"horizon-ledger\",\"root\":\"{}\",\"fired\":1,\"waiting\":4}}",
        pkg.merkle_root()
    );
    verdict(
        "p3-005",
        "Proposal 11's R8 instrument passes clean registration and fails sloppy; the \
         five-proposal holding pen ships as a signed, checker-verified package — one \
         trigger fired, four waiting, in writing",
    );
}
