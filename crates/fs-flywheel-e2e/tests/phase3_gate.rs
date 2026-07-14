//! ADDENDUM PHASE 3 — HORIZON: the terminal roadmap gate (bead
//! xpck.5). NOT a build gate but an ACTIVATION LEDGER: each horizon
//! proposal's trigger measurement is INSTRUMENTED and its current
//! verdict recorded in writing (the fixture-authenticated holding-pen package is the
//! quarterly-review artifact). Nothing opens as a broad program;
//! radical systems die of breadth more often than of ambition (R10).
#![cfg(feature = "flywheel-e2e")]

use fs_package::{Claim, EvidencePackage, Provenance};

fn with_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x7068_6173_6533_0005,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

const HORIZON_COVERAGE_SCHEMA_VERSION: u64 = 1;
const HORIZON_COVERAGE_ALGORITHM: &str = "fs-surrogate.rb-coverage/one-solve-per-mu-v2";
const HORIZON_COVERAGE_MODEL: &str = "fs-surrogate.truth-model/p1-piecewise-affine-elliptic-v1";
const HORIZON_RB_DIMENSIONS: [usize; 2] = [5, 2];
const HORIZON_TOLERANCES: [f64; 3] = [1e-2, 1e-5, 1e-8];

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-flywheel-e2e/phase3\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[derive(Debug, Clone)]
struct HorizonCoverageSpec {
    schema_version: u64,
    algorithm: &'static str,
    model: &'static str,
    truth_nodes: usize,
    mu_range: (f64, f64),
    rb_dimensions: Vec<usize>,
    concept: bool,
    mus: Vec<f64>,
    tolerances: Vec<f64>,
}

impl HorizonCoverageSpec {
    fn canonical_bytes(&self) -> Vec<u8> {
        fn push_atom(out: &mut Vec<u8>, value: &[u8]) {
            let len = u64::try_from(value.len()).expect("fixture atom length fits u64");
            out.extend_from_slice(&len.to_le_bytes());
            out.extend_from_slice(value);
        }

        let mut out = Vec::new();
        push_atom(
            &mut out,
            b"frankensim/fs-flywheel-e2e/horizon-coverage-spec/v1",
        );
        push_atom(&mut out, &self.schema_version.to_le_bytes());
        push_atom(&mut out, self.algorithm.as_bytes());
        push_atom(&mut out, self.model.as_bytes());
        push_atom(
            &mut out,
            &u64::try_from(self.truth_nodes)
                .expect("fixture truth dimension fits u64")
                .to_le_bytes(),
        );
        push_atom(&mut out, &self.mu_range.0.to_bits().to_le_bytes());
        push_atom(&mut out, &self.mu_range.1.to_bits().to_le_bytes());
        push_atom(
            &mut out,
            &u64::try_from(self.rb_dimensions.len())
                .expect("fixture rung count fits u64")
                .to_le_bytes(),
        );
        for &dimension in &self.rb_dimensions {
            push_atom(
                &mut out,
                &u64::try_from(dimension)
                    .expect("fixture rung dimension fits u64")
                    .to_le_bytes(),
            );
        }
        push_atom(&mut out, &[u8::from(self.concept)]);
        push_atom(
            &mut out,
            &u64::try_from(self.mus.len())
                .expect("fixture parameter count fits u64")
                .to_le_bytes(),
        );
        for &mu in &self.mus {
            push_atom(&mut out, &mu.to_bits().to_le_bytes());
        }
        push_atom(
            &mut out,
            &u64::try_from(self.tolerances.len())
                .expect("fixture tolerance count fits u64")
                .to_le_bytes(),
        );
        for &tolerance in &self.tolerances {
            push_atom(&mut out, &tolerance.to_bits().to_le_bytes());
        }
        out
    }

    fn canonical_text(&self) -> String {
        let dimensions = self
            .rb_dimensions
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let mus = self
            .mus
            .iter()
            .map(|value| format!("{:016x}", value.to_bits()))
            .collect::<Vec<_>>()
            .join(",");
        let tolerances = self
            .tolerances
            .iter()
            .map(|value| format!("{:016x}", value.to_bits()))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "schema={};algorithm={};model={};truth_nodes={};mu_range_bits=[{:016x},{:016x}];\
             rb_dimensions=[{dimensions}];concept={};mus_bits=[{mus}];\
             tolerance_bits=[{tolerances}]",
            self.schema_version,
            self.algorithm,
            self.model,
            self.truth_nodes,
            self.mu_range.0.to_bits(),
            self.mu_range.1.to_bits(),
            self.concept,
        )
    }

    fn hash_hex(&self) -> String {
        fs_ledger::hash_bytes(&self.canonical_bytes()).to_hex()
    }
}

fn horizon_coverage_spec() -> HorizonCoverageSpec {
    HorizonCoverageSpec {
        schema_version: HORIZON_COVERAGE_SCHEMA_VERSION,
        algorithm: HORIZON_COVERAGE_ALGORITHM,
        model: HORIZON_COVERAGE_MODEL,
        truth_nodes: 150,
        mu_range: (0.0, 4.0),
        rb_dimensions: HORIZON_RB_DIMENSIONS.to_vec(),
        concept: false,
        mus: (0..10).map(|i| 4.0 * f64::from(i) / 9.0).collect(),
        tolerances: HORIZON_TOLERANCES.to_vec(),
    }
}

#[derive(Debug)]
struct HorizonCoverageMeasurement {
    spec: HorizonCoverageSpec,
    coverage: f64,
}

fn measure_horizon_coverage() -> HorizonCoverageMeasurement {
    let spec = horizon_coverage_spec();
    let ladder = fs_surrogate::ladder::Ladder::build(
        spec.truth_nodes,
        spec.mu_range,
        &spec.rb_dimensions,
        spec.concept,
    )
    .expect("bounded horizon ladder");
    let coverage = fs_surrogate::ladder::rb_coverage(&ladder, &spec.mus, &spec.tolerances)
        .expect("bounded horizon coverage battery");
    HorizonCoverageMeasurement { spec, coverage }
}

#[test]
fn p3_000_coverage_spec_identity_binds_every_semantic_input() {
    let base = horizon_coverage_spec();
    let base_hash = base.hash_hex();
    let mut variants = Vec::new();

    let mut changed = base.clone();
    changed.schema_version += 1;
    variants.push(("schema", changed));
    let mut changed = base.clone();
    changed.algorithm = "fs-surrogate.rb-coverage/alternate";
    variants.push(("algorithm", changed));
    let mut changed = base.clone();
    changed.model = "fs-surrogate.truth-model/alternate";
    variants.push(("model", changed));
    let mut changed = base.clone();
    changed.truth_nodes += 1;
    variants.push(("truth dimension", changed));
    let mut changed = base.clone();
    changed.mu_range.1 = f64::from_bits(changed.mu_range.1.to_bits() + 1);
    variants.push(("parameter range", changed));
    let mut changed = base.clone();
    changed.rb_dimensions[0] += 1;
    variants.push(("RB dimensions", changed));
    let mut changed = base.clone();
    changed.concept = !changed.concept;
    variants.push(("concept flag", changed));
    let mut changed = base.clone();
    changed.mus[1] = f64::from_bits(changed.mus[1].to_bits() + 1);
    variants.push(("parameter battery", changed));
    let mut changed = base.clone();
    changed.tolerances[0] = f64::from_bits(changed.tolerances[0].to_bits() + 1);
    variants.push(("tolerance battery", changed));

    for (field, variant) in variants {
        assert_ne!(
            base_hash,
            variant.hash_hex(),
            "canonical coverage identity omitted {field}"
        );
    }
    assert!(base.canonical_text().contains("mus_bits=["));
    assert!(base.canonical_text().contains("tolerance_bits=["));
    verdict(
        "p3-000",
        "the domain-separated coverage-spec digest changes with every semantic model, ladder, \
         battery, algorithm, and schema input",
    );
}

#[test]
fn p3_001_proposal_a_numeric_floor_observed_but_activation_unmet() {
    // Proposal A's numeric kill-floor instrument is live, but the activation
    // gate requires certified RB coverage. This Estimated-only ladder can
    // observe the floor; it cannot satisfy that authority requirement.
    let measurement = measure_horizon_coverage();
    let coverage = measurement.coverage;
    println!(
        "{{\"metric\":\"horizon-A\",\"rb_coverage\":{coverage:.3},\"floor\":0.2,\
         \"numeric_floor_observed\":true,\"certification_activation_trigger_met\":false,\
         \"status\":\"NUMERIC KILL-FLOOR OBSERVED; CERTIFICATION/ACTIVATION TRIGGER UNMET\"}}"
    );
    assert!(
        coverage >= 0.2,
        "the numeric kill-floor is observed: {coverage}"
    );
    verdict(
        "p3-001",
        "Proposal A: the Estimated rb_coverage instrument observes the numeric 0.2 floor, \
         while certification and activation remain explicitly unmet",
    );
}

#[test]
fn p3_002_proposal_c_instrumented_awaiting_audit() {
    use fs_plan::voi::{AuditVerdict, audit_scheduling};
    // Proposal C's activation is CONDITIONAL: the machinery is live
    // (Phase-2 benchmarks), but SCHEDULING AUTHORITY requires the
    // prospective audit to show recommendations beat agent choices —
    // and with no audit evidence the verdict is Demote by design.
    let report =
        audit_scheduling("fs-flywheel-phase3-v1", &[]).expect("empty bounded audit reports safely");
    assert_eq!(
        report.verdict(),
        AuditVerdict::DemoteToReporting,
        "no evidence, no authority — the default is the safe one"
    );
    assert_eq!(report.verdict(), AuditVerdict::DemoteToReporting);
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
#[allow(clippy::too_many_lines)] // one auditable end-to-end holding-pen fixture
fn p3_005_proposal_11_r8_gate_and_the_holding_pen() -> Result<(), fs_asbuilt::RegError> {
    use fs_asbuilt::{Fiducial, Point2, register, well_posed};

    struct Phase3SignatureVerifier;

    impl fs_checker::SignatureVerifier for Phase3SignatureVerifier {
        fn verify(
            &self,
            request: &fs_checker::SignatureRequest<'_>,
        ) -> fs_checker::VerificationDecision {
            let fingerprint = fs_ledger::hash_bytes(b"fs-flywheel-e2e:phase3-signature-policy:v1");
            if request.signature
                == format!("phase3-horizon-gate:{}", request.subject_hash().to_hex())
                && request.purpose == fs_checker::SignaturePurpose::PackageRootAttestation
            {
                fs_checker::VerificationDecision::accept(fingerprint)
            } else {
                fs_checker::VerificationDecision::reject(fingerprint)
            }
        }
    }

    // Proposal 11's R8 instrument: registration must be tighter than
    // the deviations being certified. GOOD fiducials pass...
    let good = [(0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0)]
        .iter()
        .map(|&(x, y)| -> Result<Fiducial, fs_asbuilt::RegError> {
            Ok(Fiducial::new(
                Point2::new(x, y)?,
                Point2::new(x + 0.30001, y + 0.19999)?, // clean rigid shift
            ))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reg = with_cx(|cx| register(&good, cx))?;
    assert!(
        well_posed(&reg, 0.05),
        "clean registration certifies 0.05 deviations"
    );
    // ...and sloppy fiducials FAIL the same certification (the R8 kill
    // visible in the instrument).
    let sloppy = [(0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0)]
        .iter()
        .enumerate()
        .map(|(k, &(x, y))| -> Result<Fiducial, fs_asbuilt::RegError> {
            let wobble =
                0.2 * f64::from(u32::try_from(k).expect("four-point fixture index fits u32"));
            Ok(Fiducial::new(
                Point2::new(x, y)?,
                Point2::new(x + 0.3 + wobble, y + 0.2 - wobble)?,
            ))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reg = with_cx(|cx| register(&sloppy, cx))?;
    assert!(
        !well_posed(&reg, 0.05),
        "sloppy registration cannot certify what it cannot resolve"
    );
    // Recompute Proposal A's one canonical numeric kill-floor measurement for
    // the retained package. The ladder is Estimated-only, so certification and
    // activation remain unmet even though the numeric floor is observed.
    let measurement = measure_horizon_coverage();
    let coverage = measurement.coverage;
    assert!(coverage >= 0.2, "the measured numeric kill-floor holds");
    let estimator_id = measurement.spec.algorithm;
    let spec_canonical = measurement.spec.canonical_text();
    let spec_hash = measurement.spec.hash_hex();
    let coverage_statement = format!(
        "numeric kill-floor observed: rb_coverage={coverage:.17} (bits={:016x}) >= 0.2; \
         certification/activation trigger UNMET; spec_hash={spec_hash};spec={spec_canonical}",
        coverage.to_bits(),
    );

    // THE HOLDING PEN, IN WRITING: the five statuses enter a
    // fixture-authenticated package. This proves checker integration and honest
    // color separation, not external scientific certification.
    let unsigned = EvidencePackage::new(Provenance::new("phase3-horizon", "Cargo.lock"))
        .with_claim(Claim::estimated(
            "A-abstraction-ladder",
            coverage_statement,
            estimator_id,
            f64::INFINITY,
        ))
        .with_claim(Claim::estimated(
            "C-value-of-information",
            "instrumented; scheduling authority awaits the prospective audit",
            "prospective-audit-pending".to_string(),
            1.0,
        ))
        .with_claim(Claim::estimated(
            "4-spacetime-complex",
            "instrumented-but-uncontrolled; control gated on a splitting-dominated \
             paying workload",
            "workload-demand-pending".to_string(),
            1.0,
        ))
        .with_claim(Claim::estimated(
            "13b-symmetry-solver",
            "prevalence instrument live; dedicated solver waits for >=15% real-workload \
             symmetry",
            "prevalence-pending".to_string(),
            1.0,
        ))
        .with_claim(Claim::estimated(
            "11-reality-as-a-chart",
            "R8 instrument live (registration vs certified deviation); full-field \
             activation awaits metrology partnerships — point-sensor assimilation \
             ships meanwhile",
            "metrology-partnership-pending".to_string(),
            1.0,
        ));
    let root = unsigned.try_merkle_root().expect("bounded fixture root");
    let signature_subject = fs_checker::signature_subject_hash(
        root,
        fs_checker::SignaturePurpose::PackageRootAttestation,
    );
    let pkg = unsigned.signed(format!(
        "phase3-horizon-gate:{}",
        signature_subject.to_hex()
    ));
    let signature_verifier = Phase3SignatureVerifier;
    let capabilities =
        fs_checker::VerificationCapabilities::deny_all().with_signatures(&signature_verifier);
    let check = fs_checker::check_with_capabilities(&pkg, None, None, &capabilities);
    assert!(check.passed(), "the holding-pen record re-verifies");
    assert!(matches!(
        check.signature(),
        fs_checker::SignatureStatus::Authenticated(_)
    ));
    let breakdown = *check.breakdown();
    assert!(
        breakdown.verified == 0 && breakdown.estimated == 5,
        "all five horizon programs, including Proposal A's unmet certification trigger, remain \
         honestly Estimated: {breakdown:?}"
    );
    println!(
        "{{\"metric\":\"horizon-ledger\",\"root\":\"{}\",\"numeric_floor_observed\":1,\
         \"activation_trigger_met\":0,\"verified\":0,\"waiting\":5,\
         \"coverage_spec_hash\":\"{spec_hash}\"}}",
        pkg.try_merkle_root().expect("bounded fixture root")
    );
    verdict(
        "p3-005",
        "Proposal 11's R8 instrument passes clean registration and fails sloppy; the \
         five-proposal holding pen retains Proposal A's numeric floor with its certification \
         trigger unmet, alongside four other waiting programs, as Estimated claims under an \
         authenticated package-root signature",
    );
    Ok(())
}
