//! The as-built product chain end to end (bead
//! `frankensim-extreal-program-f85xj.12.2`, e2e DONE-WHEN):
//! register → bind → propagate → budget → render, with per-stage logging that
//! includes the correlation structure of the emitted geometry term.
//!
//! Each link is owned by a different crate and was previously testable only in
//! isolation:
//!
//! | stage     | crate         | what it produces                          |
//! |-----------|---------------|-------------------------------------------|
//! | register  | `fs-asbuilt`  | `CalibratedRigid3Registration` + identity  |
//! | bind      | `fs-scenario` | `PlacementBasis::AsBuilt { .. }` citation  |
//! | propagate | `fs-asbuilt`  | `GeometryPropagation` + per-QoI terms      |
//! | budget    | `fs-evidence` | the eight-term budget, geometry populated  |
//! | render    | `fs-report`   | the nominal-versus-as-built projection     |
//!
//! What this case proves is that ONE identity survives all five stages: the
//! registration's `model_identity()` is what the scenario cites, and the
//! propagation record's identity is what every geometry term and the rendered
//! correlation line agree on.
//!
//! NO-CLAIM: the "solve" here is a linear QoI evaluator, not a physics solve.
//! This case proves the seams compose and the identities survive; it does not
//! establish that any thermal result is correct, nor that a real solver
//! consumes the placement basis yet.

#![allow(clippy::float_cmp)] // Identity assertions compare exact values on purpose.

use fs_asbuilt::propagate::{
    CoveragePolicy, GeometryPropagation, QoiEvaluator, QoiSensitivity, propagate_pose_covariance,
};
use fs_asbuilt::rigid3::{
    CalibratedRigid3Registration, Covariance3, CrossFiducialModel3, Fiducial3, MetrologyModel3,
    Point3, estimate_calibrated_rigid3,
};
use fs_asbuilt::uncertainty::HuberPolicy;
use fs_blake3::ContentHash;
use fs_evidence::uncertainty::{
    EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
    TermValue, UncertaintyArtifactRef,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_report::{AsBuiltQoiDelta, nominal_vs_as_built_markdown};
use fs_scenario::FrameId;
use fs_scenario::entity::{
    EntityCatalog, EntityDeclaration, EntityId, GeometryFingerprint, PlacementBasis,
};
use std::fmt::Write as _;

/// Fixed measurement noise, in metres. Deterministic literals rather than a
/// seeded RNG: this case is about identity flow, so the scan must be byte-stable
/// without depending on any generator's own contract.
const SIGMA: f64 = 1.0e-3;

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).expect("finite fixture point")
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = fs_exec::VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x12_2E_2E_01,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(
        pool.stats().quiescent(),
        "Cx arena must be quiescent after scope: {}",
        pool.stats().to_json()
    );
    result
}

/// A stand-in for the solve: QoI deltas that are linear in the pose delta.
struct LinearEvaluator {
    gradients: Vec<[f64; 6]>,
}

impl QoiEvaluator for LinearEvaluator {
    fn evaluate(&self, pose_delta: &[f64; 6]) -> Result<Vec<f64>, String> {
        Ok(self
            .gradients
            .iter()
            .map(|gradient| (0..6).map(|axis| gradient[axis] * pose_delta[axis]).sum())
            .collect())
    }
}

/// Stage 1 — register a measured scan against the design fiducials.
fn stage_register(cx: &Cx<'_>) -> CalibratedRigid3Registration {
    let design = [
        p3(0.0, 0.0, 0.0),
        p3(2.0, 0.0, 0.0),
        p3(2.0, 3.0, 0.0),
        p3(0.0, 3.0, 0.0),
        p3(0.0, 0.0, 1.0),
        p3(2.0, 0.0, 1.0),
        p3(2.0, 3.0, 1.0),
        p3(0.5, 1.0, 0.6),
    ];
    // Deterministic per-axis perturbations, one SIGMA-scaled literal per point.
    let noise: [[f64; 3]; 8] = [
        [0.9, -0.4, 0.2],
        [-0.7, 0.3, -0.5],
        [0.1, 0.8, 0.6],
        [-0.3, -0.9, 0.4],
        [0.5, 0.2, -0.8],
        [-0.6, 0.7, 0.1],
        [0.4, -0.2, -0.3],
        [-0.1, 0.6, 0.9],
    ];
    let fiducials: Vec<Fiducial3> = design
        .iter()
        .zip(noise.iter())
        .map(|(point, delta)| {
            Fiducial3::new(
                *point,
                p3(
                    point.x() + SIGMA * delta[0],
                    point.y() + SIGMA * delta[1],
                    point.z() + SIGMA * delta[2],
                ),
            )
        })
        .collect();
    let covariance = Covariance3::new(SIGMA * SIGMA, 0.0, 0.0, SIGMA * SIGMA, 0.0, SIGMA * SIGMA)
        .expect("diagonal fiducial covariance");
    let model = MetrologyModel3::new(
        vec![covariance; fiducials.len()],
        CrossFiducialModel3::Independent,
        HuberPolicy::Disabled,
        "cmm-calibration-2026-07/cold-plate-stack@rev3",
    )
    .expect("metrology model");
    estimate_calibrated_rigid3(&fiducials, &model, cx).expect("calibrated registration")
}

/// Stage 2 — bind the registration into a scenario as an as-built placement.
///
/// Returns the catalog and the placed occurrence.
fn stage_bind(registration_ref: ContentHash) -> (EntityCatalog, EntityId) {
    let mut catalog = EntityCatalog::new();
    let assembly = catalog
        .declare(
            EntityDeclaration::assembly("assembly/cold-plate-stack")
                .with_fingerprint(GeometryFingerprint::of_bytes(b"assembly/cold-plate-stack")),
        )
        .expect("assembly");
    let plate = catalog
        .declare(
            EntityDeclaration::part(assembly, "part/cold-plate")
                .with_fingerprint(GeometryFingerprint::of_bytes(b"part/cold-plate")),
        )
        .expect("plate");
    catalog
        .declare_placement(
            plate,
            FrameId(0),
            PlacementBasis::AsBuilt { registration_ref },
        )
        .expect("as-built placement");
    (catalog, plate)
}

/// The three QoIs this fixture propagates, with their pose gradients.
fn sensitivities() -> (Vec<QoiSensitivity>, Vec<[f64; 6]>) {
    let gradients = vec![
        [0.0, 0.0, 1.0, 0.2, -0.1, 0.0],
        [0.0, 0.0, 1.0, -0.2, 0.1, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.8],
    ];
    let qois = vec![
        QoiSensitivity::new("t-junction-gap", "millimetre", gradients[0]).expect("qoi 0"),
        QoiSensitivity::new("contact-plane-gap", "millimetre", gradients[1]).expect("qoi 1"),
        QoiSensitivity::new("connector-offset-x", "millimetre", gradients[2]).expect("qoi 2"),
    ];
    (qois, gradients)
}

/// Stage 4 — assemble one complete eight-term budget for QoI `ordinal`.
fn stage_budget(
    propagation: &GeometryPropagation,
    ordinal: usize,
    qoi: &str,
) -> EngineeringUncertaintyBudget {
    let artifact = |role: &str| {
        UncertaintyArtifactRef::new(role, propagation.record_identity()).expect("artifact ref")
    };
    let terms: Vec<EngineeringUncertaintyTerm> = EngineeringUncertaintyKind::ALL
        .into_iter()
        .map(|kind| match kind {
            EngineeringUncertaintyKind::Geometry => propagation
                .geometry_term(ordinal)
                .expect("geometry term from the shared propagation record"),
            EngineeringUncertaintyKind::Roundoff
            | EngineeringUncertaintyKind::SolverAlgebraic
            | EngineeringUncertaintyKind::Discretization => EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::interval(0.0, 1e-9).expect("numerical interval"),
                artifact("numerical-certificate-placeholder"),
            )
            .expect("numerical term"),
            EngineeringUncertaintyKind::Parameters
            | EngineeringUncertaintyKind::BoundaryConditions => {
                EngineeringUncertaintyTerm::try_new(
                    kind,
                    TermValue::unknown("not propagated in this as-built chain fixture")
                        .expect("named unknown"),
                    artifact("declared-gap"),
                )
                .expect("unknown term")
            }
            EngineeringUncertaintyKind::ModelForm | EngineeringUncertaintyKind::Measurement => {
                EngineeringUncertaintyTerm::try_new(
                    kind,
                    TermValue::negligible("synthetic fixture with declared noise only")
                        .expect("named negligible"),
                    artifact("fixture-declaration"),
                )
                .expect("negligible term")
            }
            // The kind enum is non-exhaustive upstream; a source this fixture
            // does not understand is an honest evidence gap, never a silent
            // negligible.
            _ => EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::unknown("source kind unknown to this fixture").expect("named unknown"),
                artifact("declared-gap"),
            )
            .expect("wildcard term"),
        })
        .collect();
    EngineeringUncertaintyBudget::try_new(qoi, "millimetre", terms).expect("eight-term budget")
}

#[test]
fn e2e_register_bind_propagate_budget_render_keeps_one_identity() {
    let mut log = String::new();
    let (qois, gradients) = sensitivities();
    let evaluator = LinearEvaluator {
        gradients: gradients.clone(),
    };

    // ---- stage 1: register -------------------------------------------------
    let (registration, propagation) = with_default_cx(|cx| {
        let registration = stage_register(cx);
        let propagation = propagate_pose_covariance(
            &registration,
            &qois,
            CoveragePolicy::new(0.95, 2.0).expect("coverage policy"),
            Some(&evaluator),
            1e-6,
            "",
            cx,
        )
        .expect("propagation");
        (registration, propagation)
    });
    let registration_ref = registration.model_identity();
    let _ = writeln!(
        log,
        "{{\"stage\":\"register\",\"dof\":{},\"identity\":\"{registration_ref}\"}}",
        registration.degrees_of_freedom()
    );

    // ---- stage 2: bind -----------------------------------------------------
    let (catalog, plate) = stage_bind(registration_ref);
    let citations = catalog.as_built_registrations();
    assert_eq!(
        citations,
        vec![(plate, registration_ref)],
        "the scenario must cite exactly the registration it was bound to"
    );
    let _ = writeln!(
        log,
        "{{\"stage\":\"bind\",\"occurrence\":\"{plate}\",\"basis\":\"{}\",\"citation\":\"{registration_ref}\"}}",
        catalog
            .placement_of(plate)
            .expect("placement")
            .basis()
            .label()
    );

    // ---- stage 3: resolve --------------------------------------------------
    // The product layer's obligation: the cited identity must be the record it
    // is about to propagate. A wrong citation is detectable here and nowhere
    // earlier, because fs-scenario deliberately never resolves it.
    let (_, cited) = citations[0];
    assert_eq!(cited, registration.model_identity());
    let impostor = ContentHash([0x5A; 32]);
    assert_ne!(
        impostor,
        registration.model_identity(),
        "a foreign citation must not resolve to this registration"
    );
    let _ = writeln!(
        log,
        "{{\"stage\":\"resolve\",\"matches\":true,\"record\":\"{}\"}}",
        propagation.record_identity()
    );

    // ---- stage 4: propagate + budget --------------------------------------
    // One pose covariance moved everything, so the two gap QoIs that share the
    // pose translation must be genuinely correlated, not independently noisy.
    let rho = propagation.correlation(0, 1).expect("rho01");
    assert!(
        rho > 0.5,
        "shared translation must correlate the two gaps: {rho}"
    );
    let budgets: Vec<EngineeringUncertaintyBudget> = qois
        .iter()
        .enumerate()
        .map(|(ordinal, qoi)| stage_budget(&propagation, ordinal, qoi.name()))
        .collect();
    for budget in &budgets {
        let geometry = budget.term(EngineeringUncertaintyKind::Geometry);
        assert!(matches!(geometry.value(), TermValue::Distribution(_)));
        // Every geometry term cites the SAME propagation record. That shared
        // citation IS the correlation structure travelling downstream.
        assert_eq!(
            geometry.provenance().digest(),
            propagation.record_identity(),
            "each geometry term must cite the one shared propagation record"
        );
    }
    let _ = writeln!(
        log,
        "{{\"stage\":\"budget\",\"record\":\"{}\",\"method\":\"{:?}\",\"sd\":{:?},\"rho01\":{rho:.6}}}",
        propagation.record_identity(),
        propagation.method(),
        propagation.standard_deviations()
    );

    // ---- stage 5: render ---------------------------------------------------
    // Nominal values are the as-designed solve; the as-built solve is displaced
    // by the measured pose. Exactly representable literals keep the projection
    // assertions about the projection, not about decimal conversion.
    let nominal = [0.500_f64, 0.750, 1.250];
    let as_built = [0.625_f64, 0.875, 1.250];
    let geometry_terms: Vec<&EngineeringUncertaintyTerm> = budgets
        .iter()
        .map(|budget| budget.term(EngineeringUncertaintyKind::Geometry))
        .collect();
    let rows: Vec<AsBuiltQoiDelta<'_>> = qois
        .iter()
        .enumerate()
        .map(|(ordinal, qoi)| {
            AsBuiltQoiDelta::try_new(
                qoi.name(),
                qoi.unit(),
                nominal[ordinal],
                as_built[ordinal],
                geometry_terms[ordinal],
            )
            .expect("comparison row")
        })
        .collect();
    let rendered = nominal_vs_as_built_markdown(&rows);

    assert!(rendered.contains("t-junction-gap"));
    assert!(rendered.contains("shift `+0.125 millimetre`"));
    // The third QoI did not move at all; the table must still render it rather
    // than silently dropping an unchanged quantity.
    assert!(rendered.contains("shift `+0 millimetre`"));
    // The render must reach the same correlated-block conclusion the budgets
    // encode, and must name the exact record.
    assert!(
        rendered.contains("all 3 geometry terms cite one propagation record"),
        "render lost the shared-record correlation\n{rendered}"
    );
    assert!(rendered.contains(&propagation.record_identity().to_hex()));
    assert!(rendered.contains("**No-claim boundary:**"));
    // Determinism: the whole chain is replayable to the same bytes.
    assert_eq!(rendered, nominal_vs_as_built_markdown(&rows));
    let _ = writeln!(
        log,
        "{{\"stage\":\"render\",\"rows\":{},\"bytes\":{}}}",
        rows.len(),
        rendered.len()
    );

    // ---- per-stage forensic log -------------------------------------------
    for marker in [
        "\"stage\":\"register\"",
        "\"stage\":\"bind\"",
        "\"stage\":\"resolve\"",
        "\"stage\":\"budget\"",
        "\"stage\":\"render\"",
        "\"rho01\":",
    ] {
        assert!(log.contains(marker), "forensic log lost {marker:?}\n{log}");
    }
    println!("{log}{rendered}");
}

#[test]
fn e2e_chain_is_deterministic_across_independent_runs() {
    // The same fixture run twice must produce the same registration identity,
    // the same propagation record, and therefore the same citation. If any
    // stage picked up ambient state, these identities would drift.
    let (qois, gradients) = sensitivities();
    let evaluator = LinearEvaluator { gradients };
    let run = || {
        with_default_cx(|cx| {
            let registration = stage_register(cx);
            let propagation = propagate_pose_covariance(
                &registration,
                &qois,
                CoveragePolicy::new(0.95, 2.0).expect("coverage policy"),
                Some(&evaluator),
                1e-6,
                "",
                cx,
            )
            .expect("propagation");
            (registration.model_identity(), propagation.record_identity())
        })
    };
    let (first_registration, first_record) = run();
    let (second_registration, second_record) = run();
    assert_eq!(first_registration, second_registration);
    assert_eq!(first_record, second_record);

    // And the binding transports that identity unchanged.
    let (catalog, plate) = stage_bind(first_registration);
    assert_eq!(
        catalog.as_built_registrations(),
        vec![(plate, second_registration)]
    );
    println!(
        "{{\"suite\":\"fs-diffreal-e2e/as-built-chain\",\"case\":\"determinism\",\"registration\":\"{first_registration}\",\"record\":\"{first_record}\"}}"
    );
}
