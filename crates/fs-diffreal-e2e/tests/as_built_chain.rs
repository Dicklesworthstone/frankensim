//! The as-built product chain end to end (beads
//! `frankensim-extreal-program-f85xj.12.2` and `.12.10`, e2e DONE-WHEN):
//! register → bind → propagate → budget → render, plus a physics stage that
//! drives a real conduction solve from a measured thickness field, with
//! per-stage logging that includes the correlation structure of the emitted
//! geometry term.
//!
//! Each link is owned by a different crate and was previously testable only in
//! isolation:
//!
//! | stage     | crate            | what it produces                        |
//! |-----------|------------------|-----------------------------------------|
//! | register  | `fs-asbuilt`     | `CalibratedRigid3Registration` + identity|
//! | bind      | `fs-scenario`    | `PlacementBasis::AsBuilt { .. }` citation|
//! | propagate | `fs-asbuilt`     | `GeometryPropagation` + per-QoI terms    |
//! | budget    | `fs-evidence`    | the eight-term budget, geometry populated|
//! | render    | `fs-report`      | the nominal-versus-as-built projection   |
//! | solve     | `fs-conduction`  | a solve whose per-face `R''` IS the scan |
//!
//! What these cases prove is that ONE identity survives every stage: the
//! registration's `model_identity()` is what the scenario cites, what the
//! propagation record is about, and what the thickness field driving the
//! solve reports as its `registration_ref()`.
//!
//! # No-claim boundaries
//!
//! - The QoI evaluator in the propagate/budget/render stages is LINEAR, not a
//!   physics solve. Those stages prove the seams compose and the identities
//!   survive; they do not establish that any QoI value is correct.
//! - The solve stage is a real `fs-conduction` solve, but on SYNTHETIC
//!   fixtures throughout: literal-noise fiducials, a declared interface card,
//!   an injected wedge defect, and the crate's G1 two-slab geometry (metre
//!   scale, not an electronics package). It demonstrates that the seam
//!   carries a measured field into the operator and that the result moves in
//!   the physically required direction. It validates no magnitude.
//! - Every thickness value is Estimate-class and inherits `fs-asbuilt`'s
//!   no-claims unchanged: supplied correspondence (never discovered),
//!   unfiltered roughness, a first-order pose linearisation, and a composed
//!   half-width that is a declared decomposition rather than a confidence
//!   interval. `fs-conduction` does not validate a caller-supplied `R''`
//!   either — a measured map entering a solve does not make the solve
//!   validated.
//! - `f85xj.12.3`'s DONE-WHEN 2 (a real or rig-derived scan) stays open and
//!   blocked on Level-E hardware (`f85xj.4.5`). Nothing here substitutes for
//!   it.

// Identity assertions compare exact values on purpose.
#![allow(clippy::float_cmp)]
// A chain case is long because the chain is long. Splitting a stage out into a
// helper would hide exactly the thing under test — that one identity survives
// every hop — behind a function boundary, so the length is the point.
#![allow(clippy::too_many_lines)]

use fs_asbuilt::field::{ProbeModel, SurfaceSample, ThicknessField};
use fs_asbuilt::propagate::{
    CoveragePolicy, GeometryPropagation, QoiEvaluator, QoiSensitivity, propagate_pose_covariance,
};
use fs_asbuilt::rigid3::{
    CalibratedRigid3Registration, Covariance3, CrossFiducialModel3, Fiducial3, MetrologyModel3,
    Point3, estimate_calibrated_rigid3,
};
use fs_asbuilt::uncertainty::HuberPolicy;
use fs_blake3::ContentHash;
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule,
    solve_with_interfaces,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    InterfaceFacePair, InterfaceResistance, InterfaceSurface, ResistanceValueOrigin,
    ThermalInterfaces,
};
use fs_evidence::ValidityDomain;
use fs_evidence::uncertainty::{
    EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
    TermValue, UncertaintyArtifactRef,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
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

// ===========================================================================
// Stage 6 — solve (bead `frankensim-extreal-program-f85xj.12.10`)
//
// `f85xj.12.3` built the L2 as-built field extraction and `f85xj.12.9` built
// the L3 per-face consumption seam. Both were green in their own crates and
// NOTHING composed them, so "as-built fields can drive a conduction solve"
// was an inference from two passing suites rather than an observed result.
// Everything below turns it into one.
//
// The physical story: two slabs are bonded across an adhesive layer. Uneven
// clamping leaves the bond line WEDGE-SHAPED — thicker toward +y. A uniform
// `R''` cannot represent that; a per-face map extracted from the scan can.
// ===========================================================================

/// Conductivity of both slabs, W/(m·K).
const K_SOLID_W_PER_MK: f64 = 10.0;
/// The cold wall, K.
const T_COLD_K: f64 = 300.0;
/// Heat driven IN through the hot wall, W/m². `ThermalBc::neumann` takes the
/// OUTWARD flux, so this is negated at the boundary row.
///
/// Driving with flux and pinning only the cold side is what makes the
/// hot-side temperature respond to `R''` at all: with both walls pinned,
/// changing the bond line moves the heat rate instead.
const HOT_FLUX_W_PER_M2: f64 = 100.0;
/// The card's area-specific contact resistance at nominal thickness, m²·K/W.
const CARD_RESISTANCE_M2K_PER_W: f64 = 0.1;
/// The nominal bond-line thickness the card's claim describes, m.
const NOMINAL_BOND_THICKNESS_M: f64 = 2.5e-4;
/// Bond-line thickening per unit `+y` at wedge scale 1, m. Split unevenly
/// across the two faces so a recovered thickness is a genuine SUM over both,
/// not one face's growth doubled.
const BOND_WEDGE_M: f64 = 1.0e-4;
/// Coverage level and factor for both the propagation and the field.
const COVERAGE_LEVEL: f64 = 0.95;
const COVERAGE_FACTOR: f64 = 2.0;
/// How far a wedge map's conductance must exceed the conductance of a uniform
/// surface at the same MEAN thickness. See the assertion for why a bare `>`
/// does not discriminate.
const MAP_VS_AVERAGE_MARGIN: f64 = 1.005;

/// Outward normal of the bond line's face against the HOT slab.
///
/// `fs-asbuilt`'s convention is "outward means away from the solid", and the
/// solid being measured here is the ADHESIVE, not either slab. So outward at
/// the hot-side face points back into the hot slab (`-x`) and outward at the
/// cold-side face points into the cold slab (`+x`). Growth along those
/// normals is more adhesive, i.e. a thicker bond line — which is what must
/// raise the hot-side temperature.
const BOND_FACE_NORMALS: [[f64; 3]; 2] = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

fn mat3_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn coverage() -> CoveragePolicy {
    CoveragePolicy::new(COVERAGE_LEVEL, COVERAGE_FACTOR).expect("coverage policy")
}

/// Two unit slabs meeting at `x = 1`, exactly duplicated on the seam so the
/// matching-P1 contact operator admits them.
fn two_slab_mesh() -> ConductionMesh {
    let (left, mut positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let left_vertex_count = positions.len();
    let (right, right_positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let offset = u32::try_from(left_vertex_count).expect("fixture vertex count fits u32");
    let mut tets = left.tets;
    tets.extend(
        right
            .tets
            .into_iter()
            .map(|tet| tet.map(|vertex| vertex + offset)),
    );
    positions.extend(right_positions.into_iter().map(|[x, y, z]| [x + 1.0, y, z]));
    let complex = TetComplex::from_tets(positions.len(), tets);
    ConductionMesh::new(complex, positions).expect("two-slab mesh")
}

/// Flux in at `x = 0`, a pinned wall at `x = 2`, and an explicit adiabatic
/// remainder — the interface traces MUST land there, because any boundary row
/// on an interface slot is a typed refusal.
fn thermal_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "hot",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::neumann(-HOT_FLUX_W_PER_M2).expect("hot flux row"),
        )
        .expect("hot boundary")
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], 2.0),
            ThermalBc::dirichlet(T_COLD_K).expect("cold wall"),
        )
        .expect("cold boundary")
        .adiabatic_remainder()
        .finish()
        .expect("complete boundary partition")
}

fn bond_surface(material: &str, texture: &str) -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: material.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: texture.to_string(),
    }
}

/// The one material card every face of the map must cite.
fn bond_card() -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: CARD_RESISTANCE_M2K_PER_W,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic as-built bond-line fixture (f85xj.12.10)".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("bond-line resistance claim inserts");
    InterfaceSystemCard::assemble(
        bond_surface("slab-hot", "bond-normal-plus-x"),
        bond_surface("slab-cold", "bond-normal-minus-x"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-adhesive-layer".to_string()),
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .expect("bond-line card assembles")
}

/// Side A is the HOT slab, so a positive `heat_rate_a_to_b_w` is heat leaving
/// the hot side.
fn oriented_pairs(mesh: &ConductionMesh) -> Vec<InterfaceFacePair> {
    ThermalInterfaces::coincident_face_pairs(mesh)
        .expect("coincident bond-line pairs")
        .into_iter()
        .map(|pair| {
            if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                pair
            } else {
                InterfaceFacePair {
                    side_a: pair.side_b,
                    side_b: pair.side_a,
                }
            }
        })
        .collect()
}

/// Injected per-face growth at a station, m: the wedge defect.
fn injected_growth(centroid: [f64; 3], scale: f64) -> [f64; 2] {
    let wedge = scale * BOND_WEDGE_M * centroid[1];
    [wedge, 0.5 * wedge]
}

/// Everything a rung of the wedge ladder shares.
struct BondFixture<'m> {
    mesh: &'m ConductionMesh,
    boundary: &'m ThermalBoundary,
    material: &'m ConductivityModel,
    source: &'m ScalarField,
    card: &'m InterfaceSystemCard,
    pairs: &'m [InterfaceFacePair],
    registration: &'m CalibratedRigid3Registration,
}

/// What one wedge scale produced.
struct BondRung {
    scale: f64,
    field_identity: ContentHash,
    field_registration_ref: ContentHash,
    thickness: Vec<f64>,
    area_m2: f64,
    max_temperature_k: f64,
    conductance_w_per_k: f64,
    mean_jump_k: f64,
    heat_rate_a_to_b_w: f64,
    relative_closure: f64,
    flux_card_identity: ContentHash,
}

impl BondRung {
    fn min_thickness_m(&self) -> f64 {
        self.thickness.iter().copied().fold(f64::INFINITY, f64::min)
    }

    fn max_thickness_m(&self) -> f64 {
        self.thickness
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn mean_thickness_m(&self) -> f64 {
        self.thickness.iter().sum::<f64>() / self.thickness.len() as f64
    }

    /// The conductance a UNIFORM surface at this rung's mean thickness would
    /// have had, over the same total area.
    fn averaged_conductance_w_per_k(&self, k_bond_w_per_mk: f64) -> f64 {
        self.area_m2 / (self.mean_thickness_m() / k_bond_w_per_mk)
    }
}

/// Paired opposing-face stations over the bond line, one per interface face
/// pair, in the SAME order as `pairs`.
///
/// The measured point is built through the fitted pose — `R·d + t + g·(R·n)`
/// — so the injected growth `g` is recovered exactly whatever pose the
/// fiducials happened to produce. That is the point of the paired path: both
/// faces ride one pose, so the pose cancels out of their difference.
fn bond_stations(
    fixture: &BondFixture<'_>,
    scale: f64,
) -> (Vec<SurfaceSample>, Vec<SurfaceSample>, Vec<f64>) {
    let rigid = fixture.registration.registration();
    let rotation = rigid.rotation();
    let half = 0.5 * NOMINAL_BOND_THICKNESS_M;
    let mut side_a = Vec::with_capacity(fixture.pairs.len());
    let mut side_b = Vec::with_capacity(fixture.pairs.len());
    let mut nominal = Vec::with_capacity(fixture.pairs.len());
    for pair in fixture.pairs {
        let centroid = fixture.mesh.boundary()[pair.side_a].centroid;
        let growth = injected_growth(centroid, scale);
        // The adhesive occupies the gap straddling the mesh seam: its
        // hot-side face sits half a nominal thickness before the seam and its
        // cold-side face half a thickness after it.
        for (slot, along) in [-half, half].into_iter().enumerate() {
            let design = p3(centroid[0] + along, centroid[1], centroid[2]);
            let registered = rigid.apply(design).expect("registered design point");
            let normal = BOND_FACE_NORMALS[slot];
            let rotated = mat3_vec(rotation, normal);
            let measured = p3(
                registered.x() + growth[slot] * rotated[0],
                registered.y() + growth[slot] * rotated[1],
                registered.z() + growth[slot] * rotated[2],
            );
            let sample = SurfaceSample::new(design, normal, measured).expect("bond-line sample");
            if slot == 0 {
                side_a.push(sample);
            } else {
                side_b.push(sample);
            }
        }
        nominal.push(NOMINAL_BOND_THICKNESS_M);
    }
    (side_a, side_b, nominal)
}

/// Convert a measured thickness map into per-face area-specific resistances.
///
/// The card's claim is the resistance of the NOMINAL bond line, so it implies
/// the adhesive's conductivity at that thickness: `k = t_nominal / R''_card`.
/// A measured station then converts as `R''_i = t_i / k`. The card supplies
/// `k` — the MATERIAL authority — and metrology supplies `t_i`;
/// `with_measured_value` is what records that split, so a later reader can
/// see the number is not the card's own claim.
fn measured_resistances(
    base: &InterfaceResistance,
    field: &ThicknessField,
) -> Vec<InterfaceResistance> {
    let k_bond = NOMINAL_BOND_THICKNESS_M / base.value_m2_k_per_w();
    field
        .points()
        .iter()
        .map(|point| {
            base.with_measured_value(
                point.thickness() / k_bond,
                // The supplied uncertainty must describe the SUPPLIED value.
                // The field's half-width is a length; the same `k` converts
                // it. This is the caller's composition obligation, and it is
                // still the field's declared decomposition, not a confidence
                // interval.
                UncertaintyModel::HalfWidth {
                    half_width: point.half_width() / k_bond,
                    confidence: COVERAGE_LEVEL,
                },
                "as-built bond-line thickness map from a registered scan",
            )
            .expect("measured face resistance")
        })
        .collect()
}

/// One rung: extract the field at this wedge scale, map it onto the faces,
/// and solve.
fn solve_rung(fixture: &BondFixture<'_>, scale: f64) -> BondRung {
    let (side_a, side_b, nominal) = bond_stations(fixture, scale);
    let base = InterfaceResistance::from_card(
        "bondline",
        fixture.card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("card-backed bond-line resistance");
    // The probe looks along the bond normal, so neither face is grazing and
    // the station order below is the supplied order.
    let probe = ProbeModel::new([1.0, 0.0, 0.0], 0.05, 5.0e-6, 5.0e-6).expect("bond-line probe");

    with_default_cx(|cx| {
        let field = ThicknessField::extract(
            &side_a,
            &side_b,
            &nominal,
            fixture.registration,
            &probe,
            coverage(),
            cx,
        )
        .expect("thickness field extracts");
        // Load-bearing: `points()` is COMPACTED past refused stations, so the
        // station-to-face-pair index mapping is only sound when nothing was
        // refused. Assert it rather than assuming it.
        assert!(
            field.grazing().is_empty(),
            "no bond-line station may be refused: {:?}",
            field.grazing()
        );
        assert_eq!(
            field.points().len(),
            fixture.pairs.len(),
            "one thickness station per interface face pair"
        );

        let resistances = measured_resistances(&base, &field);
        // Every per-face value is caller-supplied; none may masquerade as the
        // card's own claim. Counted over the whole set rather than reported
        // one face at a time: "3 of 8 laundered" and "1 of 8 laundered" are
        // different defects and a first-failure assert cannot tell them apart.
        let laundered = resistances
            .iter()
            .filter(|resistance| {
                matches!(resistance.value_origin(), ResistanceValueOrigin::CardClaim)
            })
            .count();
        assert_eq!(
            laundered,
            0,
            "{laundered} of {} faces laundered a measurement through the card",
            resistances.len()
        );
        for (index, resistance) in resistances.iter().enumerate() {
            if let ResistanceValueOrigin::CallerSupplied { rationale } = resistance.value_origin() {
                assert!(
                    rationale.contains("as-built"),
                    "face {index} lost its measured rationale"
                );
            }
            assert_eq!(
                resistance.card_identity(),
                base.card_identity(),
                "a map varies the VALUE, not the authority"
            );
        }

        let surface =
            InterfaceSurface::new_mapped("bondline", fixture.pairs.to_vec(), resistances.clone())
                .expect("mapped bond-line surface");
        let interfaces = ThermalInterfaces::new(fixture.mesh, fixture.boundary, vec![surface])
            .expect("complete interface binding");
        assert_eq!(interfaces.surface_is_mapped("bondline"), Some(true));

        let solution = solve_with_interfaces(
            cx,
            ConductionProblem {
                mesh: fixture.mesh,
                boundary: fixture.boundary,
                material: fixture.material,
                source: fixture.source,
            },
            &interfaces,
            SolveConfig {
                nonlinearity: Nonlinearity::FixedPoint {
                    relaxation: 1.0,
                    max_backtracks: 8,
                },
                stop: StopRule {
                    residual_rtol: 1e-11,
                    residual_atol: 1e-24,
                    step_atol: 0.0,
                    max_iterations: 12,
                },
                linear: LinearConfig {
                    tolerance: 1e-13,
                    max_iterations: 60_000,
                    restart: 60,
                },
                initial: InitialGuess::DirichletMean,
            },
        )
        .expect("bond-line conduction solve");

        let flux = solution
            .report
            .interface_fluxes
            .first()
            .expect("one interface flux");
        BondRung {
            scale,
            field_identity: field.identity(),
            field_registration_ref: field.registration_ref(),
            thickness: field
                .points()
                .iter()
                .map(fs_asbuilt::field::ThicknessPoint::thickness)
                .collect(),
            area_m2: flux.area_m2,
            max_temperature_k: solution
                .temperature
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
            conductance_w_per_k: flux.conductance_w_per_k,
            mean_jump_k: flux.mean_jump_k,
            heat_rate_a_to_b_w: flux.heat_rate_a_to_b_w,
            relative_closure: solution.report.energy.relative_closure(),
            flux_card_identity: flux.card_identity,
        }
    })
}

#[test]
fn e2e_as_built_thickness_map_drives_a_conduction_solve() {
    let mut log = String::new();
    let mesh = two_slab_mesh();
    let boundary = thermal_boundary(&mesh);
    let material = ConductivityModel::isotropic_declared(K_SOLID_W_PER_MK).expect("slab material");
    let source = ScalarField::Uniform(0.0);
    let card = bond_card();
    let card_identity = card.content_hash();
    let pairs = oriented_pairs(&mesh);
    assert!(
        pairs.len() > 1,
        "a bond line with one face cannot exercise a MAP: {}",
        pairs.len()
    );

    // ---- the identity the earlier stages bound ----------------------------
    let registration = with_default_cx(stage_register);
    let registration_ref = registration.model_identity();
    let (catalog, plate) = stage_bind(registration_ref);
    assert_eq!(
        catalog.as_built_registrations(),
        vec![(plate, registration_ref)],
        "the scenario must cite exactly the registration the solve will consume"
    );

    let fixture = BondFixture {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        card: &card,
        pairs: &pairs,
        registration: &registration,
    };

    // ---- stage 6: solve, on a ladder of wedge severities -------------------
    let ladder: Vec<BondRung> = [0.0, 1.0, 2.0]
        .into_iter()
        .map(|scale| solve_rung(&fixture, scale))
        .collect();

    for rung in &ladder {
        // (c) identity/provenance: the field the solve consumed is about the
        // SAME registration stage 2 bound into the scenario. This is the link
        // that made this bead worth writing: without it, "the solve used the
        // scan" is an inference from two green crates.
        assert_eq!(
            rung.field_registration_ref, registration_ref,
            "the thickness field driving the solve must cite the bound registration"
        );
        // The interface flux still cites the one material card.
        assert_eq!(
            rung.flux_card_identity, card_identity,
            "a per-face map must not fragment the material authority"
        );
        // Conservation, independent of the bond line: in steady state every
        // watt driven in at the hot wall crosses the interface.
        assert!(
            (rung.heat_rate_a_to_b_w - HOT_FLUX_W_PER_M2).abs() < 1e-6,
            "interface heat rate {} W must equal the driven {HOT_FLUX_W_PER_M2} W at scale {}",
            rung.heat_rate_a_to_b_w,
            rung.scale
        );
        assert!(
            rung.relative_closure < 1e-9,
            "energy balance must close at scale {}: {:e}",
            rung.scale,
            rung.relative_closure
        );
        let _ = writeln!(
            log,
            "{{\"stage\":\"solve\",\"wedge_scale\":{:.1},\"stations\":{},\
             \"t_min_m\":{:e},\"t_max_m\":{:e},\"conductance_w_per_k\":{:.9},\
             \"interface_jump_k\":{:.9},\"t_hot_k\":{:.9},\"field\":\"{}\",\
             \"registration\":\"{}\"}}",
            rung.scale,
            rung.thickness.len(),
            rung.min_thickness_m(),
            rung.max_thickness_m(),
            rung.conductance_w_per_k,
            rung.mean_jump_k,
            rung.max_temperature_k,
            rung.field_identity,
            rung.field_registration_ref
        );
    }

    // ---- the zero-wedge rung reproduces the card ---------------------------
    // At scale 0 the scan measures exactly the nominal bond line, so the map
    // must land back on the card's own R''. It is NOT bitwise: the fitted
    // pose is not exactly the identity, so the recovered thickness carries a
    // second-order rotation residual. That residual being ~1e-12 relative is
    // the evidence that the pose really does cancel between paired faces.
    let nominal_rung = &ladder[0];
    for thickness in &nominal_rung.thickness {
        let relative = (thickness - NOMINAL_BOND_THICKNESS_M).abs() / NOMINAL_BOND_THICKNESS_M;
        assert!(
            relative < 1e-9,
            "an undeformed bond line must recover its nominal thickness: {thickness:e} m ({relative:e} relative)"
        );
    }
    let nominal_conductance = 1.0 / CARD_RESISTANCE_M2K_PER_W;
    assert!(
        (nominal_rung.conductance_w_per_k - nominal_conductance).abs() < 1e-6,
        "the zero-wedge map must reproduce the card's conductance: {} vs {nominal_conductance}",
        nominal_rung.conductance_w_per_k
    );

    // ---- (b) direction and monotonicity ------------------------------------
    // A thicker measured bond line conducts less, jumps more, and raises the
    // hot side. Direction and monotonicity only; no magnitude is validated.
    for window in ladder.windows(2) {
        let (thin, thick) = (&window[0], &window[1]);
        assert!(
            thick.max_thickness_m() > thin.max_thickness_m(),
            "wedge scale {} must measure thicker than {}",
            thick.scale,
            thin.scale
        );
        assert!(
            thick.conductance_w_per_k < thin.conductance_w_per_k,
            "a thicker bond line must conduct LESS: {} !< {}",
            thick.conductance_w_per_k,
            thin.conductance_w_per_k
        );
        assert!(
            thick.mean_jump_k > thin.mean_jump_k,
            "a thicker bond line must jump MORE at fixed heat rate: {} !> {}",
            thick.mean_jump_k,
            thin.mean_jump_k
        );
        assert!(
            thick.max_temperature_k > thin.max_temperature_k,
            "a thicker bond line must raise the hot side: {} K !> {} K",
            thick.max_temperature_k,
            thin.max_temperature_k
        );
    }

    // A wedge is a genuine MAP: the severe rung's own faces must disagree.
    let severe = ladder.last().expect("ladder rung");
    let (lo, hi) = (severe.min_thickness_m(), severe.max_thickness_m());
    assert!(
        hi > lo * 1.05,
        "a wedge must vary across the bond line, not shift uniformly: {lo:e} to {hi:e} m"
    );

    // ---- a MAP is not its own average -------------------------------------
    // The faces conduct in PARALLEL, so the surface conductance is the sum of
    // A_f/R''_f — an arithmetic mean of conductances, i.e. a harmonic mean of
    // resistances. By Jensen that strictly exceeds the conductance a uniform
    // surface at the same MEAN thickness would have had, for any non-constant
    // map. This is the assertion that makes the per-face map load-bearing
    // rather than decorative: collapsing the scan to its average thickness
    // gives a DIFFERENT, pessimistic answer, so the map is carrying
    // information the average destroys.
    //
    // The margin is not cosmetic. A bare `>` does NOT discriminate: if the
    // seam collapsed the map to its mean, the two sides would differ only by
    // the float noise between `Σ(A_f·x)` and `(ΣA_f)·x` and the comparison
    // would fall either way. This fixture's wedge separates them by 1.5%
    // (scale 1) and 4.1% (scale 2), so 0.5% sits far above the noise and far
    // below the real effect. Verified by falsification: replacing every face
    // value with the mean thickness passes a bare `>` and fails this.
    let k_bond = NOMINAL_BOND_THICKNESS_M / CARD_RESISTANCE_M2K_PER_W;
    for rung in &ladder[1..] {
        let averaged = rung.averaged_conductance_w_per_k(k_bond);
        assert!(
            rung.conductance_w_per_k > averaged * MAP_VS_AVERAGE_MARGIN,
            "wedge scale {}: the map must not collapse to its own average — \
             {} vs {averaged} (ratio {:.6}, need > {MAP_VS_AVERAGE_MARGIN})",
            rung.scale,
            rung.conductance_w_per_k,
            rung.conductance_w_per_k / averaged
        );
    }
    // ...and at zero wedge the map IS constant, so the two must coincide.
    let flat = nominal_rung.averaged_conductance_w_per_k(k_bond);
    assert!(
        (nominal_rung.conductance_w_per_k - flat).abs() < 1e-9,
        "a constant map must agree with its average: {} vs {flat}",
        nominal_rung.conductance_w_per_k
    );
    // Distinct maps must produce distinct field identities.
    assert_ne!(ladder[0].field_identity, ladder[1].field_identity);
    assert_ne!(ladder[1].field_identity, ladder[2].field_identity);

    for marker in [
        "\"stage\":\"solve\"",
        "\"wedge_scale\":0.0",
        "\"wedge_scale\":2.0",
        "\"conductance_w_per_k\":",
    ] {
        assert!(log.contains(marker), "forensic log lost {marker:?}\n{log}");
    }
    println!("{log}");
}

#[test]
fn e2e_solve_stage_is_deterministic_across_independent_runs() {
    // The physics stage must be replayable to the same bits: same recovered
    // thickness map, same field identity, same interface conductance, same
    // hot-side temperature. A drift here would mean the solve picked up
    // ambient state that the chain's identities do not describe.
    let mesh = two_slab_mesh();
    let boundary = thermal_boundary(&mesh);
    let material = ConductivityModel::isotropic_declared(K_SOLID_W_PER_MK).expect("slab material");
    let source = ScalarField::Uniform(0.0);
    let card = bond_card();
    let pairs = oriented_pairs(&mesh);
    let registration = with_default_cx(stage_register);
    let fixture = BondFixture {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        card: &card,
        pairs: &pairs,
        registration: &registration,
    };

    let first = solve_rung(&fixture, 1.0);
    let second = solve_rung(&fixture, 1.0);

    assert_eq!(first.field_identity, second.field_identity);
    assert_eq!(first.field_registration_ref, second.field_registration_ref);
    assert_eq!(first.thickness.len(), second.thickness.len());
    for (a, b) in first.thickness.iter().zip(second.thickness.iter()) {
        assert_eq!(a.to_bits(), b.to_bits(), "thickness station drifted");
    }
    assert_eq!(
        first.conductance_w_per_k.to_bits(),
        second.conductance_w_per_k.to_bits(),
        "interface conductance drifted"
    );
    assert_eq!(
        first.max_temperature_k.to_bits(),
        second.max_temperature_k.to_bits(),
        "hot-side temperature drifted"
    );
    println!(
        "{{\"suite\":\"fs-diffreal-e2e/as-built-chain\",\"case\":\"solve-determinism\",\
         \"field\":\"{}\",\"t_hot_k\":{:.12},\"conductance_w_per_k\":{:.12}}}",
        first.field_identity, first.max_temperature_k, first.conductance_w_per_k
    );
}
