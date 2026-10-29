//! G1/G3/G4: original failed-bond compression, residual friction and mechanics.
use fs_couple::cohesive_contact::reclosure::source::*;
use fs_couple::cohesive_contact::reclosure::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::CohesivePoint;
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Dims, Pressure};

const AREA: f64 = 1e-6;
const COMPRESSION: f64 = 2e9;

#[path = "cohesive_reclosure/native.rs"]
mod native;

fn source_card(static_mu: f64, kinetic_mu: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            COMPRESSION,
        ),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, static_mu),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, kinetic_mu),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured failed-bond friction fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured reclosure fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit interface".into(),
            third_body: None,
            environment: "fixed source query".into(),
            history: "original bond and residual surface data".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn sources() -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let card = source_card(0.6, 0.3);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let cohesive = resolve_bilinear_cohesive(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    let friction =
        resolve_residual_coulomb(&card, &point, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    (cohesive, friction)
}

#[test]
fn g0_g1_residual_source_is_explicit_fixed_and_can_describe_frictionless_failure() {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let selection = || MaterialPropertySelection::SingleClaimOnly;
    let card = source_card(0.6, 0.3);
    let source = resolve_bilinear_cohesive(&card, &query, selection(), 0.02).unwrap();
    // Cohesive strength and fracture energy do not substitute for either mu.
    assert!(ResidualCoulombSource::try_from_resolved(&card, source.resolved().clone()).is_err());
    let shifted = resolve_residual_coulomb(
        &card,
        &QueryPoint::new().with("T", 301.0).unwrap(),
        selection(),
    )
    .unwrap();
    assert!(matches!(
        ReclosureLoadLine::new(&source, &shifted, AREA, compliance(), config()),
        Err(ReclosureError::Source(_))
    ));
    let foreign = resolve_residual_coulomb(&source_card(0.7, 0.3), &query, selection()).unwrap();
    assert!(matches!(
        ReclosureLoadLine::new(&source, &foreign, AREA, compliance(), config()),
        Err(ReclosureError::Source(_))
    ));
    assert!(resolve_residual_coulomb(&source_card(0.2, 0.3), &query, selection()).is_err());
    let card = source_card(0.0, 0.0);
    let source = resolve_bilinear_cohesive(&card, &query, selection(), 0.02).unwrap();
    let friction = resolve_residual_coulomb(&card, &query, selection()).unwrap();
    let owner = ReclosureLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
    let old = [-1e-4, 0.003, -0.001];
    let end = [-2e-4, 0.0031, -0.0012];
    let expected = [AREA * COMPRESSION * 0.5 * (old[0] + end[0]), 0.0, 0.0];
    let state = owner
        .adopt(
            point(&source, old),
            receiving(&source),
            ReclosureActivation::Sliding,
        )
        .unwrap();
    let step = owner
        .trial(state, free_for(compliance(), end, expected), 0.01, || true)
        .unwrap();
    assert_eq!(step.audit.regime, ReclosureRegime::Sliding);
    assert_eq!(step.audit.friction_dissipation_j, 0.0);
    for i in 0..3 {
        close(step.audit.state.point().jump_m()[i], end[i], 1e-13);
        close(step.mean_force_n[i], expected[i], 1e-12);
    }
}
fn receiving(source: &ResolvedBilinearCohesive) -> [&MaterialStateId; 2] {
    [
        &source.surfaces()[0].material,
        &source.surfaces()[1].material,
    ]
}
fn point(source: &ResolvedBilinearCohesive, jump: [f64; 3]) -> CohesivePoint<'_> {
    let failed = source
        .law()
        .bonded_point()
        .trial([0.002, 0.0, 0.0], 1e-8)
        .unwrap()
        .point;
    assert!(failed.failed());
    failed.trial(jump, 1e-8).unwrap().point
}
fn config() -> ReclosureConfig {
    ReclosureConfig {
        max_iterations: 128,
        max_evaluations: 4096,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-10,
        maximum_jump_increment_m: 0.01,
        maximum_force_n: 100.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-9,
        material_energy_tolerance_j_m2: 1e-8,
        work_absolute_tolerance_j: 1e-13,
        work_relative_tolerance: 1e-8,
        minimum_monotonicity_margin: 0.02,
    }
}
fn compliance() -> [[f64; 3]; 3] {
    [[2e-5, 2e-6, -1e-6], [2e-6, 3e-5, 4e-6], [-1e-6, 4e-6, 4e-5]]
}
fn mul(s: [[f64; 3]; 3], f: [f64; 3]) -> [f64; 3] {
    s.map(|r| r[0] * f[0] + r[1] * f[1] + r[2] * f[2])
}
fn free_for(s: [[f64; 3]; 3], end: [f64; 3], force: [f64; 3]) -> [f64; 3] {
    let response = mul(s, force);
    core::array::from_fn(|i| end[i] + response[i])
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute,
        "{a:e} != {b:e}; tolerance {absolute:e}"
    );
}

#[test]
fn g1_static_and_occupied_kinetic_branches_use_the_same_original_pressure() {
    let (source, friction) = sources();
    let s = [[1e-5, 0.0, 0.0], [0.0, 2e-5, 0.0], [0.0, 0.0, 3e-5]];
    let owner = ReclosureLoadLine::new(&source, &friction, AREA, s, config()).unwrap();
    let old = [-1e-4, 0.003, 0.0];
    let p = point(&source, old);
    let free = free_for(s, old, [-0.2, 0.08, 0.0]);
    let static_state = owner
        .adopt(p, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let sticking = owner.trial(static_state, free, 0.01, || true).unwrap();
    assert_eq!(sticking.audit.regime, ReclosureRegime::Sticking);
    close(sticking.mean_force_n[1], 0.08, 2e-10);
    assert_eq!(sticking.audit.friction_dissipation_j, 0.0);
    let moving = owner
        .adopt(p, receiving(&source), ReclosureActivation::Sliding)
        .unwrap();
    let sliding = owner.trial(moving, free, 0.01, || true).unwrap();
    assert_eq!(sliding.audit.regime, ReclosureRegime::Sliding);
    close(sliding.mean_force_n[1], 0.3 * 0.2, 2e-10);
    close(
        sliding.audit.state.point().jump_m()[1] - old[1],
        (0.08 - 0.06) * s[1][1],
        1e-14,
    );
    assert!(sliding.audit.friction_dissipation_j > 0.0);
    let stop_free = free_for(s, sliding.audit.state.point().jump_m(), [-0.2, 0.04, 0.0]);
    let stopped = owner
        .trial(sliding.audit.state, stop_free, 0.01, || true)
        .unwrap();
    assert_eq!(stopped.audit.regime, ReclosureRegime::Sticking);
    close(
        stopped.audit.state.point().dissipated_j_m2(),
        p.dissipated_j_m2(),
        0.0,
    );
}

#[test]
fn g1_full_cross_axis_load_line_matches_manufactured_sliding_work() {
    let (source, friction) = sources();
    let s = compliance();
    let owner = ReclosureLoadLine::new(&source, &friction, AREA, s, config()).unwrap();
    let old = [-1e-4, 0.003, -0.001];
    let end = [-1.2e-4, old[1] + 3e-5, old[2] - 4e-5];
    let pressure = -AREA * COMPRESSION * 0.5 * (old[0] + end[0]);
    let force = [-pressure, 0.3 * pressure * 0.6, -0.3 * pressure * 0.8];
    let state = owner
        .adopt(
            point(&source, old),
            receiving(&source),
            ReclosureActivation::Sliding,
        )
        .unwrap();
    let result = owner
        .trial(state, free_for(s, end, force), 0.01, || true)
        .unwrap();
    assert_eq!(result.audit.regime, ReclosureRegime::Sliding);
    for i in 0..3 {
        close(result.audit.state.point().jump_m()[i], end[i], 2e-13);
        close(result.mean_force_n[i], force[i], 3e-10);
        assert!(result.compatibility_residual_m[i].abs() <= result.compatibility_tolerance_m);
    }
    close(
        result.audit.friction_dissipation_j,
        0.3 * pressure * 5e-5,
        1e-14,
    );
    let storage = 0.5 * AREA * COMPRESSION * (end[0] * end[0] - old[0] * old[0]);
    close(
        AREA * result.audit.material.total_work_j_m2(),
        storage,
        1e-17,
    );
    assert!(result.audit.work_residual_j.abs() < 1e-13);
    assert!(core::ptr::eq(owner.cohesive_source(), &source));
    assert!(core::ptr::eq(owner.friction_source(), &friction));
    assert_eq!(result.audit.state.point().damage(), 1.0);
    assert_eq!(result.audit.material.equilibrium.energy.dissipated, 0.0);
    assert_eq!(result.audit.material.viscous_dissipated_j_m2, 0.0);
}

#[test]
fn g1_contact_entry_exit_use_integrated_pressure_and_open_paths_have_no_friction() {
    let (source, friction) = sources();
    let s = compliance();
    let owner = ReclosureLoadLine::new(&source, &friction, AREA, s, config()).unwrap();
    for (start, finish) in [(-1e-4, 2e-4), (1e-4, -2e-4)] {
        let old = [start, 0.003, -0.001];
        let end = [finish, old[1] + 1e-4, old[2]];
        // Independently use the compression energy difference divided by actual opening.
        let energy = |g: f64| 0.5 * AREA * COMPRESSION * g.min(0.0).powi(2);
        let pressure = -(energy(finish) - energy(start)) / (finish - start);
        let force = [-pressure, 0.3 * pressure, 0.0];
        let state = owner
            .adopt(
                point(&source, old),
                receiving(&source),
                ReclosureActivation::Sliding,
            )
            .unwrap();
        let result = owner
            .trial(state, free_for(s, end, force), 0.01, || true)
            .unwrap();
        close(result.mean_force_n[0], -pressure, 2e-10);
        close(
            result.audit.friction_dissipation_j,
            0.3 * pressure * 1e-4,
            2e-14,
        );
        if finish > 0.0 {
            assert_eq!(result.audit.state.occupancy(), ReclosureRegime::Open);
            let mut free = result.audit.state.point().jump_m();
            free[1] += 2e-4;
            let open = owner
                .trial(result.audit.state, free, 0.01, || true)
                .unwrap();
            assert_eq!(open.mean_force_n, [0.0; 3]);
            assert_eq!(open.audit.regime, ReclosureRegime::Open);
            assert_eq!(open.audit.friction_dissipation_j, 0.0);
            close(
                open.audit.state.point().dissipated_j_m2(),
                state.point().dissipated_j_m2(),
                0.0,
            );
        }
    }
}

#[test]
fn g3_tangent_basis_rotation_preserves_full_coupling_and_fracture_history() {
    let (source, friction) = sources();
    let angle: f64 = 0.71;
    let (sin, cos) = angle.sin_cos();
    let r = [[1.0, 0.0, 0.0], [0.0, cos, -sin], [0.0, sin, cos]];
    let s = compliance();
    let rotated = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .flat_map(|k| (0..3).map(move |l| r[i][k] * s[k][l] * r[j][l]))
                .sum()
        })
    });
    let old = [-1e-4, 0.003, -0.001];
    let free = [-1.25e-4, 0.00304, -0.00107];
    let a = ReclosureLoadLine::new(&source, &friction, AREA, s, config()).unwrap();
    let b = ReclosureLoadLine::new(&source, &friction, AREA, rotated, config()).unwrap();
    let initial_a = a
        .adopt(
            point(&source, old),
            receiving(&source),
            ReclosureActivation::Sliding,
        )
        .unwrap();
    let initial_b = b
        .adopt(
            point(&source, mul(r, old)),
            receiving(&source),
            ReclosureActivation::Sliding,
        )
        .unwrap();
    let result_a = a.trial(initial_a, free, 0.01, || true).unwrap();
    let result_b = b.trial(initial_b, mul(r, free), 0.01, || true).unwrap();
    let force = mul(r, result_a.mean_force_n);
    let jump = mul(r, result_a.audit.state.point().jump_m());
    for i in 0..3 {
        close(result_b.mean_force_n[i], force[i], 4e-10);
        close(result_b.audit.state.point().jump_m()[i], jump[i], 2e-13);
    }
    close(
        result_b.audit.friction_dissipation_j,
        result_a.audit.friction_dissipation_j,
        2e-14,
    );
    close(
        result_b.audit.state.point().dissipated_j_m2(),
        result_a.audit.state.point().dissipated_j_m2(),
        0.0,
    );
}

#[test]
fn g0_g4_source_full_failure_feedback_budget_and_cancellation_refuse_atomically() {
    let (source, friction) = sources();
    let owner = ReclosureLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
    let p = point(&source, [-1e-4, 0.003, -0.001]);
    assert!(matches!(
        owner.adopt(
            source.law().bonded_point(),
            receiving(&source),
            ReclosureActivation::AtRest
        ),
        Err(ReclosureError::NotFullyFailed)
    ));
    let material = receiving(&source);
    assert!(
        owner
            .adopt(p, [material[1], material[0]], ReclosureActivation::AtRest)
            .is_err()
    );
    let other = source.clone();
    assert!(matches!(
        owner.adopt(
            point(&other, p.jump_m()),
            material,
            ReclosureActivation::AtRest
        ),
        Err(ReclosureError::StateMismatch)
    ));
    let strong = [[1e-5, 9e-6, 0.0], [9e-6, 1e-5, 0.0], [0.0, 0.0, 1e-5]];
    assert!(matches!(
        ReclosureLoadLine::new(&source, &friction, AREA, strong, config()),
        Err(ReclosureError::CouplingUnestablished { .. })
    ));
    let state = owner
        .adopt(p, material, ReclosureActivation::Sliding)
        .unwrap();
    let free = [-1.1e-4, 0.00308, -0.00104];
    let expected = owner.trial(state, free, 0.01, || true).unwrap();
    assert!(matches!(
        owner.trial(state, free, 0.01, || false),
        Err(ReclosureError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        owner.trial(state, free, 0.01, || {
            polls += 1;
            polls < 7
        }),
        Err(ReclosureError::Cancelled)
    ));
    let mut short = config();
    short.max_evaluations = 2;
    let limited = ReclosureLoadLine::new(&source, &friction, AREA, compliance(), short).unwrap();
    let limited_state = limited
        .adopt(p, material, ReclosureActivation::Sliding)
        .unwrap();
    assert!(matches!(
        limited.trial(limited_state, free, 0.01, || true),
        Err(ReclosureError::Budget)
    ));
    assert_eq!(owner.trial(state, free, 0.01, || true).unwrap(), expected);
    assert_eq!(state.point(), p);
    let mut bad = expected.audit.state.point().jump_m();
    bad[1] += 1e-5;
    assert!(
        owner
            .audit(
                state,
                bad,
                0.01,
                expected.mean_force_n,
                expected.audit.regime
            )
            .is_err()
    );
}
