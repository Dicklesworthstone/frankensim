//! Joint finite-patch damage/friction: analytical force/work and native mechanics.
use fs_couple::cohesive_contact::CohesiveContactConfig;
use fs_couple::cohesive_contact::frictional_damage::*;
use fs_couple::cohesive_contact::patch::frictional_damage::*;
use fs_couple::cohesive_contact::reclosure::source::*;
use fs_couple::cohesive_contact::reclosure::{ReclosureActivation, ResidualCoulombSource};
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
#[path = "cohesive_frictional_patch/native.rs"]
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
                    source: "manufactured partially damaged friction fixture".into(),
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
            process: "manufactured damage-friction fixture".into(),
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

fn receiving(source: &ResolvedBilinearCohesive) -> [&MaterialStateId; 2] {
    [
        &source.surfaces()[0].material,
        &source.surfaces()[1].material,
    ]
}
fn local_config() -> DamageFrictionConfig {
    DamageFrictionConfig {
        contact: CohesiveContactConfig {
            max_iterations: 48,
            max_material_evaluations: 8192,
            max_backtracks: 16,
            jump_absolute_tolerance_m: 1e-14,
            jump_relative_tolerance: 1e-10,
            maximum_jump_increment_m: 0.01,
            maximum_force_n: 100.0,
            material_energy_tolerance_j_m2: 1e-8,
        },
        max_friction_iterations: 96,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-9,
        work_absolute_tolerance_j: 1e-13,
        work_relative_tolerance: 1e-8,
    }
}

fn config() -> DamageFrictionPatchConfig {
    DamageFrictionPatchConfig {
        max_points: 8,
        friction: local_config(),
        newton_damping: 1e-14,
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}; {t:e}");
}
fn point(source: &ResolvedBilinearCohesive, maximum: f64, g: [f64; 3]) -> CohesivePoint<'_> {
    source
        .law()
        .bonded_point()
        .trial([maximum, 0.0, 0.0], 1e-8)
        .unwrap()
        .point
        .trial(g, 1e-8)
        .unwrap()
        .point
}
fn dense(n: usize) -> Vec<f64> {
    let shape = [1.0, 0.3, -0.2, 0.7, -0.1, 0.5, 0.4, 0.2, -0.6];
    (0..3 * n)
        .flat_map(|i| {
            (0..3 * n).map(move |j| 1e-4 * shape[i] * shape[j] + if i == j { 2e-4 } else { 0.0 })
        })
        .collect()
}
fn free_for(s: &[f64], target: &[[f64; 3]], force: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let dim = target.len() * 3;
    (0..target.len())
        .map(|p| {
            core::array::from_fn(|i| {
                target[p][i]
                    + (0..dim)
                        .map(|j| s[(3 * p + i) * dim + j] * force[j / 3][j % 3])
                        .sum::<f64>()
            })
        })
        .collect()
}
fn expected(point: CohesivePoint<'_>, end: [f64; 3], area: f64, mu: f64) -> [f64; 3] {
    let old = point.jump_m();
    let slip = (end[1] - old[1]).hypot(end[2] - old[2]);
    let compression = -area * COMPRESSION * (0.5 * old[0] + 0.5 * end[0]);
    let cap = mu * point.damage() * compression;
    [
        -compression,
        area * (1.0 - point.damage()) * 1e9 * 0.5 * (old[1] + end[1])
            + cap * (end[1] - old[1]) / slip,
        area * (1.0 - point.damage()) * 1e9 * 0.5 * (old[2] + end[2])
            + cap * (end[2] - old[2]) / slip,
    ]
}

#[test]
fn g1_heterogeneous_partial_bonds_use_every_cross_site_and_cross_axis_term() {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let cards = [source_card(0.6, 0.3), source_card(0.8, 0.5)];
    let sources: Vec<_> = cards
        .iter()
        .map(|c| {
            resolve_bilinear_cohesive(c, &query, MaterialPropertySelection::SingleClaimOnly, 0.02)
                .unwrap()
        })
        .collect();
    let friction: Vec<_> = cards
        .iter()
        .map(|c| {
            resolve_residual_coulomb(c, &query, MaterialPropertySelection::SingleClaimOnly).unwrap()
        })
        .collect();
    let areas = [AREA * 0.3, AREA * 0.7];
    let points = [
        point(&sources[0], 4e-4, [-2e-4, 0.0, 0.0]),
        point(&sources[1], 6e-4, [-1.5e-4, 0.0, 0.0]),
    ];
    let target = [[-2e-4, 5e-5, -2e-5], [-1.5e-4, -3e-5, 7e-5]];
    let expected = [
        expected(points[0], target[0], areas[0], 0.3),
        expected(points[1], target[1], areas[1], 0.5),
    ];
    let s = dense(2);
    let mut owner = DamageFrictionPatchLoadLine::new(
        (0..2)
            .map(|p| DamageFrictionPatchSite {
                cohesive: &sources[p],
                friction: &friction[p],
                area_m2: areas[p],
            })
            .collect(),
        s.clone(),
        config(),
    )
    .unwrap();
    let old: Vec<_> = (0..2)
        .map(|p| {
            owner
                .adopt_point(
                    p,
                    points[p],
                    receiving(&sources[p]),
                    ReclosureActivation::AtRest,
                )
                .unwrap()
        })
        .collect();
    let free = free_for(&s, &target, &expected);
    let step = owner.trial(&old, &free, 0.01, || true).unwrap();
    assert_eq!(step.kinetic_switches, 2);
    for p in 0..2 {
        assert_eq!(step.regimes[p], DamageFrictionRegime::Sliding);
        let audit = step
            .audit_point(
                p,
                old[p],
                step.jump_m[p],
                0.01,
                step.mean_force_n[p],
                step.regimes[p],
            )
            .unwrap();
        assert_eq!(audit.state.point().damage(), points[p].damage());
        assert_eq!(audit.material.equilibrium.energy.dissipated, 0.0);
        let slip =
            (target[p][1] - points[p].jump_m()[1]).hypot(target[p][2] - points[p].jump_m()[2]);
        let cap = [0.3, 0.5][p] * points[p].damage() * areas[p] * COMPRESSION * (-target[p][0]);
        close(audit.friction_dissipation_j, cap * slip, 1e-12);
        for i in 0..3 {
            close(step.jump_m[p][i], target[p][i], 1e-12);
            close(step.mean_force_n[p][i], expected[p][i], 1e-9);
        }
    }
    assert!(step.material_evaluations <= config().friction.contact.max_material_evaluations);
}

#[test]
fn g1_rank_deficient_coincident_quadrature_preserves_area_and_shared_motion() {
    let (source, friction) = sources();
    let original = point(&source, 4e-4, [-2e-4, 0.0, 0.0]);
    let local = [
        [2e-4, 3e-5, -2e-5],
        [3e-5, 3e-4, 4e-5],
        [-2e-5, 4e-5, 2.5e-4],
    ];
    for fractions in [vec![1.0], vec![0.2, 0.3, 0.5]] {
        let n = fractions.len();
        let dim = 3 * n;
        let s: Vec<_> = (0..dim)
            .flat_map(|i| (0..dim).map(move |j| local[i % 3][j % 3]))
            .collect();
        let mut owner = DamageFrictionPatchLoadLine::new(
            fractions
                .iter()
                .map(|w| DamageFrictionPatchSite {
                    cohesive: &source,
                    friction: &friction,
                    area_m2: AREA * w,
                })
                .collect(),
            s.clone(),
            config(),
        )
        .unwrap();
        let old: Vec<_> = (0..n)
            .map(|p| {
                owner
                    .adopt_point(p, original, receiving(&source), ReclosureActivation::AtRest)
                    .unwrap()
            })
            .collect();
        let total = [-0.4, 0.03, -0.02];
        let reactions: Vec<_> = fractions.iter().map(|w| total.map(|x| x * w)).collect();
        let fixed = vec![original.jump_m(); n];
        let free = free_for(&s, &fixed, &reactions);
        let step = owner.trial(&old, &free, 0.01, || true).unwrap();
        for p in 0..n {
            assert_eq!(step.regimes[p], DamageFrictionRegime::Sticking);
            for i in 0..3 {
                close(step.jump_m[p][i], original.jump_m()[i], 1e-13);
            }
            let a = step
                .audit_point(
                    p,
                    old[p],
                    step.jump_m[p],
                    0.01,
                    step.mean_force_n[p],
                    step.regimes[p],
                )
                .unwrap();
            assert_eq!(a.friction_dissipation_j, 0.0);
        }
        for i in 0..3 {
            close(step.mean_force_n.iter().map(|f| f[i]).sum(), total[i], 1e-9);
        }
        let end = [-2e-4, 8e-5, -6e-5];
        let force = expected(original, end, AREA, 0.3);
        let forces: Vec<_> = fractions.iter().map(|w| force.map(|x| w * x)).collect();
        let free = free_for(&s, &vec![end; n], &forces);
        let step = owner.trial(&old, &free, 0.01, || true).unwrap();
        let mut loss = 0.0;
        for p in 0..n {
            for i in 0..3 {
                close(step.jump_m[p][i], end[i], 1e-12);
                close(step.mean_force_n[p][i], forces[p][i], 1e-9);
            }
            loss += step
                .audit_point(
                    p,
                    old[p],
                    step.jump_m[p],
                    0.01,
                    step.mean_force_n[p],
                    step.regimes[p],
                )
                .unwrap()
                .friction_dissipation_j;
        }
        close(loss, 0.3 * 0.625 * 0.4 * 1e-4, 1e-12);
    }
}

#[test]
fn g3_one_site_reuses_existing_law_and_source_budget_refusals_are_atomic() {
    let (source, friction) = sources();
    let original = point(&source, 4e-4, [-2e-4, 0.0, 0.0]);
    let s = [
        [2e-4, 3e-5, -2e-5],
        [3e-5, 3e-4, 4e-5],
        [-2e-5, 4e-5, 2.5e-4],
    ];
    let mut single =
        DamageFrictionLoadLine::new(&source, &friction, AREA, s, local_config()).unwrap();
    let a = single
        .adopt(original, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let mut patch = DamageFrictionPatchLoadLine::new(
        vec![DamageFrictionPatchSite {
            cohesive: &source,
            friction: &friction,
            area_m2: AREA,
        }],
        s.into_iter().flatten().collect(),
        config(),
    )
    .unwrap();
    let b = patch
        .adopt_point(0, original, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let free = [[-3e-4, 1e-4, 6e-5]];
    assert!(matches!(
        patch.trial(&[b], &free, 0.01, || false),
        Err(DamageFrictionPatchError::Cancelled)
    ));
    assert!(
        patch
            .adopt_point(
                0,
                original,
                [receiving(&source)[1], receiving(&source)[0]],
                ReclosureActivation::AtRest
            )
            .is_err()
    );
    let expected = single.trial(a, free[0], 0.01, || true).unwrap();
    let actual = patch.trial(&[b], &free, 0.01, || true).unwrap();
    for i in 0..3 {
        close(
            actual.jump_m[0][i],
            expected.audit.state.point().jump_m()[i],
            2e-12,
        );
        close(actual.mean_force_n[0][i], expected.mean_force_n[i], 2e-9);
    }
    let force = actual.mean_force_n[0];
    let jump = actual.jump_m[0];
    let evaluations = actual.material_evaluations;
    let again = patch.trial(&[b], &free, 0.01, || true).unwrap();
    assert_eq!(again.mean_force_n[0], force);
    assert_eq!(again.jump_m[0], jump);
    assert_eq!(again.material_evaluations, evaluations);
    let mut exhausted = config();
    exhausted.friction.contact.max_material_evaluations = 2;
    let mut tiny = DamageFrictionPatchLoadLine::new(
        vec![DamageFrictionPatchSite {
            cohesive: &source,
            friction: &friction,
            area_m2: AREA,
        }],
        s.into_iter().flatten().collect(),
        exhausted,
    )
    .unwrap();
    let old = tiny
        .adopt_point(0, original, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    assert!(matches!(
        tiny.trial(&[old], &free, 0.01, || true),
        Err(DamageFrictionPatchError::Budget)
    ));
    let foreign = resolve_residual_coulomb(
        &source_card(0.7, 0.3),
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        DamageFrictionPatchLoadLine::new(
            vec![DamageFrictionPatchSite {
                cohesive: &source,
                friction: &foreign,
                area_m2: AREA
            }],
            s.into_iter().flatten().collect(),
            config()
        )
        .is_err()
    );
}

#[test]
fn g1_prescribed_slip_in_a_zero_compliance_direction_still_reaches_the_coulomb_boundary() {
    let (source, friction) = sources();
    let original = point(&source, 4e-4, [-2e-4, 0.0, 0.0]);
    let end = [-2e-4, 1e-4, -5e-5];
    let expected = expected(original, end, AREA, 0.3);
    for damping in [1e-14, 1e-10] {
        let mut config = config();
        config.newton_damping = damping;
        let mut owner = DamageFrictionPatchLoadLine::new(
            vec![DamageFrictionPatchSite {
                cohesive: &source,
                friction: &friction,
                area_m2: AREA,
            }],
            vec![0.0; 9],
            config,
        )
        .unwrap();
        let old = owner
            .adopt_point(0, original, receiving(&source), ReclosureActivation::AtRest)
            .unwrap();
        let step = owner.trial(&[old], &[end], 0.01, || true).unwrap();
        assert_eq!(step.jump_m[0], end);
        assert_eq!(step.regimes[0], DamageFrictionRegime::Sliding);
        for i in 0..3 {
            close(step.mean_force_n[0][i], expected[i], 1e-10);
        }
        let audit = step
            .audit_point(0, old, end, 0.01, step.mean_force_n[0], step.regimes[0])
            .unwrap();
        assert!(audit.friction_dissipation_j > 0.0);
        assert!(audit.work_residual_j.abs() < 1e-13);
    }
}
