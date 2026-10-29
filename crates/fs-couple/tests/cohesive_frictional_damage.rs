//! Independent force/work and original-state tests of the declared rigid mixture.
use fs_couple::cohesive_contact::frictional_damage::*;
use fs_couple::cohesive_contact::reclosure::source::*;
use fs_couple::cohesive_contact::reclosure::{ReclosureActivation, ResidualCoulombSource};
use fs_couple::cohesive_contact::{CohesiveContactConfig, CohesiveLoadLine};
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

#[path = "cohesive_frictional_damage/native.rs"]
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
fn config() -> DamageFrictionConfig {
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
fn compliance() -> [[f64; 3]; 3] {
    [[2e-5, 2e-6, -1e-6], [2e-6, 3e-5, 4e-6], [-1e-6, 4e-6, 4e-5]]
}
fn free_for(end: [f64; 3], force: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| end[i] + (0..3).map(|j| compliance()[i][j] * force[j]).sum::<f64>())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn partial(source: &ResolvedBilinearCohesive) -> CohesivePoint<'_> {
    let p = source
        .law()
        .bonded_point()
        .trial([4e-4, 0.0, 0.0], 1e-8)
        .unwrap()
        .point;
    assert!(p.damage() > 0.0 && p.damage() < 1.0);
    p.trial([-2e-4, 0.0, 0.0], 1e-8).unwrap().point
}

#[test]
fn g1_partial_damage_slides_with_independent_cohesive_and_damaged_area_forces() {
    let (source, friction) = sources();
    let mut owner =
        DamageFrictionLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
    let point = partial(&source);
    close(point.damage(), 0.625, 1e-14);
    let old = owner
        .adopt(point, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let end = [-2e-4, 6e-5, 8e-5];
    let pressure = AREA * COMPRESSION * 2e-4;
    let cap = 0.3 * 0.625 * pressure;
    let elastic = 0.5 * AREA * (1.0 - 0.625) * 1e9 * 1e-4;
    let expected = [-pressure, 0.6 * (elastic + cap), 0.8 * (elastic + cap)];
    let step = owner
        .trial(old, free_for(end, expected), 0.01, || true)
        .unwrap();
    assert_eq!(step.audit.regime, DamageFrictionRegime::Sliding);
    close(
        step.audit.damaged_compression_force_n,
        0.625 * pressure,
        1e-11,
    );
    close(step.audit.friction_dissipation_j, cap * 1e-4, 1e-13);
    assert_eq!(step.audit.material.equilibrium.energy.dissipated, 0.0);
    assert_eq!(
        step.audit.material.equilibrium.point.damage(),
        point.damage()
    );
    for i in 0..3 {
        close(step.audit.state.point().jump_m()[i], end[i], 1e-12);
        close(step.mean_force_n[i], expected[i], 1e-9);
    }
    let work = (0..3)
        .map(|i| expected[i] * (end[i] - point.jump_m()[i]))
        .sum::<f64>();
    close(
        work,
        AREA * step.audit.material.total_work_j_m2() + step.audit.friction_dissipation_j,
        1e-13,
    );
    let reaction = [-pressure, 0.03, -0.02];
    let stuck = owner
        .trial(old, free_for(point.jump_m(), reaction), 0.01, || true)
        .unwrap();
    assert_eq!(stuck.audit.regime, DamageFrictionRegime::Sticking);
    assert_eq!(stuck.audit.friction_dissipation_j, 0.0);
    for i in 0..3 {
        close(
            stuck.audit.state.point().jump_m()[i],
            point.jump_m()[i],
            1e-13,
        );
    }
}

#[test]
fn g1_damage_grows_while_sliding_and_force_integrates_the_entire_damage_path() {
    let (source, friction) = sources();
    let mut owner =
        DamageFrictionLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
    let point = source
        .law()
        .bonded_point()
        .trial([-2e-4, 0.0, 0.0], 1e-8)
        .unwrap()
        .point;
    let old = owner
        .adopt(point, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let length: f64 = 6e-4;
    let onset: f64 = 2e-4;
    let failure: f64 = 1e-3;
    // Analytic scalar envelope integral, independently of all vector path helpers.
    let damage_mean =
        failure / (failure - onset) * (length - onset - onset * (length / onset).ln()) / length;
    let cohesive_work = 0.5 * 1e9 * onset * onset
        + 2e5 / (failure - onset)
            * (failure * (length - onset) - 0.5 * (length * length - onset * onset));
    let pressure = AREA * COMPRESSION * 2e-4;
    let cap = 0.3 * pressure * damage_mean;
    let tangent = AREA * cohesive_work / length + cap;
    let expected = [-pressure, 0.6 * tangent, 0.8 * tangent];
    let end = [-2e-4, 0.6 * length, 0.8 * length];
    let step = owner
        .trial(old, free_for(end, expected), 0.01, || true)
        .unwrap();
    assert_eq!(step.audit.regime, DamageFrictionRegime::Sliding);
    assert!(step.audit.state.point().damage() > 0.0 && step.audit.state.point().damage() < 1.0);
    assert!(step.audit.material.equilibrium.energy.dissipated > 0.0);
    close(
        step.audit.damaged_compression_force_n,
        pressure * damage_mean,
        1e-9,
    );
    close(step.audit.friction_dissipation_j, cap * length, 1e-12);
    close(
        AREA * step.audit.material.total_work_j_m2(),
        AREA * cohesive_work,
        1e-12,
    );
    assert!(
        step.audit.damaged_compression_force_n < pressure * step.audit.state.point().damage() * 0.7
    );
    for i in 0..3 {
        close(step.mean_force_n[i], expected[i], 2e-9);
        close(step.audit.state.point().jump_m()[i], end[i], 1e-12);
    }
}

#[test]
fn g1_intact_open_and_fully_failed_limits_keep_the_original_material_owner() {
    let (source, friction) = sources();
    let mut owner =
        DamageFrictionLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
    let mut original = CohesiveLoadLine::new(AREA, compliance(), config().contact).unwrap();
    let point = source.law().bonded_point();
    let old = owner
        .adopt(point, receiving(&source), ReclosureActivation::AtRest)
        .unwrap();
    let free = [-2e-5, 3e-5, -1e-5];
    let actual = owner.trial(old, free, 0.01, || true).unwrap();
    let plain = original.trial(point, free, 0.01, || true).unwrap();
    assert_eq!(actual.audit.regime, DamageFrictionRegime::Inactive);
    assert_eq!(actual.audit.friction_dissipation_j, 0.0);
    for i in 0..3 {
        close(actual.mean_force_n[i], plain.mean_force_n[i], 1e-10);
    }
    for damage_jump in [[4e-4, 0.0, 0.0], [0.002, 0.0, 0.0]] {
        let point = source
            .law()
            .bonded_point()
            .trial(damage_jump, 1e-8)
            .unwrap()
            .point;
        let old = owner
            .adopt(point, receiving(&source), ReclosureActivation::AtRest)
            .unwrap();
        let open = owner
            .trial(old, [damage_jump[0] + 1e-5, 1e-4, -3e-5], 0.01, || true)
            .unwrap();
        assert_eq!(open.audit.regime, DamageFrictionRegime::Inactive);
        assert_eq!(open.audit.friction_dissipation_j, 0.0);
        assert!(open.audit.state.point().damage() >= point.damage());
    }
    let failed = source
        .law()
        .bonded_point()
        .trial([0.002, 0.0, 0.0], 1e-8)
        .unwrap()
        .point
        .trial([-2e-4, 0.0, 0.0], 1e-8)
        .unwrap()
        .point;
    let old = owner
        .adopt(failed, receiving(&source), ReclosureActivation::Sliding)
        .unwrap();
    let end = [-3e-4, 1e-4, -2e-4];
    let p = AREA * COMPRESSION * 2.5e-4;
    let slip = 1e-4_f64.hypot(-2e-4);
    let expected = [-p, 0.3 * p * end[1] / slip, 0.3 * p * end[2] / slip];
    let step = owner
        .trial(old, free_for(end, expected), 0.01, || true)
        .unwrap();
    assert!(step.audit.state.point().failed());
    assert_eq!(step.audit.material.equilibrium.energy.dissipated, 0.0);
    close(step.audit.friction_dissipation_j, 0.3 * p * slip, 1e-12);
    close(step.audit.damaged_compression_force_n, p, 1e-10);
    for i in 0..3 {
        close(step.mean_force_n[i], expected[i], 1e-9);
    }
}

#[test]
fn g3_coefficient_sources_controls_capacity_and_all_refusals_preserve_original_state() {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    for mu in [0.0, 0.1, 0.3, 0.55] {
        let card = source_card(0.6, mu);
        let source = resolve_bilinear_cohesive(
            &card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            0.02,
        )
        .unwrap();
        let friction =
            resolve_residual_coulomb(&card, &query, MaterialPropertySelection::SingleClaimOnly)
                .unwrap();
        let mut owner =
            DamageFrictionLoadLine::new(&source, &friction, AREA, compliance(), config()).unwrap();
        let point = partial(&source);
        let old = owner
            .adopt(point, receiving(&source), ReclosureActivation::Sliding)
            .unwrap();
        let end = [-2e-4, 1e-4, 0.0];
        let cap = mu * 0.625 * AREA * COMPRESSION * 2e-4;
        let expected = [-0.4, 0.01875 + cap, 0.0];
        let free = free_for(end, expected);
        assert!(matches!(
            owner.trial(old, free, 0.01, || false),
            Err(DamageFrictionError::Cancelled)
        ));
        assert!(matches!(
            owner.adopt(
                point,
                [receiving(&source)[1], receiving(&source)[0]],
                ReclosureActivation::AtRest
            ),
            Err(DamageFrictionError::Source(_))
        ));
        let before = old;
        let accepted = owner.trial(old, free, 0.01, || true).unwrap();
        assert_eq!(old, before);
        close(accepted.audit.friction_force_n[0], cap, 1e-9);
        let mut forged = accepted.mean_force_n;
        forged[2] += 0.01;
        assert!(matches!(
            owner.audit(
                old,
                accepted.audit.state.point().jump_m(),
                0.01,
                forged,
                DamageFrictionRegime::Sliding
            ),
            Err(DamageFrictionError::ForceMismatch { .. })
        ));
        let retry = owner.trial(old, free, 0.01, || true).unwrap();
        assert_eq!(accepted, retry);
        assert!(accepted.material_evaluations <= config().contact.max_material_evaluations);
        let mut budget = config();
        budget.contact.max_material_evaluations = 3;
        let mut tiny =
            DamageFrictionLoadLine::new(&source, &friction, AREA, compliance(), budget).unwrap();
        let state = tiny
            .adopt(point, receiving(&source), ReclosureActivation::Sliding)
            .unwrap();
        assert!(matches!(
            tiny.trial(state, free, 0.01, || true),
            Err(DamageFrictionError::Budget)
        ));
    }
    let (source, friction) = sources();
    let changed = resolve_residual_coulomb(
        &source_card(0.7, 0.3),
        &query,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        DamageFrictionLoadLine::new(&source, &changed, AREA, compliance(), config()),
        Err(DamageFrictionError::Source(_))
    ));
    let mut nonsymmetric = compliance();
    nonsymmetric[0][1] *= 2.0;
    assert!(DamageFrictionLoadLine::new(&source, &friction, AREA, nonsymmetric, config()).is_err());
}
