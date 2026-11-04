//! Heterogeneous source temperatures on the original finite native patch.
use super::*;
use fs_couple::cohesive_contact::frictional_damage::{DamageFrictionConfig, DamageFrictionRegime};
use fs_couple::cohesive_contact::patch::frictional_damage::DamageFrictionPatchConfig;
use fs_couple::cohesive_contact::reclosure::{
    ReclosureActivation, ResidualCoulombSource,
    source::{RESIDUAL_KINETIC_FRICTION_PROPERTY, RESIDUAL_STATIC_FRICTION_PROPERTY,
        resolve_residual_coulomb},
};
use fs_couple::render::schedule::force::coupled::contact::cohesive::patch::{
    CohesivePatchModalError, CohesivePatchModalFrame, CohesivePatchPointFrame,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::CohesivePoint;
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Dims, Pressure, Temperature};

const AREAS: [f64; 3] = [0.2 * AREA, 0.3 * AREA, 0.5 * AREA];
const SCALES: [f64; 3] = [0.7, 1.0, 1.6];
const TEMPERATURES: [f64; 3] = [320.0, 300.0, 310.0];

fn card(label: &str, scale: f64, hot_mu: [f64; 2]) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (COHESIVE_STIFFNESS_PROPERTY, stiffness, scale * 1e9, scale * 5e8),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, scale * 2e5, scale * 8e4),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            scale * 100.0,
            scale * 50.0,
        ),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, scale * 2e9, scale * 3e9),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, 0.6, hot_mu[0]),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, 0.3, hot_mu[1]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, cold), (320.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured temperature-friction fixture, not measured data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: format!("{label}-{name}"),
            phase: "solid".into(),
            process: "original temperature-friction fixture".into(),
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
            environment: "temperature query".into(),
            history: "original bond and residual friction".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn sources(card: &InterfaceSystemCard) -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    (
        resolve_bilinear_cohesive(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            0.02,
        )
        .unwrap(),
        resolve_residual_coulomb(card, &query, MaterialPropertySelection::SingleClaimOnly).unwrap(),
    )
}

fn receiving(source: &ResolvedBilinearCohesive) -> [&MaterialStateId; 2] {
    [&source.surfaces()[0].material, &source.surfaces()[1].material]
}

fn damage_config() -> DamageFrictionConfig {
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

fn cards() -> [InterfaceSystemCard; 3] {
    core::array::from_fn(|p| card(&format!("site-{p}"), SCALES[p], [0.3, 0.1]))
}
fn build<'a>(sources: &[&'a ResolvedBilinearCohesive]) -> CohesivePatchModalSystem<'a> {
    let sites = (0..3).map(|p| CohesiveModalSite::from_source(
        sources[p], receiving(sources[p]), AREAS[p], bar_columns(LOCATIONS[p]),
    ).unwrap()).collect();
    CohesivePatchModalSystem::new(
        bars([0.003, -0.002, 0.001, 0.04]), sites, patch_config(), &CancelGate::new(),
    ).unwrap()
}
fn modes(owner: &CohesivePatchModalSystem<'_>) -> Vec<ModalAcousticState> {
    owner.components().iter().flat_map(|body| body.states().iter().copied()).collect()
}
#[derive(Debug, PartialEq)]
struct Mechanics {
    modes: Vec<ModalAcousticState>,
    body_storage: Vec<*const ModalAcousticState>,
    sample: u64,
    frame: Option<CohesivePatchModalFrame>,
    point_frames: Option<Vec<CohesivePatchPointFrame>>,
    compliance: Vec<f64>,
    compliance_storage: *const f64,
    columns: Vec<[Vec<f64>; 3]>,
    areas: Vec<f64>,
}
fn mechanics(owner: &CohesivePatchModalSystem<'_>) -> Mechanics {
    Mechanics {
        modes: modes(owner),
        body_storage: owner.components().iter().map(|body| body.states().as_ptr()).collect(),
        sample: owner.samples_rendered(),
        frame: owner.last_frame().copied(),
        point_frames: owner.last_point_frames().map(|frames| frames.to_vec()),
        compliance: owner.compliance_m_per_n().to_vec(),
        compliance_storage: owner.compliance_m_per_n().as_ptr(),
        columns: owner.sites().iter().map(|site| site.work_columns().clone()).collect(),
        areas: owner.sites().iter().map(|site| site.area_m2()).collect(),
    }
}
#[derive(Debug, PartialEq)]
struct Accepted<'a> {
    mechanics: Mechanics,
    points: Vec<CohesivePoint<'a>>,
    sources: Vec<*const ResolvedBilinearCohesive>,
    friction: Vec<Option<*const ResidualCoulombSource>>,
    occupancy: Vec<Option<DamageFrictionRegime>>,
}
fn accepted<'a>(owner: &CohesivePatchModalSystem<'a>) -> Accepted<'a> {
    Accepted {
        mechanics: mechanics(owner),
        points: owner.points().to_vec(),
        sources: owner.sites().iter().map(|site| site.source().unwrap() as *const _).collect(),
        friction: (0..3).map(|p| owner.damage_friction_source(p).map(|s| s as *const _)).collect(),
        occupancy: (0..3).map(|p| owner.damage_friction_occupancy(p)).collect(),
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}; tolerance={tolerance:e}");
}

// The original unconnected free-coordinate bars have exact free motion v*dt.
// Prescribe a known force/jump path through EVERY native cross-compliance term.
fn drive_for(
    owner: &CohesivePatchModalSystem<'_>,
    target: &[[f64; 3]],
    forces: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let state = modes(owner);
    let compliance = owner.compliance_m_per_n();
    let n = target.len();
    (0..n).map(|p| core::array::from_fn(|i| {
        let response = (0..n).map(|q| (0..3).map(|j| {
            compliance[(3 * p + i) * 3 * n + 3 * q + j] * forces[q][j]
        }).sum::<f64>()).sum::<f64>();
        let velocity = owner.sites()[p].work_columns()[i].iter().zip(&state)
            .map(|(column, mode)| column * mode.velocity_m_sqrt_kg_per_s).sum::<f64>();
        (target[p][i] + response - owner.points()[p].jump_m()[i])
            / owner.sample_period_s() - velocity
    })).collect()
}
fn seed(owner: &mut CohesivePatchModalSystem<'_>, target: &[[f64; 3]]) {
    let forces: Vec<_> = owner.points().iter().enumerate().map(|(p, point)| {
        point.trial_over_time_with_mean_traction(target[p], owner.sample_period_s(), 1e-8)
            .unwrap().1.map(|traction| AREAS[p] * traction)
    }).collect();
    let drive = drive_for(owner, target, &forces);
    owner.step(&[0.0; 8], &drive).unwrap();
    for (point, target) in owner.points().iter().zip(target) {
        for i in 0..3 { near(point.jump_m()[i], target[i], 2e-12); }
    }
}
fn enable<'a>(
    owner: &mut CohesivePatchModalSystem<'a>,
    sources: &[&'a ResolvedBilinearCohesive],
    friction: &[&'a ResidualCoulombSource],
) {
    owner.enable_frictional_damage(
        friction, &sources.iter().map(|source| receiving(source)).collect::<Vec<_>>(),
        &[ReclosureActivation::Sliding, ReclosureActivation::AtRest, ReclosureActivation::Sliding],
        DamageFrictionPatchConfig {
            max_points: 8, friction: damage_config(), newton_damping: 1e-14,
        }, &CancelGate::new(),
    ).unwrap();
}

#[test]
fn heterogeneous_temperature_accounts_and_natural_patch_redistribution_keep_original_mechanics() {
    let cards = cards();
    let pairs: Vec<_> = cards.iter().map(sources).collect();
    let cold: Vec<_> = pairs.iter().map(|pair| &pair.0).collect();
    let hot: Vec<_> = (0..3).map(|p| {
        cold[p].at_temperature(&cards[p], "T", TEMPERATURES[p]).unwrap()
    }).collect();
    let next: Vec<_> = hot.iter().collect();
    let target = [[-2e-4, 1.8e-4, 0.0], [-1.5e-4, 2.8e-4, 0.0], [2.5e-4, -8e-5, 3e-5]];
    let mut owner = build(&cold);
    let mut unheated = build(&cold);
    seed(&mut owner, &target);
    seed(&mut unheated, &target);
    let original = owner.points().to_vec();
    let mechanical = mechanics(&owner);
    let (report, driver_work) = owner.apply_temperature_successors_with_admission(
        &next, None, "T", &CancelGate::new(),
        |report, candidate| {
            assert_eq!(candidate.len(), 3);
            for p in 0..3 {
                assert!(core::ptr::eq(candidate[p].law(), next[p].law()));
                assert_eq!(candidate[p].jump_m(), original[p].jump_m());
            }
            Ok::<f64, CohesivePatchModalError>(report.parameter_work_j)
        },
    ).unwrap();
    assert_eq!(mechanics(&owner), mechanical);
    assert_eq!(driver_work, report.parameter_work_j);
    let mut expected_work = 0.0;
    let mut expected_release = 0.0;
    for p in 0..3 {
        let point = original[p];
        let g = point.jump_m();
        let radius2 = g[0].max(0.0).powi(2) + g[1] * g[1] + g[2] * g[2];
        let compression2 = g[0].min(0.0).powi(2);
        let fraction = (TEMPERATURES[p] - 280.0) / 40.0;
        let stiffness = SCALES[p] * (1e9 - 5e8 * fraction);
        let compression = SCALES[p] * (2e9 + 1e9 * fraction);
        let onset = next[p].law().onset_separation_m();
        let failure = next[p].law().failure_separation_m();
        let envelope = if radius2 == 0.0 { 0.0 } else {
            (failure / (failure - onset) * (1.0 - onset / radius2.sqrt())).clamp(0.0, 1.0)
        };
        let damage = point.damage().max(envelope);
        let parameter = AREAS[p] * (0.5 * (1.0 - point.damage())
            * (stiffness - SCALES[p] * 1e9) * radius2
            + 0.5 * (compression - SCALES[p] * 2e9) * compression2);
        let release = AREAS[p] * 0.5 * stiffness * (damage - point.damage()) * radius2;
        near(report.points[p].parameter_work_j, parameter, 2e-16);
        near(report.points[p].damage_dissipation_j, release, 2e-16);
        near(owner.points()[p].damage(), damage, 2e-12);
        near(owner.points()[p].dissipated_j_m2(),
            point.dissipated_j_m2() + release / AREAS[p], 2e-10);
        assert!(core::ptr::eq(owner.sites()[p].source().unwrap(), next[p]));
        assert!(core::ptr::eq(owner.sites()[p].law(), next[p].law()));
        assert_eq!(report.points[p].area_m2, AREAS[p]);
        assert!(report.points[p].energy_residual_j.abs() <= report.points[p].energy_tolerance_j);
        expected_work += parameter;
        expected_release += release;
    }
    near(report.parameter_work_j, expected_work, 4e-16);
    near(report.damage_dissipation_j, expected_release, 4e-16);
    assert!(report.damage_dissipation_j > 0.0);
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    // The tolerance belongs to the interface accounts, not body kinetic energy.
    assert!(report.energy_tolerance_j < 1.1e-12);

    let before_step = owner.points().to_vec();
    let old_modes = modes(&owner);
    owner.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
    unheated.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
    assert_ne!(modes(&owner), modes(&unheated));
    let observations = owner.last_point_frames().unwrap();
    let compliance = owner.compliance_m_per_n();
    let mut cross_response = 0.0_f64;
    for p in 0..3 {
        for i in 0..3 {
            let velocity = owner.sites()[p].work_columns()[i].iter().zip(&old_modes)
                .map(|(column, mode)| column * mode.velocity_m_sqrt_kg_per_s).sum::<f64>();
            let reaction = (0..3).map(|q| (0..3).map(|j| {
                compliance[(3 * p + i) * 9 + 3 * q + j] * observations[q].mean_force_n[j]
            }).sum::<f64>()).sum::<f64>();
            near(owner.points()[p].jump_m()[i],
                before_step[p].jump_m()[i] + owner.sample_period_s() * velocity - reaction, 2e-12);
            let other_sites = (0..3).filter(|q| *q != p).map(|q| {
                (0..3).map(|j| compliance[(3 * p + i) * 9 + 3 * q + j]
                    * observations[q].mean_force_n[j]).sum::<f64>()
            }).sum::<f64>();
            cross_response = cross_response.max(other_sites.abs());
        }
    }
    assert!(cross_response > 1e-12);
    assert!(observations.iter().zip(unheated.last_point_frames().unwrap())
        .any(|(a, b)| (a.mean_force_n[0] - b.mean_force_n[0]).abs() > 1e-5));
    let heated = owner.points().to_vec();
    let mechanical = mechanics(&owner);
    let cooling = owner.apply_temperature_successors(&cold, "T", &CancelGate::new()).unwrap();
    assert_eq!(mechanics(&owner), mechanical);
    for p in 0..3 {
        assert_eq!(owner.points()[p].jump_m(), heated[p].jump_m());
        assert!(owner.points()[p].damage() >= heated[p].damage());
        assert!(owner.points()[p].dissipated_j_m2() >= heated[p].dissipated_j_m2());
    }
    assert!(cooling.energy_residual_j.abs() <= cooling.energy_tolerance_j);
}

#[test]
fn paired_patch_temperature_preserves_occupancy_and_uses_each_new_friction_law() {
    let cards = cards();
    let pairs: Vec<_> = cards.iter().map(sources).collect();
    let cold: Vec<_> = pairs.iter().map(|pair| &pair.0).collect();
    let friction: Vec<_> = pairs.iter().map(|pair| &pair.1).collect();
    let hot: Vec<_> = (0..3).map(|p| cold[p]
        .at_temperature(&cards[p], "T", TEMPERATURES[p]).unwrap()).collect();
    let hot_friction: Vec<_> = (0..3).map(|p| friction[p]
        .at_temperature(&cards[p], "T", TEMPERATURES[p]).unwrap()).collect();
    let next: Vec<_> = hot.iter().collect();
    let next_friction: Vec<_> = hot_friction.iter().collect();
    let mut owner = build(&cold);
    seed(&mut owner, &[[-2e-4, 3e-4, 0.0], [-3e-4, -2.8e-4, 1e-5], [-1e-4, 3.5e-4, -5e-5]]);
    enable(&mut owner, &cold, &friction);
    let occupancy: Vec<_> = (0..3).map(|p| owner.damage_friction_occupancy(p)).collect();
    assert_eq!(occupancy, vec![Some(DamageFrictionRegime::Sliding),
        Some(DamageFrictionRegime::Sticking), Some(DamageFrictionRegime::Sliding)]);
    let mechanical = mechanics(&owner);
    owner.apply_temperature_and_friction_successors(
        &next, &next_friction, "T", &CancelGate::new(),
    ).unwrap();
    assert_eq!(mechanics(&owner), mechanical);
    assert_eq!((0..3).map(|p| owner.damage_friction_occupancy(p)).collect::<Vec<_>>(), occupancy);
    for p in 0..3 {
        assert!(core::ptr::eq(owner.damage_friction_source(p).unwrap(), next_friction[p]));
    }
    let old = owner.points().to_vec();
    let target: Vec<_> = old.iter().map(|point| {
        let mut jump = point.jump_m(); jump[1] += 1e-6; jump
    }).collect();
    let mut expected_loss = 0.0;
    let forces: Vec<_> = (0..3).map(|p| {
        let (_, mean, pressure) = old[p].trial_over_time_with_damage_weighted_compression(
            target[p], owner.sample_period_s(), 1e-8,
        ).unwrap();
        let coefficient = next_friction[p].law().zero_speed_kinetic_coefficient().unwrap();
        let extra = coefficient * AREAS[p] * pressure;
        expected_loss += extra * 1e-6;
        [AREAS[p] * mean[0], AREAS[p] * mean[1] + extra, AREAS[p] * mean[2]]
    }).collect();
    let drive = drive_for(&owner, &target, &forces);
    let frame = *owner.step(&[0.0; 8], &drive).unwrap();
    near(frame.friction_dissipation_j, expected_loss, 2e-12);
    for p in 0..3 {
        let response = owner.last_point_frames().unwrap()[p];
        let extra = response.damage_friction.unwrap();
        assert_eq!(extra.regime, DamageFrictionRegime::Sliding);
        near(extra.force_n[0], next_friction[p].law().zero_speed_kinetic_coefficient().unwrap()
            * extra.damaged_compression_force_n, 2e-10);
        for i in 0..3 { near(response.jump_m[i], target[p][i], 2e-12); }
    }
}

#[test]
fn patch_temperature_late_receiver_refusal_and_cancellation_leave_exact_replay() {
    let cards = cards();
    let pairs: Vec<_> = cards.iter().map(sources).collect();
    let cold: Vec<_> = pairs.iter().map(|pair| &pair.0).collect();
    let friction: Vec<_> = pairs.iter().map(|pair| &pair.1).collect();
    let hot: Vec<_> = (0..3).map(|p| cold[p]
        .at_temperature(&cards[p], "T", TEMPERATURES[p]).unwrap()).collect();
    let hot_friction: Vec<_> = (0..3).map(|p| friction[p]
        .at_temperature(&cards[p], "T", TEMPERATURES[p]).unwrap()).collect();
    let next: Vec<_> = hot.iter().collect();
    let next_friction: Vec<_> = hot_friction.iter().collect();
    let mut owner = build(&cold);
    let mut replay = build(&cold);
    let target = [[-2e-4, 3e-4, 0.0]; 3];
    for system in [&mut owner, &mut replay] {
        seed(system, &target);
        enable(system, &cold, &friction);
    }
    let before = accepted(&owner);
    assert!(owner.apply_temperature_successors(&next, "T", &CancelGate::new()).is_err());
    let mut wrong_site = next.clone();
    wrong_site[2] = next[0];
    assert!(owner.apply_temperature_and_friction_successors(
        &wrong_site, &next_friction, "T", &CancelGate::new(),
    ).is_err());
    let mut wrong_query = next_friction.clone();
    wrong_query[2] = friction[2];
    assert!(owner.apply_temperature_and_friction_successors(
        &next, &wrong_query, "T", &CancelGate::new(),
    ).is_err());
    assert!(owner.apply_temperature_and_friction_successors(
        &next[..2], &next_friction, "T", &CancelGate::new(),
    ).is_err());
    assert!(owner.apply_temperature_successors_with_admission::<(), CohesivePatchModalError>(
        &next, Some(&next_friction), "T", &CancelGate::new(),
        |report, candidate| {
            assert_eq!(report.points.len(), 3);
            assert_eq!(candidate.len(), 3);
            Err(CohesivePatchModalError::Contact(
                fs_couple::cohesive_contact::patch::CohesivePatchError::Invalid("receiver refused"),
            ))
        },
    ).is_err());
    let gate = CancelGate::new();
    assert!(owner.apply_temperature_successors_with_admission(
        &next, Some(&next_friction), "T", &gate,
        |report, _| {
            assert_eq!(report.points.len(), 3);
            gate.request();
            Ok::<(), CohesivePatchModalError>(())
        },
    ).is_err());
    assert_eq!(accepted(&owner), before);
    let result = owner.apply_temperature_and_friction_successors(
        &next, &next_friction, "T", &CancelGate::new(),
    ).unwrap();
    let reference = replay.apply_temperature_and_friction_successors(
        &next, &next_friction, "T", &CancelGate::new(),
    ).unwrap();
    assert_eq!(result, reference);
    assert_eq!(owner.points(), replay.points());
    assert_eq!(modes(&owner), modes(&replay));
    assert_eq!(owner.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap(),
        replay.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap());
    assert_eq!(modes(&owner), modes(&replay));
    assert_eq!(owner.last_point_frames(), replay.last_point_frames());
}
