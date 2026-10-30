//! G1/G3/G4: original finite-rotation bodies with damaged-area tangent disks.
use fs_couple::cohesive_contact::frictional_damage::DamageFrictionRegime;
use fs_couple::cohesive_contact::objective::{
    ObjectiveCohesiveConfig, ObjectiveCohesivePoint, source::SourceBoundObjectiveCohesivePoint,
};
use fs_couple::cohesive_contact::reclosure::{
    ReclosureActivation, ResidualCoulombSource, source::*,
};
use fs_couple::cohesive_contact::rigid::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::{CohesivePoint, source::*};
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::{Dims, Pressure};

const AREA: f64 = 1e-6;
const KC: f64 = 2e9;

fn sources(mu: [f64; 2], viscous: bool) -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let mut values = vec![
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
            KC,
        ),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, mu[0]),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, mu[1]),
    ];
    if viscous {
        values.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                2e3,
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                3e3,
            ),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured original rotating damaged-area fixture".into(),
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
            process: "original rigid body".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "original bonded pair".into(),
            third_body: None,
            environment: "fixed external state".into(),
            history: "original surface history".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let cohesive = if viscous {
        resolve_viscous_bilinear_cohesive(
            &card,
            &query,
            selection.clone(),
            CohesiveEnvelope::CommonMode,
            0.05,
            10.0,
        )
        .unwrap()
    } else {
        resolve_bilinear_cohesive(&card, &query, selection.clone(), 0.05).unwrap()
    };
    (
        cohesive,
        resolve_residual_coulomb(&card, &query, selection).unwrap(),
    )
}

fn config() -> RigidCohesiveConfig {
    RigidCohesiveConfig {
        maximum_sites: 64,
        maximum_iterations: 64,
        maximum_evaluations: 16384,
        maximum_material_evaluations: 100000,
        maximum_backtracks: 24,
        maximum_step_s: 0.02,
        maximum_site_force_n: 1000.0,
        maximum_total_force_n: 1e4,
        maximum_total_torque_n_m: 1000.0,
        maximum_jump_increment_m: 0.01,
        force_absolute_tolerance_n: 2e-10,
        torque_absolute_tolerance_n_m: 2e-12,
        relative_tolerance: 2e-10,
        energy_absolute_tolerance_j: 2e-10,
        energy_relative_tolerance: 2e-9,
        linear_momentum_absolute_tolerance_n_s: 2e-11,
        angular_momentum_absolute_tolerance_n_m_s: 2e-12,
        momentum_relative_tolerance: 2e-9,
        frame_tolerance: 1e-10,
        objective: ObjectiveCohesiveConfig {
            maximum_cayley_increment: 1.0,
            kinematic_tolerance_m: 2e-12,
            material_energy_tolerance_j_m2: 2e-7,
            spatial_work_tolerance_j: 2e-12,
        },
    }
}
fn friction_config() -> RigidDamageFrictionConfig {
    RigidDamageFrictionConfig {
        maximum_iterations: 100,
        slip_absolute_tolerance_m: 2e-12,
        slip_relative_tolerance: 1e-10,
        newton_damping: 1e-12,
    }
}
fn partial(source: &ResolvedBilinearCohesive) -> CohesivePoint<'_> {
    source
        .law()
        .bonded_point()
        .trial_over_time([4e-4, 0.0, 0.0], 1.0, 1e-8)
        .unwrap()
        .equilibrium
        .point
        .trial_over_time([-2e-4, 0.0, 0.0], 1.0, 1e-8)
        .unwrap()
        .equilibrium
        .point
}
fn failed(source: &ResolvedBilinearCohesive) -> CohesivePoint<'_> {
    source
        .law()
        .bonded_point()
        .trial_over_time([0.002, 0.0, 0.0], 1.0, 1e-8)
        .unwrap()
        .equilibrium
        .point
        .trial_over_time([-2e-4, 0.0, 0.0], 1.0, 1e-8)
        .unwrap()
        .equilibrium
        .point
}
fn v(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn norm(x: Vec3) -> f64 {
    x.x.hypot(x.y).hypot(x.z)
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}; tolerance {tol:e}");
}
fn close_vec(a: Vec3, b: Vec3, tol: f64) {
    close(norm(a.sub(b)), 0.0, tol);
}

#[allow(clippy::too_many_arguments)]
fn system<'a>(
    source: &'a ResolvedBilinearCohesive,
    friction: Option<&'a ResidualCoulombSource>,
    point: CohesivePoint<'a>,
    arms: &[[f64; 3]],
    rotation: UnitQuaternion,
    velocity: [Vec3; 2],
    omega: [Vec3; 2],
    dt: f64,
    c: RigidCohesiveConfig,
) -> RigidCohesiveSystem<'a> {
    let bodies = core::array::from_fn(|i| {
        let properties = MassProperties::from_principal_frame(
            [0.4, 0.7][i],
            Vec3::ZERO,
            Vec3::new(0.002, 0.003, 0.004),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.37).unwrap(),
        )
        .unwrap();
        let position = if i == 0 {
            Vec3::ZERO
        } else {
            rotation.rotate_body_to_world(v(point.jump_m()))
        };
        let state = RigidBodyState::new(
            Pose::new(position, rotation).unwrap(),
            rotation
                .rotate_body_to_world(velocity[i])
                .scale(properties.mass()),
            properties.angular_momentum_body_checked(omega[i]).unwrap(),
        )
        .unwrap();
        RigidCohesiveBody::new(
            source.surfaces()[i].material.clone(),
            properties,
            state,
            RigidCohesiveBodyLimits {
                maximum_linear_speed_m_s: 100.0,
                maximum_angular_speed_rad_s: 1e4,
                maximum_kinetic_energy_j: 1e6,
                maximum_step_displacement_m: 0.1,
            },
        )
        .unwrap()
    });
    let sites = arms
        .iter()
        .map(|arm| {
            let point =
                ObjectiveCohesivePoint::new(point, rotation, Vec3::ZERO, AREA / arms.len() as f64)
                    .unwrap();
            let point = SourceBoundObjectiveCohesivePoint::new(
                point,
                source,
                [
                    &source.surfaces()[0].material,
                    &source.surfaces()[1].material,
                ],
            )
            .unwrap();
            RigidCohesiveSite::new(point, [v(*arm); 2], UnitQuaternion::IDENTITY).unwrap()
        })
        .collect();
    let mut system = RigidCohesiveSystem::new(
        bodies,
        sites,
        Gravity::new(Vec3::ZERO).unwrap(),
        dt,
        c,
        &CancelGate::new(),
    )
    .unwrap();
    if let Some(friction) = friction {
        system
            .enable_damaged_area_friction(
                vec![
                    RigidDamageFrictionSite {
                        source: friction,
                        activation: ReclosureActivation::AtRest,
                    };
                    arms.len()
                ],
                friction_config(),
                &CancelGate::new(),
            )
            .unwrap();
    }
    system
}

fn zero() -> [Vec3; 2] {
    [Vec3::ZERO; 2]
}
fn balanced_load(force: Vec3, middle_separation: Vec3) -> [RigidCohesiveLoad; 2] {
    [
        RigidCohesiveLoad {
            force_world_n: force.scale(-1.0),
            moment_world_n_m: middle_separation.cross(force).scale(-1.0),
        },
        RigidCohesiveLoad {
            force_world_n: force,
            moment_world_n_m: Vec3::ZERO,
        },
    ]
}

#[test]
fn g1_static_disk_holds_both_tangents_on_original_fully_failed_bodies() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = failed(&source);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let before = owner.bodies().clone();
    let force = Vec3::new(-0.4, 0.12, -0.16); // |T|=.2 < static .24; exceeds kinetic .12.
    let frame = *owner.step(balanced_load(force, v(point.jump_m()))).unwrap();
    assert_eq!(frame.friction_dissipation_j, 0.0);
    let site = owner.site_frames()[0].unwrap().friction.unwrap();
    assert_eq!(site.regime, DamageFrictionRegime::Sticking);
    close(site.local_force_n[0], 0.12, 2e-9);
    close(site.local_force_n[1], -0.16, 2e-9);
    for i in 0..2 {
        close_vec(
            owner.bodies()[i].state().pose().position_world(),
            before[i].state().pose().position_world(),
            2e-12,
        );
        close_vec(
            owner.bodies()[i].state().angular_momentum_body(),
            Vec3::ZERO,
            2e-12,
        );
    }
    assert!(owner.sites()[0].point().numerical().material().failed());
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}

#[test]
fn g1_vector_sliding_uses_only_damaged_compression_and_actual_objective_work() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = partial(&source);
    close(point.damage(), 0.625, 1e-14);
    let dt = 0.01;
    let delta = Vec3::new(0.0, 6e-5, 8e-5);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let expected = Vec3::new(-0.4, 0.05625, 0.075);
    let frame = *owner
        .step(balanced_load(
            expected,
            v(point.jump_m()).add(delta.scale(0.5)),
        ))
        .unwrap();
    let site = owner.site_frames()[0].unwrap();
    let f = site.friction.unwrap();
    assert_eq!(f.regime, DamageFrictionRegime::Sliding);
    close_vec(v(site.jump_m), v(point.jump_m()).add(delta), 1e-10);
    close(f.damaged_compression_force_n, 0.25, 1e-9);
    close(f.local_force_n[0], 0.045, 1e-9);
    close(f.local_force_n[1], 0.060, 1e-9);
    close(f.friction_dissipation_j, 7.5e-6, 1e-11);
    close(f.spatial_work_j, f.friction_dissipation_j, 2e-12);
    close(
        frame.friction_dissipation_j,
        f.friction_dissipation_j,
        1e-14,
    );
    close_vec(frame.contact_loads[0].force_world_n, expected, 2e-8);
    assert_eq!(frame.fracture_dissipation_j, 0.0);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

#[test]
fn g3_coincident_patch_partition_preserves_original_body_solution_and_total_work() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = partial(&source);
    let delta = Vec3::new(0.0, 6e-5, 8e-5);
    let load = balanced_load(
        Vec3::new(-0.4, 0.05625, 0.075),
        v(point.jump_m()).add(delta.scale(0.5)),
    );
    let mut one = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(100.0)],
        zero(),
        0.01,
        config(),
    );
    let mut split = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]; 3],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(100.0)],
        zero(),
        0.01,
        config(),
    );
    let a = *one.step(load).unwrap();
    let b = *split.step(load).unwrap();
    close(a.friction_dissipation_j, b.friction_dissipation_j, 2e-11);
    for i in 0..2 {
        close_vec(
            one.bodies()[i].state().pose().position_world(),
            split.bodies()[i].state().pose().position_world(),
            2e-10,
        );
        close_vec(
            one.bodies()[i].state().linear_momentum_world(),
            split.bodies()[i].state().linear_momentum_world(),
            2e-9,
        );
    }
    // The same redundant sites also admit static holding without physical springs.
    let mut stuck = system(
        &source,
        Some(&friction),
        failed(&source),
        &[[0.0; 3]; 3],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let frame = *stuck
        .step(balanced_load(
            Vec3::new(-0.4, 0.12, -0.16),
            Vec3::new(-2e-4, 0.0, 0.0),
        ))
        .unwrap();
    assert_eq!(frame.friction_dissipation_j, 0.0);
    assert!(
        stuck
            .site_frames()
            .iter()
            .all(|s| s.unwrap().friction.unwrap().regime == DamageFrictionRegime::Sticking)
    );
}

#[test]
fn g3_rotating_spread_patch_covaries_and_keeps_native_momentum_and_losses() {
    let (source, friction) = sources([0.6, 0.3], true);
    let arms = [
        [0.0, -0.002, 0.001],
        [0.0, 0.001, -0.002],
        [0.0, 0.002, 0.002],
    ];
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), 0.8).unwrap();
    let velocity = [Vec3::new(0.0, -0.025, 0.015), Vec3::new(0.0, 0.020, -0.010)];
    let omega = [Vec3::new(0.3, -0.2, 0.4), Vec3::new(-0.1, 0.15, 0.2)];
    let mut a = system(
        &source,
        Some(&friction),
        partial(&source),
        &arms,
        UnitQuaternion::IDENTITY,
        velocity,
        omega,
        1e-4,
        config(),
    );
    let mut b = system(
        &source,
        Some(&friction),
        partial(&source),
        &arms,
        rotation,
        velocity,
        omega,
        1e-4,
        config(),
    );
    let mut friction_loss = 0.0;
    let mut viscosity_loss = 0.0;
    for _ in 0..12 {
        let x = *a.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        let y = *b.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        close(x.friction_dissipation_j, y.friction_dissipation_j, 2e-10);
        friction_loss += x.friction_dissipation_j;
        viscosity_loss += x.viscous_dissipation_j;
        assert!(x.energy_residual_j.abs() <= x.energy_tolerance_j);
        assert!(norm(x.angular_momentum_residual_n_m_s) <= x.angular_momentum_tolerance_n_m_s);
        for i in 0..2 {
            close_vec(
                rotation.rotate_body_to_world(a.bodies()[i].state().pose().position_world()),
                b.bodies()[i].state().pose().position_world(),
                1e-9,
            );
            close_vec(
                rotation.rotate_body_to_world(a.bodies()[i].state().linear_momentum_world()),
                b.bodies()[i].state().linear_momentum_world(),
                1e-8,
            );
        }
    }
    assert!(friction_loss > 0.0 && viscosity_loss > 0.0);
    assert!(a.bodies()[0].properties().inertia_tensor_body().m[1].abs() > 0.0);
}

#[test]
fn g4_source_refusal_and_late_energy_failure_preserve_all_original_state_and_retry() {
    let (source, friction) = sources([0.6, 0.3], false);
    let (_, foreign) = sources([0.7, 0.3], false);
    let point = partial(&source);
    let mut a = system(
        &source,
        None,
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let before = a.bodies().clone();
    assert!(
        a.enable_damaged_area_friction(
            vec![RigidDamageFrictionSite {
                source: &foreign,
                activation: ReclosureActivation::AtRest
            }],
            friction_config(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert!(a.damaged_area_friction_sources().is_none());
    assert_eq!(a.bodies(), &before);
    a.enable_damaged_area_friction(
        vec![RigidDamageFrictionSite {
            source: &friction,
            activation: ReclosureActivation::AtRest,
        }],
        friction_config(),
        &CancelGate::new(),
    )
    .unwrap();
    let original = a.damaged_area_friction_states().unwrap().to_vec();
    let gate = CancelGate::new();
    gate.request();
    let load = balanced_load(Vec3::new(-0.4, 0.03, -0.02), v(point.jump_m()));
    assert_eq!(
        a.step_under_gate(load, &gate),
        Err(RigidCohesiveError::Cancelled)
    );
    assert_eq!(a.steps(), 0);
    assert_eq!(a.time_s(), 0.0);
    assert!(a.last_frame().is_none());
    assert_eq!(a.damaged_area_friction_states().unwrap(), original);
    let mut b = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    assert_eq!(a.step(load).unwrap(), b.step(load).unwrap());
    assert_eq!(a.bodies(), b.bodies());
    assert_eq!(a.site_frames(), b.site_frames());
    let mut strict = config();
    strict.energy_absolute_tolerance_j = 1e-30;
    strict.energy_relative_tolerance = 1e-30;
    let mut refused = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::new(0.0, -0.02, 0.01), Vec3::new(0.0, 0.03, -0.01)],
        zero(),
        1e-4,
        strict,
    );
    let bodies = refused.bodies().clone();
    let states = refused.damaged_area_friction_states().unwrap().to_vec();
    let error = refused.step([RigidCohesiveLoad::ZERO; 2]).unwrap_err();
    assert!(
        matches!(
            error,
            RigidCohesiveError::Balance {
                what: "actual rigid/cohesive energy",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(refused.bodies(), &bodies);
    assert_eq!(refused.damaged_area_friction_states().unwrap(), states);
    assert_eq!(refused.steps(), 0);
    assert!(refused.site_frames().iter().all(Option::is_none));
}

#[test]
fn g1_progressive_fracture_integrates_damaged_area_and_keeps_three_losses_separate() {
    let (source, friction) = sources([0.6, 0.3], true);
    let point = partial(&source);
    let end = [-2e-4, 7.2e-4, 9.6e-4];
    let delta = v(end).sub(v(point.jump_m()));
    let dt = 0.01;
    let (material, mean) = point
        .trial_over_time_with_mean_traction(end, dt, 2e-7)
        .unwrap();
    // Independent composite Simpson integral of the actual source damage law;
    // no call to the production damaged-compression path integral.
    let panels = 12000;
    let mut damage_sum = 0.0;
    for i in 0..=panels {
        let s = i as f64 / panels as f64;
        let jump = core::array::from_fn(|j| point.jump_m()[j] + s * (end[j] - point.jump_m()[j]));
        let damage = point
            .trial_over_time(jump, dt, 2e-7)
            .unwrap()
            .equilibrium
            .point
            .damage();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        damage_sum += weight * damage;
    }
    let pressure = AREA * KC * 2e-4 * damage_sum / (3.0 * panels as f64);
    let cap = 0.3 * pressure;
    let expected = v(mean.map(|x| AREA * x)).add(Vec3::new(0.0, 0.6 * cap, 0.8 * cap));
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let frame = *owner
        .step(balanced_load(
            expected,
            v(point.jump_m()).add(delta.scale(0.5)),
        ))
        .unwrap();
    let site = owner.site_frames()[0].unwrap();
    let f = site.friction.unwrap();
    assert!(owner.sites()[0].point().numerical().material().failed());
    close_vec(v(site.jump_m), v(end), 2e-10);
    close(f.damaged_compression_force_n, pressure, 2e-8);
    close(f.friction_dissipation_j, cap * 0.0012, 2e-11);
    close(
        frame.fracture_dissipation_j,
        AREA * material.equilibrium.energy.dissipated,
        2e-11,
    );
    close(
        frame.viscous_dissipation_j,
        AREA * material.viscous_dissipated_j_m2,
        2e-11,
    );
    assert!(
        frame.fracture_dissipation_j > 0.0
            && frame.viscous_dissipation_j > 0.0
            && frame.friction_dissipation_j > 0.0
    );
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}

#[test]
fn g1_open_fracture_and_zero_coefficient_recover_original_native_mechanics() {
    for opened in [false, true] {
        let (source, friction) = sources(if opened { [0.6, 0.3] } else { [0.0, 0.0] }, false);
        let point = if opened {
            failed(&source)
                .trial_over_time([2e-4, 0.0, 0.0], 1.0, 1e-8)
                .unwrap()
                .equilibrium
                .point
        } else {
            failed(&source)
        };
        let velocity = [Vec3::new(0.0, -0.02, 0.01), Vec3::new(0.0, 0.03, -0.02)];
        let omega = [Vec3::new(0.1, 0.2, 0.1), Vec3::new(0.2, -0.1, 0.1)];
        let mut a = system(
            &source,
            None,
            point,
            &[[0.0; 3]],
            UnitQuaternion::IDENTITY,
            velocity,
            omega,
            1e-4,
            config(),
        );
        let mut b = system(
            &source,
            Some(&friction),
            point,
            &[[0.0; 3]],
            UnitQuaternion::IDENTITY,
            velocity,
            omega,
            1e-4,
            config(),
        );
        let x = *a.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        let y = *b.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_eq!(y.friction_dissipation_j, 0.0);
        close(x.mechanical_energy_j, y.mechanical_energy_j, 2e-11);
        for i in 0..2 {
            close_vec(
                a.bodies()[i].state().linear_momentum_world(),
                b.bodies()[i].state().linear_momentum_world(),
                2e-10,
            );
        }
        assert_eq!(b.sites()[0].point().numerical().material().damage(), 1.0);
    }
}

#[test]
fn g4_unresolved_midpoint_capture_refuses_reflection_and_keeps_original_motion() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = failed(&source);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, Vec3::new(0.0, 3e-4, 0.0)],
        zero(),
        0.02,
        config(),
    );
    let bodies = owner.bodies().clone();
    let states = owner.damaged_area_friction_states().unwrap().to_vec();
    let error = owner
        .step(balanced_load(Vec3::new(-0.4, 0.0, 0.0), v(point.jump_m())))
        .unwrap_err();
    assert_eq!(
        error,
        RigidCohesiveError::Invalid(
            "unresolved friction stop reverses endpoint slip; subdivide the step"
        )
    );
    assert_eq!(owner.bodies(), &bodies);
    assert_eq!(owner.damaged_area_friction_states().unwrap(), states);
    assert_eq!(owner.steps(), 0);
    assert!(owner.site_frames().iter().all(Option::is_none));
}
