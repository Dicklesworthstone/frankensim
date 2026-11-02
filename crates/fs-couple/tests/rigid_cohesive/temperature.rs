//! Prescribed source temperatures must change actual original body dynamics.
use super::*;
use fs_couple::cohesive_contact::reclosure::{
    ReclosureActivation,
    source::{
        RESIDUAL_KINETIC_FRICTION_PROPERTY, RESIDUAL_STATIC_FRICTION_PROPERTY,
        resolve_residual_coulomb,
    },
};
use fs_qty::{Dims, Temperature};

fn card(history: &str) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (COHESIVE_NORMAL_STIFFNESS_PROPERTY, stiffness, 1e9, 5e8),
        (COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY, stiffness, 4e8, 3e8),
        (COHESIVE_NORMAL_STRENGTH_PROPERTY, Pressure::DIMS, 2e5, 8e4),
        (
            COHESIVE_TANGENTIAL_STRENGTH_PROPERTY,
            Pressure::DIMS,
            1.6e5,
            8e4,
        ),
        (
            COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            50.0,
        ),
        (
            COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            400.0,
            200.0,
        ),
        (COHESIVE_BK_EXPONENT_PROPERTY, Dims::NONE, 2.0, 1.5),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, 2e9, 3e9),
        (
            COHESIVE_NORMAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            2e3,
            4e3,
        ),
        (
            COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            3e3,
            6e3,
        ),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, 0.5, 0.5),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, 0.3, 0.3),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, cold), (340.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "manufactured thermal rigid cohesive fixture; not measured adhesive data"
                            .into(),
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
            process: "original native thermal fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "original finite bond".into(),
            third_body: None,
            environment: "explicit prescribed source temperature".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(card: &InterfaceSystemCard, viscous: bool) -> ResolvedBilinearCohesive {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    if viscous {
        resolve_viscous_bilinear_cohesive(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            CohesiveEnvelope::FrozenOnsetMixedMode,
            0.05,
            10.0,
        )
        .unwrap()
    } else {
        resolve_mixed_mode_bilinear_cohesive(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            0.05,
        )
        .unwrap()
    }
}

fn make<'a>(
    source: &'a ResolvedBilinearCohesive,
    jump: [f64; 3],
    cfg: RigidCohesiveConfig,
    restricted: bool,
) -> RigidCohesiveSystem<'a> {
    let mut body = bodies(
        source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [Vec3::ZERO; 2],
        [Vec3::ZERO; 2],
    );
    if restricted {
        let original = &body[0];
        let limits = RigidCohesiveBodyLimits {
            maximum_linear_speed_m_s: 1e-3,
            ..original.limits()
        };
        body[0] = RigidCohesiveBody::new(
            original.material().clone(),
            original.properties(),
            original.state(),
            limits,
        )
        .unwrap();
    }
    let point = source
        .law()
        .bonded_point()
        .trial_over_time(jump, 0.02, 1e-7)
        .unwrap()
        .equilibrium
        .point;
    let sites = [0.7e-6, 1.3e-6]
        .into_iter()
        .map(|area| site(source, &body, [Vec3::ZERO; 2], area, point))
        .collect();
    RigidCohesiveSystem::new(body, sites, Gravity::ZERO, 2e-4, cfg, &CancelGate::new()).unwrap()
}

#[test]
fn independently_heated_sites_drive_the_original_two_body_midpoint_oscillator() {
    let card = card("original bond");
    let cold = resolve(&card, false);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let warm = cold.at_temperature(&card, "T", 310.0).unwrap();
    let g0 = 3e-5;
    let mut system = make(&cold, [g0, 0.0, 0.0], config(), false);
    let initial_energy = system.total_energy_j().unwrap();
    let original_mass = system.bodies().each_ref().map(|b| b.properties());
    let h = system.step_duration_s();
    let stiffness = 0.7e-6 * 5e8 + 1.3e-6 * 7.5e8;
    let omega2 = stiffness * (1.0 / 0.4 + 1.0 / 0.7);
    let g1 = (1.0 - 0.25 * h * h * omega2) * g0 / (1.0 + 0.25 * h * h * omega2);
    let v1 = -h * omega2 * 0.5 * (g0 + g1);
    let expected_work = 0.5 * (stiffness - 2000.0) * g0 * g0;
    let frame = *system
        .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &warm], "T")
        .unwrap();
    assert_frame(&frame);
    close(frame.parameter_work_j, expected_work, 2e-15);
    close(frame.thermal_stored_change_j, expected_work, 2e-15);
    assert_eq!(frame.thermal_damage_dissipation_j, 0.0);
    assert_eq!(frame.fracture_dissipation_j, 0.0);
    assert_eq!(frame.external_work_j, 0.0);
    close(
        system.total_energy_j().unwrap(),
        initial_energy + expected_work,
        2e-14,
    );
    for (index, source) in [&hot, &warm].into_iter().enumerate() {
        assert!(core::ptr::eq(
            system.sites()[index].point().source(),
            source
        ));
        close(
            system.sites()[index]
                .point()
                .numerical()
                .material()
                .jump_m()[0],
            g1,
            2e-15,
        );
        let site_frame = system.site_frames()[index].unwrap();
        close(
            site_frame.parameter_work_j,
            0.5 * [0.7e-6, 1.3e-6][index] * ([5e8, 7.5e8][index] - 1e9) * g0 * g0,
            2e-15,
        );
    }
    let actual_v = system.bodies()[1]
        .state()
        .center_of_mass_velocity_world(original_mass[1])
        .unwrap()
        .sub(
            system.bodies()[0]
                .state()
                .center_of_mass_velocity_world(original_mass[0])
                .unwrap(),
        );
    close(actual_v.x, v1, 2e-12);
    assert_eq!(frame.material_evaluations, 2 * frame.evaluations + 2);
    // Accepted hot/warm coefficients persist, but parameter work is not replayed.
    let g2 = ((1.0 - 0.25 * h * h * omega2) * g1 + h * v1) / (1.0 + 0.25 * h * h * omega2);
    let next = *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    close(
        system.sites()[0].point().numerical().material().jump_m()[0],
        g2,
        2e-15,
    );
    assert_eq!(next.parameter_work_j, 0.0);
    assert_eq!(next.thermal_damage_dissipation_j, 0.0);
    close(
        system.total_energy_j().unwrap(),
        initial_energy + expected_work,
        2e-14,
    );
}

#[test]
fn finite_rotating_mixed_patch_preserves_sources_momenta_and_complete_thermal_energy() {
    let card = card("original bond");
    let cold = resolve(&card, true);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let warm = cold.at_temperature(&card, "T", 310.0).unwrap();
    let build = |rotation, translation| {
        let body = bodies(
            &cold,
            rotation,
            translation,
            [Vec3::new(-0.02, 0.01, -0.005), Vec3::new(0.03, -0.02, 0.01)],
            [Vec3::new(1.5, -2.0, 0.75), Vec3::new(-0.5, 1.0, -0.25)],
        );
        let history = cold
            .law()
            .bonded_point()
            .trial_over_time([4e-4, 5e-4, -2e-4], 0.02, 1e-7)
            .unwrap()
            .equilibrium
            .point;
        let sites = [
            Vec3::new(0.0, 0.004, 0.003),
            Vec3::new(0.0, -0.004, 0.003),
            Vec3::new(0.0, 0.0, -0.004),
        ]
        .into_iter()
        .zip([0.6e-6, 0.8e-6, 0.6e-6])
        .map(|(arm, area)| site(&cold, &body, [arm; 2], area, history))
        .collect();
        RigidCohesiveSystem::new(
            body,
            sites,
            Gravity::ZERO,
            2e-4,
            config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 1.2).unwrap();
    let translation = Vec3::new(0.1, -0.2, 0.3);
    let mut system = build(UnitQuaternion::IDENTITY, Vec3::ZERO);
    let mut rotated = build(rotation, translation);
    let mut fixed = build(UnitQuaternion::IDENTITY, Vec3::ZERO);
    let before = system.sites().to_vec();
    let initial_energy = system.total_energy_j().unwrap();
    let frame = *system
        .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &warm, &hot], "T")
        .unwrap();
    let covariance = *rotated
        .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &warm, &hot], "T")
        .unwrap();
    fixed.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    assert_frame(&frame);
    assert_frame(&covariance);
    assert!(frame.thermal_damage_dissipation_j > 0.0);
    assert!(frame.fracture_dissipation_j >= frame.thermal_damage_dissipation_j);
    assert!(frame.viscous_dissipation_j > 0.0);
    close(
        system.total_energy_j().unwrap() - initial_energy
            + frame.fracture_dissipation_j
            + frame.viscous_dissipation_j
            - frame.parameter_work_j,
        0.0,
        3e-11,
    );
    let (mut parameter, mut fracture) = (0.0, 0.0);
    for (index, source) in [&hot, &warm, &hot].into_iter().enumerate() {
        let old = before[index].point().numerical().material();
        let next = system.sites()[index].point().numerical().material();
        let area = before[index].point().numerical().area_m2();
        let g = old.jump_m();
        let p = source.law().mixed_mode_parameters().unwrap();
        parameter += 0.5
            * area
            * (1.0 - old.damage())
            * ((p.normal_stiffness_pa_m - 1e9) * g[0] * g[0]
                + (p.tangential_stiffness_pa_m - 4e8) * (g[1] * g[1] + g[2] * g[2]));
        fracture += area * (next.dissipated_j_m2() - old.dissipated_j_m2());
        assert!(core::ptr::eq(
            system.sites()[index].point().source(),
            source
        ));
        assert!(next.damage() >= old.damage());
        assert_eq!(
            next.mixed_mode_history().unwrap().onset_direction,
            old.mixed_mode_history().unwrap().onset_direction
        );
        close(
            system.site_frames()[index].unwrap().fracture_dissipation_j,
            area * (next.dissipated_j_m2() - old.dissipated_j_m2()),
            2e-14,
        );
    }
    close(frame.parameter_work_j, parameter, 2e-14);
    close(frame.fracture_dissipation_j, fracture, 2e-14);
    close(frame.parameter_work_j, covariance.parameter_work_j, 2e-14);
    close(
        frame.fracture_dissipation_j,
        covariance.fracture_dissipation_j,
        2e-11,
    );
    for index in 0..2 {
        let a = system.bodies()[index].state();
        let b = rotated.bodies()[index].state();
        close_vec(
            b.pose().position_world(),
            rotation
                .rotate_body_to_world(a.pose().position_world())
                .add(translation),
            3e-11,
        );
        close_vec(
            b.linear_momentum_world(),
            rotation.rotate_body_to_world(a.linear_momentum_world()),
            3e-11,
        );
        close_vec(b.angular_momentum_body(), a.angular_momentum_body(), 3e-11);
    }
    assert!(
        norm(
            system.bodies()[0]
                .state()
                .linear_momentum_world()
                .sub(fixed.bodies()[0].state().linear_momentum_world())
        ) > 1e-9
    );
}

#[test]
fn partial_source_staging_late_body_refusal_and_cancellation_leave_identical_retries() {
    let card = card("original bond");
    let cold = resolve(&card, false);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let foreign_card = self::card("different original preparation");
    let foreign = resolve(&foreign_card, false);
    let mut system = make(&cold, [4e-4, 5e-4, -2e-4], config(), true);
    let mut control = make(&cold, [4e-4, 5e-4, -2e-4], config(), true);
    system.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    control.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    let before = system.bodies().clone();
    let sites = system.sites().to_vec();
    let observed = system.site_frames().to_vec();
    let frame = *system.last_frame().unwrap();
    assert!(
        system
            .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &foreign], "T")
            .is_err()
    );
    assert!(
        system
            .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot], "T")
            .is_err()
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        system.step_with_temperature_under_gate(
            [RigidCohesiveLoad::ZERO; 2],
            &[&hot, &hot],
            "T",
            &cancelled
        ),
        Err(RigidCohesiveError::Cancelled)
    ));
    let excessive = [
        RigidCohesiveLoad {
            force_world_n: Vec3::new(100.0, 0.0, 0.0),
            ..RigidCohesiveLoad::ZERO
        },
        RigidCohesiveLoad::ZERO,
    ];
    let refusal = system.step_with_temperature(excessive, &[&hot, &hot], "T");
    assert!(
        matches!(
            refusal,
            Err(RigidCohesiveError::Range {
                what: "body linear speed",
                ..
            })
        ),
        "{refusal:?}"
    );
    assert_eq!(system.bodies(), &before);
    assert_eq!(system.sites(), sites);
    assert_eq!(system.site_frames(), observed);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.steps(), frame.step);
    assert_eq!(system.time_s(), frame.time_s);
    // A normal step must ignore all failed thermal scratch, including sources.
    assert_eq!(
        *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap(),
        *control.step([RigidCohesiveLoad::ZERO; 2]).unwrap()
    );
    assert_eq!(system.bodies(), control.bodies());
    assert_eq!(system.sites(), control.sites());
    let retry = *system
        .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &hot], "T")
        .unwrap();
    let expected = *control
        .step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &hot], "T")
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(system.bodies(), control.bodies());
    assert_eq!(system.sites(), control.sites());
    assert_eq!(system.site_frames(), control.site_frames());
}

#[test]
fn thermal_transfers_consume_original_material_budget_and_do_not_freeze_friction_sources() {
    let card = card("original bond");
    let cold = resolve(&card, false);
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    let cfg = RigidCohesiveConfig {
        maximum_material_evaluations: 2,
        ..config()
    };
    let mut system = make(&cold, [0.0; 3], cfg, false);
    let mut control = make(&cold, [0.0; 3], cfg, false);
    let before = system.bodies().clone();
    assert!(matches!(
        system.step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &hot], "T"),
        Err(RigidCohesiveError::Budget(
            "original cohesive material evaluations"
        ))
    ));
    assert_eq!(system.bodies(), &before);
    assert_eq!(system.steps(), 0);
    assert_eq!(
        *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap(),
        *control.step([RigidCohesiveLoad::ZERO; 2]).unwrap()
    );
    let friction = resolve_residual_coulomb(
        &card,
        &QueryPoint::new().with("T", 280.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut system = make(&cold, [0.0; 3], config(), false);
    system
        .enable_damaged_area_friction(
            vec![
                RigidDamageFrictionSite {
                    source: &friction,
                    activation: ReclosureActivation::AtRest
                };
                2
            ],
            RigidDamageFrictionConfig {
                maximum_iterations: 64,
                slip_absolute_tolerance_m: 1e-12,
                slip_relative_tolerance: 1e-10,
                newton_damping: 1e-10,
            },
            &CancelGate::new(),
        )
        .unwrap();
    let before = system.bodies().clone();
    let states = system.damaged_area_friction_states().unwrap().to_vec();
    assert!(matches!(
        system.step_with_temperature([RigidCohesiveLoad::ZERO; 2], &[&hot, &hot], "T"),
        Err(RigidCohesiveError::Invalid(
            "damaged-area friction requires paired temperature successors"
        ))
    ));
    assert_eq!(system.bodies(), &before);
    assert_eq!(system.damaged_area_friction_states().unwrap(), states);
    assert_eq!(system.steps(), 0);
}
