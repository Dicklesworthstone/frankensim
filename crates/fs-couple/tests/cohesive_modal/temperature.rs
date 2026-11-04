//! Source-temperature changes act on the original native interface transaction.
use super::*;
use fs_couple::cohesive_contact::frictional_damage::{DamageFrictionConfig, DamageFrictionRegime};
use fs_couple::cohesive_contact::reclosure::{
    ReclosureActivation, ReclosureConfig, ReclosureRegime, ResidualCoulombSource,
    source::{
        RESIDUAL_KINETIC_FRICTION_PROPERTY, RESIDUAL_STATIC_FRICTION_PROPERTY,
        resolve_residual_coulomb,
    },
};
use fs_couple::render::schedule::force::coupled::contact::cohesive::CohesiveModalFrame;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Dims, Pressure, Temperature};

fn card(hot_mu: [f64; 2]) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (COHESIVE_STIFFNESS_PROPERTY, stiffness, 1e9, 5e8),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5, 8e4),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            50.0,
        ),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, 2e9, 3e9),
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
            chemistry: name.into(),
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

fn reclosure_config() -> ReclosureConfig {
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

fn bound(
    source: &ResolvedBilinearCohesive,
    relative: [f64; 3],
    ceiling: f64,
) -> CohesiveModalSystem<'_> {
    CohesiveModalSystem::from_source(
        network(relative, ceiling, false),
        source,
        receiving(source),
        AREA,
        columns(),
        config(),
        &CancelGate::new(),
    )
    .unwrap()
}

// The original network consists of free coordinates. Supply the prescribed
// increment that gives this physical endpoint under the stated mean force.
fn drive_with_force(
    system: &mut CohesiveModalSystem<'_>,
    target: [f64; 3],
    force: [f64; 3],
) -> CohesiveModalFrame {
    let old = system.point().jump_m();
    let dt = system.sample_period_s();
    let velocity = relative_velocity(system);
    let compliance = system.compliance_m_per_n();
    let drive = core::array::from_fn(|i| {
        (target[i] - old[i] - dt * velocity[i]
            + (0..3).map(|j| compliance[i][j] * force[j]).sum::<f64>())
            / dt
    });
    *system.step(&[0.0; 6], drive).unwrap()
}

fn drive_to(system: &mut CohesiveModalSystem<'_>, target: [f64; 3]) -> CohesiveModalFrame {
    let (_, mean) = system
        .point()
        .trial_over_time_with_mean_traction(target, system.sample_period_s(), 1e-8)
        .unwrap();
    drive_with_force(system, target, mean.map(|traction| AREA * traction))
}

#[test]
fn held_temperature_preserves_native_motion_and_next_step_uses_the_hot_bond() {
    let card = card([0.3, 0.1]);
    let (cold, _) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let mut system = bound(&cold, [0.02, 0.005, -0.003], 1.0);
    drive_to(&mut system, [4e-4, 0.0, 0.0]);
    let old = system.point();
    close(old.damage(), 0.625, 2e-10);
    assert!(relative_velocity(&system).iter().any(|v| v.abs() > 1e-6));
    let native = states(&system);
    let sample = system.samples_rendered();
    let last = system.last_frame().copied();
    let before = system.total_energy_j().unwrap();
    let original_config = system.config();
    let original_compliance = system.compliance_m_per_n();
    let original_columns: [*const f64; 3] =
        core::array::from_fn(|i| system.work_columns()[i].as_ptr());

    let report = system
        .apply_temperature_successor(&hot, "T", &CancelGate::new())
        .unwrap();
    let squared_jump = old.jump_m().iter().map(|g| g * g).sum::<f64>();
    let radius = det::sqrt(squared_jump);
    let expected_damage =
        1.0 - (0.00016 / radius) * ((0.00125 - radius) / (0.00125 - 0.00016));
    let expected_work = 0.5 * AREA * (5e8 - 1e9) * (1.0 - old.damage()) * squared_jump;
    let expected_release = 0.5 * AREA * 5e8 * (expected_damage - old.damage()) * squared_jump;
    assert!(expected_release > 0.0);
    close(report.parameter_work_j, expected_work, 2e-13);
    close(report.damage_dissipation_j, expected_release, 2e-13);
    close(report.stored_change_j, expected_work - expected_release, 2e-13);
    close(report.energy_before_j, before, 1e-15);
    close(report.energy_after_j, system.total_energy_j().unwrap(), 1e-15);
    close(system.point().damage(), expected_damage, 2e-10);
    close(
        AREA * (system.point().dissipated_j_m2() - old.dissipated_j_m2()),
        expected_release,
        2e-13,
    );
    assert_eq!(report.sample, sample);
    assert_eq!(report.jump_m, old.jump_m());
    assert_eq!(system.point().jump_m(), old.jump_m());
    assert_eq!(states(&system), native);
    assert_eq!(system.samples_rendered(), sample);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.config(), original_config);
    assert_eq!(system.compliance_m_per_n(), original_compliance);
    for i in 0..3 {
        assert_eq!(system.work_columns()[i].as_ptr(), original_columns[i]);
    }
    assert!(core::ptr::eq(system.source().unwrap(), &hot));
    assert!(core::ptr::eq(system.law(), hot.law()));
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);

    // Analytic unloading at retained damage distinguishes the hot stiffness
    // from stale cold coefficients on the next ordinary native timestep.
    let heated = system.point();
    let target = [3e-4, 0.0, 0.0];
    let expected_force =
        AREA * 5e8 * (1.0 - heated.damage()) * 0.5 * (heated.jump_m()[0] + target[0]);
    let frame = drive_with_force(&mut system, target, [expected_force, 0.0, 0.0]);
    close(frame.mean_force_n[0], expected_force, 2e-9);
    for i in 0..3 {
        close(frame.jump_m[i], target[i], 3e-12);
    }
    close(frame.fracture_dissipation_j, 0.0, 2e-13);
    assert_eq!(system.samples_rendered(), sample + 1);
    assert!(core::ptr::eq(system.source().unwrap(), &hot));
    let before_cooling = system.point();
    let report = system
        .apply_temperature_successor(&cold, "T", &CancelGate::new())
        .unwrap();
    assert!(system.point().damage() >= before_cooling.damage());
    close(system.point().damage(), before_cooling.damage(), 2e-13);
    close(
        system.point().dissipated_j_m2(),
        before_cooling.dissipated_j_m2(),
        2e-10,
    );
    close(report.damage_dissipation_j, 0.0, 2e-13);
    assert!(report.parameter_work_j > 0.0);
    assert_eq!(system.point().jump_m(), before_cooling.jump_m());
}

#[test]
fn unchanged_temperature_source_is_zero_work_and_keeps_exact_native_replay() {
    let card = card([0.3, 0.1]);
    let (source, _) = sources(&card);
    let mut system = bound(&source, [0.015, -0.01, 0.005], 1.0);
    let mut replay = bound(&source, [0.015, -0.01, 0.005], 1.0);
    for owner in [&mut system, &mut replay] {
        drive_to(owner, [4e-4, 0.0, 0.0]);
    }
    let point = system.point();
    let last = system.last_frame().copied();
    let report = system
        .apply_temperature_successor(&source, "T", &CancelGate::new())
        .unwrap();
    assert_eq!(report.parameter_work_j, 0.0);
    assert_eq!(report.damage_dissipation_j, 0.0);
    assert_eq!(report.stored_change_j, 0.0);
    assert_eq!(report.energy_residual_j, 0.0);
    assert_eq!(report.energy_before_j, report.energy_after_j);
    assert_eq!(system.point(), point);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(states(&system), states(&replay));
    assert_eq!(system.samples_rendered(), replay.samples_rendered());
    let drive = [-0.02, 0.015, -0.007];
    assert_eq!(
        *system.step(&[0.0; 6], drive).unwrap(),
        *replay.step(&[0.0; 6], drive).unwrap()
    );
    assert_eq!(system.point(), replay.point());
    assert_eq!(states(&system), states(&replay));
}

#[test]
fn paired_temperature_continues_both_original_friction_modes_and_hot_sliding_force() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    for fully_failed in [false, true] {
        let mut system = bound(&cold, [0.01, 0.02, -0.01], 1.0);
        drive_to(&mut system, [if fully_failed { 0.002 } else { 4e-4 }, 0.0, 0.0]);
        drive_to(&mut system, [-2e-4, 0.0, 0.0]);
        if fully_failed {
            system
                .enable_frictional_reclosure(
                    &friction,
                    receiving(&cold),
                    ReclosureActivation::Sliding,
                    reclosure_config(),
                    &CancelGate::new(),
                )
                .unwrap();
        } else {
            system
                .enable_frictional_damage(
                    &friction,
                    receiving(&cold),
                    ReclosureActivation::Sliding,
                    damage_config(),
                    &CancelGate::new(),
                )
                .unwrap();
        }
        let old = system.point();
        let native = states(&system);
        let last = system.last_frame().copied();
        let sample = system.samples_rendered();
        let residual_occupancy = system.residual_contact_occupancy();
        let damaged_occupancy = system.damage_friction_occupancy();
        let report = system
            .apply_temperature_and_friction_successors(
                &hot,
                &hot_friction,
                "T",
                &CancelGate::new(),
            )
            .unwrap();
        close(report.parameter_work_j, 20e-6, 2e-13);
        close(report.damage_dissipation_j, 0.0, 2e-13);
        assert_eq!(states(&system), native);
        assert_eq!(system.last_frame().copied(), last);
        assert_eq!(system.samples_rendered(), sample);
        assert_eq!(system.point().jump_m(), old.jump_m());
        close(system.point().damage(), old.damage(), 2e-13);
        close(system.point().dissipated_j_m2(), old.dissipated_j_m2(), 2e-10);
        assert_eq!(system.residual_contact_occupancy(), residual_occupancy);
        assert_eq!(system.damage_friction_occupancy(), damaged_occupancy);
        assert!(core::ptr::eq(system.source().unwrap(), &hot));
        let actual_source = if fully_failed {
            assert_eq!(residual_occupancy, Some(ReclosureRegime::Sliding));
            system.residual_friction_source().unwrap()
        } else {
            assert_eq!(damaged_occupancy, Some(DamageFrictionRegime::Sliding));
            system.damage_friction_source().unwrap()
        };
        assert!(core::ptr::eq(actual_source, &hot_friction));

        let damage = system.point().damage();
        let slip = 1e-4;
        let capacity = 0.1 * damage * 0.6;
        let cohesive = 0.5 * AREA * 5e8 * (1.0 - damage) * slip;
        let force = [-0.6, 0.6 * (capacity + cohesive), 0.8 * (capacity + cohesive)];
        let frame = drive_with_force(&mut system, [-2e-4, 6e-5, 8e-5], force);
        for i in 0..3 {
            close(frame.mean_force_n[i], force[i], 3e-9);
        }
        close(frame.fracture_dissipation_j, 0.0, 2e-13);
        if fully_failed {
            assert_eq!(frame.residual_friction_regime, Some(ReclosureRegime::Sliding));
            close(frame.residual_friction_dissipation_j, capacity * slip, 3e-13);
        } else {
            let friction_frame = frame.damage_friction.unwrap();
            assert_eq!(friction_frame.regime, DamageFrictionRegime::Sliding);
            close(friction_frame.damaged_compression_force_n, damage * 0.6, 3e-9);
            close(friction_frame.dissipation_j, capacity * slip, 3e-13);
        }
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
}

#[test]
fn invalid_pair_receiver_refusal_and_late_cancellation_are_atomic_and_retryable() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let foreign_card = self::card([0.4, 0.1]);
    let (foreign, foreign_friction) = sources(&foreign_card);
    let mut system = bound(&cold, [0.0; 3], 1.0);
    let mut replay = bound(&cold, [0.0; 3], 1.0);
    for owner in [&mut system, &mut replay] {
        drive_to(owner, [4e-4, 0.0, 0.0]);
        drive_to(owner, [-2e-4, 0.0, 0.0]);
        owner
            .enable_frictional_damage(
                &friction,
                receiving(&cold),
                ReclosureActivation::Sliding,
                damage_config(),
                &CancelGate::new(),
            )
            .unwrap();
    }
    let original_point = system.point();
    let original_states = states(&system);
    let original_last = system.last_frame().copied();
    let original_sample = system.samples_rendered();
    let unchanged = |owner: &CohesiveModalSystem<'_>| {
        assert_eq!(owner.point(), original_point);
        assert_eq!(states(owner), original_states);
        assert_eq!(owner.last_frame().copied(), original_last);
        assert_eq!(owner.samples_rendered(), original_sample);
        assert!(core::ptr::eq(owner.source().unwrap(), &cold));
        assert!(core::ptr::eq(owner.damage_friction_source().unwrap(), &friction));
        assert_eq!(owner.damage_friction_occupancy(), Some(DamageFrictionRegime::Sliding));
    };
    assert!(system.apply_temperature_successor(&hot, "T", &CancelGate::new()).is_err());
    unchanged(&system);
    for (next, next_friction, axis) in [
        (&hot, &friction, "T"),
        (&cold, &hot_friction, "T"),
        (&foreign, &foreign_friction, "T"),
        (&hot, &hot_friction, "pressure"),
    ] {
        assert!(
            system
                .apply_temperature_and_friction_successors(
                    next,
                    next_friction,
                    axis,
                    &CancelGate::new(),
                )
                .is_err()
        );
        unchanged(&system);
    }

    let mut called = false;
    let refused = system.apply_temperature_successors_with_admission(
        &hot,
        Some(&hot_friction),
        "T",
        &CancelGate::new(),
        |report, point| {
            called = true;
            assert!(core::ptr::eq(point.law(), hot.law()));
            close(report.parameter_work_j, 20e-6, 2e-13);
            Err::<(), Box<dyn std::error::Error>>("receiving energy owner refused".into())
        },
    );
    assert!(called && refused.is_err());
    unchanged(&system);
    let gate = CancelGate::new();
    let mut called = false;
    let cancelled = system.apply_temperature_successors_with_admission(
        &hot,
        Some(&hot_friction),
        "T",
        &gate,
        |_, _| {
            called = true;
            gate.request();
            Ok::<(), CohesiveModalError>(())
        },
    );
    assert!(called && cancelled.is_err());
    unchanged(&system);

    // An absent friction mode cannot be created by a temperature operation.
    let mut plain = bound(&cold, [0.0; 3], 1.0);
    let plain_states = states(&plain);
    assert!(
        plain
            .apply_temperature_and_friction_successors(
                &hot,
                &hot_friction,
                "T",
                &CancelGate::new(),
            )
            .is_err()
    );
    assert_eq!(states(&plain), plain_states);
    assert_eq!(plain.samples_rendered(), 0);
    assert!(core::ptr::eq(plain.source().unwrap(), &cold));

    let accepted = system
        .apply_temperature_and_friction_successors(&hot, &hot_friction, "T", &CancelGate::new())
        .unwrap();
    let replayed = replay
        .apply_temperature_and_friction_successors(&hot, &hot_friction, "T", &CancelGate::new())
        .unwrap();
    assert_eq!(accepted, replayed);
    let drive = [0.0, 0.01, -0.005];
    assert_eq!(
        *system.step(&[0.0; 6], drive).unwrap(),
        *replay.step(&[0.0; 6], drive).unwrap()
    );
    assert_eq!(states(&system), states(&replay));
    assert_eq!(system.point(), replay.point());
}

#[test]
fn native_temperature_energy_ceiling_rejects_then_admits_a_smaller_update() {
    let card = card([0.3, 0.1]);
    let (cold, _) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let mild = cold.at_temperature(&card, "T", 296.0).unwrap();
    let ceiling =
        0.5 * (MA + MB) * (0.012_f64.powi(2) + 0.006_f64.powi(2) + 0.003_f64.powi(2)) + 50e-6;
    let mut system = bound(&cold, [0.0; 3], ceiling);
    drive_to(&mut system, [-2e-4, 0.0, 0.0]);
    let old = system.point();
    let native = states(&system);
    let last = system.last_frame().copied();
    let sample = system.samples_rendered();
    assert!(matches!(
        system.apply_temperature_successor(&hot, "T", &CancelGate::new()),
        Err(CohesiveModalError::Network(ModalCouplingError::Budget {
            what: "temperature-refreshed native cohesive energy",
            ..
        }))
    ));
    assert_eq!(system.point(), old);
    assert_eq!(states(&system), native);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.samples_rendered(), sample);
    assert!(core::ptr::eq(system.source().unwrap(), &cold));
    let report = system
        .apply_temperature_successor(&mild, "T", &CancelGate::new())
        .unwrap();
    close(report.parameter_work_j, 8e-6, 2e-13);
    assert!(report.energy_after_j <= ceiling);
    assert_eq!(states(&system), native);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.samples_rendered(), sample);
    assert!(core::ptr::eq(system.source().unwrap(), &mild));
}
