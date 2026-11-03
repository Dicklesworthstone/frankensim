//! G1/G3/G4: prescribed thermal weakening, real inertia and finite adhesive heat.
use super::*;
use fs_solid::cohesive::dynamics::CohesiveMotion;
use fs_solid::cohesive::temperature::TemperatureConfig;

const GAP: f64 = 0.0006;
const OFFSETS: [[f64; 3]; 1] = [[-0.001, GAP, 0.0]];
const DT: f64 = 0.01;
const SURVIVAL: f64 = 17.0 / 27.0;

fn temperature_material(environment: &str) -> (ResolvedBilinearCohesive, InterfaceSystemCard) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, cold, hot) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1000.0,
            500.0,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 1.0, 0.2),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            0.002,
            0.0004,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2000.0,
            2000.0,
        ),
        (
            COHESIVE_NORMAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            2.0,
            1.0,
        ),
        (
            COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
            4.0,
            2.0,
        ),
    ] {
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, cold), (320.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal weakening and finite-layer reference".into(),
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
            process: "original finite member".into(),
            revision: 1,
        },
        texture_frame: "fixed normal/shear axes".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "original bond".into(),
            third_body: Some(ADHESIVE.into()),
            environment: environment.into(),
            history: "initially bonded".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let resolved = resolve_interface_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = ResolvedBilinearCohesive::try_from_viscous_resolved(
        &card,
        resolved,
        CohesiveEnvelope::CommonMode,
        0.02,
        1.0,
    )
    .unwrap();
    (source, card)
}

fn source_temperature<'s, 'a>(
    target: &'s CohesiveSystem<'a>,
) -> CohesiveLayerHeatTemperature<'s, 'a> {
    CohesiveLayerHeatTemperature {
        target,
        axis: "T",
        controls: TemperatureConfig {
            local_energy_tolerance_j_m2: 1e-12,
            total_energy_tolerance_j: 1e-12,
        },
    }
}

fn loaded_motion<'s, 'a>(
    joint: &'s CohesiveSystem<'a>,
    mass: &'s Scalar,
) -> CohesiveMotion<'s, 'a> {
    // The explicit one-second preparation includes the original physical
    // dashpot force. Its completed loading precedes the motion/thermal clock.
    let state = joint
        .bonded_state()
        .solve_over_time_controlled(
            &[(1000.0 + 4.0) * GAP],
            &OFFSETS,
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state;
    assert_eq!(state.displacements_m(), &[0.0]);
    assert_eq!(state.points()[0].damage(), 0.0);
    state
        .into_dynamics(mass, vec![-0.01], 0.0, dynamics())
        .unwrap()
}

// Midpoint momentum of a retained-damage spring and parallel physical dashpot.
fn displacement_increment(u: f64, velocity: f64, stiffness: f64, viscosity: f64) -> f64 {
    (2.0 * 10.0 * velocity / DT - 1000.0 * u - stiffness * (GAP + u))
        / (2.0 * 10.0 / (DT * DT) + 0.5 * (1000.0 + stiffness) + viscosity / DT)
}

#[test]
fn thermal_weakening_changes_actual_motion_and_only_viscosity_heats_the_finite_layer() {
    let (cold, card) = temperature_material("manufactured prescribed temperature");
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &cold, &[7], true);
    let hot_joint = joint
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let motion = loaded_motion(&joint, &mass);
    let point = motion.state().points()[0];
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let u = displacement_increment(0.0, -0.01, 500.0 * SURVIVAL, 2.0 * SURVIVAL);
    let viscous = 2.0 * SURVIVAL * u * u / DT;
    let next = propose_deformable_cohesive_layer_heat_with_temperature(
        &motion,
        &[0.0],
        &OFFSETS,
        DT,
        source_temperature(&hot_joint),
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, mechanics| {
            close(
                mechanics.motion.state().points()[0].damage(),
                1.0 - SURVIVAL,
                1e-12,
            );
            close(mechanics.viscous_dissipated_j, viscous, 1e-12);
            Ok::<_, &'static str>(transport())
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let accepted = &next.heated;
    let final_motion = &accepted.mechanical.motion;
    close(final_motion.state().displacements_m()[0], u, 1e-12);
    close(final_motion.velocity_m_s()[0], 2.0 * u / DT + 0.01, 1e-11);
    assert!((u - displacement_increment(0.0, -0.01, 1000.0, 4.0)).abs() > 1e-6);
    // At held gap .0006, hot onset=.0004, failure=.004: d=10/27.
    close(next.temperature.parameter_work_j, -0.00009, 1e-13);
    close(next.temperature.damage_dissipated_j, 1.0 / 30000.0, 1e-13);
    close(next.fracture_work_j, 1.0 / 30000.0, 1e-13);
    close(accepted.fracture_work_j, 0.0, 1e-13);
    close(accepted.heat.dissipated_j, viscous, 1e-12);
    close(
        accepted.thermal.energy_j() - thermal.energy_j(),
        viscous,
        1e-11,
    );
    assert!(next.fracture_work_j > 20.0 * viscous);
    let expected = reference(DT, viscous, [300.0; 3]);
    for i in 0..3 {
        close(
            accepted.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected[i],
            2e-10,
        );
    }
    assert!(accepted.heat.layer_stored_change_j > 0.0);
    assert!(accepted.heat.layer_stored_change_j < viscous);
    close(
        accepted.heat.entropy.unwrap().friction_generation_j_k,
        viscous / expected[2],
        1e-14,
    );
    let before = 0.5 * mass.0 * 0.01_f64.powi(2) + point.stored_energy_j_m2();
    let after = 0.5 * mass.0 * final_motion.velocity_m_s()[0].powi(2)
        + 0.5 * elastic.0 * u * u
        + final_motion.state().points()[0].stored_energy_j_m2()
        + accepted.thermal.energy_j();
    close(
        after - before + next.fracture_work_j - next.temperature.parameter_work_j,
        0.0,
        2e-10,
    );
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    close(next.summed_energy_residual_j, 0.0, 2e-10);
    assert!(core::ptr::eq(
        final_motion.state().system().patches()[0].source().unwrap(),
        &hot
    ));
    assert_eq!(motion.state().points()[0], point);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn cooling_retains_thermal_damage_and_continues_the_original_inertia_and_heat_store() {
    let (cold, card) = temperature_material("manufactured prescribed temperature");
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &cold, &[7], true);
    let hot_joint = joint
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let motion = loaded_motion(&joint, &mass);
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let warm = propose_deformable_cohesive_layer_heat_with_temperature(
        &motion,
        &[0.0],
        &OFFSETS,
        DT,
        source_temperature(&hot_joint),
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let old = &warm.heated.mechanical.motion;
    let old_point = old.state().points()[0];
    let old_u = old.state().displacements_m()[0];
    let du = displacement_increment(
        old_u,
        old.velocity_m_s()[0],
        1000.0 * SURVIVAL,
        4.0 * SURVIVAL,
    );
    let viscous = 4.0 * SURVIVAL * du * du / DT;
    let expected_work = 0.5 * 500.0 * SURVIVAL * old_point.jump_m()[1].powi(2);
    let before_temperatures =
        std::array::from_fn(|i| warm.heated.thermal.temperature_k(CellId(i as u64)).unwrap());
    let cooled = propose_deformable_cohesive_layer_heat_with_temperature(
        old,
        &[0.0],
        &OFFSETS,
        DT,
        source_temperature(&joint),
        binding(&route, &card),
        &warm.heated.thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let actual = &cooled.heated.mechanical.motion;
    close(actual.state().displacements_m()[0], old_u + du, 1e-12);
    close(
        actual.velocity_m_s()[0],
        2.0 * du / DT - old.velocity_m_s()[0],
        1e-11,
    );
    close(actual.state().points()[0].damage(), 1.0 - SURVIVAL, 1e-12);
    close(
        actual.state().points()[0].dissipated_j_m2(),
        old_point.dissipated_j_m2(),
        1e-13,
    );
    close(cooled.temperature.parameter_work_j, expected_work, 1e-13);
    assert!(cooled.temperature.parameter_work_j > 0.0);
    close(cooled.temperature.damage_dissipated_j, 0.0, 1e-13);
    close(cooled.fracture_work_j, 0.0, 1e-13);
    close(cooled.heated.heat.dissipated_j, viscous, 1e-12);
    let expected = reference(DT, viscous, before_temperatures);
    for i in 0..3 {
        close(
            cooled
                .heated
                .thermal
                .temperature_k(CellId(i as u64))
                .unwrap(),
            expected[i],
            2e-10,
        );
    }
    close(
        cooled.heated.thermal.energy_j() - warm.heated.thermal.energy_j(),
        viscous,
        1e-11,
    );
    close(cooled.combined_energy_residual_j, 0.0, 2e-10);
    close(cooled.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(actual.time_s(), 2.0 * DT);
    assert_eq!(actual.steps(), 2);
    assert_eq!(cooled.heated.thermal.time_s(), actual.time_s());
    assert!(core::ptr::eq(
        actual.state().system().patches()[0].source().unwrap(),
        &cold
    ));
    assert_eq!(old.state().points()[0], old_point);
    assert_eq!(old.time_s(), DT);
}

#[test]
fn foreign_successors_and_late_source_or_cancellation_failure_preserve_both_original_owners() {
    let (cold, card) = temperature_material("manufactured prescribed temperature");
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let (foreign, foreign_card) = temperature_material("different original environment");
    let foreign_hot = foreign.at_temperature(&foreign_card, "T", 320.0).unwrap();
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &cold, &[7], true);
    let hot_joint = joint
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let stale_joint = system(&elastic, &foreign_hot, &[7], true);
    let motion = loaded_motion(&joint, &mass);
    let point = motion.state().points()[0];
    let thermal = layer_network();
    let cells = thermal.cells().clone();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let calls = Cell::new(0);
    let stale = propose_deformable_cohesive_layer_heat_with_temperature(
        &motion,
        &[0.0],
        &OFFSETS,
        DT,
        source_temperature(&stale_joint),
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(transport())
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(stale, Err(CohesiveHeatError::Mechanical(_))));
    assert_eq!(calls.get(), 0);
    for cancel in [false, true] {
        let calls = Cell::new(0);
        let gate = CancelGate::new_clock_free();
        let refused = propose_deformable_cohesive_layer_heat_with_temperature(
            &motion,
            &[0.0],
            &OFFSETS,
            DT,
            source_temperature(&hot_joint),
            binding(&route, &card),
            &thermal,
            layer_controls(),
            &[0.0; 3],
            |_, _, mechanics| {
                close(
                    mechanics.motion.state().points()[0].damage(),
                    1.0 - SURVIVAL,
                    1e-12,
                );
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("outside conductivity source");
                    }
                }
                Ok(transport())
            },
            &gate,
        );
        if cancel {
            assert!(refused.is_err());
        } else {
            assert!(matches!(
                refused,
                Err(CohesiveHeatError::Source("outside conductivity source"))
            ));
        }
        assert!(calls.get() >= 2);
        assert_eq!(motion.state().points()[0], point);
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[-0.01]);
        assert_eq!(motion.time_s(), 0.0);
        assert_eq!(motion.steps(), 0);
        assert!(core::ptr::eq(
            motion.state().system().patches()[0].source().unwrap(),
            &cold
        ));
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
}
