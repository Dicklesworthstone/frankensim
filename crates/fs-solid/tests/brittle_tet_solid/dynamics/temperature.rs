//! Original thermal coefficient/history transfer during actual brittle motion.
use super::*;

use fs_solid::brittle::solid::dynamics::temperature::BrittleTetDynamicTemperatureError;
use fs_solid::brittle::solid::dynamics::time::{BrittleTetTimeControls, BrittleTetTimeError};

fn density_curve_card(tag: &str) -> MaterialCard {
    let base = source_card(STRENGTH, tag, true);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 2000.0), (600.0, 1000.0)],
                dims: Density::DIMS,
            },
            // Typed like every query of this card: an untyped "T" axis cannot
            // answer a query that types it in kelvin.
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                600.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic temperature-dependent bulk density".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn translated_supports(
    solid: &BrittleTetSolid<'_>,
    velocity_mps: [f64; 3],
    dt_s: f64,
) -> Vec<(usize, f64)> {
    solid
        .solution()
        .displacement_m
        .as_flattened()
        .iter()
        .enumerate()
        .map(|(dof, &old)| (dof, old + dt_s * velocity_mps[dof % 3]))
        .collect()
}

#[test]
fn heating_a_translating_preloaded_solid_initiates_fracture_without_an_impulse() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = density_curve_card("moving thermal onset");
        let cold = [dynamic_material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let strain = 0.75 * STRENGTH / YOUNG;
        let preloaded = virgin
            .advance(
                &prescribed(0.03, strain),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert_eq!(preloaded.solution().element_states[0].damage(), 0.0);
        let velocity = [0.1, -0.2, 0.05];
        let old = BrittleTetMotion::try_from_solid(&preloaded, &[velocity; 4], 2.0, 3, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let dt = 1e-4;
        let supports = translated_supports(old.solid(), velocity, dt);
        let step = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let moved = &step.mechanical.motion;
        let volume = geometry.element_volumes_m3()[0];
        let hot_strength = 0.25 * STRENGTH;
        let hot_young = 0.5 * YOUNG;
        let failure =
            2.0 * FRACTURE_WORK / (hot_strength * moved.solid().band_geometry()[0].band_width_m());
        let stress = hot_strength * (failure - strain) / (failure - hot_strength / hot_young);
        let stored = 0.5 * stress * strain * volume;
        let parameter_work = 0.5 * (hot_young - YOUNG) * strain * strain * volume;
        let induced_fracture = 0.5 * (hot_young * strain * strain - stress * strain) * volume;
        assert!(induced_fracture > 0.0);
        close(step.transfer.parameter_work_j, parameter_work, 2e-9);
        close(step.transfer.fracture_dissipation_j, induced_fracture, 2e-9);
        close(step.transfer.stored_after_j, stored, 2e-9);
        close(moved.solid().stored_energy_j(), stored, 2e-9);
        close(moved.solid().solution().element_stress[0][0], stress, 2e-9);
        close(moved.solid().solution().element_strain[0][0], strain, 2e-9);
        assert!(moved.solid().solution().element_states[0].damage() > 0.0);
        assert!(step.transfer.parameter_work_j < 0.0);
        assert!(step.mechanical.energy.material_work_j.abs() < 1e-12);
        assert!(step.mechanical.energy.fracture_dissipation_j.abs() < 1e-12);
        assert!(step.mechanical.energy.inertia.kinetic_change_j.abs() < 1e-12);
        assert!(step.mechanical.energy.inertia.numerical_loss_j < 1e-20);
        assert!(step.total_balance_error_j.abs() < 1e-9);
        for v in moved.velocity_mps() {
            for axis in 0..3 {
                close(v[axis], velocity[axis], 2e-9);
            }
        }
        close(moved.kinetic_energy_j(), old.kinetic_energy_j(), 2e-9);
        assert_eq!(moved.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(moved.reference_density_kg_m3(), &[2000.0]);
        assert_eq!(
            moved.solid().materials()[0]
                .reference_density_kg_m3()
                .unwrap(),
            1000.0
        );
        assert!(std::ptr::eq(
            moved.solid().materials()[0].law(),
            hot[0].law()
        ));
        close(
            step.element_transfer[0].damage_dissipated_j_m3 * volume,
            induced_fracture,
            2e-9,
        );
        assert_eq!(moved.accepted_steps(), 1);
        assert_eq!(moved.remaining_steps(), 2);
        assert_eq!(moved.time_s(), 2.0 + dt);
        assert_eq!(old.time_s(), 2.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
    });
}

#[test]
fn new_temperature_density_cannot_change_free_body_mass_or_force_acceleration() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = density_curve_card("thermal reference inventory");
        let cold = [dynamic_material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let solid =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let velocity = [0.2, -0.1, 0.05];
        let old = BrittleTetMotion::try_from_solid(&solid, &[velocity; 4], 0.0, 3, cx).unwrap();
        let dt = 1e-6;
        let force = 2.0;
        let loads = uniform_x_force(force);
        let mass = 2000.0 * geometry.element_volumes_m3()[0];
        let change = dt * force / mass;
        let heated = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &[],
                &loads,
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let moving = &heated.mechanical.motion;
        close(moving.reference_mass_kg(), mass, 1e-13);
        assert_eq!(moving.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(moving.reference_density_kg_m3(), &[2000.0]);
        assert_eq!(hot[0].reference_density_kg_m3().unwrap(), 1000.0);
        for v in moving.velocity_mps() {
            close(v[0], velocity[0] + change, 2e-8);
            close(v[1], velocity[1], 2e-8);
            close(v[2], velocity[2], 2e-8);
        }
        close(
            heated.mechanical.energy.dead_load_endpoint_work_j,
            force * dt * (velocity[0] + change),
            2e-8,
        );
        close(
            heated.mechanical.energy.inertia.numerical_loss_j,
            0.5 * mass * change * change,
            2e-7,
        );
        assert_eq!(heated.transfer.parameter_work_j, 0.0);
        assert_eq!(heated.transfer.fracture_dissipation_j, 0.0);
        assert!(heated.total_balance_error_j.abs() < 1e-10);

        // Ordinary later steps keep that exact inventory too, although their
        // current source bundle still reports the hot bulk density.
        let continued = moving
            .advance(
                dt,
                &[],
                &loads,
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        for v in continued.motion.velocity_mps() {
            close(v[0], velocity[0] + 2.0 * change, 2e-8);
        }
        assert_eq!(
            continued.motion.reference_mass_kg(),
            old.reference_mass_kg()
        );
        assert_eq!(continued.motion.reference_density_kg_m3(), &[2000.0]);
        assert_eq!(continued.motion.accepted_steps(), 2);
        assert_eq!(continued.motion.time_s(), 2.0 * dt);
        assert_eq!(old.velocity_mps(), &[velocity; 4]);
    });
}

#[test]
fn cooling_retains_damage_fracture_work_and_the_original_moving_inventory() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = density_curve_card("moving thermal cycle");
        let cold = [dynamic_material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let cooled = [hot[0].at_temperature(&source, "T", 300.0).unwrap()];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let damaged = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let original_loss = damaged.fracture_dissipation_j();
        assert!(original_loss > 0.0);
        let velocity = [0.05, 0.0, 0.0];
        let old = BrittleTetMotion::try_from_solid(&damaged, &[velocity; 4], 0.0, 3, cx).unwrap();
        let dt = 1e-4;
        let heated = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &translated_supports(old.solid(), velocity, dt),
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let warm = &heated.mechanical.motion;
        let warm_damage = warm.solid().solution().element_states[0].damage();
        let warm_loss = warm.solid().fracture_dissipation_j();
        let warm_storage = warm.solid().stored_energy_j();
        let cooling = warm
            .advance_at_temperature(
                dt,
                &cooled,
                "T",
                &translated_supports(warm.solid(), velocity, dt),
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let next = &cooling.mechanical.motion;
        assert!(next.solid().solution().element_states[0].damage() >= warm_damage);
        close(
            next.solid().solution().element_states[0].damage(),
            warm_damage,
            2e-10,
        );
        assert!(next.solid().fracture_dissipation_j() >= original_loss);
        close(next.solid().fracture_dissipation_j(), warm_loss, 2e-9);
        close(next.solid().stored_energy_j(), 2.0 * warm_storage, 2e-9);
        close(cooling.transfer.parameter_work_j, warm_storage, 2e-9);
        assert!(cooling.transfer.parameter_work_j > 0.0);
        assert!(cooling.transfer.fracture_dissipation_j < 1e-12);
        assert!(cooling.mechanical.energy.fracture_dissipation_j.abs() < 1e-12);
        assert!(cooling.total_balance_error_j.abs() < 1e-9);
        assert_eq!(next.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(next.reference_density_kg_m3(), &[2000.0]);
        close(next.kinetic_energy_j(), old.kinetic_energy_j(), 2e-9);
        for v in next.velocity_mps() {
            close(v[0], velocity[0], 2e-9);
            assert_eq!(v[1..], [0.0; 2]);
        }
        assert_eq!(next.accepted_steps(), 2);
        assert_eq!(next.time_s(), 2.0 * dt);
        assert_eq!(warm.solid().fracture_dissipation_j(), warm_loss);
        assert_eq!(
            warm.solid().solution().element_states[0].damage(),
            warm_damage
        );
        assert_eq!(warm.time_s(), dt);
    });
}

#[test]
fn late_inertial_source_time_and_budget_refusals_keep_the_cold_motion_retryable() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = density_curve_card("thermal motion rollback");
        let cold = [dynamic_material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let preloaded = virgin
            .advance(
                &prescribed(0.03, 0.75 * STRENGTH / YOUNG),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let old = BrittleTetMotion::try_from_solid(&preloaded, &[[0.0; 3]; 4], 4.0, 1, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let storage = old.solid().stored_energy_j();
        let dt = 1e-4;
        let supports = translated_supports(old.solid(), [0.1, 0.0, 0.0], dt);
        let strict = BrittleTetDynamicControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..dynamic_controls()
        };
        // This gate occurs after the actual source transfer and native dynamic
        // candidate, including temperature-induced damage, have been computed.
        assert!(matches!(
            old.advance_at_temperature(
                dt,
                &hot,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                strict,
                cx,
            ),
            Err(BrittleTetDynamicTemperatureError::Dynamic(
                BrittleTetDynamicError::Solid(BrittleTetError::EnergyGate {
                    account: "backward-Euler numerical kinetic loss",
                    ..
                })
            ))
        ));
        let other = density_curve_card("different original material");
        let other_cold = dynamic_material(&other, &geometry);
        let wrong = [other_cold.at_temperature(&other, "T", 600.0).unwrap()];
        assert!(matches!(
            old.advance_at_temperature(
                dt,
                &wrong,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            ),
            Err(BrittleTetDynamicTemperatureError::Transfer(
                BrittleTetTemperatureError::Source { element: 0, .. }
            ))
        ));
        assert!(matches!(
            old.advance_at_temperature(
                f64::NAN,
                &hot,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            ),
            Err(BrittleTetDynamicTemperatureError::Dynamic(
                BrittleTetDynamicError::Invalid(_)
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                old.advance_at_temperature(
                    dt,
                    &hot,
                    "T",
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    dynamic_controls(),
                    cancelled,
                ),
                Err(BrittleTetDynamicTemperatureError::Dynamic(
                    BrittleTetDynamicError::Solid(BrittleTetError::Cancelled)
                ))
            ));
        });
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.solid().stored_energy_j(), storage);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.kinetic_energy_j(), 0.0);
        assert_eq!(old.time_s(), 4.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.reference_density_kg_m3(), &[2000.0]);

        let accepted = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let repeated = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.transfer, repeated.transfer);
        assert_eq!(accepted.mechanical.energy, repeated.mechanical.energy);
        assert_eq!(
            accepted.mechanical.motion.solid().solution(),
            repeated.mechanical.motion.solid().solution()
        );
        assert!(accepted.transfer.fracture_dissipation_j > 0.0);
        assert!(accepted.mechanical.energy.inertia.numerical_loss_j > 0.0);
        assert_eq!(accepted.mechanical.motion.accepted_steps(), 1);
        assert_eq!(accepted.mechanical.motion.remaining_steps(), 0);
        assert_eq!(accepted.mechanical.motion.time_s(), 4.0 + dt);
        assert!(matches!(
            accepted.mechanical.motion.advance_at_temperature(
                dt,
                &cold,
                "T",
                &supports,
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            ),
            Err(BrittleTetDynamicTemperatureError::Dynamic(
                BrittleTetDynamicError::StepBudget {
                    accepted_steps: 1,
                    maximum_steps: 1
                }
            ))
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.time_s(), 4.0);
    });
}

#[test]
fn adaptive_heated_motion_keeps_original_mass_and_refuses_partial_duration_publication() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = density_curve_card("adaptive thermal motion");
        let cold = [dynamic_material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let solid =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let velocity = [0.2, -0.1, 0.05];
        let old = BrittleTetMotion::try_from_solid(&solid, &[velocity; 4], 0.0, 8, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let duration = 4e-6;
        let force = 2.0;
        let forces: Vec<_> = (0..4).map(|node| (3 * node, [force / 4.0; 2])).collect();
        let controls = BrittleTetTimeControls {
            initial_interval_s: 2e-6,
            minimum_interval_s: 2e-6,
            maximum_interval_s: 2e-6,
            maximum_intervals: 2,
            maximum_trials: 6,
            displacement_tolerance_m: 1e-8,
            velocity_tolerance_mps: 1e-6,
            damage_tolerance: 1e-6,
            fracture_dissipation_tolerance_j: 1e-10,
            relative_tolerance: 0.0,
            maximum_damage_increment: 0.05,
            work: dynamic_controls(),
        };
        // One complete fine pair can be accepted internally, but an exhausted
        // duration budget must not publish its new source, time or velocity.
        assert!(matches!(
            old.advance_time_at_temperature(
                duration,
                &hot,
                "T",
                &[],
                &forces,
                IncrementSettings::default(),
                BrittleTetTimeControls {
                    maximum_intervals: 1,
                    ..controls
                },
                cx,
            ),
            Err(BrittleTetDynamicTemperatureError::Time(
                BrittleTetTimeError::WorkLimit {
                    accepted_intervals: 1,
                    trials: 3,
                }
            ))
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[velocity; 4]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(
            old.solid().materials()[0].material().state_coordinate("T"),
            Some(300.0)
        );

        let accepted = old
            .advance_time_at_temperature(
                duration,
                &hot,
                "T",
                &[],
                &forces,
                IncrementSettings::default(),
                controls,
                cx,
            )
            .unwrap();
        let motion = &accepted.mechanical.motion;
        let mass = 2000.0 * geometry.element_volumes_m3()[0];
        let change = force * duration / mass;
        let half_dt = 1e-6;
        let numerical_loss = 4.0 * 0.5 * mass * (force * half_dt / mass).powi(2);
        assert_eq!(accepted.mechanical.records.len(), 2);
        assert_eq!(accepted.mechanical.trials, 6);
        assert_eq!(accepted.mechanical.rejected_intervals, 0);
        assert_eq!(motion.accepted_steps(), 4);
        assert_eq!(motion.remaining_steps(), 4);
        assert_eq!(motion.time_s(), duration);
        assert_eq!(motion.reference_density_kg_m3(), &[2000.0]);
        assert_eq!(motion.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(
            motion.solid().materials()[0]
                .reference_density_kg_m3()
                .unwrap(),
            1000.0
        );
        assert!(std::ptr::eq(
            motion.solid().materials()[0].law(),
            hot[0].law()
        ));
        for v in motion.velocity_mps() {
            close(v[0], velocity[0] + change, 2e-8);
            close(v[1], velocity[1], 2e-8);
            close(v[2], velocity[2], 2e-8);
        }
        for record in &accepted.mechanical.records {
            for energy in record.energy {
                assert_eq!(energy.inertia.reference_mass_kg, old.reference_mass_kg());
            }
        }
        close(
            accepted.mechanical.energy.inertia.numerical_loss_j,
            numerical_loss,
            2e-7,
        );
        close(
            accepted.mechanical.energy.inertia.kinetic_change_j,
            mass * (velocity[0] * change + 0.5 * change * change),
            2e-8,
        );
        assert_eq!(accepted.transfer.parameter_work_j, 0.0);
        assert_eq!(accepted.transfer.fracture_dissipation_j, 0.0);
        assert!(accepted.total_balance_error_j.abs() < 1e-10);
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[velocity; 4]);
        assert_eq!(
            old.solid().materials()[0].material().state_coordinate("T"),
            Some(300.0)
        );
    });
}
