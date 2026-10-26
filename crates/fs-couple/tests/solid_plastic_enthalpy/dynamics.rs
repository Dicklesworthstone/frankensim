//! Source inventory, physical self-heating and continuing inertial trajectories.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

fn motion<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    velocity: Vec<[f64; 3]>,
    cx: &Cx<'_>,
) -> ThermoplasticDynamics<'a> {
    ThermoplasticDynamics::try_new(bind(thermal, &[fixture], &NODES, &TETS, cx), velocity, cx)
        .unwrap()
}

#[test]
fn g1_dynamic_rigid_translation_heats_the_original_material_without_artificial_loss() {
    context(|cx| {
        let cp = 500.0;
        let fixture = Fixture::new(cp);
        let duration = 0.001;
        let supplied_heat = RHO * V * cp * 3.0;
        let thermal = initial(&[&fixture], &[], &[supplied_heat / duration]);
        let dynamics = motion(&fixture, &thermal, vec![[0.1, 0.0, 0.0]; 4], cx);
        let step = dynamics
            .advance(
                &thermal,
                &(),
                duration,
                &[],
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(
            step.thermal().temperature_k(CellId(100)).unwrap(),
            303.0,
            1e-8,
        );
        close(step.report().inertia.reference_mass_kg, RHO * V, 1e-18);
        for (u, v) in step
            .dynamics()
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(step.dynamics().velocity_mps())
        {
            close(u[0], 0.1 * duration, 1e-12);
            close(v[0], 0.1, 1e-9);
        }
        close(step.report().inertia.kinetic_change_j, 0.0, 1e-12);
        close(step.report().material.physical_heat_j, 0.0, 1e-12);
        close(
            step.report().material.combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert_eq!(dynamics.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn g1_dynamic_material_density_swap_controls_consistent_inertia_and_response() {
    context(|cx| {
        let dt = 2e-6;
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let k = E * 0.8 / (1.2 * 0.6) * V / (L * L);
        let mut displacement = Vec::new();
        for density in [RHO, 2.0 * RHO] {
            let fixture = Fixture::with_density_and_constitutive_claims(
                500.0,
                None,
                scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
                scalar(BASIS.property_name(), Pressure::DIMS, H),
                density,
            );
            let thermal = initial(&[&fixture], &[], &[0.0]);
            let mut velocity = vec![[0.0; 3]; 4];
            velocity[1][0] = 1.0;
            let dynamics = motion(&fixture, &thermal, velocity, cx);
            let step = dynamics
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let mass = density * V / 10.0;
            let expected = mass * dt / (mass + k * dt * dt);
            let u = step.dynamics().solid().solution().displacement_m[1][0];
            close(u, expected, 1e-12);
            close(step.report().inertia.reference_mass_kg, density * V, 1e-18);
            close(
                step.thermal().temperature_k(CellId(100)).unwrap(),
                300.0,
                1e-8,
            );
            close(step.report().material.physical_heat_j, 0.0, 1e-12);
            assert!(step.report().inertia.numerical_loss_j > 0.0);
            assert!(step.report().material.numerical_loss_j > 0.0);
            close(
                step.report().material.combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            displacement.push(u);
        }
        assert!(displacement[1] > displacement[0]);
    });
}

#[test]
fn g1_dynamic_plastic_loss_heats_once_and_continues_accepted_motion() {
    context(|cx| {
        let fixture = Fixture::new(5.0);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let mut velocity = vec![[0.0; 3]; 4];
        velocity[1][0] = 30.0;
        let dynamics = motion(&fixture, &thermal, velocity, cx);
        let dt = 2e-6;
        let first = dynamics
            .advance(
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let report = first.report();
        assert!(report.material.physical_heat_j > 0.0);
        assert!(report.inertia.kinetic_change_j < 0.0);
        assert!(report.inertia.numerical_loss_j > 0.0);
        assert!(report.material.numerical_loss_j > 0.0);
        close(
            report.material.thermal_audit.stored_energy_change_j,
            report.material.physical_heat_j,
            config().energy_tolerance_j,
        );
        close(
            report.material.combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert!(first.thermal().temperature_k(CellId(100)).unwrap() > 300.0);
        let previous_u = first.dynamics().solid().solution().displacement_m.clone();
        let previous_history = first.dynamics().solid().solution().element_states.clone();
        let previous_kinetic = first.report().inertia.kinetic_after_j;
        let (warm, moving, observer, _) = first.into_parts();
        let second = moving
            .advance(
                &warm,
                &observer,
                2.0 * dt,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(
            second.report().inertia.kinetic_before_j,
            previous_kinetic,
            1e-12,
        );
        assert_ne!(
            second.dynamics().solid().solution().displacement_m,
            previous_u
        );
        assert_eq!(moving.solid().solution().element_states, previous_history);
        close(
            second.report().material.combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert_eq!(moving.time_s(), dt);
        assert_eq!(dynamics.time_s(), 0.0);
    });
}

fn affine_modulus(slope: f64) -> PropertyClaim {
    let mut claim = scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, E - 20.0 * slope), (320.0, E + 20.0 * slope)],
        dims: Pressure::DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    claim
}

#[test]
fn g1_dynamic_heating_refreshes_stiffness_and_reciprocal_heat_from_source_data() {
    context(|cx| {
        let dt = 2e-6;
        let cp = 500.0;
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let mut displacement = Vec::new();
        for slope in [0.0, -1e7] {
            let fixture = Fixture::with_constitutive_claims(
                cp,
                None,
                affine_modulus(slope),
                scalar(BASIS.property_name(), Pressure::DIMS, H),
            );
            let thermal = initial(&[&fixture], &[], &[RHO * V * cp * 6.0 / dt]);
            let mut velocity = vec![[0.0; 3]; 4];
            velocity[1][0] = 1.0;
            let dynamics = motion(&fixture, &thermal, velocity, cx);
            let step = dynamics
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert!(step.thermal().temperature_k(CellId(100)).unwrap() > 305.9);
            close(step.report().inertia.reference_mass_kg, RHO * V, 1e-18);
            close(step.report().material.physical_heat_j, 0.0, 1e-12);
            close(
                step.report().material.combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            if slope != 0.0 {
                assert!(step.report().material.reversible_heat_j.abs() > 0.0);
            }
            displacement.push(step.dynamics().solid().solution().displacement_m[1][0]);
        }
        assert!(displacement[1] > displacement[0]);
    });
}

#[test]
fn g4_dynamic_source_exit_and_cancel_preserve_complete_accepted_state() {
    context(|cx| {
        let fixture = Fixture::new(500.0);
        let dt = 0.001;
        let thermal = initial(&[&fixture], &[], &[RHO * V * 500.0 * 25.0 / dt]);
        let dynamics = motion(&fixture, &thermal, vec![[0.0; 3]; 4], cx);
        let old = dynamics.solid().solution().clone();
        let fixed: Vec<_> = (0..12).map(|d| (d, 0.0)).collect();
        assert!(
            dynamics
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                )
                .is_err()
        );
        assert_eq!(dynamics.solid().solution(), &old);
        assert_eq!(dynamics.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(dynamics.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        let gate = CancelGate::new();
        gate.request();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                dynamics
                    .advance(
                        &thermal,
                        &(),
                        dt,
                        &fixed,
                        &[],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        &cancelled
                    )
                    .is_err()
            );
        });
        assert_eq!(dynamics.solid().solution(), &old);
    });
}
