//! G1/G3/G4: actual native inertia and original axisymmetric surface work.
use super::*;
use fs_contact::capillary_bridge::axisymmetric::load_line::AxisymmetricLoadLineConfig;
use fs_contact::capillary_bridge::axisymmetric::*;
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    AxisymmetricCapillaryModalConfig, AxisymmetricCapillaryModalSystem,
};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection};

const MASS: [f64; 2] = [2e-4, 3e-4];
const SAMPLE_RATE: u32 = 50_000;

fn numerical() -> AxisymmetricNumerics {
    AxisymmetricNumerics {
        quadrature_relative_tolerance: 8e-14,
        gap_relative_tolerance: 3e-14,
        maximum_quadrature_panels: 4096,
        maximum_root_iterations: 80,
    }
}
fn model3d(sigma: f64) -> SourceAxisymmetricCapillaryBridge {
    let (liquid, lower) = sources(1.1, sigma, 300.0, true);
    let (_, upper) = sources(1.3, sigma, 300.0, true);
    bridge3d(liquid, [lower[0].clone(), upper[1].clone()])
}
fn bridge3d(
    liquid: ResolvedLiquidSurfaceState,
    wetting: [ResolvedLiquidWettingState; 2],
) -> SourceAxisymmetricCapillaryBridge {
    SourceAxisymmetricCapillaryBridge::new(
        liquid,
        wetting,
        VOLUME,
        AxisymmetricBranch {
            minimum_suction_parameter: 0.0,
            maximum_suction_parameter: 32.0,
            minimum_stability_margin: 0.01,
        },
        AxisymmetricLimits {
            minimum_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_wetted_radius_m: 0.1,
            maximum_suction_pa: 1e8,
            maximum_abs_force_n: 10.0,
            maximum_abs_surface_energy_j: 1.0,
            maximum_abs_step_work_j: 1.0,
            gravity_m_s2: 0.0,
            maximum_bond_number: 0.01,
        },
    )
    .unwrap()
}
fn controls() -> AxisymmetricCapillaryModalConfig {
    AxisymmetricCapillaryModalConfig {
        contact: AxisymmetricLoadLineConfig {
            mechanical: BridgeLoadLineConfig {
                max_iterations: 32,
                max_evaluations: 161,
                max_backtracks: 16,
                gap_absolute_tolerance_m: 3e-15,
                gap_relative_tolerance: 1e-10,
                maximum_gap_increment_m: 1e-4,
                minimum_compatibility_slope: 0.05,
            },
            work: AxisymmetricWorkBudget {
                equilibrium: numerical(),
                absolute_tolerance_j: 1e-18,
                relative_tolerance: 1e-5,
            },
            force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-6,
        },
        force_absolute_tolerance_n: 3e-10,
        force_relative_tolerance: 3e-6,
    }
}
fn native_network(connection: bool) -> CoupledModalSystem {
    let models = MASS
        .into_iter()
        .map(|_| {
            ModalAcousticTimeModel::try_new_with_free_coordinates(
                SAMPLE_RATE,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect();
    let connections = if connection {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / det::sqrt(MASS[0])],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0 / det::sqrt(MASS[1])],
            },
            stiffness_n_m: 0.1,
            damping_n_s_m: 1e-4,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        models,
        connections,
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 1,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 0.01,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 5e-18,
            energy_relative_tolerance: 2e-10,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn work_column() -> Vec<f64> {
    vec![-1.0 / det::sqrt(MASS[0]), 1.0 / det::sqrt(MASS[1])]
}
fn materials(source: &SourceAxisymmetricCapillaryBridge) -> [&MaterialStateId; 2] {
    [
        &source.wetting()[0].surfaces()[0].material,
        &source.wetting()[1].surfaces()[0].material,
    ]
}
fn owner(
    source: &SourceAxisymmetricCapillaryBridge,
    connection: bool,
) -> AxisymmetricCapillaryModalSystem<'_> {
    let gate = CancelGate::new_clock_free();
    let state = source.solve_at_suction(2.0, numerical(), &gate).unwrap();
    AxisymmetricCapillaryModalSystem::new_axisymmetric(
        native_network(connection),
        source,
        materials(source),
        state,
        work_column(),
        controls(),
        &gate,
    )
    .unwrap()
}

#[test]
fn g3_g4_thermal_expansion_refresh_keeps_native_motion_clock_and_observation() {
    let (original, lower) = source_cards(1.1, 0.05, true);
    let (_, upper) = source_cards(1.3, 0.05, true);
    let curves = |original: &ClaimSet| {
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims_ordered() {
            let factor = match original.key.name() {
                "density" => 0.95,
                SURFACE_TENSION_PROPERTY | EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY => 0.8,
                _ => 1.0,
            };
            let mut changed = original.clone();
            let PropertyValue::Scalar { value, dims } = original.value else {
                unreachable!()
            };
            changed.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(300.0, value), (320.0, factor * value)],
                dims,
            };
            changed.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(changed).unwrap();
        }
        claims
    };
    let material =
        MaterialCard::assemble(original.id().clone(), curves(original.claims()), vec![]).unwrap();
    let surfaces = [lower[0].clone(), upper[1].clone()].map(|card| {
        InterfaceSystemCard::assemble(
            card.surface_a().clone(),
            card.surface_b().clone(),
            card.context().clone(),
            curves(card.claims()),
            vec![],
        )
        .unwrap()
    });
    let (liquid, wetting) = resolved_sources(&material, &surfaces, 300.0, true);
    let source = bridge3d(liquid, wetting);
    let gate = CancelGate::new_clock_free();
    let mut system = owner(&source, true);
    let initial_state = system.state();
    let mut control = owner(&source, true);
    for _ in 0..3 {
        system.step(&[0.0; 2], 0.0).unwrap();
        control.step(&[0.0; 2], 0.0).unwrap();
    }
    let state = system.state();
    let gap = system.gap_m();
    let energy = system.total_energy_j().unwrap();
    let frame = *system.last_frame().unwrap();
    let bodies: Vec<_> = system
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    let hot = QueryPoint::new().with("T", 320.0).unwrap();
    let stale = source
        .refresh_at_query_point(
            initial_state,
            &material,
            surfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedMass,
            controls().contact.work,
            &gate,
        )
        .unwrap();
    let refresh = source
        .refresh_at_query_point(
            state,
            &material,
            surfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedMass,
            controls().contact.work,
            &gate,
        )
        .unwrap();
    assert!(system.apply_source_refresh(&stale, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(system.apply_source_refresh(&refresh, &cancelled).is_err());
    assert_eq!(system.state(), state);
    assert_eq!(system.total_energy_j().unwrap(), energy);
    assert_eq!(system.samples_rendered(), 3);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(
        system
            .components()
            .iter()
            .map(|m| m.states().to_vec())
            .collect::<Vec<_>>(),
        bodies
    );

    let work = system.apply_source_refresh(&refresh, &gate).unwrap();
    close(system.total_energy_j().unwrap() - energy, work, 1e-21);
    close(
        system.source().liquid().density_kg_m3() * system.state().profile().liquid_volume_m3,
        source.liquid().density_kg_m3() * VOLUME,
        1e-21,
    );
    assert_eq!(system.gap_m(), gap);
    assert_eq!(system.samples_rendered(), 3);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(
        system
            .components()
            .iter()
            .map(|m| m.states().to_vec())
            .collect::<Vec<_>>(),
        bodies
    );
    assert!(system.apply_source_refresh(&refresh, &gate).is_err());
    let next = *system.step(&[0.0; 2], 0.0).unwrap();
    let unchanged = *control.step(&[0.0; 2], 0.0).unwrap();
    assert_ne!(next.applied_mean_force_n, unchanged.applied_mean_force_n);
    let mut loss = next.network_dissipation_j;
    for _ in 0..10 {
        loss += system.step(&[0.0; 2], 0.0).unwrap().network_dissipation_j;
    }
    close(
        system.total_energy_j().unwrap() + loss,
        energy + work,
        1e-15,
    );
}

#[test]
fn g1_original_native_bodies_exchange_3d_surface_work_without_geometric_drift() {
    let source = model3d(0.05);
    let mut system = owner(&source, false);
    let initial_gap = system.gap_m();
    let initial_energy = system.total_energy_j().unwrap();
    let dt = system.sample_period_s();
    for _ in 0..200 {
        let old = [
            system.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            system.components()[1].states()[0].velocity_m_sqrt_kg_per_s,
        ];
        let frame = *system.step(&[0.0; 2], 0.0).unwrap();
        let a = system.components()[0].states()[0];
        let b = system.components()[1].states()[0];
        close(
            det::sqrt(MASS[0]) * a.velocity_m_sqrt_kg_per_s
                + det::sqrt(MASS[1]) * b.velocity_m_sqrt_kg_per_s,
            0.0,
            3e-18,
        );
        close(
            det::sqrt(MASS[0]) * (a.velocity_m_sqrt_kg_per_s - old[0]),
            -dt * frame.applied_mean_force_n,
            3e-18,
        );
        close(
            det::sqrt(MASS[1]) * (b.velocity_m_sqrt_kg_per_s - old[1]),
            dt * frame.applied_mean_force_n,
            3e-18,
        );
        let physical_gap = initial_gap - a.displacement_m_sqrt_kg / det::sqrt(MASS[0])
            + b.displacement_m_sqrt_kg / det::sqrt(MASS[1]);
        close(system.gap_m(), physical_gap, 1e-17);
        close(
            system.state().profile().gap_m,
            physical_gap,
            frame.compatibility_tolerance_m,
        );
        close(system.total_energy_j().unwrap(), initial_energy, 1e-15);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(
            frame.constitutive_residual_n.abs() + frame.estimated_constitutive_force_error_n
                <= frame.force_tolerance_n
        );
        assert!(
            frame.compatibility_residual_m.abs() + frame.estimated_compatibility_error_m
                <= frame.compatibility_tolerance_m
        );
        assert!(frame.evaluations <= controls().contact.mechanical.max_evaluations);
        assert_eq!(frame.topology_release_j, 0.0);
        assert_eq!(frame.bridge.state().profile().liquid_volume_m3, VOLUME);
    }
    assert!(system.gap_m() < 0.95 * initial_gap);
    assert!(system.last_frame().unwrap().network_energy_j > 0.0);
    assert!(std::ptr::eq(system.source(), &source));
}

#[test]
fn g1_native_connections_external_load_and_prescribed_motion_share_one_work_account() {
    let source = model3d(0.05);
    let mut system = owner(&source, true);
    let initial = system.total_energy_j().unwrap();
    let force = [3e-6 / det::sqrt(MASS[0]), 0.0];
    let prescribed_rate = 1e-5;
    let mut external = 0.0;
    let mut prescribed = 0.0;
    let mut dissipation = 0.0;
    for i in 1..=100 {
        let frame = *system.step(&force, prescribed_rate).unwrap();
        external += frame.external_work_j;
        prescribed += frame.prescribed_work_j;
        dissipation += frame.network_dissipation_j;
        assert!(frame.network_dissipation_j >= 0.0);
        close(
            system.total_energy_j().unwrap() - initial + dissipation,
            external + prescribed,
            1e-15,
        );
        let momentum: f64 = MASS
            .iter()
            .zip(system.components())
            .map(|(mass, body)| det::sqrt(*mass) * body.states()[0].velocity_m_sqrt_kg_per_s)
            .sum();
        close(
            momentum,
            3e-6 * f64::from(i) * system.sample_period_s(),
            3e-18,
        );
        close(
            frame.prescribed_work_j,
            -frame.applied_mean_force_n * prescribed_rate * system.sample_period_s(),
            1e-24,
        );
    }
    assert!(dissipation > 0.0);
    assert!(external > 0.0);
    assert!(prescribed > 0.0);
}

#[test]
fn g1_tiny_prescribed_3d_capillary_motion_carries_roundoff_separately_from_force_and_heat() {
    let source = model3d(0.05);
    let gate = CancelGate::new_clock_free();
    let state = source.solve_at_suction(2.0, numerical(), &gate).unwrap();
    let mut system = AxisymmetricCapillaryModalSystem::new_axisymmetric(
        native_network(false),
        &source,
        materials(&source),
        state,
        vec![0.0; 2],
        controls(),
        &gate,
    )
    .unwrap();
    let initial = system.total_energy_j().unwrap();
    let rate = 1e-9 * system.gap_m() / system.sample_period_s();
    let mut work = 0.0;
    for direction in [1.0, 1.0, -1.0, -1.0] {
        let frame = *system.step(&[0.0; 2], direction * rate).unwrap();
        work += frame.prescribed_work_j;
        assert!(frame.estimated_constitutive_force_error_n < 1e-3 * frame.force_tolerance_n);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert_eq!(frame.topology_release_j, 0.0);
        assert_eq!(frame.network_dissipation_j, 0.0);
        assert_eq!(frame.network_energy_j, 0.0);
        close(system.total_energy_j().unwrap() - initial, work, 1e-20);
    }
}

#[test]
fn g0_g4_source_accuracy_and_cancel_refusals_preserve_native_state_and_retry() {
    let source = model3d(0.05);
    let other = model3d(0.06);
    let gate = CancelGate::new_clock_free();
    let state = source.solve_at_suction(2.0, numerical(), &gate).unwrap();
    let wrong = other.solve_at_suction(2.0, numerical(), &gate).unwrap();
    assert!(
        AxisymmetricCapillaryModalSystem::new_axisymmetric(
            native_network(false),
            &source,
            materials(&source),
            wrong,
            work_column(),
            controls(),
            &gate
        )
        .is_err()
    );
    let assigned = materials(&source);
    assert!(
        AxisymmetricCapillaryModalSystem::new_axisymmetric(
            native_network(false),
            &source,
            [assigned[1], assigned[0]],
            state,
            work_column(),
            controls(),
            &gate
        )
        .is_err()
    );
    let make = || {
        AxisymmetricCapillaryModalSystem::new_axisymmetric(
            native_network(false),
            &source,
            materials(&source),
            state,
            vec![0.0; 2],
            controls(),
            &gate,
        )
        .unwrap()
    };
    let mut system = make();
    let mut reference = make();
    system.step(&[0.0; 2], 0.0).unwrap();
    reference.step(&[0.0; 2], 0.0).unwrap();
    assert_eq!(system.state(), state);
    let accepted = *system.last_frame().unwrap();
    let before: Vec<_> = system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect();
    let rate = (system.gap_m().next_up() - system.gap_m()) / system.sample_period_s();
    assert_eq!(
        system.gap_m() + rate * system.sample_period_s(),
        system.gap_m().next_up()
    );
    assert!(system.step(&[0.0; 2], rate).is_err());
    assert!(system.step(&[0.0; 2], -1e6).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(system.step_under_gate(&[0.0; 2], 0.0, &cancelled).is_err());
    assert_eq!(system.state(), state);
    assert_eq!(system.gap_m(), state.profile().gap_m);
    assert_eq!(*system.last_frame().unwrap(), accepted);
    assert_eq!(
        system
            .components()
            .iter()
            .flat_map(|m| m.states().iter().copied())
            .collect::<Vec<_>>(),
        before
    );
    assert_eq!(
        system.step(&[0.0; 2], -1e-4).unwrap(),
        reference.step(&[0.0; 2], -1e-4).unwrap()
    );
    let mut budget = controls();
    budget.contact.mechanical.max_evaluations = 2;
    let mut limited = AxisymmetricCapillaryModalSystem::new_axisymmetric(
        native_network(false),
        &source,
        materials(&source),
        state,
        work_column(),
        budget,
        &gate,
    )
    .unwrap();
    assert!(limited.step(&[0.0; 2], 0.0).is_err());
    assert_eq!(limited.samples_rendered(), 0);
    assert_eq!(limited.state(), state);
    assert!(limited.last_frame().is_none());
}
