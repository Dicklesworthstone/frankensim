//! Independent endpoint transport and retained-history numerical oracles.
use super::*;
use fs_lbm::material::conjugate::buoyant::MaterialEndpointConfig;

fn endpoint_controls() -> MaterialEndpointConfig {
    MaterialEndpointConfig {
        maximum_passes_per_step: 80,
        maximum_momentum_solves: 1000,
        maximum_energy_solves: 1000,
        maximum_source_refreshes: 2001,
        temperature_tolerance_k: 1e-10,
        transport_relative_tolerance: 1e-11,
        force_tolerance_n_m3: 1e-11,
        temperature_relaxation: 1.0,
    }
}

#[test]
fn nonlinear_endpoint_conductivity_solves_the_implicit_heat_equation() {
    // C=2 J/K and the six half-cell paths give 12*k(T_new) W/K.
    // The quadratic below is independent of the native Picard iteration.
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut setup = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    setup.power_w = vec![12.0];
    let unsteady = UnsteadyConfig::new(0.05, 8);
    let source = |_: usize, t: f64, _: &FluidTransportState| {
        FluidTransportState::try_new(
            1.0 - 0.001 * (t - 300.0).powi(2),
            2.0,
            1.0 + 0.1 * (t - 300.0),
            1.0 / (1.0 + 0.05 * (t - 300.0)),
        )
    };
    let explicit = fluid
        .march_natural_convection(
            &domain,
            &[],
            &flow,
            &setup,
            &EnergyConfig::default(),
            &[300.0],
            limits([0.0; 3]),
            &unsteady,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            source,
        )
        .unwrap();
    for relaxation in [1.0, 0.6] {
        let mut controls = endpoint_controls();
        controls.temperature_relaxation = relaxation;
        let mut queries = Vec::new();
        let result = fluid
            .march_natural_convection_implicit(
                &domain,
                &[],
                &flow,
                &setup,
                &EnergyConfig::default(),
                &[300.0],
                limits([0.0; 3]),
                &unsteady,
                controls,
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |cell, t, reference| {
                    queries.push(t);
                    source(cell, t, reference)
                },
            )
            .unwrap();
        let solution = &result.solution;
        let mut theta = 0.0;
        for (step, record) in solution.march.energy_records.iter().enumerate() {
            let old = theta;
            let b = 2.0 / unsteady.time_step_s + 12.0;
            let c = 2.0 / unsteady.time_step_s * old + 12.0;
            theta = 2.0 * c / (b + (b * b + 4.8 * c).sqrt());
            close(record.max_temperature_k, 300.0 + theta, 1e-11);
            close(record.stored_energy_change_j, 2.0 * (theta - old), 1e-8);
            close(record.source_j, 0.6, 1e-12);
            assert!(record.closure_j.abs() < 1e-8);
            assert!(result.records[step].passes > 1);
            assert!(result.records[step].temperature_defect_k <= controls.temperature_tolerance_k);
            assert!(
                result.records[step].transport_relative_defect
                    <= controls.transport_relative_tolerance
            );
        }
        close(solution.march.temperature[0], 300.0 + theta, 1e-11);
        assert!((solution.march.temperature[0] - explicit.march.temperature[0]).abs() > 1e-4);
        assert_eq!(result.source_refreshes, queries.len());
        assert_eq!(
            result.momentum_solves,
            result.records.iter().map(|r| r.passes).sum::<usize>()
        );
        assert_eq!(result.energy_solves, result.momentum_solves);
        assert_eq!(solution.march.flow.records.len(), unsteady.steps);
        assert_eq!(result.records.len(), unsteady.steps);
        close(solution.applied_time_s, 0.4, 1e-14);
        close(solution.endpoint_time_s, 0.4, 1e-14);
        assert_eq!(
            solution.endpoint[0].temperature_k,
            solution.march.temperature[0]
        );
        assert!(
            (solution.applied[0].temperature_k - solution.march.temperature[0]).abs()
                <= controls.temperature_tolerance_k
        );
    }
}

fn darcy_column() -> (VoxelDomain, SimpleConfig, ThermalSetup) {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(faces);
    flow.pressure_tolerance = 1e-12;
    flow.resistances.push(FlowResistance::Volume {
        lo: [0; 3],
        hi: [4, 1, 1],
        permeability_m2: [1.0; 3],
        inertial_per_m: [0.0; 3],
    });
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let mut thermal = ThermalSetup::new(faces);
    thermal.fixed_temperature = (0..4).map(|cell| (cell, 310.0)).collect();
    (domain, flow, thermal)
}

#[test]
fn endpoint_buoyancy_drives_the_first_step_without_advancing_be_or_bdf_history_twice() {
    let (domain, flow, thermal) = darcy_column();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        for expansion in [-0.001, 0.001] {
            let mut unsteady = UnsteadyConfig::new(0.25, 6);
            unsteady.scheme = scheme;
            unsteady.inner_tolerance = 1e-11;
            unsteady.inner_iterations = 1000;
            unsteady.average_from_step = 1;
            unsteady.probe = Some(1);
            let result = fluid
                .march_natural_convection_implicit(
                    &domain,
                    &[],
                    &flow,
                    &thermal,
                    &EnergyConfig::default(),
                    &[300.0; 4],
                    limits([-9.81, 0.0, 0.0]),
                    &unsteady,
                    endpoint_controls(),
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| {
                        FluidTransportState::try_new(
                            1.0 - expansion * (t - 300.0).powi(2) / 10.0,
                            2.0,
                            1.0 + 0.01 * (t - 300.0),
                            1.0 + 0.1 * (t - 300.0),
                        )
                    },
                )
                .unwrap();
            let force = 9.81 * expansion * 10.0;
            let mu = 2.0;
            let march = &result.solution.march;
            let (mut previous, mut older, mut sum) = (0.0, 0.0, 0.0);
            for (step, record) in march.flow.records.iter().enumerate() {
                let expected = if step == 0 || scheme == TimeScheme::BackwardEuler {
                    (previous + unsteady.time_step_s * force) / (1.0 + unsteady.time_step_s * mu)
                } else {
                    (4.0 * previous - older + 2.0 * unsteady.time_step_s * force)
                        / (3.0 + 2.0 * unsteady.time_step_s * mu)
                };
                close(record.probe_velocity_m_s[0], expected, 2e-8);
                older = previous;
                previous = expected;
                sum += expected;
            }
            assert!(march.flow.records[0].probe_velocity_m_s[0].abs() > 1e-3);
            assert_eq!(result.records[0].passes, 2);
            assert_eq!(result.momentum_solves, unsteady.steps + 1);
            assert_eq!(march.flow.records.len(), unsteady.steps);
            assert_eq!(march.energy_records.len(), unsteady.steps);
            assert_eq!(march.flow.averaged_steps, unsteady.steps);
            close(
                march.flow.mean_velocity_m_s[1][0],
                sum / unsteady.steps as f64,
                2e-8,
            );
            close(
                march.flow.mean_field.flux_x(2, 0, 0),
                sum / unsteady.steps as f64,
                2e-8,
            );
            close(march.flow.flow.field.flux_x(2, 0, 0), previous, 2e-8);
            close(march.energy_records[0].stored_energy_change_j, 80.0, 1e-10);
            for heat in &march.energy_records {
                assert!(heat.closure_j.abs() < 1e-8);
            }
            for source in &result.solution.applied {
                close(source.force_n_m3[0], force, 1e-12);
                close(
                    source
                        .material
                        .resolved_transport()
                        .dynamic_viscosity_pa_s(),
                    mu,
                    1e-12,
                );
            }
        }
    }
}

#[test]
fn independent_endpoint_force_and_transport_checks_cannot_be_hidden_by_temperature_tolerance() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 300.000_001)];
    let saved = thermal.clone();
    let mut endpoint = endpoint_controls();
    endpoint.maximum_passes_per_step = 1;
    endpoint.temperature_tolerance_k = 1.0;
    for force_only in [false, true] {
        let result = fluid.march_natural_convection_implicit(
            &domain,
            &[],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &EnergyConfig::default(),
            &[300.0],
            limits([0.0, -9.81, 0.0]),
            &UnsteadyConfig::new(0.1, 1),
            endpoint,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| {
                let fraction = ((t - 300.0) / 1e-6).clamp(0.0, 1.0);
                FluidTransportState::try_new(
                    if force_only {
                        1.0 - 0.01 * fraction
                    } else {
                        1.0
                    },
                    2.0,
                    if force_only { 1.0 } else { 1.0 + fraction },
                    1.0,
                )
            },
        );
        assert!(matches!(result,
            Err(MaterialBuoyancyError::Solver(ChtError::SolverNotConverged {
                system: "implicit transient material coupling", relative_residual, ..
            })) if relative_residual > 1.0));
    }
    assert_eq!(thermal, saved);
}

#[test]
fn discarded_trials_consume_each_budget_and_late_source_refusal_publishes_no_prefix() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.power_w = vec![2.0];
    let initial = [300.0];
    let saved = thermal.clone();
    let unsteady = UnsteadyConfig::new(0.1, 2);
    for kind in 0..3 {
        let mut controls = endpoint_controls();
        match kind {
            0 => controls.maximum_momentum_solves = 2,
            1 => controls.maximum_energy_solves = 2,
            _ => controls.maximum_source_refreshes = 3,
        }
        let result = fluid.march_natural_convection_implicit(
            &domain,
            &[],
            &flow,
            &thermal,
            &EnergyConfig::default(),
            &initial,
            limits([0.0; 3]),
            &unsteady,
            controls,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, _, _| Ok(reference),
        );
        assert!(matches!(
            result,
            Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput { .. }))
        ));
    }
    let mut queries = Vec::new();
    let refused = fluid.march_natural_convection_implicit(
        &domain,
        &[],
        &flow,
        &thermal,
        &EnergyConfig::default(),
        &initial,
        limits([0.0; 3]),
        &unsteady,
        endpoint_controls(),
        |_| 1.0,
        |_| 1.0,
        None,
        &CancelGate::new(),
        |_, t, _| {
            queries.push(t);
            if t > 300.15 {
                Err(MaterialError::State {
                    what: "original source support exhausted".into(),
                })
            } else {
                Ok(reference)
            }
        },
    );
    assert!(
        matches!(refused, Err(MaterialBuoyancyError::Material { temperature_k, .. }) if temperature_k > 300.15)
    );
    assert!(
        queries
            .iter()
            .filter(|&&t| (t - 300.1).abs() < 1e-9)
            .count()
            >= 2
    );
    assert_eq!(queries.iter().filter(|&&t| t > 300.15).count(), 1);
    let gate = CancelGate::new();
    gate.request();
    let mut calls = 0;
    let cancelled = fluid.march_natural_convection_implicit(
        &domain,
        &[],
        &flow,
        &thermal,
        &EnergyConfig::default(),
        &initial,
        limits([0.0; 3]),
        &unsteady,
        endpoint_controls(),
        |_| 1.0,
        |_| 1.0,
        None,
        &gate,
        |_, _, _| {
            calls += 1;
            Ok(reference)
        },
    );
    assert!(matches!(
        cancelled,
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));
    assert_eq!(calls, 0);
    assert_eq!(initial, [300.0]);
    assert_eq!(thermal, saved);
}
