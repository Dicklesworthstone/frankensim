//! Coupled native momentum and original source-caloric storage regressions.
use super::*;
use fs_lbm::material::conjugate::buoyant::MaterialEndpointConfig;

fn endpoint_controls() -> MaterialEndpointConfig {
    MaterialEndpointConfig {
        maximum_passes_per_step: 80,
        maximum_momentum_solves: 1000,
        maximum_energy_solves: 2000,
        maximum_source_refreshes: 2001,
        temperature_tolerance_k: 1e-9,
        transport_relative_tolerance: 1e-10,
        force_tolerance_n_m3: 1e-11,
        temperature_relaxation: 1.0,
    }
}

#[test]
fn implicit_caloric_feedback_retains_true_storage_and_single_be_bdf_time_levels() {
    let card = variable_card(0.04, 0.001, 300.0, 340.0);
    let law = liquid(&card);
    let chart = caloric(&card, 1e9, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
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
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let mut limits = controls().boussinesq;
    limits.gravity_m_s2 = [-9.81, 0.0, 0.0];
    limits.max_relative_density_departure = 0.05;
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let mut unsteady = UnsteadyConfig::new(0.25, 4);
        unsteady.scheme = scheme;
        unsteady.inner_tolerance = 1e-11;
        unsteady.inner_iterations = 1000;
        unsteady.average_from_step = 1;
        unsteady.probe = Some(1);
        let result = bound
            .march_natural_convection_with_enthalpy_implicit(
                &domain,
                &[],
                &flow,
                &thermal,
                &energy,
                &[300.0; 4],
                limits,
                &chart,
                transient_controls(),
                &unsteady,
                endpoint_controls(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, t, _| {
                    law.at_temperature(Temperature::new(t))
                        .map_err(|error| MaterialError::State {
                            what: error.to_string(),
                        })
                },
            )
            .unwrap();
        let solution = &result.solution;
        let march = &solution.convection.march;
        // Four reference masses each gain the integrated 22 J/kg, rather
        // than reference-Cp*10. Repeated coupling trials never add storage.
        near(march.energy_records[0].stored_energy_change_j, 88.0, 1e-9);
        near(march.energy_records[0].boundary_outflow_j, -88.0, 1e-9);
        let (mut previous, mut older, mut sum) = (0.0, 0.0, 0.0);
        for (step, record) in march.flow.records.iter().enumerate() {
            let expected = if step == 0 || scheme == TimeScheme::BackwardEuler {
                (previous + 0.25 * 0.0981) / 1.25
            } else {
                (4.0 * previous - older + 0.5 * 0.0981) / 3.5
            };
            near(record.probe_velocity_m_s[0], expected, 2e-8);
            older = previous;
            previous = expected;
            sum += expected;
        }
        near(
            march.flow.mean_velocity_m_s[1][0],
            sum / unsteady.steps as f64,
            2e-8,
        );
        assert_eq!(march.flow.averaged_steps, unsteady.steps);
        assert_eq!(march.flow.records.len(), unsteady.steps);
        assert_eq!(march.energy_records.len(), unsteady.steps);
        assert_eq!(solution.records.len(), unsteady.steps);
        assert_eq!(result.records.len(), unsteady.steps);
        assert_eq!(result.records[0].passes, 2);
        assert!(solution.records[0].passes >= 4);
        assert_eq!(
            result.energy_solves,
            solution
                .records
                .iter()
                .map(|r| r.energy_solves)
                .sum::<usize>()
        );
        assert!(result.energy_solves > result.momentum_solves);
        assert!(core::ptr::eq(solution.caloric, &chart));
        for record in &march.energy_records {
            assert!(record.closure_j.abs() <= 1e-8);
        }
        for sample in &solution.convection.applied {
            near(sample.force_n_m3[0], 0.0981, 1e-12);
            near(sample.temperature_k, 310.0, 1e-14);
            let original = sample.material.material().material();
            assert_eq!(original.card_identity(), card.content_hash());
        }
        near(
            solution.convection.applied_time_s,
            solution.convection.endpoint_time_s,
            1e-14,
        );
    }
}

#[test]
fn source_enthalpy_advection_and_outer_work_share_one_physical_heat_interval() {
    let card = variable_card(0.04, 0.001, 300.0, 340.0);
    let law = liquid(&card);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    for reverse in [false, true] {
        let (domain, flow, thermal, energy) = channel(reverse);
        let unsteady = UnsteadyConfig::new(0.5, 5);
        let mut limits = controls().boussinesq;
        limits.max_relative_density_departure = 0.05;
        let mut temperatures = Vec::new();
        for datum in [0.0, 1e9] {
            let chart = caloric(&card, datum, 340.0);
            let result = bound
                .march_natural_convection_with_enthalpy_implicit(
                    &domain,
                    &[],
                    &flow,
                    &thermal,
                    &energy,
                    &[310.0],
                    limits,
                    &chart,
                    transient_controls(),
                    &unsteady,
                    endpoint_controls(),
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| {
                        law.at_temperature(Temperature::new(t)).map_err(|error| {
                            MaterialError::State {
                                what: error.to_string(),
                            }
                        })
                    },
                )
                .unwrap();
            let mut h = 22.0;
            for record in &result.solution.convection.march.energy_records {
                let old = h;
                h = (h + unsteady.time_step_s * (2.6 + 0.1 * 22.0))
                    / (1.0 + unsteady.time_step_s * 0.1);
                assert!((record.max_temperature_k - temperature_from_h(h)).abs() < 1e-5);
                near(record.stored_energy_change_j, h - old, 1e-8);
                near(
                    record.boundary_outflow_j,
                    unsteady.time_step_s * 0.1 * (h - 22.0),
                    1e-8,
                );
                assert!(record.closure_j.abs() < 1e-8);
            }
            assert_eq!(
                result.energy_solves,
                result
                    .solution
                    .records
                    .iter()
                    .map(|r| r.energy_solves)
                    .sum::<usize>()
            );
            assert!(result.momentum_solves > unsteady.steps);
            temperatures.push(result.solution.convection.march.temperature[0]);
        }
        near(temperatures[0], temperatures[1], 1e-11);
    }
}

#[test]
fn every_caloric_and_radiation_solve_debits_both_coupled_and_original_budgets() {
    let card = card(0.04);
    let law = liquid(&card);
    let chart = caloric(&card, 1e9, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let dx = 0.01;
    let domain = VoxelDomain::from_fn(1, 1, 2, dx, |p| {
        if p[2] < dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let solids = [SolidMaterial::new("held plate", 200.0).with_heat_capacity(2.4e6)];
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 330.0)];
    let saved = thermal.clone();
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let radiation = RadiationConfig::new(vec![0.9], surroundings);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let unsteady = UnsteadyConfig::new(0.1, 1);
    let run = |outer, caloric_controls| {
        bound.march_natural_convection_with_enthalpy_implicit(
            &domain,
            &solids,
            &flow,
            &thermal,
            &energy,
            &[300.0; 2],
            controls().boussinesq,
            &chart,
            caloric_controls,
            &unsteady,
            outer,
            |_| 1.0,
            |_| 1.0,
            Some(&radiation),
            &CancelGate::new(),
            |_, t, _| {
                law.at_temperature(Temperature::new(t))
                    .map_err(|error| MaterialError::State {
                        what: error.to_string(),
                    })
            },
        )
    };
    let result = run(endpoint_controls(), transient_controls()).unwrap();
    assert_eq!(
        result.energy_solves,
        result.solution.records[0].energy_solves
    );
    assert!(result.energy_solves > result.solution.records[0].passes);
    let record = result.solution.convection.march.energy_records[0];
    let t = result.solution.convection.march.temperature[1];
    let stored = 2.4e6 * domain.volume(0) * 30.0
        + domain.volume(1) * chart.specific_enthalpy_change(300.0, t).unwrap();
    near(record.stored_energy_change_j, stored, 1e-12);
    assert!(record.closure_j.abs() < 1e-8);
    for outer_limited in [false, true] {
        let mut outer = endpoint_controls();
        let mut caloric_controls = transient_controls();
        if outer_limited {
            outer.maximum_energy_solves = result.energy_solves - 1;
        } else {
            caloric_controls.maximum_energy_solves = result.energy_solves - 1;
        }
        assert!(matches!(
            run(outer, caloric_controls),
            Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
                field: "unsteady.material_energy_solves",
                ..
            }))
        ));
    }
    assert_eq!(thermal, saved);
}
