//! Native-time source storage/advection and independent physical recurrences.
use super::*;
use fs_lbm::conjugate::{
    ChtError, FlowResistance, RadiationConfig, STEFAN_BOLTZMANN, TimeScheme, UnsteadyConfig,
};
use fs_lbm::material::conjugate::buoyant::enthalpy::MaterialEnthalpyTransientConfig;

fn transient_controls() -> MaterialEnthalpyTransientConfig {
    MaterialEnthalpyTransientConfig {
        caloric: enthalpy_controls(),
        maximum_passes_per_step: 80,
        maximum_energy_solves: 1000,
        temperature_tolerance_k: 1e-9,
        maximum_storage_defect_w: 1e-9,
        maximum_energy_closure_j: 1e-8,
        temperature_relaxation: 1.0,
    }
}

fn temperature_from_h(h: f64) -> f64 {
    // Independently integrate/invert Cp=2+.04*(T-300), h(300)=0.
    300.0 + 2.0 * h / (2.0 + (4.0 + 0.08 * h).sqrt())
}

#[test]
fn native_caloric_storage_and_open_advection_follow_exact_enthalpy_recurrences() {
    let card = card(0.04);
    let law = liquid(&card);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    for open in [false, true] {
        let mut datum_results = Vec::new();
        for datum in [0.0, 1e9] {
            let chart = caloric(&card, datum, 340.0);
            for reverse in [false, true] {
                let (domain, mut flow, mut thermal, energy) = channel(reverse);
                if !open {
                    flow = SimpleConfig::new([FvBoundary::wall(); 6]);
                    thermal.faces = [ThermalFace::Adiabatic; 6];
                }
                let unsteady = UnsteadyConfig::new(0.5, 5);
                let mut nonlinear = transient_controls();
                if datum != 0.0 {
                    // Also shift the shared numerical h datum, including both
                    // inlet and outlet, while leaving all physical heat fixed.
                    nonlinear.caloric.reference_temperature_k = 305.0;
                }
                let result = bound
                    .march_natural_convection_with_enthalpy(
                        &domain,
                        &[],
                        &flow,
                        &thermal,
                        &energy,
                        &[310.0],
                        controls().boussinesq,
                        &chart,
                        nonlinear,
                        &unsteady,
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
                let mdot = if open { 0.1 } else { 0.0 };
                for (step, record) in result.convection.march.energy_records.iter().enumerate() {
                    let old_h = h;
                    h = (h + unsteady.time_step_s * (2.6 + mdot * 22.0))
                        / (1.0 + unsteady.time_step_s * mdot);
                    assert!((record.max_temperature_k - temperature_from_h(h)).abs() < 1e-5);
                    near(record.stored_energy_change_j, h - old_h, 1e-8);
                    near(
                        record.boundary_outflow_j,
                        unsteady.time_step_s * mdot * (h - 22.0),
                        1e-8,
                    );
                    near(record.source_j, 1.3, 1e-12);
                    assert!(record.closure_j.abs() < 1e-8);
                    assert!(result.records[step].maximum_storage_defect_w <= 1e-9);
                    assert!(result.records[step].advection.maximum_cell_flux_defect_w <= 1e-10);
                }
                let t = result.convection.march.temperature[0];
                assert!((t - temperature_from_h(h)).abs() < 1e-5);
                near(result.convection.endpoint_time_s, 2.5, 1e-14);
                assert_eq!(result.convection.march.flow.records.len(), unsteady.steps);
                assert!(result.records.iter().any(|record| record.passes > 1));
                // Repeated nonlinear trials are not physical time steps.
                assert!(
                    result
                        .records
                        .iter()
                        .map(|record| record.energy_solves)
                        .sum::<usize>()
                        > unsteady.steps
                );
                assert!(core::ptr::eq(result.caloric, &chart));
                datum_results.push(t);
            }
        }
        for t in &datum_results[1..] {
            near(*t, datum_results[0], 1e-11);
        }
    }
}

#[test]
fn native_source_heat_storage_feeds_buoyancy_without_losing_momentum_history() {
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
        lo: [0, 0, 0],
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
        unsteady.probe = Some(1);
        let result = bound
            .march_natural_convection_with_enthalpy(
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
        // Four fixed reference masses each gain h(310)-h(300)=22 J/kg,
        // not Cp_reference*10=20. The thermostat supplies exactly this heat.
        near(
            result.convection.march.energy_records[0].stored_energy_change_j,
            88.0,
            1e-9,
        );
        near(
            result.convection.march.energy_records[0].boundary_outflow_j,
            -88.0,
            1e-9,
        );
        let (mut previous, mut older) = (0.0, 0.0);
        for (step, record) in result.convection.march.flow.records.iter().enumerate() {
            let expected = if step == 0 {
                0.0
            } else {
                match scheme {
                    TimeScheme::BackwardEuler => (previous + 0.25 * 0.0981) / 1.25,
                    TimeScheme::Bdf2 => (4.0 * previous - older + 0.5 * 0.0981) / 3.5,
                }
            };
            near(record.probe_velocity_m_s[0], expected, 2e-8);
            older = previous;
            previous = expected;
        }
        for record in &result.convection.march.energy_records {
            assert!(record.closure_j.abs() <= 1e-8, "{record:?}");
        }
        for sample in &result.convection.endpoint {
            near(sample.force_n_m3[0], 0.0981, 1e-12);
            near(
                sample
                    .material
                    .reference_storage_transport()
                    .density_kg_m3(),
                1.0,
                1e-14,
            );
            near(
                sample.material.resolved_transport().specific_heat_j_kg_k(),
                2.4,
                1e-12,
            );
        }
    }
}

#[test]
fn native_source_storage_retains_radiation_and_thermostat_heat_ownership() {
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
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let radiation = RadiationConfig::new(vec![0.9], surroundings);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let result = bound
        .march_natural_convection_with_enthalpy(
            &domain,
            &solids,
            &flow,
            &thermal,
            &energy,
            &[300.0; 2],
            controls().boussinesq,
            &chart,
            transient_controls(),
            &UnsteadyConfig::new(0.1, 1),
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
        .unwrap();
    let stored = 2.4e6 * domain.volume(0) * 30.0
        + domain.volume(1)
            * chart
                .specific_enthalpy_change(300.0, result.convection.march.temperature[1])
                .unwrap();
    let record = result.convection.march.energy_records[0];
    near(record.stored_energy_change_j, stored, 1e-12);
    near(record.boundary_outflow_j, -stored, 1e-10);
    near(
        record.radiated_w,
        0.9 * STEFAN_BOLTZMANN * dx * dx * (330.0_f64.powi(4) - 300.0_f64.powi(4)),
        1e-10,
    );
    assert!(record.closure_j.abs() <= 1e-8);
    // The global work ceiling counts the original radiation trial solves,
    // not just the outer caloric passes that requested them.
    assert!(result.records[0].energy_solves > result.records[0].passes);
}

#[test]
fn native_caloric_work_domain_boundary_and_cancellation_refusals_publish_no_partial_march() {
    let card = variable_card(0.04, 0.001, 300.0, 340.0);
    let law = liquid(&card);
    let chart = caloric(&card, 0.0, 340.0);
    let short_chart = caloric(&card, 0.0, 310.1);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let (domain, flow, thermal, energy) = channel(false);
    let saved = thermal.clone();
    let mut limits = controls().boussinesq;
    limits.max_relative_density_departure = 0.05;
    let unsteady = UnsteadyConfig::new(0.5, 2);
    let invoke = |caloric: &ResolvedHeatCapacityCurve, config, gate: &CancelGate| {
        bound
            .march_natural_convection_with_enthalpy(
                &domain,
                &[],
                &flow,
                &thermal,
                &energy,
                &[310.0],
                limits,
                caloric,
                config,
                &unsteady,
                |_| 1.0,
                |_| 1.0,
                None,
                gate,
                |_, t, _| {
                    law.at_temperature(Temperature::new(t))
                        .map_err(|error| MaterialError::State {
                            what: error.to_string(),
                        })
                },
            )
            .map(|result| result.convection.march.temperature)
    };
    let mut one_pass = transient_controls();
    one_pass.maximum_passes_per_step = 1;
    one_pass.temperature_tolerance_k = 100.0;
    assert!(matches!(invoke(&chart, one_pass, &CancelGate::new()),
        Err(MaterialBuoyancyError::EnthalpyStepNotConverged { storage_defect_w, .. })
            if storage_defect_w > 1e-3));
    let mut few_solves = transient_controls();
    few_solves.maximum_energy_solves = 2;
    assert!(matches!(
        invoke(&chart, few_solves, &CancelGate::new()),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "unsteady.material_energy_solves",
            ..
        }))
    ));
    assert!(matches!(
        invoke(&short_chart, transient_controls(), &CancelGate::new()),
        Err(MaterialBuoyancyError::Caloric { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        invoke(&chart, transient_controls(), &gate),
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));
    assert_eq!(thermal, saved);

    // Declared inlet/backflow support must pass the original density bound,
    // even though every cell and caloric tangent starts inside its support.
    let card = variable_card(0.0, 0.005, 280.0, 340.0);
    let law = liquid(&card);
    let chart = caloric_span(&card, 0.0, 280.0, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let mut thermal = saved;
    thermal.faces[0] = ThermalFace::Outflow {
        backflow_temperature: 280.0,
    };
    let result = bound.march_natural_convection_with_enthalpy(
        &domain,
        &[],
        &flow,
        &thermal,
        &energy,
        &[305.0],
        limits,
        &chart,
        transient_controls(),
        &unsteady,
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
    );
    assert!(matches!(
        result,
        Err(MaterialBuoyancyError::Material {
            temperature_k: 280.0,
            ..
        })
    ));
}
