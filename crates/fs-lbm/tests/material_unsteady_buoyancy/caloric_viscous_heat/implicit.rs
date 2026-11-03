//! Original molecular work closes with the same implicit caloric endpoint.
//! Synthetic laws exercise discrete balances, not material calibration.

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

// The original two closed MAC cells have no volume flow. Independently
// integrate h and apply the moving-wall half-cell stress quadrature:
// V * [h(T1)-h(T0)] = dt * (area/half_height)/2 * mu(T1) * Uwall^2.
fn wall_endpoint(old: f64, dt: f64, speed: f64, volume: f64) -> f64 {
    let residual = |t| volume * (h(t) - h(old)) - dt * 0.25 * mu(t) * speed * speed;
    let mut bracket = [
        old,
        old + dt * 0.25 * mu(old) * speed * speed / (2.0 * volume) + 1.0,
    ];
    assert!(residual(bracket[0]) < 0.0 && residual(bracket[1]) > 0.0);
    for _ in 0..80 {
        let t = 0.5 * (bracket[0] + bracket[1]);
        if residual(t) > 0.0 {
            bracket[1] = t;
        } else {
            bracket[0] = t;
        }
    }
    0.5 * (bracket[0] + bracket[1])
}

#[test]
fn g1_implicit_molecular_heat_matches_original_h_and_actual_endpoint_mu() {
    let domain = VoxelDomain::new(2, 1, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = wall_flow(0.1);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial = [300.0; 2];
    let unsteady = clock(5.0, 4);
    let caloric = Caloric;
    let gate = CancelGate::new();
    let run = |setup: &ThermalSetup| {
        bound.march_natural_convection_with_enthalpy_and_viscous_heat_implicit(
            &domain, &[], &flow, setup, &energy(), &initial, caloric_bounds(),
            &caloric, controls(), &unsteady, heat_controls(), endpoint_controls(),
            |_| 1.0, |_| 0.0, None, &gate, |_, t, _| source(t),
        )
    };
    let result = run(&thermal).unwrap();
    let heated = &result.solution;
    let march = &heated.thermal.convection.march;
    let mut expected = 300.0;
    for (index, (mechanical, heat)) in heated.mechanical_steps.iter()
        .zip(&march.energy_records).enumerate()
    {
        expected = wall_endpoint(expected, unsteady.time_step_s, 0.1, domain.volume(0));
        let expected_power = 0.25 * mu(expected) * 0.1_f64.powi(2);
        for &power in &mechanical.power.cell_heat_w {
            close(power, expected_power, 1e-9);
        }
        close(heat.source_j, 2.0 * expected_power * unsteady.time_step_s, 1e-8);
        close(heat.max_temperature_k, expected, 1e-7);
        close(heat.boundary_outflow_j, 0.0, 1e-9);
        close(mechanical.initial_kinetic_energy_j, 0.0, 1e-10);
        close(mechanical.final_kinetic_energy_j, 0.0, 1e-10);
        assert!(result.records[index].passes > 1);
    }
    for &temperature in &march.temperature {
        close(temperature, expected, 1e-7);
    }
    let stored: f64 = march.temperature.iter().enumerate()
        .map(|(cell, &t)| domain.volume(cell) * h(t)).sum();
    close(stored, heated.heating.generated_heat_j, 1e-8);
    close(heated.heating.supplied_wall_work_j, stored, 1e-8);
    assert_eq!(heated.heating.substeps, unsteady.steps);
    assert_eq!(heated.mechanical_steps.len(), unsteady.steps);
    assert_eq!(march.energy_records.len(), unsteady.steps);
    assert!(result.momentum_solves > unsteady.steps);
    close(heated.heating.duration_s, 20.0, 1e-14);

    // Repeated trials must replace physical heat even with an empty authored
    // field. Scheduling authored sources to zero cannot scale dissipation.
    let mut scheduled = thermal.clone();
    scheduled.power_w = vec![100.0; 2];
    let rescheduled = run(&scheduled).unwrap();
    assert_eq!(rescheduled.solution.thermal.convection.march.temperature, march.temperature);
    assert_eq!(rescheduled.solution.heating, heated.heating);
    let lagged = bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
        &caloric, controls(), &unsteady, heat_controls(), |_| 1.0, |_| 0.0,
        None, &gate, |_, t, _| source(t),
    ).unwrap();
    assert!(lagged.thermal.convection.march.temperature[0] > expected + 1e-3);

    // A thermostat changes the constitutive endpoint without accepting
    // another physical interval. Duplicate declarations reject heat once.
    let mut fixed = thermal.clone();
    fixed.fixed_temperature = vec![(0, 304.0), (0, 304.0), (1, 300.0)];
    let fixed_result = bound.march_natural_convection_with_enthalpy_and_viscous_heat_implicit(
        &domain, &[], &flow, &fixed, &energy(), &initial, caloric_bounds(),
        &caloric, controls(), &clock(0.5, 2), heat_controls(), endpoint_controls(),
        |_| 1.0, |_| 0.0, None, &gate, |_, t, _| source(t),
    ).unwrap();
    let fixed_heat = &fixed_result.solution;
    for (index, (mechanical, heat)) in fixed_heat.mechanical_steps.iter()
        .zip(&fixed_heat.thermal.convection.march.energy_records).enumerate()
    {
        for (&power, t) in mechanical.power.cell_heat_w.iter().zip([304.0, 300.0]) {
            close(power, 0.25 * mu(t) * 0.1_f64.powi(2), 1e-9);
        }
        let storage = if index == 0 { domain.volume(0) * h(304.0) } else { 0.0 };
        close(heat.stored_energy_change_j, storage, 1e-9);
        close(heat.source_j, mechanical.power.heat_w * 0.5, 1e-9);
        close(storage + heat.boundary_outflow_j, heat.source_j, 1e-8);
    }
    assert_eq!(fixed_heat.thermal.convection.march.temperature, [304.0, 300.0]);
    assert_eq!(initial, [300.0; 2]);
    assert_eq!(bound.transport(), &reference);
    assert!(thermal.power_w.is_empty());
}

#[test]
fn g4_actual_heat_and_reynolds_gates_survive_loose_material_tolerances() {
    let domain = VoxelDomain::new(2, 1, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = wall_flow(0.1);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial = [300.0; 2];
    let mut loose = endpoint_controls();
    loose.maximum_passes_per_step = 1;
    loose.temperature_tolerance_k = 1e6;
    loose.transport_relative_tolerance = 0.95;
    loose.force_tolerance_n_m3 = 1e6;
    let run = |law: u8| {
        bound.march_natural_convection_with_enthalpy_and_viscous_heat_implicit(
            &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
            &Caloric, controls(), &clock(5.0, 1), heat_controls(), loose,
            |_| 1.0, |_| 1.0, None, &CancelGate::new(), |_, t, _| {
                let viscosity = match law {
                    0 => 1.0,
                    1 => mu(t),
                    _ if t > 300.01 => 0.1,
                    _ => 1.0,
                };
                FluidTransportState::try_new(1.0, 2.0 + 0.05 * (t - 300.0), 0.01, viscosity)
            },
        )
    };
    // Positive control: constant mu has no endpoint mechanical/heat defect.
    assert_eq!(run(0).unwrap().records[0].passes, 1);
    // The moderate source update remains below Re=0.1 but its deposited
    // heat differs from actual endpoint stress work. Loose outer gates cannot
    // accept the deliberately incomplete single nonlinear pass.
    assert!(matches!(run(1),
        Err(MaterialBuoyancyError::Solver(ChtError::SolverNotConverged { .. }))));
    // The actual endpoint mu=0.1 gives Re=0.5. The outer relative defect
    // 0.9 fits its declared 0.95 tolerance, but native applicability must refuse.
    assert!(matches!(run(2),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "simple.viscous_heat", ..
        }))));
    assert_eq!(initial, [300.0; 2]);
    assert!(thermal.power_w.is_empty());
}

#[test]
fn g1_g4_implicit_viscous_trials_retain_native_history_and_charge_actual_work() {
    let domain = VoxelDomain::new(4, 4, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let initial = vec![300.0; domain.cell_count()];
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let flow = wall_flow(0.025);
    let unsteady = clock(0.1, 3);
    let caloric = Caloric;
    let gate = CancelGate::new();
    let run = |outer, original| {
        bound.march_natural_convection_with_enthalpy_and_viscous_heat_implicit(
            &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
            &caloric, original, &unsteady, heat_controls(), outer,
            |_| 1.0, |_| 1.0, None, &gate, |_, t, _| source(t),
        )
    };
    let result = run(endpoint_controls(), controls()).unwrap();
    let heated = &result.solution;
    assert_eq!(heated.mechanical_steps.len(), unsteady.steps);
    assert_eq!(heated.heating.substeps, unsteady.steps);
    assert_eq!(result.records.len(), unsteady.steps);
    assert!(result.momentum_solves > unsteady.steps);
    assert!(heated.heating.final_kinetic_energy_j > 1e-9);
    assert!(heated.heating.backward_euler_loss_j > 1e-9);
    for pair in heated.mechanical_steps.windows(2) {
        assert_eq!(pair[0].final_kinetic_energy_j.to_bits(), pair[1].initial_kinetic_energy_j.to_bits());
        assert_eq!(pair[0].end_s.to_bits(), pair[1].start_s.to_bits());
    }
    let march = &heated.thermal.convection.march;
    let stored: f64 = march.temperature.iter().enumerate()
        .map(|(cell, &t)| domain.volume(cell) * h(t)).sum();
    close(stored, heated.heating.generated_heat_j, 1e-8);
    close(march.energy_records.iter().map(|r| r.source_j).sum(), heated.heating.generated_heat_j, 1e-8);
    let kinetic = heated.heating.final_kinetic_energy_j - heated.heating.initial_kinetic_energy_j;
    let supplied = heated.heating.supplied_wall_work_j + heated.heating.supplied_body_work_j
        + heated.heating.pressure_constraint_work_j;
    close(supplied, stored + kinetic + heated.heating.backward_euler_loss_j, 1e-8);
    assert!((stored - heated.heating.generated_heat_j).abs()
        < 0.05 * heated.heating.backward_euler_loss_j,
        "numerical BE loss stays outside actual source enthalpy");
    assert_eq!(result.energy_solves,
        heated.thermal.records.iter().map(|r| r.energy_solves).sum::<usize>());
    let spent = result.energy_solves;
    assert!(spent > unsteady.steps);
    let mut exact = endpoint_controls();
    exact.maximum_energy_solves = spent;
    assert_eq!(run(exact, controls()).unwrap().energy_solves, spent);
    let mut outer_short = exact;
    outer_short.maximum_energy_solves -= 1;
    assert!(matches!(run(outer_short, controls()),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "unsteady.material_energy_solves", ..
        }))));
    let mut inner_short = controls();
    inner_short.maximum_energy_solves = spent - 1;
    assert!(matches!(run(endpoint_controls(), inner_short),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "unsteady.material_energy_solves", ..
        }))));
    assert_eq!(initial, vec![300.0; domain.cell_count()]);
    assert_eq!(bound.transport(), &reference);
    assert!(thermal.power_w.is_empty());
}

#[test]
fn g1_implicit_viscous_heat_and_source_solids_share_one_energy_account() {
    let fluid_card = card(0.004, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let transport = |t| {
        let state = super::super::transport(&fluid_source, t)?;
        FluidTransportState::try_new(
            state.density_kg_m3(), state.specific_heat_j_kg_k(),
            state.thermal_conductivity_w_m_k(), mu(t),
        )
    };
    let reference = transport(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let (domain, _, _, _, _) = shear_cavity(true);
    let initial = vec![300.0; domain.cell_count()];
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let unsteady = clock(0.25, 3);
    let solid_controls = SolidThermalConfig {
        temperature_tolerance_k: 1e-10,
        conductivity_relative_tolerance: 1e-11,
        storage_tolerance_j: 1e-9,
        ..SolidThermalConfig::default()
    };
    let result = bound.march_natural_convection_with_enthalpy_and_material_solids_and_viscous_heat_implicit(
        &domain, std::slice::from_ref(&solid), &wall_flow(0.01), &thermal, &energy(),
        &initial, bounds(), &caloric, controls(), solid_controls, &unsteady,
        heat_controls(), endpoint_controls(), |_| 1.0, |_| 1.0, None,
        &CancelGate::new(), |_, t, _| transport(t),
    ).unwrap();
    let heated = &result.solution;
    let joint = &heated.thermal;
    let march = &joint.fluid.convection.march;
    let stored: f64 = march.temperature.iter().enumerate().map(|(cell, &t)| {
        let slope = if domain.is_fluid(cell) { 0.004 } else { 0.2 };
        3.0 * domain.volume(cell) * (specific_energy(t, slope) - specific_energy(300.0, slope))
    }).sum();
    close(stored, heated.heating.generated_heat_j, 1e-8);
    assert_eq!(joint.solid_records.len(), unsteady.steps);
    assert_eq!(heated.mechanical_steps.len(), unsteady.steps);
    assert!(joint.solid_records.iter().map(|r| r.solid_stored_energy_change_j).sum::<f64>() > 0.0);
    assert_eq!(result.energy_solves,
        joint.fluid.records.iter().map(|r| r.energy_solves).sum::<usize>());
    for step in &heated.mechanical_steps {
        for (cell, &heat) in step.power.cell_heat_w.iter().enumerate() {
            if !domain.is_fluid(cell) {
                assert_eq!(heat, 0.0);
            }
        }
    }
    for heat in &march.energy_records {
        close(heat.stored_energy_change_j + heat.boundary_outflow_j, heat.source_j, 1e-8);
    }
}
