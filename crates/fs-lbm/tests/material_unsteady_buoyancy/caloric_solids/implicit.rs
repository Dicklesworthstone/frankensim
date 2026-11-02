//! Original solid calorics coupled to endpoint fluid transport and buoyancy.
use super::*;
use fs_lbm::material::conjugate::buoyant::MaterialEndpointConfig;

fn endpoint_controls() -> MaterialEndpointConfig {
    MaterialEndpointConfig {
        maximum_passes_per_step: 80,
        maximum_momentum_solves: 1000,
        maximum_energy_solves: 2000,
        maximum_source_refreshes: 2001,
        temperature_tolerance_k: 1e-10,
        transport_relative_tolerance: 1e-11,
        force_tolerance_n_m3: 1e-11,
        temperature_relaxation: 1.0,
    }
}

fn solid_controls() -> SolidThermalConfig {
    SolidThermalConfig {
        temperature_tolerance_k: 1e-10,
        conductivity_relative_tolerance: 1e-11,
        storage_tolerance_j: 1e-9,
        ..SolidThermalConfig::default()
    }
}

fn pair() -> VoxelDomain {
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain
}

// Independently eliminate the two original masses' energies, then bracket
// one scalar backward-Euler face balance. No native coefficient, source
// integral, caloric linearization or nonlinear iterator supplies this oracle.
fn pair_step(old: [f64; 2], fluid_slope: f64, dt: f64) -> [f64; 2] {
    let solid_h = |t| specific_energy(t, 0.2);
    let fluid_h = |t| specific_energy(t, fluid_slope);
    let total = solid_h(old[0]) + fluid_h(old[1]);
    let mut equilibrium = old;
    for _ in 0..80 {
        let t = 0.5 * (equilibrium[0] + equilibrium[1]);
        if solid_h(t) + fluid_h(t) > total {
            equilibrium[1] = t;
        } else {
            equilibrium[0] = t;
        }
    }
    let fluid_at = |solid_t: f64| {
        let target = total - solid_h(solid_t);
        let mut bracket = [solid_t, old[1]];
        for _ in 0..80 {
            let t = 0.5 * (bracket[0] + bracket[1]);
            if fluid_h(t) > target {
                bracket[1] = t;
            } else {
                bracket[0] = t;
            }
        }
        0.5 * (bracket[0] + bracket[1])
    };
    let mut bracket = [old[0], 0.5 * (equilibrium[0] + equilibrium[1])];
    for _ in 0..80 {
        let ts = 0.5 * (bracket[0] + bracket[1]);
        let tf = fluid_at(ts);
        let k_s = 1.0 + 0.04 * (ts - 300.0);
        let k_f = 1.0 + 0.02 * (tf - 300.0);
        let conductance = 1.0 / (0.5 / k_s + 0.5 / k_f);
        if 3.0 * (solid_h(ts) - solid_h(old[0])) > dt * conductance * (tf - ts) {
            bracket[1] = ts;
        } else {
            bracket[0] = ts;
        }
    }
    let ts = 0.5 * (bracket[0] + bracket[1]);
    [ts, fluid_at(ts)]
}

#[test]
fn g1_joint_implicit_solid_and_fluid_sources_match_independent_heat_balance() {
    let domain = pair();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial = [300.0, 320.0];
    let unsteady = UnsteadyConfig::new(1.0, 3);
    let solid_card = card(0.2, 0.04);
    let solid = law(&solid_card, [290.0, 400.0]);
    for fluid_slope in [0.0, 0.002] {
        let fluid_card = card(fluid_slope, 0.02);
        let source = law(&fluid_card, [290.0, 400.0]);
        let caloric = source_caloric(&fluid_card, &source, [290.0, 400.0]);
        let reference = transport(&source, 300.0).unwrap();
        let bound = MaterialChtFluid::try_new(&reference).unwrap();
        let result = bound
            .march_natural_convection_with_enthalpy_and_material_solids_implicit(
                &domain,
                std::slice::from_ref(&solid),
                &flow,
                &thermal,
                &energy(),
                &initial,
                bounds(),
                &caloric,
                controls(),
                solid_controls(),
                &unsteady,
                endpoint_controls(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, t, _| transport(&source, t),
            )
            .unwrap();
        let joint = &result.solution;
        let march = &joint.fluid.convection.march;
        let mut expected = initial;
        for (step, record) in march.energy_records.iter().enumerate() {
            let old = expected;
            expected = pair_step(old, fluid_slope, unsteady.time_step_s);
            let stored_s = 3.0 * (specific_energy(expected[0], 0.2) - specific_energy(old[0], 0.2));
            close(
                joint.solid_records[step].solid_stored_energy_change_j,
                stored_s,
                2e-8,
            );
            close(record.stored_energy_change_j, 0.0, 2e-8);
            assert!(record.closure_j.abs() <= 1e-8);
            assert_eq!(
                joint.fluid.records[step].time_s,
                joint.solid_records[step].time_s
            );
            assert_eq!(
                joint.fluid.records[step].passes,
                joint.solid_records[step].iterations
            );
            assert!(joint.fluid.records[step].passes > result.records[step].passes);
        }
        for (&actual, expected) in march.temperature.iter().zip(expected) {
            assert!((actual - expected).abs() < 2e-8);
        }
        assert_eq!(march.flow.records.len(), unsteady.steps);
        assert_eq!(joint.fluid.records.len(), unsteady.steps);
        assert_eq!(joint.solid_records.len(), unsteady.steps);
        assert_eq!(result.records.len(), unsteady.steps);
        assert!(result.momentum_solves > unsteady.steps);
        assert_eq!(
            result.energy_solves,
            joint
                .fluid
                .records
                .iter()
                .map(|r| r.energy_solves)
                .sum::<usize>()
        );
        assert!(result.energy_solves > result.momentum_solves);
        assert!(core::ptr::eq(joint.fluid.caloric, &caloric));
        assert_eq!(
            joint.solid_endpoint[0].material.card_identity(),
            solid_card.content_hash()
        );
        close(
            joint.solid_endpoint[0].temperature_k,
            march.temperature[0],
            1e-13,
        );
        close(
            joint.fluid.convection.endpoint[0].temperature_k,
            march.temperature[1],
            1e-13,
        );
        assert!(
            (joint.fluid.convection.applied[0].temperature_k - march.temperature[1]).abs()
                <= endpoint_controls().temperature_tolerance_k
        );

        // The constant-Cp specialization uses the same original solid source
        // and physical endpoint, while its independent plain heat path also
        // counts every inner solid solve against the endpoint work ceiling.
        if fluid_slope == 0.0 {
            let plain = bound
                .march_natural_convection_with_material_solids_implicit(
                    &domain,
                    std::slice::from_ref(&solid),
                    &flow,
                    &thermal,
                    &energy(),
                    &initial,
                    bounds(),
                    &unsteady,
                    solid_controls(),
                    endpoint_controls(),
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| transport(&source, t),
                )
                .unwrap();
            for (&actual, &expected) in plain
                .solution
                .fluid
                .march
                .temperature
                .iter()
                .zip(&march.temperature)
            {
                assert!((actual - expected).abs() < 2e-8);
            }
            assert_eq!(plain.solution.solid_records.len(), unsteady.steps);
            assert_eq!(
                plain.solution.fluid.march.flow.records.len(),
                unsteady.steps
            );
            assert!(plain.energy_solves > plain.momentum_solves);
            assert_eq!(
                plain.energy_solves,
                plain
                    .solution
                    .solid_records
                    .iter()
                    .map(|r| r.iterations)
                    .sum::<usize>()
            );
        }
    }
}

#[test]
fn g1_g3_heat_from_the_original_solid_drives_buoyancy_within_the_first_step() {
    let fluid_card = card(0.0, 0.02);
    let source = law(&fluid_card, [290.0, 400.0]);
    let reference = transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let domain = VoxelDomain::from_fn(4, 4, 1, 0.25, |p| {
        if p[0] < 0.25 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let flow = SimpleConfig::new(faces);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial: Vec<_> = (0..domain.cell_count())
        .map(|cell| if domain.is_fluid(cell) { 300.0 } else { 330.0 })
        .collect();
    let mut unsteady = UnsteadyConfig::new(0.1, 1);
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-10;
    let limits = limits([0.0, -9.81, 0.0]);
    let mut outcomes = Vec::new();
    for k_slope in [0.0, 0.1] {
        let solid_card = card(0.2, k_slope);
        let solid = law(&solid_card, [290.0, 400.0]);
        let result = bound
            .march_natural_convection_with_material_solids_implicit(
                &domain,
                std::slice::from_ref(&solid),
                &flow,
                &thermal,
                &energy(),
                &initial,
                limits,
                &unsteady,
                solid_controls(),
                endpoint_controls(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, t, _| transport(&source, t),
            )
            .unwrap();
        let march = &result.solution.fluid.march;
        let speed = march
            .flow
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        assert!(
            speed > 1e-8,
            "actual evolved wall heat must drive the first interval"
        );
        let mut fluid_heat = 0.0;
        let mut total_change = 0.0;
        for (cell, (&old, &new)) in initial.iter().zip(&march.temperature).enumerate() {
            let change = if domain.is_fluid(cell) {
                let heat = 6.0 * domain.volume(cell) * (new - old);
                fluid_heat += heat;
                heat
            } else {
                3.0 * domain.volume(cell) * (specific_energy(new, 0.2) - specific_energy(old, 0.2))
            };
            total_change += change;
        }
        assert!(total_change.abs() < 1e-8);
        close(
            result.solution.solid_records[0].solid_stored_energy_change_j,
            -fluid_heat,
            1e-8,
        );
        assert_eq!(result.solution.solid_records.len(), 1);
        assert_eq!(march.flow.records.len(), 1);
        assert_eq!(march.energy_records.len(), 1);
        assert_eq!(march.flow.averaged_steps, 1);
        assert!(result.records[0].passes > 1);
        close(
            result.solution.fluid.applied_time_s,
            unsteady.time_step_s,
            1e-14,
        );
        for sample in &result.solution.solid_endpoint {
            close(sample.temperature_k, march.temperature[sample.cell], 1e-13);
        }
        let explicit = bound
            .march_natural_convection_with_material_solids(
                &domain,
                std::slice::from_ref(&solid),
                &flow,
                &thermal,
                &energy(),
                &initial,
                limits,
                &unsteady,
                solid_controls(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, t, _| transport(&source, t),
            )
            .unwrap();
        assert!(
            explicit
                .fluid
                .march
                .flow
                .flow
                .velocity_m_s
                .iter()
                .flatten()
                .all(|u| u.abs() < 1e-14)
        );
        outcomes.push((fluid_heat, speed));
    }
    assert!(outcomes[1].0 > outcomes[0].0 + 1e-6);
    assert!((outcomes[1].1 - outcomes[0].1).abs() > 1e-9);
    assert_eq!(bound.transport(), &reference);
}

#[test]
fn g4_source_solid_trials_share_whole_march_heat_limits_and_publish_no_failed_endpoint() {
    let fluid_card = card(0.002, 0.02);
    let source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &source, [290.0, 400.0]);
    let reference = transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.04);
    let solid = law(&solid_card, [290.0, 400.0]);
    let domain = pair();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial = [300.0, 320.0];
    let unsteady = UnsteadyConfig::new(1.0, 2);
    let run = |joint, endpoint, caloric_controls| {
        if joint {
            bound
                .march_natural_convection_with_enthalpy_and_material_solids_implicit(
                    &domain,
                    std::slice::from_ref(&solid),
                    &flow,
                    &thermal,
                    &energy(),
                    &initial,
                    bounds(),
                    &caloric,
                    caloric_controls,
                    solid_controls(),
                    &unsteady,
                    endpoint,
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| transport(&source, t),
                )
                .map(|result| result.energy_solves)
        } else {
            bound
                .march_natural_convection_with_material_solids_implicit(
                    &domain,
                    std::slice::from_ref(&solid),
                    &flow,
                    &thermal,
                    &energy(),
                    &initial,
                    bounds(),
                    &unsteady,
                    solid_controls(),
                    endpoint,
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| transport(&source, t),
                )
                .map(|result| result.energy_solves)
        }
    };
    for joint in [false, true] {
        let spent = run(joint, endpoint_controls(), controls()).unwrap();
        assert!(spent > unsteady.steps);
        let mut outer = endpoint_controls();
        outer.maximum_energy_solves = spent - 1;
        assert!(matches!(
            run(joint, outer, controls()),
            Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
                field: "unsteady.material_energy_solves",
                ..
            }))
        ));
        if joint {
            let mut original = controls();
            original.maximum_energy_solves = spent - 1;
            assert!(matches!(
                run(joint, endpoint_controls(), original),
                Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
                    field: "unsteady.material_energy_solves",
                    ..
                }))
            ));
        }
    }
    let mut refused_queries = 0;
    let refused = bound.march_natural_convection_with_enthalpy_and_material_solids_implicit(
        &domain,
        std::slice::from_ref(&solid),
        &flow,
        &thermal,
        &energy(),
        &initial,
        bounds(),
        &caloric,
        controls(),
        solid_controls(),
        &unsteady,
        endpoint_controls(),
        |_| 1.0,
        |_| 1.0,
        None,
        &CancelGate::new(),
        |_, t, _| {
            if t < 319.99 {
                refused_queries += 1;
                Err(MaterialError::State {
                    what: "original fluid source support refused".into(),
                })
            } else {
                transport(&source, t)
            }
        },
    );
    assert!(matches!(
        refused,
        Err(MaterialBuoyancyError::Material { .. })
    ));
    assert_eq!(refused_queries, 1);
    assert_eq!(initial, [300.0, 320.0]);
    assert_eq!(solid.conductivity().reference().query().axes()["T"], 300.0);
    assert_eq!(solid.storage().temperature_span_k(), [290.0, 400.0]);
    assert_eq!(bound.transport(), &reference);
}

#[test]
fn g1_g3_implicit_sgs_and_solid_calorics_share_actual_flow_and_one_heat_ledger() {
    let fluid_card = card(0.0, 0.0);
    let source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &source, [290.0, 400.0]);
    let reference = shear_transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let (domain, flow, thermal, initial, mut unsteady) = shear_cavity(true);
    unsteady.inner_tolerance = 1e-11;
    let run = |cs| {
        bound
            .march_natural_convection_with_enthalpy_and_material_solids_and_sgs_implicit(
                &domain,
                std::slice::from_ref(&solid),
                &flow,
                &thermal,
                &energy(),
                &initial,
                bounds(),
                &caloric,
                controls(),
                solid_controls(),
                &unsteady,
                sgs(cs),
                endpoint_controls(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, t, _| shear_transport(&source, t),
            )
            .unwrap()
    };
    let active = run(0.5);
    let zero = run(0.0);
    let laminar = bound
        .march_natural_convection_with_enthalpy_and_material_solids_implicit(
            &domain,
            std::slice::from_ref(&solid),
            &flow,
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &caloric,
            controls(),
            solid_controls(),
            &unsteady,
            endpoint_controls(),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| shear_transport(&source, t),
        )
        .unwrap();
    let plain = bound
        .march_natural_convection_with_material_solids_and_sgs_implicit(
            &domain,
            std::slice::from_ref(&solid),
            &flow,
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &unsteady,
            solid_controls(),
            sgs(0.5),
            endpoint_controls(),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| shear_transport(&source, t),
        )
        .unwrap();
    let march = &active.solution.fluid.convection.march;
    let zero_march = &zero.solution.fluid.convection.march;
    let laminar_march = &laminar.solution.fluid.convection.march;
    assert!(
        march
            .flow
            .flow
            .eddy_viscosity_m2_s
            .iter()
            .any(|&nu| nu > 1e-7)
    );
    for ((&active_t, &plain_t), (&zero_t, &laminar_t)) in march
        .temperature
        .iter()
        .zip(&plain.solution.fluid.march.temperature)
        .zip(
            zero_march
                .temperature
                .iter()
                .zip(&laminar_march.temperature),
        )
    {
        assert!((active_t - plain_t).abs() < 1e-8);
        assert!((zero_t - laminar_t).abs() < 1e-10);
    }
    for (actual, expected) in zero_march
        .flow
        .flow
        .velocity_m_s
        .iter()
        .zip(&laminar_march.flow.flow.velocity_m_s)
    {
        for (&a, &b) in actual.iter().zip(expected) {
            close(a, b, 1e-10);
        }
    }
    assert!(
        march
            .temperature
            .iter()
            .zip(&zero_march.temperature)
            .any(|(a, b)| (a - b).abs() > 1e-6)
    );
    let mut total_change = 0.0;
    let mut solid_change = 0.0;
    for (cell, (&new, &old)) in march.temperature.iter().zip(&initial).enumerate() {
        let slope = if domain.is_fluid(cell) { 0.0 } else { 0.2 };
        let delta =
            3.0 * domain.volume(cell) * (specific_energy(new, slope) - specific_energy(old, slope));
        total_change += delta;
        if !domain.is_fluid(cell) {
            solid_change += delta;
        }
    }
    let supplied =
        unsteady.time_step_s * unsteady.steps as f64 * thermal.power_w.iter().sum::<f64>();
    close(total_change, supplied, 1e-8);
    close(
        active
            .solution
            .solid_records
            .iter()
            .map(|r| r.solid_stored_energy_change_j)
            .sum(),
        solid_change,
        1e-8,
    );
    assert_eq!(active.solution.solid_records.len(), unsteady.steps);
    assert_eq!(active.solution.fluid.records.len(), unsteady.steps);
    assert_eq!(march.flow.records.len(), unsteady.steps);
    assert_eq!(
        active.energy_solves,
        active
            .solution
            .fluid
            .records
            .iter()
            .map(|r| r.energy_solves)
            .sum::<usize>()
    );
    for (fluid, solid) in active
        .solution
        .fluid
        .records
        .iter()
        .zip(&active.solution.solid_records)
    {
        assert_eq!(fluid.passes, solid.iterations);
        assert!(fluid.maximum_sgs_heat_flux_defect_w <= controls().caloric.maximum_flux_defect_w);
    }
    for sample in &active.solution.fluid.convection.endpoint {
        close(sample.temperature_k, march.temperature[sample.cell], 1e-13);
        close(
            sample
                .material
                .resolved_transport()
                .dynamic_viscosity_pa_s(),
            0.03 / (1.0 + 0.15 * (sample.temperature_k - 300.0)),
            1e-13,
        );
    }
    for sample in &active.solution.solid_endpoint {
        close(sample.temperature_k, march.temperature[sample.cell], 1e-13);
        for k in sample.conductivity_w_m_k {
            close(k, 1.0 + 0.01 * (sample.temperature_k - 300.0), 1e-12);
        }
    }
}
