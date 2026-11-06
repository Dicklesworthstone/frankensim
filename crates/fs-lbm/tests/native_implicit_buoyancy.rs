//! Native endpoint buoyancy, physical heat closure and retained time histories.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    Boussinesq, ChtError, ConvectionScheme, EnergyConfig, FlowResistance, FluidProperties,
    FvBoundary, MaterialEndpointConfig, SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup,
    TimeScheme, UnsteadyConfig, Voxel, VoxelDomain, march_conjugate,
    march_conjugate_implicit_buoyancy,
};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}
fn controls() -> MaterialEndpointConfig {
    MaterialEndpointConfig {
        maximum_passes_per_step: 80,
        maximum_momentum_solves: 1000,
        maximum_energy_solves: 1000,
        maximum_source_refreshes: 2001,
        temperature_tolerance_k: 1e-9,
        transport_relative_tolerance: 1e-10,
        force_tolerance_n_m3: 1e-10,
        temperature_relaxation: 1.0,
    }
}
fn energy() -> EnergyConfig {
    EnergyConfig {
        scheme: ConvectionScheme::EnthalpyMinmod,
        tolerance: 1e-11,
        max_iterations: 2000,
        ..EnergyConfig::default()
    }
}
fn clock(dt: f64, steps: usize) -> UnsteadyConfig {
    let mut clock = UnsteadyConfig::new(dt, steps);
    clock.inner_iterations = 1000;
    clock.inner_tolerance = 1e-10;
    clock.average_from_step = 1;
    clock
}
fn buoyancy(expansion_per_k: f64, gravity_m_s2: [f64; 3]) -> Boussinesq {
    Boussinesq {
        gravity_m_s2,
        expansion_per_k,
        reference_temperature_k: 300.0,
    }
}
fn cavity() -> (
    VoxelDomain,
    SimpleConfig,
    ThermalSetup,
    Vec<SolidMaterial>,
    Vec<f64>,
) {
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
    let initial = (0..domain.cell_count())
        .map(|cell| if domain.is_fluid(cell) { 300.0 } else { 320.0 })
        .collect();
    (
        domain,
        SimpleConfig::new(faces),
        ThermalSetup::new([ThermalFace::Adiabatic; 6]),
        vec![SolidMaterial::new("conducting-wall", 1.0).with_heat_capacity(2.0)],
        initial,
    )
}

#[test]
fn g1_hot_wall_drives_the_first_interval_and_keeps_one_physical_heat_ledger() {
    let (domain, flow, thermal, solids, initial) = cavity();
    let b = buoyancy(0.001, [0.0, -9.81, 0.0]);
    let clock = clock(0.1, 1);
    let result = march_conjugate_implicit_buoyancy(
        &domain,
        &fluid(),
        &solids,
        &thermal,
        &flow,
        &clock,
        &b,
        controls(),
        &initial,
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    let march = &result.solution;
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
        "newly transferred wall heat must drive circulation"
    );
    assert!(result.records[0].passes > 1);
    assert!(result.records[0].temperature_defect_k <= controls().temperature_tolerance_k);
    assert!(result.records[0].force_defect_n_m3 <= controls().force_tolerance_n_m3);
    assert_eq!(march.flow.records.len(), 1);
    assert_eq!(march.energy_records.len(), 1);
    assert_eq!(march.flow.averaged_steps, 1);
    let stored = initial
        .iter()
        .zip(&march.temperature)
        .enumerate()
        .map(|(cell, (&old, &new))| 2.0 * domain.volume(cell) * (new - old))
        .sum::<f64>();
    assert!(stored.abs() < 1e-8, "closed cavity gained {stored} J");
    assert!(march.energy_records[0].closure_j.abs() < 1e-8);
    assert_eq!(march.energy_records[0].source_j, 0.0);
    assert_eq!(march.energy_records[0].boundary_outflow_j, 0.0);
    let lagged = march_conjugate(
        &domain,
        &fluid(),
        &solids,
        &thermal,
        &flow,
        &clock,
        Some(&b),
        &initial,
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(
        lagged
            .flow
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .all(|u| u.abs() < 1e-14)
    );

    let uniform = vec![300.0; domain.cell_count()];
    let equilibrium = march_conjugate_implicit_buoyancy(
        &domain,
        &fluid(),
        &solids,
        &thermal,
        &flow,
        &clock,
        &b,
        controls(),
        &uniform,
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(
        equilibrium
            .solution
            .temperature
            .iter()
            .all(|t| (t - 300.0).abs() < 1e-9)
    );
    assert!(
        equilibrium
            .solution
            .flow
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .all(|u| u.abs() < 1e-12)
    );
}

#[test]
fn g1_be_and_bdf2_use_each_physical_history_once_for_both_expansion_signs() {
    // Uniform column: du/dt = F - mu*u/K. This recurrence is independent
    // of the native pressure, temperature and coupling iteration machinery.
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
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        for beta in [-0.001, 0.001] {
            let mut clock = clock(0.25, 6);
            clock.scheme = scheme;
            clock.probe = Some(1);
            let result = march_conjugate_implicit_buoyancy(
                &domain,
                &fluid(),
                &[],
                &thermal,
                &flow,
                &clock,
                &buoyancy(beta, [-9.81, 0.0, 0.0]),
                controls(),
                &[300.0; 4],
                |_| 1.0,
                |_| 1.0,
                &energy(),
                None,
                &CancelGate::new(),
            )
            .unwrap();
            let force = 9.81 * beta * 10.0;
            let (mut previous, mut older, mut sum) = (0.0, 0.0, 0.0);
            for (step, record) in result.solution.flow.records.iter().enumerate() {
                let dt = clock.time_step_s;
                let expected = if step == 0 || scheme == TimeScheme::BackwardEuler {
                    (previous + dt * force) / (1.0 + dt)
                } else {
                    (4.0 * previous - older + 2.0 * dt * force) / (3.0 + 2.0 * dt)
                };
                assert!((record.probe_velocity_m_s[0] - expected).abs() < 2e-8);
                older = previous;
                previous = expected;
                sum += expected;
            }
            let march = &result.solution;
            assert_eq!(result.records[0].passes, 2);
            assert_eq!(result.momentum_solves, clock.steps + 1);
            assert_eq!(march.flow.records.len(), clock.steps);
            assert_eq!(march.energy_records.len(), clock.steps);
            assert_eq!(march.flow.averaged_steps, clock.steps);
            assert!((march.flow.mean_velocity_m_s[1][0] - sum / clock.steps as f64).abs() < 2e-8);
            assert!((march.flow.flow.field.flux_x(2, 0, 0) - previous).abs() < 2e-8);
            assert!((march.energy_records[0].stored_energy_change_j - 80.0).abs() < 1e-8);
            assert!(
                march
                    .energy_records
                    .iter()
                    .all(|record| record.closure_j.abs() < 1e-8)
            );
        }
    }
}

#[test]
fn g1_zero_buoyancy_preserves_integrated_solid_conductivity_and_native_temperature() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| {
        if p[0] < 2.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let solid = SolidMaterial::new("linear-k", 1.0)
        .with_heat_capacity(1.0)
        .with_conductivity_table(&[(300.0, 1.0), (320.0, 5.0)]);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let mut thermal = ThermalSetup::new(faces);
    thermal.fixed_temperature = vec![(1, 300.0), (2, 300.0)];
    thermal.power_w = vec![15.0, 0.0, 0.0];
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let clock = clock(1.0, 1);
    let b = buoyancy(0.0, [0.0, -9.81, 0.0]);
    let result = march_conjugate_implicit_buoyancy(
        &domain,
        &fluid(),
        std::slice::from_ref(&solid),
        &thermal,
        &flow,
        &clock,
        &b,
        controls(),
        &[300.0; 3],
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    // k=1+0.2*(T-300). Native half-cell and same-material paths have
    // reciprocal lengths 2 and 1: d + 3*(d + 0.1*d*d) = 15 J.
    let expected = 300.0 + 30.0 / (4.0 + 34.0_f64.sqrt());
    assert!((result.solution.temperature[0] - expected).abs() < 2e-8);
    assert!(result.solution.energy_records[0].closure_j.abs() < 1e-8);
    assert!(
        result
            .solution
            .flow
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .all(|u| u.abs() < 1e-14)
    );
    let native = march_conjugate(
        &domain,
        &fluid(),
        &[solid],
        &thermal,
        &flow,
        &clock,
        Some(&b),
        &[300.0; 3],
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for (implicit, old) in result.solution.temperature.iter().zip(&native.temperature) {
        assert!((implicit - old).abs() < 1e-10);
    }
}

#[test]
fn g4_cancelled_or_unclosed_coupling_publishes_no_native_step() {
    let (domain, flow, thermal, solids, initial) = cavity();
    let clock = clock(0.1, 1);
    let b = buoyancy(0.001, [0.0, -9.81, 0.0]);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        march_conjugate_implicit_buoyancy(
            &domain,
            &fluid(),
            &solids,
            &thermal,
            &flow,
            &clock,
            &b,
            controls(),
            &initial,
            |_| 1.0,
            |_| 1.0,
            &energy(),
            None,
            &gate,
        ),
        Err(ChtError::Cancelled)
    ));
    let gate = CancelGate::new();
    assert!(matches!(
        march_conjugate_implicit_buoyancy(
            &domain,
            &fluid(),
            &solids,
            &thermal,
            &flow,
            &clock,
            &b,
            controls(),
            &initial,
            |_| 1.0,
            |_| {
                gate.request();
                1.0
            },
            &energy(),
            None,
            &gate,
        ),
        Err(ChtError::Cancelled)
    ));
    let mut insufficient = controls();
    insufficient.maximum_passes_per_step = 1;
    assert!(matches!(
        march_conjugate_implicit_buoyancy(
            &domain,
            &fluid(),
            &solids,
            &thermal,
            &flow,
            &clock,
            &b,
            insufficient,
            &initial,
            |_| 1.0,
            |_| 1.0,
            &energy(),
            None,
            &CancelGate::new(),
        ),
        Err(ChtError::SolverNotConverged {
            system: "implicit transient material coupling",
            ..
        })
    ));
}
