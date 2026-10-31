//! Synthetic constitutive fixtures exercise the native transient equations.
//! They are numerical regressions, not empirical material qualifications.

#[path = "material_unsteady_buoyancy/solid.rs"]
mod solid;

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, EnergyConfig, FlowResistance, FvBoundary, RadiationConfig, STEFAN_BOLTZMANN,
    SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup, TimeScheme, UnsteadyConfig, Voxel,
    VoxelDomain, march_conjugate,
};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::conjugate::buoyant::MaterialBuoyancyError;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::{MaterialError, fluid_transport::FluidTransportState};

fn limits(gravity_m_s2: [f64; 3]) -> BoussinesqConfig<3> {
    BoussinesqConfig {
        gravity_m_s2,
        max_relative_density_departure: 0.05,
        max_relative_specific_heat_departure: 0.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_temperature_dependent_conductivity_follows_native_backward_euler_storage() {
    // Unit cube: six half-cell conduction paths give 12*k(T_n) W/K.
    // Each endpoint has the independent closed-form BE solution below.
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    setup.power_w = vec![12.0];
    let unsteady = UnsteadyConfig::new(0.05, 8);
    let mut queries = Vec::new();
    let result = fluid
        .march_natural_convection(
            &domain,
            &[],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &setup,
            &EnergyConfig::default(),
            &[300.0],
            limits([0.0; 3]),
            &unsteady,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |cell, t, _| {
                assert_eq!(cell, 0);
                queries.push(t);
                FluidTransportState::try_new(
                    1.0 - 0.001 * (t - 300.0).powi(2),
                    2.0,
                    1.0 + 0.1 * (t - 300.0),
                    1.0 / (1.0 + 0.05 * (t - 300.0)),
                )
            },
        )
        .unwrap();
    assert_eq!(queries.len(), unsteady.steps + 1);
    let mut expected = 300.0;
    for (step, &query) in queries.iter().enumerate() {
        close(query, expected, 1e-12);
        if step == unsteady.steps {
            break;
        }
        let k = 1.0 + 0.1 * (expected - 300.0);
        expected = (2.0 / unsteady.time_step_s * expected + 12.0 * k * 300.0 + 12.0)
            / (2.0 / unsteady.time_step_s + 12.0 * k);
    }
    close(result.march.temperature[0], expected, 1e-12);
    close(result.applied_time_s, 0.35, 1e-14);
    close(result.endpoint_time_s, 0.4, 1e-14);
    assert_eq!(result.applied[0].temperature_k, queries[unsteady.steps - 1]);
    assert_eq!(
        result.endpoint[0].temperature_k,
        result.march.temperature[0]
    );
    assert_eq!(result.reference, reference);
    for record in &result.march.energy_records {
        close(record.source_j, 0.6, 1e-12);
        assert!(record.closure_j.abs() < 1e-9, "{record:?}");
    }
}

#[test]
fn g1_signed_source_buoyancy_viscosity_and_darcy_drag_preserve_be_and_bdf_history() {
    // Uniform pressure-open Darcy column: rho*dU/dt + mu(T)*U/K = F(T).
    // Thermostats jump from 300 to 310 K in step one. The first force is
    // zero; every later native momentum step has the exact recurrence below.
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
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
    let mut setup = ThermalSetup::new(faces);
    setup.fixed_temperature = (0..4).map(|cell| (cell, 310.0)).collect();
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let mut unsteady = UnsteadyConfig::new(0.25, 6);
        unsteady.scheme = scheme;
        unsteady.inner_tolerance = 1e-11;
        unsteady.inner_iterations = 1000;
        unsteady.probe = Some(domain.index(1, 0, 0));
        for (expansion, viscosity_scale) in [(0.0, 1.0), (0.001, 1.0), (-0.001, 1.0), (0.001, 2.0)]
        {
            let result = fluid
                .march_natural_convection(
                    &domain,
                    &[],
                    &flow,
                    &setup,
                    &EnergyConfig::default(),
                    &[300.0; 4],
                    limits([-9.81, 0.0, 0.0]),
                    &unsteady,
                    |_| 1.0,
                    |_| 1.0,
                    None,
                    &CancelGate::new(),
                    |_, t, _| {
                        FluidTransportState::try_new(
                            1.0 - expansion * (t - 300.0).powi(2) / 10.0,
                            2.0,
                            1.0 + 0.01 * (t - 300.0),
                            viscosity_scale * (1.0 + 0.1 * (t - 300.0)),
                        )
                    },
                )
                .unwrap();
            let force = 9.81 * expansion * 10.0;
            let mu = 2.0 * viscosity_scale;
            let (mut previous, mut older) = (0.0, 0.0);
            for (step, record) in result.march.flow.records.iter().enumerate() {
                let expected = if step == 0 {
                    0.0
                } else {
                    match scheme {
                        TimeScheme::BackwardEuler => {
                            (previous + unsteady.time_step_s * force)
                                / (1.0 + unsteady.time_step_s * mu)
                        }
                        TimeScheme::Bdf2 => {
                            (4.0 * previous - older + 2.0 * unsteady.time_step_s * force)
                                / (3.0 + 2.0 * unsteady.time_step_s * mu)
                        }
                    }
                };
                close(record.probe_velocity_m_s[0], expected, 2e-8);
                assert!(record.mass_residual <= unsteady.inner_tolerance);
                older = previous;
                previous = expected;
            }
            close(result.march.flow.flow.field.flux_x(2, 0, 0), previous, 2e-8);
            for sample in &result.applied {
                close(sample.force_n_m3[0], force, 1e-13);
                close(
                    sample
                        .material
                        .resolved_transport()
                        .dynamic_viscosity_pa_s(),
                    mu,
                    1e-13,
                );
                assert_eq!(
                    sample
                        .material
                        .reference_storage_transport()
                        .density_kg_m3(),
                    1.0
                );
            }
            for record in &result.march.energy_records {
                assert!(record.closure_j.abs() < 1e-8, "{record:?}");
            }
        }
    }
}

#[test]
fn g1_radiation_and_fixed_cell_reservoir_storage_retain_original_heat_accounting() {
    let dx = 0.01;
    let domain = VoxelDomain::from_fn(1, 1, 2, dx, |p| {
        if p[2] < dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let reference = FluidTransportState::try_new(1.0, 1.0, 1e-12, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let solids = [SolidMaterial::new("held plate", 200.0).with_heat_capacity(2.4e6)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 400.0), (0, 450.0)];
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let radiation = RadiationConfig::new(vec![0.9], surroundings);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let unsteady = UnsteadyConfig::new(10.0, 2);
    let energy = EnergyConfig::default();
    let gate = CancelGate::new();
    let original = march_conjugate(
        &domain,
        fluid.properties(),
        &solids,
        &setup,
        &flow,
        &unsteady,
        None,
        &[300.0; 2],
        |_| 1.0,
        |_| 1.0,
        &energy,
        Some(&radiation),
        &gate,
    )
    .unwrap();
    let result = fluid
        .march_natural_convection(
            &domain,
            &solids,
            &flow,
            &setup,
            &energy,
            &[300.0; 2],
            limits([0.0; 3]),
            &unsteady,
            |_| 1.0,
            |_| 1.0,
            Some(&radiation),
            &gate,
            |_, _, _| Ok(reference),
        )
        .unwrap();
    assert_eq!(result.march.temperature, original.temperature);
    assert_eq!(result.march.energy_records, original.energy_records);
    let expected_radiation =
        0.9 * STEFAN_BOLTZMANN * dx * dx * (450.0_f64.powi(4) - 300.0_f64.powi(4));
    for record in &result.march.energy_records {
        close(record.radiated_w, expected_radiation, 1e-10);
        assert!(record.closure_j.abs() < 1e-7, "{record:?}");
    }
    close(
        result.march.energy_records[0].stored_energy_change_j,
        2.4e6 * domain.volume(0) * 150.0,
        1e-6,
    );
}

#[test]
fn g0_final_source_support_departure_and_cancellation_refuse_without_partial_output() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![4.0];
    let before = setup.clone();
    let initial = [300.0];
    let unsteady = UnsteadyConfig::new(0.5, 1);
    let energy = EnergyConfig::default();
    for source_refusal in [false, true] {
        let result = fluid.march_natural_convection(
            &domain,
            &[],
            &flow,
            &setup,
            &energy,
            &initial,
            limits([0.0, -9.81, 0.0]),
            &unsteady,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| {
                if t > 300.5 && source_refusal {
                    return Err(MaterialError::State {
                        what: "synthetic source ends at 300.5 K".into(),
                    });
                }
                FluidTransportState::try_new(if t > 300.5 { 0.9 } else { 1.0 }, 2.0, 1.0, 1.0)
            },
        );
        assert!(matches!(result, Err(MaterialBuoyancyError::Material {
            cell: 0, temperature_k, ..
        }) if temperature_k > 300.5));
    }
    let gate = CancelGate::new();
    let cancelled = fluid.march_natural_convection(
        &domain,
        &[],
        &flow,
        &setup,
        &energy,
        &initial,
        limits([0.0; 3]),
        &unsteady,
        |_| 1.0,
        |_| 1.0,
        None,
        &gate,
        |_, _, _| {
            gate.request();
            Ok(reference)
        },
    );
    assert!(matches!(
        cancelled,
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));
    assert_eq!(setup, before);
    assert_eq!(initial, [300.0]);
}
