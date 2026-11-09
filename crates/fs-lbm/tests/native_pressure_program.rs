//! Retained pressure forcing in the actual native flow, heat and material owners.
//! Synthetic constant/temperature-dependent laws are discrete-equation oracles.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConjugateMarch, ConvectionScheme, EnergyConfig, FlowResistance, FluidProperties,
    FvBoundary, MaterialEndpointConfig, PressureProgram, SimpleConfig, ThermalFace, ThermalSetup,
    TimeScheme, UnsteadyConfig, VoxelDomain, march_conjugate,
    simple_unsteady_with_pressure_schedule,
};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::fluid_transport::FluidTransportState;

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.02,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn channel(axis: usize, offset: f64) -> (VoxelDomain, SimpleConfig, ThermalSetup) {
    let mut dims = [1; 3];
    dims[axis] = 4;
    let domain = VoxelDomain::new(dims[0], dims[1], dims[2], 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * axis] = FvBoundary::Pressure {
        pressure_pa: offset + 1.0,
    };
    faces[2 * axis + 1] = FvBoundary::Pressure {
        pressure_pa: offset,
    };
    let mut flow = SimpleConfig::new(faces);
    flow.pressure_tolerance = 1e-12;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[2 * axis] = ThermalFace::Outflow {
        backflow_temperature: 320.0,
    };
    faces[2 * axis + 1] = ThermalFace::Outflow {
        backflow_temperature: 280.0,
    };
    (domain, flow, ThermalSetup::new(faces))
}

fn clock(scheme: TimeScheme) -> UnsteadyConfig {
    let mut time = UnsteadyConfig::new(0.125, 6);
    time.scheme = scheme;
    time.inner_iterations = 1000;
    time.inner_tolerance = 1e-11;
    time.average_from_step = 1;
    time.probe = Some(1);
    time.pressure_program = PressureProgram::new(&[(0.0, 1.0), (0.25, 1.0), (0.5, -4.0)]).unwrap();
    time
}

fn energy() -> EnergyConfig {
    EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        max_iterations: 2000,
        tolerance: 1e-12,
        ..EnergyConfig::default()
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn run(
    domain: &VoxelDomain,
    flow: &SimpleConfig,
    thermal: &ThermalSetup,
    time: &UnsteadyConfig,
    gate: &CancelGate,
) -> Result<ConjugateMarch, ChtError> {
    march_conjugate(
        domain,
        &fluid(),
        &[],
        thermal,
        flow,
        time,
        None,
        &[300.0; 4],
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        gate,
    )
}

// A fluid-only march deliberately records NaN for its absent solid maximum.
// Compare every heat diagnostic bitwise so that this sentinel cannot make
// otherwise identical pressure-program results fail floating PartialEq.
fn assert_same_thermal(a: &ConjugateMarch, b: &ConjugateMarch) {
    assert_eq!(a.temperature, b.temperature);
    assert_eq!(a.mean_temperature, b.mean_temperature);
    assert_eq!(a.energy_records.len(), b.energy_records.len());
    for (a, b) in a.energy_records.iter().zip(&b.energy_records) {
        let bits = |r: &fs_lbm::conjugate::TransientRecord| {
            [
                r.time_s,
                r.max_temperature_k,
                r.max_solid_temperature_k,
                r.stored_energy_change_j,
                r.source_j,
                r.boundary_outflow_j,
                r.closure_j,
                r.radiated_w,
            ]
            .map(f64::to_bits)
        };
        assert_eq!(bits(a), bits(b));
        assert_eq!(a.iterations, b.iterations);
    }
}

// Independent four-cell heat equation: C=2 J/K, k*A/dx=.02 W/K,
// endpoint upwind advection and no conductive exchange at pressure openings.
// Thomas elimination operates only this explicit reference tridiagonal system.
fn heat_step(previous: [f64; 4], velocity: f64, dt: f64) -> [f64; 4] {
    let capacity_rate = 2.0 / dt;
    let advective = 2.0 * velocity.abs();
    let inlet = if velocity >= 0.0 { 320.0 } else { 280.0 };
    let mut lower = [0.0; 4];
    let mut upper = [0.0; 4];
    let mut diagonal = [capacity_rate + advective; 4];
    let mut rhs = previous.map(|t| capacity_rate * t);
    for cell in 0..4 {
        if cell > 0 {
            lower[cell] = -0.02 - if velocity > 0.0 { advective } else { 0.0 };
            diagonal[cell] += 0.02;
        }
        if cell < 3 {
            upper[cell] = -0.02 - if velocity < 0.0 { advective } else { 0.0 };
            diagonal[cell] += 0.02;
        }
    }
    rhs[if velocity >= 0.0 { 0 } else { 3 }] += advective * inlet;
    for cell in 1..4 {
        let factor = lower[cell] / diagonal[cell - 1];
        diagonal[cell] -= factor * upper[cell - 1];
        rhs[cell] -= factor * rhs[cell - 1];
    }
    let mut temperature = [0.0; 4];
    temperature[3] = rhs[3] / diagonal[3];
    for cell in (0..3).rev() {
        temperature[cell] = (rhs[cell] - upper[cell] * temperature[cell + 1]) / diagonal[cell];
    }
    temperature
}

#[test]
fn reversing_pressure_drives_native_heat_and_retains_be_bdf2_history_on_all_axes() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
            let (domain, flow, thermal) = channel(axis, 0.0);
            let time = clock(scheme);
            let solved = run(&domain, &flow, &thermal, &time, &gate).unwrap();
            let mut previous = 0.0;
            let mut older = 0.0;
            let mut temperature = [300.0; 4];
            let factors = [1.0, 1.0, -1.5, -4.0, -4.0, -4.0];
            for (step, &factor) in factors.iter().enumerate() {
                // Slip-side plug: rho*dU/dt=dp/L, with rho=1 and L=4.
                let impulse = time.time_step_s * factor / 4.0;
                let velocity = if step > 0 && scheme == TimeScheme::Bdf2 {
                    (2.0 * previous - 0.5 * older + impulse) / 1.5
                } else {
                    previous + impulse
                };
                let record = &solved.flow.records[step];
                close(record.probe_velocity_m_s[axis], velocity, 2e-8);
                close(record.kinetic_energy_j, 2.0 * velocity * velocity, 2e-8);
                assert!(record.mass_residual <= time.inner_tolerance);
                let next = heat_step(temperature, velocity, time.time_step_s);
                let stored = next
                    .iter()
                    .zip(&temperature)
                    .map(|(new, old)| 2.0 * (new - old))
                    .sum();
                let (outlet, inlet) = if velocity >= 0.0 {
                    (next[3], 320.0)
                } else {
                    (next[0], 280.0)
                };
                let outward = 2.0 * velocity.abs() * time.time_step_s * (outlet - inlet);
                let heat = &solved.energy_records[step];
                close(heat.stored_energy_change_j, stored, 2e-7);
                close(heat.boundary_outflow_j, outward, 2e-7);
                close(heat.closure_j, 0.0, 2e-7);
                assert_eq!(heat.source_j, 0.0);
                temperature = next;
                older = previous;
                previous = velocity;
            }
            assert!(solved.flow.records[0].probe_velocity_m_s[axis] > 0.0);
            assert!(
                previous < 0.0,
                "the opposite thermal reservoir must become the inlet"
            );
            for (&actual, expected) in solved.temperature.iter().zip(temperature) {
                close(actual, expected, 2e-9);
            }
            assert!(solved.flow.flow.report.max_divergence_m3_s < 1e-10);
            assert!(solved.flow.mean_field.max_divergence(&domain) < 1e-10);
        }
    }
}

#[test]
fn pressure_program_preserves_identity_gauge_and_callback_composition() {
    let gate = CancelGate::new_clock_free();
    let (domain, flow, thermal) = channel(0, 0.0);
    let time = clock(TimeScheme::Bdf2);
    let driven = run(&domain, &flow, &thermal, &time, &gate).unwrap();
    let (_, shifted, _) = channel(0, 1_099_511_627_776.0);
    let offset = run(&domain, &shifted, &thermal, &time, &gate).unwrap();
    assert_eq!(driven.flow.records, offset.flow.records);
    assert_eq!(driven.flow.flow.velocity_m_s, offset.flow.flow.velocity_m_s);
    assert_eq!(driven.flow.flow.field, offset.flow.flow.field);
    assert_same_thermal(&driven, &offset);

    let mut identity = time;
    identity.pressure_program = PressureProgram::IDENTITY;
    let original = run(&domain, &flow, &thermal, &identity, &gate).unwrap();
    identity.pressure_program = PressureProgram::new(&[(0.0, 1.0), (1.0, 1.0)]).unwrap();
    let equivalent = run(&domain, &flow, &thermal, &identity, &gate).unwrap();
    assert_eq!(equivalent.flow, original.flow);
    assert_same_thermal(&equivalent, &original);

    // The existing flow callback remains independently sampled exactly once.
    let calls = core::cell::Cell::new(0);
    let combined = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &flow,
        &time,
        |_| 1.0,
        |_| {
            calls.set(calls.get() + 1);
            0.5
        },
        &gate,
    )
    .unwrap();
    let expected = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &flow,
        &identity,
        |_| 1.0,
        |t| 0.5 * time.pressure_program.sample(t).unwrap(),
        &gate,
    )
    .unwrap();
    assert_eq!(calls.get(), time.steps);
    assert_eq!(combined, expected);
}

#[test]
fn pressure_program_reaches_explicit_and_implicit_source_density_and_viscosity_feedback() {
    let (domain, mut flow, mut thermal) = channel(0, 0.0);
    flow.resistances.push(FlowResistance::Volume {
        lo: [0; 3],
        hi: [4, 1, 1],
        permeability_m2: [1.0; 3],
        inertial_per_m: [0.0; 3],
    });
    thermal.fixed_temperature = (0..4).map(|cell| (cell, 310.0)).collect();
    let reference = FluidTransportState::try_new(1.0, 2.0, 0.02, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let limits = BoussinesqConfig {
        gravity_m_s2: [-9.81, 0.0, 0.0],
        max_relative_density_departure: 0.05,
        max_relative_specific_heat_departure: 0.0,
    };
    let endpoint = MaterialEndpointConfig {
        maximum_passes_per_step: 80,
        maximum_momentum_solves: 1000,
        maximum_energy_solves: 1000,
        maximum_source_refreshes: 2001,
        temperature_tolerance_k: 1e-10,
        transport_relative_tolerance: 1e-11,
        force_tolerance_n_m3: 1e-11,
        temperature_relaxation: 1.0,
    };
    let source = |_: usize, t: f64, _: &FluidTransportState| {
        let theta = t - 300.0;
        FluidTransportState::try_new(
            1.0 - 0.001 * theta,
            2.0,
            0.02 * (1.0 + 0.01 * theta),
            1.0 + 0.1 * theta,
        )
    };
    let gate = CancelGate::new_clock_free();
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let time = clock(scheme);
        for implicit in [false, true] {
            let result = if implicit {
                bound
                    .march_natural_convection_implicit(
                        &domain,
                        &[],
                        &flow,
                        &thermal,
                        &energy(),
                        &[300.0; 4],
                        limits,
                        &time,
                        endpoint,
                        |_| 1.0,
                        |_| 1.0,
                        None,
                        &gate,
                        source,
                    )
                    .unwrap()
                    .solution
            } else {
                bound
                    .march_natural_convection(
                        &domain,
                        &[],
                        &flow,
                        &thermal,
                        &energy(),
                        &[300.0; 4],
                        limits,
                        &time,
                        |_| 1.0,
                        |_| 1.0,
                        None,
                        &gate,
                        source,
                    )
                    .unwrap()
            };
            let mut previous = 0.0;
            let mut older = 0.0;
            for (step, factor) in [1.0, 1.0, -1.5, -4.0, -4.0, -4.0].into_iter().enumerate() {
                // Thermostats change the source state in interval one. Implicit
                // coupling applies new rho(T)/mu(T) immediately; explicit uses
                // the retained start epoch. Inertia still belongs to rho_ref=1.
                let hot = implicit || step > 0;
                let mu = if hot { 2.0 } else { 1.0 };
                let force = 0.25 * factor + if hot { 0.0981 } else { 0.0 };
                let dt = time.time_step_s;
                let expected = if step > 0 && scheme == TimeScheme::Bdf2 {
                    (4.0 * previous - older + 2.0 * dt * force) / (3.0 + 2.0 * dt * mu)
                } else {
                    (previous + dt * force) / (1.0 + dt * mu)
                };
                close(
                    result.march.flow.records[step].probe_velocity_m_s[0],
                    expected,
                    3e-8,
                );
                close(result.march.energy_records[step].closure_j, 0.0, 2e-7);
                older = previous;
                previous = expected;
            }
            assert!(previous < 0.0);
            assert_eq!(result.march.flow.records.len(), time.steps);
            assert_eq!(result.march.temperature, vec![310.0; 4]);
            close(
                result.march.energy_records[0].stored_energy_change_j,
                80.0,
                1e-8,
            );
            for sample in &result.applied {
                close(sample.force_n_m3[0], 0.0981, 1e-12);
                close(
                    sample
                        .material
                        .resolved_transport()
                        .dynamic_viscosity_pa_s(),
                    2.0,
                    1e-12,
                );
            }
        }
    }
    assert_eq!(bound.transport(), &reference);
}

#[test]
fn pressure_program_refuses_unrepresentable_heads_and_cancelled_marches_atomically() {
    let (domain, flow, thermal) = channel(0, 0.0);
    let before = (flow.clone(), thermal.clone());
    let mut time = clock(TimeScheme::BackwardEuler);
    time.pressure_program =
        PressureProgram::new(&[(0.0, 1.0), (0.125, 1.0), (0.25, f64::MAX)]).unwrap();
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        run(&domain, &flow, &thermal, &time, &gate),
        Err(ChtError::InvalidInput {
            field: "unsteady.pressure_difference_pa",
            ..
        })
    ));
    assert_eq!((&flow, &thermal), (&before.0, &before.1));

    time = clock(TimeScheme::Bdf2);
    gate.request();
    assert!(matches!(
        run(&domain, &flow, &thermal, &time, &gate),
        Err(ChtError::Cancelled)
    ));
    assert_eq!((&flow, &thermal), (&before.0, &before.1));

    time.time_step_s = f64::MAX;
    assert!(matches!(
        run(
            &domain,
            &flow,
            &thermal,
            &time,
            &CancelGate::new_clock_free()
        ),
        Err(ChtError::InvalidInput {
            field: "unsteady.end_time_s",
            ..
        })
    ));
}
