use super::*;
use fs_material::gas::{GasSpec, GasState};

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.3;

fn duct(inlet: DgBoundary) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        inlet,
        DgBoundary::Closed,
    )
    .unwrap()
}
fn wall(resistance: f64) -> WallPin {
    WallPin {
        surface_density: 0.2,
        stiffness_per_area: 1e7,
        resistance,
    }
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * scale.max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

/// G1: the spatially uniform sealed tube reduces independently to
/// sigma*w'' + r*w' + (K+2B/a)*w = p0. This checks actual pressure,
/// displacement, velocity, inertia, viscosity and the physical energy maps.
#[test]
fn sealed_pressure_drives_the_analytic_material_wall_oscillator() {
    for resistance in [0.0, 300.0] {
        let dg = duct(DgBoundary::Closed)
            .with_compliant_wall(wall(resistance))
            .unwrap();
        let gas = dg.gas();
        let beta = 2.0 * gas.density * gas.sound_speed.powi(2) / RADIUS;
        let pin = wall(resistance);
        let stiffness = pin.stiffness_per_area + beta;
        let damping = resistance / (2.0 * pin.surface_density);
        let omega = (stiffness / pin.surface_density - damping.powi(2)).sqrt();
        let p0 = 10.0;
        let equilibrium = p0 / stiffness;
        let dt = 2e-6;
        let propagator = dg.propagator(dt, 1e-13).unwrap();
        let mut state = dg.project_initial_state(|_| p0, |_| 0.0);
        let initial = state.stored_energy_j();
        let mut previous = initial;
        for sample in 1..=80 {
            propagator.step(&mut state, &[]).unwrap();
            propagator.sync_nodal(&mut state).unwrap();
            let t = sample as f64 * dt;
            let envelope = (-damping * t).exp();
            let (sin, cos) = (omega * t).sin_cos();
            let displacement = equilibrium * (1.0 - envelope * (cos + damping / omega * sin));
            let velocity = p0 / (pin.surface_density * omega) * envelope * sin;
            let pressure = p0 - beta * displacement;
            let observation = dg.compliant_wall_state(&state).unwrap();
            for ((&p, &w), &v) in state
                .pressure()
                .iter()
                .zip(&observation.displacement_m)
                .zip(&observation.velocity_m_s)
            {
                close(p, pressure, p0);
                close(w, displacement, equilibrium);
                close(v, velocity, p0 / (pin.surface_density * omega));
            }
            let area = core::f64::consts::TAU * RADIUS * LENGTH;
            let physical_wall_energy = 0.5
                * area
                * (pin.surface_density * velocity.powi(2)
                    + pin.stiffness_per_area * displacement.powi(2));
            close(observation.stored_energy_j, physical_wall_energy, initial);
            close(
                observation.dissipated_power_w,
                area * resistance * velocity.powi(2),
                initial / dt,
            );
            close(
                dg.wall_dissipated_power_w(&state).unwrap(),
                observation.dissipated_power_w,
                initial / dt,
            );
            if resistance == 0.0 {
                close(state.stored_energy_j(), initial, initial);
            }
            assert!(state.stored_energy_j() <= previous + 1e-12 * initial);
            previous = state.stored_energy_j();
        }
        assert!(dg.compliant_wall_state(&state).unwrap().stored_energy_j > 0.0);
    }
}

/// G0: all reciprocal entries cancel in x^T L x; only the independent
/// physical viscous power survives, including previously bound gas memories.
#[test]
fn wall_and_gas_losses_share_the_native_operator_and_pressure_ports() {
    let gas_only = duct(DgBoundary::Port(PortId(0)))
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    let original = gas_only.total_dofs();
    let dg = gas_only.clone().with_compliant_wall(wall(300.0)).unwrap();
    let n = dg.total_dofs();
    assert_eq!(n, original + 2 * dg.field_dofs());
    for i in 0..original {
        assert_eq!(
            &dg.l_tilde[i * n..i * n + original],
            &gas_only.l_tilde[i * original..(i + 1) * original]
        );
    }
    let mut state = dg.zero_state();
    for (i, coordinate) in state.tilde.iter_mut().enumerate() {
        *coordinate = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
    }
    let mut storage_rate = 0.0;
    for i in 0..n {
        for j in 0..n {
            storage_rate += state.tilde[i] * dg.l_tilde[i * n + j] * state.tilde[j];
        }
    }
    let loss = dg.wall_dissipated_power_w(&state).unwrap();
    close(storage_rate, -loss, loss);
    assert!(loss > dg.compliant_wall_state(&state).unwrap().dissipated_power_w);
    let elastic = duct(DgBoundary::Port(PortId(0)))
        .with_compliant_wall(wall(0.0))
        .unwrap();
    let propagator = elastic.propagator(2e-6, 1e-13).unwrap();
    let mut state = elastic.zero_state();
    let mut work = 0.0;
    for _ in 0..60 {
        let frame = propagator
            .step(
                &mut state,
                &[PortInput {
                    id: PortId(0),
                    effort: 1e-7,
                }],
            )
            .unwrap();
        work += frame.ledger.port_work_step_j;
    }
    close(state.stored_energy_j(), work, work);
    assert!(
        elastic
            .compliant_wall_state(&state)
            .unwrap()
            .stored_energy_j
            > 0.0
    );
}

#[test]
fn wall_attachment_refuses_foreign_state_and_unsupported_composition() {
    let original = duct(DgBoundary::Closed);
    let dg = original.clone().with_compliant_wall(wall(0.0)).unwrap();
    assert!(dg.compliant_wall_state(&original.zero_state()).is_err());
    assert!(dg.clone().with_compliant_wall(wall(0.0)).is_err());
    assert!(
        dg.clone()
            .with_viscothermal(8, [50.0, 2000.0], 0.01)
            .is_err()
    );
    assert!(
        dg.with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 3,
            band_rad_s: [300.0, 12000.0],
            maximum_relative_error: 0.01
        })
        .is_err()
    );
    for pin in [
        WallPin {
            surface_density: 0.0,
            ..wall(0.0)
        },
        WallPin {
            stiffness_per_area: f64::NAN,
            ..wall(0.0)
        },
        wall(-1.0),
    ] {
        assert!(original.clone().with_compliant_wall(pin).is_err());
    }
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: RADIUS,
            outlet_radius: 2.0 * RADIUS,
            length: LENGTH,
        }],
    };
    let tapered = DgDuct::new(
        &cone,
        original.gas(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    assert!(tapered.with_compliant_wall(wall(0.0)).is_err());
}
