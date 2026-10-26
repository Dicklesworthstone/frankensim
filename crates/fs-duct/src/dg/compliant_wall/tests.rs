use super::*;
use fs_material::gas::{GasSpec, GasState};

mod spatial;

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

/// G1: choose a true eigenmotion of the sealed gas/Maxwell wall. With
/// a=1/tau, lambda=-a/2 and K1=K+beta+sigma*a^2/4-r*a/2, the independent
/// physical solution is w=w0*exp(lambda*t), p=-beta*w, z1=-K1*w.
/// This checks memory signs, geometry normalization, loss and real propagation.
#[test]
fn maxwell_wall_follows_an_independent_relaxing_eigenmotion() {
    let base = duct(DgBoundary::Closed);
    let gas = base.gas();
    let beta = 2.0 * gas.density * gas.sound_speed.powi(2) / RADIUS;
    let pin = wall(300.0);
    let rate: f64 = 1000.0;
    let lambda = -0.5 * rate;
    let stiffness = pin.stiffness_per_area
        + beta
        + pin.surface_density * lambda * lambda
        + pin.resistance * lambda;
    let modulus = 1e6;
    let law = fs_material::visco::GeneralizedMaxwell::new(
        modulus,
        vec![(stiffness * modulus / pin.stiffness_per_area, 1.0 / rate)],
    )
    .unwrap();
    let dg = base.with_compliant_wall_relaxation(pin, &law).unwrap();
    let m = dg.field_dofs();
    let unit = dg.project_initial_state(|_| 1.0, |_| 0.0);
    let mut state = dg.zero_state();
    let w0 = 1e-7;
    // Initialize one admissible constitutive prehistory using the independent
    // continuous energy weights, not the wall observation's inverse maps.
    for i in 0..m {
        let basis = unit.tilde[i];
        state.tilde[i] = -beta * w0 * basis;
        state.tilde[2 * m + i] = (beta * pin.stiffness_per_area).sqrt() * w0 * basis;
        state.tilde[3 * m + i] = (beta * pin.surface_density).sqrt() * lambda * w0 * basis;
        state.tilde[4 * m + i] = -(beta * stiffness).sqrt() * w0 * basis;
    }
    state.initial_energy_j = state.stored_energy_j();
    state.nodal_dirty = true;
    let initial = state.initial_energy_j();
    let area = core::f64::consts::TAU * RADIUS * LENGTH;
    let physical_initial = 0.5
        * area
        * w0
        * w0
        * (beta + pin.stiffness_per_area + pin.surface_density * lambda * lambda + stiffness);
    close(initial, physical_initial, initial);
    let dt = 1e-5;
    let propagator = dg.propagator(dt, 1e-13).unwrap();
    for step in 1..=60 {
        propagator.step(&mut state, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        let envelope = (lambda * step as f64 * dt).exp();
        let w = w0 * envelope;
        let observation = dg.compliant_wall_state(&state).unwrap();
        assert_eq!(observation.relaxation_pressure_pa.len(), 1);
        for i in 0..m {
            close(state.pressure()[i], -beta * w, beta * w0);
            close(observation.displacement_m[i], w, w0);
            close(observation.velocity_m_s[i], lambda * w, -lambda * w0);
            close(
                observation.relaxation_pressure_pa[0][i],
                -stiffness * w,
                stiffness * w0,
            );
        }
        let memory_energy = 0.5 * area * stiffness * w * w;
        close(
            observation.relaxation_stored_energy_j,
            memory_energy,
            initial,
        );
        close(
            observation.stored_energy_j,
            state.wall_stored_energy_j(),
            initial,
        );
        let energy = initial * envelope * envelope;
        close(state.stored_energy_j(), energy, initial);
        close(
            observation.dissipated_power_w,
            rate * energy,
            rate * initial,
        );
        close(
            dg.wall_dissipated_power_w(&state).unwrap(),
            rate * energy,
            rate * initial,
        );
    }
}

/// G0: gas diffusion, Kelvin viscosity and every Maxwell relaxation act on the
/// original state and each contributes its own nonnegative dissipation once.
#[test]
fn maxwell_wall_preserves_gas_memory_ports_and_total_passivity() {
    let gas = duct(DgBoundary::Port(PortId(0)))
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    let pin = wall(300.0);
    let base = gas.clone().with_compliant_wall(pin).unwrap();
    let law = fs_material::visco::GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003), (4e5, 0.0005)])
        .unwrap();
    let dg = gas.with_compliant_wall_relaxation(pin, &law).unwrap();
    let old = base.total_dofs();
    let n = dg.total_dofs();
    assert_eq!(n, old + 2 * dg.field_dofs());
    for i in 0..old {
        assert_eq!(
            &dg.l_tilde[i * n..i * n + old],
            &base.l_tilde[i * old..(i + 1) * old]
        );
    }
    for ((id, column), (old_id, old_column)) in dg.ports.iter().zip(&base.ports) {
        assert_eq!(id, old_id);
        assert_eq!(&column[..old], old_column.as_slice());
        assert!(column[old..].iter().all(|value| *value == 0.0));
    }
    let mut state = dg.zero_state();
    for (i, x) in state.tilde.iter_mut().enumerate() {
        *x = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
    }
    let mut storage_rate = 0.0;
    for i in 0..n {
        for j in 0..n {
            storage_rate += state.tilde[i] * dg.l_tilde[i * n + j] * state.tilde[j];
        }
    }
    let loss = dg.wall_dissipated_power_w(&state).unwrap();
    let wall = dg.compliant_wall_state(&state).unwrap();
    close(storage_rate, -loss, loss);
    assert!(loss > wall.dissipated_power_w);
    assert!(wall.relaxation_stored_energy_j > 0.0);
    assert!(wall.stored_energy_j > wall.relaxation_stored_energy_j);
    assert!(dg.compliant_wall_state(&base.zero_state()).is_err());
}

#[test]
fn empty_maxwell_spectrum_is_exactly_elastic_and_invalid_sources_refuse() {
    use fs_material::visco::GeneralizedMaxwell;
    let elastic = duct(DgBoundary::Closed)
        .with_compliant_wall(wall(0.0))
        .unwrap();
    for terms in [vec![], vec![(0.0, 0.001)]] {
        let law = GeneralizedMaxwell::new(1e6, terms).unwrap();
        let dg = duct(DgBoundary::Closed)
            .with_compliant_wall_relaxation(wall(0.0), &law)
            .unwrap();
        assert_eq!(dg.operator_identity, elastic.operator_identity);
        assert_eq!(dg.l_tilde, elastic.l_tilde);
        assert_eq!(dg.total_dofs(), elastic.total_dofs());
        assert!(
            dg.compliant_wall_state(&elastic.zero_state())
                .unwrap()
                .relaxation_pressure_pa
                .is_empty()
        );
    }
    for (modulus, tau) in [(-1.0, 0.001), (1e6, 0.0), (0.0, f64::NAN), (1e6, 1e-320)] {
        let law = GeneralizedMaxwell {
            e_inf: 1e6,
            terms: vec![(modulus, tau)],
        };
        assert!(
            duct(DgBoundary::Closed)
                .with_compliant_wall_relaxation(wall(0.0), &law)
                .is_err()
        );
    }
    let too_many = GeneralizedMaxwell::new(1e6, vec![(1e6, 0.001); 400]).unwrap();
    assert!(
        duct(DgBoundary::Closed)
            .with_compliant_wall_relaxation(wall(0.0), &too_many)
            .is_err()
    );
}
