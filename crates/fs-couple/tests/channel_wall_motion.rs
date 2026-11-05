//! G1/G3 axial wall drive: independent Couette–Poiseuille and startup-series
//! references, momentum reaction, signed work, pressure condensation and rollback.
//! NETL MFiX FLD05 supplies the steady continuum reference. The startup series
//! follows separation of variables for rho du/dt = mu d²u/dy², u(0)=0, u(H)=U.

use fs_couple::bernoulli_aperture::slit_channel::{
    FixedChannel, SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
    circular::{CircularChannel, CircularChannelSpec},
};
use fs_material::gas::{GasSpec, GasState};

const LENGTH: f64 = 0.015;
const HEIGHT: f64 = 0.001;
const WIDTH: f64 = 0.01;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn diffusion_time(gas: GasState) -> f64 {
    HEIGHT.powi(2) * gas.density / (gas.dynamic_viscosity * core::f64::consts::PI.powi(2))
}

fn slit(gas: GasState, cells: usize, dt: f64) -> SlitChannel {
    SlitChannel::from_gas(&gas, SlitChannelSpec {
        length_m: LENGTH, width_m: WIDTH, height_m: HEIGHT,
        transverse_cells: cells, sample_period_s: dt, exit: SlitExit::StaticPressure,
    }).unwrap()
}

fn energy(step: SlitChannelStep) {
    let scale = step.pressure_work_j.abs() + step.wall_work_j.abs()
        + step.stored_delta_j.abs() + step.viscous_dissipation_j
        + step.exit_work_j + step.numerical_dissipation_j;
    assert!(step.energy_residual_j.abs() <= 2e-12 * scale, "{step:?}");
    assert!(step.viscous_dissipation_j >= 0.0 && step.numerical_dissipation_j >= 0.0);
}

#[test]
fn couette_profile_reaction_and_heat_follow_the_actual_two_wall_speeds() {
    let gas = gas();
    let dt = diffusion_time(gas) / 10.0;
    let mut channel = slit(gas, 32, dt);
    let walls = [-0.7, 1.3];
    for _ in 0..500 { energy(channel.step_with_wall_velocities(0.0, walls).unwrap()); }
    for (i, &actual) in channel.velocities_m_s().iter().enumerate() {
        let exact = walls[0] + (walls[1] - walls[0]) * (i as f64 + 0.5) / 32.0;
        assert!((actual - exact).abs() < 1e-12, "cell {i}: {actual} versus {exact}");
    }
    let step = channel.last_step();
    let force = gas.dynamic_viscosity * LENGTH * WIDTH * (walls[1] - walls[0]) / HEIGHT;
    assert!((step.wall_impulse_n_s[0] / dt + force).abs() < 1e-12 * force);
    assert!((step.wall_impulse_n_s[1] / dt - force).abs() < 1e-12 * force);
    let heat = force * (walls[1] - walls[0]) * dt;
    assert!((step.viscous_dissipation_j / heat - 1.0).abs() < 1e-12);
    assert!((step.wall_work_j / heat - 1.0).abs() < 1e-12);
    assert_eq!(step.pressure_work_j, 0.0);
}

#[test]
fn adverse_pressure_and_wall_shear_converge_to_couette_poiseuille() {
    let gas = gas();
    let pressure = -0.4;
    let walls = [0.0, 0.2];
    let exact_mean = 0.5 * (walls[0] + walls[1])
        + pressure * HEIGHT.powi(2) / (12.0 * gas.dynamic_viscosity * LENGTH);
    let mut errors = Vec::new();
    for cells in [8, 16, 32] {
        let mut channel = slit(gas, cells, diffusion_time(gas) / 10.0);
        for _ in 0..500 {
            energy(channel.step_with_wall_velocities(pressure, walls).unwrap());
        }
        errors.push((channel.last_step().bulk_speed_m_s - exact_mean).abs());
        assert!(channel.velocities_m_s().iter().any(|v| *v < 0.0));
        assert!(channel.velocities_m_s().iter().any(|v| *v > 0.0));
    }
    assert!(errors[0] > 1e-3, "spatial axis must be live: {errors:?}");
    for pair in errors.windows(2) {
        assert!((1.95..2.05).contains(&(pair[0] / pair[1]).log2()), "{errors:?}");
    }
}

#[test]
fn wall_startup_converges_to_sine_series_and_accounts_for_fluid_impulse() {
    let gas = gas();
    let duration = 0.6 * diffusion_time(gas);
    let upper_speed = 0.8;
    let sum: f64 = (1..256).step_by(2).map(|n| {
        let n = f64::from(n);
        (-n * n * duration / diffusion_time(gas)).exp() / (n * n)
    }).sum();
    let exact = upper_speed * (0.5 - 4.0 * sum / core::f64::consts::PI.powi(2));
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let cells = 512;
        let mut channel = slit(gas, cells, duration / f64::from(steps));
        let cell_mass = gas.density * LENGTH * WIDTH * HEIGHT / cells as f64;
        for _ in 0..steps {
            let old_momentum = cell_mass * channel.velocities_m_s().iter().sum::<f64>();
            let step = channel.step_with_wall_velocities(0.0, [0.0, upper_speed]).unwrap();
            energy(step);
            let momentum = cell_mass * channel.velocities_m_s().iter().sum::<f64>();
            let impulse = step.wall_impulse_n_s.iter().sum::<f64>();
            assert!((momentum - old_momentum - impulse).abs()
                < 2e-11 * (momentum.abs() + old_momentum.abs()));
        }
        errors.push((channel.last_step().bulk_speed_m_s - exact).abs() / exact);
    }
    assert!(errors[0] > 0.01, "time axis must be live: {errors:?}");
    for pair in errors.windows(2) {
        assert!((0.9..1.1).contains(&(pair[0] / pair[1]).log2()), "{errors:?}");
    }
}

#[test]
fn moving_wall_can_recover_work_while_relative_shear_stays_dissipative() {
    let gas = gas();
    let mut channel = slit(gas, 64, diffusion_time(gas) / 20.0);
    for _ in 0..600 { channel.step(20.0).unwrap(); }
    let step = channel.step_with_wall_velocities(0.0, [0.1, 0.1]).unwrap();
    energy(step);
    assert!(step.wall_work_j < 0.0, "slower walls must receive mechanical work");
    assert!(step.viscous_dissipation_j > 0.0 && step.stored_delta_j < 0.0);
    assert!(step.wall_impulse_n_s.iter().all(|impulse| *impulse < 0.0));
}

#[test]
fn wall_drive_pays_for_the_same_implicit_kinetic_head_terminal() {
    let gas = gas();
    let mut channel = SlitChannel::from_gas(&gas, SlitChannelSpec {
        length_m: LENGTH, width_m: WIDTH, height_m: HEIGHT, transverse_cells: 32,
        sample_period_s: diffusion_time(gas) / 20.0, exit: SlitExit::UniformKineticHead,
    }).unwrap();
    for speed in [0.8; 60].into_iter().chain([-0.8; 120]) {
        let frame = channel.step_with_wall_velocities(0.0, [0.0, speed]).unwrap();
        energy(frame);
        assert_eq!(frame.pressure_work_j, 0.0);
        assert!(frame.exit_work_j >= 0.0);
        assert!(frame.viscous_dissipation_j > 0.0);
    }
    assert!(channel.last_step().flow_m3_s < 0.0);
    assert!(channel.last_step().exit_work_j > 0.0);
    assert!(channel.pressure_load_with_wall_velocities([0.0, 0.8]).is_err(),
        "a nonlinear terminal cannot silently become an affine pressure load");
}

#[test]
fn both_geometries_condense_wall_drive_into_the_same_pressure_port_and_refuse_atomically() {
    let gas = gas();
    let dt = diffusion_time(gas) / 20.0;
    let round = CircularChannel::from_gas(&gas, CircularChannelSpec {
        length_m: LENGTH, radius_m: HEIGHT / 2.0, radial_cells: 32,
        sample_period_s: dt, exit: SlitExit::StaticPressure,
    }).unwrap();
    for mut channel in [FixedChannel::from(slit(gas, 32, dt)), FixedChannel::from(round)] {
        for _ in 0..20 { channel.step(0.3).unwrap(); }
        let before = channel.velocities_m_s().to_vec();
        let last = channel.last_step();
        let walls = [0.0, 0.8];
        let flow = -1e-7;
        let load = channel.pressure_load_with_wall_velocities(walls).unwrap();
        let trial = channel.preview_at_flow_with_wall_velocities(flow, walls).unwrap();
        let step = trial.last_step();
        energy(step);
        assert!((step.flow_m3_s - flow).abs() < 1e-19);
        assert!((step.pressure_work_j - load.pressure_drop_pa(flow).unwrap()
            * step.flow_m3_s * dt).abs() < 1e-24);
        assert!(step.wall_work_j > 0.0);
        assert_eq!(channel.velocities_m_s(), before);
        assert_eq!(channel.last_step(), last);
        for walls in [[0.0, f64::NAN], [f64::INFINITY, 0.0], [0.0, f64::MAX]] {
            assert!(channel.step_with_wall_velocities(0.0, walls).is_err());
            assert_eq!(channel.velocities_m_s(), before);
            assert_eq!(channel.last_step(), last);
        }
        if channel.circular().is_some() {
            assert_eq!(step.wall_impulse_n_s[0], 0.0);
            assert!(channel.step_with_wall_velocities(0.0, [1.0, 0.0]).is_err());
        }
        let mut stationary = channel.clone();
        assert_eq!(channel.step_with_wall_velocities(0.2, [0.0; 2]).unwrap(),
            stationary.step(0.2).unwrap());
    }
}

#[test]
fn translating_circular_wall_drives_plug_translation_after_viscous_startup() {
    let gas = gas();
    let dt = diffusion_time(gas) / 20.0;
    let mut channel = CircularChannel::from_gas(&gas, CircularChannelSpec {
        length_m: LENGTH, radius_m: HEIGHT / 2.0, radial_cells: 32,
        sample_period_s: dt, exit: SlitExit::StaticPressure,
    }).unwrap();
    let first = channel.step_with_wall_velocities(0.0, [0.0, 0.5]).unwrap();
    energy(first);
    assert!(channel.velocities_m_s()[0] < 0.01 * channel.velocities_m_s()[31],
        "wall impulse must diffuse inward, never impose plug flow at startup");
    for _ in 0..1000 { channel.step_with_wall_velocities(0.0, [0.0, 0.5]).unwrap(); }
    assert!(channel.velocities_m_s().iter().all(|v| (*v - 0.5).abs() < 1e-12));
    assert!(channel.last_step().viscous_dissipation_j < 1e-25);
}
