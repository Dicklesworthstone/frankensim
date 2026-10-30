//! G0/G1 checks for material-driven unsteady plane-channel momentum.
//! The continuum reference is the independent sine-series solution, not a second
//! call to the production matrix solver or a fitted flow/inertance coefficient.

use fs_couple::bernoulli_aperture::slit_channel::{
    SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
};
use fs_material::gas::{GasSpec, GasState};

const LENGTH: f64 = 0.015;
const HEIGHT: f64 = 0.001;
const WIDTH: f64 = 0.01;

fn gas(viscosity_scale: f64, pressure_pa: f64) -> GasState {
    let mut spec = GasSpec::dry_air_ussa1976();
    // Deliberately synthetic transport intervention; it is not a named measured gas.
    spec.sutherland_beta *= viscosity_scale;
    GasState::try_new(&spec, 293.15, pressure_pa).unwrap()
}

fn channel(g: GasState, cells: usize, dt: f64, exit: SlitExit) -> SlitChannel {
    SlitChannel::from_gas(&g, SlitChannelSpec { length_m: LENGTH, width_m: WIDTH,
        height_m: HEIGHT, transverse_cells: cells, sample_period_s: dt, exit }).unwrap()
}

fn viscous_time(g: GasState) -> f64 {
    HEIGHT.powi(2) * g.density / (g.dynamic_viscosity * core::f64::consts::PI.powi(2))
}

fn continuum_mean(g: GasState, pressure: f64, time: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let sum: f64 = (1..512).step_by(2).map(|n| {
        let n = f64::from(n);
        (-n * n * time / viscous_time(g)).exp() / n.powi(4)
    }).sum();
    pressure * HEIGHT.powi(2) / (12.0 * g.dynamic_viscosity * LENGTH)
        * (1.0 - 96.0 * sum / pi.powi(4))
}

fn check_energy(frame: SlitChannelStep) {
    let scale = frame.pressure_work_j.abs() + frame.stored_delta_j.abs()
        + frame.viscous_dissipation_j + frame.exit_work_j + frame.numerical_dissipation_j;
    assert!(frame.energy_residual_j.abs() <= 1e-12 * scale,
        "unclosed channel energy {frame:?}");
    assert!(frame.kinetic_energy_j >= 0.0 && frame.viscous_dissipation_j >= 0.0
        && frame.exit_work_j >= 0.0 && frame.numerical_dissipation_j >= 0.0);
}

#[test]
fn no_slip_profile_and_flow_converge_to_continuum_poiseuille() {
    let g = gas(1.0, 101_325.0);
    let pressure = 100.0;
    let mut errors = Vec::new();
    for cells in [8, 16, 32] {
        let mut c = channel(g, cells, viscous_time(g) / 20.0, SlitExit::StaticPressure);
        for _ in 0..600 { check_energy(c.step(pressure).unwrap()); }
        let continuum = pressure * HEIGHT.powi(2) / (12.0 * g.dynamic_viscosity * LENGTH);
        errors.push((c.last_step().bulk_speed_m_s - continuum).abs() / continuum);
        assert!((c.last_step().bulk_speed_m_s - c.steady_bulk_speed_m_s(pressure).unwrap())
            .abs() / continuum < 1e-10);
        let dy = HEIGHT / cells as f64;
        for (i, &v) in c.velocities_m_s().iter().enumerate() {
            let y = (i as f64 + 0.5) * dy;
            let exact = pressure * y * (HEIGHT - y) / (2.0 * g.dynamic_viscosity * LENGTH);
            // The explicit O(dy²) wall-location error converges with the mean.
            assert!((v - exact).abs() < 0.2 * pressure * dy.powi(2)
                / (g.dynamic_viscosity * LENGTH));
        }
    }
    assert!(errors[0] > 0.02, "the h-axis must be live");
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((1.95..2.05).contains(&order), "spatial errors {errors:?}");
    }
}

#[test]
fn startup_momentum_has_independent_time_ladder_against_unsteady_sine_series() {
    let g = gas(1.0, 101_325.0);
    let duration = 0.6 * viscous_time(g);
    let exact = continuum_mean(g, 100.0, duration);
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut c = channel(g, 512, duration / f64::from(steps), SlitExit::StaticPressure);
        for _ in 0..steps { c.step(100.0).unwrap(); }
        errors.push((c.last_step().bulk_speed_m_s - exact).abs() / exact);
    }
    assert!(errors[0] > 0.01, "the time axis must be live");
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((0.9..1.1).contains(&order), "backward-Euler time errors {errors:?}");
    }
    let mut early = channel(g, 128, 0.01 * viscous_time(g), SlitExit::StaticPressure);
    let f = early.step(100.0).unwrap();
    assert!(early.velocities_m_s()[64] / f.bulk_speed_m_s < 1.15,
        "startup must resolve a flat core, not impose a parabolic profile");
    assert!(f.bulk_speed_m_s < 0.02 * early.steady_bulk_speed_m_s(100.0).unwrap(),
        "finite source inertia must survive a pressure step");
}

#[test]
fn implicit_exit_and_resolved_storage_close_through_release_and_flow_reversal() {
    let g = gas(1.0, 101_325.0);
    let mut c = channel(g, 64, 1.0 / 48_000.0, SlitExit::UniformKineticHead);
    let mut previous_energy = 0.0;
    let mut positive = false;
    let mut negative = false;
    let mut release = false;
    for n in 0..1800 {
        let p = if n < 500 { 300.0 } else if n < 800 { 0.0 }
            else if n < 1300 { -200.0 } else { 0.0 };
        let f = c.step(p).unwrap();
        check_energy(f);
        assert!((f.kinetic_energy_j - previous_energy - f.stored_delta_j).abs()
            < 1e-12 * (f.kinetic_energy_j + previous_energy + 1e-100));
        if p == 0.0 {
            assert_eq!(f.pressure_work_j, 0.0);
            assert!(f.kinetic_energy_j <= previous_energy);
            if n == 500 {
                assert!(f.bulk_speed_m_s > 1.0 && f.viscous_dissipation_j > 0.0);
                assert!(f.exit_work_j > 0.0 && f.numerical_dissipation_j > 0.0);
                release = true;
            }
        }
        positive |= f.bulk_speed_m_s > 1.0;
        negative |= f.bulk_speed_m_s < -1.0;
        previous_energy = f.kinetic_energy_j;
    }
    assert!(positive && negative && release);
}

#[test]
fn gas_data_controls_mass_and_wall_diffusion_without_channel_code_changes() {
    let a = gas(1.0, 101_325.0);
    let b = gas(2.0, 101_325.0);
    let dense = gas(1.0, 202_650.0);
    assert_eq!(a.density, b.density);
    assert_eq!(a.sound_speed, b.sound_speed);
    assert_eq!(a.dynamic_viscosity, dense.dynamic_viscosity);
    let mut low = channel(a, 64, 1.0 / 48_000.0, SlitExit::StaticPressure);
    let mut high = channel(b, 64, 1.0 / 48_000.0, SlitExit::StaticPressure);
    let mut heavy = channel(dense, 64, 1.0 / 48_000.0, SlitExit::StaticPressure);
    assert_eq!(low.gas(), &a);
    let light_first = low.step(100.0).unwrap();
    let heavy_first = heavy.step(100.0).unwrap();
    assert!(heavy_first.bulk_speed_m_s < 0.6 * light_first.bulk_speed_m_s,
        "gas density must control transient inertia");
    high.step(100.0).unwrap();
    for _ in 0..200 { low.step(100.0).unwrap(); high.step(100.0).unwrap(); }
    assert!(high.last_step().bulk_speed_m_s < 0.9 * low.last_step().bulk_speed_m_s,
        "gas transport must control the actual evolving velocity field");
    assert_eq!(low.spec(), high.spec());
}

#[test]
fn invalid_geometry_or_pressure_cannot_mutate_accepted_channel_state() {
    let g = gas(1.0, 101_325.0);
    let mut c = channel(g, 16, 1.0 / 48_000.0, SlitExit::UniformKineticHead);
    c.step(300.0).unwrap();
    let mut untouched = c.clone();
    for p in [f64::NAN, f64::INFINITY, -f64::INFINITY, f64::MAX] {
        assert!(c.step(p).is_err());
        assert_eq!(c.velocities_m_s(), untouched.velocities_m_s());
        assert_eq!(c.last_step(), untouched.last_step());
    }
    assert_eq!(c.step(100.0).unwrap(), untouched.step(100.0).unwrap());
    for cells in [0, 1, 4097] {
        let mut spec = *c.spec(); spec.transverse_cells = cells;
        assert!(SlitChannel::from_gas(&g, spec).is_err());
    }
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut spec = *c.spec(); spec.height_m = value;
        assert!(SlitChannel::from_gas(&g, spec).is_err());
    }
    let mut inconsistent = g; inconsistent.density *= 1.1;
    assert!(SlitChannel::from_gas(&inconsistent, *c.spec()).is_err());
}

#[test]
fn affine_flow_condensation_retains_the_same_profile_pressure_and_energy() {
    let g = gas(1.0, 101_325.0);
    let mut c = channel(g, 32, 1.0 / 48_000.0, SlitExit::StaticPressure);
    let initial_impedance = c.pressure_load().unwrap().impedance_pa_s_m3;
    assert!(initial_impedance > 0.0);
    assert_eq!(c.pressure_load().unwrap().history_pressure_pa, 0.0);
    for pressure in [100.0, 200.0, 0.0, -75.0, 50.0] {
        let old_profile = c.velocities_m_s().to_vec();
        let old_frame = c.last_step();
        let load = c.pressure_load().unwrap();
        assert_eq!(load.impedance_pa_s_m3, initial_impedance);
        let mut driven = c.clone();
        let expected = driven.step(pressure).unwrap();
        assert!((load.pressure_drop_pa(expected.flow_m3_s).unwrap() - pressure).abs() < 1e-10);
        let candidate = c.preview_at_flow(expected.flow_m3_s).unwrap();
        check_energy(candidate.last_step());
        assert!((candidate.last_step().flow_m3_s - expected.flow_m3_s).abs() < 1e-16);
        for (&a, &b) in candidate.velocities_m_s().iter().zip(driven.velocities_m_s()) {
            assert!((a-b).abs() < 1e-12);
        }
        assert_eq!(c.velocities_m_s(), old_profile.as_slice());
        assert_eq!(c.last_step(), old_frame);
        assert!(c.preview_at_flow(f64::NAN).is_err());
        assert_eq!(c.last_step(), old_frame);
        c = candidate;
    }
    let nonlinear = channel(g, 32, 1.0 / 48_000.0, SlitExit::UniformKineticHead);
    assert!(nonlinear.pressure_load().is_err());
    assert!(nonlinear.preview_at_flow(1e-5).is_err());
}

#[test]
fn centre_velocity_samples_the_actual_even_and_odd_mesh_and_steady_operator() {
    let g = gas(1.0, 101_325.0);
    let pressure = 100.0;
    let exact = pressure * HEIGHT.powi(2) / (8.0 * g.dynamic_viscosity * LENGTH);
    for cells in [15, 16, 31, 32] {
        for exit in [SlitExit::StaticPressure, SlitExit::UniformKineticHead] {
            let mut c = channel(g, cells, viscous_time(g) / 20.0, exit);
            assert_eq!(c.centre_speed_m_s(), 0.0);
            let steady = c.steady_centre_speed_m_s(pressure).unwrap();
            assert_eq!(c.steady_centre_speed_m_s(-pressure).unwrap(), -steady);
            assert_eq!(c.steady_centre_speed_m_s(0.0).unwrap(), 0.0);
            for _ in 0..600 { check_energy(c.step(pressure).unwrap()); }
            let field = c.velocities_m_s();
            let middle = cells / 2;
            let expected = if cells % 2 == 1 { field[middle] }
                else { 0.5 * field[middle - 1] + 0.5 * field[middle] };
            assert_eq!(c.centre_speed_m_s(), expected);
            assert!((c.centre_speed_m_s() - steady).abs() < 1e-10 * steady);
            if exit == SlitExit::StaticPressure {
                // The cell-centred no-slip boundary and reconstruction retain
                // their spatial error; no continuum ratio repairs the field.
                assert!((steady - exact).abs() / exact < 2.0 / (cells * cells) as f64);
            }
            let saved = c.velocities_m_s().to_vec();
            assert!(c.steady_centre_speed_m_s(f64::NAN).is_err());
            assert_eq!(c.velocities_m_s(), saved);
        }
    }
}

#[test]
fn centre_velocity_keeps_startup_release_and_reversal_momentum() {
    let g = gas(1.0, 101_325.0);
    for cells in [64, 65] {
        let mut c = channel(g, cells, 0.005 * viscous_time(g), SlitExit::StaticPressure);
        let first = c.step(100.0).unwrap();
        let centre = c.centre_speed_m_s();
        assert!(centre > first.bulk_speed_m_s);
        assert!(centre / first.bulk_speed_m_s < 1.15,
            "a starting flat core must not be replaced by 1.5 times bulk flow");
        let release = c.step(0.0).unwrap();
        assert!(c.centre_speed_m_s() > 0.0);
        assert!(release.kinetic_energy_j < first.kinetic_energy_j);
        assert_eq!(release.pressure_work_j, 0.0);
        check_energy(release);
        for _ in 0..64 { check_energy(c.step(-100.0).unwrap()); }
        assert!(c.centre_speed_m_s() < 0.0, "channel observation retains signed reversal");
    }
}
