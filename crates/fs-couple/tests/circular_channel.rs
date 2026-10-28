//! G1/G3 circular momentum against independent Poiseuille and Fourier–Bessel
//! solutions, with material intervention, signed work and transactional refusal.
use fs_couple::bernoulli_aperture::slit_channel::{FixedChannel, SlitExit};
use fs_couple::bernoulli_aperture::slit_channel::circular::{CircularChannel, CircularChannelSpec};
use fs_material::gas::{GasSpec, GasState};

const RADIUS: f64 = 0.0005;
const LENGTH: f64 = 0.02;
fn gas(viscosity: f64, pressure: f64) -> GasState {
    let mut source = GasSpec::dry_air_ussa1976();
    // Deliberately synthetic transport intervention, not a measured species.
    source.sutherland_beta *= viscosity;
    GasState::try_new(&source, 293.15, pressure).unwrap()
}
fn channel(gas: GasState, cells: usize, dt: f64) -> CircularChannel {
    CircularChannel::from_gas(&gas, CircularChannelSpec { length_m: LENGTH,
        radius_m: RADIUS, radial_cells: cells, sample_period_s: dt,
        exit: SlitExit::StaticPressure }).unwrap()
}
fn diffusion_time(gas: GasState) -> f64 {
    RADIUS.powi(2) * gas.density / (gas.dynamic_viscosity * 2.404_825_557_695_773f64.powi(2))
}

#[test]
fn circular_radial_profile_and_flux_converge_to_hagen_poiseuille() {
    let gas = gas(1.0, 101_325.0);
    let pressure = 10.0;
    let exact = pressure * RADIUS.powi(2) / (8.0 * gas.dynamic_viscosity * LENGTH);
    let mut errors = Vec::new();
    for cells in [8, 16, 32] {
        let mut channel = channel(gas, cells, diffusion_time(gas) / 20.0);
        for _ in 0..600 { channel.step(pressure).unwrap(); }
        let observed = channel.last_step().bulk_speed_m_s;
        assert!((observed - channel.steady_bulk_speed_m_s(pressure).unwrap()).abs() < 1e-10 * exact);
        errors.push((observed - exact).abs() / exact);
        let dr = RADIUS / cells as f64;
        let scale = pressure / (4.0 * gas.dynamic_viscosity * LENGTH);
        for (i, velocity) in channel.velocities_m_s().iter().enumerate() {
            let r = dr * (i as f64 + 0.5);
            let continuum = scale * (RADIUS.powi(2) - r.powi(2));
            assert!((velocity - continuum).abs() < 0.251 * scale * dr.powi(2));
        }
        let expected_flow = core::f64::consts::PI * RADIUS.powi(2) * observed;
        assert!((channel.last_step().flow_m3_s - expected_flow).abs() < 1e-14 * expected_flow);
    }
    assert!(errors[0] > 0.01, "radial axis must be live: {errors:?}");
    for e in errors.windows(2) {
        assert!((1.95..2.05).contains(&(e[0] / e[1]).log2()), "radial errors {errors:?}");
    }
}

#[test]
fn circular_startup_has_first_order_time_ladder_against_bessel_series() {
    let gas = gas(1.0, 101_325.0);
    let duration = 0.6 * diffusion_time(gas);
    let pressure = 10.0;
    let tau = gas.dynamic_viscosity * duration / (gas.density * RADIUS.powi(2));
    // Positive J0 zeros; the omitted tail is <1e-32 at this positive time.
    // Integrating the independent radial eigenfunction solution yields
    // U(t)/U∞ = 1 - 32 Σ exp(-j0n² νt/R²)/j0n⁴.
    let roots: [f64; 8] = [2.404_825_557_695_773, 5.520_078_110_286_311,
        8.653_727_912_911_013, 11.791_534_439_014_28, 14.930_917_708_487_787,
        18.071_063_967_910_922, 21.211_636_629_879_26, 24.352_471_530_749_302];
    let tail = roots.iter().map(|x| (-x * x * tau).exp() / x.powi(4)).sum::<f64>();
    let exact = pressure * RADIUS.powi(2) / (8.0 * gas.dynamic_viscosity * LENGTH) * (1.0 - 32.0 * tail);
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut channel = channel(gas, 512, duration / f64::from(steps));
        for _ in 0..steps { channel.step(pressure).unwrap(); }
        errors.push((channel.last_step().bulk_speed_m_s - exact).abs() / exact);
    }
    assert!(errors[0] > 0.01, "time axis must be live: {errors:?}");
    for e in errors.windows(2) {
        assert!((0.9..1.1).contains(&(e[0] / e[1]).log2()), "time errors {errors:?}");
    }
}

#[test]
fn circular_pressure_reversal_closes_separate_heat_and_numerical_work_and_refuses_atomically() {
    let gas = gas(1.0, 101_325.0);
    let mut channel = channel(gas, 32, diffusion_time(gas) / 40.0);
    let mut returned_work = false;
    for pressure in [10.0; 40].into_iter().chain([-1.0; 120]) {
        let frame = channel.step(pressure).unwrap();
        let scale = frame.pressure_work_j.abs() + frame.stored_delta_j.abs()
            + frame.viscous_dissipation_j + frame.numerical_dissipation_j;
        assert!(frame.energy_residual_j.abs() < 1e-12 * scale);
        assert!(frame.viscous_dissipation_j >= 0.0 && frame.numerical_dissipation_j >= 0.0);
        returned_work |= frame.pressure_work_j < 0.0;
    }
    assert!(returned_work, "old momentum must return work against reversed pressure");
    assert!(channel.last_step().flow_m3_s < 0.0);
    let before = channel.velocities_m_s().to_vec();
    let last = channel.last_step();
    let load = channel.pressure_load().unwrap();
    let fixed: FixedChannel = channel.clone().into();
    let candidate = fixed.preview_at_flow(1e-7).unwrap();
    let actual = candidate.last_step().flow_m3_s;
    assert!((actual - 1e-7).abs() < 1e-20);
    let dt = channel.spec().sample_period_s;
    assert!((candidate.last_step().pressure_work_j - load.pressure_drop_pa(1e-7).unwrap() * actual * dt).abs() < 1e-24);
    assert_eq!(channel.velocities_m_s(), before);
    assert_eq!(channel.last_step(), last);
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(channel.step(bad).is_err());
        assert_eq!(channel.velocities_m_s(), before);
        assert_eq!(channel.last_step(), last);
    }
}

#[test]
fn circular_material_density_changes_startup_and_viscosity_changes_steady_drag() {
    let ordinary = gas(1.0, 101_325.0);
    let dense = gas(1.0, 202_650.0);
    let viscous = gas(2.0, 101_325.0);
    let dt = diffusion_time(ordinary) / 1000.0;
    let mut a = channel(ordinary, 64, dt);
    let mut b = channel(dense, 64, dt);
    let c = channel(viscous, 64, dt);
    let ratio = a.step(10.0).unwrap().bulk_speed_m_s / b.step(10.0).unwrap().bulk_speed_m_s;
    assert!((1.98..2.0).contains(&ratio));
    let steady = a.steady_bulk_speed_m_s(10.0).unwrap();
    assert!((steady / b.steady_bulk_speed_m_s(10.0).unwrap() - 1.0).abs() < 1e-13);
    assert!((steady / c.steady_bulk_speed_m_s(10.0).unwrap() - 2.0).abs() < 1e-13);
}

#[test]
fn circular_geometry_and_incoherent_material_refuse_before_motion() {
    let gas = gas(1.0, 101_325.0);
    let base = CircularChannelSpec { length_m: LENGTH, radius_m: RADIUS,
        radial_cells: 16, sample_period_s: 1e-5, exit: SlitExit::StaticPressure };
    for spec in [CircularChannelSpec { radius_m: 0.0, ..base },
        CircularChannelSpec { radius_m: f64::NAN, ..base },
        CircularChannelSpec { radial_cells: 1, ..base },
        CircularChannelSpec { radial_cells: 4097, ..base },
        CircularChannelSpec { length_m: -1.0, ..base }] {
        assert!(CircularChannel::from_gas(&gas, spec).is_err());
    }
    let incoherent = GasState { density: 2.0 * gas.density, ..gas };
    assert!(CircularChannel::from_gas(&incoherent, base).is_err());
}
