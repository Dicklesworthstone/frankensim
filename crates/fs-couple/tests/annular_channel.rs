//! G1/G3 annular momentum: independent logarithmic Couette and annular
//! Poiseuille solutions, narrow-gap limit, material history and atomic refusal.
//! Axial wall motion only; these tests do not qualify rotating-cylinder flow.

use fs_couple::bernoulli_aperture::slit_channel::{
    SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
    annular::{AnnularChannel, AnnularChannelSpec},
};
use fs_material::gas::{GasSpec, GasState};

const INNER: f64 = 0.0002;
const OUTER: f64 = 0.0005;
const LENGTH: f64 = 0.02;

fn gas(viscosity_scale: f64, pressure: f64) -> GasState {
    let mut spec = GasSpec::dry_air_ussa1976();
    // Explicit synthetic source intervention, not a measured species.
    spec.sutherland_beta *= viscosity_scale;
    GasState::try_new(&spec, 293.15, pressure).unwrap()
}

fn diffusion_time(gas: GasState, gap: f64) -> f64 {
    gap.powi(2) * gas.density / (gas.dynamic_viscosity * core::f64::consts::PI.powi(2))
}

fn spec(gas: GasState, cells: usize) -> AnnularChannelSpec {
    AnnularChannelSpec {
        length_m: LENGTH,
        inner_radius_m: INNER,
        outer_radius_m: OUTER,
        radial_cells: cells,
        sample_period_s: diffusion_time(gas, OUTER - INNER) / 10.0,
        exit: SlitExit::StaticPressure,
    }
}

fn closes_energy(step: SlitChannelStep) {
    let scale = step.pressure_work_j.abs()
        + step.wall_work_j.abs()
        + step.stored_delta_j.abs()
        + step.viscous_dissipation_j
        + step.exit_work_j
        + step.numerical_dissipation_j;
    assert!(step.energy_residual_j.abs() <= 2e-12 * scale, "{step:?}");
    assert!(step.viscous_dissipation_j >= 0.0 && step.numerical_dissipation_j >= 0.0);
}

#[test]
fn annular_poiseuille_profile_and_flow_have_second_order_spatial_convergence() {
    let gas = gas(1.0, 101_325.0);
    let pressure = 10.0;
    let log_ratio = (OUTER / INNER).ln();
    let squared_gap = (OUTER - INNER) * (OUTER + INNER);
    let exact_mean = pressure / (8.0 * gas.dynamic_viscosity * LENGTH)
        * (INNER.powi(2) + OUTER.powi(2) - squared_gap / log_ratio);
    let mut flow_errors = Vec::new();
    let mut profile_errors = Vec::new();
    for cells in [8, 16, 32] {
        let mut channel = AnnularChannel::from_gas(&gas, spec(gas, cells)).unwrap();
        for _ in 0..500 {
            closes_energy(channel.step(pressure).unwrap());
        }
        let step = channel.last_step();
        let steady = channel.steady_bulk_speed_m_s(pressure).unwrap();
        assert!((step.bulk_speed_m_s - steady).abs() < 1e-10 * exact_mean);
        flow_errors.push((step.bulk_speed_m_s / exact_mean - 1.0).abs());
        let dr = (OUTER - INNER) / cells as f64;
        let mut maximum = 0.0_f64;
        for (i, actual) in channel.velocities_m_s().iter().enumerate() {
            let radius = INNER + (i as f64 + 0.5) * dr;
            let exact = pressure / (4.0 * gas.dynamic_viscosity * LENGTH)
                * (squared_gap * (radius / INNER).ln() / log_ratio
                    - (radius - INNER) * (radius + INNER));
            maximum = maximum.max((actual - exact).abs());
        }
        profile_errors.push(maximum);
        let area = core::f64::consts::PI * squared_gap;
        assert!((step.flow_m3_s - area * step.bulk_speed_m_s).abs() < 2e-14 * step.flow_m3_s.abs());
    }
    assert!(
        flow_errors[0] > 0.02,
        "the radial refinement axis must be live"
    );
    for errors in [&flow_errors, &profile_errors] {
        for pair in errors.windows(2) {
            assert!(
                (1.95..2.05).contains(&(pair[0] / pair[1]).log2()),
                "{errors:?}"
            );
        }
    }
}

#[test]
fn annular_axial_couette_has_logarithmic_profile_and_physical_wall_work() {
    let gas = gas(1.0, 101_325.0);
    let cells = 17;
    let mut channel = AnnularChannel::from_gas(&gas, spec(gas, cells)).unwrap();
    let walls = [-0.2, 0.8];
    for _ in 0..500 {
        closes_energy(channel.step_with_wall_velocities(0.0, walls).unwrap());
    }
    let log_ratio = (OUTER / INNER).ln();
    let dr = (OUTER - INNER) / cells as f64;
    for (i, actual) in channel.velocities_m_s().iter().enumerate() {
        let radius = INNER + (i as f64 + 0.5) * dr;
        let exact = walls[0] + (walls[1] - walls[0]) * (radius / INNER).ln() / log_ratio;
        assert!(
            (actual - exact).abs() < 2e-12,
            "cell {i}: {actual} vs {exact}"
        );
    }
    let force =
        core::f64::consts::TAU * gas.dynamic_viscosity * LENGTH * (walls[1] - walls[0]) / log_ratio;
    let step = channel.last_step();
    let dt = channel.spec().sample_period_s;
    assert!((step.wall_impulse_n_s[0] / dt + force).abs() < 2e-12 * force);
    assert!((step.wall_impulse_n_s[1] / dt - force).abs() < 2e-12 * force);
    let heat = force * (walls[1] - walls[0]) * dt;
    assert!((step.viscous_dissipation_j / heat - 1.0).abs() < 2e-12);
    assert!((step.wall_work_j / heat - 1.0).abs() < 2e-12);
    assert_eq!(step.pressure_work_j, 0.0);
}

#[test]
fn tiny_inner_cylinder_retains_its_logarithmic_drag_instead_of_becoming_an_axis() {
    let gas = gas(1.0, 101_325.0);
    let inner = 1e-12;
    let cells = 32;
    let mut shape = spec(gas, cells);
    shape.inner_radius_m = inner;
    shape.sample_period_s = diffusion_time(gas, OUTER - inner) / 10.0;
    let mut channel = AnnularChannel::from_gas(&gas, shape).unwrap();
    for _ in 0..800 {
        closes_energy(channel.step_with_wall_velocities(0.0, [1.0, 0.0]).unwrap());
    }
    let log_ratio = (OUTER / inner).ln();
    let dr = (OUTER - inner) / cells as f64;
    for (i, actual) in channel.velocities_m_s().iter().enumerate() {
        let radius = inner + (i as f64 + 0.5) * dr;
        let exact = (OUTER / radius).ln() / log_ratio;
        assert!((actual - exact).abs() < 2e-12);
    }
    let force = core::f64::consts::TAU * gas.dynamic_viscosity * LENGTH / log_ratio;
    let step = channel.last_step();
    assert!((step.wall_impulse_n_s[0] / shape.sample_period_s / force - 1.0).abs() < 2e-12);
    assert!((step.wall_impulse_n_s[1] / shape.sample_period_s / force + 1.0).abs() < 2e-12);
    assert!(
        channel
            .continuous_operator()
            .unwrap()
            .face_conductances_kg_s()[0]
            > 0.0
    );
}

#[test]
fn narrow_annulus_preserves_area_and_converges_to_the_planar_gap() {
    let gas = gas(1.0, 101_325.0);
    let inner = 1.0;
    let outer = inner + 1e-9;
    let gap = outer - inner;
    let width = core::f64::consts::TAU * (inner + 0.5 * gap);
    let dt = diffusion_time(gas, gap) / 20.0;
    let mut annular = AnnularChannel::from_gas(
        &gas,
        AnnularChannelSpec {
            length_m: LENGTH,
            inner_radius_m: inner,
            outer_radius_m: outer,
            radial_cells: 32,
            sample_period_s: dt,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap();
    let mut planar = SlitChannel::from_gas(
        &gas,
        SlitChannelSpec {
            length_m: LENGTH,
            width_m: width,
            height_m: gap,
            transverse_cells: 32,
            sample_period_s: dt,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap();
    let operator = annular.continuous_operator().unwrap();
    let area: f64 = operator.cell_areas_m2().iter().sum();
    assert!((area / (width * gap) - 1.0).abs() < 2e-14);
    assert!(
        operator
            .face_conductances_kg_s()
            .iter()
            .all(|g| g.is_finite() && *g > 0.0)
    );
    for _ in 0..60 {
        closes_energy(
            annular
                .step_with_wall_velocities(10.0, [-0.2, 0.8])
                .unwrap(),
        );
        closes_energy(planar.step_with_wall_velocities(10.0, [-0.2, 0.8]).unwrap());
    }
    for (round, flat) in annular.velocities_m_s().iter().zip(planar.velocities_m_s()) {
        assert!((round - flat).abs() < 3e-10);
    }
    assert!(
        (annular.last_step().flow_m3_s - planar.last_step().flow_m3_s).abs() < 3e-10 * width * gap
    );
}

#[test]
fn annular_pressure_port_wall_impulse_and_refusal_share_the_original_state() {
    let gas = gas(1.0, 101_325.0);
    let mut channel = AnnularChannel::from_gas(&gas, spec(gas, 24)).unwrap();
    for _ in 0..30 {
        channel.step(3.0).unwrap();
    }
    let old = channel.velocities_m_s().to_vec();
    let last = channel.last_step();
    let walls = [0.3, -0.1];
    let flow = -1e-8;
    let load = channel.pressure_load_with_wall_velocities(walls).unwrap();
    let pressure = load.pressure_drop_pa(flow).unwrap();
    let trial = channel
        .preview_at_flow_with_wall_velocities(flow, walls)
        .unwrap();
    let step = trial.last_step();
    closes_energy(step);
    assert!((step.flow_m3_s - flow).abs() < 1e-20);
    let operator = channel.continuous_operator().unwrap();
    let change: f64 = operator
        .cell_masses_kg()
        .iter()
        .zip(old.iter().zip(trial.velocities_m_s()))
        .map(|(mass, (before, after))| mass * (after - before))
        .sum();
    let pressure_impulse =
        pressure * operator.cell_areas_m2().iter().sum::<f64>() * channel.spec().sample_period_s;
    let wall_impulse: f64 = step.wall_impulse_n_s.iter().sum();
    assert!(
        (change - pressure_impulse - wall_impulse).abs()
            < 2e-12 * (change.abs() + pressure_impulse.abs() + wall_impulse.abs())
    );
    assert_eq!(channel.velocities_m_s(), old);
    assert_eq!(channel.last_step(), last);
    for walls in [[f64::NAN, 0.0], [0.0, f64::INFINITY], [f64::MAX, f64::MAX]] {
        assert!(channel.step_with_wall_velocities(0.0, walls).is_err());
        assert_eq!(channel.velocities_m_s(), old);
        assert_eq!(channel.last_step(), last);
    }
    assert!(channel.step(f64::NAN).is_err());
    let mut same = channel.clone();
    assert_eq!(
        channel.step_with_wall_velocities(-1.0, [0.0; 2]).unwrap(),
        same.step(-1.0).unwrap()
    );

    let mut shape = spec(gas, 24);
    shape.exit = SlitExit::UniformKineticHead;
    let mut terminal = AnnularChannel::from_gas(&gas, shape).unwrap();
    for speed in [0.8; 50].into_iter().chain([-0.8; 120]) {
        closes_energy(
            terminal
                .step_with_wall_velocities(0.0, [speed, 0.0])
                .unwrap(),
        );
    }
    assert!(terminal.last_step().flow_m3_s < 0.0);
    assert!(terminal.last_step().exit_work_j > 0.0);
    assert!(terminal.pressure_load_with_wall_velocities(walls).is_err());
    assert!(terminal.continuous_operator().is_err());
}

#[test]
fn annular_source_changes_preserve_momentum_and_refuse_atomically() {
    let ordinary = gas(1.0, 101_325.0);
    let dense = gas(1.0, 202_650.0);
    let viscous_dense = gas(2.0, 202_650.0);
    let mut channel = AnnularChannel::from_gas(&ordinary, spec(ordinary, 32)).unwrap();
    for _ in 0..40 {
        channel.step(4.0).unwrap();
    }
    let before = channel.velocities_m_s().to_vec();
    let last = channel.last_step();
    let energy = channel.stored_energy_j();
    let update = channel
        .rebind_gas_controlled(&dense, energy, &mut || false)
        .unwrap();
    for (old, new) in before.iter().zip(channel.velocities_m_s()) {
        assert!(
            (ordinary.density * old - dense.density * new).abs()
                <= 2e-15 * (ordinary.density * old).abs()
        );
    }
    assert!((update.parameter_work_j / energy + 0.5).abs() < 2e-14);
    assert_eq!(channel.last_step(), last);
    let before = channel.velocities_m_s().to_vec();
    let stored = channel.stored_energy_j();
    assert!(
        channel
            .rebind_gas_controlled(&ordinary, 0.0, &mut || false)
            .is_err()
    );
    assert!(
        channel
            .rebind_gas_controlled(&ordinary, energy, &mut || true)
            .is_err()
    );
    assert_eq!(channel.velocities_m_s(), before);
    assert_eq!(channel.stored_energy_j(), stored);
    assert_eq!(channel.gas(), &dense);
    let old_impedance = channel.pressure_load().unwrap().impedance_pa_s_m3;
    let update = channel
        .rebind_gas_controlled(&viscous_dense, 0.0, &mut || false)
        .unwrap();
    assert_eq!(update.parameter_work_j, 0.0);
    assert_eq!(channel.velocities_m_s(), before);
    assert!(channel.pressure_load().unwrap().impedance_pa_s_m3 > old_impedance);
}

#[test]
fn annular_invalid_geometry_and_incoherent_material_refuse() {
    let gas = gas(1.0, 101_325.0);
    let shape = spec(gas, 16);
    for bad in [
        AnnularChannelSpec {
            inner_radius_m: 0.0,
            ..shape
        },
        AnnularChannelSpec {
            inner_radius_m: OUTER,
            ..shape
        },
        AnnularChannelSpec {
            outer_radius_m: INNER / 2.0,
            ..shape
        },
        AnnularChannelSpec {
            outer_radius_m: f64::INFINITY,
            ..shape
        },
        AnnularChannelSpec {
            radial_cells: 1,
            ..shape
        },
        AnnularChannelSpec {
            radial_cells: 4097,
            ..shape
        },
        AnnularChannelSpec {
            sample_period_s: 0.0,
            ..shape
        },
    ] {
        assert!(AnnularChannel::from_gas(&gas, bad).is_err());
    }
    let incoherent = GasState {
        density: 2.0 * gas.density,
        ..gas
    };
    assert!(AnnularChannel::from_gas(&incoherent, shape).is_err());
}
