//! G1/G3 resolved harmonic pressure and moving-wall ports under exp(-i*omega*t).
//! The oscillatory planar Couette reference follows directly from mu*u''=-i*omega*rho*u:
//! u(y)=U*sinh(k*y)/sinh(k*H), k=(1-i)*sqrt(omega*rho/(2*mu)).

use fs_couple::bernoulli_aperture::slit_channel::{
    SlitChannel, SlitChannelSpec, SlitExit,
    annular::{AnnularChannel, AnnularChannelSpec},
    circular::{CircularChannel, CircularChannelSpec},
    harmonic::{
        ChannelContinuousOperator,
        walls::{ChannelWallHarmonicResponse, WALL_HARMONIC_WORK_PER_CELL},
    },
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;

const LENGTH: f64 = 0.02;
const WIDTH: f64 = 0.01;
const HEIGHT: f64 = 0.001;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn slit(cells: usize, dt: f64) -> SlitChannel {
    SlitChannel::from_gas(
        &gas(),
        SlitChannelSpec {
            length_m: LENGTH,
            width_m: WIDTH,
            height_m: HEIGHT,
            transverse_cells: cells,
            sample_period_s: dt,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap()
}
fn solve(
    operator: ChannelContinuousOperator<'_>,
    omega: f64,
    pressure: C64,
    walls: [C64; 2],
) -> ChannelWallHarmonicResponse {
    operator
        .solve_harmonic_with_walls(
            omega,
            pressure,
            walls,
            WALL_HARMONIC_WORK_PER_CELL * operator.dimension(),
            1e-11,
            &mut || false,
        )
        .unwrap()
}
fn closes_work(response: &ChannelWallHarmonicResponse, omega: f64) {
    let scale = response.viscous_dissipation_w
        + response.pressure_complex_power_w.abs()
        + response.wall_complex_power_w[0].abs()
        + response.wall_complex_power_w[1].abs()
        + 2.0 * omega * response.mean_stored_energy_j;
    assert!(
        response.power_residual_w.abs() < 2e-11 * scale,
        "{response:?}"
    );
    assert!(
        response.reactive_power_residual_w.abs() < 2e-11 * scale,
        "{response:?}"
    );
    assert!(response.viscous_dissipation_w >= 0.0 && response.mean_stored_energy_j >= 0.0);
}
fn sinh(z: C64) -> C64 {
    C64::new(z.re.sinh() * z.im.cos(), z.re.cosh() * z.im.sin())
}
fn cosh(z: C64) -> C64 {
    C64::new(z.re.cosh() * z.im.cos(), z.re.sinh() * z.im.sin())
}

#[test]
fn oscillatory_couette_profile_flow_and_wall_force_converge_to_the_continuum() {
    let gas = gas();
    let omega = 18.0 * gas.dynamic_viscosity / (gas.density * HEIGHT.powi(2));
    let k = C64::new(1.0, -1.0).scale((omega * gas.density / (2.0 * gas.dynamic_viscosity)).sqrt());
    let drive = C64::new(0.7, -0.3);
    let denominator = sinh(k.scale(HEIGHT));
    let exact_flow = drive * ((cosh(k.scale(HEIGHT)) - C64::ONE) / (k * denominator)).scale(WIDTH);
    let exact_force = (drive * k * (cosh(k.scale(HEIGHT)) / denominator))
        .scale(gas.dynamic_viscosity * LENGTH * WIDTH);
    let mut profile_errors = Vec::new();
    let mut flow_errors = Vec::new();
    let mut force_errors = Vec::new();
    for cells in [32, 64, 128] {
        let channel = slit(cells, 1e-4);
        let response = solve(
            channel.continuous_operator().unwrap(),
            omega,
            C64::ZERO,
            [C64::ZERO, drive],
        );
        closes_work(&response, omega);
        let mut maximum = 0.0_f64;
        for (i, &velocity) in response.velocities_m_s.iter().enumerate() {
            let y = (i as f64 + 0.5) * HEIGHT / cells as f64;
            let exact = drive * (sinh(k.scale(y)) / denominator);
            maximum = maximum.max((velocity - exact).abs());
        }
        profile_errors.push(maximum);
        flow_errors.push((response.flow_m3_s - exact_flow).abs());
        force_errors.push((response.wall_forces_n[1] - exact_force).abs());
        // Reversing this mechanical product preserves real work but reverses
        // reactive work. The original kinetic inventory detects that error.
        let wrong_wall_power = (drive * response.wall_forces_n[1].conj()).scale(0.5);
        let reactive = 2.0 * omega * response.mean_stored_energy_j;
        assert!((wrong_wall_power.im + reactive).abs() > reactive);
    }
    assert!(
        profile_errors[0] > 1e-3,
        "the spatial convergence axis must be live"
    );
    for errors in [&profile_errors, &flow_errors, &force_errors] {
        for pair in errors.windows(2) {
            assert!(
                (1.9..2.1).contains(&(pair[0] / pair[1]).log2()),
                "{errors:?}"
            );
        }
    }
}

#[test]
fn stationary_wall_solve_matches_original_observer_and_ignores_step_clock_and_history() {
    let mut first = slit(32, 1e-6);
    let second = slit(32, 0.25);
    for _ in 0..10 {
        first.step(3.0).unwrap();
    }
    let before = first.velocities_m_s().to_vec();
    let last = first.last_step();
    let operator = first.continuous_operator().unwrap();
    let omega = 31.0;
    let pressure = C64::new(2.0, 0.3);
    let response = solve(operator, omega, pressure, [C64::ZERO; 2]);
    let original = operator
        .observe_harmonic(
            omega,
            &response.velocities_m_s,
            pressure,
            512 * operator.dimension(),
            1e-11,
            &mut || false,
        )
        .unwrap();
    assert_eq!(response.flow_m3_s, original.flow_m3_s);
    assert_eq!(response.mean_stored_energy_j, original.mean_stored_energy_j);
    assert_eq!(
        response.viscous_dissipation_w,
        original.viscous_dissipation_w
    );
    assert_eq!(response.pressure_complex_power_w.re, original.input_power_w);
    assert_eq!(response.power_residual_w, original.power_residual_w);
    assert_eq!(
        response.solve_relative_residual,
        original.solve_relative_residual
    );
    let observed = operator
        .observe_harmonic_with_walls(
            omega,
            &response.velocities_m_s,
            pressure,
            [C64::ZERO; 2],
            WALL_HARMONIC_WORK_PER_CELL * operator.dimension(),
            1e-11,
            &mut || false,
        )
        .unwrap();
    assert_eq!(response, observed);
    assert_eq!(
        response,
        solve(
            second.continuous_operator().unwrap(),
            omega,
            pressure,
            [C64::ZERO; 2]
        )
    );
    assert_eq!(first.velocities_m_s(), before);
    assert_eq!(first.last_step(), last);
}

#[test]
fn planar_circular_and_annular_mechanical_ports_close_complex_momentum_and_work() {
    let planar = slit(32, 1e-4);
    let circular = CircularChannel::from_gas(
        &gas(),
        CircularChannelSpec {
            length_m: LENGTH,
            radius_m: HEIGHT / 2.0,
            radial_cells: 32,
            sample_period_s: 1e-4,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap();
    let annular = AnnularChannel::from_gas(
        &gas(),
        AnnularChannelSpec {
            length_m: LENGTH,
            inner_radius_m: HEIGHT / 4.0,
            outer_radius_m: HEIGHT / 2.0,
            radial_cells: 32,
            sample_period_s: 1e-4,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap();
    let omega = 47.0;
    let pressure = C64::new(0.8, -0.4);
    for operator in [
        planar.continuous_operator().unwrap(),
        circular.continuous_operator().unwrap(),
        annular.continuous_operator().unwrap(),
    ] {
        let axis = operator.face_conductances_kg_s()[0] == 0.0;
        let walls = [
            if axis { C64::ZERO } else { C64::new(-0.2, 0.1) },
            C64::new(0.3, 0.7),
        ];
        let response = solve(operator, omega, pressure, walls);
        closes_work(&response, omega);
        let momentum = operator
            .cell_masses_kg()
            .iter()
            .zip(&response.velocities_m_s)
            .fold(C64::ZERO, |total, (mass, v)| total + v.scale(*mass));
        let inertia_force = C64::new(0.0, -omega) * momentum;
        let pressure_force = pressure.scale(operator.cell_areas_m2().iter().sum());
        let net = pressure_force + response.wall_forces_n[0] + response.wall_forces_n[1];
        let scale = inertia_force.abs()
            + pressure_force.abs()
            + response.wall_forces_n[0].abs()
            + response.wall_forces_n[1].abs();
        assert!((inertia_force - net).abs() < 2e-12 * scale);
        assert_eq!(
            response.wall_complex_power_w[0],
            (response.wall_forces_n[0] * walls[0].conj()).scale(0.5)
        );
        assert_eq!(
            response.wall_complex_power_w[1],
            (response.wall_forces_n[1] * walls[1].conj()).scale(0.5)
        );
        if axis {
            assert_eq!(response.wall_forces_n[0], C64::ZERO);
        }
    }
}

#[test]
fn pressure_driven_fluid_can_return_mechanical_work_to_both_moving_walls() {
    let channel = slit(64, 1e-4);
    let gas = gas();
    let omega = 0.5 * gas.dynamic_viscosity / (gas.density * HEIGHT.powi(2));
    let response = solve(
        channel.continuous_operator().unwrap(),
        omega,
        C64::from_re(20.0),
        [C64::from_re(0.1); 2],
    );
    closes_work(&response, omega);
    assert!(response.wall_complex_power_w.iter().all(|p| p.re < 0.0));
    assert!(response.viscous_dissipation_w > 0.0);
    assert!(response.pressure_complex_power_w.re > response.viscous_dissipation_w);
}

#[test]
fn almost_rigid_low_frequency_translation_retains_inventory_without_false_refusal() {
    let channel = slit(64, 1e-4);
    let gas = gas();
    let omega = 1e-14 * gas.dynamic_viscosity / (gas.density * HEIGHT.powi(2));
    let drive = C64::new(1.0, 0.3);
    let response = solve(
        channel.continuous_operator().unwrap(),
        omega,
        C64::ZERO,
        [drive; 2],
    );
    assert!(
        response
            .velocities_m_s
            .iter()
            .all(|&v| (v - drive).abs() < 2e-11)
    );
    let expected_energy = 0.25 * gas.density * LENGTH * WIDTH * HEIGHT * drive.norm_sq();
    assert!((response.mean_stored_energy_j / expected_energy - 1.0).abs() < 2e-11);
    assert!(response.viscous_dissipation_w >= 0.0 && response.viscous_dissipation_w < 1e-24);
    assert!(
        response.power_residual_w.is_finite() && response.reactive_power_residual_w.is_finite()
    );
}

#[test]
fn invalid_wall_ports_work_budgets_candidate_fields_and_cancellation_refuse() {
    let channel = slit(32, 1e-4);
    let operator = channel.continuous_operator().unwrap();
    let work = WALL_HARMONIC_WORK_PER_CELL * operator.dimension();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            operator
                .solve_harmonic_with_walls(
                    omega,
                    C64::ONE,
                    [C64::ZERO; 2],
                    work,
                    1e-11,
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(
        operator
            .solve_harmonic_with_walls(1.0, C64::ONE, [C64::ZERO; 2], work - 1, 1e-11, &mut || {
                false
            })
            .is_err()
    );
    assert!(
        operator
            .solve_harmonic_with_walls(
                1.0,
                C64::ONE,
                [C64::new(f64::NAN, 0.0); 2],
                work,
                1e-11,
                &mut || false
            )
            .is_err()
    );
    assert!(
        operator
            .solve_harmonic_with_walls(
                1.0,
                C64::new(f64::INFINITY, 0.0),
                [C64::ZERO; 2],
                work,
                1e-11,
                &mut || false
            )
            .is_err()
    );
    for tolerance in [0.0, 1.0, f64::NAN] {
        assert!(
            operator
                .solve_harmonic_with_walls(
                    1.0,
                    C64::ONE,
                    [C64::ZERO; 2],
                    work,
                    tolerance,
                    &mut || false
                )
                .is_err()
        );
    }
    let mut calls = 0;
    assert!(
        operator
            .solve_harmonic_with_walls(1.0, C64::ONE, [C64::ZERO; 2], work, 1e-11, &mut || {
                calls += 1;
                calls == 6
            })
            .is_err()
    );
    assert_eq!(calls, 6);
    assert!(
        operator
            .observe_harmonic_with_walls(
                1.0,
                &[C64::ZERO; 31],
                C64::ONE,
                [C64::ZERO; 2],
                work,
                1e-11,
                &mut || false
            )
            .is_err()
    );
    let response = solve(operator, 1.0, C64::ONE, [C64::ZERO; 2]);
    let mut wrong = response.velocities_m_s.clone();
    wrong[10] = wrong[10] + C64::ONE;
    assert!(
        operator
            .observe_harmonic_with_walls(
                1.0,
                &wrong,
                C64::ONE,
                [C64::ZERO; 2],
                work,
                1e-11,
                &mut || false
            )
            .is_err()
    );

    let circular = CircularChannel::from_gas(
        &gas(),
        CircularChannelSpec {
            length_m: LENGTH,
            radius_m: HEIGHT / 2.0,
            radial_cells: 32,
            sample_period_s: 1e-4,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap();
    assert!(
        circular
            .continuous_operator()
            .unwrap()
            .solve_harmonic_with_walls(
                1.0,
                C64::ZERO,
                [C64::ONE, C64::ZERO],
                work,
                1e-11,
                &mut || false,
            )
            .is_err()
    );
    assert!(channel.velocities_m_s().iter().all(|v| *v == 0.0));
    assert_eq!(channel.stored_energy_j(), 0.0);
}
