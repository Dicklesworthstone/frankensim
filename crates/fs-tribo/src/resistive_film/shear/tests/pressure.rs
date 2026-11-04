use super::*;
use super::super::pressure::PressureDrivenShearedFilm;

fn supplied(alpha: f64, layered: bool) -> PressureDrivenShearedFilm {
    let native = film(&[H, H]).with_barus_viscosity(alpha).unwrap();
    let wall_rows = vec![walls(-1.0), walls(1.0)];
    let moving = if layered {
        ShearedFilm::from_stratified(
            StratifiedFilm::new(&native, vec![layers(), layers()]).unwrap(),
            wall_rows,
        )
        .unwrap()
    } else {
        ShearedFilm::new(native, wall_rows).unwrap()
    };
    PressureDrivenShearedFilm::new(moving, &[Some(1500.0), Some(-500.0)]).unwrap()
}

fn response(
    model: &PressureDrivenShearedFilm,
    q: [f64; 3],
    velocity: [f64; 3],
) -> (PressureDrivenReport, [f64; 3], f64, [ShearChannelReport; 2]) {
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default(); 2];
    let report = model
        .evaluate_into(&q, &velocity, &mut force, &mut pressure, &mut channels)
        .unwrap();
    (report, force, pressure[0], channels)
}

#[test]
fn g1_pressure_fed_couette_poiseuille_closes_flow_wall_force_and_heat() {
    let model = supplied(0.0, false);
    let (report, force, pressure, channels) =
        response(&model, [0.0; 3], [0.0, 0.02, 0.06]);
    let pressure_flow = W * H.powi(3) * 1000.0 / (12.0 * MU * L);
    let entrainment = W * H * 0.04;
    let drag = MU * W * L / H;
    near(pressure, 500.0, 1e-10);
    near(force[0], AREA * 500.0, 1e-13);
    near(force[1], -W * H * 1000.0 - 2.0 * drag * 0.04, 1e-14);
    near(force[2], -W * H * 1000.0 + 2.0 * drag * 0.04, 1e-14);
    near(channels[0].volume_flow_m3_s, -pressure_flow - entrainment, 1e-20);
    near(channels[1].volume_flow_m3_s, pressure_flow + entrainment, 1e-20);
    near(report.reservoir_power_w, 2000.0 * (pressure_flow + entrainment), 1e-16);
    near(
        report.film.dissipated_power_w,
        2000.0 * pressure_flow + 2.0 * drag * 0.04_f64.powi(2),
        1e-16,
    );
    assert!(report.mechanical_power_w < 0.0);
    near(report.power_residual_w, 0.0, 1e-16);

    // Reversed translation pumps against the imposed pressure gradient.
    let (pumping, _, _, pumped_channels) =
        response(&model, [0.0; 3], [0.0, -0.06, -0.02]);
    assert!(pumped_channels[1].volume_flow_m3_s < 0.0);
    assert!(pumping.reservoir_power_w < 0.0);
    assert!(pumping.mechanical_power_w > pumping.film.dissipated_power_w);
    near(pumping.power_residual_w, 0.0, 1e-16);
}

#[test]
fn g1_equal_pressurized_supplies_retain_barus_shear_drag_and_static_load() {
    let alpha = 3e-4;
    let model = supplied(alpha, false)
        .with_reservoir_pressures(&[Some(2000.0), Some(2000.0)])
        .unwrap();
    let (report, force, pressure, channels) =
        response(&model, [0.0; 3], [0.0, 0.02, 0.06]);
    let drag = MU * (alpha * 2000.0).exp() * W * L / H;
    assert_eq!(pressure, 2000.0);
    near(force[0], AREA * 2000.0, 1e-13);
    near(force[1], -2.0 * drag * 0.04, 1e-14);
    near(force[2], 2.0 * drag * 0.04, 1e-14);
    assert_eq!(report.reservoir_power_w, 0.0);
    for channel in channels {
        assert_eq!(channel.pressure_dissipated_power_w, 0.0);
        near(channel.shear_dissipated_power_w, drag * 0.04_f64.powi(2), 1e-16);
    }
    near(report.power_residual_w, 0.0, 1e-16);
}

#[test]
fn g1_supplied_layer_work_matches_independent_piecewise_stokes_integrals() {
    let model = supplied(0.0, true);
    let velocity = [0.0, 0.02, 0.06];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default(); 2];
    let mut heat = [0.0; 4];
    let report = model
        .evaluate_with_layers_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
            &mut heat,
        )
        .unwrap();
    let fluidity = (H / 2.0) / 0.01 + (H / 2.0) / 0.04;
    let first_moment =
        (H / 2.0).powi(2) / (2.0 * 0.01) + (H.powi(2) - (H / 2.0).powi(2)) / (2.0 * 0.04);
    for (channel, sign, boundary) in [(0, -1.0, 1500.0), (1, 1.0, -500.0)] {
        let dp = pressure[0] - boundary;
        let stress_constant =
            (sign * (velocity[2] - velocity[1]) + dp / L * first_moment) / fluidity;
        for (layer, (lower, upper, mu)) in
            [(0.0, H / 2.0, 0.01), (H / 2.0, H, 0.04)].into_iter().enumerate()
        {
            let width = upper - lower;
            let midpoint_stress = stress_constant - dp / L * (lower + upper) / 2.0;
            let expected = W * L / mu
                * (midpoint_stress.powi(2) * width + (dp / L).powi(2) * width.powi(3) / 12.0);
            near(heat[2 * channel + layer], expected, 2e-15);
        }
        near(
            heat[2 * channel] + heat[2 * channel + 1],
            channels[channel].dissipated_power_w(),
            2e-15,
        );
    }
    near(heat.iter().sum(), report.film.dissipated_power_w, 2e-15);
    near(report.power_residual_w, 0.0, 2e-15);
}

#[test]
fn g1_complete_supplied_tangent_differentiates_gap_walls_and_both_pressures() {
    for alpha in [0.0, 3e-4] {
        for layered in [false, true] {
            let model = supplied(alpha, layered);
            let q = [H / 20.0, 0.0, 0.0];
            let velocity = [0.00003, 0.002, -0.003];
            let dq = [H / 100.0, 0.0, 0.0];
            let dv = [-0.000005, 0.003, 0.002];
            let dp = [80.0, -120.0];
            let mut tangent = [0.0; 3];
            model
                .tangent_into(&q, &velocity, &dq, &dv, &dp, &mut tangent)
                .unwrap();
            let epsilon = 1e-4;
            let mut perturbed = [[0.0; 3]; 2];
            for (side, sign) in [-1.0_f64, 1.0].into_iter().enumerate() {
                let delta = sign * epsilon;
                let changed = model
                    .with_reservoir_pressures(&[
                        Some(1500.0 + delta * dp[0]),
                        Some(-500.0 + delta * dp[1]),
                    ])
                    .unwrap();
                let changed_q = std::array::from_fn(|j| q[j] + delta * dq[j]);
                let changed_velocity = std::array::from_fn(|j| velocity[j] + delta * dv[j]);
                perturbed[side] = response(&changed, changed_q, changed_velocity).1;
            }
            for j in 0..3 {
                let finite_difference = (perturbed[1][j] - perturbed[0][j]) / (2.0 * epsilon);
                near(finite_difference, tangent[j], 2e-8 * tangent[j].abs().max(1e-5));
            }
        }
    }
}

#[test]
fn g3_zero_supplies_recover_native_sliding_and_refresh_retains_supplies() {
    for alpha in [0.0, 3e-4] {
        let model = supplied(alpha, true);
        let zero = model
            .with_reservoir_pressures(&[Some(0.0), Some(0.0)])
            .unwrap();
        let q = [0.0; 3];
        let velocity = [0.00002, 0.001, -0.002];
        let (actual, force, pressure, channels) = response(&zero, q, velocity);
        let mut expected_force = [0.0; 3];
        let mut expected_pressure = [0.0];
        let mut expected_channels = [ShearChannelReport::default(); 2];
        let expected = model
            .sheared_film()
            .evaluate_with_channels_into(
                &q,
                &velocity,
                &mut expected_force,
                &mut expected_pressure,
                &mut expected_channels,
            )
            .unwrap();
        for j in 0..3 {
            near(force[j], expected_force[j], 1e-13);
        }
        near(pressure, expected_pressure[0], 1e-10);
        near(actual.film.dissipated_power_w, expected.dissipated_power_w, 1e-16);
        for j in 0..2 {
            near(channels[j].volume_flow_m3_s, expected_channels[j].volume_flow_m3_s, 1e-20);
        }

        let viscosities = [0.02, 0.08, 0.02, 0.08];
        let refreshed = model.with_layer_viscosities(&viscosities).unwrap();
        assert_eq!(refreshed.reservoir_pressures_pa(), model.reservoir_pressures_pa());
        let independent = PressureDrivenShearedFilm::new(
            model.sheared_film().with_layer_viscosities(&viscosities).unwrap(),
            model.reservoir_pressures_pa(),
        )
        .unwrap();
        assert_eq!(
            response(&refreshed, q, velocity).1,
            response(&independent, q, velocity).1,
        );
    }
}

#[test]
fn g0_bad_supplies_and_refused_solves_publish_no_partial_outputs() {
    let model = supplied(3e-4, true);
    for invalid in [
        [None, Some(0.0)],
        [Some(f64::NAN), Some(0.0)],
        [Some(20000.0), Some(0.0)],
    ] {
        assert!(model.with_reservoir_pressures(&invalid).is_err());
    }
    let mut force = [71.0; 3];
    let mut pressure = [73.0];
    let mut channels = [ShearChannelReport::default(); 2];
    let original_channels = channels;
    let mut heat = [79.0; 4];
    assert!(model
        .evaluate_with_layers_into(
            &[0.0; 3],
            &[100.0, 0.02, 0.06],
            &mut force,
            &mut pressure,
            &mut channels,
            &mut heat,
        )
        .is_err());
    assert_eq!(force, [71.0; 3]);
    assert_eq!(pressure, [73.0]);
    assert_eq!(channels, original_channels);
    assert_eq!(heat, [79.0; 4]);
    assert!(model
        .tangent_into(
            &[0.0; 3],
            &[0.0; 3],
            &[0.0; 3],
            &[0.0; 3],
            &[f64::NAN, 0.0],
            &mut force,
        )
        .is_err());
    assert_eq!(force, [71.0; 3]);
    assert!(response(&model, [0.0; 3], [0.0; 3]).0.film.dissipated_power_w > 0.0);
}
