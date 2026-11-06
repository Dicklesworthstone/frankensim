use super::*;

#[test]
fn couette_reports_signed_wall_traction_and_constant_cross_gap_stress() {
    let model = ShearedFilm::new(film(&[H, H]), vec![walls(1.0), walls(-1.0)]).unwrap();
    for (lower, upper) in [(0.002, 0.006), (0.006, 0.002), (0.004, 0.004)] {
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [ShearChannelReport::default(); 2];
        model
            .evaluate_with_channels_into(
                &[0.0; 3],
                &[0.0, lower, upper],
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        assert_eq!(pressure, [0.0]);
        for (channel, sign) in channels.into_iter().zip([1.0, -1.0]) {
            let expected = sign * MU * (upper - lower) / H;
            near(channel.lower_wall_traction_pa, -expected, 2e-15);
            near(channel.upper_wall_traction_pa, expected, 2e-15);
            for z in [0.0, 0.25, 0.5, 0.75, 1.0] {
                near(channel.shear_stress_pa(z).unwrap(), expected, 2e-15);
            }
            let work = W
                * L
                * sign
                * (channel.lower_wall_traction_pa * lower + channel.upper_wall_traction_pa * upper);
            near(channel.wall_power_w, work, 1e-19);
            near(channel.wall_power_w, channel.dissipated_power_w(), 1e-19);
        }
    }
}

#[test]
fn pressure_flow_retains_wall_shear_when_all_wall_rows_are_zero() {
    let model = PressureDrivenShearedFilm::new(
        ShearedFilm::new(film(&[H, H]), vec![walls(0.0); 2]).unwrap(),
        &[Some(1500.0), Some(-500.0)],
    )
    .unwrap();
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default(); 2];
    let report = model
        .evaluate_into(
            &[0.0; 3],
            &[0.0; 3],
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    for (channel, dp) in channels.into_iter().zip([-1000.0, 1000.0]) {
        let traction = -H * dp / (2.0 * L);
        near(channel.lower_wall_traction_pa, traction, 2e-14);
        near(channel.upper_wall_traction_pa, traction, 2e-14);
        near(channel.shear_stress_pa(0.5).unwrap(), 0.0, 2e-14);
        near(
            channel.shear_stress_pa(0.25).unwrap(),
            H * dp / (4.0 * L),
            2e-14,
        );
        near(
            W * L * (channel.lower_wall_traction_pa + channel.upper_wall_traction_pa),
            -W * H * dp,
            1e-17,
        );
        assert_eq!(channel.wall_power_w, 0.0);
        assert_eq!(channel.shear_dissipated_power_w, 0.0);
        assert!(channel.pressure_dissipated_power_w > 0.0);
    }
    near(report.mechanical_power_w, 0.0, 1e-17);
}

#[test]
fn stratified_stress_uses_fluidity_centroid_and_reproduces_wall_loads() {
    let model = pressure::supplied(0.0, true);
    let h = H * 0.9;
    let velocity = [0.0, 0.02, 0.06];
    let (_, force, pressure, channels) = pressure::response(&model, [H - h, 0.0, 0.0], velocity);
    // Integrate 1/mu and y/mu independently in physical gap coordinates.
    let fluidity = (h / 2.0) / 0.01 + (h / 2.0) / 0.04;
    let moment = (h / 2.0).powi(2) / (2.0 * 0.01) + (h.powi(2) - (h / 2.0).powi(2)) / (2.0 * 0.04);
    let mut projected = [0.0; 2];
    let mut wall_power = 0.0;
    for (channel, sign, boundary) in [(0, -1.0, 1500.0), (1, 1.0, -500.0)] {
        let dp = pressure - boundary;
        let constant = (sign * (velocity[2] - velocity[1]) + dp * moment / L) / fluidity;
        let output = channels[channel];
        near(output.lower_wall_traction_pa, -constant, 2e-14);
        near(output.upper_wall_traction_pa, constant - dp * h / L, 2e-14);
        for z in [0.0, 0.25, 0.5, 0.75, 1.0] {
            near(
                output.shear_stress_pa(z).unwrap(),
                constant - dp * h * z / L,
                2e-14,
            );
        }
        projected[0] += W * L * sign * output.lower_wall_traction_pa;
        projected[1] += W * L * sign * output.upper_wall_traction_pa;
        wall_power += output.wall_power_w;
    }
    near(projected[0], force[1], 1e-17);
    near(projected[1], force[2], 1e-17);
    near(
        wall_power,
        force[1] * velocity[1] + force[2] * velocity[2],
        1e-18,
    );
    // This case extracts mechanical power from the pressure supply; signed
    // wall work must not be misreported as negative viscous dissipation.
    assert!(wall_power < 0.0);
}

#[test]
fn translating_and_sliding_walls_close_pressure_pumping_work() {
    for alpha in [0.0, 3e-4] {
        for layered in [false, true] {
            let model = pressure::supplied(alpha, layered);
            for velocity in [[0.0, 0.02, 0.06], [0.0, 0.04, 0.04], [0.0, -0.06, -0.02]] {
                let (_, force, pressure, channels) = pressure::response(&model, [0.0; 3], velocity);
                let mut wall_input = 0.0;
                for (index, sign, boundary) in [(0, -1.0, 1500.0), (1, 1.0, -500.0)] {
                    let channel = channels[index];
                    near(
                        channel.volume_flow_m3_s * (pressure - boundary) + channel.wall_power_w,
                        channel.dissipated_power_w(),
                        1e-17,
                    );
                    near(
                        W * L
                            * sign
                            * (channel.lower_wall_traction_pa * velocity[1]
                                + channel.upper_wall_traction_pa * velocity[2]),
                        channel.wall_power_w,
                        1e-17,
                    );
                    wall_input += channel.wall_power_w;
                }
                near(
                    wall_input,
                    force[1] * velocity[1] + force[2] * velocity[2],
                    1e-17,
                );
                if velocity[1] == velocity[2] {
                    // Translation performs reversible pressure-pumping work
                    // even though the relative-wall shear contribution is zero.
                    assert_ne!(wall_input, 0.0);
                    assert_eq!(channels[0].shear_dissipated_power_w, 0.0);
                    assert_eq!(channels[1].shear_dissipated_power_w, 0.0);
                }
            }
        }
    }
}

#[test]
fn closed_passage_has_no_wall_stress_and_observation_refuses_extrapolation() {
    let model = ShearedFilm::new(film(&[H, 0.0]), vec![walls(0.0); 2]).unwrap();
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default(); 2];
    model
        .evaluate_with_channels_into(
            &[0.0; 3],
            &[0.001, 0.0, 0.0],
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert!(pressure[0] > 0.0);
    assert_eq!(channels[1].lower_wall_traction_pa, 0.0);
    assert_eq!(channels[1].upper_wall_traction_pa, 0.0);
    assert_eq!(channels[1].wall_power_w, 0.0);
    for z in [-0.1, 1.1, f64::INFINITY, f64::NAN] {
        assert!(channels[0].shear_stress_pa(z).is_err());
    }
}
