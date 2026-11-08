use super::*;

fn watts_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual >= 0.0 && expected >= 0.0);
    assert!(
        (actual - expected).abs() <= tolerance * expected,
        "{actual:e} != {expected:e}"
    );
}

// Independent volume quadrature of the original differential work law.
// Production uses positive closed-form wave integrals with no spatial grid.
fn volume_channels(
    wave: &crate::SegmentWave,
    length: f64,
    bounds: [f64; 2],
    forward: C64,
    backward: C64,
) -> [f64; 2] {
    let r = (wave.wavenumber * wave.characteristic_impedance).im;
    let g = (wave.wavenumber / wave.characteristic_impedance).im;
    let n = 1024;
    let dx = (bounds[1] - bounds[0]) / n as f64;
    let mut total = [0.0; 2];
    for i in 0..=n {
        let x = bounds[0] + i as f64 * dx;
        let [p, u] = analytic(
            wave.wavenumber, wave.characteristic_impedance, length, x, forward, backward,
        );
        let weight = if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
        total[0] += weight * 0.5 * r * u.norm_sq();
        total[1] += weight * 0.5 * g * p.norm_sq();
    }
    total.map(|value| value * dx / 3.0)
}

#[test]
fn coherent_cylinder_heat_matches_volume_work_and_cell_partitions() {
    let gas = gas();
    let (radius, length, omega) = (0.0012, 0.23, 2500.0);
    let wave = crate::segment_wave(&gas, radius, omega, LossModel::Bessel).unwrap();
    let forward = C64::new(1.4, -0.3);
    let backward = C64::new(-0.7, 0.6);
    let field = |x| analytic(
        wave.wavenumber, wave.characteristic_impedance, length, x, forward, backward,
    );
    // Deliberately author sources in the opposite order from canonical labels.
    let sources = [
        pressure_source("right", field(length)[0]),
        pressure_source("left", field(0.0)[0]),
    ];
    let mut net = cylinder_net(length);
    net.edges[0].duct.segments[0] = Segment::Cylinder { radius, length };
    let mut points = probes(length, &[0.0, 0.17, 0.53, 0.83, 1.0]);
    points[0].section = NetworkFieldSection::Inlet { segment: 0 };
    points[4].section = NetworkFieldSection::Outlet { segment: 0 };
    let fields = PreparedNetworkFields::new(
        &net, &gas, &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    ).unwrap();
    let heat = fields.prepare_cylinder_absorption(
        &[[0, 4], [0, 1], [1, 2], [2, 3], [3, 4], [2, 2]], 6, || false,
    ).unwrap();
    assert!(std::ptr::eq(heat.fields(), &fields));
    assert_eq!(heat.fields().source_clock(), "manufactured-common-clock");
    let values = heat.sample(fields.port_response().default_excitations(), || false).unwrap();
    for (value, cell) in values.iter().zip(heat.cells()) {
        let expected = volume_channels(&wave, length, [cell.start_m, cell.end_m], forward, backward);
        watts_close(value.viscous_power_w, expected[0], 2e-9);
        watts_close(value.thermal_power_w, expected[1], 2e-9);
        watts_close(value.absorbed_power_w, expected[0] + expected[1], 2e-9);
        assert_eq!(
            cell.volume_m3,
            core::f64::consts::PI * radius * radius * (cell.end_m - cell.start_m),
        );
    }
    let channels: [fn(&NetworkCylinderAbsorption) -> f64; 2] = [
        |v| v.viscous_power_w,
        |v| v.thermal_power_w,
    ];
    for channel in channels {
        watts_close(
            values[1..5].iter().map(channel).sum(), channel(&values[0]), 2e-10,
        );
    }
    assert_eq!(values[5], NetworkCylinderAbsorption::default());
    let p = |x| { let [p, u] = field(x); 0.5 * (p * u.conj()).re };
    watts_close(values[0].absorbed_power_w, p(0.0) - p(length), 2e-9);
}

#[test]
fn standing_wave_nodes_retain_cubic_heat_and_coherence_changes_both_channels() {
    let gas = gas();
    let (radius, length, omega) = (0.0008, 0.08, 1700.0);
    let wave = crate::segment_wave(&gas, radius, omega, LossModel::Bessel).unwrap();
    let r = (wave.wavenumber * wave.characteristic_impedance).im;
    let g = (wave.wavenumber / wave.characteristic_impedance).im;
    let forward = C64::new(1.3, -0.2);
    let at_middle = forward * phase(wave.wavenumber, length * 0.5);
    for sign in [-1.0, 1.0] {
        let backward = forward.scale(sign);
        let field = |x| analytic(
            wave.wavenumber, wave.characteristic_impedance, length, x, forward, backward,
        );
        let sources = [
            pressure_source("left", field(0.0)[0]),
            pressure_source("right", field(length)[0]),
        ];
        let mut net = cylinder_net(length);
        net.edges[0].duct.segments[0] = Segment::Cylinder { radius, length };
        for width in [1e-7, 3e-8] {
            let points: Vec<_> = [-0.5, 0.5].into_iter().map(|offset| NetworkFieldProbe {
                edge_id: "tube".into(),
                section: NetworkFieldSection::Cylinder {
                    segment: 0, offset_m: length * 0.5 + offset * width,
                },
            }).collect();
            let fields = PreparedNetworkFields::new(
                &net, &gas, &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
            ).unwrap();
            let heat = fields.prepare_cylinder_absorption(&[[0, 1]], 1, || false).unwrap();
            let drive = fields.port_response().default_excitations();
            let value = heat.sample(drive, || false).unwrap()[0];
            let cell = heat.cells()[0];
            let d = cell.end_m - cell.start_m;
            let (small, expected) = if sign < 0.0 {
                (
                    value.thermal_power_w,
                    g * (wave.wavenumber * at_middle).norm_sq() * d * d * d / 6.0,
                )
            } else {
                (
                    value.viscous_power_w,
                    r * (wave.wavenumber * at_middle / wave.characteristic_impedance).norm_sq()
                        * d * d * d / 6.0,
                )
            };
            assert!(small > 0.0);
            watts_close(small, expected, 2e-6);
            assert_eq!(heat.observe_coherence(None, &[drive], || false).unwrap()[0], value);
            let left = [drive[0], C64::ZERO];
            let right = [C64::ZERO, drive[1]];
            let independent = heat.observe_coherence(None, &[&left, &right], || false).unwrap()[0];
            let independent_channel = if sign < 0.0 {
                independent.thermal_power_w
            } else {
                independent.viscous_power_w
            };
            assert!(independent_channel > 1e10 * small);
            let doubled = heat.observe_coherence(Some(drive), &[drive], || false).unwrap()[0];
            watts_close(doubled.viscous_power_w, 2.0 * value.viscous_power_w, 1e-14);
            watts_close(doubled.thermal_power_w, 2.0 * value.thermal_power_w, 1e-14);
        }
    }
}

#[test]
fn next_representable_cell_retains_heat_below_surface_work_precision() {
    let gas = gas();
    let (radius, length, omega) = (0.0012, 0.23, 2500.0);
    let wave = crate::segment_wave(&gas, radius, omega, LossModel::Bessel).unwrap();
    let forward = C64::new(1.3, -0.2);
    let backward = C64::new(-0.4, 0.3);
    let field = |x| analytic(
        wave.wavenumber, wave.characteristic_impedance, length, x, forward, backward,
    );
    let sources = [
        pressure_source("left", field(0.0)[0]),
        pressure_source("right", field(length)[0]),
    ];
    let mut net = cylinder_net(length);
    net.edges[0].duct.segments[0] = Segment::Cylinder { radius, length };
    let start = 0.1_f64;
    let end = start.next_up();
    let points: Vec<_> = [start, end].into_iter().map(|offset_m| NetworkFieldProbe {
        edge_id: "tube".into(),
        section: NetworkFieldSection::Cylinder { segment: 0, offset_m },
    }).collect();
    let fields = PreparedNetworkFields::new(
        &net, &gas, &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    ).unwrap();
    let heat = fields.prepare_cylinder_absorption(&[[0, 1]], 1, || false).unwrap();
    let value = heat.sample(fields.port_response().default_excitations(), || false).unwrap()[0];
    let [p, u] = field(start);
    let viscous = 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * u.norm_sq() * (end - start);
    let thermal = 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq() * (end - start);
    assert!(value.viscous_power_w > 0.0 && value.thermal_power_w > 0.0);
    watts_close(value.viscous_power_w, viscous, 2e-9);
    watts_close(value.thermal_power_w, thermal, 2e-9);
    let face_power = (0.5 * (p * u.conj()).re).abs();
    assert!(value.absorbed_power_w < 1e-15 * face_power);
}

#[test]
fn original_liquid_heat_has_exact_zero_thermal_channel_when_cp_equals_cv() {
    use fs_material::fluid_acoustics::FluidAcousticState;
    use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    use fs_material::fluid_transport::FluidTransportState;

    let material = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(950.0, 1300.0).unwrap(),
        FluidTransportState::try_new(950.0, 4200.0, 0.6, 0.2).unwrap(),
        4200.0,
    ).unwrap();
    let (radius, length, omega) = (0.001, 0.19, 1700.0);
    let wave = crate::fluid::segment_wave_with_loss(
        &material, radius, omega, LossModel::Bessel,
    ).unwrap();
    let forward = C64::new(1.0, 0.3);
    let backward = C64::new(0.2, 0.1);
    let field = |x| analytic(
        wave.wavenumber, wave.characteristic_impedance, length, x, forward, backward,
    );
    let sources = [
        pressure_source("left", field(0.0)[0]),
        pressure_source("right", field(length)[0]),
    ];
    let mut net = cylinder_net(length);
    net.edges[0].duct.segments[0] = Segment::Cylinder { radius, length };
    let points = probes(length, &[0.0, 0.37, 1.0]);
    let fields = PreparedNetworkFields::from_material_with_loss(
        &net, &material, crate::ConeSlices::Fixed(3),
        &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    ).unwrap();
    let heat = fields.prepare_cylinder_absorption(&[[0, 1], [1, 2]], 2, || false).unwrap();
    let values = heat.sample(fields.port_response().default_excitations(), || false).unwrap();
    for (value, cell) in values.iter().zip(heat.cells()) {
        let expected = volume_channels(&wave, length, [cell.start_m, cell.end_m], forward, backward);
        assert_eq!(value.thermal_power_w, 0.0);
        assert_eq!(value.absorbed_power_w, value.viscous_power_w);
        watts_close(value.viscous_power_w, expected[0], 2e-9);
    }
    assert_eq!(fields.locations()[0].density_kg_m3, 950.0);
    let p = |x| { let [p, u] = field(x); 0.5 * (p * u.conj()).re };
    watts_close(values.iter().map(|v| v.absorbed_power_w).sum(), p(0.0) - p(length), 2e-9);
}

#[test]
fn cell_domain_budget_and_cancellation_preserve_original_cylinder_identity() {
    let gas = gas();
    let (length, omega) = (0.23, 2500.0);
    let sources = [
        pressure_source("left", C64::ONE),
        pressure_source("right", C64::new(0.2, -0.4)),
    ];
    let mut net = cylinder_net(length);
    net.edges[0].duct.segments.push(Segment::Cylinder { radius: 0.008, length });
    let points = [
        NetworkFieldProbe { edge_id: "tube".into(), section: NetworkFieldSection::Inlet { segment: 0 } },
        NetworkFieldProbe { edge_id: "tube".into(), section: NetworkFieldSection::Outlet { segment: 0 } },
        NetworkFieldProbe { edge_id: "tube".into(), section: NetworkFieldSection::Inlet { segment: 1 } },
    ];
    for loss in [LossModel::Lossless, LossModel::Bessel, LossModel::WideTube] {
        let fields = PreparedNetworkFields::new(
            &net, &gas, &request(omega, loss, &sources, &[], &points), || false,
        ).unwrap();
        assert!(fields.prepare_cylinder_absorption(&[[1, 2]], 1, || false).is_err());
        assert!(fields.prepare_cylinder_absorption(&[[1, 0]], 1, || false).is_err());
        assert!(fields.prepare_cylinder_absorption(&[[0, 3]], 1, || false).is_err());
        assert!(fields.prepare_cylinder_absorption(&[], 1, || false).is_err());
        assert!(fields.prepare_cylinder_absorption(&[[0, 1]], 0, || false).is_err());
        assert!(fields.prepare_cylinder_absorption(&[[0, 1]], MAX_FIELD_PROBES + 1, || false).is_err());
        assert!(matches!(
            fields.prepare_cylinder_absorption(&[[0, 1]], 1, || true),
            Err(DuctError::Cancelled { .. }),
        ));
        if loss == LossModel::WideTube {
            assert!(fields.prepare_cylinder_absorption(&[[0, 1]], 1, || false).is_err());
            continue;
        }
        let heat = fields.prepare_cylinder_absorption(&[[0, 1]], 1, || false).unwrap();
        let drive = fields.port_response().default_excitations();
        assert!(heat.sample(&drive[..1], || false).is_err());
        assert!(heat.sample(&[C64::new(f64::NAN, 0.0), C64::ONE], || false).is_err());
        assert!(matches!(heat.sample(drive, || true), Err(DuctError::Cancelled { .. })));
        assert!(heat.observe_coherence(None, &vec![drive; MAX_FIELD_FACTORS + 1], || false).is_err());
        assert_eq!(
            heat.observe_coherence(None, &[], || false).unwrap(),
            vec![NetworkCylinderAbsorption::default()],
        );
        if loss == LossModel::Lossless {
            assert_eq!(
                heat.sample(drive, || false).unwrap(),
                vec![NetworkCylinderAbsorption::default()],
            );
        }
    }
    let net = cone_net(0.008, 0.012, length);
    let points = cone_points(length, &[0.0, 1.0]);
    let fields = PreparedNetworkFields::new(
        &net, &gas, &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    ).unwrap();
    assert!(fields.prepare_cylinder_absorption(&[[0, 1]], 1, || false).is_err());
}
