use super::*;

fn source(label: &str, flow: C64) -> PortFlowSource {
    PortFlowSource {
        label: label.into(),
        volume_flow_m3_s: flow,
    }
}

fn star(lengths: [f64; 2]) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("left"), port("right"), port("receiver")],
        edges: vec![
            edge("left-arm", 0, 2, 0.008, lengths[0]),
            edge("right-arm", 1, 2, 0.008, lengths[1]),
        ],
    }
}

fn zc(gas: &GasState) -> f64 {
    gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008 * 0.008)
}

fn close(a: C64, b: C64, scale: f64) {
    assert!((a - b).abs() <= 3.0e-11 * scale, "{a:?} != {b:?}");
}

// Independent physical pressure/volume-flow equations for two lossless lines:
// U_out,j = U_in,j / cos(kL_j) + i tan(kL_j) p_receiver / Zc_j;
// sum U_out,j = p_receiver / Z_load. No network matrix is used by this oracle.
fn analytic_star(
    gases: [&GasState; 2],
    lengths: [f64; 2],
    omega: f64,
    flows: [C64; 2],
    load: C64,
) -> [C64; 3] {
    let phase = [
        omega * lengths[0] / gases[0].sound_speed,
        omega * lengths[1] / gases[1].sound_speed,
    ];
    let c = phase.map(f64::cos);
    let t = phase.map(f64::tan);
    let receiver = (flows[0].scale(1.0 / c[0]) + flows[1].scale(1.0 / c[1]))
        / (load.recip() - C64::new(0.0, t[0] / zc(gases[0]) + t[1] / zc(gases[1])));
    [
        receiver.scale(1.0 / c[0]) - C64::new(0.0, zc(gases[0]) * t[0]) * flows[0],
        receiver.scale(1.0 / c[1]) - C64::new(0.0, zc(gases[1]) * t[1]) * flows[1],
        receiver,
    ]
}

// G1: simultaneous sources with unequal propagation lengths retain the complete
// complex transfer and signed work. This compares physical fields, not merely
// equality to a sum of the same implementation's single-source queries.
#[test]
fn coherent_sources_match_loaded_line_equations_and_work() {
    let gas = gas();
    let lengths = [0.11, 0.17];
    let net = star(lengths);
    let load = C64::new(0.7 * zc(&gas), -0.3 * zc(&gas));
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    for omega in [713.0, 2301.0, 5107.0] {
        for angle in [0.0_f64, 0.73, core::f64::consts::PI] {
            let flows = [
                C64::new(2.0e-6, -1.0e-6),
                C64::new(angle.cos(), angle.sin()).scale(1.3e-6),
            ];
            let sources = [source("left", flows[0]), source("right", flows[1])];
            let response = network_drive_coherent(
                &net,
                &gas,
                omega,
                LossModel::Lossless,
                &sources,
                &loads,
                None,
            )
            .unwrap();
            let expected = analytic_star([&gas; 2], lengths, omega, flows, load);
            for (&actual, &wanted) in response.port_pressures.iter().zip(&expected) {
                close(actual, wanted, zc(&gas) * 3.0e-6);
            }
            assert_eq!(response.port_inflows[0], flows[0]);
            assert_eq!(response.port_inflows[1], flows[1]);
            close(response.port_inflows[2], -expected[2] / load, 3.0e-6);
            let expected_power = 0.5 * (expected[2] * (expected[2] / load).conj()).re;
            let tolerance = 3.0e-11 * zc(&gas) * 9.0e-12;
            assert!((response.external_load_power_w - expected_power).abs() < tolerance);
            assert!((response.source_power_w - expected_power).abs() < tolerance);
            assert!(response.network_power_w.abs() < tolerance);
        }
    }
}

// G1: the actual loaded system exists at its open-circuit pole. Equal-phase
// sources quadruple one source's received power; opposite phases cancel it.
// A weaker opposing source absorbs work and must not be clipped to zero.
#[test]
fn coherent_half_wave_sources_cancel_and_exchange_signed_power() {
    let gas = gas();
    let net = star([0.2; 2]);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(zc(&gas))),
    }];
    let amplitude = 1.0e-6;
    let one_source_power = 0.5 * zc(&gas) * amplitude * amplitude;
    for second in [0.0, 1.0, -1.0, -0.5] {
        let sources = [
            source("left", C64::from_re(amplitude)),
            source("right", C64::from_re(second * amplitude)),
        ];
        let response = network_drive_coherent(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &sources,
            &loads,
            None,
        )
        .unwrap();
        let total = amplitude * (1.0 + second);
        close(
            response.port_pressures[2],
            C64::from_re(-zc(&gas) * total),
            zc(&gas) * amplitude,
        );
        close(response.port_inflows[2], C64::from_re(total), amplitude);
        let expected = one_source_power * (1.0 + second).powi(2);
        assert!((response.external_load_power_w - expected).abs() < 3.0e-12 * one_source_power);
        assert!((response.source_power_w - expected).abs() < 3.0e-12 * one_source_power);
        assert!(response.network_power_w.abs() < 3.0e-12 * one_source_power);
        if second == -0.5 {
            assert!(response.port_powers_w[1] < -0.24 * one_source_power);
        }
        if second == -1.0 {
            assert!(response.port_pressures[2].abs() < 3.0e-12 * zc(&gas) * amplitude);
        }
    }
}

// G1: off-diagonal radiation acts on the total coherent source field. With
// half-wave arms the independent two-mouth load gives U=[1/4,3/4] sum(U_src)
// and p=(7/8) Zc sum(U_src), including the Re(U1 conj(U2)) power term.
#[test]
fn coherent_sources_close_mutual_radiation_at_open_circuit_pole() {
    let gas = gas();
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let net = DuctNetwork {
        nodes: vec![
            port("source-a"),
            port("source-b"),
            node("junction", NodeBoundary::Junction),
            port("mouth-a"),
            port("mouth-b"),
        ],
        edges: vec![
            edge("source-a-line", 0, 2, 0.008, 0.2),
            edge("source-b-line", 1, 2, 0.008, 0.2),
            edge("mouth-a-line", 2, 3, 0.008, 0.2),
            edge("mouth-b-line", 2, 4, 0.008, 0.2),
        ],
    };
    let flows = [C64::new(1.0e-6, 0.2e-6), C64::new(-0.4e-6, 0.7e-6)];
    let sources = [source("source-a", flows[0]), source("source-b", flows[1])];
    let coupled = CoupledTermination {
        omega,
        labels: vec!["mouth-a".into(), "mouth-b".into()],
        z: [2.0, 0.5, 0.5, 1.0]
            .map(|value| C64::from_re(value * zc(&gas)))
            .to_vec(),
        source_id: "analytic reciprocal two-mouth radiation".into(),
    };
    let response = network_drive_coherent(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &sources,
        &[],
        Some(&coupled),
    )
    .unwrap();
    let total = flows[0] + flows[1];
    for &pressure in &response.port_pressures {
        close(pressure, total.scale(0.875 * zc(&gas)), zc(&gas) * 2.0e-6);
    }
    close(response.port_inflows[2], total.scale(-0.25), 2.0e-6);
    close(response.port_inflows[3], total.scale(-0.75), 2.0e-6);
    let expected = 0.5 * 0.875 * zc(&gas) * total.norm_sq();
    assert!((response.source_power_w - expected).abs() < 3.0e-11 * expected);
    assert!((response.external_load_power_w - expected).abs() < 3.0e-11 * expected);
    assert!(response.network_power_w.abs() < 3.0e-11 * expected);
}

// G3/G5: common phase only rotates fields, all source-order permutations use
// the same sums, and the original admitted one-unit-source result keeps bits.
#[test]
fn coherent_phase_order_and_single_source_compatibility() {
    let gas = gas();
    let net = star([0.13, 0.19]);
    let load = PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(zc(&gas))),
    };
    let operator = network_port_impedance(&net, &gas, 2731.0, LossModel::Bessel).unwrap();
    let sources = [
        source("left", C64::new(1.2e-6, 0.7e-6)),
        source("right", C64::new(-0.3e-6, 0.8e-6)),
    ];
    let first = operator
        .drive_coherent(&sources, &[load.clone()], None)
        .unwrap();
    assert!(first.network_power_w > 0.0);
    let swapped = [sources[1].clone(), sources[0].clone()];
    let reordered = operator
        .drive_coherent(&swapped, &[load.clone()], None)
        .unwrap();
    for (a, b) in first.port_pressures.iter().zip(&reordered.port_pressures) {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
    let phase = C64::new(0.6, 0.8);
    let rotated = sources.map(|mut value| {
        value.volume_flow_m3_s = value.volume_flow_m3_s * phase;
        value
    });
    let second = operator
        .drive_coherent(&rotated, &[load.clone()], None)
        .unwrap();
    for (a, b) in first.port_pressures.iter().zip(&second.port_pressures) {
        close(*a * phase, *b, zc(&gas) * 2.0e-6);
    }
    assert!((second.source_power_w - first.source_power_w).abs() < 3.0e-12 * first.source_power_w);
    let loads = [
        PortTermination {
            label: "right".into(),
            load: NetworkLoad::Closed,
        },
        load,
    ];
    let old = operator.drive_bound("left", &loads).unwrap();
    let new = operator
        .drive_coherent(&[source("left", C64::ONE)], &loads, None)
        .unwrap();
    for (a, b) in old
        .port_pressures
        .iter()
        .chain(&old.port_inflows)
        .zip(new.port_pressures.iter().chain(&new.port_inflows))
    {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
}

// G1: each arm's actual medium controls its own phase and impedance, so a
// common-source cancellation in equal arms need not survive a gas change.
#[test]
fn coherent_sources_keep_regional_phase_and_impedance() {
    use fs_material::gas_species::input::GasRecord;
    let gas = gas();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 373.15, gas.pressure).unwrap();
    let lengths = [0.17; 2];
    let net = star(lengths);
    let media = NetworkMedia::new(
        &net,
        vec![
            ("left-arm".into(), GasRecord::Legacy(gas)),
            ("right-arm".into(), GasRecord::Legacy(hot)),
        ],
        GasRecord::Legacy(gas),
    )
    .unwrap();
    let flows = [C64::from_re(1.0e-6), C64::from_re(-1.0e-6)];
    let sources = [source("left", flows[0]), source("right", flows[1])];
    let load = C64::from_re(zc(&gas));
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    let omega = 2397.0;
    let response = network_drive_coherent_regional(
        &net,
        &media,
        omega,
        LossModel::Lossless,
        &sources,
        &loads,
        None,
    )
    .unwrap();
    let expected = analytic_star([&gas, &hot], lengths, omega, flows, load);
    for (a, b) in response.port_pressures.iter().zip(&expected) {
        close(*a, *b, zc(&gas) * 1.0e-6);
    }
    assert!(response.port_pressures[2].abs() > 0.01 * zc(&gas) * 1.0e-6);
}

#[test]
fn coherent_sources_require_a_complete_finite_boundary_partition() {
    let gas = gas();
    let net = star([0.13, 0.19]);
    let operator = network_port_impedance(&net, &gas, 2531.0, LossModel::Lossless).unwrap();
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Closed,
    }];
    for sources in [
        vec![],
        vec![source("left", C64::ONE)],
        vec![source("left", C64::ONE), source("left", C64::ONE)],
        vec![source("left", C64::ONE), source("absent", C64::ONE)],
        vec![source("left", C64::ONE), source("receiver", C64::ONE)],
        vec![
            source("left", C64::ONE),
            source("right", C64::new(f64::NAN, 0.0)),
        ],
        vec![
            source("left", C64::ONE),
            source("right", C64::new(0.0, f64::INFINITY)),
        ],
    ] {
        assert!(matches!(
            operator.drive_coherent(&sources, &loads, None),
            Err(DuctError::BadParameter { .. })
        ));
    }
    let all = [
        source("left", C64::ZERO),
        source("right", C64::ZERO),
        source("receiver", C64::ZERO),
    ];
    let zero = operator.drive_coherent(&all, &[], None).unwrap();
    assert!(zero.port_pressures.iter().all(|value| *value == C64::ZERO));
    assert!(zero.port_inflows.iter().all(|value| *value == C64::ZERO));
    assert_eq!(zero.source_power_w, 0.0);
    assert_eq!(zero.network_power_w, 0.0);
}

fn pressure_source(label: &str, pressure: C64, impedance: C64) -> PortSource {
    PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: pressure,
            impedance_pa_s_m3: impedance,
        },
    }
}

fn source_work_balance(response: &NetworkSourceResponse, scale_w: f64) {
    assert!(
        (response.generator_power_w
            - response.source_dissipation_w
            - response.fields.source_power_w)
            .abs()
            < 1.0e-10 * scale_w,
        "generator/internal/interface power mismatch: {response:?}"
    );
    assert!(
        response
            .source_dissipations_w
            .iter()
            .all(|&power| power >= 0.0)
    );
}

// G1: a matched lossless line has input impedance Zc at every phase. Its
// pressure advances by exp(+ikL), including exact half-wave open-port poles.
// The source divider P_open/(Zs+Zc) is independent of the graph solver.
#[test]
fn finite_source_matches_pressure_divider_and_ideal_pressure_limit() {
    let gas = gas();
    let length = 0.21;
    let net = DuctNetwork {
        nodes: vec![port("generator"), port("mouth")],
        edges: vec![edge("line", 0, 1, 0.008, length)],
    };
    let resistance = zc(&gas);
    let loads = [PortTermination {
        label: "mouth".into(),
        load: NetworkLoad::Impedance(C64::from_re(resistance)),
    }];
    let pressure = C64::new(1.7, -0.4);
    for phase in [0.37, core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
        let omega = phase * gas.sound_speed / length;
        for impedance in [C64::ZERO, C64::new(0.4, -0.7).scale(resistance)] {
            let response = network_drive_sources(
                &net,
                &gas,
                omega,
                LossModel::Lossless,
                &[pressure_source("generator", pressure, impedance)],
                &loads,
                None,
            )
            .unwrap();
            let flow = pressure / (impedance + C64::from_re(resistance));
            let pin = flow.scale(resistance);
            let pout = pin * C64::new(phase.cos(), phase.sin());
            close(
                response.fields.port_inflows[0],
                flow,
                pressure.abs() / resistance,
            );
            close(response.fields.port_pressures[0], pin, pressure.abs());
            close(response.fields.port_pressures[1], pout, pressure.abs());
            close(
                response.fields.port_inflows[1],
                -pout.scale(1.0 / resistance),
                flow.abs(),
            );
            source_work_balance(&response, pressure.norm_sq() / resistance);
            assert!(
                response.fields.network_power_w.abs() < 1.0e-11 * pressure.norm_sq() / resistance
            );
        }
    }
}

// Independent two-line elimination with each p_j+Zs_j U_j=P_j. Unlike an
// ideal flow source, U_j changes when the other generator's phase changes.
fn analytic_pressure_star(
    gases: [&GasState; 2],
    lengths: [f64; 2],
    omega: f64,
    pressures: [C64; 2],
    impedances: [C64; 2],
    load: C64,
) -> ([C64; 3], [C64; 3]) {
    let phase = [
        omega * lengths[0] / gases[0].sound_speed,
        omega * lengths[1] / gases[1].sound_speed,
    ];
    let c = phase.map(f64::cos);
    let t = phase.map(f64::tan);
    let a = [
        impedances[0] - C64::new(0.0, zc(gases[0]) * t[0]),
        impedances[1] - C64::new(0.0, zc(gases[1]) * t[1]),
    ];
    let receiver = (pressures[0] / a[0].scale(c[0]) + pressures[1] / a[1].scale(c[1]))
        / (load.recip()
            + a[0].recip().scale(1.0 / (c[0] * c[0]))
            + a[1].recip().scale(1.0 / (c[1] * c[1]))
            - C64::new(0.0, t[0] / zc(gases[0]) + t[1] / zc(gases[1])));
    let flows = [
        (pressures[0] - receiver.scale(1.0 / c[0])) / a[0],
        (pressures[1] - receiver.scale(1.0 / c[1])) / a[1],
        -receiver / load,
    ];
    (
        [
            pressures[0] - impedances[0] * flows[0],
            pressures[1] - impedances[1] * flows[1],
            receiver,
        ],
        flows,
    )
}

#[test]
fn finite_sources_match_independent_line_equations_and_react_to_each_other() {
    let gas = gas();
    let lengths = [0.11, 0.17];
    let net = star(lengths);
    let impedances = [
        C64::new(0.4, -0.3).scale(zc(&gas)),
        C64::new(0.8, 0.2).scale(zc(&gas)),
    ];
    let load = C64::new(0.7, -0.1).scale(zc(&gas));
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    for omega in [713.0, 2301.0, 5107.0] {
        let mut first_flow = C64::ZERO;
        for angle in [0.0_f64, 0.73, core::f64::consts::PI] {
            let pressures = [
                C64::new(2.0, -0.3),
                C64::new(angle.cos(), angle.sin()).scale(1.3),
            ];
            let sources = [
                pressure_source("left", pressures[0], impedances[0]),
                pressure_source("right", pressures[1], impedances[1]),
            ];
            let response = network_drive_sources(
                &net,
                &gas,
                omega,
                LossModel::Lossless,
                &sources,
                &loads,
                None,
            )
            .unwrap();
            let (expected_p, expected_u) =
                analytic_pressure_star([&gas; 2], lengths, omega, pressures, impedances, load);
            for port in 0..3 {
                close(response.fields.port_pressures[port], expected_p[port], 3.0);
                close(
                    response.fields.port_inflows[port],
                    expected_u[port],
                    3.0 / zc(&gas),
                );
            }
            source_work_balance(&response, 9.0 / zc(&gas));
            assert!(response.fields.network_power_w.abs() < 1.0e-10 * 9.0 / zc(&gas));
            if angle == 0.0 {
                first_flow = response.fields.port_inflows[0];
            } else {
                assert!((response.fields.port_inflows[0] - first_flow).abs() > 0.01 / zc(&gas));
            }
        }
    }
}

// G1: anti-phase finite generators still dissipate energy in their source
// resistances even when their external radiation cancels. A zero drive stays
// connected and absorbs; a weak generator can receive signed negative work.
#[test]
fn finite_sources_cancel_radiation_without_erasing_internal_work() {
    let gas = gas();
    let net = star([0.2; 2]);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let resistance = zc(&gas);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(resistance)),
    }];
    for second in [0.0, 1.0, -1.0, 0.1] {
        let pressures = [C64::ONE, C64::from_re(second)];
        let sources = [
            pressure_source("left", pressures[0], C64::from_re(resistance)),
            pressure_source("right", pressures[1], C64::from_re(resistance)),
        ];
        let response = network_drive_sources(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &sources,
            &loads,
            None,
        )
        .unwrap();
        let interface = C64::from_re((1.0 + second) / 3.0);
        for port in 0..2 {
            close(response.fields.port_pressures[port], interface, 1.0);
            close(
                response.fields.port_inflows[port],
                (pressures[port] - interface).scale(1.0 / resistance),
                1.0 / resistance,
            );
        }
        close(response.fields.port_pressures[2], -interface, 1.0);
        source_work_balance(&response, 1.0 / resistance);
        if second == -1.0 {
            assert!(response.fields.external_load_power_w.abs() < 1.0e-24 / resistance);
            assert!(
                (response.source_dissipation_w - 1.0 / resistance).abs() < 1.0e-11 / resistance
            );
        }
        if second == 0.0 || second == 0.1 {
            assert!(response.fields.port_powers_w[1] < 0.0);
            assert!(response.source_dissipations_w[1] > 0.0);
        }
        if second == 0.1 {
            assert!(response.generator_powers_w[1] < 0.0);
        }
    }
}

#[test]
fn finite_sources_retain_mutual_loads_and_mixed_flow_drive_at_poles() {
    let gas = gas();
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let resistance = zc(&gas);
    let net = DuctNetwork {
        nodes: vec![
            port("source-a"),
            port("source-b"),
            node("junction", NodeBoundary::Junction),
            port("mouth-a"),
            port("mouth-b"),
        ],
        edges: vec![
            edge("a", 0, 2, 0.008, 0.2),
            edge("b", 1, 2, 0.008, 0.2),
            edge("c", 2, 3, 0.008, 0.2),
            edge("d", 2, 4, 0.008, 0.2),
        ],
    };
    let coupled = CoupledTermination {
        omega,
        labels: vec!["mouth-a".into(), "mouth-b".into()],
        z: [2.0, 0.5, 0.5, 1.0]
            .map(|v| C64::from_re(v * resistance))
            .to_vec(),
        source_id: "analytic two-mouth mutual load".into(),
    };
    let fixed_flow = C64::new(0.3, 0.2).scale(1.0 / resistance);
    let pressure = C64::new(1.1, -0.3);
    let impedance = C64::new(0.7, -0.2).scale(resistance);
    let sources = [
        PortSource {
            label: "source-a".into(),
            excitation: PortExcitation::VolumeFlow(fixed_flow),
        },
        pressure_source("source-b", pressure, impedance),
    ];
    let response = network_drive_sources(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &sources,
        &[],
        Some(&coupled),
    )
    .unwrap();
    let equivalent_load = C64::from_re(0.875 * resistance);
    let expected_pressure =
        (fixed_flow + pressure / impedance) / (equivalent_load.recip() + impedance.recip());
    let total_flow = expected_pressure / equivalent_load;
    for &actual in &response.fields.port_pressures {
        close(actual, expected_pressure, 2.0);
    }
    close(
        response.fields.port_inflows[0],
        fixed_flow,
        1.0 / resistance,
    );
    close(
        response.fields.port_inflows[1],
        (pressure - expected_pressure) / impedance,
        2.0 / resistance,
    );
    close(
        response.fields.port_inflows[2],
        total_flow.scale(-0.25),
        2.0 / resistance,
    );
    close(
        response.fields.port_inflows[3],
        total_flow.scale(-0.75),
        2.0 / resistance,
    );
    source_work_balance(&response, 4.0 / resistance);
    assert!(response.fields.network_power_w.abs() < 1.0e-10 / resistance);
}

// G1/G3/G5: source impedances remain in SI through regional interfaces, and
// changing only input order cannot change the canonical solve or work sums.
#[test]
fn finite_sources_preserve_regional_media_and_source_order() {
    use fs_material::gas_species::input::GasRecord;
    let gas = gas();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 373.15, gas.pressure).unwrap();
    let net = star([0.17; 2]);
    let media = NetworkMedia::new(
        &net,
        vec![
            ("left-arm".into(), GasRecord::Legacy(gas)),
            ("right-arm".into(), GasRecord::Legacy(hot)),
        ],
        GasRecord::Legacy(gas),
    )
    .unwrap();
    let pressures = [C64::new(1.2, 0.3), C64::new(-0.7, 0.2)];
    let impedances = [
        C64::from_re(0.7 * zc(&gas)),
        C64::new(0.4, 0.6).scale(zc(&hot)),
    ];
    let sources = [
        pressure_source("left", pressures[0], impedances[0]),
        pressure_source("right", pressures[1], impedances[1]),
    ];
    let load = C64::from_re(zc(&gas));
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    let omega = 2397.0;
    let response = network_drive_sources_regional(
        &net,
        &media,
        omega,
        LossModel::Lossless,
        &sources,
        &loads,
        None,
    )
    .unwrap();
    let (expected_p, expected_u) =
        analytic_pressure_star([&gas, &hot], [0.17; 2], omega, pressures, impedances, load);
    for port in 0..3 {
        close(response.fields.port_pressures[port], expected_p[port], 2.0);
        close(
            response.fields.port_inflows[port],
            expected_u[port],
            2.0 / zc(&gas),
        );
    }
    let swapped = [sources[1].clone(), sources[0].clone()];
    let reordered = network_drive_sources_regional(
        &net,
        &media,
        omega,
        LossModel::Lossless,
        &swapped,
        &loads,
        None,
    )
    .unwrap();
    for (a, b) in response
        .fields
        .port_pressures
        .iter()
        .chain(&response.fields.port_inflows)
        .zip(
            reordered
                .fields
                .port_pressures
                .iter()
                .chain(&reordered.fields.port_inflows),
        )
    {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
    source_work_balance(&response, 4.0 / zc(&gas));
}

// The loaded boundary can have 32 unknown source flows and no passive ports.
// At half-wave the centre pressure is P_open/32 and every tip is its negative.
#[test]
fn finite_sources_allow_every_port_to_be_driven_and_refuse_invalid_boundaries() {
    let gas = gas();
    let resistance = zc(&gas);
    let mut net = DuctNetwork {
        nodes: vec![port("centre")],
        edges: vec![],
    };
    let mut sources = vec![pressure_source(
        "centre",
        C64::ONE,
        C64::from_re(resistance),
    )];
    for index in 1..MAX_NETWORK_PORTS {
        let label = format!("tip-{index}");
        net.nodes.push(port(&label));
        net.edges.push(edge(&label, 0, index, 0.008, 0.2));
        sources.push(pressure_source(&label, C64::ZERO, C64::from_re(resistance)));
    }
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let response =
        network_drive_sources(&net, &gas, omega, LossModel::Lossless, &sources, &[], None).unwrap();
    close(
        response.fields.port_pressures[0],
        C64::from_re(1.0 / 32.0),
        1.0,
    );
    for &pressure in &response.fields.port_pressures[1..] {
        close(pressure, C64::from_re(-1.0 / 32.0), 1.0);
    }
    source_work_balance(&response, 1.0 / resistance);
    assert_eq!(response.fields.external_load_power_w, 0.0);
    assert!(response.fields.network_power_w.abs() < 1.0e-10 / resistance);

    for excitation in [
        PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::new(f64::NAN, 0.0),
            impedance_pa_s_m3: C64::ONE,
        },
        PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::ONE,
            impedance_pa_s_m3: C64::from_re(-1.0),
        },
        PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::ONE,
            impedance_pa_s_m3: C64::new(0.0, f64::INFINITY),
        },
    ] {
        sources[0].excitation = excitation;
        assert!(matches!(
            network_drive_sources(&net, &gas, omega, LossModel::Lossless, &sources, &[], None),
            Err(DuctError::BadParameter { .. })
        ));
    }
}
