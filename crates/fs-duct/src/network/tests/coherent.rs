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
