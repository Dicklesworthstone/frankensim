use super::*;

pub(super) fn zc(gas: &GasState) -> f64 {
    gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008 * 0.008)
}

pub(super) fn star(lengths: [f64; 2]) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("left"), port("right"), port("receiver")],
        edges: vec![
            edge("left-arm", 0, 2, 0.008, lengths[0]),
            edge("right-arm", 1, 2, 0.008, lengths[1]),
        ],
    }
}

pub(super) fn bank(omega: f64, z: [C64; 4], pressures: [C64; 2]) -> CoupledSourceBank {
    CoupledSourceBank {
        impedance: CoupledTermination {
            omega,
            labels: vec!["left".into(), "right".into()],
            z: z.to_vec(),
            source_id: "analytic reciprocal driver pair".into(),
        },
        open_circuit_pressures_pa: pressures.to_vec(),
    }
}

pub(super) fn close(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() < 1.0e-10 * scale,
        "actual={actual:?}, expected={expected:?}"
    );
}

pub(super) fn solve_pair(a: [C64; 4], b: [C64; 2]) -> [C64; 2] {
    let determinant = a[0] * a[3] - a[1] * a[2];
    [
        (a[3] * b[0] - a[1] * b[1]) / determinant,
        (a[0] * b[1] - a[2] * b[0]) / determinant,
    ]
}

// Independent lossless line equations. Each line has
// p_i = p_r/cos(kL_i) - i Zc_i tan(kL_i) U_i, and the two outward
// line flows sum to p_r/Zload. No graph impedance or factorization is used.
fn line_fields(
    gases: [&GasState; 2],
    lengths: [f64; 2],
    omega: f64,
    u: [C64; 2],
    load: C64,
) -> ([C64; 3], [C64; 3]) {
    let phase = [
        omega * lengths[0] / gases[0].sound_speed,
        omega * lengths[1] / gases[1].sound_speed,
    ];
    let c = phase.map(f64::cos);
    let t = phase.map(f64::tan);
    let receiver = (u[0].scale(1.0 / c[0]) + u[1].scale(1.0 / c[1]))
        / (load.recip() - C64::new(0.0, t[0] / zc(gases[0]) + t[1] / zc(gases[1])));
    (
        [
            receiver.scale(1.0 / c[0]) - C64::new(0.0, zc(gases[0]) * t[0]) * u[0],
            receiver.scale(1.0 / c[1]) - C64::new(0.0, zc(gases[1]) * t[1]) * u[1],
            receiver,
        ],
        [u[0], u[1], -receiver / load],
    )
}

pub(super) fn analytic_bank(
    gases: [&GasState; 2],
    lengths: [f64; 2],
    omega: f64,
    z: [C64; 4],
    pressure: [C64; 2],
    load: C64,
) -> ([C64; 3], [C64; 3]) {
    let phase = [
        omega * lengths[0] / gases[0].sound_speed,
        omega * lengths[1] / gases[1].sound_speed,
    ];
    let c = phase.map(f64::cos);
    let t = phase.map(f64::tan);
    let receiver =
        (load.recip() - C64::new(0.0, t[0] / zc(gases[0]) + t[1] / zc(gases[1]))).recip();
    let mut a = z;
    for row in 0..2 {
        for column in 0..2 {
            a[row * 2 + column] = a[row * 2 + column] + receiver.scale(1.0 / (c[row] * c[column]));
        }
        a[row * 2 + row] = a[row * 2 + row] - C64::new(0.0, zc(gases[row]) * t[row]);
    }
    line_fields(gases, lengths, omega, solve_pair(a, pressure), load)
}

pub(super) fn balance(response: &NetworkSourceResponse, scale: f64) {
    assert!(response.source_dissipation_w >= -1.0e-11 * scale);
    assert!(
        (response.generator_power_w
            - response.source_dissipation_w
            - response.fields.external_load_power_w
            - response.fields.network_power_w)
            .abs()
            < 1.0e-10 * scale,
        "{response:?}"
    );
}

// Signed resistive mutual coupling and reactive coupling survive unequal
// propagation and independent regional sound speeds, in physical SI units.
#[test]
fn source_bank_matches_independent_two_line_elimination_and_regional_media() {
    use fs_material::gas_species::input::GasRecord;
    let gas = gas();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 373.15, gas.pressure).unwrap();
    let lengths = [0.11, 0.17];
    let net = star(lengths);
    let r = zc(&gas);
    let z = [
        C64::new(0.7, -0.2),
        C64::new(-0.35, 0.6),
        C64::new(-0.35, 0.6),
        C64::new(0.9, 0.4),
    ]
    .map(|value| value.scale(r));
    let load = C64::new(0.8, -0.3).scale(r);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    for right_gas in [gas, hot] {
        let media = NetworkMedia::new(
            &net,
            vec![
                ("left-arm".into(), GasRecord::Legacy(gas)),
                ("right-arm".into(), GasRecord::Legacy(right_gas)),
            ],
            GasRecord::Legacy(gas),
        )
        .unwrap();
        for omega in [713.0, 2301.0, 5107.0] {
            let pressure = [C64::new(1.2, -0.3), C64::new(-0.4, 0.8)];
            let bank = bank(omega, z, pressure);
            let response = network_drive_source_bank_regional(
                &net,
                &media,
                omega,
                LossModel::Lossless,
                &[],
                &bank,
                &loads,
                None,
            )
            .unwrap();
            let (p, u) = analytic_bank([&gas, &right_gas], lengths, omega, z, pressure, load);
            for port in 0..3 {
                close(response.fields.port_pressures[port], p[port], 3.0);
                close(response.fields.port_inflows[port], u[port], 3.0 / r);
            }
            balance(&response, 10.0 / r);
            assert!(response.fields.network_power_w.abs() < 1.0e-10 / r);
        }
    }
}

#[test]
fn diagonal_source_bank_recovers_independent_pressure_sources_and_permutation() {
    let gas = gas();
    let net = star([0.13, 0.19]);
    let omega = 2531.0;
    let r = zc(&gas);
    let p = [C64::new(1.3, 0.2), C64::new(-0.8, 0.1)];
    let z = [C64::new(0.7, 0.2).scale(r), C64::new(1.1, -0.3).scale(r)];
    let bank = bank(omega, [z[0], C64::ZERO, C64::ZERO, z[1]], p);
    let sources = ["left", "right"].map(|label| {
        let i = usize::from(label == "right");
        PortSource {
            label: label.into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: p[i],
                impedance_pa_s_m3: z[i],
            },
        }
    });
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let operator = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    let scalar = operator.drive_sources(&sources, &loads, None).unwrap();
    let diagonal = operator
        .drive_source_bank(&[], &bank, &loads, None)
        .unwrap();
    for port in 0..3 {
        close(
            diagonal.fields.port_pressures[port],
            scalar.fields.port_pressures[port],
            3.0,
        );
        close(
            diagonal.fields.port_inflows[port],
            scalar.fields.port_inflows[port],
            3.0 / r,
        );
    }
    assert!((diagonal.source_dissipation_w - scalar.source_dissipation_w).abs() < 1.0e-12 / r);

    let mut full = bank.clone();
    full.impedance.z[1] = C64::new(-0.3, 0.5).scale(r);
    full.impedance.z[2] = full.impedance.z[1];
    let response = operator
        .drive_source_bank(&[], &full, &loads, None)
        .unwrap();
    let mut swapped = full.clone();
    swapped.impedance.labels.reverse();
    swapped.impedance.z.reverse();
    swapped.open_circuit_pressures_pa.reverse();
    let reordered = operator
        .drive_source_bank(&[], &swapped, &loads, None)
        .unwrap();
    assert_eq!(
        response.fields.port_pressures,
        reordered.fields.port_pressures
    );
    assert_eq!(response.fields.port_inflows, reordered.fields.port_inflows);
    assert_eq!(
        response.source_dissipations_w,
        reordered.source_dissipations_w
    );
    assert_eq!(
        response.source_dissipation_w,
        reordered.source_dissipation_w
    );
    let rotation = C64::new(0.37_f64.cos(), 0.37_f64.sin());
    for pressure in &mut swapped.open_circuit_pressures_pa {
        *pressure = *pressure * rotation;
    }
    let rotated = operator
        .drive_source_bank(&[], &swapped, &loads, None)
        .unwrap();
    for port in 0..3 {
        close(
            rotated.fields.port_pressures[port],
            response.fields.port_pressures[port] * rotation,
            3.0,
        );
    }
    assert!((rotated.source_dissipation_w - response.source_dissipation_w).abs() < 1.0e-11 / r);
}

// Half-wave graph pressure cancels while the differential source mode still
// dissipates. The actual boundary regularizes the open-circuit graph pole.
#[test]
fn mutually_coupled_sources_cancel_radiation_without_cancelling_internal_loss() {
    let gas = gas();
    let net = star([0.2; 2]);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let r = zc(&gas);
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    let bank = bank(
        omega,
        [2.0, -0.5, -0.5, 2.0].map(|v| C64::from_re(v * r)),
        [C64::ONE, -C64::ONE],
    );
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let response = network_drive_source_bank(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &[],
        &bank,
        &loads,
        None,
    )
    .unwrap();
    close(
        response.fields.port_inflows[0],
        C64::from_re(1.0 / (2.5 * r)),
        1.0 / r,
    );
    close(
        response.fields.port_inflows[1],
        C64::from_re(-1.0 / (2.5 * r)),
        1.0 / r,
    );
    assert!(response.fields.external_load_power_w.abs() < 1.0e-23 / r);
    assert!((response.source_dissipation_w - 1.0 / (2.5 * r)).abs() < 1.0e-10 / r);
    balance(&response, 1.0 / r);
}

// A passive bank can transfer work between its ports. Neither signed real
// mutual impedance nor real work transferred by reactive coupling may vanish.
#[test]
fn source_bank_preserves_negative_port_exchange_and_reactive_cancellation() {
    let gas = gas();
    let r = zc(&gas);
    let lengths = [0.12, 0.16];
    let net = star(lengths);
    let omega = 1703.0;
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let cases = [
        (
            [
                C64::from_re(1.0),
                C64::from_re(2.0),
                C64::from_re(2.0),
                C64::from_re(4.0),
            ],
            [C64::from_re(1.0e-6), C64::from_re(-0.3e-6)],
            0.08e-12,
        ),
        (
            [
                C64::new(0.0, 0.3),
                C64::new(0.0, 2.0),
                C64::new(0.0, 2.0),
                C64::new(0.0, -0.7),
            ],
            [C64::from_re(1.0e-6), C64::new(0.0, 1.0e-6)],
            0.0,
        ),
        (
            [
                C64::new(1.0, 0.0),
                C64::new(0.0, 1.0e6),
                C64::new(0.0, 1.0e6),
                C64::new(1.0, 0.0),
            ],
            [C64::from_re(1.0e-6), C64::new(0.0, 1.0e-6)],
            1.0e-12,
        ),
    ];
    for (z, u, loss) in cases {
        let z = z.map(|value| value.scale(r));
        let (p, _) = line_fields([&gas; 2], lengths, omega, u, C64::from_re(r));
        let bank = bank(
            omega,
            z,
            [
                p[0] + z[0] * u[0] + z[1] * u[1],
                p[1] + z[2] * u[0] + z[3] * u[1],
            ],
        );
        let response = network_drive_source_bank(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &[],
            &bank,
            &loads,
            None,
        )
        .unwrap();
        for (port, flow) in u.iter().enumerate() {
            close(response.fields.port_inflows[port], *flow, 2.0e-6);
            let internal = z[port * 2] * u[0] + z[port * 2 + 1] * u[1];
            let expected = 0.5 * (internal * flow.conj()).re;
            assert!(
                (response.source_dissipations_w[port] - expected).abs()
                    < 1.0e-10 * (expected.abs() + 1.0e-12 * r)
            );
        }
        assert!(
            response
                .source_dissipations_w
                .iter()
                .any(|power| *power < 0.0)
        );
        assert!((response.source_dissipation_w - loss * r).abs() < 1.0e-22 * r);
        if loss == 0.0 {
            assert_eq!(response.source_dissipation_w, 0.0);
        }
    }
}

#[test]
fn coupled_source_bank_coexists_with_mixed_sources_and_full_exterior_load_at_pole() {
    let gas = gas();
    let r = zc(&gas);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let labels = ["left", "right", "scalar", "flow", "mouth-a", "mouth-b"];
    let mut nodes: Vec<NetworkNode> = labels.iter().map(|label| port(label)).collect();
    nodes.push(node("junction", NodeBoundary::Junction));
    let net = DuctNetwork {
        nodes,
        edges: labels
            .iter()
            .enumerate()
            .map(|(i, label)| edge(label, i, 6, 0.008, 0.2))
            .collect(),
    };
    let zs = [
        C64::new(1.2, 0.2),
        C64::new(-0.4, 0.3),
        C64::new(-0.4, 0.3),
        C64::new(0.8, -0.1),
    ]
    .map(|v| v.scale(r));
    let pressures = [C64::new(1.1, 0.2), C64::new(-0.3, 0.4)];
    let bank = bank(omega, zs, pressures);
    let zl = [2.0, 0.5, 0.5, 1.0].map(|value| C64::from_re(value * r));
    let exterior = CoupledTermination {
        omega,
        labels: vec!["mouth-a".into(), "mouth-b".into()],
        z: zl.to_vec(),
        source_id: "independent mutual receiver pair".into(),
    };
    let scalar_pressure = C64::new(0.7, -0.2);
    let scalar_impedance = C64::new(0.6, 0.3).scale(r);
    let flow = C64::new(0.3, 0.2).scale(1.0 / r);
    let sources = [
        PortSource {
            label: "scalar".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: scalar_pressure,
                impedance_pa_s_m3: scalar_impedance,
            },
        },
        PortSource {
            label: "flow".into(),
            excitation: PortExcitation::VolumeFlow(flow),
        },
    ];
    // All exposed ports have the same pressure after an equal half-wave arm.
    // One conservation equation closes the independent 2x2 source/load laws.
    let injected = solve_pair(zs, pressures);
    let source_y = solve_pair(zs, [C64::ONE; 2]);
    let load_y = solve_pair(zl, [C64::ONE; 2]);
    let p = (injected[0] + injected[1] + scalar_pressure / scalar_impedance + flow)
        / (source_y[0] + source_y[1] + scalar_impedance.recip() + load_y[0] + load_y[1]);
    let bank_u = solve_pair(zs, [pressures[0] - p, pressures[1] - p]);
    let expected = [
        bank_u[0],
        bank_u[1],
        (scalar_pressure - p) / scalar_impedance,
        flow,
        -load_y[0] * p,
        -load_y[1] * p,
    ];
    let response = network_drive_source_bank(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &sources,
        &bank,
        &[],
        Some(&exterior),
    )
    .unwrap();
    for (port, flow) in expected.iter().enumerate() {
        close(response.fields.port_pressures[port], p, 3.0);
        close(response.fields.port_inflows[port], *flow, 3.0 / r);
    }
    balance(&response, 10.0 / r);
}

#[test]
fn source_bank_requires_full_passivity_reciprocity_frequency_and_partition() {
    let gas = gas();
    let net = star([0.13, 0.19]);
    let omega = 2531.0;
    let operator = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    let good = bank(
        omega,
        [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
        [C64::ONE; 2],
    );
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Closed,
    }];
    for fault in 0..11 {
        let mut invalid = good.clone();
        match fault {
            0 => invalid.impedance.omega = f64::from_bits(omega.to_bits() + 1),
            1 => invalid.impedance.source_id.clear(),
            2 => {
                invalid.impedance.z[1] = C64::from_re(2.0);
                invalid.impedance.z[2] = C64::from_re(2.0);
            }
            3 => invalid.impedance.z[1] = C64::new(0.0, 0.1),
            4 => invalid.impedance.labels[1] = "left".into(),
            5 => invalid.impedance.labels[1] = "unknown".into(),
            6 => invalid.impedance.labels[1] = "receiver".into(),
            7 => invalid.impedance.z[0] = C64::new(f64::INFINITY, 0.0),
            8 => invalid.open_circuit_pressures_pa[0] = C64::new(0.0, f64::NAN),
            9 => {
                invalid.open_circuit_pressures_pa.pop();
            }
            _ => {
                invalid.impedance.z.pop();
            }
        }
        assert!(
            matches!(
                operator.drive_source_bank(&[], &invalid, &loads, None),
                Err(DuctError::BadParameter { .. })
            ),
            "fault {fault}"
        );
    }
    let duplicate = [PortSource {
        label: "left".into(),
        excitation: PortExcitation::VolumeFlow(C64::ZERO),
    }];
    assert!(matches!(
        operator.drive_source_bank(&duplicate, &good, &[], None),
        Err(DuctError::BadParameter { .. })
    ));
    assert!(matches!(
        operator.drive_source_bank(&[], &good, &[], None),
        Err(DuctError::BadParameter { .. })
    ));
    let mut oversized = good.clone();
    oversized.impedance.labels = (0..33).map(|i| format!("p-{i}")).collect();
    oversized.impedance.z = vec![C64::ZERO; 33 * 33];
    oversized.open_circuit_pressures_pa = vec![C64::ZERO; 33];
    assert!(matches!(
        operator.drive_source_bank(&[], &oversized, &[], None),
        Err(DuctError::BadParameter { .. })
    ));
    // Finite huge and tiny mutual entries must not evade full-matrix admission.
    for scale in [1.0e-300, 1.0e308] {
        let mut active = good.clone();
        active.impedance.z = [0.4, 0.8, 0.8, 0.4]
            .map(|v| C64::from_re(v * scale))
            .to_vec();
        assert!(matches!(
            operator.drive_source_bank(&[], &active, &loads, None),
            Err(DuctError::BadParameter { .. })
        ));
    }
}

#[test]
fn source_bank_does_not_regularize_a_true_loaded_differential_pole() {
    let gas = gas();
    let net = star([0.2; 2]);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    // Passive rank-one impedance has no differential-mode damping. The
    // lossless half-wave graph also has zero differential impedance here.
    let bank = bank(omega, [C64::from_re(zc(&gas)); 4], [C64::ONE, -C64::ONE]);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(zc(&gas))),
    }];
    assert!(matches!(
        network_drive_source_bank(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &[],
            &bank,
            &loads,
            None
        ),
        Err(DuctError::Singular)
    ));
}

#[test]
fn source_bank_passivity_admits_zero_subnormal_and_extreme_passive_operators() {
    // Admission, separately from solvability: a zero or rank-one source law
    // may be valid even when some physical network would leave a free mode.
    for scale in [0.0, f64::from_bits(1), 1.0e-300, f64::MAX] {
        let bank = bank(
            2531.0,
            [scale, -scale, -scale, scale].map(C64::from_re),
            [C64::ZERO; 2],
        );
        source_conditions(&["left", "right"], &[], Some(&bank), 2531.0).unwrap();
    }
    let scale = f64::from_bits(1);
    let active = bank(
        2531.0,
        [scale, 2.0 * scale, 2.0 * scale, scale].map(C64::from_re),
        [C64::ZERO; 2],
    );
    assert!(matches!(
        source_conditions(&["left", "right"], &[], Some(&active), 2531.0),
        Err(DuctError::BadParameter { .. })
    ));
}
