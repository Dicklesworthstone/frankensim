use super::*;

fn two_port_line() -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("input"), port("end")],
        edges: vec![edge("tube", 0, 1, 0.008, 0.2)],
    }
}

fn scale(gas: &GasState) -> (f64, f64) {
    (
        core::f64::consts::PI * gas.sound_speed / 0.2,
        gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008 * 0.008),
    )
}

fn same(a: &NetworkDriveResponse, b: &NetworkDriveResponse) {
    for (x, y) in a
        .port_pressures
        .iter()
        .chain(&a.port_inflows)
        .zip(b.port_pressures.iter().chain(&b.port_inflows))
    {
        assert_eq!(x.re.to_bits(), y.re.to_bits());
        assert_eq!(x.im.to_bits(), y.im.to_bits());
    }
}

// G1: a matched, complex or pressure-release termination closes the actual
// half-wave problem even though its open-circuit impedance does not exist.
#[test]
fn direct_load_closes_open_circuit_pole_without_frequency_shift() {
    let gas = gas();
    let net = two_port_line();
    let (omega, zc) = scale(&gas);
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    for load in [C64::from_re(zc), C64::new(0.3 * zc, -0.7 * zc), C64::ZERO] {
        let loads = [PortTermination {
            label: "end".into(),
            load: NetworkLoad::Impedance(load),
        }];
        let response =
            network_drive(&net, &gas, omega, LossModel::Lossless, "input", &loads).unwrap();
        assert!((response.impedance - load).abs() / zc < 2.0e-12);
        assert!((response.port_pressures[1] + load).abs() / zc < 2.0e-12);
        assert!((response.port_inflows[1] - C64::ONE).abs() < 2.0e-12);
        assert_eq!(response.port_inflows[0], C64::ONE);
    }
    let closed = [PortTermination {
        label: "end".into(),
        load: NetworkLoad::Closed,
    }];
    assert!(matches!(
        network_drive(&net, &gas, omega, LossModel::Lossless, "input", &closed),
        Err(DuctError::Singular)
    ));
}

// G1: full mutual impedance and singular passive loads must be equations,
// not inverses. The rank-one case admits a pressure-release combination.
#[test]
fn direct_mutual_load_keeps_coupled_power_and_nullspace() {
    let gas = gas();
    let (omega, zc) = scale(&gas);
    let net = DuctNetwork {
        nodes: vec![port("input"), port("left"), port("right")],
        edges: vec![
            edge("left", 0, 1, 0.008, 0.2),
            edge("right", 0, 2, 0.008, 0.2),
        ],
    };
    for (entries, impedance, flows) in [
        ([2.0, 0.5, 0.5, 1.0], 0.875, [0.25, 0.75]),
        ([1.0, 2.0, 2.0, 4.0], 0.0, [2.0, -1.0]),
    ] {
        let boundary = CoupledTermination {
            omega,
            labels: vec!["left".into(), "right".into()],
            z: entries.map(|entry| C64::from_re(entry * zc)).to_vec(),
            source_id: "analytic passive two-port load".into(),
        };
        let response = network_drive_coupled(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            "input",
            &[],
            &boundary,
        )
        .unwrap();
        assert!((response.impedance - C64::from_re(impedance * zc)).abs() / zc < 3.0e-12);
        for (i, expected) in flows.iter().enumerate() {
            assert!((response.port_inflows[i + 1] - C64::from_re(*expected)).abs() < 3.0e-12);
            let mut load_pressure = C64::ZERO;
            for j in 0..2 {
                load_pressure =
                    load_pressure + boundary.z[2 * i + j] * response.port_inflows[j + 1];
            }
            assert!((response.port_pressures[i + 1] + load_pressure).abs() / zc < 3.0e-12);
        }
        let load_power: f64 = (1..3)
            .map(|port| {
                -0.5 * (response.port_pressures[port] * response.port_inflows[port].conj()).re
            })
            .sum();
        assert!((0.5 * response.impedance.re - load_power).abs() / zc < 3.0e-12);
    }
}

// G3/G5: the direct entry point preserves admitted old results exactly;
// regional gas records retain their original propagation and exterior.
#[test]
fn direct_load_preserves_legacy_bits_and_regional_medium() {
    use fs_material::gas_species::input::GasRecord;
    let gas = gas();
    let net = two_port_line();
    let (_, zc) = scale(&gas);
    let loads = [PortTermination {
        label: "end".into(),
        load: NetworkLoad::Impedance(C64::new(0.4 * zc, -0.2 * zc)),
    }];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let omega = 2137.0;
        let old = network_port_impedance(&net, &gas, omega, loss)
            .unwrap()
            .drive("input", &loads, &gas)
            .unwrap();
        let direct = network_drive(&net, &gas, omega, loss, "input", &loads).unwrap();
        same(&old, &direct);
    }
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 343.15, gas.pressure).unwrap();
    let media = NetworkMedia::new(
        &net,
        vec![("tube".into(), GasRecord::Legacy(hot))],
        GasRecord::Legacy(gas),
    )
    .unwrap();
    let (omega, _) = scale(&hot);
    let direct =
        network_drive_regional(&net, &media, omega, LossModel::Lossless, "input", &loads).unwrap();
    let expected = match &loads[0].load {
        NetworkLoad::Impedance(z) => *z,
        _ => unreachable!(),
    };
    assert!((direct.impedance - expected).abs() / zc < 3.0e-12);
    let boundary = CoupledTermination {
        omega,
        labels: vec!["end".into()],
        z: vec![expected],
        source_id: "same regional load".into(),
    };
    let coupled = network_drive_coupled_regional(
        &net,
        &media,
        omega,
        LossModel::Lossless,
        "input",
        &[],
        &boundary,
    )
    .unwrap();
    same(&direct, &coupled);
}

#[test]
fn direct_load_refuses_incomplete_active_or_nonpassive_boundaries_at_pole() {
    let gas = gas();
    let net = two_port_line();
    let (omega, zc) = scale(&gas);
    for (active, loads) in [
        (
            "absent",
            vec![PortTermination {
                label: "end".into(),
                load: NetworkLoad::Closed,
            }],
        ),
        ("input", vec![]),
        (
            "input",
            vec![PortTermination {
                label: "input".into(),
                load: NetworkLoad::Closed,
            }],
        ),
        (
            "input",
            vec![PortTermination {
                label: "absent".into(),
                load: NetworkLoad::Closed,
            }],
        ),
        (
            "input",
            vec![PortTermination {
                label: "end".into(),
                load: NetworkLoad::Impedance(C64::from_re(-zc)),
            }],
        ),
    ] {
        assert!(matches!(
            network_drive(&net, &gas, omega, LossModel::Lossless, active, &loads),
            Err(DuctError::BadParameter { .. })
        ));
    }
    let mut boundary = CoupledTermination {
        omega: omega + 1.0,
        labels: vec!["end".into()],
        z: vec![C64::from_re(zc)],
        source_id: "wrong-frequency load".into(),
    };
    assert!(
        network_drive_coupled(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            "input",
            &[],
            &boundary,
        )
        .is_err()
    );
    boundary.omega = omega;
    boundary.z[0] = C64::from_re(-zc);
    assert!(
        network_drive_coupled(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            "input",
            &[],
            &boundary,
        )
        .is_err()
    );
}
