use super::source_bank::{analytic_bank, balance, bank, close, solve_pair, star, zc};
use super::*;

fn impedance(r: f64) -> [C64; 4] {
    [
        C64::new(0.7, -0.2),
        C64::new(-0.35, 0.6),
        C64::new(-0.35, 0.6),
        C64::new(0.9, 0.4),
    ]
    .map(|value| value.scale(r))
}

// G1/G3: response columns come from independent physical line elimination,
// and changing only coherent drive reuses the same prepared operator/buffers.
#[test]
fn prepared_source_columns_and_sweeps_match_regional_line_equations() {
    use fs_material::gas_species::input::GasRecord;
    let gas = gas();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 373.15, gas.pressure).unwrap();
    let lengths = [0.11, 0.17];
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
    let r = zc(&gas);
    let z = impedance(r);
    let load = C64::new(0.8, -0.3).scale(r);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    for omega in [713.0, 2301.0, 5107.0] {
        let bank = bank(omega, z, [C64::new(1.2, -0.3), C64::new(-0.4, 0.8)]);
        let prepared = network_prepare_sources_regional(
            &net,
            &media,
            omega,
            LossModel::Lossless,
            &[],
            Some(&bank),
            &loads,
            None,
        )
        .unwrap();
        assert_eq!(prepared.source_port_indices(), &[0, 1]);
        for column in 0..2 {
            let mut drive = [C64::ZERO; 2];
            drive[column] = C64::ONE;
            let (p, u) = analytic_bank([&gas, &hot], lengths, omega, z, drive, load);
            for port in 0..3 {
                close(
                    prepared.pressure_transfer()[port * 2 + column],
                    p[port],
                    3.0,
                );
                close(
                    prepared.flow_transfer()[port * 2 + column],
                    u[port],
                    3.0 / r,
                );
            }
        }
        let mut output = prepared.response_buffer();
        let pressure_storage = output.fields.port_pressures.as_ptr();
        for phase in [0.0_f64, 0.3, 1.1, core::f64::consts::PI] {
            let drive = [C64::new(1.2, -0.3), C64::new(phase.cos(), phase.sin())];
            prepared.drive_into(&drive, &mut output).unwrap();
            let (p, u) = analytic_bank([&gas, &hot], lengths, omega, z, drive, load);
            for port in 0..3 {
                close(output.fields.port_pressures[port], p[port], 3.0);
                close(output.fields.port_inflows[port], u[port], 3.0 / r);
            }
            balance(&output, 10.0 / r);
            assert_eq!(output.fields.port_pressures.as_ptr(), pressure_storage);
        }
    }
}

fn cross(a: [C64; 2], covariance: [C64; 4], b: [C64; 2]) -> C64 {
    let mut value = C64::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            value = value + (a[i] * covariance[i * 2 + j] * b[j].conj()).scale(0.5);
        }
    }
    value
}

// G1: independently contract the analytical line response with a declared
// complex, non-diagonal source second moment, including every cross-port p/U.
#[test]
fn partial_coherence_keeps_full_port_moments_and_source_work() {
    let gas = gas();
    let lengths = [0.11, 0.17];
    let net = star(lengths);
    let omega = 2301.0;
    let r = zc(&gas);
    let z = impedance(r);
    let load = C64::new(0.8, -0.3).scale(r);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load),
    }];
    let bank = bank(omega, z, [C64::ZERO; 2]);
    let prepared = network_prepare_sources(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &[],
        Some(&bank),
        &loads,
        None,
    )
    .unwrap();
    let rho = C64::new(0.3, -0.4);
    let amplitude = 1.7;
    let factors = [
        C64::from_re(amplitude),
        rho.scale(amplitude),
        C64::ZERO,
        C64::from_re(amplitude * (1.0 - rho.norm_sq()).sqrt()),
    ];
    let covariance =
        [C64::ONE, rho.conj(), rho, C64::ONE].map(|value| value.scale(amplitude * amplitude));
    let columns = [
        analytic_bank([&gas; 2], lengths, omega, z, [C64::ONE, C64::ZERO], load),
        analytic_bank([&gas; 2], lengths, omega, z, [C64::ZERO, C64::ONE], load),
    ];
    let p = |port| [columns[0].0[port], columns[1].0[port]];
    let u = |port| [columns[0].1[port], columns[1].1[port]];
    let mut output = prepared.statistics_buffer();
    assert!(!output.is_valid());
    prepared
        .coherence_into(&factors, 2, NetworkCoherenceBudget::default(), &mut output)
        .unwrap();
    let pp = output.pressure_second_moments().unwrap();
    let uu = output.flow_second_moments().unwrap();
    let pu = output.pressure_flow_moments().unwrap();
    for i in 0..3 {
        for j in 0..3 {
            close(pp[i * 3 + j], cross(p(i), covariance, p(j)), 20.0);
            close(uu[i * 3 + j], cross(u(i), covariance, u(j)), 20.0 / (r * r));
            close(pu[i * 3 + j], cross(p(i), covariance, u(j)), 20.0 / r);
            close(pp[i * 3 + j], pp[j * 3 + i].conj(), 20.0);
        }
    }
    let mut internal = 0.0;
    let mut generator = 0.0;
    for i in 0..2 {
        let source_coordinate = if i == 0 {
            [C64::ONE, C64::ZERO]
        } else {
            [C64::ZERO, C64::ONE]
        };
        let expected_generator = cross(source_coordinate, covariance, u(i)).re;
        let expected_internal = (z[i * 2] * cross(u(0), covariance, u(i))
            + z[i * 2 + 1] * cross(u(1), covariance, u(i)))
        .re;
        assert!((output.generator_power_w(i).unwrap() - expected_generator).abs() < 1.0e-10 / r);
        assert!(
            (output.source_internal_power_w(i).unwrap() - expected_internal).abs() < 1.0e-10 / r
        );
        generator += expected_generator;
        internal += expected_internal;
    }
    let work = output.work().unwrap();
    assert!((work.generator_power_w - generator).abs() < 1.0e-10 / r);
    assert!((work.source_dissipation_w - internal).abs() < 1.0e-10 / r);
    assert!(
        (work.generator_power_w - work.source_dissipation_w - work.external_load_power_w).abs()
            < 1.0e-10 / r
    );
    assert!(work.network_power_w.abs() < 1.0e-10 / r);
    assert!(pu[1].im.abs() > 1.0e-3 / r); // complex p/U correlation is retained.
}

// G1: fully opposed generators cancel at the receiver yet dissipate internally.
// Independent sources lie between opposed and in-phase radiation, at the exact
// open-port pole where preparation must retain the scattering boundary system.
#[test]
fn prepared_pole_coherence_and_independent_fast_path_match_common_differential_modes() {
    let gas = gas();
    let net = star([0.2; 2]);
    let r = zc(&gas);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    let bank = bank(
        omega,
        [2.0, -0.5, -0.5, 2.0].map(|value| C64::from_re(value * r)),
        [C64::ZERO; 2],
    );
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let prepared = network_prepare_sources(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &[],
        Some(&bank),
        &loads,
        None,
    )
    .unwrap();
    let mut statistics = prepared.statistics_buffer();
    prepared
        .coherence_into(
            &[C64::ONE, -C64::ONE],
            1,
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    let opposed = statistics.work().unwrap();
    assert!(opposed.external_load_power_w.abs() < 1.0e-23 / r);
    assert!((opposed.source_dissipation_w - 1.0 / (2.5 * r)).abs() < 1.0e-10 / r);
    prepared
        .coherence_into(
            &[C64::ONE; 2],
            1,
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    let together = statistics.work().unwrap();
    assert!((together.external_load_power_w - 2.0 / (3.5 * 3.5 * r)).abs() < 1.0e-10 / r);
    prepared
        .independent_into(
            &[1.0; 2],
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    let independent = statistics.work().unwrap();
    assert!((independent.external_load_power_w - 1.0 / (3.5 * 3.5 * r)).abs() < 1.0e-10 / r);
    assert!(
        (independent.source_dissipation_w
            - 0.5 * (opposed.source_dissipation_w + together.source_dissipation_w))
            .abs()
            < 1.0e-10 / r
    );
    let saved = statistics.clone();
    prepared
        .coherence_into(
            &[C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
            2,
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    assert_eq!(
        statistics.pressure_second_moments(),
        saved.pressure_second_moments()
    );
    assert_eq!(
        statistics.flow_second_moments(),
        saved.flow_second_moments()
    );
    assert_eq!(
        statistics.pressure_flow_moments(),
        saved.pressure_flow_moments()
    );
    assert_eq!(statistics.work(), saved.work());
}

#[test]
fn prepared_mixed_units_retain_generator_and_exterior_coupling_at_pole() {
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
    let zs = impedance(r);
    let bank = bank(omega, zs, [C64::ZERO; 2]);
    let zl = [2.0, 0.5, 0.5, 1.0].map(|value| C64::from_re(value * r));
    let exterior = CoupledTermination {
        omega,
        labels: vec!["mouth-a".into(), "mouth-b".into()],
        z: zl.to_vec(),
        source_id: "analytic coupled load".into(),
    };
    let scalar_impedance = C64::new(0.6, 0.3).scale(r);
    let sources = [
        PortSource {
            label: "scalar".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::ZERO,
                impedance_pa_s_m3: scalar_impedance,
            },
        },
        PortSource {
            label: "flow".into(),
            excitation: PortExcitation::VolumeFlow(C64::ZERO),
        },
    ];
    let prepared = network_prepare_sources(
        &net,
        &gas,
        omega,
        LossModel::Lossless,
        &sources,
        Some(&bank),
        &[],
        Some(&exterior),
    )
    .unwrap();
    assert_eq!(
        prepared
            .source_coordinates()
            .iter()
            .map(|value| value.label.as_str())
            .collect::<Vec<_>>(),
        ["flow", "left", "right", "scalar"]
    );
    assert_eq!(
        prepared.source_coordinates()[0].quantity,
        NetworkSourceQuantity::VolumeFlow
    );
    assert!(
        prepared.source_coordinates()[1..]
            .iter()
            .all(|value| value.quantity == NetworkSourceQuantity::OpenCircuitPressure)
    );
    let drive = [
        C64::new(0.3, 0.2).scale(1.0 / r),
        C64::new(1.1, 0.2),
        C64::new(-0.3, 0.4),
        C64::new(0.7, -0.2),
    ];
    let injected = solve_pair(zs, [drive[1], drive[2]]);
    let source_y = solve_pair(zs, [C64::ONE; 2]);
    let load_y = solve_pair(zl, [C64::ONE; 2]);
    let p = (injected[0] + injected[1] + drive[3] / scalar_impedance + drive[0])
        / (source_y[0] + source_y[1] + scalar_impedance.recip() + load_y[0] + load_y[1]);
    let bank_u = solve_pair(zs, [drive[1] - p, drive[2] - p]);
    let u = [
        bank_u[0],
        bank_u[1],
        (drive[3] - p) / scalar_impedance,
        drive[0],
        -load_y[0] * p,
        -load_y[1] * p,
    ];
    let mut coherent = prepared.response_buffer();
    prepared.drive_into(&drive, &mut coherent).unwrap();
    for (port, flow) in u.iter().enumerate() {
        close(coherent.fields.port_pressures[port], p, 3.0);
        close(coherent.fields.port_inflows[port], *flow, 3.0 / r);
    }
    balance(&coherent, 10.0 / r);
    let mut statistics = prepared.statistics_buffer();
    prepared
        .coherence_into(
            &drive,
            1,
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    for i in 0..6 {
        for j in 0..6 {
            close(
                statistics.pressure_flow_moments().unwrap()[i * 6 + j],
                p.scale(0.5) * u[j].conj(),
                5.0 / r,
            );
        }
    }
    assert!(
        (statistics.work().unwrap().generator_power_w - coherent.generator_power_w).abs()
            < 1.0e-11 / r
    );
}

#[test]
fn prepared_source_order_is_canonical_and_refusal_invalidates_only_statistics() {
    let gas = gas();
    let net = star([0.13, 0.19]);
    let r = zc(&gas);
    let omega = 2531.0;
    let original = bank(omega, impedance(r), [C64::ONE, C64::new(0.3, -0.7)]);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let operator = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    let prepared = operator
        .prepare_sources(&[], Some(&original), &loads, None)
        .unwrap();
    let mut reversed = original.clone();
    reversed.impedance.labels.reverse();
    reversed.impedance.z.reverse();
    reversed.open_circuit_pressures_pa.reverse();
    let permuted = operator
        .prepare_sources(&[], Some(&reversed), &loads, None)
        .unwrap();
    assert_eq!(prepared.pressure_transfer(), permuted.pressure_transfer());
    assert_eq!(prepared.flow_transfer(), permuted.flow_transfer());
    assert_eq!(
        prepared.default_excitations(),
        permuted.default_excitations()
    );
    let mut coherent = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut coherent)
        .unwrap();
    let original_pressure = coherent.fields.port_pressures.clone();
    for invalid in [
        vec![],
        vec![C64::ONE],
        vec![C64::new(f64::NAN, 0.0); 2],
        vec![C64::from_re(f64::MAX); 2],
    ] {
        assert!(prepared.drive_into(&invalid, &mut coherent).is_err());
        assert_eq!(coherent.fields.port_pressures, original_pressure);
    }
    let mut statistics = prepared.statistics_buffer();
    for (rank, budget) in [
        (0, NetworkCoherenceBudget::default()),
        (33, NetworkCoherenceBudget::default()),
        (
            1,
            NetworkCoherenceBudget {
                maximum_rank: 0,
                maximum_complex_products: usize::MAX,
            },
        ),
        (
            1,
            NetworkCoherenceBudget {
                maximum_rank: 32,
                maximum_complex_products: 0,
            },
        ),
    ] {
        prepared
            .coherence_into(
                prepared.default_excitations(),
                1,
                NetworkCoherenceBudget::default(),
                &mut statistics,
            )
            .unwrap();
        assert!(
            prepared
                .coherence_into(
                    prepared.default_excitations(),
                    rank,
                    budget,
                    &mut statistics
                )
                .is_err()
        );
        assert!(!statistics.is_valid());
        assert!(statistics.work().is_none());
        assert!(statistics.pressure_flow_moments().is_none());
    }
    for invalid in [[-1.0, 1.0], [f64::INFINITY, 1.0], [f64::NAN, 1.0]] {
        assert!(
            prepared
                .independent_into(&invalid, NetworkCoherenceBudget::default(), &mut statistics)
                .is_err()
        );
        assert!(statistics.pressure_second_moments().is_none());
    }
    assert!(
        prepared
            .coherence_into(
                &[C64::ONE],
                1,
                NetworkCoherenceBudget::default(),
                &mut statistics
            )
            .is_err()
    );
    assert!(
        prepared
            .coherence_into(
                &[C64::from_re(f64::MAX); 2],
                1,
                NetworkCoherenceBudget::default(),
                &mut statistics
            )
            .is_err()
    );
    assert!(statistics.work().is_none());
    prepared
        .independent_into(
            &[0.0; 2],
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )
        .unwrap();
    assert_eq!(statistics.work(), Some(NetworkSourceWork::default()));
}

// G1: reactive coupling exchanges real work between ports, while only the
// Hermitian source operator contributes to total internal dissipation.
#[test]
fn prepared_statistics_keep_signed_reactive_exchange_and_small_real_loss() {
    let gas = gas();
    let net = star([0.2; 2]);
    let r = zc(&gas);
    let omega = core::f64::consts::PI * gas.sound_speed / 0.2;
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let flows = [C64::from_re(1.0 / r), C64::new(0.0, 1.0 / r)];
    let pressure = C64::new(1.0, 1.0);
    for loss in [0.0, 1.0e-9] {
        let z = [
            C64::new(loss, 0.3),
            C64::new(0.0, 2.0),
            C64::new(0.0, 2.0),
            C64::new(loss, -0.7),
        ]
        .map(|v| v.scale(r));
        let drive = [
            pressure + z[0] * flows[0] + z[1] * flows[1],
            pressure + z[2] * flows[0] + z[3] * flows[1],
        ];
        let bank = bank(omega, z, drive);
        let prepared = network_prepare_sources(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &[],
            Some(&bank),
            &loads,
            None,
        )
        .unwrap();
        let mut statistics = prepared.statistics_buffer();
        prepared
            .coherence_into(
                &drive,
                1,
                NetworkCoherenceBudget::default(),
                &mut statistics,
            )
            .unwrap();
        assert!(statistics.generator_power_w(0).unwrap() < 0.0);
        assert!(statistics.source_internal_power_w(0).unwrap() < 0.0);
        assert!(statistics.source_internal_power_w(1).unwrap() > 0.0);
        let work = statistics.work().unwrap();
        assert!((work.external_load_power_w - 1.0 / r).abs() < 1.0e-10 / r);
        if loss == 0.0 {
            assert_eq!(work.source_dissipation_w, 0.0);
        } else {
            assert!((work.source_dissipation_w - loss / r).abs() < 1.0e-9 * loss / r);
        }
    }
}
