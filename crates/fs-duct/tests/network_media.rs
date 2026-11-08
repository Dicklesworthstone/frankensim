//! G1/G3: gas substitutions act on the original graph and boundary operators.
use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkMedia, NetworkNode,
    NodeBoundary, PortTermination, network_port_impedance, network_port_impedance_regional,
};
use fs_duct::{Duct, HoleState, LossModel, Segment, Termination, WallModel};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_math::{c64::C64, det};

fn record(text: &str) -> GasRecord {
    let mut words = text.split_ascii_whitespace();
    let record = read_gas_record(&mut words, None, None).unwrap();
    assert!(words.next().is_none());
    record
}
fn air() -> GasRecord {
    record("mixture 293.15 100000 mole 2 N2 0.79 O2 0.21")
}
fn port(id: &str) -> NetworkNode {
    NetworkNode {
        id: id.into(),
        boundary: NodeBoundary::Port { label: id.into() },
    }
}
fn edge(id: &str, from: usize, to: usize, radius: f64, length: f64) -> NetworkEdge {
    NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    }
}
fn one_edge() -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("in"), port("out")],
        edges: vec![edge("bore", 0, 1, 0.006, 0.21)],
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    let relative = (a - b).abs() / b.abs().max(1e-30);
    assert!(relative < tolerance, "{a:?} != {b:?}: {relative}");
}
fn bits(values: &[C64]) -> Vec<[u64; 2]> {
    values
        .iter()
        .map(|z| [z.re.to_bits(), z.im.to_bits()])
        .collect()
}

// Independent p/U cylinder boundary algebra, exp(-i omega t). Only the local
// constitutive wave coefficients come from their shared material/ZK owner.
fn cylinder(
    gas: &GasState,
    radius: f64,
    length: f64,
    omega: f64,
    loss: LossModel,
    load: C64,
) -> C64 {
    let wave = fs_duct::segment_wave(gas, radius, omega, loss).unwrap();
    let phase = wave.wavenumber.scale(length);
    let exp = |z: C64| {
        let (s, c) = det::sin_cos(z.re);
        C64::new(c, s).scale(det::exp(-z.im))
    };
    let (forward, backward) = (exp(phase), exp(-phase));
    let cosine = (forward + backward).scale(0.5);
    let sine = (forward - backward) * C64::new(0.0, -0.5);
    let zc = wave.characteristic_impedance;
    (cosine * load - C64::new(0.0, 1.0) * zc * sine)
        / (cosine - C64::new(0.0, 1.0) * sine * load / zc)
}

#[test]
fn g5_uniform_binding_preserves_original_bits_and_complete_records() {
    let net = one_edge();
    let gas = air();
    let media = NetworkMedia::new(&net, vec![("bore".into(), gas.clone())], gas.clone()).unwrap();
    let loads = [PortTermination {
        label: "out".into(),
        load: NetworkLoad::Radiation {
            radius: 0.006,
            termination: Termination::FlangedOpen,
        },
    }];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for hz in [239.0, 631.0, 1289.0] {
            let omega = core::f64::consts::TAU * hz;
            let old = network_port_impedance(&net, gas.state(), omega, loss).unwrap();
            let bound = network_port_impedance_regional(&net, &media, omega, loss).unwrap();
            assert_eq!(bits(&old.z), bits(&bound.z));
            assert_eq!(
                old.condition_estimate.to_bits(),
                bound.condition_estimate.to_bits()
            );
            assert_eq!(old.edge_conditions, bound.edge_conditions);
            assert_eq!(bound.media().unwrap(), &media);
            assert_eq!(bound.media().unwrap().edge_record("bore"), Some(&gas));
            assert_eq!(bound.media().unwrap().exterior_record(), &gas);
            let a = old.drive("in", &loads, gas.state()).unwrap();
            let b = bound.drive_bound("in", &loads).unwrap();
            assert_eq!(bits(&a.port_pressures), bits(&b.port_pressures));
            assert_eq!(bits(&a.port_inflows), bits(&b.port_inflows));
        }
    }
}

#[test]
fn g1_g3_distinct_branch_gases_match_serial_parallel_acoustic_solution() {
    let records = [
        air(),
        record("mixture 310 100000 mole 1 He 1"),
        record("mixture 280 100000 mole 1 CO2 1"),
    ];
    let loads = [C64::new(2e6, -1e5), C64::new(3e6, 2e5)];
    let net = DuctNetwork {
        nodes: vec![
            port("in"),
            NetworkNode {
                id: "junction".into(),
                boundary: NodeBoundary::Junction,
            },
            NetworkNode {
                id: "left-end".into(),
                boundary: NodeBoundary::Termination(NetworkLoad::Impedance(loads[0])),
            },
            NetworkNode {
                id: "right-end".into(),
                boundary: NodeBoundary::Termination(NetworkLoad::Impedance(loads[1])),
            },
        ],
        edges: vec![
            edge("stem", 0, 1, 0.006, 0.11),
            edge("left", 1, 2, 0.004, 0.17),
            edge("right", 1, 3, 0.005, 0.23),
        ],
    };
    let media = NetworkMedia::new(
        &net,
        vec![
            ("right".into(), records[2].clone()),
            ("stem".into(), records[0].clone()),
            ("left".into(), records[1].clone()),
        ],
        air(),
    )
    .unwrap();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for hz in [211.0, 537.0, 913.0] {
            let omega = core::f64::consts::TAU * hz;
            let left = cylinder(records[1].state(), 0.004, 0.17, omega, loss, loads[0]);
            let right = cylinder(records[2].state(), 0.005, 0.23, omega, loss, loads[1]);
            let at_junction = (left.recip() + right.recip()).recip();
            let expected = cylinder(records[0].state(), 0.006, 0.11, omega, loss, at_junction);
            let bound = network_port_impedance_regional(&net, &media, omega, loss).unwrap();
            close(bound.z[0], expected, 3e-12);
            let uniform = network_port_impedance(&net, records[0].state(), omega, loss).unwrap();
            assert!((bound.z[0] - uniform.z[0]).abs() > 1e-3 * bound.z[0].abs());
            assert!(bound.z[0].re >= 0.0, "passive loaded graph");
            let mut reordered = net.clone();
            reordered.edges.reverse();
            let repeat = network_port_impedance_regional(&reordered, &media, omega, loss).unwrap();
            assert_eq!(
                bits(&bound.z),
                bits(&repeat.z),
                "gas follows edge id, not order"
            );
        }
    }
}

#[test]
fn g1_exterior_load_is_independent_of_interior_and_bound_for_deferred_ports() {
    let net = one_edge();
    let interior = air();
    let outside = record("mixture 300 100000 mole 1 He 1");
    let omega = core::f64::consts::TAU * 537.0;
    let media = NetworkMedia::new(
        &net,
        vec![("bore".into(), interior.clone())],
        outside.clone(),
    )
    .unwrap();
    let radiation = NetworkLoad::Radiation {
        radius: 0.006,
        termination: Termination::FlangedOpen,
    };
    let loads = [PortTermination {
        label: "out".into(),
        load: radiation.clone(),
    }];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let result = network_port_impedance_regional(&net, &media, omega, loss).unwrap();
        let driven = result.drive_bound("in", &loads).unwrap();
        // The exterior owner retains the interior higher-mode impedances in
        // its flanged Schur load. Exterior-only rho/c is not that boundary.
        let expected = fs_duct::exterior::input_impedance(
            &net.edges[0].duct,
            interior.state(),
            outside.state(),
            omega,
            loss,
            Termination::FlangedOpen,
            WallModel::Rigid,
        )
        .unwrap()
        .impedance;
        close(driven.impedance, expected, 3e-12);
        // The same boundary directly assembled into the nodal operator agrees.
        let mut terminated = net.clone();
        terminated.nodes[1].boundary = NodeBoundary::Termination(radiation.clone());
        let assembled = network_port_impedance_regional(&terminated, &media, omega, loss).unwrap();
        close(assembled.z[0], driven.impedance, 3e-12);
        assert!(
            result.drive("in", &loads, interior.state()).is_err(),
            "no stale exterior substitution"
        );
        let uniform = network_port_impedance(&net, interior.state(), omega, loss)
            .unwrap()
            .drive("in", &loads, interior.state())
            .unwrap();
        assert!((uniform.impedance - driven.impedance).abs() > 1e-4 * driven.impedance.abs());
    }
}

#[test]
fn g0_modal_radiation_refuses_an_ambiguous_port_plane_and_single_fluid_ls_contrast() {
    let mut net = one_edge();
    net.edges.push(edge("parallel", 0, 1, 0.004, 0.33));
    let interior = air();
    let other = record("mixture 300 100000 mole 1 He 1");
    let media = NetworkMedia::new(
        &net,
        vec![
            ("bore".into(), interior.clone()),
            ("parallel".into(), other.clone()),
        ],
        interior.clone(),
    )
    .unwrap();
    let omega = core::f64::consts::TAU * 537.0;
    let radiation = NetworkLoad::Radiation {
        radius: 0.006,
        termination: Termination::FlangedOpen,
    };
    let response = network_port_impedance_regional(&net, &media, omega, LossModel::Bessel).unwrap();
    assert!(
        response
            .drive_bound(
                "in",
                &[PortTermination {
                    label: "out".into(),
                    load: radiation.clone()
                }]
            )
            .is_err()
    );
    assert!(
        response
            .drive_bound(
                "in",
                &[PortTermination {
                    label: "out".into(),
                    load: NetworkLoad::Impedance(C64::from_re(1e6))
                }]
            )
            .is_ok()
    );
    net.nodes[1].boundary = NodeBoundary::Termination(radiation);
    assert!(network_port_impedance_regional(&net, &media, omega, LossModel::Bessel).is_err());

    let net = one_edge();
    let media = NetworkMedia::new(&net, vec![("bore".into(), interior)], other).unwrap();
    let response = network_port_impedance_regional(&net, &media, omega, LossModel::Bessel).unwrap();
    let exact_ls = NetworkLoad::Radiation {
        radius: 0.006,
        termination: Termination::UnflangedThinWall {
            quadrature: fs_duct::levine_schwinger::LsQuadrature { nodes: 64 },
        },
    };
    assert!(
        response
            .drive_bound(
                "in",
                &[PortTermination {
                    label: "out".into(),
                    load: exact_ls.clone()
                }]
            )
            .is_err()
    );
    let mut terminated = net;
    terminated.nodes[1].boundary = NodeBoundary::Termination(exact_ls);
    assert!(
        network_port_impedance_regional(&terminated, &media, omega, LossModel::Bessel).is_err()
    );
}

#[test]
fn g1_embedded_chimney_and_bore_keep_their_actual_media() {
    let interior = record("mixture 300 100000 mole 1 CO2 1");
    let exterior = air();
    let omega = core::f64::consts::TAU * 731.0;
    for state in [HoleState::Open, HoleState::Closed, HoleState::Vent(0.35)] {
        let mut net = one_edge();
        net.edges[0].corrections = JunctionCorrections::ChainLocal;
        net.edges[0].duct.segments = vec![
            Segment::Cylinder {
                radius: 0.006,
                length: 0.09,
            },
            Segment::ToneHole {
                hole_radius: 0.001,
                chimney_height: 0.004,
                bore_radius: 0.006,
                state,
            },
            Segment::Cylinder {
                radius: 0.006,
                length: 0.12,
            },
        ];
        let media = NetworkMedia::new(
            &net,
            vec![("bore".into(), interior.clone())],
            exterior.clone(),
        )
        .unwrap();
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            let result = network_port_impedance_regional(&net, &media, omega, loss)
                .unwrap()
                .drive_bound(
                    "in",
                    &[PortTermination {
                        label: "out".into(),
                        load: NetworkLoad::PressureRelease,
                    }],
                )
                .unwrap();
            let direct = fs_duct::exterior::input_impedance(
                &net.edges[0].duct,
                interior.state(),
                exterior.state(),
                omega,
                loss,
                Termination::IdealOpen,
                WallModel::Rigid,
            )
            .unwrap();
            close(result.impedance, direct.impedance, 3e-12);
        }
    }
}

#[test]
fn g3_mixed_graph_remains_reciprocal_passive_and_pressure_scaling_is_physical() {
    let mut net = one_edge();
    net.edges.push(edge("parallel", 0, 1, 0.004, 0.33));
    let omega = core::f64::consts::TAU * 641.0;
    let bind = |pressure| {
        NetworkMedia::new(
            &net,
            vec![
                (
                    "bore".into(),
                    record(&format!("mixture 300 {pressure} mole 1 N2 1")),
                ),
                (
                    "parallel".into(),
                    record(&format!("mixture 310 {pressure} mole 1 He 1")),
                ),
            ],
            record(&format!("mixture 293.15 {pressure} mole 1 O2 1")),
        )
        .unwrap()
    };
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let result = network_port_impedance_regional(&net, &bind(100000), omega, loss).unwrap();
        close(result.z[1], result.z[2], 3e-12);
        for current in [
            [C64::ONE, C64::ONE],
            [C64::ONE, -C64::ONE],
            [C64::ONE, C64::new(0.0, 1.0)],
        ] {
            let pressure = [
                result.z[0] * current[0] + result.z[1] * current[1],
                result.z[2] * current[0] + result.z[3] * current[1],
            ];
            let power = (current[0].conj() * pressure[0] + current[1].conj() * pressure[1]).re;
            let scale = pressure[0].abs() + pressure[1].abs();
            assert!(power >= -1e-12 * scale, "nonnegative net incoming power");
        }
    }
    let low =
        network_port_impedance_regional(&net, &bind(50000), omega, LossModel::Lossless).unwrap();
    let high =
        network_port_impedance_regional(&net, &bind(100000), omega, LossModel::Lossless).unwrap();
    for (a, b) in low.z.iter().zip(&high.z) {
        close(a.scale(2.0), *b, 3e-12);
    }
}

#[test]
fn g0_material_binding_refuses_missing_identity_pressure_jumps_and_incoherence() {
    let mut net = one_edge();
    net.edges.push(edge("parallel", 0, 1, 0.004, 0.33));
    let original = air();
    assert!(NetworkMedia::new(&net, vec![], original.clone()).is_err());
    for names in [["bore", "bore"], ["bore", "unknown"]] {
        assert!(
            NetworkMedia::new(
                &net,
                names
                    .into_iter()
                    .map(|id| (id.into(), original.clone()))
                    .collect(),
                original.clone()
            )
            .is_err()
        );
    }
    let low = record("mixture 293.15 50000 mole 1 N2 1");
    assert!(
        NetworkMedia::new(
            &net,
            vec![("bore".into(), original.clone()), ("parallel".into(), low)],
            original.clone()
        )
        .is_err()
    );
    let mut stale = *original.state();
    stale.sound_speed *= 1.2;
    assert!(
        NetworkMedia::new(
            &net,
            vec![
                ("bore".into(), original.clone()),
                ("parallel".into(), GasRecord::Legacy(stale))
            ],
            original
        )
        .is_err()
    );
}

// N1 (MB 15606): an edge may carry the series duct's locally reacting wall.
// The walled one-edge graph is the series `WallModel::Pin` impedance on the
// same duct and load, to the tolerance the rigid comparison above uses; an
// edge without a bound wall keeps the original rigid bits.
#[test]
fn g6_edge_wall_is_the_series_pin_wall_and_unbound_edges_keep_rigid_bits() {
    let net = one_edge();
    let gas = air();
    // A test pin (not a specimen): soft enough that it moves the impedance.
    let pin = fs_phs::WallPin {
        surface_density: 1.5,
        stiffness_per_area: 2.0e5,
        resistance: 40.0,
    };
    let rigid = NetworkMedia::new(&net, vec![("bore".into(), gas.clone())], gas.clone()).unwrap();
    let walled = rigid
        .clone()
        .with_edge_walls(&net, vec![("bore".into(), pin)])
        .unwrap();
    assert_eq!(walled.edge_wall("bore"), Some(&pin));
    assert_eq!(rigid.edge_wall("bore"), None);
    let loads = [PortTermination {
        label: "out".into(),
        load: NetworkLoad::Radiation {
            radius: 0.006,
            termination: Termination::FlangedOpen,
        },
    }];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for hz in [97.0, 239.0, 631.0, 1289.0, 3517.0] {
            let omega = core::f64::consts::TAU * hz;
            let expected = fs_duct::exterior::input_impedance(
                &net.edges[0].duct,
                gas.state(),
                gas.state(),
                omega,
                loss,
                Termination::FlangedOpen,
                WallModel::Pin(&pin),
            )
            .unwrap()
            .impedance;
            let bound = network_port_impedance_regional(&net, &walled, omega, loss)
                .unwrap()
                .drive_bound("in", &loads)
                .unwrap();
            close(bound.impedance, expected, 3e-12);
            let driven = fs_duct::network::network_drive_regional(
                &net, &walled, omega, loss, "in", &loads,
            )
            .unwrap();
            close(driven.impedance, expected, 3e-12);
            let plain = network_port_impedance_regional(&net, &rigid, omega, loss).unwrap();
            let old = network_port_impedance(&net, gas.state(), omega, loss).unwrap();
            assert_eq!(bits(&plain.z), bits(&old.z));
            let stiff = plain.drive_bound("in", &loads).unwrap().impedance;
            assert!((stiff - bound.impedance).abs() > 1e-6 * stiff.abs(), "{hz} Hz");
        }
    }
    let mut two = one_edge();
    two.edges.push(edge("parallel", 0, 1, 0.004, 0.33));
    let media = NetworkMedia::new(
        &two,
        vec![("bore".into(), gas.clone()), ("parallel".into(), gas.clone())],
        gas.clone(),
    )
    .unwrap();
    assert!(media.clone().with_edge_walls(&two, vec![("missing".into(), pin)]).is_err());
    assert!(
        media
            .clone()
            .with_edge_walls(&two, vec![("bore".into(), pin), ("bore".into(), pin)])
            .is_err()
    );
    let once = media.with_edge_walls(&two, vec![("parallel".into(), pin)]).unwrap();
    assert!(once.clone().with_edge_walls(&two, vec![("bore".into(), pin)]).is_err());
    assert_eq!(once.edge_wall("bore"), None);
    assert_eq!(once.edge_wall("parallel"), Some(&pin));
}

// Network cylinder heat is the rigid-wall Bessel line. A walled edge is solved
// on the walled line, so its heat refuses instead of mixing the two lines; the
// same edge without a bound wall keeps its heat.
#[test]
fn g7_cylinder_heat_refuses_a_walled_edge_and_keeps_an_unbound_one() {
    use fs_duct::network::field::{
        NetworkFieldBudget, NetworkFieldProbe, NetworkFieldRequest, NetworkFieldSection,
        PreparedNetworkFields,
    };
    use fs_duct::network::{PortExcitation, PortSource};
    let net = one_edge();
    let gas = air();
    let pin = fs_phs::WallPin {
        surface_density: 1.5,
        stiffness_per_area: 2.0e5,
        resistance: 40.0,
    };
    let rigid = NetworkMedia::new(&net, vec![("bore".into(), gas.clone())], gas.clone()).unwrap();
    let walled = rigid
        .clone()
        .with_edge_walls(&net, vec![("bore".into(), pin)])
        .unwrap();
    let sources = [PortSource {
        label: "in".into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::ONE,
            impedance_pa_s_m3: C64::ZERO,
        },
    }];
    let loads = [PortTermination {
        label: "out".into(),
        load: NetworkLoad::Radiation {
            radius: 0.006,
            termination: Termination::FlangedOpen,
        },
    }];
    let probes: Vec<NetworkFieldProbe> = [0.0, 0.105]
        .into_iter()
        .map(|offset_m| NetworkFieldProbe {
            edge_id: "bore".into(),
            section: NetworkFieldSection::Cylinder { segment: 0, offset_m },
        })
        .collect();
    let request = NetworkFieldRequest {
        omega: core::f64::consts::TAU * 631.0,
        loss: LossModel::Bessel,
        source_clock: "n1-wall-heat",
        sources: &sources,
        source_bank: None,
        loads: &loads,
        coupled_load: None,
        probes: &probes,
        budget: NetworkFieldBudget::default(),
    };
    let plain = PreparedNetworkFields::new_regional(&net, &rigid, &request, || false).unwrap();
    assert!(plain.prepare_cylinder_absorption(&[[0, 1]], 1, || false).is_ok());
    let soft = PreparedNetworkFields::new_regional(&net, &walled, &request, || false).unwrap();
    assert!(soft.prepare_cylinder_absorption(&[[0, 1]], 1, || false).is_err());
}
