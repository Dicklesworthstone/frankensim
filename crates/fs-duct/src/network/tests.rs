use super::*;
use fs_material::gas::GasSpec;
use fs_math::det;
use std::time::Instant;

mod porous;
mod loaded;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn duct(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn node(id: &str, boundary: NodeBoundary) -> NetworkNode {
    NetworkNode {
        id: id.into(),
        boundary,
    }
}
fn port(id: &str) -> NetworkNode {
    node(id, NodeBoundary::Port { label: id.into() })
}
fn edge(id: &str, from: usize, to: usize, radius: f64, length: f64) -> NetworkEdge {
    NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: duct(radius, length),
        corrections: JunctionCorrections::Disabled,
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1.0e-300)
}
fn row(
    case: &str,
    net: &DuctNetwork,
    response: &NetworkResponse,
    error: f64,
    tol: f64,
    rung: usize,
    start: Instant,
) {
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"{case}\",\"inputs\":{{\"provenance\":\"analytic SI cylinder/graph fixtures\"}},\"axis\":\"frequency_or_case\",\"rung\":{rung},\"nodes\":{},\"edges\":{},\"cond\":{},\"observed\":{},\"oracle\":\"independent closed-form or invariant\",\"err\":{error},\"rel_err\":{error},\"tol\":{tol},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        net.nodes.len(),
        net.edges.len(),
        response.condition_estimate,
        response.omega,
        if error <= tol { "pass" } else { "fail" },
        start.elapsed().as_micros()
    );
}
fn cylinder_input(gas: &GasState, r: f64, l: f64, w: f64, loss: LossModel, load: C64) -> C64 {
    let wave = crate::segment_wave(gas, r, w, loss).unwrap();
    let phase = wave.wavenumber.scale(l);
    let exp = |z: C64| {
        let (s, c) = det::sin_cos(z.re);
        C64::new(c, s).scale(det::exp(-z.im))
    };
    let ep = exp(phase);
    let em = exp(-phase);
    let cos = (ep + em).scale(0.5);
    let sin = (ep - em) * C64::new(0.0, -0.5);
    let z = wave.characteristic_impedance;
    (cos * load - C64::new(0.0, 1.0) * z * sin) / (cos - C64::new(0.0, 1.0) * sin * load / z)
}

#[test]
fn network_t_branch_matches_parallel_admittance() {
    let gas = gas();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for (rung, f) in [211.0, 537.0, 913.0].into_iter().enumerate() {
            let start = Instant::now();
            let w = core::f64::consts::TAU * f;
            let mut net = DuctNetwork {
                nodes: vec![port("input")],
                edges: Vec::new(),
            };
            let mut sum = C64::ZERO;
            for (i, (r, l)) in [(0.004, 0.11), (0.005, 0.17), (0.006, 0.23)]
                .into_iter()
                .enumerate()
            {
                let load = C64::from_re(
                    0.7 * gas.density * gas.sound_speed / (core::f64::consts::PI * r * r),
                );
                net.nodes.push(node(
                    &format!("end-{i}"),
                    NodeBoundary::Termination(NetworkLoad::Impedance(load)),
                ));
                net.edges.push(edge(&format!("arm-{i}"), 0, i + 1, r, l));
                sum = sum + cylinder_input(&gas, r, l, w, loss, load).recip();
            }
            let response = network_port_impedance(&net, &gas, w, loss).unwrap();
            let error = relative(response.z[0], sum.recip());
            row(
                "network-t-branch",
                &net,
                &response,
                error,
                1.0e-12,
                rung,
                start,
            );
            assert!(error <= 1.0e-12, "parallel admittance {loss:?}: {error}");
        }
    }
}

// Matched pressure-wave scattering from the measured two-port impedance.
fn transmission(z: &[C64], zc: f64) -> C64 {
    let zc = C64::from_re(zc);
    (z[2] * zc).scale(2.0) / ((z[0] + zc) * (z[3] + zc) - z[1] * z[2])
}

#[test]
fn network_herschel_quincke_tl() {
    // Selamet, Dickey & Novak, JASA 96,3177 (1994), Eq.(8), not 1997:
    // https://mae.osu.edu/sites/default/files/2021-12/J11_0.pdf
    // Re-derived from p continuity and volume balance at the two ends:
    // normalized Y=i[[a,-b],[-b,a]], a=sum(Z0 cot(kLj)/Zj),
    // b=sum(Z0 csc(kLj)/Zj). |S21|²=4b²/[(1-a²+b²)²+4a²].
    let gas = gas();
    let r0 = 0.02;
    let r = r0 / det::sqrt(2.0);
    let lengths = [0.27, 0.53];
    let net = DuctNetwork {
        nodes: vec![port("in"), port("out")],
        edges: vec![
            edge("short", 0, 1, r, lengths[0]),
            edge("long", 0, 1, r, lengths[1]),
        ],
    };
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * r0 * r0);
    let zj = gas.density * gas.sound_speed / (core::f64::consts::PI * r * r);
    for rung in 0..200 {
        let start = Instant::now();
        // Include finite TL notches (matched in-phase passbands), plus a broad
        // frequency sweep. Infinite attenuation poles are checked separately.
        let f = if rung == 70 {
            gas.sound_speed / (lengths[1] - lengths[0])
        } else if rung == 140 {
            2.0 * gas.sound_speed / (lengths[1] - lengths[0])
        } else {
            97.3 + 13.173 * rung as f64
        };
        let w = core::f64::consts::TAU * f;
        let (s1, c1) = det::sin_cos(w * lengths[0] / gas.sound_speed);
        let (s2, c2) = det::sin_cos(w * lengths[1] / gas.sound_speed);
        let alpha = z0 / zj;
        let an = alpha * (c1 * s2 + c2 * s1);
        let bn = alpha * (s1 + s2);
        // Clear sine denominators analytically, including 1-cos²=sin²,
        // to retain accuracy near an individual branch's half-wave pole.
        let real = (1.0 + 2.0 * alpha * alpha) * s1 * s2 + 2.0 * alpha * alpha * (1.0 - c1 * c2);
        let analytic = 10.0 / core::f64::consts::LN_10
            * det::ln((real * real + 4.0 * an * an) / (4.0 * bn * bn));
        let response = network_port_impedance(&net, &gas, w, LossModel::Lossless).unwrap();
        let observed =
            -10.0 / core::f64::consts::LN_10 * det::ln(transmission(&response.z, z0).norm_sq());
        let error = (observed - analytic).abs();
        row(
            "network-herschel-quincke-tl-db",
            &net,
            &response,
            error,
            1.0e-10,
            rung,
            start,
        );
        assert!(
            error <= 1.0e-10,
            "HQ f={f}, TL={observed} oracle={analytic}, error={error}"
        );
    }
    let w = core::f64::consts::TAU * gas.sound_speed / (2.0 * (lengths[1] - lengths[0]));
    let start = Instant::now();
    let response = network_port_impedance(&net, &gas, w, LossModel::Lossless).unwrap();
    let zero = transmission(&response.z, z0).abs();
    row(
        "network-hq-infinite-tl-zero-transmission",
        &net,
        &response,
        zero,
        1.0e-10,
        200,
        start,
    );
    assert!(
        zero <= 1.0e-10,
        "infinite-TL pole must have vanishing transmission, not a clamped dB value"
    );
}

#[test]
fn network_single_branch_equals_tonehole() {
    let gas = gas();
    let r = 0.01;
    let b = 0.003;
    let chimney = 0.12;
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for state in [HoleState::Open, HoleState::Closed, HoleState::Vent(0.4)] {
            for (rung, f) in [121.0, 347.0, 713.0, 1231.0].into_iter().enumerate() {
                let start = Instant::now();
                let w = core::f64::consts::TAU * f;
                let load =
                    C64::from_re(gas.density * gas.sound_speed / (core::f64::consts::PI * r * r));
                let branch_load = match state {
                    HoleState::Open => NetworkLoad::Radiation {
                        radius: b,
                        termination: Termination::FlangedOpen,
                    },
                    HoleState::Closed => NetworkLoad::Closed,
                    HoleState::Vent(sigma) => {
                        // Vent is two parallel physical chimney paths whose
                        // admittances are weighted at the ideal junction below.
                        let open = crate::input_impedance(
                            &duct(b, chimney),
                            &gas,
                            w,
                            loss,
                            Termination::FlangedOpen,
                        )
                        .unwrap()
                        .impedance;
                        let closed = crate::input_impedance(
                            &duct(b, chimney),
                            &gas,
                            w,
                            loss,
                            Termination::Closed,
                        )
                        .unwrap()
                        .impedance;
                        NetworkLoad::Impedance(
                            (open.recip().scale(sigma) + closed.recip().scale(1.0 - sigma)).recip(),
                        )
                    }
                };
                let mut net = DuctNetwork {
                    nodes: vec![
                        port("in"),
                        node("junction", NodeBoundary::Junction),
                        node(
                            "end",
                            NodeBoundary::Termination(NetworkLoad::Impedance(load)),
                        ),
                        node("branch", NodeBoundary::Termination(branch_load)),
                    ],
                    edges: vec![edge("left", 0, 1, r, 0.16), edge("right", 1, 2, r, 0.19)],
                };
                if matches!(state, HoleState::Vent(_)) {
                    // The weighted branch impedance is already referred to the
                    // junction plane, so a zero-length edge is not introduced.
                    let z = if let NodeBoundary::Termination(NetworkLoad::Impedance(z)) =
                        &net.nodes[3].boundary
                    {
                        *z
                    } else {
                        unreachable!()
                    };
                    net.nodes.pop();
                    net.nodes[1].boundary = NodeBoundary::Termination(NetworkLoad::Impedance(z));
                } else {
                    net.edges.push(edge("chimney", 1, 3, b, chimney));
                }
                let chain = Duct {
                    segments: vec![
                        Segment::Cylinder {
                            radius: r,
                            length: 0.16,
                        },
                        Segment::ToneHole {
                            hole_radius: b,
                            chimney_height: chimney,
                            bore_radius: r,
                            state,
                        },
                        Segment::Cylinder {
                            radius: r,
                            length: 0.19,
                        },
                    ],
                };
                let m =
                    network_chain_two_port(&chain, &gas, w, loss, JunctionCorrections::Disabled)
                        .unwrap()
                        .abcd;
                let oracle = (m[0] * load + m[1]) / (m[2] * load + m[3]);
                let response = network_port_impedance(&net, &gas, w, loss).unwrap();
                let error = relative(response.z[0], oracle);
                row(
                    "network-single-branch-tonehole",
                    &net,
                    &response,
                    error,
                    1.0e-10,
                    rung,
                    start,
                );
                assert!(error <= 1.0e-10, "{loss:?} {state:?}: {error}");
            }
        }
    }
}

fn three_port() -> DuctNetwork {
    DuctNetwork {
        nodes: vec![
            node("junction", NodeBoundary::Junction),
            port("a"),
            port("b"),
            port("c"),
        ],
        edges: vec![
            edge("a-edge", 0, 1, 0.004, 0.13),
            edge("b-edge", 0, 2, 0.005, 0.17),
            edge("c-edge", 0, 3, 0.006, 0.23),
        ],
    }
}

#[test]
fn network_reciprocity_passivity() {
    let net = three_port();
    let gas = gas();
    let mut rng = fs_rand::StreamKey {
        seed: 0x79_2026_1016,
        kernel: 79,
        tile: 0,
    }
    .stream();
    for rung in 0..100 {
        let start = Instant::now();
        let w = core::f64::consts::TAU * (100.0 + 2900.0 * rng.next_f64());
        let response = network_port_impedance(&net, &gas, w, LossModel::Bessel).unwrap();
        let z = &response.z;
        let norm = z.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut err = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                err = err.max((z[3 * i + j] - z[3 * j + i]).abs() / norm);
            }
        }
        // All seven principal minors of a 3x3 real symmetric Hermitian part
        // are nonnegative iff it is PSD; positive diagonals alone do not prove it.
        let mut h = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                h[3 * i + j] = 0.5 * (z[3 * i + j].re + z[3 * j + i].re) / norm;
            }
        }
        let mut minors = vec![h[0], h[4], h[8]];
        for i in 0..3 {
            for j in i + 1..3 {
                minors.push(h[3 * i + i] * h[3 * j + j] - h[3 * i + j] * h[3 * j + i]);
            }
        }
        minors.push(
            h[0] * (h[4] * h[8] - h[5] * h[7]) - h[1] * (h[3] * h[8] - h[5] * h[6])
                + h[2] * (h[3] * h[7] - h[4] * h[6]),
        );
        let negative = minors.into_iter().fold(0.0_f64, f64::min);
        err = err.max(-negative);
        row(
            "network-reciprocity-passivity",
            &net,
            &response,
            err,
            1.0e-12,
            rung,
            start,
        );
        assert!(err <= 1.0e-12, "nonpassive/nonreciprocal Z: {err}");
    }
}

#[test]
fn network_node_order_invariance() {
    let net = three_port();
    let gas = gas();
    let w = 1437.0;
    let a = network_port_impedance(&net, &gas, w, LossModel::Bessel).unwrap();
    let permutation = [2, 0, 3, 1];
    let mut reverse = [0; 4];
    for (i, &old) in permutation.iter().enumerate() {
        reverse[old] = i;
    }
    let mut moved = net.clone();
    moved.nodes = permutation.iter().map(|&i| net.nodes[i].clone()).collect();
    moved.edges.reverse();
    for e in &mut moved.edges {
        e.from = reverse[e.from];
        e.to = reverse[e.to];
    }
    let b = network_port_impedance(&moved, &gas, w, LossModel::Bessel).unwrap();
    for (i, p) in a.ports.iter().enumerate() {
        for (j, q) in a.ports.iter().enumerate() {
            let pi = b.ports.iter().position(|x| x == p).unwrap();
            let qi = b.ports.iter().position(|x| x == q).unwrap();
            assert_eq!(
                (a.z[3 * i + j].re.to_bits(), a.z[3 * i + j].im.to_bits()),
                (b.z[3 * pi + qi].re.to_bits(), b.z[3 * pi + qi].im.to_bits())
            );
        }
    }
}

#[test]
fn network_bitwise_repeat() {
    let net = three_port();
    let gas = gas();
    let a = network_port_impedance(&net, &gas, 1793.0, LossModel::Bessel).unwrap();
    let b = network_port_impedance(&net, &gas, 1793.0, LossModel::Bessel).unwrap();
    let bits = |z: &[C64]| {
        z.iter()
            .map(|v| (v.re.to_bits(), v.im.to_bits()))
            .collect::<Vec<_>>()
    };
    assert_eq!(bits(&a.z), bits(&b.z));
}

#[test]
fn network_refusals() {
    let gas = gas();
    let mut net = DuctNetwork {
        nodes: vec![port("p")],
        edges: Vec::new(),
    };
    for i in 1..MAX_NETWORK_NODES {
        net.nodes.push(node(
            &format!("n-{i:03}"),
            if i == MAX_NETWORK_NODES - 1 {
                NodeBoundary::Termination(NetworkLoad::Impedance(C64::from_re(1.0e6)))
            } else {
                NodeBoundary::Junction
            },
        ));
        net.edges
            .push(edge(&format!("e-{i:03}"), i - 1, i, 0.008, 0.001));
    }
    network_port_impedance(&net, &gas, 2131.0, LossModel::Bessel).unwrap();
    let mut over = net.clone();
    over.nodes.push(node("over-cap", NodeBoundary::Junction));
    assert!(over.validate().is_err());
    let mut disconnected = three_port();
    disconnected
        .nodes
        .push(node("island", NodeBoundary::Junction));
    assert!(disconnected.validate().is_err());
    let mut zero = three_port();
    zero.edges[0].duct = duct(0.004, 0.0);
    assert!(zero.validate().is_err());
    let mut invalid = three_port();
    invalid.edges[0].to = usize::MAX;
    assert!(invalid.validate().is_err());
    let mut no_ports = three_port();
    for n in &mut no_ports.nodes {
        n.boundary = NodeBoundary::Junction;
    }
    assert!(no_ports.validate().is_err());
    let mut duplicate = three_port();
    duplicate.nodes[1].id = duplicate.nodes[0].id.clone();
    assert!(duplicate.validate().is_err());
    let mut parallel = DuctNetwork {
        nodes: vec![port("a"), port("b")],
        edges: (0..MAX_NETWORK_EDGES)
            .map(|i| edge(&format!("e-{i}"), 0, 1, 0.008, 0.03))
            .collect(),
    };
    parallel.validate().unwrap();
    parallel.edges.push(edge("over", 0, 1, 0.008, 0.03));
    assert!(parallel.validate().is_err());
    let mut segments = DuctNetwork {
        nodes: vec![port("a"), port("b")],
        edges: vec![edge("long", 0, 1, 0.008, 0.001)],
    };
    segments.edges[0].duct.segments = vec![
        Segment::Cylinder {
            radius: 0.008,
            length: 0.001
        };
        MAX_NETWORK_SEGMENTS
    ];
    segments.validate().unwrap();
    segments.edges[0].duct.segments.push(Segment::Cylinder {
        radius: 0.008,
        length: 0.001,
    });
    assert!(segments.validate().is_err());
    let mut ports = DuctNetwork {
        nodes: (0..MAX_NETWORK_PORTS)
            .map(|i| port(&format!("p{i}")))
            .collect(),
        edges: (1..MAX_NETWORK_PORTS)
            .map(|i| edge(&format!("e{i}"), 0, i, 0.008, 0.02))
            .collect(),
    };
    ports.validate().unwrap();
    ports.nodes.push(port("over"));
    ports
        .edges
        .push(edge("over", 0, MAX_NETWORK_PORTS, 0.008, 0.02));
    assert!(ports.validate().is_err());
    assert!(matches!(
        inverse(&[C64::ONE, C64::ONE, C64::ONE, C64::ONE], 2),
        Err(DuctError::Singular)
    ));
    assert!(inverse(&[C64::ONE, C64::ZERO, C64::ZERO, C64::from_re(1.0e-11)], 2).is_err());
    assert!(network_port_impedance(&three_port(), &gas, f64::NAN, LossModel::Lossless).is_err());
}

#[test]
fn network_explicit_port_loads_and_exact_ls_owner() {
    let gas = gas();
    let radius = 0.008;
    let length = 0.21;
    let w = 1649.0;
    let net = DuctNetwork {
        nodes: vec![port("source"), port("mouth")],
        edges: vec![edge("pipe", 0, 1, radius, length)],
    };
    let response = network_port_impedance(&net, &gas, w, LossModel::Bessel).unwrap();
    let termination = Termination::UnflangedThinWall {
        quadrature: crate::levine_schwinger::LsQuadrature { nodes: 64 },
    };
    let loads = vec![PortTermination {
        label: "mouth".into(),
        load: NetworkLoad::Radiation {
            radius,
            termination,
        },
    }];
    let observed = response.driving_point("source", &loads, &gas).unwrap();
    let oracle = crate::input_impedance(
        &duct(radius, length),
        &gas,
        w,
        LossModel::Bessel,
        termination,
    )
    .unwrap()
    .impedance;
    assert!(relative(observed, oracle) <= 1.0e-12);
    assert!(response.driving_point("source", &[], &gas).is_err());
    assert!(response.driving_point("unknown", &loads, &gas).is_err());
    let other = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, 101_325.0).unwrap();
    assert!(response.driving_point("source", &loads, &other).is_err());
    let bad_load = vec![PortTermination {
        label: "mouth".into(),
        load: NetworkLoad::Impedance(C64::from_re(-1.0)),
    }];
    assert!(response.driving_point("source", &bad_load, &gas).is_err());
    let table = TabulatedLoad::try_new(
        vec![w, w * 2.0],
        vec![C64::from_re(1.0e6); 2],
        "analytic flat passive load",
        0.0,
    )
    .unwrap();
    let tabulated = [PortTermination {
        label: "mouth".into(),
        load: NetworkLoad::Tabulated(table),
    }];
    let observed = response.driving_point("source", &tabulated, &gas).unwrap();
    let oracle = crate::input_impedance_load(
        &duct(radius, length),
        &gas,
        w,
        LossModel::Bessel,
        C64::from_re(1.0e6),
        None,
    )
    .unwrap()
    .impedance;
    assert!(relative(observed, oracle) <= 1.0e-12);
    let out = network_port_impedance(&net, &gas, w * 0.5, LossModel::Bessel).unwrap();
    assert!(matches!(
        out.driving_point("source", &tabulated, &gas),
        Err(DuctError::OutOfTable { .. })
    ));
}

#[test]
fn network_supplied_junction_multiport() {
    let gas = gas();
    let w = 1731.0;
    let resistance = 3.0e5;
    let conductance = C64::from_re(1.0 / resistance);
    let feature = JunctionMultiport {
        omega: w,
        source_id: "analytic reciprocal resistor; diagnostic feature seam".into(),
        terminals: vec![
            NetworkTerminal {
                edge_id: "left".into(),
                end: EdgeEnd::To,
            },
            NetworkTerminal {
                edge_id: "right".into(),
                end: EdgeEnd::From,
            },
        ],
        y: vec![conductance, -conductance, -conductance, conductance],
    };
    let mut net = DuctNetwork {
        nodes: vec![
            port("in"),
            node("feature", NodeBoundary::JunctionMultiport(feature.clone())),
            port("out"),
        ],
        edges: vec![
            edge("left", 0, 1, 0.008, 0.11),
            edge("right", 1, 2, 0.008, 0.17),
        ],
    };
    let response = network_port_impedance(&net, &gas, w, LossModel::Lossless).unwrap();
    let load = C64::from_re(1.0e6);
    let loads = [PortTermination {
        label: "out".into(),
        load: NetworkLoad::Impedance(load),
    }];
    let m1 = network_chain_two_port(
        &net.edges[0].duct,
        &gas,
        w,
        LossModel::Lossless,
        JunctionCorrections::Disabled,
    )
    .unwrap()
    .abcd;
    let m2 = network_chain_two_port(
        &net.edges[1].duct,
        &gas,
        w,
        LossModel::Lossless,
        JunctionCorrections::Disabled,
    )
    .unwrap()
    .abcd;
    let m = crate::mul2(
        crate::mul2(
            m1,
            [C64::ONE, C64::from_re(resistance), C64::ZERO, C64::ONE],
        ),
        m2,
    );
    let oracle = (m[0] * load + m[1]) / (m[2] * load + m[3]);
    assert!(relative(response.driving_point("in", &loads, &gas).unwrap(), oracle) <= 1.0e-12);
    assert!(network_port_impedance(&net, &gas, w + 1.0, LossModel::Lossless).is_err());
    if let NodeBoundary::JunctionMultiport(f) = &mut net.nodes[1].boundary {
        f.y[0] = -conductance;
    }
    assert!(net.validate().is_err());
    let mut wrong = feature;
    wrong.terminals[1].end = EdgeEnd::To;
    net.nodes[1].boundary = NodeBoundary::JunctionMultiport(wrong);
    assert!(net.validate().is_err());
}

#[test]
fn network_self_loop_matches_two_closed_half_arms() {
    let gas = gas();
    let r = 0.007;
    let length = 0.43;
    let w = 2317.0;
    let net = DuctNetwork {
        nodes: vec![port("ring-tap")],
        edges: vec![edge("ring", 0, 0, r, length)],
    };
    let response = network_port_impedance(&net, &gas, w, LossModel::Bessel).unwrap();
    let half = crate::input_impedance(
        &duct(r, length / 2.0),
        &gas,
        w,
        LossModel::Bessel,
        Termination::Closed,
    )
    .unwrap()
    .impedance;
    assert!(relative(response.z[0], half.scale(0.5)) <= 1.0e-12);
}

#[test]
fn network_expanded_unknown_cap_and_singular_edge_refuse() {
    let mut edges = vec![edge("lead", 0, 1, 0.008, 0.02)];
    for i in 0..127 {
        edges.push(edge(&format!("loop-{i:03}"), 1, 1, 0.008, 0.02));
    }
    let mut terminals = Vec::new();
    for e in &edges {
        if e.from == 1 {
            terminals.push(NetworkTerminal {
                edge_id: e.id.clone(),
                end: EdgeEnd::From,
            });
        }
        if e.to == 1 {
            terminals.push(NetworkTerminal {
                edge_id: e.id.clone(),
                end: EdgeEnd::To,
            });
        }
    }
    assert_eq!(terminals.len() + 1, MAX_NETWORK_UNKNOWNS);
    let d = terminals.len();
    let mut y = vec![C64::ZERO; d * d];
    for i in 0..d {
        y[i * d + i] = C64::from_re(1.0e-6);
    }
    let mut net = DuctNetwork {
        nodes: vec![
            port("p"),
            node(
                "feature",
                NodeBoundary::JunctionMultiport(JunctionMultiport {
                    omega: 1000.0,
                    source_id: "analytic diagonal conductance cap fixture".into(),
                    terminals,
                    y,
                }),
            ),
        ],
        edges,
    };
    net.validate().unwrap();
    net.nodes.push(node("over", NodeBoundary::Junction));
    net.edges.push(edge("over", 0, 2, 0.008, 0.02));
    assert!(net.validate().is_err());
    let singular = DuctNetwork {
        nodes: vec![port("a"), port("b")],
        edges: vec![edge("tiny-phase", 0, 1, 0.008, 0.01)],
    };
    assert!(matches!(
        network_port_impedance(&singular, &gas(), 1.0e-9, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
}

// G1: sin(kL)=0 is a pole of the edge's pressure-to-flow chart, not of
// a lossless line terminated in a finite dissipative load. Never add loss,
// perturb the frequency or subdivide the geometry to make this solve exist.
#[test]
fn network_retains_scattering_at_half_wave_and_low_frequency() {
    let gas = gas();
    let radius = 0.008;
    let length = 0.2;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let net = DuctNetwork {
        nodes: vec![
            port("input"),
            node(
                "end",
                NodeBoundary::Termination(NetworkLoad::Impedance(C64::from_re(zc))),
            ),
        ],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let half_wave = core::f64::consts::PI * gas.sound_speed / length;
    for omega in [1.0e-9, half_wave, 2.0 * half_wave, 3.0 * half_wave] {
        assert!(
            network_chain_admittance(
                &net.edges[0].duct,
                &gas,
                omega,
                LossModel::Lossless,
                JunctionCorrections::Disabled,
            )
            .is_err()
        );
        let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
        assert_eq!(response.wave_constraint_edges, 1);
        assert!(
            relative(response.z[0], C64::from_re(zc)) < 2.0e-12,
            "omega={omega}: {:?}",
            response.z
        );
        assert!(response.condition_estimate < 100.0);
    }
}

// G3: equivalent subdivisions must retain the same complex resonance loading,
// including when only one subdivision places a component exactly at a pole.
#[test]
fn network_half_wave_subdivision_preserves_loaded_response() {
    let gas = gas();
    let radius = 0.008;
    let length = 0.6;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let load = C64::new(0.3 * zc, -0.7 * zc);
    for detuning in [-1.0e-6, 0.0, 1.0e-6] {
        let omega = 3.0 * core::f64::consts::PI * gas.sound_speed / length * (1.0 + detuning);
        let expected = cylinder_input(&gas, radius, length, omega, LossModel::Lossless, load);
        for pieces in [1, 2, 3, 7] {
            let mut net = DuctNetwork {
                nodes: vec![port("input")],
                edges: Vec::new(),
            };
            for i in 0..pieces {
                let boundary = if i + 1 == pieces {
                    NodeBoundary::Termination(NetworkLoad::Impedance(load))
                } else {
                    NodeBoundary::Junction
                };
                net.nodes.push(node(&format!("node-{i}"), boundary));
                net.edges.push(edge(
                    &format!("edge-{i}"),
                    i,
                    i + 1,
                    radius,
                    length / pieces as f64,
                ));
            }
            let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
            assert!(
                relative(response.z[0], expected) < 5.0e-10,
                "pieces={pieces}, detuning={detuning}: {:?} vs {expected:?}",
                response.z
            );
        }
    }
}

// G1: an antiresonance has zero driving pressure at finite flow. A truly
// undamped closed-cavity eigenfrequency still refuses; no regularizer hides it.
#[test]
fn network_half_wave_antiresonance_is_not_a_closed_cavity_pole() {
    let gas = gas();
    let radius = 0.008;
    let length = 0.2;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    let mut net = DuctNetwork {
        nodes: vec![
            port("input"),
            node(
                "end",
                NodeBoundary::Termination(NetworkLoad::PressureRelease),
            ),
        ],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    assert_eq!(response.wave_constraint_edges, 1);
    assert!(response.z[0].abs() / zc < 1.0e-12);
    net.nodes[1].boundary = NodeBoundary::Termination(NetworkLoad::Closed);
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
}

// G1/G3: a loaded half-wave connecting two driven ports has a finite,
// reciprocal rank-one impedance. Its singular pressure chart cannot be used
// to eliminate either port, but the complete wave/KCL system is regular.
#[test]
fn network_half_wave_multiport_retains_reciprocal_mutual_loading() {
    let gas = gas();
    let radius = 0.008;
    let length = 0.2;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    let net = DuctNetwork {
        nodes: vec![
            port("left"),
            port("right"),
            node(
                "load",
                NodeBoundary::Termination(NetworkLoad::Impedance(C64::from_re(zc))),
            ),
        ],
        edges: vec![
            edge("half-wave", 0, 1, radius, length),
            edge("load-line", 1, 2, radius, 0.137),
        ],
    };
    let response = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    assert_eq!(response.wave_constraint_edges, 2);
    for (observed, sign) in response.z.iter().zip([1.0, -1.0, -1.0, 1.0]) {
        assert!(relative(*observed, C64::from_re(sign * zc)) < 2.0e-12);
    }
    assert!(relative(response.z[1], response.z[2]) < 1.0e-12);
    let repeat = network_port_impedance(&net, &gas, omega, LossModel::Lossless).unwrap();
    for (a, b) in response.z.iter().zip(&repeat.z) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
}

#[test]
fn network_coupled_termination_even_mode_and_diagonal_limit() {
    // A symmetric Y-branch: equal branches see equal flows, so a symmetric
    // coupled load [[z, m], [m, z]] acts on each branch as the single load z+m.
    let gas = gas();
    let w = 2.0 * core::f64::consts::PI * 700.0;
    let net = DuctNetwork {
        nodes: vec![
            port("in"),
            node("j", NodeBoundary::Junction),
            port("a"),
            port("b"),
        ],
        edges: vec![
            edge("trunk", 0, 1, 0.008, 0.2),
            edge("left", 1, 2, 0.005, 0.1),
            edge("right", 1, 3, 0.005, 0.1),
        ],
    };
    let response = network_port_impedance(&net, &gas, w, LossModel::Bessel).unwrap();
    let (z, m) = (C64::new(1.0e5, -4.0e5), C64::new(5.0e4, -1.0e5));
    let coupled = |matrix: Vec<C64>, labels: &[&str]| CoupledTermination {
        omega: w,
        labels: labels.iter().map(|s| (*s).to_string()).collect(),
        z: matrix,
        source_id: "analytic symmetric coupled load".into(),
    };
    let single = |label: &str, value: C64| PortTermination {
        label: label.into(),
        load: NetworkLoad::Impedance(value),
    };
    let observed = response
        .drive_coupled("in", &[], &coupled(vec![z, m, m, z], &["a", "b"]), &gas)
        .unwrap()
        .impedance;
    let even = response
        .driving_point("in", &[single("a", z + m), single("b", z + m)], &gas)
        .unwrap();
    assert!(
        relative(observed, even) <= 1.0e-12,
        "{observed:?} vs {even:?}"
    );
    // The mutual term is live: the uncoupled answer differs.
    let uncoupled = response
        .driving_point("in", &[single("a", z), single("b", z)], &gas)
        .unwrap();
    assert!(relative(observed, uncoupled) > 1.0e-3);
    // Zero mutual terms reproduce the per-port loads bitwise, also when only
    // one port is coupled.
    let diagonal = response
        .drive_coupled(
            "in",
            &[],
            &coupled(vec![z, C64::ZERO, C64::ZERO, z], &["a", "b"]),
            &gas,
        )
        .unwrap()
        .impedance;
    assert_eq!(
        (diagonal.re.to_bits(), diagonal.im.to_bits()),
        (uncoupled.re.to_bits(), uncoupled.im.to_bits())
    );
    let mixed = response
        .drive_coupled("in", &[single("b", z)], &coupled(vec![z], &["a"]), &gas)
        .unwrap()
        .impedance;
    assert_eq!(
        (mixed.re.to_bits(), mixed.im.to_bits()),
        (uncoupled.re.to_bits(), uncoupled.im.to_bits())
    );
    // Refusals: active, duplicated, unknown or omitted ports; another
    // frequency; an active Hermitian part; a missing source.
    let ok = coupled(vec![z, m, m, z], &["a", "b"]);
    assert!(response.drive_coupled("a", &[], &ok, &gas).is_err());
    assert!(
        response
            .drive_coupled("in", &[single("a", z)], &ok, &gas)
            .is_err()
    );
    assert!(
        response
            .drive_coupled("in", &[], &coupled(vec![z, m, m, z], &["a", "x"]), &gas)
            .is_err()
    );
    assert!(
        response
            .drive_coupled("in", &[], &coupled(vec![z], &["a"]), &gas)
            .is_err()
    );
    let mut shifted = ok.clone();
    shifted.omega = w * (1.0 + 1e-9);
    assert!(response.drive_coupled("in", &[], &shifted, &gas).is_err());
    let active = coupled(
        vec![z, C64::new(3.0e5, 0.0), C64::new(3.0e5, 0.0), z],
        &["a", "b"],
    );
    assert!(response.drive_coupled("in", &[], &active, &gas).is_err());
    let mut anonymous = ok;
    anonymous.source_id = " ".into();
    assert!(response.drive_coupled("in", &[], &anonymous, &gas).is_err());
}
