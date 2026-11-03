use super::*;
use crate::network::{
    JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary,
    network_drive_sources, network_port_impedance,
};
use crate::{Duct, Segment};
use fs_material::gas::GasSpec;
use fs_material::gas_species::input::GasRecord;
use std::f64::consts::PI;

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap()
}

fn port(label: &str) -> NetworkNode {
    NetworkNode {
        id: label.into(),
        boundary: NodeBoundary::Port {
            label: label.into(),
        },
    }
}

fn edge(label: &str, from: usize, to: usize, radius: f64, length: f64) -> NetworkEdge {
    NetworkEdge {
        id: label.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    }
}

fn reference(label: &str, value: f64) -> WavePort {
    WavePort {
        label: label.into(),
        reference_impedance_pa_s_m3: value,
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 3e-10 * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn exp_i(phase: f64) -> C64 {
    let (s, c) = det::sin_cos(phase);
    C64::new(c, s)
}

fn unitary(s: &[C64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            let inner = coherent_sum((0..n).map(|k| s[k * n + i].conj() * s[k * n + j]));
            close(inner, if i == j { C64::ONE } else { C64::ZERO });
            close(s[i * n + j], s[j * n + i]);
        }
    }
}

// G1: matched lossless pipe has no reflection and exactly its propagation
// phase in both transmission entries, including the undefined Z-matrix pole.
#[test]
fn matched_pipe_scattering_survives_open_circuit_poles() {
    let gas = gas(293.15);
    let length = 0.2;
    let radius = 0.008;
    let net = DuctNetwork {
        nodes: vec![port("right"), port("left")],
        edges: vec![edge("tube", 1, 0, radius, length)],
    };
    let zc = gas.density * gas.sound_speed / (PI * radius * radius);
    let ports = [reference("right", zc), reference("left", zc)];
    let pole = PI * gas.sound_speed / length;
    assert!(matches!(
        network_port_impedance(&net, &gas, pole, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    for omega in [0.73 * pole, pole, 2.0 * pole] {
        let response = PreparedNetworkScattering::new(
            &net,
            &gas,
            omega,
            LossModel::Lossless,
            &ports,
            &[],
            None,
        )
        .unwrap();
        assert_eq!(response.ports()[0].label, "left");
        assert_eq!(response.omega(), omega);
        let delay = exp_i(omega * length / gas.sound_speed);
        for (actual, expected) in
            response
                .pressure_matrix()
                .iter()
                .zip([C64::ZERO, delay, delay, C64::ZERO])
        {
            close(*actual, expected);
        }
        unitary(response.power_normalized_matrix(), 2);
    }
}

// G1/G3: unequal reference impedances require power normalization. Pressure
// transmission can exceed one in a passive system; that is not power gain.
#[test]
fn unequal_references_preserve_reciprocity_and_coherent_power() {
    let gas = gas(293.15);
    let length = 0.2;
    let radius = 0.008;
    let net = DuctNetwork {
        nodes: vec![port("a"), port("b")],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let zc = gas.density * gas.sound_speed / (PI * radius * radius);
    let ports = [reference("a", zc), reference("b", 4.0 * zc)];
    let response = PreparedNetworkScattering::new(
        &net,
        &gas,
        PI * gas.sound_speed / length,
        LossModel::Lossless,
        &ports,
        &[],
        None,
    )
    .unwrap();
    for (actual, expected) in response
        .pressure_matrix()
        .iter()
        .zip([0.6, -0.4, -1.6, -0.6])
    {
        close(*actual, C64::from_re(expected));
    }
    for (actual, expected) in response
        .power_normalized_matrix()
        .iter()
        .zip([0.6, -0.8, -0.8, -0.6])
    {
        close(*actual, C64::from_re(expected));
    }
    unitary(response.power_normalized_matrix(), 2);
    let incident = [C64::new(0.7, -0.3), C64::new(-0.4, 0.8)];
    let mut outgoing = [C64::ZERO; 2];
    response
        .scatter_power_normalized_into(&incident, &mut outgoing)
        .unwrap();
    let incoming_power: f64 = incident.iter().map(|z| z.norm_sq()).sum();
    let outgoing_power: f64 = outgoing.iter().map(|z| z.norm_sq()).sum();
    assert!((incoming_power - outgoing_power).abs() < 2e-10);
}

fn star(radii: [f64; 3], lengths: [f64; 3]) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![
            port("a"),
            port("b"),
            port("c"),
            NetworkNode {
                id: "junction".into(),
                boundary: NodeBoundary::Junction,
            },
        ],
        edges: (0..3)
            .map(|i| {
                edge(
                    ["a-edge", "b-edge", "c-edge"][i],
                    i,
                    3,
                    radii[i],
                    lengths[i],
                )
            })
            .collect(),
    }
}

// G1: independent travelling-wave elimination gives
// Sij=(2Yj/sum(Y)-delta_ij)*exp(i*(ki Li+kj Lj)).
#[test]
fn regional_three_way_reflection_matches_analytic_junction() {
    let gases = [gas(293.15), gas(343.15), gas(393.15)];
    let radii = [0.006, 0.009, 0.012];
    let lengths = [0.11, 0.13, 0.17];
    let net = star(radii, lengths);
    let media = NetworkMedia::new(
        &net,
        (0..3)
            .map(|i| {
                (
                    ["a-edge", "b-edge", "c-edge"][i].into(),
                    GasRecord::Legacy(gases[i]),
                )
            })
            .collect(),
        GasRecord::Legacy(gases[0]),
    )
    .unwrap();
    let z: [f64; 3] = std::array::from_fn(|i| {
        gases[i].density * gases[i].sound_speed / (PI * radii[i] * radii[i])
    });
    let ports: [WavePort; 3] = std::array::from_fn(|i| reference(["a", "b", "c"][i], z[i]));
    let omega = 2301.0;
    let response = PreparedNetworkScattering::regional(
        &net,
        &media,
        omega,
        LossModel::Lossless,
        &ports,
        &[],
        None,
    )
    .unwrap();
    let admittance: f64 = z.iter().map(|r| r.recip()).sum();
    for i in 0..3 {
        for j in 0..3 {
            let phase =
                omega * (lengths[i] / gases[i].sound_speed + lengths[j] / gases[j].sound_speed);
            let expected =
                exp_i(phase).scale(2.0 / z[j] / admittance - if i == j { 1.0 } else { 0.0 });
            close(response.pressure_matrix()[i * 3 + j], expected);
        }
    }
    unitary(response.power_normalized_matrix(), 3);
    let mut reordered = net.clone();
    reordered.edges.reverse();
    let reverse_ports = [ports[2].clone(), ports[0].clone(), ports[1].clone()];
    let again = PreparedNetworkScattering::regional(
        &reordered,
        &media,
        omega,
        LossModel::Lossless,
        &reverse_ports,
        &[],
        None,
    )
    .unwrap();
    assert_eq!(response.pressure_matrix(), again.pressure_matrix());
}

// G1/G3: original losses and explicit exterior loads remain in the coherent
// solve. Interference is combined BEFORE comparing incident/outgoing power.
#[test]
fn lossy_loaded_reflections_match_original_physical_pressure_sources() {
    let gas = gas(293.15);
    let radii = [0.006, 0.009, 0.012];
    let net = star(radii, [0.11, 0.13, 0.17]);
    let z: [f64; 3] = radii.map(|r| gas.density * gas.sound_speed / (PI * r * r));
    let ports = [reference("b", z[1]), reference("a", z[0])];
    let loads = [PortTermination {
        label: "c".into(),
        load: NetworkLoad::Impedance(C64::new(0.7 * z[2], -0.3 * z[2])),
    }];
    let omega = 1807.0;
    let response =
        PreparedNetworkScattering::new(&net, &gas, omega, LossModel::Bessel, &ports, &loads, None)
            .unwrap();
    for phase in [0.0, 0.7, PI] {
        let incoming = [C64::new(0.8, -0.2), exp_i(phase).scale(1.2)];
        let sources: Vec<_> = (0..2)
            .map(|i| PortSource {
                label: ["a", "b"][i].into(),
                excitation: PortExcitation::Pressure {
                    open_circuit_pressure_pa: incoming[i].scale(2.0),
                    impedance_pa_s_m3: C64::from_re(z[i]),
                },
            })
            .collect();
        let original =
            network_drive_sources(&net, &gas, omega, LossModel::Bessel, &sources, &loads, None)
                .unwrap();
        let mut outgoing = [C64::ZERO; 2];
        response.scatter_into(&incoming, &mut outgoing).unwrap();
        for i in 0..2 {
            close(
                outgoing[i],
                (original.fields.port_pressures[i] - original.fields.port_inflows[i].scale(z[i]))
                    .scale(0.5),
            );
        }
        let net_power: f64 = (0..2)
            .map(|i| 0.5 * (incoming[i].norm_sq() - outgoing[i].norm_sq()) / z[i])
            .sum();
        assert!(net_power >= 0.0);
        let physical_power =
            original.fields.network_power_w + original.fields.external_load_power_w;
        assert!((net_power - physical_power).abs() < 2e-10 * physical_power.abs().max(1e-9));
        let mut normalized = [C64::ZERO; 2];
        response
            .scatter_power_normalized_into(
                &std::array::from_fn::<_, 2, _>(|i| incoming[i].scale(1.0 / det::sqrt(z[i]))),
                &mut normalized,
            )
            .unwrap();
        for i in 0..2 {
            close(normalized[i].scale(det::sqrt(z[i])), outgoing[i]);
        }
    }
}

#[test]
fn invalid_ports_and_refused_application_publish_nothing() {
    let gas = gas(293.15);
    let net = DuctNetwork {
        nodes: vec![port("a"), port("b")],
        edges: vec![edge("tube", 0, 1, 0.008, 0.2)],
    };
    let zc = gas.density * gas.sound_speed / (PI * 0.008 * 0.008);
    let ports = [reference("a", zc), reference("b", zc)];
    for invalid in [
        vec![],
        vec![ports[0].clone(), ports[0].clone()],
        vec![reference("a", 0.0), ports[1].clone()],
        vec![reference("a", f64::NAN), ports[1].clone()],
        vec![reference("absent", zc), ports[1].clone()],
        vec![ports[0].clone()],
    ] {
        assert!(
            PreparedNetworkScattering::new(
                &net,
                &gas,
                2000.0,
                LossModel::Lossless,
                &invalid,
                &[],
                None
            )
            .is_err()
        );
    }
    let response =
        PreparedNetworkScattering::new(&net, &gas, 2000.0, LossModel::Lossless, &ports, &[], None)
            .unwrap();
    let sentinel = [C64::new(8.0, -9.0); 2];
    let mut output = sentinel;
    assert!(response.scatter_into(&[C64::ONE], &mut output).is_err());
    assert_eq!(output, sentinel);
    assert!(
        response
            .scatter_into(&[C64::ONE, C64::new(f64::NAN, 0.0)], &mut output)
            .is_err()
    );
    assert_eq!(output, sentinel);
    let mut short = [C64::ONE];
    assert!(response.scatter_into(&[C64::ONE; 2], &mut short).is_err());
    assert_eq!(short, [C64::ONE]);
}
