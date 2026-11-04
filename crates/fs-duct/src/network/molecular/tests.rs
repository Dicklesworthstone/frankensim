use super::*;
use crate::network::{
    JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary, PortExcitation,
};
use crate::{Duct, HoleState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::GasRecord;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101325.0).unwrap()
}

// Manufactured reservoirs isolate constitutive effects; they are not empirical
// dry-air relaxation measurements or an experimental qualification.
fn law(gas: GasState, count: usize, rate_scale: f64) -> RelaxationState {
    let id = SpeciesId::new("synthetic").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|j| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * j as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (2e-5 + j as f64 * 3e-4) * rate_scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
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
fn records(net: &DuctNetwork, states: &[RelaxationState]) -> NetworkMedia {
    NetworkMedia::new(
        net,
        net.edges
            .iter()
            .zip(states)
            .map(|(edge, state)| (edge.id.clone(), GasRecord::Legacy(*state.gas())))
            .collect(),
        GasRecord::Legacy(gas(293.15)),
    )
    .unwrap()
}
fn bindings(
    net: &DuctNetwork,
    states: &[RelaxationState],
    wall: Option<MolecularWallResolution>,
) -> Vec<MolecularEdge> {
    net.edges
        .iter()
        .zip(states)
        .map(|(edge, state)| MolecularEdge {
            edge_id: edge.id.clone(),
            state: state.clone(),
            wall,
            cones: ConeSlices::Fixed(16),
        })
        .collect()
}
fn model(
    net: DuctNetwork,
    states: &[RelaxationState],
    wall: Option<MolecularWallResolution>,
) -> MolecularNetwork {
    let media = records(&net, states);
    let edges = bindings(&net, states, wall);
    MolecularNetwork::new(net, media, edges, 4096).unwrap()
}
fn near(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-10 * scale,
        "actual={actual:?}, expected={expected:?}, scale={scale:e}"
    );
}
fn sin_cos(z: C64) -> (C64, C64) {
    let (s, c) = z.re.sin_cos();
    (
        C64::new(s * z.im.cosh(), c * z.im.sinh()),
        C64::new(c * z.im.cosh(), -s * z.im.sinh()),
    )
}
// Independent physical p/U transfer from the material's dynamic bulk law.
// This deliberately uses no duct transfer, scattering or network routine.
fn line(state: &RelaxationState, omega: f64, radius: f64, length: f64) -> [C64; 4] {
    let k = state.relaxing_wavenumber(omega).unwrap();
    let z =
        C64::from_re(omega * state.gas().density / (core::f64::consts::PI * radius * radius)) / k;
    let (s, c) = sin_cos(k.scale(length));
    let minus_i = C64::new(0.0, -1.0);
    [c, minus_i * z * s, minus_i * s / z, c]
}
fn drive(prepared: &PreparedNetworkSources) -> crate::network::NetworkSourceResponse {
    let mut result = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut result)
        .unwrap();
    result
}

#[test]
fn molecular_network_matches_complex_line_storage_loss_and_finite_source_work() {
    let state = law(gas(293.15), 2, 1.0);
    let radius = 0.008;
    let length = 0.27;
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let original = model(net.clone(), &[state.clone()], None);
    let changed = model(net, &[law(gas(293.15), 2, 0.1)], None);
    let reference =
        state.gas().characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let source_z = C64::new(0.3, -0.2).scale(reference);
    let load_z = C64::new(0.8, 0.4).scale(reference);
    let amplitude = C64::new(1.2, -0.3);
    let sources = [PortSource {
        label: "drive".into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: amplitude,
            impedance_pa_s_m3: source_z,
        },
    }];
    let loads = [PortTermination {
        label: "load".into(),
        load: NetworkLoad::Impedance(load_z),
    }];
    let mut effect = 0.0_f64;
    for omega in [100.0, 1000.0, 7000.0, 30000.0] {
        let actual = drive(
            &original
                .prepare_sources(omega, &sources, None, &loads, None)
                .unwrap(),
        );
        let [a, b, c, d] = line(&state, omega, radius, length);
        let outlet_flow = amplitude / (a * load_z + b + source_z * (c * load_z + d));
        let p = [(a * load_z + b) * outlet_flow, load_z * outlet_flow];
        let u = [(c * load_z + d) * outlet_flow, -outlet_flow];
        for i in 0..2 {
            near(actual.fields.port_pressures[i], p[i], 3.0);
            near(actual.fields.port_inflows[i], u[i], 3.0 / reference);
        }
        let work = 0.5 * (amplitude * u[0].conj()).re;
        let source_loss = 0.5 * source_z.re * u[0].norm_sq();
        let gas_loss = 0.5 * ((p[0] * u[0].conj()).re + (p[1] * u[1].conj()).re);
        assert!((actual.generator_power_w - work).abs() < 1e-10 / reference);
        assert!((actual.source_dissipation_w - source_loss).abs() < 1e-10 / reference);
        assert!((actual.fields.network_power_w - gas_loss).abs() < 1e-10 / reference);
        assert!(actual.fields.network_power_w > 0.0);
        let alternative = drive(
            &changed
                .prepare_sources(omega, &sources, None, &loads, None)
                .unwrap(),
        );
        effect = effect.max((alternative.fields.port_pressures[1] - p[1]).abs());
    }
    assert!(
        effect > 1e-3,
        "actual transmitted fields must retain relaxation times"
    );
    assert_eq!(original.edge("tube").unwrap().state(), &state);
}

#[test]
fn regional_molecular_branches_retain_full_source_bank_and_coherent_interference() {
    let states = [law(gas(293.15), 2, 1.0), law(gas(343.15), 2, 0.2)];
    let lengths = [0.11, 0.17];
    let radius = 0.008;
    let net = DuctNetwork {
        nodes: vec![port("left"), port("right"), port("receiver")],
        edges: vec![
            edge("left-arm", 0, 2, radius, lengths[0]),
            edge("right-arm", 1, 2, radius, lengths[1]),
        ],
    };
    let media = records(&net, &states);
    let mut reordered = net.clone();
    reordered.edges.reverse();
    for edge in &mut reordered.edges {
        edge.corrections = JunctionCorrections::ChainLocal;
    }
    let mut reordered_bindings = bindings(&net, &states, None);
    reordered_bindings.reverse();
    let reordered = MolecularNetwork::new(reordered, media, reordered_bindings, 16).unwrap();
    let original = model(net, &states, None);
    let r = states[0].gas().characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let source_z = [
        C64::new(0.7, -0.2),
        C64::new(-0.35, 0.6),
        C64::new(-0.35, 0.6),
        C64::new(0.9, 0.4),
    ]
    .map(|z| z.scale(r));
    let load_z = C64::new(0.8, -0.3).scale(r);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load_z),
    }];
    for omega in [713.0, 5107.0, 22000.0] {
        let bank = CoupledSourceBank {
            impedance: CoupledTermination {
                omega,
                labels: vec!["left".into(), "right".into()],
                z: source_z.to_vec(),
                source_id: "synthetic reciprocal source pair".into(),
            },
            open_circuit_pressures_pa: vec![C64::ZERO; 2],
        };
        let prepared = original
            .prepare_sources(omega, &[], Some(&bank), &loads, None)
            .unwrap();
        let reordered = reordered
            .prepare_sources(omega, &[], Some(&bank), &loads, None)
            .unwrap();
        assert_eq!(prepared.pressure_transfer(), reordered.pressure_transfer());
        assert_eq!(prepared.flow_transfer(), reordered.flow_transfer());
        let m = [
            line(&states[0], omega, radius, lengths[0]),
            line(&states[1], omega, radius, lengths[1]),
        ];
        // Eliminate the two physical line outlet currents and their common
        // receiver pressure, then close the FULL pressure-generator matrix.
        let receiver_z = (load_z.recip() + m[0][2] / m[0][3] + m[1][2] / m[1][3]).recip();
        let mut matrix = source_z;
        for i in 0..2 {
            for j in 0..2 {
                matrix[2 * i + j] = matrix[2 * i + j] + receiver_z / (m[i][3] * m[j][3]);
            }
            matrix[2 * i + i] = matrix[2 * i + i] + m[i][1] / m[i][3];
        }
        let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
        let mut powers = Vec::new();
        let mut actual = prepared.response_buffer();
        for phase in [0.0_f64, 0.8, core::f64::consts::PI] {
            let amplitude = [C64::new(1.2, -0.3), C64::new(phase.cos(), phase.sin())];
            let u = [
                (matrix[3] * amplitude[0] - matrix[1] * amplitude[1]) / determinant,
                (matrix[0] * amplitude[1] - matrix[2] * amplitude[0]) / determinant,
            ];
            let receiver = receiver_z * (u[0] / m[0][3] + u[1] / m[1][3]);
            let p = [
                (receiver + m[0][1] * u[0]) / m[0][3],
                (receiver + m[1][1] * u[1]) / m[1][3],
                receiver,
            ];
            prepared.drive_into(&amplitude, &mut actual).unwrap();
            for i in 0..3 {
                near(actual.fields.port_pressures[i], p[i], 4.0);
            }
            near(actual.fields.port_inflows[0], u[0], 4.0 / r);
            near(actual.fields.port_inflows[1], u[1], 4.0 / r);
            near(actual.fields.port_inflows[2], -receiver / load_z, 4.0 / r);
            assert!(actual.fields.network_power_w > 0.0);
            assert!(actual.source_dissipation_w > 0.0);
            assert!(
                (actual.generator_power_w
                    - actual.source_dissipation_w
                    - actual.fields.network_power_w
                    - actual.fields.external_load_power_w)
                    .abs()
                    < 1e-10 / r
            );
            powers.push(actual.fields.external_load_power_w);
        }
        assert!(
            (powers[0] - powers[2]).abs() > 1e-4 / r,
            "source phase must affect received work"
        );
    }
}

#[test]
fn joint_molecular_wall_spectra_survive_graph_subdivision_and_frequency_band_refusal() {
    let state = law(gas(293.15), 2, 1.0);
    let wall = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [0.1, 7000.0],
        maximum_relative_error: 0.02,
    };
    let net = DuctNetwork {
        nodes: vec![
            port("drive"),
            NetworkNode {
                id: "join".into(),
                boundary: NodeBoundary::Junction,
            },
            port("load"),
        ],
        edges: vec![
            edge("first", 0, 1, 0.001, 0.03),
            edge("second", 1, 2, 0.0016, 0.04),
        ],
    };
    let duct = Duct {
        segments: net
            .edges
            .iter()
            .flat_map(|edge| edge.duct.segments.clone())
            .collect(),
    };
    let chain =
        MolecularDuct::new(&duct, state.clone(), Some(wall), ConeSlices::Fixed(16), 16).unwrap();
    let graph = model(net, &[state.clone(), state], Some(wall));
    let flow = C64::new(1e-8, -0.3e-8);
    let load = C64::new(1e8, -2e7);
    let sources = [PortSource {
        label: "drive".into(),
        excitation: PortExcitation::VolumeFlow(flow),
    }];
    let loads = [PortTermination {
        label: "load".into(),
        load: NetworkLoad::Impedance(load),
    }];
    for omega in [10.0, 1000.0, 6000.0] {
        let expected = chain
            .response(omega, Some(load), ComposeOptions::default())
            .unwrap();
        let actual = drive(
            &graph
                .prepare_sources(omega, &sources, None, &loads, None)
                .unwrap(),
        );
        let transmitted = expected.u_out_over_u_in * flow;
        near(actual.fields.port_pressures[0], expected.z_in * flow, 10.0);
        near(actual.fields.port_pressures[1], load * transmitted, 10.0);
        near(actual.fields.port_inflows[1], -transmitted, 1e-7);
        assert!(actual.fields.network_power_w > 0.0);
    }
    assert!(
        graph
            .prepare_sources(7500.0, &sources, None, &loads, None)
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        graph.prepare_sources_controlled(1000.0, &sources, None, &loads, None, &mut || {
            calls += 1;
            calls == 3
        }),
        Err(MolecularDuctError::Duct(DuctError::Cancelled { .. }))
    ));
}

#[test]
fn empty_reservoir_graph_keeps_the_loaded_solution_at_a_half_wave_pole() {
    let gas = gas(293.15);
    let radius = 0.008;
    let length = 0.2;
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let graph = model(net, &[law(gas, 0, 1.0)], None);
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    let r = gas.characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let sources = [PortSource {
        label: "drive".into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: C64::ONE,
            impedance_pa_s_m3: C64::from_re(r),
        },
    }];
    let loads = [PortTermination {
        label: "load".into(),
        load: NetworkLoad::Impedance(C64::from_re(r)),
    }];
    let actual = drive(
        &graph
            .prepare_sources(omega, &sources, None, &loads, None)
            .unwrap(),
    );
    near(actual.fields.port_pressures[0], C64::from_re(0.5), 1.0);
    near(actual.fields.port_pressures[1], C64::from_re(-0.5), 1.0);
    near(
        actual.fields.port_inflows[1],
        C64::from_re(0.5 / r),
        1.0 / r,
    );
    assert!(actual.fields.network_power_w.abs() < 1e-10 / r);
}

#[test]
fn exact_material_identity_aggregate_resolution_and_explicit_chimneys_are_required() {
    let states = [law(gas(293.15), 2, 1.0), law(gas(293.15), 2, 0.1)];
    let net = DuctNetwork {
        nodes: vec![port("a"), port("b"), port("c")],
        edges: vec![
            edge("first", 0, 1, 0.008, 0.1),
            edge("second", 1, 2, 0.008, 0.1),
        ],
    };
    let media = records(&net, &states);
    let good = bindings(&net, &states, None);
    let mut wrong = good.clone();
    wrong[0].state = law(gas(303.15), 2, 1.0);
    assert!(MolecularNetwork::new(net.clone(), media.clone(), wrong, 16).is_err());
    let mut duplicate = good.clone();
    duplicate[1].edge_id = duplicate[0].edge_id.clone();
    assert!(MolecularNetwork::new(net.clone(), media.clone(), duplicate, 16).is_err());
    let mut missing = good.clone();
    missing[0].edge_id = "unknown".into();
    assert!(MolecularNetwork::new(net.clone(), media.clone(), missing, 16).is_err());
    assert!(MolecularNetwork::new(net.clone(), media.clone(), good.clone(), 1).is_err());
    assert!(matches!(
        MolecularNetwork::new_controlled(net.clone(), media.clone(), good.clone(), 16, &mut || {
            true
        }),
        Err(MolecularDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        }))
    ));
    let mut compact = net;
    compact.edges[0].duct.segments.push(Segment::ToneHole {
        hole_radius: 0.002,
        chimney_height: 0.003,
        bore_radius: 0.008,
        state: HoleState::Closed,
    });
    assert!(MolecularNetwork::new(compact, media, good, 16).is_err());
}
