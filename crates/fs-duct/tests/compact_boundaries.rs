//! Compact independent mouths and explicit signed multi-opening boundaries.
use fs_duct::compose::{ComposeOptions, TwoPort, compose_chain, duct_chain};
use fs_duct::network::{
    CoupledTermination, DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode,
    NodeBoundary, PortTermination, network_port_impedance,
};
use fs_duct::{Duct, HoleState, LossModel, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn close(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-11 * b.abs().max(1.0), "{a:?} != {b:?}");
}
fn mul(a: [C64; 4], b: [C64; 4]) -> [C64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

#[test]
fn compact_chimneys_keep_their_local_linear_stamps_for_every_actual_terminal() {
    let gas = gas();
    let (a, b, h, spacing) = (0.0075, 0.0035, 0.006, 0.024);
    let hole = Segment::ToneHole {
        hole_radius: b,
        chimney_height: h,
        bore_radius: a,
        state: HoleState::Open,
    };
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: a,
                length: 0.22,
            },
            hole.clone(),
            Segment::Cylinder {
                radius: a,
                length: spacing,
            },
            hole,
            Segment::Cylinder {
                radius: a,
                length: 0.12,
            },
        ],
    };
    let options = ComposeOptions::default();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for hz in [350.0, 950.0, 1150.0] {
            let omega = core::f64::consts::TAU * hz;
            let chain =
                duct_chain(&duct, &gas, omega, loss, None, options, true, || false).unwrap();
            // The declared compact approximation is a local, linear T junction:
            // its original chimney self load, plus the existing signed bore
            // series correction. No downstream bell is used in these stamps.
            let chimney =
                fs_duct::tone_hole_shunt(&gas, b, h, HoleState::Open, omega, loss, a).unwrap();
            let length = fs_phs::side_hole_series_length(b, a)
                + fs_phs::side_hole_mutual_length(b, b, a, spacing);
            let half = C64::new(
                0.0,
                -omega * gas.density * length / (core::f64::consts::PI * a * a),
            )
            .scale(0.5);
            let series = [C64::ONE, half, C64::ZERO, C64::ONE];
            let shunt = [C64::ONE, C64::ZERO, chimney.recip(), C64::ONE];
            let expected =
                TwoPort::from_abcd(mul(mul(series, shunt), series), chain.reference_impedance)
                    .unwrap();
            for index in [1, 3] {
                for (actual, target) in chain.blocks[index]
                    .scattering
                    .iter()
                    .zip(expected.scattering)
                {
                    close(*actual, target);
                }
            }
            let flanged =
                fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, a, omega)
                    .unwrap();
            for load in [
                None,
                Some(C64::ZERO),
                Some(flanged),
                Some(C64::new(1e6, 2e6)),
            ] {
                let mut reference_blocks = chain.blocks.clone();
                reference_blocks[1] = expected;
                reference_blocks[3] = expected;
                let reference = compose_chain(&reference_blocks, load, options).unwrap();
                let actual = match load {
                    None => fs_duct::input_impedance(&duct, &gas, omega, loss, Termination::Closed),
                    Some(z) => fs_duct::input_impedance_load(&duct, &gas, omega, loss, z, None),
                }
                .unwrap();
                close(actual.impedance, reference.z_in);
                close(actual.u_mouth_over_u_in, reference.u_out_over_u_in);
            }
        }
    }
}

#[test]
fn explicit_passive_matrix_preserves_negative_mutual_work_and_signed_reactance() {
    let gas = gas();
    let omega = core::f64::consts::TAU * 700.0;
    let port = |name: &str| NetworkNode {
        id: name.into(),
        boundary: NodeBoundary::Port { label: name.into() },
    };
    let edge = |id: &str, from, to, radius, length| NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    };
    let net = DuctNetwork {
        nodes: vec![
            port("in"),
            NetworkNode {
                id: "junction".into(),
                boundary: NodeBoundary::Junction,
            },
            port("left"),
            port("right"),
        ],
        edges: vec![
            edge("stem", 0, 1, 0.008, 0.2),
            edge("left", 1, 2, 0.005, 0.1),
            edge("right", 1, 3, 0.005, 0.1),
        ],
    };
    let response = network_port_impedance(&net, &gas, omega, LossModel::Bessel).unwrap();
    // A frequency-specific passive matrix: Hermitian eigenvalues 0.8e5,3.2e5.
    // Its even mode has LOWER resistance and opposite-signed reactance to z.
    // This is analytic matrix data, not a measured or geometry-baked boundary.
    let z = C64::new(2e5, -1e5);
    let mutual = C64::new(-1.2e5, 4e5);
    let boundary = CoupledTermination {
        omega,
        labels: vec!["left".into(), "right".into()],
        z: vec![z, mutual, mutual, z],
        source_id: "analytic passive signed-coupling fixture".into(),
    };
    let coupled = response.drive_coupled("in", &[], &boundary, &gas).unwrap();
    let loads = |value| {
        [
            PortTermination {
                label: "left".into(),
                load: NetworkLoad::Impedance(value),
            },
            PortTermination {
                label: "right".into(),
                load: NetworkLoad::Impedance(value),
            },
        ]
    };
    let even = response.drive("in", &loads(z + mutual), &gas).unwrap();
    close(coupled.impedance, even.impedance);
    for (actual, expected) in coupled.port_inflows.iter().zip(&even.port_inflows) {
        close(*actual, *expected);
    }
    let uncoupled = response.drive("in", &loads(z), &gas).unwrap();
    assert!((coupled.impedance - uncoupled.impedance).abs() > 1e-3 * coupled.impedance.abs());
    let (left, right) = (-coupled.port_inflows[1], -coupled.port_inflows[2]);
    close(left, right);
    let self_power = 0.5 * z.re * (left.norm_sq() + right.norm_sq());
    let mutual_power = 0.5 * (left.conj() * mutual * right + right.conj() * mutual * left).re;
    assert!(mutual_power < 0.0);
    assert!(self_power + mutual_power > 0.0);
    close(coupled.port_pressures[1], (z + mutual) * left);
    close(coupled.port_pressures[2], (z + mutual) * right);
}
