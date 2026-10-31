//! Two 400 Hz peak-flow sources feeding a common matched receiver.
//! Run `cargo run -p fs-duct --example coherent_network -- 0` and repeat with
//! `180` degrees: pressure interferes before mean power is calculated.
//! The equal 0.2 m, 8 mm-radius rigid arms contain explicitly specified dry air.
//! This example uses one scalar impedance load. Multiple exterior mouths can
//! instead take their full sourced `CoupledTermination` through the final
//! argument; a diagonal scalar load does not stand in for mutual radiation.

use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary,
    PortFlowSource, PortTermination, network_drive_coherent,
};
use fs_duct::{Duct, LossModel, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 1 {
        return Err("usage: coherent_network <relative-source-phase-degrees>".into());
    }
    let degrees: f64 = args[0].parse()?;
    if !degrees.is_finite() {
        return Err("source phase must be finite".into());
    }
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)?;
    let radius = 0.008;
    let network = DuctNetwork {
        nodes: ["left", "right", "receiver"]
            .map(|label| NetworkNode {
                id: label.into(),
                boundary: NodeBoundary::Port {
                    label: label.into(),
                },
            })
            .to_vec(),
        edges: [0, 1]
            .map(|from| NetworkEdge {
                id: format!("arm-{from}"),
                from,
                to: 2,
                duct: Duct {
                    segments: vec![Segment::Cylinder {
                        radius,
                        length: 0.2,
                    }],
                },
                corrections: JunctionCorrections::Disabled,
            })
            .to_vec(),
    };
    let (sin, cos) = det::sin_cos(degrees.rem_euclid(360.0).to_radians());
    let sources = [
        PortFlowSource {
            label: "left".into(),
            volume_flow_m3_s: C64::from_re(1.0e-6),
        },
        PortFlowSource {
            label: "right".into(),
            volume_flow_m3_s: C64::new(cos, sin).scale(1.0e-6),
        },
    ];
    let characteristic = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(characteristic)),
    }];
    let response = network_drive_coherent(
        &network,
        &gas,
        core::f64::consts::TAU * 400.0,
        LossModel::Lossless,
        &sources,
        &loads,
        None,
    )?;
    println!("frequency_hz=400,phase_degrees={degrees},peak_source_flow_m3_s=1e-6");
    println!(
        "port,pressure_re_pa,pressure_im_pa,flow_re_m3_s,flow_im_m3_s,mean_power_into_graph_w"
    );
    for (index, label) in response.ports.iter().enumerate() {
        let p = response.port_pressures[index];
        let u = response.port_inflows[index];
        println!(
            "{label},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            p.re, p.im, u.re, u.im, response.port_powers_w[index]
        );
    }
    println!(
        "source_power_w={:.12e},external_load_power_w={:.12e},net_graph_work_w={:.12e}",
        response.source_power_w, response.external_load_power_w, response.network_power_w
    );
    Ok(())
}
