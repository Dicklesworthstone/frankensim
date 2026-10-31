//! Two 400 Hz peak-flow sources feeding a common matched receiver.
//! Run `cargo run -p fs-duct --example coherent_network -- 0` and repeat with
//! `180` degrees: pressure interferes before mean power is calculated.
//! The equal 0.2 m, 8 mm-radius rigid arms contain explicitly specified dry air.
//! This example uses one scalar impedance load. Multiple exterior mouths can
//! instead take their full sourced `CoupledTermination` through the final
//! argument; a diagonal scalar load does not stand in for mutual radiation.
//! An optional second argument sets each source resistance as a multiple of
//! Zc. This selects pressure generators with |p_open| = Zc * 1e-6 Pa, whose
//! actual flows react to each other; e.g. `-- 180 0.5` cancels receiver pressure
//! while retaining work dissipated inside the two source resistances.

use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary,
    PortExcitation, PortSource, PortTermination, network_drive_sources,
};
use fs_duct::{Duct, LossModel, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(1..=2).contains(&args.len()) {
        return Err(
            "usage: coherent_network <relative-source-phase-degrees> [source-resistance-over-Zc]"
                .into(),
        );
    }
    let degrees: f64 = args[0].parse()?;
    if !degrees.is_finite() {
        return Err("source phase must be finite".into());
    }
    let source_resistance = args.get(1).map(|value| value.parse::<f64>()).transpose()?;
    if source_resistance.is_some_and(|value| !value.is_finite() || value < 0.0) {
        return Err("source resistance / Zc must be finite and nonnegative".into());
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
    let characteristic = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let excitation = |phase: C64| match source_resistance {
        Some(ratio) => PortExcitation::Pressure {
            open_circuit_pressure_pa: phase.scale(characteristic * 1.0e-6),
            impedance_pa_s_m3: C64::from_re(ratio * characteristic),
        },
        None => PortExcitation::VolumeFlow(phase.scale(1.0e-6)),
    };
    let sources = [
        PortSource {
            label: "left".into(),
            excitation: excitation(C64::ONE),
        },
        PortSource {
            label: "right".into(),
            excitation: excitation(C64::new(cos, sin)),
        },
    ];
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(C64::from_re(characteristic)),
    }];
    let driven = network_drive_sources(
        &network,
        &gas,
        core::f64::consts::TAU * 400.0,
        LossModel::Lossless,
        &sources,
        &loads,
        None,
    )?;
    if let Some(ratio) = source_resistance {
        println!(
            "frequency_hz=400,phase_degrees={degrees},source_resistance_over_Zc={ratio},peak_open_pressure_pa={:.12e}",
            characteristic * 1.0e-6
        );
    } else {
        println!("frequency_hz=400,phase_degrees={degrees},peak_source_flow_m3_s=1e-6");
    }
    let response = &driven.fields;
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
    println!(
        "generator_power_w={:.12e},source_internal_dissipation_w={:.12e}",
        driven.generator_power_w, driven.source_dissipation_w
    );
    Ok(())
}
