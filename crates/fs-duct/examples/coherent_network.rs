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
//! Optional third/fourth arguments add reciprocal mutual source resistance
//! and reactance / Zc, e.g. `-- 180 0.5 -0.25 0.1`. The complete internal
//! matrix is admitted as passive; mutual entries and port work retain signs.
//! An optional fifth argument supplies coherence magnitude in [0,1], retaining
//! the declared phase: `-- 180 0.5 -0.25 0.1 0.6`. The prepared boundary serves
//! both the deterministic vector and its explicitly declared source statistics.

use fs_duct::network::{
    CoupledSourceBank, CoupledTermination, DuctNetwork, JunctionCorrections,
    NetworkCoherenceBudget, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary, PortExcitation,
    PortSource, PortTermination, network_prepare_sources,
};
use fs_duct::{Duct, LossModel, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(1..=5).contains(&args.len()) {
        return Err(
            "usage: coherent_network <relative-source-phase-degrees> [source-resistance-over-Zc [mutual-resistance-over-Zc [mutual-reactance-over-Zc [coherence-magnitude]]]]"
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
    let mutual_resistance = args.get(2).map(|value| value.parse::<f64>()).transpose()?;
    let mutual_reactance = args
        .get(3)
        .map(|value| value.parse::<f64>())
        .transpose()?
        .unwrap_or(0.0);
    let coherence = args.get(4).map(|value| value.parse::<f64>()).transpose()?;
    if coherence.is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value)) {
        return Err("coherence magnitude must be finite and in [0,1]".into());
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
    let omega = core::f64::consts::TAU * 400.0;
    let prepared = if let Some(mutual) = mutual_resistance {
        let resistance = source_resistance.expect("mutual arguments follow source resistance");
        let diagonal = C64::from_re(resistance * characteristic);
        let off_diagonal = C64::new(mutual, mutual_reactance).scale(characteristic);
        let bank = CoupledSourceBank {
            impedance: CoupledTermination {
                omega,
                labels: vec!["left".into(), "right".into()],
                z: vec![diagonal, off_diagonal, off_diagonal, diagonal],
                source_id: "explicit reciprocal two-generator example".into(),
            },
            open_circuit_pressures_pa: [C64::ONE, C64::new(cos, sin)]
                .map(|phase| phase.scale(characteristic * 1.0e-6))
                .to_vec(),
        };
        network_prepare_sources(
            &network,
            &gas,
            omega,
            LossModel::Lossless,
            &[],
            Some(&bank),
            &loads,
            None,
        )?
    } else {
        network_prepare_sources(
            &network,
            &gas,
            omega,
            LossModel::Lossless,
            &sources,
            None,
            &loads,
            None,
        )?
    };
    let mut driven = prepared.response_buffer();
    prepared.drive_into(prepared.default_excitations(), &mut driven)?;
    if let Some(ratio) = source_resistance {
        println!(
            "frequency_hz=400,phase_degrees={degrees},source_resistance_over_Zc={ratio},peak_open_pressure_pa={:.12e}",
            characteristic * 1.0e-6
        );
    } else {
        println!("frequency_hz=400,phase_degrees={degrees},peak_source_flow_m3_s=1e-6");
    }
    if let Some(mutual) = mutual_resistance {
        println!(
            "mutual_source_resistance_over_Zc={mutual},mutual_source_reactance_over_Zc={mutual_reactance}"
        );
    }
    let response = &driven.fields;
    println!(
        "port,pressure_re_pa,pressure_im_pa,flow_re_m3_s,flow_im_m3_s,mean_power_into_graph_w,generator_power_w,signed_internal_source_work_w"
    );
    for (index, label) in response.ports.iter().enumerate() {
        let p = response.port_pressures[index];
        let u = response.port_inflows[index];
        println!(
            "{label},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            p.re,
            p.im,
            u.re,
            u.im,
            response.port_powers_w[index],
            driven.generator_powers_w[index],
            driven.source_dissipations_w[index]
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
    if let Some(magnitude) = coherence {
        let defaults = prepared.default_excitations();
        let factors = [
            defaults[0],
            defaults[1].scale(magnitude),
            C64::ZERO,
            defaults[1].scale((1.0 - magnitude * magnitude).sqrt()),
        ];
        let mut statistics = prepared.statistics_buffer();
        prepared.coherence_into(
            &factors,
            2,
            NetworkCoherenceBudget::default(),
            &mut statistics,
        )?;
        let work = statistics
            .work()
            .expect("completed statistical contraction");
        let receiver_pressure_ms = statistics
            .pressure_second_moments()
            .expect("completed pressure moments")[8]
            .re;
        println!(
            "coherence_magnitude={magnitude},receiver_mean_square_pressure_pa2={receiver_pressure_ms:.12e},statistical_generator_power_w={:.12e},statistical_source_loss_w={:.12e},statistical_external_power_w={:.12e},statistical_graph_power_w={:.12e}",
            work.generator_power_w,
            work.source_dissipation_w,
            work.external_load_power_w,
            work.network_power_w
        );
    }
    Ok(())
}
