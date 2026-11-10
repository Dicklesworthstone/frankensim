//! The actual file selector reaches the original sampled graph and heat owner.
use super::*;
use crate::bernoulli_aperture::viscothermal::{
    ViscothermalTimeOptions, with_regional_viscothermal_sections_options,
};
use fs_material::gas::GasSpec;

const DT: f64 = 1.0 / 96_000.0;

fn input(explicit_network: bool) -> Input {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    Input {
        graph: TubeNetworkSpec {
            nodes: vec![
                NetworkNode::Inlet,
                NetworkNode::Termination { reflection: 0.0 },
            ],
            sections: vec![TubeSection {
                nodes: [0, 1],
                length_m: 0.25,
                radius_m: 0.007,
                max_length_error_m: 1e-12,
            }],
            sound_speed_m_s: gas.sound_speed,
            max_wave_memory_bytes: 8 << 20,
        },
        gases: vec![gas],
        inner: Vec::new(),
        closed: Vec::new(),
        explicit_network,
    }
}

fn parse(
    input: &Input,
    text: &str,
) -> Result<(Vec<ViscothermalSectionSpec>, RelaxationSampling), FlueInputError> {
    let mut rows = Rows {
        lines: text.lines(),
        line: 0,
    };
    let result = input.read_losses(&mut rows)?;
    assert!(rows.lines.next().is_none());
    Ok(result)
}

#[test]
fn all_original_loss_headers_admit_one_explicit_sampling_choice() {
    for (network, header, tail, count) in [
        (false, "none", "", 0),
        (false, "16 50 1500", "", 1),
        (true, "sections 1", "loss_section 0 16 50 1500\n", 1),
    ] {
        let source = input(network);
        let original = parse(&source, &format!("viscothermal {header}\n{tail}")).unwrap();
        assert_eq!(original.1, RelaxationSampling::BilinearMidpoint);
        assert_eq!(original.0.len(), count);
        for (name, sampling) in [
            ("bilinear-midpoint", RelaxationSampling::BilinearMidpoint),
            ("exact-held-flow", RelaxationSampling::ExactHeldFlow),
        ] {
            let selected =
                parse(&source, &format!("viscothermal {header} {name}\n{tail}")).unwrap();
            assert_eq!(selected.1, sampling);
            assert_eq!(selected.0, original.0);
        }
        for suffix in [
            "exact",
            "backward-euler",
            "exact-held-flow ignored",
            "bilinear-midpoint exact-held-flow",
        ] {
            assert!(parse(&source, &format!("viscothermal {header} {suffix}\n{tail}")).is_err());
        }
    }
    assert!(parse(&input(true), "viscothermal 16 50 1500 exact-held-flow\n").is_err());
    assert!(parse(&input(true), "viscothermal sections 2 exact-held-flow\n").is_err());
}

#[test]
fn parsed_exact_loss_reaches_every_runtime_leg_and_retains_physical_work() {
    let gate = CancelGate::new_clock_free();
    let policy = FractionalTubePolicy {
        maximum_frequency_hz: 1500.0,
        maximum_summed_phase_error_rad: 0.05,
    };
    let mut outputs = Vec::new();
    for suffix in ["", " bilinear-midpoint", " exact-held-flow"] {
        let source = input(true);
        let (selections, sampling) = parse(
            &source,
            &format!("viscothermal sections 1{suffix}\nloss_section 0 16 50 1500\n"),
        )
        .unwrap();
        let (graph, gases, reports) = with_regional_viscothermal_sections_options(
            source.graph,
            source.gases,
            DT,
            &selections,
            ViscothermalTimeOptions {
                sampling,
                fractional: Some(policy),
            },
            &gate,
        )
        .unwrap();
        assert_eq!(reports[0].loss.sampling(), sampling);
        assert!(reports[0].max_scattering_error <= 0.03);
        let (mut viscous, mut thermal) = (0, 0);
        for node in &graph.nodes {
            match node {
                NetworkNode::Series { load } => {
                    assert_eq!(load.sampling(), sampling);
                    viscous += 1;
                }
                NetworkNode::ShuntAdmittance { load } => {
                    assert_eq!(load.sampling(), sampling);
                    thermal += 1;
                }
                _ => {}
            }
        }
        assert!(viscous > 0 && thermal > 0);
        let (mut network, _, _) = graph
            .realize_fractional_waveguide(&gases, DT, policy, &gate)
            .unwrap();
        let (mut heat, mut work) = (0.0, 0.0);
        let mut signal = Vec::new();
        for sample in 0..2048 {
            let before = network.stored_energy_j();
            let forcing = (core::f64::consts::TAU * f64::from(sample) / 192.0).sin();
            let frame = network.step(forcing).unwrap();
            let scale = before
                + frame.stored_energy_j
                + frame.inlet_work_j.abs()
                + frame.terminal_loss_j
                + frame.interior_loss_j;
            assert!(frame.balance_residual_j().abs() <= 3e-11 * scale.max(1e-30));
            assert!(frame.interior_loss_j >= 0.0);
            heat += frame.terminal_loss_j + frame.interior_loss_j;
            work += frame.inlet_work_j;
            signal.push(network.node_frame(1).unwrap().pressure_pa);
        }
        assert!(heat > 0.0 && signal.iter().any(|x| x.abs() > 1e-6));
        assert!((network.stored_energy_j() + heat - work).abs() < 3e-10 * work.abs());
        outputs.push(signal);
    }
    assert_eq!(
        outputs[0], outputs[1],
        "explicit midpoint must preserve the original file default"
    );
    assert_ne!(
        outputs[0], outputs[2],
        "exact selection must change actual memory evolution"
    );
}
