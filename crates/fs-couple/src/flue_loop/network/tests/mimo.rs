use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::{
    MimoRadiationLimits, MimoRadiationModel,
};
use crate::flue_loop::network::performance::{
    FlueOutletConfig, FlueOutletPerformance, FlueRadiatingOutlet, FlueRadiationInput,
};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctPressure, DuctRadiationSpec, DuctReconstruction, DuctView,
};
use crate::pcm_wav::observation::PressureRenderer;
use fs_math::c64::C64;
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};
use fs_vfit::mimo::MimoModel;
use fs_vfit::model::PoleTerm;

fn branched() -> FlueNetworkSpec {
    let g = gas();
    let mut source = source(g, -0.9);
    source.graph.nodes = vec![
        NetworkNode::Inlet,
        NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.9 },
        NetworkNode::Termination { reflection: -0.9 },
    ];
    source.graph.sections = vec![
        TubeSection {
            nodes: [0, 1],
            length_m: 0.12,
            radius_m: 0.006,
            max_length_error_m: 0.005,
        },
        TubeSection {
            nodes: [1, 2],
            length_m: 0.18,
            radius_m: 0.006,
            max_length_error_m: 0.005,
        },
        TubeSection {
            nodes: [1, 3],
            length_m: 0.23,
            radius_m: 0.005,
            max_length_error_m: 0.005,
        },
    ];
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, g.pressure).unwrap();
    source.gases = vec![g, g, hot];
    source
}
fn observation() -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: vec![
            DuctOpening {
                network: 0,
                node: 2,
                center_m: [0.0, 0.0],
            },
            DuctOpening {
                network: 0,
                node: 3,
                center_m: [0.06, 0.0],
            },
        ],
        exterior: Some(gas()),
        receiver: CircularOutletReceiver {
            position_m: [0.1, 0.0, 1.0],
            radial_rings: 4,
            angular_points: 16,
            maximum_frequency_hz: 1000.0,
        },
        reconstruction: DuctReconstruction::FiniteDifference {
            sample_rate_hz: 48_000,
        },
    }
}
fn radiation(mutual: f64) -> MimoRadiationModel {
    // An explicitly synthetic reciprocal passive matrix fixture. It is not a
    // measured flute or an independently validated radiation approximation.
    let z = gas().characteristic_impedance / (core::f64::consts::PI * 0.006f64.powi(2));
    let d = vec![0.15 * z, mutual * z, mutual * z, 0.2 * z];
    MimoRadiationModel {
        declaration_identity: [1; 32],
        realization_identity: [2; 32],
        exterior: gas(),
        apertures: vec![
            CircularOpening {
                center_m: [0.0, 0.0],
                radius_m: 0.006,
            },
            CircularOpening {
                center_m: [0.06, 0.0],
                radius_m: 0.005,
            },
        ],
        band_hz: [1.0, 24_000.0],
        impedance: MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real {
                pole: -3000.0,
                residue: 0.0,
            }],
            residues: vec![d.iter().map(|v| C64::from_re(100.0 * v)).collect()],
            d,
            e: vec![0.0; 4],
            reciprocal: true,
        },
    }
}
fn limits() -> MimoRadiationLimits {
    MimoRadiationLimits {
        max_states: 32,
        max_memory_bytes: 1 << 20,
        sampled_passivity_points: 257,
    }
}
fn coupled(mutual: f64) -> FlueVoice {
    voice(branched())
        .with_mimo_environment(
            radiation(mutual),
            observation(),
            limits(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}
fn schedule() -> GestureSchedule {
    GestureSchedule::try_new(
        700,
        vec![GestureTrack {
            id: "blow".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(300.0),
            events: vec![],
        }],
    )
    .unwrap()
}
fn config() -> FlueOutletConfig {
    FlueOutletConfig {
        sample_rate_hz: 48_000,
        samples: 2048,
        max_block: 2048,
        max_controls: 64,
        max_compile_work: 10_000,
        outlet_node: 2,
        receiver: observation().receiver,
        maximum_abs_pressure_pa: 1e6,
    }
}
fn performance(
    voice: FlueVoice,
    exterior: GasState,
) -> Result<FlueOutletPerformance, crate::render::schedule::GestureCompileError> {
    FlueOutletPerformance::with_radiation(
        voice,
        schedule(),
        config(),
        FlueRadiationInput {
            exterior,
            additional_outlets: vec![FlueRadiatingOutlet {
                node: 3,
                center_m: [0.06, 0.0],
            }],
        },
    )
}

#[test]
fn mutual_operator_changes_the_actual_jet_feedback_and_common_microphone() {
    let original = branched();
    let mut with = coupled(0.04);
    let mut without = coupled(0.0);
    assert_eq!(with.network().unwrap().source().graph, original.graph);
    assert_eq!(with.network().unwrap().section_gases(), original.gases);
    let environment = with.network().unwrap().mimo_environment().unwrap();
    assert_eq!(environment.source().declaration_identity, [1; 32]);
    assert_eq!(environment.source().realization_identity, [2; 32]);
    assert_eq!(environment.source().exterior, gas());
    let (continuous, sampled) = with.network().unwrap().mimo_passivity().unwrap();
    assert!(continuous.passive && sampled.passive);
    assert_eq!(
        sampled.class,
        fs_vfit::passivity::CertificateClass::GridOnly
    );
    let mut observer = DuctPressure::new(
        &[DuctView {
            graph: with.network().unwrap().graph(),
            gases: with.network().unwrap().section_gases(),
        }],
        observation(),
    )
    .unwrap();
    for voice in [&mut with, &mut without] {
        voice.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    }
    let mut a = vec![0.0; 2048];
    let mut b = vec![0.0; 2048];
    let mut expected_pressure = Vec::new();
    let mut coupled_work = 0.0;
    with.step_block_observed(&mut a, |line| {
        let FlueResonator::Network(network) = line else {
            panic!("physical graph");
        };
        let frame = network.last_frame().unwrap();
        coupled_work += frame.coupled_terminal_work_j.abs();
        assert!(frame.balance_residual_j().abs() < 1e-10);
        for (i, node) in [2, 3].iter().enumerate() {
            assert_eq!(
                network.node_frame(*node).unwrap().pressure_pa,
                network.mimo_port_frames().unwrap()[i].pressure_pa
            );
        }
        expected_pressure.push(
            observer
                .push_frame(|_, node| network.node_frame(node).map(|f| f.net_flow_into_node_m3_s))
                .unwrap()
                .unwrap(),
        );
        Ok(())
    })
    .unwrap();
    without.step_block(&mut b).unwrap();
    assert!(coupled_work > 0.0);
    assert!(a.iter().zip(&b).any(|(x, y)| (x - y).abs() > 1e-8));
    assert!(expected_pressure.iter().any(|p| p.abs() > 1e-12));
    let mut playback = performance(coupled(0.04), gas()).unwrap();
    let mut actual = vec![0.0; 2048];
    for block in actual.chunks_mut(37) {
        playback.block(block).unwrap();
    }
    assert_eq!(actual, expected_pressure);
    let mut independent = performance(coupled(0.0), gas()).unwrap();
    let mut other_pressure = vec![0.0; 2048];
    independent.block(&mut other_pressure).unwrap();
    assert!(
        actual
            .iter()
            .zip(&other_pressure)
            .any(|(a, b)| (a - b).abs() > 1e-10)
    );
}

#[test]
fn the_installed_environment_refuses_a_second_gas_or_aperture_identity() {
    let different = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas().pressure).unwrap();
    assert!(performance(coupled(0.04), different).is_err());
    let mut input = radiation(0.04);
    input.exterior = different;
    assert!(
        voice(branched())
            .with_mimo_environment(
                input,
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut input = radiation(0.04);
    input.apertures.swap(0, 1);
    assert!(
        voice(branched())
            .with_mimo_environment(
                input,
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut input = radiation(0.04);
    input.band_hz[1] *= 1.0 - 1e-9;
    assert!(
        voice(branched())
            .with_mimo_environment(
                input,
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut input = radiation(0.04);
    input.impedance.e[0] = 0.1;
    assert!(
        voice(branched())
            .with_mimo_environment(
                input,
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut input = radiation(0.04);
    input.impedance.d[1] *= 0.5;
    assert!(
        voice(branched())
            .with_mimo_environment(
                input,
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut running = coupled(0.04);
    assert!(running.set_network_terminal_reflection(2, -0.5).is_err());
    running.step_block(&mut [0.0]).unwrap();
    assert!(
        running
            .with_mimo_environment(
                radiation(0.04),
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
}

#[test]
fn authored_reactive_outlets_are_not_erased_by_mutual_loading() {
    let mut source = branched();
    source.graph.nodes[2] = NetworkNode::Impedance {
        load: fs_vfit::impedance::SeriesImpedanceSpec {
            resistance_pa_s_m3: 1e5,
            inertance_pa_s2_m3: 10.0,
            compliance_m3_pa: None,
        },
    };
    assert!(
        voice(source)
            .with_mimo_environment(
                radiation(0.04),
                observation(),
                limits(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
}
