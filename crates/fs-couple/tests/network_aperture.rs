//! Public composition regressions. Synthetic valve geometry is explicitly a
//! solver fixture, not an instrument-fidelity or microphone validation claim.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::network_line::performance::NetworkLinePerformance;
use fs_couple::bernoulli_aperture::network_line::{LoadedApertureNetworkLine, NetworkLineLoad};
use fs_couple::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use fs_couple::bernoulli_aperture::performance::AperturePerformanceConfig;
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_couple::bernoulli_aperture::two_sided::LoadedAperture;
use fs_couple::driving_point::LineLadder;
use fs_couple::driving_point::network::NetworkDrivingPoint;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_dcontact::Obstacle;
use fs_duct::network::{
    CoupledTermination, DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkMedia,
    NetworkNode, NodeBoundary, PortTermination,
};
use fs_duct::{Duct, LossModel, Segment};
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_math::c64::C64;
use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortSchema, PortTimestamp, StableId};
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};
use fs_vfit::discretize::DelayedFilter;

const RATE: u32 = 96_000;
const SAMPLES: u64 = 256;
const AIR: &str = "mixture 293.15 100000 mole 2 N2 0.79 O2 0.21";
const CO2: &str = "mixture 293.15 100000 mole 1 CO2 1";

fn gas(text: &str) -> GasRecord {
    read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap()
}
fn schema() -> PortSchema {
    let id = |s| StableId::new(s).unwrap();
    PortKind::FluidPressureFlux
        .scalar_seed_schema(
            id("network-valve-test"),
            CoordinateBinding::new(
                id("pressure-volume-flow"),
                id("inlet-plane"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(id("sample-clock"), 0),
        )
        .unwrap()
}
fn graph() -> DuctNetwork {
    let edge = |id: &str, from, to| NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.006,
                length: 0.06,
            }],
        },
        corrections: JunctionCorrections::Disabled,
    };
    DuctNetwork {
        nodes: vec![
            NetworkNode {
                id: "source".into(),
                boundary: NodeBoundary::Port {
                    label: "inlet".into(),
                },
            },
            NetworkNode {
                id: "joint".into(),
                boundary: NodeBoundary::Junction,
            },
            NetworkNode {
                id: "end".into(),
                boundary: NodeBoundary::Port {
                    label: "outlet".into(),
                },
            },
        ],
        edges: vec![edge("source-edge", 0, 1), edge("remote-edge", 1, 2)],
    }
}
fn media(net: &DuctNetwork, remote: &str) -> NetworkMedia {
    NetworkMedia::new(
        net,
        vec![
            ("remote-edge".into(), gas(remote)),
            ("source-edge".into(), gas(AIR)),
        ],
        gas(AIR),
    )
    .unwrap()
}
fn terminal_z() -> f64 {
    let record = gas(AIR);
    1.7 * record.state().density * record.state().sound_speed
        / (core::f64::consts::PI * 0.006 * 0.006)
}
fn point(net: DuctNetwork) -> NetworkDrivingPoint {
    NetworkDrivingPoint::new(
        net,
        "inlet".into(),
        vec![PortTermination {
            label: "outlet".into(),
            load: NetworkLoad::Impedance(C64::new(terminal_z(), 0.0)),
        }],
        schema(),
    )
    .unwrap()
}
fn load(remote: &str, reversed: bool) -> NetworkLineLoad {
    let mut net = graph();
    let regions = media(&net, remote);
    if reversed {
        net.edges.reverse();
    }
    NetworkLineLoad::new(
        point(net),
        regions,
        LossModel::Lossless,
        RATE,
        &LineLadder::default(),
        None,
    )
    .unwrap()
}
fn contact() -> Obstacle {
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        1e8,
        2.0,
        "synthetic source-plane valve contact".into(),
    )
    .unwrap()
}
fn valve(load: &NetworkLineLoad) -> DynamicAperture {
    DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: 4e-4,
                width_m: 0.013,
                closing_pressure_pa: 6000.0,
            },
            mass_kg: 1e-5,
            stiffness_n_m: 500.0,
            damping_ratio: 0.35,
            density_kg_m3: load.density_kg_m3(),
            impedance_pa_s_m3: load.effective_impedance_pa_s_m3(),
            time_step_s: 1.0 / f64::from(RATE),
            max_steps: SAMPLES,
        },
        ApertureState {
            opening_m: 4e-4,
            opening_velocity_m_s: 0.0,
        },
        contact(),
        BladePropagator::ExactZoh,
    )
    .unwrap()
}
fn system(remote: &str, reversed: bool) -> LoadedApertureNetworkLine {
    let l = load(remote, reversed);
    LoadedApertureNetworkLine::new(valve(&l).into(), l).unwrap()
}
fn drive(n: u64) -> TubeDrive {
    TubeDrive {
        upstream_pressure_pa: if n < 100 { 50.0 } else { 0.0 },
        body_flow_m3_s: 0.0,
    }
}

#[test]
fn uniform_line_and_original_valve_retain_the_same_sampled_source_equations() {
    let l = load(AIR, false);
    let zc = l.characteristic_impedance_pa_s_m3();
    let (uniform, uniform_report) = l
        .driving_point()
        .realize_line(
            l.source_record().state(),
            LossModel::Lossless,
            RATE,
            zc,
            &LineLadder::default(),
            None,
        )
        .unwrap();
    let (mut reference_line, report) = l
        .driving_point()
        .realize_line_regional(
            l.media(),
            LossModel::Lossless,
            RATE,
            zc,
            &LineLadder::default(),
            None,
        )
        .unwrap();
    assert_eq!(&report, l.realization());
    assert_eq!(uniform_report.n_fft, report.n_fft);
    let taps = reference_line.quiescent_fir().unwrap().1.to_vec();
    for (a, b) in taps.iter().zip(uniform.quiescent_fir().unwrap().1) {
        assert!((a - b).abs() < 1e-12);
    }
    assert_eq!(taps.len(), l.tap_count());
    assert_eq!(taps[0], l.direct_reflection());
    let mut native = valve(&l);
    let mut coupled = LoadedApertureNetworkLine::new(valve(&l).into(), l).unwrap();
    let mut outgoing = Vec::new();
    for n in 0..SAMPLES {
        let h: f64 = taps
            .iter()
            .enumerate()
            .skip(1)
            .take(outgoing.len())
            .map(|(k, tap)| tap * outgoing[outgoing.len() - k])
            .sum();
        let input = drive(n);
        let expected = native
            .step(ApertureDrive {
                upstream_pressure_pa: input.upstream_pressure_pa,
                incoming_pressure_pa: h / (1.0 - taps[0]),
                body_flow_m3_s: input.body_flow_m3_s,
            })
            .unwrap();
        let actual = coupled.step(input).unwrap();
        assert!((actual.aperture.bore_flow_m3_s - expected.bore_flow_m3_s).abs() < 1e-16);
        assert!((actual.aperture.bore_pressure_pa - expected.bore_pressure_pa).abs() < 1e-9);
        assert!((actual.aperture.state.opening_m - expected.state.opening_m).abs() < 1e-14);
        let reflected = reference_line
            .push(actual.aperture.outgoing_pressure_pa)
            .unwrap();
        assert!((reflected - actual.incoming_pressure_pa).abs() < 1e-11);
        assert!(
            (actual.aperture.bore_pressure_pa - actual.aperture.outgoing_pressure_pa - reflected)
                .abs()
                < 1e-10
        );
        assert!(
            ((actual.aperture.outgoing_pressure_pa - reflected) / zc
                - actual.aperture.bore_flow_m3_s)
                .abs()
                < 1e-16
        );
        outgoing.push(actual.aperture.outgoing_pressure_pa);
    }

    let net = graph();
    let regions = media(&net, AIR);
    let model = NetworkDrivingPoint::new_coupled(
        net,
        "inlet".into(),
        vec![],
        vec!["outlet".into()],
        schema(),
    )
    .unwrap();
    let mut boundary = |omega| {
        Ok(CoupledTermination {
            omega,
            labels: vec!["outlet".into()],
            z: vec![C64::new(terminal_z(), 0.0)],
            source_id: "retained physical boundary source".into(),
        })
    };
    let coupled = NetworkLineLoad::new(
        model,
        regions,
        LossModel::Lossless,
        RATE,
        &LineLadder::default(),
        Some(&mut boundary),
    )
    .unwrap();
    assert_eq!(
        coupled.coupled_source_ids(),
        &["retained physical boundary source"]
    );
    assert_eq!(coupled.source_record(), &gas(AIR));
}

#[test]
fn regional_edge_identity_survives_reordering_and_a_remote_gas_changes_the_actual_trajectory() {
    let (mut air, mut changed, mut reordered) =
        (system(AIR, false), system(CO2, false), system(CO2, true));
    assert_eq!(changed.load().source_record(), &gas(AIR));
    assert_eq!(
        changed.load().media().edge_record("remote-edge"),
        Some(&gas(CO2))
    );
    assert_eq!(changed.load().media().exterior_record(), &gas(AIR));
    let mut changed_signal = false;
    for n in 0..SAMPLES {
        let a = air.step(drive(n)).unwrap();
        let b = changed.step(drive(n)).unwrap();
        let c = reordered.step(drive(n)).unwrap();
        assert!((b.aperture.bore_pressure_pa - c.aperture.bore_pressure_pa).abs() < 1e-9);
        changed_signal |= (a.aperture.bore_pressure_pa - b.aperture.bore_pressure_pa).abs() > 1e-4;
    }
    assert!(
        changed_signal,
        "a retained remote gas must affect feedback, not just metadata"
    );
}

#[test]
fn rejected_sample_preserves_source_line_and_last_observation() {
    let (mut a, mut reference) = (system(AIR, false), system(AIR, false));
    assert_eq!(a.step(drive(0)).unwrap(), reference.step(drive(0)).unwrap());
    let before = a.load().outgoing_history();
    let last = *a.last_frame().unwrap();
    assert!(
        a.step(TubeDrive {
            upstream_pressure_pa: f64::NAN,
            body_flow_m3_s: 0.0
        })
        .is_err()
    );
    assert!(
        a.step_with_forces_and_pressures(drive(1), [0.0, 1.0], [0.0; 2])
            .is_err()
    );
    assert!(
        a.step_with_forces_and_pressures(drive(1), [0.0; 2], [1.0, 0.0])
            .is_err()
    );
    assert_eq!(a.aperture().accepted_steps(), 1);
    assert_eq!(a.load().outgoing_history(), before);
    assert_eq!(a.last_frame(), Some(&last));
    assert_eq!(a.step(drive(1)).unwrap(), reference.step(drive(1)).unwrap());
}

#[test]
fn independent_signed_modes_keep_their_original_force_and_pressure_ports() {
    let l = load(AIR, false);
    let spec = ModalApertureSpec {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        modes: vec![
            ApertureMode {
                mass_kg: 1e-5,
                stiffness_n_m: 500.0,
                damping_n_s_m: 0.03,
                pressure_area_m2: 1.25e-4,
                gap_projection: 1.0,
                force_projection: [1.0, 0.0],
            },
            ApertureMode {
                mass_kg: 8e-6,
                stiffness_n_m: 2500.0,
                damping_n_s_m: 0.02,
                pressure_area_m2: -4e-5,
                gap_projection: -0.2,
                force_projection: [-0.3, 0.7],
            },
            ApertureMode {
                mass_kg: 6e-6,
                stiffness_n_m: 6000.0,
                damping_n_s_m: 0.01,
                pressure_area_m2: 7e-5,
                gap_projection: 0.0,
                force_projection: [0.2, -0.4],
            },
        ],
        density_kg_m3: l.density_kg_m3(),
        impedance_pa_s_m3: l.effective_impedance_pa_s_m3(),
        time_step_s: 1.0 / f64::from(RATE),
        max_steps: SAMPLES,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 1.0,
            max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13,
            relative_tolerance: 1e-11,
        },
    };
    let make = || {
        ModalAperture::new_with_propagator(
            spec.clone(),
            vec![BladeState::default(); 3],
            contact(),
            BladePropagator::ExactZoh,
        )
        .unwrap()
    };
    let mut native = make();
    let mut system = LoadedApertureNetworkLine::new(make().into(), l).unwrap();
    for n in 0..64 {
        let incoming =
            system.load().history_pressure_pa() / (1.0 - system.load().direct_reflection());
        let forces = [0.001, -0.0002];
        let expected = native
            .step(
                ApertureDrive {
                    upstream_pressure_pa: 50.0,
                    incoming_pressure_pa: incoming,
                    body_flow_m3_s: 0.0,
                },
                forces,
            )
            .unwrap();
        let actual = system
            .step_with_forces_and_pressures(drive(n), forces, [0.0; 2])
            .unwrap();
        let LoadedAperture::Modal(owner) = system.aperture() else {
            panic!("original modal owner was lost")
        };
        assert_eq!(owner.modal_states(), &expected.modes);
        assert_eq!(
            actual.aperture.mechanical_work_j,
            expected.aperture.mechanical_work_j
        );
        assert_eq!(
            actual.aperture.bore_pressure_pa,
            expected.aperture.bore_pressure_pa
        );
    }
    let LoadedAperture::Modal(owner) = system.aperture() else {
        unreachable!()
    };
    assert!(
        owner
            .modal_states()
            .iter()
            .all(|mode| mode.displacement_m.abs() > 1e-12)
    );
    assert_eq!(
        owner.spec().modes[2].gap_projection,
        0.0,
        "the swept-volume mode was not collapsed away"
    );
}

fn schedule() -> GestureSchedule {
    GestureSchedule::try_new(
        RATE,
        vec![GestureTrack {
            id: "blow".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.0),
            events: vec![
                GestureEvent {
                    time_s: 0.0,
                    transition_s: 0.0,
                    value: GestureValue::PressurePa(50.0),
                },
                GestureEvent {
                    time_s: 100.0 / f64::from(RATE),
                    transition_s: 0.0,
                    value: GestureValue::PressurePa(0.0),
                },
            ],
        }],
    )
    .unwrap()
}

#[test]
fn canonical_source_plane_pressure_is_invariant_to_callback_partition() {
    let mut reference = system(AIR, false);
    let expected: Vec<u64> = (0..SAMPLES)
        .map(|n| {
            reference
                .step(drive(n))
                .unwrap()
                .aperture
                .bore_pressure_pa
                .to_bits()
        })
        .collect();
    for block in [1, 17, 256] {
        let mut render = NetworkLinePerformance::new(
            system(AIR, false),
            schedule(),
            AperturePerformanceConfig {
                sample_rate_hz: RATE,
                samples: SAMPLES,
                max_block: block,
                max_compile_work: 4096,
                max_controls: 512,
            },
        )
        .unwrap();
        let mut actual = vec![0.0; SAMPLES as usize];
        for chunk in actual.chunks_mut(block) {
            render.block(chunk).unwrap();
        }
        assert_eq!(
            actual.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
            expected
        );
        assert_eq!(render.samples_rendered(), SAMPLES);
        assert_eq!(render.pending_controls(), 0);
        assert!(
            actual[110..].iter().any(|p| p.abs() > 1e-8),
            "release retains mechanical and wave history"
        );
    }
}

#[test]
fn full_quiescent_fir_inspection_has_no_probe_horizon_and_cannot_drop_hidden_state() {
    let mut taps = vec![0.0; 9000];
    taps[0] = 0.05;
    taps[8999] = 0.2;
    let mut line =
        DelayedFilter::from_impulse_response(1.0 / f64::from(RATE), taps.clone()).unwrap();
    let (dt, actual) = line.quiescent_fir().unwrap();
    assert_eq!(dt, 1.0 / f64::from(RATE));
    assert_eq!(actual, taps);
    assert!(
        line.clone()
            .with_dispersion(0.2, 1)
            .unwrap()
            .quiescent_fir()
            .is_none()
    );
    line.push(1.0).unwrap();
    assert!(line.quiescent_fir().is_none());
    for _ in 1..8999 {
        line.push(0.0).unwrap();
    }
    assert_eq!(line.push(0.0).unwrap(), 0.2);
}

#[test]
fn an_ambiguous_source_plane_or_self_loop_cannot_supply_an_arbitrary_gas_or_radius() {
    let mut net = graph();
    net.edges[1].from = 0;
    let regions = media(&net, AIR);
    let p = point(net);
    assert!(p.active_port_radius().is_err());
    assert!(p.active_port_record(&regions).is_err());
    assert!(
        NetworkLineLoad::new(
            p,
            regions,
            LossModel::Lossless,
            RATE,
            &LineLadder::default(),
            None
        )
        .is_err()
    );

    let net = DuctNetwork {
        nodes: vec![NetworkNode {
            id: "loop".into(),
            boundary: NodeBoundary::Port {
                label: "inlet".into(),
            },
        }],
        edges: vec![NetworkEdge {
            id: "self-loop".into(),
            from: 0,
            to: 0,
            duct: Duct {
                segments: vec![Segment::Cylinder {
                    radius: 0.006,
                    length: 0.12,
                }],
            },
            corrections: JunctionCorrections::Disabled,
        }],
    };
    let regions = NetworkMedia::new(&net, vec![("self-loop".into(), gas(AIR))], gas(AIR)).unwrap();
    let p = NetworkDrivingPoint::new(net, "inlet".into(), vec![], schema()).unwrap();
    assert!(p.active_port_radius().is_err());
    assert!(p.active_port_record(&regions).is_err());
}
