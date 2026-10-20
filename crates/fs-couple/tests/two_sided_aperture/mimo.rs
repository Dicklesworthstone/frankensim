//! Source-bound synthetic exterior data on the original three mechanical owners.
//! These integration controls are not measured instrument radiation models.
use super::*;
use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::radiation::environment::multiport::{
    MimoRadiationLimits, MimoRadiationModel,
};
use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use fs_couple::bernoulli_aperture::viscothermal::ViscothermalSectionSpec;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::baffled::openings::CircularOpening;
use fs_couple::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctPressure, DuctRadiationSpec, DuctReconstruction, DuctView,
};
use fs_couple::pcm_wav::observation::{ObservationPolicy, PressureRenderer};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::model::PoleTerm;

const SIDES: [ValveSide; 2] = [ValveSide::Upstream, ValveSide::Downstream];
#[derive(Clone, Copy)]
enum Kind { Single, Paired, Modal }
const KINDS: [Kind; 3] = [Kind::Single, Kind::Paired, Kind::Modal];

fn branched_sides() -> [ApertureDuctSpec; 2] {
    let g = gas();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, g.pressure).unwrap();
    [0, 4].map(|extra| ApertureDuctSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
                NetworkNode::Termination { reflection: -0.8 },
                NetworkNode::Termination { reflection: -0.7 }],
            sections: vec![
                TubeSection { nodes: [0, 1], length_m: (8 + extra) as f64 * g.sound_speed * DT,
                    radius_m: 0.007, max_length_error_m: 1e-12 },
                TubeSection { nodes: [1, 2], length_m: 12.0 * g.sound_speed * DT,
                    radius_m: 0.007, max_length_error_m: 1e-12 },
                TubeSection { nodes: [1, 3], length_m: 20.0 * hot.sound_speed * DT,
                    radius_m: 0.006, max_length_error_m: 1e-12 },
            ],
            sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
        },
        gases: vec![g, g, hot], viscothermal: vec![],
    })
}

fn valve(kind: Kind, sides: &[ApertureDuctSpec; 2]) -> LoadedAperture {
    let (g, z) = combined_inlet_state(sides).unwrap();
    match kind {
        Kind::Single => DynamicAperture::new_with_propagator(DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013,
                closing_pressure_pa: 1000.0 },
            mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.05,
            density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: 1024,
        }, ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 }, contact(),
            BladePropagator::ExactZoh).unwrap().into(),
        Kind::Paired => paired(z, g, 1024, 1.5).into(),
        Kind::Modal => ModalAperture::new_with_propagator(ModalApertureSpec {
            rest_opening_m: 4e-4, width_m: 0.013,
            modes: vec![
                ApertureMode { mass_kg: 1e-5, stiffness_n_m: 500.0, damping_n_s_m: 0.03,
                    pressure_area_m2: 1.25e-4, gap_projection: 1.0, force_projection: [1.0, 0.0] },
                ApertureMode { mass_kg: 8e-6, stiffness_n_m: 2500.0, damping_n_s_m: 0.02,
                    pressure_area_m2: -4e-5, gap_projection: -0.2, force_projection: [-0.3, 0.7] },
                ApertureMode { mass_kg: 6e-6, stiffness_n_m: 6000.0, damping_n_s_m: 0.01,
                    pressure_area_m2: 7e-5, gap_projection: 0.0, force_projection: [0.2, -0.4] },
            ],
            density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: 1024,
            limits: PairedApertureLimits { maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
                max_evaluations: 256, velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
        }, vec![BladeState::default(); 3], contact(), BladePropagator::ExactZoh).unwrap().into(),
    }
}

fn local_observation(side: ValveSide) -> DuctRadiationSpec {
    let centers = match side {
        ValveSide::Upstream => [[-0.09, 0.0], [-0.03, 0.0]],
        ValveSide::Downstream => [[0.03, 0.0], [0.09, 0.0]],
    };
    DuctRadiationSpec {
        openings: [2, 3].into_iter().zip(centers).map(|(node, center_m)|
            DuctOpening { network: 0, node, center_m }).collect(),
        exterior: Some(gas()),
        receiver: CircularOutletReceiver { position_m: [0.1, 0.0, 0.3], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        reconstruction: DuctReconstruction::FiniteDifference { sample_rate_hz: 48_000 },
    }
}

fn observation() -> DuctRadiationSpec {
    let mut s = local_observation(ValveSide::Upstream);
    let down = local_observation(ValveSide::Downstream);
    // Interleaving the groups must preserve each group's exact matrix order.
    s.openings = (0..2).flat_map(|i| [s.openings[i], DuctOpening {
        network: 1, ..down.openings[i]
    }]).collect();
    s
}

fn source(side: ValveSide, mutual: f64) -> MimoRadiationModel {
    let z = gas().characteristic_impedance / (core::f64::consts::PI * 0.007f64.powi(2));
    let d = vec![0.15 * z, mutual * z, mutual * z, 0.2 * z];
    let geometry = local_observation(side);
    let identity = match side { ValveSide::Upstream => 1, ValveSide::Downstream => 3 };
    MimoRadiationModel {
        declaration_identity: [identity; 32], realization_identity: [identity + 1; 32], exterior: gas(),
        apertures: geometry.openings.iter().zip([0.007, 0.006]).map(|(o, radius_m)|
            CircularOpening { center_m: o.center_m, radius_m }).collect(),
        band_hz: [1.0, 24_000.0],
        impedance: MimoModel { ports: 2,
            poles: vec![PoleTerm::Real { pole: -3000.0, residue: 0.0 }],
            residues: vec![d.iter().map(|v| C64::from_re(3000.0 * v)).collect()],
            d, e: vec![0.0; 4], reciprocal: true },
    }
}

fn radiation_limits() -> MimoRadiationLimits {
    MimoRadiationLimits { max_states: 32, max_memory_bytes: 1 << 20,
        sampled_passivity_points: 257 }
}
fn attach(mut system: TwoSidedAperture, mutual: f64) -> TwoSidedAperture {
    for side in SIDES {
        system = system.with_mimo_environment(side, source(side, mutual),
            local_observation(side), radiation_limits(), &CancelGate::new_clock_free()).unwrap();
    }
    system
}
fn coupled(kind: Kind, mutual: f64) -> TwoSidedAperture {
    let sides = branched_sides();
    let system = TwoSidedAperture::new(valve(kind, &sides), sides, limits(),
        &CancelGate::new_clock_free()).unwrap();
    attach(system, mutual)
}

#[test]
fn mutual_reaction_reaches_all_original_valves_and_complete_work_balance() {
    let gate = CancelGate::new_clock_free();
    for kind in KINDS {
        let mut with = coupled(kind, 0.04);
        let mut diagonal = coupled(kind, 0.0);
        let mut feedback_changed = false;
        let mut mechanics_changed = false;
        let mut work_balance = -with.stored_energy_j();
        let mut exchanged = 0.0;
        for side in SIDES {
            let original = source(side, 0.04);
            let environment = with.mimo_environment(side).unwrap();
            assert_eq!(environment.source().declaration_identity, original.declaration_identity);
            assert_eq!(environment.source().realization_identity, original.realization_identity);
            assert_eq!(environment.source().impedance.d, original.impedance.d);
            assert_eq!(environment.radiation().openings, local_observation(side).openings);
            assert_eq!(environment.source().exterior, gas());
            let (continuous, sampled) = with.mimo_passivity(side).unwrap();
            assert!(continuous.passive && sampled.passive);
            assert_eq!(sampled.class, fs_vfit::passivity::CertificateClass::GridOnly);
        }
        for n in 0..512 {
            let input = drive(if n < 128 { 50.0 } else { 0.0 });
            let a = with.step_under_gate(input, &gate).unwrap();
            let b = diagonal.step_under_gate(input, &gate).unwrap();
            assert_eq!(a.coupled_terminal_work_j,
                a.upstream.coupled_terminal_work_j + a.downstream.coupled_terminal_work_j);
            assert_eq!(a.dissipated_energy_j, a.aperture.dissipated_energy_j
                + a.upstream.terminal_loss_j + a.upstream.interior_loss_j
                + a.downstream.terminal_loss_j + a.downstream.interior_loss_j);
            assert!(a.balance_residual_j().abs() < 1e-12);
            assert!(a.upstream.balance_residual_j().abs() < 1e-12);
            assert!(a.downstream.balance_residual_j().abs() < 1e-12);
            for side in SIDES {
                for (node, port) in [2, 3].iter().zip(with.mimo_port_frames(side).unwrap()) {
                    assert_eq!(with.node_frame(side, *node).unwrap().pressure_pa, port.pressure_pa);
                    assert_eq!(with.node_frame(side, *node).unwrap().net_flow_into_node_m3_s,
                        port.flow_m3_s);
                }
            }
            exchanged += a.coupled_terminal_work_j.abs();
            work_balance += a.dissipated_energy_j + a.coupled_terminal_work_j
                - a.source_work_j - a.aperture.mechanical_work_j;
            feedback_changed |= (a.aperture.bore_pressure_pa - b.aperture.bore_pressure_pa).abs() > 1e-8;
            mechanics_changed |= (a.aperture.state.opening_m - b.aperture.state.opening_m).abs() > 1e-12;
        }
        assert!(feedback_changed && mechanics_changed && exchanged > 0.0);
        assert!((with.stored_energy_j() + work_balance).abs() < 1e-10);
    }
}

#[test]
fn shared_multirate_microphone_matches_actual_flows_and_callback_partitions() {
    let mut radiation = observation();
    radiation.reconstruction = DuctReconstruction::Differentiating {
        mechanics_rate_hz: 48_000, output_rate_hz: 24_000,
    };
    let gate = CancelGate::new_clock_free();
    for kind in KINDS {
        let mut direct = coupled(kind, 0.04);
        let mut receiver = DuctPressure::new(&[
            DuctView { graph: direct.graph(SIDES[0]), gases: direct.section_gases(SIDES[0]) },
            DuctView { graph: direct.graph(SIDES[1]), gases: direct.section_gases(SIDES[1]) },
        ], radiation.clone()).unwrap();
        let mut expected = Vec::new();
        let mut last = None;
        for _ in 0..512 {
            last = Some(direct.step_under_gate(drive(50.0), &gate).unwrap());
            if let Some(p) = receiver.push_frame(|network, node|
                direct.node_frame(SIDES[network], node).map(|f| f.net_flow_into_node_m3_s)).unwrap()
            { expected.push(p); }
        }
        assert_eq!(expected.len(), 256);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        for chunk in [1, 37, 128] {
            let mut renderer = TwoSidedPressure::with_radiation(coupled(kind, 0.04),
                radiation.clone(), 128).unwrap();
            assert_eq!(renderer.observation_policy(), ObservationPolicy::Diagnostic);
            renderer.set_drive(drive(50.0)).unwrap();
            let mut out = vec![0.0; 256];
            for block in out.chunks_mut(chunk) { renderer.block(block).unwrap(); }
            assert_eq!(out, expected);
            assert_eq!(renderer.last_frame(), last.as_ref());
            assert_eq!(renderer.mechanical_samples_rendered(), 512);
            for side in SIDES {
                assert_eq!(renderer.system().mimo_pole_states(side), direct.mimo_pole_states(side));
            }
        }
    }
}

#[test]
fn global_receiver_cannot_relabel_omit_or_reorder_either_local_source() {
    let mut bad = Vec::new();
    let mut changed = observation();
    changed.exterior = Some(GasState::try_new(&GasSpec::dry_air_ussa1976(), 303.15, 101325.0).unwrap());
    bad.push(changed);
    let mut changed = observation(); changed.openings.swap(0, 2); bad.push(changed);
    let mut changed = observation(); changed.openings.swap(1, 3); bad.push(changed);
    let mut changed = observation(); changed.openings.remove(1); bad.push(changed);
    let mut changed = observation(); changed.openings[1].center_m[0] += 0.01; bad.push(changed);
    let mut changed = observation();
    changed.reconstruction = DuctReconstruction::FiniteDifference { sample_rate_hz: 96_000 };
    bad.push(changed);
    for radiation in bad {
        assert!(TwoSidedPressure::with_radiation(coupled(Kind::Paired, 0.04), radiation, 64).is_err());
    }
    assert!(TwoSidedPressure::with_observation_policy(coupled(Kind::Paired, 0.04),
        observation(), 64, ObservationPolicy::Acceptance, &[]).is_err());
}

#[test]
fn attachment_keeps_regional_gases_fractional_travel_and_physical_loss_maps() {
    let mut sides = branched_sides();
    for side in &mut sides {
        for (section, gas) in side.graph.sections.iter_mut().zip(&side.gases) {
            section.length_m += 0.25 * gas.sound_speed * DT;
        }
        side.viscothermal = vec![ViscothermalSectionSpec { section: 2, cells: 2,
            minimum_frequency_hz: 200.0, maximum_frequency_hz: 400.0 }];
    }
    let policy = FractionalTubePolicy { maximum_frequency_hz: 1000.0,
        maximum_summed_phase_error_rad: 0.1 };
    let original = sides.clone();
    let system = TwoSidedAperture::with_fractional_transits(valve(Kind::Modal, &sides), sides,
        limits(), [policy; 2], &CancelGate::new_clock_free()).unwrap();
    let graphs = SIDES.map(|s| system.graph(s).clone());
    let gases = SIDES.map(|s| system.section_gases(s).to_vec());
    let represented = SIDES.map(|s| system.represented_sections(s).to_vec());
    let fractional = SIDES.map(|s| system.fractional_sections(s).unwrap().1.to_vec());
    let mut attached = attach(system, 0.04);
    for (i, side) in SIDES.into_iter().enumerate() {
        assert_eq!(attached.source(side).graph, original[i].graph);
        assert_eq!(attached.source(side).gases, original[i].gases);
        assert_eq!(attached.source(side).viscothermal, original[i].viscothermal);
        assert_eq!(attached.graph(side), &graphs[i]);
        assert_eq!(attached.section_gases(side), gases[i]);
        assert_eq!(attached.represented_sections(side), represented[i]);
        assert_eq!(attached.fractional_sections(side).unwrap(), (policy, fractional[i].as_slice()));
        assert_eq!(attached.viscothermal_sections(side).len(), 1);
        assert_eq!(attached.viscothermal_sections(side)[0].source, original[i].viscothermal[0]);
    }
    let gate = CancelGate::new_clock_free();
    for _ in 0..128 {
        let f = attached.step_under_gate(drive(20.0), &gate).unwrap();
        assert!(f.balance_residual_j().abs() < 1e-12);
    }
}

#[test]
fn refused_samples_and_terminal_changes_preserve_all_original_histories() {
    let gate = CancelGate::new_clock_free();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    for kind in KINDS {
        let mut a = coupled(kind, 0.04);
        let mut b = coupled(kind, 0.04);
        for _ in 0..128 {
            assert_eq!(a.step_under_gate(drive(50.0), &gate).unwrap(),
                b.step_under_gate(drive(50.0), &gate).unwrap());
        }
        for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(a.step_under_gate(drive(bad), &gate).is_err());
        }
        assert!(a.step_under_gate(drive(50.0), &cancelled).is_err());
        for side in SIDES {
            assert!(a.set_terminal_reflection(side, 2, -0.5).is_err());
            assert_eq!(a.mimo_pole_states(side), b.mimo_pole_states(side));
            assert!(a.mimo_pole_states(side).unwrap().iter().any(|s| s.re != 0.0 || s.im != 0.0));
            assert_eq!(a.mimo_port_frames(side), b.mimo_port_frames(side));
            for node in 0..4 { assert_eq!(a.node_frame(side, node), b.node_frame(side, node)); }
        }
        assert_eq!(a.valve().accepted_steps(), 128);
        assert_eq!(a.stored_energy_j(), b.stored_energy_j());
        for _ in 0..128 {
            assert_eq!(a.step_under_gate(drive(0.0), &gate).unwrap(),
                b.step_under_gate(drive(0.0), &gate).unwrap());
        }
    }
}

#[test]
fn whole_system_energy_refusal_does_not_publish_either_rational_candidate() {
    let gate = CancelGate::new_clock_free();
    let mut reference = coupled(Kind::Paired, 0.04);
    let mut probe = coupled(Kind::Paired, 0.04);
    let mut peak = 0.0_f64;
    for _ in 0..128 {
        let frame = reference.step_under_gate(drive(1.0), &gate).unwrap();
        assert_eq!(frame, probe.step_under_gate(drive(1.0), &gate).unwrap());
        peak = peak.max(frame.stored_energy_j);
    }
    let next = probe.step_under_gate(drive(50.0), &gate).unwrap();
    assert!(next.stored_energy_j > peak);
    let mut budget = limits();
    budget.maximum_energy_j = 0.5 * (peak + next.stored_energy_j);
    let sides = branched_sides();
    let mut limited = attach(TwoSidedAperture::new(valve(Kind::Paired, &sides), sides,
        budget, &gate).unwrap(), 0.04);
    for _ in 0..128 { limited.step_under_gate(drive(1.0), &gate).unwrap(); }
    let error = limited.step_under_gate(drive(50.0), &gate).unwrap_err();
    assert!(matches!(error, fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription {
        what: "two-sided valve refused the complete storage, work or flow balance",
    }));
    assert_eq!(limited.valve().accepted_steps(), reference.valve().accepted_steps());
    assert_eq!(limited.stored_energy_j(), reference.stored_energy_j());
    for side in SIDES {
        assert_eq!(limited.mimo_pole_states(side), reference.mimo_pole_states(side));
        assert_eq!(limited.mimo_port_frames(side), reference.mimo_port_frames(side));
        for node in 0..4 { assert_eq!(limited.node_frame(side, node), reference.node_frame(side, node)); }
    }
    for _ in 0..64 {
        assert_eq!(limited.step_under_gate(drive(0.0), &gate).unwrap(),
            reference.step_under_gate(drive(0.0), &gate).unwrap());
    }
}

#[test]
fn a_second_or_late_group_and_mismatched_source_geometry_refuse() {
    let gate = CancelGate::new_clock_free();
    let side = ValveSide::Downstream;
    assert!(coupled(Kind::Paired, 0.04).with_mimo_environment(side, source(side, 0.04),
        local_observation(side), radiation_limits(), &gate).is_err());
    let sides = branched_sides();
    let new = || TwoSidedAperture::new(valve(Kind::Paired, &sides), sides.clone(), limits(), &gate).unwrap();
    let mut advanced = new(); advanced.step_under_gate(drive(50.0), &gate).unwrap();
    assert!(advanced.with_mimo_environment(side, source(side, 0.04),
        local_observation(side), radiation_limits(), &gate).is_err());
    let mut wrong = source(side, 0.04); wrong.apertures[0].radius_m *= 1.1;
    assert!(new().with_mimo_environment(side, wrong,
        local_observation(side), radiation_limits(), &gate).is_err());
    let mut wrong = local_observation(side); wrong.openings[0].network = 1;
    assert!(new().with_mimo_environment(side, source(side, 0.04), wrong,
        radiation_limits(), &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(new().with_mimo_environment(side, source(side, 0.04), local_observation(side),
        radiation_limits(), &cancelled).is_err());
}

mod keyed;
