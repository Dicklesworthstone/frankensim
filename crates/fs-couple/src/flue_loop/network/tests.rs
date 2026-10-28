use super::*;
use crate::flue_loop::FlueControl;
use crate::bernoulli_aperture::network::{NetworkNode, TubeSection};
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_material::gas::GasSpec;

mod mimo;

pub(super) fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
pub(super) fn geometry() -> FlueGeometry {
    FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m: 0.004,
        labium_offset_m: 0.0, jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * 0.006 * 0.006, mouth_loss_rel: 0.02 }
}
pub(super) fn source(gas: GasState, reflection: f64) -> FlueNetworkSpec {
    FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
            sections: vec![TubeSection { nodes: [0, 1], length_m: 0.3, radius_m: 0.006,
                max_length_error_m: 0.005 }],
            sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3, max_jet_history_samples: 4096,
    }
}
pub(super) fn voice(source: FlueNetworkSpec) -> FlueVoice {
    let card = mint_tonal_interim_card().unwrap();
    FlueVoice::from_network(Some(&card), 90.0, 1.0e-6, geometry(), source,
        48_000, 3, 0, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn jet_and_whole_tube_network_close_the_same_work_window() {
    let mut v = voice(source(gas(), -0.9));
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut launched = vec![0.0; 2048];
    let mut input = 0.0;
    let mut loss = 0.0;
    v.step_block_observed(&mut launched, |line| {
        let FlueResonator::Network(n) = line else { panic!("physical graph") };
        let f = n.last_frame().unwrap();
        assert!(f.balance_residual_j().abs() < 1.0e-12);
        input += f.inlet_work_j;
        loss += f.terminal_loss_j + f.interior_loss_j;
        Ok(())
    }).unwrap();
    let d = v.diagnostics().last().unwrap();
    let stored = v.network().unwrap().stored_energy_j();
    let scale = 1.0 + d.source_work_j.abs();
    assert!((input - d.pipe_work_j).abs() < 1.0e-12 * scale);
    assert!((stored + loss - input).abs() < 1.0e-12 * scale);
    assert!((d.source_work_j - stored - loss - d.stored_delta_j
        - d.numerical_dissipation_j - d.radiation_loss_j).abs() < 1.0e-12 * scale);
    assert!(launched.iter().any(|p| p.abs() > 1.0e-6));
    assert_eq!(v.samples_rendered(), 2048);
}

#[test]
fn the_actual_outlet_flow_arrives_after_geometry_transit_and_reacts_back() {
    let mut a = voice(source(gas(), -0.9));
    let mut b = voice(source(gas(), 0.0));
    let delay = a.network().unwrap().represented_sections()[0].one_way_samples;
    for v in [&mut a, &mut b] { v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap(); }
    let mut pa = vec![0.0; 1024];
    let mut pb = vec![0.0; 1024];
    let mut flows = Vec::new();
    a.step_block_observed(&mut pa, |line| {
        let FlueResonator::Network(n) = line else { panic!("physical graph") };
        flows.push(n.node_frame(1).unwrap().net_flow_into_node_m3_s);
        Ok(())
    }).unwrap();
    b.step_block(&mut pb).unwrap();
    assert!(flows[..delay].iter().all(|q| *q == 0.0));
    assert!(flows[delay].abs() > 0.0);
    // A far-end load cannot affect the source before a complete round trip.
    assert_eq!(&pa[..2*delay], &pb[..2*delay]);
    assert!(pa[2*delay..].iter().zip(&pb[2*delay..]).any(|(x,y)| (x-y).abs() > 1.0e-8));
    let energy = a.network().unwrap().stored_energy_j();
    a.set_network_terminal_reflection(1, 0.2).unwrap();
    assert_eq!(a.network().unwrap().stored_energy_j(), energy);
    assert!(a.apply(FlueControl::Fingering { index: 0, fade_samples: 32 }).is_err());
}

#[test]
fn data_only_gas_substitution_changes_jet_impedance_timing_and_actual_waveform() {
    let g = gas();
    // Synthetic gas parameters test causal substitution, not measured gas fidelity.
    let alternate = GasState::try_new(&GasSpec {
        molar_mass: 0.020, gamma: 1.5,
        conductivity: fs_material::gas::ConductivityModel::Eucken,
        ..GasSpec::dry_air_ussa1976()
    }, 293.15, 101_325.0).unwrap();
    let mut a = voice(source(g, -0.9));
    let mut b = voice(source(alternate, -0.9));
    for v in [&mut a, &mut b] { v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap(); }
    assert_ne!(a.jet_speed_m_s(), b.jet_speed_m_s());
    assert_ne!(a.reynolds(), b.reynolds());
    assert_ne!(a.transit_s(), b.transit_s());
    assert_ne!(a.network().unwrap().represented_sections()[0], b.network().unwrap().represented_sections()[0]);
    let mut pa = [0.0; 512];
    let mut pb = [0.0; 512];
    a.step_block(&mut pa).unwrap();
    b.step_block(&mut pb).unwrap();
    assert_ne!(pa, pb);
    assert_eq!(b.network().unwrap().section_gases(), &[alternate]);
}

#[test]
fn selected_material_viscothermal_losses_join_the_same_jet_network() {
    let mut spec = source(gas(), -0.9);
    spec.viscothermal = vec![ViscothermalSectionSpec {
        section: 0, cells: 8, minimum_frequency_hz: 300.0, maximum_frequency_hz: 500.0,
    }];
    let mut v = voice(spec);
    let n = v.network().unwrap();
    assert_eq!(n.viscothermal_sections().len(), 1);
    assert!(n.graph().sections.len() > 1);
    assert_eq!(n.viscothermal_sections()[0].loss.dynamic_viscosity_pa_s(), gas().dynamic_viscosity);
    let represented = n.viscothermal_sections()[0].one_way_samples;
    assert_eq!(n.represented_sections().iter().map(|s| s.one_way_samples).sum::<usize>(), represented);
    v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
    let mut wave = [0.0; 512];
    v.step_block(&mut wave).unwrap();
    let f = v.network().unwrap().last_frame().unwrap();
    assert!(f.interior_loss_j > 0.0);
    assert!(f.balance_residual_j().abs() < 1.0e-12);
}

#[test]
fn network_preserves_block_partition_and_rejects_stale_geometry_or_jet_limits() {
    let mut one = voice(source(gas(), -0.9));
    let mut many = voice(source(gas(), -0.9));
    for v in [&mut one, &mut many] { v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap(); }
    let mut a = [0.0; 513];
    let mut b = [0.0; 513];
    one.step_block(&mut a).unwrap();
    for block in b.chunks_mut(37) { many.step_block(block).unwrap(); }
    assert_eq!(a, b);
    assert_eq!(one.network().unwrap().last_frame(), many.network().unwrap().last_frame());
    let old = one.jet_speed_m_s();
    assert!(one.apply(FlueControl::SetBlowingPressure(1.0e9)).is_err());
    assert_eq!(one.jet_speed_m_s(), old);
    // Ruling 8132: a jet below the declared minimum speed is admitted. Past the onset, its
    // lookup beyond the retained history takes the declared state instead of refusing.
    one.apply(FlueControl::SetBlowingPressure(1.0e-10)).unwrap();
    one.step_block(&mut a[..1]).unwrap();
    assert_eq!(one.diagnostics().last().unwrap().below_speed_edge_samples, 1);
    one.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
    one.step_block(&mut a[..1]).unwrap();
    assert!(one.network().unwrap().stored_energy_j() > 0.0);
    let card = mint_tonal_interim_card().unwrap();
    let mut geom = geometry();
    geom.bore_area_m2 *= 2.0;
    assert!(FlueVoice::from_network(Some(&card), 90.0, 0.0, geom, source(gas(), -0.9),
        48_000, 3, 0, &CancelGate::new_clock_free()).is_err());
}
