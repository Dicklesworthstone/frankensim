//! Shared physical loss/propagation integration; synthetic geometry, not measured timbre.
use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::viscothermal::{
    ViscothermalSectionSpec, with_regional_fractional_viscothermal_sections,
    with_regional_viscothermal_sections, with_viscothermal_sections,
};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_vfit::waveguide::fractional::FractionalTransit;

const DT: f64 = 1.0 / 48_000.0;
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
fn policy() -> FractionalTubePolicy {
    FractionalTubePolicy { maximum_frequency_hz: 400.0,
        maximum_summed_phase_error_rad: 0.01 }
}
fn selection(section: usize) -> ViscothermalSectionSpec {
    ViscothermalSectionSpec { section, cells: 8,
        minimum_frequency_hz: 200.0, maximum_frequency_hz: 400.0 }
}
fn graph(g: GasState, length: f64) -> TubeNetworkSpec {
    TubeNetworkSpec { nodes: vec![NetworkNode::Inlet,
        NetworkNode::Termination { reflection: -0.8 }],
        sections: vec![TubeSection { nodes: [0, 1], length_m: length,
            radius_m: 0.007, max_length_error_m: 1e-12 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 8 << 20 }
}

#[test]
fn phase_queries_match_the_actual_complex_allpass_not_an_ideal_delay() {
    for d in [0.0, 0.01, 0.25, 0.5, 0.99] {
        let transit = FractionalTransit::new(d).unwrap();
        let a = (1.0-d)/(1.0+d);
        for i in 0..=128 {
            let w = 3.0*f64::from(i)/128.0;
            // Direct complex division of (a+exp(-iw))/(1+a*exp(-iw)).
            let (s,c) = w.sin_cos();
            let den = (1.0+a*c).powi(2)+(a*s).powi(2);
            let re = ((a+c)*(1.0+a*c)+a*s*s)/den;
            let im = (-s*(1.0+a*c)+(a+c)*a*s)/den;
            let lag = transit.phase_lag_rad(w).unwrap();
            assert!((lag+im.atan2(re)).abs() < 2e-13);
            assert!((re*re+im*im-1.0).abs() < 2e-13);
            assert!((lag-d*w).abs() <= transit.phase_error_bound_rad(w).unwrap()+2e-15);
        }
    }
    assert!(FractionalTransit::new(0.5).unwrap().phase_lag_rad(f64::NAN).is_err());
    assert!(FractionalTransit::new(0.5).unwrap().phase_lag_rad(core::f64::consts::PI).is_err());
}

#[test]
fn loss_subdivision_preserves_subsample_geometry_and_its_actual_wave_response() {
    let g = gas(293.15);
    let gate = CancelGate::new_clock_free();
    let mut outputs = Vec::new();
    for length in [0.4001, 0.4002] {
        let original = graph(g, length);
        let (lowered, gases, reports) = with_regional_fractional_viscothermal_sections(
            original, vec![g], DT, &[selection(0)], policy(), &gate).unwrap();
        assert_eq!(lowered.nodes[0], NetworkNode::Inlet);
        assert_eq!(lowered.nodes[1], NetworkNode::Termination { reflection: -0.8 });
        assert_eq!(reports[0].original_nodes, [0,1]);
        assert_eq!(reports[0].section_range, [0,32]);
        assert_eq!(reports[0].requested_length_m, length);
        assert!((reports[0].represented_length_m-length).abs() < 1e-12);
        assert!((lowered.sections.iter().map(|s|s.length_m).sum::<f64>()-length).abs() < 1e-14);
        assert!(reports[0].max_scattering_error <= 0.03);
        let (mut network, represented, fractional) = lowered.realize_fractional_waveguide(
            &gases, DT, policy(), &gate).unwrap();
        assert!((represented.iter().map(|s|s.represented_length_m).sum::<f64>()-length).abs() < 1e-12);
        assert!(fractional.iter().any(|s|s.transit.samples()>0.0));
        let mut trace = Vec::new();
        let mut work = 0.0;
        let mut loss = 0.0;
        for k in 0..2048 {
            let f = network.step(if k < 24 { 0.1 } else { 0.0 }).unwrap();
            assert!(f.balance_residual_j().abs() < 1e-18);
            assert!(f.stored_energy_j >= 0.0);
            work += f.inlet_work_j;
            loss += f.terminal_loss_j+f.interior_loss_j;
            trace.push(network.node_frame(1).unwrap().pressure_pa);
        }
        assert!(loss > 0.0);
        assert!((network.stored_energy_j()+loss-work).abs() < 1e-17);
        assert!(trace.iter().any(|p|p.abs()>1e-6));
        outputs.push(trace);
    }
    assert!(outputs[0].iter().zip(&outputs[1]).any(|(a,b)|(a-b).abs()>1e-7));
}

#[test]
fn regional_temperature_changes_transit_and_loss_without_changing_physical_geometry() {
    let inlet = gas(293.15);
    let gate = CancelGate::new_clock_free();
    let mut snapshots = Vec::new();
    for t in [293.15, 313.15] {
        let remote = gas(t);
        let mut original = graph(inlet, 0.4001);
        original.nodes.push(NetworkNode::Termination { reflection: -0.7 });
        original.nodes[1] = NetworkNode::Junction;
        original.sections.push(TubeSection { nodes: [1,2], length_m: 0.4001,
            radius_m: 0.007, max_length_error_m: 1e-12 });
        let (lowered, gases, reports) = with_regional_fractional_viscothermal_sections(
            original, vec![inlet, remote], DT, &[selection(0),selection(1)], policy(), &gate).unwrap();
        let (_, _, transit) = lowered.realize_fractional_waveguide(&gases, DT, policy(), &gate).unwrap();
        assert_eq!(reports[1].original_nodes, [1,2]);
        assert_eq!(gases[0].sound_speed, inlet.sound_speed);
        assert_eq!(gases[32].sound_speed, remote.sound_speed);
        let physical: Vec<_> = lowered.sections.iter().map(|s|s.length_m).collect();
        let samples: f64 = transit[32..].iter().map(|s|s.requested_samples).sum();
        assert!((samples-0.4001/(remote.sound_speed*DT)).abs() < 1e-12);
        let viscosity = reports[1].loss.dynamic_viscosity_pa_s();
        assert_eq!(viscosity, remote.dynamic_viscosity);
        snapshots.push((physical,samples,viscosity));
    }
    assert_eq!(snapshots[0].0, snapshots[1].0);
    assert!(snapshots[1].1 < snapshots[0].1);
    assert_ne!(snapshots[0].2, snapshots[1].2);
}

#[test]
fn legacy_integer_image_and_empty_selection_keep_their_original_paths() {
    let g = gas(293.15);
    let gate = CancelGate::new_clock_free();
    let mut original = graph(g,0.4001);
    original.sections[0].max_length_error_m = 0.01;
    let (uniform, a) = with_viscothermal_sections(original.clone(), &g, DT, &[selection(0)], &gate).unwrap();
    let (regional, _, b) = with_regional_viscothermal_sections(original.clone(), vec![g], DT, &[selection(0)], &gate).unwrap();
    assert_eq!(uniform.nodes, regional.nodes);
    assert_eq!(a[0].represented_length_m, b[0].represented_length_m);
    assert_eq!(a[0].max_scattering_error, b[0].max_scattering_error);
    for (a,b) in uniform.sections.iter().zip(&regional.sections) {
        assert_eq!(a.length_m,b.length_m);
    }
    assert!((a[0].represented_length_m-0.4001).abs()>1e-5);
    let (empty, gases, reports) = with_regional_fractional_viscothermal_sections(
        original.clone(),vec![g],DT,&[],policy(),&gate).unwrap();
    assert!(reports.is_empty());
    assert_eq!(empty.nodes,original.nodes);
    assert_eq!(empty.sections[0].length_m,original.sections[0].length_m);
    assert_eq!(gases.len(),1);
}

#[test]
fn unresolved_geometry_phase_band_and_memory_are_refused_without_relaxing_gates() {
    let g=gas(293.15);
    let gate=CancelGate::new_clock_free();
    let try_build=|graph, specs: &[ViscothermalSectionSpec], p| {
        with_regional_fractional_viscothermal_sections(graph,vec![g],DT,specs,p,&gate)
    };
    let mut p=policy(); p.maximum_summed_phase_error_rad=0.0;
    assert!(try_build(graph(g,0.4001),&[selection(0)],p).is_err());
    p=policy();p.maximum_frequency_hz=300.0;
    assert!(try_build(graph(g,0.4001),&[selection(0)],p).is_err());
    assert!(try_build(graph(g,0.1),&[selection(0)],policy()).is_err());
    assert!(try_build(graph(g,0.4001),&[selection(0),selection(0)],policy()).is_err());
    let mut small=graph(g,0.4001);small.max_wave_memory_bytes=1;
    assert!(try_build(small,&[selection(0)],policy()).is_err());
    let mut s=selection(0);s.cells=0;
    assert!(try_build(graph(g,0.4001),&[s],policy()).is_err());
    p=policy();p.maximum_frequency_hz=f64::NAN;
    assert!(try_build(graph(g,0.4001),&[],p).is_err());
}

#[path = "fractional_viscothermal/reeds.rs"]
mod reeds;
