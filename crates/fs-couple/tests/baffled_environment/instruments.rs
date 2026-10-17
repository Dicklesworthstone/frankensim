use super::*;
use fs_couple::bernoulli_aperture::paired::{
    ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits,
    PairedApertureSpec,
};
use fs_couple::bernoulli_aperture::two_sided::{
    ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive,
    TwoSidedLimits, ValveSide, combined_inlet_state,
};
use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::baffled::openings::duct::{DuctPressure,DuctView};
use fs_dcontact::Obstacle;

const STEPS:u64=1024;
fn limits()->TwoSidedLimits {
    TwoSidedLimits {maximum_energy_j:1.0,energy_absolute_tolerance_j:1e-16,
        energy_relative_tolerance:1e-7,flow_absolute_tolerance_m3_s:1e-16,
        flow_relative_tolerance:1e-8}
}
fn reed()->(LoadedAperture,[ApertureDuctSpec;2]) {
    let g=gas(293.15);
    let sources=std::array::from_fn(|_|ApertureDuctSpec {graph:graph(g,0.007),
        gases:vec![g;2],viscothermal:vec![]});
    let (_,z)=combined_inlet_state(&sources).unwrap();
    let b=ApertureBlade {mass_kg:1e-5,stiffness_n_m:500.0,damping_n_s_m:0.03,pressure_area_m2:1.25e-4};
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic exterior-coupled reed fixture".into()).unwrap();
    // Outside the midpoint one-cent band at this rate (DEF-04): exact blades.
    let valve=PairedAperture::new_with_propagator(PairedApertureSpec {rest_opening_m:4e-4,width_m:0.013,
        blades:[b,ApertureBlade {mass_kg:1.3e-5,..b}],density_kg_m3:g.density,
        impedance_pa_s_m3:z,time_step_s:1.0/f64::from(RATE),max_steps:STEPS,
        limits:PairedApertureLimits {maximum_abs_displacement_m:0.01,
            maximum_abs_velocity_m_s:100.0,maximum_energy_j:1.0,max_evaluations:256,
            velocity_tolerance_m_s:1e-13,relative_tolerance:1e-11}},
        [BladeState::default();2],contact,BladePropagator::ExactZoh).unwrap();
    (valve.into(),sources)
}
fn observation(exterior:GasState,ratio:u32)->DuctRadiationSpec {
    let mut input=radiation(Some(exterior),1);
    input.openings[0].network=1;
    input.reconstruction=DuctReconstruction::Differentiating {
        mechanics_rate_hz:RATE,output_rate_hz:RATE/ratio};
    input
}
fn performer(exterior:GasState,ratio:u32)->TwoSidedPressure {
    let (v,s)=reed();
    TwoSidedPressure::from_baffled_environment(v,s,[None,None],limits(),None,
        observation(exterior,ratio),64,&CancelGate::new_clock_free()).unwrap().0
}
fn drive(p:f64)->TwoSidedDrive {
    TwoSidedDrive {source_pressure_pa:p,blade_forces_n:[0.0;2]}
}
fn blades(p:&TwoSidedPressure)->[BladeState;2] {
    let LoadedAperture::Paired(a)=p.system().valve() else {panic!("paired fixture")};
    *a.blade_states()
}

/// G3: the actual accepted terminal flow, not an independently rebuilt source.
#[test]
fn two_sided_pressure_matches_owner() {
    use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium};
    let exterior = gas(333.15);
    let mut actual = performer(exterior, 1);
    let declaration = actual.radiation();
    let port = declaration.openings[0];
    let side = if port.network == 0 {
        ValveSide::Upstream
    } else {
        ValveSide::Downstream
    };
    let section = actual
        .system()
        .graph(side)
        .sections
        .iter()
        .find(|s| s.nodes.contains(&port.node))
        .unwrap();
    let mut owner = BaffledPressure::circular_outlet_multirate_differentiating(
        section.radius_m,
        RATE,
        RATE,
        declaration.receiver,
        RayleighMedium {
            density: exterior.density,
            sound_speed: exterior.sound_speed,
        },
    )
    .unwrap();
    let mut peak = 0.0_f64;
    let start = std::time::Instant::now();
    for sample in 0..STEPS {
        actual
            .set_drive(drive(if sample < 173 { 50.0 } else { 0.0 }))
            .unwrap();
        let mut out = [0.0];
        actual.block(&mut out).unwrap();
        let flow = actual
            .system()
            .node_frame(side, port.node)
            .unwrap()
            .net_flow_into_node_m3_s;
        assert_eq!(out[0].to_bits(), owner.step(&[flow]).unwrap().to_bits());
        peak = peak.max(out[0].abs());
    }
    assert!(peak > 0.0);
    println!(
        "{{\"case\":\"two_sided_pressure_matches_owner\",\"inputs\":\"authored paired valve and two physical gas-bound ducts; terminal volume flow\",\"observed\":0,\"err\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
        start.elapsed().as_secs_f64()
    );
}

/// G0: load provenance must agree with receiver media before publication.
#[test]
fn reaction_receiver_gas_mismatch_refuses() {
    use fs_couple::pcm_wav::observation::ObservationPolicy;
    let (valve, [up, down]) = reed();
    let mut graphs = [up.graph, down.graph];
    let environment = bind_baffled_environment(
        &mut graphs,
        &[&up.gases, &down.gases],
        observation(gas(293.15), 1),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let (mut receiver, bindings) = environment.into_parts();
    receiver.exterior = Some(gas(333.15));
    let [a, b] = graphs;
    let system = TwoSidedAperture::new(
        valve,
        [
            ApertureDuctSpec { graph: a, ..up },
            ApertureDuctSpec { graph: b, ..down },
        ],
        limits(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let error = TwoSidedPressure::with_observation_policy(
        system,
        receiver,
        64,
        ObservationPolicy::Acceptance,
        &bindings,
    )
    .err()
    .unwrap();
    assert!(
        error.to_string().contains("reaction/receiver gas mismatch"),
        "{error}"
    );
}

/// G5: reconstruction, mechanics and control history replay across callback sizes.
#[test]
fn observer_bitwise_repeat() {
    let (mut a, mut b) = (performer(gas(333.15), 1), performer(gas(333.15), 1));
    a.set_drive(drive(50.0)).unwrap();
    b.set_drive(drive(50.0)).unwrap();
    let (mut x, mut y) = ([0.0; 1024], [0.0; 1024]);
    for part in x.chunks_mut(64) {
        a.block(part).unwrap();
    }
    for part in y.chunks_mut(7) {
        b.block(part).unwrap();
    }
    assert!(x.iter().any(|v| *v != 0.0));
    assert!(x.iter().zip(y).all(|(a, b)| a.to_bits() == b.to_bits()));
}

#[test]
fn one_environment_matches_manual_original_owners_at_each_multirate_sample() {
    for ratio in [1,2,4] {
        let ext=gas(333.15);let mut actual=performer(ext,ratio);
        let (v,mut sources)=reed();
        let load=BaffledRadiationLoad::new(0.007,ext.density,ext.sound_speed,
            1.0/f64::from(RATE),400.0,&CancelGate::new_clock_free()).unwrap();
        sources[1].graph.nodes[2]=load.termination();
        let model=TwoSidedAperture::new(v,sources,limits(),&CancelGate::new_clock_free()).unwrap();
        let mut expected=TwoSidedPressure::with_radiation(model,observation(ext,ratio),64).unwrap();
        for n in 0..STEPS/u64::from(ratio) {
            let held=drive(if n<64 {50.0}else{0.0});
            actual.set_drive(held).unwrap();expected.set_drive(held).unwrap();
            let (mut a,mut b)=([0.0],[0.0]);actual.block(&mut a).unwrap();expected.block(&mut b).unwrap();
            assert_eq!(a,b);assert_eq!(actual.last_frame(),expected.last_frame());
            assert_eq!(blades(&actual),blades(&expected));
        }
        assert!(actual.system().stored_energy_j()>0.0);
        assert_eq!(actual.mechanical_samples_rendered(),STEPS);
    }
}

#[test]
fn exterior_only_changes_reed_mechanics_after_return_not_just_the_microphone() {
    let mut cold=performer(gas(293.15),1);let mut warm=performer(gas(333.15),1);
    let delay=2*cold.system().represented_sections(ValveSide::Downstream)
        .iter().map(|s|s.one_way_samples).sum::<usize>();
    let mut changed=false;let mut heard=false;
    for n in 0..STEPS {
        for p in [&mut cold,&mut warm] {p.set_drive(drive(if n<173 {50.0}else{0.0})).unwrap();}
        let (mut a,mut b)=([0.0],[0.0]);cold.block(&mut a).unwrap();warm.block(&mut b).unwrap();
        if (n as usize)<delay {assert_eq!(blades(&cold),blades(&warm));}
        else {changed|=blades(&cold)!=blades(&warm);}
        heard|=a!=b;
        assert_eq!(cold.system().section_gases(ValveSide::Downstream),warm.system().section_gases(ValveSide::Downstream));
        for p in [&cold,&warm] {assert!(p.last_frame().unwrap().balance_residual_j().abs()<1e-14);}
    }
    assert!(changed && heard);assert!(cold.system().stored_energy_j()>0.0);
}

#[test]
fn jet_source_and_receiver_share_the_exterior_while_internal_gas_and_card_stay_fixed() {
    use fs_couple::flue_loop::{FlueControl,FlueGeometry,FlueVoice};
    use fs_couple::flue_loop::network::FlueNetworkSpec;
    let g=gas(293.15);let card=fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    let geometry=FlueGeometry {flue_width_m:0.010,flue_height_m:0.001,cut_up_m:0.004,
        labium_offset_m:0.0,jet_angle_rad:0.0,bore_area_m2:std::f64::consts::PI*0.007*0.007,mouth_loss_rel:0.02};
    let make=|exterior| {
        let source=FlueNetworkSpec {graph:graph(g,0.007),gases:vec![g;2],viscothermal:vec![],
            minimum_jet_speed_m_s:1.0,maximum_jet_mach:0.3,max_jet_history_samples:4096};
        let (voice,environment)=FlueVoice::from_baffled_environment(Some(&card),90.0,1e-6,
            geometry,source,None,None,radiation(Some(exterior),1),
            3,0,&CancelGate::new_clock_free()).unwrap();
        let net=voice.network().unwrap();
        let mic=DuctPressure::new(&[DuctView {graph:net.graph(),gases:net.section_gases()}],
            environment.radiation().clone()).unwrap();
        (voice,mic)
    };
    let (mut a,mut ma)=make(gas(293.15));let (mut b,mut mb)=make(gas(333.15));
    let delay=2*a.network().unwrap().represented_sections().iter().map(|s|s.one_way_samples).sum::<usize>();
    let mut feedback=false;let mut heard=false;
    for n in 0..STEPS {
        for v in [&mut a,&mut b] {v.apply(FlueControl::SetBlowingPressure(if n<256 {300.0}else{0.0})).unwrap();}
        let (mut x,mut y)=([0.0],[0.0]);a.step_block(&mut x).unwrap();b.step_block(&mut y).unwrap();
        if (n as usize)<delay {assert_eq!(x,y);}else{feedback|=x!=y;}
        let pa=ma.push_frame(|_,node|a.network().unwrap().node_frame(node).map(|f|f.net_flow_into_node_m3_s)).unwrap();
        let pb=mb.push_frame(|_,node|b.network().unwrap().node_frame(node).map(|f|f.net_flow_into_node_m3_s)).unwrap();
        heard|=pa!=pb;
        assert_eq!(a.network().unwrap().section_gases(),b.network().unwrap().section_gases());
        assert_eq!(a.island(),b.island());
    }
    assert!(feedback && heard);assert!(a.network().unwrap().stored_energy_j()>0.0);
}

#[test]
fn coupled_environment_rejects_clock_mismatch_and_retains_callback_independence() {
    let (v,s)=reed();let mut wrong=observation(gas(293.15),1);
    wrong.reconstruction=DuctReconstruction::FiniteDifference {sample_rate_hz:RATE/2};
    assert!(TwoSidedPressure::from_baffled_environment(v,s,[None,None],limits(),None,wrong,64,
        &CancelGate::new_clock_free()).is_err());
    let mut a=performer(gas(333.15),2);let mut b=performer(gas(333.15),2);
    a.set_drive(drive(50.0)).unwrap();b.set_drive(drive(50.0)).unwrap();
    let mut expected=[0.0;128];let mut actual=expected;
    for p in &mut expected {a.block(core::slice::from_mut(p)).unwrap();}
    for chunk in actual.chunks_mut(31) {b.block(chunk).unwrap();}
    assert_eq!(actual,expected);assert_eq!(a.last_frame(),b.last_frame());
    assert_eq!(blades(&a),blades(&b));
}

#[test]
fn fractional_reed_selection_preserves_original_geometry_and_declared_radiation_band() {
    use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
    let (v,s)=reed();let original=format!("{:?}",s[1].graph.sections);
    let policy=FractionalTubePolicy {maximum_frequency_hz:400.0,maximum_summed_phase_error_rad:0.01};
    let (mut p,loads)=TwoSidedPressure::from_baffled_environment(v,s,[None,None],limits(),
        Some([policy;2]),observation(gas(333.15),2),64,&CancelGate::new_clock_free()).unwrap();
    assert_eq!(format!("{:?}",p.system().source(ValveSide::Downstream).graph.sections),original);
    assert_eq!(loads[0].load().maximum_frequency_hz(),400.0);
    p.set_drive(drive(50.0)).unwrap();p.block(&mut [0.0;64]).unwrap();
    assert_eq!(p.mechanical_samples_rendered(),128);
    let (v,s)=reed();let narrow=FractionalTubePolicy {maximum_frequency_hz:200.0,..policy};
    assert!(TwoSidedPressure::from_baffled_environment(v,s,[None,None],limits(),Some([narrow;2]),
        observation(gas(333.15),2),64,&CancelGate::new_clock_free()).is_err());
}
