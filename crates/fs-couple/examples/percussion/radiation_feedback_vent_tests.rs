//! Actual drum/neck/BEM integration. Declared passive loads below are explicit
//! mechanics oracles, never a fallback for a refused exterior boundary fit.
use super::*;
use fs_couple::render::plate::impact::{ImpactSystem,ImpactSubstepConfig,radiation::Pole};

fn drum(frames:u64,internal:bool)->Experiment {
    let neck=crate::cavity::NeckOptions {exterior_load:internal,radius_m:0.003,
        effective_length_m:0.008,resistance_pa_s_m3:1000.,
        azimuth_rad:core::f64::consts::PI/8.,axial_position_m:0.08};
    let wires=crate::snare::SnareSet {strands:2,modes_per_strand:2,contact_cells:4,
        ..crate::snare::SnareSet::reference(false)};
    let specimen=crate::drum_spec::Spec {radial_intervals:2,azimuths:8,
        ..crate::drum_spec::Spec::reference()};
    crate::drum_with_radiation(frames*SUBSTEPS as u64,MECHANICAL_DT,true,false,
        Some(wires),true,crate::Stroke {speed_m_s:0.8,position_m:Some([0.06,0.01])},
        true,Some(neck),Some(specimen),
        Some(crate::Stroke {speed_m_s:0.6,position_m:Some([-0.05,0.02])}),
        &[],25.,None,true).unwrap()
}
fn prepared(mut e:Experiment)->Experiment {
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig {max_depth:8,max_attempts:511}).unwrap();e
}
fn owner(s:&Mechanics)->&ImpactSystem {
    match s {
        Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,Mechanics::Substepped(s)=>s,
        Mechanics::Driven {inner,..}=>owner(inner),Mechanics::Prepared(_)=>panic!("joint time owner required"),
    }
}
fn receivers()->[Receiver;2] {
    [Receiver::FinitePoint([0.08,0.05,0.35]),Receiver::FinitePoint([-0.12,0.05,0.4])]
}
fn band()->radiation_spec::Spec {
    // Same resolved 500 Hz observation window as the existing small-drum
    // prescribed-vent BEM fixture. Do not bypass the default-band TooCoarse gate.
    radiation_spec::Spec::parse(include_str!("coupled-vent.radiation")).unwrap()
}
fn declared_load(sources:usize)->Model {
    Model {ports:sources,poles:vec![
        Pole {omega:1200.,zeta:0.2,coupling:(0..sources).map(|i|
            if i+1==sources {180.}else if i%2==0 {-120.}else {70.}).collect()},
        Pole {omega:2400.,zeta:0.15,coupling:(0..sources).map(|i|
            if i+1==sources {-60.}else {40.}).collect()},
    ]}
}

#[test]
fn vent_feedback_rejects_effective_impedance_missing_mouth_and_interior_source_before_playing() {
    let effective=drum(4,false);let before=effective.system.state().to_vec();
    assert!(admit_instrument(&effective).is_err());
    assert_eq!(effective.system.state(),before);assert_eq!(owner(&effective.system).samples(),0);
    let mut e=drum(4,true);admit_instrument(&e).unwrap();
    let before=e.system.state().to_vec();
    let ids=e.acoustics.as_ref().unwrap().state_modes.clone();
    let mouth_row=e.acoustics.as_mut().unwrap().weights.pop().unwrap();
    let mouth=e.acoustics.as_mut().unwrap().state_modes.pop().unwrap();
    assert!(admit_instrument(&e).is_err());
    e.acoustics.as_mut().unwrap().state_modes.push(mouth);
    e.acoustics.as_mut().unwrap().weights.push(mouth_row);
    e.acoustics.as_mut().unwrap().state_modes.push(mouth);
    assert!(admit_instrument(&e).is_err());
    e.acoustics.as_mut().unwrap().state_modes.pop();
    let interior=e.air.as_ref().unwrap().coupling.structural_modes();
    e.acoustics.as_mut().unwrap().state_modes[0]=interior;
    assert!(admit_instrument(&e).is_err());
    e.acoustics.as_mut().unwrap().state_modes=ids.clone();
    admit_instrument(&e).unwrap();
    assert_eq!(e.system.state(),before);assert_eq!(owner(&e.system).samples(),0);
    let mut wrong_order=ids;wrong_order.reverse();
    assert!(attach_model(e,&declared_load(wrong_order.len()),&wrong_order).is_err());

    // This is the integrated preparation front door, not only a helper test.
    // 8000 Hz violates the declared 8 mm slug's chart before any BEM fitting.
    let wide=radiation_spec::Spec {band_hz:[40.,8000.],..band()};
    let error=prepare(drum(4,true),4,20.,&receivers(),wide,&CancelGate::new_clock_free())
        .err().expect("broad neck band must refuse");
    assert!(error.to_string().contains("compact-neck"));
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(prepare(drum(4,true),4,20.,&receivers(),band(),&cancelled).is_err());
}

#[test]
fn played_vent_reaction_changes_neck_flow_and_keeps_total_energy_and_exact_retry() {
    let original=drum(32,true);let prefix=original.system.state().to_vec();
    let sources=original.acoustics.as_ref().unwrap().state_modes.clone();
    let model=declared_load(sources.len());
    let mut loaded=prepared(attach_model(original,&model,&sources).unwrap());
    let mut clean=prepared(attach_model(drum(32,true),&model,&sources).unwrap());
    let mut unloaded=prepared(drum(32,true));
    assert_eq!(&loaded.system.state()[..prefix.len()],prefix);
    assert_eq!(loaded.system.state().len(),prefix.len()+2*model.poles.len());
    assert_eq!(loaded.force,unloaded.force);assert_eq!(loaded.observer_a,unloaded.observer_a);
    assert_eq!(loaded.observer_b,unloaded.observer_b);
    let port=loaded.air.as_ref().unwrap().coupling.neck_radiation_port(0).unwrap();
    assert_eq!(sources.last(),Some(&port.coordinate));
    assert!(!sources.contains(&0) && !sources.contains(&loaded.second_stick.unwrap().coordinate));
    let initial=owner(&loaded.system).stored_energy_j();let gate=CancelGate::new_clock_free();
    let (mut loss,mut supplied,mut changed,mut pressure_peak)=(0.,0.,0.0_f64,0.0_f64);
    for tick in 1..=512 {
        if tick==160 {
            let state=loaded.system.state().to_vec();
            let cancelled=CancelGate::new_clock_free();cancelled.request();
            assert!(loaded.system.step(&loaded.force,&cancelled).is_err());
            let mut bad=loaded.force.clone();bad[0]=f64::NAN;
            assert!(loaded.system.step(&bad,&gate).is_err());
            assert_eq!(loaded.system.state(),state);
        }
        let f=loaded.system.step(&loaded.force,&gate).unwrap();
        clean.system.step(&clean.force,&gate).unwrap();unloaded.system.step(&unloaded.force,&gate).unwrap();
        assert_eq!(loaded.system.state(),clean.system.state(),"same-path refusal/retry");
        assert_eq!(f.time_s,tick as f64*MECHANICAL_DT);
        loss+=f.dissipated_energy_j;supplied+=f.supplied_work_j;
        assert!(f.balance_residual_j.abs()<1e-7);
        assert!((f.stored_energy_j+loss-initial-supplied).abs()<1e-6);
        let system=owner(&loaded.system);let air=&loaded.air.as_ref().unwrap().coupling;
        let neck=air.neck_observation(system.state(),0).unwrap();
        let pressure=air.neck_exterior_pressure(system,0).unwrap();
        let reaction=system.radiation_force(port.coordinate).unwrap();
        let power=pressure*neck.volume_flow_m3_s;
        assert!((power+reaction*system.state()[2*port.coordinate+1]).abs()<1e-12*(1.+power.abs()));
        pressure_peak=pressure_peak.max(pressure.abs());
        changed=changed.max((system.state()[2*port.coordinate+1]
            -unloaded.system.state()[2*port.coordinate+1]).abs());
    }
    assert!(changed>1e-12 && pressure_peak>0. && loss>0.);
    assert!(observation(&loaded.system).unwrap().stored_energy_j>0.);
    assert!(loaded.system.membrane_observation(1).unwrap().stretching_energy_j>0.);
    assert_ne!(&loaded.system.state()[..prefix.len()],unloaded.system.state());
}

#[test]
fn mouth_bem_pressure_transposes_into_the_same_head_and_neck_work_coordinates() {
    let e=drum(4,true);let boundary=e.acoustics.as_ref().unwrap();
    let port=e.air.as_ref().unwrap().coupling.neck_radiation_port(0).unwrap();
    let surface=SpherePanels::from_triangles(boundary.triangles.clone()).unwrap();
    let mouth=boundary.state_modes.iter().position(|&i|i==port.coordinate).unwrap();
    let row=&boundary.weights[mouth];
    let volume_weight=row.iter().zip(surface.areas()).map(|(b,a)|b*a).sum::<f64>();
    assert!((volume_weight-port.volume_weight_m2_per_sqrt_kg).abs()<1e-12*volume_weight);
    let omega=core::f64::consts::TAU*250.;let medium=Medium::air();
    let fields=acceleration_fields(&boundary.weights,omega);
    let refs:Vec<_>=fields.iter().map(Vec::as_slice).collect();
    let geometry=reciprocal_geometry(&surface,omega/medium.sound_speed,band(),refs.len(),1,
        &CancelGate::new_clock_free()).unwrap();
    let solutions=geometry.solve_batch(omega/medium.sound_speed,medium,&refs,||false).unwrap();
    let traces:Vec<&[C64]>=solutions.iter().map(|s|s.pressure.as_slice()).collect();
    let z=project_pressures(&surface,&boundary.weights,&traces,omega).unwrap();
    let n=boundary.weights.len();let mut cross_reaction=false;
    let scale=z.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    for i in 0..n {for j in 0..n {assert!((z[i*n+j]-z[j*n+i]).abs()<1e-9*scale);}}
    for (source,solution) in solutions.iter().enumerate() {
        assert!(solution.flux_balance<1e-8);
        let mouth_pressure=solution.pressure.iter().zip(row).zip(surface.areas())
            .fold(C64::ZERO,|sum,((p,b),a)|sum+p.scale(b*a/volume_weight))
            *C64::new(0.,-omega);
        let expected=mouth_pressure.scale(volume_weight);
        assert!((z[mouth*n+source]-expected).abs()<1e-11*(1.+expected.abs()));
        if source!=mouth {cross_reaction|=expected.abs()>1e-12;}
    }
    assert!(cross_reaction,"head-driven exterior pressure must also react on the mouth");
}

#[test]
fn vent_feedback_bakes_close_directional_stereo_from_one_loaded_mechanical_trajectory() {
    let frames=96;let original=drum(frames as u64,true);
    let prefix=original.system.state().to_vec();let gate=CancelGate::new_clock_free();
    let microphones=crate::acoustics::receivers::input::Spec::parse(
        include_str!("coupled-vent-microphones.frm")).unwrap();
    microphones.admit_command("snare-mic",false,false).unwrap();
    let microphones=microphones.into_receivers();
    assert!(matches!(microphones[0],Receiver::NearField(m) if m.pressure_fraction()==0.5));
    assert!(matches!(microphones[1],Receiver::NearField(m) if m.pressure_fraction()==1.));
    let (loaded,observer)=prepare(original,frames,20.,&microphones,band(),&gate).unwrap();
    assert_eq!(&loaded.system.state()[..prefix.len()],prefix);
    assert_eq!(observer.baked.len(),2);
    assert!(observer.baked.iter().all(|b|b.propagation_delay_s==0. && b.pressure_gain==1.));
    assert_eq!(observer.sources.last(),Some(&loaded.air.as_ref().unwrap().coupling
        .neck_radiation_port(0).unwrap().coordinate));
    let initial=owner(&loaded.system).stored_energy_j();let mut loaded=prepared(loaded);
    let mut unloaded=prepared(drum(frames as u64,true));
    // Identical admitted observer bank; only the mechanics lacks the feedback.
    // This is an explicitly unloaded comparison, not a supported coupled export.
    let one_way=render_baked_with_gate(&mut unloaded,frames,20.,&observer.baked,&gate).unwrap();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    let before=loaded.system.state().to_vec();
    assert!(render_baked_with_gate(&mut loaded,frames,20.,&observer.baked,&cancelled).is_err());
    assert_eq!(loaded.system.state(),before);
    let wav=observer.render(&mut loaded,frames,20.,&gate).unwrap();
    assert_eq!(&wav[..4],b"RIFF");assert_eq!(&wav[8..12],b"WAVE");
    assert_eq!(u16::from_le_bytes(wav[22..24].try_into().unwrap()),2);
    assert_eq!(wav.len(),44+frames*4);assert_eq!(one_way.len(),wav.len());
    assert!(wav[44..].iter().any(|b|*b!=0));
    assert_ne!(&loaded.system.state()[..prefix.len()],unloaded.system.state());
    assert!(observation(&loaded.system).unwrap().stored_energy_j>0.);
    assert!(owner(&loaded.system).stored_energy_j()<initial);
    // The complete mechanical budget is consumed once, not once per receiver.
    let before=loaded.system.state().to_vec();
    assert!(loaded.system.step(&loaded.force,&gate).is_err());
    assert_eq!(loaded.system.state(),before);
}
