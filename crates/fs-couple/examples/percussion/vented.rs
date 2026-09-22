//! Explicit observer-only vent radiation, never an implicit physical upgrade.
use super::{Error,cavity::{self,NeckOptions}};

pub fn admit(neck:Option<NeckOptions>,distributed:bool,command:&str,prescribed:bool)->Result<(),Error> {
    if !prescribed {return cavity::admit_neck_command(neck,distributed,command);}
    if neck.is_none() || !distributed || !matches!(command,
        "drum-wav"|"drum-mic"|"drum-modal-wav"|"drum-modal-mic"|
        "drum-stretch-wav"|"drum-stretch-mic"|"snare-wav"|"snare-mic"|"snare-off-wav"|"snare-off-mic") {
        return Err("--prescribed-vent-radiation requires a drum/snare audio command, --cavity-modes and an explicit --cavity-neck; it does not add radiation feedback".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{acoustics,drum_with_radiation,drum_spec,Experiment,Stroke,snare::SnareSet};
    use fs_exec::CancelGate;
    use fs_couple::render::plate::impact::ImpactSubstepConfig;

    fn neck()->NeckOptions {NeckOptions {radius_m:0.003,effective_length_m:0.008,
        resistance_pa_s_m3:1000.0,azimuth_rad:std::f64::consts::PI/8.0,axial_position_m:0.08}}
    fn build(audio:bool,nonlinear:bool,frames:u64)->Experiment {
        let wires=SnareSet {strands:2,modes_per_strand:2,contact_cells:4,..SnareSet::reference(false)};
        let spec=drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()};
        let mut e=drum_with_radiation(frames*acoustics::SUBSTEPS as u64,acoustics::MECHANICAL_DT,
            audio,!nonlinear,Some(wires),nonlinear,
            Stroke {speed_m_s:2.0,position_m:Some([0.06,0.01])},true,Some(neck()),Some(spec),
            Some(Stroke {speed_m_s:1.6,position_m:Some([-0.05,0.02])}),&[],25.0,None,audio).unwrap();
        if nonlinear {
            e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
                ImpactSubstepConfig {max_depth:8,max_attempts:511}).unwrap();
        }
        e
    }

    #[test]
    fn selecting_the_approximation_never_silently_promotes_a_bare_vent_request() {
        for command in ["drum-wav","drum-mic","drum-modal-wav","drum-modal-mic",
            "drum-stretch-wav","drum-stretch-mic","snare-wav","snare-mic","snare-off-wav","snare-off-mic"] {
            admit(Some(neck()),true,command,true).unwrap();
            assert!(admit(Some(neck()),true,command,false).is_err());
            assert!(admit(None,true,command,true).is_err());
            assert!(admit(Some(neck()),false,command,true).is_err());
        }
        for command in ["drum","drum-modal","drum-stretch","snare","snare-off"] {
            admit(Some(neck()),true,command,false).unwrap();
            assert!(admit(Some(neck()),true,command,true).is_err());
        }
        for command in ["splash","splash-mic","splash-wav","unknown"] {
            assert!(admit(Some(neck()),true,command,true).is_err());
        }
    }

    #[test]
    fn vent_observation_retains_actual_neck_and_two_stick_nonlinear_snare_state_work_and_clocks() {
        let gate=CancelGate::new_clock_free();
        for nonlinear in [false,true] {
            let mut observed=build(true,nonlinear,16);let mut mechanics=build(false,nonlinear,16);
            let port=observed.air.as_ref().unwrap().coupling.neck_radiation_port(0).unwrap();
            let modes=observed.acoustics.as_ref().unwrap().state_modes();
            assert_eq!(modes.last(),Some(&port.coordinate));
            let second=observed.second_stick.unwrap();
            assert!(!modes.contains(&second.coordinate));assert!(!modes.contains(&0));
            assert_eq!(modes.len(),second.coordinate); // heads plus one neck, neither stick
            assert_eq!(observed.force,mechanics.force);assert_eq!(observed.system.state(),mechanics.system.state());
            let mut peak_flow=0.0_f64;let mut peak_work=0.0_f64;
            for _ in 0..256 {
                let f=observed.system.step(&observed.force,&gate).unwrap();
                let g=mechanics.system.step(&mechanics.force,&gate).unwrap();
                assert_eq!(observed.system.state(),mechanics.system.state());
                assert_eq!(f.time_s,g.time_s);assert_eq!(f.stored_energy_j,g.stored_energy_j);
                assert_eq!(f.dissipated_energy_j,g.dissipated_energy_j);assert_eq!(f.supplied_work_j,g.supplied_work_j);
                let neck=observed.air.as_ref().unwrap().coupling.neck_observation(observed.system.state(),0).unwrap();
                peak_flow=peak_flow.max(neck.volume_flow_m3_s.abs());peak_work=peak_work.max(f.stored_energy_j);
            }
            assert!(peak_flow>0.0 && peak_work>0.0);
            assert_eq!(observed.system.membrane_observation(1).is_some(),nonlinear);
        }
    }

    #[test]
    fn prescribed_vent_runs_through_real_bem_fits_stereo_and_pcm_without_advancing_mechanics_twice() {
        let frames=192;let mut e=build(true,false,frames);let mut reference=build(false,false,frames);
        let wav=acoustics::stereo::render_receivers(&mut e,frames as usize,20.0,&[
            acoustics::Receiver::FinitePoint([0.08,0.05,0.35]),
            acoustics::Receiver::FinitePoint([-0.12,0.05,0.4])]).unwrap();
        assert_eq!(&wav[..4],b"RIFF");assert_eq!(&wav[8..12],b"WAVE");
        assert_eq!(u16::from_le_bytes(wav[22..24].try_into().unwrap()),2);
        assert_eq!(u32::from_le_bytes(wav[24..28].try_into().unwrap()),48000);
        assert_eq!(wav.len(),44+frames as usize*4);
        assert!(wav[44..].iter().any(|b|*b!=0),"actual aperture/head pressure must reach PCM");
        let gate=CancelGate::new_clock_free();
        for _ in 0..frames*acoustics::SUBSTEPS as u64 {reference.system.step(&reference.force,&gate).unwrap();}
        assert_eq!(e.system.state(),reference.system.state());
    }
}
