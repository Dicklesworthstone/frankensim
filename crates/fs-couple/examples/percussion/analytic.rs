//! Numerical realization only. The existing constructors own all physics.
use super::{Mechanics,Error};
pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let count=args.iter().filter(|s|s.as_str()=="--analytic-newton").count();
    if count>1 {return Err("--analytic-newton may be supplied only once".into());}
    args.retain(|s|s!="--analytic-newton");Ok(count==1)
}
impl Mechanics {
    pub fn into_analytic_nonlinear(self)->Result<Self,Error> {
        match self {
            Self::Reference(s)=>Ok(Self::Nonlinear(s.prepare_analytic()?)),
            Self::Nonlinear(mut s)=>{s.set_analytic_newton(true);Ok(Self::Nonlinear(s))},
            Self::Substepped(mut s)=>{s.set_analytic_newton(true);Ok(Self::Substepped(s))},
            Self::Prepared(_)=>Err("--analytic-newton needs nonlinear splash/drum/drum-stretch mechanics; no conversion of modal/snare physics".into()),
            Self::Driven {inner,drive}=>Ok(Self::Driven {inner:Box::new((*inner).into_analytic_nonlinear()?),drive}),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::{drum,splash_with_stroke,drum_with_mufflers,Stroke,muffling,acoustics};
    use fs_exec::CancelGate;
    #[test]
    fn numerical_option_preserves_physical_arguments_and_rejects_duplicates_or_modal_conversion() {
        let mut args=vec!["drum-stretch".into(),"--analytic-newton".into(),"--strike-speed-m-s".into(),"3".into()];
        assert!(option(&mut args).unwrap());assert_eq!(args,["drum-stretch","--strike-speed-m-s","3"]);
        assert!(!option(&mut args).unwrap());
        assert!(option(&mut vec!["--analytic-newton".into();2]).is_err());
        assert!(drum(4,2e-6,false,true).unwrap().system.into_analytic_nonlinear().is_err());
    }
    #[test]
    fn actual_hard_splash_substeps_keep_nonlinear_felt_memory_and_pressure_surface() {
        let mut e=splash_with_stroke(128,acoustics::MECHANICAL_DT,true,Stroke {speed_m_s:4.0,..Stroke::default()}).unwrap();
        let before=e.system.state().to_vec();assert!(before.len()>2*e.force.len());
        assert!(e.acoustics.is_some());e.system=e.system.into_analytic_nonlinear().unwrap();
        // Native fixed-step execution stalls at this unchanged 4 m/s impact.
        // Exercise the explicitly selected recovery path, not a gentler strike
        // or relaxed energy tolerance. The pressure geometry is untouched.
        e.system=e.system.with_impact_substeps(fs_couple::render::plate::impact::ImpactSubstepConfig {
            max_depth:8,max_attempts:511 }).unwrap();
        assert_eq!(e.system.state(),before);
        let gate=CancelGate::new_clock_free();let mut refined=false;
        for tick in 1..=128 {
            let f=e.system.step(&e.force,&gate).unwrap();
            let Mechanics::Substepped(s)=&e.system else {panic!("substepped shell")};
            refined|=s.last_substeps().accepted_substeps>1;
            assert_eq!(f.time_s,tick as f64*acoustics::MECHANICAL_DT);
            assert!(f.stored_energy_j.is_finite() && f.balance_residual_j.abs()<1e-7);
        }
        let Mechanics::Substepped(s)=e.system else {panic!("nonlinear shell")};
        assert!(refined,"hard impact must exercise recovery rather than just the unchanged one-leaf path");
        assert!(s.felt_history(0).is_some());assert_eq!(s.samples(),128);
        assert!(s.state()[1]<before[1],"the physical stick has entered contact and decelerated");
    }
    #[test]
    fn two_stick_stretching_drums_cavity_and_spatial_loss_use_the_analytic_joint_solve() {
        let first=Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])};
        let second=Some(Stroke {speed_m_s:2.5,position_m:Some([-0.05,0.02])});
        let muffler=muffling::Muffler {surface:muffling::Surface::Batter,position_m:[0.08,0.01],resistance_n_s_m:0.4};
        let mut e=drum_with_mufflers(128,2e-6,false,false,None,true,first,true,None,None,second,&[muffler]).unwrap();
        assert!(e.air.is_some() && e.second_stick.is_some());
        let before=e.system.state().to_vec();e.system=e.system.into_analytic_nonlinear().unwrap();
        assert_eq!(e.system.state(),before);let gate=CancelGate::new_clock_free();
        for _ in 0..128 {
            let f=e.system.step(&e.force,&gate).unwrap();
            assert!(f.dissipated_energy_j>=0.0 && f.balance_residual_j.abs()<1e-7);
        }
        assert!(e.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
    }
}
