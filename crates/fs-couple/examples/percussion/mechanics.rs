//! Choose a numerical image only AFTER the shared geometric construction.
//! No mesh, material, force, radiation or output code is duplicated here.
use super::{Error, config};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::render::plate::impact::{ImpactBody, ImpactError, ImpactSystem, PreparedImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::linear::{LinearImpactSystem, LinearImpactConfig, VolumeConnection};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::{ModalContactConfig, multiple::MultiContactConfig}};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

pub enum Mechanics {
    Reference(ImpactSystem),
    Prepared(LinearImpactSystem),
    Nonlinear(PreparedImpactSystem),
}
/// Only the diagnostics shared by both images. In particular, a normal-force
/// residual is not exposed under the reference solver's different residual unit.
pub struct Frame {
    pub time_s: f64,
    pub stored_energy_j: f64,
    pub felt_crush_loss_j: f64,
    pub dissipated_energy_j: f64,
    pub balance_residual_j: f64,
}
fn prepared_config(steps: u64, dt_s: f64) -> Result<LinearImpactConfig, Error> {
    let rate = dt_s.recip().round();
    if !dt_s.is_finite() || dt_s <= 0.0 || !rate.is_finite()
        || rate < 1.0 || rate > f64::from(u32::MAX)
        || (1.0 / rate).to_bits() != dt_s.to_bits()
    { return Err("prepared image requires an exact integer-Hz clock; no retiming is permitted".into()); }
    let original = config(steps, dt_s);
    Ok(LinearImpactConfig {
        sample_rate_hz: rate as u32, max_steps: steps,
        maximum_generalized_force: original.maximum_generalized_force,
        component: ModalAcousticTimeBudget { maximum_total_energy_j: original.maximum_energy_j,
            ..ModalAcousticTimeBudget::audible_reference() },
        coupling: ModalCouplingConfig {
            max_modes: 64, max_connections: 8, max_setup_terms: 100000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: original.maximum_energy_j,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: original.energy_absolute_tolerance_j,
            energy_relative_tolerance: original.energy_relative_tolerance,
        },
        contact: ModalContactConfig {
            max_iterations: 100, maximum_force_n: 1e4, maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1e-9, force_relative_tolerance: 1e-9,
        },
        multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 100000 },
    })
}
/// Extract only the numerical option; the existing playing parser owns physics.
pub fn prepared_option(args: &mut Vec<String>) -> Result<bool, Error> {
    let count = args.iter().filter(|arg| arg.as_str() == "--prepared-nonlinear").count();
    if count > 1 { return Err("--prepared-nonlinear may be supplied only once".into()); }
    args.retain(|arg| arg != "--prepared-nonlinear");
    Ok(count == 1)
}
impl Mechanics {
    /// Prepare the exact nonlinear model after construction, with no reset,
    /// retiming, modal truncation, or replacement by the linear-contact image.
    pub fn into_prepared_nonlinear(self) -> Result<Self, Error> {
        match self {
            Self::Reference(system) => Ok(Self::Nonlinear(system.prepare()?)),
            Self::Nonlinear(system) => Ok(Self::Nonlinear(system)),
            Self::Prepared(_) => Err("--prepared-nonlinear needs splash, drum or drum-stretch; the modal/snare image is a different physical admission".into()),
        }
    }

    pub fn prepared(bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volume: VolumeSpring,
        reference_area_m2: f64, steps: u64, dt_s: f64) -> Result<Self, Error>
    {
        Self::with_configuration(bodies,contacts,volume,reference_area_m2,prepared_config(steps,dt_s)?)
    }
    /// Separate explicit work envelope for the multi-strand example. Existing
    /// drum and cymbal commands keep their original limits and numerical image.
    pub fn prepared_snares(bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volume: VolumeSpring,
        reference_area_m2: f64, steps: u64, dt_s: f64) -> Result<Self, Error>
    {
        let mut configuration=prepared_config(steps,dt_s)?;
        configuration.coupling.max_modes=256;
        configuration.multiple.max_contacts=512;
        configuration.multiple.max_setup_terms=50_000_000;
        Self::with_configuration(bodies,contacts,volume,reference_area_m2,configuration)
    }
    fn with_configuration(bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volume: VolumeSpring,
        reference_area_m2: f64, configuration: LinearImpactConfig) -> Result<Self, Error>
    {
        let system = LinearImpactSystem::new(bodies, contacts,
            vec![VolumeConnection { spring: volume, reference_area_m2 }], configuration,
            &CancelGate::new_clock_free())?;
        eprintln!("mechanical image: prepared exact-ZOH bodies plus simultaneous volume/contact reactions; unchanged physical cards; native real-time performance is unqualified");
        Ok(Self::Prepared(system))
    }
    pub fn membrane_observation(&self,body:usize) -> Option<fs_couple::render::plate::impact::membrane::MembraneObservation> {
        match self { Self::Reference(s)=>s.membrane_observation(body), Self::Prepared(_)=>None,
            Self::Nonlinear(s)=>s.membrane_observation(body) }
    }
    pub fn state(&self) -> &[f64] {
        match self { Self::Reference(s) => s.state(), Self::Prepared(s) => s.state(), Self::Nonlinear(s) => s.state() }
    }
    pub fn step(&mut self, external: &[f64], gate: &CancelGate) -> Result<Frame, ImpactError> {
        Ok(match self {
            Self::Reference(s) => {
                let f = s.step(external, gate)?;
                Frame { time_s: f.time_s, stored_energy_j: f.stored_energy_j,
                    felt_crush_loss_j: f.felt_crush_loss_j, dissipated_energy_j: f.dissipated_energy_j,
                    balance_residual_j: f.balance_residual_j }
            }
            Self::Nonlinear(s) => {
                let f = s.step(external, gate)?;
                Frame { time_s: f.time_s, stored_energy_j: f.stored_energy_j,
                    felt_crush_loss_j: f.felt_crush_loss_j, dissipated_energy_j: f.dissipated_energy_j,
                    balance_residual_j: f.balance_residual_j }
            }
            Self::Prepared(s) => {
                let f = s.step(external, gate)?;
                Frame { time_s: f.time_s, stored_energy_j: f.stored_energy_j,
                    felt_crush_loss_j: 0.0, dissipated_energy_j: f.dissipated_energy_j,
                    balance_residual_j: f.balance_residual_j }
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn prepared_images_keep_both_existing_clocks_without_retiming() {
        for dt in [2e-6, super::super::acoustics::MECHANICAL_DT] {
            let c = prepared_config(123, dt).unwrap();
            assert_eq!((1.0/f64::from(c.sample_rate_hz)).to_bits(), dt.to_bits());
            assert_eq!(c.max_steps, 123);
            assert_eq!(c.coupling.maximum_total_energy_j, config(123,dt).maximum_energy_j);
        }
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0/48000.25] {
            assert!(prepared_config(10,dt).is_err());
        }
    }
    #[test]
    fn prepared_drum_consumes_the_same_real_geometry_and_initial_physical_state() {
        // Actual existing film assembly/eigensolve, not an authored frequency bank.
        // Native execution remains required to establish this regression passes.
        let mut old = super::super::drum(128, 2e-6, false, false).unwrap();
        let mut new = super::super::drum(128, 2e-6, false, true).unwrap();
        assert_eq!(old.system.state(), new.system.state());
        assert_eq!(old.observer_a, new.observer_a);
        assert_eq!(old.observer_b, new.observer_b);
        assert_eq!(old.pressure.as_ref().unwrap().areas, new.pressure.as_ref().unwrap().areas);
        assert!(matches!(&new.system, Mechanics::Prepared(_)));
        let gate = CancelGate::new_clock_free();
        let mut delta = 0.0_f64;
        for _ in 0..128 {
            old.system.step(&old.force, &gate).unwrap();
            new.system.step(&new.force, &gate).unwrap();
            for (i, (&a,&b)) in old.system.state().iter().zip(new.system.state()).enumerate() {
                delta = delta.max((a-b).abs() * if i%2==0 {4000.0} else {1.0});
            }
        }
        assert!(delta < 1e-3, "prepared/reference onset discrepancy {delta}");
    }
}

#[cfg(test)]
mod nonlinear_tests {
    use super::*;
    #[test]
    fn numerical_flag_does_not_change_physical_arguments() {
        let mut args=vec!["drum-stretch-mic".into(),"100".into(),"--prepared-nonlinear".into(),
            "--strike-speed-m-s".into(),"4.0".into()];
        assert!(prepared_option(&mut args).unwrap());
        let (positional,stroke)=super::super::playing::parse(args).unwrap();
        assert_eq!(positional,vec!["drum-stretch-mic","100"]); assert_eq!(stroke.speed_m_s,4.0);
        assert!(prepared_option(&mut vec!["--prepared-nonlinear".into();2]).is_err());
    }
    #[test]
    fn actual_stretching_drum_prepares_without_changing_geometry_or_initial_motion() {
        let mut experiment=super::super::drum_with_playing(64,2e-6,false,false,None,true,
            super::super::Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])}).unwrap();
        let state=experiment.system.state().to_vec(); let force=experiment.force.clone();
        let observation=experiment.observer_a.clone();
        experiment.system=experiment.system.into_prepared_nonlinear().unwrap();
        assert_eq!(experiment.system.state(),state);
        assert_eq!(experiment.force,force); assert_eq!(experiment.observer_a,observation);
        assert!(matches!(&experiment.system,Mechanics::Nonlinear(_)));
        let gate=CancelGate::new_clock_free(); let mut stretch=0.0_f64;
        for _ in 0..64 {
            let f=experiment.system.step(&experiment.force,&gate).unwrap();
            assert!(f.balance_residual_j.abs()<1e-7);
            let head=experiment.system.membrane_observation(1).unwrap();
            stretch=stretch.max(head.stretching_energy_j); assert!(head.maximum_slope<=0.2);
        }
        assert!(stretch>0.0,"the prepared executable must retain actual nonlinear head storage");
    }
    #[test]
    fn actual_splash_preserves_shell_and_each_felt_history_when_prepared() {
        let mut experiment=super::super::splash(192,2e-6,false).unwrap();
        let state=experiment.system.state().to_vec();
        let Mechanics::Reference(reference)=&experiment.system else {panic!("reference construction");};
        let histories:Vec<_>=(0..6).map(|i|reference.felt_history(i)).collect();
        let energy=reference.stored_energy_j();
        experiment.system=experiment.system.into_prepared_nonlinear().unwrap();
        let Mechanics::Nonlinear(prepared)=&experiment.system else {panic!("prepared construction");};
        assert_eq!(prepared.state(),state); assert_eq!(prepared.stored_energy_j().to_bits(),energy.to_bits());
        for (i,h) in histories.iter().enumerate() {assert_eq!(&prepared.felt_history(i),h);}
        let gate=CancelGate::new_clock_free();
        for _ in 0..192 {experiment.system.step(&experiment.force,&gate).unwrap();}
    }
}
