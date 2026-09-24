//! Passive multiport acoustic storage in the loaded board's velocity basis.
//!
//! Each row l realizes Z(s) = l^T l s/(s^2 + 2*zeta*omega*s + omega^2).
//! Acoustic coordinates have unit modal mass. Their forcing is l*v_board;
//! the opposing board force is -l^T*v_air. Thus the interconnection does
//! exactly zero continuous work. No entrywise SISO passivity assumption.
//!
//! Free acoustic motion is lowered from the SAME exact oscillator owner as
//! the piano bank. The power interconnection is the collective exact flow of
//! the complete rectangular coupling, prepared by the shared fs-phs owner. The
//! coupling/free-air/piano/free-air/coupling composition is second-order, not
//! an exact full-system propagator. Each coupling flow conserves
//! combined kinetic energy to its checked decomposition/roundoff tolerance;
//! accuracy still requires rate convergence. No allocation in these flows.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_math::c64::C64;
use fs_phs::{PortExchangeBudget, PreparedPortExchange};

pub const MAX_AIR_STATES: usize = 1024;

// Shared positive-real model; collective exchange does not alter its impedance.
pub use fs_couple::render::plate::impact::radiation::{Model,Pole};

#[derive(Clone, Copy, Default)]
struct Free { xx:f64, xp:f64, px:f64, pp:f64 }

pub struct Prepared {
    ports: usize,
    // x = omega*q_air, p = v_air, so storage is their Euclidean norm / 2.
    state: Vec<[f64;2]>,
    saved: Vec<[f64;2]>,
    free: Vec<Free>,
    exchange: PreparedPortExchange,
    velocities: Vec<f64>,
}
impl Prepared {
    pub fn new(model:&Model, rate:u32, ports:usize) -> Result<Self,String> {
        model.validate()?;
        if rate < 8_000 || ports != model.ports
            || model.poles.iter().any(|p| p.omega >= 0.9*std::f64::consts::PI*f64::from(rate)) {
            return Err("acoustic model does not cover the loaded basis or mechanical rate".into());
        }
        let half_rate=rate.checked_mul(2).ok_or("acoustic half-step rate overflow")?;
        let modes:Vec<_>=model.poles.iter().map(|p| ModalAcousticMode {
            angular_frequency_rad_s:p.omega, damping_ratio:p.zeta,
            pressure_per_modal_velocity:C64::ZERO,
        }).collect();
        let mut oscillator=ModalAcousticTimeModel::try_new(half_rate,modes,
            ModalAcousticTimeBudget::audible_reference()).map_err(|e|e.to_string())?;
        let zero=vec![0.;model.poles.len()];
        let mut states:Vec<_>=model.poles.iter().map(|p|ModalAcousticState {
            displacement_m_sqrt_kg:1./p.omega,velocity_m_sqrt_kg_per_s:0.,
        }).collect();
        oscillator.restore_states(&states).map_err(|e|e.to_string())?;
        oscillator.step(&zero).map_err(|e|e.to_string())?;
        let mut free=vec![Free::default();states.len()];
        for ((t,s),p) in free.iter_mut().zip(oscillator.states()).zip(&model.poles) {
            t.xx=p.omega*s.displacement_m_sqrt_kg;t.px=s.velocity_m_sqrt_kg_per_s;
        }
        states.fill(ModalAcousticState {displacement_m_sqrt_kg:0.,velocity_m_sqrt_kg_per_s:1.});
        oscillator.restore_states(&states).map_err(|e|e.to_string())?;
        oscillator.step(&zero).map_err(|e|e.to_string())?;
        for ((t,s),p) in free.iter_mut().zip(oscillator.states()).zip(&model.poles) {
            t.xp=p.omega*s.displacement_m_sqrt_kg;t.pp=s.velocity_m_sqrt_kg_per_s;
        }
        let coupling:Vec<_>=model.poles.iter().flat_map(|p|p.coupling.iter().copied()).collect();
        // One collective half-step, not separate physical-port rotations.
        // The old full-step infinity-norm limit .25 remains .125 per half-step.
        let exchange=PreparedPortExchange::new(ports,model.poles.len(),&coupling,
            0.5/f64::from(rate),PortExchangeBudget {max_left:32,max_right:MAX_AIR_STATES,
                max_setup_terms:64*1024*1024,maximum_dt_coupling:0.125})
            .map_err(|e|e.to_string())?;
        if free.iter().any(|t|[t.xx,t.xp,t.px,t.pp].iter().any(|x|!x.is_finite())) {
            return Err("nonfinite prepared acoustic transition".into());
        }
        Ok(Self {ports,state:vec![[0.;2];free.len()],saved:vec![[0.;2];free.len()],
            velocities:vec![0.;free.len()],free,exchange})
    }
    pub fn energy(&self)->f64 {0.5*self.state.iter().map(|s|s[0]*s[0]+s[1]*s[1]).sum::<f64>()}
    pub fn checkpoint(&mut self) {self.saved.copy_from_slice(&self.state);}
    pub fn restore(&mut self) {self.state.copy_from_slice(&self.saved);}
    fn exchange(&mut self,board:&mut[f64])->Result<(),&'static str> {
        if board.len()!=self.ports {return Err("radiation exchange has a different board basis");}
        for (v,s) in self.velocities.iter_mut().zip(&self.state) {*v=s[1];}
        self.exchange.apply(board,&mut self.velocities)
            .map_err(|_|"collective radiation exchange failed finite/energy admission")?;
        for (s,&v) in self.state.iter_mut().zip(&self.velocities) {s[1]=v;}
        Ok(())
    }
    fn free_half(&mut self)->f64 {
        let before=self.energy();
        for (s,t) in self.state.iter_mut().zip(&self.free) {
            let [x,p]=*s;*s=[t.xx*x+t.xp*p,t.px*x+t.pp*p];
        }
        // Signed floating-point loss, not a clipped energy-balance correction.
        before-self.energy()
    }
    pub fn before(&mut self,board:&mut[f64])->Result<f64,&'static str> {
        self.exchange(board)?;Ok(self.free_half())
    }
    pub fn after(&mut self,board:&mut[f64])->Result<f64,&'static str> {
        let loss=self.free_half();self.exchange(board)?;Ok(loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn model(zeta:f64)->Model {Model {ports:2,poles:vec![
        Pole {omega:1200.,zeta,coupling:vec![180.,-90.]},
        Pole {omega:3200.,zeta,coupling:vec![50.,220.]},
    ]}}
    #[test]
    fn complete_complex_load_is_reciprocal_passive_and_mass_like_at_low_frequency() {
        let m=model(0.2);
        for w in [10.,500.,1200.,3200.,10000.] {
            let z=m.impedance(w).unwrap();assert_eq!(z[1],z[2]);
            for u in [[1.,0.],[0.,1.],[1.,1.],[1.,-2.]] {
                let power=(0..2).map(|i|(0..2).map(|j|u[i]*z[2*i+j].re*u[j]).sum::<f64>()).sum::<f64>();
                assert!(power>=-1e-14);
            }
        }
        assert!(m.impedance(10.).unwrap()[0].im<0.);
        assert!(m.impedance(500.).unwrap()[1].abs()>1e-6);
    }
    #[test]
    fn free_body_and_acoustic_storage_close_work_without_artificial_loss() {
        for damping in [0.,0.2] {
            let mut air=Prepared::new(&model(damping),192_000,2).unwrap();
            let mut v=[0.03,-0.02];let initial=0.5*(v[0]*v[0]+v[1]*v[1]);let mut loss=0.;
            for _ in 0..10000 {loss+=air.before(&mut v).unwrap();loss+=air.after(&mut v).unwrap();}
            let final_energy=air.energy()+0.5*(v[0]*v[0]+v[1]*v[1]);
            assert!((final_energy+loss-initial).abs()<2e-13);
            assert!(air.energy()>1e-10);assert!(loss>=-1e-12);
            if damping>0. {assert!(loss>1e-8);}
        }
    }
    #[test]
    fn refused_frame_can_restore_air_history_exactly() {
        let mut air=Prepared::new(&model(0.2),192_000,2).unwrap();let mut v=[0.1,0.2];
        air.before(&mut v).unwrap();air.after(&mut v).unwrap();air.checkpoint();let old=air.state.clone();
        let saved=v;let loss=air.before(&mut v).unwrap()+air.after(&mut v).unwrap();let candidate=air.state.clone();
        air.restore();assert_eq!(air.state,old);v=saved;
        let retried=air.before(&mut v).unwrap()+air.after(&mut v).unwrap();
        assert_eq!(candidate,air.state);assert_eq!(loss.to_bits(),retried.to_bits());
    }
    #[test]
    fn bad_or_under_resolved_models_refuse_without_changing_a_valid_model() {
        let good=model(0.2);assert!(Prepared::new(&good,192_000,3).is_err());
        for bad in [f64::NAN,-1.,1.] {let mut m=good.clone();m.poles[0].zeta=bad;assert!(m.validate().is_err());}
        let mut m=good.clone();m.poles[0].coupling[0]=1e9;
        assert!(Prepared::new(&m,192_000,2).is_err());
        let mut m=good;m.poles[0].coupling.pop();assert!(m.validate().is_err());
    }
}
