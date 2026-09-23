//! Passive multiport acoustic storage in the loaded board's velocity basis.
//!
//! Each row l realizes Z(s) = l^T l s/(s^2 + 2*zeta*omega*s + omega^2).
//! Acoustic coordinates have unit modal mass. Their forcing is l*v_board;
//! the opposing board force is -l^T*v_air. Thus the interconnection does
//! exactly zero continuous work. No entrywise SISO passivity assumption.
//!
//! Free acoustic motion is lowered from the SAME exact oscillator owner as
//! the piano bank. The power interconnection is a symmetric composition of
//! exact two-coordinate rotations. Coupling/free-air/piano/free-air/coupling
//! is second-order, not an exact full-system propagator. All factors conserve
//! or dissipate the combined mechanical/acoustic energy, at any step size;
//! accuracy still requires rate convergence. No allocation in these flows.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_math::{c64::C64, det};

pub const MAX_AIR_STATES: usize = 1024;

#[derive(Clone, Debug)]
pub struct Pole {
    pub omega: f64,
    pub zeta: f64,
    /// Loaded-board velocity -> acoustic modal force. Signed entries matter.
    pub coupling: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct Model {
    pub ports: usize,
    pub poles: Vec<Pole>,
}
impl Model {
    pub fn validate(&self) -> Result<(), String> {
        if !(1..=super::super::linear::MAX_BOARD_MODES).contains(&self.ports)
            || self.poles.is_empty() || self.poles.len() > MAX_AIR_STATES
            || self.poles.iter().any(|p| !p.omega.is_finite() || p.omega <= 0.
                || !p.zeta.is_finite() || !(0.0..1.0).contains(&p.zeta)
                || p.coupling.len() != self.ports
                || p.coupling.iter().any(|x| !x.is_finite())) {
            return Err("invalid complete passive acoustic pole/port model".into());
        }
        Ok(())
    }
    /// Continuous impedance in the BEM's exp(-i omega t) convention. This is
    /// the physical model being split, not its finite-step numerical transfer.
    pub fn impedance(&self, omega: f64) -> Result<Vec<C64>, String> {
        self.validate()?;
        if !omega.is_finite() || omega <= 0. {return Err("invalid acoustic frequency".into());}
        let mut z = vec![C64::ZERO; self.ports*self.ports];
        for p in &self.poles {
            let denominator = C64::new(p.omega*p.omega-omega*omega, -2.*p.zeta*p.omega*omega);
            if denominator.abs() == 0. {return Err("unresolved lossless acoustic pole".into());}
            let h = C64::new(0., -omega)/denominator;
            for i in 0..self.ports {for j in 0..self.ports {
                z[i*self.ports+j] = z[i*self.ports+j] + h.scale(p.coupling[i]*p.coupling[j]);
            }}
        }
        if z.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err("acoustic impedance overflow".into());
        }
        Ok(z)
    }
}

#[derive(Clone, Copy, Default)]
struct Free { xx:f64, xp:f64, px:f64, pp:f64 }
#[derive(Clone, Copy)]
struct Rotation { air:usize, board:usize, c:f64, s:f64 }

pub struct Prepared {
    ports: usize,
    // x = omega*q_air, p = v_air, so storage is their Euclidean norm / 2.
    state: Vec<[f64;2]>,
    saved: Vec<[f64;2]>,
    free: Vec<Free>,
    rotations: Vec<Rotation>,
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
        let mut rotations=Vec::new();
        let mut row_sum=vec![0.;ports];
        let mut largest_air=0.0_f64;
        for (air,p) in model.poles.iter().enumerate() {
            largest_air=largest_air.max(p.coupling.iter().map(|x|x.abs()).sum());
            for (board,&l) in p.coupling.iter().enumerate() {
                row_sum[board]+=l.abs();
                if l==0. {continue;}
                // Forward and reverse each rotate dt/4: one coupling half-step.
                let angle=0.25*l/f64::from(rate);
                rotations.push(Rotation {air,board,c:det::cos(angle),s:det::sin(angle)});
            }
        }
        let norm_bound=largest_air.max(row_sum.into_iter().fold(0.0_f64,f64::max));
        if !norm_bound.is_finite() || norm_bound/f64::from(rate)>0.25 {
            return Err("acoustic coupling exceeds the explicit mechanical-rate resolution budget".into());
        }
        if free.iter().any(|t|[t.xx,t.xp,t.px,t.pp].iter().any(|x|!x.is_finite())) {
            return Err("nonfinite prepared acoustic transition".into());
        }
        Ok(Self {ports,state:vec![[0.;2];free.len()],saved:vec![[0.;2];free.len()],free,rotations})
    }
    pub fn energy(&self)->f64 {0.5*self.state.iter().map(|s|s[0]*s[0]+s[1]*s[1]).sum::<f64>()}
    pub fn checkpoint(&mut self) {self.saved.copy_from_slice(&self.state);}
    pub fn restore(&mut self) {self.state.copy_from_slice(&self.saved);}
    fn exchange(&mut self,board:&mut[f64]) {
        assert_eq!(board.len(),self.ports);
        for r in self.rotations.iter().chain(self.rotations.iter().rev()) {
            let b=board[r.board];let a=self.state[r.air][1];
            board[r.board]=r.c*b-r.s*a;
            self.state[r.air][1]=r.s*b+r.c*a;
        }
    }
    fn free_half(&mut self)->f64 {
        let before=self.energy();
        for (s,t) in self.state.iter_mut().zip(&self.free) {
            let [x,p]=*s;*s=[t.xx*x+t.xp*p,t.px*x+t.pp*p];
        }
        // Signed floating-point loss, not a clipped energy-balance correction.
        before-self.energy()
    }
    pub fn before(&mut self,board:&mut[f64])->f64 {self.exchange(board);self.free_half()}
    pub fn after(&mut self,board:&mut[f64])->f64 {let loss=self.free_half();self.exchange(board);loss}
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
            for _ in 0..10000 {loss+=air.before(&mut v);loss+=air.after(&mut v);}
            let final_energy=air.energy()+0.5*(v[0]*v[0]+v[1]*v[1]);
            assert!((final_energy+loss-initial).abs()<2e-13);
            assert!(air.energy()>1e-10);assert!(loss>=-1e-12);
            if damping>0. {assert!(loss>1e-8);}
        }
    }
    #[test]
    fn refused_frame_can_restore_air_history_exactly() {
        let mut air=Prepared::new(&model(0.2),192_000,2).unwrap();let mut v=[0.1,0.2];
        air.before(&mut v);air.after(&mut v);air.checkpoint();let old=air.state.clone();
        let saved=v;let loss=air.before(&mut v)+air.after(&mut v);let candidate=air.state.clone();
        air.restore();assert_eq!(air.state,old);v=saved;
        let retried=air.before(&mut v)+air.after(&mut v);
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
