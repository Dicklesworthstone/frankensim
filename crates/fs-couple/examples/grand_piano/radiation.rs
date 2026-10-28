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
    ModalAcousticTimeBudget, ModalAcousticTimeModel, MAX_TIME_DOMAIN_ACOUSTIC_MODES};
use fs_math::c64::C64;
use fs_phs::{PortExchangeBudget, PreparedPortExchange};

/// Every acoustic oscillator a complete passive model may carry.
pub const MAX_AIR_STATES: usize = MAX_POLES;

// Shared positive-real model; collective exchange does not alter its impedance.
pub use fs_couple::render::plate::impact::radiation::{Model,Pole,MAX_POLES,MAX_PORTS};
use fs_couple::render::plate::impact::radiation::fit::{Fit,fit};

#[path = "radiation_ports.rs"]
mod ports;
pub use ports::PortBasis;
pub(crate) use ports::{geometry_basis, validate_complete_fit};

pub(crate) fn fit_projected(omega:&[f64],data:&[Vec<C64>],basis:&PortBasis)->Result<Fit,String> {
    let projected=data.iter().map(|row|basis.project_impedance(row)).collect::<Result<Vec<_>,_>>()?;
    let mut fitted=fit(omega,&projected,basis.ports())?;
    validate_complete_fit(omega,data,basis,&mut fitted)?;
    Ok(fitted)
}

#[derive(Clone, Copy, Default)]
struct Free { xx:f64, xp:f64, px:f64, pp:f64 }

struct Projected {
    basis: PortBasis,
    incoming: Vec<f64>,
    outgoing: Vec<f64>,
    board: Vec<f64>,
    air: Vec<[f64; 2]>,
}

fn stored_energy(state:&[[f64;2]])->f64 {
    0.5*state.iter().map(|s|s[0]*s[0]+s[1]*s[1]).sum::<f64>()
}
fn free_motion(state:&mut[[f64;2]],free:&[Free])->f64 {
    let before=stored_energy(state);
    for (s,t) in state.iter_mut().zip(free) {
        let [x,p]=*s;*s=[t.xx*x+t.xp*p,t.px*x+t.pp*p];
    }
    // Signed floating-point loss, not a clipped energy-balance correction.
    before-stored_energy(state)
}

pub struct Prepared {
    ports: usize,
    // x = omega*q_air, p = v_air, so storage is their Euclidean norm / 2.
    state: Vec<[f64;2]>,
    saved: Vec<[f64;2]>,
    free: Vec<Free>,
    exchange: PreparedPortExchange,
    velocities: Vec<f64>,
    projected: Option<Projected>,
}
impl Prepared {
    pub fn new(model:&Model, rate:u32, ports:usize) -> Result<Self,String> {
        model.validate()?;
        if rate < 8_000 || ports != model.ports
            || model.poles.iter().any(|p| p.omega >= 0.9*std::f64::consts::PI*f64::from(rate)) {
            return Err("acoustic model does not cover the loaded basis or mechanical rate".into());
        }
        let half_rate=rate.checked_mul(2).ok_or("acoustic half-step rate overflow")?;
        // Each oscillator's exact transition is independent of the others, so they are
        // taken from the shared owner in groups within its own mode ceiling.
        let mut free=vec![Free::default();model.poles.len()];
        for (poles,free) in model.poles.chunks(MAX_TIME_DOMAIN_ACOUSTIC_MODES)
            .zip(free.chunks_mut(MAX_TIME_DOMAIN_ACOUSTIC_MODES)) {
            let modes:Vec<_>=poles.iter().map(|p| ModalAcousticMode {
                angular_frequency_rad_s:p.omega, damping_ratio:p.zeta,
                pressure_per_modal_velocity:C64::ZERO,
            }).collect();
            let mut oscillator=ModalAcousticTimeModel::try_new(half_rate,modes,
                ModalAcousticTimeBudget::audible_reference()).map_err(|e|e.to_string())?;
            let zero=vec![0.;poles.len()];
            let mut states:Vec<_>=poles.iter().map(|p|ModalAcousticState {
                displacement_m_sqrt_kg:1./p.omega,velocity_m_sqrt_kg_per_s:0.,
            }).collect();
            oscillator.restore_states(&states).map_err(|e|e.to_string())?;
            oscillator.step(&zero).map_err(|e|e.to_string())?;
            for ((t,s),p) in free.iter_mut().zip(oscillator.states()).zip(poles) {
                t.xx=p.omega*s.displacement_m_sqrt_kg;t.px=s.velocity_m_sqrt_kg_per_s;
            }
            states.fill(ModalAcousticState {displacement_m_sqrt_kg:0.,velocity_m_sqrt_kg_per_s:1.});
            oscillator.restore_states(&states).map_err(|e|e.to_string())?;
            oscillator.step(&zero).map_err(|e|e.to_string())?;
            for ((t,s),p) in free.iter_mut().zip(oscillator.states()).zip(poles) {
                t.xp=p.omega*s.displacement_m_sqrt_kg;t.pp=s.velocity_m_sqrt_kg_per_s;
            }
        }
        let coupling:Vec<_>=model.poles.iter().flat_map(|p|p.coupling.iter().copied()).collect();
        // One collective half-step, not separate physical-port rotations. Its
        // largest rotation, sigma_max(L) dt/2, stays within .125 per half-step.
        let exchange=PreparedPortExchange::new(ports,model.poles.len(),&coupling,
            0.5/f64::from(rate),PortExchangeBudget {max_left:MAX_PORTS,max_right:MAX_AIR_STATES,
                max_setup_terms:60*MAX_PORTS*MAX_PORTS*MAX_AIR_STATES,maximum_dt_coupling:0.125})
            .map_err(|e|e.to_string())?;
        if free.iter().any(|t|[t.xx,t.xp,t.px,t.pp].iter().any(|x|!x.is_finite())) {
            return Err("nonfinite prepared acoustic transition".into());
        }
        Ok(Self {ports,state:vec![[0.;2];free.len()],saved:vec![[0.;2];free.len()],
            velocities:vec![0.;free.len()],free,exchange,projected:None})
    }
    /// Attach a passive model through explicit acoustic combinations of the
    /// complete board. The existing exchange still sees at most 32 ports;
    /// omitted acoustic directions retain their original mechanical motion.
    /// Both board and air are staged before each half-flow is published.
    pub fn new_projected(model:&Model,rate:u32,basis:&PortBasis)->Result<Self,String> {
        if model.ports!=basis.ports() {
            return Err("radiation model and selected acoustic basis have different port counts".into());
        }
        let mut prepared=Self::new(model,rate,basis.ports())?;
        prepared.ports=basis.board_ports();
        prepared.projected=Some(Projected {basis:basis.clone(),
            incoming:vec![0.;basis.ports()],outgoing:vec![0.;basis.ports()],
            board:vec![0.;basis.board_ports()],air:vec![[0.;2];model.poles.len()]});
        Ok(prepared)
    }
    pub fn energy(&self)->f64 {stored_energy(&self.state)}
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
        free_motion(&mut self.state,&self.free)
    }
    fn projected_half(&mut self,board:&mut[f64],free_first:bool)->Result<f64,&'static str> {
        if board.len()!=self.ports || board.iter().any(|v|!v.is_finite()) {
            return Err("projected radiation requires the complete finite board velocity");
        }
        let work=self.projected.as_mut().ok_or("missing radiation port basis")?;
        work.air.copy_from_slice(&self.state);
        let mut loss=if free_first {free_motion(&mut work.air,&self.free)}else{0.};
        for (u,row) in work.incoming.iter_mut().zip(work.basis.vectors()) {
            *u=row.iter().zip(board.iter()).map(|(q,v)|q*v).sum();
        }
        work.outgoing.copy_from_slice(&work.incoming);
        for (v,s) in self.velocities.iter_mut().zip(&work.air) {*v=s[1];}
        let before=0.5*(board.iter().map(|v|v*v).sum::<f64>()
            +self.velocities.iter().map(|v|v*v).sum::<f64>());
        if !before.is_finite() || !loss.is_finite() {
            return Err("projected radiation energy overflow");
        }
        self.exchange.apply(&mut work.outgoing,&mut self.velocities)
            .map_err(|_|"collective projected radiation exchange failed finite/energy admission")?;
        // Lift the CHANGE only: (I-QQ^T)v is never discarded or re-integrated.
        for (i,(candidate,&original)) in work.board.iter_mut().zip(board.iter()).enumerate() {
            let mut change=0.;
            for ((row,&outgoing),&incoming) in work.basis.vectors().iter()
                .zip(&work.outgoing).zip(&work.incoming) {
                change+=row[i]*(outgoing-incoming);
            }
            *candidate=if change==0. {original}else{original+change};
        }
        let after=0.5*(work.board.iter().map(|v|v*v).sum::<f64>()
            +self.velocities.iter().map(|v|v*v).sum::<f64>());
        let tolerance=(512.*f64::EPSILON
            *(self.ports+self.velocities.len()+work.basis.ports()+1) as f64
            +4.*work.basis.orthogonality_error)*(before+after);
        if !after.is_finite() || !tolerance.is_finite() || (after-before).abs()>tolerance {
            return Err("projected radiation fails full board/air energy admission");
        }
        for (s,&v) in work.air.iter_mut().zip(&self.velocities) {s[1]=v;}
        if !free_first {loss=free_motion(&mut work.air,&self.free);}
        if !loss.is_finite() || !stored_energy(&work.air).is_finite() {
            return Err("projected acoustic free motion is nonfinite");
        }
        board.copy_from_slice(&work.board);
        self.state.copy_from_slice(&work.air);
        Ok(loss)
    }
    pub fn before(&mut self,board:&mut[f64])->Result<f64,&'static str> {
        if self.projected.is_some() {return self.projected_half(board,false);}
        self.exchange(board)?;Ok(self.free_half())
    }
    pub fn after(&mut self,board:&mut[f64])->Result<f64,&'static str> {
        if self.projected.is_some() {return self.projected_half(board,true);}
        let loss=self.free_half();self.exchange(board)?;Ok(loss)
    }
}

#[cfg(test)]
#[path = "radiation_projection_tests.rs"]
mod projection_tests;

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
        for bad in [f64::NAN,-1.,f64::INFINITY] {let mut m=good.clone();m.poles[0].zeta=bad;assert!(m.validate().is_err());}
        let mut m=good.clone();m.poles[0].coupling[0]=1e9;
        assert!(Prepared::new(&m,192_000,2).is_err());
        let mut m=good;m.poles[0].coupling.pop();assert!(m.validate().is_err());
    }
    // Independent real-pole impulse integrals. These do not obtain loss by
    // subtracting endpoint energies as the runtime's free flow does.
    fn exponential_integral(rate:f64,t:f64)->f64 { -(-rate*t).exp_m1()/rate }
    fn free_oracle(omega:f64,zeta:f64,t:f64)->([f64;2],f64) {
        if zeta==1. {
            let a=omega*t;let e=(-a).exp();
            return ([e*(1.+a),-a*e],0.5*(1.-(-2.*a).exp()*(1.+2.*a+2.*a*a)));
        }
        let sum=zeta+(zeta*zeta-1.).sqrt();
        let (slow,fast)=(omega/sum,omega*sum);
        let (a,b)=((-slow*t).exp(),(-fast*t).exp());
        let factor=omega/(fast-slow);
        let loss=2.*zeta*omega*factor*factor
            *(exponential_integral(2.*slow,t)-2.*exponential_integral(slow+fast,t)
                +exponential_integral(2.*fast,t));
        ([(fast*a-slow*b)/(fast-slow),factor*(b-a)],loss)
    }

    #[test]
    fn critical_and_overdamped_free_memory_match_motion_and_independent_heat() {
        // The fastest real pole in the last case is far above the sample rate.
        // It is a decay, not an oscillatory frequency; exact propagation keeps
        // the physical slow tail without an artificial Nyquist refusal.
        let (omega,rate)=(1200.,8000);
        for zeta in [1.,2.,1e6] {
            let m=Model {ports:1,poles:vec![Pole {omega,zeta,coupling:vec![0.]}]};
            let mut air=Prepared::new(&m,rate,1).unwrap();
            air.state[0]=[1.,0.];
            let loss=air.free_half();
            let (expected,heat)=free_oracle(omega,zeta,0.5/f64::from(rate));
            for (a,b) in air.state[0].iter().zip(expected) {
                assert!((a-b).abs()<2e-11,"zeta={zeta}: {a} != {b}");
            }
            assert!(loss>=0.);
            assert!((loss-heat).abs()<2e-11*heat.max(1e-4),"{zeta}: {loss} != {heat}");
        }
    }

    #[test]
    fn nonresonant_radiation_impedance_keeps_signed_reciprocal_passivity() {
        for zeta in [1.,2.,1000.] {
            let m=model(zeta);
            for w in [10.,500.,1200.,3200.,10000.] {
                let z=m.impedance(w).unwrap();assert_eq!(z[1],z[2]);
                let v=[C64::new(0.3,0.8),C64::new(-0.1,0.4)];
                let power=(0..2).map(|i|(0..2)
                    .map(|j|(v[i].conj()*z[2*i+j]*v[j]).re).sum::<f64>()).sum::<f64>();
                let expected=m.poles.iter().map(|p| {
                    let drive=p.coupling.iter().zip(v).fold(C64::ZERO,|s,(l,v)|s+v.scale(*l));
                    let resistance=2.*zeta*p.omega*w*w
                        /((p.omega*p.omega-w*w).powi(2)+(2.*zeta*p.omega*w).powi(2));
                    drive.abs().powi(2)*resistance
                }).sum::<f64>();
                assert!(power>=0.);
                assert!((power-expected).abs()<2e-13*expected.max(1e-20));
            }
        }
    }

    #[test]
    fn critical_and_overdamped_reaction_converge_to_the_complete_physical_solution() {
        let (omega,g,v0,t)=(1200.0_f64,200.0_f64,0.1,0.004);
        let loaded=(omega*omega+g*g).sqrt();
        for ratio in [1.,2.] {
            let zeta=ratio*loaded/omega;
            let force=g*v0;
            let (q,p,heat)=if ratio==1. {
                let a=loaded*t;let e=(-a).exp();
                (force/(loaded*loaded)*(1.-e*(1.+a)),force*t*e,
                    force*force/(2.*loaded*loaded)
                        *(1.-(-2.*a).exp()*(1.+2.*a+2.*a*a)))
            } else {
                let sum=ratio+(ratio*ratio-1.0_f64).sqrt();
                let (slow,fast)=(loaded/sum,loaded*sum);
                let (a,b)=((-slow*t).exp(),(-fast*t).exp());
                let factor=force/(fast-slow);
                (force/(loaded*loaded)*(1.-(fast*a-slow*b)/(fast-slow)),
                    factor*(a-b),2.*zeta*omega*factor*factor
                        *(exponential_integral(2.*slow,t)-2.*exponential_integral(slow+fast,t)
                            +exponential_integral(2.*fast,t)))
            };
            let mut previous=None;
            for rate in [48_000,96_000,192_000] {
                let m=Model {ports:1,poles:vec![Pole {omega,zeta,coupling:vec![g]}]};
                let mut air=Prepared::new(&m,rate,1).unwrap();
                let mut board=[v0];let mut loss=0.;
                for _ in 0..(t*f64::from(rate)).round() as usize {
                    loss+=air.before(&mut board).unwrap();
                    loss+=air.after(&mut board).unwrap();
                }
                let initial=0.5*v0*v0;
                assert!((air.energy()+0.5*board[0]*board[0]+loss-initial).abs()<2e-14);
                assert!(loss>0.);
                let error=((air.state[0][0]-omega*q).powi(2)
                    +(air.state[0][1]-p).powi(2)+(board[0]-(v0-g*q)).powi(2)).sqrt()/v0
                    +(loss-heat).abs()/initial;
                if let Some(old)=previous {
                    assert!((3.8..4.2).contains(&(old/error)),"zeta={zeta}, rate={rate}: {old} -> {error}");
                }
                eprintln!("radiation zeta={zeta:.6}, rate={rate}: response/work error={error:e}");
                previous=Some(error);
            }
            assert!(previous.unwrap()<1e-5);
        }
    }

}
