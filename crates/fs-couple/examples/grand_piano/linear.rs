//! Reciprocal moving-boundary string/soundboard mechanics.
//!
//! With y(x)=sum(phi_n(x) q_n)+x/L*b, the string has fixed-interface
//! modal mass 1, cross inertia beta_n=sqrt(2 mu L)(-1)^(n+1)/(n pi),
//! and endpoint mass mu L/3. Completing the kinetic square with
//! z_n=q_n+beta_n*b leaves positive residual endpoint mass
//! mu L/3-sum(beta_n^2). fs-modal mass-normalizes the small loaded board.
//!
//! The resulting potential is sum(w_n^2 (z_n-beta_n*b)^2)/2 plus the
//! bare board and T/L endpoint stiffness. Thus string and board forces
//! are RECIPROCAL, not a gyroscopic/velocity-coupling substitute.
//!
//! Existing fs-couple exact-ZOH transitions advance each diagonal oscillator.
//! The cross-potential uses average-displacement forces; a small Schur solve
//! closes their work exactly. Linear component damping is dissipative. The
//! coupling is second-order consistent, NOT the exact full coupled propagator.
//! No new oscillator, eigensolver or matrix factorization is implemented here.

use std::f64::consts::{PI, TAU};
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::{c64::C64, det};
use super::geometry::Course;

#[derive(Clone, Debug)]
pub struct BoardMode {
    pub frequency_hz: f64,
    pub damping_ratio: f64,
    /// Mass-normalized displacement at each key's bridge, 1/sqrt(kg).
    pub bridge: [f64; 88],
    /// Integral of the surface mode, m^2/sqrt(kg); not a pressure transfer.
    pub volume: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct StringMode {
    pub omega: f64,
    pub beta: f64,
    pub a: f64,
    pub hammer_shape: f64,
    pub damper_shape: f64,
    pub string: usize,
}

#[derive(Clone, Debug)]
pub struct StringPort {
    pub course: usize,
    pub modes: std::ops::Range<usize>,
    pub bridge: Vec<f64>,
    pub hammer_lift: f64,
    pub damper_lift: f64,
    pub contact: Option<usize>,
}

#[derive(Clone, Copy, Debug, Default)]
struct Transition { qq: f64, qv: f64, vq: f64, vv: f64, bq: f64, bv: f64 }

/// Cold lowering of the EXISTING public exact-ZOH implementation to a reusable
/// linear map. Unit-energy basis probes avoid extreme-frequency energy limits.
fn transition_chunk(rate: u32, modes: &[ModalAcousticMode]) -> Result<Vec<Transition>, String> {
    let mut model = ModalAcousticTimeModel::try_new(rate, modes.to_vec(),
        ModalAcousticTimeBudget::audible_reference()).map_err(|e| e.to_string())?;
    let zeros = vec![0.0; modes.len()];
    let mut result = vec![Transition::default(); modes.len()];
    let mut states: Vec<ModalAcousticState> = modes.iter().map(|m| ModalAcousticState {
        displacement_m_sqrt_kg: 1.0 / m.angular_frequency_rad_s,
        velocity_m_sqrt_kg_per_s: 0.0,
    }).collect();
    model.restore_states(&states).map_err(|e| e.to_string())?;
    model.step(&zeros).map_err(|e| e.to_string())?;
    for ((t, state), mode) in result.iter_mut().zip(model.states()).zip(modes) {
        t.qq = state.displacement_m_sqrt_kg * mode.angular_frequency_rad_s;
        t.vq = state.velocity_m_sqrt_kg_per_s * mode.angular_frequency_rad_s;
    }
    states.fill(ModalAcousticState { displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: 1.0 });
    model.restore_states(&states).map_err(|e| e.to_string())?;
    model.step(&zeros).map_err(|e| e.to_string())?;
    for (t, state) in result.iter_mut().zip(model.states()) {
        t.qv = state.displacement_m_sqrt_kg;
        t.vv = state.velocity_m_sqrt_kg_per_s;
    }
    states.fill(ModalAcousticState::default());
    model.restore_states(&states).map_err(|e| e.to_string())?;
    model.step(&vec![1.0; modes.len()]).map_err(|e| e.to_string())?;
    for (t, state) in result.iter_mut().zip(model.states()) {
        t.bq = state.displacement_m_sqrt_kg;
        t.bv = state.velocity_m_sqrt_kg_per_s;
        if !t.bq.is_finite() || t.bq <= 0.0 { return Err("nonpositive ZOH displacement compliance".into()); }
    }
    Ok(result)
}

fn transitions(rate: u32, modes: &[ModalAcousticMode]) -> Result<Vec<Transition>, String> {
    let mut result = Vec::with_capacity(modes.len());
    for chunk in modes.chunks(4096) { result.extend(transition_chunk(rate, chunk)?); }
    Ok(result)
}

fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n { a[i * n + i] = 1.0; }
    a
}

/// All dense inversion is cold, through the existing symmetric eigenfacility.
fn inverse_spd(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let modes = fs_modal::eigh_gen_dense(a, &identity(n), n).map_err(|e| e.to_string())?;
    let mut inverse = vec![0.0; n * n];
    for m in modes {
        if !m.lambda.is_finite() || m.lambda <= 0.0 || m.residual > 1e-7 * m.lambda {
            return Err("bridge Schur complement is not resolved positive definite".into());
        }
        for i in 0..n { for j in 0..n { inverse[i*n+j] += m.phi[i]*m.phi[j]/m.lambda; } }
    }
    Ok(inverse)
}

pub struct Bank {
    pub strings: Vec<StringPort>,
    pub modes: Vec<StringMode>,
    pub contact_strings: Vec<usize>,
    pub board_count: usize,
    pub q: Vec<f64>,
    pub v: Vec<f64>,
    pub next_q: Vec<f64>,
    pub next_v: Vec<f64>,
    /// Symmetric contact displacement compliance, m/N, row-major.
    pub contact_compliance: Vec<f64>,
    pub free_contact: Vec<f64>,
    pub rate: u32,
    pub omitted_duplex_modes: usize,
    transition: Vec<Transition>,
    diagonal_omega2: Vec<f64>,
    pub last_modal_loss_j: f64,
    physical_board_k: Vec<f64>,
    board_volume: Vec<f64>,
    /// Columns map loaded coordinates to the supplied bare-board coordinates.
    board_basis: Vec<f64>,
    schur_inverse: Vec<f64>,
    contact_board: Vec<f64>,
    free_q: Vec<f64>,
    free_v: Vec<f64>,
    r_string: Vec<f64>,
    board_rhs: Vec<f64>,
    board_end: Vec<f64>,
}

impl Bank {
    /// `rate` is the mechanics rate; `band_hz` is the OUTPUT band, so
    /// oversampling never silently admits inaudible/aliased retained modes.
    pub fn new(courses: &[Course], board: &[BoardMode], rate: u32, band_hz: f64,
        max_modes: usize, damping: bool) -> Result<Self, String> {
        let r = board.len();
        if courses.is_empty() || courses.len() > 88 || !(1..=32).contains(&r)
            || !(1..=128).contains(&max_modes) || rate < 8_000
            || !band_hz.is_finite() || band_hz <= 0.0 || band_hz > 0.45*f64::from(rate) {
            return Err("invalid course, modal, frequency or sample-rate budget".into());
        }
        for b in board {
            if !b.frequency_hz.is_finite() || b.frequency_hz <= 0.0
                || !b.damping_ratio.is_finite() || b.damping_ratio < 0.0
                || !b.volume.is_finite() || b.bridge.iter().any(|x| !x.is_finite()) {
                return Err("invalid mass-normalized soundboard mode".into());
            }
        }
        let mut strings = Vec::new();
        let mut modes = Vec::new();
        let mut contact_strings = Vec::new();
        let mut oscillator = Vec::new();
        let mut mass = identity(r);
        let mut bare_k = vec![0.0; r*r];
        for (i, b) in board.iter().enumerate() { bare_k[i*r+i] = (TAU*b.frequency_hz).powi(2); }
        let mut endpoint_k = bare_k;
        let mut loaded_add = vec![0.0; r*r];
        // Authored viscoelastic bending spectrum, not measured music wire.
        // Its loss is weighted by the actual bending/tension energy fraction.
        let bending = GeneralizedMaxwell::new(200e9, vec![(8e9, 0.0004), (2e9, 0.02)])
            .map_err(|e| e.to_string())?;
        let mut omitted_duplex_modes = 0;
        for (ci, c) in courses.iter().enumerate() {
            c.validate()?;
            for member in 0..c.unison {
                let cents = (member as f64 - 0.5*(c.unison-1) as f64)*c.detune_cents;
                let tension = c.tension_at_cents(cents)?;
                for duplex in [false, true] {
                    if duplex && c.duplex_length_m == 0.0 { continue; }
                    let card = if duplex { Course { length_m: c.duplex_length_m, ..*c } } else { *c };
                    let start = modes.len();
                    let si = strings.len();
                    let mut beta2 = 0.0;
                    let mut hammer_lift = c.strike_fraction;
                    let mut damper_lift = 0.35; // authored station, replace for measured dampers
                    let bridge: Vec<f64> = board.iter().map(|b| b.bridge[usize::from(c.midi-21)]).collect();
                    for n in 1..=max_modes {
                        let f = card.partial_hz(n, tension);
                        if f > band_hz { break; }
                        let omega = TAU*f;
                        let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
                        let beta = det::sqrt(2.0*card.linear_density_kg_m*card.length_m)*sign/(n as f64*PI);
                        let hammer_shape = card.strike_shape(n);
                        let damper_shape = det::sin(n as f64*PI*0.35)/det::sqrt(card.modal_mass_kg());
                        let k = n as f64*PI/card.length_m;
                        let bend_fraction = card.flexural_rigidity_nm2*k*k/(tension + card.flexural_rigidity_nm2*k*k);
                        let zeta = if damping { 0.30/omega + 0.5*bending.loss_factor(omega)*bend_fraction } else { 0.0 };
                        oscillator.push(ModalAcousticMode { angular_frequency_rad_s: omega,
                            damping_ratio: zeta, pressure_per_modal_velocity: C64::new(0.0,0.0) });
                        modes.push(StringMode { omega, beta, a: omega*omega*beta,
                            hammer_shape, damper_shape, string: si });
                        beta2 += beta*beta;
                        hammer_lift -= hammer_shape*beta;
                        damper_lift -= damper_shape*beta;
                        for i in 0..r { for j in 0..r {
                            loaded_add[i*r+j] += omega*omega*beta*beta*bridge[i]*bridge[j];
                        } }
                    }
                    if modes.len() == start {
                        if !duplex { return Err(format!("key {} fundamental exceeds the output band",c.midi)); }
                        omitted_duplex_modes += 1;
                    }
                    let residual_mass = card.linear_density_kg_m*card.length_m/3.0 - beta2;
                    if !residual_mass.is_finite() || residual_mass <= 0.0 { return Err("invalid residual string endpoint mass".into()); }
                    for i in 0..r { for j in 0..r {
                        mass[i*r+j] += residual_mass*bridge[i]*bridge[j];
                        endpoint_k[i*r+j] += tension/card.length_m*bridge[i]*bridge[j];
                    } }
                    let contact = if duplex { None } else {
                        let index = contact_strings.len(); contact_strings.push(si); Some(index)
                    };
                    strings.push(StringPort { course: ci, modes: start..modes.len(), bridge,
                        hammer_lift, damper_lift, contact });
                }
            }
        }
        let loaded: Vec<f64> = endpoint_k.iter().zip(&loaded_add).map(|(x,y)| x+y).collect();
        let eig = fs_modal::eigh_gen_dense(&loaded, &mass, r).map_err(|e| e.to_string())?;
        for e in &eig {
            if !e.lambda.is_finite() || e.lambda <= 0.0 || e.residual > 1e-7*e.lambda
                || det::sqrt(e.lambda)/TAU > band_hz {
                return Err("loaded soundboard mode unresolved or above output band".into());
            }
        }
        let mut board_basis = vec![0.0; r*r];
        for i in 0..r { for j in 0..r { board_basis[i*r+j] = eig[j].phi[i]; } }
        for s in &mut strings {
            s.bridge = eig.iter().map(|e| s.bridge.iter().zip(&e.phi).map(|(g,p)| g*p).sum()).collect();
        }
        let mut physical_board_k = vec![0.0; r*r];
        for i in 0..r { for j in 0..r { for a in 0..r { for b in 0..r {
            physical_board_k[i*r+j] += eig[i].phi[a]*endpoint_k[a*r+b]*eig[j].phi[b];
        } } } }
        let board_volume: Vec<f64> = eig.iter().map(|e| board.iter().zip(&e.phi).map(|(b,p)| b.volume*p).sum()).collect();
        for e in &eig {
            // Loaded-coordinate modal damping is an authored reduction. It
            // is not an exact transform of the original nonproportional C.
            let zeta = if damping { board.iter().map(|b| b.damping_ratio).sum::<f64>()/r as f64 } else { 0.0 };
            oscillator.push(ModalAcousticMode { angular_frequency_rad_s: det::sqrt(e.lambda),
                damping_ratio: zeta, pressure_per_modal_velocity: C64::new(0.0,0.0) });
        }
        let transition = transitions(rate, &oscillator)?;
        let n = modes.len();
        let nc = contact_strings.len();
        let mut schur = vec![0.0; r*r];
        for i in 0..r { schur[i*r+i] = 1.0/transition[n+i].bq; }
        for (m,t) in modes.iter().zip(&transition) {
            let g = &strings[m.string].bridge;
            for i in 0..r { for j in 0..r { schur[i*r+j] -= 0.25*m.a*m.a*t.bq*g[i]*g[j]; } }
        }
        let schur_inverse = inverse_spd(&schur,r)?;
        let mut contact_board = vec![0.0; nc*r];
        let mut contact_compliance = vec![0.0; nc*nc];
        for (c,&si) in contact_strings.iter().enumerate() {
            let s = &strings[si];
            let mut lift = s.hammer_lift;
            for k in s.modes.clone() {
                lift += 0.5*modes[k].a*transition[k].bq*modes[k].hammer_shape;
                contact_compliance[c*nc+c] += transition[k].bq*modes[k].hammer_shape.powi(2);
            }
            for j in 0..r { contact_board[c*r+j] = lift*s.bridge[j]; }
        }
        for i in 0..nc { for j in 0..nc { for a in 0..r { for b in 0..r {
            contact_compliance[i*nc+j] += contact_board[i*r+a]*schur_inverse[a*r+b]*contact_board[j*r+b];
        } } } }
        Ok(Self { strings,modes,contact_strings,board_count:r,q:vec![0.0;n+r],v:vec![0.0;n+r],
            next_q:vec![0.0;n+r],next_v:vec![0.0;n+r],contact_compliance,
            free_contact:vec![0.0;nc],rate,omitted_duplex_modes,transition,physical_board_k,
            diagonal_omega2:oscillator.iter().map(|m|m.angular_frequency_rad_s.powi(2)).collect(),last_modal_loss_j:0.0,
            board_volume,board_basis,schur_inverse,contact_board,free_q:vec![0.0;n+r],free_v:vec![0.0;n+r],
            r_string:vec![0.0;n],board_rhs:vec![0.0;r],board_end:vec![0.0;r] })
    }

    /// Cold projection of a physical bare-board shape into the SAME loaded
    /// coordinates used for bridge work. This preserves modal normalization for
    /// spatial microphones; a unit-mass or frequency-only observer is not used.
    pub fn project_board_shape(&self, shape: &[f64]) -> Result<Vec<f64>, String> {
        let r=self.board_count;
        if shape.len()!=r || shape.iter().any(|x| !x.is_finite()) {
            return Err("surface projection must cover every finite bare-board mode".into());
        }
        let out:Vec<f64>=(0..r).map(|j|(0..r).map(|i|shape[i]*self.board_basis[i*r+j]).sum()).collect();
        if out.iter().any(|x| !x.is_finite()) { return Err("surface projection overflow".into()); }
        Ok(out)
    }

    pub fn contact_position(&self, contact: usize, q: &[f64]) -> f64 {
        let s = &self.strings[self.contact_strings[contact]];
        let n = self.modes.len();
        s.modes.clone().map(|k| self.modes[k].hammer_shape*q[k]).sum::<f64>()
            + s.hammer_lift*s.bridge.iter().zip(&q[n..]).map(|(g,x)| g*x).sum::<f64>()
    }

    /// Prepare the unforced coupled end positions. No allocations.
    pub fn predict(&mut self) {
        let n = self.modes.len(); let r = self.board_count;
        for k in 0..self.q.len() {
            let t=self.transition[k];
            self.free_q[k]=t.qq*self.q[k]+t.qv*self.v[k];
            self.free_v[k]=t.vq*self.q[k]+t.vv*self.v[k];
        }
        for j in 0..r { self.board_rhs[j]=self.free_q[n+j]/self.transition[n+j].bq; }
        for (k,m) in self.modes.iter().enumerate() {
            let g=&self.strings[m.string].bridge;
            let b0=g.iter().zip(&self.q[n..]).map(|(g,q)| g*q).sum::<f64>();
            let rs=self.free_q[k]+0.5*self.transition[k].bq*m.a*b0;
            self.r_string[k]=rs;
            for j in 0..r { self.board_rhs[j]+=0.5*m.a*g[j]*(self.q[k]+rs); }
        }
        for j in 0..r {
            self.board_end[j]=(0..r).map(|k| self.schur_inverse[j*r+k]*self.board_rhs[k]).sum();
            self.next_q[n+j]=self.board_end[j];
        }
        for (k,m) in self.modes.iter().enumerate() {
            let g=&self.strings[m.string].bridge;
            self.next_q[k]=self.r_string[k]+0.5*self.transition[k].bq*m.a
                *g.iter().zip(&self.board_end).map(|(g,q)| g*q).sum::<f64>();
        }
        for c in 0..self.contact_strings.len() { self.free_contact[c]=self.contact_position(c,&self.next_q); }
    }

    /// Finish the SAME prediction with held contact forces [N]. Candidate
    /// arrays are not committed until the nonlinear island accepts the step.
    pub fn finish(&mut self, forces: &[f64]) {
        let n=self.modes.len();let r=self.board_count;
        self.last_modal_loss_j=0.0;
        for j in 0..r {
            let mut response=0.0;
            for a in 0..r {
                let rhs=(0..forces.len()).map(|c|self.contact_board[c*r+a]*forces[c]).sum::<f64>();
                response+=self.schur_inverse[j*r+a]*rhs;
            }
            self.next_q[n+j]=self.board_end[j]+response;
        }
        for (k,m) in self.modes.iter().enumerate() {
            let s=&self.strings[m.string];
            let bbar=s.bridge.iter().enumerate().map(|(j,g)|g*0.5*(self.q[n+j]+self.next_q[n+j])).sum::<f64>();
            let contact=s.contact.map_or(0.0,|c|forces[c]*m.hammer_shape);
            let f=m.a*bbar+contact;
            self.next_q[k]=self.free_q[k]+self.transition[k].bq*f;
            self.next_v[k]=self.free_v[k]+self.transition[k].bv*f;
            self.last_modal_loss_j+=f*(self.next_q[k]-self.q[k])-0.5*(
                self.next_v[k].powi(2)-self.v[k].powi(2)
                +self.diagonal_omega2[k]*(self.next_q[k].powi(2)-self.q[k].powi(2)));
        }
        self.board_rhs.fill(0.0);
        for (k,m) in self.modes.iter().enumerate() {
            let g=&self.strings[m.string].bridge;
            for j in 0..r { self.board_rhs[j]+=m.a*g[j]*0.5*(self.q[k]+self.next_q[k]); }
        }
        for (c,&si) in self.contact_strings.iter().enumerate() {
            let s=&self.strings[si];
            for j in 0..r { self.board_rhs[j]+=s.hammer_lift*s.bridge[j]*forces[c]; }
        }
        for j in 0..r {
            let k=n+j;let f=self.board_rhs[j];
            self.next_v[k]=self.free_v[k]+self.transition[k].bv*f;
            self.last_modal_loss_j+=f*(self.next_q[k]-self.q[k])-0.5*(
                self.next_v[k].powi(2)-self.v[k].powi(2)
                +self.diagonal_omega2[k]*(self.next_q[k].powi(2)-self.q[k].powi(2)));
        }
    }

    pub fn commit(&mut self) { self.q.copy_from_slice(&self.next_q);self.v.copy_from_slice(&self.next_v); }

    pub fn energy(&self) -> f64 { self.energy_at(&self.q,&self.v) }
    pub fn energy_at(&self,q:&[f64],v:&[f64])->f64 {
        let n=self.modes.len();let r=self.board_count;
        let mut energy=0.5*v.iter().map(|x|x*x).sum::<f64>();
        for (k,m) in self.modes.iter().enumerate() {
            let b=self.strings[m.string].bridge.iter().zip(&q[n..]).map(|(g,q)|g*q).sum::<f64>();
            energy+=0.5*m.omega.powi(2)*(q[k]-m.beta*b).powi(2);
        }
        for i in 0..r { for j in 0..r {energy+=0.5*q[n+i]*self.physical_board_k[i*r+j]*q[n+j];} }
        energy
    }

    pub fn volume_velocity(&self)->f64 {
        self.board_volume.iter().zip(&self.v[self.modes.len()..]).map(|(a,v)|a*v).sum()
    }

    /// Exact dissipative rank-one velocity map at an actual displacement port.
    /// This is a viscous damper approximation, not a measured wool-pad model.
    /// `drag` is N s/m; damper/rank-one flow removes kinetic energy only.
    pub fn damp_string(&mut self,si:usize,drag:f64,dt:f64)->f64 {
        let s=&self.strings[si];let n=self.modes.len();
        let mut norm=0.0;let mut speed=0.0;
        for k in s.modes.clone(){let g=self.modes[k].damper_shape;norm+=g*g;speed+=g*self.v[k];}
        for j in 0..self.board_count {let g=s.damper_lift*s.bridge[j];norm+=g*g;speed+=g*self.v[n+j];}
        if norm==0.0 || drag==0.0 {return 0.0;}
        let change=det::expm1(-drag*norm*dt)*speed/norm;
        for k in s.modes.clone(){self.v[k]+=self.modes[k].damper_shape*change;}
        for j in 0..self.board_count{self.v[n+j]+=s.damper_lift*s.bridge[j]*change;}
        -change*speed-0.5*change*change*norm
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn bank(damped:bool)->Bank{
        let scale=super::super::geometry::demonstration_scale().unwrap();
        Bank::new(&[scale[48]],&super::super::board::demonstration(),192_000,21_600.0,12,damped).unwrap()
    }
    #[test]
    fn reciprocal_bridge_conserves_energy_without_gyroscopic_coupling(){
        let mut b=bank(false);
        for k in 0..b.q.len(){b.q[k]=1e-6*det::sin(k as f64);b.v[k]=0.01*det::cos(k as f64);}
        let initial=b.energy();let forces=vec![0.0;b.contact_strings.len()];
        for _ in 0..1000{b.predict();b.finish(&forces);b.commit();}
        assert!((b.energy()/initial-1.0).abs()<1e-8);
    }
    #[test]
    fn contact_compliance_is_reciprocal_and_force_work_is_physical(){
        let mut b=bank(true);let nc=b.contact_strings.len();
        for i in 0..nc{assert!(b.contact_compliance[i*nc+i]>0.0);for j in 0..nc{
            assert!((b.contact_compliance[i*nc+j]-b.contact_compliance[j*nc+i]).abs()<1e-15);
        }}
        let forces:Vec<f64>=(1..=nc).map(|x|x as f64).collect();
        let before=b.energy();b.predict();b.finish(&forces);
        let work=(0..nc).map(|i|forces[i]*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum::<f64>();
        let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
        assert!(defect.abs()<1e-12,"{defect:e}");
    }
    #[test]
    fn damper_port_removes_exactly_its_reported_energy(){
        let mut b=bank(false);b.v.fill(0.03);let before=b.energy();
        let loss=b.damp_string(0,0.4,1.0/192_000.0);
        assert!(loss>0.0);assert!((b.energy()+loss-before).abs()<1e-12);
    }
    #[test]
    fn surface_projection_uses_the_same_loaded_basis_as_volume_velocity(){
        let mut b=bank(false);
        let shape:Vec<f64>=super::super::board::demonstration().iter().map(|m|m.volume).collect();
        let projected=b.project_board_shape(&shape).unwrap();
        assert_eq!(projected,b.board_volume);
        b.v.fill(0.03);
        let observed=projected.iter().zip(&b.v[b.modes.len()..]).map(|(g,v)|g*v).sum::<f64>();
        assert_eq!(observed,b.volume_velocity());
        assert!(b.project_board_shape(&[]).is_err());
    }
}
