//! Frequency-domain bridge response of the existing moving-boundary bank.
//! exp(-i omega t), force amplitudes are peak (not RMS). No time stepper,
//! fitted resonance, invented string frequency, or altered playback state.
//! The same endpoint inertia completion and reciprocal cross-potential are
//! retained. Regular string coordinates are eliminated into a board-sized Schur
//! complement. Near fixed-interface poles stay in a bounded, pivoted border;
//! the full recovered equation and power balance are audited in both cases.
use super::{geometry::Course, linear::{Bank, BoardMode, MAX_BOARD_MODES}};
use fs_la::eigen_complex::lu_complex;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

/// Numerical border capacity, not a physical string or retained-mode cutoff.
const MAX_RETAINED_STRING_POLES: usize = 128;

fn finite(z:C64)->bool {z.re.is_finite() && z.im.is_finite()}
fn product(row:&[f64],x:&[C64])->C64 {
    row.iter().zip(x).fold(C64::ZERO,|s,(a,b)|s+b.scale(*a))
}

/// Cold continuous-operator image. It never advances the bank or its controls.
/// Damping is the bank's authored bending spectrum and physical bare-board
/// viscous matrix, not a diagonal approximation in the loaded board basis.
/// Hammer contact, key dampers, and nonlinear soundboard motion are excluded.
pub struct BridgeResponse {
    bank: Bank,
    endpoint_k: Vec<f64>,
    board_c: Vec<f64>,
    string_c: Vec<f64>,
    keys: Vec<u8>,
    band_hz: f64,
}
#[derive(Debug)]
pub struct Response {
    /// Loaded mass-normalized generalized displacements [m sqrt(kg)].
    pub board_displacement: Vec<C64>,
    pub string_displacement: Vec<C64>,
    /// All admitted course velocities, in the original scale order [m/s].
    pub bridge_velocity: Vec<C64>,
    /// Cycle-average powers [W]; amplitudes are peak phasors.
    pub input_w: f64,
    pub board_loss_w: f64,
    pub string_loss_w: f64,
    pub radiation_w: f64,
    pub power_defect_w: f64,
    /// Normwise backward error of the FULL recovered block equation.
    pub backward_error: f64,
    /// String coordinates retained in the pivoted border, not removed modes.
    pub retained_string_poles: usize,
}
impl BridgeResponse {
    pub fn new(courses:&[Course],board:&[BoardMode],rate:u32,band_hz:f64,
        max_modes:usize,damping:bool)->Result<Self,String> {
        let bank=Bank::new(courses,board,rate,band_hz,max_modes,damping)?;
        let r=bank.board_count;
        let mut endpoint_k=vec![0.;r*r];let mut board_c=vec![0.;r*r];
        // A bare coordinate's unit row is mapped by the bank's actual Phi.
        // This retains the physical wood damping after string mass loading.
        for (a,mode) in board.iter().enumerate() {
            let mut row=vec![0.;r];row[a]=1.;
            let row=bank.project_board_shape(&row)?;
            let omega=TAU*mode.frequency_hz;
            for i in 0..r {for j in 0..r {
                endpoint_k[i*r+j]+=omega*omega*row[i]*row[j];
                if damping {board_c[i*r+j]+=2.*mode.damping_ratio*omega*row[i]*row[j];}
            }}
        }
        // Same established material image as Bank::new, NOT a new loss law.
        // The parity regression below discriminates the complete generator.
        let bending=GeneralizedMaxwell::new(200e9,vec![(8e9,0.0004),(2e9,0.02)])
            .map_err(|e|e.to_string())?;
        let mut string_c=vec![0.;bank.modes.len()];let mut si=0;
        for (ci,course) in courses.iter().enumerate() {
            for member in 0..course.unison {
                let cents=(member as f64-0.5*(course.unison-1) as f64)*course.detune_cents;
                let tension=course.tension_at_cents(cents)?;
                for duplex in [false,true] {
                    if duplex && course.duplex_length_m==0. {continue;}
                    let length=if duplex {course.duplex_length_m}else{course.length_m};
                    let port=&bank.strings[si];
                    if port.course!=ci {return Err("harmonic string order disagrees with bank".into());}
                    for i in 0..r {for j in 0..r {
                        endpoint_k[i*r+j]+=tension/length*port.bridge[i]*port.bridge[j];
                    }}
                    for (index,k) in port.modes.clone().enumerate() {
                        let mode=bank.modes[k];let wave=(index+1) as f64*PI/length;
                        let fraction=course.flexural_rigidity_nm2*wave*wave/
                            (tension+course.flexural_rigidity_nm2*wave*wave);
                        if damping {string_c[k]=2.*mode.omega*(0.30/mode.omega
                            +0.5*bending.loss_factor(mode.omega)*fraction);}
                    }
                    si+=1;
                }
            }
        }
        if si!=bank.strings.len() || endpoint_k.iter().chain(&board_c).chain(&string_c)
            .any(|v|!v.is_finite()) {return Err("nonfinite or incomplete harmonic bank".into());}
        Ok(Self {bank,endpoint_k,board_c,string_c,keys:courses.iter().map(|c|c.midi).collect(),band_hz})
    }
    pub fn bank(&self)->&Bank {&self.bank}
    pub fn keys(&self)->&[u8] {&self.keys}
    pub fn bridge_row(&self,key:u8)->Result<&[f64],String> {
        let course=self.keys.iter().position(|k|*k==key).ok_or("drive key absent from physical scale")?;
        self.bank.strings.iter().find(|s|s.course==course).map(|s|s.bridge.as_slice())
            .ok_or_else(||"course has no physical bridge port".into())
    }
    /// Z is a row-major force/velocity radiation impedance in the SAME loaded
    /// board coordinates. None is vacuum, not a fallback after a failed solve.
    /// The applied force is at one physical bridge station and is never divided
    /// by the number of unison strings. Strings and duplexes respond passively.
    pub fn solve(&self,hz:f64,key:u8,force_n:C64,z:Option<&[C64]>)->Result<Response,String> {
        let r=self.bank.board_count;let n=self.bank.modes.len();
        if !(1..=MAX_BOARD_MODES).contains(&r) || !hz.is_finite() || hz<=0.
            || hz>self.band_hz || !finite(force_n)
            || z.is_some_and(|z|z.len()!=r*r || z.iter().any(|v|!finite(*v))) {
            return Err("invalid harmonic frequency, force or complete radiation matrix".into());
        }
        let g=self.bridge_row(key)?;let w=TAU*hz;let iw=C64::new(0.,-w);
        let force:Vec<C64>=g.iter().map(|v|force_n.scale(*v)).collect();
        let mut matrix=vec![C64::ZERO;r*r];
        for i in 0..r {for j in 0..r {
            matrix[i*r+j]=C64::new(self.endpoint_k[i*r+j]-if i==j {w*w}else{0.},
                -w*self.board_c[i*r+j]);
            if let Some(z)=z {matrix[i*r+j]=matrix[i*r+j]+iw*z[i*r+j];}
        }}
        let mut divisors=Vec::with_capacity(n);let mut poles=Vec::new();
        for (k,mode) in self.bank.modes.iter().enumerate() {
            let d=C64::new(mode.omega.powi(2)-w*w,-w*self.string_c[k]);
            if !finite(d) {return Err(format!("nonfinite string divisor at {hz} Hz, mode {k}"));}
            // A fixed-interface pole is NOT necessarily a resonance of the
            // coupled piano. Do not divide by its small diagonal; retain that
            // original coordinate and let the existing pivoted LU solve it.
            // This threshold changes only the partition, not any coefficient.
            if d.abs()<=8.*f64::EPSILON.sqrt()*(mode.omega.powi(2)+w*w) {
                if poles.len()==MAX_RETAINED_STRING_POLES {
                    return Err("harmonic response exceeds the 128-coordinate string-pole border budget".into());
                }
                poles.push(k);
            }
            divisors.push(d);
        }
        for port in &self.bank.strings {
            let mut weight=C64::ZERO;
            for k in port.modes.clone() {
                let mode=self.bank.modes[k];
                // omega_s^2 beta^2 - a^2/D, evaluated without catastrophic
                // low-frequency cancellation of the two static terms.
                weight=weight+if poles.binary_search(&k).is_ok() {
                    // Keep this coordinate's full board self-term; its cross
                    // reaction is solved explicitly in the bordered equation.
                    C64::new(mode.omega.powi(2)*mode.beta.powi(2),0.)
                } else {C64::new(-w*w,-w*self.string_c[k])
                    .scale(mode.omega.powi(2)*mode.beta.powi(2))/divisors[k]};
            }
            for i in 0..r {for j in 0..r {
                matrix[i*r+j]=matrix[i*r+j]+weight.scale(port.bridge[i]*port.bridge[j]);
            }}
        }
        if matrix.iter().any(|v|!finite(*v)) {return Err("harmonic Schur overflow".into());}
        let (q,pole_displacements)=if poles.is_empty() {
            // Preserve the original operation order away from string poles.
            let lu=lu_complex(&matrix,r).map_err(|_|"singular radiation-loaded bridge equation")?;
            let mut q=force.clone();lu.solve(&mut q);(q,Vec::new())
        } else {self.solve_pole_border(&matrix,&force,&divisors,&poles)?};
        if q.iter().chain(&pole_displacements).any(|v|!finite(*v)) {
            return Err("nonfinite harmonic board/string response".into());
        }
        let mut strings=vec![C64::ZERO;n];
        let mut reconstructed=vec![C64::ZERO;r];let mut row_norm=vec![0.;r];
        let mut string_loss_w=0.;let mut worst_residual=0.0_f64;let mut operator_norm=0.0_f64;
        for port in &self.bank.strings {
            let b=product(&port.bridge,&q);let abs_g:f64=port.bridge.iter().map(|v|v.abs()).sum();
            for k in port.modes.clone() {
                let m=self.bank.modes[k];let displacement=match poles.binary_search(&k) {
                    Ok(slot)=>pole_displacements[slot],Err(_)=>b.scale(m.a)/divisors[k],
                };
                if !finite(displacement) {return Err("nonfinite recovered string response".into());}
                strings[k]=displacement;
                worst_residual=worst_residual.max((divisors[k]*displacement-b.scale(m.a)).abs());
                operator_norm=operator_norm.max(divisors[k].abs()+m.a.abs()*abs_g);
                string_loss_w+=0.5*w*w*self.string_c[k]*displacement.abs().powi(2);
                let reaction=(b.scale(m.beta)-displacement).scale(m.a);
                for i in 0..r {
                    reconstructed[i]=reconstructed[i]+reaction.scale(port.bridge[i]);
                    row_norm[i]+=m.a.abs()*port.bridge[i].abs()*(1.+m.beta.abs()*abs_g);
                }
            }
        }
        let mut board_loss_w=0.;let mut radiation_w=0.;
        for i in 0..r {
            let mut loss=C64::ZERO;let mut radiation=C64::ZERO;
            for j in 0..r {
                let entry=C64::new(self.endpoint_k[i*r+j]-if i==j {w*w}else{0.},-w*self.board_c[i*r+j]);
                reconstructed[i]=reconstructed[i]+entry*q[j];row_norm[i]+=entry.abs();
                loss=loss+q[j].scale(self.board_c[i*r+j]);
                if let Some(z)=z {radiation=radiation+z[i*r+j]*q[j];row_norm[i]+=w*z[i*r+j].abs();}
            }
            reconstructed[i]=reconstructed[i]+iw*radiation;
            worst_residual=worst_residual.max((reconstructed[i]-force[i]).abs());
            operator_norm=operator_norm.max(row_norm[i]);
            board_loss_w+=0.5*w*w*(q[i].conj()*loss).re;
            radiation_w+=0.5*w*w*(q[i].conj()*radiation).re;
        }
        let bridge_velocity:Vec<_>=self.keys.iter().map(|k|
            self.bridge_row(*k).map(|g|iw*product(g,&q))).collect::<Result<_,_>>()?;
        let input_w=0.5*(force_n.conj()*iw*product(g,&q)).re;
        let power_defect_w=input_w-board_loss_w-string_loss_w-radiation_w;
        let maximum_q=q.iter().chain(&strings).map(|v|v.abs()).fold(0.0_f64,f64::max);
        let maximum_force=force.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
        let scale=operator_norm*maximum_q+maximum_force;
        if !scale.is_finite() {return Err("harmonic residual scale overflow".into());}
        let backward_error=if scale==0. {worst_residual}else{worst_residual/scale};
        let power_scale=input_w.abs()+board_loss_w.abs()+string_loss_w.abs()+radiation_w.abs();
        let tolerance=1e-10+1e-7*power_scale;
        if !backward_error.is_finite() || backward_error>1e-9
            || [input_w,board_loss_w,string_loss_w,radiation_w,power_defect_w].iter().any(|v|!v.is_finite())
            || radiation_w < -tolerance || board_loss_w < -tolerance || string_loss_w < -tolerance
            || power_defect_w.abs()>tolerance {
            return Err(format!("bridge response failed full-equation/power admission: backward={backward_error}, defect={power_defect_w} W, radiation={radiation_w} W"));
        }
        Ok(Response {board_displacement:q,string_displacement:strings,bridge_velocity,input_w,
            board_loss_w,string_loss_w,radiation_w,power_defect_w,backward_error,
            retained_string_poles:poles.len()})
    }

    /// [board, retained string coordinates]. Regular strings are already
    /// eliminated. The exact symmetric cross-potential stays on BOTH sides;
    /// supplied radiation retains its original entries (including reciprocity
    /// or lack of it). No pseudoinverse, pole shift or artificial loss.
    fn solve_pole_border(&self, board:&[C64], force:&[C64], divisors:&[C64],
        poles:&[usize])->Result<(Vec<C64>,Vec<C64>),String> {
        let r=self.bank.board_count;let size=r+poles.len();
        let mut matrix=vec![C64::ZERO;size*size];let mut rhs=vec![C64::ZERO;size];
        for i in 0..r {matrix[i*size..i*size+r].copy_from_slice(&board[i*r..(i+1)*r]);}
        rhs[..r].copy_from_slice(force);
        for (slot,&k) in poles.iter().enumerate() {
            let row=r+slot;let mode=self.bank.modes[k];
            matrix[row*size+row]=divisors[k];
            for (i,g) in self.bank.strings[mode.string].bridge.iter().enumerate() {
                let cross=C64::new(-mode.a*g,0.);
                matrix[i*size+row]=cross;matrix[row*size+i]=cross;
            }
        }
        if matrix.iter().any(|v|!finite(*v)) {return Err("harmonic string-pole border overflow".into());}
        let lu=lu_complex(&matrix,size).map_err(|_|"singular coupled bridge/string-pole equation")?;
        lu.solve(&mut rhs);
        Ok((rhs[..r].to_vec(),rhs[r..].to_vec()))
    }
}

#[cfg(test)]
#[path="bridge_pole_tests.rs"]
mod pole_tests;

#[cfg(test)]
mod tests {
    use super::*;
    fn inputs()->(Vec<Course>,Vec<BoardMode>) {
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let courses=[48,51].map(|i|Course {unison:2,duplex_length_m:0.2,..scale[i]}).to_vec();
        let mut board=vec![BoardMode {frequency_hz:170.,damping_ratio:0.012,bridge:[0.;88],volume:0.1},
            BoardMode {frequency_hz:310.,damping_ratio:0.018,bridge:[0.;88],volume:-0.03}];
        board[0].bridge[48]=0.08;board[0].bridge[51]=0.06;
        board[1].bridge[48]=-0.04;board[1].bridge[51]=0.09;
        (courses,board)
    }
    fn model(damping:bool)->BridgeResponse {
        let (c,b)=inputs();BridgeResponse::new(&c,&b,192_000,21_600.,4,damping).unwrap()
    }
    fn load(w:f64)->Vec<C64> {
        // Authored passive two-port: non-diagonal resistance + added mass.
        vec![C64::new(3.,-w*0.01),C64::new(0.5,-w*0.002),
             C64::new(0.5,-w*0.002),C64::new(2.,-w*0.015)]
    }
    #[test]
    fn condensed_solution_matches_full_energy_hessian_and_closes_power() {
        let model=model(true);let bank=&model.bank;let n=bank.modes.len();let r=bank.board_count;let size=n+r;
        let hz=277.;let w=TAU*hz;let z=load(w);
        let answer=model.solve(hz,69,C64::new(1.,0.3),Some(&z)).unwrap();
        // Independent dense oracle: polarize the ACTUAL time bank energy, not
        // a second copy of this module's condensed stiffness construction.
        let zero=vec![0.;size];let mut q=zero.clone();let mut energies=vec![0.;size];
        for i in 0..size {q[i]=1.;energies[i]=bank.energy_at(&q,&zero);q[i]=0.;}
        let mut full=vec![C64::ZERO;size*size];
        for i in 0..size {for j in i..size {
            let value=if i==j {2.*energies[i]}else {
                q[i]=1.;q[j]=1.;let value=bank.energy_at(&q,&zero)-energies[i]-energies[j];q[i]=0.;q[j]=0.;value
            };
            full[i*size+j]=C64::new(value,0.);full[j*size+i]=C64::new(value,0.);
        }}
        for i in 0..size {full[i*size+i]=full[i*size+i]-C64::new(w*w,if i<n {w*model.string_c[i]}else{0.});}
        for i in 0..r {for j in 0..r {
            let at=(n+i)*size+n+j;
            full[at]=full[at]+C64::new(0.,-w*model.board_c[i*r+j])+C64::new(0.,-w)*z[i*r+j];
        }}
        let mut rhs=vec![C64::ZERO;size];
        for (out,g) in rhs[n..].iter_mut().zip(model.bridge_row(69).unwrap()) {*out=C64::new(1.,0.3).scale(*g);}
        lu_complex(&full,size).unwrap().solve(&mut rhs);
        let expected:Vec<_>=answer.string_displacement.iter().chain(&answer.board_displacement).copied().collect();
        let maximum=rhs.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
        for (actual,expected) in rhs.iter().zip(expected) {assert!((*actual-expected).abs()<1e-7*maximum);}
        assert!(answer.input_w>0. && answer.board_loss_w>0. && answer.string_loss_w>0. && answer.radiation_w>0.);
        assert!(answer.power_defect_w.abs()<1e-9*answer.input_w);
        assert!(answer.backward_error<1e-12);
    }
    #[test]
    fn cross_bridge_response_is_reciprocal_and_silent_strings_remain_present() {
        let model=model(true);let z=load(TAU*277.);
        let a=model.solve(277.,69,C64::ONE,Some(&z)).unwrap();
        let b=model.solve(277.,72,C64::ONE,Some(&z)).unwrap();
        assert!((a.bridge_velocity[1]-b.bridge_velocity[0]).abs()<1e-10*a.bridge_velocity[1].abs());
        let (courses,board)=inputs();
        let one=BridgeResponse::new(&courses[..1],&board,192_000,21_600.,4,true).unwrap();
        let first=model.bank.strings.iter().position(|s|s.course==1).unwrap();
        // Probe the actual silent course's fundamental, not an unrelated
        // frequency where its correctly retained reaction can be tiny.
        let hz=model.bank.modes[model.bank.strings[first].modes.start].omega/TAU;
        let one=one.solve(hz,69,C64::ONE,None).unwrap();
        let both=model.solve(hz,69,C64::ONE,None).unwrap();
        assert!((one.bridge_velocity[0]-both.bridge_velocity[0]).abs()>1e-6*both.bridge_velocity[0].abs());
        assert!(both.string_displacement[model.bank.strings[first].modes.clone()].iter().any(|v|v.abs()>0.));
        assert!(model.bank.q.iter().chain(&model.bank.v).all(|v|*v==0.));
    }
    #[test]
    fn material_loss_matches_the_existing_time_bank_generator() {
        let (mut courses,board)=inputs();courses.truncate(1);courses[0].duplex_length_m=0.;
        let mut errors=Vec::new();
        for rate in [96_000,192_000] {
            let mut m=BridgeResponse::new(&courses,&board,rate,21_600.,4,true).unwrap();
            let n=m.bank.modes.len();
            for i in 0..n {m.bank.v[i]=0.01*((i+1) as f64).sin();}
            m.bank.v[n]=0.006;m.bank.v[n+1]=-0.004;
            let mut expected:f64=(0..n).map(|i|m.string_c[i]*m.bank.v[i].powi(2)).sum();
            for i in 0..2 {for j in 0..2 {expected+=m.board_c[2*i+j]*m.bank.v[n+i]*m.bank.v[n+j];}}
            m.bank.predict();m.bank.finish(&vec![0.;m.bank.contact_strings.len()]);
            let measured=m.bank.last_modal_loss_j*f64::from(rate);
            errors.push((measured-expected).abs()/expected);
        }
        assert!(errors[1]<0.004,"continuous damping disagrees with time bank: {errors:?}");
        assert!(errors[1]<errors[0],"generator comparison must improve at finer step: {errors:?}");
    }
    #[test]
    fn refusals_do_not_shift_poles_or_invent_damping_and_zero_force_is_zero() {
        let m=model(false);let pole=m.bank.modes[0].omega/TAU;
        let at_pole=m.solve(pole,69,C64::ONE,None).unwrap();
        assert!(at_pole.retained_string_poles>0);assert!(at_pole.backward_error<1e-9);
        for hz in [0.,f64::NAN,21_601.] {assert!(m.solve(hz,69,C64::ONE,None).is_err());}
        assert!(m.solve(277.,60,C64::ONE,None).is_err());
        assert!(m.solve(277.,69,C64::ONE,Some(&[C64::ZERO])).is_err());
        let active=vec![C64::new(-1e5,0.),C64::ZERO,C64::ZERO,C64::new(-1e5,0.)];
        assert!(m.solve(277.,69,C64::ONE,Some(&active)).is_err());
        let zero=m.solve(277.,69,C64::ZERO,None).unwrap();
        assert_eq!(zero.input_w,0.);assert_eq!(zero.backward_error,0.);
        assert!(zero.bridge_velocity.iter().all(|v|*v==C64::ZERO));
    }
}
