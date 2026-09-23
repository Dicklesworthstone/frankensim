//! Cold positive-real matrix approximation of the ACTUAL modal BEM load.
//! A fixed damped-resonator dictionary and PSD residue matrices preserve all
//! signed cross-port terms. This is a restricted passive fit, not general MIMO
//! vector fitting. Inadequate data/basis refuses; no entrywise repair or gain EQ.
use super::{engine::radiation::{Model,Pole}, exterior_geometry::{Boundary,Specification,Samples}, exterior_loading};
use fs_math::{c64::C64,det};
const TERMS:usize=32;
const MAX_PORTS:usize=32;
const SWEEPS:usize=256;

pub struct Fit {
    pub model:Model,
    pub peak_error:f64,
    pub rms_error:f64,
    pub resistance_peak_error:f64,
    pub resistance_rms_error:f64,
}
fn dictionary(low:f64,high:f64)->Vec<(f64,f64)> {
    let mut poles:Vec<_>=(0..TERMS-1).map(|j| {
        let w=det::exp(det::ln(low/4.)+(det::ln(8.*high)-det::ln(low/4.))*j as f64/(TERMS-2) as f64);
        (w,0.2)
    }).collect();
    // Off-band lossless storage permits compact dipole loads without imposing
    // an artificial minimum resistance. Its pole must still resolve at runtime.
    poles.push((8.*high,0.));poles
}
fn normalization(p:(f64,f64))->f64 {if p.1==0. {p.0} else {2.*p.1*p.0}}
fn kernel(w:f64,p:(f64,f64))->C64 {
    let (o,z)=p;C64::new(0.,-normalization(p)*w)/C64::new(o*o-w*w,-2.*z*o*w)
}
fn norm(row:&[C64])->f64 {row.iter().fold(0.0_f64,|s,z|s.hypot(z.abs()))}
fn identity(n:usize)->Vec<f64> {let mut a=vec![0.;n*n];for i in 0..n {a[i*n+i]=1.;}a}
fn psd(a:&[f64],n:usize,id:&[f64])->Result<Vec<f64>,String> {
    if a.iter().any(|x|!x.is_finite()) {return Err("passive load fit overflow".into());}
    if n==1 {return Ok(vec![a[0].max(0.)]);}
    let scale=a.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
    if scale==0. {return Ok(vec![0.;n*n]);}
    let pairs=fs_modal::eigh_gen_dense(a,id,n).map_err(|e|e.to_string())?;
    let mut out=vec![0.;n*n];
    for p in pairs {
        if !p.lambda.is_finite() || !p.residual.is_finite() || p.residual>1e-8*scale*n as f64
            || p.phi.len()!=n || p.phi.iter().any(|x|!x.is_finite()) {
            return Err("passive fit PSD projection did not resolve its eigenproblem".into());
        }
        // Projection of a FIT ITERATE onto the PSD cone, not alteration of BEM
        // evidence. The original complex matrices remain the validation target.
        let l=p.lambda.max(0.);
        for i in 0..n {for j in 0..n {out[i*n+j]+=l*p.phi[i]*p.phi[j];}}
    }
    Ok(out)
}

/// Fit EVEN samples only. Odd samples are never used for coefficients,
/// convergence stopping, dictionary choice or tolerance selection.
pub fn fit(omega:&[f64],data:&[Vec<C64>],ports:usize)->Result<Fit,String> {
    if !(1..=MAX_PORTS).contains(&ports) || !(33..=257).contains(&omega.len()) || omega.len()%2==0
        || data.len()!=omega.len() || data.iter().any(|r|r.len()!=ports*ports
            || r.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()))
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
            || (i>0 && *w<=omega[i-1])) {
        return Err("passive load fit requires a complete 1..32-port matrix on an odd 33..257-frequency grid".into());
    }
    let poles=dictionary(omega[0],omega[omega.len()-1]);let size=ports*ports;
    let scale=data.iter().step_by(2).map(|r|norm(r)).fold(0.0_f64,f64::max);
    if !scale.is_finite() {return Err("acoustic load norm overflow".into());}
    if scale==0. {
        if data.iter().flatten().any(|z|z.abs()!=0.) {return Err("zero training load has a nonzero held-out response".into());}
        return Ok(Fit {model:Model {ports,poles:vec![Pole {omega:poles[0].0,zeta:poles[0].1,coupling:vec![0.;ports]}]},
            peak_error:0.,rms_error:0.,resistance_peak_error:0.,resistance_rms_error:0.});
    }
    // Balance resistive and reactive error using TRAINING data only. A very
    // small physical radiation resistance must not disappear in an inertance-
    // dominated complex least-squares objective. Blocks still have a scalar
    // Hessian, so the PSD projection remains their exact minimizer.
    let real_scale=data.iter().step_by(2).map(|row|row.iter().fold(0.0_f64,|s,z|s.hypot(z.re)))
        .fold(0.0_f64,f64::max);
    let real_weight=if real_scale==0. {1.} else {scale/real_scale};
    if !real_weight.is_finite() || real_weight>1e12 {
        return Err("radiation resistance is below the load fit's relative resolution budget".into());
    }
    let mut gram=vec![0.;TERMS*TERMS];let mut rhs=vec![vec![0.;size];TERMS];
    for f in (0..omega.len()).step_by(2) {
        let h:Vec<_>=poles.iter().map(|&p|kernel(omega[f],p)).collect();
        for j in 0..TERMS {
            for k in 0..TERMS {
                gram[j*TERMS+k]+=h[j].re*real_weight*h[k].re*real_weight+h[j].im*h[k].im;
            }
            for (a,value) in rhs[j].iter_mut().enumerate() {
                let y=data[f][a].scale(1./scale);
                *value+=h[j].re*real_weight*y.re*real_weight+h[j].im*y.im;
            }
        }
    }
    if gram.iter().any(|x|!x.is_finite()) || (0..TERMS).any(|j|gram[j*TERMS+j]<=0.) {
        return Err("unresolved passive fit dictionary".into());
    }
    let id=identity(ports);let mut residue=vec![vec![0.;size];TERMS];
    // Exact block minimization: for fixed other blocks, the quadratic Hessian
    // is a positive scalar times identity. Symmetric PSD projection solves
    // this block's constrained least-squares problem, using fs-modal's owner.
    for _ in 0..SWEEPS {
        let mut change=0.0_f64;let mut magnitude=0.0_f64;
        for j in 0..TERMS {
            let mut candidate=vec![0.;size];
            for (a,value) in candidate.iter_mut().enumerate() {
                *value=(rhs[j][a]-(0..TERMS).filter(|k|*k!=j)
                    .map(|k|gram[j*TERMS+k]*residue[k][a]).sum::<f64>())/gram[j*TERMS+j];
            }
            for i in 0..ports {for k in 0..i {
                let v=0.5*(candidate[i*ports+k]+candidate[k*ports+i]);
                candidate[i*ports+k]=v;candidate[k*ports+i]=v;
            }}
            let projected=psd(&candidate,ports,&id)?;
            for (old,new) in residue[j].iter_mut().zip(projected) {
                change=change.max((*old-new).abs());magnitude=magnitude.max(new.abs());*old=new;
            }
        }
        if change<=1e-10*(1.+magnitude) {break;}
    }
    let mut model=Model {ports,poles:Vec::new()};
    for (matrix,&(o,z)) in residue.iter().zip(&poles) {
        let pairs=fs_modal::eigh_gen_dense(matrix,&id,ports).map_err(|e|e.to_string())?;
        let tolerance=1e-8*ports as f64*matrix.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
        for p in pairs {
            if !p.lambda.is_finite() || !p.residual.is_finite() || p.residual>tolerance
                || p.lambda < -tolerance || p.phi.len()!=ports || p.phi.iter().any(|x|!x.is_finite()) {
                return Err("passive residue factorization failed".into());
            }
            if p.lambda<=0. {continue;}
            let gain=det::sqrt(p.lambda*scale*normalization((o,z)));
            model.poles.push(Pole {omega:o,zeta:z,coupling:p.phi.iter().map(|x|x*gain).collect()});
        }
    }
    model.validate()?;
    let mut peak=0.0_f64;let mut rms=0.0_f64;let mut rp=0.0_f64;let mut rr=0.0_f64;
    for offset in 0..2 {
        let mut max_reference=0.0_f64;let mut max_error=0.0_f64;
        let mut references=0.0_f64;let mut errors=0.0_f64;
        let mut real_max=0.0_f64;let mut real_error=0.0_f64;
        let mut real_references=0.0_f64;let mut real_errors=0.0_f64;
        for f in (offset..omega.len()).step_by(2) {
            let fitted=model.impedance(omega[f])?;
            let error:Vec<_>=fitted.iter().zip(&data[f]).map(|(a,b)|*a-*b).collect();
            let reference=norm(&data[f]);let delta=norm(&error);
            max_reference=max_reference.max(reference);max_error=max_error.max(delta);
            references=references.hypot(reference);errors=errors.hypot(delta);
            // Hermitian part controls power for COMPLEX multiport velocities.
            // Keeping its own bound prevents the much larger inertive part
            // from hiding an inaccurate radiation resistance.
            let mut h=0.0_f64;let mut dh=0.0_f64;
            for i in 0..ports {for j in 0..ports {
                h=h.hypot((data[f][i*ports+j]+data[f][j*ports+i].conj()).scale(0.5).abs());
                dh=dh.hypot((error[i*ports+j]+error[j*ports+i].conj()).scale(0.5).abs());
            }}
            real_max=real_max.max(h);real_error=real_error.max(dh);
            real_references=real_references.hypot(h);real_errors=real_errors.hypot(dh);
        }
        let ratio=|error:f64,reference:f64|if error==0. {0.} else {error/reference};
        peak=peak.max(ratio(max_error,max_reference));rms=rms.max(ratio(errors,references));
        rp=rp.max(ratio(real_error,real_max));rr=rr.max(ratio(real_errors,real_references));
    }
    if ![peak,rms,rp,rr].iter().all(|x|x.is_finite()) || peak>0.05 || rms>0.02 || rp>0.10 || rr>0.05 {
        return Err(format!("passive BEM load fit refused: complex peak/RMS={peak:.5}/{rms:.5}, power peak/RMS={rp:.5}/{rr:.5}; limits .05/.02/.10/.05"));
    }
    Ok(Fit {model,peak_error:peak,rms_error:rms,resistance_peak_error:rp,resistance_rms_error:rr})
}

/// Solve one modal BEM batch per frequency. Those SAME fields feed the load
/// and receiver fits. Pressure/acceleration = pressure/velocity * i/omega.
pub fn prepare(boundary:&Boundary,spec:&Specification,mechanics_rate:u32)->Result<(Fit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if r>MAX_PORTS || omega.len()<33 {return Err("loaded playback requires <=32 complete board modes and >=33 frequencies".into());}
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=exterior_loading::sample(boundary,spec,w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fitted=fit(&omega,&matrices,r)?;
    let delays_s=spec.receivers.iter().map(|p| {
        let radius=(0..3).fold(0.0_f64,|a,j|a.hypot(p[j]-boundary.center[j]));
        (radius-boundary.radius)/spec.medium.sound_speed
    }).collect();
    Ok((fitted,Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn grid()->Vec<f64> {(0..65).map(|i|std::f64::consts::TAU*(40.+360.*i as f64/64.)).collect()}
    #[test]
    fn signed_cross_port_resonances_fit_without_entrywise_passivity() {
        let omega=grid();let dictionary=dictionary(omega[0],omega[64]);
        let rows=[[1.,-0.7],[0.3,1.]];
        let data:Vec<Vec<_>>=omega.iter().map(|&w|(0..4).map(|a| {
            kernel(w,dictionary[8]).scale(rows[0][a/2]*rows[0][a%2])
                +kernel(w,dictionary[20]).scale(rows[1][a/2]*rows[1][a%2])
        }).collect()).collect();
        let fitted=fit(&omega,&data,2).unwrap();assert!(fitted.rms_error<0.01);
        assert!(fitted.model.poles.iter().any(|p|p.coupling[0]*p.coupling[1]<0.));
        for &w in &[1.,200.,20000.,1e6] {
            let z=fitted.model.impedance(w).unwrap();assert!((z[1]-z[2]).abs()<1e-12);
            assert!(z[0].re>=0. && z[3].re>=0.);
            assert!(z[0].re*z[3].re-z[1].re*z[2].re>=-1e-12);
        }
    }
    #[test]
    fn analytic_spherical_relaxation_keeps_inertance_and_resistance() {
        let w=grid();let data:Vec<_>=w.iter().map(|&w|vec![C64::new(0.,-w)/C64::new(6800.,-w)]).collect();
        let fit=fit(&w,&data,1).unwrap();assert!(fit.rms_error<0.02);assert!(fit.resistance_rms_error<0.05);
        assert!(fit.model.impedance(w[0]).unwrap()[0].im<0.);
    }
    #[test]
    fn weak_dipole_resistance_is_not_hidden_by_the_reactive_norm() {
        let w=grid();
        // Passive rational dipole coupon, not a fitted instrument recording.
        for a in [6800.,34000.] {
            let data:Vec<_>=w.iter().map(|&w| {
                let s=C64::new(0.,-w);
                vec![s*(s+C64::new(a,0.))/(s*s+s.scale(2.*a)+C64::new(2.*a*a,0.))]
            }).collect();
            let fit=fit(&w,&data,1).unwrap();
            assert!(fit.resistance_rms_error<0.05);assert!(fit.rms_error<0.02);
            assert!(fit.model.poles.iter().any(|p|p.zeta==0.));
        }
    }
    #[test]
    fn held_out_active_and_missing_loads_cannot_be_silently_repaired() {
        let w=grid();let mut data:Vec<_>=w.iter().map(|&w|vec![C64::new(0.,-w)/C64::new(6800.,-w)]).collect();
        for f in (1..data.len()).step_by(2) {data[f][0]=data[f][0].scale(2.);}
        assert!(fit(&w,&data,1).is_err());
        let active=vec![vec![C64::new(-1.,0.)];w.len()];assert!(fit(&w,&active,1).is_err());
        let active_cross=vec![vec![C64::ONE,C64::new(2.,0.),C64::new(2.,0.),C64::ONE];w.len()];
        assert!(fit(&w,&active_cross,2).is_err()); // positive diagonals do NOT establish multiport passivity
        assert!(fit(&w,&data,2).is_err());data[0][0]=C64::new(f64::NAN,0.);assert!(fit(&w,&data,1).is_err());
        let zero=vec![vec![C64::ZERO];w.len()];assert_eq!(fit(&w,&zero,1).unwrap().rms_error,0.);
    }
}
