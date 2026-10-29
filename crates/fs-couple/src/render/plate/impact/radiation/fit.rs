//! Cold positive-real matrix approximation of the ACTUAL modal BEM load.
//! A fixed damped-resonator dictionary and PSD residue matrices preserve all
//! signed cross-port terms. This is a restricted passive fit, not general MIMO
//! vector fitting. Inadequate data/basis refuses; no entrywise repair or gain EQ.
//!
//! [`fit`] keeps the shared 32-term dictionary. [`fit_resolved`] derives its
//! dictionary from the load's longest travel time instead: a radiator's kernel
//! is a superposition of delays up to tau = D/c, so its impedance has detail on a
//! scale of 1/tau in FREQUENCY at every frequency, which a fixed log-spaced set
//! cannot follow over a wide band. Same solver, same unchanged error limits.
use super::{Model,Pole,MAX_PORTS};
use fs_math::{c64::C64,det};
const TERMS:usize=32;
const SWEEPS:usize=512;
/// A resolved dictionary that needs more terms than this refuses: the band is
/// too wide for the radiator's travel time to be represented.
pub const MAX_RESOLVED_TERMS:usize=512;
/// Resolved resonators continue this many spacings past the band top, so the
/// fitted endpoint is not a truncation edge of the representation.
const GUARD_TERMS:usize=6;

pub struct Fit {
    pub model:Model,
    pub peak_error:f64,
    pub rms_error:f64,
    pub resistance_peak_error:f64,
    pub resistance_rms_error:f64,
}

/// The shared passive fit with frequency-local discrepancy limits. These
/// are sampled Frobenius-norm errors, not bounds between sampled frequencies
/// or on every multiport velocity direction.
pub struct RelativeFit {
    pub fit:Fit,
    pub complex_peak_relative_error:f64,
    pub resistance_peak_relative_error:f64,
}
fn dictionary(low:f64,high:f64)->Vec<(f64,f64)> {
    let mut poles=Vec::with_capacity(TERMS);
    for &factor in &[0.35,0.55,0.80] {
        poles.push((low*factor,0.25));
    }
    let in_band=22usize;
    for j in 0..in_band {
        let w=det::exp(det::ln(low)+(det::ln(high)-det::ln(low))*j as f64/(in_band-1) as f64);
        poles.push((w,0.12));
    }
    let high_band=TERMS-1-poles.len();
    for j in 1..=high_band {
        let w=det::exp(det::ln(high)+(det::ln(8.*high)-det::ln(high))*j as f64/high_band as f64);
        poles.push((w,0.25));
    }
    // Off-band lossless storage permits compact dipole loads without imposing
    // an artificial minimum resistance. Its pole must still resolve at runtime.
    poles.push((8.*high,0.));poles
}
/// In-band resonators every `spacing` rad/s from the band bottom, each with
/// half-power bandwidth 2*spacing (zeta = spacing/omega, at most 0.5), continued
/// GUARD_TERMS past the top; storage just above the band for the frequency-
/// dependent added mass; and the shared sub-band and high-band terms. The
/// storage terms carry the same travel-time bandwidth as the in-band
/// resonators: the kernel holds no delay longer than D/c, so no term may ring
/// longer than they do. Lossless storage did: above a board whose modes stop
/// at the band top, each such pole was an undamped mode of the whole system,
/// heard as steady lines at 1.15..8x the top in every air-loaded piano render.
fn resolved_dictionary(low:f64,high:f64,spacing:f64)->Vec<(f64,f64)> {
    let mut poles=Vec::new();
    for &factor in &[0.35,0.55,0.80] {poles.push((low*factor,0.25));}
    let in_band=((high-low)/spacing).floor() as usize;
    for j in 0..=in_band+GUARD_TERMS {
        let w=low+spacing*j as f64;poles.push((w,(spacing/w).min(0.5)));
    }
    for j in 1..=6 {
        poles.push((det::exp(det::ln(high)+(det::ln(8.*high)-det::ln(high))*f64::from(j)/6.),0.25));
    }
    for &factor in &[1.15,1.3,1.6,2.,3.,5.,8.] {let w=factor*high;poles.push((w,(spacing/w).min(0.5)));}
    poles
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
/// The same projection through Householder/QL, for the resolved dictionary's
/// hundred-odd terms per sweep; the shared fixed dictionary keeps [`psd`].
fn psd_householder(a:&[f64],n:usize)->Result<Vec<f64>,String> {
    if a.iter().any(|x|!x.is_finite()) {return Err("passive load fit overflow".into());}
    if n==1 {return Ok(vec![a[0].max(0.)]);}
    let scale=a.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
    if scale==0. {return Ok(vec![0.;n*n]);}
    let (values,vectors)=fs_la::eigen::householder_eigh(a,n);
    let mut out=vec![0.;n*n];
    for (k,&l) in values.iter().enumerate() {
        let mut residual=0.0_f64;
        for i in 0..n {
            let mut s=0.;for j in 0..n {s+=a[i*n+j]*vectors[j*n+k];}
            residual=residual.max((s-l*vectors[i*n+k]).abs());
        }
        if !l.is_finite() || !(residual<=1e-8*scale*n as f64) {
            return Err("passive fit PSD projection did not resolve its eigenproblem".into());
        }
        if l<=0. {continue;}
        for i in 0..n {for j in 0..n {out[i*n+j]+=l*vectors[i*n+k]*vectors[j*n+k];}}
    }
    Ok(out)
}

fn quadratic_objective(residue:&[Vec<f64>],gram:&[f64],rhs:&[Vec<f64>])->f64 {
    let terms=residue.len();
    let mut total=0.0_f64;
    for j in 0..terms {
        let mut dot_rhs=0.0_f64;
        for (&r,&b) in residue[j].iter().zip(&rhs[j]) {dot_rhs+=r*b;}
        let mut quad=0.0_f64;
        for k in 0..terms {
            let g=gram[j*terms+k];
            let mut inner=0.0_f64;
            for (&rj,&rk) in residue[j].iter().zip(&residue[k]) {inner+=rj*rk;}
            quad+=g*inner;
        }
        total+=0.5*quad-dot_rhs;
    }
    total
}

fn block_psd_sweep(
    start:&[Vec<f64>],
    gram:&[f64],
    rhs:&[Vec<f64>],
    ports:usize,
    id:&[f64],
    householder:bool,
)->Result<(Vec<Vec<f64>>,f64,f64),String> {
    let size=ports*ports;let terms=start.len();
    let mut state=start.to_vec();
    let mut change=0.0_f64;
    let mut magnitude=0.0_f64;
    for j in 0..terms {
        let mut candidate=vec![0.;size];
        for (a,value) in candidate.iter_mut().enumerate() {
            *value=(rhs[j][a]-(0..terms).filter(|k|*k!=j)
                .map(|k|gram[j*terms+k]*state[k][a]).sum::<f64>())/gram[j*terms+j];
        }
        for i in 0..ports {for k in 0..i {
            let v=0.5*(candidate[i*ports+k]+candidate[k*ports+i]);
            candidate[i*ports+k]=v;candidate[k*ports+i]=v;
        }}
        let projected=if householder {psd_householder(&candidate,ports)?} else {psd(&candidate,ports,id)?};
        for (old,new) in state[j].iter_mut().zip(projected) {
            change=change.max((*old-new).abs());magnitude=magnitude.max(new.abs());*old=new;
        }
    }
    Ok((state,change,magnitude))
}

fn admit(omega:&[f64],data:&[Vec<C64>],ports:usize)->Result<(),String> {
    if !(1..=MAX_PORTS).contains(&ports) || !(33..=257).contains(&omega.len()) || omega.len()%2==0
        || data.len()!=omega.len() || data.iter().any(|r|r.len()!=ports*ports
            || r.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()))
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
            || (i>0 && *w<=omega[i-1])) {
        return Err(format!("passive load fit requires a complete 1..{MAX_PORTS}-port matrix on an odd 33..257-frequency grid"));
    }
    Ok(())
}

/// Fit EVEN samples only. Odd samples are never used for coefficients,
/// convergence stopping, dictionary choice or tolerance selection.
pub fn fit(omega:&[f64],data:&[Vec<C64>],ports:usize)->Result<Fit,String> {
    admit(omega,data,ports)?;
    fit_dictionary(omega,data,ports,&dictionary(omega[0],omega[omega.len()-1]),false,false)
}

/// Fit with a dictionary that resolves the load's longest travel time. `spacing`
/// is the in-band resonator spacing in rad/s; a radiator of largest extent D in a
/// medium of sound speed c needs spacing <= 2 pi c / (6 D) (a sixth of 1/tau). The
/// training samples (even indices) must be no further apart than `spacing`, or
/// the fit refuses rather than interpolate detail it never saw.
pub fn fit_resolved(omega:&[f64],data:&[Vec<C64>],ports:usize,spacing:f64)->Result<Fit,String> {
    let poles=admitted_resolved_dictionary(omega,data,ports,spacing)?;
    fit_dictionary(omega,data,ports,&poles,true,false)
}

/// Resolve weak low-frequency radiation resistance as well as the larger
/// reactive response. Training remains EVEN samples only. The same PSD block
/// solver minimizes frequency-local complex and real-part errors; column
/// scaling changes numerical conditioning, not the objective or residue cone.
/// Eight logarithmic resonators fill the first linear dictionary interval,
/// where a travel-time spacing alone may not resolve the small resistance.
///
/// Both training and held-out samples must have at most 2% relative complex
/// and 2% relative Hermitian-part error at EVERY supplied frequency, in addition
/// to the original global peak/RMS gates. A zero local reference or more than
/// twelve decades of required weighting refuses. No recorded datum, frequency
/// band or tolerance is modified to make a fit pass; no new solver is selected.
pub fn fit_resolved_relative(omega:&[f64],data:&[Vec<C64>],ports:usize,spacing:f64)
    ->Result<RelativeFit,String>
{
    let mut poles=admitted_resolved_dictionary(omega,data,ports,spacing)?;
    let low=omega[0];let first=low+spacing;
    let extra=(1..=8).map(|i| (det::exp(det::ln(low)
        +(det::ln(first)-det::ln(low))*f64::from(i)/9.),0.5));
    poles.splice(3..3,extra);
    if poles.len()>MAX_RESOLVED_TERMS {
        return Err("relative load fit exceeds the resolved dictionary ceiling".into());
    }
    let fit=fit_dictionary(omega,data,ports,&poles,true,true)?;
    let mut complex=0.0_f64;let mut resistance=0.0_f64;
    for (&w,reference) in omega.iter().zip(data) {
        let actual=fit.model.impedance(w)?;
        let error:Vec<_>=actual.iter().zip(reference).map(|(a,b)|*a-*b).collect();
        let mut h=0.0_f64;let mut dh=0.0_f64;
        for i in 0..ports {for j in 0..ports {
            h=h.hypot((reference[i*ports+j]+reference[j*ports+i].conj()).scale(0.5).abs());
            dh=dh.hypot((error[i*ports+j]+error[j*ports+i].conj()).scale(0.5).abs());
        }}
        if h==0. || norm(reference)==0. {
            return Err("relative load fit needs nonzero complex and resistive references at every frequency".into());
        }
        complex=complex.max(norm(&error)/norm(reference));
        resistance=resistance.max(dh/h);
    }
    if !complex.is_finite() || !resistance.is_finite() || complex>0.02 || resistance>0.02 {
        return Err(format!("relative passive load fit refused: pointwise complex/resistance={complex:.5}/{resistance:.5}; limits .02/.02"));
    }
    Ok(RelativeFit {fit,complex_peak_relative_error:complex,resistance_peak_relative_error:resistance})
}

fn admitted_resolved_dictionary(omega:&[f64],data:&[Vec<C64>],ports:usize,spacing:f64)
    ->Result<Vec<(f64,f64)>,String>
{
    admit(omega,data,ports)?;
    // At most the band itself, so the guard terms stay below the 8x off-band storage.
    if !spacing.is_finite() || spacing<=0. || spacing>omega[omega.len()-1]-omega[0] {
        return Err("resolved load fit needs a finite positive resonator spacing no wider than its band".into());
    }
    let widest=(2..omega.len()).step_by(2).map(|f|omega[f]-omega[f-2]).fold(0.0_f64,f64::max);
    if widest>spacing {
        return Err(format!("resolved load fit: training samples are up to {:.1} Hz apart, wider than the {:.1} Hz the load's travel time needs; sample the band more densely",
            widest/std::f64::consts::TAU,spacing/std::f64::consts::TAU));
    }
    let poles=resolved_dictionary(omega[0],omega[omega.len()-1],spacing);
    if poles.len()>MAX_RESOLVED_TERMS {
        return Err(format!("resolved load fit needs {} terms, above its {MAX_RESOLVED_TERMS}-term ceiling",poles.len()));
    }
    Ok(poles)
}

fn fit_dictionary(omega:&[f64],data:&[Vec<C64>],ports:usize,poles:&[(f64,f64)],householder:bool,
    relative:bool)
    ->Result<Fit,String> {
    let size=ports*ports;let terms=poles.len();
    let scale=data.iter().step_by(2).map(|r|norm(r)).fold(0.0_f64,f64::max);
    if !scale.is_finite() {return Err("acoustic load norm overflow".into());}
    if scale==0. {
        if data.iter().flatten().any(|z|z.abs()!=0.) {return Err("zero training load has a nonzero held-out response".into());}
        return Ok(Fit {model:Model {ports,poles:vec![Pole {omega:poles[0].0,zeta:poles[0].1,coupling:vec![0.;ports]}]},
            peak_error:0.,rms_error:0.,resistance_peak_error:0.,resistance_rms_error:0.});
    }
    // Balance resistive and reactive error using TRAINING data only. A very
    // small physical radiation resistance must not disappear in an inertance-
    // dominated complex least-squares objective, while respecting that the
    // complex RMS limit (.02) is tighter than the resistance RMS limit (.05).
    // Blocks still have a scalar Hessian, so the PSD projection remains their
    // exact minimizer.
    let real_scale=data.iter().step_by(2).map(|row|row.iter().fold(0.0_f64,|s,z|s.hypot(z.re)))
        .fold(0.0_f64,f64::max);
    let real_weight=if real_scale==0. {1.} else {(0.35*(scale/real_scale)).max(1.0)};
    if !real_weight.is_finite() || real_weight>1e12 {
        return Err("radiation resistance is below the load fit's relative resolution budget".into());
    }
    let mut gram=vec![0.;terms*terms];let mut rhs=vec![vec![0.;size];terms];
    for f in (0..omega.len()).step_by(2) {
        let (real_weight,imag_weight)=if relative {
            let real=data[f].iter().fold(0.0_f64,|s,z|s.hypot(z.re));
            let imaginary=norm(&data[f]);
            let weights=(scale/real,scale/imaginary);
            if !weights.0.is_finite() || !weights.1.is_finite()
                || weights.0>1e12 || weights.1>1e12 {
                return Err("relative radiation response exceeds the fit's twelve-decade resolution budget".into());
            }
            weights
        } else {(real_weight,1.)};
        let h:Vec<_>=poles.iter().map(|&p|kernel(omega[f],p)).collect();
        for j in 0..terms {
            for k in 0..terms {
                gram[j*terms+k]+=h[j].re*real_weight*h[k].re*real_weight
                    +h[j].im*imag_weight*h[k].im*imag_weight;
            }
            for (a,value) in rhs[j].iter_mut().enumerate() {
                let y=data[f][a].scale(1./scale);
                *value+=h[j].re*real_weight*y.re*real_weight+h[j].im*imag_weight*y.im*imag_weight;
            }
        }
    }
    if gram.iter().any(|x|!x.is_finite()) || (0..terms).any(|j|gram[j*terms+j]<=0.) {
        return Err("unresolved passive fit dictionary".into());
    }
    let columns:Vec<_>=(0..terms).map(|j|if relative {det::sqrt(gram[j*terms+j])} else {1.}).collect();
    if relative {
        for j in 0..terms {
            for k in 0..terms {gram[j*terms+k]=gram[j*terms+k]/columns[j]/columns[k];}
            for b in &mut rhs[j] {*b/=columns[j];}
        }
    }
    let id=identity(ports);let mut residue=vec![vec![0.;size];terms];
    let mut y=residue.clone();
    let mut obj=quadratic_objective(&residue,&gram,&rhs);
    let mut t=1.0_f64;
    // Exact block minimization: for fixed other blocks, the quadratic Hessian
    // is a positive scalar times identity. Symmetric PSD projection solves
    // this block's constrained least-squares problem, using fs-modal's owner.
    // Monotone Nesterov extrapolation accelerates ill-conditioned dictionary
    // coupling while keeping every accepted iterate strictly in the PSD cone.
    for _ in 0..SWEEPS {
        let (mut next,mut change,mut magnitude)=block_psd_sweep(&y,&gram,&rhs,ports,&id,householder)?;
        let mut next_obj=quadratic_objective(&next,&gram,&rhs);
        if !(next_obj<=obj) {
            (next,change,magnitude)=block_psd_sweep(&residue,&gram,&rhs,ports,&id,householder)?;
            next_obj=quadratic_objective(&next,&gram,&rhs);
            t=1.0;
            y.clone_from(&next);
        } else {
            let next_t=0.5*(1.0+det::sqrt(1.0+4.0*t*t));
            let beta=(t-1.0)/next_t;
            for j in 0..terms {
                for a in 0..size {
                    y[j][a]=next[j][a]+beta*(next[j][a]-residue[j][a]);
                }
            }
            t=next_t;
        }
        residue=next;
        obj=next_obj;
        if change<=1e-10*(1.+magnitude) {break;}
    }
    let mut model=Model {ports,poles:Vec::new()};
    for (index,(matrix,&(o,z))) in residue.iter().zip(poles).enumerate() {
        let pairs=fs_modal::eigh_gen_dense(matrix,&id,ports).map_err(|e|e.to_string())?;
        let tolerance=1e-8*ports as f64*matrix.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
        for p in pairs {
            if !p.lambda.is_finite() || !p.residual.is_finite() || p.residual>tolerance
                || p.lambda < -tolerance || p.phi.len()!=ports || p.phi.iter().any(|x|!x.is_finite()) {
                return Err("passive residue factorization failed".into());
            }
            if p.lambda<=0. {continue;}
            let gain=det::sqrt(p.lambda/columns[index]*scale*normalization((o,z)));
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

#[cfg(test)]
mod tests {
    use super::*;
    fn grid()->Vec<f64> {(0..65).map(|i|std::f64::consts::TAU*(40.+360.*i as f64/64.)).collect()}
    #[test]
    fn fit_and_model_share_the_complete_port_capacity() {
        let omega = grid();
        for ports in [0, MAX_PORTS, MAX_PORTS + 1] {
            let data = vec![vec![C64::ZERO; ports * ports]; omega.len()];
            let result = fit(&omega, &data, ports);
            let model = Model { ports, poles: vec![Pole {
                omega: omega[0], zeta: 0.2, coupling: vec![0.0; ports],
            }] };
            assert_eq!(model.validate().is_ok(), ports == MAX_PORTS,
                "model capacity at {ports} complete ports");
            if ports == MAX_PORTS {
                result.unwrap().model.validate().unwrap();
            } else {
                assert!(result.is_err(), "accepted {ports} complete ports");
            }
        }
    }
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
    fn resolved_dictionary_follows_a_rayleigh_load_the_fixed_one_cannot() {
        use fs_bem::helmholtz::{Medium, rayleigh_modal::{AreaSample, PreparedSurface}};
        use std::f64::consts::{PI, TAU};
        // A baffled 2.0 m x 1.5 m plate in 5 cm cells carrying 25 fields
        // cos(m pi x / L) cos(n pi y / W): its kernel holds delays up to about 2.4 m / c,
        // the scale of a grand piano's soundboard (whose own 128 fields need the same).
        let (nx, ny, h, side) = (40usize, 30usize, 0.05, 5usize);
        let cells: Vec<_> = (0..nx * ny).map(|k| {
            let (x, y) = (((k % nx) as f64 + 0.5) * h, ((k / nx) as f64 + 0.5) * h);
            AreaSample { position_m: [x, y], area_m2: h * h, mode_shape: (0..side * side)
                .map(|f| (PI * (f / side) as f64 * x / 2.0).cos() * (PI * (f % side) as f64 * y / 1.5).cos()).collect() }
        }).collect();
        let ports = side * side;
        let surface = PreparedSurface::new(&cells).unwrap();
        let medium = Medium { density: 1.2, sound_speed: 343.0 };
        let omega: Vec<f64> = (0..129).map(|i| TAU * (90.0 + 900.0 * i as f64 / 128.0)).collect();
        let data: Vec<Vec<C64>> = omega.iter()
            .map(|&w| surface.power_consistent_impedance(w, medium).unwrap()).collect();
        assert!(fit(&omega, &data, ports).is_err(), "the fixed 32-term dictionary should not resolve a 2.4 m travel time");
        let spacing = TAU * medium.sound_speed / (6.0 * surface.maximum_distance_m());
        let resolved = fit_resolved(&omega, &data, ports, spacing).unwrap();
        assert!(resolved.rms_error <= 0.02 && resolved.resistance_rms_error <= 0.05);
        // No term at or above the band outlasts the travel time: a lossless off-band pole rang for ever.
        assert!(resolved.model.poles.iter().filter(|p| p.omega >= omega[0])
            .all(|p| p.zeta * p.omega >= spacing * (1.0 - 1e-12)));
        // Too few training samples for that spacing refuses instead of interpolating.
        let sparse: Vec<f64> = (0..33).map(|i| TAU * (90.0 + 900.0 * i as f64 / 32.0)).collect();
        let sparse_data: Vec<Vec<C64>> = sparse.iter()
            .map(|&w| surface.power_consistent_impedance(w, medium).unwrap()).collect();
        assert!(fit_resolved(&sparse, &sparse_data, ports, spacing).is_err());
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
