//! Geometry sampling uses the shared passive multiport fitter, not a second fit.
use super::{engine::radiation::PortBasis,exterior_geometry::{Boundary,Specification,Samples},exterior_loading};
use fs_math::c64::C64;
pub use fs_couple::render::plate::impact::radiation::fit::{Fit,fit};
const MAX_PORTS:usize=32;

/// Solve one modal BEM batch per frequency. Those SAME fields feed the load
/// and receiver fits. Pressure/acceleration = pressure/velocity * i/omega.
pub fn prepare(boundary:&Boundary,spec:&Specification,mechanics_rate:u32)->Result<(Fit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if r>MAX_PORTS || omega.len()<33 {return Err("loaded playback requires <=32 complete board modes and >=33 frequencies".into());}
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let prepared=exterior_loading::LoadingSweep::new(boundary,spec,&omega)?;
    let delays_s=prepared.delays_s().to_vec();
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=prepared.sample(w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fitted=fit(&omega,&matrices,r)?;
    Ok((fitted,Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}

/// The fitted acoustic coordinates are independent of the complete structural
/// basis. Error fields in `fit` refer to the lifted COMPLETE BEM load, including
/// projection error, not just the retained acoustic matrix.
pub struct ProjectedFit {
    pub fit:Fit,
    pub basis:PortBasis,
    /// Relative area-weighted normal-velocity projection error. This geometric
    /// diagnostic does not replace the full complex and power acceptance gates.
    pub surface_rms_error:f64,
}

/// Select radiating combinations from geometry alone, before seeing any BEM
/// frequencies or receiver positions. Rows of Q^T are Euclidean-orthonormal in
/// the bank's mass-normalized board coordinates. The area-weighted Gram picks
/// the dominant normal-velocity fields, not a new structural eigenspace.
fn geometry_basis(weights:&[Vec<f64>],areas:&[f64],maximum:usize)
    ->Result<(PortBasis,f64),String> {
    let n=weights.len();let panels=areas.len();
    if !(1..=super::linear::MAX_BOARD_MODES).contains(&n)
        || !(1..=MAX_PORTS).contains(&maximum)
        || !(1..=super::exterior_geometry::MAX_PANELS).contains(&panels)
        || areas.iter().any(|a|!a.is_finite() || *a<=0.)
        || weights.iter().any(|row|row.len()!=panels || row.iter().any(|v|!v.is_finite())) {
        return Err("radiation port selection needs complete finite board fields, positive panel areas and 1..32 ports".into());
    }
    let mut weighted=vec![vec![0.;panels];n];let mut scale=0.0_f64;
    for (out,row) in weighted.iter_mut().zip(weights) {
        for ((value,&shape),&area) in out.iter_mut().zip(row).zip(areas) {
            *value=shape*area.sqrt();
            if !value.is_finite() || (shape!=0. && *value==0.) {
                return Err("radiation surface weighting exceeds finite numerical range".into());
            }
            scale=scale.max(value.abs());
        }
    }
    if scale==0. {
        // A physically stationary exterior has exactly zero load. Keep one
        // named coordinate for the shared fitter's existing zero-load image.
        let mut row=vec![0.;n];row[0]=1.;
        return Ok((PortBasis::new(n,vec![row])?,0.));
    }
    for value in weighted.iter_mut().flatten() {
        let old=*value;*value/=scale;
        if old!=0. && *value==0. {
            return Err("radiation surface scale hides a nonzero field".into());
        }
    }
    let mut gram=vec![0.;n*n];let mut identity=vec![0.;n*n];
    for i in 0..n {
        identity[i*n+i]=1.;
        for j in 0..=i {
            let value=weighted[i].iter().zip(&weighted[j]).map(|(a,b)|a*b).sum();
            gram[i*n+j]=value;gram[j*n+i]=value;
        }
    }
    let trace=(0..n).map(|i|gram[i*n+i]).sum::<f64>();
    let mut pairs=fs_modal::eigh_gen_dense(&gram,&identity,n).map_err(|e|e.to_string())?;
    let tolerance=256.*f64::EPSILON*n as f64*trace;
    if !trace.is_finite() || trace<=0. || pairs.len()!=n
        || pairs.iter().any(|p|!p.lambda.is_finite() || p.lambda< -tolerance
            || !p.residual.is_finite() || p.residual>tolerance
            || p.phi.len()!=n || p.phi.iter().any(|v|!v.is_finite())) {
        return Err("radiation surface Gram eigenproblem did not resolve its positive semidefinite field".into());
    }
    // Stable ordering retains the eigensolver's order for equal eigenvalues.
    // No pressure sample or fit failure changes this geometry-only selection.
    pairs.sort_by(|a,b|b.lambda.total_cmp(&a.lambda));
    let mut vectors=Vec::new();
    for pair in pairs.into_iter().filter(|p|p.lambda>tolerance).take(maximum) {
        let mut row=pair.phi;let mut pivot=0;
        for i in 1..n {if row[i].abs()>row[pivot].abs() {pivot=i;}}
        if row[pivot]<0. {for value in &mut row {*value= -*value;}}
        vectors.push(row);
    }
    if vectors.is_empty() {return Err("radiation surface has no resolved nonzero port".into());}
    let basis=PortBasis::new(n,vectors)?;
    let mut error=0.0_f64;let mut reference=0.0_f64;
    let mut projected=vec![0.;n];
    for panel in 0..panels {
        projected.fill(0.);
        for row in basis.vectors() {
            let amplitude=row.iter().zip(&weighted).map(|(q,w)|q*w[panel]).sum::<f64>();
            for (value,q) in projected.iter_mut().zip(row) {*value+=q*amplitude;}
        }
        for (value,row) in projected.iter().zip(&weighted) {
            reference=reference.hypot(row[panel]);error=error.hypot(value-row[panel]);
        }
    }
    let surface_rms_error=error/reference;
    if !surface_rms_error.is_finite() {return Err("nonfinite radiation surface projection error".into());}
    Ok((basis,surface_rms_error))
}

/// Apply the same whole-matrix tolerances as the existing shared fit to the
/// lifted response. Both even training and odd held-out samples are checked;
/// no held-out value selects the geometry basis or any fit coefficient.
fn validate_complete_fit(omega:&[f64],data:&[Vec<C64>],basis:&PortBasis,
    fitted:&mut Fit)->Result<(),String> {
    let ports=basis.board_ports();
    if data.len()!=omega.len() || data.iter().any(|row|row.len()!=ports*ports
        || row.iter().any(|z|!z.re.is_finite() || !z.im.is_finite())) {
        return Err("projected radiation validation requires the complete original BEM matrices".into());
    }
    let mut peak=0.0_f64;let mut rms=0.0_f64;let mut rp=0.0_f64;let mut rr=0.0_f64;
    for offset in 0..2 {
        let mut max_reference=0.0_f64;let mut max_error=0.0_f64;
        let mut references=0.0_f64;let mut errors=0.0_f64;
        let mut real_max=0.0_f64;let mut real_error=0.0_f64;
        let mut real_references=0.0_f64;let mut real_errors=0.0_f64;
        for f in (offset..omega.len()).step_by(2) {
            let lifted=basis.lift_impedance(&fitted.model.impedance(omega[f])?)?;
            let error:Vec<_>=lifted.iter().zip(&data[f]).map(|(a,b)|*a-*b).collect();
            let reference=data[f].iter().fold(0.0_f64,|s,z|s.hypot(z.abs()));
            let delta=error.iter().fold(0.0_f64,|s,z|s.hypot(z.abs()));
            max_reference=max_reference.max(reference);max_error=max_error.max(delta);
            references=references.hypot(reference);errors=errors.hypot(delta);
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
        return Err(format!("projected BEM load fit refused on the complete board: complex peak/RMS={peak:.5}/{rms:.5}, power peak/RMS={rp:.5}/{rr:.5}; limits .05/.02/.10/.05; supply more radiation ports or retain the complete load"));
    }
    fitted.peak_error=peak;fitted.rms_error=rms;
    fitted.resistance_peak_error=rp;fitted.resistance_rms_error=rr;
    Ok(())
}

fn fit_projected(omega:&[f64],data:&[Vec<C64>],basis:&PortBasis)->Result<Fit,String> {
    let projected=data.iter().map(|row|basis.project_impedance(row)).collect::<Result<Vec<_>,_>>()?;
    let mut fitted=fit(omega,&projected,basis.ports())?;
    validate_complete_fit(omega,data,basis,&mut fitted)?;
    Ok(fitted)
}

/// Keep every structural coordinate and receiver input, fitting only the
/// explicitly bounded geometry-derived radiating subspace. Each frequency is
/// solved once by the existing complete modal BEM batch. An insufficient rank
/// refuses after full-load validation; it never silently drops a mechanical
/// mode, relaxes the resistance gate, or substitutes an output-only load.
pub fn prepare_projected(boundary:&Boundary,spec:&Specification,mechanics_rate:u32,
    max_ports:usize)->Result<(ProjectedFit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if !(33..=257).contains(&omega.len()) || omega.len()%2==0
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
            || (i>0 && *w<=omega[i-1])) {
        return Err("projected loaded playback requires an odd 33..257-frequency increasing positive grid".into());
    }
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let (basis,surface_rms_error)=geometry_basis(&boundary.weights,boundary.surface.areas(),max_ports)?;
    let prepared=exterior_loading::LoadingSweep::new(boundary,spec,&omega)?;
    let delays_s=prepared.delays_s().to_vec();
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=prepared.sample(w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fit=fit_projected(&omega,&matrices,&basis)?;
    Ok((ProjectedFit {fit,basis,surface_rms_error},
        Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}

#[cfg(test)]
#[path="radiation_port_fit_tests.rs"]
mod port_tests;
