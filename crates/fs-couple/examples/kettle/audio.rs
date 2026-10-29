//! Closed-boundary BEM -> reciprocal load + fixed-point pressure -> shared PCM.
//! No compact monopole, material-name gain, new integrator or WAV implementation.
use super::{boundary::Boundary, model::Model, spec::Spec, Error};
use fs_bem::{helmholtz::{Medium, Formulation, PreparedCbieGeometry, solve_radiation_batch,
    exterior_pressure_at_points_integrated}, panel3d::SpherePanels};
use fs_couple::{pcm_wav::{decimate::Decimator, encode_pcm16_wav_interleaved},
    render::plate::impact::{audio::ImpactSource, radiation::{self, fit}}};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_plate::{ModePair, shell::head::TensionedDisk};
use fs_vfit::{vf::{FitOptions, WeightPreset, vector_fit},
    discretize::{DiscreteStateSpace, DelayedFilter, DigitalFilter, bilinear_state_space}};
use core::f64::consts::{PI, TAU};

/// Same mechanical source with variational work-conjugate exterior acoustics.
#[path="reciprocal_air.rs"]
pub mod reciprocal;

pub const OUTPUT_RATE: u32 = 48_000;
#[derive(Clone, Copy)]
pub struct Settings {
    pub full_scale_pa: f64,
    pub band_hz: [f64;2],
    /// An even/odd train/audit lattice with 2*intervals+1 frequencies.
    pub intervals: usize,
    pub order: usize,
    pub receiver_m: [f64;3],
}
impl Settings {
    pub fn read(args: &[String], spec: &Spec) -> Result<Self, Error> {
        if args.len()!=9 || args[0]!="--wav" {
            return Err("audio syntax: --wav FULL_SCALE_PA MIN_HZ MAX_HZ INTERVALS ORDER X_M Y_M Z_M".into());
        }
        let s=Self { full_scale_pa:args[1].parse()?, band_hz:[args[2].parse()?,args[3].parse()?],
            intervals:args[4].parse()?,order:args[5].parse()?,
            receiver_m:[args[6].parse()?,args[7].parse()?,args[8].parse()?] };
        s.validate(spec)?; Ok(s)
    }
    pub fn validate(self, spec: &Spec) -> Result<(), Error> {
        spec.validate_exterior()?;
        if !self.full_scale_pa.is_finite() || self.full_scale_pa<=0.0
            || self.band_hz.iter().chain(&self.receiver_m).any(|v|!v.is_finite())
            || self.band_hz[0]<=0.0 || self.band_hz[1]<=self.band_hz[0]
            || self.band_hz[1]>=0.45*f64::from(OUTPUT_RATE)
            || !(16..=128).contains(&self.intervals) || !(2..=32).contains(&self.order)
            || self.intervals<2*self.order || self.order%2!=0
            || spec.rate_hz%OUTPUT_RATE!=0 || !(2..=32).contains(&(spec.rate_hz/OUTPUT_RATE))
            // The shared passive dictionary extends to eight times the band.
            || 8.0*TAU*self.band_hz[1]/f64::from(spec.rate_hz)>=0.9*PI
        { return Err("audio needs finite physical controls, an overdetermined 33..257-point grid, resolved load poles and integer 2..32x mechanical oversampling".into()); }
        Ok(())
    }
}
pub struct Bake { filters:Vec<DiscreteStateSpace>, delay_s:f64 }
fn phase(value:C64,angle:f64)->C64 { value*C64::new(det::cos(angle),det::sin(angle)) }

/// Fixed declared order: odd samples never choose normalization or poles.
fn observer(omega:&[f64], values:&[C64], order:usize)->Result<DiscreteStateSpace,Error> {
    let dt=1.0/f64::from(OUTPUT_RATE);
    if omega.len()!=values.len() || omega.len()<5 || omega.len()%2!=1
        || !(2..=32).contains(&order) || (omega.len()+1)/2<2*order+1
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.0 || *w*dt>=PI
            || i>0 && *w<=omega[i-1])
        || values.iter().any(|v|!v.re.is_finite() || !v.im.is_finite()) {
        return Err("invalid fixed pressure-fit basis".into());
    }
    let scale=values.iter().step_by(2).map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale==0.0 {
        if values.iter().any(|v|v.re!=0.0 || v.im!=0.0) { return Err("zero training pressure has nonzero audit pressure".into()); }
        return Ok(DiscreteStateSpace {n:0,a:vec![],b:vec![],c:vec![],d:0.0,e_leftover:0.0,t_s:dt});
    }
    if !scale.is_finite() { return Err("nonfinite pressure normalization".into()); }
    let anchor=values.iter().step_by(2).max_by(|a,b|a.abs().total_cmp(&b.abs())).ok_or("empty pressure training")?;
    let scale=scale.copysign(if anchor.re.abs()>=anchor.im.abs(){anchor.re}else{anchor.im});
    // Acoustic owner uses exp(-i*w*t), fit owner evaluates at s=+i*w.
    // Prewarp TRAINING frequencies for the intended output-clock realization.
    let frequencies:Vec<_>=omega.iter().step_by(2).map(|w|2.0/dt*det::tan(*w*dt/2.0)).collect();
    let response:Vec<_>=values.iter().step_by(2).map(|v|C64::new(v.re/scale,-v.im/scale)).collect();
    let result=vector_fit(&frequencies,&response,&FitOptions {
        order,iterations:12,weights:WeightPreset::Uniform,fit_d:true,fit_e:false })?;
    if !result.model.is_stable() || result.model.e!=0.0 || !result.report.weighted_rms.is_finite() { return Err("pressure fit is not stable and proper".into()); }
    let mut filter=bilinear_state_space(&result.model,dt,0.0)?;
    for c in &mut filter.c { *c*=scale; } filter.d*=scale;
    filter.try_runtime()?;
    let mut peak=0.0_f64; let mut square=0.0; let mut count=0;
    for i in (1..omega.len()).step_by(2) {
        let error=(filter.eval(omega[i])?.conj()-values[i]).abs()/scale.abs();
        if !error.is_finite() { return Err("nonfinite independent pressure audit".into()); }
        peak=peak.max(error); square+=error*error; count+=1;
    }
    let rms=(square/count as f64).sqrt();
    if peak>0.15 || rms>0.05 { return Err(format!("independent pressure audit refused: peak={peak}, rms={rms}; fixed order {order}").into()); }
    eprintln!("pressure transfer: order={order}, independent_audit_peak={peak}, independent_audit_rms={rms}");
    Ok(filter)
}
/// Test the whole Hermitian power form, not only each diagonal's power.
fn admit_power(z:&[C64], n:usize)->Result<(),Error> {
    if n==0 || n>32 || z.len()!=n*n || z.iter().any(|v|!v.re.is_finite() || !v.im.is_finite()) {
        return Err("incomplete or nonfinite exterior impedance".into());
    }
    let m=2*n; let mut h=vec![0.0;m*m];
    for i in 0..n { for j in 0..n {
        let v=(z[i*n+j]+z[j*n+i].conj()).scale(0.5);
        h[i*m+j]=v.re; h[(i+n)*m+j+n]=v.re;
        h[i*m+j+n]=-v.im; h[(i+n)*m+j]=v.im;
    }}
    fs_phs::QuadraticStorage::new(h,m).map_err(|e|format!("exterior multiport power refused: {e}"))?;
    Ok(())
}
/// Reciprocity makes the exact modal radiation impedance complex symmetric
/// (Z = Z^T for real mode shapes); collocation breaks it by its discretisation
/// error. Keep the reciprocal part (Z+Z^T)/2 in place and return the stated
/// asymmetry |Z-Z^T|/|Z| (Frobenius), refused above the load fit's own 0.02
/// complex RMS limit, with the lowest eigenvalue of the power form Re Z_r over
/// |Z|. That eigenvalue must be indistinguishable from zero at the stated
/// error, at least -asymmetry*|Z|; nothing is clipped.
fn reciprocal_part(z:&mut [C64], n:usize)->Result<(f64,f64),String> {
    if n==0 || n>32 || z.len()!=n*n || z.iter().any(|v|!v.re.is_finite() || !v.im.is_finite()) {
        return Err("incomplete or nonfinite exterior impedance".into());
    }
    let (mut norm,mut skew)=(0.0_f64,0.0_f64);
    for i in 0..n { for j in 0..n {
        norm=norm.hypot(z[i*n+j].abs()); skew=skew.hypot((z[i*n+j]-z[j*n+i]).abs());
    }}
    if norm==0.0 { return Ok((0.0,0.0)); }
    let asymmetry=skew/norm;
    if !(asymmetry<=0.02) {
        return Err(format!("exterior impedance is not reciprocal: |Z-Z^T|/|Z|={asymmetry:.3e} exceeds 0.02"));
    }
    for i in 0..n { for j in 0..i {
        let v=(z[i*n+j]+z[j*n+i]).scale(0.5); z[i*n+j]=v; z[j*n+i]=v;
    }}
    let power:Vec<f64>=z.iter().map(|v|v.re).collect();
    let lowest=fs_la::eigen::householder_eigh(&power,n).0.into_iter().fold(f64::INFINITY,f64::min)/norm;
    if !lowest.is_finite() || lowest < -asymmetry {
        return Err(format!("exterior multiport power refused: lowest eigenvalue of Re((Z+Z^T)/2) is {lowest:.3e}|Z|, below -{asymmetry:.3e}|Z|"));
    }
    Ok((asymmetry,lowest))
}
/// One set of BEM source solutions supplies BOTH acoustic work and observation.
pub fn bake(spec:&Spec,film:&TensionedDisk,modes:&[ModePair],settings:Settings,gate:&CancelGate)
    ->Result<(radiation::Model,Bake),Error> {
    settings.validate(spec)?;
    if gate.is_requested() { return Err("exterior preparation cancelled".into()); }
    if modes.iter().any(|m|m.lambda.sqrt()/TAU<settings.band_hz[0] || m.lambda.sqrt()/TAU>settings.band_hz[1]) {
        return Err("the acoustic band must cover the complete retained head basis".into());
    }
    let boundary=Boundary::new(spec,film,modes)?;
    let n=modes.len(); let panels=boundary.triangles.len(); let samples=2*settings.intervals+1;
    let work=(panels as u64).checked_pow(3).and_then(|w|w.checked_mul(samples as u64))
        .ok_or("exterior work budget overflow")?;
    if work>100_000_000_000 { return Err("exterior sweep exceeds 100 billion dense matrix work units".into()); }
    let medium=Medium {density:spec.air[0],sound_speed:spec.air[1]};
    let range=settings.receiver_m[0].hypot(settings.receiver_m[1]).hypot(settings.receiver_m[2]);
    let delay_s=(range-boundary.radius_m)/medium.sound_speed;
    if !delay_s.is_finite() || delay_s<2.0/f64::from(OUTPUT_RATE) || delay_s>2.0 {
        return Err("finite receiver must be outside the enclosing sphere, with 2 samples..2 seconds of explicit propagation".into());
    }
    let omega:Vec<_>=(0..samples).map(|i|TAU*(settings.band_hz[0]
        +(settings.band_hz[1]-settings.band_hz[0])*i as f64/(samples-1) as f64)).collect();
    let k_max=omega[samples-1]/medium.sound_speed;
    // Slab Dirichlet bound excludes CBIE irregular frequencies across this band.
    let formulation=if (k_max*boundary.minimum_box_width_m).next_up()<PI.next_down() {
        Formulation::PlainCbie
    } else { Formulation::BurtonMiller };
    // Largest extent D of the closed surface: its load has detail every c/D in frequency.
    let extent=boundary.triangles.iter().flatten().flat_map(|a|boundary.triangles.iter().flatten()
        .map(move |b|(a[0]-b[0]).hypot(a[1]-b[1]).hypot(a[2]-b[2]))).fold(0.0_f64,f64::max);
    let surface=SpherePanels::from_triangles(boundary.triangles)?;
    let geometry=if formulation==Formulation::PlainCbie {
        Some(PreparedCbieGeometry::new_with_cancel(&surface,k_max,
            2048*2048*core::mem::size_of::<(f64,f64)>(),||gate.is_requested())?)
    } else { None };
    let mut values=vec![vec![C64::ZERO;samples];n];
    let mut impedance=vec![vec![C64::ZERO;n*n];samples];
    let mut ppw=f64::INFINITY;
    // Largest asymmetry, and the power eigenvalue with the least headroom above -asymmetry.
    let (mut asymmetry,mut worst)=(0.0_f64,(f64::INFINITY,0.0,0.0));
    for index in (0..samples).rev() {
        if gate.is_requested() { return Err("exterior frequency sweep cancelled".into()); }
        let w=omega[index]; let k=w/medium.sound_speed;
        let fields:Vec<Vec<_>>=boundary.weights.iter().map(|row|row.iter().map(|b|C64::new(0.0,*b/w)).collect()).collect();
        let refs:Vec<&[C64]>=fields.iter().map(Vec::as_slice).collect();
        let solutions=if let Some(geometry)=&geometry {geometry.solve_batch(k,medium,&refs)?}
            else {solve_radiation_batch(&surface,k,medium,&refs,formulation)?};
        if solutions.len()!=n { return Err("exterior solve lost a head source".into()); }
        for (j,solution) in solutions.iter().enumerate() {
            if !solution.panels_per_wavelength.is_finite() || solution.panels_per_wavelength<6.0
                || !solution.radiated_power_roundoff_interval.1.is_finite() || solution.radiated_power_roundoff_interval.1<0.0
                || solution.pressure.len()!=panels {
                return Err("exterior resolution, source power or pressure trace refused".into());
            }
            ppw=ppw.min(solution.panels_per_wavelength);
            let pressure=exterior_pressure_at_points_integrated(&surface,solution,medium,&[settings.receiver_m])?[0];
            // Finite-point pressure ALREADY includes propagation and spreading.
            // Peel only this lower-bound delay and restore it once at playback.
            values[j][index]=phase(pressure,-w*delay_s);
            for (i,row) in boundary.weights.iter().enumerate() {
                let mut force=C64::ZERO;
                for ((p,area),b) in solution.pressure.iter().zip(surface.areas()).zip(row) {force=force+p.scale(area*b);}
                impedance[index][i*n+j]=force*C64::new(0.0,-w);
            }
        }
        let (skew,lowest)=reciprocal_part(&mut impedance[index],n).map_err(|e|format!("{e} at {:.2} Hz",w/TAU))?;
        asymmetry=asymmetry.max(skew);
        if lowest+skew<worst.0 { worst=(lowest+skew,lowest,w/TAU); }
    }
    eprintln!("exterior reciprocity: max |Z-Z^T|/|Z|={asymmetry:.3e} (refused above 0.02); reciprocal part Re((Z+Z^T)/2) lowest eigenvalue {:.3e}|Z| at {:.2} Hz, {:.3e}|Z| above that frequency's -|Z-Z^T|/|Z| refusal",
        worst.1,worst.2,worst.0);
    if gate.is_requested() { return Err("exterior fit cancelled".into()); }
    // The fixed log-spaced dictionary cannot follow that detail over the band; the
    // travel-time dictionary does, under the same unchanged error limits.
    let spacing=(TAU*medium.sound_speed/(6.0*extent)).min(omega[samples-1]-omega[0]);
    let fitted=fit::fit_resolved(&omega,&impedance,n,spacing)?;
    let mut filters=Vec::with_capacity(n);
    for row in &values {
        if gate.is_requested() { return Err("pressure fit cancelled".into()); }
        filters.push(observer(&omega,row,settings.order)?);
    }
    eprintln!("closed kettle BEM: panels={panels}, sources={n}, samples={samples}, min_panels_per_wavelength={ppw}, formulation={formulation:?}, band_hz={:?}; extent_m={extent}, load_resonator_spacing_hz={}, passive_load_poles={}, complex_peak/RMS={}/{}, resistance_peak/RMS={}/{}",
        settings.band_hz,spacing/TAU,fitted.model.poles.len(),fitted.peak_error,fitted.rms_error,
        fitted.resistance_peak_error,fitted.resistance_rms_error);
    Ok((fitted.model,Bake {filters,delay_s}))
}
/// Keep all nonlinear, felt, air and exterior-memory state between samples.
/// No bytes are published until the complete offline render and PCM admission.
pub fn render(model:&mut Model,bake:&Bake,spec:&Spec,settings:Settings,frames:usize,gate:&CancelGate)
    ->Result<Vec<u8>,Error> {
    if frames==0 || frames>480_000 || bake.filters.len()!=model.modes.len() { return Err("invalid pressure render extent or source basis".into()); }
    settings.validate(spec)?;
    if model.mechanics.samples()!=0
        || model.mechanics.sample_period_s().to_bits()!=(1.0/f64::from(spec.rate_hz)).to_bits() { return Err("pressure playback requires quiescent prior source history".into()); }
    let count=model.modes.len(); let factor=(spec.rate_hz/OUTPUT_RATE) as usize;
    let mut filters=bake.filters.iter().map(DiscreteStateSpace::try_runtime).collect::<Result<Vec<_>,_>>()?;
    let mut delay=if bake.delay_s==0.0 {None}else{Some(DelayedFilter::new(bake.delay_s*f64::from(OUTPUT_RATE),DigitalFilter {
        sections:vec![],direct:1.0,t_s:1.0/f64::from(OUTPUT_RATE),prewarp:0.0 })?)};
    let mut decimator=Decimator::new(factor,count)?;
    let mut block=vec![0.0;factor*count]; let mut previous=vec![0.0;count];
    let mut pressure=Vec::with_capacity(frames);
    for _ in 0..frames {
        for frame in block.chunks_exact_mut(count) {
            model.step(gate)?;
            for (j,(value,previous)) in frame.iter_mut().zip(&mut previous).enumerate() {
                let v=model.mechanics.state()[2*(j+1)+1];
                *value=(v-*previous)*f64::from(spec.rate_hz); *previous=v;
            }
        }
        let acceleration=decimator.preview(&block)?;
        let mut p=0.0;
        for (filter,&a) in filters.iter_mut().zip(acceleration) { p+=filter.step(a)?; }
        if let Some(delay)=&mut delay {p=delay.push(p)?;}
        if !p.is_finite() { return Err("exterior pressure overflow".into()); }
        decimator.commit(); pressure.push(p);
    }
    let (wav,clips)=encode_pcm16_wav_interleaved(&pressure,OUTPUT_RATE,1,settings.full_scale_pa)?;
    if clips!=0 { return Err(format!("{clips} PCM clips: increase the explicit full-scale Pa; no clipped WAV is published").into()); }
    eprintln!("kettle WAV: frames={frames}, rate_hz={OUTPUT_RATE}, full_scale_pa={}, peak_pa={}, propagation_delay_s={}, anti_alias_delay_frames={}; stationary finite-point pressure, not measured-instrument or real-time certification",
        settings.full_scale_pa,pressure.iter().map(|p|p.abs()).fold(0.0_f64,f64::max),bake.delay_s,decimator.delay_output_frames());
    Ok(wav)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_acceleration_and_peeled_delay_preserve_units_phase_and_gain() {
        for w in [20.0,300.0,1500.0] {
            assert!((C64::new(0.0,2.0/w)*C64::new(0.0,-w)-C64::from_re(2.0)).abs()<1e-14);
            let original=C64::new(0.3,-0.7);
            assert!((phase(phase(original,-w*0.003),w*0.003)-original).abs()<1e-14);
        }
    }
    #[test]
    fn power_admission_checks_cross_mode_combinations() {
        let positive=vec![C64::from_re(2.0),C64::from_re(1.0),C64::from_re(1.0),C64::from_re(2.0)];
        assert!(admit_power(&positive,2).is_ok());
        let active=vec![C64::from_re(1.0),C64::from_re(2.0),C64::from_re(2.0),C64::from_re(1.0)];
        assert!(admit_power(&active,2).is_err()); // both individual diagonals are positive
    }
    #[test]
    fn reciprocal_part_admits_collocation_skew_and_refuses_beyond_it() {
        // A passive symmetric pair; collocation adds a small skew reactance.
        let pair=|r22:f64,skew:f64| vec![C64::new(2.0,-30.0),C64::new(0.0,-4.0+skew),
            C64::new(0.0,-4.0-skew),C64::new(r22,-20.0)];
        let mut z=pair(1e-6,0.1);
        assert!(admit_power(&z,2).is_err()); // the skew alone makes the Hermitian part indefinite
        let (skew,lowest)=reciprocal_part(&mut z,2).unwrap();
        assert!((0.007..0.008).contains(&skew) && lowest>0.0);
        assert!(z.iter().zip(pair(1e-6,0.0)).all(|(a,b)|(*a-b).abs()<1e-15));
        // A negative eigenvalue within the stated error is admitted and reported, not clipped.
        let (_,lowest)=reciprocal_part(&mut pair(-1e-3,0.1),2).unwrap();
        assert!(lowest<0.0 && (lowest*36.55+1e-3).abs()<1e-5);
        assert!(reciprocal_part(&mut pair(-1.0,0.1),2).is_err());
        assert!(reciprocal_part(&mut pair(1e-6,1.0),2).is_err());
    }
    #[test]
    fn zero_training_cannot_hide_unseen_audit_pressure() {
        let omega:Vec<_>=(0..33).map(|i|TAU*(40.0+10.0*i as f64)).collect();
        let mut values=vec![C64::ZERO;33];
        assert!(observer(&omega,&values,2).is_ok());
        values[1]=C64::from_re(1.0);
        assert!(observer(&omega,&values,2).is_err());
    }
    #[test]
    fn known_stable_pressure_transfer_keeps_phase_and_rejects_audit_corruption() {
        let original=DiscreteStateSpace {n:2,a:vec![0.95,0.0,0.0,0.7],b:vec![0.2,0.1],
            c:vec![0.3,-0.4],d:0.1,e_leftover:0.0,t_s:1.0/f64::from(OUTPUT_RATE)};
        let omega:Vec<_>=(0..65).map(|i|TAU*(40.0+20.0*i as f64)).collect();
        let mut values:Vec<_>=omega.iter().map(|w|original.eval(*w).unwrap().conj()).collect();
        let f=observer(&omega,&values,2).unwrap();
        for w in [300.0,1900.0,7000.0] { assert!((f.eval(w).unwrap()-original.eval(w).unwrap()).abs()<1e-5); }
        values[1]=values[1].scale(100.0);
        assert!(observer(&omega,&values,2).is_err());
    }

}
