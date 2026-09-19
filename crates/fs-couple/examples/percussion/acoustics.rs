//! Cold geometry -> BEM -> fixed-receiver filters, then observed acceleration -> Pa.
//!
//! This example is composition, not another BEM, constitutive law or integrator.
//! The fixed direction is projected BEFORE vector fitting: one proper filter per
//! generalized acceleration, not an entire SH bank for a stationary microphone.
//! A moving receiver needs the existing broadband directional source instead.
//!
//! Exterior acoustics is linear, one-way and evaluated on undeformed geometry.
//! The far-field approximation is explicit; a nearby drum microphone, radiation
//! loading, air absorption, room reflections and stand/stick radiation are NOT
//! modeled. Fit error is checked only in the declared 40..1640 Hz band. Neither
//! that check nor successful WAV export certifies a full-band physical instrument.
use super::{Error, Experiment};
use fs_bem::helmholtz::{Formulation, Medium, far_field, solve_radiation_batch};
use fs_bem::panel3d::SpherePanels;
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_plate::{ModePair, shell::head::TensionedDisk};
use fs_plate::shell::reduction::radiation::ShellRadiationSurface;
use fs_vfit::discretize::{DelayedFilter, DigitalFilter, DiscreteStateSpace, DiscreteStateSpaceRuntime, bilinear_state_space};
use fs_vfit::vf::{FitOptions, WeightPreset, vector_fit};
use std::collections::BTreeMap;

pub const OUTPUT_RATE: u32 = 48_000;
pub const SUBSTEPS: usize = 16;
pub const MECHANICAL_DT: f64 = 1.0 / (OUTPUT_RATE as f64 * SUBSTEPS as f64);
const MAX_PANELS: usize = 2048;
const MAX_INPUTS: usize = 63;
const MAX_RELATIVE_ERROR: f64 = 0.15;
const MAX_RMS_ERROR: f64 = 0.05;

/// Arbitrary closed exterior triangles, prescribed modal normal velocities and
/// explicit addresses in ImpactSystem's interleaved q/v state. No instrument
/// name enters the acoustic solve. The striker is deliberately not a source.
pub struct Boundary {
    triangles: Vec<[[f64; 3]; 3]>,
    weights: Vec<Vec<f64>>,
    state_modes: Vec<usize>,
}
impl Boundary {
    pub fn shell(surface: &ShellRadiationSurface, first_mode: usize) -> Result<Self, Error> {
        let count = surface.normal_velocity_weights().len();
        let end=first_mode.checked_add(count).ok_or("shell acoustic mode address overflow")?;
        if count==0 || count>MAX_INPUTS || end>64 || surface.triangles().len()>MAX_PANELS {
            return Err("shell acoustic boundary exceeds its mode/panel limits".into());
        }
        Ok(Self { triangles: surface.triangles().to_vec(), weights: surface.normal_velocity_weights().to_vec(),
            state_modes: (first_mode..end).collect() })
    }

    /// Two actual film meshes, rigid bearing-edge annuli and a rigid outer shell.
    /// Both mechanical head coordinates are positive DOWNWARD; global z is up.
    /// Thus top normal velocity is -w_dot, bottom normal velocity is +w_dot.
    /// Meshes must match at this reference seam; no remeshing/interpolation hides
    /// a different physical surface. Rim and shell motion is explicitly zero.
    pub fn drum(films: &[TensionedDisk], modes: &[Vec<ModePair>], depth: f64, outer_radius: f64) -> Result<Self, Error> {
        if films.len()!=2 || modes.len()!=2 || films[0].mesh.nodes!=films[1].mesh.nodes
            || films[0].mesh.tris!=films[1].mesh.tris || !depth.is_finite() || depth<=0.0
            || !outer_radius.is_finite() || outer_radius<=films[0].spec.radius_m
        { return Err("drum exterior needs two matching films, positive depth and an outer radius beyond their span".into()); }
        let mesh=&films[0].mesh; let n=mesh.nodes.len();
        let count=modes[0].len().checked_add(modes[1].len()).ok_or("drum mode count overflow")?;
        if count==0 || count>MAX_INPUTS || modes.iter().zip(films).any(|(ms,film)| ms.iter().any(|m|
            m.phi.len()!=film.model.free || m.phi.iter().any(|v| !v.is_finite()))) {
            return Err("drum radiation modes must match their film pencils and input budget".into());
        }
        let mut uses=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
        for tri in &mesh.tris { for a in 0..3 {
            let (i,j)=(tri[a],tri[(a+1)%3]); uses.entry((i.min(j),i.max(j))).or_default().push((i,j));
        }}
        let edges:Vec<_>=uses.values().filter(|e|e.len()==1).map(|e|e[0]).collect();
        if edges.len()<3 { return Err("drum film has no boundary loop".into()); }
        let mut rim:Vec<_>=edges.iter().flat_map(|&(a,b)|[a,b]).collect(); rim.sort_unstable(); rim.dedup();
        // Long axial strips would underresolve BEM even when the films are fine.
        let requested_axial=(depth/(2.0*std::f64::consts::PI*outer_radius/rim.len() as f64)).ceil().max(1.0);
        if !requested_axial.is_finite() || requested_axial>MAX_PANELS as f64 {
            return Err("drum axial subdivision exceeds its explicit panel budget".into());
        }
        let axial=requested_axial as usize;
        let panels=mesh.tris.len().checked_mul(2).and_then(|x|edges.len().checked_mul(4+2*axial).and_then(|r|x.checked_add(r)))
            .ok_or("drum exterior panel count overflow")?;
        if panels>MAX_PANELS { return Err("drum exterior exceeds its explicit panel budget".into()); }
        let mut points:Vec<_>=mesh.nodes.iter().map(|&(x,y)|[x,y,0.5*depth])
            .chain(mesh.nodes.iter().map(|&(x,y)|[x,y,-0.5*depth])).collect();
        let mut rings=Vec::new();
        for level in 0..=axial {
            let mut ring=BTreeMap::new();
            for &i in &rim {
                let (x,y)=mesh.nodes[i]; let r=x.hypot(y);
                if (r-films[0].spec.radius_m).abs()>1e-10*films[0].spec.radius_m {
                    return Err("drum exterior requires one circular film boundary".into());
                }
                ring.insert(i,points.len());
                points.push([outer_radius*x/r,outer_radius*y/r,depth*(0.5-level as f64/axial as f64)]);
            }
            rings.push(ring);
        }
        let mut indices=Vec::with_capacity(panels);
        for &tri in &mesh.tris { indices.push(tri); }
        for &tri in &mesh.tris { indices.push([tri[2]+n,tri[1]+n,tri[0]+n]); }
        for &(a,b) in &edges {
            let (ta,tb)=(rings[0][&a],rings[0][&b]);
            indices.extend([[a,ta,tb],[a,tb,b]]);
            let (ba,bb)=(rings[axial][&a],rings[axial][&b]);
            indices.extend([[a+n,bb,ba],[a+n,b+n,bb]]);
            for level in 0..axial {
                let (ta,tb,ba,bb)=(rings[level][&a],rings[level][&b],rings[level+1][&a],rings[level+1][&b]);
                indices.extend([[ta,ba,bb],[ta,bb,tb]]);
            }
        }
        check_closed(&indices)?;
        let triangles:Vec<_>=indices.iter().map(|t|t.map(|i|points[i])).collect();
        let mut weights=vec![vec![0.0;panels];count]; let mut at=0;
        for head in 0..2 {
            for mode in &modes[head] {
                for (f,tri) in mesh.tris.iter().enumerate() {
                    let shape=tri.iter().map(|&i|films[head].model.dof_map[3*i].map_or(0.0,|k|mode.phi[k])/3.0).sum::<f64>();
                    weights[at][head*mesh.tris.len()+f]=if head==0 {-shape}else{shape};
                }
                at+=1;
            }
        }
        Ok(Self {triangles,weights,state_modes:(1..=count).collect()})
    }
}
fn check_closed(tris:&[[usize;3]])->Result<(),Error> {
    let mut edges=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
    for t in tris { for i in 0..3 {
        let (a,b)=(t[i],t[(i+1)%3]); edges.entry((a.min(b),a.max(b))).or_default().push((a,b));
    }}
    if edges.values().any(|e| e.len()!=2 || e[0]!=(e[1].1,e[1].0)) {
        return Err("drum acoustic boundary is not consistently closed".into());
    }
    Ok(())
}

/// e^{-i omega t}: acceleration=-i omega velocity, so velocity=i acceleration/omega.
fn acceleration_fields(weights:&[Vec<f64>],omega:f64)->Vec<Vec<C64>> {
    weights.iter().map(|row|row.iter().map(|&b|C64::new(0.0,b/omega)).collect()).collect()
}
fn zero_filter(dt:f64)->DiscreteStateSpace {
    DiscreteStateSpace{n:0,a:vec![],b:vec![],c:vec![],d:0.0,e_leftover:0.0,t_s:dt}
}
/// Fit only even-index samples. Odd frequencies are never supplied to fs-vfit.
fn fit_observer(omega:&[f64],values:&[C64],dt:f64,order:usize)->Result<(DiscreteStateSpace,f64,f64),Error> {
    if omega.len()!=values.len() || omega.len()<5 || !dt.is_finite() || dt<=0.0
        || omega.iter().enumerate().any(|(i,w)| !w.is_finite() || *w<=0.0 || *w*dt>=std::f64::consts::PI
            || (i>0 && *w<=omega[i-1]))
        || values.iter().any(|v|!v.re.is_finite() || !v.im.is_finite()) {
        return Err("observer samples must be finite, ordered and below Nyquist".into());
    }
    let scale=values.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale==0.0 { return Ok((zero_filter(dt),0.0,0.0)); }
    if !scale.is_finite() {return Err("observer response scale overflow".into());}
    let frequencies:Vec<_>=omega.iter().step_by(2).map(|w|2.0/dt*det::tan(w*dt/2.0)).collect();
    // The fit owner uses s=+i omega. Conjugation and frequency warping are both
    // required; fitting raw negative-time phasors reverses the physical phase.
    let response:Vec<_>=values.iter().step_by(2).map(|v|C64::new(v.re/scale,-v.im/scale)).collect();
    let fit=vector_fit(&frequencies,&response,&FitOptions{order,iterations:12,weights:WeightPreset::Uniform,fit_d:true,fit_e:false})?;
    if !fit.model.is_stable() || fit.model.e!=0.0 || !fit.report.weighted_rms.is_finite() {
        return Err("observer fit must be stable, finite and proper".into());
    }
    let mut filter=bilinear_state_space(&fit.model,dt,0.0)?;
    for c in &mut filter.c {*c*=scale;} filter.d*=scale;
    // Realization admission also checks finite scaling and any improper leftover.
    filter.try_runtime()?;
    let mut maximum=0.0_f64; let mut square=0.0; let mut held=0;
    for i in (1..omega.len()).step_by(2) {
        let error=(filter.eval(omega[i])?.conj()-values[i]).abs()/scale;
        if !error.is_finite() {return Err("observer held-out error is nonfinite".into());}
        maximum=maximum.max(error); square+=error*error; held+=1;
    }
    let rms=(square/held as f64).sqrt();
    if maximum>MAX_RELATIVE_ERROR || rms>MAX_RMS_ERROR {
        return Err(format!("observer holdout refused: max={maximum}, rms={rms}; limits={MAX_RELATIVE_ERROR}/{MAX_RMS_ERROR}; refine the declared bake, do not bypass the gate").into());
    }
    Ok((filter,maximum,rms))
}

struct Bake {
    filters:Vec<DiscreteStateSpace>,
    range_m:f64,
    medium:Medium,
    source_delay_s:f64,
}
// In negative-time phasors this factor DELAYS the transfer by radius/c.
// F alone is referenced to the origin and can contain advance from the near
// side of the surface. Delay it to the enclosing sphere before a causal fit;
// the explicit propagation line then travels only (range-radius)/c.
fn shift_to_enclosing_sphere(amplitude:C64,omega:f64,source_delay_s:f64)->C64 {
    let phase=omega*source_delay_s;
    amplitude*C64::new(det::cos(phase),det::sin(phase))
}
fn bake(boundary:&Boundary,observer:[f64;3])->Result<Bake,Error> {
    let count=boundary.weights.len(); let panels=boundary.triangles.len();
    if count==0 || count>MAX_INPUTS || count!=boundary.state_modes.len() || panels==0 || panels>MAX_PANELS
        || boundary.weights.iter().any(|b|b.len()!=panels || b.iter().any(|v|!v.is_finite())) {
        return Err("exterior radiation exceeds its shape/work contract".into());
    }
    let range_m=observer.iter().map(|x|x*x).sum::<f64>().sqrt();
    let radius=boundary.triangles.iter().flatten().map(|p|p.iter().map(|x|x*x).sum::<f64>().sqrt()).fold(0.0_f64,f64::max);
    if !range_m.is_finite() || !radius.is_finite() || radius<=0.0 || range_m<10.0*radius {
        return Err("fixed far-field observer must be finite and at least ten source radii away (a screening rule, not a near-field certificate)".into());
    }
    let medium=Medium::air();
    let source_delay_s=radius/medium.sound_speed;
    let surface=SpherePanels::from_triangles(boundary.triangles.clone())?;
    let omega:Vec<_>=(0..41).map(|i|2.0*std::f64::consts::PI*(40.0+40.0*i as f64)).collect();
    let mut values=vec![Vec::with_capacity(omega.len());count];
    let mut ppw=f64::INFINITY; let mut condition=0.0_f64;
    for &w in &omega {
        let k=w/medium.sound_speed;
        let fields=acceleration_fields(&boundary.weights,w);
        let formulation=if k*radius<0.5 {Formulation::PlainCbie}else{Formulation::BurtonMiller};
        let field_refs:Vec<&[C64]>=fields.iter().map(Vec::as_slice).collect();
        let solutions=solve_radiation_batch(&surface,k,medium,&field_refs,formulation)?;
        for (row,solution) in values.iter_mut().zip(&solutions) {
            // A materially negative power result is not repaired into audible data.
            if solution.radiated_power_roundoff_interval.1<0.0 {
                return Err("BEM reports negative radiation power beyond roundoff; refine the acoustic solve".into());
            }
            ppw=ppw.min(solution.panels_per_wavelength); condition=condition.max(solution.condition_lower_bound);
            row.push(shift_to_enclosing_sphere(far_field(&surface,solution,medium,&[observer])[0],w,source_delay_s));
        }
    }
    let mut filters=Vec::with_capacity(count); let mut maximum=0.0_f64; let mut rms=0.0_f64;
    for row in &values {
        let (filter,m,r)=fit_observer(&omega,row,1.0/f64::from(OUTPUT_RATE),8)?;
        maximum=maximum.max(m); rms=rms.max(r); filters.push(filter);
    }
    eprintln!("fixed receiver BEM bake: panels={panels}, inputs={count}, band_hz=40..1640, training=21, held_out=20, min_panels_per_wavelength={ppw}, condition_lower_bound_max={condition}, max_error={maximum}, worst_input_rms={rms}");
    eprintln!("observer_m={observer:?}; one causal filter per modal acceleration; linear undeformed one-way FAR FIELD, not a close microphone or radiation loading");
    Ok(Bake{filters,range_m,medium,source_delay_s})
}

struct Observer<'a> { filters:Vec<DiscreteStateSpaceRuntime<'a>>, delay:DelayedFilter, inverse_range:f64 }
impl Bake {
    fn runtime(&self)->Result<Observer<'_>,Error> {
        let dt=self.filters.first().ok_or("empty observer bank")?.t_s;
        let delay=DelayedFilter::new((self.range_m/self.medium.sound_speed-self.source_delay_s)/dt,DigitalFilter {
            sections:vec![],direct:1.0,t_s:dt,prewarp:0.0,
        })?;
        Ok(Observer {filters:self.filters.iter().map(DiscreteStateSpace::try_runtime).collect::<Result<_,_>>()?,
            delay,inverse_range:1.0/self.range_m})
    }
}
impl Observer<'_> {
    // A rendering error is terminal for this example: discard the candidate WAV
    // and this observer. No cross-owner rollback/retry claim is made here.
    fn step(&mut self,accelerations:&[f64])->Result<f64,Error> {
        if accelerations.len()!=self.filters.len() || accelerations.iter().any(|a|!a.is_finite()) {
            return Err("observer acceleration input mismatch".into());
        }
        let mut amplitude=0.0;
        for (filter,&a) in self.filters.iter_mut().zip(accelerations) {amplitude+=filter.step(a)?;}
        if !amplitude.is_finite() {return Err("observer modal sum overflow".into());}
        Ok(self.delay.push(amplitude*self.inverse_range)?)
    }
}

/// Offline reference export; all expensive geometry/BEM/fitting precedes the
/// observed samples. Mechanics still uses the allocating implicit fs-phs owner.
/// No native wall-clock or callback allocation bound is asserted.
pub fn render(experiment:&mut Experiment,frames:usize,full_scale_pa:f64)->Result<Vec<u8>,Error> {
    let boundary=experiment.acoustics.as_ref().ok_or("missing acoustic boundary")?;
    if frames==0 || frames>480000 || !full_scale_pa.is_finite() || full_scale_pa<=0.0
        || boundary.state_modes.iter().any(|&k|k>=experiment.force.len()) {
        return Err("pressure export requires bounded frames, finite positive full-scale and valid state addresses".into());
    }
    let baked=bake(boundary,[1.5,0.7,1.5])?;
    let mut observer=baked.runtime()?;
    let count=boundary.state_modes.len();
    let mut decimator=Decimator::new(SUBSTEPS,count)?;
    let mut block=vec![0.0;SUBSTEPS*count];
    let x=experiment.system.state();
    let mut previous:Vec<f64>=boundary.state_modes.iter().map(|&k|x[2*k+1]).collect();
    let mut pressure=Vec::with_capacity(frames); let gate=CancelGate::new_clock_free();
    for _ in 0..frames {
        for frame in block.chunks_exact_mut(count) {
            experiment.system.step(&experiment.force,&gate)?;
            let x=experiment.system.state();
            for ((sample,previous),&k) in frame.iter_mut().zip(&mut previous).zip(&boundary.state_modes) {
                let velocity=x[2*k+1]; *sample=(velocity-*previous)/MECHANICAL_DT; *previous=velocity;
            }
        }
        // Interval-average mechanical accelerations, not a hand-authored force
        // pulse. Filter only this observation before its lower-rate transfer;
        // the nonlinear contact, felt, membrane and energy state are untouched.
        let acceleration=decimator.preview(&block)?;
        let p=observer.step(acceleration)?; decimator.commit(); pressure.push(p);
    }
    let (wav,clips)=encode_pcm16_wav(&pressure,OUTPUT_RATE,full_scale_pa)?;
    let peak=pressure.iter().map(|p|p.abs()).fold(0.0_f64,f64::max);
    eprintln!("pressure WAV: frames={frames}, rate_hz={OUTPUT_RATE}, full_scale_pa={full_scale_pa}, clips={clips}, peak_pa={peak}, observer_flight_s={}, decimator_delay_frames={}; no peak normalization; acceleration is a step-average tagged at step end",baked.range_m/baked.medium.sound_speed,decimator.delay_output_frames());
    Ok(wav)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_plate::shell::{head::TensionedDiskSpec,profile::ProfileBudget};
    fn film()->TensionedDisk {
        TensionedDisk::new(TensionedDiskSpec{radius_m:0.1,thickness_m:0.0002,young_pa:4e9,poisson:0.38,
            density_kg_m3:1390.,tension_n_m:1000.,radial_intervals:2,azimuths:8},
            ProfileBudget{max_nodes:100,max_triangles:200,max_feature_evaluations:0}).unwrap()
    }
    fn shape(f:&TensionedDisk)->ModePair {
        let mut phi=vec![0.0;f.model.free];
        for (i,&(x,y)) in f.mesh.nodes.iter().enumerate() {if let Some(k)=f.model.dof_map[3*i] {
            phi[k]=1.0-(x*x+y*y)/f.spec.radius_m.powi(2);
        }}
        // A synthetic kinematic fixture, deliberately not an eigensolve claim.
        ModePair{lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}
    }
    #[test]
    fn drum_exterior_uses_opposite_head_normals_and_real_outer_shell() {
        let films=[film(),film()]; let modes=vec![vec![shape(&films[0])],vec![shape(&films[1])]];
        let b=Boundary::drum(&films,&modes,0.12,0.11).unwrap();
        let s=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
        assert_eq!(b.state_modes,vec![1,2]);
        for head in 0..2 {
            let flux:f64=b.weights[head].iter().zip(s.areas()).map(|(b,a)|b*a).sum();
            let expected=(if head==0 {-1.0}else{1.0})*films[head].modal_area(&modes[head][0].phi).unwrap();
            assert!((flux-expected).abs()<1e-14);
        }
        let moving=2*films[0].mesh.tris.len();
        assert!(b.weights.iter().all(|row|row[moving..].iter().all(|v|*v==0.0)));
        let volume:f64=s.centroids().iter().zip(s.normals()).zip(s.areas())
            .map(|((p,n),a)|p.iter().zip(n).map(|(p,n)|p*n).sum::<f64>()*a/3.0).sum();
        let expected=4.0*0.11_f64.powi(2)*(std::f64::consts::PI/4.0).sin()*0.12;
        assert!((volume-expected).abs()<1e-14);
    }
    #[test]
    fn acceleration_to_bem_velocity_preserves_the_time_convention() {
        let fields=acceleration_fields(&[vec![2.,-3.,0.]],10.0);
        for (&v,w) in fields[0].iter().zip([2.,-3.,0.]) {
            assert!((v*C64::new(0.,-10.)-C64::from_re(w)).abs()<1e-15);
        }
    }
    #[test]
    fn fixed_observer_fit_retains_phase_and_direct_feedthrough() {
        let original=DiscreteStateSpace{n:2,a:vec![0.95,0.,0.,0.7],b:vec![0.2,0.1],c:vec![0.3,-0.4],
            d:0.1,e_leftover:0.,t_s:1./48000.};
        let omega:Vec<_>=(0..41).map(|i|2.*std::f64::consts::PI*(40.+40.*i as f64)).collect();
        let values:Vec<_>=omega.iter().map(|w|original.eval(*w).unwrap().conj()).collect();
        let (f,max,rms)=fit_observer(&omega,&values,original.t_s,2).unwrap();
        assert!(max<1e-5 && rms<1e-5);
        for w in [300.,1900.,7000.] {assert!((f.eval(w).unwrap()-original.eval(w).unwrap()).abs()<1e-5);}
        let mut changed=values.clone(); changed[1]=changed[1].scale(100.0);
        assert!(fit_observer(&omega,&changed,original.t_s,2).is_err(),"held-out data must reject, never enter fitting");
    }
    #[test]
    fn enclosing_sphere_shift_and_remaining_propagation_preserve_total_phase() {
        let source=C64::new(0.3,-0.7);
        for omega in [1.0,200.0,9000.0] {
            let shifted=shift_to_enclosing_sphere(source,omega,0.0006);
            let propagated=shift_to_enclosing_sphere(shifted,omega,0.006-0.0006);
            let expected=shift_to_enclosing_sphere(source,omega,0.006);
            assert!((propagated-expected).abs()<1e-13);
        }
    }
    #[test]
    fn pressure_delay_and_range_are_applied_exactly_once() {
        let mut filter=zero_filter(0.001); filter.d=2.0;
        let baked=Bake{filters:vec![filter],range_m:2.0,medium:Medium{density:1.2,sound_speed:100.0},source_delay_s:0.0};
        let mut r=baked.runtime().unwrap();
        for i in 0..30 {
            let y=r.step(&[if i==0 {1.0}else{0.0}]).unwrap();
            assert_eq!(y,if i==20 {1.0}else{0.0});
        }
    }
}
