//! Maxwell stress on the SAME retained moderate-slope strain as elastic storage.
//! Each original layer/arm retains mean and cosine relaxed strains. Their
//! midpoint update and mechanical pullback close discrete work exactly for the
//! quadratic epsilon=u_x+(y_x^2+z_x^2)/2; no fitted modal damping is substituted.
use super::*;
mod maxwell;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::history::{
    MappedModalHistoryBinding, ModalMaterialHistory,
};
fn constitutive_invalid(what:&'static str) -> ModalCouplingError {ModalCouplingError::Invalid {what}}

#[derive(Clone, Debug)]
struct Arm {
    /// Delta(EA)*L. Mean storage is C*r^2/2, cosine storage C*r_j^2/4.
    stiffness_j: f64,
    time_s: f64,
}

#[derive(Clone, Debug)]
struct StrainHistory {
    string: MaterialLongitudinalString,
    length_m: f64,
    arms: Vec<Arm>,
}
impl StrainHistory {
    fn width(&self) -> usize {1+self.string.model.axial_modes().len()}
    fn stiffness_bound(&self,memory:Option<&[f64]>)->f64 {
        // Projection cannot increase the spatial L2 norm. The combined slope
        // envelope bounds both planes; no extra plane multiplier belongs in
        // B^T C B. The stress-geometric Hessian acts on each plane separately.
        let axial=0.5*self.string.model.axial_slope_scales().iter().map(|s|s*s).sum::<f64>();
        let transverse=0.5*self.string.slope_scales.iter().map(|s|s*s).sum::<f64>();
        let strain=self.string.resolution.maximum_strain;
        self.arms.iter().enumerate().map(|(i,arm)| {
            // For dt/tau<=2 each updated z is a convex combination of its
            // accepted value and the mean endpoint strain. Thus this original
            // continuous-span bound remains valid through the entire history.
            let z_bound=memory.map_or(strain,|m|m[i*self.width()..(i+1)*self.width()]
                .iter().map(|v|v.abs()).sum::<f64>().max(strain));
            arm.stiffness_j*(axial+transverse*(self.string.resolution.maximum_slope.powi(2)
                +self.string.planes as f64*(strain+z_bound)))
        }).sum()
    }
    fn strain(&self,q:&[f64]) -> Result<Vec<f64>,ModalCouplingError> {
        self.string.validate_coordinates(q)?;
        let mut strain=vec![0.0;self.width()];
        strain[0]=self.string.model.geometric_strain(&self.string.slopes(q),&mut strain[1..]).map_err(owner)?;
        for (j,&scale) in self.string.model.axial_slope_scales().iter().enumerate() {
            strain[j+1]+=scale*q[self.string.transverse_count()+j];
        }
        Ok(strain)
    }
    fn validate_memory(&self,strain:&[f64],memory:&[f64]) -> Result<(),ModalCouplingError> {
        if memory.len()!=self.memory_count() || memory.iter().any(|v|!v.is_finite()) {
            return Err(constitutive_invalid("full-strain material history requires complete finite relaxed strains"));
        }
        let mut tension:Vec<_>=strain.iter().map(|v|self.string.axial_rigidity_n*v).collect();
        tension[0]+=self.string.prestress_n;
        for (arm,z) in self.arms.iter().zip(memory.chunks_exact(self.width())) {
            for ((force,epsilon),z) in tension.iter_mut().zip(strain).zip(z) {
                *force+=arm.stiffness_j/self.length_m*(epsilon-z);
            }
        }
        if tension.iter().any(|v|!v.is_finite()) || tension[0]-tension[1..].iter().map(|v|v.abs()).sum::<f64>()<=0.0 {
            return Err(constitutive_invalid("hereditary string stress loses local tension or finite storage"));
        }
        Ok(())
    }
    fn pull(&self,q:&[f64],stress:&[f64],out:&mut[f64]) -> Result<(),ModalCouplingError> {
        if out.len()!=self.coordinate_count() {return Err(constitutive_invalid("full-strain material force shape"));}
        let slopes=self.string.slopes(q);let mut gradient=vec![0.0;slopes.len()];
        self.string.model.geometric_strain_gradient(&slopes,stress[0],&stress[1..],&mut gradient).map_err(owner)?;
        for plane in 0..self.string.planes {for mode in 0..self.string.transverse_modes {
            out[plane*self.string.transverse_modes+mode]=-self.string.slope_scales[mode]
                *gradient[plane*(self.string.transverse_modes+1)+mode+1];
        }}
        for (j,&scale) in self.string.model.axial_slope_scales().iter().enumerate() {
            out[self.string.transverse_count()+j]=-scale*stress[j+1];
        }
        if out.iter().any(|v|!v.is_finite()) {return Err(constitutive_invalid("hereditary material force overflows"));}
        Ok(())
    }
}
impl ModalMaterialHistory for StrainHistory {
    fn coordinate_count(&self) -> usize {self.string.coordinate_count()}
    fn memory_count(&self) -> usize {self.width()*self.arms.len()}
    fn stiffness_upper_bound_per_s2(&self,memory:&[f64])->Result<f64,ModalCouplingError> {
        if memory.len()!=self.memory_count() || memory.iter().any(|v|!v.is_finite()) {
            return Err(constitutive_invalid("complete initial history is required for stiffness admission"));
        }
        let bound=self.stiffness_bound(Some(memory));
        if !bound.is_finite() {return Err(constitutive_invalid("hereditary stiffness bound overflows"));}
        Ok(bound)
    }
    fn relaxed_memory(&self,q:&[f64]) -> Result<Vec<f64>,ModalCouplingError> {
        Ok(self.strain(q)?.repeat(self.arms.len()))
    }
    fn energy(&self,q:&[f64],memory:&[f64]) -> Result<f64,ModalCouplingError> {
        let strain=self.strain(q)?;self.validate_memory(&strain,memory)?;
        let mut energy=0.0;
        for (arm,z) in self.arms.iter().zip(memory.chunks_exact(self.width())) {
            for (j,(&epsilon,&z)) in strain.iter().zip(z).enumerate() {
                energy+=0.5*arm.stiffness_j*if j==0 {1.0} else {0.5}*(epsilon-z).powi(2);
            }
        }
        if !energy.is_finite() {return Err(constitutive_invalid("hereditary strain energy overflows"));}
        Ok(energy)
    }
    fn reaction(&self,q:&[f64],memory:&[f64],out:&mut[f64]) -> Result<(),ModalCouplingError> {
        let strain=self.strain(q)?;self.validate_memory(&strain,memory)?;
        let mut stress=vec![0.0;self.width()];
        for (arm,z) in self.arms.iter().zip(memory.chunks_exact(self.width())) {
            for (j,((force,epsilon),z)) in stress.iter_mut().zip(&strain).zip(z).enumerate() {
                *force+=arm.stiffness_j*if j==0 {1.0} else {0.5}*(epsilon-z);
            }
        }
        self.pull(q,&stress,out)
    }
    fn stage(&self,before:&[f64],after:&[f64],memory:&[f64],dt:f64,
        candidate:&mut[f64],reaction:&mut[f64]) -> Result<f64,ModalCouplingError> {
        let a=self.strain(before)?;let b=self.strain(after)?;self.validate_memory(&a,memory)?;
        if candidate.len()!=memory.len() || !dt.is_finite() || dt<=0.0 {
            return Err(constitutive_invalid("hereditary strain update needs complete state and positive physical duration"));
        }
        let mut stress=vec![0.0;self.width()];let mut loss=0.0;
        for ((arm,z0),z1) in self.arms.iter().zip(memory.chunks_exact(self.width()))
            .zip(candidate.chunks_exact_mut(self.width())) {
            let rate=dt/arm.time_s;
            if !rate.is_finite() || rate>2.0 {return Err(constitutive_invalid("physical clock does not resolve material relaxation"));}
            for j in 0..self.width() {
                // Use mean ENDPOINT strains here, not strain at midpoint q.
                let c=arm.stiffness_j*if j==0 {1.0} else {0.5};
                let (next,effort,heat)=maxwell::step(a[j],b[j],z0[j],rate,c)
                    .ok_or_else(||constitutive_invalid("material strain-history step overflows"))?;
                z1[j]=next;stress[j]+=effort;loss+=heat;
            }
        }
        self.validate_memory(&b,candidate)?;
        let midpoint:Vec<_>=before.iter().zip(after).map(|(&a,&b)|f64::midpoint(a,b)).collect();
        self.pull(&midpoint,&stress,reaction)?;
        if !loss.is_finite() {return Err(constitutive_invalid("hereditary dissipation overflows"));}
        Ok(loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law()->StrainHistory {
        let (length,solid,transverse,ea,tension)=(0.5,0.002,0.003,20_000.0,20.0);
        let indices=vec![1,2,3,4];
        let model=fs_nlmodal::longitudinal::LongitudinalString::new(
            fs_nlmodal::KcStringParams {length,lin_density:solid,tension,ea},2,2,&indices).unwrap();
        StrainHistory {length_m:length,arms:vec![Arm {stiffness_j:2500.0,time_s:0.001},
            Arm {stiffness_j:1000.0,time_s:0.003}],string:MaterialLongitudinalString {
                model,resolution:LongitudinalStringResolution {axial_modes:indices,maximum_slope:0.2,maximum_strain:0.01},
                transverse_modes:2,planes:2,slope_scales:(1..=2).map(|n|
                    (2.0_f64/(transverse*length)).sqrt()*n as f64*core::f64::consts::PI/length).collect(),
                solid_linear_density_kg_m:solid,transverse_linear_density_kg_m:transverse,
                axial_rigidity_n:ea,prestress_n:tension,axial_component:1,elastic_frequency_band_rad_s:None,
                axial_material_loss:None,finite_strain_viscosity:false,
            }}
    }
    fn spatial(law:&StrainHistory,q:&[f64],s:f64)->(f64,Vec<f64>) {
        let k=core::f64::consts::PI/law.length_m;
        let mut slopes=[0.0;2];let mut b=vec![0.0;q.len()];let mut strain=0.0;
        for plane in 0..2 {for mode in 0..2 {
            let shape=(2.0/(law.string.transverse_linear_density_kg_m*law.length_m)).sqrt()
                *(mode+1) as f64*k*((mode+1) as f64*core::f64::consts::PI*s).cos();
            b[2*plane+mode]=shape;slopes[plane]+=shape*q[2*plane+mode];
        }}
        for plane in 0..2 {
            strain+=0.5*slopes[plane]*slopes[plane];
            for mode in 0..2 {b[2*plane+mode]*=slopes[plane];}
        }
        for mode in 0..4 {
            let shape=(2.0/(law.string.solid_linear_density_kg_m*law.length_m)).sqrt()
                *(mode+1) as f64*k*((mode+1) as f64*core::f64::consts::PI*s).cos();
            b[4+mode]=shape;strain+=shape*q[4+mode];
        }
        (strain,b)
    }
    #[test]
    fn g1_history_force_and_loss_match_independent_spatial_work() {
        let law=law();let a=vec![2e-5,-1e-5,1e-5,0.5e-5,1e-8,-2e-8,0.0,1e-8];
        let d=[2e-6,-1e-6,-1e-6,1e-6,2e-9,-1e-9,1e-9,-2e-9];
        let b:Vec<_>=a.iter().zip(d).map(|(a,d)|a+d).collect();let dt=1e-4;
        let mut memory=law.relaxed_memory(&a).unwrap();memory[0]-=2e-6;memory[6]+=1e-7;
        let old=memory.clone();let mut next=vec![0.0;memory.len()];let mut force=vec![0.0;a.len()];
        let loss=law.stage(&a,&b,&memory,dt,&mut next,&mut force).unwrap();
        assert_eq!(memory,old,"a trial cannot advance accepted history");
        let midpoint:Vec<_>=a.iter().zip(&b).map(|(&a,&b)|f64::midpoint(a,b)).collect();
        let (mut spatial_loss,mut spatial_force)=(0.0,vec![0.0;a.len()]);
        for cell in 0..512 {
            let s=(cell as f64+0.5)/512.0;
            let ea=spatial(&law,&a,s).0;let eb=spatial(&law,&b,s).0;let gradient=spatial(&law,&midpoint,s).1;
            for (arm,z) in law.arms.iter().zip(memory.chunks_exact(law.width())) {
                let z=z[0]+z[1..].iter().enumerate().map(|(i,z)|
                    z*((i+1) as f64*core::f64::consts::PI*s).cos()).sum::<f64>();
                let rate=dt/arm.time_s;let r=(0.5*(ea+eb)-z)/(1.0+0.5*rate);
                spatial_loss+=arm.stiffness_j*rate*r*r/512.0;
                for (force,g) in spatial_force.iter_mut().zip(&gradient) {*force-=arm.stiffness_j*r*g/512.0;}
            }
        }
        assert!((loss-spatial_loss).abs()<1e-12*loss);
        for (a,b) in force.iter().zip(spatial_force) {assert!((a-b).abs()<1e-10*a.abs().max(1e-12));}
        let defect=law.energy(&b,&next).unwrap()-law.energy(&a,&memory).unwrap()+loss
            +force.iter().zip(d).map(|(f,d)|f*d).sum::<f64>();
        assert!(defect.abs()<1e-15);
        let mut repeat=vec![0.0;next.len()];let mut again=vec![0.0;force.len()];
        assert_eq!(loss,law.stage(&a,&b,&memory,dt,&mut repeat,&mut again).unwrap());
        assert_eq!(repeat,next);assert_eq!(again,force);
    }
    #[test]
    fn g1_stiffness_admission_covers_bent_geometry_and_retained_material_stress() {
        let law=law();let q=vec![2e-5,-1e-5,1e-5,0.5e-5,1e-8,-2e-8,0.0,1e-8];
        let mut memory=law.relaxed_memory(&q).unwrap();memory[0]-=2e-4;
        let bound=law.stiffness_bound(Some(&memory));
        let axial: f64=law.arms.iter().map(|arm|arm.stiffness_j*0.5
            *law.string.model.axial_slope_scales().iter().map(|s|s*s).sum::<f64>()).sum();
        assert!(bound>axial,"a bent string's added stiffness includes transverse reactions");
        let direction=[0.2,-0.3,0.1,0.4,0.15,-0.25,0.35,-0.1];let h=1e-10;
        let plus:Vec<_>=q.iter().zip(direction).map(|(q,d)|q+h*d).collect();
        let minus:Vec<_>=q.iter().zip(direction).map(|(q,d)|q-h*d).collect();
        let mut fp=vec![0.0;q.len()];let mut fm=fp.clone();
        law.reaction(&plus,&memory,&mut fp).unwrap();law.reaction(&minus,&memory,&mut fm).unwrap();
        let tangent:f64=fp.iter().zip(fm).zip(direction).map(|((p,m),d)|-d*(p-m)/(2.0*h)).sum();
        assert!(tangent.abs()<=bound*direction.iter().map(|d|d*d).sum::<f64>());
        let mut larger=memory.clone();larger[0]-=0.1;
        assert!(law.stiffness_bound(Some(&larger))>bound,"explicit initial stress must widen admission");
    }
}

impl ReciprocalStringBank {
    /// Source-defined axial Maxwell spectra acting on the complete retained
    /// strain, including its mean geometric stretch and both transverse planes.
    /// The original joint single/multi-contact solve stages all memory, force,
    /// storage, dissipation and endpoint pressure; publication commits it once.
    ///
    /// Bind after longitudinal waves and static initialization, before playback.
    /// One state per mean/cosine strain, in selection/layer/nonzero-arm order.
    /// Explicit memory is relaxed strain (dimensionless); None is relaxed at
    /// the actual initial q. The maintained reference prestress does not relax.
    /// This replaces the linear axial-memory approximation. Parallel bending,
    /// torsion and independently selected Kelvin--Voigt mechanisms remain intact.
    /// No live gas/topology remap, self-heating, slip or finite-strain geometry
    /// beyond the original moderate-slope model is claimed.
    pub fn with_material_finite_strain_longitudinal_relaxation(mut self,
        selections:Vec<LongitudinalRelaxationSelection>,initial_memory:Option<Vec<f64>>,
        max_memory_states:usize,gate:&CancelGate) -> Result<Self,ImpactError> {
        let (image,bindings)=bindings(&self.source,&self.cards,&self.system,selections,
            initial_memory,max_memory_states,gate)?;
        if !bindings.is_empty() {
            self.system=self.system.with_mapped_material_histories(bindings,max_memory_states,gate)
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
        }
        retain_source(&mut self.source,image);
        validate_sources(&self.source,&self.system,gate)?;
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok(self)
    }
}

pub(in crate::bowed_string::runtime::reciprocal::bank) fn retain_source(
    source:&mut MaterialStringBankSource,image:StringBankLongitudinalRelaxation) {
    source.longitudinal_waves.as_mut().expect("validated original waves").material_relaxation=Some(image);
}

/// Re-admit original material participants against the SAME expanded band.
pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_sources(
    source:&MaterialStringBankSource,system:&FrictionModalSystem,gate:&CancelGate) -> Result<(),ImpactError> {
    source.validate_body_memory_frequencies(system,gate).map_err(bank_error)?;
    source.validate_material_bridge_frequencies(system,gate).map_err(bank_error)?;
    let axial=source.longitudinal_waves.as_ref().expect("validated original waves");
    super::super::bearings::validate_bands(source,&axial.strings,&axial.bearings,system,gate)?;
    let high=system.linear_frequency_upper_bound_rad_s();
    for string in &axial.strings {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        let loss=string.axial_material_loss.as_ref().map(|law|law.omega_band_rad_s());
        if string.elastic_frequency_band_rad_s.into_iter().chain(loss).any(|(_,hi)|high>hi) {
            return Err(invalid("hereditary coupled spectrum exceeds original axial elasticity or viscosity support"));
        }
    }
    Ok(())
}

/// Common source projection for single and already-preloaded multi-contact owners.
pub(in crate::bowed_string::runtime::reciprocal::bank) fn bindings(
        source:&MaterialStringBankSource,cards:&[crate::bowed_string::BowedStringCard],
        system:&FrictionModalSystem,selections:Vec<LongitudinalRelaxationSelection>,
        initial_memory:Option<Vec<f64>>,max_memory_states:usize,gate:&CancelGate)
    -> Result<(StringBankLongitudinalRelaxation,Vec<MappedModalHistoryBinding>),ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if system.samples_rendered()!=0 || source.longitudinal_waves.as_ref()
            .is_none_or(|axial|axial.material_relaxation.is_some()) {
            return Err(invalid("full-strain memory requires one cold binding to original axial waves"));
        }
        let mut image=StringBankLongitudinalRelaxation::bind_source(source,cards,system,selections,max_memory_states,gate)?;
        image.finite_strain=true;
        let axial=source.longitudinal_waves.as_ref().expect("checked original waves");
        let mut laws=Vec::new();let mut count=0usize;
        for (selection,section) in image.selections.iter().zip(&image.sections) {
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            let mut arms=Vec::new();
            let length=cards[selection.string].length_m;
            for (spectrum,geometry) in section.spectra().iter().zip(section.elastic().layer_integrals()) {
                for &(delta,time_s) in &spectrum.terms {
                    if delta==0.0 {continue;}
                    let stiffness_j=delta*geometry.area_m2()*length;
                    if !stiffness_j.is_finite() || stiffness_j<=0.0 {return Err(invalid("finite-strain memory rigidity overflows"));}
                    arms.push(Arm {stiffness_j,time_s});
                }
            }
            if arms.is_empty() {continue;}
            let law=StrainHistory {string:axial.strings[selection.string].clone(),length_m:length,arms};
            count=count.checked_add(law.memory_count()).ok_or_else(||invalid("strain history count overflow"))?;
            if count>max_memory_states {return Err(invalid("mean and harmonic histories exceed the full-strain memory budget"));}
            laws.push((selection.string,law));
        }
        if initial_memory.as_ref().is_some_and(|m|m.len()!=count || m.iter().any(|v|!v.is_finite())) {
            return Err(invalid("complete finite mean/harmonic initial material history required"));
        }
        let mut start=0;let mut bindings=Vec::with_capacity(laws.len());
        for (string,law) in laws {
            let end=start+law.memory_count();
            image.finite_strain_stiffness_bound_per_s2+=law.stiffness_bound(
                initial_memory.as_ref().map(|m|&m[start..end]));
            let coordinates=(0..law.string.transverse_count()).map(|mode|ModalCoordinate {component:string,mode})
                .chain((0..law.string.model.axial_modes().len()).map(|mode|ModalCoordinate {
                    component:law.string.axial_component,mode})).collect();
            bindings.push(MappedModalHistoryBinding {coordinates,
                initial_memory:initial_memory.as_ref().map(|m|m[start..end].to_vec()),law:Box::new(law)});
            start=end;
        }
        // Complete source/clock admission follows binding, using the shared
        // owner's expanded frequency report exactly once for every material.
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok((image,bindings))
}
