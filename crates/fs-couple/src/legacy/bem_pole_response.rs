//! Force participation and analytic gas/material sensitivities of loaded modes.
//!
//! The existing BEM model owns the entire pressure/displacement/Maxwell pencil.
//! This adapter evaluates it ONCE at a resolved pole and uses fs-la's prepared
//! bordered factor for all forcing columns and parameter directions. No extra
//! frequency solve, material fit, pole extraction or boundary discretization.
//!
//! Derivatives hold the source basis, geometry and pressure-work maps fixed.
//! A density/sound-speed direction can be composed with the caller's actual
//! temperature, pressure or composition derivatives; no atmosphere or equation
//! of state is inferred. Coefficient directions must come from the same source
//! coordinates. These are local derivatives, not finite-change predictions,
//! mode tracking at crossings or derivatives of omitted structural modes.
use crate::bem_loaded_modes::{BemLoadedModes,BemModePole,BemModesError};
use fs_la::nonlinear_eig::PoleEvaluation;
use fs_la::pole_response::{PoleResponseError,PoleResponseOptions,PreparedPole};
use fs_math::c64::C64;

#[derive(Debug,Clone,Copy)]
pub struct BemPoleResponseOptions {
    pub numerical: PoleResponseOptions,
    /// Original mechanical, pressure and memory residual gate at preparation.
    pub physical_residual_tolerance: f64,
    /// Additional original-operator/action/output payload, excluding the source
    /// model and separately budgeted numerical factor, allocator overhead/RSS.
    pub max_operator_bytes: usize,
    /// Per-batch action/observation units 32*max_vectors*dimension^2,
    /// additional to the numerical solve budget. Not a measured FLOP count.
    pub max_action_work: u64,
}

#[derive(Debug)]
pub enum BemPoleResponseError {
    Invalid { what: &'static str },
    Model(BemModesError),
    Numerical(PoleResponseError),
    PhysicalResidual { mechanical:f64,boundary:f64,memory:f64 },
    Cancelled,
}
impl core::fmt::Display for BemPoleResponseError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"FS-COUPLE-BEM-POLE-RESPONSE: {self:?}")
    }
}
impl std::error::Error for BemPoleResponseError {}
impl From<BemModesError> for BemPoleResponseError {fn from(e:BemModesError)->Self{Self::Model(e)}}
impl From<PoleResponseError> for BemPoleResponseError {fn from(e:PoleResponseError)->Self{Self::Numerical(e)}}
fn invalid(what:&'static str)->BemPoleResponseError{BemPoleResponseError::Invalid{what}}
fn finite(z:C64)->bool{z.re.is_finite()&&z.im.is_finite()}
fn poll(c:&mut impl FnMut()->bool)->Result<(),BemPoleResponseError>{
    if c(){Err(BemPoleResponseError::Cancelled)}else{Ok(())}
}
fn add(sum:&mut C64,correction:&mut C64,term:C64){
    let y=term-*correction;let next=*sum+y;*correction=(next-*sum)-y;*sum=next;
}

/// Derivatives with respect to ONE explicitly chosen real parameter theta.
/// Empty slices mean exactly zero change, never missing source data. Nonempty
/// matrices/rows must be complete and finite. Entries are per unit theta.
#[derive(Debug,Clone,Copy,Default)]
pub struct BemParameterDirection<'a> {
    pub density_kg_m3: f64,
    pub sound_speed_m_s: f64,
    /// Derivative of mass in the FIXED original unit-mass coordinates; m*m.
    pub mass: &'a [f64],
    pub stiffness_per_s2: &'a [f64],
    pub viscosity_per_s: &'a [f64],
    /// Row-major derivative of c_j=sqrt(stiffness_j)*projection_j; branches*m.
    /// For a log branch-stiffness change with fixed projection, dc_j=c_j/2.
    pub relaxation_rows_per_s: &'a [f64],
    pub relaxation_times_s: &'a [f64],
}

/// Laurent coefficient in omega: response(omega) ~ residue/(omega-omega_p).
/// Every field therefore has its ordinary physical units TIMES rad/s. These
/// fields are not RMS amplitudes, energy, power, or a standalone causal model.
#[derive(Debug,Clone)]
pub struct BemForceResidue {
    pub modal_displacement_residue: Vec<C64>,
    pub surface_pressure_residue: Vec<C64>,
    pub elastic_memory_residue: Vec<C64>,
    pub surface_normal_velocity_residue: Vec<C64>,
    pub linear_residual: f64,
}

/// Pole and original-coordinate mode-shape derivative per unit theta.
/// Gauge: one coordinate of the prepared reference mode is fixed at one.
#[derive(Debug,Clone)]
pub struct BemPoleDerivative {
    pub angular_frequency_rad_s: C64,
    pub modal_displacement: Vec<C64>,
    pub surface_pressure: Vec<C64>,
    pub elastic_memory: Vec<C64>,
    pub linear_residual: f64,
}
impl BemPoleDerivative {
    pub fn frequency_hz(&self)->f64{self.angular_frequency_rad_s.re/core::f64::consts::TAU}
    pub fn decay_rate_per_s(&self)->f64{-self.angular_frequency_rad_s.im}
}

/// Borrows the immutable physical model: the factor cannot be cross-wired with
/// another source or medium. Source geometry/material validity remains upstream.
pub struct BemPoleResponse<'a> {
    model:&'a BemLoadedModes,
    raw:PoleEvaluation,
    prepared:PreparedPole,
    maximum_vectors:usize,
    physical_residuals:[f64;3],
}

impl BemLoadedModes {
    /// Prepare forcing residues and gas/material derivatives of this model's
    /// isolated pole. Reported fields/residuals in the supplied pole are NOT
    /// trusted: the full original equations are evaluated and checked again.
    pub fn prepare_pole_response<'a>(
        &'a self,pole:&BemModePole,options:BemPoleResponseOptions,
        mut cancelled:impl FnMut()->bool,
    )->Result<BemPoleResponse<'a>,BemPoleResponseError>{
        poll(&mut cancelled)?;
        let (m,p,r,n)=(self.mode_count(),self.panel_count(),self.memory_count(),self.dimension());
        let count=options.numerical.max_vectors as u128;
        let size=n as u128;
        if pole.modal_displacement.len()!=m||pole.surface_pressure_pa.len()!=p
            ||pole.elastic_memory_sqrt_j.len()!=r||!finite(pole.angular_frequency_rad_s)
            ||pole.angular_frequency_rad_s.im>0.0
            ||!options.physical_residual_tolerance.is_finite()
            ||options.physical_residual_tolerance<=0.0||options.physical_residual_tolerance>=1.0
            ||64*size*size+128*size*(count+4)>options.max_operator_bytes as u128
            ||32*count*size*size>options.max_action_work as u128
        {return Err(invalid("complete stable pole, physical tolerance and operator/action budgets required"));}
        let mut vector=pole.modal_displacement.clone();
        vector.extend_from_slice(&pole.surface_pressure_pa);
        vector.extend_from_slice(&pole.elastic_memory_sqrt_j);
        let raw=self.original_pencil(pole.angular_frequency_rad_s,&mut cancelled)?;
        let prepared=PreparedPole::new(pole.angular_frequency_rad_s,&vector,&raw,
            options.numerical,&mut cancelled)?;
        let residuals=physical_rows(self,&raw,prepared.right_vector(),pole.angular_frequency_rad_s,&mut cancelled)?;
        if residuals.iter().any(|&v|v>options.physical_residual_tolerance){
            return Err(BemPoleResponseError::PhysicalResidual{
                mechanical:residuals[0],boundary:residuals[1],memory:residuals[2]});
        }
        poll(&mut cancelled)?;
        Ok(BemPoleResponse{model:self,raw,prepared,maximum_vectors:options.numerical.max_vectors,
            physical_residuals:residuals})
    }
}

fn physical_rows(model:&BemLoadedModes,raw:&PoleEvaluation,vector:&[C64],w:C64,
    cancelled:&mut impl FnMut()->bool)->Result<[f64;3],BemPoleResponseError>{
    let (m,p,n)=(model.mode_count(),model.panel_count(),model.dimension());
    let mut worst=[0.0_f64;3];
    for i in 0..n {
        poll(cancelled)?;
        let(mut sum,mut correction)=(C64::ZERO,C64::ZERO);let mut denominator=0.0;
        for j in 0..n {
            add(&mut sum,&mut correction,raw.matrix[i*n+j]*vector[j]);
            // Do not normalize a near-zero scalar stiffness by itself. Keep
            // stiffness/inertia/viscosity and the memory diagonal separated.
            let magnitude=if i<m&&j<m {
                model.stiffness_per_s2()[i*m+j].abs()+w.abs()*model.viscosity_per_s()[i*m+j].abs()
                    +if i==j{w.abs()*w.abs()}else{0.0}
            }else if i>=m+p&&j==i{
                1.0+w.abs()*model.relaxation_times_s()[i-m-p]
            }else{raw.matrix[i*n+j].abs()};
            denominator+=magnitude*vector[j].abs();
        }
        if !finite(sum)||!denominator.is_finite(){return Err(invalid("original pole row overflow"));}
        let relative=if denominator==0.0{if sum.abs()==0.0{0.0}else{f64::INFINITY}}else{sum.abs()/denominator};
        let block=if i<m{0}else if i<m+p{1}else{2};
        worst[block]=worst[block].max(relative);
    }
    Ok(worst)
}

impl BemPoleResponse<'_> {
    pub fn angular_frequency_rad_s(&self)->C64{self.prepared.value()}
    /// Full original coordinate order [modal displacement, panel pressure, memory].
    pub fn reference_mode(&self)->&[C64]{self.prepared.right_vector()}
    pub fn normalization_pivot(&self)->usize{self.prepared.normalization_pivot()}
    pub fn physical_residuals(&self)->[f64;3]{self.physical_residuals}

    fn admit_count(&self,count:usize)->Result<(),BemPoleResponseError>{
        if count==0||count>self.maximum_vectors{return Err(invalid("nonempty batch within declared vector ceiling required"));}
        Ok(())
    }
    /// Original modal force columns, evaluated at this pole. No extra BEM solve
    /// or LU; all pressure and material residues follow the same coupled state.
    pub fn force_residues(&self,forces:&[&[C64]],mut cancelled:impl FnMut()->bool)
        ->Result<Vec<BemForceResidue>,BemPoleResponseError>{
        poll(&mut cancelled)?;self.admit_count(forces.len())?;
        let(m,p,n)=(self.model.mode_count(),self.model.panel_count(),self.model.dimension());
        let mut full=Vec::with_capacity(forces.len());
        for &f in forces{
            poll(&mut cancelled)?;
            if f.len()!=m||f.iter().any(|&z|!finite(z)){return Err(invalid("complete finite source-modal force required"));}
            let mut rhs=vec![C64::ZERO;n];rhs[..m].copy_from_slice(f);full.push(rhs);
        }
        let refs:Vec<_>=full.iter().map(Vec::as_slice).collect();
        let results=self.prepared.forcing_residues(&refs,&mut cancelled)?;
        let mut out=Vec::with_capacity(results.len());
        let minus_i_w=C64::new(self.prepared.value().im,-self.prepared.value().re);
        for result in results{
            let mut state=result.vector;
            let memory=state.split_off(m+p);let pressure=state.split_off(m);
            let mut normal=vec![C64::ZERO;p];
            for i in 0..p{
                poll(&mut cancelled)?;let mut correction=C64::ZERO;
                for j in 0..m{add(&mut normal[i],&mut correction,
                    (minus_i_w*state[j]).scale(self.model.normal_velocity_weights()[j][i]));}
            }
            if normal.iter().any(|&z|!finite(z)){return Err(invalid("surface-velocity residue overflow"));}
            out.push(BemForceResidue{modal_displacement_residue:state,surface_pressure_residue:pressure,
                elastic_memory_residue:memory,surface_normal_velocity_residue:normal,linear_residual:result.linear_residual});
        }
        poll(&mut cancelled)?;Ok(out)
    }

    /// Directional local derivatives, with all cross-mode and material-history
    /// equations retained. No recomputation of the source modal basis is implied.
    pub fn parameter_derivatives(&self,directions:&[BemParameterDirection<'_>],mut cancelled:impl FnMut()->bool)
        ->Result<Vec<BemPoleDerivative>,BemPoleResponseError>{
        poll(&mut cancelled)?;self.admit_count(directions.len())?;
        let mut actions=Vec::with_capacity(directions.len());
        for direction in directions{actions.push(self.parameter_action(direction,&mut cancelled)?);}
        let refs:Vec<_>=actions.iter().map(Vec::as_slice).collect();
        let values=self.prepared.parameter_derivatives(&refs,&mut cancelled)?;
        let(m,p)=(self.model.mode_count(),self.model.panel_count());
        let mut out=Vec::with_capacity(values.len());
        for value in values{
            poll(&mut cancelled)?;
            let mut vector=value.vector;let memory=vector.split_off(m+p);let pressure=vector.split_off(m);
            out.push(BemPoleDerivative{angular_frequency_rad_s:value.value,modal_displacement:vector,
                surface_pressure:pressure,elastic_memory:memory,linear_residual:value.linear_residual});
        }
        poll(&mut cancelled)?;Ok(out)
    }

    fn parameter_action(&self,d:&BemParameterDirection<'_>,cancelled:&mut impl FnMut()->bool)
        ->Result<Vec<C64>,BemPoleResponseError>{
        let(m,p,r,n)=(self.model.mode_count(),self.model.panel_count(),self.model.memory_count(),self.model.dimension());
        for (slice,size) in [(d.mass,m*m),(d.stiffness_per_s2,m*m),(d.viscosity_per_s,m*m),
            (d.relaxation_rows_per_s,r*m),(d.relaxation_times_s,r)]{
            if (!slice.is_empty()&&slice.len()!=size)||slice.iter().any(|v|!v.is_finite()){
                return Err(invalid("parameter direction needs complete finite original-coordinate coefficients"));
            }
        }
        if !d.density_kg_m3.is_finite()||!d.sound_speed_m_s.is_finite(){return Err(invalid("finite medium direction required"));}
        let get=|slice:&[f64],i:usize|slice.get(i).copied().unwrap_or(0.0);
        let w=self.prepared.value();let v=self.prepared.right_vector();
        let iw=C64::new(-w.im,w.re);
        let mut out=vec![C64::ZERO;n];
        for i in 0..m{
            poll(cancelled)?;let mut correction=C64::ZERO;
            for j in 0..m{
                let coefficient=C64::from_re(get(d.stiffness_per_s2,i*m+j))
                    -(w*w).scale(get(d.mass,i*m+j))-iw.scale(get(d.viscosity_per_s,i*m+j));
                add(&mut out[i],&mut correction,coefficient*v[j]);
            }
            for j in 0..r{add(&mut out[i],&mut correction,v[m+p+j].scale(get(d.relaxation_rows_per_s,j*m+i)));}
        }
        let medium=self.model.medium();
        for i in m..m+p{
            poll(cancelled)?;let mut correction=C64::ZERO;
            for j in 0..m{
                let t=self.raw.matrix[i*n+j];let dt=self.raw.derivative[i*n+j];
                // -rho*w^2 B(w/c)W: d/dc = (2*T-w*T_w)/c.
                let coefficient=t.scale(d.density_kg_m3/medium.density)
                    +(t.scale(2.0)-w*dt).scale(d.sound_speed_m_s/medium.sound_speed);
                add(&mut out[i],&mut correction,coefficient*v[j]);
            }
            for j in m..m+p{
                let coefficient=(w*self.raw.derivative[i*n+j]).scale(-d.sound_speed_m_s/medium.sound_speed);
                add(&mut out[i],&mut correction,coefficient*v[j]);
            }
        }
        for j in 0..r{
            poll(cancelled)?;let mut correction=C64::ZERO;let at=m+p+j;
            let tau=self.model.relaxation_times_s()[j];let dtau=get(d.relaxation_times_s,j);
            for i in 0..m{
                let dc=get(d.relaxation_rows_per_s,j*m+i);let c=self.model.relaxation_rows_per_s()[j][i];
                add(&mut out[at],&mut correction,iw.scale(dtau*c+tau*dc)*v[i]);
            }
            add(&mut out[at],&mut correction,-iw.scale(dtau)*v[at]);
        }
        if out.iter().any(|&z|!finite(z)){return Err(invalid("parameter action overflow"));}
        poll(cancelled)?;Ok(out)
    }
}
