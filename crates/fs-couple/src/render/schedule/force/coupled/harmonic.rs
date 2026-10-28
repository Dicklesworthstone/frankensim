//! Continuous forced response of the original coupled linear system.
//! Peak phasors use exp(-i omega t); averages therefore carry a factor 1/2.
use super::*;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, det};

/// Dense response capacity; no retained mode or acoustic pole is truncated.
pub const MAX_HARMONIC_UNKNOWNS: usize = 1024;

/// AC particular solution, independent of accepted histories and DC preload.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalHarmonicResponse {
    pub omega_rad_s: f64,
    /// Original flattened modal coordinates, including any internal components.
    pub displacement_m_sqrt_kg: Vec<C64>,
    pub velocity_m_sqrt_kg_per_s: Vec<C64>,
    /// Energy-normalized Maxwell internal strains, in original branch order.
    pub material_memory: Vec<C64>,
    /// Acoustic displacement/momentum pairs, in original radiation-pole order.
    pub radiation_state: Vec<C64>,
    /// AC force on each connection's positive side, in connection order [N].
    pub connection_forces_n: Vec<C64>,
    /// Mean quadratic storage of the AC solution only [J].
    pub mean_stored_energy_j: f64,
    /// Mean work rate of the supplied peak-force phasors [W].
    pub external_power_w: f64,
    pub component_dissipation_w: f64,
    pub connection_dissipation_w: f64,
    pub material_dissipation_w: f64,
    pub radiation_dissipation_w: f64,
    /// Total physical mean loss minus external mean power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the ORIGINAL equations.
    pub solve_relative_residual: f64,
}
impl ModalHarmonicResponse {
    pub fn total_dissipation_w(&self) -> f64 {
        self.component_dissipation_w + self.connection_dissipation_w
            + self.material_dissipation_w + self.radiation_dissipation_w
    }
}

fn complex_finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }
fn lu_error(error: LuComplexError) -> ModalCouplingError {
    match error {
        LuComplexError::Cancelled => ModalCouplingError::Cancelled,
        LuComplexError::Singular { .. } => invalid("harmonic dynamic stiffness is singular"),
    }
}
fn projection(column: &[f64], q: &[C64]) -> C64 {
    column.iter().zip(q).fold(C64::ZERO, |sum,(&b,&x)| sum+x.scale(b))
}
fn add_rank(matrix: &mut [C64], dimension: usize, column: &[f64], gain: C64,
    gate: &CancelGate) -> Result<(),ModalCouplingError>
{
    for (i,&a) in column.iter().enumerate() {
        poll(Some(gate))?;
        for (j,&b) in column.iter().enumerate() {
            matrix[i*dimension+j] = matrix[i*dimension+j] + gain.scale(a*b);
        }
    }
    Ok(())
}

impl CoupledModalSystem {
    /// Continuous linear response to peak generalized-force phasors in
    /// exp(-i omega t), with ALL original connections, Maxwell arms and air.
    /// This is the AC particular solution: initial transients, DC rest lengths,
    /// nonlinear contacts and the discrete held-force transfer are not included.
    /// Material/fit frequency validity remains with the coefficient producer.
    ///
    /// Maxwell states are eliminated algebraically, without approximating their
    /// frequency dependence. Acoustic displacement unknowns remain explicit, so
    /// a bare lossless air pole can produce its physical antiresonance rather
    /// than a spurious refusal from dividing by an isolated-pole impedance.
    ///
    /// `maximum_work` bounds d^3+d^2*(connections+material_arms+3), with
    /// d=mechanical_modes+acoustic_poles, before allocation; d is also capped at
    /// MAX_HARMONIC_UNKNOWNS. The existing solve-relative tolerance gates the
    /// original equations. Constant-loss prepared connections require their own
    /// basis law and explicitly refuse; they are never treated as dashpots here.
    /// No accepted state, material history, sample clock or ledger is changed.
    pub fn harmonic_response(&self, omega_rad_s: f64, force: &[C64], maximum_work: usize,
        gate: &CancelGate) -> Result<ModalHarmonicResponse,ModalCouplingError>
    {
        poll(Some(gate))?;
        let n=self.mode_count();
        let poles=self.radiation.as_ref().map_or(0,|r|r.model.poles.len());
        let dimension=n.checked_add(poles).ok_or_else(||invalid("harmonic dimension overflow"))?;
        let square=dimension.checked_mul(dimension).ok_or_else(||invalid("harmonic matrix overflow"))?;
        let work=dimension.checked_add(self.connections.len())
            .and_then(|v|v.checked_add(self.relaxation_branches.len()))
            .and_then(|v|v.checked_add(3)).and_then(|v|square.checked_mul(v))
            .ok_or_else(||invalid("harmonic work overflow"))?;
        if dimension>MAX_HARMONIC_UNKNOWNS || work>maximum_work {
            return Err(invalid("complete harmonic solve exceeds its explicit dense work budget"));
        }
        let omega_sq=omega_rad_s*omega_rad_s;
        if !omega_rad_s.is_finite() || omega_rad_s<=0.0 || !omega_sq.is_finite() || omega_sq==0.0
            || force.len()!=n || force.iter().any(|&f|!complex_finite(f)) {
            return Err(invalid("harmonic response requires positive representable frequency and one finite force per mode"));
        }
        if self.connections.iter().any(|c|c.constant_loss_reference_omega_rad_s.is_some()) {
            return Err(invalid("harmonic constant-loss prepared connections are not Kelvin-Voigt"));
        }
        let s=C64::new(0.0,-omega_rad_s);
        let modes:Vec<_>=self.models.iter().flat_map(|m|m.modes()).collect();
        let mut matrix=vec![C64::ZERO;square];
        for (i,mode) in modes.iter().enumerate() {
            poll(Some(gate))?;
            let w=mode.angular_frequency_rad_s;
            matrix[i*dimension+i]=C64::new(w*w-omega_sq,-omega_rad_s*(2.0*mode.damping_ratio*w));
        }
        for (column,law) in self.columns.iter().zip(&self.connections) {
            add_rank(&mut matrix,dimension,column,C64::new(law.stiffness_n_m,
                -omega_rad_s*law.damping_n_s_m),gate)?;
        }
        let mut column=vec![0.0;n];
        for arm in &self.relaxation_branches {
            let root=det::sqrt(arm.stiffness);
            for (b,&p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) { *b=root*p; }
            let gain=s/(s+C64::from_re(arm.relaxation_time_s.recip()));
            add_rank(&mut matrix,dimension,&column,gain,gate)?;
        }
        if let Some(air)=&self.radiation {
            for (j,pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let index=n+j;
                matrix[index*dimension+index]=C64::new(pole.omega*pole.omega-omega_sq,
                    -omega_rad_s*(2.0*pole.zeta*pole.omega));
                for (&i,&g) in air.source_modes.iter().zip(&pole.coupling) {
                    matrix[i*dimension+index]=s.scale(g);
                    matrix[index*dimension+i]=-s.scale(g);
                }
            }
        }
        if matrix.iter().any(|&a|!complex_finite(a)) { return Err(invalid("harmonic operator overflow")); }
        let scale=matrix.iter().map(|a|a.abs()).fold(0.0_f64,f64::max);
        if !scale.is_finite() || scale==0.0 { return Err(invalid("harmonic operator has no finite scale")); }
        let scaled:Vec<_>=matrix.iter().map(|a|C64::new(a.re/scale,a.im/scale)).collect();
        let factor=lu_complex_with_cancel(&scaled,dimension,||gate.is_requested()).map_err(lu_error)?;
        let mut state=vec![C64::ZERO;dimension];
        for (x,f) in state.iter_mut().zip(force) { *x=C64::new(f.re/scale,f.im/scale); }
        if state.iter().any(|&x|!complex_finite(x)) { return Err(invalid("harmonic force scaling overflow")); }
        factor.solve_with_cancel(&mut state,||gate.is_requested()).map_err(lu_error)?;
        if state.iter().any(|&x|!complex_finite(x)) { return Err(invalid("harmonic response is not representable")); }
        let mut residual_norm=0.0_f64;let mut operator_norm=0.0_f64;
        let state_norm=state.iter().map(|x|x.abs()).fold(0.0_f64,f64::max);
        let force_norm=force.iter().map(|x|x.abs()).fold(0.0_f64,f64::max);
        for i in 0..dimension {
            poll(Some(gate))?;
            let rhs=force.get(i).copied().unwrap_or(C64::ZERO);
            let mut residual=-rhs;let mut row_norm=0.0;
            for (j,&x) in state.iter().enumerate() {
                residual=residual+matrix[i*dimension+j]*x;
                row_norm+=matrix[i*dimension+j].abs();
            }
            if !complex_finite(residual) || !row_norm.is_finite() { return Err(invalid("harmonic residual overflow")); }
            residual_norm=residual_norm.max(residual.abs());operator_norm=operator_norm.max(row_norm);
        }
        let denominator=operator_norm*state_norm+force_norm;
        if !denominator.is_finite() { return Err(invalid("harmonic residual scale overflow")); }
        let relative=if denominator==0.0 {0.0} else {residual_norm/denominator};
        if relative>self.config.solve_relative_tolerance {
            return Err(ModalCouplingError::SolveResidual {relative,tolerance:self.config.solve_relative_tolerance});
        }
        let mut result=ModalHarmonicResponse {omega_rad_s,
            displacement_m_sqrt_kg:state[..n].to_vec(),
            velocity_m_sqrt_kg_per_s:state[..n].iter().map(|&q|s*q).collect(),
            material_memory:Vec::with_capacity(self.relaxation_branches.len()),
            radiation_state:Vec::with_capacity(2*poles),connection_forces_n:Vec::with_capacity(self.connections.len()),
            mean_stored_energy_j:0.0,external_power_w:0.0,component_dissipation_w:0.0,
            connection_dissipation_w:0.0,material_dissipation_w:0.0,radiation_dissipation_w:0.0,
            power_residual_w:0.0,solve_relative_residual:relative};
        for (i,mode) in modes.iter().enumerate() {
            let q=state[i];let v=result.velocity_m_sqrt_kg_per_s[i];let w=mode.angular_frequency_rad_s;
            result.mean_stored_energy_j+=0.25*(q.scale(w).norm_sq()+v.norm_sq());
            result.component_dissipation_w+=mode.damping_ratio*w*v.norm_sq();
            result.external_power_w+=0.5*(force[i]*v.conj()).re;
        }
        for (column,law) in self.columns.iter().zip(&self.connections) {
            poll(Some(gate))?;
            let extension=projection(column,&state[..n]);let velocity=s*extension;
            result.connection_forces_n.push(-(extension.scale(law.stiffness_n_m)+velocity.scale(law.damping_n_s_m)));
            result.mean_stored_energy_j+=0.25*law.stiffness_n_m*extension.norm_sq();
            result.connection_dissipation_w+=0.5*law.damping_n_s_m*velocity.norm_sq();
        }
        for arm in &self.relaxation_branches {
            poll(Some(gate))?;
            let root=det::sqrt(arm.stiffness);let rate=arm.relaxation_time_s.recip();
            for (b,&p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) { *b=root*p; }
            let strain=projection(&column,&state[..n]);let denominator=s+C64::from_re(rate);
            let memory=(C64::from_re(rate)/denominator)*strain;
            let mismatch=(s/denominator)*strain;
            result.material_memory.push(memory);
            result.mean_stored_energy_j+=0.25*mismatch.norm_sq();
            result.material_dissipation_w+=0.5*rate*mismatch.norm_sq();
        }
        if let Some(air)=&self.radiation {
            for (j,pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let q=state[n+j];let p=s*q;
                result.radiation_state.extend([q,p]);
                result.mean_stored_energy_j+=0.25*(q.scale(pole.omega).norm_sq()+p.norm_sq());
                result.radiation_dissipation_w+=pole.zeta*pole.omega*p.norm_sq();
            }
        }
        result.power_residual_w=result.total_dissipation_w()-result.external_power_w;
        if result.velocity_m_sqrt_kg_per_s.iter().chain(&result.material_memory)
            .chain(&result.radiation_state).chain(&result.connection_forces_n).any(|&z|!complex_finite(z))
            || [result.mean_stored_energy_j,result.external_power_w,result.total_dissipation_w(),
                result.power_residual_w].iter().any(|v|!v.is_finite()) {
            return Err(invalid("harmonic physical power or storage overflow"));
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
