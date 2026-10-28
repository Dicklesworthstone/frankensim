//! Radial thermal diffusion in an explicitly cylindrical, compact cavity.
//!
//! The cylindrical sidewall is isothermal and both endcaps are insulated.
//! The existing shared-cavity pressure remains spatially uniform. These are
//! physical boundary conditions, not a thermal radius inferred from volume.
//! For the same Bessel law as the duct owner, C(s)/C_ad = 1 + sum g/(s+p).
//! With u_j=sqrt(g_j), inversion uses diag(p)+u*u^T. Its positive spectral
//! measure gives K(s)=K_eq + sum k_j*s/(s+lambda_j). Thus genuine first-order
//! thermal memory joins the existing Maxwell/exact owner without fake masses.

use super::*;
use fs_la::eigen::jacobi_eigh;

/// Geometry and diffusion accuracy for an isothermal cylindrical sidewall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylindricalCavityThermalOptions {
    /// Actual cavity radius [m]. Length follows from the original cavity volume.
    pub radius_m: f64,
    /// Retained radial diffusion poles, an explicit convergence axis.
    pub poles: usize,
    /// Sampled angular-frequency accuracy band [rad/s].
    pub band_rad_s: [f64; 2],
    /// Maximum sampled relative complex and real-loss Bessel error in (0,1).
    pub maximum_relative_error: f64,
}

/// One positive thermal compression relaxation in the neck-displacement chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalCompressionArm {
    /// Incremental high-frequency stiffness [N/m], acting on summed volume/area.
    pub stiffness_n_m: f64,
    /// Positive thermal relaxation time [s].
    pub relaxation_time_s: f64,
    root_stiffness: f64,
}

/// Admitted radial heat-diffusion compression law and its original coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct CylindricalThermalCompression {
    options: CylindricalCavityThermalOptions,
    length_m: f64,
    diffusion: fs_phs::ZkDiffusive,
    equilibrium_stiffness_n_m: f64,
    adiabatic_stiffness_n_m: f64,
    arms: Vec<ThermalCompressionArm>,
}

impl CylindricalThermalCompression {
    fn new(source: &HelmholtzModalCavity, options: CylindricalCavityThermalOptions,
        gate: &CancelGate) -> Result<Self, ModalCouplingError>
    {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        let gas=source.gas();
        let diffusion=fs_phs::ZkDiffusive::new(options.radius_m,gas.density,
            gas.dynamic_viscosity,gas.gamma,gas.prandtl,options.poles,
            options.band_rad_s,options.maximum_relative_error)
            .map_err(ModalCouplingError::Propagation)?;
        let length_m=source.geometry().volume_m3
            /(core::f64::consts::PI*options.radius_m*options.radius_m);
        if !length_m.is_finite() || length_m<=0.0 {
            return Err(network_invalid("cylindrical cavity thermal geometry is not representable"));
        }
        let terms=diffusion.shunt_terms();
        let excess: f64=terms.iter().map(|t|t.gain_per_s/t.rate_per_s).sum();
        let adiabatic=source.compression_stiffness_n_m();
        let equilibrium=adiabatic/(1.0+excess);
        if !excess.is_finite() || !equilibrium.is_finite() || equilibrium<=0.0 {
            return Err(network_invalid("cylindrical cavity thermal compliance is not representable"));
        }
        let n=terms.len();
        let mut arms=Vec::with_capacity(n);
        if n>0 {
            let scale=terms.iter().flat_map(|t|[t.gain_per_s,t.rate_per_s])
                .fold(0.0_f64,f64::max);
            let u:Vec<_>=terms.iter().map(|t|det::sqrt(t.gain_per_s/scale)).collect();
            let mut matrix=vec![0.0;n*n];
            for i in 0..n {
                for j in 0..n { matrix[i*n+j]=u[i]*u[j]; }
                matrix[i*n+i]+=terms[i].rate_per_s/scale;
            }
            let (values,vectors)=jacobi_eigh(&matrix,n);
            let norm=matrix.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
            let mut stiffness_fraction=0.0;
            for j in 0..n {
                if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
                let lambda=values[j];
                let projection:f64=(0..n).map(|i|u[i]*vectors[i*n+j]).sum();
                let weight=projection*projection/lambda;
                let rate=lambda*scale;
                let stiffness=adiabatic*weight;
                let tau=rate.recip();
                if ![lambda,weight,rate,stiffness,tau].iter().all(|v|v.is_finite() && *v>0.0) {
                    return Err(network_invalid("cylindrical cavity thermal spectrum is unresolved"));
                }
                let mut residual=0.0_f64;
                let mut unit=0.0;
                for i in 0..n {
                    let av:f64=(0..n).map(|k|matrix[i*n+k]*vectors[k*n+j]).sum();
                    residual=residual.max((av-lambda*vectors[i*n+j]).abs());
                    unit+=vectors[i*n+j]*vectors[i*n+j];
                }
                if residual>2e-12*norm || (unit-1.0).abs()>2e-12 {
                    return Err(network_invalid("cylindrical cavity thermal eigensolve is unresolved"));
                }
                stiffness_fraction+=weight;
                arms.push(ThermalCompressionArm {stiffness_n_m:stiffness,
                    relaxation_time_s:tau,root_stiffness:det::sqrt(stiffness)});
            }
            let target=excess/(1.0+excess);
            if (stiffness_fraction-target).abs()>2e-11*target {
                return Err(network_invalid("cylindrical cavity thermal spectrum lost static compliance"));
            }
        }
        Ok(Self {options,length_m,diffusion,equilibrium_stiffness_n_m:equilibrium,
            adiabatic_stiffness_n_m:adiabatic,arms})
    }

    /// Original physical radius and sampled approximation controls.
    pub const fn options(&self)->CylindricalCavityThermalOptions {self.options}
    /// Actual cylinder length from the original cavity volume and supplied radius.
    pub const fn length_m(&self)->f64 {self.length_m}
    /// Original Bessel coefficients and sampled error evidence.
    pub fn diffusion(&self)->&fs_phs::ZkDiffusive {&self.diffusion}
    /// Fully relaxed isothermal compression stiffness [N/m].
    pub const fn equilibrium_stiffness_n_m(&self)->f64 {self.equilibrium_stiffness_n_m}
    /// Original adiabatic compression stiffness [N/m].
    pub const fn adiabatic_stiffness_n_m(&self)->f64 {self.adiabatic_stiffness_n_m}
    /// The admitted positive thermal arms, in accepted-memory order.
    pub fn arms(&self)->&[ThermalCompressionArm] {&self.arms}

    /// Complex stiffness [N/m] in the +i omega convention, inside the admitted band.
    pub fn stiffness_at(&self,omega:f64)->Result<C64,ModalCouplingError> {
        if !omega.is_finite() || !(self.options.band_rad_s[0]..=self.options.band_rad_s[1]).contains(&omega) {
            return Err(network_invalid("cavity thermal response is outside its admitted band"));
        }
        let s=C64::new(0.0,omega);
        Ok(self.arms.iter().fold(C64::from_re(self.equilibrium_stiffness_n_m),|sum,arm|
            sum+(s/(s+C64::from_re(arm.relaxation_time_s.recip()))).scale(arm.stiffness_n_m)))
    }

    pub(super) fn observation(&self,extension:f64,memory:&[f64])->Result<(f64,f64),ModalCouplingError> {
        if memory.len()!=self.arms.len() {
            return Err(network_invalid("cavity thermal history must retain every arm"));
        }
        let mut force=self.equilibrium_stiffness_n_m*extension;
        let mut energy=0.5*self.equilibrium_stiffness_n_m*extension*extension;
        for (arm,&z) in self.arms.iter().zip(memory) {
            let mismatch=arm.root_stiffness*extension-z;
            force+=arm.root_stiffness*mismatch;
            energy+=0.5*mismatch*mismatch;
        }
        Ok((force,energy))
    }
}

impl SharedCavity {
    /// Retain radial thermal wall diffusion in an explicitly cylindrical cavity.
    ///
    /// `radius_m` is the actual cavity radius; original volume fixes its length.
    /// The sidewall is isothermal, endcaps thermally insulated, and the existing
    /// compact-cavity pressure is uniform. This is not an arbitrary-body thermal
    /// approximation, finite wall conductivity, or distributed acoustic solver.
    /// Neck viscous memory may be selected separately; exterior radiation attaches
    /// afterwards through the same existing builder. No baseline loss is doubled.
    ///
    /// `None` initializes thermal arms fully relaxed at the supplied displacement.
    /// `Some(z)` supplies energy-normalized thermal history [sqrt(J)] in arm order;
    /// all zeros represent an instantaneous adiabatic displacement from equilibrium.
    /// Every arm is a true first-order coordinate in the SAME exact propagator.
    /// Its acoustic availability and dissipation use the existing constitutive
    /// memory account, not an additional heat deposit or a second energy ledger.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_cylindrical_thermal_wall(
        source:HelmholtzModalCavity,surfaces:Vec<CavitySurface>,sample_rate_hz:u32,
        neck_initial:ModalAcousticState,neck_budget:ModalAcousticTimeBudget,
        config:ModalCouplingConfig,neck_options:Option<ViscousNeckOptions>,
        options:CylindricalCavityThermalOptions,initial_thermal_memory:Option<Vec<f64>>,
        gate:&CancelGate,
    )->Result<Self,ModalCouplingError> {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        if options.band_rad_s[1]>core::f64::consts::PI*sample_rate_hz as f64*config.nyquist_guard_fraction {
            return Err(network_invalid("cavity thermal band exceeds the original Nyquist guard"));
        }
        let thermal=CylindricalThermalCompression::new(&source,options,gate)?;
        let viscous=neck_options.map(|options|Self::prepare_viscous_neck(&source,
            sample_rate_hz,config,options,gate)).transpose()?;
        Self::new_inner(source,surfaces,sample_rate_hz,neck_initial,neck_budget,
            config,viscous,Some((thermal,initial_thermal_memory)),gate)
    }
}

#[cfg(test)]
mod tests;
