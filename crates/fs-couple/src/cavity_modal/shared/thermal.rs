//! Radial thermal diffusion in an explicitly cylindrical, compact cavity.
//!
//! The cylindrical sidewall is isothermal and both endcaps are insulated.
//! The existing shared-cavity pressure remains spatially uniform. These are
//! physical boundary conditions, not a thermal radius inferred from volume.
//! For the same Bessel law as the duct owner, C(s)/C_ad = 1 + sum g/(s+p).
//! With u_j=sqrt(g_j), inversion uses diag(p)+u*u^T. Its positive spectral
//! measure gives K(s)=K_eq + sum k_j*s/(s+lambda_j). Thus genuine first-order
//! thermal memory joins the existing Maxwell/exact owner without fake masses.
//!
//! The molecular constructor couples the same diffusion coordinates to the
//! original material heat reservoirs BEFORE forming compression arms. It retains
//! the isothermal static limit and the frozen molecular instantaneous limit in
//! one positive storage/heat account, used by both time and harmonic response.

use super::*;
use fs_la::eigen::jacobi_eigh;
use fs_material::relaxation::RelaxationState;

mod molecular;
mod spectrum;
pub mod material_wall;
pub use material_wall::{CylindricalMaterialWallCompression, CylindricalMaterialWallOptions};

// A single runtime compression owner, preserving each public source-specific
// accessor. Material walls and prescribed-temperature walls never run in parallel.
pub(super) enum Compression {
    Isothermal(CylindricalThermalCompression),
    Material(CylindricalMaterialWallCompression),
}
impl Compression {
    pub(super) fn band_rad_s(&self)->[f64;2] { match self {
        Self::Isothermal(law)=>law.options().band_rad_s,
        Self::Material(law)=>law.options().band_rad_s,
    }}
    pub(super) fn equilibrium_stiffness_n_m(&self)->f64 { match self {
        Self::Isothermal(law)=>law.equilibrium_stiffness_n_m(),
        Self::Material(law)=>law.equilibrium_stiffness_n_m(),
    }}
    pub(super) fn arms(&self)->&[ThermalCompressionArm] { match self {
        Self::Isothermal(law)=>law.arms(), Self::Material(law)=>law.arms(),
    }}
    pub(super) fn stiffness_at(&self,omega:f64)->Result<C64,ModalCouplingError> { match self {
        Self::Isothermal(law)=>law.stiffness_at(omega),
        Self::Material(law)=>law.stiffness_at(omega),
    }}
    pub(super) fn observation(&self,extension:f64,memory:&[f64])->Result<(f64,f64),ModalCouplingError> { match self {
        Self::Isothermal(law)=>law.observation(extension,memory),
        Self::Material(law)=>law.observation(extension,memory),
    }}
}

fn compression_stiffness(equilibrium:f64,arms:&[ThermalCompressionArm],band:[f64;2],omega:f64)
    ->Result<C64,ModalCouplingError>
{
    if !omega.is_finite() || !(band[0]..=band[1]).contains(&omega) {
        return Err(network_invalid("cavity thermal response is outside its admitted band"));
    }
    let s=C64::new(0.0,omega);
    Ok(arms.iter().fold(C64::from_re(equilibrium),|sum,arm|
        sum+(s/(s+C64::from_re(arm.relaxation_time_s.recip()))).scale(arm.stiffness_n_m)))
}

fn compression_observation(equilibrium:f64,arms:&[ThermalCompressionArm],extension:f64,memory:&[f64])
    ->Result<(f64,f64),ModalCouplingError>
{
    if memory.len()!=arms.len() {
        return Err(network_invalid("cavity thermal history must retain every arm"));
    }
    let mut force=equilibrium*extension;
    let mut energy=0.5*equilibrium*extension*extension;
    for (arm,&z) in arms.iter().zip(memory) {
        let mismatch=arm.root_stiffness*extension-z;
        force+=arm.root_stiffness*mismatch;
        energy+=0.5*mismatch*mismatch;
    }
    Ok((force,energy))
}

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
    instantaneous_stiffness_n_m: f64,
    arms: Vec<ThermalCompressionArm>,
    molecular_relaxation: Option<RelaxationState>,
}

impl CylindricalThermalCompression {
    fn new(source: &HelmholtzModalCavity, options: CylindricalCavityThermalOptions,
        gate: &CancelGate) -> Result<Self, ModalCouplingError>
    {
        Self::new_inner(source, options, None, gate)
    }

    /// Couple the original molecular reservoirs to radial wall heat diffusion.
    ///
    /// Both mechanisms act on the same gas temperatures. Independent molecular
    /// and thermal springs would double-count storage and miss their coupled
    /// transition. The source gas must match the complete relaxation snapshot.
    /// The isothermal limit remains p A_neck^2/V; the instantaneous limit uses
    /// the molecularly frozen heat capacity.
    ///
    /// Preparation uses one joint symmetric eigensolve, bounded by
    /// maximum_coordinates (at most 128). Its dimension is the thermal pole
    /// count times one plus the number of nonzero molecular reservoirs.
    /// Every resulting arm joins the original exact mechanical propagator.
    /// Algebraic accuracy of complex stiffness AND its loss is checked at
    /// 33 frequencies to 1e-8 relative error inside the supplied band.
    /// The separate equilibrium Bessel admission is not a bound on the spatial
    /// error of the joint molecular model; thermal pole refinement remains
    /// necessary. Local Debye reservoirs and constant Fourier conductivity are
    /// retained, without separate diffusion of molecular internal energy.
    pub fn new_with_molecular_relaxation(
        source: &HelmholtzModalCavity,
        options: CylindricalCavityThermalOptions,
        relaxation: RelaxationState,
        maximum_coordinates: usize,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        Self::new_inner(source, options, Some((relaxation, maximum_coordinates)), gate)
    }

    fn new_inner(source: &HelmholtzModalCavity, options: CylindricalCavityThermalOptions,
        molecular: Option<(RelaxationState, usize)>,
        gate: &CancelGate) -> Result<Self, ModalCouplingError>
    {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        let gas=source.gas();
        if let Some((law, _)) = &molecular {
            if law.gas() != gas {
                return Err(network_invalid("cavity thermal and molecular owners need the same gas snapshot"));
            }
        }
        let diffusion=fs_phs::ZkDiffusive::new(options.radius_m,gas.density,
            gas.dynamic_viscosity,gas.gamma,gas.prandtl,options.poles,
            options.band_rad_s,options.maximum_relative_error)
            .map_err(ModalCouplingError::Propagation)?;
        let length_m=source.geometry().volume_m3
            /(core::f64::consts::PI*options.radius_m*options.radius_m);
        if !length_m.is_finite() || length_m<=0.0 {
            return Err(network_invalid("cylindrical cavity thermal geometry is not representable"));
        }
        if let Some((law, maximum_coordinates)) = molecular {
            let modes: Vec<_> = diffusion.shunt_terms().iter().map(|term|
                molecular::ThermalMode {
                    weight: term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0)),
                    rate_per_s: term.rate_per_s,
                }).collect();
            let spectrum = molecular::realize(&law, &modes, options.band_rad_s,
                maximum_coordinates, &mut || gate.is_requested()).map_err(|reason| {
                    if gate.is_requested() { ModalCouplingError::Cancelled }
                    else { network_invalid(reason) }
                })?;
            let area = source.neck_area_m2();
            let scale = area * area / source.geometry().volume_m3;
            let equilibrium = gas.pressure * scale;
            let instantaneous = spectrum.frozen_bulk_modulus_pa * scale;
            if ![scale,equilibrium,instantaneous].iter().all(|v|v.is_finite() && *v>0.0) {
                return Err(network_invalid("molecular cavity compression geometry is unrepresentable"));
            }
            let mut arms = Vec::with_capacity(spectrum.arms.len());
            for arm in spectrum.arms {
                let stiffness = arm.increment_pa * scale;
                if !stiffness.is_finite() || stiffness<=0.0 {
                    return Err(network_invalid("molecular cavity compression arm is unrepresentable"));
                }
                arms.push(ThermalCompressionArm {stiffness_n_m:stiffness,
                    relaxation_time_s:arm.relaxation_time_s,
                    root_stiffness:det::sqrt(stiffness)});
            }
            return Ok(Self {options,length_m,diffusion,
                equilibrium_stiffness_n_m:equilibrium,
                adiabatic_stiffness_n_m:source.compression_stiffness_n_m(),
                instantaneous_stiffness_n_m:instantaneous,
                arms,molecular_relaxation:Some(law)});
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
            adiabatic_stiffness_n_m:adiabatic,instantaneous_stiffness_n_m:adiabatic,
            arms,molecular_relaxation:None})
    }

    /// Original physical radius and sampled approximation controls.
    pub const fn options(&self)->CylindricalCavityThermalOptions {self.options}
    /// Actual cylinder length from the original cavity volume and supplied radius.
    pub const fn length_m(&self)->f64 {self.length_m}
    /// Original Bessel coefficients and sampled error evidence.
    pub fn diffusion(&self)->&fs_phs::ZkDiffusive {&self.diffusion}
    /// Fully relaxed isothermal compression stiffness [N/m].
    pub const fn equilibrium_stiffness_n_m(&self)->f64 {self.equilibrium_stiffness_n_m}
    /// Original equilibrium-gas adiabatic stiffness [N/m], with molecular
    /// reservoirs equilibrated. Use instantaneous_stiffness_n_m for the frozen limit.
    pub const fn adiabatic_stiffness_n_m(&self)->f64 {self.adiabatic_stiffness_n_m}
    /// Instantaneous compression stiffness [N/m], including frozen reservoirs.
    pub const fn instantaneous_stiffness_n_m(&self)->f64 {self.instantaneous_stiffness_n_m}
    /// Original material law, when molecular relaxation is attached.
    pub fn molecular_relaxation(&self)->Option<&RelaxationState> {self.molecular_relaxation.as_ref()}
    /// The admitted positive thermal arms, in accepted-memory order.
    pub fn arms(&self)->&[ThermalCompressionArm] {&self.arms}

    /// Complex stiffness [N/m] in the +i omega convention, inside the admitted band.
    pub fn stiffness_at(&self,omega:f64)->Result<C64,ModalCouplingError> {
        compression_stiffness(self.equilibrium_stiffness_n_m,&self.arms,self.options.band_rad_s,omega)
    }

    pub(super) fn observation(&self,extension:f64,memory:&[f64])->Result<(f64,f64),ModalCouplingError> {
        compression_observation(self.equilibrium_stiffness_n_m,&self.arms,extension,memory)
    }
}

impl SharedCavity {
    /// Couple molecular relaxation and a cylindrical thermal sidewall to the
    /// SAME moving surfaces, neck, exterior loading and accepted time history.
    ///
    /// The original molecular source and wall diffusion form joint positive
    /// compression arms. None initializes them fully relaxed; supplied history
    /// must contain every admitted joint arm in the returned spectrum's order.
    /// Neck viscous memory is an independent option. Constant conductivity,
    /// isothermal sidewall, insulated endcaps and uniform acoustic pressure are
    /// the same physical restrictions as the thermal constructor below.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_cylindrical_molecular_wall(
        source:HelmholtzModalCavity,surfaces:Vec<CavitySurface>,sample_rate_hz:u32,
        neck_initial:ModalAcousticState,neck_budget:ModalAcousticTimeBudget,
        config:ModalCouplingConfig,neck_options:Option<ViscousNeckOptions>,
        options:CylindricalCavityThermalOptions,relaxation:RelaxationState,
        maximum_coordinates:usize,initial_thermal_memory:Option<Vec<f64>>,
        gate:&CancelGate,
    )->Result<Self,ModalCouplingError> {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        if options.band_rad_s[1]>core::f64::consts::PI*sample_rate_hz as f64*config.nyquist_guard_fraction {
            return Err(network_invalid("cavity thermal band exceeds the original Nyquist guard"));
        }
        let thermal=CylindricalThermalCompression::new_with_molecular_relaxation(
            &source,options,relaxation,maximum_coordinates,gate)?;
        let viscous=neck_options.map(|options|Self::prepare_viscous_neck(&source,
            sample_rate_hz,config,options,gate)).transpose()?;
        Self::new_inner(source,surfaces,sample_rate_hz,neck_initial,neck_budget,
            config,viscous,Some((Compression::Isothermal(thermal),initial_thermal_memory)),None,gate)
    }

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
            config,viscous,Some((Compression::Isothermal(thermal),initial_thermal_memory)),None,gate)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod molecular_tests;
