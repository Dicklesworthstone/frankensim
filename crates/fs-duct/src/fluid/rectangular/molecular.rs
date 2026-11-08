//! Joint molecular and wall heat storage in the original rectangular section.
//!
//! The linear ideal-gas energy equation is
//! `rho cp(s) s T' - conductivity Laplacian(T') = s p'`.
//! The SAME material reservoirs supply cp(s), so its section mean is
//! `mean(T') = p' H(z)/(rho cp(s))`, z=s*rho*cp(s)/conductivity.
//! Continuity therefore uses kappa(s)=[1-R H(z)/cp(s)]/(rho R T).
//! This is the thermodynamic law of the original molecular duct owner with
//! the actual rectangular diffusion operator in place of a circular spectrum.
//! It retains isothermal storage as frequency vanishes, while the original
//! internal heat capacities equilibrate or freeze on their own time scales.
//!
//! No molecular attenuation coefficient, equilibrium thermal shunt, bulk
//! viscosity or Stokes--Kirchhoff term is added separately. Walls are rigid,
//! no-slip and isothermal; the original plane-pressure envelope, no-axial-heat
//! approximation and material model domains still apply. This evaluates a
//! frozen background, not a thermal/species evolution or transient solver.

use super::{
    FluidDuctError, PI, RectangularAdmission, RectangularDiffusion, RectangularSection,
    RectangularWave, SegmentWave, admit_tail, bad, finite, poll, positive,
};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};

impl RectangularSection {
    /// Molecular dispersion and rectangular wall diffusion from one original gas.
    /// The returned native wave and TwoPort use the same density, actual area,
    /// caloric response and thermal equation at this angular frequency.
    pub fn evaluate_molecular(
        &self,
        material: &RelaxationState,
        omega: f64,
        admission: RectangularAdmission,
    ) -> Result<RectangularWave, FluidDuctError> {
        self.evaluate_molecular_controlled(material, omega, admission, || false)
    }

    /// Allocation-free frequency query, bounded by the existing Fourier and
    /// molecular-reservoir caps. Poll every 32 transverse terms and around
    /// the original caloric evaluation. No partially prepared wave escapes.
    ///
    /// The analytic tail bounds include Re(z)/N^3 as well as Im(z)^2/N^5.
    /// In addition to the original diffusion admission, propagate those bounds
    /// through the joint compliance and admit both storage and physical loss.
    /// These exclude floating-point rounding and are not interval certificates.
    pub fn evaluate_molecular_controlled(
        &self,
        material: &RelaxationState,
        omega: f64,
        admission: RectangularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RectangularWave, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        // Preserve the exact original response for an explicitly non-relaxing
        // material. In particular there is no extra absorption or new EOS.
        if material.processes().is_empty() {
            return self.evaluate_controlled(material.gas(), omega, admission, cancelled);
        }
        if !positive(omega)
            || !positive(admission.maximum_relative_tail)
            || admission.maximum_relative_tail >= 1.0
            || !positive(admission.maximum_transverse_phase)
            || admission.maximum_transverse_phase >= PI
        {
            return Err(
                bad("rectangular molecular frequency, tail or phase envelope is invalid").into(),
            );
        }
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let acoustic = state.acoustic();
        let transport = state.transport();
        let k0 = omega / acoustic.sound_speed_m_s();
        if !positive(k0) || k0 * self.long_m > admission.maximum_transverse_phase {
            return Err(
                bad("rectangular molecular gas exceeds its plane-pressure envelope").into(),
            );
        }
        // The material law uses exp(-i omega t); the diffusion/telegraph owner
        // uses +i omega. cp_plus is its conjugate, with nonpositive imaginary
        // part, so s*cp_plus has NONNEGATIVE real and positive imaginary parts.
        let cp = material.caloric_response(omega)?.cp_j_kg_k;
        poll(&mut cancelled, 0)?;
        let viscous =
            self.diffusion(omega, transport.kinematic_viscosity_m2_s(), &mut cancelled)?;
        admit_tail(viscous, admission.maximum_relative_tail)?;
        let argument = C64::new(cp.im, cp.re)
            .scale(omega * acoustic.density_kg_m3() / transport.thermal_conductivity_w_m_k());
        let thermal = self.diffusion_argument(argument, &mut cancelled)?;
        admit_tail(thermal, admission.maximum_relative_tail)?;
        let r_over_cp = cp.conj().recip().scale(gas.specific_gas_constant);
        let factor = C64::ONE - thermal.mean * r_over_cp;
        admit_compression(factor, r_over_cp, thermal, admission.maximum_relative_tail)?;
        let compliance =
            factor.scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
        let s = C64::new(0.0, omega);
        let series = s.scale(gas.density / self.area_m2) / viscous.mean;
        let shunt = s.scale(self.area_m2) * compliance;
        let (wavenumber, characteristic_impedance) =
            crate::helmholtz_wave_from_telegraph(series, shunt);
        let specific_impedance = characteristic_impedance.scale(self.area_m2);
        let transverse_phase = k0.max(wavenumber.abs()) * self.long_m;
        if [wavenumber, characteristic_impedance, specific_impedance]
            .iter()
            .any(|v| !finite(*v))
            || wavenumber.re <= 0.0
            || wavenumber.im < 0.0
            || characteristic_impedance.re <= 0.0
            || !transverse_phase.is_finite()
            || transverse_phase > admission.maximum_transverse_phase
        {
            return Err(bad(
                "rectangular molecular wave is unresolved or outside its phase envelope",
            )
            .into());
        }
        poll(&mut cancelled, self.modes.len())?;
        Ok(RectangularWave {
            wave: SegmentWave {
                wavenumber,
                characteristic_impedance,
                specific_impedance,
                shear_number: 0.5
                    * self.short_m
                    * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
            },
            viscous,
            thermal: Some(thermal),
            transverse_phase,
            area_m2: self.area_m2,
        })
    }
}

fn admit_compression(
    factor: C64,
    r_over_cp: C64,
    thermal: RectangularDiffusion,
    tolerance: f64,
) -> Result<(), super::DuctError> {
    let error_real =
        r_over_cp.re.abs() * thermal.tail_real + r_over_cp.im.abs() * thermal.tail_imaginary;
    let error_imaginary =
        r_over_cp.re.abs() * thermal.tail_imaginary + r_over_cp.im.abs() * thermal.tail_real;
    if !finite(factor)
        || factor.re <= 0.0
        || factor.im >= 0.0
        || !error_real.is_finite()
        || !error_imaginary.is_finite()
        || error_real > tolerance * (factor.re - error_real)
        || error_imaginary > tolerance * -factor.im
    {
        return Err(bad(
            "rectangular molecular storage or loss has an unresolved Fourier tail",
        ));
    }
    Ok(())
}
