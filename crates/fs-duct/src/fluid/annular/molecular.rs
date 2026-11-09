//! Molecular relaxation and both annular wall heat layers in one energy equation.
//!
//! In +i omega bookkeeping, the original reservoirs supply cp(s) and
//! `rho s cp(s) T' - conductivity Laplacian(T') = s p'`.
//! Thus `mean(T') = p' H(z)/(rho cp(s))`, z=s*rho*cp(s)/conductivity,
//! and the original gas equation of state gives the section compliance
//! `[1 - R H(z)/cp(s)]/(rho R T)`. The SAME weighted radial P1 mass,
//! stiffness, no-slip/isothermal walls and complex tridiagonal solve remain
//! in use. A second equilibrium thermal shunt, bulk viscosity or separately
//! added molecular/Stokes--Kirchhoff attenuation would double-count storage.
//!
//! This is the complex-caloric specialization of Stinson's arbitrary-section
//! energy/continuity equations, JASA 89 (1991), 550-558, equations 57-69,
//! <https://doi.org/10.1121/1.400379>, using the original material's Debye
//! reservoirs. It has the same fully developed, linear, continuum and
//! transverse-pressure assumptions as the constant-Cp annular owner.
//! Molecular rates and heat capacities retain their source validity;
//! neither radial refinement nor algebraic work residuals establish that
//! validity or bound entrance, rarefaction, axial-heat or mean-flow error.

use super::field::PreparedAnnularFields;
use super::{
    AnnularAdmission, AnnularSection, AnnularWave, FluidDuctError, RadialSolution, SegmentWave,
    bad, finite, poll, positive,
};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};

struct MolecularSolution {
    state: ThermoviscousFluidState,
    cp: C64,
    wave: AnnularWave,
    viscous: RadialSolution,
    thermal: RadialSolution,
}

impl AnnularSection {
    /// Joint wall diffusion and molecular caloric response of one original gas.
    /// No equivalent circular radius or additive attenuation replaces the
    /// actual annular geometry, local pressure-dependent rates or heat storage.
    pub fn evaluate_molecular(
        &self,
        material: &RelaxationState,
        omega: f64,
        admission: AnnularAdmission,
    ) -> Result<AnnularWave, FluidDuctError> {
        self.evaluate_molecular_controlled(material, omega, admission, || false)
    }

    /// Each active field uses one native O(N) radial solve under the original
    /// payload and cell-phase budgets. Poll before preparation, around each
    /// solve and every 64 cells. No partial response is published on refusal.
    ///
    /// The thermal `gradient_work` remains the independently integrated radial
    /// identity Im(H), with scale Im(z)/|z|². It does not include the molecular
    /// relaxation work, which is already present in the complex compliance.
    pub fn evaluate_molecular_controlled(
        &self,
        material: &RelaxationState,
        omega: f64,
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AnnularWave, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        if material.processes().is_empty() {
            return self.evaluate_controlled(material.gas(), omega, admission, cancelled);
        }
        Ok(self
            .solve_molecular(material, omega, admission, &mut cancelled)?
            .wave)
    }

    /// Retain the original complex-Cp radial fields, source and average wave.
    /// Fourier and molecular availability losses remain separate, with the
    /// original coherent finite-segment excitation available on the result.
    pub fn prepare_molecular_fields<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        admission: AnnularAdmission,
    ) -> Result<PreparedAnnularFields<'g, 'm, RelaxationState>, FluidDuctError> {
        self.prepare_molecular_fields_controlled(material, omega, admission, || false)
    }

    /// The same two native O(N) solves supply fields and propagation. All
    /// original reservoirs and both isothermal walls participate in one energy
    /// equation; no summed attenuation or equilibrium thermal shunt is added.
    /// Subsequent point, ring and finite-interval queries allocate nothing and
    /// neither reevaluate the material nor solve another radial problem.
    pub fn prepare_molecular_fields_controlled<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAnnularFields<'g, 'm, RelaxationState>, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        if material.processes().is_empty() {
            return self
                .prepare_fields_controlled(material.gas(), omega, admission, cancelled)
                .map(|fields| fields.with_source(material));
        }
        let solved = self.solve_molecular(material, omega, admission, &mut cancelled)?;
        Ok(PreparedAnnularFields::from_solutions(
            self,
            material,
            solved.state,
            omega,
            solved.wave,
            solved.viscous,
            Some(solved.thermal),
            Some((solved.cp, material.gas().temperature)),
            &mut cancelled,
        )?)
    }

    fn solve_molecular(
        &self,
        material: &RelaxationState,
        omega: f64,
        admission: AnnularAdmission,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularSolution, FluidDuctError> {
        self.admit(omega, admission)?;
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let acoustic = state.acoustic();
        let transport = state.transport();
        let k0 = omega / acoustic.sound_speed_m_s();
        if !positive(k0) || k0 * self.outer_m > admission.maximum_transverse_phase {
            return Err(bad("annular molecular gas exceeds its plane-pressure envelope").into());
        }
        // Original calorics use exp(-i omega t). Conjugate into the diffusion
        // convention so s*cp(s) has nonnegative real and positive imaginary parts.
        let cp = material.caloric_response(omega)?.cp_j_kg_k;
        poll(cancelled, 0)?;
        let viscous = self.diffusion_solution(
            omega,
            transport.kinematic_viscosity_m2_s(),
            admission,
            cancelled,
        )?;
        let argument = C64::new(cp.im, cp.re).scale(
            (omega * acoustic.density_kg_m3() / transport.thermal_conductivity_w_m_k())
                * self.gap_m
                * self.gap_m,
        );
        let magnitude = argument.abs();
        let thermal = self.diffusion_solution_scaled(
            argument,
            det::sqrt(magnitude),
            det::sqrt(argument.im) / magnitude,
            admission,
            cancelled,
        )?;
        let factor =
            C64::ONE - thermal.response.mean * cp.conj().recip().scale(gas.specific_gas_constant);
        let compliance =
            factor.scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
        if !finite(compliance) || compliance.re <= 0.0 || compliance.im >= 0.0 {
            return Err(
                bad("annular molecular compression must retain positive storage and loss").into(),
            );
        }
        let s = C64::new(0.0, omega);
        let series = s.scale(gas.density / self.area_m2) / viscous.response.mean;
        let shunt = s.scale(self.area_m2) * compliance;
        let (wavenumber, characteristic_impedance) =
            crate::helmholtz_wave_from_telegraph(series, shunt);
        let specific_impedance = characteristic_impedance.scale(self.area_m2);
        let transverse_phase = k0.max(wavenumber.abs()) * self.outer_m;
        if ![wavenumber, characteristic_impedance, specific_impedance]
            .into_iter()
            .all(finite)
            || wavenumber.re <= 0.0
            || wavenumber.im < 0.0
            || characteristic_impedance.re <= 0.0
            || !transverse_phase.is_finite()
            || transverse_phase > admission.maximum_transverse_phase
        {
            return Err(
                bad("annular molecular wave is unresolved or outside its phase envelope").into(),
            );
        }
        poll(cancelled, self.elements())?;
        let wave = AnnularWave {
            wave: SegmentWave {
                wavenumber,
                characteristic_impedance,
                specific_impedance,
                shear_number: 0.5
                    * self.gap_m
                    * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
            },
            viscous: viscous.response,
            thermal: Some(thermal.response),
            transverse_phase,
            area_m2: self.area_m2,
        };
        Ok(MolecularSolution {
            state,
            cp,
            wave,
            viscous,
            thermal,
        })
    }
}
