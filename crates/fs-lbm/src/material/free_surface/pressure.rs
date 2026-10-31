//! Absolute external gas pressure enters the existing single-liquid VOF step.

use super::*;
use crate::d3q19::freesurface3::FreeSurfaceError3;

impl MaterialFreeSurfaceScaling<'_> {
    /// Assign the positive absolute pressure [Pa] represented by density one,
    /// before constructing a trajectory. Without an explicit reference the
    /// original solver reports only pressure perturbations and cannot admit
    /// an absolute external-pressure field.
    ///
    /// This choice fixes the pressure gauge; it does not requery the frozen
    /// liquid's source coefficients or claim physical compressibility.
    pub fn with_pressure_reference(
        mut self,
        reference_pressure_pa: f64,
    ) -> Result<Self, MaterialLbmError> {
        positive(
            reference_pressure_pa,
            "gas pressure reference must be finite and positive",
        )?;
        self.reference_pressure_pa = Some(reference_pressure_pa);
        Ok(self)
    }

    /// Fixed absolute-pressure reference, if explicitly assigned [Pa].
    #[must_use]
    pub const fn reference_pressure_pa(&self) -> Option<f64> {
        self.reference_pressure_pa
    }

    /// Reconstruct absolute liquid pressure from the fixed declared reference
    /// and the existing incompressible pressure perturbation [Pa]. This does
    /// not turn the numerical density into a thermodynamic equation of state.
    pub fn absolute_pressure_si(&self, lattice_density: f64) -> Result<f64, MaterialLbmError> {
        let reference = self.reference_pressure_pa.ok_or(MaterialLbmError::Invalid(
            "absolute liquid pressure requires an explicit fixed pressure reference",
        ))?;
        positive(
            reference + self.pressure_perturbation_si(lattice_density)?,
            "absolute liquid pressure is nonpositive or unrepresentable",
        )
    }
}

impl MaterialFreeSurface3<'_> {
    /// Advance under a spatial external gas pressure in absolute Pascals.
    ///
    /// The scaling's declared pressure reference is the fixed physical pressure
    /// corresponding to liquid lattice density one. The callback supplies
    /// positive absolute gas pressure at exposed liquid cell centers [m], with the domain origin at
    /// zero. Both are explicit: the numerical lattice baseline is never used
    /// as an absolute gas pressure. Adding the same constant to the reference
    /// and all imposed pressures leaves the numerical trajectory unchanged.
    ///
    /// The original liquid source, viscosity, capillarity and native mass
    /// routes remain the owners. The combined gas and capillary pressure must
    /// obey the supplied fractional density-excursion bound in (0,1).
    /// Pressure may vary between calls, but one call describes a frozen
    /// boundary field over the admitted physical timestep.
    ///
    /// This is prescribed normal traction on a frozen isothermal liquid, not
    /// an evolved gas reservoir or thermodynamic pressure-work coupling. The
    /// caller owns applicability of the frozen liquid coefficients and spatial
    /// and temporal convergence. Sampling reports pressure relative to this
    /// fixed reference, including the imposed gas load and capillary jump.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite absolute pressures, unrepresentable
    /// conversion, excess numerical density excursion, or native cancellation
    /// before publishing any liquid populations, mass or conversion counters.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 3]) -> f64,
    ) -> Result<(), MaterialLbmError> {
        let reference_pressure_pa =
            self.scaling
                .reference_pressure_pa
                .ok_or(MaterialLbmError::Invalid(
                    "absolute gas loading requires an explicit fixed pressure reference",
                ))?;
        let dx = self.scaling.cell_size_m;
        let pressure_unit = self.scaling.pressure_unit_pa;
        self.lattice.step_with_gas_pressure(
            gate,
            maximum_density_excursion,
            |coordinates| {
                let position = coordinates.map(|c| (c as f64 + 0.5) * dx);
                let absolute = pressure_pa(position);
                let perturbation = (absolute - reference_pressure_pa) / pressure_unit;
                if position.iter().any(|x| !x.is_finite())
                    || !absolute.is_finite()
                    || absolute <= 0.0
                    || !perturbation.is_finite()
                {
                    return Err(FreeSurfaceError3::InvalidGasPressure {
                        coordinates: Some(coordinates),
                        reason: "absolute gas pressure must be positive with finite physical coordinates and lattice conversion",
                    });
                }
                Ok(perturbation)
            },
        )?;
        Ok(())
    }
}
