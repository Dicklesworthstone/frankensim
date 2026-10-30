//! Source-owned thermal expansion for the existing Boussinesq solvers.

use fs_material::fluid_transport::{
    FluidTransport, ThermalFluidTransport, ThermalFluidTransportState,
};

use super::{
    MaterialLbmError, MaterialScaling, MaterialSlab, MaterialSlab3, MaterialThermalSlab,
    MaterialThermalSlab3,
};

/// Physical walls and a declared limit on the frozen-density approximation.
///
/// Expansion and reference temperature come from the SAME material that
/// supplies viscosity, conductivity, density and heat capacity. The limit
/// checks `abs(beta * (T - T_material))`; it is a necessary local Boussinesq
/// check, not an EOS support, phase-stability or solver-stability certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialBuoyancyBoundary {
    /// Absolute bottom-wall temperature [K].
    pub bottom_temperature_k: f64,
    /// Absolute top-wall temperature [K].
    pub top_temperature_k: f64,
    /// Nonnegative gravitational acceleration [m/s²], directed downwards.
    pub gravity_m_s2: f64,
    /// Positive fractional density-departure budget, strictly below one.
    pub max_relative_density_change: f64,
}

impl<F: ThermalFluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// Construct 2-D coupled flow and heat with material-owned buoyancy.
    ///
    /// `fluid_cells` excludes the two wall rows. No separately authored
    /// expansion coefficient can become stale after a material substitution.
    /// The initial linear temperature profile satisfies the declared density
    /// budget; subsequent forcing must keep the evolving field in that regime.
    ///
    /// # Errors
    /// Refuses a stale/incoherent material snapshot, unknown expansion,
    /// excessive temperature departures and the existing lattice admission errors.
    pub fn thermal_slab_from_material(
        &self,
        fluid_cells: [usize; 2],
        boundary: MaterialBuoyancyBoundary,
        max_cells: usize,
    ) -> Result<MaterialSlab, MaterialLbmError> {
        let (thermal, beta) = self.admit_material_buoyancy(boundary)?;
        let reference = normalized_reference(boundary, thermal.temperature_k())?;
        let mut slab = self.thermal_slab(MaterialThermalSlab {
            nx: fluid_cells[0],
            ny: fluid_cells[1],
            bottom_temperature_k: boundary.bottom_temperature_k,
            top_temperature_k: boundary.top_temperature_k,
            gravity_m_s2: boundary.gravity_m_s2,
            expansion_per_k: beta,
            max_cells,
        })?;
        slab.lattice.t_ref = reference;
        Ok(slab)
    }

    /// Three-dimensional counterpart using the same thermal material state.
    /// `fluid_cells[2]` excludes the two z-wall layers. D3Q7 heat and D3Q19
    /// momentum retain their original, distinct lattice diffusivity conversions.
    pub fn thermal_slab3_from_material(
        &self,
        fluid_cells: [usize; 3],
        boundary: MaterialBuoyancyBoundary,
        max_cells: usize,
    ) -> Result<MaterialSlab3, MaterialLbmError> {
        let (thermal, beta) = self.admit_material_buoyancy(boundary)?;
        let reference = normalized_reference(boundary, thermal.temperature_k())?;
        let mut slab = self.thermal_slab3(MaterialThermalSlab3 {
            nx: fluid_cells[0],
            ny: fluid_cells[1],
            nz: fluid_cells[2],
            bottom_temperature_k: boundary.bottom_temperature_k,
            top_temperature_k: boundary.top_temperature_k,
            gravity_m_s2: boundary.gravity_m_s2,
            expansion_per_k: beta,
            max_cells,
        })?;
        slab.lattice.t_ref = reference;
        Ok(slab)
    }

    /// Refresh live 2-D transport AND buoyancy from one new material state.
    ///
    /// Every current cell and both walls must satisfy the new state's density
    /// budget. The existing conservative rebind preserves populations and
    /// refuses density, heat-capacity or lattice-unit changes. All admission
    /// precedes mutation; a refused update leaves the previous solver intact.
    pub fn rebind_thermal_slab_from_material<G: FluidTransport + ?Sized>(
        &self,
        slab: &mut MaterialSlab,
        previous: &MaterialScaling<'_, G>,
        gravity_m_s2: f64,
        max_relative_density_change: f64,
    ) -> Result<(), MaterialLbmError> {
        let boundary = MaterialBuoyancyBoundary {
            bottom_temperature_k: slab.top_temperature_k + slab.temperature_difference_k,
            top_temperature_k: slab.top_temperature_k,
            gravity_m_s2,
            max_relative_density_change,
        };
        let (thermal, beta) = self.admit_material_buoyancy(boundary)?;
        let reference = normalized_reference(boundary, thermal.temperature_k())?;
        for population in &slab.lattice.gpop {
            admit_temperature(
                slab.temperature_si(population.iter().sum())?,
                thermal.temperature_k(),
                beta,
                max_relative_density_change,
            )?;
        }
        self.rebind_thermal_slab(slab, previous, gravity_m_s2, beta)?;
        slab.lattice.t_ref = reference;
        Ok(())
    }

    /// Three-dimensional conservative refresh with the same source and field checks.
    pub fn rebind_thermal_slab3_from_material<G: FluidTransport + ?Sized>(
        &self,
        slab: &mut MaterialSlab3,
        previous: &MaterialScaling<'_, G>,
        gravity_m_s2: f64,
        max_relative_density_change: f64,
    ) -> Result<(), MaterialLbmError> {
        let boundary = MaterialBuoyancyBoundary {
            bottom_temperature_k: slab.top_temperature_k + slab.temperature_difference_k,
            top_temperature_k: slab.top_temperature_k,
            gravity_m_s2,
            max_relative_density_change,
        };
        let (thermal, beta) = self.admit_material_buoyancy(boundary)?;
        let reference = normalized_reference(boundary, thermal.temperature_k())?;
        for population in &slab.lattice.gpop {
            admit_temperature(
                slab.temperature_si(population.iter().sum())?,
                thermal.temperature_k(),
                beta,
                max_relative_density_change,
            )?;
        }
        self.rebind_thermal_slab3(slab, previous, gravity_m_s2, beta)?;
        slab.lattice.t_ref = reference;
        Ok(())
    }

    fn admit_material_buoyancy(
        &self,
        boundary: MaterialBuoyancyBoundary,
    ) -> Result<(ThermalFluidTransportState, f64), MaterialLbmError> {
        if !boundary.gravity_m_s2.is_finite()
            || boundary.gravity_m_s2 < 0.0
            || !boundary.max_relative_density_change.is_finite()
            || !(0.0..1.0).contains(&boundary.max_relative_density_change)
            || boundary.max_relative_density_change == 0.0
        {
            return Err(MaterialLbmError::Invalid(
                "material buoyancy needs finite nonnegative gravity and a density budget in (0,1)",
            ));
        }
        let thermal = self.material.thermal_transport_state()?;
        if thermal.transport() != self.transport {
            return Err(MaterialLbmError::Invalid(
                "thermal material must match the momentum/heat transport snapshot; rebuild the scaling",
            ));
        }
        let beta = thermal
            .volumetric_expansion_per_k()
            .ok_or(MaterialLbmError::Invalid(
                "material buoyancy requires known isobaric volumetric expansion",
            ))?;
        for temperature in [boundary.bottom_temperature_k, boundary.top_temperature_k] {
            admit_temperature(
                temperature,
                thermal.temperature_k(),
                beta,
                boundary.max_relative_density_change,
            )?;
        }
        Ok((thermal, beta))
    }
}

fn normalized_reference(
    boundary: MaterialBuoyancyBoundary,
    temperature_k: f64,
) -> Result<f64, MaterialLbmError> {
    let difference = boundary.bottom_temperature_k - boundary.top_temperature_k;
    let reference = (temperature_k - boundary.top_temperature_k) / difference;
    if !difference.is_finite() || difference == 0.0 || !reference.is_finite() {
        return Err(MaterialLbmError::Invalid(
            "material reference temperature must be representable in the slab's temperature mapping",
        ));
    }
    Ok(reference)
}

fn admit_temperature(
    temperature_k: f64,
    reference_k: f64,
    beta: f64,
    limit: f64,
) -> Result<(), MaterialLbmError> {
    let departure = (beta * (temperature_k - reference_k)).abs();
    if !temperature_k.is_finite()
        || temperature_k <= 0.0
        || !departure.is_finite()
        || departure > limit
    {
        return Err(MaterialLbmError::Invalid(
            "temperature exceeds the material's declared Boussinesq density-departure budget",
        ));
    }
    Ok(())
}
