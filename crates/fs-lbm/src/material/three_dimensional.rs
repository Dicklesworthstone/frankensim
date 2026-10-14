//! The existing 3-D solvers, driven by the same frozen material as 2-D.

use fs_material::fluid_transport::FluidTransport;

use super::{MaterialLbmError, MaterialScaling, MaterialSlab};
use crate::d3q19::{PlatesGrid3, Q3, ThermalLbm3};

/// Physical boundary data for a D3Q19/D3Q7 thermal slab.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialThermalSlab3 {
    /// Fluid cells along the periodic x axis.
    pub nx: usize,
    /// Fluid cells along the periodic y axis.
    pub ny: usize,
    /// Fluid layers along z, excluding the two wall layers.
    pub nz: usize,
    /// Positive absolute bottom-wall temperature [K].
    pub bottom_temperature_k: f64,
    /// Positive absolute top-wall temperature [K].
    pub top_temperature_k: f64,
    /// Nonnegative gravity magnitude [m/s²], directed downwards.
    pub gravity_m_s2: f64,
    /// Signed volumetric expansion [1/K] at the reference material state.
    /// Must be supplied by the material model; never defaulted to air's 1/T.
    pub expansion_per_k: f64,
    /// Explicit cell budget, including both wall layers.
    pub max_cells: usize,
}

/// The existing 3-D thermal solver with the shared absolute-temperature map.
pub type MaterialSlab3 = MaterialSlab<ThermalLbm3>;

impl<F: FluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// Construct the existing D3Q19 plate channel at rest, periodic in x/y.
    /// `nz_fluid` excludes the two rigid z-wall layers. Unlike the tiled duct,
    /// this reference grid does not require dimensions divisible by four.
    /// The material reference density maps to lattice density one, just as in
    /// [`Self::flow_grid`]; all SI conversions use this same scaling object.
    ///
    /// # Errors
    /// Refuses empty/overflowing dimensions, fewer than two fluid layers and
    /// requests exceeding the explicit cell budget, before solver allocation.
    pub fn flow_plates3(
        &self,
        nx: usize,
        ny: usize,
        nz_fluid: usize,
        max_cells: usize,
    ) -> Result<PlatesGrid3, MaterialLbmError> {
        admit_plate_cells(nx, ny, nz_fluid, max_cells)?;
        Ok(PlatesGrid3::plates(nx, ny, nz_fluid, self.flow.tau))
    }

    /// Thermal relaxation for the existing D3Q7 temperature lattice.
    ///
    /// Its c_s² is 1/4, not the D2Q5 value 1/3. Both dimensions represent
    /// exactly the same physical alpha=k/(rho cp) and share dx and dt; their
    /// relaxation times must therefore differ. The scalar is a lattice
    /// property, not a material-dependent constant or an air default.
    ///
    /// # Errors
    /// Refuses overflow or a relaxation indistinguishable from one half.
    pub fn thermal_relaxation_time3(&self) -> Result<f64, MaterialLbmError> {
        let tau = 0.5 + self.thermal_diffusivity_lattice / 0.25;
        if !tau.is_finite() || tau <= 0.5 {
            return Err(MaterialLbmError::Invalid(
                "D3Q7 thermal relaxation must be representably greater than one half",
            ));
        }
        Ok(tau)
    }

    /// Construct the actual coupled 3-D flow/heat solver in physical units.
    ///
    /// Uses the same material snapshot, signed buoyancy and temperature map
    /// as [`Self::thermal_slab`], but the D3Q7 thermal relaxation. Material
    /// substitution changes momentum AND heat diffusion together. No source
    /// is re-queried and no independent fluid or thermal material is created.
    /// This remains the existing frozen-property Boussinesq model: it does not
    /// assert compressible-gas, melting, or evolving-property validity.
    ///
    /// # Errors
    /// Refuses invalid boundaries, relaxation, dimensions or allocation budget
    /// before constructing the solver (whose raw constructor may panic).
    pub fn thermal_slab3(
        &self,
        input: MaterialThermalSlab3,
    ) -> Result<MaterialSlab3, MaterialLbmError> {
        admit_plate_cells(input.nx, input.ny, input.nz, input.max_cells)?;
        let (difference, gbeta) = self.slab_boundary(
            input.bottom_temperature_k,
            input.top_temperature_k,
            input.gravity_m_s2,
            input.expansion_per_k,
        )?;
        let tau_g = self.thermal_relaxation_time3()?;
        Ok(MaterialSlab {
            lattice: ThermalLbm3::slab(
                input.nx, input.ny, input.nz, self.flow.tau, tau_g, gbeta,
            ),
            top_temperature_k: input.top_temperature_k,
            temperature_difference_k: difference,
        })
    }
}

fn admit_plate_cells(
    nx: usize,
    ny: usize,
    nz_fluid: usize,
    max_cells: usize,
) -> Result<(), MaterialLbmError> {
    let cells = nz_fluid
        .checked_add(2)
        .and_then(|nz| nx.checked_mul(ny)?.checked_mul(nz))
        .ok_or(MaterialLbmError::Invalid("3-D plate cell count overflows"))?;
    // The D3Q19 population is the largest element allocated by either solver.
    // Match Vec's representable-size bound, even for an unbounded caller budget.
    if nx == 0 || ny == 0 || nz_fluid < 2 || cells > max_cells
        || cells > isize::MAX as usize / core::mem::size_of::<[f64; Q3]>()
    {
        return Err(MaterialLbmError::Invalid(
            "3-D plates need positive x/y extents, at least two fluid layers and an admissible cell budget",
        ));
    }
    Ok(())
}
