//! Material-driven physical units for the existing flow and thermal solvers.
//!
//! One frozen transport snapshot and one (dx, dt) pair supply both relaxation
//! times. Lattice density is relative to this material's reference density.
//! This is an incompressible/Boussinesq initialization, not a compressible gas
//! EOS, phase transition, or conservative replacement of an evolving material.
//! A successful plan checks necessary numerical bounds, not solver stability.

pub mod conjugate;

use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

mod forcing;
mod buoyancy;
pub use buoyancy::MaterialBuoyancyBoundary;
pub mod free_surface;
pub use forcing::{electrothermal, transport};
mod runtime_rebind;
mod three_dimensional;
pub use three_dimensional::{MaterialSlab3, MaterialThermalSlab3};

use crate::{CS2, Grid, MACH_LIMIT, ScalingPlan, ThermalLbm};

/// A material, physical-unit, numerical-range or allocation-budget refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialLbmError {
    /// The original material transport owner refused its snapshot.
    Transport(MaterialError),
    /// An external thermal-flow force field was refused before stepping.
    Forcing(crate::thermal::ThermalForcingError),
    /// The existing sparse free-surface owner refused construction or evolution.
    FreeSurface(crate::d3q19::freesurface3::FreeSurfaceError3),
    /// Invalid or unrepresentable physical/numerical input.
    Invalid(&'static str),
}
impl core::fmt::Display for MaterialLbmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Transport(error) => write!(f, "material LBM: {error}"),
            Self::Forcing(error) => write!(f, "material LBM: {error}"),
            Self::FreeSurface(error) => write!(f, "material LBM: {error}"),
            Self::Invalid(message) => write!(f, "material LBM: {message}"),
        }
    }
}
impl std::error::Error for MaterialLbmError {}
impl From<crate::thermal::ThermalForcingError> for MaterialLbmError {
    fn from(error: crate::thermal::ThermalForcingError) -> Self { Self::Forcing(error) }
}
impl From<MaterialError> for MaterialLbmError {
    fn from(error: MaterialError) -> Self { Self::Transport(error) }
}
impl From<crate::d3q19::freesurface3::FreeSurfaceError3> for MaterialLbmError {
    fn from(error: crate::d3q19::freesurface3::FreeSurfaceError3) -> Self { Self::FreeSurface(error) }
}

/// A frozen source binding and its common momentum/thermal lattice conversion.
///
/// Source receipts remain on `material()`. This numerical adapter does not
/// turn those receipts into an experimental accuracy or stability certificate.
pub struct MaterialScaling<'a, F: FluidTransport + ?Sized> {
    material: &'a F,
    transport: FluidTransportState,
    cell_size_m: f64,
    timestep_s: f64,
    velocity_scale: f64,
    flow: ScalingPlan,
    thermal_diffusivity_lattice: f64,
    thermal_relaxation_time: f64,
}

/// Physical boundary data for the existing normalized-temperature slab.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialThermalSlab {
    /// Fluid cells horizontally.
    pub nx: usize,
    /// Fluid rows vertically, excluding the two wall rows.
    pub ny: usize,
    /// Positive absolute bottom-wall temperature [K].
    pub bottom_temperature_k: f64,
    /// Positive absolute top-wall temperature [K].
    pub top_temperature_k: f64,
    /// Nonnegative gravity magnitude [m/s²], directed downwards.
    pub gravity_m_s2: f64,
    /// Signed volumetric expansion [1/K] at the material reference state.
    /// Negative expansion is meaningful; it is never replaced with 1/T.
    pub expansion_per_k: f64,
    /// Explicit allocation budget, including the two wall rows.
    pub max_cells: usize,
}

/// Compile a gas or liquid evaluation into the existing LBM scaling contract.
///
/// `max_speed_m_s` is an explicit scenario bound. It must be maintained by the
/// caller during evolution; this constructor cannot predict future velocities.
/// The lattice Mach gate does NOT establish physical gas incompressibility or
/// equate numerical lattice sound speed to the material's acoustic sound speed.
///
/// # Errors
/// Refuses invalid transport, dx/dt/speed, unrepresentable conversions,
/// relaxation times indistinguishable from 1/2, and excessive lattice Mach.
pub fn plan_scaling_from_material<F: FluidTransport + ?Sized>(
    material: &F,
    cell_size_m: f64,
    timestep_s: f64,
    max_speed_m_s: f64,
) -> Result<MaterialScaling<'_, F>, MaterialLbmError> {
    let transport = material.transport_state()?;
    let (flow, velocity_scale) = plan_momentum_scaling(
        transport.kinematic_viscosity_m2_s(), cell_size_m, timestep_s, max_speed_m_s,
    )?;
    let diffusivity_scale = positive((timestep_s / cell_size_m) / cell_size_m, "diffusivity conversion is unrepresentable")?;
    let thermal_diffusivity_lattice = positive(transport.thermal_diffusivity_m2_s() * diffusivity_scale, "lattice thermal diffusivity is unrepresentable")?;
    let thermal_relaxation_time = relaxation(thermal_diffusivity_lattice)?;
    Ok(MaterialScaling {
        material, transport, cell_size_m, timestep_s, velocity_scale,
        flow,
        thermal_diffusivity_lattice, thermal_relaxation_time,
    })
}

impl<F: FluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// Original material owner, including receipts for source-resolved inputs.
    #[must_use]
    pub const fn material(&self) -> &F { self.material }
    /// The single snapshot used by both momentum and heat populations.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Existing momentum scaling representation and its necessary bounds.
    #[must_use]
    pub const fn flow(&self) -> &ScalingPlan { &self.flow }
    /// D2Q5 thermal relaxation time from that same snapshot and physical scaling.
    #[must_use]
    pub const fn thermal_relaxation_time(&self) -> f64 { self.thermal_relaxation_time }
    /// Thermal diffusivity in lattice units.
    #[must_use]
    pub const fn thermal_diffusivity_lattice(&self) -> f64 { self.thermal_diffusivity_lattice }
    /// Physical grid spacing [m].
    #[must_use]
    pub const fn cell_size_m(&self) -> f64 { self.cell_size_m }
    /// Physical integration step [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 { self.timestep_s }

    /// Build the actual D2Q9 grid at rest; no second flow solver is introduced.
    /// Its numerical reference density is one. Retain this scaling for output.
    ///
    /// # Errors
    /// Refuses empty/overflowing dimensions or an exceeded allocation budget.
    pub fn flow_grid(&self, nx: usize, ny: usize, max_cells: usize) -> Result<Grid, MaterialLbmError> {
        admit_cells(nx, ny, max_cells)?;
        Ok(Grid::uniform(nx, ny, self.flow.tau))
    }

    /// Build the actual D2Q9/D2Q5 conducting slab and its temperature mapping.
    ///
    /// The existing solver uses theta=(T-T_top)/(T_bottom-T_top), with walls
    /// at 1 and 0. Hence its force coefficient is g beta DeltaT dt²/dx. Signed
    /// DeltaT and expansion are retained: heating from above or a negative
    /// expansion coefficient reverses buoyancy instead of taking an absolute.
    /// The caller must establish the small-density-variation Boussinesq regime.
    ///
    /// # Errors
    /// Refuses invalid wall temperatures, zero contrast, gravity, expansion,
    /// numerical conversions and allocation requests before creating the grid.
    pub fn thermal_slab(&self, input: MaterialThermalSlab) -> Result<MaterialSlab, MaterialLbmError> {
        if input.ny < 2 || !input.gravity_m_s2.is_finite() || input.gravity_m_s2 < 0.0
            || !input.expansion_per_k.is_finite()
        {
            return Err(MaterialLbmError::Invalid("slab needs at least two fluid rows, finite expansion and nonnegative gravity"));
        }
        let (difference, gbeta) = self.slab_boundary(
            input.bottom_temperature_k, input.top_temperature_k,
            input.gravity_m_s2, input.expansion_per_k,
        )?;
        let rows = input.ny.checked_add(2).ok_or(MaterialLbmError::Invalid("slab row count overflows"))?;
        admit_cells(input.nx, rows, input.max_cells)?;
        Ok(MaterialSlab {
            lattice: ThermalLbm::slab(input.nx, input.ny, self.flow.tau, self.thermal_relaxation_time, gbeta),
            top_temperature_k: input.top_temperature_k,
            temperature_difference_k: difference,
            force_units: forcing::ForceUnits::from_scale(self),
        })
    }

    // Shared by the 2-D and 3-D thermal solvers. Neither dimension may silently
    // change the sign, absolute-temperature convention or material expansion.
    fn slab_boundary(
        &self,
        bottom_temperature_k: f64,
        top_temperature_k: f64,
        gravity_m_s2: f64,
        expansion_per_k: f64,
    ) -> Result<(f64, f64), MaterialLbmError> {
        positive(bottom_temperature_k, "bottom absolute temperature must be positive and finite")?;
        positive(top_temperature_k, "top absolute temperature must be positive and finite")?;
        if !gravity_m_s2.is_finite() || gravity_m_s2 < 0.0 || !expansion_per_k.is_finite() {
            return Err(MaterialLbmError::Invalid("slab needs finite expansion and nonnegative finite gravity"));
        }
        let difference = bottom_temperature_k - top_temperature_k;
        if !difference.is_finite() || difference == 0.0 {
            return Err(MaterialLbmError::Invalid("normalized slab requires a finite nonzero wall temperature contrast"));
        }
        let gbeta = self.acceleration_lattice(gravity_m_s2 * expansion_per_k * difference)?;
        Ok((difference, gbeta))
    }

    /// Numerical reference-one density to physical density [kg/m³].
    pub fn density_si(&self, lattice_density: f64) -> Result<f64, MaterialLbmError> {
        positive(lattice_density, "lattice density must be finite and positive")?;
        positive(lattice_density * self.transport.density_kg_m3(), "physical density is unrepresentable")
    }

    /// A signed velocity component, converted to lattice units.
    /// The full vector's speed must still respect the admitted scenario bound.
    pub fn velocity_lattice(&self, velocity_m_s: f64) -> Result<f64, MaterialLbmError> {
        finite(velocity_m_s, "physical velocity must be finite")?;
        finite(velocity_m_s / self.velocity_scale, "lattice velocity is unrepresentable")
    }

    /// A signed lattice velocity component to m/s.
    pub fn velocity_si(&self, lattice_velocity: f64) -> Result<f64, MaterialLbmError> {
        finite(lattice_velocity, "lattice velocity must be finite")?;
        finite(lattice_velocity * self.velocity_scale, "physical velocity is unrepresentable")
    }

    /// Signed physical acceleration to lattice acceleration, NOT force density.
    /// A Guo force density additionally multiplies by numerical cell density.
    pub fn acceleration_lattice(&self, acceleration_m_s2: f64) -> Result<f64, MaterialLbmError> {
        finite(acceleration_m_s2, "physical acceleration must be finite")?;
        finite((acceleration_m_s2 * self.timestep_s) / self.velocity_scale, "lattice acceleration is unrepresentable")
    }

    /// Incompressible pressure perturbation [Pa] about the material reference.
    /// This is not absolute gas pressure or a physical acoustic equation of state.
    pub fn pressure_perturbation_si(&self, lattice_density: f64) -> Result<f64, MaterialLbmError> {
        positive(lattice_density, "lattice density must be finite and positive")?;
        finite(CS2 * (lattice_density - 1.0) * self.transport.density_kg_m3()
            * self.velocity_scale * self.velocity_scale, "physical pressure perturbation is unrepresentable")
    }
}

/// Existing thermal solver together with its dimension-independent temperature mapping.
/// Defaults to the D2Q9/D2Q5 solver; [`MaterialSlab3`] uses D3Q19/D3Q7.
/// Numerical fields remain public through `lattice`; manually changing its
/// normalized wall values invalidates this frozen physical-unit mapping.
pub struct MaterialSlab<L = ThermalLbm> {
    /// Actual solver; step it through `lattice.step()`.
    pub lattice: L,
    top_temperature_k: f64,
    temperature_difference_k: f64,
    force_units: forcing::ForceUnits,
}
impl<L> MaterialSlab<L> {
    /// Convert any finite normalized temperature to absolute K.
    /// This also works for negative bottom-minus-top temperature contrasts.
    pub fn temperature_si(&self, normalized: f64) -> Result<f64, MaterialLbmError> {
        finite(normalized, "normalized temperature must be finite")?;
        positive(self.top_temperature_k + normalized *self.temperature_difference_k,
            "physical temperature must remain finite and positive")
    }
    /// Physical temperature scale used by both wall mapping and buoyancy.
    #[must_use]
    pub const fn temperature_difference_k(&self) -> f64 { self.temperature_difference_k }
}

// Both caloric and isothermal liquid consumers use exactly the same numerical
// momentum contract. Heat capacity/conductivity are not momentum prerequisites.
fn plan_momentum_scaling(
    kinematic_viscosity_m2_s: f64,
    cell_size_m: f64,
    timestep_s: f64,
    max_speed_m_s: f64,
) -> Result<(ScalingPlan, f64), MaterialLbmError> {
    positive(cell_size_m, "cell size must be finite and positive")?;
    positive(timestep_s, "timestep must be finite and positive")?;
    positive(kinematic_viscosity_m2_s, "kinematic viscosity must be finite and positive")?;
    if !max_speed_m_s.is_finite() || max_speed_m_s < 0.0 {
        return Err(MaterialLbmError::Invalid("maximum speed must be finite and nonnegative"));
    }
    let velocity_scale = positive(cell_size_m / timestep_s, "velocity unit is unrepresentable")?;
    let diffusivity_scale = positive((timestep_s / cell_size_m) / cell_size_m, "diffusivity conversion is unrepresentable")?;
    let viscosity = positive(kinematic_viscosity_m2_s * diffusivity_scale, "lattice viscosity is unrepresentable")?;
    let tau = relaxation(viscosity)?;
    let u_lattice = max_speed_m_s / velocity_scale;
    let mach = u_lattice / CS2.sqrt();
    if !mach.is_finite() || mach >= MACH_LIMIT {
        return Err(MaterialLbmError::Invalid("scenario speed exceeds the existing lattice Mach limit"));
    }
    Ok((ScalingPlan { tau, viscosity, u_lattice, mach, tau_margin: tau - 0.5, stable: true }, velocity_scale))
}

fn relaxation(diffusivity: f64) -> Result<f64, MaterialLbmError> {
    let tau = 0.5 + diffusivity / CS2;
    if !tau.is_finite() || tau <= 0.5 {
        Err(MaterialLbmError::Invalid("relaxation must be representably greater than one half"))
    } else { Ok(tau) }
}
fn admit_cells(nx: usize, ny: usize, max_cells: usize) -> Result<(), MaterialLbmError> {
    let cells = nx.checked_mul(ny).ok_or(MaterialLbmError::Invalid("grid cell count overflows"))?;
    // Bound the largest per-cell vector before Vec allocation, even when the
    // caller supplies usize::MAX as its budget. Allocation itself may still OOM.
    if nx == 0 || ny == 0 || cells > max_cells
        || cells > isize::MAX as usize / core::mem::size_of::<[f64; crate::Q]>()
    {
        return Err(MaterialLbmError::Invalid("grid is empty, oversized or exceeds its explicit cell budget"));
    }
    Ok(())
}
fn positive(value: f64, message: &'static str) -> Result<f64, MaterialLbmError> {
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(MaterialLbmError::Invalid(message)) }
}
fn finite(value: f64, message: &'static str) -> Result<f64, MaterialLbmError> {
    if value.is_finite() { Ok(value) } else { Err(MaterialLbmError::Invalid(message)) }
}
