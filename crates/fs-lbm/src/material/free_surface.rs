//! Material-driven physical units for the existing sparse 3-D free surface.
//!
//! This is the frozen, isothermal, single-liquid rung: the existing VOF owner
//! evolves momentum, fill and its conservative mass ledger. The complete source
//! point supplies rho, mu (or nu), and optionally uniform sigma; no caloric
//! property is fabricated. Geometry, timestep, gravity and the numerical wall
//! ghost are explicit scenario data. No thermal, moving-wall, gas-flow, wetting
//! calibration or conservative surface-energy evolution is introduced here.
//! The `enthalpy` owner adds thermal evolution; its `rheology` path alone may
//! consume source bundles that explicitly require evolving physical shear rate.

use fs_exec::CancelGate;
use fs_material::liquid::surface::{ResolvedLiquidSurfaceState, ResolvedLiquidWettingState};
use fs_material::state_point::MaterialStatePointError;

use super::{MaterialLbmError, finite, plan_momentum_scaling, positive};
use crate::d3q19::freesurface3::{Cell3, ContactModel3, FreeSurface3, PlanarWall3};
use crate::d3q19::sparse::MORTON_COORD_BITS;
use crate::d3q19::{Q3, TILE};
use crate::{CS2, ScalingPlan};

mod initial;
pub use initial::MaterialFreeSurfaceInitialCell3;
/// Geometric, moment-preserving ingress of an already-liquid tetrahedral body.
pub mod transfer;
/// Conservative moving-liquid enthalpy and source-temperature viscosity.
pub mod enthalpy;
/// Source-pinned hysteresis on the actual planar VOF contact line.
pub mod wetting;

/// One reference source state and its physical-to-lattice conversion.
///
/// A necessary numerical admission, not a stability or accuracy certificate.
/// The caller must maintain its maximum-speed bound and establish adequate
/// interface resolution and small density excursions during the actual run.
#[derive(Clone, Copy, Debug)]
pub struct MaterialFreeSurfaceScaling<'a> {
    material: &'a ResolvedLiquidSurfaceState,
    flow: ScalingPlan,
    cell_size_m: f64,
    timestep_s: f64,
    velocity_scale: f64,
    pressure_unit_pa: f64,
    mass_unit_kg: f64,
    surface_tension_lattice: f64,
    planar_wetting: Option<(PlanarWall3, &'a ResolvedLiquidWettingState)>,
}

/// Stationary voxel domain and explicitly chosen numerical wall model.
#[derive(Clone, Copy, Debug)]
pub struct MaterialFreeSurfaceDomain3<'a> {
    /// Global cell counts, each a positive multiple of the native tile edge.
    pub cells: [usize; 3],
    /// Gravity [m/s²]. Native forcing is the constant reference-density force
    /// rho0*g, not a compressible density-dependent force law.
    pub gravity_m_s2: [f64; 3],
    /// Existing neutral/wetting fill ghost, explicitly chosen by the scenario.
    /// Neither choice infers a measured contact angle from material names.
    pub contact: ContactModel3,
    /// Permanently solid tile coordinates, in the native 4-cell tile convention.
    pub solid_tiles: &'a [(u32, u32, u32)],
    /// Bound on the complete domain's cells, including currently inactive space.
    /// This bounds possible later activation, not just the initial wet region.
    pub max_cells: usize,
}

/// A wet-cell sample in the admitted physical units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceSample3 {
    /// Native fluid/interface classification.
    pub cell: Cell3,
    /// Native liquid fill fraction.
    pub fill_fraction: f64,
    /// Numerical liquid density mapped about the sourced reference [kg/m³].
    pub density_kg_m3: f64,
    /// Guo-corrected liquid velocity [m/s].
    pub velocity_m_s: [f64; 3],
    /// Incompressible pressure perturbation relative to the imposed gas-side
    /// reference [Pa], not absolute pressure or a physical acoustic EOS.
    pub pressure_perturbation_pa: f64,
}

/// The actual sparse VOF solver with its original material and physical scale.
///
/// Step through [`MaterialFreeSurface3::step`]; the original mass/conversion accounting and
/// error semantics are retained. Source coefficients and geometry units cannot
/// drift through this wrapper. A different source point starts a different
/// frozen scenario; it is not an energy-accounted runtime replacement.
pub struct MaterialFreeSurface3<'a> {
    /// Existing solver, with no alternate collision or free-surface integrator.
    lattice: FreeSurface3,
    scaling: MaterialFreeSurfaceScaling<'a>,
}

/// Compile a liquid source bundle's reference momentum/capillary scaling.
/// Frozen constructors refuse a source with required rheology evolution;
/// the dedicated enthalpy/rheology constructor resolves it at every flow state.
///
/// sigma_lattice = sigma_SI / (rho0 * (dx/dt)^2 * dx).
/// The native bounded-curvature pressure reconstruction must remain strictly
/// positive at either sign of its curvature limit. This prevents a known
/// inadmissible boundary density; it does not establish low compressibility,
/// capillary-wave convergence or a contact-angle law.
///
/// # Errors
/// Refuses the shared viscosity/relaxation/Mach/dx/dt gates, unrepresentable SI
/// units, and a surface coefficient exceeding the native reconstruction domain.
pub fn plan_free_surface_from_material(
    material: &ResolvedLiquidSurfaceState,
    cell_size_m: f64,
    timestep_s: f64,
    max_speed_m_s: f64,
) -> Result<MaterialFreeSurfaceScaling<'_>, MaterialLbmError> {
    let (flow, velocity_scale) = plan_momentum_scaling(
        material.kinematic_viscosity_m2_s(), cell_size_m, timestep_s, max_speed_m_s,
    )?;
    let pressure_unit_pa = positive(
        (material.density_kg_m3() * velocity_scale) * velocity_scale,
        "free-surface pressure unit is unrepresentable",
    )?;
    let mass_unit_kg = positive(
        ((material.density_kg_m3() * cell_size_m) * cell_size_m) * cell_size_m,
        "free-surface cell mass unit is unrepresentable",
    )?;
    let surface_tension_lattice = match material.surface_tension_n_m() {
        None => 0.0,
        Some(sigma) => positive(
            (sigma / pressure_unit_pa) / cell_size_m,
            "lattice surface tension is unrepresentable",
        )?,
    };
    let density_excursion = surface_tension_lattice * FreeSurface3::CURVATURE_LIMIT / CS2;
    if !density_excursion.is_finite() || density_excursion >= 1.0 {
        return Err(MaterialLbmError::Invalid(
            "surface tension and scaling can produce nonpositive native capillary boundary density",
        ));
    }
    Ok(MaterialFreeSurfaceScaling {
        material, flow, cell_size_m, timestep_s, velocity_scale,
        pressure_unit_pa, mass_unit_kg, surface_tension_lattice,
        planar_wetting: None,
    })
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    fn require_frozen_source(&self) -> Result<(), MaterialLbmError> {
        if self.material.rheology_axes().is_some() {
            return Err(MaterialLbmError::Invalid(
                "shear-rate material requires the evolving free-surface rheology constructor",
            ));
        }
        Ok(())
    }

    /// Original source state, property receipts and interface environment.
    #[must_use]
    pub const fn material(&self) -> &'a ResolvedLiquidSurfaceState { self.material }

    /// The same momentum scaling representation as the thermal/flow adapters.
    #[must_use]
    pub const fn flow(&self) -> &ScalingPlan { &self.flow }

    /// Physical grid spacing [m].
    #[must_use]
    pub const fn cell_size_m(&self) -> f64 { self.cell_size_m }

    /// Physical integration step [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 { self.timestep_s }

    /// Native uniform surface-tension coefficient; zero means explicitly omitted.
    #[must_use]
    pub const fn surface_tension_lattice(&self) -> f64 { self.surface_tension_lattice }

    /// Worst capillary reference-density excursion from the native curvature
    /// clamp. This is a numerical bound, not a predicted actual interface value.
    #[must_use]
    pub fn capillary_density_excursion_bound(&self) -> f64 {
        self.surface_tension_lattice * FreeSurface3::CURVATURE_LIMIT / CS2
    }

    /// Mass [kg] represented by one lattice mass unit, rho0*dx³.
    #[must_use]
    pub const fn mass_unit_kg(&self) -> f64 { self.mass_unit_kg }

    /// Signed velocity component from lattice units to m/s.
    pub fn velocity_si(&self, lattice_velocity: f64) -> Result<f64, MaterialLbmError> {
        finite(lattice_velocity, "free-surface lattice velocity must be finite")?;
        finite(lattice_velocity * self.velocity_scale, "free-surface velocity is unrepresentable")
    }

    /// Constant rho0*g forcing in native reference-density-one units.
    /// The actual solver's velocity includes the Guo half-force correction.
    pub fn gravity_lattice(&self, gravity_m_s2: f64) -> Result<f64, MaterialLbmError> {
        finite(gravity_m_s2, "free-surface gravity must be finite")?;
        finite((gravity_m_s2 * self.timestep_s) / self.velocity_scale,
            "free-surface lattice gravity is unrepresentable")
    }

    /// Incompressible pressure perturbation [Pa] about the reference-one density.
    pub fn pressure_perturbation_si(&self, lattice_density: f64) -> Result<f64, MaterialLbmError> {
        positive(lattice_density, "free-surface lattice density must be finite and positive")?;
        finite(CS2 * (lattice_density - 1.0) * self.pressure_unit_pa,
            "free-surface pressure perturbation is unrepresentable")
    }

    /// Assign the sourced equilibrium angle to one planar domain wall before
    /// constructing any flow state. The angle is measured through the liquid;
    /// other domain walls and internal solid tiles retain the neutral ghost.
    ///
    /// The caller assigns the substrate and named environment to this geometry.
    /// The retained interface does not validate an ambient gas EOS or claim a
    /// calibrated moving contact line. Sigma remains owned by this SAME liquid
    /// surface source; the wetting interface contributes only the angle.
    ///
    /// # Errors
    /// Refuses a different liquid phase/process/revision, complete query point
    /// or quantity descriptors before evaluating geometry or allocating flow.
    pub fn with_planar_wetting(
        mut self,
        wall: PlanarWall3,
        wetting: &'a ResolvedLiquidWettingState,
    ) -> Result<Self, MaterialStatePointError> {
        wetting.require_liquid_state(self.material)?;
        self.planar_wetting = Some((wall, wetting));
        Ok(self)
    }

    /// Explicit geometry assignment and source-bearing equilibrium-angle owner.
    #[must_use]
    pub const fn planar_wetting(&self) -> Option<(PlanarWall3, &'a ResolvedLiquidWettingState)> {
        self.planar_wetting
    }

    /// Construct the existing sparse solver after admitting the complete domain.
    ///
    /// The classifier is a fixed predicate on global cell coordinates, as for
    /// the native constructor. Physical fill follows its half-filled initial
    /// interface cells; binary input occupancy is not an exact voxel mass.
    ///
    /// # Errors
    /// Refuses invalid/oversized dimensions, out-of-domain or contradictory solid
    /// tiles, a legacy Wetting ghost combined with a sourced planar angle,
    /// invalid force conversion and native free-surface construction.
    pub fn free_surface3(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        fluid: impl Fn(usize, usize, usize) -> bool,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        self.require_frozen_source()?;
        if self.planar_wetting.is_some() && domain.contact != ContactModel3::Neutral {
            return Err(MaterialLbmError::Invalid(
                "a sourced planar angle requires Neutral on all other walls",
            ));
        }
        admit_domain(&domain, &fluid)?;
        let mut gravity = [0.0; 3];
        for (component, physical) in gravity.iter_mut().zip(domain.gravity_m_s2) {
            *component = self.gravity_lattice(physical)?;
        }
        let [nx, ny, nz] = domain.cells;
        let mut lattice = FreeSurface3::with_solid_tiles(
            nx, ny, nz, self.flow.tau, gravity, self.surface_tension_lattice,
            domain.contact, domain.solid_tiles, fluid,
        )?;
        if let Some((wall, wetting)) = self.planar_wetting {
            lattice.set_planar_contact_angle(wall, wetting.contact_angle_radians())?;
        }
        Ok(MaterialFreeSurface3 { lattice, scaling: self })
    }
}

impl<'a> MaterialFreeSurface3<'a> {
    /// Read-only access to the actual sparse VOF owner.
    #[must_use]
    pub const fn lattice(&self) -> &FreeSurface3 { &self.lattice }

    /// Advance the existing native owner by the admitted physical timestep.
    /// Its original refusal and mass-conversion semantics are unchanged.
    ///
    /// # Errors
    /// Propagates the native collision/grid refusal; this wrapper does not add
    /// transactional guarantees beyond those of the native step.
    pub fn step(&mut self) -> Result<(), MaterialLbmError> {
        self.lattice.step()?;
        Ok(())
    }

    /// Advance with the native cancellation polls before physical publication.
    /// An observed cancellation preserves populations, mass and conversions;
    /// the derived fill cache is refreshed on a retry.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), MaterialLbmError> {
        self.lattice.step_with_gate(gate)?;
        Ok(())
    }

    /// Retire fully gas tiles through the native owner; the liquid mass ledger
    /// is unchanged. Returns its retired tile count.
    pub fn retire_gas_tiles(&mut self) -> Result<usize, MaterialLbmError> {
        Ok(self.lattice.retire_gas_tiles()?)
    }

    /// Frozen source and physical unit owner used to construct this trajectory.
    #[must_use]
    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'a> { &self.scaling }

    /// Conservative native ledger in kg, including its unregioned carry.
    pub fn ledger_mass_kg(&self) -> Result<f64, MaterialLbmError> {
        finite(self.lattice.ledger_mass() * self.scaling.mass_unit_kg,
            "free-surface mass ledger is unrepresentable")
    }

    /// Region's native fluid/interface mass in kg; excludes the global carry.
    pub fn region_mass_kg(
        &self, region: impl Fn(usize, usize, usize) -> bool,
    ) -> Result<f64, MaterialLbmError> {
        finite(self.lattice.region_mass(region) * self.scaling.mass_unit_kg,
            "free-surface regional mass is unrepresentable")
    }

    /// Signed unregioned conservation carry in kg, not an extra fluid cell.
    pub fn carry_mass_kg(&self) -> Result<f64, MaterialLbmError> {
        finite(self.lattice.carry() * self.scaling.mass_unit_kg,
            "free-surface conservation carry is unrepresentable")
    }

    /// Sample a wet cell in physical units. Gas/wall/out-of-domain stays None.
    ///
    /// # Errors
    /// Refuses inadmissible or unrepresentable numerical output. Sampling does
    /// not certify the state at other cells or the entire run's Mach bound.
    pub fn sample_at_si(
        &self, x: usize, y: usize, z: usize,
    ) -> Result<Option<MaterialFreeSurfaceSample3>, MaterialLbmError> {
        let Some(sample) = self.lattice.sample_at(x, y, z) else { return Ok(None); };
        let density_kg_m3 = positive(
            sample.density * self.scaling.material.density_kg_m3(),
            "free-surface physical density must be finite and positive",
        )?;
        let mut velocity_m_s = [0.0; 3];
        for (component, lattice) in velocity_m_s.iter_mut().zip(sample.velocity) {
            *component = self.scaling.velocity_si(lattice)?;
        }
        Ok(Some(MaterialFreeSurfaceSample3 {
            cell: sample.cell,
            fill_fraction: finite(sample.fill_fraction, "free-surface fill must be finite")?,
            density_kg_m3, velocity_m_s,
            pressure_perturbation_pa: self.scaling.pressure_perturbation_si(sample.density)?,
        }))
    }
}

fn admit_domain(
    domain: &MaterialFreeSurfaceDomain3<'_>,
    fluid: &impl Fn(usize, usize, usize) -> bool,
) -> Result<(), MaterialLbmError> {
    let cells = domain.cells.into_iter().try_fold(1usize, |count, n| count.checked_mul(n))
        .ok_or(MaterialLbmError::Invalid("free-surface domain cell count overflows"))?;
    if domain.cells.into_iter().any(|n| n == 0 || n % TILE != 0
        || n / TILE > (1usize << MORTON_COORD_BITS))
        || cells > domain.max_cells
        || cells > isize::MAX as usize / core::mem::size_of::<[f64; Q3]>()
    {
        return Err(MaterialLbmError::Invalid(
            "free-surface domain must be bounded, tile-aligned and within its explicit cell budget",
        ));
    }
    for &(tx, ty, tz) in domain.solid_tiles {
        let tile = [tx as usize, ty as usize, tz as usize];
        if (0..3).any(|axis| tile[axis] >= domain.cells[axis] / TILE) {
            return Err(MaterialLbmError::Invalid("free-surface solid tile lies outside the domain"));
        }
        for lz in 0..TILE {
            for ly in 0..TILE {
                for lx in 0..TILE {
                    if fluid(tile[0] * TILE + lx, tile[1] * TILE + ly, tile[2] * TILE + lz) {
                        return Err(MaterialLbmError::Invalid(
                            "free-surface initial liquid intersects a permanently solid tile",
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}
