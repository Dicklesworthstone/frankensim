//! Shrinking spherical liquid geometry on the shared evaporation/heat owner.
//!
//! An isolated sphere in a quasi-steady stagnant binary gas has diameter-based
//! Sherwood number 2: `k_m = D_AB / r`. The original Stefan logarithm and any
//! selected interface kinetics therefore remain in `surface_flux`; this module
//! supplies geometry, not a second phase or caloric law. The classical radius
//! scaling is described by Langmuir, *The Evaporation of Small Spheres*, Phys.
//! Rev. 12, 368 (1918), <https://doi.org/10.1103/PhysRev.12.368>.
//! The stagnant `Sh=2` limit also follows from equations (2.12)--(2.16) in
//! Majee et al., JFM 956 A17 (2023), <https://doi.org/10.1017/jfm.2023.30>.
//! Concentration and Stefan-flow conventions remain those of the parent law.
//!
//! Mass and original specific enthalpy are the only retained liquid state.
//! Radius follows mass and the original material's current liquid density.
//! Uniform liquid temperature, spherical shape, continuum radial diffusion,
//! and negligible forced/buoyant convection are explicit assumptions. This is
//! not a Kelvin, rarefied, breakup, coalescence, or complete-dryout model.

use super::{
    EvaporationError, EvaporationLaw, GasFilm, LiquidSurfaceState, SolverConfig, SurfaceStep,
    advance_surface_with_heat, nonnegative, positive, properties,
};

/// A single positive pure-liquid inventory, with no independent cached radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalDropletState {
    pub mass_kg: f64,
    pub specific_enthalpy_j_kg: f64,
}

impl SphericalDropletState {
    /// Derive geometry from this inventory and the original liquid density.
    pub fn geometry<L: EvaporationLaw + ?Sized>(
        self,
        law: &L,
        temperature_k: f64,
    ) -> Result<SphericalDropletGeometry, EvaporationError> {
        if !self.specific_enthalpy_j_kg.is_finite() {
            return Err(EvaporationError::Invalid(
                "finite droplet enthalpy required",
            ));
        }
        SphericalDropletGeometry::from_mass(
            self.mass_kg,
            properties(law, temperature_k)?.liquid_density_kg_m3,
        )
    }

    fn surface(self, geometry: SphericalDropletGeometry) -> LiquidSurfaceState {
        LiquidSurfaceState {
            mass_kg: self.mass_kg,
            specific_enthalpy_j_kg: self.specific_enthalpy_j_kg,
            area_m2: geometry.area_m2,
        }
    }
}

/// Checked spherical geometry; every coordinate derives from one mass/density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalDropletGeometry {
    radius_m: f64,
    area_m2: f64,
    volume_m3: f64,
    liquid_density_kg_m3: f64,
}

impl SphericalDropletGeometry {
    pub fn from_mass(mass_kg: f64, liquid_density_kg_m3: f64) -> Result<Self, EvaporationError> {
        if !positive(mass_kg) || !positive(liquid_density_kg_m3) {
            return Err(EvaporationError::Invalid(
                "positive finite spherical liquid mass and density required",
            ));
        }
        let volume_m3 = mass_kg / liquid_density_kg_m3;
        if !positive(volume_m3) {
            return Err(EvaporationError::Unrepresentable("spherical liquid volume"));
        }
        // Factor after the root so a representable subnormal volume is not
        // lost in multiplication by 3/(4*pi). This is the geometric cube root,
        // not a fitted transfer coefficient.
        const RADIUS_FACTOR: f64 = 0.620_350_490_899_4; // cbrt(3/(4*pi))
        let radius_m = fs_math::det::pow(volume_m3, 1.0 / 3.0) * RADIUS_FACTOR;
        let area_m2 = (4.0 * core::f64::consts::PI * radius_m) * radius_m;
        if !positive(radius_m) || !positive(area_m2) {
            return Err(EvaporationError::Unrepresentable(
                "spherical liquid radius or area",
            ));
        }
        Ok(Self {
            radius_m,
            area_m2,
            volume_m3,
            liquid_density_kg_m3,
        })
    }

    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }
    pub const fn diameter_m(self) -> f64 {
        2.0 * self.radius_m
    }
    pub const fn area_m2(self) -> f64 {
        self.area_m2
    }
    pub const fn volume_m3(self) -> f64 {
        self.volume_m3
    }
    pub const fn liquid_density_kg_m3(self) -> f64 {
        self.liquid_density_kg_m3
    }

    /// Exact stagnant-sphere geometry factor for an authored binary diffusivity.
    /// Zero diffusivity seals mass transfer; heat can still change density/size.
    pub fn mass_transfer_coefficient_m_s(
        self,
        binary_diffusivity_m2_s: f64,
    ) -> Result<f64, EvaporationError> {
        if !nonnegative(binary_diffusivity_m2_s) {
            return Err(EvaporationError::Invalid(
                "finite nonnegative binary molecular diffusivity required",
            ));
        }
        let coefficient = binary_diffusivity_m2_s / self.radius_m;
        if !nonnegative(coefficient) || (binary_diffusivity_m2_s > 0.0 && coefficient == 0.0) {
            return Err(EvaporationError::Unrepresentable(
                "spherical molecular conductance",
            ));
        }
        Ok(coefficient)
    }

    /// Convert an explicitly authored diameter-based Sherwood number to `k_m`.
    ///
    /// `D_AB` remains molecular diffusivity; it is never replaced with an
    /// effective diffusivity. This geometric conversion supplies no correlation
    /// or validity claim. The caller must admit the selected transport model.
    /// `Sh=2` uses the original stagnant-sphere arithmetic exactly. A zero
    /// molecular diffusivity seals mass transfer; Sherwood must be positive.
    pub fn mass_transfer_coefficient_at_sherwood_m_s(
        self,
        binary_diffusivity_m2_s: f64,
        sherwood_number: f64,
    ) -> Result<f64, EvaporationError> {
        if !positive(sherwood_number) || !nonnegative(binary_diffusivity_m2_s) {
            return Err(EvaporationError::Invalid(
                "positive finite Sherwood number and nonnegative molecular diffusivity required",
            ));
        }
        if sherwood_number == 2.0 {
            return self.mass_transfer_coefficient_m_s(binary_diffusivity_m2_s);
        }
        if binary_diffusivity_m2_s == 0.0 {
            return Ok(0.0);
        }
        let diameter = self.diameter_m();
        let coefficient = (binary_diffusivity_m2_s / diameter) * sherwood_number;
        // A representable result can survive a different multiplication order
        // even when the first positive intermediate overflowed/underflowed.
        let coefficient = if positive(coefficient) {
            coefficient
        } else {
            (binary_diffusivity_m2_s * sherwood_number) / diameter
        };
        if !positive(coefficient) {
            return Err(EvaporationError::Unrepresentable(
                "spherical Sherwood conductance",
            ));
        }
        Ok(coefficient)
    }
}

/// Prescribed far-field gas coordinates and molecular binary diffusivity.
/// The finite-gas owner supplies its own actual coordinates instead.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalGasFilm {
    pub temperature_k: f64,
    pub total_pressure_pa: f64,
    pub vapor_partial_pressure_pa: f64,
    pub binary_diffusivity_m2_s: f64,
}

impl SphericalGasFilm {
    pub fn at_geometry(
        self,
        geometry: SphericalDropletGeometry,
    ) -> Result<GasFilm, EvaporationError> {
        let film = GasFilm {
            temperature_k: self.temperature_k,
            total_pressure_pa: self.total_pressure_pa,
            vapor_partial_pressure_pa: self.vapor_partial_pressure_pa,
            mass_transfer_coefficient_m_s: geometry
                .mass_transfer_coefficient_m_s(self.binary_diffusivity_m2_s)?,
        };
        film.validate()?;
        Ok(film)
    }
}

/// Accepted state and both actual spherical geometries around one heat step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalDropletStep {
    pub droplet: SphericalDropletState,
    pub geometry: [SphericalDropletGeometry; 2],
    /// Original heat/phase receipt using the interval's initial exposed area.
    /// Its planar `recession_m` is not the spherical radius change.
    pub phase: SurfaceStep,
}

/// Advance a sphere against a prescribed bath using the original heat balance.
///
/// Radius, area and D/r are frozen over this interval; the original phase/heat
/// solve is implicit in liquid enthalpy. Accepted geometry is reconstructed
/// from final mass AND final material density, then used by the next interval.
/// This geometry split is first order: refine the whole timestep. It is not an
/// exact d-squared integrator, and finite humidity/latent cooling/kinetics need
/// not follow a constant d-squared law. Supplied heat is counted only once.
/// An exhausted positive liquid seed refuses through the original phase owner.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplet_with_heat<L: EvaporationLaw + ?Sized>(
    law: &L,
    film: SphericalGasFilm,
    droplet: SphericalDropletState,
    duration_s: f64,
    supplied_heat_j: f64,
    config: SolverConfig,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    mut checkpoint: impl FnMut() -> Result<(), EvaporationError>,
) -> Result<SphericalDropletStep, EvaporationError> {
    checkpoint()?;
    let initial_t = temperature(droplet.specific_enthalpy_j_kg)?;
    let initial = droplet.geometry(law, initial_t)?;
    let phase = advance_surface_with_heat(
        law,
        film.at_geometry(initial)?,
        droplet.surface(initial),
        duration_s,
        supplied_heat_j,
        config,
        temperature,
        &mut checkpoint,
    )?;
    checkpoint()?;
    let droplet = SphericalDropletState {
        mass_kg: phase.liquid.mass_kg,
        specific_enthalpy_j_kg: phase.liquid.specific_enthalpy_j_kg,
    };
    let final_geometry = droplet.geometry(law, phase.surface_temperature_k)?;
    checkpoint()?;
    Ok(SphericalDropletStep {
        droplet,
        geometry: [initial, final_geometry],
        phase,
    })
}

#[cfg(test)]
mod tests;
