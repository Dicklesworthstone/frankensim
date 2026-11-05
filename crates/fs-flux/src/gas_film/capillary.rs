//! Finite-conductance, stationary no-slip vents for an isothermal gas film.
//!
//! The long-capillary limit is equation (5) of Berg, Metrologia 42 (2005),
//! doi:10.1088/0026-1394/42/1/002. Only that ideal-gas, creeping-flow limit is
//! used here, not the paper's entrance, slip, thermal or inertial corrections.

use super::{GasFilmBoundaryTopology, GasFilmError, GasFilmInput, derived, positive};

/// Parallel straight circular capillaries connecting one film cell to a bath.
///
/// Walls are stationary and no-slip, independently of the film wall policy.
/// Flow is quasi-steady, isothermal and fully developed; neither gas storage
/// nor inertia inside the capillaries is evolved. The conservative model
/// envelope is `radius/length <= 0.01`, diameter-based `Re <= 1`, radius-based
/// `Kn <= 0.01`, and centreline `Mach <= 0.1`. Stricter caller Kn/Mach limits
/// still apply. These gates are not an error bound on omitted end effects or
/// a check that the caller's physical forcing varies slowly enough.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CapillaryVent {
    /// Internal radius of each identical capillary [m].
    pub radius_m: f64,
    /// Fully developed flow length of each capillary [m].
    pub length_m: f64,
    /// Number of parallel capillaries per out-of-plane width [m^-1].
    /// This converts each tube's mass flow to the film's per-width account.
    pub capillaries_per_m: f64,
}

/// Accepted capillary mechanics, separate from dissipation inside the film.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CapillaryVentReport {
    /// Pressure-gradient force work, equal to viscous dissipation [W m^-1]:
    /// `m_dot*R*T*ln(p_film/p_bath)`. Nonnegative for either flow direction.
    /// It is not an extra enthalpy flux or a solved thermal-energy source.
    pub viscous_dissipation_w_per_m: f64,
    /// Diameter-based Reynolds number using one tube's mean velocity.
    pub reynolds_number: f64,
    /// Maximum `lambda(p)/radius`, attained at the lower-pressure endpoint.
    pub maximum_knudsen_number: f64,
    /// Maximum centreline speed divided by the declared gas sound speed.
    pub maximum_gas_mach_number: f64,
}

impl CapillaryVent {
    pub(super) fn validate(self) -> Result<(), GasFilmError> {
        positive(self.radius_m, "capillary.radius_m")?;
        positive(self.length_m, "capillary.length_m")?;
        positive(self.capillaries_per_m, "capillary.capillaries_per_m")?;
        if self.radius_m / self.length_m > 0.01 {
            return Err(GasFilmError::Unavailable {
                reason: "capillary-outside-long-tube-envelope",
            });
        }
        Ok(())
    }
}

// Compiled once per step. The exact positive pressure derivative keeps the vent
// implicit even when its discharge time is shorter than the requested timestep.
pub(super) struct CapillaryConnection {
    pub cell: usize,
    pressure_pa: f64,
    geometry: CapillaryVent,
    mass_coefficient: f64,
    area_m2: f64,
    rt: f64,
}

impl CapillaryConnection {
    pub fn from_input(input: &GasFilmInput, rt: f64) -> Result<Option<Self>, GasFilmError> {
        let GasFilmBoundaryTopology::CapillaryVented {
            cell_index,
            absolute_pressure_pa,
            capillary,
        } = input.boundary
        else {
            return Ok(None);
        };
        let area_m2 = derived(
            core::f64::consts::PI * capillary.radius_m * capillary.radius_m,
            "capillary.area_m2",
        )?;
        positive(area_m2, "capillary.area_m2")?;
        let mass_coefficient = derived(
            capillary.capillaries_per_m * area_m2 * capillary.radius_m.powi(2)
                / (16.0 * input.gas.dynamic_viscosity_pa_s * capillary.length_m * rt),
            "capillary.mass_conductance_kg_per_m_s_pa2",
        )?;
        positive(
            mass_coefficient,
            "capillary.mass_conductance_kg_per_m_s_pa2",
        )?;
        Ok(Some(Self {
            cell: cell_index,
            pressure_pa: absolute_pressure_pa,
            geometry: capillary,
            mass_coefficient,
            area_m2,
            rt,
        }))
    }

    /// Signed outward mass flow [kg m^-1 s^-1] and its pressure derivative.
    pub fn exchange(&self, film_pressure_pa: f64) -> Result<(f64, f64), GasFilmError> {
        let flow = derived(
            self.mass_coefficient
                * (film_pressure_pa - self.pressure_pa)
                * (film_pressure_pa + self.pressure_pa),
            "capillary.outward_mass_flux_kg_per_m_s",
        )?;
        let derivative = derived(
            2.0 * self.mass_coefficient * film_pressure_pa,
            "capillary.mass_flux_pressure_derivative",
        )?;
        Ok((flow, derivative))
    }

    pub fn observe(
        &self,
        film_pressure_pa: f64,
        input: &GasFilmInput,
    ) -> Result<CapillaryVentReport, GasFilmError> {
        let lower_pressure = film_pressure_pa.min(self.pressure_pa);
        let knudsen = derived(
            (input.applicability.mean_free_path_m / self.geometry.radius_m)
                * (input.initial_absolute_pressure_pa / lower_pressure),
            "capillary.knudsen_number",
        )?;
        if knudsen > input.applicability.maximum_knudsen_number.min(0.01) {
            return Err(GasFilmError::Unavailable {
                reason: "capillary-knudsen-outside-no-slip-envelope",
            });
        }
        let (flow, _) = self.exchange(film_pressure_pa)?;
        let tube_flow = flow.abs() / self.geometry.capillaries_per_m;
        let reynolds = derived(
            2.0 * tube_flow
                / (core::f64::consts::PI
                    * self.geometry.radius_m
                    * input.gas.dynamic_viscosity_pa_s),
            "capillary.reynolds_number",
        )?;
        if reynolds > 1.0 {
            return Err(GasFilmError::Unavailable {
                reason: "capillary-reynolds-outside-creeping-envelope",
            });
        }
        // Mass flow is constant; density and therefore the parabolic peak
        // speed attain their extrema at opposite endpoints. Check that peak,
        // not the mean speed or the speed at the high-pressure inlet.
        let mach = derived(
            2.0 * tube_flow * self.rt
                / (lower_pressure * self.area_m2 * input.applicability.speed_of_sound_m_per_s),
            "capillary.maximum_gas_mach_number",
        )?;
        if mach > input.applicability.maximum_mach_number.min(0.1) {
            return Err(GasFilmError::Unavailable {
                reason: "capillary-mach-outside-quasisteady-envelope",
            });
        }
        // log1p keeps near-equilibrium pressure work, including exact zero,
        // without subtracting nearly equal logarithms.
        let dissipation = derived(
            flow * self.rt * ((film_pressure_pa - self.pressure_pa) / self.pressure_pa).ln_1p(),
            "capillary.viscous_dissipation_w_per_m",
        )?;
        Ok(CapillaryVentReport {
            viscous_dissipation_w_per_m: dissipation,
            reynolds_number: reynolds,
            maximum_knudsen_number: knudsen,
            maximum_gas_mach_number: mach,
        })
    }
}
