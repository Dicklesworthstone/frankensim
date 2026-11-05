//! Finite-conductance, stationary no-slip vents for an isothermal gas film.
//!
//! The long-capillary limit is equation (5) of Berg, Metrologia 42 (2005),
//! doi:10.1088/0026-1394/42/1/002. Only that ideal-gas, creeping-flow limit is
//! used here, not the paper's entrance, slip, thermal or inertial corrections.

use super::{
    GasFilmBoundaryTopology, GasFilmCheckpoint, GasFilmError, GasFilmInput, derived, positive,
};

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

/// Finite rigid reservoir containing the same isothermal ideal gas as the film.
/// Its mass is advanced through the capillary; its pressure is never prescribed
/// after initialization. Temperature remains the film's declared fixed value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmReservoir {
    /// Reservoir volume per out-of-plane width [m^3/m = m^2].
    pub volume_per_width_m2: f64,
    /// Initial absolute pressure [Pa], used only without a checkpoint.
    pub initial_absolute_pressure_pa: f64,
}

/// Storage and transfer of the finite reservoir at the accepted endpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmReservoirReport {
    /// Absolute pressure from the accepted reservoir mass and ideal-gas EOS [Pa].
    pub absolute_pressure_pa: f64,
    /// Accepted reservoir gas mass per out-of-plane width [kg m^-1].
    pub gas_mass_kg_per_m: f64,
    /// Represented mass-storage rate [kg m^-1 s^-1].
    pub storage_rate_kg_per_m_s: f64,
    /// Signed enthalpy entering the reservoir from the film [W m^-1].
    /// This is the film's outward vent enthalpy with the receiving convention.
    pub inward_enthalpy_flux_w_per_m: f64,
    /// Heat into the reservoir required to maintain its fixed temperature [W m^-1].
    /// Equals `-R*T*m_dot`; filling rejects heat and draining requires heat.
    /// This is a thermal requirement, not a solved wall heat-transfer law.
    pub required_isothermal_heat_to_gas_w_per_m: f64,
    /// Film plus reservoir storage rate [kg m^-1 s^-1]. Their capillary
    /// exchange is internal; both exterior film ends are sealed.
    pub combined_mass_closure_residual_kg_per_m_s: f64,
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

#[derive(Clone, Copy)]
struct ReservoirStorage {
    old_mass_kg_per_m: f64,
    capacity_kg_per_m_pa: f64,
    storage_coefficient: f64,
    duration_s: f64,
}

#[derive(Clone, Copy)]
pub(super) struct CapillaryExchange {
    pub flow_kg_per_m_s: f64,
    pub derivative: f64,
    exterior_pressure_pa: f64,
    pressure_drop_pa: f64,
    pub reservoir_mass_kg_per_m: Option<f64>,
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
    reservoir: Option<ReservoirStorage>,
}

impl CapillaryConnection {
    pub fn from_input(
        input: &GasFilmInput,
        rt: f64,
        checkpoint: Option<&GasFilmCheckpoint>,
    ) -> Result<Option<Self>, GasFilmError> {
        let (cell_index, absolute_pressure_pa, capillary, reservoir) = match input.boundary {
            GasFilmBoundaryTopology::CapillaryVented {
                cell_index,
                absolute_pressure_pa,
                capillary,
            } => (cell_index, absolute_pressure_pa, capillary, None),
            GasFilmBoundaryTopology::CapillaryReservoir {
                cell_index,
                capillary,
                reservoir,
            } => {
                let capacity = derived(
                    reservoir.volume_per_width_m2 / rt,
                    "reservoir.mass_pressure_capacity",
                )?;
                positive(capacity, "reservoir.mass_pressure_capacity")?;
                let mass = match checkpoint {
                    Some(state) => state.capillary_reservoir_mass_kg_per_m.ok_or(
                        GasFilmError::CheckpointMismatch {
                            field: "capillary_reservoir_mass_kg_per_m",
                        },
                    )?,
                    None => derived(
                        capacity * reservoir.initial_absolute_pressure_pa,
                        "reservoir.initial_mass_kg_per_m",
                    )?,
                };
                positive(mass, "reservoir.gas_mass_kg_per_m")?;
                let pressure = derived(mass / capacity, "reservoir.absolute_pressure_pa")?;
                positive(pressure, "reservoir.absolute_pressure_pa")?;
                let storage_coefficient =
                    derived(capacity / input.timestep_s, "reservoir.storage_coefficient")?;
                positive(storage_coefficient, "reservoir.storage_coefficient")?;
                (
                    cell_index,
                    pressure,
                    capillary,
                    Some(ReservoirStorage {
                        old_mass_kg_per_m: mass,
                        capacity_kg_per_m_pa: capacity,
                        storage_coefficient,
                        duration_s: input.timestep_s,
                    }),
                )
            }
            _ => return Ok(None),
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
            reservoir,
        }))
    }

    /// Signed outward mass flow [kg m^-1 s^-1] and its pressure derivative.
    fn prescribed_exchange(&self, film_pressure_pa: f64) -> Result<(f64, f64), GasFilmError> {
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

    pub fn exchange(&self, film_pressure_pa: f64) -> Result<CapillaryExchange, GasFilmError> {
        let Some(reservoir) = self.reservoir else {
            let (flow, derivative) = self.prescribed_exchange(film_pressure_pa)?;
            return Ok(CapillaryExchange {
                flow_kg_per_m_s: flow,
                derivative,
                exterior_pressure_pa: self.pressure_pa,
                pressure_drop_pa: film_pressure_pa - self.pressure_pa,
                reservoir_mass_kg_per_m: None,
            });
        };
        // Eliminate the reservoir's backward-Euler storage equation exactly:
        // G*p_r^2 + B*p_r = G*p_f^2 + B*p_r_old, B = V/(R*T*dt).
        // Scale pressures and coefficients before the positive quadratic root.
        // No outer iteration, explicit discharge step or pressure clamp is used.
        let pressure_scale = self.pressure_pa.max(film_pressure_pa);
        let conductance_scale = derived(
            self.mass_coefficient * pressure_scale,
            "reservoir.conductance_scale",
        )?;
        let coefficient_scale = reservoir.storage_coefficient.max(conductance_scale);
        let a = conductance_scale / coefficient_scale;
        let b = reservoir.storage_coefficient / coefficient_scale;
        let film = film_pressure_pa / pressure_scale;
        let old = self.pressure_pa / pressure_scale;
        let c = a * film * film + b * old;
        let denominator = 0.5 * b + (0.5 * b).hypot(a.sqrt() * c.sqrt());
        let root_pressure = derived(
            pressure_scale * (c / denominator),
            "reservoir.implicit_pressure_pa",
        )?;
        positive(root_pressure, "reservoir.implicit_pressure_pa")?;
        // The series form avoids subtracting p_f-p_r when a tiny reservoir
        // equilibrates nearly completely in one step. It also retains a small
        // flux when a large reservoir's pressure change rounds away.
        let tube = derived(
            self.mass_coefficient * (film_pressure_pa + root_pressure),
            "reservoir.tube_conductance",
        )?;
        let flow = derived(
            series_conductance(tube, reservoir.storage_coefficient)
                * (film_pressure_pa - self.pressure_pa),
            "capillary.outward_mass_flux_kg_per_m_s",
        )?;
        let mass = derived(
            reservoir.old_mass_kg_per_m + reservoir.duration_s * flow,
            "reservoir.accepted_mass_kg_per_m",
        )?;
        positive(mass, "reservoir.accepted_mass_kg_per_m")?;
        let pressure = derived(
            mass / reservoir.capacity_kg_per_m_pa,
            "reservoir.accepted_pressure_pa",
        )?;
        positive(pressure, "reservoir.accepted_pressure_pa")?;
        let slope = derived(
            2.0 * self.mass_coefficient * pressure,
            "reservoir.pressure_derivative",
        )?;
        let derivative = derived(
            series_conductance(slope, reservoir.storage_coefficient)
                * (film_pressure_pa / pressure),
            "capillary.mass_flux_pressure_derivative",
        )?;
        Ok(CapillaryExchange {
            flow_kg_per_m_s: flow,
            derivative,
            exterior_pressure_pa: pressure,
            pressure_drop_pa: derived(flow / tube, "reservoir.capillary_pressure_drop_pa")?,
            reservoir_mass_kg_per_m: Some(mass),
        })
    }

    pub fn admit_initial(
        &self,
        film_pressure_pa: f64,
        input: &GasFilmInput,
    ) -> Result<(), GasFilmError> {
        let (flow, _) = self.prescribed_exchange(film_pressure_pa)?;
        self.observe_flow(
            film_pressure_pa,
            self.pressure_pa,
            film_pressure_pa - self.pressure_pa,
            flow,
            input,
        )?;
        Ok(())
    }

    pub fn observe(
        &self,
        film_pressure_pa: f64,
        input: &GasFilmInput,
    ) -> Result<CapillaryVentReport, GasFilmError> {
        let exchange = self.exchange(film_pressure_pa)?;
        self.observe_flow(
            film_pressure_pa,
            exchange.exterior_pressure_pa,
            exchange.pressure_drop_pa,
            exchange.flow_kg_per_m_s,
            input,
        )
    }

    fn observe_flow(
        &self,
        film_pressure_pa: f64,
        exterior_pressure_pa: f64,
        pressure_drop_pa: f64,
        flow: f64,
        input: &GasFilmInput,
    ) -> Result<CapillaryVentReport, GasFilmError> {
        let lower_pressure = film_pressure_pa.min(exterior_pressure_pa);
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
            flow * self.rt * (pressure_drop_pa / exterior_pressure_pa).ln_1p(),
            "capillary.viscous_dissipation_w_per_m",
        )?;
        Ok(CapillaryVentReport {
            viscous_dissipation_w_per_m: dissipation,
            reynolds_number: reynolds,
            maximum_knudsen_number: knudsen,
            maximum_gas_mach_number: mach,
        })
    }

    pub fn reservoir_residual(&self, exchange: CapillaryExchange) -> Result<f64, GasFilmError> {
        match (self.reservoir, exchange.reservoir_mass_kg_per_m) {
            (Some(storage), Some(mass)) => derived(
                (mass - storage.old_mass_kg_per_m) / storage.duration_s - exchange.flow_kg_per_m_s,
                "reservoir.mass_closure_residual_kg_per_m_s",
            ),
            _ => Ok(0.0),
        }
    }

    pub fn reservoir_report(
        &self,
        exchange: CapillaryExchange,
        input: &GasFilmInput,
        film_storage_rate: f64,
    ) -> Result<Option<GasFilmReservoirReport>, GasFilmError> {
        let (Some(storage), Some(mass)) = (self.reservoir, exchange.reservoir_mass_kg_per_m) else {
            return Ok(None);
        };
        let rate = derived(
            (mass - storage.old_mass_kg_per_m) / storage.duration_s,
            "reservoir.storage_rate_kg_per_m_s",
        )?;
        Ok(Some(GasFilmReservoirReport {
            absolute_pressure_pa: exchange.exterior_pressure_pa,
            gas_mass_kg_per_m: mass,
            storage_rate_kg_per_m_s: rate,
            inward_enthalpy_flux_w_per_m: derived(
                exchange.flow_kg_per_m_s * input.gas.declared_specific_enthalpy_j_kg,
                "reservoir.inward_enthalpy_flux_w_per_m",
            )?,
            required_isothermal_heat_to_gas_w_per_m: derived(
                -exchange.flow_kg_per_m_s * self.rt,
                "reservoir.required_isothermal_heat_to_gas_w_per_m",
            )?,
            combined_mass_closure_residual_kg_per_m_s: derived(
                film_storage_rate + rate,
                "reservoir.combined_mass_closure_residual_kg_per_m_s",
            )?,
        }))
    }
}

// ab/(a+b), without multiplying or adding potentially disparate conductances.
fn series_conductance(a: f64, b: f64) -> f64 {
    let low = a.min(b);
    let high = a.max(b);
    low / (1.0 + low / high)
}
