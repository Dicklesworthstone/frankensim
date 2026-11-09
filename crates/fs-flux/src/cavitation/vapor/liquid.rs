//! Resolved liquid heat transport outside a spherical phase-changing interface.
//!
//! Shell faces retain `r_i^3-R^3`, so every cell has fixed liquid volume and
//! mass while the radial mesh follows the moving bubble. The relative mass
//! flux through EVERY face is `F_m=-dm_v/dt`, outward-positive. The same
//! original saturation-path liquid enthalpy supplies storage, interior upwind
//! transport, and both exact Dirichlet boundary traces. Interface conductive
//! heat and transferred enthalpy therefore cancel their native vapor accounts.
//!
//! The existing material geometry supplies the spherical resistance integral.
//! Cell conductivities come from the actual original IF97/transport state;
//! neighboring half-shell resistances combine in series. Boundary conductivity
//! is the mean of its original face/cell endpoint values. This two-point
//! finite-volume approximation requires independent spatial refinement.
//! An outflow Dirichlet trace is admitted only when its half-cell Peclet scale
//! fits the declared bound; a coarse outflow boundary must be refined.
//!
//! Storage follows h_l,sat(T), an explicitly near-saturation incompressible
//! caloric approximation, not an assignment of mechanical pressure to p_sat(T).
//! Actual source-property drift, saturation-path versus isobaric capacity and
//! estimated omitted reversible pressure heating remain checked. The native
//! spherical velocity field supplies its resolved viscous heating exactly.
//! The unrepresented exterior liquid remains an explicit thermal reservoir.
//! This spatial operator owns no time integrator or separately accepted clock.

use super::{
    VaporBubbleError, VaporBubbleObservation, VaporBubbleState, VaporBubbleThermalState, bound,
};
use core::f64::consts::PI;
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::cavitation::vapor::liquid::SaturatedWaterLiquidThermal;
use fs_material::thermal::conduction::ReservoirLayerGeometry;

pub const MAX_LIQUID_THERMAL_CELLS: usize = 64;
const MAX_FACES: usize = MAX_LIQUID_THERMAL_CELLS + 1;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(what: &'static str) -> VaporBubbleError {
    VaporBubbleError::InvalidInput(what)
}
fn finite(value: f64, what: &'static str) -> Result<f64, VaporBubbleError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(what))
    }
}

/// Fixed liquid-volume coordinates with explicit reference physical faces.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalLiquidGrid {
    cells: usize,
    reference_radius_m: f64,
    cubic_offsets_m3: [f64; MAX_FACES],
    volumes_m3: [f64; MAX_LIQUID_THERMAL_CELLS],
}
impl SphericalLiquidGrid {
    /// N+1 strictly increasing radii, starting at the bubble's reference radius.
    /// Arbitrary positive grading permits resolution of the thermal layer.
    pub fn from_face_radii(radii: &[f64]) -> Result<Self, VaporBubbleError> {
        if !(2..=MAX_FACES).contains(&radii.len())
            || radii.iter().any(|&r| !positive(r))
            || radii.windows(2).any(|p| p[0] >= p[1])
        {
            return Err(invalid(
                "ordered positive spherical liquid faces, 1..=64 cells",
            ));
        }
        let mut result = Self {
            cells: radii.len() - 1,
            reference_radius_m: radii[0],
            cubic_offsets_m3: [0.0; MAX_FACES],
            volumes_m3: [0.0; MAX_LIQUID_THERMAL_CELLS],
        };
        let base = radii[0];
        for i in 1..radii.len() {
            let r = radii[i];
            let previous = radii[i - 1];
            let offset = (r - base) * (r * r + r * base + base * base);
            let volume =
                (4.0 * PI / 3.0) * (r - previous) * (r * r + r * previous + previous * previous);
            if !positive(offset) || !positive(volume) || offset <= result.cubic_offsets_m3[i - 1] {
                return Err(invalid("representable spherical liquid volume coordinates"));
            }
            result.cubic_offsets_m3[i] = offset;
            result.volumes_m3[i - 1] = volume;
        }
        Ok(result)
    }
    pub fn equal_volume(inner: f64, outer: f64, cells: usize) -> Result<Self, VaporBubbleError> {
        if !positive(inner)
            || !positive(outer)
            || outer <= inner
            || !(1..=MAX_LIQUID_THERMAL_CELLS).contains(&cells)
        {
            return Err(invalid("positive spherical liquid extent and 1..=64 cells"));
        }
        let span = (outer - inner) * (outer * outer + outer * inner + inner * inner);
        let base = inner * inner * inner;
        if !positive(span) || !positive(base) {
            return Err(invalid("represented liquid extent"));
        }
        let mut faces = [0.0; MAX_FACES];
        faces[0] = inner;
        for (i, face) in faces.iter_mut().enumerate().take(cells).skip(1) {
            *face = (base + span * (i as f64 / cells as f64)).cbrt();
        }
        faces[cells] = outer;
        Self::from_face_radii(&faces[..=cells])
    }
    pub const fn cell_count(&self) -> usize {
        self.cells
    }
    pub const fn reference_radius_m(&self) -> f64 {
        self.reference_radius_m
    }
    pub fn volumes_m3(&self) -> &[f64] {
        &self.volumes_m3[..self.cells]
    }

    fn geometry(&self, radius: f64) -> Result<Geometry, VaporBubbleError> {
        if !positive(radius) {
            return Err(invalid("positive current liquid interface radius"));
        }
        let base = radius * radius * radius;
        if !positive(base) {
            return Err(invalid("represented liquid interface volume"));
        }
        let mut result = Geometry {
            faces: [0.0; MAX_FACES],
            centres: [0.0; MAX_LIQUID_THERMAL_CELLS],
        };
        result.faces[0] = radius;
        for i in 1..=self.cells {
            result.faces[i] = (base + self.cubic_offsets_m3[i]).cbrt();
            let lo = result.faces[i - 1];
            let hi = result.faces[i];
            if !positive(hi) || hi <= lo {
                return Err(invalid("resolved positive current liquid cell widths"));
            }
            // Exact volume centroid, factoring both power differences.
            let centre = 0.75 * (lo + hi) * ((lo * lo + hi * hi) / (lo * lo + lo * hi + hi * hi));
            if !centre.is_finite() || centre <= lo || centre >= hi {
                return Err(invalid("represented spherical liquid volume centroid"));
            }
            result.centres[i - 1] = centre;
        }
        Ok(result)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidThermalLimits {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Original rho, cp, k, mu and c drift from the retained reference, <=0.05.
    pub maximum_source_relative_change: f64,
    /// |dh_sat/dT-cp|/cp, <=0.05; the source pressure path is not mechanical p.
    pub maximum_pressure_path_capacity_fraction: f64,
    /// Estimated |beta T Delta p|/(rho cp), in K, for omitted pressure heating.
    pub maximum_pressure_heating_temperature_k: f64,
    /// Outflow half-cell |mdot| (dh/dT)/G, in (0,1].
    pub maximum_outflow_peclet: f64,
}
impl LiquidThermalLimits {
    pub fn validate(self) -> Result<(), VaporBubbleError> {
        if !self.minimum_temperature_k.is_finite()
            || !self.maximum_temperature_k.is_finite()
            || self.minimum_temperature_k < 273.16
            || self.maximum_temperature_k > 623.15
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || !positive(self.maximum_pressure_heating_temperature_k)
        {
            return Err(invalid(
                "original liquid thermal domain and pressure-heating budget",
            ));
        }
        for (value, maximum) in [
            (self.maximum_source_relative_change, 0.05),
            (self.maximum_pressure_path_capacity_fraction, 0.05),
            (self.maximum_outflow_peclet, 1.0),
        ] {
            if !positive(value) || value > maximum {
                return Err(invalid("bounded liquid thermal approximation limits"));
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SphericalLiquidField {
    source: SaturatedWaterCavitation,
    grid: SphericalLiquidGrid,
    reference: SaturatedWaterLiquidThermal,
    limits: LiquidThermalLimits,
}
impl SphericalLiquidField {
    pub fn new(
        source: SaturatedWaterCavitation,
        grid: SphericalLiquidGrid,
        limits: LiquidThermalLimits,
    ) -> Result<Self, VaporBubbleError> {
        limits.validate()?;
        let reference = SaturatedWaterLiquidThermal::try_from_saturation(*source.saturation())?;
        let result = Self {
            source,
            grid,
            reference,
            limits,
        };
        result.admit(reference)?;
        result.grid.geometry(result.grid.reference_radius_m)?;
        Ok(result)
    }
    pub const fn source(&self) -> &SaturatedWaterCavitation {
        &self.source
    }
    pub const fn grid(&self) -> &SphericalLiquidGrid {
        &self.grid
    }
    pub const fn limits(&self) -> LiquidThermalLimits {
        self.limits
    }
    pub fn reference_enthalpy_j_kg(&self) -> f64 {
        self.reference.enthalpy_j_kg()
    }

    fn admit(&self, point: SaturatedWaterLiquidThermal) -> Result<f64, VaporBubbleError> {
        let t = point.temperature_k();
        if t < self.limits.minimum_temperature_k || t > self.limits.maximum_temperature_k {
            return Err(VaporBubbleError::TemperatureOutsideDomain {
                temperature_k: t,
                minimum_k: self.limits.minimum_temperature_k,
                maximum_k: self.limits.maximum_temperature_k,
            });
        }
        bound(
            "liquid saturation-path heat-capacity correction",
            point.pressure_path_capacity_fraction(),
            self.limits.maximum_pressure_path_capacity_fraction,
        )?;
        let mut change: f64 = 0.0;
        for (value, reference) in [
            (point.density_kg_m3(), self.reference.density_kg_m3()),
            (point.cp_j_kg_k(), self.reference.cp_j_kg_k()),
            (
                point.conductivity_w_m_k(),
                self.reference.conductivity_w_m_k(),
            ),
            (
                point.dynamic_viscosity_pa_s(),
                self.reference.dynamic_viscosity_pa_s(),
            ),
            (point.sound_speed_m_s(), self.reference.sound_speed_m_s()),
        ] {
            change = change.max((value / reference - 1.0).abs());
        }
        bound(
            "resolved liquid source variation",
            change,
            self.limits.maximum_source_relative_change,
        )?;
        Ok(change)
    }
    fn sample(&self, temperature: f64) -> Result<SaturatedWaterLiquidThermal, VaporBubbleError> {
        if !temperature.is_finite()
            || temperature < self.limits.minimum_temperature_k
            || temperature > self.limits.maximum_temperature_k
        {
            return Err(VaporBubbleError::TemperatureOutsideDomain {
                temperature_k: temperature,
                minimum_k: self.limits.minimum_temperature_k,
                maximum_k: self.limits.maximum_temperature_k,
            });
        }
        let state = if temperature == self.reference.temperature_k() {
            self.reference
        } else {
            SaturatedWaterLiquidThermal::at_temperature_k(temperature)?
        };
        self.admit(state)?;
        Ok(state)
    }

    /// Sample original sources once and prepare heat before native phase flux.
    pub fn prepare(
        &self,
        thermal: &VaporBubbleThermalState,
        temperatures: &[f64; MAX_LIQUID_THERMAL_CELLS],
        bath_temperature_k: f64,
    ) -> Result<LiquidThermalStage, VaporBubbleError> {
        if thermal.state.temperature_k() != thermal.vapor.temperature_k() {
            return Err(invalid("coherent vapor/liquid interface temperature"));
        }
        let interface = if thermal.state.temperature_k() == self.reference.temperature_k() {
            self.reference
        } else {
            SaturatedWaterLiquidThermal::try_from_saturation(*thermal.vapor.saturation())?
        };
        self.admit(interface)?;
        let bath = self.sample(bath_temperature_k)?;
        let geometry = self.grid.geometry(thermal.state.radius_m())?;
        let mut result = LiquidThermalStage {
            state: thermal.state,
            cells: self.grid.cells,
            geometry,
            reference_h: self.reference.enthalpy_j_kg(),
            temperature: [0.0; MAX_LIQUID_THERMAL_CELLS],
            enthalpy: [0.0; MAX_LIQUID_THERMAL_CELLS],
            capacity: [0.0; MAX_LIQUID_THERMAL_CELLS],
            mass: [0.0; MAX_LIQUID_THERMAL_CELLS],
            conductivity: [0.0; MAX_LIQUID_THERMAL_CELLS],
            pressure_heating_factor: [0.0; MAX_LIQUID_THERMAL_CELLS],
            compressibility: [0.0; MAX_LIQUID_THERMAL_CELLS],
            saturation_pressure: [0.0; MAX_LIQUID_THERMAL_CELLS],
            conductance: [0.0; MAX_FACES],
            conductive_power: [0.0; MAX_FACES],
            interface_h: interface.enthalpy_j_kg(),
            bath_h: bath.enthalpy_j_kg(),
            interface_capacity: interface.enthalpy_temperature_derivative_j_kg_k(),
            bath_capacity: bath.enthalpy_temperature_derivative_j_kg_k(),
            stored_energy_j: 0.0,
            maximum_source_relative_change: 0.0,
            maximum_pressure_path_capacity_fraction: 0.0,
        };
        for (i, &temperature) in temperatures.iter().enumerate().take(self.grid.cells) {
            let point = self.sample(temperature)?;
            let mass = self.source.state().liquid_density_kg_m3() * self.grid.volumes_m3[i];
            let capacity = mass * point.enthalpy_temperature_derivative_j_kg_k();
            if !positive(mass) || !positive(capacity) {
                return Err(invalid("represented liquid thermal inventory"));
            }
            result.temperature[i] = temperature;
            result.enthalpy[i] = point.enthalpy_j_kg();
            result.mass[i] = mass;
            result.capacity[i] = capacity;
            result.conductivity[i] = point.conductivity_w_m_k();
            result.pressure_heating_factor[i] = (point.thermal_expansion_per_k() * temperature)
                .abs()
                / (point.density_kg_m3() * point.cp_j_kg_k());
            result.compressibility[i] = point.isothermal_compressibility_per_pa();
            result.saturation_pressure[i] = point.saturation_pressure_pa();
            result.maximum_source_relative_change = result
                .maximum_source_relative_change
                .max(self.admit(point)?);
            result.maximum_pressure_path_capacity_fraction = result
                .maximum_pressure_path_capacity_fraction
                .max(point.pressure_path_capacity_fraction());
            result.stored_energy_j += mass * (point.enthalpy_j_kg() - result.reference_h);
        }
        let n = self.grid.cells;
        result.conductance[0] = metric(result.geometry.faces[0], result.geometry.centres[0])?
            * f64::midpoint(interface.conductivity_w_m_k(), result.conductivity[0]);
        result.conductive_power[0] =
            result.conductance[0] * (interface.temperature_k() - result.temperature[0]);
        for i in 1..n {
            let face = result.geometry.faces[i];
            let left = metric(face, result.geometry.centres[i - 1])? * result.conductivity[i - 1];
            let right = metric(face, result.geometry.centres[i])? * result.conductivity[i];
            // Harmonic combination, avoiding product overflow.
            let lower = left.min(right);
            result.conductance[i] = lower / (1.0 + lower / left.max(right));
            result.conductive_power[i] =
                result.conductance[i] * (result.temperature[i - 1] - result.temperature[i]);
        }
        result.conductance[n] = metric(result.geometry.faces[n], result.geometry.centres[n - 1])?
            * f64::midpoint(bath.conductivity_w_m_k(), result.conductivity[n - 1]);
        result.conductive_power[n] =
            result.conductance[n] * (result.temperature[n - 1] - bath.temperature_k());
        if !result.stored_energy_j.is_finite()
            || result.conductance[..=n].iter().any(|&v| !positive(v))
            || result.conductive_power[..=n].iter().any(|v| !v.is_finite())
        {
            return Err(invalid("represented resolved liquid heat powers"));
        }
        Ok(result)
    }

    /// Add the SAME native phase mass flux and radial mechanical dissipation.
    pub fn finish(
        &self,
        stage: LiquidThermalStage,
        observation: &VaporBubbleObservation,
        far_pressure_pa: f64,
    ) -> Result<LiquidThermalRates, VaporBubbleError> {
        if stage.state != observation.state
            || observation.heat_power_w != stage.interface_heat_power_w()
            || !far_pressure_pa.is_finite()
            || stage.cells != self.grid.cells
        {
            return Err(invalid(
                "same liquid/vapor thermal stage and physical pressure",
            ));
        }
        let mass_flux = -observation.vapor_mass_rate_kg_s;
        if !mass_flux.is_finite() {
            return Err(invalid("finite liquid Stefan mass flux"));
        }
        let n = stage.cells;
        let mut advective = [0.0; MAX_FACES];
        advective[0] = mass_flux * (stage.interface_h - stage.reference_h);
        advective[n] = mass_flux * (stage.bath_h - stage.reference_h);
        for (i, value) in advective.iter_mut().enumerate().take(n).skip(1) {
            let upstream = if mass_flux >= 0.0 { i - 1 } else { i };
            *value = mass_flux * (stage.enthalpy[upstream] - stage.reference_h);
        }
        let outflow_peclet = if mass_flux < 0.0 {
            -mass_flux
                * stage
                    .interface_capacity
                    .max(stage.capacity[0] / stage.mass[0])
                / stage.conductance[0]
        } else {
            mass_flux
                * stage
                    .bath_capacity
                    .max(stage.capacity[n - 1] / stage.mass[n - 1])
                / stage.conductance[n]
        };
        bound(
            "outflow Dirichlet liquid half-cell Peclet number",
            outflow_peclet,
            self.limits.maximum_outflow_peclet,
        )?;
        let rho = self.source.state().liquid_density_kg_m3();
        let viscosity = self.source.state().dynamic_viscosity_pa_s();
        let radial_flux = observation.liquid_volume_flow_m3_s / (4.0 * PI);
        let radial_flux_rate = observation.liquid_volume_flow_rate_m3_s2 / (4.0 * PI);
        let viscous_factor = 16.0 * PI * viscosity * radial_flux * radial_flux;
        let mut result = LiquidThermalRates {
            temperature_rate_k_s: [0.0; MAX_LIQUID_THERMAL_CELLS],
            interface_heat_power_w: stage.interface_heat_power_w(),
            interface_enthalpy_power_w: -advective[0],
            outer_conductive_power_w: stage.conductive_power[n],
            outer_advected_power_w: advective[n],
            stored_energy_j: stage.stored_energy_j,
            energy_rate_residual_w: 0.0,
            maximum_outflow_peclet: outflow_peclet,
            maximum_pressure_heating_temperature_k: 0.0,
            maximum_source_relative_change: stage.maximum_source_relative_change,
            maximum_pressure_path_capacity_fraction: stage.maximum_pressure_path_capacity_fraction,
            minimum_transport_time_s: f64::INFINITY,
            viscous_heat_power_w: 0.0,
            exterior_viscous_heat_power_w: 0.0,
            cells: n,
            enthalpy: stage.enthalpy,
            mass: stage.mass,
            reference_h: stage.reference_h,
        };
        let mut stored_power = 0.0;
        for i in 0..n {
            let radius = stage.geometry.centres[i];
            let pressure = far_pressure_pa
                + rho
                    * (radial_flux_rate / radius - 0.5 * (radial_flux / (radius * radius)).powi(2));
            let pressure_departure = (pressure - stage.saturation_pressure[i]).abs();
            let temperature_scale = stage.pressure_heating_factor[i]
                * (pressure - self.reference.saturation_pressure_pa())
                    .abs()
                    .max(pressure_departure);
            bound(
                "omitted liquid pressure-heating temperature scale",
                temperature_scale,
                self.limits.maximum_pressure_heating_temperature_k,
            )?;
            bound(
                "liquid pressure-density correction",
                stage.compressibility[i] * pressure_departure,
                self.limits.maximum_source_relative_change,
            )?;
            result.maximum_pressure_heating_temperature_k = result
                .maximum_pressure_heating_temperature_k
                .max(temperature_scale);
            let inner = stage.geometry.faces[i];
            let outer = stage.geometry.faces[i + 1];
            let inverse_cube_difference = ((outer - inner) / outer)
                * (1.0 + inner / outer + (inner / outer).powi(2))
                / (inner * inner * inner);
            let viscous_heat = viscous_factor * inverse_cube_difference;
            if !viscous_heat.is_finite() || viscous_heat < 0.0 {
                return Err(invalid("resolved native liquid viscous heating"));
            }
            let power = stage.conductive_power[i] - stage.conductive_power[i + 1] + advective[i]
                - advective[i + 1]
                + viscous_heat;
            result.temperature_rate_k_s[i] = finite(
                power / stage.capacity[i],
                "resolved liquid temperature rate",
            )?;
            result.viscous_heat_power_w += viscous_heat;
            stored_power += stage.capacity[i] * result.temperature_rate_k_s[i];
            let decay = (stage.conductance[i] + stage.conductance[i + 1]) / stage.capacity[i]
                + mass_flux.abs() / stage.mass[i];
            if !positive(decay) {
                return Err(invalid("liquid thermal transport time scale"));
            }
            result.minimum_transport_time_s = result.minimum_transport_time_s.min(1.0 / decay);
        }
        let outer = stage.geometry.faces[n];
        result.exterior_viscous_heat_power_w = finite(
            viscous_factor / (outer * outer * outer),
            "exterior native viscous heating",
        )?;
        result.energy_rate_residual_w = finite(
            stored_power
                + result.interface_heat_power_w
                + result.interface_enthalpy_power_w
                + result.outer_conductive_power_w
                + result.outer_advected_power_w
                - result.viscous_heat_power_w,
            "resolved liquid differential first law",
        )?;
        if advective[..=n].iter().any(|v| !v.is_finite()) {
            return Err(invalid("represented liquid advective power"));
        }
        Ok(result)
    }
}

fn metric(surface: f64, trace: f64) -> Result<f64, VaporBubbleError> {
    ReservoirLayerGeometry::spherical_face(surface, trace)
        .map(|g| g.conductance_multiplier_m())
        .map_err(|_| invalid("original spherical liquid conduction geometry"))
}

struct Geometry {
    faces: [f64; MAX_FACES],
    centres: [f64; MAX_LIQUID_THERMAL_CELLS],
}

/// Prepared original material responses reused after the native phase solve.
pub struct LiquidThermalStage {
    state: VaporBubbleState,
    cells: usize,
    geometry: Geometry,
    reference_h: f64,
    temperature: [f64; MAX_LIQUID_THERMAL_CELLS],
    enthalpy: [f64; MAX_LIQUID_THERMAL_CELLS],
    capacity: [f64; MAX_LIQUID_THERMAL_CELLS],
    mass: [f64; MAX_LIQUID_THERMAL_CELLS],
    conductivity: [f64; MAX_LIQUID_THERMAL_CELLS],
    pressure_heating_factor: [f64; MAX_LIQUID_THERMAL_CELLS],
    compressibility: [f64; MAX_LIQUID_THERMAL_CELLS],
    saturation_pressure: [f64; MAX_LIQUID_THERMAL_CELLS],
    conductance: [f64; MAX_FACES],
    conductive_power: [f64; MAX_FACES],
    interface_h: f64,
    bath_h: f64,
    interface_capacity: f64,
    bath_capacity: f64,
    stored_energy_j: f64,
    maximum_source_relative_change: f64,
    maximum_pressure_path_capacity_fraction: f64,
}
impl LiquidThermalStage {
    pub fn interface_heat_power_w(&self) -> f64 {
        -self.conductive_power[0]
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
}

/// Instantaneous rates and separate physical boundary/storage accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidThermalRates {
    pub temperature_rate_k_s: [f64; MAX_LIQUID_THERMAL_CELLS],
    pub interface_heat_power_w: f64,
    /// Positive into vapor, using the fixed field reference enthalpy.
    pub interface_enthalpy_power_w: f64,
    /// Positive leaving the resolved liquid at its actual outer face.
    pub outer_conductive_power_w: f64,
    /// Positive leaving the resolved liquid, using the same enthalpy reference.
    pub outer_advected_power_w: f64,
    pub stored_energy_j: f64,
    pub energy_rate_residual_w: f64,
    pub maximum_outflow_peclet: f64,
    pub maximum_pressure_heating_temperature_k: f64,
    pub maximum_source_relative_change: f64,
    pub maximum_pressure_path_capacity_fraction: f64,
    /// Inverse local outgoing conduction/advection rate; not an error estimate.
    pub minimum_transport_time_s: f64,
    pub viscous_heat_power_w: f64,
    /// Native viscous loss beyond this finite thermal mesh, not deposited here.
    pub exterior_viscous_heat_power_w: f64,
    cells: usize,
    enthalpy: [f64; MAX_LIQUID_THERMAL_CELLS],
    mass: [f64; MAX_LIQUID_THERMAL_CELLS],
    reference_h: f64,
}
impl LiquidThermalRates {
    /// Factored original caloric change, preserving small sensible increments.
    pub fn energy_change_from(&self, initial: &Self) -> Result<f64, VaporBubbleError> {
        if self.cells != initial.cells
            || self.mass != initial.mass
            || self.reference_h != initial.reference_h
        {
            return Err(invalid(
                "same liquid cells, mass and caloric reference for energy change",
            ));
        }
        let mut value = 0.0;
        let mut correction = 0.0;
        for i in 0..self.cells {
            let term = self.mass[i] * (self.enthalpy[i] - initial.enthalpy[i]);
            let adjusted = term - correction;
            let next = value + adjusted;
            correction = (next - value) - adjusted;
            value = next;
        }
        finite(value, "original liquid caloric energy change")
    }
}

#[cfg(test)]
#[path = "liquid/tests.rs"]
mod tests;
