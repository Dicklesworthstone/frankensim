//! Read-only flash-temperature comparison against resolved thermal motion.
//!
//! Accepted, conduction-derived body heat supplies piecewise-constant surface
//! fluxes. Their on/off responses are superposed through the ORIGINAL
//! `fs_tribo::flash_temperature_candidate` kernel. Signed transfers and cooling
//! after a pulse are retained. This observer neither routes heat nor creates a
//! thermal store, and its estimate never replaces an accepted network state.
//!
//! The comparison is the stationary, one-dimensional, constant-property,
//! semi-infinite uniform-flux limit. Each body must initially be uniform; the
//! original fixed geometry/materials and source thermal properties are retained.
//! Finite depth, lateral diffusion, changing properties, third-body storage and
//! a moving footprint can invalidate this limit. In particular this is NOT a
//! two-percent accuracy claim for arbitrary moving-source Peclet number.
//! See the classical sliding-contact context in J. C. Jaeger (1943),
//! "Moving sources of heat and the temperature at sliding contacts",
//! https://doi.org/10.5962/p.360338. Spatial and temporal convergence of the
//! actual network must establish the comparison's useful domain.

use super::{FrictionHeatReport, FrictionHeatRoute, finite, invalid};
use fs_math::det;
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, MaterialId, ThermalLink,
};
use fs_tribo::{
    FlashTemperatureEstimate, FlashTemperatureInput, HeatPartition, InterfaceSystemRef,
    SurfaceThermalProperties, TriboError, flash_temperature_candidate,
};
use std::collections::BTreeMap;

#[derive(Debug)]
pub enum FlashOracleError {
    Thermal(DepositionError),
    Tribology(TriboError),
    HistoryBudget { maximum: usize },
}
impl From<DepositionError> for FlashOracleError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<TriboError> for FlashOracleError {
    fn from(error: TriboError) -> Self {
        Self::Tribology(error)
    }
}
impl core::fmt::Display for FlashOracleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "flash comparison: {self:?}")
    }
}
impl std::error::Error for FlashOracleError {}

#[derive(Clone, Debug)]
struct CellBinding {
    id: CellId,
    body: u64,
    material: MaterialId,
    geometry: CellGeometry,
    mass: f64,
}
#[derive(Clone, Copy, Debug)]
struct Pulse {
    start_s: f64,
    end_s: f64,
    power_w: [f64; 2],
}

/// Diagnostic error at the same time as an accepted resolved interface state.
#[derive(Clone, Debug)]
pub struct FlashTemperatureComparison {
    pub time_s: f64,
    pub resolved_interface_temperature_k: f64,
    pub oracle_surface_temperature_k: [f64; 2],
    pub absolute_error_k: [f64; 2],
    /// Relative to each oracle rise from its initial uniform temperature.
    /// None denotes exactly zero predicted rise; no denominator floor is added.
    pub relative_rise_error: [Option<f64>; 2],
    pub retained_pulses: usize,
}

/// Bounded flux-history observer, not an energy or temperature state owner.
/// Clone preserves its exact comparison history for checkpoint/replay.
#[derive(Clone, Debug)]
pub struct FlashHeatHistory {
    interface: InterfaceSystemRef,
    route: FrictionHeatRoute,
    properties: [SurfaceThermalProperties; 2],
    bindings: Vec<CellBinding>,
    links: BTreeMap<u64, ThermalLink>,
    initial_temperature_k: [f64; 2],
    time_s: f64,
    pulses: Vec<Pulse>,
    maximum_pulses: usize,
    energy_tolerance_j: f64,
}

impl FlashHeatHistory {
    /// Bind the two original, initially uniform bodies and constant material
    /// thermal properties. The route supplies actual area; there is no caller
    /// heat fraction. Bounds are explicit, and the oracle retains caller/source
    /// authority of k and alpha rather than inferring material calibration.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &DepositionSimulation<'_>,
        route: &FrictionHeatRoute,
        interface: InterfaceSystemRef,
        properties: [SurfaceThermalProperties; 2],
        maximum_cells: usize,
        maximum_pulses: usize,
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, FlashOracleError> {
        checkpoint()?;
        if maximum_cells == 0
            || initial.cells().len() > maximum_cells
            || maximum_pulses == 0
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(invalid(
                "flash comparison needs explicit finite cell/history/energy budgets",
            )
            .into());
        }
        route.validate(initial)?;
        FrictionHeatRoute::from_link(
            initial,
            route.patch.link.id,
            properties.map(|p| p.conductivity_w_per_m_k()),
        )?;
        let ids = route.cells();
        let body = ids.map(|id| initial.cells()[&id].body_id);
        if body[0] == body[1] {
            return Err(invalid("half-space comparison needs two distinct bodies").into());
        }
        let initial_temperature_k = [
            initial.temperature_k(ids[0])?,
            initial.temperature_k(ids[1])?,
        ];
        let mut bindings = Vec::with_capacity(initial.cells().len());
        for (i, (&id, cell)) in initial.cells().iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let side = body.iter().position(|&b| b == cell.body_id).ok_or(invalid(
                "flash comparison network contains an unrelated body",
            ))?;
            if initial.temperature_k(id)? != initial_temperature_k[side] {
                return Err(invalid(
                    "half-space comparison requires initially uniform body temperatures",
                )
                .into());
            }
            bindings.push(CellBinding {
                id,
                body: cell.body_id,
                material: cell.material_id,
                geometry: cell.geometry,
                mass: cell.mass_kg,
            });
        }
        let observer = Self {
            interface,
            route: route.clone(),
            properties,
            bindings,
            links: initial.links().clone(),
            initial_temperature_k,
            time_s: initial.time_s(),
            pulses: Vec::new(),
            maximum_pulses,
            energy_tolerance_j,
        };
        observer.validate(initial, &mut checkpoint)?;
        // Invoke the original owner even for zero power, so interface and input
        // admission are identical to a nonzero flash-temperature evaluation.
        observer.response(
            Pulse {
                start_s: 0.0,
                end_s: 1.0,
                power_w: [0.0; 2],
            },
            1.0,
        )?;
        checkpoint()?;
        Ok(observer)
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn retained_pulses(&self) -> usize {
        self.pulses.len()
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.route.validate(thermal)?;
        if thermal.cells().len() != self.bindings.len() || thermal.links() != &self.links {
            return Err(invalid(
                "flash comparison changed its original topology or conductance sources",
            ));
        }
        for (i, binding) in self.bindings.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let cell = thermal
                .cells()
                .get(&binding.id)
                .ok_or(invalid("inactive flash body cell"))?;
            if cell.body_id != binding.body
                || cell.material_id != binding.material
                || cell.mass_kg != binding.mass
                || cell.geometry != binding.geometry
            {
                return Err(invalid(
                    "flash comparison changed original body/material/geometry/mass",
                ));
            }
            if cell.boundary.source_w != 0.0
                || cell.boundary.htc_w_m2_k != 0.0
                || cell.boundary.reservoir.is_some()
            {
                return Err(invalid(
                    "half-space pulse comparison excludes other boundary heat",
                ));
            }
        }
        for (i, link) in thermal.links().values().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if thermal.cells()[&link.a].body_id != thermal.cells()[&link.b].body_id
                && link.id != self.route.patch.link.id
            {
                return Err(invalid(
                    "flash comparison has an unobserved body-to-body heat path",
                ));
            }
        }
        Ok(())
    }

    fn response(&self, pulse: Pulse, target_s: f64) -> Result<[f64; 2], FlashOracleError> {
        let duration = finite(pulse.end_s - pulse.start_s)?;
        let age = finite(target_s - pulse.start_s)?;
        let cool = finite(target_s - pulse.end_s)?;
        if duration <= 0.0 || age < duration || cool < 0.0 {
            return Err(invalid(
                "flash pulse must have a positive resolved duration before observation",
            )
            .into());
        }
        let total = finite(pulse.power_w[0].abs() + pulse.power_w[1].abs())?;
        // These are DERIVED magnitudes of measured body flux. Signed directions
        // are restored below; this partition does not route production heat.
        let partition = if total == 0.0 {
            HeatPartition::new(0.0, 0.0, 1.0)?
        } else {
            HeatPartition::new(
                pulse.power_w[0].abs() / total,
                pulse.power_w[1].abs() / total,
                0.0,
            )?
        };
        let input = FlashTemperatureInput::for_duration(
            total,
            partition,
            self.route.patch.link.area_m2,
            duration,
            Some(self.properties[0]),
            Some(self.properties[1]),
        )?;
        let FlashTemperatureEstimate::Candidate(candidate) =
            flash_temperature_candidate(&self.interface, input)?
        else {
            return Err(invalid("flash source has no thermal-property closure").into());
        };
        // Stable on/off superposition: (sqrt(age)-sqrt(cool))/sqrt(duration)
        // = sqrt(duration/age)/(1+sqrt(cool/age)). No subtraction of nearly
        // equal old pulses, nor evaluation at an artificial traverse speed.
        let weight = finite(det::sqrt(duration / age) / (1.0 + det::sqrt(cool / age)))?;
        Ok([
            finite(pulse.power_w[0].signum() * candidate.surface_a_rise_k() * weight)?,
            finite(pulse.power_w[1].signum() * candidate.surface_b_rise_k() * weight)?,
        ])
    }

    /// Observe one accepted thermal interval without modifying either thermal
    /// endpoint. The physical receipt is independently checked against the
    /// actual cell storage and endpoint interface conduction before retention.
    /// A refusal leaves this observer's time/history unchanged as well.
    pub fn observe(
        &mut self,
        before: &DepositionSimulation<'_>,
        after: &DepositionSimulation<'_>,
        report: &FrictionHeatReport,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FlashTemperatureComparison, FlashOracleError> {
        checkpoint()?;
        self.validate(before, &mut checkpoint)?;
        self.validate(after, &mut checkpoint)?;
        let start = before.time_s();
        let end = after.time_s();
        let duration = finite(end - start)?;
        if start != self.time_s
            || duration <= 0.0
            || !report.interface_temperature_k.is_finite()
            || report.interface_temperature_k <= 0.0
        {
            return Err(
                invalid("flash observations need contiguous accepted physical intervals").into(),
            );
        }
        let check = |a: f64, b: f64| -> Result<(), DepositionError> {
            if finite(a - b)?.abs() > self.energy_tolerance_j {
                Err(invalid("flash receipt differs from accepted thermal work"))
            } else {
                Ok(())
            }
        };
        check(report.external_heat_j, 0.0)?;
        check(report.interface_residual_j, 0.0)?;
        check(report.total_energy_residual_j, 0.0)?;
        let mut stored = 0.0;
        for (i, binding) in self.bindings.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            stored = finite(
                stored
                    + finite(
                        binding.mass
                            * finite(
                                after.cells()[&binding.id].specific_h_j_kg
                                    - before.cells()[&binding.id].specific_h_j_kg,
                            )?,
                    )?,
            )?;
        }
        check(stored, report.thermal.stored_energy_change_j)?;
        check(stored, report.dissipated_j)?;
        let expected_g = self
            .route
            .conductances(self.properties.map(|p| p.conductivity_w_per_m_k()))?;
        let mut power_w = [0.0; 2];
        for side in 0..2 {
            if !report.conductances_w_k[side].is_finite()
                || (report.conductances_w_k[side] - expected_g[side]).abs()
                    > 64.0 * f64::EPSILON * expected_g[side].abs()
            {
                return Err(invalid(
                    "flash oracle conductivity differs from the accepted material source",
                )
                .into());
            }
            let actual_heat = finite(
                duration
                    * finite(
                        expected_g[side]
                            * (report.interface_temperature_k
                                - after.temperature_k(self.route.cells()[side])?),
                    )?,
            )?;
            check(actual_heat, report.body_heat_j[side])?;
            power_w[side] = finite(report.body_heat_j[side] / duration)?;
        }
        check(
            finite(report.body_heat_j[0] + report.body_heat_j[1])?,
            report.dissipated_j,
        )?;
        let mut pulse = Pulse {
            start_s: start,
            end_s: end,
            power_w,
        };
        let active = power_w != [0.0; 2];
        let merge = active
            && self
                .pulses
                .last()
                .is_some_and(|last| last.end_s == start && last.power_w == power_w);
        if merge {
            pulse.start_s = self.pulses.last().expect("checked last pulse").start_s;
        }
        if active && !merge && self.pulses.len() >= self.maximum_pulses {
            return Err(FlashOracleError::HistoryBudget {
                maximum: self.maximum_pulses,
            });
        }
        let mut rise = [0.0; 2];
        for (i, previous) in self.pulses.iter().copied().enumerate() {
            if i % 64 == 0 {
                checkpoint()?;
            }
            if merge && i + 1 == self.pulses.len() {
                continue;
            }
            let response = self.response(previous, end)?;
            for side in 0..2 {
                rise[side] = finite(rise[side] + response[side])?;
            }
        }
        if active {
            let response = self.response(pulse, end)?;
            for side in 0..2 {
                rise[side] = finite(rise[side] + response[side])?;
            }
        }
        let mut temperature = [0.0; 2];
        let mut absolute_error = [0.0; 2];
        let mut relative_error = [None; 2];
        for side in 0..2 {
            temperature[side] = finite(self.initial_temperature_k[side] + rise[side])?;
            if temperature[side] <= 0.0 {
                return Err(
                    invalid("flash comparison left its positive-temperature domain").into(),
                );
            }
            absolute_error[side] =
                finite(report.interface_temperature_k - temperature[side])?.abs();
            if rise[side] != 0.0 {
                relative_error[side] = Some(finite(absolute_error[side] / rise[side].abs())?);
            }
        }
        let retained = self.pulses.len() + usize::from(active && !merge);
        let comparison = FlashTemperatureComparison {
            time_s: end,
            resolved_interface_temperature_k: report.interface_temperature_k,
            oracle_surface_temperature_k: temperature,
            absolute_error_k: absolute_error,
            relative_rise_error: relative_error,
            retained_pulses: retained,
        };
        checkpoint()?;
        if merge {
            *self.pulses.last_mut().expect("checked last pulse") = pulse;
        } else if active {
            self.pulses.push(pulse);
        }
        self.time_s = end;
        Ok(comparison)
    }
}
