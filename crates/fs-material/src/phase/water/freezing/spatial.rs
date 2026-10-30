//! Spatial freezing on the original conservative thermal timeline.
//!
//! A stationary plane slab is divided into equal reference volumes. Its initial
//! source density fixes each cell's mass. Every later temperature, phase mass,
//! and phase volume comes from the original ice/liquid enthalpy chart; fusion
//! is never added as a separate heat source. The shared conduction observer
//! evaluates the selected phase conductivity every actual thermal substep.
//! Surface reservoirs include the boundary half-cell resistance.
//!
//! This is a one-dimensional conduction model on fixed REFERENCE geometry.
//! Pressure is externally maintained, and the reported free phase volumes and
//! p*Delta V do not move that geometry. Convection, mechanical compatibility,
//! nucleation, supercooling, phase-front kinetics and a deforming free surface
//! require additional models. Spatial and timestep refinement remain separate
//! from the numerical energy gates.

use super::thermal::FreezingWaterThermalLaw;
use super::transport::{
    FreezingConductivityError, FreezingWaterConductivity, FreezingWaterImplicitLaw,
};
use super::{FreezingWaterError, FreezingWaterState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionError,
    ImplicitConductionLaw, ImplicitConductionReport, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{
    ConductionObserver, ConductiveLink, ConductiveReservoir, ConductivityLaw,
};
use fs_thermal::enthalpy::entropy::EntropyAudit;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, EnergyAudit, Event, MaterialId, Reservoir, StepConfig, ThermalLink,
    TimedEvent,
};

/// Plane reference geometry. The two transverse directions are adiabatic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabGeometry {
    pub thickness_m: f64,
    pub area_m2: f64,
    /// At least two cells; the two face reservoirs bind distinct boundary cells.
    pub cells: usize,
}

/// A face bath, at its actual reservoir temperature rather than cell temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabReservoir {
    pub temperature_k: f64,
    /// Additional face resistance [m² K/W]. Zero imposes temperature at the
    /// surface; 1/htc gives a Robin bath. The half-cell resistance is added once.
    pub surface_resistance_m2_k_w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlabFace {
    Left,
    Right,
}

/// Change an installed face bath without changing its area or resistance.
/// Events are supplied in nondecreasing time order and hit by the original
/// thermal event owner. A face declared insulated cannot acquire a bath here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabTemperatureEvent {
    pub time_s: f64,
    pub face: SlabFace,
    pub temperature_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreezingSlabError {
    Thermal(DepositionError),
    Material(FreezingWaterError),
    Transport(FreezingConductivityError),
    Implicit(ImplicitConductionError),
    Invalid(&'static str),
    EnergyBalance {
        maximum_residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<DepositionError> for FreezingSlabError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<FreezingWaterError> for FreezingSlabError {
    fn from(error: FreezingWaterError) -> Self {
        Self::Material(error)
    }
}
impl From<FreezingConductivityError> for FreezingSlabError {
    fn from(error: FreezingConductivityError) -> Self {
        Self::Transport(error)
    }
}
impl From<ImplicitConductionError> for FreezingSlabError {
    fn from(error: ImplicitConductionError) -> Self {
        Self::Implicit(error)
    }
}
impl core::fmt::Display for FreezingSlabError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial water freezing: {self:?}")
    }
}
impl std::error::Error for FreezingSlabError {}

/// Resolved state of one reference control volume, in left-to-right order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabCell {
    pub cell_id: CellId,
    pub reference_centre_m: f64,
    pub reference_mass_kg: f64,
    /// Actual timeline coordinate; the source inverse has its own small tolerance.
    pub stored_specific_enthalpy_j_kg: f64,
    pub material: FreezingWaterState,
}

/// Current inventory. These free volumes are not the frozen transport geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabInventory {
    pub mass_kg: f64,
    pub ice_mass_kg: f64,
    pub liquid_mass_kg: f64,
    pub free_ice_volume_m3: f64,
    pub free_liquid_volume_m3: f64,
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Thickness times ice MASS fraction, because reference density is uniform.
    /// A phase-inventory equivalent length, not an identified sharp interface or
    /// the physical ice thickness. Inspect the cell profile to locate a front.
    pub equivalent_frozen_reference_length_m: f64,
}
impl FreezingSlabInventory {
    pub fn free_volume_m3(self) -> f64 {
        self.free_ice_volume_m3 + self.free_liquid_volume_m3
    }
}

/// Independently checked source calorimetry for one accepted advance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabAudit {
    pub thermal: EnergyAudit,
    pub entropy: EntropyAudit,
    /// Difference of original, unshifted source enthalpies, summed by cell.
    pub source_enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub pressure_work_j: f64,
    pub source_enthalpy_residual_j: f64,
    /// Delta U + p*Delta V - Q, independently evaluated from source states.
    pub internal_energy_residual_j: f64,
}

/// Independent material calorimetry and original implicit solver diagnostics.
/// Physical energy appears in both nested reports for comparison, not addition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingSlabImplicitAudit {
    pub physical: FreezingSlabAudit,
    pub solve: ImplicitConductionReport,
}

/// A reusable physical slab assembled from the existing caloric and transport
/// owners. All advances return a candidate: failure leaves this state untouched.
#[derive(Clone)]
pub struct FreezingSlab<'a> {
    thermal: DepositionSimulation<'a>,
    conduction: ConductionObserver<'a>,
    caloric: &'a FreezingWaterThermalLaw,
    transport: &'a FreezingWaterConductivity,
    geometry: FreezingSlabGeometry,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn bath_boundary(area: f64, temperature: Option<f64>) -> Boundary {
    Boundary {
        ambient_k: 273.15,
        htc_w_m2_k: 0.0,
        // The geometric conduction observer replaces this authored zero before
        // heat rates and CFL. No artificial large surface coefficient is used.
        reservoir: temperature.map(|temperature_k| Reservoir {
            area_m2: area,
            conductance_w_k: 0.0,
            temperature_k,
        }),
        source_w: 0.0,
    }
}

impl<'a> FreezingSlab<'a> {
    /// Create and activate an initially uniform, closed-mass reference slab.
    /// Density and mass are derived from the supplied original source state.
    /// Both material bindings must retain exactly the same pressure chart.
    pub fn new(
        caloric: &'a FreezingWaterThermalLaw,
        transport: &'a FreezingWaterConductivity,
        geometry: FreezingSlabGeometry,
        initial: FreezingWaterState,
        reservoirs: [Option<FreezingSlabReservoir>; 2],
        temperature_events: &[FreezingSlabTemperatureEvent],
        budget: DepositionBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, FreezingSlabError> {
        checkpoint()?;
        let n = geometry.cells;
        let event_count = n
            .checked_mul(2)
            .and_then(|v| v.checked_sub(1))
            .and_then(|v| v.checked_add(temperature_events.len()));
        if n < 2
            || n > budget.max_cells
            || n - 1 > budget.max_links
            || event_count.is_none_or(|v| v > budget.max_events)
            || !positive(geometry.thickness_m)
            || !positive(geometry.area_m2)
            || initial.law() != caloric.chart()
            || transport.chart() != caloric.chart()
        {
            return Err(FreezingSlabError::Invalid(
                "positive plane geometry, bounded 2+ cells and identical caloric/transport charts required",
            ));
        }
        for bath in reservoirs.into_iter().flatten() {
            if !positive(bath.temperature_k) || !nonnegative(bath.surface_resistance_m2_k_w) {
                return Err(FreezingSlabError::Invalid(
                    "finite positive bath temperature and nonnegative face resistance required",
                ));
            }
        }
        let dx = geometry.thickness_m / n as f64;
        let volume = geometry.area_m2 * dx;
        let mass = volume / initial.specific_volume_m3_kg();
        if !positive(dx) || !positive(volume) || !positive(mass) {
            return Err(FreezingSlabError::Invalid(
                "unrepresentable cell length, volume or source mass",
            ));
        }
        let mut events = Vec::with_capacity(event_count.expect("checked event count"));
        for i in 0..n {
            checkpoint()?;
            let bath = if i == 0 {
                reservoirs[0]
            } else if i + 1 == n {
                reservoirs[1]
            } else {
                None
            };
            let start = dx * i as f64;
            let end = if i + 1 == n {
                geometry.thickness_m
            } else {
                dx * (i + 1) as f64
            };
            if !(end > start) {
                return Err(FreezingSlabError::Invalid(
                    "cell coordinates are not representably ordered",
                ));
            }
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: 0,
                    body_id: 0,
                    material_id: MaterialId(0),
                    mass_kg: mass,
                    specific_h_j_kg: initial.specific_enthalpy_j_kg(),
                    geometry: CellGeometry {
                        start_m: [start, 0.0, 0.0],
                        end_m: [end, 0.0, 0.0],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * geometry.area_m2,
                    },
                    boundary: bath_boundary(geometry.area_m2, bath.map(|b| b.temperature_k)),
                }),
            });
        }
        let mut links = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            checkpoint()?;
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: CellId(i as u64),
                    b: CellId((i + 1) as u64),
                    area_m2: geometry.area_m2,
                    conductance_w_k: 0.0,
                }),
            });
            links.push(ConductiveLink {
                link_id: i as u64,
                half_distances_m: [0.5 * dx; 2],
                contact_resistance_m2_k_w: 0.0,
            });
        }
        for event in temperature_events {
            checkpoint()?;
            let side = usize::from(event.face == SlabFace::Right);
            if reservoirs[side].is_none() || !positive(event.temperature_k) {
                return Err(FreezingSlabError::Invalid(
                    "temperature event needs an installed face bath and positive temperature",
                ));
            }
            events.push(TimedEvent {
                time_s: event.time_s,
                event: Event::Boundary(
                    CellId(if side == 0 { 0 } else { (n - 1) as u64 }),
                    bath_boundary(geometry.area_m2, Some(event.temperature_k)),
                ),
            });
        }
        let thermal = DepositionSimulation::new(events, &[(MaterialId(0), caloric)], budget)?;
        let (thermal, _) = thermal.advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            &mut checkpoint,
        )?;
        let faces: Vec<_> = reservoirs
            .iter()
            .enumerate()
            .filter_map(|(side, bath)| {
                bath.map(|bath| ConductiveReservoir {
                    cell_id: CellId(if side == 0 { 0 } else { (n - 1) as u64 }),
                    half_distance_m: 0.5 * dx,
                    surface_resistance_m2_k_w: bath.surface_resistance_m2_k_w,
                })
            })
            .collect();
        let conduction = ConductionObserver::new(
            &thermal,
            &[(MaterialId(0), transport as &dyn ConductivityLaw)],
            &links,
            &faces,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(Self {
            thermal,
            conduction,
            caloric,
            transport,
            geometry,
        })
    }

    pub fn time_s(&self) -> f64 {
        self.thermal.time_s()
    }
    pub fn mass_kg(&self) -> f64 {
        self.thermal.mass_kg()
    }
    pub const fn geometry(&self) -> FreezingSlabGeometry {
        self.geometry
    }
    /// Original accepted thermal owner, including its source enthalpies and IDs.
    pub const fn thermal(&self) -> &DepositionSimulation<'a> {
        &self.thermal
    }

    /// Resolve the current material field without modifying physical history.
    pub fn profile(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreezingSlabCell>, FreezingSlabError> {
        checkpoint()?;
        let mut cells = Vec::with_capacity(self.geometry.cells);
        for cell in self.thermal.cells().values() {
            checkpoint()?;
            let material = self
                .caloric
                .chart()
                .state_at_specific_enthalpy(cell.specific_h_j_kg, || false)?
                .into_state();
            cells.push(FreezingSlabCell {
                cell_id: cell.id,
                reference_centre_m: cell.geometry.start_m[0]
                    + 0.5 * (cell.geometry.end_m[0] - cell.geometry.start_m[0]),
                reference_mass_kg: cell.mass_kg,
                stored_specific_enthalpy_j_kg: cell.specific_h_j_kg,
                material,
            });
        }
        checkpoint()?;
        Ok(cells)
    }

    pub fn inventory(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FreezingSlabInventory, FreezingSlabError> {
        let cells = self.profile(&mut checkpoint)?;
        let mut mass = Sum::default();
        let mut ice_mass = Sum::default();
        let mut liquid_mass = Sum::default();
        let mut ice_volume = Sum::default();
        let mut liquid_volume = Sum::default();
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        for cell in cells {
            checkpoint()?;
            let state = cell.material;
            mass.add(cell.reference_mass_kg);
            ice_mass.add(cell.reference_mass_kg * state.ice_mass_fraction());
            liquid_mass.add(cell.reference_mass_kg * state.liquid_mass_fraction());
            let volume = cell.reference_mass_kg * state.specific_volume_m3_kg();
            ice_volume.add(volume * state.ice_volume_fraction());
            liquid_volume.add(volume * state.liquid_volume_fraction());
            minimum = minimum.min(state.temperature_k());
            maximum = maximum.max(state.temperature_k());
        }
        let mass = mass.total()?;
        let ice_mass = ice_mass.total()?;
        checkpoint()?;
        Ok(FreezingSlabInventory {
            mass_kg: mass,
            ice_mass_kg: ice_mass,
            liquid_mass_kg: liquid_mass.total()?,
            free_ice_volume_m3: ice_volume.total()?,
            free_liquid_volume_m3: liquid_volume.total()?,
            minimum_temperature_k: minimum,
            maximum_temperature_k: maximum,
            equivalent_frozen_reference_length_m: self.geometry.thickness_m * (ice_mass / mass),
        })
    }

    /// Advance on the shared CFL/event grid and independently check original
    /// source Delta H = Q and Delta U + p*Delta V = Q. The same joule tolerance
    /// applies to all three accounts; it is not an accuracy or model-error bound.
    pub fn advance_to(
        &self,
        target_s: f64,
        controls: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FreezingSlabAudit), FreezingSlabError> {
        checkpoint()?;
        let (thermal, conduction, energy, entropy) = self.thermal.advance_observed_with_entropy(
            &self.conduction,
            target_s,
            controls,
            &mut checkpoint,
        )?;
        self.finish_advance(
            thermal,
            conduction,
            energy,
            entropy,
            controls.energy_tolerance_j,
            checkpoint,
        )
    }

    /// Advance the SAME spatial material with endpoint-evaluated implicit heat
    /// transport. Analytic source Cp and phase-volume conductivity tangents
    /// retain the exact fusion plateau. The original thermal residual, local
    /// and global joule gates, entropy and event owner determine acceptance.
    ///
    /// The reusable workspace has O(cells) storage and a pivoted tridiagonal
    /// linear solve. Maximum timestep and refinement remain caller choices:
    /// implicit stepping removes the explicit diffusion CFL restriction, not
    /// time-discretization error or the fixed-reference-geometry assumptions.
    /// Failure can change numerical scratch only; this physical slab is intact.
    pub fn advance_implicit_to(
        &self,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FreezingSlabImplicitAudit), FreezingSlabError> {
        checkpoint()?;
        let law = FreezingWaterImplicitLaw::new(self.caloric, self.transport, || false)?;
        let chain = ImplicitConductionChain::new(
            &self.thermal,
            &self.conduction,
            &[(MaterialId(0), &law as &dyn ImplicitConductionLaw)],
            &mut checkpoint,
        )?;
        let (thermal, solve) = chain.advance_to(
            &self.thermal,
            target_s,
            controls,
            workspace,
            &mut checkpoint,
        )?;
        let (next, physical) = self.finish_advance(
            thermal,
            self.conduction.clone(),
            solve.energy,
            solve.entropy,
            controls.total_tolerance_j,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok((next, FreezingSlabImplicitAudit { physical, solve }))
    }

    fn finish_advance(
        &self,
        thermal: DepositionSimulation<'a>,
        conduction: ConductionObserver<'a>,
        energy: EnergyAudit,
        entropy: EntropyAudit,
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FreezingSlabAudit), FreezingSlabError> {
        let initial = self.profile(&mut checkpoint)?;
        let next = Self {
            thermal,
            conduction,
            caloric: self.caloric,
            transport: self.transport,
            geometry: self.geometry,
        };
        let final_cells = next.profile(&mut checkpoint)?;
        let mut h = Sum::default();
        let mut u = Sum::default();
        let mut work = Sum::default();
        for (a, b) in initial.iter().zip(&final_cells) {
            checkpoint()?;
            let m = a.reference_mass_kg;
            h.add(m * (b.material.specific_enthalpy_j_kg() - a.material.specific_enthalpy_j_kg()));
            u.add(
                m * (b.material.specific_internal_energy_j_kg()
                    - a.material.specific_internal_energy_j_kg()),
            );
            work.add(
                m * self.caloric.chart().pressure_pa()
                    * (b.material.specific_volume_m3_kg() - a.material.specific_volume_m3_kg()),
            );
        }
        let h = h.total()?;
        let u = u.total()?;
        let work = work.total()?;
        let dh = h - energy.external_heat_j;
        let du = u + work - energy.external_heat_j;
        let residual = energy.residual_j.abs().max(dh.abs()).max(du.abs());
        if !residual.is_finite() || residual > energy_tolerance_j {
            return Err(FreezingSlabError::EnergyBalance {
                maximum_residual_j: residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        checkpoint()?;
        Ok((
            next,
            FreezingSlabAudit {
                thermal: energy,
                entropy,
                source_enthalpy_change_j: h,
                internal_energy_change_j: u,
                pressure_work_j: work,
                source_enthalpy_residual_j: dh,
                internal_energy_residual_j: du,
            },
        ))
    }
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn total(self) -> Result<f64, FreezingSlabError> {
        let value = self.sum + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DepositionError::Arithmetic.into())
        }
    }
}

#[cfg(test)]
mod tests;
