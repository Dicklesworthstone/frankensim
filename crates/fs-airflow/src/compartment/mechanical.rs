//! Pressure-dependent ventilation from the original fan and loss owners.
//!
//! `FanCurve` stores AMCA fan *static* pressure. Matching it directly to a
//! total-pressure duct loss would omit the fan outlet velocity pressure. This
//! adapter retains the original bank curve and adds that velocity pressure from
//! the declared per-fan outlet area before matching the original `LossNetwork`.
//! Fan pressure and quadratic loss both scale with actual donor density relative
//! to their common reference density. At fixed speed this leaves volume flow
//! unchanged when the complete pressure load scales with density, while mass
//! flow and shaft power change. The source is ASHRAE Handbook 2016, HVAC Systems
//! and Equipment, chapter 21, sections 5 and 10 and the nomenclature:
//! <https://handbook.ashrae.org/Handbooks/S16/SI/s16_ch21/s16_ch21_si.aspx>.
//!
//! The operating point is a bounded nominal estimate. Supplied curve/loss
//! tolerances and validity cards remain attached to their original owners;
//! this adapter does not turn them into physical uncertainty bounds. Density
//! similarity, low pressure change, frozen inlet state, negligible endpoint
//! kinetic storage and equal-elevation ports are explicit approximations.
//! Reverse flow, fan stall, surge, duct storage and fan-off windmilling refuse
//! or require another model. Closing a connection explicitly stops its fan and
//! closes an ideal damper; it is not inferred from a failed operating point.

use super::{Endpoint, Layer, NetworkError, Room, ZoneNetwork, poll, positive, vent_layer};
use crate::{AirflowError, FanArrangement, FanBank, LossNetwork, SourceProvenance};
use crate::graph::thermal::medium::fan::density_affinity_ratio;
use fs_qty::{Area, Density, Length, MassFlowRate, Power, Pressure, VolumetricFlowRate};

/// Bound on original fan samples admitted by this low-order connection.
pub const MAX_FAN_POINTS: usize = 256;
/// Bound on original series/parallel loss nodes admitted at construction.
pub const MAX_LOSS_NODES: usize = 512;

/// Explicit applicability of the incompressible density-affinity approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanDomain {
    /// Actual donor density / reference density, inclusive positive bounds.
    pub density_ratio: [f64; 2],
    /// Maximum fan rise, duct drop or endpoint difference / minimum absolute
    /// endpoint pressure. Must lie in (0, 0.1]; this is a caller's applicability
    /// limit, not a certified compressibility-error bound.
    pub maximum_pressure_fraction: f64,
}

/// Sourced power conversion for the complete fan bank at its operating point.
#[derive(Clone, Debug, PartialEq)]
pub struct FanEnergy {
    /// Total aerodynamic power / fan shaft input power, in (0, 1]. This is
    /// total efficiency, not static or motor electrical efficiency.
    pub total_efficiency: f64,
    /// Fraction of shaft-minus-aerodynamic power rejected to surroundings,
    /// in [0, 1]. The rest remains in the gas. Zero gives an adiabatic fan/duct
    /// assembly with all shaft work eventually retained by the gas.
    pub rejected_loss_fraction: f64,
    /// Authority for efficiency and heat-rejection assumptions.
    pub source: SourceProvenance,
}

/// Bounded operating-point numerical work, separate from physical uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanControl {
    /// At most 128 bracket refinements.
    pub maximum_iterations: usize,
    /// Absolute fan-minus-system pressure residual [Pa].
    pub pressure_tolerance_pa: f64,
    /// Absolute width of the retained volume-flow bracket [m^3/s].
    pub flow_tolerance_m3_s: f64,
}
impl Default for FanControl {
    fn default() -> Self {
        Self {
            maximum_iterations: 64,
            pressure_tolerance_pa: 1e-7,
            flow_tolerance_m3_s: 1e-12,
        }
    }
}

/// An original fan bank and total-pressure loss network at one declared density.
#[derive(Clone, Debug, PartialEq)]
pub struct FanModel {
    bank: FanBank,
    losses: LossNetwork,
    reference_density: Density,
    outlet_area_per_fan: Area,
    energy: FanEnergy,
    domain: FanDomain,
    resistance: f64,
    velocity_resistance: f64,
    flow_domain: [f64; 2],
}

/// Actual nominal operating point; no clipping or curve extrapolation occurs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanOperatingPoint {
    /// Actual donor-side volume flow through the complete bank.
    pub volume_flow: VolumetricFlowRate,
    /// Actual donor density times volume flow.
    pub mass_rate: MassFlowRate,
    /// Original density-scaled static pressure plus outlet velocity pressure.
    pub total_pressure_rise: Pressure,
    /// Density-scaled total-pressure loss of the original network.
    pub duct_pressure_drop: Pressure,
    /// Fan shaft input; electrical motor loss is outside this declaration.
    pub shaft_power: Power,
    /// Shaft power retained as gas total enthalpy after declared heat rejection.
    pub gas_power: Power,
    /// Shaft loss rejected directly to the surroundings.
    pub rejected_power: Power,
    /// Fan rise minus duct loss minus receiver/donor pressure difference [Pa].
    pub pressure_residual_pa: f64,
    /// Actual refinements used by the bounded root solve.
    pub iterations: usize,
}

/// Original owner, applicability or bounded operating-point refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FanError {
    /// Invalid geometry, model assumptions, controls or nonfinite arithmetic.
    Invalid(&'static str),
    /// Original fan-curve evaluation refused its source domain.
    Airflow(AirflowError),
    /// Available forward curve cannot support the endpoint pressure difference.
    NoOperatingPoint {
        /// Fan-minus-system residual at the admitted low flow [Pa].
        low_residual_pa: f64,
        /// Fan-minus-system residual at the admitted high flow [Pa].
        high_residual_pa: f64,
    },
    /// Numerical bracket/residual gates did not both close within the budget.
    IterationLimit,
    /// Cancellation before a complete result; no owner was modified.
    Cancelled,
}
impl core::fmt::Display for FanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compartment fan: {self:?}")
    }
}
impl std::error::Error for FanError {}
impl From<AirflowError> for FanError {
    fn from(value: AirflowError) -> Self {
        Self::Airflow(value)
    }
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), FanError> {
    if cancelled() {
        Err(FanError::Cancelled)
    } else {
        Ok(())
    }
}
fn source_present(source: &SourceProvenance) -> bool {
    !source.citation.trim().is_empty() && !source.identifier.trim().is_empty()
}
fn finite(value: f64) -> Result<f64, FanError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FanError::Invalid("nonfinite fan arithmetic"))
    }
}
fn admit_loss(
    node: &LossNetwork,
    depth: usize,
    count: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FanError> {
    check(cancelled)?;
    *count += 1;
    if depth > 16 || *count > MAX_LOSS_NODES {
        return Err(FanError::Invalid("fan loss topology budget"));
    }
    match node {
        LossNetwork::Element(e) => {
            if e.name.trim().is_empty()
                || !positive(e.resistance.value())
                || !(0.0..1.0).contains(&e.uncertainty_rel)
                || !source_present(&e.source)
            {
                return Err(FanError::Invalid("fan loss declaration"));
            }
        }
        LossNetwork::Series(children) | LossNetwork::Parallel(children) => {
            if children.is_empty() {
                return Err(FanError::Invalid("empty fan loss group"));
            }
            for child in children {
                admit_loss(child, depth + 1, count, cancelled)?;
            }
        }
    }
    Ok(())
}

impl FanModel {
    /// Retain the exact source owners and precompute only immutable geometric
    /// factors. The original static curve plus velocity correction minus the
    /// original loss must decrease strictly over every admitted curve segment;
    /// unstable or multiple-root possibilities require a different fan model.
    /// `losses` must include the inlet, outlet/mixing and installation total-
    /// pressure losses appropriate to the retained fan's rating configuration.
    /// Those losses are not inferred from the supplied outlet area.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bank: FanBank,
        losses: LossNetwork,
        reference_density: Density,
        outlet_area_per_fan: Area,
        energy: FanEnergy,
        domain: FanDomain,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FanError> {
        check(&mut cancelled)?;
        if bank.has_material_density() {
            return Err(FanError::Invalid(
                "evolving compartment fan requires the original reference bank, not a frozen material binding",
            ));
        }
        let [rho_lo, rho_hi] = domain.density_ratio;
        if bank.curve().points().len() > MAX_FAN_POINTS
            || !source_present(bank.curve().source())
            || !source_present(&energy.source)
            || !positive(reference_density.value())
            || !positive(outlet_area_per_fan.value())
            || !positive(energy.total_efficiency)
            || energy.total_efficiency > 1.0
            || !(0.0..=1.0).contains(&energy.rejected_loss_fraction)
            || !positive(rho_lo)
            || !rho_hi.is_finite()
            || rho_hi < rho_lo
            || !positive(domain.maximum_pressure_fraction)
            || domain.maximum_pressure_fraction > 0.1
        {
            return Err(FanError::Invalid(
                "fan source, geometry, efficiency or applicability",
            ));
        }
        admit_loss(&losses, 0, &mut 0, &mut cancelled)?;
        let resistance = losses.equivalent_resistance().value();
        let (series, parallel) = match bank.arrangement() {
            FanArrangement::Series => (bank.count() as f64, 1.0),
            FanArrangement::Parallel => (1.0, bank.count() as f64),
        };
        let outlet = finite(parallel * outlet_area_per_fan.value())?;
        let velocity_resistance =
            finite(0.5 * reference_density.value() * series / outlet / outlet)?;
        let flow_factor = bank.flow_factor();
        let points = bank.curve().points();
        let flow_domain = [
            finite(bank.curve().admissible_min_flow().value() * flow_factor)?,
            finite(points[points.len() - 1].flow.value() * flow_factor)?,
        ];
        if !positive(resistance)
            || !positive(velocity_resistance)
            || !positive(flow_factor)
            || !positive(flow_domain[1])
            || flow_domain[1] <= flow_domain[0]
        {
            return Err(FanError::Invalid("fan scaling or loss is unrepresentable"));
        }
        let effective_resistance = finite(resistance - velocity_resistance)?;
        for pair in points.windows(2) {
            check(&mut cancelled)?;
            let lo = (pair[0].flow.value() * flow_factor).max(flow_domain[0]);
            let hi = pair[1].flow.value() * flow_factor;
            if hi <= lo {
                continue;
            }
            let slope = finite(
                (pair[1].pressure.value() - pair[0].pressure.value()) * bank.pressure_factor()
                    / ((pair[1].flow.value() - pair[0].flow.value()) * flow_factor),
            )?;
            let derivative_lo = finite(slope - 2.0 * effective_resistance * lo)?;
            let derivative_hi = finite(slope - 2.0 * effective_resistance * hi)?;
            if derivative_lo.max(derivative_hi) > 0.0
                || (slope == 0.0 && effective_resistance == 0.0)
            {
                return Err(FanError::Invalid(
                    "fan-system curve is not strictly decreasing",
                ));
            }
        }
        check(&mut cancelled)?;
        Ok(Self {
            bank,
            losses,
            reference_density,
            outlet_area_per_fan,
            energy,
            domain,
            resistance,
            velocity_resistance,
            flow_domain,
        })
    }
    /// Original speed-scaled bank and its complete source/uncertainty declaration.
    pub const fn bank(&self) -> &FanBank {
        &self.bank
    }
    /// Original total-pressure loss topology and its sources/validity cards.
    pub const fn losses(&self) -> &LossNetwork {
        &self.losses
    }
    /// Density at which both original fan and loss data were declared.
    pub const fn reference_density(&self) -> Density {
        self.reference_density
    }
    /// Original explicit efficiency and heat-rejection assumptions.
    pub const fn energy(&self) -> &FanEnergy {
        &self.energy
    }
    /// Caller-declared applicability of the density-affinity approximation.
    pub const fn domain(&self) -> FanDomain {
        self.domain
    }

    fn pressures(&self, q: f64, ratio: f64) -> Result<(f64, f64), FanError> {
        let velocity = finite(self.velocity_resistance * q * q)?;
        let rise = finite(
            ratio * (self.bank.pressure_at(VolumetricFlowRate::new(q))?.value() + velocity),
        )?;
        let loss = finite(ratio * self.resistance * q * q)?;
        Ok((rise, loss))
    }

    /// Solve against actual absolute endpoint pressures and donor density.
    /// No operation below the original stall limit or beyond its last sample
    /// is inferred. The complete forward fan/duct residual is bracketed, and
    /// both the actual residual and final bracket width must meet their gates.
    pub fn operating_point(
        &self,
        donor_pressure: Pressure,
        receiver_pressure: Pressure,
        donor_density: Density,
        controls: FanControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FanOperatingPoint, FanError> {
        check(&mut cancelled)?;
        if !positive(donor_pressure.value())
            || !positive(receiver_pressure.value())
            || !positive(donor_density.value())
            || !(1..=128).contains(&controls.maximum_iterations)
            || !positive(controls.pressure_tolerance_pa)
            || !positive(controls.flow_tolerance_m3_s)
        {
            return Err(FanError::Invalid("fan endpoint or numerical controls"));
        }
        let ratio = density_affinity_ratio(self.reference_density.value(),
            donor_density.value(), self.domain.density_ratio)
            .map_err(|_| FanError::Invalid("fan density-affinity domain"))?;
        let head = finite(receiver_pressure.value() - donor_pressure.value())?;
        let min_pressure = donor_pressure.value().min(receiver_pressure.value());
        if head.abs() / min_pressure > self.domain.maximum_pressure_fraction {
            return Err(FanError::Invalid("fan endpoint pressure-ratio domain"));
        }
        let mut bounds = self.flow_domain;
        let residual = |q: f64| -> Result<f64, FanError> {
            let (rise, loss) = self.pressures(q, ratio)?;
            finite(finite(rise - loss)? - head)
        };
        let low_residual = residual(bounds[0])?;
        check(&mut cancelled)?;
        let high_residual = residual(bounds[1])?;
        if low_residual < 0.0 || high_residual > 0.0 {
            return Err(FanError::NoOperatingPoint {
                low_residual_pa: low_residual,
                high_residual_pa: high_residual,
            });
        }
        for iterations in 0..=controls.maximum_iterations {
            check(&mut cancelled)?;
            let q = if low_residual == 0.0 {
                bounds[0]
            } else if high_residual == 0.0 {
                bounds[1]
            } else {
                bounds[0] + 0.5 * (bounds[1] - bounds[0])
            };
            let r = residual(q)?;
            let exact_endpoint = low_residual == 0.0 || high_residual == 0.0;
            if r.abs() <= controls.pressure_tolerance_pa
                && (exact_endpoint || bounds[1] - bounds[0] <= controls.flow_tolerance_m3_s)
            {
                if q <= 0.0 {
                    return Err(FanError::Invalid(
                        "zero-flow fan requires a shutoff power model",
                    ));
                }
                let (rise, loss) = self.pressures(q, ratio)?;
                if rise.max(loss) / min_pressure > self.domain.maximum_pressure_fraction {
                    return Err(FanError::Invalid("fan internal pressure-ratio domain"));
                }
                let aerodynamic = finite(q * rise)?;
                let shaft = finite(aerodynamic / self.energy.total_efficiency)?;
                let rejected = finite((shaft - aerodynamic) * self.energy.rejected_loss_fraction)?;
                let gas_power = finite(shaft - rejected)?;
                let mass_rate = finite(donor_density.value() * q)?;
                if !positive(mass_rate) || !positive(aerodynamic) {
                    return Err(FanError::Invalid("fan mass or aerodynamic power underflow"));
                }
                check(&mut cancelled)?;
                return Ok(FanOperatingPoint {
                    volume_flow: VolumetricFlowRate::new(q),
                    mass_rate: MassFlowRate::new(mass_rate),
                    total_pressure_rise: Pressure::new(rise),
                    duct_pressure_drop: Pressure::new(loss),
                    shaft_power: Power::new(shaft),
                    gas_power: Power::new(gas_power),
                    rejected_power: Power::new(rejected),
                    pressure_residual_pa: r,
                    iterations,
                });
            }
            if iterations == controls.maximum_iterations {
                break;
            }
            if r > 0.0 {
                bounds[0] = q;
            } else {
                bounds[1] = q;
            }
        }
        Err(FanError::IterationLimit)
    }
}

/// A point port in world vertical coordinates. The present layer at this point
/// is selected from accepted geometry for each attempted interval; no fixed
/// layer label can continue drawing an absent inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanPort {
    /// Room or the declared ambient reservoir.
    pub endpoint: Endpoint,
    /// Physical point height; finite room ports must lie within that room.
    pub height: Length,
}

/// One forward fan connection, with explicitly equal port elevations.
#[derive(Clone, Debug, PartialEq)]
pub struct FanVent {
    /// Donor point; actual donor composition and enthalpy are transported.
    pub from: FanPort,
    /// Receiving point; fan work enters this gas owner's energy account.
    pub to: FanPort,
    /// Source-bound fan bank, total-pressure losses and energy assumptions.
    pub model: FanModel,
    /// Explicit operating-point work and numerical gates.
    pub controls: FanControl,
    /// True means the fan is stopped and an ideal isolating damper is closed.
    pub closed: bool,
}

/// Operating point used by the most recently accepted network interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanTransfer {
    /// Actual donor layer selected at the point port; arbitrary for ambient.
    pub from_layer: Layer,
    /// Actual receiver layer selected at the point port; arbitrary for ambient.
    pub to_layer: Layer,
    /// None for an explicitly stopped and isolated connection.
    pub operating_point: Option<FanOperatingPoint>,
}

fn port_layer(rooms: &[Room], port: FanPort) -> Layer {
    match port.endpoint {
        Endpoint::Ambient => Layer::Lower,
        Endpoint::Room(i) => {
            let room = &rooms[i];
            let gas = &room.gas;
            if gas.layer(Layer::Lower).mass_kg() == 0.0 {
                Layer::Upper
            } else if gas.layer(Layer::Upper).mass_kg() == 0.0 {
                Layer::Lower
            } else if port.height.value()
                < room.floor_elevation.value() + gas.interface_height().value()
            {
                Layer::Lower
            } else {
                Layer::Upper
            }
        }
    }
}

impl ZoneNetwork {
    /// Atomically install pressure-dependent fan connections between accepted
    /// steps. Sources and point-port geometry are admitted before mutation.
    /// Previous operating-point diagnostics are cleared on topology change.
    pub fn set_fan_vents(
        &mut self,
        fans: Vec<FanVent>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), NetworkError> {
        poll(&mut cancelled)?;
        if self
            .walls
            .len()
            .saturating_add(self.ceilings.len())
            .saturating_add(self.mechanical.len())
            .saturating_add(fans.len())
            > super::MAX_CONNECTIONS
        {
            return Err(NetworkError::Invalid("fan connection budget"));
        }
        for fan in &fans {
            poll(&mut cancelled)?;
            super::admit_endpoints(self.rooms.len(), fan.from.endpoint, fan.to.endpoint)?;
            if !fan.from.height.value().is_finite() || fan.from.height != fan.to.height {
                return Err(NetworkError::Invalid(
                    "fan ports require one finite common elevation",
                ));
            }
            if !(1..=128).contains(&fan.controls.maximum_iterations)
                || !positive(fan.controls.pressure_tolerance_pa)
                || !positive(fan.controls.flow_tolerance_m3_s)
            {
                return Err(NetworkError::Invalid("fan operating controls"));
            }
            for port in [fan.from, fan.to] {
                if let Endpoint::Room(i) = port.endpoint {
                    let room = &self.rooms[i];
                    if port.height.value() < room.floor_elevation.value()
                        || port.height.value()
                            > room.floor_elevation.value() + room.gas.height().value()
                    {
                        return Err(NetworkError::Invalid("fan port outside original room"));
                    }
                }
            }
        }
        let flows = vec![None; fans.len()];
        let trials = flows.clone();
        poll(&mut cancelled)?;
        self.fans = fans;
        self.fan_flows = flows;
        self.fan_trials = trials;
        Ok(())
    }
    /// Original authored fan models, in connection order.
    pub fn fan_vents(&self) -> &[FanVent] {
        &self.fans
    }
    /// Last accepted interval's actual transport operating points. Entries are
    /// absent before the first accepted step or after topology replacement.
    pub fn last_fan_transfers(&self) -> &[Option<FanTransfer>] {
        &self.fan_flows
    }

    pub(super) fn evaluate_fans(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), NetworkError> {
        for (index, fan) in self.fans.iter().enumerate() {
            poll(cancelled)?;
            let from_layer = port_layer(&self.rooms, fan.from);
            let to_layer = port_layer(&self.rooms, fan.to);
            let operating_point = if fan.closed {
                None
            } else {
                let from = super::column(
                    &self.rooms,
                    &self.ambient,
                    self.ambient_density,
                    &self.perturbation,
                    fan.from.endpoint,
                )?;
                let to = super::column(
                    &self.rooms,
                    &self.ambient,
                    self.ambient_density,
                    &self.perturbation,
                    fan.to.endpoint,
                )?;
                let point = fan.model.operating_point(
                    from.pressure_at(fan.from.height, self.gravity)?,
                    to.pressure_at(fan.to.height, self.gravity)?,
                    from.density(vent_layer(from_layer)),
                    fan.controls,
                    &mut *cancelled,
                )?;
                let carried = super::carried(
                    &self.rooms,
                    &self.ambient,
                    fan.from.endpoint,
                    from_layer,
                    point.mass_rate,
                )?;
                super::transfer(
                    &mut self.rates,
                    &mut self.external,
                    fan.from.endpoint,
                    from_layer,
                    fan.to.endpoint,
                    to_layer,
                    carried,
                )?;
                if let Endpoint::Room(i) = fan.to.endpoint {
                    self.rates[i][to_layer.index()].add_heat(point.gas_power)?;
                    self.external.add_heat(point.gas_power)?;
                }
                // Total enthalpy already carries pressure-flow work. Duct
                // dissipation is an internal conversion, not another heat source.
                // Work deposited downstream into ambient does not enter a room.
                Some(point)
            };
            self.fan_trials[index] = Some(FanTransfer {
                from_layer,
                to_layer,
                operating_point,
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
