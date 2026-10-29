//! Through-thickness wall storage on the original `fs-thermal` enthalpy owner.
//!
//! Uniform finite volumes share the original series-conductance law. Surface
//! half-cell resistance is retained when balancing area-weighted gas contacts,
//! an outer reservoir, and explicitly supplied incident radiant power. The
//! latter is EXTERNAL wall input: it is never silently debited from gas energy.
//! Signed internal radiation from a coupled enclosure has a separate input and
//! report account; the caller applies the opposite gas/surface exchanges once.
//! No gas emissivity, soot absorption, material heat capacity, or radiation
//! temperature is inferred. A single profile averages the declared face area;
//! lateral variations require separate wall patches owned by the caller.

use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_qty::{Area, Density, Energy, Length, Power, Qty, Temperature, Time};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, StepCoefficients,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, EnergyAudit, Event, MaterialId, Reservoir, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent, series_conductance_w_k,
};
use std::collections::BTreeMap;

/// Bounded cell axis for the shared dense Newton solve.
pub const MAX_WALL_CELLS: usize = 128;
/// Thermal conductivity [W/(m K)].
pub type ThermalConductivity = Qty<1, 1, -3, -1, 0, 0>;
/// Convective heat-transfer coefficient [W/(m² K)].
pub type HeatTransferCoefficient = Qty<0, 1, -3, -1, 0, 0>;
/// Original material-specific enthalpy [J/kg].
pub type SpecificEnthalpy = Qty<2, 0, -2, 0, 0, 0>;

/// Frozen homogeneous wall geometry and explicitly supplied material properties.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallGeometry {
    /// Area of each through-thickness face.
    pub area: Area,
    /// Total wall thickness.
    pub thickness: Length,
    /// Stored reference density; no thermal expansion is inferred.
    pub density: Density,
    /// Positive conductivity retained over this wall's declared source domain.
    pub conductivity: ThermalConductivity,
    /// Number of equal finite volumes through the thickness.
    pub cells: usize,
}
/// One layer's contact with the common inner surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasContact {
    /// Nonnegative contacted area. The two areas may not exceed the wall area.
    pub area: Area,
    /// Prescribed gas endpoint temperature, held during this candidate solve.
    pub temperature: Temperature,
    /// Nonnegative explicitly resolved convective coefficient; zero insulates.
    pub htc: HeatTransferCoefficient,
}
/// Conditions frozen over one implicit candidate interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallBoundary {
    /// Contacts ordered `[lower gas layer, upper gas layer]`.
    pub gas: [GasContact; 2],
    /// Prescribed temperature on the outer side.
    pub ambient_temperature: Temperature,
    /// Outer Robin coefficient over the complete face; zero insulates.
    pub ambient_htc: HeatTransferCoefficient,
    /// Nonnegative prescribed radiant input at the INNER surface [W]. This is
    /// an external source, e.g. an explicitly allocated fire radiant fraction.
    pub external_radiant_power: Power,
}
/// Explicit nonlinear budget and original local/aggregate energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallSolve {
    /// Maximum Newton iterations, 1 through 64.
    pub maximum_iterations: usize,
    /// Maximum backtracking probes per Newton update, 1 through 32.
    pub maximum_backtracks: usize,
    /// Minimum finite-difference energy increment [J], adjusted upward only
    /// when the original stored enthalpy coordinate needs a larger increment.
    pub energy_difference_j: f64,
    /// Maximum absolute original per-cell energy residual [J].
    pub cell_energy_tolerance_j: f64,
    /// Maximum absolute original whole-wall energy residual [J].
    pub total_energy_tolerance_j: f64,
}
impl Default for WallSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 24,
            maximum_backtracks: 24,
            energy_difference_j: 1e-3,
            cell_energy_tolerance_j: 1e-6,
            total_energy_tolerance_j: 1e-5,
        }
    }
}
/// Candidate diagnostics from exactly the accepted thermal endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallReport {
    /// Accepted duration; no automatic substep reduction is performed here.
    pub elapsed: Time,
    /// Newton evaluations of the unperturbed candidate.
    pub iterations: usize,
    /// Original largest per-cell residual before the original acceptance gate.
    pub maximum_cell_residual_j: f64,
    /// Gas-layer powers [W], positive INTO gas, opposite to their wall contacts.
    pub gas_heat_power: [Power; 2],
    /// Gas-layer energies [J], the identical powers integrated over the step.
    pub gas_heat_energy: [Energy; 2],
    /// Signed energy [J] entering the wall from the outer ambient.
    pub ambient_heat_energy: Energy,
    /// Prescribed external radiant energy entering the wall, without a gas debit.
    pub external_radiant_energy: Energy,
    /// Signed coupled radiant energy entering the inner face. This is an
    /// internal exchange with the explicitly coupled enclosure, not fire input.
    pub coupled_radiant_energy: Energy,
    /// Inner and outer FACE temperatures, including half-cell resistance.
    pub surface_temperatures: [Temperature; 2],
    /// Original storage/external-heat audit, not a newly reconstructed account.
    pub audit: EnergyAudit,
}
/// Physical candidate and its coupled exchange; the original wall is unchanged.
pub struct WallCandidate<'a> {
    /// The original thermal owner advanced to the candidate endpoint.
    pub wall: Wall1d<'a>,
    /// Exact endpoint exchanges and original acceptance observations.
    pub report: WallReport,
}
/// A wall retains the original borrowed caloric law and original thermal history.
#[derive(Clone)]
pub struct Wall1d<'a> {
    state: DepositionSimulation<'a>,
    geometry: WallGeometry,
    half_cell_conductance: f64,
}
/// Reusable numerical storage; it carries no physical wall state.
pub struct WallWorkspace {
    energy: Vec<f64>,
    probe: Vec<f64>,
    residual: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    delta: Vec<f64>,
    lu: LuWorkspace,
}
/// Atomic source, numerical, work-limit or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum WallError {
    /// Invalid geometry, coefficient, boundary, or solver control.
    Invalid(&'static str),
    /// Original material/storage/energy admission refused.
    Thermal(DepositionError),
    /// Original prepared dense solve refused the Jacobian.
    Linear(LuWorkspaceError),
    /// Bounded Newton/backtracking failed to close the actual energy residual.
    Iterations,
    /// No candidate was published.
    Cancelled,
}
impl core::fmt::Display for WallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compartment wall: {self:?}")
    }
}
impl std::error::Error for WallError {}
impl From<DepositionError> for WallError {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl From<LuWorkspaceError> for WallError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Linear(e)
    }
}
fn finite(x: f64) -> Result<f64, WallError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

impl<'a> Wall1d<'a> {
    /// Construct the actual cell masses, thermal contacts and original enthalpy
    /// timeline. Initial enthalpy is explicit because `ThermalLaw` owns T(h);
    /// the adapter never substitutes a Cp or guesses the inverse chart.
    pub fn new(
        law: &'a dyn ThermalLaw,
        geometry: WallGeometry,
        initial_specific_enthalpy: SpecificEnthalpy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallError> {
        checkpoint(&mut cancelled)?;
        let n = geometry.cells;
        if !(1..=MAX_WALL_CELLS).contains(&n)
            || ![
                geometry.area.value(),
                geometry.thickness.value(),
                geometry.density.value(),
                geometry.conductivity.value(),
            ]
            .into_iter()
            .all(positive)
            || !initial_specific_enthalpy.value().is_finite()
        {
            return Err(WallError::Invalid(
                "positive finite wall properties, bounded cells and finite initial enthalpy required",
            ));
        }
        let dx = geometry.thickness.value() / n as f64;
        let volume = finite(geometry.area.value() * dx)?;
        let mass = finite(geometry.density.value() * volume)?;
        let k = geometry.conductivity.value();
        let conductance = series_conductance_w_k(geometry.area.value(), [0.5 * dx; 2], [k; 2])?;
        let half_cell_conductance =
            series_conductance_w_k(geometry.area.value(), [0.25 * dx; 2], [k; 2])?;
        let initial_temperature = law.temperature_k(initial_specific_enthalpy.value())?;
        let mut schedule = Vec::with_capacity(2 * n - 1);
        for i in 0..n {
            checkpoint(&mut cancelled)?;
            schedule.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64 + 1),
                    segment_id: 1,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: mass,
                    specific_h_j_kg: initial_specific_enthalpy.value(),
                    geometry: CellGeometry {
                        start_m: [i as f64 * dx, 0.0, 0.0],
                        end_m: [(i + 1) as f64 * dx, 0.0, 0.0],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * geometry.area.value(),
                    },
                    // All changing reservoir contacts enter through the original
                    // endpoint coefficient callback; nothing persists twice.
                    boundary: Boundary {
                        ambient_k: initial_temperature,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
        for i in 0..n - 1 {
            schedule.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64 + 1,
                    a: CellId(i as u64 + 1),
                    b: CellId(i as u64 + 2),
                    area_m2: geometry.area.value(),
                    conductance_w_k: conductance,
                }),
            });
        }
        let state = DepositionSimulation::new(
            schedule,
            &[(MaterialId(1), law)],
            DepositionBudget {
                max_cells: n,
                max_links: n - 1,
                max_events: 2 * n - 1,
                max_substeps: 1,
            },
        )?
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || checkpoint(&mut cancelled),
        )?
        .0;
        checkpoint(&mut cancelled)?;
        Ok(Self {
            state,
            geometry,
            half_cell_conductance,
        })
    }
    /// Immutable access to the original cell/contact/material history.
    pub const fn simulation(&self) -> &DepositionSimulation<'a> {
        &self.state
    }
    /// Original uniform wall geometry and explicit frozen transport declarations.
    pub const fn geometry(&self) -> WallGeometry {
        self.geometry
    }
    /// Accepted physical clock.
    pub fn time(&self) -> Time {
        Time::new(self.state.time_s())
    }
    /// Temperature of a zero-based through-thickness cell.
    pub fn cell_temperature(&self, index: usize) -> Result<Temperature, WallError> {
        if index >= self.geometry.cells {
            return Err(WallError::Invalid("wall cell index"));
        }
        Ok(Temperature::new(
            self.state.temperature_k(CellId(index as u64 + 1))?,
        ))
    }
    /// Prepare bounded numerical work once for this cell count.
    pub fn workspace(&self) -> Result<WallWorkspace, WallError> {
        let n = self.geometry.cells;
        Ok(WallWorkspace {
            energy: vec![0.0; n],
            probe: vec![0.0; n],
            residual: vec![0.0; n],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            delta: vec![0.0; n],
            lu: LuWorkspace::new(n)?,
        })
    }
    /// Solve one backward-Euler candidate through the original thermal residual.
    /// Failed Newton trials change only numerical workspace. Native thermal
    /// candidates retain their original map allocations and material-domain
    /// checks; this adapter does not claim an allocation-free thermal solve.
    pub fn candidate(
        &self,
        duration: Time,
        boundary: WallBoundary,
        controls: WallSolve,
        work: &mut WallWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<WallCandidate<'a>, WallError> {
        self.candidate_with_exchange(
            duration,
            boundary,
            Power::new(0.0),
            controls,
            work,
            cancelled,
        )
    }

    /// Solve with signed inner-face radiation supplied by a coupled enclosure.
    /// Positive power enters this wall; negative power leaves it. The same
    /// power must be paired with the enclosure's other surface/gas exchanges.
    /// It is held over this candidate interval and recorded separately from
    /// prescribed external fire radiation. Conduction and convective contact
    /// use the actual backward-Euler wall endpoint and half-cell resistance.
    ///
    /// This is a frozen-flux coupling seam, not an unconditionally stable
    /// implicit radiation integrator. An outer coupled iteration or timestep
    /// refinement must control temperature-dependent radiation error.
    pub fn candidate_with_exchange(
        &self,
        duration: Time,
        boundary: WallBoundary,
        coupled_radiant_power: Power,
        controls: WallSolve,
        work: &mut WallWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WallCandidate<'a>, WallError> {
        checkpoint(&mut cancelled)?;
        self.validate_boundary(boundary)?;
        if !coupled_radiant_power.value().is_finite() {
            return Err(WallError::Invalid(
                "finite signed coupled radiation required",
            ));
        }
        let n = self.geometry.cells;
        if !positive(duration.value())
            || work.energy.len() != n
            || !(1..=64).contains(&controls.maximum_iterations)
            || !(1..=32).contains(&controls.maximum_backtracks)
            || !positive(controls.energy_difference_j)
            || !nonnegative(controls.cell_energy_tolerance_j)
            || !nonnegative(controls.total_energy_tolerance_j)
        {
            return Err(WallError::Invalid(
                "wall duration, workspace or nonlinear controls",
            ));
        }
        let step =
            BackwardEulerStep::new(&self.state, self.state.time_s() + duration.value(), || {
                checkpoint(&mut cancelled)
            })?;
        work.energy.fill(0.0);
        for iteration in 0..controls.maximum_iterations {
            checkpoint(&mut cancelled)?;
            let trial = self.evaluate(
                &step,
                &work.energy,
                boundary,
                coupled_radiant_power,
                &mut cancelled,
            )?;
            let norm = trial.maximum_cell_residual_j();
            if norm <= controls.cell_energy_tolerance_j
                && trial.audit().residual_j.abs() <= controls.total_energy_tolerance_j
            {
                let (state, exchange, audit) = trial.accept(
                    controls.cell_energy_tolerance_j,
                    controls.total_energy_tolerance_j,
                    || checkpoint(&mut cancelled),
                )?;
                let dt = step.duration_s();
                let report = WallReport {
                    elapsed: Time::new(dt),
                    iterations: iteration + 1,
                    maximum_cell_residual_j: norm,
                    gas_heat_power: exchange.gas_into_wall.map(|q| Power::new(-q)),
                    gas_heat_energy: exchange.gas_into_wall.map(|q| Energy::new(-q * dt)),
                    ambient_heat_energy: Energy::new(exchange.ambient_into_wall * dt),
                    external_radiant_energy: Energy::new(
                        boundary.external_radiant_power.value() * dt,
                    ),
                    coupled_radiant_energy: Energy::new(coupled_radiant_power.value() * dt),
                    surface_temperatures: exchange.surfaces.map(Temperature::new),
                    audit,
                };
                if !report
                    .gas_heat_energy
                    .iter()
                    .map(|q| q.value())
                    .chain([
                        report.ambient_heat_energy.value(),
                        report.external_radiant_energy.value(),
                        report.coupled_radiant_energy.value(),
                    ])
                    .all(f64::is_finite)
                {
                    return Err(DepositionError::Arithmetic.into());
                }
                checkpoint(&mut cancelled)?;
                return Ok(WallCandidate {
                    wall: Self {
                        state,
                        geometry: self.geometry,
                        half_cell_conductance: self.half_cell_conductance,
                    },
                    report,
                });
            }
            work.residual.copy_from_slice(trial.residual_j());
            drop(trial);
            for j in 0..n {
                checkpoint(&mut cancelled)?;
                let cell = &self.state.cells()[&CellId(j as u64 + 1)];
                let stored_scale = finite(cell.mass_kg * cell.specific_h_j_kg.abs())?;
                let delta = controls
                    .energy_difference_j
                    .max(64.0 * f64::EPSILON * stored_scale.max(work.energy[j].abs()));
                let mut derivative = None;
                for sign in [1.0, -1.0] {
                    work.probe.copy_from_slice(&work.energy);
                    work.probe[j] = finite(work.energy[j] + sign * delta)?;
                    let actual_delta = work.probe[j] - work.energy[j];
                    if actual_delta == 0.0 {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    match self.evaluate(
                        &step,
                        &work.probe,
                        boundary,
                        coupled_radiant_power,
                        &mut cancelled,
                    ) {
                        Ok(probe) => {
                            derivative = Some((probe, actual_delta));
                            break;
                        }
                        Err(WallError::Cancelled) => return Err(WallError::Cancelled),
                        Err(WallError::Thermal(
                            DepositionError::MaterialDomain | DepositionError::Arithmetic,
                        )) => {}
                        Err(e) => return Err(e),
                    }
                }
                let (probe, delta) = derivative.ok_or(WallError::Iterations)?;
                for i in 0..n {
                    work.matrix[i * n + j] =
                        finite((probe.residual_j()[i] - work.residual[i]) / delta)?;
                }
            }
            for i in 0..n {
                work.rhs[i] = -work.residual[i];
            }
            checkpoint(&mut cancelled)?;
            work.lu
                .solve_into(&work.matrix, &work.rhs, &mut work.delta)?;
            let mut scale = 1.0;
            let mut accepted = false;
            for _ in 0..controls.maximum_backtracks {
                checkpoint(&mut cancelled)?;
                for i in 0..n {
                    work.probe[i] = finite(work.energy[i] + scale * work.delta[i])?;
                }
                match self.evaluate(
                    &step,
                    &work.probe,
                    boundary,
                    coupled_radiant_power,
                    &mut cancelled,
                ) {
                    Ok(probe) if probe.maximum_cell_residual_j() < norm => {
                        accepted = true;
                        break;
                    }
                    Err(WallError::Cancelled) => return Err(WallError::Cancelled),
                    _ => {}
                }
                scale *= 0.5;
            }
            if !accepted {
                return Err(WallError::Iterations);
            }
            work.energy.copy_from_slice(&work.probe);
        }
        Err(WallError::Iterations)
    }
    fn validate_boundary(&self, boundary: WallBoundary) -> Result<(), WallError> {
        if !positive(boundary.ambient_temperature.value())
            || !nonnegative(boundary.ambient_htc.value())
            || !nonnegative(boundary.external_radiant_power.value())
            || boundary.gas.iter().any(|c| {
                !positive(c.temperature.value())
                    || !nonnegative(c.htc.value())
                    || !nonnegative(c.area.value())
            })
        {
            return Err(WallError::Invalid(
                "wall contact temperature, area, coefficient or radiant input",
            ));
        }
        let area = finite(boundary.gas[0].area.value() + boundary.gas[1].area.value())?;
        if area > self.geometry.area.value() * (1.0 + 8.0 * f64::EPSILON) {
            return Err(WallError::Invalid("gas contact areas exceed the wall face"));
        }
        Ok(())
    }
    fn evaluate(
        &self,
        step: &BackwardEulerStep<'_, 'a>,
        energy: &[f64],
        boundary: WallBoundary,
        coupled_radiant_power: Power,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BackwardEulerTrial<'a, Exchange>, WallError> {
        step.evaluate(
            energy,
            |state| self.coefficients(state, boundary, coupled_radiant_power),
            || checkpoint(cancelled),
        )
    }
    fn coefficients(
        &self,
        state: &DepositionSimulation<'a>,
        boundary: WallBoundary,
        coupled_radiant_power: Power,
    ) -> Result<StepCoefficients<Exchange>, WallError> {
        let first = CellId(1);
        let last = CellId(self.geometry.cells as u64);
        let inner = state.temperature_k(first)?;
        let outer = state.temperature_k(last)?;
        let mut gas = [Reservoir {
            area_m2: 0.0,
            conductance_w_k: 0.0,
            temperature_k: 300.0,
        }; 2];
        for (i, c) in boundary.gas.iter().enumerate() {
            gas[i] = Reservoir {
                area_m2: c.area.value(),
                conductance_w_k: finite(c.area.value() * c.htc.value())?,
                temperature_k: c.temperature.value(),
            };
        }
        let ambient = [Reservoir {
            area_m2: self.geometry.area.value(),
            conductance_w_k: finite(self.geometry.area.value() * boundary.ambient_htc.value())?,
            temperature_k: boundary.ambient_temperature.value(),
        }];
        let inner_radiant_power =
            finite(boundary.external_radiant_power.value() + coupled_radiant_power.value())?;
        let (inner_surface, gas_into_wall) =
            surface(inner, self.half_cell_conductance, gas, inner_radiant_power)?;
        let (outer_surface, ambient_into_wall) =
            surface(outer, self.half_cell_conductance, ambient, 0.0)?;
        let mut sources_w = BTreeMap::new();
        sources_w.insert(
            first,
            finite(finite(gas_into_wall[0] + gas_into_wall[1])? + inner_radiant_power)?,
        );
        let source = sources_w.entry(last).or_insert(0.0);
        *source = finite(*source + ambient_into_wall[0])?;
        Ok(StepCoefficients {
            sources_w,
            conductances_w_k: BTreeMap::new(),
            payload: Exchange {
                gas_into_wall,
                ambient_into_wall: ambient_into_wall[0],
                surfaces: [inner_surface, outer_surface],
            },
        })
    }
}
#[derive(Clone, Copy)]
struct Exchange {
    gas_into_wall: [f64; 2],
    ambient_into_wall: f64,
    surfaces: [f64; 2],
}

fn surface<const N: usize>(
    cell_temperature: f64,
    bulk: f64,
    contacts: [Reservoir; N],
    external_power: f64,
) -> Result<(f64, [f64; N]), WallError> {
    let mut conductance = bulk;
    let mut drive = external_power;
    for contact in contacts {
        conductance = finite(conductance + contact.conductance_w_k)?;
        drive = finite(
            drive + finite(contact.conductance_w_k * (contact.temperature_k - cell_temperature))?,
        )?;
    }
    let surface_change = finite(drive / conductance)?;
    let temperature = finite(cell_temperature + surface_change)?;
    if !positive(temperature) {
        return Err(WallError::Invalid("nonpositive wall surface temperature"));
    }
    let mut powers = [0.0; N];
    for (i, contact) in contacts.iter().enumerate() {
        // Use the small relative difference even when its absolute surface
        // temperature rounds onto the adjacent cell temperature.
        powers[i] = finite(
            contact.conductance_w_k * ((contact.temperature_k - cell_temperature) - surface_change),
        )?;
    }
    Ok((temperature, powers))
}

#[cfg(test)]
mod tests;
