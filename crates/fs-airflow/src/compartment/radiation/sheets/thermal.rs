//! Implicit heating of existing finite thermal cells by transmitting sheets.
//!
//! The two optical faces of each sheet MUST refer to the same stored CellId.
//! Their true absorption-minus-emission powers enter that cell once, without
//! another mass, heat-capacity law or temperature update. Every original cell
//! and conductive link remains in the original backward-Euler energy residual,
//! including cells without an optical face. Prescribed opaque faces are baths.
//!
//! Gas temperatures, geometry and supplied optical bands stay fixed during the
//! interval. The returned gas/bath powers must be paired with their own owners;
//! this is not a fully implicit gas update. Uniform-cell temperatures require
//! the caller's isothermal/Biot justification. Spectral calculations retain
//! the original finite window, source domains and omitted-emission diagnostics.

use std::cell::RefCell;
use std::collections::BTreeMap;

use fs_conduction::ConductionError;
use fs_conduction::radiation::SurfaceEmissivity;
use fs_gas::multicomponent::compartment::LayerRates;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_math::dd::Dd;
use fs_qty::{Power, Temperature, Time};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit,
    backward_euler::{
        BackwardEulerStep, BackwardEulerTrial, StepCoefficients, resolvable_increment_j,
    },
};

use super::super::{
    LayerAbsorption, ParticipatingEnclosure, RadiationError, RadiationSolve, RadiationWorkspace,
    spectral::{SpectralOptics, SpectralSolve, SpectralWorkspace},
    walls::WallRadiation,
};
use super::SheetPair;

/// Deliberate small-system ceiling for the dense coupled thermal Jacobian.
pub const MAX_SHEET_THERMAL_CELLS: usize = 256;

/// Temperature and thermal ownership of one original optical surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SheetThermalSurface {
    /// Absorbed-minus-emitted power heats this original active thermal cell.
    Stored(CellId),
    /// An opaque surface held at this temperature by an external bath.
    Prescribed(Temperature),
}

/// Source-bearing optics held over one thermal candidate interval.
#[derive(Clone, Copy)]
pub enum SheetThermalOptics<'a, 's> {
    /// Total-gray sheets and original gas absorption.
    Gray {
        absorption: &'a [LayerAbsorption; 2],
        sheets: &'a [SheetPair],
        solve: RadiationSolve,
    },
    /// Original finite-window solver, including its declared sheet bands.
    Spectral {
        optics: &'a SpectralOptics<'s>,
        solve: SpectralSolve,
    },
}

impl SheetThermalOptics<'_, '_> {
    fn sheets(&self) -> &[SheetPair] {
        match self {
            Self::Gray { sheets, .. } => sheets,
            Self::Spectral { optics, .. } => optics
                .sheet_bands()
                .and_then(|bands| bands.first())
                .map_or(&[], Vec::as_slice),
        }
    }
}

/// Work ceilings and independent physical energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SheetThermalSolve {
    /// Outer Newton residual evaluations, 1 through 64.
    pub maximum_iterations: usize,
    /// Globalization probes per update AND probe halvings per Jacobian column,
    /// 1 through 32. Every probe is also charged to maximum_evaluations.
    pub maximum_backtracks: usize,
    /// ALL thermal trials, including Jacobian probes and rejected steps.
    /// Each optical trial also retains its original gray/spectral work limits.
    pub maximum_evaluations: usize,
    /// Nominal finite-difference increment in the original cell energy [J].
    /// Raised for coordinate resolution, halved at material/optical support edges.
    pub energy_difference_j: f64,
    /// Original maximum absolute individual-cell residual [J].
    pub cell_energy_tolerance_j: f64,
    /// Original complete thermal inventory residual [J].
    pub total_energy_tolerance_j: f64,
    /// Thermal residual plus the actual gas/bath optical closure defect [J].
    pub combined_energy_tolerance_j: f64,
}

impl Default for SheetThermalSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 32,
            maximum_backtracks: 24,
            maximum_evaluations: 16_384,
            energy_difference_j: 1e-3,
            cell_energy_tolerance_j: 1e-6,
            total_energy_tolerance_j: 1e-5,
            combined_energy_tolerance_j: 1e-4,
        }
    }
}

/// Accepted exchanges evaluated at the actual returned thermal endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetThermalReport<'s> {
    pub elapsed: Time,
    pub iterations: usize,
    /// Counts every original thermal trial, including unsuccessful probes.
    pub evaluations: usize,
    /// Original gray/spectral report, with sheet sources and spectral coverage.
    pub radiation: WallRadiation<'s>,
    /// Endpoint temperatures in the original optical surface order.
    pub surface_temperatures: Vec<Temperature>,
    /// Powers already applied to stored cells. Multiple faces are summed once.
    /// Do NOT add these powers to the returned thermal state again.
    pub cell_radiant_power: BTreeMap<CellId, Power>,
    /// Watts INTO prescribed baths, in optical order; zero for stored faces.
    /// Transmission between faces is not assigned to the sheet's heat storage.
    pub prescribed_surface_heat_power: Vec<Power>,
    /// Original material-backed storage audit. Its external heat includes the
    /// radiant input above plus all original authored sources and boundaries.
    pub audit: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    /// Signed sum of the ACTUALLY applied cell, gas and prescribed-bath powers.
    /// Includes any rounding introduced when several faces share one cell.
    pub applied_power_residual_w: f64,
    /// Original thermal residual plus dt times the applied power residual.
    /// This is the combined storage/gas/bath balance after excluding authored
    /// thermal external heat. No defect is redistributed into a physical owner.
    pub combined_energy_residual_j: f64,
}

impl SheetThermalReport<'_> {
    /// Add the identical radiative gas powers once. Gas is not advanced here;
    /// install the thermal candidate only with acceptance of its coupled owners.
    pub fn add_gas_heat(&self, rates: &mut [LayerRates; 2]) -> Result<(), RadiationError> {
        match &self.radiation {
            WallRadiation::Gray(report) => report.add_gas_heat(rates),
            WallRadiation::Spectral(report) => report.add_gas_heat(rates),
        }
    }
}

/// Original thermal owner at its accepted endpoint; the input is unchanged.
pub struct SheetThermalCandidate<'m, 's> {
    pub thermal: DepositionSimulation<'m>,
    pub report: SheetThermalReport<'s>,
}

/// Reusable numerical storage; contains no accepted physical state.
pub struct SheetThermalWorkspace {
    thermal: ThermalWorkspace,
    optical: OpticalWorkspace,
}

/// Numerical buffers shared by the sheet and non-isothermal slab adapters.
pub(crate) struct ThermalWorkspace {
    ids: Vec<CellId>,
    energy: Vec<f64>,
    probe: Vec<f64>,
    residual: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    delta: Vec<f64>,
    lu: LuWorkspace,
}

impl ThermalWorkspace {
    pub(crate) fn new(
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Self, SheetThermalError> {
        let n = thermal.cells().len();
        if !(1..=MAX_SHEET_THERMAL_CELLS).contains(&n) {
            return Err(SheetThermalError::Invalid(
                "sheet thermal storage needs 1..=256 cells",
            ));
        }
        Ok(Self {
            ids: thermal.cells().keys().copied().collect(),
            energy: vec![0.0; n],
            probe: vec![0.0; n],
            residual: vec![0.0; n],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            delta: vec![0.0; n],
            lu: LuWorkspace::new(n)?,
        })
    }
}

struct OpticalWorkspace {
    gray: RadiationWorkspace,
    spectral: SpectralWorkspace,
}

/// A refusal never modifies the original thermal inventory.
#[derive(Debug)]
pub enum SheetThermalError {
    Invalid(&'static str),
    Thermal(DepositionError),
    Radiation(RadiationError),
    Emissivity {
        surface: usize,
        error: ConductionError,
    },
    /// Sourced volumetric optics or their original angular preparation refused.
    LayerOptics {
        band: Option<usize>,
        slab: usize,
        layer: Option<usize>,
        error: ConductionError,
    },
    Linear(LuWorkspaceError),
    Budget {
        evaluations: usize,
        maximum_evaluations: usize,
    },
    Iterations,
    EnergyResidual {
        residual_j: f64,
        allowed_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SheetThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transmitting-sheet thermal coupling: {self:?}")
    }
}
impl std::error::Error for SheetThermalError {}
impl From<DepositionError> for SheetThermalError {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl From<RadiationError> for SheetThermalError {
    fn from(error: RadiationError) -> Self {
        if error == RadiationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Radiation(error)
        }
    }
}
impl From<LuWorkspaceError> for SheetThermalError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Linear(error)
    }
}

impl ParticipatingEnclosure {
    /// Prepare a bounded Jacobian for ALL active original thermal cells.
    /// Original thermal trials retain their normal candidate allocations.
    pub fn sheet_thermal_workspace(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<SheetThermalWorkspace, SheetThermalError> {
        Ok(SheetThermalWorkspace {
            thermal: ThermalWorkspace::new(thermal)?,
            optical: OpticalWorkspace {
                gray: self.workspace()?,
                spectral: self.spectral_workspace()?,
            },
        })
    }

    /// Advance original stored cells with endpoint sheet/opaque-surface radiation.
    ///
    /// Both faces of each sheet must map to the SAME Stored cell; their area and
    /// reciprocal optical assumptions remain enforced by the original solver.
    /// Other surfaces may map to cells or prescribed baths. A cell may own several
    /// patches, but its inventory is never duplicated. Original sources, boundary
    /// convection and links stay additive; do not author this radiation twice.
    ///
    /// The gray emissivity callback is queried ONLY for unpaired opaque faces,
    /// at each actual trial temperature. Spectral surface bands skip it entirely.
    /// Fixed optics must cover every trial temperature in their source domains.
    /// All trials use immutable old history. Acceptance uses the original local
    /// and total joule residuals, not a small temperature/Newton update.
    ///
    /// Process pending start-time events through the original thermal owner first.
    /// An event strictly inside this interval refuses; target-time events remain
    /// pending in the accepted left-limit state. Process those events only after
    /// all coupled owners accept. Elapsed time is the original owner's actual
    /// representable interval. Supplied gas/bath states must describe that same
    /// interval; this method has no independent gas clock to check.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_sheet_thermal<'m, 's>(
        &self,
        thermal: &DepositionSimulation<'m>,
        duration: Time,
        surfaces: &[SheetThermalSurface],
        gas_temperatures: [Temperature; 2],
        optics: SheetThermalOptics<'_, 's>,
        mut emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        controls: SheetThermalSolve,
        work: &mut SheetThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SheetThermalCandidate<'m, 's>, SheetThermalError> {
        checkpoint(&mut cancelled)?;
        if surfaces.len() != self.geometry().len()
            || gas_temperatures.iter().any(|t| !positive(t.value()))
        {
            return Err(SheetThermalError::Invalid("sheet thermal optical axes or gas temperature"));
        }
        let paired = admit_surfaces(thermal, surfaces, optics.sheets())?;
        let optical = &mut work.optical;
        candidate_thermal(
            thermal, duration, controls, &mut work.thermal,
            |state, cancel| resolve_sheet(
                self, state, surfaces, &paired, gas_temperatures, optics,
                &mut emissivity, optical, cancel,
            ),
            cancelled,
        )
    }
}

/// One original-history thermal solve; only the endpoint optical resolver varies.
pub(crate) fn candidate_thermal<'m, 's>(
    thermal: &DepositionSimulation<'m>,
    duration: Time,
    controls: SheetThermalSolve,
    work: &mut ThermalWorkspace,
    mut resolve: impl FnMut(
        &DepositionSimulation<'m>, &mut dyn FnMut() -> bool,
    ) -> Result<Applied<'s>, SheetThermalError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SheetThermalCandidate<'m, 's>, SheetThermalError> {
    checkpoint(&mut cancelled)?;
    let n = thermal.cells().len();
    if !(1..=MAX_SHEET_THERMAL_CELLS).contains(&n)
        || work.ids.len() != n
        || !work.ids.iter().copied().eq(thermal.cells().keys().copied())
        || !positive(duration.value())
        || !(1..=64).contains(&controls.maximum_iterations)
        || !(1..=32).contains(&controls.maximum_backtracks)
        || !(1..=65_536).contains(&controls.maximum_evaluations)
        || !positive(controls.energy_difference_j)
        || !nonnegative(controls.cell_energy_tolerance_j)
        || !nonnegative(controls.total_energy_tolerance_j)
        || !nonnegative(controls.combined_energy_tolerance_j)
    {
        return Err(SheetThermalError::Invalid("sheet thermal axes, workspace or controls"));
    }
    let step = BackwardEulerStep::new(
        thermal,
        finite(thermal.time_s() + duration.value())?,
        || checkpoint(&mut cancelled),
    )?;
    work.energy.fill(0.0);
    let mut evaluations = 0;
    for iteration in 0..controls.maximum_iterations {
        checkpoint(&mut cancelled)?;
        let trial = evaluate(
            &step,
            &work.energy,
            controls,
            &mut resolve,
            &mut evaluations,
            &mut cancelled,
        )?;
        let maximum_cell = trial.maximum_cell_residual_j();
        if maximum_cell <= controls.cell_energy_tolerance_j
            && trial.audit().residual_j.abs() <= controls.total_energy_tolerance_j
        {
            return finish(
                trial,
                step.duration_s(),
                iteration + 1,
                evaluations,
                maximum_cell,
                controls,
                &mut cancelled,
            );
        }
        let norm = maximum_cell.max(trial.audit().residual_j.abs());
        work.residual.copy_from_slice(trial.residual_j());
        drop(trial);
        for column in 0..n {
            checkpoint(&mut cancelled)?;
            let cell = &thermal.cells()[&work.ids[column]];
            let scale = finite(cell.mass_kg * cell.specific_h_j_kg.abs())?;
            let mut increment = controls
                .energy_difference_j
                .max(64.0 * f64::EPSILON * scale.max(work.energy[column].abs()));
            let mut derivative = None;
            for _ in 0..controls.maximum_backtracks {
                let mut representable = false;
                for sign in [1.0, -1.0] {
                    work.probe.copy_from_slice(&work.energy);
                    let value = work.energy[column] + sign * increment;
                    if !value.is_finite() {
                        continue;
                    }
                    work.probe[column] =
                        resolvable_increment_j(cell.specific_h_j_kg, value, cell.mass_kg);
                    let difference = finite(work.probe[column] - work.energy[column])?;
                    if difference == 0.0 {
                        continue;
                    }
                    representable = true;
                    match evaluate(
                        &step,
                        &work.probe,
                        controls,
                        &mut resolve,
                        &mut evaluations,
                        &mut cancelled,
                    ) {
                        Ok(probe) => {
                            derivative = Some((probe, difference));
                            break;
                        }
                        Err(error) if trial_domain_error(&error) => {}
                        Err(error) => return Err(error),
                    }
                }
                if derivative.is_some() || !representable {
                    break;
                }
                increment *= 0.5;
            }
            let (probe, difference) = derivative.ok_or(SheetThermalError::Iterations)?;
            for row in 0..n {
                work.matrix[row * n + column] =
                    finite((probe.residual_j()[row] - work.residual[row]) / difference)?;
            }
        }
        for row in 0..n {
            work.rhs[row] = -work.residual[row];
        }
        checkpoint(&mut cancelled)?;
        work.lu
            .solve_into(&work.matrix, &work.rhs, &mut work.delta)?;
        let mut scale = 1.0;
        let mut improved = false;
        for _ in 0..controls.maximum_backtracks {
            checkpoint(&mut cancelled)?;
            let mut representable = true;
            for (i, id) in work.ids.iter().enumerate() {
                let cell = &thermal.cells()[id];
                let value = work.energy[i] + scale * work.delta[i];
                if !value.is_finite() {
                    representable = false;
                    break;
                }
                work.probe[i] =
                    resolvable_increment_j(cell.specific_h_j_kg, value, cell.mass_kg);
            }
            if !representable {
                scale *= 0.5;
                continue;
            }
            match evaluate(
                &step,
                &work.probe,
                controls,
                &mut resolve,
                &mut evaluations,
                &mut cancelled,
            ) {
                Ok(probe) => {
                    let candidate_norm = probe
                        .maximum_cell_residual_j()
                        .max(probe.audit().residual_j.abs());
                    if candidate_norm < norm {
                        improved = true;
                        break;
                    }
                }
                Err(error) if trial_domain_error(&error) => {}
                Err(error) => return Err(error),
            }
            scale *= 0.5;
        }
        if !improved {
            return Err(SheetThermalError::Iterations);
        }
        work.energy.copy_from_slice(&work.probe);
    }
    Err(SheetThermalError::Iterations)
}

fn admit_surfaces(
    thermal: &DepositionSimulation<'_>,
    surfaces: &[SheetThermalSurface],
    sheets: &[SheetPair],
) -> Result<Vec<bool>, SheetThermalError> {
    for surface in surfaces {
        match *surface {
            SheetThermalSurface::Stored(id) if thermal.cells().contains_key(&id) => {}
            SheetThermalSurface::Prescribed(t) if positive(t.value()) => {}
            _ => {
                return Err(SheetThermalError::Invalid(
                    "inactive cell or invalid bath temperature",
                ));
            }
        }
    }
    let mut paired = vec![false; surfaces.len()];
    for sheet in sheets {
        let [a, b] = sheet.surfaces;
        if a == b || a >= surfaces.len() || b >= surfaces.len() || paired[a] || paired[b] {
            return Err(SheetThermalError::Invalid(
                "invalid or repeated sheet faces",
            ));
        }
        match (surfaces[a], surfaces[b]) {
            (SheetThermalSurface::Stored(first), SheetThermalSurface::Stored(second))
                if first == second => {}
            _ => {
                return Err(SheetThermalError::Invalid(
                    "both sheet faces must share one original stored thermal cell",
                ));
            }
        }
        paired[a] = true;
        paired[b] = true;
    }
    Ok(paired)
}

pub(crate) struct Applied<'s> {
    radiation: WallRadiation<'s>,
    temperatures: Vec<Temperature>,
    cell_power: BTreeMap<CellId, Power>,
    bath_power: Vec<Power>,
    power_residual_w: f64,
}

fn evaluate<'m, 's>(
    step: &BackwardEulerStep<'_, 'm>,
    energy: &[f64],
    controls: SheetThermalSolve,
    resolve: &mut impl FnMut(
        &DepositionSimulation<'m>, &mut dyn FnMut() -> bool,
    ) -> Result<Applied<'s>, SheetThermalError>,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BackwardEulerTrial<'m, Applied<'s>>, SheetThermalError> {
    checkpoint(cancelled)?;
    if *evaluations >= controls.maximum_evaluations {
        return Err(SheetThermalError::Budget {
            evaluations: *evaluations,
            maximum_evaluations: controls.maximum_evaluations,
        });
    }
    *evaluations += 1;
    // The original thermal and optical callbacks execute serially.
    let cancellation = RefCell::new(cancelled);
    step.evaluate(
        energy,
        |state| {
            let applied = resolve(state, &mut || (*cancellation.borrow_mut())())?;
            Ok(StepCoefficients {
                sources_w: applied.cell_power.iter()
                    .map(|(&id, &power)| (id, power.value())).collect(),
                conductances_w_k: BTreeMap::new(),
                payload: applied,
            })
        },
        || checkpoint(&mut || (*cancellation.borrow_mut())()),
    )
}

#[allow(clippy::too_many_arguments)]
fn resolve_sheet<'s>(
    enclosure: &ParticipatingEnclosure,
    state: &DepositionSimulation<'_>,
    surfaces: &[SheetThermalSurface],
    paired: &[bool],
    gas_temperatures: [Temperature; 2],
    optics: SheetThermalOptics<'_, 's>,
    emissivity: &mut impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
    work: &mut OpticalWorkspace,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Applied<'s>, SheetThermalError> {
    let mut temperatures = Vec::with_capacity(surfaces.len());
    for surface in surfaces {
        checkpoint(cancelled)?;
        temperatures.push(match *surface {
            SheetThermalSurface::Stored(id) => Temperature::new(state.temperature_k(id)?),
            SheetThermalSurface::Prescribed(t) => t,
        });
    }
    let spectral_sources = matches!(
        optics, SheetThermalOptics::Spectral { optics, .. }
            if optics.surface_emissivity_bands().is_some()
    );
    let mut cards = Vec::new();
    if !spectral_sources {
        for (surface, &temperature) in temperatures.iter().enumerate() {
            if paired[surface] {
                continue;
            }
            checkpoint(cancelled)?;
            let card = emissivity(surface, temperature)
                .map_err(|error| SheetThermalError::Emissivity { surface, error })?;
            if card.temperature_k() != temperature.value() {
                return Err(SheetThermalError::Invalid(
                    "emissivity must query the actual endpoint temperature",
                ));
            }
            cards.push(card);
        }
    }
    let radiation = match optics {
        SheetThermalOptics::Gray {
            absorption,
            sheets,
            solve,
        } => WallRadiation::Gray(enclosure.evaluate_with_sheets(
            &temperatures,
            &cards,
            sheets,
            gas_temperatures,
            absorption,
            solve,
            &mut work.gray,
            &mut *cancelled,
        )?),
        SheetThermalOptics::Spectral { optics, solve } => {
            WallRadiation::Spectral(enclosure.evaluate_spectral(
                &temperatures,
                &cards,
                gas_temperatures,
                optics,
                solve,
                &mut work.spectral,
                &mut *cancelled,
            )?)
        }
    };
    let mut sums = BTreeMap::new();
    let mut bath_power = vec![Power::new(0.0); surfaces.len()];
    for (i, (&surface, &power)) in surfaces
        .iter()
        .zip(radiation.surface_heat_power())
        .enumerate()
    {
        checkpoint(cancelled)?;
        match surface {
            SheetThermalSurface::Stored(id) => {
                let entry = sums.entry(id).or_insert(Dd::from_f64(0.0));
                *entry = *entry + Dd::from_f64(power.value());
            }
            SheetThermalSurface::Prescribed(_) => bath_power[i] = power,
        }
    }
    pack_applied(radiation, temperatures, sums, bath_power)
}

/// Finalize the exact sources consumed by the native thermal residual once.
pub(crate) fn pack_applied<'s>(
    radiation: WallRadiation<'s>,
    temperatures: Vec<Temperature>,
    sums: BTreeMap<CellId, Dd>,
    bath_power: Vec<Power>,
) -> Result<Applied<'s>, SheetThermalError> {
    let cell_power: BTreeMap<CellId, Power> = sums
        .into_iter()
        .map(|(id, sum)| Ok((id, Power::new(finite(sum.to_f64())?))))
        .collect::<Result<_, DepositionError>>()?;
    let mut balance = Dd::from_f64(0.0);
    for power in cell_power.values() {
        balance = balance + Dd::from_f64(power.value());
    }
    for power in bath_power.iter().chain(radiation.gas_heat_power().iter()) {
        balance = balance + Dd::from_f64(power.value());
    }
    let power_residual_w = finite(balance.to_f64())?;
    Ok(Applied { radiation, temperatures, cell_power, bath_power, power_residual_w })
}

#[allow(clippy::too_many_arguments)]
fn finish<'m, 's>(
    trial: BackwardEulerTrial<'m, Applied<'s>>,
    dt: f64,
    iterations: usize,
    evaluations: usize,
    maximum_cell_residual_j: f64,
    controls: SheetThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SheetThermalCandidate<'m, 's>, SheetThermalError> {
    let applied_power_residual_w = trial.payload().power_residual_w;
    // Preserve the native compensated storage residual instead of subtracting
    // large rounded energy totals or assigning the optical defect to a gas layer.
    let combined = finite(trial.audit().residual_j + finite(dt * applied_power_residual_w)?)?;
    if combined.abs() > controls.combined_energy_tolerance_j {
        return Err(SheetThermalError::EnergyResidual {
            residual_j: combined,
            allowed_j: controls.combined_energy_tolerance_j,
        });
    }
    let (thermal, applied, audit) = trial.accept(
        controls.cell_energy_tolerance_j,
        controls.total_energy_tolerance_j,
        || checkpoint(cancelled),
    )?;
    Ok(SheetThermalCandidate {
        thermal,
        report: SheetThermalReport {
            elapsed: Time::new(dt),
            iterations,
            evaluations,
            radiation: applied.radiation,
            surface_temperatures: applied.temperatures,
            cell_radiant_power: applied.cell_power,
            prescribed_surface_heat_power: applied.bath_power,
            audit,
            maximum_cell_residual_j,
            applied_power_residual_w,
            combined_energy_residual_j: combined,
        },
    })
}

fn trial_domain_error(error: &SheetThermalError) -> bool {
    matches!(
        error,
        SheetThermalError::Thermal(DepositionError::MaterialDomain | DepositionError::Arithmetic)
            | SheetThermalError::Radiation(
                RadiationError::Invalid(_) | RadiationError::SpectralTruncation { .. }
            )
            | SheetThermalError::Emissivity { .. }
    )
}
pub(crate) fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
pub(crate) fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
pub(crate) fn checkpoint(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
