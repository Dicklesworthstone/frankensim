//! Implicit radiation feedback to the original stored wall profiles.
//!
//! For trial surface powers `q`, the original backward-Euler wall solves give
//! faces `T_w(q)`. The original enclosure then gives `R(T_w(q), T_g)`. This
//! adapter solves `q - R = 0`; it does not duplicate conduction, surface
//! half-cell resistance, radiosity, participating transfer or Planck integration.
//! A Jacobian column perturbs only one wall response before reevaluating the
//! enclosure. The remaining wall profiles need not be solved again for that
//! column. Prepared wall, optical and dense-LU storage is reused.
//!
//! At acceptance the identical enclosure powers are applied to fresh original
//! wall candidates and the gas account. The actual returned wall faces are
//! reevaluated to gate both temperature and radiative-power consistency. These
//! nonlinear residuals remain separate from the original wall, radiation and
//! spectral quadrature residuals; no imbalance is redistributed to a gas layer.
//! Each returned optical report names the temperatures at which its applied
//! powers were evaluated, which agree with the returned faces to the declared
//! tolerance. Material emissivity is resolved at every queried temperature.
//!
//! Gas temperatures, absorption, geometry, convection and external boundaries
//! are supplied for this interval. The gas itself is not advanced here: use the
//! result as one original network source candidate and publish all wall owners
//! only after network acceptance. This gives first-order implicit wall/radiation
//! coupling, not a fully implicit reacting gas solve, an entropy theorem, a
//! spectral database or an unconditional convergence/accuracy guarantee.

use super::spectral::{SpectralOptics, SpectralRadiationReport, SpectralSolve, SpectralWorkspace};
use super::{
    LayerAbsorption, ParticipatingEnclosure, RadiationError, RadiationReport, RadiationSolve,
    RadiationWorkspace,
};
use crate::compartment::wall::{
    Wall1d, WallBoundary, WallCandidate, WallError, WallSolve, WallWorkspace,
};
use fs_conduction::ConductionError;
use fs_conduction::radiation::SurfaceEmissivity;
use fs_gas::multicomponent::compartment::{Layer, LayerRates, TwoZone};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_qty::{Energy, Power, Temperature, Time};

/// Original, source-bearing optical closure, held over the candidate interval.
#[derive(Clone, Copy)]
pub enum WallOptics<'a, 's> {
    /// Original gray participating-medium solve.
    Gray {
        /// Explicit `[lower, upper]` absorption at the supplied gas state.
        absorption: &'a [LayerAbsorption; 2],
        /// Original algebraic gates.
        solve: RadiationSolve,
    },
    /// Original finite-window wavelength-dependent solve.
    Spectral {
        /// Sourced absorption bands and optional original soot populations.
        optics: &'a SpectralOptics<'s>,
        /// Original algebraic, quadrature, coverage and wavelength-work gates.
        solve: SpectralSolve,
    },
}

/// The unmodified original report whose powers were applied to all owners.
#[derive(Clone, Debug, PartialEq)]
pub enum WallRadiation<'s> {
    /// Original gray solve, retaining surface and absorption sources.
    Gray(RadiationReport),
    /// Original spectral solve, retaining sources, coverage and error estimates.
    Spectral(SpectralRadiationReport<'s>),
}
impl WallRadiation<'_> {
    /// Signed powers into the original enclosure's surfaces, in surface order.
    pub fn surface_heat_power(&self) -> &[Power] {
        match self {
            Self::Gray(r) => &r.surface_heat_power,
            Self::Spectral(r) => &r.surface_heat_power,
        }
    }
    /// Signed radiative powers into the `[lower, upper]` gas owners.
    pub fn gas_heat_power(&self) -> [Power; 2] {
        match self {
            Self::Gray(r) => r.gas_heat_power,
            Self::Spectral(r) => r.gas_heat_power,
        }
    }
}

/// Bounded nonlinear work and actual endpoint consistency gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnclosureWallSolve {
    /// Original wall storage/nonlinear controls, applied to every trial.
    pub wall: WallSolve,
    /// Maximum outer Newton iterations, 1 through 64.
    pub maximum_iterations: usize,
    /// Maximum globalization probes per Newton update, 1 through 32.
    pub maximum_backtracks: usize,
    /// Total gray/spectral evaluations, including Jacobians and rejected trials.
    pub maximum_radiation_evaluations: usize,
    /// Total original wall candidate evaluations, including rejected trials.
    pub maximum_wall_evaluations: usize,
    /// Nominal finite-difference surface-power increment [W].
    pub power_difference_w: f64,
    /// Actual applied-versus-endpoint radiative-power gate for every owner [W].
    pub power_tolerance_w: f64,
    /// Applied optical-query versus returned wall-face temperature gate [K].
    pub temperature_tolerance_k: f64,
    /// Actual all-wall storage plus gas heat minus external heat gate [J].
    pub combined_energy_tolerance_j: f64,
}
impl Default for EnclosureWallSolve {
    fn default() -> Self {
        Self {
            wall: WallSolve::default(),
            maximum_iterations: 32,
            maximum_backtracks: 24,
            maximum_radiation_evaluations: 4096,
            maximum_wall_evaluations: 16_384,
            power_difference_w: 0.01,
            power_tolerance_w: 1e-4,
            temperature_tolerance_k: 1e-7,
            combined_energy_tolerance_j: 1e-4,
        }
    }
}

/// Actual work, including discarded numerical probes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EnclosureWallWork {
    /// Calls to the original wall candidate solver.
    pub wall_evaluations: usize,
    /// Calls to the original gray or spectral enclosure solver.
    pub radiation_evaluations: usize,
    /// Wavelength evaluations across all successful spectral trials.
    pub wavelength_evaluations: usize,
    /// Optical matrix factorizations across all successful spectral trials.
    pub spectral_factorizations: usize,
}

/// Exact applied exchanges and measured endpoint consistency of one candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct EnclosureWallReport<'s> {
    /// Actual common interval from the original wall time owner.
    pub elapsed: Time,
    /// Number of outer Newton residual evaluations through acceptance.
    pub iterations: usize,
    /// Unmodified original optical report. Do not separately add its gas powers
    /// after calling this report's `add_gas_heat`, which includes them already.
    pub radiation: WallRadiation<'s>,
    /// Face temperatures queried for the applied radiation powers.
    pub radiation_surface_temperatures: Vec<Temperature>,
    /// Returned face minus applied optical-query temperature [K], by surface.
    pub surface_temperature_residual_k: Vec<f64>,
    /// Radiation reevaluated at returned faces minus applied radiation [W].
    pub surface_power_residual_w: Vec<f64>,
    /// The same actual endpoint mismatch for `[lower, upper]` gas radiation [W].
    pub gas_power_residual_w: [f64; 2],
    /// Largest actual nonlinear surface-power or endpoint power residual [W].
    pub maximum_power_residual_w: f64,
    /// Largest actual optical-query versus returned face mismatch [K].
    pub maximum_temperature_residual_k: f64,
    /// TOTAL power into each gas layer: applied radiation plus original signed
    /// convection from every returned wall candidate. Add exactly once.
    pub gas_heat_power: [Power; 2],
    /// Those identical total gas powers integrated over the accepted interval.
    pub gas_heat_energy: [Energy; 2],
    /// Sum of original wall storage changes [J].
    pub wall_storage_change: Energy,
    /// Sum of original outer ambient exchanges, positive into walls [J].
    pub ambient_heat_energy: Energy,
    /// Sum of separately supplied external radiant input, positive into walls [J].
    pub external_radiant_energy: Energy,
    /// Wall storage plus total gas heat minus ambient/external input [J].
    pub combined_energy_residual_j: f64,
    /// Bounded actual work used by the complete nonlinear call.
    pub work: EnclosureWallWork,
}
impl EnclosureWallReport<'_> {
    /// Add total radiation-plus-convection gas heat atomically. Matching wall
    /// candidates already contain the opposite internal exchanges.
    pub fn add_gas_heat(&self, rates: &mut [LayerRates; 2]) -> Result<(), RadiationError> {
        let mut next = *rates;
        for (rate, heat) in next.iter_mut().zip(self.gas_heat_power) {
            rate.add_heat(heat)?;
        }
        *rates = next;
        Ok(())
    }
}

/// All original wall candidates and their matching common gas exchange.
pub struct EnclosureWallCandidate<'m, 's> {
    /// Same ordering as the admitted enclosure surfaces. Install only after
    /// the enclosing gas/network step has accepted this candidate's interval.
    pub walls: Vec<WallCandidate<'m>>,
    /// Applied optical/thermal accounts and numerical consistency checks.
    pub report: EnclosureWallReport<'s>,
}

/// Reusable numerical state; it owns no accepted thermal or gas history.
pub struct EnclosureWallWorkspace {
    walls: Vec<WallWorkspace>,
    gray: RadiationWorkspace,
    spectral: SpectralWorkspace,
    power: Vec<f64>,
    probe: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    delta: Vec<f64>,
    lu: LuWorkspace,
}

/// Atomic source, geometry, numerical-work or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum EnclosureWallError {
    /// Invalid axes, geometry, clocks, contact temperatures or controls.
    Invalid(&'static str),
    /// Original wall/material refusal at the named enclosure surface.
    Wall { surface: usize, error: WallError },
    /// Original material-card emissivity query refused at this surface.
    Emissivity {
        surface: usize,
        error: ConductionError,
    },
    /// Original gray or spectral radiation owner refused.
    Radiation(RadiationError),
    /// Original prepared dense solver refused the coupling Jacobian.
    Linear(LuWorkspaceError),
    /// Explicit total evaluation budget exhausted.
    Budget {
        owner: &'static str,
        evaluations: usize,
        maximum: usize,
    },
    /// Bounded Newton/globalization did not close the actual residual.
    Iterations,
    /// Actual combined storage/external account exceeded its supplied gate.
    EnergyResidual { residual_j: f64, allowed_j: f64 },
    /// No candidate was published.
    Cancelled,
}
impl core::fmt::Display for EnclosureWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled compartment walls/radiation: {self:?}")
    }
}
impl std::error::Error for EnclosureWallError {}
impl From<RadiationError> for EnclosureWallError {
    fn from(error: RadiationError) -> Self {
        if error == RadiationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Radiation(error)
        }
    }
}
impl From<LuWorkspaceError> for EnclosureWallError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Linear(error)
    }
}

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), EnclosureWallError> {
    if cancelled() {
        Err(EnclosureWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, EnclosureWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnclosureWallError::Invalid(
            "non-finite coupling arithmetic",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn fatal_probe(error: &EnclosureWallError) -> bool {
    matches!(
        error,
        EnclosureWallError::Cancelled
            | EnclosureWallError::Budget { .. }
            | EnclosureWallError::Radiation(RadiationError::SpectralBudget { .. })
    )
}
fn charge(
    count: &mut usize,
    maximum: usize,
    owner: &'static str,
) -> Result<(), EnclosureWallError> {
    if *count >= maximum {
        return Err(EnclosureWallError::Budget {
            owner,
            evaluations: *count,
            maximum,
        });
    }
    *count += 1;
    Ok(())
}

impl ParticipatingEnclosure {
    /// Prepare the original per-wall and enclosure buffers plus the bounded
    /// coupling Jacobian. Wall reports and native thermal trials retain their
    /// original source-bearing allocations; this is not an allocation-free solve.
    pub fn wall_workspace(
        &self,
        walls: &[Wall1d<'_>],
    ) -> Result<EnclosureWallWorkspace, EnclosureWallError> {
        self.admit_walls(walls)?;
        let n = walls.len();
        let wall_work = walls
            .iter()
            .enumerate()
            .map(|(surface, wall)| {
                wall.workspace()
                    .map_err(|error| EnclosureWallError::Wall { surface, error })
            })
            .collect::<Result<_, _>>()?;
        Ok(EnclosureWallWorkspace {
            walls: wall_work,
            gray: self.workspace()?,
            spectral: self.spectral_workspace()?,
            power: vec![0.0; n],
            probe: vec![0.0; n],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            delta: vec![0.0; n],
            lu: LuWorkspace::new(n)?,
        })
    }

    /// Use actual room temperatures and clocks, refusing optical/convective
    /// paths through an absent layer. Gas state itself remains unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_room_walls<'m, 's>(
        &self,
        room: &TwoZone,
        walls: &[Wall1d<'m>],
        duration: Time,
        boundaries: &[WallBoundary],
        optics: WallOptics<'_, 's>,
        emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        initial_surface_power: Option<&[Power]>,
        controls: EnclosureWallSolve,
        work: &mut EnclosureWallWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnclosureWallCandidate<'m, 's>, EnclosureWallError> {
        check(&mut cancelled)?;
        if walls.iter().any(|w| w.time() != room.time()) {
            return Err(EnclosureWallError::Invalid("wall and room clocks differ"));
        }
        for layer in [Layer::Lower, Layer::Upper] {
            if room.layer(layer).volume().value() == 0.0 {
                if self.paths().iter().any(|p| {
                    p.segments
                        .iter()
                        .any(|s| s.layer == layer && s.length.value() > 0.0)
                }) || boundaries.iter().any(|b| {
                    b.gas[layer.index()].area.value() > 0.0
                        && b.gas[layer.index()].htc.value() > 0.0
                }) {
                    return Err(EnclosureWallError::Invalid(
                        "radiative or convective contact crosses an absent gas layer",
                    ));
                }
            }
        }
        self.candidate_walls(
            walls,
            duration,
            boundaries,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            optics,
            emissivity,
            initial_surface_power,
            controls,
            work,
            cancelled,
        )
    }

    /// Solve all wall faces and their participating radiation together.
    ///
    /// One original wall must own every enclosure surface, with the same area
    /// and accepted clock. Nonzero convective contacts must use the supplied
    /// gas temperatures. The emissivity callback resolves the original card at
    /// the requested temperature and must be deterministic for unchanged input.
    /// It is called for every optical trial, with bounded cancellation between
    /// surfaces. Opacity/bands, gas temperatures and geometry remain frozen.
    ///
    /// An optional initial power vector is only a nonlinear guess, never a
    /// prescribed extra source. `None` starts from zero radiation. A useful
    /// previous-step guess can admit cases whose isolated zero-radiation wall
    /// endpoint lies outside its material domain. Failed calls mutate scratch
    /// only; subsequent calls always start from the caller's fresh guess.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_walls<'m, 's>(
        &self,
        walls: &[Wall1d<'m>],
        duration: Time,
        boundaries: &[WallBoundary],
        gas_temperatures: [Temperature; 2],
        optics: WallOptics<'_, 's>,
        mut emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        initial_surface_power: Option<&[Power]>,
        controls: EnclosureWallSolve,
        work: &mut EnclosureWallWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnclosureWallCandidate<'m, 's>, EnclosureWallError> {
        check(&mut cancelled)?;
        self.admit_walls(walls)?;
        let n = walls.len();
        if boundaries.len() != n
            || work.walls.len() != n
            || work.power.len() != n
            || !positive(duration.value())
            || gas_temperatures.iter().any(|t| !positive(t.value()))
            || !(1..=64).contains(&controls.maximum_iterations)
            || !(1..=32).contains(&controls.maximum_backtracks)
            || !(1..=65_536).contains(&controls.maximum_radiation_evaluations)
            || !(1..=65_536).contains(&controls.maximum_wall_evaluations)
            || !positive(controls.power_difference_w)
            || !nonnegative(controls.power_tolerance_w)
            || !nonnegative(controls.temperature_tolerance_k)
            || !nonnegative(controls.combined_energy_tolerance_j)
            || initial_surface_power
                .is_some_and(|p| p.len() != n || p.iter().any(|q| !q.value().is_finite()))
        {
            return Err(EnclosureWallError::Invalid(
                "wall coupling axes, duration, workspace, guess or controls",
            ));
        }
        for boundary in boundaries {
            for (contact, gas_t) in boundary.gas.iter().zip(gas_temperatures) {
                if contact.area.value() > 0.0
                    && contact.htc.value() > 0.0
                    && contact.temperature != gas_t
                {
                    return Err(EnclosureWallError::Invalid(
                        "convective and radiative gas temperatures differ",
                    ));
                }
            }
        }
        for (i, q) in work.power.iter_mut().enumerate() {
            *q = initial_surface_power.map_or(0.0, |p| p[i].value());
        }
        let mut counts = EnclosureWallWork::default();
        for iteration in 0..controls.maximum_iterations {
            check(&mut cancelled)?;
            let base_walls = solve_walls(
                walls,
                duration,
                boundaries,
                &work.power,
                controls,
                &mut work.walls,
                &mut counts,
                &mut cancelled,
            )?;
            let temperatures: Vec<_> = base_walls
                .iter()
                .map(|w| w.report.surface_temperatures[0])
                .collect();
            let radiation = optical(
                self,
                &temperatures,
                gas_temperatures,
                optics,
                &mut emissivity,
                controls,
                &mut work.gray,
                &mut work.spectral,
                &mut counts,
                &mut cancelled,
            )?;
            let mut norm = 0.0_f64;
            for i in 0..n {
                work.rhs[i] = finite(radiation.surface_heat_power()[i].value() - work.power[i])?;
                norm = norm.max(work.rhs[i].abs());
            }
            if norm <= controls.power_tolerance_w {
                // Apply this ONE conservative optical result to both owners.
                for i in 0..n {
                    work.probe[i] = radiation.surface_heat_power()[i].value();
                }
                let accepted = solve_walls(
                    walls,
                    duration,
                    boundaries,
                    &work.probe,
                    controls,
                    &mut work.walls,
                    &mut counts,
                    &mut cancelled,
                )?;
                let endpoint: Vec<_> = accepted
                    .iter()
                    .map(|w| w.report.surface_temperatures[0])
                    .collect();
                let endpoint_radiation = optical(
                    self,
                    &endpoint,
                    gas_temperatures,
                    optics,
                    &mut emissivity,
                    controls,
                    &mut work.gray,
                    &mut work.spectral,
                    &mut counts,
                    &mut cancelled,
                )?;
                let temperature_residual: Vec<_> = endpoint
                    .iter()
                    .zip(&temperatures)
                    .map(|(a, b)| a.value() - b.value())
                    .collect();
                let surface_residual: Vec<_> = endpoint_radiation
                    .surface_heat_power()
                    .iter()
                    .zip(radiation.surface_heat_power())
                    .map(|(a, b)| a.value() - b.value())
                    .collect();
                let gas_residual = std::array::from_fn(|i| {
                    endpoint_radiation.gas_heat_power()[i].value()
                        - radiation.gas_heat_power()[i].value()
                });
                let max_temperature = temperature_residual
                    .iter()
                    .map(|r| r.abs())
                    .fold(0.0, f64::max);
                let max_power = surface_residual
                    .iter()
                    .chain(gas_residual.iter())
                    .map(|r| r.abs())
                    .fold(norm, f64::max);
                if max_temperature <= controls.temperature_tolerance_k
                    && max_power <= controls.power_tolerance_w
                {
                    return finish(
                        accepted,
                        radiation,
                        temperatures,
                        temperature_residual,
                        surface_residual,
                        gas_residual,
                        max_temperature,
                        max_power,
                        iteration + 1,
                        counts,
                        controls,
                        &mut cancelled,
                    );
                }
            }
            // Each radiation-power column needs only its own perturbed wall.
            for j in 0..n {
                check(&mut cancelled)?;
                let increment = controls
                    .power_difference_w
                    .max(64.0 * f64::EPSILON * work.power[j].abs());
                let mut column = None;
                for sign in [1.0, -1.0] {
                    let q = finite(work.power[j] + sign * increment)?;
                    let actual_increment = q - work.power[j];
                    if actual_increment == 0.0 {
                        return Err(EnclosureWallError::Invalid(
                            "unrepresentable surface-power increment",
                        ));
                    }
                    let trial = solve_wall(
                        j,
                        &walls[j],
                        duration,
                        boundaries[j],
                        q,
                        controls,
                        &mut work.walls[j],
                        &mut counts,
                        &mut cancelled,
                    );
                    let trial = match trial {
                        Ok(t) => t,
                        Err(e) if fatal_probe(&e) => return Err(e),
                        Err(_) => continue,
                    };
                    let mut probe_temperatures = temperatures.clone();
                    probe_temperatures[j] = trial.report.surface_temperatures[0];
                    match optical(
                        self,
                        &probe_temperatures,
                        gas_temperatures,
                        optics,
                        &mut emissivity,
                        controls,
                        &mut work.gray,
                        &mut work.spectral,
                        &mut counts,
                        &mut cancelled,
                    ) {
                        Ok(r) => {
                            column = Some((r, actual_increment));
                            break;
                        }
                        Err(e) if fatal_probe(&e) => return Err(e),
                        Err(_) => {}
                    }
                }
                let (probe, increment) = column.ok_or(EnclosureWallError::Iterations)?;
                for i in 0..n {
                    let derivative = finite(
                        (probe.surface_heat_power()[i].value()
                            - radiation.surface_heat_power()[i].value())
                            / increment,
                    )?;
                    work.matrix[i * n + j] = finite(if i == j {
                        1.0 - derivative
                    } else {
                        -derivative
                    })?;
                }
            }
            check(&mut cancelled)?;
            work.lu
                .solve_into(&work.matrix, &work.rhs, &mut work.delta)?;
            let mut scale = 1.0;
            let mut found = false;
            for _ in 0..controls.maximum_backtracks {
                check(&mut cancelled)?;
                for i in 0..n {
                    work.probe[i] = finite(work.power[i] + scale * work.delta[i])?;
                }
                let trial = solve_walls(
                    walls,
                    duration,
                    boundaries,
                    &work.probe,
                    controls,
                    &mut work.walls,
                    &mut counts,
                    &mut cancelled,
                );
                let trial = match trial {
                    Ok(t) => t,
                    Err(e) if fatal_probe(&e) => return Err(e),
                    Err(_) => {
                        scale *= 0.5;
                        continue;
                    }
                };
                let t: Vec<_> = trial
                    .iter()
                    .map(|w| w.report.surface_temperatures[0])
                    .collect();
                match optical(
                    self,
                    &t,
                    gas_temperatures,
                    optics,
                    &mut emissivity,
                    controls,
                    &mut work.gray,
                    &mut work.spectral,
                    &mut counts,
                    &mut cancelled,
                ) {
                    Ok(r) => {
                        let next_norm = r
                            .surface_heat_power()
                            .iter()
                            .zip(&work.probe)
                            .map(|(r, q)| (q - r.value()).abs())
                            .fold(0.0, f64::max);
                        if next_norm < norm {
                            found = true;
                            break;
                        }
                    }
                    Err(e) if fatal_probe(&e) => return Err(e),
                    Err(_) => {}
                }
                scale *= 0.5;
            }
            if !found {
                return Err(EnclosureWallError::Iterations);
            }
            work.power.copy_from_slice(&work.probe);
        }
        Err(EnclosureWallError::Iterations)
    }

    fn admit_walls(&self, walls: &[Wall1d<'_>]) -> Result<(), EnclosureWallError> {
        if walls.len() != self.geometry().len() || walls.is_empty() {
            return Err(EnclosureWallError::Invalid(
                "one original wall required per enclosure surface",
            ));
        }
        let time = walls[0].time();
        for (wall, area) in walls.iter().zip(self.geometry().areas_m2()) {
            if wall.time() != time
                || (wall.geometry().area.value() - area).abs() > 32.0 * f64::EPSILON * area
            {
                return Err(EnclosureWallError::Invalid(
                    "wall and enclosure areas or wall clocks differ",
                ));
            }
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_wall<'m>(
    surface: usize,
    wall: &Wall1d<'m>,
    duration: Time,
    boundary: WallBoundary,
    power: f64,
    controls: EnclosureWallSolve,
    work: &mut WallWorkspace,
    counts: &mut EnclosureWallWork,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WallCandidate<'m>, EnclosureWallError> {
    check(cancelled)?;
    charge(
        &mut counts.wall_evaluations,
        controls.maximum_wall_evaluations,
        "wall",
    )?;
    wall.candidate_with_exchange(
        duration,
        boundary,
        Power::new(power),
        controls.wall,
        work,
        cancelled,
    )
    .map_err(|error| {
        if error == WallError::Cancelled {
            EnclosureWallError::Cancelled
        } else {
            EnclosureWallError::Wall { surface, error }
        }
    })
}

#[allow(clippy::too_many_arguments)]
fn solve_walls<'m>(
    walls: &[Wall1d<'m>],
    duration: Time,
    boundaries: &[WallBoundary],
    powers: &[f64],
    controls: EnclosureWallSolve,
    work: &mut [WallWorkspace],
    counts: &mut EnclosureWallWork,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<WallCandidate<'m>>, EnclosureWallError> {
    walls
        .iter()
        .enumerate()
        .map(|(i, w)| {
            solve_wall(
                i,
                w,
                duration,
                boundaries[i],
                powers[i],
                controls,
                &mut work[i],
                counts,
                cancelled,
            )
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn optical<'s>(
    enclosure: &ParticipatingEnclosure,
    temperatures: &[Temperature],
    gas: [Temperature; 2],
    optics: WallOptics<'_, 's>,
    emissivity: &mut impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
    controls: EnclosureWallSolve,
    gray: &mut RadiationWorkspace,
    spectral: &mut SpectralWorkspace,
    counts: &mut EnclosureWallWork,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WallRadiation<'s>, EnclosureWallError> {
    check(cancelled)?;
    charge(
        &mut counts.radiation_evaluations,
        controls.maximum_radiation_evaluations,
        "radiation",
    )?;
    let mut emissivities = Vec::with_capacity(temperatures.len());
    for (surface, &temperature) in temperatures.iter().enumerate() {
        check(cancelled)?;
        let value = emissivity(surface, temperature)
            .map_err(|error| EnclosureWallError::Emissivity { surface, error })?;
        if value.temperature_k() != temperature.value() {
            return Err(EnclosureWallError::Invalid(
                "emissivity callback must query the actual requested temperature",
            ));
        }
        emissivities.push(value);
    }
    match optics {
        WallOptics::Gray { absorption, solve } => Ok(WallRadiation::Gray(enclosure.evaluate(
            temperatures,
            &emissivities,
            gas,
            absorption,
            solve,
            gray,
            cancelled,
        )?)),
        WallOptics::Spectral { optics, solve } => {
            let result = enclosure.evaluate_spectral(
                temperatures,
                &emissivities,
                gas,
                optics,
                solve,
                spectral,
                cancelled,
            )?;
            counts.wavelength_evaluations += result.wavelength_evaluations;
            counts.spectral_factorizations += result.factorizations;
            Ok(WallRadiation::Spectral(result))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn finish<'m, 's>(
    walls: Vec<WallCandidate<'m>>,
    radiation: WallRadiation<'s>,
    temperatures: Vec<Temperature>,
    temperature_residual: Vec<f64>,
    surface_residual: Vec<f64>,
    gas_residual: [f64; 2],
    max_temperature: f64,
    max_power: f64,
    iterations: usize,
    counts: EnclosureWallWork,
    controls: EnclosureWallSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EnclosureWallCandidate<'m, 's>, EnclosureWallError> {
    let elapsed = walls[0].report.elapsed;
    let dt = elapsed.value();
    let mut gas = radiation.gas_heat_power().map(Power::value);
    let mut storage = 0.0;
    let mut ambient = 0.0;
    let mut external = 0.0;
    for wall in &walls {
        check(cancelled)?;
        if wall.report.elapsed != elapsed {
            return Err(EnclosureWallError::Invalid(
                "wall candidate intervals differ",
            ));
        }
        for (power, convection) in gas.iter_mut().zip(wall.report.gas_heat_power) {
            *power = finite(*power + convection.value())?;
        }
        storage = finite(storage + wall.report.audit.stored_energy_change_j)?;
        ambient = finite(ambient + wall.report.ambient_heat_energy.value())?;
        external = finite(external + wall.report.external_radiant_energy.value())?;
    }
    let gas_energy = [finite(dt * gas[0])?, finite(dt * gas[1])?];
    let residual = finite(storage + gas_energy[0] + gas_energy[1] - ambient - external)?;
    if residual.abs() > controls.combined_energy_tolerance_j {
        return Err(EnclosureWallError::EnergyResidual {
            residual_j: residual,
            allowed_j: controls.combined_energy_tolerance_j,
        });
    }
    check(cancelled)?;
    Ok(EnclosureWallCandidate {
        walls,
        report: EnclosureWallReport {
            elapsed,
            iterations,
            radiation,
            radiation_surface_temperatures: temperatures,
            surface_temperature_residual_k: temperature_residual,
            surface_power_residual_w: surface_residual,
            gas_power_residual_w: gas_residual,
            maximum_power_residual_w: max_power,
            maximum_temperature_residual_k: max_temperature,
            gas_heat_power: gas.map(Power::new),
            gas_heat_energy: gas_energy.map(Energy::new),
            wall_storage_change: Energy::new(storage),
            ambient_heat_energy: Energy::new(ambient),
            external_radiant_energy: Energy::new(external),
            combined_energy_residual_j: residual,
            work: counts,
        },
    })
}

#[cfg(test)]
mod tests;
