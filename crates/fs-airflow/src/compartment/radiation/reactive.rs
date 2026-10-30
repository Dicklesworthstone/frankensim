//! Participating radiation coupled to original inert walls and receding fuel.
//!
//! The original enclosure feedback iteration is shared with the wall-only API.
//! An inert wall's unknown is its incident net radiation. A reactive slab's
//! unknown is TOTAL front-face heat: radiation plus convection and explicitly
//! supplied external irradiation. Its original implicit reaction/conduction
//! and conservative remap supply the endpoint; current half-cell resistance
//! recovers the actual moving-face temperature. The shared iteration closes
//! that temperature against both radiation and convection.
//!
//! Accepted radiation and convection are applied with identical signed powers
//! to condensed owners and gas. Their query-versus-endpoint mismatch is checked
//! explicitly, never repaired through a gas source. Volatile mass and enthalpy
//! remain in the ORIGINAL slab report; transfer them through `VolatileMixture`
//! separately, exactly once. `report.add_gas_heat` adds only radiation/convection.
//! The combined audit includes exported volatile enthalpy as an outgoing stream.
//!
//! The slab back is insulated. Gas state and supplied optical paths/opacity are
//! frozen over each interval. Geometry recession does not imply a change of gas
//! volume or pressure work. Original condensed/source domains and first-order
//! time integration remain; this is not a flame, soot-production or gas-moving-
//! boundary model. Candidates publish only alongside network acceptance.

use super::ParticipatingEnclosure;
use super::walls::{
    AppliedExchange, EnclosureWallError, EnclosureWallReport, EnclosureWallSolve,
    FeedbackWorkspace, SurfaceTrial, WallOptics, admit_contact_temperatures, wall_error,
    wall_trial,
};
use crate::compartment::wall::{GasContact, Wall1d, WallBoundary, WallCandidate, WallWorkspace};
use fs_conduction::ConductionError;
use fs_conduction::radiation::SurfaceEmissivity;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_qty::{Energy, Power, Temperature, Time};
use fs_thermal::enthalpy::pyrolysis::PyrolysisControls;
use fs_thermal::enthalpy::pyrolysis::implicit::{ImplicitSolve, ImplicitWorkspace};
use fs_thermal::enthalpy::pyrolysis::slab::{RecedingSlab, SlabError, SlabReport, SlabSolve};
use fs_thermal::enthalpy::{Boundary, DepositionError};

/// Explicit front contacts and external input. The slab back is insulated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabBoundary {
    /// Original lower/upper gas contacts on the common moving front face.
    pub gas: [GasContact; 2],
    /// Prescribed absorbed external irradiation, positive into the front [W].
    /// This is separate from participating enclosure radiation.
    pub external_radiant_power: Power,
}

/// Original physical and nonlinear budgets, applied to every slab trial.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SlabSurfaceSolve {
    /// Original time, mass and energy controls.
    pub pyrolysis: PyrolysisControls,
    /// Original implicit conduction/reaction work controls.
    pub implicit: ImplicitSolve,
    /// Original bounded recession/remapping controls.
    pub geometry: SlabSolve,
}

/// Exactly one original condensed owner per enclosure surface.
#[derive(Clone, Copy)]
pub enum ThermalSurface<'a, 'm> {
    /// Original finite wall with inner gas contacts and outer ambient.
    Wall {
        /// Accepted original wall, unchanged by candidate evaluation.
        wall: &'a Wall1d<'m>,
        /// Original wall boundary inputs.
        boundary: WallBoundary,
    },
    /// Original homogeneous-material slab with optional char and an insulated back.
    Slab {
        /// Accepted original receding solid.
        slab: &'a RecedingSlab<'m>,
        /// Explicit front-face contacts and absorbed external irradiation.
        boundary: SlabBoundary,
        /// Original thermal, reaction and geometry controls.
        solve: SlabSurfaceSolve,
    },
}
impl ThermalSurface<'_, '_> {
    fn area(&self) -> f64 {
        match self {
            Self::Wall { wall, .. } => wall.geometry().area.value(),
            Self::Slab { slab, .. } => slab.geometry().area_m2,
        }
    }
    fn time(&self) -> Time {
        match self {
            Self::Wall { wall, .. } => wall.time(),
            Self::Slab { slab, .. } => Time::new(slab.time_s()),
        }
    }
    fn gas(&self) -> &[GasContact; 2] {
        match self {
            Self::Wall { boundary, .. } => &boundary.gas,
            Self::Slab { boundary, .. } => &boundary.gas,
        }
    }
}

/// Accepted original slab and its distinct heat and volatile accounts.
pub struct SlabSurfaceCandidate<'m> {
    /// Install only with the matching accepted gas/network interval.
    pub slab: RecedingSlab<'m>,
    /// Unmodified original integrated release, storage and remap diagnostics.
    pub report: SlabReport,
    /// Actual moving face temperature from returned state and half-cell law.
    pub surface_temperature: Temperature,
    /// Applied convection into gas [W], at the accepted optical-query face.
    /// Its difference from actual-endpoint convection enters the shared gate.
    pub gas_heat_power: [Power; 2],
    /// Identical enclosure surface power integrated over the actual interval.
    pub coupled_radiant_energy: Energy,
    /// Separately prescribed external radiant energy, positive into solid.
    pub external_radiant_energy: Energy,
}

/// Original owner candidates, retaining the supplied enclosure ordering.
pub enum ThermalSurfaceCandidate<'m> {
    /// Original inert wall candidate.
    Wall(WallCandidate<'m>),
    /// Original receding reactive solid candidate.
    Slab(SlabSurfaceCandidate<'m>),
}

/// One jointly accepted optical/thermal candidate; gas itself is unchanged.
pub struct EnclosureSurfaceCandidate<'m, 's> {
    /// Matching candidate for each supplied original surface owner.
    pub surfaces: Vec<ThermalSurfaceCandidate<'m>>,
    /// Shared gas heat, volatile export and original conservation diagnostics.
    pub report: EnclosureWallReport<'s>,
    /// Numerical next-call guess in the SAME surface ordering. For walls this
    /// is radiation; for slabs it is total front heat. It is never an extra source.
    pub input_powers: Vec<Power>,
}

enum SurfaceWorkspace {
    Wall(WallWorkspace),
    Slab(ImplicitWorkspace),
}

/// Prepared original per-owner workspaces and the one shared feedback solver.
pub struct EnclosureSurfaceWorkspace {
    surfaces: Vec<SurfaceWorkspace>,
    feedback: FeedbackWorkspace,
}

impl ParticipatingEnclosure {
    /// Prepare original bounded wall/slab and optical workspaces once.
    pub fn surface_workspace(
        &self,
        surfaces: &[ThermalSurface<'_, '_>],
    ) -> Result<EnclosureSurfaceWorkspace, EnclosureWallError> {
        self.admit_surfaces(surfaces)?;
        let mut work = Vec::with_capacity(surfaces.len());
        for (i, surface) in surfaces.iter().enumerate() {
            work.push(match surface {
                ThermalSurface::Wall { wall, .. } => {
                    SurfaceWorkspace::Wall(wall.workspace().map_err(|e| wall_error(i, e))?)
                }
                ThermalSurface::Slab { slab, .. } => SurfaceWorkspace::Slab(
                    slab.state()
                        .implicit_workspace()
                        .map_err(|e| slab_error(i, SlabError::from(e)))?,
                ),
            });
        }
        Ok(EnclosureSurfaceWorkspace {
            surfaces: work,
            feedback: FeedbackWorkspace::new(self)?,
        })
    }

    /// Use actual original room temperatures and clocks; absent layers may not
    /// own nonzero optical paths or active convective contact.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_room_surfaces<'m, 's>(
        &self,
        room: &TwoZone,
        surfaces: &[ThermalSurface<'_, 'm>],
        duration: Time,
        optics: WallOptics<'_, 's>,
        emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        initial_input_powers: Option<&[Power]>,
        controls: EnclosureWallSolve,
        work: &mut EnclosureSurfaceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnclosureSurfaceCandidate<'m, 's>, EnclosureWallError> {
        poll(&mut cancelled)?;
        if surfaces.iter().any(|s| s.time() != room.time()) {
            return Err(EnclosureWallError::Invalid(
                "surface and room clocks differ",
            ));
        }
        for layer in [Layer::Lower, Layer::Upper] {
            if room.layer(layer).volume().value() == 0.0
                && (self.paths().iter().any(|p| {
                    p.segments
                        .iter()
                        .any(|s| s.layer == layer && s.length.value() > 0.0)
                }) || surfaces.iter().any(|s| {
                    let c = s.gas()[layer.index()];
                    c.area.value() > 0.0 && c.htc.value() > 0.0
                }))
            {
                return Err(EnclosureWallError::Invalid(
                    "radiative or convective contact crosses an absent gas layer",
                ));
            }
        }
        self.candidate_surfaces(
            surfaces,
            duration,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            optics,
            emissivity,
            initial_input_powers,
            controls,
            work,
            cancelled,
        )
    }

    /// Coupled original condensed owners at explicitly supplied frozen gas
    /// temperatures. Each candidate uses the common actual duration and signed
    /// exchanges; all original sources and physical domains remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_surfaces<'m, 's>(
        &self,
        surfaces: &[ThermalSurface<'_, 'm>],
        duration: Time,
        gas_temperatures: [Temperature; 2],
        optics: WallOptics<'_, 's>,
        emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        initial_input_powers: Option<&[Power]>,
        controls: EnclosureWallSolve,
        work: &mut EnclosureSurfaceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EnclosureSurfaceCandidate<'m, 's>, EnclosureWallError> {
        poll(&mut cancelled)?;
        self.admit_surfaces(surfaces)?;
        if work.surfaces.len() != surfaces.len() {
            return Err(EnclosureWallError::Invalid("surface workspace axis"));
        }
        for surface in surfaces {
            admit_contact_temperatures(surface.gas(), gas_temperatures)?;
        }
        let result = self.candidate_feedback(
            duration,
            gas_temperatures,
            optics,
            emissivity,
            initial_input_powers,
            controls,
            &mut work.feedback,
            |i, power, applied, cancel| match (&surfaces[i], &mut work.surfaces[i]) {
                (ThermalSurface::Wall { wall, boundary }, SurfaceWorkspace::Wall(work)) => {
                    let candidate = wall
                        .candidate_with_exchange(
                            duration,
                            *boundary,
                            Power::new(power),
                            controls.wall,
                            work,
                            cancel,
                        )
                        .map_err(|e| wall_error(i, e))?;
                    let trial = wall_trial(candidate);
                    Ok(SurfaceTrial {
                        candidate: ThermalSurfaceCandidate::Wall(trial.candidate),
                        temperature: trial.temperature,
                        elapsed: trial.elapsed,
                        gas_heat_power: trial.gas_heat_power,
                        storage_change: trial.storage_change,
                        volatile_enthalpy: trial.volatile_enthalpy,
                        ambient_heat_energy: trial.ambient_heat_energy,
                        external_radiant_energy: trial.external_radiant_energy,
                        feedback_power: trial.feedback_power,
                    })
                }
                (
                    ThermalSurface::Slab {
                        slab,
                        boundary,
                        solve,
                    },
                    SurfaceWorkspace::Slab(work),
                ) => slab_trial(
                    i, slab, *boundary, *solve, duration, power, applied, work, cancel,
                ),
                _ => Err(EnclosureWallError::Invalid("surface workspace owner kind")),
            },
            cancelled,
        )?;
        Ok(EnclosureSurfaceCandidate {
            surfaces: result.surfaces,
            report: result.report,
            input_powers: result.input_powers,
        })
    }

    fn admit_surfaces(
        &self,
        surfaces: &[ThermalSurface<'_, '_>],
    ) -> Result<(), EnclosureWallError> {
        if surfaces.is_empty() || surfaces.len() != self.geometry().len() {
            return Err(EnclosureWallError::Invalid(
                "one original owner per enclosure surface",
            ));
        }
        let time = surfaces[0].time();
        for (surface, area) in surfaces.iter().zip(self.geometry().areas_m2()) {
            if surface.time() != time || (surface.area() - area).abs() > 32.0 * f64::EPSILON * area
            {
                return Err(EnclosureWallError::Invalid(
                    "surface areas or clocks differ",
                ));
            }
            if let ThermalSurface::Slab { slab, boundary, .. } = surface {
                let area_sum: f64 = boundary.gas.iter().map(|c| c.area.value()).sum();
                if !boundary.external_radiant_power.value().is_finite()
                    || boundary.external_radiant_power.value() < 0.0
                    || !area_sum.is_finite()
                    || area_sum > area * (1.0 + 32.0 * f64::EPSILON)
                    || boundary.gas.iter().any(|c| {
                        !c.area.value().is_finite()
                            || c.area.value() < 0.0
                            || !c.htc.value().is_finite()
                            || c.htc.value() < 0.0
                            || !c.temperature.value().is_finite()
                            || c.temperature.value() <= 0.0
                    })
                {
                    return Err(EnclosureWallError::Invalid(
                        "slab contact area, coefficient or irradiation",
                    ));
                }
                if slab.front_cell_id() != slab.back_cell_id() {
                    let back = slab.state().cells()[&slab.back_cell_id()].boundary;
                    if back.htc_w_m2_k != 0.0
                        || back.source_w != 0.0
                        || back.reservoir.is_some_and(|r| r.conductance_w_k != 0.0)
                    {
                        return Err(EnclosureWallError::Invalid(
                            "radiating slab requires insulated back",
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EnclosureWallError> {
    if cancelled() {
        Err(EnclosureWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn slab_error(surface: usize, error: SlabError) -> EnclosureWallError {
    if error == SlabError::Cancelled {
        EnclosureWallError::Cancelled
    } else {
        EnclosureWallError::Slab { surface, error }
    }
}

fn convective_into_slab(boundary: SlabBoundary, temperature: Temperature) -> [f64; 2] {
    boundary
        .gas
        .map(|c| c.area.value() * c.htc.value() * (c.temperature.value() - temperature.value()))
}

#[allow(clippy::too_many_arguments)]
fn slab_trial<'m>(
    surface: usize,
    slab: &RecedingSlab<'m>,
    boundary: SlabBoundary,
    solve: SlabSurfaceSolve,
    duration: Time,
    power: f64,
    applied: Option<AppliedExchange>,
    work: &mut ImplicitWorkspace,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<SurfaceTrial<ThermalSurfaceCandidate<'m>>, EnclosureWallError> {
    let checkpoint = |cancelled: &mut dyn FnMut() -> bool| {
        if cancelled() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    let mut prepared = slab.clone();
    let front = slab.front_cell_id();
    let temperature = slab
        .state()
        .temperature_k(front)
        .map_err(|e| slab_error(surface, e.into()))?;
    prepared
        .set_boundary(
            front,
            Boundary {
                ambient_k: temperature,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: power,
            },
            || checkpoint(cancelled),
        )
        .map_err(|e| slab_error(surface, e))?;
    let (next, report) = prepared
        .advance_implicit_to(
            slab.time_s() + duration.value(),
            solve.pyrolysis,
            solve.implicit,
            solve.geometry,
            work,
            || checkpoint(cancelled),
        )
        .map_err(|e| slab_error(surface, e))?;
    let elapsed = Time::new(next.time_s() - slab.time_s());
    let dt = elapsed.value();
    let conductance = next
        .front_half_cell_conductance_w_k()
        .map_err(|e| slab_error(surface, e))?;
    let cell_temperature = next
        .state()
        .temperature_k(next.front_cell_id())
        .map_err(|e| slab_error(surface, e.into()))?;
    let face = Temperature::new(cell_temperature + power / conductance);
    if !face.value().is_finite() || face.value() <= 0.0 {
        return Err(EnclosureWallError::Invalid(
            "reactive face temperature outside positive domain",
        ));
    }
    let endpoint_convection = convective_into_slab(boundary, face);
    let external = boundary.external_radiant_power.value();
    let feedback = endpoint_convection[0] + endpoint_convection[1] + external;
    let applied_convection =
        convective_into_slab(boundary, applied.map_or(face, |a| a.temperature));
    let radiation = applied.map_or(power - feedback, |a| a.radiation.value());
    let gas_heat_power = applied_convection.map(|p| Power::new(-p));
    let storage_change = Energy::new(report.pyrolysis.thermal.stored_energy_change_j);
    let volatile_enthalpy = Energy::new(report.pyrolysis.volatile_enthalpy_j);
    let external_radiant_energy = Energy::new(dt * external);
    if ![
        feedback,
        radiation,
        applied_convection[0],
        applied_convection[1],
        dt * radiation,
        external_radiant_energy.value(),
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(EnclosureWallError::Invalid(
            "unrepresentable reactive face exchange",
        ));
    }
    Ok(SurfaceTrial {
        candidate: ThermalSurfaceCandidate::Slab(SlabSurfaceCandidate {
            slab: next,
            report,
            surface_temperature: face,
            gas_heat_power,
            coupled_radiant_energy: Energy::new(dt * radiation),
            external_radiant_energy,
        }),
        temperature: face,
        elapsed,
        gas_heat_power,
        storage_change,
        volatile_enthalpy,
        ambient_heat_energy: Energy::new(0.0),
        external_radiant_energy,
        feedback_power: feedback,
    })
}

#[cfg(test)]
mod tests;
