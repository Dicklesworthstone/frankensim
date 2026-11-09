//! Finite reacting discharge through a converging restriction to fixed back pressure.
//!
//! The original material reservoir owns its immutable birth isentrope, remaining
//! inventory and integrated outgoing stagnation enthalpy. This consumer supplies
//! the actual reacting nozzle rate to the same midpoint mass solver as rocket
//! blowdown. It admits choked and subsonic discharge on the original single-sonic
//! branch; no frozen-gamma rate, new reservoir EOS or prescribed mass flux enters.
//!
//! A step crossing the external pressure refuses without clipping. This is a
//! one-direction, well-mixed, rigid adiabatic tank and an external pressure
//! boundary, not a finite receiver, reverse-flow model or resolved jet.

use super::{
    EquilibriumNozzleControl, EquilibriumNozzleError, EquilibriumNozzleFlow,
    publish_on_path, resolve_on_isentrope, validate_flow_inputs,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::{
    EquilibriumBlowdownError, EquilibriumBlowdownSolve, EquilibriumBlowdownStep,
    EquilibriumGasReservoir,
    midpoint::{self, OutflowMidpointError, OutflowMidpointSolve, OutflowRate},
};
use fs_qty::{Area, Pressure};
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumExpansionErrorV1,
    isentrope::density::EquilibriumDensityErrorV1,
};

/// Original source/finite-energy policies plus physical-time mass limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumNozzleBlowdownControl {
    pub reservoir: EquilibriumBlowdownSolve,
    /// Gibbs/entropy policy must exactly match the immutable reservoir birth.
    pub nozzle: EquilibriumNozzleControl,
    /// Maximum fraction of OLD tank mass discharged, in (0,0.5]. No clipping.
    pub maximum_outgoing_fraction: f64,
    /// One to 128 bounded midpoint mass-rate evaluations.
    pub maximum_mass_iterations: usize,
}

/// Complete finite interval; carried H already contains flow work and jet energy.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumNozzleBlowdownStep {
    pub duration_s: f64,
    /// Existing owner's actual mass loss and integrated stagnation enthalpy.
    /// Molecular species inventories may react; their differences are not fluxes.
    pub gas: EquilibriumBlowdownStep,
    /// Full source-bearing physical midpoint inlet and discharge station.
    pub midpoint_nozzle: EquilibriumNozzleFlow,
    pub initial_choked: bool,
    pub endpoint_choked: bool,
    /// Actual expelled mass minus dt times the final published midpoint rate.
    pub mass_quadrature_residual_kg: f64,
    pub mass_iterations: usize,
    /// All density inversion pressure queries, including final reconstruction.
    pub density_pressure_queries: usize,
    /// All initial/trial/midpoint/final nozzle SP queries; preparation is separate.
    pub nozzle_pressure_queries: usize,
    /// One original birth-isentrope preparation for the whole finite interval.
    pub isentrope_preparations: usize,
}

/// The original tank is unchanged until its caller accepts this complete value.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumNozzleBlowdownCandidate {
    endpoint: EquilibriumGasReservoir,
    step: EquilibriumNozzleBlowdownStep,
}
impl EquilibriumNozzleBlowdownCandidate {
    pub const fn endpoint(&self) -> &EquilibriumGasReservoir {
        &self.endpoint
    }
    pub const fn step(&self) -> &EquilibriumNozzleBlowdownStep {
        &self.step
    }
    pub fn into_parts(self) -> (EquilibriumGasReservoir, EquilibriumNozzleBlowdownStep) {
        (self.endpoint, self.step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumNozzleBlowdownError {
    Invalid(&'static str),
    Reservoir(EquilibriumBlowdownError),
    Nozzle(EquilibriumNozzleError),
    /// The requested physical interval exceeds the original outgoing mass cap.
    InventoryLimit {
        maximum_duration_s: f64,
    },
    /// The endpoint would reverse the admitted pressure drive. Refine the
    /// interval or resolve its pressure-equilibrium event; no state is published.
    BackPressureCrossing {
        endpoint_pressure_pa: f64,
        back_pressure_pa: f64,
    },
    NoConvergence {
        mass_iterations: usize,
        mass_residual_kg: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EquilibriumNozzleBlowdownError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite reacting converging-nozzle discharge: {self:?}")
    }
}
impl std::error::Error for EquilibriumNozzleBlowdownError {}
type Error = EquilibriumNozzleBlowdownError;

impl From<EquilibriumBlowdownError> for Error {
    fn from(error: EquilibriumBlowdownError) -> Self {
        match error {
            EquilibriumBlowdownError::Cancelled
            | EquilibriumBlowdownError::Thermochemistry(
                EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled),
            )
            | EquilibriumBlowdownError::Density(EquilibriumDensityErrorV1::Thermochemistry(
                EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled),
            )) => Self::Cancelled,
            other => Self::Reservoir(other),
        }
    }
}
impl From<EquilibriumExpansionErrorV1> for Error {
    fn from(error: EquilibriumExpansionErrorV1) -> Self {
        EquilibriumBlowdownError::Thermochemistry(error).into()
    }
}
impl From<EquilibriumNozzleError> for Error {
    fn from(error: EquilibriumNozzleError) -> Self {
        if matches!(error, EquilibriumNozzleError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Nozzle(error)
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// Advance finite equilibrium gas through a pressure/area-driven restriction.
///
/// Initial, midpoint and endpoint nozzle regimes must all pass the original
/// reacting flow law. Every finite energy debit comes from the existing
/// reservoir's direct source integral and U=m*h-p*V balance, not dm*h_mid and
/// not an extra kinetic-energy sink. Elements per unit mass retain the original
/// well-mixed feed ratios; molecular composition re-equilibrates at every trial.
///
/// Smooth time accuracy is second order. The outgoing-fraction limit and actual
/// mass quadrature residual are not temporal error estimates. A fixed ambient
/// boundary cannot absorb an unreported receiver inventory. Source, work,
/// physical endpoint and final cancellation failures leave the old tank intact.
pub fn equilibrium_nozzle_blowdown(
    reservoir: &EquilibriumGasReservoir,
    back_pressure: Pressure,
    effective_area: Area,
    duration_s: f64,
    controls: EquilibriumNozzleBlowdownControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumNozzleBlowdownCandidate, Error> {
    poll(&mut cancelled)?;
    let initial = reservoir.isentropic_sample();
    let back = back_pressure.value();
    let area = effective_area.value();
    validate_flow_inputs(initial.pressure.value(), back, area, controls.nozzle)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !controls.maximum_outgoing_fraction.is_finite()
        || controls.maximum_outgoing_fraction <= 0.0
        || controls.maximum_outgoing_fraction > 0.5
        || !(1..=128).contains(&controls.maximum_mass_iterations)
        || controls.nozzle.isentropic != reservoir.thermochemical_controls()
    {
        return Err(Error::Invalid(
            "positive time, bounded withdrawal/work and unchanged birth thermochemistry required",
        ));
    }
    let mut prepared = reservoir.prepare_depletion(controls.reservoir, &mut cancelled)?;
    let first = resolve_on_isentrope(
        prepared.isentrope(),
        initial,
        reservoir.caloric_anchor(),
        back,
        area,
        controls.nozzle,
        &mut cancelled,
    )?;
    let initial_choked = first.bracket.is_some();
    let mut nozzle_queries = first.work.pressure;
    let mut density_queries = 0;

    if first.mass_flow_kg_s == 0.0 {
        let midpoint_nozzle = publish_on_path(
            prepared.isentrope(),
            reservoir.thermodynamics().clone(),
            initial,
            first,
            back,
            area,
            &mut cancelled,
        )?;
        let candidate =
            reservoir.exhaust_mass_adiabatic(0.0, controls.reservoir, &mut cancelled)?;
        let (endpoint, gas) = candidate.into_parts();
        poll(&mut cancelled)?;
        return Ok(EquilibriumNozzleBlowdownCandidate {
            endpoint,
            step: EquilibriumNozzleBlowdownStep {
                duration_s,
                gas,
                midpoint_nozzle,
                initial_choked: false,
                endpoint_choked: false,
                mass_quadrature_residual_kg: 0.0,
                mass_iterations: 0,
                density_pressure_queries: 0,
                nozzle_pressure_queries: nozzle_queries,
                isentrope_preparations: 1,
            },
        });
    }

    let quadrature = midpoint::solve(
        reservoir.mass_kg(),
        first.mass_flow_kg_s,
        OutflowMidpointSolve {
            duration_s,
            maximum_outgoing_fraction: controls.maximum_outgoing_fraction,
            mass_tolerance_kg: controls.reservoir.mass_tolerance_kg,
            maximum_iterations: controls.maximum_mass_iterations,
        },
        |mass, mut cancelled| -> Result<OutflowRate, Error> {
            let chamber = prepared.sample_remaining_mass(mass, &mut cancelled)?;
            density_queries += chamber.pressure_queries;
            // A non-driving trial supplies only the lower outward-flux bound
            // for root bracketing. No endpoint beyond ambient can be accepted.
            if chamber.sample.pressure.value() <= back {
                return Ok(OutflowRate {
                    mass_flow_kg_s: 0.0,
                    derivative_per_s: None,
                });
            }
            let anchor = prepared.isentrope().caloric_anchor(&mut cancelled)?;
            let flow = resolve_on_isentrope(
                prepared.isentrope(),
                chamber.sample,
                &anchor,
                back,
                area,
                controls.nozzle,
                &mut cancelled,
            )?;
            nozzle_queries += flow.work.pressure;
            let denominator = if flow.bracket.is_some() {
                flow.accepted.sample.response.sound_speed_m_s
            } else {
                flow.velocity_m_s
            };
            // Along the birth isentrope dh_ch/dm=c_ch^2/m. At fixed subsonic
            // back pressure rho_exit is fixed, so d(mdot)/dh_ch=mdot/v_exit^2.
            // At an ordinary sonic maximum it is mdot/c_throat^2.
            let ratio = chamber.sample.response.sound_speed_m_s / denominator;
            let derivative = (flow.mass_flow_kg_s / mass) * ratio * ratio;
            Ok(OutflowRate {
                mass_flow_kg_s: flow.mass_flow_kg_s,
                derivative_per_s: derivative.is_finite().then_some(derivative),
            })
        },
        &mut cancelled,
    )
    .map_err(|error| match error {
        OutflowMidpointError::Invalid(reason) => Error::Invalid(reason),
        OutflowMidpointError::Source(source) => source,
        OutflowMidpointError::InventoryLimit { maximum_duration_s } => {
            Error::InventoryLimit { maximum_duration_s }
        }
        OutflowMidpointError::NoConvergence {
            iterations,
            mass_residual_kg,
        } => Error::NoConvergence {
            mass_iterations: iterations,
            mass_residual_kg,
        },
        OutflowMidpointError::Cancelled => Error::Cancelled,
    })?;

    let midpoint = prepared.sample_remaining_mass(quadrature.midpoint_mass_kg, &mut cancelled)?;
    density_queries += midpoint.pressure_queries;
    if midpoint.sample.pressure.value() <= back {
        return Err(Error::Invalid(
            "positive withdrawal has no outward midpoint pressure drive",
        ));
    }
    let anchor = prepared.isentrope().caloric_anchor(&mut cancelled)?;
    let inlet = prepared.isentrope().publish_current(&mut cancelled)?;
    let resolved = resolve_on_isentrope(
        prepared.isentrope(),
        midpoint.sample,
        &anchor,
        back,
        area,
        controls.nozzle,
        &mut cancelled,
    )?;
    nozzle_queries += resolved.work.pressure;
    let midpoint_nozzle = publish_on_path(
        prepared.isentrope(),
        inlet,
        midpoint.sample,
        resolved,
        back,
        area,
        &mut cancelled,
    )?;
    let residual = quadrature.expelled_mass_kg - duration_s * midpoint_nozzle.mass_flow.value();
    if !residual.is_finite() || residual.abs() > controls.reservoir.mass_tolerance_kg {
        return Err(Error::NoConvergence {
            mass_iterations: quadrature.iterations,
            mass_residual_kg: residual,
        });
    }
    let candidate = prepared.finish_remaining_mass(quadrature.remaining_mass_kg, &mut cancelled)?;
    density_queries += candidate.step().pressure_queries;
    let endpoint = candidate.endpoint();
    let endpoint_sample = endpoint.isentropic_sample();
    if endpoint_sample.pressure.value() < back {
        return Err(Error::BackPressureCrossing {
            endpoint_pressure_pa: endpoint_sample.pressure.value(),
            back_pressure_pa: back,
        });
    }
    let final_flow = resolve_on_isentrope(
        prepared.isentrope(),
        endpoint_sample,
        endpoint.caloric_anchor(),
        back,
        area,
        controls.nozzle,
        &mut cancelled,
    )?;
    nozzle_queries += final_flow.work.pressure;
    let endpoint_choked = final_flow.bracket.is_some();
    let (endpoint, gas) = candidate.into_parts();
    poll(&mut cancelled)?;
    Ok(EquilibriumNozzleBlowdownCandidate {
        endpoint,
        step: EquilibriumNozzleBlowdownStep {
            duration_s,
            gas,
            midpoint_nozzle,
            initial_choked,
            endpoint_choked,
            mass_quadrature_residual_kg: residual,
            mass_iterations: quadrature.iterations,
            density_pressure_queries: density_queries,
            nozzle_pressure_queries: nozzle_queries,
            isentrope_preparations: 1,
        },
    })
}
