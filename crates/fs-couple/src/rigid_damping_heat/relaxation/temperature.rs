//! Calibrated temperature kinetics on the original spatial Maxwell graph.
//!
//! The existing material owner supplies a_T = tau(T)/tau_ref. Every original
//! Kelvin form scales by a_T and every original branch time by the SAME a_T.
//! Elastic factors, memory coordinates and storage remain unchanged. Actual
//! finite-cell temperatures are sampled on the existing thermal event grid;
//! only the native Kelvin/Maxwell dissipation is deposited as heat.
//!
//! Kinetics are frozen at the start of each actual thermal interval. This is
//! first-order thermal feedback; refine max_dt_s independently of energy
//! conservation. No vertical modulus shift, thermal expansion, entropy model,
//! equilibrium spring, moisture calibration or source validity is inferred.

use super::super::{
    DampingSource, Result, RigidDampingThermalError, RigidDampingThermalReport, ThermalControls,
    ThermalDriver, advance_with_source_driver, require_implicit_report,
};
use super::{
    HeldRelaxationKinetics, RelaxationDampingImplicitConfig, RelaxationDampingThermalCandidate,
    RelaxationDampingThermalConfig, RelaxationDampingThermalInput, RelaxationMotion,
    RelaxationNetworkMotionControls, SpatialRelaxationBranch, prepare_links,
};
use fs_material::visco::{ViscoError, temperature::TemperatureShift};
use fs_mbd::{
    RigidBodyState,
    spatial_loads::internal_impulse::viscous::network::{
        ViscousNetworkLink, relaxation::RelaxationNetworkLink,
    },
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionObserver, DepositionSimulation};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionReport, ImplicitConductionWorkspace,
};
use std::{borrow::Cow, collections::BTreeMap};

/// Calibration or an attempted reinterpretation of original source data.
#[derive(Clone, Debug, PartialEq)]
pub enum RelaxationTemperatureError {
    Invalid(&'static str),
    Reference { connector: usize },
    Shift(ViscoError),
}
impl core::fmt::Display for RelaxationTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "relaxation temperature calibration: {message}"),
            Self::Reference { connector } => {
                write!(
                    f,
                    "connector {connector} changed its original relaxation source"
                )
            }
            Self::Shift(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for RelaxationTemperatureError {}

/// One original spatial source and explicit material calibration per connector.
///
/// All branch times and the instantaneous Kelvin form are declared at that
/// connector's reference temperature. Retain this binding and those original
/// input laws across continuation calls. Factors, branch order and reference
/// times are checked against the binding before any candidate advances.
/// Different connectors may have different references and calibrated support.
///
/// This is thermorheological simplicity within each connector: its complete
/// relaxation spectrum shares one shift, while its elastic storage stays fixed.
/// Calibration is authored source data; no material-name constants are supplied.
#[derive(Clone, Copy, Debug)]
pub struct CalibratedRelaxationDamping<'a> {
    reference: &'a [ViscousNetworkLink],
    branches: &'a [Vec<SpatialRelaxationBranch>],
    shifts: &'a [TemperatureShift],
}
impl<'a> CalibratedRelaxationDamping<'a> {
    pub fn try_new(
        reference: &'a [ViscousNetworkLink],
        branches: &'a [Vec<SpatialRelaxationBranch>],
        shifts: &'a [TemperatureShift],
    ) -> core::result::Result<Self, RelaxationTemperatureError> {
        if reference.len() != branches.len() || reference.len() != shifts.len() {
            return Err(RelaxationTemperatureError::Invalid(
                "one complete original branch bank and temperature calibration per connector required",
            ));
        }
        Ok(Self {
            reference,
            branches,
            shifts,
        })
    }

    pub fn reference_links(&self) -> &[ViscousNetworkLink] {
        self.reference
    }
    pub fn reference_branches(&self) -> &[Vec<SpatialRelaxationBranch>] {
        self.branches
    }
    pub fn shifts(&self) -> &[TemperatureShift] {
        self.shifts
    }
}

/// Jointly accepted native bodies, elastic histories, finite heat and kinetics.
#[must_use]
pub struct TemperatureRelaxationDampingThermalCandidate<'m, O> {
    native: RelaxationDampingThermalCandidate<'m, O>,
    temperatures_k: Vec<f64>,
    factors: Vec<f64>,
    evaluations: usize,
}
impl<'m, O> TemperatureRelaxationDampingThermalCandidate<'m, O> {
    pub fn states(&self) -> &[RigidBodyState] {
        self.native.states()
    }
    pub fn memory(&self) -> &[Vec<[f64; 6]>] {
        self.native.memory()
    }
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        self.native.thermal()
    }
    pub fn observer(&self) -> &O {
        self.native.observer()
    }
    pub fn report(&self) -> &RigidDampingThermalReport {
        self.native.report()
    }
    /// Actual endpoint caloric temperature for each original connector [K].
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    /// Endpoint tau(T)/tau_ref and C(T)/C_ref in original connector order.
    /// These are outputs; continuation keeps the original calibration binding.
    pub fn shift_factors(&self) -> &[f64] {
        &self.factors
    }
    /// Successful material evaluations, including the final endpoint refresh.
    /// Unchanged temperatures reuse the previous admitted coefficients.
    pub const fn law_evaluations(&self) -> usize {
        self.evaluations
    }
    /// Move accepted native mechanics/history plus endpoint temperature/shift.
    pub fn into_parts(self) -> (RelaxationDampingThermalCandidate<'m, O>, Vec<f64>, Vec<f64>) {
        (self.native, self.temperatures_k, self.factors)
    }
}

/// Advance original material history and finite heat with calibrated feedback.
///
/// One explicit CellId per connector supplies its temperature independently of
/// heat-allocation routes. Shared temperature cells are sampled once per source
/// synchronization. Same-time births complete before sampling, and all final
/// native events complete before endpoint coefficients are admitted.
///
/// Native held kinetics are always evaluated from the original laws. A rate
/// change leaves the original energy-amplitude basis and stored energy intact;
/// no memory reset, factor rebinding or extra parameter-work/heat term occurs.
/// Full Kelvin and Maxwell losses follow the same native coupled trajectory.
///
/// A source-domain, effective-coefficient, caloric, observer, budget or energy
/// refusal discards the entire candidate. Original inputs remain borrowed.
/// Additional caches are linear in the admitted connector/branch counts.
pub fn advance_candidate<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    calibration: &CalibratedRelaxationDamping<'_>,
    temperature_cells: &[CellId],
    target_s: f64,
    config: RelaxationDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<TemperatureRelaxationDampingThermalCandidate<'m, O>, O::Error, RelaxationTemperatureError>
{
    advance_with_driver(
        input,
        thermal,
        observer,
        calibration,
        temperature_cells,
        target_s,
        config.mechanical,
        config.into(),
        ThermalDriver::Explicit,
        cancelled,
    )
    .map(|(candidate, _)| candidate)
}

/// Advance calibrated Maxwell history with native implicit thermal conduction.
///
/// The original Kelvin forms and relaxation times receive the same calibrated
/// shift at each accepted physical interval start. Both SDIRK stages hold the
/// resulting whole-interval heat proposal. Rejected thermal intervals retain
/// accepted bodies, memory and source coefficients; the proposed mechanical
/// interval is recomputed without resetting physical history or work counts.
/// Final caloric temperatures and source-domain admission are synchronized
/// after all native events before the joint candidate is returned.
///
/// BE/SDIRK conduction does not change the first-order lagged temperature
/// feedback. Refine maximum_dt_s independently for kinetics and moving geometry.
/// The original branch factors and memory basis stay fixed, so no elastic
/// parameter work or extra heat is introduced. Graph/source/event admission is
/// the same as [`super::super::advance_implicit_candidate`].
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    calibration: &CalibratedRelaxationDamping<'_>,
    temperature_cells: &[CellId],
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: RelaxationDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (TemperatureRelaxationDampingThermalCandidate<'m, O>, ImplicitConductionReport),
    O::Error,
    RelaxationTemperatureError,
> {
    let (candidate, report) = advance_with_driver(
        input,
        thermal,
        observer,
        calibration,
        temperature_cells,
        target_s,
        config.mechanical,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

#[allow(clippy::too_many_arguments)]
fn advance_with_driver<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    calibration: &CalibratedRelaxationDamping<'_>,
    temperature_cells: &[CellId],
    target_s: f64,
    mechanical: RelaxationNetworkMotionControls,
    config: ThermalControls,
    driver: ThermalDriver<'_, '_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (TemperatureRelaxationDampingThermalCandidate<'m, O>, Option<ImplicitConductionReport>),
    O::Error,
    RelaxationTemperatureError,
> {
    let links = prepare_links(input, mechanical, &mut cancelled)?;
    let n = links.len();
    if temperature_cells.len() != n
        || calibration.reference.len() != n
        || calibration.branches.len() != n
    {
        return Err(RigidDampingThermalError::Invalid(
            "one explicit caloric coordinate and original calibration per connector required",
        ));
    }
    for (connector, link) in links.iter().enumerate() {
        if cancelled() {
            return Err(RigidDampingThermalError::Cancelled);
        }
        if link.viscous != calibration.reference[connector]
            || link.branches != calibration.branches[connector]
        {
            return Err(RigidDampingThermalError::Source {
                connector,
                error: RelaxationTemperatureError::Reference { connector },
            });
        }
    }
    let source = TemperatureRelaxationSource {
        reference: input.mechanical.links,
        links: &links,
        shifts: calibration.shifts,
        cells: temperature_cells,
        temperatures_k: vec![f64::NAN; n],
        factors: vec![1.0; n],
        time_scales: links
            .iter()
            .map(|link| vec![1.0; link.branches.len()])
            .collect(),
        evaluations: 0,
    };
    let mechanics = RelaxationMotion {
        links: &links,
        initial_memory: input.memory,
        memory: Cow::Borrowed(input.memory),
        efforts: input.efforts,
        controls: mechanical,
    };
    let (native, source, mechanics, report) = advance_with_source_driver(
        input.mechanical,
        thermal,
        observer,
        source,
        mechanics,
        target_s,
        config,
        driver,
        cancelled,
    )?;
    Ok((
        TemperatureRelaxationDampingThermalCandidate {
            native: RelaxationDampingThermalCandidate {
                native,
                memory: mechanics.memory.into_owned(),
            },
            temperatures_k: source.temperatures_k,
            factors: source.factors,
            evaluations: source.evaluations,
        },
        report,
    ))
}

#[derive(Clone)]
struct TemperatureRelaxationSource<'a> {
    reference: &'a [ViscousNetworkLink],
    links: &'a [RelaxationNetworkLink],
    shifts: &'a [TemperatureShift],
    cells: &'a [CellId],
    temperatures_k: Vec<f64>,
    factors: Vec<f64>,
    time_scales: Vec<Vec<f64>>,
    evaluations: usize,
}
impl TemperatureRelaxationSource<'_> {
    fn kinetics(&self) -> HeldRelaxationKinetics<'_> {
        HeldRelaxationKinetics {
            viscosity_scales: &self.factors,
            relaxation_time_scales: &self.time_scales,
        }
    }
}
impl DampingSource for TemperatureRelaxationSource<'_> {
    type Error = RelaxationTemperatureError;
    fn links(&self) -> &[ViscousNetworkLink] {
        self.reference
    }
    fn relaxation_kinetics(&self) -> Option<HeldRelaxationKinetics<'_>> {
        Some(self.kinetics())
    }
    fn synchronize<E>(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        requested: &dyn Fn() -> bool,
    ) -> Result<(), E, Self::Error> {
        let mut sampled = BTreeMap::new();
        let mut changed = false;
        for (connector, &cell) in self.cells.iter().enumerate() {
            if requested() {
                return Err(RigidDampingThermalError::Cancelled);
            }
            let temperature = match sampled.entry(cell) {
                std::collections::btree_map::Entry::Occupied(entry) => *entry.get(),
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let temperature = thermal.temperature_k(cell)?;
                    if !temperature.is_finite() || temperature <= 0.0 {
                        return Err(RigidDampingThermalError::Invalid(
                            "finite positive original caloric temperature required",
                        ));
                    }
                    *entry.insert(temperature)
                }
            };
            if self.temperatures_k[connector] == temperature {
                continue;
            }
            let state = self.shifts[connector]
                .evaluate(temperature)
                .map_err(|error| RigidDampingThermalError::Source {
                    connector,
                    error: RelaxationTemperatureError::Shift(error),
                })?;
            self.factors[connector] = state.factor();
            self.time_scales[connector].fill(state.factor());
            self.temperatures_k[connector] = temperature;
            self.evaluations = self
                .evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget(
                    "temperature-dependent relaxation evaluations",
                ))?;
            changed = true;
        }
        if changed {
            // Includes final endpoint refresh: no dense preparation is required
            // to refuse an unrepresentable original Kelvin form or branch rate.
            self.kinetics().validate(self.links)?;
        }
        Ok(())
    }
}
