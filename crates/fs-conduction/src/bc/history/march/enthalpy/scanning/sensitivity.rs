//! Source-control derivatives through the actual implicit scan endpoint.
//!
//! The existing enthalpy/radiation owner revalidates the primal and builds the
//! complete physical tangent. One transpose solve supplies every source's
//! center/width/watts pullback, radiative reservoir/emissivity controls and the
//! previous enthalpy seed for a reverse sweep on a FIXED time grid. No solver
//! iterations, adaptive decisions, event times, material selections or
//! source-support admission rules are differentiated.

/// Full fixed-grid inverse objectives with authored-schedule gradients.
pub mod objective;
mod fraction;

use super::{ScanningEnthalpyRun, ScanStep, ScanTransport, events, patches_before, poll};
use super::super::super::moving::sensitivity::{MovingHeatGradient, MovingHeatLinearization};
use crate::{ConductionError, LinearConfig};
use crate::interface::ContactResistanceGradient;
use crate::interface::temperature::enthalpy::{
    ContactEnthalpyLinearization, ContactEnthalpyRadiationLinearization,
};
use crate::transient::enthalpy::EnthalpyStepConfig;
use crate::transient::enthalpy::adjoint::{EnthalpyAdjointError, EnthalpyStepLinearization};
use crate::transient::enthalpy::heterogeneous::{HeterogeneousEnthalpyError,
    HeterogeneousEnthalpyStepGradient, HeterogeneousEnthalpyStepLinearization,
    HeterogeneousEnthalpyRadiationStepLinearization};
use fs_exec::Cx;

/// Extra retained-source and radiation-tangent caps for one linearized step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanSensitivityBudget {
    /// Maximum sources times mesh vertices. Each binding retains at most three
    /// nodal scalar arrays; this count is not a byte/allocator-overhead bound.
    pub max_source_vertices: usize,
    /// Original enthalpy adjoint's total radiation-feedback scalar allowance.
    pub max_radiation_feedback_entries: usize,
}

/// Preserve the original typed refusal, including material corners and budgets.
#[derive(Debug)]
pub enum ScanSensitivityError {
    /// Primal/history/source assembly or phase admission refused.
    History(HeterogeneousEnthalpyError),
    /// Original checked implicit-function or transpose solver refused.
    Adjoint(EnthalpyAdjointError),
    /// Moving-source derivative, scalar admission, or cancellation refused.
    Source(ConductionError),
}
impl std::fmt::Display for ScanSensitivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { Self::History(e) => std::fmt::Display::fmt(e, f),
            Self::Adjoint(e) => std::fmt::Display::fmt(e, f), Self::Source(e) => std::fmt::Display::fmt(e, f) }
    }
}
impl std::error::Error for ScanSensitivityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::History(e) => Some(e), Self::Adjoint(e) => Some(e), Self::Source(e) => Some(e) }
    }
}
impl From<HeterogeneousEnthalpyError> for ScanSensitivityError {
    fn from(e: HeterogeneousEnthalpyError) -> Self { Self::History(e) }
}
impl From<EnthalpyAdjointError> for ScanSensitivityError {
    fn from(e: EnthalpyAdjointError) -> Self { Self::Adjoint(e) }
}
impl From<ConductionError> for ScanSensitivityError {
    fn from(e: ConductionError) -> Self { Self::Source(e) }
}

/// Independent physical radiation controls at this fixed scan endpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanRadiationGradient {
    /// Radiating regions in the original patch declaration order.
    pub radiation_regions: Vec<String>,
    /// Derivative per kelvin of each black-reservoir temperature, including
    /// coefficient and explicit reference-temperature changes.
    pub reservoir_temperatures: Vec<f64>,
    /// Derivative per unit absolute hemispherical emissivity; the source claim
    /// and query support remain fixed. At epsilon=1 only feasible directions apply.
    pub emissivities: Vec<f64>,
}

/// All sources and radiation controls share the SAME physical transpose.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanStepGradient {
    /// Original source-density, history, reference-density and true-residual data.
    /// `transport.previous_specific_enthalpy` seeds the preceding fixed-grid step.
    pub storage: HeterogeneousEnthalpyStepGradient,
    /// Independent instantaneous controls, in source declaration order.
    pub sources: Vec<MovingHeatGradient>,
    /// Original radiative parameters, or None when no radiation is attached.
    pub radiation: Option<ScanRadiationGradient>,
    /// Whole-law resistance multipliers for the attached live contact model,
    /// in its original canonical interface order. Empty without that model.
    pub contact_resistance_multipliers: Vec<ContactResistanceGradient>,
}

// Preserve each native owner's physical parameters alongside its common
// transport tangent. In particular, do not discard radiation metadata when
// binding a scan endpoint or solve the transpose again for its controls.
#[derive(Debug)]
enum ScanTangent<'m> {
    Storage(HeterogeneousEnthalpyStepLinearization<'m>),
    Radiation(HeterogeneousEnthalpyRadiationStepLinearization<'m>),
    Contact { binding: ContactEnthalpyLinearization<'m>, max_face_pairs: usize },
    ContactRadiation { binding: ContactEnthalpyRadiationLinearization<'m>, max_face_pairs: usize },
}
impl<'m> ScanTangent<'m> {
    fn transport(&self) -> &EnthalpyStepLinearization<'_> {
        match self {
            Self::Storage(binding) => binding.transport(),
            Self::Radiation(binding) => binding.transport(),
            Self::Contact { binding, .. } => binding.transport().transport(),
            Self::ContactRadiation { binding, .. } => binding.transport().transport(),
        }
    }

    fn pullback(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<(HeterogeneousEnthalpyStepGradient, Option<ScanRadiationGradient>,
            Vec<ContactResistanceGradient>), EnthalpyAdjointError>
    {
        match self {
            Self::Storage(binding) => Ok((binding.pullback(cx, seed, linear)?, None, Vec::new())),
            Self::Radiation(binding) => {
                let gradient = binding.pullback(cx, seed, linear)?;
                let radiation = ScanRadiationGradient {
                    radiation_regions: gradient.radiation_regions,
                    reservoir_temperatures: gradient.reservoir_temperatures,
                    emissivities: gradient.emissivities,
                };
                Ok((gradient.storage, Some(radiation), Vec::new()))
            }
            Self::Contact { binding, max_face_pairs } => {
                let gradient = binding.pullback_with_resistance_controls(cx, seed, linear, *max_face_pairs)?;
                Ok((gradient.storage, None, gradient.resistance_multipliers))
            }
            Self::ContactRadiation { binding, max_face_pairs } => {
                let gradient = binding.pullback_with_resistance_controls(cx, seed, linear, *max_face_pairs)?;
                let radiation = ScanRadiationGradient {
                    radiation_regions: gradient.radiation.radiation_regions,
                    reservoir_temperatures: gradient.radiation.reservoir_temperatures,
                    emissivities: gradient.radiation.emissivities,
                };
                Ok((gradient.radiation.storage, Some(radiation), gradient.resistance_multipliers))
            }
        }
    }
}

/// Primal report and immutable derivative binding for one accepted scan step.
/// Only the original mesh is borrowed; callers may retain a bounded collection
/// for a reverse sweep. The run itself never grows a trajectory tape.
#[derive(Debug)]
pub struct ScanStepLinearization<'m> {
    report: ScanStep,
    transport: ScanTangent<'m>,
    sources: Vec<MovingHeatLinearization<'m>>,
    // Retain each node's original chart derivative without retaining chart
    // borrows. A fraction-only corner must not disable temperature objectives.
    fraction_derivatives: Vec<fraction::FractionDerivative>,
    start_time_s: f64,
    end_time_s: f64,
}
impl ScanStepLinearization<'_> {
    /// Original full physical report, retaining source and radiation observations.
    #[must_use]
    pub const fn report(&self) -> &ScanStep { &self.report }
    /// Physical beginning of this fixed differentiation interval, seconds.
    #[must_use]
    pub const fn start_time_s(&self) -> f64 { self.start_time_s }
    /// Physical endpoint, held fixed by every returned derivative, seconds.
    #[must_use]
    pub const fn end_time_s(&self) -> f64 { self.end_time_s }
    /// Rechecked physical tangent, including any nonlinear radiation feedback.
    #[must_use]
    pub fn transport(&self) -> &EnthalpyStepLinearization<'_> { self.transport.transport() }

    /// Pull back an objective seed in nodal specific-enthalpy coordinates.
    /// Complete contact and radiation feedback, whole-law resistance multipliers
    /// and reservoir/emissivity controls use the same transpose as all moving
    /// sources, with no additional primal solve or frozen contact approximation.
    /// Prescribed nodal bath-temperature gradients remain in `storage.transport`;
    /// their constraint multipliers never enter physical heat-input pullbacks.
    /// Geometry, phase charts, conductivity, time grid and event side remain
    /// fixed; instantaneous boundary controls require a separate schedule chain
    /// rule for authored values.
    pub fn pullback(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<ScanStepGradient, ScanSensitivityError>
    {
        let (storage, radiation, contact_resistance_multipliers) =
            self.transport.pullback(cx, seed, linear)?;
        let mut sources = Vec::with_capacity(self.sources.len());
        for (i, binding) in self.sources.iter().enumerate() {
            poll(cx, i)?;
            sources.push(binding.pullback(cx, &storage.transport.source_density)?);
        }
        poll(cx, sources.len())?;
        Ok(ScanStepGradient { storage, sources, radiation, contact_resistance_multipliers })
    }

    /// Pull back a nodal temperature objective through the original dT/dh.
    /// Plateau interiors correctly have zero temperature sensitivity even when
    /// their enthalpy and liquid-fraction state can still change with heating.
    pub fn pullback_temperature(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<ScanStepGradient, ScanSensitivityError>
    {
        let enthalpy_seed = self.transport.transport().temperature_pullback(cx, seed)?;
        self.pullback(cx, &enthalpy_seed, linear)
    }

    /// Convert a dimensionless liquid-mass-fraction seed to specific-enthalpy
    /// coordinates using each node's ORIGINAL constitutive chart. No inversion
    /// of plateau temperature and no common/averaged latent heat is used.
    /// A nonzero seed at a fraction-only chart corner refuses even when the
    /// transport and temperature have unique derivatives there. Zero local
    /// seeds do not require a derivative of an unobserved fraction.
    pub fn liquid_mass_fraction_pullback(&self, cx: &Cx<'_>, seed: &[f64])
        -> Result<Vec<f64>, ScanSensitivityError>
    {
        fraction::pullback(cx, &self.fraction_derivatives,
            &self.report.solution().specific_enthalpy_j_kg, seed)
    }

    /// Differentiate a liquid-fraction objective through the same complete
    /// thermal/radiative transpose solve and normalized moving-source laws.
    /// Solid/liquid interiors have zero fraction derivative; a melting plateau
    /// generally does not. Returned previous-state seeds are in J/kg coordinates.
    pub fn pullback_liquid_mass_fraction(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<ScanStepGradient, ScanSensitivityError>
    {
        let enthalpy_seed = self.liquid_mass_fraction_pullback(cx, seed)?;
        self.pullback(cx, &enthalpy_seed, linear)
    }
}

impl<'h, 'm, 'c> ScanningEnthalpyRun<'h, 'm, 'c> {
    /// Advance and bind derivatives at one explicit right-continuous endpoint.
    /// The original primal solver runs once. Publication occurs only after
    /// source binding and physical residual/energy/smoothness checks all pass.
    /// A valid primal at a phase kink may therefore refuse this derivative API;
    /// ordinary advance_to remains available and unchanged.
    ///
    /// # Errors
    /// Original history/motion/solver refusal, derivative memory admission,
    /// unsupported classical derivative, or cancellation. Every failure leaves
    /// the accepted caloric state, clock and step counter unchanged.
    pub fn advance_linearized_to(&mut self, cx: &Cx<'_>, time_s: f64,
        config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        self.admit(cx, time_s)?;
        if self.motion_endpoint(cx, self.run.time_s(), time_s)? < time_s {
            return Err(EnthalpyAdjointError::InvalidInput("linearized scan endpoint exceeds its selected motion limit").into());
        }
        let binding = self.linearize_endpoint(cx, time_s, false, config, budget)?;
        self.run.publish(cx, time_s, binding.report.solution())?;
        Ok(binding)
    }

    /// Advance to the next event/motion-limited endpoint and bind that interval.
    /// Held loads and their derivatives use the identical pre-switch values.
    /// Returned gradients hold the selected endpoint and event side FIXED; they
    /// are not derivatives of which event or motion limit selected the interval.
    pub fn advance_linearized_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        self.admit(cx, target_s)?;
        let start = self.run.time_s();
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.signals(), start, target_s, max_step_s)?;
        let end = self.motion_endpoint(cx, start, end)?;
        let binding = self.linearize_endpoint(cx, end, true, config, budget)?;
        self.run.publish(cx, end, binding.report.solution())?;
        Ok(binding)
    }

    fn linearize_endpoint(&self, cx: &Cx<'_>, end: f64, left: bool,
        config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let mesh = self.run.problem.boundary.mesh();
        if self.sources.len().checked_mul(mesh.vertex_count())
            .is_none_or(|n| n > budget.max_source_vertices)
        {
            return Err(EnthalpyAdjointError::InvalidInput("scan source tangent storage exceeds its explicit source-vertex cap").into());
        }
        let mut sources = Vec::with_capacity(self.sources.len());
        for (i, source) in self.sources.iter().enumerate() {
            poll(cx, i)?;
            sources.push(if left { source.linearize_before(cx, mesh, end)? }
                else { source.linearize_at(cx, mesh, end)? });
        }
        let start = self.run.time_s();
        let report = self.solve(cx, &self.run.specific_h, start, end, left, config)?;
        // Immutable declarations and physical time are identical to the primal.
        // The existing binder independently reassembles the complete endpoint.
        let boundary = if left { events::boundary_before(cx, self.run.problem.boundary, end)? }
            else { self.run.problem.boundary.at(cx, end)? };
        self.run.storage.admit_prescribed_boundary(cx, &boundary)?;
        let (source, _) = self.source_at(cx, end, left)?;
        let problem = self.run.conduction_problem(&boundary, &source)?;
        let patches = if let Some(radiation) = &self.radiation {
            Some(if left { patches_before(cx, &radiation.patches, end)? } else {
                let mut patches = Vec::with_capacity(radiation.patches.len());
                for (i, patch) in radiation.patches.iter().enumerate() {
                    poll(cx, i)?;
                    patches.push(patch.at(end)?);
                }
                patches
            })
        } else { None };
        let transport = if let Some(contacts) = &self.contacts {
            let ScanTransport::Contact(accepted) = &report.transport else {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "scan contact tangent requires its original accepted contact report").into());
            };
            let max_face_pairs = contacts.max_face_pairs;
            let residual = contacts.physical_residual_tolerance_j();
            match patches.as_deref() {
                None => ScanTangent::Contact {
                    binding: contacts.contacts.linearize_accepted_enthalpy(cx,
                        &self.run.storage, problem, &self.run.specific_h, end - start,
                        config, residual, accepted.as_ref().clone(), max_face_pairs)?,
                    max_face_pairs,
                },
                Some(patches) => ScanTangent::ContactRadiation {
                    binding: contacts.contacts.linearize_accepted_enthalpy_with_ambient_radiation(cx,
                        &self.run.storage, problem, &self.run.specific_h, end - start,
                        config, residual, patches, accepted.as_ref().clone(),
                        budget.max_radiation_feedback_entries, max_face_pairs)?,
                    max_face_pairs,
                },
            }
        } else {
            let accepted = report.solution().clone();
            match patches.as_deref() {
                None => ScanTangent::Storage(self.run.storage.linearize_accepted(cx, problem,
                    self.run.problem.interfaces, &self.run.specific_h, end - start, config, accepted)?),
                Some(patches) => ScanTangent::Radiation(
                    self.run.storage.linearize_accepted_parameters_with_ambient_radiation(
                        cx, problem, self.run.problem.interfaces, &self.run.specific_h, end - start,
                        config, patches, accepted, budget.max_radiation_feedback_entries)?),
            }
        };
        let fraction_derivatives = fraction::bind(cx, &self.run.storage,
            &report.solution().specific_enthalpy_j_kg)?;
        poll(cx, self.sources.len())?;
        Ok(ScanStepLinearization { report, transport, sources, fraction_derivatives,
            start_time_s: start, end_time_s: end })
    }
}
