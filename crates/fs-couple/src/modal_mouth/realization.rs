//! One geometric modal mouth's complete rational radiation memory.
//!
//! All modes share the same physical reference plane. There is no separately
//! inferred inter-aperture propagation delay; geometric diffraction and modal
//! mutual memory are approximated together by the native shared poles. The
//! original BEM work matrix is already in Pa s/m³: conjugation from exp(-iwt)
//! to the fitter's s=+iw axis is the only coordinate conversion.
//!
//! The shared physical-source fitter owns the three-order ladder, independent
//! training/withheld/causality sweeps, bounded and reported residue repair,
//! continuous passivity and actual interval-mean ZOH passivity. These checks
//! do not certify the continuum mesh, unsampled frequencies or closed-loop
//! resonance/decay errors. In particular, GridOnly remains GridOnly.

use super::{
    ApertureRoute, BemModalMouth, BemModalMouthBake, CircularMouthPatch, ModalMouthBemBudget,
    ModalMouthBemError,
};
use crate::bernoulli_aperture::radiation::environment::multiport::{
    MimoRadiationLimits,
    bake::{
        BaffledMimoBakeError, BaffledMimoBakeOptions, BaffledMimoBakeReport, admit_common,
        fit_source,
    },
};
use fs_duct::modal::ModeCount;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::mimo::load::{MimoLoadError, MimoLoadLimits, PreparedMimoLoad};
use std::cell::RefCell;

/// Explicit physical domain and the original shared fitter/runtime controls.
#[derive(Clone, Copy, Debug)]
pub struct ModalMouthRealizationOptions {
    /// Physical load clock, independent of any later microphone output rate.
    pub mechanics_rate_hz: u32,
    /// Authored source support [Hz]. Every original BEM query, including the
    /// dispersion test's lower probes, must lie inside this declared band.
    /// The selected native BEM route retains all of its separate refusals.
    pub source_band_hz: [f64; 2],
    /// Rational approximation's declared band [Hz], reaching mechanics Nyquist.
    /// Fidelity below its positive lower endpoint, including DC, is not inferred.
    pub operator_band_hz: [f64; 2],
    /// Existing shared-pole order ladder, fit accuracy, memory and causality
    /// assumptions. A residue repair is explicit in the returned rung report,
    /// bounded to one tenth of its error allowance, and never preferred to an
    /// otherwise admitted unchanged fit.
    pub fit: BaffledMimoBakeOptions,
    /// Existing complete MIMO state, retained-memory and sampled-passivity limits.
    pub load: MimoRadiationLimits,
    /// Raw source transpose-defect allowance, at most one tenth of fit error.
    /// The reciprocal fit still measures error against the unchanged raw data.
    pub maximum_reciprocity_defect: f64,
}

/// Complete original physical samples and approximation declarations.
/// No source matrix or higher-mode net-flow diagnostic is rescaled.
pub struct BemModalMouthRealizationRecord {
    patch: CircularMouthPatch,
    exterior: GasState,
    modes: ModeCount,
    route: ApertureRoute,
    source_budget: ModalMouthBemBudget,
    options: ModalMouthRealizationOptions,
    training: BemModalMouthBake,
    held_out: BemModalMouthBake,
    causality: BemModalMouthBake,
    report: BaffledMimoBakeReport,
}

impl BemModalMouthRealizationRecord {
    /// Original circular section and its cap on the retained exterior mesh.
    pub fn patch(&self) -> &CircularMouthPatch {
        &self.patch
    }
    /// Original coherent exterior gas for every physical sample.
    pub const fn exterior(&self) -> &GasState {
        &self.exterior
    }
    /// Complete signed radial coordinate count, including the plane mode.
    pub const fn modes(&self) -> ModeCount {
        self.modes
    }
    /// Original selected BEM route and its numerical controls, without fallback.
    pub const fn route(&self) -> &ApertureRoute {
        &self.route
    }
    /// Original modal source binding and aggregate sampling allowances.
    pub const fn source_budget(&self) -> ModalMouthBemBudget {
        self.source_budget
    }
    /// Authored support, physical clock and unchanged numerical declarations.
    pub const fn options(&self) -> &ModalMouthRealizationOptions {
        &self.options
    }
    /// Original log-spaced training samples, including both operator endpoints.
    pub fn training(&self) -> &BemModalMouthBake {
        &self.training
    }
    /// Independently solved interleaved samples, never used in identification.
    pub fn held_out(&self) -> &BemModalMouthBake {
        &self.held_out
    }
    /// Original uniform half-step sweep used by the dispersion consistency test.
    pub fn causality(&self) -> &BemModalMouthBake {
        &self.causality
    }
    /// All three fit attempts, source-causality diagnostic and original classes
    /// of continuous and actual mean-ZOH passivity evidence.
    pub const fn report(&self) -> &BaffledMimoBakeReport {
        &self.report
    }
}

/// Quiescent native load plus its immutable source and realization record.
pub struct BemModalMouthRealization {
    load: PreparedMimoLoad,
    record: BemModalMouthRealizationRecord,
}

impl BemModalMouthRealization {
    /// Original complete MIMO load, suitable for existing continuous operator
    /// or exact interval-mean held-flow connections. Port order is radial order.
    pub fn load(&self) -> &PreparedMimoLoad {
        &self.load
    }
    /// Physical geometry, both bands, original samples and fit/passivity reports.
    pub fn record(&self) -> &BemModalMouthRealizationRecord {
        &self.record
    }
    /// Transfer the original runtime into a physical interconnection while
    /// retaining its source record; no refit, state reset or port conversion.
    pub fn into_parts(self) -> (PreparedMimoLoad, BemModalMouthRealizationRecord) {
        (self.load, self.record)
    }
}

/// A failed realization never publishes a partial runtime or guessed memory.
#[derive(Debug)]
pub enum ModalMouthRealizationError {
    /// Source support, raw reciprocity or aggregate BEM work declaration refused.
    Invalid(&'static str),
    /// Original physical geometry, source solve or full matrix admission refusal.
    Source(ModalMouthBemError),
    /// Original shared fitting, causality, accuracy or passivity refusal.
    Fit(BaffledMimoBakeError),
    /// Original complete MIMO runtime admission refusal.
    Runtime(MimoLoadError),
    /// Cancellation before publishing the complete result.
    Cancelled,
}

impl core::fmt::Display for ModalMouthRealizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "modal mouth realization: {self:?}")
    }
}
impl std::error::Error for ModalMouthRealizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source(error) => Some(error),
            Self::Fit(error) => Some(error),
            Self::Runtime(error) => Some(error),
            _ => None,
        }
    }
}
impl From<BaffledMimoBakeError> for ModalMouthRealizationError {
    fn from(error: BaffledMimoBakeError) -> Self {
        match error {
            BaffledMimoBakeError::Cancelled => Self::Cancelled,
            other => Self::Fit(other),
        }
    }
}
impl From<ModalMouthBemError> for ModalMouthRealizationError {
    fn from(error: ModalMouthBemError) -> Self {
        match error {
            ModalMouthBemError::Cancelled => Self::Cancelled,
            other => Self::Source(other),
        }
    }
}
fn invalid(what: &'static str) -> ModalMouthRealizationError {
    ModalMouthRealizationError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), ModalMouthRealizationError> {
    if gate.is_requested() {
        Err(ModalMouthRealizationError::Cancelled)
    } else {
        Ok(())
    }
}

impl BemModalMouth<'_> {
    /// Realize every original modal mouth coordinate through the common passive
    /// MIMO owner, retaining the geometry-derived signed self and mutual matrix.
    ///
    /// All three original BEM sweeps count against this binding's aggregate
    /// matrix/projection budgets before the first source solve. The source band
    /// must also cover the causality test's low-frequency probes. No table
    /// extrapolation or fabricated high-band samples fill an unsupported clock.
    /// Raw reciprocity error is bounded before reciprocal identification and
    /// remains in the final error comparison. Only the native explicitly
    /// reported, error-budgeted residue repair is available; no gain/damping fix.
    ///
    /// The prepared load retains its exact physical clock, all shared-pole
    /// states and original passivity reports. Its continuous view can join the
    /// existing linear coupled assemblies; each consumer must preserve this
    /// signed modal work pairing and own its complete interconnection audit.
    ///
    /// # Errors
    /// Incompatible declared support or budgets, any native BEM/source failure,
    /// unresolved causality/fit/passivity, runtime admission, or cancellation.
    pub fn realize(
        &self,
        route: &ApertureRoute,
        options: ModalMouthRealizationOptions,
        gate: &CancelGate,
    ) -> Result<BemModalMouthRealization, ModalMouthRealizationError> {
        poll(gate)?;
        let n = self.basis.modes().0;
        admit_common(
            &self.exterior,
            n,
            options.mechanics_rate_hz,
            options.operator_band_hz,
            options.fit,
            options.load,
        )?;
        let source = options.source_band_hz;
        let band = options.operator_band_hz;
        let count = options.fit.frequency_samples;
        let lowest = band[0].min(band[1] / (2 * count - 1) as f64);
        if source.iter().any(|f| !f.is_finite() || *f <= 0.0)
            || source[0] >= source[1]
            || source[0] > lowest
            || source[1] < band[1]
            || !options.maximum_reciprocity_defect.is_finite()
            || options.maximum_reciprocity_defect < 0.0
            || options.maximum_reciprocity_defect > 0.1 * options.fit.maximum_relative_error
        {
            return Err(invalid(
                "authored source band must cover every physical probe and reciprocity must fit its error budget",
            ));
        }
        let frequencies = (3 * count - 1) as u128;
        let entries = 2 * n as u128 * n as u128 * frequencies;
        let terms = self.surface.areas().len() as u128 * n as u128 * n as u128 * (1 + frequencies);
        if entries > self.budget.maximum_matrix_entries as u128
            || terms > u128::from(self.budget.maximum_projection_terms)
        {
            return Err(invalid(
                "complete training, withheld and causality sweeps exceed the original mouth budget",
            ));
        }
        let sampled = RefCell::new(Vec::with_capacity(3));
        let (model, report) = fit_source(
            n,
            options.mechanics_rate_hz,
            band,
            options.fit,
            options.load,
            gate,
            |hz, gate| {
                poll(gate)?;
                if hz
                    .iter()
                    .any(|f| !f.is_finite() || *f < source[0] || *f > source[1])
                {
                    return Err(invalid(
                        "original BEM query is outside its authored source band",
                    ));
                }
                let omega: Vec<_> = hz.iter().map(|f| core::f64::consts::TAU * f).collect();
                let bake = self.bake(&omega, route, gate)?;
                if bake.source().rows.iter().any(|row| {
                    !row.reciprocity_defect.is_finite()
                        || row.reciprocity_defect > options.maximum_reciprocity_defect
                }) {
                    return Err(invalid(
                        "original modal mouth reciprocity exceeds its declared fit budget",
                    ));
                }
                // These are already pressure/modal-volume-flow coordinates.
                // A higher mode's net flow can be exactly zero; never divide it.
                let rows = bake
                    .source()
                    .rows
                    .iter()
                    .map(|row| row.impedance.iter().map(|z| z.conj()).collect())
                    .collect();
                sampled.borrow_mut().push(bake);
                Ok(rows)
            },
        )?;
        poll(gate)?;
        let load = PreparedMimoLoad::new(
            &model,
            1.0 / f64::from(options.mechanics_rate_hz),
            MimoLoadLimits {
                band_rad_s: [
                    core::f64::consts::TAU * band[0],
                    core::f64::consts::TAU * band[1],
                ],
                sampled_passivity_points: options.load.sampled_passivity_points,
                max_states: options.load.max_states,
                max_memory_bytes: options.load.max_memory_bytes,
            },
        )
        .map_err(ModalMouthRealizationError::Runtime)?;
        poll(gate)?;
        let [training, held_out, causality]: [BemModalMouthBake; 3] = sampled
            .into_inner()
            .try_into()
            .map_err(|_| invalid("three complete original BEM source sweeps are required"))?;
        Ok(BemModalMouthRealization {
            load,
            record: BemModalMouthRealizationRecord {
                patch: self.patch.clone(),
                exterior: self.exterior,
                modes: self.basis.modes(),
                route: route.clone(),
                source_budget: self.budget,
                options,
                training,
                held_out,
                causality,
                report,
            },
        })
    }
}

#[cfg(test)]
mod tests;
