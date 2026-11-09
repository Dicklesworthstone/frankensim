//! Cold, source-bound construction of a complete interior/exterior time load.
//!
//! The original FE boundary shapes define both the velocity work columns and
//! the original BEM aperture matrix. The homogeneous pressure operator is
//! assembled from that same coherent source gas, not an unrelated coefficient
//! field. Its immutable owner supplies the field, mode order and fitted load
//! together to the existing full-field time closure.
//!
//! Native BEM reactions remain in N s/m: no area or net-flow normalization is
//! made, including for signed modes with zero volume gain. exp(-i omega t)
//! BEM samples are conjugated exactly once before the native exp(+i omega t)
//! causality and shared-pole fit. The same fit approximates all self and mutual
//! entries. No prescribed first-arrival delay is inferred from aperture centres
//! or bounding boxes; phase propagation is present in the finite-band data.
//!
//! Three increasing proper rational orders and independent interleaved BEM
//! rows are required. This is bounded finite-band realization evidence, not a
//! proof of a uniform error bound, spatial convergence, first-arrival accuracy,
//! exterior stored energy, or coupled decay accuracy. Native continuous and
//! sampled passivity reports retain their original certificate classes. A
//! GridOnly report stays GridOnly. No passivity repair or global gain scaling
//! is applied. BEM, fit and exact-ZOH refusals propagate.

use super::super::causal::{
    CausalExteriorLoad, CausalInteriorTime, CausalTimeError, CausalTimeLimits,
};
use super::{GasState, HelmholtzFe, RadiationBudget, RadiationError, bind_openings, poll};
use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::aperture::{ApertureImpedance, ApertureRoute, aperture_impedance};
use fs_bem::panel3d::SpherePanels;
use fs_exec::CancelGate;
use fs_feec::acoustics::helmholtz::AcousticFeError;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::time::exact::ExactVelocityPort;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousHelmholtzFe, VolumeVelocitySource,
};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_math::{c64::C64, det};
use fs_vfit::causality::{CausalityError, CausalityOptions, CausalityReport, check_mimo_causality};
use fs_vfit::mimo::load::{MimoLoadError, MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::mimo::{
    MimoFitReport, MimoModel, MimoPassivityReport, mimo_vector_fit,
    shared_pole_relocation_scratch_bytes,
};
use fs_vfit::vf::{FitOptions, VfError};

/// Explicit cold bake and native fitting choices. There is no inferred band.
#[derive(Clone, Debug)]
pub struct RadiationRealizationOptions {
    /// Uniform training samples at (j + 1/2) * angular_step_rad_s.
    /// The native causality checker requires at least 8*(subtractions+1).
    pub samples: usize,
    pub angular_step_rad_s: f64,
    /// Exactly three increasing orders; every rung uses the original shared
    /// denominator fitter, including its original matrix-relative report.
    pub fits: [FitOptions; 3],
    /// Original finite-band dispersion assumptions, retained with the result.
    /// In particular the declared tail bound is not inferred from finite data.
    pub causality: CausalityOptions,
    /// Maximum matrix-relative training AND independent validation error.
    /// Must lie in (0,1). The three measured maxima must not increase beyond
    /// floating-point roundoff. This does not bound unsampled frequencies.
    pub relative_tolerance: f64,
    /// Admission for retained BEM rows, conjugated fitting data and validation
    /// evaluations. Excludes BEM native work, native pole eigensolver scratch,
    /// passivity scratch, original FE forms and the final load's own allowance.
    pub maximum_sample_bytes: usize,
    /// Original shared-denominator QR scratch ceiling, checked before baking.
    pub maximum_relocation_scratch_bytes: usize,
    /// Conservative sum of dense QR arithmetic terms across all three orders.
    /// This is a work admission, not a wall-clock promise or eigensolver bound.
    pub maximum_relocation_work: u64,
    /// Conservative entry-by-entry finite-band dispersion reconstruction work.
    /// Includes both conjugate frequency halves and subtraction products.
    pub maximum_causality_work: u64,
}

impl RadiationRealizationOptions {
    /// Exact endpoints required by the native load's declared source band.
    /// new() admits finiteness, ordering and this clock's Nyquist coverage.
    pub fn band_rad_s(&self) -> [f64; 2] {
        [
            0.5 * self.angular_step_rad_s,
            (self.samples as f64 - 0.5) * self.angular_step_rad_s,
        ]
    }
}

/// An unmodified native fit report plus independent physical-row validation.
#[derive(Clone, Debug)]
pub struct RadiationFitRung {
    pub options: FitOptions,
    pub fit: MimoFitReport,
    pub validation_max_relative_error: f64,
    /// Change from the preceding fitted order at independent BEM frequencies,
    /// normalized by the original BEM matrix, never by the preceding fit.
    pub previous_validation_change: Option<f64>,
}

/// Refusal never returns an unbound or partially admitted time owner.
#[derive(Debug)]
pub enum CausalRadiationError {
    Invalid(&'static str),
    Geometry(RadiationError),
    Field(AcousticFeError),
    Causality(CausalityError),
    Fit(VfError),
    Load(MimoLoadError),
    Time(CausalTimeError),
    Accuracy {
        reason: &'static str,
        tolerance: f64,
        /// Preserve the actual native diagnostics, including an unsuccessful
        /// measured order ladder, rather than report an invented certificate.
        rungs: Vec<RadiationFitRung>,
    },
}
impl core::fmt::Display for CausalRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound causal radiation: {self:?}")
    }
}
impl std::error::Error for CausalRadiationError {}
impl From<RadiationError> for CausalRadiationError {
    fn from(error: RadiationError) -> Self {
        Self::Geometry(error)
    }
}
fn invalid(what: &'static str) -> CausalRadiationError {
    CausalRadiationError::Invalid(what)
}

/// One immutable physical field/mode definition, original exterior samples,
/// measured rational-order ladder, and zero-history native exact-ZOH load.
///
/// Construct volume sources, quadratic boundary sources and probes from field();
/// preparation binds those original owners to this same complete pressure field.
/// Separate prepared simulations have independent native load histories.
pub struct BoundCausalRadiation<'f, 'c> {
    field: HeterogeneousHelmholtzFe<'f, 'c>,
    ports: Vec<ExactVelocityPort>,
    exterior_ports: Vec<usize>,
    source: ApertureImpedance,
    validation_source: ApertureImpedance,
    options: RadiationRealizationOptions,
    causality: Vec<CausalityReport>,
    rungs: Vec<RadiationFitRung>,
    load: PreparedMimoLoad,
}

impl<'f, 'c> BoundCausalRadiation<'f, 'c> {
    /// Bind the original planar P1 openings, bake native full matrices, run
    /// causality BEFORE fitting, and admit the highest of three proper orders.
    /// Other original modes become explicit interior feeds. All original
    /// pressure DOFs remain; no FE modal truncation or alternate engine occurs.
    ///
    /// The supplied load band must equal the two actual training endpoints and
    /// already reach this clock's Nyquist frequency. Below the positive first
    /// sample, the native sampled/DC passivity check still runs, but finite-band
    /// fit accuracy is not claimed. The selected BEM route owns its mesh,
    /// resolution and work limits; neither this function nor its fit coarsens it.
    ///
    /// This constructor intentionally assembles homogeneous source-gas forms.
    /// An arbitrary heterogeneous operator's geometry.medium() is not evidence
    /// that its actual quadrature coefficients equal this exterior gas.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fe: &'f HelmholtzFe<'c>,
        ports: &[BoundaryMode],
        exterior_ports: &[usize],
        surface: &SpherePanels,
        gas: &GasState,
        quadrature_order: usize,
        spatial_budget: RadiationBudget,
        route: &ApertureRoute,
        options: RadiationRealizationOptions,
        duration_s: f64,
        load_limits: MimoLoadLimits,
        gate: &CancelGate,
    ) -> Result<Self, CausalRadiationError> {
        poll(gate)?;
        // Reject oversized fitting/sample work before any BEM solve or bank.
        let (omega, validation_omega) =
            admit(&options, exterior_ports.len(), duration_s, load_limits)?;
        let (_, apertures) = bind_openings(
            fe,
            ports,
            exterior_ports,
            surface,
            gas,
            spatial_budget,
            gate,
        )?;
        let field = HeterogeneousHelmholtzFe::new::<AcousticFeError>(
            fe,
            quadrature_order,
            &mut |_, _| Ok(fe.medium()),
            &mut || gate.is_requested(),
        )
        .map_err(CausalRadiationError::Field)?;
        let medium = Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        let source = aperture_impedance(surface, &apertures, &omega, medium, route, gate)
            .map_err(|error| CausalRadiationError::Geometry(RadiationError::Exterior(error)))?;
        let data = conjugated_rows(&source)?;
        // Zero extraction is deliberate: reference-point distance is not the
        // first arrival of an extended mode. Keep the native checker and its
        // explicitly declared omitted-tail assumptions without asserting one.
        let zero_delays = vec![0.0; apertures.len() * apertures.len()];
        let causality = check_mimo_causality(
            &omega,
            &data,
            apertures.len(),
            &zero_delays,
            &options.causality,
        )
        .map_err(CausalRadiationError::Causality)?;
        poll(gate)?;
        let validation_source =
            aperture_impedance(surface, &apertures, &validation_omega, medium, route, gate)
                .map_err(|error| CausalRadiationError::Geometry(RadiationError::Exterior(error)))?;
        let validation = conjugated_rows(&validation_source)?;
        let (model, rungs) = fit_ladder(
            &omega,
            &data,
            &validation_omega,
            &validation,
            apertures.len(),
            &options,
            gate,
        )?;
        poll(gate)?;
        let load = PreparedMimoLoad::new(&model, duration_s, load_limits)
            .map_err(CausalRadiationError::Load)?;
        poll(gate)?;
        Ok(Self {
            field,
            ports: ports
                .iter()
                .map(|port| ExactVelocityPort {
                    profile: port.traces.clone(),
                })
                .collect(),
            exterior_ports: exterior_ports.to_vec(),
            source,
            validation_source,
            options,
            causality,
            rungs,
            load,
        })
    }

    /// Original complete pressure field, assembled from the bound source gas.
    pub fn field(&self) -> &HeterogeneousHelmholtzFe<'f, 'c> {
        &self.field
    }
    pub fn ports(&self) -> &[ExactVelocityPort] {
        &self.ports
    }
    pub fn exterior_ports(&self) -> &[usize] {
        &self.exterior_ports
    }
    /// Original exp(-i omega t) BEM rows and native diagnostics, unchanged.
    pub fn source(&self) -> &ApertureImpedance {
        &self.source
    }
    pub fn validation_source(&self) -> &ApertureImpedance {
        &self.validation_source
    }
    pub fn options(&self) -> &RadiationRealizationOptions {
        &self.options
    }
    pub fn causality_reports(&self) -> &[CausalityReport] {
        &self.causality
    }
    pub fn fit_rungs(&self) -> &[RadiationFitRung] {
        &self.rungs
    }
    /// Original admitted N s/m model, limits and native reports. This immutable
    /// template starts at zero history; each prepare() owns its own clone.
    pub fn load(&self) -> &PreparedMimoLoad {
        &self.load
    }
    pub fn passivity_reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.load.reports()
    }

    /// Reuse the existing complete-field exact step and original source work.
    /// No model, gas, aperture shape, mode index map or fit can be substituted.
    /// Local resistances, release faces and original volume sources retain the
    /// exact FE owner's original source identity and overlap refusals.
    pub fn prepare(
        &self,
        boundaries: &[ResistiveBoundary],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: CausalTimeLimits,
        gate: &CancelGate,
    ) -> Result<CausalInteriorTime<'_, 'f, 'c>, CausalRadiationError> {
        self.prepare_with_boundary_sources(
            boundaries,
            &[],
            volume_sources,
            pressure_release_faces,
            limits,
            gate,
        )
    }

    /// Drive the original BEM-loaded time field with complete quadratic inward
    /// boundary motion and volume injection, without rebuilding the exterior.
    ///
    /// Prepare every source through `field()`. The native exact FE owner keeps
    /// its original cached spatial projection, including midside-only and signed
    /// zero-net-flow shapes; no corner interpolation or area normalization is
    /// performed. Held real boundary amplitudes have units m/s, and their signed
    /// work is already included once in `CausalTimeWork::boundary_drive_work_j`.
    ///
    /// Original unselected linear ports remain first in `driven_ports()` order,
    /// followed by these sources in supplied order. BEM mode indices, native
    /// force/velocity units and fitted load states remain unchanged. A supplied
    /// boundary source is a prescribed interior feed, so it cannot share any
    /// exterior opening face, even when its full integral or amplitude is zero.
    /// Pressure-release, owner-identity, resource and cancellation refusals pass
    /// through the original exact-time closure. A successful preparation starts
    /// a separate zero-history load; advancing it never changes this template.
    pub fn prepare_with_boundary_sources(
        &self,
        boundaries: &[ResistiveBoundary],
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: CausalTimeLimits,
        gate: &CancelGate,
    ) -> Result<CausalInteriorTime<'_, 'f, 'c>, CausalRadiationError> {
        poll(gate)?;
        CausalInteriorTime::new_with_boundary_sources(
            &self.field,
            self.load.time_step_s(),
            boundaries,
            &self.ports,
            boundary_sources,
            &self.exterior_ports,
            volume_sources,
            pressure_release_faces,
            CausalExteriorLoad::Rational(self.load.clone()),
            limits,
            &mut || gate.is_requested(),
        )
        .map_err(CausalRadiationError::Time)
    }
}

fn admit(
    options: &RadiationRealizationOptions,
    ports: usize,
    duration_s: f64,
    load: MimoLoadLimits,
) -> Result<(Vec<f64>, Vec<f64>), CausalRadiationError> {
    let samples = options.samples;
    let step = options.angular_step_rad_s;
    if !(1..=64).contains(&ports)
        || !(16..=4096).contains(&samples)
        || !(1..=8).contains(&options.causality.subtractions)
        || !options.causality.tail_bound_ratio.is_finite()
        || options.causality.tail_bound_ratio < 1.0
        || !options.causality.relative_data_tolerance.is_finite()
        || options.causality.relative_data_tolerance < 0.0
        || samples < 8 * (options.causality.subtractions + 1)
        || !step.is_finite()
        || step <= 0.0
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !options.relative_tolerance.is_finite()
        || !(0.0 < options.relative_tolerance && options.relative_tolerance < 1.0)
        || options.fits.windows(2).any(|p| p[0].order >= p[1].order)
        || load.max_states > 256
        || !(17..=65537).contains(&load.sampled_passivity_points)
    {
        return Err(invalid(
            "bounded samples, increasing orders, physical clock and finite tolerance required",
        ));
    }
    let band = options.band_rad_s();
    if !band.iter().all(|w| w.is_finite() && *w > 0.0)
        || band[0].to_bits() != load.band_rad_s[0].to_bits()
        || band[1].to_bits() != load.band_rad_s[1].to_bits()
        || band[1] < core::f64::consts::PI / duration_s
    {
        return Err(invalid(
            "native source band must equal the sampled endpoints and cover Nyquist",
        ));
    }
    // Two source artifacts, once-conjugated fit/validation data and concurrent
    // previous/current validation responses plus row/report metadata.
    let entries = ports as u128 * ports as u128;
    let subtraction_terms = options.causality.subtractions as u128 + 1;
    let causality_work = 32
        * entries
        * samples as u128
        * (2 * samples as u128 + 1024)
        * subtraction_terms
        * subtraction_terms;
    if causality_work > u128::from(options.maximum_causality_work) {
        return Err(invalid(
            "native dispersion reconstruction exceeds its declared work",
        ));
    }
    let sample_bytes = samples as u128 * (entries * 160 + 256) + options.fits.len() as u128 * 512;
    if sample_bytes > options.maximum_sample_bytes as u128 {
        return Err(invalid(
            "original sample and validation banks exceed their byte budget",
        ));
    }
    let responses = ports * (ports + 1) / 2;
    let mut work = 0_u128;
    for fit in &options.fits {
        if fit.order == 0
            || fit.order > samples - 2
            || fit
                .order
                .checked_mul(ports)
                .is_none_or(|n| n > load.max_states)
            || !(1..=64).contains(&fit.iterations)
            || fit.fit_e
            || shared_pole_relocation_scratch_bytes(samples, responses, fit)
                .is_none_or(|bytes| bytes > options.maximum_relocation_scratch_bytes)
        {
            return Err(invalid(
                "proper native fit exceeds its order, state or relocation scratch limit",
            ));
        }
        let columns = fit.order as u128 + 2;
        let rows = 2 * samples as u128 * responses as u128 + 1;
        work = work.saturating_add(32 * rows * columns * columns * fit.iterations as u128);
    }
    if work > u128::from(options.maximum_relocation_work) {
        return Err(invalid(
            "three original pole relocations exceed the declared QR work",
        ));
    }
    let omega = (0..samples).map(|j| (j as f64 + 0.5) * step).collect();
    let validation = (1..samples).map(|j| j as f64 * step).collect();
    Ok((omega, validation))
}

fn conjugated_rows(source: &ApertureImpedance) -> Result<Vec<Vec<C64>>, CausalRadiationError> {
    let entries = source.apertures * source.apertures;
    source
        .rows
        .iter()
        .map(|row| {
            if row.impedance.len() != entries
                || row
                    .impedance
                    .iter()
                    .any(|z| !z.re.is_finite() || !z.im.is_finite())
                || !row.passivity_margin.is_finite()
                || row.passivity_margin < -1e-12
                || !row.reciprocity_defect.is_finite()
            {
                return Err(invalid(
                    "original BEM row is nonfinite, malformed or measurably nonpassive",
                ));
            }
            Ok(row.impedance.iter().map(|z| z.conj()).collect())
        })
        .collect()
}

/// Stable Frobenius error in original matrix units. A zero or nonfinite
/// reference refuses rather than inventing a relative-error floor.
fn relative_gap(got: &[C64], reference: &[C64]) -> Result<f64, CausalRadiationError> {
    relative_difference(got, reference, reference)
}

fn relative_difference(
    left: &[C64],
    right: &[C64],
    reference: &[C64],
) -> Result<f64, CausalRadiationError> {
    if left.len() != reference.len() || right.len() != reference.len() || reference.is_empty() {
        return Err(invalid("original matrix dimensions disagree"));
    }
    let scale = reference
        .iter()
        .map(|z| z.re.abs().max(z.im.abs()))
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid(
            "a finite nonzero original matrix is required for relative error",
        ));
    }
    let (mut numerator, mut denominator) = (0.0, 0.0);
    for ((a, b), r) in left.iter().zip(right).zip(reference) {
        let (re, im) = (r.re / scale, r.im / scale);
        let (dr, di) = (a.re / scale - b.re / scale, a.im / scale - b.im / scale);
        numerator += dr * dr + di * di;
        denominator += re * re + im * im;
    }
    let error = det::sqrt(numerator / denominator);
    if !error.is_finite() {
        return Err(invalid("native fit error left the finite set"));
    }
    Ok(error)
}

#[allow(clippy::too_many_arguments)]
fn fit_ladder(
    omega: &[f64],
    data: &[Vec<C64>],
    validation_omega: &[f64],
    validation: &[Vec<C64>],
    ports: usize,
    options: &RadiationRealizationOptions,
    gate: &CancelGate,
) -> Result<(MimoModel, Vec<RadiationFitRung>), CausalRadiationError> {
    let mut rungs = Vec::with_capacity(3);
    let mut previous: Option<Vec<Vec<C64>>> = None;
    let mut selected = None;
    for fit in &options.fits {
        poll(gate)?;
        let (model, report) =
            mimo_vector_fit(omega, data, ports, fit, true).map_err(CausalRadiationError::Fit)?;
        poll(gate)?;
        let mut worst = 0.0_f64;
        let mut change = 0.0_f64;
        let mut evaluations = Vec::with_capacity(validation.len());
        for (j, (&w, reference)) in validation_omega.iter().zip(validation).enumerate() {
            poll(gate)?;
            let got = model.eval_iw(w);
            worst = worst.max(relative_gap(&got, reference)?);
            if let Some(prior) = &previous {
                // Normalize by the original BEM matrix, not either fitted model.
                change = change.max(relative_difference(&got, &prior[j], reference)?);
            }
            evaluations.push(got);
        }
        rungs.push(RadiationFitRung {
            options: *fit,
            fit: report,
            validation_max_relative_error: worst,
            previous_validation_change: previous.as_ref().map(|_| change),
        });
        previous = Some(evaluations);
        selected = Some(model);
    }
    let measured = |rung: &RadiationFitRung| {
        rung.fit
            .max_relative_error
            .max(rung.validation_max_relative_error)
    };
    let nonfinite = rungs.iter().any(|r| {
        !r.fit.max_relative_error.is_finite()
            || !r.fit.rms_relative_error.is_finite()
            || !r.fit.data_reciprocity_defect.is_finite()
    });
    let increasing = rungs.windows(2).any(|r| {
        measured(&r[1]) > measured(&r[0]) * (1.0 + 256.0 * f64::EPSILON) + 256.0 * f64::EPSILON
    });
    if nonfinite || increasing || measured(&rungs[2]) > options.relative_tolerance {
        return Err(CausalRadiationError::Accuracy {
            reason: if nonfinite {
                "native fit diagnostics are nonfinite"
            } else if increasing {
                "measured original-row error increased with rational order"
            } else {
                "highest order exceeds original training or independent BEM tolerance"
            },
            tolerance: options.relative_tolerance,
            rungs,
        });
    }
    Ok((selected.expect("three declared orders"), rungs))
}

#[cfg(test)]
mod tests;
