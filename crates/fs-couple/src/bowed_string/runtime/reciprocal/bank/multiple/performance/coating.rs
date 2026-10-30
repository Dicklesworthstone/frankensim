//! Local film-strain applicability on the original accepted performance state.
//!
//! This observer uses every retained body mode, every coated triangle and both
//! film faces through the original MITC owner. It neither changes the mechanics
//! nor rescales pressure. Checks cover the initial state and each mechanical
//! endpoint, including thermal playback, before pressure reaches resampling.
//! They do not bound extrema between endpoints, omitted modes, cure strain or
//! mean strain absent from the retained coordinates.

use super::{MultiBowPerformance, MultiBowStringBank, render_invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::shear_shell::coating::strain::{
    CoatingStrainSource, OwnedShellCoatingStrain, ShellCoatingFilmSource, ShellCoatingStrainBudget,
    ShellCoatingStrainResources,
};
use fs_exec::CancelGate;
use std::{mem::size_of, sync::Arc};

/// Additional local-observation allowances; original mechanics limits stand.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShellCoatingMonitorOptions {
    /// Cached observer, per-film running maxima and reusable real modal buffer.
    /// Borrowed original images, allocator overhead and original core stack
    /// frames are excluded, as in [`ShellCoatingStrainBudget`].
    pub maximum_bytes: usize,
    /// Preparation plus the initial state and every scheduled endpoint.
    pub maximum_work: usize,
}

/// The original observation's refusal, retained after an output failure.
#[derive(Debug)]
pub enum ShellCoatingMonitorError {
    /// Invalid source, clock, extent or resource declaration.
    Invalid(&'static str),
    /// The original strain owner refused at this local mechanical endpoint.
    Observation {
        sample: u64,
        source: AcousticRealizeError,
    },
    /// A sufficient spatial bound cannot establish the authored envelope.
    /// This is not evidence that the actual film has failed or been damaged.
    Envelope {
        sample: u64,
        component_index: Option<usize>,
        film_index: usize,
        bound: f64,
        limit: f64,
        ratio: f64,
    },
}
impl core::fmt::Display for ShellCoatingMonitorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shell coating strain monitor: {self:?}")
    }
}
impl std::error::Error for ShellCoatingMonitorError {}
fn invalid(what: &'static str) -> ShellCoatingMonitorError {
    ShellCoatingMonitorError::Invalid(what)
}

/// Per-film maxima across all completely checked endpoints, including zero.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellCoatingFilmMaximum {
    pub covered_elements: usize,
    pub frobenius_upper_bound: f64,
    /// Separately maximized absolute tensor components [xx, yy, zz, xy, yz, zx].
    /// Their maxima need not share a sample, location or phase.
    pub component_abs_upper_bounds: [f64; 6],
    pub limit_ratio: f64,
    /// Earliest local endpoint attaining the largest Frobenius bound.
    pub sample: u64,
    /// These identify a Bernstein coefficient, not a physical maximum site.
    pub ruling_element: usize,
    pub ruling_control_point: usize,
    pub ruling_z_face: usize,
}

/// A numeric maximum together with its exact immutable material/region source.
#[derive(Clone, Copy, Debug)]
pub struct ShellCoatingMonitorFilm<'a> {
    pub source: ShellCoatingFilmSource<'a>,
    pub maximum: &'a ShellCoatingFilmMaximum,
}

/// Borrowed running observation, available even after a refused output callback.
pub struct ShellCoatingMonitorReport<'a> {
    pub options: ShellCoatingMonitorOptions,
    /// Includes both the original observer and this adapter's retained buffers
    /// and complete planned query work. No successful query is omitted to fit.
    pub resources: &'a ShellCoatingStrainResources,
    pub initial_absolute_sample: u64,
    /// Includes the initial state. A failed numerical query is not counted.
    pub checked_endpoints: usize,
    /// Local endpoint index; zero is the state before this performance starts.
    pub last_checked_sample: u64,
    pub failure: Option<&'a ShellCoatingMonitorError>,
    source: CoatingStrainSource<'a>,
    maxima: &'a [ShellCoatingFilmMaximum],
}
impl<'a> ShellCoatingMonitorReport<'a> {
    /// Original component/film order; all film regions remain represented.
    pub fn films(&self) -> impl Iterator<Item = ShellCoatingMonitorFilm<'a>> + '_ {
        self.source
            .film_sources()
            .zip(self.maxima)
            .map(|(source, maximum)| ShellCoatingMonitorFilm { source, maximum })
    }
}

pub(super) struct BoundShellCoatingMonitor {
    options: ShellCoatingMonitorOptions,
    resources: ShellCoatingStrainResources,
    observer: OwnedShellCoatingStrain,
    gate: Arc<CancelGate>,
    displacement: Vec<f64>,
    maxima: Vec<ShellCoatingFilmMaximum>,
    initial_absolute_sample: u64,
    last_checked_sample: u64,
    failure: Option<ShellCoatingMonitorError>,
}

fn source(bank: &MultiBowStringBank) -> Result<CoatingStrainSource<'_>, ShellCoatingMonitorError> {
    let source = bank.source().bank();
    match (
        source.material_shell_image(),
        source.material_shell_assembly_image(),
    ) {
        (Some(image), None) => Ok(CoatingStrainSource::Shell(image)),
        (None, Some(image)) => Ok(CoatingStrainSource::Assembly(image)),
        _ => Err(invalid(
            "coating monitoring requires one original material shell image",
        )),
    }
}

impl MultiBowPerformance {
    /// Observe local film strain before rendering and at every accepted step.
    ///
    /// The shared cancellation gate is retained so large local-field queries
    /// remain cancellable within callbacks. Binding is allowed once before this
    /// performance advances, after all preload/state-changing bindings. The original accepted q includes any retained mean
    /// displacement; no absent prestress/cure strain is invented. An envelope
    /// refusal stops output without reverting the already accepted mechanics.
    pub fn with_material_shell_coating_monitor(
        mut self,
        options: ShellCoatingMonitorOptions,
        gate: Arc<CancelGate>,
    ) -> Result<Self, ShellCoatingMonitorError> {
        self.validate_sample_count(0)
            .map_err(|_| invalid("coating monitoring requires a healthy performance"))?;
        if self.samples_rendered() != 0 || self.coating_monitor.is_some() {
            return Err(invalid(
                "bind coating monitoring once before performance playback",
            ));
        }
        let source = source(&self.bank)?;
        let modes = match source {
            CoatingStrainSource::Shell(image) => image.physical_modes().len(),
            CoatingStrainSource::Assembly(image) => image
                .physical_modes()
                .first()
                .ok_or_else(|| invalid("coating monitoring requires an assembled modal basis"))?
                .len(),
        };
        let films = source.film_sources().count();
        let evaluations = self
            .config
            .samples
            .checked_add(1)
            .and_then(|count| usize::try_from(count).ok())
            .ok_or_else(|| invalid("complete coating-monitor endpoint count overflows"))?;
        // Includes the source-independent running maxima and one reusable real
        // q buffer, with no per-sample allocation or copied material image.
        let overhead_bytes = modes
            .checked_mul(size_of::<f64>())
            .and_then(|q| {
                films
                    .checked_mul(size_of::<ShellCoatingFilmMaximum>())
                    .and_then(|peaks| q.checked_add(peaks))
            })
            .and_then(|n| n.checked_add(size_of::<BoundShellCoatingMonitor>()))
            .ok_or_else(|| invalid("coating-monitor storage overflows"))?;
        let per_evaluation_work = films
            .checked_mul(64)
            .and_then(|n| modes.checked_mul(4).and_then(|q| n.checked_add(q)))
            .and_then(|n| n.checked_add(32))
            .ok_or_else(|| invalid("coating-monitor endpoint work overflows"))?;
        let overhead_work = per_evaluation_work
            .checked_mul(evaluations)
            .ok_or_else(|| invalid("complete coating-monitor work overflows"))?;
        let maximum_bytes = options
            .maximum_bytes
            .checked_sub(overhead_bytes)
            .ok_or_else(|| invalid("coating-monitor state exceeds its byte allowance"))?;
        let maximum_work = options
            .maximum_work
            .checked_sub(overhead_work)
            .ok_or_else(|| {
                invalid("complete coating-monitor queries exceed their work allowance")
            })?;
        let observer = source
            .prepare_snapshot(
                ShellCoatingStrainBudget {
                    maximum_bytes,
                    maximum_work,
                    evaluations,
                },
                &gate,
            )
            .map_err(|source| ShellCoatingMonitorError::Observation { sample: 0, source })?;
        let mut resources = *observer.resources();
        resources.retained_bytes = resources
            .retained_bytes
            .checked_add(overhead_bytes)
            .ok_or_else(|| invalid("complete coating-monitor storage overflows"))?;
        resources.per_evaluation_work = resources
            .per_evaluation_work
            .checked_add(per_evaluation_work)
            .ok_or_else(|| invalid("complete coating-monitor endpoint work overflows"))?;
        resources.total_work = resources
            .total_work
            .checked_add(overhead_work)
            .ok_or_else(|| invalid("complete coating-monitor work overflows"))?;
        let mut monitor = BoundShellCoatingMonitor {
            options,
            resources,
            observer,
            gate,
            displacement: vec![0.0; modes],
            maxima: vec![ShellCoatingFilmMaximum::default(); films],
            initial_absolute_sample: self.bank.system().samples_rendered(),
            last_checked_sample: 0,
            failure: None,
        };
        monitor.observe(&self.bank, 0)?;
        self.coating_monitor = Some(monitor);
        Ok(self)
    }

    /// Read running film maxima and the first refusal without advancing state.
    pub fn coating_monitor_report(&self) -> Option<ShellCoatingMonitorReport<'_>> {
        self.coating_monitor.as_ref().and_then(|monitor| {
            source(&self.bank)
                .ok()
                .map(|source| ShellCoatingMonitorReport {
                    options: monitor.options,
                    resources: &monitor.resources,
                    initial_absolute_sample: monitor.initial_absolute_sample,
                    checked_endpoints: monitor.observer.completed_evaluations(),
                    last_checked_sample: monitor.last_checked_sample,
                    failure: monitor.failure.as_ref(),
                    source,
                    maxima: &monitor.maxima,
                })
        })
    }
}

impl BoundShellCoatingMonitor {
    fn stop(&mut self, error: ShellCoatingMonitorError) -> RenderError {
        self.failure = Some(error);
        render_invalid("local coating strain observation refused; inspect coating_monitor_report")
    }

    pub(super) fn before_step(
        &mut self,
        bank: &MultiBowStringBank,
        sample: u64,
    ) -> Result<(), RenderError> {
        if self.gate.is_requested() {
            return Err(self.stop(ShellCoatingMonitorError::Observation {
                sample,
                source: AcousticRealizeError::Cancelled,
            }));
        }
        if self.last_checked_sample != sample || self.failure.is_some() {
            return Err(self.stop(invalid("coating monitor lost its accepted endpoint clock")));
        }
        // Later setup may settle a prescribed load without advancing time.
        // Such state changes must precede monitor binding so sample zero is
        // actually checked, not silently replaced before the first endpoint.
        let component = bank.source().bank().strings().len();
        let states = bank
            .system()
            .components()
            .get(component)
            .map(|body| body.states());
        if states.is_none_or(|states| {
            states.len() != self.displacement.len()
                || states
                    .iter()
                    .zip(&self.displacement)
                    .any(|(state, q)| state.displacement_m_sqrt_kg.to_bits() != q.to_bits())
        }) {
            return Err(self.stop(invalid(
                "bind coating monitoring after all preload or body-state changes",
            )));
        }
        Ok(())
    }

    pub(super) fn accepted_endpoint(
        &mut self,
        bank: &MultiBowStringBank,
        sample: u64,
    ) -> Result<(), RenderError> {
        self.observe(bank, sample).map_err(|error| self.stop(error))
    }

    fn observe(
        &mut self,
        bank: &MultiBowStringBank,
        sample: u64,
    ) -> Result<(), ShellCoatingMonitorError> {
        if self.initial_absolute_sample.checked_add(sample)
            != Some(bank.system().samples_rendered())
            || (self.observer.completed_evaluations() != 0
                && self.last_checked_sample.checked_add(1) != Some(sample))
        {
            return Err(invalid(
                "coating monitor must observe every accepted body endpoint",
            ));
        }
        let component = bank.source().bank().strings().len();
        let states = bank
            .system()
            .components()
            .get(component)
            .ok_or_else(|| invalid("coating monitor lost the original body component"))?
            .states();
        if states.len() != self.displacement.len() {
            return Err(invalid(
                "coating monitor lost the complete original body mode extent",
            ));
        }
        for (q, state) in self.displacement.iter_mut().zip(states) {
            *q = state.displacement_m_sqrt_kg;
        }
        let source = source(bank)?;
        let summary = self
            .observer
            .evaluate_snapshot(source, &self.displacement, &self.gate)
            .map_err(|source| ShellCoatingMonitorError::Observation { sample, source })?;
        let initial = self.observer.completed_evaluations() == 1;
        for (peak, value) in self.maxima.iter_mut().zip(self.observer.values()) {
            if initial || value.frobenius_upper_bound > peak.frobenius_upper_bound {
                peak.covered_elements = value.covered_elements;
                peak.frobenius_upper_bound = value.frobenius_upper_bound;
                peak.limit_ratio = value.limit_ratio;
                peak.sample = sample;
                peak.ruling_element = value.ruling_element;
                peak.ruling_control_point = value.ruling_control_point;
                peak.ruling_z_face = value.ruling_z_face;
            }
            for (maximum, component) in peak
                .component_abs_upper_bounds
                .iter_mut()
                .zip(value.component_amplitudes)
            {
                *maximum = maximum.max(component);
            }
        }
        self.last_checked_sample = sample;
        if !summary.within_authored_envelope {
            for (film, value) in source.film_sources().zip(self.observer.values()) {
                if value.limit_ratio > 1.0 {
                    return Err(ShellCoatingMonitorError::Envelope {
                        sample,
                        component_index: film.component_index,
                        film_index: film.film_index,
                        bound: value.frobenius_upper_bound,
                        limit: film.source.input().material.strain_limit(),
                        ratio: value.limit_ratio,
                    });
                }
            }
            return Err(invalid(
                "coating monitor received an inconsistent envelope status",
            ));
        }
        Ok(())
    }
}
