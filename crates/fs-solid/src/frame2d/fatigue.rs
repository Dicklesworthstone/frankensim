//! Strain-life observation on the original cyclic fiber-section frame solver.
//!
//! Each gauge follows one existing member/Lobatto-section/metal fiber. Accepted
//! mechanical strain is `epsilon_0 - y*kappa`; stress comes from the ORIGINAL
//! committed plastic history, never from an elastic reconstruction `E*epsilon`.
//! Static load, displacement control and inertial Newmark blocks share this
//! observer. Newton trials never enter rainflow, and a refused observation
//! publishes neither its mechanical block nor any counter. Unresolved cycles
//! survive calls, clones and changes between the native loading routes.
//!
//! Calibration is fixed-condition and rate-independent. The original frame
//! does not retain provenance on raw fiber laws: an optional `FiberMaterial`
//! is an EXPLICIT caller-supplied binding whose exact law and source conditions
//! are checked, not provenance inferred from a matching modulus. Thermal/creep
//! fibers and concrete fatigue are unsupported. Damage is an empirical Miner
//! estimate; it does not remove stiffness, deposit heat or certify service life.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fatigue::strain_life::{
    PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainLifeSummary, StrainStressSample,
};
use fs_material::fatigue::{FatigueError, FatigueSummary, RainflowLimits, SnCurve};
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::state_point::{ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{DynamicHistory, Frame2d, Rayleigh, SectionModel, SolidError, StaticStep};
use crate::fiber::chaboche::{ChabocheFiberLaw, ChabocheFiberState};
use crate::fiber::material::FiberMaterial;
use crate::fiber::{Fiber, FiberLaw, FiberState};

mod elastic;
use elastic::ElasticGaugeState;
pub use elastic::FrameElasticFatigueGauge;

/// Both original material families observed on the same accepted native path.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameFatigueSummary {
    pub fibers: Vec<StrainLifeSummary>,
    pub elastic: Vec<FatigueSummary>,
}

/// One original integration point and an explicitly declared axial calibration.
#[derive(Debug, Clone)]
pub struct FrameFiberFatigueGauge {
    pub element: usize,
    pub section: usize,
    pub fiber: usize,
    /// Caller-supplied original binding, required by a source-backed CMB curve.
    /// Its exact native law must match the existing frame fiber. An authored
    /// law or curve never acquires measured authority from this association.
    pub material_source: Option<Arc<FiberMaterial>>,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
}

/// Bounds on observation storage, work and the small axial-strain domain.
#[derive(Debug, Clone, Copy)]
pub struct FrameFatigueConfig {
    /// Includes the initial accepted native state, observed exactly once.
    pub history: RainflowLimits,
    pub maximum_gauges: usize,
    pub maximum_total_residue: usize,
    /// Positive finite bound, no larger than 0.05 for this small-strain model.
    /// Native J2/Voce/Chaboche full-tensor strain limits remain in force too.
    pub maximum_absolute_axial_strain: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FrameFatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch { gauge: usize, what: &'static str },
    Gauge { gauge: usize, error: FatigueError },
    Mechanics(SolidError),
    Cancelled,
    Finished,
}
impl core::fmt::Display for FrameFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native fiber-frame fatigue refused: {self:?}")
    }
}
impl std::error::Error for FrameFatigueError {}
impl From<SolidError> for FrameFatigueError {
    fn from(error: SolidError) -> Self {
        match error {
            SolidError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

#[derive(Debug, Clone)]
enum Counter {
    Scalar(StrainLifeMiner),
    Paired(PairedStrainLifeMiner),
}
impl Counter {
    fn new(gauge: &FrameFiberFatigueGauge, limits: RainflowLimits) -> Result<Self, FatigueError> {
        if gauge.mean_stress == StrainLifeMeanStress::Uncorrected {
            Ok(Self::Scalar(StrainLifeMiner::new(
                gauge.curve.clone(),
                gauge.mean_stress,
                limits,
            )?))
        } else {
            Ok(Self::Paired(PairedStrainLifeMiner::new(
                gauge.curve.clone(),
                gauge.mean_stress,
                limits,
            )?))
        }
    }
    fn push(&mut self, sample: StrainStressSample, cx: &Cx<'_>) -> Result<(), FatigueError> {
        let mut cancelled = || cx.checkpoint().is_err();
        match self {
            Self::Scalar(counter) => {
                counter.push_with_cancel(sample.mechanical_strain, &mut cancelled)
            }
            Self::Paired(counter) => counter.push_with_cancel(sample, &mut cancelled),
        }
    }
    fn summary(&self) -> StrainLifeSummary {
        match self {
            Self::Scalar(counter) => counter.summary(),
            Self::Paired(counter) => counter.summary(),
        }
    }
    fn finish(&mut self, cx: &Cx<'_>) -> Result<StrainLifeSummary, FatigueError> {
        let mut cancelled = || cx.checkpoint().is_err();
        match self {
            Self::Scalar(counter) => counter.finish_with_cancel(&mut cancelled),
            Self::Paired(counter) => counter.finish_with_cancel(&mut cancelled),
        }
    }
    fn preview_finish(&self, cx: &Cx<'_>) -> Result<StrainLifeSummary, FatigueError> {
        let mut cancelled = || cx.checkpoint().is_err();
        match self {
            Self::Scalar(counter) => counter.preview_finish_with_cancel(&mut cancelled),
            Self::Paired(counter) => counter.preview_finish_with_cancel(&mut cancelled),
        }
    }
}

#[derive(Debug, Clone)]
struct GaugeState {
    binding: FrameFiberFatigueGauge,
    counter: Counter,
    sample: StrainStressSample,
}

/// Native mechanics and bounded empirical observation, advanced atomically.
/// Cloning preserves every original constitutive state and unresolved reversal.
/// The supplied initial state is the start of observation; prior damage is not
/// reconstructed from a preloaded frame. Its immutable access prevents missed
/// native steps or an unannounced law/source switch during observation.
#[derive(Debug, Clone)]
pub struct FrameFatiguePath {
    frame: Frame2d,
    gauges: Vec<GaugeState>,
    elastic: Vec<ElasticGaugeState>,
    config: FrameFatigueConfig,
    samples: u64,
    finished: bool,
}

impl FrameFatiguePath {
    pub fn new(
        frame: Frame2d,
        gauges: Vec<FrameFiberFatigueGauge>,
        config: FrameFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameFatigueError> {
        Self::new_with_elastic(frame, gauges, Vec::new(), config, cx)
    }

    /// Observe original metal fibers and source-bound elastic material layers
    /// together. At least one gauge across the two families is required.
    pub fn new_with_elastic(
        frame: Frame2d,
        gauges: Vec<FrameFiberFatigueGauge>,
        elastic_gauges: Vec<FrameElasticFatigueGauge>,
        config: FrameFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameFatigueError> {
        poll(cx)?;
        frame.require_untimed_materials()?;
        let count = gauges
            .len()
            .checked_add(elastic_gauges.len())
            .ok_or(FrameFatigueError::Budget)?;
        if count == 0
            || count > config.maximum_gauges
            || count
                .checked_mul(config.history.maximum_residue)
                .is_none_or(|count| count > config.maximum_total_residue)
        {
            return Err(FrameFatigueError::Budget);
        }
        if !config.maximum_absolute_axial_strain.is_finite()
            || !(0.0..=0.05).contains(&config.maximum_absolute_axial_strain)
            || config.maximum_absolute_axial_strain == 0.0
        {
            return Err(FrameFatigueError::Invalid(
                "declare a positive axial strain bound no greater than 0.05",
            ));
        }
        let mut prepared = Vec::with_capacity(gauges.len());
        for (index, binding) in gauges.into_iter().enumerate() {
            poll(cx)?;
            let (fiber, _, _) = point(&frame, &binding)?;
            bind(index, &binding, &fiber.law)?;
            let sample = sample(&frame, &binding, config.maximum_absolute_axial_strain)?;
            let mut counter = Counter::new(&binding, config.history)
                .map_err(|error| gauge_error(index, error))?;
            counter
                .push(sample, cx)
                .map_err(|error| gauge_error(index, error))?;
            prepared.push(GaugeState {
                binding,
                counter,
                sample,
            });
        }
        let mut elastic = Vec::with_capacity(elastic_gauges.len());
        for (index, gauge) in elastic_gauges.into_iter().enumerate() {
            elastic.push(ElasticGaugeState::new(
                &frame,
                gauge,
                prepared.len() + index,
                config.history,
                config.maximum_absolute_axial_strain,
                cx,
            )?);
        }
        poll(cx)?;
        Ok(Self {
            frame,
            gauges: prepared,
            elastic,
            config,
            samples: 1,
            finished: false,
        })
    }

    pub const fn frame(&self) -> &Frame2d {
        &self.frame
    }
    pub const fn sample_count(&self) -> u64 {
        self.samples
    }
    pub const fn is_finished(&self) -> bool {
        self.finished
    }
    pub fn samples(&self) -> Vec<StrainStressSample> {
        self.gauges.iter().map(|g| g.sample).collect()
    }
    pub fn elastic_samples(&self) -> Vec<StrainStressSample> {
        self.elastic.iter().map(ElasticGaugeState::sample).collect()
    }
    /// Already counted damage; unresolved endpoint cycles remain live.
    pub fn summaries(&self) -> FrameFatigueSummary {
        FrameFatigueSummary {
            fibers: self.gauges.iter().map(|g| g.counter.summary()).collect(),
            elastic: self
                .elastic
                .iter()
                .map(ElasticGaugeState::summary)
                .collect(),
        }
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.gauges.get(gauge).map(|g| g.binding.curve.as_ref())
    }
    pub fn material_source(&self, gauge: usize) -> Option<&FiberMaterial> {
        self.gauges.get(gauge)?.binding.material_source.as_deref()
    }
    pub fn elastic_curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.elastic.get(gauge).map(ElasticGaugeState::curve)
    }
    pub fn elastic_material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.elastic
            .get(gauge)
            .map(ElasticGaugeState::material_source)
    }

    /// Native load ramp: observes EVERY accepted `lambda = step/steps`.
    /// As in `Frame2d::static_load`, each call ramps the supplied load from zero.
    /// Use `steps=1` for each explicitly sampled variable-amplitude load target.
    pub fn static_load(
        &mut self,
        load: &[f64],
        steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Vec<StaticStep>, FrameFatigueError> {
        self.advance(steps, cx, |frame, accepted| {
            frame.static_load_observed(load, steps, &mut || checkpoint(cx), accepted)
        })
    }

    /// Observe each original converged displacement target, including reversals.
    pub fn displacement_control(
        &mut self,
        pattern: &[f64],
        control: usize,
        targets: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<StaticStep>, FrameFatigueError> {
        self.advance(targets.len(), cx, |frame, accepted| {
            frame.displacement_control_observed(
                pattern,
                control,
                targets,
                &mut || checkpoint(cx),
                accepted,
            )
        })
    }

    /// One original Newmark block, preparing the native operator once. The
    /// first ground sample is the existing endpoint; only subsequent accepted
    /// motion steps enter rainflow. Adjacent blocks must repeat that endpoint
    /// load as required by the native solver; observation never counts it twice.
    /// Native time/error and damping semantics are unchanged. Sampling must
    /// resolve actual stress/strain extrema; endpoint observation is no proof
    /// of accuracy for an underresolved oscillation.
    #[allow(clippy::too_many_arguments)]
    pub fn newmark_with_load(
        &mut self,
        ground: &[f64],
        dt: f64,
        dir: (f64, f64),
        damping: Rayleigh,
        sustained_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicHistory, FrameFatigueError> {
        self.advance(ground.len().saturating_sub(1), cx, |frame, accepted| {
            frame.newmark_with_load_observed(
                ground,
                dt,
                dir,
                damping,
                sustained_load,
                &mut || checkpoint(cx),
                accepted,
            )
        })
    }

    /// Preview end-of-record half cycles without advancing or closing the run.
    /// These hypothetical finalized summaries must not be added across windows.
    pub fn preview_finish(&self, cx: &Cx<'_>) -> Result<FrameFatigueSummary, FrameFatigueError> {
        poll(cx)?;
        let mut fibers = Vec::with_capacity(self.gauges.len());
        for (index, gauge) in self.gauges.iter().enumerate() {
            poll(cx)?;
            fibers.push(
                gauge
                    .counter
                    .preview_finish(cx)
                    .map_err(|error| gauge_error(index, error))?,
            );
        }
        let mut elastic = Vec::with_capacity(self.elastic.len());
        for gauge in &self.elastic {
            elastic.push(gauge.preview_finish(cx)?);
        }
        poll(cx)?;
        Ok(FrameFatigueSummary { fibers, elastic })
    }

    /// Close the observation once, committing remaining half cycles atomically.
    /// Repeated successful finalization returns the same summaries.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<FrameFatigueSummary, FrameFatigueError> {
        poll(cx)?;
        let mut proposed = self.gauges.clone();
        let mut proposed_elastic = self.elastic.clone();
        let fibers = finish_gauges(&mut proposed, cx)?;
        let mut elastic = Vec::with_capacity(proposed_elastic.len());
        for gauge in &mut proposed_elastic {
            elastic.push(gauge.finish(cx)?);
        }
        poll(cx)?;
        self.gauges = proposed;
        self.elastic = proposed_elastic;
        self.finished = true;
        Ok(FrameFatigueSummary { fibers, elastic })
    }

    fn advance<T>(
        &mut self,
        steps: usize,
        cx: &Cx<'_>,
        run: impl FnOnce(
            &mut Frame2d,
            &mut dyn FnMut(&Frame2d) -> Result<(), SolidError>,
        ) -> Result<T, SolidError>,
    ) -> Result<T, FrameFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(FrameFatigueError::Finished);
        }
        let additional = u64::try_from(steps).map_err(|_| FrameFatigueError::Budget)?;
        let samples = self
            .samples
            .checked_add(additional)
            .filter(|&n| n <= self.config.history.maximum_samples)
            .ok_or(FrameFatigueError::Budget)?;
        let mut frame = self.frame.clone();
        let mut gauges = self.gauges.clone();
        let mut elastic = self.elastic.clone();
        let mut failure = None;
        let mut observed = 0_usize;
        let result = run(&mut frame, &mut |accepted| {
            let result = observe(
                accepted,
                &mut gauges,
                self.config.maximum_absolute_axial_strain,
                cx,
            )
            .and_then(|()| {
                for gauge in &mut elastic {
                    gauge.observe(accepted, cx)?;
                }
                Ok(())
            });
            match result {
                Ok(()) => {
                    observed += 1;
                    Ok(())
                }
                Err(error) => {
                    failure = Some(error);
                    Err(SolidError::MaterialRefused {
                        what: "native frame fatigue observation refused".into(),
                    })
                }
            }
        });
        let result = result.map_err(|error| failure.unwrap_or_else(|| error.into()))?;
        poll(cx)?;
        if observed != steps {
            return Err(FrameFatigueError::Invalid(
                "native accepted-step count differs from observation",
            ));
        }
        self.frame = frame;
        self.gauges = gauges;
        self.elastic = elastic;
        self.samples = samples;
        Ok(result)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameFatigueError> {
    checkpoint(cx).map_err(Into::into)
}
fn gauge_error(gauge: usize, error: FatigueError) -> FrameFatigueError {
    if error == FatigueError::Cancelled {
        FrameFatigueError::Cancelled
    } else {
        FrameFatigueError::Gauge { gauge, error }
    }
}
fn finish_gauges(
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<StrainLifeSummary>, FrameFatigueError> {
    let mut summaries = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter_mut().enumerate() {
        poll(cx)?;
        summaries.push(
            gauge
                .counter
                .finish(cx)
                .map_err(|error| gauge_error(index, error))?,
        );
    }
    poll(cx)?;
    Ok(summaries)
}
fn observe(
    frame: &Frame2d,
    gauges: &mut [GaugeState],
    maximum_strain: f64,
    cx: &Cx<'_>,
) -> Result<(), FrameFatigueError> {
    for (index, gauge) in gauges.iter_mut().enumerate() {
        poll(cx)?;
        let next = sample(frame, &gauge.binding, maximum_strain)?;
        gauge
            .counter
            .push(next, cx)
            .map_err(|error| gauge_error(index, error))?;
        gauge.sample = next;
    }
    Ok(())
}

fn point<'a>(
    frame: &'a Frame2d,
    gauge: &FrameFiberFatigueGauge,
) -> Result<(&'a Fiber, &'a FiberState, f64), FrameFatigueError> {
    let invalid = || {
        FrameFatigueError::Invalid("gauge must name an original member, fiber section and fiber")
    };
    let element = frame.elements.get(gauge.element).ok_or_else(invalid)?;
    let SectionModel::Fiber(section) = element.sections.get(gauge.section).ok_or_else(invalid)?
    else {
        return Err(invalid());
    };
    let fiber = section.fibers.get(gauge.fiber).ok_or_else(invalid)?;
    let state = section.states.get(gauge.fiber).ok_or_else(invalid)?;
    let &(eps0, curvature) = element.commit.d.get(gauge.section).ok_or_else(invalid)?;
    Ok((fiber, state, eps0 - fiber.y * curvature))
}

fn chaboche_modulus<const N: usize>(law: &UniaxialChaboche<N>) -> f64 {
    law.backbone().backbone().backbone().elastic.youngs
}
fn modulus(law: &FiberLaw) -> Option<f64> {
    Some(match law {
        FiberLaw::J2(law) => law.backbone().elastic.youngs,
        FiberLaw::Voce(law) => law.backbone().backbone().elastic.youngs,
        FiberLaw::Tabulated(law) => law.backbone().backbone().elastic.youngs,
        FiberLaw::Chaboche(law) => chaboche_modulus(law),
        FiberLaw::ChabocheBranches(law) => match law {
            ChabocheFiberLaw::Two(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Three(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Four(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Five(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Six(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Seven(l) => chaboche_modulus(l),
            ChabocheFiberLaw::Eight(l) => chaboche_modulus(l),
        },
        FiberLaw::Steel(law) => law.e0,
        _ => return None,
    })
}

fn sample(
    frame: &Frame2d,
    gauge: &FrameFiberFatigueGauge,
    maximum_strain: f64,
) -> Result<StrainStressSample, FrameFatigueError> {
    let (fiber, state, strain) = point(frame, gauge)?;
    let e = modulus(&fiber.law).ok_or(FrameFatigueError::Invalid(
        "fatigue gauges require fixed-condition metal fibers; concrete and timed/thermal laws are unsupported",
    ))?;
    if !strain.is_finite() || strain.abs() > maximum_strain {
        return Err(FrameFatigueError::Invalid(
            "accepted axial strain exceeds the observation's small-strain domain",
        ));
    }
    let stress = match (&fiber.law, state) {
        (FiberLaw::J2(_), FiberState::J2(state))
        | (FiberLaw::Voce(_), FiberState::Voce(state))
        | (FiberLaw::Tabulated(_), FiberState::Tabulated(state)) => {
            e * (strain - state.plastic_strain())
        }
        (FiberLaw::Chaboche(_), FiberState::Chaboche(state)) => {
            e * (strain - state.plastic_strain())
        }
        (FiberLaw::ChabocheBranches(law), FiberState::ChabocheBranches(state)) => {
            let plastic = match (law, state) {
                (ChabocheFiberLaw::Two(_), ChabocheFiberState::Two(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Three(_), ChabocheFiberState::Three(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Four(_), ChabocheFiberState::Four(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Five(_), ChabocheFiberState::Five(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Six(_), ChabocheFiberState::Six(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Seven(_), ChabocheFiberState::Seven(s)) => s.plastic_strain(),
                (ChabocheFiberLaw::Eight(_), ChabocheFiberState::Eight(s)) => s.plastic_strain(),
                _ => {
                    return Err(FrameFatigueError::Invalid(
                        "native Chaboche law/history branch counts differ",
                    ));
                }
            };
            e * (strain - plastic)
        }
        (FiberLaw::Steel(_), FiberState::Steel(state)) if state.eps_prev == strain => {
            state.sig_prev
        }
        _ => {
            return Err(FrameFatigueError::Invalid(
                "native metal law and accepted fiber history differ",
            ));
        }
    };
    if !stress.is_finite() {
        return Err(FrameFatigueError::Invalid(
            "accepted native fiber stress is nonfinite",
        ));
    }
    Ok(StrainStressSample {
        mechanical_strain: strain,
        normal_stress_pa: stress,
    })
}

fn same_law(left: &FiberLaw, right: &FiberLaw) -> bool {
    match (left, right) {
        (FiberLaw::J2(a), FiberLaw::J2(b)) => a == b,
        (FiberLaw::Voce(a), FiberLaw::Voce(b)) => a == b,
        (FiberLaw::Tabulated(a), FiberLaw::Tabulated(b)) => a == b,
        (FiberLaw::Chaboche(a), FiberLaw::Chaboche(b)) => a == b,
        (FiberLaw::ChabocheBranches(a), FiberLaw::ChabocheBranches(b)) => a == b,
        (FiberLaw::Steel(a), FiberLaw::Steel(b)) => {
            [a.e0, a.fy, a.b, a.r0, a.a1, a.a2] == [b.e0, b.fy, b.b, b.r0, b.a1, b.a2]
        }
        _ => false,
    }
}
fn bind(
    index: usize,
    gauge: &FrameFiberFatigueGauge,
    law: &FiberLaw,
) -> Result<(), FrameFatigueError> {
    let mismatch = |what| FrameFatigueError::SourceMismatch { gauge: index, what };
    if modulus(law) != Some(gauge.curve.parameters().youngs_modulus_pa) {
        return Err(mismatch(
            "CMB modulus must equal the original metal fiber modulus",
        ));
    }
    let Some(binding) = gauge.material_source.as_ref() else {
        return if gauge.curve.source().is_none() {
            Ok(())
        } else {
            Err(mismatch(
                "source-backed CMB needs the original caller-supplied fiber material binding",
            ))
        };
    };
    if !same_law(binding.law(), law) {
        return Err(mismatch(
            "supplied material law differs from the original frame fiber law",
        ));
    }
    let Some(material) = binding.source() else {
        return if gauge.curve.source().is_none() {
            Ok(())
        } else {
            Err(mismatch(
                "authored fiber material cannot establish a CMB source binding",
            ))
        };
    };
    let scope = gauge.curve.scope();
    let axis = gauge
        .curve
        .source()
        .map_or("T", |s| s.request().temperature_axis.as_str());
    let absolute_temperature = material
        .axis_quantities()
        .get(axis)
        .is_some_and(|quantity| {
            *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                || quantity.semantic_type().is_some_and(|semantic| {
                    semantic.kind() == QuantityKind::AbsoluteTemperature
                        && semantic.form() == ValueForm::Static
                })
        });
    if scope.material != material.material().chemistry
        || scope.process != material.material().process
        || material.state_coordinate(axis) != Some(scope.temperature_k)
        || !absolute_temperature
    {
        return Err(mismatch(
            "CMB material/process/absolute-temperature quantity differs from the supplied native source",
        ));
    }
    if let Some(source) = gauge.curve.source() {
        let youngs = material
            .property(YOUNG_MODULUS_PROPERTY)
            .ok_or_else(|| mismatch("original native modulus receipt is missing or ambiguous"))?;
        if !source.matches_material_state(
            material.card_identity(),
            material.material(),
            &material.query(),
            youngs.answer(),
        ) {
            return Err(mismatch(
                "CMB card, state, operating context or modulus claim differs from the native source",
            ));
        }
    }
    Ok(())
}
