//! State-independent compilation of physical loads in an existing modal basis.
//!
//! The original force renderer and coupled material performers share this
//! compiler. It never constructs, clones, preloads or advances a mechanical
//! model. The caller supplies its actual conjugate displacement/load columns.

use super::{ForceRenderConfig, ModalForceEvent, charge, control, project, sizing};
use crate::render::RenderError;

/// Physical actuator columns and initial held newtons for one structural image.
/// A column maps newtons to mass-normalized modal force and modal velocity back
/// to the work-conjugate point velocity. Source/geometry authority stays with
/// the producer of these columns, not this numerical compiler.
#[derive(Clone, Debug)]
pub struct PhysicalForceInput {
    columns: Vec<Vec<f64>>,
    initial_force_n: Vec<f64>,
}
impl PhysicalForceInput {
    /// Require aligned, finite columns and physical inputs. Signed columns and
    /// exact zeros are preserved; no gain or modal-mass normalization is added.
    pub fn new(columns: Vec<Vec<f64>>, initial_force_n: Vec<f64>) -> Result<Self, RenderError> {
        let modes = columns.first().map_or(0, Vec::len);
        if modes == 0 || columns.len() != initial_force_n.len()
            || columns.iter().any(|b| b.len() != modes || b.iter().any(|x| !x.is_finite()))
            || initial_force_n.iter().any(|x| !x.is_finite()) {
            return Err(control("physical force input requires aligned nonempty finite columns and newtons"));
        }
        Ok(Self { columns, initial_force_n })
    }
}

/// One compiled held-force assignment, not an impulse or contact reaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProjectedForceEvent {
    /// Apply before this sample on the caller's performance clock.
    pub sample: u64,
    /// Original input-image order.
    pub voice: usize,
    /// Original mode order inside that image.
    pub mode: usize,
    /// Complete external modal force, after summing every physical actuator.
    pub force_n_per_sqrt_kg: f64,
}

/// Bounded, preprojected controls reusable by a retained mechanical owner.
/// Initial forces are separate from events at sample zero: a static preload
/// uses the former, then zero-time assignments may release it before stepping.
#[derive(Clone, Debug)]
pub struct ProjectedForceProgram {
    initial: Vec<Vec<f64>>,
    events: Vec<ProjectedForceEvent>,
}
impl ProjectedForceProgram {
    /// Compile all authored inputs before any mechanical state is touched.
    ///
    /// Unsorted inputs are legal. Equal-sample assignments retain input order;
    /// the last assignment to one port wins. All changes at that sample apply
    /// together BEFORE summation in original mode-then-port order. Only changed
    /// modal values emit controls, using the original bitwise-change criterion.
    /// Event, projection-work and emitted-control limits are the existing
    /// force compiler's. Time stepping and preload are deliberately absent.
    pub fn compile(inputs: Vec<PhysicalForceInput>, mut events: Vec<ModalForceEvent>,
        config: ForceRenderConfig) -> Result<Self, RenderError>
    {
        if inputs.is_empty() || config.sample_rate_hz == 0 || config.max_block == 0 {
            return Err(control("physical performance needs inputs, a positive sample rate and block capacity"));
        }
        if events.len() > config.max_events {
            return Err(sizing("physical force event budget exceeded"));
        }
        for event in &events {
            let input = inputs.get(event.voice)
                .ok_or(RenderError::UnknownVoice { index: event.voice })?;
            if event.sample == u64::MAX || event.port >= input.columns.len() || !event.force_n.is_finite() {
                return Err(control("force event needs a renderable sample, an existing port and finite newtons"));
            }
        }
        events.sort_by_key(|event| event.sample);
        let mut remaining = config.max_projection_terms;
        let mut initial = Vec::new();
        initial.try_reserve(inputs.len()).map_err(|_| sizing("cannot reserve initial force images"))?;
        let mut prepared = Vec::new();
        prepared.try_reserve(inputs.len()).map_err(|_| sizing("cannot reserve force inputs"))?;
        for input in inputs {
            charge(&mut remaining, &input.columns)?;
            let generalized = project(&input.columns, &input.initial_force_n)?;
            initial.push(generalized.clone());
            prepared.push((input.columns, input.initial_force_n, generalized));
        }
        let mut controls = Vec::new();
        let mut changed = Vec::new();
        let mut at = 0;
        while at < events.len() {
            let sample = events[at].sample;
            changed.clear();
            while at < events.len() && events[at].sample == sample {
                let event = events[at];
                prepared[event.voice].1[event.port] = event.force_n;
                changed.push(event.voice);
                at += 1;
            }
            changed.sort_unstable();
            changed.dedup();
            for &index in &changed {
                let (columns, force_n, previous) = &mut prepared[index];
                charge(&mut remaining, columns)?;
                let projected = project(columns, force_n)?;
                let count = projected.iter().zip(previous.iter())
                    .filter(|(a, b)| a.to_bits() != b.to_bits()).count();
                let total = controls.len().checked_add(count)
                    .ok_or_else(|| sizing("projected control count overflow"))?;
                if total > config.max_controls { return Err(sizing("projected modal control budget exceeded")); }
                controls.try_reserve(count).map_err(|_| sizing("cannot reserve projected controls"))?;
                for (mode, (&next, &old)) in projected.iter().zip(previous.iter()).enumerate() {
                    if next.to_bits() != old.to_bits() {
                        controls.push(ProjectedForceEvent { sample, voice: index, mode,
                            force_n_per_sqrt_kg: next });
                    }
                }
                *previous = projected;
            }
        }
        Ok(Self { initial, events: controls })
    }
    /// Complete initial external modal loads, in original image/mode order.
    #[must_use] pub fn initial_forces(&self) -> &[Vec<f64>] { &self.initial }
    /// Chronological assignments; all equal-sample entries apply as one group.
    #[must_use] pub fn events(&self) -> &[ProjectedForceEvent] { &self.events }
    /// Transfer prepared values without retaining another compiler or model.
    #[must_use]
    pub fn into_parts(self) -> (Vec<Vec<f64>>, Vec<ProjectedForceEvent>) { (self.initial, self.events) }
}
