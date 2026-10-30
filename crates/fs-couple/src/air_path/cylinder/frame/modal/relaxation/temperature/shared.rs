//! Original temperature claims for several bodies on one mechanical timeline.
//! Source evaluation is staged for ALL bindings before ONE native rate rebuild.
//! A later body's refusal must not leave an earlier body at new conditions.
use super::*;

/// Change every selected body's layer temperatures at one accepted boundary.
///
/// Bindings and input slices share declaration order; each input slice retains
/// its own original member/layer order. Native memory order may differ. Source
/// blocks must be disjoint, and every source must still match the complete
/// mechanical owner. This does not permit duplicated physical dissipation.
///
/// All selected source queries finish before one ordered native rate update.
/// Source snapshots publish only after native preparation/cancellation succeeds.
/// No displacement, fluid momentum, internal strain, clock or stored energy is
/// changed. Original scalar-storage and fixed-geometry restrictions apply.
/// For different fluid source types, the original models may borrow
/// `&dyn FluidTransport`; no transport coefficients are flattened or averaged.
/// This prescribes conditions; it does not solve heat transport.
pub fn set_layer_temperatures<F>(
    bindings: &mut [FrameRelaxationBinding<'_, '_, F>],
    mechanical: &mut CoupledModalSystem,
    layers: &[&[FrameRelaxationLayerTemperature<'_>]],
    gate: &CancelGate,
) -> Result<MaterialRelaxationTimeReport> {
    poll(gate)?;
    if layers.len() != bindings.len() {
        return Err(invalid("shared conditions require every original source binding"));
    }
    admit_bindings(bindings, mechanical, gate)?;
    let sources = bindings.iter().zip(layers).map(|(binding, input)|
        binding.prepare_layer_temperatures(mechanical, input, gate)
    ).collect::<Result<Vec<_>>>()?;
    let prepared = PreparedSharedTemperatures::new(sources, mechanical.material_memory().len(), gate)?;
    let report = mechanical.set_material_relaxation_times(&prepared.updates, gate)?;
    // The native owner made its final cancellation check. Only infallible
    // source moves may follow; do not introduce a second publication boundary.
    prepared.install(bindings);
    Ok(report)
}

/// Source snapshots stay in declaration order; the native batch is independently
/// sorted once, not applied body by body. Reused by finite caloric composition.
pub(super) struct PreparedSharedTemperatures {
    sources: Vec<PreparedLayerTemperatures>,
    pub(super) updates: Vec<MaterialRelaxationTime>,
}
impl PreparedSharedTemperatures {
    pub(super) fn new(sources: Vec<PreparedLayerTemperatures>, maximum_branches: usize,
        gate: &CancelGate) -> Result<Self> {
        poll(gate)?;
        let count = sources.iter().try_fold(0usize, |n, source|
            n.checked_add(source.updates.len())
                .ok_or_else(|| invalid("shared temperature update count overflow")))?;
        if count > maximum_branches {
            return Err(invalid("shared updates exceed the original complete material history"));
        }
        let mut updates = Vec::with_capacity(count);
        for source in &sources {
            poll(gate)?;
            updates.extend_from_slice(&source.updates);
        }
        updates.sort_unstable_by_key(|update| update.branch);
        if updates.iter().any(|update| update.branch >= maximum_branches)
            || updates.windows(2).any(|pair| pair[0].branch == pair[1].branch) {
            return Err(invalid("shared temperatures require disjoint original branch addresses"));
        }
        poll(gate)?;
        Ok(Self { sources, updates })
    }

    pub(super) fn install<F>(self, bindings: &mut [FrameRelaxationBinding<'_, '_, F>]) {
        // Only the constructors above and the thermal child can create this
        // proposal; both retain the unchanged binding count until publication.
        for (binding, source) in bindings.iter_mut().zip(self.sources) {
            binding.install_layer_temperatures(source);
        }
    }
}

/// O(B log B) interval admission, leaving original source declaration order
/// unchanged. Empty zero-increment blocks own no native memory and may coincide.
pub(super) fn admit_bindings<F>(bindings: &[FrameRelaxationBinding<'_, '_, F>],
    mechanical: &CoupledModalSystem, gate: &CancelGate) -> Result<()> {
    poll(gate)?;
    if bindings.is_empty() {
        return Err(invalid("shared material temperatures require at least one source binding"));
    }
    let mut ranges = Vec::with_capacity(bindings.len());
    for binding in bindings {
        poll(gate)?;
        let range = binding.memory_range();
        if range.start > range.end || range.end > mechanical.material_memory().len() {
            return Err(invalid("shared source lost its complete original memory range"));
        }
        if !range.is_empty() { ranges.push(range); }
    }
    ranges.sort_unstable_by_key(|range| range.start);
    if ranges.windows(2).any(|pair| pair[0].end > pair[1].start) {
        return Err(invalid("shared source memory ranges overlap"));
    }
    poll(gate)?;
    Ok(())
}
