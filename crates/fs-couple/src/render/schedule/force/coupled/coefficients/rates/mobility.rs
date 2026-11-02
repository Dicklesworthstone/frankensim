//! Native diagonal resistance changes composed with original Maxwell rates.
use super::*;

/// Source-derived replacement of one original modal resistance, after the
/// geometry owner has retained all unchanged material/ambient contributions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ModalDampingRatio {
    pub mode: usize,
    pub previous: f64,
    pub next: f64,
}

impl CoupledModalSystem {
    pub(super) fn prepare_modal_mobility(
        &self,
        updates: &[ModalDampingRatio],
        gate: &CancelGate,
    ) -> Result<
        Option<(Vec<ModalAcousticTimeModel>, Vec<ModalAcousticTimeModel>)>,
        ModalCouplingError,
    > {
        if updates.is_empty() {
            return Ok(None);
        }
        let modes: Vec<_> = self.models.iter().flat_map(|m| m.modes()).collect();
        let mut previous = None;
        let mut changed = false;
        for update in updates {
            poll(Some(gate))?;
            let original = modes
                .get(update.mode)
                .ok_or_else(|| invalid("material mobility names no original modal coordinate"))?;
            if previous.is_some_and(|p| update.mode <= p)
                || original.damping_ratio.to_bits() != update.previous.to_bits()
                || !update.next.is_finite()
                || update.next < 0.0
            {
                return Err(invalid(
                    "material mobility must retain ordered original modes and nonnegative finite resistance",
                ));
            }
            previous = Some(update.mode);
            changed |= update.previous != update.next;
        }
        if !changed {
            return Ok(None);
        }
        let mut accepted = self.models.clone();
        let mut staged = self.candidates.clone();
        for component in 0..self.models.len() {
            poll(Some(gate))?;
            let begin = self.offsets[component];
            let end = self.offsets[component + 1];
            let selected: Vec<_> = updates
                .iter()
                .filter(|u| (begin..end).contains(&u.mode))
                .collect();
            if selected.is_empty() {
                continue;
            }
            let mut damping: Vec<_> = self.models[component]
                .modes()
                .iter()
                .map(|m| m.damping_ratio)
                .collect();
            for update in selected {
                damping[update.mode - begin] = update.next;
            }
            accepted[component] = self.models[component]
                .with_damping_ratios(&damping)
                .map_err(|source| ModalCouplingError::Component { component, source })?;
            staged[component] = self.candidates[component]
                .with_damping_ratios(&damping)
                .map_err(|source| ModalCouplingError::Component { component, source })?;
        }
        Ok(Some((accepted, staged)))
    }
}
