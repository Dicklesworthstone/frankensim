//! Model-independent, right-continuous histories of authored nuclide movement.
//!
//! The history admits rate values, clock ordering and fixed directed paths.
//! Its physical owner additionally binds original nuclides, compartments and
//! volumes. Frames retain their input vectors once; choosing a frame allocates
//! nothing and never averages rates across a discontinuity.

use super::NuclideTransport;
use crate::nonnegative;

/// Replace all path rates at this absolute time, retaining original path order.
#[derive(Clone, Debug, PartialEq)]
pub struct NuclideTransportChange {
    pub time_s: f64,
    pub transport: NuclideTransport,
}

/// A shared piecewise-constant movement history, independent of a nuclear model.
#[derive(Clone, Debug, PartialEq)]
pub struct NuclideTransportTimeline {
    initial: NuclideTransport,
    changes: Vec<NuclideTransportChange>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NuclideTransportTimelineError {
    Invalid(&'static str),
    Budget(&'static str),
}

impl core::fmt::Display for NuclideTransportTimelineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nuclide transport timeline refused: {self:?}")
    }
}

impl std::error::Error for NuclideTransportTimelineError {}

impl NuclideTransportTimeline {
    pub fn new(
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        maximum_changes: usize,
    ) -> Result<Self, NuclideTransportTimelineError> {
        Self::new_with_checkpoint(initial, changes, maximum_changes, || {
            Ok::<(), NuclideTransportTimelineError>(())
        })
    }

    /// Poll each frame and at most 256 coefficients apart. The callback's error
    /// type is retained, so physical owners keep their original cancellation or
    /// admission reason without a second complete rate-validation pass.
    pub fn new_with_checkpoint<E: From<NuclideTransportTimelineError>>(
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        maximum_changes: usize,
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<Self, E> {
        checkpoint()?;
        if changes.len() > maximum_changes {
            return Err(
                NuclideTransportTimelineError::Budget("nuclide movement schedule changes").into(),
            );
        }
        admit_rows(&initial, &mut checkpoint)?;
        let mut previous = 0.0;
        for change in &changes {
            checkpoint()?;
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(NuclideTransportTimelineError::Invalid(
                    "movement changes need strictly increasing positive times",
                )
                .into());
            }
            admit_frame(&initial, &change.transport, &mut checkpoint)?;
            previous = change.time_s;
        }
        checkpoint()?;
        Ok(Self { initial, changes })
    }

    /// A physical owner has already admitted these constant rate vectors.
    #[cfg(feature = "thermal")]
    pub(crate) fn admitted_constant(initial: NuclideTransport) -> Self {
        Self {
            initial,
            changes: Vec::new(),
        }
    }

    pub fn initial(&self) -> &NuclideTransport {
        &self.initial
    }

    pub fn changes(&self) -> &[NuclideTransportChange] {
        &self.changes
    }

    /// Zero selects the initial frame; k>0 selects change k-1.
    pub fn position_at(&self, time_s: f64) -> Result<usize, NuclideTransportTimelineError> {
        if !nonnegative(time_s) {
            return Err(NuclideTransportTimelineError::Invalid(
                "movement schedule time must be finite and nonnegative",
            ));
        }
        Ok(self
            .changes
            .partition_point(|change| change.time_s <= time_s))
    }

    pub fn at_position(&self, position: usize) -> Option<&NuclideTransport> {
        if position == 0 {
            Some(&self.initial)
        } else {
            self.changes
                .get(position - 1)
                .map(|change| &change.transport)
        }
    }

    /// Right-continuous rates at a finite nonnegative absolute time [s].
    pub fn transport_at(
        &self,
        time_s: f64,
    ) -> Result<&NuclideTransport, NuclideTransportTimelineError> {
        // A partition point is always an admitted frame position.
        let position = self.position_at(time_s)?;
        Ok(if position == 0 {
            &self.initial
        } else {
            &self.changes[position - 1].transport
        })
    }

    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        self.changes
            .get(self.position_at(time_s).ok()?)
            .map(|change| change.time_s)
    }
}

fn rows(transport: &NuclideTransport) -> impl Iterator<Item = &[f64]> {
    transport
        .transfers
        .iter()
        .map(|path| path.rates_per_s.as_slice())
        .chain(
            transport
                .outflows
                .iter()
                .map(|path| path.rates_per_s.as_slice()),
        )
}

fn admit_rows<E: From<NuclideTransportTimelineError>>(
    transport: &NuclideTransport,
    checkpoint: &mut impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    for row in rows(transport) {
        admit_row(row, checkpoint)?;
    }
    Ok(())
}

fn admit_frame<E: From<NuclideTransportTimelineError>>(
    original: &NuclideTransport,
    current: &NuclideTransport,
    checkpoint: &mut impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    let invalid = || {
        NuclideTransportTimelineError::Invalid(
            "movement histories require fixed ordered paths and coefficient shapes",
        )
    };
    if original.transfers.len() != current.transfers.len()
        || original.outflows.len() != current.outflows.len()
    {
        return Err(invalid().into());
    }
    for (before, after) in original.transfers.iter().zip(&current.transfers) {
        checkpoint()?;
        if (before.from_compartment, before.to_compartment)
            != (after.from_compartment, after.to_compartment)
            || before.rates_per_s.len() != after.rates_per_s.len()
        {
            return Err(invalid().into());
        }
        admit_row(&after.rates_per_s, checkpoint)?;
    }
    for (before, after) in original.outflows.iter().zip(&current.outflows) {
        checkpoint()?;
        if before.compartment != after.compartment
            || before.rates_per_s.len() != after.rates_per_s.len()
        {
            return Err(invalid().into());
        }
        admit_row(&after.rates_per_s, checkpoint)?;
    }
    Ok(())
}

fn admit_row<E: From<NuclideTransportTimelineError>>(
    row: &[f64],
    checkpoint: &mut impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    checkpoint()?;
    for chunk in row.chunks(256) {
        checkpoint()?;
        if chunk.iter().any(|&rate| !nonnegative(rate)) {
            return Err(NuclideTransportTimelineError::Invalid(
                "movement coefficients must be finite and nonnegative",
            )
            .into());
        }
    }
    Ok(())
}
