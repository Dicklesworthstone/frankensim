//! Lowering of the M6.1-W2 measured generator (`fs_scenario`'s
//! `ArchingSurface`). Measured arching and graduation have one owner: these
//! constructors take its own degree-one graphs and its declared reference
//! state. Nothing is refitted here, and the caller declares no load state.

use super::{ArchingParameters, ArchingSurface, ReferenceLoadState, ThicknessField, bad};
use crate::PlateError;
use fs_scenario::assembly::ReferenceState;
use fs_scenario::assembly::generators::ArchingSurface as MeasuredArching;

impl ArchingSurface {
    /// The generator's height graph, with its source's load state. A loaded
    /// state is described by the source label and the declared state.
    pub fn from_measured(source: &MeasuredArching) -> Result<Self, PlateError> {
        let reference_state = match source.reference_state() {
            ReferenceState::Unloaded => ReferenceLoadState::Unloaded,
            state => ReferenceLoadState::Loaded {
                description: format!("{} ({state:?})", source.label()),
            },
        };
        Self::new(
            ArchingParameters::NurbsHeightMap(source.surface().clone()),
            reference_state,
        )
    }
}

impl ThicknessField {
    /// The generator's measured graduation graph. A source without a measured
    /// thickness table refuses; no thickness is defaulted.
    pub fn from_measured(source: &MeasuredArching) -> Result<Self, PlateError> {
        let chart = source
            .thickness_surface()
            .ok_or(bad("measured arching source has no thickness table"))?;
        Self::new(ArchingParameters::NurbsHeightMap(chart.clone()))
    }
}
