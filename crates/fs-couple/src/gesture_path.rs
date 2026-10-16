//! Exact typed-path bindings to existing pressure and physical bow consumers.
//!
//! Every source path must have an explicit destination. Unsupported inputs
//! return NoConsumer by path. V2 bytes, paths and hash remain intact through
//! lowering; only legacy entry points interpret numeric alias metadata.
use crate::render::schedule::force::coupled::render::friction::{
    CompiledFrictionTrack, FrictionGestureConfig, compile_friction_track_admitted,
};
use crate::render::schedule::{
    GestureCompileError, PressureGestureBinding, ScheduledControl,
    compile_pressure_gestures_admitted,
};
use fs_scenario::assembly::ParameterPath;
use fs_scenario::gesture::path::{legacy_member, path_name};
use fs_scenario::gesture::{
    GestureError, GestureTarget, GestureTrack, GestureValue, ResolvedSchedule,
};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathInput {
    BlowingPressure {
        voice: usize,
    },
    /// Actual material-string member, not an encoded MIDI or u8 alias index.
    Bow {
        member: usize,
    },
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathBinding {
    pub path: ParameterPath,
    pub input: PathInput,
}
#[derive(Clone, Debug)]
pub struct BoundPaths {
    source: ResolvedSchedule,
    bindings: Vec<PathBinding>,
}
fn no_consumer(path: &ParameterPath) -> GestureError {
    GestureError::NoConsumer {
        path: path_name(path),
    }
}
/// Admission only; typed callers supply the exact description-admitted path.
pub(crate) fn matches_input(
    track: &GestureTrack,
    input: PathInput,
    path: Option<&ParameterPath>,
) -> bool {
    if let Some(path) = path {
        matches!(&track.target,GestureTarget::Parameter(p) if p==path)
            && std::iter::once(&track.initial)
                .chain(track.events.iter().map(|e| &e.value))
                .all(|v| match input {
                    PathInput::BlowingPressure { .. } => matches!(v, GestureValue::PressurePa(_)),
                    PathInput::Bow { .. } => matches!(v, GestureValue::Bow { .. }),
                })
    } else {
        match input {
            PathInput::BlowingPressure { .. } => {
                legacy_member(&track.target, "blowing_pressure") == Some(0)
            }
            PathInput::Bow { .. } => legacy_member(&track.target, "bow").is_some(),
        }
    }
}
/// Validate the COMPLETE binding before compiling any input or touching state.
pub fn bind_paths(
    source: &ResolvedSchedule,
    bindings: &[PathBinding],
) -> Result<BoundPaths, GestureError> {
    let mut by_path = BTreeMap::new();
    for b in bindings {
        source.track(&b.path)?;
        if by_path.insert(b.path.clone(), b.input).is_some() {
            return Err(GestureError::DuplicateParameter {
                path: path_name(&b.path),
            });
        }
    }
    for track in source.schedule().tracks() {
        let GestureTarget::Parameter(path) = &track.target else {
            return Err(GestureError::Invalid {
                what: "resolved schedule contains a legacy alias",
            });
        };
        let input = by_path.get(path).ok_or_else(|| no_consumer(path))?;
        if !matches_input(track, *input, Some(path)) {
            return Err(no_consumer(path));
        }
    }
    Ok(BoundPaths {
        source: source.clone(),
        bindings: bindings.to_vec(),
    })
}
impl BoundPaths {
    pub fn source(&self) -> &ResolvedSchedule {
        &self.source
    }
    pub fn bindings(&self) -> &[PathBinding] {
        &self.bindings
    }
    /// Lower all pressure paths to the original sample-clock compiler. Bow
    /// bindings require their bow consumer and refuse here rather than vanish.
    pub fn pressure_controls(
        &self,
        audio_rate_hz: u32,
        samples: u64,
        max_work: u64,
    ) -> Result<Vec<ScheduledControl>, GestureCompileError> {
        let mut bindings = Vec::new();
        for b in &self.bindings {
            let PathInput::BlowingPressure { voice } = b.input else {
                return Err(GestureCompileError::Gesture(no_consumer(&b.path)));
            };
            let track = self
                .source
                .track(&b.path)
                .map_err(GestureCompileError::Gesture)?;
            bindings.push(PressureGestureBinding {
                track: track.id.clone(),
                voice,
            });
        }
        compile_pressure_gestures_admitted(
            self.source.schedule(),
            &bindings,
            audio_rate_hz,
            samples,
            max_work,
            Some(&self.bindings),
        )
    }
    /// Lower one explicitly bound fixed physical contact through the existing
    /// bow compiler. The caller executes each bound contact in its joint model.
    pub fn friction_controls(
        &self,
        path: &ParameterPath,
        config: FrictionGestureConfig,
    ) -> Result<CompiledFrictionTrack, GestureCompileError> {
        if !self
            .bindings
            .iter()
            .any(|b| &b.path == path && matches!(b.input, PathInput::Bow { .. }))
        {
            return Err(GestureCompileError::Gesture(no_consumer(path)));
        }
        let track = self
            .source
            .track(path)
            .map_err(GestureCompileError::Gesture)?;
        compile_friction_track_admitted(self.source.schedule(), &track.id, config, Some(path))
    }
}
