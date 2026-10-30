//! Bounded waves of the existing calibrated paths on the native tile executor.

use core::mem::size_of;
use core::ops::ControlFlow;
use std::sync::Mutex;

use fs_alloc::OperationMemoryLease;
use fs_exec::{Cancelled, Cx, TileKernel, TilePlan};

use crate::photometry::PhotometryError;
use crate::tracer::{RenderExecutionConfig, RenderExecutionError, RenderWorkerPool, TracerError};

/// Raw samples and private moments remain bounded independently of image size.
const MAX_WAVE_PATHS: usize = 65_536;
const KERNEL: &str = "fs-render/photometric-path-wave-v1";

#[derive(Clone, Copy)]
pub(super) struct PathSample {
    pub xyz: [f64; 3],
    /// Only the first absolute sample of an uncached probe carries this value.
    pub direct_lux: Option<f64>,
}

impl PathSample {
    const EMPTY: Self = Self {
        xyz: [0.0; 3],
        direct_lux: None,
    };
}

struct WaveKernel<'a, M, F> {
    moments: &'a [M],
    sample: &'a F,
    tiles: &'a [Mutex<Vec<PathSample>>],
    failure: &'a Mutex<Option<(u64, PhotometryError)>>,
    first: u64,
    paths: usize,
    tile_paths: usize,
    active_tiles: usize,
}

impl<M, F> WaveKernel<'_, M, F> {
    fn fail(&self, address: u64, error: PhotometryError) -> ControlFlow<Cancelled, ()> {
        let mut failure = self
            .failure
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if failure
            .as_ref()
            .is_none_or(|(previous, _)| address < *previous)
        {
            *failure = Some((address, error));
        }
        ControlFlow::Break(Cancelled)
    }
}

impl<M, F> TileKernel for WaveKernel<'_, M, F>
where
    M: Sync,
    F: Fn(&Cx<'_>, u32, u32, &M) -> Result<PathSample, PhotometryError> + Sync,
{
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(KERNEL, self.active_tiles as u64)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Ok(tile) = usize::try_from(tile) else {
            return self.fail(self.first, internal("unrepresentable photometric tile"));
        };
        if tile >= self.active_tiles {
            return self.fail(self.first, internal("photometric tile outside wave"));
        }
        let offset = tile * self.tile_paths;
        let paths = self.tile_paths.min(self.paths - offset);
        // Each tile owns one slot for the complete dispatch. Holding this
        // tile-local mutex during tracing introduces no cross-worker lock.
        let mut output = self.tiles[tile]
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let outputs = self.moments.len() as u64;
        for local in 0..paths {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let address = self.first + (offset + local) as u64;
            let index = (address % outputs) as usize;
            let sample = (address / outputs) as u32;
            output[local] = match (self.sample)(cx, index as u32, sample, &self.moments[index]) {
                Ok(value) => value,
                Err(PhotometryError::Cancelled) => return ControlFlow::Break(Cancelled),
                Err(error) => return self.fail(address, error),
            };
        }
        if cx.checkpoint().is_err() {
            ControlFlow::Break(Cancelled)
        } else {
            ControlFlow::Continue(())
        }
    }
}

fn internal(reason: &'static str) -> PhotometryError {
    RenderExecutionError::Internal(reason).into()
}

fn bytes<T>(count: usize) -> Result<u64, PhotometryError> {
    u64::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<T>() as u64))
        .ok_or_else(|| internal("photometric execution storage size overflow"))
}

/// Run independent absolute path addresses, then replay their moment updates
/// in logical order. `sample` reads only the accepted pre-wave moment, while
/// `accumulate` mutates a private copy. A wave publishes only after every path
/// and every candidate moment succeeds and the original Cx remains live.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance<M>(
    cx: &Cx<'_>,
    max_path_samples: usize,
    execution: &RenderExecutionConfig,
    moments: &mut [M],
    retained_capacity: usize,
    completed: &mut u64,
    target: u64,
    sample: impl Fn(&Cx<'_>, u32, u32, &M) -> Result<PathSample, PhotometryError> + Sync,
    accumulate: impl Fn(&mut M, PathSample) -> Result<(), PhotometryError>,
) -> Result<usize, PhotometryError>
where
    M: Copy + Send + Sync,
{
    cx.checkpoint()?;
    let start = *completed;
    let limit = start
        .saturating_add(u64::try_from(max_path_samples).unwrap_or(u64::MAX))
        .min(target);
    if limit == start {
        return Ok(0);
    }
    if moments.is_empty() || retained_capacity < moments.len() || start > target {
        return Err(internal("invalid photometric wave state"));
    }
    let tile_paths =
        usize::try_from(u64::from(execution.tile_width()) * u64::from(execution.tile_height()))
            .map_err(|_| internal("unrepresentable photometric tile size"))?;
    let wave_limit = tile_paths
        .saturating_mul(execution.workers())
        .saturating_mul(2)
        .min(MAX_WAVE_PATHS);
    let capacity = (limit - start).min(wave_limit as u64) as usize;
    let tile_count = 1 + (capacity - 1) / tile_paths;
    let candidate_capacity = capacity.min(moments.len());
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes());
    // Borrowed scene/fixture/field assets are not copied or claimed by this
    // operation. Retained result payload, wave buffers and native executor
    // metadata are charged; allocator bookkeeping and thread stacks are not.
    let _retained = lease
        .reserve(
            "photometric-retained-moments",
            bytes::<M>(retained_capacity)?,
        )
        .map_err(RenderExecutionError::Memory)?;
    let candidate_bytes = bytes::<M>(candidate_capacity)?;
    let staging_bytes = bytes::<PathSample>(capacity)?
        .checked_add(bytes::<Mutex<Vec<PathSample>>>(tile_count)?)
        .and_then(|value| value.checked_add(candidate_bytes))
        .ok_or_else(|| internal("photometric wave storage size overflow"))?;
    let _staging = lease
        .reserve("photometric-wave-staging", staging_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let mut tiles = Vec::new();
    tiles
        .try_reserve_exact(tile_count)
        .map_err(|_| PhotometryError::AllocationRefused)?;
    for tile in 0..tile_count {
        cx.checkpoint()?;
        let count = tile_paths.min(capacity - tile * tile_paths);
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        values.resize(count, PathSample::EMPTY);
        tiles.push(Mutex::new(values));
    }
    let mut candidates = Vec::new();
    candidates
        .try_reserve_exact(candidate_capacity)
        .map_err(|_| PhotometryError::AllocationRefused)?;
    candidates.resize(candidate_capacity, moments[0]);
    cx.checkpoint()?;

    let pool = RenderWorkerPool::new(execution, cx.mode(), execution.run_id().0);
    pool.with_parked_crew_local(|crew| {
        while *completed < limit {
            cx.checkpoint()?;
            let first = *completed;
            let paths = (limit - first).min(capacity as u64) as usize;
            let active_tiles = 1 + (paths - 1) / tile_paths;
            let failure = Mutex::new(None);
            let outcome = {
                let kernel = WaveKernel {
                    moments: &*moments,
                    sample: &sample,
                    tiles: &tiles,
                    failure: &failure,
                    first,
                    paths,
                    tile_paths,
                    active_tiles,
                };
                crew.run_photometric_kernel(cx, &kernel, execution, &lease)
            };
            let failure = failure
                .into_inner()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            match outcome {
                Err(RenderExecutionError::Tracer(TracerError::Cancelled)) => {
                    return Err(failure.map_or(PhotometryError::Cancelled, |(_, error)| error));
                }
                Err(error) => return Err(error.into()),
                Ok(_) if failure.is_some() => {
                    return Err(internal("photometric tile failure disagreed with executor"));
                }
                Ok(_) => {}
            }
            cx.checkpoint()?;
            let outputs = moments.len();
            let first_index = (first % outputs as u64) as usize;
            let changed = paths.min(outputs);
            for (offset, candidate) in candidates[..changed].iter_mut().enumerate() {
                if offset % 256 == 0 {
                    cx.checkpoint()?;
                }
                *candidate = moments[(first_index + offset) % outputs];
            }
            let mut offset = 0;
            for tile in &mut tiles[..active_tiles] {
                cx.checkpoint()?;
                let values = tile
                    .get_mut()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                let count = tile_paths.min(paths - offset);
                for &value in &values[..count] {
                    if offset % 256 == 0 {
                        cx.checkpoint()?;
                    }
                    accumulate(&mut candidates[offset % outputs], value)?;
                    offset += 1;
                }
            }
            // The two bounded copies are the publication point. They cannot
            // fail and contain no checkpoints that could expose half a wave.
            cx.checkpoint()?;
            let first_count = changed.min(outputs - first_index);
            moments[first_index..first_index + first_count]
                .copy_from_slice(&candidates[..first_count]);
            moments[..changed - first_count].copy_from_slice(&candidates[first_count..changed]);
            *completed += paths as u64;
        }
        Ok((*completed - start) as usize)
    })
}
