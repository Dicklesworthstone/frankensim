//! Resolve retained material filaments without changing their present field.

use super::{FarCell, FarField, cell_from_segs_with_cores, dist};
use crate::filament3d::{Refusal, preserves_straight_segment};

/// Geometric resolution change; no physical time or diffusion has elapsed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FarSegmentRefinement {
    /// Retained physical segments before refinement, including zero strengths.
    pub segments_before: usize,
    /// Retained physical segments after refinement.
    pub segments_after: usize,
    /// Cells whose source arrays and expansions were rebuilt.
    pub cells_refined: usize,
    /// Longest retained segment before refinement [m].
    pub longest_segment_before_m: f64,
    /// Longest retained segment after refinement [m].
    pub longest_segment_after_m: f64,
}

impl FarField {
    /// Subdivide stretched retained filaments into material segments no longer
    /// than `max_segment_length_m`, within 64 machine epsilons of the target.
    /// New vertices participate in the existing coupled near/far advection,
    /// allowing subsequent nonuniform flow to bend previously long segments.
    ///
    /// Every child keeps its parent's signed circulation and BASE core. The
    /// accumulated diffusion offset is retained separately, without resetting
    /// elapsed core growth. Original endpoints are copied exactly, preserving
    /// shared near/far and inter-cell nodes. Interior children meet at identical
    /// coordinates and carry equal circulation, preserving nodal closure.
    /// Exact velocity, its Jacobian and both geometric impulses are unchanged
    /// to roundoff at refinement; no missing historical curvature is recovered.
    ///
    /// Only changed cells allocate replacement sources and rebuild their
    /// multipole caches. Unchanged cells retain their buffers and bits. Work is
    /// linear in original plus refined source count. `max_segments` bounds the
    /// COMPLETE output source count and separately the cell count, including
    /// empty cells. All replacements remain private until the final checkpoint;
    /// any refusal leaves this owner unchanged. Rebuild cached image snapshots
    /// and prepared trees after success, as for other geometry changes.
    ///
    /// # Errors
    /// Invalid target, mode, interval or sources; insufficient source capacity;
    /// allocation failure; unrepresentable subdivision/moments; cancellation.
    pub fn refine_segments(
        &mut self,
        max_segment_length_m: f64,
        max_segments: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FarSegmentRefinement, Refusal> {
        checkpoint(&mut cancelled)?;
        self.mode.admit()?;
        if !max_segment_length_m.is_finite()
            || max_segment_length_m <= 0.0
            || !self.shed_dt_s.is_finite()
            || self.shed_dt_s <= 0.0
        {
            return Err(refusal(
                "far-refinement-invalid",
                "segment target and retained shedding interval must be finite and positive",
            ));
        }
        if self.cells.len() > max_segments {
            return Err(capacity());
        }
        let mut before = 0_usize;
        for cell in &self.cells {
            checkpoint(&mut cancelled)?;
            before = before.checked_add(cell.segs.len()).ok_or_else(capacity)?;
            if before > max_segments {
                return Err(capacity());
            }
        }
        let mut after = before;
        let mut changed = 0;
        let mut longest_before = 0.0_f64;
        let mut longest_after = 0.0_f64;
        for cell in &self.cells {
            checkpoint(&mut cancelled)?;
            if cell.segs.len() != cell.source_core2_m2.len()
                || !cell.diffusion_offset_m2.is_finite()
                || cell.diffusion_offset_m2 < 0.0
            {
                return Err(invalid_sources());
            }
            let mut cell_changed = false;
            for (i, (&(a, b, gamma), &core)) in
                cell.segs.iter().zip(&cell.source_core2_m2).enumerate()
            {
                if i % 256 == 0 {
                    checkpoint(&mut cancelled)?;
                }
                if !gamma.is_finite()
                    || !a.iter().chain(&b).all(|v| v.is_finite())
                    || !core.is_finite()
                    || core < 0.0
                    || !(core + cell.diffusion_offset_m2).is_finite()
                {
                    return Err(invalid_sources());
                }
                let length = checked_length(a, b)?;
                longest_before = longest_before.max(length);
                let count = pieces(length, max_segment_length_m, max_segments - after)?;
                after = after.checked_add(count - 1).ok_or_else(capacity)?;
                if count > 1 {
                    cell_changed = true;
                } else {
                    longest_after = longest_after.max(length);
                }
            }
            changed += usize::from(cell_changed);
        }
        let mut report = FarSegmentRefinement {
            segments_before: before,
            segments_after: after,
            cells_refined: changed,
            longest_segment_before_m: longest_before,
            longest_segment_after_m: longest_after,
        };
        if changed == 0 {
            checkpoint(&mut cancelled)?;
            return Ok(report);
        }
        let mut replacements: Vec<(usize, FarCell)> = Vec::new();
        reserve(&mut replacements, changed)?;
        for (cell_index, cell) in self.cells.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            let mut cell_count = 0_usize;
            for (i, &(a, b, _)) in cell.segs.iter().enumerate() {
                if i % 256 == 0 {
                    checkpoint(&mut cancelled)?;
                }
                cell_count = cell_count
                    .checked_add(pieces(dist(a, b), max_segment_length_m, max_segments)?)
                    .ok_or_else(capacity)?;
            }
            if cell_count == cell.segs.len() {
                continue;
            }
            let mut segments = Vec::new();
            let mut cores = Vec::new();
            reserve(&mut segments, cell_count)?;
            reserve(&mut cores, cell_count)?;
            for (&(a, b, gamma), &core) in cell.segs.iter().zip(&cell.source_core2_m2) {
                let count = pieces(dist(a, b), max_segment_length_m, max_segments)?;
                let mut start = a;
                for part in 1..=count {
                    if part % 256 == 1 {
                        checkpoint(&mut cancelled)?;
                    }
                    // Preserve the original endpoint bit-for-bit. A shared
                    // vertex must not acquire two rounded coordinate identities.
                    let end = if part == count {
                        b
                    } else {
                        let fraction = part as f64 / count as f64;
                        std::array::from_fn(|axis| a[axis] + fraction * (b[axis] - a[axis]))
                    };
                    let length = checked_length(start, end)?;
                    if !within_target(length / max_segment_length_m)
                        || (count > 1 && !preserves_straight_segment(a, b, start, end))
                    {
                        return Err(unrepresentable());
                    }
                    longest_after = longest_after.max(length);
                    segments.push((start, end, gamma));
                    cores.push(core);
                    start = end;
                }
            }
            checkpoint(&mut cancelled)?;
            let mut rebuilt = cell_from_segs_with_cores(segments, cores);
            rebuilt.diffusion_offset_m2 = cell.diffusion_offset_m2;
            rebuilt.core2_m2 += cell.diffusion_offset_m2;
            if !rebuilt.expansion().is_finite() {
                return Err(unrepresentable());
            }
            replacements.push((cell_index, rebuilt));
        }
        report.longest_segment_after_m = longest_after;
        checkpoint(&mut cancelled)?;
        for (index, cell) in replacements {
            self.cells[index] = cell;
        }
        Ok(report)
    }
}

fn pieces(length: f64, target: f64, additional: usize) -> Result<usize, Refusal> {
    let ratio = length / target;
    if within_target(ratio) {
        return Ok(1);
    }
    let count = ratio.ceil();
    if !count.is_finite() || count >= usize::MAX as f64 || count - 1.0 > additional as f64 {
        return Err(capacity());
    }
    let count = count as usize;
    if count.checked_sub(1).is_none_or(|n| n > additional) {
        return Err(capacity());
    }
    Ok(count)
}

fn within_target(ratio: f64) -> bool {
    ratio <= 1.0 + 64.0 * f64::EPSILON
}

fn checked_length(a: [f64; 3], b: [f64; 3]) -> Result<f64, Refusal> {
    let length = dist(a, b);
    if length.is_finite() {
        Ok(length)
    } else {
        Err(unrepresentable())
    }
}

fn reserve<T>(buffer: &mut Vec<T>, count: usize) -> Result<(), Refusal> {
    buffer.try_reserve_exact(count).map_err(|_| {
        refusal(
            "far-refinement-allocation",
            "replacement filament storage could not be reserved",
        )
    })
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "far-refinement-cancelled",
            "refinement was cancelled",
        ))
    } else {
        Ok(())
    }
}

fn capacity() -> Refusal {
    refusal(
        "far-refinement-budget",
        "complete retained source or cell count exceeds the declared capacity",
    )
}

fn invalid_sources() -> Refusal {
    refusal(
        "far-refinement-invalid",
        "retained source positions, circulation and per-source cores must be finite and valid",
    )
}

fn unrepresentable() -> Refusal {
    refusal(
        "far-refinement-unrepresentable",
        "target subdivision, finite source length or rebuilt moments are not representable",
    )
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "retain finite source geometry and cores, a resolvable target and sufficient capacity"
                .into(),
        ],
    }
}
