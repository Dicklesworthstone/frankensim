//! Continue a bounded near wake by retaining its aged physical sources.

use super::{FarField, cell_from_segs_with_cores};
use crate::filament3d::{FilamentWake, Refusal};

/// Bounds on the complete retained far field after a transfer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FarRetirementLimits {
    /// Includes existing and newly appended cells, including empty cells.
    pub max_cells: usize,
    /// Includes all existing and newly retained physical source records.
    pub max_segments: usize,
}

/// A storage transfer at one physical instant; no time or diffusion advances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FarRetirement {
    pub rows_transferred: usize,
    pub near_rows_remaining: usize,
    pub cells_added: usize,
    pub far_cells: usize,
    pub segments_added: usize,
    pub far_segments: usize,
}

impl FarField {
    /// Transfer the oldest near rows into this existing retained far field.
    /// Call before the next shedding interval would exhaust the near-row cap.
    /// Positive `n_rows` must be a multiple of positive `rows_per_cell`, and
    /// at least one near row must remain. Zero rows is an admitted identity.
    ///
    /// The original physical-edge visitor captures each source exactly once,
    /// including its current evolved core and the trailing edges ending at
    /// the first remaining near row. Its captured downstream circulation
    /// becomes the new interface. Existing far sources, diffusion offsets,
    /// source-buffer allocations and expansion caches remain unchanged.
    /// New cells append in original oldest-to-newest source order. No vortex
    /// is discarded, circulation averaged, core history reset, or time elapsed.
    ///
    /// This permits repeated retirement while retaining the near certificate's
    /// row cap. Far storage still grows and must fit the explicit total limits;
    /// this is not a pruning or coarsening approximation. Both exact fields
    /// and geometric impulses are invariant to roundoff. Grouping newly moved
    /// sources can change multipole approximation error, so prepared source
    /// trees and separately stored image snapshots must be rebuilt.
    ///
    /// Source validation is linear in the existing inventory; only the newly
    /// retired sources allocate arrays or rebuild moments. No full near/far
    /// clone or global endpoint map is needed. All fallible preparation and
    /// cancellation precede publication. This preserves any existing nodal
    /// closure; it does not certify earlier pruning or repair an open wake.
    ///
    /// # Errors
    /// Invalid rows/grouping, incompatible shedding clocks, invalid sources,
    /// exceeded total storage caps, unrepresentable moments, allocation failure
    /// or cancellation. Both physical inputs remain unchanged on every error.
    pub fn retire_near_rows(
        &mut self,
        near: &mut FilamentWake,
        n_rows: usize,
        rows_per_cell: usize,
        limits: FarRetirementLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FarRetirement, Refusal> {
        checkpoint(&mut cancelled)?;
        if rows_per_cell == 0
            || n_rows % rows_per_cell != 0
            || (n_rows != 0 && n_rows >= near.rows.len())
        {
            return Err(refusal(
                "far-retirement-invalid",
                "retire complete positive-size cell groups while retaining at least one near row",
            ));
        }
        let new_cells = n_rows / rows_per_cell;
        let final_cells = self
            .cells
            .len()
            .checked_add(new_cells)
            .ok_or_else(capacity)?;
        if final_cells > limits.max_cells {
            return Err(capacity());
        }
        let mut final_segments = 0_usize;
        for cell in &self.cells {
            checkpoint(&mut cancelled)?;
            final_segments = final_segments
                .checked_add(cell.segs.len())
                .ok_or_else(capacity)?;
            if final_segments > limits.max_segments {
                return Err(capacity());
            }
        }
        self.validate_kinematics_sources(near)?;
        let growth = self.near_core_growth_m2().ok_or_else(|| {
            refusal(
                "far-retirement-invalid",
                "invalid retained core law or shedding interval",
            )
        })?;
        checkpoint(&mut cancelled)?;

        let before_segments = final_segments;
        let mut counts = Vec::new();
        reserve(&mut counts, new_cells)?;
        for first in (0..n_rows).step_by(rows_per_cell) {
            let mut count = 0_usize;
            for row in first..first + rows_per_cell {
                checkpoint(&mut cancelled)?;
                near.visit_row_segments(row, growth, &mut |_| count += 1);
            }
            final_segments = final_segments.checked_add(count).ok_or_else(capacity)?;
            if final_segments > limits.max_segments {
                return Err(capacity());
            }
            counts.push(count);
        }
        let report = FarRetirement {
            rows_transferred: n_rows,
            near_rows_remaining: near.rows.len() - n_rows,
            cells_added: new_cells,
            far_cells: final_cells,
            segments_added: final_segments - before_segments,
            far_segments: final_segments,
        };
        if n_rows == 0 {
            checkpoint(&mut cancelled)?;
            return Ok(report);
        }

        let mut prepared = Vec::new();
        reserve(&mut prepared, new_cells)?;
        for (index, count) in counts.into_iter().enumerate() {
            let mut segments = Vec::new();
            let mut cores = Vec::new();
            reserve(&mut segments, count)?;
            reserve(&mut cores, count)?;
            let first = index * rows_per_cell;
            for row in first..first + rows_per_cell {
                checkpoint(&mut cancelled)?;
                near.visit_row_segments(row, growth, &mut |s| {
                    segments.push((s.start, s.end, s.circulation));
                    // Bake only this source's elapsed history into its base.
                    // Future common growth starts from a zero new-cell offset.
                    cores.push(s.core2_m2);
                });
            }
            checkpoint(&mut cancelled)?;
            let cell = cell_from_segs_with_cores(segments, cores);
            if !cell.expansion().is_finite() {
                return Err(refusal(
                    "far-retirement-unrepresentable",
                    "retired source geometry or core moments exceed the numerical range",
                ));
            }
            prepared.push(cell);
        }
        let mut interface = Vec::new();
        reserve(&mut interface, near.cert.n_stations)?;
        interface.extend_from_slice(&near.rows[n_rows - 1].gamma);
        checkpoint(&mut cancelled)?;
        // Reserve the destination only after the final cancellation point.
        // Failure changes neither source inventory; success leaves only
        // infallible moves and bounded near-row header compaction below.
        self.cells
            .try_reserve(new_cells)
            .map_err(|_| allocation())?;
        self.cells.append(&mut prepared);
        near.rows.drain(..n_rows);
        near.downstream_gamma = Some(interface);
        Ok(report)
    }
}

fn reserve<T>(values: &mut Vec<T>, count: usize) -> Result<(), Refusal> {
    values.try_reserve_exact(count).map_err(|_| allocation())
}

fn capacity() -> Refusal {
    refusal(
        "far-retirement-capacity",
        "complete retained far field exceeds the declared capacity",
    )
}

fn allocation() -> Refusal {
    refusal(
        "far-retirement-allocation",
        "retired source storage is not representable or available",
    )
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "far-retirement-cancelled",
            "near-to-far transfer cancelled before publication",
        ))
    } else {
        Ok(())
    }
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "retain a valid common-clock wake and sufficient explicit storage".into(),
        ],
    }
}
