//! Compact immutable optical metadata on the original simulation cell order.

use super::{VolumeTransportError, checkpoint};
use fs_exec::Cx;
use std::sync::Arc;

/// Explicit limits for optical model storage and the cell-index scan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CellPaletteBudget {
    /// Maximum retained optical models, including unreferenced entries.
    pub max_models: usize,
    /// Maximum cell indices scanned and retained.
    pub max_cells: usize,
}

/// Shared optical models and one immutable `u32` model index per physical cell.
///
/// Values and indices use `Arc`; cloning or sharing a map never copies density
/// or temperature fields. Indices retain the native grid's x-major order.
/// Models are not deduplicated or reordered: the caller owns their declaration
/// order, and the transport owner validates their physical wavelength support.
#[derive(Debug, PartialEq)]
pub struct CellPalette<T> {
    values: Arc<[T]>,
    cell_indices: Arc<[u32]>,
    largest_index: u32,
}

impl<T> Clone for CellPalette<T> {
    fn clone(&self) -> Self {
        Self {
            values: Arc::clone(&self.values),
            cell_indices: Arc::clone(&self.cell_indices),
            largest_index: self.largest_index,
        }
    }
}

impl<T> CellPalette<T> {
    /// Admit nonempty values and a nonempty map whose every index is valid.
    /// Storage/work limits are checked before scanning; cancellation is polled
    /// at most every 1024 indices. `cell_indices` is retained without copying.
    pub fn try_new(
        cx: &Cx<'_>,
        values: Vec<T>,
        cell_indices: Arc<[u32]>,
        budget: CellPaletteBudget,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        validate_lengths(values.len(), cell_indices.len(), budget)?;
        let mut largest_index = 0;
        for chunk in cell_indices.chunks(1024) {
            checkpoint(cx)?;
            for &index in chunk {
                if index as usize >= values.len() {
                    return Err(VolumeTransportError::Palette);
                }
                largest_index = largest_index.max(index);
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            values: values.into(),
            cell_indices,
            largest_index,
        })
    }

    /// Bind another set of optical models to this exact immutable cell map.
    ///
    /// The cached largest validated index admits the new values without
    /// rescanning or copying the map. Entries beyond that index may differ in
    /// number, but every referenced entry must exist. This is useful when
    /// absorption and scattering share one composition map.
    pub fn with_values<U>(
        &self,
        cx: &Cx<'_>,
        values: Vec<U>,
        budget: CellPaletteBudget,
    ) -> Result<CellPalette<U>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_lengths(values.len(), self.cell_indices.len(), budget)?;
        if self.largest_index as usize >= values.len() {
            return Err(VolumeTransportError::Palette);
        }
        Ok(CellPalette {
            values: values.into(),
            cell_indices: Arc::clone(&self.cell_indices),
            largest_index: self.largest_index,
        })
    }

    /// Original optical models in declaration order.
    #[must_use]
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// Shared native-order model indices; `Arc::clone` preserves the same map.
    #[must_use]
    pub const fn cell_indices(&self) -> &Arc<[u32]> {
        &self.cell_indices
    }

    /// Model index for one physical cell; out-of-range cell access is explicit.
    #[must_use]
    pub fn index_at_cell(&self, cell: usize) -> Option<u32> {
        self.cell_indices.get(cell).copied()
    }

    /// Optical model for one physical cell, without interpolation or allocation.
    #[must_use]
    pub fn at_cell(&self, cell: usize) -> Option<&T> {
        self.index_at_cell(cell)
            .map(|index| &self.values[index as usize])
    }
}

fn validate_lengths(
    models: usize,
    cells: usize,
    budget: CellPaletteBudget,
) -> Result<(), VolumeTransportError> {
    if models > budget.max_models || cells > budget.max_cells {
        return Err(VolumeTransportError::Budget);
    }
    if models == 0 || cells == 0 || models > u32::MAX as usize {
        return Err(VolumeTransportError::Palette);
    }
    Ok(())
}
