//! Original normal-pressure data for a moving shear brush.
//!
//! Cells are held in the physical footprint coordinates, independent of the
//! material/thermal mesh. Their actual pressure is never normalized to repair
//! a load mismatch. Zero-pressure gaps require a separate contact topology;
//! this owner represents one completely loaded, connected footprint.

use super::{MovingPartialSlipError, NormalPatchView, close, finite};
use fs_math::dd::Dd;

/// One original constant-pressure interval across the full strip width.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingPressureCell {
    /// Physical interval in the route's footprint coordinates [m].
    pub span_m: [f64; 2],
    pub normal_pressure_pa: f64,
}

/// Supplied normal-pressure discretization with its original load authority.
///
/// Retaining a caller/source identity does not establish measured calibration.
/// Spatial pressure resolution and thermal resolution are independent. The
/// caller must refine a sampled smooth pressure field before making a resolved
/// edge-traction claim; no interpolation or extrapolation is supplied here.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingPressureProfile {
    normal: NormalPatchView,
    cells: Vec<MovingPressureCell>,
    line_load_n_per_m: f64,
}

impl MovingPressureProfile {
    /// Bind finite positive, ordered, exactly contiguous pressure cells to the
    /// original normal patch. Footprint coverage and total force are checked
    /// when attached to the original moving route, whose width supplies area.
    pub fn new(
        normal: &NormalPatchView,
        cells: Vec<MovingPressureCell>,
    ) -> Result<Self, MovingPartialSlipError> {
        if cells.is_empty()
            || cells.iter().any(|cell| {
                !cell.span_m.iter().all(|v| v.is_finite())
                    || cell.span_m[0] >= cell.span_m[1]
                    || !cell.normal_pressure_pa.is_finite()
                    || cell.normal_pressure_pa <= 0.0
            })
            || cells
                .windows(2)
                .any(|pair| pair[0].span_m[1] != pair[1].span_m[0])
        {
            return Err(MovingPartialSlipError::Invalid(
                "finite positive contiguous pressure cells in original footprint order required",
            ));
        }
        let line_load = cells.iter().fold(Dd::ZERO, |sum, cell| {
            sum + (Dd::from_f64(cell.span_m[1]) - Dd::from_f64(cell.span_m[0]))
                * Dd::from_f64(cell.normal_pressure_pa)
        });
        let line_load_n_per_m = finite(line_load.to_f64())?;
        if line_load_n_per_m <= 0.0 {
            return Err(MovingPartialSlipError::Invalid(
                "positive representable integrated normal pressure required",
            ));
        }
        Ok(Self {
            normal: normal.clone(),
            cells,
            line_load_n_per_m,
        })
    }

    pub fn cells(&self) -> &[MovingPressureCell] {
        &self.cells
    }

    pub const fn normal_patch(&self) -> &NormalPatchView {
        &self.normal
    }

    pub(super) fn validate_for(
        &self,
        normal: &NormalPatchView,
        footprint: [f64; 2],
        width_m: f64,
    ) -> Result<(), MovingPartialSlipError> {
        if &self.normal != normal {
            return Err(MovingPartialSlipError::Invalid(
                "pressure profile must retain the exact original normal patch and source authority",
            ));
        }
        if self.cells[0].span_m[0] != footprint[0]
            || self.cells.last().unwrap().span_m[1] != footprint[1]
        {
            return Err(MovingPartialSlipError::Invalid(
                "original pressure cells must exactly cover the physical footprint",
            ));
        }
        let load = finite(self.line_load_n_per_m * width_m)?;
        if !close(load, normal.normal_load_n()) {
            return Err(MovingPartialSlipError::Invalid(
                "supplied pressure integral must match the original normal load without rescaling",
            ));
        }
        Ok(())
    }

    // The existing local return operates in reference-patch normalization.
    // Multiplication by the region's original physical area fraction recovers
    // its actual force and extensive work. Pressure and stiffness density are
    // not changed by this representation, and all normal provenance survives.
    pub(super) fn normalized_patch<E>(
        &self,
        pressure_pa: f64,
        footprint_area_m2: f64,
    ) -> Result<NormalPatchView, MovingPartialSlipError<E>> {
        Ok(NormalPatchView::new(
            self.normal.patch_id(),
            self.normal.normal_model_id(),
            self.normal.source_id(),
            self.normal.authority(),
            finite(pressure_pa * footprint_area_m2)?,
            self.normal.semi_axis_longitudinal_m(),
            self.normal.semi_axis_lateral_m(),
            self.normal.pressure_second_moment_m2(),
        )?)
    }
}
