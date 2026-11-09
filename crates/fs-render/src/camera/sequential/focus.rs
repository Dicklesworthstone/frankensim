//! Explicit paraxial focusing by translating the original lens prescription.
//!
//! The sensor and requested object plane stay fixed. Every authored surface
//! receives the same axial translation; glass thicknesses, curvatures, aspheres,
//! stops, coatings and ghost policy retain their original owners. The caller
//! supplies permitted mechanical travel and an object-focus residual tolerance.
//! This does not choose a scene target or optimize finite-aperture image quality.

use super::paraxial::ParaxialMatrix;
use super::{CameraError, SequentialLens};
use fs_math::dd::Dd;

/// An admitted axial repositioning of the original sequential prescription.
#[derive(Clone, Debug, PartialEq)]
pub struct SequentialLensFocus {
    lens: SequentialLens,
    axial_shift_m: f64,
    target_distance_m: f64,
    achieved_distance_m: f64,
}

impl SequentialLensFocus {
    /// The translated prescription, focused at its new reference wavelength.
    pub const fn lens(&self) -> &SequentialLens {
        &self.lens
    }

    /// Move the prescription into the original physical-camera owner.
    pub fn into_lens(self) -> SequentialLens {
        self.lens
    }

    /// Common surface displacement; positive moves the lens toward the object.
    pub const fn axial_shift_m(&self) -> f64 {
        self.axial_shift_m
    }

    /// Requested axial object-plane distance from the fixed sensor, in metres.
    pub const fn target_distance_m(&self) -> f64 {
        self.target_distance_m
    }

    /// Focus reread from the translated prescription's original paraxial model.
    pub const fn achieved_distance_m(&self) -> f64 {
        self.achieved_distance_m
    }
}

impl SequentialLens {
    /// Reposition the complete lens to image a prescribed finite axial plane.
    ///
    /// Distances are SI metres, measured from the unchanged sensor at z=0;
    /// wavelength is vacuum nm. `shift_range_m` is a finite inclusive mechanical
    /// travel interval. The accepted absolute residual in the paraxial object
    /// distance is the caller's finite nonnegative `focus_tolerance_m`.
    ///
    /// The two conjugate positions are solved analytically from the original
    /// reduced-angle matrix. Among positions inside the travel interval and
    /// wholly between sensor and object, choose the smallest absolute shift;
    /// an exact tie chooses the lower signed shift. Source support, positive
    /// focal power, ordered cap bounds and the requested residual are checked
    /// before returning a complete prescription. The input remains unchanged.
    ///
    /// The returned lens uses `wavelength_nm` as its reference, so
    /// `PhysicalCamera::with_sequential_lens` consumes that same focus directly.
    /// Every other material, surface and sampling setting is retained. This
    /// moves the whole prescription, not a selected internal focusing group.
    /// Finite-aperture aberrations, coating phase and diffraction retain their
    /// existing models; a small paraxial residual is not a spot-size guarantee.
    pub fn refocused_for_object(
        &self,
        wavelength_nm: f64,
        object_distance_m: f64,
        shift_range_m: [f64; 2],
        focus_tolerance_m: f64,
    ) -> Result<SequentialLensFocus, CameraError> {
        if !object_distance_m.is_finite() || object_distance_m <= 0.0 {
            return Err(CameraError::InvalidFocusDistance);
        }
        if !shift_range_m.into_iter().all(f64::is_finite)
            || shift_range_m[0] > shift_range_m[1]
            || !focus_tolerance_m.is_finite()
            || focus_tolerance_m < 0.0
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let matrix = self.paraxial_matrix_at(wavelength_nm)?;
        if matrix.c >= 0.0 {
            return Err(CameraError::SequentialLensOutsideModel);
        }

        // Preserve an exactly satisfied geometry, including a zero-width
        // allowed travel interval. Changing the reference wavelength remains
        // explicit even when no mechanical motion is required.
        if shift_range_m[0] <= 0.0 && shift_range_m[1] >= 0.0 {
            if let Ok((_, focus)) = self.paraxial_focus_at(wavelength_nm) {
                if focus == object_distance_m
                    && self.surfaces().iter().all(|surface| {
                        let bounds = surface.bounds();
                        bounds[0] > 0.0 && bounds[1] < object_distance_m
                    })
                {
                    let mut lens = self.clone();
                    lens.reference_wavelength_nm = wavelength_nm;
                    return Ok(SequentialLensFocus {
                        lens,
                        axial_shift_m: 0.0,
                        target_distance_m: object_distance_m,
                        achieved_distance_m: focus,
                    });
                }
            }
        }

        let mut shifts = focus_shifts(matrix, object_distance_m)
            .ok_or(CameraError::SequentialLensOutsideModel)?;
        shifts.sort_by(|left, right| {
            left.abs()
                .total_cmp(&right.abs())
                .then_with(|| left.total_cmp(right))
        });
        'candidate: for shift in shifts {
            if !shift.is_finite() || shift < shift_range_m[0] || shift > shift_range_m[1] {
                continue;
            }
            let mut surfaces = Vec::with_capacity(self.surfaces().len());
            let mut previous_high = 0.0;
            for original in self.surfaces() {
                let mut surface = *original;
                surface.vertex_z_m += shift;
                let [low, high] = surface.bounds();
                if ![surface.vertex_z_m, low, high].into_iter().all(f64::is_finite)
                    || low <= previous_high
                    || high >= object_distance_m
                {
                    continue 'candidate;
                }
                previous_high = high;
                surfaces.push(surface);
            }
            let mut lens = self.clone();
            lens.surfaces = surfaces.into();
            lens.reference_wavelength_nm = wavelength_nm;
            // Reevaluate the actual rounded surface coordinates, rather than
            // reporting the requested distance or trusting the quadratic alone.
            let Ok((_, achieved)) = lens.paraxial_focus_at(wavelength_nm) else {
                continue;
            };
            if (achieved - object_distance_m).abs() > focus_tolerance_m {
                continue;
            }
            return Ok(SequentialLensFocus {
                lens,
                axial_shift_m: if shift == 0.0 { 0.0 } else { shift },
                target_distance_m: object_distance_m,
                achieved_distance_m: achieved,
            });
        }
        Err(CameraError::SequentialLensOutsideModel)
    }
}

// M=[a,b;c,d] maps (height,n*angle) from sensor to the original last vertex.
// Translation by x gives M*T(x/ns), then the object gap is S-z_front-x.
// Its imaging condition B_total=0 is
// -c*x² + [nt*a+c*(S-z_front)-ns*d]*x
//       + ns*[nt*b+(S-z_front)*d] = 0.
// The existing nominal matrix, not a thin-lens approximation, supplies M.
// Retained products and stable roots protect the two conjugate branches near
// their coalescence. This is numerical evaluation, not an interval certificate.
fn focus_shifts(matrix: ParaxialMatrix, object_distance_m: f64) -> Option<[f64; 2]> {
    let ns = Dd::from_f64(matrix.sensor_index);
    let nt = Dd::from_f64(matrix.scene_index);
    let gap = Dd::from_f64(object_distance_m) - Dd::from_f64(matrix.last_vertex_z_m);
    let mut a = -Dd::from_f64(matrix.c);
    let mut b = nt * Dd::from_f64(matrix.a) + Dd::from_f64(matrix.c) * gap
        - ns * Dd::from_f64(matrix.d);
    let mut c = ns * (nt * Dd::from_f64(matrix.b) + gap * Dd::from_f64(matrix.d));
    if [a, b, c]
        .into_iter()
        .any(|value| !value.hi.is_finite() || !value.lo.is_finite())
    {
        return None;
    }
    let scale = a.hi.abs().max(b.hi.abs()).max(c.hi.abs());
    if scale == 0.0 {
        return None;
    }
    let scale = Dd::from_f64(scale);
    a = a / scale;
    b = b / scale;
    c = c / scale;
    if a == Dd::ZERO {
        return None;
    }
    let discriminant = b * b - Dd::from_f64(4.0) * a * c;
    if !discriminant.hi.is_finite()
        || !discriminant.lo.is_finite()
        || discriminant.lt(Dd::ZERO)
    {
        return None;
    }
    let root = discriminant.sqrt();
    let q = Dd::from_f64(-0.5) * (b + if b.lt(Dd::ZERO) { -root } else { root });
    if q == Dd::ZERO {
        let root = (-b / (Dd::from_f64(2.0) * a)).to_f64();
        return Some([root, root]);
    }
    Some([(q / a).to_f64(), (c / q).to_f64()])
}
