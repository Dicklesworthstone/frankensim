//! Paraxial conjugates and analytic chromatic sensitivities of original media.

use super::{CameraError, DielectricGlass, SequentialLens};

/// Direct paraxial focus and its local derivative with respect to the original
/// vacuum wavelength. All positions are measured from the fixed sensor plane.
/// This is an analytic gradient for chromatic optical design, not a new material
/// fit or a prediction of finite-aperture/diffraction best focus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParaxialFocusSensitivity {
    /// Sensor-side effective focal length, in metres.
    pub focal_length_m: f64,
    /// Real axial object focus measured from the sensor, in metres.
    pub object_focus_m: f64,
    /// Derivative of sensor-side effective focal length, in metres per nm.
    pub focal_length_derivative_m_per_nm: f64,
    /// Derivative of the axial object focus, in metres per nm.
    pub object_focus_derivative_m_per_nm: f64,
}

/// Sensor-to-last-interface map of (height, n*angle), including the first gap.
/// The explicit indices keep physical conjugate distances distinct from the
/// reduced-angle matrix entries when either exterior is an immersion medium.
#[derive(Clone, Copy, Debug)]
pub(super) struct ParaxialMatrix {
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) c: f64,
    pub(super) d: f64,
    pub(super) sensor_index: f64,
    pub(super) scene_index: f64,
    pub(super) last_vertex_z_m: f64,
}

impl ParaxialMatrix {
    fn focus(self) -> Result<(f64, f64), CameraError> {
        let focal = -self.sensor_index / self.c;
        let beyond_front = -self.scene_index * (self.b / self.d);
        let focus = self.last_vertex_z_m + beyond_front;
        if ![focal, beyond_front, focus]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        Ok((focal, focus))
    }
}

struct MatrixDerivative {
    b: f64,
    c: f64,
    d: f64,
    sensor_index: f64,
    scene_index: f64,
}

impl SequentialLens {
    /// Sensor-side effective focal length and axial object focus, in metres,
    /// at an explicitly admitted vacuum wavelength. The object focus is
    /// measured from the sensor, not from the last vertex or a principal plane.
    /// Differences between wavelengths expose longitudinal chromatic aberration
    /// without changing the prescription's reference focus or tracing probe rays.
    ///
    /// This uses the original phase indices of every bulk medium, including
    /// both exteriors, with the reduced-angle coordinates `(height, n*angle)`.
    /// It needs no group index or spectral derivative, so an admitted tabulated
    /// phase-index knot remains usable even when its derivative is undefined.
    /// The complete lens wavelength band, including attached coatings, applies.
    ///
    /// This is the paraxial direct-image limit. Aperture clipping, finite-angle
    /// aberrations, coating phase, ghosts and diffraction do not alter this
    /// geometric focus. Afocal/diverging systems, virtual object focus and
    /// unrepresentable results refuse with `SequentialLensOutsideModel`, as
    /// they do for [`Self::reference_focus`]. Exact rays retain their own model.
    pub fn paraxial_focus_at(&self, wavelength_nm: f64) -> Result<(f64, f64), CameraError> {
        self.paraxial_matrix_at(wavelength_nm)?.focus()
    }

    /// Differentiate the real paraxial object conjugate and sensor-side focal
    /// length at fixed surface positions. Every original bulk law supplies
    /// its phase index and group index, giving dn/dlambda = (n - n_g)/lambda
    /// in vacuum-nanometre units. Both immersion exteriors participate.
    ///
    /// The matrix and its derivative propagate together in one allocation-free
    /// traversal. There is no finite-difference step, probe wavelength, fitted
    /// dispersion or numerical differentiation of traced rays. Use the signed
    /// derivatives in chromatic focus/design objectives; this query does not
    /// move the lens or alter its reference wavelength.
    ///
    /// The complete original lens band applies. A phase-index knot with no
    /// derivative refuses with SequentialLensOutsideModel even when ordinary
    /// paraxial focus at that knot remains valid. At a source-band endpoint,
    /// the original law's admitted one-sided derivative is retained. Group
    /// index supplies only the spectral derivative, never the refractive index
    /// or an additional positive-transit-time restriction.
    ///
    /// These are direct paraxial derivatives: coating phase, finite-aperture
    /// aberrations, ghosts and diffraction remain separate optical models.
    pub fn paraxial_focus_sensitivity_at(
        &self,
        wavelength_nm: f64,
    ) -> Result<ParaxialFocusSensitivity, CameraError> {
        let (matrix, derivative) = self.paraxial_matrix_with_derivative_at::<true>(wavelength_nm)?;
        let (focal_length_m, object_focus_m) = matrix.focus()?;
        // Quotient rules use ratios rather than squaring small C or D. This
        // preserves finite derivatives that an intermediate square can lose.
        let focal_length_derivative_m_per_nm = -derivative.sensor_index / matrix.c
            - focal_length_m * (derivative.c / matrix.c);
        let conjugate_ratio = matrix.b / matrix.d;
        let object_focus_derivative_m_per_nm = -derivative.scene_index * conjugate_ratio
            - matrix.scene_index
                * (derivative.b / matrix.d - conjugate_ratio * (derivative.d / matrix.d));
        if ![
            focal_length_derivative_m_per_nm,
            object_focus_derivative_m_per_nm,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        Ok(ParaxialFocusSensitivity {
            focal_length_m,
            object_focus_m,
            focal_length_derivative_m_per_nm,
            object_focus_derivative_m_per_nm,
        })
    }

    // The mechanical refocus solve needs the nominal map even when this sensor
    // position has a virtual/afocal conjugate. Only focus() imposes that gate.
    pub(super) fn paraxial_matrix_at(
        &self,
        wavelength_nm: f64,
    ) -> Result<ParaxialMatrix, CameraError> {
        Ok(self.paraxial_matrix_with_derivative_at::<false>(wavelength_nm)?.0)
    }

    // Compile-time selection keeps the phase-only path free of group-index
    // evaluation and derivative arithmetic, including at admitted phase kinks.
    fn paraxial_matrix_with_derivative_at<const CHROMATIC: bool>(
        &self,
        wavelength_nm: f64,
    ) -> Result<(ParaxialMatrix, MatrixDerivative), CameraError> {
        let (lo, hi) = self.wavelength_bounds_nm;
        if !wavelength_nm.is_finite()
            || wavelength_nm <= 0.0
            || wavelength_nm < lo
            || wavelength_nm > hi
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let (sensor_index, sensor_derivative) =
            phase_and_derivative::<CHROMATIC>(self.sensor_medium, wavelength_nm)?;
        let (mut a, mut b, mut c, mut d) = (1.0, 0.0, 0.0, 1.0);
        let (mut da, mut db, mut dc, mut dd) = (0.0, 0.0, 0.0, 0.0);
        let (mut z, mut n, mut dn) = (0.0, sensor_index, sensor_derivative);
        for surface in self.surfaces() {
            let travel = (surface.vertex_z_m - z) / n;
            if CHROMATIC {
                let dtravel = -travel * (dn / n);
                da += dtravel * c + travel * dc;
                db += dtravel * d + travel * dd;
            }
            a += travel * c;
            b += travel * d;
            let (next, dnext) =
                phase_and_derivative::<CHROMATIC>(surface.medium_after, wavelength_nm)?;
            let power = (n - next) * surface.curvature_per_m;
            if CHROMATIC {
                let dpower = (dn - dnext) * surface.curvature_per_m;
                dc += dpower * a + power * da;
                dd += dpower * b + power * db;
            }
            c += power * a;
            d += power * b;
            z = surface.vertex_z_m;
            n = next;
            dn = dnext;
        }
        if ![a, b, c, d].iter().all(|value| value.is_finite())
            || (CHROMATIC && ![da, db, dc, dd].iter().all(|value| value.is_finite()))
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        Ok((
            ParaxialMatrix {
                a,
                b,
                c,
                d,
                sensor_index,
                scene_index: n,
                last_vertex_z_m: z,
            },
            MatrixDerivative {
                b: db,
                c: dc,
                d: dd,
                sensor_index: sensor_derivative,
                scene_index: dn,
            },
        ))
    }
}

fn phase_and_derivative<const CHROMATIC: bool>(
    medium: Option<DielectricGlass>,
    wavelength_nm: f64,
) -> Result<(f64, f64), CameraError> {
    let Some(glass) = medium else {
        return Ok((1.0, 0.0));
    };
    let phase = glass
        .ior()
        .eval(wavelength_nm)
        .map_err(|_| CameraError::SequentialLensOutsideModel)?;
    let derivative = if CHROMATIC {
        let group = glass
            .ior()
            .group_index(wavelength_nm)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        let derivative = (phase - group) / wavelength_nm;
        if !derivative.is_finite() {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        derivative
    } else {
        0.0
    };
    Ok((phase, derivative))
}
