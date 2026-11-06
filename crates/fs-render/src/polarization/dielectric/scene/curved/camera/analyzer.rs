//! Passive polarization measurements through the original reciprocal camera.
//!
//! The physical detector row a acts on light arriving along -ray.dir. Its
//! horizontal axis is the supplied world axis projected onto that ray's
//! transverse plane. In opposite travel frames sharing this axis the native
//! reciprocity convention is A=diag(1,1,-1,1), so the camera seed is A*a.
//! With M_reverse=A*M_forward^T*A, unpolarized emission gives
//! e0^T*M_reverse*A*a = a^T*M_forward*e0. V is NOT negated, and the seed's
//! intensity retains the physical analyzer throughput instead of being one.
//!
//! This is an ideal local analyzer at the geometric camera launch plane, with
//! a wavelength-independent declared row. It supplies no finite-thickness or
//! angle-dependent polarizer, polarized source law, bulk birefringence, or
//! propagation through a resolved camera lens. The existing ray geometry,
//! complete Stokes Fresnel response, Beer path and radiance index conversion
//! remain authoritative. All-zero rows refuse because the original traversal
//! cannot distinguish extinguished power from floating-point underflow.

use super::{CameraRadianceError, CameraRadianceSample, CurvedCameraScene};
use crate::charts::Ray;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{CurvedSceneError, CurvedTraceSettings};
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, PolarizedBeam};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

#[cfg(feature = "volumes")]
mod molecular;

/// Passive intensity row in the actual arriving-light detector frame.
///
/// For every physical Stokes input its response lies between zero and the
/// original intensity: |(aQ,aU,aV)| <= min(a0,1-a0), with 0<a0<=1.
/// The axis is projected independently for each physical camera ray. Reversing
/// its sign is a pi frame rotation and does not change the measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CameraAnalyzer {
    row: Stokes,
    reference_axis: Vec3,
}

impl CameraAnalyzer {
    /// Retain the declared passive row [a0,aQ,aU,aV] without normalization or
    /// fitting. A0 is the transmission of unpolarized light, not peak gain.
    /// The original frame owner refuses an axis parallel to a queried ray.
    pub fn try_new(reference_axis: Vec3, row: Stokes) -> Result<Self, BeamError> {
        if ![reference_axis.x, reference_axis.y, reference_axis.z]
            .into_iter()
            .all(f64::is_finite)
            || reference_axis
                .x
                .abs()
                .max(reference_axis.y.abs())
                .max(reference_axis.z.abs())
                == 0.0
        {
            return Err(BeamError::InvalidFrame);
        }
        if ![row.i, row.q, row.u, row.v].into_iter().all(f64::is_finite)
            || row.i <= 0.0
            || row.i > 1.0
        {
            return Err(BeamError::InvalidStokes);
        }
        let degree = det::hypot(det::hypot(row.q, row.u), row.v);
        if !degree.is_finite() || degree > row.i.min(1.0 - row.i) {
            return Err(BeamError::InvalidStokes);
        }
        Ok(Self {
            row,
            reference_axis,
        })
    }

    /// Ideal linear analyzer at the stated angle in the incoming detector
    /// frame. Its unpolarized transmission is one half. Generated trigonometric
    /// roundoff is directed inward to the passive cone; authored rows are never
    /// repaired by [`Self::try_new`].
    pub fn try_linear(reference_axis: Vec3, angle_rad: f64) -> Result<Self, BeamError> {
        if !angle_rad.is_finite() || !(2.0 * angle_rad).is_finite() {
            return Err(BeamError::InvalidStokes);
        }
        let (s, c) = det::sin_cos(2.0 * angle_rad);
        let mut row = Stokes {
            i: 0.5,
            q: 0.5 * c,
            u: 0.5 * s,
            v: 0.0,
        };
        // A unit-circle evaluation can lie one ulp outside the passive cone.
        // Move only the larger generated component toward zero; its orientation
        // changes by roundoff rather than admitting a negative response.
        for _ in 0..4 {
            if det::hypot(row.q, row.u) <= 0.5 {
                break;
            }
            let value = if row.q.abs() >= row.u.abs() {
                &mut row.q
            } else {
                &mut row.u
            };
            *value = value.abs().next_down().copysign(*value);
        }
        Self::try_new(reference_axis, row)
    }

    /// Ideal circular analyzer selecting the given sign of incoming Stokes V.
    /// Reciprocal launch preserves that V sign, unlike a geometric reversal.
    pub fn try_circular(reference_axis: Vec3, positive_v: bool) -> Result<Self, BeamError> {
        Self::try_new(
            reference_axis,
            Stokes {
                i: 0.5,
                q: 0.0,
                u: 0.0,
                v: if positive_v { 0.5 } else { -0.5 },
            },
        )
    }

    pub const fn row(&self) -> Stokes {
        self.row
    }
    pub const fn reference_axis(&self) -> Vec3 {
        self.reference_axis
    }

    fn launch(
        self,
        camera: &CurvedCameraScene<'_, '_>,
        ray: Ray,
        wavelength_nm: f64,
        medium: usize,
    ) -> Result<PolarizedBeam, CameraRadianceError> {
        if self.row.q == 0.0 && self.row.u == 0.0 && self.row.v == 0.0 {
            // An isotropic response has no distinguished detector axis.
            // Preserve the original basis and every bit for an identity row.
            let original = camera.launch(ray, wavelength_nm, medium)?;
            if self.row.i == 1.0 {
                return Ok(original);
            }
            return PolarizedBeam::try_new(
                original.ray().origin,
                original.frame(),
                wavelength_nm,
                original.medium(),
                self.row,
            )
            .map_err(CurvedSceneError::from)
            .map_err(Into::into);
        }
        let frame = PolarizationFrame::try_new(ray.dir, self.reference_axis)
            .map_err(CurvedSceneError::from)?;
        let seed = Stokes {
            u: -self.row.u,
            ..self.row
        };
        Ok(camera
            .scene
            .launch_beam(ray.origin, frame, wavelength_nm, medium, seed)?)
    }
}

impl CurvedCameraScene<'_, '_> {
    /// Measure the original unpolarized source scene through a declared passive
    /// analyzer at this camera ray's launch. Native counter streams, per-path
    /// budgets, thermal/surface emission and endpoint n² conversion are shared
    /// with [`Self::sample_ray`]; no post-hoc scalar attenuation approximates
    /// the interaction of the analyzer with intervening polarizing interfaces.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// The same physical analyzer with original [facet_u,facet_v,branch]
    /// variates. Native scalar proposal probabilities remain independent of
    /// the analyzer; every Stokes response retains their division exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples(cx, beam, initial_medium, settings, samples)
    }
}

#[cfg(test)]
mod tests;
