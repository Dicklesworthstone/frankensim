//! Causal ideal-wall images of the original physical stress histories.
//!
//! An affine orthogonal image H with pressure character s carries the full
//! stress tensor s H Q H^T. Evaluate the original Q at virtual receiver H^-1 x
//! instead, and return velocity through (H^-1)^T. Each image retains its own
//! r/c: a transient cannot use the harmonic owner's folded transfer or one
//! effective delay. All geometry and signs are prepared once; observation
//! reuses the native retarded stress kernel and source-history callback.

use super::*;
use crate::curle3d::array::IdealAcousticPlane;
use crate::curle3d::array::wedge::{WedgeImages, multiply, reflection};

pub(super) const ORTHOGONAL_SCOPE: &str = "Prescribed stationary transient symmetric 3D \
stress histories in one homogeneous lossless intersection of one to three perpendicular \
infinite ideal rigid or pressure-release half-spaces; each full tensor image retains its \
own causal delay, pressure and vector velocity before instantaneous energy flux. \
No source-model or bandwidth validation, finite cell size, finite wall impedance, edge \
diffraction, closed-room reverberation, motion, absorption, feedback or measured-SPL qualification.";
pub(super) const WEDGE_SCOPE: &str = "Prescribed stationary transient symmetric 3D stress \
histories in one homogeneous lossless ideal pi/m wedge; every full tensor image retains \
its own causal delay, pressure and vector velocity before instantaneous energy flux. \
Equal rigid/release walls, or mixed walls for even m. No source-model or bandwidth validation, \
finite cell size, finite wall impedance, general diffraction, closed-room reverberation, \
motion, absorption, feedback or measured-SPL qualification.";

type Matrix = [[f64; 3]; 3];
const IDENTITY: Matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

impl PreparedTransientQuadrupoleObserver {
    /// Prepare every exact image of one to three perpendicular ideal walls.
    ///
    /// `maximum_source_paths` includes the direct group and bounds 2^m times
    /// the physical source count. Preparation and storage are O(paths).
    /// Sources must be strictly inside all half-spaces; the receiver may lie
    /// on their walls. Previously attached boundaries are refused.
    pub fn with_orthogonal_planes(
        self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_orthogonal_planes_cancellable(planes, maximum_source_paths, || false)
    }

    /// Bounded preparation; cancellation or refusal publishes no observer.
    pub fn with_orthogonal_planes_cancellable(
        mut self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !(1..=3).contains(&planes.len()) {
            return Err(invalid(
                "transient quadrupole boundaries need one to three planes",
            ));
        }
        for (index, plane) in planes.iter().enumerate() {
            for other in &planes[..index] {
                let dot = (0..3).fold(Dd::ZERO, |sum, axis| {
                    sum + Dd::from_f64(plane.inward_normal()[axis])
                        * Dd::from_f64(other.inward_normal()[axis])
                });
                if dot.to_f64().abs() > 64.0 * f64::EPSILON {
                    return Err(invalid(
                        "transient quadrupole image planes must be perpendicular",
                    ));
                }
            }
        }
        let count = 1usize << planes.len();
        self.admit_images(planes, count, maximum_source_paths, &mut cancelled)?;
        for mask in 1..count {
            let mut receiver = self.observer_m;
            let mut inverse = IDENTITY;
            let mut sign = 1.0;
            for (index, plane) in planes.iter().enumerate() {
                if mask & (1 << index) != 0 {
                    receiver = plane.reflected_position(receiver)?;
                    inverse = multiply(reflection(plane.inward_normal()), inverse);
                    sign *= plane.boundary().sign();
                }
            }
            self.append_image(receiver, inverse, sign, &mut cancelled)?;
        }
        self.reflection_planes.extend_from_slice(planes);
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Prepare the exact 2m dihedral image group of a wedge opening by pi/m.
    ///
    /// Reuses native geometry admission: m is in 2..=128, and mixed rigid/
    /// pressure-release conditions require even m. Source/path ordering is
    /// direct group first, followed by the native dihedral image order.
    pub fn with_wedge(
        self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_wedge_cancellable(planes, order, maximum_source_paths, || false)
    }

    /// Poll for each image/source; all paths retain the original source clock.
    pub fn with_wedge_cancellable(
        mut self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let images = WedgeImages::new(planes, order)?;
        self.admit_images(
            &planes,
            images.path_count(),
            maximum_source_paths,
            &mut cancelled,
        )?;
        images.visit_inverse(|index, inverse, sign| {
            poll(&mut cancelled)?;
            if index != 0 {
                self.append_image(
                    inverse.position(self.observer_m)?,
                    inverse.matrix,
                    sign,
                    &mut cancelled,
                )?;
            }
            Ok(())
        })?;
        self.reflection_planes.extend_from_slice(&planes);
        self.reflection_wedge_order = Some(order);
        poll(&mut cancelled)?;
        Ok(self)
    }

    fn admit_images(
        &mut self,
        planes: &[IdealAcousticPlane],
        image_count: usize,
        maximum_source_paths: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        if !self.reflection_planes.is_empty() {
            return Err(invalid(
                "transient quadrupole boundaries are already prepared",
            ));
        }
        let paths = self
            .source_count()
            .checked_mul(image_count)
            .filter(|paths| *paths <= maximum_source_paths)
            .ok_or_else(|| invalid("transient quadrupole image paths exceed preparation budget"))?;
        for plane in planes {
            poll(cancelled)?;
            if plane.distance(self.observer_m)?.to_f64() < 0.0 {
                return Err(invalid(
                    "transient quadrupole receiver is outside a fluid half-space",
                ));
            }
            for transfer in &self.transfers {
                poll(cancelled)?;
                if plane.distance(transfer.source_m)?.to_f64() <= 0.0 {
                    return Err(invalid(
                        "transient quadrupole sources must lie strictly inside every half-space",
                    ));
                }
            }
        }
        self.transfers
            .try_reserve_exact(paths - self.source_count())
            .map_err(|_| invalid("transient quadrupole image paths exceed available memory"))?;
        self.reflection_planes
            .try_reserve_exact(planes.len())
            .map_err(|_| {
                invalid("transient quadrupole boundary allocation exceeds available memory")
            })?;
        Ok(())
    }

    fn append_image(
        &mut self,
        receiver: [f64; 3],
        inverse: Matrix,
        sign: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        for source in 0..self.source_count() {
            poll(cancelled)?;
            let mut transfer =
                Transfer::new(self.acoustic, receiver, self.transfers[source].source_m)?;
            transfer.pressure_scales = transfer.pressure_scales.map(|value| sign * value);
            transfer.velocity_scales = transfer.velocity_scales.map(|value| sign * value);
            transfer.inverse_image = Some(inverse);
            self.transfers.push(transfer);
        }
        Ok(())
    }
}
