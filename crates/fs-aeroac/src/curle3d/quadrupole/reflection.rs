//! Exact ideal images of the original symmetric volume-stress tensor.
//!
//! For an orthogonal affine image H with pressure sign s, a source image
//! carries s H Q H^T. Equivalently evaluate the ORIGINAL tensor at the virtual
//! receiver H^-1 x: p(x)=s p_free(H^-1 x), u(x)=s H u_free(H^-1 x).
//! The latter form retains every mixed stress without a separate six-by-six
//! tensor map. The velocity pullback is the transpose of the inverse map,
//! including rotations in a wedge. All 24 transfer coefficients are folded
//! during preparation; source clocks, coherence and O(N) evaluation survive.

use super::*;
use crate::curle3d::array::wedge::{WedgeImages, multiply, reflection};
use fs_math::dd::Dd;

pub(super) const ORTHOGONAL_SCOPE: &str = "Prescribed symmetric 3D Lighthill stress moments \
[N m] in a homogeneous stationary lossless intersection of one to three perpendicular \
infinite ideal rigid or pressure-release half-spaces; all image pressure/velocity paths \
interfere coherently. No source-size validation, finite wall impedance, edge diffraction, \
enclosed-room reverberation, absorption, motion, source feedback, unresolved SGS sound or \
measured-SPL qualification.";
pub(super) const WEDGE_SCOPE: &str = "Prescribed symmetric 3D Lighthill stress moments [N m] \
in a homogeneous stationary lossless ideal pi/m wedge; all 2m tensor image paths interfere \
coherently. Equal rigid/release walls, or mixed walls for even m. No source-size validation, \
finite wall impedance, general wedge diffraction, enclosed-room reverberation, absorption, \
motion, source feedback, unresolved SGS sound or measured-SPL qualification.";

type Matrix = [[f64; 3]; 3];
type FieldSum = [[ComplexSum; 6]; 4];
const IDENTITY: Matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

impl PreparedQuadrupoleObserver {
    /// Fold the exact 2^m image group of one to three perpendicular ideal walls.
    ///
    /// Sources must lie strictly inside every half-space; receivers may lie
    /// on a wall. The complete group includes the existing direct path and is
    /// bounded by `maximum_source_paths`. Previously attached walls are refused.
    pub fn with_orthogonal_planes(
        self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_orthogonal_planes_cancellable(planes, maximum_source_paths, || false)
    }

    /// Poll before each path/source; refusal publishes no partial observer.
    pub fn with_orthogonal_planes_cancellable(
        mut self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !(1..=3).contains(&planes.len()) {
            return Err(invalid("quadrupole boundaries need one to three planes"));
        }
        for (index, plane) in planes.iter().enumerate() {
            for other in &planes[..index] {
                let dot = (0..3).fold(Dd::ZERO, |sum, axis| {
                    sum + Dd::from_f64(plane.inward_normal()[axis])
                        * Dd::from_f64(other.inward_normal()[axis])
                });
                if dot.to_f64().abs() > 64.0 * f64::EPSILON {
                    return Err(invalid("quadrupole image planes must be perpendicular"));
                }
            }
        }
        let count = 1usize << planes.len();
        let mut sums = self.admit_images(planes, count, maximum_source_paths, &mut cancelled)?;
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
            self.fold_image(receiver, inverse, sign, &mut sums, &mut cancelled)?;
        }
        self.finish_images(planes, None, sums, &mut cancelled)?;
        Ok(self)
    }

    /// Fold the exact 2m image group of two walls with opening pi/m.
    ///
    /// Reuses the original wedge geometry admission (m in 2..=128). Mixed
    /// rigid/release walls require even m; other angles need diffraction and
    /// are refused. Every source remains one coherent physical tensor.
    pub fn with_wedge(
        self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_wedge_cancellable(planes, order, maximum_source_paths, || false)
    }

    /// Bounded preparation using the existing dihedral inverse-image visitor.
    pub fn with_wedge_cancellable(
        mut self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let images = WedgeImages::new(planes, order)?;
        let mut sums = self.admit_images(
            &planes,
            images.path_count(),
            maximum_source_paths,
            &mut cancelled,
        )?;
        images.visit_inverse(|index, inverse, sign| {
            poll(&mut cancelled)?;
            if index != 0 {
                self.fold_image(
                    inverse.position(self.observer_m)?,
                    inverse.matrix,
                    sign,
                    &mut sums,
                    &mut cancelled,
                )?;
            }
            Ok(())
        })?;
        self.finish_images(&planes, Some(order), sums, &mut cancelled)?;
        Ok(self)
    }

    fn admit_images(
        &self,
        planes: &[IdealAcousticPlane],
        image_count: usize,
        maximum_source_paths: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FieldSum>, AeroacError> {
        if !self.reflection_planes.is_empty() {
            return Err(invalid("quadrupole boundaries are already prepared"));
        }
        if self
            .source_count()
            .checked_mul(image_count)
            .is_none_or(|n| n > maximum_source_paths)
        {
            return Err(invalid("quadrupole image paths exceed preparation budget"));
        }
        for plane in planes {
            poll(cancelled)?;
            if plane.distance(self.observer_m)?.to_f64() < 0.0 {
                return Err(invalid("quadrupole receiver is outside a fluid half-space"));
            }
            for transfer in &self.transfers {
                poll(cancelled)?;
                if plane.distance(transfer.source_m)?.to_f64() <= 0.0 {
                    return Err(invalid(
                        "quadrupole sources must lie strictly inside every half-space",
                    ));
                }
            }
        }
        let mut sums = Vec::new();
        sums.try_reserve_exact(self.source_count())
            .map_err(|_| invalid("quadrupole image sums exceed available memory"))?;
        for transfer in &self.transfers {
            poll(cancelled)?;
            let mut sum = [[ComplexSum::default(); 6]; 4];
            for (out, value) in sum
                .iter_mut()
                .flatten()
                .zip(transfer.field.iter().flatten())
            {
                out.add(*value)?;
            }
            sums.push(sum);
        }
        Ok(sums)
    }

    fn fold_image(
        &self,
        receiver: [f64; 3],
        inverse: Matrix,
        sign: f64,
        sums: &mut [FieldSum],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        for (transfer, sum) in self.transfers.iter().zip(sums) {
            poll(cancelled)?;
            let displacement =
                core::array::from_fn(|axis| receiver[axis] - transfer.source_m[axis]);
            let field = kernel::transfer(
                self.acoustic.density_kg_m3(),
                self.acoustic.sound_speed_m_s(),
                self.frequency_hz,
                displacement,
            )
            .ok_or_else(|| {
                invalid("quadrupole image separation or transfer is not representable")
            })?;
            for component in 0..6 {
                sum[0][component].add(field[0][component].scale(sign))?;
                for axis in 0..3 {
                    // grad_x p(H^-1 x) = (H^-1)^T grad p. For a wedge,
                    // using the inverse instead of its transpose is WRONG.
                    let mut velocity = ComplexSum::default();
                    for image_axis in 0..3 {
                        velocity.add(
                            field[image_axis + 1][component].scale(inverse[image_axis][axis]),
                        )?;
                    }
                    sum[axis + 1][component].add(velocity.value().scale(sign))?;
                }
            }
        }
        Ok(())
    }

    fn finish_images(
        &mut self,
        planes: &[IdealAcousticPlane],
        order: Option<usize>,
        sums: Vec<FieldSum>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        self.reflection_planes
            .try_reserve_exact(planes.len())
            .map_err(|_| invalid("quadrupole boundary allocation exceeds available memory"))?;
        for (transfer, sum) in self.transfers.iter_mut().zip(sums) {
            poll(cancelled)?;
            transfer.field = sum.map(|row| row.map(ComplexSum::value));
            if transfer.field.iter().flatten().any(|value| !finite(*value)) {
                return Err(nonfinite("coherent quadrupole image transfer"));
            }
        }
        self.reflection_planes.extend_from_slice(planes);
        self.reflection_wedge_order = order;
        poll(cancelled)
    }
}
