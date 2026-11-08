//! Ideal images for uniform flow parallel to stationary reflecting planes.
//!
//! Reflection leaves both |x-y| and U.(x-y) invariant when U.n=0. Thus the
//! original convected Green function commutes with each admitted reflection.
//! Its even extension has zero normal particle velocity; its odd extension
//! has zero pressure. This symmetry acts on the full potential, including
//! convection and induction, not just a pressure coefficient or travel delay.

use super::*;
use crate::curle3d::array::{IdealAcousticPlane, wedge::WedgeImages};

pub(super) const CONVECTED_IDEAL_BOUNDARY_SCOPE: &str = "Prescribed compact volume injection \
at the uniform mean velocity in a homogeneous subsonic lossless fluid with one to three \
perpendicular infinite ideal rigid or pressure-release planes; flow is tangent to every \
wall and all images remain coherent. No shear/wall layers, finite impedance, nonuniform \
flow, diffraction, closed-room reverberation, source feedback or measured-SPL qualification.";

pub(super) const CONVECTED_WEDGE_SCOPE: &str = "Prescribed compact volume injection \
at the uniform mean velocity in a homogeneous subsonic lossless fluid within an infinite \
ideal pi/m wedge; flow is parallel to the edge and all 2m images remain coherent. Equal \
rigid/release walls, or mixed walls for even m. No shear/wall layers, finite impedance, \
nonuniform flow, general wedge diffraction, closed-room reverberation, source feedback \
or measured-SPL qualification.";

impl PreparedConvectedMonopoleObserver {
    /// Add the complete ideal image group for one to three perpendicular walls.
    ///
    /// The original mean flow must be tangent to EVERY plane. A normal flow
    /// component would violate the stationary boundary and the image symmetry;
    /// it is refused, never projected away. Perpendicularity and tangency are
    /// admitted within 64 machine epsilons in normalized directions, without
    /// changing the supplied geometry or velocity. Sources lie strictly inside
    /// every fluid half-space; receivers may lie on the boundaries.
    ///
    /// Preparation evaluates all 2^m paths per physical source using the same
    /// material, flow, frequency and clearance limits. `maximum_source_paths`
    /// bounds this complete work, including the direct paths. Compensated sums
    /// fold four p/u coefficients per source, so evaluation and statistical
    /// factors remain allocation-free O(N) with the original source count.
    /// Direct-path geometry/delay accessors keep their original meaning.
    pub fn with_orthogonal_planes(
        self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_orthogonal_planes_cancellable(planes, maximum_source_paths, || false)
    }

    /// Poll per source/path; cancellation or refusal publishes no partial owner.
    pub fn with_orthogonal_planes_cancellable(
        mut self,
        planes: &[IdealAcousticPlane],
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !self.reflection_planes.is_empty() || !(1..=3).contains(&planes.len()) {
            return Err(invalid(
                "convected boundaries need one original set of one to three planes",
            ));
        }
        let image_count = 1usize << planes.len();
        if self
            .source_count()
            .checked_mul(image_count)
            .is_none_or(|n| n > maximum_source_paths)
        {
            return Err(invalid("convected image paths exceed preparation budget"));
        }
        for (index, plane) in planes.iter().enumerate() {
            poll(&mut cancelled)?;
            let normal = plane.inward_normal();
            if planes[..index]
                .iter()
                .any(|other| dot(normal, other.inward_normal()).abs() > 64.0 * f64::EPSILON)
            {
                return Err(invalid("convected image planes must be perpendicular"));
            }
        }
        self.admit_reflecting_half_spaces(planes, &mut cancelled)?;
        self.reflection_planes
            .try_reserve_exact(planes.len())
            .map_err(|_| invalid("convected boundary allocation exceeds available memory"))?;
        let mut sums = Vec::new();
        sums.try_reserve_exact(self.source_count())
            .map_err(|_| invalid("convected image sum allocation exceeds available memory"))?;
        for path in &self.paths {
            poll(&mut cancelled)?;
            let mut values = [ComplexSum::default(); 4];
            for (sum, field) in values.iter_mut().zip(path.field) {
                sum.add(field)?;
            }
            sums.push(values);
        }
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(self.source_count())
            .map_err(|_| invalid("convected image geometry exceeds available memory"))?;
        for mask in 1..image_count {
            positions.clear();
            let mut sign = 1.0;
            for (index, plane) in planes.iter().enumerate() {
                if mask & (1 << index) != 0 {
                    sign *= plane.boundary().sign();
                }
            }
            for path in &self.paths {
                poll(&mut cancelled)?;
                let mut position = path.source_m;
                for (index, plane) in planes.iter().enumerate() {
                    if mask & (1 << index) != 0 {
                        position = plane.image(position)?;
                    }
                }
                positions.push(position);
            }
            self.fold_image(&positions, sign, &mut sums, &mut cancelled)?;
        }
        for (path, values) in self.paths.iter_mut().zip(sums) {
            poll(&mut cancelled)?;
            path.field = values.map(ComplexSum::value);
            if path.field.iter().any(|&value| !finite(value)) {
                return Err(nonfinite("coherent convected image transfer"));
            }
        }
        self.reflection_planes.extend_from_slice(planes);
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Add the complete 2m image group of an ideal wedge opening by pi/m.
    ///
    /// `order` is m in 2..=128. The original flow must be tangent to BOTH
    /// walls, hence parallel to their edge (or zero). Equal rigid/release
    /// conditions admit every supported m; mixed walls require even m.
    /// The original plane geometry and flow are never snapped or projected.
    /// Sources must lie strictly inside the wedge; receivers may lie on it.
    ///
    /// `maximum_source_paths` bounds 2m times the physical source count,
    /// including direct paths. Preparation uses O(mN) work and O(N) scratch;
    /// four compensated p/u coefficients per original source preserve O(N),
    /// allocation-free field evaluation and the original source statistics,
    /// clock, material and direct-path delay metadata.
    pub fn with_wedge(
        self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_wedge_cancellable(planes, order, maximum_source_paths, || false)
    }

    /// Poll per source/image; cancellation publishes no partial reflection bank.
    pub fn with_wedge_cancellable(
        mut self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !self.reflection_planes.is_empty() {
            return Err(invalid(
                "convected acoustic boundaries are already prepared",
            ));
        }
        let images = WedgeImages::new(planes, order)?;
        if self
            .source_count()
            .checked_mul(images.path_count())
            .is_none_or(|n| n > maximum_source_paths)
        {
            return Err(invalid("convected wedge paths exceed preparation budget"));
        }
        self.admit_reflecting_half_spaces(&planes, &mut cancelled)?;
        self.reflection_planes
            .try_reserve_exact(2)
            .map_err(|_| invalid("convected boundary allocation exceeds available memory"))?;
        let mut sums = Vec::new();
        sums.try_reserve_exact(self.source_count())
            .map_err(|_| invalid("convected image sum allocation exceeds available memory"))?;
        for path in &self.paths {
            poll(&mut cancelled)?;
            let mut values = [ComplexSum::default(); 4];
            for (sum, field) in values.iter_mut().zip(path.field) {
                sum.add(field)?;
            }
            sums.push(values);
        }
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(self.source_count())
            .map_err(|_| invalid("convected image geometry exceeds available memory"))?;
        // The finite dihedral group is closed under inversion and its real
        // +/-1 boundary character obeys chi(H^-1)=chi(H). Thus the existing
        // inverse-map visitor also enumerates the complete SOURCE image group.
        images.visit_inverse(|index, image, sign| {
            poll(&mut cancelled)?;
            if index == 0 {
                return Ok(()); // direct coefficients are already included
            }
            positions.clear();
            for path in &self.paths {
                poll(&mut cancelled)?;
                positions.push(image.position(path.source_m)?);
            }
            self.fold_image(&positions, sign, &mut sums, &mut cancelled)
        })?;
        for (path, values) in self.paths.iter_mut().zip(sums) {
            poll(&mut cancelled)?;
            path.field = values.map(ComplexSum::value);
            if path.field.iter().any(|&value| !finite(value)) {
                return Err(nonfinite("coherent convected wedge transfer"));
            }
        }
        self.reflection_planes.extend_from_slice(&planes);
        self.reflection_wedge_order = Some(order);
        poll(&mut cancelled)?;
        Ok(self)
    }

    fn admit_reflecting_half_spaces(
        &self,
        planes: &[IdealAcousticPlane],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        // Scale before normalization so finite, very large or tiny velocities
        // retain the SAME directional tangency test without squaring overflow.
        let scale = self
            .flow_velocity_m_s
            .iter()
            .fold(0.0_f64, |s, v| s.max(v.abs()));
        let direction = if scale == 0.0 {
            [0.0; 3]
        } else {
            let scaled = self.flow_velocity_m_s.map(|v| v / scale);
            let length = norm(scaled);
            scaled.map(|v| v / length)
        };
        for plane in planes {
            poll(cancelled)?;
            if dot(direction, plane.inward_normal()).abs() > 64.0 * f64::EPSILON {
                return Err(invalid(
                    "uniform flow must be tangent to every reflecting plane",
                ));
            }
            if plane.distance(self.observer_m)?.to_f64() < 0.0 {
                return Err(invalid(
                    "convected receiver is outside an acoustic half-space",
                ));
            }
            for source in &self.paths {
                poll(cancelled)?;
                if plane.distance(source.source_m)?.to_f64() <= 0.0 {
                    return Err(invalid(
                        "convected sources must lie strictly inside every half-space",
                    ));
                }
            }
        }
        Ok(())
    }

    fn fold_image(
        &self,
        positions: &[[f64; 3]],
        sign: f64,
        sums: &mut [[ComplexSum; 4]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        let image = Self::new_cancellable(
            self.acoustic,
            self.flow_velocity_m_s,
            self.frequency_hz,
            self.observer_m,
            positions,
            self.source_count(),
            self.limits,
            &mut *cancelled,
        )?;
        for (path, values) in image.paths.iter().zip(sums) {
            poll(cancelled)?;
            // The image potential's gradient already uses the observer's
            // original axes. Do not reflect particle velocity a second time.
            for (sum, field) in values.iter_mut().zip(path.field) {
                sum.add(field.scale(sign))?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
