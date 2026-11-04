//! Finite dihedral image groups for an ideal acoustic wedge of opening pi/m.
//!
//! With the two physical wall reflections S1 and S2, the group consists of
//! R^k and R^k S1, where R = S2 S1 and 0 <= k < m. Its pressure character is
//! (s2 s1)^k and s1 (s2 s1)^k, respectively. The relation R^m = I requires
//! (s2 s1)^m = +1: unlike equal wall conditions, mixed walls need even m.
//! Every image carries the FULL Cartesian dipole transform of its group
//! element. Summing just reflected monopoles cannot enforce the walls for a
//! prescribed force, especially in its tangential reactive near field.

use super::*;
use fs_math::dd::Dd;

pub(super) const LOSSLESS_SCOPE: &str = "Prescribed compact 3D Q/F sources in a homogeneous \
stationary lossless ideal pi/m wedge; all 2m dihedral image paths interfere coherently. \
Equal rigid/release walls, or mixed walls for even m. No finite wall impedance, \
general wedge diffraction, enclosed-room reverberation, motion, source feedback or measured-SPL qualification.";
pub(super) const MOLECULAR_SCOPE: &str = "Prescribed compact 3D Q/F sources with material-owned \
molecular bulk relaxation and inviscid momentum in an ideal pi/m wedge; all 2m dihedral \
image paths interfere coherently. Equal rigid/release walls, or mixed walls for even m. \
No wall boundary layers, finite wall impedance, general wedge diffraction, enclosed-room \
reverberation, motion, source feedback or measured-SPL qualification.";

type Matrix = [[f64; 3]; 3];
type FieldSum = [[ComplexSum; 4]; 4];
const IDENTITY: Matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const MAXIMUM_ORDER: usize = 128;

impl PreparedCompactSourceObserver {
    /// Add the exact finite image group of two ideal walls opening by pi/m.
    ///
    /// `order` is m, in 2..=128. Inward normals must subtend pi - pi/m,
    /// within floating-point roundoff; geometry is never snapped. Both walls
    /// may be rigid or pressure release. Mixed conditions require even m;
    /// an odd mixed wedge needs a diffracted field beyond a finite image sum.
    /// Sources must be strictly inside both walls; receivers may be on them.
    ///
    /// `maximum_source_paths` bounds 2m times the physical source count,
    /// including direct paths. Preparation takes O(mN) work and O(N) scratch;
    /// it folds all sixteen pressure/velocity coefficients into each original
    /// source, preserving O(N), allocation-free evaluation and source clock,
    /// material and coherence semantics. Previously attached walls are refused.
    pub fn with_wedge(
        self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
    ) -> Result<Self, AeroacError> {
        self.with_wedge_cancellable(planes, order, maximum_source_paths, || false)
    }

    /// Poll during geometry admission, every image path and coefficient fold.
    pub fn with_wedge_cancellable(
        mut self,
        planes: [IdealAcousticPlane; 2],
        order: usize,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        plane::poll(&mut cancelled)?;
        if !self.reflection_planes.is_empty() {
            return Err(invalid("compact acoustic boundaries are already prepared"));
        }
        if !(2..=MAXIMUM_ORDER).contains(&order) {
            return Err(invalid("ideal acoustic wedge order must be in 2..=128"));
        }
        let step_sign = planes[0].boundary().sign() * planes[1].boundary().sign();
        if step_sign < 0.0 && order % 2 != 0 {
            return Err(invalid("mixed ideal wedge boundaries require even order"));
        }
        if self
            .source_count()
            .checked_mul(2 * order)
            .is_none_or(|count| count > maximum_source_paths)
        {
            return Err(invalid("wedge source paths exceed preparation budget"));
        }
        check_opening(planes, order)?;
        for wall in planes {
            plane::poll(&mut cancelled)?;
            if wall.distance(self.observer_m)?.to_f64() < 0.0 {
                return Err(invalid("compact receiver is outside the acoustic wedge"));
            }
            for source in &self.transfers {
                plane::poll(&mut cancelled)?;
                if wall.distance(source.source_m)?.to_f64() <= 0.0 {
                    return Err(invalid("compact sources must be strictly inside the wedge"));
                }
            }
        }
        self.reflection_planes
            .try_reserve_exact(2)
            .map_err(|_| invalid("acoustic boundary allocation exceeds available memory"))?;
        let mut sums = Vec::new();
        let mut rotations = Vec::new();
        let mut mirrors = Vec::new();
        sums.try_reserve_exact(self.source_count())
            .map_err(|_| invalid("acoustic image sum allocation exceeds available memory"))?;
        for positions in [&mut rotations, &mut mirrors] {
            positions
                .try_reserve_exact(self.source_count())
                .map_err(|_| invalid("acoustic image geometry exceeds available memory"))?;
        }
        for transfer in &self.transfers {
            plane::poll(&mut cancelled)?;
            let mut sum = [[ComplexSum::default(); 4]; 4];
            add_field(&mut sum, transfer.field)?;
            sums.push(sum);
            rotations.push(transfer.source_m);
            mirrors.push(planes[0].image(transfer.source_m)?);
        }
        let first_reflection = reflection(planes[0].inward_normal());
        let step = multiply(reflection(planes[1].inward_normal()), first_reflection);
        let mut rotation_matrix = IDENTITY;
        let mut mirror_matrix = first_reflection;
        let mut rotation_sign = 1.0;
        for index in 0..order {
            plane::poll(&mut cancelled)?;
            if index > 0 {
                fold_image(
                    &self,
                    &rotations,
                    planes,
                    rotation_matrix,
                    rotation_sign,
                    &mut sums,
                    &mut cancelled,
                )?;
            }
            fold_image(
                &self,
                &mirrors,
                planes,
                mirror_matrix,
                rotation_sign * planes[0].boundary().sign(),
                &mut sums,
                &mut cancelled,
            )?;
            if index + 1 < order {
                for positions in [&mut rotations, &mut mirrors] {
                    for position in positions {
                        plane::poll(&mut cancelled)?;
                        *position = planes[1]
                            .reflected_position(planes[0].reflected_position(*position)?)?;
                    }
                }
                rotation_matrix = multiply(step, rotation_matrix);
                mirror_matrix = multiply(step, mirror_matrix);
                rotation_sign *= step_sign;
            }
        }
        for (transfer, sum) in self.transfers.iter_mut().zip(sums) {
            plane::poll(&mut cancelled)?;
            transfer.field = sum.map(|row| row.map(ComplexSum::value));
            if transfer.field.iter().flatten().any(|value| !finite(*value)) {
                return Err(nonfinite("coherent acoustic wedge transfer"));
            }
        }
        self.reflection_planes.extend_from_slice(&planes);
        self.reflection_wedge_order = Some(order);
        plane::poll(&mut cancelled)?;
        Ok(self)
    }
}

fn check_opening(planes: [IdealAcousticPlane; 2], order: usize) -> Result<(), AeroacError> {
    let a = planes[0].inward_normal();
    let b = planes[1].inward_normal();
    let dot = (0..3).fold(Dd::ZERO, |sum, axis| {
        sum + Dd::from_f64(a[axis]) * Dd::from_f64(b[axis])
    });
    let cross: [f64; 3] = core::array::from_fn(|axis| {
        let j = (axis + 1) % 3;
        let k = (axis + 2) % 3;
        (Dd::from_f64(a[j]) * Dd::from_f64(b[k]) - Dd::from_f64(a[k]) * Dd::from_f64(b[j])).to_f64()
    });
    let sine = det::hypot(det::hypot(cross[0], cross[1]), cross[2]);
    let (expected_sine, expected_cosine) = det::sin_cos(core::f64::consts::PI / order as f64);
    let tolerance = 64.0 * f64::EPSILON;
    if (dot.to_f64() + expected_cosine).abs() > tolerance
        || (sine - expected_sine).abs() > tolerance * expected_sine
    {
        return Err(invalid(
            "acoustic wedge planes must open by exactly pi/order",
        ));
    }
    Ok(())
}

fn reflection(normal: [f64; 3]) -> Matrix {
    core::array::from_fn(|row| {
        core::array::from_fn(|col| {
            let identity = if row == col { 1.0 } else { 0.0 };
            (Dd::from_f64(identity) - Dd::from_f64(2.0 * normal[row]) * Dd::from_f64(normal[col]))
                .to_f64()
        })
    })
}

fn multiply(left: Matrix, right: Matrix) -> Matrix {
    core::array::from_fn(|row| {
        core::array::from_fn(|col| {
            (0..3)
                .fold(Dd::ZERO, |sum, axis| {
                    sum + Dd::from_f64(left[row][axis]) * Dd::from_f64(right[axis][col])
                })
                .to_f64()
        })
    })
}

fn add_field(sum: &mut FieldSum, field: [[C64; 4]; 4]) -> Result<(), AeroacError> {
    for (row, values) in sum.iter_mut().zip(field) {
        for (entry, value) in row.iter_mut().zip(values) {
            entry.add(value)?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn fold_image(
    owner: &PreparedCompactSourceObserver,
    positions: &[[f64; 3]],
    planes: [IdealAcousticPlane; 2],
    transform: Matrix,
    sign: f64,
    sums: &mut [FieldSum],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AeroacError> {
    for &position in positions {
        plane::poll(cancelled)?;
        if planes[0].distance(position)?.to_f64() >= 0.0
            && planes[1].distance(position)?.to_f64() >= 0.0
        {
            return Err(invalid(
                "acoustic wedge image is not representably outside the fluid",
            ));
        }
    }
    let mut image = PreparedCompactSourceObserver::prepare_cancellable(
        owner.acoustic,
        owner.frequency_hz,
        owner.observer_m,
        positions,
        owner.source_count(),
        owner.wavenumber_per_m,
        owner.molecular.clone(),
        cancelled,
    )?;
    for (transfer, sum) in image.transfers.iter_mut().zip(sums) {
        plane::poll(cancelled)?;
        plane::transform_source_components(transfer, transform, sign)?;
        add_field(sum, transfer.field)?;
    }
    Ok(())
}
