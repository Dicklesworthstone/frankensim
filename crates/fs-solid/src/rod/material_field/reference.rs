//! Intrinsic geometry for the existing source-bound Cosserat equilibrium.
//!
//! Stress resultants depend on Gamma-Gamma0 and kappa-kappa0, using the SAME
//! midpoint directors and SO(3) logarithm on reference and current geometry.
//! This is the intrinsic-strain constitutive relation, e.g. Gazzola et al.,
//! R. Soc. Open Sci. 5 (2018) 171628, doi:10.1098/rsos.171628; see also
//! https://www.cosseratrods.org/cosserat_rods/theory/ (Linear Elasticity).
//! Curvature of an unloaded part is not elastic pre-stress. A loaded scan must
//! NOT be passed as its stress-free shape; this API cannot infer residual stress.
//!
//! Nodes discretize explicit MATERIAL-coordinate intervals, not inferred spatial
//! chords. The original constructor retains uniform intervals; the nonuniform
//! entry retains every supplied length. Refine centerline, frames and sections
//! together. This is the existing slender-rod model, not a 3D solid, a scan
//! reconstruction, contact solve, or dynamics. No material names or target
//! resonances select its geometry.

use fs_ga::So3;
use fs_material::section::ConcentricShearSection;
use std::sync::Arc;

use super::{MAX_MATERIAL_ROD_SEGMENTS, Rod, RodSection, SpatialMaterialRod, SpatialRodError};

/// The same material-coordinate intervals drive strain, storage and inventory.
/// Uniform inputs retain the original node-multiplication arithmetic.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct ReferencePartition {
    pub(super) length_m: f64,
    pub(super) lengths_m: Vec<f64>,
    pub(super) coordinates_m: Vec<f64>,
}

impl ReferencePartition {
    pub(super) fn uniform(length_m: f64, count: usize) -> Result<Self, SpatialRodError> {
        if !length_m.is_finite()
            || length_m <= 0.0
            || !(1..=MAX_MATERIAL_ROD_SEGMENTS).contains(&count)
        {
            return Err(SpatialRodError::Invalid(
                "finite positive length and 1..=64 material segments are required",
            ));
        }
        let length = length_m / count as f64;
        if !length.is_finite() || length <= 0.0 {
            return Err(SpatialRodError::Invalid(
                "reference segment length is unrepresentable",
            ));
        }
        let coordinates_m = (0..=count)
            .map(|node| node as f64 * length)
            .collect::<Vec<_>>();
        if coordinates_m.iter().any(|x| !x.is_finite())
            || coordinates_m.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(SpatialRodError::Invalid(
                "reference material coordinates are unrepresentable",
            ));
        }
        Ok(Self {
            length_m,
            lengths_m: vec![length; count],
            coordinates_m,
        })
    }

    pub(super) fn nonuniform(lengths_m: Vec<f64>) -> Result<Self, SpatialRodError> {
        if !(1..=MAX_MATERIAL_ROD_SEGMENTS).contains(&lengths_m.len())
            || lengths_m
                .iter()
                .any(|length| !length.is_finite() || *length <= 0.0)
        {
            return Err(SpatialRodError::Invalid(
                "1..=64 finite positive reference segment lengths are required",
            ));
        }
        let mut coordinates_m = Vec::with_capacity(lengths_m.len() + 1);
        coordinates_m.push(0.0);
        let mut length_m = 0.0;
        for &length in &lengths_m {
            let next = length_m + length;
            if !next.is_finite() || next <= length_m {
                return Err(SpatialRodError::Invalid(
                    "reference segment boundaries must be finite and distinctly representable",
                ));
            }
            // Distinct stations alone can still shorten/stretch a tiny interval
            // by O(1) through cancellation in the actual nodal chord. Bound the
            // resulting dimensionless axial strain at the same roundoff scale
            // as reference-chord admission below, without changing authored L.
            if (((next - length_m) / length) - 1.0).abs() > 64.0 * f64::EPSILON {
                return Err(SpatialRodError::Invalid(
                    "reference material stations cannot represent the supplied segment length accurately",
                ));
            }
            coordinates_m.push(next);
            length_m = next;
        }
        Ok(Self {
            length_m,
            lengths_m,
            coordinates_m,
        })
    }
}

/// An explicit, immutable declaration of the unloaded reference configuration.
#[derive(Clone, Debug)]
pub struct RodReferenceGeometry {
    length_m: f64,
    pub(super) partition: ReferencePartition,
    positions_m: Vec<[f64; 3]>,
    frames: Vec<So3>,
    pub(super) strains: Vec<([f64; 3], [f64; 3])>,
}

impl RodReferenceGeometry {
    /// Admit a centerline and body-to-world material frames in SI units.
    ///
    /// `length_m` is the reference material length; `n` segments each occupy
    /// `length_m/n`. Supplied frames fix the material directions and intrinsic
    /// twist; they are not guessed from an exterior outline. The caller must
    /// resolve neighboring rotations below pi and supply noncoincident nodes.
    /// Every chord must fit within its material interval; an inconsistent
    /// length would otherwise change reference mass and compliance silently.
    /// No self-intersection or slenderness certificate is implied.
    pub fn try_new(
        length_m: f64,
        positions_m: Vec<[f64; 3]>,
        frames: Vec<So3>,
    ) -> Result<Self, SpatialRodError> {
        let partition = ReferencePartition::uniform(length_m, positions_m.len().saturating_sub(1))?;
        Self::from_partition(partition, positions_m, frames)
    }

    /// Admit explicit nonuniform material intervals and their unloaded geometry.
    /// Lengths are physical reference arc intervals, not inferred chord lengths.
    /// Their order is the node/section order and no interval is averaged or rescaled.
    /// Cumulative material stations must preserve each length to the same
    /// roundoff precision as the reference-chord admission.
    pub fn try_nonuniform(
        segment_lengths_m: Vec<f64>,
        positions_m: Vec<[f64; 3]>,
        frames: Vec<So3>,
    ) -> Result<Self, SpatialRodError> {
        let partition = ReferencePartition::nonuniform(segment_lengths_m)?;
        Self::from_partition(partition, positions_m, frames)
    }

    fn from_partition(
        partition: ReferencePartition,
        positions_m: Vec<[f64; 3]>,
        frames: Vec<So3>,
    ) -> Result<Self, SpatialRodError> {
        let invalid = SpatialRodError::Invalid;
        if positions_m.len() != partition.lengths_m.len() + 1
            || frames.len() != positions_m.len()
            || positions_m.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "reference rod requires finite positive length and matching 2..=65 nodes/frames",
            ));
        }
        let count = positions_m.len() - 1;
        // The kinematic owner does not consult stiffness. Actual coefficients
        // enter only when a material field is bound below.
        let rod = Rod {
            positions: positions_m,
            frames,
            l0: partition.lengths_m[0],
            section: RodSection {
                ea: 0.0,
                ga: 0.0,
                gj: 0.0,
                ei: 0.0,
            },
        };
        let mut strains = Vec::with_capacity(count);
        for segment in 0..count {
            let l0 = partition.lengths_m[segment];
            let chord =
                [0, 1, 2].map(|c| rod.positions[segment + 1][c] - rod.positions[segment][c]);
            let chord_length =
                fs_math::det::hypot(fs_math::det::hypot(chord[0], chord[1]), chord[2]);
            if !chord_length.is_finite() || chord_length <= 0.0 {
                return Err(invalid(
                    "reference centerline has a coincident or unrepresentable segment",
                ));
            }
            if chord_length - l0 > 64.0 * f64::EPSILON * l0 {
                return Err(invalid(
                    "reference chord exceeds its material-coordinate interval; supply the actual unloaded material length",
                ));
            }
            let (g, k) = crate::rod::derivatives::analytic::Kinematics::from_nodes(
                &rod.positions,
                &rod.frames,
                l0,
                segment,
            )?
            .strains();
            let angle = fs_math::det::hypot(fs_math::det::hypot(k[0] * l0, k[1] * l0), k[2] * l0);
            if !g.iter().chain(&k).all(|v| v.is_finite())
                || !angle.is_finite()
                || angle >= core::f64::consts::PI - 1e-8
            {
                return Err(invalid(
                    "reference directors need finite strains and neighboring rotation strictly below pi; refine the geometry",
                ));
            }
            strains.push((g, k));
        }
        Ok(Self {
            length_m: partition.length_m,
            partition,
            positions_m: rod.positions,
            frames: rod.frames,
            strains,
        })
    }

    /// Reference material length [m], independently of chord lengths.
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Original ordered material interval lengths [m].
    pub fn segment_lengths_m(&self) -> &[f64] {
        &self.partition.lengths_m
    }
    /// Original material-coordinate node stations [m], starting at zero.
    pub fn material_coordinates_m(&self) -> &[f64] {
        &self.partition.coordinates_m
    }
    /// Unloaded centerline nodes [m].
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }
    /// Unloaded body-to-world material frames.
    pub fn frames(&self) -> &[So3] {
        &self.frames
    }
    /// Reference strain measures from the same discretization as the solve.
    pub fn intrinsic_strains(&self) -> &[([f64; 3], [f64; 3])] {
        &self.strains
    }
}

impl SpatialMaterialRod {
    /// Bind actual material sections to an intrinsically curved/twisted part.
    /// All segments share the original equilibrium, thermal, gravity, support
    /// and reaction owners. The initial geometry is exactly stress free.
    pub fn try_from_reference(
        reference: RodReferenceGeometry,
        sections: Vec<ConcentricShearSection>,
    ) -> Result<Self, SpatialRodError> {
        if sections.len() != reference.strains.len() {
            return Err(SpatialRodError::Invalid(
                "one material section is required per reference segment",
            ));
        }
        let mut rod =
            Self::try_straight_with_partition(reference.partition.clone(), sections, None)?;
        rod.rod.positions.clone_from(&reference.positions_m);
        rod.rod.frames.clone_from(&reference.frames);
        rod.reference_geometry = Some(Arc::new(reference));
        rod.energy_j()?;
        Ok(rod)
    }

    /// Explicit intrinsic geometry; None retains the original straight datum.
    pub fn reference_geometry(&self) -> Option<&RodReferenceGeometry> {
        self.reference_geometry.as_deref()
    }
}
