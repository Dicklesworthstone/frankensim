//! Reference geometry for fixed-plane, energy-regularized brittle solids.
//!
//! The existing [`crate::plastic3::SmallStrainTetProblem`] consumes
//! [`fs_material::brittle::BrittleCrackBand::small_strain_law`] directly. This
//! module supplies the element length needed to turn fracture energy per area
//! into dissipated energy per reference volume; it adds no constitutive or
//! equilibrium solver.
//!
//! For a convex tetrahedron and an explicitly selected crack normal `n`, the
//! projected area is `A = sum_faces |n . area_vector| / 2`, and the band width
//! is `h = V / A`. Consequently complete pure-mode failure dissipates
//! `V * Gc / h = Gc * A`. This is the tetrahedron's projected footprint, not
//! the area of a particular plane cutting its interior. The fixed normal is
//! a modeling input. Directional length regularization does not establish
//! mesh-objective crack paths, automatic crack initiation, separated topology,
//! or contact between newly created surfaces.

/// Retained source-bound brittle solids on the native P1 equilibrium path.
pub mod solid;

use crate::loads3::TetLoadGeometry3;
use fs_material::brittle::BrittleCrackBandSettings;

/// A bounded refusal while deriving a crack band from admitted reference geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TetCrackBandGeometryError {
    /// The index does not identify an element in the supplied original mesh.
    ElementOutOfRange { element: usize },
    /// A crack plane needs a finite, nonzero normal.
    InvalidNormal,
    /// The projected area or band width cannot be represented as positive SI data.
    UnrepresentableGeometry,
}

impl core::fmt::Display for TetCrackBandGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ElementOutOfRange { element } => {
                write!(f, "brittle crack-band element {element} is outside the reference mesh")
            }
            Self::InvalidNormal => f.write_str("brittle crack-band normal must be finite and nonzero"),
            Self::UnrepresentableGeometry => {
                f.write_str("brittle crack-band projected area or width is unrepresentable")
            }
        }
    }
}

impl std::error::Error for TetCrackBandGeometryError {}

/// Original tetrahedral volume, explicit crack direction, and directional length.
///
/// Values are immutable after admission. Original connectivity is retained to
/// make per-element assignment inspectable. Keep this record with the matching
/// material and element history; moving it to another mesh is not an admitted
/// history-transfer operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetCrackBandGeometry {
    element: usize,
    nodes: [usize; 4],
    normal: [f64; 3],
    reference_volume_m3: f64,
    projected_area_m2: f64,
    band_width_m: f64,
}

impl TetCrackBandGeometry {
    /// Derive a band from one element of the original admitted reference mesh.
    ///
    /// [`TetLoadGeometry3`] already owns finite coordinates, positive volume,
    /// connectivity, quality, work limits and face-incidence admission. Both
    /// element windings give the same geometric length. Mechanical consumers
    /// still apply their own orientation and constitutive admission rules.
    /// This bounded four-face query allocates no memory and builds no matrices.
    ///
    /// # Errors
    /// An unknown element, a nonfinite/zero normal, or unrepresentable positive
    /// projected area or width. No input geometry or history is modified.
    pub fn from_element(
        geometry: &TetLoadGeometry3<'_>,
        element: usize,
        normal: [f64; 3],
    ) -> Result<Self, TetCrackBandGeometryError> {
        let nodes = *geometry
            .tetrahedra()
            .get(element)
            .ok_or(TetCrackBandGeometryError::ElementOutOfRange { element })?;
        if normal.iter().any(|v| !v.is_finite()) {
            return Err(TetCrackBandGeometryError::InvalidNormal);
        }
        // Scaling first admits finite normals near either floating-point limit
        // without overflowing a squared norm or underflowing it to zero.
        let scale = normal.iter().fold(0.0_f64, |largest, v| largest.max(v.abs()));
        if scale == 0.0 {
            return Err(TetCrackBandGeometryError::InvalidNormal);
        }
        let scaled_normal = normal.map(|v| v / scale);
        let length = fs_math::det::sqrt(dot(scaled_normal, scaled_normal));
        let normal = scaled_normal.map(|v| v / length);
        let points = nodes.map(|node| geometry.nodes_m()[node]);
        let mut edge_scale = 0.0_f64;
        for a in 0..4 {
            for b in a + 1..4 {
                let edge = sub(points[b], points[a]);
                edge_scale = edge_scale.max(edge[0].hypot(edge[1]).hypot(edge[2]));
            }
        }
        // Scale all cross products as in the original load-geometry owner.
        // Translating the reference mesh does not enter any area product.
        let mut projection = 0.0;
        for [a, b, c] in [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] {
            let ab = sub(points[b], points[a]).map(|v| v / edge_scale);
            let ac = sub(points[c], points[a]).map(|v| v / edge_scale);
            projection += dot(normal, cross(ab, ac)).abs();
        }
        // One half converts a cross product to triangle area, the other half
        // removes the equal front/back contribution of the convex boundary.
        let projected_area_m2 = ((0.25 * projection) * edge_scale) * edge_scale;
        let reference_volume_m3 = geometry.element_volumes_m3()[element];
        let band_width_m = reference_volume_m3 / projected_area_m2;
        if !(projected_area_m2.is_finite()
            && projected_area_m2 > 0.0
            && band_width_m.is_finite()
            && band_width_m > 0.0)
        {
            return Err(TetCrackBandGeometryError::UnrepresentableGeometry);
        }
        Ok(Self {
            element,
            nodes,
            normal,
            reference_volume_m3,
            projected_area_m2,
            band_width_m,
        })
    }

    /// Original element index, in the supplied mesh's connectivity order.
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }

    /// Original global node IDs, preserving their local connectivity order.
    #[must_use]
    pub const fn nodes(&self) -> [usize; 4] {
        self.nodes
    }

    /// Unit crack-plane normal in the original reference-coordinate frame.
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }

    /// Undeformed element volume [m³], retained from the original geometry owner.
    #[must_use]
    pub const fn reference_volume_m3(&self) -> f64 {
        self.reference_volume_m3
    }

    /// Area [m²] of the tetrahedron's orthogonal projection onto the crack plane.
    #[must_use]
    pub const fn projected_area_m2(&self) -> f64 {
        self.projected_area_m2
    }

    /// Directional crack-band width `V / A_projected` [m].
    #[must_use]
    pub const fn band_width_m(&self) -> f64 {
        self.band_width_m
    }

    /// Supply the original length/direction to the material constructor.
    ///
    /// The caller specifies the physical small-strain domain and work tolerance.
    /// The material constructor still validates both; no defaults or calibrated
    /// strength/fracture data are inferred from mesh geometry.
    #[must_use]
    pub const fn material_settings(
        &self,
        maximum_strain: f64,
        energy_tolerance_j_m3: f64,
    ) -> BrittleCrackBandSettings {
        BrittleCrackBandSettings {
            band_width_m: self.band_width_m,
            normal: self.normal,
            maximum_strain,
            energy_tolerance_j_m3,
        }
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
