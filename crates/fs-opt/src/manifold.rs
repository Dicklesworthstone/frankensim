//! Manifold metadata for variables: the gradient stack consumes
//! retraction + tangent projection so "optimize an orientation" never
//! degenerates into "optimize 9 numbers and renormalize when it
//! explodes". Euclidean, sphere, and SO(3) (unit quaternion) are
//! implemented; Stiefel and fixed-volume level sets are representable
//! metadata whose retractions land with their consumer beads.

/// Where a variable lives.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Manifold {
    /// Flat ℝⁿ.
    Euclidean(usize),
    /// Unit sphere S^{n−1} embedded in ℝⁿ.
    Sphere(usize),
    /// Rotations as unit quaternions (S³ with antipodal identification —
    /// the sign ambiguity is the caller's to canonicalize).
    So3,
    /// Stiefel(n, p): orthonormal n×p frames. METADATA-ONLY in v1.
    Stiefel(usize, usize),
    /// Fixed-volume level-set fields (topology optimization).
    /// METADATA-ONLY in v1.
    FixedVolumeLevelSet(usize),
}

impl Manifold {
    /// Ambient (storage) dimension.
    #[must_use]
    pub fn ambient_dim(&self) -> usize {
        match *self {
            Manifold::Euclidean(n) | Manifold::Sphere(n) | Manifold::FixedVolumeLevelSet(n) => n,
            Manifold::So3 => 4,
            Manifold::Stiefel(n, p) => n * p,
        }
    }

    /// Project an ambient vector onto the tangent space at `x`.
    /// Panics (structured) for metadata-only manifolds — using them in a
    /// descent before their consumer bead lands is a modeling error.
    #[must_use]
    pub fn tangent_project(&self, x: &[f64], v: &[f64]) -> Vec<f64> {
        match self {
            Manifold::Euclidean(_) => v.to_vec(),
            Manifold::Sphere(_) | Manifold::So3 => {
                // v − (v·x)·x on the unit sphere (x assumed unit).
                let dot: f64 = v.iter().zip(x).map(|(a, b)| a * b).sum();
                v.iter().zip(x).map(|(vi, xi)| (-dot).mul_add(*xi, *vi)).collect()
            }
            Manifold::Stiefel(..) | Manifold::FixedVolumeLevelSet(_) => {
                panic!("{self:?} is metadata-only in v1 (its consumer bead supplies the retraction)")
            }
        }
    }

    /// Retract from `x` along `step` (metric projection retractions).
    /// Same metadata-only panic policy as [`Manifold::tangent_project`].
    #[must_use]
    pub fn retract(&self, x: &[f64], step: &[f64]) -> Vec<f64> {
        match self {
            Manifold::Euclidean(_) => x.iter().zip(step).map(|(a, b)| a + b).collect(),
            Manifold::Sphere(_) | Manifold::So3 => {
                let moved: Vec<f64> = x.iter().zip(step).map(|(a, b)| a + b).collect();
                let nrm = fs_math::det::sqrt(moved.iter().map(|t| t * t).sum());
                moved.iter().map(|t| t / nrm).collect()
            }
            Manifold::Stiefel(..) | Manifold::FixedVolumeLevelSet(_) => {
                panic!("{self:?} is metadata-only in v1 (its consumer bead supplies the retraction)")
            }
        }
    }
}
