//! Outgoing Helmholtz point kernels, with time dependence exp(-i omega t).
//! These are off-diagonal kernels, not panel quadrature or finite-part self terms.

use fs_math::{c64::C64, det};

/// A complex-valued kernel smooth away from x=y. `Sync`, so FMM passes may
/// evaluate it on declared workers.
pub trait ComplexKernel: Sync {
    /// Evaluate K(target, source). The FMM excludes self interactions by index.
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64;

    /// Helmholtz wavenumber in inverse coordinate units, for regime admission.
    /// Generic non-Helmholtz kernels leave this as `None`.
    fn wavenumber(&self) -> Option<f64> {
        None
    }
    /// Whether K(x+a,y+a)=K(x,y), enabling shared translation couplings.
    /// Custom kernels must opt in; declaring a wavenumber alone is insufficient.
    fn translation_invariant(&self) -> bool {
        false
    }

    /// Exact analytic identity, when this kernel is a Cartesian component of
    /// the outgoing Helmholtz gradient. Directional interpolation can then
    /// retain its coordinate-difference factor exactly. Implementors returning
    /// `Some(g)` promise that `eval` and `wavenumber` agree with `g`.
    fn helmholtz_gradient(&self) -> Option<HelmholtzGradient> {
        None
    }
}

/// Outgoing single layer G=exp(ikr)/(4 pi r), zero on the diagonal.
#[derive(Debug, Clone, Copy)]
pub struct HelmholtzSingle {
    /// Nonnegative real wavenumber (inverse coordinate units).
    pub k: f64,
}

impl ComplexKernel for HelmholtzSingle {
    fn translation_invariant(&self) -> bool {
        true
    }
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64 {
        radial(self.k, x, y).map_or(C64::ZERO, |v| v.g)
    }

    fn wavenumber(&self) -> Option<f64> {
        Some(self.k)
    }
}

/// Source-normal derivative dG/dn_y; normals belong to the source cloud.
#[derive(Debug, Clone, Copy)]
pub struct HelmholtzDouble {
    /// Nonnegative real wavenumber.
    pub k: f64,
}

impl HelmholtzDouble {
    /// Evaluate dG/dn_y, excluding the diagonal. No jump or self term is added.
    #[must_use]
    pub fn eval_with_normal(self, x: [f64; 3], y: [f64; 3], ny: [f64; 3]) -> C64 {
        radial(self.k, x, y).map_or(C64::ZERO, |v| v.gp.scale(-dot(ny, v.unit)))
    }
}

/// Target-normal derivative dG/dn_x (adjoint double layer).
#[derive(Debug, Clone, Copy)]
pub struct HelmholtzAdjointDouble {
    /// Nonnegative real wavenumber.
    pub k: f64,
}

impl HelmholtzAdjointDouble {
    /// Evaluate dG/dn_x, excluding the diagonal.
    #[must_use]
    pub fn eval_with_normal(self, x: [f64; 3], y: [f64; 3], nx: [f64; 3]) -> C64 {
        radial(self.k, x, y).map_or(C64::ZERO, |v| v.gp.scale(dot(nx, v.unit)))
    }
}

/// Mixed normal derivative d²G/(dn_x dn_y), for the off-diagonal BM operator.
/// The caller supplies panel quadrature, finite-part self terms and jumps.
#[derive(Debug, Clone, Copy)]
pub struct HelmholtzHypersingular {
    /// Nonnegative real wavenumber.
    pub k: f64,
}

impl HelmholtzHypersingular {
    /// Evaluate the mixed derivative away from x=y.
    #[must_use]
    pub fn eval_with_normals(self, x: [f64; 3], y: [f64; 3], nx: [f64; 3], ny: [f64; 3]) -> C64 {
        radial(self.k, x, y).map_or(C64::ZERO, |v| {
            -(v.gpp - v.gp.scale(v.inv_r)).scale(dot(nx, v.unit) * dot(ny, v.unit))
                - v.gp.scale(v.inv_r * dot(nx, ny))
        })
    }
}

/// A Cartesian component of grad_x G. Three scalar applications give the
/// gradient used in regularized Burton–Miller formulations.
#[derive(Debug, Clone, Copy)]
pub struct HelmholtzGradient {
    pub(super) k: f64,
    pub(super) axis: usize,
}

impl HelmholtzGradient {
    /// Select a Cartesian component (0, 1, or 2).
    pub fn new(k: f64, axis: usize) -> Result<Self, crate::FmmError> {
        if !k.is_finite() || k < 0.0 || axis >= 3 {
            return Err(crate::FmmError::ComplexRequest {
                reason: "gradient requires finite k >= 0 and Cartesian axis < 3",
            });
        }
        Ok(Self { k, axis })
    }
}

impl ComplexKernel for HelmholtzGradient {
    fn helmholtz_gradient(&self) -> Option<HelmholtzGradient> {
        Some(*self)
    }
    fn translation_invariant(&self) -> bool {
        true
    }
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64 {
        radial(self.k, x, y).map_or(C64::ZERO, |v| v.gp.scale(v.unit[self.axis]))
    }

    fn wavenumber(&self) -> Option<f64> {
        Some(self.k)
    }
}

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) struct Radial {
    pub unit: [f64; 3],
    pub inv_r: f64,
    pub g: C64,
    pub gp: C64,
    pub gpp: C64,
}

/// G'=(ik-1/r)G; G''=(ik-1/r)G'+G/r². Scaled norm avoids
/// underflowing r² for otherwise representable, small separations.
pub(super) fn radial(k: f64, x: [f64; 3], y: [f64; 3]) -> Option<Radial> {
    let d = [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let scale = d[0].abs().max(d[1].abs()).max(d[2].abs());
    if scale == 0.0 {
        return None;
    }
    let scaled = d.map(|v| v / scale);
    let r = scale * det::sqrt(dot(scaled, scaled));
    let inv_r = 1.0 / r;
    let (s, c) = det::sin_cos(k * r);
    let g = C64::new(c, s).scale(inv_r / (4.0 * core::f64::consts::PI));
    let gp = g * C64::new(-inv_r, k);
    Some(Radial {
        unit: d.map(|v| v * inv_r),
        inv_r,
        g,
        gp,
        gpp: gp * C64::new(-inv_r, k) + g.scale(inv_r * inv_r),
    })
}
