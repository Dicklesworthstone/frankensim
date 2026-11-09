//! Existing physical source projections on the exact lossy fluid/wall owner.

use super::{
    BoundaryVelocitySource, HeterogeneousError, PreparedEquivalentFluid, PressureSourceOwner,
    VolumeVelocitySource,
};

impl<'f, 'c> PreparedEquivalentFluid<'f, 'c> {
    /// Prepare an original distributed volume-injection shape [m^-3]. A complex
    /// peak volume-flow amplitude [m^3/s] supplies the continuity source [1/s].
    ///
    /// The original tetrahedron index and physical position reach the callback;
    /// its signed shape is integrated without normalization. The SAME volume
    /// projection as HeterogeneousHelmholtzFe retains every pressure-basis test
    /// and the original affine/quadratic geometry. Its separate quadrature order
    /// must lie between pressure order + 2 and 16; source and mesh convergence
    /// remain caller obligations. No effective density, porosity, extra volume
    /// factor, heat-release model or mass-injection model is inferred.
    ///
    /// The source borrows this exact prepared frequency/material/wall owner.
    /// Subsequent drives reuse its load and the original LU without invoking
    /// this callback again. Original source errors propagate unchanged; failure
    /// or cancellation returns no partially projected source.
    pub fn volume_velocity_source<'p, E: From<HeterogeneousError>>(
        &'p self,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VolumeVelocitySource<'p, 'f, 'c>, E> {
        VolumeVelocitySource::prepare(
            PressureSourceOwner::Equivalent(self),
            quadrature_order,
            density_per_m3,
            cancelled,
        )
    }

    /// Prepare original quadratic inward boundary motion without discarding
    /// midside values. Each triangle supplies dimensionless values at
    /// `corner0, corner1, corner2, edge01, edge12, edge20`, in its given vertex
    /// order. A complex peak inward velocity [m/s] multiplies the cached load.
    ///
    /// Uses the SAME boundary source owner and original physical face rule as
    /// the real heterogeneous solver. The default rule resolves at least
    /// pressure order + 2 and the geometry's declared face quadrature order.
    /// Signed shapes may have zero area integral while still exciting the
    /// pressure field. No area normalization or nodal-to-corner reduction occurs.
    /// Sources retain this exact preparation; solve-time pressure-release and
    /// source-count admission applies even at zero amplitude.
    pub fn quadratic_velocity_source<'p>(
        &'p self,
        traces: &[([u32; 3], [f64; 6])],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'p, 'f, 'c>, HeterogeneousError> {
        self.quadratic_velocity_source_with_quadrature(
            traces,
            self.geometry.face_quadrature_order(self.geometry.order + 2),
            cancelled,
        )
    }

    /// The same physical source with explicit source/curved-metric quadrature,
    /// from pressure order + 2 through 16 Gauss points per face Duffy axis.
    /// This changes integration only; basis, geometry and material stay original.
    pub fn quadratic_velocity_source_with_quadrature<'p>(
        &'p self,
        traces: &[([u32; 3], [f64; 6])],
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'p, 'f, 'c>, HeterogeneousError> {
        BoundaryVelocitySource::prepare(
            PressureSourceOwner::Equivalent(self),
            traces,
            quadrature_order,
            cancelled,
        )
    }
}
