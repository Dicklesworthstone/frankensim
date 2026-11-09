//! The original master-space constraints survive a native material refresh.
use super::*;
use crate::elastic3::refresh::{MaterialRefreshOwner3, refresh_material_owner};

impl MaterialRefreshOwner3 for AdaptiveElasticity3 {
    fn material_refresh_raw(&self) -> &CutElasticity3 {
        &self.raw
    }

    fn material_refresh_raw_mut(&mut self) -> &mut CutElasticity3 {
        &mut self.raw
    }
}

impl AdaptiveElasticity3 {
    /// Atomically replace material-dependent raw blocks and stabilization while
    /// retaining the SAME hanging reconstruction, masters, supports and geometry.
    /// Reference density and design scales remain fixed. All prepared native
    /// mass/body/eigenstrain/prestress maps become stale only after success.
    /// See [`CutElasticity3::refresh_material_field`] for work and physical scope.
    pub fn refresh_material_field(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.refresh_material_field_with(field, limits, |_| Ok(()), checkpoint)
    }

    /// Accept a native material field only after dependent work on this complete
    /// adaptive owner succeeds. The immutable candidate retains the original
    /// hanging reconstruction, masters and homogeneous supports, so existing
    /// reduced solves, recovery, mass and prestress preparation remain available.
    ///
    /// Error, final cancellation or an unwinding callback/checkpoint panic
    /// restores the previous raw caches and exact preparation binding. Geometry,
    /// master rows, scales and reference masses never change. Owned preparations
    /// returned by a successful attempt remain valid against the accepted owner.
    /// Staging uses the same buffers for rollback; geometry is never cloned.
    ///
    /// [CutElasticity3::refresh_material_field_with] defines the physical and
    /// transaction scope: callback work keeps its own solve/admission budgets;
    /// external mutations and independently nested transactions are not rolled
    /// back. Return proposed peer state or observations and publish after success.
    pub fn refresh_material_field_with<R>(
        &mut self,
        field: ElasticMaterialField3<'_>,
        limits: MaterialRefreshLimits3,
        attempt: impl FnOnce(&Self) -> Result<R, ElasticityError3>,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<R, ElasticityError3> {
        refresh_material_owner(self, field, limits, attempt, checkpoint)
    }
}
