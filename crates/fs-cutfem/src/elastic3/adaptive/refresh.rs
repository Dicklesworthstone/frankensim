//! The original master-space constraints survive a native material refresh.
use super::*;

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
        self.raw.refresh_material_field(field, limits, checkpoint)
    }
}
