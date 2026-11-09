//! Atomic source-temperature equilibrium on the original master coordinates.
use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::surface::ReferenceLoad3;
use crate::elastic3::thermal::equilibrium::{
    EquilibriumOwner3, admit_settings, resolve_candidate,
};
use crate::elastic3::thermal::{
    ThermalCoordinateMap3, ThermoelasticEquilibrium3, ThermoelasticEquilibriumSettings3,
};
use crate::elastic3::{ElasticMaterialField3, ElasticityError3};

impl AdaptiveElasticity3 {
    /// Solve the current native thermal state through the original hanging map.
    /// Both physical and independent DOF allowances precede native allocations.
    /// Load, physical admission, force gates, inertia and optional prestress use
    /// the same owners as the atomic refresh entry; no material token changes.
    ///
    /// This supports evaluation inside a material-decision callback. Acceptance
    /// remains with the enclosing transaction; rejected native preparations are
    /// stale against the restored operator. See
    /// [crate::elastic3::CutElasticity3::solve_thermoelastic_equilibrium] for scope and settings.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_thermoelastic_equilibrium(
        &self,
        load: ReferenceLoad3<'_>,
        nodal_force: Option<&[f64]>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        settings: ThermoelasticEquilibriumSettings3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticEquilibrium3, ElasticityError3> {
        admit_settings(self.n(), self.raw.n(), nodal_force, settings, &mut checkpoint)?;
        resolve_candidate(
            EquilibriumOwner3::Adaptive(self),
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            load, nodal_force, prescribed, settings, checkpoint,
        )
    }

    /// Atomically refresh native source materials and solve compatible thermal
    /// equilibrium through the SAME original hanging-node map and transpose.
    /// Independent force/displacement entries use original master coordinates;
    /// physical source and domain checks use original cell reconstruction.
    ///
    /// All physical scopes, result bindings, independent force/residual gates
    /// and rollback guarantees match
    /// [`crate::elastic3::CutElasticity3::refresh_thermoelastic_equilibrium`].
    /// The constraint map, reference geometry and original mass inventory are
    /// neither copied nor rebuilt; no alternative solver or projection is used.
    #[allow(clippy::too_many_arguments)]
    pub fn refresh_thermoelastic_equilibrium(
        &mut self,
        field: ElasticMaterialField3<'_>,
        load: ReferenceLoad3<'_>,
        nodal_force: Option<&[f64]>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        settings: ThermoelasticEquilibriumSettings3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticEquilibrium3, ElasticityError3> {
        admit_settings(self.n(), self.raw.n(), nodal_force, settings, &mut checkpoint)?;
        let completed = Cell::new(0);
        let checkpoint = RefCell::new(checkpoint);
        self.refresh_material_field_with(field, settings.refresh, |candidate| {
            resolve_candidate(EquilibriumOwner3::Adaptive(candidate),
                ThermalCoordinateMap3 {
                    rows: Some(&candidate.rows), fixed: &candidate.fixed,
                }, load, nodal_force, prescribed, settings, |iterations| {
                    completed.set(iterations);
                    (checkpoint.borrow_mut())(iterations)
                })
        }, || (checkpoint.borrow_mut())(completed.get()))
    }
}
