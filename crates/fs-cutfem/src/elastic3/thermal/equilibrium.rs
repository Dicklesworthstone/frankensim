//! Atomic source-temperature changes through the original equilibrium owners.
//!
//! The material transaction owns publication. This module composes the native
//! complete RHS, controlled solve, physical strain admission and original
//! reference-density inertia; it introduces no heat equation or time scheme.

use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::ThermalCoordinateMap3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::eigenstrain::{PreparedEigenstrain3, physical_response};
use crate::elastic3::equilibrium::EquilibriumResponse3;
use crate::elastic3::mass::MaterialMass3;
use crate::elastic3::surface::ReferenceLoad3;
use crate::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticitySolution3,
    MaterialRefreshLimits3, PreparedPrestress3,
};

/// Independent bounds for coefficient refresh, equilibrium and retained results.
/// Nodal forces and residual allowances are in newtons; displacements are in
/// metres. No unbounded default is supplied. Limits apply before the
/// corresponding native allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticEquilibriumSettings3 {
    pub refresh: MaterialRefreshLimits3,
    /// Maximum BOTH original physical and independent displacement dimensions,
    /// including clamps. Adaptive load owners may use physical-coordinate scratch.
    pub maximum_native_dofs: usize,
    /// True relative Euclidean residual gate of the original controlled solve.
    pub relative_residual_tolerance: f64,
    pub maximum_iterations: usize,
    pub poll_iterations: usize,
    /// Independent absolute gate on the recovered free-coordinate force norm.
    /// Strong-support reactions are not counted as an equilibrium defect.
    pub maximum_free_residual_n: f64,
    /// Actual and B-bar constitutive TOTAL-strain Frobenius ceiling in (0,1).
    /// Original source component and ELASTIC-strain domains remain separate.
    pub maximum_kinematic_strain: f64,
    /// Bulk, selected embedded-support and ghost points in the original grid.
    /// No all-point response vector is retained by this operation.
    pub maximum_points: usize,
    /// Upper-triangular consistent-mass scalar slots, before xyz replication.
    pub maximum_mass_terms: usize,
    /// None skips geometric stiffness. Some bounds its retained scalar entries.
    /// Signed compression is retained; G and K+G need not be positive definite.
    pub maximum_prestress_entries: Option<usize>,
    /// Reference origin for original strong-support moment recovery.
    pub moment_origin_m: [f64; 3],
}

/// Compatible equilibrium and preparations on the queried native material owner.
/// A surrounding material transaction owns acceptance of proposed source states.
/// Moving this packet does not copy geometry or constraints.
///
/// The complete RHS includes source free strain, independent external forces
/// and optional embedded lifting. Strong reactions have sign K_full u-rhs_full;
/// embedded reactions remain with the separate native support observer.
/// The algebraic residual and retained-point physical admission do not establish
/// uniqueness, coercivity, continuum accuracy or nonlinear stability.
///
/// Native eigenstrain/mass/prestress bindings must still be validated after any
/// subsequent material or geometry change. G additionally freezes exact scales.
/// No caloric storage, boundary work, heat exchange or total energy balance is
/// inferred here. Source thermodynamic observations keep their original owner.
#[derive(Debug)]
pub struct ThermoelasticEquilibrium3 {
    pub solution: ElasticitySolution3,
    pub equilibrium: EquilibriumResponse3,
    pub eigenstrain: PreparedEigenstrain3,
    pub mass: MaterialMass3,
    pub prestress: Option<PreparedPrestress3>,
}

#[derive(Clone, Copy)]
pub(in crate::elastic3) enum EquilibriumOwner3<'a> {
    Cartesian(&'a CutElasticity3),
    Adaptive(&'a AdaptiveElasticity3),
}

macro_rules! native {
    ($owner:expr, $method:ident $(, $argument:expr)* $(,)?) => {
        match $owner {
            EquilibriumOwner3::Cartesian(op) => op.$method($($argument),*),
            EquilibriumOwner3::Adaptive(op) => op.$method($($argument),*),
        }
    };
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

fn poll(
    iterations: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    if checkpoint(iterations).is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

pub(in crate::elastic3) fn admit_settings(
    n: usize,
    physical_n: usize,
    nodal_force: Option<&[f64]>,
    settings: ThermoelasticEquilibriumSettings3,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    poll(0, checkpoint)?;
    if n == 0 || n > settings.maximum_native_dofs
        || physical_n > settings.maximum_native_dofs
        || !settings.relative_residual_tolerance.is_finite()
        || settings.relative_residual_tolerance <= 0.0
        || settings.relative_residual_tolerance >= 1.0
        || settings.poll_iterations == 0
        || !settings.maximum_free_residual_n.is_finite()
        || settings.maximum_free_residual_n < 0.0
        || !settings.maximum_kinematic_strain.is_finite()
        || settings.maximum_kinematic_strain <= 0.0
        || settings.maximum_kinematic_strain >= 1.0
        || settings.moment_origin_m.iter().any(|v| !v.is_finite())
        || nodal_force.is_some_and(|force| force.len() != n)
    {
        return Err(invalid("invalid atomic thermoelastic equilibrium controls or force shape"));
    }
    if let Some(force) = nodal_force {
        for chunk in force.chunks(192) {
            poll(0, checkpoint)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("nonfinite atomic thermoelastic external force"));
            }
        }
    }
    poll(0, checkpoint)
}

// Bulk admission comes from the native equilibrium response below. Selected
// support strains must also be admitted when G is not requested: a numerical
// reaction functional by itself does not perform physical strain admission.
pub(super) fn admit_selected_supports(
    raw: &CutElasticity3,
    map: ThermalCoordinateMap3<'_>,
    prepared: &PreparedEigenstrain3,
    displacement: &[f64],
    maximum_kinematic_strain: f64,
    completed: usize,
    checkpoint: &mut impl FnMut(usize) -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    let Some(support) = &raw.embedded else { return poll(completed, checkpoint); };
    for (id, (cell, free)) in raw.cells.iter().zip(&prepared.cells).enumerate() {
        poll(completed, checkpoint)?;
        let selected = support.selected_points(id);
        if selected.is_empty() { continue; }
        if selected.len() != free.surface_strains.len() {
            return Err(invalid("thermoelastic support/free-strain alignment mismatch"));
        }
        let surface = cell.rules.surface()
            .ok_or(invalid("thermoelastic support rule is absent"))?;
        let local = std::array::from_fn(|i|
            map.value(displacement, cell.nodes[i / 3], i % 3));
        for (&index, &free_strain) in selected.iter().zip(&free.surface_strains) {
            poll(completed, checkpoint)?;
            let point = surface.points().get(index)
                .ok_or(invalid("thermoelastic selected support point is absent"))?;
            physical_response(cell, &local, point.position, free_strain, free.mean_trace,
                raw.scales[id], maximum_kinematic_strain)?;
        }
    }
    poll(completed, checkpoint)
}

/// Shared equilibrium evaluation on either original displacement space.
/// A surrounding material transaction, when present, owns publication of this
/// result and its source after its final cancellation check.
#[allow(clippy::too_many_arguments)]
pub(in crate::elastic3) fn resolve_candidate(
    owner: EquilibriumOwner3<'_>,
    map: ThermalCoordinateMap3<'_>,
    load: ReferenceLoad3<'_>,
    nodal_force: Option<&[f64]>,
    prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
    settings: ThermoelasticEquilibriumSettings3,
    mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
) -> Result<ThermoelasticEquilibrium3, ElasticityError3> {
    poll(0, &mut checkpoint)?;
    let raw = match owner {
        EquilibriumOwner3::Cartesian(op) => op,
        EquilibriumOwner3::Adaptive(op) => &op.raw,
    };
    // The existing source/point preflight includes selected support and ghost
    // work and refuses even one nonthermal cell before external load callbacks.
    raw.admit_thermal_temperature_points(settings.maximum_points, || checkpoint(0))?;
    let eigenstrain = native!(owner, prepare_thermoelastic,
        settings.maximum_points, || checkpoint(0))?;
    let mut rhs = native!(owner, eigenstrain_equilibrium_load,
        &eigenstrain, load, prescribed, || checkpoint(0))?;
    if let Some(external) = nodal_force {
        for (i, (force, additional)) in rhs.iter_mut().zip(external).enumerate() {
            if i % 192 == 0 { poll(0, &mut checkpoint)?; }
            *force += additional;
            if !force.is_finite() {
                return Err(invalid("atomic thermoelastic complete force overflow"));
            }
        }
    }
    let solution = native!(owner, solve_controlled, &rhs,
        settings.relative_residual_tolerance, settings.maximum_iterations,
        settings.poll_iterations, &mut checkpoint)?;
    let completed = solution.iterations();
    let equilibrium = native!(owner, eigenstrain_equilibrium_response,
        &eigenstrain, solution.coefficients(), &rhs, settings.maximum_kinematic_strain,
        settings.maximum_points, settings.moment_origin_m, || checkpoint(completed))?;
    if equilibrium.free_residual_norm > settings.maximum_free_residual_n {
        return Err(invalid("thermoelastic equilibrium exceeds the absolute force allowance"));
    }
    // Preparation admits the SAME bulk and selected support state when G is
    // wanted. Otherwise retain surface admission without assembling/storing G.
    let prestress = match settings.maximum_prestress_entries {
        Some(entries) => Some(native!(owner, prepare_eigenstrain_prestress,
            &eigenstrain, solution.coefficients(), settings.maximum_kinematic_strain,
            settings.maximum_points, entries, || checkpoint(completed))?),
        None => {
            admit_selected_supports(raw, map, &eigenstrain, solution.coefficients(),
                settings.maximum_kinematic_strain, completed, &mut checkpoint)?;
            None
        }
    };
    let mass = native!(owner, prepare_material_mass,
        settings.maximum_mass_terms, || checkpoint(completed))?;
    poll(completed, &mut checkpoint)?;
    Ok(ThermoelasticEquilibrium3 { solution, equilibrium, eigenstrain, mass, prestress })
}

impl CutElasticity3 {
    /// Solve the current source-owned thermal equilibrium without changing it.
    ///
    /// Reuses the original complete load, controlled solve, absolute free-force
    /// gate, total/elastic strain admission and reference mass/prestress owners.
    /// Returned preparations bind to this current material state; no source
    /// query, material refresh, ownership rotation or geometry copy occurs.
    ///
    /// This entry can evaluate a candidate inside a material-decision callback.
    /// The enclosing transaction owns acceptance; a packet returned from a
    /// rejected candidate has stale native preparations against the restored
    /// operator. The settings' refresh limits are unused by this nonmutating
    /// solve; its independent solve, point, DOF and preparation bounds apply.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_thermoelastic_equilibrium(
        &self,
        load: ReferenceLoad3<'_>,
        nodal_force: Option<&[f64]>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        settings: ThermoelasticEquilibriumSettings3,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ThermoelasticEquilibrium3, ElasticityError3> {
        admit_settings(self.n(), self.n(), nodal_force, settings, &mut checkpoint)?;
        resolve_candidate(
            EquilibriumOwner3::Cartesian(self),
            ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            load, nodal_force, prescribed, settings, checkpoint,
        )
    }

    /// Change the complete source-temperature field and solve its compatible
    /// equilibrium as one atomic native material operation.
    ///
    /// Uses original retained geometry, scales, full reference dead loads,
    /// optional independent nodal forces (including supported entries), native
    /// free strain and optional embedded displacement lifting. Strong supports
    /// remain homogeneous. Prescribed motion is a pure held reference law;
    /// it does not supply a support-work or follower-load model.
    ///
    /// The original controlled solver must pass its true relative residual gate
    /// and recovered free forces must pass the independent absolute allowance.
    /// TOTAL and ELASTIC source domains are admitted at original bulk and
    /// selected support points. No bulk-response array is allocated.
    ///
    /// Success commits source materials and returns matching owned displacement,
    /// eigenstrain, original-density consistent mass, reactions and optional G.
    /// Any material/solve/domain/budget/force refusal or cancellation restores
    /// the prior operator and ALL prior preparations. Geometry and reference
    /// inventory are never copied, reintegrated or changed. This is a fixed
    /// reference quasistatic update, not an atomic heat/mechanical macrostep.
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
        admit_settings(self.n(), self.n(), nodal_force, settings, &mut checkpoint)?;
        let completed = Cell::new(0);
        let checkpoint = RefCell::new(checkpoint);
        self.refresh_material_field_with(field, settings.refresh, |candidate| {
            resolve_candidate(EquilibriumOwner3::Cartesian(candidate),
                ThermalCoordinateMap3 { rows: None, fixed: &candidate.fixed },
                load, nodal_force, prescribed, settings, |iterations| {
                    completed.set(iterations);
                    (checkpoint.borrow_mut())(iterations)
                })
        }, || (checkpoint.borrow_mut())(completed.get()))
    }
}
