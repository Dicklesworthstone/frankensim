//! Prescribed free strain in the existing cut-cell elastic weak form.
//!
//! At fixed temperature/composition, `sigma = C:(epsilon(u)-epsilon_free)`.
//! The independent load is `integral B^T C epsilon_free dV`. Symmetric Nitsche
//! supports additionally require `-integral N^T (C epsilon_free)n dA`; omitting
//! this term creates artificial restraint even in a stress-free expansion.
//! See Bleyer, *Linear thermoelasticity (weak coupling)*, Numerical Tours of
//! Computational Mechanics with FEniCSx. No heat equation or expansion law is
//! invented here: callers supply already integrated, global-frame free strain.
//!
//! B-bar projects BOTH total and free dilatation with the retained physical
//! volume rule. The deviator remains local. Thus local pressure is constant
//! within a cell, as in the same condensed constant-pressure formulation.
//! Preparation retains sampled strains and reference load columns. Changing
//! stiffness scales reuses these columns without re-evaluating the field or
//! rebuilding geometry. Material identity, density and reference mass stay
//! with the original operator. These are quadrature-level numerical responses,
//! not certified field extrema, finite deformation or thermal back-reaction.

use std::ops::ControlFlow;
use std::sync::Arc;

use fs_material::tensor::{contract, from_mandel, to_mandel};
use fs_solver::op::LinearOp;

use super::material::{ReferenceMaterial3, engineering_strain};
use super::strain_domain::FrobeniusStrainBall;
use super::{Cell3, CutElasticity3, ElasticityError3, q1};

mod sensitivity;

pub(super) struct CellEigenstrain3 {
    /// Unscaled bulk and Nitsche load, including strongly constrained entries.
    pub(super) force: [f64; 24],
    /// Original tensor-shear samples, in exactly the retained bulk-rule order.
    pub(super) bulk_strains: Vec<[f64; 6]>,
    /// Physical-volume average of the original free dilatation.
    pub(super) mean_trace: f64,
}

/// Immutable free strain and reference loads bound to one original operator.
///
/// A rebuilt/refined operator requires its own preparation even when dimensions
/// match. Current stiffness scales are deliberately excluded from this binding.
/// Callback samples are numerical input, not source or interpolation authority.
pub struct PreparedEigenstrain3 {
    pub(super) binding: Arc<()>,
    pub(super) cells: Vec<CellEigenstrain3>,
    points: usize,
}

impl core::fmt::Debug for PreparedEigenstrain3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedEigenstrain3")
            .field("cells", &self.cells.len())
            .field("bulk_points", &self.points)
            .finish_non_exhaustive()
    }
}

impl PreparedEigenstrain3 {
    /// Number of retained physical bulk samples, excluding boundary samples.
    #[must_use]
    pub const fn bulk_points(&self) -> usize {
        self.points
    }
}

/// Physical elastic response at one original bulk quadrature point.
/// All strains and stresses use tensor shear `[xx, yy, zz, xy, yz, zx]`.
#[derive(Debug, Clone, PartialEq)]
pub struct EigenstrainResponse3 {
    /// Original active-cell index, also indexing material scales.
    pub cell: usize,
    /// Reference Cartesian position in metres.
    pub position: [f64; 3],
    /// Original positive reference-volume quadrature weight in m³.
    pub weight: f64,
    /// Actual total displacement strain, before any B-bar projection.
    pub strain: [f64; 6],
    /// Prescribed stress-free strain, before any B-bar projection.
    pub free_strain: [f64; 6],
    /// Actual total minus actual free strain.
    pub elastic_strain: [f64; 6],
    /// Strain used by the material; B-bar projects the elastic dilatation.
    pub constitutive_elastic_strain: [f64; 6],
    /// Original unscaled material stress in Pa.
    pub reference_stress: [f64; 6],
    /// Current material stress, `scale[cell] * reference_stress`, in Pa.
    pub stress: [f64; 6],
    /// Recoverable elastic energy per reference volume in J/m³.
    /// Excludes ghost/Nitsche terms and thermal or transformation energy.
    pub energy_density: f64,
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn finite_strain(strain: [f64; 6]) -> Result<[f64; 6], ElasticityError3> {
    if strain.iter().any(|v| !v.is_finite()) {
        return Err(ElasticityError3::Invalid(
            "nonfinite prescribed eigenstrain",
        ));
    }
    // Mandel conversion can itself overflow for a finite tensor shear input.
    if to_mandel(&strain).iter().any(|v| !v.is_finite()) {
        return Err(ElasticityError3::Invalid(
            "unrepresentable prescribed eigenstrain",
        ));
    }
    Ok(strain)
}

fn project_free(cell: &Cell3, mut strain: [f64; 6], mean_trace: f64) -> [f64; 6] {
    if cell.kinematics.is_projected() {
        let correction = (mean_trace - strain[..3].iter().sum::<f64>()) / 3.0;
        for value in &mut strain[..3] {
            *value += correction;
        }
    }
    strain
}

// Apply the SAME reference tensor as the bulk owner. Separating spherical and
// deviatoric modes avoids subtracting large hydrostatic stresses after thermal
// expansion. No native elastic-domain check is applied to a free-strain load;
// that domain belongs to the resulting elastic strain, not expansion alone.
fn material_stress(cell: &Cell3, strain: [f64; 6]) -> [f64; 6] {
    let isotropic = cell
        .kinematics
        .projected_moduli()
        .map(|(bulk, shear, _)| (bulk, shear))
        .or_else(|| match &cell.material {
            ReferenceMaterial3::Legacy {
                bulk,
                lame: [_, shear],
                ..
            } => Some((*bulk, *shear)),
            ReferenceMaterial3::Source(_) => None,
        });
    if let Some((bulk, shear)) = isotropic {
        let trace = strain[..3].iter().sum::<f64>();
        let mean = trace / 3.0;
        let pressure = bulk * trace;
        return std::array::from_fn(|i| {
            if i < 3 {
                (2.0 * shear).mul_add(strain[i] - mean, pressure)
            } else {
                core::f64::consts::SQRT_2 * shear * strain[i]
            }
        });
    }
    let c = cell
        .material
        .source_state()
        .expect("native tensor")
        .stiffness_mandel_pa();
    let mandel = std::array::from_fn(|i| c[i].iter().zip(strain).map(|(c, e)| c * e).sum());
    from_mandel(&mandel)
}

fn material_domain(cell: &Cell3, strain: [f64; 6]) -> Result<(), ElasticityError3> {
    let admitted = if let Some(limit) = cell.material.frobenius_strain_limit() {
        FrobeniusStrainBall::new(strain, limit).is_some()
    } else {
        let limit = cell
            .material
            .source_state()
            .expect("native tensor")
            .maximum_absolute_mandel_strain();
        strain.iter().all(|v| v.is_finite() && v.abs() <= limit)
    };
    if admitted {
        Ok(())
    } else {
        Err(ElasticityError3::Invalid(
            "eigenstrain response exceeds retained elastic material domain",
        ))
    }
}

impl CutElasticity3 {
    /// Integrate a prescribed thermal/transformation strain on the original
    /// geometry, material field, B-bar projection and embedded support patch.
    ///
    /// The pure callback receives active-cell index and reference position and
    /// returns tensor-shear strain in the global frame. For thermal expansion,
    /// integrate the original expansion law over temperature first; this API
    /// never assumes `alpha * delta_temperature` for a variable coefficient.
    /// `max_points` bounds all bulk plus selected boundary evaluations before
    /// any callback or sample allocation. Cancellation publishes no candidate.
    pub fn prepare_eigenstrain(
        &self,
        field: &dyn Fn(usize, [f64; 3]) -> Result<[f64; 6], ElasticityError3>,
        max_points: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedEigenstrain3, ElasticityError3> {
        poll(&mut checkpoint)?;
        let mut evaluations = 0usize;
        let mut bulk_points = 0usize;
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let surface = self
                .embedded
                .as_ref()
                .map_or(0, |data| data.selected_points(id).len());
            bulk_points = bulk_points.checked_add(cell.rules.bulk().len()).ok_or(
                ElasticityError3::Invalid("eigenstrain point count overflow"),
            )?;
            evaluations = evaluations
                .checked_add(cell.rules.bulk().len())
                .and_then(|n| n.checked_add(surface))
                .ok_or(ElasticityError3::Invalid(
                    "eigenstrain point count overflow",
                ))?;
            if evaluations > max_points {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain point allowance exhausted",
                ));
            }
        }
        let mut cells = Vec::with_capacity(self.cells.len());
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let mut bulk_strains = Vec::with_capacity(cell.rules.bulk().len());
            let mut mean_trace = 0.0;
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let strain = finite_strain(field(id, position)?)?;
                poll(&mut checkpoint)?;
                mean_trace += (weight / cell.rules.volume()) * strain[..3].iter().sum::<f64>();
                bulk_strains.push(strain);
            }
            if !mean_trace.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain mean dilatation overflow",
                ));
            }
            let mut force = [0.0; 24];
            for (&(position, weight), strain) in cell.rules.bulk().iter().zip(&bulk_strains) {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let free = project_free(cell, to_mandel(strain), mean_trace);
                let stress = material_stress(cell, free);
                for i in 3..24 {
                    let (a, component) = (i / 3, i % 3);
                    let mut basis = engineering_strain(gradients[a], component);
                    if let Some((_, _, mean_gradients)) = cell.kinematics.projected_moduli() {
                        let correction =
                            (mean_gradients[a][component] - gradients[a][component]) / 3.0;
                        for value in &mut basis[..3] {
                            *value += correction;
                        }
                    }
                    force[i] += weight * basis.iter().zip(stress).map(|(e, s)| e * s).sum::<f64>();
                }
            }
            // Bulk virtual strains annihilate translations exactly, using the
            // same u_a-u_0 convention as native stress observations.
            for component in 0..3 {
                force[component] = -(1..8).map(|a| force[3 * a + component]).sum::<f64>();
            }
            if let Some(data) = &self.embedded {
                let surface = cell.rules.surface().expect("retained Nitsche surface");
                for &index in data.selected_points(id) {
                    poll(&mut checkpoint)?;
                    let point = &surface.points()[index];
                    let strain = finite_strain(field(id, point.position)?)?;
                    poll(&mut checkpoint)?;
                    let stress =
                        material_stress(cell, project_free(cell, to_mandel(&strain), mean_trace));
                    let n = point.normal;
                    let traction = [
                        stress[0] * n[0] + stress[3] * n[1] + stress[5] * n[2],
                        stress[3] * n[0] + stress[1] * n[1] + stress[4] * n[2],
                        stress[5] * n[0] + stress[4] * n[1] + stress[2] * n[2],
                    ];
                    let (basis, _) = q1(cell.bounds, point.position);
                    for i in 0..24 {
                        force[i] -= point.weight * basis[i / 3] * traction[i % 3];
                    }
                }
            }
            if force.iter().any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "eigenstrain reference load overflow",
                ));
            }
            cells.push(CellEigenstrain3 {
                force,
                bulk_strains,
                mean_trace,
            });
        }
        poll(&mut checkpoint)?;
        Ok(PreparedEigenstrain3 {
            binding: Arc::clone(&self.material_body_binding),
            cells,
            points: bulk_points,
        })
    }

    /// Require the exact geometry/material/support owner used in preparation.
    pub fn validate_eigenstrain(
        &self,
        prepared: &PreparedEigenstrain3,
    ) -> Result<(), ElasticityError3> {
        if !Arc::ptr_eq(&self.material_body_binding, &prepared.binding) {
            return Err(ElasticityError3::Invalid(
                "eigenstrain belongs to another operator",
            ));
        }
        Ok(())
    }

    /// Current `sum scale[cell] * b_free[cell]`, in independent displacement
    /// coordinates. Add external and prescribed-displacement loads separately.
    /// A stiffness-design change therefore changes this load as well as K.
    pub fn eigenstrain_load(
        &self,
        prepared: &PreparedEigenstrain3,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        let mut force = vec![0.0; self.n()];
        for (id, (cell, free)) in self.cells.iter().zip(&prepared.cells).enumerate() {
            poll(&mut checkpoint)?;
            for i in 0..24 {
                let node = cell.nodes[i / 3];
                if !self.fixed[node] {
                    force[3 * node + i % 3] += self.scales[id] * free.force[i];
                }
            }
        }
        for chunk in force.chunks(256) {
            poll(&mut checkpoint)?;
            if chunk.iter().any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "scaled eigenstrain load overflow",
                ));
            }
        }
        poll(&mut checkpoint)?;
        Ok(force)
    }

    /// Observe total, free and elastic strain, stress and recoverable storage.
    ///
    /// `maximum_kinematic_strain` must lie in `(0,1)` and bounds the Frobenius
    /// norm of total actual and B-bar strain, separately from the original
    /// material's ELASTIC domain. A tighter linearization budget remains the
    /// caller's modeling choice; this ceiling is not a finite-strain error bound.
    /// Thus compatible free expansion may exceed a narrow elastic calibration
    /// while remaining within the user's small-deformation budget. Both actual
    /// and projected elastic strains must pass the original material envelope.
    /// These checks apply at retained bulk points; no between-point bound is
    /// asserted for an arbitrary sampled eigenstrain field.
    ///
    /// Use `solve_controlled` for equilibrium, then this method for physical
    /// recovery. `solve_physical_controlled` and `bulk_stress` describe zero
    /// eigenstrain and must not be used to interpret a thermally loaded state.
    pub fn eigenstrain_response(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<EigenstrainResponse3>, ElasticityError3> {
        self.eigenstrain_response_mapped(
            prepared,
            u,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn eigenstrain_response_mapped(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        n: usize,
        maximum_kinematic_strain: f64,
        max_points: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<EigenstrainResponse3>, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        if u.len() != n
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(ElasticityError3::Invalid(
                "invalid eigenstrain response displacement or kinematic bound",
            ));
        }
        if prepared.points > max_points {
            return Err(ElasticityError3::Invalid(
                "eigenstrain response point allowance exhausted",
            ));
        }
        for chunk in u.chunks(256) {
            poll(&mut checkpoint)?;
            if chunk.iter().any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "nonfinite eigenstrain response displacement",
                ));
            }
        }
        let mut response = Vec::with_capacity(prepared.points);
        for (id, (cell, free)) in self.cells.iter().zip(&prepared.cells).enumerate() {
            poll(&mut checkpoint)?;
            let local = local(cell, u);
            for (&(position, weight), &free_strain) in
                cell.rules.bulk().iter().zip(&free.bulk_strains)
            {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let strain = super::stress::mandel_strain(&gradients, &local);
                let assumed = cell.project_mandel_strain(strain, &local);
                for total in [strain, assumed] {
                    if FrobeniusStrainBall::new(total, maximum_kinematic_strain).is_none() {
                        return Err(ElasticityError3::Invalid(
                            "total strain exceeds small-deformation budget; finite-deformation model required",
                        ));
                    }
                }
                let eigen = to_mandel(&free_strain);
                let projected_free = project_free(cell, eigen, free.mean_trace);
                let elastic = std::array::from_fn(|i| strain[i] - eigen[i]);
                let constitutive = std::array::from_fn(|i| assumed[i] - projected_free[i]);
                material_domain(cell, elastic)?;
                material_domain(cell, constitutive)?;
                let reference_stress = material_stress(cell, constitutive);
                let stress = reference_stress.map(|s| self.scales[id] * s);
                let constitutive_elastic_strain = from_mandel(&constitutive);
                let energy_density = 0.5 * contract(&constitutive_elastic_strain, &stress);
                if reference_stress
                    .iter()
                    .chain(&stress)
                    .any(|v| !v.is_finite())
                    || !energy_density.is_finite()
                    || energy_density < 0.0
                {
                    return Err(ElasticityError3::Invalid(
                        "unrepresentable eigenstrain response",
                    ));
                }
                response.push(EigenstrainResponse3 {
                    cell: id,
                    position,
                    weight,
                    strain: from_mandel(&strain),
                    free_strain,
                    elastic_strain: from_mandel(&elastic),
                    constitutive_elastic_strain,
                    reference_stress,
                    stress,
                    energy_density,
                });
            }
        }
        poll(&mut checkpoint)?;
        Ok(response)
    }
}
