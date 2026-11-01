//! Full affine load-set assessment on two native adaptive grids.
//!
//! Both grids remaximize their complete quadratic, including the common
//! material body force. Two frozen load directions use the existing DWR owner;
//! neither a smooth envelope derivative nor a continuum bound is assumed.
use super::grid::{combine_fields, evaluate_grid, finite, normalize, relative_difference};
use super::{AffineLoadSpace3, AffineLoadSpaceGrid3, RobustOptions3};
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3, Sdf3Elasticity};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3, prepare_material_goal3};
use crate::{SolveControl, SolveWork};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::{ReferenceLoad3, SurfaceForce3};
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::reaction::{
    AffineGoalForm3, MaterialBodyLoads3, estimate_affine_assembled_goal3,
};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
use fs_dwr::elasticity3::{
    GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3, dorfler3,
};
use std::collections::BTreeMap;
use std::ops::ControlFlow;

/// Transfer and numerical identity limits beyond the shared solve budget.
/// Each supplied operator owns its preparation policy, including configured
/// AdaptiveSolveSpace3 Jacobi/two-level/multilevel preparation and accounting.
#[derive(Debug, Clone, Copy)]
pub struct AffineLoadSpaceRefinementOptions3 {
    pub max_transfer_terms: usize,
    pub numerical: GoalOptions3,
}
impl Default for AffineLoadSpaceRefinementOptions3 {
    fn default() -> Self {
        Self {
            max_transfer_terms: 2_000_000,
            numerical: GoalOptions3::default(),
        }
    }
}

/// Dilated material inventory on the two actual native grids. Fine occupancy
/// is copied from source parents, independently of eroded stiffness/self-weight.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceInventory3 {
    pub objective: MinimumMaterialObjective3,
    pub volume_scale: f64,
    pub coarse_volume_fraction: f64,
    pub fine_volume_fraction: f64,
    /// Absolute change in numerical dilated volume fraction.
    pub volume_change: f64,
    pub coarse_mass_kg: Option<f64>,
    pub fine_mass_kg: Option<f64>,
    /// Inventory divided by volume_scale or the selected MassKg scale_kg.
    pub coarse_objective: f64,
    pub fine_objective: f64,
    /// Absolute physical inventory change divided by its objective scale.
    pub objective_change: f64,
    /// Absolute dilated inventory change in kg, present only for MassKg.
    pub mass_change_kg: Option<f64>,
}

/// Complete successful two-grid experiment at one raw source design.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceRefinement3 {
    pub policy: AffineLoadSpace3,
    pub robustness: RobustOptions3,
    pub rho: Vec<f64>,
    pub coarse_eroded_rho: Vec<f64>,
    /// Literal parent copies; no fine filtering or reconstruction from SIMP.
    pub fine_eroded_rho: Vec<f64>,
    pub coarse_eroded_scales: Vec<f64>,
    pub fine_eroded_scales: Vec<f64>,
    pub coarse: AffineLoadSpaceGrid3,
    pub fine: AffineLoadSpaceGrid3,
    /// Optional DILATED inventory, evaluated without another filter or solve.
    pub inventory: Option<AffineLoadSpaceInventory3>,
    /// ||G_f-G_c||_F + 2||h_f-h_c||_2 + |c_f-c_c|, in normalized units.
    /// In exact arithmetic this bounds the quadratic change on the unit ball.
    /// This computed numerical metric is not an interval or continuum bound.
    pub normalized_model_change: f64,
    /// Difference of the independently remaximized reduced quadratic values.
    pub normalized_maximum_change: f64,
    /// Difference of actual physical endpoint work, divided by the reference.
    pub normalized_witness_change: f64,
    /// Physical compliance DWR for the coarse maximizing coefficient vector.
    pub coarse_witness: GoalEstimate3,
    /// Physical compliance DWR for the fine maximizing coefficient vector.
    pub fine_witness: GoalEstimate3,
    /// [J_f(a_f)-J_f(a_c)]/reference. Signed; never clipped or localized.
    pub fine_switch_correction: f64,
    /// [J_c(a_c)-J_c(a_f)]/reference. Signed; never clipped or localized.
    pub coarse_switch_correction: f64,
    /// Scaled accounting defect, including comparison to the oracle maxima.
    /// Scaling includes endpoint values and cancelling DWR terms; this is not
    /// relative accuracy of a nearly zero grid difference.
    pub identity_relative_defect: f64,
    /// Mean of the two frozen witnesses' absolute native local masses, divided
    /// by the reference. Diagnostic marking ONLY, not whole-envelope error.
    pub marking_mass: BTreeMap<Octant3, f64>,
    pub work: SolveWork,
}

fn inventory_scale(
    objective: MinimumMaterialObjective3,
    volume_scale: f64,
) -> Result<f64, GoalRefinementError3> {
    let scale = match objective {
        MinimumMaterialObjective3::VolumeFraction => volume_scale,
        MinimumMaterialObjective3::MassKg { scale_kg } => scale_kg,
    };
    // Mass mode carries volume_scale as metadata only, matching design admission.
    if !scale.is_finite() || scale <= 0.0 || !(1.0 / scale).is_finite() || 1.0 / scale == 0.0 {
        return Err(GoalRefinementError3::Invalid(
            "affine inventory requires positive invertible scales",
        ));
    }
    Ok(scale)
}

#[allow(clippy::too_many_arguments)]
fn assess_inventory<O: Sdf3Elasticity, F: Sdf3Elasticity>(
    coarse: &CutDensityStudy3<O>,
    fine: &F,
    projected: &[f64],
    parents: &[usize],
    objective: MinimumMaterialObjective3,
    volume_scale: f64,
    control: &mut SolveControl<'_>,
) -> Result<AffineLoadSpaceInventory3, GoalRefinementError3> {
    let scale = inventory_scale(objective, volume_scale)?;
    if parents.len() != fine.cells() || projected.len() != coarse.cells() {
        return Err(GoalRefinementError3::Invalid(
            "affine inventory parent-map shape mismatch",
        ));
    }
    let mut fine_projected = Vec::with_capacity(parents.len());
    for (i, &parent) in parents.iter().enumerate() {
        if i % 64 == 0 {
            control.checkpoint("sdf3-affine-goal-dilated-occupancy")?;
        }
        fine_projected.push(*projected.get(parent).ok_or(GoalRefinementError3::Invalid(
            "invalid affine inventory parent",
        ))?);
    }
    let coarse_volume_fraction =
        crate::sdf3::inventory::weighted_occupancy(&coarse.mass, projected, control)?;
    let fine_volume_fraction =
        crate::sdf3::inventory::volume_fraction(fine, &fine_projected, control)?;
    let volume_change = (fine_volume_fraction - coarse_volume_fraction).abs();
    let (coarse_mass_kg, fine_mass_kg, coarse_value, fine_value, mass_change_kg) = match objective {
        MinimumMaterialObjective3::VolumeFraction => (
            None,
            None,
            coarse_volume_fraction,
            fine_volume_fraction,
            None,
        ),
        MinimumMaterialObjective3::MassKg { .. } => {
            let a = crate::sdf3::inventory::weighted_occupancy(
                crate::sdf3::inventory::mass_weights(&coarse.operator)?,
                projected,
                control,
            )?;
            let b = crate::sdf3::inventory::weighted_occupancy(
                crate::sdf3::inventory::mass_weights(fine)?,
                &fine_projected,
                control,
            )?;
            (Some(a), Some(b), a, b, Some((b - a).abs()))
        }
    };
    let coarse_objective = normalize(coarse_value, scale, "sdf3-affine-goal-inventory")?;
    let fine_objective = normalize(fine_value, scale, "sdf3-affine-goal-inventory")?;
    let objective_change = normalize(
        (fine_value - coarse_value).abs(),
        scale,
        "sdf3-affine-goal-inventory-change",
    )?;
    Ok(AffineLoadSpaceInventory3 {
        objective,
        volume_scale,
        coarse_volume_fraction,
        fine_volume_fraction,
        volume_change,
        coarse_mass_kg,
        fine_mass_kg,
        coarse_objective,
        fine_objective,
        objective_change,
        mass_change_kg,
    })
}
impl AffineLoadSpaceRefinement3 {
    /// Symmetric reconstruction of the normalized actual endpoint difference.
    /// The signed load-switch terms remain separate from the DWR cell masses.
    pub fn correction(&self) -> f64 {
        0.5 * (self.coarse_witness.correction() / self.policy.reference())
            + 0.5 * (self.fine_witness.correction() / self.policy.reference())
            + 0.5 * self.fine_switch_correction
            - 0.5 * self.coarse_switch_correction
    }

    /// Mark source leaves from the two frozen-direction diagnostics. Empty
    /// marks never establish that the complete uncertain load set is resolved;
    /// consult normalized_model_change as well as the full signed identities.
    pub fn mark(
        &self,
        theta: f64,
        max_marks: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, theta, max_marks, checkpoint)
    }

    /// Interpret the policy's reference as a PHYSICAL compliance cap, as in a
    /// minimum-material study. Both reduced maxima and actual witness work
    /// must satisfy that cap. The tolerance applies ONLY to full-model change;
    /// it never waives cap excess. A reference used solely for normalization
    /// need not be a design cap, so this interpretation is explicitly opt-in.
    /// This is numerical two-grid evidence, not continuum feasibility.
    pub fn stable_and_feasible(&self, change_tolerance: f64) -> bool {
        change_tolerance.is_finite()
            && change_tolerance >= 0.0
            && self.normalized_model_change.is_finite()
            && self.normalized_model_change >= 0.0
            && self.normalized_model_change <= change_tolerance
            && [&self.coarse, &self.fine].iter().all(|grid| {
                grid.maximum.is_finite()
                    && (0.0..=1.0).contains(&grid.maximum)
                    && grid.worst_compliance.is_finite()
                    && grid.worst_compliance >= 0.0
                    && grid.worst_compliance <= self.policy.reference()
            })
    }
}

fn norm_difference(a: &[f64], b: &[f64]) -> Result<f64, GoalRefinementError3> {
    let differences: Vec<_> = a.iter().zip(b).map(|(a, b)| b - a).collect();
    if differences.iter().any(|v| !v.is_finite()) {
        return Err(GoalRefinementError3::Invalid(
            "affine model difference overflow",
        ));
    }
    let scale = differences.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let norm = if scale == 0.0 {
        0.0
    } else {
        scale * fs_math::det::sqrt(differences.iter().map(|v| (v / scale) * (v / scale)).sum())
    };
    finite(norm, "sdf3-affine-goal-model-change").map_err(Into::into)
}

// Each reference law is integrated independently on each actual grid. Frozen
// DWR laws sum the SAME physical densities; native B*rho is supplied separately.
#[allow(clippy::too_many_arguments)]
fn frozen_goal(
    transfer: &AdaptiveTransfer3<'_>,
    material: Option<MaterialBodyLoads3<'_>>,
    dead: ReferenceLoad3<'_>,
    basis: &[ReferenceLoad3<'_>],
    coefficients: &[f64],
    forces: (&[f64], &[f64]),
    fields: GoalFields3<'_>,
    options: GoalOptions3,
    control: &mut SolveControl<'_>,
) -> Result<GoalEstimate3, GoalRefinementError3> {
    let body = |p| {
        let mut value = [0.0; 3];
        for (law, weight) in std::iter::once((dead, 1.0))
            .chain(basis.iter().copied().zip(coefficients.iter().copied()))
        {
            if weight == 0.0 {
                continue;
            }
            if let Some(force) = law.body {
                let f = force(p);
                for c in 0..3 {
                    value[c] += weight * f[c];
                }
            }
        }
        value
    };
    let surface = |p, n: [f64; 3]| {
        let mut value = [0.0; 3];
        for (law, weight) in std::iter::once((dead, 1.0))
            .chain(basis.iter().copied().zip(coefficients.iter().copied()))
        {
            if weight == 0.0 {
                continue;
            }
            let f = match law.surface {
                Some(SurfaceForce3::Traction(force)) => force(p, n),
                Some(SurfaceForce3::Pressure(pressure)) => {
                    let pressure = pressure(p);
                    n.map(|v| -pressure * v)
                }
                None => continue,
            };
            for c in 0..3 {
                value[c] += weight * f[c];
            }
        }
        value
    };
    let external = ReferenceLoad3 {
        body: if dead.body.is_some() || basis.iter().any(|law| law.body.is_some()) {
            Some(&body)
        } else {
            None
        },
        surface: if dead.surface.is_some() || basis.iter().any(|law| law.surface.is_some()) {
            Some(SurfaceForce3::Traction(&surface))
        } else {
            None
        },
    };
    // The quadratic's c is already contained in total force dot displacement;
    // it is NOT a constant offset of this frozen linear compliance goal.
    let coarse = AffineGoalForm3 {
        gradient: forces.0.to_vec(),
        offset: 0.0,
    };
    let fine = AffineGoalForm3 {
        gradient: forces.1.to_vec(),
        offset: 0.0,
    };
    controlled_goal3(control, "sdf3-affine-goal-witness", |poll| {
        estimate_affine_assembled_goal3(
            transfer,
            EquilibriumLoad3 {
                external,
                prescribed: None,
            },
            material,
            &coarse,
            &fine,
            fields,
            options,
            poll,
        )
        .map(|report| report.linearized)
    })
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Reintegrate and independently remaximize the entire affine load set on
    /// the source and enriched grid. Baseline and basis laws are fixed, pure
    /// reference loads, with amplitudes/correlations already encoded in them.
    /// The same implicit domain, supports and material assignment are required.
    /// Full native material identity is checked before load callbacks; source
    /// physical occupancy and stiffness are inherited by fine children exactly.
    /// A configured body load requires a retained material acceleration recipe.
    ///
    /// One source filter and one preparation plus dimension+1 direct solves on
    /// each grid suffice. Frozen-witness DWR adds no primal, adjoint, or filter
    /// solve. Each configured backend's preparation/work policy is retained.
    /// Incoming source scales are restored on every returned outcome. The fine
    /// operator is only a probe and is not installed as a refitted study.
    ///
    /// Hard/flat maxima are admitted without asserting a differentiable
    /// envelope. DWR is physical self-adjoint compliance under homogeneous
    /// supports. This does not implement imposed motion, adaptive optimization,
    /// shape derivatives, or an error certificate for the uncertain continuum.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_affine_load_space_enrichment<F: AdaptiveSdf3Elasticity>(
        &mut self,
        enriched: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        options: AffineLoadSpaceRefinementOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefinement3, GoalRefinementError3> {
        self.affine_load_space_enrichment(
            enriched, rho, dead, basis, robustness, policy, None, options, control,
        )
    }

    /// The same complete affine experiment with DILATED inventory assessment.
    /// Source eroded and dilated projections share one forward filter. Fine
    /// dilated occupancy is copied from the same checked source parents, then
    /// integrated against that grid's own volume and native mass weights.
    /// No fine filter, inventory adjoint or additional elasticity solve occurs.
    /// VolumeFraction does not require density; MassKg refuses either grid
    /// without complete native material density. The existing load-set/DWR
    /// identities and eroded fields retain their original meaning.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_affine_load_space_enrichment_with_inventory<F: AdaptiveSdf3Elasticity>(
        &mut self,
        enriched: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        objective: MinimumMaterialObjective3,
        volume_scale: f64,
        options: AffineLoadSpaceRefinementOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefinement3, GoalRefinementError3> {
        self.affine_load_space_enrichment(
            enriched,
            rho,
            dead,
            basis,
            robustness,
            policy,
            Some((objective, volume_scale)),
            options,
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn affine_load_space_enrichment<F: AdaptiveSdf3Elasticity>(
        &mut self,
        mut enriched: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        inventory_policy: Option<(MinimumMaterialObjective3, f64)>,
        options: AffineLoadSpaceRefinementOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefinement3, GoalRefinementError3> {
        control.checkpoint("sdf3-affine-goal-start")?;
        self.require_material_body_recipe("sdf3-affine-goal-material-recipe")?;
        if !robustness.valid(self.params)
            || policy.oracle_options().validate().is_err()
            || rho.len() != self.cells()
            || rho
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || basis.len() != policy.dimension()
            || options.max_transfer_terms == 0
            || ![
                options.numerical.residual_tolerance,
                options.numerical.identity_tolerance,
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0 && *v < 1.0)
        {
            return Err(GoalRefinementError3::Invalid(
                "invalid affine enrichment design/load family/options",
            ));
        }
        if let Some((objective, volume_scale)) = inventory_policy {
            inventory_scale(objective, volume_scale)?;
            if matches!(objective, MinimumMaterialObjective3::MassKg { .. }) {
                crate::sdf3::inventory::mass_weights(&self.operator)?;
                crate::sdf3::inventory::mass_weights(&enriched)?;
            }
        }
        let filtered = self.filter_apply(rho, control)?;
        let eroded = self.design_from_filtered(
            &filtered,
            self.params.eta + robustness.threshold_offset,
            control,
        )?;
        let dilated = inventory_policy
            .map(|_| {
                self.design_from_filtered(
                    &filtered,
                    self.params.eta - robustness.threshold_offset,
                    control,
                )
            })
            .transpose()?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&eroded.scales)?;
        let result = (|| {
            let initial = controlled_goal3(control, "sdf3-affine-goal-transfer", |poll| {
                AdaptiveTransfer3::new(
                    self.operator.adaptive(),
                    enriched.adaptive(),
                    options.max_transfer_terms,
                    poll,
                )
                .map_err(Into::into)
            })?;
            if !initial
                .fine()
                .leaves()
                .iter()
                .zip(initial.parents())
                .any(|(leaf, &parent)| leaf.level() > initial.coarse().leaves()[parent].level())
            {
                return Err(GoalRefinementError3::Invalid(
                    "affine goal requires active-space enrichment",
                ));
            }
            let fine_eroded_scales = initial.inherited_scales();
            let mut fine_eroded_rho = Vec::with_capacity(initial.parents().len());
            for (i, &parent) in initial.parents().iter().enumerate() {
                if i % 64 == 0 {
                    control.checkpoint("sdf3-affine-goal-occupancy")?;
                }
                fine_eroded_rho.push(eroded.projected[parent]);
            }
            let inventory = inventory_policy
                .map(|(objective, volume_scale)| {
                    assess_inventory(
                        self,
                        &enriched,
                        &dilated
                            .as_ref()
                            .expect("inventory projection was prepared")
                            .projected,
                        initial.parents(),
                        objective,
                        volume_scale,
                        control,
                    )
                })
                .transpose()?;
            drop(initial);
            enriched.set_scales(&fine_eroded_scales)?;
            let transfer = controlled_goal3(control, "sdf3-affine-goal-transfer", |poll| {
                AdaptiveTransfer3::new(
                    self.operator.adaptive(),
                    enriched.adaptive(),
                    options.max_transfer_terms,
                    poll,
                )
                .map_err(Into::into)
            })?;
            let material = prepare_material_goal3(
                self,
                &enriched,
                transfer.parents(),
                Some(&eroded.projected),
                control,
            )?;
            let mut coarse_offset = controlled_goal3(control, "sdf3-affine-goal-load", |poll| {
                transfer
                    .coarse()
                    .reference_load(dead, poll)
                    .map_err(Into::into)
            })?;
            let mut fine_offset = controlled_goal3(control, "sdf3-affine-goal-load", |poll| {
                transfer
                    .fine()
                    .reference_load(dead, poll)
                    .map_err(Into::into)
            })?;
            crate::sdf3::body::add_force(
                &mut coarse_offset,
                material.as_ref().map(|m| m.coarse_force.as_slice()),
                control,
            )?;
            crate::sdf3::body::add_force(
                &mut fine_offset,
                material.as_ref().map(|m| m.fine_force.as_slice()),
                control,
            )?;
            let mut coarse_basis = Vec::with_capacity(basis.len());
            let mut fine_basis = Vec::with_capacity(basis.len());
            for &law in basis {
                coarse_basis.push(controlled_goal3(
                    control,
                    "sdf3-affine-goal-load",
                    |poll| {
                        transfer
                            .coarse()
                            .reference_load(law, poll)
                            .map_err(Into::into)
                    },
                )?);
                fine_basis.push(controlled_goal3(
                    control,
                    "sdf3-affine-goal-load",
                    |poll| {
                        transfer
                            .fine()
                            .reference_load(law, poll)
                            .map_err(Into::into)
                    },
                )?);
            }
            let coarse_basis: Vec<_> = coarse_basis.iter().map(Vec::as_slice).collect();
            let fine_basis: Vec<_> = fine_basis.iter().map(Vec::as_slice).collect();
            let coarse = {
                let prepared = self.operator.prepare_elasticity(control)?;
                evaluate_grid(
                    &self.operator,
                    &prepared,
                    coarse_offset,
                    &coarse_basis,
                    policy,
                    control,
                )?
            };
            let fine = {
                let prepared = enriched.prepare_elasticity(control)?;
                evaluate_grid(
                    &enriched,
                    &prepared,
                    fine_offset,
                    &fine_basis,
                    policy,
                    control,
                )?
            };
            let normalized_model_change = finite(
                norm_difference(&coarse.gram, &fine.gram)?
                    + 2.0 * norm_difference(&coarse.linear, &fine.linear)?
                    + (fine.offset - coarse.offset).abs(),
                "sdf3-affine-goal-model-change",
            )?;
            let normalized_maximum_change = finite(
                fine.maximum - coarse.maximum,
                "sdf3-affine-goal-maximum-change",
            )?;
            let coarse_witness = {
                let (force, displacement) = combine_fields(
                    enriched.fixed(),
                    &fine.offset_force,
                    &fine.offset_displacement,
                    &fine_basis,
                    &fine.basis_displacements,
                    &coarse.worst_coefficients,
                    control,
                )?;
                frozen_goal(
                    &transfer,
                    material.as_ref().map(|m| m.loads()),
                    dead,
                    basis,
                    &coarse.worst_coefficients,
                    (&coarse.worst_force, &force),
                    GoalFields3::compliance(&coarse.worst_displacement, &displacement),
                    options.numerical,
                    control,
                )?
            };
            let fine_witness = {
                let (force, displacement) = combine_fields(
                    self.operator.fixed(),
                    &coarse.offset_force,
                    &coarse.offset_displacement,
                    &coarse_basis,
                    &coarse.basis_displacements,
                    &fine.worst_coefficients,
                    control,
                )?;
                frozen_goal(
                    &transfer,
                    material.as_ref().map(|m| m.loads()),
                    dead,
                    basis,
                    &fine.worst_coefficients,
                    (&force, &fine.worst_force),
                    GoalFields3::compliance(&displacement, &fine.worst_displacement),
                    options.numerical,
                    control,
                )?
            };
            let reference = policy.reference();
            let normalized = |value| normalize(value, reference, "sdf3-affine-goal-normalization");
            let [a, b, c, d] = [
                coarse_witness.coarse_value,
                coarse_witness.fine_value,
                fine_witness.coarse_value,
                fine_witness.fine_value,
            ]
            .map(normalized);
            let [a, b, c, d] = [a?, b?, c?, d?];
            let normalized_witness_change =
                normalized(fine_witness.fine_value - coarse_witness.coarse_value)?;
            let fine_switch_correction =
                normalized(fine_witness.fine_value - coarse_witness.fine_value)?;
            let coarse_switch_correction =
                normalized(coarse_witness.coarse_value - fine_witness.coarse_value)?;
            let dc = normalized(coarse_witness.correction())?;
            let df = normalized(fine_witness.correction())?;
            let reconstructed = finite(
                0.5 * dc + 0.5 * df + 0.5 * fine_switch_correction - 0.5 * coarse_switch_correction,
                "sdf3-affine-goal-identity",
            )?;
            let mut scale = [
                a,
                b,
                c,
                d,
                dc,
                df,
                reconstructed,
                normalized_witness_change,
                normalized_maximum_change,
                fine_switch_correction,
                coarse_switch_correction,
                coarse.maximum,
                fine.maximum,
            ]
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
            for goal in [&coarse_witness, &fine_witness] {
                for value in [
                    goal.dwr,
                    goal.coarse_space,
                    goal.goal_transfer,
                    goal.algebraic,
                ] {
                    scale = scale.max(normalized(value)?.abs());
                }
            }
            let mut identity_relative_defect = coarse_witness
                .identity_relative_defect
                .max(fine_witness.identity_relative_defect);
            if scale > 0.0 {
                identity_relative_defect = identity_relative_defect
                    .max(
                        ((d / scale - a / scale) - (dc / scale + fine_switch_correction / scale))
                            .abs(),
                    )
                    .max(
                        ((d / scale - a / scale) - (df / scale - coarse_switch_correction / scale))
                            .abs(),
                    )
                    .max(relative_difference(
                        normalized_witness_change,
                        reconstructed,
                        scale,
                    ))
                    .max(relative_difference(
                        normalized_maximum_change,
                        normalized_witness_change,
                        scale,
                    ))
                    .max(relative_difference(coarse.maximum, a, scale))
                    .max(relative_difference(fine.maximum, d, scale));
            }
            if !identity_relative_defect.is_finite()
                || identity_relative_defect > options.numerical.identity_tolerance
            {
                return Err(GoalError3::Identity {
                    relative_defect: identity_relative_defect,
                }
                .into());
            }
            let mut marking_mass = BTreeMap::new();
            for goal in [&coarse_witness, &fine_witness] {
                for (i, (&leaf, cell)) in goal.cells.iter().enumerate() {
                    if i % 64 == 0 {
                        control.checkpoint("sdf3-affine-goal-marking")?;
                    }
                    let value = normalized(cell.marking_mass)?;
                    let half = 0.5 * value;
                    if value < 0.0 || (value != 0.0 && half == 0.0) {
                        return Err(GoalRefinementError3::Invalid(
                            "affine marking mass is not representable",
                        ));
                    }
                    let entry = marking_mass.entry(leaf).or_insert(0.0);
                    *entry = finite(*entry + half, "sdf3-affine-goal-marking")?;
                }
            }
            drop(material);
            control.checkpoint("sdf3-affine-goal-publish")?;
            Ok(AffineLoadSpaceRefinement3 {
                policy: policy.clone(),
                robustness,
                rho: rho.to_vec(),
                coarse_eroded_rho: eroded.projected,
                fine_eroded_rho,
                coarse_eroded_scales: eroded.scales,
                fine_eroded_scales,
                coarse,
                fine,
                inventory,
                normalized_model_change,
                normalized_maximum_change,
                normalized_witness_change,
                coarse_witness,
                fine_witness,
                fine_switch_correction,
                coarse_switch_correction,
                identity_relative_defect,
                marking_mass,
                work: control.work(),
            })
        })();
        self.operator
            .set_scales(&previous)
            .expect("previous affine enrichment scales remain valid");
        result
    }
}
