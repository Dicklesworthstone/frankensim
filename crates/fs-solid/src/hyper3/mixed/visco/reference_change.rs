//! Held-configuration reference changes on the original mixed energy owner.
//!
//! Pressure remains an algebraic variable. The caller supplies its compatible
//! endpoint values; this transaction neither invents a pressure projection nor
//! advances displacement to repair an incompatible exact-volume constraint.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::MixedFiniteStrainPotential;

use super::super::pressure::cofactor;
use super::{
    MixedUpError, MixedViscoLawField, MixedViscoTetError, MixedViscoTetSolid, checkpoint, gate,
    invalid,
};

/// Independent weak-volume and pressure-work allowances for a reference change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedViscoReferenceChangeControls {
    /// Native normalized norm `||r_p / original_reference_volume||₂`.
    pub maximum_relative_volume_residual: f64,
    /// Maximum absolute `∫p*r_p` at either endpoint and of its change [J].
    /// This numerical constraint-work allowance is never viscous heat.
    pub maximum_pressure_constraint_work_j: f64,
}

/// Complete proposed reference change, including the mixed constraint defect.
#[derive(Debug)]
pub struct MixedViscoTetReferenceChange<'a> {
    pub solid: MixedViscoTetSolid<'a>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed native bulk-free exchange, with the same retained viscous metrics.
    pub remainder_reference_work_j: f64,
    /// Change of the original mixed bulk saddle potential at held total F.
    /// At exact weak pressure equilibrium this equals the bulk storage change.
    pub bulk_reference_work_j: f64,
    /// Signed external reference work: remainder plus bulk saddle change.
    pub reference_work_j: f64,
    /// `Δ∫p*r_p`, reported and gated separately from physical dissipation.
    /// `Δstored = reference_work - pressure_constraint_defect` to roundoff.
    pub pressure_constraint_defect_j: f64,
    /// Actual native weak-volume norm at the proposed endpoint.
    pub relative_volume_residual: f64,
}

impl<'a> MixedViscoTetSolid<'a> {
    /// Change stress-free references at fixed displacement, time and memory.
    ///
    /// Pressure is supplied in [`Self::pressure_nodes`] order. The original
    /// material owner permits only reference/relaxation-time changes and retains
    /// total F, natural-coordinate viscous metrics, mechanical work and loss.
    /// Pressure must remain exactly unchanged if no reference map changed.
    ///
    /// Every new weak volume equation must satisfy the declared allowance.
    /// Exact-incompressible volume changes that need mechanical motion refuse;
    /// the subsequent motion/equilibrium problem must solve that displacement.
    /// The original pressure dual supplies physical bulk storage. Its saddle
    /// value equals storage plus `p*r_p`, so the pressure constraint defect is
    /// separated from the signed reference work and independently gated at both
    /// endpoints. Neither account is mechanical loading work or viscous heat.
    ///
    /// Success shares original geometry/displacement and publishes a new
    /// accepted identity. Failure leaves the entire previous solid untouched.
    pub fn change_reference<'b>(
        &self,
        target_material: MixedViscoLawField<'b>,
        pressure_pa: &[f64],
        controls: MixedViscoReferenceChangeControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoTetReferenceChange<'b>, MixedViscoTetError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        if pressure_pa.len() != self.pressure_pa.len()
            || pressure_pa.iter().any(|p| !p.is_finite())
            || [
                controls.maximum_relative_volume_residual,
                controls.maximum_pressure_constraint_work_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "reference change needs finite pressure at every corner and positive finite volume/work allowances",
            ).into());
        }
        let expected = match target_material {
            MixedViscoLawField::Uniform(_) => None,
            MixedViscoLawField::PerElement(laws) => Some((laws.len(), self.mesh.elements.len())),
            MixedViscoLawField::PerQuadrature(laws) => Some((laws.len(), self.states.len())),
        };
        if expected.is_some_and(|(got, wanted)| got != wanted) {
            return Err(invalid("reference-change laws must cover every admitted point").into());
        }
        let points = self.mesh.points_per_element;
        let mut states = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        let mut remainder_after_j = 0.0;
        let mut bulk_after_j = 0.0;
        let mut remainder_reference_work_j = 0.0;
        let mut old_constraint_work_j = 0.0;
        let mut new_constraint_work_j = 0.0;
        let mut volume_residual = vec![0.0; self.mesh.pressure_nodes.len()];
        let mut reference_changed = false;
        for (e, element) in self.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = e * points + q;
                let previous = self.material.at(e, q, points);
                let target = target_material.at(e, q, points);
                reference_changed |=
                    previous.reference_deformation() != target.reference_deformation();
                let material_error = |error| MixedUpError::Material {
                    element: e,
                    point: q,
                    error,
                };
                let change = target
                    .change_remainder_reference_from(previous, &self.states[index])
                    .map_err(material_error)?;
                let old_pressure = (0..4)
                    .map(|a| {
                        self.pressure_pa[element.pressure_indices[a]] * point.pressure_shapes[a]
                    })
                    .sum();
                let new_pressure = (0..4)
                    .map(|a| pressure_pa[element.pressure_indices[a]] * point.pressure_shapes[a])
                    .sum();
                // The reference update already prepared the new remainder.
                // Reuse it; pressure reconstruction is only the scalar dual.
                let f = &self.states[index].deformation_gradient;
                let cof = cofactor(f);
                let volume_ratio = f[0] * cof[0] + f[1] * cof[1] + f[2] * cof[2];
                let old_bulk = previous
                    .freeze_remainder(&self.states[index], 0.0)
                    .map_err(material_error)?
                    .mixed_pressure_potential(f, old_pressure, self.volumetric)
                    .map_err(material_error)?;
                let new_bulk = target
                    .freeze_remainder(&change.state, 0.0)
                    .map_err(material_error)?
                    .mixed_pressure_potential(f, new_pressure, self.volumetric)
                    .map_err(material_error)?;
                // Preserve the native pressure dual's resolved near-reference
                // increment rather than subtracting two rounded volume ratios.
                let old_residual = (volume_ratio - 1.0) - old_bulk.volume_increment;
                let new_residual = (volume_ratio - 1.0) - new_bulk.volume_increment;
                let total_piola: [f64; 9] =
                    core::array::from_fn(|i| change.first_piola_pa[i] + new_pressure * cof[i]);
                if total_piola.iter().any(|value| !value.is_finite()) {
                    return Err(material_error(fs_material::MaterialError::State {
                        what: "mixed reference-change Piola stress is unrepresentable".into(),
                    })
                    .into());
                }
                let w = point.weight_m3;
                remainder_after_j += w * change.state.stored;
                bulk_after_j += w * new_bulk.bulk_stored_pa;
                remainder_reference_work_j += w * change.reference_work_pa;
                old_constraint_work_j += w * old_pressure * old_residual;
                new_constraint_work_j += w * new_pressure * new_residual;
                for a in 0..4 {
                    volume_residual[element.pressure_indices[a]] +=
                        w * point.pressure_shapes[a] * new_residual;
                }
                first_piola_pa.push(total_piola);
                states.push(change.state);
            }
        }
        if !reference_changed && pressure_pa != self.pressure_pa.as_ref() {
            return Err(invalid(
                "unchanged material references require unchanged physical pressure",
            )
            .into());
        }
        for residual in &mut volume_residual {
            *residual /= self.mesh.total_volume_m3;
        }
        let relative_volume_residual = fs_solver::norm2(&volume_residual);
        if !relative_volume_residual.is_finite()
            || relative_volume_residual > controls.maximum_relative_volume_residual
        {
            return Err(MixedViscoTetError::InitialVolumeGate {
                relative_residual: relative_volume_residual,
                maximum_relative_residual: controls.maximum_relative_volume_residual,
            });
        }
        let pressure_constraint_defect_j = new_constraint_work_j - old_constraint_work_j;
        for (quantity, value) in [
            ("previous pressure constraint work", old_constraint_work_j),
            ("target pressure constraint work", new_constraint_work_j),
            (
                "reference-change pressure constraint defect",
                pressure_constraint_defect_j,
            ),
        ] {
            gate(quantity, value, controls.maximum_pressure_constraint_work_j)?;
        }
        // A kinetic-only change is an exact history/storage no-op. New identity
        // still prevents using a split receipt prepared with different kinetics.
        if !reference_changed {
            remainder_after_j = self.remainder_energy_j;
            bulk_after_j = self.bulk_energy_j;
        }
        let bulk_reference_work_j =
            (bulk_after_j - self.bulk_energy_j) + pressure_constraint_defect_j;
        let reference_work_j = remainder_reference_work_j + bulk_reference_work_j;
        let stored_after_j = remainder_after_j + bulk_after_j;
        if !stored_after_j.is_finite()
            || !remainder_reference_work_j.is_finite()
            || !bulk_reference_work_j.is_finite()
            || !reference_work_j.is_finite()
        {
            return Err(invalid("mixed reference-change energy is unrepresentable").into());
        }
        checkpoint(cx)?;
        Ok(MixedViscoTetReferenceChange {
            solid: MixedViscoTetSolid {
                nodes_m: self.nodes_m,
                elements: self.elements,
                material: target_material,
                volumetric: self.volumetric,
                mesh: Arc::clone(&self.mesh),
                budget: self.budget,
                quadrature_points_per_dir: self.quadrature_points_per_dir,
                maximum_quadrature_points: self.maximum_quadrature_points,
                displacement_m: Arc::clone(&self.displacement_m),
                pressure_pa: if reference_changed {
                    Arc::from(pressure_pa)
                } else {
                    Arc::clone(&self.pressure_pa)
                },
                states: if reference_changed {
                    states.into()
                } else {
                    Arc::clone(&self.states)
                },
                first_piola_pa: if reference_changed {
                    first_piola_pa.into()
                } else {
                    Arc::clone(&self.first_piola_pa)
                },
                remainder_energy_j: remainder_after_j,
                bulk_energy_j: bulk_after_j,
                time_s: self.time_s,
                accepted_identity: Arc::new(()),
            },
            stored_before_j: self.stored_energy_j(),
            stored_after_j,
            remainder_reference_work_j,
            bulk_reference_work_j,
            reference_work_j,
            pressure_constraint_defect_j,
            relative_volume_residual,
        })
    }
}
