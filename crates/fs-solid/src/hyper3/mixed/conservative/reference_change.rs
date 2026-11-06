//! Prescribed natural-strain and constitutive-reference events.
//!
//! Rebinding is performed at the accepted total deformation, velocity and time.
//! Its signed input is the change of physical storage, integrated pointwise in
//! the ORIGINAL reference. Retained finite-bulk volumes are carried unchanged.
//! Original geometry, consistent mass, density and accumulated energy baseline
//! are shared; neither motion nor a pressure-space volume is reconstructed.
//!
//! An exact constraint may change natural volume only as part of an atomic
//! physical interval. For P = P_remainder + p_alg Gd, constraint actuation supplies
//! -sum_i p_alg_i integral N_i Delta(det(G)) dV. This sign follows by subtracting
//! the old and new weak constraints from the actual discrete pressure work.
//! The target law/reference is held during that interval: this is a prescribed
//! endpoint split, not a continuous thermal trajectory or a heat closure.

use std::sync::Arc;

use fs_exec::Cx;

use super::{
    ConservativeFormulation, ConservativeReferenceField, ConservativeTetDefinition,
    ConservativeTetError, ConservativeTetSolid, ConservativeTetStep, ConstraintReferenceChange,
    MixedLawField, MixedUpSettings, checkpoint, gate, invalid,
};

/// Fixed-deformation physical storage and its signed source input [J].
/// Bulk reference work is a subaccount already INCLUDED in reference_work_j.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConservativeReferenceEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub bulk_stored_before_j: f64,
    pub bulk_stored_after_j: f64,
    pub reference_work_j: f64,
    pub bulk_reference_work_j: f64,
    /// Pointwise integrated reference work minus the accepted storage change.
    pub balance_error_j: f64,
    /// All prior input plus this event, relative to the original total energy.
    pub cumulative_balance_error_j: f64,
}

#[derive(Debug)]
pub struct ConservativeTetReferenceChange<'a> {
    pub solid: ConservativeTetSolid<'a>,
    pub energy: ConservativeReferenceEnergy,
}

/// One atomic reference event followed by one native mechanical interval.
/// `reference` reports fixed-F source input; `step.energy` reports the interval
/// AFTER that event, including any exact-volume constraint actuation. Its
/// cumulative balance includes both source inputs exactly once.
#[derive(Debug)]
pub struct ConservativeTetReferenceStep<'a> {
    pub step: ConservativeTetStep<'a>,
    pub reference: ConservativeReferenceEnergy,
}

/// Only the combined operation can hold an off-constraint intermediate; it
/// cannot be returned as an accepted ConservativeTetSolid.
struct PreparedReferenceChange<'a> {
    solid: ConservativeTetSolid<'a>,
    energy: ConservativeReferenceEnergy,
    natural_volume_change_m3: Option<Vec<f64>>,
}

fn add_compensated(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

impl<'a> ConservativeTetSolid<'a> {
    /// Change the prescribed law/natural strain at FIXED accepted F, v and time.
    ///
    /// FullMaterial and OriginalMaterial admit the new elastic state. The
    /// latter retains each original TOTAL bulk volume and evaluates the new
    /// J_ref U(Jbar_total/J_ref); no volume or endpoint pressure is projected.
    /// ExactIncompressible also requires the new weak natural-volume constraint
    /// to hold at the unchanged displacement. Use the atomic interval entry for
    /// a volume-changing actuation. Accepted centered velocities are retained
    /// without reapplying preparation's initial tangent-rate restriction.
    ///
    /// The predictor is in sorted pressure-node order. An empty slice retains
    /// the last interval multiplier; a nonempty slice must cover every pressure
    /// node. It is only a Newton predictor, never a physical pressure reset.
    /// A changed finite-bulk branch may require an explicit admitted predictor.
    ///
    /// All material, work, volume and cancellation refusals leave self intact.
    pub fn change_reference<'b>(
        &self,
        material: MixedLawField<'b>,
        references: ConservativeReferenceField<'b>,
        pressure_predictor_pa: &[f64],
        maximum_work_error_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetReferenceChange<'b>, ConservativeTetError>
    where
        'a: 'b,
    {
        let prepared = self.prepare_reference_change(
            material,
            references,
            pressure_predictor_pa,
            maximum_work_error_j,
            false,
            cx,
        )?;
        Ok(ConservativeTetReferenceChange {
            solid: prepared.solid,
            energy: prepared.energy,
        })
    }

    /// Prescribe the new law/natural strain and advance one positive interval
    /// as a single accepted result. The fixed-F reference work is separated from
    /// the native interval's load/support and exact-volume actuation work.
    ///
    /// For exact incompressibility, the target natural volume may differ from
    /// the old one. Only the final native endpoint is published; its pressure
    /// rows, physical motion, original work identities and cumulative balance
    /// must all pass. The off-constraint fixed-F intermediate is private.
    ///
    /// This endpoint split assigns the new law/reference before the mechanical
    /// interval. No second-order accuracy in a varying temperature/reference
    /// path is claimed; callers control and refine that prescribed path. Both
    /// the old F evaluated with the target law/reference and the solved endpoint
    /// must lie in the target constitutive domain; refinement can be necessary.
    #[allow(clippy::too_many_arguments)]
    pub fn change_reference_and_advance<'b>(
        &self,
        material: MixedLawField<'b>,
        references: ConservativeReferenceField<'b>,
        pressure_predictor_pa: &[f64],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: MixedUpSettings,
        maximum_work_error_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetReferenceStep<'b>, ConservativeTetError>
    where
        'a: 'b,
    {
        let prepared = self.prepare_reference_change(
            material,
            references,
            pressure_predictor_pa,
            maximum_work_error_j,
            true,
            cx,
        )?;
        let reference_change =
            prepared
                .natural_volume_change_m3
                .as_deref()
                .map(|change| ConstraintReferenceChange {
                    pressure_residual_before_m3: &self.pressure_residual_m3,
                    natural_volume_change_m3: change,
                });
        let step = prepared.solid.advance_with_reference_change(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            maximum_work_error_j,
            reference_change,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(ConservativeTetReferenceStep {
            step,
            reference: prepared.energy,
        })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepare_reference_change<'b>(
        &self,
        material: MixedLawField<'b>,
        references: ConservativeReferenceField<'b>,
        pressure_predictor_pa: &[f64],
        maximum_work_error_j: f64,
        allow_constraint_change: bool,
        cx: &Cx<'_>,
    ) -> Result<PreparedReferenceChange<'b>, ConservativeTetError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        if !maximum_work_error_j.is_finite() || maximum_work_error_j < 0.0 {
            return Err(invalid(
                "finite nonnegative reference-work allowance required",
            ));
        }
        let points = self.points.len();
        let count = match material {
            MixedLawField::Uniform(_) => None,
            MixedLawField::PerElement(laws) => Some((laws.len(), self.mesh.elements.len())),
            MixedLawField::PerQuadrature(laws) => Some((laws.len(), points)),
        };
        if count.is_some_and(|(got, expected)| got != expected) {
            return Err(invalid(
                "target elastic laws must cover the original quadrature",
            ));
        }
        references.validate(self.mesh.elements.len(), points)?;
        let np = self.algorithmic_pressure_pa.len();
        if (!pressure_predictor_pa.is_empty() && pressure_predictor_pa.len() != np)
            || pressure_predictor_pa.iter().any(|p| !p.is_finite())
        {
            return Err(invalid(
                "finite reference-event predictor must cover the pressure nodes",
            ));
        }
        let algorithmic_pressure_pa = if pressure_predictor_pa.is_empty() {
            self.algorithmic_pressure_pa.clone()
        } else {
            pressure_predictor_pa.to_vec()
        };
        let definition: ConservativeTetDefinition<'b> = ConservativeTetDefinition {
            material,
            ..self.definition
        };
        let frozen = definition.freeze_points(
            &self.mesh,
            &self.displacement_m,
            definition
                .formulation
                .bulk_settings()
                .map(|_| self.bulk_volume_ratios.as_slice()),
            &algorithmic_pressure_pa,
            references,
            cx,
        )?;
        let exact = matches!(
            definition.formulation,
            ConservativeFormulation::ExactIncompressible { .. }
        );
        if !exact || !allow_constraint_change {
            definition.admit_volume(&self.mesh, &frozen.pressure_residual_m3)?;
        }

        let mut material_work_j = 0.0;
        let mut material_correction_j = 0.0;
        let mut bulk_reference_work_j = 0.0;
        let mut bulk_correction_j = 0.0;
        let mut natural_volume_change_m3 =
            (exact && allow_constraint_change).then(|| vec![0.0; np]);
        for (element_index, element) in self.mesh.elements.iter().enumerate() {
            for (point_index, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element_index * self.mesh.points_per_element + point_index;
                add_compensated(
                    &mut material_work_j,
                    &mut material_correction_j,
                    point.weight_m3
                        * (frozen.points[index].energy_before_pa()
                            - self.points[index].energy_before_pa()),
                );
                if !self.bulk_points.is_empty() {
                    add_compensated(
                        &mut bulk_reference_work_j,
                        &mut bulk_correction_j,
                        point.weight_m3
                            * (frozen.bulk_points[index].energy_before_pa()
                                - self.bulk_points[index].energy_before_pa()),
                    );
                }
                if let Some(change) = &mut natural_volume_change_m3 {
                    let old_volume = self
                        .references
                        .at(element_index, point_index, self.mesh.points_per_element)
                        .map_or(1.0, |reference| reference.volume_ratio());
                    let new_volume = references
                        .at(element_index, point_index, self.mesh.points_per_element)
                        .map_or(1.0, |reference| reference.volume_ratio());
                    // Integrate natural-volume change directly. Subtracting two
                    // geometric residuals can erase a small prescribed change.
                    let delta = new_volume - old_volume;
                    for a in 0..4 {
                        change[element.pressure_indices[a]] +=
                            point.weight_m3 * point.pressure_shapes[a] * delta;
                    }
                }
            }
        }
        if natural_volume_change_m3
            .as_ref()
            .is_some_and(|change| change.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid("natural-volume actuation is unrepresentable"));
        }
        let reference_work_j = material_work_j + bulk_reference_work_j;
        let increment = reference_work_j - self.work_correction_j;
        let cumulative_work_j = self.cumulative_work_j + increment;
        let work_correction_j = (cumulative_work_j - self.cumulative_work_j) - increment;
        let energy = ConservativeReferenceEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: frozen.stored_energy_j,
            bulk_stored_before_j: self.bulk_stored_energy_j,
            bulk_stored_after_j: frozen.bulk_stored_energy_j,
            reference_work_j,
            bulk_reference_work_j,
            balance_error_j: reference_work_j - (frozen.stored_energy_j - self.stored_energy_j),
            cumulative_balance_error_j: cumulative_work_j
                - (frozen.stored_energy_j + self.kinetic_energy_j - self.initial_total_energy_j),
        };
        for (quantity, value) in [
            ("reference storage work", energy.balance_error_j),
            (
                "reference original bulk work",
                bulk_reference_work_j - (frozen.bulk_stored_energy_j - self.bulk_stored_energy_j),
            ),
            (
                "accumulated reference work",
                energy.cumulative_balance_error_j,
            ),
        ] {
            gate(quantity, value, maximum_work_error_j)?;
        }
        checkpoint(cx)?;
        Ok(PreparedReferenceChange {
            solid: ConservativeTetSolid {
                definition,
                references,
                mesh: Arc::clone(&self.mesh),
                mass: Arc::clone(&self.mass),
                corner_reference_density: self.corner_reference_density.clone(),
                points: frozen.points,
                bulk_points: frozen.bulk_points,
                bulk_volume_ratios: frozen.bulk_volume_ratios,
                displacement_m: self.displacement_m.clone(),
                velocity_mps: self.velocity_mps.clone(),
                algorithmic_pressure_pa,
                pressure_residual_m3: frozen.pressure_residual_m3,
                stored_energy_j: frozen.stored_energy_j,
                bulk_stored_energy_j: frozen.bulk_stored_energy_j,
                kinetic_energy_j: self.kinetic_energy_j,
                initial_total_energy_j: self.initial_total_energy_j,
                cumulative_work_j,
                work_correction_j,
                time_s: self.time_s,
            },
            energy,
            natural_volume_change_m3,
        })
    }
}
