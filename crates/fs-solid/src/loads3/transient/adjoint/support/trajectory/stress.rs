//! Stress objectives on the original accepted full-displacement trajectory.
use super::{
    MaterialScaleGradient3, NodalMotionCotangent3, NodalMotionObservation3,
    PrescribedDesignGradient3, RecordedPrescribedTransient3, ReverseVisit3, TransientAdjointError3,
    add_material, check, filled, number, reserved,
};
use fs_exec::Cx;

/// A scalar objective's derivative with respect to one element's stress.
/// State zero denotes initialization; later states are accepted endpoints.
/// Values are objective/Pa in the mesh-frame Mandel order [xx,yy,zz,xy,yz,zx],
/// with sqrt(2)-scaled shear. They are stress cotangents, not physical stresses.
/// Include integration, quadrature or least-squares weights explicitly.
/// Duplicate element/state rows add in caller order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementStressObservation3 {
    pub state: usize,
    pub element: usize,
    pub values: [f64; 6],
}

impl RecordedPrescribedTransient3<'_, '_, '_, '_> {
    /// Differentiate mixed element-stress and nodal-motion objectives through
    /// every original accepted endpoint, including the initial acceleration
    /// solve. Enable material recording before the first finite step.
    ///
    /// Stress is C(Bq - free_thermal_strain). Both the full displacement seed
    /// B' C' w and the direct stiffness-scale term w' stress are retained.
    /// Their implicit and direct contributions are accumulated in one reverse
    /// sweep; density and damping enter through the original trajectory.
    /// Support displacements and applied loads remain independent controls.
    ///
    /// maximum_strain is the same finite positive total/elastic Mandel-strain
    /// ceiling as stress recovery, checked at states with stress observations.
    /// The source mesh frame, geometry, material model, free thermal strain,
    /// restraint axes and accepted grid remain fixed. Frame transformations
    /// and state-dependent material/thermal laws require their own chain rules.
    /// Refusal and cancellation publish no partial result and leave the
    /// recorded trajectory reusable.
    pub fn reverse_design_with_stress(
        &self,
        terminal: &NodalMotionCotangent3,
        motion_observations: &[NodalMotionObservation3],
        stress_observations: &[ElementStressObservation3],
        maximum_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<PrescribedDesignGradient3, TransientAdjointError3> {
        check(cx)?;
        if !maximum_strain.is_finite() || maximum_strain <= 0.0 {
            return Err(TransientAdjointError3::Invalid(
                "stress design requires a finite positive strain ceiling",
            ));
        }
        let tape = self
            .material
            .as_ref()
            .ok_or(TransientAdjointError3::Invalid(
                "stress design needs motion recording enabled before advancement",
            ))?;
        let count = motion_observations
            .len()
            .checked_add(stress_observations.len())
            .ok_or(TransientAdjointError3::Budget)?;
        if count > self.budget.maximum_observations {
            return Err(TransientAdjointError3::Budget);
        }
        let elements = self.run.geometry.tets.len();
        for row in stress_observations {
            check(cx)?;
            if row.state >= self.sample_count() || row.element >= elements {
                return Err(TransientAdjointError3::Invalid(
                    "stress observation must name an accepted state and existing element",
                ));
            }
            for value in row.values {
                number(value)?;
            }
        }
        let mut heads = filled(self.sample_count(), usize::MAX)?;
        let mut next = filled(stress_observations.len(), usize::MAX)?;
        for (index, row) in stress_observations.iter().enumerate().rev() {
            check(cx)?;
            next[index] = heads[row.state];
            heads[row.state] = index;
        }
        let mut weights = filled(elements, [0.0; 6])?;
        let mut rows = reserved(elements)?;
        let mut material = MaterialScaleGradient3 {
            element_scales: filled(elements, [0.0; 2])?,
            rayleigh: [0.0; 2],
        };
        let controls = self.reverse_accumulating(
            terminal,
            motion_observations,
            |visit| {
                match visit {
                    ReverseVisit3::Observation { state, seed } => {
                        if heads[state] == usize::MAX {
                            return Ok(());
                        }
                        weights.fill([0.0; 6]);
                        rows.clear();
                        let mut index = heads[state];
                        while index != usize::MAX {
                            check(cx)?;
                            let row = stress_observations[index];
                            for c in 0..6 {
                                weights[row.element][c] =
                                    number(weights[row.element][c] + row.values[c])?;
                            }
                            index = next[index];
                        }
                        for (element, &value) in weights.iter().enumerate() {
                            check(cx)?;
                            if value.iter().any(|v| *v != 0.0) {
                                rows.push((element, value));
                            }
                        }
                        let direct = self.run.stress_pullback_nodal_at(
                            &tape.motion[state][0],
                            &rows,
                            maximum_strain,
                            tape.budget.gradient,
                            cx,
                        )?;
                        for (sum, value) in seed.q.iter_mut().zip(direct.motion.q) {
                            check(cx)?;
                            for c in 0..3 {
                                sum[c] = number(sum[c] + value[c])?;
                            }
                        }
                        add_material(&mut material, direct.material, cx)?;
                    }
                    ReverseVisit3::Residual { state, force } => {
                        let implicit = self.run.material_scale_pullback_nodal_at(
                            &tape.motion[state],
                            force,
                            tape.budget.gradient,
                            cx,
                        )?;
                        add_material(&mut material, implicit, cx)?;
                    }
                }
                Ok(())
            },
            cx,
        )?;
        check(cx)?;
        Ok(PrescribedDesignGradient3 { controls, material })
    }
}
