//! Time-dependent thermal strain on the existing sparse Newmark trajectory.
//!
//! Endpoint constitutive tensors give exact discrete material work:
//! W_T = -sum_e V_e * (sigma_old + sigma_next)/2 : delta(epsilon_free),
//! W_C = sum_e V_e/2 * elastic_strain_old : delta(C) : elastic_strain_next.
//! Both are signed work supplied by a prescribed thermal/material-state path.
//! They include constrained storage changes with zero displacement work;
//! no caloric reservoir or reciprocal thermoelastic heat is inferred here.

use super::*;
use crate::linear3::{TetElasticMaterial, TetMaterialField};

/// The original numerical audit plus the physical shifted elastic-energy audit.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalTransientStep3 {
    pub native: TransientStep3,
    /// Kinetic plus V/2 (epsilon-epsilon_free):C:(epsilon-epsilon_free) [J].
    pub kinetic_plus_elastic_before_j: f64,
    pub kinetic_plus_elastic_after_j: f64,
    /// Original shifted elastic-storage change, evaluated element by element.
    pub elastic_energy_change_j: f64,
    /// Work supplied TO mechanics by the changing stress-free strain [J].
    pub thermal_material_work_j: f64,
    /// Signed work from the changing constitutive tensor at fixed mass [J].
    /// This differs from native unshifted stiffness work when eigenstrain is present.
    pub stiffness_material_work_j: f64,
    /// Delta(K+Psi) + damping - external force/support work - both material work terms [J].
    pub physical_balance_error_j: f64,
}

fn material_at(field: TetMaterialField<'_>, element: usize) -> Option<&TetElasticMaterial> {
    match field {
        TetMaterialField::Uniform(material) => Some(material),
        TetMaterialField::PerElement(materials) => materials.get(element),
    }
}

fn free_strain(field: Option<&response::FrozenThermal3>, element: usize) -> [f64; 6] {
    match field {
        None => [0.0; 6],
        Some(response::FrozenThermal3::Uniform(state)) => state.free_strain_mandel(),
        Some(response::FrozenThermal3::PerElement(states)) => states[element].free_strain_mandel(),
    }
}

impl<'g, 'm> TetTransientRun3<'g, 'm, '_, '_> {
    /// Exact accepted thermal field, including its material/law identities.
    pub fn thermal_strain_field(&self) -> Option<TetThermalStrainField<'_>> {
        self.thermal_strains.as_ref().map(|field| match field {
            response::FrozenThermal3::Uniform(state) => TetThermalStrainField::Uniform(state),
            response::FrozenThermal3::PerElement(states) => TetThermalStrainField::PerElement(states),
        })
    }

    /// Advance the SAME accepted motion using one new absolute thermal state.
    /// The caller retains the refreshed material problem; its complete source
    /// identities are retained by the returned run and native stress recovery.
    /// Reference density stays bit-identical. Explicit nodal inventory remains
    /// authoritative even when the refreshed card's current density changes.
    /// A changed Mandel tensor rebuilds
    /// the native endpoint stiffness; an unchanged tensor reuses its operator.
    /// Changed stiffness uses Cbar = eta*M + zeta*(K_old+K_next)/2 at the
    /// midpoint velocity, keeping Rayleigh loss nonnegative for the admitted K.
    /// Endpoint acceleration is recovered with the actual endpoint operators;
    /// both solves share the original iteration limit. This entry point requires
    /// stationary zero supports; the explicit thermal/support candidates admit
    /// support motion at unchanged stiffness. Mass changes refuse. This remains
    /// the explicitly supplied Rayleigh model, not a newly inferred broadband
    /// or temperature-dependent loss law.
    ///
    /// Eigenstrain follows the endpoint trapezoidal/Newmark convention, as do
    /// external forces. Resolve jumps and short heating pulses explicitly and
    /// check timestep refinement. Both numerical and physical energy gates use
    /// the existing run's tolerances. Refusal changes no accepted state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermal_candidate<'p, 's>(
        &self, problem: &'p TetLinearElasticProblem<'s>, thermal: TetThermalStrainField<'_>,
        load_next: &TetLoadCase3<'_, '_>, h_s: f64, mesh_frame: ContentHash,
        maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), TransientError3> {
        self.advance_thermal_endpoint(
            problem, thermal, load_next, None, h_s, mesh_frame, maximum_strain, cx,
        )
    }

    /// Advance source-refreshed eigenstrain and explicit support displacements
    /// together on the original full consistent-mass equations. The support
    /// convention is exactly `advance_with_supports`: omitted fixed components
    /// have zero endpoint displacement; Newmark supplies their velocity and
    /// acceleration. The accepted support and free states publish together.
    ///
    /// Every element's complete stiffness and reference inertia must remain
    /// unchanged. Source identities and thermal strain may change. Reactions use
    /// the actual endpoint thermal load; support work, external force work and
    /// signed thermal material work remain separate. No heat feedback is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermal_with_supports_candidate<'p, 's>(
        &self, problem: &'p TetLinearElasticProblem<'s>, thermal: TetThermalStrainField<'_>,
        load_next: &TetLoadCase3<'_, '_>, prescribed_dofs_m: &[(usize, f64)],
        h_s: f64, mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), TransientError3> {
        let supports = self.support_endpoint_candidate(
            prescribed_dofs_m, motion::SupportEndpoint::Displacement, h_s, cx,
        )?;
        self.advance_thermal_endpoint(
            problem, thermal, load_next, Some(supports), h_s, mesh_frame, maximum_strain, cx,
        )
    }

    /// The same coupled thermal/support candidate with prescribed endpoint
    /// accelerations [m/s²]. Omitted fixed components have zero acceleration;
    /// their accepted displacement and velocity continue through the original
    /// Newmark map. This retains the direct-acceleration convention without
    /// differencing sampled support displacements. Stiffness and inertia stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermal_with_support_accelerations_candidate<'p, 's>(
        &self, problem: &'p TetLinearElasticProblem<'s>, thermal: TetThermalStrainField<'_>,
        load_next: &TetLoadCase3<'_, '_>, prescribed_dofs_m_s2: &[(usize, f64)],
        h_s: f64, mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), TransientError3> {
        let supports = self.support_endpoint_candidate(
            prescribed_dofs_m_s2, motion::SupportEndpoint::Acceleration, h_s, cx,
        )?;
        self.advance_thermal_endpoint(
            problem, thermal, load_next, Some(supports), h_s, mesh_frame, maximum_strain, cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_thermal_endpoint<'p, 's>(
        &self, problem: &'p TetLinearElasticProblem<'s>, thermal: TetThermalStrainField<'_>,
        load_next: &TetLoadCase3<'_, '_>, mut next_supports: Option<motion::PrescribedMotion3>,
        h_s: f64, mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), TransientError3> {
        check(cx)?;
        if self.prescribed.is_some() != next_supports.is_some()
            || problem.nodes_m != self.problem.nodes_m
            || problem.tetrahedra != self.problem.tetrahedra
            || problem.fixed_dofs != self.problem.fixed_dofs {
            return Err(TransientError3::Invalid("thermal refresh requires the same mesh and an explicit endpoint for prescribed supports"));
        }
        if self.state.step >= self.config.maximum_steps { return Err(TransientError3::Budget); }
        let time = self.state.time_s + h_s;
        if !h_s.is_finite() || h_s <= 0.0 || !time.is_finite() || time <= self.state.time_s {
            return Err(TransientError3::Invalid("thermal endpoint must reach a distinct finite later time"));
        }
        let mut stiffness_changed = false;
        for element in 0..problem.tetrahedra.len() {
            check(cx)?;
            let before = material_at(self.problem.materials, element).ok_or(
                TransientError3::Invalid("missing accepted element material"))?;
            let after = material_at(problem.materials, element).ok_or(
                TransientError3::Invalid("missing refreshed element material"))?;
            if self.nodal_reference_density.is_none()
                && before.density_kg_m3.to_bits() != after.density_kg_m3.to_bits() {
                return Err(TransientError3::Invalid("thermal material refresh cannot change reference density"));
            }
            stiffness_changed |= before.stiffness_mandel_pa().iter().flatten()
                .zip(after.stiffness_mandel_pa().iter().flatten())
                .any(|(a,b)| a.to_bits() != b.to_bits());
        }
        if stiffness_changed && next_supports.is_some() {
            return Err(TransientError3::Invalid(
                "thermal refresh with prescribed supports requires unchanged stiffness",
            ));
        }
        // Existing assembly checks full source identity, counts, constraints,
        // geometry/quality and bounds before any source is retained.
        let load = problem.assemble_thermal_load(thermal, cx)?;
        if load.free_dofs != self.operators.free_dofs {
            return Err(TransientError3::Invalid("thermal refresh changed the original free-DOF map"));
        }
        let operators = if stiffness_changed {
            let refreshed = problem.assemble_prescribed_dynamics(
                self.nodal_reference_density.as_deref(), cx)?;
            if refreshed.free_dofs != self.operators.free_dofs
                || refreshed.element_volumes_m3 != self.operators.element_volumes_m3
                || refreshed.total_mass_kg.to_bits() != self.operators.total_mass_kg.to_bits() {
                return Err(TransientError3::Invalid("stiffness refresh changed the reference inventory or free-DOF map"));
            }
            for row in 0..refreshed.mass.nrows() {
                check(cx)?;
                let (old_columns, old_values) = self.operators.mass.row(row);
                let (new_columns, new_values) = refreshed.mass.row(row);
                if old_columns != new_columns
                    || old_values.iter().zip(new_values).any(|(a,b)| a.to_bits() != b.to_bits())
                    || refreshed.stiffness.row(row).1.iter().any(|v| !v.is_finite()) {
                    return Err(TransientError3::Invalid("stiffness refresh requires identical finite mass and finite stiffness"));
                }
            }
            Arc::new(refreshed)
        } else { self.operators.clone() };
        let before = self.recover_response(mesh_frame, maximum_strain, cx)?;
        let captured = response::FrozenThermal3::capture(thermal, cx)?;
        let mut candidate = TetTransientRun3 {
            geometry: self.geometry, problem, operators,
            nodal_reference_density: self.nodal_reference_density.clone(),
            state: self.state.clone(), thermal: self.thermal.clone(),
            thermal_strains: self.thermal_strains.clone(), force: self.force.clone(),
            full_force: self.full_force.clone(), config: self.config,
            prescribed: self.prescribed.clone(),
        };
        if let Some(supports) = &mut next_supports {
            supports.thermal = Arc::from(load.full_equivalent_force_n);
        }
        let native = candidate.advance_with_thermal_impl(
            load_next, next_supports, Some(&load.reduced_force_n),
            stiffness_changed.then_some(&self.operators.stiffness), h_s, cx)?;
        candidate.thermal_strains = Some(captured);
        let after = candidate.recover_response(mesh_frame, maximum_strain, cx)?;
        let mut work = Dd::default();
        let mut stiffness_work = Dd::default();
        let mut elastic_change = Dd::default();
        for (element, &volume) in self.operators.element_volumes_m3.iter().enumerate() {
            check(cx)?;
            let old_strain = free_strain(self.thermal_strains.as_ref(), element);
            let next_strain = free_strain(candidate.thermal_strains.as_ref(), element);
            let old_stress = before.stress.elements()[element].stress_mandel_pa();
            let next_stress = after.stress.elements()[element].stress_mandel_pa();
            for component in 0..6 {
                let average = Dd::from_f64(0.5 * old_stress[component])
                    + Dd::from_f64(0.5 * next_stress[component]);
                let delta = Dd::from_f64(next_strain[component]) - Dd::from_f64(old_strain[component]);
                work = work - Dd::from_f64(volume) * average * delta;
            }
            let old_material = material_at(self.problem.materials, element).ok_or(
                TransientError3::Invalid("missing accepted element material"))?;
            let next_material = material_at(problem.materials, element).ok_or(
                TransientError3::Invalid("missing refreshed element material"))?;
            let old_elastic = before.stress.elements()[element].elastic_strain_mandel();
            let next_elastic = after.stress.elements()[element].elastic_strain_mandel();
            // Evaluate the original quadratic independently of stress work:
            // delta(Psi) = e_bar' C0 delta(e) + e1' delta(C) e1 / 2.
            // Coefficient differences are formed before common strain products.
            let mut fixed_leading = Dd::default();
            let mut fixed_tail = Dd::default();
            let mut changed = Dd::default();
            for i in 0..6 { for j in 0..6 {
                let mean = Dd::from_f64(0.5 * old_elastic[i]) + Dd::from_f64(0.5 * next_elastic[i]);
                let delta = Dd::from_f64(next_elastic[j]) - Dd::from_f64(old_elastic[j]);
                let coefficient = Dd::from_f64(old_material.stiffness_mandel_pa()[i][j]);
                fixed_leading = fixed_leading + Dd::from_f64(mean.hi) * coefficient * Dd::from_f64(delta.hi);
                fixed_tail = fixed_tail + Dd::from_f64(mean.lo) * coefficient * delta
                    + Dd::from_f64(mean.hi) * coefficient * Dd::from_f64(delta.lo);
                let difference = Dd::from_f64(next_material.stiffness_mandel_pa()[i][j]) - coefficient;
                changed = changed + Dd::from_f64(0.5 * next_elastic[i]) * difference
                    * Dd::from_f64(next_elastic[j]);
            }}
            elastic_change = elastic_change + Dd::from_f64(volume) * (fixed_leading + fixed_tail + changed);
            if stiffness_changed {
                for i in 0..6 { for j in 0..6 {
                    let delta = Dd::from_f64(next_material.stiffness_mandel_pa()[i][j])
                        - Dd::from_f64(old_material.stiffness_mandel_pa()[i][j]);
                    stiffness_work = stiffness_work + Dd::from_f64(0.5 * volume)
                        * Dd::from_f64(old_elastic[i]) * delta * Dd::from_f64(next_elastic[j]);
                }}
            }
        }
        let old_energy = before.kinetic_plus_elastic_energy_j;
        let next_energy = after.kinetic_plus_elastic_energy_j;
        let elastic_change = energy::retained(elastic_change)?;
        let work = energy::retained(work)?;
        let stiffness_work = energy::retained(stiffness_work)?;
        let defect = energy::retained(native.retained_kinetic_energy_change_j() + elastic_change
            + native.retained_damping_loss_j() - native.retained_external_work_j()
            - work - stiffness_work)?.to_f64();
        let scale = native.kinetic_energy_change_j.abs().max(elastic_change.to_f64().abs())
            .max(work.to_f64().abs()).max(stiffness_work.to_f64().abs())
            .max(native.damping_loss_j.abs()).max(native.mechanical_work_j.abs())
            .max(native.support_work_j.abs());
        let allowance = number(self.config.absolute_energy_tolerance_j
            + self.config.relative_energy_tolerance*scale)?;
        if defect.abs() > allowance {
            return Err(TransientError3::Energy { defect_j: defect, allowance_j: allowance });
        }
        check(cx)?;
        Ok((candidate, ThermalTransientStep3 { native,
            kinetic_plus_elastic_before_j: old_energy,
            kinetic_plus_elastic_after_j: next_energy,
            elastic_energy_change_j: elastic_change.to_f64(),
            thermal_material_work_j: work.to_f64(), stiffness_material_work_j: stiffness_work.to_f64(),
            physical_balance_error_j: defect }))
    }
}
