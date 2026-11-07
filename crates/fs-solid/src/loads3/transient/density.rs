//! Spatial reference inventory shared by every native transient operator.
use super::*;

fn capture(
    problem: &TetLinearElasticProblem<'_>,
    density: &[f64],
    cx: &Cx<'_>,
) -> Result<Arc<[f64]>, TransientError3> {
    check(cx)?;
    if density.len() != problem.nodes_m.len()
        || density.len() > problem.budget.maximum_nodes
    {
        return Err(TransientError3::Invalid(
            "reference density requires one sample per admitted mesh node",
        ));
    }
    let mut retained = Vec::with_capacity(density.len());
    for &value in density {
        check(cx)?;
        if !value.is_finite() || value <= 0.0 {
            return Err(TransientError3::Invalid(
                "reference density must be finite and positive",
            ));
        }
        retained.push(value);
    }
    check(cx)?;
    Ok(Arc::from(retained))
}

impl<'g, 'm> TetLoadCase3<'g, 'm> {
    /// Initialize with a positive P1 density [kg/m³] on the ORIGINAL mesh.
    /// The exact cubic integral of rho_0 N_a N_b supplies consistent inertia;
    /// material cards still supply stiffness. The field is copied once, shared
    /// by forks/checkpoints, and retained through material/thermal refreshes.
    /// It is not resampled at the current temperature or deformed geometry.
    ///
    /// Full reactions, impulses, mass damping, adjoints and acceleration-history
    /// forcing use this same inventory. Initial displacements and velocities
    /// have the stationary-support semantics of `begin_transient`.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient_with_nodal_reference_density<'p, 's>(
        &self,
        problem: &'p TetLinearElasticProblem<'s>,
        reference_density_kg_m3: &[f64],
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        config: TransientConfig3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        config.validate()?;
        let density = capture(problem, reference_density_kg_m3, cx)?;
        self.begin_transient_impl(problem, thermal, time_s, displacement_m,
            velocity_m_s, None, None, Some(density), config, cx)
    }

    /// Graded reference inertia with the full free/fixed mass coupling of
    /// `begin_transient_with_supports`. Prescribed acceleration is full-node
    /// xyz and zero on free components; subsequent support APIs are unchanged.
    /// Fully prescribed solids retain their full physical mass and support work.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient_with_nodal_reference_density_and_supports<'p, 's>(
        &self,
        problem: &'p TetLinearElasticProblem<'s>,
        reference_density_kg_m3: &[f64],
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: &[[f64; 3]],
        config: TransientConfig3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        config.validate()?;
        let density = capture(problem, reference_density_kg_m3, cx)?;
        self.begin_transient_impl(problem, thermal, time_s, displacement_m,
            velocity_m_s, Some(prescribed_acceleration_m_s2), None,
            Some(density), config, cx)
    }
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Explicit reference field, if supplied. `None` uses cellwise material
    /// density. The returned immutable values are the original mass inventory.
    pub fn nodal_reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.nodal_reference_density.as_deref()
    }

    /// Physical mass [kg], including every constrained component's share.
    pub fn reference_mass_kg(&self) -> f64 {
        self.operators.total_mass_kg
    }
}
