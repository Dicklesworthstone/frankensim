//! Discrete force balance before strong homogeneous support rows are eliminated.
use super::surface::ReferenceLoad3;
use super::*;

/// Physical force balance of one supplied displacement and complete variational RHS.
/// No equilibrium solve or continuum traction recovery is performed.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumResponse3 {
    /// `K_full u - rhs_full`, in the operator's node-major xyz ordering.
    /// Fixed entries are forces exerted ON the solid by strong supports; free
    /// entries are the remaining equilibrium imbalance, not support reactions.
    /// Includes the original ghost and embedded Nitsche terms, when present.
    pub residual: Vec<f64>,
    /// Euclidean norm of free entries only, in force units. This is a measured
    /// absolute imbalance, not a solver tolerance or continuum-error certificate.
    pub free_residual_norm: f64,
    /// Sum of strong-support forces only. Embedded reactions remain available
    /// separately through `embedded_reaction` and are not counted twice here.
    pub support_force: [f64; 3],
    /// Sum `(node - moment_origin) cross reaction` at strong supports.
    pub support_moment: [f64; 3],
    /// Caller-selected reference point for the reported support moment.
    pub moment_origin: [f64; 3],
}

pub(super) fn poll(
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

pub(super) fn admit_load(
    rhs: &[f64],
    n: usize,
    origin: [f64; 3],
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    poll(checkpoint)?;
    if rhs.len() != n || origin.iter().any(|v| !v.is_finite()) {
        return Err(ElasticityError3::Invalid(
            "invalid equilibrium load or moment origin",
        ));
    }
    for chunk in rhs.chunks(192) {
        poll(checkpoint)?;
        if chunk.iter().any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid("nonfinite equilibrium load"));
        }
    }
    Ok(())
}

pub(super) fn response(
    mut force: Vec<f64>,
    rhs: &[f64],
    nodes: &[[f64; 3]],
    fixed: &[bool],
    origin: [f64; 3],
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<EquilibriumResponse3, ElasticityError3> {
    let mut free_residual_norm = 0.0_f64;
    let mut support_force = [0.0; 3];
    let mut support_moment = [0.0; 3];
    for (node, (&position, &clamped)) in nodes.iter().zip(fixed).enumerate() {
        poll(checkpoint)?;
        let dof = 3 * node;
        for component in 0..3 {
            force[dof + component] -= rhs[dof + component];
            if !force[dof + component].is_finite() {
                return Err(ElasticityError3::Invalid("equilibrium residual overflow"));
            }
        }
        if clamped {
            let arm: [f64; 3] = std::array::from_fn(|i| position[i] - origin[i]);
            for component in 0..3 {
                support_force[component] += force[dof + component];
                support_moment[component] += arm[(component + 1) % 3]
                    * force[dof + (component + 2) % 3]
                    - arm[(component + 2) % 3] * force[dof + (component + 1) % 3];
            }
        } else {
            for &value in &force[dof..dof + 3] {
                free_residual_norm = free_residual_norm.hypot(value);
            }
        }
    }
    if !free_residual_norm.is_finite()
        || support_force
            .iter()
            .chain(&support_moment)
            .any(|v| !v.is_finite())
    {
        return Err(ElasticityError3::Invalid(
            "support force or moment overflow",
        ));
    }
    poll(checkpoint)?;
    Ok(EquilibriumResponse3 {
        residual: force,
        free_residual_norm,
        support_force,
        support_moment,
        moment_origin: origin,
    })
}

impl CutElasticity3 {
    /// Assemble body/surface forces and optional embedded displacement lifting,
    /// retaining entries at strong supports. The same vector can be passed to
    /// `solve_physical_controlled`, which masks supported rows only for solving,
    /// and then to `equilibrium_response` without losing supported loading.
    ///
    /// Pure reference laws and the original retained integration rules are used.
    /// Add any independent nodal forces to this COMPLETE vector, including forces
    /// at supports. Reassemble after changing scales when `prescribed` is present.
    /// The lifting is variational data, not applied traction or actuator work.
    /// Existing masked load vectors, including prepared material-body loads, do
    /// not contain enough information to reconstruct their supported entries.
    pub fn equilibrium_load(
        &self,
        load: ReferenceLoad3<'_>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if prescribed.is_some() && self.embedded.is_none() {
            return Err(ElasticityError3::Invalid(
                "motion requires embedded Dirichlet support",
            ));
        }
        let mut rhs = self.reference_load_impl(load, true, &mut checkpoint)?;
        if let Some(g) = prescribed {
            let lifting = self.prescribed_displacement_load_impl(g, true, &mut checkpoint)?;
            for (i, (force, lift)) in rhs.iter_mut().zip(lifting).enumerate() {
                if i % 192 == 0 {
                    poll(&mut checkpoint)?;
                }
                *force += lift;
                if !force.is_finite() {
                    return Err(ElasticityError3::Invalid("equilibrium load overflow"));
                }
            }
        }
        poll(&mut checkpoint)?;
        Ok(rhs)
    }

    /// Recover the full discrete residual and strong-support force/moment.
    /// `full_rhs` must retain supported external forces AND any embedded lifting;
    /// use `equilibrium_load`, not the masked `body_load`/`reference_load` paths.
    /// Supports impose zero displacement exactly as in the solver. Original
    /// material strain domains are checked throughout every active Q1 cell.
    ///
    /// The supplied field need not be equilibrated: `free_residual_norm` exposes
    /// its imbalance. For a solved field, fixed residual entries are reactions.
    /// Current cell scales, B-bar stiffness, Nitsche terms and ghost forces are
    /// retained. These discrete constraint forces are not raw sampled stress.
    /// Cancellation or invalid arithmetic returns no partial response.
    pub fn equilibrium_response(
        &self,
        u: &[f64],
        full_rhs: &[f64],
        moment_origin: [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EquilibriumResponse3, ElasticityError3> {
        admit_load(full_rhs, self.n(), moment_origin, &mut checkpoint)?;
        self.admit_material_domain(u, &mut checkpoint)?;
        let force = self.uneliminated_force(
            self.n(),
            |node, component| {
                if self.fixed[node] {
                    0.0
                } else {
                    u[3 * node + component]
                }
            },
            |node, component, value, out| out[3 * node + component] += value,
            &mut checkpoint,
        )?;
        response(
            force,
            full_rhs,
            &self.nodes,
            &self.fixed,
            moment_origin,
            &mut checkpoint,
        )
    }

    // Share the original retained element/ghost action across Cartesian and
    // adaptive response. Only the output vector is allocated. Local reconstruction
    // masks input supports; the force transpose deliberately retains ALL rows.
    pub(super) fn uneliminated_force(
        &self,
        n: usize,
        coefficient: impl Fn(usize, usize) -> f64,
        scatter: impl Fn(usize, usize, f64, &mut [f64]),
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        let mut force = vec![0.0; n];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(checkpoint)?;
            let local: [f64; 24] = std::array::from_fn(|j| coefficient(cell.nodes[j / 3], j % 3));
            for i in 0..24 {
                let value: f64 = cell.stiffness[i]
                    .iter()
                    .zip(local)
                    .map(|(k, u)| k * u)
                    .sum();
                scatter(
                    cell.nodes[i / 3],
                    i % 3,
                    self.scales[id] * value,
                    &mut force,
                );
            }
        }
        for face in &self.ghosts {
            poll(checkpoint)?;
            let weight = face.scaled_weight(&self.scales);
            for component in 0..3 {
                let jump: f64 = face
                    .nodes
                    .iter()
                    .zip(&face.jump)
                    .map(|(&node, &j)| j * coefficient(node, component))
                    .sum();
                for (&node, &j) in face.nodes.iter().zip(&face.jump) {
                    scatter(node, component, weight * j * jump, &mut force);
                }
            }
        }
        poll(checkpoint)?;
        Ok(force)
    }
}
