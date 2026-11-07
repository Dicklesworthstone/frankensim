//! Parameter contractions of the original P1 stiffness, consistent mass and B.
use super::super::*;

impl TetLinearElasticProblem<'_> {
    /// -lambda' d(M a + (eta M + zeta K) v + K u - fT).
    /// Columns are independent multiplicative stiffness/density scales at one;
    /// the final pair differentiates eta and zeta themselves. Thermal eigenstrain
    /// is fixed, so the stiffness derivative includes its equivalent-load term.
    /// The accepted transient owner embeds free equilibrium multipliers here.
    /// Reaction observations instead use negative constrained-row weights to
    /// obtain their direct positive residual/material contraction.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn transient_material_contraction(
        &self,
        motion: &[Vec<[f64; 3]>; 3],
        lambda: &[[f64; 3]],
        thermal: Option<TetThermalStrainField<'_>>,
        eta: f64,
        zeta: f64,
        nodal_reference_density: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, [f64; 2]), TetElasticError>{
        let invalid = |what| TetElasticError::InvalidStressObservation { what };
        let finite = |x: f64| {
            if x.is_finite() {
                Ok(x)
            } else {
                Err(invalid("nonfinite transient material-gradient arithmetic"))
            }
        };
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        if let Some(field) = thermal {
            field.validate(self.tetrahedra.len())?;
        }
        if motion.iter().any(|v| v.len() != self.nodes_m.len())
            || lambda.len() != self.nodes_m.len()
            || [eta, zeta].iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid("invalid material-gradient motion, seed or damping"));
        }
        for field in motion
            .iter()
            .map(Vec::as_slice)
            .chain(std::iter::once(lambda))
        {
            for node in field {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                for &value in node {
                    finite(value)?;
                }
            }
        }
        let mut elements = Vec::with_capacity(self.tetrahedra.len());
        let mut damping = [0.0; 2];
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let material = self.materials.at(element);
            let free = match thermal {
                Some(field) => {
                    let state = field.at(element);
                    if state.elastic_material_state_identity != material.material_state_identity {
                        return Err(TetElasticError::ThermalElasticStateMismatch { element });
                    }
                    state.free_strain_mandel
                }
                None => [0.0; 6],
            };
            let me = p1_tet_reference_mass(geometry.volume, material.density_kg_m3,
                *tet, nodal_reference_density)?;
            let b = strain_displacement_matrix(&geometry);
            let strain = |field: &[[f64; 3]]| {
                relative_mandel_strain(
                    tet.len(),
                    |d| core::array::from_fn(|i| b[i][d]),
                    |d| field[tet[d / 3]][d % 3],
                )
            };
            let total = strain(&motion[0]);
            let rate = strain(&motion[1]);
            let elastic: [f64; 6] = core::array::from_fn(|i| total[i] - free[i]);
            let stress = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(elastic[j], sum)
                })
            });
            let rate_stress = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(rate[j], sum)
                })
            });
            let force = relative_mandel_force(&b, &stress, geometry.volume);
            let kv = relative_mandel_force(&b, &rate_stress, geometry.volume);
            let mut scales = [0.0; 2];
            for d in 0..12 {
                let seed = lambda[tet[d / 3]][d % 3];
                let mut ma = 0.0;
                let mut mv = 0.0;
                for j in 0..12 {
                    let node = tet[j / 3];
                    let c = j % 3;
                    ma = me[d][j].mul_add(motion[2][node][c], ma);
                    mv = me[d][j].mul_add(motion[1][node][c], mv);
                }
                let mass = finite(finite(ma)? + eta * finite(mv)?)?;
                scales[1] = finite((-seed).mul_add(mass, scales[1]))?;
                damping[0] = finite((-seed).mul_add(mv, damping[0]))?;
                if d >= 3 {
                    // The anchor receives minus the other forces. Contract
                    // that transpose with a relative seed before multiplying,
                    // so a large common cotangent cannot create elastic work.
                    let relative_seed = finite(seed - lambda[tet[0]][d % 3])?;
                    let stiffness = finite(finite(force[d])? + zeta * finite(kv[d])?)?;
                    scales[0] = finite((-relative_seed).mul_add(stiffness, scales[0]))?;
                    damping[1] = finite((-relative_seed).mul_add(kv[d], damping[1]))?;
                }
            }
            elements.push(scales);
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok((elements, damping))
    }
}
