//! Pointwise P1 stress observations on the original B and material tensors.
use super::super::*;

impl TetLinearElasticProblem<'_> {
    /// Transpose sigma=C(Bq-epsilon_free) in mesh-frame Mandel coordinates
    /// [xx,yy,zz,xy,yz,zx]. A pointwise stress seed has NO element-volume factor.
    /// Geometry and free strain are fixed; stiffness scales multiply each C.
    pub(crate) fn stress_observation_contraction(
        &self,
        displacement: &[[f64; 3]],
        rows: &[(usize, [f64; 6])],
        thermal: Option<TetThermalStrainField<'_>>,
        maximum_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, Vec<f64>), TetElasticError> {
        let invalid = |what| TetElasticError::InvalidStressObservation { what };
        let finite = |value: f64| {
            if value.is_finite() {
                Ok(value)
            } else {
                Err(invalid("nonfinite stress-observation derivative"))
            }
        };
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        if let Some(field) = thermal {
            field.validate(self.tetrahedra.len())?;
        }
        if displacement.len() != self.nodes_m.len()
            || !maximum_strain.is_finite()
            || maximum_strain <= 0.0
        {
            return Err(invalid(
                "full displacement and finite positive strain ceiling required",
            ));
        }
        for node in displacement {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for &value in node {
                finite(value)?;
            }
        }
        let mut previous = None;
        for &(element, weights) in rows {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            if element >= self.tetrahedra.len() || previous.is_some_and(|p| p >= element) {
                return Err(invalid(
                    "stress rows must name sorted unique existing elements",
                ));
            }
            for value in weights {
                finite(value)?;
            }
            previous = Some(element);
        }
        let mut displacement_bar = vec![[0.0; 3]; self.nodes_m.len()];
        let mut stiffness_bar = vec![0.0; self.tetrahedra.len()];
        let mut row = 0;
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|n| self.nodes_m[n]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let b = strain_displacement_matrix(&geometry);
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
            let total: [f64; 6] = core::array::from_fn(|i| {
                (0..12).fold(0.0, |sum, d| {
                    b[i][d].mul_add(displacement[tet[d / 3]][d % 3], sum)
                })
            });
            let elastic: [f64; 6] = core::array::from_fn(|i| total[i] - free[i]);
            if total
                .iter()
                .chain(&elastic)
                .any(|v| !v.is_finite() || v.abs() > maximum_strain)
            {
                return Err(invalid(
                    "stress pullback exceeds the finite total or elastic small-strain ceiling",
                ));
            }
            if row >= rows.len() || rows[row].0 != element {
                continue;
            }
            let weights = rows[row].1;
            row += 1;
            let stress: [f64; 6] = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(elastic[j], sum)
                })
            });
            // Transpose the actual canonical entries, preserving every
            // normal/shear coupling in an admitted anisotropic tensor.
            let elastic_bar: [f64; 6] = core::array::from_fn(|j| {
                (0..6).fold(0.0, |sum, i| {
                    material.stiffness_mandel_pa[i][j].mul_add(weights[i], sum)
                })
            });
            for &value in stress.iter().chain(&elastic_bar) {
                finite(value)?;
            }
            stiffness_bar[element] =
                finite((0..6).fold(0.0, |sum, i| weights[i].mul_add(stress[i], sum)))?;
            for d in 0..12 {
                let value =
                    finite((0..6).fold(0.0, |sum, i| b[i][d].mul_add(elastic_bar[i], sum)))?;
                let accumulated = &mut displacement_bar[tet[d / 3]][d % 3];
                *accumulated = finite(*accumulated + value)?;
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok((displacement_bar, stiffness_bar))
    }
}
