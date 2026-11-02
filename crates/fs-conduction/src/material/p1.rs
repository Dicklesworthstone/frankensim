//! Constitutive average and consistent nodal derivatives for a P1 tetrahedron.

use fs_exec::Cx;

use super::{ConductivityKind, ConductivityModel, ConductivityTable, TemperatureSpan};
use crate::ConductionError;

mod moments;

pub(crate) struct ElementConductivity {
    /// Volume-normalized integral of K(T_h).
    pub(crate) tensor: [[f64; 3]; 3],
    /// Volume-normalized integral of lambda_c K'(T_h), in local vertex order.
    pub(crate) derivatives: [[[f64; 3]; 3]; 4],
}

fn scalar(
    cx: Option<&Cx<'_>>,
    table: &ConductivityTable,
    nodal: [f64; 4],
    with_derivatives: bool,
) -> Result<moments::ScalarMoments, ConductionError> {
    let (average, derivatives) = integrate_scalar_path(cx, table.knots(), nodal, with_derivatives)?;
    Ok(moments::ScalarMoments {
        average,
        derivatives,
    })
}

/// Integrate a positive piecewise-linear scalar law over a P1 coordinate.
/// The coordinate may be temperature or specific enthalpy; units and complete
/// source support are owned by the caller. No temperature table is synthesized.
pub(crate) fn integrate_scalar_path(
    cx: Option<&Cx<'_>>,
    knots: &[(f64, f64)],
    nodal: [f64; 4],
    with_derivatives: bool,
) -> Result<(f64, [f64; 4]), ConductionError> {
    moments::integrate_table(nodal, knots, with_derivatives, |interval| {
        if let Some(cx) = cx {
            cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                stage: "integrate-conductivity-knots",
                at: interval,
            })?;
        }
        Ok(())
    })
    .map(|moments| (moments.average, moments.derivatives))
    .map_err(|error| match error {
        moments::MomentError::Interrupted(error) => error,
        moments::MomentError::Numerical(what) => {
            ConductionError::Conductivity { what: what.into() }
        }
    })
}

impl ConductivityModel {
    /// Integrate the admitted scalar components, then rotate once into the
    /// original fixed principal frame. Every nodal value must be supported;
    /// an admissible mean cannot hide unsupported material elsewhere in a tet.
    /// The scalar rule is exact for the retained piecewise-linear table (up to
    /// roundoff), including any crossed source knots. Tangents are optional so
    /// operator/flux assembly never computes four unused derivative moments.
    pub(crate) fn p1_conductivity(
        &self,
        cx: Option<&Cx<'_>>,
        nodal: [f64; 4],
        with_derivatives: bool,
    ) -> Result<ElementConductivity, ConductionError> {
        let span = self.temperature_span();
        for t in nodal {
            crate::require_finite("temperature iterate", t)?;
            span.check(t)?;
        }
        let mut result = ElementConductivity {
            tensor: [[0.0; 3]; 3],
            derivatives: [[[0.0; 3]; 3]; 4],
        };
        match &self.kind {
            ConductivityKind::ConstantTensor { tensor, .. } => result.tensor = **tensor,
            ConductivityKind::Isotropic(table) => {
                let values = scalar(cx, table, nodal, with_derivatives)?;
                for i in 0..3 {
                    result.tensor[i][i] = values.average;
                    if with_derivatives {
                        for c in 0..4 {
                            result.derivatives[c][i][i] = values.derivatives[c];
                        }
                    }
                }
            }
            ConductivityKind::Orthotropic { axes, tables } => {
                for (axis, table) in axes.iter().zip(tables.iter()) {
                    let values = scalar(cx, table, nodal, with_derivatives)?;
                    for i in 0..3 {
                        for j in 0..3 {
                            result.tensor[i][j] =
                                (values.average * axis[i]).mul_add(axis[j], result.tensor[i][j]);
                            if with_derivatives {
                                for c in 0..4 {
                                    result.derivatives[c][i][j] = (values.derivatives[c] * axis[i])
                                        .mul_add(axis[j], result.derivatives[c][i][j]);
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    /// Interior knots have zero volume in a nonuniform P1 field, so crossing
    /// one does not destroy differentiability of its material integral. A
    /// uniform element exactly at a true kink has only a semismooth tangent;
    /// source endpoints also prohibit an unrestricted two-sided perturbation.
    pub(crate) fn p1_is_smooth(&self, nodal: [f64; 4]) -> Result<bool, ConductionError> {
        let span = self.temperature_span();
        for t in nodal {
            crate::require_finite("temperature iterate", t)?;
            span.check(t)?;
        }
        if !self.is_temperature_dependent() {
            return Ok(true);
        }
        if let TemperatureSpan::Sampled { low, high } = span {
            if nodal.iter().any(|&t| t == low || t == high) {
                return Ok(false);
            }
        }
        if nodal.iter().any(|&t| t != nodal[0]) {
            return Ok(true);
        }
        let slope = self.tensor_derivative_at(nodal[0])?;
        for row in slope {
            for value in row {
                crate::require_finite("conductivity derivative", value)?;
            }
        }
        for neighbour in [nodal[0].next_down(), nodal[0].next_up()] {
            match self.tensor_derivative_at(neighbour) {
                Ok(other) => {
                    for row in other {
                        for value in row {
                            crate::require_finite("conductivity derivative", value)?;
                        }
                    }
                    if slope != other {
                        return Ok(false);
                    }
                }
                Err(ConductionError::OutsideTemperatureSpan { .. }) => return Ok(false),
                Err(error) => return Err(error),
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests;
