//! Thermal strain at the same integration points as the quadratic solid.

use crate::linear3::{TetElasticError, TetThermalStrainField, TetThermalStrainState};

/// Integrated thermal eigenstrain on TET10, HEX20 or HEX27 elements.
///
/// The elementwise variant retains the original constant-strain fast path.
/// Pointwise input resolves thermal gradients without replacing them by an
/// element mean: loads, stored energy, stress and initial-stress stiffness all
/// consume the SAME samples. Elastic tensors remain constant per element.
/// This field prescribes thermal loading; it does not solve heat transport or
/// infer an expansion law from temperatures.
#[derive(Debug, Clone, Copy)]
pub enum QuadraticThermalStrainField<'a> {
    /// Original uniform or per-element thermal states.
    Elementwise(TetThermalStrainField<'a>),
    /// One state per stiffness quadrature point, in element-major order.
    ///
    /// Within each element, use `duffy_quadrature` order for TET10 and
    /// `HexQuadrature::points` order for HEX20/HEX27. Samples must be recomputed
    /// when changing the integration rule. Each state must reference its
    /// element's exact elastic material identity and already contain the
    /// source-integrated mesh-frame Mandel strain at that point.
    AtQuadrature {
        /// Flat array with `element_count * points_per_element` states.
        states: &'a [TetThermalStrainState],
        /// Exact number of stiffness points in the selected integration rule.
        points_per_element: usize,
    },
}

impl<'a> From<TetThermalStrainField<'a>> for QuadraticThermalStrainField<'a> {
    fn from(field: TetThermalStrainField<'a>) -> Self {
        Self::Elementwise(field)
    }
}

impl QuadraticThermalStrainField<'_> {
    pub(super) fn validate(
        &self,
        element_count: usize,
        quadrature_points: usize,
    ) -> Result<(), TetElasticError> {
        match self {
            Self::Elementwise(field) => field.validate(element_count),
            Self::AtQuadrature {
                states,
                points_per_element,
            } => {
                if *points_per_element != quadrature_points || quadrature_points == 0 {
                    return Err(TetElasticError::InvalidThermalStrain {
                        what: "thermal samples must use the exact stiffness quadrature rule",
                    });
                }
                let expected = element_count
                    .checked_mul(quadrature_points)
                    .ok_or(TetElasticError::SizeOverflow)?;
                if states.len() != expected {
                    return Err(TetElasticError::ThermalStrainCountMismatch {
                        expected,
                        got: states.len(),
                    });
                }
                Ok(())
            }
        }
    }

    pub(super) fn at(&self, element: usize, point: usize) -> &TetThermalStrainState {
        match self {
            Self::Elementwise(field) => field.at(element),
            Self::AtQuadrature {
                states,
                points_per_element,
            } => &states[element * points_per_element + point],
        }
    }

    pub(super) const fn is_elementwise(&self) -> bool {
        matches!(self, Self::Elementwise(_))
    }
}
