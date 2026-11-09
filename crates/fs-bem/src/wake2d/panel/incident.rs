//! One prescribed external field shared by boundary solves and wake transport.

use super::{BemError, checkpoint};

type Variation<'a> = dyn FnMut([f64; 2], f64) -> Result<[f64; 2], BemError> + 'a;

pub(super) struct IncidentField<'a> {
    variation: Option<&'a mut Variation<'a>>,
}

impl<'a> IncidentField<'a> {
    pub(super) fn uniform() -> Self {
        Self { variation: None }
    }

    pub(super) fn new(variation: &'a mut Variation<'a>) -> Self {
        Self {
            variation: Some(variation),
        }
    }

    pub(super) fn velocity(
        &mut self,
        reference: [f64; 2],
        point: [f64; 2],
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], BemError> {
        let Some(variation) = self.variation.as_mut() else {
            return Ok(reference);
        };
        if !time.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "panel wake incident time",
            });
        }
        checkpoint(cancelled)?;
        let change = variation(point, time)?;
        checkpoint(cancelled)?;
        if change.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "panel wake incident variation",
            });
        }
        // Preserve the original operation tree, including signed zero, when
        // a caller supplies an identically zero variation.
        let velocity = core::array::from_fn(|axis| {
            if change[axis] == 0.0 {
                reference[axis]
            } else {
                reference[axis] + change[axis]
            }
        });
        if velocity.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "panel wake total incident velocity",
            });
        }
        Ok(velocity)
    }
}
