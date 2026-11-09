//! One prescribed external field shared by boundary solves and wake transport.

use super::motion::{MovingField, PanelFrame2d};
use super::{BemError, checkpoint};

type Variation<'a> = dyn FnMut([f64; 2], f64) -> Result<[f64; 2], BemError> + 'a;

pub(super) struct IncidentField<'a> {
    variation: Option<&'a mut Variation<'a>>,
    moving: Option<MovingField<'a>>,
}

impl<'a> IncidentField<'a> {
    pub(super) fn uniform() -> Self {
        Self {
            variation: None,
            moving: None,
        }
    }

    pub(super) fn new(variation: &'a mut Variation<'a>) -> Self {
        Self {
            variation: Some(variation),
            moving: None,
        }
    }

    pub(super) fn moving(
        motion: &'a mut dyn FnMut(f64) -> Result<PanelFrame2d, BemError>,
        air: &'a mut Variation<'a>,
    ) -> Self {
        Self {
            variation: None,
            moving: Some(MovingField::new(motion, air)),
        }
    }

    pub(super) fn begin(
        &mut self,
        time: f64,
        frame: Option<PanelFrame2d>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), BemError> {
        match (&mut self.moving, frame) {
            (None, None) => Ok(()),
            (Some(moving), Some(frame)) => moving.begin(time, frame, cancelled),
            _ => Err(BemError::InvalidScalar {
                name: "panel wake frame mode",
                value: if frame.is_some() { 1.0 } else { 0.0 },
                requirement: "use the stepping API corresponding to the constructed frame mode",
            }),
        }
    }

    pub(super) fn accepted_frame(&self) -> Option<PanelFrame2d> {
        self.moving.as_ref().and_then(MovingField::accepted_frame)
    }

    pub(super) fn velocity(
        &mut self,
        reference: [f64; 2],
        point: [f64; 2],
        time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], BemError> {
        if let Some(moving) = &mut self.moving {
            return moving.velocity(point, time, cancelled);
        }
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
