//! Exact commuting subflows for spherical and aligned axisymmetric gyrostats.

use super::{LieStepError, So3, So3Tangent, Vec3};

/// Inputs have already passed the original gyrostat admission. Only exact
/// inertia equalities and exact axis alignment select a symmetric-top flow.
pub(super) fn step(
    rotation: So3,
    momentum: Vec3,
    inertia: Vec3,
    rotor: Vec3,
    duration: f64,
) -> Result<Option<(So3, Vec3)>, LieStepError> {
    if momentum == rotor {
        return Ok(Some((rotation, momentum)));
    }
    let moments = [inertia.x, inertia.y, inertia.z];
    let stored = [rotor.x, rotor.y, rotor.z];
    let total = [momentum.x, momentum.y, momentum.z];
    let mut selected = None;
    for axis in 0..3 {
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        if moments[a] == moments[b] && stored[a] == 0.0 && stored[b] == 0.0 {
            let carrier_rate = (total[axis] - stored[axis]) / moments[axis];
            // Pure axial motion must not subtract two large rotor phases.
            if total[a] == 0.0 && total[b] == 0.0 {
                let mut tangent = [0.0; 3];
                tangent[axis] = duration * carrier_rate;
                let increment = So3::exp(So3Tangent::new(Vec3::new(
                    tangent[0], tangent[1], tangent[2],
                )))?;
                return Ok(Some((rotation.compose(increment)?, momentum)));
            }
            let mut residual = [0.0; 3];
            residual[axis] = carrier_rate - total[axis] / moments[a];
            selected = Some((moments[a], Vec3::new(residual[0], residual[1], residual[2])));
            break;
        }
    }
    let (reference, residual) = match selected {
        Some(value) => value,
        None if moments[0] == moments[1] && moments[1] == moments[2] => (
            moments[0],
            Vec3::new(
                -rotor.x / moments[0],
                -rotor.y / moments[0],
                -rotor.z / moments[0],
            ),
        ),
        None => return Ok(None),
    };

    // H=|L|²/(2A)+H_residual. The spherical Casimir flow commutes with
    // the residual. Here that residual is one exact constant-rate rotation:
    // about the symmetry axis, or about -h/A for a spherical carrier.
    let increment = So3::exp(So3Tangent::new(residual.scale(duration)))?;
    let next_momentum = increment.inverse().rotate(momentum)?;
    let next_rotation = super::super::so3_body_exp_step(
        rotation.compose(increment)?,
        Vec3::new(
            next_momentum.x / reference,
            next_momentum.y / reference,
            next_momentum.z / reference,
        ),
        duration,
    )?;
    Ok(Some((next_rotation, next_momentum)))
}
