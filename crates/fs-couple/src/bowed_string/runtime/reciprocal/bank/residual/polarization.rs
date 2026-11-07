//! Fixed point contacts through both independently condensed bridge axes.
//! The bank admits only orthogonal bridge axes, so their omitted-mode Green
//! matrix is diagonal. Sum each axis's Schur term without changing the retained
//! modes, contact force law or time integrator.

use super::super::polarization::StringBodyTranslationConnection;
use super::*;

pub(crate) struct PointCorrection {
    /// Subtract from the contact's original string row; torsion stays zero.
    pub string_shapes: Vec<f64>,
    /// Add each body row to the complete map; retain the actual opposite surface.
    pub bodies: Vec<ModalAttachment>,
}

/// One original bridge and any admitted orthogonal directed bridge per string.
/// Keep every support's original receiving body, combining rows only within
/// the same body basis. The existing contact owner assembles their reciprocal
/// global column once, with no additional state or integration path.
pub(crate) fn point_residual(
    strings: &[MaterialStringMember],
    cards: &[BowedStringCard],
    directed: &[StringBodyTranslationConnection],
    points: &[ResidualStop],
    gate: &CancelGate,
) -> Result<(Vec<PointCorrection>, Vec<f64>), ReciprocalStringError> {
    poll(gate)?;
    let bridges = residual_bridges(strings, cards)?;
    let (cross, mut compliance) = stop_residual(&bridges, points)?;
    let mut corrections = Vec::with_capacity(points.len());
    for (point, factor) in points.iter().zip(cross) {
        poll(gate)?;
        let member = &strings[point.string];
        let connection = &member.connections[0];
        let bridge = member
            .source
            .bridge_shape(&cards[point.string], connection.station_fraction)?;
        let body = ModalAttachment {
            component: cards.len() + connection.body.component,
            shapes: connection.body.shapes.iter().map(|b| factor * b).collect(),
        };
        corrections.push(PointCorrection {
            string_shapes: bridge.into_iter().map(|b| factor * b).collect(),
            bodies: if body.shapes.iter().any(|b| *b != 0.0) {
                vec![body]
            } else {
                Vec::new()
            },
        });
    }
    let q = points.len();
    let mut toward = vec![0.0; q];
    for connection in directed {
        poll(gate)?;
        if !points.iter().any(|p| p.string == connection.string) {
            continue;
        }
        let member = strings
            .get(connection.string)
            .ok_or_else(|| invalid("residual directed bridge names an unknown string"))?;
        let card = &cards[connection.string];
        let axis = connection.direction.components();
        let station = connection.station_fraction;
        let condensed = condense(
            dot2(axis, axis) * residual_compliance_m_n(card, station, station)?,
            connection.stiffness_n_m,
            connection.damping_n_s_m,
        )?;
        let bridge = member
            .source
            .transverse_shape_along(card, station, connection.direction)?;
        for (i, point) in points.iter().enumerate() {
            if point.string != connection.string {
                toward[i] = 0.0;
                continue;
            }
            toward[i] = dot2(point.direction, axis)
                * residual_compliance_m_n(card, point.station_fraction, station)?;
            let factor = condensed.stiffness_n_m * toward[i];
            if !factor.is_finite() {
                return Err(invalid(
                    "point/bridge residual correction is not representable",
                ));
            }
            if factor == 0.0 {
                continue;
            }
            let correction = &mut corrections[i];
            if correction.string_shapes.len() != bridge.len() {
                return Err(invalid(
                    "residual bridge maps must use the same material string basis",
                ));
            }
            for (shape, b) in correction.string_shapes.iter_mut().zip(&bridge) {
                *shape += factor * b;
            }
            let component = cards.len() + connection.body.component;
            if let Some(body) = correction
                .bodies
                .iter_mut()
                .find(|body| body.component == component)
            {
                if body.shapes.len() != connection.body.shapes.len() {
                    return Err(invalid(
                        "residual body corrections must use the same receiving body basis",
                    ));
                }
                for (shape, b) in body.shapes.iter_mut().zip(&connection.body.shapes) {
                    *shape += factor * b;
                }
            } else {
                correction.bodies.push(ModalAttachment {
                    component,
                    shapes: connection.body.shapes.iter().map(|b| factor * b).collect(),
                });
            }
        }
        // G_ij - sum_b k_eff,b G_ib G_jb. Compute each pair once so symmetry
        // is exact even after multiple axes and differently signed contacts.
        for i in 0..q {
            if points[i].string != connection.string {
                continue;
            }
            for j in i..q {
                if points[j].string != connection.string {
                    continue;
                }
                let g = compliance[i * q + j] - condensed.stiffness_n_m * (toward[i] * toward[j]);
                if !g.is_finite() {
                    return Err(invalid("point residual compliance is not representable"));
                }
                compliance[i * q + j] = g;
                compliance[j * q + i] = g;
            }
        }
    }
    for correction in &mut corrections {
        correction
            .bodies
            .sort_unstable_by_key(|body| body.component);
    }
    Ok((corrections, compliance))
}
