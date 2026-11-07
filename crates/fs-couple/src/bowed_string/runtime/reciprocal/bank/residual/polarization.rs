//! Fixed point contacts through both independently condensed bridge axes.
//! The bank admits only orthogonal bridge axes, so their omitted-mode Green
//! matrix is diagonal. Sum each axis's Schur term without changing the retained
//! modes, contact force law or time integrator.

use super::super::polarization::StringBodyTranslationConnection;
use super::*;

pub(crate) struct PointCorrection {
    /// Subtract from the contact's original string row; torsion stays zero.
    pub string_shapes: Vec<f64>,
    /// Add to its complete contact map, retaining the actual opposite surface.
    pub body: ModalAttachment,
}

/// One original bridge and any admitted orthogonal directed bridge per string.
/// The contact owner accepts one extra body map: combine all nonzero support
/// corrections in that SAME body basis, or refuse distinct receiving bodies.
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
        corrections.push(PointCorrection {
            string_shapes: bridge.into_iter().map(|b| factor * b).collect(),
            body: ModalAttachment {
                component: cards.len() + connection.body.component,
                shapes: connection.body.shapes.iter().map(|b| factor * b).collect(),
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
            if component != correction.body.component {
                if correction.body.shapes.iter().any(|b| *b != 0.0) {
                    return Err(invalid(
                        "point residual corrections require one receiving body per contact",
                    ));
                }
                correction.body = ModalAttachment {
                    component,
                    shapes: vec![0.0; connection.body.shapes.len()],
                };
            }
            if correction.body.shapes.len() != connection.body.shapes.len() {
                return Err(invalid(
                    "residual body corrections must use the same receiving body basis",
                ));
            }
            for (shape, b) in correction
                .body
                .shapes
                .iter_mut()
                .zip(&connection.body.shapes)
            {
                *shape += factor * b;
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
    Ok((corrections, compliance))
}
