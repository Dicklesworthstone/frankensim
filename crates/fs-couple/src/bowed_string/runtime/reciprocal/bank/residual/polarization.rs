//! Fixed point contacts through the full two-support Schur complement. Use
//! the SAME passive elastic factors as the retained string/body network, so
//! nonorthogonal support reactions preserve reciprocal maps and compliance.

use super::super::polarization::StringBodyTranslationConnection;
use super::*;
use fs_math::dd::Dd;

pub(crate) struct PointCorrection {
    /// Subtract from the contact's original string row; torsion stays zero.
    pub string_shapes: Vec<f64>,
    /// Add each body row to the complete map; retain the actual opposite surface.
    pub bodies: Vec<ModalAttachment>,
    /// Fixed point displacement from bridge rest extensions: g_i^T K_eff r.
    /// Positive along this contact's declared direction; it is not actuator work.
    pub displacement_offset_m: f64,
}

/// One original bridge and one independently directed bridge per string.
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
    let pairs = coupled_bridge::pairs(strings, cards, directed, gate)?;
    let bridges = residual_bridges(strings, cards)?;
    let (cross, mut compliance) = stop_residual(&bridges, points)?;
    let mut corrections = Vec::with_capacity(points.len());
    for (point, &factor) in points.iter().zip(&cross) {
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
        let displacement_offset_m = factor * connection.rest_extension_m;
        if !displacement_offset_m.is_finite() {
            return Err(invalid(
                "residual bridge rest displacement is not representable",
            ));
        }
        corrections.push(PointCorrection {
            string_shapes: bridge.into_iter().map(|b| factor * b).collect(),
            bodies: if body.shapes.iter().any(|b| *b != 0.0) {
                vec![body]
            } else {
                Vec::new()
            },
            displacement_offset_m,
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
        let pair = pairs[connection.string]
            .as_ref()
            .ok_or_else(|| invalid("residual directed bridge lost its condensed pair"))?
            .law;
        let original = &member.connections[0];
        // The same e1-h e0 rest coordinate used by the native elastic factor.
        let rest = if pair.mix == 0.0 {
            connection.rest_extension_m
        } else {
            connection.rest_extension_m - pair.mix * original.rest_extension_m
        };
        if !rest.is_finite() {
            return Err(invalid(
                "residual bridge rest displacement is not representable",
            ));
        }
        let first_axis = directions(&member.source).1;
        let first_bridge = member
            .source
            .bridge_shape(card, original.station_fraction)?;
        let condensed = pair.second;
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
            if pair.mix != 0.0 {
                // The second elastic coordinate is e1-h e0, not e1 alone.
                toward[i] -= pair.mix
                    * dot2(point.direction, first_axis)
                    * residual_compliance_m_n(
                        card,
                        point.station_fraction,
                        original.station_fraction,
                    )?;
            }
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
            // Sum before rounding so opposing support-rest contributions do
            // not erase a representable weak point displacement. Constant rest
            // coordinates are already included in the network's spring energy.
            correction.displacement_offset_m = (Dd::from_f64(cross[i])
                * Dd::from_f64(original.rest_extension_m)
                + Dd::from_f64(factor) * Dd::from_f64(rest))
            .to_f64();
            if !correction.displacement_offset_m.is_finite() {
                return Err(invalid(
                    "residual bridge rest displacement is not representable",
                ));
            }
            if correction.string_shapes.len() != bridge.len() {
                return Err(invalid(
                    "residual bridge maps must use the same material string basis",
                ));
            }
            for ((shape, b), a) in correction
                .string_shapes
                .iter_mut()
                .zip(&bridge)
                .zip(&first_bridge)
            {
                let projected = if pair.mix == 0.0 {
                    *b
                } else {
                    b - pair.mix * a
                };
                *shape += factor * projected;
                if !shape.is_finite() {
                    return Err(invalid(
                        "point/bridge residual correction is not representable",
                    ));
                }
            }
            add_body(
                &mut correction.bodies,
                cards.len(),
                &connection.body,
                factor,
            )?;
            if pair.mix != 0.0 {
                add_body(
                    &mut correction.bodies,
                    cards.len(),
                    &original.body,
                    -pair.mix * factor,
                )?;
            }
        }
        // Subtract the second elastic factor from the already-condensed first.
        // This is G_ij-g_i^T K_eff g_j, with exact matrix symmetry.
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

fn add_body(
    bodies: &mut Vec<ModalAttachment>,
    strings: usize,
    original: &ModalAttachment,
    factor: f64,
) -> Result<(), ReciprocalStringError> {
    if factor == 0.0 {
        return Ok(());
    }
    let component = strings + original.component;
    if let Some(body) = bodies.iter_mut().find(|body| body.component == component) {
        if body.shapes.len() != original.shapes.len() {
            return Err(invalid(
                "residual body corrections must use the same receiving body basis",
            ));
        }
        for (shape, b) in body.shapes.iter_mut().zip(&original.shapes) {
            *shape += factor * b;
        }
        if body.shapes.iter().any(|shape| !shape.is_finite()) {
            return Err(invalid(
                "residual receiving-body correction is not representable",
            ));
        }
    } else {
        let shapes: Vec<_> = original.shapes.iter().map(|b| factor * b).collect();
        if shapes.iter().any(|shape| !shape.is_finite()) {
            return Err(invalid(
                "residual receiving-body correction is not representable",
            ));
        }
        bodies.push(ModalAttachment { component, shapes });
    }
    Ok(())
}
