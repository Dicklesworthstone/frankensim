//! Version-5 friction declarations; normal records keep their existing parser.
use super::{Reader, ModalPerformanceError, ModalForceVoice, input, read_attachment};
use crate::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use crate::render::schedule::force::coupled::contact::multiple::friction::ModalFriction;

pub(super) fn read(reader: &mut Reader<'_>, voices: &[ModalForceVoice],
    contacts: &[(ModalContact, ModalContactConfig)], weights: &mut usize)
    -> Result<Vec<Option<ModalFriction>>, ModalPerformanceError>
{
    let count: usize = reader.one("frictions")?;
    if count != contacts.len() {
        return Err(input(reader.line, "frictions count must equal the complete normal-contact count"));
    }
    let mut result = Vec::with_capacity(count);
    for (contact, _) in contacts {
        let mut row = reader.row("friction")?;
        match row.word()? {
            "none" => { row.finish()?; result.push(None); }
            "regularized-coulomb" => {
                let coefficient = row.scalar()?;
                let regularization_speed_m_s = row.scalar()?;
                let maximum_force_n = row.scalar()?;
                let source = row.word()?;
                if source.len() > 1024 { return Err(input(row.line, "friction source label exceeds 1024 bytes")); }
                let source = source.to_string();
                row.finish()?;
                // Reuse bounded, finite attachment decoding. Friction is a
                // second direction of THIS contact pair, never another body pair.
                let left = read_attachment(reader, "friction_left", voices, weights)?;
                if left.component != contact.left.component {
                    return Err(input(reader.line, "friction left body must match its normal contact"));
                }
                let right = read_attachment(reader, "friction_right", voices, weights)?;
                if right.component != contact.right.component {
                    return Err(input(reader.line, "friction right body must match its normal contact"));
                }
                result.push(Some(ModalFriction { left_shapes: left.shapes, right_shapes: right.shapes,
                    coefficient, regularization_speed_m_s, maximum_force_n, source }));
            }
            _ => return Err(input(row.line, "expected explicit none or regularized-coulomb friction model")),
        }
    }
    Ok(result)
}
