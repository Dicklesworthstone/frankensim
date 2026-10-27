//! Project physical nodal translations on the body's original material basis.
//! The same signed row measures displacement and applies the reaction force;
//! a lateral bridge never borrows a normal-force row or a separately solved basis.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::StringBodyTranslationConnection;
use fs_couple::render::schedule::force::coupled::ModalAttachment;

pub(super) fn bind(
    bank: ReciprocalStringBank,
    declarations: &[input::TranslationalBridge],
    gate: &CancelGate,
) -> Result<ReciprocalStringBank> {
    if declarations.is_empty() { return Ok(bank); }
    let mut joints = Vec::with_capacity(declarations.len());
    for joint in declarations {
        poll(gate)?;
        let source = bank.source();
        let modes = if let Some(image) = source.material_shell_assembly_image() {
            image.physical_modes().get(joint.shell_component)
                .ok_or("translational bridge names an absent shell component")?.as_slice()
        } else if let Some(image) = source.material_shell_image() {
            if joint.shell_component != 0 {
                return Err("a single-shell body has only component zero".into());
            }
            image.physical_modes()
        } else {
            return Err("translational bridge geometry requires a shell or assembled-shell --body-profile".into());
        };
        let offset = joint.node.checked_mul(6).ok_or("shell node offset overflow")?;
        let end = offset.checked_add(3).ok_or("shell translation extent overflow")?;
        let norm2 = joint.body_direction.iter().map(|x| x*x).sum::<f64>();
        if !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12 {
            return Err("translational bridge body direction must remain a finite unit vector".into());
        }
        let shapes = modes.iter().map(|mode| -> Result<f64> {
            let translation = mode.get(offset..end)
                .ok_or("translational bridge names an absent physical shell node")?;
            Ok(translation.iter().zip(joint.body_direction).map(|(x, d)| x*d).sum())
        }).collect::<Result<Vec<_>>>()?;
        joints.push(StringBodyTranslationConnection {
            string: joint.string, station_fraction: joint.station_fraction,
            direction: joint.direction,
            // Shell assembly parts share ONE assembled modal body at index zero.
            body: ModalAttachment { component: 0, shapes },
            stiffness_n_m: joint.stiffness_n_m, damping_n_s_m: joint.damping_n_s_m,
            rest_extension_m: joint.rest_extension_m,
        });
    }
    Ok(bank.with_translational_body_connections(joints, gate)?)
}
