//! Keep the existing material loss, but move its GAS contribution to the
//! absolute-velocity image. Do not replace either source with a damping preset.
use super::{NonlinearMaterialString, ImpactError, CancelGate, invalid, poll};
use crate::bowed_string::runtime::reciprocal::bank::MaterialStringSource;
use crate::air_path::oscillating_cylinder_air_resistance_per_length;
use crate::render::plate::impact::linear::wire::WireSpan;
use crate::render::plate::impact::supported::AbsoluteStringDrag;

pub(super) fn split_loss(source: &NonlinearMaterialString, remainder_hz: f64,
    gate: &CancelGate) -> Result<(WireSpan, AbsoluteStringDrag), ImpactError>
{
    poll(gate)?;
    let original=match source.source() {
        MaterialStringSource::WithPolarizations(p) => p.source(),
        other => other,
    };
    let bending=match original {
        MaterialStringSource::WithTorsion(t) => t.bending(),
        other => other,
    };
    let radius=match bending {
        MaterialStringSource::Homogeneous(s) => 0.5*s.string().width_m,
        MaterialStringSource::Measured(s) => s.radius_m(),
        MaterialStringSource::Wound(s) => s.section().outer_radius_m(),
        MaterialStringSource::Concentric { section, .. } => section.elastic().outer_radius_m(),
        MaterialStringSource::ConcentricProny { section, .. } => section.elastic().outer_radius_m(),
        MaterialStringSource::ConcentricModalLoss { section, .. } => section.elastic().outer_radius_m(),
        MaterialStringSource::WithTorsion(_) => return Err(invalid("nested torsion is not a physical gas-drag source")),
        MaterialStringSource::WithPolarizations(_) => return Err(invalid("nested polarizations are not a physical gas-drag source")),
    };
    let card=source.card();
    let mut wire=source.wire_span().clone();
    let mut modal_per_s=Vec::with_capacity(card.mode_count);
    for (i,total) in wire.damping_per_s.iter_mut().enumerate() {
        poll(gate)?;
        let omega=card.mode_omega_rad_s(i+1);
        let resistance=oscillating_cylinder_air_resistance_per_length(radius,omega,source.gas())
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        // Reproduce the original gas contribution's arithmetic, not a second
        // approximation of intrinsic Kelvin-Voigt loss. With zero material
        // loss the subtraction is exact. Unresolved positive source loss
        // refuses: no epsilon floor or sign clipping is used.
        let zeta=0.5*resistance/card.linear_density_kg_m/omega;
        let gas_rate=2.0*zeta*omega;
        let material_rate=*total-gas_rate;
        if !gas_rate.is_finite() || gas_rate<=0.0 || !material_rate.is_finite()
            || material_rate<0.0 || (card.viscous_bending_n_m2_s>0.0 && material_rate==0.0) {
            return Err(invalid("material and gas modal losses cannot be separated at the retained precision"));
        }
        *total=material_rate;modal_per_s.push(gas_rate);
    }
    if matches!(source.source(),MaterialStringSource::WithPolarizations(_)) {
        // The same gas resistance acts on two velocity components. The
        // orthogonal support-lift remainder is still one physical field.
        modal_per_s.extend_from_within(..);
    }
    let resistance=oscillating_cylinder_air_resistance_per_length(radius,
        core::f64::consts::TAU*remainder_hz,source.gas())
        .map_err(|e| ImpactError::Owner(e.to_string()))?;
    let remainder_per_s=resistance/card.linear_density_kg_m;
    if !remainder_per_s.is_finite() || remainder_per_s<=0.0 {
        return Err(invalid("stationary gas remainder resistance is unrepresentable"));
    }
    Ok((wire,AbsoluteStringDrag {modal_per_s,remainder_per_s}))
}
