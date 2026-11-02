//! Exact original layer viscosity in the native string mass coordinates.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::torsion::SHEAR_VISCOSITY_REFERENCE_OMEGA;

impl MultiBowStringBank {
    pub(super) fn string_viscous_heat_source(
        &self,
        string: usize,
        layer: usize,
        gate: &CancelGate,
    ) -> Result<Vec<MaterialViscousHeatSource>, ReciprocalStringError> {
        let member = &self.source.bank.strings[string];
        let card = &self.cards[string];
        let bending = match bending_source(&member.source)? {
            MaterialStringSource::Concentric { section, .. } => *section
                .layer_viscous_stiffness_n_m2_s()
                .get(layer)
                .ok_or_else(|| invalid("viscous heat lost the original bending layer"))?,
            MaterialStringSource::Homogeneous(specimen) if layer == 0 => specimen
                .string()
                .kelvin_voigt_bending
                .map_or(0.0, |law| law.viscous_stiffness_n_m2_s),
            MaterialStringSource::ConcentricProny { .. } => 0.0,
            _ => {
                return Err(invalid(
                    "viscous heat requires its original source-resolved bending operator",
                ));
            }
        };
        let component = &self.system.components()[string];
        let offset: usize = self.system.components()[..string]
            .iter()
            .map(|m| m.modes().len())
            .sum();
        let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) {
            2
        } else {
            1
        };
        let mut rates = Vec::new();
        if bending > 0.0 {
            // This card retains the ACTUAL native transverse inertia, including
            // admitted potential-flow added mass. The caloric inventory remains
            // the original solid annulus; ambient resistance is never included.
            for plane in 0..planes {
                for mode in 0..card.mode_count {
                    poll(gate)?;
                    let k = (mode + 1) as f64 * core::f64::consts::PI / card.length_m;
                    rates.push((
                        offset + plane * card.mode_count + mode,
                        bending * k.powi(4) / card.linear_density_kg_m,
                    ));
                }
            }
        }
        if let Some(twist) = member.source.torsion_source() {
            let section = twist.section();
            let shear = *section
                .layer_viscous_stiffness_n_m2_s()
                .get(layer)
                .ok_or_else(|| invalid("viscous heat lost the original torsional layer"))?;
            if shear > 0.0 {
                if section.shear().elastic().layers()[layer]
                    .material
                    .property(SHEAR_VISCOSITY_REFERENCE_OMEGA)
                    .is_some()
                {
                    return Err(invalid(
                        "constant-loss shear proxies require their own original heat operator",
                    ));
                }
                let twist_offset = component
                    .modes()
                    .len()
                    .checked_sub(twist.mode_count())
                    .ok_or_else(|| {
                        invalid("viscous heat lost the original torsional coordinates")
                    })?;
                for mode in 0..twist.mode_count() {
                    poll(gate)?;
                    let k = (mode + 1) as f64 * core::f64::consts::PI / card.length_m;
                    rates.push((
                        offset + twist_offset + mode,
                        shear * k.powi(2) / section.rotary_inertia_kg_m(),
                    ));
                }
            }
        }
        if rates
            .iter()
            .any(|(_, rate)| !rate.is_finite() || *rate <= 0.0)
        {
            return Err(invalid(
                "source layer viscosity is unrepresentable in the original modal mass basis",
            ));
        }
        Ok(if rates.is_empty() {
            Vec::new()
        } else {
            vec![MaterialViscousHeatSource::from_material_region(
                string, layer, rates,
            )]
        })
    }
}
