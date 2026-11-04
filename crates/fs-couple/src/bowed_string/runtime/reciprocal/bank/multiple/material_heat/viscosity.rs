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
        let (bending, moment, property) = match bending_source(&member.source)? {
            MaterialStringSource::Concentric { section, .. } => {
                let bending = *section.layer_viscous_stiffness_n_m2_s().get(layer)
                    .ok_or_else(|| invalid("viscous heat lost the original bending layer"))?;
                let property = section.elastic().layers()[layer].material.property(section.viscosity_property())
                    .ok_or_else(|| invalid("viscous heat lost its original bending property selector"))?;
                (bending, section.elastic().layer_integrals()[layer].second_moment_m4(), Some(property))
            }
            MaterialStringSource::Homogeneous(specimen) if layer == 0 => {
                let law = specimen.string().kelvin_voigt_bending;
                let property = law.map(|_| specimen.material().property(
                    crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
                    .ok_or_else(|| invalid("viscous heat lost its homogeneous source selector"))).transpose()?;
                (law.map_or(0.0, |law| law.viscous_stiffness_n_m2_s), specimen.second_moment_m4(), property)
            }
            MaterialStringSource::ConcentricProny { .. } => (0.0, 0.0, None),
            _ => {
                return Err(invalid(
                    "viscous heat requires its original source-resolved bending operator",
                ));
            }
        };
        let varies = |property: &fs_material::state_point::ResolvedScalarProperty|
            property.source_axes().iter().any(|(_, dims)| *dims == fs_qty::Temperature::DIMS);
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
        let mut properties = Vec::new();
        if let Some(property) = property.filter(|p| bending > 0.0 || varies(p)) {
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
                    properties.push((property.requirement().name().to_owned(),
                        moment * k.powi(4) / card.linear_density_kg_m, property.value_si()));
                }
            }
        }
        if let Some(twist) = member.source.torsion_source() {
            let section = twist.section()
                .ok_or_else(|| invalid("viscous heat requires material torsion; measured whole-string torsion has no layers"))?;
            let shear = *section
                .layer_viscous_stiffness_n_m2_s()
                .get(layer)
                .ok_or_else(|| invalid("viscous heat lost the original torsional layer"))?;
            let property = section.shear().elastic().layers()[layer].material.property(section.viscosity_property())
                .ok_or_else(|| invalid("viscous heat lost its original shear property selector"))?;
            if shear > 0.0 || varies(property) {
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
                    properties.push((property.requirement().name().to_owned(),
                        2.0 * section.shear().elastic().layer_integrals()[layer].second_moment_m4()
                            * k.powi(2) / section.rotary_inertia_kg_m(), property.value_si()));
                }
            }
        }
        if rates
            .iter()
            .any(|(_, rate)| !rate.is_finite() || *rate < 0.0)
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
            ).with_property_rows(properties).map_err(ModalFrictionError::from)
                .map_err(ReciprocalStringError::Friction)?]
        })
    }
}
