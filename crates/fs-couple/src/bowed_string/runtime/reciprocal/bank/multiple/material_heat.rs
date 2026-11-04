//! Original string-layer bending/torsion loss in its actual finite material volume.
//! The source producer owns arm/plane/mode order; no damping fraction, modal
//! mass or invented heat partition replaces the physical annular inventory.
use super::*;
use crate::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
use crate::render::schedule::force::coupled::coefficients::heat::viscosity::MaterialViscousHeatSource;
use crate::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    material::MaterialRelaxationHeatDestination, rates::MaterialRelaxationTemperature,
};
use core::ops::Range;
use fs_matdb::MaterialCard;
use fs_material::{
    phase::EquilibriumEnthalpyPhaseCurve,
    section::{ConcentricElasticSection, relaxation::LayerRelaxationProperties},
    state_point::ResolvedMaterialStatePoint,
    visco::GeneralizedMaxwell,
};
use fs_thermal::enthalpy::CellId;

mod viscosity;

struct Source {
    material: ResolvedMaterialStatePoint,
    volume_m3: f64,
    branches: Vec<Range<usize>>,
    times: Vec<(usize, String)>,
    viscosity: Vec<MaterialViscousHeatSource>,
}

// Unwrap the same physical section. Torsional loss is collected separately in
// its original native block; axial/nonlinear storage needs its own continuation.
fn bending_source(
    source: &MaterialStringSource,
) -> Result<&MaterialStringSource, ReciprocalStringError> {
    match source {
        MaterialStringSource::WithPolarizations(p) => bending_source(p.source()),
        MaterialStringSource::WithTorsion(t) => bending_source(t.bending()),
        MaterialStringSource::Homogeneous(_)
        | MaterialStringSource::Concentric { .. }
        | MaterialStringSource::ConcentricProny { .. } => Ok(source),
        _ => Err(invalid(
            "string material heat requires original homogeneous or bonded layers",
        )),
    }
}

#[allow(clippy::too_many_arguments)]
fn layer_source(
    section: &ConcentricElasticSection,
    spectra: &[GeneralizedMaxwell],
    properties: &[LayerRelaxationProperties],
    layer: usize,
    length_m: f64,
    per_arm: usize,
    mut start: usize,
    gate: &CancelGate,
) -> Result<Source, ReciprocalStringError> {
    if layer >= section.layers().len()
        || spectra.len() != section.layers().len()
        || properties.len() != spectra.len()
        || per_arm == 0
    {
        return Err(invalid(
            "string material heat requires one original physical layer and complete spectra",
        ));
    }
    for spectrum in &spectra[..layer] {
        poll(gate)?;
        start += per_arm
            * spectrum
                .terms
                .iter()
                .filter(|(increment, _)| *increment > 0.0)
                .count();
    }
    let spectrum = &spectra[layer];
    let roles = &properties[layer];
    if spectrum.terms.len() != roles.branches.len() {
        return Err(invalid(
            "string heat lost the original source arm selectors",
        ));
    }
    let mut times = Vec::new();
    for (arm, &(increment, _)) in spectrum.terms.iter().enumerate() {
        poll(gate)?;
        if increment == 0.0 {
            continue;
        }
        for _ in 0..per_arm {
            times.push((start + times.len(), roles.branches[arm].1.clone()));
        }
    }
    let volume_m3 = section.layer_integrals()[layer].area_m2() * length_m;
    if !volume_m3.is_finite() || volume_m3 <= 0.0 {
        return Err(invalid(
            "string heat requires an actual finite layer volume",
        ));
    }
    Ok(Source {
        material: section.layers()[layer].material.clone(),
        volume_m3,
        branches: if times.is_empty() {
            Vec::new()
        } else {
            vec![start..start + times.len()]
        },
        times,
        viscosity: Vec::new(),
    })
}

fn elastic_source(
    source: &MaterialStringSource,
    layer: usize,
    length_m: f64,
) -> Result<Source, ReciprocalStringError> {
    let (material, area) = match source {
        MaterialStringSource::Homogeneous(specimen) if layer == 0 => {
            (specimen.material().clone(), specimen.area_m2())
        }
        MaterialStringSource::Concentric { section, .. } => {
            let original = section
                .elastic()
                .layers()
                .get(layer)
                .ok_or_else(|| invalid("string heat names no original material layer"))?;
            (
                original.material.clone(),
                section.elastic().layer_integrals()[layer].area_m2(),
            )
        }
        _ => {
            return Err(invalid(
                "string heat needs its original homogeneous or bonded loss section",
            ));
        }
    };
    let volume_m3 = area * length_m;
    if !volume_m3.is_finite() || volume_m3 <= 0.0 {
        return Err(invalid("string heat requires its physical layer volume"));
    }
    Ok(Source {
        material,
        volume_m3,
        branches: Vec::new(),
        times: Vec::new(),
        viscosity: Vec::new(),
    })
}

impl MultiBowStringBank {
    fn string_heat_source(
        &self,
        string: usize,
        layer: usize,
        gate: &CancelGate,
    ) -> Result<Source, ReciprocalStringError> {
        let mut source = self.string_bending_heat_source(string, layer, gate)?;
        let bank = &self.source.bank;
        if let Some(twist) = bank.strings[string].source.torsion_source() {
            // The torsion constructor retains this same section. Check the
            // full original bundle and volume rather than infer a heat split
            // from GJ/EI or from the relative vibration amplitudes.
            let section = twist.section()
                .ok_or_else(|| invalid("string material heat requires material torsion; measured whole-string torsion has no layers"))?
                .shear().elastic();
            if section
                .layers()
                .get(layer)
                .is_none_or(|original| original.material != source.material)
                || (section.layer_integrals()[layer].area_m2() * self.cards[string].length_m)
                    .to_bits()
                    != source.volume_m3.to_bits()
            {
                return Err(invalid(
                    "string bending and torsion heat must retain the same original layer and volume",
                ));
            }
            if let Some(law) = twist.relaxation() {
                // ALL intrinsic bending rows precede torsion; explicit bending
                // rows are appended after source/body histories. Preserve both
                // orders and leave intervening strings' material states alone.
                let mut start = 0;
                for (member, card) in bank.strings.iter().zip(&self.cards) {
                    poll(gate)?;
                    if let Some(bending) = member.source.bending_relaxation() {
                        let planes =
                            if matches!(&member.source, MaterialStringSource::WithPolarizations(_))
                            {
                                2
                            } else {
                                1
                            };
                        start += planes
                            * card.mode_count
                            * bending
                                .branches
                                .iter()
                                .filter(|arm| arm.relaxing_stiffness_n_m2 > 0.0)
                                .count();
                    }
                }
                for member in &bank.strings[..string] {
                    poll(gate)?;
                    if let Some(previous) = member.source.torsion_source() {
                        if let Some(law) = previous.relaxation() {
                            start += previous.mode_count()
                                * law
                                    .branches()
                                    .iter()
                                    .filter(|arm| arm.relaxing_stiffness_n_m2 > 0.0)
                                    .count();
                        }
                    }
                }
                let torsion = layer_source(
                    section,
                    law.spectra(),
                    law.properties(),
                    layer,
                    self.cards[string].length_m,
                    twist.mode_count(),
                    start,
                    gate,
                )?;
                source.branches.extend(torsion.branches);
                source.times.extend(torsion.times);
            }
        }
        source.viscosity = self.string_viscous_heat_source(string, layer, gate)?;
        if source.branches.is_empty() && source.viscosity.is_empty() {
            return Err(invalid(
                "string heat requires positive original material loss channels",
            ));
        }
        source.branches.sort_unstable_by_key(|range| range.start);
        source.times.sort_unstable_by_key(|(branch, _)| *branch);
        Ok(source)
    }

    fn string_bending_heat_source(
        &self,
        string: usize,
        layer: usize,
        gate: &CancelGate,
    ) -> Result<Source, ReciprocalStringError> {
        poll(gate)?;
        let bank = &self.source.bank;
        let member = bank
            .strings
            .get(string)
            .ok_or_else(|| invalid("string heat names no original material string"))?;
        let source = bending_source(&member.source)?;
        if bank.geometric_stretching.is_some() || bank.longitudinal_waves.is_some() {
            return Err(invalid(
                "string body heat needs a complete axial/nonlinear heat continuation before including stretching",
            ));
        }
        let string_model = source.resolved_transverse_string()?;
        if string_model.damping_ratio != 0.0
            || string_model.rayleigh.is_some()
            || string_model
                .tension_loss_factor
                .is_some_and(|loss| loss != 0.0)
            || string_model
                .bending_loss_factor
                .is_some_and(|loss| loss != 0.0)
        {
            return Err(invalid(
                "string body heat cannot omit parallel intrinsic material loss",
            ));
        }
        let card = &self.cards[string];
        let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) {
            2
        } else {
            1
        };
        let per_arm = planes * card.mode_count;
        if let MaterialStringSource::ConcentricProny { section, .. } = source {
            // All original bending histories precede torsion/body/additional
            // histories, exactly as source_memory_branches constructs them.
            let mut start = 0;
            for (previous, card) in bank.strings[..string].iter().zip(&self.cards[..string]) {
                poll(gate)?;
                if let Some(law) = previous.source.bending_relaxation() {
                    let planes =
                        if matches!(&previous.source, MaterialStringSource::WithPolarizations(_)) {
                            2
                        } else {
                            1
                        };
                    start += planes
                        * card.mode_count
                        * law
                            .branches
                            .iter()
                            .filter(|arm| arm.relaxing_stiffness_n_m2 > 0.0)
                            .count();
                }
            }
            return layer_source(
                section.elastic(),
                section.spectra(),
                section.properties(),
                layer,
                card.length_m,
                per_arm,
                start,
                gate,
            );
        }
        if source.bending_relaxation().is_some() {
            return Err(invalid(
                "homogeneous intrinsic spectra need their complete original layer-role binding",
            ));
        }
        let Some(image) = bank.string_bending_relaxation.as_ref() else {
            return elastic_source(source, layer, card.length_m);
        };
        let Some(selected) = image
            .layered_selections()
            .iter()
            .position(|selection| selection.string == string)
        else {
            return elastic_source(source, layer, card.length_m);
        };
        let prefix = image
            .channels()
            .iter()
            .position(|channel| channel.string == string)
            .ok_or_else(|| invalid("selected string layer has no physical bending dissipation"))?;
        let body_count = bank
            .material_chart_relaxation
            .as_ref()
            .map_or(0, |image| image.coupling_rows_per_s().len())
            + bank
                .material_chart_tensor_relaxation
                .as_ref()
                .map_or(0, |image| image.coupling_rows_per_s().len())
            + bank
                .material_shell_coating_relaxation
                .as_ref()
                .map_or(0, |image| image.coupling_rows_per_s().len());
        let start = bank.source_memory_count(&self.cards, gate)? + body_count + prefix;
        layer_source(
            &image.sections()[selected],
            &image.spectra()[selected],
            &image.layered_selections()[selected].layers,
            layer,
            card.length_m,
            per_arm,
            start,
            gate,
        )
    }

    /// Bind every original bending/torsion-memory arm/plane/mode of one layer
    /// to that layer's finite uniform-temperature caloric body. The annular
    /// area and full supported span supply its volume; the existing heat owner
    /// integrates actual resistance loss along the original mechanical path.
    ///
    /// Original homogeneous/bonded Kelvin-Voigt sources, bonded Prony sections,
    /// and explicit layered bending spectra are supported, including torsional
    /// Prony memory and native shear viscosity. Their original source operators
    /// and disjoint bending/torsion blocks heat the SAME annular
    /// cell once. Axial histories, geometric stretching, measured/wound sections
    /// and modal/tension loss or constant-loss shear proxies refuse
    /// rather than silently losing part of one material body's heat. Gas drag,
    /// body/contact loss and radiation retain their independent owners. Stored
    /// coefficients and reference geometry remain fixed; all source and caloric
    /// domain checks belong to the original thermal transaction.
    #[allow(clippy::too_many_arguments)]
    pub fn string_material_heat_destination<'a>(
        &self,
        string: usize,
        layer: usize,
        cell: CellId,
        material: &'a MaterialCard,
        caloric: &'a EquilibriumEnthalpyPhaseCurve,
        temperature_axis: &'a str,
        budget: MaterialHeatBudget,
        gate: &CancelGate,
    ) -> Result<MaterialRelaxationHeatDestination<'a>, ReciprocalStringError> {
        let source = self.string_heat_source(string, layer, gate)?;
        let channel = self
            .system
            .prepare_material_heat_with_viscosity(&source.branches, &source.viscosity, budget, gate)
            .map_err(ReciprocalStringError::Friction)?;
        Ok(MaterialRelaxationHeatDestination {
            cell,
            material,
            caloric,
            temperature_axis,
            density_property: "density",
            volume_m3: source.volume_m3,
            source: source.material,
            channel: Some(channel),
            nonlinear_channels: Vec::new(),
        })
    }

    /// The same complete material layer with its original source mobility roles.
    /// The existing accepted-rate/endpoint-temperature transaction changes only
    /// relaxation times and original bending/shear viscosities, retaining material strain/storage and the full coupled
    /// state. This is the same first-order temperature partition as fiber heat;
    /// it does not re-resolve geometry, modulus, density or prestress. Purely
    /// Kelvin-Voigt layers need no additional material-memory coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn string_material_rate_heat_destination<'a>(
        &self,
        string: usize,
        layer: usize,
        cell: CellId,
        material: &'a MaterialCard,
        caloric: &'a EquilibriumEnthalpyPhaseCurve,
        temperature_axis: &'a str,
        budget: MaterialHeatBudget,
        gate: &CancelGate,
    ) -> Result<
        (
            MaterialRelaxationHeatDestination<'a>,
            MaterialRelaxationTemperature,
        ),
        ReciprocalStringError,
    > {
        let source = self.string_heat_source(string, layer, gate)?;
        if source.times.is_empty() && source.viscosity.iter().all(|v| v.property_names().next().is_none()) {
            return Err(invalid(
                "source-temperature mobility feedback needs original Maxwell or viscosity selectors",
            ));
        }
        let channel = self
            .system
            .prepare_material_heat_with_viscosity(&source.branches, &source.viscosity, budget, gate)
            .map_err(ReciprocalStringError::Friction)?;
        Ok((
            MaterialRelaxationHeatDestination {
                cell,
                material,
                caloric,
                temperature_axis,
                density_property: "density",
                volume_m3: source.volume_m3,
                source: source.material,
                channel: Some(channel),
                nonlinear_channels: Vec::new(),
            },
            MaterialRelaxationTemperature {
                cell,
                times: source.times,
            },
        ))
    }
}
