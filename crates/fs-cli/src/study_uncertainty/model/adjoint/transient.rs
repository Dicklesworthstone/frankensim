//! Bind native mean controls to the existing accepted-trajectory adjoint.
//! The probability samples still run the unchanged physical project. This
//! adapter selects the requested final or sampled-peak report and converts its
//! physical coordinates. Fixed peak controls may use a deterministic tied branch.

use fs_project::{ConductionTransient, ProjectSpec, ThermalBoundaryCondition, spec::OutputRequest};

use super::{J, NominalCalibration, Result, Target, UniformParameter, array, invalid};

#[path = "transient/peak.rs"]
mod peak;
#[path = "transient/schedule.rs"]
mod schedule;

const FINAL: &str = "temperature-final-adjoint";
const FINAL_MATERIAL: &str = "temperature-final-material-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const PEAK_MATERIAL: &str = "temperature-sampled-peak-material-adjoint";

pub(super) fn policy(spec: &ProjectSpec) -> Option<&ConductionTransient> {
    spec.cooling
        .as_ref()?
        .conduction
        .as_ref()?
        .transient
        .as_ref()
}

fn admit_controls(spec: &ProjectSpec, parameters: &[UniformParameter]) -> Result<()> {
    let setup = spec
        .cooling
        .as_ref()
        .and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| invalid("transient calibration has no conduction setup"))?;
    let time = policy(spec).ok_or_else(|| invalid("transient calibration has no time policy"))?;
    for parameter in parameters.iter().filter(|p| p.low != p.high) {
        let supported = match parameter.target {
            Target::Power => setup.regions.iter().any(|r| r.region == parameter.entity),
            Target::PowerScheduleScale => time.power_schedules.iter()
                .filter(|row| row.region == parameter.entity).count() == 1
                && setup.regions.iter().any(|row| row.region == parameter.entity),
            Target::ConvectionCoefficient | Target::ConvectionTemperature => true,
            Target::AirInletTemperature
            | Target::FanSpeedRatio
            | Target::NaturalConvectionAmbient
            | Target::RadiationReservoirTemperature => true,
            Target::FixedTemperature | Target::HeatFlux => setup.boundaries.iter().any(|boundary| {
                boundary.target == parameter.entity && matches!(
                    (&boundary.condition, parameter.target),
                    (ThermalBoundaryCondition::FixedTemperature { .. }, Target::FixedTemperature)
                        | (ThermalBoundaryCondition::HeatFlux { .. }, Target::HeatFlux)
                )
            }),
            Target::ContactPressure
            | Target::ContactThickness
            | Target::ContactGap
            | Target::ContactTorque => spec.interface_cards.as_deref().unwrap_or(&[])
                .iter().any(|binding| binding.interface == parameter.entity),
        };
        if !supported {
            return Err(invalid(format!(
                "transient nominal adjoint does not supply the declared coordinate {}; use support secants or a supported trajectory control",
                parameter.name,
            )));
        }
    }
    Ok(())
}

/// Preserve report intent and storage physics; do not replace a peak request
/// with a different objective or add a conflicting second adjoint.
pub(super) fn configure(
    spec: &mut ProjectSpec, parameters: &[UniformParameter], sampled_peak: bool,
) -> Result<()> {
    admit_controls(spec, parameters)?;
    let (ordinary, material) = if sampled_peak { (PEAK, PEAK_MATERIAL) }
        else { (FINAL, FINAL_MATERIAL) };
    let needs_contacts = parameters.iter()
        .any(|parameter| parameter.target.contact_axis().is_some() && parameter.low != parameter.high);
    let outputs = spec.outputs.as_deref().unwrap_or(&[]);
    let requests: Vec<_> = outputs
        .iter()
        .filter(|row| {
            matches!(
                row.name.as_str(),
                FINAL
                    | FINAL_MATERIAL
                    | PEAK
                    | PEAK_MATERIAL
                    | super::OUTPUT
                    | super::output::CONTACT
                    | super::output::COMBINED
                    | "temperature-max-boundary-adjoint"
                    | "temperature-volume-mean-adjoint"
            )
        })
        .collect();
    if requests.len() > 1
        || requests.first().is_some_and(|row| {
            (row.name != ordinary && row.name != material)
                || row.kind != "report"
                || row.region.is_some()
        })
    {
        return Err(invalid(
            "transient uncertainty calibration requires one adjoint report for its selected final or sampled-peak objective, without an output region",
        ));
    }
    if requests.is_empty() {
        spec.outputs
            .get_or_insert_with(Vec::new)
            .push(OutputRequest {
                name: if needs_contacts { material } else { ordinary }.into(),
                kind: "report".into(),
                region: None,
            });
    } else if needs_contacts {
        // Keep the selected functional while explicitly admitting the
        // material contractions and their work on the same native trajectory.
        for output in spec.outputs.iter_mut().flatten() {
            if output.name == ordinary { output.name = material.into(); }
        }
    }
    Ok(())
}

fn number(row: &J, key: &str) -> Result<f64> {
    row.f64_field(key)
        .filter(|v| v.is_finite())
        .ok_or_else(|| invalid(format!("transient nominal adjoint has no finite {key}")))
}

fn named<'a>(owner: &'a J, array_name: &str, key: &str, entity: &str) -> Result<&'a J> {
    let mut rows = array(owner, array_name, 64)?
        .iter()
        .filter(|row| row.str_field(key) == Some(entity));
    let row = rows.next().ok_or_else(|| {
        invalid(format!(
            "transient nominal adjoint has no {array_name} control for {entity}"
        ))
    })?;
    if rows.next().is_some() {
        return Err(invalid(format!(
            "ambiguous transient {array_name} control for {entity}"
        )));
    }
    Ok(row)
}

fn family<'a>(report: &'a J, name: &str) -> Result<&'a J> {
    report
        .get(name)
        .ok_or_else(|| invalid(format!("transient nominal adjoint lacks {name} controls")))
}

/// Total-enthalpy controls are absolute original-boundary derivatives. Their
/// physical point and selected prefix must agree with this calibration, rather
/// than a frozen radiative Robin row or an unrelated same-unit temperature.
fn enthalpy_convection(
    report: &J,
    spec: &ProjectSpec,
    entity: &str,
    target: Target,
    fine_steps: usize,
) -> Result<f64> {
    let controls = family(report, "prescribed_convection")?;
    if controls.f64_field("control_steps") != Some(fine_steps as f64)
        || controls.str_field("parameter_scope")
            != Some("independent-original-reference-and-positive-coefficient-held-throughout-selected-trajectory")
    {
        return Err(invalid("prescribed-convection controls do not span the final trajectory"));
    }
    let row = named(controls, "boundaries", "target", entity)?;
    let setup = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
        .ok_or_else(|| invalid("missing sampled prescribed-convection setup"))?;
    let mut laws = setup.boundaries.iter().filter_map(|boundary| {
        match &boundary.condition {
            ThermalBoundaryCondition::Convection { coefficient, reference_temperature }
                if boundary.target == entity => Some((coefficient.value, reference_temperature.value)),
            _ => None,
        }
    });
    let (coefficient, reference) = laws.next()
        .ok_or_else(|| invalid("transient calibration lost its prescribed-convection boundary"))?;
    if laws.next().is_some()
        || !(coefficient.is_finite() && coefficient > 0.0
            && reference.is_finite() && reference > 0.0)
        || number(row, "coefficient_w_m2_k")?.to_bits() != coefficient.to_bits()
        || number(row, "reference_temperature_k")?.to_bits() != reference.to_bits()
    {
        return Err(invalid("prescribed-convection control differs from the sampled original boundary"));
    }
    number(row, match target {
        Target::ConvectionCoefficient => "dtemperature_dcoefficient_k_per_w_m2_k",
        Target::ConvectionTemperature => "dtemperature_dreference_temperature_k_per_k",
        _ => return Err(invalid("requested coordinate is not prescribed convection")),
    })
}

/// Prescribed bath and outward-flux coordinates retain the same meaning in
/// capacity and enthalpy storage. Bind them to their original native law,
/// sampled value and complete reverse prefix before consuming the derivative.
fn prescribed_boundary(
    report: &J,
    spec: &ProjectSpec,
    entity: &str,
    target: Target,
    fine_steps: usize,
) -> Result<f64> {
    let (name, scope, value_key, derivative_key) = match target {
        Target::FixedTemperature => (
            "prescribed_temperature",
            "independent-uniform-bath-temperatures-held-throughout-selected-trajectory",
            "temperature_k", "dtemperature_dprescribed_temperature_k_per_k",
        ),
        Target::HeatFlux => (
            "prescribed_heat_flux",
            "independent-uniform-outward-heat-fluxes-held-throughout-selected-trajectory",
            "flux_w_m2", "dtemperature_doutward_flux_k_per_w_m2",
        ),
        _ => return Err(invalid("requested coordinate is not a prescribed boundary")),
    };
    let controls = family(report, name)?;
    if controls.f64_field("control_steps") != Some(fine_steps as f64)
        || controls.str_field("parameter_scope") != Some(scope)
    {
        return Err(invalid(format!("{name} controls do not span the final trajectory")));
    }
    let row = named(controls, "boundaries", "target", entity)?;
    if target == Target::FixedTemperature && row.str_field("status") != Some("available") {
        return Err(invalid(format!("transient prescribed-temperature control {entity} is unavailable: {}",
            row.str_field("reason").unwrap_or("no independent bath derivative"))));
    }
    let setup = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
        .ok_or_else(|| invalid("missing sampled prescribed-boundary setup"))?;
    let mut values = setup.boundaries.iter().filter_map(|boundary| {
        if boundary.target != entity { return None; }
        match (&boundary.condition, target) {
            (ThermalBoundaryCondition::FixedTemperature { temperature }, Target::FixedTemperature) =>
                Some(temperature.value),
            (ThermalBoundaryCondition::HeatFlux { outward_flux }, Target::HeatFlux) =>
                Some(outward_flux.value),
            _ => None,
        }
    });
    let sampled = values.next()
        .ok_or_else(|| invalid("transient calibration lost its original prescribed boundary"))?;
    if values.next().is_some() || !sampled.is_finite()
        || number(row, value_key)?.to_bits() != sampled.to_bits()
    {
        return Err(invalid(format!("{name} control differs from the sampled original boundary")));
    }
    number(row, derivative_key)
}

/// The caller first verifies the sealed run, its project and its selected QoI.
/// Bind the derivative to the same selected state, storage law and fine grid.
/// A peak branch supplies fixed controls even when the maximum is not unique;
/// its active state, vertex and tie count remain explicit calibration data.
#[cfg(test)]
pub(super) fn coefficients(
    receipt: &J, spec: &ProjectSpec, parameters: &[UniformParameter],
    nominal: f64, region: &str, sampled_peak: bool,
) -> Result<NominalCalibration> {
    coefficients_with_reference(receipt, spec, parameters, nominal, region, sampled_peak, None)
}

/// Retain the immutable original schedule, including its nonzero watts at a
/// zero-amplitude calibration point. Other coordinates keep their existing
/// sampled-law contractions. No additional trajectory or transpose is solved.
pub(super) fn coefficients_with_reference(
    receipt: &J,
    spec: &ProjectSpec,
    parameters: &[UniformParameter],
    nominal: f64,
    region: &str,
    sampled_peak: bool,
    reference: Option<(&ProjectSpec, &[f64])>,
) -> Result<NominalCalibration> {
    admit_controls(spec, parameters)?;
    let time = policy(spec).ok_or_else(|| invalid("missing transient calibration policy"))?;
    let transient = family(receipt, "transient")?;
    let report = family(transient, "adjoint")?;
    let enthalpy = time.enthalpy.is_some();
    let schema = if enthalpy {
        "fs-cli-native-enthalpy-adjoint-v1"
    } else {
        "fs-cli-native-transient-adjoint-v1"
    };
    let time_schema = if enthalpy {
        "fs-cli-transient-enthalpy-v1"
    } else {
        "fs-cli-transient-conduction-v1"
    };
    let fine_steps = 2 * time
        .coarse_step_ends_s()
        .map_err(|e| invalid(format!("{}: {}", e.code, e.what)))?
        .len();
    let goal = family(transient, "temporal_error")?;
    // The enthalpy receipt predates an explicit goal-region field. Its region
    // is bound by the sealed project and QoI above; capacity also records it.
    let region_matches = match goal.str_field("region") {
        Some(recorded) => recorded == region,
        None => enthalpy,
    };
    let selected = sampled_peak.then(|| peak::bind(receipt, nominal, region, fine_steps, enthalpy))
        .transpose().map_err(invalid)?;
    let prefix_steps = selected.as_ref().map_or(fine_steps, |selected| selected.state);
    let selected_time = selected.as_ref().map_or(time.horizon.value, |selected| selected.time_s);
    let output_matches = if sampled_peak {
        matches!(report.str_field("output"), Some(PEAK | PEAK_MATERIAL))
    } else {
        matches!(report.str_field("output"), Some(FINAL | FINAL_MATERIAL))
    };
    if transient.str_field("schema") != Some(time_schema)
        || transient.str_field("status") != Some("completed")
        || transient.str_field("qoi_time") != Some("final")
        || number(transient, "final_time_s")?.to_bits() != time.horizon.value.to_bits()
        || transient.f64_field("fine_steps") != Some(fine_steps as f64)
        || report.str_field("schema") != Some(schema)
        || !output_matches
        || number(report, "value_k")?.to_bits() != nominal.to_bits()
        || number(report, "time_s")?.to_bits() != selected_time.to_bits()
        || report.f64_field("state_index") != Some(prefix_steps as f64)
        || report.f64_field("max_reverse_steps") != Some(fine_steps as f64)
        || number(report, "maximum_dual_relative_residual")? < 0.0
        || goal.str_field("qoi") != Some("temperature-max")
        || (!sampled_peak && number(goal, "fine_final_k")?.to_bits() != nominal.to_bits())
        || !region_matches
    {
        return Err(invalid(
            "transient adjoint storage, selected state, grid, functional or residual differs from the sealed calibration",
        ));
    }
    let gradient = parameters.iter().enumerate().map(|(index, parameter)| {
        if parameter.low == parameter.high { return Ok(0.0); }
        let entity = parameter.entity.as_str();
        let derivative = match parameter.target {
            Target::PowerScheduleScale => {
                let (original, coordinates) = reference
                    .ok_or_else(|| invalid("schedule sensitivity requires its immutable original workload"))?;
                if coordinates.len() != parameters.len() {
                    return Err(invalid("schedule sensitivity calibration coordinate count differs"));
                }
                let amplitude = coordinates[index];
                if !amplitude.is_finite() || amplitude < parameter.low || amplitude > parameter.high {
                    return Err(invalid("schedule sensitivity point is outside its declared amplitude support"));
                }
                schedule::coefficient(report, original, spec, entity, amplitude, prefix_steps)?
            },
            Target::Power => {
                let power = spec.power.as_deref().unwrap_or(&[]).iter()
                    .find(|row| row.region == entity)
                    .ok_or_else(|| invalid("transient calibration lost its declared power"))?;
                // A schedule supplies absolute delivered watts for the whole
                // horizon. It completely overrides this static declaration.
                if time.power_schedules.iter().any(|row| row.region == entity) { 0.0 }
                else {
                    number(named(report, "regional_power", "region", entity)?,
                        "dtemperature_duniform_added_watts_k_per_w")? * power.duty
                }
            }
            Target::ConvectionCoefficient | Target::ConvectionTemperature if enthalpy =>
                enthalpy_convection(report, spec, entity, parameter.target, prefix_steps)?,
            Target::ConvectionCoefficient | Target::ConvectionTemperature => {
                let row = named(report, "robin", "region", entity)?;
                if parameter.target == Target::ConvectionTemperature {
                    number(row, "dtemperature_dreference_k_per_k")?
                } else {
                    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
                        .ok_or_else(|| invalid("missing sampled convection setup"))?;
                    let h = setup.boundaries.iter().find_map(|boundary| {
                        match &boundary.condition {
                            ThermalBoundaryCondition::Convection { coefficient, .. }
                                if boundary.target == entity => Some(coefficient.value),
                            _ => None,
                        }
                    }).filter(|h| h.is_finite() && *h > 0.0)
                        .ok_or_else(|| invalid("transient log-htc control needs the sampled positive coefficient"))?;
                    number(row, "dtemperature_dlog_htc_k")? / h
                }
            }
            Target::NaturalConvectionAmbient => number(named(family(report, "convection")?,
                "natural_reservoirs", "target", entity)?, "dtemperature_dambient_temperature_k_per_k")?,
            Target::AirInletTemperature => number(named(family(report, "convection")?,
                "air_inlets", "branch", entity)?, "dtemperature_dinlet_temperature_k_per_k")?,
            Target::FanSpeedRatio => {
                let row = named(family(report, "convection")?, "fan_banks", "bank", entity)?;
                if row.str_field("status") != Some("available") {
                    return Err(invalid(format!("transient fan control {entity} is unavailable: {}",
                        row.str_field("reason").unwrap_or("no complete physical derivative"))));
                }
                number(row, "dtemperature_dspeed_ratio_k")?
            }
            Target::FixedTemperature | Target::HeatFlux =>
                prescribed_boundary(report, spec, entity, parameter.target, prefix_steps)?,
            Target::RadiationReservoirTemperature => number(named(family(report, "radiation")?,
                "surfaces", "name", entity)?, "dtemperature_dreservoir_temperature_k_per_k")?,
            Target::ContactPressure | Target::ContactThickness | Target::ContactGap
            | Target::ContactTorque => super::contact::transient_coefficient(report, entity)?,
        };
        if !derivative.is_finite() {
            return Err(invalid(format!("nonfinite transient nominal coefficient for {}", parameter.name)));
        }
        // These coordinates act only on future transport. Do not extend this
        // check to initial-state or constitutive-chart coordinates: those must
        // retain their genuine native initial-observation contribution.
        let future_only = matches!(parameter.target,
            Target::Power | Target::PowerScheduleScale | Target::ConvectionCoefficient | Target::ConvectionTemperature
                | Target::AirInletTemperature | Target::FanSpeedRatio
                | Target::NaturalConvectionAmbient | Target::RadiationReservoirTemperature
                | Target::FixedTemperature | Target::HeatFlux | Target::ContactPressure
                | Target::ContactThickness | Target::ContactGap | Target::ContactTorque);
        if sampled_peak && prefix_steps == 0 && future_only && derivative != 0.0 {
            return Err(invalid("an initial peak cannot depend on a future thermal control"));
        }
        Ok(derivative)
    }).collect::<Result<Vec<_>>>()?;
    Ok(NominalCalibration { gradient, peak_branch: selected.map(|selected| selected.json()) })
}

#[cfg(test)]
#[path = "transient/tests.rs"]
mod tests;
