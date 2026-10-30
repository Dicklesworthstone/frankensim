//! Whole-string observations enter the same reciprocal renderer without an
//! invented bulk material, winding construction, or intrinsic damping law.
use super::*;
use fs_couple::material_assembly::measured_string::{
    MeasuredStringProperties, ResolvedMeasuredStringSpecimen,
};

pub(super) fn resolve(
    wire: &input::Wire,
    axes: &[(String, f64)],
    properties: &MeasuredStringProperties,
    catalog: &mut Catalog,
    loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
    spec: &input::Spec,
    gate: &CancelGate,
) -> Result<MaterialStringSource> {
    poll(gate)?;
    let StringPrestress::FixedTension(tension) = wire.prestress else {
        return Err("measured whole strings require their measured fixed-tension setup; axial stiffness is unknown".into());
    };
    if wire.twist_modes != 0
        || wire.thermal_reference_kelvin.is_some()
        || wire.torsion_memory.is_some()
        || wire.tension_loss.is_some()
        || !matches!(&wire.bending_loss, input::BendingLoss::KelvinVoigt)
        || spec.longitudinal_waves.is_some()
        || spec.geometric_stretching.is_some()
    {
        return Err("measured whole strings supply transverse mass and bending stiffness, not torsion, thermal, axial, or intrinsic-loss laws".into());
    }
    load_pack(&wire.pack, catalog, loader)?;
    let card = &catalog[&wire.pack].card;
    // A source's temperature_known=0 is not a measurement at the body's T.
    // Every measurement coordinate is declared explicitly and the original
    // source owns quantity semantics, admissibility, uncertainty and receipts.
    let query = fs_couple::material_source::MaterialRecord {
        pack: wire.pack.clone(),
        state: card.id().process.clone(),
        axes: axes.to_vec(),
    }
    .query_point(card)?;
    let measurements = resolve_material_state_point(
        card,
        &query,
        &properties.requirements()?,
        Pick::SingleClaimOnly,
    )?;
    let specimen = ResolvedMeasuredStringSpecimen::try_new(
        wire.bending_modes,
        &measurements,
        properties.clone(),
    )?;
    let string = specimen.string();
    for (name, declared, measured) in [
        ("length [m]", wire.length, string.length_m),
        ("radius [m]", wire.radius, specimen.radius_m()),
        ("tension [N]", tension, string.tension_n),
    ] {
        if (declared - measured).abs() > 8.0 * f64::EPSILON * declared.abs().max(measured.abs()) {
            return Err(format!("string {name} disagrees with its measured setup: declared {declared}, measured {measured}").into());
        }
    }
    let mut source = MaterialStringSource::Measured(specimen);
    if let Some(directions) = wire.polarizations {
        source = source.with_transverse_polarizations(directions)?;
    }
    poll(gate)?;
    Ok(source)
}
