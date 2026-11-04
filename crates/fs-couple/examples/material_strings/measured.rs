//! Whole-string observations enter the same reciprocal renderer without an
//! invented bulk material, winding construction, or intrinsic damping law. A
//! damping law exists only as the source's own rows, selected by name.
//! Torsion exists only as the string's own measured torsion, selected by a
//! `torsion-measured` row bound to this specimen; without it twist refuses.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::MeasuredStringTorsion;
use fs_couple::material_assembly::measured_string::ResolvedMeasuredStringSpecimen;

pub(super) fn resolve(
    wire: &input::Wire,
    catalog: &mut Catalog,
    loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
    spec: &input::Spec,
    gate: &CancelGate,
) -> Result<MaterialStringSource> {
    poll(gate)?;
    let input::WireSection::Measured { axes, properties, loss, length } = &wire.section else {
        return Err("measured resolution requires a measured section".into());
    };
    let StringPrestress::FixedTension(tension) = wire.prestress else {
        return Err("measured whole strings require their measured fixed-tension setup; axial stiffness is unknown".into());
    };
    if (wire.twist_modes != 0 && wire.measured_torsion.is_none())
        || wire.thermal_reference_kelvin.is_some()
        || wire.torsion_memory.is_some()
        || wire.tension_loss.is_some()
        || !matches!(&wire.bending_loss, input::BendingLoss::KelvinVoigt)
        || spec.longitudinal_waves.is_some()
        || spec.geometric_stretching.is_some()
    {
        return Err("measured whole strings supply transverse mass and bending stiffness, not torsion, thermal, axial, or intrinsic-loss laws from outside their own source rows".into());
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
    let mut requirements = properties.requirements()?;
    if let Some(loss) = loss {
        requirements.extend(loss.requirements()?);
    }
    let measurements = resolve_material_state_point(card, &query, &requirements, Pick::SingleClaimOnly)?;
    let specimen = ResolvedMeasuredStringSpecimen::try_new(
        wire.bending_modes,
        &measurements,
        properties.clone(),
    )?;
    let specimen = match loss {
        None => specimen,
        Some(loss) => specimen.with_loss_factors(loss.clone())?,
    };
    // A speaking-length source measured nut to saddle: bind it there, on this string's own span and bridge station.
    let measured_length = specimen.string().length_m;
    let specimen = match length {
        input::MeasuredLength::TestSpan => specimen,
        input::MeasuredLength::Speaking => specimen
            .on_span_from_speaking_length(wire.length, wire.bridge.station_fraction)
            .map_err(|error| format!("string speaking length [m] disagrees with its measured setup: declared span \
                {} x bridge station {} = {}, measured {measured_length}: {error}", wire.length,
                wire.bridge.station_fraction, wire.length * wire.bridge.station_fraction))?,
    };
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
    let speaking_length = wire.length * wire.bridge.station_fraction;
    let mut source = MaterialStringSource::Measured(specimen);
    if let Some(directions) = wire.polarizations {
        source = source.with_transverse_polarizations(directions)?;
    }
    if let Some(row) = &wire.measured_torsion {
        // The row's numbers belong to one specimen at one speaking length; no other string borrows them.
        if wire.twist_modes == 0 {
            return Err(format!("torsion-measured {} needs the string row's twist mode count, at least 1", row.id).into());
        }
        let pack = catalog[&wire.pack].pack.pack_id();
        if pack != row.pack || properties.linear_density != row.linear_density {
            return Err(format!("torsion-measured {} was measured on {} ({}), not on this string's {} ({})",
                row.id, row.pack, row.linear_density, pack, properties.linear_density).into());
        }
        if (speaking_length - row.speaking_length_m).abs()
            > 8.0 * f64::EPSILON * speaking_length.abs().max(row.speaking_length_m) {
            return Err(format!("torsion-measured {} belongs to a {} m speaking length; span {} x bridge station {} = {speaking_length}",
                row.id, row.speaking_length_m, wire.length, wire.bridge.station_fraction).into());
        }
        let torsion = MeasuredStringTorsion::try_new(row.fundamental_hz, row.impedance_kg_s, row.quality_factor,
            wire.bridge.station_fraction)?;
        source = source.with_measured_torsion(torsion, wire.twist_modes)?;
    }
    poll(gate)?;
    Ok(source)
}

/// The bank's own measured twist of one string, and its twist mode count.
fn measured_twist(source: &MaterialStringSource) -> Option<(&MeasuredStringTorsion, usize)> {
    match source {
        MaterialStringSource::WithPolarizations(p) => measured_twist(p.source()),
        MaterialStringSource::WithTorsion(t) => t.measured().map(|m| (m, t.mode_count())),
        _ => None,
    }
}

/// One receipt line per string with measured torsion, read from the source the bank holds.
pub(super) fn describe_torsion(ready: &Ready) {
    let bank = ready.performance.bank();
    for (string, (wire, member)) in ready.spec.strings.iter().zip(bank.source().bank().strings()).enumerate() {
        let (Some(row), Some((twist, modes))) = (&wire.measured_torsion, measured_twist(&member.source)) else { continue };
        let length = bank.cards()[string].length_m;
        let speed = twist.wave_speed_m_s(length);
        eprintln!("measured_torsion_string={string}; row={}; twist_modes={modes}; f1_Hz={}; Z_kg_s={}; Q={} on every twist mode \
            (zeta={}), not a material law; span_m={} (nut to bridge station {}; the afterlength does not twist); c_m_s={speed}; \
            surface_inertia_kg_m={}; surface_stiffness_N={}; gas_drag=none (the measured Q is the string's own loss in air); source={:?}",
            row.id, twist.fundamental_hz(), twist.impedance_kg_s(), twist.quality_factor(), 0.5 / twist.quality_factor(),
            twist.span_m(length), twist.support_fraction(), twist.surface_inertia_kg_m(length),
            twist.impedance_kg_s() * speed, row.source);
    }
}

/// The bank's own measured specimen of one string, under any polarization or twist extension.
fn measured_specimen(source: &MaterialStringSource) -> Option<&ResolvedMeasuredStringSpecimen> {
    match source {
        MaterialStringSource::WithPolarizations(p) => measured_specimen(p.source()),
        MaterialStringSource::WithTorsion(t) => measured_specimen(t.bending()),
        MaterialStringSource::Measured(specimen) => Some(specimen),
        _ => None,
    }
}

/// The tension and bending loss factors a measured string's own `loss-factors` rows put on the string the bank holds,
/// with each row's name and the citation of its selected claim. `None` when the string is not measured or is
/// `no-intrinsic-loss`; with no tension-loss table row either, such a string has no tension loss at all.
pub(super) fn loss_receipt(ready: &Ready, string: usize) -> Option<String> {
    let member = ready.performance.bank().source().bank().strings().get(string)?;
    let specimen = measured_specimen(&member.source)?;
    let rows = specimen.loss_properties()?;
    let entry = ready.catalog.get(&ready.spec.strings.get(string)?.pack)?;
    let bound = specimen.string();
    let cited = |row: &str| specimen.measurements().property(row)
        .and_then(|property| entry.card.claims().claim(property.answer().receipt.selected))
        .map_or_else(|| "unresolved".to_owned(), |claim| claim.provenance.source.clone());
    let factor = |value: Option<f64>| value.map_or_else(|| "none".to_owned(), |value| value.to_string());
    Some(format!("eta_F={} (the string's own section measured loss-factors row {}, not the tension-loss table; source={:?}); \
        eta_B={} (row {}; source={:?}); pack={}",
        factor(bound.tension_loss_factor), rows.friction, cited(&rows.friction),
        factor(bound.bending_loss_factor), rows.bending, cited(&rows.bending), entry.pack.pack_id()))
}
