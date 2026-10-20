//! Native source selection over the original two-sided radiation owners.

use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInput;
use crate::bernoulli_aperture::two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture};

pub(super) enum Selection {
    Compact,
    Mimo(Vec<(ValveSide, usize, MimoRadiationInput)>),
}
impl Selection {
    pub(super) fn read(rows: &mut Reader<'_>, spec: &DuctRadiationSpec)
        -> Result<Self, PlateValveInputError>
    {
        if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
            != Some("radiation") {
            return Ok(Self::Compact);
        }
        let mut row = rows.row("radiation")?;
        match row.word()? {
            "compact-baffled" => { row.finish()?; Ok(Self::Compact) }
            "mimo" => {
                let count = row.count(2)?;
                row.finish()?;
                if count == 0 { return Err(bad(rows.line, "MIMO radiation needs one or two explicitly selected sides")); }
                let mut selected = [false; 2];
                let mut groups = Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row = rows.row("radiation_group")?;
                    let (side, network) = match row.word()? {
                        "upstream" => (ValveSide::Upstream, 0),
                        "downstream" => (ValveSide::Downstream, 1),
                        _ => return Err(bad(row.line, "MIMO radiation group must name upstream or downstream")),
                    };
                    if selected[network] || !spec.openings.iter().any(|o| o.network == network) {
                        return Err(bad(row.line, "each radiation group needs a distinct side with observed physical openings"));
                    }
                    selected[network] = true;
                    // Four source fields, at most one extra to trigger refusal;
                    // no unbounded token collection from a malformed header.
                    let declaration = row.fields.by_ref().take(5).collect::<Vec<_>>().join(" ");
                    let source = MimoRadiationInput::parse(&declaration)
                        .map_err(|source| PlateValveInputError::RadiationSource { line: row.line, source })?;
                    let line = row.line;
                    row.finish()?;
                    groups.push((side, line, source));
                }
                if spec.openings.iter().any(|o| !selected[o.network]) {
                    return Err(bad(rows.line, "MIMO mode must bind every observed side; mixed unbound reaction is not implicit"));
                }
                Ok(Self::Mimo(groups))
            }
            _ => Err(bad(row.line, "radiation must be compact-baffled or mimo GROUP_COUNT")),
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn build(
    selection: Selection,
    valve: LoadedAperture,
    sides: [ApertureDuctSpec; 2],
    channel: Option<SlitChannel>,
    coupling: TwoSidedLimits,
    policies: Option<[FractionalTubePolicy; 2]>,
    spec: DuctRadiationSpec,
    max_block: usize,
    gas_sources: Option<&ThermochemicalSources<'_>>,
    load: &mut impl FnMut(&str) -> Result<Vec<u8>, String>,
    gate: &CancelGate,
) -> Result<(TwoSidedPressure, Vec<BaffledOutletBinding>, Vec<(ValveSide, MimoRadiationRecord)>), PlateValveInputError> {
    let Selection::Mimo(groups) = selection else {
        let (source, bindings) = match channel {
            Some(channel) => TwoSidedPressure::from_baffled_environment_with_supply(
                valve, sides, [None, None], channel, coupling, policies, spec, max_block, gate),
            None => TwoSidedPressure::from_baffled_environment(
                valve, sides, [None, None], coupling, policies, spec, max_block, gate),
        }.map_err(|e| PlateValveInputError::Gesture(GestureCompileError::Render(e)))?;
        return Ok((source, bindings, Vec::new()));
    };
    let mut system = match channel {
        Some(channel) => TwoSidedAperture::with_material_walls_and_supply_channel(
            valve, sides, [None, None], channel, coupling, policies, gate),
        None => TwoSidedAperture::with_material_walls(
            valve, sides, [None, None], coupling, policies, gate),
    }.map_err(PlateValveInputError::Physics)?;
    let mut records = Vec::with_capacity(groups.len());
    for (side, line, input) in groups {
        let (model, record) = input.load(gas_sources, load, gate)
            .map_err(|source| PlateValveInputError::RadiationSource { line, source })?;
        let network = match side { ValveSide::Upstream => 0, ValveSide::Downstream => 1 };
        let mut local = spec.clone();
        local.openings.retain(|o| o.network == network);
        for opening in &mut local.openings { opening.network = 0; }
        system = system.with_mimo_environment(side, model, local, input.limits(), gate)
            .map_err(PlateValveInputError::Physics)?;
        records.push((side, record));
    }
    let source = TwoSidedPressure::with_radiation(system, spec, max_block)
        .map_err(|e| PlateValveInputError::Gesture(GestureCompileError::Render(e)))?;
    Ok((source, Vec::new(), records))
}
