//! Native source selection over the original two-sided radiation owners.

use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInput;
use crate::bernoulli_aperture::radiation::environment::multiport::source::{
    ComputedBaffledInput, MimoSourceInput, MimoSourceRecord,
};
use crate::bernoulli_aperture::two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture};
use crate::bernoulli_aperture::wall::material::MaterialWallInput;
use crate::bernoulli_aperture::radiation::environment::bind_baffled_environment;
use crate::pcm_wav::observation::ObservationPolicy;

pub(super) enum Selection {
    Compact,
    Mimo(Vec<(ValveSide, usize, MimoSourceInput)>),
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
            kind @ ("mimo" | "computed-baffled") => {
                let computed = kind == "computed-baffled";
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
                    // Bounded source fields, at most one extra to trigger refusal;
                    // no unbounded token collection from a malformed header.
                    let fields = if computed { ComputedBaffledInput::FIELDS } else { 4 };
                    let declaration = row.fields.by_ref().take(fields + 1).collect::<Vec<_>>().join(" ");
                    let source = if computed {
                        ComputedBaffledInput::parse(&declaration).map(MimoSourceInput::ComputedBaffled)
                    } else {
                        MimoRadiationInput::parse(&declaration).map(MimoSourceInput::Artifact)
                    }
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
            _ => Err(bad(row.line, "radiation must be compact-baffled, mimo GROUP_COUNT or computed-baffled GROUP_COUNT")),
        }
    }
}

/// Cold source construction precedes material-key attachment; the original
/// receiver is admitted only after every moving boundary is present.
pub(super) struct Prepared {
    pub(super) system: TwoSidedAperture,
    spec: DuctRadiationSpec,
    policy: ObservationPolicy,
    bindings: Vec<BaffledOutletBinding>,
    records: Vec<(ValveSide, MimoRadiationRecord)>,
    computed: Vec<(ValveSide, ComputedBaffledRecord)>,
}
impl Prepared {
    #[allow(clippy::type_complexity)]
    pub(super) fn observe(self, max_block: usize) -> Result<(TwoSidedPressure,
        Vec<BaffledOutletBinding>, Vec<(ValveSide, MimoRadiationRecord)>,
        Vec<(ValveSide, ComputedBaffledRecord)>), PlateValveInputError>
    {
        let source = TwoSidedPressure::with_observation_policy(
            self.system, self.spec, max_block, self.policy, &self.bindings)
            .map_err(|e| PlateValveInputError::Gesture(GestureCompileError::Render(e)))?;
        Ok((source, self.bindings, self.records, self.computed))
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    selection: Selection,
    valve: LoadedAperture,
    sides: [ApertureDuctSpec; 2],
    walls: [Option<MaterialWallInput>; 2],
    channel: Option<FixedChannel>,
    coupling: TwoSidedLimits,
    policies: Option<[FractionalTubePolicy; 2]>,
    spec: DuctRadiationSpec,
    exterior: &GasRecord,
    max_block: usize,
    gas_sources: Option<&ThermochemicalSources<'_>>,
    load: &mut impl FnMut(&str) -> Result<Vec<u8>, String>,
    gate: &CancelGate,
) -> Result<Prepared, PlateValveInputError> {
    checkpoint(gate)?;
    let rate = spec.reconstruction.mechanics_rate_hz();
    if rate == 0 || max_block == 0 || valve.accepted_steps() != 0
        || valve.time_step_s().to_bits() != (1.0 / f64::from(rate)).to_bits() {
        return Err(bad(0, "radiation needs the original unadvanced valve and matching mechanical clock"));
    }
    let (sides, spec, bindings, policy, groups) = match selection {
        Selection::Compact => {
            if policies.is_some_and(|ps| ps.iter().any(|p|
                p.maximum_frequency_hz < spec.receiver.maximum_frequency_hz)) {
                return Err(bad(0, "fractional propagation must cover the exterior radiation band"));
            }
            let [up, down] = sides;
            let mut graphs = [up.graph, down.graph];
            let environment = bind_baffled_environment(
                &mut graphs, &[&up.gases, &down.gases], spec, gate)
                .map_err(PlateValveInputError::Physics)?;
            let [up_graph, down_graph] = graphs;
            let sides = [ApertureDuctSpec { graph: up_graph, ..up },
                ApertureDuctSpec { graph: down_graph, ..down }];
            let (spec, bindings) = environment.into_parts();
            (sides, spec, bindings, ObservationPolicy::Acceptance, Vec::new())
        }
        Selection::Mimo(groups) => (sides, spec, Vec::new(), ObservationPolicy::Diagnostic, groups),
    };
    let mut system = match channel {
        Some(channel) => TwoSidedAperture::with_material_walls_and_supply_channel(
            valve, sides, walls, channel, coupling, policies, gate),
        None => TwoSidedAperture::with_material_walls(
            valve, sides, walls, coupling, policies, gate),
    }.map_err(PlateValveInputError::Physics)?;
    let mut records = Vec::with_capacity(groups.len());
    let mut computed = Vec::with_capacity(groups.len());
    for (side, line, input) in groups {
        let network = match side { ValveSide::Upstream => 0, ValveSide::Downstream => 1 };
        let mut local = spec.clone();
        local.openings.retain(|o| o.network == network);
        for opening in &mut local.openings { opening.network = 0; }
        let original = system.source(side);
        let (model, record) = input.build(&original.graph, &original.gases, &local,
            exterior, gas_sources, load, gate)
            .map_err(|source| PlateValveInputError::RadiationSource { line, source })?;
        system = system.with_mimo_environment(side, model, local, input.limits(), gate)
            .map_err(PlateValveInputError::Physics)?;
        match record {
            MimoSourceRecord::Artifact(record) => records.push((side, record)),
            MimoSourceRecord::ComputedBaffled(record) => computed.push((side, record)),
        }
    }
    checkpoint(gate)?;
    Ok(Prepared { system, spec, policy, bindings, records, computed })
}
