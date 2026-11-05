//! Original source-ordinate controls of local-temperature ambient emissivity.
//! Shared card/claim identities own one control vector; each bath stays local.

use fs_blake3::ContentHash;
use fs_matdb::{ClaimId, PropertyValue};
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;

use super::{Cx, EnthalpyStepLinearization, Instant, SolveRefusal, bad, failed, finite,
    json_string, number, poll, radiation, zeros};

fn parameter_count(lowering: &radiation::LoweredRadiation) -> Result<usize, SolveRefusal> {
    let patches = lowering.surface_count();
    if !lowering.is_material_ambient() || !lowering.is_pointwise() || lowering.is_enclosure()
        || !(1..=64).contains(&patches) || lowering.patches.len() != patches
    {
        return Err(bad("material emissivity controls lost their complete ambient source layout"));
    }
    let mut count = 0usize;
    for patch in 0..patches {
        let claim = lowering.material_ambient_claim(patch)
            .ok_or_else(|| bad("material emissivity control lost its original selected claim"))?;
        let entries = match &claim.value {
            PropertyValue::Scalar { .. } => 1,
            PropertyValue::Curve { abscissa, knots, .. }
                if abscissa == fs_project::TEMPERATURE_AXIS && knots.len() >= 2 => knots.len(),
            _ => return Err(bad("local-temperature controls require original scalar or linear T-curve emissivity")),
        };
        count = count.checked_add(entries)
            .ok_or_else(|| bad("material emissivity source-control count overflows"))?;
    }
    Ok(count)
}

/// Use the native material owner's complete face blocks and source bases.
/// Counting every boundary face and each patch's ordinates is conservative
/// even when the same immutable source is shared by several disjoint traces.
pub(super) fn workspace_entries(
    nodes: usize, faces: usize, lowering: &radiation::LoweredRadiation,
) -> Result<usize, SolveRefusal> {
    if nodes == 0 || faces == 0 {
        return Err(bad("material radiation adjoints require nonempty original geometry"));
    }
    let parameters = parameter_count(lowering)?;
    nodes.checked_mul(8)
        .and_then(|n| faces.checked_mul(128).and_then(|f| n.checked_add(f)))
        .and_then(|n| parameters.checked_mul(4).and_then(|p| n.checked_add(p)))
        .and_then(|n| lowering.surface_count().checked_mul(64).and_then(|p| n.checked_add(p)))
        .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<usize>())
        .ok_or_else(|| bad("material radiation adjoint workspace overflows"))
}

/// Additional reverse storage: original scalar/ordinate controls, complete
/// material face contractions, and copied nodal source receipts in fine-grid
/// checkpoints. Forward trajectory reports have their own admission charge.
/// As elsewhere, this is an allocation allowance rather than an RSS bound.
pub(in crate::solve::transient) fn storage_bytes(
    nodes: usize, faces: usize, lowering: &radiation::LoweredRadiation, fine_steps: usize,
) -> Result<usize, SolveRefusal> {
    let base = super::storage_bytes(nodes, faces, lowering.surface_count(), fine_steps, false)?;
    let native = workspace_entries(nodes, faces, lowering)?;
    let parameters = parameter_count(lowering)?;
    base.checked_add(native.checked_mul(std::mem::size_of::<usize>())
        .ok_or_else(|| bad("material radiation adjoint storage overflows"))?)
        .and_then(|n| parameters.checked_mul(512).and_then(|p| n.checked_add(p)))
        .and_then(|n| faces.checked_mul(3).and_then(|p| p.checked_mul(fine_steps))
            .and_then(|p| p.checked_mul(4096)).and_then(|p| n.checked_add(p)))
        .filter(|&n| isize::try_from(n).is_ok())
        .ok_or_else(|| bad("material emissivity history storage allowance overflows"))
}

#[derive(PartialEq)]
enum Parameters {
    Scalar { value: f64, derivative: f64 },
    Curve { axis: String, knots: Vec<(f64, f64)>, derivatives: Vec<f64> },
}

impl Parameters {
    fn new(cx: &Cx<'_>, value: &PropertyValue, deadline: Option<(Instant, f64)>)
        -> Result<Self, SolveRefusal>
    {
        match value {
            PropertyValue::Scalar { value, .. } =>
                Ok(Self::Scalar { value: finite(*value)?, derivative: 0.0 }),
            PropertyValue::Curve { abscissa, knots, .. }
                if abscissa == fs_project::TEMPERATURE_AXIS && knots.len() >= 2 =>
            {
                let mut original = Vec::new();
                original.try_reserve_exact(knots.len())
                    .map_err(|_| bad("material emissivity source-control allocation refused"))?;
                for (index, &(temperature, emissivity)) in knots.iter().enumerate() {
                    if index % 256 == 0 { poll(cx, deadline)?; }
                    original.push((finite(temperature)?, finite(emissivity)?));
                }
                Ok(Self::Curve { axis: abscissa.clone(), knots: original,
                    derivatives: zeros(knots.len())? })
            }
            _ => Err(bad("material emissivity controls lost their original scalar or T curve")),
        }
    }

    fn add(&mut self, cx: &Cx<'_>, gradient: &MaterialEmissivityGradient,
        deadline: Option<(Instant, f64)>) -> Result<(), SolveRefusal>
    {
        match (self, gradient) {
            (Self::Scalar { derivative, .. }, MaterialEmissivityGradient::Scalar { derivative: value }) => {
                *derivative = finite(*derivative + finite(*value)?)?;
            }
            (Self::Curve { axis, knots, derivatives }, MaterialEmissivityGradient::Curve {
                axis: original_axis, knot_abscissae, knot_value_derivatives,
            }) if axis == original_axis && knots.len() == knot_abscissae.len()
                && derivatives.len() == knot_value_derivatives.len() =>
            {
                for index in 0..knots.len() {
                    if index % 256 == 0 { poll(cx, deadline)?; }
                    if knots[index].0.to_bits() != knot_abscissae[index].to_bits() {
                        return Err(bad("material emissivity pullback changed an original knot abscissa"));
                    }
                    derivatives[index] = finite(derivatives[index] + finite(knot_value_derivatives[index])?)?;
                }
            }
            _ => return Err(bad("material emissivity pullback changed its original source-control layout")),
        }
        Ok(())
    }

    fn fragment(&self, cx: &Cx<'_>, deadline: Option<(Instant, f64)>) -> Result<String, SolveRefusal> {
        match self {
            Self::Scalar { value, derivative } => Ok(format!(
                "\"model\":\"scalar\",\"emissivity\":{},\"dtemperature_dsource_emissivity_k\":{}",
                number(*value)?, number(*derivative)?,
            )),
            Self::Curve { axis, knots, derivatives } => {
                let rows = knots.iter().zip(derivatives).enumerate().map(|(index, (&(temperature, value), &derivative))| {
                    if index % 256 == 0 { poll(cx, deadline)?; }
                    Ok(format!("{{\"index\":{index},\"temperature_k\":{},\"emissivity\":{},\"dtemperature_demissivity_ordinate_k\":{}}}",
                        number(temperature)?, number(value)?, number(derivative)?)) })
                    .collect::<Result<Vec<_>, SolveRefusal>>()?;
                Ok(format!("\"model\":\"temperature-curve\",\"axis\":{},\"knots\":[{}]",
                    json_string(axis), rows.join(",")))
            }
        }
    }
}

struct Source {
    card: ContentHash,
    claim: ClaimId,
    material_state: String,
    surfaces: usize,
    parameters: Parameters,
}

struct Surface {
    target: String,
    source: usize,
    fields: String,
    reservoir: f64,
}

pub(super) struct Controls {
    sources: Vec<Source>,
    surfaces: Vec<Surface>,
    steps: usize,
}

impl Controls {
    pub(super) fn new(cx: &Cx<'_>, lowering: &radiation::LoweredRadiation,
        deadline: Option<(Instant, f64)>) -> Result<Self, SolveRefusal>
    {
        parameter_count(lowering)?;
        let mut sources: Vec<Source> = Vec::new();
        let mut surfaces = Vec::with_capacity(lowering.surface_count());
        for (index, ((name, target), patch)) in lowering.surface_labels().zip(&lowering.patches).enumerate() {
            poll(cx, deadline)?;
            if target != patch.region() {
                return Err(bad("material emissivity controls changed the declared trace order"));
            }
            let emissivity = patch.emissivity();
            let (card, claim) = (emissivity.card_identity(), emissivity.receipt().selected);
            let original = lowering.material_ambient_claim(index)
                .ok_or_else(|| bad("material emissivity controls lost the pinned original source"))?;
            let parameters = Parameters::new(cx, &original.value, deadline)?;
            let source = if let Some(source) = sources.iter().position(|s| s.card == card && s.claim == claim) {
                if sources[source].parameters != parameters || sources[source].material_state != emissivity.material_state() {
                    return Err(bad("one material emissivity identity acquired different original controls"));
                }
                sources[source].surfaces += 1;
                source
            } else {
                sources.push(Source { card, claim, material_state: emissivity.material_state().into(),
                    surfaces: 1, parameters });
                sources.len() - 1
            };
            let fields = format!(
                "\"name\":{},\"target\":{},\"source_index\":{source},\"card\":{},\"claim\":{},\"material_state\":{},\"source_query_emissivity\":{},\"source_query_temperature_k\":{},\"emissivity_receipt\":{},\"reservoir_temperature_k\":{}",
                json_string(name), json_string(target), json_string(&card.to_hex()), json_string(&claim.0.to_hex()),
                json_string(emissivity.material_state()), number(emissivity.value())?, number(emissivity.temperature_k())?,
                json_string(&emissivity.receipt().content_hash().to_hex()), number(patch.ambient_temperature_k())?,
            );
            surfaces.push(Surface { target: target.into(), source, fields, reservoir: 0.0 });
        }
        Ok(Self { sources, surfaces, steps: 0 })
    }

    pub(super) fn add(&mut self, cx: &Cx<'_>, step: &EnthalpyStepLinearization<'_>,
        nodal_load_adjoint: &[f64], deadline: Option<(Instant, f64)>) -> Result<(), SolveRefusal>
    {
        // P_free*dt*lambda is already supplied by the original transport owner.
        // Outward radiative heat has the opposite sign; never introduce dt again.
        let mut outward = zeros(nodal_load_adjoint.len())?;
        for (vertex, (&load, weight)) in nodal_load_adjoint.iter().zip(&mut outward).enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            *weight = finite(-load)?;
        }
        let gradient = step.material_ambient_parameter_pullback(cx, &outward).map_err(failed)?;
        if gradient.patches.len() != self.surfaces.len() {
            return Err(bad("material emissivity pullback lost its complete patch layout"));
        }
        for (surface, patch) in self.surfaces.iter_mut().zip(&gradient.patches) {
            poll(cx, deadline)?;
            let source = &mut self.sources[surface.source];
            if patch.region != surface.target || patch.material_identity != source.card || patch.claim != source.claim {
                return Err(bad("material emissivity pullback changed the original source or trace identity"));
            }
            source.parameters.add(cx, &patch.emissivity, deadline)?;
            surface.reservoir = finite(surface.reservoir + finite(patch.reservoir_temperature)?)?;
        }
        self.steps = self.steps.checked_add(1).ok_or_else(|| bad("material radiation reverse prefix overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn fragment(&self, cx: &Cx<'_>, expected_steps: usize, max_radiation_entries: usize,
        workspace_bytes: usize, deadline: Option<(Instant, f64)>) -> Result<String, SolveRefusal>
    {
        if self.steps != expected_steps {
            return Err(bad("material emissivity controls do not cover the selected reverse prefix"));
        }
        let sources = self.sources.iter().map(|source| Ok(format!(
            "{{\"card\":{},\"claim\":{},\"material_state\":{},\"shared_surface_count\":{},{} }}",
            json_string(&source.card.to_hex()), json_string(&source.claim.0.to_hex()),
            json_string(&source.material_state), source.surfaces, source.parameters.fragment(cx, deadline)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let surfaces = self.surfaces.iter().map(|surface| Ok(format!(
            "{{{},\"dtemperature_dreservoir_temperature_k_per_k\":{}}}", surface.fields, number(surface.reservoir)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        poll(cx, deadline)?;
        Ok(format!(
            ",\"radiation\":{{\"method\":\"complete-pointwise-p1-material-ambient-enthalpy-tangent\",\"spatial_model\":\"pointwise\",\"emissivity_response\":\"local-temperature\",\"patch_count\":{},\"control_steps\":{},\"max_radiation_workspace_entries\":{max_radiation_entries},\"additional_workspace_bytes\":{workspace_bytes},\"parameters\":\"original-shared-emissivity-source-values-and-independent-reservoir-temperatures-throughout-selected-trajectory\",\"parameter_pullback\":\"prepared-material-pointwise-outward-load-contraction\",\"parameter_krylov_iterations\":0,\"parameter_linear_solves\":0,\"residual_target\":\"original-accepted-forward-physical-tolerance\",\"source_feasible_directions\":\"preserve-all-queried-local-emissivities-in-zero-exclusive-one-inclusive\",\"emissivity_sources\":[{}],\"surfaces\":[{}],\"scope\":{}}}",
            self.surfaces.len(), self.steps, sources.join(","), surfaces.join(","), json_string(SCOPE),
        ))
    }
}

pub(super) const SCOPE: &str = "Local nodal temperature queries the original pinned scalar or linear T-curve emissivity before P1 interpolation. The complete emission and source-slope feedback enter the same accepted enthalpy transpose as storage, conductivity, live contact and physical convection. Original scalar values or curve ordinates are independent controls at fixed abscissae; contributions sharing the exact material-card and selected-claim identity are summed once across physical surfaces and the selected history prefix. Per-surface black-bath temperatures remain independent. Parameter contractions use minus the already masked time-scaled thermal multiplier, with no additional dt or parameter solve. Initial selected peaks retain original identities and zero future derivatives. Source selection, source abscissae, nonthermal context, caloric chart, geometry and time grid stay fixed. Every local physical emissivity remains in (0,1]; source or chart corners and unsupported states refuse classical derivatives. No per-patch frozen-emissivity derivative, phase-dependent emissivity source, enclosure response, uncertainty propagation or physical-validation certificate is supplied by this local-temperature ambient adapter.";
