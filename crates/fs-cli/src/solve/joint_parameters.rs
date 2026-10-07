//! Declared manufactured-joint bands through the ordinary native thermal solve.
use super::{BindingRequirements, BindingTarget, CardPackSet, EvidenceWork, ProjectSpec,
    PropagatedTerm, SolveRefusal, conduction_error, resolve_bindings};
use fs_matdb::{ClaimId, ClaimSelection, SelectionPolicy};
use fs_project::InterfaceState;
use fs_project::interface_state::band::{ResistanceBand, resistance_band};

const MAX_VERTICES: usize = 64;
const MAX_SOURCE_CANDIDATES: usize = 1024;

fn width(state: &InterfaceState) -> f64 {
    match state {
        InterfaceState::DryContact { pressure_half_width, .. } => pressure_half_width.value,
        InterfaceState::Tim { thickness_half_width, .. }
        | InterfaceState::Adhesive { thickness_half_width, .. } => thickness_half_width.value,
        InterfaceState::GapWithFluid { gap_half_width, .. } => gap_half_width.value,
        InterfaceState::BoltedWithPattern { torque_half_width, .. } => torque_half_width.value,
    }
}
fn gap(reason: impl Into<String>) -> PropagatedTerm { PropagatedTerm::Unmeasured { reason: reason.into() } }
fn poll(work: EvidenceWork<'_>) -> Result<(), SolveRefusal> {
    if work.is_requested() { Err(conduction_error("cli-solve-cancelled",
        "manufactured-joint parameter propagation interrupted", "resume the retained pipeline prefix")) }
    else { Ok(()) }
}
fn corner_count(varying: usize, conductivity: bool) -> Option<usize> {
    u32::try_from(varying).ok().and_then(|n| 1usize.checked_shl(n))
        .and_then(|n| n.checked_mul(if conductivity { 2 } else { 1 }))
        .filter(|n| *n <= MAX_VERTICES)
}

/// Extend the legacy conductivity calculation only when a nonzero joint band
/// is explicitly declared. Enumerate ALL combinations of resistance extrema,
/// crossed with the existing common conductivity lower/upper perturbations.
/// Admit the complete joint design BEFORE additional physical work; no partial maximum.
/// A card query/physical failure leaves the term unmeasured; cancellation and
/// invocation-work errors propagate through the caller's existing vertex seam.
/// These deterministic corners assume monotone thermal response in resistance,
/// not in manufactured coordinates. They do not bound unprovided card errors,
/// geometry, other inputs or parameter/ambient interactions, or define a PDF.
#[allow(clippy::too_many_arguments)]
pub(super) fn propagate(
    base: &ProjectSpec, cards: &CardPackSet, nominal: f64, conductivity: PropagatedTerm,
    work: EvidenceWork<'_>,
    mut evaluate: impl FnMut(String, &ProjectSpec, f64) -> Result<Result<(String, f64), String>, SolveRefusal>,
) -> Result<PropagatedTerm, SolveRefusal> {
    let declarations = base.interface_cards.as_deref().unwrap_or(&[]);
    if !declarations.iter().any(|row| width(&row.state) != 0.0) { return Ok(conductivity); }
    poll(work)?;
    if declarations.len() > MAX_VERTICES { return Ok(gap("joint-band declaration allowance exhausted")); }
    let has_conductivity = base.materials.iter().flatten().any(|b| b.conductivity_tolerance.is_some());
    // A failed declared conductivity solve cannot be repaired by reporting only
    // the successful contact subproblem. With no conductivity tolerance, this
    // is explicitly a joint-only sensitivity at the nominal material law.
    if has_conductivity && conductivity.half_width().is_none() { return Ok(conductivity); }
    let library = cards.library();
    let requirements = BindingRequirements::thermal_steady_v1();
    let resolved = resolve_bindings(base, &library, &requirements);
    poll(work)?;
    if !resolved.admissible() { return Ok(gap("joint-band propagation cannot reproduce the nominal material bindings")); }
    let mut bands: Vec<(usize, ResistanceBand)> = Vec::new();
    let mut detail = Vec::new();
    for (i, declared) in declarations.iter().enumerate() {
        poll(work)?;
        if width(&declared.state) == 0.0 { continue; }
        let Some(binding) = resolved.bindings.iter().find(|b|
            matches!(&b.target, BindingTarget::Interface(name) if name == &declared.interface))
            else { return Ok(gap("joint-band propagation lost its interface binding")); };
        let Some(card) = library.interface(&declared.card)
            else { return Ok(gap("joint-band propagation lost its original card")); };
        let selection = match declared.claim.as_deref() {
            Some(text) => match fs_blake3::ContentHash::from_hex(text) {
                Some(pin) => ClaimSelection::Pinned(ClaimId(pin)),
                None => return Ok(gap("invalid manufactured-joint claim pin")),
            },
            None => ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        };
        let evaluated = resistance_band(card.claims(), &declared.state,
            &requirements.temperature_axis, binding.range_lo, binding.range_hi,
            selection, MAX_SOURCE_CANDIDATES);
        poll(work)?;
        let band = match evaluated {
            Ok(band) => band,
            Err(error) => return Ok(gap(format!("joint `{}` band refused: {}", declared.interface, error.detail))),
        };
        let Some(property) = binding.properties.iter().find(|p| p.property == fs_project::CONTACT_RESISTANCE_PROPERTY)
            else { return Ok(gap("missing resolved nominal contact resistance")); };
        if band.claim.0.to_hex() != property.selected_claim
            || band.nominal_resistance.to_bits() != property.value_lo.to_bits()
        { return Ok(gap("joint band does not reproduce its nominal source and resistance")); }
        detail.push(format!("`{}` {} [{}, {}] {} -> R'' [{}, {}] m2 K/W",
            declared.interface, band.axis, band.coordinate_low, band.coordinate_high, band.unit,
            band.minimum_resistance, band.maximum_resistance));
        if band.minimum_resistance != band.maximum_resistance { bands.push((i, band)); }
    }
    poll(work)?;
    let Some(count) = corner_count(bands.len(), has_conductivity)
        else { return Ok(gap("joint/conductivity corner design exceeds 64 physical re-solves; no partial parameter bound was produced")); };
    let (mut half_width, mut vertices, conductivity_detail) = match conductivity {
        PropagatedTerm::Measured { half_width_k, detail, vertices, .. } => (half_width_k, vertices, detail),
        PropagatedTerm::Unmeasured { .. } => (0.0, Vec::new(),
            "no conductivity tolerance declared: material laws held at their nominal values; unprovided material uncertainty is not bounded".into()),
    };
    if !nominal.is_finite() || !half_width.is_finite() || half_width < 0.0 {
        return Ok(gap("nonfinite nominal value or conductivity width in joint-band propagation"));
    }
    // Flat selected laws cannot change the contact operator. Preserve legacy
    // conductivity results without spending duplicate physical solves.
    if !bands.is_empty() {
        let sides: &[f64] = if has_conductivity { &[-1.0, 1.0] } else { &[0.0] };
        for mask in 0..count / sides.len() {
            poll(work)?;
            let mut vertex = base.clone();
            let mut labels = Vec::new();
            for (bit, (i, band)) in bands.iter().enumerate() {
                let high = mask & (1usize << bit) != 0;
                let row = &mut vertex.interface_cards.as_mut().expect("declared interfaces")[*i];
                row.state = if high { band.maximum_state.clone() } else { band.minimum_state.clone() };
                labels.push(format!("{} {}", row.interface, row.state.render()));
            }
            for &side in sides {
                poll(work)?;
                let label = format!("joint corner {mask}: {}; conductivity side {side}", labels.join("; "));
                match evaluate(label, &vertex, side)? {
                    Ok(row) if row.1.is_finite() && (row.1 - nominal).is_finite() => {
                        half_width = half_width.max((row.1 - nominal).abs());
                        vertices.push(row);
                    }
                    Ok(_) => return Ok(gap("joint parameter vertex produced nonfinite thermal deviation")),
                    Err(reason) => return Ok(gap(format!("joint parameter corner refused; incomplete design is not a bound: {reason}"))),
                }
            }
        }
    }
    poll(work)?;
    Ok(PropagatedTerm::Measured { half_width_k: half_width, method: "joint-state-resistance-corner-resolve",
        detail: format!("{}; {conductivity_detail}; all independent contact resistance extrema crossed with the existing common conductivity perturbations; source knots included when finding resistance extrema; internal vertices are point states, not bands applied twice; nominal project and cards unchanged; assumes monotone thermal response in each resistance; no probability law, validated interval, missing-source bound or parameter/ambient interaction certificate",
            detail.join("; ")), vertices })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complete_corner_design_is_admitted_before_any_subset_can_be_published() {
        assert_eq!(corner_count(0, false), Some(1));
        assert_eq!(corner_count(5, true), Some(64));
        assert_eq!(corner_count(6, false), Some(64));
        assert_eq!(corner_count(6, true), None);
        assert_eq!(corner_count(usize::MAX, false), None);
        let mut choices = std::collections::BTreeSet::new();
        for mask in 0..corner_count(3, false).unwrap() {
            choices.insert((0..3).map(|bit| mask & (1usize << bit) != 0).collect::<Vec<_>>());
        }
        assert_eq!(choices.len(), 8);
        assert!(choices.contains(&vec![false, true, false]), "mixed corners must not be omitted");
    }
}
