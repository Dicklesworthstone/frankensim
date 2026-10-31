//! Source-derived reciprocal links on the original full inductance matrix.

use fs_material::field_storage::{
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor, linear_magnetic_permeability,
};
use fs_phs::circuit::MutualInductance;

use super::{MaterialCircuitError, MaterialStorage, MaterialStorageElement, invalid};

const MAX_BRANCHES: usize = super::super::MAX_MATERIAL_CIRCUIT_BRANCHES;
const MAX_LINKS: usize = MAX_BRANCHES * (MAX_BRANCHES - 1) / 2;

/// Authored shared magnetic geometry on two original inductor branches.
///
/// Both winding materials must resolve the same absolute linear permeability
/// source at the same physical query. The mutual coefficient is
/// `M = geometry_m * mu + constant_h`; this metadata does not add another field
/// energy store. The existing descriptor owns the complete coupled energy and
/// admits positive definiteness of the full matrix, including every held link.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialMagneticCoupling {
    /// Original graph branch indices, with polarity inherited from terminals.
    pub branches: (usize, usize),
    /// Signed shared-path factor [m], e.g. oriented `N1 * N2 * A / length`.
    /// Fixed during source refresh; it is not inferred from equal material cards.
    pub geometry_m: f64,
    /// Authored temperature-independent contribution to this same link [H].
    pub constant_h: f64,
}

fn pair(branches: (usize, usize)) -> Result<(usize, usize), MaterialCircuitError> {
    let (a, b) = branches;
    if a == b || a >= MAX_BRANCHES || b >= MAX_BRANCHES {
        return Err(invalid(
            "magnetic coupling needs two distinct bounded branch indices",
        ));
    }
    Ok((a.min(b), a.max(b)))
}

pub(super) fn canonical_couplings(
    mut couplings: Vec<MaterialMagneticCoupling>,
) -> Result<Vec<MaterialMagneticCoupling>, MaterialCircuitError> {
    if couplings.len() > MAX_LINKS {
        return Err(invalid(
            "source-bound magnetic links exceed bounded pair capacity",
        ));
    }
    for coupling in &mut couplings {
        coupling.branches = pair(coupling.branches)?;
        if !coupling.geometry_m.is_finite() || !coupling.constant_h.is_finite() {
            return Err(invalid(
                "magnetic geometry and constant mutual offset must be finite",
            ));
        }
    }
    couplings.sort_by_key(|coupling| coupling.branches);
    if couplings.windows(2).any(|p| p[0].branches == p[1].branches) {
        return Err(invalid("duplicate source-bound magnetic pair"));
    }
    Ok(couplings)
}

pub(super) fn complete_mutual(
    mut held: Vec<MutualInductance>,
    couplings: &[MaterialMagneticCoupling],
    storage: &[MaterialStorage],
) -> Result<Vec<MutualInductance>, MaterialCircuitError> {
    if held.len() > MAX_LINKS || couplings.len() > MAX_LINKS - held.len() {
        return Err(invalid(
            "complete magnetic links exceed bounded pair capacity",
        ));
    }
    held.extend(couplings.iter().map(|coupling| MutualInductance {
        branches: coupling.branches,
        henries: 0.0,
    }));
    // The shared refresh admission also rejects a held link that duplicates a
    // governed pair, including an opposite ordering of the same two branches.
    refresh_mutual(&held, couplings, storage)
}

fn common_permeability(
    a: &ResolvedInductor,
    b: &ResolvedInductor,
) -> Result<f64, MaterialCircuitError> {
    let a = a.material();
    let b = b.material();
    let mu = linear_magnetic_permeability(a)?;
    let other = linear_magnetic_permeability(b)?;
    let first = a
        .property(MAGNETIC_PERMEABILITY_PROPERTY)
        .ok_or_else(|| invalid("first winding lost its permeability source"))?;
    let second = b
        .property(MAGNETIC_PERMEABILITY_PROPERTY)
        .ok_or_else(|| invalid("second winding lost its permeability source"))?;
    if a.card_identity() != b.card_identity()
        || a.material() != b.material()
        || a.query_point().len() != b.query_point().len()
        || a.query_point()
            .iter()
            .zip(b.query_point())
            .any(|((na, va), (nb, vb))| na != nb || va.to_bits() != vb.to_bits())
        || a.axis_quantities() != b.axis_quantities()
        || first.requirement() != second.requirement()
        || first.source_axes() != second.source_axes()
        || first.answer().receipt.selected != second.answer().receipt.selected
        || mu.to_bits() != other.to_bits()
    {
        return Err(invalid(
            "coupled windings must share the actual permeability source, card and typed physical query",
        ));
    }
    // Receipt selection policies may honestly differ after a pinned refresh.
    // Other properties retained by either winding do not define this shared law.
    Ok(mu)
}

/// Refresh only governed values, preserving complete link order and orientation.
/// The caller may supply just the participating storage owners; unrelated
/// capacitors, inductors and held mutual links need no material replay here.
pub(crate) fn refresh_mutual(
    links: &[MutualInductance],
    couplings: &[MaterialMagneticCoupling],
    storage: &[MaterialStorage],
) -> Result<Vec<MutualInductance>, MaterialCircuitError> {
    if links.len() > MAX_LINKS || couplings.len() > MAX_LINKS || storage.len() > MAX_BRANCHES {
        return Err(invalid(
            "magnetic refresh exceeds bounded link or source capacity",
        ));
    }
    if couplings.is_empty() {
        return Ok(links.to_vec());
    }
    let couplings = canonical_couplings(couplings.to_vec())?;
    let mut indices = Vec::with_capacity(links.len());
    for (index, link) in links.iter().enumerate() {
        if !link.henries.is_finite() {
            return Err(invalid("mutual inductance must be finite"));
        }
        indices.push((pair(link.branches)?, index));
    }
    indices.sort_unstable_by_key(|(branches, _)| *branches);
    if indices.windows(2).any(|p| p[0].0 == p[1].0) {
        return Err(invalid("duplicate held or governed mutual inductance pair"));
    }
    let mut sources: Vec<_> = storage.iter().collect();
    sources.sort_unstable_by_key(|source| source.branch);
    if sources.iter().any(|source| source.branch >= MAX_BRANCHES)
        || sources.windows(2).any(|p| p[0].branch == p[1].branch)
    {
        return Err(invalid(
            "magnetic refresh requires unique bounded storage branches",
        ));
    }
    let inductor = |branch| -> Result<&ResolvedInductor, MaterialCircuitError> {
        let index = sources
            .binary_search_by_key(&branch, |source| source.branch)
            .map_err(|_| invalid("governed mutual link is missing a winding source"))?;
        match &sources[index].element {
            MaterialStorageElement::Inductor(value) => Ok(value),
            MaterialStorageElement::Capacitor(_) => Err(invalid(
                "governed magnetic links must name source-bound inductors",
            )),
        }
    };
    let mut refreshed = links.to_vec();
    for coupling in couplings {
        let index = indices
            .binary_search_by_key(&coupling.branches, |(branches, _)| *branches)
            .map_err(|_| invalid("governed magnetic pair is absent from complete mutual links"))?;
        let mu = common_permeability(
            inductor(coupling.branches.0)?,
            inductor(coupling.branches.1)?,
        )?;
        let henries = coupling.geometry_m.mul_add(mu, coupling.constant_h);
        if !henries.is_finite()
            || coupling.geometry_m != 0.0 && coupling.constant_h == 0.0 && henries == 0.0
        {
            return Err(invalid(
                "source-derived mutual inductance is unrepresentable",
            ));
        }
        refreshed[indices[index].1].henries = henries;
    }
    Ok(refreshed)
}
