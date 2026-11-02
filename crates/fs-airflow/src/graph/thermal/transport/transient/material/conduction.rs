//! Source-backed axial fluid conduction on the original finite pipe inventory.
//!
//! The two neighbouring cell centres exchange heat through their actual half
//! lengths and cross-sectional areas: `R_i = L_i / (2 k_i A_i)`. The original
//! source state at each cell's own bulk temperature supplies `k_i`; it is held
//! within that finite cell, not averaged into a new material or sampled at a
//! fictitious common temperature. Native reciprocal links own both heat debits,
//! storage, cancellation and integration. Advection remains its separate owner.
//!
//! The declared link must cover the whole shared section, with area equal to
//! the smaller adjacent cross section. A stepped-area one-dimensional model
//! assumes no additional constriction/spreading/contact resistance. Smaller
//! partial apertures, junctions and wall links need their own explicit geometry
//! and law and cannot borrow this adapter. No turbulent dispersion is inferred.

use super::*;
use std::collections::BTreeSet;

/// An already authored native thermal link across one whole internal pipe face.
/// Both endpoints must be consecutive cells of this original hydraulic branch.
/// The link's authored conductance must be zero; its source has one owner.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GraphConductiveLink {
    pub link: u64,
    pub branch: usize,
}

/// Actual geometry, bulk transport and signed heat at one thermal candidate.
/// Endpoint arrays follow the native link's `a,b` ordering, independently of
/// branch flow direction. Source owners are retained by `source_sections()`.
#[derive(Clone, Debug)]
pub struct GraphConductionSample {
    pub link: u64,
    pub branch: usize,
    pub cells: [CellId; 2],
    /// Original section indices corresponding to these two native endpoints.
    pub sections: [usize; 2],
    pub temperatures_k: [f64; 2],
    pub lengths_m: [f64; 2],
    pub areas_m2: [f64; 2],
    pub link_area_m2: f64,
    pub transport: [FluidTransportState; 2],
    pub half_resistances_k_w: [f64; 2],
    pub conductance_w_k: f64,
    /// Positive from native endpoint a to b, never an external heat source.
    pub heat_a_to_b_w: f64,
}

#[derive(Debug)]
pub enum GraphConductionError {
    Source(SourceGraphError),
    Link { link: u64, reason: &'static str },
    Thermal { link: u64, source: DepositionError },
}
impl core::fmt::Display for GraphConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite graph axial conduction: {self:?}")
    }
}
impl std::error::Error for GraphConductionError {}

/// Native coefficient overrides sharing the exact immutable source epoch.
/// Sharing retains original liquid/gas receipts without re-querying or cloning
/// them once per face, and permits `(sources, conduction)` in a native payload.
#[derive(Debug)]
pub struct ResolvedGraphConduction<M> {
    conductances_w_k: BTreeMap<u64, f64>,
    samples: Vec<GraphConductionSample>,
    sections: Arc<Vec<Vec<DistributedMaterialSection<M>>>>,
}
impl<M> ResolvedGraphConduction<M> {
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }
    /// Samples are in native link-ID order, independent of binding order.
    pub fn samples(&self) -> &[GraphConductionSample] {
        &self.samples
    }
    /// The exact source section owners used by `GraphSourceState`, shared by
    /// reference count. Sample branch/section indices address this same epoch.
    pub fn source_sections(&self) -> &[Vec<DistributedMaterialSection<M>>] {
        &self.sections
    }

    /// Apply once to the native thermal residual's link conductance map.
    /// Conflicting ownership or cancellation leaves `destination` unchanged.
    /// Do not also add the diagnostic signed powers as cell heat sources.
    pub fn apply_to(
        &self,
        cx: &Cx<'_>,
        destination: &mut BTreeMap<u64, f64>,
    ) -> Result<(), GraphConductionError> {
        for &link in self.conductances_w_k.keys() {
            poll(cx).map_err(GraphConductionError::Source)?;
            if destination.contains_key(&link) {
                return Err(GraphConductionError::Link {
                    link,
                    reason: "axial link already has a conductance owner",
                });
            }
        }
        poll(cx).map_err(GraphConductionError::Source)?;
        destination.extend(self.conductances_w_k.iter().map(|(&id, &g)| (id, g)));
        Ok(())
    }
}

/// Resolve held-cell axial conductivity at one exact native thermal candidate.
/// `maximum_links` bounds all active links inspected for duplicate whole-face
/// ownership, as well as the selected bindings. The native thermal timeline
/// retains its own independent cell/link and nonlinear/step work limits.
pub(super) fn evaluate_links<M>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    sources: &GraphSourceState<M>,
    bindings: &[GraphConductiveLink],
    maximum_links: usize,
) -> Result<ResolvedGraphConduction<M>, GraphConductionError> {
    let source_error = GraphConductionError::Source;
    poll(cx).map_err(source_error)?;
    if !(1..=1_000_000).contains(&maximum_links)
        || bindings.len() > maximum_links
        || thermal.links().len() > maximum_links
    {
        return Err(source_error(SourceGraphError::Invalid(
            "bounded active and selected axial link counts required",
        )));
    }
    core.admit_sources(cx, thermal, sources)
        .map_err(source_error)?;
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.link);
    let mut ids = BTreeSet::new();
    let mut owned_pairs = BTreeMap::new();
    let mut samples = Vec::with_capacity(bindings.len());
    let mut conductances_w_k = BTreeMap::new();
    for binding in ordered {
        poll(cx).map_err(source_error)?;
        let invalid = |reason| GraphConductionError::Link {
            link: binding.link,
            reason,
        };
        if !ids.insert(binding.link) {
            return Err(invalid("duplicate axial conduction binding"));
        }
        let link = thermal
            .links()
            .get(&binding.link)
            .ok_or_else(|| invalid("inactive axial conduction link"))?;
        if link.conductance_w_k != 0.0 {
            return Err(invalid(
                "axial source links require zero authored conductance",
            ));
        }
        let mut sections = [0; 2];
        for (side, cell) in [link.a, link.b].into_iter().enumerate() {
            let Some(Some((branch, section))) = core.locations.get(&cell) else {
                return Err(invalid("axial endpoints must be original pipe cells"));
            };
            if *branch != binding.branch {
                return Err(invalid("axial link must remain inside its original branch"));
            }
            sections[side] = *section;
        }
        if sections[0].abs_diff(sections[1]) != 1 {
            return Err(invalid(
                "axial link must join consecutive original pipe cells",
            ));
        }
        if owned_pairs
            .insert(pair(link.a, link.b), binding.link)
            .is_some()
        {
            return Err(invalid("a whole axial face cannot have two link owners"));
        }
        let lengths = sections.map(|section| core.section_lengths[binding.branch][section]);
        let areas = sections.map(|section| core.section_areas[binding.branch][section]);
        let face = link.area_m2;
        if !near_positive(face, areas[0].min(areas[1])) {
            return Err(invalid(
                "axial link area must cover the whole shared pipe section",
            ));
        }
        let transport = sections.map(|section| {
            sources.sections[binding.branch][section]
                .momentum
                .evaluation
                .transport
        });
        let conductivity = transport.map(|state| state.thermal_conductivity_w_m_k());
        // Reuse the native common-area kernel with shape-equivalent distances:
        // d*_i/A_face = L_i/(2 A_i). Actual source conductivity stays untouched.
        // The face-area gate above prevents a partial aperture from cancelling
        // out of the one-dimensional resistance without a spreading model.
        let distances = std::array::from_fn(|i| (0.5 * lengths[i]) * (face / areas[i]));
        let conductance =
            fs_thermal::enthalpy::series_conductance_w_k(face, distances, conductivity).map_err(
                |source| GraphConductionError::Thermal {
                    link: binding.link,
                    source,
                },
            )?;
        let half_resistances = std::array::from_fn(|i| (distances[i] / conductivity[i]) / face);
        let temperatures = [thermal.temperature_k(link.a), thermal.temperature_k(link.b)];
        let mut t = [0.0; 2];
        for (side, temperature) in temperatures.into_iter().enumerate() {
            t[side] = temperature.map_err(|source| GraphConductionError::Thermal {
                link: binding.link,
                source,
            })?;
        }
        let difference = t[0] - t[1];
        let heat = conductance * difference;
        if !half_resistances.into_iter().all(positive)
            || !difference.is_finite()
            || !heat.is_finite()
            || (difference != 0.0 && heat == 0.0)
        {
            return Err(invalid("axial half resistance or heat is unrepresentable"));
        }
        conductances_w_k.insert(link.id, conductance);
        samples.push(GraphConductionSample {
            link: link.id,
            branch: binding.branch,
            cells: [link.a, link.b],
            sections,
            temperatures_k: t,
            lengths_m: lengths,
            areas_m2: areas,
            link_area_m2: face,
            transport,
            half_resistances_k_w: half_resistances,
            conductance_w_k: conductance,
            heat_a_to_b_w: heat,
        });
    }
    // One linear scan, not a link-by-link quadratic search, prevents a second
    // native contact from silently doubling a claimed whole axial face.
    for link in thermal.links().values() {
        poll(cx).map_err(source_error)?;
        if owned_pairs
            .get(&pair(link.a, link.b))
            .is_some_and(|&owner| owner != link.id)
        {
            return Err(GraphConductionError::Link {
                link: link.id,
                reason: "another native link duplicates a whole axial face",
            });
        }
    }
    poll(cx).map_err(source_error)?;
    Ok(ResolvedGraphConduction {
        conductances_w_k,
        samples,
        sections: sources.sections.clone(),
    })
}

fn pair(a: CellId, b: CellId) -> (CellId, CellId) {
    if a < b { (a, b) } else { (b, a) }
}
