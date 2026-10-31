//! Causal solid-winding skin/proximity currents on the native nonlinear circuit.
//!
//! The original A_z history supplies magnetic diffusion memory. Regional
//! conductivity supplies the SAME normalized source, uniform electric drive,
//! DC resistance and total-field heat. A refreshed source F(T) contributes
//! terminal voltage F_new^T delta(A)/dt; delta(F^T A)/dt would invent an EMF
//! when a stationary biased field is merely observed with different weights.

use super::*;
pub use crate::planar::eddy::solid::SolidWindingRegion;
use crate::planar::eddy::{EddyError, solid::SolidRegionLayout};
use std::collections::BTreeMap;

pub(super) struct PreparedSolid {
    pub profiles: Rc<Vec<Vec<f64>>>,
    pub layout: SolidRegionLayout,
}

fn induction(error: EddyError) -> NonlinearWindingError {
    NonlinearEddyError::from(error).into()
}

/// Source geometry is already bounded by the accepted native field. Original
/// material/claim continuity is independently rechecked by that field on every
/// solve; this preparation cannot publish a candidate on its own.
pub(super) fn prepare(
    field: &PlanarNonlinearEddy<'_>,
    sources: &[EddyConductor],
    profiles: &[Vec<f64>],
    bindings: &[SolidWindingRegion],
    admit_original: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedSolid> {
    let mesh = field.problem().mesh();
    if sources.len() != field.conductors().len()
        || profiles.iter().any(|p| p.len() != mesh.faces().len())
    {
        return Err(NonlinearWindingError::Invalid(
            "complete original solid source geometry required",
        ));
    }
    let mut order: Vec<_> = (0..sources.len()).collect();
    order.sort_by_key(|&i| sources[i].triangle);
    let mut weights = BTreeMap::<u32, f64>::new();
    let mut cells = Vec::with_capacity(sources.len());
    for (slot, &index) in order.iter().enumerate() {
        poll(cancelled)?;
        let source = &sources[index];
        let original = &field.conductors()[slot];
        if source.triangle != original.triangle
            || source.region != original.region
            || !source
                .conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
        {
            return Err(NonlinearWindingError::Invalid(
                "solid regions retain original triangles and static conductivity",
            ));
        }
        let sigma = source.conductor.resistivity_ohm_m().recip();
        let volume = mesh
            .face_measure(source.triangle)
            .ok_or(NonlinearWindingError::Invalid(
                "missing original solid volume",
            ))?;
        let weight = finite(sigma * volume)?;
        if !positive(sigma) || !positive(weight) {
            return Err(NonlinearWindingError::Invalid(
                "positive representable solid conductivity volume required",
            ));
        }
        let total = weights.entry(source.region).or_default();
        *total = finite(*total + weight)?;
        cells.push((source.triangle, source.region, sigma));
    }
    let regions: Vec<_> = weights.keys().copied().collect();
    let weights: Vec<_> = weights.values().copied().collect();
    let layout = SolidRegionLayout::new(
        &regions,
        &weights,
        profiles.len(),
        field.problem().thickness_m(),
        bindings.to_vec(),
        cancelled,
    )
    .map_err(induction)?;
    let mut next = profiles.to_vec();
    for (triangle, region, sigma) in cells {
        poll(cancelled)?;
        let group = regions
            .binary_search(&region)
            .map_err(|_| NonlinearWindingError::Invalid("missing prepared solid region"))?;
        if let Some(index) = layout.region_bindings[group] {
            let binding = layout.bindings[index];
            let expected = finite(sigma * layout.electric_per_current_v_m_a[group])?;
            let original = profiles[binding.winding][triangle];
            if expected == 0.0
                || original == 0.0
                || expected.signum() != original.signum()
                || profiles
                    .iter()
                    .enumerate()
                    .any(|(w, p)| w != binding.winding && p[triangle] != 0.0)
                || admit_original
                    && (original - expected).abs()
                        > 256.0 * f64::EPSILON * original.abs().max(expected.abs())
            {
                return Err(NonlinearWindingError::Invalid(
                    "solid profile must match its original conductivity-normalized region and polarity",
                ));
            }
            next[binding.winding][triangle] = expected;
        } else if profiles.iter().any(|p| p[triangle] != 0.0) {
            return Err(NonlinearWindingError::Invalid(
                "active resolved conductor requires an explicit solid current binding",
            ));
        }
    }
    Ok(PreparedSolid {
        profiles: Rc::new(next),
        layout,
    })
}

impl<'m> NonlinearWindingCircuit<'m> {
    /// Bind source-normalized solid sections, with default current-input native
    /// ports. The supplied initial profiles must already represent their source
    /// conductivity; initial field, native storage and current are never reset.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_solid_regions(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        sources: Vec<WindingConductor>,
        solid_regions: Vec<SolidWindingRegion>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let count = windings.profiles_per_m2.len();
        if count == 0 || count > config.max_windings || config.max_windings > MAX_WINDINGS {
            return Err(NonlinearWindingError::Invalid(
                "invalid winding count or capacity",
            ));
        }
        Self::new_complete(
            field,
            graph,
            windings,
            vec![NonlinearWindingPort::CurrentInput; count],
            sources,
            solid_regions,
            initial,
            config,
            cancelled,
        )
    }

    /// Solid and stranded sections may share a winding, with explicit native
    /// current/voltage-input causality. Every solid return section is explicit;
    /// unbound resolved regions remain floating. Source cells are disjoint.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_solid_regions_and_ports(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        sources: Vec<WindingConductor>,
        solid_regions: Vec<SolidWindingRegion>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_complete(
            field,
            graph,
            windings,
            terminal_ports,
            sources,
            solid_regions,
            initial,
            config,
            cancelled,
        )
    }

    /// Original connected-region winding assignment and signed return geometry.
    pub fn solid_regions(&self) -> &[SolidWindingRegion] {
        &self.solid_regions
    }
}
