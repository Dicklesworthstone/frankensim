//! Fixed, cell-resolved substrate chemistry on the original planar boundary.
//!
//! Each patch owns one original source bundle, shared by its cells. The cell
//! index is immutable and shared by thermal/native trial clones. No averaging
//! creates a synthetic material claim at a chemical boundary.

use super::*;

/// One explicit stationary substrate patch, resolved at native wall-cell scale.
///
/// Surface A and its finish/history may differ from the base substrate. Surface
/// B, liquid medium and named gas environment must remain the same. Each patch
/// selects its own original equilibrium angle or reversible adhesion-work claim;
/// the thermal owner supplies the actual local liquid temperature and sigma.
/// Cells are full wall footprints, not subcell area fractions. Refine the grid
/// to resolve patch edges; no Cassie/Wenzel average or exact pinning is inferred.
#[derive(Clone, Debug)]
pub struct EquilibriumWettingPatch<'m> {
    pub cells: Vec<[usize; 3]>,
    pub interface: &'m InterfaceSystemCard,
    pub origin: EquilibriumAngleOrigin,
    pub selection: MaterialPropertySelection,
}

#[derive(Clone)]
pub(super) struct WettingPattern<'m> {
    sources: Vec<WettingPhysics<'m>>,
    cells: Arc<BTreeMap<[usize; 3], usize>>,
}

impl<'m> WettingBinding<'m> {
    pub(super) fn physics_at(&self, position: [usize; 3]) -> &WettingPhysics<'m> {
        self.pattern
            .as_ref()
            .and_then(|pattern| pattern.cells.get(&position).map(|&i| &pattern.sources[i]))
            .unwrap_or(&self.physics)
    }

    pub(super) fn physics_mut(&mut self) -> impl Iterator<Item = &mut WettingPhysics<'m>> {
        std::iter::once(&mut self.physics).chain(
            self.pattern
                .iter_mut()
                .flat_map(|pattern| pattern.sources.iter_mut()),
        )
    }

    pub(in crate::material::free_surface) fn with_equilibrium_patches(
        mut self,
        patches: &[EquilibriumWettingPatch<'m>],
        domain_cells: [usize; 3],
        scaling: &MaterialFreeSurfaceScaling<'_>,
        temperature_axis: &str,
        evolving_surface_tension: bool,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        if self.pattern.is_some()
            || self.physics.equilibrium_hysteresis.is_some()
            || self
                .physics
                .equilibrium
                .as_ref()
                .is_none_or(|equilibrium| equilibrium.dynamic)
        {
            return Err(Error::Invalid(
                "substrate patches require an unpatterned equilibrium binding before hysteresis/rate-law selection",
            ));
        }
        let maximum_cells = self.physics.controls.maximum_wall_cells;
        if patches.is_empty() || patches.len() > maximum_cells {
            return Err(Error::Invalid(
                "substrate patches require a nonempty bounded source set",
            ));
        }
        let (axis, upper) = match self.wall {
            PlanarWall3::XMin => (0, false),
            PlanarWall3::XMax => (0, true),
            PlanarWall3::YMin => (1, false),
            PlanarWall3::YMax => (1, true),
            PlanarWall3::ZMin => (2, false),
            PlanarWall3::ZMax => (2, true),
        };
        let edge = if upper { domain_cells[axis] - 1 } else { 0 };
        let mut cells = BTreeMap::new();
        let mut sources = Vec::with_capacity(patches.len());
        for patch in patches {
            checkpoint(gate)?;
            if patch.cells.is_empty() || patch.cells.len() > maximum_cells - cells.len() {
                return Err(Error::Invalid(
                    "substrate patch cells exceed the original wall-cell budget",
                ));
            }
            // A new finish or coating is local to the substrate. It cannot
            // silently change this liquid or the shared surrounding atmosphere.
            if patch.interface.surface_b() != self.physics.card.surface_b()
                || patch.interface.medium() != self.physics.card.medium()
                || patch.interface.environment() != self.physics.card.environment()
            {
                return Err(Error::Invalid(
                    "substrate patches must retain the original liquid and named gas environment",
                ));
            }
            let index = sources.len();
            for &position in &patch.cells {
                checkpoint(gate)?;
                if position[axis] != edge || (0..3).any(|d| position[d] >= domain_cells[d]) {
                    return Err(Error::Invalid(
                        "substrate patch cell is outside the assigned planar wall",
                    ));
                }
                if cells.insert(position, index).is_some() {
                    return Err(Error::Invalid(
                        "substrate patches overlap or repeat a wall cell",
                    ));
                }
            }
            sources.push(WettingPhysics::try_equilibrium(
                scaling,
                patch.interface,
                patch.origin,
                patch.selection.clone(),
                temperature_axis,
                evolving_surface_tension,
                maximum_cells,
                gate,
            )?);
        }
        checkpoint(gate)?;
        self.pattern = Some(WettingPattern {
            sources,
            cells: Arc::new(cells),
        });
        // The thermal caller resolves all currently active cells at their actual
        // temperatures before installing this candidate. Unlisted cells keep
        // the explicit base substrate; inactive declared cells retain their
        // patch source if a native VOF tile activates later.
        Ok(self)
    }
}
