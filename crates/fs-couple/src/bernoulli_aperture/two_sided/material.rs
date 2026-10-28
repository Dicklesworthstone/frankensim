//! Original solid material cards on both sides of the original moving valve.
//!
//! This is a cold composition of `wall::material` and `TwoSidedAperture`, not
//! another wall or reed equation. Wall flow reacts at the same junction whose
//! pressure drives it; its storage and loss already belong to the network audit.

use super::{ApertureDuctSpec, FractionalTubePolicy, LoadedAperture, TwoSidedAperture,
    TwoSidedLimits, ValveSide, invalid, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::wall::WallState;
use crate::bernoulli_aperture::wall::material::{
    BoundMaterialWall, MaterialWallInput, bind_material_walls,
};
use fs_exec::CancelGate;

// Receipts only. The sole evolving wall state is inside the acoustic network.
pub(super) struct MaterialDuctWalls {
    input: MaterialWallInput,
    walls: Vec<BoundMaterialWall>,
}

impl TwoSidedAperture {
    /// Bind independent structural materials in the upstream tract and downstream
    /// bore before either duct is subdivided or any valve/wave state is advanced.
    ///
    /// Entries in `walls` and `policies` are ordered [upstream, downstream]. An
    /// absent wall entry leaves that side's existing loads unchanged. `policies`
    /// selects the same optional fractional transit image as
    /// [`Self::with_fractional_transits`]; `None` retains original integer travel.
    /// Empty/absent wall selections do not add a physical load or change arithmetic.
    ///
    /// All geometry, regional gas and wall source records are retained, including
    /// generalized-Maxwell spectra. The existing single-/independent-blade solver,
    /// pressure/flow coupling, full energy audit and atomic two-network publication
    /// remain the only runtime. Ordinary two-sided pressure performers accept the
    /// returned object unchanged; there is no downstream-only audio filter.
    ///
    /// This admits only the existing fixed-state, thin isotropic local membrane
    /// wall image. It does not add carved/orthotropic shells, moving geometry,
    /// static-pressure prestress, thermal transport or exterior shell radiation.
    /// Source applicability bands are not filters on nonlinear generated sound.
    ///
    /// # Errors
    /// Original material/topology/load refusals, a wall band at/above Nyquist,
    /// unresolved wall time scales, exceeded source/numerical limits, an advanced
    /// valve, inconsistent gas/load declarations or cancellation. A failed cold
    /// construction publishes no partially coupled system.
    pub fn with_material_walls(
        valve: LoadedAperture,
        sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2],
        limits: TwoSidedLimits,
        policies: Option<[FractionalTubePolicy; 2]>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::build_material_walls(valve, sources, walls, limits, policies, gate, None)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn build_material_walls(
        valve: LoadedAperture,
        mut sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2],
        limits: TwoSidedLimits,
        policies: Option<[FractionalTubePolicy; 2]>,
        gate: &CancelGate,
        supply_channel: Option<crate::bernoulli_aperture::slit_channel::FixedChannel>,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let dt = valve.time_step_s();
        if valve.accepted_steps() != 0 {
            return Err(invalid("material ducts require an unadvanced valve"));
        }
        let mut original_graphs = [None, None];
        let mut retained = [None, None];
        for (i, input) in walls.into_iter().enumerate() {
            poll(gate)?;
            let Some(input) = input else { continue; };
            // Bound duplicated source receipts before cloning any assignment.
            // The physical owner below validates the full material band/domain.
            if input.max_walls > 1024 || input.assignments.len() > input.max_walls
                || !dt.is_finite() || dt <= 0.0 || input.band_hz[1] >= 0.5 / dt {
                return Err(invalid("material ducts exceed their wall-count, clock or Nyquist limits"));
            }
            let (graph, bound) = bind_material_walls(sources[i].graph.clone(),
                input.assignments.clone(), input.band_hz, input.max_walls, gate)?;
            for wall in &bound {
                poll(gate)?;
                wall.validate_time_step(dt)?;
            }
            // Keep original physical inputs rather than reporting lowered shunt
            // coefficients as if the caller had authored them. Gas and loss
            // selections are never replaced; original node indices are unchanged.
            original_graphs[i] = Some(core::mem::replace(&mut sources[i].graph, graph));
            retained[i] = Some(MaterialDuctWalls { input, walls: bound });
        }
        let mut model = Self::build(valve, sources, limits, gate, policies, supply_channel)?;
        for (side, original) in model.sides.iter_mut().zip(original_graphs) {
            if let Some(original) = original { side.source.graph = original; }
        }
        model.material_ducts = retained;
        poll(gate)?;
        Ok(model)
    }

    /// Exact original solid assignments and applicability limits, if supplied.
    /// Raw acoustic load coefficients never acquire synthetic material receipts.
    #[must_use]
    pub fn wall_input(&self, side: ValveSide) -> Option<&MaterialWallInput> {
        self.material_ducts[side.index()].as_ref().map(|w| &w.input)
    }

    /// Original source-bound walls; no second copy of their evolving states.
    #[must_use]
    pub fn material_walls(&self, side: ValveSide) -> &[BoundMaterialWall] {
        self.material_ducts[side.index()].as_ref().map_or(&[], |w| w.walls.as_slice())
    }

    fn material_wall(&self, side: ValveSide, node: usize)
        -> Result<&BoundMaterialWall, AcousticRealizeError>
    {
        self.material_walls(side).iter().find(|wall| wall.source().node == node)
            .ok_or_else(|| invalid("wall observation requires a source-bound node on the selected duct side"))
    }

    /// Accepted displacement, velocity and base mechanical energy at a physical
    /// wall node. Hereditary energy is separate; both are already in system energy.
    /// Reading this value cannot step or reset a material history.
    pub fn wall_state(&self, side: ValveSide, node: usize)
        -> Result<WallState, AcousticRealizeError>
    {
        self.material_wall(side, node)?.observe_waveguide(&self.sides[side.index()].network)
    }

    /// Accepted hereditary spring energy [J] at the selected material node.
    /// Together with `wall_state` this gives complete wall storage; do not add it
    /// again to `stored_energy_j`, which already includes every acoustic load.
    pub fn wall_relaxation_energy_j(&self, side: ValveSide, node: usize)
        -> Result<f64, AcousticRealizeError>
    {
        self.material_wall(side, node)?.relaxation_stored_energy_j(&self.sides[side.index()].network)
    }
}
