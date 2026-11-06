//! Frozen-damage vibration of the original thermally evolved brittle solid.
//!
//! Current source coefficients and accepted MECHANICAL strain choose the
//! native fixed-damage stiffness. Thermal eigenstrain is held, so it has no
//! incremental strain derivative; total displacement strain must not choose
//! the opening/compression branch. The original consistent reference mass
//! survives temperature-dependent density changes exactly as in native motion.
//!
//! This is a local conservative pencil with frozen temperature and fracture
//! history, not an evolving damage tangent or a thermoelastic damping model.

use fs_exec::Cx;

use super::ThermalBrittleTetMotion;
pub use crate::brittle::solid::dynamics::modes::{FrozenBrittleModesControls, FrozenBrittlePencil};
use crate::brittle::solid::dynamics::modes::{assemble_frozen_pencil, require_open_released_faces};
use crate::brittle::solid::{BrittleTetError, checkpoint};

impl ThermalBrittleTetMotion<'_> {
    /// Assemble the native `(K, M)` pencil at the accepted thermal state.
    ///
    /// `fixed_dofs` holds incremental displacements, independently of the
    /// loading that produced the current motion. Stiffness comes from the
    /// current source laws and their ORIGINAL accepted mechanical histories;
    /// neither the source temperature nor damage advances. Reference inertia
    /// is reduced from the same cached consistent matrix used by this motion.
    ///
    /// A damaged material at its normal closure switch refuses because it
    /// has no unique two-sided frozen tangent. Released face pairs must be
    /// strictly open beyond their original gap tolerance; closed or ambiguous
    /// pairs require an explicit contact tangent before extracting their modes.
    /// No geometric stiffness, thermal feedback, growth of damage, modal loss
    /// or radiation is added. A moving snapshot need not be an equilibrium
    /// state. Refusal and successful extraction both leave the motion intact.
    pub fn frozen_vibration_pencil(
        &self,
        fixed_dofs: &[usize],
        controls: FrozenBrittleModesControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittlePencil, BrittleTetError> {
        checkpoint(cx)?;
        require_open_released_faces(
            self.solid.geometry.nodes_m().len(),
            &self.solid.released_faces,
            &self.solid.solution.displacement_m,
            cx,
        )?;
        assemble_frozen_pencil(
            self.solid.geometry,
            self.solid.budget,
            self.solid
                .materials
                .iter()
                .map(|source| source.band().law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.mass,
            fixed_dofs,
            controls,
            cx,
        )
    }
}
