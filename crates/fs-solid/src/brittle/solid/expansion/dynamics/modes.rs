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
use crate::brittle::solid::dynamics::fracture::contact::friction::BrittleFaceFriction;
use crate::brittle::solid::dynamics::modes::contact::assemble_contact_modes;
pub use crate::brittle::solid::dynamics::modes::contact::friction::FrozenBrittleStickingControls;
use crate::brittle::solid::dynamics::modes::contact::friction::assemble_sticking_contact_modes;
pub use crate::brittle::solid::dynamics::modes::contact::{
    FrozenBrittleContactControls, FrozenBrittleContactError, FrozenBrittleContactModes,
};
pub use crate::brittle::solid::dynamics::modes::{FrozenBrittleModesControls, FrozenBrittlePencil};
use crate::brittle::solid::dynamics::modes::{assemble_frozen_pencil, require_open_released_faces};
use crate::brittle::solid::{BrittleTetError, checkpoint};

impl ThermalBrittleTetMotion<'_> {
    /// Observe a held static-friction branch using the native directional ellipse.
    ///
    /// Each face has three preload vectors in its declared normal/tangent frame.
    /// Loaded, resting forces strictly inside a positive friction ellipse
    /// impose exact normal and both tangential constraints. A zero ellipse
    /// preserves frictionless tangent freedom; a singular one-axis ellipse or
    /// force at the slip boundary refuses. Coefficients are explicit model
    /// inputs, not new source evidence or tangential penalty stiffnesses.
    ///
    /// Thermal coefficients, mechanical histories, original reference mass
    /// and time stay with this motion. Contact membership and temperature are
    /// held. Equilibrium and finite-amplitude sticking are caller obligations;
    /// no frictional heat, damping or evolving slip is supplied by these modes.
    pub fn frozen_sticking_contact_modes(
        &self,
        friction: &[BrittleFaceFriction],
        preload_n: &[[[f64; 3]; 3]],
        fixed_dofs: &[usize],
        controls: FrozenBrittleStickingControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittleContactModes, FrozenBrittleContactError> {
        assemble_sticking_contact_modes(
            self.solid.geometry,
            self.solid.budget,
            self.solid
                .materials
                .iter()
                .map(|source| source.band().law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.velocity_mps,
            &self.solid.released_faces,
            &self.mass,
            friction,
            preload_n,
            fixed_dofs,
            controls,
            cx,
        )
    }

    /// Observe conditional modes with explicitly preloaded normal crack contact.
    ///
    /// The three entries of each `normal_preload_n` row are normal forces at
    /// the retained face's three paired vertices. Loaded, matched, normally
    /// resting pairs impose exact incremental normal constraints; tangential
    /// motion remains free. Touching without strictly positive preload is a
    /// unilateral branch switch and refuses. The caller holds the active set
    /// and establishes equilibrium if equilibrium natural frequencies are
    /// required; an interval-mean force receipt alone does not establish it.
    ///
    /// This shares the mechanical contact reduction and original reference
    /// inertia. Current thermal coefficients and accepted mechanical strain
    /// determine bulk stiffness. Temperature, eigenstrain, fracture history,
    /// contact geometry and normal direction remain fixed. No contact spring,
    /// contact inertia, heat, damping or finite-amplitude guarantee is added.
    pub fn frozen_contact_modes(
        &self,
        normal_preload_n: &[[f64; 3]],
        fixed_dofs: &[usize],
        controls: FrozenBrittleContactControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittleContactModes, FrozenBrittleContactError> {
        assemble_contact_modes(
            self.solid.geometry,
            self.solid.budget,
            self.solid
                .materials
                .iter()
                .map(|source| source.band().law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.velocity_mps,
            &self.solid.released_faces,
            &self.mass,
            normal_preload_n,
            fixed_dofs,
            controls,
            cx,
        )
    }

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
