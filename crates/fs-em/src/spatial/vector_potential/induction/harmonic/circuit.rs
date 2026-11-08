//! Original electrical networks driving closed three-dimensional AC windings.
//!
//! The shared native circuit descriptor and terminal closure retain the same
//! branch, source and transformer constraints as planar winding circuits. The
//! recovered 3D field supplies its own spatial copper heat and signed motion
//! work. Prescribed mechanical motion may return electrical power; it is never
//! replaced by a positive equivalent loss or an inferred network residual.

use super::winding::{
    SpatialHarmonicWindingError, SpatialHarmonicWindingFrequency, SpatialHarmonicWindingResponse,
};
use crate::planar::eddy::harmonic::circuit::{
    HarmonicCircuitField, HarmonicCircuitResponse, circuit_response,
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::c64::C64;
use fs_phs::circuit::{CircuitGraph, MutualInductance};

pub use crate::planar::eddy::harmonic::circuit::{
    HarmonicCircuitBranch, HarmonicCircuitControls, HarmonicCircuitError, HarmonicCircuitSource,
    HarmonicCircuitTransformer, HarmonicSourceKind,
};

/// Physical original network and its resolved three-dimensional field response.
pub type SpatialHarmonicCircuitResponse<'s, 'm, C = ResolvedConductorState> =
    HarmonicCircuitResponse<SpatialHarmonicWindingResponse<'s, 'm, C>>;

impl<C: SpatialConductorMaterial> HarmonicCircuitField
    for SpatialHarmonicWindingResponse<'_, '_, C>
{
    fn currents_a(&self) -> &[C64] {
        &self.currents_a
    }
    fn terminal_voltages_v(&self) -> &[C64] {
        &self.terminal_voltages_v
    }
    fn mean_dissipated_power_w(&self) -> f64 {
        self.mean_dissipated_power_w()
    }
    fn mean_mechanical_power_w(&self) -> f64 {
        self.field.balance().mean_mechanical_power_w
    }
    fn absorbed_power_va(&self) -> C64 {
        self.supplied_power_va
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicWindingFrequency<'_, 's, '_, 'm, C> {
    /// Close the original RLC/transformer/source graph on every field winding.
    /// `winding_branches` lists distinct oriented inductor placeholders in the
    /// original winding order. Source-bound winding resistance, resolved mutual
    /// flux and declared omitted end storage replace those numeric inductances.
    /// External mutual links may connect only remaining ordinary inductors.
    ///
    /// The shared native descriptor is factored once and its small terminal
    /// interface is solved once. Cached field cochains seed the actual drive,
    /// whose original sparse field equations and physical work are rechecked.
    /// Source indices, peak exp(-i omega t) phases, transformer power and KCL
    /// remain in the original graph convention. Numerical reference resistor
    /// heat is removed before reporting any physical loss.
    ///
    /// Mean field heat and signed mechanical work are separate outputs. Their
    /// sum, including possible electrical generation, closes the real circuit
    /// supply; reactive storage and motional power close its imaginary part.
    /// This remains a linear harmonic response at prescribed motion and fixed
    /// material state, without advancing temperature, geometry or shaft speed.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_response(
        &self,
        graph: &CircuitGraph,
        winding_branches: &[usize],
        external_mutual: &[MutualInductance],
        source_phasors: &[C64],
        controls: HarmonicCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicCircuitResponse<'s, 'm, C>, SpatialHarmonicWindingError> {
        circuit_response(
            self.omega_rad_s(),
            self.impedance_ohm(),
            None,
            graph,
            winding_branches,
            external_mutual,
            source_phasors,
            controls,
            cancelled,
            |currents, cancel| self.current_response(currents, cancel),
        )
    }
}
