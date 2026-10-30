//! Three-dimensional original-species transport on the shared molecular kernel.

use super::*;
use crate::multicomponent::{MAX_SPATIAL_COMPONENTS, SpatialCell, SpatialPrimitive};

/// Prepared molecular transport in an orthonormal three-dimensional face frame.
///
/// This owner reuses the exact axial/planar Maxwell--Stefan solver, caloric
/// references, source identity and prepared face workspace. Both transverse
/// velocity derivatives enter the three-dimensional Newtonian stress. Species
/// enthalpy, Fourier heat and all three components of stress work accompany the
/// same conservative species/momentum/energy flux; no extra heating source is
/// needed. Soret/Dufour effects and nonideal transport remain separate closures.
pub struct SpatialTransport<L: TransportLaw> {
    prepared: PlanarTransport<L>,
}

impl<L: TransportLaw> SpatialTransport<L> {
    /// Bind the exact original caloric species and one transport law. Numerical
    /// scratch is prepared once and reused by every successful face evaluation.
    pub fn new(
        gas: MixtureGas,
        law: L,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        Ok(Self {
            prepared: PlanarTransport::new(gas, law, cancelled)?,
        })
    }

    /// Original caloric source, species order and admitted temperature interval.
    pub const fn gas(&self) -> &MixtureGas {
        self.prepared.gas()
    }

    /// Original transport law and its source and validity declarations.
    pub const fn law(&self) -> &L {
        self.prepared.law()
    }

    /// Prepared numerical heap capacity, including the original Maxwell--Stefan
    /// workspace. Immutable gas storage, bounded stack scratch and the supplied
    /// law's own storage/evaluation allocations retain their original contracts.
    pub fn workspace_bytes(&self) -> usize {
        self.prepared.workspace_bytes()
    }

    /// Conservative normal flux `[j_i..., -tau_nn, -tau_t1n, -tau_t2n, E_flux]`.
    ///
    /// `distance` is the positive center-to-center separation. The gradient rows
    /// are derivatives along the two face tangents, and columns are velocities
    /// `[u_n, u_t1, u_t2]` in that same orthonormal frame. Normal derivatives use
    /// the two supplied traces. Output is scratch until success; inputs remain
    /// unchanged. The caller owns spatial divergence and transactional updates.
    #[allow(clippy::too_many_arguments)]
    pub fn face(
        &mut self,
        left: SpatialCell<'_>,
        right: SpatialCell<'_>,
        distance: f64,
        tangential_velocity_gradient: [[f64; 3]; 2],
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let left = conservative_trace(&self.prepared.gas, left, &mut cancelled)?;
        let right = conservative_trace(&self.prepared.gas, right, &mut cancelled)?;
        self.interior(
            left,
            right,
            distance,
            tangential_velocity_gradient,
            controls,
            out,
            &mut cancelled,
        )
    }

    /// Use already reconstructed primitives while retaining the original gas's
    /// complete three-momentum caloric admission, without another T inversion.
    #[allow(clippy::too_many_arguments)]
    pub fn face_from_primitives(
        &mut self,
        left: SpatialPrimitive<'_>,
        right: SpatialPrimitive<'_>,
        distance: f64,
        tangential_velocity_gradient: [[f64; 3]; 2],
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let left = primitive_trace(&self.prepared.gas, left, &mut cancelled)?;
        let right = primitive_trace(&self.prepared.gas, right, &mut cancelled)?;
        self.interior(
            left,
            right,
            distance,
            tangential_velocity_gradient,
            controls,
            out,
            &mut cancelled,
        )
    }

    /// Stationary adiabatic impermeable slip wall. The transverse derivatives
    /// are `[d_t1 u_t1, d_t2 u_t2]` in the same frame as `d_n u_n`. Normal
    /// traction retains both transverse dilatations; every species, both shear
    /// tractions and energy flux are exactly zero. A spatial owner supplies the
    /// appropriate lower/upper sign for the normal wall velocity derivative.
    pub fn wall(
        &mut self,
        cell: SpatialCell<'_>,
        normal_velocity_gradient: f64,
        transverse_tangential_gradients: [f64; 2],
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let trace = conservative_trace(&self.prepared.gas, cell, &mut cancelled)?;
        self.wall_trace(
            trace,
            normal_velocity_gradient,
            transverse_tangential_gradients,
            out,
            &mut cancelled,
        )
    }

    /// The same slip-wall closure from source-admitted reconstructed primitives.
    pub fn wall_from_primitive(
        &mut self,
        cell: SpatialPrimitive<'_>,
        normal_velocity_gradient: f64,
        transverse_tangential_gradients: [f64; 2],
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let trace = primitive_trace(&self.prepared.gas, cell, &mut cancelled)?;
        self.wall_trace(
            trace,
            normal_velocity_gradient,
            transverse_tangential_gradients,
            out,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn interior(
        &mut self,
        left: VectorTrace<'_, 2>,
        right: VectorTrace<'_, 2>,
        distance: f64,
        tangential_velocity_gradient: [[f64; 3]; 2],
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let prepared = &mut self.prepared;
        interior_components(
            &prepared.gas,
            &mut prepared.law,
            &prepared.diffusion,
            &mut prepared.workspace,
            left,
            right,
            distance,
            distance,
            tangential_velocity_gradient,
            prepared.mass_ratio,
            &mut prepared.binary,
            out,
            controls,
            cancelled,
        )
    }

    fn wall_trace(
        &mut self,
        trace: VectorTrace<'_, 2>,
        normal_velocity_gradient: f64,
        transverse_tangential_gradients: [f64; 2],
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let prepared = &mut self.prepared;
        wall_components(
            &prepared.gas,
            &mut prepared.law,
            trace,
            normal_velocity_gradient,
            transverse_tangential_gradients,
            &mut prepared.binary,
            out,
            cancelled,
        )
    }
}

fn conservative_trace<'a>(
    gas: &MixtureGas,
    cell: SpatialCell<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<VectorTrace<'a, 2>, DissipationError> {
    let primitive = gas.spatial_primitive(cell, cancelled)?;
    let cv_density = heat_capacity_density(gas, cell.partial_densities, primitive)?;
    Ok(VectorTrace {
        densities: cell.partial_densities,
        primitive,
        transverse_velocity: cell.transverse_momentum.map(|p| p / primitive.density),
        cv_density,
    })
}

fn primitive_trace<'a>(
    gas: &MixtureGas,
    cell: SpatialPrimitive<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<VectorTrace<'a, 2>, DissipationError> {
    poll(cancelled)?;
    let mut q = [0.0; MAX_SPATIAL_COMPONENTS];
    let primitive = gas.spatial_conserved_into(cell, &mut q[..gas.species_count() + 4])?;
    let cv_density = heat_capacity_density(gas, cell.partial_densities, primitive)?;
    Ok(VectorTrace {
        densities: cell.partial_densities,
        primitive,
        transverse_velocity: cell.transverse_velocity,
        cv_density,
    })
}
