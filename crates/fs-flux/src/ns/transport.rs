//! Directional transport on the same sign-split faces as momentum convection.

use super::{FluxSystem, convection_edge_quadrature};

/// Integrated volume rates per unit depth on one original mesh edge.
///
/// The global normal points outward from `edge.tris.0`, toward `edge.tris.1`
/// on an interior edge. Both directional rates are nonnegative; a reversing
/// affine normal trace can carry material BOTH ways even when its net is zero.
/// Multiplying each rate by depth [m] and density [kg/m³] gives mass flow [kg/s].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EdgeVolumeTransport {
    along_normal_m2_s: f64,
    against_normal_m2_s: f64,
    net_m2_s: f64,
}

impl EdgeVolumeTransport {
    /// Integral of `max(u.n, 0)` [m²/s].
    #[must_use]
    pub const fn along_normal_m2_s(&self) -> f64 {
        self.along_normal_m2_s
    }

    /// Integral of `max(-u.n, 0)` [m²/s].
    #[must_use]
    pub const fn against_normal_m2_s(&self) -> f64 {
        self.against_normal_m2_s
    }

    /// Signed integral of `u.n` [m²/s], positive along the global normal.
    /// Retains the original mean moment instead of subtracting potentially
    /// large opposing rates. Their difference agrees up to quadrature roundoff.
    #[must_use]
    pub const fn net_m2_s(&self) -> f64 {
        self.net_m2_s
    }
}

/// Invalid input or unrepresentable directional face transport.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FluxTransportError {
    /// Supply only the original velocity DOFs, without appended pressures.
    VelocityLength { expected: usize, actual: usize },
    /// One original velocity moment is not finite.
    NonFiniteVelocity { dof: usize },
    /// The original mesh edge length is not finite and positive.
    InvalidEdgeLength { edge: usize },
    /// A nonzero rate underflows, or an integrated rate overflows.
    UnrepresentableRate { edge: usize },
}

impl core::fmt::Display for FluxTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "BDM1 face transport: {self:?}")
    }
}

impl std::error::Error for FluxTransportError {}

// Multiplying the smallest and largest factors first avoids intermediate
// overflow/underflow when their final positive product is representable.
fn dimensional_rate(normalized: f64, scale: f64, length: f64) -> Option<f64> {
    if normalized == 0.0 {
        return Some(0.0);
    }
    let mut factors = [normalized, scale, length];
    factors.sort_unstable_by(f64::total_cmp);
    let value = (factors[0] * factors[2]) * factors[1];
    (value.is_finite() && value > 0.0).then_some(value)
}

impl FluxSystem<'_> {
    /// Resolve both physical directions on every original BDM1 mesh edge.
    ///
    /// Uses the SAME affine normal trace and zero-split quadrature as the
    /// momentum convection operator, evaluated in its unit face parameter.
    /// No cell interpolation, density, timestep, upwind enthalpy, or continuity
    /// correction is introduced. The caller owns the carried material state
    /// and must check mass balance when using these rates on fixed inventories.
    ///
    /// Results follow `mesh.edges` order. Work is O(edges), with at most four
    /// quadrature sites per edge and one output allocation. This is transport
    /// through the fixed mesh; ALE callers must supply relative normal moments.
    ///
    /// # Errors
    /// Refuses a vector other than exactly `n_u` finite velocity moments,
    /// invalid edge lengths, and unrepresentable rates. Inputs are unchanged.
    pub fn edge_volume_transport(
        &self,
        velocity: &[f64],
    ) -> Result<Vec<EdgeVolumeTransport>, FluxTransportError> {
        if velocity.len() != self.n_u {
            return Err(FluxTransportError::VelocityLength {
                expected: self.n_u,
                actual: velocity.len(),
            });
        }
        if let Some(dof) = velocity.iter().position(|value| !value.is_finite()) {
            return Err(FluxTransportError::NonFiniteVelocity { dof });
        }
        let mut result = Vec::with_capacity(self.mesh.edges.len());
        for (index, edge) in self.mesh.edges.iter().enumerate() {
            if !edge.len.is_finite() || edge.len <= 0.0 {
                return Err(FluxTransportError::InvalidEdgeLength { edge: index });
            }
            let moments = [velocity[2 * index], velocity[2 * index + 1]];
            let scale = moments[0].abs().max(moments[1].abs());
            if scale == 0.0 {
                result.push(EdgeVolumeTransport::default());
                continue;
            }
            let mean = moments[0] / scale;
            let linear = moments[1] / scale;
            let mut rates = [0.0; 2];
            let mut corrections = [0.0; 2];
            for (point, weight, _) in convection_edge_quadrature([0.0, 0.0], [1.0, 0.0], moments) {
                let normal = mean + 12.0 * linear * (point[0] - 0.5);
                let direction = usize::from(normal < 0.0);
                let corrected = weight * normal.abs() - corrections[direction];
                let next = rates[direction] + corrected;
                corrections[direction] = (next - rates[direction]) - corrected;
                rates[direction] = next;
            }
            let refusal = || FluxTransportError::UnrepresentableRate { edge: index };
            let along_normal_m2_s =
                dimensional_rate(rates[0], scale, edge.len).ok_or_else(refusal)?;
            let against_normal_m2_s =
                dimensional_rate(rates[1], scale, edge.len).ok_or_else(refusal)?;
            let net_m2_s = moments[0] * edge.len;
            if !net_m2_s.is_finite() || (moments[0] != 0.0 && net_m2_s == 0.0) {
                return Err(refusal());
            }
            result.push(EdgeVolumeTransport {
                along_normal_m2_s,
                against_normal_m2_s,
                net_m2_s,
            });
        }
        Ok(result)
    }
}
