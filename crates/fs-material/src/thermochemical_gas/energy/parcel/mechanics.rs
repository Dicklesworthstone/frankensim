//! Work-conjugate pressure loads and matrix-free mechanical tangents.
//!
//! The geometry owner supplies V(q), its gradient and Hessian actions in its
//! existing coordinate chart. This module pulls the SAME material pressure and
//! acoustic bulk modulus back to that chart; it owns no geometry or integrator.
//! For a constant exterior pressure pe, F=(p-pe) grad(V), and
//! DF[dq]=(p-pe) Hess(V)[dq]-(rho*c^2/V) grad(V) (grad(V).dq).
//! Both terms matter: the geometric term does not vanish in curved chambers.
//! Changing gas data changes force, stiffness, transport and sound together.

use super::{ClosedGasParcel, GasEntropyError, ResolvedThermochemicalGasState, positive};

/// Default work/allocation bound for a small local mechanical response.
/// Mesh consumers may supply their own admitted coordinate limit and cancellation.
pub const MAX_VOLUME_COORDINATES: usize = 256;

/// Local adiabatic pressure response from the original acoustic material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasPressureResponse {
    /// Absolute gas pressure [Pa].
    pub pressure_pa: f64,
    /// Adiabatic bulk modulus rho*c^2 [Pa].
    pub bulk_modulus_pa: f64,
    /// dp/dV at fixed mass and entropy [Pa/m^3].
    pub pressure_volume_derivative_pa_m3: f64,
}

/// Immutable borrowed material/geometry epoch. The parcel cannot be heated or
/// moved while this view is used, preventing stale thermodynamic coefficients.
/// Geometry derivatives and coordinate/frame correctness remain caller authority.
#[derive(Debug)]
pub struct GasMechanicalLoad<'a, S = ResolvedThermochemicalGasState> {
    state: &'a S,
    gradient: &'a [f64],
    response: GasPressureResponse,
    net_pressure_pa: f64,
}

fn invalid(message: &'static str) -> GasEntropyError { GasEntropyError::Invalid(message) }
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasEntropyError> {
    if cancelled() { Err(GasEntropyError::Cancelled) } else { Ok(()) }
}
fn finite_vector(
    values: &[f64], length: usize, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasEntropyError> {
    if values.len() != length {
        return Err(invalid("volume response vector has wrong length or nonfinite entries"));
    }
    for block in values.chunks(MAX_VOLUME_COORDINATES) {
        checkpoint(cancelled)?;
        if block.iter().any(|v| !v.is_finite()) {
            return Err(invalid("volume response vector has wrong length or nonfinite entries"));
        }
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64, GasEntropyError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("mechanical pressure response overflow")) }
}
fn dot(
    a: &[f64], b: &[f64], cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasEntropyError> {
    // Fixed-order compensated accumulation, independent of worker scheduling.
    let (mut sum, mut correction) = (0.0, 0.0);
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        if i % MAX_VOLUME_COORDINATES == 0 { checkpoint(cancelled)?; }
        let adjusted = finite(x * y)? - correction;
        let next = finite(sum + adjusted)?;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    finite(sum)
}

impl ClosedGasParcel {
    /// Reuse frozen-composition acoustic c for mechanical compressibility.
    /// No nominal air constant or fixed-gamma approximation is introduced.
    pub fn pressure_response(&self) -> Result<GasPressureResponse, GasEntropyError> {
        let gas = self.state().state();
        let bulk = gas.density * gas.sound_speed * gas.sound_speed;
        let derivative = -bulk / self.volume_m3();
        if !positive(bulk) || !derivative.is_finite() || derivative >= 0.0 {
            return Err(invalid("adiabatic pressure-volume response is unrepresentable"));
        }
        Ok(GasPressureResponse { pressure_pa: gas.pressure, bulk_modulus_pa: bulk,
            pressure_volume_derivative_pa_m3: derivative })
    }

    /// Pull pressure back through a caller-owned volume gradient, at the EXACT
    /// accepted parcel volume. Exterior pressure is explicit, absolute and
    /// constant during linearization; zero means a vacuum, never a default.
    ///
    /// `gradient[i]` is dV/dq_i. q remains the caller's existing chart: e.g.
    /// metres -> m^2 and newtons, radians -> m^3/rad and torque. Rate/direction
    /// and Hessian actions must use that same chart, units, frame and ordering.
    /// No coordinate registry, derivative estimator or volume law is invented.
    pub fn mechanical_load<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
    ) -> Result<GasMechanicalLoad<'a>, GasEntropyError> {
        GasMechanicalLoad::new(
            self.state(),
            self.volume_m3(),
            volume_m3,
            gradient,
            exterior_pressure_pa,
            self.pressure_response()?,
        )
    }

    /// Borrow the same pressure response for a caller-budgeted spatial chart.
    ///
    /// The admitted mesh owns maximum_coordinates; it is not a new constitutive
    /// limit. This view borrows all derivatives without allocating. Admission
    /// and cancellable actions poll at most every 256 coordinates. Use the
    /// cancellable actions when applying the load to a large chart.
    pub fn mechanical_load_with_limit<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'a>, GasEntropyError> {
        GasMechanicalLoad::new_with_limit(
            self.state(),
            self.volume_m3(),
            volume_m3,
            gradient,
            exterior_pressure_pa,
            self.pressure_response()?,
            maximum_coordinates,
            cancelled,
        )
    }
}

impl<'a, S> GasMechanicalLoad<'a, S> {
    // Shared pullback for ideal and nonideal material states. The borrowed state
    // fixes the physical epoch; each EOS supplies its own pressure/bulk response.
    pub(in crate::thermochemical_gas) fn new(
        state: &'a S,
        current_volume_m3: f64,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        response: GasPressureResponse,
    ) -> Result<Self, GasEntropyError> {
        Self::new_with_limit(
            state, current_volume_m3, volume_m3, gradient, exterior_pressure_pa,
            response, MAX_VOLUME_COORDINATES, || false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::thermochemical_gas) fn new_with_limit(
        state: &'a S,
        current_volume_m3: f64,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        response: GasPressureResponse,
        maximum_coordinates: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEntropyError> {
        checkpoint(&mut cancelled)?;
        if volume_m3 != current_volume_m3 || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0 || gradient.is_empty()
            || maximum_coordinates == 0 || gradient.len() > maximum_coordinates
        {
            return Err(invalid("mechanical load needs current volume, nonnegative exterior pressure and coordinates within the declared limit"));
        }
        finite_vector(gradient, gradient.len(), &mut cancelled)?;
        let net_pressure_pa = finite(response.pressure_pa - exterior_pressure_pa)?;
        checkpoint(&mut cancelled)?;
        Ok(Self {
            state,
            gradient,
            response,
            net_pressure_pa,
        })
    }

    /// Original source-bearing medium, also consumed by thermal/acoustic paths.
    #[must_use]
    pub const fn state(&self) -> &S {
        self.state
    }
    /// Absolute gas pressure and its adiabatic derivative, not net pressure.
    #[must_use]
    pub const fn pressure_response(&self) -> GasPressureResponse { self.response }
    /// Gas minus the caller's fixed exterior pressure [Pa].
    #[must_use]
    pub const fn net_pressure_pa(&self) -> f64 { self.net_pressure_pa }

    /// Work-conjugate generalized force, (p-pe)*grad(V), in the caller's chart.
    pub fn generalized_forces(&self) -> Result<Vec<f64>, GasEntropyError> {
        self.generalized_forces_cancellable(|| false)
    }

    /// The same force with bounded cancellation latency and atomic publication.
    pub fn generalized_forces_cancellable(
        &self, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, GasEntropyError> {
        checkpoint(&mut cancelled)?;
        let mut force = Vec::with_capacity(self.gradient.len());
        for (i, &g) in self.gradient.iter().enumerate() {
            if i % MAX_VOLUME_COORDINATES == 0 { checkpoint(&mut cancelled)?; }
            force.push(finite(self.net_pressure_pa * g)?);
        }
        checkpoint(&mut cancelled)?;
        Ok(force)
    }

    /// Actual force Jacobian action at fixed mass/entropy. `volume_hessian_action`
    /// is Hess(V)[direction], NOT a dense matrix or a directional second scalar.
    /// This O(n) operation includes the geometric pressure term and the shared
    /// acoustic compressibility term. It does not assert positive stiffness,
    /// stability or symmetry of unverified caller-supplied geometry derivatives.
    pub fn tangent_action(
        &self, direction: &[f64], volume_hessian_action: &[f64],
    ) -> Result<Vec<f64>, GasEntropyError> {
        self.tangent_action_cancellable(direction, volume_hessian_action, || false)
    }

    /// Matrix-free Jacobian action on a caller-budgeted chart, with the same
    /// source acoustic compressibility and geometric pressure terms.
    pub fn tangent_action_cancellable(
        &self,
        direction: &[f64],
        volume_hessian_action: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, GasEntropyError> {
        checkpoint(&mut cancelled)?;
        finite_vector(direction, self.gradient.len(), &mut cancelled)?;
        finite_vector(volume_hessian_action, self.gradient.len(), &mut cancelled)?;
        let dp = finite(self.response.pressure_volume_derivative_pa_m3
            * dot(self.gradient, direction, &mut cancelled)?)?;
        let mut action = Vec::with_capacity(self.gradient.len());
        for (i, (&g, &h)) in self.gradient.iter().zip(volume_hessian_action).enumerate() {
            if i % MAX_VOLUME_COORDINATES == 0 { checkpoint(&mut cancelled)?; }
            action.push(finite(finite(self.net_pressure_pa * h)? + finite(dp * g)?)?);
        }
        checkpoint(&mut cancelled)?;
        Ok(action)
    }

    /// Net mechanical power OUT of gas plus the fixed-pressure exterior [W].
    /// F.q_dot equals (p-pe)*V_dot; energy integration remains with the caller.
    pub fn mechanical_power_w(&self, coordinate_rates: &[f64]) -> Result<f64, GasEntropyError> {
        self.mechanical_power_w_cancellable(coordinate_rates, || false)
    }

    /// Work-conjugate power on a spatial chart with bounded cancellation checks.
    pub fn mechanical_power_w_cancellable(
        &self, coordinate_rates: &[f64], mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, GasEntropyError> {
        checkpoint(&mut cancelled)?;
        finite_vector(coordinate_rates, self.gradient.len(), &mut cancelled)?;
        let power = finite(self.net_pressure_pa
            * dot(self.gradient, coordinate_rates, &mut cancelled)?)?;
        checkpoint(&mut cancelled)?;
        Ok(power)
    }
}

#[cfg(test)]
mod tests;
