//! Endpoint pressure response from the original two-mode midpoint operators.

use super::*;

/// End of the reference channel; outward points toward -x on the left.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChannelEnd {
    Left,
    Right,
}

/// Exact discrete pressure relation `p_mid = p_blocked - R * U_outward`.
/// Both components share the impermeable end velocity; the other end is fixed.
/// Pressure is a perturbation about a balanced reference pressure. The half-cell
/// fluid mass belongs in the mechanical composition's momentum and energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpermeableEndResponse {
    pub blocked_pressure_pa: f64,
    pub pressure_velocity_impedance_pa_s_m: f64,
    pub boundary_half_mass_kg: f64,
}

impl TwoFluidChannel {
    /// Query a reciprocal impermeable-end response without publishing any state.
    ///
    /// This reuses the original closed wave/friction proposal and its admitted
    /// tridiagonal factors. The two endpoint inverse diagonals are cached for
    /// each time step. No new allocation or alternative wave solver is used.
    /// The blocked proposal must satisfy the original physical admission; the
    /// final moving-boundary step must independently satisfy it too.
    ///
    /// Friction leaves this pressure response affine: its first half-step is
    /// independent of the new boundary velocity, while the second half-step
    /// changes no pressure or barycentric momentum. The complete coupled step
    /// still applies the original friction law to its actual moving solution.
    pub fn impermeable_end_response(
        &mut self,
        dt_s: f64,
        end: ChannelEnd,
        friction: MutualFriction,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ImpermeableEndResponse, ChannelError> {
        self.stage_step(
            dt_s,
            BoundaryVelocities::CLOSED,
            friction,
            BodyForce::None,
            &mut cancelled,
        )?;
        let slot = match end {
            ChannelEnd::Left => 0,
            ChannelEnd::Right => 1,
        };
        let boundary_cell = if slot == 0 { 0 } else { self.cells() - 1 };
        let old = self.modes.cell(
            &self.material,
            [
                self.accepted[0].cells[boundary_cell],
                self.accepted[1].cells[boundary_cell],
            ],
        );
        let new = self.modes.cell(
            &self.material,
            [
                self.pending[0].cells[boundary_cell],
                self.pending[1].cells[boundary_cell],
            ],
        );
        let pressure = 0.5 * old.pressure_perturbation_pa + 0.5 * new.pressure_perturbation_pa;
        let bulk_column = self.modes.face_to_modes(1.0, 0.0);
        let mut resistance = 0.0;
        for mode in 0..2 {
            let inverse = endpoint_inverse(
                &mut self.workspace[mode],
                slot,
                boundary_cell,
                &mut cancelled,
            )?;
            // In energy coordinates p=sum(c_m*k_m*P_m), k_m is the modal
            // velocity produced by unit barycentric speed. A prescribed end
            // adds -lambda*k_m*U to the midpoint RHS; midpoint pressure has
            // the additional factor 1/2. This is the original Schur response.
            let c = self.modes.speed[mode];
            let k = bulk_column[mode];
            resistance += (0.5 * (c * dt_s / self.dx) * k) * (c * k) * inverse;
        }
        let half_mass = (0.5 * self.cell_volume) * self.material.density_kg_m3();
        if !pressure.is_finite() || !positive(resistance) || !positive(half_mass) {
            return Err(ChannelError::Invalid(
                "impermeable-end response is not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ImpermeableEndResponse {
            blocked_pressure_pa: pressure,
            pressure_velocity_impedance_pa_s_m: resistance,
            boundary_half_mass_kg: half_mass,
        })
    }
}

fn endpoint_inverse(
    work: &mut SolveWorkspace,
    slot: usize,
    end: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ChannelError> {
    poll(cancelled)?;
    if let Some(value) = work.boundary_inverse[slot] {
        return Ok(value);
    }
    for chunk in work.rhs.chunks_mut(128) {
        poll(cancelled)?;
        chunk.fill(0.0);
    }
    work.rhs[end] = 1.0;
    work.factor
        .solve_factored_in_place(&mut work.rhs, &mut *cancelled)?;
    let mut error: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for i in 0..work.rhs.len() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let diagonal = work.diagonal[i] * work.rhs[i];
        let left = if i == 0 {
            0.0
        } else {
            work.off_diagonal[i - 1] * work.rhs[i - 1]
        };
        let right = if i + 1 == work.rhs.len() {
            0.0
        } else {
            work.off_diagonal[i] * work.rhs[i + 1]
        };
        let target = if i == end { 1.0 } else { 0.0 };
        error = error.max((diagonal + left + right - target).abs());
        scale = scale.max(diagonal.abs() + left.abs() + right.abs() + target);
    }
    let value = work.rhs[end];
    if !positive(value)
        || !error.is_finite()
        || !scale.is_finite()
        || error > 4096.0 * f64::EPSILON * scale
    {
        return Err(ChannelError::Invalid(
            "impermeable-end pressure response is unresolved",
        ));
    }
    poll(cancelled)?;
    work.boundary_inverse[slot] = Some(value);
    Ok(value)
}
