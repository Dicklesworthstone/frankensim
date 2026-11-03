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

/// Reciprocal response of both impermeable ends about a supplied midpoint
/// velocity: `p_i = p_reference_i - sum_j R_ij * (U_j - U_reference_j)`.
/// All velocities and half-cell forces are outward, ordered [left, right].
/// The matrix is symmetric; its cross term retains communication between ends.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpermeableEndsResponse {
    pub reference_pressure_pa: [f64; 2],
    pub reference_velocity_outward_m_s: [f64; 2],
    pub pressure_velocity_impedance_pa_s_m: [[f64; 2]; 2],
    pub boundary_half_mass_kg: f64,
    /// Authored body force on each endpoint fluid half-cell, excluding pressure.
    pub boundary_body_force_outward_n: [f64; 2],
}

impl TwoFluidChannel {
    /// Fluid inertia excluded from the resolved interior at EACH end [kg].
    /// An impermeable moving boundary must carry this momentum and energy.
    pub fn boundary_half_mass_kg(&self) -> f64 {
        (0.5 * self.cell_volume) * self.material.density_kg_m3()
    }

    /// Authored total component force on the two endpoint half dual volumes.
    /// The returned forces point outward, with left negative in channel +x.
    pub fn boundary_body_force_outward_n(
        &self,
        force: BodyForce<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 2], ChannelError> {
        force.validate(self.cells() + 1, &self.material, &mut cancelled)?;
        let volume = 0.5 * self.cell_volume;
        let values = [
            -volume * force.components(0, &self.material).2,
            volume * force.components(self.cells(), &self.material).2,
        ];
        if !values.into_iter().all(f64::is_finite) {
            return Err(ChannelError::Invalid(
                "endpoint body force is not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(values)
    }

    /// Query the exact coupled endpoint response on the original two-mode
    /// midpoint operator, including authored force kicks and mutual friction.
    ///
    /// A force-only mechanical midpoint is a useful reference: a rigidly
    /// translating or freely accelerating fluid then needs no cancellation of
    /// large artificial blocked pressures. This algebraic reference is never
    /// published and is not subject to physical amplitude limits. Finiteness,
    /// the Courant policy, original linear residuals, and work/momentum checks
    /// still apply. The actual `step_forced_with_heat` independently admits all
    /// physical limits before the coupled owners publish.
    ///
    /// Endpoint columns use the existing tridiagonal factors and buffers. They
    /// are cached per time step; there is no new allocation or wave solver.
    pub fn impermeable_ends_response(
        &mut self,
        dt_s: f64,
        reference_velocity_outward_m_s: [f64; 2],
        force: BodyForce<'_>,
        friction: MutualFriction,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ImpermeableEndsResponse, ChannelError> {
        let [left, right] = reference_velocity_outward_m_s;
        let boundary = BoundaryVelocities {
            left: EndVelocities {
                normal_m_s: -left,
                superfluid_m_s: -left,
            },
            right: EndVelocities {
                normal_m_s: right,
                superfluid_m_s: right,
            },
        };
        self.stage_step_with_limits(dt_s, boundary, friction, force, false, &mut cancelled)?;
        let mut pressure = [0.0; 2];
        for (slot, cell) in [0, self.cells() - 1].into_iter().enumerate() {
            let old = self.modes.cell(
                &self.material,
                [self.accepted[0].cells[cell], self.accepted[1].cells[cell]],
            );
            let new = self.modes.cell(
                &self.material,
                [self.pending[0].cells[cell], self.pending[1].cells[cell]],
            );
            pressure[slot] =
                0.5 * old.pressure_perturbation_pa + 0.5 * new.pressure_perturbation_pa;
        }
        let bulk_column = self.modes.face_to_modes(1.0, 0.0);
        let mut resistance = [[0.0; 2]; 2];
        let last = self.cells() - 1;
        for mode in 0..2 {
            let c = self.modes.speed[mode];
            let k = bulk_column[mode];
            let multiplier = (0.5 * (c * dt_s / self.dx) * k) * (c * k);
            for (column, end) in [0, last].into_iter().enumerate() {
                let inverse =
                    endpoint_inverse(&mut self.workspace[mode], column, end, &mut cancelled)?;
                for row in 0..2 {
                    resistance[row][column] += multiplier * inverse[row];
                }
            }
        }
        let scale = resistance[0][0].abs().max(resistance[1][1].abs());
        let asymmetry = (resistance[0][1] - resistance[1][0]).abs();
        let cross = 0.5 * resistance[0][1] + 0.5 * resistance[1][0];
        let half_mass = self.boundary_half_mass_kg();
        if !pressure.into_iter().all(f64::is_finite)
            || !resistance.into_iter().flatten().all(f64::is_finite)
            || !positive(half_mass)
            || !cross.is_finite()
            || cross < 0.0
            || resistance[0][0] <= cross
            || resistance[1][1] <= cross
            || asymmetry > 8192.0 * f64::EPSILON * scale
        {
            return Err(ChannelError::Invalid(
                "coupled impermeable-end response is unresolved",
            ));
        }
        // Enforce reciprocal arithmetic only after checking the independently
        // solved original endpoint columns. A distant cross term may underflow.
        resistance[0][1] = cross;
        resistance[1][0] = cross;
        let body_force = self.boundary_body_force_outward_n(force, &mut cancelled)?;
        Ok(ImpermeableEndsResponse {
            reference_pressure_pa: pressure,
            reference_velocity_outward_m_s,
            pressure_velocity_impedance_pa_s_m: resistance,
            boundary_half_mass_kg: half_mass,
            boundary_body_force_outward_n: body_force,
        })
    }

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
            resistance += (0.5 * (c * dt_s / self.dx) * k) * (c * k) * inverse[slot];
        }
        let half_mass = self.boundary_half_mass_kg();
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
) -> Result<[f64; 2], ChannelError> {
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
    let value = [work.rhs[0], work.rhs[work.rhs.len() - 1]];
    if !positive(value[slot])
        || !value[1 - slot].is_finite()
        || value[1 - slot] < 0.0
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
