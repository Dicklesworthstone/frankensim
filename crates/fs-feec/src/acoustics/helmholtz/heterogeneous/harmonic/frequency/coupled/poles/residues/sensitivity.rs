//! Original wall-material directions through the complete pole linearization.

use super::super::{MechanicalPoleError, MechanicalPoleLimits, cap, finite, poll};
use super::{PreparedCoupledResonance, response_error};
use fs_math::{c64::C64, dd::Dd};

/// Derivatives of one ORIGINAL positive Maxwell branch per caller parameter.
/// These are physical stiffness and relaxation-time directions, not derivatives
/// of a frequency-fitted loss factor or of the internal square-root coordinate.
#[derive(Clone, Copy, Debug, Default)]
pub struct MechanicalRelaxationParameterDirection {
    pub stiffness_pa_m: f64,
    pub relaxation_time_s: f64,
}

/// One patch's physical coefficient derivatives per caller parameter.
/// Geometry, trace mass, bulk fluid and all unselected patches remain fixed.
/// Signed directions are local derivatives; no finite parameter step or new
/// material validity domain is admitted by supplying them.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalWallParameterDirection<'a> {
    pub wall: usize,
    pub surface_density_kg_m2: f64,
    pub stiffness_pa_m: f64,
    pub resistance_pa_s_m: f64,
    /// Exactly one entry per original branch, in the original branch order.
    pub relaxation: &'a [MechanicalRelaxationParameterDirection],
}

/// Complete physical wall-coordinate derivatives in the fixed modal gauge.
#[derive(Clone, Debug)]
pub struct MechanicalWallParameterDerivative {
    pub displacement_m_per_parameter: Vec<C64>,
    pub normal_velocity_m_s_per_parameter: Vec<C64>,
    /// Includes the derivative of the source's sqrt(stiffness) coordinate.
    pub relaxation_memory_per_parameter: Vec<Vec<C64>>,
    pub relaxation_extension_per_parameter: Vec<Vec<C64>>,
}

/// Derivative of an isolated simple damped pole and its COMPLETE original mode.
/// Re(s') changes growth/decay; Im(s') changes signed oscillation frequency.
/// Neither component is a driven peak or a half-power bandwidth derivative.
#[derive(Clone, Debug)]
pub struct MechanicalResonanceParameterDerivative {
    pub rate_per_s_per_parameter: C64,
    /// Same complete free-coordinate order as the prepared owner's mode().
    pub coordinate_derivative: Vec<C64>,
    pub pressure_pa_per_parameter: Vec<C64>,
    pub walls: Vec<MechanicalWallParameterDerivative>,
    /// This physical coordinate is held at its BASE mode() value. The largest
    /// modal amplitude is one at the base, but that norm does not move with the
    /// parameter. This avoids a nonsmooth max-norm/phase derivative.
    pub normalization_pivot: usize,
    pub normalization_value: C64,
    pub bordered_relative_residual: f64,
    /// Independently checked T*x' + T_s*x*s' + T_parameter*x in original rows.
    pub physical_derivative_residual: f64,
}

impl PreparedCoupledResonance<'_> {
    pub(crate) fn wall_parameter_derivative(
        &self,
        directions: &[MechanicalWallParameterDirection<'_>],
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonanceParameterDerivative, MechanicalPoleError> {
        poll(cancelled)?;
        let source = self.source;
        let full_count = source.blocks[0].nrows();
        let n = source.free_dofs.len();
        if limits.max_evaluations == 0
            || !limits.physical_residual_tolerance.is_finite()
            || limits.physical_residual_tolerance <= 0.0
            || directions.len() > source.walls.len()
        {
            return Err(bad(
                "bounded complete wall direction and positive audit limits required",
            ));
        }
        // One original derivative-row audit; the already prepared shared owner
        // separately bounds its single RHS solve. All additional storage is O(N).
        cap(
            "wall parameter derivative payload bytes",
            320 * full_count as u128 + 160 * n as u128 + 32 * source.walls.len() as u128,
            limits.max_payload_bytes as u128,
        )?;
        let mut selected = vec![None; source.walls.len()];
        for direction in directions {
            poll(cancelled)?;
            if direction.wall >= source.walls.len()
                || selected[direction.wall].is_some()
                || direction.relaxation.len() != source.walls[direction.wall].arms.len()
                || [
                    direction.surface_density_kg_m2,
                    direction.stiffness_pa_m,
                    direction.resistance_pa_s_m,
                ]
                .iter()
                .any(|x| !x.is_finite())
                || direction
                    .relaxation
                    .iter()
                    .any(|x| !x.stiffness_pa_m.is_finite() || !x.relaxation_time_s.is_finite())
            {
                return Err(bad(
                    "finite unique wall directions with every original branch required",
                ));
            }
            selected[direction.wall] = Some(direction);
        }
        let s = self.numerical.value();
        let mut full = vec![C64::ZERO; full_count];
        let mut amplitude = 0.0_f64;
        for ((&row, &scale), &y) in source
            .free_dofs
            .iter()
            .zip(&source.scales)
            .zip(self.numerical.right_vector())
        {
            poll(cancelled)?;
            full[row] = y.scale(scale);
            amplitude = amplitude.max(full[row].abs());
        }
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(bad("physical pole normalization is unrepresentable"));
        }
        let mut action = vec![[Dd::ZERO; 2]; full_count];
        for (wall_index, wall) in source.walls.iter().enumerate() {
            let Some(direction) = selected[wall_index] else {
                continue;
            };
            let offset = source.offsets[wall_index];
            let count = wall.dofs.len();
            for row in 0..count {
                poll(cancelled)?;
                let (columns, weights) = wall.mass.row(row);
                for (&col, &weight) in columns.iter().zip(weights) {
                    add(
                        &mut action[offset + count + row],
                        C64::from_re(weight),
                        full[offset + col].scale(direction.stiffness_pa_m),
                    );
                    add(
                        &mut action[offset + count + row],
                        C64::from_re(weight),
                        (C64::from_re(direction.resistance_pa_s_m)
                            + s.scale(direction.surface_density_kg_m2))
                            * full[offset + count + col],
                    );
                    for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                        let (dc, _) = branch_direction(coupling, rate, direction.relaxation[arm])?;
                        add(
                            &mut action[offset + count + row],
                            C64::from_re(weight * dc),
                            full[offset + (2 + arm) * count + col],
                        );
                    }
                }
                for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                    let (dc, dr) = branch_direction(coupling, rate, direction.relaxation[arm])?;
                    let memory = offset + (2 + arm) * count + row;
                    add(
                        &mut action[memory],
                        C64::from_re(-dc),
                        full[offset + count + row],
                    );
                    add(&mut action[memory], C64::from_re(dr), full[memory]);
                }
            }
        }
        // Numerical congruence is fixed: T_hat=D*T*D, x=D*y, hence
        // T_hat_parameter*y=D*(T_parameter*x). No equilibration is differentiated.
        let mut rhs = Vec::with_capacity(n);
        for (&row, &scale) in source.free_dofs.iter().zip(&source.scales) {
            poll(cancelled)?;
            let value = C64::new(
                (action[row][0] * Dd::from_f64(scale)).to_f64(),
                (action[row][1] * Dd::from_f64(scale)).to_f64(),
            );
            if !finite(value) {
                return Err(bad("original wall parameter action is unrepresentable"));
            }
            rhs.push(value);
        }
        let tangent = self
            .numerical
            .parameter_derivatives(&[rhs.as_slice()], &mut *cancelled)
            .map_err(response_error)?
            .pop()
            .ok_or_else(|| bad("pole tangent was not returned"))?;
        let mut derivative = vec![C64::ZERO; full_count];
        let mut coordinates = Vec::with_capacity(n);
        for ((&row, &scale), &y) in source
            .free_dofs
            .iter()
            .zip(&source.scales)
            .zip(&tangent.vector)
        {
            poll(cancelled)?;
            let value = divide(y.scale(scale), amplitude);
            if !finite(value) {
                return Err(bad("original physical mode tangent is unrepresentable"));
            }
            derivative[row] = value;
            coordinates.push(value);
        }
        for (value, force) in full.iter_mut().zip(&mut action) {
            *value = divide(*value, amplitude);
            force[0] = force[0] / Dd::from_f64(amplitude);
            force[1] = force[1] / Dd::from_f64(amplitude);
        }
        let derivative_norm = coordinates.iter().map(|x| x.abs()).fold(0.0, f64::max);
        let mut physical_residual = 0.0_f64;
        for &row in &source.free_dofs {
            poll(cancelled)?;
            let mut sum = action[row];
            let mut envelope = C64::new(sum[0].to_f64(), sum[1].to_f64()).abs();
            for (block, multiplier, slope) in [
                (&source.blocks[0], C64::ONE, C64::ZERO),
                (&source.blocks[2], s, C64::ONE),
                (&source.blocks[1], s * s, s.scale(2.0)),
            ] {
                let (columns, values) = block.row(row);
                for (&col, &value) in columns.iter().zip(values) {
                    add(&mut sum, multiplier.scale(value), derivative[col]);
                    add(&mut sum, (slope * tangent.value).scale(value), full[col]);
                    envelope += value.abs()
                        * (multiplier.abs() * derivative_norm + slope.abs() * tangent.value.abs());
                }
            }
            let error = C64::new(sum[0].to_f64(), sum[1].to_f64()).abs();
            let residual = if envelope == 0.0 && error == 0.0 {
                0.0
            } else {
                error / envelope
            };
            if !envelope.is_finite()
                || !residual.is_finite()
                || residual > limits.physical_residual_tolerance
            {
                return Err(MechanicalPoleError::Residual {
                    row,
                    observed: residual,
                    maximum: limits.physical_residual_tolerance,
                });
            }
            physical_residual = physical_residual.max(residual);
        }
        let mut walls = Vec::with_capacity(source.walls.len());
        for (wall_index, wall) in source.walls.iter().enumerate() {
            poll(cancelled)?;
            let count = wall.dofs.len();
            let offset = source.offsets[wall_index];
            let mut memories = Vec::with_capacity(wall.arms.len());
            let mut extensions = Vec::with_capacity(wall.arms.len());
            for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                let direction = selected[wall_index]
                    .map_or(MechanicalRelaxationParameterDirection::default(), |d| {
                        d.relaxation[arm]
                    });
                let (dc, dr) = branch_direction(coupling, rate, direction)?;
                let start = offset + (2 + arm) * count;
                let mut memory = Vec::with_capacity(count);
                for j in 0..count {
                    poll(cancelled)?;
                    let mut value = [Dd::ZERO; 2];
                    if s.abs() > rate {
                        // Differentiate z=(rate/s)e to retain very small physical
                        // memory at high rate. No s+rate denominator appears.
                        let z = (full[start + j] / s).scale(rate);
                        add(&mut value, C64::from_re(rate), derivative[start + j] / s);
                        add(&mut value, C64::from_re(dr), full[start + j] / s);
                        add(&mut value, -tangent.value, z / s);
                    } else {
                        // z=c*w-e, with c=sqrt(k). The dc*w term is essential
                        // when the material stiffness, rather than c, is varied.
                        add(&mut value, C64::from_re(coupling), derivative[offset + j]);
                        add(&mut value, C64::from_re(dc), full[offset + j]);
                        add(&mut value, C64::from_re(-1.0), derivative[start + j]);
                    }
                    let value = C64::new(value[0].to_f64(), value[1].to_f64());
                    if !finite(value) {
                        return Err(bad("original Maxwell memory tangent is unrepresentable"));
                    }
                    memory.push(value);
                }
                memories.push(memory);
                extensions.push(derivative[start..start + count].to_vec());
            }
            walls.push(MechanicalWallParameterDerivative {
                displacement_m_per_parameter: derivative[offset..offset + count].to_vec(),
                normal_velocity_m_s_per_parameter: derivative[offset + count..offset + 2 * count]
                    .to_vec(),
                relaxation_memory_per_parameter: memories,
                relaxation_extension_per_parameter: extensions,
            });
        }
        let pivot = self.numerical.normalization_pivot();
        poll(cancelled)?;
        Ok(MechanicalResonanceParameterDerivative {
            rate_per_s_per_parameter: tangent.value,
            coordinate_derivative: coordinates,
            pressure_pa_per_parameter: derivative[..source.fluid_dofs].to_vec(),
            walls,
            normalization_pivot: pivot,
            normalization_value: self.mode.coordinates[pivot],
            bordered_relative_residual: tangent.linear_residual,
            physical_derivative_residual: physical_residual,
        })
    }
}

fn branch_direction(
    coupling: f64,
    rate: f64,
    d: MechanicalRelaxationParameterDirection,
) -> Result<(f64, f64), MechanicalPoleError> {
    // Halving a subnormal physical stiffness direction first can erase it,
    // although division by the source's small sqrt(k) recovers a normal value.
    let dc = if d.stiffness_pa_m.abs() < f64::MIN_POSITIVE {
        (d.stiffness_pa_m / coupling) * 0.5
    } else {
        (0.5 * d.stiffness_pa_m) / coupling
    };
    let dr = -(d.relaxation_time_s * rate) * rate;
    if !dc.is_finite() || !dr.is_finite() {
        Err(bad(
            "original Maxwell parameter direction is unrepresentable",
        ))
    } else {
        Ok((dc, dr))
    }
}

fn bad(message: &'static str) -> MechanicalPoleError {
    MechanicalPoleError::Invalid(message)
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn add(sum: &mut [Dd; 2], a: C64, b: C64) {
    sum[0] =
        sum[0] + Dd::from_f64(a.re) * Dd::from_f64(b.re) - Dd::from_f64(a.im) * Dd::from_f64(b.im);
    sum[1] =
        sum[1] + Dd::from_f64(a.re) * Dd::from_f64(b.im) + Dd::from_f64(a.im) * Dd::from_f64(b.re);
}

#[cfg(test)]
mod tests;
