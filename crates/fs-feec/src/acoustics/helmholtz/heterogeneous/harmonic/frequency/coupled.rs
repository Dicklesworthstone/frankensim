//! Original pressure, distributed wall motion and every Maxwell memory in ONE
//! harmonic pencil. No reciprocal wall impedance or frequency-dependent Schur
//! complement is formed. The invertible memory coordinate e=sqrt(k) w-z keeps
//! small physical relaxation extensions explicit, as in the exact time owner.

use super::{
    Forcing, PressureError, PressurePencil, PressureSolution, finite, poll, quadratic_nonnegative,
    real_form,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::{Coo, Csr};

pub(crate) mod poles;

/// Dense admission counts every free pressure, wall displacement/velocity and
/// Maxwell coordinate. Original sparse forms and lifted fields are additional.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalHarmonicLimits {
    pub max_states: usize,
    /// Peak dense input plus LU, 32*n^2 bytes; original sparse forms, O(n)
    /// vectors and a configured parallel LU worker's workspace are additional.
    pub max_dense_bytes: usize,
    /// Componentwise backward error of the complete original equations.
    pub residual_tolerance: f64,
}

impl Default for MechanicalHarmonicLimits {
    fn default() -> Self {
        Self {
            max_states: 1024,
            max_dense_bytes: 64 * 1024 * 1024,
            residual_tolerance: 5e-13,
        }
    }
}

/// Complex peak phasors in the patch's original sorted wall-DOF order.
#[derive(Clone, Debug)]
pub struct MechanicalWallPhasors {
    pub displacement_m: Vec<C64>,
    pub normal_velocity_m_s: Vec<C64>,
    /// Original energy-normalized z_j, with e_j=sqrt(k_j) w-z_j.
    pub relaxation_memory: Vec<Vec<C64>>,
}

/// Analytic derivatives of the same original wall phasors per rad/s.
#[derive(Clone, Debug)]
pub struct MechanicalWallFrequencyDerivative {
    pub displacement_m_per_rad_s: Vec<C64>,
    pub normal_velocity_m_s_per_rad_s: Vec<C64>,
    pub relaxation_memory_per_rad_s: Vec<Vec<C64>>,
}

/// Time-averaged wall storage, irreversible heat and signed fluid-to-wall work.
/// Energies use peak-phasor factor 1/4; powers use factor 1/2.
#[derive(Clone, Copy, Debug, Default)]
pub struct MechanicalWallAccounting {
    pub kinetic_energy_j: f64,
    pub elastic_energy_j: f64,
    pub relaxation_energy_j: f64,
    pub resistance_dissipation_w: f64,
    pub relaxation_dissipation_w: f64,
    pub fluid_to_wall_power_w: f64,
    /// Signed external mechanical work, independent of fluid source work and
    /// internal pressure exchange: 1/2 Re((W*traction)^H v).
    pub wall_traction_power_w: f64,
}

impl MechanicalWallAccounting {
    #[must_use]
    pub fn total_energy_j(self) -> f64 {
        self.kinetic_energy_j + self.elastic_energy_j + self.relaxation_energy_j
    }
    #[must_use]
    pub fn dissipation_w(self) -> f64 {
        self.resistance_dissipation_w + self.relaxation_dissipation_w
    }
    /// Fluid exchange plus external mechanical work minus actual wall heat.
    #[must_use]
    pub fn power_balance_residual_w(self) -> f64 {
        self.fluid_to_wall_power_w + self.wall_traction_power_w - self.dissipation_w()
    }
}

pub(crate) struct Wall {
    pub(crate) dofs: Vec<usize>,
    pub(crate) mass: Csr,
    pub(crate) surface_density: f64,
    pub(crate) stiffness: f64,
    pub(crate) resistance: f64,
    /// Original (sqrt(k), 1/tau), admitted through the wall's material owner.
    pub(crate) arms: Vec<(f64, f64)>,
}

/// One cached, signed consistent W*shape force column in its original wall.
/// The public prepared source additionally binds this to its exact owner.
pub(crate) struct WallTraction {
    wall: usize,
    force: Vec<Dd>,
}

pub(crate) struct CoupledPencil {
    omega: f64,
    fluid_dofs: usize,
    blocks: [Csr; 4], // A0, A2, A1, B in A0-omega² A2+i omega A1-omega B.
    free_dofs: Vec<usize>,
    scales: Vec<f64>,
    factor: LuComplex,
    walls: Vec<Wall>,
    offsets: Vec<usize>,
    tolerance: f64,
}

impl CoupledPencil {
    pub(crate) fn prepare(
        omega: f64,
        fluid: [&Csr; 4],
        free_pressure: &[usize],
        walls: Vec<Wall>,
        limits: MechanicalHarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PressureError> {
        poll(cancelled)?;
        let fluid_dofs = fluid[0].nrows();
        if fluid_dofs == 0
            || fluid
                .iter()
                .any(|a| a.nrows() != fluid_dofs || a.ncols() != fluid_dofs)
            || free_pressure.iter().any(|&i| i >= fluid_dofs)
            || free_pressure.windows(2).any(|ij| ij[0] >= ij[1])
            || !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
            || limits.max_states == 0
            || limits.max_dense_bytes == 0
            || !limits.residual_tolerance.is_finite()
            || limits.residual_tolerance <= 0.0
        {
            return Err(PressureError::Invalid(
                "invalid complete harmonic pencil or limits",
            ));
        }
        let mut dimension = fluid_dofs;
        let mut states = free_pressure.len();
        let mut offsets = Vec::with_capacity(walls.len());
        for wall in &walls {
            poll(cancelled)?;
            let n = wall.dofs.len();
            if n == 0
                || wall.mass.nrows() != n
                || wall.mass.ncols() != n
                || wall.dofs.iter().any(|&i| i >= fluid_dofs)
                || wall.dofs.windows(2).any(|ij| ij[0] >= ij[1])
                || ![wall.surface_density, wall.stiffness]
                    .iter()
                    .all(|x| x.is_finite() && *x > 0.0)
                || !wall.resistance.is_finite()
                || wall.resistance < 0.0
                || wall.arms.iter().any(|&(c, rate)| {
                    !c.is_finite() || c <= 0.0 || !rate.is_finite() || rate <= 0.0
                })
            {
                return Err(PressureError::Invalid(
                    "invalid original mechanical trace or coefficients",
                ));
            }
            offsets.push(dimension);
            let count = wall
                .arms
                .len()
                .checked_add(2)
                .and_then(|m| n.checked_mul(m))
                .ok_or(PressureError::Invalid("coupled wall state extent overflow"))?;
            dimension = dimension
                .checked_add(count)
                .ok_or(PressureError::Invalid("coupled state extent overflow"))?;
            states = states
                .checked_add(count)
                .ok_or(PressureError::Invalid("coupled free-state extent overflow"))?;
        }
        if states == 0 {
            return Err(PressureError::Invalid(
                "complete harmonic pencil has no free states",
            ));
        }
        if states > limits.max_states {
            return Err(PressureError::Limit {
                what: "complete pressure and mechanical wall states",
                requested: states,
                maximum: limits.max_states,
            });
        }
        let entries = states
            .checked_mul(states)
            .ok_or(PressureError::Invalid("coupled dense extent overflow"))?;
        let bytes = entries
            .checked_mul(2 * core::mem::size_of::<C64>())
            .ok_or(PressureError::Invalid("coupled dense byte extent overflow"))?;
        if bytes > limits.max_dense_bytes {
            return Err(PressureError::Limit {
                what: "complete harmonic dense matrix bytes",
                requested: bytes,
                maximum: limits.max_dense_bytes,
            });
        }
        let mut coo: [Coo; 4] = core::array::from_fn(|_| Coo::new(dimension, dimension));
        for (destination, original) in coo.iter_mut().zip(fluid) {
            for row in 0..fluid_dofs {
                poll(cancelled)?;
                let (cols, values) = original.row(row);
                for (&col, &value) in cols.iter().zip(values) {
                    destination.push(row, col, value);
                }
            }
        }
        for (wall, &offset) in walls.iter().zip(&offsets) {
            let n = wall.dofs.len();
            let velocity = offset + n;
            for (row, &fluid_row) in wall.dofs.iter().enumerate() {
                poll(cancelled)?;
                let (columns, values) = wall.mass.row(row);
                for (&col, &w) in columns.iter().zip(values) {
                    // Reciprocal original-face work: +W v in continuity,
                    // +W p forcing the mechanical momentum equation.
                    coo[2].push(fluid_row, velocity + col, w);
                    coo[0].push(velocity + row, wall.dofs[col], -w);
                    coo[0].push(velocity + row, offset + col, wall.stiffness * w);
                    coo[0].push(velocity + row, velocity + col, wall.resistance * w);
                    coo[2].push(velocity + row, velocity + col, wall.surface_density * w);
                    for (arm, &(coupling, _)) in wall.arms.iter().enumerate() {
                        coo[0].push(velocity + row, offset + (2 + arm) * n + col, coupling * w);
                    }
                }
                coo[0].push(offset + row, velocity + row, 1.0);
                coo[2].push(offset + row, offset + row, -1.0); // v=i omega w.
                for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                    let memory = offset + (2 + arm) * n + row;
                    coo[0].push(memory, velocity + row, -coupling);
                    coo[0].push(memory, memory, rate);
                    coo[2].push(memory, memory, 1.0); // e_dot=c v-rate e.
                }
            }
        }
        let blocks = coo.map(|a| a.assemble());
        let mut free_dofs = free_pressure.to_vec();
        free_dofs.extend(fluid_dofs..dimension);
        let mut inverse = vec![usize::MAX; dimension];
        for (row, &original) in free_dofs.iter().enumerate() {
            inverse[original] = row;
        }
        let mut dense = vec![C64::ZERO; entries];
        let mut scales = vec![0.0; states];
        for (row, &original) in free_dofs.iter().enumerate() {
            poll(cancelled)?;
            for (matrix, multiplier, imaginary) in [
                (&blocks[0], 1.0, false),
                (&blocks[1], -omega * omega, false),
                (&blocks[3], -omega, false),
                (&blocks[2], omega, true),
            ] {
                let (columns, values) = matrix.row(original);
                for (&col, &value) in columns.iter().zip(values) {
                    if inverse[col] != usize::MAX {
                        let entry = &mut dense[row * states + inverse[col]];
                        if imaginary {
                            entry.im += multiplier * value;
                        } else {
                            entry.re += multiplier * value;
                        }
                    }
                }
            }
            let values = &dense[row * states..(row + 1) * states];
            if values.iter().any(|&v| !finite(v)) {
                return Err(PressureError::Invalid(
                    "complete harmonic operator overflow",
                ));
            }
            let maximum = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if maximum == 0.0 {
                return Err(LuComplexError::Singular { index: row }.into());
            }
            scales[row] = 1.0 / det::sqrt(maximum);
        }
        for row in 0..states {
            poll(cancelled)?;
            for col in 0..states {
                dense[row * states + col] = dense[row * states + col]
                    .scale(scales[row])
                    .scale(scales[col]);
            }
        }
        if dense.iter().any(|&v| !finite(v)) || scales.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return Err(PressureError::Invalid(
                "complete harmonic equilibration overflow",
            ));
        }
        let factor = lu_complex_with_cancel(&dense, states, &mut *cancelled)?;
        poll(cancelled)?;
        Ok(Self {
            omega,
            fluid_dofs,
            blocks,
            free_dofs,
            scales,
            factor,
            walls,
            offsets,
            tolerance: limits.residual_tolerance,
        })
    }

    fn pencil(&self) -> PressurePencil<'_> {
        PressurePencil {
            omega: self.omega,
            stiffness: &self.blocks[0],
            mass: &self.blocks[1],
            conductance: &self.blocks[2],
            susceptance: &self.blocks[3],
            boundary_frequency_derivative: None,
            factor: &self.factor,
            scales: &self.scales,
            free_dofs: &self.free_dofs,
            residual_tolerance: self.tolerance,
        }
    }

    fn load(&self, load: &[C64]) -> Result<Vec<C64>, PressureError> {
        if load.len() != self.fluid_dofs {
            return Err(PressureError::Invalid(
                "original acoustic load extent required",
            ));
        }
        let mut full = load.to_vec();
        full.resize(self.blocks[0].nrows(), C64::ZERO);
        Ok(full)
    }

    pub(crate) fn solve(
        &self,
        load: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.pencil().solve(&self.load(load)?, cancelled)
    }

    pub(crate) fn frequency_derivative(
        &self,
        state: &[C64],
        load: &[C64],
        derivative: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.pencil().frequency_derivative(
            state,
            &self.load(load)?,
            &self.load(derivative)?,
            cancelled,
        )
    }

    pub(crate) fn wall_dofs(&self, wall: usize) -> Option<&[usize]> {
        self.walls.get(wall).map(|w| w.dofs.as_slice())
    }

    pub(crate) fn prepare_traction(
        &self,
        wall: usize,
        coefficients: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallTraction, PressureError> {
        let original = self
            .walls
            .get(wall)
            .ok_or(PressureError::Invalid("traction wall index"))?;
        if coefficients.len() != original.dofs.len() || coefficients.iter().any(|v| !v.is_finite())
        {
            return Err(PressureError::Invalid(
                "original finite wall traction shape required",
            ));
        }
        let mut force = vec![Dd::ZERO; coefficients.len()];
        for (row, output) in force.iter_mut().enumerate() {
            poll(cancelled)?;
            let (columns, values) = original.mass.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                *output = *output + Dd::from_f64(value) * Dd::from_f64(coefficients[column]);
            }
            if !output.hi.is_finite() || !output.lo.is_finite() {
                return Err(PressureError::Invalid("wall traction projection overflow"));
            }
        }
        poll(cancelled)?;
        Ok(WallTraction { wall, force })
    }

    /// Full physical RHS: fluid i*omega*b and wall W*traction. Supplying a
    /// fluid slope requests rhs' = i*b + i*omega*b' and W*traction', without
    /// dividing mechanical forcing by frequency or differentiating that ratio.
    fn traction_forcing(
        &self,
        load: &[C64],
        slope: Option<&[C64]>,
        tractions: &[(&WallTraction, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Forcing, PressureError> {
        poll(cancelled)?;
        if load.len() != self.fluid_dofs
            || load.iter().any(|&v| !finite(v))
            || slope.is_some_and(|s| s.len() != load.len() || s.iter().any(|&v| !finite(v)))
        {
            return Err(PressureError::Invalid(
                "original finite acoustic source and slope required",
            ));
        }
        let n = self.blocks[0].nrows();
        let mut forcing = Forcing {
            values: vec![[Dd::ZERO; 2]; n],
            magnitude: vec![0.0; n],
        };
        for (row, &value) in load.iter().enumerate() {
            poll(cancelled)?;
            let (rhs, magnitude) = if let Some(slope) = slope {
                let w = Dd::from_f64(self.omega);
                (
                    [
                        -Dd::from_f64(value.im) - w * Dd::from_f64(slope[row].im),
                        Dd::from_f64(value.re) + w * Dd::from_f64(slope[row].re),
                    ],
                    value.abs() + self.omega * slope[row].abs(),
                )
            } else {
                let rhs = C64::new(-self.omega * value.im, self.omega * value.re);
                ([Dd::from_f64(rhs.re), Dd::from_f64(rhs.im)], rhs.abs())
            };
            forcing.values[row] = rhs;
            forcing.magnitude[row] = magnitude;
        }
        for &(traction, amplitude) in tractions {
            poll(cancelled)?;
            if !finite(amplitude) {
                return Err(PressureError::Invalid(
                    "finite wall traction amplitude required",
                ));
            }
            let start = self.offsets[traction.wall] + traction.force.len();
            for (local, &force) in traction.force.iter().enumerate() {
                poll(cancelled)?;
                let row = start + local;
                forcing.values[row][0] =
                    forcing.values[row][0] + force * Dd::from_f64(amplitude.re);
                forcing.values[row][1] =
                    forcing.values[row][1] + force * Dd::from_f64(amplitude.im);
                forcing.magnitude[row] += force.to_f64().abs() * amplitude.abs();
            }
        }
        Ok(forcing)
    }

    pub(crate) fn solve_with_tractions(
        &self,
        load: &[C64],
        tractions: &[(&WallTraction, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PressureSolution, f64), PressureError> {
        if tractions.is_empty() {
            return Ok((self.solve(load, cancelled)?, 0.0));
        }
        let forcing = self.traction_forcing(load, None, tractions, cancelled)?;
        let solution = self.pencil().solve_rhs(&forcing, cancelled)?;
        let mut power = Dd::ZERO;
        for (wall, &offset) in self.walls.iter().zip(&self.offsets) {
            let n = wall.dofs.len();
            for row in offset + n..offset + 2 * n {
                poll(cancelled)?;
                let v = solution.pressure[row];
                power = power
                    + Dd::from_f64(0.5)
                        * (forcing.values[row][0] * Dd::from_f64(v.re)
                            + forcing.values[row][1] * Dd::from_f64(v.im));
            }
        }
        if !power.hi.is_finite() || !power.lo.is_finite() {
            return Err(PressureError::Invalid("wall traction work overflow"));
        }
        poll(cancelled)?;
        Ok((solution, power.to_f64()))
    }

    pub(crate) fn frequency_derivative_with_tractions(
        &self,
        state: &[C64],
        load: &[C64],
        slope: &[C64],
        traction_slopes: &[(&WallTraction, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        if traction_slopes.is_empty() {
            return self.frequency_derivative(state, load, slope, cancelled);
        }
        let forcing = self.traction_forcing(load, Some(slope), traction_slopes, cancelled)?;
        self.pencil()
            .frequency_derivative_rhs(state, forcing, cancelled)
    }

    pub(crate) fn wall_fields(
        &self,
        state: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<MechanicalWallPhasors>, MechanicalWallAccounting), PressureError> {
        if state.len() != self.blocks[0].nrows() || state.iter().any(|&v| !finite(v)) {
            return Err(PressureError::Invalid(
                "complete finite wall state required",
            ));
        }
        let mut fields = Vec::with_capacity(self.walls.len());
        let mut account = [Dd::from_f64(0.0); 6];
        for (wall, &offset) in self.walls.iter().zip(&self.offsets) {
            poll(cancelled)?;
            let n = wall.dofs.len();
            let displacement = &state[offset..offset + n];
            let velocity = &state[offset + n..offset + 2 * n];
            let pressure: Vec<_> = wall.dofs.iter().map(|&i| state[i]).collect();
            let vv = quadratic_nonnegative(&wall.mass, velocity, cancelled)?;
            let ww = quadratic_nonnegative(&wall.mass, displacement, cancelled)?;
            account[0] = account[0] + Dd::from_f64(0.25 * wall.surface_density) * Dd::from_f64(vv);
            account[1] = account[1] + Dd::from_f64(0.25 * wall.stiffness) * Dd::from_f64(ww);
            account[3] = account[3] + Dd::from_f64(0.5 * wall.resistance) * Dd::from_f64(vv);
            account[5] = account[5]
                + Dd::from_f64(0.5 * real_form(&wall.mass, &pressure, velocity, cancelled)?.0);
            let mut memories = Vec::with_capacity(wall.arms.len());
            for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                let start = offset + (2 + arm) * n;
                let extension = &state[start..start + n];
                let ee = quadratic_nonnegative(&wall.mass, extension, cancelled)?;
                account[2] = account[2] + Dd::from_f64(0.25 * ee);
                account[4] = account[4] + Dd::from_f64(0.5 * rate) * Dd::from_f64(ee);
                memories.push(if self.omega > rate {
                    extension
                        .iter()
                        .map(|&e| memory_from_extension(e, rate, self.omega))
                        .collect()
                } else {
                    original_memory(displacement, extension, coupling)
                });
            }
            fields.push(MechanicalWallPhasors {
                displacement_m: displacement.to_vec(),
                normal_velocity_m_s: velocity.to_vec(),
                relaxation_memory: memories,
            });
        }
        let a = account.map(Dd::to_f64);
        if a.iter().any(|v| !v.is_finite())
            || fields
                .iter()
                .any(|w| w.relaxation_memory.iter().flatten().any(|&v| !finite(v)))
        {
            return Err(PressureError::Invalid(
                "unrepresentable mechanical wall storage, heat or memory",
            ));
        }
        poll(cancelled)?;
        Ok((
            fields,
            MechanicalWallAccounting {
                kinetic_energy_j: a[0],
                elastic_energy_j: a[1],
                relaxation_energy_j: a[2],
                resistance_dissipation_w: a[3],
                relaxation_dissipation_w: a[4],
                fluid_to_wall_power_w: a[5],
                wall_traction_power_w: 0.0,
            },
        ))
    }

    pub(crate) fn wall_derivatives(
        &self,
        state: &[C64],
        tangent: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MechanicalWallFrequencyDerivative>, PressureError> {
        if tangent.len() != self.blocks[0].nrows()
            || state.len() != tangent.len()
            || state.iter().chain(tangent).any(|&v| !finite(v))
        {
            return Err(PressureError::Invalid(
                "complete finite wall tangent required",
            ));
        }
        let mut fields = Vec::with_capacity(self.walls.len());
        for (wall, &offset) in self.walls.iter().zip(&self.offsets) {
            poll(cancelled)?;
            let n = wall.dofs.len();
            let displacement = &tangent[offset..offset + n];
            let mut memories = Vec::with_capacity(wall.arms.len());
            for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                let start = offset + (2 + arm) * n;
                let memory = if self.omega > rate {
                    tangent[start..start + n]
                        .iter()
                        .zip(&state[start..start + n])
                        .map(|(&de, &e)| {
                            let dz = memory_from_extension(de, rate, self.omega);
                            let z = memory_from_extension(e, rate, self.omega);
                            let omega = Dd::from_f64(self.omega);
                            // Differentiate i*omega*z=rate*e, including the explicit omega.
                            C64::new(
                                (Dd::from_f64(dz.re) - Dd::from_f64(z.re) / omega).to_f64(),
                                (Dd::from_f64(dz.im) - Dd::from_f64(z.im) / omega).to_f64(),
                            )
                        })
                        .collect()
                } else {
                    original_memory(displacement, &tangent[start..start + n], coupling)
                };
                if memory.iter().any(|&v| !finite(v)) {
                    return Err(PressureError::Invalid("wall memory tangent overflow"));
                }
                memories.push(memory);
            }
            fields.push(MechanicalWallFrequencyDerivative {
                displacement_m_per_rad_s: displacement.to_vec(),
                normal_velocity_m_s_per_rad_s: tangent[offset + n..offset + 2 * n].to_vec(),
                relaxation_memory_per_rad_s: memories,
            });
        }
        poll(cancelled)?;
        Ok(fields)
    }
}

// For omega*tau>1, z=c*w-e cancels nearly equal phasors. The original
// memory equation instead gives z=rate*e/(i*omega), without cancellation.
// Divide by omega/rate when representable, otherwise divide then multiply:
// do not underflow rate/omega before multiplying a large finite extension.
fn memory_from_extension(extension: C64, rate: f64, omega: f64) -> C64 {
    let rotated = C64::new(extension.im, -extension.re);
    let denominator = omega / rate;
    if denominator.is_finite() {
        rotated / C64::from_re(denominator)
    } else {
        (rotated / C64::from_re(omega)).scale(rate)
    }
}

fn original_memory(displacement: &[C64], extension: &[C64], coupling: f64) -> Vec<C64> {
    let c = Dd::from_f64(coupling);
    displacement
        .iter()
        .zip(extension)
        .map(|(w, e)| {
            C64::new(
                (c * Dd::from_f64(w.re) - Dd::from_f64(e.re)).to_f64(),
                (c * Dd::from_f64(w.im) - Dd::from_f64(e.im)).to_f64(),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests;
