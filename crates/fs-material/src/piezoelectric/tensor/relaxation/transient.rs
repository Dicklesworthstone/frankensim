//! Retained passive midpoint memories of the original Maxwell/Debye material.
//!
//! This is a local constitutive companion, not a structural or electric solver.
//! The caller supplies delta_S=(S1-S0)/2 and E_mid from ONE coupled midpoint
//! solve. With h=dt/2, a=tau/(tau+h), b=h/(tau+h), the eliminated memories are
//! q_mid=a(q0+B delta_S), p_mid=a p0+b delta_epsilon (n.E_mid).
//! Thus a*C_j and b*delta_epsilon*n*n^T enter the original field operator.
//! Accepted endpoint fields include the retained old memories; the mechanical
//! history correction is +b*B^T q0, and the polarization increment independent
//! of delta_E is b*(delta_epsilon*n.E0-p0)*n.
//!
//! Every update is a candidate. Native Gram coordinates and oriented physical
//! polarizations remain bound to the exact source/model/frame; no public setter,
//! history projection, coefficient refresh or alternate time integrator exists.
//! Only unrestricted source bands admit time-domain use. A bounded harmonic fit
//! cannot establish applicability of an arbitrary transient spectrum.
//!
//! Physical heat uses the ACTUAL represented midpoint memory: dt*|q_mid|^2/tau
//! and dt*(delta_epsilon*n.E_mid-p_mid)^2/(tau*delta_epsilon). Stored-energy
//! changes use actual endpoint increments, preserving changes below the spacing
//! of a large stored energy. Work and the independent residual remain visible
//! for the enclosing spatial owner's energy tolerance; heat is never adjusted
//! to close that residual. Midpoint is second order and passive for this linear
//! constitutive law, but is not L-stable: unresolved fast memories can alternate.

use super::{RelaxingTensorPiezoelectric, Result, invalid, scalar};
use fs_blake3::ContentHash;
use fs_math::det;

/// Immutable accepted material memory. Mechanical components are the native
/// Gram stresses q [sqrt(Pa)]; polarization components p are in C/m² along each
/// original directed pole. Only the matching material can advance this state.
#[derive(Clone, Debug, PartialEq)]
pub struct PiezoelectricMemory {
    model_identity: ContentHash,
    mechanical: Vec<[f64; 6]>,
    dielectric: Vec<f64>,
}

impl PiezoelectricMemory {
    /// Exact identity of the original source, frame and memory representation.
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }
    /// Native six-component Gram stresses [sqrt(Pa)], in source branch order.
    pub fn mechanical_memory_sqrt_pa(&self) -> &[[f64; 6]] {
        &self.mechanical
    }
    /// Signed scalar polarizations [C/m²], in source pole order.
    pub fn polarization_c_m2(&self) -> &[f64] {
        &self.dielectric
    }
    /// Scalar count charged to the consuming spatial owner's memory budget.
    pub fn internal_variables(&self) -> usize {
        6 * self.mechanical.len() + self.dielectric.len()
    }
}

/// Additional real tensors for the same midpoint field increment. They contain
/// no equilibrium/base stiffness or instantaneous/base permittivity.
#[derive(Clone, Copy, Debug)]
pub struct MemoryMidpointCoefficients {
    /// Additional constant-field Mandel stiffness [Pa] acting on delta_S.
    pub stiffness_increment_mandel_pa: [[f64; 6]; 6],
    /// Additional fully clamped permittivity [F/m] acting on delta_E.
    pub permittivity_increment_f_m: [[f64; 3]; 3],
}

/// History terms when the accepted field already contains its physical memory.
#[derive(Clone, Copy, Debug)]
pub struct MemoryHistoryCorrection {
    /// Add +integral B_strain^T * this stress to the mechanical increment RHS.
    pub mechanical_rhs_stress_mandel_pa: [f64; 6],
    /// Delta_P independent of delta_E; the electric RHS gets -integral gradN.Delta_P.
    pub polarization_history_increment_c_m2: [f64; 3],
}

/// Physical additional fields and memory storage at one represented state.
/// The two named energies are subaccounts of memory_energy_j_m3.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MemoryResponse {
    /// Additional tensile Maxwell stress in the original Mandel order [Pa].
    pub additional_stress_mandel_pa: [f64; 6],
    /// Additional directed electric displacement in the current frame [C/m²].
    pub polarization_c_m2: [f64; 3],
    /// Mechanical memory free-energy density [J/m³].
    pub mechanical_memory_energy_j_m3: f64,
    /// Dielectric polarization free-energy density [J/m³].
    pub dielectric_memory_energy_j_m3: f64,
    /// Sum of the two memory energy densities, counted once [J/m³].
    pub memory_energy_j_m3: f64,
}

/// Local work/storage/heat account for one original coupled midpoint interval.
#[derive(Clone, Copy, Debug)]
pub struct PiezoelectricMemoryStep {
    /// Actual represented midpoint memories and their fields/storage.
    pub midpoint: MemoryResponse,
    /// Proposed endpoint memories and their fields/storage.
    pub endpoint: MemoryResponse,
    /// Original accepted memory free-energy density [J/m³].
    pub memory_energy_before_j_m3: f64,
    /// Proposed endpoint memory free-energy density [J/m³].
    pub memory_energy_after_j_m3: f64,
    /// Computed from state increments, not a subtraction of endpoint energies.
    pub memory_energy_change_j_m3: f64,
    /// Nonnegative mechanical relaxation heat released in the interval [J/m³].
    pub mechanical_heat_j_m3: f64,
    /// Nonnegative dielectric polarization heat released in the interval [J/m³].
    pub dielectric_heat_j_m3: f64,
    /// Integral of additional mechanical stress against actual strain change.
    pub mechanical_work_j_m3: f64,
    /// Midpoint field times actual directed polarization change.
    pub dielectric_work_j_m3: f64,
    /// Delta(memory energy) + both heat channels - both work channels.
    pub energy_residual_j_m3: f64,
}

/// Proposed history plus its independent account; no accepted input is mutated.
#[derive(Clone, Debug)]
pub struct PiezoelectricMemoryCandidate {
    memory: PiezoelectricMemory,
    step: PiezoelectricMemoryStep,
}

impl PiezoelectricMemoryCandidate {
    /// Proposed history, still separate from the accepted input.
    pub const fn memory(&self) -> &PiezoelectricMemory {
        &self.memory
    }
    /// Independent material work/storage/heat account.
    pub const fn step(&self) -> &PiezoelectricMemoryStep {
        &self.step
    }
    /// Publish only together with the enclosing accepted mechanical/electric state.
    pub fn into_parts(self) -> (PiezoelectricMemory, PiezoelectricMemoryStep) {
        (self.memory, self.step)
    }
}

fn half_step(dt: f64) -> Result<f64> {
    let half = 0.5 * dt;
    if !dt.is_finite() || dt <= 0.0 || half == 0.0 {
        return Err(invalid(
            "material midpoint requires a representable positive duration",
        ));
    }
    Ok(half)
}

// Neither tau+h nor h/tau needs to overflow. Preserve the small term instead
// of computing it by subtracting the large term from one.
fn weights(half: f64, tau: f64) -> Result<(f64, f64)> {
    let (a, b) = if half <= tau {
        let ratio = half / tau;
        if ratio == 0.0 {
            return Err(invalid(
                "material midpoint relaxation rate ratio underflows",
            ));
        }
        (1.0 / (1.0 + ratio), ratio / (1.0 + ratio))
    } else {
        let ratio = tau / half;
        if ratio == 0.0 {
            return Err(invalid(
                "material midpoint inverse relaxation rate ratio underflows",
            ));
        }
        (ratio / (1.0 + ratio), 1.0 / (1.0 + ratio))
    };
    Ok((a, b))
}

fn square_energy(value: f64) -> Result<f64> {
    scalar((0.5 * value) * value)
}

fn finish_response(response: &mut MemoryResponse) -> Result<()> {
    response.memory_energy_j_m3 =
        scalar(response.mechanical_memory_energy_j_m3 + response.dielectric_memory_energy_j_m3)?;
    if response
        .additional_stress_mandel_pa
        .iter()
        .chain(&response.polarization_c_m2)
        .any(|x| !x.is_finite())
    {
        return Err(invalid("material memory field is not representable"));
    }
    Ok(())
}

// Neumaier's fixed-order compensated addition preserves weak work/storage
// increments when different branches carry much larger amounts of energy.
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let value = scalar(value)?;
        let next = scalar(self.sum + value)?;
        self.correction = scalar(
            self.correction
                + if self.sum.abs() >= value.abs() {
                    (self.sum - next) + value
                } else {
                    (value - next) + self.sum
                },
        )?;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> Result<f64> {
        scalar(self.sum + self.correction)
    }
}

impl RelaxingTensorPiezoelectric {
    /// Scalar memory count before any state allocation. Native rank-deficient
    /// mechanical branches retain all six Gram coordinates without projection.
    pub fn internal_variables(&self) -> usize {
        6 * self.mechanical.branches().len() + self.dielectric.len()
    }

    fn admit_transient(&self) -> Result<()> {
        if self.frequency_band_rad_s() != [0.0, f64::MAX] {
            return Err(invalid(
                "frequency-limited relaxation sources cannot authorize unrestricted transients",
            ));
        }
        Ok(())
    }

    fn admit_memory(&self, memory: &PiezoelectricMemory) -> Result<()> {
        self.admit_transient()?;
        if memory.model_identity != self.identity()
            || memory.mechanical.len() != self.mechanical.branches().len()
            || memory.dielectric.len() != self.dielectric.len()
            || memory
                .mechanical
                .iter()
                .flatten()
                .chain(&memory.dielectric)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "material memory must retain its exact original model, sources and frame",
            ));
        }
        Ok(())
    }

    /// Zero memory is compatible with the original instantaneous dielectric
    /// field and relaxed mechanical stress. No charge/strain projection occurs.
    pub fn zero_memory(&self) -> Result<PiezoelectricMemory> {
        self.admit_transient()?;
        Ok(PiezoelectricMemory {
            model_identity: self.identity(),
            mechanical: vec![[0.0; 6]; self.mechanical.branches().len()],
            dielectric: vec![0.0; self.dielectric.len()],
        })
    }

    /// Allocation-free physical response of an accepted memory.
    pub fn memory_response(&self, memory: &PiezoelectricMemory) -> Result<MemoryResponse> {
        self.admit_memory(memory)?;
        let mut response = MemoryResponse::default();
        for ((branch, factor), state) in self
            .mechanical
            .branches()
            .iter()
            .zip(self.mechanical.factors())
            .zip(&memory.mechanical)
        {
            let gram = branch.resolved_gram_rows(factor);
            for row in 0..6 {
                response.mechanical_memory_energy_j_m3 =
                    scalar(response.mechanical_memory_energy_j_m3 + square_energy(state[row])?)?;
                for axis in 0..6 {
                    response.additional_stress_mandel_pa[axis] = gram[row][axis]
                        .mul_add(state[row], response.additional_stress_mandel_pa[axis]);
                }
            }
        }
        for (pole, &state) in self.dielectric.iter().zip(&memory.dielectric) {
            response.dielectric_memory_energy_j_m3 = scalar(
                response.dielectric_memory_energy_j_m3
                    + square_energy(state / det::sqrt(pole.permittivity_increment_f_m))?,
            )?;
            for axis in 0..3 {
                response.polarization_c_m2[axis] =
                    pole.direction[axis].mul_add(state, response.polarization_c_m2[axis]);
            }
        }
        finish_response(&mut response)?;
        Ok(response)
    }

    /// Read-only additional fields and storage of the same accepted memory.
    pub fn response(&self, memory: &PiezoelectricMemory) -> Result<MemoryResponse> {
        self.memory_response(memory)
    }

    /// Exact midpoint companion tensors, with no source query or allocation.
    pub fn midpoint_coefficients(&self, duration_s: f64) -> Result<MemoryMidpointCoefficients> {
        self.admit_transient()?;
        let half = half_step(duration_s)?;
        let mut stiffness = [[0.0; 6]; 6];
        let mut permittivity = [[0.0; 3]; 3];
        for branch in self.mechanical.branches() {
            let (a, _) = weights(half, branch.tau_s)?;
            for row in 0..6 {
                for column in row..6 {
                    let value =
                        scalar(a.mul_add(branch.c_pa[row][column], stiffness[row][column]))?;
                    stiffness[row][column] = value;
                    stiffness[column][row] = value;
                }
            }
        }
        for pole in self.dielectric.iter() {
            let (_, b) = weights(half, pole.relaxation_time_s)?;
            let coefficient = b * pole.permittivity_increment_f_m;
            for row in 0..3 {
                for column in row..3 {
                    let value = scalar(
                        (coefficient * pole.direction[row])
                            .mul_add(pole.direction[column], permittivity[row][column]),
                    )?;
                    permittivity[row][column] = value;
                    permittivity[column][row] = value;
                }
            }
        }
        Ok(MemoryMidpointCoefficients {
            stiffness_increment_mandel_pa: stiffness,
            permittivity_increment_f_m: permittivity,
        })
    }

    /// History corrections for an increment around an accepted field that
    /// already includes its memory stress and polarization. The signs here
    /// must not be reused for a field whose old memory was omitted.
    pub fn history_correction(
        &self,
        memory: &PiezoelectricMemory,
        duration_s: f64,
        electric_field_before_v_m: [f64; 3],
    ) -> Result<MemoryHistoryCorrection> {
        self.admit_memory(memory)?;
        let half = half_step(duration_s)?;
        if electric_field_before_v_m.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "memory history correction needs a finite original electric field",
            ));
        }
        let mut stress = [0.0; 6];
        let mut polarization = [0.0; 3];
        for ((branch, factor), state) in self
            .mechanical
            .branches()
            .iter()
            .zip(self.mechanical.factors())
            .zip(&memory.mechanical)
        {
            let (_, b) = weights(half, branch.tau_s)?;
            let gram = branch.resolved_gram_rows(factor);
            for row in 0..6 {
                for axis in 0..6 {
                    stress[axis] = gram[row][axis].mul_add(b * state[row], stress[axis]);
                }
            }
        }
        for (pole, &state) in self.dielectric.iter().zip(&memory.dielectric) {
            let (_, b) = weights(half, pole.relaxation_time_s)?;
            let along = pole
                .direction
                .into_iter()
                .zip(electric_field_before_v_m)
                .fold(0.0, |sum, (n, e)| n.mul_add(e, sum));
            let increment = (b * pole.permittivity_increment_f_m).mul_add(along, -b * state);
            for axis in 0..3 {
                polarization[axis] = pole.direction[axis].mul_add(increment, polarization[axis]);
            }
        }
        if stress
            .iter()
            .chain(&polarization)
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "material memory history correction is not representable",
            ));
        }
        Ok(MemoryHistoryCorrection {
            mechanical_rhs_stress_mandel_pa: stress,
            polarization_history_increment_c_m2: polarization,
        })
    }

    /// Advance material memories inside the caller's ONE coupled midpoint
    /// solve. delta_strain_mid is HALF the actual endpoint strain increment.
    /// The new vectors are the only runtime allocations; all observations and
    /// coefficients use fixed-size fields and the original native factors.
    pub fn midpoint_candidate(
        &self,
        memory: &PiezoelectricMemory,
        duration_s: f64,
        delta_strain_mid: [f64; 6],
        electric_field_mid_v_m: [f64; 3],
    ) -> Result<PiezoelectricMemoryCandidate> {
        self.admit_memory(memory)?;
        let half = half_step(duration_s)?;
        if delta_strain_mid
            .iter()
            .chain(&electric_field_mid_v_m)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "material midpoint needs finite strain increments and electric field",
            ));
        }
        let mut mechanical = Vec::with_capacity(memory.mechanical.len());
        let mut dielectric = Vec::with_capacity(memory.dielectric.len());
        let mut midpoint = MemoryResponse::default();
        let mut endpoint = MemoryResponse::default();
        let mut before = Sum::default();
        let mut change = Sum::default();
        let mut mechanical_heat = Sum::default();
        let mut dielectric_heat = Sum::default();
        let mut mechanical_work = Sum::default();
        let mut dielectric_work = Sum::default();
        let root_dt = det::sqrt(duration_s);
        for ((branch, factor), state) in self
            .mechanical
            .branches()
            .iter()
            .zip(self.mechanical.factors())
            .zip(&memory.mechanical)
        {
            let (a, b) = weights(half, branch.tau_s)?;
            let gram = branch.resolved_gram_rows(factor);
            let heat_scale = root_dt / det::sqrt(branch.tau_s);
            let mut next = [0.0; 6];
            for row in 0..6 {
                let drive = gram[row]
                    .into_iter()
                    .zip(delta_strain_mid)
                    .fold(0.0, |sum, (coefficient, strain)| {
                        coefficient.mul_add(strain, sum)
                    });
                // The increment form retains weak evolution under a large bias.
                let increment = scalar(2.0 * a.mul_add(drive, -b * state[row]))?;
                next[row] = scalar(state[row] + increment)?;
                let actual_delta = scalar(next[row] - state[row])?;
                let mid = f64::midpoint(state[row], next[row]);
                let offset = 0.5 * (state[row] - mid) + 0.5 * (next[row] - mid);
                before.add(square_energy(state[row])?)?;
                change.add(actual_delta * mid)?;
                change.add(actual_delta * offset)?;
                mechanical_work.add((2.0 * drive) * mid)?;
                let loss_amplitude = heat_scale * mid;
                mechanical_heat.add(loss_amplitude * loss_amplitude)?;
                midpoint.mechanical_memory_energy_j_m3 =
                    scalar(midpoint.mechanical_memory_energy_j_m3 + square_energy(mid)?)?;
                endpoint.mechanical_memory_energy_j_m3 =
                    scalar(endpoint.mechanical_memory_energy_j_m3 + square_energy(next[row])?)?;
                for axis in 0..6 {
                    midpoint.additional_stress_mandel_pa[axis] =
                        gram[row][axis].mul_add(mid, midpoint.additional_stress_mandel_pa[axis]);
                    endpoint.additional_stress_mandel_pa[axis] = gram[row][axis]
                        .mul_add(next[row], endpoint.additional_stress_mandel_pa[axis]);
                }
            }
            mechanical.push(next);
        }
        for (pole, &state) in self.dielectric.iter().zip(&memory.dielectric) {
            let (_, b) = weights(half, pole.relaxation_time_s)?;
            let root_epsilon = det::sqrt(pole.permittivity_increment_f_m);
            let along = pole
                .direction
                .into_iter()
                .zip(electric_field_mid_v_m)
                .fold(0.0, |sum, (n, e)| n.mul_add(e, sum));
            let increment =
                scalar(2.0 * (b * pole.permittivity_increment_f_m).mul_add(along, -b * state))?;
            let next = scalar(state + increment)?;
            let actual_delta = scalar(next - state)?;
            let mid = f64::midpoint(state, next);
            let before_amplitude = scalar(state / root_epsilon)?;
            let after_amplitude = scalar(next / root_epsilon)?;
            let mid_amplitude = scalar(mid / root_epsilon)?;
            let delta_amplitude = scalar(actual_delta / root_epsilon)?;
            let offset = (0.5 * (state - mid) + 0.5 * (next - mid)) / root_epsilon;
            before.add(square_energy(before_amplitude)?)?;
            change.add(delta_amplitude * mid_amplitude)?;
            change.add(delta_amplitude * offset)?;
            dielectric_work.add(along * actual_delta)?;
            // Actual represented p_mid, independent of the solved p increment.
            // If a positive physical loss is hidden by rounded state evolution,
            // it remains visible in the independent energy residual.
            let friction = root_epsilon.mul_add(along, -mid_amplitude);
            let loss_amplitude = (root_dt / det::sqrt(pole.relaxation_time_s)) * friction;
            dielectric_heat.add(loss_amplitude * loss_amplitude)?;
            midpoint.dielectric_memory_energy_j_m3 =
                scalar(midpoint.dielectric_memory_energy_j_m3 + square_energy(mid_amplitude)?)?;
            endpoint.dielectric_memory_energy_j_m3 =
                scalar(endpoint.dielectric_memory_energy_j_m3 + square_energy(after_amplitude)?)?;
            for axis in 0..3 {
                midpoint.polarization_c_m2[axis] =
                    pole.direction[axis].mul_add(mid, midpoint.polarization_c_m2[axis]);
                endpoint.polarization_c_m2[axis] =
                    pole.direction[axis].mul_add(next, endpoint.polarization_c_m2[axis]);
            }
            dielectric.push(next);
        }
        finish_response(&mut midpoint)?;
        finish_response(&mut endpoint)?;
        let mechanical_heat_j_m3 = mechanical_heat.value()?;
        let dielectric_heat_j_m3 = dielectric_heat.value()?;
        let mechanical_work_j_m3 = mechanical_work.value()?;
        let dielectric_work_j_m3 = dielectric_work.value()?;
        let memory_energy_change_j_m3 = change.value()?;
        let mut balance = Sum::default();
        for value in [
            memory_energy_change_j_m3,
            mechanical_heat_j_m3,
            dielectric_heat_j_m3,
            -mechanical_work_j_m3,
            -dielectric_work_j_m3,
        ] {
            balance.add(value)?;
        }
        let step = PiezoelectricMemoryStep {
            midpoint,
            endpoint,
            memory_energy_before_j_m3: before.value()?,
            memory_energy_after_j_m3: endpoint.memory_energy_j_m3,
            memory_energy_change_j_m3,
            mechanical_heat_j_m3,
            dielectric_heat_j_m3,
            mechanical_work_j_m3,
            dielectric_work_j_m3,
            energy_residual_j_m3: balance.value()?,
        };
        Ok(PiezoelectricMemoryCandidate {
            memory: PiezoelectricMemory {
                model_identity: self.identity(),
                mechanical,
                dielectric,
            },
            step,
        })
    }
}
