//! Exact held-force structural response for the existing scalar slit solve.
//!
//! Diagonal modes use the same ExactBlade as the one-/two-blade owners. Full
//! material C retains every entry in one energy-scaled M/C/K exponential. The
//! pressure/contact solver consumes displacement INCREMENTS over dt, so swept
//! volume and mechanical work remain conjugate. End velocity is not inferred
//! from a midpoint identity. No contact potential or channel law lives here.

use super::{AcousticRealizeError, BladeState, ModalApertureSpec, invalid};
use crate::bernoulli_aperture::paired::{ExactArm, ExactBlade};
use fs_math::det;
use fs_phs::linear_exact::{MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, exponential};

/// One positive material-memory channel in physical modal coordinates.
/// Its storage is (coupling^T q - z)^2 / 2; z has units sqrt(J).
#[derive(Clone, Debug)]
pub(super) struct CoupledArm {
    pub coupling: Vec<f64>,
    pub relaxation_time_s: f64,
}

#[derive(Clone)]
pub(super) struct ExactModes {
    modes: ModeMaps,
    pressure: Vec<f64>,
    contact: Vec<f64>,
    force_response: Vec<f64>,
    arms: Vec<Vec<ExactArm>>,
}

#[derive(Clone)]
enum ModeMaps {
    Diagonal(Vec<ExactBlade>),
    Coupled(CoupledModes),
}

pub(super) struct ExactEndpoint {
    pub state: Vec<BladeState>,
    pub mean_displacement: Vec<f64>,
    pub viscous_loss_j: f64,
    pub memory: Vec<Vec<f64>>,
    pub relaxation_loss_j: f64,
    pub relaxation_force_n: Vec<f64>,
}

impl ExactModes {
    pub(super) fn new(
        spec: &ModalApertureSpec,
        viscosity: Option<&[f64]>,
    ) -> Result<Self, AcousticRealizeError> {
        let dt = spec.time_step_s;
        let modes = match viscosity {
            Some(c) if c.iter().any(|x| *x != 0.0) => {
                ModeMaps::Coupled(CoupledModes::new(spec, c)?)
            }
            _ => ModeMaps::Diagonal(
                spec.modes
                    .iter()
                    .map(|m| {
                        ExactBlade::from_mechanics(m.mass_kg, m.damping_n_s_m, m.stiffness_n_m, dt)
                    })
                    .collect::<Result<_, _>>()?,
            ),
        };
        Self::from_maps(spec, modes, Vec::new())
    }

    pub(super) fn with_arms(
        spec: &ModalApertureSpec,
        arms: Vec<Vec<ExactArm>>,
    ) -> Result<Self, AcousticRealizeError> {
        if arms.len() != spec.modes.len() {
            return Err(invalid("modal relaxation must cover every retained mode"));
        }
        let modes = ModeMaps::Diagonal(
            spec.modes
                .iter()
                .zip(&arms)
                .map(|(m, arms)| {
                    ExactBlade::with_arms(
                        m.mass_kg,
                        m.damping_n_s_m,
                        m.stiffness_n_m,
                        arms,
                        spec.time_step_s,
                    )
                })
                .collect::<Result<_, _>>()?,
        );
        Self::from_maps(spec, modes, arms)
    }

    /// Retain complete regional Maxwell forms instead of diagonalizing away
    /// their cross-mode reactions. All channels share one accepted history row.
    pub(super) fn with_coupled_arms(
        spec: &ModalApertureSpec,
        arms: Vec<CoupledArm>,
    ) -> Result<Self, AcousticRealizeError> {
        let n = spec.modes.len();
        CoupledModes::dimensions(n, arms.len())?;
        let mut c = vec![0.0; n * n];
        for (i, mode) in spec.modes.iter().enumerate() {
            c[i * n + i] = mode.damping_n_s_m;
        }
        Self::with_coupled_memory(spec, &c, arms)
    }

    /// Full Kelvin matrix and full Maxwell memory share one physical block.
    /// C is the complete viscosity, including any admitted diagonal entries.
    pub(super) fn with_coupled_memory(
        spec: &ModalApertureSpec,
        c: &[f64],
        arms: Vec<CoupledArm>,
    ) -> Result<Self, AcousticRealizeError> {
        let block = CoupledModes::with_memory(spec, c, Some(arms))?;
        Self::from_maps(spec, ModeMaps::Coupled(block), Vec::new())
    }

    fn from_maps(
        spec: &ModalApertureSpec,
        modes: ModeMaps,
        arms: Vec<Vec<ExactArm>>,
    ) -> Result<Self, AcousticRealizeError> {
        let dt = spec.time_step_s;
        let n = spec.modes.len();
        let mut pressure = vec![0.0; n];
        let mut contact = vec![0.0; n];
        let mut force_response = vec![0.0; n * n];
        match &modes {
            ModeMaps::Diagonal(blocks) => {
                for (i, block) in blocks.iter().enumerate() {
                    let zero_memory = vec![0.0; arms.get(i).map_or(0, Vec::len)];
                    let (_, mobility) =
                        block.free_and_mobility(BladeState::default(), 0.0, dt, &zero_memory);
                    force_response[n * i + i] = mobility;
                    pressure[i] = mobility * spec.modes[i].pressure_area_m2;
                    contact[i] = mobility * spec.modes[i].gap_projection;
                }
            }
            ModeMaps::Coupled(block) => {
                for i in 0..n {
                    for j in 0..n {
                        let mobility = block.map[(block.core + i) * block.side + 2 * n + j]
                            / block.root_stiffness[i]
                            / block.root_stiffness[j]
                            / dt;
                        force_response[n * i + j] = mobility;
                        pressure[i] += mobility * spec.modes[j].pressure_area_m2;
                        contact[i] += mobility * spec.modes[j].gap_projection;
                    }
                }
            }
        }
        if pressure.iter().chain(&contact).chain(&force_response).any(|v| !v.is_finite()) {
            return Err(invalid(
                "exact modal pressure/contact response is unrepresentable",
            ));
        }
        Ok(Self {
            modes,
            pressure,
            contact,
            force_response,
            arms,
        })
    }

    pub(super) fn pressure_response(&self) -> &[f64] {
        &self.pressure
    }

    /// Change only the work-conjugate acoustic geometry. The original complete
    /// mechanical propagator, viscosity, Maxwell arms and force map survive.
    pub(super) fn rebind_pressure(&mut self, spec: &ModalApertureSpec)
        -> Result<(), AcousticRealizeError>
    {
        let n = spec.modes.len();
        let pressure: Vec<f64> = self.force_response.chunks_exact(n).map(|row| {
            row.iter().zip(&spec.modes).map(|(b, mode)| b * mode.pressure_area_m2).sum()
        }).collect();
        if pressure.len() != n || pressure.iter().any(|p| !p.is_finite()) {
            return Err(invalid("wetted geometry gives an unrepresentable exact pressure response"));
        }
        self.pressure = pressure;
        Ok(())
    }
    pub(super) fn contact_response(&self) -> &[f64] {
        &self.contact
    }

    // Every spatial contact column uses the SAME cached full structural map,
    // including off-diagonal viscosity and regional material-memory reactions.
    pub(super) fn force_velocity(&self, force: &[f64], out: &mut [f64]) {
        for (row, value) in self.force_response.chunks_exact(force.len()).zip(out) {
            *value = row.iter().zip(force).map(|(b, f)| b * f).sum();
        }
    }

    pub(super) fn free_velocity(
        &self,
        state: &[BladeState],
        external: &[f64],
        dt: f64,
        memory: Option<&[Vec<f64>]>,
        out: &mut [f64],
    ) -> Result<(), AcousticRealizeError> {
        match &self.modes {
            ModeMaps::Diagonal(blocks) => {
                for (i, block) in blocks.iter().enumerate() {
                    let memory = memory.map_or(&[][..], |states| states[i].as_slice());
                    out[i] = block.free_and_mobility(state[i], external[i], dt, memory).0;
                }
            }
            ModeMaps::Coupled(block) => {
                block.validate_memory(memory)?;
                let z = block.scaled(state, external, memory);
                for (i, value) in out.iter_mut().enumerate() {
                    *value = block.row(block.core + i, &z) / block.root_stiffness[i] / dt;
                }
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("exact modal free response is unrepresentable"));
        }
        Ok(())
    }

    pub(super) fn endpoint(
        &self,
        old: &[BladeState],
        held: &[f64],
        mean_velocity: &[f64],
        dt: f64,
        memory: Option<&[Vec<f64>]>,
    ) -> ExactEndpoint {
        let mut state = Vec::with_capacity(old.len());
        let mut mean_displacement = Vec::with_capacity(old.len());
        let mut viscous_loss_j = 0.0;
        let mut next_memory = Vec::new();
        let mut relaxation_loss_j = 0.0;
        let mut relaxation_force_n = vec![0.0; old.len()];
        match &self.modes {
            ModeMaps::Diagonal(blocks) => {
                for (i, block) in blocks.iter().enumerate() {
                    let memory = memory.map_or(&[][..], |states| states[i].as_slice());
                    state.push(BladeState {
                        displacement_m: old[i].displacement_m + dt * mean_velocity[i],
                        velocity_m_s: block.end_velocity(old[i], held[i], memory),
                    });
                    mean_displacement.push(block.mean_displacement(old[i], held[i], dt, memory));
                    viscous_loss_j += block.viscous_loss_j(old[i], held[i], memory);
                    if !self.arms.is_empty() {
                        let next = block.end_memory(old[i], held[i], memory);
                        relaxation_loss_j += block.arm_loss_j(old[i], held[i], memory);
                        // dz/dt=(c*q-z)/tau: the independently integrated arm
                        // force is c*tau*Delta(z)/dt, never an energy remainder.
                        relaxation_force_n[i] = self.arms[i]
                            .iter()
                            .zip(memory)
                            .zip(&next)
                            .map(|((arm, old), next)| {
                                arm.coupling * arm.relaxation_time_s * (next - old) / dt
                            })
                            .sum();
                        next_memory.push(next);
                    }
                }
            }
            ModeMaps::Coupled(block) => {
                let z = block.scaled(old, held, memory);
                let n = old.len();
                for i in 0..n {
                    state.push(BladeState {
                        displacement_m: old[i].displacement_m + dt * mean_velocity[i],
                        velocity_m_s: block.row(n + i, &z) / block.root_mass[i],
                    });
                    mean_displacement.push(
                        block.row(block.core + n + i, &z)
                            / block.omega[i]
                            / block.root_stiffness[i]
                            / dt,
                    );
                }
                for i in 0..z.len() {
                    for j in 0..z.len() {
                        viscous_loss_j += z[i] * block.viscous[i * z.len() + j] * z[j];
                        relaxation_loss_j += z[i] * block.relaxation[i * z.len() + j] * z[j];
                    }
                }
                if let Some(arms) = &block.arms {
                    let next: Vec<_> = (0..arms.len()).map(|a| block.row(3 * n + a, &z)).collect();
                    for (a, arm) in arms.iter().enumerate() {
                        // Integral(c^T q-z) dt = tau * Delta(z). This is a
                        // constitutive-force integral, not an energy remainder.
                        let integral = arm.relaxation_time_s * (next[a] - z[3 * n + a]) / dt;
                        for (force, coupling) in relaxation_force_n.iter_mut().zip(&arm.coupling) {
                            *force += coupling * integral;
                        }
                    }
                    next_memory.push(next);
                }
            }
        }
        ExactEndpoint {
            state,
            mean_displacement,
            viscous_loss_j,
            memory: next_memory,
            relaxation_loss_j,
            relaxation_force_n,
        }
    }
}

/// z = [sqrt(K) q, sqrt(M) v, F/sqrt(K), memory] has energy-scaled components.
/// Two extra rows per mode integrate displacement increments and mean position
/// without cancellation. This is the original material pencil, not new modes.
#[derive(Clone)]
struct CoupledModes {
    root_mass: Vec<f64>,
    root_stiffness: Vec<f64>,
    omega: Vec<f64>,
    map: Vec<f64>,
    viscous: Vec<f64>,
    relaxation: Vec<f64>,
    arms: Option<Vec<CoupledArm>>,
    core: usize,
    side: usize,
}

impl CoupledModes {
    fn new(spec: &ModalApertureSpec, c: &[f64]) -> Result<Self, AcousticRealizeError> {
        Self::with_memory(spec, c, None)
    }

    fn dimensions(n: usize, arms: usize) -> Result<(usize, usize), AcousticRealizeError> {
        let core = n.checked_mul(3).and_then(|v| v.checked_add(arms));
        let side = n.checked_mul(5).and_then(|v| v.checked_add(arms));
        match (core, side) {
            (Some(core), Some(side))
                if n > 0
                    && side <= MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
                    && core <= MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 =>
            {
                Ok((core, side))
            }
            _ => Err(invalid(
                "complete coupled modal memory exceeds the exact exponential capacity",
            )),
        }
    }

    fn with_memory(
        spec: &ModalApertureSpec,
        c: &[f64],
        arms: Option<Vec<CoupledArm>>,
    ) -> Result<Self, AcousticRealizeError> {
        let n = spec.modes.len();
        let (core, side) = Self::dimensions(n, arms.as_ref().map_or(0, Vec::len))?;
        let dt = spec.time_step_s;
        if !dt.is_finite()
            || dt <= 0.0
            || arms.as_ref().is_some_and(|arms| {
                arms.iter().any(|a| {
                    a.coupling.len() != n
                        || a.coupling.iter().any(|v| !v.is_finite())
                        || !a.relaxation_time_s.is_finite()
                        || a.relaxation_time_s <= 0.0
                })
            })
        {
            return Err(invalid(
                "coupled material arms need complete finite rows and positive times",
            ));
        }
        let root_mass: Vec<_> = spec.modes.iter().map(|m| det::sqrt(m.mass_kg)).collect();
        let root_stiffness: Vec<_> = spec
            .modes
            .iter()
            .map(|m| det::sqrt(m.stiffness_n_m))
            .collect();
        let omega: Vec<_> = root_stiffness
            .iter()
            .zip(&root_mass)
            .map(|(k, m)| k / m)
            .collect();
        if c.len() != n * n || omega.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(invalid("exact modal material pencil is unrepresentable"));
        }
        let mut generator = vec![0.0; side * side];
        let mut weights = vec![0.0; core * core];
        for i in 0..n {
            let w = dt * omega[i];
            generator[i * side + n + i] = w;
            generator[(n + i) * side + i] = -w;
            generator[(n + i) * side + 2 * n + i] = w;
            generator[(core + i) * side + n + i] = w;
            generator[(core + n + i) * side + i] = w;
            for j in 0..n {
                let damping = c[i * n + j] / root_mass[i] / root_mass[j];
                generator[(n + i) * side + n + j] = -dt * damping;
                weights[(n + i) * core + n + j] = damping;
            }
        }
        let mut memory_weights = vec![0.0; core * core];
        if let Some(arms) = &arms {
            for (a, arm) in arms.iter().enumerate() {
                let at = 3 * n + a;
                let rate = 1.0 / arm.relaxation_time_s;
                generator[at * side + at] = -dt * rate;
                memory_weights[at * core + at] = rate;
                for i in 0..n {
                    let strain_i = arm.coupling[i] / root_stiffness[i];
                    let force_i = arm.coupling[i] / root_mass[i];
                    generator[(n + i) * side + at] = dt * force_i;
                    generator[at * side + i] = dt * rate * strain_i;
                    memory_weights[i * core + at] = -rate * strain_i;
                    memory_weights[at * core + i] = -rate * strain_i;
                    for j in 0..n {
                        let strain_j = arm.coupling[j] / root_stiffness[j];
                        generator[(n + i) * side + j] -= dt * force_i * strain_j;
                        memory_weights[i * core + j] += rate * strain_i * strain_j;
                    }
                }
            }
        }
        let map = exponential(&generator, side).map_err(super::nonlinear)?;
        let viscous = integrated_loss(&generator, side, core, dt, &weights)?;
        let relaxation = integrated_loss(&generator, side, core, dt, &memory_weights)?;
        if map
            .iter()
            .chain(&viscous)
            .chain(&relaxation)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "exact modal transition or material loss is unrepresentable",
            ));
        }
        Ok(Self {
            root_mass,
            root_stiffness,
            omega,
            map,
            viscous,
            relaxation,
            arms,
            core,
            side,
        })
    }

    fn validate_memory(&self, memory: Option<&[Vec<f64>]>) -> Result<(), AcousticRealizeError> {
        if let Some(arms) = &self.arms {
            if !memory.is_some_and(|rows| {
                rows.len() == 1
                    && rows[0].len() == arms.len()
                    && rows[0].iter().all(|v| v.is_finite())
            }) {
                return Err(invalid(
                    "coupled material memory requires one complete finite history row",
                ));
            }
        }
        Ok(())
    }

    fn scaled(&self, state: &[BladeState], held: &[f64], memory: Option<&[Vec<f64>]>) -> Vec<f64> {
        let n = state.len();
        let mut z = vec![0.0; self.core];
        for i in 0..n {
            z[i] = self.root_stiffness[i] * state[i].displacement_m;
            z[n + i] = self.root_mass[i] * state[i].velocity_m_s;
            z[2 * n + i] = held[i] / self.root_stiffness[i];
        }
        if self.arms.is_some() {
            for (a, value) in z[3 * n..].iter_mut().enumerate() {
                *value = memory
                    .and_then(|rows| rows.first())
                    .and_then(|row| row.get(a))
                    .copied()
                    .unwrap_or(f64::NAN);
            }
        }
        z
    }

    fn row(&self, row: usize, z: &[f64]) -> f64 {
        z.iter()
            .enumerate()
            .map(|(j, v)| self.map[row * self.side + j] * v)
            .sum()
    }
}

fn integrated_loss(
    generator: &[f64],
    side: usize,
    core: usize,
    dt: f64,
    weights: &[f64],
) -> Result<Vec<f64>, AcousticRealizeError> {
    if weights.iter().all(|v| *v == 0.0) {
        return Ok(vec![0.0; core * core]);
    }
    // Independent integral of v^T C v on the same held-force trajectory.
    // Short-interval Van Loan and Gramian doubling avoid unstable inverse
    // exponentials when the original material has very rapid viscous decay.
    let mut scales = 0;
    let mut norm = 0.0_f64;
    for row in generator[..core * side].chunks_exact(side) {
        norm = norm.max(row[..core].iter().map(|v| v.abs()).sum());
    }
    while norm > 0.5 {
        norm *= 0.5;
        scales += 1;
        if scales > 48 {
            return Err(invalid(
                "exact modal loss integral exceeds its bounded scaling work",
            ));
        }
    }
    let fraction = 2.0_f64.powi(-(scales as i32));
    let block = 2 * core;
    let mut van_loan = vec![0.0; block * block];
    for i in 0..core {
        for j in 0..core {
            van_loan[i * block + j] = -fraction * generator[j * side + i];
            van_loan[(core + i) * block + core + j] = fraction * generator[i * side + j];
            van_loan[i * block + core + j] = fraction * dt * weights[i * core + j];
        }
    }
    let exp = exponential(&van_loan, block).map_err(super::nonlinear)?;
    let mut transition = vec![0.0; core * core];
    let mut viscous = vec![0.0; core * core];
    for i in 0..core {
        for j in 0..core {
            transition[i * core + j] = exp[(core + i) * block + core + j];
            viscous[i * core + j] = (0..core)
                .map(|k| exp[(core + k) * block + core + i] * exp[k * block + core + j])
                .sum();
        }
    }
    for _ in 0..scales {
        let mut product = vec![0.0; core * core];
        for i in 0..core {
            for k in 0..core {
                for j in 0..core {
                    product[i * core + j] += viscous[i * core + k] * transition[k * core + j];
                }
            }
        }
        let mut next = vec![0.0; core * core];
        for i in 0..core {
            for j in 0..core {
                for k in 0..core {
                    viscous[i * core + j] += transition[k * core + i] * product[k * core + j];
                    next[i * core + j] += transition[i * core + k] * transition[k * core + j];
                }
            }
        }
        transition = next;
    }
    if viscous.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "exact modal transition or material loss is unrepresentable",
        ));
    }
    Ok(viscous)
}

#[cfg(test)]
mod tests;
