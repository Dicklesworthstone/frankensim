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
use fs_phs::linear_exact::exponential;

pub(super) struct ExactModes {
    modes: ModeMaps,
    pressure: Vec<f64>,
    contact: Vec<f64>,
    arms: Vec<Vec<ExactArm>>,
}

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

    fn from_maps(
        spec: &ModalApertureSpec,
        modes: ModeMaps,
        arms: Vec<Vec<ExactArm>>,
    ) -> Result<Self, AcousticRealizeError> {
        let dt = spec.time_step_s;
        let n = spec.modes.len();
        let mut pressure = vec![0.0; n];
        let mut contact = vec![0.0; n];
        match &modes {
            ModeMaps::Diagonal(blocks) => {
                for (i, block) in blocks.iter().enumerate() {
                    let zero_memory = vec![0.0; arms.get(i).map_or(0, Vec::len)];
                    let (_, mobility) =
                        block.free_and_mobility(BladeState::default(), 0.0, dt, &zero_memory);
                    pressure[i] = mobility * spec.modes[i].pressure_area_m2;
                    contact[i] = mobility * spec.modes[i].gap_projection;
                }
            }
            ModeMaps::Coupled(block) => {
                for i in 0..n {
                    for j in 0..n {
                        let mobility = block.map[(3 * n + i) * 5 * n + 2 * n + j]
                            / block.root_stiffness[i]
                            / block.root_stiffness[j]
                            / dt;
                        pressure[i] += mobility * spec.modes[j].pressure_area_m2;
                        contact[i] += mobility * spec.modes[j].gap_projection;
                    }
                }
            }
        }
        if pressure.iter().chain(&contact).any(|v| !v.is_finite()) {
            return Err(invalid(
                "exact modal pressure/contact response is unrepresentable",
            ));
        }
        Ok(Self {
            modes,
            pressure,
            contact,
            arms,
        })
    }

    pub(super) fn pressure_response(&self) -> &[f64] {
        &self.pressure
    }
    pub(super) fn contact_response(&self) -> &[f64] {
        &self.contact
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
                let z = block.scaled(state, external);
                let n = state.len();
                for (i, value) in out.iter_mut().enumerate() {
                    *value = block.row(3 * n + i, &z) / block.root_stiffness[i] / dt;
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
                let z = block.scaled(old, held);
                let n = old.len();
                for i in 0..n {
                    state.push(BladeState {
                        displacement_m: old[i].displacement_m + dt * mean_velocity[i],
                        velocity_m_s: block.row(n + i, &z) / block.root_mass[i],
                    });
                    mean_displacement.push(
                        block.row(4 * n + i, &z) / block.omega[i] / block.root_stiffness[i] / dt,
                    );
                }
                for i in 0..z.len() {
                    for j in 0..z.len() {
                        viscous_loss_j += z[i] * block.viscous[i * z.len() + j] * z[j];
                    }
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

/// z = [sqrt(K) q, sqrt(M) v, F/sqrt(K)] has energy-scaled components.
/// Two extra rows per mode integrate displacement increments and mean position
/// without cancellation. This is the original material pencil, not new modes.
struct CoupledModes {
    root_mass: Vec<f64>,
    root_stiffness: Vec<f64>,
    omega: Vec<f64>,
    map: Vec<f64>,
    viscous: Vec<f64>,
}

impl CoupledModes {
    fn new(spec: &ModalApertureSpec, c: &[f64]) -> Result<Self, AcousticRealizeError> {
        let n = spec.modes.len();
        let (core, side) = (3 * n, 5 * n);
        let dt = spec.time_step_s;
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
            generator[(3 * n + i) * side + n + i] = w;
            generator[(4 * n + i) * side + i] = w;
            for j in 0..n {
                let damping = c[i * n + j] / root_mass[i] / root_mass[j];
                generator[(n + i) * side + n + j] = -dt * damping;
                weights[(n + i) * core + n + j] = damping;
            }
        }
        let map = exponential(&generator, side).map_err(super::nonlinear)?;
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
        if map.iter().chain(&viscous).any(|v| !v.is_finite()) {
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
        })
    }

    fn scaled(&self, state: &[BladeState], held: &[f64]) -> Vec<f64> {
        let n = state.len();
        let mut z = vec![0.0; 3 * n];
        for i in 0..n {
            z[i] = self.root_stiffness[i] * state[i].displacement_m;
            z[n + i] = self.root_mass[i] * state[i].velocity_m_s;
            z[2 * n + i] = held[i] / self.root_stiffness[i];
        }
        z
    }

    fn row(&self, row: usize, z: &[f64]) -> f64 {
        let side = 5 * self.root_mass.len();
        z.iter()
            .enumerate()
            .map(|(j, v)| self.map[row * side + j] * v)
            .sum()
    }
}
