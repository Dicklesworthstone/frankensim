//! Continuous forced response of the complete already selected causal image.
//!
//! Peak phasors use exp(-i omega t). The same original shell, material and load
//! coordinates as the free-mode analysis enter one shared physical LU. Bare
//! component impedances are never inverted. Quadratic observations here are
//! cycle means; the separate pole API keeps its complex-envelope convention.

use super::{
    PreparedShellResonance, ShellResonanceError, cap, finite, memories, poll, projection, quadratic,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::harmonic::{HarmonicFactor, HarmonicResponse};

/// Additional complete physical observations, beyond the prepared harmonic
/// owner's existing factor/response budgets. No source coordinate is truncated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellForcedBudget {
    pub maximum_payload_bytes: usize,
    pub maximum_observation_work: u64,
}

/// Forced load amplitudes and independent cycle-mean storage and work.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellForcedLoadResponse {
    pub state: Vec<C64>,
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    /// Re(conj(flow_p)*effort_p)/2, positive into the original load.
    pub absorbed_port_power_w: Vec<f64>,
    /// Re(z^H Q z)/4, using the load's declared physical storage.
    pub mean_stored_energy_j: f64,
    /// Re((Qz)^H R (Qz))/2, independently evaluated original load heat.
    pub dissipated_power_w: f64,
    /// Load heat minus absorbed port power, independently evaluated.
    pub power_residual_w: f64,
    pub relative_equation_residual: f64,
}

/// Complete source-backed response with actual external and internal work.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellForcedResponse {
    /// Original joint state, output, mean storage and total physical heat.
    /// State order matches PreparedShellResonance, including every load state.
    pub physical: HarmonicResponse,
    /// Original force inputs [N], followed by signed exterior and cavity
    /// pressure inputs [Pa]. Pressure sign follows the authored surface maps.
    pub input: Vec<C64>,
    pub displacement_m_sqrt_kg: Vec<C64>,
    pub velocity_m_sqrt_kg_per_s: Vec<C64>,
    pub force_port_displacement_m: Vec<C64>,
    pub mean_mechanical_storage_j: f64,
    pub background_dissipation_w: f64,
    pub interface_dissipation_w: f64,
    /// Original selected substrate-row order; these are actual branch channels.
    pub substrate_branch_mean_storage_j: Vec<f64>,
    pub substrate_branch_dissipation_w: Vec<f64>,
    /// Original selected coating-row order, not inferred per-film heat shares.
    pub coating_branch_mean_storage_j: Vec<f64>,
    pub coating_branch_dissipation_w: Vec<f64>,
    pub load: Option<ShellForcedLoadResponse>,
    /// Shell heat + signed power into the load - external supplied power.
    pub shell_power_residual_w: f64,
    /// Maximum independent source/load energy, heat and work discrepancy.
    pub physical_relative_residual: f64,
    /// Actual source mechanical, kinematic, material and load equation checks.
    pub original_equation_relative_residual: f64,
}

impl ShellForcedResponse {
    pub fn force_port_velocity_m_s(&self) -> &[C64] {
        &self.physical.output[..self.force_port_displacement_m.len()]
    }
    pub fn exterior_volume_velocity_m3_s(&self) -> C64 {
        self.physical.output[self.force_port_displacement_m.len()]
    }
    pub fn cavity_volume_velocity_m3_s(&self) -> C64 {
        self.physical.output[self.force_port_displacement_m.len() + 1]
    }
}

/// One full-state physical factorization, reusable across all external drives.
pub struct ShellForcedFactor<'model, 'source> {
    model: &'model PreparedShellResonance<'source>,
    factor: HarmonicFactor<'model>,
    omega_rad_s: f64,
}

impl<'source> PreparedShellResonance<'source> {
    /// Original force ports, then exterior and cavity pressure ports. Internal
    /// attachment ports do not become independent external excitation slots.
    pub fn physical_port_count(&self) -> usize {
        self.image.force_ports.len() + 2
    }

    /// Factor the joint continuous response within the common authored source
    /// band. All source memory and independent load states remain in the solve,
    /// so a bare load resonance/antiresonance needs no special inversion path.
    ///
    /// This evaluates the already selected causal Prony/substrate continuation;
    /// it does not replace it with a frozen or newly fitted source law. The
    /// original fixed shell basis and source-law approximation limits remain.
    /// This is a linear image response. Finite-amplitude coating strain
    /// admission remains the source image's separate physical strain check.
    /// No accepted time history, sample clock or free-mode normalization is used.
    pub fn harmonic_factor(
        &self,
        omega_rad_s: f64,
        budget: ShellForcedBudget,
        gate: &CancelGate,
    ) -> Result<ShellForcedFactor<'_, 'source>, ShellResonanceError> {
        poll(gate)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || omega_rad_s < self.band.0
            || omega_rad_s > self.band.1
        {
            return Err(ShellResonanceError::Invalid(
                "positive forced frequency within the common shell/load source band required",
            ));
        }
        let (bytes, work) = self.observation_budget(1);
        cap(
            "forced physical observation bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "forced physical observation work",
            work,
            budget.maximum_observation_work as u128,
        )?;
        let factor = self
            .model
            .factor(omega_rad_s, || gate.is_requested())
            .map_err(ShellResonanceError::Preparation)?;
        poll(gate)?;
        Ok(ShellForcedFactor {
            model: self,
            factor,
            omega_rad_s,
        })
    }
}

fn relative(defect: f64, scale: f64) -> f64 {
    defect.abs() / scale.max(f64::MIN_POSITIVE)
}

impl ShellForcedFactor<'_, '_> {
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// One finite peak input per original external port. Multiple drive vectors
    /// reuse this same factor. Refusal publishes no partial response and leaves
    /// the factor, source images and accepted histories unchanged.
    pub fn response(
        &self,
        input: &[C64],
        gate: &CancelGate,
    ) -> Result<ShellForcedResponse, ShellResonanceError> {
        poll(gate)?;
        let model = self.model;
        if input.len() != model.physical_port_count() || input.iter().any(|&u| !finite(u)) {
            return Err(ShellResonanceError::Invalid(
                "one finite input per original shell port required",
            ));
        }
        let mut physical = self
            .factor
            .response(input, || gate.is_requested())
            .map_err(ShellResonanceError::Preparation)?;
        let m = model.mode_count();
        let s = C64::new(0.0, -self.omega_rad_s);
        let displacement: Vec<_> = (0..m).map(|i| physical.state[2 * i]).collect();
        let velocity: Vec<_> = (0..m).map(|i| physical.state[2 * i + 1]).collect();
        let mechanical: f64 = model
            .image
            .modes
            .iter()
            .enumerate()
            .map(|(i, body)| {
                0.25 * (velocity[i].norm_sq() + displacement[i].scale(body.omega).norm_sq())
            })
            .sum();
        let background = 0.5 * quadratic(model.image.background, &velocity, gate)?;
        let mut interface = 0.0;
        for &(row, c) in &model.image.interfaces {
            poll(gate)?;
            interface += 0.5 * c * projection(row, &velocity).norm_sq();
        }
        let mut stores = [Vec::new(), Vec::new()];
        let mut losses = [Vec::new(), Vec::new()];
        let mut at = 2 * m;
        for (family, (rows, times)) in memories(&model.image).into_iter().enumerate() {
            for (row, &tau) in rows.iter().zip(times) {
                poll(gate)?;
                let square = (projection(row, &displacement) - physical.state[at]).norm_sq();
                stores[family].push(0.25 * square);
                losses[family].push(0.5 * square / tau);
                at += 1;
            }
        }
        // The existing load observer checks its driven original state equation
        // for arbitrary s and supplied shell velocity, independently of poles.
        // Convert only its quadratic envelope forms to harmonic cycle means.
        let load = model
            .load
            .as_ref()
            .map(|load| {
                let raw = load.observe(&physical.state[at..], &velocity, s, gate)?;
                Ok::<_, ShellResonanceError>(ShellForcedLoadResponse {
                    state: raw.state,
                    flow: raw.flow,
                    effort: raw.effort,
                    absorbed_port_power_w: raw
                        .absorbed_port_power_w
                        .into_iter()
                        .map(|p| 0.5 * p)
                        .collect(),
                    mean_stored_energy_j: 0.5 * raw.stored_energy_j,
                    dissipated_power_w: 0.5 * raw.dissipated_power_w,
                    power_residual_w: 0.5 * raw.decay_power_residual_w,
                    relative_equation_residual: raw.relative_equation_residual,
                })
            })
            .transpose()?;
        let original = self.original_rows(
            input,
            &physical.state,
            &displacement,
            &velocity,
            load.as_ref(),
            gate,
        )?;
        physical.relative_residual = physical.relative_residual.max(original);
        let shell_energy = mechanical + stores.iter().flatten().sum::<f64>();
        let shell_heat = background + interface + losses.iter().flatten().sum::<f64>();
        let load_energy = load.as_ref().map_or(0.0, |l| l.mean_stored_energy_j);
        let load_heat = load.as_ref().map_or(0.0, |l| l.dissipated_power_w);
        let load_power = load
            .as_ref()
            .map_or(0.0, |l| l.absorbed_port_power_w.iter().sum::<f64>());
        let interface_scale = load.as_ref().map_or(0.0, |l| {
            l.flow
                .iter()
                .zip(&l.effort)
                .map(|(u, y)| 0.5 * u.abs() * y.abs())
                .sum::<f64>()
        });
        let supplied_scale: f64 = input
            .iter()
            .zip(&physical.output)
            .map(|(u, y)| 0.5 * u.abs() * y.abs())
            .sum();
        let energy = shell_energy + load_energy;
        let heat = shell_heat + load_heat;
        let shell_residual = shell_heat + load_power - physical.supplied_power_w;
        let power_scale = 2.0 * self.omega_rad_s * energy + heat.abs() + supplied_scale;
        let mut physical_relative = relative(
            energy - physical.mean_stored_energy_j,
            energy.abs() + physical.mean_stored_energy_j.abs(),
        )
        .max(relative(
            heat - physical.dissipated_power_w,
            heat.abs() + physical.dissipated_power_w.abs(),
        ))
        .max(relative(heat - physical.supplied_power_w, power_scale))
        .max(relative(
            shell_residual,
            2.0 * self.omega_rad_s * shell_energy
                + shell_heat.abs()
                + interface_scale
                + supplied_scale,
        ));
        if let Some(load) = &load {
            physical_relative = physical_relative
                .max(relative(
                    load.power_residual_w,
                    2.0 * self.omega_rad_s * load_energy + load_heat.abs() + interface_scale,
                ))
                .max(load.relative_equation_residual);
        }
        if [energy, heat, shell_residual, physical_relative]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(ShellResonanceError::Invalid(
                "forced shell physical accounting overflow",
            ));
        }
        if physical_relative > model.tolerance {
            return Err(ShellResonanceError::PhysicalAccounting {
                relative: physical_relative,
                tolerance: model.tolerance,
            });
        }
        let mut force_displacement = Vec::with_capacity(model.image.force_ports.len());
        for row in model.image.force_ports {
            poll(gate)?;
            force_displacement.push(projection(row, &displacement));
        }
        let [substrate_store, coating_store] = stores;
        let [substrate_loss, coating_loss] = losses;
        poll(gate)?;
        Ok(ShellForcedResponse {
            physical,
            input: input.to_vec(),
            displacement_m_sqrt_kg: displacement,
            velocity_m_sqrt_kg_per_s: velocity,
            force_port_displacement_m: force_displacement,
            mean_mechanical_storage_j: mechanical,
            background_dissipation_w: background,
            interface_dissipation_w: interface,
            substrate_branch_mean_storage_j: substrate_store,
            substrate_branch_dissipation_w: substrate_loss,
            coating_branch_mean_storage_j: coating_store,
            coating_branch_dissipation_w: coating_loss,
            load,
            shell_power_residual_w: shell_residual,
            physical_relative_residual: physical_relative,
            original_equation_relative_residual: original,
        })
    }

    fn original_rows(
        &self,
        input: &[C64],
        state: &[C64],
        displacement: &[C64],
        velocity: &[C64],
        load: Option<&ShellForcedLoadResponse>,
        gate: &CancelGate,
    ) -> Result<f64, ShellResonanceError> {
        let model = self.model;
        let m = model.mode_count();
        let s = C64::new(0.0, -self.omega_rad_s);
        // Project each original row once. Repeating these projections inside
        // every mechanical equation would introduce an unbudgeted cubic cost.
        let mut interface_velocity = Vec::with_capacity(model.image.interfaces.len());
        for &(row, _) in &model.image.interfaces {
            poll(gate)?;
            interface_velocity.push(projection(row, velocity));
        }
        let mut branch_drive = Vec::with_capacity(model.substrate_rows + model.coating_rows);
        for (rows, _) in memories(&model.image) {
            for row in rows {
                poll(gate)?;
                branch_drive.push(projection(row, displacement));
            }
        }
        let mut maximum = load.map_or(0.0, |l| l.relative_equation_residual);
        let mut check = |value: C64, scale: f64| -> Result<(), ShellResonanceError> {
            if !finite(value) || !scale.is_finite() {
                return Err(ShellResonanceError::Invalid(
                    "forced original shell equation overflow",
                ));
            }
            maximum = maximum.max(relative(value.abs(), scale));
            Ok(())
        };
        for i in 0..m {
            poll(gate)?;
            check(
                s * displacement[i] - velocity[i],
                s.abs() * displacement[i].abs() + velocity[i].abs(),
            )?;
            let mut value = s * velocity[i];
            let mut scale = value.abs();
            let mut add = |v: C64| {
                value = value + v;
                scale += v.abs();
            };
            add(displacement[i].scale(model.image.modes[i].omega.powi(2)));
            for j in 0..m {
                add(velocity[j].scale(model.image.background[i * m + j]));
            }
            for (&(row, c), &flow) in model.image.interfaces.iter().zip(&interface_velocity) {
                add(flow.scale(c * row[i]));
            }
            let mut at = 2 * m;
            for (rows, times) in memories(&model.image) {
                for (row, _) in rows.iter().zip(times) {
                    add((branch_drive[at - 2 * m] - state[at]).scale(row[i]));
                    at += 1;
                }
            }
            if let (Some(prepared), Some(load)) = (&model.load, load) {
                for (p, &effort) in load.effort.iter().enumerate() {
                    add(effort.scale(prepared.rows[p * m + i]));
                }
            }
            for (p, row) in model.image.force_ports.iter().enumerate() {
                add(input[p].scale(-row[i]));
            }
            let p = model.image.force_ports.len();
            add(input[p].scale(-model.image.exterior[i]));
            add(input[p + 1].scale(-model.image.cavity[i]));
            check(value, scale)?;
        }
        let mut at = 2 * m;
        for (rows, times) in memories(&model.image) {
            for (_, &tau) in rows.iter().zip(times) {
                poll(gate)?;
                let drive = branch_drive[at - 2 * m];
                let derivative = s.scale(tau) * state[at];
                check(
                    derivative + state[at] - drive,
                    derivative.abs() + state[at].abs() + drive.abs(),
                )?;
                at += 1;
            }
        }
        if maximum > model.tolerance {
            return Err(ShellResonanceError::PhysicalAccounting {
                relative: maximum,
                tolerance: model.tolerance,
            });
        }
        Ok(maximum)
    }
}
