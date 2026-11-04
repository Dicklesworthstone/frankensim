//! Causal locally reacting walls in the original full acoustic trace space.
//!
//! The original [`WallPin`] supplies surface mass `sigma` [kg/m²], stiffness
//! `K` [Pa/m] and resistance `r` [Pa s/m]. Positive displacement `w` and
//! velocity `v` point OUT of the fluid: `sigma v_dot + r v + K w = p`.
//! Every active original boundary coordinate retains both states. With the
//! original face mass `W = integral(phi_i phi_j)`, wall energy is
//! `(sigma v^T W v + K w^T W w)/2`, and fluid feedback is `+W v`.
//! No patch-average piston, diagonal mass or extra area factor is substituted.
//!
//! One midpoint elimination gives
//! `v_mid = (p_mid + 2 sigma v_0/h - K w_0)/(2 sigma/h + r + h K/2)`.
//! The positive effective admittance enters the ORIGINAL sparse acoustic
//! factor; the history enters its original projected load. The same solve
//! reconstructs all wall states, then checks the uncondensed equations and
//! physical work. This effective admittance is never reported as wall heat.
//!
//! Patches have constant coefficients and disjoint original faces, including
//! disjointness from real Robin patches. Adjacent mechanical patches retain
//! independent states at their shared edges. Imposed velocity may share a
//! wall face. This is the declared local wall law, not shell bending, lateral
//! elastic waves, a hereditary material fit, or a moving fluid mesh.

use super::{
    BoundaryAdmittance, C64, Coo, Csr, Dd, HeterogeneousError, HeterogeneousHelmholtzFe,
    HeterogeneousTimeEnergy, HeterogeneousTimeLimits, HeterogeneousTimeState,
    PreparedHeterogeneousTime, RealVelocityTrace, ResistiveBoundary, energy_change, invalid, poll,
    positive_form, relative,
};
use fs_phs::WallPin;
use std::collections::BTreeSet;

/// Original surface-specific mechanical law on a complete boundary patch.
#[derive(Clone, Debug)]
pub struct MechanicalBoundary {
    pub faces: Vec<[u32; 3]>,
    pub wall: WallPin,
}

/// Existing acoustic admission plus an explicit full wall-state extent.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTimeLimits {
    pub acoustic: HeterogeneousTimeLimits,
    /// Sum of active boundary DOFs across patches; each has displacement and
    /// velocity. Patch CSR and O(ndof) work vectors are additional to factor bytes.
    pub max_wall_dofs: usize,
}

impl Default for MechanicalTimeLimits {
    fn default() -> Self {
        Self {
            acoustic: HeterogeneousTimeLimits::default(),
            max_wall_dofs: 4096,
        }
    }
}

struct PreparedWall {
    boundary: MechanicalBoundary,
    dofs: Vec<usize>,
    mass: Csr,
    admittance: f64,
    inertia_weight: f64,
    spring_weight: f64,
}

/// One original acoustic factor with every declared distributed wall state.
pub struct PreparedMechanicalTime<'m, 'f, 'c> {
    condensed: PreparedHeterogeneousTime<'m, 'f, 'c>,
    conductance: Csr,
    walls: Vec<PreparedWall>,
}

/// Candidate fields bound to this exact preparation, including its walls.
#[derive(Clone)]
pub struct MechanicalTimeState<'p, 'm, 'f, 'c> {
    owner: &'p PreparedMechanicalTime<'m, 'f, 'c>,
    acoustic: HeterogeneousTimeState<'m, 'f, 'c>,
    displacements: Vec<Vec<f64>>,
    velocities: Vec<Vec<f64>>,
}

impl core::fmt::Debug for MechanicalTimeState<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MechanicalTimeState")
            .field("acoustic", &self.acoustic)
            .field("wall_patches", &self.displacements.len())
            .finish_non_exhaustive()
    }
}

impl MechanicalTimeState<'_, '_, '_, '_> {
    /// Original hierarchical pressure-potential coefficients [Pa s].
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        self.acoustic.pressure_potential_pa_s()
    }
    /// Original hierarchical pressure coefficients [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        self.acoustic.pressure_pa()
    }
    /// Displacement [m] in the patch's original `wall_dofs` order.
    #[must_use]
    pub fn wall_displacement_m(&self, wall: usize) -> Option<&[f64]> {
        self.displacements.get(wall).map(Vec::as_slice)
    }
    /// Outward normal velocity [m/s] in the same original trace order.
    #[must_use]
    pub fn wall_normal_velocity_m_s(&self, wall: usize) -> Option<&[f64]> {
        self.velocities.get(wall).map(Vec::as_slice)
    }
}

/// Independently evaluated physical acoustic and wall storage [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicalTimeEnergy {
    pub acoustic: HeterogeneousTimeEnergy,
    pub wall_kinetic_j: f64,
    pub wall_elastic_j: f64,
}

impl MechanicalTimeEnergy {
    #[must_use]
    pub fn total_j(self) -> f64 {
        self.acoustic.total_j() + self.wall_kinetic_j + self.wall_elastic_j
    }
}

/// Original source work, reversible exchange and irreversible losses [J].
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub acoustic_delta_energy_j: f64,
    pub wall_delta_energy_j: f64,
    pub input_work_j: f64,
    /// Dissipation of the original REAL Robin patches only.
    pub boundary_dissipation_j: f64,
    /// `h sum r v_mid^T W v_mid`, never effective-admittance loss.
    pub wall_dissipation_j: f64,
    /// Signed fluid-to-wall exchange `h sum p_mid^T W v_mid`.
    pub wall_work_j: f64,
    /// Original uncondensed fluid and wall equations, not the stored factor.
    pub backward_error: f64,
}

impl MechanicalTimeWork {
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.boundary_dissipation_j + self.wall_dissipation_j
            - self.input_work_j
    }
    #[must_use]
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_delta_energy_j + self.boundary_dissipation_j + self.wall_work_j
            - self.input_work_j
    }
    #[must_use]
    pub fn wall_balance_residual_j(self) -> f64 {
        self.wall_delta_energy_j + self.wall_dissipation_j - self.wall_work_j
    }
}

/// Fully checked candidate; the caller owns acceptance and the clock.
#[derive(Clone, Debug)]
pub struct MechanicalTimeStep<'p, 'm, 'f, 'c> {
    pub state: MechanicalTimeState<'p, 'm, 'f, 'c>,
    pub work: MechanicalTimeWork,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare reciprocal full-trace mechanical walls on the original time owner.
    /// All wall/Robin faces must be disjoint. No material data is inferred.
    pub fn prepare_time_with_walls<'m>(
        &'m self,
        duration_s: f64,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        limits: MechanicalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 || limits.max_wall_dofs == 0 {
            return Err(
                invalid("positive mechanical time step and wall DOF limit required").into(),
            );
        }
        if self.ndof() > limits.acoustic.max_dofs {
            return Err(HeterogeneousError::Limit {
                what: "transient pressure degrees of freedom",
                requested: self.ndof(),
                maximum: limits.acoustic.max_dofs,
            });
        }
        let mut seen = BTreeSet::new();
        for faces in resistive
            .iter()
            .map(|b| &b.faces)
            .chain(boundaries.iter().map(|b| &b.faces))
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid("mechanical and real Robin faces must be disjoint").into());
                }
            }
        }
        let real = resistive
            .iter()
            .map(|b| BoundaryAdmittance {
                faces: b.faces.clone(),
                admittance_m_pa_s: C64::from_re(b.admittance_m_pa_s),
            })
            .collect::<Vec<_>>();
        let (conductance, _) = self.admittance_matrices(&real, cancelled)?;
        let mut total = 0usize;
        let mut walls = Vec::with_capacity(boundaries.len());
        let mut effective = resistive.to_vec();
        for boundary in boundaries {
            poll(cancelled)?;
            // Reuse the original WallPin admission, not a fitted frequency sample.
            fs_phs::wall_specific_impedance(&boundary.wall, 1.0)
                .map_err(HeterogeneousError::Wall)?;
            let pin = boundary.wall;
            let inertia = 2.0 * (pin.surface_density / duration_s);
            let impedance = inertia + pin.resistance + (0.5 * duration_s) * pin.stiffness_per_area;
            let admittance = 1.0 / impedance;
            let inertia_weight = inertia * admittance;
            let spring_weight = pin.stiffness_per_area * admittance;
            if ![impedance, admittance, inertia_weight, spring_weight]
                .iter()
                .all(|v| v.is_finite() && *v > 0.0)
            {
                return Err(
                    invalid("mechanical wall midpoint coefficients are unrepresentable").into(),
                );
            }
            // This is the SAME symmetric original face quadrature as Robin loading.
            let (global, _) = self.admittance_matrices(
                &[BoundaryAdmittance {
                    faces: boundary.faces.clone(),
                    admittance_m_pa_s: C64::ONE,
                }],
                cancelled,
            )?;
            let dofs = (0..self.ndof())
                .filter(|&row| {
                    let (columns, values) = global.row(row);
                    columns.binary_search(&row).is_ok_and(|i| values[i] > 0.0)
                })
                .collect::<Vec<_>>();
            total = total
                .checked_add(dofs.len())
                .ok_or_else(|| invalid("mechanical trace extent overflow"))?;
            if total > limits.max_wall_dofs {
                return Err(HeterogeneousError::Limit {
                    what: "distributed wall trace degrees of freedom",
                    requested: total,
                    maximum: limits.max_wall_dofs,
                });
            }
            let mut local = Coo::new(dofs.len(), dofs.len());
            for (row, &global_row) in dofs.iter().enumerate() {
                poll(cancelled)?;
                let (cols, entries) = global.row(global_row);
                for (&col, &entry) in cols.iter().zip(entries) {
                    if entry != 0.0 {
                        let col = dofs
                            .binary_search(&col)
                            .map_err(|_| invalid("incomplete original wall trace support"))?;
                        local.push(row, col, entry);
                    }
                }
            }
            effective.push(ResistiveBoundary {
                faces: boundary.faces.clone(),
                admittance_m_pa_s: admittance,
            });
            walls.push(PreparedWall {
                boundary: boundary.clone(),
                dofs,
                mass: local.assemble(),
                admittance,
                inertia_weight,
                spring_weight,
            });
        }
        let condensed = self.prepare_time(duration_s, &effective, limits.acoustic, cancelled)?;
        Ok(PreparedMechanicalTime {
            condensed,
            conductance,
            walls,
        })
    }
}

impl<'m, 'f, 'c> PreparedMechanicalTime<'m, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.condensed.duration_s()
    }
    #[must_use]
    pub fn wall_count(&self) -> usize {
        self.walls.len()
    }
    /// Active original hierarchical trace coordinates, sorted by global DOF.
    #[must_use]
    pub fn wall_dofs(&self, wall: usize) -> Option<&[usize]> {
        self.walls.get(wall).map(|w| w.dofs.as_slice())
    }
    #[must_use]
    pub fn wall(&self, wall: usize) -> Option<&MechanicalBoundary> {
        self.walls.get(wall).map(|w| &w.boundary)
    }
    #[must_use]
    pub fn zero_state(&self) -> MechanicalTimeState<'_, 'm, 'f, 'c> {
        MechanicalTimeState {
            owner: self,
            acoustic: self.condensed.zero_state(),
            displacements: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
            velocities: self.walls.iter().map(|w| vec![0.0; w.dofs.len()]).collect(),
        }
    }
    /// Import the original acoustic field; every wall starts at rest.
    pub fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        let acoustic = self.condensed.state(potential, pressure)?;
        let mut state = self.zero_state();
        state.acoustic = acoustic;
        Ok(state)
    }
    /// Import every physical history in `wall_dofs` order, without normalization.
    pub fn state_with_walls(
        &self,
        potential: &[f64],
        pressure: &[f64],
        displacements_m: &[Vec<f64>],
        velocities_m_s: &[Vec<f64>],
    ) -> Result<MechanicalTimeState<'_, 'm, 'f, 'c>, HeterogeneousError> {
        if displacements_m.len() != self.walls.len()
            || velocities_m_s.len() != self.walls.len()
            || self.walls.iter().enumerate().any(|(index, wall)| {
                displacements_m[index].len() != wall.dofs.len()
                    || velocities_m_s[index].len() != wall.dofs.len()
                    || displacements_m[index]
                        .iter()
                        .chain(&velocities_m_s[index])
                        .any(|v| !v.is_finite())
            })
        {
            return Err(invalid(
                "finite complete original wall displacement/velocity fields required",
            )
            .into());
        }
        Ok(MechanicalTimeState {
            owner: self,
            acoustic: self.condensed.state(potential, pressure)?,
            displacements: displacements_m.to_vec(),
            velocities: velocities_m_s.to_vec(),
        })
    }
    fn validate_state(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self, state.owner) {
            return Err(
                invalid("mechanical state belongs to another prepared wall configuration").into(),
            );
        }
        Ok(())
    }
    pub fn energy(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
    ) -> Result<MechanicalTimeEnergy, HeterogeneousError> {
        self.validate_state(state)?;
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }
    fn energy_with_roundoff(
        &self,
        state: &MechanicalTimeState<'_, '_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(MechanicalTimeEnergy, f64, f64), HeterogeneousError> {
        let (acoustic, acoustic_roundoff) = self
            .condensed
            .energy_with_roundoff(&state.acoustic, cancelled)?;
        let mut energy = MechanicalTimeEnergy {
            acoustic,
            wall_kinetic_j: 0.0,
            wall_elastic_j: 0.0,
        };
        let mut roundoff = 0.0;
        for (index, wall) in self.walls.iter().enumerate() {
            let (kinetic, rk) = positive_form(&wall.mass, &state.velocities[index], cancelled)?;
            let (elastic, re) = positive_form(&wall.mass, &state.displacements[index], cancelled)?;
            let pin = wall.boundary.wall;
            energy.wall_kinetic_j += 0.5 * pin.surface_density * kinetic;
            energy.wall_elastic_j += 0.5 * pin.stiffness_per_area * elastic;
            roundoff += 0.5 * (pin.surface_density * rk + pin.stiffness_per_area * re);
        }
        if !energy.total_j().is_finite() || !roundoff.is_finite() {
            return Err(invalid("mechanical acoustic storage is unrepresentable").into());
        }
        Ok((energy, acoustic_roundoff, roundoff))
    }

    /// Solve fluid and every wall together; no participant is advanced explicitly.
    /// Sources are original inward normal velocities at the step midpoint.
    pub fn step<'p>(
        &'p self,
        initial: &MechanicalTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate_state(initial)?;
        let core = &self.condensed;
        let h = self.duration_s();
        let source = core.model.geometry.boundary_linear_load(inward_velocity)?;
        let mut load = source.clone();
        for (index, wall) in self.walls.iter().enumerate() {
            poll(cancelled)?;
            let history = initial.velocities[index]
                .iter()
                .zip(&initial.displacements[index])
                .map(|(&v, &w)| wall.inertia_weight.mul_add(v, -wall.spring_weight * w))
                .collect::<Vec<_>>();
            let mut projected = vec![0.0; wall.dofs.len()];
            wall.mass.spmv(&history, &mut projected);
            for (&dof, value) in wall.dofs.iter().zip(projected) {
                load[dof] -= value;
            }
        }
        // Only this existing owner factors/solves. Its condensed Robin work is
        // an algebraic check, not an observable physical dissipation account.
        let acoustic = core.step_projected(&initial.acoustic, &load, cancelled)?;
        let midpoint_p = initial
            .pressure_pa()
            .iter()
            .zip(acoustic.state.pressure_pa())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let midpoint_psi = initial
            .pressure_potential_pa_s()
            .iter()
            .zip(acoustic.state.pressure_potential_pa_s())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let mut next = MechanicalTimeState {
            owner: self,
            acoustic: acoustic.state,
            displacements: Vec::with_capacity(self.walls.len()),
            velocities: Vec::with_capacity(self.walls.len()),
        };
        let mut wall_load = vec![0.0; source.len()];
        let mut wall_change = 0.0;
        let mut wall_loss = 0.0;
        let mut wall_loss_roundoff = 0.0;
        let mut wall_work = 0.0;
        let mut wall_work_roundoff = 0.0;
        let mut backward_error: f64 = 0.0;
        for (index, wall) in self.walls.iter().enumerate() {
            poll(cancelled)?;
            let pin = wall.boundary.wall;
            let previous_w = &initial.displacements[index];
            let previous_v = &initial.velocities[index];
            let mut displacement = Vec::with_capacity(wall.dofs.len());
            let mut velocity = Vec::with_capacity(wall.dofs.len());
            let mut middle_v = Vec::with_capacity(wall.dofs.len());
            let mut pressure = Vec::with_capacity(wall.dofs.len());
            for (i, &dof) in wall.dofs.iter().enumerate() {
                let p = midpoint_p[dof];
                let history = wall
                    .inertia_weight
                    .mul_add(previous_v[i], -wall.spring_weight * previous_w[i]);
                let v = wall.admittance.mul_add(p, history);
                let end_v = v + (v - previous_v[i]);
                let end_w = h.mul_add(v, previous_w[i]);
                if !end_v.is_finite() || !end_w.is_finite() {
                    return Err(invalid("mechanical wall candidate overflow").into());
                }
                let mid_v = f64::midpoint(previous_v[i], end_v);
                let mid_w = f64::midpoint(previous_w[i], end_w);
                let momentum = Dd::from_f64(pin.surface_density)
                    * (Dd::from_f64(end_v) - Dd::from_f64(previous_v[i]));
                let drag = Dd::from_f64(h) * Dd::from_f64(pin.resistance) * Dd::from_f64(mid_v);
                let spring =
                    Dd::from_f64(h) * Dd::from_f64(pin.stiffness_per_area) * Dd::from_f64(mid_w);
                let force = Dd::from_f64(h) * Dd::from_f64(p);
                let defect = (momentum + drag + spring - force).to_f64();
                let scale = [momentum, drag, spring, force]
                    .iter()
                    .map(|v| v.to_f64().abs())
                    .sum::<f64>();
                let kinematics = (Dd::from_f64(end_w)
                    - Dd::from_f64(previous_w[i])
                    - Dd::from_f64(h) * Dd::from_f64(mid_v))
                .to_f64();
                let kinematic_scale = end_w.abs() + previous_w[i].abs() + h * mid_v.abs();
                if !defect.is_finite()
                    || !scale.is_finite()
                    || !kinematics.is_finite()
                    || !kinematic_scale.is_finite()
                {
                    return Err(invalid("mechanical wall physical residual overflow").into());
                }
                backward_error = backward_error
                    .max(relative(defect.abs(), scale))
                    .max(relative(kinematics.abs(), kinematic_scale));
                displacement.push(end_w);
                velocity.push(end_v);
                middle_v.push(mid_v);
                pressure.push(p);
            }
            let mut projected = vec![0.0; wall.dofs.len()];
            wall.mass.spmv(&middle_v, &mut projected);
            for (&dof, &value) in wall.dofs.iter().zip(&projected) {
                wall_load[dof] += value;
            }
            let (work, rw) = paired_work(h, &pressure, &projected);
            wall_work += work;
            wall_work_roundoff += rw;
            let (loss, rl) = positive_form(&wall.mass, &middle_v, cancelled)?;
            wall_loss += h * pin.resistance * loss;
            wall_loss_roundoff += h * pin.resistance * rl;
            wall_change += pin.surface_density
                * energy_change(&wall.mass, previous_v, &velocity, cancelled)?
                + pin.stiffness_per_area
                    * energy_change(&wall.mass, previous_w, &displacement, cancelled)?;
            next.displacements.push(displacement);
            next.velocities.push(velocity);
        }
        let (_, fluid_error) = core.physical_residual(
            &initial.acoustic,
            &next.acoustic,
            &midpoint_p,
            &midpoint_psi,
            &self.conductance,
            &source,
            Some(&wall_load),
            cancelled,
        )?;
        backward_error = backward_error.max(fluid_error);
        if backward_error > core.limits.residual_tolerance {
            return Err(HeterogeneousError::Residual {
                observed: backward_error,
                limit: core.limits.residual_tolerance,
            });
        }
        let (initial_energy, initial_fluid_roundoff, initial_wall_roundoff) =
            self.energy_with_roundoff(initial, cancelled)?;
        let (final_energy, final_fluid_roundoff, final_wall_roundoff) =
            self.energy_with_roundoff(&next, cancelled)?;
        let (boundary_loss, boundary_roundoff) =
            positive_form(&self.conductance, &midpoint_p, cancelled)?;
        let (input, input_roundoff) = paired_work(h, &midpoint_p, &source);
        let work = MechanicalTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: acoustic.work.delta_energy_j + wall_change,
            acoustic_delta_energy_j: acoustic.work.delta_energy_j,
            wall_delta_energy_j: wall_change,
            input_work_j: input,
            boundary_dissipation_j: h * boundary_loss,
            wall_dissipation_j: wall_loss,
            wall_work_j: wall_work,
            backward_error,
        };
        let fluid_roundoff =
            initial_fluid_roundoff + final_fluid_roundoff + h * boundary_roundoff + input_roundoff;
        let wall_roundoff = initial_wall_roundoff + final_wall_roundoff + wall_loss_roundoff;
        let tolerance = core.limits.energy_balance_tolerance;
        balance(
            work.acoustic_balance_residual_j(),
            work.acoustic_delta_energy_j.abs()
                + work.boundary_dissipation_j
                + wall_work.abs()
                + input.abs(),
            fluid_roundoff + wall_work_roundoff,
            tolerance,
        )?;
        balance(
            work.wall_balance_residual_j(),
            wall_change.abs() + wall_loss + wall_work.abs(),
            wall_roundoff + wall_work_roundoff,
            tolerance,
        )?;
        balance(
            work.balance_residual_j(),
            work.delta_energy_j.abs() + work.boundary_dissipation_j + wall_loss + input.abs(),
            fluid_roundoff + wall_roundoff,
            tolerance,
        )?;
        poll(cancelled)?;
        Ok(MechanicalTimeStep { state: next, work })
    }
}

fn paired_work(h: f64, pressure: &[f64], load: &[f64]) -> (f64, f64) {
    let mut result = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (&p, &b) in pressure.iter().zip(load) {
        let term = Dd::from_f64(h) * Dd::from_f64(p) * Dd::from_f64(b);
        result = result + term;
        absolute += term.to_f64().abs();
    }
    (result.to_f64(), 64.0 * f64::EPSILON * absolute)
}

fn balance(
    residual: f64,
    scale: f64,
    roundoff: f64,
    tolerance: f64,
) -> Result<(), HeterogeneousError> {
    let limit = tolerance * scale + roundoff;
    if !residual.is_finite() || !limit.is_finite() {
        return Err(invalid("mechanical acoustic work is unrepresentable").into());
    }
    if residual.abs() > limit {
        return Err(HeterogeneousError::EnergyBalance {
            observed_j: residual,
            limit_j: limit,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
