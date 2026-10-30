//! Joint damaged-area friction on the existing finite native cohesive patch.
//! Original bodies, clock, local source identities and all damage histories
//! publish together only after every actual-motion and complete energy gate.
use super::super::frictional_damage::DamageFrictionFrame;
use super::*;
use crate::cohesive_contact::frictional_damage::{DamageFrictionRegime, DamageFrictionState};
use crate::cohesive_contact::patch::frictional_damage::{
    DamageFrictionPatchConfig, DamageFrictionPatchError, DamageFrictionPatchLoadLine,
    DamageFrictionPatchSite,
};
use crate::cohesive_contact::reclosure::{ReclosureActivation, ResidualCoulombSource};

pub(super) struct DamageFrictionPatchBinding<'a> {
    load_line: DamageFrictionPatchLoadLine<'a>,
    states: Vec<DamageFrictionState<'a>>,
    previous: Vec<DamageFrictionState<'a>>,
    staged: Vec<DamageFrictionState<'a>>,
}
impl<'a> CohesivePatchModalSystem<'a> {
    /// Select the same damaged-area rigid friction law at EVERY original site.
    /// Reference quadrature and body compliance are retained. Sources and actual
    /// receiving-material assignments are matched independently at each site.
    /// This changes no accepted point, body, clock or previous observation.
    pub fn enable_frictional_damage(
        &mut self,
        friction: &[&'a ResidualCoulombSource],
        receiving: &[[&MaterialStateId; 2]],
        activation: &[ReclosureActivation],
        mut config: DamageFrictionPatchConfig,
        gate: &CancelGate,
    ) -> Result<(), CohesivePatchModalError> {
        poll(Some(gate))?;
        let n = self.points.len();
        if self.damage_friction.is_some()
            || friction.len() != n
            || receiving.len() != n
            || activation.len() != n
        {
            return Err(invalid("unbound complete patch friction selection required").into());
        }
        if config.friction.contact.max_material_evaluations < 3 * n {
            return Err(
                invalid("patch friction budget must reserve every actual-motion audit").into(),
            );
        }
        config.friction.contact.max_material_evaluations -= n;
        let mut sites = Vec::with_capacity(n);
        for p in 0..n {
            sites.push(DamageFrictionPatchSite {
                cohesive: self.sites[p].source.ok_or_else(|| {
                    invalid("every frictional patch site needs its original cohesive source")
                })?,
                friction: friction[p],
                area_m2: self.sites[p].area_m2,
            });
        }
        let load_line =
            DamageFrictionPatchLoadLine::new(sites, self.compliance_m_per_n().to_vec(), config)?;
        let mut states = Vec::with_capacity(n);
        for p in 0..n {
            poll(Some(gate))?;
            states.push(load_line.adopt_point(p, self.points[p], receiving[p], activation[p])?);
        }
        let binding = DamageFrictionPatchBinding {
            previous: states.clone(),
            staged: states.clone(),
            states,
            load_line,
        };
        poll(Some(gate))?;
        self.damage_friction = Some(binding);
        Ok(())
    }
    pub fn damage_friction_source(&self, point: usize) -> Option<&'a ResidualCoulombSource> {
        self.damage_friction
            .as_ref()
            .and_then(|b| b.load_line.friction_source(point))
    }
    pub fn damage_friction_occupancy(&self, point: usize) -> Option<DamageFrictionRegime> {
        self.damage_friction
            .as_ref()
            .and_then(|b| b.states.get(point).map(|s| s.occupancy()))
    }
    pub(super) fn step_damage_friction(
        &mut self,
        external: &[f64],
        prescribed_rate: &[[f64; 3]],
        gate: Option<&CancelGate>,
    ) -> Result<&CohesivePatchModalFrame, CohesivePatchModalError> {
        poll(gate)?;
        if prescribed_rate.len() != self.points.len() {
            return Err(invalid("prescribed cohesive velocities must match all sites").into());
        }
        let before = self.total_energy_j()?;
        let dt = self.network.dt;
        let binding = self
            .damage_friction
            .as_mut()
            .ok_or_else(|| invalid("patch friction selection absent"))?;
        for p in 0..self.points.len() {
            if binding.states[p].point() != self.points[p] {
                return Err(
                    DamageFrictionPatchError::Invalid("canonical patch history mismatch").into(),
                );
            }
            binding.previous[p] = binding.states[p]
                .with_endpoint_tangential_velocity(tangential_velocity(
                    &self.sites[p],
                    self.network.components(),
                    prescribed_rate[p],
                )?)
                .map_err(|source| DamageFrictionPatchError::Point { point: p, source })?;
        }
        for (delta, rate) in self.prescribed_delta.iter_mut().zip(prescribed_rate) {
            for i in 0..3 {
                delta[i] = finite(rate[i] * dt)?;
            }
        }
        self.network.prepare_prediction(external, gate)?;
        for p in 0..self.points.len() {
            poll(gate)?;
            let old = self.points[p].jump_m();
            for i in 0..3 {
                self.free_jump[p][i] = finite(
                    old[i]
                        + self.prescribed_delta[p][i]
                        + dot(&self.sites[p].columns[i], &self.network.free_delta)?,
                )?;
            }
        }
        let solved = binding
            .load_line
            .trial(&binding.previous, &self.free_jump, dt, || {
                !gate.is_some_and(CancelGate::is_requested)
            })?;
        self.forces.copy_from_slice(external);
        for p in 0..self.points.len() {
            poll(gate)?;
            for i in 0..3 {
                for (force, &weight) in self.forces.iter_mut().zip(&self.sites[p].columns[i]) {
                    *force = finite(*force - weight * solved.mean_force_n[p][i])?;
                }
            }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut cohesive_energy = 0.0;
        let mut fracture = 0.0;
        let mut viscous = 0.0;
        let mut friction = 0.0;
        let mut cohesive_work = 0.0;
        let mut prescribed_work = 0.0;
        let mut failed_points = 0;
        for p in 0..self.points.len() {
            poll(gate)?;
            let old = self.points[p].jump_m();
            let mut actual_jump = [0.0; 3];
            let mut compatibility_residual = [0.0; 3];
            for i in 0..3 {
                actual_jump[i] = finite(
                    old[i]
                        + self.prescribed_delta[p][i]
                        + dot(&self.sites[p].columns[i], &self.network.free_delta)?,
                )?;
                compatibility_residual[i] = finite(actual_jump[i] - solved.jump_m[p][i])?;
                if compatibility_residual[i].abs() > solved.compatibility_tolerance_m[p] {
                    return Err(CohesivePatchModalError::JumpMismatch {
                        point: p,
                        component: i,
                        residual_m: compatibility_residual[i],
                        tolerance_m: solved.compatibility_tolerance_m[p],
                    });
                }
                limit(
                    "actual cohesive patch jump increment",
                    (actual_jump[i] - old[i]).abs(),
                    self.config.patch.contact.maximum_jump_increment_m,
                )?;
            }
            let audited = solved.audit_point(
                p,
                binding.previous[p],
                actual_jump,
                dt,
                solved.mean_force_n[p],
                solved.regimes[p],
            )?;
            binding.staged[p] = audited
                .state
                .with_endpoint_tangential_velocity(tangential_velocity(
                    &self.sites[p],
                    self.network.staged_components(),
                    prescribed_rate[p],
                )?)
                .map_err(|source| DamageFrictionPatchError::Point { point: p, source })?;
            let actual = audited.material;
            let area = self.sites[p].area_m2;
            let mut force_residual = [0.0; 3];
            let mut force_tolerance = [0.0; 3];
            for i in 0..3 {
                let expected = audited.expected_force_n[i];
                limit(
                    "actual cohesive patch mean force",
                    expected.abs(),
                    self.config.patch.contact.maximum_force_n,
                )?;
                force_residual[i] = finite(solved.mean_force_n[p][i] - expected)?;
                force_tolerance[i] = finite(
                    self.config.force_absolute_tolerance_n
                        + self.config.force_relative_tolerance
                            * expected.abs().max(solved.mean_force_n[p][i].abs()),
                )?
                .min(audited.force_tolerance_n[i]);
                if force_residual[i].abs() > force_tolerance[i] {
                    return Err(CohesivePatchModalError::ForceMismatch {
                        point: p,
                        component: i,
                        residual_n: force_residual[i],
                        tolerance_n: force_tolerance[i],
                    });
                }
            }
            let point_energy = finite(area * actual.equilibrium.point.stored_energy_j_m2())?;
            let point_fracture = finite(area * actual.equilibrium.energy.dissipated)?;
            let point_viscous = finite(area * actual.viscous_dissipated_j_m2)?;
            let point_work = finite(area * actual.total_work_j_m2())?;
            if point_fracture < 0.0 || point_viscous < 0.0 {
                return Err(invalid(
                    "cohesive patch fracture and viscosity must not create energy",
                )
                .into());
            }
            cohesive_energy = finite(cohesive_energy + point_energy)?;
            fracture = finite(fracture + point_fracture)?;
            viscous = finite(viscous + point_viscous)?;
            friction = finite(friction + audited.friction_dissipation_j)?;
            cohesive_work = finite(cohesive_work + point_work)?;
            prescribed_work =
                finite(prescribed_work + dot(&solved.mean_force_n[p], &self.prescribed_delta[p])?)?;
            failed_points += usize::from(actual.equilibrium.point.failed());
            self.staged_points[p] = actual.equilibrium.point;
            self.staged_observations[p] = CohesivePatchPointFrame {
                area_m2: area,
                jump_m: actual_jump,
                damage: actual.equilibrium.point.damage(),
                mean_force_n: solved.mean_force_n[p],
                endpoint_traction_pa: actual.traction_pa,
                cohesive_energy_j: point_energy,
                fracture_dissipation_j: point_fracture,
                viscous_dissipation_j: point_viscous,
                damage_friction: Some(DamageFrictionFrame {
                    regime: audited.regime,
                    damaged_compression_force_n: audited.damaged_compression_force_n,
                    force_n: audited.friction_force_n,
                    dissipation_j: audited.friction_dissipation_j,
                    work_defect_j: audited.work_residual_j,
                }),
                cohesive_work_j: point_work,
                compatibility_residual_m: compatibility_residual,
                compatibility_tolerance_m: solved.compatibility_tolerance_m[p],
                constitutive_residual_n: force_residual,
                force_tolerance_n: force_tolerance,
            };
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let after = finite(network_energy + cohesive_energy)?;
        limit(
            "total energy including finite cohesive patch",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let loss = finite(network_loss + fracture + viscous + friction)?;
        let residual = finite((after - before) + loss - external_work - prescribed_work)?;
        let scale = before
            .max(after)
            .max(loss.abs())
            .max(external_work.abs())
            .max(prescribed_work.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = CohesivePatchModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            cohesive_energy_j: cohesive_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            fracture_dissipation_j: fracture,
            viscous_dissipation_j: viscous,
            friction_dissipation_j: friction,
            cohesive_work_j: cohesive_work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            failed_points,
            iterations: solved.iterations,
            kinetic_switches: solved.kinetic_switches,
            material_evaluations: solved.material_evaluations + self.points.len(),
        };
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        self.network.publish_staged();
        core::mem::swap(&mut binding.states, &mut binding.staged);
        std::mem::swap(&mut self.points, &mut self.staged_points);
        std::mem::swap(&mut self.observations, &mut self.staged_observations);
        Ok(self.last.insert(accepted))
    }
}

fn tangential_velocity(
    site: &CohesiveModalSite<'_>,
    components: &[ModalAcousticTimeModel],
    prescribed: [f64; 3],
) -> Result<[f64; 2], ModalCouplingError> {
    let mut result = [prescribed[1], prescribed[2]];
    for (j, mode) in components.iter().flat_map(|m| m.states()).enumerate() {
        for i in 0..2 {
            result[i] = finite(result[i] + site.columns[i + 1][j] * mode.velocity_m_sqrt_kg_per_s)?;
        }
    }
    Ok(result)
}
