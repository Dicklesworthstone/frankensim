//! Intact adhesive relaxation on the original native modal network.
//!
//! The original ordered interface supplies thickness and independent normal /
//! shear spectra. Its Prony history follows actual accepted body motion through
//! all three work columns and the complete native compliance. One affine 3-by-3
//! load-line solve uses cached exact ramp coefficients; a second constitutive
//! evaluation audits the actual propagated coordinates before publication.
//!
//! The bond is massless, small-strain, and intact in a fixed local frame. A
//! straight jump ramp is the constitutive approximation; refine the mechanical
//! timestep for trajectory convergence. Physical relaxation heat and native
//! numerical loss remain separate. No damage, healing, thermal-expansion strain
//! or automatic heat destination is inferred.
use super::{
    CohesiveLoadLine, CohesiveModalConfig, CohesiveModalError, CoupledModalSystem,
    ModalAcousticTimeModel, ModalCouplingError, dot, finite, invalid, limit, network_response,
    poll,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_material::cohesive::relaxation::{
    ViscoelasticBondPoint, ViscoelasticBondRamp, source::ResolvedViscoelasticBond,
};
use fs_math::dd::Dd;

/// Accepted native motion and original-area adhesive work for one timestep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticBondModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub bond_energy_j: f64,
    /// Increment from the original constitutive history, not a subtraction of
    /// large cumulative storage or work inventories.
    pub bond_stored_change_j: f64,
    pub bond_work_j: f64,
    pub relaxation_heat_j: f64,
    pub network_dissipation_j: f64,
    pub numerical_dissipation_j: f64,
    pub external_work_j: f64,
    pub prescribed_work_j: f64,
    pub bond_work_on_network_j: f64,
    /// Small bond account, independent of unrelated native stored energy.
    pub bond_work_residual_j: f64,
    pub bond_work_tolerance_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub jump_m: [f64; 3],
    pub mean_force_n: [f64; 3],
    pub endpoint_traction_pa: [f64; 3],
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub constitutive_residual_n: [f64; 3],
    pub force_tolerance_n: [f64; 3],
}

/// One originally bonded patch retaining the complete native mechanical owner.
/// Setup allocates the response/force buffers and Prony coefficients. Timesteps
/// reuse those buffers and the existing small LU workspace; only bounded
/// constitutive successors are allocated, never copies of the native network.
pub struct ViscoelasticBondModalSystem<'a> {
    network: CoupledModalSystem,
    source: &'a ResolvedViscoelasticBond,
    point: ViscoelasticBondPoint<'a>,
    prepared: ViscoelasticBondRamp<'a>,
    columns: [Vec<f64>; 3],
    load_line: CohesiveLoadLine,
    config: CohesiveModalConfig,
    forces: Vec<f64>,
    last: Option<ViscoelasticBondModalFrame>,
}

impl<'a> ViscoelasticBondModalSystem<'a> {
    /// Explicitly bond the sample-zero configuration at zero jump and zero
    /// internal branch strain. Columns retain every native coordinate in
    /// component-then-mode order: jump = B^T delta_q, reaction = -B F.
    /// Surface identities are ordered, and area is original bonded area [m2].
    /// Initial body velocity is preserved; no preload equilibrium is invented.
    pub fn from_source(
        network: CoupledModalSystem,
        source: &'a ResolvedViscoelasticBond,
        surface_materials: [&MaterialStateId; 2],
        area_m2: f64,
        columns: [Vec<f64>; 3],
        config: CohesiveModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CohesiveModalError> {
        poll(Some(gate))?;
        if source.surfaces()[0].material != *surface_materials[0]
            || source.surfaces()[1].material != *surface_materials[1]
        {
            return Err(
                invalid("relaxing bond requires its original ordered surface materials").into(),
            );
        }
        if network.samples_rendered() != 0 {
            return Err(invalid("relaxing bond requires a sample-zero native network").into());
        }
        if columns.iter().any(|column| {
            column.len() != network.mode_count() || column.iter().any(|x| !x.is_finite())
        }) || !config.force_absolute_tolerance_n.is_finite()
            || config.force_absolute_tolerance_n <= 0.0
            || !config.force_relative_tolerance.is_finite()
            || !(0.0 < config.force_relative_tolerance && config.force_relative_tolerance < 1.0)
            || config.contact.max_material_evaluations < 2
        {
            return Err(invalid("complete bond work columns, positive force gates and two material evaluations required").into());
        }
        let mut compliance = [[0.0; 3]; 3];
        for j in 0..3 {
            let response = network_response(&network, &columns[j], gate)?;
            for i in 0..3 {
                compliance[i][j] = dot(&columns[i], &response)?;
            }
        }
        let mut solve_config = config.contact;
        solve_config.max_material_evaluations -= 1;
        let load_line = CohesiveLoadLine::new(area_m2, compliance, solve_config)?;
        let prepared = source.law().prepare(network.sample_period_s())?;
        let system = Self {
            forces: vec![0.0; network.mode_count()],
            point: source.law().bonded_point(),
            prepared,
            network,
            source,
            columns,
            load_line,
            config,
            last: None,
        };
        limit(
            "total initial energy including relaxing bond",
            system.total_energy_j()?,
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
    }

    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn point(&self) -> &ViscoelasticBondPoint<'a> {
        &self.point
    }
    pub const fn source(&self) -> &'a ResolvedViscoelasticBond {
        self.source
    }
    pub fn area_m2(&self) -> f64 {
        self.load_line.area_m2()
    }
    pub fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.load_line.compliance_m_per_n()
    }
    pub const fn work_columns(&self) -> &[Vec<f64>; 3] {
        &self.columns
    }
    pub const fn config(&self) -> CohesiveModalConfig {
        self.config
    }
    pub fn last_frame(&self) -> Option<&ViscoelasticBondModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, CohesiveModalError> {
        Ok(finite(
            self.network.total_energy_j()? + self.area_m2() * self.point.stored_energy_j_m2(),
        )?)
    }

    /// Constant authored modal forces and prescribed local jump rate [m/s].
    /// Use a zero prescribed rate for a wholly internal two-body bond.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: [f64; 3],
    ) -> Result<&ViscoelasticBondModalFrame, CohesiveModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, None)
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<&ViscoelasticBondModalFrame, CohesiveModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, Some(gate))
    }

    /// Admit a receiving work/heat account against staged original bodies.
    /// The callback returns a candidate; install it infallibly only after this
    /// method succeeds. Refusal and late cancellation publish no body state,
    /// Prony history, time or frame.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: [f64; 3],
        gate: &CancelGate,
        admit: impl FnOnce(&ViscoelasticBondModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<CohesiveModalError>,
    {
        let (point, frame) = self
            .prepare_step(external, prescribed_jump_rate_m_s, Some(gate))
            .map_err(E::from)?;
        let receiving = admit(&frame, self.network.staged_components())?;
        poll(Some(gate))
            .map_err(CohesiveModalError::from)
            .map_err(E::from)?;
        self.network.publish_staged();
        self.point = point;
        self.last = Some(frame);
        Ok(receiving)
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        rate: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<&ViscoelasticBondModalFrame, CohesiveModalError> {
        let (point, frame) = self.prepare_step(external, rate, gate)?;
        self.network.publish_staged();
        self.point = point;
        Ok(self.last.insert(frame))
    }

    fn prepare_step(
        &mut self,
        external: &[f64],
        rate: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<(ViscoelasticBondPoint<'a>, ViscoelasticBondModalFrame), CohesiveModalError> {
        poll(gate)?;
        let old = self.point.jump_m();
        let mut prescribed = [0.0; 3];
        for i in 0..3 {
            prescribed[i] = finite(rate[i] * self.network.dt)?;
        }
        self.network.prepare_prediction(external, gate)?;
        let mut free = [0.0; 3];
        for i in 0..3 {
            free[i] =
                finite(old[i] + prescribed[i] + dot(&self.columns[i], &self.network.free_delta)?)?;
        }
        let solved = self
            .load_line
            .trial_relaxation(&self.point, &self.prepared, free, || {
                !gate.is_some_and(CancelGate::is_requested)
            })?;
        let held = solved.mean_force_n;
        for (j, (force, input)) in self.forces.iter_mut().zip(external).enumerate() {
            if j % 256 == 0 {
                poll(gate)?;
            }
            let mut sum = Dd::from_f64(*input);
            for i in 0..3 {
                sum = sum - Dd::from_f64(self.columns[i][j]) * Dd::from_f64(held[i]);
            }
            *force = finite(sum.to_f64())?;
        }
        self.network.stage_inner(&self.forces, gate)?;
        // Audit propagated original coordinates, not merely the load-line
        // prediction or an assumed meaning of native propagation scratch.
        for (j, ((change, end), start)) in self
            .network
            .free_delta
            .iter_mut()
            .zip(self.network.candidates.iter().flat_map(|m| m.states()))
            .zip(self.network.models.iter().flat_map(|m| m.states()))
            .enumerate()
        {
            if j % 256 == 0 {
                poll(gate)?;
            }
            *change = finite(end.displacement_m_sqrt_kg - start.displacement_m_sqrt_kg)?;
        }
        let mut body = [0.0; 3];
        let mut travel = [0.0; 3];
        let mut jump = [0.0; 3];
        let mut compatibility = [0.0; 3];
        for i in 0..3 {
            body[i] = dot(&self.columns[i], &self.network.free_delta)?;
            travel[i] = finite(prescribed[i] + body[i])?;
            limit(
                "actual relaxing bond jump increment",
                travel[i].abs(),
                self.config.contact.maximum_jump_increment_m,
            )?;
            jump[i] = finite(old[i] + travel[i])?;
            compatibility[i] = finite(jump[i] - solved.material.point.jump_m()[i])?;
            if compatibility[i].abs() > solved.compatibility_tolerance_m {
                return Err(CohesiveModalError::JumpMismatch {
                    component: i,
                    residual_m: compatibility[i],
                    tolerance_m: solved.compatibility_tolerance_m,
                });
            }
        }
        poll(gate)?;
        let actual = self.point.trial(
            &self.prepared,
            jump,
            self.config.contact.material_energy_tolerance_j_m2,
        )?;
        let area = self.area_m2();
        let mut force_residual = [0.0; 3];
        let mut force_tolerance = [0.0; 3];
        for i in 0..3 {
            let expected = finite(area * actual.mean_traction_pa[i])?;
            limit(
                "actual relaxing bond held force",
                expected.abs(),
                self.config.contact.maximum_force_n,
            )?;
            force_residual[i] = finite(held[i] - expected)?;
            force_tolerance[i] = finite(
                self.config.force_absolute_tolerance_n
                    + self.config.force_relative_tolerance * held[i].abs().max(expected.abs()),
            )?;
            if force_residual[i].abs() > force_tolerance[i] {
                return Err(CohesiveModalError::ForceMismatch {
                    component: i,
                    residual_n: force_residual[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let native = &self.network.staged_frame;
        let network_energy = finite(native.total_energy_j())?;
        let bond_energy = finite(area * actual.point.stored_energy_j_m2())?;
        limit(
            "total energy including relaxing bond",
            finite(network_energy + bond_energy)?,
            self.network.config.maximum_total_energy_j,
        )?;
        let storage = finite(area * actual.stored_change_j_m2)?;
        let heat = finite(area * actual.viscous_dissipated_j_m2)?;
        if heat < 0.0 {
            return Err(invalid("adhesive relaxation cannot create energy").into());
        }
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = dot(&held, &prescribed)?;
        let mut represented = Dd::ZERO;
        for (j, ((force, input), dq)) in self
            .forces
            .iter()
            .zip(external)
            .zip(&self.network.free_delta)
            .enumerate()
        {
            if j % 256 == 0 {
                poll(gate)?;
            }
            represented =
                represented + (Dd::from_f64(*force) - Dd::from_f64(*input)) * Dd::from_f64(*dq);
        }
        let work_on_network = finite(represented.to_f64())?;
        let work = finite(area * actual.work_j_m2)?;
        let scale = storage
            .abs()
            .max(heat.abs())
            .max(work.abs())
            .max(work_on_network.abs())
            .max(prescribed_work.abs());
        let mut tolerance = Dd::from_f64(area * self.config.contact.material_energy_tolerance_j_m2)
            + Dd::from_f64(128.0 * f64::EPSILON * scale);
        let mut transpose = represented;
        for i in 0..3 {
            tolerance =
                tolerance + Dd::from_f64(force_tolerance[i]) * Dd::from_f64(travel[i].abs());
            transpose = transpose + Dd::from_f64(held[i]) * Dd::from_f64(body[i]);
        }
        let bond_tolerance = finite(tolerance.to_f64())?;
        check_work(finite(transpose.to_f64())?, bond_tolerance)?;
        check_work(finite(area * actual.defect_j_m2)?, bond_tolerance)?;
        let bond_residual = finite(
            (Dd::from_f64(storage) + Dd::from_f64(heat) + represented
                - Dd::from_f64(prescribed_work))
            .to_f64(),
        )?;
        check_work(bond_residual, bond_tolerance)?;
        let port_residual =
            Dd::from_f64(native.external_work_j) - Dd::from_f64(external_work) - represented;
        let energy_residual = finite(
            (Dd::from_f64(native.energy_residual_j) + Dd::from_f64(bond_residual) + port_residual)
                .to_f64(),
        )?;
        let energy_tolerance = finite(
            native.energy_tolerance_j
                + bond_tolerance
                + 128.0 * f64::EPSILON * native.external_work_j.abs().max(external_work.abs()),
        )?;
        check_work(energy_residual, energy_tolerance)?;
        let frame = ViscoelasticBondModalFrame {
            sample: native.sample,
            observer_pressure_pa: native.observer_pressure_pa,
            network_energy_j: network_energy,
            bond_energy_j: bond_energy,
            bond_stored_change_j: storage,
            bond_work_j: work,
            relaxation_heat_j: heat,
            network_dissipation_j: native.total_dissipation_j(),
            numerical_dissipation_j: native.numerical_dissipation_j,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            bond_work_on_network_j: work_on_network,
            bond_work_residual_j: bond_residual,
            bond_work_tolerance_j: bond_tolerance,
            energy_residual_j: energy_residual,
            energy_tolerance_j: energy_tolerance,
            jump_m: jump,
            mean_force_n: held,
            endpoint_traction_pa: actual.traction_pa,
            compatibility_residual_m: compatibility,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            constitutive_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
        };
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        Ok((actual.point, frame))
    }
}

fn check_work(residual_j: f64, tolerance_j: f64) -> Result<(), CohesiveModalError> {
    if residual_j.abs() > tolerance_j {
        Err(ModalCouplingError::EnergyBalance {
            residual_j,
            tolerance_j,
        }
        .into())
    } else {
        Ok(())
    }
}
