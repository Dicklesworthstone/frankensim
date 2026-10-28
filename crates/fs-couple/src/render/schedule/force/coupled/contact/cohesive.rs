//! Irreversible opening/shear debonding on the original native modal network.
//!
//! Each supplied work column maps every retained modal coordinate into one
//! fixed local jump component. Its transpose applies the same physical force.
//! The full held-force compliance, native propagation, original component
//! limits and actual-state energy transaction remain with that network.
//! Construction explicitly bonds the current configuration at zero jump.
//!
//! Mean path traction supplies body work; endpoint traction remains a separate
//! observation. Fracture and physical bond-viscosity loss are reported separately
//! and are not automatically assigned to heat. Straight jump paths and a fixed
//! local frame are the constitutive approximation: resolve reversals by refining
//! the physical timestep. No evolving geometry, rebonding or thermal refresh is
//! inferred from this fixed-source interface.

use super::{
    CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, network_response, poll,
};
use crate::cohesive_contact::{CohesiveContactConfig, CohesiveContactError, CohesiveLoadLine};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::cohesive::{BilinearCohesive, CohesiveError, CohesivePoint};
use fs_material::state_point::MaterialStatePointError;

/// A finite bonded interface with coupled opening, shear and progressive peeling.
pub mod patch;

/// Bounded constitutive solve plus the independent actual-body force gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveModalConfig {
    /// The material-evaluation ceiling includes the final actual-state audit
    /// and must be at least two; one evaluation is reserved from the root solve.
    pub contact: CohesiveContactConfig,
    /// Positive absolute mismatch allowed in each physical force component [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative mismatch in (0,1), scaled by actual/applied force, not its ceiling.
    pub force_relative_tolerance: f64,
}

/// No refusal publishes modal state, irreversible damage or the last frame.
#[derive(Debug)]
pub enum CohesiveModalError {
    Network(ModalCouplingError),
    Contact(CohesiveContactError),
    Material(CohesiveError),
    Source(MaterialStatePointError),
    JumpMismatch {
        component: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    ForceMismatch {
        component: usize,
        residual_n: f64,
        tolerance_n: f64,
    },
}
impl core::fmt::Display for CohesiveModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native cohesive trial refused: {self:?}")
    }
}
impl std::error::Error for CohesiveModalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Network(e) => Some(e),
            Self::Contact(e) => Some(e),
            Self::Material(e) => Some(e),
            Self::Source(e) => Some(e),
            _ => None,
        }
    }
}
impl From<ModalCouplingError> for CohesiveModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<CohesiveContactError> for CohesiveModalError {
    fn from(error: CohesiveContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<CohesiveError> for CohesiveModalError {
    fn from(error: CohesiveError) -> Self {
        Self::Material(error)
    }
}
impl From<MaterialStatePointError> for CohesiveModalError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}

/// Last accepted complete mechanical and cohesive work account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    /// Recoverable bond plus independent compression-penalty storage [J].
    pub cohesive_energy_j: f64,
    /// Authored modal-force work, excluding the internal cohesive reaction [J].
    pub external_work_j: f64,
    /// Applied interface force dotted with prescribed jump increment [J].
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    /// Intrinsic irreversible fracture, not automatically frictional heat [J].
    pub fracture_dissipation_j: f64,
    /// Separately integrated physical surviving-bond dashpot loss [J].
    pub viscous_dissipation_j: f64,
    /// Actual constitutive path work on all three local jump components [J].
    pub cohesive_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Held force actually applied through the complete work columns [N].
    pub mean_force_n: [f64; 3],
    /// Timed endpoint traction, distinct from the force driving the bodies [Pa].
    pub endpoint_traction_pa: [f64; 3],
    pub jump_m: [f64; 3],
    pub damage: f64,
    /// Actual body jump minus the nonlinear load-line solution [m].
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    /// Applied force minus the constitutive mean on actual body motion [N].
    pub constitutive_residual_n: [f64; 3],
    pub force_tolerance_n: [f64; 3],
    pub iterations: usize,
    /// Includes the final independent actual-state material evaluation.
    pub material_evaluations: usize,
}

/// One explicitly bonded patch in a fixed local frame. The original network
/// and source remain in place; stepping clones neither mechanical nor material
/// owners. Only the original propagator's existing scratch/allocation is used.
pub struct CohesiveModalSystem<'a> {
    network: CoupledModalSystem,
    point: CohesivePoint<'a>,
    source: Option<&'a ResolvedBilinearCohesive>,
    columns: [Vec<f64>; 3],
    load_line: CohesiveLoadLine,
    config: CohesiveModalConfig,
    forces: Vec<f64>,
    last: Option<CohesiveModalFrame>,
}

impl<'a> CohesiveModalSystem<'a> {
    /// Explicitly bond a sample-zero network in its current configuration.
    /// Columns are positive opening, shear 1 and shear 2, in the complete
    /// component-then-mode order. `B_i^T delta_q` is a physical jump [m]; the
    /// reciprocal modal reaction is `-sum B_i F_i`. Area is original bonded m2.
    /// The initial jump and cohesive storage are zero even for supplied initial
    /// vibration; no preload equilibrium or bond creation from proximity occurs.
    pub fn new(
        network: CoupledModalSystem,
        law: &'a BilinearCohesive,
        area_m2: f64,
        columns: [Vec<f64>; 3],
        config: CohesiveModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CohesiveModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0 {
            return Err(invalid("cohesive bonding requires a sample-zero network").into());
        }
        if columns.iter().any(|column| {
            column.len() != network.mode_count() || column.iter().any(|value| !value.is_finite())
        }) || !config.force_absolute_tolerance_n.is_finite()
            || config.force_absolute_tolerance_n <= 0.0
            || !config.force_relative_tolerance.is_finite()
            || !(0.0 < config.force_relative_tolerance && config.force_relative_tolerance < 1.0)
            || config.contact.max_material_evaluations < 2
        {
            return Err(invalid("cohesive work columns must retain every mode and force gates must be positive finite").into());
        }
        let mut compliance = [[0.0; 3]; 3];
        for j in 0..3 {
            let response = network_response(&network, &columns[j], gate)?;
            for i in 0..3 {
                compliance[i][j] = dot(&columns[i], &response)?;
            }
        }
        // Reserve one evaluation for the independent actual-state audit so the
        // native transaction never exceeds the caller's complete material budget.
        let mut solve_config = config.contact;
        solve_config.max_material_evaluations -= 1;
        let load_line = CohesiveLoadLine::new(area_m2, compliance, solve_config)?;
        let system = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            point: law.bonded_point(),
            source: None,
            columns,
            load_line,
            config,
            last: None,
        };
        limit(
            "total initial energy including cohesive bond",
            system.total_energy_j()?,
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
    }

    /// Retain the original ordered interface source at its fixed external query.
    /// The caller supplies the material identities actually assigned to the two
    /// receiving surfaces; reversed or foreign assignments are refused. Frame,
    /// geometry and physical bonding are the caller's explicit declarations.
    pub fn from_source(
        network: CoupledModalSystem,
        source: &'a ResolvedBilinearCohesive,
        surface_materials: [&MaterialStateId; 2],
        area_m2: f64,
        columns: [Vec<f64>; 3],
        config: CohesiveModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CohesiveModalError> {
        poll(Some(gate))?;
        source.require_surface_materials(surface_materials[0], surface_materials[1])?;
        let mut system = Self::new(network, source.law(), area_m2, columns, config, gate)?;
        system.source = Some(source);
        Ok(system)
    }

    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    #[must_use]
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    #[must_use]
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    #[must_use]
    pub const fn point(&self) -> CohesivePoint<'a> {
        self.point
    }
    #[must_use]
    pub const fn source(&self) -> Option<&'a ResolvedBilinearCohesive> {
        self.source
    }
    #[must_use]
    pub fn law(&self) -> &'a BilinearCohesive {
        self.point.law()
    }
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        self.load_line.area_m2()
    }
    #[must_use]
    pub fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.load_line.compliance_m_per_n()
    }
    #[must_use]
    pub const fn work_columns(&self) -> &[Vec<f64>; 3] {
        &self.columns
    }
    #[must_use]
    pub const fn config(&self) -> CohesiveModalConfig {
        self.config
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&CohesiveModalFrame> {
        self.last.as_ref()
    }

    /// Accepted original mechanical storage plus recoverable interface storage.
    pub fn total_energy_j(&self) -> Result<f64, CohesiveModalError> {
        Ok(finite(
            self.network.total_energy_j()? + self.area_m2() * self.point.stored_energy_j_m2(),
        )?)
    }

    /// Advance held external modal forces and a constant prescribed local jump
    /// rate [m/s]. Use `[0.0; 3]` for a wholly internal two-body interface.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: [f64; 3],
    ) -> Result<&CohesiveModalFrame, CohesiveModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, None)
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<&CohesiveModalFrame, CohesiveModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, Some(gate))
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        prescribed_rate: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<&CohesiveModalFrame, CohesiveModalError> {
        poll(gate)?;
        let before = self.total_energy_j()?;
        let old = self.point.jump_m();
        let dt = self.network.dt;
        let mut prescribed_delta = [0.0; 3];
        for i in 0..3 {
            prescribed_delta[i] = finite(prescribed_rate[i] * dt)?;
        }
        self.network.prepare_prediction(external, gate)?;
        let mut free_jump = [0.0; 3];
        for i in 0..3 {
            free_jump[i] = finite(
                old[i] + prescribed_delta[i] + dot(&self.columns[i], &self.network.free_delta)?,
            )?;
        }
        let solved = self.load_line.trial(self.point, free_jump, dt, || {
            !gate.is_some_and(CancelGate::is_requested)
        })?;
        self.forces.copy_from_slice(external);
        for i in 0..3 {
            for (force, &weight) in self.forces.iter_mut().zip(&self.columns[i]) {
                *force = finite(*force - weight * solved.mean_force_n[i])?;
            }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut actual_jump = [0.0; 3];
        let mut compatibility_residual = [0.0; 3];
        let solved_jump = solved.material.equilibrium.point.jump_m();
        for i in 0..3 {
            actual_jump[i] = finite(
                old[i] + prescribed_delta[i] + dot(&self.columns[i], &self.network.free_delta)?,
            )?;
            compatibility_residual[i] = finite(actual_jump[i] - solved_jump[i])?;
            if compatibility_residual[i].abs() > solved.compatibility_tolerance_m {
                return Err(CohesiveModalError::JumpMismatch {
                    component: i,
                    residual_m: compatibility_residual[i],
                    tolerance_m: solved.compatibility_tolerance_m,
                });
            }
            limit(
                "actual cohesive jump increment",
                (actual_jump[i] - old[i]).abs(),
                self.config.contact.maximum_jump_increment_m,
            )?;
        }
        poll(gate)?;
        let (actual, mean) = self.point.trial_over_time_with_mean_traction(
            actual_jump,
            dt,
            self.config.contact.material_energy_tolerance_j_m2,
        )?;
        let area = self.area_m2();
        let mut force_residual = [0.0; 3];
        let mut force_tolerance = [0.0; 3];
        for i in 0..3 {
            let expected = finite(area * mean[i])?;
            limit(
                "actual cohesive held force",
                expected.abs(),
                self.config.contact.maximum_force_n,
            )?;
            force_residual[i] = finite(solved.mean_force_n[i] - expected)?;
            force_tolerance[i] = finite(
                self.config.force_absolute_tolerance_n
                    + self.config.force_relative_tolerance
                        * expected.abs().max(solved.mean_force_n[i].abs()),
            )?;
            if force_residual[i].abs() > force_tolerance[i] {
                return Err(CohesiveModalError::ForceMismatch {
                    component: i,
                    residual_n: force_residual[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let cohesive_energy = finite(area * actual.equilibrium.point.stored_energy_j_m2())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let fracture = finite(area * actual.equilibrium.energy.dissipated)?;
        let viscous = finite(area * actual.viscous_dissipated_j_m2)?;
        if fracture < 0.0 || viscous < 0.0 {
            return Err(invalid("cohesive fracture and viscosity must not create energy").into());
        }
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = dot(&solved.mean_force_n, &prescribed_delta)?;
        let after = finite(network_energy + cohesive_energy)?;
        limit(
            "total energy including cohesive interface",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let loss = finite(network_loss + fracture + viscous)?;
        let residual = finite((after - before) + loss - external_work - prescribed_work)?;
        let scale = before
            .max(after)
            .max(external_work.abs())
            .max(prescribed_work.abs())
            .max(loss.abs());
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
        let accepted = CohesiveModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            cohesive_energy_j: cohesive_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            fracture_dissipation_j: fracture,
            viscous_dissipation_j: viscous,
            cohesive_work_j: finite(area * actual.total_work_j_m2())?,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            mean_force_n: solved.mean_force_n,
            endpoint_traction_pa: actual.traction_pa,
            jump_m: actual_jump,
            damage: actual.equilibrium.point.damage(),
            compatibility_residual_m: compatibility_residual,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            constitutive_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
            iterations: solved.iterations,
            material_evaluations: solved.material_evaluations + 1,
        };
        poll(gate)?;
        self.network.publish_staged();
        self.point = actual.equilibrium.point;
        Ok(self.last.insert(accepted))
    }
}
