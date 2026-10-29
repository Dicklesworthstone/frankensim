//! Two-axis retained shear on the original coupled modal network.
//!
//! The caller supplies both displacement projections in one fixed orthonormal
//! physical tangent frame. Their full response is C_ij = b_i^T D_network b_j;
//! cross-axis compliance is retained even when the constitutive law is isotropic.
//! The shared matrix adapter queries that same law and accepted history. One
//! combined force then advances the existing network, whose actual displacement
//! and work must agree before either mechanics or contact state is published.
//! Normal load is prescribed. Frame transport and normal evolution require
//! their own owners. The optional heat child couples the original finite-body
//! thermal network without adding a reservoir or replacing retained shear.

use super::super::super::network_response;
use super::super::{
    CoupledModalSystem, ModalAttachment, ModalCouplingError, ModalFrictionConfig,
    ModalFrictionError, dot, finite, invalid, limit, poll, scalar,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_tribo::tangential::matrix::{
    MatrixTangentialConfig, MatrixTangentialError, advance_with_matrix_compliance,
};
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// Both contact directions and original finite-body heat accept together.
#[path = "planar/heat.rs"]
pub mod heat;

/// Original source-bound layer with affine elastic temperature dependence.
#[path = "planar/material_layer.rs"]
pub mod material_layer;

/// Held driver motion in the same two physical axes as the attachment shapes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarTangentialModalDrive {
    pub speed_m_s: [f64; 2],
    pub normal_force_n: f64,
    pub contact_temperature_k: Option<f64>,
    pub nominal_pressure_pa: Option<f64>,
}

/// A complete accepted mechanical and constitutive work interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarTangentialModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub contact_energy_j: f64,
    pub contact_stored_change_j: f64,
    pub network_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub external_work_j: f64,
    pub drive_work_j: f64,
    /// Held work effort on the positive attachments [N].
    pub traction_n: [f64; 2],
    /// Driver advance minus actual accepted relative body advance [m].
    pub slip_distance_m: [f64; 2],
    pub compatibility_residual_m: [f64; 2],
    pub contact_work_residual_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub regime: FrictionRegime,
    pub coupling_iterations: usize,
    pub law_evaluations: usize,
}

/// Original network/source error or a refused shared matrix constitutive solve.
#[derive(Debug)]
pub enum PlanarTangentialModalError {
    Mechanical(ModalFrictionError),
    Interface(MatrixTangentialError),
}

impl core::fmt::Display for PlanarTangentialModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Mechanical(error) => write!(f, "planar modal friction: {error}"),
            Self::Interface(error) => write!(f, "planar modal friction: {error}"),
        }
    }
}
impl std::error::Error for PlanarTangentialModalError {}
impl From<ModalFrictionError> for PlanarTangentialModalError {
    fn from(error: ModalFrictionError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<ModalCouplingError> for PlanarTangentialModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Mechanical(error.into())
    }
}
impl From<MatrixTangentialError> for PlanarTangentialModalError {
    fn from(error: MatrixTangentialError) -> Self {
        Self::Interface(error)
    }
}

/// Two work-conjugate contact axes on one unchanged modal network.
///
/// Attachments on each side must name the same component and retain their
/// original modal basis. Setup prepares two response columns once; stepping
/// reuses them and the force workspace. The law owns all shear memory, static
/// friction, plasticity and constitutive applicability.
pub struct PlanarTangentialModalSystem<L: TangentialLaw> {
    network: CoupledModalSystem,
    law: L,
    state: L::State,
    config: ModalFrictionConfig,
    matrix_config: MatrixTangentialConfig,
    columns: [Vec<f64>; 2],
    compliance: [[f64; 2]; 2],
    forces: Vec<f64>,
    last: Option<PlanarTangentialModalFrame>,
}

struct PlanarTangentialModalCandidate<S> {
    frame: PlanarTangentialModalFrame,
    state: S,
}

impl<L: TangentialLaw> PlanarTangentialModalSystem<L> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        left: [ModalAttachment; 2],
        right: Option<[ModalAttachment; 2]>,
        law: L,
        config: ModalFrictionConfig,
        matrix_config: MatrixTangentialConfig,
        gate: &CancelGate,
    ) -> Result<Self, PlanarTangentialModalError> {
        poll(Some(gate))?;
        matrix_config.validate()?;
        if left[0].component != left[1].component
            || right
                .as_ref()
                .is_some_and(|r| r[0].component != r[1].component)
        {
            return Err(invalid("two tangent axes must retain each original contact body").into());
        }
        let columns = [
            scalar::admit_attachments(&network, &left[0], right.as_ref().map(|r| &r[0]), config)?,
            scalar::admit_attachments(&network, &left[1], right.as_ref().map(|r| &r[1]), config)?,
        ];
        let responses = [
            network_response(&network, &columns[0], gate)?,
            network_response(&network, &columns[1], gate)?,
        ];
        let mut compliance = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                compliance[i][j] = dot(&columns[i], &responses[j])?;
            }
        }
        let state = law.initial_state();
        let contact_energy = finite(law.stored_energy_j(&state))?;
        if contact_energy < 0.0 {
            return Err(invalid("planar contact requires nonnegative reversible storage").into());
        }
        limit(
            "initial planar network and contact storage",
            finite(network.total_energy_j()? + contact_energy)?,
            network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(Self {
            forces: vec![0.0; network.mode_count()],
            network,
            law,
            state,
            config,
            matrix_config,
            columns,
            compliance,
            last: None,
        })
    }

    #[must_use]
    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    #[must_use]
    pub const fn law(&self) -> &L {
        &self.law
    }
    #[must_use]
    pub const fn state(&self) -> &L::State {
        &self.state
    }
    /// Complete original response [m/N], including cross-axis terms.
    #[must_use]
    pub const fn compliance_m_per_n(&self) -> [[f64; 2]; 2] {
        self.compliance
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&PlanarTangentialModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        finite(self.network.total_energy_j()? + self.law.stored_energy_j(&self.state))
    }

    pub fn step(
        &mut self,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
    ) -> Result<&PlanarTangentialModalFrame, PlanarTangentialModalError> {
        self.step_inner(external, drive, None)
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        gate: &CancelGate,
    ) -> Result<&PlanarTangentialModalFrame, PlanarTangentialModalError> {
        self.step_inner(external, drive, Some(gate))
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        gate: Option<&CancelGate>,
    ) -> Result<&PlanarTangentialModalFrame, PlanarTangentialModalError> {
        let candidate = self.stage_candidate(external, drive, gate)?;
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        Ok(self.publish_candidate(candidate))
    }

    fn stage_candidate(
        &mut self,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        gate: Option<&CancelGate>,
    ) -> Result<PlanarTangentialModalCandidate<L::State>, PlanarTangentialModalError> {
        self.stage_candidate_with_law(external, drive, gate, None)
    }

    // Thermal continuation supplies an admitted endpoint law and the SAME
    // retained strain rebound by its original source owner. Accepted fields
    // remain untouched throughout nonlinear/Jacobian probes.
    fn stage_candidate_with_law(
        &mut self,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        gate: Option<&CancelGate>,
        constitutive: Option<(&L, &L::State)>,
    ) -> Result<PlanarTangentialModalCandidate<L::State>, PlanarTangentialModalError> {
        poll(gate)?;
        let (law, state) = constitutive.unwrap_or((&self.law, &self.state));
        if drive.speed_m_s.iter().any(|v| !v.is_finite())
            || !drive.normal_force_n.is_finite()
            || drive.normal_force_n < 0.0
            || drive
                .contact_temperature_k
                .is_some_and(|t| !t.is_finite() || t <= 0.0)
            || drive
                .nominal_pressure_pa
                .is_some_and(|p| !p.is_finite() || p < 0.0)
        {
            return Err(invalid("planar drive requires finite physical contact inputs").into());
        }
        let network_before = self.network.total_energy_j()?;
        let contact_before = finite(law.stored_energy_j(state))?;
        self.network.prepare_prediction(external, gate)?;
        let mut drive_distance = [0.0; 2];
        let mut free_slip = [0.0; 2];
        for i in 0..2 {
            drive_distance[i] = finite(drive.speed_m_s[i] * self.network.dt)?;
            free_slip[i] =
                finite(drive_distance[i] - dot(&self.columns[i], &self.network.free_delta)?)?;
        }
        let result = advance_with_matrix_compliance(
            law,
            state,
            &TangentialStepInput {
                normal_n: drive.normal_force_n,
                slip_increment_m: free_slip,
                dt_s: self.network.dt,
                contact_temperature_k: drive.contact_temperature_k,
                nominal_pressure_pa: drive.nominal_pressure_pa,
                compliance_m_per_n: 0.0,
                solve: TangentialSolveConfig {
                    max_iterations: self.config.max_iterations,
                    maximum_force_n: self.config.maximum_force_n,
                    force_absolute_tolerance_n: self.config.force_absolute_tolerance_n,
                    force_relative_tolerance: self.config.force_relative_tolerance,
                },
            },
            self.compliance,
            self.matrix_config,
            || gate.is_some_and(CancelGate::is_requested),
        )?;
        poll(gate)?;
        let step = result.contact;
        let magnitude = finite(det::hypot(step.traction_n[0], step.traction_n[1]))?;
        limit(
            "planar contact traction",
            magnitude,
            self.config.maximum_force_n,
        )?;
        self.forces.copy_from_slice(external);
        for ((force, bx), by) in self
            .forces
            .iter_mut()
            .zip(&self.columns[0])
            .zip(&self.columns[1])
        {
            *force = finite(*force + bx * step.traction_n[0] + by * step.traction_n[1])?;
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut slip = [0.0; 2];
        let mut compatibility = [0.0; 2];
        for i in 0..2 {
            slip[i] = finite(drive_distance[i] - dot(&self.columns[i], &self.network.free_delta)?)?;
            compatibility[i] = finite(slip[i] - step.slip_increment_m[i])?;
            for residual in [compatibility[i], step.compatibility_residual_m[i]] {
                limit(
                    "actual planar contact compatibility",
                    finite(residual)?.abs(),
                    self.config.sticking_tolerance_m,
                )?;
            }
        }
        limit(
            "planar contact slip speed",
            finite(det::hypot(slip[0], slip[1]) / self.network.dt)?,
            self.config.maximum_slip_speed_m_s,
        )?;
        let contact_after = finite(law.stored_energy_j(&step.state))?;
        let loss = finite(step.dissipated_j)?;
        if contact_before < 0.0 || contact_after < 0.0 || loss < 0.0 {
            return Err(invalid("planar storage and friction heat must be nonnegative").into());
        }
        let contact_change = finite(contact_after - contact_before)?;
        let frame = &self.network.staged_frame;
        let network_after = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let drive_work = dot(&step.traction_n, &drive_distance)?;
        let before = finite(network_before + contact_before)?;
        let after = finite(network_after + contact_after)?;
        limit(
            "planar network and contact storage",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let scale = before
            .max(after)
            .max(network_loss)
            .max(loss)
            .max(external_work.abs())
            .max(drive_work.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        let contact_work_residual = finite(contact_change + loss - dot(&step.traction_n, &slip)?)?;
        let residual = finite((after - before) + network_loss + loss - external_work - drive_work)?;
        for checked in [
            residual,
            contact_work_residual,
            finite(contact_change - step.stored_change_j)?,
            finite(step.closure_residual_j)?,
        ] {
            if checked.abs() > tolerance {
                return Err(ModalCouplingError::EnergyBalance {
                    residual_j: checked,
                    tolerance_j: tolerance,
                }
                .into());
            }
        }
        let accepted = PlanarTangentialModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_after,
            contact_energy_j: contact_after,
            contact_stored_change_j: contact_change,
            network_dissipation_j: network_loss,
            friction_dissipation_j: loss,
            external_work_j: external_work,
            drive_work_j: drive_work,
            traction_n: step.traction_n,
            slip_distance_m: slip,
            compatibility_residual_m: compatibility,
            contact_work_residual_j: contact_work_residual,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            regime: step.regime,
            coupling_iterations: result.iterations,
            law_evaluations: result.law_evaluations,
        };
        poll(gate)?;
        Ok(PlanarTangentialModalCandidate {
            frame: accepted,
            state: step.state,
        })
    }

    fn publish_candidate(
        &mut self,
        candidate: PlanarTangentialModalCandidate<L::State>,
    ) -> &PlanarTangentialModalFrame {
        self.network.publish_staged();
        self.state = candidate.state;
        self.last.insert(candidate.frame)
    }
}
