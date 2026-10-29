//! Work-conjugate stateful contact on the original exact modal network.
//!
//! The constitutive owner receives the released relative increment and the
//! actual network compliance. Its held effort drives that same network; the
//! accepted displacement, contact storage and physical heat close one work
//! window. No endpoint force is substituted for a law's discrete work effort.
//! Each law integrates one declared displacement segment; time refinement must
//! resolve intrastep reversals. Exact linear propagation and energy closure do
//! not establish exact nonlinear contact trajectories.
//! The single retained tangent and prescribed normal load do not claim normal
//! contact evolution or two-dimensional slip. The optional [`heat`] composition
//! solves the interface temperature and deposits physical heat into the original
//! finite thermal network without replacing the constitutive law or its history.

use super::{
    CoupledModalSystem, ModalAttachment, ModalCouplingError, ModalFrictionConfig,
    ModalFrictionError, dot, effective_compliance, finite, invalid, limit, poll, scalar,
};
use fs_exec::CancelGate;
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// Stateful shear and the original finite thermal network accept together.
pub mod heat;

/// Two retained tangent directions with the full original network response.
pub mod planar;

/// Simultaneous retained shear at several contacts on the same modal bodies.
pub mod multiple;

/// One held driver motion and its explicitly supplied contact environment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentialModalDrive {
    /// Signed driver-minus-receiver speed [m/s]. Zero gives an internal brake.
    pub speed_m_s: f64,
    /// Prescribed compressive load [N]. The constitutive law owns opening rules.
    pub normal_force_n: f64,
    /// Explicit thermal boundary [K]; absence is passed to the original law.
    pub contact_temperature_k: Option<f64>,
    /// Explicit nominal pressure [Pa]; never inferred from the normal load.
    pub nominal_pressure_pa: Option<f64>,
}

/// Accepted mechanics and constitutive work on one shared sample clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentialModalFrame {
    /// Original network's accepted one-based sample number.
    pub sample: u64,
    /// Original coupled observer [Pa], without extra scaling or synthesis.
    pub observer_pressure_pa: f64,
    /// Original network's reversible storage, excluding contact storage [J].
    pub network_energy_j: f64,
    /// Reversible contact storage in the accepted constitutive state [J].
    pub contact_energy_j: f64,
    /// Change of that contact storage during this sample [J].
    pub contact_stored_change_j: f64,
    /// Original network's material, radiation and connection dissipation [J].
    pub network_dissipation_j: f64,
    /// Exported physical contact heat; reversible contact work is excluded [J].
    pub friction_dissipation_j: f64,
    /// Work of original external generalized forces [J].
    pub external_work_j: f64,
    /// Work supplied by prescribed driver motion [J].
    pub drive_work_j: f64,
    /// Held work-conjugate traction on the positive attachment [N].
    pub traction_n: f64,
    /// Actual driver advance minus accepted attachment increment [m].
    pub slip_distance_m: f64,
    /// Actual increment minus the constitutive solve's returned increment [m].
    pub compatibility_residual_m: f64,
    /// Actual contact storage change + physical heat - actual traction work [J].
    pub contact_work_residual_j: f64,
    /// Complete storage change + losses - original external supplies [J].
    pub energy_residual_j: f64,
    /// Original network's absolute plus scaled relative energy budget [J].
    pub energy_tolerance_j: f64,
    /// Original law's branch; a compliant sticking branch can move elastically.
    pub regime: FrictionRegime,
    /// Constitutive solver iterations, without an additional coupling root loop.
    pub iterations: usize,
}

/// A single stateful tangential port on an already admitted coupled network.
///
/// Both attachments, all bilateral connections, material memory and prepared
/// acoustic loads participate in the original response. Configure those on the
/// network before construction. The law remains immutable; its state, mechanics
/// and frame publish together only after the actual work checks succeed.
pub struct TangentialModalSystem<L: TangentialLaw> {
    network: CoupledModalSystem,
    law: L,
    state: L::State,
    config: ModalFrictionConfig,
    column: Vec<f64>,
    compliance: f64,
    forces: Vec<f64>,
    last: Option<TangentialModalFrame>,
}

struct TangentialModalCandidate<S> {
    frame: TangentialModalFrame,
    state: S,
}

impl<L: TangentialLaw> TangentialModalSystem<L> {
    /// Attach an unloaded law without replacing or resetting any body state.
    /// Source identity and constitutive applicability remain the law's concern.
    pub fn new(
        network: CoupledModalSystem,
        left: ModalAttachment,
        right: Option<ModalAttachment>,
        law: L,
        config: ModalFrictionConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        let column = scalar::admit_attachments(&network, &left, right.as_ref(), config)?;
        let compliance = effective_compliance(&network, &column, gate)?;
        let state = law.initial_state();
        let contact_energy = finite(law.stored_energy_j(&state))?;
        if contact_energy < 0.0 {
            return Err(invalid("tangential law requires nonnegative reversible storage").into());
        }
        limit(
            "initial network and contact storage",
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
            column,
            compliance,
            last: None,
        })
    }

    /// Original accepted network; mutable access cannot bypass contact history.
    #[must_use]
    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    /// Unmodified original constitutive law.
    #[must_use]
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Accepted checkpoint; failed trials do not replace it.
    #[must_use]
    pub const fn state(&self) -> &L::State {
        &self.state
    }
    /// Last complete accepted mechanical/constitutive work window.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&TangentialModalFrame> {
        self.last.as_ref()
    }
    /// Current reversible network plus interface storage [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        finite(self.network.total_energy_j()? + self.law.stored_energy_j(&self.state))
    }

    /// Advance one physical interval using the original exact held-force response.
    pub fn step(
        &mut self,
        external: &[f64],
        drive: TangentialModalDrive,
    ) -> Result<&TangentialModalFrame, ModalFrictionError> {
        self.step_inner(external, drive, None)
    }

    /// Same atomic operation with cancellation at prediction, law and commit boundaries.
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        drive: TangentialModalDrive,
        gate: &CancelGate,
    ) -> Result<&TangentialModalFrame, ModalFrictionError> {
        self.step_inner(external, drive, Some(gate))
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        drive: TangentialModalDrive,
        gate: Option<&CancelGate>,
    ) -> Result<&TangentialModalFrame, ModalFrictionError> {
        let candidate = self.stage_candidate(external, drive, gate)?;
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        Ok(self.publish_candidate(candidate))
    }

    // Every thermal trial starts from the same accepted body and contact state.
    // Only the final jointly admitted trial may cross the publication boundary.
    fn stage_candidate(
        &mut self,
        external: &[f64],
        drive: TangentialModalDrive,
        gate: Option<&CancelGate>,
    ) -> Result<TangentialModalCandidate<L::State>, ModalFrictionError> {
        poll(gate)?;
        if !drive.speed_m_s.is_finite()
            || !drive.normal_force_n.is_finite()
            || drive.normal_force_n < 0.0
            || drive
                .contact_temperature_k
                .is_some_and(|t| !t.is_finite() || t <= 0.0)
            || drive
                .nominal_pressure_pa
                .is_some_and(|p| !p.is_finite() || p < 0.0)
        {
            return Err(invalid(
                "tangential drive requires finite motion and physical contact inputs",
            )
            .into());
        }
        let network_before = self.network.total_energy_j()?;
        let contact_before = finite(self.law.stored_energy_j(&self.state))?;
        self.network.prepare_prediction(external, gate)?;
        let drive_distance = finite(drive.speed_m_s * self.network.dt)?;
        let free_slip = finite(drive_distance - dot(&self.column, &self.network.free_delta)?)?;
        let step = self.law.advance(
            &self.state,
            &TangentialStepInput {
                normal_n: drive.normal_force_n,
                slip_increment_m: [free_slip, 0.0],
                dt_s: self.network.dt,
                contact_temperature_k: drive.contact_temperature_k,
                nominal_pressure_pa: drive.nominal_pressure_pa,
                compliance_m_per_n: self.compliance,
                solve: TangentialSolveConfig {
                    max_iterations: self.config.max_iterations,
                    maximum_force_n: self.config.maximum_force_n,
                    force_absolute_tolerance_n: self.config.force_absolute_tolerance_n,
                    force_relative_tolerance: self.config.force_relative_tolerance,
                },
            },
        )?;
        poll(gate)?;
        // This owner has one work coordinate. A transverse response must not
        // disappear just because its body attachment was never assembled.
        if step.traction_n[1] != 0.0 || step.slip_increment_m[1] != 0.0 {
            return Err(invalid("single-tangent contact refuses transverse response").into());
        }
        let traction = finite(step.traction_n[0])?;
        limit(
            "tangential contact traction",
            traction.abs(),
            self.config.maximum_force_n,
        )?;
        for residual in step.compatibility_residual_m {
            limit(
                "constitutive contact compatibility",
                finite(residual)?.abs(),
                self.config.sticking_tolerance_m,
            )?;
        }
        self.forces.copy_from_slice(external);
        for (force, b) in self.forces.iter_mut().zip(&self.column) {
            *force = finite(*force + b * traction)?;
        }
        self.network.stage_inner(&self.forces, gate)?;
        let slip = finite(drive_distance - dot(&self.column, &self.network.free_delta)?)?;
        let compatibility = finite(slip - step.slip_increment_m[0])?;
        limit(
            "accepted tangential contact compatibility",
            compatibility.abs(),
            self.config.sticking_tolerance_m,
        )?;
        limit(
            "tangential contact slip speed",
            finite(slip / self.network.dt)?.abs(),
            self.config.maximum_slip_speed_m_s,
        )?;
        let contact_after = finite(self.law.stored_energy_j(&step.state))?;
        let loss = finite(step.dissipated_j)?;
        if contact_before < 0.0 || contact_after < 0.0 || loss < 0.0 {
            return Err(invalid("tangential storage and physical heat must be nonnegative").into());
        }
        let contact_change = finite(contact_after - contact_before)?;
        let frame = &self.network.staged_frame;
        let network_after = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let drive_work = finite(traction * drive_distance)?;
        let before = finite(network_before + contact_before)?;
        let after = finite(network_after + contact_after)?;
        limit(
            "network and contact storage",
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
        let contact_work_residual = finite(contact_change + loss - traction * slip)?;
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
        let accepted = TangentialModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_after,
            contact_energy_j: contact_after,
            contact_stored_change_j: contact_change,
            network_dissipation_j: network_loss,
            friction_dissipation_j: loss,
            external_work_j: external_work,
            drive_work_j: drive_work,
            traction_n: traction,
            slip_distance_m: slip,
            compatibility_residual_m: compatibility,
            contact_work_residual_j: contact_work_residual,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            regime: step.regime,
            iterations: step.iterations,
        };
        poll(gate)?;
        Ok(TangentialModalCandidate {
            frame: accepted,
            state: step.state,
        })
    }

    fn publish_candidate(
        &mut self,
        candidate: TangentialModalCandidate<L::State>,
    ) -> &TangentialModalFrame {
        self.network.publish_staged();
        self.state = candidate.state;
        self.last.insert(candidate.frame)
    }
}
