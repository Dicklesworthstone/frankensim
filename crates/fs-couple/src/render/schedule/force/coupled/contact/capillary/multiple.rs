//! Simultaneous circular wetting bridges on the original mechanical network.
//!
//! All bridges solve `g = g_free + S F(g)` together, with the complete native
//! displacement response `S_ij = B_i^T D B_j`. Every probe starts from the SAME
//! accepted footprints over the SAME physical timestep. The original bridge
//! supplies its two independent pinning/sliding histories, surface free energy,
//! line dissipation and analytic force tangent. Mechanical iteration never
//! advances a liquid history twice or freezes another bridge's reaction.
//!
//! A preallocated dense Newton solve retains every cross-bridge term. The
//! final inverse-Jacobian gain screens unresolved local continuation and
//! propagates the original force-error estimates into gap estimates. These
//! floating-point estimates are not interval certificates or global stability
//! and uniqueness proofs. The original quasistatic, connected, acute-angle,
//! nonpositive-pressure liquid regime remains unchanged. Bridges do not merge,
//! exchange liquid, rupture or acquire liquid inertia through this adapter.
//!
//! One original native step applies the sum of all mean forces. Every bridge
//! is then replayed at its ACTUAL mechanical gap. Native motion, all footprint
//! histories, gap coordinates, diagnostics and clock publish together after
//! independent compatibility, work, source-accuracy and cancellation checks.

use super::{
    AsymmetricCircularContactLineState, AsymmetricCircularContactLineStep, CapillaryModalError,
    CoupledModalSystem, MaterialStateId, ModalAcousticTimeModel, ModalCouplingError,
    SourceAsymmetricCircularContactLineBridge, dot, finite, invalid, limit, network_response, poll,
};
use crate::render::schedule::force::coupled::CoupledModalFrame;
use fs_exec::CancelGate;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};

mod solve;

pub const MAX_CIRCULAR_WETTING_BRIDGES: usize = 32;

/// One original bridge with its complete reciprocal mechanical gap column.
/// Ordered substrate identities are checked before the input is adopted.
pub struct CircularWettingBridge<'a> {
    pub source: &'a SourceAsymmetricCircularContactLineBridge,
    pub substrate_materials: [&'a MaterialStateId; 2],
    pub state: AsymmetricCircularContactLineState,
    pub gap_column: Vec<f64>,
}

/// All outer source calls share ONE ceiling, including rejected line-search
/// probes and the final replay of every bridge. Each source call retains its
/// original bounded geometry/footprint work, so total work is bounded by the
/// product of the outer ceiling and those original source ceilings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleCircularWettingModalConfig {
    pub maximum_bridges: usize,
    pub maximum_setup_terms: usize,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_evaluations: usize,
    pub gap_absolute_tolerance_m: f64,
    pub gap_relative_tolerance: f64,
    pub maximum_gap_increment_m: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    /// Upper admitted infinity norm of the inverse joint compatibility
    /// Jacobian. This is a numerical continuation screen, not a stability proof.
    pub maximum_compatibility_gain: f64,
}

#[derive(Debug)]
pub enum MultipleCircularWettingModalError {
    Native(CapillaryModalError),
    Linear(LuWorkspaceError),
    Budget,
    NotConverged,
    CompatibilityGain {
        gain: f64,
        maximum: f64,
    },
    Bridge {
        index: usize,
        error: CapillaryModalError,
    },
}
impl core::fmt::Display for MultipleCircularWettingModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous native circular wetting refused: {self:?}")
    }
}
impl std::error::Error for MultipleCircularWettingModalError {}
impl From<ModalCouplingError> for MultipleCircularWettingModalError {
    fn from(e: ModalCouplingError) -> Self {
        Self::Native(e.into())
    }
}
impl From<CapillaryModalError> for MultipleCircularWettingModalError {
    fn from(e: CapillaryModalError) -> Self {
        Self::Native(e)
    }
}
impl From<LuWorkspaceError> for MultipleCircularWettingModalError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Linear(e)
    }
}
type Error = MultipleCircularWettingModalError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularWettingBridgeFrame {
    pub applied_mean_force_n: f64,
    /// Independent original proposal on the actual staged body displacement.
    pub bridge: AsymmetricCircularContactLineStep,
    pub actual_gap_m: f64,
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    pub estimated_compatibility_error_m: f64,
    pub constitutive_residual_n: f64,
    pub force_tolerance_n: f64,
    pub estimated_constitutive_force_error_n: f64,
    pub prescribed_work_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleCircularWettingModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    /// Sum of original signed surface energies, with their original references.
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    pub contact_line_dissipation_j: f64,
    pub viscous_line_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub estimated_constitutive_energy_error_j: f64,
    pub compatibility_gain: f64,
    pub iterations: usize,
    /// Total original advance calls, including rejected probes and final replay.
    pub evaluations: usize,
}

struct Point<'a> {
    source: &'a SourceAsymmetricCircularContactLineBridge,
    state: AsymmetricCircularContactLineState,
    gap: f64,
    column: Vec<f64>,
    last: Option<CircularWettingBridgeFrame>,
}

/// Fixed bridge attachments share every original native mode and connection.
/// All nonlinear and candidate buffers are allocated at admission. The existing
/// native propagator retains its own allocation and physical model behavior.
pub struct MultipleCircularWettingModalSystem<'a> {
    network: CoupledModalSystem,
    points: Vec<Point<'a>>,
    config: MultipleCircularWettingModalConfig,
    compliance: Vec<f64>,
    forces: Vec<f64>,
    free_gaps: Vec<f64>,
    gaps: Vec<f64>,
    candidate_gaps: Vec<f64>,
    residual: Vec<f64>,
    correction: Vec<f64>,
    jacobian: Vec<f64>,
    inverse: Vec<f64>,
    gap_errors: Vec<f64>,
    steps: Vec<Option<AsymmetricCircularContactLineStep>>,
    candidate_steps: Vec<Option<AsymmetricCircularContactLineStep>>,
    staged_points: Vec<Option<CircularWettingBridgeFrame>>,
    lu: LuWorkspace,
    last: Option<MultipleCircularWettingModalFrame>,
}

impl<'a> MultipleCircularWettingModalSystem<'a> {
    pub fn new(
        network: CoupledModalSystem,
        bridges: Vec<CircularWettingBridge<'a>>,
        config: MultipleCircularWettingModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(Some(gate))?;
        let p = bridges.len();
        let n = network.mode_count();
        if network.samples_rendered() != 0
            || !(1..=MAX_CIRCULAR_WETTING_BRIDGES).contains(&config.maximum_bridges)
            || p == 0
            || p > config.maximum_bridges
            || config.maximum_iterations == 0
            || config.maximum_backtracks == 0
            || config.maximum_evaluations < 2 * p
        {
            return Err(invalid("multiple wetting needs a sample-zero network and bounded nonempty bridges/iteration budgets").into());
        }
        for x in [
            config.gap_absolute_tolerance_m,
            config.maximum_gap_increment_m,
            config.force_absolute_tolerance_n,
            config.maximum_compatibility_gain,
        ] {
            if !x.is_finite() || x <= 0.0 {
                return Err(invalid(
                    "multiple wetting needs positive finite physical and accuracy limits",
                )
                .into());
            }
        }
        for x in [
            config.gap_relative_tolerance,
            config.force_relative_tolerance,
        ] {
            if !x.is_finite() || !(0.0 < x && x < 1.0) {
                return Err(
                    invalid("multiple wetting relative tolerances must be in (0,1)").into(),
                );
            }
        }
        let k = network.columns.len();
        let terms = n
            .checked_mul(k + 2)
            .and_then(|x| (k + 1).checked_mul(k + 1).and_then(|y| x.checked_add(y)))
            .and_then(|x| x.checked_mul(p))
            .and_then(|x| n.checked_mul(p * p).and_then(|y| x.checked_add(y)))
            .ok_or_else(|| invalid("multiple wetting setup work overflow"))?;
        if terms > config.maximum_setup_terms {
            return Err(Error::Budget);
        }
        let mut points = Vec::with_capacity(p);
        for (i, bridge) in bridges.into_iter().enumerate() {
            poll(Some(gate))?;
            if bridge.gap_column.len() != n || bridge.gap_column.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "every wetting column must cover the original finite native basis",
                )
                .into());
            }
            for side in 0..2 {
                if &bridge.source.original().wetting()[side].surfaces()[0].material
                    != bridge.substrate_materials[side]
                {
                    return Err(invalid(
                        "wetting substrate order differs from original source identities",
                    )
                    .into());
                }
            }
            bridge
                .source
                .evaluate(bridge.state)
                .map_err(|e| Error::Bridge {
                    index: i,
                    error: e.into(),
                })?;
            points.push(Point {
                source: bridge.source,
                state: bridge.state,
                gap: bridge.state.gap_m(),
                column: bridge.gap_column,
                last: None,
            });
        }
        let mut compliance = vec![0.0; p * p];
        for j in 0..p {
            let response = network_response(&network, &points[j].column, gate)?;
            for i in 0..p {
                poll(Some(gate))?;
                compliance[i * p + j] = dot(&points[i].column, &response)?;
            }
            if compliance[j * p + j] <= 0.0 {
                return Err(invalid("each wetting bridge needs positive representable original mechanical compliance").into());
            }
        }
        let owner = Self {
            network,
            points,
            config,
            compliance,
            forces: vec![0.0; n],
            free_gaps: vec![0.0; p],
            gaps: vec![0.0; p],
            candidate_gaps: vec![0.0; p],
            residual: vec![0.0; p],
            correction: vec![0.0; p],
            jacobian: vec![0.0; p * p],
            inverse: vec![0.0; p * p],
            gap_errors: vec![0.0; p],
            steps: vec![None; p],
            candidate_steps: vec![None; p],
            staged_points: vec![None; p],
            lu: LuWorkspace::new(p)?,
            last: None,
        };
        limit(
            "initial native/multiple-wetting energy magnitude",
            owner.total_energy_j()?.abs(),
            owner.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(owner)
    }

    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub fn bridge_count(&self) -> usize {
        self.points.len()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    /// Original body-network account; its external work includes bridge work.
    pub fn native_frame(&self) -> Option<&CoupledModalFrame> {
        self.network.last_frame()
    }
    pub const fn config(&self) -> MultipleCircularWettingModalConfig {
        self.config
    }
    /// Complete row-major cross-bridge mechanical compliance [m/N].
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    pub fn source(&self, index: usize) -> Option<&'a SourceAsymmetricCircularContactLineBridge> {
        self.points.get(index).map(|p| p.source)
    }
    pub fn state(&self, index: usize) -> Option<AsymmetricCircularContactLineState> {
        self.points.get(index).map(|p| p.state)
    }
    pub fn gap_m(&self, index: usize) -> Option<f64> {
        self.points.get(index).map(|p| p.gap)
    }
    pub fn work_column(&self, index: usize) -> Option<&[f64]> {
        self.points.get(index).map(|p| p.column.as_slice())
    }
    pub fn bridge_frame(&self, index: usize) -> Option<&CircularWettingBridgeFrame> {
        self.points.get(index).and_then(|p| p.last.as_ref())
    }
    pub const fn last_frame(&self) -> Option<&MultipleCircularWettingModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, Error> {
        let mut energy = self.network.total_energy_j()?;
        for (index, point) in self.points.iter().enumerate() {
            let profile = point
                .source
                .evaluate(point.state)
                .map_err(|e| Error::Bridge {
                    index,
                    error: e.into(),
                })?;
            energy = finite(energy + profile.surface_free_energy_j)?;
        }
        Ok(energy)
    }

    /// External modal forces and one prescribed additional gap velocity [m/s]
    /// per bridge are held over this sample. Use zero rates for internal bonds.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_gap_rates_m_s: &[f64],
    ) -> Result<&MultipleCircularWettingModalFrame, Error> {
        self.step_under_gate(
            external,
            prescribed_gap_rates_m_s,
            &CancelGate::new_clock_free(),
        )
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_gap_rates_m_s: &[f64],
        gate: &CancelGate,
    ) -> Result<&MultipleCircularWettingModalFrame, Error> {
        poll(Some(gate))?;
        let p = self.points.len();
        if prescribed_gap_rates_m_s.len() != p
            || prescribed_gap_rates_m_s.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "one finite prescribed gap velocity per original bridge is required",
            )
            .into());
        }
        let before = self.total_energy_j()?;
        let dt = self.network.dt;
        self.network.prepare_prediction(external, Some(gate))?;
        for i in 0..p {
            let point = &self.points[i];
            self.free_gaps[i] = finite(
                point.gap
                    + finite(prescribed_gap_rates_m_s[i] * dt)?
                    + dot(&point.column, &self.network.free_delta)?,
            )?;
            self.gaps[i] = point.gap;
        }
        let (iterations, mut evaluations, gain) = self.solve(gate)?;
        self.forces.copy_from_slice(external);
        for i in 0..p {
            let step = self.steps[i]
                .as_ref()
                .expect("successful joint solve fills every bridge");
            for (f, &b) in self.forces.iter_mut().zip(&self.points[i].column) {
                *f = finite(*f + b * step.mean_normal_force_n)?;
            }
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        let mut surface_energy = 0.0;
        let mut prescribed_work = 0.0;
        let mut viscous_loss = 0.0;
        let mut pinning_loss = 0.0;
        let mut energy_error = 0.0;
        for i in 0..p {
            poll(Some(gate))?;
            let point = &self.points[i];
            let solved = self.steps[i]
                .as_ref()
                .expect("successful joint solve fills every bridge");
            let actual_gap = finite(
                point.gap
                    + finite(prescribed_gap_rates_m_s[i] * dt)?
                    + dot(&point.column, &self.network.free_delta)?,
            )?;
            let tolerance = self.gap_tolerance(i, actual_gap)?;
            let residual = finite(actual_gap - self.gaps[i])?;
            if finite(residual.abs() + self.gap_errors[i])? > tolerance {
                return Err(Error::Bridge {
                    index: i,
                    error: CapillaryModalError::GapAccuracy {
                        residual_m: residual,
                        estimated_error_m: self.gap_errors[i],
                        tolerance_m: tolerance,
                    },
                });
            }
            limit(
                "actual multiple-wetting gap increment",
                (actual_gap - point.gap).abs(),
                self.config.maximum_gap_increment_m,
            )?;
            if evaluations == self.config.maximum_evaluations {
                return Err(Error::Budget);
            }
            evaluations += 1;
            let actual = point
                .source
                .advance(point.state, actual_gap, dt, gate)
                .map_err(|e| Error::Bridge {
                    index: i,
                    error: e.into(),
                })?;
            let geometric_residual = finite(actual_gap - actual.state().gap_m())?;
            if geometric_residual.abs() > tolerance {
                return Err(Error::Bridge {
                    index: i,
                    error: CapillaryModalError::GapMismatch {
                        residual_m: geometric_residual,
                        tolerance_m: tolerance,
                    },
                });
            }
            let force_residual = finite(solved.mean_normal_force_n - actual.mean_normal_force_n)?;
            let force_tolerance = self.force_tolerance(
                solved
                    .mean_normal_force_n
                    .abs()
                    .max(actual.mean_normal_force_n.abs()),
            )?;
            let force_error =
                finite(solved.estimated_force_error_n + actual.estimated_force_error_n)?;
            if finite(force_residual.abs() + force_error)? > force_tolerance {
                return Err(Error::Bridge {
                    index: i,
                    error: CapillaryModalError::ForceAccuracy {
                        residual_n: force_residual,
                        estimated_error_n: force_error,
                        tolerance_n: force_tolerance,
                    },
                });
            }
            let work = finite(-solved.mean_normal_force_n * prescribed_gap_rates_m_s[i] * dt)?;
            surface_energy = finite(surface_energy + actual.profile.surface_free_energy_j)?;
            prescribed_work = finite(prescribed_work + work)?;
            viscous_loss = finite(viscous_loss + actual.viscous_line_dissipation_j)?;
            pinning_loss = finite(pinning_loss + actual.pinning_dissipation_j)?;
            energy_error = finite(energy_error + actual.estimated_energy_error_j)?;
            self.staged_points[i] = Some(CircularWettingBridgeFrame {
                applied_mean_force_n: solved.mean_normal_force_n,
                bridge: actual,
                actual_gap_m: actual_gap,
                compatibility_residual_m: residual,
                compatibility_tolerance_m: tolerance,
                estimated_compatibility_error_m: self.gap_errors[i],
                constitutive_residual_n: force_residual,
                force_tolerance_n: force_tolerance,
                estimated_constitutive_force_error_n: force_error,
                prescribed_work_j: work,
            });
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let after = finite(network_energy + surface_energy)?;
        limit(
            "native/multiple-wetting energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let external_work = dot(external, &self.network.free_delta)?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let contact_loss = finite(viscous_loss + pinning_loss)?;
        if viscous_loss < 0.0 || pinning_loss < 0.0 || energy_error < 0.0 {
            return Err(
                invalid("original wetting loss and uncertainty must be nonnegative").into(),
            );
        }
        let loss = finite(network_loss + contact_loss)?;
        let residual = finite(after - before + loss - external_work - prescribed_work)?;
        let scale = before
            .abs()
            .max(after.abs())
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
        if finite(residual.abs() + energy_error)? > tolerance {
            return Err(CapillaryModalError::EnergyAccuracy {
                residual_j: residual,
                estimated_error_j: energy_error,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = MultipleCircularWettingModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            surface_energy_j: surface_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            contact_line_dissipation_j: contact_loss,
            viscous_line_dissipation_j: viscous_loss,
            pinning_dissipation_j: pinning_loss,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            estimated_constitutive_energy_error_j: energy_error,
            compatibility_gain: gain,
            iterations,
            evaluations,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        for (point, staged) in self.points.iter_mut().zip(&self.staged_points) {
            let frame = staged.expect("all original bridges audited before publication");
            point.state = frame.bridge.state();
            point.gap = frame.actual_gap_m;
            point.last = Some(frame);
        }
        Ok(self.last.insert(accepted))
    }
}
