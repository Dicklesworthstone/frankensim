//! Low-Tabor adhesive release, recontact and microslip on the native network.
//!
//! The source retains DMT normal potential and Hertz--Mindlin shear memory.
//! Two complete work columns retain the FULL reciprocal normal/shear response.
//! For each candidate normal approach, the original analytic Mindlin path owner
//! solves shear recoil; one bounded scalar bracket closes normal compatibility.
//! No contact history is allocated until replay on the actual native endpoint.
//!
//! Coordinates and source conditions are fixed over a step, with affine local
//! paths and one proportional shear direction. Signed indentation admits gaps,
//! release and reattachment; the original path owner integrates only the active
//! contact interval and retains total slip during flight. The caller supplies
//! the signed root bracket and resolves loading reversals. A bracketed root is
//! not a uniqueness or snap-stability certificate.

use super::{
    CancelGate, CoupledModalSystem, Dd, MaterialStateId, ModalAcousticTimeModel,
    ModalCouplingError, dot, finite, invalid, limit, network_response, poll,
};
use fs_contact::interface_binding::mindlin::loading::path::SourceMindlinPathInput;
use fs_contact::normal_patch::dmt::mindlin::{
    DmtMindlinError, DmtMindlinPathTrial, DmtMindlinState, SourceDmtMindlinContact,
};
use fs_contact::normal_patch::dmt::{DmtError, DmtResponse};
use fs_tribo::TangentialSolveConfig;

/// Physical path and numerical budgets. The shear force ceiling also bounds
/// the magnitude of the signed normal effort. No material coefficient lives here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveDmtMindlinConfig {
    /// Total original constitutive probes, including one actual-motion replay.
    pub maximum_evaluations: usize,
    pub maximum_step_s: f64,
    pub maximum_coordinate_increment_m: f64,
    pub maximum_coordinate_speed_m_s: f64,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    /// Contact-sized work allowance, independent of unrelated native storage.
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
    pub shear: TangentialSolveConfig,
}

#[derive(Debug)]
pub enum AdhesiveDmtMindlinError {
    Network(ModalCouplingError),
    Contact(DmtMindlinError),
    Normal(DmtError),
    Bracket {
        indentation_m: [f64; 2],
        residual_m: [f64; 2],
    },
    NotConverged {
        evaluations: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    Coordinate {
        component: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    Force {
        component: usize,
        residual_n: f64,
        tolerance_n: f64,
    },
    Work {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for AdhesiveDmtMindlinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native DMT--Mindlin motion: {self:?}")
    }
}
impl std::error::Error for AdhesiveDmtMindlinError {}
impl From<ModalCouplingError> for AdhesiveDmtMindlinError {
    fn from(e: ModalCouplingError) -> Self {
        Self::Network(e)
    }
}
impl From<DmtMindlinError> for AdhesiveDmtMindlinError {
    fn from(e: DmtMindlinError) -> Self {
        Self::Contact(e)
    }
}
impl From<DmtError> for AdhesiveDmtMindlinError {
    fn from(e: DmtError) -> Self {
        Self::Normal(e)
    }
}

/// Actual native and original constitutive endpoint. Tangential contact storage,
/// microslip heat, normal adhesion potential and native numerical loss stay separate.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveDmtMindlinFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub normal_energy_j: f64,
    pub tangential_energy_j: f64,
    pub normal_stored_change_j: f64,
    pub tangential_stored_change_j: f64,
    pub microslip_dissipation_j: f64,
    pub network_dissipation_j: f64,
    pub numerical_dissipation_j: f64,
    pub external_work_j: f64,
    pub prescribed_work_j: f64,
    pub contact_work_on_network_j: f64,
    pub contact_work_residual_j: f64,
    pub contact_work_tolerance_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Work-average [normal, shear] efforts actually applied to native bodies.
    pub mean_force_n: [f64; 2],
    pub position_m: [f64; 2],
    pub endpoint_shear_force_n: f64,
    pub history_intervals: usize,
    pub normal_response: DmtResponse,
    /// Fraction of this affine step with positive contact indentation.
    pub contact_fraction: f64,
    /// Normalized active-contact start/end; [0,0] for a wholly separated path.
    pub contact_interval_fraction: [f64; 2],
    pub compatibility_residual_m: [f64; 2],
    pub compatibility_tolerance_m: [f64; 2],
    pub constitutive_residual_n: [f64; 2],
    pub evaluations: usize,
    /// Inner shear-recoil iterations at the selected normal root.
    pub shear_iterations: usize,
}

/// One source-bound adhesive contact on the complete unchanged body network.
pub struct AdhesiveDmtMindlinModalSystem<'a> {
    network: CoupledModalSystem,
    source: &'a SourceDmtMindlinContact,
    state: DmtMindlinState,
    columns: [Vec<f64>; 2],
    compliance: [[f64; 2]; 2],
    forces: Vec<f64>,
    config: AdhesiveDmtMindlinConfig,
    last: Option<AdhesiveDmtMindlinFrame>,
}

struct Probe {
    indentation: f64,
    path: DmtMindlinPathTrial,
    residual: f64,
    tolerance: f64,
}

impl<'a> AdhesiveDmtMindlinModalSystem<'a> {
    /// Bind a virgin shear history at the supplied actual initial indentation.
    /// Nonzero initial modal position is allowed; the attachment coordinates
    /// are explicit, not reconstructed from an invented absolute modal origin.
    /// The columns map dq to [approach, shear] and apply reaction -B*mean_force.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        source: &'a SourceDmtMindlinContact,
        surface_materials: [&MaterialStateId; 2],
        initial_indentation_m: f64,
        columns: [Vec<f64>; 2],
        config: AdhesiveDmtMindlinConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveDmtMindlinError> {
        poll(Some(gate))?;
        validate_config(config)?;
        let materials = source
            .normal_source()
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if network.samples_rendered() != 0
            || surface_materials[0] != &materials[0]
            || surface_materials[1] != &materials[1]
            || columns
                .iter()
                .any(|c| c.len() != network.mode_count() || c.iter().any(|x| !x.is_finite()))
        {
            return Err(invalid("sample-zero DMT--Mindlin network, original ordered surfaces and complete work columns").into());
        }
        limit(
            "DMT--Mindlin sample duration",
            network.sample_period_s(),
            config.maximum_step_s,
        )?;
        let mut compliance = [[0.0; 2]; 2];
        for j in 0..2 {
            let response = network_response(&network, &columns[j], gate)?;
            for i in 0..2 {
                compliance[i][j] = dot(&columns[i], &response)?;
            }
        }
        let scale = compliance
            .iter()
            .flatten()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if scale > 0.0 {
            let c = compliance.map(|row| row.map(|x| x / scale));
            let tolerance = 128.0 * f64::EPSILON + 8.0 * network.config.solve_relative_tolerance;
            if c[0][0] < 0.0
                || c[1][1] < 0.0
                || (c[0][1] - c[1][0]).abs() > tolerance
                || c[0][0] * c[1][1] - c[0][1] * c[1][0] < -tolerance
            {
                return Err(
                    invalid("reciprocal positive-semidefinite native contact compliance").into(),
                );
            }
        }
        let state = source.initialize(initial_indentation_m)?;
        let system = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            source,
            state,
            columns,
            compliance,
            config,
            last: None,
        };
        limit(
            "native DMT--Mindlin initial storage magnitude",
            system.total_energy_j()?.abs(),
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
    }

    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub const fn source(&self) -> &'a SourceDmtMindlinContact {
        self.source
    }
    pub const fn state(&self) -> &DmtMindlinState {
        &self.state
    }
    pub const fn compliance_m_per_n(&self) -> [[f64; 2]; 2] {
        self.compliance
    }
    pub const fn work_columns(&self) -> &[Vec<f64>; 2] {
        &self.columns
    }
    pub const fn config(&self) -> AdhesiveDmtMindlinConfig {
        self.config
    }
    pub fn last_frame(&self) -> Option<&AdhesiveDmtMindlinFrame> {
        self.last.as_ref()
    }
    /// Recover the actual native bodies and retained original contact history.
    pub fn into_parts(self) -> (CoupledModalSystem, DmtMindlinState) {
        (self.network, self.state)
    }
    pub fn total_energy_j(&self) -> Result<f64, AdhesiveDmtMindlinError> {
        Ok(finite(
            self.network.total_energy_j()?
                + self
                    .source
                    .normal_source()
                    .evaluate(self.state.indentation_m())?
                    .reversible_energy_j
                + self.state.tangential_energy_j(),
        )?)
    }

    /// The explicit signed bracket selects a source-admitted root. Both bodies,
    /// all connections and the normal/shear cross response participate. Actual
    /// source and mechanical path limits apply to every constitutive proposal.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 2],
        indentation_bracket_m: [f64; 2],
    ) -> Result<&AdhesiveDmtMindlinFrame, AdhesiveDmtMindlinError> {
        self.step_under_gate(
            external,
            prescribed_rate_m_s,
            indentation_bracket_m,
            &CancelGate::new_clock_free(),
        )
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 2],
        indentation_bracket_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<&AdhesiveDmtMindlinFrame, AdhesiveDmtMindlinError> {
        let (state, frame) =
            self.prepare(external, prescribed_rate_m_s, indentation_bracket_m, gate)?;
        self.network.publish_staged();
        self.state = state;
        Ok(self.last.insert(frame))
    }

    /// Receiver refusal or late cancellation publishes neither native motion,
    /// irreversible shear memory, the clock nor the last accepted frame.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 2],
        indentation_bracket_m: [f64; 2],
        gate: &CancelGate,
        admit: impl FnOnce(&AdhesiveDmtMindlinFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveDmtMindlinError>,
    {
        let (state, frame) = self
            .prepare(external, prescribed_rate_m_s, indentation_bracket_m, gate)
            .map_err(E::from)?;
        let receiving = admit(&frame, self.network.staged_components())?;
        poll(Some(gate))
            .map_err(AdhesiveDmtMindlinError::from)
            .map_err(E::from)?;
        self.network.publish_staged();
        self.state = state;
        self.last = Some(frame);
        Ok(receiving)
    }

    fn probe(
        &self,
        indentation: f64,
        free: [f64; 2],
        gate: &CancelGate,
    ) -> Result<Probe, AdhesiveDmtMindlinError> {
        poll(Some(gate))?;
        let mean = self
            .source
            .normal_source()
            .mean_force_between(self.state.indentation_m(), indentation)?;
        limit(
            "DMT--Mindlin normal effort",
            mean.abs(),
            self.config.shear.maximum_force_n,
        )?;
        let free_shear = finite(free[1] - self.compliance[1][0] * mean)?;
        let mut cancelled = || gate.is_requested();
        let path = self.source.trial_coupled_path(
            &self.state,
            &SourceMindlinPathInput {
                indentation_m: indentation,
                free_slip_increment_m: free_shear,
                dt_s: self.network.dt,
                compliance_m_per_n: self.compliance[1][1],
                solve: self.config.shear,
            },
            &mut cancelled,
        )?;
        self.check_travel([
            indentation - self.state.indentation_m(),
            path.tangential.slip_increment_m,
        ])?;
        let residual = finite(
            (Dd::from_f64(indentation) - Dd::from_f64(free[0])
                + Dd::from_f64(self.compliance[0][0]) * Dd::from_f64(path.mean_normal_force_n)
                + Dd::from_f64(self.compliance[0][1])
                    * Dd::from_f64(path.tangential.average_force_n))
            .to_f64(),
        )?;
        let tolerance = self.coordinate_tolerance(indentation, free[0])?;
        Ok(Probe {
            indentation,
            path,
            residual,
            tolerance,
        })
    }

    fn solve(
        &self,
        free: [f64; 2],
        bracket: [f64; 2],
        gate: &CancelGate,
    ) -> Result<(Probe, usize), AdhesiveDmtMindlinError> {
        if bracket.iter().any(|x| !x.is_finite()) || bracket[0] > bracket[1] {
            return Err(invalid("finite ordered signed indentation bracket").into());
        }
        let mut a = self.probe(bracket[0], free, gate)?;
        if a.residual.abs() <= a.tolerance {
            return Ok((a, 1));
        }
        let mut b = self.probe(bracket[1], free, gate)?;
        if b.residual.abs() <= b.tolerance {
            return Ok((b, 2));
        }
        if a.residual.is_sign_positive() == b.residual.is_sign_positive() {
            return Err(AdhesiveDmtMindlinError::Bracket {
                indentation_m: bracket,
                residual_m: [a.residual, b.residual],
            });
        }
        // Bracket-preserving interpolation uses the already evaluated original
        // constitutive paths. The central fallback guarantees interval progress.
        let mut used = 2;
        while used + 1 < self.config.maximum_evaluations {
            let width = b.indentation - a.indentation;
            let weight = a.residual / (a.residual - b.residual);
            let fraction = if weight.is_finite() && (0.1..=0.9).contains(&weight) {
                weight
            } else {
                0.5
            };
            let x = a.indentation + fraction * width;
            if x == a.indentation || x == b.indentation {
                break;
            }
            let p = self.probe(x, free, gate)?;
            used += 1;
            if p.residual.abs() <= p.tolerance {
                return Ok((p, used));
            }
            if p.residual.is_sign_positive() == a.residual.is_sign_positive() {
                a = p;
            } else {
                b = p;
            }
        }
        let p = if a.residual.abs() < b.residual.abs() {
            a
        } else {
            b
        };
        Err(AdhesiveDmtMindlinError::NotConverged {
            evaluations: used,
            residual_m: p.residual,
            tolerance_m: p.tolerance,
        })
    }

    fn coordinate_tolerance(&self, a: f64, b: f64) -> Result<f64, ModalCouplingError> {
        finite(
            self.config.coordinate_absolute_tolerance_m
                + self.config.coordinate_relative_tolerance * a.abs().max(b.abs()),
        )
    }
    fn check_travel(&self, travel: [f64; 2]) -> Result<(), ModalCouplingError> {
        for x in travel {
            limit(
                "DMT--Mindlin coordinate increment",
                finite(x)?.abs(),
                self.config.maximum_coordinate_increment_m,
            )?;
            limit(
                "DMT--Mindlin coordinate speed",
                finite(x / self.network.dt)?.abs(),
                self.config.maximum_coordinate_speed_m_s,
            )?;
        }
        Ok(())
    }

    fn prepare(
        &mut self,
        external: &[f64],
        rate: [f64; 2],
        bracket: [f64; 2],
        gate: &CancelGate,
    ) -> Result<(DmtMindlinState, AdhesiveDmtMindlinFrame), AdhesiveDmtMindlinError> {
        poll(Some(gate))?;
        let old = [self.state.indentation_m(), self.state.slip_m()];
        let prescribed = [
            finite(rate[0] * self.network.dt)?,
            finite(rate[1] * self.network.dt)?,
        ];
        self.network.prepare_prediction(external, Some(gate))?;
        let free = [
            finite(old[0] + prescribed[0] + dot(&self.columns[0], &self.network.free_delta)?)?,
            finite(prescribed[1] + dot(&self.columns[1], &self.network.free_delta)?)?,
        ];
        let (solved, evaluations) = self.solve(free, bracket, gate)?;
        let held = [
            solved.path.mean_normal_force_n,
            solved.path.tangential.average_force_n,
        ];
        for (j, (force, input)) in self.forces.iter_mut().zip(external).enumerate() {
            if j % 256 == 0 {
                poll(Some(gate))?;
            }
            *force = finite(
                (Dd::from_f64(*input)
                    - Dd::from_f64(self.columns[0][j]) * Dd::from_f64(held[0])
                    - Dd::from_f64(self.columns[1][j]) * Dd::from_f64(held[1]))
                .to_f64(),
            )?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        for (j, ((change, end), start)) in self
            .network
            .free_delta
            .iter_mut()
            .zip(self.network.candidates.iter().flat_map(|m| m.states()))
            .zip(self.network.models.iter().flat_map(|m| m.states()))
            .enumerate()
        {
            if j % 256 == 0 {
                poll(Some(gate))?;
            }
            *change = finite(end.displacement_m_sqrt_kg - start.displacement_m_sqrt_kg)?;
        }
        let body = [
            dot(&self.columns[0], &self.network.free_delta)?,
            dot(&self.columns[1], &self.network.free_delta)?,
        ];
        let travel = [
            finite(prescribed[0] + body[0])?,
            finite(prescribed[1] + body[1])?,
        ];
        self.check_travel(travel)?;
        let position = [finite(old[0] + travel[0])?, finite(old[1] + travel[1])?];
        let compatibility = [
            finite(position[0] - solved.indentation)?,
            finite(travel[1] - solved.path.tangential.slip_increment_m)?,
        ];
        let tolerance = [
            self.coordinate_tolerance(position[0], solved.indentation)?,
            self.coordinate_tolerance(travel[1], solved.path.tangential.slip_increment_m)?,
        ];
        for i in 0..2 {
            if compatibility[i].abs() > tolerance[i] {
                return Err(AdhesiveDmtMindlinError::Coordinate {
                    component: i,
                    residual_m: compatibility[i],
                    tolerance_m: tolerance[i],
                });
            }
        }
        let mut cancelled = || gate.is_requested();
        let actual = self.source.advance_coupled_path(
            &self.state,
            &SourceMindlinPathInput {
                indentation_m: position[0],
                free_slip_increment_m: travel[1],
                dt_s: self.network.dt,
                compliance_m_per_n: 0.0,
                solve: self.config.shear,
            },
            &mut cancelled,
        )?;
        let material_force = [actual.mean_normal_force_n, actual.tangential.traction_n[0]];
        let mut constitutive = [0.0; 2];
        let mut force_tolerance = [0.0; 2];
        for i in 0..2 {
            constitutive[i] = finite(held[i] - material_force[i])?;
            force_tolerance[i] = finite(
                self.config.shear.force_absolute_tolerance_n
                    + self.config.shear.force_relative_tolerance
                        * held[i].abs().max(material_force[i].abs()),
            )?;
            if constitutive[i].abs() > force_tolerance[i] {
                return Err(AdhesiveDmtMindlinError::Force {
                    component: i,
                    residual_n: constitutive[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let native = &self.network.staged_frame;
        let normal_energy = actual.normal_response.reversible_energy_j;
        let shear_energy = actual.tangential.state.tangential_energy_j();
        let network_energy = finite(native.total_energy_j())?;
        limit(
            "native DMT--Mindlin storage magnitude",
            finite(network_energy + normal_energy + shear_energy)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
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
                poll(Some(gate))?;
            }
            represented =
                represented + (Dd::from_f64(*force) - Dd::from_f64(*input)) * Dd::from_f64(*dq);
        }
        let contact_work = finite(represented.to_f64())?;
        let normal_change = actual.normal_stored_change_j;
        let shear_change = actual.tangential.stored_change_j;
        let heat = actual.tangential.dissipated_j;
        let scale = normal_change
            .abs()
            .max(shear_change.abs())
            .max(heat.abs())
            .max(contact_work.abs())
            .max(prescribed_work.abs());
        let work_tolerance = finite(
            self.config.work_absolute_tolerance_j
                + self.config.work_relative_tolerance * scale
                + force_tolerance[0] * travel[0].abs()
                + force_tolerance[1] * travel[1].abs()
                + 128.0 * f64::EPSILON * scale,
        )?;
        let transpose = represented
            + Dd::from_f64(held[0]) * Dd::from_f64(body[0])
            + Dd::from_f64(held[1]) * Dd::from_f64(body[1]);
        check_work(finite(transpose.to_f64())?, work_tolerance)?;
        check_work(actual.tangential.closure_residual_j, work_tolerance)?;
        check_work(
            finite((Dd::from_f64(normal_change) - Dd::from_f64(actual.normal_work_j)).to_f64())?,
            work_tolerance,
        )?;
        let contact_residual = finite(
            (Dd::from_f64(normal_change)
                + Dd::from_f64(shear_change)
                + Dd::from_f64(heat)
                + represented
                - Dd::from_f64(prescribed_work))
            .to_f64(),
        )?;
        check_work(contact_residual, work_tolerance)?;
        // Compose the native increment account with the original contact work,
        // retaining small adhesive increments beside a large stored preload.
        let port_residual =
            Dd::from_f64(native.external_work_j) - Dd::from_f64(external_work) - represented;
        let energy_residual = finite(
            (Dd::from_f64(native.energy_residual_j)
                + Dd::from_f64(contact_residual)
                + port_residual)
                .to_f64(),
        )?;
        let energy_tolerance = finite(
            native.energy_tolerance_j
                + work_tolerance
                + 128.0 * f64::EPSILON * native.external_work_j.abs().max(external_work.abs()),
        )?;
        if energy_residual.abs() > energy_tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: energy_residual,
                tolerance_j: energy_tolerance,
            }
            .into());
        }
        let frame = AdhesiveDmtMindlinFrame {
            sample: native.sample,
            observer_pressure_pa: native.observer_pressure_pa,
            network_energy_j: network_energy,
            normal_energy_j: normal_energy,
            tangential_energy_j: shear_energy,
            normal_stored_change_j: normal_change,
            tangential_stored_change_j: shear_change,
            microslip_dissipation_j: heat,
            network_dissipation_j: native.total_dissipation_j(),
            numerical_dissipation_j: native.numerical_dissipation_j,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            contact_work_on_network_j: contact_work,
            contact_work_residual_j: contact_residual,
            contact_work_tolerance_j: work_tolerance,
            energy_residual_j: energy_residual,
            energy_tolerance_j: energy_tolerance,
            mean_force_n: held,
            position_m: position,
            endpoint_shear_force_n: actual.tangential.state.endpoint_force_n(),
            history_intervals: actual.tangential.state.history_intervals(),
            normal_response: actual.normal_response,
            contact_fraction: actual.contact_fraction,
            contact_interval_fraction: actual.contact_interval_fraction,
            compatibility_residual_m: compatibility,
            compatibility_tolerance_m: tolerance,
            constitutive_residual_n: constitutive,
            evaluations: evaluations + 1,
            shear_iterations: solved.path.tangential.iterations,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok((actual.tangential.state, frame))
    }
}

fn check_work(residual_j: f64, tolerance_j: f64) -> Result<(), AdhesiveDmtMindlinError> {
    if residual_j.abs() > tolerance_j {
        Err(AdhesiveDmtMindlinError::Work {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

fn validate_config(c: AdhesiveDmtMindlinConfig) -> Result<(), ModalCouplingError> {
    if c.maximum_evaluations < 3
        || !(1..=128).contains(&c.shear.max_iterations)
        || c.shear.force_relative_tolerance >= 1.0
        || c.shear.force_relative_tolerance <= 0.0
        || c.shear.force_absolute_tolerance_n <= 0.0
        || [
            c.maximum_step_s,
            c.maximum_coordinate_increment_m,
            c.maximum_coordinate_speed_m_s,
            c.shear.maximum_force_n,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || [
            c.coordinate_absolute_tolerance_m,
            c.coordinate_relative_tolerance,
            c.work_absolute_tolerance_j,
            c.work_relative_tolerance,
            c.shear.force_absolute_tolerance_n,
            c.shear.force_relative_tolerance,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        Err(invalid(
            "finite positive DMT--Mindlin path budgets and nonnegative accuracy gates",
        ))
    } else {
        Ok(())
    }
}
