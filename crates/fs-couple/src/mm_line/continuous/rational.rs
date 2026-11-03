//! Complete modal radiation memory coupled to the original continuous duct.
//!
//! The original Galerkin state [q; M*P] and every original MIMO pole coordinate
//! evolve together. Mouth pressure is continuous and outward duct volume flow
//! drives the impedance. No flow is held across an interval: only the incident
//! pressure waves are held. Exactly zero and reciprocal PSD mixed-rank direct
//! impedance retain their original null pressure constraints, differentiated
//! once by the existing owner without artificial resistance or state removal.
//!
//! The original source band, clock and passivity reports survive attachment.
//! Generic rational coordinates have no declared physical storage; their work
//! remains signed exchange. The independently known duct energy is not total
//! coupled energy. Spatial, radial, source-fit and clock accuracy remain
//! separate from exact evolution of this admitted finite spatial operator.

use super::{GasState, LinearSplitBudget, ModalGalerkinSource, ModalTimeError, ModeCount};
use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{TerminalSolve, eliminate_mixed_direct, eliminate_zero_direct};
use crate::modal_mouth::realization::{BemModalMouthRealization, BemModalMouthRealizationRecord};
use fs_duct::modal::time_domain::{ModalEnd, assemble_lossless, assemble_lossless_regional};
use fs_duct::{LossModel, RadiusProfile, segment_wave};
use fs_exec::CancelGate;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{DenseExponentialBudget, ScatteringPortRecord, split_quadratic};
use fs_phs::{PhsError, PortHamiltonian, PreparedStepError};
use fs_vfit::mimo::load::{MimoLoadError, PreparedMimoLoad};

/// Continuous response and finite resonances of the original joint descriptor.
#[path = "rational/analysis.rs"]
pub mod analysis;

/// Source, numerical, constraint or transactional refusal; no history is reset.
#[derive(Debug)]
pub enum RationalModalError {
    Invalid(&'static str),
    Duct(ModalTimeError),
    Linear(PhsError),
    Load(MimoLoadError),
    Terminal(AcousticRealizeError),
    Cancelled,
}
impl core::fmt::Display for RationalModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous rational modal line: {self:?}")
    }
}
impl std::error::Error for RationalModalError {}
impl From<ModalTimeError> for RationalModalError {
    fn from(value: ModalTimeError) -> Self {
        Self::Duct(value)
    }
}
impl From<PhsError> for RationalModalError {
    fn from(value: PhsError) -> Self {
        Self::Linear(value)
    }
}
impl From<MimoLoadError> for RationalModalError {
    fn from(value: MimoLoadError) -> Self {
        Self::Load(value)
    }
}
fn invalid(what: &'static str) -> RationalModalError {
    RationalModalError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), RationalModalError> {
    if gate.is_requested() {
        Err(RationalModalError::Cancelled)
    } else {
        Ok(())
    }
}

/// A proper impedance in the original signed radial pressure/volume-flow basis.
/// All ports share the same circular mouth plane; no propagation delay or
/// scalar net-flow normalization is inserted. The original load may already
/// have accepted pole history, which attachment preserves.
pub struct ModalRationalLoad {
    load: PreparedMimoLoad,
    radius_m: f64,
    exterior: GasState,
    record: Option<BemModalMouthRealizationRecord>,
}
impl ModalRationalLoad {
    /// Bind a producer-declared P=Z*U source, in Pa and m³/s. The complete port
    /// count sets the radial basis size. This declaration does not qualify a
    /// fit's physical accuracy or supply a storage matrix for its coordinates.
    pub fn new(
        load: PreparedMimoLoad,
        radius_m: f64,
        exterior: GasState,
    ) -> Result<Self, RationalModalError> {
        ModeCount(load.ports())
            .validate()
            .map_err(ModalTimeError::from)?;
        segment_wave(&exterior, radius_m, 1.0, LossModel::Lossless)
            .map_err(ModalTimeError::from)?;
        Ok(Self {
            load,
            radius_m,
            exterior,
            record: None,
        })
    }
    /// Consume the actual BEM producer, retaining its geometry, signed basis,
    /// independent samples, authored bands and original fit/admission reports.
    pub fn from_bem(realization: BemModalMouthRealization) -> Self {
        let (load, record) = realization.into_parts();
        Self {
            load,
            radius_m: record.patch().radius_m,
            exterior: *record.exterior(),
            record: Some(record),
        }
    }
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    pub const fn exterior(&self) -> &GasState {
        &self.exterior
    }
    pub fn modes(&self) -> ModeCount {
        ModeCount(self.load.ports())
    }
    pub fn load(&self) -> &PreparedMimoLoad {
        &self.load
    }
    pub fn radiation_record(&self) -> Option<&BemModalMouthRealizationRecord> {
        self.record.as_ref()
    }
}

/// Complete state, cached exact work and algebraic accuracy allowances.
#[derive(Clone, Copy, Debug)]
pub struct RationalModalLineBudget {
    /// Bounds every original duct and rational coordinate; none are discarded.
    pub linear: LinearSplitBudget,
    /// Includes the explicit held-wave coordinates and original unused input
    /// slots of the exact owner, plus every independent observation cache.
    pub exact: DenseExponentialBudget,
    /// Relative pressure/work agreement, finite and in (0, 1e-8]. No physical
    /// damping or gain is changed to meet this algebraic tolerance.
    pub relative_tolerance: f64,
}

/// One mode's independent observations of the same continuous mouth exchange.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModalTerminalAccount {
    pub mean_duct_pressure_pa: f64,
    pub mean_load_pressure_pa: f64,
    pub mean_outward_flow_m3_s: f64,
    /// Integral P_duct(t)*U_out(t) dt; positive into the rational load [J].
    pub duct_work_j: f64,
    /// Independently integrated original (C*z+D*U_out)*U_out [J].
    /// Signed exchange, not heat or an inferred rational stored energy.
    pub load_work_j: f64,
}

/// Independently measured duct storage and boundary exchanges for one interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RationalModalStep {
    pub delta_duct_energy_j: f64,
    pub inlet_work_j: f64,
    pub incident_energy_j: f64,
    pub reflected_energy_j: f64,
    pub terminal_work_j: f64,
    /// delta(H_duct) - inlet_work + terminal_work [J].
    pub balance_residual_j: f64,
}

/// One exact joint trajectory and the complete accepted physical history.
/// Preparation uses the existing exponential/Gramian owner. Each warm push
/// performs bounded matrix contractions and a private original-load candidate;
/// this path does not claim allocation-free stepping.
pub struct RationalModalLine {
    source: ModalGalerkinSource,
    duct: PortHamiltonian,
    boundary: ModalRationalLoad,
    flow: PreparedQuadraticFlow,
    // Only duct and original pole states. Held incident coordinates are NEVER
    // stored here and are not assigned physical energy or source history.
    state: Vec<f64>,
    flow_rows: Vec<Vec<f64>>,
    duct_pressure_rows: Vec<Vec<f64>>,
    load_pressure_rows: Vec<Vec<f64>>,
    constraints: Vec<Vec<f64>>,
    mean_inlet_rows: Vec<Vec<f64>>,
    inlet_accounts: Vec<ScatteringPortRecord>,
    terminal_accounts: Vec<ModalTerminalAccount>,
    reflected_means: Vec<f64>,
    tolerance: f64,
}

impl RationalModalLine {
    /// Attach every original signed mouth mode and preserve all pole history.
    /// The interval is exactly the original load clock. Inlet waves use the
    /// positive real reference rho*c/S, including evanescent radial modes.
    ///
    /// Invertible, exactly zero and reciprocal PSD mixed-rank direct D are
    /// supported. Mixed-rank flow uses the complete original positive range
    /// and differentiated null pressure constraints. Unresolved rank or singular
    /// differentiated constraints refuse. The source history must already
    /// satisfy the retained pressure manifold with the quiescent duct;
    /// inconsistent accepted history is never projected or discarded.
    pub fn new(
        profile: &RadiusProfile,
        gas: &GasState,
        elements: usize,
        boundary: ModalRationalLoad,
        budget: RationalModalLineBudget,
        gate: &CancelGate,
    ) -> Result<Self, RationalModalError> {
        Self::prepare_source(profile, gas, None, elements, boundary, budget, gate)
    }

    /// Attach the original rational mouth to frozen regional gases, with one
    /// original state per uniform axial element at a common static pressure.
    /// Temperature and composition interfaces remain on that grid; refining it
    /// requires repeating each gas over the same physical interval.
    ///
    /// Each region retains its own density, sound speed and transverse cutoff
    /// in the original modal mass and stiffness. Inlet waves use the inlet
    /// gas's rho*c/S, while the mouth keeps the original exterior source and
    /// every accepted load coordinate. No gas, impedance or history is averaged.
    /// The retained regions do not add mean flow, mixing or thermal evolution.
    ///
    /// # Errors
    /// Original regional-source refusals and the same complete-state, clock,
    /// pressure-manifold, capacity and cancellation refusals as [Self::new].
    pub fn new_regional(
        profile: &RadiusProfile,
        gases: &[GasState],
        boundary: ModalRationalLoad,
        budget: RationalModalLineBudget,
        gate: &CancelGate,
    ) -> Result<Self, RationalModalError> {
        poll(gate)?;
        let gas = gases.first().ok_or(ModalTimeError::Invalid(
            "regional modal source needs one gas state per axial element",
        ))?;
        Self::prepare_source(
            profile,
            gas,
            Some(gases),
            gases.len(),
            boundary,
            budget,
            gate,
        )
    }

    fn prepare_source(
        profile: &RadiusProfile,
        gas: &GasState,
        regional_gases: Option<&[GasState]>,
        elements: usize,
        boundary: ModalRationalLoad,
        budget: RationalModalLineBudget,
        gate: &CancelGate,
    ) -> Result<Self, RationalModalError> {
        poll(gate)?;
        let p = boundary.load.ports();
        let radius = profile.outlet_radius();
        if (boundary.radius_m - radius).abs() > 128.0 * f64::EPSILON * radius
            || boundary.exterior.pressure.to_bits() != gas.pressure.to_bits()
            || !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance <= 0.0
            || budget.relative_tolerance > 1e-8
        {
            return Err(invalid(
                "mouth radius/static pressure or explicit algebraic tolerance mismatch",
            ));
        }
        let view = boundary.load.continuous_operator();
        let nt = view.dimension();
        let nd = elements
            .checked_add(1)
            .and_then(|v| v.checked_mul(2 * p))
            .ok_or_else(|| invalid("complete modal spatial dimension overflow"))?;
        let physical = nd
            .checked_add(nt)
            .ok_or_else(|| invalid("complete modal state overflow"))?;
        let n = physical
            .checked_add(p)
            .ok_or_else(|| invalid("held-wave dimension overflow"))?;
        let d = n
            .checked_add(2 * p + 1)
            .ok_or_else(|| invalid("exact-flow dimension overflow"))?;
        if physical > budget.linear.max_modes || d > 256 || d > budget.exact.max_dimension {
            return Err(invalid(
                "complete duct, load and held waves exceed authored exact-flow dimension",
            ));
        }
        // Nine independent forms per mode, full original source matrices and
        // existing owner's conservative Gramian workspace; admit before forms.
        let bytes = d as u128 * d as u128 * 8 * (18 * p + 48) as u128;
        let work = (d as u128).pow(3) * (224 * 9 * p + 512) as u128;
        if bytes > budget.exact.max_dense_bytes as u128
            || work > budget.exact.max_multiply_terms as u128
        {
            return Err(invalid(
                "complete duct, load, held waves and observations exceed authored exact-flow capacity",
            ));
        }
        let (source, duct) = match regional_gases {
            Some(gases) => assemble_lossless_regional(
                profile,
                gases,
                ModeCount(p),
                ModalEnd::Ports,
                budget.linear.max_modes,
                || gate.is_requested(),
            ),
            None => assemble_lossless(
                profile,
                gas,
                ModeCount(p),
                elements,
                ModalEnd::Ports,
                budget.linear.max_modes,
                || gate.is_requested(),
            ),
        }?
        .into_parts();
        let split = split_quadratic(&duct, budget.linear)?;
        let model = view.state_space(budget.exact.max_dense_bytes, &mut || gate.is_requested())?;
        let q = split.quadratic_matrix();
        let (j, r, g) = duct.structure();
        let mut generator = vec![0.0; n * n];
        let mut inlet_pressure = vec![vec![0.0; n]; p];
        let mut duct_pressure_rows = vec![vec![0.0; n]; p];
        let mut inputs = vec![vec![0.0; n]; p];
        let mut constraints = vec![vec![0.0; n]; p];
        let z0 = source.inlet_reference_impedance();
        for i in 0..nd {
            poll(gate)?;
            for col in 0..nd {
                generator[i * n + col] = (0..nd)
                    .map(|k| (j[i * nd + k] - r[i * nd + k]) * q[k * nd + col])
                    .sum();
            }
            for k in 0..p {
                inlet_pressure[k][i] = (0..nd).map(|a| g[a * 2 * p + k] * q[a * nd + i]).sum();
                duct_pressure_rows[k][i] =
                    (0..nd).map(|a| g[a * 2 * p + p + k] * q[a * nd + i]).sum();
                inputs[k][i] = -g[i * 2 * p + p + k]; // outward mouth flow
            }
        }
        for i in 0..nd {
            for k in 0..p {
                for col in 0..nd {
                    generator[i * n + col] -= g[i * 2 * p + k] * inlet_pressure[k][col] / z0;
                }
                // Explicit constant coordinates hold ALL incident modes. This
                // lets constrained elimination react to the complete forcing at once.
                generator[i * n + physical + k] = 2.0 * g[i * 2 * p + k] / z0;
            }
        }
        for i in 0..nt {
            generator[(nd + i) * n + nd..(nd + i) * n + physical]
                .copy_from_slice(&model.a[i * nt..(i + 1) * nt]);
            for k in 0..p {
                inputs[k][nd + i] = model.b[i * p + k];
            }
        }
        for k in 0..p {
            constraints[k].copy_from_slice(&duct_pressure_rows[k]);
            for i in 0..nt {
                constraints[k][nd + i] = -model.c[k * nt + i];
            }
        }
        let flow_rows = if model.d.iter().all(|v| *v == 0.0) {
            let original =
                eliminate_zero_direct(&generator, &inputs, &constraints, &vec![0.0; n], gate)
                    .map_err(RationalModalError::Terminal)?;
            generator = original.generator;
            constraints = original.constraints;
            original.flow_columns
        } else if let Ok(solve) = TerminalSolve::new(&model.d, p) {
            let mut rows = vec![vec![0.0; n]; p];
            for col in 0..n {
                poll(gate)?;
                let rhs: Vec<_> = constraints.iter().map(|row| row[col]).collect();
                for (k, value) in solve
                    .apply(&rhs)
                    .map_err(RationalModalError::Terminal)?
                    .into_iter()
                    .enumerate()
                {
                    rows[k][col] = value;
                    for i in 0..n {
                        generator[i * n + col] += inputs[k][i] * value;
                    }
                }
            }
            constraints.clear();
            rows
        } else {
            let mixed_bytes = bytes + 256u128 * (n as u128).pow(2);
            if mixed_bytes > budget.exact.max_dense_bytes as u128 {
                return Err(invalid(
                    "original mixed-rank terminal workspace exceeds authored capacity",
                ));
            }
            let original = eliminate_mixed_direct(
                &generator,
                &inputs,
                &constraints,
                &model.d,
                &vec![0.0; n],
                gate,
            )
            .map_err(RationalModalError::Terminal)?;
            generator = original.generator;
            constraints = original.constraints;
            original.flow_columns
        };
        let mut load_pressure_rows = vec![vec![0.0; n]; p];
        for k in 0..p {
            for col in 0..n {
                load_pressure_rows[k][col] = (0..p)
                    .map(|a| model.d[k * p + a] * flow_rows[a][col])
                    .sum::<f64>()
                    + if (nd..physical).contains(&col) {
                        model.c[k * nt + col - nd]
                    } else {
                        0.0
                    };
            }
        }
        let dt = boundary.load.time_step_s();
        let mut extra_scales = Vec::with_capacity(nt + p);
        let mut at = 0;
        for term in &boundary.load.source_model().poles {
            let width = term.state_dim() * p;
            let norm = (0..p)
                .flat_map(|row| model.c[row * nt + at..row * nt + at + width].iter())
                .map(|v| v.abs())
                .fold(0.0, f64::max);
            extra_scales.extend(core::iter::repeat_n(
                if norm == 0.0 {
                    1.0
                } else {
                    fs_math::det::sqrt(norm)
                },
                width,
            ));
            at += width;
        }
        // Numerical conditioning of constant forcing, not additional storage.
        extra_scales.extend(core::iter::repeat_n(fs_math::det::sqrt(dt / z0), p));
        let mut quadratic = Vec::with_capacity(5 * p);
        let mut linear = Vec::with_capacity(4 * p);
        for k in 0..p {
            let mut incident = vec![0.0; n];
            incident[physical + k] = 1.0;
            let mut inlet_flow: Vec<_> = inlet_pressure[k].iter().map(|v| -v / z0).collect();
            inlet_flow[physical + k] += 2.0 / z0;
            let mut outgoing = inlet_pressure[k].clone();
            outgoing[physical + k] -= 1.0;
            quadratic.push(product_form(&inlet_pressure[k], &inlet_flow, 1.0));
            quadratic.push(product_form(&incident, &incident, 1.0 / z0));
            quadratic.push(product_form(&outgoing, &outgoing, 1.0 / z0));
            quadratic.push(product_form(&duct_pressure_rows[k], &flow_rows[k], 1.0));
            quadratic.push(product_form(&load_pressure_rows[k], &flow_rows[k], 1.0));
            linear.extend([
                inlet_pressure[k].clone(),
                duct_pressure_rows[k].clone(),
                load_pressure_rows[k].clone(),
                flow_rows[k].clone(),
            ]);
        }
        let mut state = vec![0.0; physical];
        state[nd..].copy_from_slice(&view.real_coordinates());
        check_constraints(&constraints, &state, budget.relative_tolerance)?;
        let flow = PreparedQuadraticFlow::new_extended_controlled(
            &duct,
            &split,
            &generator,
            &extra_scales,
            dt,
            &quadratic,
            &linear,
            budget.exact,
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            PreparedStepError::Cancelled => RationalModalError::Cancelled,
            PreparedStepError::Solver(error) => RationalModalError::Linear(error),
        })?;
        let mut mean_inlet_rows = Vec::with_capacity(p);
        for k in 0..p {
            let (mut row, _) = flow.linear_integral_map(5 * p + 4 * k)?;
            for value in &mut row {
                *value /= dt;
            }
            if row.iter().any(|value| !value.is_finite()) {
                return Err(invalid("modal mean reflection map is not representable"));
            }
            mean_inlet_rows.push(row);
        }
        poll(gate)?;
        Ok(Self {
            source,
            duct,
            boundary,
            flow,
            state,
            flow_rows,
            duct_pressure_rows,
            load_pressure_rows,
            constraints,
            mean_inlet_rows,
            inlet_accounts: vec![ScatteringPortRecord::default(); p],
            terminal_accounts: vec![ModalTerminalAccount::default(); p],
            reflected_means: vec![0.0; p],
            tolerance: budget.relative_tolerance,
        })
    }

    pub fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }
    pub fn boundary(&self) -> &ModalRationalLoad {
        &self.boundary
    }
    pub fn load(&self) -> &PreparedMimoLoad {
        &self.boundary.load
    }
    pub fn duration(&self) -> f64 {
        self.boundary.load.time_step_s()
    }
    /// Accepted original [duct q; duct M*P; all original real pole coordinates].
    /// Held incident waves are inputs, never physical history in this slice.
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    pub fn duct_state(&self) -> &[f64] {
        &self.state[..self.duct.state_dim()]
    }
    pub fn load_state(&self) -> &[f64] {
        &self.state[self.duct.state_dim()..]
    }
    /// Original duct storage only; no rational storage matrix is inferred.
    pub fn stored_duct_energy_j(&self) -> f64 {
        self.duct.hamiltonian(self.duct_state())
    }
    pub fn port_count(&self) -> usize {
        self.boundary.load.ports()
    }
    /// Accepted interval-mean reflected pressure waves, distinct from energy.
    pub fn incoming(&self) -> &[f64] {
        &self.reflected_means
    }
    pub fn port_accounts(&self) -> &[ScatteringPortRecord] {
        &self.inlet_accounts
    }
    pub fn terminal_accounts(&self) -> &[ModalTerminalAccount] {
        &self.terminal_accounts
    }

    /// Same-interval b_bar=free+direct*a in every original inlet mode. The native
    /// cached integral map supplies the complete mutual relation; querying it
    /// changes no history. Both caller buffers are unchanged on refusal.
    pub fn reflection_relation_into(
        &self,
        free: &mut [f64],
        direct: &mut [f64],
    ) -> Result<(), RationalModalError> {
        let p = self.port_count();
        let n = self.state.len();
        if free.len() != p || direct.len() != p * p {
            return Err(invalid("modal reflection output shapes"));
        }
        let candidate: Vec<_> = self
            .mean_inlet_rows
            .iter()
            .map(|row| dot(row, &self.state))
            .collect();
        if candidate.iter().any(|v| !v.is_finite()) {
            return Err(invalid("modal mean reflection overflow"));
        }
        for k in 0..p {
            for j in 0..p {
                direct[k * p + j] = self.mean_inlet_rows[k][n + j] - if k == j { 1.0 } else { 0.0 };
            }
        }
        free.copy_from_slice(&candidate);
        Ok(())
    }

    /// Evolve one original-clock interval with held incident pressure waves.
    /// All mouth flows, pressures and pole coordinates vary jointly within it.
    /// Pressure continuity, two independent mouth work forms and original duct
    /// storage balance are checked before publishing any state or output.
    pub fn push(
        &mut self,
        incident: &[f64],
        gate: &CancelGate,
    ) -> Result<RationalModalStep, RationalModalError> {
        poll(gate)?;
        let p = self.port_count();
        let physical = self.state.len();
        let nd = self.duct.state_dim();
        if incident.len() != p || incident.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "one finite incident wave is required per original inlet mode",
            ));
        }
        let mut start = self.state.clone();
        start.extend_from_slice(incident);
        check_constraints(&self.constraints, &start, self.tolerance)?;
        let next = self.flow.evaluate(&start, &vec![0.0; 2 * p])?;
        check_constraints(&self.constraints, &next.state, self.tolerance)?;
        for k in 0..p {
            let a = &self.duct_pressure_rows[k];
            let b = &self.load_pressure_rows[k];
            let scale = a
                .iter()
                .zip(b)
                .zip(&next.state)
                .map(|((a, b), x)| (a * x).abs() + (b * x).abs())
                .sum();
            check_equal(
                dot(a, &next.state),
                dot(b, &next.state),
                scale,
                self.tolerance,
                "rational modal endpoint pressure continuity",
            )?;
        }
        let dt = self.duration();
        let before = self.stored_duct_energy_j();
        let after = self.duct.hamiltonian(&next.state[..nd]);
        let mut record = RationalModalStep {
            delta_duct_energy_j: after - before,
            ..RationalModalStep::default()
        };
        let mut ports = Vec::with_capacity(p);
        let mut terminal = Vec::with_capacity(p);
        for k in 0..p {
            let q = &next.integrals[5 * k..5 * k + 5];
            let l = &next.integrals[5 * p + 4 * k..5 * p + 4 * k + 4];
            let scale = before.abs() + after.abs() + q.iter().map(|v| v.abs()).sum::<f64>();
            check_equal(
                q[0],
                q[1] - q[2],
                scale,
                self.tolerance,
                "modal incident/reflected work mismatch",
            )?;
            check_equal(
                q[3],
                q[4],
                scale,
                self.tolerance,
                "modal duct/load work mismatch",
            )?;
            if q[1] < -self.tolerance * scale || q[2] < -self.tolerance * scale {
                return Err(invalid("modal wave energy is materially negative"));
            }
            // Original signed values survive roundoff checks; never clip energy.
            ports.push(ScatteringPortRecord {
                boundary_work: q[0],
                incident_energy: q[1],
                outgoing_energy: q[2],
                outgoing_wave: l[0] / dt - incident[k],
            });
            terminal.push(ModalTerminalAccount {
                mean_duct_pressure_pa: l[1] / dt,
                mean_load_pressure_pa: l[2] / dt,
                mean_outward_flow_m3_s: l[3] / dt,
                duct_work_j: q[3],
                load_work_j: q[4],
            });
            record.inlet_work_j += q[0];
            record.incident_energy_j += q[1];
            record.reflected_energy_j += q[2];
            record.terminal_work_j += q[3];
        }
        record.balance_residual_j =
            record.delta_duct_energy_j - record.inlet_work_j + record.terminal_work_j;
        let scale = before.abs()
            + after.abs()
            + record.incident_energy_j.abs()
            + record.reflected_energy_j.abs()
            + terminal.iter().map(|v| v.duct_work_j.abs()).sum::<f64>();
        check_equal(
            record.delta_duct_energy_j,
            record.inlet_work_j - record.terminal_work_j,
            scale,
            self.tolerance,
            "original modal duct energy/port work mismatch",
        )?;
        let endpoint_input: Vec<_> = self
            .flow_rows
            .iter()
            .map(|row| dot(row, &next.state))
            .collect();
        let load = self
            .boundary
            .load
            .continuous_operator()
            .candidate_from_coordinates(&next.state[nd..physical], &endpoint_input)?;
        if ports.iter().any(|v| !v.outgoing_wave.is_finite())
            || terminal.iter().any(|v| {
                !v.mean_duct_pressure_pa.is_finite()
                    || !v.mean_load_pressure_pa.is_finite()
                    || !v.mean_outward_flow_m3_s.is_finite()
            })
        {
            return Err(invalid("modal interval mean overflows"));
        }
        poll(gate)?;
        self.state.copy_from_slice(&next.state[..physical]);
        self.boundary.load = load;
        for (out, port) in self.reflected_means.iter_mut().zip(&ports) {
            *out = port.outgoing_wave;
        }
        self.inlet_accounts = ports;
        self.terminal_accounts = terminal;
        Ok(record)
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn product_form(a: &[f64], b: &[f64], factor: f64) -> Vec<f64> {
    let n = a.len();
    let mut form = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            form[i * n + j] = 0.5 * factor * (a[i] * b[j] + b[i] * a[j]);
        }
    }
    form
}
fn check_equal(
    a: f64,
    b: f64,
    scale: f64,
    tolerance: f64,
    what: &'static str,
) -> Result<(), RationalModalError> {
    if !a.is_finite()
        || !b.is_finite()
        || !scale.is_finite()
        || (a - b).abs() > tolerance * scale.max(f64::MIN_POSITIVE)
    {
        Err(invalid(what))
    } else {
        Ok(())
    }
}
fn check_constraints(
    rows: &[Vec<f64>],
    state: &[f64],
    tolerance: f64,
) -> Result<(), RationalModalError> {
    for row in rows {
        let scale = row.iter().zip(state).map(|(h, x)| (h * x).abs()).sum();
        check_equal(
            dot(row, state),
            0.0,
            scale,
            tolerance,
            "original modal null-pressure constraint is inconsistent; state projection is forbidden",
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
