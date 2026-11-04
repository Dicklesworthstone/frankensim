//! Developed Newtonian sections with one exactly common incompressible flow.
//!
//! Each original pipe retains its geometry, source, viscosity and every radial
//! diffusion arm. Only the duplicate MEAN flows are eliminated analytically:
//! `L=sum(L_i)`, `L Q' = dp-sum(R_i)Q-sum(L_i*g_ij*(Q-z_ij))`,
//! `z'_ij=p_ij*(Q-z_ij)`. The energy is the sum of the original pipe energies.
//! This is the exact series connection of their admitted rational Stokes laws,
//! with no intermediate pressure storage or time-discrete continuity repair.
//!
//! Every section's irreversible heat is integrated from its ORIGINAL positive
//! viscous source squares along that same complete coupled trajectory. It is
//! not a share of total heat, an endpoint power, or memory-pressure work. Thus
//! a thermal caller can deposit heat into the correct material inventory.
//!
//! Sections share exactly one fixed reference density. Their different radii
//! describe ideal developed series components; entrance, area-change, convective
//! momentum and junction losses require their own owner. Density expansion,
//! radial temperature gradients, mixing, turbulence and thermal evolution are
//! not inferred. All single-pipe frequency-band and Reynolds limits still apply.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fluid_transport::FluidTransport;
use fs_math::{c64::C64, det};
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh, PreparedSourceLoss,
};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage, ResistanceWitness};
use fs_qty::{Pressure, VolumetricFlowRate};

use super::{
    GraphError, LaminarPipe, LaminarPipeError, LaminarPipeState, LaminarPipeStep, checkpoint,
    energy, positive,
};

/// Complete source-loss integration retains `2*n` coordinates in the original
/// native 256-coordinate bound. No material-memory arm is truncated to fit.
pub const MAX_LAMINAR_SERIES_STATES: usize = 128;

/// The existing energy-normalized profile state, bound to THIS series owner.
/// Its single mean coordinate stores the sum of the section mean energies;
/// subsequent coordinates retain each section's original profile arms.
pub type LaminarSeriesState = LaminarPipeState;

/// Existing whole-path flow/energy record. `reynolds_envelope` is the largest
/// section Reynolds number; every section's own ceiling is admitted separately.
pub type LaminarSeriesStep = LaminarPipeStep;

/// Original source-bearing developed sections and their exact common-flow law.
pub struct LaminarSeries<M> {
    sections: Vec<LaminarPipe<M>>,
    offsets: Vec<usize>,
    inertance: f64,
    steady_resistance: f64,
    flow_scales: Vec<f64>,
    source_rows: Vec<Vec<f64>>,
    system: PortHamiltonian,
    identity: Arc<()>,
}

impl<M: FluidTransport> LaminarSeries<M> {
    /// Connect admitted sections in their authored positive direction. Source
    /// snapshots are retained without resampling. The complete state count and
    /// exact density agreement are checked before dense allocation.
    pub fn try_new(cx: &Cx<'_>, sections: Vec<LaminarPipe<M>>) -> Result<Self, LaminarPipeError> {
        checkpoint(cx)?;
        let Some(first) = sections.first() else {
            return Err(LaminarPipeError::Invalid(
                "a laminar series needs a section",
            ));
        };
        let density = first.transport.density_kg_m3();
        let mut n = 1_usize;
        let mut inertance = 0.0;
        let mut steady_resistance = 0.0;
        for pipe in &sections {
            checkpoint(cx)?;
            n = n.checked_add(pipe.memory.series_terms().len()).ok_or(
                LaminarPipeError::Invalid("laminar series state extent overflow"),
            )?;
            if n > MAX_LAMINAR_SERIES_STATES || pipe.transport.density_kg_m3() != density {
                return Err(LaminarPipeError::Invalid(
                    "bounded complete memory and exactly common reference density required",
                ));
            }
            inertance += pipe.inertance;
            steady_resistance += pipe.steady_resistance;
        }
        let mean_scale = det::sqrt(inertance);
        let mean_drag = steady_resistance / inertance;
        if ![
            inertance,
            steady_resistance,
            mean_scale,
            1.0 / mean_scale,
            mean_drag,
        ]
        .into_iter()
        .all(positive)
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable series inertia or drag",
            ));
        }
        let mut resistance = vec![0.0; n * n];
        let mut squares = vec![0.0; n * n];
        let mut q = vec![0.0; n * n];
        let mut g = vec![0.0; 2 * n];
        let mut flow_scales = Vec::with_capacity(n);
        let mut offsets = Vec::with_capacity(sections.len() + 1);
        let mut source_rows = Vec::with_capacity(sections.len());
        flow_scales.push(mean_scale);
        resistance[0] = mean_drag;
        squares[0] = det::sqrt(mean_drag);
        g[0] = 1.0 / mean_scale;
        g[1] = -g[0];
        q[0] = 1.0;
        for pipe in &sections {
            checkpoint(cx)?;
            let start = flow_scales.len();
            offsets.push(start);
            let weight = pipe.inertance / inertance;
            let rows = 1 + pipe.memory.series_terms().len();
            let mut local = vec![0.0; rows * n];
            local[0] = det::sqrt(pipe.steady_resistance / inertance);
            if !positive(weight) || !positive(local[0]) {
                return Err(LaminarPipeError::Invalid(
                    "unrepresentable section resistance share",
                ));
            }
            for (j, term) in pipe.memory.series_terms().iter().enumerate() {
                let k = flow_scales.len();
                let gain = weight * term.gain_per_s;
                let a = det::sqrt(gain);
                let b = det::sqrt(term.rate_per_s);
                let coupling = a * b;
                let scale = pipe.flow_scales[j + 1];
                if ![gain, a, b, coupling, scale, 1.0 / scale]
                    .into_iter()
                    .all(positive)
                {
                    return Err(LaminarPipeError::Invalid(
                        "unrepresentable series profile scale",
                    ));
                }
                resistance[0] += gain;
                resistance[k] = -coupling;
                resistance[k * n] = -coupling;
                resistance[k * n + k] = term.rate_per_s;
                squares[k * n] = a;
                squares[k * n + k] = -b;
                local[(j + 1) * n] = a;
                local[(j + 1) * n + k] = -b;
                q[k * n + k] = 1.0;
                flow_scales.push(scale);
            }
            source_rows.push(local);
        }
        offsets.push(n);
        let witness = ResistanceWitness::try_from_source_rows(n, resistance.clone(), squares)?;
        let system = PortHamiltonian::new(
            n,
            2,
            vec![0.0; n * n],
            resistance,
            g,
            Box::new(QuadraticStorage::new(q, n)?),
        )?
        .with_resistance_witness(witness)?;
        checkpoint(cx)?;
        Ok(Self {
            sections,
            offsets,
            inertance,
            steady_resistance,
            flow_scales,
            source_rows,
            system,
            identity: Arc::new(()),
        })
    }

    /// Continue every original section under its own new viscosity. Each source
    /// is sampled once by the original pipe; density, geometry and arm order
    /// remain fixed, and every band admission is renewed. The accepted coupled
    /// profile and clock remain compatible, with no reset or injected energy.
    pub fn with_sources<N: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        sources: Vec<N>,
    ) -> Result<LaminarSeries<N>, LaminarPipeError> {
        checkpoint(cx)?;
        if sources.len() != self.sections.len() {
            return Err(LaminarPipeError::Invalid(
                "one new source per original section required",
            ));
        }
        let sections = self
            .sections
            .iter()
            .zip(sources)
            .map(|(pipe, source)| pipe.with_source(cx, source))
            .collect::<Result<Vec<_>, _>>()?;
        let mut next = LaminarSeries::try_new(cx, sections)?;
        next.identity = Arc::clone(&self.identity);
        Ok(next)
    }

    pub fn sections(&self) -> &[LaminarPipe<M>] {
        &self.sections
    }
    pub const fn mean_inertance_pa_s2_m3(&self) -> f64 {
        self.inertance
    }
    pub const fn steady_resistance_pa_s_m3(&self) -> f64 {
        self.steady_resistance
    }
    pub const fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Sum of the original complex pressure/flow impedances, with every local
    /// frequency band checked by the original wall-memory source.
    pub fn impedance_at(&self, omega_rad_s: f64) -> Result<C64, LaminarPipeError> {
        let mut impedance = C64::new(0.0, 0.0);
        for section in &self.sections {
            impedance = impedance + section.impedance_at(omega_rad_s)?;
        }
        if !positive(impedance.re) || !positive(impedance.im) {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable series impedance",
            ));
        }
        Ok(impedance)
    }

    /// Common steady profile; each section includes its parabolic-profile energy.
    pub fn steady_state(
        &self,
        volume_flow: VolumetricFlowRate,
        accepted_time_s: f64,
    ) -> Result<LaminarSeriesState, LaminarPipeError> {
        if !accepted_time_s.is_finite() || accepted_time_s < 0.0 {
            return Err(LaminarPipeError::Invalid(
                "finite nonnegative series clock required",
            ));
        }
        self.admit_flow(volume_flow.value().abs())?;
        let coordinates: Vec<_> = self
            .flow_scales
            .iter()
            .map(|s| s * volume_flow.value())
            .collect();
        if coordinates.iter().any(|x| !x.is_finite())
            || !energy(&coordinates).is_finite()
            || (volume_flow.value() != 0.0
                && (coordinates.iter().any(|x| *x == 0.0) || energy(&coordinates) == 0.0))
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable initial series state",
            ));
        }
        Ok(LaminarPipeState {
            coordinates,
            mean_flow_scale: self.flow_scales[0],
            accepted_time_s,
            identity: Arc::clone(&self.identity),
        })
    }

    /// Physical kinetic inventory of one section, including its profile memory.
    pub fn section_kinetic_energy_j(
        &self,
        state: &LaminarSeriesState,
        section: usize,
    ) -> Result<f64, LaminarPipeError> {
        self.admit_owner(state)?;
        let pipe = self
            .sections
            .get(section)
            .ok_or(LaminarPipeError::Invalid("laminar series section index"))?;
        Ok(0.5 * pipe.inertance * state.volume_flow_m3_s().powi(2)
            + energy(&state.coordinates[self.offsets[section]..self.offsets[section + 1]]))
    }

    /// Cold preparation of the exact complete trajectory and all original local
    /// heat integrals. Warm `advance` reuses every buffer and allocates nothing.
    pub fn prepare(
        &self,
        cx: &Cx<'_>,
        duration_s: f64,
    ) -> Result<PreparedLaminarSeries<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        if !positive(duration_s) {
            return Err(LaminarPipeError::Invalid(
                "positive finite series interval required",
            ));
        }
        let n = self.flow_scales.len();
        let native = PreparedExactZoh::new(
            &self.system,
            duration_s,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: 1.0,
            },
        )?;
        // These are energy coordinates: Q=I, J=0, so the exact physical
        // generator is -R. Each original section's source squares observe it.
        let (_, r, g) = self.system.structure();
        let generator: Vec<_> = r.iter().map(|r| -r).collect();
        let mut heat = Vec::with_capacity(self.sections.len());
        for rows in &self.source_rows {
            heat.push(
                PreparedSourceLoss::new_controlled(
                    n,
                    2,
                    &generator,
                    g,
                    rows,
                    duration_s,
                    DenseExponentialBudget {
                        max_dimension: 2 * n,
                        max_dense_bytes: usize::MAX,
                        max_multiply_terms: u64::MAX,
                    },
                    || cx.checkpoint().is_err(),
                )
                .map_err(prepared_error)?,
            );
        }
        checkpoint(cx)?;
        Ok(PreparedLaminarSeries {
            series: self,
            native,
            heat,
            pending: vec![0.0; n],
            heat_scratch: vec![0.0; 2 * n],
            pending_heat: vec![0.0; self.sections.len()],
        })
    }

    fn admit_owner(&self, state: &LaminarSeriesState) -> Result<(), LaminarPipeError> {
        if !Arc::ptr_eq(&state.identity, &self.identity) {
            return Err(LaminarPipeError::Invalid(
                "viscous history belongs to a different series",
            ));
        }
        Ok(())
    }
    fn admit_flow(&self, maximum: f64) -> Result<f64, LaminarPipeError> {
        let mut reynolds = 0.0_f64;
        for pipe in &self.sections {
            reynolds = reynolds.max(pipe.admit_flow(maximum)?);
        }
        Ok(reynolds)
    }
}

/// Fixed-clock original native trajectory with one heat cache per section.
pub struct PreparedLaminarSeries<'a, M> {
    series: &'a LaminarSeries<M>,
    native: PreparedExactZoh<'a>,
    heat: Vec<PreparedSourceLoss>,
    pending: Vec<f64>,
    heat_scratch: Vec<f64>,
    pending_heat: Vec<f64>,
}
impl<M: FluidTransport> PreparedLaminarSeries<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }

    /// Actual step-average common flow as `history + mobility*(p_from-p_to)`.
    pub fn average_flow_relation(
        &self,
        state: &LaminarSeriesState,
    ) -> Result<(f64, f64), LaminarPipeError> {
        self.series.admit_owner(state)?;
        let relation = self
            .native
            .average_relation(&state.coordinates, &[0.0, 0.0], 0)?;
        if !positive(relation.1) {
            return Err(LaminarPipeError::Invalid("unrepresentable series mobility"));
        }
        Ok(relation)
    }

    /// Advance the COMMON physical trajectory, publishing state and local heat
    /// together. `section_heat_j` must retain one slot per authored section.
    /// All outputs remain unchanged on source/regime, arithmetic or cancellation
    /// refusal. The cooperative Q/z comparison bound checks every section's
    /// Reynolds support throughout this held-pressure interval, even in reversal.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarSeriesState,
        pressures: [Pressure; 2],
        section_heat_j: &mut [f64],
    ) -> Result<LaminarSeriesStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.series.admit_owner(state)?;
        let duration_s = self.native.duration();
        let next_time = state.accepted_time_s + duration_s;
        let pressures = pressures.map(|p| p.value());
        let dp = pressures[0] - pressures[1];
        if section_heat_j.len() != self.heat.len()
            || pressures.iter().any(|p| !p.is_finite())
            || !dp.is_finite()
            || !next_time.is_finite()
            || next_time <= state.accepted_time_s
        {
            return Err(LaminarPipeError::Invalid(
                "complete heat slots, finite pressure and advancing clock required",
            ));
        }
        let initial_max = state
            .coordinates
            .iter()
            .zip(&self.series.flow_scales)
            .map(|(x, s)| (x / s).abs())
            .fold(0.0_f64, f64::max);
        let maximum = initial_max
            .max(dp.abs() / self.series.steady_resistance)
            .min(initial_max + (dp.abs() / self.series.inertance) * duration_s);
        if !maximum.is_finite() || (initial_max == 0.0 && dp != 0.0 && maximum == 0.0) {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable series flow envelope",
            ));
        }
        let reynolds_envelope = self.series.admit_flow(maximum)?;
        let input = [dp, 0.0];
        let mut output = [0.0; 2];
        let native =
            self.native
                .step_into(&state.coordinates, &input, &mut self.pending, &mut output)?;
        for (i, heat) in self.heat.iter().enumerate() {
            self.pending_heat[i] = heat
                .evaluate_into_controlled(
                    &state.coordinates,
                    &input,
                    &mut self.heat_scratch,
                    || cx.checkpoint().is_err(),
                )
                .map_err(prepared_error)?;
        }
        let sum: f64 = self.pending_heat.iter().sum();
        let allowance = 4096.0
            * f64::EPSILON
            * self.pending.len() as f64
            * sum.abs().max(native.dissipated.abs());
        if !sum.is_finite() || (sum - native.dissipated).abs() > allowance {
            return Err(LaminarPipeError::Invalid(
                "original section heat does not close the native viscous account",
            ));
        }
        let final_volume_flow_m3_s = self.pending[0] / self.series.flow_scales[0];
        let final_mean_kinetic_energy_j = 0.5 * self.pending[0].powi(2);
        let final_profile_memory_energy_j = energy(&self.pending[1..]);
        let displaced_volume_m3 = duration_s * output[0];
        if ![
            final_volume_flow_m3_s,
            final_mean_kinetic_energy_j,
            final_profile_memory_energy_j,
            displaced_volume_m3,
        ]
        .into_iter()
        .all(f64::is_finite)
            || (output[0] != 0.0 && displaced_volume_m3 == 0.0)
            || output[1] != -output[0]
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable series flow, energy or continuity",
            ));
        }
        let report = LaminarPipeStep {
            duration_s,
            initial_volume_flow_m3_s: state.volume_flow_m3_s(),
            final_volume_flow_m3_s,
            average_volume_flow_m3_s: output[0],
            displaced_volume_m3,
            reynolds_envelope,
            initial_mean_kinetic_energy_j: state.mean_kinetic_energy_j(),
            final_mean_kinetic_energy_j,
            initial_profile_memory_energy_j: state.profile_memory_energy_j(),
            final_profile_memory_energy_j,
            native,
        };
        checkpoint(cx)?;
        state.coordinates.copy_from_slice(&self.pending);
        state.accepted_time_s = next_time;
        section_heat_j.copy_from_slice(&self.pending_heat);
        Ok(report)
    }
}

fn prepared_error(error: PreparedStepError) -> LaminarPipeError {
    match error {
        PreparedStepError::Solver(error) => LaminarPipeError::Native(error),
        PreparedStepError::Cancelled => LaminarPipeError::Execution(GraphError::Interrupted),
    }
}

#[cfg(test)]
mod tests;
