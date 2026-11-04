//! Free modes of the complete frozen modal/material/DG/terminal assembly.
//!
//! A real polynomial pencil retains every original Maxwell coordinate, duct
//! memory, finite branch, terminal realization, gap cell and algebraic flow.
//! Only mechanical velocities are represented analytically as v=s*q. The
//! native nonlinear-eigenvalue owner never inverts a singular descriptor mass,
//! a terminal direct matrix or an isolated acoustic impedance.
use super::*;
use core::ops::Range;
pub use fs_la::nonlinear_eig::PoleOptions;
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, NonlinearPole, newton_pole_with_cancel,
};

mod pencil;
use pencil::Pencil;
mod spectrum;
pub use spectrum::{BeynOptions, DgLoadedModalSpectrum, DgLoadedSpectrumError, ResonanceRegion};
mod response;
pub use response::{
    DgLoadedModalForceResidue, DgLoadedModalPoleDerivative, DgPoleParameterDirection,
    DgPoleResponseError, DgPoleResponseOptions, PoleResponseOptions, PreparedDgPoleResponse,
};

/// Cold complete-pencil resources, excluding already owned physical sources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgPolePreparation {
    /// Positive numerical row-equilibration scale [1/s], not a physical loss.
    pub reference_rate_per_s: f64,
    pub maximum_dimension: usize,
    pub maximum_payload_bytes: usize,
    pub maximum_work: u64,
}

/// Original finite gap coordinates, or an empty range for a rigid closure.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPoleChannelCoordinates {
    pub port: PortId,
    pub velocity_cells: Range<usize>,
    pub pressure_drop: Option<usize>,
}

/// Exact state ordering of the prepared pencil and every returned eigenvector.
/// Mechanical displacement is mass-normalized; each other owner retains its
/// original units. A unit Euclidean vector is not unit physical energy.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPoleCoordinates {
    pub modal_displacement: Range<usize>,
    pub radiation_displacement: Range<usize>,
    pub duct_energy: Range<usize>,
    pub branch_energy: Vec<(PortId, Range<usize>)>,
    pub terminal_state: Range<usize>,
    pub terminal_flow: Vec<(PortId, usize)>,
    pub junction_flow: Vec<(PortId, usize)>,
    pub channels: Vec<DgPoleChannelCoordinates>,
    /// Energy-normalized internal strains in original Maxwell branch order.
    pub material_memory: Range<usize>,
}

/// One local free mode of the complete assembly, in exp(s*t).
#[derive(Clone, Debug, PartialEq)]
pub struct DgLoadedModalPole {
    pub laplace_pole_per_s: C64,
    /// omega=i*s under the AC owner's exp(-i*omega*t) convention.
    pub angular_frequency_rad_s: C64,
    /// Original mixed-unit coordinates, with the native numerical normalization.
    pub state: Vec<C64>,
    pub modal_velocity_m_sqrt_kg_per_s: Vec<C64>,
    pub radiation_velocity: Vec<C64>,
    pub inlet_flow_m3_s: C64,
    pub inlet_pressure_pa: C64,
    /// Maximum original-row backward residual, scaled by each split coefficient
    /// row envelope times the full state infinity norm. This is coordinate
    /// dependent numerical evidence, not an error bound on physical mode shape.
    pub relative_residual: f64,
    pub iterations: usize,
}

#[derive(Debug)]
pub enum DgLoadedPoleError {
    Source(DgHarmonicCouplingError),
    Solve(NonlinearEigError),
    Residual { relative: f64, tolerance: f64 },
}
impl From<DgHarmonicCouplingError> for DgLoadedPoleError {
    fn from(error: DgHarmonicCouplingError) -> Self {
        Self::Source(error)
    }
}
impl From<ModalCouplingError> for DgLoadedPoleError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Source(error.into())
    }
}
impl From<NonlinearEigError> for DgLoadedPoleError {
    fn from(error: NonlinearEigError) -> Self {
        Self::Solve(error)
    }
}
impl core::fmt::Display for DgLoadedPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete modal/DG pole refused: {self:?}")
    }
}
impl std::error::Error for DgLoadedPoleError {}

/// Immutable complete polynomial pencil. Repeated local refinements reuse all
/// assembled real coefficients and fixed row scales. No source history changes.
pub struct PreparedModalDgPoles<'a> {
    source: &'a CoupledModalSystem,
    load: &'a DgPortLoad,
    pencil: Pencil,
    coordinates: DgPoleCoordinates,
    flow_column: Vec<f64>,
    pressure_column: Vec<f64>,
}

impl operator::Entries for Pencil {
    fn dimension(&self) -> usize {
        self.dimension
    }
    fn set(&mut self, i: usize, j: usize, constant: f64, linear: f64) {
        self.constant[i * self.dimension + j] = constant;
        self.linear[i * self.dimension + j] = linear;
    }
    fn clear_row(&mut self, row: usize) {
        let span = row * self.dimension..(row + 1) * self.dimension;
        self.constant[span.clone()].fill(0.0);
        self.linear[span].fill(0.0);
        self.quadratic_diagonal[row] = 0.0;
    }
}

impl CoupledModalSystem {
    /// Prepare true free resonance/relaxation of the original connected owners.
    ///
    /// Acoustic blocks and constraints share the AC assembler. Material arms
    /// remain explicit, including an unobservable arm at s=-1/tau. The full
    /// nonproportional mechanical damping, prepared anchored-loss columns and
    /// direct acoustic oscillator load are retained. Rest extensions set the
    /// static equilibrium; they do not alter this frozen linear tangent.
    ///
    /// No dense full-state exponential or sample clock is needed. Preparation
    /// reserves 64*d^2+512*d bytes and d^2*(64+8*(links+arms+ports+branches+gaps))
    /// work units before coefficient allocation, with d<=1024. The returned
    /// model has its own explicitly budgeted local queries. Local plane-stress
    /// strain admission and driven mouths require their nonlinear tangent owners
    /// and refuse here. This is not a nonlinear operating-point solver.
    pub fn prepare_resonance_poles_with_load<'a>(
        &'a self,
        load: &'a DgPortLoad,
        flow_column: &[f64],
        budget: DgPolePreparation,
        gate: &CancelGate,
    ) -> Result<PreparedModalDgPoles<'a>, DgLoadedPoleError> {
        poll(Some(gate))?;
        if load.has_driven_mouth() || !self.plane_stress_material_points.is_empty() {
            return Err(invalid("joint poles require the original linear owners without driven-mouth or local plane-stress tangents").into());
        }
        if !budget.reference_rate_per_s.is_finite() || budget.reference_rate_per_s <= 0.0 {
            return Err(invalid("joint pole reference rate must be positive and finite").into());
        }
        let n = self.mode_count();
        let air = self.radiation.as_ref().map_or(0, |r| r.model.poles.len());
        let offset = n
            .checked_add(air)
            .ok_or_else(|| invalid("joint pole mode extent overflow"))?;
        let view = load.duct().continuous_operator();
        let branches = load.branches();
        let branch_views: Vec<_> = branches
            .iter()
            .map(|branch| branch.duct().continuous_operator())
            .collect();
        let terminal = load.terminal().map(|owner| (owner, load.terminal_ports()));
        let binding = operator::DgBinding::new(
            n,
            &view,
            load.inlet(),
            flow_column,
            terminal,
            branches,
            &branch_views,
            load.terminal_channels(),
            gate,
        )?;
        let memory_offset = offset
            .checked_add(binding.additional_dimension()?)
            .ok_or_else(|| invalid("joint pole material offset overflow"))?;
        let dimension = memory_offset
            .checked_add(self.relaxation_branches.len())
            .ok_or_else(|| invalid("joint pole complete dimension overflow"))?;
        let d = dimension as u128;
        let sources = (self.connections.len() as u128)
            + (self.relaxation_branches.len() as u128)
            + (binding.terminal_ports.len() as u128)
            + (branches.len() as u128)
            + (binding.channels.len() as u128);
        let bytes = 64 * d * d + 512 * d;
        let work = d * d * (64 + 8 * sources);
        if dimension == 0
            || dimension > MAX_NONLINEAR_POLE_DIMENSION
            || dimension > budget.maximum_dimension
            || bytes > budget.maximum_payload_bytes as u128
            || work > budget.maximum_work as u128
        {
            return Err(invalid("complete joint pole preparation exceeds its explicit dimension, byte or work budget").into());
        }
        let mut pencil = Pencil::new(dimension, budget.reference_rate_per_s);
        for (i, (mode, &damping)) in self
            .models
            .iter()
            .flat_map(|m| m.modes().iter().zip(m.damping_rates_per_s()))
            .enumerate()
        {
            poll(Some(gate))?;
            let w = mode.angular_frequency_rad_s;
            pencil.constant[i * dimension + i] = w * w;
            pencil.linear[i * dimension + i] = damping;
            pencil.quadratic_diagonal[i] = 1.0;
        }
        for (index, (column, law)) in self.columns.iter().zip(&self.connections).enumerate() {
            let damping = self.harmonic_damping_column(index)?;
            for i in 0..n {
                poll(Some(gate))?;
                for j in 0..n {
                    pencil.constant[i * dimension + j] += law.stiffness_n_m * column[i] * column[j];
                    pencil.linear[i * dimension + j] += law.damping_n_s_m * damping[i] * damping[j];
                }
            }
        }
        for (index, arm) in self.relaxation_branches.iter().enumerate() {
            let memory = memory_offset + index;
            let root = det::sqrt(arm.stiffness);
            let rate = arm.relaxation_time_s.recip();
            pencil.constant[memory * dimension + memory] = rate;
            pencil.linear[memory * dimension + memory] = 1.0;
            for i in 0..n {
                poll(Some(gate))?;
                let b = root * arm.projection[2 * i];
                pencil.constant[i * dimension + memory] = -b;
                pencil.constant[memory * dimension + i] = -rate * b;
                for j in 0..n {
                    pencil.constant[i * dimension + j] += b * (root * arm.projection[2 * j]);
                }
            }
        }
        if let Some(air) = &self.radiation {
            for (j, pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let row = n + j;
                pencil.constant[row * dimension + row] = pole.omega * pole.omega;
                pencil.linear[row * dimension + row] = 2.0 * pole.zeta * pole.omega;
                pencil.quadratic_diagonal[row] = 1.0;
                for (&i, &g) in air.source_modes.iter().zip(&pole.coupling) {
                    pencil.linear[i * dimension + row] = g;
                    pencil.linear[row * dimension + i] = -g;
                }
            }
        }
        let appended = binding.append(
            offset,
            flow_column,
            budget.maximum_payload_bytes,
            &mut pencil,
            gate,
        )?;
        pencil.condition(&mut || gate.is_requested())?;
        let channels = binding
            .channels
            .iter()
            .zip(&appended.channel_offsets)
            .map(|(channel, &at)| {
                let width = channel
                    .channel()
                    .map(|owner| owner.continuous_operator())
                    .transpose()?
                    .map_or(0, |view| view.dimension());
                Ok(DgPoleChannelCoordinates {
                    port: channel.input().port,
                    velocity_cells: at..at + width,
                    pressure_drop: if width == 0 { None } else { Some(at + width) },
                })
            })
            .collect::<Result<Vec<_>, DgHarmonicCouplingError>>()?;
        let coordinates = DgPoleCoordinates {
            modal_displacement: 0..n,
            radiation_displacement: n..offset,
            duct_energy: offset..offset + view.dimension(),
            branch_energy: branches
                .iter()
                .enumerate()
                .map(|(i, branch)| {
                    (
                        branch.port(),
                        offset + binding.branch_offsets[i]
                            ..offset + binding.branch_offsets[i] + branch_views[i].dimension(),
                    )
                })
                .collect(),
            terminal_state: appended.terminal_offset..appended.flow_offset,
            terminal_flow: binding
                .terminal_ports
                .iter()
                .enumerate()
                .map(|(i, &port)| (port, appended.flow_offset + i))
                .collect(),
            junction_flow: branches
                .iter()
                .enumerate()
                .map(|(i, branch)| (branch.port(), appended.junction_offset + i))
                .collect(),
            channels,
            material_memory: memory_offset..dimension,
        };
        poll(Some(gate))?;
        Ok(PreparedModalDgPoles {
            source: self,
            load,
            pencil,
            coordinates,
            flow_column: flow_column.to_vec(),
            pressure_column: binding.pressure_column.to_vec(),
        })
    }
}

impl PreparedModalDgPoles<'_> {
    pub fn dimension(&self) -> usize {
        self.pencil.dimension
    }
    pub fn coordinates(&self) -> &DgPoleCoordinates {
        &self.coordinates
    }

    fn admit_oscillation(&self, s: C64, gate: &CancelGate) -> Result<(), DgLoadedPoleError> {
        poll(Some(gate))?;
        if !complex_finite(s) {
            return Err(invalid("joint pole must be finite").into());
        }
        // A real relaxation coordinate has no oscillation frequency. Its rate
        // is retained without assigning a resonance Q or pretending it is AC.
        if s.im != 0.0 {
            let omega = s.im.abs();
            self.load
                .duct()
                .continuous_operator()
                .admit_frequency(omega)
                .map_err(|e| duct_error(e, gate))?;
            for branch in self.load.branches() {
                poll(Some(gate))?;
                branch
                    .duct()
                    .continuous_operator()
                    .admit_frequency(omega)
                    .map_err(|e| duct_error(e, gate))?;
            }
            if let Some(owner) = self.load.terminal() {
                owner
                    .continuous_operator()
                    .admit_frequency(omega)
                    .map_err(|e| invalid(e.0))?;
            }
        }
        Ok(())
    }

    /// Build a seed from the existing joint AC result, retaining all owner
    /// fields. This is only a Newton starting vector; the original pole rows
    /// independently decide acceptance even when the response came from a
    /// different parameter snapshot with the same coordinate layout.
    pub fn seed_from_harmonic_response(
        &self,
        response: &DgLoadedModalHarmonicResponse,
    ) -> Result<Vec<C64>, DgLoadedPoleError> {
        let layout = &self.coordinates;
        let mut state = vec![C64::ZERO; self.dimension()];
        let copy = |state: &mut [C64], span: Range<usize>, values: &[C64]| {
            if span.len() != values.len() || values.iter().any(|&v| !complex_finite(v)) {
                return Err(invalid(
                    "joint pole seed requires every original finite owner coordinate",
                ));
            }
            state[span].copy_from_slice(values);
            Ok(())
        };
        copy(
            &mut state,
            layout.modal_displacement.clone(),
            &response.modal.displacement_m_sqrt_kg,
        )?;
        copy(
            &mut state,
            layout.material_memory.clone(),
            &response.modal.material_memory,
        )?;
        copy(
            &mut state,
            layout.duct_energy.clone(),
            &response.duct.energy_coordinates,
        )?;
        if response.modal.radiation_state.len() != 2 * layout.radiation_displacement.len()
            || response.branches.len() != layout.branch_energy.len()
            || response.channels.len() != layout.channels.len()
        {
            return Err(invalid("joint pole seed has a different owner layout").into());
        }
        for (i, pair) in response.modal.radiation_state.chunks_exact(2).enumerate() {
            state[layout.radiation_displacement.start + i] = pair[0];
        }
        for ((port, span), branch) in layout.branch_energy.iter().zip(&response.branches) {
            if *port != branch.port {
                return Err(invalid("joint pole branch seed order differs").into());
            }
            copy(&mut state, span.clone(), &branch.duct.energy_coordinates)?;
        }
        for ((port, index), branch) in layout.junction_flow.iter().zip(&response.branches) {
            if *port != branch.port {
                return Err(invalid("joint pole junction seed order differs").into());
            }
            state[*index] = branch.junction_flow_m3_s;
        }
        match (&response.terminal, self.load.terminal()) {
            (Some(terminal), Some(_)) => {
                if terminal.ports.len() != layout.terminal_flow.len()
                    || terminal.outward_flow_m3_s.len() != layout.terminal_flow.len()
                {
                    return Err(invalid("joint pole terminal seed is incomplete").into());
                }
                copy(&mut state, layout.terminal_state.clone(), &terminal.state)?;
                for (i, &(port, index)) in layout.terminal_flow.iter().enumerate() {
                    if port != terminal.ports[i] {
                        return Err(invalid("joint pole terminal order differs").into());
                    }
                    state[index] = terminal.outward_flow_m3_s[i];
                }
            }
            (None, None) => (),
            _ => return Err(invalid("joint pole seed has a different terminal owner").into()),
        }
        for (coordinate, channel) in layout.channels.iter().zip(&response.channels) {
            if coordinate.port != channel.port {
                return Err(invalid("joint pole channel seed order differs").into());
            }
            match (&channel.channel, coordinate.pressure_drop) {
                (Some(field), Some(index)) => {
                    copy(
                        &mut state,
                        coordinate.velocity_cells.clone(),
                        &field.velocities_m_s,
                    )?;
                    state[index] = channel.pressure_drop_pa;
                }
                (None, None) => (),
                _ => {
                    return Err(
                        invalid("joint pole seed changes a finite gap or rigid closure").into(),
                    );
                }
            }
        }
        if state.iter().any(|&v| !complex_finite(v)) {
            return Err(invalid("joint pole seed contains a nonfinite field or multiplier").into());
        }
        Ok(state)
    }

    /// Refine one locally selected, isolated finite pole, retaining growth and
    /// pure real relaxation. This is not a complete spectrum search, a forward
    /// error enclosure, a transfer-observability claim or evidence of attraction.
    ///
    /// Oscillatory seeds/results obey each original duct and terminal band.
    /// A real-axis source band does not certify empirical continuation to complex
    /// frequency. Normalized eigenvectors have arbitrary amplitude, so neither
    /// finite-amplitude material/Mach admission nor terminal heat is invented.
    /// The already prepared cache is excluded from this query's resource budget;
    /// callback work and returned physical observations are included.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalPole, DgLoadedPoleError> {
        self.admit_oscillation(laplace_seed_per_s, gate)?;
        let n = self.dimension() as u128;
        let bytes = 128 * n;
        let work = (128 + 16 * options.max_iterations as u128) * n * n;
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(
                invalid("joint pole callbacks and observations exceed query resources").into(),
            );
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let pole = newton_pole_with_cancel(
            self.dimension(),
            laplace_seed_per_s,
            state_seed,
            &numerical,
            |s, cancel| self.pencil.evaluate(s, cancel),
            &mut || gate.is_requested(),
        )?;
        self.observe_pole(pole, options.residual_tolerance, gate)
    }

    // Local and contour solves use exactly the same source admission,
    // original-equation residual and physical observations.
    fn observe_pole(
        &self,
        pole: NonlinearPole,
        residual_tolerance: f64,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalPole, DgLoadedPoleError> {
        self.admit_oscillation(pole.value, gate)?;
        let relative = self
            .pencil
            .relative_residual(pole.value, &pole.vector, &mut || gate.is_requested())?;
        let tolerance = residual_tolerance.min(self.source.config.solve_relative_tolerance);
        if relative > tolerance {
            return Err(DgLoadedPoleError::Residual {
                relative,
                tolerance,
            });
        }
        let layout = &self.coordinates;
        let s = pole.value;
        let modal_velocity: Vec<_> = pole.vector[layout.modal_displacement.clone()]
            .iter()
            .map(|&q| s * q)
            .collect();
        let radiation_velocity: Vec<_> = pole.vector[layout.radiation_displacement.clone()]
            .iter()
            .map(|&q| s * q)
            .collect();
        let inlet_flow = projection(&self.flow_column, &modal_velocity);
        let inlet_pressure = projection(
            &self.pressure_column,
            &pole.vector[layout.duct_energy.clone()],
        );
        if modal_velocity
            .iter()
            .chain(&radiation_velocity)
            .any(|&v| !complex_finite(v))
            || !complex_finite(inlet_flow)
            || !complex_finite(inlet_pressure)
        {
            return Err(invalid("joint pole physical observation is not representable").into());
        }
        poll(Some(gate))?;
        Ok(DgLoadedModalPole {
            laplace_pole_per_s: s,
            angular_frequency_rad_s: C64::new(-s.im, s.re),
            state: pole.vector,
            modal_velocity_m_sqrt_kg_per_s: modal_velocity,
            radiation_velocity,
            inlet_flow_m3_s: inlet_flow,
            inlet_pressure_pa: inlet_pressure,
            relative_residual: relative,
            iterations: pole.iterations,
        })
    }
}
