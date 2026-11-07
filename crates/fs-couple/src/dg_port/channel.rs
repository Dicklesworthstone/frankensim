//! Fixed finite gap channels and rigid closures in the shared terminal group.
//!
//! The original transverse no-slip momentum owner supplies the affine drop
//! `p_duct - p_exterior = Z Q - H`. Its actual field, material viscosity,
//! kinetic storage and separately measured viscous/numerical losses survive
//! the simultaneous DG/rational solve. A rigid closure constrains Q exactly;
//! it does not delete a row from the exterior's retained physical history.

use super::*;
use crate::bernoulli_aperture::slit_channel::{
    ChannelGasUpdate, FixedChannel, SlitChannel, SlitChannelSpec, SlitChannelStep, SlitExit,
    annular::{AnnularChannel, AnnularChannelSpec},
    circular::{CircularChannel, CircularChannelSpec},
};
use fs_material::gas::GasState;

/// Sum of independently resolved transverse cells admitted to one DG source.
pub const MAX_DG_TERMINAL_TRANSVERSE_CELLS: usize = 4096;

/// A fixed port reduction. No pad motion, swept volume or contact is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DgTerminalChannelKind {
    /// Fixed rigid closure. Interior and exterior pressures may differ.
    Closed,
    /// Fixed, parallel, no-slip walls with their original transverse fluid
    /// momentum field. End pressures are static; entrance and exit jets,
    /// lateral walls and compressible channel storage are not resolved.
    Slit {
        width_m: f64,
        height_m: f64,
        length_m: f64,
        transverse_cells: usize,
        /// Maximum absolute resolved cell speed / local sound speed, in (0,.3].
        maximum_mach: f64,
    },
    /// Fixed circular no-slip passage with an independently resolved radial
    /// momentum field. This is additional geometry between the duct pressure
    /// plane and its exterior load, not a second loss on an existing chimney.
    /// Axial development, compressible storage and entrance jets are absent.
    Circular {
        radius_m: f64,
        length_m: f64,
        radial_cells: usize,
        /// Maximum absolute resolved cell speed / local sound speed, in (0,.3].
        maximum_mach: f64,
    },
    /// Fixed concentric no-slip passage with both material walls and its
    /// original annular momentum field. Flow is axial through the supplied
    /// length; eccentricity, radial pad openings and entrance jets are absent.
    Annular {
        inner_radius_m: f64,
        outer_radius_m: f64,
        length_m: f64,
        radial_cells: usize,
        /// Maximum absolute resolved cell speed / local sound speed, in (0,.3].
        maximum_mach: f64,
    },
    /// The quasi-steady sharp-edged orifice law of an open tone hole, in series with
    /// the terminal's exterior load: `p_duct - p_exterior = K |Q| Q` with
    /// `K = rho / (2 C_vc^2 A^2)`, `A` the actual chimney outlet area and `rho` the
    /// terminal's own gas (the live reed card's `open_hole_orifice` law, which takes
    /// the same form). The coefficient uses the PREVIOUS accepted sample's `|Q|`, so
    /// each sample solves a linear terminal group whose series resistance
    /// `K |Q_prev|` is nonnegative: the group stays passive and an absent channel
    /// changes no bit. The lag is one mechanical step and is measured per render
    /// (`DgOrificeReport`). No storage; the law is inertialess.
    Orifice {
        /// Vena-contracta coefficient `C_vc` of the card, finite and in (0, 2].
        vena_contracta: f64,
    },
}

impl DgTerminalChannelKind {
    pub(super) fn section(self) -> Option<(f64, usize)> {
        match self {
            Self::Closed | Self::Orifice { .. } => None,
            Self::Slit { width_m, height_m, transverse_cells, .. } =>
                Some((width_m * height_m, transverse_cells)),
            Self::Circular { radius_m, radial_cells, .. } =>
                Some((core::f64::consts::PI * radius_m * radius_m, radial_cells)),
            Self::Annular { inner_radius_m, outer_radius_m, radial_cells, .. } =>
                Some((core::f64::consts::PI * (outer_radius_m - inner_radius_m)
                    * (outer_radius_m + inner_radius_m), radial_cells)),
        }
    }

    fn maximum_mach(self) -> Option<f64> {
        match self {
            Self::Closed | Self::Orifice { .. } => None,
            Self::Slit { maximum_mach, .. } | Self::Circular { maximum_mach, .. }
            | Self::Annular { maximum_mach, .. } => Some(maximum_mach),
        }
    }
}

/// Original external port and an explicit finite gap or fixed closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgTerminalChannelInput {
    pub port: PortId,
    pub kind: DgTerminalChannelKind,
}

impl DgTerminalChannelInput {
    /// Geometry-only admission, usable before source preparation/material I/O.
    /// The planar reduction requires height <= width/10 and length/10;
    /// the circular reduction requires radius <= length/10, and the annular
    /// reduction requires a positive radial gap <= length/10.
    pub fn validate(&self) -> Result<(), AcousticRealizeError> {
        if let DgTerminalChannelKind::Slit {
            width_m, height_m, length_m, transverse_cells, maximum_mach,
        } = self.kind {
            if ![width_m, height_m, length_m, maximum_mach].iter()
                .all(|value| value.is_finite() && *value > 0.0)
                || height_m > 0.1 * width_m || height_m > 0.1 * length_m
                || !(2..=MAX_DG_TERMINAL_TRANSVERSE_CELLS).contains(&transverse_cells)
                || maximum_mach > 0.3
                || !(width_m * height_m).is_finite()
            {
                return Err(invalid("DG terminal slit requires a finite thin planar gap, 2..=4096 cells and Mach in (0,.3]"));
            }
        }
        if let DgTerminalChannelKind::Circular {
            radius_m, length_m, radial_cells, maximum_mach,
        } = self.kind {
            if ![radius_m, length_m, maximum_mach].iter()
                .all(|value| value.is_finite() && *value > 0.0)
                || radius_m > 0.1 * length_m
                || !(2..=MAX_DG_TERMINAL_TRANSVERSE_CELLS).contains(&radial_cells)
                || maximum_mach > 0.3
                || !(core::f64::consts::PI * radius_m * radius_m).is_finite()
            {
                return Err(invalid("DG terminal circular passage requires a finite slender cylinder, 2..=4096 cells and Mach in (0,.3]"));
            }
        }
        if let DgTerminalChannelKind::Annular {
            inner_radius_m, outer_radius_m, length_m, radial_cells, maximum_mach,
        } = self.kind {
            let gap = outer_radius_m - inner_radius_m;
            let area = core::f64::consts::PI * gap * (outer_radius_m + inner_radius_m);
            if ![inner_radius_m, outer_radius_m, length_m, maximum_mach, gap, area].iter()
                .all(|value| value.is_finite() && *value > 0.0)
                || gap > 0.1 * length_m
                || !(2..=MAX_DG_TERMINAL_TRANSVERSE_CELLS).contains(&radial_cells)
                || maximum_mach > 0.3
            {
                return Err(invalid("DG terminal annulus requires finite positive ordered radii, area and a slender gap, 2..=4096 cells and Mach in (0,.3]"));
            }
        }
        if let DgTerminalChannelKind::Orifice { vena_contracta } = self.kind {
            if !(vena_contracta.is_finite() && vena_contracta > 0.0 && vena_contracta <= 2.0) {
                return Err(invalid("DG terminal orifice requires a finite vena-contracta coefficient in (0, 2]"));
            }
        }
        Ok(())
    }
}

/// The lagged sharp-edged orifice of one terminal: its coefficient, area and the
/// running account of how far the solved drop (`K |Q_prev| Q`) has stayed from the
/// law it stands for (`K |Q| Q`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct OrificeState {
    vena_contracta: f64,
    area_m2: f64,
    last_flow_m3_s: f64,
    last_step_work_j: f64,
    max_lag_drop_pa: f64,
    lag_work_j: f64,
    law_work_j: f64,
    steps: u64,
}

impl OrificeState {
    /// `K = rho / (2 C_vc^2 A^2)` in the terminal's own gas [Pa s^2/m^6].
    fn coefficient(&self, density_kg_m3: f64) -> f64 {
        density_kg_m3 / (2.0 * self.vena_contracta * self.vena_contracta * self.area_m2 * self.area_m2)
    }

    /// The affine load of the next sample: `dp = K |Q_prev| Q`, no offset. Zero at rest (a
    /// resistance of exactly zero adds nothing to the terminal matrix).
    fn load(&self, density_kg_m3: f64) -> crate::bernoulli_aperture::slit_channel::SlitChannelLoad {
        crate::bernoulli_aperture::slit_channel::SlitChannelLoad {
            impedance_pa_s_m3: self.coefficient(density_kg_m3) * self.last_flow_m3_s.abs(),
            history_pressure_pa: 0.0,
        }
    }
}

/// Measured account of one orifice terminal after the last accepted sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgOrificeReport {
    pub vena_contracta: f64,
    /// Actual chimney outlet area the law is written on [m^2].
    pub area_m2: f64,
    /// Outward flow of the last accepted sample [m^3/s].
    pub last_flow_m3_s: f64,
    /// Accepted samples.
    pub steps: u64,
    /// Largest `|K (|Q| - |Q_prev|) Q|` over the accepted samples [Pa]: how far the solved drop was from the law.
    pub max_lag_drop_pa: f64,
    /// `sum |K (|Q| - |Q_prev|) Q| |Q| dt` [J] and `sum K |Q|^3 dt` [J]: their ratio is the share of the law's
    /// dissipation the lag misattributed over the render.
    pub lag_work_j: f64,
    pub law_work_j: f64,
    /// Work the solved resistance took, `sum K |Q_prev| Q^2 dt` [J]; all of it is dissipated.
    pub dissipated_j: f64,
}

/// Bound material and the unchanged original channel owner, without a mutable
/// escape. Closed ports retain their source gas but no fabricated fluid state.
#[derive(Clone, Debug)]
pub struct DgTerminalChannelLoad {
    input: DgTerminalChannelInput,
    gas: GasState,
    channel: Option<FixedChannel>,
    orifice: Option<OrificeState>,
    cumulative_pressure_work_j: f64,
    cumulative_viscous_dissipation_j: f64,
    cumulative_numerical_dissipation_j: f64,
    cumulative_parameter_work_j: f64,
}

/// One terminal's original physical gas transfer. Rigid closures have no
/// fluid storage and return no fabricated channel transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgTerminalChannelGasUpdate {
    pub port: PortId,
    pub transfer: Option<ChannelGasUpdate>,
}

impl DgTerminalChannelLoad {
    pub const fn input(&self) -> &DgTerminalChannelInput { &self.input }
    pub const fn gas(&self) -> &GasState { &self.gas }
    pub const fn channel(&self) -> Option<&FixedChannel> { self.channel.as_ref() }
    pub fn stored_energy_j(&self) -> f64 {
        self.channel.as_ref().map_or(0.0, FixedChannel::stored_energy_j)
    }
    /// Signed gas-intervention work, separate from pressure work and heat [J].
    pub const fn cumulative_parameter_work_j(&self) -> f64 {
        self.cumulative_parameter_work_j
    }
    pub const fn cumulative_pressure_work_j(&self) -> f64 { self.cumulative_pressure_work_j }
    pub const fn cumulative_viscous_dissipation_j(&self) -> f64 {
        self.cumulative_viscous_dissipation_j
    }
    pub const fn cumulative_numerical_dissipation_j(&self) -> f64 {
        self.cumulative_numerical_dissipation_j
    }
    pub(super) fn cumulative_dissipated_j(&self) -> f64 {
        if self.orifice.is_some() {
            // The orifice stores nothing: every joule its resistance takes is dissipated.
            self.cumulative_pressure_work_j
        } else {
            self.cumulative_viscous_dissipation_j + self.cumulative_numerical_dissipation_j
        }
    }
    /// The orifice account of this terminal after the last accepted sample, if it is one.
    pub fn orifice_report(&self) -> Option<DgOrificeReport> {
        self.orifice.map(|orifice| DgOrificeReport {
            vena_contracta: orifice.vena_contracta,
            area_m2: orifice.area_m2,
            last_flow_m3_s: orifice.last_flow_m3_s,
            steps: orifice.steps,
            max_lag_drop_pa: orifice.max_lag_drop_pa,
            lag_work_j: orifice.lag_work_j,
            law_work_j: orifice.law_work_j,
            dissipated_j: self.cumulative_pressure_work_j,
        })
    }
    pub(super) fn pressure_load(&self)
        -> Result<Option<crate::bernoulli_aperture::slit_channel::SlitChannelLoad>, AcousticRealizeError>
    {
        if let Some(orifice) = &self.orifice {
            return Ok(Some(orifice.load(self.gas.density)));
        }
        self.channel.as_ref().map(FixedChannel::pressure_load).transpose()
    }

    pub(super) fn prepare_gas(
        &self, gas: &GasState, maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, DgTerminalChannelGasUpdate), AcousticRealizeError> {
        if cancelled() { return Err(AcousticRealizeError::Cancelled); }
        gas.validate_coherence().map_err(|_| invalid("live terminal channel requires a coherent complete gas"))?;
        let mut candidate = self.clone();
        let transfer = if let Some(channel) = &mut candidate.channel {
            let transfer = channel.rebind_gas_controlled(
                gas, maximum_absolute_parameter_work_j, cancelled,
            )?;
            let maximum_mach = self.input.kind.maximum_mach()
                .ok_or_else(|| invalid("live terminal channel lost its physical kind"))?;
            if channel.velocities_m_s().iter().any(|speed|
                speed.abs() > maximum_mach * gas.sound_speed)
            {
                return Err(invalid("live terminal gas transfer exceeds the retained resolved-speed Mach limit"));
            }
            candidate.cumulative_parameter_work_j += transfer.parameter_work_j;
            if !candidate.cumulative_parameter_work_j.is_finite() {
                return Err(invalid("live terminal gas work cannot be represented"));
            }
            Some(transfer)
        } else { None };
        candidate.gas = *gas;
        if cancelled() { return Err(AcousticRealizeError::Cancelled); }
        Ok((candidate, DgTerminalChannelGasUpdate { port: self.input.port, transfer }))
    }
}

/// Original channel exchange; rigid closure has no fluid step and zero work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgTerminalChannelFrame {
    pub port: PortId,
    pub outward_flow_m3_s: f64,
    /// Independently observed duct-side minus exterior-side mean pressure [Pa].
    pub mean_pressure_drop_pa: f64,
    /// Original channel's next-step constitutive drop [Pa]. For a rigid
    /// closure this is its unconstrained reaction pressure instead.
    pub constitutive_pressure_drop_pa: f64,
    /// Work entering the original channel [J], including kinetic storage.
    pub pressure_work_j: f64,
    /// Observed pressure drop times held volume minus original channel work [J].
    pub pressure_work_residual_j: f64,
    /// Original accepted FV flow minus the requested shared held flow [m³/s].
    pub flow_residual_m3_s: f64,
    /// Complete original storage, wall heat and numerical-loss account.
    pub step: Option<SlitChannelStep>,
}

impl DgTerminalChannelFrame {
    // Retain the independent original storage/loss terms when a requested
    // zero flow cancels inertial history to roundoff. Scaling by Q alone
    // would demand a bit-exact FV cancellation instead of the fixed work gate.
    pub(super) fn work_scale_j(&self) -> f64 {
        // A rigid closure has no work (0.0); an orifice has no step but does work.
        self.step.map_or(self.pressure_work_j.abs(), |step| step.pressure_work_j.abs()
            + step.stored_delta_j.abs() + step.viscous_dissipation_j
            + step.numerical_dissipation_j)
    }
}

impl DgPortLoad {
    /// Insert fixed finite slit/circular/annular channels or rigid closures before the existing
    /// full exterior impedance. Unspecified terminals retain their original
    /// direct connection. Material density, viscosity and sound speed come from
    /// the ACTUAL local terminal, including the end of a finite DG chimney.
    ///
    /// All supplied geometry is cold; runtime opening changes would require
    /// moving-boundary work, displaced volume and conservative field transfer
    /// and are deliberately not exposed by this fixed-geometry constructor.
    /// All geometries use the original laminar unsteady Stokes owner. This
    /// does not change the exterior surface geometry or resolve a radial pad
    /// opening. Each channel is an additional supplied passage; it must not
    /// duplicate geometry or losses already resolved by a DG chimney.
    pub fn with_terminal_channels(
        mut self, inputs: Vec<DgTerminalChannelInput>, maximum_transverse_cells: usize,
    ) -> Result<Self, AcousticRealizeError> {
        if inputs.is_empty() { return Ok(self); }
        if self.state.step_count() != 0
            || self.branches.iter().any(|branch| branch.state.step_count() != 0)
            || !self.terminal_channels.is_empty()
            || self.terminal.as_ref().is_none_or(|load|
                load.pole_states().iter().any(|state| state.re != 0.0 || state.im != 0.0))
            || inputs.len() > self.terminal_ports().len()
            || maximum_transverse_cells > MAX_DG_TERMINAL_TRANSVERSE_CELLS
        {
            return Err(invalid("DG terminal channels require one cold complete terminal group and a bounded cell allowance"));
        }
        let mut cells = 0usize;
        for (index, input) in inputs.iter().enumerate() {
            input.validate()?;
            if !self.terminal_ports().contains(&input.port)
                || inputs[..index].iter().any(|previous| previous.port == input.port)
            {
                return Err(invalid("DG terminal channels require distinct actual terminal ports"));
            }
            if matches!(input.kind, DgTerminalChannelKind::Orifice { .. })
                && (!self.branches.iter().any(|branch| branch.port() == input.port)
                    || !self.terminal_area_m2(input.port)?.is_some_and(|area| area.is_finite() && area > 0.0)
                    || self.has_driven_mouth())
            {
                return Err(invalid("DG terminal orifice is a tone-hole law: it needs a finite DG chimney with an outlet area and no driven exterior mouth"));
            }
            if let Some((channel_area, transverse_cells)) = input.kind.section() {
                cells = cells.checked_add(transverse_cells)
                    .ok_or_else(|| invalid("DG terminal transverse state count overflow"))?;
                // A flush side port has no aperture area of its own; its passage
                // still cannot exceed the bore section it is cut into.
                let attached_area = match self.terminal_area_m2(input.port)? {
                    Some(area) => area,
                    None => self.duct.side_port_source(input.port)
                        .map_err(|_| invalid("DG terminal lost its original side-port geometry"))?
                        .bore_area_m2,
                };
                if channel_area > attached_area * (1.0 + 64.0 * f64::EPSILON) {
                    return Err(invalid("DG terminal channel area exceeds its actual attached duct area"));
                }
            }
        }
        if cells > maximum_transverse_cells {
            return Err(invalid("DG terminal transverse fields exceed their declared combined cell allowance"));
        }
        let mut channels = Vec::with_capacity(inputs.len());
        for input in inputs {
            let gas = self.terminal_gas(input.port)?;
            let orifice = if let DgTerminalChannelKind::Orifice { vena_contracta } = input.kind {
                let area_m2 = self.terminal_area_m2(input.port)?
                    .ok_or_else(|| invalid("DG terminal orifice lost its chimney outlet area"))?;
                Some(OrificeState {
                    vena_contracta, area_m2, last_flow_m3_s: 0.0, last_step_work_j: 0.0,
                    max_lag_drop_pa: 0.0, lag_work_j: 0.0, law_work_j: 0.0, steps: 0,
                })
            } else { None };
            let channel = match input.kind {
                DgTerminalChannelKind::Closed | DgTerminalChannelKind::Orifice { .. } => None,
                DgTerminalChannelKind::Slit { width_m, height_m, length_m, transverse_cells, .. } => {
                    Some(SlitChannel::from_gas(&gas, SlitChannelSpec {
                        length_m, width_m, height_m, transverse_cells,
                        sample_period_s: self.propagator.dt(), exit: SlitExit::StaticPressure,
                    })?.into())
                }
                DgTerminalChannelKind::Circular { radius_m, length_m, radial_cells, .. } => {
                    Some(CircularChannel::from_gas(&gas, CircularChannelSpec {
                        length_m, radius_m, radial_cells,
                        sample_period_s: self.propagator.dt(), exit: SlitExit::StaticPressure,
                    })?.into())
                }
                DgTerminalChannelKind::Annular {
                    inner_radius_m, outer_radius_m, length_m, radial_cells, ..
                } => {
                    Some(AnnularChannel::from_gas(&gas, AnnularChannelSpec {
                        length_m, inner_radius_m, outer_radius_m, radial_cells,
                        sample_period_s: self.propagator.dt(), exit: SlitExit::StaticPressure,
                    })?.into())
                }
            };
            channels.push(DgTerminalChannelLoad { input, gas, channel, orifice,
                cumulative_pressure_work_j: 0.0, cumulative_viscous_dissipation_j: 0.0,
                cumulative_numerical_dissipation_j: 0.0, cumulative_parameter_work_j: 0.0 });
        }
        let terminal = self.terminal.as_ref().expect("admitted complete terminal group");
        let response = self.multiport_response()?;
        let coupling = TerminalCoupling::with_channels(&response, self.inlet,
            self.terminal_ports().to_vec(), terminal, self.has_driven_mouth(), &channels)?;
        let impedance = coupling.relation(response, terminal, &channels)?.inlet.impedance_pa_s_m3;
        self.terminal_channels = channels;
        self.terminal_coupling = Some(coupling);
        self.impedance_pa_s_m3 = impedance;
        Ok(self)
    }

    pub fn terminal_channels(&self) -> &[DgTerminalChannelLoad] { &self.terminal_channels }
    pub fn last_terminal_channel_frames(&self) -> &[DgTerminalChannelFrame] {
        &self.last_terminal_channel_frames
    }

    pub(super) fn step_terminal_channels(
        &self, output: &DgStepOutput, outward_flows: &[f64],
        load_pressure: &[f64], offset: usize,
    ) -> Result<(Vec<DgTerminalChannelLoad>, Vec<DgTerminalChannelFrame>), AcousticRealizeError> {
        let mut candidates = Vec::with_capacity(self.terminal_channels.len());
        let mut frames = Vec::with_capacity(self.terminal_channels.len());
        for owner in &self.terminal_channels {
            let index = self.terminal_ports().iter().position(|port| *port == owner.input.port)
                .ok_or_else(|| invalid("DG terminal channel lost its original load row"))?;
            let port = output.port_outputs.iter().find(|port| port.id == owner.input.port)
                .ok_or_else(|| invalid("DG terminal channel lost its original pressure observation"))?;
            let flow = outward_flows[index];
            let drop = port.conjugate_rate - load_pressure[index + offset];
            let mut candidate = owner.clone();
            let (constitutive_drop, work, flow_residual, step) = if let Some(channel) = &owner.channel {
                let relation = channel.pressure_load()?;
                let constitutive_drop = relation.pressure_drop_pa(flow)?;
                let next = channel.preview_at_flow(flow)?;
                let step = next.last_step();
                let maximum_mach = owner.input.kind.maximum_mach()
                    .ok_or_else(|| invalid("DG terminal channel lost its physical kind"))?;
                if next.velocities_m_s().iter().any(|speed|
                    speed.abs() > maximum_mach * owner.gas.sound_speed)
                {
                    return Err(invalid("DG terminal channel exceeded its declared local resolved-speed Mach limit"));
                }
                let flow_residual = step.flow_m3_s - flow;
                let scale = step.flow_m3_s.abs() + flow.abs()
                    + relation.history_pressure_pa.abs() / relation.impedance_pa_s_m3;
                if !scale.is_finite() || flow_residual.abs() > 1e-10 * scale.max(1e-300) {
                    return Err(invalid("DG terminal channel disagrees with the shared held flow"));
                }
                candidate.cumulative_pressure_work_j += step.pressure_work_j;
                candidate.cumulative_viscous_dissipation_j += step.viscous_dissipation_j;
                candidate.cumulative_numerical_dissipation_j += step.numerical_dissipation_j;
                candidate.channel = Some(next);
                (constitutive_drop, step.pressure_work_j, flow_residual, Some(step))
            } else if let Some(orifice) = &owner.orifice {
                // The sample was solved with the resistance Z = K |Q_prev|: that drop, not the law at the new
                // flow, is what the same-step pressure and work audits must see. The law's own value and the
                // difference are kept for the render's report.
                let k = orifice.coefficient(owner.gas.density);
                let used = k * orifice.last_flow_m3_s.abs();
                let constitutive_drop = used * flow;
                let work = constitutive_drop * (flow * self.propagator.dt());
                let law_drop = k * flow.abs() * flow;
                let lag_drop = (law_drop - constitutive_drop).abs();
                let mut next = *orifice;
                next.last_flow_m3_s = flow;
                next.last_step_work_j = work;
                next.max_lag_drop_pa = next.max_lag_drop_pa.max(lag_drop);
                next.lag_work_j += lag_drop * flow.abs() * self.propagator.dt();
                next.law_work_j += k * flow.abs() * flow * flow * self.propagator.dt();
                next.steps += 1;
                candidate.orifice = Some(next);
                candidate.cumulative_pressure_work_j += work;
                (constitutive_drop, work, 0.0, None)
            } else {
                if flow != 0.0 { return Err(invalid("rigid DG terminal closure admitted nonzero flow")); }
                (drop, 0.0, 0.0, None)
            };
            let residual = drop * (flow * self.propagator.dt()) - work;
            if ![drop, constitutive_drop, work, flow_residual, residual,
                candidate.cumulative_pressure_work_j, candidate.cumulative_dissipated_j()]
                .iter().all(|value| value.is_finite())
            {
                return Err(invalid("DG terminal channel candidate left the finite work domain"));
            }
            frames.push(DgTerminalChannelFrame { port: owner.input.port,
                outward_flow_m3_s: flow, mean_pressure_drop_pa: drop,
                constitutive_pressure_drop_pa: constitutive_drop, pressure_work_j: work,
                pressure_work_residual_j: residual, flow_residual_m3_s: flow_residual, step });
            candidates.push(candidate);
        }
        Ok((candidates, frames))
    }
}

/// Extend the original physical ledger once with each original channel.
/// Its detailed frame keeps numerical dissipation separate from wall heat.
pub(super) fn include_channel_ledgers(
    output: &mut DgStepOutput, channels: &[DgTerminalChannelLoad],
) -> Result<(), AcousticRealizeError> {
    for channel in channels {
        if let Some(orifice) = &channel.orifice {
            // An orifice stores nothing and returns nothing: the work its resistance took is dissipated.
            output.ledger.port_work_step_j += orifice.last_step_work_j;
            output.ledger.cumulative_port_work_j += channel.cumulative_pressure_work_j;
            output.ledger.cumulative_dissipated_j += channel.cumulative_dissipated_j();
            output.ledger.energy_balance_residual_j += channel.cumulative_dissipated_j()
                - channel.cumulative_pressure_work_j;
            continue;
        }
        let Some(owner) = &channel.channel else { continue; };
        let step = owner.last_step();
        output.ledger.stored_energy_j += step.kinetic_energy_j;
        output.ledger.port_work_step_j += step.pressure_work_j;
        output.ledger.cumulative_port_work_j += channel.cumulative_pressure_work_j;
        output.ledger.cumulative_dissipated_j += channel.cumulative_dissipated_j();
        output.ledger.cumulative_parameter_work_j += channel.cumulative_parameter_work_j;
        // Original DG cumulative sign is E - W + D, opposite the channel's
        // per-step pressure-work residual. Keep each owner's convention.
        output.ledger.energy_balance_residual_j += step.kinetic_energy_j
            - channel.cumulative_pressure_work_j + channel.cumulative_dissipated_j()
            - channel.cumulative_parameter_work_j;
    }
    if !channels.is_empty() {
        output.ledger.relative_energy_drift = output.ledger.energy_balance_residual_j.abs()
            / (output.ledger.cumulative_port_work_j.abs()
                + output.ledger.cumulative_parameter_work_j.abs()).max(1e-300);
    }
    let ledger = &output.ledger;
    if ![ledger.stored_energy_j, ledger.port_work_step_j, ledger.cumulative_port_work_j,
        ledger.cumulative_dissipated_j, ledger.energy_balance_residual_j,
        ledger.relative_energy_drift, ledger.cumulative_parameter_work_j]
        .iter().all(|value| value.is_finite())
    {
        return Err(invalid("DG aggregate channel storage/work ledger is not finite"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod orifice_tests;
