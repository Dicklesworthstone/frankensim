//! Retained source walls and native material work on a continuing port load.

use super::*;
use crate::bernoulli_aperture::wall::material::{
    LocalMaterialWall, MaterialWallDg, RegionalMaterialWallDg,
};
use fs_duct::dg::DgWallMaterialUpdate;
use fs_matdb::MaterialCard;
use fs_phs::linear_exact::DenseExponentialBudget;

/// Original main-bore wall sources beside the load's sole native DG operator.
/// Regional entries follow original elements; homogeneous entries retain their
/// declared regional identity rather than manufacturing an averaged material.
#[derive(Clone, Debug, PartialEq)]
pub enum DgPortWallSource {
    /// One complete original law for the whole main bore.
    Uniform(LocalMaterialWall),
    /// One complete original law per native main-bore element.
    Regional(Vec<LocalMaterialWall>),
}

impl DgPortWallSource {
    /// Original whole-bore source, when the binding was explicitly uniform.
    #[must_use]
    pub const fn uniform_wall(&self) -> Option<&LocalMaterialWall> {
        match self {
            Self::Uniform(wall) => Some(wall),
            Self::Regional(_) => None,
        }
    }

    /// Original per-element sources in physical inlet-to-outlet order.
    #[must_use]
    pub fn regional_walls(&self) -> Option<&[LocalMaterialWall]> {
        match self {
            Self::Uniform(_) => None,
            Self::Regional(walls) => Some(walls),
        }
    }

    fn element_walls(&self, elements: usize) -> Vec<LocalMaterialWall> {
        match self {
            Self::Uniform(wall) => vec![wall.clone(); elements],
            Self::Regional(walls) => walls.clone(),
        }
    }

    fn at_temperature_pinned(
        &self,
        card: &MaterialCard,
        temperature_axis: &str,
        initial_temperature_k: f64,
        temperature_k: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AcousticRealizeError> {
        let mut refresh = |wall: &LocalMaterialWall| {
            wall.at_temperature_pinned(
                card,
                temperature_axis,
                initial_temperature_k,
                temperature_k,
                cancelled,
            )
        };
        match self {
            Self::Uniform(wall) => Ok(Self::Uniform(refresh(wall)?)),
            Self::Regional(walls) => walls
                .iter()
                .map(&mut refresh)
                .collect::<Result<Vec<_>, _>>()
                .map(Self::Regional),
        }
    }
}

/// Accepted main-wall temperature intervention in the complete port load.
/// Exterior rational storage is not inferred; the reported represented storage
/// includes the main bore, original finite branches and finite terminal channels.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPortWallUpdate {
    /// Complete refreshed original sources, including pinned material receipts.
    pub source: DgPortWallSource,
    /// Native wall-state transfer and its signed parameter-work account.
    pub native: DgWallMaterialUpdate,
    /// Complete represented acoustic/channel storage before the change [J].
    pub energy_before_j: f64,
    /// Complete represented acoustic/channel storage after the change [J].
    pub energy_after_j: f64,
    /// Signed native work INTO the main acoustic/wall state [J]. The finite
    /// thermal or prescribed background owner must supply this work once.
    pub parameter_work_j: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticRealizeError> {
    if cancelled() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn duct_error(error: fs_duct::DuctError) -> AcousticRealizeError {
    match error {
        fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
        error => AcousticRealizeError::Duct(error),
    }
}

impl DgPortLoad {
    /// Prepare a quiescent load while retaining the complete uniform source law.
    /// The original operator moves into this load exactly once, with all native
    /// gas/wall memory coordinates intact. No parallel acoustic owner is kept.
    pub fn new_with_material_wall(
        wall: MaterialWallDg,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
    ) -> Result<Self, AcousticRealizeError> {
        let source = DgPortWallSource::Uniform(wall.wall().clone());
        let mut load = Self::new(
            wall.into_duct(),
            time_step_s,
            exponential_tolerance,
            maximum_dofs,
        )?;
        load.wall_source = Some(source);
        Ok(load)
    }

    /// Prepare a quiescent load retaining one original source per main-bore
    /// element, including independently sourced Kelvin and Maxwell laws.
    pub fn new_with_regional_material_wall(
        wall: RegionalMaterialWallDg,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
    ) -> Result<Self, AcousticRealizeError> {
        let source = DgPortWallSource::Regional(wall.walls().to_vec());
        let mut load = Self::new(
            wall.into_duct(),
            time_step_s,
            exponential_tolerance,
            maximum_dofs,
        )?;
        load.wall_source = Some(source);
        Ok(load)
    }

    /// Retain a uniform original wall while attaching the existing causal
    /// axial terminal. Both owners keep the original same-step mean-flow law.
    pub fn new_with_material_wall_and_terminal(
        wall: MaterialWallDg,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        terminal: PreparedMimoLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let source = DgPortWallSource::Uniform(wall.wall().clone());
        let mut load = Self::new_with_terminal(
            wall.into_duct(),
            time_step_s,
            exponential_tolerance,
            maximum_dofs,
            terminal,
        )?;
        load.wall_source = Some(source);
        Ok(load)
    }

    /// Retain per-element original wall sources with the existing causal
    /// axial terminal, including every native wall and exterior pole history.
    pub fn new_with_regional_material_wall_and_terminal(
        wall: RegionalMaterialWallDg,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        terminal: PreparedMimoLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let source = DgPortWallSource::Regional(wall.walls().to_vec());
        let mut load = Self::new_with_terminal(
            wall.into_duct(),
            time_step_s,
            exponential_tolerance,
            maximum_dofs,
            terminal,
        )?;
        load.wall_source = Some(source);
        Ok(load)
    }

    /// Original source laws attached at material-aware construction, retained
    /// through accepted acoustic steps and subsequent native gas refreshes.
    /// A raw native duct does not invent missing material receipts.
    #[must_use]
    pub const fn material_wall_source(&self) -> Option<&DgPortWallSource> {
        self.wall_source.as_ref()
    }

    /// Refresh the retained main-bore wall at a supported new temperature.
    ///
    /// Every original source must identify the supplied immutable card and
    /// accepted initial temperature on the same absolute-temperature axis.
    /// Selected claims, non-temperature coordinates, thickness, restraint,
    /// material axes, fixed density and original loss-property roles survive.
    /// Regional sources retain their distinct local queries and receipts.
    ///
    /// The existing native transfer preserves radial displacement, conjugate
    /// momentum, gas/diffusion coordinates and admitted hereditary wall history.
    /// Its signed storage change is parameter work; this method consumes no
    /// sample and supplies no heat. A finite thermal owner must settle the work.
    /// Phase change, density/geometry changes and generated normal-plane
    /// viscosity still require their separate physical transfer model.
    ///
    /// Finite branches, finite terminal channels and exterior pole histories
    /// remain on their original sources and states. Their branch and terminal
    /// elimination is rebuilt against the refreshed main operator, including a
    /// driven mouth and sealed branches. All accepted interval frames survive.
    /// Source, native state, same-clock propagator and algebraic connections
    /// publish together after every source, work and cancellation check.
    #[allow(clippy::too_many_arguments)]
    pub fn rebind_wall_temperature_controlled(
        &mut self,
        card: &MaterialCard,
        temperature_axis: &str,
        initial_temperature_k: f64,
        temperature_k: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPortWallUpdate, AcousticRealizeError> {
        poll(cancelled)?;
        if !maximum_absolute_parameter_work_j.is_finite() || maximum_absolute_parameter_work_j < 0.0
        {
            return Err(invalid(
                "live wall temperature refresh requires a finite nonnegative parameter-work allowance",
            ));
        }
        let previous = self.wall_source.as_ref().ok_or_else(|| {
            invalid("live wall temperature refresh requires retained original material sources")
        })?;
        let source = previous.at_temperature_pinned(
            card,
            temperature_axis,
            initial_temperature_k,
            temperature_k,
            cancelled,
        )?;
        let elements = self.duct.mesh().elements;
        let previous_walls = previous.element_walls(elements);
        let walls = source.element_walls(elements);
        let retained_band_hz = previous_walls
            .iter()
            .map(LocalMaterialWall::band_hz)
            .reduce(|band, next| [band[0].max(next[0]), band[1].min(next[1])])
            .ok_or_else(|| invalid("live wall refresh lost its original source band"))?;
        let energy_before_j = self.stored_energy_j();
        let (duct, state, propagator, native) = LocalMaterialWall::prepare_dg_update(
            &previous_walls,
            &walls,
            retained_band_hz,
            &self.duct,
            &self.state,
            &self.propagator,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        if duct.total_dofs() != self.duct.total_dofs() {
            return Err(invalid(
                "live wall temperature refresh cannot change the admitted native state count",
            ));
        }
        poll(cancelled)?;
        let response = propagator.multiport_response(&state).map_err(duct_error)?;
        let branch_coupling = if self.branches.is_empty() {
            None
        } else {
            Some(BranchCoupling::new(&response, &self.branches)?)
        };
        let response = match &branch_coupling {
            Some(coupling) => coupling.response(response, &self.branches)?,
            None => response,
        };
        let terminal_coupling = match (&self.terminal, &self.terminal_coupling) {
            (Some(terminal), Some(previous)) => Some(TerminalCoupling::with_channels(
                &response,
                self.inlet,
                previous.ports.clone(),
                terminal,
                previous.driven_mouth,
                &self.terminal_channels,
            )?),
            (None, None) => None,
            _ => {
                return Err(invalid(
                    "live wall refresh lost the original terminal owner",
                ));
            }
        };
        let impedance = match (&self.terminal, &terminal_coupling) {
            (Some(terminal), Some(coupling)) => {
                coupling
                    .relation(response, terminal, &self.terminal_channels)?
                    .inlet
                    .impedance_pa_s_m3
            }
            _ if branch_coupling.is_none() => {
                propagator
                    .port_response(&state, self.inlet)
                    .map_err(duct_error)?
                    .impedance_pa_s_m3
            }
            _ => {
                if response.ports != [self.inlet] {
                    return Err(invalid(
                        "live wall inlet-only response lost sealed branch ownership",
                    ));
                }
                response.impedance_pa_s_m3[0]
            }
        };
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid(
                "live wall refresh requires strictly positive exact mean-port impedance",
            ));
        }
        let duct_energy = self
            .branches
            .iter()
            .fold(state.stored_energy_j(), |sum, branch| {
                sum + branch.state().stored_energy_j()
            });
        let energy_after_j = self
            .terminal_channels
            .iter()
            .fold(duct_energy, |sum, channel| sum + channel.stored_energy_j());
        if !energy_before_j.is_finite() || !energy_after_j.is_finite() {
            return Err(invalid(
                "live wall acoustic/channel storage cannot be represented",
            ));
        }
        let report = DgPortWallUpdate {
            source: source.clone(),
            native,
            energy_before_j,
            energy_after_j,
            parameter_work_j: native.parameter_work_j,
        };
        poll(cancelled)?;
        self.duct = duct;
        self.state = state;
        self.propagator = propagator;
        self.wall_source = Some(source);
        self.branch_coupling = branch_coupling;
        self.terminal_coupling = terminal_coupling;
        self.impedance_pa_s_m3 = impedance;
        Ok(report)
    }
}

impl DgPortTrial<'_> {
    /// Original MAIN-bore native loss, before branch/channel aggregation.
    /// The thermal transaction compares its main-wall heat against this owner
    /// with its explicit energy tolerance; unrelated losses cannot fund heat.
    pub(crate) fn main_dissipated_j(&self) -> f64 {
        self.state.cumulative_dissipated_j() - self.owner.state.cumulative_dissipated_j()
    }

    /// Original solid-material heat on this already solved held-flow interval.
    /// The exact inlet, loaded outlet and finite-branch junction inputs used by
    /// the native main step are retained in the trial. Gas, terminal, branch
    /// and numerical loss are not partitioned into this physical source.
    /// Continuous trajectories require their own integral and explicitly refuse.
    pub(crate) fn wall_dissipated_heat_j(
        &self,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AcousticRealizeError> {
        poll(cancelled)?;
        let inputs = self.main_inputs.as_deref().ok_or_else(|| {
            invalid("solid wall heat requires the original held-flow DG trajectory")
        })?;
        if self.owner.wall_source.is_none() {
            return Err(invalid(
                "solid wall heat requires retained original material sources",
            ));
        }
        let heat = self
            .owner
            .duct
            .prepare_compliant_wall_heat(&self.owner.propagator, budget, cancelled)
            .map_err(duct_error)?
            .integrated_heat_j(&self.owner.state, inputs, cancelled)
            .map_err(duct_error)?;
        if !heat.is_finite() || heat < 0.0 {
            return Err(invalid(
                "original solid wall heat must be finite and nonnegative",
            ));
        }
        poll(cancelled)?;
        Ok(heat)
    }
}
