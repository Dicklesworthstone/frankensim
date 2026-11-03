//! CFD-clock OU evolution of the existing full-tensor Mann inlet.
//!
//! The atmosphere owns both temporal decorrelation and advection. A successful
//! apply publishes its advanced state together with the original reservoir
//! profile. Native fluid advancement remains a separate operation, as for the
//! frozen adapter. OU correlation times are authored model inputs, not a
//! resolved eddy-turnover or wall-layer model.

use super::*;
use fs_atmo::mann_field::fixed_step::FixedStepMannField;

/// Both inlet policies use the same physical-frame and boundary refusals.
pub type NsfEvolvingMannError = NsfFrozenMannError;

/// Explicit physical reservoir and work controls. The field's own advection
/// supplies the only temporal phase displacement; no second convection is added.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfEvolvingMannConfig {
    pub pressure_pa: f64,
    pub temperature_k: f64,
    pub mean_velocity_m_s: [f64; 3],
    pub grid_origin_in_field_m: [f64; 3],
    /// Proper orthonormal field-to-grid rotation.
    pub field_to_grid: [[f64; 3]; 3],
    /// Aperture tiles times represented wavevectors per apply.
    pub maximum_mode_samples: usize,
}

impl NsfEvolvingMannConfig {
    fn frame(self) -> NsfFrozenMannConfig {
        NsfFrozenMannConfig {
            pressure_pa: self.pressure_pa,
            temperature_k: self.temperature_k,
            mean_velocity_m_s: self.mean_velocity_m_s,
            convection_velocity_m_s: [0.0; 3],
            grid_origin_in_field_m: self.grid_origin_in_field_m,
            field_to_grid: self.field_to_grid,
            reference_time_s: 0.0,
            maximum_mode_samples: self.maximum_mode_samples,
        }
    }
}

/// The two explicit clocks and work actually used by a successful update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfEvolvingMannApplied {
    pub sampled_time_s: f64,
    pub atmosphere_step: u64,
    pub atmosphere_elapsed_time_s: f64,
    /// Original snapshot phase time plus the fixed-step elapsed time.
    pub atmosphere_time_s: f64,
    pub dt_s: f64,
    pub correlation_time_range_s: [f64; 2],
    pub ou_mode_transitions: usize,
    pub tiles: usize,
    pub mode_samples: usize,
    pub spatial_phase_bound_rad: f64,
    /// Bounds advection only. OU innovations are not temporally bandlimited.
    pub advective_phase_bound_rad: f64,
}

#[derive(Clone, Debug)]
struct GridBinding {
    dims: [usize; 3],
    dx_m: f64,
    dt_s: f64,
    aperture: usize,
    axis: usize,
    high: bool,
    cfd_step: usize,
    atmosphere_step: u64,
}

/// Reusable accepted/trial atmosphere states and one canonical profile buffer.
///
/// Construction binds the CURRENT grid step to the CURRENT atmosphere step.
/// Every later CFD step corresponds to exactly one atmosphere step at the same
/// admitted dt; no prehistory is invented for an already-running grid. A reused
/// grid must retain its dimensions, spacing, dt and selected aperture tiles.
///
/// Refused or cancelled applies change neither the accepted atmosphere nor the
/// grid. After a successful apply, a later failed native fluid step retains that
/// accepted boundary control. Reapplying at the same CFD step consumes no noise.
/// The original profile rebind and fluid solver may allocate their own work;
/// this adapter reuses its temporal state and sample buffer without spectral
/// cloning on each update.
#[derive(Debug)]
pub struct NsfEvolvingMannInflow {
    field: FixedStepMannField,
    trial: FixedStepMannField,
    config: NsfEvolvingMannConfig,
    binding: GridBinding,
    profile: Vec<NsfReservoirTileState>,
    averaging: Option<SpatialAverage>,
    wall: Option<MannPlaneWall>,
    maximum_wavenumber_rad_m: f64,
    advection_speed_m_s: f64,
    maximum_workspace_bytes: usize,
}

fn admit_bytes(required: usize, maximum: usize) -> Result<(), NsfEvolvingMannError> {
    if required > maximum {
        Err(NsfEvolvingMannError::WorkspaceBytes { required, maximum })
    } else {
        Ok(())
    }
}

impl NsfEvolvingMannInflow {
    /// Unfiltered numerical workspace for both temporal owners and the retained
    /// profile. Shared atmosphere parameter arrays are counted conservatively
    /// in both owners. Immutable spectral geometry and fixed metadata are
    /// excluded; prepared spatial-filter gains are added to `workspace_bytes`.
    pub fn required_workspace_bytes(
        field: &FixedStepMannField,
        tiles: usize,
    ) -> Result<usize, NsfEvolvingMannError> {
        field
            .workspace_bytes()
            .checked_mul(2)
            .and_then(|bytes| {
                tiles
                    .checked_mul(core::mem::size_of::<NsfReservoirTileState>())
                    .and_then(|profile| bytes.checked_add(profile))
            })
            .ok_or(NsfEvolvingMannError::Input(
                "evolving inlet workspace size overflows",
            ))
    }

    /// Bind the two current clock epochs without rounding or advancing either.
    /// The spectral source must already use this grid's exact physical dt.
    pub fn new<L: Lattice3<27>>(
        field: FixedStepMannField,
        config: NsfEvolvingMannConfig,
        grid: &NsfGrid3<L>,
        aperture: usize,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, NsfEvolvingMannError> {
        admit_config(config.frame())?;
        if field.dt_s() != grid.law.dt_s {
            return Err(NsfEvolvingMannError::Input(
                "atmosphere and CFD steps must have identical dt",
            ));
        }
        let source = grid
            .reservoir_aperture_at(aperture)
            .ok_or(NsfEvolvingMannError::Input(
                "evolving Mann inlet requires an existing aperture",
            ))?;
        let tiles = source.cells.len();
        let work = tiles
            .checked_mul(field.model().mode_count())
            .ok_or(NsfEvolvingMannError::Input("inlet sample work overflows"))?;
        if work > config.maximum_mode_samples {
            return Err(NsfEvolvingMannError::Work {
                requested: work,
                maximum: config.maximum_mode_samples,
            });
        }
        admit_bytes(
            Self::required_workspace_bytes(&field, tiles)?,
            maximum_workspace_bytes,
        )?;
        let maximum_wavenumber_rad_m =
            field.model().quadrature().max_kl / field.model().parameters().length_m;
        let advection_speed_m_s = speed(field.advection_mps());
        if !maximum_wavenumber_rad_m.is_finite()
            || maximum_wavenumber_rad_m <= 0.0
            || !advection_speed_m_s.is_finite()
        {
            return Err(NsfEvolvingMannError::Input(
                "inlet wavenumber or advection speed is not representable",
            ));
        }
        let mut profile = Vec::new();
        profile
            .try_reserve_exact(tiles)
            .map_err(|_| NsfEvolvingMannError::Allocation)?;
        profile.extend(
            grid.reservoir_aperture_tile_states(aperture)
                .expect("checked aperture"),
        );
        let binding = GridBinding {
            dims: grid.spec.dims,
            dx_m: grid.law.dx_m,
            dt_s: grid.law.dt_s,
            aperture,
            axis: source.axis,
            high: source.high,
            cfd_step: grid.steps(),
            atmosphere_step: field.step(),
        };
        let result = Self {
            trial: field.clone(),
            field,
            config,
            binding,
            profile,
            averaging: None,
            wall: None,
            maximum_wavenumber_rad_m,
            advection_speed_m_s,
            maximum_workspace_bytes,
        };
        admit_bytes(result.workspace_bytes(), maximum_workspace_bytes)?;
        Ok(result)
    }

    pub fn workspace_bytes(&self) -> usize {
        self.field.workspace_bytes()
            + self.trial.workspace_bytes()
            + self.profile.capacity() * core::mem::size_of::<NsfReservoirTileState>()
            + self
                .averaging
                .as_ref()
                .map_or(0, |a| a.filter.workspace_bytes())
    }

    /// Accepted atmosphere only; use its checkpoint for exact source restart.
    pub const fn field(&self) -> &FixedStepMannField {
        &self.field
    }
    pub const fn config(&self) -> NsfEvolvingMannConfig {
        self.config
    }
    pub fn spatial_filter(&self) -> Option<&MannBoxFilter> {
        self.averaging.as_ref().map(|a| &a.filter)
    }

    /// Observe the accepted/trial Mann realization through a harmonic slip
    /// plane in its own frame. The mean must be tangent; every ghost centre
    /// must be in the exact upper half-space. OU draws and clocks are unchanged.
    /// Combining this with box/face averaging refuses in either order.
    pub fn with_plane_wall(mut self, wall: MannPlaneWall) -> Result<Self, NsfEvolvingMannError> {
        admit_plane_wall(self.config.frame(), self.averaging.is_some())?;
        self.wall = Some(wall);
        Ok(self)
    }

    #[must_use]
    pub const fn plane_wall(&self) -> Option<MannPlaneWall> {
        self.wall
    }

    pub fn with_cell_average(
        self,
        maximum_filter_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NsfEvolvingMannError> {
        self.prepare_average(None, maximum_filter_bytes, cancelled)
    }

    pub fn with_face_average(
        self,
        maximum_filter_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NsfEvolvingMannError> {
        let axis = self.binding.axis;
        self.prepare_average(Some(axis), maximum_filter_bytes, cancelled)
    }

    fn prepare_average(
        mut self,
        face_axis: Option<usize>,
        maximum_filter_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NsfEvolvingMannError> {
        if self.wall.is_some() {
            return Err(NsfEvolvingMannError::Input(
                "Mann plane blocking currently requires point observations",
            ));
        }
        let remaining = self
            .maximum_workspace_bytes
            .saturating_sub(self.workspace_bytes());
        self.averaging = Some(prepare_spatial_average(
            self.field.model(),
            self.config.field_to_grid,
            self.binding.dx_m,
            face_axis,
            maximum_filter_bytes.min(remaining),
            cancelled,
        )?);
        admit_bytes(self.workspace_bytes(), self.maximum_workspace_bytes)?;
        Ok(self)
    }

    pub fn apply<L: Lattice3<27>>(
        &mut self,
        grid: &mut NsfGrid3<L>,
        maximum_ou_mode_transitions: usize,
    ) -> Result<NsfEvolvingMannApplied, NsfEvolvingMannError> {
        self.apply_cancellable(grid, maximum_ou_mode_transitions, || false)
    }

    /// Evolve a trial source, sample the actual ghost centres, and atomically
    /// rebind the existing gas/Guo-force/energy owner. No fallible operation
    /// follows the successful rebind. Cancellation is polled during temporal
    /// evolution and sampling, and immediately before the synchronous rebind.
    /// The transition budget counts four quadrature innovations per mode per
    /// elapsed step; zero permits an idempotent same-time apply.
    pub fn apply_cancellable<L: Lattice3<27>>(
        &mut self,
        grid: &mut NsfGrid3<L>,
        maximum_ou_mode_transitions: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NsfEvolvingMannApplied, NsfEvolvingMannError> {
        if cancelled() {
            return Err(NsfEvolvingMannError::Cancelled);
        }
        let b = &self.binding;
        let source = grid
            .reservoir_aperture_at(b.aperture)
            .ok_or(NsfEvolvingMannError::Input("bound Mann aperture is absent"))?;
        if grid.spec.dims != b.dims
            || grid.law.dx_m != b.dx_m
            || grid.law.dt_s != b.dt_s
            || source.axis != b.axis
            || source.high != b.high
            || source.cells.len() != self.profile.len()
            || !grid
                .reservoir_aperture_tile_states(b.aperture)
                .expect("checked aperture")
                .zip(&self.profile)
                .all(|(tile, old)| tile.cell == old.cell)
        {
            return Err(NsfEvolvingMannError::Input(
                "evolving Mann inlet requires its bound grid and aperture geometry",
            ));
        }
        let elapsed_steps = grid
            .steps()
            .checked_sub(b.cfd_step)
            .and_then(|steps| u64::try_from(steps).ok())
            .ok_or(NsfEvolvingMannError::Input(
                "CFD clock precedes the bound inlet epoch",
            ))?;
        let target = b
            .atmosphere_step
            .checked_add(elapsed_steps)
            .ok_or(NsfEvolvingMannError::Input("atmosphere step overflows"))?;
        let span = target
            .checked_sub(self.field.step())
            .and_then(|span| usize::try_from(span).ok())
            .ok_or(NsfEvolvingMannError::Input(
                "CFD clock precedes the accepted atmosphere",
            ))?;
        let transitions = span
            .checked_mul(self.field.model().mode_count())
            .and_then(|work| work.checked_mul(4))
            .ok_or(NsfEvolvingMannError::Input(
                "atmosphere transition work overflows",
            ))?;
        self.trial.copy_temporal_state_from(&self.field)?;
        self.trial
            .advance_to(target, maximum_ou_mode_transitions, &mut cancelled)
            .map_err(|error| {
                if error.code == "mann-fixed-step-cancelled" {
                    NsfEvolvingMannError::Cancelled
                } else {
                    NsfEvolvingMannError::Atmosphere(error)
                }
            })?;
        let sampled = sample_profile(
            grid,
            b.aperture,
            SamplingPlan {
                config: self.config.frame(),
                maximum_wavenumber_rad_m: self.maximum_wavenumber_rad_m,
                convection_speed_m_s: self.advection_speed_m_s,
                mode_count: self.trial.model().mode_count(),
                averaging: self.averaging.as_ref(),
            },
            &mut self.profile,
            |point, filter| match self.wall {
                Some(wall) => self.trial.sample_wall_blocked(point, wall),
                None => match filter {
                    Some(filter) => self.trial.sample_filtered(point, filter),
                    None => self.trial.sample(point),
                },
            },
            &mut cancelled,
        )?;
        let times = self.trial.correlation_times_s();
        let report = NsfEvolvingMannApplied {
            sampled_time_s: sampled.time_s,
            atmosphere_step: self.trial.step(),
            atmosphere_elapsed_time_s: self.trial.elapsed_time_s(),
            atmosphere_time_s: self.trial.time_s(),
            dt_s: b.dt_s,
            correlation_time_range_s: [
                times.iter().copied().fold(f64::INFINITY, f64::min),
                times.iter().copied().fold(0.0, f64::max),
            ],
            ou_mode_transitions: transitions,
            tiles: sampled.tiles,
            mode_samples: sampled.mode_samples,
            spatial_phase_bound_rad: sampled.spatial_phase_bound_rad,
            advective_phase_bound_rad: sampled.temporal_phase_bound_rad,
        };
        grid.set_reservoir_aperture_profile(b.aperture, &self.profile)?;
        core::mem::swap(&mut self.field, &mut self.trial);
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
