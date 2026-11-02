//! Full-tensor homogeneous Mann turbulence on an original NSF aperture.
//!
//! A retained snapshot is convected continuously using Taylor's frozen-field
//! hypothesis. Sampling at the actual CFD time avoids rounding to fs-atmo's
//! 120 Hz OU clock. The snapshot's own tick and phase stay fixed; only the
//! displacement relative to the caller's CFD reference time changes.
//!
//! Spatial synthesis is the existing fs-atmo Mann tensor quadrature (Mann,
//! 1998, doi:10.1016/S0266-8920(97)00036-2). Convected synthetic inflow is also
//! used by Nathan et al. (2018), doi:10.5194/wes-3-905-2018. This adapter adds
//! neither wall conditioning nor eddy turnover, mean subtraction, intensity
//! rescaling, thermodynamic fluctuations or nonreflecting characteristics.

use super::{NsfError, NsfGrid3, NsfOpenState, NsfReservoirTileState};
use crate::lattice3::Lattice3;
use fs_atmo::mann_field::HomogeneousMannField;
use fs_math::det;

/// Every flow, frame, clock and work choice is caller-authored.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfFrozenMannConfig {
    /// Uniform absolute reservoir pressure [Pa], resolved by the grid's gas.
    pub pressure_pa: f64,
    /// Uniform absolute reservoir temperature [K], resolved by the same gas.
    pub temperature_k: f64,
    /// Mean physical velocity in grid Cartesian axes [m/s].
    pub mean_velocity_m_s: [f64; 3],
    /// Frozen-field convection velocity in grid Cartesian axes [m/s].
    /// Usually the mean velocity under Taylor's hypothesis; no default is used.
    pub convection_velocity_m_s: [f64; 3],
    /// Coordinates of the grid's low domain corner in the Mann field frame [m].
    /// Cell centres lie at (index + 1/2)*dx relative to this corner.
    pub grid_origin_in_field_m: [f64; 3],
    /// Proper orthonormal rotation from Mann field components to grid components.
    /// Positions use its transpose; velocity fluctuations use this matrix.
    /// No NED/altitude sign convention or axis rotation is inferred.
    pub field_to_grid: [[f64; 3]; 3],
    /// CFD time [s] at which the supplied snapshot occupies its original coordinates.
    /// The snapshot's internal tick is retained, including its existing phase.
    pub reference_time_s: f64,
    /// Maximum aperture-tile count times represented wavevector count per update.
    pub maximum_mode_samples: usize,
}

/// Resolution and provenance of a successfully applied inlet state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfFrozenMannApplied {
    /// Actual accepted grid clock, steps*dt [s].
    pub sampled_time_s: f64,
    /// Unchanged source snapshot tick, not a rounded CFD time.
    pub frozen_tick: i64,
    pub tiles: usize,
    pub mode_samples: usize,
    /// Conservative radial-band phase bounds; each must be at most pi.
    pub spatial_phase_bound_rad: f64,
    pub temporal_phase_bound_rad: f64,
}

/// A refused update changes no aperture, fluid, source snapshot or clock.
#[derive(Clone, Debug, PartialEq)]
pub enum NsfFrozenMannError {
    Input(&'static str),
    Allocation,
    Cancelled,
    Work {
        requested: usize,
        maximum: usize,
    },
    /// k_max*dx > pi, with k_max=max_kl/L in rad/m.
    SpatialNyquist {
        phase_bound_rad: f64,
        maximum_rad: f64,
    },
    /// k_max*|U_conv|*dt > pi, a conservative bound on |k dot U_conv|*dt.
    TemporalNyquist {
        phase_bound_rad: f64,
        maximum_rad: f64,
    },
    Atmosphere(fs_atmo::Refusal),
    Boundary(NsfError),
}
impl core::fmt::Display for NsfFrozenMannError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frozen Mann inlet: {self:?}")
    }
}
impl std::error::Error for NsfFrozenMannError {}
impl From<NsfError> for NsfFrozenMannError {
    fn from(error: NsfError) -> Self {
        Self::Boundary(error)
    }
}
impl From<fs_atmo::Refusal> for NsfFrozenMannError {
    fn from(error: fs_atmo::Refusal) -> Self {
        Self::Atmosphere(error)
    }
}

/// Immutable full-tensor snapshot sampled into the original prescribed-state
/// boundary. No extra mass, energy or SGS kinetic-energy owner is introduced.
///
/// Call apply once at each accepted time boundary, then call the original
/// NsfGrid3::step. States remain held during that step, as for other explicit
/// reservoir controls. Repeating an update at the same time is deterministic
/// and changes no source RNG state. Updating and advancing the fluid are two
/// explicit operations; a later failed fluid step retains the already accepted
/// boundary control. This is not a second-order time-dependent boundary claim.
///
/// Only existing localized apertures are used. Their transverse wall topology
/// avoids falsely wrapping the arbitrary Mann wavevectors into a periodic
/// inlet plane. The finite aperture's sample mean need not vanish: changing it
/// would alter the represented covariance. The original accounts retain the
/// actual signed mass and total-energy exchanges, including fluctuations.
#[derive(Clone, Debug, PartialEq)]
pub struct NsfFrozenMannInflow {
    snapshot: HomogeneousMannField,
    config: NsfFrozenMannConfig,
    maximum_wavenumber_rad_m: f64,
    convection_speed_m_s: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl NsfFrozenMannInflow {
    /// Retain an already initialized source snapshot without cloning its
    /// spectral model or advancing its OU paths. Grid-specific material,
    /// spatial/temporal resolution and work bounds are admitted on apply.
    pub fn new(
        snapshot: HomogeneousMannField,
        config: NsfFrozenMannConfig,
    ) -> Result<Self, NsfFrozenMannError> {
        if !config.pressure_pa.is_finite()
            || config.pressure_pa <= 0.0
            || !config.temperature_k.is_finite()
            || config.temperature_k <= 0.0
            || !config.reference_time_s.is_finite()
            || config.maximum_mode_samples == 0
            || !config
                .mean_velocity_m_s
                .iter()
                .chain(&config.convection_velocity_m_s)
                .chain(&config.grid_origin_in_field_m)
                .chain(config.field_to_grid.iter().flatten())
                .all(|value| value.is_finite())
        {
            return Err(NsfFrozenMannError::Input(
                "invalid inlet state, frame, time or work bound",
            ));
        }
        let rotation = config.field_to_grid;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (dot(rotation[i], rotation[j]) - expected).abs() > 256.0 * f64::EPSILON {
                    return Err(NsfFrozenMannError::Input(
                        "field-to-grid rotation must be orthonormal",
                    ));
                }
            }
        }
        let determinant = rotation[0][0]
            * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
            - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
            + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]);
        if (determinant - 1.0).abs() > 512.0 * f64::EPSILON {
            return Err(NsfFrozenMannError::Input(
                "field-to-grid rotation must preserve handedness",
            ));
        }
        let maximum_wavenumber_rad_m =
            snapshot.model().quadrature().max_kl / snapshot.model().parameters().length_m;
        let scale = config
            .convection_velocity_m_s
            .iter()
            .fold(0.0_f64, |m, &u| m.max(u.abs()));
        let convection_speed_m_s = if scale == 0.0 {
            0.0
        } else {
            let normalized = config.convection_velocity_m_s.map(|u| u / scale);
            scale * det::sqrt(dot(normalized, normalized))
        };
        if !maximum_wavenumber_rad_m.is_finite()
            || maximum_wavenumber_rad_m <= 0.0
            || !convection_speed_m_s.is_finite()
        {
            return Err(NsfFrozenMannError::Input(
                "inlet wavenumber or convection speed is not representable",
            ));
        }
        Ok(Self {
            snapshot,
            config,
            maximum_wavenumber_rad_m,
            convection_speed_m_s,
        })
    }

    #[must_use]
    pub const fn config(&self) -> NsfFrozenMannConfig {
        self.config
    }

    /// The original immutable field, for checkpoint and represented-covariance queries.
    #[must_use]
    pub const fn snapshot(&self) -> &HomogeneousMannField {
        &self.snapshot
    }

    /// Sample at each original exterior ghost centre and the accepted CFD time,
    /// then atomically rebind one complete aperture profile. The density and
    /// thermal energy come from the grid's original gas at the authored p/T;
    /// every sampled velocity contributes its own resolved kinetic energy.
    pub fn apply<L: Lattice3<27>>(
        &self,
        grid: &mut NsfGrid3<L>,
        aperture: usize,
    ) -> Result<NsfFrozenMannApplied, NsfFrozenMannError> {
        self.apply_cancellable(grid, aperture, || false)
    }

    /// Bounded counterpart polling before allocation, between tile samples and
    /// before the original atomic boundary rebind. No random draws, spectral
    /// copies or source-clock changes occur. The original boundary rebind is
    /// synchronous; cancellation does not interrupt its internal material calls.
    pub fn apply_cancellable<L: Lattice3<27>>(
        &self,
        grid: &mut NsfGrid3<L>,
        aperture: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NsfFrozenMannApplied, NsfFrozenMannError> {
        if cancelled() {
            return Err(NsfFrozenMannError::Cancelled);
        }
        let source = grid
            .reservoir_aperture_at(aperture)
            .ok_or(NsfFrozenMannError::Input(
                "frozen Mann inlet requires an existing localized reservoir aperture",
            ))?;
        let axis = source.axis;
        let high = source.high;
        let count = source.cells.len();
        let mode_samples = count
            .checked_mul(self.snapshot.model().mode_count())
            .ok_or(NsfFrozenMannError::Input("inlet sample work overflows"))?;
        if mode_samples > self.config.maximum_mode_samples {
            return Err(NsfFrozenMannError::Work {
                requested: mode_samples,
                maximum: self.config.maximum_mode_samples,
            });
        }
        let spatial_phase_bound_rad = self.maximum_wavenumber_rad_m * grid.law.dx_m;
        if !spatial_phase_bound_rad.is_finite() || spatial_phase_bound_rad > core::f64::consts::PI {
            return Err(NsfFrozenMannError::SpatialNyquist {
                phase_bound_rad: spatial_phase_bound_rad,
                maximum_rad: core::f64::consts::PI,
            });
        }
        let temporal_phase_bound_rad =
            self.maximum_wavenumber_rad_m * (self.convection_speed_m_s * grid.law.dt_s);
        if !temporal_phase_bound_rad.is_finite() || temporal_phase_bound_rad > core::f64::consts::PI
        {
            return Err(NsfFrozenMannError::TemporalNyquist {
                phase_bound_rad: temporal_phase_bound_rad,
                maximum_rad: core::f64::consts::PI,
            });
        }
        // Beyond this, conversion of consecutive integer ticks can alias even
        // though the modal Nyquist bound is satisfied in real arithmetic.
        if grid.steps() as u128 >= (1u128 << 53) {
            return Err(NsfFrozenMannError::Input(
                "CFD step clock exceeds exact f64 integer range",
            ));
        }
        let time_s = grid.steps() as f64 * grid.law.dt_s;
        let elapsed_s = time_s - self.config.reference_time_s;
        if !time_s.is_finite()
            || !elapsed_s.is_finite()
            || (grid.steps() as f64 + 1.0) * grid.law.dt_s <= time_s
        {
            return Err(NsfFrozenMannError::Input(
                "CFD inlet time is not representable",
            ));
        }
        let gas = grid
            .spec
            .gas
            .at(self.config.temperature_k, self.config.pressure_pa)?;
        let mut states = Vec::new();
        states
            .try_reserve_exact(count)
            .map_err(|_| NsfFrozenMannError::Allocation)?;
        for tile in grid
            .reservoir_aperture_tile_states(aperture)
            .expect("checked aperture")
        {
            if cancelled() {
                return Err(NsfFrozenMannError::Cancelled);
            }
            let mut point_grid = tile.cell.map(|x| (x as f64 + 0.5) * grid.law.dx_m);
            point_grid[axis] += if high { grid.law.dx_m } else { -grid.law.dx_m };
            let displacement: [f64; 3] = core::array::from_fn(|a| {
                point_grid[a] - self.config.convection_velocity_m_s[a] * elapsed_s
            });
            let point_field = core::array::from_fn(|a| {
                self.config.grid_origin_in_field_m[a]
                    + (0..3)
                        .map(|i| self.config.field_to_grid[i][a] * displacement[i])
                        .sum::<f64>()
            });
            let fluctuation = self.snapshot.sample(point_field)?.u;
            let velocity_m_s = core::array::from_fn(|a| {
                self.config.mean_velocity_m_s[a] + dot(self.config.field_to_grid[a], fluctuation)
            });
            states.push(NsfReservoirTileState {
                cell: tile.cell,
                state: NsfOpenState {
                    density_kg_m3: gas.density,
                    velocity_m_s,
                    temperature_k: gas.temperature,
                },
            });
        }
        if cancelled() {
            return Err(NsfFrozenMannError::Cancelled);
        }
        grid.set_reservoir_aperture_profile(aperture, &states)?;
        Ok(NsfFrozenMannApplied {
            sampled_time_s: time_s,
            frozen_tick: self.snapshot.tick(),
            tiles: count,
            mode_samples,
            spatial_phase_bound_rad,
            temporal_phase_bound_rad,
        })
    }
}

#[cfg(test)]
mod tests;
