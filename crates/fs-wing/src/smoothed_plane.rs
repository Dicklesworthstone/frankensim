//! SmoothedTangentPlane for rolling terrain ground effect (bead `frankensim-wf-root-guzez.5.11`, E4.4b).
//!
//! Filtered global ground plane with hysteresis and slope limits for non-flat terrain (V-06c).
//! Receipts carry plane velocity, filter state, and the ARTIFICIAL BOUNDARY-POWER residual.
//! Filter history is checkpointable for deterministic replay.
//!
//! NEVER inherits 'exact' (claim class is always EstimateOnly / Approximate).

use crate::{refuse, Refusal};

/// Default maximum allowable terrain slope [rad] (~8.5 degrees).
pub const DEFAULT_MAX_SLOPE_RAD: f64 = 0.15;
/// Default maximum vertical rate of the smoothed plane [m/s].
pub const DEFAULT_MAX_PLANE_RATE_M_S: f64 = 3.0;
/// Default hysteresis deadband [m].
pub const DEFAULT_HYSTERESIS_BAND_M: f64 = 0.05;
/// Default filter cutoff frequency [Hz].
pub const DEFAULT_FILTER_CUTOFF_HZ: f64 = 1.0;

/// Configuration for the smoothed tangent plane filter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmoothedTangentPlaneConfig {
    /// Maximum terrain slope admitted [rad].
    pub max_slope_rad: f64,
    /// Maximum allowable plane normal velocity [m/s].
    pub max_plane_rate_m_s: f64,
    /// Hysteresis deadband [m] to prevent plane chatter.
    pub hysteresis_band_m: f64,
    /// Low-pass filter cutoff frequency [Hz].
    pub filter_cutoff_hz: f64,
}

impl Default for SmoothedTangentPlaneConfig {
    fn default() -> Self {
        Self {
            max_slope_rad: DEFAULT_MAX_SLOPE_RAD,
            max_plane_rate_m_s: DEFAULT_MAX_PLANE_RATE_M_S,
            hysteresis_band_m: DEFAULT_HYSTERESIS_BAND_M,
            filter_cutoff_hz: DEFAULT_FILTER_CUTOFF_HZ,
        }
    }
}

/// Checkpointable internal filter state for the smoothed plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmoothedTangentPlaneState {
    /// Filtered plane height z [m] (positive down in NED).
    pub plane_z_m: f64,
    /// Filtered plane normal vector [nx, ny, nz] (unit vector pointing up: nz < 0).
    pub plane_normal: [f64; 3],
    /// Velocity of the plane height dz/dt [m/s].
    pub plane_z_dot_m_s: f64,
    /// Rate of change of plane normal d(normal)/dt [1/s].
    pub normal_dot: [f64; 3],
    /// Cumulative artificial work introduced by the moving boundary [J].
    pub cumulative_artificial_work_j: f64,
    /// Discrete time step counter.
    pub step_count: u64,
}

impl SmoothedTangentPlaneState {
    /// Create initial state at a given ground height and horizontal level normal.
    #[must_use]
    pub fn new(initial_z_m: f64) -> Self {
        Self {
            plane_z_m: initial_z_m,
            plane_normal: [0.0, 0.0, -1.0], // pointing up in NED (+z down)
            plane_z_dot_m_s: 0.0,
            normal_dot: [0.0, 0.0, 0.0],
            cumulative_artificial_work_j: 0.0,
            step_count: 0,
        }
    }
}

/// Receipt emitted on every plane update step (V-06c).
#[derive(Clone, Debug, PartialEq)]
pub struct SmoothedTangentPlaneReceipt {
    /// Filtered plane elevation z [m].
    pub plane_z_m: f64,
    /// Filtered unit normal vector.
    pub plane_normal: [f64; 3],
    /// Plane velocity vector [m/s].
    pub plane_velocity_m_s: [f64; 3],
    /// Instantaneous slope of the tangent plane [rad].
    pub slope_rad: f64,
    /// Instantaneous artificial boundary power residual [W].
    pub artificial_boundary_power_w: f64,
    /// Cumulative artificial boundary work [J].
    pub cumulative_artificial_work_j: f64,
    /// Authority/Claim class (always EstimateOnly; never Exact).
    pub claim_class: &'static str,
}

/// Dynamic filter engine for the smoothed tangent ground plane.
#[derive(Clone, Debug)]
pub struct SmoothedTangentPlane {
    config: SmoothedTangentPlaneConfig,
    state: SmoothedTangentPlaneState,
}

impl SmoothedTangentPlane {
    /// Construct a new filter engine with config and initial ground elevation.
    #[must_use]
    pub fn new(config: SmoothedTangentPlaneConfig, initial_z_m: f64) -> Self {
        Self {
            config,
            state: SmoothedTangentPlaneState::new(initial_z_m),
        }
    }

    /// Access current filter state (for checkpointing/snapshots).
    #[must_use]
    pub const fn state(&self) -> &SmoothedTangentPlaneState {
        &self.state
    }

    /// Restore filter state from a checkpoint.
    pub fn restore_state(&mut self, state: SmoothedTangentPlaneState) {
        self.state = state;
    }

    /// Update the smoothed ground plane from raw terrain probe coordinates and normal.
    ///
    /// # Errors
    /// [`Refusal`] if terrain slope exceeds `config.max_slope_rad` or inputs are non-finite.
    pub fn update(
        &mut self,
        raw_terrain_z_m: f64,
        raw_terrain_normal: [f64; 3],
        aircraft_induced_force_n: f64,
        dt_s: f64,
    ) -> Result<SmoothedTangentPlaneReceipt, Refusal> {
        if !raw_terrain_z_m.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(refuse(
                "smoothed-plane-non-finite",
                "terrain height or dt non-finite/non-positive".into(),
                "provide valid terrain probe and positive dt",
            ));
        }

        // Validate unit normal
        let n_sq = raw_terrain_normal[0] * raw_terrain_normal[0]
            + raw_terrain_normal[1] * raw_terrain_normal[1]
            + raw_terrain_normal[2] * raw_terrain_normal[2];
        if !n_sq.is_finite() || n_sq < 1e-6 {
            return Err(refuse(
                "smoothed-plane-degenerate-normal",
                "raw terrain normal has near-zero magnitude".into(),
                "provide valid unit normal vector",
            ));
        }
        let inv_n = 1.0 / n_sq.sqrt();
        let target_nx = raw_terrain_normal[0] * inv_n;
        let target_ny = raw_terrain_normal[1] * inv_n;
        let target_nz = raw_terrain_normal[2] * inv_n;

        // Calculate slope relative to horizontal plane
        let slope_rad = (target_nx * target_nx + target_ny * target_ny).sqrt().atan2(target_nz.abs());
        if slope_rad > self.config.max_slope_rad {
            return Err(refuse(
                "smoothed-plane-slope-exceeded",
                format!(
                    "terrain slope {:.4} rad exceeds maximum limit {:.4} rad",
                    slope_rad, self.config.max_slope_rad
                ),
                "terrain beyond smoothed tangent plane domain; use full 3D boundary element",
            ));
        }

        // Hysteresis deadband on target height
        let dz = raw_terrain_z_m - self.state.plane_z_m;
        let effective_target_z = if dz.abs() < self.config.hysteresis_band_m {
            self.state.plane_z_m // Keep current height inside deadband
        } else {
            raw_terrain_z_m - dz.signum() * self.config.hysteresis_band_m
        };

        // Low-pass filter dynamics (alpha = dt / (RC + dt))
        let tau = 1.0 / (2.0 * std::f64::consts::PI * self.config.filter_cutoff_hz);
        let alpha = (dt_s / (tau + dt_s)).clamp(0.0, 1.0);

        // Filter height with velocity rate limiting
        let desired_z = self.state.plane_z_m + alpha * (effective_target_z - self.state.plane_z_m);
        let raw_rate = (desired_z - self.state.plane_z_m) / dt_s;
        let clamped_rate = raw_rate.clamp(
            -self.config.max_plane_rate_m_s,
            self.config.max_plane_rate_m_s,
        );
        let new_z = self.state.plane_z_m + clamped_rate * dt_s;

        // Filter normal vector
        let new_nx = self.state.plane_normal[0] + alpha * (target_nx - self.state.plane_normal[0]);
        let new_ny = self.state.plane_normal[1] + alpha * (target_ny - self.state.plane_normal[1]);
        let new_nz = self.state.plane_normal[2] + alpha * (target_nz - self.state.plane_normal[2]);
        let new_norm = (new_nx * new_nx + new_ny * new_ny + new_nz * new_nz).sqrt();
        let unit_normal = [new_nx / new_norm, new_ny / new_norm, new_nz / new_norm];

        let normal_dot = [
            (unit_normal[0] - self.state.plane_normal[0]) / dt_s,
            (unit_normal[1] - self.state.plane_normal[1]) / dt_s,
            (unit_normal[2] - self.state.plane_normal[2]) / dt_s,
        ];

        // Artificial boundary power = induced force on plane * plane normal velocity
        let plane_vel = [0.0, 0.0, clamped_rate];
        let artificial_boundary_power_w = (aircraft_induced_force_n * clamped_rate).abs();
        let artificial_work_step_j = artificial_boundary_power_w * dt_s;

        // Update internal state
        self.state.plane_z_m = new_z;
        self.state.plane_normal = unit_normal;
        self.state.plane_z_dot_m_s = clamped_rate;
        self.state.normal_dot = normal_dot;
        self.state.cumulative_artificial_work_j += artificial_work_step_j;
        self.state.step_count += 1;

        Ok(SmoothedTangentPlaneReceipt {
            plane_z_m: new_z,
            plane_normal: unit_normal,
            plane_velocity_m_s: plane_vel,
            slope_rad,
            artificial_boundary_power_w,
            cumulative_artificial_work_j: self.state.cumulative_artificial_work_j,
            claim_class: "EstimateOnly",
        })
    }
}
