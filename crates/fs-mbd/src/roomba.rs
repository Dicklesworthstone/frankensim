//! Source-bounded optical redirect composition for US 6,594,844.
//!
//! The patented lane is the intersecting emitter/detector field and the
//! resulting change of travel direction when a floor surface is absent or a
//! nearby wall enters the modeled field. Chassis motion composes the generic
//! [`crate::planar_drive`] constant-twist owner. Room and furniture boundaries
//! use deterministic kinematic non-penetration projection for the museum
//! display; they are not a force, friction, impact, or contact-impulse model.

use core::fmt;

use crate::planar_drive::{
    DifferentialDriveStep, PlanarDriveError, PlanarDriveState, step_differential_drive,
};

/// Display bumper radius used by the shared Classic Patents room, in metres.
pub const ROOMBA_BUMPER_RADIUS_M: f64 = 0.17;
/// Distance between the two driven wheels in the procedural chassis, in metres.
pub const ROOMBA_TRACK_WIDTH_M: f64 = 0.24;
/// Effective procedural wheel radius used for visible rotation, in metres.
pub const ROOMBA_WHEEL_RADIUS_M: f64 = 0.035;
const INCHES_PER_METER: f64 = 39.370_078_740_157_48;
/// Maximum low-solid collider count admitted by the bounded museum kernel.
pub const ROOMBA_MAX_COLLIDERS: usize = 64;

/// Discrete contextual cleaning-path state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoombaMode {
    /// Expanding constant-speed spiral used as contextual navigation motion.
    Spiral,
    /// Straight contextual travel between redirects.
    Straight,
    /// In-place turn after a wall or bumper redirect.
    Turn,
    /// Short reverse interval after a cliff or bumper redirect.
    Backup,
}

impl RoombaMode {
    /// Stable browser-facing mode label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Spiral => "spiral",
            Self::Straight => "straight",
            Self::Turn => "turn",
            Self::Backup => "backup",
        }
    }
}

/// Why the patented optical subsystem requested a redirect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoombaRedirectReason {
    /// No optical redirect is currently requested.
    None,
    /// The downward fields no longer overlap on an expected floor surface.
    SurfaceAbsent,
    /// A nearby wall lies within the lateral field-intersection range.
    WallDetected,
}

impl RoombaRedirectReason {
    /// Stable browser-facing reason label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::SurfaceAbsent => "surface-absent",
            Self::WallDetected => "wall-detected",
        }
    }
}

/// One low rectangular obstacle footprint in the shared display room.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RectCollider {
    /// World x coordinate of the footprint center, in metres.
    pub x_m: f64,
    /// World y coordinate of the footprint center, in metres.
    pub y_m: f64,
    /// Footprint extent along world x, in metres.
    pub width_m: f64,
    /// Footprint extent along world y, in metres.
    pub height_m: f64,
}

/// Reader controls and declared environment for one Roomba step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoombaStepParams {
    /// Contextual chassis speed magnitude, in m/s.
    pub wheel_speed_mps: f64,
    /// Contextual in-place redirect rate, in rad/s.
    pub turn_rate_rad_s: f64,
    /// Shared display room width, in metres.
    pub room_width_m: f64,
    /// Shared display room height, in metres.
    pub room_height_m: f64,
    /// Downward optical sensor height above the expected surface, in inches.
    pub sensor_height_inches: f64,
    /// Optional caller-declared lateral wall distance, in inches. When absent,
    /// the room and collider geometry determine the range.
    pub wall_distance_inches: Option<f64>,
    /// Whether the Claim 1 emitter, detector, and redirect circuit are present.
    pub optical_sensor_enabled: bool,
}

/// Complete deterministic state shared by the museum's 2D and 3D faces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoombaState {
    /// Axle-midpoint world x coordinate, in metres.
    pub x_m: f64,
    /// Axle-midpoint world y coordinate, in metres.
    pub y_m: f64,
    /// Chassis heading, in radians.
    pub heading_rad: f64,
    /// Current contextual motion mode.
    pub mode: RoombaMode,
    /// Time committed in the current mode, in seconds.
    pub time_in_mode_s: f64,
    /// Deterministic pseudo-random state used only to vary turn duration.
    pub random_seed: u32,
    /// Whether the Claim 1 optical subsystem is present.
    pub optical_sensor_enabled: bool,
    /// Normalized overlap of the downward emitter and detector fields.
    pub surface_overlap_fraction: f64,
    /// Whether the overlap exceeds the display's source-facing threshold.
    pub surface_present: bool,
    /// Whether the lateral field sees a wall within its modeled range.
    pub wall_present: bool,
    /// Current patented optical redirect cause.
    pub redirect_reason: RoombaRedirectReason,
    /// `-1` for clear, `-2` for a room wall, otherwise the collider slice index.
    pub contact_index: i32,
    /// World x component of the kinematic projection normal.
    pub contact_normal_x: f64,
    /// World y component of the kinematic projection normal.
    pub contact_normal_y: f64,
    /// Current left prescribed wheel speed, in m/s.
    pub left_wheel_speed_mps: f64,
    /// Current right prescribed wheel speed, in m/s.
    pub right_wheel_speed_mps: f64,
    /// Left visible wheel rotation coordinate, in radians.
    pub left_wheel_angle_rad: f64,
    /// Right visible wheel rotation coordinate, in radians.
    pub right_wheel_angle_rad: f64,
    /// Side-brush display rotation coordinate, in radians.
    pub side_brush_angle_rad: f64,
}

/// Canonical initial state for the shared Roomba tape.
#[must_use]
pub const fn initial_roomba_state() -> RoombaState {
    RoombaState {
        x_m: 0.0,
        y_m: 0.0,
        heading_rad: 0.0,
        mode: RoombaMode::Spiral,
        time_in_mode_s: 0.0,
        random_seed: 42,
        optical_sensor_enabled: true,
        surface_overlap_fraction: 1.0,
        surface_present: true,
        wall_present: false,
        redirect_reason: RoombaRedirectReason::None,
        contact_index: -1,
        contact_normal_x: 0.0,
        contact_normal_y: 0.0,
        left_wheel_speed_mps: 0.0,
        right_wheel_speed_mps: 0.0,
        left_wheel_angle_rad: 0.0,
        right_wheel_angle_rad: 0.0,
        side_brush_angle_rad: 0.0,
    }
}

/// Typed refusal from the US 6,594,844 composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoombaError {
    /// A named input lies outside the admitted finite geometric domain.
    InvalidInput(&'static str),
    /// The bounded museum kernel refuses an unbounded obstacle list.
    TooManyColliders,
    /// The generic differential-drive owner refused the prescribed step.
    PlanarDrive(PlanarDriveError),
    /// Finite admitted input produced an unrepresentable output.
    UnrepresentableOutput,
}

impl From<PlanarDriveError> for RoombaError {
    fn from(value: PlanarDriveError) -> Self {
        Self::PlanarDrive(value)
    }
}

impl fmt::Display for RoombaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(field) => write!(formatter, "invalid Roomba input: {field}"),
            Self::TooManyColliders => {
                write!(
                    formatter,
                    "at most {ROOMBA_MAX_COLLIDERS} colliders are admitted"
                )
            }
            Self::PlanarDrive(error) => write!(formatter, "planar-drive refusal: {error}"),
            Self::UnrepresentableOutput => {
                formatter.write_str("Roomba output is not representable as finite f64")
            }
        }
    }
}

impl std::error::Error for RoombaError {}

#[derive(Debug, Clone, Copy)]
struct Projection {
    x_m: f64,
    y_m: f64,
    hit: bool,
    normal_x: f64,
    normal_y: f64,
}

fn valid_collider(collider: RectCollider) -> bool {
    [
        collider.x_m,
        collider.y_m,
        collider.width_m,
        collider.height_m,
    ]
    .into_iter()
    .all(f64::is_finite)
        && collider.width_m > 0.0
        && collider.height_m > 0.0
}

fn project_outside_collider(x_m: f64, y_m: f64, collider: RectCollider) -> Projection {
    let min_x = collider.x_m - collider.width_m / 2.0;
    let max_x = collider.x_m + collider.width_m / 2.0;
    let min_y = collider.y_m - collider.height_m / 2.0;
    let max_y = collider.y_m + collider.height_m / 2.0;
    let closest_x = x_m.clamp(min_x, max_x);
    let closest_y = y_m.clamp(min_y, max_y);
    let dx = x_m - closest_x;
    let dy = y_m - closest_y;
    let distance = dx.hypot(dy);

    if distance >= ROOMBA_BUMPER_RADIUS_M - 1.0e-9 {
        return Projection {
            x_m,
            y_m,
            hit: false,
            normal_x: 0.0,
            normal_y: 0.0,
        };
    }
    if distance > 1.0e-9 {
        let normal_x = dx / distance;
        let normal_y = dy / distance;
        return Projection {
            x_m: closest_x + normal_x * ROOMBA_BUMPER_RADIUS_M,
            y_m: closest_y + normal_y * ROOMBA_BUMPER_RADIUS_M,
            hit: true,
            normal_x,
            normal_y,
        };
    }

    let exits = [
        (x_m - min_x, min_x - ROOMBA_BUMPER_RADIUS_M, y_m, -1.0, 0.0),
        (max_x - x_m, max_x + ROOMBA_BUMPER_RADIUS_M, y_m, 1.0, 0.0),
        (y_m - min_y, x_m, min_y - ROOMBA_BUMPER_RADIUS_M, 0.0, -1.0),
        (max_y - y_m, x_m, max_y + ROOMBA_BUMPER_RADIUS_M, 0.0, 1.0),
    ];
    let mut exit = exits[0];
    for candidate in &exits[1..] {
        if candidate.0.total_cmp(&exit.0).is_lt() {
            exit = *candidate;
        }
    }
    Projection {
        x_m: exit.1,
        y_m: exit.2,
        hit: true,
        normal_x: exit.3,
        normal_y: exit.4,
    }
}

fn ray_distance_to_expanded_collider(
    x_m: f64,
    y_m: f64,
    direction_x: f64,
    direction_y: f64,
    collider: RectCollider,
) -> f64 {
    let min_x = collider.x_m - collider.width_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    let max_x = collider.x_m + collider.width_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let min_y = collider.y_m - collider.height_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    let max_y = collider.y_m + collider.height_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let mut entry: f64 = 0.0;
    let mut exit = f64::INFINITY;

    for (origin, direction, min, max) in [
        (x_m, direction_x, min_x, max_x),
        (y_m, direction_y, min_y, max_y),
    ] {
        if direction.abs() < 1.0e-12 {
            if origin < min || origin > max {
                return f64::INFINITY;
            }
            continue;
        }
        let first = (min - origin) / direction;
        let second = (max - origin) / direction;
        entry = entry.max(first.min(second));
        exit = exit.min(first.max(second));
        if exit < entry {
            return f64::INFINITY;
        }
    }
    if exit < 0.0 {
        f64::INFINITY
    } else {
        entry.max(0.0)
    }
}

fn modeled_wall_distance_inches(
    state: RoombaState,
    params: RoombaStepParams,
    colliders: &[RectCollider],
) -> f64 {
    let direction_x = state.heading_rad.cos();
    let direction_y = state.heading_rad.sin();
    let min_x = -params.room_width_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let max_x = params.room_width_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    let min_y = -params.room_height_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let max_y = params.room_height_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    let wall_x = if direction_x > 1.0e-12 {
        (max_x - state.x_m) / direction_x
    } else if direction_x < -1.0e-12 {
        (min_x - state.x_m) / direction_x
    } else {
        f64::INFINITY
    };
    let wall_y = if direction_y > 1.0e-12 {
        (max_y - state.y_m) / direction_y
    } else if direction_y < -1.0e-12 {
        (min_y - state.y_m) / direction_y
    } else {
        f64::INFINITY
    };
    let mut nearest = wall_x.min(wall_y);
    for collider in colliders {
        nearest = nearest.min(ray_distance_to_expanded_collider(
            state.x_m,
            state.y_m,
            direction_x,
            direction_y,
            *collider,
        ));
    }
    nearest.max(0.0) * INCHES_PER_METER
}

fn validate_inputs(
    params: RoombaStepParams,
    state: RoombaState,
    colliders: &[RectCollider],
    dt_s: f64,
) -> Result<(), RoombaError> {
    if colliders.len() > ROOMBA_MAX_COLLIDERS {
        return Err(RoombaError::TooManyColliders);
    }
    if !colliders.iter().copied().all(valid_collider) {
        return Err(RoombaError::InvalidInput("collider geometry"));
    }
    if !params.wheel_speed_mps.is_finite() || params.wheel_speed_mps < 0.0 {
        return Err(RoombaError::InvalidInput("wheel_speed_mps"));
    }
    if !params.turn_rate_rad_s.is_finite() || params.turn_rate_rad_s < 0.0 {
        return Err(RoombaError::InvalidInput("turn_rate_rad_s"));
    }
    if !params.room_width_m.is_finite()
        || params.room_width_m <= 2.0 * ROOMBA_BUMPER_RADIUS_M
        || !params.room_height_m.is_finite()
        || params.room_height_m <= 2.0 * ROOMBA_BUMPER_RADIUS_M
    {
        return Err(RoombaError::InvalidInput("room dimensions"));
    }
    if !params.sensor_height_inches.is_finite() || params.sensor_height_inches < 0.0 {
        return Err(RoombaError::InvalidInput("sensor_height_inches"));
    }
    if params
        .wall_distance_inches
        .is_some_and(|distance| !distance.is_finite() || distance < 0.0)
    {
        return Err(RoombaError::InvalidInput("wall_distance_inches"));
    }
    if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > 0.25 {
        return Err(RoombaError::InvalidInput("dt_s"));
    }
    if ![
        state.x_m,
        state.y_m,
        state.heading_rad,
        state.time_in_mode_s,
        state.left_wheel_angle_rad,
        state.right_wheel_angle_rad,
        state.side_brush_angle_rad,
    ]
    .into_iter()
    .all(f64::is_finite)
        || state.time_in_mode_s < 0.0
    {
        return Err(RoombaError::InvalidInput("state"));
    }
    Ok(())
}

/// Step the patented optical redirect and contextual differential-drive tape.
///
/// # Errors
/// Refuses invalid or non-finite controls/state/geometry, more than 64 low
/// colliders, a generic planar-drive refusal, or non-finite derived output.
pub fn step_roomba(
    params: RoombaStepParams,
    mut state: RoombaState,
    colliders: &[RectCollider],
    dt_s: f64,
) -> Result<RoombaState, RoombaError> {
    validate_inputs(params, state, colliders, dt_s)?;

    let wall_distance_inches = params
        .wall_distance_inches
        .unwrap_or_else(|| modeled_wall_distance_inches(state, params, colliders));
    let surface_overlap_fraction =
        (1.0 - (params.sensor_height_inches - 0.5).abs() / 0.5).clamp(0.0, 1.0);
    let surface_present = surface_overlap_fraction > 0.2;
    let wall_present = params.optical_sensor_enabled && wall_distance_inches <= 2.95;
    state.optical_sensor_enabled = params.optical_sensor_enabled;
    state.surface_overlap_fraction = surface_overlap_fraction;
    state.surface_present = surface_present;
    state.wall_present = wall_present;
    state.redirect_reason = if !params.optical_sensor_enabled {
        RoombaRedirectReason::None
    } else if !surface_present {
        RoombaRedirectReason::SurfaceAbsent
    } else if wall_present {
        RoombaRedirectReason::WallDetected
    } else {
        RoombaRedirectReason::None
    };
    state.time_in_mode_s += dt_s;

    let outside_room = state.x_m > params.room_width_m / 2.0 - ROOMBA_BUMPER_RADIUS_M
        || state.x_m < -params.room_width_m / 2.0 + ROOMBA_BUMPER_RADIUS_M
        || state.y_m > params.room_height_m / 2.0 - ROOMBA_BUMPER_RADIUS_M
        || state.y_m < -params.room_height_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let embedded_in_collider = colliders
        .iter()
        .copied()
        .any(|collider| project_outside_collider(state.x_m, state.y_m, collider).hit);
    if ((params.optical_sensor_enabled && !surface_present) || outside_room || embedded_in_collider)
        && !matches!(state.mode, RoombaMode::Backup | RoombaMode::Turn)
    {
        state.mode = RoombaMode::Backup;
        state.time_in_mode_s = 0.0;
    }
    if wall_present && !matches!(state.mode, RoombaMode::Backup | RoombaMode::Turn) {
        state.mode = RoombaMode::Turn;
        state.time_in_mode_s = 0.0;
    }
    if state.mode == RoombaMode::Backup && state.time_in_mode_s > 0.4 {
        state.mode = RoombaMode::Turn;
        state.time_in_mode_s = 0.0;
        state.random_seed = state
            .random_seed
            .wrapping_mul(1_103_515_245)
            .wrapping_add(12_345)
            & 0x7fff_ffff;
    } else if state.mode == RoombaMode::Turn {
        let turn_duration = 0.4 + f64::from(state.random_seed % 100) / 100.0 * 1.2;
        if state.time_in_mode_s > turn_duration {
            state.mode = RoombaMode::Straight;
            state.time_in_mode_s = 0.0;
        }
    }

    let (mut left_speed_mps, mut right_speed_mps) = match state.mode {
        RoombaMode::Spiral => {
            let radius_m = 0.12 + state.time_in_mode_s * 0.045;
            let yaw_rate_rad_s = params.wheel_speed_mps / radius_m;
            (
                params.wheel_speed_mps - yaw_rate_rad_s * ROOMBA_TRACK_WIDTH_M / 2.0,
                params.wheel_speed_mps + yaw_rate_rad_s * ROOMBA_TRACK_WIDTH_M / 2.0,
            )
        }
        RoombaMode::Straight => (params.wheel_speed_mps, params.wheel_speed_mps),
        RoombaMode::Backup => (-params.wheel_speed_mps, -params.wheel_speed_mps),
        RoombaMode::Turn => (
            -params.turn_rate_rad_s * ROOMBA_TRACK_WIDTH_M / 2.0,
            params.turn_rate_rad_s * ROOMBA_TRACK_WIDTH_M / 2.0,
        ),
    };
    let drive = step_differential_drive(
        PlanarDriveState {
            x_m: state.x_m,
            y_m: state.y_m,
            heading_rad: state.heading_rad,
            left_wheel_angle_rad: state.left_wheel_angle_rad,
            right_wheel_angle_rad: state.right_wheel_angle_rad,
        },
        DifferentialDriveStep {
            left_speed_mps,
            right_speed_mps,
            track_width_m: ROOMBA_TRACK_WIDTH_M,
            wheel_radius_m: ROOMBA_WHEEL_RADIUS_M,
            dt_s,
        },
    )?;
    state.x_m = drive.x_m;
    state.y_m = drive.y_m;
    state.heading_rad = drive.heading_rad;
    state.left_wheel_angle_rad = drive.left_wheel_angle_rad;
    state.right_wheel_angle_rad = drive.right_wheel_angle_rad;

    let min_x = -params.room_width_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let max_x = params.room_width_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    let min_y = -params.room_height_m / 2.0 + ROOMBA_BUMPER_RADIUS_M;
    let max_y = params.room_height_m / 2.0 - ROOMBA_BUMPER_RADIUS_M;
    state.contact_index = -1;
    state.contact_normal_x = 0.0;
    state.contact_normal_y = 0.0;
    if state.x_m < min_x {
        state.x_m = min_x;
        state.contact_index = -2;
        state.contact_normal_x = 1.0;
    } else if state.x_m > max_x {
        state.x_m = max_x;
        state.contact_index = -2;
        state.contact_normal_x = -1.0;
    }
    if state.y_m < min_y {
        state.y_m = min_y;
        state.contact_index = -2;
        state.contact_normal_x = 0.0;
        state.contact_normal_y = 1.0;
    } else if state.y_m > max_y {
        state.y_m = max_y;
        state.contact_index = -2;
        state.contact_normal_x = 0.0;
        state.contact_normal_y = -1.0;
    }

    for (index, collider) in colliders.iter().copied().enumerate() {
        let projection = project_outside_collider(state.x_m, state.y_m, collider);
        if !projection.hit {
            continue;
        }
        state.x_m = projection.x_m;
        state.y_m = projection.y_m;
        state.contact_index = i32::try_from(index).map_err(|_| RoombaError::TooManyColliders)?;
        state.contact_normal_x = projection.normal_x;
        state.contact_normal_y = projection.normal_y;
    }
    if state.contact_index != -1 && !matches!(state.mode, RoombaMode::Backup | RoombaMode::Turn) {
        state.mode = RoombaMode::Backup;
        state.time_in_mode_s = 0.0;
        left_speed_mps = -params.wheel_speed_mps;
        right_speed_mps = -params.wheel_speed_mps;
    }

    state.left_wheel_speed_mps = left_speed_mps;
    state.right_wheel_speed_mps = right_speed_mps;
    state.side_brush_angle_rad +=
        60.0 * 0.5 * (left_speed_mps.abs() + right_speed_mps.abs()) * dt_s;
    if [
        state.x_m,
        state.y_m,
        state.heading_rad,
        state.time_in_mode_s,
        state.surface_overlap_fraction,
        state.contact_normal_x,
        state.contact_normal_y,
        state.left_wheel_speed_mps,
        state.right_wheel_speed_mps,
        state.left_wheel_angle_rad,
        state.right_wheel_angle_rad,
        state.side_brush_angle_rad,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        Ok(state)
    } else {
        Err(RoombaError::UnrepresentableOutput)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> RoombaStepParams {
        RoombaStepParams {
            wheel_speed_mps: 0.3,
            turn_rate_rad_s: 1.5,
            room_width_m: 4.0,
            room_height_m: 4.0,
            sensor_height_inches: 0.5,
            wall_distance_inches: None,
            optical_sensor_enabled: true,
        }
    }

    #[test]
    fn nominal_step_composes_the_generic_drive_and_optical_overlap() {
        let next = step_roomba(params(), initial_roomba_state(), &[], 1.0 / 120.0)
            .expect("valid Roomba step");
        assert_eq!(next.mode, RoombaMode::Spiral);
        assert!(next.x_m > 0.0);
        assert!(next.heading_rad > 0.0);
        assert_eq!(next.surface_overlap_fraction, 1.0);
        assert!(next.surface_present);
        assert_eq!(next.redirect_reason, RoombaRedirectReason::None);
    }

    #[test]
    fn absent_surface_triggers_the_patented_backup_redirect() {
        let next = step_roomba(
            RoombaStepParams {
                sensor_height_inches: 2.0,
                ..params()
            },
            initial_roomba_state(),
            &[],
            1.0 / 120.0,
        )
        .expect("valid cliff probe");
        assert!(!next.surface_present);
        assert_eq!(next.redirect_reason, RoombaRedirectReason::SurfaceAbsent);
        assert_eq!(next.mode, RoombaMode::Backup);
        assert!(next.left_wheel_speed_mps < 0.0);
        assert!(next.right_wheel_speed_mps < 0.0);
    }

    #[test]
    fn claim_inversion_removes_optical_redirect_but_not_room_projection() {
        let mut state = initial_roomba_state();
        state.x_m = 1.9;
        state.mode = RoombaMode::Straight;
        let next = step_roomba(
            RoombaStepParams {
                sensor_height_inches: 2.0,
                optical_sensor_enabled: false,
                ..params()
            },
            state,
            &[],
            1.0 / 120.0,
        )
        .expect("valid claim inversion");
        assert_eq!(next.redirect_reason, RoombaRedirectReason::None);
        assert!(!next.wall_present);
        assert_eq!(next.contact_index, -2);
        assert!(next.x_m <= 2.0 - ROOMBA_BUMPER_RADIUS_M);
        assert_eq!(next.mode, RoombaMode::Backup);
    }

    #[test]
    fn embedded_chassis_exits_the_nearest_collider_face() {
        let collider = RectCollider {
            x_m: 0.0,
            y_m: 0.0,
            width_m: 0.1,
            height_m: 0.1,
        };
        let mut state = initial_roomba_state();
        state.mode = RoombaMode::Straight;
        let next = step_roomba(params(), state, &[collider], 1.0 / 120.0)
            .expect("valid obstacle projection");
        assert_eq!(next.contact_index, 0);
        assert!((next.contact_normal_x.hypot(next.contact_normal_y) - 1.0).abs() < 1.0e-12);
        assert_eq!(next.mode, RoombaMode::Backup);
        let second = project_outside_collider(next.x_m, next.y_m, collider);
        assert!(!second.hit);
    }

    #[test]
    fn invalid_geometry_time_and_unbounded_collider_lists_refuse() {
        assert_eq!(
            step_roomba(
                RoombaStepParams {
                    room_width_m: 0.2,
                    ..params()
                },
                initial_roomba_state(),
                &[],
                1.0 / 120.0,
            ),
            Err(RoombaError::InvalidInput("room dimensions"))
        );
        assert_eq!(
            step_roomba(params(), initial_roomba_state(), &[], 0.5),
            Err(RoombaError::InvalidInput("dt_s"))
        );
        let colliders = vec![
            RectCollider {
                x_m: 0.0,
                y_m: 0.0,
                width_m: 0.1,
                height_m: 0.1
            };
            65
        ];
        assert_eq!(
            step_roomba(params(), initial_roomba_state(), &colliders, 1.0 / 120.0),
            Err(RoombaError::TooManyColliders)
        );
    }

    #[test]
    fn replay_is_bit_deterministic() {
        let mut first = initial_roomba_state();
        let mut second = initial_roomba_state();
        for _ in 0..1_000 {
            first = step_roomba(params(), first, &[], 1.0 / 120.0).expect("first replay");
            second = step_roomba(params(), second, &[], 1.0 / 120.0).expect("second replay");
        }
        assert_eq!(first, second);
    }
}
