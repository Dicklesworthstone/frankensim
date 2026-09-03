//! Browser boundary for the US 6,594,844 optical redirect composition.
//!
//! The generic differential-drive law and source-bounded Roomba composition
//! remain in `fs-mbd`. This crate admits one numeric packet and serializes one
//! stable success or typed-refusal JSON envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;

use fs_mbd::roomba::{
    ROOMBA_MAX_COLLIDERS, RectCollider, RoombaError, RoombaMode, RoombaRedirectReason, RoombaState,
    RoombaStepParams, step_roomba,
};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

const PACKET_VERSION: f64 = 1.0;
const HEADER_WORDS: usize = 18;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

fn parse_mode(value: f64) -> Option<RoombaMode> {
    match value {
        0.0 => Some(RoombaMode::Spiral),
        1.0 => Some(RoombaMode::Straight),
        2.0 => Some(RoombaMode::Turn),
        3.0 => Some(RoombaMode::Backup),
        _ => None,
    }
}

fn parse_bool(value: f64) -> Option<bool> {
    match value {
        0.0 => Some(false),
        1.0 => Some(true),
        _ => None,
    }
}

fn parse_seed(value: f64) -> Option<u32> {
    if value.is_finite() && value.fract() == 0.0 && (0.0..=f64::from(u32::MAX)).contains(&value) {
        Some(value as u32)
    } else {
        None
    }
}

fn serialize_success(result: RoombaState) -> String {
    let mut output = String::with_capacity(1_200);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"x_m\":{},\
         \"y_m\":{},\
         \"heading_rad\":{},\
         \"mode\":\"{}\",\
         \"time_in_mode_s\":{},\
         \"random_seed\":{},\
         \"optical_sensor_enabled\":{},\
         \"surface_overlap_fraction\":{},\
         \"surface_present\":{},\
         \"wall_present\":{},\
         \"redirect_reason\":\"{}\",\
         \"contact_index\":{},\
         \"contact_normal_x\":{},\
         \"contact_normal_y\":{},\
         \"left_wheel_speed_mps\":{},\
         \"right_wheel_speed_mps\":{},\
         \"left_wheel_angle_rad\":{},\
         \"right_wheel_angle_rad\":{},\
         \"side_brush_angle_rad\":{}\
         }}}}",
        result.x_m,
        result.y_m,
        result.heading_rad,
        result.mode.as_str(),
        result.time_in_mode_s,
        result.random_seed,
        result.optical_sensor_enabled,
        result.surface_overlap_fraction,
        result.surface_present,
        result.wall_present,
        result.redirect_reason.as_str(),
        result.contact_index,
        result.contact_normal_x,
        result.contact_normal_y,
        result.left_wheel_speed_mps,
        result.right_wheel_speed_mps,
        result.left_wheel_angle_rad,
        result.right_wheel_angle_rad,
        result.side_brush_angle_rad,
    );
    output
}

/// Step one explicit Roomba packet through the authoritative `fs-mbd` owner.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn roomba_step(packet: &[f64]) -> String {
    if packet.len() < HEADER_WORDS || !(packet.len() - HEADER_WORDS).is_multiple_of(4) {
        return refusal_json(
            "malformed-packet",
            "packet must contain 18 header words followed by collider quadruples",
            "Rebuild the version-1 packet using the documented field order",
        );
    }
    if packet[0] != PACKET_VERSION {
        return refusal_json(
            "unsupported-version",
            "only Roomba packet version 1 is supported",
            "Set packet word 0 to 1",
        );
    }
    let collider_count = (packet.len() - HEADER_WORDS) / 4;
    if collider_count > ROOMBA_MAX_COLLIDERS {
        return refusal_json(
            "resource-bound",
            "the bounded Roomba owner admits at most 64 colliders",
            "Reduce the low-solid collider receipt to at most 64 entries",
        );
    }
    if !packet.iter().all(|value| value.is_finite()) {
        return refusal_json(
            "non-finite-input",
            "every packet coordinate must be finite",
            "Replace NaN and infinity with finite declared values",
        );
    }
    let Some(optical_sensor_enabled) = parse_bool(packet[8]) else {
        return refusal_json(
            "malformed-packet",
            "optical subsystem flag must be exactly 0 or 1",
            "Encode the boolean as 0 or 1",
        );
    };
    let Some(mode) = parse_mode(packet[12]) else {
        return refusal_json(
            "malformed-packet",
            "mode must be 0 spiral, 1 straight, 2 turn, or 3 backup",
            "Encode one documented mode value",
        );
    };
    let Some(random_seed) = parse_seed(packet[14]) else {
        return refusal_json(
            "malformed-packet",
            "random seed must be an exact unsigned 32-bit integer",
            "Preserve the prior successful random_seed value",
        );
    };
    let wall_distance_inches = if packet[7] == -1.0 {
        None
    } else if packet[7] >= 0.0 {
        Some(packet[7])
    } else {
        return refusal_json(
            "malformed-packet",
            "wall distance must be nonnegative or exactly -1 for geometry-derived range",
            "Use -1 or a nonnegative distance in inches",
        );
    };

    let (collider_words, remainder) = packet[HEADER_WORDS..].as_chunks::<4>();
    debug_assert!(remainder.is_empty());
    let colliders: Vec<RectCollider> = collider_words
        .iter()
        .map(|values| RectCollider {
            x_m: values[0],
            y_m: values[1],
            width_m: values[2],
            height_m: values[3],
        })
        .collect();
    let result = step_roomba(
        RoombaStepParams {
            wheel_speed_mps: packet[2],
            turn_rate_rad_s: packet[3],
            room_width_m: packet[4],
            room_height_m: packet[5],
            sensor_height_inches: packet[6],
            wall_distance_inches,
            optical_sensor_enabled,
        },
        RoombaState {
            x_m: packet[9],
            y_m: packet[10],
            heading_rad: packet[11],
            mode,
            time_in_mode_s: packet[13],
            random_seed,
            optical_sensor_enabled,
            surface_overlap_fraction: 0.0,
            surface_present: false,
            wall_present: false,
            redirect_reason: RoombaRedirectReason::None,
            contact_index: -1,
            contact_normal_x: 0.0,
            contact_normal_y: 0.0,
            left_wheel_speed_mps: 0.0,
            right_wheel_speed_mps: 0.0,
            left_wheel_angle_rad: packet[15],
            right_wheel_angle_rad: packet[16],
            side_brush_angle_rad: packet[17],
        },
        &colliders,
        packet[1],
    );
    match result {
        Ok(result) => serialize_success(result),
        Err(RoombaError::InvalidInput(_)) => refusal_json(
            "input-outside-domain",
            "the source-bounded Roomba owner refused the declared state or geometry",
            "Use finite nonnegative controls, positive room/collider geometry, and dt no greater than 0.25 s",
        ),
        Err(RoombaError::TooManyColliders) => refusal_json(
            "resource-bound",
            "the bounded Roomba owner admits at most 64 colliders",
            "Reduce the low-solid collider receipt to at most 64 entries",
        ),
        Err(RoombaError::PlanarDrive(_)) => refusal_json(
            "multibody-refusal",
            "the generic planar differential-drive owner refused the step",
            "Inspect the fs-mbd planar-drive contract and admitted fixed-step geometry",
        ),
        Err(RoombaError::UnrepresentableOutput) => refusal_json(
            "unrepresentable-output",
            "the admitted step did not produce finite representable telemetry",
            "Reduce coordinate, speed, angle, or elapsed-time magnitude",
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_mbd::planar_drive::{DifferentialDriveStep, PlanarDriveState, step_differential_drive};

    fn packet() -> Vec<f64> {
        vec![
            1.0,
            1.0 / 120.0,
            0.3,
            1.5,
            4.0,
            4.0,
            0.5,
            -1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            42.0,
            0.0,
            0.0,
            0.0,
        ]
    }

    #[test]
    fn nominal_packet_returns_complete_success_tape() {
        let output = roomba_step(&packet());
        assert!(output.contains("\"ok\""));
        assert!(output.contains("\"mode\":\"spiral\""));
        assert!(output.contains("\"surface_overlap_fraction\":1"));
        assert!(output.contains("\"contact_index\":-1"));
        assert!(!output.contains("\"refusal\""));
    }

    #[test]
    fn optical_cliff_and_claim_inversion_are_distinct() {
        let mut cliff = packet();
        cliff[6] = 2.0;
        let redirected = roomba_step(&cliff);
        assert!(redirected.contains("\"redirect_reason\":\"surface-absent\""));
        assert!(redirected.contains("\"mode\":\"backup\""));

        cliff[8] = 0.0;
        let removed = roomba_step(&cliff);
        assert!(removed.contains("\"redirect_reason\":\"none\""));
        assert!(removed.contains("\"optical_sensor_enabled\":false"));
    }

    #[test]
    fn collider_packet_projects_an_embedded_chassis() {
        let mut input = packet();
        input.extend_from_slice(&[0.0, 0.0, 0.1, 0.1]);
        let output = roomba_step(&input);
        assert!(output.contains("\"contact_index\":0"));
        assert!(output.contains("\"mode\":\"backup\""));
    }

    #[test]
    fn malformed_non_finite_and_invalid_domain_packets_refuse() {
        assert!(roomba_step(&packet()[..17]).contains("\"code\":\"malformed-packet\""));
        let mut non_finite = packet();
        non_finite[2] = f64::NAN;
        assert!(roomba_step(&non_finite).contains("\"code\":\"non-finite-input\""));
        let mut invalid_room = packet();
        invalid_room[4] = 0.2;
        assert!(roomba_step(&invalid_room).contains("\"code\":\"input-outside-domain\""));
    }

    #[test]
    fn identical_packets_return_identical_bytes() {
        assert_eq!(roomba_step(&packet()), roomba_step(&packet()));
    }

    #[test]
    fn composed_generic_drive_matches_the_closed_form_arc_oracle() {
        let next = step_differential_drive(
            PlanarDriveState {
                x_m: 0.0,
                y_m: 0.0,
                heading_rad: 0.0,
                left_wheel_angle_rad: 0.0,
                right_wheel_angle_rad: 0.0,
            },
            DifferentialDriveStep {
                left_speed_mps: 0.1,
                right_speed_mps: 0.3,
                track_width_m: 0.2,
                wheel_radius_m: 0.05,
                dt_s: 0.2,
            },
        )
        .expect("valid analytic arc");
        assert!((next.heading_rad - 0.2).abs() < 1.0e-14);
        assert!((next.x_m - 0.2 * 0.2_f64.sin()).abs() < 1.0e-14);
        assert!((next.y_m - 0.2 * (1.0 - 0.2_f64.cos())).abs() < 1.0e-14);
    }

    #[test]
    fn oversized_collider_packet_refuses_before_owner_allocation() {
        let mut input = packet();
        for _ in 0..=ROOMBA_MAX_COLLIDERS {
            input.extend_from_slice(&[0.0, 0.0, 0.1, 0.1]);
        }
        assert!(roomba_step(&input).contains("\"code\":\"resource-bound\""));
    }
}
