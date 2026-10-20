//! Shared transaction and work composition for single and paired reed sources.
//! The original aperture and graph owners still perform every physical solve.

use super::boundaries::{BoundaryAperture, BoundaryApertures, BoundaryDrive, BoundaryTrial};
use super::{ApertureNetworkFrame, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureFrame;
use crate::bernoulli_aperture::tube::TubeDrive;
use fs_vfit::waveguide::network::WaveguideNetwork;

pub(in crate::bernoulli_aperture) struct NetworkTrial {
    pub frame: ApertureNetworkFrame,
    boundaries: Option<BoundaryTrial>,
}

pub(in crate::bernoulli_aperture) fn admit_horizon(
    inputs: &[BoundaryAperture],
    max_apertures: usize,
    max_steps: u64,
) -> Result<(), AcousticRealizeError> {
    if !(1..=64).contains(&max_apertures)
        || inputs.is_empty()
        || inputs.len() > max_apertures
        || inputs
            .iter()
            .any(|b| b.aperture.spec().max_steps != max_steps)
    {
        return Err(invalid(
            "reed and moving boundaries require one shared sample horizon and bounded aperture count",
        ));
    }
    Ok(())
}

pub(in crate::bernoulli_aperture) fn preview(
    aperture: ApertureFrame,
    drive: TubeDrive,
    dt: f64,
    network: &mut WaveguideNetwork,
    boundaries: Option<&BoundaryApertures>,
    drives: &[BoundaryDrive],
) -> Result<NetworkTrial, AcousticRealizeError> {
    let boundary_trial = match boundaries {
        Some(boundaries) => {
            Some(boundaries.preview_step(network, aperture.outgoing_pressure_pa, drives)?)
        }
        None if drives.is_empty() => None,
        None => {
            return Err(invalid(
                "boundary drives require previously bound moving apertures",
            ));
        }
    };
    let boundary_frame = boundary_trial.as_ref().map(|trial| trial.frame);
    let network = match boundary_frame {
        Some(frame) => frame.network,
        None => network
            .preview_step(aperture.outgoing_pressure_pa)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?,
    };
    // Preserve the original arithmetic exactly when no moving keys are bound.
    // The boundary frame already owns graph storage and dissipation; do not add
    // either a second time when composing the primary reed with those keys.
    let (stored_energy_j, storage_change_j, dissipated_energy_j) = match boundary_frame {
        Some(frame) => (
            aperture.stored_energy_j + frame.stored_energy_j,
            aperture.storage_change_j + frame.storage_change_j,
            aperture.dissipated_energy_j + frame.dissipated_energy_j,
        ),
        None => (
            aperture.stored_energy_j + network.stored_energy_j,
            aperture.storage_change_j + network.storage_change_j,
            aperture.dissipated_energy_j + network.terminal_loss_j + network.interior_loss_j,
        ),
    };
    let frame = ApertureNetworkFrame {
        aperture,
        network,
        boundaries: boundary_frame,
        stored_energy_j,
        storage_change_j,
        dissipated_energy_j,
        upstream_work_j: drive.upstream_pressure_pa
            * ((aperture.bore_flow_m3_s - drive.body_flow_m3_s) * dt),
        body_work_j: aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt),
    };
    if ![
        frame.stored_energy_j,
        frame.storage_change_j,
        frame.dissipated_energy_j,
        frame.upstream_work_j,
        frame.body_work_j,
        frame.balance_residual_j(),
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(invalid("coupled network observation left the finite set"));
    }
    Ok(NetworkTrial {
        frame,
        boundaries: boundary_trial,
    })
}

// The caller completes source energy/cancellation checks before this function.
// Graph validation is the last fallible physical operation; key and primary
// aperture acceptance are infallible and reuse their original candidate states.
pub(in crate::bernoulli_aperture) fn accept(
    network: &mut WaveguideNetwork,
    boundaries: Option<&mut BoundaryApertures>,
    trial: NetworkTrial,
) -> Result<ApertureNetworkFrame, AcousticRealizeError> {
    match (boundaries, trial.boundaries) {
        (Some(boundaries), Some(boundary_trial)) => {
            boundaries.accept_trial(network, boundary_trial)?;
        }
        (None, None) => {
            network
                .step(trial.frame.aperture.outgoing_pressure_pa)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        }
        _ => {
            return Err(invalid(
                "moving-boundary trial must retain its original owner",
            ));
        }
    }
    Ok(trial.frame)
}

pub(in crate::bernoulli_aperture) fn validate_step_budget(
    boundaries: Option<&BoundaryApertures>,
    total: u64,
) -> Result<(), AcousticRealizeError> {
    if let Some(boundaries) = boundaries {
        boundaries.validate_step_budget(total)?;
    }
    Ok(())
}

pub(in crate::bernoulli_aperture) fn extend_step_budget(
    boundaries: Option<&mut BoundaryApertures>,
    total: u64,
) -> Result<(), AcousticRealizeError> {
    if let Some(boundaries) = boundaries {
        boundaries.extend_step_budget(total)?;
    }
    Ok(())
}
