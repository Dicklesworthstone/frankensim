//! Whole-trajectory resistance multipliers from the original contact duals.
//!
//! One scalar multiplies the complete source-backed law of each named
//! interface at every endpoint. Native temperature/phase feedback stays in
//! the transpose; its endpoint contractions already include temporal weights.

use super::{ScanSensitivityError, ScanningEnthalpyRun, finite, invalid, poll};
use crate::interface::ContactResistanceGradient;
use fs_exec::Cx;

pub(super) fn control_count(
    cx: &Cx<'_>,
    run: &ScanningEnthalpyRun<'_, '_, '_>,
) -> Result<usize, ScanSensitivityError> {
    poll(cx, 0)?;
    // Count named controls once, including fixed interfaces in the live model
    // and nonmatching interfaces with delegated exact-matching pieces.
    Ok(run
        .contacts()
        .map_or(0, |contacts| contacts.surface_count()))
}

pub(super) fn accumulate(
    cx: &Cx<'_>,
    run: &ScanningEnthalpyRun<'_, '_, '_>,
    total: &mut Vec<ContactResistanceGradient>,
    values: Vec<ContactResistanceGradient>,
) -> Result<(), ScanSensitivityError> {
    if values.len() != control_count(cx, run)? {
        return Err(invalid(
            "scan contact-gradient layout changed from its complete native model",
        ));
    }
    let mut previous: Option<&str> = None;
    for (i, value) in values.iter().enumerate() {
        poll(cx, i)?;
        if previous.is_some_and(|name| name >= value.interface.as_str()) {
            return Err(invalid(
                "scan contact controls lost their native canonical order",
            ));
        }
        finite(value.derivative)?;
        previous = Some(&value.interface);
    }
    if total.is_empty() {
        *total = values;
        return Ok(());
    }
    if total.len() != values.len() {
        return Err(invalid("scan contact layout changed during reverse sweep"));
    }
    for (i, (a, b)) in total.iter_mut().zip(values).enumerate() {
        poll(cx, i)?;
        if a.interface != b.interface
            || a.card_identity != b.card_identity
            || a.mapped != b.mapped
            || a.face_pairs != b.face_pairs
        {
            return Err(invalid(
                "scan contact source or trace identity changed during reverse sweep",
            ));
        }
        // This is dJ/ds at s=1 for the WHOLE resistance law, not dJ/dR at a
        // single frozen face. Do not multiply by dt or solve another transpose.
        a.derivative = finite(a.derivative + b.derivative)?;
    }
    Ok(())
}
