//! Spatial viscous attachments in the unchanged mechanical coordinates.
//!
//! A port has physical velocity b.v and reaction -c*b*(b.v). Thus the
//! momentum resistance is c*b*b^T and its power loss is c*(b.v)^2 >= 0.
//! Keeping cross terms is essential: damping a point is NOT independent
//! damping of each participating mode. A mode with a node at the attachment
//! is unaffected. Opposite-signed body weights give reciprocal relative drag.
//!
//! These are fixed, bilateral, ideal dashpots, optionally grounded. No finger
//! mass, preload, friction, contact detection or measured material is implied.
//! They change mechanics before radiation, never PCM gain or modal frequency.
use super::{ImpactError, invalid};

/// Bounded number of localized dissipative ports in the impact description.
/// Each execution image also applies its existing matrix/connection budgets.
pub const MAX_VISCOUS_DAMPERS: usize = 16;

/// One physical viscous port, e.g. an idealized held muffler on a head or shell.
#[derive(Debug, Clone)]
pub struct ViscousDamper {
    /// Signed physical displacement participation [1/sqrt(kg)], in complete
    /// body/mode order. Interpolate the actual mode shapes at the attachment;
    /// do not normalize this row. Unattached coordinates have exact zeros.
    pub weights: Vec<f64>,
    /// Nonnegative physical resistance [N s/m]. Zero disables the port without
    /// introducing a numerical connection or altering the undamped trajectory.
    pub damping_n_s_m: f64,
}

pub(crate) fn validate(dampers: &[ViscousDamper], modes: usize) -> Result<(), ImpactError> {
    if dampers.len() > MAX_VISCOUS_DAMPERS {
        return Err(invalid("localized damping exceeds the physical port budget"));
    }
    for damper in dampers {
        if damper.weights.len() != modes || modes == 0
            || damper.weights.iter().any(|b| !b.is_finite())
            || damper.weights.iter().all(|b| *b == 0.0)
            || !damper.damping_n_s_m.is_finite() || damper.damping_n_s_m < 0.0 {
            return Err(invalid("viscous damper needs a finite moving port and nonnegative SI resistance"));
        }
        if damper.damping_n_s_m == 0.0 { continue; }
        let root = damper.damping_n_s_m.sqrt();
        for &weight in &damper.weights {
            let scaled = root*weight;
            let diagonal = scaled*scaled;
            if !diagonal.is_finite() || (weight != 0.0 && diagonal == 0.0) {
                return Err(invalid("viscous damper resistance is not representable in this basis"));
            }
        }
    }
    Ok(())
}

/// Add the symmetric momentum block only; Kelvin and position rows stay intact.
/// Called during construction, never while a step has partially advanced.
pub(super) fn add_resistance(dampers: &[ViscousDamper], modes: usize,
    dimension: usize, resistance: &mut [f64]) -> Result<(), ImpactError>
{
    validate(dampers, modes)?;
    for damper in dampers {
        if damper.damping_n_s_m == 0.0 { continue; }
        let root = damper.damping_n_s_m.sqrt();
        for i in 0..modes {
            for j in 0..=i {
                let index = (2*i+1)*dimension+2*j+1;
                let value = resistance[index]+(root*damper.weights[i])*(root*damper.weights[j]);
                if !value.is_finite() {
                    return Err(invalid("combined viscous resistance overflows the mechanical basis"));
                }
                resistance[index] = value;
                resistance[(2*j+1)*dimension+2*i+1] = value;
            }
        }
    }
    Ok(())
}

/// Preserve solid ports when appending acoustic/neck inertia. A muffler cannot
/// act directly on gas coordinates merely because they share a state vector.
pub(crate) fn extend(mut dampers: Vec<ViscousDamper>, structural: usize,
    total: usize) -> Result<Vec<ViscousDamper>, ImpactError>
{
    validate(&dampers, structural)?;
    if total < structural { return Err(invalid("damper extension cannot discard structural coordinates")); }
    for damper in &mut dampers { damper.weights.resize(total, 0.0); }
    Ok(dampers)
}

#[cfg(test)]
mod tests;
