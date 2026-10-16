//! Admission of physical observations; this does not certify instrument fidelity.
use crate::pcm_wav::baffled::RayleighMedium;

/// Acceptance excludes internal-pressure/jet proxies and unmatched reaction loads.
/// Diagnostic retains those explicitly selected observations for comparisons.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObservationPolicy {
    Acceptance,
    Diagnostic,
}

/// Named cold-construction refusals, before either physical history advances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObservationError {
    LegacyListenerInAcceptance,
    UnmatchedReactionLoad,
    ReactionReceiverGasMismatch { field: &'static str },
    ReactionBinding { what: &'static str },
}

impl core::fmt::Display for ObservationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::LegacyListenerInAcceptance => write!(
                f,
                "legacy listener_m is diagnostic-only; acceptance requires a physical outlet receiver"
            ),
            Self::UnmatchedReactionLoad => write!(
                f,
                "acceptance requires the actual matched reaction-load binding"
            ),
            Self::ReactionReceiverGasMismatch { field } => write!(
                f,
                "reaction/receiver gas mismatch: {field} differs by more than 1e-12 relative"
            ),
            Self::ReactionBinding { what } => {
                write!(f, "reaction/receiver binding mismatch: {what}")
            }
        }
    }
}
impl std::error::Error for ObservationError {}

/// Compare the actual load and receiver media, never the internal propagation gas.
/// Invalid values refuse too; an unknown medium cannot establish a match.
pub fn validate_reaction_medium(
    reaction: RayleighMedium,
    receiver: RayleighMedium,
) -> Result<(), ObservationError> {
    for (field, a, b) in [
        ("density_kg_m3", reaction.density, receiver.density),
        (
            "sound_speed_m_s",
            reaction.sound_speed,
            receiver.sound_speed,
        ),
    ] {
        if !a.is_finite()
            || !b.is_finite()
            || a <= 0.0
            || b <= 0.0
            || (a - b).abs() > 1e-12 * a.abs().max(b.abs())
        {
            return Err(ObservationError::ReactionReceiverGasMismatch { field });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reaction_receiver_gas_mismatch_refuses() {
        let original = RayleighMedium {
            density: 1.0,
            sound_speed: 1.0,
        };
        assert_eq!(validate_reaction_medium(original, original), Ok(()));
        for field in ["density_kg_m3", "sound_speed_m_s"] {
            // Exactly representable values bracket the 1e-12 relative band.
            let cap = 1.0 + 4503.0 * f64::EPSILON;
            let over = 1.0 + 4504.0 * f64::EPSILON;
            let changed = |value| {
                if field == "density_kg_m3" {
                    RayleighMedium {
                        density: value,
                        ..original
                    }
                } else {
                    RayleighMedium {
                        sound_speed: value,
                        ..original
                    }
                }
            };
            assert_eq!(validate_reaction_medium(original, changed(cap)), Ok(()));
            assert_eq!(
                validate_reaction_medium(original, changed(over)),
                Err(ObservationError::ReactionReceiverGasMismatch { field })
            );
            for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert!(validate_reaction_medium(original, changed(invalid)).is_err());
            }
        }
    }
}
