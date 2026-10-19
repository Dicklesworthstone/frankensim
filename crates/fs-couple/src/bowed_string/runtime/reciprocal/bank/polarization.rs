//! Two transverse planes of ONE circular material string.
//!
//! Both planes solve mu*y_tt = T*y_xx - EI*y_xxxx with the same material/gas
//! losses. For a unit transverse direction d, a physical port has row
//! [d_x*phi, d_y*phi]; its transpose applies force to those same coordinates.
//! Thus rotation changes neither mass nor work. A bridge spring stores
//! k*(d dot y - y_body)^2/2, so body compliance couples the planes and shifts
//! their assembled resonances without an authored detuning or output gain.
//!
//! See Mansour, Woodhouse & Scavone (2016), "Enhanced wave-based modelling of
//! musical strings. Part 2: Bowed strings", for two-polarization/body coupling:
//! https://api.repository.cam.ac.uk/server/api/core/bitstreams/bc2b1bfc-8b61-4499-999e-3de6dd494927/content
//!
//! The section is circular and frozen-state, with identical flexural laws in
//! both planes. Translational bridge ports follow their supplied direction;
//! optional torque/rotation connections use the original material torsion in
//! the same joint solve. The explicit `MaterialStringSource::into_nonlinear`
//! path retains both planes in one Kirchhoff--Carrier elongation/tension channel;
//! the linear bank continues to use its original quadratic storage. A six-DOF
//! termination and a resolved bow-normal contact law remain outside this image.
//! Optional twist carries no fictitious direct acoustic output.

use super::*;
use fs_material::state_point::MaterialStatePointError;

/// Physical directions expressed in one common orthonormal transverse frame.
/// No modal weights, gains, frequency offsets or implicit normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseDirections {
    bow: [f64; 2],
    bridge: [f64; 2],
}
impl TransverseDirections {
    /// Admit the tangent of the bow stroke and the string/body connector axis.
    /// Each vector must have unit squared norm within 1e-12. This tolerance is
    /// only floating-point direction admission; supplied vectors are retained.
    pub fn try_new(bow: [f64; 2], bridge: [f64; 2]) -> Result<Self, MaterialStatePointError> {
        for direction in [bow, bridge] {
            let norm2 = direction[0] * direction[0] + direction[1] * direction[1];
            if direction.iter().any(|x| !x.is_finite())
                || !norm2.is_finite()
                || (norm2 - 1.0).abs() > 1e-12
            {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "transverse bow and bridge directions must be finite unit vectors",
                });
            }
        }
        Ok(Self { bow, bridge })
    }
    /// Bow tangent; centreline load ports use this same direction without torque.
    #[must_use]
    pub const fn bow(&self) -> [f64; 2] {
        self.bow
    }
    /// Direction of the physical scalar connection to the common receiver.
    #[must_use]
    pub const fn bridge(&self) -> [f64; 2] {
        self.bridge
    }
}

/// Original circular section, optional original torsion and physical directions.
/// Coordinate order is first transverse plane, second plane, then any twist.
#[derive(Clone, Debug)]
pub struct PolarizedStringSource {
    source: MaterialStringSource,
    directions: TransverseDirections,
    total_modes: usize,
}
impl PolarizedStringSource {
    /// Original source; no repeated material query or copied material identity.
    #[must_use]
    pub const fn source(&self) -> &MaterialStringSource {
        &self.source
    }
    /// Original directions, independent of the arbitrary transverse basis choice.
    #[must_use]
    pub const fn directions(&self) -> TransverseDirections {
        self.directions
    }
    pub(super) const fn total_modes(&self) -> usize {
        self.total_modes
    }

    pub(super) fn add_torsion(
        self,
        shear: &str,
        viscosity: &str,
        modes: usize,
    ) -> Result<MaterialStringSource, MaterialStatePointError> {
        self.source
            .with_torsion(shear, viscosity, modes)?
            .with_transverse_polarizations(self.directions)
    }
}
impl MaterialStringSource {
    /// Retain both transverse displacement components from the SAME source.
    /// No second material, mass, damping curve or artificial pitch split exists.
    /// The body/friction network acts on the declared directions and determines
    /// its own polarization mixing. Adding torsion before or after this operation
    /// gives the same coordinate order and physical source.
    ///
    /// # Errors
    /// Already polarized/nested sources, zero or oversized transverse basis,
    /// or mode-count overflow. Full network and material-band limits still apply
    /// when the source is lowered into the bank.
    pub fn with_transverse_polarizations(
        self,
        directions: TransverseDirections,
    ) -> Result<Self, MaterialStatePointError> {
        let bending = match &self {
            Self::WithTorsion(t) => t.bending(),
            other => other,
        };
        let count = match bending {
            Self::Homogeneous(s) => s.string().n_modes,
            Self::Concentric { template, .. } => template.n_modes,
            Self::WithPolarizations(_) | Self::WithTorsion(_) => {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "two transverse planes require one original unpolarized source",
                });
            }
        };
        if count == 0 || count > 4096 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "each transverse plane requires 1..=4096 modes",
            });
        }
        let total_modes =
            self.modes()
                .checked_add(count)
                .ok_or(MaterialStatePointError::InvalidDerived {
                    quantity: "combined transverse and torsional dimension overflow",
                })?;
        Ok(Self::WithPolarizations(Box::new(PolarizedStringSource {
            source: self,
            directions,
            total_modes,
        })))
    }

    pub(super) fn torsion_source(&self) -> Option<&torsion::TorsionalStringSource> {
        match self {
            Self::WithPolarizations(p) => p.source.torsion_source(),
            Self::WithTorsion(t) => Some(t),
            Self::Homogeneous(_) | Self::Concentric { .. } => None,
        }
    }

    pub(super) fn append_extra_modes(
        &self,
        card: &BowedStringCard,
        modes: &mut Vec<ModalAcousticMode>,
        gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        if matches!(self, Self::WithPolarizations(_)) {
            // An independent displacement component has the SAME physical mass,
            // stiffness and loss; this does not double the specimen's line mass.
            for k in 0..card.mode_count {
                poll(gate)?;
                modes.push(modes[k]);
            }
        }
        if let Some(twist) = self.torsion_source() {
            twist.append_modes(card, modes, gate)?;
        }
        Ok(())
    }

    pub(super) fn transverse_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        self.projected_shape(card, station, false)
    }
    pub(super) fn bridge_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        self.projected_shape(card, station, true)
    }
    fn projected_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
        bridge: bool,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        let mut row = shape(card, station);
        if let Self::WithPolarizations(p) = self {
            let direction = if bridge {
                p.directions.bridge
            } else {
                p.directions.bow
            };
            for k in 0..card.mode_count {
                row.push(direction[1] * row[k]);
            }
            for value in &mut row[..card.mode_count] {
                *value *= direction[0];
            }
        }
        row.resize(self.modes(), 0.0);
        if row.iter().any(|x| !x.is_finite()) || row.iter().all(|x| *x == 0.0) {
            return Err(invalid(
                "transverse string direction has no finite nonzero physical port",
            ));
        }
        Ok(row)
    }
}
