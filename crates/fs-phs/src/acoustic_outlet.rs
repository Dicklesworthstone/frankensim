//! Physical outward volume flow from the original LC-chain coordinates.
//!
//! The native charges and fluxes alternate orientation along the chain. The
//! final storage effort is therefore multiplied by (-1)^N. Radiation and
//! viscothermal inertia already belong to that storage; an observer must not
//! replace it with rho*dx/A. The spherical psi-line refuses because its loaded
//! terminal omits an Euler near-field flow needed by a physical microphone.

use crate::{
    AcousticSection, AcousticTap, ChainOutlet, PhsError, PortHamiltonian, ViscothermalPin, WallPin,
    acoustic_chain_outlet_wall, spherical_chain_section,
};

/// Immutable physical flow row belonging to one returned acoustic chain.
///
/// Positive flow leaves the terminal, in m^3/s. A rigid terminal has an exactly
/// zero row, including when its disconnected retained flux coordinate is nonzero.
/// This observes the loaded terminal only; side openings need their own physical
/// receivers. It neither adds a port nor changes the system's power balance.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticOutletFlow {
    coefficients: Vec<f64>,
}

impl AcousticOutletFlow {
    /// Original-state linear functional; retained for composition and harmonic
    /// response without inferring an acoustic state layout in another crate.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Observe one accepted native state without allocation or state changes.
    ///
    /// # Errors
    /// A different state extent, nonfinite coordinate or unrepresentable flow.
    pub fn volume_flow_m3_s(&self, state: &[f64]) -> Result<f64, PhsError> {
        if state.len() != self.coefficients.len() || state.iter().any(|x| !x.is_finite()) {
            return Err(PhsError::Dimension {
                what: "outlet observation requires the complete finite original chain state",
            });
        }
        let flow = self
            .coefficients
            .iter()
            .zip(state)
            .map(|(a, x)| a * x)
            .sum::<f64>();
        if !flow.is_finite() {
            return Err(PhsError::Unadmitted {
                what: "physical outlet flow is unrepresentable",
            });
        }
        Ok(flow)
    }
}

/// Assemble the unchanged acoustic chain and its physical terminal-flow row.
///
/// All geometry, material, boundary, loss and wall admission remains with
/// [`acoustic_chain_outlet_wall`]. The row uses its actual complete quadratic
/// storage, including terminal radiation mass and resolved viscothermal inertia.
/// Cylinder/frustum chains retain their original coordinate conventions. The
/// automatically selected single widening spherical cone refuses: its compact
/// loaded terminal neglects the Euler near-field term in physical volume flow,
/// so scaled terminal storage effort alone is not a physical outlet history.
/// The original cone remains available through the source-chain API. Spatial
/// and constitutive approximations remain those of the original admitted chain;
/// this is not a new acoustic solver.
///
/// # Errors
/// Original chain refusals or an unrepresentable physical observation map.
#[allow(clippy::too_many_arguments)]
pub fn acoustic_chain_with_outlet_flow(
    sections: &[AcousticSection],
    density: f64,
    sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    taps: &[AcousticTap],
    viscothermal: Option<&ViscothermalPin>,
    wall: Option<&WallPin>,
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
    if spherical_chain_section(sections, outlet).is_some() {
        return Err(PhsError::Unadmitted {
            what: "physical spherical outlet requires the unresolved Euler near-field flow",
        });
    }
    let system = acoustic_chain_outlet_wall(
        sections,
        density,
        sound_speed,
        outlet,
        inlets,
        taps,
        viscothermal,
        wall,
    )?;
    let n = system.state_dim();
    let mut coefficients = vec![0.0; n];
    if outlet != ChainOutlet::Rigid {
        let cells: usize = sections.iter().map(|section| section.cells).sum();
        let orientation = if cells % 2 == 0 { 1.0 } else { -1.0 };
        let q = system
            .storage
            .exact_quadratic_matrix()
            .ok_or(PhsError::Unadmitted {
                what: "native acoustic outlet requires its declared quadratic storage",
            })?;
        let terminal = 2 * cells - 1;
        for (coefficient, value) in coefficients
            .iter_mut()
            .zip(&q[terminal * n..(terminal + 1) * n])
        {
            *coefficient = orientation * value;
        }
        if coefficients.iter().any(|value| !value.is_finite()) || coefficients[terminal] == 0.0 {
            return Err(PhsError::Unadmitted {
                what: "native acoustic outlet flow map is unrepresentable",
            });
        }
    }
    Ok((system, AcousticOutletFlow { coefficients }))
}
