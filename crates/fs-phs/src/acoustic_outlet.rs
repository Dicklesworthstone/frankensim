//! Physical outward volume flow from the original LC-chain coordinates.
//!
//! The native charges and fluxes alternate orientation along the chain. The
//! final storage effort is therefore multiplied by (-1)^N. Radiation and
//! viscothermal inertia already belong to that storage; an observer must not
//! replace it with rho*dx/A. This image explicitly selects the original
//! compression-volume/volume-flow frustum owner, including for a single cone.

use crate::{
    AcousticSection, AcousticTap, ChainOutlet, PhsError, PortHamiltonian, ViscothermalPin, WallPin,
    acoustic_frustum_chain_outlet_wall, acoustic_frustum_chain_outlet_wall_in_exterior,
};

/// Local numeric gas coefficients for one authored frustum section.
///
/// The caller supplies a coherent material state at the common static pressure.
/// Geometry and region boundaries remain those of [`AcousticSection`]; each
/// cell belongs wholly to its section, with no averaged interface density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticSectionMedium {
    /// Local equilibrium mass density [kg/m3].
    pub density: f64,
    /// Local adiabatic sound speed [m/s].
    pub sound_speed: f64,
    /// Same region's viscothermal coefficients and numerical memory budget.
    pub viscothermal: Option<ViscothermalPin>,
}

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

/// Assemble the original frustum chain and its physical terminal-flow row.
///
/// All geometry, material, boundary, loss and wall admission remains with
/// [`acoustic_frustum_chain_outlet_wall`]. The row uses its actual complete quadratic
/// storage, including terminal radiation mass and resolved viscothermal inertia.
/// A single widening cone uses that same mass-conserving owner rather than the
/// legacy spherical psi-line, whose compact loaded boundary omits the Euler
/// near-field contribution to physical flow. Callers must admit the returned
/// frustum model's spatial error and advance this same model when observing it.
/// The legacy source-chain API keeps its existing automatic model choice.
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
    let system = acoustic_frustum_chain_outlet_wall(
        sections,
        density,
        sound_speed,
        outlet,
        inlets,
        taps,
        viscothermal,
        wall,
    )?;
    with_outlet_flow(system, sections, outlet)
}

/// The same physical chain with an independently supplied exterior at its end.
///
/// Cell storage, viscothermal response, wall loading and the acoustic pin
/// frequency retain the interior `density` and `sound_speed`. Only the original
/// compact/Foster terminal radiation load uses the exterior density and speed.
/// Pressure and outward volume flow are continuous across this frozen boundary;
/// neither receives a density-ratio multiplier. The caller must supply coherent
/// media at a common static pressure: this numerical layer has no gas-state EOS.
/// Distinct media with side taps refuse because their exterior loads need their
/// own owners. Equal media preserve the original assembly arithmetic.
///
/// # Errors
/// Original chain refusals, nonpositive/nonfinite exterior properties, distinct
/// exterior with side taps, or an unrepresentable physical observation map.
#[allow(clippy::too_many_arguments)]
pub fn acoustic_chain_with_outlet_flow_in_exterior(
    sections: &[AcousticSection],
    density: f64,
    sound_speed: f64,
    exterior_density: f64,
    exterior_sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    taps: &[AcousticTap],
    viscothermal: Option<&ViscothermalPin>,
    wall: Option<&WallPin>,
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
    let system = acoustic_frustum_chain_outlet_wall_in_exterior(
        sections,
        density,
        sound_speed,
        exterior_density,
        exterior_sound_speed,
        outlet,
        inlets,
        taps,
        viscothermal,
        wall,
        None,
    )?;
    with_outlet_flow(system, sections, outlet)
}

/// Original frustum chain with one local medium per section and an exterior.
///
/// Each complete series cell retains its exact integral of rho/area, and each
/// compression cell its volume/(rho*c^2). The original incidence joins their
/// pressure/volume-flow ports without a density ratio or extra interface mass.
/// Local Bessel/Foster terms use that same cell's medium. Their shared pin is
/// the quarter-wave scale from the sum of section travel times. The terminal
/// radiation owner uses the exterior and the terminal cell's memory budget;
/// the authored numerical wall retains its original geometry and law.
///
/// Identical section media use the uniform assembler's original arithmetic.
/// Gas EOS, common static pressure, frozen-region validity and spatial-error
/// admission remain caller responsibilities. Side openings are not represented
/// by this regional entry point.
///
/// # Errors
/// Missing/misaligned or invalid regional coefficients, original chain
/// admission, or an unrepresentable physical flow map.
pub fn acoustic_chain_with_regional_outlet_flow(
    sections: &[AcousticSection],
    media: &[AcousticSectionMedium],
    exterior_density: f64,
    exterior_sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    wall: Option<&WallPin>,
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
    let Some(first) = media.first() else {
        return Err(PhsError::Unadmitted {
            what: "one acoustic medium per section required",
        });
    };
    if media.len() != sections.len() {
        return Err(PhsError::Unadmitted {
            what: "one acoustic medium per section required",
        });
    }
    let regional = (!media.iter().all(|medium| medium == first)).then_some(media);
    let system = acoustic_frustum_chain_outlet_wall_in_exterior(
        sections,
        first.density,
        first.sound_speed,
        exterior_density,
        exterior_sound_speed,
        outlet,
        inlets,
        &[],
        first.viscothermal.as_ref(),
        wall,
        regional,
    )?;
    with_outlet_flow(system, sections, outlet)
}

fn with_outlet_flow(
    system: PortHamiltonian,
    sections: &[AcousticSection],
    outlet: ChainOutlet,
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
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
