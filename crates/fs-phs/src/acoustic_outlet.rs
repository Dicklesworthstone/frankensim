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

/// Geometry of one original frustum cell, in native cell order.
/// Material adapters use this query instead of reconstructing a second mesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticCellGeometry {
    /// Index of the authored section that owns this complete cell.
    pub section: usize,
    /// Original midpoint radius [m].
    pub radius: f64,
    /// Axial integration length [m] used by gas storage and series inertia.
    /// Surface slant never multiplies this quantity.
    pub axial_length: f64,
    /// Wetted strip length [m], including the frustum's surface slant.
    /// Its area is `2*pi*radius*wall_strip_length`.
    pub wall_strip_length: f64,
    pub(crate) x0: f64,
    pub(crate) inertance_area: f64,
    pub(crate) compression_area: f64,
}

impl AcousticCellGeometry {
    /// Actual wetted cell area [m2], distinct from its axial cross section.
    #[must_use]
    pub fn wall_area_m2(&self) -> f64 {
        2.0 * core::f64::consts::PI * self.radius * self.wall_strip_length
    }
}

/// One original Maxwell arm in local surface-displacement units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticWallRelaxation {
    /// Positive relaxing stiffness [Pa/m], excluding equilibrium stiffness.
    pub stiffness_per_area: f64,
    /// Positive relaxation time [s].
    pub relaxation_time_s: f64,
}

/// Complete local wall law for one cell returned by [`acoustic_cell_geometry`].
/// Material source and frequency validity remain with the caller.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticCellWall {
    /// Original surface mass, equilibrium stiffness and Kelvin resistance.
    pub wall: WallPin,
    /// Every positive original hereditary arm, without a frequency fit.
    pub relaxation: Vec<AcousticWallRelaxation>,
}

/// Query the same cell geometry used by the original frustum assembler.
///
/// # Errors
/// Empty/nonphysical sections or an unrepresentable cell count/geometry.
pub fn acoustic_cell_geometry(
    sections: &[AcousticSection],
) -> Result<Vec<AcousticCellGeometry>, PhsError> {
    let mut count = 0usize;
    for section in sections {
        if ![section.length, section.radius, section.outlet()]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            || section.cells == 0
        {
            return Err(PhsError::Unadmitted {
                what: "acoustic section parameters",
            });
        }
        count = count
            .checked_add(section.cells)
            .ok_or(PhsError::Unadmitted {
                what: "acoustic cell count overflow",
            })?;
    }
    if count == 0 {
        return Err(PhsError::Unadmitted {
            what: "acoustic chain parameters",
        });
    }
    let mut cells = Vec::with_capacity(count);
    let mut x_acc = 0.0;
    for (section, s) in sections.iter().enumerate() {
        let dx = s.length / s.cells as f64;
        let r_in = s.radius;
        let r_out = s.outlet();
        let n_s = s.cells as f64;
        for i in 0..s.cells {
            let ra = r_in + (r_out - r_in) * (i as f64 / n_s);
            let rb = r_in + (r_out - r_in) * ((i + 1) as f64 / n_s);
            let slope = (rb - ra) / dx;
            cells.push(AcousticCellGeometry {
                section,
                radius: f64::midpoint(ra, rb),
                wall_strip_length: (1.0 + slope * slope).sqrt(),
                x0: x_acc,
                axial_length: dx,
                inertance_area: core::f64::consts::PI * ra * rb,
                compression_area: core::f64::consts::PI * (ra * ra + ra * rb + rb * rb) / 3.0,
            });
            x_acc += dx;
        }
    }
    for index in 0..cells.len() {
        let end = cells.get(index + 1).map_or(x_acc, |cell| cell.x0);
        let strip = cells[index].wall_strip_length * (end - cells[index].x0);
        cells[index].wall_strip_length = strip;
        if !cells[index].wall_strip_length.is_finite() || cells[index].wall_strip_length <= 0.0 {
            return Err(PhsError::Unadmitted {
                what: "acoustic wall strip geometry",
            });
        }
    }
    Ok(cells)
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

/// Outward flow from an active original side opening, in native tap order.
///
/// A compact opening uses its single neck flow; a distributed chimney uses
/// the terminal flux after its final compression cell. Closed pads have no
/// emitting row. The effective area already includes the authored open fraction:
/// no additional square-root or area multiplier belongs on `flow`.
#[derive(Clone, Debug, PartialEq)]
pub struct AcousticTapFlow {
    /// Index in the original `AcousticTap` slice, including closed entries.
    pub tap_index: usize,
    /// Authored fraction of the chain's total axial length.
    pub station: f64,
    /// Radius of the original emitting aperture: neck radius * sqrt(open fraction).
    pub effective_radius_m: f64,
    /// Complete-state outward volume flow [m3/s], with native orientation removed.
    pub flow: AcousticOutletFlow,
}

// Filled only where the original assembler selects a compact or long branch.
// Coupling code never reconstructs state offsets or estimates terminal inertia.
pub(crate) struct AcousticTapTerminal {
    pub tap_index: usize,
    pub station: f64,
    pub effective_radius_m: f64,
    pub coordinate: usize,
    pub orientation: f64,
}

/// Assemble one original uniform frustum chain and every radiating mouth's flow.
///
/// The returned terminal row and active side rows belong to the same complete
/// system, including original viscothermal, pad, wall and radiation memories.
/// Source geometry, aperture corrections and radiation loading are unchanged;
/// these immutable observations add no resistance or power port. Each side
/// opening emits into the same gas used by its native neck/chimney owner.
/// Separate exterior media and regional gas/wall partitions need their own
/// complete tap ownership and are not introduced by this entry point.
///
/// # Errors
/// Original chain admission or an unrepresentable complete-state flow row.
#[allow(clippy::too_many_arguments)]
pub fn acoustic_chain_with_opening_flows(
    sections: &[AcousticSection],
    density: f64,
    sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    taps: &[AcousticTap],
    viscothermal: Option<&ViscothermalPin>,
    wall: Option<&WallPin>,
) -> Result<(PortHamiltonian, AcousticOutletFlow, Vec<AcousticTapFlow>), PhsError> {
    let mut terminals = Vec::new();
    let system = acoustic_frustum_chain_outlet_wall_in_exterior(
        sections,
        density,
        sound_speed,
        density,
        sound_speed,
        outlet,
        inlets,
        taps,
        viscothermal,
        wall,
        None,
        None,
        Some(&mut terminals),
    )?;
    let (system, bell) = with_outlet_flow(system, sections, outlet)?;
    let q = declared_quadratic(&system)?;
    let openings = terminals
        .into_iter()
        .map(|terminal| {
            Ok(AcousticTapFlow {
                tap_index: terminal.tap_index,
                station: terminal.station,
                effective_radius_m: terminal.effective_radius_m,
                flow: effort_flow_row(
                    &q,
                    system.state_dim(),
                    terminal.coordinate,
                    terminal.orientation,
                )?,
            })
        })
        .collect::<Result<Vec<_>, PhsError>>()?;
    drop(q);
    Ok((system, bell, openings))
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
        None,
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
    regional_outlet(
        sections,
        media,
        exterior_density,
        exterior_sound_speed,
        outlet,
        inlets,
        wall,
        None,
    )
}

/// Bind complete cell walls to the same regional gas/frustum owner.
///
/// Wall pressure and volume flow share the original cell junction. Equilibrium
/// storage and Kelvin resistance use `apply_cell_wall`; hereditary arms use
/// the original `PortHamiltonian::with_relaxation_branches`. The outlet row is
/// derived only after that complete system exists. Source walls use the actual
/// slanted strip area; the existing uniform numerical-wall API is unchanged.
/// No axial shell coupling, pressure prestress or exterior wall radiation is
/// inferred. The caller retains material-source and excited-band admission.
///
/// # Errors
/// Missing/misaligned cell walls, invalid base or hereditary coefficients,
/// original gas/geometry admission, or an unrepresentable flow map.
#[allow(clippy::too_many_arguments)]
pub fn acoustic_chain_with_regional_outlet_flow_and_cell_walls(
    sections: &[AcousticSection],
    media: &[AcousticSectionMedium],
    exterior_density: f64,
    exterior_sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    walls: &[AcousticCellWall],
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
    regional_outlet(
        sections,
        media,
        exterior_density,
        exterior_sound_speed,
        outlet,
        inlets,
        None,
        Some(walls),
    )
}

#[allow(clippy::too_many_arguments)]
fn regional_outlet(
    sections: &[AcousticSection],
    media: &[AcousticSectionMedium],
    exterior_density: f64,
    exterior_sound_speed: f64,
    outlet: ChainOutlet,
    inlets: usize,
    wall: Option<&WallPin>,
    cell_walls: Option<&[AcousticCellWall]>,
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
        cell_walls,
        None,
    )?;
    with_outlet_flow(system, sections, outlet)
}

fn declared_quadratic(system: &PortHamiltonian) -> Result<std::borrow::Cow<'_, [f64]>, PhsError> {
    let q = match system.storage.exact_quadratic_matrix() {
        Some(q) => std::borrow::Cow::Borrowed(q),
        None => {
            // Original Maxwell storage supplies an owned exact quadratic
            // image. Derive immutable rows without gradient probes or a second
            // state model; retain the borrowed path for native quadratic storage.
            let part = system
                .storage
                .quadratic_part()
                .filter(|part| part.exact)
                .ok_or(PhsError::Unadmitted {
                    what: "native acoustic outlet requires its declared quadratic storage",
                })?;
            std::borrow::Cow::Owned(part.matrix)
        }
    };
    if Some(q.len()) != system.state_dim().checked_mul(system.state_dim()) {
        return Err(PhsError::Dimension {
            what: "native acoustic outlet quadratic storage",
        });
    }
    Ok(q)
}

fn effort_flow_row(
    q: &[f64],
    n: usize,
    terminal: usize,
    orientation: f64,
) -> Result<AcousticOutletFlow, PhsError> {
    if terminal >= n {
        return Err(PhsError::Dimension {
            what: "native acoustic outlet coordinate",
        });
    }
    let coefficients: Vec<_> = q[terminal * n..(terminal + 1) * n]
        .iter()
        .map(|value| orientation * value)
        .collect();
    if coefficients.iter().any(|value| !value.is_finite()) || coefficients[terminal] == 0.0 {
        return Err(PhsError::Unadmitted {
            what: "native acoustic outlet flow map is unrepresentable",
        });
    }
    Ok(AcousticOutletFlow { coefficients })
}

fn with_outlet_flow(
    system: PortHamiltonian,
    sections: &[AcousticSection],
    outlet: ChainOutlet,
) -> Result<(PortHamiltonian, AcousticOutletFlow), PhsError> {
    let flow = if outlet == ChainOutlet::Rigid {
        AcousticOutletFlow {
            coefficients: vec![0.0; system.state_dim()],
        }
    } else {
        let cells: usize = sections.iter().map(|section| section.cells).sum();
        let orientation = if cells % 2 == 0 { 1.0 } else { -1.0 };
        effort_flow_row(
            &declared_quadratic(&system)?,
            system.state_dim(),
            2 * cells - 1,
            orientation,
        )?
    };
    Ok((system, flow))
}
