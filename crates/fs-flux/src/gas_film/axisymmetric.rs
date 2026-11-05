//! Radial geometry and total-unit output for the shared Reynolds solver.
//!
//! The no-slip continuum equation is d(ph)/dt =
//! (1/r) d[r*h^3*p*dp/dr/(12*mu)]/dr. The existing Maxwell wall law changes
//! its conductance without changing annular storage or circular face measures.

use std::f64::consts::PI;

use super::{
    FilmFaceMetric, FilmMetrics, GasFilmBoundaryTopology, GasFilmCheckpoint, GasFilmError,
    GasFilmInput, GasFilmInputAuthority, GasFilmUncertainty, derived, nonnegative, positive,
    solve_film,
};

/// Parallel no-slip capillaries distributed over one complete annular cell.
///
/// The connection is azimuthally homogenized: it resolves the annular pressure,
/// not individual holes or junction flow. Tube storage, inertia and end effects
/// are excluded under the same envelope as [`super::CapillaryVent`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricCapillaryVent {
    /// Internal radius of each identical capillary [m].
    pub radius_m: f64,
    /// Fully developed flow length of each capillary [m].
    pub length_m: f64,
    /// Total number of identical tubes connected to the annular cell; nonzero.
    /// This is neither a count per circumference nor a count per radial width.
    pub capillary_count: u32,
}

/// Finite rigid volume containing the same isothermal ideal gas as the film.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricGasFilmReservoir {
    /// Total reservoir volume [m^3], independent of the film's radius or mesh.
    pub volume_m3: f64,
    /// Initial absolute pressure [Pa], used only without a checkpoint.
    pub initial_absolute_pressure_pa: f64,
}

/// Total capillary work and single-tube applicability at the accepted endpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricCapillaryVentReport {
    /// Total pressure-gradient work/viscous dissipation [W], `M*R*T*ln(p_f/p_r)`.
    /// Separate from film heat and from transported enthalpy.
    pub viscous_dissipation_w: f64,
    /// Diameter-based Reynolds number in one tube, not the complete bundle.
    pub reynolds_number: f64,
    /// Maximum radius-based Knudsen number at either tube endpoint.
    pub maximum_knudsen_number: f64,
    /// Maximum centreline Mach number in one tube.
    pub maximum_gas_mach_number: f64,
}

/// Actual finite gas inventory and its exchange with the complete radial film.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricGasFilmReservoirReport {
    /// Absolute pressure from the accepted mass, fixed temperature and volume [Pa].
    pub absolute_pressure_pa: f64,
    /// Total reservoir gas mass [kg].
    pub gas_mass_kg: f64,
    /// Represented total mass-storage rate [kg s^-1].
    pub storage_rate_kg_per_s: f64,
    /// Signed enthalpy entering the reservoir from the film [W].
    pub inward_enthalpy_power_w: f64,
    /// Heat required to maintain fixed temperature [W], `-R*T*M`.
    /// Filling rejects heat; this is not a solved heat-transfer law.
    pub required_isothermal_heat_to_gas_w: f64,
    /// Film plus reservoir storage [kg s^-1]; their exchange is internal.
    pub combined_mass_closure_residual_kg_per_s: f64,
}

/// Conservative totals for a complete axisymmetric film, not a strip per width.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisymmetricGasFilmReceipt {
    /// Shared Reynolds model identity; checkpoint identity also binds radial geometry.
    pub model_id: String,
    /// Caller input authority retained without promotion.
    pub input_authority: GasFilmInputAuthority,
    /// Caller-declared uncertainty, without propagation or certification.
    pub input_uncertainty: GasFilmUncertainty,
    /// Deterministic Newton corrections used; zero at initial equilibrium.
    pub iterations: u32,
    /// Maximum local mass-balance residual per annular area [kg m^-2 s^-1].
    pub max_mass_residual_kg_m2_s: f64,
    /// Gas mass integrated over the complete active disk/annulus [kg].
    pub gas_mass_kg: f64,
    /// Total gas storage rate [kg s^-1].
    pub storage_rate_kg_per_s: f64,
    /// Outward mass rate across the inner circular boundary [kg s^-1].
    /// Exactly zero at a disk's sealed symmetry axis.
    pub inner_boundary_outward_mass_flow_kg_per_s: f64,
    /// Outward mass rate across the outer circular boundary [kg s^-1].
    pub outer_boundary_outward_mass_flow_kg_per_s: f64,
    /// Total outward exchange through the annular vent/capillary connection [kg s^-1].
    pub vent_outward_mass_flow_kg_per_s: f64,
    /// Signed inner-boundary enthalpy transport using caller-declared enthalpy [W].
    pub inner_boundary_outward_enthalpy_power_w: f64,
    /// Signed outer-boundary enthalpy transport using caller-declared enthalpy [W].
    pub outer_boundary_outward_enthalpy_power_w: f64,
    /// Signed annular vent/capillary enthalpy transport [W].
    pub vent_outward_enthalpy_power_w: f64,
    /// Finite capillary mechanics, when present, outside the film's dissipation.
    pub capillary_vent: Option<AxisymmetricCapillaryVentReport>,
    /// Total finite reservoir inventory, when present. Both film ends are sealed.
    pub finite_reservoir: Option<AxisymmetricGasFilmReservoirReport>,
    /// Total storage plus outward boundary/vent mass rate [kg s^-1].
    pub mass_closure_residual_kg_per_s: f64,
    /// Integral of absolute pressure over annular cell areas [N].
    pub absolute_pressure_load_n: f64,
    /// Integral of pressure minus the declared gauge reference [N].
    pub gauge_pressure_load_n: f64,
    /// Local gas-on-wall radial shear at each cell [Pa], positive radially outward.
    /// Equal on both walls since in-plane wall motion is excluded. Wall-on-gas
    /// traction is its negative; this is a scalar radial field, not a net
    /// Cartesian force (the azimuthally integrated radial force cancels).
    /// Normal gas-on-upper-wall traction is the step's absolute pressure field.
    pub gas_on_wall_radial_traction_pa: Vec<Option<f64>>,
    /// Normal squeeze work into the gas, -integral(p*dh/dt*dA) [W].
    /// This is all wall work for the admitted stationary in-plane walls.
    pub normal_gap_power_to_gas_w: f64,
    /// Discrete pressure-gradient force work [W], separate from enthalpy transport.
    pub pressure_gradient_power_to_gas_w: f64,
    /// Bulk plus Maxwell-interface dissipation [W]; no temperature is evolved.
    pub viscous_heat_w: f64,
    /// Nonnegative interface part of `viscous_heat_w` [W].
    pub wall_slip_dissipation_w: f64,
    /// Largest admitted pressure-scaled film/bath Knudsen number.
    pub maximum_knudsen_number: f64,
    /// Largest radial gas Mach number in reconstructed face velocity profiles.
    pub maximum_gas_mach_number: f64,
    /// Largest adjacent-cell |dh/dr|.
    pub maximum_gap_slope: f64,
}

/// Accepted axisymmetric pressure field, total-unit receipt and restart state.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisymmetricGasFilmStep {
    /// Cell-centered absolute pressure [Pa]; excluded outer cells have no gas state.
    pub absolute_pressure_pa: Vec<Option<f64>>,
    /// Annular storage, circular mass rates, total load/work, and local radial shear.
    pub receipt: AxisymmetricGasFilmReceipt,
    /// Shared pressure/gap checkpoint bound to the radial geometry. A finite
    /// reservoir uses `capillary_reservoir_mass_kg`; the per-width field is absent.
    pub checkpoint: GasFilmCheckpoint,
}

/// Advance a fixed disk (`inner_radius_m == 0`) or annulus using the shared
/// isothermal Reynolds solver. `input.grid.length_m` is its radial span; cell
/// samples lie at `inner_radius_m + (i+0.5)*dr`. The active prefix ends at its
/// last circular face, so contact-excluded outer cells contribute no storage/load.
///
/// A disk requires a sealed inner face. Use `MixedEnds { left: None, right:
/// Some(p) }` (with the full field names) for an ambient-open rim. Normal gap
/// motion/profiles, both film wall policies and ideal annular pressure vents are
/// supported. `AxisymmetricCapillaryVented`/`AxisymmetricCapillaryReservoir`
/// supply total tube counts and total finite volumes, with both radial ends
/// sealed. The selected positive-area cell homogenizes the exchange around its
/// circumference, including a central disk cell; it is not a point source at r=0.
/// Radial/azimuthal wall motion and per-width capillary/reservoir connections
/// are refused; they require different physical geometry/semantics.
///
/// Cell areas are exact annuli. Face resistance integrates 1/r between pressure
/// samples, preserving the constant-gap steady logarithmic pressure potential.
/// The axis has zero circumference and is never a fictitious pressure reservoir.
/// All inherited continuum, slope, Mach, positive-pressure and restart gates apply.
pub fn solve_isothermal_axisymmetric_gas_film_1d(
    input: &GasFilmInput,
    inner_radius_m: f64,
    checkpoint: Option<&GasFilmCheckpoint>,
) -> Result<AxisymmetricGasFilmStep, GasFilmError> {
    let step = solve_film(input, checkpoint, Some(inner_radius_m))?;
    let receipt = step.receipt;
    Ok(AxisymmetricGasFilmStep {
        absolute_pressure_pa: step.absolute_pressure_pa,
        checkpoint: step.checkpoint,
        receipt: AxisymmetricGasFilmReceipt {
            model_id: receipt.model_id,
            input_authority: receipt.input_authority,
            input_uncertainty: receipt.input_uncertainty,
            iterations: receipt.iterations,
            max_mass_residual_kg_m2_s: receipt.max_mass_residual_kg_m2_s,
            gas_mass_kg: receipt.gas_mass_kg_per_m,
            storage_rate_kg_per_s: receipt.storage_rate_kg_per_m_s,
            inner_boundary_outward_mass_flow_kg_per_s: receipt
                .left_boundary_outward_mass_flux_kg_per_m_s,
            outer_boundary_outward_mass_flow_kg_per_s: receipt
                .right_boundary_outward_mass_flux_kg_per_m_s,
            vent_outward_mass_flow_kg_per_s: receipt.vent_outward_mass_flux_kg_per_m_s,
            inner_boundary_outward_enthalpy_power_w: receipt
                .left_boundary_outward_enthalpy_flux_w_per_m,
            outer_boundary_outward_enthalpy_power_w: receipt
                .right_boundary_outward_enthalpy_flux_w_per_m,
            vent_outward_enthalpy_power_w: receipt.vent_outward_enthalpy_flux_w_per_m,
            capillary_vent: receipt
                .capillary_vent
                .map(|vent| AxisymmetricCapillaryVentReport {
                    viscous_dissipation_w: vent.viscous_dissipation_w_per_m,
                    reynolds_number: vent.reynolds_number,
                    maximum_knudsen_number: vent.maximum_knudsen_number,
                    maximum_gas_mach_number: vent.maximum_gas_mach_number,
                }),
            finite_reservoir: receipt.finite_reservoir.map(|reservoir| {
                AxisymmetricGasFilmReservoirReport {
                    absolute_pressure_pa: reservoir.absolute_pressure_pa,
                    gas_mass_kg: reservoir.gas_mass_kg_per_m,
                    storage_rate_kg_per_s: reservoir.storage_rate_kg_per_m_s,
                    inward_enthalpy_power_w: reservoir.inward_enthalpy_flux_w_per_m,
                    required_isothermal_heat_to_gas_w: reservoir
                        .required_isothermal_heat_to_gas_w_per_m,
                    combined_mass_closure_residual_kg_per_s: reservoir
                        .combined_mass_closure_residual_kg_per_m_s,
                }
            }),
            mass_closure_residual_kg_per_s: receipt.mass_closure_residual_kg_per_m_s,
            absolute_pressure_load_n: receipt.absolute_pressure_load_n_per_m,
            gauge_pressure_load_n: receipt.gauge_pressure_load_n_per_m,
            gas_on_wall_radial_traction_pa: receipt.gas_on_upper_wall_tangential_traction_pa,
            normal_gap_power_to_gas_w: receipt.normal_gap_power_to_gas_w_per_m,
            pressure_gradient_power_to_gas_w: receipt.pressure_gradient_power_to_gas_w_per_m,
            viscous_heat_w: receipt.viscous_heat_w_per_m,
            wall_slip_dissipation_w: receipt.wall_slip_dissipation_w_per_m,
            maximum_knudsen_number: receipt.maximum_knudsen_number,
            maximum_gas_mach_number: receipt.maximum_gas_mach_number,
            maximum_gap_slope: receipt.maximum_gap_slope,
        },
    })
}

pub(super) fn compile_metrics(
    input: &GasFilmInput,
    active: usize,
    spacing: f64,
    inner_radius_m: f64,
) -> Result<FilmMetrics, GasFilmError> {
    nonnegative(inner_radius_m, "axisymmetric.inner_radius_m")?;
    let outer_radius = derived(
        inner_radius_m + input.grid.length_m,
        "axisymmetric.outer_radius_m",
    )?;
    if outer_radius <= inner_radius_m {
        return Err(GasFilmError::InvalidInput {
            field: "axisymmetric.unresolved_radial_span",
        });
    }
    if input.wall_motion.lower_tangential_velocity_m_per_s != 0.0
        || input.wall_motion.upper_tangential_velocity_m_per_s != 0.0
    {
        return Err(GasFilmError::Unavailable {
            reason: "axisymmetric-in-plane-wall-motion",
        });
    }
    if matches!(
        input.boundary,
        GasFilmBoundaryTopology::CapillaryVented { .. }
            | GasFilmBoundaryTopology::CapillaryReservoir { .. }
    ) {
        return Err(GasFilmError::Unavailable {
            reason: "axisymmetric-per-width-capillary-connection",
        });
    }
    if inner_radius_m == 0.0 && input.boundary.end_pressures()[0].is_some() {
        return Err(GasFilmError::Unavailable {
            reason: "axisymmetric-axis-must-be-sealed",
        });
    }
    let mut areas = Vec::with_capacity(active);
    let mut faces = Vec::with_capacity(active + 1);
    for cell in 0..active {
        let radius = inner_radius_m + (cell as f64 + 0.5) * spacing;
        // Factored annular area avoids subtracting nearly equal squared radii.
        let area = derived(2.0 * PI * radius * spacing, "axisymmetric.cell_area_m2")?;
        positive(area, "axisymmetric.cell_area_m2")?;
        areas.push(area);
    }
    for face in 0..=active {
        let radius = inner_radius_m + face as f64 * spacing;
        let width = derived(2.0 * PI * radius, "axisymmetric.face_circumference_m")?;
        let distance_m = if radius == 0.0 {
            // No logarithm at the symmetry axis. Its flux and dual weight are zero.
            0.5 * spacing
        } else {
            let left = if face == 0 {
                radius
            } else {
                radius - 0.5 * spacing
            };
            let span = if face == 0 || face == active {
                0.5 * spacing
            } else {
                spacing
            };
            positive(left, "axisymmetric.left_pressure_radius_m")?;
            if left + span <= left {
                return Err(GasFilmError::InvalidInput {
                    field: "axisymmetric.unresolved_radial_cell",
                });
            }
            derived(
                radius * (span / left).ln_1p(),
                "axisymmetric.face_resistance_distance_m",
            )?
        };
        positive(distance_m, "axisymmetric.face_resistance_distance_m")?;
        faces.push(FilmFaceMetric { width, distance_m });
    }
    Ok(FilmMetrics::Axisymmetric {
        inner_radius_m,
        cell_area_m2: areas,
        faces,
    })
}
