//! Reduced, deterministic circular-capillary flow screening.
//!
//! This module owns the Hagen–Poiseuille relation for a declared Newtonian,
//! incompressible fluid in a straight, rigid, circular capillary with
//! fully-developed laminar flow and no slip. It deliberately does not turn
//! that screening law into a polymer-rheology claim: entrance contraction,
//! shear thinning, viscoelasticity, wall slip and heating are outside this rung.
//! Lucas–Washburn meniscus models admit horizontal and inclined capillary
//! imbibition/drainage with an explicit equilibrium angle and energy account.

use core::f64::consts::PI;

/// Conservative triangle-surface forces with exact endpoint area work.
pub mod surface;

/// SI inputs admitted by the circular-capillary screening law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularCapillaryInput {
    /// Dynamic viscosity in pascal-seconds.
    pub dynamic_viscosity_pa_s: f64,
    /// Straight capillary land length in metres.
    pub length_m: f64,
    /// Circular capillary radius in metres.
    pub radius_m: f64,
    /// Volumetric flow rate in cubic metres per second.
    pub volumetric_flow_m3_s: f64,
}

/// Deterministic output of the admitted reduced law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularCapillaryStep {
    /// Required pressure drop in pascals.
    pub pressure_drop_pa: f64,
    /// Newtonian wall shear rate, `4Q/(πR³)`, in reciprocal seconds.
    pub wall_shear_rate_per_s: f64,
    /// Hydraulic power `ΔP Q` in watts.
    pub hydraulic_power_w: f64,
}

/// Typed refusal at the reduced-law boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CircularCapillaryError {
    /// A named input was NaN or infinite.
    NonFinite {
        /// Stable input field name.
        field: &'static str,
    },
    /// A named input that must be strictly positive was zero or negative.
    NonPositive {
        /// Stable input field name.
        field: &'static str,
    },
    /// Reverse flow is outside this scalar browser-facing screening rung.
    NegativeFlow,
    /// A finite value lies outside its stated physical or geometric domain.
    OutsideDomain {
        /// Stable input field name.
        field: &'static str,
    },
    /// The admitted calculation overflowed or otherwise became non-finite.
    NonFiniteResult,
}

/// Evaluate Hagen–Poiseuille flow within the declared reduced-law boundary.
///
/// Zero flow is admitted and returns zero pressure, shear rate, and power.
pub fn step_newtonian_circular_capillary(
    input: CircularCapillaryInput,
) -> Result<CircularCapillaryStep, CircularCapillaryError> {
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("length_m", input.length_m),
        ("radius_m", input.radius_m),
        ("volumetric_flow_m3_s", input.volumetric_flow_m3_s),
    ] {
        if !value.is_finite() {
            return Err(CircularCapillaryError::NonFinite { field });
        }
    }
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("length_m", input.length_m),
        ("radius_m", input.radius_m),
    ] {
        if value <= 0.0 {
            return Err(CircularCapillaryError::NonPositive { field });
        }
    }
    if input.volumetric_flow_m3_s < 0.0 {
        return Err(CircularCapillaryError::NegativeFlow);
    }

    let radius_cubed = input.radius_m.powi(3);
    let pressure_drop_pa =
        8.0 * input.dynamic_viscosity_pa_s * input.length_m * input.volumetric_flow_m3_s
            / (PI * input.radius_m.powi(4));
    let wall_shear_rate_per_s = 4.0 * input.volumetric_flow_m3_s / (PI * radius_cubed);
    let hydraulic_power_w = pressure_drop_pa * input.volumetric_flow_m3_s;

    if !pressure_drop_pa.is_finite()
        || !wall_shear_rate_per_s.is_finite()
        || !hydraulic_power_w.is_finite()
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    Ok(CircularCapillaryStep {
        pressure_drop_pa,
        wall_shear_rate_per_s,
        hydraulic_power_w,
    })
}

/// SI inputs for a horizontal, partially filled, circular tube.
///
/// The liquid is Newtonian and incompressible; flow is quasisteady, laminar,
/// and no-slip. Gas is vented at its prescribed pressure. Gravity, inertia,
/// entrance losses, evaporation, gas resistance and dynamic-angle hysteresis
/// are excluded. The equilibrium angle belongs to the actual solid/liquid/gas
/// system; no material name selects a coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularMeniscusInput {
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Tube radius [m].
    pub radius_m: f64,
    /// Available tube length [m].
    pub tube_length_m: f64,
    /// Current liquid-column length [m], in [0, tube_length_m].
    pub wet_length_m: f64,
    /// Nonnegative liquid/gas surface tension [N/m]; zero explicitly omits it.
    pub surface_tension_n_m: f64,
    /// Cosine of the equilibrium angle measured through the liquid, in [-1, 1].
    /// Positive values wet the wall; negative values resist entry.
    pub contact_angle_cosine: f64,
    /// Reservoir liquid pressure minus downstream gas pressure [Pa], signed.
    pub applied_pressure_pa: f64,
    /// Requested nonnegative integration interval [s].
    pub duration_s: f64,
}

/// A meniscus reaching an open tube end requires a new boundary model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CircularMeniscusEnd {
    /// The complete requested interval was integrated without an end event.
    TimeLimit,
    /// The liquid column emptied; no negative liquid volume is introduced.
    Inlet,
    /// The column reached the outlet; throughflow is a separate model.
    Outlet,
}

/// Accepted moving-meniscus state and its signed physical transfers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularMeniscusStep {
    /// Accepted liquid-column length [m].
    pub wet_length_m: f64,
    /// Time actually advanced [s]; an end event can stop before duration_s.
    pub elapsed_time_s: f64,
    /// Reason the advance ended.
    pub end: CircularMeniscusEnd,
    /// Signed volume transferred from the reservoir into the tube [m³].
    pub transferred_volume_m3: f64,
    /// Gas minus liquid pressure at the meniscus, 2 sigma cos(theta)/R [Pa].
    pub capillary_pressure_pa: f64,
    /// Work by the imposed reservoir-to-gas pressure difference [J].
    pub pressure_work_j: f64,
    /// Change of solid/liquid minus solid/gas surface free energy [J].
    /// The translating meniscus has constant area; end topology is excluded.
    pub wall_surface_energy_change_j: f64,
    /// Integral of positive viscous loss over the accepted advance [J].
    pub viscous_dissipation_j: f64,
}

/// Integrate the horizontal Lucas–Washburn law exactly for constant inputs.
///
/// With drive = applied pressure + 2 sigma cos(theta)/R,
/// d(length²)/dt = R² drive/(4 mu). Integrating length squared avoids
/// dividing by zero at an initially empty tube. It does not resolve the
/// inertial startup omitted by the quasisteady model. Reverse flow and
/// nonwetting entry pressure are retained. End events consume their actual
/// travel time; residual requested time is never silently throughflow.
///
/// Young's law gives delta E_wall = -2 pi R sigma cos(theta) delta length;
/// imposed pressure work and viscous loss close W - delta E = D >= 0.
/// Initial meniscus formation and final end topology are outside this
/// translating-meniscus energy account.
///
/// Reference: E. W. Washburn, Phys. Rev. 17, 273 (1921),
/// <https://doi.org/10.1103/PhysRev.17.273> (no-slip horizontal limit).
///
/// # Errors
/// Refuses nonfinite/out-of-domain inputs and unrepresentable geometry,
/// mobility or transfers; no caller state is mutated on refusal.
pub fn step_newtonian_circular_meniscus(
    input: CircularMeniscusInput,
) -> Result<CircularMeniscusStep, CircularCapillaryError> {
    validate_circular_meniscus(input)?;
    let area = PI * input.radius_m * input.radius_m;
    let length_squared = input.wet_length_m * input.wet_length_m;
    let tube_squared = input.tube_length_m * input.tube_length_m;
    let capillary_pressure_pa = if input.contact_angle_cosine == 0.0 {
        0.0
    } else {
        (input.surface_tension_n_m / input.radius_m) * (2.0 * input.contact_angle_cosine)
    };
    let drive = input.applied_pressure_pa + capillary_pressure_pa;
    let rate = (0.25 * input.radius_m) * input.radius_m * (drive / input.dynamic_viscosity_pa_s);
    if !area.is_finite()
        || area <= 0.0
        || !tube_squared.is_finite()
        || tube_squared <= 0.0
        || (input.wet_length_m > 0.0 && length_squared == 0.0)
        || !capillary_pressure_pa.is_finite()
        || !drive.is_finite()
        || !rate.is_finite()
        || (drive != 0.0 && rate == 0.0)
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    let (wet_length_m, elapsed_time_s, end) = if rate == 0.0 || input.duration_s == 0.0 {
        (
            input.wet_length_m,
            input.duration_s,
            CircularMeniscusEnd::TimeLimit,
        )
    } else {
        let (distance_squared, boundary, event) = if rate > 0.0 {
            // Factored subtraction retains the small remaining distance.
            (
                (input.tube_length_m - input.wet_length_m)
                    * (input.tube_length_m + input.wet_length_m),
                input.tube_length_m,
                CircularMeniscusEnd::Outlet,
            )
        } else {
            (length_squared, 0.0, CircularMeniscusEnd::Inlet)
        };
        let travel_time = distance_squared / rate.abs();
        if distance_squared > 0.0 && travel_time == 0.0 {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        if travel_time <= input.duration_s {
            (boundary, travel_time, event)
        } else {
            (
                (rate.mul_add(input.duration_s, length_squared)).sqrt(),
                input.duration_s,
                CircularMeniscusEnd::TimeLimit,
            )
        }
    };
    let transferred_volume_m3 = area * (wet_length_m - input.wet_length_m);
    let pressure_work_j = input.applied_pressure_pa * transferred_volume_m3;
    let wall_surface_energy_change_j = -capillary_pressure_pa * transferred_volume_m3;
    let viscous_dissipation_j = drive * transferred_volume_m3;
    if !wet_length_m.is_finite()
        || !elapsed_time_s.is_finite()
        || !transferred_volume_m3.is_finite()
        || !pressure_work_j.is_finite()
        || !wall_surface_energy_change_j.is_finite()
        || !viscous_dissipation_j.is_finite()
        || viscous_dissipation_j < 0.0
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    Ok(CircularMeniscusStep {
        wet_length_m,
        elapsed_time_s,
        end,
        transferred_volume_m3,
        capillary_pressure_pa,
        pressure_work_j,
        wall_surface_energy_change_j,
        viscous_dissipation_j,
    })
}

/// SI inputs for capillary motion in a straight tube at a fixed inclination.
///
/// The tube inlet is at the reservoir free-surface elevation. The projected
/// acceleration is positive when gravity opposes advancing liquid, negative
/// for a downward tube, and zero for a horizontal tube. Gas is vented at a
/// prescribed uniform pressure. The same quasisteady Newtonian/no-slip model
/// as `CircularMeniscusInput` applies; inertia and dynamic wetting are excluded.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InclinedCircularMeniscusInput {
    /// Geometry, liquid properties, pressure forcing, state and requested time.
    pub meniscus: CircularMeniscusInput,
    /// Liquid density [kg/m³], held constant during the advance.
    pub density_kg_m3: f64,
    /// Opposing gravitational acceleration along the tube [m/s²], signed.
    /// For upward inclination alpha from horizontal, this is g sin(alpha).
    pub opposing_gravity_m_s2: f64,
}

/// Gravity-aware meniscus transport and its mechanical energy account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InclinedCircularMeniscusStep {
    /// Accepted state, pressure work, wall energy and actual viscous loss.
    pub meniscus: CircularMeniscusStep,
    /// Change of the column's gravitational potential relative to the inlet [J].
    pub gravitational_energy_change_j: f64,
    /// Hydrostatic pressure opposing advance at the accepted meniscus [Pa].
    pub hydrostatic_pressure_pa: f64,
}

/// Integrate capillary rise or drainage with hydrostatic feedback.
///
/// With A = delta_p + 2 sigma cos(theta)/R, b = rho g_parallel and
/// K = R²/(8 mu), the equation is x dx/dt = K (A - b x). For upward
/// wetting, it approaches Jurin's height A/b without stepping across it.
/// Pressure-driven reversal and gravity-assisted downward flow are retained.
///
/// The separated analytic travel-time integral is inverted monotonically
/// with a bounded safeguarded solve; there is no explicit time-discretization
/// error. Endpoint weights and log1p avoid cancellation in the horizontal
/// limit and near equilibrium. The zero-gravity path is exactly the existing
/// horizontal solver. The computation allocates nothing.
///
/// The account closes W_pressure - delta E_wall - delta E_gravity = D >= 0.
/// As in the horizontal model, meniscus formation, changing end topology,
/// liquid/gas inertia and dynamic-contact-angle dissipation are excluded.
/// At A=0 the continuous quasisteady branch has dx/dt=-K b; it does not
/// resolve the inertial start of a downward initially empty column.
///
/// Reference: Fries & Dreyer, J. Colloid Interface Sci. 320, 259–263 (2008),
/// <https://doi.org/10.1016/j.jcis.2008.01.009>.
///
/// # Errors
/// Refuses invalid input or unrepresentable coefficients/transfers. No caller
/// state is mutated. The finite scalar solve has no internal cancellation point.
pub fn step_newtonian_inclined_circular_meniscus(
    input: InclinedCircularMeniscusInput,
) -> Result<InclinedCircularMeniscusStep, CircularCapillaryError> {
    validate_circular_meniscus(input.meniscus)?;
    for (field, value) in [
        ("density_kg_m3", input.density_kg_m3),
        ("opposing_gravity_m_s2", input.opposing_gravity_m_s2),
    ] {
        if !value.is_finite() {
            return Err(CircularCapillaryError::NonFinite { field });
        }
    }
    if input.density_kg_m3 <= 0.0 {
        return Err(CircularCapillaryError::NonPositive {
            field: "density_kg_m3",
        });
    }
    if input.opposing_gravity_m_s2 == 0.0 {
        return Ok(InclinedCircularMeniscusStep {
            meniscus: step_newtonian_circular_meniscus(input.meniscus)?,
            gravitational_energy_change_j: 0.0,
            hydrostatic_pressure_pa: 0.0,
        });
    }
    let m = input.meniscus;
    let area = PI * m.radius_m * m.radius_m;
    let mobility = (m.radius_m / m.dynamic_viscosity_pa_s) * (m.radius_m / 8.0);
    let b = input.density_kg_m3 * input.opposing_gravity_m_s2;
    let capillary = if m.contact_angle_cosine == 0.0 {
        0.0
    } else {
        (m.surface_tension_n_m / m.radius_m) * (2.0 * m.contact_angle_cosine)
    };
    let drive = m.applied_pressure_pa + capillary;
    let initial_drive = (-b).mul_add(m.wet_length_m, drive);
    if !area.is_finite()
        || area <= 0.0
        || !mobility.is_finite()
        || mobility <= 0.0
        || !b.is_finite()
        || b == 0.0
        || !capillary.is_finite()
        || !drive.is_finite()
        || !initial_drive.is_finite()
        || !(b * m.wet_length_m).is_finite()
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    let (length, elapsed, end) = if m.duration_s == 0.0 {
        (m.wet_length_m, 0.0, CircularMeniscusEnd::TimeLimit)
    } else if drive == 0.0 {
        // Here x cancels from the ODE, including the limiting arrival at x=0.
        let velocity = -mobility * b;
        if !velocity.is_finite() || velocity == 0.0 {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        let (boundary, event) = if velocity > 0.0 {
            (m.tube_length_m, CircularMeniscusEnd::Outlet)
        } else {
            (0.0, CircularMeniscusEnd::Inlet)
        };
        let travel = (boundary - m.wet_length_m) / velocity;
        if travel <= m.duration_s {
            (boundary, travel, event)
        } else {
            (
                velocity.mul_add(m.duration_s, m.wet_length_m),
                m.duration_s,
                CircularMeniscusEnd::TimeLimit,
            )
        }
    } else if initial_drive == 0.0 || drive / b == m.wet_length_m {
        (m.wet_length_m, m.duration_s, CircularMeniscusEnd::TimeLimit)
    } else {
        advance_inclined_meniscus(m, mobility, b, drive, initial_drive)?
    };

    let volume = area * (length - m.wet_length_m);
    let final_drive = (-b).mul_add(length, drive);
    let pressure_work = m.applied_pressure_pa * volume;
    let wall_energy = -capillary * volume;
    let mean_length = m.wet_length_m + 0.5 * (length - m.wet_length_m);
    let gravity_energy = (b * mean_length) * volume;
    // The drive is affine in x: its exact displacement average gives the
    // dissipation without subtracting two nearly equal large energy accounts.
    let dissipation = (0.5 * initial_drive + 0.5 * final_drive) * volume;
    let hydrostatic = b * length;
    if [
        length,
        elapsed,
        volume,
        pressure_work,
        wall_energy,
        gravity_energy,
        dissipation,
        hydrostatic,
    ]
    .into_iter()
    .any(|v| !v.is_finite())
        || !(0.0..=m.tube_length_m).contains(&length)
        || dissipation < 0.0
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    Ok(InclinedCircularMeniscusStep {
        meniscus: CircularMeniscusStep {
            wet_length_m: length,
            elapsed_time_s: elapsed,
            end,
            transferred_volume_m3: volume,
            capillary_pressure_pa: capillary,
            pressure_work_j: pressure_work,
            wall_surface_energy_change_j: wall_energy,
            viscous_dissipation_j: dissipation,
        },
        gravitational_energy_change_j: gravity_energy,
        hydrostatic_pressure_pa: hydrostatic,
    })
}

fn advance_inclined_meniscus(
    m: CircularMeniscusInput,
    mobility: f64,
    b: f64,
    drive: f64,
    initial_drive: f64,
) -> Result<(f64, f64, CircularMeniscusEnd), CircularCapillaryError> {
    let direction = initial_drive.signum();
    let (boundary, event) = if direction > 0.0 {
        (m.tube_length_m, CircularMeniscusEnd::Outlet)
    } else {
        (0.0, CircularMeniscusEnd::Inlet)
    };
    if boundary == m.wet_length_m {
        return Ok((boundary, 0.0, event));
    }
    let equilibrium = drive / b;
    let has_barrier = b > 0.0
        && equilibrium >= 0.0
        && direction * (equilibrium - m.wet_length_m) > 0.0
        && direction * (boundary - equilibrium) >= 0.0;
    let target = if has_barrier { equilibrium } else { boundary };
    let limit = (target - m.wet_length_m).abs();
    let time_at = |distance: f64| {
        inclined_meniscus_travel_time(
            m.wet_length_m,
            direction * distance,
            mobility,
            b,
            initial_drive,
        )
    };
    if !has_barrier {
        let travel = time_at(limit)?;
        if travel <= m.duration_s {
            return Ok((boundary, travel, event));
        }
    }

    // A constant force equal to the largest force on this segment bounds the
    // displacement. It also keeps tiny timesteps local instead of bisecting
    // the entire tube to locate a microscopic displacement.
    let maximum_drive = initial_drive.abs().max((-b).mul_add(target, drive).abs());
    let increment = (2.0 * mobility * maximum_drive) * m.duration_s;
    let initial_squared = m.wet_length_m * m.wet_length_m;
    let estimate = if increment.is_finite() && initial_squared.is_finite() {
        if direction > 0.0 {
            increment / ((initial_squared + increment).sqrt() + m.wet_length_m)
        } else if increment >= initial_squared {
            m.wet_length_m
        } else {
            increment / ((initial_squared - increment).sqrt() + m.wet_length_m)
        }
    } else {
        limit
    };
    let mut high = if estimate.is_finite() && estimate > 0.0 {
        estimate.min(limit)
    } else {
        limit
    };
    if high < limit && time_at(high)? < m.duration_s {
        high = limit; // Preserve the bracket if the analytic bound rounded down.
    }
    let mut low = 0.0;
    let mut distance = 0.5 * high;
    for iteration in 0..80 {
        let travel = time_at(distance)?;
        if (travel - m.duration_s).abs() <= 8.0 * f64::EPSILON * m.duration_s {
            break;
        }
        if travel < m.duration_s {
            low = distance;
        } else {
            high = distance;
        }
        let midpoint = low + 0.5 * (high - low);
        if midpoint == low
            || midpoint == high
            || m.wet_length_m + direction * low == m.wet_length_m + direction * high
        {
            distance = midpoint;
            break;
        }
        // Newton accelerates normal cases; the remaining 64 iterations are
        // bisection, so a poorly conditioned equilibrium cannot stall progress.
        let position = m.wet_length_m + direction * distance;
        let force = (-b).mul_add(direction * distance, initial_drive).abs();
        let candidate = distance - (travel - m.duration_s) * (mobility * force / position);
        distance = if iteration < 16 && candidate.is_finite() && candidate > low && candidate < high
        {
            candidate
        } else {
            midpoint
        };
    }
    let length = (m.wet_length_m + direction * distance)
        .clamp(m.wet_length_m.min(target), m.wet_length_m.max(target));
    Ok((length, m.duration_s, CircularMeniscusEnd::TimeLimit))
}

/// Integral of x/(K (A-bx)) from x0 to x0+delta, without subtracting
/// antiderivatives. Both endpoint weights are positive for every z<1.
fn inclined_meniscus_travel_time(
    initial: f64,
    delta: f64,
    mobility: f64,
    b: f64,
    initial_drive: f64,
) -> Result<f64, CircularCapillaryError> {
    if delta == 0.0 {
        return Ok(0.0);
    }
    let ratio = delta / initial_drive;
    let z = b * ratio;
    if z >= 1.0 {
        return Ok(f64::INFINITY);
    }
    if !z.is_finite() || !ratio.is_finite() {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    let (left, right) = if z.abs() <= 0.125 {
        // Integral weights of (1-u)/(1-zu) and u/(1-zu). At z=0 both
        // are 1/2, giving the exact horizontal squared-length relation.
        let mut left = 0.5;
        let mut right = 0.5;
        let mut power = 1.0;
        for k in 1..=20 {
            power *= z;
            right += power / f64::from(k + 2);
            left += power / f64::from((k + 1) * (k + 2));
        }
        (left, right)
    } else {
        let logarithm = (-z).ln_1p();
        (
            (1.0 + (1.0 / z - 1.0) * logarithm) / z,
            (-logarithm / z - 1.0) / z,
        )
    };
    let time = ratio * (initial * left + (initial + delta) * right) / mobility;
    if time.is_nan() || time < 0.0 {
        Err(CircularCapillaryError::NonFiniteResult)
    } else {
        Ok(time)
    }
}

fn validate_circular_meniscus(input: CircularMeniscusInput) -> Result<(), CircularCapillaryError> {
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("radius_m", input.radius_m),
        ("tube_length_m", input.tube_length_m),
        ("wet_length_m", input.wet_length_m),
        ("surface_tension_n_m", input.surface_tension_n_m),
        ("contact_angle_cosine", input.contact_angle_cosine),
        ("applied_pressure_pa", input.applied_pressure_pa),
        ("duration_s", input.duration_s),
    ] {
        if !value.is_finite() {
            return Err(CircularCapillaryError::NonFinite { field });
        }
    }
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("radius_m", input.radius_m),
        ("tube_length_m", input.tube_length_m),
    ] {
        if value <= 0.0 {
            return Err(CircularCapillaryError::NonPositive { field });
        }
    }
    for (field, admitted) in [
        (
            "wet_length_m",
            (0.0..=input.tube_length_m).contains(&input.wet_length_m),
        ),
        ("surface_tension_n_m", input.surface_tension_n_m >= 0.0),
        (
            "contact_angle_cosine",
            (-1.0..=1.0).contains(&input.contact_angle_cosine),
        ),
        ("duration_s", input.duration_s >= 0.0),
    ] {
        if !admitted {
            return Err(CircularCapillaryError::OutsideDomain { field });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meniscus() -> CircularMeniscusInput {
        CircularMeniscusInput {
            dynamic_viscosity_pa_s: 0.1,
            radius_m: 1e-3,
            tube_length_m: 0.1,
            wet_length_m: 0.0,
            surface_tension_n_m: 0.072,
            contact_angle_cosine: 0.5,
            applied_pressure_pa: 0.0,
            duration_s: 0.01,
        }
    }

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-20),
            "actual={actual}, expected={expected}"
        );
    }

    fn inclined() -> InclinedCircularMeniscusInput {
        InclinedCircularMeniscusInput {
            meniscus: meniscus(),
            density_kg_m3: 1000.0,
            opposing_gravity_m_s2: 9.81,
        }
    }

    fn capillary_equilibrium(input: InclinedCircularMeniscusInput) -> (f64, f64) {
        let m = input.meniscus;
        let pressure = m.applied_pressure_pa
            + 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m;
        let b = input.density_kg_m3 * input.opposing_gravity_m_s2;
        let height = pressure / b;
        let time_scale =
            8.0 * m.dynamic_viscosity_pa_s * pressure / (m.radius_m * m.radius_m * b * b);
        (height, time_scale)
    }

    fn gravity_energy_closes(result: InclinedCircularMeniscusStep) {
        let m = result.meniscus;
        near(
            m.pressure_work_j
                - m.wall_surface_energy_change_j
                - result.gravitational_energy_change_j,
            m.viscous_dissipation_j,
        );
        assert!(m.viscous_dissipation_j >= 0.0);
    }

    #[test]
    fn inclined_g1_rise_matches_independent_dimensionless_washburn_solution() {
        let input = inclined();
        let (height, scale) = capillary_equilibrium(input);
        // The published upward-rise solution is t/tau=-y-ln(1-y), y=h/h_J.
        for fraction in [0.01_f64, 0.25, 0.5, 0.9, 0.99] {
            let step = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    duration_s: scale * (-fraction - (-fraction).ln_1p()),
                    ..input.meniscus
                },
                ..input
            })
            .unwrap();
            near(step.meniscus.wet_length_m, fraction * height);
            assert_eq!(step.meniscus.end, CircularMeniscusEnd::TimeLimit);
            near(
                step.hydrostatic_pressure_pa,
                input.density_kg_m3 * input.opposing_gravity_m_s2 * step.meniscus.wet_length_m,
            );
            let potential = 0.5
                * input.density_kg_m3
                * input.opposing_gravity_m_s2
                * PI
                * input.meniscus.radius_m.powi(2)
                * step.meniscus.wet_length_m.powi(2);
            near(step.gravitational_energy_change_j, potential);
            gravity_energy_closes(step);
        }
    }

    #[test]
    fn inclined_g1_jurin_arrests_rise_and_drainage_without_overshoot() {
        let input = inclined();
        let (height, scale) = capillary_equilibrium(input);
        for initial_fraction in [0.0, 0.5, 1.0, 1.5, 3.0] {
            let initial = initial_fraction * height;
            let step = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    wet_length_m: initial,
                    duration_s: 100.0 * scale,
                    ..input.meniscus
                },
                ..input
            })
            .unwrap();
            near(step.meniscus.wet_length_m, height);
            assert!(
                (initial.min(height)..=initial.max(height)).contains(&step.meniscus.wet_length_m)
            );
            assert_eq!(step.meniscus.end, CircularMeniscusEnd::TimeLimit);
            if initial_fraction == 1.0 {
                assert_eq!(step.meniscus.transferred_volume_m3, 0.0);
                assert_eq!(step.meniscus.viscous_dissipation_j, 0.0);
            } else {
                gravity_energy_closes(step);
            }
        }
        // Above the equilibrium, G(y)=-y-ln(y-1) gives drainage travel time.
        let initial = 1.6_f64;
        let final_fraction = 1.2_f64;
        let time =
            scale * (initial - final_fraction + ((initial - 1.0) / (final_fraction - 1.0)).ln());
        let drained = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
            meniscus: CircularMeniscusInput {
                wet_length_m: initial * height,
                duration_s: time,
                ..input.meniscus
            },
            ..input
        })
        .unwrap();
        near(drained.meniscus.wet_length_m, final_fraction * height);
        assert!(drained.meniscus.transferred_volume_m3 < 0.0);
        assert!(drained.gravitational_energy_change_j < 0.0);
        gravity_energy_closes(drained);
    }

    #[test]
    fn inclined_g3_horizontal_limit_and_timestep_composition() {
        let base = inclined();
        for gravity in [0.0, 1e-12, -1e-12] {
            let input = InclinedCircularMeniscusInput {
                opposing_gravity_m_s2: gravity,
                ..base
            };
            let result = step_newtonian_inclined_circular_meniscus(input).unwrap();
            let horizontal = step_newtonian_circular_meniscus(input.meniscus).unwrap();
            if gravity == 0.0 {
                assert_eq!(result.meniscus, horizontal);
                assert_eq!(result.gravitational_energy_change_j, 0.0);
            } else {
                near(result.meniscus.wet_length_m, horizontal.wet_length_m);
            }
        }
        for gravity in [9.81, -9.81] {
            let input = InclinedCircularMeniscusInput {
                opposing_gravity_m_s2: gravity,
                meniscus: CircularMeniscusInput {
                    wet_length_m: 0.002,
                    applied_pressure_pa: 17.0,
                    duration_s: 0.08,
                    ..base.meniscus
                },
                ..base
            };
            let full = step_newtonian_inclined_circular_meniscus(input).unwrap();
            let first = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    duration_s: input.meniscus.duration_s / 2.0,
                    ..input.meniscus
                },
                ..input
            })
            .unwrap();
            let second = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    wet_length_m: first.meniscus.wet_length_m,
                    duration_s: input.meniscus.duration_s / 2.0,
                    ..input.meniscus
                },
                ..input
            })
            .unwrap();
            near(second.meniscus.wet_length_m, full.meniscus.wet_length_m);
            near(
                first.gravitational_energy_change_j + second.gravitational_energy_change_j,
                full.gravitational_energy_change_j,
            );
            near(
                first.meniscus.viscous_dissipation_j + second.meniscus.viscous_dissipation_j,
                full.meniscus.viscous_dissipation_j,
            );
            gravity_energy_closes(full);
        }
    }

    #[test]
    fn inclined_end_event_uses_actual_time_and_inclination_changes_height() {
        let input = inclined();
        let (height, scale) = capillary_equilibrium(input);
        let result = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
            meniscus: CircularMeniscusInput {
                tube_length_m: 0.5 * height,
                duration_s: 20.0,
                ..input.meniscus
            },
            ..input
        })
        .unwrap();
        assert_eq!(result.meniscus.end, CircularMeniscusEnd::Outlet);
        near(
            result.meniscus.elapsed_time_s,
            scale * (-0.5 - 0.5_f64.ln()),
        );
        near(result.meniscus.wet_length_m, 0.5 * height);
        gravity_energy_closes(result);
        let half_gravity =
            step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                opposing_gravity_m_s2: 0.5 * input.opposing_gravity_m_s2,
                meniscus: CircularMeniscusInput {
                    duration_s: 1000.0,
                    ..input.meniscus
                },
                ..input
            })
            .unwrap();
        near(half_gravity.meniscus.wet_length_m, 2.0 * height);
    }

    #[test]
    fn inclined_nonwetting_downhill_entry_threshold_and_gravity_only_limit() {
        let base = inclined();
        let downhill = InclinedCircularMeniscusInput {
            opposing_gravity_m_s2: -9.81,
            meniscus: CircularMeniscusInput {
                contact_angle_cosine: -0.5,
                duration_s: 100.0,
                ..base.meniscus
            },
            ..base
        };
        let (threshold, _) = capillary_equilibrium(downhill);
        for (initial, expected_end) in [
            (0.5 * threshold, CircularMeniscusEnd::Inlet),
            (1.5 * threshold, CircularMeniscusEnd::Outlet),
        ] {
            let step = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    wet_length_m: initial,
                    ..downhill.meniscus
                },
                ..downhill
            })
            .unwrap();
            assert_eq!(step.meniscus.end, expected_end);
            assert!(step.meniscus.elapsed_time_s < downhill.meniscus.duration_s);
            gravity_energy_closes(step);
        }
        for gravity in [9.81, -9.81] {
            let input = InclinedCircularMeniscusInput {
                opposing_gravity_m_s2: gravity,
                meniscus: CircularMeniscusInput {
                    surface_tension_n_m: 0.0,
                    wet_length_m: 0.02,
                    duration_s: 0.1,
                    ..base.meniscus
                },
                ..base
            };
            let step = step_newtonian_inclined_circular_meniscus(input).unwrap();
            let velocity = -input.density_kg_m3 * gravity * input.meniscus.radius_m.powi(2)
                / (8.0 * input.meniscus.dynamic_viscosity_pa_s);
            near(
                step.meniscus.wet_length_m,
                input.meniscus.wet_length_m + velocity * input.meniscus.duration_s,
            );
            assert_eq!(step.meniscus.pressure_work_j, 0.0);
            assert_eq!(step.meniscus.wall_surface_energy_change_j, 0.0);
            gravity_energy_closes(step);
        }
    }

    #[test]
    fn inclined_invalid_density_gravity_and_duration_refuse() {
        let input = inclined();
        for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                    density_kg_m3: density,
                    ..input
                },)
                .is_err()
            );
        }
        for gravity in [f64::NAN, f64::INFINITY] {
            assert!(
                step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                    opposing_gravity_m_s2: gravity,
                    ..input
                },)
                .is_err()
            );
        }
        for duration in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
                    meniscus: CircularMeniscusInput {
                        duration_s: duration,
                        ..input.meniscus
                    },
                    ..input
                },)
                .is_err()
            );
        }
    }

    #[test]
    fn washburn_g1_capillary_penetration_and_surface_energy() {
        let input = meniscus();
        let step = step_newtonian_circular_meniscus(input).unwrap();
        let expected_squared = input.surface_tension_n_m
            * input.radius_m
            * input.contact_angle_cosine
            * input.duration_s
            / (2.0 * input.dynamic_viscosity_pa_s);
        near(step.wet_length_m.powi(2), expected_squared);
        near(step.capillary_pressure_pa, 72.0);
        near(
            step.transferred_volume_m3,
            PI * input.radius_m.powi(2) * step.wet_length_m,
        );
        let surface_change = -2.0
            * PI
            * input.radius_m
            * input.surface_tension_n_m
            * input.contact_angle_cosine
            * step.wet_length_m;
        near(step.wall_surface_energy_change_j, surface_change);
        near(step.viscous_dissipation_j, -surface_change);
        assert_eq!(step.pressure_work_j, 0.0);
        assert_eq!(step.end, CircularMeniscusEnd::TimeLimit);
    }

    #[test]
    fn washburn_g3_split_steps_and_parameter_scaling() {
        let input = meniscus();
        let full = step_newtonian_circular_meniscus(input).unwrap();
        let half = step_newtonian_circular_meniscus(CircularMeniscusInput {
            duration_s: input.duration_s / 2.0,
            ..input
        })
        .unwrap();
        let second = step_newtonian_circular_meniscus(CircularMeniscusInput {
            wet_length_m: half.wet_length_m,
            duration_s: input.duration_s / 2.0,
            ..input
        })
        .unwrap();
        near(second.wet_length_m, full.wet_length_m);
        near(
            half.transferred_volume_m3 + second.transferred_volume_m3,
            full.transferred_volume_m3,
        );
        near(
            half.viscous_dissipation_j + second.viscous_dissipation_j,
            full.viscous_dissipation_j,
        );
        for changed in [
            CircularMeniscusInput {
                radius_m: 4.0 * input.radius_m,
                ..input
            },
            CircularMeniscusInput {
                surface_tension_n_m: 4.0 * input.surface_tension_n_m,
                ..input
            },
            CircularMeniscusInput {
                dynamic_viscosity_pa_s: input.dynamic_viscosity_pa_s / 4.0,
                ..input
            },
        ] {
            near(
                step_newtonian_circular_meniscus(changed)
                    .unwrap()
                    .wet_length_m,
                2.0 * full.wet_length_m,
            );
        }
    }

    #[test]
    fn washburn_nonwetting_drainage_and_pressure_work_close() {
        let input = CircularMeniscusInput {
            wet_length_m: 0.01,
            contact_angle_cosine: -0.5,
            applied_pressure_pa: 20.0,
            ..meniscus()
        };
        let step = step_newtonian_circular_meniscus(input).unwrap();
        assert!(step.wet_length_m < input.wet_length_m);
        assert!(step.transferred_volume_m3 < 0.0 && step.pressure_work_j < 0.0);
        assert!(step.wall_surface_energy_change_j < 0.0 && step.viscous_dissipation_j > 0.0);
        near(
            step.pressure_work_j - step.wall_surface_energy_change_j,
            step.viscous_dissipation_j,
        );
        let forced = step_newtonian_circular_meniscus(CircularMeniscusInput {
            applied_pressure_pa: 144.0,
            ..input
        })
        .unwrap();
        assert!(forced.wet_length_m > input.wet_length_m);
        assert!(forced.pressure_work_j > 0.0 && forced.wall_surface_energy_change_j > 0.0);
        near(
            forced.pressure_work_j - forced.wall_surface_energy_change_j,
            forced.viscous_dissipation_j,
        );
    }

    #[test]
    fn washburn_end_events_consume_only_travel_time() {
        let input = CircularMeniscusInput {
            tube_length_m: 0.002,
            duration_s: 1.0,
            ..meniscus()
        };
        let filled = step_newtonian_circular_meniscus(input).unwrap();
        assert_eq!(filled.wet_length_m, input.tube_length_m);
        assert_eq!(filled.end, CircularMeniscusEnd::Outlet);
        near(
            filled.elapsed_time_s,
            2.0 * input.dynamic_viscosity_pa_s * input.tube_length_m.powi(2)
                / (input.surface_tension_n_m * input.radius_m * input.contact_angle_cosine),
        );
        assert!(filled.elapsed_time_s < input.duration_s);
        let drained = step_newtonian_circular_meniscus(CircularMeniscusInput {
            wet_length_m: input.tube_length_m,
            contact_angle_cosine: -0.5,
            ..input
        })
        .unwrap();
        assert_eq!(drained.wet_length_m, 0.0);
        assert_eq!(drained.end, CircularMeniscusEnd::Inlet);
        near(drained.elapsed_time_s, filled.elapsed_time_s);
        near(drained.transferred_volume_m3, -filled.transferred_volume_m3);
    }

    #[test]
    fn washburn_zero_drive_time_and_empty_nonwetting_do_no_work() {
        let input = CircularMeniscusInput {
            wet_length_m: 0.01,
            ..meniscus()
        };
        for changed in [
            CircularMeniscusInput {
                duration_s: 0.0,
                ..input
            },
            CircularMeniscusInput {
                surface_tension_n_m: 0.0,
                ..input
            },
            CircularMeniscusInput {
                contact_angle_cosine: 0.0,
                ..input
            },
            CircularMeniscusInput {
                applied_pressure_pa: -72.0,
                ..input
            },
        ] {
            let step = step_newtonian_circular_meniscus(changed).unwrap();
            assert_eq!(step.wet_length_m, changed.wet_length_m);
            assert_eq!(step.transferred_volume_m3, 0.0);
            assert_eq!(step.pressure_work_j, 0.0);
            assert_eq!(step.wall_surface_energy_change_j, 0.0);
            assert_eq!(step.viscous_dissipation_j, 0.0);
        }
        let empty = step_newtonian_circular_meniscus(CircularMeniscusInput {
            contact_angle_cosine: -0.5,
            ..meniscus()
        })
        .unwrap();
        assert_eq!(empty.end, CircularMeniscusEnd::Inlet);
        assert_eq!(empty.elapsed_time_s, 0.0);
        assert_eq!(empty.transferred_volume_m3, 0.0);
    }

    #[test]
    fn washburn_pressure_only_recovers_integrated_poiseuille() {
        let input = CircularMeniscusInput {
            surface_tension_n_m: 0.0,
            applied_pressure_pa: 100.0,
            wet_length_m: 0.02,
            ..meniscus()
        };
        let step = step_newtonian_circular_meniscus(input).unwrap();
        let midpoint_law = step_newtonian_circular_capillary(CircularCapillaryInput {
            dynamic_viscosity_pa_s: input.dynamic_viscosity_pa_s,
            radius_m: input.radius_m,
            length_m: (input.wet_length_m + step.wet_length_m) / 2.0,
            volumetric_flow_m3_s: step.transferred_volume_m3 / step.elapsed_time_s,
        })
        .unwrap();
        near(midpoint_law.pressure_drop_pa, input.applied_pressure_pa);
        assert_eq!(step.wall_surface_energy_change_j, 0.0);
        near(step.pressure_work_j, step.viscous_dissipation_j);
    }

    #[test]
    fn washburn_refuses_invalid_input_without_publishing_a_step() {
        let input = meniscus();
        for (changed, field) in [
            (
                CircularMeniscusInput {
                    contact_angle_cosine: 1.01,
                    ..input
                },
                "contact_angle_cosine",
            ),
            (
                CircularMeniscusInput {
                    surface_tension_n_m: -0.1,
                    ..input
                },
                "surface_tension_n_m",
            ),
            (
                CircularMeniscusInput {
                    duration_s: -1.0,
                    ..input
                },
                "duration_s",
            ),
            (
                CircularMeniscusInput {
                    wet_length_m: 2.0 * input.tube_length_m,
                    ..input
                },
                "wet_length_m",
            ),
        ] {
            assert_eq!(
                step_newtonian_circular_meniscus(changed),
                Err(CircularCapillaryError::OutsideDomain { field })
            );
        }
        assert_eq!(
            step_newtonian_circular_meniscus(CircularMeniscusInput {
                radius_m: 0.0,
                ..input
            }),
            Err(CircularCapillaryError::NonPositive { field: "radius_m" })
        );
        assert_eq!(
            step_newtonian_circular_meniscus(CircularMeniscusInput {
                applied_pressure_pa: f64::NAN,
                ..input
            }),
            Err(CircularCapillaryError::NonFinite {
                field: "applied_pressure_pa"
            })
        );
        assert_eq!(
            step_newtonian_circular_meniscus(CircularMeniscusInput {
                radius_m: 1e200,
                ..input
            }),
            Err(CircularCapillaryError::NonFiniteResult)
        );
    }

    #[test]
    fn matches_the_closed_form_and_quartic_radius_scaling() {
        let base = step_newtonian_circular_capillary(CircularCapillaryInput {
            dynamic_viscosity_pa_s: 280.0,
            length_m: 1.6e-3,
            radius_m: 0.2e-3,
            volumetric_flow_m3_s: 4.05e-9,
        })
        .expect("admitted capillary");
        let wide = step_newtonian_circular_capillary(CircularCapillaryInput {
            radius_m: 0.4e-3,
            ..CircularCapillaryInput {
                dynamic_viscosity_pa_s: 280.0,
                length_m: 1.6e-3,
                radius_m: 0.2e-3,
                volumetric_flow_m3_s: 4.05e-9,
            }
        })
        .expect("wider admitted capillary");

        assert!((base.pressure_drop_pa - 2_887_707.287_459_349).abs() < 1.0e-6);
        assert!((base.wall_shear_rate_per_s - 644.577_519_522_176).abs() < 1.0e-9);
        assert!((base.pressure_drop_pa / wide.pressure_drop_pa - 16.0).abs() < 1.0e-12);
        assert!((base.hydraulic_power_w - base.pressure_drop_pa * 4.05e-9).abs() < 1.0e-15);
    }

    #[test]
    fn zero_flow_is_a_valid_stationary_state() {
        let result = step_newtonian_circular_capillary(CircularCapillaryInput {
            dynamic_viscosity_pa_s: 1.0,
            length_m: 1.0,
            radius_m: 0.1,
            volumetric_flow_m3_s: 0.0,
        })
        .expect("zero-flow state");
        assert_eq!(result.pressure_drop_pa, 0.0);
        assert_eq!(result.wall_shear_rate_per_s, 0.0);
        assert_eq!(result.hydraulic_power_w, 0.0);
    }

    #[test]
    fn refuses_invalid_physical_inputs() {
        let input = CircularCapillaryInput {
            dynamic_viscosity_pa_s: 1.0,
            length_m: 1.0,
            radius_m: 0.1,
            volumetric_flow_m3_s: 1.0,
        };
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                radius_m: 0.0,
                ..input
            }),
            Err(CircularCapillaryError::NonPositive { field: "radius_m" })
        );
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                volumetric_flow_m3_s: -1.0,
                ..input
            }),
            Err(CircularCapillaryError::NegativeFlow)
        );
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                length_m: f64::NAN,
                ..input
            }),
            Err(CircularCapillaryError::NonFinite { field: "length_m" })
        );
    }
}
