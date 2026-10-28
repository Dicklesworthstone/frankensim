//! Reduced, deterministic circular-capillary flow screening.
//!
//! This module owns the Hagen–Poiseuille relation for a declared Newtonian,
//! incompressible fluid in a straight, rigid, circular capillary with
//! fully-developed laminar flow and no slip. It deliberately does not turn
//! that screening law into a polymer-rheology claim: entrance contraction,
//! shear thinning, viscoelasticity, wall slip and heating are outside this rung.
//! A separate horizontal Lucas–Washburn meniscus model admits capillary
//! imbibition and drainage with an explicit equilibrium contact angle.

use core::f64::consts::PI;

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
