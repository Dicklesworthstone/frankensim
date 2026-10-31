//! Rate-dependent wetting and static hysteresis in the inclined two-fluid tube.
//!
//! A linear contact-line friction `zeta` adds the constant series resistance
//! `zeta*R/4` to `mu_l*x + mu_d*(L-x)`. Advancing and receding thresholds select
//! an affine-force branch, so the original analytic gravity travel integral and
//! positive phase-loss moments still apply. The equilibrium Young angle owns
//! wall free energy; threshold pinning and line friction dissipate separately.
//!
//! This is an explicitly bounded, quasisteady constitutive law. It does not
//! resolve inertial startup, hydrodynamic contact-line slip, films, temperature
//! evolution or a molecular-kinetic/Cox--Voinov law. The line-power identity is
//! `D_line/dt = 2*pi*R*zeta*u²`; see Fricke et al., Physica D (2023),
//! <https://doi.org/10.1016/j.physd.2023.133895>, equations (11), (23)--(25).

use super::super::{CircularCapillaryError, CircularMeniscusEnd, CircularMeniscusStep};
use super::TwoFluidCircularMeniscusStep;
use super::gravity::{
    InclinedTwoFluidCircularMeniscusInput, InclinedTwoFluidCircularMeniscusStep, Law, phase_losses,
    validate_input,
};
use core::f64::consts::PI;

/// The force-selected contact-line branch during the accepted interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactLineMotion {
    /// The applied force lies inside the static hysteresis interval.
    Pinned,
    /// The invading liquid advances against its advancing threshold.
    Advancing,
    /// The invading liquid recedes against its receding threshold.
    Receding,
}

/// Explicit linear-friction/hysteresis closure for the original tube geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineCircularMeniscusInput {
    /// Original phase, geometry, reservoir-pressure and gravity inputs. Its
    /// contact-angle cosine is the equilibrium Young value, not a moving angle.
    pub meniscus: InclinedTwoFluidCircularMeniscusInput,
    /// Nonnegative linear contact-line friction coefficient [Pa s].
    pub contact_line_friction_pa_s: f64,
    /// Static advancing angle [rad], at least the equilibrium Young angle.
    pub advancing_angle_radians: f64,
    /// Static receding angle [rad], at most the equilibrium Young angle.
    pub receding_angle_radians: f64,
    /// Strictly positive admitted magnitude of line speed [m/s].
    pub max_speed_m_s: f64,
    /// Admitted absolute apparent-angle departure from the Young angle [rad].
    pub max_angle_departure_radians: f64,
}

/// Accepted motion and independently resolved mechanical loss channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineCircularMeniscusStep {
    /// Original gravity/pressure/surface accounts. The nested capillary pressure
    /// is the final apparent meniscus pressure, while wall free energy uses the
    /// equilibrium angle. The nested viscous loss remains the two bulk phases'
    /// sum; the complete mechanical balance uses `total_dissipation_j` below.
    pub meniscus: InclinedTwoFluidCircularMeniscusStep,
    /// Integrated positive contact-line friction loss [J].
    pub contact_line_viscous_dissipation_j: f64,
    /// Integrated positive static advancing/receding threshold loss [J].
    pub pinning_dissipation_j: f64,
    /// Sum of both bulk losses, contact-line friction and pinning [J].
    pub total_dissipation_j: f64,
    /// Constitutive branch selected by the initial force.
    pub motion: ContactLineMotion,
    /// Signed instantaneous speed at the beginning of the accepted interval.
    pub initial_speed_m_s: f64,
    /// Signed instantaneous speed at the accepted endpoint.
    pub final_speed_m_s: f64,
    /// Apparent angle through the invading liquid at the initial state [rad].
    pub initial_apparent_angle_radians: f64,
    /// Apparent angle through the invading liquid at the accepted state [rad].
    pub final_apparent_angle_radians: f64,
}

/// Advance a two-fluid capillary with linear contact-line friction and pinning.
///
/// Moving branches satisfy
/// `8*(mu_l*x + mu_d*(L-x) + zeta*R/4)*u = R²*F_branch(x)`.
/// The branch angle is `theta_A` for `u>0` and `theta_R` for `u<0`, with
/// `cos(theta_apparent) = cos(theta_branch) - zeta*u/sigma`.
/// Inside the static interval, the apparent angle balances the pressure and
/// hydrostatic head at zero speed. No work or dissipation is invented there.
///
/// The complete balance is
/// `W_pressure - delta_E_Young_wall - delta_E_gravity - E_outlet_potential`
/// `= D_liquid + D_displaced + D_line + D_pinning`.
/// The original gravity solver owns branch travel and endpoint events. It
/// never crosses a stable branch equilibrium. Speed is a ratio of affine
/// functions with positive resistance, so checking both interval endpoints
/// bounds speed and apparent-angle departure throughout the accepted motion.
///
/// # Errors
/// Refuses original tube-domain errors, nonpositive tension, malformed wetting
/// limits, impossible apparent angles, or any accepted-interval source-limit
/// excursion. It returns no partial step when a limit fails. The finite speed
/// bound excludes the dry, zero-resistance singularity when both downstream
/// viscosity and line friction are zero.
pub fn step_newtonian_contact_line_circular_meniscus(
    input: ContactLineCircularMeniscusInput,
) -> Result<ContactLineCircularMeniscusStep, CircularCapillaryError> {
    validate_input(input.meniscus)?;
    let m = input.meniscus.meniscus.meniscus;
    let young = m.contact_angle_cosine.acos();
    let (advancing, receding) = validate_contact_line(input)?;
    let sigma = m.surface_tension_n_m;
    let g = input.meniscus.opposing_gravity_m_s2;
    let rho = input.meniscus.density_kg_m3;
    let rho_d = input.meniscus.displaced_density_kg_m3;
    let mu_d = input.meniscus.meniscus.displaced_dynamic_viscosity_pa_s;
    let zeta = input.contact_line_friction_pa_s;
    let hydrostatic = |x: f64| g * rho.mul_add(x, rho_d * (m.tube_length_m - x));
    let area = PI * m.radius_m * m.radius_m;
    let capillary_scale = (sigma / m.radius_m) * 2.0;
    let young_pressure = capillary_scale * m.contact_angle_cosine;
    let force_young = (m.applied_pressure_pa - hydrostatic(m.wet_length_m)) + young_pressure;
    let advancing_threshold = capillary_scale * (m.contact_angle_cosine - advancing);
    let receding_threshold = capillary_scale * (receding - m.contact_angle_cosine);
    let (motion, branch_cosine, pinning_pressure, force0) = if force_young > advancing_threshold {
        (
            ContactLineMotion::Advancing,
            advancing,
            advancing_threshold,
            force_young - advancing_threshold,
        )
    } else if force_young < -receding_threshold {
        (
            ContactLineMotion::Receding,
            receding,
            -receding_threshold,
            force_young + receding_threshold,
        )
    } else {
        (ContactLineMotion::Pinned, m.contact_angle_cosine, 0.0, 0.0)
    };
    let line_viscosity_scale = (zeta * 0.25) * (m.radius_m / m.tube_length_m);
    let mu_scale = m.dynamic_viscosity_pa_s.max(mu_d).max(line_viscosity_scale);
    let law = Law {
        initial: m.wet_length_m,
        length: m.tube_length_m,
        liquid: m.dynamic_viscosity_pa_s / mu_scale,
        displaced: mu_d / mu_scale,
        contact_line: (zeta / mu_scale) * (m.radius_m / 4.0),
        mobility: (m.radius_m / mu_scale) * (m.radius_m / 8.0),
        gravity_slope: (rho - rho_d) * g,
        force0,
    };
    let outlet_head = (rho_d * g) * m.tube_length_m;
    if [
        area,
        capillary_scale,
        young_pressure,
        force_young,
        advancing_threshold,
        receding_threshold,
        line_viscosity_scale,
        law.contact_line,
        law.mobility,
        law.gravity_slope,
        law.force0,
        outlet_head,
        hydrostatic(0.0),
        hydrostatic(law.length),
    ]
    .into_iter()
    .any(|v| !v.is_finite())
        || area <= 0.0
        || capillary_scale <= 0.0
        || law.mobility <= 0.0
        || law.liquid <= 0.0
        || (mu_d > 0.0 && law.displaced <= 0.0)
        || (zeta > 0.0 && law.contact_line <= 0.0)
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    let pinned_angle = if motion == ContactLineMotion::Pinned {
        // The known static interval only corrects roundoff in the reconstruction.
        (m.contact_angle_cosine - force_young / capillary_scale)
            .clamp(advancing, receding)
            .acos()
    } else {
        young
    };
    let diagnostics = |x: f64| -> Result<(f64, f64), CircularCapillaryError> {
        let (speed, angle) = if motion == ContactLineMotion::Pinned {
            (0.0, pinned_angle)
        } else {
            let resistance = law.resistance(x);
            let speed = law.mobility * law.force(x - law.initial) / resistance;
            if resistance <= 0.0 || !speed.is_finite() {
                return Err(CircularCapillaryError::OutsideDomain {
                    field: "contact_line_speed_m_s",
                });
            }
            let cosine = branch_cosine - (zeta / sigma) * speed;
            if !(-1.0..=1.0).contains(&cosine) {
                return Err(CircularCapillaryError::OutsideDomain {
                    field: "apparent_contact_angle",
                });
            }
            (speed, cosine.acos())
        };
        if speed.abs() > input.max_speed_m_s {
            return Err(CircularCapillaryError::OutsideDomain {
                field: "contact_line_speed_m_s",
            });
        }
        if (angle - young).abs() > input.max_angle_departure_radians {
            return Err(CircularCapillaryError::OutsideDomain {
                field: "contact_angle_departure_radians",
            });
        }
        Ok((speed, angle))
    };
    let (initial_speed, initial_angle) = diagnostics(law.initial)?;
    let (length, elapsed, end) = if motion == ContactLineMotion::Pinned || m.duration_s == 0.0 {
        (law.initial, m.duration_s, CircularMeniscusEnd::TimeLimit)
    } else {
        law.advance(m.duration_s)?
    };
    let (final_speed, final_angle) = diagnostics(length)?;
    let delta = length - law.initial;
    let volume = area * delta;
    let dynamic_loss = (0.5 * law.force0 + 0.5 * law.force(delta)) * volume;
    let (liquid_loss, displaced_loss, line_loss) = phase_losses(law, length, dynamic_loss)?;
    let bulk_loss = liquid_loss + displaced_loss;
    let pinning_loss = pinning_pressure * volume;
    let total_loss = (bulk_loss + line_loss) + pinning_loss;
    let pressure_work = m.applied_pressure_pa * volume;
    let wall_energy = -young_pressure * volume;
    let gravity_storage = (law.gravity_slope * (law.initial + 0.5 * delta)) * volume;
    let gravity_export = outlet_head * volume;
    let final_head = hydrostatic(length);
    let final_capillary = capillary_scale * final_angle.cos();
    if [
        length,
        elapsed,
        volume,
        dynamic_loss,
        liquid_loss,
        displaced_loss,
        line_loss,
        bulk_loss,
        pinning_loss,
        total_loss,
        pressure_work,
        wall_energy,
        gravity_storage,
        gravity_export,
        final_head,
        final_capillary,
    ]
    .into_iter()
    .any(|v| !v.is_finite())
        || !(0.0..=law.length).contains(&length)
        || elapsed < 0.0
        || [
            dynamic_loss,
            liquid_loss,
            displaced_loss,
            line_loss,
            pinning_loss,
        ]
        .into_iter()
        .any(|v| v < 0.0)
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    Ok(ContactLineCircularMeniscusStep {
        meniscus: InclinedTwoFluidCircularMeniscusStep {
            meniscus: TwoFluidCircularMeniscusStep {
                meniscus: CircularMeniscusStep {
                    wet_length_m: length,
                    elapsed_time_s: elapsed,
                    end,
                    transferred_volume_m3: volume,
                    capillary_pressure_pa: final_capillary,
                    pressure_work_j: pressure_work,
                    wall_surface_energy_change_j: wall_energy,
                    viscous_dissipation_j: bulk_loss,
                },
                liquid_viscous_dissipation_j: liquid_loss,
                displaced_viscous_dissipation_j: displaced_loss,
            },
            gravitational_storage_change_j: gravity_storage,
            outlet_potential_energy_j: gravity_export,
            hydrostatic_pressure_pa: final_head,
        },
        contact_line_viscous_dissipation_j: line_loss,
        pinning_dissipation_j: pinning_loss,
        total_dissipation_j: total_loss,
        motion,
        initial_speed_m_s: initial_speed,
        final_speed_m_s: final_speed,
        initial_apparent_angle_radians: initial_angle,
        final_apparent_angle_radians: final_angle,
    })
}

fn validate_contact_line(
    input: ContactLineCircularMeniscusInput,
) -> Result<(f64, f64), CircularCapillaryError> {
    for (field, value) in [
        (
            "contact_line_friction_pa_s",
            input.contact_line_friction_pa_s,
        ),
        ("advancing_angle_radians", input.advancing_angle_radians),
        ("receding_angle_radians", input.receding_angle_radians),
        ("max_speed_m_s", input.max_speed_m_s),
        (
            "max_angle_departure_radians",
            input.max_angle_departure_radians,
        ),
    ] {
        if !value.is_finite() {
            return Err(CircularCapillaryError::NonFinite { field });
        }
    }
    for (field, value) in [
        (
            "surface_tension_n_m",
            input.meniscus.meniscus.meniscus.surface_tension_n_m,
        ),
        ("max_speed_m_s", input.max_speed_m_s),
    ] {
        if value <= 0.0 {
            return Err(CircularCapillaryError::NonPositive { field });
        }
    }
    for (field, admitted) in [
        (
            "contact_line_friction_pa_s",
            input.contact_line_friction_pa_s >= 0.0,
        ),
        (
            "advancing_angle_radians",
            (0.0..=PI).contains(&input.advancing_angle_radians),
        ),
        (
            "receding_angle_radians",
            (0.0..=PI).contains(&input.receding_angle_radians),
        ),
        (
            "max_angle_departure_radians",
            (0.0..=PI).contains(&input.max_angle_departure_radians),
        ),
    ] {
        if !admitted {
            return Err(CircularCapillaryError::OutsideDomain { field });
        }
    }
    let young = input.meniscus.meniscus.meniscus.contact_angle_cosine;
    let advancing = input.advancing_angle_radians.cos();
    let receding = input.receding_angle_radians.cos();
    // Angle/cosine round trips may differ by one ulp at a zero-width interval.
    let tolerance = 8.0 * f64::EPSILON;
    if advancing > young + tolerance
        || receding < young - tolerance
        || input.receding_angle_radians > input.advancing_angle_radians
    {
        return Err(CircularCapillaryError::OutsideDomain {
            field: "contact_angle_hysteresis",
        });
    }
    let advancing = if (advancing - young).abs() <= tolerance {
        young
    } else {
        advancing
    };
    let receding = if (receding - young).abs() <= tolerance {
        young
    } else {
        receding
    };
    Ok((advancing, receding))
}

#[cfg(test)]
mod tests {
    use super::super::super::CircularMeniscusInput;
    use super::super::TwoFluidCircularMeniscusInput;
    use super::super::gravity::step_newtonian_inclined_two_fluid_circular_meniscus;
    use super::*;

    fn input() -> ContactLineCircularMeniscusInput {
        ContactLineCircularMeniscusInput {
            meniscus: InclinedTwoFluidCircularMeniscusInput {
                meniscus: TwoFluidCircularMeniscusInput {
                    meniscus: CircularMeniscusInput {
                        dynamic_viscosity_pa_s: 0.1,
                        radius_m: 0.001,
                        tube_length_m: 0.02,
                        wet_length_m: 0.004,
                        surface_tension_n_m: 0.072,
                        contact_angle_cosine: 0.5,
                        applied_pressure_pa: 160.0,
                        duration_s: 0.25,
                    },
                    displaced_dynamic_viscosity_pa_s: 0.2,
                },
                density_kg_m3: 1000.0,
                displaced_density_kg_m3: 800.0,
                opposing_gravity_m_s2: 9.81,
            },
            contact_line_friction_pa_s: 0.6,
            advancing_angle_radians: 70.0_f64.to_radians(),
            receding_angle_radians: 50.0_f64.to_radians(),
            max_speed_m_s: 1.0,
            max_angle_departure_radians: 0.7,
        }
    }

    fn near(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-20),
            "{a} vs {b}"
        );
    }

    fn balance(step: ContactLineCircularMeniscusStep) {
        let g = step.meniscus;
        let m = g.meniscus.meniscus;
        near(
            m.pressure_work_j
                - m.wall_surface_energy_change_j
                - g.gravitational_storage_change_j
                - g.outlet_potential_energy_j,
            step.total_dissipation_j,
            2e-12,
        );
        near(
            g.meniscus.liquid_viscous_dissipation_j + g.meniscus.displaced_viscous_dissipation_j,
            m.viscous_dissipation_j,
            2e-15,
        );
        near(
            m.viscous_dissipation_j
                + step.contact_line_viscous_dissipation_j
                + step.pinning_dissipation_j,
            step.total_dissipation_j,
            2e-15,
        );
        for loss in [
            g.meniscus.liquid_viscous_dissipation_j,
            g.meniscus.displaced_viscous_dissipation_j,
            step.contact_line_viscous_dissipation_j,
            step.pinning_dissipation_j,
        ] {
            assert!(loss >= 0.0);
        }
    }

    #[test]
    fn contact_line_g3_zero_friction_hysteresis_recovers_original_tube() {
        for (gravity, mu_d, rho_d) in [(0.0, 0.2, 800.0), (9.81, 0.2, 800.0), (9.81, 0.0, 0.0)] {
            let mut i = input();
            i.contact_line_friction_pa_s = 0.0;
            i.advancing_angle_radians = i.meniscus.meniscus.meniscus.contact_angle_cosine.acos();
            i.receding_angle_radians = i.advancing_angle_radians;
            i.meniscus.opposing_gravity_m_s2 = gravity;
            i.meniscus.displaced_density_kg_m3 = rho_d;
            i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = mu_d;
            let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
            let original = step_newtonian_inclined_two_fluid_circular_meniscus(i.meniscus).unwrap();
            near(
                s.meniscus.meniscus.meniscus.wet_length_m,
                original.meniscus.meniscus.wet_length_m,
                2e-13,
            );
            near(
                s.total_dissipation_j,
                original.meniscus.meniscus.viscous_dissipation_j,
                2e-13,
            );
            near(
                s.meniscus.meniscus.liquid_viscous_dissipation_j,
                original.meniscus.liquid_viscous_dissipation_j,
                2e-13,
            );
            assert_eq!(s.contact_line_viscous_dissipation_j, 0.0);
            assert_eq!(s.pinning_dissipation_j, 0.0);
            balance(s);
        }
    }

    #[test]
    fn contact_line_g1_friction_regularizes_dry_entry_with_exact_travel_and_loss() {
        let mut i = input();
        i.meniscus.opposing_gravity_m_s2 = 0.0;
        i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 0.0;
        i.contact_line_friction_pa_s = 1.0;
        i.advancing_angle_radians = 0.5_f64.acos();
        i.receding_angle_radians = i.advancing_angle_radians;
        i.meniscus.meniscus.meniscus.wet_length_m = 0.0;
        i.meniscus.meniscus.meniscus.applied_pressure_pa = 0.0;
        i.meniscus.meniscus.meniscus.duration_s = 0.01;
        let m = i.meniscus.meniscus.meniscus;
        let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        let pressure = 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m;
        let line = i.contact_line_friction_pa_s * m.radius_m / 4.0;
        let advance = m.radius_m.powi(2) * pressure * m.duration_s / 8.0;
        let expected = 2.0 * advance
            / (line + (line * line + 2.0 * m.dynamic_viscosity_pa_s * advance).sqrt());
        near(s.meniscus.meniscus.meniscus.wet_length_m, expected, 2e-13);
        near(
            s.initial_speed_m_s,
            m.radius_m.powi(2) * pressure / (8.0 * line),
            2e-15,
        );
        let expected_line_loss = PI * m.radius_m.powi(2) * pressure * line
            / m.dynamic_viscosity_pa_s
            * (m.dynamic_viscosity_pa_s * expected / line).ln_1p();
        near(
            s.contact_line_viscous_dissipation_j,
            expected_line_loss,
            2e-13,
        );
        near(s.initial_apparent_angle_radians, PI / 2.0, 2e-15);
        assert!(s.final_speed_m_s < s.initial_speed_m_s);
        assert_eq!(s.meniscus.meniscus.displaced_viscous_dissipation_j, 0.0);
        balance(s);
    }

    #[test]
    fn contact_line_g1_pinning_and_both_direction_thresholds_keep_young_wall_energy() {
        let mut i = input();
        i.meniscus.opposing_gravity_m_s2 = 0.0;
        let m = i.meniscus.meniscus.meniscus;
        let scale = 2.0 * m.surface_tension_n_m / m.radius_m;
        let a = scale * (m.contact_angle_cosine - i.advancing_angle_radians.cos());
        let r = scale * (i.receding_angle_radians.cos() - m.contact_angle_cosine);
        for force in [0.0, 0.5 * a, -0.5 * r] {
            i.meniscus.meniscus.meniscus.applied_pressure_pa =
                force - scale * m.contact_angle_cosine;
            let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
            assert_eq!(s.motion, ContactLineMotion::Pinned);
            assert_eq!(s.meniscus.meniscus.meniscus.wet_length_m, m.wet_length_m);
            assert_eq!(s.total_dissipation_j, 0.0);
            assert_eq!(s.initial_speed_m_s, 0.0);
            near(
                s.meniscus.meniscus.meniscus.capillary_pressure_pa,
                -i.meniscus.meniscus.meniscus.applied_pressure_pa,
                2e-15,
            );
        }
        for (force, motion, threshold) in [
            (a + 10.0, ContactLineMotion::Advancing, a),
            (-r - 10.0, ContactLineMotion::Receding, -r),
        ] {
            i.meniscus.meniscus.meniscus.applied_pressure_pa =
                force - scale * m.contact_angle_cosine;
            let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
            assert_eq!(s.motion, motion);
            let result = s.meniscus.meniscus.meniscus;
            near(
                s.pinning_dissipation_j,
                threshold * result.transferred_volume_m3,
                2e-13,
            );
            near(
                result.wall_surface_energy_change_j,
                -scale * m.contact_angle_cosine * result.transferred_volume_m3,
                2e-15,
            );
            assert!(s.pinning_dissipation_j > 0.0);
            balance(s);
        }
    }

    #[test]
    fn contact_line_g1_all_three_viscous_channels_match_spatial_quadrature() {
        for pressure in [160.0, 30.0] {
            let mut i = input();
            i.meniscus.meniscus.meniscus.applied_pressure_pa = pressure;
            i.meniscus.meniscus.meniscus.duration_s = 1.0;
            let m = i.meniscus.meniscus.meniscus;
            let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
            let end = s.meniscus.meniscus.meniscus.wet_length_m;
            let angle = if end > m.wet_length_m {
                i.advancing_angle_radians
            } else {
                i.receding_angle_radians
            };
            let mut integrals = [0.0; 3];
            let count = 2048;
            for k in 0..=count {
                let x = m.wet_length_m + (end - m.wet_length_m) * f64::from(k) / f64::from(count);
                let resistance = [
                    m.dynamic_viscosity_pa_s * x,
                    i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s * (m.tube_length_m - x),
                    i.contact_line_friction_pa_s * m.radius_m / 4.0,
                ];
                let total = resistance.iter().sum::<f64>();
                let force = pressure + 2.0 * m.surface_tension_n_m * angle.cos() / m.radius_m
                    - i.meniscus.opposing_gravity_m_s2
                        * (i.meniscus.density_kg_m3 * x
                            + i.meniscus.displaced_density_kg_m3 * (m.tube_length_m - x));
                let weight = if k == 0 || k == count {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for channel in 0..3 {
                    integrals[channel] += weight * force * resistance[channel] / total;
                }
            }
            let factor =
                PI * m.radius_m.powi(2) * (end - m.wet_length_m) / (3.0 * f64::from(count));
            for (actual, integral) in [
                s.meniscus.meniscus.liquid_viscous_dissipation_j,
                s.meniscus.meniscus.displaced_viscous_dissipation_j,
                s.contact_line_viscous_dissipation_j,
            ]
            .into_iter()
            .zip(integrals)
            {
                near(actual, factor * integral, 2e-12);
            }
            balance(s);
        }
    }

    #[test]
    fn contact_line_g3_gravity_exponential_composition_and_equilibrium_barrier() {
        let mut i = input();
        i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 0.1;
        i.meniscus.meniscus.meniscus.applied_pressure_pa = 130.0;
        i.meniscus.meniscus.meniscus.duration_s = 5.0;
        let m = i.meniscus.meniscus.meniscus;
        let slope = (i.meniscus.density_kg_m3 - i.meniscus.displaced_density_kg_m3)
            * i.meniscus.opposing_gravity_m_s2;
        let inlet = m.applied_pressure_pa
            + 2.0 * m.surface_tension_n_m * i.advancing_angle_radians.cos() / m.radius_m
            - i.meniscus.displaced_density_kg_m3
                * i.meniscus.opposing_gravity_m_s2
                * m.tube_length_m;
        let equilibrium = inlet / slope;
        let resistance = m.dynamic_viscosity_pa_s * m.tube_length_m
            + i.contact_line_friction_pa_s * m.radius_m / 4.0;
        let decay = m.radius_m.powi(2) * slope / (8.0 * resistance);
        let full = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        near(
            full.meniscus.meniscus.meniscus.wet_length_m,
            equilibrium + (m.wet_length_m - equilibrium) * (-decay * m.duration_s).exp(),
            2e-13,
        );
        i.meniscus.meniscus.meniscus.duration_s = 0.5 * m.duration_s;
        let first = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        i.meniscus.meniscus.meniscus.wet_length_m = first.meniscus.meniscus.meniscus.wet_length_m;
        let second = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        near(
            second.meniscus.meniscus.meniscus.wet_length_m,
            full.meniscus.meniscus.meniscus.wet_length_m,
            2e-13,
        );
        near(
            first.total_dissipation_j + second.total_dissipation_j,
            full.total_dissipation_j,
            2e-13,
        );
        near(
            first.contact_line_viscous_dissipation_j + second.contact_line_viscous_dissipation_j,
            full.contact_line_viscous_dissipation_j,
            2e-13,
        );
        i.meniscus.meniscus.meniscus.duration_s = 10000.0;
        let settled = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        near(
            settled.meniscus.meniscus.meniscus.wet_length_m,
            equilibrium,
            2e-13,
        );
        assert!(settled.final_speed_m_s >= 0.0);
        assert_eq!(
            settled.meniscus.meniscus.meniscus.end,
            CircularMeniscusEnd::TimeLimit
        );
        balance(full);
        balance(settled);
    }

    #[test]
    fn contact_line_g1_both_end_events_retain_actual_time_and_line_power() {
        for pressure in [180.0, -180.0] {
            let mut i = input();
            i.meniscus.opposing_gravity_m_s2 = 0.0;
            i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 0.1;
            i.meniscus.meniscus.meniscus.applied_pressure_pa = pressure;
            i.meniscus.meniscus.meniscus.duration_s = 10.0;
            let m = i.meniscus.meniscus.meniscus;
            let s = step_newtonian_contact_line_circular_meniscus(i).unwrap();
            let result = s.meniscus.meniscus.meniscus;
            let (boundary, event) = if pressure > 0.0 {
                (m.tube_length_m, CircularMeniscusEnd::Outlet)
            } else {
                (0.0, CircularMeniscusEnd::Inlet)
            };
            assert_eq!(result.end, event);
            assert_eq!(result.wet_length_m, boundary);
            near(
                result.elapsed_time_s,
                (boundary - m.wet_length_m) / s.initial_speed_m_s,
                2e-14,
            );
            near(
                s.contact_line_viscous_dissipation_j,
                2.0 * PI
                    * m.radius_m
                    * i.contact_line_friction_pa_s
                    * s.initial_speed_m_s.powi(2)
                    * result.elapsed_time_s,
                2e-14,
            );
            assert!(result.elapsed_time_s < m.duration_s);
            balance(s);
        }
    }

    #[test]
    fn contact_line_g0_endpoint_limits_and_original_input_domain_refuse_atomically() {
        let mut i = input();
        i.meniscus.opposing_gravity_m_s2 = 0.0;
        i.meniscus.meniscus.meniscus.duration_s = 1.0;
        let full = step_newtonian_contact_line_circular_meniscus(i).unwrap();
        assert!(full.final_speed_m_s > full.initial_speed_m_s);
        let mut limited = i;
        limited.max_speed_m_s = 0.5 * (full.initial_speed_m_s + full.final_speed_m_s);
        assert_eq!(
            step_newtonian_contact_line_circular_meniscus(limited),
            Err(CircularCapillaryError::OutsideDomain {
                field: "contact_line_speed_m_s"
            })
        );
        limited.meniscus.meniscus.meniscus.duration_s = 0.0;
        assert!(step_newtonian_contact_line_circular_meniscus(limited).is_ok());
        limited = i;
        limited.max_angle_departure_radians = 0.5
            * (full.initial_apparent_angle_radians + full.final_apparent_angle_radians)
            - i.meniscus.meniscus.meniscus.contact_angle_cosine.acos();
        assert_eq!(
            step_newtonian_contact_line_circular_meniscus(limited),
            Err(CircularCapillaryError::OutsideDomain {
                field: "contact_angle_departure_radians"
            })
        );
        limited.meniscus.meniscus.meniscus.duration_s = 0.0;
        assert!(step_newtonian_contact_line_circular_meniscus(limited).is_ok());
        limited = i;
        limited.contact_line_friction_pa_s = -1.0;
        assert!(step_newtonian_contact_line_circular_meniscus(limited).is_err());
        limited = i;
        limited.meniscus.density_kg_m3 = 0.0;
        assert_eq!(
            step_newtonian_contact_line_circular_meniscus(limited),
            Err(CircularCapillaryError::NonPositive {
                field: "density_kg_m3"
            })
        );
        limited = i;
        limited.advancing_angle_radians = limited.receding_angle_radians;
        assert!(step_newtonian_contact_line_circular_meniscus(limited).is_err());
    }
}
