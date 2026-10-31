//! Exact horizontal meniscus transport with both phases' viscous resistance.
//!
//! The series Poiseuille resistance is proportional to
//! `mu_liquid*x + mu_displaced*(L-x)`. Integrating that affine resistance
//! gives a quadratic travel coordinate, inverted without subtracting roots.
//! This is the straight-tube limit of Budaraju et al., Langmuir (2016),
//! <https://doi.org/10.1021/acs.langmuir.6b02788>. The two phases are immiscible,
//! Newtonian and incompressible, or the displaced gas is explicitly treated
//! as a uniform-property, negligible-compression approximation. Both ends
//! communicate with prescribed-pressure reservoirs. No trapped-gas pressure,
//! film deposition, inertia, entry loss or dynamic contact angle is inferred.

use super::{
    CircularCapillaryError, CircularMeniscusEnd, CircularMeniscusInput, CircularMeniscusStep,
    step_newtonian_circular_meniscus, validate_circular_meniscus,
};
use core::f64::consts::PI;

/// A translating meniscus with a viscous downstream phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidCircularMeniscusInput {
    /// Invading phase, tube geometry, current state and signed pressure drive.
    /// Surface tension and angle must describe this actual pair of phases
    /// and substrate. Applied pressure is inlet minus outlet pressure.
    pub meniscus: CircularMeniscusInput,
    /// Viscosity [Pa s] of the displaced phase occupying the remaining tube.
    /// Zero explicitly recovers the original negligible-gas-resistance law.
    pub displaced_dynamic_viscosity_pa_s: f64,
}

/// Accepted two-phase displacement with separately resolved viscous losses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidCircularMeniscusStep {
    /// State, signed exchanged volume, pressure work and surface-energy change.
    /// Its viscous dissipation is the sum of both phase losses below.
    pub meniscus: CircularMeniscusStep,
    /// Poiseuille dissipation inside the invading liquid column [J].
    pub liquid_viscous_dissipation_j: f64,
    /// Poiseuille dissipation inside the displaced fluid column [J].
    pub displaced_viscous_dissipation_j: f64,
}

/// Advance the two-fluid Washburn equation with constant coefficients.
///
/// The motion satisfies `8*(mu_l*x+mu_d*(L-x))*dx/dt = R²*drive`, with
/// `drive = applied_pressure + 2*sigma*cos(theta)/R`. Positive downstream
/// viscosity gives finite initial velocity at x=0, even in this inertia-free
/// model. Equal viscosities give constant speed; zero downstream viscosity
/// delegates exactly to the original single-phase law.
///
/// End events consume only their exact travel time. The translating meniscus
/// has constant area; the changing wall energy supplies capillary work.
/// The two viscous losses are analytic integrals of each phase's share of
/// the instantaneous resistance, not a split evaluated at one endpoint.
/// This bounded scalar kernel allocates no memory and mutates no caller state.
///
/// # Errors
/// Refuses malformed inputs and unrepresentable geometry, mobility or work.
pub fn step_newtonian_two_fluid_circular_meniscus(
    input: TwoFluidCircularMeniscusInput,
) -> Result<TwoFluidCircularMeniscusStep, CircularCapillaryError> {
    let m = input.meniscus;
    validate_circular_meniscus(m)?;
    let mu_d = input.displaced_dynamic_viscosity_pa_s;
    if !mu_d.is_finite() {
        return Err(CircularCapillaryError::NonFinite {
            field: "displaced_dynamic_viscosity_pa_s",
        });
    }
    if mu_d < 0.0 {
        return Err(CircularCapillaryError::OutsideDomain {
            field: "displaced_dynamic_viscosity_pa_s",
        });
    }
    if mu_d == 0.0 {
        let meniscus = step_newtonian_circular_meniscus(m)?;
        return Ok(TwoFluidCircularMeniscusStep {
            liquid_viscous_dissipation_j: meniscus.viscous_dissipation_j,
            displaced_viscous_dissipation_j: 0.0,
            meniscus,
        });
    }

    // Normalize both viscosity and position before the quadratic solve. Its
    // discriminant then uses numbers of order one even for microscopic tubes.
    let mu_scale = m.dynamic_viscosity_pa_s.max(mu_d);
    let liquid = m.dynamic_viscosity_pa_s / mu_scale;
    let displaced = mu_d / mu_scale;
    let area = PI * m.radius_m * m.radius_m;
    let aspect = m.radius_m / m.tube_length_m;
    let capillary = if m.contact_angle_cosine == 0.0 {
        0.0
    } else {
        (m.surface_tension_n_m / m.radius_m) * (2.0 * m.contact_angle_cosine)
    };
    let drive = m.applied_pressure_pa + capillary;
    let rate = (aspect * aspect) * (drive / mu_scale) / 8.0;
    let x0 = m.wet_length_m / m.tube_length_m;
    let resistance0 = liquid * x0 + displaced * (1.0 - x0);
    if !area.is_finite()
        || area <= 0.0
        || liquid <= 0.0
        || displaced <= 0.0
        || !capillary.is_finite()
        || !drive.is_finite()
        || !rate.is_finite()
        || (drive != 0.0 && rate == 0.0)
        || resistance0 <= 0.0
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    let (length, elapsed, end) = if rate == 0.0 || m.duration_s == 0.0 {
        (m.wet_length_m, m.duration_s, CircularMeniscusEnd::TimeLimit)
    } else {
        let (distance, boundary, endpoint_resistance, event) = if rate > 0.0 {
            (
                (m.tube_length_m - m.wet_length_m) / m.tube_length_m,
                m.tube_length_m,
                liquid,
                CircularMeniscusEnd::Outlet,
            )
        } else {
            (-x0, 0.0, displaced, CircularMeniscusEnd::Inlet)
        };
        let travel_time = distance * (0.5 * resistance0 + 0.5 * endpoint_resistance) / rate;
        if distance != 0.0 && travel_time == 0.0 {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        if travel_time <= m.duration_s {
            (boundary, travel_time, event)
        } else {
            let advance = rate * m.duration_s;
            let discriminant =
                (2.0 * (liquid - displaced)).mul_add(advance, resistance0 * resistance0);
            if discriminant < 0.0 || !discriminant.is_finite() {
                return Err(CircularCapillaryError::NonFiniteResult);
            }
            let delta = 2.0 * advance / (resistance0 + discriminant.sqrt());
            let length = delta.mul_add(m.tube_length_m, m.wet_length_m);
            if length == m.wet_length_m {
                return Err(CircularCapillaryError::NonFiniteResult);
            }
            (length, m.duration_s, CircularMeniscusEnd::TimeLimit)
        }
    };

    let delta = length - m.wet_length_m;
    let volume = area * delta;
    let pressure_work = m.applied_pressure_pa * volume;
    let wall_energy = -capillary * volume;
    let dissipation = drive * volume;
    let x1 = length / m.tube_length_m;
    let resistance1 = liquid * x1 + displaced * (1.0 - x1);
    let (w0, w1) = resistance_weights(resistance0, resistance1);
    let liquid_share = w0 * (liquid * x0 / resistance0) + w1 * (liquid * x1 / resistance1);
    let displaced_share =
        w0 * (displaced * (1.0 - x0) / resistance0) + w1 * (displaced * (1.0 - x1) / resistance1);
    // Evaluate the smaller loss directly, then retain exact aggregate closure
    // without subtracting a small loss from two large phase integrals.
    let (liquid_loss, displaced_loss) = if liquid_share <= displaced_share {
        let loss = dissipation * liquid_share;
        (loss, dissipation - loss)
    } else {
        let loss = dissipation * displaced_share;
        (dissipation - loss, loss)
    };
    if [
        length,
        elapsed,
        volume,
        pressure_work,
        wall_energy,
        dissipation,
        liquid_loss,
        displaced_loss,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || !(0.0..=m.tube_length_m).contains(&length)
        || elapsed < 0.0
        || dissipation < 0.0
        || liquid_loss < 0.0
        || displaced_loss < 0.0
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    Ok(TwoFluidCircularMeniscusStep {
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
        liquid_viscous_dissipation_j: liquid_loss,
        displaced_viscous_dissipation_j: displaced_loss,
    })
}

// Integrate the endpoint barycentric weights a*(1-t)/R(t), b*t/R(t).
// Positive weights sum to one. Interpolating each endpoint's phase resistance
// share with these weights gives its exact displacement-averaged share.
fn resistance_weights(a: f64, b: f64) -> (f64, f64) {
    let ratio = a.min(b) / a.max(b);
    let d = 1.0 - ratio;
    let small = if d < 0.125 {
        // (-d - log(1-d))/d², with no removable-singularity cancellation.
        let mut sum: f64 = 1.0 / 18.0;
        for k in (2..18).rev() {
            sum = sum.mul_add(d, 1.0 / f64::from(k));
        }
        ratio * sum
    } else if ratio == 0.0 {
        0.0
    } else {
        ratio * (-d - ratio.ln()) / (d * d)
    };
    if a <= b {
        (small, 1.0 - small)
    } else {
        (1.0 - small, small)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input() -> TwoFluidCircularMeniscusInput {
        TwoFluidCircularMeniscusInput {
            meniscus: CircularMeniscusInput {
                dynamic_viscosity_pa_s: 0.002,
                radius_m: 0.0002,
                tube_length_m: 0.1,
                wet_length_m: 0.0,
                surface_tension_n_m: 0.072,
                contact_angle_cosine: 0.5,
                applied_pressure_pa: 0.0,
                duration_s: 0.02,
            },
            displaced_dynamic_viscosity_pa_s: 0.001,
        }
    }
    fn near(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-30),
            "{a} != {b}"
        );
    }
    fn balance(step: TwoFluidCircularMeniscusStep) {
        let s = step.meniscus;
        near(
            s.pressure_work_j - s.wall_surface_energy_change_j,
            s.viscous_dissipation_j,
            2e-13,
        );
        near(
            step.liquid_viscous_dissipation_j + step.displaced_viscous_dissipation_j,
            s.viscous_dissipation_j,
            2e-15,
        );
        assert!(step.liquid_viscous_dissipation_j >= 0.0);
        assert!(step.displaced_viscous_dissipation_j >= 0.0);
    }

    #[test]
    fn two_fluid_g1_integrated_poiseuille_and_finite_startup() {
        let mut i = input();
        i.meniscus.duration_s = 1e-7;
        let s = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        let m = i.meniscus;
        let drive = 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m;
        let initial_velocity = m.radius_m.powi(2) * drive
            / (8.0 * i.displaced_dynamic_viscosity_pa_s * m.tube_length_m);
        near(
            s.meniscus.wet_length_m / m.duration_s,
            initial_velocity,
            1e-6,
        );
        for mu_d in [1e-5, 0.001, 0.002, 0.1] {
            i = input();
            i.displaced_dynamic_viscosity_pa_s = mu_d;
            let s = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
            let x = s.meniscus.wet_length_m;
            let integrated = mu_d * i.meniscus.tube_length_m * x
                + 0.5 * (i.meniscus.dynamic_viscosity_pa_s - mu_d) * x * x;
            near(
                integrated,
                i.meniscus.radius_m.powi(2) * drive * i.meniscus.duration_s / 8.0,
                2e-13,
            );
            balance(s);
        }
    }

    #[test]
    fn two_fluid_equal_viscosities_give_constant_speed_and_exact_loss_split() {
        let mut i = input();
        i.displaced_dynamic_viscosity_pa_s = i.meniscus.dynamic_viscosity_pa_s;
        i.meniscus.wet_length_m = 0.03;
        let s = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        let m = i.meniscus;
        let drive = s.meniscus.capillary_pressure_pa;
        let speed = m.radius_m.powi(2) * drive / (8.0 * m.dynamic_viscosity_pa_s * m.tube_length_m);
        near(
            s.meniscus.wet_length_m - m.wet_length_m,
            speed * m.duration_s,
            2e-13,
        );
        let mean_fraction = 0.5 * (m.wet_length_m + s.meniscus.wet_length_m) / m.tube_length_m;
        near(
            s.liquid_viscous_dissipation_j,
            mean_fraction * s.meniscus.viscous_dissipation_j,
            2e-13,
        );
        balance(s);
    }

    #[test]
    fn two_fluid_g3_zero_displaced_viscosity_is_original_law_exactly() {
        let mut i = input();
        i.displaced_dynamic_viscosity_pa_s = 0.0;
        let s = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(
            s.meniscus,
            step_newtonian_circular_meniscus(i.meniscus).unwrap()
        );
        assert_eq!(s.displaced_viscous_dissipation_j, 0.0);
        balance(s);
    }

    #[test]
    fn two_fluid_g3_split_steps_reverse_flow_and_phase_losses() {
        for mu_d in [1e-6, 0.001, 0.002, 0.00200000001, 0.01, 1.0] {
            for drive_sign in [-1.0, 1.0] {
                let mut i = input();
                i.displaced_dynamic_viscosity_pa_s = mu_d;
                i.meniscus.wet_length_m = 0.03;
                i.meniscus.contact_angle_cosine *= drive_sign;
                let full = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
                i.meniscus.duration_s *= 0.5;
                let first = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
                i.meniscus.wet_length_m = first.meniscus.wet_length_m;
                let second = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
                near(
                    full.meniscus.wet_length_m,
                    second.meniscus.wet_length_m,
                    2e-13,
                );
                near(
                    full.liquid_viscous_dissipation_j,
                    first.liquid_viscous_dissipation_j + second.liquid_viscous_dissipation_j,
                    2e-10,
                );
                near(
                    full.displaced_viscous_dissipation_j,
                    first.displaced_viscous_dissipation_j + second.displaced_viscous_dissipation_j,
                    2e-10,
                );
                balance(full);
            }
        }
    }

    #[test]
    fn two_fluid_end_events_stop_before_invalid_quadratic_branches() {
        let mut i = input();
        i.meniscus.duration_s = 1e6;
        i.displaced_dynamic_viscosity_pa_s = 1.0;
        let filled = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(filled.meniscus.end, CircularMeniscusEnd::Outlet);
        assert_eq!(filled.meniscus.wet_length_m, i.meniscus.tube_length_m);
        let m = i.meniscus;
        let time = 4.0
            * (m.dynamic_viscosity_pa_s + i.displaced_dynamic_viscosity_pa_s)
            * m.tube_length_m.powi(2)
            / (m.radius_m.powi(2) * filled.meniscus.capillary_pressure_pa);
        near(filled.meniscus.elapsed_time_s, time, 2e-14);
        i.meniscus.wet_length_m = m.tube_length_m;
        i.meniscus.contact_angle_cosine *= -1.0;
        let drained = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(drained.meniscus.end, CircularMeniscusEnd::Inlet);
        assert_eq!(drained.meniscus.wet_length_m, 0.0);
        near(drained.meniscus.elapsed_time_s, time, 2e-14);
        near(
            drained.liquid_viscous_dissipation_j,
            filled.liquid_viscous_dissipation_j,
            2e-14,
        );
        near(
            drained.displaced_viscous_dissipation_j,
            filled.displaced_viscous_dissipation_j,
            2e-14,
        );
        balance(filled);
        balance(drained);
    }

    #[test]
    fn two_fluid_phase_loss_matches_independent_spatial_quadrature() {
        let mut i = input();
        i.meniscus.wet_length_m = 0.01;
        i.meniscus.duration_s = 0.1;
        for mu_d in [0.001, 0.002, 0.01] {
            i.displaced_dynamic_viscosity_pa_s = mu_d;
            let step = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
            let m = i.meniscus;
            let dx = (step.meniscus.wet_length_m - m.wet_length_m) / 4096.0;
            let mut integral = 0.0;
            for k in 0..4096 {
                let x = m.wet_length_m + (f64::from(k) + 0.5) * dx;
                integral += m.dynamic_viscosity_pa_s * x
                    / (m.dynamic_viscosity_pa_s * x + mu_d * (m.tube_length_m - x));
            }
            near(
                step.liquid_viscous_dissipation_j,
                step.meniscus.viscous_dissipation_j * integral / 4096.0,
                2e-9,
            );
        }
    }

    #[test]
    fn two_fluid_stationary_and_invalid_inputs() {
        let mut i = input();
        i.meniscus.contact_angle_cosine = 0.0;
        let s = step_newtonian_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(s.meniscus.wet_length_m, 0.0);
        assert_eq!(s.meniscus.elapsed_time_s, i.meniscus.duration_s);
        assert_eq!(s.meniscus.viscous_dissipation_j, 0.0);
        for mu in [-1.0, f64::NAN, f64::INFINITY] {
            i.displaced_dynamic_viscosity_pa_s = mu;
            assert!(step_newtonian_two_fluid_circular_meniscus(i).is_err());
        }
    }
}
