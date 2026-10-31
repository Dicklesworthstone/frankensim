//! Quasisteady displacement of two fluids through an inclined circular tube.
//!
//! Both phases contribute Poiseuille resistance and hydrostatic head. The
//! separated travel-time integral is analytic; a bounded monotone inversion
//! retains finite end events and never crosses a stable equilibrium. Reservoir
//! pressures refer to the actual inlet and outlet elevations. Consequently the
//! gravitational account includes both tube storage and potential energy carried
//! through the outlet. Omitting the latter loses the displaced phase's head.
//!
//! Properties and contact angle are constant during one step. The model excludes
//! inertia, trapped-gas compression, deposited films and dynamic contact-line
//! dissipation, just as the horizontal two-fluid owner does.

use super::super::{
    CircularCapillaryError, CircularMeniscusEnd, CircularMeniscusStep,
    InclinedCircularMeniscusInput, step_newtonian_inclined_circular_meniscus,
    validate_circular_meniscus,
};
use super::{
    TwoFluidCircularMeniscusInput, TwoFluidCircularMeniscusStep,
    step_newtonian_two_fluid_circular_meniscus,
};
use core::f64::consts::PI;

/// Two-fluid displacement with independently specified phase densities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InclinedTwoFluidCircularMeniscusInput {
    /// Geometry, both viscosities, wetting and prescribed reservoir pressures.
    pub meniscus: TwoFluidCircularMeniscusInput,
    /// Density of the invading liquid [kg/m³], strictly positive.
    pub density_kg_m3: f64,
    /// Density of the displaced phase [kg/m³]. Zero explicitly neglects its
    /// hydrostatic head; it is not inferred from its viscosity or phase name.
    pub displaced_density_kg_m3: f64,
    /// Signed opposing gravity along the tube [m/s²], g sin(inclination).
    pub opposing_gravity_m_s2: f64,
}

/// Accepted displacement with separate open-system gravitational accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InclinedTwoFluidCircularMeniscusStep {
    /// Original motion/pressure/surface account and both viscous losses.
    pub meniscus: TwoFluidCircularMeniscusStep,
    /// Change of gravitational energy stored inside the tube [J], using zero
    /// potential at the inlet: g*(rho_l-rho_d)*mean(x)*transferred_volume.
    pub gravitational_storage_change_j: f64,
    /// Signed potential energy leaving through the outlet [J]:
    /// rho_d*g*L*transferred_volume. Negative flow imports this energy.
    pub outlet_potential_energy_j: f64,
    /// Total hydrostatic pressure across both columns at the accepted state [Pa].
    pub hydrostatic_pressure_pa: f64,
}

/// Advance the gravity-aware two-fluid Washburn equation.
///
/// `8*(mu_l*x + mu_d*(L-x))*dx/dt = R²*F(x)`, where
/// `F(x) = delta_p + 2*sigma*cos(theta)/R - g*(rho_l*x + rho_d*(L-x))`.
/// Positive `g*(rho_l-rho_d)` admits a stable equilibrium; negative density
/// contrast reverses its stability. Equal densities retain the full constant
/// head `rho*g*L`, even though tube gravitational storage does not change.
///
/// The mechanical account is
/// `W_pressure - delta_E_wall - delta_E_tube_gravity - E_outlet_potential = D`.
/// Both phase losses integrate their actual resistance share over the motion.
/// Zero gravity delegates bit-for-bit to the horizontal owner. The massless,
/// inviscid displaced-phase limit delegates to the original inclined owner.
/// No caller state is mutated and the bounded scalar kernel allocates nothing.
///
/// # Errors
/// Refuses malformed inputs or coefficients/transfers not representable in f64.
/// Reaching an equilibrium consumes the requested time, not an invented finite
/// arrival time. Inlet/outlet events consume only their actual travel time.
pub fn step_newtonian_inclined_two_fluid_circular_meniscus(
    input: InclinedTwoFluidCircularMeniscusInput,
) -> Result<InclinedTwoFluidCircularMeniscusStep, CircularCapillaryError> {
    validate_input(input)?;
    let m = input.meniscus.meniscus;
    let mu_d = input.meniscus.displaced_dynamic_viscosity_pa_s;
    let g = input.opposing_gravity_m_s2;
    if g == 0.0 {
        return Ok(InclinedTwoFluidCircularMeniscusStep {
            meniscus: step_newtonian_two_fluid_circular_meniscus(input.meniscus)?,
            gravitational_storage_change_j: 0.0,
            outlet_potential_energy_j: 0.0,
            hydrostatic_pressure_pa: 0.0,
        });
    }
    if mu_d == 0.0 && input.displaced_density_kg_m3 == 0.0 {
        let step = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
            meniscus: m,
            density_kg_m3: input.density_kg_m3,
            opposing_gravity_m_s2: g,
        })?;
        return Ok(InclinedTwoFluidCircularMeniscusStep {
            meniscus: TwoFluidCircularMeniscusStep {
                meniscus: step.meniscus,
                liquid_viscous_dissipation_j: step.meniscus.viscous_dissipation_j,
                displaced_viscous_dissipation_j: 0.0,
            },
            gravitational_storage_change_j: step.gravitational_energy_change_j,
            outlet_potential_energy_j: 0.0,
            hydrostatic_pressure_pa: step.hydrostatic_pressure_pa,
        });
    }

    let mu_scale = m.dynamic_viscosity_pa_s.max(mu_d);
    let area = PI * m.radius_m * m.radius_m;
    let capillary = if m.contact_angle_cosine == 0.0 {
        0.0
    } else {
        (m.surface_tension_n_m / m.radius_m) * (2.0 * m.contact_angle_cosine)
    };
    let hydrostatic = |x: f64| {
        g * input
            .density_kg_m3
            .mul_add(x, input.displaced_density_kg_m3 * (m.tube_length_m - x))
    };
    let outlet_head = (input.displaced_density_kg_m3 * g) * m.tube_length_m;
    let inlet_drive = (m.applied_pressure_pa - outlet_head) + capillary;
    let mut law = Law {
        initial: m.wet_length_m,
        length: m.tube_length_m,
        liquid: m.dynamic_viscosity_pa_s / mu_scale,
        displaced: mu_d / mu_scale,
        contact_line: 0.0,
        mobility: (m.radius_m / mu_scale) * (m.radius_m / 8.0),
        gravity_slope: (input.density_kg_m3 - input.displaced_density_kg_m3) * g,
        force0: (m.applied_pressure_pa - hydrostatic(m.wet_length_m)) + capillary,
    };
    // The only zero-resistance endpoint is an empty tube with mu_d=0. If the
    // inlet drive also vanishes, x cancels analytically before taking x -> 0.
    let constant_velocity = law.displaced == 0.0 && inlet_drive == 0.0;
    if law.gravity_slope == 0.0 {
        // Equal densities leave a constant full-tube head. Preserve the same
        // pressure arithmetic as the horizontal delegate, including balance.
        law.force0 = inlet_drive;
    } else if constant_velocity {
        law.force0 = -law.gravity_slope * law.initial;
    }
    if [
        area,
        capillary,
        outlet_head,
        inlet_drive,
        law.mobility,
        law.gravity_slope,
        law.force0,
        hydrostatic(0.0),
        hydrostatic(law.length),
    ]
    .into_iter()
    .any(|v| !v.is_finite())
        || area <= 0.0
        || law.mobility <= 0.0
        || law.liquid <= 0.0
        || (mu_d > 0.0 && law.displaced <= 0.0)
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    let (length, elapsed, end) = if m.duration_s == 0.0 {
        (law.initial, 0.0, CircularMeniscusEnd::TimeLimit)
    } else if constant_velocity && law.gravity_slope != 0.0 {
        let speed = -law.mobility * law.gravity_slope / law.liquid;
        if !speed.is_finite() || speed == 0.0 {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        let (boundary, event) = if speed > 0.0 {
            (law.length, CircularMeniscusEnd::Outlet)
        } else {
            (0.0, CircularMeniscusEnd::Inlet)
        };
        let travel = (boundary - law.initial) / speed;
        if travel <= m.duration_s {
            (boundary, travel, event)
        } else {
            (
                speed.mul_add(m.duration_s, law.initial),
                m.duration_s,
                CircularMeniscusEnd::TimeLimit,
            )
        }
    } else if law.force0 == 0.0 {
        (law.initial, m.duration_s, CircularMeniscusEnd::TimeLimit)
    } else if law.gravity_slope == 0.0 {
        let horizontal =
            step_newtonian_two_fluid_circular_meniscus(TwoFluidCircularMeniscusInput {
                meniscus: super::super::CircularMeniscusInput {
                    applied_pressure_pa: m.applied_pressure_pa - outlet_head,
                    ..m
                },
                ..input.meniscus
            })?
            .meniscus;
        (
            horizontal.wet_length_m,
            horizontal.elapsed_time_s,
            horizontal.end,
        )
    } else {
        law.advance(m.duration_s)?
    };

    let delta = length - law.initial;
    let volume = area * delta;
    let final_drive = law.force(delta);
    let pressure_work = m.applied_pressure_pa * volume;
    let wall_energy = -capillary * volume;
    let mean_length = law.initial + 0.5 * delta;
    let gravity_storage = (law.gravity_slope * mean_length) * volume;
    let gravity_export = outlet_head * volume;
    let dissipation = (0.5 * law.force0 + 0.5 * final_drive) * volume;
    let (liquid_loss, displaced_loss, _) = phase_losses(law, length, dissipation)?;
    let final_head = hydrostatic(length);
    if [
        length,
        elapsed,
        volume,
        pressure_work,
        wall_energy,
        gravity_storage,
        gravity_export,
        dissipation,
        liquid_loss,
        displaced_loss,
        final_head,
    ]
    .into_iter()
    .any(|v| !v.is_finite())
        || !(0.0..=law.length).contains(&length)
        || elapsed < 0.0
        || dissipation < 0.0
        || liquid_loss < 0.0
        || displaced_loss < 0.0
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    Ok(InclinedTwoFluidCircularMeniscusStep {
        meniscus: TwoFluidCircularMeniscusStep {
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
        },
        gravitational_storage_change_j: gravity_storage,
        outlet_potential_energy_j: gravity_export,
        hydrostatic_pressure_pa: final_head,
    })
}

pub(super) fn validate_input(
    input: InclinedTwoFluidCircularMeniscusInput,
) -> Result<(), CircularCapillaryError> {
    validate_circular_meniscus(input.meniscus.meniscus)?;
    let mu_d = input.meniscus.displaced_dynamic_viscosity_pa_s;
    for (field, value) in [
        ("density_kg_m3", input.density_kg_m3),
        ("displaced_density_kg_m3", input.displaced_density_kg_m3),
        ("opposing_gravity_m_s2", input.opposing_gravity_m_s2),
        ("displaced_dynamic_viscosity_pa_s", mu_d),
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
    for (field, value) in [
        ("displaced_density_kg_m3", input.displaced_density_kg_m3),
        ("displaced_dynamic_viscosity_pa_s", mu_d),
    ] {
        if value < 0.0 {
            return Err(CircularCapillaryError::OutsideDomain { field });
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(super) struct Law {
    pub(super) initial: f64,
    pub(super) length: f64,
    pub(super) liquid: f64,
    pub(super) displaced: f64,
    /// Constant extra resistance, normalized by the same viscosity scale.
    pub(super) contact_line: f64,
    pub(super) mobility: f64,
    pub(super) gravity_slope: f64,
    pub(super) force0: f64,
}

impl Law {
    pub(super) fn resistance(self, x: f64) -> f64 {
        self.liquid.mul_add(x, self.displaced * (self.length - x)) + self.contact_line
    }

    pub(super) fn force(self, delta: f64) -> f64 {
        (-self.gravity_slope).mul_add(delta, self.force0)
    }

    fn travel_time(self, delta: f64) -> Result<f64, CircularCapillaryError> {
        if delta == 0.0 {
            return Ok(0.0);
        }
        let ratio = delta / self.force0;
        let z = self.gravity_slope * ratio;
        if z >= 1.0 {
            return Ok(f64::INFINITY);
        }
        if !z.is_finite() || !ratio.is_finite() {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        let (left, right) = if z.abs() <= 0.125 {
            let (mut left, mut right, mut power) = (0.5, 0.5, 1.0);
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
        let time = ratio
            * (self.resistance(self.initial) * left
                + self.resistance(self.initial + delta) * right)
            / self.mobility;
        if time.is_nan() || time < 0.0 {
            Err(CircularCapillaryError::NonFiniteResult)
        } else {
            Ok(time)
        }
    }

    pub(super) fn advance(
        self,
        duration: f64,
    ) -> Result<(f64, f64, CircularMeniscusEnd), CircularCapillaryError> {
        let direction = self.force0.signum();
        let (boundary, event) = if direction > 0.0 {
            (self.length, CircularMeniscusEnd::Outlet)
        } else {
            (0.0, CircularMeniscusEnd::Inlet)
        };
        if boundary == self.initial {
            return Ok((boundary, 0.0, event));
        }
        let to_boundary = (boundary - self.initial).abs();
        let equilibrium_offset = self.force0 / self.gravity_slope;
        let barrier = self.gravity_slope > 0.0 && equilibrium_offset.abs() <= to_boundary;
        let limit = if barrier {
            equilibrium_offset.abs()
        } else {
            to_boundary
        };
        if !barrier {
            let travel = self.travel_time(direction * limit)?;
            if travel <= duration {
                return Ok((boundary, travel, event));
            }
        }
        let target = direction.mul_add(limit, self.initial);
        let minimum_resistance = self.resistance(self.initial).min(self.resistance(target));
        let maximum_force = self.force0.abs().max(self.force(direction * limit).abs());
        let bound = (self.mobility * maximum_force / minimum_resistance) * duration;
        let mut high = if bound.is_finite() && bound > 0.0 {
            bound.min(limit)
        } else {
            limit
        };
        if high < limit && self.travel_time(direction * high)? < duration {
            high = limit;
        }
        let mut low = 0.0;
        let initial_resistance = self.resistance(self.initial);
        let estimate = if initial_resistance > 0.0 {
            (self.mobility * self.force0.abs() / initial_resistance) * duration
        } else {
            ((2.0 * self.mobility * self.force0.abs()) * duration / self.liquid).sqrt()
        };
        let mut distance = if estimate.is_finite() && estimate > 0.0 && estimate < high {
            estimate
        } else {
            0.5 * high
        };
        for iteration in 0..80 {
            let travel = self.travel_time(direction * distance)?;
            if (travel - duration).abs() <= 8.0 * f64::EPSILON * duration {
                break;
            }
            if travel < duration {
                low = distance;
            } else {
                high = distance;
            }
            let midpoint = low + 0.5 * (high - low);
            if midpoint == low
                || midpoint == high
                || direction.mul_add(low, self.initial) == direction.mul_add(high, self.initial)
            {
                distance = midpoint;
                break;
            }
            let position = direction.mul_add(distance, self.initial);
            let inverse_derivative =
                self.mobility * self.force(direction * distance).abs() / self.resistance(position);
            let candidate = distance - (travel - duration) * inverse_derivative;
            distance =
                if iteration < 16 && candidate.is_finite() && candidate > low && candidate < high {
                    candidate
                } else {
                    midpoint
                };
        }
        let mut length = direction
            .mul_add(distance, self.initial)
            .clamp(self.initial.min(target), self.initial.max(target));
        // A rounded representation of the equilibrium may lie one ulp beyond
        // the root of the anchored force. Retain the admissible side explicitly.
        if barrier && direction * self.force(length - self.initial) < 0.0 {
            length = if direction > 0.0 {
                length.next_down()
            } else {
                length.next_up()
            };
        }
        if direction * self.force(length - self.initial) < 0.0 {
            return Err(CircularCapillaryError::NonFiniteResult);
        }
        Ok((length, duration, CircularMeniscusEnd::TimeLimit))
    }
}

// Integrate F(u)*R_phase(u)/R_total(u). Three positive quadratic Bernstein
// weights retain small phase losses without subtracting large antiderivatives.
pub(super) fn phase_losses(
    law: Law,
    end: f64,
    total: f64,
) -> Result<(f64, f64, f64), CircularCapillaryError> {
    if total == 0.0 {
        return Ok((0.0, 0.0, 0.0));
    }
    if law.displaced == 0.0 && law.contact_line == 0.0 {
        return Ok((total, 0.0, 0.0));
    }
    let mut resistance = [law.resistance(law.initial), law.resistance(end)];
    let scale = resistance[0].max(resistance[1]);
    let mut force = [law.force0.abs(), law.force(end - law.initial).abs()];
    let force_scale = force[0].max(force[1]);
    force = force.map(|f| f / force_scale);
    let mut liquid = [law.liquid * law.initial / scale, law.liquid * end / scale];
    let mut displaced = [
        law.displaced * (law.length - law.initial) / scale,
        law.displaced * (law.length - end) / scale,
    ];
    if resistance[0] < resistance[1] {
        resistance.swap(0, 1);
        force.swap(0, 1);
        liquid.swap(0, 1);
        displaced.swap(0, 1);
    }
    let ratio = resistance[1] / resistance[0];
    if !ratio.is_finite() || ratio <= 0.0 || ratio > 1.0 {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    let d = 1.0 - ratio;
    let (left, cross, right) = if d < 0.125 {
        let (mut left, mut cross, mut right, mut power) = (1.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0, 1.0);
        for k in 1..=20 {
            power *= d;
            left += 2.0 * power / f64::from((k + 1) * (k + 2) * (k + 3));
            cross += power / f64::from((k + 2) * (k + 3));
            right += power / f64::from(k + 3);
        }
        (left, cross, right)
    } else {
        let logarithm = -ratio.ln();
        let cube = d * d * d;
        (
            (0.5 * d * d - ratio * d + ratio * ratio * logarithm) / cube,
            (0.5 * d * d - ratio * (logarithm - d)) / cube,
            (logarithm - d - 0.5 * d * d) / cube,
        )
    };
    let integral = |phase: [f64; 2]| {
        force[0] * phase[0] * left
            + (force[0] * phase[1] + force[1] * phase[0]) * cross
            + force[1] * phase[1] * right
    };
    let a = integral(liquid);
    let b = integral(displaced);
    let c = if law.contact_line == 0.0 {
        0.0
    } else {
        integral([law.contact_line / scale; 2])
    };
    let mean_force = 0.5 * force[0] + 0.5 * force[1];
    if [a, b, c].into_iter().any(|v| !v.is_finite() || v < 0.0) || mean_force <= 0.0 {
        return Err(CircularCapillaryError::NonFiniteResult);
    }
    // Compute small channels directly. Recover the largest by closure so tiny
    // bulk or line losses are not erased by subtracting two large integrals.
    if a >= b && a >= c {
        let displaced_loss = total * (b / mean_force);
        let line_loss = total * (c / mean_force);
        Ok((
            total - displaced_loss - line_loss,
            displaced_loss,
            line_loss,
        ))
    } else if b >= c {
        let liquid_loss = total * (a / mean_force);
        let line_loss = total * (c / mean_force);
        Ok((liquid_loss, total - liquid_loss - line_loss, line_loss))
    } else {
        let liquid_loss = total * (a / mean_force);
        let displaced_loss = total * (b / mean_force);
        Ok((
            liquid_loss,
            displaced_loss,
            total - liquid_loss - displaced_loss,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::CircularMeniscusInput;
    use super::*;

    fn input() -> InclinedTwoFluidCircularMeniscusInput {
        InclinedTwoFluidCircularMeniscusInput {
            meniscus: TwoFluidCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    dynamic_viscosity_pa_s: 0.1,
                    tube_length_m: 0.02,
                    radius_m: 0.001,
                    wet_length_m: 0.004,
                    surface_tension_n_m: 0.072,
                    contact_angle_cosine: 0.5,
                    applied_pressure_pa: 100.0,
                    duration_s: 1.0,
                },
                displaced_dynamic_viscosity_pa_s: 0.2,
            },
            density_kg_m3: 1000.0,
            displaced_density_kg_m3: 800.0,
            opposing_gravity_m_s2: 9.81,
        }
    }
    fn near(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-20),
            "{a} vs {b}"
        );
    }
    fn balance(s: InclinedTwoFluidCircularMeniscusStep) {
        let m = s.meniscus.meniscus;
        near(
            m.pressure_work_j
                - m.wall_surface_energy_change_j
                - s.gravitational_storage_change_j
                - s.outlet_potential_energy_j,
            m.viscous_dissipation_j,
            2e-12,
        );
        near(
            s.meniscus.liquid_viscous_dissipation_j + s.meniscus.displaced_viscous_dissipation_j,
            m.viscous_dissipation_j,
            2e-15,
        );
        assert!(s.meniscus.liquid_viscous_dissipation_j >= 0.0);
        assert!(s.meniscus.displaced_viscous_dissipation_j >= 0.0);
    }

    #[test]
    fn gravity_two_fluid_g3_exact_horizontal_and_negligible_outer_limits() {
        let mut i = input();
        for mu in [0.0, 0.001, 0.2] {
            i.meniscus.displaced_dynamic_viscosity_pa_s = mu;
            i.opposing_gravity_m_s2 = 0.0;
            let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            assert_eq!(
                s.meniscus,
                step_newtonian_two_fluid_circular_meniscus(i.meniscus).unwrap()
            );
            assert_eq!(s.gravitational_storage_change_j, 0.0);
            assert_eq!(s.outlet_potential_energy_j, 0.0);
        }
        i = input();
        i.meniscus.displaced_dynamic_viscosity_pa_s = 0.0;
        i.displaced_density_kg_m3 = 0.0;
        let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
        let one = step_newtonian_inclined_circular_meniscus(InclinedCircularMeniscusInput {
            meniscus: i.meniscus.meniscus,
            density_kg_m3: i.density_kg_m3,
            opposing_gravity_m_s2: i.opposing_gravity_m_s2,
        })
        .unwrap();
        assert_eq!(s.meniscus.meniscus, one.meniscus);
        assert_eq!(
            s.gravitational_storage_change_j,
            one.gravitational_energy_change_j
        );
        assert_eq!(s.meniscus.displaced_viscous_dissipation_j, 0.0);
    }

    #[test]
    fn gravity_two_fluid_g1_equal_density_retains_downstream_hydrostatic_head() {
        let mut i = input();
        i.displaced_density_kg_m3 = i.density_kg_m3;
        i.meniscus.meniscus.duration_s = 0.1;
        let head = i.density_kg_m3 * i.opposing_gravity_m_s2 * i.meniscus.meniscus.tube_length_m;
        let mut horizontal = i.meniscus;
        horizontal.meniscus.applied_pressure_pa -= head;
        let expected = step_newtonian_two_fluid_circular_meniscus(horizontal).unwrap();
        let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(
            s.meniscus.meniscus.wet_length_m,
            expected.meniscus.wet_length_m
        );
        assert!(s.meniscus.meniscus.transferred_volume_m3 < 0.0);
        assert_eq!(s.gravitational_storage_change_j, 0.0);
        near(
            s.outlet_potential_energy_j,
            head * s.meniscus.meniscus.transferred_volume_m3,
            2e-15,
        );
        near(s.hydrostatic_pressure_pa, head, 2e-15);
        balance(s);
    }

    #[test]
    fn gravity_two_fluid_g1_equal_viscosity_matches_exponential_and_unstable_branch() {
        for rho_d in [800.0, 1200.0] {
            let mut i = input();
            i.displaced_density_kg_m3 = rho_d;
            i.meniscus.displaced_dynamic_viscosity_pa_s =
                i.meniscus.meniscus.dynamic_viscosity_pa_s;
            i.meniscus.meniscus.applied_pressure_pa = if rho_d > 1000.0 { 145.0 } else { 100.0 };
            let m = i.meniscus.meniscus;
            let b = (i.density_kg_m3 - rho_d) * i.opposing_gravity_m_s2;
            let c = m.applied_pressure_pa
                + 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m
                - rho_d * i.opposing_gravity_m_s2 * m.tube_length_m;
            let lambda =
                m.radius_m.powi(2) * b / (8.0 * m.dynamic_viscosity_pa_s * m.tube_length_m);
            let delta = (c - b * m.wet_length_m) / b * -(-lambda * m.duration_s).exp_m1();
            let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            near(
                s.meniscus.meniscus.wet_length_m,
                m.wet_length_m + delta,
                2e-13,
            );
            balance(s);
            if b < 0.0 {
                assert!(s.meniscus.meniscus.wet_length_m < m.wet_length_m);
            }
        }
    }

    #[test]
    fn gravity_two_fluid_g1_stable_equilibrium_and_split_time_do_not_overshoot() {
        let base = input();
        let m = base.meniscus.meniscus;
        let equilibrium = (m.applied_pressure_pa
            + 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m
            - base.displaced_density_kg_m3 * base.opposing_gravity_m_s2 * m.tube_length_m)
            / ((base.density_kg_m3 - base.displaced_density_kg_m3) * base.opposing_gravity_m_s2);
        for initial in [0.002, 0.015] {
            let mut i = base;
            i.meniscus.meniscus.wet_length_m = initial;
            i.meniscus.meniscus.duration_s = 1e5;
            let long = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            near(long.meniscus.meniscus.wet_length_m, equilibrium, 2e-14);
            assert_eq!(long.meniscus.meniscus.end, CircularMeniscusEnd::TimeLimit);
            balance(long);
            let settled = step_newtonian_inclined_two_fluid_circular_meniscus(
                InclinedTwoFluidCircularMeniscusInput {
                    meniscus: TwoFluidCircularMeniscusInput {
                        meniscus: CircularMeniscusInput {
                            wet_length_m: long.meniscus.meniscus.wet_length_m,
                            ..i.meniscus.meniscus
                        },
                        ..i.meniscus
                    },
                    ..i
                },
            )
            .unwrap();
            near(settled.meniscus.meniscus.wet_length_m, equilibrium, 2e-14);
            assert_eq!(
                settled.meniscus.meniscus.end,
                CircularMeniscusEnd::TimeLimit
            );
            i.meniscus.meniscus.duration_s = 2.0;
            let whole = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            i.meniscus.meniscus.duration_s = 1.0;
            let a = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            i.meniscus.meniscus.wet_length_m = a.meniscus.meniscus.wet_length_m;
            let b = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            near(
                whole.meniscus.meniscus.wet_length_m,
                b.meniscus.meniscus.wet_length_m,
                3e-13,
            );
            near(
                whole.meniscus.liquid_viscous_dissipation_j,
                a.meniscus.liquid_viscous_dissipation_j + b.meniscus.liquid_viscous_dissipation_j,
                3e-12,
            );
            near(
                whole.meniscus.displaced_viscous_dissipation_j,
                a.meniscus.displaced_viscous_dissipation_j
                    + b.meniscus.displaced_viscous_dissipation_j,
                3e-12,
            );
        }
    }

    #[test]
    fn gravity_two_fluid_g1_phase_losses_match_independent_spatial_quadrature() {
        for mu_d in [0.001, 0.1, 10.0] {
            for pressure in [-100.0, 100.0] {
                let mut i = input();
                i.meniscus.displaced_dynamic_viscosity_pa_s = mu_d;
                i.meniscus.meniscus.applied_pressure_pa = pressure;
                i.meniscus.meniscus.duration_s = 0.01;
                let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
                let m = i.meniscus.meniscus;
                let delta = s.meniscus.meniscus.wet_length_m - m.wet_length_m;
                let mut loss = [0.0; 2];
                for k in 0..2048 {
                    let x = m.wet_length_m + delta * (f64::from(k) + 0.5) / 2048.0;
                    let r = [m.dynamic_viscosity_pa_s * x, mu_d * (m.tube_length_m - x)];
                    let force = m.applied_pressure_pa
                        + 2.0 * m.surface_tension_n_m * m.contact_angle_cosine / m.radius_m
                        - i.opposing_gravity_m_s2
                            * (i.density_kg_m3 * x
                                + i.displaced_density_kg_m3 * (m.tube_length_m - x));
                    for p in 0..2 {
                        loss[p] +=
                            PI * m.radius_m.powi(2) * delta / 2048.0 * force * r[p] / (r[0] + r[1]);
                    }
                }
                near(s.meniscus.liquid_viscous_dissipation_j, loss[0], 3e-7);
                near(s.meniscus.displaced_viscous_dissipation_j, loss[1], 3e-7);
                balance(s);
            }
        }
    }

    #[test]
    fn gravity_two_fluid_end_events_and_inviscid_displaced_hydrostatic_limit() {
        for gravity in [-9.81, 9.81] {
            let mut i = input();
            i.opposing_gravity_m_s2 = gravity;
            i.meniscus.meniscus.applied_pressure_pa = if gravity > 0.0 { -100.0 } else { 100.0 };
            i.meniscus.meniscus.duration_s = 1000.0;
            let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
            assert!(s.meniscus.meniscus.elapsed_time_s < 1000.0);
            assert_eq!(
                s.meniscus.meniscus.wet_length_m,
                if gravity > 0.0 { 0.0 } else { 0.02 }
            );
            balance(s);
        }
        let mut i = input();
        i.meniscus.displaced_dynamic_viscosity_pa_s = 0.0;
        i.meniscus.meniscus.contact_angle_cosine = 0.0;
        i.meniscus.meniscus.applied_pressure_pa =
            i.displaced_density_kg_m3 * i.opposing_gravity_m_s2 * i.meniscus.meniscus.tube_length_m;
        let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
        let m = i.meniscus.meniscus;
        let speed = -m.radius_m.powi(2)
            * (i.density_kg_m3 - i.displaced_density_kg_m3)
            * i.opposing_gravity_m_s2
            / (8.0 * m.dynamic_viscosity_pa_s);
        near(
            s.meniscus.meniscus.wet_length_m,
            m.wet_length_m + speed * m.duration_s,
            2e-14,
        );
        assert_eq!(s.meniscus.displaced_viscous_dissipation_j, 0.0);
        balance(s);
    }

    #[test]
    fn gravity_two_fluid_stationary_duration_and_invalid_inputs() {
        let mut i = input();
        i.meniscus.meniscus.duration_s = 0.0;
        let s = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(
            s.meniscus.meniscus.wet_length_m,
            i.meniscus.meniscus.wet_length_m
        );
        assert_eq!(s.meniscus.meniscus.transferred_volume_m3, 0.0);
        balance(s);
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let mut i = input();
            i.displaced_density_kg_m3 = bad;
            assert!(step_newtonian_inclined_two_fluid_circular_meniscus(i).is_err());
        }
        // A hydrostatically balanced equal-density tube stays stationary even
        // when the displaced phase is inviscid; the canceled velocity is zero.
        i = input();
        i.displaced_density_kg_m3 = i.density_kg_m3;
        i.meniscus.displaced_dynamic_viscosity_pa_s = 0.0;
        i.meniscus.meniscus.contact_angle_cosine = 0.0;
        i.meniscus.meniscus.applied_pressure_pa =
            i.density_kg_m3 * i.opposing_gravity_m_s2 * i.meniscus.meniscus.tube_length_m;
        let stationary = step_newtonian_inclined_two_fluid_circular_meniscus(i).unwrap();
        assert_eq!(
            stationary.meniscus.meniscus.wet_length_m,
            i.meniscus.meniscus.wet_length_m
        );
        assert_eq!(
            stationary.meniscus.meniscus.elapsed_time_s,
            i.meniscus.meniscus.duration_s
        );
        balance(stationary);
        i = input();
        i.density_kg_m3 = 0.0;
        assert!(step_newtonian_inclined_two_fluid_circular_meniscus(i).is_err());
        i = input();
        i.opposing_gravity_m_s2 = f64::NAN;
        assert!(step_newtonian_inclined_two_fluid_circular_meniscus(i).is_err());
    }
}
