//! Gravity-resolved, quasi-steady capillary rise and drainage.
//!
//! The original pipe owner resolves Q(x) under P0 - rho g x at every sample.
//! A logarithmic distance to a stable hydrostatic equilibrium regularizes the
//! infinite physical travel time. A finite step never jumps across or snaps
//! to that equilibrium. Density and projected gravity are explicit SI inputs;
//! this module adds no density, wetting or viscosity calibration authority.

use super::{
    CarreauCapillaryTube, CarreauMeniscusControls, CarreauMeniscusEnd, CarreauMeniscusError,
    CarreauMeniscusInput, CarreauMeniscusStep, PathSample, Travel,
};

/// Accepted column transport with constant density and signed axial gravity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauInclinedMeniscusStep {
    pub meniscus: CarreauMeniscusStep,
    pub density_kg_m3: f64,
    /// Positive gravity opposes increasing wet length (an upward-pointing tube).
    pub opposing_gravity_m_s2: f64,
    /// Signed mass entering the column, rho times transferred volume.
    pub transferred_mass_kg: f64,
    /// Signed hydrostatic pressure across the final liquid column, rho g x.
    pub hydrostatic_pressure_pa: f64,
}

#[derive(Clone, Copy)]
enum Coordinate {
    Distance {
        boundary: f64,
        maximum: f64,
    },
    /// y = -ln(|x - x_eq| / |x0 - x_eq|).
    Equilibrium {
        distance: f64,
        maximum: f64,
        equilibrium: f64,
    },
}

impl Coordinate {
    fn maximum(self) -> f64 {
        match self {
            Self::Distance { maximum, .. } | Self::Equilibrium { maximum, .. } => maximum,
        }
    }

    fn point(
        self,
        coordinate: f64,
        x0: f64,
        direction: f64,
        pressure0: f64,
        drive0: f64,
        gravity_gradient: f64,
    ) -> Result<(f64, f64, f64), CarreauMeniscusError> {
        let point = match self {
            Self::Distance { boundary, maximum } => {
                let x = if coordinate == maximum {
                    boundary
                } else {
                    x0 + direction * coordinate
                };
                (x, (-gravity_gradient).mul_add(x, pressure0), 1.0)
            }
            Self::Equilibrium {
                distance,
                equilibrium,
                ..
            } => {
                let decay = fs_math::det::exp(-coordinate);
                let x = x0 + direction * distance * -fs_math::det::expm1(-coordinate);
                if x == equilibrium {
                    return Err(CarreauMeniscusError::Unrepresentable);
                }
                // Do not subtract two nearly equal hydrostatic pressures.
                (x, drive0 * decay, distance * decay)
            }
        };
        if !point.0.is_finite()
            || !point.1.is_finite()
            || !point.2.is_finite()
            || point.2 <= 0.0
            || point.1 == 0.0
            || point.1.signum() != direction
            || (coordinate > 0.0 && point.0 == x0)
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        Ok(point)
    }
}

// The equilibrium coordinate regularizes A dx/|Q|, but makes the heat
// density decay as exp(-2z). Resolve time in z and heat in physical x:
// both integrands are smooth even after many relaxation times. The same
// positive Simpson rule, two-agreement gate and caller work bounds apply.
fn integrate_positive(
    span: f64,
    endpoints: [f64; 2],
    controls: CarreauMeniscusControls,
    mut sample: impl FnMut(f64) -> Result<f64, CarreauMeniscusError>,
) -> Result<(f64, f64), CarreauMeniscusError> {
    let ([value], change) = super::integrate_positive(
        span,
        [[endpoints[0]], [endpoints[1]]],
        controls,
        |fraction| Ok([sample(fraction)?]),
    )?;
    Ok((value, change))
}

impl CarreauCapillaryTube {
    /// Resolve constant-density gravity, capillary pressure and imposed
    /// pressure together with the original thermal Carreau radial profile.
    ///
    /// The opposing gravitational acceleration is g sin(inclination) for a
    /// straight tube; either sign is allowed. All fluid remains at the pipe's
    /// declared temperature. There is no inertia, entrance correction,
    /// displaced-fluid column, variable density or dynamic contact angle.
    ///
    /// Both endpoints of the declared wet interval must admit the original
    /// stress/rate support. A stable equilibrium is approached in logarithmic
    /// distance, never crossed in a finite step. Travel below representable
    /// position resolution is refused, rather than clipped into equilibrium.
    pub fn step_inclined(
        &self,
        input: CarreauMeniscusInput,
        density_kg_m3: f64,
        opposing_gravity_m_s2: f64,
        controls: CarreauMeniscusControls,
    ) -> Result<CarreauInclinedMeniscusStep, CarreauMeniscusError> {
        controls.admit()?;
        if !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
            || !opposing_gravity_m_s2.is_finite()
            || !input.wet_length_m.is_finite()
            || !(self.minimum_wet_length_m..=self.pipe.length_m()).contains(&input.wet_length_m)
            || !input.applied_pressure_pa.is_finite()
            || !input.duration_s.is_finite()
            || input.duration_s < 0.0
        {
            return Err(CarreauMeniscusError::Input(
                "invalid inclined meniscus inputs",
            ));
        }
        let gradient = density_kg_m3 * opposing_gravity_m_s2;
        let pressure0 = input.applied_pressure_pa + self.capillary_pressure_pa;
        if !gradient.is_finite() || !pressure0.is_finite() {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let wrap = |mut meniscus: CarreauMeniscusStep| {
            let transferred_mass_kg = density_kg_m3 * meniscus.transferred_volume_m3;
            let hydrostatic_pressure_pa = gradient * meniscus.wet_length_m;
            let mean_length = 0.5 * input.wet_length_m + 0.5 * meniscus.wet_length_m;
            let gravitational_energy_change_j = if gradient == 0.0 {
                0.0
            } else {
                meniscus.transferred_volume_m3 * (gradient * mean_length)
            };
            meniscus.gravitational_energy_change_j = gravitational_energy_change_j;
            meniscus.constitutive_work_defect_j = meniscus.reservoir_work_j
                - meniscus.wall_energy_change_j
                - gravitational_energy_change_j
                - meniscus.viscous_dissipation_j;
            if ![
                transferred_mass_kg,
                hydrostatic_pressure_pa,
                gravitational_energy_change_j,
                meniscus.constitutive_work_defect_j,
            ]
            .iter()
            .all(|v| v.is_finite())
            {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            Ok(CarreauInclinedMeniscusStep {
                meniscus,
                density_kg_m3,
                opposing_gravity_m_s2,
                transferred_mass_kg,
                hydrostatic_pressure_pa,
            })
        };
        // Preserve the existing horizontal computation exactly.
        if opposing_gravity_m_s2 == 0.0 {
            return wrap(self.step_horizontal(input, controls)?);
        }
        if gradient == 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let drive0 = (-gradient).mul_add(input.wet_length_m, pressure0);
        if !drive0.is_finite() {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut evaluations = 0;
        let start_profile = self.profile(input.wet_length_m, drive0, controls, &mut evaluations)?;
        // |P0/x - rho*g| has its maximum at an interval endpoint. Admit the
        // entire declared interval even for an initially stationary column.
        let mut boundary_profiles = [start_profile; 2];
        for (i, x) in [self.minimum_wet_length_m, self.pipe.length_m()]
            .into_iter()
            .enumerate()
        {
            if x != input.wet_length_m {
                let drive = (-gradient).mul_add(x, pressure0);
                boundary_profiles[i] = self.profile(x, drive, controls, &mut evaluations)?;
            }
        }
        if drive0 == 0.0 || input.duration_s == 0.0 {
            let end = if drive0 == 0.0 {
                CarreauMeniscusEnd::Stationary
            } else {
                CarreauMeniscusEnd::Duration
            };
            return wrap(self.finish(
                input,
                input.duration_s,
                end,
                Travel {
                    time_s: input.duration_s,
                    heat_j: 0.0,
                    relative_change: 0.0,
                    endpoint: start_profile,
                },
                evaluations,
                0,
            )?);
        }
        let direction = drive0.signum();
        let boundary_index = if direction > 0.0 { 1 } else { 0 };
        let boundary_profile = boundary_profiles[boundary_index];
        let boundary = boundary_profile.pipe().length_m();
        let boundary_kind = if direction > 0.0 {
            CarreauMeniscusEnd::Full
        } else {
            CarreauMeniscusEnd::MinimumWetLength
        };
        let maximum_distance = (boundary - input.wet_length_m).abs();
        if maximum_distance == 0.0 {
            return wrap(self.finish(
                input,
                0.0,
                boundary_kind,
                Travel {
                    time_s: 0.0,
                    heat_j: 0.0,
                    relative_change: 0.0,
                    endpoint: start_profile,
                },
                evaluations,
                0,
            )?);
        }
        let equilibrium_distance = drive0.abs() / gradient;
        let coordinate = if gradient > 0.0 && equilibrium_distance <= maximum_distance {
            let equilibrium = input.wet_length_m + direction * equilibrium_distance;
            // Keep a finite representable distance from equilibrium. This is
            // an admissibility limit, never an altered physical endpoint.
            let resolution = 4.0 * f64::EPSILON * equilibrium.abs().max(input.wet_length_m.abs());
            let maximum = fs_math::det::ln(equilibrium_distance) - fs_math::det::ln(resolution);
            if !maximum.is_finite() || maximum <= 0.0 {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            Coordinate::Equilibrium {
                distance: equilibrium_distance,
                maximum,
                equilibrium,
            }
        } else {
            Coordinate::Distance {
                boundary,
                maximum: maximum_distance,
            }
        };
        let make_sample = |value: f64, evaluations: &mut u32| {
            let (x, drive, length_jacobian) = coordinate.point(
                value,
                input.wet_length_m,
                direction,
                pressure0,
                drive0,
                gradient,
            )?;
            let profile = if value == 0.0 {
                start_profile
            } else {
                self.profile(x, drive, controls, evaluations)?
            };
            Ok(PathSample {
                profile,
                length_jacobian,
            })
        };
        let start = make_sample(0.0, &mut evaluations)?;
        let at_coordinate = |value: f64, evaluations: &mut u32| {
            let endpoint = make_sample(value, evaluations)?;
            if let Coordinate::Equilibrium { distance, .. } = coordinate {
                let (time_s, time_change) = integrate_positive(
                    value,
                    [self.integrands(&start)?[0], self.integrands(&endpoint)?[0]],
                    controls,
                    |fraction| {
                        Ok(self.integrands(&make_sample(value * fraction, evaluations)?)?[0])
                    },
                )?;
                let wet_distance = distance * -fs_math::det::expm1(-value);
                let decay = fs_math::det::exp(-value);
                let (heat_j, heat_change) = integrate_positive(
                    wet_distance,
                    [
                        self.area_m2 * start.profile.flow().pressure_drop_pa.abs(),
                        self.area_m2 * endpoint.profile.flow().pressure_drop_pa.abs(),
                    ],
                    controls,
                    |fraction| {
                        let x = input.wet_length_m + direction * wet_distance * fraction;
                        // Stable affine physical pressure, including near the
                        // equilibrium endpoint. The original pipe still owns
                        // its accepted pressure and constitutive work residual.
                        let drive = drive0 * ((1.0 - fraction) + fraction * decay);
                        let profile = self.profile(x, drive, controls, evaluations)?;
                        Ok(self.area_m2 * profile.flow().pressure_drop_pa.abs())
                    },
                )?;
                return Ok(Travel {
                    time_s,
                    heat_j,
                    relative_change: time_change.max(heat_change),
                    endpoint: endpoint.profile,
                });
            }
            self.integrate_travel(start, endpoint, value, controls, |fraction| {
                make_sample(value * fraction, evaluations)
            })
        };
        let maximum = coordinate.maximum();
        let mut upper = match coordinate {
            // Use the physical tangent estimate directly; imposing a unit
            // coordinate cap repeatedly reintegrates prefixes on long steps.
            Coordinate::Equilibrium { .. } => {
                (input.duration_s / self.integrands(&start)?[0]).min(maximum)
            }
            Coordinate::Distance { .. } => {
                let speed = start_profile
                    .flow()
                    .volumetric_flow_m3_s
                    .abs()
                    .max(boundary_profile.flow().volumetric_flow_m3_s.abs())
                    / self.area_m2;
                if !speed.is_finite() || speed <= 0.0 {
                    return Err(CarreauMeniscusError::Unrepresentable);
                }
                if input.duration_s >= maximum / speed {
                    maximum
                } else {
                    input.duration_s * speed
                }
            }
        };
        if !upper.is_finite() || upper <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut high = at_coordinate(upper, &mut evaluations)?;
        while high.time_s < input.duration_s && upper < maximum {
            let next = (2.0 * upper).min(maximum);
            if next <= upper {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            upper = next;
            high = at_coordinate(upper, &mut evaluations)?;
        }
        if upper == maximum && high.time_s <= input.duration_s {
            return match coordinate {
                Coordinate::Distance { .. } => {
                    wrap(self.finish(input, high.time_s, boundary_kind, high, evaluations, 0)?)
                }
                Coordinate::Equilibrium { .. } => Err(CarreauMeniscusError::Unrepresentable),
            };
        }
        if (high.time_s - input.duration_s).abs()
            <= controls.relative_time_tolerance * input.duration_s
        {
            return wrap(self.finish(
                input,
                input.duration_s,
                CarreauMeniscusEnd::Duration,
                high,
                evaluations,
                0,
            )?);
        }
        let mut lower = 0.0;
        let mut value = upper * (input.duration_s / high.time_s);
        for iteration in 1..=controls.maximum_iterations {
            let trial = at_coordinate(value, &mut evaluations)?;
            let residual = trial.time_s - input.duration_s;
            if residual.abs() <= controls.relative_time_tolerance * input.duration_s {
                return wrap(self.finish(
                    input,
                    input.duration_s,
                    CarreauMeniscusEnd::Duration,
                    trial,
                    evaluations,
                    iteration,
                )?);
            }
            if residual < 0.0 {
                lower = value;
            } else {
                upper = value;
            }
            let (_, _, length_jacobian) = coordinate.point(
                value,
                input.wet_length_m,
                direction,
                pressure0,
                drive0,
                gradient,
            )?;
            let derivative = self.integrands(&PathSample {
                profile: trial.endpoint,
                length_jacobian,
            })?[0];
            let newton = value - residual / derivative;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
            if next == value || next <= 0.0 {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            value = next;
        }
        Err(CarreauMeniscusError::IterationBudget)
    }
}
