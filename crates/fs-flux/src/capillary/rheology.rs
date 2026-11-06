//! Capillary filling and drainage with the original Carreau law.
//!
//! At liquid length x, the fully developed circular-pipe owner resolves
//! Q(x) from dp = p_in - p_out + 2 sigma cos(theta)/R. The moving front obeys
//! A dx/dt = Q(x). We invert its travel-time integral A integral(dx/Q),
//! evaluating the radial shear-dependent viscosity at every spatial sample.
//! The same positive quadrature integrates the accepted viscous work.
//!
//! The prescribed pressure and temperature are frozen throughout a step.
//! A strictly positive minimum wet length bounds the declared developed-flow
//! interval; there is no dry-entry singularity, entrance model, inertia,
//! evolving contact angle, displaced-fluid resistance or thermal feedback.
//! The horizontal owner excludes gravity; [`inclined`] adds its resolved
//! hydrostatic load and potential energy. A step admits its load over the
//! entire declared wet-length interval.
//! This numerical consumer adds no empirical authority to the supplied law.

pub mod inclined;

use core::f64::consts::PI;
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::pipe::{
    CarreauPipeControls, CarreauPipeProfile, LaminarCarreauPipe,
};

/// Hard work bounds and estimated numerical tolerances, distinct from
/// uncertainty in the authored viscosity and wetting data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauMeniscusControls {
    pub pipe: CarreauPipeControls,
    /// Relative travel-time inversion residual and positive-integral changes.
    pub relative_time_tolerance: f64,
    /// Maximum safeguarded Newton corrections, in 1..=128.
    pub maximum_iterations: u32,
    /// Maximum Simpson panel doublings per travel integral, in 3..=16.
    pub maximum_quadrature_refinements: u32,
    /// Total original pipe-profile queries across the entire step.
    pub maximum_profile_evaluations: u32,
}

impl Default for CarreauMeniscusControls {
    fn default() -> Self {
        Self {
            pipe: CarreauPipeControls::default(),
            relative_time_tolerance: 1e-8,
            maximum_iterations: 64,
            maximum_quadrature_refinements: 12,
            maximum_profile_evaluations: 8192,
        }
    }
}

impl CarreauMeniscusControls {
    fn admit(self) -> Result<(), CarreauMeniscusError> {
        if !self.relative_time_tolerance.is_finite()
            || self.relative_time_tolerance <= 0.0
            || self.relative_time_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_iterations)
            || !(3..=16).contains(&self.maximum_quadrature_refinements)
            || !(1..=1_000_000).contains(&self.maximum_profile_evaluations)
        {
            return Err(CarreauMeniscusError::Input("invalid numerical controls"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CarreauMeniscusError {
    Input(&'static str),
    /// The original constitutive support or pipe solver refused the load.
    Material(MaterialError),
    ProfileBudget,
    QuadratureBudget,
    IterationBudget,
    Unrepresentable,
}

impl From<MaterialError> for CarreauMeniscusError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for CarreauMeniscusError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Carreau capillary: {self:?}")
    }
}
impl std::error::Error for CarreauMeniscusError {}

/// Original rheology/temperature/radius and declared total tube length.
/// Surface tension and equilibrium angle are explicit SI inputs; a material
/// name does not supply them. The model is immutable across tentative solves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauCapillaryTube {
    pipe: LaminarCarreauPipe,
    minimum_wet_length_m: f64,
    surface_tension_n_m: f64,
    contact_angle_radians: f64,
    capillary_pressure_pa: f64,
    area_m2: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauMeniscusInput {
    pub wet_length_m: f64,
    /// Inlet minus uniform vented downstream pressure [Pa].
    pub applied_pressure_pa: f64,
    pub duration_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarreauMeniscusEnd {
    Duration,
    Stationary,
    Full,
    MinimumWetLength,
}

/// Accepted travel and physical work, with numerical defects kept separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauMeniscusStep {
    pub input: CarreauMeniscusInput,
    pub wet_length_m: f64,
    /// An endpoint event consumes only its arrival time.
    pub elapsed_s: f64,
    pub end: CarreauMeniscusEnd,
    pub transferred_volume_m3: f64,
    pub reservoir_work_j: f64,
    /// Change in wall wetting energy, -2 pi R sigma cos(theta) delta_x.
    pub wall_energy_change_j: f64,
    /// Change in gravitational potential energy; zero for horizontal travel.
    pub gravitational_energy_change_j: f64,
    /// Integral of the ORIGINAL accepted pipe pressure times signed dV.
    pub viscous_dissipation_j: f64,
    /// Reservoir work minus wall and gravitational energy changes minus heat.
    /// Constitutive pressure residuals and quadrature error are not extra heat.
    pub constitutive_work_defect_j: f64,
    /// Evaluated travel time minus elapsed_s, not a physical time correction.
    pub travel_time_residual_s: f64,
    pub quadrature_relative_change: f64,
    pub profile_evaluations: u32,
    pub iterations: u32,
    /// Actual endpoint liquid length and original radial constitutive profile.
    /// At an end event this is the approach flow; post-event flow is not evolved.
    pub endpoint_profile: CarreauPipeProfile,
}

#[derive(Clone, Copy)]
struct Travel {
    time_s: f64,
    heat_j: f64,
    relative_change: f64,
    endpoint: CarreauPipeProfile,
}

#[derive(Clone, Copy)]
struct PathSample {
    profile: CarreauPipeProfile,
    /// Absolute physical length change per integration coordinate.
    length_jacobian: f64,
}

impl CarreauCapillaryTube {
    pub fn new(
        pipe: LaminarCarreauPipe,
        minimum_wet_length_m: f64,
        surface_tension_n_m: f64,
        contact_angle_radians: f64,
    ) -> Result<Self, CarreauMeniscusError> {
        if !minimum_wet_length_m.is_finite()
            || minimum_wet_length_m <= 0.0
            || minimum_wet_length_m > pipe.length_m()
            || !surface_tension_n_m.is_finite()
            || surface_tension_n_m < 0.0
            || !contact_angle_radians.is_finite()
            || !(0.0..=PI).contains(&contact_angle_radians)
        {
            return Err(CarreauMeniscusError::Input(
                "invalid wet interval or wetting inputs",
            ));
        }
        let capillary_pressure_pa = 2.0
            * (surface_tension_n_m * fs_math::det::cos(contact_angle_radians))
            / pipe.radius_m();
        let area_m2 = (PI * pipe.radius_m()) * pipe.radius_m();
        if !capillary_pressure_pa.is_finite() || !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        Ok(Self {
            pipe,
            minimum_wet_length_m,
            surface_tension_n_m,
            contact_angle_radians,
            capillary_pressure_pa,
            area_m2,
        })
    }

    pub const fn pipe(&self) -> &LaminarCarreauPipe {
        &self.pipe
    }
    pub const fn minimum_wet_length_m(&self) -> f64 {
        self.minimum_wet_length_m
    }
    pub const fn surface_tension_n_m(&self) -> f64 {
        self.surface_tension_n_m
    }
    pub const fn contact_angle_radians(&self) -> f64 {
        self.contact_angle_radians
    }
    pub const fn capillary_pressure_pa(&self) -> f64 {
        self.capillary_pressure_pa
    }

    fn profile(
        &self,
        length: f64,
        pressure: f64,
        controls: CarreauMeniscusControls,
        evaluations: &mut u32,
    ) -> Result<CarreauPipeProfile, CarreauMeniscusError> {
        if *evaluations >= controls.maximum_profile_evaluations {
            return Err(CarreauMeniscusError::ProfileBudget);
        }
        *evaluations += 1;
        Ok(self
            .pipe
            .with_length(length)?
            .profile_at_pressure(pressure, controls.pipe)?)
    }

    fn integrands(&self, sample: &PathSample) -> Result<[f64; 2], CarreauMeniscusError> {
        let flow = sample.profile.flow();
        let values = [
            if sample.length_jacobian == 1.0 {
                self.area_m2 / flow.volumetric_flow_m3_s.abs()
            } else {
                self.area_m2 * (sample.length_jacobian / flow.volumetric_flow_m3_s.abs())
            },
            (self.area_m2 * flow.pressure_drop_pa.abs()) * sample.length_jacobian,
        ];
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        Ok(values)
    }

    fn travel(
        &self,
        start: &CarreauPipeProfile,
        endpoint: CarreauPipeProfile,
        pressure: f64,
        controls: CarreauMeniscusControls,
        evaluations: &mut u32,
    ) -> Result<Travel, CarreauMeniscusError> {
        let x0 = start.pipe().length_m();
        let x1 = endpoint.pipe().length_m();
        let signed_distance = x1 - x0;
        self.integrate_travel(
            PathSample {
                profile: *start,
                length_jacobian: 1.0,
            },
            PathSample {
                profile: endpoint,
                length_jacobian: 1.0,
            },
            signed_distance.abs(),
            controls,
            |fraction| {
                Ok(PathSample {
                    profile: self.profile(
                        x0 + signed_distance * fraction,
                        pressure,
                        controls,
                        evaluations,
                    )?,
                    length_jacobian: 1.0,
                })
            },
        )
    }

    /// Shared positive quadrature in physical distance or a regularizing
    /// coordinate. The original profile owns the accepted work at each sample.
    fn integrate_travel(
        &self,
        start: PathSample,
        endpoint: PathSample,
        span: f64,
        controls: CarreauMeniscusControls,
        mut sample: impl FnMut(f64) -> Result<PathSample, CarreauMeniscusError>,
    ) -> Result<Travel, CarreauMeniscusError> {
        if !span.is_finite() || span <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let a = self.integrands(&start)?;
        let b = self.integrands(&endpoint)?;
        let mut trapezoid = core::array::from_fn::<_, 2, _>(|k| 0.5 * a[k] + 0.5 * b[k]);
        let mut previous: Option<[f64; 2]> = None;
        let mut consecutive = 0;
        let mut panels = 1_u32;
        for _ in 0..controls.maximum_quadrature_refinements {
            let mut sums = [0.0; 2];
            let mut compensation = [0.0; 2];
            for i in 0..panels {
                let fraction = (i as f64 + 0.5) / panels as f64;
                let values = self.integrands(&sample(fraction)?)?;
                for k in 0..2 {
                    // Average before accumulation to avoid a panel-count overflow.
                    let term = values[k] / panels as f64 - compensation[k];
                    let next = sums[k] + term;
                    compensation[k] = (next - sums[k]) - term;
                    sums[k] = next;
                }
            }
            let refined = core::array::from_fn::<_, 2, _>(|k| 0.5 * trapezoid[k] + 0.5 * sums[k]);
            let integral = core::array::from_fn::<_, 2, _>(|k| {
                span * (refined[k] + (refined[k] - trapezoid[k]) / 3.0)
            });
            if integral.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            if let Some(old) = previous {
                let change = (0..2)
                    .map(|k| (integral[k] - old[k]).abs() / integral[k].max(old[k]))
                    .fold(0.0_f64, f64::max);
                consecutive = if change <= controls.relative_time_tolerance {
                    consecutive + 1
                } else {
                    0
                };
                if consecutive >= 2 {
                    return Ok(Travel {
                        time_s: integral[0],
                        heat_j: integral[1],
                        relative_change: change,
                        endpoint: endpoint.profile,
                    });
                }
            }
            previous = Some(integral);
            trapezoid = refined;
            panels *= 2;
        }
        Err(CarreauMeniscusError::QuadratureBudget)
    }

    /// Advance the complete requested interval or stop at a declared geometric
    /// endpoint. Monotone bounded Newton in distance inverts a positive
    /// travel-time quadrature; its derivative is the actual endpoint A/|Q|.
    /// No timestep-sized constitutive lag or apparent-viscosity fit is used.
    ///
    /// The fixed drive must be supported at the minimum declared wet length,
    /// which admits the whole operating interval before any result is returned.
    /// Refusal leaves this immutable model and the caller's input unchanged.
    pub fn step_horizontal(
        &self,
        input: CarreauMeniscusInput,
        controls: CarreauMeniscusControls,
    ) -> Result<CarreauMeniscusStep, CarreauMeniscusError> {
        controls.admit()?;
        if !input.wet_length_m.is_finite()
            || !(self.minimum_wet_length_m..=self.pipe.length_m()).contains(&input.wet_length_m)
            || !input.applied_pressure_pa.is_finite()
            || !input.duration_s.is_finite()
            || input.duration_s < 0.0
        {
            return Err(CarreauMeniscusError::Input(
                "invalid meniscus state, pressure or duration",
            ));
        }
        let pressure = input.applied_pressure_pa + self.capillary_pressure_pa;
        if !pressure.is_finite() {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut evaluations = 0;
        let start = self.profile(input.wet_length_m, pressure, controls, &mut evaluations)?;
        if pressure == 0.0 || input.duration_s == 0.0 {
            let end = if pressure == 0.0 {
                CarreauMeniscusEnd::Stationary
            } else {
                CarreauMeniscusEnd::Duration
            };
            return self.finish(
                input,
                input.duration_s,
                end,
                Travel {
                    time_s: input.duration_s,
                    heat_j: 0.0,
                    relative_change: 0.0,
                    endpoint: start,
                },
                evaluations,
                0,
            );
        }
        let minimum = if input.wet_length_m == self.minimum_wet_length_m {
            start
        } else {
            self.profile(
                self.minimum_wet_length_m,
                pressure,
                controls,
                &mut evaluations,
            )?
        };
        let direction = pressure.signum();
        let boundary = if pressure > 0.0 {
            self.pipe.length_m()
        } else {
            self.minimum_wet_length_m
        };
        let boundary_kind = if pressure > 0.0 {
            CarreauMeniscusEnd::Full
        } else {
            CarreauMeniscusEnd::MinimumWetLength
        };
        let maximum_distance = (boundary - input.wet_length_m).abs();
        if maximum_distance == 0.0 {
            return self.finish(
                input,
                0.0,
                boundary_kind,
                Travel {
                    time_s: 0.0,
                    heat_j: 0.0,
                    relative_change: 0.0,
                    endpoint: start,
                },
                evaluations,
                0,
            );
        }
        // |Q| decreases with length: the original stress curve is strictly
        // monotone. This upper speed bounds travel, shortening the first
        // integration interval for small requested durations.
        let peak_flow = if pressure > 0.0 {
            start.flow()
        } else {
            minimum.flow()
        }
        .volumetric_flow_m3_s
        .abs();
        let peak_speed = peak_flow / self.area_m2;
        if !peak_speed.is_finite() || peak_speed <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut upper = if input.duration_s >= maximum_distance / peak_speed {
            maximum_distance
        } else {
            input.duration_s * peak_speed
        };
        let at_distance = |distance: f64, evaluations: &mut u32| {
            let x = if distance == maximum_distance {
                boundary
            } else {
                input.wet_length_m + direction * distance
            };
            let endpoint = self.profile(x, pressure, controls, evaluations)?;
            self.travel(&start, endpoint, pressure, controls, evaluations)
        };
        let mut high = at_distance(upper, &mut evaluations)?;
        // Rounding in the speed bound may undershoot very small motions.
        // Enlarge the actual bracket; never extrapolate its constitutive law.
        while high.time_s < input.duration_s && upper < maximum_distance {
            let next = (2.0 * upper).min(maximum_distance);
            if next <= upper {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            upper = next;
            high = at_distance(upper, &mut evaluations)?;
        }
        if upper == maximum_distance && high.time_s <= input.duration_s {
            return self.finish(input, high.time_s, boundary_kind, high, evaluations, 0);
        }
        let mut lower = 0.0;
        let mut distance = upper * (input.duration_s / high.time_s);
        for iteration in 1..=controls.maximum_iterations {
            let trial = at_distance(distance, &mut evaluations)?;
            let residual = trial.time_s - input.duration_s;
            if residual.abs() <= controls.relative_time_tolerance * input.duration_s {
                return self.finish(
                    input,
                    input.duration_s,
                    CarreauMeniscusEnd::Duration,
                    trial,
                    evaluations,
                    iteration,
                );
            }
            if residual < 0.0 {
                lower = distance;
            } else {
                upper = distance;
            }
            let speed = trial.endpoint.flow().volumetric_flow_m3_s.abs() / self.area_m2;
            let newton = distance - residual * speed;
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
            if next == distance || next <= 0.0 {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            distance = next;
        }
        Err(CarreauMeniscusError::IterationBudget)
    }

    fn finish(
        &self,
        input: CarreauMeniscusInput,
        elapsed_s: f64,
        end: CarreauMeniscusEnd,
        travel: Travel,
        profile_evaluations: u32,
        iterations: u32,
    ) -> Result<CarreauMeniscusStep, CarreauMeniscusError> {
        let wet_length_m = travel.endpoint.pipe().length_m();
        let transferred_volume_m3 = self.area_m2 * (wet_length_m - input.wet_length_m);
        let reservoir_work_j = input.applied_pressure_pa * transferred_volume_m3;
        let wall_energy_change_j = -self.capillary_pressure_pa * transferred_volume_m3;
        let constitutive_work_defect_j = reservoir_work_j - wall_energy_change_j - travel.heat_j;
        if ![
            transferred_volume_m3,
            reservoir_work_j,
            wall_energy_change_j,
            constitutive_work_defect_j,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        Ok(CarreauMeniscusStep {
            input,
            wet_length_m,
            elapsed_s,
            end,
            transferred_volume_m3,
            reservoir_work_j,
            wall_energy_change_j,
            gravitational_energy_change_j: 0.0,
            viscous_dissipation_j: travel.heat_j,
            constitutive_work_defect_j,
            travel_time_residual_s: travel.time_s - elapsed_s,
            quadrature_relative_change: travel.relative_change,
            profile_evaluations,
            iterations,
            endpoint_profile: travel.endpoint,
        })
    }
}

#[cfg(test)]
mod tests;
