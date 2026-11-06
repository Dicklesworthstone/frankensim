//! Horizontal and inclined displacement with both original Carreau fluids.
//!
//! The two fully developed pipe columns carry one volume flux. At each
//! front position, a bounded pressure solve resolves the sum of their actual
//! radial constitutive drops. Time and the two accepted viscous powers are
//! integrated together along the interface path, preserving each phase's heat.
//! The positive minimum column lengths are explicit fully developed operating
//! limits; this owner does not invent an entrance law when a column disappears.

use super::inclined::Coordinate;

use super::{
    CarreauCapillaryTube, CarreauMeniscusControls, CarreauMeniscusError, CarreauMeniscusInput,
    CarreauPipeProfile, LaminarCarreauPipe, MaterialError, integrate_positive,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauTwoFluidControls {
    pub meniscus: CarreauMeniscusControls,
    /// Relative series-pressure residual AND relative flux Newton correction.
    pub relative_pressure_tolerance: f64,
    /// Maximum pressure-sum corrections per spatial sample, in 1..=128.
    pub maximum_pressure_iterations: u32,
}

impl Default for CarreauTwoFluidControls {
    fn default() -> Self {
        Self {
            meniscus: CarreauMeniscusControls::default(),
            relative_pressure_tolerance: 1e-9,
            maximum_pressure_iterations: 64,
        }
    }
}

impl CarreauTwoFluidControls {
    fn admit(self) -> Result<(), CarreauMeniscusError> {
        self.meniscus.admit()?;
        if !self.relative_pressure_tolerance.is_finite()
            || self.relative_pressure_tolerance <= 0.0
            || self.relative_pressure_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_pressure_iterations)
        {
            return Err(CarreauMeniscusError::Input(
                "invalid two-fluid pressure controls",
            ));
        }
        Ok(())
    }
}

/// Two immiscible incompressible fluids in one straight circular tube.
/// Radius and total length must agree exactly; each law retains its own frozen
/// temperature and original physical wall-rate support. The capillary pressure
/// and contact angle on `invading` must describe this actual pair and substrate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauTwoFluidCapillary {
    invading: CarreauCapillaryTube,
    displaced: LaminarCarreauPipe,
    maximum_wet_length_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarreauTwoFluidEnd {
    Duration,
    Stationary,
    MinimumInvadingLength,
    MinimumDisplacedLength,
}

/// Coupled accepted radial states, in invading/displaced order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauTwoFluidProfile {
    pub wet_length_m: f64,
    /// Common signed transport flux used to advance the interface [m³/s].
    pub volumetric_flow_m3_s: f64,
    pub phase_profiles: [CarreauPipeProfile; 2],
    /// Actual original profile flow minus the common transport flux [m³/s].
    pub phase_flow_residuals_m3_s: [f64; 2],
    /// Sum of accepted phase pressures minus imposed total drive [Pa].
    pub pressure_residual_pa: f64,
}

/// Physical phase heat and front motion, with unresolved work kept separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauTwoFluidStep {
    pub input: CarreauMeniscusInput,
    pub wet_length_m: f64,
    pub elapsed_s: f64,
    pub end: CarreauTwoFluidEnd,
    /// Invading inventory gains this volume; displaced inventory loses it.
    /// The same signed volume enters/leaves the prescribed-pressure reservoirs.
    pub transferred_volume_m3: f64,
    pub reservoir_work_j: f64,
    pub wall_energy_change_j: f64,
    /// Each original profile's dissipated power integrated over actual time.
    pub phase_viscous_dissipation_j: [f64; 2],
    pub viscous_dissipation_j: f64,
    /// Pressure/wall work minus both actual heats, including numerical pressure,
    /// flux-continuity and quadrature defects; never extra physical heat.
    pub constitutive_work_defect_j: f64,
    pub travel_time_residual_s: f64,
    pub quadrature_relative_change: f64,
    /// Original pipe solve calls, including source-envelope preparation and
    /// rejected pressure/time iterates. Each call has its own fixed work cap.
    pub profile_evaluations: u32,
    pub iterations: u32,
    pub endpoint: CarreauTwoFluidProfile,
}

/// Two-phase inventory and gravity accounting for an inclined straight tube.
/// Potential energy uses the inlet elevation as zero. Reservoir transfer is
/// signed outward: displaced liquid exits at the outlet elevation during
/// invasion, and enters there during drainage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauInclinedTwoFluidStep {
    pub meniscus: CarreauTwoFluidStep,
    pub phase_densities_kg_m3: [f64; 2],
    /// Positive acceleration opposes increasing wet length.
    pub opposing_gravity_m_s2: f64,
    /// Signed changes of invading/displaced mass stored inside the tube.
    pub phase_mass_changes_kg: [f64; 2],
    pub hydrostatic_pressure_pa: f64,
    pub stored_potential_energy_change_j: f64,
    pub reservoir_potential_energy_out_j: f64,
    /// Stored change plus the net potential energy transported to reservoirs.
    pub gravitational_work_j: f64,
}

#[derive(Clone, Copy)]
struct Envelope {
    flow_ceiling: f64,
    /// Zero-flow differential resistance per unit length [Pa s/m⁴].
    rest_resistance_gradient: [f64; 2],
}

struct EvaluationWork<'a> {
    count: u32,
    check: &'a mut dyn FnMut() -> Result<(), CarreauMeniscusError>,
}

#[derive(Clone, Copy)]
struct TwoFluidTravel {
    integrals: [f64; 3],
    relative_change: f64,
    endpoint: CarreauTwoFluidProfile,
}

fn count(
    controls: CarreauTwoFluidControls,
    evaluations: &mut EvaluationWork<'_>,
) -> Result<(), CarreauMeniscusError> {
    (evaluations.check)()?;
    if evaluations.count >= controls.meniscus.maximum_profile_evaluations {
        return Err(CarreauMeniscusError::ProfileBudget);
    }
    evaluations.count += 1;
    Ok(())
}

fn material_state(message: &str) -> CarreauMeniscusError {
    CarreauMeniscusError::Material(MaterialError::State {
        what: message.into(),
    })
}

impl CarreauTwoFluidCapillary {
    pub fn new(
        invading: CarreauCapillaryTube,
        displaced: LaminarCarreauPipe,
        minimum_displaced_length_m: f64,
    ) -> Result<Self, CarreauMeniscusError> {
        let length = invading.pipe.length_m();
        let maximum_wet_length_m = length - minimum_displaced_length_m;
        if displaced.radius_m() != invading.pipe.radius_m()
            || displaced.length_m() != length
            || !minimum_displaced_length_m.is_finite()
            || minimum_displaced_length_m <= 0.0
            || !maximum_wet_length_m.is_finite()
            || maximum_wet_length_m < invading.minimum_wet_length_m
            || maximum_wet_length_m >= length
        {
            return Err(CarreauMeniscusError::Input(
                "two-fluid geometry or wet interval disagrees",
            ));
        }
        Ok(Self {
            invading,
            displaced,
            maximum_wet_length_m,
        })
    }

    pub const fn invading(&self) -> &CarreauCapillaryTube {
        &self.invading
    }
    pub const fn displaced(&self) -> &LaminarCarreauPipe {
        &self.displaced
    }
    pub const fn maximum_wet_length_m(&self) -> f64 {
        self.maximum_wet_length_m
    }

    fn phases(&self) -> [LaminarCarreauPipe; 2] {
        [self.invading.pipe, self.displaced]
    }

    fn at_flow(
        &self,
        phase: LaminarCarreauPipe,
        length: f64,
        flow: f64,
        controls: CarreauTwoFluidControls,
        evaluations: &mut EvaluationWork<'_>,
    ) -> Result<CarreauPipeProfile, CarreauMeniscusError> {
        count(controls, evaluations)?;
        let pipe = controls.meniscus.pipe;
        Ok(*phase
            .with_length(length)?
            .profile_at_flow(
                flow,
                pipe,
                pipe.relative_constitutive_tolerance,
                pipe.maximum_constitutive_iterations,
            )?
            .profile())
    }

    fn envelope(
        &self,
        boundary_drives: [f64; 2],
        controls: CarreauTwoFluidControls,
        evaluations: &mut EvaluationWork<'_>,
    ) -> Result<Envelope, CarreauMeniscusError> {
        let phases = self.phases();
        let mut flow_ceiling = f64::INFINITY;
        for phase in phases {
            count(controls, evaluations)?;
            let cap = phase.profile_at_wall_rate(
                phase.maximum_laminar_wall_rate_per_s(),
                controls.meniscus.pipe,
            )?;
            flow_ceiling = flow_ceiling.min(cap.flow().volumetric_flow_m3_s);
        }
        if !flow_ceiling.is_finite() || flow_ceiling <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut ceiling_gradient = [0.0; 2];
        let mut rest_resistance_gradient = [0.0; 2];
        for (i, phase) in phases.into_iter().enumerate() {
            let cap = self.at_flow(phase, phase.length_m(), flow_ceiling, controls, evaluations)?;
            ceiling_gradient[i] = cap.flow().pressure_drop_pa / phase.length_m();
            let rest = self.at_flow(phase, phase.length_m(), 0.0, controls, evaluations)?;
            rest_resistance_gradient[i] =
                1.0 / (rest.pressure_conductance_m3_s_pa()? * phase.length_m());
        }
        if ceiling_gradient
            .iter()
            .chain(&rest_resistance_gradient)
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        // At fixed common flux, each original pressure is proportional to its
        // column length. Thus the limiting total pressure is affine in x;
        // With gravity the drive is affine too. The absolute drive is convex;
        // checking both endpoints admits the whole declared displacement range.
        for (x, drive) in [
            self.invading.minimum_wet_length_m,
            self.maximum_wet_length_m,
        ]
        .into_iter()
        .zip(boundary_drives)
        {
            let available =
                ceiling_gradient[0] * x + ceiling_gradient[1] * (self.displaced.length_m() - x);
            if !available.is_finite() {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            if drive.abs() > available {
                return Err(material_state(
                    "two-fluid drive exceeds common-flow source support over the wet interval",
                ));
            }
        }
        Ok(Envelope {
            flow_ceiling,
            rest_resistance_gradient,
        })
    }

    fn profile(
        &self,
        x: f64,
        drive: f64,
        envelope: Envelope,
        controls: CarreauTwoFluidControls,
        evaluations: &mut EvaluationWork<'_>,
    ) -> Result<CarreauTwoFluidProfile, CarreauMeniscusError> {
        let phases = self.phases();
        let lengths = [x, self.displaced.length_m() - x];
        let target = drive.abs();
        let resistance = envelope.rest_resistance_gradient[0] * lengths[0]
            + envelope.rest_resistance_gradient[1] * lengths[1];
        let mut q = (target / resistance).min(envelope.flow_ceiling);
        if !resistance.is_finite()
            || resistance <= 0.0
            || !q.is_finite()
            || target > 0.0 && q <= 0.0
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let (mut low, mut high) = (0.0, envelope.flow_ceiling);
        for _ in 0..controls.maximum_pressure_iterations {
            let signed_flow = q.copysign(drive);
            let profiles = [
                self.at_flow(phases[0], lengths[0], signed_flow, controls, evaluations)?,
                self.at_flow(phases[1], lengths[1], signed_flow, controls, evaluations)?,
            ];
            let pressure =
                profiles[0].flow().pressure_drop_pa + profiles[1].flow().pressure_drop_pa;
            let residual = pressure - drive;
            let derivative = 1.0 / profiles[0].pressure_conductance_m3_s_pa()?
                + 1.0 / profiles[1].pressure_conductance_m3_s_pa()?;
            let correction = (pressure.abs() - target) / derivative;
            if ![pressure, residual, derivative, correction]
                .iter()
                .all(|v| v.is_finite())
                || derivative <= 0.0
            {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            if target == 0.0
                || residual.abs() <= controls.relative_pressure_tolerance * target
                    && correction.abs() <= controls.relative_pressure_tolerance * q
            {
                return Ok(CarreauTwoFluidProfile {
                    wet_length_m: x,
                    volumetric_flow_m3_s: signed_flow,
                    phase_flow_residuals_m3_s: profiles
                        .map(|p| p.flow().volumetric_flow_m3_s - signed_flow),
                    phase_profiles: profiles,
                    pressure_residual_pa: residual,
                });
            }
            if pressure.abs() < target {
                low = q;
            } else {
                high = q;
            }
            let candidate = q - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                low + 0.5 * (high - low)
            };
            if next == q || next <= 0.0 {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            q = next;
        }
        Err(CarreauMeniscusError::SeriesPressureBudget)
    }

    fn integrands(
        &self,
        profile: &CarreauTwoFluidProfile,
    ) -> Result<[f64; 3], CarreauMeniscusError> {
        let inverse_speed = self.invading.area_m2 / profile.volumetric_flow_m3_s.abs();
        let values = [
            inverse_speed,
            inverse_speed * profile.phase_profiles[0].flow().dissipated_power_w,
            inverse_speed * profile.phase_profiles[1].flow().dissipated_power_w,
        ];
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        Ok(values)
    }

    /// Advance horizontal two-fluid displacement under held reservoir pressure.
    /// Both phase profiles contribute to the resistance at EVERY spatial sample.
    /// Heat integrates each original profile's power over the resolved travel
    /// time; it is not apportioned using one endpoint viscosity or pressure.
    ///
    /// Temperatures are frozen independently. There is no interphase mixing,
    /// compressibility, gravity, inertia, deposited film, dynamic contact angle
    /// or thermal evolution. A gas may be used only under the caller's explicit
    /// uniform-property incompressible approximation, with vented reservoirs.
    /// Checked positive column limits avoid an unmodeled entrance singularity.
    /// Immutable solves allocate no per-sample storage and publish no partial step.
    pub fn step_horizontal(
        &self,
        input: CarreauMeniscusInput,
        controls: CarreauTwoFluidControls,
    ) -> Result<CarreauTwoFluidStep, CarreauMeniscusError> {
        self.step_horizontal_with_check(input, controls, || Ok(()))
    }

    /// Same physical solve with a caller-owned check before each original
    /// bounded pipe solve and before final publication. No state is changed on
    /// cancellation; an identical input can be retried without recovery work.
    pub fn step_horizontal_with_check(
        &self,
        input: CarreauMeniscusInput,
        controls: CarreauTwoFluidControls,
        check: impl FnMut() -> Result<(), CarreauMeniscusError>,
    ) -> Result<CarreauTwoFluidStep, CarreauMeniscusError> {
        let pressure0 = input.applied_pressure_pa + self.invading.capillary_pressure_pa;
        self.step_with_affine_load(input, pressure0, 0.0, controls, check)
    }

    /// Advance both original fluid columns with explicit constant densities
    /// and signed axial gravity. Reservoir pressures are measured at the two
    /// tube ends, so BOTH columns contribute to hydrostatic resistance.
    ///
    /// A stable density-contrast equilibrium is approached logarithmically,
    /// never crossed or reached in finite time. The original source support is
    /// checked over the entire authored wet interval. Source temperatures stay
    /// fixed; no buoyant mixing, inertia or variable density is inferred.
    pub fn step_inclined(
        &self,
        input: CarreauMeniscusInput,
        phase_densities_kg_m3: [f64; 2],
        opposing_gravity_m_s2: f64,
        controls: CarreauTwoFluidControls,
    ) -> Result<CarreauInclinedTwoFluidStep, CarreauMeniscusError> {
        self.step_inclined_with_check(
            input,
            phase_densities_kg_m3,
            opposing_gravity_m_s2,
            controls,
            || Ok(()),
        )
    }

    /// Checked inclined advance; a refused step leaves the immutable owner
    /// available for an identical retry, including late cancellation.
    pub fn step_inclined_with_check(
        &self,
        input: CarreauMeniscusInput,
        phase_densities_kg_m3: [f64; 2],
        opposing_gravity_m_s2: f64,
        controls: CarreauTwoFluidControls,
        mut check: impl FnMut() -> Result<(), CarreauMeniscusError>,
    ) -> Result<CarreauInclinedTwoFluidStep, CarreauMeniscusError> {
        if phase_densities_kg_m3
            .iter()
            .any(|rho| !rho.is_finite() || *rho <= 0.0)
            || !opposing_gravity_m_s2.is_finite()
        {
            return Err(CarreauMeniscusError::Input(
                "invalid two-fluid gravity inputs",
            ));
        }
        let [rho0, rho1] = phase_densities_kg_m3;
        let length = self.displaced.length_m();
        let pressure = input.applied_pressure_pa + self.invading.capillary_pressure_pa;
        let gradient = opposing_gravity_m_s2 * (rho0 - rho1);
        let pressure0 = (-opposing_gravity_m_s2 * rho1).mul_add(length, pressure);
        let mut meniscus = if opposing_gravity_m_s2 == 0.0 {
            self.step_horizontal_with_check(input, controls, &mut check)?
        } else {
            self.step_with_affine_load(input, pressure0, gradient, controls, &mut check)?
        };
        let mean_length = 0.5 * input.wet_length_m + 0.5 * meniscus.wet_length_m;
        let volume = meniscus.transferred_volume_m3;
        let phase_mass_changes_kg = [rho0 * volume, -rho1 * volume];
        let hydrostatic_pressure_pa = opposing_gravity_m_s2
            * (rho0 * meniscus.wet_length_m + rho1 * (length - meniscus.wet_length_m));
        let stored_potential_energy_change_j = volume * (gradient * mean_length);
        let reservoir_potential_energy_out_j = volume * ((opposing_gravity_m_s2 * rho1) * length);
        // Sum column weights before multiplying by transferred volume. This
        // avoids cancellation between stored and advected potential energy.
        let gravitational_work_j =
            volume * (opposing_gravity_m_s2 * (rho0 * mean_length + rho1 * (length - mean_length)));
        meniscus.constitutive_work_defect_j = meniscus.reservoir_work_j
            - meniscus.wall_energy_change_j
            - gravitational_work_j
            - meniscus.viscous_dissipation_j;
        if !phase_mass_changes_kg
            .iter()
            .chain(&[
                hydrostatic_pressure_pa,
                stored_potential_energy_change_j,
                reservoir_potential_energy_out_j,
                gravitational_work_j,
                meniscus.constitutive_work_defect_j,
            ])
            .all(|v| v.is_finite())
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        check()?;
        Ok(CarreauInclinedTwoFluidStep {
            meniscus,
            phase_densities_kg_m3,
            opposing_gravity_m_s2,
            phase_mass_changes_kg,
            hydrostatic_pressure_pa,
            stored_potential_energy_change_j,
            reservoir_potential_energy_out_j,
            gravitational_work_j,
        })
    }

    fn step_with_affine_load(
        &self,
        input: CarreauMeniscusInput,
        pressure0: f64,
        gradient: f64,
        controls: CarreauTwoFluidControls,
        mut check: impl FnMut() -> Result<(), CarreauMeniscusError>,
    ) -> Result<CarreauTwoFluidStep, CarreauMeniscusError> {
        check()?;
        controls.admit()?;
        if !input.wet_length_m.is_finite()
            || !(self.invading.minimum_wet_length_m..=self.maximum_wet_length_m)
                .contains(&input.wet_length_m)
            || !input.applied_pressure_pa.is_finite()
            || !input.duration_s.is_finite()
            || input.duration_s < 0.0
        {
            return Err(CarreauMeniscusError::Input(
                "invalid two-fluid meniscus state, pressure or duration",
            ));
        }
        let drive = (-gradient).mul_add(input.wet_length_m, pressure0);
        let boundary_drives = [
            self.invading.minimum_wet_length_m,
            self.maximum_wet_length_m,
        ]
        .map(|x| (-gradient).mul_add(x, pressure0));
        if ![drive, pressure0, gradient]
            .iter()
            .chain(&boundary_drives)
            .all(|v| v.is_finite())
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let mut evaluations = EvaluationWork {
            count: 0,
            check: &mut check,
        };
        let envelope = self.envelope(boundary_drives, controls, &mut evaluations)?;
        let start = self.profile(
            input.wet_length_m,
            drive,
            envelope,
            controls,
            &mut evaluations,
        )?;
        if drive == 0.0 || input.duration_s == 0.0 {
            let end = if drive == 0.0 {
                CarreauTwoFluidEnd::Stationary
            } else {
                CarreauTwoFluidEnd::Duration
            };
            return self.finish(
                input,
                input.duration_s,
                end,
                TwoFluidTravel {
                    integrals: [input.duration_s, 0.0, 0.0],
                    relative_change: 0.0,
                    endpoint: start,
                },
                evaluations,
                0,
            );
        }
        let direction = drive.signum();
        let boundary = if direction > 0.0 {
            self.maximum_wet_length_m
        } else {
            self.invading.minimum_wet_length_m
        };
        let end = if direction > 0.0 {
            CarreauTwoFluidEnd::MinimumDisplacedLength
        } else {
            CarreauTwoFluidEnd::MinimumInvadingLength
        };
        let maximum = (boundary - input.wet_length_m).abs();
        if maximum == 0.0 {
            return self.finish(
                input,
                0.0,
                end,
                TwoFluidTravel {
                    integrals: [0.0; 3],
                    relative_change: 0.0,
                    endpoint: start,
                },
                evaluations,
                0,
            );
        }
        let coordinate = Coordinate::new(input.wet_length_m, drive, gradient, boundary)?;
        let make_sample = |value: f64, evaluations: &mut EvaluationWork<'_>| {
            let (x, local_drive, jacobian) = coordinate.point(
                value,
                input.wet_length_m,
                direction,
                pressure0,
                drive,
                gradient,
            )?;
            let profile = if value == 0.0 {
                start
            } else {
                self.profile(x, local_drive, envelope, controls, evaluations)?
            };
            Ok::<_, CarreauMeniscusError>((profile, jacobian))
        };
        let (_, start_jacobian) = make_sample(0.0, &mut evaluations)?;
        let at_coordinate = |value: f64, evaluations: &mut EvaluationWork<'_>| {
            let (endpoint, jacobian) = make_sample(value, evaluations)?;
            if endpoint.wet_length_m == input.wet_length_m {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            if let Coordinate::Equilibrium { distance, .. } = coordinate {
                // Integrate time in log distance, and each actual source heat
                // in physical distance. Heat decays exponentially in the log
                // coordinate, making a shared time/heat grid inefficient.
                let ([time], time_change) = integrate_positive(
                    value,
                    [
                        [self.integrands(&start)?[0] * start_jacobian],
                        [self.integrands(&endpoint)?[0] * jacobian],
                    ],
                    controls.meniscus,
                    |fraction| {
                        let (profile, jacobian) = make_sample(value * fraction, evaluations)?;
                        Ok([self.integrands(&profile)?[0] * jacobian])
                    },
                )?;
                let wet_distance = distance * -fs_math::det::expm1(-value);
                let decay = fs_math::det::exp(-value);
                let heat_density = |profile: &CarreauTwoFluidProfile| {
                    let integrands = self.integrands(profile)?;
                    Ok::<_, CarreauMeniscusError>([integrands[1], integrands[2]])
                };
                let (heat, heat_change) = integrate_positive(
                    wet_distance,
                    [heat_density(&start)?, heat_density(&endpoint)?],
                    controls.meniscus,
                    |fraction| {
                        let x = input.wet_length_m + direction * wet_distance * fraction;
                        // Resolve the small remaining drive without subtracting
                        // nearly equal reservoir and hydrostatic pressures.
                        let local_drive = drive * ((1.0 - fraction) + fraction * decay);
                        heat_density(&self.profile(
                            x,
                            local_drive,
                            envelope,
                            controls,
                            evaluations,
                        )?)
                    },
                )?;
                return Ok(TwoFluidTravel {
                    integrals: [time, heat[0], heat[1]],
                    relative_change: time_change.max(heat_change),
                    endpoint,
                });
            }
            let (integrals, relative_change) = integrate_positive(
                value,
                [self.integrands(&start)?, self.integrands(&endpoint)?],
                controls.meniscus,
                |fraction| self.integrands(&make_sample(value * fraction, evaluations)?.0),
            )?;
            Ok(TwoFluidTravel {
                integrals,
                relative_change,
                endpoint,
            })
        };
        // A local physical speed estimate keeps small steps local. Bracket
        // growth resolves either ordering of the two nonlinear resistance laws.
        let speed = (start.volumetric_flow_m3_s.abs() / self.invading.area_m2) / start_jacobian;
        if !speed.is_finite() || speed <= 0.0 {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        let maximum = coordinate.maximum();
        let mut upper = if input.duration_s >= maximum / speed {
            maximum
        } else {
            input.duration_s * speed
        };
        let mut high = at_coordinate(upper, &mut evaluations)?;
        while high.integrals[0] < input.duration_s && upper < maximum {
            let next = (2.0 * upper).min(maximum);
            if next <= upper {
                return Err(CarreauMeniscusError::Unrepresentable);
            }
            upper = next;
            high = at_coordinate(upper, &mut evaluations)?;
        }
        if upper == maximum && high.integrals[0] <= input.duration_s {
            return match coordinate {
                Coordinate::Distance { .. } => {
                    self.finish(input, high.integrals[0], end, high, evaluations, 0)
                }
                Coordinate::Equilibrium { .. } => Err(CarreauMeniscusError::Unrepresentable),
            };
        }
        if (high.integrals[0] - input.duration_s).abs()
            <= controls.meniscus.relative_time_tolerance * input.duration_s
        {
            return self.finish(
                input,
                input.duration_s,
                CarreauTwoFluidEnd::Duration,
                high,
                evaluations,
                0,
            );
        }
        let mut lower = 0.0;
        let mut distance = upper * (input.duration_s / high.integrals[0]);
        for iteration in 1..=controls.meniscus.maximum_iterations {
            let trial = at_coordinate(distance, &mut evaluations)?;
            let residual = trial.integrals[0] - input.duration_s;
            if residual.abs() <= controls.meniscus.relative_time_tolerance * input.duration_s {
                return self.finish(
                    input,
                    input.duration_s,
                    CarreauTwoFluidEnd::Duration,
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
            let (_, _, jacobian) = coordinate.point(
                distance,
                input.wet_length_m,
                direction,
                pressure0,
                drive,
                gradient,
            )?;
            let endpoint_speed =
                (trial.endpoint.volumetric_flow_m3_s.abs() / self.invading.area_m2) / jacobian;
            let candidate = distance - residual * endpoint_speed;
            let next = if candidate.is_finite() && candidate > lower && candidate < upper {
                candidate
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
        end: CarreauTwoFluidEnd,
        travel: TwoFluidTravel,
        evaluations: EvaluationWork<'_>,
        iterations: u32,
    ) -> Result<CarreauTwoFluidStep, CarreauMeniscusError> {
        let wet_length_m = travel.endpoint.wet_length_m;
        let transferred_volume_m3 = self.invading.area_m2 * (wet_length_m - input.wet_length_m);
        let reservoir_work_j = input.applied_pressure_pa * transferred_volume_m3;
        let wall_energy_change_j = -self.invading.capillary_pressure_pa * transferred_volume_m3;
        let phase_viscous_dissipation_j = [travel.integrals[1], travel.integrals[2]];
        let viscous_dissipation_j = travel.integrals[1] + travel.integrals[2];
        let constitutive_work_defect_j =
            reservoir_work_j - wall_energy_change_j - viscous_dissipation_j;
        if ![
            transferred_volume_m3,
            reservoir_work_j,
            wall_energy_change_j,
            viscous_dissipation_j,
            constitutive_work_defect_j,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(CarreauMeniscusError::Unrepresentable);
        }
        (evaluations.check)()?;
        Ok(CarreauTwoFluidStep {
            input,
            wet_length_m,
            elapsed_s,
            end,
            transferred_volume_m3,
            reservoir_work_j,
            wall_energy_change_j,
            phase_viscous_dissipation_j,
            viscous_dissipation_j,
            constitutive_work_defect_j,
            travel_time_residual_s: travel.integrals[0] - elapsed_s,
            quadrature_relative_change: travel.relative_change,
            profile_evaluations: evaluations.count,
            iterations,
            endpoint: travel.endpoint,
        })
    }
}
