//! Coherent acoustic/elastic waves on both sides of a fluid-solid interface.
//!
//! The original fluid occupies z < 0 and the original isotropic solid z > 0.
//! Peak phasors use `exp(i omega (s.x - t))`, as in the elastic interface owner.
//! Normal velocity and normal traction are continuous; both solid shear
//! tractions vanish. Tangential velocity may slip. These three P/SV equations
//! share one native complex LU; SH reflects independently with zero traction.
//!
//! An incoming solid wave can therefore radiate into the fluid, including
//! coherent illumination from the fluid at the same tangential slowness.
//! Outgoing evanescence retains pressure, velocity and traction without being
//! counted as absorption. External incoming evanescent waves are not inferred.
//! The local, nondispersive, lossless small-signal model adds no viscosity,
//! radiation-reaction source model, finite aperture or material calibration.

use super::*;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

/// Channel order is fluid P, solid P, solid SV, solid SH.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FluidElasticChannel {
    Fluid,
    Solid(ElasticMode),
}

impl FluidElasticChannel {
    const fn index(self) -> usize {
        match self {
            Self::Fluid => 0,
            Self::Solid(mode) => 1 + mode.index(),
        }
    }
}

/// Explicit fluid small-signal limits, separate from the original solid's
/// strain ceiling. These caller-declared limits do not authenticate a fluid
/// closure or establish its experimental validity at the chosen amplitudes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidDriveLimits {
    maximum_relative_density: f64,
    maximum_particle_mach: f64,
}

impl FluidDriveLimits {
    pub fn try_new(
        maximum_relative_density: f64,
        maximum_particle_mach: f64,
    ) -> Result<Self, MaterialError> {
        if [maximum_relative_density, maximum_particle_mach]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0 || x >= 1.0)
        {
            return Err(invalid(
                "fluid wave limits must be finite and strictly between zero and one",
            ));
        }
        Ok(Self {
            maximum_relative_density,
            maximum_particle_mach,
        })
    }
    pub const fn maximum_relative_density(self) -> f64 {
        self.maximum_relative_density
    }
    pub const fn maximum_particle_mach(self) -> f64 {
        self.maximum_particle_mach
    }
}

/// Original fluid/EOS owner and explicitly isotropic solid. Fluid coefficients
/// are evaluated once from the same state; no air or material-name fallback.
#[derive(Debug)]
pub struct FluidElasticInterface<'a, F: FluidAcoustics + ?Sized> {
    fluid_source: &'a F,
    fluid_state: FluidAcousticState,
    solid: &'a ElasticAcousticState,
    speeds: [f64; 2],
    impedances: [f64; 2],
}

/// Four coherent channels at a fixed tangential slowness. Repeated finite
/// drives allocate nothing and reuse the prepared scattering coefficients.
#[derive(Debug)]
pub struct PreparedFluidElasticInterface<'a, F: FluidAcoustics + ?Sized> {
    fluid_source: &'a F,
    fluid_state: FluidAcousticState,
    solid: &'a ElasticAcousticState,
    tangential_slowness_s_m: f64,
    incoming: [ElasticWaveTrace; 4],
    outgoing: [ElasticWaveTrace; 4],
    scattering: [[C64; 4]; 4],
}

/// Physical interface fields and independently reconstructed normal power.
#[derive(Clone, Debug)]
pub struct FluidElasticResponse {
    pub fluid_state: FluidAcousticState,
    pub solid_state_identity: ContentHash,
    pub incoming_velocity_m_s: [C64; 4],
    pub outgoing_velocity_m_s: [C64; 4],
    pub incoming_normal_power_w_m2: [f64; 4],
    pub outgoing_normal_power_w_m2: [f64; 4],
    pub total_incoming_power_w_m2: f64,
    pub total_outgoing_power_w_m2: f64,
    /// Fluid first, solid second. Their tangential velocities need not match.
    pub boundary_velocity_m_s: [[C64; 3]; 2],
    /// Positive fluid compression equals minus the solid's normal traction.
    pub fluid_pressure_pa: C64,
    pub solid_traction_pa: [C64; 3],
    /// Fluid `Re(p conj(v_z))/2` and solid `-Re(t.conj(v))/2`, both on +z.
    pub boundary_normal_flux_w_m2: [f64; 2],
    /// Conservative bounds over each complete half-space, including phases
    /// away from the contact and decaying evanescent outgoing fields.
    pub fluid_relative_density_bound: f64,
    pub fluid_particle_mach_bound: f64,
    pub solid_strain_norm_bound: f64,
}

impl<'a, F: FluidAcoustics + ?Sized> FluidElasticInterface<'a, F> {
    pub fn new(
        fluid_source: &'a F,
        solid: &'a ElasticAcousticState,
    ) -> Result<Self, MaterialError> {
        if !solid.is_explicitly_isotropic() {
            return Err(invalid(
                "fluid-solid reflection requires an explicitly isotropic solid",
            ));
        }
        let fluid_state = fluid_source.acoustic_state()?;
        let modes = solid.directional_waves([0.0, 0.0, 1.0])?.modes;
        Ok(Self {
            fluid_source,
            fluid_state,
            solid,
            speeds: [modes[2].phase_speed_m_s, modes[0].phase_speed_m_s],
            impedances: [modes[2].impedance_pa_s_m, modes[0].impedance_pa_s_m],
        })
    }

    pub fn fluid_source(&self) -> &'a F {
        self.fluid_source
    }
    pub const fn fluid_state(&self) -> FluidAcousticState {
        self.fluid_state
    }
    pub const fn solid_state(&self) -> &'a ElasticAcousticState {
        self.solid
    }

    /// One propagating incident channel fixes conserved tangential slowness.
    /// Every other propagating channel at that slowness may illuminate it
    /// coherently. Exact grazing incidence has no incoming normal power.
    pub fn prepare_incidence(
        &self,
        channel: FluidElasticChannel,
        cosine: f64,
    ) -> Result<PreparedFluidElasticInterface<'a, F>, MaterialError> {
        self.prepare_incidence_cancellable(channel, cosine, &mut || false)
    }

    pub fn prepare_incidence_cancellable(
        &self,
        channel: FluidElasticChannel,
        cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedFluidElasticInterface<'a, F>, MaterialError> {
        poll(cancelled)?;
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(invalid("fluid-solid incident cosine must be in (0,1]"));
        }
        let incident_speed = match channel {
            FluidElasticChannel::Fluid => self.fluid_state.sound_speed_m_s(),
            FluidElasticChannel::Solid(mode) => self.speeds[mode.index().min(1)],
        };
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let tangential = sine / incident_speed;
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) {
            return Err(invalid(
                "fluid-solid tangential slowness is unrepresentable",
            ));
        }
        let mut incoming = [empty_trace(); 4];
        let mut outgoing = incoming;
        for i in 0..4 {
            let speed = if i == 0 {
                self.fluid_state.sound_speed_m_s()
            } else {
                self.speeds[(i - 1).min(1)]
            };
            let ratio = speed / incident_speed;
            let normal = normal_ratio(incident_speed, speed, cosine)?;
            let wave_cosine = normal.scale(ratio);
            let wave_sine = sine * ratio;
            let normal_slowness = divided(normal, incident_speed);
            if !finite(wave_cosine)
                || !wave_sine.is_finite()
                || (normal != C64::ZERO && normal_slowness == C64::ZERO)
            {
                return Err(invalid("fluid-solid wave direction is unrepresentable"));
            }
            for (slot, direction) in [
                (&mut outgoing[i], if i == 0 { -1.0 } else { 1.0 }),
                (&mut incoming[i], if i == 0 { 1.0 } else { -1.0 }),
            ] {
                *slot = if i == 0 {
                    fluid_trace(
                        self.fluid_state,
                        direction,
                        wave_sine,
                        wave_cosine,
                        tangential,
                        normal_slowness,
                    )?
                } else {
                    trace(
                        self.solid,
                        self.impedances,
                        i - 1,
                        direction,
                        wave_sine,
                        wave_cosine,
                        tangential,
                        normal_slowness,
                    )?
                };
            }
        }
        let mut result = PreparedFluidElasticInterface {
            fluid_source: self.fluid_source,
            fluid_state: self.fluid_state,
            solid: self.solid,
            tangential_slowness_s_m: tangential,
            incoming,
            outgoing,
            scattering: [[C64::ZERO; 4]; 4],
        };
        let mut matrix = [C64::ZERO; 9];
        let mut scales = [0.0_f64; 3];
        for row in 0..3 {
            for column in 0..3 {
                let value = contact_entry(result.outgoing[column], column, row);
                matrix[3 * row + column] = value;
                scales[row] = scales[row].max(value.abs());
            }
            if !scales[row].is_finite() || scales[row] == 0.0 {
                return Err(invalid(
                    "fluid-solid contact row is singular or unrepresentable",
                ));
            }
            for column in 0..3 {
                matrix[3 * row + column] = divided(matrix[3 * row + column], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 3, &mut *cancelled).map_err(solve_error)?;
        for incident in 0..3 {
            if incoming[incident].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let rhs: [C64; 3] = core::array::from_fn(|r| {
                -divided(contact_entry(incoming[incident], incident, r), scales[r])
            });
            let mut solution = rhs;
            factor
                .solve_with_cancel(&mut solution, &mut *cancelled)
                .map_err(solve_error)?;
            if solution.iter().any(|&v| !finite(v)) {
                return Err(invalid("fluid-solid contact solution is unrepresentable"));
            }
            for row in 0..3 {
                let terms: [C64; 3] = core::array::from_fn(|c| matrix[3 * row + c] * solution[c]);
                let residual = (complex_sum(terms) - rhs[row]).abs();
                let magnitude = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                if !residual.is_finite() || (magnitude > 0.0 && residual / magnitude > 1e-11) {
                    return Err(invalid("fluid-solid contact continuity did not converge"));
                }
            }
            for row in 0..3 {
                result.scattering[row][incident] = solution[row];
            }
        }
        // SH has equal velocity polarization and opposite traction on its two
        // branches; a free shear boundary therefore reflects it with +1.
        if incoming[3].normal_power_factor_pa_s_m > 0.0 {
            result.scattering[3][3] = C64::ONE;
        }
        poll(cancelled)?;
        Ok(result)
    }
}

impl<'a, F: FluidAcoustics + ?Sized> PreparedFluidElasticInterface<'a, F> {
    pub fn fluid_source(&self) -> &'a F {
        self.fluid_source
    }
    pub const fn fluid_state(&self) -> FluidAcousticState {
        self.fluid_state
    }
    pub const fn solid_state(&self) -> &'a ElasticAcousticState {
        self.solid
    }
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential_slowness_s_m
    }
    pub fn incoming_mode(&self, channel: FluidElasticChannel) -> Option<ElasticWaveTrace> {
        let wave = self.incoming[channel.index()];
        (wave.normal_power_factor_pa_s_m > 0.0).then_some(wave)
    }
    pub fn outgoing_mode(&self, channel: FluidElasticChannel) -> ElasticWaveTrace {
        self.outgoing[channel.index()]
    }
    pub fn velocity_coefficients(&self, channel: FluidElasticChannel) -> Option<[C64; 4]> {
        self.incoming_mode(channel)?;
        Some(core::array::from_fn(|i| {
            self.scattering[i][channel.index()]
        }))
    }

    /// Apply peak modal velocities in channel order fluid P, solid P/SV/SH.
    /// All inputs share one frequency, tangential wave vector and boundary
    /// phase reference. No phase averaging or power renormalization is used.
    pub fn drive(
        &self,
        incoming: [C64; 4],
        limits: FluidDriveLimits,
    ) -> Result<FluidElasticResponse, MaterialError> {
        for (i, &a) in incoming.iter().enumerate() {
            if !finite(a) || (a != C64::ZERO && self.incoming[i].normal_power_factor_pa_s_m == 0.0)
            {
                return Err(invalid(
                    "fluid-solid drive requires finite propagating incident amplitudes",
                ));
            }
        }
        let outgoing: [C64; 4] = core::array::from_fn(|i| {
            complex_sum((0..4).map(|j| self.scattering[i][j] * incoming[j]))
        });
        if outgoing.iter().any(|&a| !finite(a)) {
            return Err(invalid(
                "fluid-solid outgoing amplitudes are unrepresentable",
            ));
        }
        let mut input_power = [0.0; 4];
        let mut output_power = input_power;
        for i in 0..4 {
            input_power[i] = power(self.incoming[i], incoming[i])?;
            output_power[i] = power(self.outgoing[i], outgoing[i])?;
        }
        let total_input = real_sum(input_power);
        let total_output = real_sum(output_power);
        let amplitude_sum = real_sum([incoming[0].abs(), outgoing[0].abs()]);
        let density_bound = amplitude_sum / self.fluid_state.sound_speed_m_s();
        let mach_bound = real_sum(
            [
                (self.incoming[0], incoming[0]),
                (self.outgoing[0], outgoing[0]),
            ]
            .map(|(w, a)| {
                let norm = w.velocity.iter().fold(0.0, |n, v| det::hypot(n, v.abs()));
                norm * (a.abs() / self.fluid_state.sound_speed_m_s())
            }),
        );
        if !density_bound.is_finite()
            || !mach_bound.is_finite()
            || density_bound > limits.maximum_relative_density
            || mach_bound > limits.maximum_particle_mach
        {
            return Err(invalid(
                "fluid-solid drive exceeds the declared fluid small-signal limits",
            ));
        }
        let strain_bound = real_sum(
            (1..4)
                .flat_map(|i| {
                    [
                        (self.incoming[i], incoming[i]),
                        (self.outgoing[i], outgoing[i]),
                    ]
                })
                .map(|(w, a)| {
                    w.strain_mandel_s_m
                        .iter()
                        .fold(0.0, |n, e| det::hypot(n, e.abs()))
                        * a.abs()
                }),
        );
        if !strain_bound.is_finite() || strain_bound > self.solid.maximum_absolute_mandel_strain() {
            return Err(invalid(
                "fluid-solid drive exceeds original solid strain domain",
            ));
        }
        let field = |solid: bool, traction: bool, component: usize| {
            complex_sum((0..4).filter(|&i| (i > 0) == solid).flat_map(|i| {
                let value = |w: ElasticWaveTrace| {
                    if traction {
                        w.traction_pa_s_m[component]
                    } else {
                        w.velocity[component]
                    }
                };
                [
                    value(self.incoming[i]) * incoming[i],
                    value(self.outgoing[i]) * outgoing[i],
                ]
            }))
        };
        let velocity: [[C64; 3]; 2] =
            core::array::from_fn(|s| core::array::from_fn(|i| field(s == 1, false, i)));
        let traction: [C64; 3] = core::array::from_fn(|i| field(true, true, i));
        let pressure = -field(false, true, 2);
        let flux = [
            (pressure * velocity[0][2].conj()).re * 0.5,
            -0.5 * real_sum((0..3).map(|i| (traction[i] * velocity[1][i].conj()).re)),
        ];
        if velocity
            .iter()
            .flatten()
            .chain(&traction)
            .any(|&v| !finite(v))
            || !finite(pressure)
            || [total_input, total_output]
                .into_iter()
                .chain(flux)
                .any(|v| !v.is_finite())
            || !balanced(total_input, total_output, total_input)
            || !balanced(flux[0], input_power[0] - output_power[0], total_input)
            || !balanced(
                flux[1],
                real_sum(output_power[1..].iter().copied())
                    - real_sum(input_power[1..].iter().copied()),
                total_input,
            )
        {
            return Err(invalid(
                "fluid-solid interface power is unrepresentable or does not close",
            ));
        }
        Ok(FluidElasticResponse {
            fluid_state: self.fluid_state,
            solid_state_identity: self.solid.material_state_identity(),
            incoming_velocity_m_s: incoming,
            outgoing_velocity_m_s: outgoing,
            incoming_normal_power_w_m2: input_power,
            outgoing_normal_power_w_m2: output_power,
            total_incoming_power_w_m2: total_input,
            total_outgoing_power_w_m2: total_output,
            boundary_velocity_m_s: velocity,
            fluid_pressure_pa: pressure,
            solid_traction_pa: traction,
            boundary_normal_flux_w_m2: flux,
            fluid_relative_density_bound: density_bound,
            fluid_particle_mach_bound: mach_bound,
            solid_strain_norm_bound: strain_bound,
        })
    }
}

fn contact_entry(wave: ElasticWaveTrace, channel: usize, row: usize) -> C64 {
    match row {
        0 => wave.velocity[2].scale(if channel == 0 { 1.0 } else { -1.0 }),
        1 => wave.traction_pa_s_m[2].scale(if channel == 0 { 1.0 } else { -1.0 }),
        _ => {
            if channel == 0 {
                C64::ZERO
            } else {
                wave.traction_pa_s_m[0]
            }
        }
    }
}

fn fluid_trace(
    state: FluidAcousticState,
    direction: f64,
    sine: f64,
    cosine: C64,
    tangential: f64,
    normal: C64,
) -> Result<ElasticWaveTrace, MaterialError> {
    let normal = normal.scale(direction);
    let velocity = [C64::from_re(sine), C64::ZERO, cosine.scale(direction)];
    let strain = [
        velocity[0].scale(-tangential),
        C64::ZERO,
        -(normal * velocity[2]),
        C64::ZERO,
        C64::ZERO,
        (normal * velocity[0] + velocity[2].scale(tangential))
            .scale(-core::f64::consts::FRAC_1_SQRT_2),
    ];
    let power = 0.5 * state.characteristic_impedance_pa_s_m() * cosine.re;
    if velocity.iter().chain(&strain).any(|&v| !finite(v))
        || !power.is_finite()
        || power < 0.0
        || (cosine.re > 0.0 && power == 0.0)
    {
        return Err(invalid("fluid wave trace is unrepresentable"));
    }
    Ok(ElasticWaveTrace {
        normal_slowness_s_m: normal,
        velocity,
        traction_pa_s_m: [
            C64::ZERO,
            C64::ZERO,
            C64::from_re(-state.characteristic_impedance_pa_s_m()),
        ],
        strain_mandel_s_m: strain,
        normal_power_factor_pa_s_m: power,
    })
}

#[cfg(test)]
mod tests;
