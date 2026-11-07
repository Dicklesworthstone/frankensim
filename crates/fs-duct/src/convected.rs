//! Exact plane sound in a rigid cylinder with uniform subsonic mean flow.
//!
//! The original [`FluidAcoustics`] owner supplies one coherent, frozen density
//! and adiabatic sound speed. The signed mean velocity is independently
//! declared along +x. In exp(-i omega t), the two characteristics have
//! `k+ = omega/(c+U0)` and `k- = -omega/(c-U0)`, with `u' = (p+ - p-)/(rho*c)`.
//! Forward and backward always mean +x and -x, even when the mean flow reverses.
//!
//! Physical disturbance power is **not** stationary pressure/volume-flow work.
//! With `rho' = p'/c²`, the perturbation mass flow and stagnation enthalpy are
//! `m' = A*(rho*u' + U0*rho')` and `H' = p'/rho + U0*u'`. Their cycle-mean
//! product `Re(H' conj(m'))/2` is Myers' disturbance power. A plane pressure
//! wave therefore carries `+A*(1+M)²*|p+|²/(2*rho*c)` or
//! `-A*(1-M)²*|p-|²/(2*rho*c)`. These are peak phasors.
//! The simpler linearized-Euler norm flux `Re(p' conj(u'))/2 + U0*E0`
//! is another conserved quantity for this uniform model; it is not the same
//! disturbance power and is not substituted for it.
//!
//! Derivation: S. W. Rienstra, *Fundamentals of Duct Acoustics*, sections 5.1
//! and A.3, <https://sjoerdr.win.tue.nl/papers/VKI_Rienstra.pdf>; the physical
//! mass-flux/enthalpy pairing is also equation (3) of Rienstra (2022),
//! <https://doi.org/10.1177/1475472X221107363>. Those references use +i omega t;
//! the fields here use the crate's opposite time convention.
//!
//! Admission is below the first transverse hard-wall cut-on,
//! `omega < j'11*c*sqrt(1-M²)/radius`. Sources and observation sections must
//! also be far enough from discontinuities for evanescent transverse modes to
//! be negligible. Only a straight cylinder is admitted. There is no shear,
//! swirl, nonuniform mean state, entropy mode, wall/bulk absorption, compliant
//! wall, end radiation or vortex shedding. The caller owns the small-signal
//! and uniform-flow assumptions and the original material's source limits.
//! No stationary [`crate::compose::TwoPort`] is exported for these moving-
//! medium fields. In particular, a sealed end cannot sustain the prescribed
//! nonzero through-flow. Prescribed incident waves remain available, while
//! [`termination`] supplies conjugate enthalpy/mass-flow end loads, coherent
//! sources and their exact multiple-reflection solution in the maintained flow.

use crate::fluid::FluidDuctError;
use crate::{DuctError, Segment};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, det};

/// Finite-length quasi-one-dimensional sound in a maintained isentropic flow.
pub mod profile;
/// Conjugate end loads, coherent driving and complete reflected fields.
pub mod termination;

/// First positive zero of J1', the first non-plane circular hard-wall mode.
const FIRST_TRANSVERSE_ROOT: f64 = 1.841_183_781_340_659_3;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn phase(angle: f64) -> C64 {
    let (sin, cos) = det::sin_cos(angle);
    C64::new(cos, sin)
}

/// One source-bound, frequency-specific uniform-flow cylinder.
/// Preparation and every subsequent field query use fixed work and no heap
/// allocation. The source is evaluated once; its identity and receipts stay
/// on the retained original owner, not an invented ideal-gas replacement.
#[derive(Debug)]
pub struct UniformFlowCylinder<'source, M: FluidAcoustics + ?Sized> {
    source: &'source M,
    state: FluidAcousticState,
    radius_m: f64,
    length_m: f64,
    area_m2: f64,
    omega: f64,
    mean_velocity_m_s: f64,
    mach_factors: [f64; 2],
    wavenumbers: [f64; 2],
    travel_times: [f64; 2],
    transmission: [C64; 2],
    transverse_cutoff: f64,
    pressure_power_weight: f64,
    disturbance_power_weights: [f64; 2],
    disturbance_energy_weights: [f64; 2],
    mass_flow_scale: f64,
}

/// Boundary result for a coherent pair of incoming pressure waves.
/// Left and right are the x=0 and x=L sections of the same uniform flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformFlowPower {
    /// Outgoing pressure waves [Pa], ordered [left-going at 0, right-going at L].
    pub outgoing_pressure_pa: [C64; 2],
    /// Sum of the two incoming Myers disturbance powers [W], both positive.
    pub incident_disturbance_power_w: f64,
    /// Sum of the two outgoing Myers disturbance powers [W], both positive.
    pub outgoing_disturbance_power_w: f64,
    /// Incoming minus outgoing power [W]. This lossless model has zero exact
    /// balance; the signed floating-point residual is not an absorbed-heat source.
    pub conservation_residual_w: f64,
    /// Signed disturbance power along +x [W], constant through the cylinder.
    pub axial_disturbance_power_w: f64,
}

/// Coherent interior observables from one pair of traveling acoustic waves.
/// All complex quantities are peak perturbation phasors, excluding the mean.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformFlowField {
    /// Acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Acoustic particle velocity along +x [m/s].
    pub particle_velocity_m_s: C64,
    /// Acoustic volume flow along +x [m3/s].
    pub volume_flow_m3_s: C64,
    /// Total mass-flow perturbation A*(rho*u' + U0*p'/c²) [kg/s].
    pub mass_flow_perturbation_kg_s: C64,
    /// Stagnation-enthalpy perturbation p'/rho + U0*u' [J/kg].
    pub stagnation_enthalpy_perturbation_j_kg: C64,
    /// Cycle-mean Re(p' conj(volume flow))/2 [W]. This omits mean-flow
    /// contributions and is deliberately distinct from disturbance power.
    pub pressure_work_power_w: f64,
    /// Cycle-mean Myers energy density [J/m3], including
    /// U0*Re(p' conj(u'))/(2*c²). Positive for admitted subsonic flow.
    pub disturbance_energy_density_j_m3: f64,
    /// Cycle-mean Re(H' conj(m'))/2 [W], signed along +x.
    pub disturbance_power_w: f64,
}

/// Prepared two-sided excitation. The original material and mean-flow state
/// remain bound to every observation; fields need no re-query or spatial solve.
#[derive(Debug)]
pub struct UniformFlowExcitation<'duct, 'source, M: FluidAcoustics + ?Sized> {
    duct: &'duct UniformFlowCylinder<'source, M>,
    incoming_pressure_pa: [C64; 2],
    power: UniformFlowPower,
    pressure_work_power_w: f64,
    disturbance_energy_density_j_m3: f64,
}

impl<'source, M: FluidAcoustics + ?Sized> UniformFlowCylinder<'source, M> {
    /// Prepare a rigid cylinder from one original acoustic material state.
    /// `mean_velocity_m_s` is signed along +x, with |U0| strictly below c.
    /// `omega` is positive angular frequency [rad/s]; the first transverse
    /// cut-on is a strict upper bound. Other segment kinds refuse.
    pub fn try_new(
        source: &'source M,
        segment: &Segment,
        omega: f64,
        mean_velocity_m_s: f64,
    ) -> Result<Self, FluidDuctError> {
        segment.validate()?;
        let Segment::Cylinder { radius, length } = *segment else {
            return Err(bad("uniform-flow plane acoustics requires a straight cylinder").into());
        };
        if !positive(omega) || !mean_velocity_m_s.is_finite() {
            return Err(
                bad("uniform-flow frequency must be positive finite and velocity finite").into(),
            );
        }
        let state = source.acoustic_state()?;
        let c = state.sound_speed_m_s();
        if mean_velocity_m_s.abs() >= c {
            return Err(
                bad("uniform-flow plane acoustics requires strictly subsonic mean flow").into(),
            );
        }
        let speeds = [c + mean_velocity_m_s, c - mean_velocity_m_s];
        let mach_factors = speeds.map(|speed| speed / c);
        let wavenumbers = [omega / speeds[0], -omega / speeds[1]];
        let travel_times = speeds.map(|speed| length / speed);
        let angles = [wavenumbers[0] * length, -wavenumbers[1] * length];
        let transverse_cutoff =
            FIRST_TRANSVERSE_ROOT * (c / radius) * det::sqrt(mach_factors[0] * mach_factors[1]);
        let area_m2 = core::f64::consts::PI * radius * radius;
        let pressure_power_weight = 0.5 * area_m2 / state.characteristic_impedance_pa_s_m();
        let disturbance_power_weights =
            mach_factors.map(|factor| pressure_power_weight * factor * factor);
        let disturbance_energy_weights =
            mach_factors.map(|factor| 0.5 * factor / state.adiabatic_bulk_modulus_pa());
        let mass_flow_scale = area_m2 / c;
        if !speeds
            .into_iter()
            .chain(mach_factors)
            .chain(travel_times)
            .chain(angles)
            .chain(disturbance_power_weights)
            .chain(disturbance_energy_weights)
            .chain([
                area_m2,
                transverse_cutoff,
                pressure_power_weight,
                mass_flow_scale,
                wavenumbers[0],
                -wavenumbers[1],
            ])
            .all(positive)
        {
            return Err(bad(
                "uniform-flow geometry, propagation or power weights are unrepresentable",
            )
            .into());
        }
        if omega >= transverse_cutoff {
            return Err(
                bad("uniform-flow frequency reaches the first transverse duct cut-on").into(),
            );
        }
        let transmission = angles.map(phase);
        if !transmission.into_iter().all(finite) {
            return Err(bad("uniform-flow propagation phase is unrepresentable").into());
        }
        Ok(Self {
            source,
            state,
            radius_m: radius,
            length_m: length,
            area_m2,
            omega,
            mean_velocity_m_s,
            mach_factors,
            wavenumbers,
            travel_times,
            transmission,
            transverse_cutoff,
            pressure_power_weight,
            disturbance_power_weights,
            disturbance_energy_weights,
            mass_flow_scale,
        })
    }

    /// The retained original material/EOS owner, including its source receipts.
    pub const fn source(&self) -> &'source M {
        self.source
    }
    /// The immutable coherent state used throughout this prepared solve.
    pub const fn acoustic_state(&self) -> &FluidAcousticState {
        &self.state
    }
    /// Cylinder radius [m].
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Cylinder length [m].
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Full cross-sectional area [m2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Angular frequency [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Authored uniform mean velocity along +x [m/s].
    pub const fn mean_velocity_m_s(&self) -> f64 {
        self.mean_velocity_m_s
    }
    /// Signed mean-flow Mach number.
    pub fn mach_number(&self) -> f64 {
        self.mean_velocity_m_s / self.state.sound_speed_m_s()
    }
    /// Signed axial wavenumbers [rad/m], ordered [+x wave, -x wave].
    pub const fn wavenumbers_rad_m(&self) -> [f64; 2] {
        self.wavenumbers
    }
    /// Positive travel times [s], ordered [left to right, right to left].
    pub const fn travel_times_s(&self) -> [f64; 2] {
        self.travel_times
    }
    /// Pressure transmission phasors, ordered [left to right, right to left].
    /// These are directional characteristics, not stationary p/U power waves.
    pub const fn pressure_transmission(&self) -> [C64; 2] {
        self.transmission
    }
    /// First non-plane hard-wall cut-on [rad/s], including mean-flow shift.
    pub const fn first_transverse_cutoff_rad_s(&self) -> f64 {
        self.transverse_cutoff
    }

    /// Resolve coherent incident pressures [right-going at x=0, left-going at
    /// x=L] into outgoing waves and reusable interior fields. Amplitudes must
    /// be finite and small enough for the caller's linear-acoustic regime.
    pub fn prepare_excitation(
        &self,
        incoming_pressure_pa: [C64; 2],
    ) -> Result<UniformFlowExcitation<'_, 'source, M>, DuctError> {
        if !incoming_pressure_pa.into_iter().all(finite) {
            return Err(bad("uniform-flow incident pressure must be finite"));
        }
        let squared = incoming_pressure_pa.map(C64::norm_sq);
        let outgoing_pressure_pa = [
            incoming_pressure_pa[1] * self.transmission[1],
            incoming_pressure_pa[0] * self.transmission[0],
        ];
        let forward_power = self.disturbance_power_weights[0] * squared[0];
        let backward_power = self.disturbance_power_weights[1] * squared[1];
        let incident_disturbance_power_w = forward_power + backward_power;
        let outgoing_disturbance_power_w = self.disturbance_power_weights[1]
            * outgoing_pressure_pa[0].norm_sq()
            + self.disturbance_power_weights[0] * outgoing_pressure_pa[1].norm_sq();
        let power = UniformFlowPower {
            outgoing_pressure_pa,
            incident_disturbance_power_w,
            outgoing_disturbance_power_w,
            conservation_residual_w: incident_disturbance_power_w - outgoing_disturbance_power_w,
            axial_disturbance_power_w: forward_power - backward_power,
        };
        let pressure_work_power_w = self.pressure_power_weight * (squared[0] - squared[1]);
        let disturbance_energy_density_j_m3 = self.disturbance_energy_weights[0] * squared[0]
            + self.disturbance_energy_weights[1] * squared[1];
        if !outgoing_pressure_pa.into_iter().all(finite)
            || ![
                incident_disturbance_power_w,
                outgoing_disturbance_power_w,
                power.conservation_residual_w,
                power.axial_disturbance_power_w,
                pressure_work_power_w,
                disturbance_energy_density_j_m3,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(bad(
                "uniform-flow incident fields or disturbance powers are unrepresentable",
            ));
        }
        Ok(UniformFlowExcitation {
            duct: self,
            incoming_pressure_pa,
            power,
            pressure_work_power_w,
            disturbance_energy_density_j_m3,
        })
    }
}

impl<'duct, 'source, M: FluidAcoustics + ?Sized> UniformFlowExcitation<'duct, 'source, M> {
    /// The original prepared cylinder and retained material owner.
    pub const fn duct(&self) -> &'duct UniformFlowCylinder<'source, M> {
        self.duct
    }
    /// Exact characteristic boundary solution and its floating-point balance.
    pub const fn boundary_power(&self) -> &UniformFlowPower {
        &self.power
    }

    /// Observe at x [m] measured from the left section, with 0 <= x <= L.
    /// The +x wave is anchored at 0 and the -x wave at L. Conservation-law
    /// quantities are evaluated once from their characteristic weights, so
    /// interference nodes and observation order cannot invent energy changes.
    pub fn at_position(&self, x_m: f64) -> Result<UniformFlowField, DuctError> {
        let duct = self.duct;
        if !x_m.is_finite() || x_m < 0.0 || x_m > duct.length_m {
            return Err(bad("uniform-flow observation must lie inside the cylinder"));
        }
        let forward = self.incoming_pressure_pa[0]
            * if x_m == duct.length_m {
                duct.transmission[0]
            } else {
                phase(duct.wavenumbers[0] * x_m)
            };
        let backward = self.incoming_pressure_pa[1]
            * if x_m == 0.0 {
                duct.transmission[1]
            } else {
                phase(-duct.wavenumbers[1] * (duct.length_m - x_m))
            };
        let pressure_pa = forward + backward;
        let particle_velocity_m_s =
            (forward - backward) / C64::from_re(duct.state.characteristic_impedance_pa_s_m());
        let volume_flow_m3_s = particle_velocity_m_s.scale(duct.area_m2);
        let enthalpy_forward = forward.scale(duct.mach_factors[0]);
        let enthalpy_backward = backward.scale(duct.mach_factors[1]);
        let mass_flow_perturbation_kg_s =
            (enthalpy_forward - enthalpy_backward).scale(duct.mass_flow_scale);
        let stagnation_enthalpy_perturbation_j_kg =
            (enthalpy_forward + enthalpy_backward) / C64::from_re(duct.state.density_kg_m3());
        if ![
            pressure_pa,
            particle_velocity_m_s,
            volume_flow_m3_s,
            mass_flow_perturbation_kg_s,
            stagnation_enthalpy_perturbation_j_kg,
        ]
        .into_iter()
        .all(finite)
        {
            return Err(bad(
                "uniform-flow interior pressure, flow or enthalpy is unrepresentable",
            ));
        }
        Ok(UniformFlowField {
            pressure_pa,
            particle_velocity_m_s,
            volume_flow_m3_s,
            mass_flow_perturbation_kg_s,
            stagnation_enthalpy_perturbation_j_kg,
            pressure_work_power_w: self.pressure_work_power_w,
            disturbance_energy_density_j_m3: self.disturbance_energy_density_j_m3,
            disturbance_power_w: self.power.axial_disturbance_power_w,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compose::TwoPort;
    use fs_material::gas::{GasSpec, GasState};

    fn gas() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
    }

    fn cylinder() -> Segment {
        Segment::Cylinder {
            radius: 0.012,
            length: 0.83,
        }
    }

    fn omega() -> f64 {
        core::f64::consts::TAU * 723.0
    }

    fn close(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 2e-13 * scale.max(1e-30),
            "actual {actual:e}, expected {expected:e}, scale {scale:e}"
        );
    }

    fn complex_close(actual: C64, expected: C64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 2e-13 * scale.max(1e-30),
            "actual {actual:?}, expected {expected:?}, scale {scale:e}"
        );
    }

    #[test]
    fn zero_flow_recovers_the_original_cylinder_basis_and_power() {
        let gas = gas();
        let segment = cylinder();
        let duct = UniformFlowCylinder::try_new(&gas, &segment, omega(), 0.0).unwrap();
        let wave = crate::fluid::segment_wave(&gas, duct.radius_m(), omega()).unwrap();
        let stationary =
            TwoPort::from_wave(&segment, &wave, wave.characteristic_impedance.re).unwrap();
        let incoming = [C64::new(1.2, -0.4), C64::new(-0.3, 0.7)];
        let excitation = duct.prepare_excitation(incoming).unwrap();
        let boundary = excitation.boundary_power();
        for i in 0..2 {
            let expected = stationary.scattering[2 * i] * incoming[0]
                + stationary.scattering[2 * i + 1] * incoming[1];
            complex_close(boundary.outgoing_pressure_pa[i], expected, expected.abs());
        }
        let matrix = crate::segment_matrix(&segment, &wave).unwrap();
        let left = excitation.at_position(0.0).unwrap();
        let right = excitation.at_position(duct.length_m()).unwrap();
        complex_close(
            left.pressure_pa,
            matrix[0] * right.pressure_pa + matrix[1] * right.volume_flow_m3_s,
            incoming[0].abs() + incoming[1].abs(),
        );
        complex_close(
            left.volume_flow_m3_s,
            matrix[2] * right.pressure_pa + matrix[3] * right.volume_flow_m3_s,
            (incoming[0].abs() + incoming[1].abs()) / wave.characteristic_impedance.abs(),
        );
        close(
            left.disturbance_power_w,
            left.pressure_work_power_w,
            boundary.incident_disturbance_power_w,
        );
        close(
            boundary.conservation_residual_w,
            0.0,
            boundary.incident_disturbance_power_w,
        );
    }

    #[test]
    fn counter_propagating_waves_carry_the_correct_myers_power_and_delays() {
        let gas = gas();
        let area = core::f64::consts::PI * 0.012 * 0.012;
        let amplitude = C64::new(2.0, -0.5);
        let stationary_power = 0.5 * area * amplitude.norm_sq() / (gas.density * gas.sound_speed);
        for mach in [-0.6, 0.0, 0.6] {
            let duct =
                UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), mach * gas.sound_speed)
                    .unwrap();
            let delays = duct.travel_times_s();
            close(
                delays[0] * (gas.sound_speed + duct.mean_velocity_m_s()),
                0.83,
                0.83,
            );
            close(
                delays[1] * (gas.sound_speed - duct.mean_velocity_m_s()),
                0.83,
                0.83,
            );
            for (slot, direction) in [(0, 1.0), (1, -1.0)] {
                let mut incoming = [C64::ZERO; 2];
                incoming[slot] = amplitude;
                let excitation = duct.prepare_excitation(incoming).unwrap();
                let factor = 1.0 + direction * mach;
                let expected_power = direction * stationary_power * factor * factor;
                let field = excitation.at_position(0.27).unwrap();
                close(field.disturbance_power_w, expected_power, stationary_power);
                // The wave's energy travels at its signed characteristic speed.
                let transport_speed = duct.mean_velocity_m_s() + direction * gas.sound_speed;
                close(
                    field.disturbance_energy_density_j_m3 * area * transport_speed,
                    expected_power,
                    stationary_power,
                );
                let trace_power = 0.5
                    * (field.stagnation_enthalpy_perturbation_j_kg
                        * field.mass_flow_perturbation_kg_s.conj())
                    .re;
                close(trace_power, expected_power, stationary_power);
            }
        }
    }

    #[test]
    fn coherent_fields_satisfy_linearized_euler_and_physical_power_pairing() {
        let gas = gas();
        let rho = gas.density;
        let c = gas.sound_speed;
        for mach in [-0.5, 0.0, 0.5] {
            let mean = mach * c;
            let duct = UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), mean).unwrap();
            let excitation = duct
                .prepare_excitation([C64::new(1.3, -0.2), C64::new(-0.4, 0.9)])
                .unwrap();
            for x in [0.1, 0.37, 0.72] {
                let field = excitation.at_position(x).unwrap();
                let h = 1e-5;
                let lo = excitation.at_position(x - h).unwrap();
                let hi = excitation.at_position(x + h).unwrap();
                let dp = (hi.pressure_pa - lo.pressure_pa).scale(0.5 / h);
                let du = (hi.particle_velocity_m_s - lo.particle_velocity_m_s).scale(0.5 / h);
                let dt_p = C64::new(0.0, -omega()) * field.pressure_pa;
                let dt_u = C64::new(0.0, -omega()) * field.particle_velocity_m_s;
                let pressure_residual = dt_p + dp.scale(mean) + du.scale(rho * c * c);
                let momentum_residual = (dt_u + du.scale(mean)).scale(rho) + dp;
                let pressure_scale = dt_p.abs() + mean.abs() * dp.abs() + rho * c * c * du.abs();
                let momentum_scale = rho * (dt_u.abs() + mean.abs() * du.abs()) + dp.abs();
                assert!(pressure_residual.abs() < 1e-7 * pressure_scale);
                assert!(momentum_residual.abs() < 1e-7 * momentum_scale);

                let density_perturbation = field.pressure_pa.scale(1.0 / (c * c));
                let mass_flow = field.volume_flow_m3_s.scale(rho)
                    + density_perturbation.scale(mean * duct.area_m2());
                let enthalpy =
                    field.pressure_pa.scale(1.0 / rho) + field.particle_velocity_m_s.scale(mean);
                complex_close(
                    field.mass_flow_perturbation_kg_s,
                    mass_flow,
                    mass_flow.abs(),
                );
                complex_close(
                    field.stagnation_enthalpy_perturbation_j_kg,
                    enthalpy,
                    enthalpy.abs(),
                );
                let physical_power = 0.5 * (enthalpy * mass_flow.conj()).re;
                let pressure_work = 0.5 * (field.pressure_pa * field.volume_flow_m3_s.conj()).re;
                let energy = 0.25
                    * (field.pressure_pa.norm_sq() / (rho * c * c)
                        + rho * field.particle_velocity_m_s.norm_sq())
                    + 0.5 * mean / (c * c)
                        * (field.pressure_pa * field.particle_velocity_m_s.conj()).re;
                let scale = excitation.boundary_power().incident_disturbance_power_w;
                close(field.disturbance_power_w, physical_power, scale);
                close(field.pressure_work_power_w, pressure_work, scale);
                close(field.disturbance_energy_density_j_m3, energy, energy);
            }
        }
    }

    #[test]
    fn equal_pressure_waves_have_mean_flow_power_even_when_pressure_work_cancels() {
        let gas = gas();
        let mach = 0.4;
        let duct = UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), mach * gas.sound_speed)
            .unwrap();
        let excitation = duct.prepare_excitation([C64::ONE, C64::ONE]).unwrap();
        let field = excitation.at_position(0.2).unwrap();
        let expected = 2.0 * duct.area_m2() * mach / (gas.density * gas.sound_speed);
        assert_eq!(field.pressure_work_power_w, 0.0);
        close(field.disturbance_power_w, expected, expected);
        assert!(field.disturbance_power_w > 0.0);
    }

    #[test]
    fn reversing_flow_sources_and_coordinates_reverses_the_physical_flux() {
        let gas = gas();
        let forward =
            UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), 0.37 * gas.sound_speed)
                .unwrap();
        let reverse =
            UniformFlowCylinder::try_new(&gas, &cylinder(), omega(), -0.37 * gas.sound_speed)
                .unwrap();
        let incoming = [C64::new(0.8, -0.9), C64::new(1.2, 0.3)];
        let a = forward.prepare_excitation(incoming).unwrap();
        let b = reverse
            .prepare_excitation([incoming[1], incoming[0]])
            .unwrap();
        assert_eq!(
            forward.travel_times_s(),
            [reverse.travel_times_s()[1], reverse.travel_times_s()[0]]
        );
        for x in [0.0, 0.23, 0.71, 0.83] {
            let original = a.at_position(forward.length_m() - x).unwrap();
            let reversed = b.at_position(x).unwrap();
            complex_close(reversed.pressure_pa, original.pressure_pa, 3.0);
            complex_close(
                reversed.particle_velocity_m_s,
                -original.particle_velocity_m_s,
                3.0 / (gas.density * gas.sound_speed),
            );
            complex_close(
                reversed.mass_flow_perturbation_kg_s,
                -original.mass_flow_perturbation_kg_s,
                original.mass_flow_perturbation_kg_s.abs(),
            );
            complex_close(
                reversed.stagnation_enthalpy_perturbation_j_kg,
                original.stagnation_enthalpy_perturbation_j_kg,
                original.stagnation_enthalpy_perturbation_j_kg.abs(),
            );
            assert_eq!(reversed.disturbance_power_w, -original.disturbance_power_w);
            close(
                reversed.disturbance_energy_density_j_m3,
                original.disturbance_energy_density_j_m3,
                original.disturbance_energy_density_j_m3,
            );
        }
    }

    #[test]
    fn original_material_is_retained_and_queried_only_once() {
        use std::cell::Cell;
        struct Counted {
            state: FluidAcousticState,
            queries: Cell<usize>,
        }
        impl FluidAcoustics for Counted {
            fn acoustic_state(&self) -> Result<FluidAcousticState, fs_material::MaterialError> {
                self.queries.set(self.queries.get() + 1);
                Ok(self.state)
            }
        }
        let source = Counted {
            state: gas().acoustic_state().unwrap(),
            queries: Cell::new(0),
        };
        let duct = UniformFlowCylinder::try_new(&source, &cylinder(), omega(), 20.0).unwrap();
        assert!(core::ptr::eq(duct.source(), &source));
        let excitation = duct
            .prepare_excitation([C64::ONE, C64::new(0.2, -0.3)])
            .unwrap();
        for i in 0..20 {
            excitation
                .at_position(duct.length_m() * f64::from(i) / 20.0)
                .unwrap();
        }
        assert_eq!(source.queries.get(), 1);
        let x = 0.23;
        let original = excitation.at_position(x).unwrap();
        excitation.at_position(0.67).unwrap();
        assert_eq!(excitation.at_position(x).unwrap(), original);
        let mut stale = gas();
        stale.sound_speed *= 1.1;
        assert!(matches!(
            UniformFlowCylinder::try_new(&stale, &cylinder(), omega(), 0.0),
            Err(FluidDuctError::Material(_))
        ));
    }

    #[test]
    fn sonic_flow_non_cylinders_and_transverse_cut_on_refuse() {
        let gas = gas();
        let segment = cylinder();
        for velocity in [
            gas.sound_speed,
            -gas.sound_speed,
            1.01 * gas.sound_speed,
            f64::NAN,
            f64::INFINITY,
        ] {
            assert!(UniformFlowCylinder::try_new(&gas, &segment, omega(), velocity).is_err());
        }
        for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(UniformFlowCylinder::try_new(&gas, &segment, frequency, 0.0).is_err());
        }
        let cone = Segment::Cone {
            inlet_radius: 0.01,
            outlet_radius: 0.02,
            length: 0.83,
        };
        assert!(UniformFlowCylinder::try_new(&gas, &cone, omega(), 0.0).is_err());
        let invalid = Segment::Cylinder {
            radius: 0.0,
            length: 0.83,
        };
        assert!(UniformFlowCylinder::try_new(&gas, &invalid, omega(), 0.0).is_err());
        let still = UniformFlowCylinder::try_new(&gas, &segment, omega(), 0.0).unwrap();
        let flowing =
            UniformFlowCylinder::try_new(&gas, &segment, omega(), 0.6 * gas.sound_speed).unwrap();
        let cutoff = still.first_transverse_cutoff_rad_s();
        close(
            flowing.first_transverse_cutoff_rad_s(),
            0.8 * cutoff,
            cutoff,
        );
        assert!(UniformFlowCylinder::try_new(&gas, &segment, 0.9 * cutoff, 0.0).is_ok());
        assert!(
            UniformFlowCylinder::try_new(&gas, &segment, 0.9 * cutoff, 0.6 * gas.sound_speed)
                .is_err()
        );
        assert!(UniformFlowCylinder::try_new(&gas, &segment, cutoff, 0.0).is_err());
        let excitation = flowing.prepare_excitation([C64::ONE, C64::ZERO]).unwrap();
        for x in [-0.01, 0.84, f64::NAN, f64::INFINITY] {
            assert!(excitation.at_position(x).is_err());
        }
        assert!(
            flowing
                .prepare_excitation([C64::new(f64::NAN, 0.0), C64::ZERO])
                .is_err()
        );
        assert!(
            flowing
                .prepare_excitation([C64::from_re(f64::MAX), C64::ZERO])
                .is_err()
        );
    }
}
