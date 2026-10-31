//! Original radial constitutive owner plus one jointly solved inertia matrix.

use super::*;
use core::f64::consts::PI;
use fs_la::LuWorkspace;

pub(super) struct Workspace {
    solver: LuWorkspace,
    matrix: Vec<f64>,
    rhs: [f64; MAX_INTERACTING_BUBBLES],
    solution: [f64; MAX_INTERACTING_BUBBLES],
    native: [Option<BubbleObservation>; MAX_INTERACTING_BUBBLES],
}
impl Workspace {
    pub(super) fn new(count: usize) -> Result<Self, BubbleClusterError> {
        Ok(Self {
            solver: LuWorkspace::new(count).map_err(|_| BubbleClusterError::LinearSolve)?,
            matrix: vec![0.0; count * count],
            rhs: [0.0; MAX_INTERACTING_BUBBLES],
            solution: [0.0; MAX_INTERACTING_BUBBLES],
            native: [None; MAX_INTERACTING_BUBBLES],
        })
    }
}

pub(super) struct Frame {
    pub accelerations: [f64; MAX_INTERACTING_BUBBLES],
    effective_pressures: [f64; MAX_INTERACTING_BUBBLES],
    pub viscous_dissipation_w: f64,
    pub drive_power_w: f64,
}

impl BubbleCluster {
    pub(super) fn evaluate(
        &self,
        time_s: f64,
        values: &[f64; DIMENSION],
        pressure_pa: f64,
        work: &mut Workspace,
    ) -> Result<Frame, BubbleClusterError> {
        let count = self.members.len();
        let mut volume_rate = 0.0;
        let mut dissipation = 0.0;
        let mut roots = [0.0; MAX_INTERACTING_BUBBLES];
        for (i, member) in self.members.iter().enumerate() {
            // Raw RK stages use the ORIGINAL model's domain checks so an
            // out-of-radius trial remains refinable by the shared stepper.
            let state = BubbleState {
                time_s,
                radius_m: values[2 * i],
                wall_velocity_m_s: values[2 * i + 1],
            };
            let native = member
                .bubble
                .observe(&state, pressure_pa)
                .map_err(|cause| super::member(i, cause))?;
            let radius = state.radius_m();
            let velocity = state.wall_velocity_m_s();
            roots[i] = radius.sqrt();
            work.rhs[i] = finite(
                roots[i] * (radius * native.wall_acceleration_m_s2 + 2.0 * velocity * velocity),
                "normalized cluster acceleration forcing",
            )?;
            volume_rate += native.volume_rate_m3_s;
            dissipation += native.viscous_dissipation_w;
            work.native[i] = Some(native);
        }
        let mut acceleration = [0.0; MAX_INTERACTING_BUBBLES];
        let mut sources = [0.0; MAX_INTERACTING_BUBBLES];
        if count == 1 {
            // Preserve the exact original single-bubble arithmetic and source.
            let native = work.native[0].expect("populated native member");
            acceleration[0] = native.wall_acceleration_m_s2;
            sources[0] = native.volume_acceleration_m3_s2 / (4.0 * PI);
        } else {
            for i in 0..count {
                for j in 0..count {
                    let value = if i == j {
                        1.0
                    } else {
                        roots[i] * roots[j] * self.inverse_distances[i * count + j]
                    };
                    if !positive(value) {
                        return Err(BubbleClusterError::OutsideDomain(
                            "cluster inertia is unrepresentable",
                        ));
                    }
                    work.matrix[i * count + j] = value;
                }
            }
            work.solver
                .solve_into(
                    &work.matrix,
                    &work.rhs[..count],
                    &mut work.solution[..count],
                )
                .map_err(|_| BubbleClusterError::LinearSolve)?;
            for i in 0..count {
                let radius = values[2 * i];
                let velocity = values[2 * i + 1];
                sources[i] = finite(roots[i] * work.solution[i], "coupled monopole source")?;
                acceleration[i] = finite(
                    (sources[i] - 2.0 * radius * velocity * velocity) / (radius * radius),
                    "coupled wall acceleration",
                )?;
            }
        }
        let sound_speed = self.material.liquid_sound_speed_m_s();
        for i in 0..count {
            let radius = values[2 * i];
            let velocity = values[2 * i + 1];
            let pressure_ratio = finite(
                ((radius * acceleration[i] + 1.5 * velocity * velocity) / sound_speed).abs()
                    / sound_speed,
                "coupled RP pressure/bulk-modulus ratio",
            )?;
            if pressure_ratio > self.limits.maximum_liquid_pressure_ratio {
                return Err(super::member(
                    i,
                    CavitationError::PressureRatioExceeded {
                        pressure_ratio,
                        maximum_pressure_ratio: self.limits.maximum_liquid_pressure_ratio,
                    },
                ));
            }
        }
        let density = self.material.liquid_density_kg_m3();
        let mut effective_pressures = [0.0; MAX_INTERACTING_BUBBLES];
        for i in 0..count {
            let mut pressure = pressure_pa;
            for (j, source) in sources[..count].iter().enumerate() {
                if i != j {
                    pressure += density * source * self.inverse_distances[i * count + j];
                }
            }
            effective_pressures[i] = finite(pressure, "effective cluster pressure")?;
        }
        Ok(Frame {
            accelerations: acceleration,
            effective_pressures,
            viscous_dissipation_w: finite(dissipation, "cluster viscous loss rate")?,
            drive_power_w: finite(
                -(pressure_pa - self.reference_pressure_pa) * volume_rate,
                "cluster pressure work rate",
            )?,
        })
    }

    pub(super) fn observe_values(
        &self,
        time_s: f64,
        values: &[f64; DIMENSION],
        pressure_pa: f64,
        work: &mut Workspace,
    ) -> Result<BubbleClusterObservation, BubbleClusterError> {
        let frame = self.evaluate(time_s, values, pressure_pa, work)?;
        let count = self.members.len();
        let mut bubbles = Vec::with_capacity(count);
        let mut volume = 0.0;
        let mut rate = 0.0;
        let mut acceleration = 0.0;
        let mut kinetic = 0.0;
        let mut potential = 0.0;
        for (i, member) in self.members.iter().enumerate() {
            let native = work.native[i].expect("populated native member");
            let current = super::super::observation(
                self.material,
                &native.state,
                frame.effective_pressures[i],
                native.gas_partial_pressure_pa,
                frame.accelerations[i],
            )
            .map_err(|cause| super::member(i, cause))?;
            let gas = member
                .bubble
                .gas()
                .volume_response(current.volume_m3)
                .map_err(|cause| super::member(i, cause))?;
            volume += current.volume_m3;
            rate += current.volume_rate_m3_s;
            acceleration += current.volume_acceleration_m3_s2;
            kinetic += current.liquid_kinetic_energy_j;
            potential += gas.mechanical_potential_j
                + current.surface_energy_j
                + (self.reference_pressure_pa - self.material.saturation_pressure_pa())
                    * current.volume_m3;
            bubbles.push(current);
        }
        for i in 0..count {
            for j in i + 1..count {
                // Algebraically rho*Q_i*Q_j/(4pi*d), evaluated through each
                // native self energy to avoid overflowing an intermediate Q_i*Q_j.
                let coupling = values[2 * i].sqrt()
                    * values[2 * j].sqrt()
                    * self.inverse_distances[i * count + j];
                kinetic += 2.0
                    * coupling
                    * bubbles[i].liquid_kinetic_energy_j.sqrt()
                    * bubbles[j].liquid_kinetic_energy_j.sqrt()
                    * bubbles[i].volume_rate_m3_s.signum()
                    * bubbles[j].volume_rate_m3_s.signum();
            }
        }
        if !kinetic.is_finite() || kinetic < 0.0 {
            return Err(BubbleClusterError::OutsideDomain(
                "cluster kinetic energy is unrepresentable",
            ));
        }
        Ok(BubbleClusterObservation {
            time_s,
            external_absolute_pressure_pa: pressure_pa,
            bubbles,
            total_volume_m3: finite(volume, "cluster volume")?,
            total_volume_rate_m3_s: finite(rate, "cluster volume rate")?,
            total_volume_acceleration_m3_s2: finite(acceleration, "cluster volume acceleration")?,
            liquid_kinetic_energy_j: kinetic,
            bubble_potential_j: finite(potential, "cluster mechanical potential")?,
            total_mechanical_energy_j: finite(kinetic + potential, "cluster mechanical storage")?,
            viscous_dissipation_w: frame.viscous_dissipation_w,
            drive_power_w: frame.drive_power_w,
        })
    }
}
