//! Sound-driven original mechanics and coherent incident/motion interference.
use super::*;
use fs_bem::helmholtz::{
    bc::{Admittance, IncidentField},
    molecular::MolecularScatteringOptions,
};

/// Additional illumination-stage resources. The retained source solution bank,
/// original BEM assembly and cached mechanical solve keep their existing limits.
#[derive(Debug, Clone, Copy)]
pub struct MolecularIncidentLimits {
    pub minimum_source_clearance_m: f64,
    /// All analytic source terms times (panels + receivers), over the batch.
    pub maximum_source_point_products: u64,
    /// Conservative live logical bytes, including blocked solutions, force
    /// columns, mechanical responses, total boundary/receiver fields and scratch:
    /// 16*[cases*(12*p+8*m+8*r+2*d+64)+16*p+4*r+4*m+2*t]. Here p/m/r are panel,
    /// mode and receiver counts; d=modes+memories, t=all incident source terms.
    /// The last term covers finite-source position/clearance scratch.
    /// Adaptive receivers replace cases*8*r+4*r with cases*12*r+6*r, retaining
    /// motion and total estimates plus shared paired-field scratch. Their
    /// separately admitted adaptive operator rows are outside this reservation.
    /// Allocator overhead and the
    /// separately admitted original operator/source bank are excluded.
    pub maximum_live_bytes: usize,
    /// Projection/reconstruction and boundary-observation reservation:
    /// cases*p*(4*m+4*r+24). The original owner admits BEM assembly separately.
    pub maximum_panel_products: u64,
}

/// Cycle-average work for peak phasors. `total_exterior_w` is SIGNED total
/// outward boundary work under illumination; it is not exterior gas heat or
/// outgoing radiation alone. Sound can transfer energy INTO the source.
/// Only the motion-only exterior term has the usual outgoing-field meaning.
///
/// The blocked, motion and interference channels are independently contracted.
/// Total channels are independently contracted from summed physical traces;
/// neither a cross term nor a small heat channel is obtained by subtraction.
#[derive(Debug, Clone, Copy)]
pub struct MolecularIncidentPower {
    pub external_mechanical_input_w: f64,
    pub blocked_acoustic_input_w: f64,
    pub kelvin_w: f64,
    pub relaxation_w: f64,
    pub blocked_exterior_w: f64,
    pub motion_exterior_w: f64,
    pub interference_exterior_w: f64,
    pub total_exterior_w: f64,
    pub blocked_wall_w: f64,
    pub motion_wall_w: f64,
    pub interference_wall_w: f64,
    pub total_wall_w: f64,
    /// External mechanical input minus structural heats minus TOTAL outward
    /// boundary work minus TOTAL wall absorption. Incident emitter work enters
    /// through that signed total boundary field, not as an added radiation loss.
    pub balance_w: f64,
    pub apparent_power_scale_w: f64,
    pub relative_balance: f64,
}

/// Actual sound-driven motion, total traction and matched total observations.
#[derive(Debug, Clone)]
pub struct MolecularIncidentResponse {
    /// Mechanics solved with external force plus blocked-body acoustic force.
    /// This member's receivers/reaction/power describe the motion-radiated field
    /// and its effective drive. Use the explicit total fields below for sound.
    pub motion: MolecularLoadedResponse,
    /// Force ON the original modal coordinates from the stationary illuminated
    /// body's TOTAL pressure, including its coherent scattered pressure.
    pub blocked_acoustic_force: Vec<C64>,
    /// Direct projection of total pressure traction ON the original source.
    /// It equals blocked_acoustic_force minus Z*modal_velocity, up to roundoff.
    pub total_fluid_reaction: Vec<C64>,
    /// Same declared order/locations as the prepared harmonic receiver bank.
    pub total_receivers: Vec<HarmonicField>,
    /// Absolute blocked-field plus coherent-motion quadrature estimates.
    /// They exclude source, solve and boundary discretization error and do not
    /// bound relative error near interference nulls. None denotes centroid use.
    pub total_receiver_quadrature_error_estimates: Option<Vec<MolecularReceiverError>>,
    pub total_surface_pressure_pa: Vec<C64>,
    pub total_surface_normal_velocity_m_s: Vec<C64>,
    pub power: MolecularIncidentPower,
}

/// Componentwise Neumaier sums keep quiet coherent terms between cancelling
/// loud terms. This is numerical accumulation, not a passivity correction.
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn finish(self) -> f64 {
        self.value + self.correction
    }
}

struct BoundaryField {
    pressure: Vec<C64>,
    velocity: Vec<C64>,
}

fn sum_receiver(
    pressure_pa: C64,
    particle_velocity_m_s: [C64; 3],
    radiated: &HarmonicField,
) -> Result<HarmonicField, MolecularModesError> {
    let total = HarmonicField {
        pressure_pa: pressure_pa + radiated.pressure_pa,
        particle_velocity_m_s: std::array::from_fn(|axis| {
            particle_velocity_m_s[axis] + radiated.particle_velocity_m_s[axis]
        }),
    };
    if !finite(total.pressure_pa) || total.particle_velocity_m_s.iter().any(|&u| !finite(u)) {
        return Err(bad("molecular illuminated total receiver overflow").into());
    }
    Ok(total)
}

impl PreparedMolecularLoadedModes<'_> {
    /// Illuminate the ORIGINAL source, solve its resulting loaded motion, and
    /// observe the coherent incident + blocked scattering + motion field.
    /// Every illumination has one complete external mechanical-force column;
    /// use an explicit zero vector for purely sound-driven motion.
    ///
    /// One blocked-body BEM assembly/LU is shared by the incident batch. The
    /// previously cached mechanical/memory LU and motion boundary/receiver bank
    /// are reused. The source bank's earlier BEM LU is not retained or claimed
    /// to be reused. All original material arms and mutual loads survive.
    /// If adaptive receivers were selected, one additional receiver row bank
    /// serves the complete blocked incident batch. Motion and total estimates
    /// retain the same original receiver choice without new motion quadrature.
    ///
    /// The same immutable molecular exterior, wall laws, formulation, original
    /// source modes and receiver declarations used at preparation are enforced.
    /// Finite emitters must remain outside the body. Original receiver clearance
    /// applies; neither field accuracy nor spatial convergence is inferred from
    /// clearance. Failed/cancelled batches publish no partial responses.
    pub fn harmonic_incident_response_batch(
        &self,
        incidents: &[IncidentField<'_>],
        mechanical_forces: &[&[C64]],
        limits: MolecularIncidentLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularIncidentResponse>, MolecularModesError> {
        poll(&mut cancelled)?;
        self.admit_force_batch(mechanical_forces)?;
        if incidents.len() != mechanical_forces.len() {
            return Err(bad(
                "molecular illumination needs one complete force column per incident field",
            )
            .into());
        }
        let (m, p, d) = (
            self.source.mode_count(),
            self.source.panel_count(),
            self.rows.len(),
        );
        let receivers = self.radiation.harmonic().receivers();
        let (c, r) = (incidents.len(), receivers.len());
        let source_terms: u128 = incidents
            .iter()
            .map(|incident| match incident {
                IncidentField::PlaneWave(_) => 1,
                IncidentField::Monopoles(sources) => sources.len() as u128,
                IncidentField::ForceDipoles(sources) => sources.len() as u128,
                IncidentField::Sources(sources) => sources.len() as u128,
            })
            .sum();
        let (case_receiver_entries, shared_receiver_entries) = if self.near_field_options.is_some()
        {
            (12, 6)
        } else {
            (8, 4)
        };
        let live_bytes = 16
            * (c as u128
                * (12 * p as u128
                    + 8 * m as u128
                    + case_receiver_entries * r as u128
                    + 2 * d as u128
                    + 64)
                + 16 * p as u128
                + shared_receiver_entries * r as u128
                + 4 * m as u128
                + 2 * source_terms);
        let panel_products = c as u128 * p as u128 * (4 * m as u128 + 4 * r as u128 + 24);
        if !limits.minimum_source_clearance_m.is_finite()
            || limits.minimum_source_clearance_m <= 0.0
            || live_bytes > limits.maximum_live_bytes as u128
            || panel_products > limits.maximum_panel_products as u128
            || limits.maximum_source_point_products == 0
            || source_terms * (p as u128 + r as u128) > limits.maximum_source_point_products as u128
        {
            return Err(bad("molecular illumination needs complete paired finite forces and admitted source, field and projection resources").into());
        }
        let omega = self.radiation.harmonic().angular_frequency_rad_s();
        let beta = self
            .boundary
            .as_ref()
            .map(|boundary| boundary.normalized_admittance(omega, &mut cancelled))
            .transpose()?;
        let admittance = beta
            .as_ref()
            .map_or(Admittance::Rigid, |beta| Admittance::Local {
                beta,
                declared_active: false,
            });
        let blocked = self.exterior.solve_scattering_batch(
            &self.source.surface,
            omega,
            incidents,
            admittance,
            self.source.settings.formulation,
            &MolecularScatteringOptions {
                boundary: self.source.settings.boundary,
                minimum_source_clearance_m: limits.minimum_source_clearance_m,
                maximum_source_point_products: limits.maximum_source_point_products,
            },
            &mut cancelled,
        )?;
        let points: Vec<_> = receivers.iter().map(|r| r.position_m).collect();
        // The entire incident batch belongs to this same immutable exterior.
        // Prepare once before consuming it; each case only contracts its traces.
        let blocked_observer = self
            .near_field_options
            .map(|options| {
                blocked[0].solution().scattered.prepare_near_field(
                    &self.source.surface,
                    &points,
                    self.minimum_clearance_m,
                    options,
                    &mut cancelled,
                )
            })
            .transpose()?;
        let mut acoustic_forces = Vec::with_capacity(c);
        let mut effective_forces = Vec::with_capacity(c);
        for (field, &force) in blocked.iter().zip(mechanical_forces) {
            poll(&mut cancelled)?;
            let solution = field.solution();
            let acoustic = self.project_pressure_force(
                |panel| {
                    solution.incident_pressure_pa[panel] + solution.scattered.pressure_pa()[panel]
                },
                &mut cancelled,
            )?;
            let effective: Vec<_> = force.iter().zip(&acoustic).map(|(&a, &b)| a + b).collect();
            if effective.iter().any(|&force| !finite(force)) {
                return Err(bad("molecular incident effective force overflow").into());
            }
            acoustic_forces.push(acoustic);
            effective_forces.push(effective);
        }
        let force_columns: Vec<_> = effective_forces.iter().map(Vec::as_slice).collect();
        let motion = self.harmonic_response_batch(&force_columns, &mut cancelled)?;
        let mut output = Vec::with_capacity(c);
        for (((field, acoustic_force), motion), &external_force) in blocked
            .into_iter()
            .zip(acoustic_forces)
            .zip(motion)
            .zip(mechanical_forces)
        {
            poll(&mut cancelled)?;
            let (total_receivers, total_receiver_quadrature_error_estimates) =
                if let Some(observer) = &blocked_observer {
                    let samples =
                        field.total_pressure_velocity_with_near_field(observer, &mut cancelled)?;
                    // Allocate retained totals only after the observer releases
                    // its temporary analytic-incident sample vector.
                    let mut totals = Vec::with_capacity(r);
                    let mut estimates = Vec::with_capacity(r);
                    let motion_estimates = motion
                        .receiver_quadrature_error_estimates
                        .as_ref()
                        .ok_or_else(|| {
                            bad("adaptive molecular motion has no receiver estimates")
                        })?;
                    for i in 0..r {
                        totals.push(sum_receiver(
                            samples.scalar.pressure[i],
                            samples.particle_velocity_m_s[i],
                            &motion.receivers[i],
                        )?);
                        let estimate = MolecularReceiverError {
                            pressure_pa: samples.scalar.quadrature_error_estimate_pa[i]
                                + motion_estimates[i].pressure_pa,
                            particle_velocity_m_s: std::array::from_fn(|axis| {
                                samples.quadrature_error_estimate_m_s[i][axis]
                                    + motion_estimates[i].particle_velocity_m_s[axis]
                            }),
                        };
                        if !estimate.pressure_pa.is_finite()
                            || estimate
                                .particle_velocity_m_s
                                .iter()
                                .any(|v| !v.is_finite())
                        {
                            return Err(bad(
                                "molecular total receiver quadrature estimate overflow",
                            )
                            .into());
                        }
                        estimates.push(estimate);
                    }
                    (totals, Some(estimates))
                } else {
                    let samples = field.total_pressure_velocity_at(
                        &points,
                        self.minimum_clearance_m,
                        &mut cancelled,
                    )?;
                    let mut totals = Vec::with_capacity(r);
                    for (blocked, radiated) in samples.into_iter().zip(&motion.receivers) {
                        totals.push(sum_receiver(
                            blocked.pressure_pa,
                            blocked.velocity_m_s,
                            radiated,
                        )?);
                    }
                    (totals, None)
                };
            let blocked = BoundaryField {
                pressure: field.solution().total_pressure_pa(),
                velocity: field.solution().total_normal_velocity_m_s(),
            };
            let mut radiated = BoundaryField {
                pressure: vec![C64::ZERO; p],
                velocity: vec![C64::ZERO; p],
            };
            self.radiation.evaluate_boundary_coherent(
                &motion.modal_velocity,
                &mut radiated.pressure,
                &mut radiated.velocity,
                &mut cancelled,
            )?;
            let mut total = BoundaryField {
                pressure: Vec::with_capacity(p),
                velocity: Vec::with_capacity(p),
            };
            for panel in 0..p {
                if panel % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                let pressure = blocked.pressure[panel] + radiated.pressure[panel];
                let velocity = blocked.velocity[panel] + radiated.velocity[panel];
                if !finite(pressure) || !finite(velocity) {
                    return Err(bad("molecular illuminated total boundary overflow").into());
                }
                total.pressure.push(pressure);
                total.velocity.push(velocity);
            }
            let reaction =
                self.project_pressure_force(|panel| total.pressure[panel], &mut cancelled)?;
            let power = self.incident_power(
                external_force,
                &acoustic_force,
                &motion,
                &blocked,
                &radiated,
                &total,
                beta.as_deref(),
                &mut cancelled,
            )?;
            output.push(MolecularIncidentResponse {
                motion,
                blocked_acoustic_force: acoustic_force,
                total_fluid_reaction: reaction,
                total_receivers,
                total_receiver_quadrature_error_estimates,
                total_surface_pressure_pa: total.pressure,
                total_surface_normal_velocity_m_s: total.velocity,
                power,
            });
        }
        poll(&mut cancelled)?;
        Ok(output)
    }

    fn project_pressure_force(
        &self,
        pressure: impl Fn(usize) -> C64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, MolecularModesError> {
        let mut forces = Vec::with_capacity(self.source.mode_count());
        for weights in &self.source.weights {
            let mut sums = [Sum::default(); 2];
            for (panel, (&weight, &area)) in
                weights.iter().zip(self.source.surface.areas()).enumerate()
            {
                if panel % 128 == 0 {
                    poll(cancelled)?;
                }
                let force = pressure(panel).scale(-weight * area);
                sums[0].add(force.re);
                sums[1].add(force.im);
            }
            let force = C64::new(sums[0].finish(), sums[1].finish());
            if !finite(force) {
                return Err(bad("molecular incident pressure-force projection overflow").into());
            }
            forces.push(force);
        }
        Ok(forces)
    }

    #[allow(clippy::too_many_arguments)]
    fn incident_power(
        &self,
        external: &[C64],
        acoustic: &[C64],
        motion: &MolecularLoadedResponse,
        blocked: &BoundaryField,
        radiated: &BoundaryField,
        total: &BoundaryField,
        beta: Option<&[C64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularIncidentPower, MolecularModesError> {
        // Normalize by all independent forces and fields, never their possibly
        // cancelling effective mechanical force. Restore physical watts once.
        let mut amplitude: f64 = 0.0;
        for values in [
            external,
            acoustic,
            &motion.modal_velocity,
            &motion.elastic_memory_sqrt_j,
            &blocked.pressure,
            &blocked.velocity,
            &radiated.pressure,
            &radiated.velocity,
        ] {
            for &value in values {
                amplitude = amplitude.max(value.abs());
            }
        }
        if !amplitude.is_finite() {
            return Err(bad("molecular incident work scale overflow").into());
        }
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        let normalize = |value: C64| -> Result<C64, MolecularModesError> {
            let scaled = divide(value, amplitude);
            if !finite(scaled)
                || (value.re != 0.0 && scaled.re == 0.0)
                || (value.im != 0.0 && scaled.im == 0.0)
                || (scaled != C64::ZERO && scaled.norm_sq() == 0.0)
            {
                Err(bad("molecular incident work normalization is unrepresentable").into())
            } else {
                Ok(scaled)
            }
        };
        let velocity: Vec<_> = motion
            .modal_velocity
            .iter()
            .map(|&v| normalize(v))
            .collect::<Result<_, _>>()?;
        // external, acoustic, Kelvin, Maxwell, blocked/motion/cross/total
        // exterior, blocked/motion/cross/total wall, apparent scale.
        let mut sums = [Sum::default(); 13];
        let m = velocity.len();
        for i in 0..m {
            poll(cancelled)?;
            let (f, a, v) = (
                normalize(external[i])?,
                normalize(acoustic[i])?,
                velocity[i],
            );
            sums[0].add(0.5 * (f * v.conj()).re);
            sums[1].add(0.5 * (a * v.conj()).re);
            sums[12].add(0.5 * (f.abs() + a.abs()) * v.abs());
            for (j, &other) in velocity.iter().enumerate() {
                let c = self.source.viscosity[i * m + j];
                sums[2].add(0.5 * (v.conj() * other).re * c);
                sums[12].add(0.5 * v.abs() * other.abs() * c.abs());
            }
        }
        for (&memory, &tau) in motion
            .elastic_memory_sqrt_j
            .iter()
            .zip(&self.source.memory.times)
        {
            poll(cancelled)?;
            let heat = 0.5 * normalize(memory)?.norm_sq() / tau;
            sums[3].add(heat);
            sums[12].add(heat);
        }
        let rho_c = self.source.medium.density * self.source.medium.sound_speed;
        for (panel, &area) in self.source.surface.areas().iter().enumerate() {
            if panel % 128 == 0 {
                poll(cancelled)?;
            }
            let pb = normalize(blocked.pressure[panel])?;
            let vb = normalize(blocked.velocity[panel])?;
            let pm = normalize(radiated.pressure[panel])?;
            let vm = normalize(radiated.velocity[panel])?;
            let pt = normalize(total.pressure[panel])?;
            let vt = normalize(total.velocity[panel])?;
            let y = beta.map_or(C64::ZERO, |beta| divide(beta[panel], rho_c));
            let (yb, ym, yt) = (y * pb, y * pm, y * pt);
            let half_area = 0.5 * area;
            for (channel, term) in [
                (4, pb * vb.conj()),
                (5, pm * vm.conj()),
                (6, pb * vm.conj()),
                (6, pm * vb.conj()),
                (7, pt * vt.conj()),
                (8, pb * yb.conj()),
                (9, pm * ym.conj()),
                (10, pb * ym.conj()),
                (10, pm * yb.conj()),
                (11, pt * yt.conj()),
            ] {
                if !finite(term) {
                    return Err(bad("molecular incident boundary work overflow").into());
                }
                sums[channel].add(half_area * term.re);
            }
            sums[12].add(
                half_area * (pb.abs() + pm.abs()) * (vb.abs() + vm.abs() + yb.abs() + ym.abs()),
            );
        }
        let values = sums.map(Sum::finish);
        let balance = values[0] - values[2] - values[3] - values[7] - values[11];
        let relative = if values[12] == 0.0 {
            balance.abs()
        } else {
            balance.abs() / values[12]
        };
        if values.iter().any(|value| !value.is_finite()) || !relative.is_finite() {
            return Err(bad("molecular incident total work overflow").into());
        }
        if relative > self.options.power_relative_tolerance {
            return Err(BemModesError::PowerBalance {
                relative,
                tolerance: self.options.power_relative_tolerance,
            }
            .into());
        }
        let restore = |value: f64| -> Result<f64, MolecularModesError> {
            let scaled = (value * amplitude) * amplitude;
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                Err(bad("physical molecular incident power is unrepresentable").into())
            } else {
                Ok(scaled)
            }
        };
        poll(cancelled)?;
        Ok(MolecularIncidentPower {
            external_mechanical_input_w: restore(values[0])?,
            blocked_acoustic_input_w: restore(values[1])?,
            kelvin_w: restore(values[2])?,
            relaxation_w: restore(values[3])?,
            blocked_exterior_w: restore(values[4])?,
            motion_exterior_w: restore(values[5])?,
            interference_exterior_w: restore(values[6])?,
            total_exterior_w: restore(values[7])?,
            blocked_wall_w: restore(values[8])?,
            motion_wall_w: restore(values[9])?,
            interference_wall_w: restore(values[10])?,
            total_wall_w: restore(values[11])?,
            balance_w: restore(balance)?,
            apparent_power_scale_w: restore(values[12])?,
            relative_balance: relative,
        })
    }
}

#[cfg(test)]
mod sphere_tests;
#[cfg(test)]
mod tests;
