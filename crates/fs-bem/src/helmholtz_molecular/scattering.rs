//! Coherent illumination and diffraction in the original relaxing exterior.
//!
//! Finite sources and receiver fields use the same complex Green function as
//! the boundary operator. The source strengths retain their physical volume
//! flow [m³/s] and force on fluid [N], and momentum remains grad(p)=i*rho*omega*u.
//! The original lossless incident sampler shares these analytic kernels.
use super::*;
use crate::ComplexSum;
use crate::helmholtz::bc::incident::{IncidentField, IncidentSample};

/// Boundary assembly limits plus explicit analytic-source work and clearance.
#[derive(Debug, Clone, Copy)]
pub struct MolecularScatteringOptions {
    pub boundary: ComplexOptions,
    /// Minimum distance from each finite emitter to the closed body [m].
    pub minimum_source_clearance_m: f64,
    /// Maximum source terms times sampled points. The whole incident batch
    /// counts against this cap on the boundary; each later receiver query is
    /// admitted against the same cap. This excludes original BEM assembly work.
    pub maximum_source_point_products: u64,
}

/// Outgoing scattered traces and coherent total boundary work.
#[derive(Debug, Clone)]
pub struct MolecularScatteringSolution {
    /// Scattered field only; never substitute total traces in an outgoing
    /// Green representation, since the incident field is regular in the body.
    pub scattered: MolecularRadiation,
    pub incident_pressure_pa: Vec<C64>,
    pub incident_normal_velocity_m_s: Vec<C64>,
    /// Absorption from TOTAL pressure, including incident/scattered cross terms.
    pub absorbed_power_w: f64,
    /// Signed total outward boundary work. A stationary passive surface gives
    /// minus the absorbed wall power, up to floating-point cancellation.
    pub total_boundary_work_w: f64,
    pub declared_active: bool,
}

impl MolecularScatteringSolution {
    pub fn total_pressure_pa(&self) -> Vec<C64> {
        self.scattered
            .pressure
            .iter()
            .zip(&self.incident_pressure_pa)
            .map(|(&scattered, &incident)| scattered + incident)
            .collect()
    }

    pub fn total_normal_velocity_m_s(&self) -> Vec<C64> {
        self.scattered
            .velocity
            .iter()
            .zip(&self.incident_normal_velocity_m_s)
            .map(|(&scattered, &incident)| scattered + incident)
            .collect()
    }
}

/// One immutable body, incident field, original material and solved frequency.
/// Incident sources remain borrowed, so their phases cannot change under a
/// retained solution. All closed components are solved together.
pub struct MolecularScatteredField<'s, 'i> {
    surface: &'s SpherePanels,
    incident: IncidentField<'i>,
    solution: MolecularScatteringSolution,
    maximum_source_point_products: u64,
}

fn admit_source_work(terms: usize, points: usize, maximum: u64) -> Result<(), MolecularError> {
    if maximum == 0
        || (terms as u64)
            .checked_mul(points as u64)
            .is_none_or(|work| work > maximum)
    {
        return Err(bad(
            "molecular incident source-point work exceeds its declared budget",
        ));
    }
    Ok(())
}

impl MolecularExterior {
    fn incident_wave(&self, omega: f64) -> Result<DampedWaveNumber, MolecularError> {
        if !omega.is_finite() || omega < self.band_rad_s[0] || omega > self.band_rad_s[1] {
            return Err(bad(
                "molecular incident frequency is outside its declared band",
            ));
        }
        let k = self.source.relaxing_wavenumber(omega)?;
        let momentum = self.source.gas().density * omega;
        if !finite(k) || k.re <= 0.0 || k.im < 0.0 || !momentum.is_finite() || momentum <= 0.0 {
            return Err(bad(
                "molecular incident wave or momentum scale is unrepresentable",
            ));
        }
        Ok(DampedWaveNumber {
            k,
            phase_origin_m: 0.0,
        })
    }

    /// Exact coherent free-space pressure and velocity from the original
    /// plane-wave, monopole and force-dipole source declarations. A plane wave's
    /// amplitude is its pressure at the coordinate origin; it decays along its
    /// propagation direction and grows upstream in this unbounded idealization.
    /// Singular/unrepresentable receivers refuse. Body admission belongs to
    /// `solve_scattering_batch` below.
    pub fn incident_pressure_velocity_at(
        &self,
        omega: f64,
        field: IncidentField<'_>,
        point_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IncidentSample, MolecularError> {
        poll(&mut cancelled)?;
        field.validate()?;
        let wave = self.incident_wave(omega)?;
        let out = field.sample_wave(
            wave,
            self.source.gas().density * omega,
            None,
            point_m,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(out)
    }

    /// Scatter coherent plane waves, finite monopoles and compact force sources
    /// off rigid or local/nonlocal impedance bodies in the SAME molecular gas.
    ///
    /// A stationary body satisfies v_total=-Y*p_total. The original Robin
    /// owner therefore receives v_drive=-v_incident-Y*p_incident; both of its
    /// boundary equations keep q=i*rho*omega*v. All illuminations share ONE
    /// assembly and LU. Source geometry and cardinality, total source-panel
    /// work, the original complex operator budgets and cancellation are admitted
    /// before solving. A mutable nonlocal-law fingerprint refuses publication.
    ///
    /// Total wall absorption includes coherent cross terms. Scattered outward
    /// work also feeds gas absorption; it is not a far-field scattering cross
    /// section, nor the incident power removed from an absorbing exterior.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_scattering_batch<'s, 'i>(
        &self,
        surface: &'s SpherePanels,
        omega: f64,
        incidents: &[IncidentField<'i>],
        admittance: Admittance<'_>,
        formulation: Formulation,
        options: &MolecularScatteringOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularScatteredField<'s, 'i>>, MolecularError> {
        poll(&mut cancelled)?;
        if incidents.is_empty() || incidents.len() > MAX_RADIATION_FIELDS_PER_BATCH {
            return Err(bad(
                "molecular scattering requires a nonempty bounded incident batch",
            ));
        }
        if !options.minimum_source_clearance_m.is_finite()
            || options.minimum_source_clearance_m <= 0.0
        {
            return Err(bad(
                "molecular incident source clearance must be positive and finite",
            ));
        }
        complex::admit_values_storage(surface, formulation, &options.boundary)?;
        super::super::bc::validate_closed(surface)?;
        let n = surface.areas().len();
        let mut terms = 0usize;
        for &field in incidents {
            field.validate()?;
            terms = terms
                .checked_add(field.source_count())
                .ok_or_else(|| bad("molecular incident source count overflow"))?;
        }
        admit_source_work(terms, n, options.maximum_source_point_products)?;
        admittance.validate(n)?;
        let stamp = admittance.stamp();
        let wave = self.incident_wave(omega)?;
        let gas = self.source.gas();
        let momentum = gas.density * omega;
        let rho_c = gas.density * gas.sound_speed;
        if !rho_c.is_finite() || rho_c <= 0.0 {
            return Err(bad(
                "molecular incident admittance normalization is unrepresentable",
            ));
        }
        // Source locations must be regular in every closed component. These
        // are the original winding/triangle-clearance checks, not a sphere proxy.
        for &field in incidents {
            poll(&mut cancelled)?;
            let positions = field.source_positions();
            if !positions.is_empty() {
                crate::near_field::admit_exterior_points(
                    surface,
                    &positions,
                    options.minimum_source_clearance_m,
                    MAX_DENSE_PANELS,
                )?;
            }
        }
        let mut traces = Vec::with_capacity(incidents.len());
        let mut drives = Vec::with_capacity(incidents.len());
        for &field in incidents {
            let mut pressure = Vec::with_capacity(n);
            let mut velocity = Vec::with_capacity(n);
            for (&point, &normal) in surface.centroids().iter().zip(surface.normals()) {
                poll(&mut cancelled)?;
                let sample = field.sample_wave(wave, momentum, None, point, &mut cancelled)?;
                pressure.push(sample.pressure_pa);
                let mut normal_velocity = ComplexSum::default();
                for (value, direction) in sample.velocity_m_s.iter().zip(normal) {
                    normal_velocity.add(value.scale(direction));
                }
                velocity.push(normal_velocity.finish());
            }
            let wall = admittance.image(&pressure)?;
            let drive: Vec<_> = velocity
                .iter()
                .zip(wall)
                .map(|(&incident, absorbed)| -incident - absorbed.scale(1.0 / rho_c))
                .collect();
            if drive.iter().any(|&v| !finite(v)) {
                return Err(bad("molecular incident Robin drive is nonfinite"));
            }
            traces.push((pressure, velocity));
            drives.push(drive);
        }
        if admittance.stamp() != stamp {
            return Err(bad(
                "molecular scattering admittance changed during source preparation",
            ));
        }
        let fields: Vec<_> = drives.iter().map(Vec::as_slice).collect();
        let solved = self.solve_boundary_batch(
            surface,
            omega,
            &fields,
            admittance,
            formulation,
            &options.boundary,
            &mut cancelled,
        )?;
        let mut out = Vec::with_capacity(incidents.len());
        for ((solved, (pressure, velocity)), &incident) in
            solved.into_iter().zip(traces).zip(incidents)
        {
            poll(&mut cancelled)?;
            let scattered = solved.radiation;
            let total: Vec<_> = pressure
                .iter()
                .zip(&scattered.pressure)
                .map(|(&a, &b)| a + b)
                .collect();
            let wall = admittance.image(&total)?;
            let mut absorbed_work = ComplexSum::default();
            let mut total_work = ComplexSum::default();
            let mut apparent = 0.0;
            for i in 0..n {
                if i % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                let p = total[i];
                let v = velocity[i] + scattered.velocity[i];
                let inward = wall[i].scale(1.0 / rho_c);
                if !finite(p) || !finite(v) || !finite(inward) {
                    return Err(bad("nonfinite molecular total boundary field"));
                }
                let half_area = 0.5 * surface.areas()[i];
                absorbed_work.add((p * inward.conj()).scale(half_area));
                apparent += p.abs() * inward.abs() * half_area;
                total_work.add((p * v.conj()).scale(half_area));
            }
            let absorbed_work = absorbed_work.finish();
            let total_work = total_work.finish();
            let absorbed_power_w = absorbed_work.re;
            let total_boundary_work_w = total_work.re;
            if !finite(absorbed_work)
                || !finite(total_work)
                || !apparent.is_finite()
                || (!admittance.declared_active()
                    && absorbed_power_w < -1e-10 * apparent.max(f64::MIN_POSITIVE))
            {
                return Err(bad(
                    "molecular total wall work is nonfinite or undeclared active",
                ));
            }
            out.push(MolecularScatteredField {
                surface,
                incident,
                solution: MolecularScatteringSolution {
                    scattered,
                    incident_pressure_pa: pressure,
                    incident_normal_velocity_m_s: velocity,
                    absorbed_power_w,
                    total_boundary_work_w,
                    declared_active: admittance.declared_active(),
                },
                maximum_source_point_products: options.maximum_source_point_products,
            });
        }
        if admittance.stamp() != stamp {
            return Err(bad("molecular scattering admittance changed during solve"));
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

impl MolecularScatteredField<'_, '_> {
    pub fn solution(&self) -> &MolecularScatteringSolution {
        &self.solution
    }

    /// Total incident plus scattered pressure/particle velocity using the
    /// original adaptive triangle receiver. Analytic incident fields are added
    /// once at the same physical points. Returned componentwise discrepancies
    /// cover scattered-field quadrature only, not solve or material error.
    pub fn total_pressure_velocity_adaptive_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_receiver_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::near_field::VectorObservation, MolecularError> {
        let incident = self.near_field_incident(points_m, &mut cancelled)?;
        let prepared = self.solution.scattered.prepare_near_field(
            self.surface,
            points_m,
            minimum_receiver_clearance_m,
            options,
            &mut cancelled,
        )?;
        self.combine_near_field(&prepared, incident, &mut cancelled)
    }

    /// Reuse one source-bound adaptive row bank across a coherent incident
    /// batch. Its retained receiver positions determine both scattered and
    /// analytic fields; callers cannot supply mismatched observation points.
    pub fn total_pressure_velocity_with_near_field(
        &self,
        prepared: &PreparedMolecularNearField<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::near_field::VectorObservation, MolecularError> {
        let incident = self.near_field_incident(prepared.points_m(), &mut cancelled)?;
        self.combine_near_field(prepared, incident, &mut cancelled)
    }

    fn near_field_incident(
        &self,
        points: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<IncidentSample>, MolecularError> {
        poll(cancelled)?;
        if points.is_empty() || points.len() > crate::near_field::MAX_RECEIVERS {
            return Err(bad(
                "molecular adaptive receiver array exceeds its cardinality limit",
            ));
        }
        admit_source_work(
            self.incident.source_count(),
            points.len(),
            self.maximum_source_point_products,
        )?;
        let radiation = &self.solution.scattered;
        let wave = DampedWaveNumber {
            k: radiation.k,
            phase_origin_m: 0.0,
        };
        let momentum = radiation.source.gas().density * radiation.omega;
        points
            .iter()
            .map(|&point| {
                self.incident
                    .sample_wave(wave, momentum, None, point, cancelled)
                    .map_err(MolecularError::from)
            })
            .collect()
    }

    fn combine_near_field(
        &self,
        prepared: &PreparedMolecularNearField<'_>,
        incident: Vec<IncidentSample>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<crate::near_field::VectorObservation, MolecularError> {
        let mut output = prepared.evaluate(&self.solution.scattered, &mut *cancelled)?;
        for (i, source) in incident.into_iter().enumerate() {
            poll(cancelled)?;
            let pressure = output.scalar.pressure[i] + source.pressure_pa;
            let velocity = std::array::from_fn(|axis| {
                output.particle_velocity_m_s[i][axis] + source.velocity_m_s[axis]
            });
            if !finite(pressure) || velocity.iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite adaptive molecular total receiver"));
            }
            output.scalar.pressure[i] = pressure;
            output.particle_velocity_m_s[i] = velocity;
        }
        poll(cancelled)?;
        Ok(output)
    }

    /// Coherent total p [Pa] and u [m/s] at finite exterior points. The incident
    /// analytic gradient and the original scattered Green derivative both enter
    /// inviscid momentum. Retains angular/near-field force-source velocity even
    /// at a pressure null. No plane-wave pressure/velocity approximation occurs.
    ///
    /// The existing molecular receiver uses centroid quadrature; its clearance
    /// gate is not a near-field spatial convergence certificate. Molecular gas
    /// absorption is already present, so no post-hoc attenuation is applied.
    pub fn total_pressure_velocity_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_receiver_clearance_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<IncidentSample>, MolecularError> {
        poll(&mut cancelled)?;
        if points_m.is_empty() || points_m.len() > crate::near_field::MAX_RECEIVERS {
            return Err(bad(
                "molecular scattering receiver array exceeds its cardinality limit",
            ));
        }
        admit_source_work(
            self.incident.source_count(),
            points_m.len(),
            self.maximum_source_point_products,
        )?;
        let radiation = &self.solution.scattered;
        let momentum = radiation.source.gas().density * radiation.omega;
        let wave = DampedWaveNumber {
            k: radiation.k,
            phase_origin_m: 0.0,
        };
        let mut out = Vec::with_capacity(points_m.len());
        // Singular source/receiver collisions refuse before BEM evaluation.
        for &point in points_m {
            out.push(
                self.incident
                    .sample_wave(wave, momentum, None, point, &mut cancelled)?,
            );
        }
        let scattered = radiation.observe(
            self.surface,
            points_m,
            minimum_receiver_clearance_m,
            &mut cancelled,
        )?;
        let inverse_momentum = C64::new(0.0, -1.0 / momentum);
        for (sample, scattered) in out.iter_mut().zip(scattered) {
            poll(&mut cancelled)?;
            sample.pressure_pa = sample.pressure_pa + scattered.pressure_pa;
            for (u, gradient) in sample
                .velocity_m_s
                .iter_mut()
                .zip(scattered.gradient_pa_per_m)
            {
                *u = *u + gradient * inverse_momentum;
            }
            if !finite(sample.pressure_pa) || sample.velocity_m_s.iter().any(|&u| !finite(u)) {
                return Err(bad(
                    "nonfinite molecular total receiver pressure or velocity",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

#[cfg(test)]
#[path = "scattering/sphere_tests.rs"]
mod sphere_tests;
#[cfg(test)]
#[path = "scattering/tests.rs"]
mod tests;
