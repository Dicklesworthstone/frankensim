//! Source-owned molecular compression across one penetrable gas interface.
//!
//! The original direct transmission combination and dense factorization are
//! retained. Both outgoing wavenumbers come from their original caloric laws;
//! pressure and physical normal velocity remain continuous. The real reference
//! k_ref=omega/c_eq,out only scales the unknown q_out/k_ref=i*rho_out*c_eq,out*v.
//! Replacing this momentum factor with complex k*c_eq would change the problem.
//!
//! This is a stationary, source-free, homogeneous gas interior, inviscid
//! momentum and relaxing bulk compression. Both mean pressures must agree.
//! No molecular mixing, mean heat transport, shear viscosity, solid interface
//! or multiple-body solve is implied. Finite receivers use the original
//! centroid Green derivative, not the lossless owner's adaptive quadrature.
use super::*;
use crate::ComplexSum;
use crate::helmholtz::bc::incident::IncidentSample;
use crate::helmholtz::laplace::{DampedWaveNumber, fields_from_traces};
use fs_material::relaxation::RelaxationState;

mod kernels;

/// One original molecular gas on each side, with one reusable factorization.
pub struct MolecularFluidInterface<'a> {
    surface: &'a SpherePanels,
    exterior: &'a RelaxationState,
    interior: &'a RelaxationState,
    omega: f64,
    k_o: DampedWaveNumber,
    k_i: DampedWaveNumber,
    contrast: Vec<C64>,
    factor: Option<LuComplex>,
    condition_lower_bound: f64,
    ppw: f64,
}

impl<'a> MolecularFluidInterface<'a> {
    /// Evaluate each original reservoir law once at the actual frequency.
    /// Admit both |k| wavelengths and the original worker-dependent dense
    /// payload BEFORE allocating the contrast matrix or LU. Identical density
    /// and complex k are exactly transparent, including molecular attenuation.
    pub fn new(
        surface: &'a SpherePanels,
        exterior: &'a RelaxationState,
        interior: &'a RelaxationState,
        omega_rad_s: f64,
        options: Options,
        gate: &CancelGate,
    ) -> Result<Self, TransmissionError> {
        poll(gate)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || options.maximum_panels == 0
            || options.maximum_panels > MAX_TRANSMISSION_PANELS
        {
            return Err(bad("invalid molecular transmission frequency or panel budget").into());
        }
        let n = surface.areas().len();
        if n > options.maximum_panels {
            return Err(HelmholtzError::WorkCap { panels: n }.into());
        }
        surface
            .validate()
            .map_err(|_| bad("invalid molecular interface surface"))?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("molecular interfaces require retained triangles"))?;
        if near_field::components(triangles)?.len() != 1 {
            return Err(bad("molecular interface requires one closed outward component").into());
        }
        let workers = admit_storage(n, options)?;
        exterior.gas().validate_coherence()?;
        interior.gas().validate_coherence()?;
        if exterior.gas().pressure != interior.gas().pressure {
            return Err(
                bad("a stationary molecular gas interface needs equal mean pressures").into(),
            );
        }
        poll(gate)?;
        let k_o = DampedWaveNumber {
            k: exterior.relaxing_wavenumber(omega_rad_s)?,
            phase_origin_m: 0.0,
        };
        poll(gate)?;
        let k_i = DampedWaveNumber {
            k: interior.relaxing_wavenumber(omega_rad_s)?,
            phase_origin_m: 0.0,
        };
        for k in [k_o.k, k_i.k] {
            if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
                return Err(
                    bad("molecular transmission needs finite passive outgoing waves").into(),
                );
            }
        }
        let ppw =
            admit_resolution(surface, k_o.k.abs())?.min(admit_resolution(surface, k_i.k.abs())?);
        let outside = exterior.gas();
        let inside = interior.gas();
        let eta = inside.density / outside.density;
        let delta_rho = (inside.density - outside.density) / outside.density;
        let k_reference = omega_rad_s / outside.sound_speed;
        if [
            eta,
            k_reference,
            outside.density * outside.sound_speed,
            omega_rad_s * outside.density,
            omega_rad_s * inside.density,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
            || !delta_rho.is_finite()
        {
            return Err(bad("unrepresentable molecular transmission momentum or contrast").into());
        }
        let (contrast, factor, condition_lower_bound) =
            if outside.density == inside.density && k_o.k == k_i.k {
                (Vec::new(), None, 1.0)
            } else {
                let matrix = assemble_wave(
                    surface,
                    k_o,
                    k_i,
                    k_reference,
                    eta,
                    delta_rho,
                    gate,
                    kernels::difference,
                    kernels::self_contrast,
                )?;
                let (matrix, factor, condition) = factor_contrast(matrix, 2 * n, workers, gate)?;
                (matrix, Some(factor), condition)
            };
        poll(gate)?;
        Ok(Self {
            surface,
            exterior,
            interior,
            omega: omega_rad_s,
            k_o,
            k_i,
            contrast,
            factor,
            condition_lower_bound,
            ppw,
        })
    }

    pub fn exterior_source(&self) -> &'a RelaxationState {
        self.exterior
    }
    pub fn interior_source(&self) -> &'a RelaxationState {
        self.interior
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn exterior_wavenumber_per_m(&self) -> C64 {
        self.k_o.k
    }
    pub fn interior_wavenumber_per_m(&self) -> C64 {
        self.k_i.k
    }
    /// Original mesh admission diagnostic for the smaller complex wavelength.
    pub fn panels_per_wavelength(&self) -> f64 {
        self.ppw
    }
    /// Original probe-based condition lower bound; not a certificate.
    pub fn condition_lower_bound(&self) -> f64 {
        self.condition_lower_bound
    }

    /// Transmit a coherent field through the SAME molecular interface. All
    /// sources must remain exterior and regular throughout the bounded gas.
    /// The explicit source-point budget covers this preparation and each later
    /// total-field query; the factorization is reused for subsequent incidents.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<MolecularTransmissionField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate()?;
        let n = self.surface.areas().len();
        admit_products(incident.source_count(), n, maximum_source_point_products)?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad(
                "molecular transmission source clearance must be positive",
            ));
        }
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new(self.surface, &positions, minimum_source_clearance_m)?;
        }
        let mut incoming_pressure = Vec::with_capacity(n);
        let mut incoming_velocity = Vec::with_capacity(n);
        let momentum = self.exterior.gas().density * self.omega;
        for (&point, &normal) in self.surface.centroids().iter().zip(self.surface.normals()) {
            poll(gate)?;
            let sample = incident
                .sample_wave(self.k_o, momentum, None, point, &mut || gate.is_requested())?;
            incoming_pressure.push(sample.pressure_pa);
            let mut v = ComplexSum::default();
            for (component, direction) in sample.velocity_m_s.into_iter().zip(normal) {
                v.add(component.scale(direction));
            }
            incoming_velocity.push(v.finish());
        }
        let z = self.exterior.gas().density * self.exterior.gas().sound_speed;
        let mut incoming = incoming_pressure.clone();
        incoming.extend(incoming_velocity.iter().map(|&v| v * C64::new(0.0, z)));
        if incoming.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite molecular transmission incident trace"));
        }
        let mut scattered = vec![C64::ZERO; 2 * n];
        if let Some(factor) = &self.factor {
            for (i, value) in scattered.iter_mut().enumerate() {
                poll(gate)?;
                let mut sum = ComplexSum::default();
                for (j, &b) in incoming.iter().enumerate() {
                    sum.add(self.contrast[i * 2 * n + j] * b);
                }
                *value = -sum.finish();
            }
            if scattered.iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite molecular transmission contrast drive"));
            }
            factor
                .solve_with_cancel(&mut scattered, || gate.is_requested())
                .map_err(lu_error)?;
        }
        let scattered_pressure = scattered[..n].to_vec();
        let transparent = self.factor.is_none();
        let scattered_velocity = if transparent {
            vec![C64::ZERO; n]
        } else {
            scattered[n..]
                .iter()
                .map(|&q| C64::new(q.im / z, -q.re / z))
                .collect()
        };
        let total_pressure = if transparent {
            incoming_pressure
        } else {
            scattered_pressure
                .iter()
                .zip(incoming_pressure)
                .map(|(&a, b)| a + b)
                .collect()
        };
        let total_velocity = if transparent {
            incoming_velocity
        } else {
            scattered_velocity
                .iter()
                .zip(incoming_velocity)
                .map(|(&a, b)| a + b)
                .collect()
        };
        if scattered_pressure
            .iter()
            .chain(&scattered_velocity)
            .chain(&total_pressure)
            .chain(&total_velocity)
            .any(|&z| !finite(z))
        {
            return Err(bad("nonfinite molecular transmission traces"));
        }
        let (_, interval, interface_work) = evaluate_radiated_power(
            &total_pressure,
            &total_velocity,
            self.surface.areas(),
            Medium::from(self.exterior.gas()),
        );
        if !finite(interface_work.surface_power) || !interface_work.apparent_power.is_finite() {
            return Err(bad("nonfinite molecular transmission interface work"));
        }
        poll(gate)?;
        Ok(MolecularTransmissionField {
            interface: self,
            incident,
            maximum_source_point_products,
            scattered_pressure,
            scattered_velocity,
            total_pressure,
            total_velocity,
            interface_work,
            interface_work_roundoff_interval: interval,
        })
    }
}

/// Scattered outgoing traces and continuous physical traces, all in exp(-iwt).
pub struct MolecularTransmissionField<'p, 's> {
    interface: &'p MolecularFluidInterface<'s>,
    incident: IncidentField<'p>,
    maximum_source_point_products: u64,
    scattered_pressure: Vec<C64>,
    scattered_velocity: Vec<C64>,
    total_pressure: Vec<C64>,
    total_velocity: Vec<C64>,
    /// Signed OUTWARD work. Its negative real part is the continuum molecular
    /// dissipation inside the source-free gas. Numerical defects remain signed;
    /// this is neither clipped absorption nor remote exterior radiation power.
    pub interface_work: RadiationPowerDiagnostics,
    /// Panel-dot roundoff interval only, excluding spatial/material errors.
    pub interface_work_roundoff_interval: (f64, f64),
}

impl MolecularTransmissionField<'_, '_> {
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.total_pressure
    }
    /// One physical normal velocity, shared by both sides of the interface.
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.total_velocity
    }
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.scattered_velocity
    }

    /// Outgoing scattered p/u with the original centroid Green derivative.
    /// Budget counts panels times receivers (four kernel evaluations per pair).
    pub fn scattered_pressure_velocity_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        self.observe(
            false,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            gate,
        )
    }

    /// Total exterior pressure and velocity, retaining source interference and
    /// reactive/angular velocity. Molecular attenuation is already in k.
    pub fn total_pressure_velocity_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        poll(gate)?;
        admit_products(
            self.incident.source_count(),
            points.len(),
            self.maximum_source_point_products,
        )?;
        let owner = self.interface;
        let mut fields = self.scattered_pressure_velocity_at(
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            gate,
        )?;
        for (&point, field) in points.iter().zip(&mut fields) {
            let incident = self.incident.sample_wave(
                owner.k_o,
                owner.exterior.gas().density * owner.omega,
                None,
                point,
                &mut || gate.is_requested(),
            )?;
            field.pressure_pa = field.pressure_pa + incident.pressure_pa;
            for (v, u) in field.velocity_m_s.iter_mut().zip(incident.velocity_m_s) {
                *v = *v + u;
            }
            if !finite(field.pressure_pa) || field.velocity_m_s.iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite total molecular transmission receiver"));
            }
        }
        poll(gate)?;
        Ok(fields)
    }

    /// Transmitted interior p/u from S_i q_i-D_i p, with the original
    /// interior winding/triangle-clearance admission and actual rho_i*omega.
    /// No gas loss or pressure-to-velocity plane-wave approximation is added.
    pub fn interior_pressure_velocity_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        self.observe(
            true,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            gate,
        )
    }

    fn observe(
        &self,
        interior: bool,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        poll(gate)?;
        let owner = self.interface;
        admit_products(
            owner.surface.areas().len(),
            points.len(),
            maximum_panel_point_products,
        )?;
        if points.is_empty() || points.len() > near_field::MAX_RECEIVERS {
            return Err(bad(
                "molecular transmission receiver count exceeds admission",
            ));
        }
        if interior {
            near_field::Geometry::new_interior(owner.surface, points, minimum_clearance_m)?;
        } else {
            near_field::Geometry::new(owner.surface, points, minimum_clearance_m)?;
        }
        let (wave, gas, p, v, sign) = if interior {
            (
                owner.k_i,
                owner.interior.gas(),
                &self.total_pressure,
                &self.total_velocity,
                -1.0,
            )
        } else {
            (
                owner.k_o,
                owner.exterior.gas(),
                &self.scattered_pressure,
                &self.scattered_velocity,
                1.0,
            )
        };
        if interior && owner.factor.is_none() {
            admit_products(
                self.incident.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
            return points
                .iter()
                .map(|&point| {
                    self.incident.sample_wave(
                        wave,
                        gas.density * owner.omega,
                        None,
                        point,
                        &mut || gate.is_requested(),
                    )
                })
                .collect();
        }
        let momentum = gas.density * owner.omega;
        let mut fields = Vec::with_capacity(points.len());
        for &point in points {
            poll(gate)?;
            let field =
                fields_from_traces(owner.surface, p, v, wave, C64::new(0.0, momentum), &[point])?
                    [0];
            let sample = IncidentSample {
                pressure_pa: field.pressure_pa.scale(sign),
                velocity_m_s: field
                    .gradient_pa_per_m
                    .map(|g| g * C64::new(0.0, -sign / momentum)),
            };
            if !finite(sample.pressure_pa) || sample.velocity_m_s.iter().any(|&v| !finite(v)) {
                return Err(bad(
                    "nonfinite physical molecular transmission receiver velocity",
                ));
            }
            fields.push(sample);
        }
        poll(gate)?;
        Ok(fields)
    }
}

fn admit_products(terms: usize, points: usize, maximum: u64) -> Result<(), HelmholtzError> {
    if maximum == 0
        || (terms as u64)
            .checked_mul(points as u64)
            .is_none_or(|n| n > maximum)
    {
        return Err(bad(
            "molecular transmission source/receiver work exceeds its declared budget",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
