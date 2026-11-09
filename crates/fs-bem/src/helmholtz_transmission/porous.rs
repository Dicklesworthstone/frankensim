//! Finite rigid-frame porous bodies on the original closed transmission mesh.
//!
//! The original JCA/JCAL owner supplies complex dynamic density and bulk
//! modulus, including optional coupled molecular/pore thermal relaxation.
//! Conjugate its exp(+iwt) response ONCE into this owner's exp(-iwt) convention.
//! Pressure and superficial normal velocity are continuous: the porosity is
//! already included in the effective coefficients and is never applied twice.
//! The full closed surface resolves front, rear and rim interactions together;
//! no locally reacting or infinite-planar transmission coefficient is fitted.
//!
//! The parent two-trace formulation retains its original P0 discretization.
//! Its density ratio is now complex, so q_i=i*omega*rho_eff*v differs physically
//! from a molecular gas's real inertial density. Adaptive receivers use this
//! same momentum factor for both Green traces and the pressure gradient.
//!
//! The single-interface owner admits one homogeneous isotropic, stationary,
//! rigid-frame inclusion; `multiple` couples distinct bodies using the existing
//! multiple-scattering engine. No Biot/frame displacement, mean flow, pore-
//! resolved field, or thin-layer accuracy is implied. Refinement must resolve BOTH the
//! physical thickness and complex wavelength; the wavelength gate alone does
//! not certify close-panel or thickness accuracy. Source validity remains
//! with the original material owner. Equivalent-fluid model reference:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.182.html>.
use super::*;
use crate::ComplexSum;
use crate::helmholtz::laplace::DampedWaveNumber;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::porous::{PorousFluidLaw, PorousFluidResponse};
use fs_material::relaxation::RelaxationState;

#[path = "porous/receivers.rs"]
mod receivers;
pub use receivers::PorousAbsorptionObservation;

/// Mutually interacting finite porous, lossless-fluid and molecular-gas bodies.
#[path = "porous/multiple.rs"]
pub mod multiple;

/// Original saturating-fluid owner, retained for the lifetime of the solve.
/// A relaxing gas uses the porous owner's coupled thermal equation; its free
/// gas attenuation is not added a second time.
#[derive(Clone, Copy)]
pub enum PoreFluid<'a> {
    Fluid(&'a dyn ThermoviscousFluid),
    RelaxingGas(&'a RelaxationState),
}

/// One original porous skeleton, coherent fluid, geometry and frequency.
pub struct PorousFluidInterface<'a> {
    surface: &'a SpherePanels,
    exterior_source: &'a dyn FluidAcoustics,
    exterior: FluidAcousticState,
    law: &'a PorousFluidLaw,
    filling: PoreFluid<'a>,
    response: PorousFluidResponse,
    omega: f64,
    k_o: DampedWaveNumber,
    k_i: DampedWaveNumber,
    density: C64,
    bulk: C64,
    contrast: Vec<C64>,
    factor: LuComplex,
    condition_lower_bound: f64,
    ppw: f64,
}

impl<'a> PorousFluidInterface<'a> {
    /// Prepare the material and one reusable coupled factorization. The
    /// exterior is the original lossless acoustic fluid; the interior is the
    /// original passive equivalent fluid. Source-band and constitutive
    /// refusals propagate before dense allocation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        exterior_source: &'a dyn FluidAcoustics,
        law: &'a PorousFluidLaw,
        filling: PoreFluid<'a>,
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
            return Err(bad("invalid porous transmission frequency or panel budget").into());
        }
        let n = surface.areas().len();
        if n > options.maximum_panels {
            return Err(HelmholtzError::WorkCap { panels: n }.into());
        }
        surface
            .validate()
            .map_err(|_| bad("invalid porous interface surface"))?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("porous interfaces require retained triangles"))?;
        if near_field::components(triangles)?.len() != 1 {
            return Err(bad("porous interface requires one closed outward component").into());
        }
        let workers = admit_storage(n, options)?;
        poll(gate)?;
        let exterior = exterior_source.acoustic_state()?;
        poll(gate)?;
        let response = match filling {
            PoreFluid::Fluid(fluid) => law.effective_fluid(fluid, omega_rad_s)?,
            PoreFluid::RelaxingGas(state) => law.effective_relaxing(state, omega_rad_s)?,
        };
        poll(gate)?;
        let density = response.dynamic_density_kg_m3().conj();
        let bulk = response.dynamic_bulk_modulus_pa().conj();
        let k_reference = omega_rad_s / exterior.sound_speed_m_s();
        let k_o = DampedWaveNumber {
            k: C64::from_re(k_reference),
            phase_origin_m: 0.0,
        };
        let k_i = DampedWaveNumber {
            k: response.characteristic().1.conj(),
            phase_origin_m: 0.0,
        };
        let rho_o = exterior.density_kg_m3();
        let eta = density.scale(1.0 / rho_o);
        // Preserve a weak contrast before adding/subtracting a rounded unity.
        let delta_rho = (density - C64::from_re(rho_o)).scale(1.0 / rho_o);
        let momentum = density.scale(omega_rad_s);
        let thermal_loss = bulk.recip().im * (0.5 * omega_rad_s);
        if !finite(density)
            || density.re <= 0.0
            || density.im < 0.0
            || !finite(bulk)
            || bulk.re <= 0.0
            || bulk.im > 0.0
            || !finite(k_i.k)
            || k_i.k.re <= 0.0
            || k_i.k.im < 0.0
            || !finite(eta)
            || !finite(delta_rho)
            || !finite(momentum)
            || momentum.re <= 0.0
            || !thermal_loss.is_finite()
            || thermal_loss < 0.0
            || [
                k_reference,
                rho_o * exterior.sound_speed_m_s(),
                omega_rad_s * rho_o,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(bad("unrepresentable or nonpassive porous transmission response").into());
        }
        let ppw =
            admit_resolution(surface, k_reference)?.min(admit_resolution(surface, k_i.k.abs())?);
        let matrix = assemble_wave(
            surface,
            k_o,
            k_i,
            k_reference,
            eta,
            delta_rho,
            gate,
            molecular::kernels::difference,
            molecular::kernels::self_contrast,
        )?;
        let (contrast, factor, condition_lower_bound) =
            factor_contrast(matrix, 2 * n, workers, gate)?;
        poll(gate)?;
        Ok(Self {
            surface,
            exterior_source,
            exterior,
            law,
            filling,
            response,
            omega: omega_rad_s,
            k_o,
            k_i,
            density,
            bulk,
            contrast,
            factor,
            condition_lower_bound,
            ppw,
        })
    }

    pub fn exterior_source(&self) -> &'a dyn FluidAcoustics {
        self.exterior_source
    }
    pub fn porous_source(&self) -> &'a PorousFluidLaw {
        self.law
    }
    pub fn pore_fluid_source(&self) -> PoreFluid<'a> {
        self.filling
    }
    /// Original material sample, with its ORIGINAL exp(+iwt) convention.
    pub fn material_response(&self) -> &PorousFluidResponse {
        &self.response
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Outgoing exp(-iwt) interior wavenumber [rad/m].
    pub fn interior_wavenumber_per_m(&self) -> C64 {
        self.k_i.k
    }
    /// Effective exp(-iwt) inertia, already including the open-area convention.
    pub fn interior_dynamic_density_kg_m3(&self) -> C64 {
        self.density
    }
    pub fn interior_dynamic_bulk_modulus_pa(&self) -> C64 {
        self.bulk
    }
    pub fn panels_per_wavelength(&self) -> f64 {
        self.ppw
    }
    /// Original probe-based diagnostic, not a condition-number certificate.
    pub fn condition_lower_bound(&self) -> f64 {
        self.condition_lower_bound
    }

    /// Reuse the factorization for a coherent analytic incident field. Finite
    /// sources must be exterior and regular throughout the porous volume.
    /// No source strength or interference term is normalized away.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<PorousTransmissionField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate()?;
        let n = self.surface.areas().len();
        molecular::admit_products(incident.source_count(), n, maximum_source_point_products)?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad("porous transmission source clearance must be positive"));
        }
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new(self.surface, &positions, minimum_source_clearance_m)?;
        }
        let mut incoming_pressure = Vec::with_capacity(n);
        let mut incoming_velocity = Vec::with_capacity(n);
        let momentum = self.exterior.density_kg_m3() * self.omega;
        for (&point, &normal) in self.surface.centroids().iter().zip(self.surface.normals()) {
            poll(gate)?;
            let sample = incident
                .sample_wave(self.k_o, momentum, None, point, &mut || gate.is_requested())?;
            incoming_pressure.push(sample.pressure_pa);
            let mut velocity = ComplexSum::default();
            for (component, direction) in sample.velocity_m_s.into_iter().zip(normal) {
                velocity.add(component.scale(direction));
            }
            incoming_velocity.push(velocity.finish());
        }
        let z = self.exterior.density_kg_m3() * self.exterior.sound_speed_m_s();
        let mut incoming = incoming_pressure.clone();
        incoming.extend(incoming_velocity.iter().map(|&v| v * C64::new(0.0, z)));
        if incoming.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite porous transmission incident trace"));
        }
        let mut scattered = vec![C64::ZERO; 2 * n];
        for (i, value) in scattered.iter_mut().enumerate() {
            poll(gate)?;
            let mut sum = ComplexSum::default();
            for (j, &b) in incoming.iter().enumerate() {
                sum.add(self.contrast[i * 2 * n + j] * b);
            }
            *value = -sum.finish();
        }
        if scattered.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite porous transmission contrast drive"));
        }
        self.factor
            .solve_with_cancel(&mut scattered, || gate.is_requested())
            .map_err(lu_error)?;
        let scattered_pressure = scattered[..n].to_vec();
        let scattered_velocity: Vec<_> = scattered[n..]
            .iter()
            .map(|&q| C64::new(q.im / z, -q.re / z))
            .collect();
        let total_pressure: Vec<_> = scattered_pressure
            .iter()
            .zip(incoming_pressure)
            .map(|(&a, b)| a + b)
            .collect();
        let total_velocity: Vec<_> = scattered_velocity
            .iter()
            .zip(incoming_velocity)
            .map(|(&a, b)| a + b)
            .collect();
        if scattered_pressure
            .iter()
            .chain(&scattered_velocity)
            .chain(&total_pressure)
            .chain(&total_velocity)
            .any(|&value| !finite(value))
        {
            return Err(bad("nonfinite porous transmission traces"));
        }
        let (_, interval, interface_work) = evaluate_radiated_power(
            &total_pressure,
            &total_velocity,
            self.surface.areas(),
            medium(self.exterior),
        );
        if !finite(interface_work.surface_power) || !interface_work.apparent_power.is_finite() {
            return Err(bad("nonfinite porous transmission boundary work"));
        }
        poll(gate)?;
        Ok(PorousTransmissionField {
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

/// A coherent solution retaining its original body, material and illumination.
pub struct PorousTransmissionField<'p, 's> {
    interface: &'p PorousFluidInterface<'s>,
    incident: IncidentField<'p>,
    maximum_source_point_products: u64,
    scattered_pressure: Vec<C64>,
    scattered_velocity: Vec<C64>,
    total_pressure: Vec<C64>,
    total_velocity: Vec<C64>,
    /// Signed outward work through the actual interface. In the continuum its
    /// negative real part equals irreversible loss within this source-free
    /// porous volume. Discretization defects remain signed and are not clipped
    /// into absorption or equated to remote scattered radiation power.
    pub interface_work: RadiationPowerDiagnostics,
    /// Roundoff in the panel dot product only, excluding spatial/material error.
    pub interface_work_roundoff_interval: (f64, f64),
}
impl PorousTransmissionField<'_, '_> {
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.total_pressure
    }
    /// Physical velocity per total exposed area: superficial velocity inside
    /// the porous body, matching exterior normal particle velocity.
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.total_velocity
    }
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.scattered_velocity
    }
}

#[cfg(test)]
#[path = "porous/tests.rs"]
mod tests;
