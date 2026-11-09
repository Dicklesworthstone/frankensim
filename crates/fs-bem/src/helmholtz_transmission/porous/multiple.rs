//! Mutually interacting porous and fluid bodies in one lossless exterior.
//!
//! Each original material supplies its own local transmission contrast. The
//! existing multiple-body engine solves every mutual reflection in one cached
//! factorization; there is no independent-body sum or reflection-order cutoff.
//! Complex porous inertia is retained in both interface continuity and local
//! receivers. Fluid and molecular-gas bodies can share the same collection.
//!
//! Every body is one closed outward component. The original separated-box
//! admission excludes touching, nested and overlapping bodies. P0/centroid
//! cross-body integration requires spatial and separation refinement. Finite
//! thickness, complex wavelength and inter-body gaps must all be resolved.
//! The background is stationary and lossless; mechanical frame motion, Biot
//! coupling, mean flow and a general intersection solver are outside this API.
use super::*;
use crate::helmholtz::transmission::multiple::engine;

#[path = "multiple/receivers.rs"]
mod receivers;

/// The original constitutive owner for one body's interior.
#[derive(Clone, Copy)]
pub enum BodyMaterial<'a> {
    /// Rigid-frame equivalent fluid, with all porosity factors already included.
    Porous {
        law: &'a PorousFluidLaw,
        filling: PoreFluid<'a>,
    },
    /// Original lossless acoustic fluid. Repeated identical trait-object views
    /// share their first snapshot, following the existing fluid-body owner.
    Fluid(&'a dyn FluidAcoustics),
    /// Free molecular gas, with real inertia and its original complex bulk law.
    /// Its molecular attenuation is not added to a porous response separately.
    RelaxingGas(&'a RelaxationState),
}

/// Original closed mesh and source-owned interior material.
#[derive(Clone, Copy)]
pub struct Body<'a> {
    pub surface: &'a SpherePanels,
    pub material: BodyMaterial<'a>,
}

struct PreparedBody<'a> {
    original: Body<'a>,
    acoustic_state: Option<FluidAcousticState>,
    porous_response: Option<PorousFluidResponse>,
    density: C64,
    bulk: C64,
    wave: DampedWaveNumber,
    first_panel: usize,
    first_unknown: usize,
    contrast: Vec<C64>,
}
impl engine::BoundaryBlock for PreparedBody<'_> {
    fn surface(&self) -> &SpherePanels {
        self.original.surface
    }
    fn first_panel(&self) -> usize {
        self.first_panel
    }
    fn first_unknown(&self) -> usize {
        self.first_unknown
    }
    fn contrast(&self) -> &[C64] {
        &self.contrast
    }
}

/// Source-bound finite bodies with one reusable fully coupled factorization.
/// Traces retain body order, then the original panel order within each body.
pub struct PorousFluidBodies<'a> {
    bodies: Vec<PreparedBody<'a>>,
    surface: SpherePanels,
    exterior_source: &'a dyn FluidAcoustics,
    exterior: FluidAcousticState,
    omega: f64,
    k_o: f64,
    ppw: f64,
    transfer: Vec<C64>,
    factor: Option<LuComplex>,
    condition_lower_bound: f64,
}
impl<'a> PorousFluidBodies<'a> {
    /// Prepare original material responses at one shared positive frequency.
    /// A porous material's authored frequency range and thermodynamic admission
    /// remain with its original owner. The caller supplies a stationary mean
    /// equilibrium; generic acoustic fluid views do not expose mean pressure.
    ///
    /// The inherited peak-payload admission includes local contrast matrices,
    /// mutual transfer, coupled LU, copied geometry and these descriptors.
    /// Allocator/runtime metadata, stacks and receiver scratch remain excluded.
    pub fn new(
        bodies: &[Body<'a>],
        exterior_source: &'a dyn FluidAcoustics,
        omega_rad_s: f64,
        minimum_body_clearance_m: f64,
        options: Options,
        gate: &CancelGate,
    ) -> Result<Self, TransmissionError> {
        poll(gate)?;
        if bodies.is_empty()
            || !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || !minimum_body_clearance_m.is_finite()
            || minimum_body_clearance_m <= 0.0
            || options.maximum_panels == 0
            || options.maximum_panels > MAX_TRANSMISSION_PANELS
            || bodies.len() > options.maximum_panels
        {
            return Err(bad("invalid porous multiple-body frequency, clearance or budget").into());
        }
        let (count, workers) = engine::admit(
            bodies.iter().map(|body| body.surface),
            bodies.len(),
            core::mem::size_of::<PreparedBody<'_>>(),
            minimum_body_clearance_m,
            options,
            gate,
        )?;
        let exterior = exterior_source.acoustic_state()?;
        let rho_o = exterior.density_kg_m3();
        let k_o = omega_rad_s / exterior.sound_speed_m_s();
        if [k_o, rho_o * omega_rad_s, rho_o * exterior.sound_speed_m_s()]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(bad("unrepresentable porous multiple-body exterior scale").into());
        }
        let outside_wave = DampedWaveNumber {
            k: C64::from_re(k_o),
            phase_origin_m: 0.0,
        };
        let mut prepared: Vec<PreparedBody<'a>> = Vec::with_capacity(bodies.len());
        let mut triangles = Vec::with_capacity(count);
        let mut first_panel = 0;
        let mut ppw = f64::INFINITY;
        for &body in bodies {
            poll(gate)?;
            let (acoustic_state, porous_response, density, bulk, k) = match body.material {
                BodyMaterial::Fluid(source) => {
                    let state = acoustic_snapshot_for_view(
                        source,
                        core::iter::once((exterior_source, exterior)).chain(
                            prepared.iter().filter_map(|previous| {
                                match previous.original.material {
                                    BodyMaterial::Fluid(original) => Some((
                                        original,
                                        previous.acoustic_state.expect("fluid snapshot"),
                                    )),
                                    _ => None,
                                }
                            }),
                        ),
                    )?;
                    (
                        Some(state),
                        None,
                        C64::from_re(state.density_kg_m3()),
                        C64::from_re(
                            state.density_kg_m3()
                                * state.sound_speed_m_s()
                                * state.sound_speed_m_s(),
                        ),
                        C64::from_re(omega_rad_s / state.sound_speed_m_s()),
                    )
                }
                BodyMaterial::Porous { law, filling } => {
                    let response = match filling {
                        PoreFluid::Fluid(fluid) => law.effective_fluid(fluid, omega_rad_s)?,
                        PoreFluid::RelaxingGas(state) => {
                            law.effective_relaxing(state, omega_rad_s)?
                        }
                    };
                    (
                        None,
                        Some(response),
                        response.dynamic_density_kg_m3().conj(),
                        response.dynamic_bulk_modulus_pa().conj(),
                        response.characteristic().1.conj(),
                    )
                }
                BodyMaterial::RelaxingGas(state) => {
                    state.gas().validate_coherence()?;
                    (
                        None,
                        None,
                        C64::from_re(state.gas().density),
                        state.caloric_response(omega_rad_s)?.bulk_modulus_pa,
                        state.relaxing_wavenumber(omega_rad_s)?,
                    )
                }
            };
            let eta = density.scale(1.0 / rho_o);
            let delta = (density - C64::from_re(rho_o)).scale(1.0 / rho_o);
            let momentum = density.scale(omega_rad_s);
            let compression_loss = 0.5 * omega_rad_s * bulk.recip().im;
            if !finite(density)
                || density.re <= 0.0
                || density.im < 0.0
                || !finite(bulk)
                || bulk.re <= 0.0
                || bulk.im > 0.0
                || !finite(k)
                || k.re <= 0.0
                || k.im < 0.0
                || !finite(eta)
                || !finite(delta)
                || !finite(momentum)
                || momentum.re <= 0.0
                || !compression_loss.is_finite()
                || compression_loss < 0.0
            {
                return Err(
                    bad("unrepresentable or nonpassive multiple-body material response").into(),
                );
            }
            let wave = DampedWaveNumber {
                k,
                phase_origin_m: 0.0,
            };
            ppw = ppw.min(
                admit_resolution(body.surface, k_o)?.min(admit_resolution(body.surface, k.abs())?),
            );
            let contrast = if let Some(state) = acoustic_state {
                if medium(state) == medium(exterior) {
                    Vec::new()
                } else {
                    // Preserve the original real-density assembly and arithmetic.
                    assemble(
                        body.surface,
                        k_o,
                        k.re,
                        state.density_kg_m3() / rho_o,
                        (state.density_kg_m3() - rho_o) / rho_o,
                        gate,
                    )?
                }
            } else if density == C64::from_re(rho_o) && k == outside_wave.k {
                Vec::new()
            } else {
                assemble_wave(
                    body.surface,
                    outside_wave,
                    wave,
                    k_o,
                    eta,
                    delta,
                    gate,
                    molecular::kernels::difference,
                    molecular::kernels::self_contrast,
                )?
            };
            prepared.push(PreparedBody {
                original: body,
                acoustic_state,
                porous_response,
                density,
                bulk,
                wave,
                first_panel,
                first_unknown: 2 * first_panel,
                contrast,
            });
            first_panel += body.surface.areas().len();
            triangles.extend_from_slice(
                body.surface
                    .triangles()
                    .expect("admitted original triangles"),
            );
        }
        let surface = SpherePanels::from_triangles(triangles)
            .map_err(|_| bad("invalid assembled porous multiple-body surface"))?;
        // The background is real: keep the existing real-k mutual arithmetic.
        let (transfer, factor, condition_lower_bound) =
            engine::prepare(&prepared, k_o, k_o, 2 * count, workers, gate)?;
        poll(gate)?;
        Ok(Self {
            bodies: prepared,
            surface,
            exterior_source,
            exterior,
            omega: omega_rad_s,
            k_o,
            ppw,
            transfer,
            factor,
            condition_lower_bound,
        })
    }

    pub fn surface(&self) -> &SpherePanels {
        &self.surface
    }
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    pub fn body(&self, index: usize) -> Option<Body<'a>> {
        self.bodies.get(index).map(|body| body.original)
    }
    pub fn exterior_source(&self) -> &'a dyn FluidAcoustics {
        self.exterior_source
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn exterior_wavenumber_per_m(&self) -> f64 {
        self.k_o
    }
    pub fn interior_wavenumber_per_m(&self, index: usize) -> Option<C64> {
        self.bodies.get(index).map(|body| body.wave.k)
    }
    pub fn interior_dynamic_density_kg_m3(&self, index: usize) -> Option<C64> {
        self.bodies.get(index).map(|body| body.density)
    }
    pub fn interior_dynamic_bulk_modulus_pa(&self, index: usize) -> Option<C64> {
        self.bodies.get(index).map(|body| body.bulk)
    }
    /// Original porous response in its ORIGINAL exp(+iwt) convention, or None
    /// when the selected body is a free fluid/gas or the index does not exist.
    pub fn porous_material_response(&self, index: usize) -> Option<&PorousFluidResponse> {
        self.bodies.get(index)?.porous_response.as_ref()
    }
    pub fn panels_per_wavelength(&self) -> f64 {
        self.ppw
    }
    /// Original probe-based lower bound, not a condition-number certificate.
    pub fn condition_lower_bound(&self) -> f64 {
        self.condition_lower_bound
    }

    /// Reuse the coupled factor for a coherent original plane/volume/force
    /// source field. Finite sources must lie outside every body's volume.
    /// Source-point admission remains distinct from panel receiver work.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<PorousMultipleField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate()?;
        let n = self.surface.areas().len();
        molecular::admit_products(incident.source_count(), n, maximum_source_point_products)?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad(
                "porous multiple-body source clearance must be positive",
            ));
        }
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new(&self.surface, &positions, minimum_source_clearance_m)?;
        }
        let mut pressure = Vec::with_capacity(n);
        let mut velocity = Vec::with_capacity(n);
        let momentum = self.exterior.density_kg_m3() * self.omega;
        for (&point, &normal) in self.surface.centroids().iter().zip(self.surface.normals()) {
            poll(gate)?;
            let sample = incident
                .sample_wave(self.k_o, momentum, None, point, &mut || gate.is_requested())?;
            pressure.push(sample.pressure_pa);
            let mut normal_velocity = ComplexSum::default();
            for (component, direction) in sample.velocity_m_s.into_iter().zip(normal) {
                normal_velocity.add(component.scale(direction));
            }
            velocity.push(normal_velocity.finish());
        }
        let traces = engine::solve(
            &self.bodies,
            &self.surface,
            medium(self.exterior),
            self.exterior.density_kg_m3() * self.exterior.sound_speed_m_s(),
            &self.transfer,
            self.factor.as_ref(),
            &pressure,
            &velocity,
            gate,
        )?;
        poll(gate)?;
        Ok(PorousMultipleField {
            owner: self,
            incident,
            maximum_source_point_products,
            traces,
        })
    }
}

/// Complete coherent solution retaining original body/material/source identity.
pub struct PorousMultipleField<'p, 's> {
    owner: &'p PorousFluidBodies<'s>,
    incident: IncidentField<'p>,
    maximum_source_point_products: u64,
    traces: engine::Traces,
}
impl PorousMultipleField<'_, '_> {
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.traces.total_pressure
    }
    /// Superficial velocity on porous faces, matching the exterior velocity
    /// per gross face area. Porosity is never multiplied a second time.
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.total_velocity
    }
    /// Complete physical scattered traces include every neighbor's field.
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.traces.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.scattered_velocity
    }
    /// Only these LOCAL outgoing coefficients belong in a summed Green field.
    pub fn outgoing_pressure_pa(&self) -> &[C64] {
        &self.traces.pressure
    }
    pub fn outgoing_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.velocity
    }
    /// Signed complete scattered surface work, retaining all mutual terms.
    pub fn scattered_work(&self) -> &RadiationPowerDiagnostics {
        &self.traces.scattered_work
    }
    pub fn scattered_work_roundoff_interval(&self) -> (f64, f64) {
        self.traces.scattered_work_interval
    }
    /// Signed outward total work for each original body in its original order.
    /// Negative real work is inward absorption in the continuum. Panel defects
    /// stay signed; the intervals cover dot-product roundoff only, not mesh,
    /// constitutive-model or close-interaction error.
    pub fn interface_work(&self) -> &[(RadiationPowerDiagnostics, (f64, f64))] {
        &self.traces.work
    }
}

#[cfg(test)]
#[path = "multiple/tests.rs"]
mod tests;
