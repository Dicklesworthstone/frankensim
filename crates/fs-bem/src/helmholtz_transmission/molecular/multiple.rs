//! Mutual scattering by distinct molecular-gas bodies in one relaxing gas.
//!
//! The original multiple-body engine retains every reflection order, each
//! source-owned caloric law and one reusable coupled factorization. Complex
//! propagation changes the Green kernels, never the real rho*omega momentum.
//! Complete scattered traces include neighbors; local outgoing coefficients
//! alone reconstruct the field without counting those neighbors repeatedly.
//! Bodies are stationary, disjoint and non-nested with equal mean pressures.
//! Centroid receiver quadrature requires refinement; no far-field power in an
//! absorbing exterior or general close-interaction accuracy is inferred.
use super::*;
use crate::helmholtz::transmission::multiple::engine;

/// One original closed body and its actual molecular caloric source.
#[derive(Clone, Copy)]
pub struct Body<'a> {
    pub surface: &'a SpherePanels,
    pub material: &'a RelaxationState,
}

struct PreparedBody<'a> {
    original: Body<'a>,
    first_panel: usize,
    first_unknown: usize,
    wave: DampedWaveNumber,
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

/// Source-owned gas interfaces with all mutual reflections solved together.
/// Original material laws are evaluated once per body at the shared frequency.
pub struct MolecularFluidBodies<'a> {
    bodies: Vec<PreparedBody<'a>>,
    surface: SpherePanels,
    exterior: &'a RelaxationState,
    omega: f64,
    wave: DampedWaveNumber,
    ppw: f64,
    transfer: Vec<C64>,
    factor: Option<LuComplex>,
    condition_lower_bound: f64,
}
impl<'a> MolecularFluidBodies<'a> {
    /// Admit the original total-panel/peak-payload bounds and separated vertex
    /// boxes before quadratic allocation. Both complex wavelengths constrain
    /// every mesh; rho and mean pressure remain the original gas state's.
    /// The logical numerical payload excludes receiver geometry scratch,
    /// allocator/runtime metadata and thread stacks.
    pub fn new(
        bodies: &[Body<'a>],
        exterior: &'a RelaxationState,
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
            return Err(
                bad("invalid molecular multiple-body frequency, clearance or budget").into(),
            );
        }
        let (count, workers) = engine::admit(
            bodies.iter().map(|body| body.surface),
            bodies.len(),
            core::mem::size_of::<PreparedBody<'_>>(),
            minimum_body_clearance_m,
            options,
            gate,
        )?;
        exterior.gas().validate_coherence()?;
        let outside = exterior.gas();
        let wave = outgoing_wave(exterior, omega_rad_s, gate)?;
        let k_reference = omega_rad_s / outside.sound_speed;
        if [
            k_reference,
            outside.density * outside.sound_speed,
            outside.density * omega_rad_s,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
        {
            return Err(bad("unrepresentable molecular multiple-body exterior momentum").into());
        }
        let mut prepared = Vec::with_capacity(bodies.len());
        let mut triangles = Vec::with_capacity(count);
        let mut first_panel = 0;
        let mut ppw = f64::INFINITY;
        for &body in bodies {
            poll(gate)?;
            body.material.gas().validate_coherence()?;
            let inside = body.material.gas();
            if inside.pressure != outside.pressure {
                return Err(bad("stationary molecular bodies need the same mean pressure").into());
            }
            let inside_wave = outgoing_wave(body.material, omega_rad_s, gate)?;
            let eta = inside.density / outside.density;
            let delta_rho = (inside.density - outside.density) / outside.density;
            if [eta, inside.density * omega_rad_s]
                .into_iter()
                .any(|x| !x.is_finite() || x <= 0.0)
                || !delta_rho.is_finite()
            {
                return Err(bad("unrepresentable molecular body momentum or contrast").into());
            }
            ppw = ppw.min(
                admit_resolution(body.surface, wave.k.abs())?
                    .min(admit_resolution(body.surface, inside_wave.k.abs())?),
            );
            let contrast = if inside.density == outside.density && inside_wave.k == wave.k {
                Vec::new()
            } else {
                assemble_wave(
                    body.surface,
                    wave,
                    inside_wave,
                    k_reference,
                    eta,
                    delta_rho,
                    gate,
                    kernels::difference,
                    kernels::self_contrast,
                )?
            };
            prepared.push(PreparedBody {
                original: body,
                first_panel,
                first_unknown: 2 * first_panel,
                wave: inside_wave,
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
            .map_err(|_| bad("invalid assembled molecular multiple-body surface"))?;
        let (transfer, factor, condition_lower_bound) =
            engine::prepare(&prepared, wave, k_reference, 2 * count, workers, gate)?;
        poll(gate)?;
        Ok(Self {
            bodies: prepared,
            surface,
            exterior,
            omega: omega_rad_s,
            wave,
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
    pub fn exterior_source(&self) -> &'a RelaxationState {
        self.exterior
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn exterior_wavenumber_per_m(&self) -> C64 {
        self.wave.k
    }
    pub fn interior_wavenumber_per_m(&self, index: usize) -> Option<C64> {
        self.bodies.get(index).map(|body| body.wave.k)
    }
    pub fn panels_per_wavelength(&self) -> f64 {
        self.ppw
    }
    /// Probe-based lower bound from the original coupled factorization.
    pub fn condition_lower_bound(&self) -> f64 {
        self.condition_lower_bound
    }

    /// Sample the original coherent plane/volume/force field once and solve all
    /// reflections using the cached factor. Sources must be exterior to every
    /// body. The source-point budget applies here and to each total-field query.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<MolecularMultipleField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate()?;
        let n = self.surface.areas().len();
        admit_products(incident.source_count(), n, maximum_source_point_products)?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad(
                "molecular multiple-body source clearance must be positive",
            ));
        }
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new(&self.surface, &positions, minimum_source_clearance_m)?;
        }
        let mut pressure = Vec::with_capacity(n);
        let mut velocity = Vec::with_capacity(n);
        let momentum = self.exterior.gas().density * self.omega;
        for (&point, &normal) in self.surface.centroids().iter().zip(self.surface.normals()) {
            poll(gate)?;
            let sample = incident.sample_wave(self.wave, momentum, None, point, &mut || {
                gate.is_requested()
            })?;
            pressure.push(sample.pressure_pa);
            let mut normal_velocity = ComplexSum::default();
            for (v, direction) in sample.velocity_m_s.into_iter().zip(normal) {
                normal_velocity.add(v.scale(direction));
            }
            velocity.push(normal_velocity.finish());
        }
        let traces = engine::solve(
            &self.bodies,
            &self.surface,
            Medium::from(self.exterior.gas()),
            self.exterior.gas().density * self.exterior.gas().sound_speed,
            &self.transfer,
            self.factor.as_ref(),
            &pressure,
            &velocity,
            gate,
        )?;
        poll(gate)?;
        Ok(MolecularMultipleField {
            owner: self,
            incident,
            maximum_source_point_products,
            traces,
        })
    }
}

fn outgoing_wave(
    source: &RelaxationState,
    omega: f64,
    gate: &CancelGate,
) -> Result<DampedWaveNumber, TransmissionError> {
    poll(gate)?;
    let k = source.relaxing_wavenumber(omega)?;
    if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
        return Err(bad("molecular bodies need finite passive outgoing waves").into());
    }
    Ok(DampedWaveNumber {
        k,
        phase_origin_m: 0.0,
    })
}

/// One coherent field, local outgoing coefficients and complete physical traces.
pub struct MolecularMultipleField<'p, 's> {
    owner: &'p MolecularFluidBodies<'s>,
    incident: IncidentField<'p>,
    maximum_source_point_products: u64,
    traces: engine::Traces,
}
impl MolecularMultipleField<'_, '_> {
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.traces.total_pressure
    }
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.total_velocity
    }
    /// Complete scattered traces include every other body's regular field.
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.traces.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.scattered_velocity
    }
    /// Local outgoing coefficients exclude other bodies; use these only in
    /// the summed Green representation, never as the physical scattered work.
    pub fn outgoing_pressure_pa(&self) -> &[C64] {
        &self.traces.pressure
    }
    pub fn outgoing_normal_velocity_m_s(&self) -> &[C64] {
        &self.traces.velocity
    }
    /// Signed complete scattered flux including interference in a lossy gas.
    /// This is surface work, not power propagated to infinity.
    pub fn scattered_work(&self) -> &RadiationPowerDiagnostics {
        &self.traces.scattered_work
    }
    pub fn scattered_work_roundoff_interval(&self) -> (f64, f64) {
        self.traces.scattered_work_interval
    }
    /// Signed total work for each original body. An absorbing interior has
    /// inward work; numerical real defects remain visible, without clipping.
    /// The intervals enclose panel-dot roundoff, not discretization errors.
    pub fn interface_work(&self) -> &[(RadiationPowerDiagnostics, (f64, f64))] {
        &self.traces.work
    }

    /// Exterior scattered p/u from the original centroid Green derivatives.
    /// The bound counts panel-point products; each pair needs four kernels.
    pub fn scattered_pressure_velocity_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        admit_receivers(
            self.owner.surface.areas().len(),
            points,
            maximum_panel_point_products,
            gate,
        )?;
        near_field::Geometry::new(&self.owner.surface, points, minimum_clearance_m)?;
        self.outgoing_at(points, gate)
    }

    /// Coherent exterior p/u with the original incident field added once.
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
        let mut fields = self.scattered_pressure_velocity_at(
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            gate,
        )?;
        self.add_incident(points, &mut fields, gate)?;
        Ok(fields)
    }

    /// Transmitted p/u inside the specified body with its original rho*omega.
    /// An exactly transparent body still receives all active neighbors. Its
    /// field is evaluated in the common exterior gas without a fictitious
    /// local source or an approximate total-boundary representation.
    pub fn interior_pressure_velocity_at(
        &self,
        body: usize,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        poll(gate)?;
        let body = self
            .owner
            .bodies
            .get(body)
            .ok_or_else(|| bad("unknown molecular body"))?;
        let transparent = body.contrast.is_empty();
        let panel_count = if transparent {
            self.owner.surface.areas().len()
        } else {
            body.original.surface.areas().len()
        };
        admit_receivers(panel_count, points, maximum_panel_point_products, gate)?;
        near_field::Geometry::new_interior(body.original.surface, points, minimum_clearance_m)?;
        if transparent {
            admit_products(
                self.incident.source_count(),
                points.len(),
                self.maximum_source_point_products,
            )?;
            for other in &self.owner.bodies {
                poll(gate)?;
                if other.first_panel != body.first_panel {
                    near_field::Geometry::new(other.original.surface, points, minimum_clearance_m)?;
                }
            }
            let mut fields = self.outgoing_at(points, gate)?;
            self.add_incident(points, &mut fields, gate)?;
            return Ok(fields);
        }
        let range = body.first_panel..body.first_panel + body.original.surface.areas().len();
        samples_from_traces(
            body.original.surface,
            &self.traces.total_pressure[range.clone()],
            &self.traces.total_velocity[range],
            body.wave,
            body.original.material.gas().density * self.owner.omega,
            -1.0,
            points,
            gate,
        )
    }

    fn outgoing_at(
        &self,
        points: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<Vec<IncidentSample>, HelmholtzError> {
        if self.owner.factor.is_none() {
            poll(gate)?;
            return Ok(vec![
                IncidentSample {
                    pressure_pa: C64::ZERO,
                    velocity_m_s: [C64::ZERO; 3]
                };
                points.len()
            ]);
        }
        samples_from_traces(
            &self.owner.surface,
            &self.traces.pressure,
            &self.traces.velocity,
            self.owner.wave,
            self.owner.exterior.gas().density * self.owner.omega,
            1.0,
            points,
            gate,
        )
    }

    fn add_incident(
        &self,
        points: &[[f64; 3]],
        fields: &mut [IncidentSample],
        gate: &CancelGate,
    ) -> Result<(), HelmholtzError> {
        let owner = self.owner;
        for (&point, field) in points.iter().zip(fields) {
            poll(gate)?;
            let incident = self.incident.sample_wave(
                owner.wave,
                owner.exterior.gas().density * owner.omega,
                None,
                point,
                &mut || gate.is_requested(),
            )?;
            // Exact transparency retains the incident values, including zero
            // components; no divide/re-multiply through the pressure scale.
            if owner.factor.is_none() {
                *field = incident;
            } else {
                field.pressure_pa = field.pressure_pa + incident.pressure_pa;
                for (v, u) in field.velocity_m_s.iter_mut().zip(incident.velocity_m_s) {
                    *v = *v + u;
                }
                if !finite(field.pressure_pa) || field.velocity_m_s.iter().any(|&v| !finite(v)) {
                    return Err(bad("nonfinite total molecular multiple-body receiver"));
                }
            }
        }
        poll(gate)
    }
}

fn admit_receivers(
    panels: usize,
    points: &[[f64; 3]],
    maximum: u64,
    gate: &CancelGate,
) -> Result<(), HelmholtzError> {
    poll(gate)?;
    admit_products(panels, points.len(), maximum)?;
    if points.is_empty() || points.len() > near_field::MAX_RECEIVERS {
        return Err(bad(
            "molecular multiple-body receiver count exceeds admission",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn samples_from_traces(
    surface: &SpherePanels,
    pressure: &[C64],
    velocity: &[C64],
    wave: DampedWaveNumber,
    momentum: f64,
    sign: f64,
    points: &[[f64; 3]],
    gate: &CancelGate,
) -> Result<Vec<IncidentSample>, HelmholtzError> {
    let mut result = Vec::with_capacity(points.len());
    for &point in points {
        poll(gate)?;
        let field = fields_from_traces(
            surface,
            pressure,
            velocity,
            wave,
            C64::new(0.0, momentum),
            &[point],
        )?[0];
        let sample = IncidentSample {
            pressure_pa: field.pressure_pa.scale(sign),
            velocity_m_s: field
                .gradient_pa_per_m
                .map(|g| g * C64::new(0.0, -sign / momentum)),
        };
        if !finite(sample.pressure_pa) || sample.velocity_m_s.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite physical molecular multiple-body receiver"));
        }
        result.push(sample);
    }
    poll(gate)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
