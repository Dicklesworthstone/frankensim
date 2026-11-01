//! Mutually reflecting, distinct fluid bodies in one homogeneous exterior.
//!
//! For body i, C_i is the ORIGINAL single-body contrast operator, so its
//! outgoing traces s_i obey (I+C_i)s_i=-C_i h_i. The exciting field h_i is
//! the original incident field plus every other body's outgoing field:
//! h_i=inc_i+sum E_ij s_j. With pressure-scaled normal derivative q/k_o,
//! E_ij=[D_o,-k_o S_o; N_o/k_o,-D'_o]. One coupled factorization therefore
//! solves (I+C_i)s_i+sum C_i E_ij s_j=-C_i inc_i. All reflection orders are
//! retained; no independent-body sum or scattering-order truncation is used.
//!
//! This scattering formulation keeps exact material transparency: C_i=0
//! makes s_i exactly zero even beside a strongly scattering body. It also
//! retains the original weak-contrast forcing. Applying a global M-I to
//! sampled incident traces would instead leave discrete cancellation errors
//! from other bodies' regular incident representations.
//!
//! Each body has one closed outward, non-self-intersecting triangle component.
//! Pairwise separated vertex bounding boxes conservatively admit disjoint,
//! non-nested bodies. Overlapping boxes refuse, even if their surfaces might
//! be disjoint. This is a bounded geometry domain, not an intersection solver.
//! Cross-body kernels use the original constant-panel centroid quadrature;
//! close-interaction accuracy and continuum accuracy need mesh refinement.

use super::*;

/// Original surface and inviscid fluid source for one bounded body.
#[derive(Clone, Copy)]
pub struct Body<'a> {
    pub surface: &'a SpherePanels,
    pub material: &'a dyn FluidAcoustics,
}

struct PreparedBody<'a> {
    original: Body<'a>,
    state: FluidAcousticState,
    first_panel: usize,
    first_unknown: usize,
    k: f64,
    ppw: f64,
    contrast: Vec<C64>,
}

/// One complete dense multiple-scattering solve and retained source identities.
///
/// Unknowns are grouped by body, then pressure and q/k_o. The exterior Green
/// transfer and factorization are prepared once and reused across incidents.
/// Interior states are sampled once each; the shared exterior is sampled once.
pub struct FluidBodies<'a> {
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

impl<'a> FluidBodies<'a> {
    /// Prepare the coupled system after admitting all original wavelengths,
    /// the total panel count, geometric separation and peak numerical payload.
    /// A positive box clearance is a lower bound on surface separation. The
    /// logical byte budget includes copied surface arrays and original LU
    /// worker scratch, but excludes allocator/runtime metadata and stacks.
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
            return Err(
                bad("invalid multiple-fluid body, frequency, clearance or panel budget").into(),
            );
        }
        let mut count = 0usize;
        let mut local_entries = 0usize;
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(bodies.len());
        for body in bodies {
            poll(gate)?;
            let n = body.surface.areas().len();
            count = count
                .checked_add(n)
                .ok_or_else(|| bad("multiple-fluid panel overflow"))?;
            if count > options.maximum_panels {
                return Err(HelmholtzError::WorkCap { panels: count }.into());
            }
            body.surface
                .validate()
                .map_err(|_| bad("invalid multiple-fluid surface"))?;
            let triangles = body
                .surface
                .triangles()
                .ok_or_else(|| bad("multiple-fluid bodies require retained triangles"))?;
            if near_field::components(triangles)?.len() != 1 {
                return Err(bad("each fluid body requires one closed outward component").into());
            }
            let mut lower = [f64::INFINITY; 3];
            let mut upper = [f64::NEG_INFINITY; 3];
            for vertex in triangles.iter().flatten() {
                for axis in 0..3 {
                    lower[axis] = lower[axis].min(vertex[axis]);
                    upper[axis] = upper[axis].max(vertex[axis]);
                }
            }
            for &(lo, hi) in &boxes {
                let mut gap2 = 0.0;
                for axis in 0..3 {
                    let gap = (lower[axis] - hi[axis])
                        .max(lo[axis] - upper[axis])
                        .max(0.0);
                    gap2 += gap * gap;
                }
                let scale = lower
                    .iter()
                    .chain(&upper)
                    .chain(&lo)
                    .chain(&hi)
                    .fold(0.0_f64, |a, &b| a.max(b.abs()));
                if gap2.sqrt() < minimum_body_clearance_m.max(128.0 * f64::EPSILON * scale) {
                    return Err(
                        bad("fluid body bounding boxes overlap or violate separation").into(),
                    );
                }
            }
            boxes.push((lower, upper));
            local_entries += 4 * n * n;
        }
        drop(boxes);
        let m = 2 * count;
        let workers = DENSE_WORKERS.with(core::cell::Cell::get);
        let scratch = lu_complex_blocked_workspace_bytes(m, workers)
            .ok_or_else(|| bad("multiple-fluid LU workspace byte count overflow"))?;
        let complex_bytes = core::mem::size_of::<C64>();
        let transfer_entries = if bodies.len() == 1 { 0 } else { m * m };
        // Factorization clones A; retained C_i and E are live throughout.
        // The solve peak includes incident, solution, exciting, output traces.
        let fixed = (local_entries + transfer_entries) * complex_bytes
            + m * core::mem::size_of::<usize>()
            + count * 16 * core::mem::size_of::<f64>()
            + bodies.len() * core::mem::size_of::<PreparedBody<'_>>();
        let factor_peak = fixed
            .checked_add((2 * m * m + m) * complex_bytes)
            .and_then(|x| x.checked_add(scratch))
            .ok_or_else(|| bad("multiple-fluid factor byte count overflow"))?;
        let solve_peak = fixed
            .checked_add((m * m + 8 * m) * complex_bytes)
            .and_then(|x| {
                x.checked_add(
                    bodies.len() * core::mem::size_of::<(RadiationPowerDiagnostics, (f64, f64))>(),
                )
            })
            .ok_or_else(|| bad("multiple-fluid solve byte count overflow"))?;
        let required_bytes = factor_peak.max(solve_peak);
        if required_bytes > options.maximum_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes,
                maximum_bytes: options.maximum_bytes,
            }
            .into());
        }
        let exterior = exterior_source.acoustic_state()?;
        let k_o = omega_rad_s / exterior.sound_speed_m_s();
        if !k_o.is_finite()
            || k_o <= 0.0
            || !(omega_rad_s * exterior.density_kg_m3()).is_finite()
            || omega_rad_s * exterior.density_kg_m3() <= 0.0
        {
            return Err(bad("unrepresentable multiple-fluid exterior scale").into());
        }
        let mut prepared = Vec::with_capacity(bodies.len());
        let mut triangles = Vec::with_capacity(count);
        let mut first_panel = 0;
        let mut ppw = f64::INFINITY;
        for &body in bodies {
            poll(gate)?;
            let state = body.material.acoustic_state()?;
            let k = omega_rad_s / state.sound_speed_m_s();
            let eta = state.density_kg_m3() / exterior.density_kg_m3();
            let delta_rho =
                (state.density_kg_m3() - exterior.density_kg_m3()) / exterior.density_kg_m3();
            if [k, eta, omega_rad_s * state.density_kg_m3()]
                .into_iter()
                .any(|v| !v.is_finite() || v <= 0.0)
                || !delta_rho.is_finite()
            {
                return Err(
                    bad("unrepresentable multiple-fluid interior contrast or scale").into(),
                );
            }
            let body_ppw =
                admit_resolution(body.surface, k_o)?.min(admit_resolution(body.surface, k)?);
            ppw = ppw.min(body_ppw);
            let contrast = if medium(state) == medium(exterior) {
                Vec::new()
            } else {
                assemble(body.surface, k_o, k, eta, delta_rho, gate)?
            };
            prepared.push(PreparedBody {
                original: body,
                state,
                first_panel,
                first_unknown: 2 * first_panel,
                k,
                ppw: body_ppw,
                contrast,
            });
            first_panel += body.surface.areas().len();
            triangles.extend_from_slice(body.surface.triangles().expect("admitted triangles"));
        }
        let surface = SpherePanels::from_triangles(triangles)
            .map_err(|_| bad("invalid assembled multiple-fluid surface"))?;
        let transparent = prepared.iter().all(|body| body.contrast.is_empty());
        let (transfer, factor, condition_lower_bound) = if transparent {
            (Vec::new(), None, 1.0)
        } else {
            let transfer = if prepared.len() == 1 {
                Vec::new()
            } else {
                assemble_transfer(&prepared, k_o, m, gate)?
            };
            let mut matrix = zero_matrix(m * m, "multiple-fluid matrix")?;
            for body in &prepared {
                let mi = 2 * body.original.surface.areas().len();
                for row in 0..mi {
                    poll(gate)?;
                    let global_row = body.first_unknown + row;
                    matrix[global_row * m + global_row] = C64::ONE;
                    if body.contrast.is_empty() {
                        continue;
                    }
                    for col in 0..mi {
                        matrix[global_row * m + body.first_unknown + col] = matrix
                            [global_row * m + body.first_unknown + col]
                            + body.contrast[row * mi + col];
                    }
                    for other in &prepared {
                        // An identical-fluid body's exact equation is s_j=0.
                        // Eliminate those zero columns algebraically; its
                        // interior still receives every active body's field.
                        if other.first_unknown == body.first_unknown || other.contrast.is_empty() {
                            continue;
                        }
                        let mj = 2 * other.original.surface.areas().len();
                        for col in 0..mj {
                            if col % 64 == 0 {
                                poll(gate)?;
                            }
                            let global_col = other.first_unknown + col;
                            let mut value = C64::ZERO;
                            for inner in 0..mi {
                                value = value
                                    + body.contrast[row * mi + inner]
                                        * transfer[(body.first_unknown + inner) * m + global_col];
                            }
                            matrix[global_row * m + global_col] = value;
                        }
                    }
                }
            }
            if matrix.iter().any(|&z| !finite(z)) {
                return Err(bad("nonfinite multiple-fluid coupled operator").into());
            }
            let factor =
                lu_complex_blocked_with_cancel(&matrix, m, workers, || gate.is_requested())
                    .map_err(lu_error)?;
            poll(gate)?;
            let condition = condition_lower_bound(&matrix, m, &factor);
            if !condition.is_finite() {
                return Err(bad("nonfinite multiple-fluid condition diagnostic").into());
            }
            (transfer, Some(factor), condition)
        };
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

    /// Combined original triangles, in input-body order, for existing observers.
    pub fn surface(&self) -> &SpherePanels {
        &self.surface
    }
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    pub fn exterior_source(&self) -> &'a dyn FluidAcoustics {
        self.exterior_source
    }
    pub fn exterior_state(&self) -> FluidAcousticState {
        self.exterior
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn exterior_wavenumber_rad_m(&self) -> f64 {
        self.k_o
    }
    pub fn body(&self, index: usize) -> Option<Body<'a>> {
        self.bodies.get(index).map(|body| body.original)
    }
    pub fn interior_state(&self, index: usize) -> Option<FluidAcousticState> {
        self.bodies.get(index).map(|body| body.state)
    }

    /// Solve every mutual reflection together, reusing the prepared factor.
    /// Caller-authored traces must extend regularly through every body.
    pub fn scatter<'p>(
        &'p self,
        incident: &IncidentTraces,
        gate: &CancelGate,
    ) -> Result<Solution<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate_binding(&self.surface, self.k_o, medium(self.exterior))?;
        let n = self.surface.areas().len();
        let m = 2 * n;
        let z = self.exterior.characteristic_impedance_pa_s_m();
        let mut incoming = vec![C64::ZERO; m];
        for body in &self.bodies {
            let ni = body.original.surface.areas().len();
            for j in 0..ni {
                incoming[body.first_unknown + j] = incident.pressure_pa()[body.first_panel + j];
                incoming[body.first_unknown + ni + j] =
                    incident.normal_velocity_m_s()[body.first_panel + j] * C64::new(0.0, z);
            }
        }
        if incoming.iter().any(|&value| !finite(value)) {
            return Err(bad(
                "nonfinite pressure-scaled multiple-fluid incident field",
            ));
        }
        let mut outgoing = vec![C64::ZERO; m];
        if let Some(factor) = &self.factor {
            for body in &self.bodies {
                if body.contrast.is_empty() {
                    continue;
                }
                let mi = 2 * body.original.surface.areas().len();
                for row in 0..mi {
                    poll(gate)?;
                    for col in 0..mi {
                        outgoing[body.first_unknown + row] = outgoing[body.first_unknown + row]
                            - body.contrast[row * mi + col] * incoming[body.first_unknown + col];
                    }
                }
            }
            if outgoing.iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite multiple-fluid contrast drive"));
            }
            factor
                .solve_with_cancel(&mut outgoing, || gate.is_requested())
                .map_err(lu_error)?;
        }
        let mut total_scattered = outgoing.clone();
        if !self.transfer.is_empty() {
            for row in 0..m {
                poll(gate)?;
                for (col, &value) in outgoing.iter().enumerate() {
                    total_scattered[row] =
                        total_scattered[row] + self.transfer[row * m + col] * value;
                }
            }
        }
        let mut pressure = vec![C64::ZERO; n];
        let mut velocity = vec![C64::ZERO; n];
        let mut scattered_pressure = vec![C64::ZERO; n];
        let mut scattered_velocity = vec![C64::ZERO; n];
        let mut total_pressure = vec![C64::ZERO; n];
        let mut total_velocity = vec![C64::ZERO; n];
        let mut work = Vec::with_capacity(self.bodies.len());
        for body in &self.bodies {
            poll(gate)?;
            let ni = body.original.surface.areas().len();
            for j in 0..ni {
                let (p, q) = (body.first_unknown + j, body.first_unknown + ni + j);
                let target = body.first_panel + j;
                pressure[target] = outgoing[p];
                velocity[target] = if body.contrast.is_empty() {
                    C64::ZERO
                } else {
                    C64::new(outgoing[q].im / z, -outgoing[q].re / z)
                };
                scattered_pressure[target] = total_scattered[p];
                scattered_velocity[target] = if self.factor.is_none() {
                    C64::ZERO
                } else {
                    C64::new(total_scattered[q].im / z, -total_scattered[q].re / z)
                };
                if self.factor.is_none() {
                    total_pressure[target] = incident.pressure_pa()[target];
                    total_velocity[target] = incident.normal_velocity_m_s()[target];
                } else {
                    total_pressure[target] = total_scattered[p] + incoming[p];
                    let total_q = total_scattered[q] + incoming[q];
                    total_velocity[target] = C64::new(total_q.im / z, -total_q.re / z);
                }
            }
            let range = body.first_panel..body.first_panel + ni;
            let (_, interval, diagnostics) = evaluate_radiated_power(
                &total_pressure[range.clone()],
                &total_velocity[range],
                body.original.surface.areas(),
                medium(self.exterior),
            );
            if !finite(diagnostics.surface_power) || !diagnostics.apparent_power.is_finite() {
                return Err(bad("nonfinite multiple-fluid interface work"));
            }
            work.push((diagnostics, interval));
        }
        if pressure
            .iter()
            .chain(&velocity)
            .chain(&total_pressure)
            .chain(&total_velocity)
            .chain(&scattered_pressure)
            .chain(&scattered_velocity)
            .any(|&v| !finite(v))
        {
            return Err(bad("nonfinite multiple-fluid solution traces"));
        }
        // Local outgoing densities reconstruct the field without adding each
        // neighbor's regular representation repeatedly. Physical flux uses
        // COMPLETE scattered boundary traces, including interference terms.
        let outgoing = representation(
            &self.surface,
            self.k_o,
            medium(self.exterior),
            self.ppw,
            self.condition_lower_bound,
            pressure,
            velocity,
        )?;
        let (_, scattered_work_interval, scattered_work) = evaluate_radiated_power(
            &scattered_pressure,
            &scattered_velocity,
            self.surface.areas(),
            medium(self.exterior),
        );
        if !finite(scattered_work.surface_power) || !scattered_work.apparent_power.is_finite() {
            return Err(bad("nonfinite complete multiple-fluid scattered work"));
        }
        poll(gate)?;
        Ok(Solution {
            interface: self,
            outgoing,
            scattered_pressure,
            scattered_velocity,
            scattered_work,
            scattered_work_interval,
            total_pressure,
            total_velocity,
            work,
        })
    }

    /// Sample and retain one coherent incident field regular in every body.
    pub fn scatter_field<'p>(
        &'p self,
        field: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        gate: &CancelGate,
    ) -> Result<ScatteredField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        let traces = IncidentTraces::from_field(
            &self.surface,
            self.k_o,
            medium(self.exterior),
            field,
            minimum_source_clearance_m,
        )?;
        let solution = self.scatter(&traces, gate)?;
        Ok(ScatteredField {
            incident: field,
            solution,
        })
    }
}

/// Outgoing contributions and total physical traces, ordered by original body.
pub struct Solution<'p, 's> {
    interface: &'p FluidBodies<'s>,
    outgoing: RadiationSolution,
    scattered_pressure: Vec<C64>,
    scattered_velocity: Vec<C64>,
    scattered_work: RadiationPowerDiagnostics,
    scattered_work_interval: (f64, f64),
    total_pressure: Vec<C64>,
    total_velocity: Vec<C64>,
    work: Vec<(RadiationPowerDiagnostics, (f64, f64))>,
}
impl Solution<'_, '_> {
    /// Complete scattered field on every body, including all other bodies.
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.scattered_velocity
    }
    /// Per-body outgoing representation coefficients. These exclude the
    /// other bodies' fields and are NOT the complete scattered boundary data.
    pub fn outgoing_pressure_pa(&self) -> &[C64] {
        &self.outgoing.pressure
    }
    pub fn outgoing_normal_velocity_m_s(&self) -> &[C64] {
        &self.outgoing.velocity
    }
    /// Signed physical scattered flux from COMPLETE boundary traces. Includes
    /// mutual interference; no sum of isolated outgoing powers or clipping.
    pub fn scattered_work(&self) -> &RadiationPowerDiagnostics {
        &self.scattered_work
    }
    pub fn scattered_work_roundoff_interval(&self) -> (f64, f64) {
        self.scattered_work_interval
    }
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.total_pressure
    }
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.total_velocity
    }
    /// Signed total work per body; real discretization defects are not absorption.
    /// Each interval bounds panel-dot roundoff only, not continuum error.
    pub fn interface_work(&self) -> &[(RadiationPowerDiagnostics, (f64, f64))] {
        &self.work
    }

    /// Original retained-triangle far-field quadrature of the summed outgoing
    /// representation. Multiply by exp(i*k*r)/r once for far-field pressure.
    pub fn far_field(
        &self,
        directions: &[[f64; 3]],
        maximum_kernel_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<Vec<C64>, HelmholtzError> {
        super::super::far_field_integrated_cancellable(
            &self.interface.surface,
            &self.outgoing,
            medium(self.interface.exterior),
            directions,
            maximum_kernel_evaluations,
            || gate.is_requested(),
        )
    }

    /// Scattered pressure outside every body, with original receiver controls.
    pub fn scattered_pressure_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let interface = self.interface;
        near_field::Geometry::new(&interface.surface, points_m, minimum_clearance_m)?
            .prepare(interface.k_o, medium(interface.exterior), options)?
            .evaluate(&self.outgoing)
    }

    /// Transmitted pressure in the specified original body, from all mutual
    /// interactions. Original winding/clearance and adaptive quadrature apply.
    pub fn interior_pressure_at(
        &self,
        body: usize,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let body = self
            .interface
            .bodies
            .get(body)
            .ok_or_else(|| bad("unknown fluid body"))?;
        let geometry = near_field::Geometry::new_interior(
            body.original.surface,
            points_m,
            minimum_clearance_m,
        )?;
        let range = body.first_panel..body.first_panel + body.original.surface.areas().len();
        let traces = representation(
            body.original.surface,
            body.k,
            medium(body.state),
            body.ppw,
            self.interface.condition_lower_bound,
            self.total_pressure[range.clone()]
                .iter()
                .map(|&p| -p)
                .collect(),
            self.total_velocity[range].iter().map(|&v| -v).collect(),
        )?;
        geometry
            .prepare(body.k, medium(body.state), options)?
            .evaluate(&traces)
    }
}

/// A coherent incident field and the coupled outgoing solution.
pub struct ScatteredField<'p, 's> {
    incident: IncidentField<'p>,
    solution: Solution<'p, 's>,
}
impl ScatteredField<'_, '_> {
    pub fn solution(&self) -> &Solution<'_, '_> {
        &self.solution
    }
    pub fn total_pressure_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let interface = self.solution.interface;
        let mut observed =
            self.solution
                .scattered_pressure_at(points_m, minimum_clearance_m, options)?;
        for (&point, pressure) in points_m.iter().zip(&mut observed.pressure) {
            *pressure = *pressure
                + self
                    .incident
                    .pressure_velocity_at(interface.k_o, medium(interface.exterior), point)?
                    .pressure_pa;
            if !finite(*pressure) {
                return Err(bad("nonfinite multiple-fluid exterior pressure"));
            }
        }
        Ok(observed)
    }
}

fn zero_matrix(entries: usize, what: &'static str) -> Result<Vec<C64>, HelmholtzError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(entries)
        .map_err(|_| HelmholtzError::AllocationFailed { what })?;
    values.resize(entries, C64::ZERO);
    Ok(values)
}

fn assemble_transfer(
    bodies: &[PreparedBody<'_>],
    k: f64,
    m: usize,
    gate: &CancelGate,
) -> Result<Vec<C64>, HelmholtzError> {
    let mut transfer = zero_matrix(m * m, "multiple-fluid exterior transfer")?;
    for target in bodies {
        let a = target.original.surface;
        let ni = a.areas().len();
        for source in bodies {
            if source.first_unknown == target.first_unknown || source.contrast.is_empty() {
                continue;
            }
            let b = source.original.surface;
            let nj = b.areas().len();
            for i in 0..ni {
                poll(gate)?;
                for j in 0..nj {
                    let (s, d, dp, n) = kernels(
                        k,
                        a.centroids()[i],
                        a.normals()[i],
                        b.centroids()[j],
                        b.normals()[j],
                    );
                    let (r, c) = (target.first_unknown + i, source.first_unknown + j);
                    let area = b.areas()[j];
                    transfer[r * m + c] = d.scale(area);
                    transfer[r * m + c + nj] = s.scale(-k * area);
                    transfer[(r + ni) * m + c] = n.scale(area / k);
                    transfer[(r + ni) * m + c + nj] = dp.scale(-area);
                }
            }
        }
    }
    if transfer.iter().any(|&v| !finite(v)) {
        return Err(bad("nonfinite multiple-fluid exterior transfer"));
    }
    Ok(transfer)
}

fn representation(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    ppw: f64,
    condition: f64,
    pressure: Vec<C64>,
    velocity: Vec<C64>,
) -> Result<RadiationSolution, HelmholtzError> {
    let (power, interval, diagnostics) =
        evaluate_radiated_power(&pressure, &velocity, surface.areas(), medium);
    if !power.is_finite()
        || !finite(diagnostics.surface_power)
        || !diagnostics.apparent_power.is_finite()
    {
        return Err(bad("nonfinite multiple-fluid representation work"));
    }
    Ok(RadiationSolution {
        pressure,
        velocity,
        k,
        medium,
        surface_fingerprint: surface_fingerprint(surface),
        panels_per_wavelength: ppw,
        radiated_power: power,
        radiated_power_roundoff_interval: interval,
        power_diagnostics: diagnostics,
        condition_lower_bound: condition,
        dense_cap_utilization: 2.0 * surface.areas().len() as f64 / MAX_DENSE_PANELS as f64,
    })
}
