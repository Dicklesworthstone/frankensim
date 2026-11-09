//! Complete pressure FE fields coupled to finite-opening exterior radiation.
//!
//! Every original pressure DOF is retained. The same declared real P1 opening
//! traces supply the FE work columns and their exact dimensionless panel means
//! `psi` supply the original BEM velocity ports. Signed and zero-net-flow modes are
//! retained, with every self and mutual reaction. No plane-port conversion by
//! net flow, piston replacement, artificial damping or passivity repair occurs.
//!
//! The FE convention is exp(+i omega t); the native outgoing BEM matrix uses
//! exp(-i omega t). It is conjugated once and divided by the two geometric
//! support areas to obtain Pa s/m^3. FE flow Q points into the cavity,
//! exterior flow U points out, so U=-Q and the interface is P+Z Q=0. Positive
//! real power-wave references are numerical coordinates, not physical loads.
//! The native Robin-loaded FE factorization and a small port closure retain
//! the original loaded field, including at regular loaded-cavity resonances.
//! Original distributed volume injection and quadratic inward boundary motion
//! enter that same loaded field and its full exterior reaction coherently.
//! Their signed work is retained separately from driven-port work and losses.
//!
//! This bounded harmonic model has a rigid exterior except at the selected
//! openings. Other FE ports are explicitly driven interior feeds, not exterior
//! radiators. Each opening is planar and must use identical physical triangles
//! on both sides. The exterior uses their P0 velocity projection; pressure FE,
//! signed trace and BEM mesh/route convergence remain separate requirements.
//! The [`causal`] child constructs source-bound full-field exact-ZOH loads
//! from these same original opening traces and native BEM matrices.

/// Source-bound BEM bake, native rational realization and full-field time closure.
pub mod causal;

use super::{BoundaryLayerGeometry, GasState, HelmholtzFe, boundary_layer_transport};
use fs_bem::helmholtz::aperture::{
    Aperture, ApertureError, ApertureImpedance, ApertureRoute, aperture_impedance,
};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_exec::CancelGate;
use fs_feec::acoustics::helmholtz::AcousticFeError;
use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerPower;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousError, VolumeVelocitySource,
};
use fs_feec::acoustics::helmholtz::multiport::{
    BoundaryMode, InteriorMultiport, InteriorPortResponse, MAX_PORT_MODES,
    PreparedInteriorScattering,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;
use fs_math::det;
use std::collections::{BTreeMap, BTreeSet};

/// Default relative backward-error ceiling for the original physical port
/// constraints, measured in the declared real-reference power coordinates.
pub const DEFAULT_PORT_RELATIVE_TOLERANCE: f64 = 1e-8;

/// Explicit complete-field and exterior-work ceilings; exceeding one refuses.
#[derive(Clone, Copy, Debug)]
pub struct RadiationBudget {
    /// All original pressure DOFs, bounded by the FE owner's dense ceiling.
    pub maximum_pressure_dofs: usize,
    /// Actual closed exterior panels, before geometry-map allocation.
    pub maximum_panels: usize,
    /// Complete real opening-mode by exterior-panel trace bank.
    pub maximum_trace_entries: usize,
    /// `panels * radiating_modes^2` for one full pressure projection. Native
    /// BEM assembly/factorization retains its selected route's own limits.
    pub maximum_projection_terms: u64,
}

/// Original physical-owner refusals, invalid binding or cooperative cancellation.
#[derive(Debug)]
pub enum RadiationError {
    Invalid(&'static str),
    Interior(AcousticFeError),
    Source(HeterogeneousError),
    Geometry(HelmholtzError),
    Exterior(ApertureError),
    Closure(LuComplexError),
    Cancelled,
}
impl core::fmt::Display for RadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "interior/exterior radiation: {self:?}")
    }
}
impl std::error::Error for RadiationError {}
impl From<AcousticFeError> for RadiationError {
    fn from(error: AcousticFeError) -> Self {
        Self::Interior(error)
    }
}
impl From<HeterogeneousError> for RadiationError {
    fn from(error: HeterogeneousError) -> Self {
        Self::Source(error)
    }
}
fn bad(what: &'static str) -> RadiationError {
    RadiationError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), RadiationError> {
    if gate.is_requested() {
        Err(RadiationError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

type TriangleKey = [[u64; 3]; 3];
fn triangle_key(points: [[f64; 3]; 3]) -> TriangleKey {
    let mut key = points.map(|p| p.map(|x| if x == 0.0 { 0 } else { x.to_bits() }));
    key.sort_unstable();
    key
}
fn face_key(mut face: [u32; 3]) -> [u32; 3] {
    face.sort_unstable();
    face
}

/// Source-bound complete FE pressure field and its physical exterior openings.
pub struct RadiatingInterior<'surface> {
    interior: InteriorMultiport,
    surface: &'surface SpherePanels,
    gas: GasState,
    radiating_modes: Vec<usize>,
    driven_modes: Vec<usize>,
    apertures: Vec<Aperture>,
}

/// One original loaded FE factor, one full BEM row and one small port closure.
/// Reuse for arbitrary coherent source flows at this same frequency.
pub struct PreparedRadiatingInterior<'model, 'surface> {
    model: &'model RadiatingInterior<'surface>,
    interior: PreparedInteriorScattering<'model>,
    source: ApertureImpedance,
    reaction_impedance: Vec<C64>,
    roots: Vec<f64>,
    closure: LuComplex,
    frequency_hz: f64,
}

/// Original complete pressure field and independently evaluated radiation work.
#[derive(Debug, Clone)]
pub struct RadiationResponse {
    /// Amplitudes contain EVERY original `HelmholtzFe` pressure DOF in its
    /// original order. Port pressures/flows retain the declared mode order.
    pub interior: InteriorPortResponse,
    /// Outward modal volume-flow coordinates, in `radiating_modes()` order.
    /// These need not equal net volume flow for a signed mode.
    pub outward_flow_m3_s: Vec<C64>,
    /// Exterior pressure reaction Z U, in the FE exp(+i omega t) convention.
    pub exterior_pressure_pa: Vec<C64>,
    /// Cycle-mean work supplied at the explicitly driven interior ports [W].
    pub source_power_w: f64,
    /// Cycle-mean exterior work Re(U^H Z U)/2 [W], without clipping.
    pub radiation_power_w: f64,
    /// Native wall operator's separate viscous/thermal dissipation, evaluated
    /// on the original complete field; not inferred from a power residual.
    pub boundary_layer_power: BoundaryLayerPower,
    /// Maximum original driven-flow/opening-pressure constraint residual,
    /// scaled by the largest pressure/flow power-coordinate magnitude.
    pub port_relative_residual: f64,
    /// Total driven-port, boundary and volume source work minus exterior work
    /// and native wall losses [W].
    /// This measured solver residual is never relabeled as physical heat.
    pub power_balance_residual_w: f64,
}

impl RadiationResponse {
    /// All signed physical excitation work, counted once on the final field.
    /// This is an energy-balance observation, not an additional heat source.
    pub fn total_source_power_w(&self) -> f64 {
        self.source_power_w
            + self.interior.boundary_source_power_w
            + self.interior.volume_source_power_w
    }
}

/// Shared original geometry/shape admission for harmonic and causal owners.
#[allow(clippy::too_many_arguments)]
fn bind_openings(
    fe: &HelmholtzFe<'_>,
    ports: &[BoundaryMode],
    radiating_modes: &[usize],
    surface: &SpherePanels,
    gas: &GasState,
    budget: RadiationBudget,
    gate: &CancelGate,
) -> Result<(Vec<usize>, Vec<Aperture>), RadiationError> {
    poll(gate)?;
    boundary_layer_transport(gas)?;
    let medium = fe.medium();
    if medium.density_kg_m3.to_bits() != gas.density.to_bits()
        || medium.sound_speed_m_s.to_bits() != gas.sound_speed.to_bits()
    {
        return Err(bad(
            "original interior and exterior must use the same coherent source gas",
        ));
    }
    let panels = surface.areas().len();
    let count = radiating_modes.len();
    let terms = panels as u128 * count as u128 * count as u128;
    if ports.is_empty()
        || ports.len() > MAX_PORT_MODES
        || count == 0
        || count > ports.len()
        || panels == 0
        || panels > budget.maximum_panels
        || panels
            .checked_mul(count)
            .is_none_or(|n| n > budget.maximum_trace_entries)
        || terms > u128::from(budget.maximum_projection_terms)
    {
        return Err(bad(
            "complete opening trace bank exceeds its declared limits",
        ));
    }
    let mut selected = vec![false; ports.len()];
    for &mode in radiating_modes {
        if mode >= ports.len() || selected[mode] {
            return Err(bad("unique original radiating mode indices required"));
        }
        selected[mode] = true;
    }
    let driven_modes: Vec<_> = (0..ports.len()).filter(|&i| !selected[i]).collect();
    let radiation_faces: BTreeSet<_> = radiating_modes
        .iter()
        .flat_map(|&i| ports[i].traces.iter().map(|(face, _)| face_key(*face)))
        .collect();
    if driven_modes.iter().any(|&i| {
        ports[i]
            .traces
            .iter()
            .any(|(face, _)| radiation_faces.contains(&face_key(*face)))
    }) {
        return Err(bad(
            "driven interior feed and radiation modes share an opening face",
        ));
    }

    if fe.ndof() > budget.maximum_pressure_dofs {
        return Err(bad(
            "original pressure field exceeds its declared DOF limit",
        ));
    }
    // Validate every original trace before indexing or creating a panel bank.
    // Harmonic and causal construction use the same physical opening binding.
    for port in ports {
        if port.traces.is_empty() {
            return Err(bad("each original opening mode needs boundary traces"));
        }
        let mut faces = BTreeSet::new();
        for &(face, values) in &port.traces {
            poll(gate)?;
            if !faces.insert(face_key(face)) || values.iter().any(|x| !x.is_finite()) {
                return Err(bad("opening traces need unique faces and finite values"));
            }
            fe.boundary_geometry(face)?;
        }
    }
    GeometryPolicy::new(surface).map_err(RadiationError::Geometry)?;
    poll(gate)?;
    let triangles = surface
        .triangles()
        .ok_or_else(|| bad("exterior requires original retained triangles"))?;
    let mut panel_by_triangle = BTreeMap::new();
    for (i, &triangle) in triangles.iter().enumerate() {
        poll(gate)?;
        if panel_by_triangle
            .insert(triangle_key(triangle), i)
            .is_some()
        {
            return Err(bad("duplicate physical exterior triangle"));
        }
    }
    let mut apertures = Vec::with_capacity(count);
    for &mode in radiating_modes {
        poll(gate)?;
        let first = ports[mode].traces[0].0;
        let (_, normal, _) = fe.boundary_geometry(first)?;
        let reference_point = fe.positions()[first[0] as usize];
        let mut shape = vec![0.0; panels];
        for &(face, values) in &ports[mode].traces {
            poll(gate)?;
            let points = face.map(|v| fe.positions()[v as usize]);
            let panel = *panel_by_triangle
                .get(&triangle_key(points))
                .ok_or_else(|| bad("opening FE triangle is absent from original exterior"))?;
            let (_, face_normal, _) = fe.boundary_geometry(face)?;
            if face_normal
                .iter()
                .zip(normal)
                .any(|(a, b)| (a - b).abs() > 1e-12)
                || surface.normals()[panel]
                    .iter()
                    .zip(normal)
                    .any(|(a, b)| (a - b).abs() > 1e-12)
            {
                return Err(bad(
                    "opening must share its planar outward normal with the exterior",
                ));
            }
            for point in points {
                let offset: f64 = (0..3)
                    .map(|c| normal[c] * (point[c] - reference_point[c]))
                    .sum();
                let scale = point
                    .iter()
                    .chain(&reference_point)
                    .map(|x| x.abs())
                    .fold(f64::MIN_POSITIVE, f64::max);
                if !offset.is_finite() || offset.abs() > 64.0 * f64::EPSILON * scale {
                    return Err(bad("each radiating mode must lie on one physical plane"));
                }
            }
            // Exact area mean of the original real P1 trace, retaining its
            // sign. Never divide by the mode's possibly zero net volume.
            shape[panel] = values[0] / 3.0 + values[1] / 3.0 + values[2] / 3.0;
        }
        if shape.iter().any(|v| !v.is_finite()) || shape.iter().all(|&v| v == 0.0) {
            return Err(bad(
                "exterior P0 mesh does not resolve a supplied signed opening mode",
            ));
        }
        apertures.push(Aperture {
            shape,
            reference_point,
            normal,
        });
    }

    poll(gate)?;
    Ok((driven_modes, apertures))
}

impl<'surface> RadiatingInterior<'surface> {
    /// Bind original pressure coordinates, source gas and identical opening
    /// triangles. All modes supported on a radiating face must be radiating;
    /// interior feeds cannot share those faces. Different radiation modes may
    /// overlap, as required for transverse opening modes. `walls`, when present,
    /// uses the existing source-gas boundary-layer operator on disjoint faces.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fe: &HelmholtzFe<'_>,
        ports: &[BoundaryMode],
        radiating_modes: &[usize],
        surface: &'surface SpherePanels,
        gas: &GasState,
        walls: Option<(&[[u32; 3]], BoundaryLayerGeometry)>,
        budget: RadiationBudget,
        gate: &CancelGate,
    ) -> Result<Self, RadiationError> {
        let (driven_modes, apertures) =
            bind_openings(fe, ports, radiating_modes, surface, gas, budget, gate)?;
        let transport = boundary_layer_transport(gas)?;
        let cancelled = || gate.is_requested();
        let interior = match walls {
            Some((faces, geometry)) => InteriorMultiport::new_full_with_boundary_layers(
                fe,
                ports,
                budget.maximum_pressure_dofs,
                faces,
                transport,
                geometry,
                &cancelled,
            )?,
            None => {
                InteriorMultiport::new_full(fe, ports, budget.maximum_pressure_dofs, &cancelled)?
            }
        };
        poll(gate)?;
        Ok(Self {
            interior,
            surface,
            gas: *gas,
            radiating_modes: radiating_modes.to_vec(),
            driven_modes,
            apertures,
        })
    }

    /// Complete original pressure operator; immutable so the exterior cannot
    /// become detached from the source geometry or pressure work columns.
    pub fn interior(&self) -> &InteriorMultiport {
        &self.interior
    }
    pub fn radiating_modes(&self) -> &[usize] {
        &self.radiating_modes
    }
    /// Remaining original mode indices, in original order. A solve requires
    /// one explicit complex inward flow for each; zero means a sealed mode.
    pub fn driven_modes(&self) -> &[usize] {
        &self.driven_modes
    }
    /// Actual dimensionless signed P0 velocity traces submitted to the BEM.
    /// Physical outward velocities are `shape * U / geometric_support_area`.
    pub fn apertures(&self) -> &[Aperture] {
        &self.apertures
    }

    /// Assemble the complete native radiation matrix, convert its convention
    /// once, and close it against the original Robin-loaded pressure solve.
    /// The selected BEM route's resolution/band refusals are propagated. A
    /// measured passivity defect beyond roundoff refuses without matrix repair.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        route: &ApertureRoute,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingInterior<'_, 'surface>, RadiationError> {
        poll(gate)?;
        let cancelled = || gate.is_requested();
        let interior =
            self.interior
                .prepare_scattering(frequency_hz, reference_resistance, &cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let source = aperture_impedance(
            self.surface,
            &self.apertures,
            &[omega],
            Medium {
                density: self.gas.density,
                sound_speed: self.gas.sound_speed,
            },
            route,
            gate,
        )
        .map_err(RadiationError::Exterior)?;
        let row = &source.rows[0];
        if !row.passivity_margin.is_finite() || row.passivity_margin < -1e-12 {
            return Err(bad(
                "original exterior radiation matrix has a resolved passivity deficit",
            ));
        }
        let nr = self.radiating_modes.len();
        let areas = self.interior.port_areas_m2();
        // Preserve native dimensionless shape and N s/m artifact units. The
        // work-conjugate FE pressure/volume coordinates require TWO geometric
        // support areas; zero signed net-flow gains play no role here.
        let reaction_impedance: Vec<_> = row
            .impedance
            .iter()
            .enumerate()
            .map(|(entry, z)| {
                let r = self.radiating_modes[entry / nr];
                let s = self.radiating_modes[entry % nr];
                z.conj().scale(areas[r].recip()).scale(areas[s].recip())
            })
            .collect();
        if reaction_impedance.iter().any(|&z| !finite(z)) {
            return Err(bad("geometric opening-area conversion is unrepresentable"));
        }
        let scattering = interior.scattering(&cancelled)?;
        let roots: Vec<_> = reference_resistance.iter().map(|&r| det::sqrt(r)).collect();
        let n = self.interior.port_mode_count();
        let mut matrix = vec![C64::ZERO; n * n];
        for &i in &self.driven_modes {
            poll(gate)?;
            for j in 0..n {
                matrix[i * n + j] =
                    C64::from_re(if i == j { 1.0 } else { 0.0 }) - scattering[i * n + j];
            }
        }
        for (r, &i) in self.radiating_modes.iter().enumerate() {
            poll(gate)?;
            for j in 0..n {
                let mut value =
                    C64::from_re(if i == j { 1.0 } else { 0.0 }) + scattering[i * n + j];
                for (s, &k) in self.radiating_modes.iter().enumerate() {
                    let flow = C64::from_re(if k == j { 1.0 } else { 0.0 }) - scattering[k * n + j];
                    value = value
                        + reaction_impedance[r * nr + s].scale(roots[i].recip() * roots[k].recip())
                            * flow;
                }
                matrix[i * n + j] = value;
            }
        }
        if matrix.iter().any(|&z| !finite(z)) {
            return Err(bad(
                "radiation closure is unrepresentable in the supplied power references",
            ));
        }
        let closure =
            lu_complex_with_cancel(&matrix, n, &cancelled).map_err(RadiationError::Closure)?;
        poll(gate)?;
        Ok(PreparedRadiatingInterior {
            model: self,
            interior,
            source,
            reaction_impedance,
            roots,
            closure,
            frequency_hz,
        })
    }
}

impl PreparedRadiatingInterior<'_, '_> {
    /// Unmodified native BEM artifact, including signed net-flow diagnostics,
    /// every matrix entry, route, source declaration and passivity/reciprocity.
    pub fn source(&self) -> &ApertureImpedance {
        &self.source
    }
    /// Full modal pressure/flow impedance in the FE exp(+i omega t) convention.
    pub fn reaction_impedance(&self) -> &[C64] {
        &self.reaction_impedance
    }

    /// Solve simultaneous prescribed inward flows in `driven_modes()` order.
    /// Both the complete FE field and full exterior reaction use the cached
    /// original factors; source changes require no radiation reassembly.
    pub fn response(
        &self,
        driven_inward_flow_m3_s: &[C64],
        gate: &CancelGate,
    ) -> Result<RadiationResponse, RadiationError> {
        self.response_with_tolerance(
            driven_inward_flow_m3_s,
            DEFAULT_PORT_RELATIVE_TOLERANCE,
            gate,
        )
    }

    /// Same original solve with an explicit positive relative port residual
    /// ceiling below one. Refuses a failed physical interface/source constraint
    /// rather than hiding it in recovered flows or inferred dissipation.
    pub fn response_with_tolerance(
        &self,
        driven_inward_flow_m3_s: &[C64],
        maximum_port_relative_residual: f64,
        gate: &CancelGate,
    ) -> Result<RadiationResponse, RadiationError> {
        self.response_with_sources_and_tolerance(
            driven_inward_flow_m3_s,
            &[],
            &[],
            maximum_port_relative_residual,
            gate,
        )
    }

    /// Close original distributed injection and inward boundary motion against
    /// the same full self/mutual exterior load as the prescribed port flows.
    ///
    /// Prepare each source through this model's `interior()` and its original
    /// `HelmholtzFe`. Volume amplitudes are peak m^3/s; quadratic boundary
    /// amplitudes are peak inward m/s, all on the FE exp(+i omega t) clock.
    /// Native exact-owner, signed-shape, support and source-count checks apply.
    /// Sources are admitted even when their requested amplitude is zero.
    ///
    /// All source columns combine before the retained loaded solve. One affine
    /// outgoing source response enters the existing small port closure, then
    /// the same loaded factor reconstructs the complete final pressure field.
    /// There is no per-source factorization or exterior reassembly. Signed
    /// source work is evaluated on that final field, including interference.
    ///
    /// Boundary motion represents an interior actuator/feed on non-opening
    /// faces. The exterior remains rigid outside its declared openings; this
    /// does not model direct exterior radiation from a moving structural wall.
    pub fn response_with_sources(
        &self,
        driven_inward_flow_m3_s: &[C64],
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        gate: &CancelGate,
    ) -> Result<RadiationResponse, RadiationError> {
        self.response_with_sources_and_tolerance(
            driven_inward_flow_m3_s,
            boundary,
            volume,
            DEFAULT_PORT_RELATIVE_TOLERANCE,
            gate,
        )
    }

    /// Source-driven response with the original explicit physical port
    /// residual ceiling. Refusal/cancellation publishes no partial response;
    /// the source projections and both original factors remain reusable.
    pub fn response_with_sources_and_tolerance(
        &self,
        driven_inward_flow_m3_s: &[C64],
        boundary: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        maximum_port_relative_residual: f64,
        gate: &CancelGate,
    ) -> Result<RadiationResponse, RadiationError> {
        poll(gate)?;
        if !maximum_port_relative_residual.is_finite()
            || maximum_port_relative_residual <= 0.0
            || maximum_port_relative_residual >= 1.0
        {
            return Err(bad(
                "finite positive physical port residual tolerance below one required",
            ));
        }
        if driven_inward_flow_m3_s.len() != self.model.driven_modes.len()
            || driven_inward_flow_m3_s.iter().any(|&q| !finite(q))
        {
            return Err(bad(
                "one finite complex flow per original driven mode required",
            ));
        }
        let cancelled = || gate.is_requested();
        let mut incident = vec![C64::ZERO; self.model.interior.port_mode_count()];
        for (&i, &q) in self.model.driven_modes.iter().zip(driven_inward_flow_m3_s) {
            incident[i] = q.scale(self.roots[i]);
        }
        let has_sources = !boundary.is_empty() || !volume.is_empty();
        if has_sources {
            let outgoing = self
                .interior
                .source_outgoing(boundary, volume, &cancelled)?;
            // b = S a + e, P = sqrt(R) (a+b), Q = (a-b)/sqrt(R).
            // Driven rows require (I-S)a = sqrt(R)Q + e. Radiating rows
            // require [(I+S)+Zbar(I-S)]a = -e + Zbar e, Zbar=R^-1/2 Z R^-1/2.
            // These are affine loads on the original factor, not new impedance.
            for &i in &self.model.driven_modes {
                poll(gate)?;
                incident[i] = incident[i] + outgoing[i];
            }
            let nr = self.model.radiating_modes.len();
            for (r, &i) in self.model.radiating_modes.iter().enumerate() {
                poll(gate)?;
                let mut value = -outgoing[i];
                for (s, &k) in self.model.radiating_modes.iter().enumerate() {
                    value = value
                        + self.reaction_impedance[r * nr + s]
                            .scale(self.roots[i].recip() * self.roots[k].recip())
                            * outgoing[k];
                }
                incident[i] = value;
            }
        }
        if incident.iter().any(|&z| !finite(z)) {
            return Err(bad("source power coordinates overflow"));
        }
        self.closure
            .solve_with_cancel(&mut incident, &cancelled)
            .map_err(RadiationError::Closure)?;
        let interior = if has_sources {
            self.interior
                .response_with_sources(&incident, boundary, volume, &cancelled)?
        } else {
            self.interior.response(&incident, &cancelled)?
        };
        let outward_flow_m3_s: Vec<_> = self
            .model
            .radiating_modes
            .iter()
            .map(|&i| interior.inward_flow_m3_s[i].scale(-1.0))
            .collect();
        let nr = outward_flow_m3_s.len();
        let exterior_pressure_pa: Vec<_> = self
            .reaction_impedance
            .chunks(nr)
            .map(|row| {
                row.iter()
                    .zip(&outward_flow_m3_s)
                    .fold(C64::ZERO, |sum, (&z, &q)| sum + z * q)
            })
            .collect();
        let radiation_power_w = 0.5
            * exterior_pressure_pa
                .iter()
                .zip(&outward_flow_m3_s)
                .map(|(&p, &q)| (q.conj() * p).re)
                .sum::<f64>();
        let source_power_w = 0.5
            * self
                .model
                .driven_modes
                .iter()
                .zip(driven_inward_flow_m3_s)
                .map(|(&i, &prescribed)| (interior.pressure_pa[i].conj() * prescribed).re)
                .sum::<f64>();
        if exterior_pressure_pa.iter().any(|&z| !finite(z))
            || !radiation_power_w.is_finite()
            || !source_power_w.is_finite()
        {
            return Err(bad("radiation reaction or physical work overflow"));
        }
        let mut scale = 0.0_f64;
        let mut residual = 0.0_f64;
        for (i, (&p, &q)) in interior
            .pressure_pa
            .iter()
            .zip(&interior.inward_flow_m3_s)
            .enumerate()
        {
            let pressure = p.abs() / self.roots[i];
            let flow = q.abs() * self.roots[i];
            if !pressure.is_finite() || !flow.is_finite() {
                return Err(bad("physical port residual scaling overflow"));
            }
            scale = scale.max(pressure).max(flow);
        }
        for (&i, &prescribed) in self.model.driven_modes.iter().zip(driven_inward_flow_m3_s) {
            let error = (interior.inward_flow_m3_s[i] - prescribed).abs() * self.roots[i];
            let input = prescribed.abs() * self.roots[i];
            if !error.is_finite() || !input.is_finite() {
                return Err(bad("physical driven-flow residual overflow"));
            }
            residual = residual.max(error);
            scale = scale.max(input);
        }
        for (&i, &pressure) in self.model.radiating_modes.iter().zip(&exterior_pressure_pa) {
            let error = (interior.pressure_pa[i] - pressure).abs() / self.roots[i];
            let exterior = pressure.abs() / self.roots[i];
            if !error.is_finite() || !exterior.is_finite() {
                return Err(bad("physical opening-pressure residual overflow"));
            }
            residual = residual.max(error);
            scale = scale.max(exterior);
        }
        let port_relative_residual = if scale == 0.0 { 0.0 } else { residual / scale };
        if port_relative_residual > maximum_port_relative_residual {
            return Err(bad(
                "original interior/exterior port constraints exceed residual tolerance",
            ));
        }
        let boundary_layer_power = self
            .model
            .interior
            .boundary_layer_power(self.frequency_hz, &interior.amplitudes)?;
        let power_balance_residual_w =
            source_power_w + interior.boundary_source_power_w + interior.volume_source_power_w
                - radiation_power_w
                - boundary_layer_power.viscous_w
                - boundary_layer_power.thermal_w;
        if !power_balance_residual_w.is_finite() {
            return Err(bad("physical power-balance residual overflow"));
        }
        poll(gate)?;
        Ok(RadiationResponse {
            interior,
            outward_flow_m3_s,
            exterior_pressure_pa,
            source_power_w,
            radiation_power_w,
            boundary_layer_power,
            port_relative_residual,
            power_balance_residual_w,
        })
    }
}

#[cfg(test)]
mod tests;
