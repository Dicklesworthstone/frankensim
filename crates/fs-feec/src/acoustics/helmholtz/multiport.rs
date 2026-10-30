//! Spatial boundary modes of the original three-dimensional Helmholtz problem.
//!
//! A real, piecewise-linear shape psi on a port of area A prescribes inward
//! velocity `v = psi Q/A`. Its conjugate pressure is `P = integral(psi p)/A`,
//! so mean input power is exactly `Re(sum(conj(P) Q))/2`, including overlapping
//! transverse shapes. Q is a modal volume-flow coordinate, not the net flow of
//! a zero-mean mode. No projection onto a face average discards those modes.
//!
//! The original volume/trace quadrature, Galerkin basis, static enrichment and
//! thermoviscous wall operator are reused. One complex LU factorization serves
//! every excitation at a frequency. Traces are a convergent P1 representation
//! of supplied cross-section modes; geometry, trace, volume and basis refinement
//! remain independent obligations. No analytic duct mode or cutoff is invented.
//!
//! Time convention is exp(+i omega t), as in the parent FE owner. Impedance uses
//! inward flows. Scattering uses positive REAL reference resistances and power
//! coordinates, including for reactive near-field traces. These are not the
//! traveling-wave amplitudes of evanescent modes; a full evanescent scattering
//! matrix must not be tested for Euclidean contractivity under that convention
//! (see Svendsen, Haakestad & Skaar, Phys. Rev. A 87, 013838, 2013).

use super::boundary_layer::thermal_wall::ThermalWall;
use super::boundary_layer::{BoundaryLayerGeometry, BoundaryLayerPower, BoundaryLayerTransport};
use super::{AcousticFeError, HelmholtzFe, HelmholtzModes, ReducedHelmholtz};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;
use fs_math::det;

/// Explicit cap on retained boundary modes. Exceeding it refuses, never truncates.
pub const MAX_PORT_MODES: usize = 64;

/// One real normal-velocity shape, in caller-specified vertex order on each face.
/// Faces within a mode must be unique. Different modes may share faces.
/// Values are dimensionless; the caller owns modal normalization and orientation.
#[derive(Debug, Clone)]
pub struct BoundaryMode {
    /// `(boundary face, psi at its three vertices)`.
    pub traces: Vec<([u32; 3], [f64; 3])>,
}

impl BoundaryMode {
    /// The uniform shape psi=1, recovering ordinary volume flow and mean pressure.
    #[must_use]
    pub fn uniform(faces: &[[u32; 3]]) -> Self {
        Self {
            traces: faces.iter().map(|&face| (face, [1.0; 3])).collect(),
        }
    }
}

/// Prepared spatial port projection on the existing Helmholtz reduction.
pub struct InteriorMultiport {
    pub(super) reduced: ReducedHelmholtz,
    pub(crate) physical_modes: Vec<Vec<PhysicalTrace>>,
}

/// Original source traces retained for an exactly identified hybrid interface.
/// Coordinates use normalized-zero bits, so vertex numbering and winding may
/// differ across the interface without accepting a substituted trace shape.
#[derive(Clone)]
pub(crate) struct PhysicalTrace {
    pub vertices: [[u64; 3]; 3],
    pub values: [f64; 3],
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
}

pub(crate) fn physical_trace(
    fe: &HelmholtzFe<'_>,
    mode: &BoundaryMode,
) -> Result<Vec<PhysicalTrace>, AcousticFeError> {
    mode.traces
        .iter()
        .map(|&(face, values)| {
            let (_, outward_normal, area_m2) = fe.boundary_geometry(face)?;
            let mut vertices = std::array::from_fn::<_, 3, _>(|i| {
                let point = fe.positions()[face[i] as usize].map(|x| {
                    if x == 0.0 {
                        0.0_f64.to_bits()
                    } else {
                        x.to_bits()
                    }
                });
                (point, values[i])
            });
            vertices.sort_by_key(|&(point, _)| point);
            Ok(PhysicalTrace {
                vertices: vertices.map(|(point, _)| point),
                values: vertices.map(|(_, value)| value),
                outward_normal,
                area_m2,
            })
        })
        .collect()
}

/// One factorization, reusable for arbitrary simultaneous modal excitations.
/// Borrowing the model prevents changing its material/wall configuration midway
/// through a frequency solve. No partially factored result is returned.
pub struct PreparedInteriorMultiport<'a> {
    model: &'a InteriorMultiport,
    frequency_hz: f64,
    omega_rho: f64,
    factor: LuComplex,
}

/// Original pressure FE operator with real-reference Robin loads, factored
/// once for arbitrary simultaneous incident power waves. This remains usable
/// at a pole of the unloaded impedance whenever the loaded operator is regular.
pub struct PreparedInteriorScattering<'a> {
    model: &'a InteriorMultiport,
    frequency_hz: f64,
    omega_rho: f64,
    reference: Vec<f64>,
    roots: Vec<f64>,
    factor: LuComplex,
    bulk_loss_factor: f64,
}

/// Complex pressure and flow coordinates, with their original FE amplitudes.
#[derive(Debug, Clone)]
pub struct InteriorPortResponse {
    /// Conjugate pressure `integral(psi p)/A` [Pa], in declared mode order.
    pub pressure_pa: Vec<C64>,
    /// Prescribed inward modal flow coordinates [m^3/s].
    pub inward_flow_m3_s: Vec<C64>,
    /// Original reduced pressure amplitudes `(re, im)`; same basis as the model.
    pub amplitudes: Vec<[f64; 2]>,
}

impl InteriorPortResponse {
    /// Mean total acoustic input work [W], including modal interference through
    /// the solved pressure. There is no independent area or impedance scaling.
    #[must_use]
    pub fn input_power_w(&self) -> f64 {
        0.5 * self
            .pressure_pa
            .iter()
            .zip(&self.inward_flow_m3_s)
            .map(|(p, q)| (p.conj() * *q).re)
            .sum::<f64>()
    }
}

fn invalid(what: &'static str) -> AcousticFeError {
    AcousticFeError::InvalidInput { what }
}

fn poll(cancelled: &dyn Fn() -> bool) -> Result<(), AcousticFeError> {
    if cancelled() {
        Err(invalid("interior multiport cancelled"))
    } else {
        Ok(())
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn linear_error(error: LuComplexError, frequency_hz: f64) -> AcousticFeError {
    match error {
        LuComplexError::Cancelled => invalid("interior multiport cancelled"),
        LuComplexError::Singular { .. } => AcousticFeError::SingularAtFrequency { frequency_hz },
    }
}

impl InteriorMultiport {
    /// Project supplied spatial modes onto the original rigid-wall FE problem.
    /// Static enrichment, when requested, uses each ACTUAL mode load.
    pub fn new(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[BoundaryMode],
        static_vectors: bool,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, AcousticFeError> {
        Self::build(fe, modes, ports, static_vectors, None, cancelled)
    }

    /// Same projection with the existing spatial viscous/isothermal wall layers.
    /// Port modes may overlap one another, but none may overlap a wall face.
    pub fn new_with_boundary_layers(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[BoundaryMode],
        static_vectors: bool,
        wall_faces: &[[u32; 3]],
        transport: BoundaryLayerTransport,
        geometry: BoundaryLayerGeometry,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, AcousticFeError> {
        Self::build(
            fe,
            modes,
            ports,
            static_vectors,
            Some((wall_faces, transport, geometry)),
            cancelled,
        )
    }

    fn build(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[BoundaryMode],
        static_vectors: bool,
        wall: Option<(&[[u32; 3]], BoundaryLayerTransport, BoundaryLayerGeometry)>,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, AcousticFeError> {
        poll(cancelled)?;
        if ports.is_empty() || ports.len() > MAX_PORT_MODES {
            return Err(invalid(
                "interior port mode count outside 1..=MAX_PORT_MODES",
            ));
        }
        if modes.shapes.is_empty()
            || modes
                .shapes
                .iter()
                .any(|p| p.len() != fe.ndof() || p.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid(
                "interior port reduction needs a finite nonempty FE basis",
            ));
        }
        let mut faces = Vec::with_capacity(ports.len());
        let mut loads = Vec::with_capacity(ports.len());
        for port in ports {
            poll(cancelled)?;
            // The original trace owner checks boundary ownership, duplicates,
            // winding and all nodal values before any projection is published.
            let mut load = fe.boundary_linear_load(&port.traces)?;
            let selected: Vec<_> = port.traces.iter().map(|(face, _)| *face).collect();
            let (_, area) = fe.boundary_load(&selected)?;
            if !area.is_finite() || area <= 0.0 || load.iter().all(|&v| v == 0.0) {
                return Err(invalid(
                    "interior port has no area or unresolved zero trace",
                ));
            }
            for value in &mut load {
                *value /= area;
            }
            faces.push(selected);
            loads.push(load);
        }
        let face_sets: Vec<_> = faces.iter().map(Vec::as_slice).collect();
        let reduced = ReducedHelmholtz::build_with_port_loads(
            fe,
            modes,
            &face_sets,
            static_vectors,
            wall,
            Some(&loads),
        )?;
        if reduced.b_r.iter().any(|b| b.iter().all(|&v| v == 0.0)) {
            return Err(invalid("interior reduction discards a supplied port mode"));
        }
        poll(cancelled)?;
        let physical_modes = ports
            .iter()
            .map(|port| physical_trace(fe, port))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            reduced,
            physical_modes,
        })
    }

    /// Bind the original conductive/adiabatic thermal wall to this same port
    /// projection. Wall matrices and spatial port loads are not rebuilt.
    pub fn with_thermal_wall(mut self, wall: ThermalWall) -> Result<Self, AcousticFeError> {
        self.reduced = self.reduced.with_thermal_wall(wall)?;
        Ok(self)
    }

    /// Number of spatial mode coordinates, including zero-mean near-field modes.
    #[must_use]
    pub fn port_mode_count(&self) -> usize {
        self.reduced.b_r.len()
    }

    /// Actual reduced FE dimension, including retained static enrichment.
    #[must_use]
    pub const fn dim(&self) -> usize {
        self.reduced.n
    }

    /// Physical area [m^2] supporting each mode, without automatic renormalization.
    #[must_use]
    pub fn port_areas_m2(&self) -> &[f64] {
        &self.reduced.area
    }

    /// Uniform medium used by the original pressure operator.
    #[must_use]
    pub const fn medium(&self) -> super::Medium {
        self.reduced.medium
    }

    /// Original gas coefficients used by the optional pressure-wall operator.
    /// Source adapters can verify transport as well as volume rho/c without
    /// replacing or reconstructing the original spatial wall matrices.
    #[must_use]
    pub fn boundary_layer_transport(&self) -> Option<BoundaryLayerTransport> {
        self.reduced
            .boundary_layer
            .as_ref()
            .map(|wall| wall.transport)
    }

    pub(super) fn operator(
        &self,
        frequency_hz: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(Vec<C64>, f64), AcousticFeError> {
        poll(cancelled)?;
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(invalid("interior multiport frequency"));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        let omega_rho = omega * self.reduced.medium.density_kg_m3;
        let k = omega / self.reduced.medium.sound_speed_m_s;
        if !omega_rho.is_finite() || !(k * k).is_finite() {
            return Err(invalid("interior multiport frequency overflow"));
        }
        let mut a: Vec<_> = self
            .reduced
            .k_r
            .iter()
            .zip(&self.reduced.m_r)
            .map(|(stiffness, mass)| [stiffness - k * k * mass, 0.0])
            .collect();
        if let Some(wall) = &self.reduced.boundary_layer {
            wall.add_operator(self.reduced.medium, frequency_hz, &mut a)?;
        }
        let a: Vec<_> = a.into_iter().map(|z| C64::new(z[0], z[1])).collect();
        if a.iter().any(|&z| !finite(z)) {
            return Err(invalid("interior multiport operator overflow"));
        }
        Ok((a, omega_rho))
    }

    /// Factor once for imposed modal flows. A lossless closed-cavity pole may
    /// be singular; [`Self::scattering`] instead factors a loaded operator.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedInteriorMultiport<'_>, AcousticFeError> {
        let (a, omega_rho) = self.operator(frequency_hz, cancelled)?;
        self.prepare_operator(a, omega_rho, frequency_hz, cancelled)
    }

    pub(super) fn prepare_operator(
        &self,
        a: Vec<C64>,
        omega_rho: f64,
        frequency_hz: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedInteriorMultiport<'_>, AcousticFeError> {
        let factor = lu_complex_with_cancel(&a, self.dim(), cancelled)
            .map_err(|e| linear_error(e, frequency_hz))?;
        Ok(PreparedInteriorMultiport {
            model: self,
            frequency_hz,
            omega_rho,
            factor,
        })
    }

    /// Row-major reciprocal impedance Z [Pa s/m^3], P=Z Q, using one LU and
    /// one right-hand side per spatial mode. No symmetry repair is applied.
    pub fn impedance(
        &self,
        frequency_hz: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, AcousticFeError> {
        let prepared = self.prepare(frequency_hz, cancelled)?;
        let ports = self.port_mode_count();
        let mut matrix = vec![C64::ZERO; ports * ports];
        let mut flow = vec![C64::ZERO; ports];
        for col in 0..ports {
            flow[col] = C64::ONE;
            let response = prepared.response(&flow, cancelled)?;
            flow[col] = C64::ZERO;
            for (row, value) in response.pressure_pa.into_iter().enumerate() {
                matrix[row * ports + col] = value;
            }
        }
        poll(cancelled)?;
        Ok(matrix)
    }

    /// Row-major power scattering b=S a for positive real R [Pa s/m^3]:
    /// `a=(P+R Q)/(2 sqrt(R))`, `b=(P-R Q)/(2 sqrt(R))`.
    /// Mean input work is `(||a||^2-||b||^2)/2`. Reactive/evanescent traces
    /// remain in the solve; R is a coordinate choice, not their characteristic
    /// impedance. Direct Robin loading avoids first forming a pole in Z.
    pub fn scattering(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, AcousticFeError> {
        self.validate_references(reference_resistance)?;
        let (a, omega_rho) = self.operator(frequency_hz, cancelled)?;
        self.scattering_operator(a, omega_rho, frequency_hz, reference_resistance, cancelled)
    }

    /// Retain the loaded factorization for field recovery and network closure.
    /// Incident coordinates are a=(P+RQ)/(2 sqrt(R)); flow is inward.
    pub fn prepare_scattering(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedInteriorScattering<'_>, AcousticFeError> {
        self.validate_references(reference_resistance)?;
        let (a, omega_rho) = self.operator(frequency_hz, cancelled)?;
        self.prepare_scattering_operator(
            a,
            omega_rho,
            frequency_hz,
            reference_resistance,
            cancelled,
        )
    }

    /// Prepare the original pressure operator with a source-evaluated uniform
    /// bulk modulus K [Pa] in exp(+iwt). Passive K has Re(K)>0 and Im(K)>=0.
    /// The mass coefficient is rho omega^2/K; momentum remains inviscid.
    /// This is a frequency sample, not a broadband causal material realization.
    ///
    /// An existing boundary-layer overlay refuses: its equilibrium thermal
    /// term does not represent the coupled molecular/thermal response with
    /// complex K. The joint pressure/temperature molecular owner handles that
    /// separate case. No additional dense matrix or empirical attenuation is
    /// introduced here; the original basis, mass and loaded LU are reused.
    pub fn prepare_scattering_with_bulk_modulus(
        &self,
        frequency_hz: f64,
        bulk_modulus_pa: C64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedInteriorScattering<'_>, AcousticFeError> {
        poll(cancelled)?;
        self.validate_references(reference_resistance)?;
        if self.reduced.boundary_layer.is_some() {
            return Err(invalid(
                "molecular bulk modulus cannot overlay Helmholtz boundary layers",
            ));
        }
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !finite(bulk_modulus_pa)
            || bulk_modulus_pa.re <= 0.0
            || bulk_modulus_pa.im < 0.0
        {
            return Err(invalid(
                "interior molecular frequency or passive bulk modulus",
            ));
        }
        let compressibility = C64::ONE / bulk_modulus_pa;
        let omega = core::f64::consts::TAU * frequency_hz;
        let omega_rho = omega * self.reduced.medium.density_kg_m3;
        let k2 = compressibility.scale(omega_rho * omega);
        let bulk_loss_factor = -0.5 * omega * compressibility.im;
        if !finite(compressibility)
            || compressibility.re <= 0.0
            || !finite(k2)
            || !omega_rho.is_finite()
            || !bulk_loss_factor.is_finite()
        {
            return Err(invalid("interior molecular operator overflow"));
        }
        let n = self.dim();
        let mut a = vec![C64::ZERO; n * n];
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                let i = row * n + col;
                a[i] = C64::from_re(self.reduced.k_r[i]) - k2.scale(self.reduced.m_r[i]);
            }
        }
        let mut prepared = self.prepare_scattering_operator(
            a,
            omega_rho,
            frequency_hz,
            reference_resistance,
            cancelled,
        )?;
        prepared.bulk_loss_factor = bulk_loss_factor;
        Ok(prepared)
    }

    pub(super) fn scattering_operator(
        &self,
        a: Vec<C64>,
        omega_rho: f64,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, AcousticFeError> {
        self.prepare_scattering_operator(
            a,
            omega_rho,
            frequency_hz,
            reference_resistance,
            cancelled,
        )?
        .scattering(cancelled)
    }

    fn prepare_scattering_operator(
        &self,
        mut a: Vec<C64>,
        omega_rho: f64,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedInteriorScattering<'_>, AcousticFeError> {
        self.validate_references(reference_resistance)?;
        let n = self.dim();
        for (b, &resistance) in self.reduced.b_r.iter().zip(reference_resistance) {
            poll(cancelled)?;
            let scale = omega_rho / resistance;
            for i in 0..n {
                for j in 0..n {
                    a[i * n + j].im += scale * b[i] * b[j];
                }
            }
        }
        if a.iter().any(|&z| !finite(z)) {
            return Err(invalid("interior power-wave load overflow"));
        }
        let factor =
            lu_complex_with_cancel(&a, n, cancelled).map_err(|e| linear_error(e, frequency_hz))?;
        let roots: Vec<_> = reference_resistance.iter().map(|r| det::sqrt(*r)).collect();
        poll(cancelled)?;
        Ok(PreparedInteriorScattering {
            model: self,
            frequency_hz,
            omega_rho,
            reference: reference_resistance.to_vec(),
            roots,
            factor,
            bulk_loss_factor: 0.0,
        })
    }

    pub(super) fn validate_references(
        &self,
        reference_resistance: &[f64],
    ) -> Result<(), AcousticFeError> {
        let ports = self.port_mode_count();
        if reference_resistance.len() != ports
            || reference_resistance
                .iter()
                .any(|r| !r.is_finite() || *r <= 0.0)
        {
            return Err(invalid("interior power-wave reference resistance"));
        }
        Ok(())
    }

    /// Use the original wall work owner on amplitudes returned by this model.
    pub fn boundary_layer_power(
        &self,
        frequency_hz: f64,
        amplitudes: &[[f64; 2]],
    ) -> Result<BoundaryLayerPower, AcousticFeError> {
        self.reduced.boundary_layer_power(frequency_hz, amplitudes)
    }
}

fn project(load: &[f64], amplitudes: &[C64]) -> C64 {
    load.iter()
        .zip(amplitudes)
        .fold(C64::ZERO, |sum, (b, p)| sum + p.scale(*b))
}

impl PreparedInteriorMultiport<'_> {
    /// Solve arbitrary simultaneous complex modal flows without refactorization.
    pub fn response(
        &self,
        inward_flow_m3_s: &[C64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<InteriorPortResponse, AcousticFeError> {
        poll(cancelled)?;
        if inward_flow_m3_s.len() != self.model.port_mode_count()
            || inward_flow_m3_s.iter().any(|&q| !finite(q))
        {
            return Err(invalid("interior modal flow shape or nonfinite value"));
        }
        let mut rhs = vec![C64::ZERO; self.model.dim()];
        for (b, &flow) in self.model.reduced.b_r.iter().zip(inward_flow_m3_s) {
            let scale = C64::new(-self.omega_rho * flow.im, self.omega_rho * flow.re);
            for (entry, &weight) in rhs.iter_mut().zip(b) {
                *entry = *entry + scale.scale(weight);
            }
        }
        self.factor
            .solve_with_cancel(&mut rhs, cancelled)
            .map_err(|e| linear_error(e, self.frequency_hz))?;
        let pressure_pa: Vec<_> = self
            .model
            .reduced
            .b_r
            .iter()
            .map(|b| project(b, &rhs))
            .collect();
        if rhs.iter().chain(&pressure_pa).any(|&z| !finite(z)) {
            return Err(invalid("interior modal response overflow"));
        }
        poll(cancelled)?;
        Ok(InteriorPortResponse {
            pressure_pa,
            inward_flow_m3_s: inward_flow_m3_s.to_vec(),
            amplitudes: rhs.into_iter().map(|z| [z.re, z.im]).collect(),
        })
    }
}

impl PreparedInteriorScattering<'_> {
    /// Declared positive real resistances [Pa s/m^3], in original mode order.
    #[must_use]
    pub fn reference_resistance(&self) -> &[f64] {
        &self.reference
    }

    /// Independent volume dissipation [W] in the original pressure basis:
    /// -omega Im(1/K) (p* M_r p)/2. M_r is the existing geometric mass
    /// matrix, not a fitted loss or boundary power defect. The ordinary real
    /// medium path returns zero. No loss clipping or extra dense array occurs.
    pub fn bulk_dissipation_w(
        &self,
        amplitudes: &[[f64; 2]],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<f64, AcousticFeError> {
        poll(cancelled)?;
        let n = self.model.dim();
        if amplitudes.len() != n || amplitudes.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("interior bulk work pressure coordinates"));
        }
        if self.bulk_loss_factor == 0.0 {
            return Ok(0.0);
        }
        let mut quadratic = 0.0;
        for row in 0..n {
            poll(cancelled)?;
            let mut mass_re = 0.0;
            let mut mass_im = 0.0;
            for (col, p) in amplitudes.iter().enumerate() {
                let mass = self.model.reduced.m_r[row * n + col];
                mass_re += mass * p[0];
                mass_im += mass * p[1];
            }
            quadratic += amplitudes[row][0] * mass_re + amplitudes[row][1] * mass_im;
        }
        let work = self.bulk_loss_factor * quadratic;
        if !work.is_finite() {
            return Err(invalid("interior bulk work overflow"));
        }
        poll(cancelled)?;
        Ok(work)
    }

    /// Resolve simultaneous incident waves with the retained original LU.
    /// Both pressure and inward flow are recovered; no impedance inversion
    /// or separate unloaded solve is needed to recover the original FE field.
    pub fn response(
        &self,
        incident: &[C64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<InteriorPortResponse, AcousticFeError> {
        poll(cancelled)?;
        if incident.len() != self.reference.len() || incident.iter().any(|&a| !finite(a)) {
            return Err(invalid("interior incident wave shape or nonfinite value"));
        }
        let mut rhs = vec![C64::ZERO; self.model.dim()];
        for (port, &wave) in incident.iter().enumerate() {
            poll(cancelled)?;
            let scale = C64::new(-wave.im, wave.re).scale(2.0 * self.omega_rho / self.roots[port]);
            for (entry, &weight) in rhs.iter_mut().zip(&self.model.reduced.b_r[port]) {
                *entry = *entry + scale.scale(weight);
            }
        }
        self.factor
            .solve_with_cancel(&mut rhs, cancelled)
            .map_err(|e| linear_error(e, self.frequency_hz))?;
        let pressure_pa = self
            .model
            .reduced
            .b_r
            .iter()
            .map(|b| project(b, &rhs))
            .collect::<Vec<_>>();
        let inward_flow_m3_s = pressure_pa
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                (incident[i].scale(2.0 * self.roots[i]) - p).scale(1.0 / self.reference[i])
            })
            .collect::<Vec<_>>();
        if rhs
            .iter()
            .chain(&pressure_pa)
            .chain(&inward_flow_m3_s)
            .any(|&z| !finite(z))
        {
            return Err(invalid("interior loaded response overflow"));
        }
        poll(cancelled)?;
        Ok(InteriorPortResponse {
            pressure_pa,
            inward_flow_m3_s,
            amplitudes: rhs.into_iter().map(|z| [z.re, z.im]).collect(),
        })
    }

    /// Full power scattering from the same loaded factorization.
    pub fn scattering(&self, cancelled: &dyn Fn() -> bool) -> Result<Vec<C64>, AcousticFeError> {
        let m = self.reference.len();
        let mut scattering = vec![C64::ZERO; m * m];
        let mut incident = vec![C64::ZERO; m];
        for col in 0..m {
            incident[col] = C64::ONE;
            let response = self.response(&incident, cancelled)?;
            incident[col] = C64::ZERO;
            for row in 0..m {
                scattering[row * m + col] = response.pressure_pa[row].scale(1.0 / self.roots[row])
                    - if row == col { C64::ONE } else { C64::ZERO };
            }
        }
        poll(cancelled)?;
        Ok(scattering)
    }
}

#[cfg(test)]
mod tests;
