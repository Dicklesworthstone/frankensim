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
        Ok(Self { reduced })
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

    pub(super) fn scattering_operator(
        &self,
        mut a: Vec<C64>,
        omega_rho: f64,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, AcousticFeError> {
        self.validate_references(reference_resistance)?;
        let ports = self.port_mode_count();
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
        let mut scattering = vec![C64::ZERO; ports * ports];
        for col in 0..ports {
            poll(cancelled)?;
            let mut rhs: Vec<_> = self.reduced.b_r[col]
                .iter()
                .map(|b| C64::new(0.0, 2.0 * omega_rho / roots[col] * b))
                .collect();
            factor
                .solve_with_cancel(&mut rhs, cancelled)
                .map_err(|e| linear_error(e, frequency_hz))?;
            for row in 0..ports {
                let pressure = project(&self.reduced.b_r[row], &rhs);
                let reflected = pressure.scale(1.0 / roots[row])
                    - if row == col { C64::ONE } else { C64::ZERO };
                if !finite(reflected) {
                    return Err(invalid("interior scattering overflow"));
                }
                scattering[row * ports + col] = reflected;
            }
        }
        poll(cancelled)?;
        Ok(scattering)
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

#[cfg(test)]
mod tests;
