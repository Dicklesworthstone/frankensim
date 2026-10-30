//! Joint molecular calorics and resolved heat conduction on the original FE.
//!
//! With exp(+i omega t), the material supplies passive cp(omega), Im(cp)<=0.
//! Ideal-gas continuity and the heat equation give, in the original bases,
//! `K p - omega² rho/p0 M p + omega² rho/T0 H^T T = i omega rho B Q` and
//! `[kappa K_T + i omega rho cp M_T + Y_wall C_T] T = i omega H p`.
//! The pressure/temperature equation follows the linearized ideal-gas energy
//! and continuity equations, with the material's dynamic caloric law:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_thermo.09.79.html>.
//! Temperature is eliminated once per frequency; all pressure ports reuse the
//! same heat elimination and original cancellable pressure LU. No independent
//! complex bulk modulus or gas thermal wall term is added a second time.
//!
//! Heat is resolved throughout the original volume, including overlapping
//! thermal layers. Optional viscosity still uses the original THIN tangential
//! wall layer; this is not resolved viscous LNS. Pressure uses the explicitly
//! supplied basis without automatic static enrichment. Temperature retains all
//! original FE coordinates except exact isothermal trace constraints. Mesh,
//! pressure basis and thermal penetration resolution need separate convergence.
//! Mean state is uniform, stationary and fixed; acoustic mouths are thermally
//! insulated. A conductive solid retains its original local normal-profile
//! approximation; no mean heating, flow, species evolution or time step is added.

use super::boundary_layer::thermal_wall::ThermalWall;
use super::boundary_layer::{BoundaryLayerGeometry, BoundaryLayerTransport};
use super::multiport::{
    BoundaryMode, InteriorMultiport, InteriorPortResponse, PreparedInteriorMultiport,
};
use super::{AcousticFeError, HelmholtzFe, HelmholtzModes, tri_quadrature};
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;
use std::collections::BTreeSet;

/// Explicit dense thermal cap; larger domains need a separate reduced heat basis.
pub const MAX_THERMAL_DIM: usize = 640;

/// Source-evaluated mean calorics and transport. The FE supplies rho and c.
#[derive(Clone, Copy, Debug)]
pub struct CaloricGas {
    /// Equilibrium pressure [Pa].
    pub pressure_pa: f64,
    /// Common equilibrium gas/solid temperature [K].
    pub temperature_k: f64,
    /// Equilibrium cp/gamma and original shear/conduction coefficients.
    pub transport: BoundaryLayerTransport,
}

/// One original pressure multiport and its resolved, work-conjugate heat field.
pub struct CaloricInterior {
    ports: InteriorMultiport,
    gas: CaloricGas,
    wall: ThermalWall,
    geometry: BoundaryLayerGeometry,
    thermal_dofs: Vec<usize>,
    original_dofs: usize,
    heat_k: Vec<f64>,
    heat_m: Vec<f64>,
    heat_c: Vec<f64>,
    overlap: Vec<f64>,
}

/// Frequency-local factor and heat elimination, shared by all modal excitations.
pub struct PreparedCaloricInterior<'a> {
    model: &'a CaloricInterior,
    pressure: PreparedInteriorMultiport<'a>,
    temperature_from_pressure: Vec<C64>,
    frequency_hz: f64,
    cp: C64,
    wall_admittance: C64,
}

/// Positive cycle-averaged dissipation, each mechanism counted once [W].
#[derive(Clone, Copy, Debug, Default)]
pub struct CaloricPower {
    /// Internal molecular reservoir relaxation in the solved temperature field.
    pub molecular_w: f64,
    /// Fourier conduction within the resolved gas volume.
    pub gas_conduction_w: f64,
    /// Fourier conduction within the source solid normal profile.
    pub solid_conduction_w: f64,
    /// Original tangential viscous boundary-layer work.
    pub viscous_wall_w: f64,
}

impl CaloricPower {
    #[must_use]
    pub fn total_w(self) -> f64 {
        self.molecular_w + self.gas_conduction_w + self.solid_conduction_w + self.viscous_wall_w
    }
}

/// Spatial acoustic response plus the actual temperature perturbation field.
pub struct CaloricResponse {
    /// Original work-conjugate port pressures/flows and pressure amplitudes.
    pub acoustic: InteriorPortResponse,
    /// Temperature [K] in the original FE coordinate order; fixed entries are zero.
    pub temperature_k: Vec<C64>,
    /// Separate loss terms whose sum equals mean input acoustic work.
    pub dissipation: CaloricPower,
}

fn invalid(what: &'static str) -> AcousticFeError {
    AcousticFeError::InvalidInput { what }
}
fn poll(cancelled: &dyn Fn() -> bool) -> Result<(), AcousticFeError> {
    if cancelled() {
        Err(invalid("caloric interior cancelled"))
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn linear(error: LuComplexError, frequency_hz: f64) -> AcousticFeError {
    match error {
        LuComplexError::Cancelled => invalid("caloric interior cancelled"),
        LuComplexError::Singular { .. } => AcousticFeError::SingularAtFrequency { frequency_hz },
    }
}
fn quadratic(matrix: &[f64], x: &[C64]) -> Result<f64, AcousticFeError> {
    let n = x.len();
    let mut value = 0.0;
    let mut absolute = 0.0;
    for i in 0..n {
        for j in 0..n {
            let term = matrix[i * n + j] * (x[i].conj() * x[j]).re;
            value += term;
            absolute += term.abs();
        }
    }
    if !value.is_finite() || !absolute.is_finite() || value < -128.0 * f64::EPSILON * absolute {
        return Err(invalid("caloric dissipation form overflow or nonpassivity"));
    }
    Ok(value)
}

impl CaloricInterior {
    /// Prepare the original FE pressure/temperature geometry. All selected wall
    /// faces share one declared thermal condition. Other faces have zero heat
    /// flux. Acoustic modes may overlap each other but cannot overlap these walls.
    /// The positive heat-capacity response is supplied separately per frequency.
    pub fn new(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[BoundaryMode],
        wall_faces: &[[u32; 3]],
        gas: CaloricGas,
        wall: ThermalWall,
        geometry: BoundaryLayerGeometry,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, AcousticFeError> {
        poll(cancelled)?;
        let medium = fe.medium();
        let t = gas.transport;
        if !gas.pressure_pa.is_finite()
            || gas.pressure_pa <= 0.0
            || !gas.temperature_k.is_finite()
            || gas.temperature_k <= 0.0
            || !t.cp_j_kg_k.is_finite()
            || t.cp_j_kg_k <= 0.0
            || !t.heat_capacity_ratio.is_finite()
            || t.heat_capacity_ratio <= 1.0
            || !t.thermal_conductivity_w_m_k.is_finite()
            || t.thermal_conductivity_w_m_k < 0.0
        {
            return Err(invalid("caloric interior mean gas"));
        }
        let r = t.cp_j_kg_k * (1.0 - 1.0 / t.heat_capacity_ratio);
        for ratio in [
            gas.pressure_pa / (medium.density_kg_m3 * r * gas.temperature_k),
            medium.sound_speed_m_s.powi(2) / (t.heat_capacity_ratio * r * gas.temperature_k),
        ] {
            if !ratio.is_finite() || (ratio - 1.0).abs() > 128.0 * f64::EPSILON {
                return Err(invalid(
                    "caloric gas does not match original FE equilibrium",
                ));
            }
        }
        // The original pressure wall retains ONLY its viscous form. Temperature
        // and molecular loss are accounted by the coupled volume equation below.
        let viscous = BoundaryLayerTransport {
            thermal_conductivity_w_m_k: 0.0,
            ..t
        };
        let pressure = InteriorMultiport::new_with_boundary_layers(
            fe, modes, ports, false, wall_faces, viscous, geometry, cancelled,
        )?;
        let mut fixed = vec![false; fe.ndof()];
        let quad = tri_quadrature(fe.order + 1);
        let mut seen = BTreeSet::new();
        for &face in wall_faces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !seen.insert(key) {
                return Err(invalid("duplicate caloric wall face"));
            }
            let rule = fe.face_rule(face, &quad)?;
            if wall == ThermalWall::Isothermal {
                for dof in rule.dofs {
                    fixed[dof] = true;
                }
            }
        }
        let thermal_dofs: Vec<_> = (0..fe.ndof()).filter(|&i| !fixed[i]).collect();
        let nt = thermal_dofs.len();
        if nt == 0 || nt > MAX_THERMAL_DIM {
            return Err(invalid(
                "caloric thermal dimension outside 1..=MAX_THERMAL_DIM",
            ));
        }
        let mut index = vec![usize::MAX; fe.ndof()];
        for (i, &global) in thermal_dofs.iter().enumerate() {
            index[global] = i;
        }
        let restrict = |matrix: &fs_sparse::Csr| {
            let mut out = vec![0.0; nt * nt];
            for (i, &global) in thermal_dofs.iter().enumerate() {
                let (columns, values) = matrix.row(global);
                for (&column, &value) in columns.iter().zip(values) {
                    if index[column] != usize::MAX {
                        out[i * nt + index[column]] = value;
                    }
                }
            }
            out
        };
        let heat_k = restrict(&fe.stiffness);
        let heat_m = restrict(&fe.mass);
        let heat_c = restrict(&fe.boundary_mass(wall_faces)?);
        let np = pressure.dim();
        let mut overlap = vec![0.0; nt * np];
        let mut product = vec![0.0; fe.ndof()];
        for (j, shape) in modes.shapes.iter().enumerate() {
            poll(cancelled)?;
            fe.mass.spmv(shape, &mut product);
            for (i, &global) in thermal_dofs.iter().enumerate() {
                overlap[i * np + j] = product[global];
            }
        }
        poll(cancelled)?;
        Ok(Self {
            ports: pressure,
            gas,
            wall,
            geometry,
            thermal_dofs,
            original_dofs: fe.ndof(),
            heat_k,
            heat_m,
            heat_c,
            overlap,
        })
    }

    #[must_use]
    pub fn port_mode_count(&self) -> usize {
        self.ports.port_mode_count()
    }
    #[must_use]
    pub fn pressure_dimension(&self) -> usize {
        self.ports.dim()
    }
    #[must_use]
    pub fn thermal_dimension(&self) -> usize {
        self.thermal_dofs.len()
    }

    fn operator(
        &self,
        frequency_hz: f64,
        cp: C64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(Vec<C64>, f64, Vec<C64>, C64), AcousticFeError> {
        poll(cancelled)?;
        let t = self.gas.transport;
        let r = t.cp_j_kg_k * (1.0 - 1.0 / t.heat_capacity_ratio);
        if !finite(cp) || cp.re <= r || cp.im > 0.0 {
            return Err(invalid(
                "caloric cp requires passive exp(+i omega t) response and positive cv",
            ));
        }
        let (mut a, omega_rho) = self.ports.operator(frequency_hz, cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let medium = self.ports.reduced.medium;
        let y = match self.wall {
            ThermalWall::Conductive(wall) => {
                let y = wall.normal_admittance(medium, self.geometry, frequency_hz)?;
                C64::new(y[0], y[1])
            }
            _ => C64::ZERO,
        };
        let nt = self.thermal_dimension();
        let np = self.pressure_dimension();
        let rate = C64::new(-omega_rho * cp.im, omega_rho * cp.re);
        let mut heat = Vec::with_capacity(nt * nt);
        for i in 0..nt * nt {
            heat.push(
                C64::from_re(t.thermal_conductivity_w_m_k * self.heat_k[i])
                    + rate.scale(self.heat_m[i])
                    + y.scale(self.heat_c[i]),
            );
        }
        if heat.iter().any(|&z| !finite(z)) {
            return Err(invalid("caloric heat operator overflow"));
        }
        let factor =
            lu_complex_with_cancel(&heat, nt, cancelled).map_err(|e| linear(e, frequency_hz))?;
        drop(heat);
        let mut map = vec![C64::ZERO; nt * np];
        for j in 0..np {
            poll(cancelled)?;
            let mut rhs: Vec<_> = (0..nt)
                .map(|i| C64::new(0.0, omega * self.overlap[i * np + j]))
                .collect();
            factor
                .solve_with_cancel(&mut rhs, cancelled)
                .map_err(|e| linear(e, frequency_hz))?;
            for i in 0..nt {
                map[i * np + j] = rhs[i];
            }
        }
        let mass_shift =
            (omega / medium.sound_speed_m_s).powi(2) - omega * omega_rho / self.gas.pressure_pa;
        let thermal_scale = omega * omega_rho / self.gas.temperature_k;
        for i in 0..np {
            poll(cancelled)?;
            for j in 0..np {
                let thermal = (0..nt).fold(C64::ZERO, |s, k| {
                    s + map[k * np + j].scale(self.overlap[k * np + i])
                });
                a[i * np + j] = a[i * np + j]
                    + C64::from_re(mass_shift * self.ports.reduced.m_r[i * np + j])
                    + thermal.scale(thermal_scale);
            }
        }
        if a.iter().chain(&map).any(|&z| !finite(z)) {
            return Err(invalid("caloric pressure operator overflow"));
        }
        poll(cancelled)?;
        Ok((a, omega_rho, map, y))
    }

    /// Evaluate cp with the material owner at this SAME frequency, in exp(+iwt).
    /// No total plane-wave attenuation may be substituted for this caloric law.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        cp_j_kg_k: C64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedCaloricInterior<'_>, AcousticFeError> {
        let (a, omega_rho, map, y) = self.operator(frequency_hz, cp_j_kg_k, cancelled)?;
        let pressure = self
            .ports
            .prepare_operator(a, omega_rho, frequency_hz, cancelled)?;
        Ok(PreparedCaloricInterior {
            model: self,
            pressure,
            temperature_from_pressure: map,
            frequency_hz,
            cp: cp_j_kg_k,
            wall_admittance: y,
        })
    }

    /// Original real-reference power scattering with the complete coupled heat
    /// response. Loaded pressure assembly avoids intermediate closed-cavity poles.
    pub fn scattering(
        &self,
        frequency_hz: f64,
        cp_j_kg_k: C64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, AcousticFeError> {
        self.ports.validate_references(reference_resistance)?;
        let (a, omega_rho, _, _) = self.operator(frequency_hz, cp_j_kg_k, cancelled)?;
        self.ports
            .scattering_operator(a, omega_rho, frequency_hz, reference_resistance, cancelled)
    }
}

impl PreparedCaloricInterior<'_> {
    /// Simultaneous complex modal flow with original pressure/flow coordinates.
    /// Heat dissipation follows directly from the solved thermal field:
    /// `[kappa T*K_T*T - omega rho Im(cp) T*M_T*T + Re(Yw) T*C_T*T]/(2 T0)`.
    /// The molecular and conductive terms are positive forms, not separate
    /// lossy corrections to an already complex bulk modulus.
    pub fn response(
        &self,
        inward_flow_m3_s: &[C64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<CaloricResponse, AcousticFeError> {
        let acoustic = self.pressure.response(inward_flow_m3_s, cancelled)?;
        let model = self.model;
        let np = model.pressure_dimension();
        let mut temperature = vec![C64::ZERO; model.thermal_dimension()];
        for (i, value) in temperature.iter_mut().enumerate() {
            poll(cancelled)?;
            for (j, p) in acoustic.amplitudes.iter().enumerate() {
                *value = *value + self.temperature_from_pressure[i * np + j] * C64::new(p[0], p[1]);
            }
        }
        let scale = 0.5 / model.gas.temperature_k;
        let omega_rho =
            core::f64::consts::TAU * self.frequency_hz * model.ports.reduced.medium.density_kg_m3;
        let dissipation = CaloricPower {
            molecular_w: scale
                * (-omega_rho * self.cp.im)
                * quadratic(&model.heat_m, &temperature)?,
            gas_conduction_w: scale
                * model.gas.transport.thermal_conductivity_w_m_k
                * quadratic(&model.heat_k, &temperature)?,
            solid_conduction_w: scale
                * self.wall_admittance.re
                * quadratic(&model.heat_c, &temperature)?,
            viscous_wall_w: model
                .ports
                .boundary_layer_power(self.frequency_hz, &acoustic.amplitudes)?
                .viscous_w,
        };
        if !dissipation.total_w().is_finite() {
            return Err(invalid("caloric work overflow"));
        }
        let mut temperature_k = vec![C64::ZERO; model.original_dofs];
        for (&global, &value) in model.thermal_dofs.iter().zip(&temperature) {
            temperature_k[global] = value;
        }
        poll(cancelled)?;
        Ok(CaloricResponse {
            acoustic,
            temperature_k,
            dissipation,
        })
    }
}

#[cfg(test)]
mod tests;
