//! Exact normal-incidence characteristics of the original full elastic tensor.
//!
//! A checked frame maps surface Cartesian coordinates into the material's
//! original stiffness frame. Its third column is the positive wave normal.
//! The original Christoffel owner supplies all three real acoustic modes;
//! no isotropic projection, branch relabeling or second eigensolver is used.
//!
//! Quasi-mode indices are the original ascending-speed indices. The SAME
//! native polarization describes velocity for both propagation directions.
//! Thus a free surface reflects +1 and a clamp -1 for every mode. These are
//! not the isotropic interface's direction-dependent P/SV/SH polarizations.
//! At a repeated eigenvalue only the subspace and Cartesian response are
//! physical; individual degenerate mode labels carry no branch identity.
//!
//! On the wave-normal plane, t_+=-D v_+ and t_-=+D v_-, with
//! D=sum rho*c_m e_m e_m^T. The same modes give D^-1 without another inverse.
//! Energy velocity can walk off the normal in an anisotropic material; the
//! original directional-wave result retains that transport information.
//! This is a frozen lossless normal-wave owner, not oblique anisotropic
//! reflection, attenuation, prestress, finite-source or lateral-edge physics.

use super::super::DirectionalElasticWaves;
use super::surface::tensor::complex_dot;
use super::*;

mod surface;
mod welded;
pub use surface::{NormalElasticSurfaceResponse, PreparedNormalElasticSurface};
pub use welded::{NormalElasticInterfaceResponse, PreparedNormalElasticInterface};

/// Source-bound characteristic in one explicitly checked surface frame.
#[derive(Clone, Debug)]
pub struct NormalElasticCharacteristic<'a> {
    state: &'a ElasticAcousticState,
    surface_to_material: [[f64; 3]; 3],
    waves: DirectionalElasticWaves,
    polarizations: [[f64; 3]; 3],
    impedances: [f64; 3],
    impedance: [[f64; 3]; 3],
    admittance: [[f64; 3]; 3],
    incoming: [ElasticWaveTrace; 3],
}

impl<'a> NormalElasticCharacteristic<'a> {
    /// Columns are surface x, y and +normal in the ORIGINAL material frame.
    /// The finite proper orthonormal frame is checked to 1e-12 and retained
    /// without repair. The native wave owner normalizes its queried direction.
    pub fn new(
        state: &'a ElasticAcousticState,
        surface_to_material: [[f64; 3]; 3],
    ) -> Result<Self, MaterialError> {
        Self::new_cancellable(state, surface_to_material, &mut || false)
    }

    /// Bounded original Christoffel preparation, with no partial publication.
    pub fn new_cancellable(
        state: &'a ElasticAcousticState,
        surface_to_material: [[f64; 3]; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialError> {
        poll(cancelled)?;
        check_frame(surface_to_material)?;
        let normal = core::array::from_fn(|i| surface_to_material[i][2]);
        let waves = state.directional_waves(normal)?;
        poll(cancelled)?;
        let normal = waves.direction;
        let mut incoming = [empty_trace(); 3];
        let mut polarizations = [[0.0; 3]; 3];
        let impedances = waves.modes.map(|wave| wave.impedance_pa_s_m);
        for mode in 0..3 {
            poll(cancelled)?;
            let wave = waves.modes[mode];
            let p = wave.polarization;
            let c = wave.phase_speed_m_s;
            polarizations[mode] = to_surface(surface_to_material, p);
            // -sym(n tensor p)/c is the actual strain per modal velocity in
            // exp(i omega(s.x-t)). Its Mandel components stay in material axes.
            let mut strain = [0.0; 6];
            for i in 0..3 {
                strain[i] = -(normal[i] * p[i]) / c;
            }
            for (i, (a, b)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
                strain[3 + i] = -dot([normal[a], normal[b], 0.0], [p[b], p[a], 0.0])
                    * core::f64::consts::FRAC_1_SQRT_2 / c;
            }
            let stress = crate::tensor::from_mandel(&state.apply_stiffness_mandel(&strain));
            let tensor = [
                [stress[0], stress[3], stress[5]],
                [stress[3], stress[1], stress[4]],
                [stress[5], stress[4], stress[2]],
            ];
            let material_traction = tensor.map(|row| dot(row, normal));
            let traction = to_surface(surface_to_material, material_traction);
            let slowness = 1.0 / c;
            let power_factor = 0.5 * impedances[mode];
            let defect = (0..3).fold(0.0, |norm, i| {
                det::hypot(norm, traction[i] + impedances[mode] * polarizations[mode][i])
            });
            if strain.iter().chain(&traction).chain(&polarizations[mode]).any(|v| !v.is_finite())
                || !slowness.is_finite() || slowness == 0.0
                || !power_factor.is_finite() || power_factor <= 0.0
                || !defect.is_finite() || defect / impedances[mode] > 1e-10
            {
                return Err(invalid("original normal elastic trace is unresolved or unrepresentable"));
            }
            incoming[mode] = ElasticWaveTrace {
                normal_slowness_s_m: C64::from_re(slowness),
                velocity: polarizations[mode].map(C64::from_re),
                traction_pa_s_m: traction.map(C64::from_re),
                strain_mandel_s_m: strain.map(C64::from_re),
                normal_power_factor_pa_s_m: power_factor,
            };
        }
        let impedance = characteristic_tensor(polarizations, impedances)?;
        let admittance = characteristic_tensor(polarizations, impedances.map(f64::recip))?;
        poll(cancelled)?;
        Ok(Self {
            state, surface_to_material, waves, polarizations, impedances,
            impedance, admittance, incoming,
        })
    }

    /// Retained original material, including its identity and strain domain.
    pub const fn state(&self) -> &'a ElasticAcousticState { self.state }

    /// Unmodified checked surface-to-original-material frame.
    pub const fn surface_to_material(&self) -> [[f64; 3]; 3] { self.surface_to_material }

    /// Original Christoffel modes, polarizations and energy velocities in material axes.
    pub const fn directional_waves(&self) -> &DirectionalElasticWaves { &self.waves }

    /// Outer index is native ascending-speed MODE; inner index is surface x/y/z.
    pub const fn mode_polarizations_surface(&self) -> [[f64; 3]; 3] { self.polarizations }

    /// Original positive rho*c for each native mode, without a cycle-average factor.
    pub const fn mode_impedances_pa_s_m(&self) -> [f64; 3] { self.impedances }

    /// Symmetric Cartesian traction/velocity tensor [Pa s/m]. Native modal
    /// factors remain authoritative for weak directions below matrix roundoff.
    pub const fn impedance_pa_s_m(&self) -> [[f64; 3]; 3] { self.impedance }

    /// Same original modal inverse, in Cartesian surface axes [m/(Pa s)].
    pub const fn admittance_m_pa_s(&self) -> [[f64; 3]; 3] { self.admittance }

    /// Original energy-transport velocity for propagation along +normal.
    pub fn mode_energy_velocities_surface(&self) -> [[f64; 3]; 3] {
        self.waves.modes.map(|m| to_surface(self.surface_to_material, m.energy_velocity_m_s))
    }

    /// Native mode toward +normal. Velocity/traction are surface Cartesian;
    /// Mandel strain remains in ORIGINAL MATERIAL axes. Index is not ElasticMode.
    pub fn incoming_mode(&self, index: usize) -> Option<ElasticWaveTrace> {
        self.incoming.get(index).copied()
    }

    /// Same native polarization toward -normal. Traction, slowness and strain
    /// reverse sign; velocity polarization and positive power factor do not.
    pub fn reflected_mode(&self, index: usize) -> Option<ElasticWaveTrace> {
        let mut wave = self.incoming_mode(index)?;
        wave.normal_slowness_s_m = -wave.normal_slowness_s_m;
        wave.traction_pa_s_m = wave.traction_pa_s_m.map(|v| -v);
        wave.strain_mandel_s_m = wave.strain_mandel_s_m.map(|v| -v);
        Some(wave)
    }

    /// Reconstruct surface Cartesian velocity from native mode amplitudes.
    pub fn modal_velocity_to_surface(&self, modal: [C64; 3]) -> [C64; 3] {
        core::array::from_fn(|axis| complex_dot((0..3).map(|mode| {
            (C64::from_re(self.polarizations[mode][axis]), modal[mode])
        })))
    }

    /// Project surface Cartesian velocity onto the retained native modes.
    pub fn surface_velocity_to_modal(&self, surface: [C64; 3]) -> [C64; 3] {
        core::array::from_fn(|mode| complex_dot((0..3).map(|axis| {
            (C64::from_re(self.polarizations[mode][axis]), surface[axis])
        })))
    }

    fn powers(&self, amplitudes: [C64; 3]) -> Result<[f64; 3], MaterialError> {
        let mut result = [0.0; 3];
        for i in 0..3 {
            if !finite(amplitudes[i]) {
                return Err(invalid("normal elastic amplitudes must be finite"));
            }
            result[i] = power(self.incoming[i], amplitudes[i])?;
        }
        Ok(result)
    }

    // Both sides of a welded interface share +normal. The second side's
    // incoming waves therefore occupy negative, not positive.
    fn field(&self, positive: [C64; 3], negative: [C64; 3]) -> Result<NormalField, MaterialError> {
        let values: [C64; 12] = core::array::from_fn(|component| {
            complex_dot((0..3).flat_map(|m| {
                let (coefficient, parity) = match component {
                    0..3 => (self.incoming[m].velocity[component], 1.0),
                    3..6 => (self.incoming[m].traction_pa_s_m[component - 3], -1.0),
                    _ => (self.incoming[m].strain_mandel_s_m[component - 6], -1.0),
                };
                [(coefficient, positive[m]), (coefficient.scale(parity), negative[m])]
            }))
        });
        let strain_bound = real_sum((0..3).map(|m| {
            let norm = self.incoming[m].strain_mandel_s_m.iter()
                .fold(0.0, |norm, v| det::hypot(norm, v.abs()));
            norm * (positive[m].abs() + negative[m].abs())
        }));
        if !strain_bound.is_finite() || strain_bound > self.state.maximum_absolute_mandel_strain()
            || values.iter().any(|&v| !finite(v))
        {
            return Err(invalid("normal elastic field exceeds original strain domain or numerical range"));
        }
        let velocity: [C64; 3] = core::array::from_fn(|i| values[i]);
        let traction: [C64; 3] = core::array::from_fn(|i| values[3 + i]);
        let flux = -0.5 * complex_dot((0..3).map(|i| (velocity[i].conj(), traction[i]))).re;
        if !flux.is_finite() {
            return Err(invalid("normal elastic flux is unrepresentable"));
        }
        Ok(NormalField {
            velocity, traction, strain: core::array::from_fn(|i| values[6 + i]),
            flux, strain_bound,
        })
    }
}

#[derive(Clone, Copy)]
struct NormalField {
    velocity: [C64; 3],
    traction: [C64; 3],
    strain: [C64; 6],
    flux: f64,
    strain_bound: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).fold(Dd::ZERO, |sum, i| sum + Dd::from_f64(a[i]) * Dd::from_f64(b[i])).to_f64()
}

fn to_surface(frame: [[f64; 3]; 3], material: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| dot(core::array::from_fn(|j| frame[j][i]), material))
}

fn check_frame(frame: [[f64; 3]; 3]) -> Result<(), MaterialError> {
    if frame.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid("normal elastic surface frame must be finite"));
    }
    let columns: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| frame[j][i]));
    for i in 0..3 {
        for j in 0..3 {
            let value = dot(columns[i], columns[j]);
            if !value.is_finite() || (value - if i == j { 1.0 } else { 0.0 }).abs() > 1e-12 {
                return Err(invalid("normal elastic surface frame must be orthonormal"));
            }
        }
    }
    let cross = core::array::from_fn(|i| {
        columns[0][(i + 1) % 3] * columns[1][(i + 2) % 3]
            - columns[0][(i + 2) % 3] * columns[1][(i + 1) % 3]
    });
    if dot(cross, columns[2]) <= 0.0 {
        return Err(invalid("normal elastic surface frame must be right-handed"));
    }
    Ok(())
}

fn characteristic_tensor(p: [[f64; 3]; 3], values: [f64; 3]) -> Result<[[f64; 3]; 3], MaterialError> {
    if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(invalid("normal elastic characteristic is unrepresentable"));
    }
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            let value = (0..3).fold(Dd::ZERO, |sum, m| sum
                + Dd::from_f64(values[m]) * Dd::from_f64(p[m][i]) * Dd::from_f64(p[m][j])).to_f64();
            if !value.is_finite() || (i == j && value <= 0.0) {
                return Err(invalid("normal elastic characteristic tensor is unrepresentable"));
            }
            result[i][j] = value;
            result[j][i] = value;
        }
    }
    Ok(result)
}

fn check_solution(matrix: &[C64], solution: &[C64], rhs: &[C64]) -> Result<(), MaterialError> {
    let n = rhs.len();
    if solution.iter().any(|&v| !finite(v)) {
        return Err(invalid("normal elastic boundary solution is unrepresentable"));
    }
    for i in 0..n {
        let residual = (complex_dot((0..n).map(|j| (matrix[n * i + j], solution[j]))) - rhs[i]).abs();
        let magnitude = real_sum((0..n).map(|j| (matrix[n * i + j] * solution[j]).abs())) + rhs[i].abs();
        if !residual.is_finite() || !magnitude.is_finite()
            || (magnitude > 0.0 && residual / magnitude > 1e-11)
        {
            return Err(invalid("normal elastic boundary solve did not converge"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
