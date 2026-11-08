//! Resolved Cauchy stress and wall shear on the original P2/P1 fluid fields.
//!
//! The constitutive law is the same Newtonian law as the mixed volume operator:
//! `sigma = -p I + 2 mu dev(sym(grad v)) + zeta div(v) I`. See Kampinga,
//! *Viscothermal acoustics using finite elements* (2010), section 2.4,
//! <https://doi.org/10.3990/1.9789036530507>.
//! Preparation retains the original affine/quadratic physical gradients. Warm
//! observations read only ten velocity and four pressure coefficients, without
//! allocations, finite differences or another solve. Interface stress is the
//! trace from the selected element; no stress smoothing or mesh-error claim is
//! implied. Pointwise tractions do not replace the weak, work-conjugate reaction
//! used by the moving-wall coupling.

use super::time::ThermoviscousTimeState;
use super::{C64, ThermoviscousError, ThermoviscousFe, ThermoviscousField, finite, poll};

/// Original pressure, velocity and temperature fields and their energy flux.
pub mod field;

/// Physical stress with tensor shear in `[xx, yy, zz, xy, yz, zx]` order.
/// `C64` denotes peak phasors with exp(+i omega t); `f64` denotes an instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousStress<T> {
    /// Thermodynamic pressure perturbation [Pa], positive in compression.
    pub pressure_pa: T,
    /// Trace-free viscous stress [Pa]; includes normal deviatoric components.
    pub deviatoric_pa: [T; 6],
    /// Independent spherical viscous stress zeta div(v) [Pa].
    pub bulk_pa: T,
    shear_power: f64,
    bulk_power: f64,
}

/// Cauchy traction `sigma n` on the side with outward unit normal `n`.
/// The force on the other side has the opposite sign.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StressTraction<T> {
    /// Total traction [Pa].
    pub total_pa: [T; 3],
    /// Signed normal stress [Pa], positive in tension.
    pub normal_pa: T,
    /// Tangential traction [Pa], formed without subtracting pressure-sized terms.
    pub shear_pa: [T; 3],
}

impl ThermoviscousStress<C64> {
    fn from_gradient(
        pressure_pa: C64,
        gradient: [[C64; 3]; 3],
        mu: f64,
        bulk: f64,
    ) -> Result<Self, ThermoviscousError> {
        if !finite(pressure_pa) || gradient.iter().flatten().any(|&v| !finite(v)) {
            return Err(ThermoviscousError::NonFinite);
        }
        let normal = [gradient[0][0], gradient[1][1], gradient[2][2]];
        // Differences preserve a small deviator on a large spherical rate and
        // make pure isotropic expansion exactly shear-free.
        let rate = [
            (normal[0] - normal[1]).scale(1.0 / 3.0) + (normal[0] - normal[2]).scale(1.0 / 3.0),
            (normal[1] - normal[0]).scale(1.0 / 3.0) + (normal[1] - normal[2]).scale(1.0 / 3.0),
            (normal[2] - normal[0]).scale(1.0 / 3.0) + (normal[2] - normal[1]).scale(1.0 / 3.0),
            gradient[0][1].scale(0.5) + gradient[1][0].scale(0.5),
            gradient[1][2].scale(0.5) + gradient[2][1].scale(0.5),
            gradient[2][0].scale(0.5) + gradient[0][2].scale(0.5),
        ];
        let divergence = normal[0] + normal[1] + normal[2];
        let deviatoric_pa = rate.map(|v| v.scale(mu).scale(2.0));
        let bulk_pa = divergence.scale(bulk);
        // Mean power is Re(tau : conj(D))/2, with tensor shear counted twice.
        // Multiply stress by rate instead of squaring the rate before mu.
        let shear_power = (0..6)
            .map(|i| {
                let weight = if i < 3 { 0.5 } else { 1.0 };
                weight * (deviatoric_pa[i].re * rate[i].re + deviatoric_pa[i].im * rate[i].im)
            })
            .sum::<f64>();
        let bulk_power = 0.5 * (bulk_pa.re * divergence.re + bulk_pa.im * divergence.im);
        if deviatoric_pa.iter().any(|&v| !finite(v))
            || !finite(bulk_pa)
            || !shear_power.is_finite()
            || !bulk_power.is_finite()
        {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(Self {
            pressure_pa,
            deviatoric_pa,
            bulk_pa,
            shear_power,
            bulk_power,
        })
    }

    /// Cycle-averaged physical shear heating [W/m^3], including the 1/2 factor.
    #[must_use]
    pub const fn mean_shear_heating_w_m3(&self) -> f64 {
        self.shear_power
    }

    /// Cycle-averaged independent bulk-viscous heating [W/m^3].
    #[must_use]
    pub const fn mean_bulk_heating_w_m3(&self) -> f64 {
        self.bulk_power
    }

    /// Resolve traction for a finite nonzero normal direction (normalized here).
    /// Pressure and bulk viscosity contribute only to normal traction.
    pub fn traction(&self, direction: [f64; 3]) -> Result<StressTraction<C64>, ThermoviscousError> {
        let n = unit_normal(direction)?;
        let s = &self.deviatoric_pa;
        let tensor = [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]];
        let dev_traction: [C64; 3] =
            std::array::from_fn(|i| (0..3).fold(C64::ZERO, |a, j| a + tensor[i][j].scale(n[j])));
        let dev_normal = (0..3).fold(C64::ZERO, |a, i| a + dev_traction[i].scale(n[i]));
        let spherical = self.bulk_pa - self.pressure_pa;
        let normal_pa = dev_normal + spherical;
        let shear_pa = std::array::from_fn(|i| dev_traction[i] - dev_normal.scale(n[i]));
        let total_pa = std::array::from_fn(|i| dev_traction[i] + spherical.scale(n[i]));
        if !finite(normal_pa) || total_pa.iter().chain(&shear_pa).any(|&v| !finite(v)) {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(StressTraction {
            total_pa,
            normal_pa,
            shear_pa,
        })
    }
}

impl ThermoviscousStress<f64> {
    /// Instantaneous physical shear heating [W/m^3]; no harmonic 1/2 factor.
    #[must_use]
    pub const fn shear_heating_w_m3(&self) -> f64 {
        self.shear_power
    }

    /// Instantaneous independent bulk-viscous heating [W/m^3].
    #[must_use]
    pub const fn bulk_heating_w_m3(&self) -> f64 {
        self.bulk_power
    }

    /// Instantaneous traction for a finite nonzero normal direction.
    pub fn traction(&self, direction: [f64; 3]) -> Result<StressTraction<f64>, ThermoviscousError> {
        let complex = ThermoviscousStress {
            pressure_pa: C64::from_re(self.pressure_pa),
            deviatoric_pa: self.deviatoric_pa.map(C64::from_re),
            bulk_pa: C64::from_re(self.bulk_pa),
            shear_power: self.shear_power,
            bulk_power: self.bulk_power,
        }
        .traction(direction)?;
        Ok(StressTraction {
            total_pa: complex.total_pa.map(|v| v.re),
            normal_pa: complex.normal_pa.re,
            shear_pa: complex.shear_pa.map(|v| v.re),
        })
    }
}

fn unit_normal(direction: [f64; 3]) -> Result<[f64; 3], ThermoviscousError> {
    let scale = direction.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if direction.iter().any(|x| !x.is_finite()) || scale == 0.0 {
        return Err(ThermoviscousError::InvalidInput(
            "finite nonzero traction normal required",
        ));
    }
    let scaled = direction.map(|x| x / scale);
    let norm = fs_math::det::sqrt(scaled.iter().map(|x| x * x).sum());
    Ok(scaled.map(|x| x / norm))
}

/// One physical point bound to the original mixed mesh and uniform transport.
/// The supplied field/state must use this model's original coefficient ordering.
pub struct PreparedThermoviscousStress {
    pressure_vertices: [usize; 4],
    barycentric: [f64; 4],
    velocity_dofs: [usize; 10],
    physical_gradients: [[f64; 3]; 10],
    pressure_count: usize,
    velocity_count: usize,
    mu: f64,
    bulk: f64,
}

impl ThermoviscousFe<'_> {
    /// Prepare an interior or boundary point of one original tetrahedron.
    /// Curved meshes use the admitted quadratic Jacobian at that same point.
    /// Preparation is cold; repeated field observations need no allocation.
    pub fn prepare_stress_at(
        &self,
        element: usize,
        barycentric: [f64; 4],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousStress, ThermoviscousError> {
        poll(cancelled)?;
        if element >= self.space.complex.tets.len()
            || barycentric
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0 || *x > 1.0)
            || (barycentric.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(ThermoviscousError::InvalidInput(
                "stress point must lie in an original element",
            ));
        }
        let (_, gradients) = self.volume_geometry(element, barycentric)?;
        let dofs = self.space.element_dofs(element);
        let mut velocity_dofs = [0; 10];
        let mut physical_gradients = [[0.0; 3]; 10];
        for (a, &(index, basis)) in dofs.iter().enumerate() {
            velocity_dofs[a] = index;
            let derivative = basis.d_lambda(barycentric, 2);
            physical_gradients[a] =
                std::array::from_fn(|j| (0..4).map(|i| derivative[i] * gradients[i][j]).sum());
        }
        if physical_gradients.iter().flatten().any(|x| !x.is_finite()) {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(PreparedThermoviscousStress {
            pressure_vertices: self.space.complex.tets[element].map(|v| v as usize),
            barycentric,
            velocity_dofs,
            physical_gradients,
            pressure_count: self.space.complex.vertex_count,
            velocity_count: self.space.ndof,
            mu: self.medium.dynamic_viscosity_pa_s,
            bulk: self.medium.bulk_viscosity_pa_s,
        })
    }
}

impl PreparedThermoviscousStress {
    fn evaluate(
        &self,
        pressure_count: usize,
        velocity_count: usize,
        pressure: impl Fn(usize) -> C64,
        velocity: impl Fn(usize) -> [C64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousStress<C64>, ThermoviscousError> {
        poll(cancelled)?;
        if pressure_count != self.pressure_count || velocity_count != self.velocity_count {
            return Err(ThermoviscousError::InvalidInput(
                "stress observation requires original mixed field dimensions",
            ));
        }
        let mut p = C64::ZERO;
        for i in 0..4 {
            p = p + pressure(self.pressure_vertices[i]).scale(self.barycentric[i]);
        }
        let mut gradient = [[C64::ZERO; 3]; 3];
        for a in 0..10 {
            let v = velocity(self.velocity_dofs[a]);
            for i in 0..3 {
                for j in 0..3 {
                    gradient[i][j] = gradient[i][j] + v[i].scale(self.physical_gradients[a][j]);
                }
            }
        }
        ThermoviscousStress::from_gradient(p, gradient, self.mu, self.bulk)
    }

    /// Observe peak stress phasors and cycle-averaged shear/bulk heating.
    pub fn harmonic(
        &self,
        field: &ThermoviscousField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousStress<C64>, ThermoviscousError> {
        self.evaluate(
            field.pressure_pa.len(),
            field.velocity_m_s.len(),
            |i| field.pressure_pa[i],
            |i| field.velocity_m_s[i],
            cancelled,
        )
    }

    /// Observe instantaneous stress and heating from the full original time state.
    /// Stored moving-wall velocity is already in the physical P2 coefficients.
    pub fn instantaneous(
        &self,
        state: &ThermoviscousTimeState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousStress<f64>, ThermoviscousError> {
        let value = self.evaluate(
            state.pressure_pa.len(),
            state.velocity_m_s.len(),
            |i| C64::from_re(state.pressure_pa[i]),
            |i| state.velocity_m_s[i].map(C64::from_re),
            cancelled,
        )?;
        let shear_power = 2.0 * value.shear_power;
        let bulk_power = 2.0 * value.bulk_power;
        if !shear_power.is_finite() || !bulk_power.is_finite() {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(ThermoviscousStress {
            pressure_pa: value.pressure_pa.re,
            deviatoric_pa: value.deviatoric_pa.map(|v| v.re),
            bulk_pa: value.bulk_pa.re,
            shear_power,
            bulk_power,
        })
    }
}

#[cfg(test)]
mod tests;
