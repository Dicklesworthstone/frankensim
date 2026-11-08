//! Physical point fields and disturbance-energy transport in a fluid at rest.
//!
//! The original P1 pressure and hierarchical P2 velocity/temperature are
//! evaluated at one admitted affine or curved element point. No plane-wave
//! impedance, recovered pressure gradient, smoothing or extra solve is used.
//! The outward disturbance-energy flux is
//! `J = p v - tau v + (T'/T0) q`, with `q = -kappa grad(T')` and viscous
//! stress `tau` excluding pressure. Harmonic products use peak phasors and
//! `Re(a conj(b))/2`. This is the linear thermoviscous energy corollary, not
//! the first-order heat current or a second-order streaming model; see
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_thermo.09.74.html>.
//! Local fluxes do not replace the weak work-conjugate terminal reactions.

use super::{
    C64, PreparedThermoviscousStress, ThermoviscousError, ThermoviscousFe, ThermoviscousField,
    ThermoviscousStress, ThermoviscousTimeState, finite, poll, unit_normal,
};
use fs_math::dd::Dd;

/// A physical field sample. Temperature is a perturbation, not absolute kelvin.
/// `C64` denotes peak exp(+i omega t) phasors; `f64` denotes an instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousPoint<T> {
    pub pressure_pa: T,
    pub velocity_m_s: [T; 3],
    pub temperature_k: T,
    pub temperature_gradient_k_m: [T; 3],
    /// First-order Fourier heat current, distinct from disturbance-energy flux.
    pub heat_flux_w_m2: [T; 3],
    pub stress: ThermoviscousStress<T>,
    equilibrium_temperature_k: f64,
}

/// Signed Cartesian energy-transport contributions, all in W/m^2.
/// Keeping viscous work separate avoids subtracting pressure-sized stresses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousEnergyFlux {
    pub pressure_work_w_m2: [f64; 3],
    pub viscous_work_w_m2: [f64; 3],
    pub conductive_work_w_m2: [f64; 3],
}

impl ThermoviscousEnergyFlux {
    /// Sum signed mechanisms with compensated arithmetic; no clipping.
    pub fn total_w_m2(&self) -> Result<[f64; 3], ThermoviscousError> {
        let total = std::array::from_fn(|i| {
            (Dd::from_f64(self.pressure_work_w_m2[i])
                + Dd::from_f64(self.viscous_work_w_m2[i])
                + Dd::from_f64(self.conductive_work_w_m2[i]))
            .to_f64()
        });
        if total.iter().any(|x| !x.is_finite()) {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(total)
    }

    /// Signed outward flux through a finite nonzero direction, normalized here.
    pub fn normal_w_m2(&self, direction: [f64; 3]) -> Result<f64, ThermoviscousError> {
        let normal = unit_normal(direction)?;
        let total = self.total_w_m2()?;
        let value = (0..3)
            .fold(Dd::ZERO, |a, i| {
                a + Dd::from_f64(total[i]) * Dd::from_f64(normal[i])
            })
            .to_f64();
        if !value.is_finite() {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(value)
    }
}

fn real_product(a: C64, b: C64) -> Dd {
    Dd::from_f64(a.re) * Dd::from_f64(b.re) + Dd::from_f64(a.im) * Dd::from_f64(b.im)
}

fn energy_flux(
    pressure: C64,
    velocity: [C64; 3],
    temperature_ratio: C64,
    heat_flux: [C64; 3],
    deviator: [C64; 6],
    bulk: C64,
    factor: f64,
) -> Result<ThermoviscousEnergyFlux, ThermoviscousError> {
    let tensor = [
        [deviator[0], deviator[3], deviator[5]],
        [deviator[3], deviator[1], deviator[4]],
        [deviator[5], deviator[4], deviator[2]],
    ];
    let scale = Dd::from_f64(factor);
    let value = ThermoviscousEnergyFlux {
        pressure_work_w_m2: velocity.map(|v| (scale * real_product(pressure, v)).to_f64()),
        viscous_work_w_m2: std::array::from_fn(|i| {
            // Bulk and deviatoric stress remain separate until the final sum.
            let work = (0..3).fold(real_product(bulk, velocity[i]), |a, j| {
                a + real_product(tensor[i][j], velocity[j])
            });
            (-scale * work).to_f64()
        }),
        conductive_work_w_m2: heat_flux
            .map(|q| (scale * real_product(temperature_ratio, q)).to_f64()),
    };
    value.total_w_m2()?;
    Ok(value)
}

impl ThermoviscousPoint<C64> {
    /// Cycle-averaged disturbance-energy flux, including the peak-phasor 1/2.
    pub fn mean_energy_flux(&self) -> Result<ThermoviscousEnergyFlux, ThermoviscousError> {
        energy_flux(
            self.pressure_pa,
            self.velocity_m_s,
            self.temperature_k
                .scale(1.0 / self.equilibrium_temperature_k),
            self.heat_flux_w_m2,
            self.stress.deviatoric_pa,
            self.stress.bulk_pa,
            0.5,
        )
    }
}

impl ThermoviscousPoint<f64> {
    /// Instantaneous disturbance-energy flux; no harmonic 1/2 factor.
    pub fn energy_flux(&self) -> Result<ThermoviscousEnergyFlux, ThermoviscousError> {
        energy_flux(
            C64::from_re(self.pressure_pa),
            self.velocity_m_s.map(C64::from_re),
            C64::from_re(self.temperature_k / self.equilibrium_temperature_k),
            self.heat_flux_w_m2.map(C64::from_re),
            self.stress.deviatoric_pa.map(C64::from_re),
            C64::from_re(self.stress.bulk_pa),
            1.0,
        )
    }
}

/// Fixed-size physical observation on the original mixed finite-element basis.
/// The supplied fields must belong to this model's original coefficient order.
/// Element-interface gradients and stress are traces from the selected cell.
pub struct PreparedThermoviscousPoint {
    stress: PreparedThermoviscousStress,
    basis_values: [f64; 10],
    equilibrium_temperature_k: f64,
    conductivity_w_m_k: f64,
}

impl ThermoviscousFe<'_> {
    /// Prepare pressure, velocity, temperature, stress and heat-current readout.
    /// Retains the original curved Jacobian and material transport coefficients.
    /// Warm evaluations use fixed-size storage and do not allocate.
    pub fn prepare_field_at(
        &self,
        element: usize,
        barycentric: [f64; 4],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousPoint, ThermoviscousError> {
        let stress = self.prepare_stress_at(element, barycentric, cancelled)?;
        let dofs = self.space.element_dofs(element);
        let mut basis_values = [0.0; 10];
        for (a, &(_, basis)) in dofs.iter().enumerate() {
            basis_values[a] = basis.eval(barycentric, 2);
        }
        poll(cancelled)?;
        Ok(PreparedThermoviscousPoint {
            stress,
            basis_values,
            equilibrium_temperature_k: self.medium.temperature_k,
            conductivity_w_m_k: self.medium.thermal_conductivity_w_m_k,
        })
    }
}

impl PreparedThermoviscousPoint {
    fn evaluate(
        &self,
        counts: [usize; 3],
        pressure: impl Fn(usize) -> C64,
        velocity: impl Fn(usize) -> [C64; 3],
        temperature: impl Fn(usize) -> C64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousPoint<C64>, ThermoviscousError> {
        poll(cancelled)?;
        if counts
            != [
                self.stress.pressure_count,
                self.stress.velocity_count,
                self.stress.velocity_count,
            ]
        {
            return Err(ThermoviscousError::InvalidInput(
                "point observation requires original pressure, velocity and temperature dimensions",
            ));
        }
        let stress = self
            .stress
            .evaluate(counts[0], counts[1], pressure, &velocity, cancelled)?;
        let mut v = [C64::ZERO; 3];
        let mut t = C64::ZERO;
        let mut gradient = [C64::ZERO; 3];
        for a in 0..10 {
            let index = self.stress.velocity_dofs[a];
            let va = velocity(index);
            let ta = temperature(index);
            if !finite(ta) || va.iter().any(|&x| !finite(x)) {
                return Err(ThermoviscousError::NonFinite);
            }
            t = t + ta.scale(self.basis_values[a]);
            for j in 0..3 {
                v[j] = v[j] + va[j].scale(self.basis_values[a]);
                gradient[j] = gradient[j] + ta.scale(self.stress.physical_gradients[a][j]);
            }
        }
        let heat_flux = gradient.map(|g| g.scale(-self.conductivity_w_m_k));
        if !finite(t)
            || v.iter()
                .chain(&gradient)
                .chain(&heat_flux)
                .any(|&x| !finite(x))
        {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(ThermoviscousPoint {
            pressure_pa: stress.pressure_pa,
            velocity_m_s: v,
            temperature_k: t,
            temperature_gradient_k_m: gradient,
            heat_flux_w_m2: heat_flux,
            stress,
            equilibrium_temperature_k: self.equilibrium_temperature_k,
        })
    }

    /// Observe the original complete harmonic field with its physical phase.
    pub fn harmonic(
        &self,
        field: &ThermoviscousField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousPoint<C64>, ThermoviscousError> {
        self.evaluate(
            [
                field.pressure_pa.len(),
                field.velocity_m_s.len(),
                field.temperature_k.len(),
            ],
            |i| field.pressure_pa[i],
            |i| field.velocity_m_s[i],
            |i| field.temperature_k[i],
            cancelled,
        )
    }

    /// Observe the original complete time state, including physical wall motion.
    pub fn instantaneous(
        &self,
        state: &ThermoviscousTimeState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousPoint<f64>, ThermoviscousError> {
        let value = self.evaluate(
            [
                state.pressure_pa.len(),
                state.velocity_m_s.len(),
                state.temperature_k.len(),
            ],
            |i| C64::from_re(state.pressure_pa[i]),
            |i| state.velocity_m_s[i].map(C64::from_re),
            |i| C64::from_re(state.temperature_k[i]),
            cancelled,
        )?;
        let shear_power = value.stress.shear_power * 2.0;
        let bulk_power = value.stress.bulk_power * 2.0;
        if !shear_power.is_finite() || !bulk_power.is_finite() {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(ThermoviscousPoint {
            pressure_pa: value.pressure_pa.re,
            velocity_m_s: value.velocity_m_s.map(|x| x.re),
            temperature_k: value.temperature_k.re,
            temperature_gradient_k_m: value.temperature_gradient_k_m.map(|x| x.re),
            heat_flux_w_m2: value.heat_flux_w_m2.map(|x| x.re),
            stress: ThermoviscousStress {
                pressure_pa: value.pressure_pa.re,
                deviatoric_pa: value.stress.deviatoric_pa.map(|x| x.re),
                bulk_pa: value.stress.bulk_pa.re,
                shear_power,
                bulk_power,
            },
            equilibrium_temperature_k: self.equilibrium_temperature_k,
        })
    }
}

#[cfg(test)]
mod tests;
