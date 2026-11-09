//! Spatial physical dissipation of the original resolved mixed fields.
//!
//! The disturbance-energy corollary includes both viscous caloric heating and
//! destruction of thermal availability. Only the former is directly a mean
//! heat source: Fourier conduction and complex-Cp relaxation require their own
//! mean-energy closure before coupling them to a caloric thermal equation.

use super::{
    PreparedThermoviscous, PreparedThermoviscousScattering, ThermoviscousError, finite, poll,
};
use crate::highorder::simplex::duffy_quadrature;
use fs_math::{c64::C64, dd::Dd};

/// Cycle-mean volume integrals for ONE original tetrahedron, in SI units.
/// Peak phasors use `exp(+i omega t)`. Geometry, density, transport, equilibrium
/// temperature and complex heat capacity belong to the immutable prepared
/// owner. No boundary or structural loss is included.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousElementPower {
    /// Original physical volume, using the retained geometry quadrature [m^3].
    pub volume_m3: f64,
    /// Original equilibrium fluid mass `rho * volume_m3` [kg].
    pub mass_kg: f64,
    /// `integral(mu * |dev(sym(grad v))|^2)` [W].
    pub shear_heating_w: f64,
    /// `integral(zeta * |div(v)|^2 / 2)` [W].
    pub bulk_heating_w: f64,
    /// `integral(kappa * |grad(T')|^2 / (2*T0^2))` [W/K].
    pub conductive_entropy_w_k: f64,
    /// `T0 * conductive_entropy_w_k` [W], distinct from caloric heating.
    pub thermal_availability_loss_w: f64,
    /// `integral(-omega*rho*Im(Cp)*|T'|^2/(2*T0))` [W]. This is
    /// disturbance-energy loss, not an independently closed mean heat source.
    pub caloric_availability_loss_w: f64,
}

impl ThermoviscousElementPower {
    /// Direct irreversible shear-plus-bulk caloric heating [W]. Thermal
    /// availability destruction must not be added to this mean heat source.
    #[must_use]
    pub fn viscous_heat_w(self) -> f64 {
        self.shear_heating_w + self.bulk_heating_w
    }

    /// Full local loss in the native acoustic disturbance-energy balance [W].
    #[must_use]
    pub fn disturbance_dissipation_w(self) -> f64 {
        self.viscous_heat_w() + self.thermal_availability_loss_w + self.caloric_availability_loss_w
    }
}

/// Preserve the original full-domain accumulation order, rather than summing
/// rounded per-element totals when assembling an existing response.
#[derive(Default)]
pub(super) struct VolumeDissipation {
    pub(super) viscous_w: f64,
    pub(super) bulk_w: f64,
    pub(super) conductive_w: f64,
    pub(super) caloric_w: f64,
}

impl PreparedThermoviscous<'_, '_> {
    /// Original volume-cell count, available before allocating or integrating
    /// a spatial loss observation. This includes every retained tetrahedron.
    #[must_use]
    pub fn element_count(&self) -> usize {
        self.model.space.complex.tets.len()
    }

    /// Independently integrate the declared full P2 velocity and temperature
    /// coefficients on every original tetrahedron, in `TetComplex::tets` order.
    /// Vectors use the native hierarchical ordering: vertices followed by edge
    /// bubbles; temperature is the actual perturbation in K, not `T'/T0`.
    /// Pass a response from THIS prepared owner to account for that solution.
    /// Arbitrary finite coefficients are also admitted as a field observation;
    /// this does not certify that they solve the mixed equations or boundaries.
    ///
    /// Uses original affine/quadratic geometry and the SAME quadrature kernel
    /// as `response`, including any already-lifted moving-wall velocity. No
    /// factorization, attenuation inference, pressure-exchange subtraction or
    /// clipped balance residual is used. Original model work limits apply;
    /// cancellation or invalid coefficients return no partial vector.
    pub fn element_dissipation(
        &self,
        velocity_m_s: &[[C64; 3]],
        temperature_k: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ThermoviscousElementPower>, ThermoviscousError> {
        let mut elements = Vec::with_capacity(self.model.space.complex.tets.len());
        self.volume_dissipation(velocity_m_s, temperature_k, cancelled, |power| {
            elements.push(power);
        })?;
        Ok(elements)
    }

    pub(super) fn volume_dissipation(
        &self,
        velocity_m_s: &[[C64; 3]],
        temperature_k: &[C64],
        cancelled: &mut impl FnMut() -> bool,
        mut element: impl FnMut(ThermoviscousElementPower),
    ) -> Result<VolumeDissipation, ThermoviscousError> {
        poll(cancelled)?;
        let model = self.model;
        if velocity_m_s.len() != model.space.ndof || temperature_k.len() != model.space.ndof {
            return Err(ThermoviscousError::InvalidInput(
                "full P2 velocity/temperature dimensions",
            ));
        }
        for (i, (velocity, &temperature)) in velocity_m_s.iter().zip(temperature_k).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if velocity.iter().any(|&v| !finite(v)) || !finite(temperature) {
                return Err(ThermoviscousError::NonFinite);
            }
        }
        let mut total = VolumeDissipation::default();
        let quad = duffy_quadrature(model.volume_quadrature_order());
        for t in 0..model.space.complex.tets.len() {
            poll(cancelled)?;
            let dofs = model.space.element_dofs(t);
            let mut local = VolumeDissipation::default();
            let mut volume = Dd::ZERO;
            for &(bary, w) in &quad {
                poll(cancelled)?;
                let (jac, barycentric_gradients) = model.volume_geometry(t, bary)?;
                let mut grad_v = [[C64::ZERO; 3]; 3];
                let mut grad_t = [C64::ZERO; 3];
                let mut temp = C64::ZERO;
                for &(i, f) in &dofs {
                    let dl = f.d_lambda(bary, 2);
                    let gradient = std::array::from_fn::<_, 3, _>(|c| {
                        (0..4)
                            .map(|a| dl[a] * barycentric_gradients[a][c])
                            .sum::<f64>()
                    });
                    temp = temp + temperature_k[i].scale(f.eval(bary, 2));
                    for a in 0..3 {
                        grad_t[a] = grad_t[a] + temperature_k[i].scale(gradient[a]);
                        for b in 0..3 {
                            grad_v[a][b] = grad_v[a][b] + velocity_m_s[i][a].scale(gradient[b]);
                        }
                    }
                }
                let divergence = grad_v[0][0] + grad_v[1][1] + grad_v[2][2];
                let mut dev_norm = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        let dev = (grad_v[a][b] + grad_v[b][a]).scale(0.5)
                            - if a == b {
                                divergence.scale(1.0 / 3.0)
                            } else {
                                C64::ZERO
                            };
                        dev_norm += dev.norm_sq();
                    }
                }
                let volume_weight = jac * w;
                let gas = model.medium;
                let shear = volume_weight * gas.dynamic_viscosity_pa_s * dev_norm;
                let bulk = 0.5 * volume_weight * gas.bulk_viscosity_pa_s * divergence.norm_sq();
                let conduction = 0.5 * volume_weight * gas.thermal_conductivity_w_m_k
                    / gas.temperature_k
                    * grad_t.iter().map(|z| z.norm_sq()).sum::<f64>();
                let caloric = -0.5
                    * volume_weight
                    * std::f64::consts::TAU
                    * self.frequency_hz
                    * gas.density_kg_m3
                    * self.cp.im
                    / gas.temperature_k
                    * temp.norm_sq();
                // These four additions retain the original response's order
                // and point arithmetic. Per-element accumulation is separate.
                total.viscous_w += shear;
                total.bulk_w += bulk;
                total.conductive_w += conduction;
                total.caloric_w += caloric;
                local.viscous_w += shear;
                local.bulk_w += bulk;
                local.conductive_w += conduction;
                local.caloric_w += caloric;
                volume = volume + Dd::from_f64(volume_weight);
            }
            let volume_m3 = volume.to_f64();
            let mass_kg = model.medium.density_kg_m3 * volume_m3;
            let power = ThermoviscousElementPower {
                volume_m3,
                mass_kg,
                shear_heating_w: local.viscous_w,
                bulk_heating_w: local.bulk_w,
                conductive_entropy_w_k: local.conductive_w / model.medium.temperature_k,
                thermal_availability_loss_w: local.conductive_w,
                caloric_availability_loss_w: local.caloric_w,
            };
            if [
                power.volume_m3,
                power.mass_kg,
                power.shear_heating_w,
                power.bulk_heating_w,
                power.conductive_entropy_w_k,
                power.thermal_availability_loss_w,
                power.caloric_availability_loss_w,
                power.disturbance_dissipation_w(),
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(ThermoviscousError::NonFinite);
            }
            if volume_m3 <= 0.0 || mass_kg <= 0.0 {
                return Err(ThermoviscousError::InvalidInput(
                    "unrepresentable positive element volume or mass",
                ));
            }
            element(power);
        }
        if [
            total.viscous_w,
            total.bulk_w,
            total.conductive_w,
            total.caloric_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(total)
    }
}

impl PreparedThermoviscousScattering<'_, '_> {
    /// Original element loss observation for fields of this loaded scattering
    /// owner. Uses its retained frequency, complex Cp and original geometry;
    /// the same full P2 ordering and heat/availability distinction apply as in
    /// [`PreparedThermoviscous::element_dissipation`].
    pub fn element_dissipation(
        &self,
        velocity_m_s: &[[C64; 3]],
        temperature_k: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ThermoviscousElementPower>, ThermoviscousError> {
        self.inner
            .element_dissipation(velocity_m_s, temperature_k, cancelled)
    }
}
