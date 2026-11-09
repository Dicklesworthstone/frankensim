//! Passive support reflection in the original Christoffel mode basis.

use super::*;
use crate::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedanceTensor,
};

/// Fixed normal-incidence reflection and wall-motion maps. Preparation uses
/// one native 3x3 complex factorization; coherent drives allocate nothing.
#[derive(Debug)]
pub struct PreparedNormalElasticSurface<'a> {
    characteristic: NormalElasticCharacteristic<'a>,
    condition: ElasticSurfaceCondition,
    support: Option<ElasticSurfaceImpedanceTensor>,
    scattering: [[C64; 3]; 3],
    wall_modal_velocity: [[C64; 3]; 3],
}

/// Native ascending-speed quasi-mode amplitudes and physical surface fields.
/// Peak phasors use exp(-i omega t); powers are cycle means per surface area.
#[derive(Clone, Debug)]
pub struct NormalElasticSurfaceResponse {
    pub material_state_identity: ContentHash,
    pub incoming_velocity_m_s: [C64; 3],
    pub reflected_velocity_m_s: [C64; 3],
    pub incoming_normal_power_w_m2: [f64; 3],
    pub reflected_normal_power_w_m2: [f64; 3],
    pub total_incoming_power_w_m2: f64,
    pub total_reflected_power_w_m2: f64,
    /// Positive original support work, retaining its authored Gram channels.
    /// Reactive and ideal supports have exactly zero loss.
    pub absorbed_normal_power_w_m2: f64,
    /// Physical surface Cartesian velocity. Solved independently of reflection
    /// so weak motion at a stiff support survives cancelling incident waves.
    pub boundary_velocity_m_s: [C64; 3],
    /// Traction on +normal, in the same surface Cartesian axes [Pa].
    pub boundary_traction_pa: [C64; 3],
    /// Independently summed original-wave -Re(v*.traction)/2 [W/m²], before
    /// imposing the retained wall coordinates. It can lose relative accuracy
    /// for tiny absorption beside nearly complete reflection.
    pub boundary_normal_flux_w_m2: f64,
    /// Coherent ORIGINAL MATERIAL-frame Mandel strain at the surface. Like
    /// the independent flux, direct wave cancellation limits weak residuals.
    pub boundary_strain_mandel: [C64; 6],
    /// Conservative norm bound over the entire original solid half-space.
    pub strain_norm_bound: f64,
}

impl<'a> NormalElasticCharacteristic<'a> {
    /// Bind a free, clamped or passive tensor support in these surface axes.
    pub fn prepare_surface(
        &self,
        condition: ElasticSurfaceCondition,
    ) -> Result<PreparedNormalElasticSurface<'a>, MaterialError> {
        self.prepare_surface_cancellable(condition, &mut || false)
    }

    /// Bounded preparation. No partial reflection map escapes cancellation.
    pub fn prepare_surface_cancellable(
        &self,
        condition: ElasticSurfaceCondition,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedNormalElasticSurface<'a>, MaterialError> {
        poll(cancelled)?;
        let support = match condition {
            ElasticSurfaceCondition::TractionFree | ElasticSurfaceCondition::Clamped => None,
            ElasticSurfaceCondition::TensorImpedance(value) => Some(value),
            ElasticSurfaceCondition::Impedance(value) => {
                let mut matrix = [[C64::ZERO; 3]; 3];
                matrix[0][0] = value.tangential_pa_s_m();
                matrix[1][1] = value.tangential_pa_s_m();
                matrix[2][2] = value.normal_pa_s_m();
                Some(ElasticSurfaceImpedanceTensor::new(matrix)?)
            }
        };
        let mut result = PreparedNormalElasticSurface {
            characteristic: self.clone(),
            condition,
            support,
            scattering: [[C64::ZERO; 3]; 3],
            wall_modal_velocity: [[C64::ZERO; 3]; 3],
        };
        if let Some(support) = support {
            result.prepare_support(support, cancelled)?;
        } else {
            for i in 0..3 {
                if condition == ElasticSurfaceCondition::Clamped {
                    result.scattering[i][i] = -C64::ONE;
                } else {
                    result.scattering[i][i] = C64::ONE;
                    result.wall_modal_velocity[i][i] = C64::from_re(2.0);
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }
}

impl<'a> PreparedNormalElasticSurface<'a> {
    /// Original retained characteristic; it does not borrow the input owner.
    pub const fn characteristic(&self) -> &NormalElasticCharacteristic<'a> {
        &self.characteristic
    }

    pub const fn condition(&self) -> ElasticSurfaceCondition { self.condition }

    /// Reflected quasi-mode velocities per unit incident native mode velocity.
    pub fn velocity_coefficients(&self, index: usize) -> Option<[C64; 3]> {
        (index < 3).then(|| core::array::from_fn(|i| self.scattering[i][index]))
    }

    /// Drive a physical incoming Cartesian velocity in the checked surface
    /// frame. The returned modal arrays retain the original eigenbasis.
    pub fn drive_cartesian(
        &self,
        incoming: [C64; 3],
    ) -> Result<NormalElasticSurfaceResponse, MaterialError> {
        self.drive(self.characteristic.surface_velocity_to_modal(incoming))
    }

    /// Coherently apply native quasi-mode peak velocities. Source strain,
    /// positive support work and independent physical flux are all checked.
    pub fn drive(
        &self,
        incoming: [C64; 3],
    ) -> Result<NormalElasticSurfaceResponse, MaterialError> {
        let input_power = self.characteristic.powers(incoming)?;
        let reflected = core::array::from_fn(|i| complex_dot((0..3).map(|j| {
            (self.scattering[i][j], incoming[j])
        })));
        let output_power = self.characteristic.powers(reflected)?;
        let field = self.characteristic.field(incoming, reflected)?;
        let wall_modal = core::array::from_fn(|i| complex_dot((0..3).map(|j| {
            (self.wall_modal_velocity[i][j], incoming[j])
        })));
        let velocity = self.characteristic.modal_velocity_to_surface(wall_modal);
        let (traction, absorbed) = if let Some(support) = self.support {
            (support.traction(velocity), support.absorbed_power(velocity)?)
        } else if self.condition == ElasticSurfaceCondition::TractionFree {
            ([C64::ZERO; 3], 0.0)
        } else {
            (field.traction, 0.0)
        };
        let total_incoming = real_sum(input_power);
        let total_reflected = real_sum(output_power);
        if velocity.iter().chain(&traction).any(|&v| !finite(v))
            || !total_incoming.is_finite() || !total_reflected.is_finite()
            || !balanced(total_incoming, total_reflected + absorbed, total_incoming)
            || !balanced(field.flux, absorbed, total_incoming)
        {
            return Err(invalid("normal elastic support work is unrepresentable or does not close"));
        }
        Ok(NormalElasticSurfaceResponse {
            material_state_identity: self.characteristic.state.material_state_identity(),
            incoming_velocity_m_s: incoming,
            reflected_velocity_m_s: reflected,
            incoming_normal_power_w_m2: input_power,
            reflected_normal_power_w_m2: output_power,
            total_incoming_power_w_m2: total_incoming,
            total_reflected_power_w_m2: total_reflected,
            absorbed_normal_power_w_m2: absorbed,
            boundary_velocity_m_s: velocity,
            boundary_traction_pa: traction,
            boundary_normal_flux_w_m2: field.flux,
            boundary_strain_mandel: field.strain,
            strain_norm_bound: field.strain_bound,
        })
    }

    fn prepare_support(
        &mut self,
        support: ElasticSurfaceImpedanceTensor,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MaterialError> {
        let modes = self.characteristic.polarizations;
        let z = self.characteristic.impedances;
        let cartesian = support.matrix_pa_s_m();
        // E^T Z E uses compensated original products, not a second passive
        // admission. The original support factor remains the work owner.
        let mut modal = [[C64::ZERO; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
                for a in 0..3 {
                    for b in 0..3 {
                        let weight = Dd::from_f64(modes[i][a]) * Dd::from_f64(modes[j][b]);
                        re = re + weight * Dd::from_f64(cartesian[a][b].re);
                        im = im + weight * Dd::from_f64(cartesian[a][b].im);
                    }
                }
                let value = C64::new(re.to_f64(), im.to_f64());
                if !finite(value) {
                    return Err(invalid("normal elastic modal support is unrepresentable"));
                }
                modal[i][j] = value;
                modal[j][i] = value;
            }
        }
        let scales: [f64; 3] = core::array::from_fn(|i| modal[i].iter()
            .fold(z[i], |s, v| s.max(v.re.abs()).max(v.im.abs())));
        let matrix: [C64; 9] = core::array::from_fn(|entry| {
            let (i, j) = (entry / 3, entry % 3);
            divided(modal[i][j], scales[i])
                + C64::from_re(if i == j { z[i] / scales[i] } else { 0.0 })
        });
        let factor = lu_complex_with_cancel(&matrix, 3, &mut *cancelled).map_err(solve_error)?;
        for incident in 0..3 {
            // Solve reflected waves and actual wall velocity independently.
            // Subtracting a solved wall field from identity would destroy weak
            // reflection at matching; adding opposite waves loses stiff motion.
            for wall in [false, true] {
                let rhs: [C64; 3] = core::array::from_fn(|i| {
                    let incident_term = if i == incident { z[i] / scales[i] } else { 0.0 };
                    if wall {
                        C64::from_re(2.0 * incident_term)
                    } else {
                        let difference = (Dd::from_f64(if i == incident { z[i] } else { 0.0 })
                            - Dd::from_f64(modal[i][incident].re)).to_f64();
                        divided(C64::new(difference, -modal[i][incident].im), scales[i])
                    }
                });
                let mut solution = rhs;
                factor.solve_with_cancel(&mut solution, &mut *cancelled).map_err(solve_error)?;
                check_solution(&matrix, &solution, &rhs)?;
                for i in 0..3 {
                    if wall {
                        self.wall_modal_velocity[i][incident] = solution[i];
                    } else {
                        self.scattering[i][incident] = solution[i];
                    }
                }
            }
        }
        Ok(())
    }
}
