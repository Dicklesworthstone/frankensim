//! Compatible thermal feedback from the original accepted frame velocity.
//!
//! A force-based member does not generally have constant axial strain: thermal
//! gradients and axial/bending coupling produce internal strain redistribution.
//! Differentiating its ORIGINAL equilibrium and Lobatto compatibility gives
//! `d_dot_i = K_i^-1 (B_i q_dot - sum_l f_T,il T_dot,il)` and
//! `q_dot = K_basic (v_dot + sum_i L*w_i B_i^T K_i^-1 sum_l f_T,il T_dot,il)`.
//! The current native basic tangent, interpolation and reference geometric
//! columns supply every mechanical operator. No displacement interpolation,
//! thermal eigenstrain rate or extra stiffness derivative is fitted here.
//!
//! Only original linear-geometry Euler–Bernoulli members are admitted. Their
//! accepted velocities may come from any supported native motion owner. This
//! instantaneous derivative does not advance that owner, solve heat transport,
//! establish positive total caloric capacity or infer finite-interval work.

use super::{
    Frame2d, FrameLaminateThermoelasticLayerHeatRate, FrameLaminateThermoelasticResponse,
    checkpoint, finite, invalid, layer_heat_rate,
};
use crate::SolidError;
use crate::frame2d::{FrameElement, ductile_damage};
use fs_exec::Cx;
use fs_exec::reduce::Compensated;

/// One original station's compatible total-strain rates and signed layer heat.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameLaminateStationMotionHeatRate {
    pub station: usize,
    pub station_coordinate: f64,
    pub reference_length_weight_m: f64,
    /// [total axial strain rate (1/s), native curvature rate (1/(m s))].
    pub generalized_strain_rate: [f64; 2],
    /// Per reference length; multiply each heat term by the station weight
    /// once to obtain its original member quadrature contribution in watts.
    pub layers: Vec<FrameLaminateThermoelasticLayerHeatRate>,
}

/// Instantaneous response to native accepted motion and original layer Tdot.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameLaminateMotionHeatRate {
    pub element: usize,
    /// Axial extension rate [m/s], then two relative end-rotation rates [1/s].
    pub basic_deformation_rate: [f64; 3],
    /// Native [N_dot, M_i_dot, M_j_dot] in [N/s, N m/s, N m/s].
    pub basic_force_rate: [f64; 3],
    pub stations: Vec<FrameLaminateStationMotionHeatRate>,
}

struct Station {
    source: FrameLaminateThermoelasticResponse,
    flexibility: [[f64; 2]; 2],
    thermal_force_rate: [f64; 2],
}

fn action2(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> Result<[f64; 2], SolidError> {
    Ok([
        finite(matrix[0][0].mul_add(vector[0], matrix[0][1] * vector[1]))?,
        finite(matrix[1][0].mul_add(vector[0], matrix[1][1] * vector[1]))?,
    ])
}

impl Frame2d {
    /// Recover compatible total-strain rates and reversible layer heat directly
    /// from this frame's accepted displacement/velocity and supplied physical
    /// temperature rates in `[station][original layer]` order.
    ///
    /// The original linear reference-coordinate owner supplies vdot; common
    /// translation, including base translation, contributes no strain rate.
    /// The existing force-based basic tangent and scaled sectional flexibility
    /// then recover thermal and axial/bending redistribution at every station.
    /// Layer temperature derivatives already include changing elastic moduli,
    /// integrated expansion and lateral relaxation; none are added again.
    /// Each original layer's thermodynamic state is evaluated once and reused
    /// for its mechanical derivative and separated reversible heat terms.
    ///
    /// Corotational, transverse-shear and separately viscous member operators
    /// refuse. Reference length, geometry, mass and inertia remain unchanged.
    /// All failures/cancellation leave accepted motion and sources untouched.
    pub fn laminate_thermoelastic_motion_heat_rate(
        &self,
        element: usize,
        temperature_rates_k_s: &[Vec<f64>],
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateMotionHeatRate, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("thermoelastic motion member index is out of range"))?;
        let original = member
            .laminate_thermal
            .as_ref()
            .ok_or_else(|| invalid("thermoelastic motion requires its original laminate member"))?;
        if original.thermoelastic.is_none()
            || member.transverse_shear.is_some()
            || member.viscous_section.is_some()
            || temperature_rates_k_s.len() != member.sections.len()
            || temperature_rates_k_s.iter().any(|rates| {
                rates.len() != original.source.layers().len()
                    || rates.iter().any(|v| !v.is_finite())
            })
        {
            return Err(invalid(
                "thermoelastic motion requires every original prepared elastic station and finite layer temperature rates",
            ));
        }
        let (dofs, columns, length) = self.reference_strain_traces(element)?;
        // Apply the original geometric columns after removing common
        // translation, retaining the actual two nodal angular velocities.
        let velocity = [
            0.0,
            0.0,
            self.vel[dofs[2]],
            self.vel[dofs[3]] - self.vel[dofs[0]],
            self.vel[dofs[4]] - self.vel[dofs[1]],
            self.vel[dofs[5]],
        ];
        if velocity.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite accepted native member velocities required"));
        }
        let mut coordinate_rates = [0.0; 3];
        for row in 0..3 {
            let sum = columns[row]
                .iter()
                .zip(velocity)
                .fold(Compensated::zero(), |sum, (&coefficient, rate)| {
                    sum.accumulate(coefficient * rate)
                });
            coordinate_rates[row] = finite(sum.value())?;
        }
        let basic_deformation_rate = [
            coordinate_rates[0],
            finite((0.5 * coordinate_rates[1] + 0.5 * coordinate_rates[2]) / length)?,
            finite((0.5 * coordinate_rates[1] - 0.5 * coordinate_rates[2]) / length)?,
        ];
        let mut right = basic_deformation_rate.map(|v| Compensated::zero().accumulate(v));
        let mut stations = Vec::with_capacity(member.sections.len());
        for (station, temperature_rates) in temperature_rates_k_s.iter().enumerate() {
            checkpoint(cx)?;
            let source = self.laminate_thermoelastic_response(element, station, cx)?;
            let section = self
                .laminate_thermal_section(element, station)
                .ok_or_else(|| invalid("original thermal station is absent"))?;
            let flexibility = ductile_damage::section_flexibility(section.prepared().tangent())?;
            let mut thermal = [Compensated::zero(); 2];
            for (layer, &temperature_rate) in source.layers.iter().zip(temperature_rates) {
                checkpoint(cx)?;
                for row in 0..2 {
                    thermal[row] = thermal[row].accumulate(finite(
                        layer.source.generalized_force_temperature_derivative[row]
                            * temperature_rate,
                    )?);
                }
            }
            let thermal_force_rate = [finite(thermal[0].value())?, finite(thermal[1].value())?];
            let thermal_strain_rate = action2(flexibility, thermal_force_rate)?;
            let b = FrameElement::b(member.xi[station]);
            for row in 0..3 {
                right[row] = right[row].accumulate(finite(
                    source.reference_length_weight_m
                        * b[0][row]
                            .mul_add(thermal_strain_rate[0], b[1][row] * thermal_strain_rate[1]),
                )?);
            }
            stations.push(Station {
                source,
                flexibility,
                thermal_force_rate,
            });
        }
        // This is the original accepted-pose tangent owner, not another
        // inversion or a stiffness inferred from averaged layer constants.
        let basic_tangent = self.element_trial(member, &self.u)?.kb;
        checkpoint(cx)?;
        let right = [
            finite(right[0].value())?,
            finite(right[1].value())?,
            finite(right[2].value())?,
        ];
        let mut basic_force_rate = [0.0; 3];
        for row in 0..3 {
            let value = basic_tangent[row]
                .iter()
                .zip(right)
                .fold(Compensated::zero(), |sum, (&coefficient, rate)| {
                    sum.accumulate(coefficient * rate)
                });
            basic_force_rate[row] = finite(value.value())?;
        }
        let mut response = Vec::with_capacity(stations.len());
        for (station, data) in stations.into_iter().enumerate() {
            checkpoint(cx)?;
            let b = FrameElement::b(member.xi[station]);
            let force_rate = [
                basic_force_rate[0] - data.thermal_force_rate[0],
                b[1][1].mul_add(basic_force_rate[1], b[1][2] * basic_force_rate[2])
                    - data.thermal_force_rate[1],
            ];
            let generalized_strain_rate = action2(data.flexibility, force_rate)?;
            let mut layers = Vec::with_capacity(data.source.layers.len());
            for (layer, &temperature_rate) in data
                .source
                .layers
                .into_iter()
                .zip(&temperature_rates_k_s[station])
            {
                checkpoint(cx)?;
                layers.push(layer_heat_rate(
                    layer.layer,
                    layer.source,
                    layer.temperature_k,
                    generalized_strain_rate,
                    temperature_rate,
                )?);
            }
            response.push(FrameLaminateStationMotionHeatRate {
                station,
                station_coordinate: data.source.station_coordinate,
                reference_length_weight_m: data.source.reference_length_weight_m,
                generalized_strain_rate,
                layers,
            });
        }
        checkpoint(cx)?;
        Ok(FrameLaminateMotionHeatRate {
            element,
            basic_deformation_rate,
            basic_force_rate,
            stations: response,
        })
    }
}
