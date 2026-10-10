//! Convection controls of the checked enthalpy residual, including a direct
//! enthalpy cotangent. Never invert dT/dh: it is zero on a latent plateau.

use super::{
    EnthalpyAdjointError, EnthalpyStepGradient, EnthalpyStepLinearization, admit_linear, finite,
    poll, vector,
};
use crate::adjoint::robin::{RobinPort, bind_boundary_ports};
use crate::assemble::ASSEMBLY_TILE;
use crate::{LinearConfig, RobinFlux};
use fs_exec::Cx;

/// Complete endpoint history/source derivatives and convection controls.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpyRobinGradient {
    /// Joule-residual multiplier, previous h and consistent P1 source VJPs.
    pub transport: EnthalpyStepGradient,
    /// Derivative per kelvin of each ORIGINAL convective reference.
    pub references: Vec<f64>,
    /// Derivative per ln(HTC), with radiation properties held fixed.
    pub log_htc: Vec<f64>,
    /// Derivative per watt of assembled nodal load: dt times the multiplier.
    pub nodal_load: Vec<f64>,
}

/// A convection-only view of a privately checked physical enthalpy tangent.
/// The borrowed endpoint may contain ambient-radiation feedback; its original
/// convection boundary is retained by that producer, never supplied here.
pub struct EnthalpyRobinResponse<'a, 'm> {
    step: &'a EnthalpyStepLinearization<'m>,
    ports: Vec<RobinPort>,
    fluxes: Vec<RobinFlux>,
    linear: LinearConfig,
}

impl<'m> EnthalpyStepLinearization<'m> {
    /// Bind selected uniform convective ports in the requested order. The
    /// checked endpoint, complete tangent and immutable history are unchanged.
    /// Radiation's combined Robin coefficient/reference is never an air port.
    ///
    /// # Errors
    /// Duplicate/unknown/nonuniform ports, invalid linear work, cancellation
    /// and unrepresentable convection quadrature refuse without a response.
    pub fn robin_response<'a>(
        &'a self,
        cx: &Cx<'_>,
        regions: &[&str],
        linear: LinearConfig,
    ) -> Result<EnthalpyRobinResponse<'a, 'm>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        admit_linear(self.masses.len(), linear)?;
        let ports = bind_boundary_ports(cx, self.mesh, &self.convection_boundary, regions)?;
        let mut response = EnthalpyRobinResponse {
            step: self,
            ports,
            fluxes: Vec::new(),
            linear,
        };
        let means = response.wall_means(cx, response.temperature())?;
        for (port, mean) in response.ports.iter().zip(means) {
            poll(cx, response.fluxes.len())?;
            response.fluxes.push(RobinFlux {
                region: port.name.clone(),
                faces: port.faces.len(),
                area_m2: port.area_m2,
                mean_htc_w_per_m2_k: port.htc_w_m2_k,
                mean_wall_temperature_k: mean,
                mean_reference_temperature_k: port.reference_k,
                heat_rate_w: finite(port.htc_w_m2_k * port.area_m2 * (mean - port.reference_k))?,
            });
        }
        poll(cx, regions.len())?;
        Ok(response)
    }
}

impl<'m> EnthalpyRobinResponse<'_, 'm> {
    /// Original accepted endpoint with the complete physical state Jacobian.
    #[must_use]
    pub const fn linearization(&self) -> &EnthalpyStepLinearization<'m> {
        self.step
    }

    /// Chart-resolved accepted temperatures, not the physical history field.
    #[must_use]
    pub fn temperature(&self) -> &[f64] {
        &self.step.primal().temperature
    }

    /// Original convective ports in response order.
    #[must_use]
    pub fn ports(&self) -> &[RobinPort] {
        &self.ports
    }

    /// Convection-only heat rows, recomputed from the checked temperature field.
    #[must_use]
    pub fn robin_fluxes(&self) -> &[RobinFlux] {
        &self.fluxes
    }

    /// Area-average any nodal field on each selected convection trace.
    pub fn wall_means(&self, cx: &Cx<'_>, field: &[f64]) -> Result<Vec<f64>, EnthalpyAdjointError> {
        vector(cx, field, self.temperature().len())?;
        self.ports
            .iter()
            .map(|port| {
                let mut mean = 0.0;
                for (vertices, area) in &port.faces {
                    poll(cx, 0)?;
                    for &vertex in vertices {
                        mean = finite(mean + (area / port.area_m2 / 3.0) * field[vertex])?;
                    }
                }
                Ok(mean)
            })
            .collect()
    }

    /// Pull back direct h carry plus nodal temperature, wall-mean and outward
    /// CONVECTIVE heat-rate seeds. h carry is added after multiplying the
    /// temperature seed by dT/dh, retaining latent energy sensitivity exactly.
    /// One complete physical transpose solve owns all returned derivatives.
    /// No boundary, Krylov or nonlinear iteration is differentiated.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        h_weights: &[f64],
        nodal_temperature_weights: &[f64],
        wall_weights: &[f64],
        heat_weights: &[f64],
    ) -> Result<EnthalpyRobinGradient, EnthalpyAdjointError> {
        let n = self.temperature().len();
        vector(cx, h_weights, n)?;
        vector(cx, nodal_temperature_weights, n)?;
        vector(cx, wall_weights, self.ports.len())?;
        vector(cx, heat_weights, self.ports.len())?;
        let mut seed = nodal_temperature_weights.to_vec();
        for (index, port) in self.ports.iter().enumerate() {
            for (vertices, area) in &port.faces {
                poll(cx, index)?;
                let weight = finite(
                    (area / 3.0)
                        * (wall_weights[index] / port.area_m2
                            + heat_weights[index] * port.htc_w_m2_k),
                )?;
                for &vertex in vertices {
                    seed[vertex] = finite(seed[vertex] + weight)?;
                }
            }
        }
        let mut seed = self.step.temperature_pullback(cx, &seed)?;
        for (index, (value, &carry)) in seed.iter_mut().zip(h_weights).enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            *value = finite(*value + carry)?;
        }
        let transport = self.step.pullback(cx, &seed, self.linear)?;
        let mut nodal_load = Vec::with_capacity(n);
        for (index, &lambda) in transport.adjoint.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            nodal_load.push(finite(self.step.dt * lambda)?);
        }
        let mut references = Vec::with_capacity(self.ports.len());
        let mut log_htc = Vec::with_capacity(self.ports.len());
        for (index, port) in self.ports.iter().enumerate() {
            let conductance = finite(port.htc_w_m2_k * port.area_m2)?;
            let mut reference = finite(-heat_weights[index] * conductance)?;
            let mut coefficient = finite(heat_weights[index] * self.fluxes[index].heat_rate_w)?;
            for (vertices, area) in &port.faces {
                poll(cx, index)?;
                for (a, &row) in vertices.iter().enumerate() {
                    reference =
                        finite(reference + nodal_load[row] * port.htc_w_m2_k * (area / 3.0))?;
                    for (b, &column) in vertices.iter().enumerate() {
                        let mass = port.htc_w_m2_k * (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        coefficient = finite(
                            coefficient
                                + nodal_load[row]
                                    * mass
                                    * (port.reference_k - self.temperature()[column]),
                        )?;
                    }
                }
            }
            references.push(reference);
            log_htc.push(coefficient);
        }
        poll(cx, n)?;
        Ok(EnthalpyRobinGradient {
            transport,
            references,
            log_htc,
            nodal_load,
        })
    }
}
