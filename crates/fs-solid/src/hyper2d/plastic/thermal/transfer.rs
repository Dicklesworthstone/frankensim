//! Conservative transfer between original sheet nodes and native material points.
//!
//! Temperature uses the same cached P1/Q1 shape functions that define the sheet
//! geometry. Caloric energy returns through their transpose, weighted by the
//! ORIGINAL reference volume. It includes signed reversible exchange and the
//! native physical loss once; endpoint-work defects never enter this transfer.
//! This is a spatial transfer, not a caloric storage or conduction solver.

use fs_exec::Cx;
use fs_material::plastic::finite::thermal::FiniteThermalJ2Material;

use super::super::{checkpoint, invalid};
use super::{
    FiniteThermoplasticSheet, FiniteThermoplasticSheetError, FiniteThermoplasticSheetStep,
};
use crate::hyper3::plastic::FinitePlasticStepControls;
use crate::plastic3::IncrementSettings;

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticSheet<'a, L> {
    /// Interpolate positive absolute nodal temperatures to the native points.
    /// The cell-minimum form preserves constant temperatures exactly and avoids
    /// subtracting large positive and negative weighted temperatures. Individual
    /// source-family domains remain enforced by the normal endpoint advance.
    pub fn quadrature_temperatures_from_nodes(
        &self,
        nodal_temperatures_k: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        if nodal_temperatures_k.len() != self.mesh.nodes.len()
            || nodal_temperatures_k
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid(
                "one positive finite absolute temperature per original sheet node required",
            )
            .into());
        }
        let mut temperatures = Vec::with_capacity(self.states.len());
        for element in self.elements.iter() {
            let minimum = element.nodes[..element.arity]
                .iter()
                .map(|&node| nodal_temperatures_k[node])
                .fold(f64::INFINITY, f64::min);
            for point in &element.points {
                checkpoint(cx)?;
                let shape = &point.shape[..element.arity];
                let mut temperature = minimum;
                for (a, &weight) in shape.iter().enumerate() {
                    temperature = weight.mul_add(
                        nodal_temperatures_k[element.nodes[a]] - minimum,
                        temperature,
                    );
                }
                if !temperature.is_finite() || temperature <= 0.0 {
                    return Err(invalid(
                        "interpolated sheet temperature is not positive and representable",
                    )
                    .into());
                }
                temperatures.push(temperature);
            }
        }
        checkpoint(cx)?;
        Ok(temperatures)
    }

    /// Advance the existing thermal sheet with temperatures from original
    /// nodes. This uses the same source admission, actual duration, nonlinear
    /// solve and atomic physical-account gates as the point-temperature API.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_from_nodal_temperatures(
        &self,
        nodal_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetStep<'a, L>, FiniteThermoplasticSheetError> {
        let temperatures = self.quadrature_temperatures_from_nodes(nodal_temperatures_k, cx)?;
        self.advance_with_quadrature_temperatures(
            &temperatures,
            prescribed_m,
            nodal_forces_n,
            dt_s,
            settings,
            controls,
            cx,
        )
    }
}

impl<L: FiniteThermalJ2Material> FiniteThermoplasticSheetStep<'_, L> {
    /// Signed caloric energy in joules for every original sheet node.
    ///
    /// Apply once to the receiving thermal owner. The transfer is the transpose
    /// of the native temperature interpolation: `Q_a = sum_q N_a(q) V0(q) h_q`.
    /// Reversible cooling remains negative; neither absolute values, clipping,
    /// current-volume factors nor numerical work defects are introduced.
    pub fn nodal_caloric_heat_j(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        if self.quadrature_energy.len() != self.sheet.states.len() {
            return Err(invalid(
                "sheet caloric transfer requires the complete native point account",
            )
            .into());
        }
        let mut heat = vec![0.0; self.sheet.mesh.nodes.len()];
        for element in self.sheet.elements.iter() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let shape = &point.shape[..element.arity];
                let density = self.quadrature_energy[element.offset + q].caloric_heat_j_m3();
                let point_heat = point.volume_m3 * density;
                if !density.is_finite() || !point_heat.is_finite() {
                    return Err(invalid("native sheet caloric energy is not representable").into());
                }
                for (a, &weight) in shape.iter().enumerate() {
                    let node = element.nodes[a];
                    heat[node] = weight.mul_add(point_heat, heat[node]);
                    if !heat[node].is_finite() {
                        return Err(
                            invalid("nodal sheet caloric energy is not representable").into()
                        );
                    }
                }
            }
        }
        checkpoint(cx)?;
        Ok(heat)
    }
}
