//! Net depth heating from exterior irradiation and original thermal layers.

use super::{
    AbsorbingLayer, AbsorbingSlab, Arc, ConductionError, GL8, GL16, SlabQuadrature, add,
    checkpoint, difference_product, finite, invalid, layered, nonnegative,
};

/// Numerical angular diagnostics, not certified continuum error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabDepthReport {
    pub controls: SlabQuadrature,
    pub angular_evaluations: usize,
    pub accepted_panels: usize,
    /// Largest sum of accepted GL8/GL16 discrepancies over retained couplings.
    pub maximum_estimated_coefficient_error: f64,
}

/// Prepared net radiative heating within one original uniform-temperature layer.
///
/// Positive coefficients couple an interval to both exterior faces and every
/// other original layer. Detailed balance removes the target layer's own
/// emission and reabsorption analytically. Uniform equilibrium therefore gives
/// exactly zero, without subtracting independently integrated emission fields.
/// The full conserved Snell domain includes all propagating and trapped modes
/// with the thermal n² measure. Sources and physical interfaces remain shared.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabDepthResponse {
    layers: Arc<[AbsorbingLayer]>,
    layer: usize,
    interval_m: [f64; 2],
    face_coupling: [f64; 2],
    layer_coupling: Arc<[f64]>,
    report: SlabDepthReport,
}

impl SlabDepthResponse {
    pub fn layers(&self) -> &[AbsorbingLayer] {
        &self.layers
    }
    pub const fn layer_index(&self) -> usize {
        self.layer
    }
    pub const fn interval_m(&self) -> [f64; 2] {
        self.interval_m
    }
    pub const fn face_coupling(&self) -> [f64; 2] {
        self.face_coupling
    }
    /// Interval exchange coefficients in original layer order; self is zero.
    pub fn layer_coupling(&self) -> &[f64] {
        &self.layer_coupling
    }
    pub const fn quadrature(&self) -> SlabDepthReport {
        self.report
    }

    /// Signed net heat INTO the interval [W/m² of slab area]. Divide by positive
    /// interval thickness for cell-average volumetric heating [W/m³].
    ///
    /// Irradiances and original layer blackbody exitances use the same gray,
    /// band or wavelength scope. Finite signed deviations from ONE common
    /// background are also admitted. Each original layer retains its uniform
    /// temperature; this does not invent independent temperatures inside it.
    /// Applications cost O(N), without allocation, quadrature or a ray solve.
    pub fn net_flux(
        &self,
        irradiation_w_m2: [f64; 2],
        blackbody_exitance_w_m2: &[f64],
    ) -> Result<f64, ConductionError> {
        if blackbody_exitance_w_m2.len() != self.layers.len()
            || irradiation_w_m2
                .iter()
                .chain(blackbody_exitance_w_m2)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "net depth heating needs finite irradiation and one exitance per original layer",
            ));
        }
        let local = blackbody_exitance_w_m2[self.layer];
        let mut heat = 0.0;
        let mut correction = 0.0;
        for side in 0..2 {
            add(
                &mut heat,
                &mut correction,
                difference_product(self.face_coupling[side], irradiation_w_m2[side], local)?,
            );
        }
        for (&coefficient, &source) in self.layer_coupling.iter().zip(blackbody_exitance_w_m2) {
            add(
                &mut heat,
                &mut correction,
                difference_product(coefficient, source, local)?,
            );
        }
        finite(heat)
    }
}

impl AbsorbingSlab {
    /// Prepare net interval heating from diffuse exterior light AND thermal
    /// emission of all original layers. Reuses the original positive reflection
    /// elimination, without splitting or retabulating any physical material.
    pub fn depth_response(
        &self,
        layer: usize,
        interval_m: [f64; 2],
        controls: SlabQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SlabDepthResponse, ConductionError> {
        checkpoint(&mut cancelled, 0)?;
        let material = self
            .layers
            .get(layer)
            .ok_or_else(|| invalid("net depth heating needs an original layer index"))?;
        if !controls.is_valid()
            || interval_m.iter().any(|x| !nonnegative(*x))
            || interval_m[0] > interval_m[1]
            || interval_m[1] > material.thickness_m()
        {
            return Err(invalid(
                "net depth heating needs an ordered interval within the layer and bounded quadrature controls",
            ));
        }
        let channels = 2 + self.layers.len();
        let mut total = vec![0.0; channels];
        let mut correction = vec![0.0; channels];
        let mut errors = vec![0.0; channels];
        let mut coarse = vec![0.0; channels];
        let mut fine = vec![0.0; channels];
        let mut rule_correction = vec![0.0; channels];
        let mut report = SlabDepthReport {
            controls,
            angular_evaluations: 0,
            accepted_panels: 0,
            maximum_estimated_coefficient_error: 0.0,
        };
        let mut edges = vec![0.0, 1.0];
        edges.extend(self.layers.iter().map(AbsorbingLayer::refractive_index));
        edges.sort_by(f64::total_cmp);
        edges.dedup();
        // Modes above the target's index cannot propagate into this interval.
        let bands: Vec<_> = edges
            .windows(2)
            .take_while(|edge| edge[0] < material.refractive_index())
            .map(|edge| [edge[0], edge[1]])
            .collect();
        let mut pending = Vec::new();
        if material.absorption_m_inv() > 0.0 && interval_m[0] < interval_m[1] {
            for band in (0..bands.len()).rev() {
                pending.push((0.0, 1.0, 0_usize, band));
            }
        }
        while let Some((lo, hi, depth, band)) = pending.pop() {
            checkpoint(&mut cancelled, report.angular_evaluations)?;
            if report.angular_evaluations + 24 > controls.maximum_evaluations {
                return Err(invalid("net depth angular evaluation budget exhausted"));
            }
            for (nodes, output) in [(&GL8[..], &mut coarse), (&GL16[..], &mut fine)] {
                layered::rule(
                    &self.layers,
                    bands[band],
                    lo,
                    hi,
                    nodes,
                    output,
                    &mut rule_correction,
                    &mut cancelled,
                    &mut report.angular_evaluations,
                    Some((layer, interval_m)),
                )?;
            }
            if coarse.iter().zip(&fine).all(|(&a, &b)| {
                (b - a).abs()
                    <= controls.absolute_coefficient_error * (hi - lo) / bands.len() as f64
                        + controls.relative_coefficient_error * b.abs()
            }) {
                for i in 0..channels {
                    add(&mut total[i], &mut correction[i], fine[i]);
                    errors[i] += (fine[i] - coarse[i]).abs();
                }
                report.accepted_panels += 1;
            } else {
                if depth >= controls.maximum_depth {
                    return Err(invalid("net depth angular subdivision depth exhausted"));
                }
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    return Err(invalid("net depth angular panel is unrepresentable"));
                }
                pending.push((mid, hi, depth + 1, band));
                pending.push((lo, mid, depth + 1, band));
            }
        }
        if total.iter().any(|&x| !nonnegative(x)) {
            return Err(invalid("unrepresentable net depth coupling"));
        }
        report.maximum_estimated_coefficient_error = errors.into_iter().fold(0.0_f64, f64::max);
        checkpoint(&mut cancelled, report.angular_evaluations)?;
        Ok(SlabDepthResponse {
            layers: self.layers.clone(),
            layer,
            interval_m,
            face_coupling: [total[0], total[1]],
            layer_coupling: total[2..].into(),
            report,
        })
    }
}
