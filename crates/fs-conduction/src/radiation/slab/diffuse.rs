//! Prepared depth-resolved absorption of diffuse external irradiation.

use super::{
    AbsorbingLayer, AbsorbingSlab, Arc, ConductionError, GL8, GL16, SlabQuadrature, add,
    checkpoint, finite, invalid, nonnegative,
};

/// Bounded preparation work and embedded angular discrepancies.
/// These are numerical diagnostics, not certified continuum error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabDiffuseDepthReport {
    pub controls: SlabQuadrature,
    /// Counts every ordinate, including rejected panels.
    pub angular_evaluations: usize,
    pub accepted_panels: usize,
    /// Sum of accepted GL8/GL16 discrepancies for [front, back] absorptance.
    pub estimated_absorptance_error: [f64; 2],
}

/// Immutable diffuse absorption in one interval of an original physical layer.
///
/// Uniform unpolarized external radiance has the hemispherical measure 2μ dμ.
/// Each direction retains Snell refraction and every incoherent Fresnel return
/// through the original stack. Depth integration uses positive exact Beer loss;
/// no optical layer is split and no angular redistribution is introduced.
///
/// Repeated heat applications cost O(1), without allocation or quadrature.
/// This is incident heating only; it does not assign a depth profile to thermal
/// emission or interlayer thermal exchange. The caller retains the original
/// gray, band or wavelength scope of the sourced optical material.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabDiffuseDepthResponse {
    layers: Arc<[AbsorbingLayer]>,
    exterior_refractive_indices: [f64; 2],
    exterior_index_sources: [Option<Arc<str>>; 2],
    layer: usize,
    interval_m: [f64; 2],
    absorptance: [f64; 2],
    report: SlabDiffuseDepthReport,
}

impl SlabDiffuseDepthResponse {
    pub fn layers(&self) -> &[AbsorbingLayer] {
        &self.layers
    }

    /// Index at the FRONT face; use the plural accessor for both exterior media.
    pub const fn exterior_refractive_index(&self) -> f64 {
        self.exterior_refractive_indices[0]
    }

    pub const fn exterior_refractive_indices(&self) -> [f64; 2] {
        self.exterior_refractive_indices
    }

    pub fn exterior_index_source(&self) -> Option<&str> {
        self.exterior_index_sources[0].as_deref()
    }

    pub fn exterior_index_sources(&self) -> [Option<&str>; 2] {
        self.exterior_index_sources
            .each_ref()
            .map(|source| source.as_deref())
    }

    pub const fn layer_index(&self) -> usize {
        self.layer
    }

    /// Depths [lo, hi] measured from the front face of the original layer [m].
    pub const fn interval_m(&self) -> [f64; 2] {
        self.interval_m
    }

    /// Absorption fractions of incident diffuse flux from [front, back].
    pub const fn absorptance(&self) -> [f64; 2] {
        self.absorptance
    }

    pub const fn quadrature(&self) -> SlabDiffuseDepthReport {
        self.report
    }

    /// Positive absorbed power per physical slab area [W/m²]. Inputs are the
    /// already hemispherically integrated [front, back] irradiances [W/m²].
    /// Divide by positive interval thickness for cell-average heating [W/m³].
    /// A zero-width interval absorbs zero; no second cosine factor is applied.
    pub fn absorbed_flux(&self, irradiation_w_m2: [f64; 2]) -> Result<f64, ConductionError> {
        if irradiation_w_m2.iter().any(|x| !nonnegative(*x)) {
            return Err(invalid(
                "diffuse depth heating needs nonnegative finite irradiation",
            ));
        }
        finite(
            self.absorptance[0] * irradiation_w_m2[0] + self.absorptance[1] * irradiation_w_m2[1],
        )
    }
}

impl AbsorbingSlab {
    /// Prepare diffuse incident heating for a spatial cell within one original
    /// layer, retaining its source and the complete original reflecting stack.
    /// Preparation costs O(N) per angular ordinate and is bounded separately
    /// from the owner's whole-layer thermal preparation. Cache the result for
    /// repeated heat steps; changing irradiance needs no further angular work.
    pub fn diffuse_depth_response(
        &self,
        layer: usize,
        interval_m: [f64; 2],
        controls: SlabQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SlabDiffuseDepthResponse, ConductionError> {
        checkpoint(&mut cancelled, 0)?;
        let material = self
            .layers
            .get(layer)
            .ok_or_else(|| invalid("diffuse depth heating needs an original layer index"))?;
        if !controls.is_valid()
            || interval_m.iter().any(|x| !nonnegative(*x))
            || interval_m[0] > interval_m[1]
            || interval_m[1] > material.thickness_m()
        {
            return Err(invalid(
                "diffuse depth heating needs an ordered interval within the layer and bounded quadrature controls",
            ));
        }
        let mut report = SlabDiffuseDepthReport {
            controls,
            angular_evaluations: 0,
            accepted_panels: 0,
            estimated_absorptance_error: [0.0; 2],
        };
        let mut total = [0.0; 2];
        let mut correction = [0.0; 2];
        let mut error_correction = [0.0; 2];
        let mut pending = Vec::new();
        let mut bands = Vec::new();
        if material.absorption_m_inv() > 0.0 && interval_m[0] < interval_m[1] {
            // Conserved p=n*sin(theta), split at each accessible critical
            // direction. Sine-squared panel coordinates remove endpoint cusps
            // and preserve positive Snell discriminants through the beam owner.
            let maximum_exterior =
                self.exterior_refractive_indices[0].max(self.exterior_refractive_indices[1]);
            let mut edges = vec![
                0.0,
                self.exterior_refractive_indices[0],
                self.exterior_refractive_indices[1],
            ];
            edges.extend(
                self.layers
                    .iter()
                    .map(AbsorbingLayer::refractive_index)
                    .filter(|&index| index < maximum_exterior),
            );
            edges.sort_by(f64::total_cmp);
            edges.dedup();
            bands.extend(edges.windows(2).map(|edge| [edge[0], edge[1]]));
            for band in (0..bands.len()).rev() {
                pending.push((0.0, 1.0, 0_usize, band));
            }
        }
        while let Some((lo, hi, depth, band)) = pending.pop() {
            checkpoint(&mut cancelled, report.angular_evaluations)?;
            if report.angular_evaluations + 24 > controls.maximum_evaluations {
                return Err(invalid(format!(
                    "diffuse depth angular evaluation budget exhausted: {} of {}",
                    report.angular_evaluations, controls.maximum_evaluations,
                )));
            }
            let [lower, upper] = bands[band];
            let band_mass = self.exterior_refractive_indices.map(|index| {
                if index > lower {
                    ((upper - lower) / index) * ((upper + lower) / index)
                } else {
                    0.0
                }
            });
            let mut estimates = [[0.0; 2]; 2];
            for (nodes, estimate) in [&GL8[..], &GL16[..]].into_iter().zip(&mut estimates) {
                // The beam owner uses the same positive Fresnel elimination and
                // original endpoint fields as collimated depth queries. Interior
                // nodes avoid evaluating the exactly tangent external ray.
                let midpoint = 0.5 * (lo + hi);
                let half_width = 0.5 * (hi - lo);
                let mut rule_correction = [0.0; 2];
                for &(node, weight) in nodes {
                    for sign in [-1.0, 1.0] {
                        checkpoint(&mut cancelled, report.angular_evaluations)?;
                        let ordinate = midpoint + sign * half_width * node;
                        let (sin, cos) =
                            fs_math::det::sin_cos(core::f64::consts::FRAC_PI_2 * ordinate);
                        let beam = self.beam_response_in_band(bands[band], cos)?;
                        let fractions = beam.layer_interval_absorptance(layer, interval_m)?;
                        let angular_weight =
                            core::f64::consts::PI * sin * cos * half_width * weight;
                        for side in 0..2 {
                            add(
                                &mut estimate[side],
                                &mut rule_correction[side],
                                band_mass[side] * angular_weight * fractions[side],
                            );
                        }
                        report.angular_evaluations += 1;
                    }
                }
            }
            let [coarse, fine] = estimates;
            let errors: [f64; 2] = core::array::from_fn(|side| (fine[side] - coarse[side]).abs());
            if (0..2).all(|side| {
                errors[side]
                    <= controls.absolute_coefficient_error * band_mass[side] * (hi - lo)
                        + controls.relative_coefficient_error * fine[side]
            }) {
                for side in 0..2 {
                    add(&mut total[side], &mut correction[side], fine[side]);
                    add(
                        &mut report.estimated_absorptance_error[side],
                        &mut error_correction[side],
                        errors[side],
                    );
                }
                report.accepted_panels += 1;
            } else {
                if depth >= controls.maximum_depth {
                    return Err(invalid("diffuse depth angular subdivision depth exhausted"));
                }
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    return Err(invalid("diffuse depth angular panel is unrepresentable"));
                }
                pending.push((mid, hi, depth + 1, band));
                pending.push((lo, mid, depth + 1, band));
            }
        }
        let allowance = 256.0 * f64::EPSILON * self.layers.len() as f64;
        if total
            .iter()
            .any(|x| !nonnegative(*x) || *x > 1.0 + allowance)
        {
            return Err(invalid("unrepresentable diffuse depth absorptance"));
        }
        checkpoint(&mut cancelled, report.angular_evaluations)?;
        Ok(SlabDiffuseDepthResponse {
            layers: self.layers.clone(),
            exterior_refractive_indices: self.exterior_refractive_indices,
            exterior_index_sources: self.exterior_index_sources.clone(),
            layer,
            interval_m,
            absorptance: total,
            report,
        })
    }
}
