//! A narrow, free-sided Euler--Bernoulli strip from the original thermal laminate.
//!
//! Each original layer relaxes `sigma_yy = tau_xy = 0` locally. Its axial
//! modulus is the Schur complement of its retained plane-stress stiffness;
//! axial free strain remains that layer's original integrated thermal strain.
//! The beam convention is `epsilon_x(z) = eps0 - (z-axis_offset)*kappa`, so
//! the conjugate moment is `-width * integral((z-axis_offset)*sigma_x dz)`.
//!
//! This is an explicit one-dimensional narrow-strip approximation, not a
//! compatible finite-width plate solution. Lateral contraction can differ
//! between layers; edge boundary layers, transverse/twist inertia, warping,
//! interlaminar tractions and torsion are omitted. Significant layerwise
//! extension/shear coupling or thermal shear refuses instead of silently
//! entering a planar frame. Temperature, heat/work accounting and finite-strain
//! admission stay with the existing source and structural owners.

use super::{LaminateThermalError, LayerThermalStress, PreparedLaminateThermal};
use fs_math::det;
use std::sync::Arc;

/// Reversible source-bound potential with the same free-sided relaxation.
pub mod thermoelastic;

#[derive(Clone, Copy, Debug, PartialEq)]
struct LayerReduction {
    axial_pa: f64,
    // Elastic [epsilon_y, gamma_xy] = -relaxation * elastic epsilon_x.
    relaxation: [f64; 2],
}

/// Exact section resultants and positive elastic storage for the declared
/// narrow-strip closure, per unit original beam length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateBeamResponse {
    pub axial_force_n: f64,
    pub bending_moment_n_m: f64,
    /// Includes incompatible layer self-stress even when N and M both vanish.
    pub strain_energy_j_m: f64,
}

/// Prepared reduction of the same immutable laminate, source profiles and mass.
/// Layer preparation is O(layers). Cached resultants/tangent allocate nothing;
/// energy uses original Gauss samples and stress recovery uses original faces.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedLaminateBeam {
    source: Arc<PreparedLaminateThermal>,
    width_m: f64,
    axis_offset_m: f64,
    layers: Vec<LayerReduction>,
    tangent: [[f64; 2]; 2],
    thermal: [f64; 2],
    free: [f64; 2],
    mass_moments: [f64; 3],
}

impl PreparedLaminateBeam {
    /// Select the free-sided narrow-strip closure explicitly. `axis_offset_m`
    /// is measured from the original laminate reference toward its top face.
    /// Choosing I1/I0 places the beam axis at the original density centroid.
    /// Neither geometry nor density is refreshed from the current temperature.
    pub fn try_new(
        source: Arc<PreparedLaminateThermal>,
        width_m: f64,
        axis_offset_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateThermalError> {
        poll(&mut cancelled)?;
        if !width_m.is_finite() || width_m <= 0.0 || !axis_offset_m.is_finite() {
            return Err(invalid(
                "positive finite strip width and finite reference axis required",
            ));
        }
        let mut layers = Vec::with_capacity(source.samples.len());
        let mut tangent = [[0.0; 2]; 2];
        let mut thermal = [0.0; 2];
        let mut mass_moments = [width_m * source.section().i0_kg_m2(), 0.0, 0.0];
        for (index, sample) in source.samples.iter().enumerate() {
            poll(&mut cancelled)?;
            let q = sample.q;
            let shear_scale = det::sqrt(q[8]);
            for normal in [0, 1] {
                let scale = det::sqrt(q[3 * normal + normal]) * shear_scale;
                if !scale.is_finite() || q[3 * normal + 2].abs() > 128.0 * f64::EPSILON * scale {
                    return Err(invalid(
                        "planar strip requires each layer's axial/shear alignment",
                    ));
                }
            }
            for free in sample.face_free {
                let scale = free[0].abs().max(free[1].abs());
                if free[2].abs() > 128.0 * f64::EPSILON * scale {
                    return Err(invalid("planar strip cannot admit a layer's thermal shear"));
                }
            }
            let relaxation = solve_spd2([[q[4], q[5]], [q[7], q[8]]], [q[3], q[6]])?;
            let axial_pa = (-q[2]).mul_add(relaxation[1], (-q[1]).mul_add(relaxation[0], q[0]));
            if !axial_pa.is_finite() || axial_pa <= 0.0 {
                return Err(invalid(
                    "positive representable layer axial Schur complement required",
                ));
            }
            let reduction = LayerReduction {
                axial_pa,
                relaxation,
            };
            let stiffness_weight = width_m * sample.weight_m * axial_pa;
            for p in 0..2 {
                let z = sample.z[p] - axis_offset_m;
                tangent[0][0] += stiffness_weight;
                tangent[0][1] -= stiffness_weight * z;
                tangent[1][1] += (stiffness_weight * z) * z;
                thermal[0] += stiffness_weight * sample.free[p][0];
                thermal[1] -= (stiffness_weight * sample.free[p][0]) * z;
            }
            // Integrate original mass about the requested axis without the
            // cancellation in I2 - 2*offset*I1 + offset^2*I0.
            let integral = &source.section().layer_integrals()[index];
            let layer = &source.section().layers()[index];
            let h = integral.z_top_m - integral.z_bottom_m;
            let center = 0.5 * integral.z_bottom_m + 0.5 * integral.z_top_m - axis_offset_m;
            let mass = width_m * layer.density_kg_m3 * h;
            mass_moments[1] += mass * center;
            mass_moments[2] += (mass * center) * center + (mass * h) * (h / 12.0);
            layers.push(reduction);
        }
        tangent[1][0] = tangent[0][1];
        finite(
            tangent
                .iter()
                .flatten()
                .chain(&thermal)
                .chain(&mass_moments),
        )?;
        if mass_moments[0] <= 0.0 || mass_moments[2] <= 0.0 {
            return Err(invalid(
                "positive representable reference strip mass and rotary inertia required",
            ));
        }
        let free = solve_spd2(tangent, thermal)?;
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            width_m,
            axis_offset_m,
            layers,
            tangent,
            thermal,
            free,
            mass_moments,
        })
    }

    /// Original thermal section, including every layer law, frame and receipt.
    pub fn source(&self) -> &PreparedLaminateThermal {
        &self.source
    }
    pub const fn width_m(&self) -> f64 {
        self.width_m
    }
    pub const fn axis_offset_m(&self) -> f64 {
        self.axis_offset_m
    }
    /// d[N,M]/d[eps0,kappa], with units [[N,N m],[N m,N m²]].
    pub const fn tangent(&self) -> [[f64; 2]; 2] {
        self.tangent
    }
    /// Equivalent thermal loading t in [N,M] = K*[eps0,kappa] - t.
    pub const fn thermal_resultants(&self) -> [f64; 2] {
        self.thermal
    }
    /// Generalized strain at zero net axial force and bending moment.
    /// This need not remove incompatible layer stress or elastic storage.
    pub const fn free_generalized_strain(&self) -> [f64; 2] {
        self.free
    }
    /// [integral rho dA, integral rho*z_axis dA, integral rho*z_axis² dA],
    /// in [kg/m, kg, kg m]. A nonzero first moment remains explicit.
    pub const fn reference_mass_moments(&self) -> [f64; 3] {
        self.mass_moments
    }

    /// Evaluate the original layer samples after local lateral relaxation.
    /// Squaring each elastic strain retains positive self-stress energy;
    /// subtracting aggregate thermal work from total storage would lose it.
    pub fn evaluate(
        &self,
        eps0: f64,
        kappa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateBeamResponse, LaminateThermalError> {
        finite([eps0, kappa].iter())?;
        let mut result = LaminateBeamResponse {
            axial_force_n: 0.0,
            bending_moment_n_m: 0.0,
            strain_energy_j_m: 0.0,
        };
        for (sample, reduced) in self.source.samples.iter().zip(&self.layers) {
            poll(&mut cancelled)?;
            let weight = self.width_m * sample.weight_m;
            for p in 0..2 {
                let z = sample.z[p] - self.axis_offset_m;
                let elastic = (-z).mul_add(kappa, eps0) - sample.free[p][0];
                let stress = reduced.axial_pa * elastic;
                let force = weight * stress;
                result.axial_force_n += force;
                result.bending_moment_n_m -= z * force;
                result.strain_energy_j_m += 0.5 * force * elastic;
            }
        }
        finite(
            [
                result.axial_force_n,
                result.bending_moment_n_m,
                result.strain_energy_j_m,
            ]
            .iter(),
        )?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Total engineering [epsilon_x, epsilon_y, gamma_xy] at this layer's
    /// original faces. Transverse contraction is recovered separately per layer;
    /// a common affine laminate strain cannot represent this strip closure.
    pub fn layer_strain_faces(
        &self,
        layer: usize,
        eps0: f64,
        kappa: f64,
    ) -> Result<[[f64; 3]; 2], LaminateThermalError> {
        finite([eps0, kappa].iter())?;
        let reduced = self
            .layers
            .get(layer)
            .ok_or(LaminateThermalError::LayerIndex {
                index: layer,
                count: self.layers.len(),
            })?;
        let sample = &self.source.samples[layer];
        let integral = &self.source.section().layer_integrals()[layer];
        let positions = [integral.z_bottom_m, integral.z_top_m];
        let strains = core::array::from_fn(|p| {
            let total = (self.axis_offset_m - positions[p]).mul_add(kappa, eps0);
            let elastic = total - sample.face_free[p][0];
            [
                total,
                (-reduced.relaxation[0]).mul_add(elastic, sample.face_free[p][1]),
                (-reduced.relaxation[1]).mul_add(elastic, sample.face_free[p][2]),
            ]
        });
        finite(strains.iter().flatten())?;
        Ok(strains)
    }

    /// Actual relaxed in-plane face stresses. The admitted closure sets the
    /// transverse and shear stresses to zero, not their total strains.
    pub fn layer_stress(
        &self,
        layer: usize,
        eps0: f64,
        kappa: f64,
    ) -> Result<LayerThermalStress, LaminateThermalError> {
        let strains = self.layer_strain_faces(layer, eps0, kappa)?;
        let sample = &self.source.samples[layer];
        let axial = self.layers[layer].axial_pa;
        let stress: [[f64; 3]; 2] =
            core::array::from_fn(|p| [axial * (strains[p][0] - sample.face_free[p][0]), 0.0, 0.0]);
        finite(stress.iter().flatten())?;
        Ok(LayerThermalStress {
            bottom_pa: stress[0],
            top_pa: stress[1],
        })
    }
}

// Diagonal scaling admits stiffnesses with unlike physical units without a
// dimensional pivot threshold. The positive normalized determinant is checked;
// singular or unrepresentable reductions refuse rather than being regularized.
fn solve_spd2(a: [[f64; 2]; 2], rhs: [f64; 2]) -> Result<[f64; 2], LaminateThermalError> {
    finite(a.iter().flatten().chain(&rhs))?;
    if a[0][0] <= 0.0 || a[1][1] <= 0.0 {
        return Err(invalid("positive section stiffness required"));
    }
    let sx = det::sqrt(a[0][0]);
    let sy = det::sqrt(a[1][1]);
    let r = (a[0][1] / sx) / sy;
    let d = (1.0 - r) * (1.0 + r);
    if !d.is_finite() || d <= 0.0 {
        return Err(invalid(
            "positive representable axial/bending compliance required",
        ));
    }
    let u = rhs[0] / sx;
    let v = (-r).mul_add(u, rhs[1] / sy) / d;
    let solution = [(-r).mul_add(v, u) / sx, v / sy];
    finite(solution.iter())?;
    Ok(solution)
}

fn invalid(quantity: &'static str) -> LaminateThermalError {
    LaminateThermalError::Nonfinite { quantity }
}
fn finite<'a>(values: impl Iterator<Item = &'a f64>) -> Result<(), LaminateThermalError> {
    if values.into_iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(invalid("representable laminate beam arithmetic required"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LaminateThermalError> {
    if cancelled() {
        Err(LaminateThermalError::Cancelled)
    } else {
        Ok(())
    }
}
