//! Bounded torsional observations and exact resultants for concentric layers.
//!
//! Point recovery delegates to the original [`super::torsion::TorsionalStressPoint`]
//! law. An explicit outer-radius engineering-shear limit applies to the whole
//! selected layer, so querying near the axis cannot hide surface overstrain.
//! The exact layer peak uses that canonical law at the analytical outer-radius
//! extremum; torque and storage use the SAME cached source G and polar moment
//! as native circular torsion. No averaged G, spatial search, quadrature or
//! transverse-shear correction enters these observations.
//!
//! This is the torsional contribution in the material frame, not the full
//! three-dimensional stress of a bent/sheared rod. End effects, warping of
//! noncircular sections, interlayer slip, plasticity and failure are absent.
//! Large frame rotations are compatible with small local elastic shear;
//! `maximum_engineering_shear` is the caller's explicit applicability bound.
//! Reference: MIT OpenCourseWare, 16.20 Structural Mechanics, torsion of
//! circular shafts, Unit 10:
//! <https://ocw.mit.edu/courses/16-20-structural-mechanics-fall-2002/e888339011c834eb2a02bd365fa25052_unit10.pdf>.

use super::{ConcentricShearSection, MaterialStatePointError, invalid};
use super::torsion::TorsionalStressPoint;

/// Exact radial extremum and integrals of the torsional contribution in a layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularTorsionLayer {
    /// Maximum |tau| at this layer's outer radius [Pa], including its own G.
    pub maximum_shear_pa: f64,
    /// Signed torque carried by this layer [N m].
    pub torque_n_m: f64,
    /// Reversible energy per original axial length [J/m].
    pub energy_per_length_j_m: f64,
}

impl ConcentricShearSection {
    /// Bound the original torsional point law by an explicit small-shear domain.
    ///
    /// Both sides of a shared radius are admitted with their OWN layer index:
    /// bonded layers share twist, while different G values give different
    /// axial-cut shear stresses. The cylindrical interface traction is zero.
    /// The bore and exterior are not material. The strain bound is checked at
    /// the layer's outer radius, so a centerline query cannot hide overstrain.
    /// The returned point is the SAME `TorsionalStressPoint` as the unbounded
    /// observation; stress/strain conventions and representability refusals are
    /// delegated unchanged. The original source lookup remains allocation-free.
    pub fn torsional_stress_at_twist_bounded(
        &self,
        layer: usize,
        coordinates_yz_m: [f64; 2],
        twist_per_m: f64,
        maximum_engineering_shear: f64,
    ) -> Result<TorsionalStressPoint, MaterialStatePointError> {
        self.admit_torsional_layer(layer, twist_per_m, maximum_engineering_shear)?;
        self.torsional_stress_at_twist(layer, coordinates_yz_m, twist_per_m)
    }

    /// Recover one layer's exact outer-radius peak through the original point
    /// law and analytical torque/energy from its original J=2I and cached G.
    /// Work is independent of the layer count; no search, quadrature or solve runs.
    pub fn torsional_layer_response(
        &self,
        layer: usize,
        twist_per_m: f64,
        maximum_engineering_shear: f64,
    ) -> Result<CircularTorsionLayer, MaterialStatePointError> {
        let modulus = self.admit_torsional_layer(layer, twist_per_m, maximum_engineering_shear)?;
        let radius = self.elastic.layers[layer].outer_radius_m;
        // Match the admitted native GJ multiplication order. Forming 2I first
        // could overflow even when the original GJ is representable.
        let stiffness = 2.0 * (modulus * self.elastic.layer_integrals[layer].second_moment_m4());
        let maximum_shear_pa = self
            .torsional_stress_at_twist(layer, [radius, 0.0], twist_per_m)?
            .shear_stress_pa[1]
            .abs();
        let torque_n_m = stiffness * twist_per_m;
        let energy_per_length_j_m = 0.5 * torque_n_m * twist_per_m;
        if !maximum_shear_pa.is_finite()
            || !torque_n_m.is_finite()
            || !energy_per_length_j_m.is_finite()
            || (twist_per_m != 0.0
                && (maximum_shear_pa == 0.0 || torque_n_m == 0.0 || energy_per_length_j_m == 0.0))
        {
            return Err(invalid(
                "layer torsional stress or resultants are unrepresentable",
            ));
        }
        Ok(CircularTorsionLayer {
            maximum_shear_pa,
            torque_n_m,
            energy_per_length_j_m,
        })
    }

    fn admit_torsional_layer(
        &self,
        layer: usize,
        twist_per_m: f64,
        maximum_engineering_shear: f64,
    ) -> Result<f64, MaterialStatePointError> {
        let geometry = self
            .elastic
            .layers
            .get(layer)
            .ok_or_else(|| invalid("torsional stress layer is out of range"))?;
        if !twist_per_m.is_finite()
            || !maximum_engineering_shear.is_finite()
            || maximum_engineering_shear <= 0.0
            || maximum_engineering_shear >= 1.0
            || twist_per_m.abs() * geometry.outer_radius_m > maximum_engineering_shear
        {
            return Err(invalid(
                "torsion requires finite twist within an explicit small-shear domain",
            ));
        }
        self.layer_shear_modulus_pa(layer)
    }
}
