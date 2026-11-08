//! Original, immutable tensor conductivity laws and a declared-frame adapter.
//!
//! This layer owns neither spatial discretization nor physical time/history.
//! A law supplies both endpoint tensors and its actual interval integral.
//! Point sampling at an averaged temperature is not an interval integral.

use super::super::ConductivityTensor;
use crate::conjugate::{ChtError, SolidMaterial};

mod source;
pub use source::RotatedSourceConductivity;

/// Original conductivity in Cartesian grid coordinates, W/(m K).
///
/// Implementations must be deterministic for one solve and must not advance
/// material history, refit a source, or substitute another conductivity law.
/// Temperatures are positive finite Kelvin values. A source-backed law must
/// refuse unsupported endpoints **and every unsupported part of an interval**,
/// including interior gaps; two supported endpoints do not establish support.
///
/// The nonlinear spatial owner uses the fundamental theorem of calculus:
/// derivatives of an interval integral are the original endpoint tensors.
/// It therefore never adds a second conductivity-temperature derivative.
/// Neither this contract nor positive point tensors establishes a discrete
/// maximum principle, entropy inequality, or global nonlinear convergence.
pub trait TensorConductivityLaw {
    /// Original symmetric positive-definite K(T).
    ///
    /// # Errors
    /// Invalid temperature, unavailable original source support, or a tensor
    /// that cannot be represented as finite symmetric positive definite.
    fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError>;

    /// Actual interval average, integral_from^to K(T) dT / (to - from).
    ///
    /// The result is symmetric in the two arguments. At equal endpoints it
    /// is K(T). Implementations must integrate every crossed source segment;
    /// evaluating K at a mean temperature is generally incorrect. Returning
    /// the average avoids subtracting large absolute Kirchhoff potentials
    /// and permits geometry scaling before multiplication by a temperature
    /// difference. This is an integral representation, not a frozen law.
    ///
    /// # Errors
    /// Invalid temperature, unsupported interval, or an unrepresentable
    /// positive-definite average. Never extrapolate without source authority.
    fn mean_conductivity(&self, from_k: f64, to_k: f64) -> Result<ConductivityTensor, ChtError>;
}

impl TensorConductivityLaw for ConductivityTensor {
    fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError> {
        temperature(temperature_k)?;
        Ok(*self)
    }

    fn mean_conductivity(&self, from_k: f64, to_k: f64) -> Result<ConductivityTensor, ChtError> {
        temperature(from_k)?;
        temperature(to_k)?;
        Ok(*self)
    }
}

/// An explicitly oriented native solid with its original directional laws.
///
/// The columns of `rotation` are the principal material directions expressed
/// in Cartesian grid coordinates: K(T) = R diag(k_0(T),k_1(T),k_2(T)) R^T.
/// No orientation is inferred from labels, geometry or scalar conductivity.
///
/// Constant, common-factor tabulated and independent directional tabulated
/// solids retain the native `SolidMaterial` semantics. In particular, native
/// tables explicitly declare constant extensions beyond their end knots.
/// Independent directional tables are absolute laws, without a second scalar
/// multiplier. Every crossed knot is integrated by the native curve owner.
/// Borrowing the material prevents mutation through this adapter during a solve.
#[derive(Debug, Clone, Copy)]
pub struct RotatedSolidConductivity<'a> {
    material: &'a SolidMaterial,
    rotation: [[f64; 3]; 3],
}

impl<'a> RotatedSolidConductivity<'a> {
    /// Bind the original solid to a declared proper orthogonal frame.
    ///
    /// # Errors
    /// Invalid native tables or conductivity references, or a non-finite,
    /// non-orthogonal or reflected frame. The frame is checked, never fitted,
    /// renormalized or projected.
    pub fn new(material: &'a SolidMaterial, rotation: [[f64; 3]; 3]) -> Result<Self, ChtError> {
        material.validate_conductivity_tables()?;
        if !material.conductivity_w_m_k.is_finite()
            || material.conductivity_w_m_k <= 0.0
            || material
                .axis_conductivity()
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "native conductivity references must be positive and finite",
            ));
        }
        if !proper_rotation(rotation) {
            return Err(invalid(
                "principal axes need a finite proper orthogonal rotation",
            ));
        }
        Ok(Self { material, rotation })
    }

    fn rotated(&self, principal: [f64; 3]) -> Result<ConductivityTensor, ChtError> {
        let matrix = rotated_diagonal(principal, self.rotation)
            .ok_or_else(|| invalid("directional conductivity is not positive and finite"))?;
        ConductivityTensor::new(matrix)
    }
}

impl TensorConductivityLaw for RotatedSolidConductivity<'_> {
    fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError> {
        temperature(temperature_k)?;
        self.rotated(std::array::from_fn(|axis| {
            self.material.conductivity_curve(axis).at(temperature_k)
        }))
    }

    fn mean_conductivity(&self, from_k: f64, to_k: f64) -> Result<ConductivityTensor, ChtError> {
        temperature(from_k)?;
        temperature(to_k)?;
        self.rotated(std::array::from_fn(|axis| {
            self.material.conductivity_curve(axis).mean(from_k, to_k)
        }))
    }
}

fn temperature(value: f64) -> Result<(), ChtError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "tensor-law temperature must be finite and positive Kelvin",
        ))
    }
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "thermal.tensor_conduction.material",
        reason: reason.into(),
    }
}

fn proper_rotation(rotation: [[f64; 3]; 3]) -> bool {
    if rotation.iter().flatten().any(|v| !v.is_finite()) {
        return false;
    }
    // A declared floating-point frame gets a small orthogonality allowance;
    // this is not an eigenvalue floor or an adjustment to authored axes.
    const TOLERANCE: f64 = 1.0e-12;
    for a in 0..3 {
        for b in 0..=a {
            let dot = (0..3).fold(0.0, |sum, i| rotation[i][a].mul_add(rotation[i][b], sum));
            let expected = if a == b { 1.0 } else { 0.0 };
            if (dot - expected).abs() > TOLERANCE {
                return false;
            }
        }
    }
    let determinant = rotation[0][0]
        * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]);
    (determinant - 1.0).abs() <= 3.0 * TOLERANCE
}

fn rotated_diagonal(principal: [f64; 3], rotation: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    if principal.iter().any(|k| !k.is_finite() || *k <= 0.0) {
        return None;
    }
    // A truly isotropic native law has no frame dependence, even to roundoff.
    if principal[0] == principal[1] && principal[1] == principal[2] {
        return Some([
            [principal[0], 0.0, 0.0],
            [0.0, principal[0], 0.0],
            [0.0, 0.0, principal[0]],
        ]);
    }
    let mut tensor = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let value = (0..3).fold(0.0, |sum, a| {
                (rotation[i][a] * principal[a]).mul_add(rotation[j][a], sum)
            });
            if !value.is_finite() {
                return None;
            }
            tensor[i][j] = value;
            tensor[j][i] = value;
        }
    }
    Some(tensor)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];

    fn close(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) {
        for i in 0..3 {
            for j in 0..3 {
                assert!((a[i][j] - b[i][j]).abs() < 2e-13 * b[i][j].abs().max(1.0));
                assert_eq!(a[i][j], a[j][i]);
            }
        }
    }

    #[test]
    fn native_independent_tables_keep_all_knots_and_declared_rotation() {
        let material = SolidMaterial::new("independent", 99.0)
            .with_orthotropic([71.0, 83.0, 97.0])
            .with_axis_conductivity_tables([
                &[(300.0, 2.0), (350.0, 8.0), (400.0, 4.0)],
                &[(250.0, 1.0), (375.0, 6.0), (500.0, 3.0)],
                &[(300.0, 7.0)],
            ]);
        let law = RotatedSolidConductivity::new(&material, ROTATION).unwrap();
        // Independent hand integrals over 200..550 include both constant tails
        // and all differing axis knots: x=1350, y=1200, z=2450 W/m.
        let expected = rotated_diagonal([1350.0 / 350.0, 1200.0 / 350.0, 7.0], ROTATION).unwrap();
        close(
            law.mean_conductivity(200.0, 550.0).unwrap().as_array(),
            expected,
        );
        close(
            law.mean_conductivity(550.0, 200.0).unwrap().as_array(),
            expected,
        );
        let point = law.conductivity_at_temperature(350.0).unwrap();
        close(
            law.mean_conductivity(350.0, 350.0).unwrap().as_array(),
            point.as_array(),
        );
        assert_ne!(
            law.conductivity_at_temperature(375.0).unwrap().as_array(),
            expected,
        );
    }

    #[test]
    fn common_factor_has_exact_kirchhoff_integral_without_double_scaling() {
        let material = SolidMaterial::new("common", 2.0)
            .with_orthotropic([6.0, 2.0, 10.0])
            .with_conductivity_table(&[(300.0, 2.0), (400.0, 6.0), (500.0, 2.0)]);
        let law = RotatedSolidConductivity::new(&material, ROTATION).unwrap();
        close(
            law.mean_conductivity(300.0, 500.0).unwrap().as_array(),
            rotated_diagonal([12.0, 4.0, 20.0], ROTATION).unwrap(),
        );
        close(
            law.conductivity_at_temperature(400.0).unwrap().as_array(),
            rotated_diagonal([18.0, 6.0, 30.0], ROTATION).unwrap(),
        );
    }

    #[test]
    fn isotropic_law_is_exactly_frame_invariant_and_constant_tensor_is_a_law() {
        let material = SolidMaterial::new("isotropic", 3.25);
        let law = RotatedSolidConductivity::new(&material, ROTATION).unwrap();
        let expected = [[3.25, 0.0, 0.0], [0.0, 3.25, 0.0], [0.0, 0.0, 3.25]];
        assert_eq!(
            law.mean_conductivity(100.0, 900.0).unwrap().as_array(),
            expected
        );
        let constant = ConductivityTensor::new(expected).unwrap();
        assert_eq!(constant.mean_conductivity(100.0, 900.0).unwrap(), constant);
        assert_eq!(
            constant.conductivity_at_temperature(450.0).unwrap(),
            constant
        );
    }

    #[test]
    fn refuses_bad_axes_bad_native_sources_and_nonphysical_temperatures() {
        let material = SolidMaterial::new("solid", 1.0);
        for frame in [
            [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ] {
            assert!(RotatedSolidConductivity::new(&material, frame).is_err());
        }
        let invalid = material
            .clone()
            .with_conductivity_table(&[(300.0, 1.0), (300.0, 2.0)]);
        assert!(RotatedSolidConductivity::new(&invalid, ROTATION).is_err());
        let competing = material
            .clone()
            .with_conductivity_table(&[(300.0, 1.0)])
            .with_axis_conductivity_tables([&[(300.0, 1.0)]; 3]);
        assert!(RotatedSolidConductivity::new(&competing, ROTATION).is_err());
        let law = RotatedSolidConductivity::new(&material, ROTATION).unwrap();
        for t in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(law.conductivity_at_temperature(t).is_err());
            assert!(law.mean_conductivity(300.0, t).is_err());
            assert!(law.mean_conductivity(t, 300.0).is_err());
        }
    }
}
