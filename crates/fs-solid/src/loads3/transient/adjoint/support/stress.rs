//! Elastic stress and smooth von-Mises strength objectives at accepted states.
//! All six-component vectors use mesh-frame Mandel [xx,yy,zz,xy,yz,zx]: the
//! final three stress entries are sqrt(2) times the physical tensor shears.
use super::super::material::{MaterialGradientBudget3, MaterialScaleGradient3};
use super::super::{TransientAdjointError3, check, number};
use super::NodalMotionCotangent3;
use crate::linear3::TetThermalStrainField;
use crate::loads3::transient::{TetTransientRun3, TransientError3, response};
use fs_exec::Cx;

#[derive(Debug, Clone, PartialEq)]
pub struct StressPullback3 {
    /// Full displacement seed, including fixed DOFs; velocity and acceleration
    /// entries are zero. Chain this through the matching endpoint's reverse.
    pub motion: NodalMotionCotangent3,
    /// Direct stiffness-scale derivative at fixed q and free thermal strain.
    /// Direct density/Rayleigh entries are zero. Add implicit trajectory terms
    /// for total design derivatives; direct thermal-stress changes are included.
    pub material: MaterialScaleGradient3,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Reverse pointwise element stresses sigma=C(Bq-epsilon_free) at this
    /// exact accepted endpoint, using the retained material and thermal field.
    /// Rows are sorted unique element IDs and finite dJ/dsigma [objective/Pa]
    /// in mesh-frame Mandel [xx,yy,zz,xy,yz,zx]. Omitted rows have zero seed.
    /// No volume or time-quadrature factor is inserted: callers include those
    /// weights in their objective seeds. Supports may be moving or stationary.
    ///
    /// Geometry, material frame and free thermal strain stay fixed. Both total
    /// and elastic Mandel strains must meet the finite positive ceiling, as in
    /// recover_response. This is a small-strain elastic observation derivative,
    /// not a yield/plasticity law or derivative of its admission decisions.
    /// Budget refusal/cancellation returns no partial gradient and changes no
    /// accepted state. Direct material terms must accompany the motion reverse.
    pub fn stress_pullback(
        &self,
        rows: &[(usize, [f64; 6])],
        maximum_strain: f64,
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<StressPullback3, TransientAdjointError3> {
        self.admit_stress_pullback(rows, maximum_strain, budget, cx)?;
        let motion = self.nodal_motion(cx)?;
        self.stress_pullback_nodal_at(&motion[0], rows, maximum_strain, budget, cx)
    }

    /// The recorder owns association of this exact q with the immutable run's
    /// constitutive and thermal inputs. All rows use mesh-frame Mandel stress.
    pub(in crate::loads3::transient::adjoint) fn stress_pullback_nodal_at(
        &self,
        displacement: &[[f64; 3]],
        rows: &[(usize, [f64; 6])],
        maximum_strain: f64,
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<StressPullback3, TransientAdjointError3> {
        self.admit_stress_pullback(rows, maximum_strain, budget, cx)?;
        let thermal = self.thermal_strains.as_ref().map(|field| match field {
            response::FrozenThermal3::Uniform(v) => TetThermalStrainField::Uniform(v),
            response::FrozenThermal3::PerElement(v) => TetThermalStrainField::PerElement(v),
        });
        let (q, stiffness) = self
            .problem
            .stress_observation_contraction(displacement, rows, thermal, maximum_strain, cx)
            .map_err(TransientError3::from)?;
        let nodes = self.problem.nodes_m.len();
        let result = StressPullback3 {
            motion: NodalMotionCotangent3 {
                q,
                v: vec![[0.0; 3]; nodes],
                a: vec![[0.0; 3]; nodes],
            },
            material: MaterialScaleGradient3 {
                element_scales: stiffness.into_iter().map(|s| [s, 0.0]).collect(),
                rayleigh: [0.0; 2],
            },
        };
        check(cx)?;
        Ok(result)
    }

    fn admit_stress_pullback(
        &self,
        rows: &[(usize, [f64; 6])],
        maximum_strain: f64,
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        if !maximum_strain.is_finite() || maximum_strain <= 0.0 {
            return Err(TransientAdjointError3::Invalid(
                "stress needs a finite positive strain ceiling",
            ));
        }
        let full = self
            .problem
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TransientAdjointError3::Budget)?;
        if full > budget.maximum_full_dofs
            || self.problem.tetrahedra.len() > budget.maximum_elements
        {
            return Err(TransientAdjointError3::Budget);
        }
        let mut previous = None;
        for &(element, weights) in rows {
            check(cx)?;
            if element >= self.problem.tetrahedra.len() || previous.is_some_and(|p| p >= element) {
                return Err(TransientAdjointError3::Invalid(
                    "stress rows must name sorted unique existing elements",
                ));
            }
            for value in weights {
                number(value)?;
            }
            previous = Some(element);
        }
        Ok(())
    }
}

/// Smooth dimensionless strength constraint g=(sigma_vm/limit)^2-1 and its
/// d(g)/d(stress) seed in Mandel [xx,yy,zz,xy,yz,zx], including sqrt(2) shears.
/// `g <= 0` is equivalent to the ordinary von-Mises limit for positive limit.
/// The squared formulation has a defined zero derivative at hydrostatic/zero
/// stress. The supplied limit is a caller-owned criterion, not a material-law
/// admission, physical failure probability, or yield authority.
pub fn von_mises_limit(
    stress_mandel_pa: [f64; 6],
    limit_pa: f64,
) -> Result<(f64, [f64; 6]), TransientAdjointError3> {
    if !limit_pa.is_finite() || limit_pa <= 0.0 {
        return Err(TransientAdjointError3::Invalid(
            "von-Mises limit must be finite and positive",
        ));
    }
    let mut s = [0.0; 6];
    for i in 0..6 {
        s[i] = number(number(stress_mandel_pa[i])? / limit_pa)?;
    }
    let xy = number(s[0] - s[1])?;
    let yz = number(s[1] - s[2])?;
    let zx = number(s[2] - s[0])?;
    let normal = number(0.5 * number(xy * xy + yz * yz + zx * zx)?)?;
    let shear = number(1.5 * number(s[3] * s[3] + s[4] * s[4] + s[5] * s[5])?)?;
    let value = number(normal + shear - 1.0)?;
    let mut seed = [0.0; 6];
    seed[0] = number(number(xy - zx)? / limit_pa)?;
    seed[1] = number(number(yz - xy)? / limit_pa)?;
    seed[2] = number(number(zx - yz)? / limit_pa)?;
    for i in 3..6 {
        seed[i] = number(number(3.0 * s[i])? / limit_pa)?;
    }
    Ok((value, seed))
}
