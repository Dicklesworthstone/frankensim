//! Element fields of the original finite-strain P1 tetrahedral problem.
//!
//! First Piola stress P remains conjugate to the TOTAL deformation gradient F
//! on the original reference mesh. Cauchy stress is the full push-forward
//! sigma = P F^T / det(F), so current-area tractions retain the same resultant
//! as reference-area Piola tractions. Material interfaces remain discontinuous;
//! no nodal smoothing or stress symmetrization alters the original response.
//!
//! Total Green--Lagrange strain is geometric, including any thermal or plastic
//! deformation represented by a material adapter. It is not that adapter's
//! elastic strain. Likewise the original material potential equals stored
//! elastic energy only for a hyperelastic law; a frozen inelastic owner may
//! supply an incremental potential and must report storage and heat separately.
//! Recovery evaluates the same combined stress/potential seam as native force
//! assembly, without extracting a bare law or forming a constitutive tangent.
//!
//! The supplied displacement need not be an equilibrium solution. These fields
//! describe its admitted P1 state; they do not certify equilibrium, a continuum
//! peak stress, material failure or spatial convergence.

use fs_exec::Cx;
use fs_material::hyper::cards::FiniteStrainPotential;

use super::{Element, HyperTetError, HyperTetProblem, checkpoint};

/// Constant material and geometric fields on one original P1 tetrahedron.
///
/// All tensors are row-major 3-by-3 arrays. Material stress retains the source
/// law and its original reference configuration, including reference-deformation
/// adapters. Element and node indices bind the fields to the borrowed problem.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperTetElementResponse {
    /// Element index in the original connectivity/material-field order.
    pub element: usize,
    /// Original global node IDs, in original local-node order.
    pub nodes: [usize; 4],
    /// Total deformation gradient F, dimensionless.
    pub deformation_gradient: [f64; 9],
    /// Positive det(F), the current/reference volume ratio.
    pub jacobian: f64,
    /// First Piola stress P = d(potential density)/dF, Pa.
    pub first_piola_pa: [f64; 9],
    /// Current spatial stress P F^T / det(F), Pa.
    ///
    /// Objective material laws give a symmetric tensor to roundoff. Every
    /// original component is retained rather than averaging away asymmetry.
    pub cauchy_pa: [f64; 9],
    /// Total Green--Lagrange strain (F^T F - I)/2, dimensionless.
    ///
    /// This is total geometric strain, not an inferred elastic strain or
    /// an engineering-shear/Mandel vector.
    pub green_lagrange_strain: [f64; 9],
    /// Original reference volume, m^3.
    pub reference_volume_m3: f64,
    /// Current volume det(F) times original reference volume, m^3.
    pub current_volume_m3: f64,
    /// Original material potential per reference volume, J/m^3.
    ///
    /// Stored energy for hyperelastic laws; possibly a reduced incremental
    /// potential for frozen inelastic adapters. Negative potential datums are
    /// retained. This field must not be interpreted as dissipated heat.
    pub potential_density_j_m3: f64,
    /// Original reference volume times material potential density, J.
    ///
    /// Has the same stored-versus-incremental semantics as the density above.
    pub potential_j: f64,
}

impl<L: FiniteStrainPotential> HyperTetProblem<'_, L> {
    /// Recover original element stress, total strain, volume and potential.
    ///
    /// Admits the reference mesh once and returns one constant P1 response per
    /// original tetrahedron. Repeated observations can reuse `Self::prepare`
    /// and the prepared owner's `element_responses` method.
    ///
    /// # Errors
    /// Native input/budget/geometry refusal, malformed displacement, material
    /// domain exit, unrepresentable response or cancellation. No partial field
    /// is returned and no constitutive history is advanced.
    pub fn element_responses(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<HyperTetElementResponse>, HyperTetError> {
        self.prepare(cx)?.element_responses(displacement_m, cx)
    }
}

pub(super) fn recover_all<L: FiniteStrainPotential>(
    problem: &HyperTetProblem<'_, L>,
    elements: &[Element],
    displacement_m: &[f64],
    cx: &Cx<'_>,
) -> Result<Vec<HyperTetElementResponse>, HyperTetError> {
    checkpoint(cx)?;
    problem.check_len(displacement_m)?;
    let mut responses = Vec::with_capacity(elements.len());
    for (element, geometry) in elements.iter().enumerate() {
        if element % 256 == 0 {
            checkpoint(cx)?;
        }
        let deformation = HyperTetProblem::<L>::deformation_gradient(geometry, displacement_m);
        responses.push(recover(problem, element, geometry, deformation)?);
    }
    checkpoint(cx)?;
    Ok(responses)
}

pub(super) fn recover<L: FiniteStrainPotential>(
    problem: &HyperTetProblem<'_, L>,
    element: usize,
    geometry: &Element,
    deformation_gradient: [f64; 9],
) -> Result<HyperTetElementResponse, HyperTetError> {
    let refused = |what: String| HyperTetError::MaterialRefused { element, what };
    let jacobian =
        HyperTetProblem::<L>::deformation_jacobian(&deformation_gradient).map_err(refused)?;
    let (first_piola_pa, potential_density_j_m3) = problem
        .material
        .at(element)
        .stress_and_potential(&deformation_gradient)
        .map_err(|error| refused(error.to_string()))?;
    if first_piola_pa.iter().any(|value| !value.is_finite()) || !potential_density_j_m3.is_finite()
    {
        return Err(refused(
            "native stress or material potential is not finite".into(),
        ));
    }

    let cauchy_pa: [f64; 9] = core::array::from_fn(|index| {
        let (i, j) = (index / 3, index % 3);
        (0..3).fold(0.0, |sum, k| {
            first_piola_pa[3 * i + k].mul_add(deformation_gradient[3 * j + k], sum)
        }) / jacobian
    });

    // E = sym(H) + H^T H/2, H = F-I, retains small strain without subtracting
    // nearly equal metric diagonals. Mirror strain only, never material stress.
    let mut h = deformation_gradient;
    for i in 0..3 {
        h[3 * i + i] -= 1.0;
    }
    let mut green_lagrange_strain = [0.0; 9];
    for i in 0..3 {
        for j in i..3 {
            let quadratic = (0..3).fold(0.0, |sum, k| h[3 * k + i].mul_add(h[3 * k + j], sum));
            let strain = 0.5 * (h[3 * i + j] + h[3 * j + i] + quadratic);
            green_lagrange_strain[3 * i + j] = strain;
            green_lagrange_strain[3 * j + i] = strain;
        }
    }
    let current_volume_m3 = geometry.volume * jacobian;
    let potential_j = geometry.volume * potential_density_j_m3;
    if cauchy_pa
        .iter()
        .chain(&green_lagrange_strain)
        .any(|value| !value.is_finite())
        || !(current_volume_m3.is_finite() && current_volume_m3 > 0.0)
        || !potential_j.is_finite()
    {
        return Err(refused(
            "recovered stress, strain, current volume or integrated potential is not representable"
                .into(),
        ));
    }
    Ok(HyperTetElementResponse {
        element,
        nodes: geometry.nodes,
        deformation_gradient,
        jacobian,
        first_piola_pa,
        cauchy_pa,
        green_lagrange_strain,
        reference_volume_m3: geometry.volume,
        current_volume_m3,
        potential_density_j_m3,
        potential_j,
    })
}
