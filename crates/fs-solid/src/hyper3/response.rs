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
//!
//! Motion recovery adds the actual referential velocity gradient, objective
//! strain rates and stress power on that same state. Rates are total geometric
//! rates, not an inferred elastic/plastic/thermal split. Stress power is not
//! heat or a stored-energy rate when a material adapter has evolving history.

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

/// Native finite-deformation fields and their instantaneous motion.
///
/// Every tensor is row-major, with tensor shear components rather than an
/// engineering-shear or Mandel convention. No material history is advanced.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperTetElementMotionResponse {
    /// Original stress, strain, volume and material-potential response.
    pub response: HyperTetElementResponse,
    /// Total referential velocity gradient Fdot = Grad_X(v), 1/s.
    pub deformation_gradient_rate_per_s: [f64; 9],
    /// Current spatial velocity gradient L = Fdot F^-1, 1/s.
    pub spatial_velocity_gradient_per_s: [f64; 9],
    /// Objective spatial deformation rate D = sym(L), 1/s.
    /// Under a time-dependent superposed rotation Q, D becomes Q D Q^T.
    pub spatial_strain_rate_per_s: [f64; 9],
    /// Material total Green--Lagrange rate Edot = sym(F^T Fdot), 1/s.
    /// This is invariant under a superposed spatial rigid motion.
    pub green_lagrange_strain_rate_per_s: [f64; 9],
    /// Jdot = J tr(L), 1/s, with positive values denoting dilation.
    pub jacobian_rate_per_s: f64,
    /// Instantaneous current-volume rate, V_reference Jdot, m^3/s.
    pub current_volume_rate_m3_s: f64,
    /// Original reference stress power V_reference P:Fdot, W.
    pub reference_stress_power_w: f64,
    /// Independently reconstructed current-volume stress power V_current sigma:L, W.
    /// The FULL stress and velocity gradient are retained. For objective laws,
    /// symmetry of sigma also gives sigma:D; asymmetry is never hidden here.
    pub spatial_stress_power_w: f64,
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

    /// Recover native finite-strain fields, objective rates and stress power.
    ///
    /// Displacement and velocity use the same original node-major xyz order,
    /// in metres and metres/second respectively. Repeated observations should
    /// reuse the prepared owner; no equilibrium solve or history update occurs.
    ///
    /// # Errors
    /// Native mesh/budget/material refusal, malformed motion, unrepresentable
    /// rate or power, or cancellation. No partial field is published.
    pub fn element_motion_responses(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<HyperTetElementMotionResponse>, HyperTetError> {
        self.prepare(cx)?
            .element_motion_responses(displacement_m, velocity_m_s, cx)
    }
}

pub(super) fn recover_motion_all<L: FiniteStrainPotential>(
    problem: &HyperTetProblem<'_, L>,
    elements: &[Element],
    displacement_m: &[f64],
    velocity_m_s: &[f64],
    cx: &Cx<'_>,
) -> Result<Vec<HyperTetElementMotionResponse>, HyperTetError> {
    checkpoint(cx)?;
    problem.check_len(displacement_m)?;
    if velocity_m_s.len() != displacement_m.len()
        || velocity_m_s.iter().any(|value| !value.is_finite())
    {
        return Err(HyperTetError::InvalidInput {
            what: "one finite xyz velocity per original node required".into(),
        });
    }
    let mut responses = Vec::with_capacity(elements.len());
    for (element, geometry) in elements.iter().enumerate() {
        if element % 256 == 0 {
            checkpoint(cx)?;
        }
        let (f, rate) = motion_gradients_with::<L>(geometry, |node| {
            (
                core::array::from_fn(|i| displacement_m[3 * node + i]),
                core::array::from_fn(|i| velocity_m_s[3 * node + i]),
            )
        })?;
        responses.push(with_motion(recover(problem, element, geometry, f)?, rate)?);
    }
    checkpoint(cx)?;
    Ok(responses)
}

pub(super) fn motion_gradients_with<L: FiniteStrainPotential>(
    geometry: &Element,
    mut motion_at_node: impl FnMut(usize) -> ([f64; 3], [f64; 3]),
) -> Result<([f64; 9], [f64; 9]), HyperTetError> {
    let mut rate = [0.0; 9];
    let mut first_velocity = [0.0; 3];
    let mut local_node = 0;
    let mut finite = true;
    // Reuse the native F owner and its original callback order. Differentiate
    // velocity directly from zero: forming I+Grad(v) would erase tiny rates.
    let f = HyperTetProblem::<L>::deformation_gradient_with(geometry, |node| {
        let (displacement, velocity) = motion_at_node(node);
        finite &= displacement.iter().chain(&velocity).all(|v| v.is_finite());
        if local_node == 0 {
            first_velocity = velocity;
        } else {
            for i in 0..3 {
                let relative = velocity[i] - first_velocity[i];
                for j in 0..3 {
                    rate[3 * i + j] =
                        relative.mul_add(geometry.gradients[local_node][j], rate[3 * i + j]);
                }
            }
        }
        local_node += 1;
        displacement
    });
    if !finite || rate.iter().any(|v| !v.is_finite()) {
        return Err(HyperTetError::InvalidInput {
            what: "finite nodal motion and representable referential velocity gradient required"
                .into(),
        });
    }
    Ok((f, rate))
}

pub(super) fn with_motion(
    response: HyperTetElementResponse,
    rate: [f64; 9],
) -> Result<HyperTetElementMotionResponse, HyperTetError> {
    let fields = motion_fields(
        response.deformation_gradient,
        rate,
        response.jacobian,
        response.reference_volume_m3,
        response.current_volume_m3,
        &response.first_piola_pa,
        &response.cauchy_pa,
    )
    .map_err(|what| HyperTetError::MaterialRefused {
        element: response.element,
        what: what.into(),
    })?;
    Ok(HyperTetElementMotionResponse {
        response,
        deformation_gradient_rate_per_s: rate,
        spatial_velocity_gradient_per_s: fields.spatial_velocity_gradient_per_s,
        spatial_strain_rate_per_s: fields.spatial_strain_rate_per_s,
        green_lagrange_strain_rate_per_s: fields.green_lagrange_strain_rate_per_s,
        jacobian_rate_per_s: fields.jacobian_rate_per_s,
        current_volume_rate_m3_s: fields.current_volume_rate_m3_s,
        reference_stress_power_w: fields.reference_stress_power_w,
        spatial_stress_power_w: fields.spatial_stress_power_w,
    })
}

// One tensor calculation shared by native P1 and mixed-quadrature observers.
// Geometry, material admission and pressure/history ownership stay with callers.
pub(super) struct MotionFields {
    pub(super) spatial_velocity_gradient_per_s: [f64; 9],
    pub(super) spatial_strain_rate_per_s: [f64; 9],
    pub(super) green_lagrange_strain_rate_per_s: [f64; 9],
    pub(super) jacobian_rate_per_s: f64,
    pub(super) current_volume_rate_m3_s: f64,
    pub(super) reference_stress_power_w: f64,
    pub(super) spatial_stress_power_w: f64,
}

pub(super) fn motion_fields(
    f: [f64; 9],
    rate: [f64; 9],
    jacobian: f64,
    reference_volume_m3: f64,
    current_volume_m3: f64,
    first_piola_pa: &[f64; 9],
    cauchy_pa: &[f64; 9],
) -> Result<MotionFields, &'static str> {
    let inverse: [f64; 9] = core::array::from_fn(|ij| {
        let (i, j) = (ij / 3, ij % 3);
        let (a, b, c, d) = ((j + 1) % 3, (j + 2) % 3, (i + 1) % 3, (i + 2) % 3);
        (f[3 * a + c] * f[3 * b + d] - f[3 * a + d] * f[3 * b + c]) / jacobian
    });
    let spatial_velocity_gradient_per_s: [f64; 9] = core::array::from_fn(|ij| {
        (0..3).fold(0.0, |sum, k| {
            rate[3 * (ij / 3) + k].mul_add(inverse[3 * k + ij % 3], sum)
        })
    });
    let material_rate: [f64; 9] = core::array::from_fn(|ij| {
        (0..3).fold(0.0, |sum, k| {
            f[3 * k + ij / 3].mul_add(rate[3 * k + ij % 3], sum)
        })
    });
    let symmetric =
        |a: &[f64; 9]| core::array::from_fn(|ij| f64::midpoint(a[ij], a[3 * (ij % 3) + ij / 3]));
    let spatial_strain_rate_per_s = symmetric(&spatial_velocity_gradient_per_s);
    let green_lagrange_strain_rate_per_s = symmetric(&material_rate);
    let jacobian_rate_per_s = jacobian
        * (spatial_velocity_gradient_per_s[0]
            + spatial_velocity_gradient_per_s[4]
            + spatial_velocity_gradient_per_s[8]);
    let current_volume_rate_m3_s = reference_volume_m3 * jacobian_rate_per_s;
    let contraction =
        |a: &[f64; 9], b: &[f64; 9]| a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(*y, sum));
    let reference_stress_power_w = reference_volume_m3 * contraction(first_piola_pa, &rate);
    let spatial_stress_power_w =
        current_volume_m3 * contraction(cauchy_pa, &spatial_velocity_gradient_per_s);
    if inverse
        .iter()
        .chain(&spatial_velocity_gradient_per_s)
        .chain(&spatial_strain_rate_per_s)
        .chain(&green_lagrange_strain_rate_per_s)
        .any(|value| !value.is_finite())
        || [
            jacobian_rate_per_s,
            current_volume_rate_m3_s,
            reference_stress_power_w,
            spatial_stress_power_w,
        ]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err("native motion strain rate, volume rate or stress power is not representable");
    }
    Ok(MotionFields {
        spatial_velocity_gradient_per_s,
        spatial_strain_rate_per_s,
        green_lagrange_strain_rate_per_s,
        jacobian_rate_per_s,
        current_volume_rate_m3_s,
        reference_stress_power_w,
        spatial_stress_power_w,
    })
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

    let (cauchy_pa, green_lagrange_strain) =
        stress_strain_fields(&deformation_gradient, &first_piola_pa, jacobian);
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

// Preserve the full original stress while recovering total geometric strain.
// This uses the same arithmetic for P1 and mixed quadrature, without smoothing.
pub(super) fn stress_strain_fields(
    deformation_gradient: &[f64; 9],
    first_piola_pa: &[f64; 9],
    jacobian: f64,
) -> ([f64; 9], [f64; 9]) {
    let cauchy_pa: [f64; 9] = core::array::from_fn(|index| {
        let (i, j) = (index / 3, index % 3);
        (0..3).fold(0.0, |sum, k| {
            first_piola_pa[3 * i + k].mul_add(deformation_gradient[3 * j + k], sum)
        }) / jacobian
    });

    (cauchy_pa, green_lagrange_strain(deformation_gradient))
}

pub(super) fn green_lagrange_strain(deformation_gradient: &[f64; 9]) -> [f64; 9] {
    // E = sym(H) + H^T H/2, H = F-I, retains small strain without subtracting
    // nearly equal metric diagonals. Mirror strain only, never material stress.
    let mut h = *deformation_gradient;
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
    green_lagrange_strain
}
