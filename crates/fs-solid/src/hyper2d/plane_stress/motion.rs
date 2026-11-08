//! Native P1/Q1 motion at fixed temperature, reference and material history.
//!
//! A finite-compliance branch differentiates all three transverse equilibrium
//! equations. F22dot is not zero and is not inferred from incompressibility.
//! An exact-volume branch differentiates the constrained thickness and relaxed
//! transverse shear, including the original pressure multiplier. A plane-strain
//! branch instead holds the transverse deformation at identity.
//! All reuse the original geometry/material-point ordering and report total
//! geometric strain rates, without an elastic/plastic/thermal decomposition.

use super::{
    Element, Geometry, HyperProblem, PlaneStressSettings, Point, checkpoint,
    incompressible_plane_stress, invalid, material_error, plane_stress, reference_seed,
};
use crate::SolidError;
use fs_exec::Cx;

/// Original P1/Q1 material point and its held-material directional motion.
/// All tensors are row-major, with physical tensor shear components.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneMaterialPointMotion {
    /// Original mesh element index.
    pub element: usize,
    /// Native quadrature index: one centroid for P1, four Gauss sites for Q1.
    pub local_point: usize,
    /// Original quadrature area [m²]; multiply by original thickness for volume.
    pub reference_area_m2: f64,
    /// Total three-dimensional F, including relaxed transverse deformation.
    pub deformation_gradient: [f64; 9],
    /// Total Fdot at fixed temperature, natural reference and internal history [1/s].
    pub deformation_gradient_rate_per_s: [f64; 9],
    /// Original first Piola stress, including the actual transverse residual [Pa].
    pub first_piola_pa: [f64; 9],
    /// Full spatial Cauchy stress P F^T/J [Pa]; no stress symmetrization.
    pub cauchy_pa: [f64; 9],
    /// Total material Green--Lagrange strain (F^T F-I)/2.
    pub green_lagrange_strain: [f64; 9],
    /// Material rate sym(F^T Fdot), invariant under superposed in-plane motion [1/s].
    pub green_lagrange_strain_rate_per_s: [f64; 9],
    /// Current spatial velocity gradient L=Fdot F^-1 [1/s].
    pub spatial_velocity_gradient_per_s: [f64; 9],
    /// Objective spatial deformation rate D=sym(L) [1/s].
    pub spatial_strain_rate_per_s: [f64; 9],
    /// Positive current/original volume ratio J.
    pub jacobian: f64,
    /// Jdot=J tr(L) [1/s].
    pub jacobian_rate_per_s: f64,
    /// P:Fdot per original volume [W/m³]. This is stress power, not heat.
    pub reference_stress_power_density_w_m3: f64,
    /// sigma:L per current volume [W/m³]. Multiplication by J recovers
    /// original-volume power. Full stress/L retain any source asymmetry.
    pub spatial_stress_power_density_w_m3: f64,
}

/// Exact-volume P1/Q1 motion, retaining the physical pressure multiplier.
/// The multiplier follows the native remainder's gauge and need not equal
/// minus the mean Cauchy stress. Material sources and history are held fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct IncompressiblePlaneMaterialPointMotion {
    /// Total kinematics and full physical stress/power at the native point.
    pub motion: PlaneMaterialPointMotion,
    /// Compression-positive multiplier in P=Pr-p cof(F) [Pa].
    pub pressure_pa: f64,
    /// Held-material directional rate of that same multiplier [Pa/s].
    pub pressure_rate_pa_s: f64,
}

impl HyperProblem<'_> {
    /// Native plane-strain fields and motion with F22=1 and transverse rates zero.
    ///
    /// Both input buffers contain one finite xy pair per node, in m and m/s.
    /// Temperature, material coefficients, reference mapping and history are held.
    /// No history advances; source-induced stress rates or heat are not inferred.
    /// Reference/current geometry and original material-domain gates remain active.
    pub fn plane_strain_held_material_motion(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<PlaneMaterialPointMotion>, SolidError> {
        self.held_material_motion(displacement_m, velocity_m_s, None, cx)
    }

    /// Native finite-compliance plane-stress motion on the stable relaxed branch.
    ///
    /// Differentiates P[:,2]=0 using the SAME eliminated tangent columns as the
    /// native constitutive response. This includes thickness and transverse shear
    /// rates. It is isothermal and holds the natural reference, coefficients and
    /// internal history fixed. Independent thermal/plastic/relaxation rates are
    /// not supplied by this API; a frozen incremental potential gives only its
    /// held-potential directional rate, not complete evolving-specimen motion.
    ///
    /// Exact-volume plane stress, mean dilatation, shell bending and wrinkling are
    /// different formulations and are not represented by this observation. Loads
    /// are not applied and equilibrium of the supplied displacement is not claimed.
    pub fn plane_stress_held_material_motion(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<PlaneMaterialPointMotion>, SolidError> {
        self.held_material_motion(displacement_m, velocity_m_s, Some(local), cx)
    }

    /// Native exact-incompressibility plane-stress motion at held material state.
    ///
    /// The supplied xy displacements and velocities determine the in-plane
    /// motion. Thickness obeys det(F)=J_reference exactly; both transverse shear
    /// rates follow their consistently differentiated stationary branch. Full
    /// stress includes the original pressure multiplier, whose rate is returned.
    /// No bulk penalty is evaluated, including for finite-bulk native laws.
    ///
    /// Temperature, coefficients, natural reference (including its volume) and
    /// accepted internal history are fixed. A frozen inelastic potential gives
    /// only that directional branch, not an evolving specimen's complete rate.
    /// Loads, history advancement, shell bending and wrinkling are not supplied.
    pub fn incompressible_plane_stress_held_material_motion(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<IncompressiblePlaneMaterialPointMotion>, SolidError> {
        checkpoint(cx)?;
        self.validate_materials_with(|| checkpoint(cx))?;
        let geometry = Geometry::new(self.mesh, cx)?;
        if velocity_m_s.len() != geometry.ndof || velocity_m_s.iter().any(|v| !v.is_finite()) {
            return Err(invalid("one finite xy velocity per original node required"));
        }
        let mut responses = Vec::new();
        geometry.for_points(
            displacement_m,
            0.0,
            cx,
            |index, local_point, element, point, in_plane| {
                let material = self.point_material(index, local_point)?;
                let seed = self.material.transverse_seed_at(index, local_point);
                let response = if let Some(seed) = seed {
                    incompressible_plane_stress::incompressible_plane_stress_with_kinematics(
                        material,
                        &in_plane,
                        Some([seed[0], seed[1]]),
                        local,
                    )
                } else if let Some(native) = material.native_hyperelastic() {
                    native.incompressible_plane_stress_with_kinematics(&in_plane, local)
                } else {
                    incompressible_plane_stress::incompressible_plane_stress_with_kinematics(
                        material, &in_plane, None, local,
                    )
                }
                .map_err(|error| material_error(index, error))?;
                let rate = velocity_gradient(element, point, velocity_m_s);
                let fdot = response
                    .deformation_gradient_rate(rate)
                    .map_err(|error| material_error(index, error))?;
                let pressure_rate_pa_s = response
                    .pressure_rate(rate)
                    .map_err(|error| material_error(index, error))?;
                let state = response.response;
                responses.push(IncompressiblePlaneMaterialPointMotion {
                    motion: recover(
                        index,
                        local_point,
                        point.weight,
                        state.deformation_gradient,
                        fdot,
                        state.full_piola,
                    )?,
                    pressure_pa: state.pressure_pa,
                    pressure_rate_pa_s,
                });
                Ok(())
            },
        )?;
        Ok(responses)
    }

    fn held_material_motion(
        &self,
        u: &[f64],
        velocity: &[f64],
        local: Option<PlaneStressSettings>,
        cx: &Cx<'_>,
    ) -> Result<Vec<PlaneMaterialPointMotion>, SolidError> {
        checkpoint(cx)?;
        self.validate_materials_with(|| checkpoint(cx))?;
        let geometry = Geometry::new(self.mesh, cx)?;
        if velocity.len() != geometry.ndof || velocity.iter().any(|v| !v.is_finite()) {
            return Err(invalid("one finite xy velocity per original node required"));
        }
        let mut responses = Vec::new();
        geometry.for_points(
            u,
            0.0,
            cx,
            |index, local_point, element, point, in_plane| {
                let material = self.point_material(index, local_point)?;
                let rate = velocity_gradient(element, point, velocity);
                let (f, fdot) = if let Some(local) = local {
                    let seed = self.material.transverse_seed_at(index, local_point);
                    let response = if let Some(seed) = seed {
                        plane_stress::plane_stress_with_kinematics(
                            material,
                            &in_plane,
                            Some(seed),
                            local,
                        )
                    } else if let Some(native) = material.native_hyperelastic() {
                        native.plane_stress_with_kinematics(&in_plane, local)
                    } else {
                        plane_stress::plane_stress_with_kinematics(
                            material,
                            &in_plane,
                            Some(reference_seed(material, &in_plane)),
                            local,
                        )
                    }
                    .map_err(|error| material_error(index, error))?;
                    (
                        response.response.deformation_gradient,
                        response
                            .deformation_gradient_rate(rate)
                            .map_err(|error| material_error(index, error))?,
                    )
                } else {
                    (
                        [
                            in_plane[0],
                            in_plane[1],
                            0.0,
                            in_plane[2],
                            in_plane[3],
                            0.0,
                            0.0,
                            0.0,
                            1.0,
                        ],
                        [rate[0], rate[1], 0.0, rate[2], rate[3], 0.0, 0.0, 0.0, 0.0],
                    )
                };
                // Preserve the full original stress and its remaining transverse
                // residual. Setting relaxed entries to zero would hide power error.
                let (stress, potential) = material
                    .stress_and_potential(&f)
                    .map_err(|error| material_error(index, error))?;
                if !potential.is_finite() {
                    return Err(SolidError::MaterialRefused {
                        what: "native motion material potential is nonfinite".into(),
                    });
                }
                responses.push(recover(index, local_point, point.weight, f, fdot, stress)?);
                Ok(())
            },
        )?;
        Ok(responses)
    }
}

fn velocity_gradient(element: &Element, point: &Point, velocity: &[f64]) -> [f64; 4] {
    // Differentiate from zero, before adding any reference axes. Common velocity
    // is removed before contraction, as in the native displacement calculation.
    let mut rate = [0.0_f64; 4];
    for a in 1..element.arity {
        for i in 0..2 {
            let relative = velocity[2 * element.nodes[a] + i] - velocity[2 * element.nodes[0] + i];
            for j in 0..2 {
                rate[2 * i + j] = relative.mul_add(point.gradients[a][j], rate[2 * i + j]);
            }
        }
    }
    rate
}

fn recover(
    element: usize,
    local_point: usize,
    reference_area_m2: f64,
    f: [f64; 9],
    rate: [f64; 9],
    first_piola_pa: [f64; 9],
) -> Result<PlaneMaterialPointMotion, SolidError> {
    let jacobian = (f[0] * f[4] - f[1] * f[3]) * f[8];
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
    let sym =
        |a: &[f64; 9]| core::array::from_fn(|ij| f64::midpoint(a[ij], a[3 * (ij % 3) + ij / 3]));
    let spatial_strain_rate_per_s = sym(&spatial_velocity_gradient_per_s);
    let green_lagrange_strain_rate_per_s = sym(&material_rate);
    let cauchy_pa: [f64; 9] = core::array::from_fn(|ij| {
        (0..3).fold(0.0, |sum, k| {
            first_piola_pa[3 * (ij / 3) + k].mul_add(f[3 * (ij % 3) + k], sum)
        }) / jacobian
    });
    let mut h = f;
    for i in 0..3 {
        h[3 * i + i] -= 1.0;
    }
    let green_lagrange_strain: [f64; 9] = core::array::from_fn(|ij| {
        let (i, j) = (ij / 3, ij % 3);
        let quadratic = (0..3).fold(0.0, |sum, k| h[3 * k + i].mul_add(h[3 * k + j], sum));
        0.5 * (h[ij] + h[3 * j + i] + quadratic)
    });
    let jacobian_rate_per_s = jacobian
        * (spatial_velocity_gradient_per_s[0]
            + spatial_velocity_gradient_per_s[4]
            + spatial_velocity_gradient_per_s[8]);
    let contract =
        |a: &[f64; 9], b: &[f64; 9]| a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(*y, sum));
    let reference_stress_power_density_w_m3 = contract(&first_piola_pa, &rate);
    let spatial_stress_power_density_w_m3 = contract(&cauchy_pa, &spatial_velocity_gradient_per_s);
    if !jacobian.is_finite()
        || jacobian <= 0.0
        || f.iter()
            .chain(&rate)
            .chain(&first_piola_pa)
            .chain(&inverse)
            .chain(&cauchy_pa)
            .chain(&green_lagrange_strain)
            .chain(&spatial_velocity_gradient_per_s)
            .chain(&spatial_strain_rate_per_s)
            .chain(&green_lagrange_strain_rate_per_s)
            .any(|v| !v.is_finite())
        || [
            jacobian_rate_per_s,
            reference_stress_power_density_w_m3,
            spatial_stress_power_density_w_m3,
        ]
        .iter()
        .any(|v| !v.is_finite())
    {
        return Err(SolidError::MaterialRefused {
            what: format!(
                "element {element}, point {local_point}: motion strain, volume rate or stress power is unrepresentable"
            ),
        });
    }
    Ok(PlaneMaterialPointMotion {
        element,
        local_point,
        reference_area_m2,
        deformation_gradient: f,
        deformation_gradient_rate_per_s: rate,
        first_piola_pa,
        cauchy_pa,
        green_lagrange_strain,
        green_lagrange_strain_rate_per_s,
        spatial_velocity_gradient_per_s,
        spatial_strain_rate_per_s,
        jacobian,
        jacobian_rate_per_s,
        reference_stress_power_density_w_m3,
        spatial_stress_power_density_w_m3,
    })
}
