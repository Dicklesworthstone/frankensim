//! Complete original tensor memory for contiguous equivalent-single-layer
//! stacks. Each layer keeps its own normal-traction constraint, orientation,
//! shear correction, material band and history. No branchwise plane-stress
//! condensation, scalar coating fit or drilling loss is introduced.
use super::*;
use fs_math::dd::Dd;

impl TensorShellRelaxation {
    /// Original source layers on every physical MITC facet. Each stack is
    /// nonempty, ordered from its lower to upper face and exactly contiguous.
    /// Faces are signed coordinates about the ORIGINAL shell reference surface.
    ///
    /// Every layer contributes its own seven area and two thickness points.
    /// The sum must reproduce the original equilibrium A/B/D, transverse shear
    /// and all three density moments; no reference mass or stiffness is changed.
    /// The supported equilibrium section has no membrane/shear or bending/shear
    /// cross block. Layer cross terms may cancel only through their actual
    /// integrated section, which is checked before projection.
    ///
    /// This retains the native equivalent-single-layer kinematics with an
    /// authored work-conjugate shear correction in each layer. It does not
    /// resolve interlaminar shear/normal stress, slip, delamination or warped
    /// through-thickness displacement. Those need a resolved spatial owner.
    pub fn from_layered_mitc_shell(
        model: &MitcShellModel,
        facets: &[Vec<TensorShellMaterial>],
        modes: &[Vec<f64>],
        limits: TensorShellRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if facets.len() > limits.max_material_points.min(4096) / 14 {
            return Err(invalid(
                "layered tensor shell material points exceed their bound",
            ));
        }
        let stacks: Vec<_> = facets.iter().map(Vec::as_slice).collect();
        Self::from_material_stacks(model, &stacks, modes, limits, gate)
    }

    pub(super) fn from_material_stacks(
        model: &MitcShellModel,
        stacks: &[&[TensorShellMaterial]],
        modes: &[Vec<f64>],
        limits: TensorShellRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let (mesh, directors, sections) = model.original_geometry();
        let n = modes.len();
        let materials = stacks
            .iter()
            .try_fold(0usize, |sum, stack| sum.checked_add(stack.len()))
            .ok_or_else(|| invalid("tensor shell material extent overflows"))?;
        let count = materials
            .checked_mul(14)
            .ok_or_else(|| invalid("tensor shell material point extent overflows"))?;
        if !(1..=64).contains(&limits.max_modes)
            || n == 0
            || n > limits.max_modes
            || !(1..=4096).contains(&limits.max_material_points)
            || limits.max_memory_states > 8192
            || limits.max_projection_terms == 0
            || stacks.len() != mesh.tris.len()
            || stacks.is_empty()
            || stacks.iter().any(|stack| stack.is_empty())
            || count > limits.max_material_points
            || modes.iter().any(|v| v.len() != model.pencil().free)
        {
            return Err(invalid(
                "tensor shell needs complete original facets, modes and bounded point/history/work capacities",
            ));
        }
        for stack in stacks {
            poll(gate)?;
            if stack
                .windows(2)
                .any(|layers| layers[0].faces_m[1] != layers[1].faces_m[0])
            {
                return Err(invalid(
                    "tensor shell layers must be ordered and exactly contiguous without gaps or overlaps",
                ));
            }
        }
        let states = stacks
            .iter()
            .flat_map(|stack| stack.iter())
            .try_fold(0usize, |sum, layer| {
                let normal = usize::from(
                    layer
                        .law
                        .kv_pa_s()
                        .is_some_and(|d| d[2].iter().any(|v| *v != 0.0)),
                );
                layer
                    .law
                    .branches()
                    .len()
                    .checked_mul(6)
                    .and_then(|r| r.checked_add(normal))
                    .and_then(|r| r.checked_mul(14))
                    .and_then(|r| sum.checked_add(r))
            })
            .ok_or_else(|| invalid("tensor shell material history extent overflows"))?;
        // Original geometry and condensed fields are prepared once per FACET;
        // each layer owns only its constitutive image and thickness quadrature.
        let setup = stacks.len() as u128
            * (PreparedElementStrain::PREPARATION_WORK + n * PreparedElementStrain::REAL_FIELD_WORK)
                as u128
            + model.pencil().m.nnz() as u128 * n as u128
            + 16 * model.pencil().k.nnz() as u128 * n as u128
            + 17 * model.pencil().free as u128 * (n as u128).pow(2)
            + mesh.nodes.len() as u128 * 6 * n as u128
            + count as u128 * n as u128 * 256
            + count as u128 * (n as u128).pow(2)
            + 64 * (n as u128).pow(3)
            + if materials == stacks.len() {
                0
            } else {
                materials as u128 * 256
            };
        if setup >= limits.max_projection_terms as u128 {
            return Err(invalid(
                "tensor shell original geometry and modal pencil exceed projection work",
            ));
        }
        let work = limits.max_projection_terms - setup as usize;
        let factor_work = work / 3;
        if states.checked_mul(n).is_none_or(|v| v > factor_work) {
            return Err(invalid(
                "tensor shell original histories exceed their explicit projection work",
            ));
        }
        let point_work = (work - factor_work) / 2 / (count + materials);
        if point_work == 0 {
            return Err(invalid(
                "tensor shell constitutive projection requires positive work",
            ));
        }
        validate_mass_basis(model.pencil(), modes, gate)?;
        let equilibrium_omega_squared = project_original_stiffness(model, modes, gate)?;
        let physical = modes
            .iter()
            .map(|v| {
                poll(gate)?;
                model.expand(v).map_err(plate_error)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut rows = Vec::with_capacity(states);
        let mut times = Vec::with_capacity(states);
        let mut points = Vec::with_capacity(count);
        let mut viscosity = vec![Dd::ZERO; n * n];
        let mut band = (0.0_f64, f64::MAX);
        let mut catalog = Vec::with_capacity(materials);
        let mut facet_materials = Vec::with_capacity(stacks.len());
        let rule = integration_rule();
        for (e, stack) in stacks.iter().enumerate() {
            poll(gate)?;
            let first = catalog.len();
            let mut constitutive = Vec::with_capacity(stack.len());
            for layer in *stack {
                poll(gate)?;
                let unit: [[f64; 5]; 5] = core::array::from_fn(|i| {
                    core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 })
                });
                constitutive.push(
                    PlaneStressPoint::new_shell(
                        &layer.law,
                        &unit,
                        1.0,
                        TensorProjectionBudget {
                            max_coordinates: 6,
                            max_scalar_terms: point_work,
                        },
                        gate,
                    )
                    .map_err(tensor_error)?
                    .viscoelastic_image(point_work, gate)
                    .map_err(tensor_error)?,
                );
            }
            if stack.len() == 1 {
                // Preserve the established homogeneous admission and arithmetic.
                validate_section(
                    model,
                    e,
                    &sections[e],
                    &stack[0],
                    constitutive[0].equilibrium_stiffness(),
                )?;
            } else {
                validate_layered_section(model, e, &sections[e], stack, &constitutive, gate)?;
            }
            let observer =
                PreparedElementStrain::new(mesh, e, directors, &sections[e], &mut || {
                    gate.is_requested()
                })
                .map_err(plate_error)?;
            let tri = mesh.tris[e];
            let fields = physical
                .iter()
                .map(|u| {
                    observer
                        .real_field(
                            &core::array::from_fn(|i| u[6 * tri[i / 6] + i % 6]),
                            &mut || gate.is_requested(),
                        )
                        .map_err(plate_error)
                })
                .collect::<Result<Vec<_>, _>>()?;
            let area = mesh.facet(e).map_err(plate_error)?.area_m2;
            // Reuse all original in-plane samples for every layer.
            let strains = rule
                .iter()
                .map(|&(lambda, _)| {
                    fields
                        .iter()
                        .map(|field| field.at_barycentric(lambda).map_err(plate_error))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()?;
            for (layer, constitutive) in stack.iter().zip(&constitutive) {
                poll(gate)?;
                let material = catalog.len();
                catalog.push(layer.clone());
                let source_band = layer.law.band_rad_s();
                band = (band.0.max(source_band.0), band.1.min(source_band.1));
                let h = layer.faces_m[1] - layer.faces_m[0];
                let center = layer.faces_m[0] + 0.5 * h;
                let offset = 0.5 * h / det::sqrt(3.0);
                let shear = det::sqrt(layer.shear_correction);
                for (&(lambda, weight), strains) in rule.iter().zip(&strains) {
                    for z in [center - offset, center + offset] {
                        poll(gate)?;
                        let columns: Vec<[f64; 5]> = strains
                            .iter()
                            .map(|s| {
                                [
                                    s.membrane[0] + z * s.gradient_per_m[0],
                                    s.membrane[1] + z * s.gradient_per_m[1],
                                    s.membrane[2] + z * s.gradient_per_m[2],
                                    shear * s.transverse_shear[1],
                                    shear * s.transverse_shear[0],
                                ]
                            })
                            .collect();
                        let image = constitutive
                            .project_shell(
                                &columns,
                                area * weight * h * 0.5,
                                TensorProjectionBudget {
                                    max_coordinates: n + 1,
                                    max_scalar_terms: point_work,
                                },
                                gate,
                            )
                            .map_err(tensor_error)?;
                        let start = rows.len();
                        for (total, &value) in
                            viscosity.iter_mut().zip(image.feedthrough_viscosity())
                        {
                            *total = *total + Dd::from_f64(value);
                        }
                        rows.extend_from_slice(image.coupling_rows());
                        times.extend(
                            image
                                .relaxation_rates_per_s()
                                .iter()
                                .map(|rate| rate.recip()),
                        );
                        points.push(TensorShellMaterialPoint {
                            facet: e,
                            material,
                            barycentric: lambda,
                            thickness_coordinate_m: z,
                            memory: start..rows.len(),
                            memory_map: None,
                            image,
                        });
                    }
                }
            }
            facet_materials.push(first..catalog.len());
        }
        if band.0 != 0.0 || band.1 < band.0 || rows.len() != states || points.len() != count {
            return Err(invalid(
                "tensor shell needs complete zero-to-transient source support and original histories",
            ));
        }
        if equilibrium_omega_squared
            .iter()
            .any(|&value| det::sqrt(value) > band.1)
        {
            return Err(invalid(
                "original shell modal frequencies exceed source support",
            ));
        }
        let shared = RelaxationCoordinates::new(
            &rows,
            &times,
            n,
            RelaxationCoordinateBudget {
                maximum_states: limits.max_memory_states,
                maximum_scalar_terms: factor_work,
            },
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_material::visco_tensor::TensorPronyError::Cancelled => {
                AcousticRealizeError::Cancelled
            }
            other => owner(other),
        })?;
        let viscosity = ModalPlateViscosity::from_psd_modal_matrix(
            viscosity.into_iter().map(|value| value.to_f64()).collect(),
            n,
            1.0,
            gate,
        )?;
        for point in &mut points {
            poll(gate)?;
            point.memory_map = Some(shared.memory_map(point.memory.clone()).map_err(owner)?);
        }
        rows = shared.coupling_rows().to_vec();
        times = shared.relaxation_times_s().to_vec();
        poll(gate)?;
        Ok(Self {
            facets: catalog,
            facet_materials,
            points,
            band_rad_s: band,
            shared,
            viscosity,
            equilibrium_omega_squared,
            modal: ModalShellRelaxation {
                projected: Vec::new(),
                rows,
                times,
                modes: n,
                limits: PlateRelaxationLimits {
                    max_modes: limits.max_modes,
                    max_memory_states: limits.max_memory_states,
                    max_projection_terms: limits.max_projection_terms,
                },
            },
        })
    }
}

// Generalized section coordinates are [membrane(3), curvature(3), gamma_xz,
// gamma_yz]. The constrained material order is [xx, yy, xy, yz, xz].
// Integrating one congruence retains signed B and both shear cross blocks.
fn add_layer_resultants(
    section: &mut [Dd; 64],
    density: &mut [Dd; 3],
    c: &[f64],
    moments: [f64; 3],
    shear_correction: f64,
    density_kg_m3: f64,
) {
    let shear = det::sqrt(shear_correction);
    let axes = [0, 1, 2, 0, 1, 2, 4, 3];
    let degree = [0, 0, 0, 1, 1, 1, 0, 0];
    let weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, shear, shear];
    for i in 0..8 {
        for j in 0..8 {
            let value = Dd::from_f64(c[5 * axes[i] + axes[j]])
                * Dd::from_f64(moments[degree[i] + degree[j]])
                * Dd::from_f64(weights[i])
                * Dd::from_f64(weights[j]);
            section[8 * i + j] = section[8 * i + j] + value;
        }
    }
    for i in 0..3 {
        density[i] = density[i] + Dd::from_f64(moments[i]) * Dd::from_f64(density_kg_m3);
    }
}

fn validate_layered_section(
    model: &MitcShellModel,
    facet: usize,
    section: &ShearShellSection,
    layers: &[TensorShellMaterial],
    images: &[PlaneStressRelaxation],
    gate: &CancelGate,
) -> Result<(), AcousticRealizeError> {
    let mut expected = [Dd::ZERO; 64];
    let mut density = [Dd::ZERO; 3];
    for (layer, image) in layers.iter().zip(images) {
        poll(gate)?;
        add_layer_resultants(
            &mut expected,
            &mut density,
            image.equilibrium_stiffness(),
            layer.moments()?,
            layer.shear_correction,
            layer.source.density_kg_m3(),
        );
    }
    let expected = expected.map(|v| v.to_f64());
    let density = density.map(|v| v.to_f64());
    let close = |actual: f64, expected: f64, scale: f64| {
        actual.is_finite()
            && expected.is_finite()
            && scale.is_finite()
            && scale > 0.0
            && (actual - expected).abs() <= 1e-9 * scale
    };
    let (a, b, d, shear) = section.coupled_matrices();
    for i in 0..8 {
        for j in 0..8 {
            let actual = match (i, j) {
                (0..=2, 0..=2) => a[3 * i + j],
                (0..=2, 3..=5) => b[3 * i + j - 3],
                (3..=5, 0..=2) => b[3 * j + i - 3],
                (3..=5, 3..=5) => d[3 * (i - 3) + j - 3],
                (6..=7, 6..=7) => shear[i - 6][j - 6],
                _ => 0.0,
            };
            let scale = det::sqrt(expected[8 * i + i]) * det::sqrt(expected[8 * j + j]);
            if !close(actual, expected[8 * i + j], scale) {
                return Err(invalid(
                    "tensor shell layers must reproduce complete original equilibrium and shear section work",
                ));
            }
        }
    }
    let scales = [
        density[0],
        det::sqrt(density[0]) * det::sqrt(density[2]),
        density[2],
    ];
    for inertia in model
        .reference_inertia_moments(facet)
        .ok_or_else(|| invalid("original layered shell inertia is absent"))?
    {
        let actual = [
            inertia.areal_mass_kg_m2(),
            inertia.mass_first_moment_kg_m(),
            inertia.rotary_inertia_kg(),
        ];
        if (0..3).any(|i| !close(actual[i], density[i], scales[i])) {
            return Err(invalid(
                "tensor shell layers must reproduce every original reference inertia moment",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::dd::Dd;
    fn moments(lo: f64, hi: f64) -> [f64; 3] {
        let h = hi - lo;
        let center = lo + 0.5 * h;
        [h, h * center, h * (center * center + h * h / 12.0)]
    }
    fn law() -> [f64; 25] {
        let a: [[f64; 5]; 5] = [
            [2.0, 0.2, -0.3, 0.8, -0.4],
            [0.0, 3.0, 0.5, -0.2, 0.7],
            [0.0, 0.0, 1.5, 0.1, 0.4],
            [0.0, 0.0, 0.0, 2.7, -0.3],
            [0.0, 0.0, 0.0, 0.0, 1.1],
        ];
        core::array::from_fn(|k| (0..5).map(|r| a[r][k / 5] * a[r][k % 5]).sum())
    }
    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 3e-14 * a.abs().max(b.abs()).max(1e-100),
            "{a:e} vs {b:e}"
        );
    }
    #[test]
    fn g1_layer_resultants_match_independent_original_strain_quadrature() {
        let c = law();
        let layers = [
            (-0.017, -0.003, 0.6, 1000.0),
            (-0.003, 0.004, 0.9, 1700.0),
            (0.004, 0.025, 0.7, 800.0),
        ];
        let mut section = [Dd::ZERO; 64];
        let mut density = [Dd::ZERO; 3];
        for &(lo, hi, kappa, rho) in &layers {
            add_layer_resultants(&mut section, &mut density, &c, moments(lo, hi), kappa, rho);
        }
        let csec = section.map(|x| x.to_f64());
        let d = density.map(|x| x.to_f64());
        for k in 0..129 {
            let q: [f64; 8] = core::array::from_fn(|i| ((7 * i + 11 * k) % 19) as f64 / 7.0 - 1.0);
            let actual: (f64, f64) = (
                (0..8)
                    .flat_map(|i| (0..8).map(move |j| q[i] * csec[8 * i + j] * q[j]))
                    .sum(),
                d[0] * q[0] * q[0] + 2.0 * d[1] * q[0] * q[3] + d[2] * q[3] * q[3],
            );
            let mut work = 0.0;
            let mut inertia = 0.0;
            // Three-point Gauss integration uses the ORIGINAL five-dimensional
            // source, not the section block or its thickness moment formula.
            for &(lo, hi, kappa, rho) in &layers {
                for (xi, weight) in [
                    (-(3.0_f64 / 5.0).sqrt(), 5.0 / 9.0),
                    (0.0, 8.0 / 9.0),
                    ((3.0_f64 / 5.0).sqrt(), 5.0 / 9.0),
                ] {
                    let z = 0.5 * (lo + hi) + 0.5 * (hi - lo) * xi;
                    let e = [
                        q[0] + z * q[3],
                        q[1] + z * q[4],
                        q[2] + z * q[5],
                        kappa.sqrt() * q[7],
                        kappa.sqrt() * q[6],
                    ];
                    let sigma: [f64; 5] =
                        core::array::from_fn(|i| (0..5).map(|j| c[5 * i + j] * e[j]).sum());
                    work += 0.5
                        * (hi - lo)
                        * weight
                        * e.iter().zip(sigma).map(|(a, b)| a * b).sum::<f64>();
                    inertia += 0.5 * (hi - lo) * weight * rho * (q[0] + z * q[3]).powi(2);
                }
            }
            near(actual.0, work);
            near(actual.1, inertia);
            assert!(actual.0 >= 0.0 && actual.1 >= 0.0);
        }
    }
    #[test]
    fn g3_splitting_identical_layers_preserves_signed_work_and_inertia() {
        let c = law();
        let mut one = [Dd::ZERO; 64];
        let mut od = [Dd::ZERO; 3];
        add_layer_resultants(&mut one, &mut od, &c, moments(-0.013, 0.021), 0.7, 1400.0);
        let mut split = [Dd::ZERO; 64];
        let mut sd = [Dd::ZERO; 3];
        for (lo, hi) in [(-0.013, -0.007), (-0.007, 0.002), (0.002, 0.021)] {
            add_layer_resultants(&mut split, &mut sd, &c, moments(lo, hi), 0.7, 1400.0);
        }
        for (a, b) in one.iter().zip(split) {
            near(a.to_f64(), b.to_f64());
        }
        for (a, b) in od.iter().zip(sd) {
            near(a.to_f64(), b.to_f64());
        }
        assert!(one[3].to_f64() > 0.0 && od[1].to_f64() > 0.0);
    }
}
