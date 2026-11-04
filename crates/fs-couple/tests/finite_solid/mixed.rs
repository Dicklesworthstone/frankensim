//! Original-volume mixed thermal potentials and native compatible prestress.
use super::*;
use fs_couple::finite_solid::thermal::MeanThermalReference;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpProblem, MixedUpSettings};
use fs_solid::hyper3::mixed::vibration::MixedVibrationSettings;
use fs_solid::hyper3::vibration::VibrationSettings;
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

const BULK: f64 = 40.0;
const SHEAR: f64 = 15.0;
const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const Q: [[f64; 3]; 3] = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];

fn det(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3]).sum()
    })
}
fn thermal_tensor(stretches: [f64; 3]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).map(|k| Q[ij / 3][k] * stretches[k] * Q[ij % 3][k]).sum()
    })
}
fn potential(
    f: &[f64; 9],
    inverse_thermal: &[f64; 9],
    jth: f64,
    p: f64,
    finite: bool,
) -> f64 {
    let fe = multiply(f, inverse_thermal);
    let remainder = 0.5 * SHEAR
        * (fe.iter().map(|v| v * v).sum::<f64>() * det(&fe).powf(-2.0 / 3.0) - 3.0);
    jth * remainder + p * (det(f) - jth)
        - if finite { jth * p * p / (2.0 * BULK) } else { 0.0 }
}
fn mean_at(card: &MaterialCard, temperature: f64) -> FiniteSolidMaterial {
    let quantity = QuantitySpec::dimensional(Temperature::DIMS);
    let query = QueryPoint::new()
        .with_quantity("T_reference", quantity, 300.0).unwrap()
        .with_quantity("T_endpoint", quantity, temperature).unwrap();
    let thermal = MeanThermalReference::resolve(
        card,
        &PropertyKey::new(
            "mean_engineering_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
        ),
        &query,
        "T_reference",
        "T_endpoint",
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    FiniteSolidMaterial::resolve_mean_thermal(
        card, 0, thermal, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
fn mixed_problem<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: MixedLawField<'a, FiniteSolidMaterial>,
    fixed: &'a [(usize, f64)],
    treatment: VolumetricTreatment,
) -> MixedUpProblem<'a, FiniteSolidMaterial> {
    MixedUpProblem {
        nodes_m: nodes,
        elements,
        material,
        volumetric: treatment,
        prescribed_m: fixed,
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        pressure_gauge: None,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 10_000,
        settings: MixedUpSettings {
            load_steps: 1,
            maximum_newton_iterations: 24,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            length_scale_m: 1.0,
            stress_scale_pa: 40.0,
            maximum_assembly_bytes: 16 * 1024 * 1024,
            maximum_factor_bytes: 16 * 1024 * 1024,
        },
    }
}
fn modal_settings() -> MixedVibrationSettings {
    MixedVibrationSettings {
        spectrum: VibrationSettings {
            maximum_free_dofs: 32,
            relative_residual_tolerance: 1e-8,
        },
        maximum_mass_assembly_bytes: 8 * 1024 * 1024,
        maximum_condensation_bytes: 8 * 1024 * 1024,
        pressure_residual_tolerance: 1e-10,
    }
}

#[test]
fn oriented_thermal_mixed_potential_has_all_work_conjugate_derivatives() {
    let stretches = [1.10, 0.96, 1.04];
    let card = super::expansion::card(stretches, "T", false);
    let material = super::expansion::bind_expansion(&card, Q, 310.0);
    material.validate_mixed().unwrap();
    let fth = thermal_tensor(stretches);
    let inverse = thermal_tensor(stretches.map(|s| 1.0 / s));
    let jth = stretches.into_iter().product::<f64>();
    let fe = [1.07, 0.13, -0.04, 0.05, 0.94, 0.08, 0.0, -0.03, 1.02];
    let f = multiply(&fe, &fth);
    let pressure = 3.4;
    for treatment in [
        VolumetricTreatment::OriginalMaterial,
        VolumetricTreatment::ExactIncompressible,
    ] {
        let finite = treatment == VolumetricTreatment::OriginalMaterial;
        let response = material.mixed_response(&f, pressure, treatment).unwrap();
        let scalar = material.mixed_pressure_potential(&f, pressure, treatment).unwrap();
        near(scalar.compliance_per_pa, -response.pressure_tangent_per_pa, 2e-15);
        near(scalar.volume_ratio, response.state.bulk_volume_ratio, 2e-14);
        near(scalar.volume_increment, scalar.volume_ratio - 1.0, 2e-14);
        near(scalar.bulk_stored_pa, response.state.bulk_stored_pa, 2e-14);
        near(scalar.complementary_pa, response.state.bulk_complementary_pa, 2e-14);
        near(
            response.state.remainder_potential_pa + pressure * (det(&f) - 1.0)
                - scalar.complementary_pa,
            response.state.potential_pa,
            3e-12,
        );
        near(response.state.potential_pa, potential(&f, &inverse, jth, pressure, finite), 3e-12);
        near(response.state.volume_ratio, det(&f), 2e-14);
        near(
            response.state.bulk_volume_ratio,
            jth * if finite { 1.0 + pressure / BULK } else { 1.0 },
            2e-14,
        );
        near(response.state.pressure_residual, det(&f) - response.state.bulk_volume_ratio, 2e-14);
        near(response.pressure_tangent_per_pa, if finite { -jth / BULK } else { 0.0 }, 2e-15);
        let epsilon = 1e-6;
        for column in 0..9 {
            let mut plus = f;
            plus[column] += epsilon;
            let mut minus = f;
            minus[column] -= epsilon;
            near(
                (potential(&plus, &inverse, jth, pressure, finite)
                    - potential(&minus, &inverse, jth, pressure, finite)) / (2.0 * epsilon),
                response.state.first_piola_pa[column],
                5e-8,
            );
            let a = material.mixed_state(&plus, pressure, treatment).unwrap();
            let b = material.mixed_state(&minus, pressure, treatment).unwrap();
            near(
                (a.pressure_residual - b.pressure_residual) / (2.0 * epsilon),
                response.pressure_coupling[column],
                5e-9,
            );
            near(
                (det(&plus) - det(&minus)) / (2.0 * epsilon),
                response.pressure_coupling[column],
                5e-9,
            );
            for row in 0..9 {
                near(
                    (a.first_piola_pa[row] - b.first_piola_pa[row]) / (2.0 * epsilon),
                    response.tangent_pa[row][column],
                    3e-6,
                );
                near(response.tangent_pa[row][column], response.tangent_pa[column][row], 2e-12);
            }
        }
        let a = material.mixed_state(&f, pressure + epsilon, treatment).unwrap();
        let b = material.mixed_state(&f, pressure - epsilon, treatment).unwrap();
        let ap = material.mixed_pressure_potential(&f, pressure + epsilon, treatment).unwrap();
        let bp = material.mixed_pressure_potential(&f, pressure - epsilon, treatment).unwrap();
        near(
            (ap.complementary_pa - bp.complementary_pa) / (2.0 * epsilon),
            scalar.volume_increment, 2e-9,
        );
        near(
            (ap.volume_increment - bp.volume_increment) / (2.0 * epsilon),
            scalar.compliance_per_pa, 2e-9,
        );
        near(
            (a.potential_pa - b.potential_pa) / (2.0 * epsilon),
            response.state.pressure_residual,
            2e-8,
        );
        near(
            (a.pressure_residual - b.pressure_residual) / (2.0 * epsilon),
            response.pressure_tangent_per_pa,
            2e-9,
        );
        for i in 0..9 {
            near(
                (a.first_piola_pa[i] - b.first_piola_pa[i]) / (2.0 * epsilon),
                response.pressure_coupling[i],
                4e-9,
            );
        }
    }
    // Eliminate the original physical pressure at p=U'(Je). The complete
    // Schur tangent must recover the existing displacement material exactly.
    let p = BULK * (det(&fe) - 1.0);
    let response = material.mixed_response(&f, p, VolumetricTreatment::OriginalMaterial).unwrap();
    let tangent = material.consistent_tangent(&f).unwrap();
    let stress = material.first_piola(&f).unwrap();
    near(response.state.potential_pa, material.energy_density(&f), 3e-12);
    for i in 0..9 {
        near(response.state.first_piola_pa[i], stress[i], 2e-12);
        for j in 0..9 {
            near(
                response.tangent_pa[i][j]
                    - response.pressure_coupling[i] * response.pressure_coupling[j]
                        / response.pressure_tangent_per_pa,
                tangent[i][j],
                3e-11,
            );
        }
    }
    let remainder = material.remainder_response(&f).unwrap();
    let (remainder_stress, remainder_energy) = material.remainder_stress_and_potential(&f).unwrap();
    near(remainder.potential_pa, potential(&f, &inverse, jth, 0.0, false), 3e-12);
    near(remainder.potential_pa, remainder_energy, 3e-12);
    for i in 0..9 {
        near(remainder.first_piola_pa[i], remainder_stress[i], 2e-12);
    }
}

#[test]
fn thermal_mixed_reference_preserves_original_domain_and_nonthermal_parity() {
    let card = source_card(6.0, 1.0, false, false);
    let material = bind(&card);
    let f = [1.1, 0.1, 0.0, 0.0, 0.95, 0.0, 0.0, 0.0, 1.03];
    let native = material.law().mixed_response(&f, 2.0, VolumetricTreatment::OriginalMaterial).unwrap();
    let actual = material.mixed_response(&f, 2.0, VolumetricTreatment::OriginalMaterial).unwrap();
    assert_eq!(actual.tangent_pa, native.tangent_pa);
    assert_eq!(actual.pressure_coupling, native.pressure_coupling);
    assert_eq!(actual.pressure_tangent_per_pa, native.pressure_tangent_per_pa);
    assert_eq!(actual.state.first_piola_pa, native.state.first_piola_pa);
    assert_eq!(actual.state.potential_pa, native.state.potential_pa);

    let card = super::thermal::thermal_card(0.3, 1.0, false);
    let material = super::thermal::prepared(&card);
    let fth = IDENTITY.map(|x| 4.0 * x);
    assert!(material.law().mixed_state(&fth, 0.0, VolumetricTreatment::OriginalMaterial).is_err());
    for treatment in [
        VolumetricTreatment::OriginalMaterial,
        VolumetricTreatment::ExactIncompressible,
    ] {
        let response = material.mixed_response(&fth, 0.0, treatment).unwrap();
        let scalar = material.mixed_pressure_potential(&fth, 0.0, treatment).unwrap();
        near(scalar.volume_ratio, 64.0, 1e-12);
        near(scalar.compliance_per_pa, -response.pressure_tangent_per_pa, 2e-14);
        assert!(material.mixed_pressure_potential(&IDENTITY, 0.0, treatment).is_err());
        near(response.state.volume_ratio, 64.0, 1e-12);
        near(response.state.bulk_volume_ratio, 64.0, 1e-12);
        near(response.state.pressure_residual, 0.0, 1e-12);
        for p in response.state.first_piola_pa {
            near(p, 0.0, 1e-10);
        }
        // Domain checks belong to Fe. A valid bare total F cannot hide an
        // inadmissible elastically compressed thermal reference.
        assert!(material.mixed_state(&IDENTITY, 0.0, treatment).is_err());
        let mut inverted = fth;
        inverted[0] *= -1.0;
        assert!(material.mixed_response(&inverted, 0.0, treatment).is_err());
    }
    assert!(material.mixed_response(&fth, -BULK, VolumetricTreatment::OriginalMaterial).is_err());
    near(material.reference_density_kg_m3(), 6.0, 0.0);
    assert_eq!(material.density().card_identity(), card.content_hash());
}

#[test]
fn unequal_original_temperatures_produce_native_compatible_pressure_and_reactions() {
    with_cx(|_, cx| {
        let card = super::thermal::thermal_card(0.004, 1.0, false);
        let materials = [mean_at(&card, 310.0), mean_at(&card, 304.0)];
        let p1 = [
            [0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 1.0, 1.0],
        ];
        let cells = [[0, 1, 2, 3], [4, 1, 3, 2]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &cells, None, cx).unwrap();
        let fixed: Vec<_> = (0..3 * nodes.len()).map(|d| (d, 0.0)).collect();
        let mut problem = mixed_problem(
            &nodes, &elements, MixedLawField::PerElement(&materials),
            &fixed, VolumetricTreatment::OriginalMaterial,
        );
        problem.settings.maximum_newton_iterations = 2;
        let prepared = problem.prepare(cx).unwrap();
        let solution = prepared.solve(cx).unwrap();
        // A quadratic source gives a linear pressure equation. The native
        // all-prescribed specialization must use the same Jth-scaled scalar
        // compliance as the full mixed Hessian.
        let jth = [1.04_f64.powi(3), 1.016_f64.powi(3)];
        let volumes = [1.0 / 6.0, 1.0 / 3.0];
        // Independent P1 mass projection: He=Ve*Jth/(20K)*(I+11^T),
        // be=Ve*(1-Jth)/4. Shared-face corner pressures are equal by symmetry.
        let h = [
            volumes[0] * jth[0] / (20.0 * BULK),
            volumes[1] * jth[1] / (20.0 * BULK),
        ];
        let b = [
            volumes[0] * (1.0 - jth[0]) / 4.0,
            volumes[1] * (1.0 - jth[1]) / 4.0,
        ];
        let shared = (b[0] + b[1]) / (5.0 * (h[0] + h[1]));
        let expected = [
            0.5 * b[0] / h[0] - 1.5 * shared, shared, shared, shared,
            0.5 * b[1] / h[1] - 1.5 * shared,
        ];
        assert_eq!(solution.pressure_nodes, [0, 1, 2, 3, 4]);
        for (&actual, expected) in solution.pressure_pa.iter().zip(expected) {
            near(actual, expected, 2e-9);
        }
        assert!(solution.reactions_n.iter().any(|(_, force)| force.abs() > 1e-3));
        assert!(solution.displacement_m.iter().flatten().all(|u| *u == 0.0));
        near(solution.min_det_f, 1.0, 2e-14);
        let energy = 0.5
            * (b[0] * expected[0] + b[1] * expected[4] + 3.0 * (b[0] + b[1]) * shared);
        near(solution.potential_j, energy, 2e-10);
        let state = prepared.evaluate(
            &solution.displacement_m, &solution.pressure_pa, true, cx,
        ).unwrap();
        for residual in state.pressure_residual_m3 {
            near(residual, 0.0, 2e-11);
        }
        for (material, temperature) in materials.iter().zip([310.0, 304.0]) {
            assert_eq!(material.law().context(), &point(temperature));
            assert_eq!(material.density().query(), point(300.0));
            assert_eq!(material.density().card_identity(), card.content_hash());
            near(material.reference_density_kg_m3(), 6.0, 0.0);
        }
    });
}

#[test]
fn restrained_thermal_prestress_reaches_native_modes_with_original_density() {
    with_cx(|_, cx| {
        let card = super::thermal::thermal_card(0.01, 1.0, false);
        let material = super::thermal::prepared(&card);
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let problem = mixed_problem(
            &nodes, &elements, MixedLawField::Uniform(&material),
            &fixed, VolumetricTreatment::OriginalMaterial,
        );
        let rho = material.reference_density_kg_m3();
        let spectrum = problem.solve_vibration(&[rho], modal_settings(), cx).unwrap();
        let pressure = BULK * (1.0 / 1.1_f64.powi(3) - 1.0);
        for p in &spectrum.equilibrium().pressure_pa {
            near(*p, pressure, 2e-9);
        }
        let tangent = material.consistent_tangent(&IDENTITY).unwrap();
        let expected_k = 3.0 * VOLUME / 5.0 * tangent[0][0];
        let expected_m = rho * VOLUME / 70.0;
        near(spectrum.stiffness_n_m()[0], expected_k, 2e-9);
        near(spectrum.mass_kg()[0], expected_m, 2e-14);
        near(spectrum.modes()[0].displacement.lambda, expected_k / expected_m, 2e-7);
        assert!(spectrum.equilibrium().reactions_n.iter().any(|(_, force)| force.abs() > 0.01));
        let MixedLawField::Uniform(retained) = spectrum.problem().material else {
            panic!("source material lost")
        };
        assert!(core::ptr::eq(retained, &material));
        assert_eq!(spectrum.reference_density_kg_m3(), &[rho]);
    });
}

#[test]
fn free_thermal_volume_survives_exact_incompressible_vibration() {
    with_cx(|_, cx| {
        let card = super::thermal::thermal_card(0.01, 1.0, false);
        let material = super::thermal::prepared(&card);
        let stretch = 1.1;
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let initial: Vec<_> = nodes.iter().map(|x| x.map(|v| (stretch - 1.0) * v)).collect();
        let free = [0, 3, 7, 11];
        let fixed: Vec<_> = (0..3 * nodes.len()).filter(|d| !free.contains(d))
            .map(|d| (d, initial[d / 3][d % 3])).collect();
        let mut problem = mixed_problem(
            &nodes, &elements, MixedLawField::Uniform(&material),
            &fixed, VolumetricTreatment::ExactIncompressible,
        );
        problem.pressure_gauge = Some((0, 0.0));
        let rho = material.reference_density_kg_m3();
        let spectrum = problem.solve_incompressible_vibration_from(
            &initial, &[0.0; 4], &[rho], modal_settings(), cx,
        ).unwrap();
        assert_eq!(spectrum.dimension(), 1);
        assert_eq!(spectrum.free_dofs(), &free);
        near(spectrum.equilibrium().min_det_f, stretch.powi(3), 2e-11);
        for (actual, expected) in spectrum.equilibrium().displacement_m.iter().zip(&initial) {
            for i in 0..3 {
                near(actual[i], expected[i], 2e-10);
            }
        }
        for pressure in &spectrum.equilibrium().pressure_pa {
            near(*pressure, 0.0, 2e-9);
        }
        // At Fe=I and p=0, A0=s*Ae while reference mass is unchanged.
        // The unit-tet divergence-free corner shape [1,-1,-1,-1] has
        // K=26*mu*V/5 and M=11*rho*V/210, independently of the pressure gauge.
        let mode = &spectrum.modes()[0];
        near(mode.coordinates.lambda, stretch * 1092.0 * SHEAR / (11.0 * rho), 3e-7);
        let mut mass_norm = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                mass_norm += mode.displacement_shape[i]
                    * spectrum.mass_kg()[4 * i + j] * mode.displacement_shape[j];
            }
        }
        near(mass_norm, 1.0, 3e-11);
        assert!(mode.constraint_relative_residual <= modal_settings().pressure_residual_tolerance);
        assert!(mode.momentum_relative_residual <= modal_settings().pressure_residual_tolerance);
        assert_eq!(spectrum.reference_density_kg_m3(), &[rho]);
        let MixedLawField::Uniform(retained) = spectrum.problem().material else {
            panic!("thermal source material lost")
        };
        assert!(core::ptr::eq(retained, &material));
    });
}
