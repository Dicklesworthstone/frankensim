//! Independent narrow-beam thermoelastic oracles on native frame equilibrium.
//! Synthetic material curves establish numerical behavior, not measured fidelity.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, SolidError};
use std::sync::Arc;

#[derive(Clone, Copy)]
struct Ply {
    h: f64,
    rho: f64,
    young: f64,
    poisson: f64,
    alpha: f64,
    softening: f64,
}

fn source(plies: &[Ply], reference: ReferenceSurface) -> Arc<LaminateThermalSource> {
    let inputs = plies
        .iter()
        .enumerate()
        .map(|(index, p)| {
            let mut claims = ClaimSet::new();
            for (key, dims, value, slope) in [
                (DENSITY_PROPERTY, Density::DIMS, p.rho, 0.0),
                (
                    "young_modulus",
                    Pressure::DIMS,
                    p.young,
                    -p.young * p.softening / 100.0,
                ),
                ("poisson_ratio", Dims::NONE, p.poisson, 0.0),
                (ALPHA, INVERSE_TEMPERATURE_DIMS, p.alpha, p.alpha / 100.0),
            ] {
                // Density is intentionally unavailable at every hot temperature:
                // reference inventory, not endpoint density, must drive inertia.
                let validity = if key == DENSITY_PROPERTY {
                    ValidityDomain::unconstrained().with("T", 300.0, 300.0)
                } else {
                    ValidityDomain::unconstrained().with("T", 260.0, 400.0)
                };
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(key, dims),
                        value: if slope == 0.0 {
                            PropertyValue::Scalar { value, dims }
                        } else {
                            PropertyValue::Curve {
                                abscissa: "T".into(),
                                abscissa_dims: Temperature::DIMS,
                                knots: vec![
                                    (260.0, value - 40.0 * slope),
                                    (400.0, value + 100.0 * slope),
                                ],
                                dims,
                            }
                        },
                        interpolation: if slope == 0.0 {
                            InterpolationPolicy::ConstantWithinValidity
                        } else {
                            InterpolationPolicy::LinearInside
                        },
                        validity,
                        uncertainty: UncertaintyModel::Unstated,
                        observations: vec![],
                        provenance: Provenance {
                            source: "synthetic native laminate beam".into(),
                            license: "test-only".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
            LaminaSourceInput {
                card: Arc::new(
                    MaterialCard::assemble(
                        MaterialStateId {
                            chemistry: format!("thermal ply {index}"),
                            phase: "solid".into(),
                            process: "numerical oracle".into(),
                            revision: 0,
                        },
                        claims,
                        vec![],
                    )
                    .unwrap(),
                ),
                reference: QueryPoint::new().with("T", 300.0).unwrap(),
                thickness_m: p.h,
                orientation_rad: 0.0,
                elasticity: LaminaElasticSource::Isotropic(
                    IsotropicElasticParameters::YoungPoisson {
                        young_modulus: ScalarPropertyRequirement::try_new(
                            "young_modulus",
                            Pressure::DIMS,
                            ScalarAdmissibility::StrictlyPositive,
                        )
                        .unwrap(),
                        poisson_ratio: ScalarPropertyRequirement::try_new(
                            "poisson_ratio",
                            Dims::NONE,
                            ScalarAdmissibility::Finite,
                        )
                        .unwrap(),
                    },
                ),
                expansion: LaminaExpansionSource::Isotropic,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        })
        .collect();
    Arc::new(LaminateThermalSource::try_new(inputs, reference, || false).unwrap())
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * scale.abs().max(expected.abs()).max(1e-30),
        "{actual:e} != {expected:e} (scale {scale:e})"
    );
}
fn free(p: Ply, t: f64) -> f64 {
    let dt = t - 300.0;
    p.alpha * (dt + dt * dt / 200.0)
}
fn young(p: Ply, t: f64) -> f64 {
    p.young * (1.0 - p.softening * (t - 300.0) / 100.0)
}
fn member(
    source: Arc<LaminateThermalSource>,
    width: f64,
    length: f64,
    angle: f64,
    temperatures: &[Vec<f64>],
    cx: &Cx<'_>,
) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(length * angle.cos(), length * angle.sin());
    frame.fix(0, [true; 3]);
    frame
        .add_laminate_thermal_element(0, 1, source, width, temperatures, 0.02, cx)
        .unwrap();
    frame
}
fn bimetal() -> [Ply; 2] {
    [
        Ply {
            h: 0.003,
            rho: 2700.0,
            young: 2e9,
            poisson: 0.19,
            alpha: 1e-5,
            softening: 0.1,
        },
        Ply {
            h: 0.002,
            rho: 7800.0,
            young: 5e9,
            poisson: 0.37,
            alpha: 2.5e-5,
            softening: 0.3,
        },
    ]
}

// Independent exact layer integrals about the reference density centroid.
fn oracle(
    plies: &[Ply],
    temperatures: &[f64],
    width: f64,
) -> ([[f64; 2]; 2], [f64; 2], Vec<[f64; 2]>) {
    let mut z = 0.0;
    let mut mass = 0.0;
    let mut first_mass = 0.0;
    for p in plies {
        mass += p.rho * p.h;
        first_mass += p.rho * p.h * (z + 0.5 * p.h);
        z += p.h;
    }
    let center = first_mass / mass;
    let mut k = [[0.0; 2]; 2];
    let mut thermal = [0.0; 2];
    let mut faces = Vec::new();
    z = -center;
    for (&p, &t) in plies.iter().zip(temperatures) {
        let mid = z + 0.5 * p.h;
        let a = width * young(p, t) * p.h;
        k[0][0] += a;
        k[0][1] -= a * mid;
        k[1][1] += a * (mid * mid + p.h * p.h / 12.0);
        thermal[0] += a * free(p, t);
        thermal[1] -= a * mid * free(p, t);
        faces.push([z, z + p.h]);
        z += p.h;
    }
    k[1][0] = k[0][1];
    (k, thermal, faces)
}
fn solve(k: [[f64; 2]; 2], force: [f64; 2]) -> [f64; 2] {
    let det = k[0][0] * k[1][1] - k[0][1] * k[1][0];
    [
        (k[1][1] * force[0] - k[0][1] * force[1]) / det,
        (k[0][0] * force[1] - k[1][0] * force[0]) / det,
    ]
}

#[test]
fn heterogeneous_poisson_bimetal_bends_with_coupling_and_residual_layer_storage() {
    with_cx(|cx| {
        let plies = bimetal();
        let temperatures = [340.0, 360.0];
        let width = 0.0005;
        let length = 0.1;
        let angle: f64 = 0.63;
        let original = source(&plies, ReferenceSurface::Offset(0.011));
        let mut frame = member(
            Arc::clone(&original),
            width,
            length,
            angle,
            &vec![temperatures.to_vec(); 3],
            cx,
        );
        frame.static_load(&[0.0; 6], 1).unwrap();
        let (k, thermal, faces) = oracle(&plies, &temperatures, width);
        let [eps, curvature] = solve(k, thermal);
        assert!(k[0][1].abs() > 1.0);
        assert!(curvature < 0.0);
        let axial = eps * length;
        let transverse = 0.5 * curvature * length * length;
        near(
            frame.displacements()[3],
            axial * angle.cos() - transverse * angle.sin(),
            0.001,
        );
        near(
            frame.displacements()[4],
            axial * angle.sin() + transverse * angle.cos(),
            0.001,
        );
        near(frame.displacements()[5], curvature * length, 0.01);
        let response = frame.laminate_thermal_response(0, 1).unwrap();
        near(response.axial_strain, eps, 0.001);
        near(response.curvature_per_m, curvature, 1.0);
        near(response.axial_force_n, 0.0, thermal[0]);
        near(response.bending_moment_n_m, 0.0, thermal[1]);
        let mut energy = 0.0;
        for (layer, ((&p, &t), &[bottom, top])) in
            plies.iter().zip(&temperatures).zip(&faces).enumerate()
        {
            let mid = 0.5 * (bottom + top);
            let elastic = eps - mid * curvature - free(p, t);
            energy += 0.5
                * width
                * young(p, t)
                * (p.h * elastic * elastic + curvature * curvature * p.h.powi(3) / 12.0);
            for (y, stress) in [bottom, top].into_iter().zip([
                response.layer_stress_pa[layer].bottom_pa,
                response.layer_stress_pa[layer].top_pa,
            ]) {
                near(
                    stress[0],
                    young(p, t) * (eps - y * curvature - free(p, t)),
                    1e6,
                );
                near(stress[1], 0.0, 1e6);
                near(stress[2], 0.0, 1e6);
            }
        }
        near(response.strain_energy_j_m, energy, energy);
        assert!(energy > 0.0);
        let bound = frame.laminate_thermal_section(0, 1).unwrap();
        for row in 0..2 {
            for col in 0..2 {
                near(
                    bound.prepared().tangent()[row][col],
                    k[row][col],
                    k[row][col],
                );
            }
        }
        assert!(Arc::ptr_eq(
            &bound.source().source().layers()[0].input().card,
            &original.layers()[0].input().card
        ));
    });
}

#[test]
fn full_restraint_recovers_actual_thermal_force_moment_and_layer_stresses() {
    with_cx(|cx| {
        let plies = bimetal();
        let temperatures = [350.0, 370.0];
        let width = 0.0005;
        let cold = vec![vec![300.0; 2]; 3];
        let mut frame = member(
            source(&plies, ReferenceSurface::Midsurface),
            width,
            0.1,
            0.0,
            &cold,
            cx,
        );
        frame.fix(1, [true; 3]);
        let mut hot = frame
            .with_laminate_temperatures(0, &vec![temperatures.to_vec(); 3], cx)
            .unwrap();
        hot.static_load(&[0.0; 6], 1).unwrap();
        let (_, thermal, _) = oracle(&plies, &temperatures, width);
        let response = hot.laminate_thermal_response(0, 1).unwrap();
        near(response.axial_force_n, -thermal[0], thermal[0]);
        near(response.bending_moment_n_m, -thermal[1], thermal[1]);
        near(hot.element_forces(0)[1], thermal[1], thermal[1]);
        near(hot.element_forces(0)[2], -thermal[1], thermal[1]);
        assert_eq!(hot.displacements(), &[0.0; 6]);
        let mut energy = 0.0;
        for (layer, (&p, &t)) in plies.iter().zip(&temperatures).enumerate() {
            let expected = -young(p, t) * free(p, t);
            near(
                response.layer_stress_pa[layer].bottom_pa[0],
                expected,
                expected,
            );
            near(
                response.layer_stress_pa[layer].top_pa[0],
                expected,
                expected,
            );
            energy += 0.5 * width * young(p, t) * p.h * free(p, t).powi(2);
        }
        near(response.strain_energy_j_m, energy, energy);
    });
}

#[test]
fn spatial_thermal_gradient_satisfies_native_force_equilibrium_and_compatibility() {
    with_cx(|cx| {
        let plies = bimetal();
        let width = 0.0005;
        let length = 0.1;
        let temperatures = vec![vec![310.0, 330.0], vec![330.0, 360.0], vec![350.0, 380.0]];
        let mut frame = member(
            source(&plies, ReferenceSurface::Midsurface),
            width,
            length,
            0.0,
            &temperatures,
            cx,
        );
        frame.static_load(&[0.0; 6], 1).unwrap();
        let mut extension = 0.0;
        let mut rotation = 0.0;
        let mut deflection = 0.0;
        for (station, ((t, weight), xi)) in temperatures
            .iter()
            .zip([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0])
            .zip([0.0, 0.5, 1.0])
            .enumerate()
        {
            let (k, thermal, _) = oracle(&plies, t, width);
            let [eps, curvature] = solve(k, thermal);
            let response = frame.laminate_thermal_response(0, station).unwrap();
            near(response.axial_strain, eps, 0.001);
            near(response.curvature_per_m, curvature, 1.0);
            near(response.axial_force_n, 0.0, thermal[0]);
            near(response.bending_moment_n_m, 0.0, thermal[1]);
            extension += weight * length * eps;
            rotation += weight * length * curvature;
            deflection += weight * length * length * (1.0 - xi) * curvature;
        }
        near(frame.displacements()[3], extension, extension);
        near(frame.displacements()[4], deflection, deflection);
        near(frame.displacements()[5], rotation, rotation);
    });
}

#[test]
fn temperature_refresh_keeps_source_reference_mass_spin_and_exact_current_modulus() {
    with_cx(|cx| {
        let ply = Ply {
            h: 0.005,
            rho: 4300.0,
            young: 2e9,
            poisson: 0.31,
            alpha: 1e-5,
            softening: 0.5,
        };
        let width = 0.0005;
        let original_source = source(&[ply], ReferenceSurface::Offset(0.023));
        let original = member(
            Arc::clone(&original_source),
            width,
            0.1,
            0.0,
            &vec![vec![300.0]; 3],
            cx,
        );
        let acceleration = [0.1, -9.81, 0.7, 0.3, -9.81, -0.2];
        let original_mass = original.inertial_forces(&acceleration).unwrap();
        let original_frequencies = original.natural_frequencies().unwrap();
        let [mass, first, spin] = original
            .laminate_thermal_section(0, 0)
            .unwrap()
            .prepared()
            .reference_mass_moments();
        near(mass, width * ply.rho * ply.h, mass);
        near(first, 0.0, mass * ply.h);
        near(spin, width * ply.rho * ply.h.powi(3) / 12.0, spin);
        let mut current = original.clone();
        for t in [350.0, 280.0, 370.0, 300.0] {
            current = current
                .with_laminate_temperatures(0, &vec![vec![t]; 3], cx)
                .unwrap();
            assert_eq!(
                current.inertial_forces(&acceleration).unwrap(),
                original_mass
            );
            current.static_load(&[0.0; 6], 1).unwrap();
            near(current.displacements()[3], 0.1 * free(ply, t), 0.001);
            near(current.displacements()[4], 0.0, 0.001);
            for (&before, after) in original_frequencies
                .iter()
                .zip(current.natural_frequencies().unwrap())
            {
                near(after, before * (young(ply, t) / ply.young).sqrt(), before);
            }
            let bound = current.laminate_thermal_section(0, 1).unwrap();
            assert_eq!(
                bound.source().source().layers()[0].reference_material(),
                original_source.layers()[0].reference_material()
            );
        }
        assert_eq!(original.displacements(), &[0.0; 6]);
    });
}

#[test]
fn thin_thermal_sections_keep_physical_flexibility_without_ea_based_floor() {
    with_cx(|cx| {
        for h in [0.01, 1e-8] {
            let ply = Ply {
                h,
                rho: 4000.0,
                young: 2e9,
                poisson: 0.31,
                alpha: 0.0,
                softening: 0.0,
            };
            let width = 10.0 * h;
            let length = 1000.0 * h;
            let mut frame = member(
                source(&[ply], ReferenceSurface::Midsurface),
                width,
                length,
                0.0,
                &vec![vec![300.0]; 3],
                cx,
            );
            // Constrain axial motion so a global LU tolerance comparing unlike
            // axial and rotational units is not this sectional regression.
            frame.fix(1, [true, false, false]);
            let target = 0.01 * length;
            let mut pattern = [0.0; 6];
            pattern[4] = 1.0;
            let steps = frame.displacement_control(&pattern, 4, &[target]).unwrap();
            let ei = ply.young * width * h.powi(3) / 12.0;
            let force = 3.0 * ei * target / length.powi(3);
            near(steps[0].lambda, force, force);
            near(frame.displacements()[5], 1.5 * target / length, 0.015);
            near(
                frame
                    .laminate_thermal_section(0, 0)
                    .unwrap()
                    .prepared()
                    .tangent()[1][1],
                ei,
                ei,
            );
        }
    });
}

#[test]
fn source_layer_refinement_converges_to_exact_affine_gradient_curvature() {
    with_cx(|cx| {
        let h = 0.005;
        let length = 0.1;
        let gradient = 0.05;
        let center_free = 0.0005;
        let mut previous_error = None;
        let mut previous_energy = f64::INFINITY;
        for count in [2, 4, 8, 16] {
            let ply = Ply {
                h: h / count as f64,
                rho: 4000.0,
                young: 2e9,
                poisson: 0.31,
                alpha: 1e-5,
                softening: 0.0,
            };
            let plies = vec![ply; count];
            // Invert the known instantaneous-alpha integral so the physical
            // free strain, rather than alpha(T)*DeltaT, is affine in y.
            let temperatures: Vec<f64> = (0..count)
                .map(|layer| {
                    let y = -0.5 * h + (layer as f64 + 0.5) * ply.h;
                    let target = center_free + gradient * y;
                    200.0 + (10000.0 + 200.0 * target / ply.alpha).sqrt()
                })
                .collect();
            let mut frame = member(
                source(&plies, ReferenceSurface::Midsurface),
                0.0005,
                length,
                0.0,
                &vec![temperatures; 3],
                cx,
            );
            frame.static_load(&[0.0; 6], 1).unwrap();
            let response = frame.laminate_thermal_response(0, 1).unwrap();
            let expected = -gradient * (1.0 - 1.0 / (count * count) as f64);
            near(response.curvature_per_m, expected, gradient);
            near(
                frame.displacements()[4],
                0.5 * expected * length * length,
                0.001,
            );
            near(response.axial_strain, center_free, center_free);
            let error = (response.curvature_per_m + gradient).abs();
            if let Some(previous) = previous_error {
                near(previous / error, 4.0, 4.0);
            }
            assert!(response.strain_energy_j_m < previous_energy);
            previous_error = Some(error);
            previous_energy = response.strain_energy_j_m;
        }
    });
}

#[test]
fn invalid_temperature_bounds_and_cancellation_publish_no_partial_member_or_refresh() {
    with_cx(|cx| {
        let plies = bimetal();
        let src = source(&plies, ReferenceSurface::Midsurface);
        let original = member(
            Arc::clone(&src),
            0.0005,
            0.1,
            0.0,
            &vec![vec![300.0; 2]; 3],
            cx,
        );
        for temperatures in [
            vec![vec![350.0; 2]; 2],
            vec![vec![401.0; 2]; 3],
            vec![vec![350.0]; 3],
        ] {
            assert!(
                original
                    .with_laminate_temperatures(0, &temperatures, cx)
                    .is_err()
            );
        }
        assert_eq!(original.displacements(), &[0.0; 6]);
        assert_eq!(original.element_forces(0), [0.0; 3]);
        let mut empty = Frame2d::new(Geometry::Linear);
        empty.add_node(0.0, 0.0);
        empty.add_node(0.1, 0.0);
        assert!(
            empty
                .add_laminate_thermal_element(
                    0,
                    1,
                    Arc::clone(&src),
                    0.0005,
                    &vec![vec![390.0; 2]; 3],
                    1e-6,
                    cx
                )
                .is_err()
        );
        assert!(empty.laminate_thermal_section(0, 0).is_none());
        assert_eq!(empty.inertial_forces(&[1.0; 6]).unwrap(), vec![0.0; 6]);

        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 21,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                original.with_laminate_temperatures(0, &vec![vec![350.0; 2]; 3], &cancelled),
                Err(SolidError::Cancelled)
            ));
            assert!(matches!(
                empty.add_laminate_thermal_element(
                    0,
                    1,
                    src,
                    0.0005,
                    &vec![vec![350.0; 2]; 3],
                    0.02,
                    &cancelled
                ),
                Err(SolidError::Cancelled)
            ));
        });
        assert!(empty.laminate_thermal_section(0, 0).is_none());
    });
}
