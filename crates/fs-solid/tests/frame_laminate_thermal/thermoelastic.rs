//! Native accepted-state storage and reversible heat, with independent oracles.

use super::*;
use fs_material::laminate::strength::{LaminaStrengthProperties, LaminaStrengthSelection};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

#[path = "thermoelastic/motion.rs"]
mod motion;

#[path = "thermoelastic/dynamics.rs"]
mod dynamics;

fn kelvin() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

// Authored test cards explicitly declare absolute temperature. The production
// attachment must never add that meaning to an existing legacy source.
fn typed_source(plies: &[Ply], angle: f64, strengths: bool) -> Arc<LaminateThermalSource> {
    let original = source(plies, ReferenceSurface::Offset(0.011));
    let inputs = original
        .layers()
        .iter()
        .map(|layer| {
            let mut input = layer.input().clone();
            let mut claims = ClaimSet::new();
            for (_, old) in input.card.claims().claims_ordered() {
                let mut claim = old.clone();
                let &(lo, hi) = claim.validity.bounds().get("T").unwrap();
                claim.validity = claim.validity.with_quantity("T", kelvin(), lo, hi);
                claims.insert_claim(claim).unwrap();
            }
            if strengths {
                for name in ["tension", "compression", "shear"] {
                    claims
                        .insert_claim(PropertyClaim {
                            key: PropertyKey::new(name, Pressure::DIMS),
                            value: PropertyValue::Curve {
                                abscissa: "T".into(),
                                abscissa_dims: Temperature::DIMS,
                                knots: vec![(260.0, 12e6), (400.0, 5e6)],
                                dims: Pressure::DIMS,
                            },
                            interpolation: InterpolationPolicy::LinearInside,
                            validity: ValidityDomain::unconstrained().with_quantity(
                                "T",
                                kelvin(),
                                260.0,
                                400.0,
                            ),
                            uncertainty: UncertaintyModel::Unstated,
                            observations: vec![],
                            provenance: Provenance {
                                source: "synthetic thermal strength".into(),
                                license: "test-only".into(),
                                artifact: None,
                            },
                        })
                        .unwrap();
                }
            }
            input.card =
                Arc::new(MaterialCard::assemble(input.card.id().clone(), claims, vec![]).unwrap());
            input.reference = QueryPoint::new()
                .with_quantity("T", kelvin(), 300.0)
                .unwrap();
            input.orientation_rad = angle;
            input
        })
        .collect();
    Arc::new(
        LaminateThermalSource::try_new(inputs, ReferenceSurface::Offset(0.011), || false).unwrap(),
    )
}

fn strength_selection() -> LaminaStrengthSelection {
    LaminaStrengthSelection {
        properties: LaminaStrengthProperties {
            tensile: ["tension".into(), "tension".into()],
            compressive: ["compression".into(), "compression".into()],
            shear: "shear".into(),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

// Integrate the independent one-dimensional potential E(T)*(eps-y*kappa-f(T))²/2.
// Returns N, M, Psi, S, U, C, N_T, M_T per original member length.
fn potential(p: Ply, t: f64, width: f64, mid: f64, eps: f64, kappa: f64) -> [f64; 8] {
    let area = width * p.h;
    let e = young(p, t);
    let e_t = -p.young * p.softening / 100.0;
    let alpha = p.alpha * (1.0 + (t - 300.0) / 100.0);
    let alpha_t = p.alpha / 100.0;
    let elastic = eps - mid * kappa - free(p, t);
    let square = elastic * elastic + kappa * kappa * p.h * p.h / 12.0;
    let moment = mid * elastic - kappa * p.h * p.h / 12.0;
    let psi = 0.5 * area * e * square;
    let entropy = area * (e * alpha * elastic - 0.5 * e_t * square);
    let capacity = t * area * ((2.0 * e_t * alpha + e * alpha_t) * elastic - e * alpha * alpha);
    [
        area * e * elastic,
        -area * e * moment,
        psi,
        entropy,
        psi + t * entropy,
        capacity,
        area * (e_t * elastic - e * alpha),
        -area * (e_t * moment - e * alpha * mid),
    ]
}

#[test]
fn actual_bimetal_stations_integrate_source_energy_entropy_and_capacity() {
    with_cx(|cx| {
        let plies = bimetal();
        let temperatures = vec![vec![310.0, 330.0], vec![330.0, 360.0], vec![350.0, 380.0]];
        let width = 0.0005;
        let length = 0.1;
        let mut frame = member(
            typed_source(&plies, 0.0, false),
            width,
            length,
            0.61,
            &temperatures,
            cx,
        )
        .with_laminate_thermoelasticity(0, &[[260.0, 400.0]; 2], cx)
        .unwrap();
        frame.static_load(&[0.0; 6], 1).unwrap();
        let mut total = [0.0; 4];
        let mut expected_heat = 0.0;
        for (station, (&weight, t)) in [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]
            .iter()
            .zip(&temperatures)
            .enumerate()
        {
            let response = frame
                .laminate_thermoelastic_response(0, station, cx)
                .unwrap();
            let (k, thermal, faces) = oracle(&plies, t, width);
            let [eps, kappa] = solve(k, thermal);
            near(response.axial_strain, eps, 0.001);
            near(response.curvature_per_m, kappa, 1.0);
            near(response.reference_length_weight_m, weight * length, length);
            let rate = [1e-4 * (station + 1) as f64, -0.2 + 0.1 * station as f64];
            let t_rate = [2.0 + station as f64, -1.0];
            let heat = frame
                .laminate_thermoelastic_heat_rate(0, station, rate, &t_rate, cx)
                .unwrap();
            for layer in 0..2 {
                let mid = 0.5 * (faces[layer][0] + faces[layer][1]);
                let expected = potential(plies[layer], t[layer], width, mid, eps, kappa);
                let actual = response.layers[layer].source;
                for (a, b) in [
                    actual.axial_force_n,
                    actual.bending_moment_n_m,
                    actual.free_energy_j_m,
                    actual.entropy_j_m_k,
                    actual.internal_energy_j_m,
                    actual.heat_capacity_j_m_k,
                    actual.generalized_force_temperature_derivative[0],
                    actual.generalized_force_temperature_derivative[1],
                ]
                .into_iter()
                .zip(expected)
                {
                    near(a, b, b.abs().max(1e-10));
                }
                assert!(
                    response.layers[layer].minimum_elastic_heat_capacity_j_m_k
                        <= actual.heat_capacity_j_m_k
                );
                for axis in 0..4 {
                    total[axis] += length * weight * expected[axis + 2];
                }
                let deform = t[layer] * (expected[6] * rate[0] + expected[7] * rate[1]);
                let thermal = -expected[5] * t_rate[layer];
                near(
                    heat.layers[layer].deformation_heat_w_m,
                    deform,
                    deform.abs().max(1e-10),
                );
                near(
                    heat.layers[layer].temperature_heat_w_m,
                    thermal,
                    thermal.abs().max(1e-10),
                );
                expected_heat += length * weight * (deform + thermal);
            }
        }
        let energy = frame.laminate_thermoelastic_energy(0, cx).unwrap();
        for (actual, expected) in [
            energy.free_energy_j,
            energy.entropy_j_k,
            energy.internal_energy_j,
            energy.elastic_heat_capacity_j_k,
        ]
        .into_iter()
        .zip(total)
        {
            near(actual, expected, expected.abs());
        }
        assert!(energy.free_energy_j > 0.0);
        let actual_heat: f64 = (0..3)
            .map(|station| {
                let heat = frame
                    .laminate_thermoelastic_heat_rate(
                        0,
                        station,
                        [1e-4 * (station + 1) as f64, -0.2 + 0.1 * station as f64],
                        &[2.0 + station as f64, -1.0],
                        cx,
                    )
                    .unwrap();
                heat.reference_length_weight_m
                    * heat
                        .layers
                        .iter()
                        .map(|p| p.reversible_heat_w_m)
                        .sum::<f64>()
            })
            .sum();
        near(actual_heat, expected_heat, expected_heat.abs());
    });
}

#[test]
fn free_expansion_heat_terms_cancel_but_fixed_member_heating_changes_elastic_storage() {
    with_cx(|cx| {
        for sign in [1.0, -1.0] {
            let p = Ply {
                h: 0.005,
                rho: 4300.0,
                young: 2e9,
                poisson: 0.31,
                alpha: sign * 1e-5,
                softening: 0.4,
            };
            let width = 0.0005;
            let t = 350.0;
            let mut free_frame = member(
                typed_source(&[p], 0.0, false),
                width,
                0.1,
                0.0,
                &vec![vec![t]; 3],
                cx,
            )
            .with_laminate_thermoelasticity(0, &[[260.0, 400.0]], cx)
            .unwrap();
            free_frame.static_load(&[0.0; 6], 1).unwrap();
            let alpha = p.alpha * 1.5;
            for temperature_rate in [4.0, -3.0] {
                let heat = free_frame
                    .laminate_thermoelastic_heat_rate(
                        0,
                        1,
                        [alpha * temperature_rate, 0.0],
                        &[temperature_rate],
                        cx,
                    )
                    .unwrap();
                let layer = heat.layers[0];
                assert!(layer.elastic_heat_capacity_j_m_k < 0.0);
                assert!(layer.deformation_heat_w_m.abs() > 0.0);
                near(layer.reversible_heat_w_m, 0.0, layer.deformation_heat_w_m);
                let expected = -t * width * p.h * young(p, t) * alpha * alpha;
                near(layer.elastic_heat_capacity_j_m_k, expected, expected);
            }
            let mut fixed = member(
                typed_source(&[p], 0.0, false),
                width,
                0.1,
                0.0,
                &vec![vec![300.0]; 3],
                cx,
            );
            fixed.fix(1, [true; 3]);
            let fixed = fixed
                .with_laminate_thermoelasticity(0, &[[260.0, 400.0]], cx)
                .unwrap()
                .with_laminate_temperatures(0, &vec![vec![t]; 3], cx)
                .unwrap();
            let expected = potential(p, t, width, 0.0, 0.0, 0.0);
            let heat = fixed
                .laminate_thermoelastic_heat_rate(0, 1, [0.0; 2], &[4.0], cx)
                .unwrap();
            assert_eq!(heat.layers[0].deformation_heat_w_m, 0.0);
            near(
                heat.layers[0].temperature_heat_w_m,
                -4.0 * expected[5],
                expected[5],
            );
            assert_eq!(fixed.displacements(), &[0.0; 6]);
        }
    });
}

#[test]
fn source_potentials_and_strength_refresh_together_without_rebinding_or_remassing() {
    with_cx(|cx| {
        let plies = bimetal();
        let source = typed_source(&plies, 0.0, true);
        let original = member(
            Arc::clone(&source),
            0.0005,
            0.1,
            0.0,
            &vec![vec![300.0; 2]; 3],
            cx,
        );
        assert!(original.laminate_thermoelastic_energy(0, cx).is_err());
        let selected = vec![strength_selection(); 2];
        let ranges = [[280.0, 370.0]; 2];
        for thermal_first in [false, true] {
            let bound = if thermal_first {
                original
                    .with_laminate_thermoelasticity(0, &ranges, cx)
                    .unwrap()
                    .with_laminate_strength(0, &selected, cx)
                    .unwrap()
            } else {
                original
                    .with_laminate_strength(0, &selected, cx)
                    .unwrap()
                    .with_laminate_thermoelasticity(0, &ranges, cx)
                    .unwrap()
            };
            assert!(
                bound
                    .with_laminate_thermoelasticity(0, &ranges, cx)
                    .is_err()
            );
            let mut hot = bound
                .with_laminate_temperatures(0, &vec![vec![350.0, 360.0]; 3], cx)
                .unwrap();
            assert_eq!(
                hot.inertial_forces(&[1.0; 6]).unwrap(),
                original.inertial_forces(&[1.0; 6]).unwrap()
            );
            hot.static_load(&[0.0; 6], 1).unwrap();
            assert!(
                hot.assess_laminate_element_strength(0, cx)
                    .unwrap()
                    .utilization()
                    > 0.0
            );
            for station in 0..3 {
                let section = hot.laminate_thermal_section(0, station).unwrap();
                let potential = section.thermoelastic().unwrap();
                assert!(std::ptr::eq(
                    potential.source().prepared().as_ref(),
                    section.prepared().source()
                ));
                assert!(std::ptr::eq(
                    section.strength().unwrap().thermal().prepared().as_ref(),
                    section.prepared().source()
                ));
                assert!(Arc::ptr_eq(
                    &potential.source().source().layers()[0].input().card,
                    &source.layers()[0].input().card
                ));
                assert_eq!(potential.layer_temperature_k(1), Some(360.0));
            }
            let before = hot.laminate_thermoelastic_energy(0, cx).unwrap();
            let copied = hot.laminate_thermal_section(0, 0).unwrap().clone();
            let mut raw = Frame2d::new(Geometry::Linear);
            raw.add_node(0.0, 0.0);
            raw.add_node(0.1, 0.0);
            raw.add_element(0, 1, 3, 9.0, &|| {
                fs_solid::SectionModel::LaminateThermal(copied.clone())
            })
            .unwrap();
            assert!(raw.laminate_thermoelastic_response(0, 0, cx).is_err());
            assert!(raw.laminate_thermoelastic_energy(0, cx).is_err());
            assert!(
                raw.laminate_thermoelastic_heat_rate(0, 0, [0.0; 2], &[0.0; 2], cx)
                    .is_err()
            );
            assert!(
                hot.with_laminate_temperatures(0, &vec![vec![380.0; 2]; 3], cx)
                    .is_err()
            );
            assert_eq!(hot.laminate_thermoelastic_energy(0, cx).unwrap(), before);
            assert!(
                hot.laminate_thermoelastic_heat_rate(0, 1, [f64::NAN, 0.0], &[0.0; 2], cx)
                    .is_err()
            );
            assert!(
                hot.laminate_thermoelastic_heat_rate(0, 1, [0.0; 2], &[0.0], cx)
                    .is_err()
            );
        }
        assert!(
            original
                .with_laminate_thermoelasticity(0, &ranges[..1], cx)
                .is_err()
        );
        let legacy = member(
            super::source(&plies, ReferenceSurface::Midsurface),
            0.0005,
            0.1,
            0.0,
            &vec![vec![300.0; 2]; 3],
            cx,
        );
        assert!(
            legacy
                .with_laminate_thermoelasticity(0, &ranges, cx)
                .is_err()
        );
        assert!(
            original
                .laminate_thermal_section(0, 0)
                .unwrap()
                .thermoelastic()
                .is_none()
        );
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
                original.with_laminate_thermoelasticity(0, &ranges, &cancelled),
                Err(SolidError::Cancelled)
            ));
        });
    });
}

#[test]
fn rotated_isotropic_principal_engineering_shear_retains_native_tensor_envelope() {
    with_cx(|cx| {
        let p = Ply {
            h: 0.005,
            rho: 4300.0,
            young: 2e9,
            poisson: 0.31,
            alpha: 0.0,
            softening: 0.0,
        };
        let width = 0.0005;
        let length = 0.1;
        let eps: f64 = 0.018;
        let mut frame = member(
            typed_source(&[p], std::f64::consts::FRAC_PI_4, false),
            width,
            length,
            0.0,
            &vec![vec![300.0]; 3],
            cx,
        );
        let mut load = [0.0; 6];
        load[3] = 1.0;
        frame
            .displacement_control(&load, 3, &[eps * length])
            .unwrap();
        assert!(eps.hypot(p.poisson * eps) < 0.02);
        assert!((1.0 + p.poisson) * eps > 0.02);
        let bound = frame
            .with_laminate_thermoelasticity(0, &[[280.0, 370.0]], cx)
            .unwrap();
        let energy = bound.laminate_thermoelastic_energy(0, cx).unwrap();
        near(
            energy.free_energy_j,
            0.5 * width * p.h * length * p.young * eps * eps,
            energy.free_energy_j,
        );
        near(energy.entropy_j_k, 0.0, 1e-10);
    });
}
