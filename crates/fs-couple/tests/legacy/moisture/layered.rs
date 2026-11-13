//! Unlike-material equilibrium, conserved partition and nonlinear dry uptake.
//! These authored laws verify implementation, not measured material behavior.
use super::*;
use fs_couple::moisture::layered::*;
use fs_material::state_point::{
    MonotoneMoistureSorption, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

#[path = "layered/adaptive.rs"]
mod adaptive;
#[path = "layered/atmosphere.rs"]
mod atmosphere;
#[path = "layered/radial.rs"]
mod radial;
#[path = "layered/sdirk.rs"]
mod sdirk;
#[path = "layered/solution.rs"]
mod solution;
#[path = "layered/transfer.rs"]
mod transfer;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}
fn material(rho: f64, sorption: Vec<(f64, f64)>, diffusion: [f64; 2]) -> MaterialCard {
    material_profile(
        rho,
        sorption,
        vec![(0.0, diffusion[0]), (0.6, diffusion[1])],
    )
}
fn material_profile(
    rho: f64,
    sorption: Vec<(f64, f64)>,
    diffusion: Vec<(f64, f64)>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, axis, axis_quantity, knots, lo, hi) in [
        (
            PropertyKey::new("density", Density::DIMS),
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            vec![(0.0, 1.1 * rho), (0.6, 1.1 * rho)],
            0.1,
            0.1,
        ),
        (
            PropertyKey::new("moisture_diffusivity", MOISTURE_DIFFUSIVITY_DIMS),
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            diffusion,
            0.0,
            0.6,
        ),
        (
            PropertyKey::with_quantity(
                "sorption",
                qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ),
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            sorption,
            0.0,
            1.0,
        ),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: Dims::NONE,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        qty(QuantityKind::AbsoluteTemperature),
                        293.15,
                        293.15,
                    )
                    .with_quantity(axis, axis_quantity, lo, hi),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored heterogeneous moisture verification".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "verification sorbing solid".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn isotherm(card: &MaterialCard) -> MonotoneMoistureSorption {
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    isotherm_at(card, point)
}
fn isotherm_at(card: &MaterialCard, point: QueryPoint) -> MonotoneMoistureSorption {
    let source = resolve_material_state_point(
        card,
        &point,
        &[ScalarPropertyRequirement::try_with_quantity(
            "sorption",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MonotoneMoistureSorption::from_binding(card, &source, "sorption", [0.0, 1.0]).unwrap()
}
fn layer(card: MaterialCard, length: f64, initial: Vec<f64>, cx: &Cx<'_>) -> MoistureLayer {
    let sorption = isotherm(&card);
    let slab = MoistureSlab::new_kirchhoff(
        card,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureSlabGeometry {
            length_m: length,
            cross_section_m2: 0.001,
            maximum_moisture: 0.6,
            maximum_cells: 256,
        },
        initial,
        nonlinear(),
        cx,
    )
    .unwrap();
    MoistureLayer::new(slab, sorption, cx).unwrap()
}
fn sealed() -> [MoistureLayeredFace; 2] {
    [MoistureLayeredFace::Sealed; 2]
}

#[test]
fn g0_layer_sorption_and_diffusion_hold_identical_external_source_context() {
    with_cx(|cx| {
        let card = material(500.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]);
        let dry_point = point(0.1)
            .with_quantity("pressure", qty(QuantityKind::Pressure), 100_000.0)
            .unwrap();
        let slab = MoistureSlab::new_kirchhoff(
            card.clone(),
            dry_point,
            "moisture_diffusivity",
            MaterialPropertySelection::SingleClaimOnly,
            MoistureSlabGeometry {
                length_m: 0.2,
                cross_section_m2: 0.001,
                maximum_moisture: 0.6,
                maximum_cells: 256,
            },
            vec![0.1],
            nonlinear(),
            cx,
        )
        .unwrap();
        let base = isotherm(&card).source_point().clone();
        let matching = base
            .clone()
            .with_quantity("pressure", qty(QuantityKind::Pressure), 100_000.0)
            .unwrap();
        assert!(MoistureLayer::new(slab.clone(), isotherm_at(&card, matching.clone()), cx).is_ok());
        for point in [
            matching.clone().with("pressure", 200_000.0).unwrap(),
            // An equal number without the original pressure quantity is not
            // an equal source context; neither is an omitted held coordinate.
            base.clone().with("pressure", 100_000.0).unwrap(),
            base,
            // The other owner's evolving coordinate must not become a hidden
            // held condition on the isotherm side of this material join.
            matching
                .with_quantity(
                    "moisture_content",
                    qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
                    0.1,
                )
                .unwrap(),
        ] {
            assert!(matches!(
                MoistureLayer::new(slab.clone(), isotherm_at(&card, point), cx),
                Err(MoistureError::Invalid(
                    "layer sorption and transport must hold identical external coordinates and quantities"
                ))
            ));
        }
        assert_eq!(slab.moisture(), &[0.1]);
        assert_eq!(slab.time_s(), 0.0);
    });
}

#[test]
fn g1_unlike_materials_keep_different_moisture_at_exact_stiff_activity_equilibrium() {
    with_cx(|cx| {
        let left = layer(
            material(800.0, vec![(0.0, 0.0), (1.0, 0.1)], [1e100; 2]),
            0.2,
            vec![0.05; 2],
            cx,
        );
        let right = layer(
            material(300.0, vec![(0.0, 0.0), (1.0, 0.4)], [1e100; 2]),
            0.3,
            vec![0.2; 3],
            cx,
        );
        let old = MoistureLayeredSlab::new(vec![left, right], 256, cx).unwrap();
        let mut step = config();
        step.linear.max_iterations = 0;
        let (next, report) = old.advance(1e50, &sealed(), step, nonlinear(), cx).unwrap();
        assert_eq!(next.layers()[0].slab().moisture(), &[0.05; 2]);
        assert_eq!(next.layers()[1].slab().moisture(), &[0.2; 3]);
        assert_eq!(next.activity(), &[0.5; 5]);
        assert_eq!(report.boundary_water_kg, [0.0; 2]);
        assert_eq!(report.interfaces[0].water_left_to_right_kg, 0.0);
        assert_eq!(report.interfaces[0].moisture, [0.05, 0.2]);
        assert_eq!(report.nonlinear.tridiagonal_solves, 0);
        assert_eq!(next.water_mass_kg(), old.water_mass_kg());
        assert_eq!(next.dry_mass_kg(), old.dry_mass_kg());
        assert_eq!(next.layers()[0].slab().time_s(), 1e50);
    });
}

#[test]
fn g1_linear_two_material_transient_matches_capacity_weighted_partition() {
    with_cx(|cx| {
        let old = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(800.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]),
                    0.2,
                    vec![0.16],
                    cx,
                ),
                layer(
                    material(300.0, vec![(0.0, 0.0), (1.0, 0.4)], [0.02; 2]),
                    0.1,
                    vec![0.2],
                    cx,
                ),
            ],
            256,
            cx,
        )
        .unwrap();
        let dt = 2.0;
        let capacity = [0.16 * 0.2, 0.03 * 0.4];
        let half = [0.016, 0.048];
        let conductance = 0.012;
        // Water moves from lower M=.16 into higher M=.2 because the source
        // activities, .8 and .5, establish the thermodynamic driving direction.
        let mean = (capacity[0] * 0.8 + capacity[1] * 0.5) / (capacity[0] + capacity[1]);
        let difference = 0.3 / (1.0 + dt * conductance * (1.0 / capacity[0] + 1.0 / capacity[1]));
        let expected = [
            mean + capacity[1] / (capacity[0] + capacity[1]) * difference,
            mean - capacity[0] / (capacity[0] + capacity[1]) * difference,
        ];
        let (next, report) = old
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        for (&actual, expected) in next.activity().iter().zip(expected) {
            near(actual, expected, 3e-13);
        }
        let surface = (half[0] * expected[0] + half[1] * expected[1]) / (half[0] + half[1]);
        let received = dt * conductance * difference;
        let face = &report.interfaces[0];
        near(face.relative_activity, surface, 3e-13);
        near(face.moisture[0], 0.2 * surface, 3e-13);
        near(face.moisture[1], 0.4 * surface, 3e-13);
        near(face.water_left_to_right_kg, received, 3e-13);
        near(
            old.layers()[0].slab().water_mass_kg() - next.layers()[0].slab().water_mass_kg(),
            received,
            3e-13,
        );
        near(
            next.layers()[1].slab().water_mass_kg() - old.layers()[1].slab().water_mass_kg(),
            received,
            3e-13,
        );
        near(next.water_mass_kg(), old.water_mass_kg(), 3e-13);
        assert_eq!(report.nonlinear.iterations, 1);
    });
}

#[test]
fn g1_degenerate_unlike_material_interface_wets_the_dry_layer_with_balanced_flux() {
    with_cx(|cx| {
        let old = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(800.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.0, 0.06]),
                    0.2,
                    vec![0.16],
                    cx,
                ),
                layer(
                    material(300.0, vec![(0.0, 0.0), (1.0, 0.4)], [0.0, 0.12]),
                    0.1,
                    vec![0.0],
                    cx,
                ),
            ],
            256,
            cx,
        )
        .unwrap();
        let dt: f64 = 0.1;
        let capacity: [f64; 2] = [0.032, 0.012];
        // Integrated half-cell flux = B_k*(a_cell²-a_surface²).
        let half: [f64; 2] = [0.016, 0.096];
        let b = half[0] * half[1] / (half[0] + half[1]);
        let qa = dt * b * (1.0 / capacity[0].powi(2) - 1.0 / capacity[1].powi(2));
        let qb = -1.0 - 2.0 * dt * b * 0.8 / capacity[0];
        let qc = dt * b * 0.8_f64.powi(2);
        let water = 2.0 * qc / (-qb + (qb * qb - 4.0 * qa * qc).sqrt());
        let expected = [0.8 - water / capacity[0], water / capacity[1]];
        let (next, report) = old
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        for (&actual, expected) in next.activity().iter().zip(expected) {
            near(actual, expected, 3e-12);
        }
        let face = &report.interfaces[0];
        let surface = ((half[0] * expected[0].powi(2) + half[1] * expected[1].powi(2))
            / (half[0] + half[1]))
            .sqrt();
        near(face.relative_activity, surface, 3e-12);
        near(face.water_left_to_right_kg, water, 3e-13);
        near(
            dt * half[0] * (next.activity()[0].powi(2) - face.relative_activity.powi(2)),
            water,
            3e-13,
        );
        near(
            dt * half[1] * (face.relative_activity.powi(2) - next.activity()[1].powi(2)),
            water,
            3e-13,
        );
        assert!(next.layers()[1].slab().moisture()[0] > 0.0);
        near(next.water_mass_kg(), old.water_mass_kg(), 3e-13);
        assert_eq!(old.layers()[1].slab().moisture(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
        let limited = MoistureNonlinearConfig {
            maximum_iterations: 1,
            ..nonlinear()
        };
        assert!(matches!(
            old.advance(dt, &sealed(), config(), limited, cx),
            Err(MoistureError::NonlinearNotConverged(_))
        ));
        assert_eq!(old.layers()[1].slab().moisture(), &[0.0]);
    });
}

#[test]
fn g1_equal_material_layers_reduce_to_the_native_nonlinear_slab() {
    with_cx(|cx| {
        let material = material(500.0, vec![(0.0, 0.0), (0.5, 0.1), (1.0, 0.4)], [0.0, 0.06]);
        let initial = vec![0.2, 0.15, 0.1, 0.05, 0.0, 0.0, 0.0, 0.0];
        let whole = layer(material.clone(), 1.0, initial.clone(), cx);
        let old = MoistureLayeredSlab::new(
            vec![
                layer(material.clone(), 0.5, initial[..4].to_vec(), cx),
                layer(material, 0.5, initial[4..].to_vec(), cx),
            ],
            256,
            cx,
        )
        .unwrap();
        let faces = [
            MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
            MoistureLayeredFace::PrescribedRelativeHumidity(0.2),
        ];
        let native = [
            MoistureFace::PrescribedDryBasis(whole.sorption().moisture(0.8).unwrap()),
            MoistureFace::PrescribedDryBasis(whole.sorption().moisture(0.2).unwrap()),
        ];
        let (next, report) = old.advance(0.2, &faces, config(), nonlinear(), cx).unwrap();
        let (expected, native_report) = whole.slab().advance(0.2, &native, config(), cx).unwrap();
        let actual = next
            .layers()
            .iter()
            .flat_map(|l| l.slab().moisture().iter().copied())
            .collect::<Vec<_>>();
        for (&actual, &expected) in actual.iter().zip(expected.moisture()) {
            near(actual, expected, 5e-12);
        }
        for side in 0..2 {
            near(
                report.boundary_water_kg[side],
                native_report.boundary_water_kg[side],
                3e-13,
            );
        }
        near(next.water_mass_kg(), expected.water_mass_kg(), 3e-13);
    });
}

#[test]
fn g0_impermeable_layer_blocks_internal_and_external_transport_without_losing_water() {
    with_cx(|cx| {
        let old = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(800.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]),
                    0.2,
                    vec![0.04],
                    cx,
                ),
                layer(
                    material(300.0, vec![(0.0, 0.0), (1.0, 0.4)], [0.0; 2]),
                    0.1,
                    vec![0.2],
                    cx,
                ),
            ],
            256,
            cx,
        )
        .unwrap();
        let faces = [
            MoistureLayeredFace::PrescribedRelativeHumidity(0.9),
            MoistureLayeredFace::PrescribedRelativeHumidity(0.0),
        ];
        let (next, report) = old.advance(1.0, &faces, config(), nonlinear(), cx).unwrap();
        assert!(next.layers()[0].slab().moisture()[0] > 0.04);
        assert_eq!(next.layers()[1].slab().moisture(), &[0.2]);
        assert_eq!(report.interfaces[0].water_left_to_right_kg, 0.0);
        assert_eq!(report.boundary_water_kg[1], 0.0);
        near(
            next.water_mass_kg() - old.water_mass_kg(),
            report.boundary_water_kg[0],
            3e-13,
        );
    });
}
