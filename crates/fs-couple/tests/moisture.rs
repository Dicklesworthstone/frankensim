//! G1 transport limits and actual native mechanical response to local moisture.
//! Diffusivity/elasticity fixtures are AUTHORED, not measured wood validation.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::moisture::plate::*;
use fs_couple::moisture::*;
use fs_couple::thin_plate::PlateMaterialModel;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, SliceOptions, modes};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};
use fs_solver::CheckedCgConfig;

#[path = "moisture/solid.rs"]
mod solid;
#[path = "moisture/atmosphere.rs"]
mod atmosphere;
#[path = "moisture/surface.rs"]
mod surface;
#[path = "moisture/implicit.rs"]
mod implicit;

fn qty(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(m: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            m,
        )
        .unwrap()
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn card(diffusivity: [f64; 2], elastic_upper: f64) -> MaterialCard {
    card_with_swelling(diffusivity, elastic_upper, [0.01, 0.04, 0.02])
}
fn card_with_swelling(diffusivity: [f64; 2], elastic_upper: f64, beta: [f64; 3]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims: Dims, a: f64, b: f64, lower: f64, upper: f64| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "moisture_content".into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, a), (0.3, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        qty(QuantityKind::AbsoluteTemperature),
                        293.15,
                        293.15,
                    )
                    .with_quantity(
                        "moisture_content",
                        qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
                        lower,
                        upper,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored Fick/elastic verification specimen".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        "moisture_diffusivity",
        MOISTURE_DIFFUSIVITY_DIMS,
        diffusivity[0],
        diffusivity[1],
        0.0,
        0.3,
    );
    // Density is known ONLY at reference moisture. Transport must not re-query
    // this curve at the wet endpoint or manufacture another dry inventory.
    add("density", Density::DIMS, 550.0, 550.0, 0.1, 0.1);
    add(
        "young_modulus",
        Pressure::DIMS,
        1.0e9,
        0.7e9,
        0.0,
        elastic_upper,
    );
    add("poisson_ratio", Dims::NONE, 0.25, 0.25, 0.0, elastic_upper);
    add(
        "young_modulus_1",
        Pressure::DIMS,
        2.0e9,
        1.4e9,
        0.0,
        elastic_upper,
    );
    add(
        "young_modulus_2",
        Pressure::DIMS,
        0.5e9,
        0.35e9,
        0.0,
        elastic_upper,
    );
    add("poisson_ratio_12", Dims::NONE, 0.2, 0.2, 0.0, elastic_upper);
    add(
        "shear_modulus_12",
        Pressure::DIMS,
        0.3e9,
        0.21e9,
        0.0,
        elastic_upper,
    );
    add("young_modulus_3", Pressure::DIMS, 0.4e9, 0.28e9, 0.0, elastic_upper);
    add("poisson_ratio_13", Dims::NONE, 0.1, 0.1, 0.0, elastic_upper);
    add("poisson_ratio_23", Dims::NONE, 0.2, 0.2, 0.0, elastic_upper);
    add("shear_modulus_23", Pressure::DIMS, 0.2e9, 0.14e9, 0.0, elastic_upper);
    add("shear_modulus_31", Pressure::DIMS, 0.25e9, 0.175e9, 0.0, elastic_upper);
    for (key,value) in fs_material::state_point::ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES.into_iter().zip(beta) {
        add(key, Dims::NONE, value, value, 0.0, 0.3);
    }
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::with_quantity(fs_material::state_point::HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY,
            qty(QuantityKind::Moisture(MoistureBasis::DryMass))),
        value: PropertyValue::Scalar { value: 0.1, dims: Dims::NONE },
        validity: ValidityDomain::unconstrained().with_quantity("moisture_content",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.3)
            .with_quantity("temperature", qty(QuantityKind::AbsoluteTemperature),293.15,293.15),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "authored stress-free moisture datum".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored sorbing solid".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn geometry(length: f64) -> MoistureSlabGeometry {
    MoistureSlabGeometry {
        length_m: length,
        cross_section_m2: 0.001,
        maximum_moisture: 0.3,
        maximum_cells: 256,
    }
}
fn slab(material: MaterialCard, length: f64, contents: Vec<f64>, cx: &Cx<'_>) -> MoistureSlab {
    MoistureSlab::new(
        material,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        geometry(length),
        contents,
        cx,
    )
    .unwrap()
}
fn config() -> MoistureStepConfig {
    MoistureStepConfig {
        linear: CheckedCgConfig {
            tolerance: 1.0e-13,
            max_iterations: 2048,
            max_corrections: 2,
        },
        water_balance_tolerance_kg: 1.0e-11,
        maximum_content_increment: 0.3,
    }
}
fn faces(m: f64) -> [MoistureFace; 2] {
    [
        MoistureFace::PrescribedDryBasis(m),
        MoistureFace::PrescribedDryBasis(m),
    ]
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b} (tol {tol})");
}

#[test]
fn g1_one_cell_has_half_cell_boundary_resistance_and_conserves_water() {
    with_cx(|cx| {
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let (next, report) = old.advance(2.0, &faces(0.2), config(), cx).unwrap();
        // m'=(m+2*dt*D/L²*(left+right))/(1+4*dt*D/L²).
        let expected = (0.1 + 0.04 * 0.4) / 1.08;
        near(next.moisture()[0], expected, 1e-14);
        near(next.dry_mass_kg(), 0.5, 1e-14);
        near(
            next.wet_mass_kg() - old.wet_mass_kg(),
            report.boundary_water_kg.iter().sum(),
            1e-14,
        );
        near(report.water_balance_residual_kg, 0.0, 1e-14);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.moisture(), &[0.1]);
        assert!(next.moisture()[0] < 0.2);
    });
}

#[test]
fn g0_sealed_transport_smooths_the_field_without_changing_water() {
    with_cx(|cx| {
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.2, 0.1, 0.05, 0.1], cx);
        let (next, report) = old
            .advance(
                1.0,
                &[MoistureFace::Sealed, MoistureFace::Sealed],
                config(),
                cx,
            )
            .unwrap();
        near(old.water_mass_kg(), next.water_mass_kg(), 1e-13);
        assert!(next.moisture()[0] < old.moisture()[0]);
        assert!(next.moisture()[2] > old.moisture()[2]);
        assert_eq!(report.boundary_water_kg, [0.0; 2]);
        let replay = old
            .advance(
                1.0,
                &[MoistureFace::Sealed, MoistureFace::Sealed],
                config(),
                cx,
            )
            .unwrap()
            .0;
        assert_eq!(next.moisture(), replay.moisture());
    });
}

#[test]
fn g1_slab_diffusion_approaches_the_independent_fourier_series() {
    with_cx(|cx| {
        let run = |dt: f64, steps: usize| {
            let mut state = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 64], cx);
            for _ in 0..steps {
                state = state.advance(dt, &faces(0.2), config(), cx).unwrap().0;
            }
            state.moisture().iter().sum::<f64>() / 64.0
        };
        // Exact mean of a unit slab with both surfaces stepped from .1 to .2.
        let exact = 0.2
            - 0.1
                * (0..100)
                    .map(|i| {
                        let k = (2 * i + 1) as f64;
                        8.0 / (std::f64::consts::PI.powi(2) * k * k)
                            * (-0.01 * k * k * std::f64::consts::PI.powi(2) * 2.0).exp()
                    })
                    .sum::<f64>();
        let coarse = run(0.2, 10);
        let fine = run(0.05, 40);
        assert!((fine - exact).abs() < (coarse - exact).abs());
        near(fine, exact, 2e-4);
    });
}

#[test]
fn g3_length_squared_time_scaling_and_data_only_diffusivity_changes() {
    with_cx(|cx| {
        let small = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 8], cx);
        let large = slab(card([0.01; 2], 0.3), 2.0, vec![0.1; 8], cx);
        let fast = slab(card([0.02; 2], 0.3), 1.0, vec![0.1; 8], cx);
        let a = small.advance(1.0, &faces(0.2), config(), cx).unwrap().0;
        let b = large.advance(4.0, &faces(0.2), config(), cx).unwrap().0;
        let c = fast.advance(0.5, &faces(0.2), config(), cx).unwrap().0;
        assert_eq!(a.moisture(), b.moisture());
        assert_eq!(a.moisture(), c.moisture());
        assert!(a.moisture()[0] > a.moisture()[3]);
        near(b.dry_mass_kg(), 2.0 * a.dry_mass_kg(), 1e-14);
    });
}

#[test]
fn changing_diffusivity_is_refreshed_locally_and_changes_the_next_interval() {
    with_cx(|cx| {
        let old = slab(card([0.01, 0.04], 0.3), 1.0, vec![0.1; 4], cx);
        let (next, report) = old.advance(1.0, &faces(0.2), config(), cx).unwrap();
        let property = |s: &fs_material::state_point::ResolvedMaterialStatePoint| {
            s.property("moisture_diffusivity").unwrap().value_si()
        };
        for (state, &m) in report.transport_end.iter().zip(next.moisture()) {
            near(property(state), 0.01 + 0.1 * m, 1e-14);
        }
        let (_, subsequent) = next.advance(1.0, &faces(0.2), config(), cx).unwrap();
        assert_eq!(report.transport_end, subsequent.transport_start);
        assert!(property(&report.transport_end[0]) > property(&report.transport_start[0]));
    });
}

#[test]
fn real_fpl_sorption_changes_the_surface_without_assigning_bulk_equilibrium() {
    with_cx(|cx| {
        let pack = NormalizedPack::from_bytes(include_bytes!(
            "../examples/grand_piano/packs/wood-emc-hailwood-horrobin-fpl-gtr282.fsmatpk"
        ))
        .unwrap();
        let sorption = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "FPL sorption".into(),
                phase: "solid".into(),
                process: "handbook".into(),
                revision: 0,
            },
            pack.claims().clone(),
            vec![],
        )
        .unwrap();
        let ambient = MoistureFace::WoodEquilibrium {
            card: sorption,
            temperature_k: 293.15,
            relative_humidity: 0.8,
        };
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 16], cx);
        let (next, report) = old
            .advance(0.1, &[ambient.clone(), ambient], config(), cx)
            .unwrap();
        let equilibrium = report.surface_moisture[0].unwrap();
        assert!(equilibrium > 0.15 && equilibrium < 0.18);
        assert!(next.moisture()[0] > next.moisture()[7]);
        assert!(next.moisture().iter().all(|m| *m < equilibrium));
    });
}

fn plate(slab: &MoistureSlab, model: PlateMaterialModel, cx: &Cx<'_>) -> MoisturePlate {
    let n = slab.moisture().len();
    MoisturePlate::new(
        slab,
        0.001,
        4,
        model,
        0.0,
        MaterialPropertySelection::SingleClaimOnly,
        PlateMesh::rectangle_boundary(n, 4),
        AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported,
        },
        MoisturePlateBudget {
            maximum_nodes: 1024,
            maximum_triangles: 2048,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn native_plate_modes_follow_local_source_stiffness_and_transported_wet_mass() {
    with_cx(|cx| {
        let material = card([0.01; 2], 0.3);
        let initial = slab(material.clone(), 1.0, vec![0.1; 4], cx);
        let wet = slab(material, 1.0, vec![0.2; 4], cx);
        let binding = plate(&initial, PlateMaterialModel::Isotropic, cx);
        let before = binding.resolve(&initial, cx).unwrap();
        let after = binding.resolve(&wet, cx).unwrap();
        let base = modes(before.model(), (20.0, 90.0), &SliceOptions::default()).unwrap();
        let hydrated = modes(after.model(), (15.0, 90.0), &SliceOptions::default()).unwrap();
        assert_eq!(base.modes.len(), 1);
        assert_eq!(hydrated.modes.len(), 1);
        let predicted_ratio = (0.8 / 0.9) / (1.2 / 1.1);
        near(
            hydrated.modes[0].lambda / base.modes[0].lambda,
            predicted_ratio,
            1e-7,
        );
        near(before.wet_mass_kg(), 0.55, 1e-13);
        near(after.wet_mass_kg(), 0.6, 1e-13);
        let step = binding
            .advance(&initial, 1.0, &faces(0.2), config(), cx)
            .unwrap();
        assert!(
            step.plate.elastic_states()[0]
                .property("young_modulus")
                .unwrap()
                .value_si()
                < step.plate.elastic_states()[1]
                    .property("young_modulus")
                    .unwrap()
                    .value_si()
        );
        assert_eq!(step.plate.time_s(), step.moisture.time_s());
        near(
            step.plate.wet_mass_kg() - before.wet_mass_kg(),
            step.transport.boundary_water_kg.iter().sum(),
            1e-12,
        );
        // Orthotropy uses the same owner and full local directional constants.
        let orthotropic = plate(&initial, PlateMaterialModel::Orthotropic21, cx);
        let ortho = orthotropic
            .advance(&initial, 1.0, &faces(0.2), config(), cx)
            .unwrap();
        assert!(ortho.plate.sections()[0].d[4] > ortho.plate.sections()[0].d[0]);
    });
}

#[test]
fn unsupported_mechanical_endpoint_and_budget_leave_both_accepted_domains_unchanged() {
    with_cx(|cx| {
        let old = slab(card([0.01; 2], 0.105), 1.0, vec![0.1; 4], cx);
        let binding = plate(&old, PlateMaterialModel::Isotropic, cx);
        assert!(
            binding
                .advance(&old, 1.0, &faces(0.2), config(), cx)
                .is_err()
        );
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.moisture(), &[0.1; 4]);
        let mut exhausted = config();
        exhausted.linear.max_iterations = 0;
        assert!(matches!(
            old.advance(1.0, &faces(0.2), exhausted, cx),
            Err(MoistureError::NotConverged(_))
        ));
        let wet_basis = QueryPoint::new()
            .with_quantity(
                "temperature",
                qty(QuantityKind::AbsoluteTemperature),
                293.15,
            )
            .unwrap()
            .with_quantity(
                "moisture_content",
                qty(QuantityKind::Moisture(MoistureBasis::WetMass)),
                0.1,
            )
            .unwrap();
        assert!(matches!(
            MoistureSlab::new(
                card([0.01; 2], 0.3),
                wet_basis,
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                geometry(1.0),
                vec![0.1],
                cx
            ),
            Err(MoistureError::Invalid(_))
        ));
    });
}

#[test]
fn cancelled_transport_cannot_publish_a_moisture_or_mechanical_endpoint() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 4], &cx);
        let binding = plate(&old, PlateMaterialModel::Isotropic, &cx);
        gate.request();
        assert!(matches!(
            binding.advance(&old, 1.0, &faces(0.2), config(), &cx),
            Err(MoistureError::Cancelled)
        ));
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.moisture(), &[0.1; 4]);
    });
}

