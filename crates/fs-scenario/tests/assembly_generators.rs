//! G0 source preservation/refusals; G1 independent analytic geometry oracles;
//! G3 direct reuse of the .55 mesher and arching lift. All shapes are fixtures.
use fs_exec::{CancelGate, Cx};
use fs_math::det;
use fs_qty::{Dims, QtyAny};
use fs_scenario::ToleranceSource;
use fs_scenario::assembly::generators::*;
use fs_scenario::assembly::{DimensionBound, GeometrySource, ParameterKey, ReferenceState};

const B: GeneratorBudget = GeneratorBudget {
    max_dimensions: 10000,
    max_samples: 10000,
    max_work: 10_000_000,
};
const LABEL: &str = "test-fixture:analytic-geometry";
fn dimension(name: &str, value: f64, error: f64, dims: Dims) -> MeasuredDimension {
    let row = format!("test-fixture:geometry-table-{name}");
    MeasuredDimension {
        value: QtyAny::new(value, dims),
        bound: DimensionBound {
            property: ParameterKey::new(name, dims),
            lower: QtyAny::new(value - error, dims),
            upper: QtyAny::new(value + error, dims),
            source: ToleranceSource::Measured {
                ledger_row: row.clone(),
                method: "analytic test fixture, not a measurement".into(),
            },
        },
        ledger_row: Some(row),
    }
}
fn length(name: &str, v: f64) -> MeasuredDimension {
    dimension(name, v, 0.0, LENGTH)
}
fn stations() -> Vec<(MeasuredDimension, MeasuredDimension)> {
    [(0.0, 0.01), (0.25, 0.018), (0.75, 0.012), (1.0, 0.02)]
        .iter()
        .enumerate()
        .map(|(i, (x, r))| (length(&format!("x{i}"), *x), length(&format!("r{i}"), *r)))
        .collect()
}
fn bore() -> BoreProfile {
    BoreProfile::new(
        LABEL.into(),
        stations(),
        BoreInterpolant::PiecewiseConical,
        B,
    )
    .unwrap()
}
fn loop_points(points: &[[f64; 2]], name: &str) -> Vec<[MeasuredDimension; 2]> {
    points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            [
                length(&format!("{name}{i}x"), p[0]),
                length(&format!("{name}{i}y"), p[1]),
            ]
        })
        .collect()
}
fn outer() -> Vec<[MeasuredDimension; 2]> {
    loop_points(
        &[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]],
        "outer",
    )
}
fn circle_hole() -> Vec<[MeasuredDimension; 2]> {
    loop_points(
        &(0..16)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 16.0;
                match i {
                    0 => [0.25, 0.0],
                    4 => [0.0, 0.25],
                    8 => [-0.25, 0.0],
                    12 => [0.0, -0.25],
                    _ => [0.25 * det::cos(a), 0.25 * det::sin(a)],
                }
            })
            .collect::<Vec<_>>(),
        "hole",
    )
}
fn arch(n: usize, f: impl Fn(f64, f64) -> f64) -> ArchingSurface {
    let x: Vec<_> = (0..=n)
        .map(|i| length(&format!("x{i}"), -1.0 + 2.0 * i as f64 / n as f64))
        .collect();
    let y: Vec<_> = (0..=n)
        .map(|i| length(&format!("y{i}"), -1.0 + 2.0 * i as f64 / n as f64))
        .collect();
    let mut z = Vec::new();
    for (j, yy) in y.iter().enumerate() {
        for (i, xx) in x.iter().enumerate() {
            z.push(length(
                &format!("z{i}_{j}"),
                f(xx.value.value, yy.value.value),
            ));
        }
    }
    ArchingSurface::new(LABEL.into(), x, y, z, None, ReferenceState::Unloaded, B).unwrap()
}
fn assert_relative(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs(),
        "{a:.17e} != {b:.17e}, relative {}",
        (a - b).abs() / b.abs()
    );
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 179,
                kernel_id: 179,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
#[test]
fn ag_001_bore_retains_stations_bit_exact() {
    let profile = bore();
    let source = profile.geometry_source();
    let Generator::BoreProfile(decoded) = Generator::from_source(&source, B).unwrap() else {
        panic!("wrong generator");
    };
    assert_eq!(source, decoded.geometry_source());
    assert_eq!(profile.stations(), decoded.stations());
    for divisions in [1, 2, 4] {
        let sampled = decoded
            .sample(
                SamplingRung {
                    subdivisions_per_span: divisions,
                },
                B,
                || false,
            )
            .unwrap();
        assert_eq!(sampled.stations.len(), 3 * divisions + 1);
        for (i, (x, r)) in decoded.stations().iter().enumerate() {
            assert_eq!(
                sampled.stations[i * divisions][0].value.value.to_bits(),
                x.value.value.to_bits()
            );
            assert_eq!(
                sampled.stations[i * divisions][1].value.value.to_bits(),
                r.value.value.to_bits()
            );
            assert_eq!(
                decoded
                    .radius_at(x.value.value)
                    .unwrap()
                    .value
                    .value
                    .to_bits(),
                r.value.value.to_bits()
            );
        }
    }
    assert!(decoded.radius_at(-1.0).is_err());
}
#[test]
fn ag_002_piecewise_conical_volume_closed_form() {
    let p = bore();
    // Independent frustum formula using inlet/outlet diameters.
    let expected = p
        .stations()
        .windows(2)
        .map(|w| {
            let d0 = 2.0 * w[0].1.value.value;
            let d1 = 2.0 * w[1].1.value.value;
            std::f64::consts::PI
                * (w[1].0.value.value - w[0].0.value.value)
                * (d0 * d0 + d0 * d1 + d1 * d1)
                / 12.0
        })
        .sum::<f64>();
    assert_relative(p.volume().unwrap().value.value, expected, 1e-12);
    for n in [1, 2, 4] {
        let sample = p
            .sample(
                SamplingRung {
                    subdivisions_per_span: n,
                },
                B,
                || false,
            )
            .unwrap();
        let sum = sample
            .segments
            .iter()
            .map(|s| {
                let a = s.value.radius_m;
                let b = s.value.outlet_radius_m;
                std::f64::consts::PI * s.value.length_m * (a * a + a * b + b * b) / 3.0
            })
            .sum();
        assert_relative(sum, expected, 1e-12);
    }
}
#[test]
fn ag_003_smooth_bore_interpolant_converges_three_rungs() {
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let mut table = Vec::new();
        let mut slopes = Vec::new();
        for i in 0..=n {
            let x = i as f64 / n as f64;
            let r = 0.01 * det::exp(0.6 * x);
            table.push((length(&format!("x{i}"), x), length(&format!("r{i}"), r)));
            slopes.push(dimension(&format!("s{i}"), 0.6 * r, 0.0, Dims::NONE));
        }
        let p = BoreProfile::new(
            LABEL.into(),
            table,
            BoreInterpolant::CubicHermite { slopes },
            B,
        )
        .unwrap();
        assert_eq!(p.interpolant().formal_order(), 4);
        let source = p.geometry_source();
        let Generator::BoreProfile(rebuilt) = Generator::from_source(&source, B).unwrap() else {
            panic!("wrong generator");
        };
        assert_eq!(source, rebuilt.geometry_source());
        let error = (0..n)
            .map(|i| {
                let x = (i as f64 + 0.5) / n as f64;
                (p.radius_at(x).unwrap().value.value - 0.01 * det::exp(0.6 * x)).abs()
            })
            .fold(0.0, f64::max);
        errors.push(error);
        let sampled = p
            .sample(
                SamplingRung {
                    subdivisions_per_span: 2,
                },
                B,
                || false,
            )
            .unwrap();
        assert_eq!(sampled.stations.len(), 2 * n + 1);
        assert!(
            p.volume().is_err(),
            "no silent conical reduction of cubic geometry"
        );
    }
    for w in errors.windows(2) {
        let order = det::ln(w[0] / w[1]) / det::ln(2.0);
        eprintln!("exponential horn errors={errors:?} order={order}");
        assert!((order - 4.0).abs() <= 0.2);
    }
}
#[test]
fn ag_004_outline_holes_validate_and_area_exact() {
    let outline = OutlineWithHoles::new(LABEL.into(), outer(), vec![circle_hole()], B).unwrap();
    let expected = 4.0 - 0.5 * 16.0 * 0.25 * 0.25 * det::sin(2.0 * std::f64::consts::PI / 16.0);
    assert_relative(outline.area().unwrap().value.value, expected, 1e-12);
    for (source, vertices) in outline.contours().iter().zip(outline.vertices_m()) {
        for (p, q) in source.iter().zip(vertices) {
            for axis in 0..2 {
                assert_eq!(p[axis].value.value.to_bits(), q[axis].to_bits());
            }
        }
    }
    let source = outline.geometry_source();
    let Generator::OutlineWithHoles(rebuilt) = Generator::from_source(&source, B).unwrap() else {
        panic!("wrong type");
    };
    assert_eq!(source, rebuilt.geometry_source());
    use fs_rep_nurbs::{Classification, Rat};
    assert_eq!(
        outline
            .trimmed_patch()
            .classify([Rat::int(0), Rat::int(0)])
            .unwrap(),
        Classification::Outside
    );
    assert_eq!(
        outline
            .trimmed_patch()
            .classify([Rat::new(1, 2), Rat::int(0)])
            .unwrap(),
        Classification::Inside
    );
}
#[test]
fn ag_005_hole_outside_outline_refuses() {
    let hole = loop_points(&[[0.5, 0.0], [1.5, 0.0], [1.5, 0.5], [0.5, 0.5]], "outside");
    assert!(matches!(
        OutlineWithHoles::new(LABEL.into(), outer(), vec![hole], B),
        Err(GeneratorError::HoleOutsideOutline { hole: 0 })
    ));
    let hole = circle_hole();
    assert!(matches!(
        OutlineWithHoles::new(LABEL.into(), outer(), vec![hole.clone(), hole], B),
        Err(GeneratorError::OverlappingHoles)
    ));
}
#[test]
fn ag_006_self_intersecting_loop_refuses() {
    let points = loop_points(
        &[[-1.0, -1.0], [1.0, 1.0], [-1.0, 1.0], [1.0, -1.0]],
        "crossed",
    );
    assert!(matches!(
        OutlineWithHoles::new(LABEL.into(), points, vec![], B),
        Err(GeneratorError::SelfIntersection { loop_index: 0 })
    ));
}
#[test]
fn ag_007_arching_retains_contour_heights() {
    let a = arch(4, |x, y| 0.01 * (1.0 - x * x) * (1.0 - y * y));
    let (x, y) = a.axes();
    for (j, yy) in y.iter().enumerate() {
        for (i, xx) in x.iter().enumerate() {
            let p = [xx.value.value, yy.value.value];
            let expected = a.contour_heights()[j * x.len() + i].value.value;
            assert_eq!(
                a.height_at(p).unwrap().value.value.to_bits(),
                expected.to_bits()
            );
            assert_eq!(
                a.surface().eval(p[0], p[1]).unwrap()[2].to_bits(),
                expected.to_bits()
            );
        }
    }
    assert!(a.height_at([2.0, 0.0]).is_err());
    assert!(a.thickness_at([0.0, 0.0]).is_err());
    let source = a.geometry_source();
    let Generator::ArchingSurface(rebuilt) = Generator::from_source(&source, B).unwrap() else {
        panic!("wrong generator");
    };
    assert_eq!(source, rebuilt.geometry_source());
    assert_eq!(rebuilt.reference_state(), ReferenceState::Unloaded);
}
#[test]
fn ag_008_uncertainty_propagates_to_derived_dims() {
    let mut table = stations();
    for (x, r) in &mut table {
        *x = dimension(&x.bound.property.name, x.value.value, 1e-5, LENGTH);
        *r = dimension(&r.bound.property.name, r.value.value, 1e-4, LENGTH);
    }
    let p = BoreProfile::new(
        LABEL.into(),
        table.clone(),
        BoreInterpolant::PiecewiseConical,
        B,
    )
    .unwrap();
    let volume = p.volume().unwrap();
    assert_eq!(volume.inputs.len(), 8);
    // Exhaust independent endpoint perturbations; no probabilistic independence
    // assumption is used in the interval propagation.
    for mask in 0..256 {
        let values: Vec<_> = table
            .iter()
            .enumerate()
            .map(|(i, (x, r))| {
                let choose = |d: &MeasuredDimension, bit| {
                    if (mask >> bit) & 1 == 0 {
                        d.bound.lower.value
                    } else {
                        d.bound.upper.value
                    }
                };
                (choose(x, 2 * i), choose(r, 2 * i + 1))
            })
            .collect();
        let v = values
            .windows(2)
            .map(|w| {
                std::f64::consts::PI
                    * (w[1].0 - w[0].0)
                    * (w[0].1 * w[0].1 + w[0].1 * w[1].1 + w[1].1 * w[1].1)
                    / 3.0
            })
            .sum::<f64>();
        assert!(v >= volume.lower.value && v <= volume.upper.value);
    }
    for sample in p
        .sample(
            SamplingRung {
                subdivisions_per_span: 2,
            },
            B,
            || false,
        )
        .unwrap()
        .segments
    {
        for dim in [sample.length, sample.inlet_radius, sample.outlet_radius] {
            assert_eq!(dim.inputs, volume.inputs);
            assert!(dim.lower.value <= dim.value.value && dim.value.value <= dim.upper.value);
        }
    }
    let mut polygon = outer();
    for p in &mut polygon {
        for d in p {
            *d = dimension(&d.bound.property.name, d.value.value, 1e-3, LENGTH);
        }
    }
    let area = OutlineWithHoles::new(LABEL.into(), polygon, vec![], B)
        .unwrap()
        .area()
        .unwrap();
    assert!(area.lower.value < 4.0 && area.upper.value > 4.0);
    assert_eq!(area.inputs.len(), 8);
    let a = arch(2, |x, y| 0.02 + x * x + y * y);
    let (x, y) = a.axes();
    let heights = a
        .contour_heights()
        .iter()
        .map(|h| dimension(&h.bound.property.name, h.value.value, 0.001, LENGTH))
        .collect::<Vec<_>>();
    let thickness = heights
        .iter()
        .enumerate()
        .map(|(i, _)| dimension(&format!("t{i}"), 0.003, 0.0002, LENGTH))
        .collect();
    let a = ArchingSurface::new(
        LABEL.into(),
        x.to_vec(),
        y.to_vec(),
        heights,
        Some(thickness),
        ReferenceState::Assembled,
        B,
    )
    .unwrap();
    for d in [
        a.height_at([0.25, 0.25]).unwrap(),
        a.thickness_at([0.25, 0.25]).unwrap(),
    ] {
        assert!(d.lower.value < d.value.value && d.value.value < d.upper.value);
        assert_eq!(d.inputs.len(), 15);
    }
    let source = a.geometry_source();
    let Generator::ArchingSurface(rebuilt) = Generator::from_source(&source, B).unwrap() else {
        panic!();
    };
    assert_eq!(source, rebuilt.geometry_source());
}
#[test]
fn ag_009_unsourced_dimension_requires_fixture_label() {
    // These citation tokens exercise structural admission, not source authority.
    let cite = |d: &mut MeasuredDimension| {
        let row = "M8-W1:structural-test-row".to_owned();
        d.ledger_row = Some(row.clone());
        d.bound.source = ToleranceSource::Measured {
            ledger_row: row,
            method: "structural admission test only".into(),
        };
    };
    let mut sourced = stations();
    for (x, r) in &mut sourced {
        cite(x);
        cite(r);
    }
    let measured = BoreProfile::new(
        "declared-table:irregular-bore".into(),
        sourced,
        BoreInterpolant::PiecewiseConical,
        B,
    )
    .unwrap();
    assert!(
        Generator::from_source(&measured.geometry_source(), B)
            .unwrap()
            .require_specimen()
            .is_ok()
    );
    let mut rectangle = outer();
    for p in &mut rectangle {
        for d in p {
            cite(d);
        }
    }
    assert!(matches!(
        OutlineWithHoles::new("declared-table:rectangle".into(), rectangle, vec![], B),
        Err(GeneratorError::FixtureRequired)
    ));
    let mut axis_x = vec![length("x0", 0.0), length("x1", 1.0)];
    let mut axis_y = vec![length("y0", 0.0), length("y1", 1.0)];
    let mut heights = vec![
        length("z0", 0.0),
        length("z1", 1.0),
        length("z2", 2.0),
        length("z3", 3.0),
    ];
    for d in axis_x.iter_mut().chain(&mut axis_y).chain(&mut heights) {
        cite(d);
    }
    assert!(matches!(
        ArchingSurface::new(
            "declared-table:tilted-plane".into(),
            axis_x,
            axis_y,
            heights,
            None,
            ReferenceState::Unloaded,
            B
        ),
        Err(GeneratorError::FixtureRequired)
    ));
    let mut table = stations();
    for (x, r) in &mut table {
        for d in [x, r] {
            d.ledger_row = None;
            d.bound.source = ToleranceSource::Assumed {
                rationale: "analytic test fixture".into(),
            };
        }
    }
    assert!(matches!(
        BoreProfile::new(
            "specimen:unsourced".into(),
            table.clone(),
            BoreInterpolant::PiecewiseConical,
            B
        ),
        Err(GeneratorError::MissingSource { .. })
    ));
    let p = BoreProfile::new(LABEL.into(), table, BoreInterpolant::PiecewiseConical, B).unwrap();
    assert!(
        Generator::from_source(&p.geometry_source(), B)
            .unwrap()
            .require_specimen()
            .is_err()
    );
    // Even a cited row cannot turn an ideal cone into specimen geometry.
    let mut cone = vec![
        (length("x0", 0.0), length("r0", 0.01)),
        (length("x1", 1.0), length("r1", 0.02)),
    ];
    for (x, r) in &mut cone {
        for d in [x, r] {
            let row = "M8-W1:declared-cone".to_owned();
            d.ledger_row = Some(row.clone());
            d.bound.source = ToleranceSource::Measured {
                ledger_row: row,
                method: "dimension table".into(),
            };
        }
    }
    assert!(matches!(
        BoreProfile::new(
            "specimen:ideal-cone".into(),
            cone,
            BoreInterpolant::PiecewiseConical,
            B
        ),
        Err(GeneratorError::FixtureRequired)
    ));
    let mut source = p.geometry_source();
    if let GeometrySource::Parametric { params, .. } = &mut source {
        params[4].1 = QtyAny::dimensionless(9.0);
    }
    assert!(Generator::from_source(&source, B).is_err());
}
#[test]
fn ag_010_station_budget_cap_and_cap_plus_one() {
    let table = stations();
    let exact = GeneratorBudget {
        max_dimensions: 8,
        ..B
    };
    assert!(
        BoreProfile::new(
            LABEL.into(),
            table.clone(),
            BoreInterpolant::PiecewiseConical,
            exact
        )
        .is_ok()
    );
    assert!(matches!(
        BoreProfile::new(
            LABEL.into(),
            table,
            BoreInterpolant::PiecewiseConical,
            GeneratorBudget {
                max_dimensions: 7,
                ..B
            }
        ),
        Err(GeneratorError::BudgetExceeded {
            resource: "dimensions",
            cap: 7
        })
    ));
    let p = bore();
    let sample = GeneratorBudget {
        max_samples: 7,
        ..B
    };
    assert!(
        p.sample(
            SamplingRung {
                subdivisions_per_span: 2
            },
            sample,
            || false
        )
        .is_ok()
    );
    assert!(matches!(
        p.sample(
            SamplingRung {
                subdivisions_per_span: 2
            },
            GeneratorBudget {
                max_samples: 6,
                ..B
            },
            || false
        ),
        Err(GeneratorError::BudgetExceeded {
            resource: "samples",
            cap: 6
        })
    ));
    assert_eq!(
        p.sample(
            SamplingRung {
                subdivisions_per_span: 1
            },
            B,
            || true
        ),
        Err(GeneratorError::Cancelled)
    );
    assert!(
        p.sample(
            SamplingRung {
                subdivisions_per_span: 0
            },
            B,
            || false
        )
        .is_err()
    );
    assert!(
        OutlineWithHoles::new(
            LABEL.into(),
            outer(),
            vec![],
            GeneratorBudget { max_work: 15, ..B }
        )
        .is_err()
    );
    assert!(
        OutlineWithHoles::new(
            LABEL.into(),
            outer(),
            vec![],
            GeneratorBudget { max_work: 16, ..B }
        )
        .is_ok()
    );
    let a = arch(2, |x, y| 1.0 + x * x + y * y);
    let (x, y) = a.axes();
    assert!(
        ArchingSurface::new(
            LABEL.into(),
            x.to_vec(),
            y.to_vec(),
            a.contour_heights().to_vec(),
            None,
            ReferenceState::Unloaded,
            GeneratorBudget {
                max_dimensions: 15,
                ..B
            }
        )
        .is_ok()
    );
    assert!(matches!(
        ArchingSurface::new(
            LABEL.into(),
            x.to_vec(),
            y.to_vec(),
            a.contour_heights().to_vec(),
            None,
            ReferenceState::Unloaded,
            GeneratorBudget {
                max_dimensions: 14,
                ..B
            }
        ),
        Err(GeneratorError::BudgetExceeded {
            resource: "dimensions",
            cap: 14
        })
    ));
}
#[test]
fn ag_011_arching_second_order_and_existing_mesh_lift() {
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let a = arch(n, |x, y| 0.01 * (x * x + y * y));
        let p = [-1.0 + 1.0 / n as f64, -1.0 + 1.0 / n as f64];
        errors
            .push((a.height_at(p).unwrap().value.value - 0.01 * (p[0] * p[0] + p[1] * p[1])).abs());
    }
    for w in errors.windows(2) {
        let order = det::ln(w[0] / w[1]) / det::ln(2.0);
        assert!((order - ArchingSurface::FORMAL_ORDER as f64).abs() <= 0.2);
    }
    let a = arch(4, |x, y| 0.01 * (1.0 - x * x) * (1.0 - y * y));
    let outline = OutlineWithHoles::new(LABEL.into(), outer(), vec![circle_hole()], B).unwrap();
    let loops = outline.vertices_m();
    use fs_mesh::planar::{Loop, PlanarBudget, PlanarOptions, PlanarOutline, mesh_outline};
    let mesh = with_cx(|cx| {
        mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&loops[0]),
                holes: vec![Loop::Polyline(&loops[1])],
            },
            &PlanarOptions {
                maximum_edge_m: 10.0,
                minimum_angle_deg: 0.0,
                size_field: None,
                budget: PlanarBudget {
                    max_vertices: 10000,
                    max_triangles: 20000,
                    max_refinement_steps: 10000,
                },
            },
            cx,
        )
        .unwrap()
    });
    assert_relative(
        mesh.receipt.mesh_area_m2,
        outline.area().unwrap().value.value,
        1e-12,
    );
    assert_eq!(mesh.receipt.boundary_loops.len(), 2);
    let plate = fs_plate::PlateMesh::from_unstructured(
        mesh.vertices.iter().map(|p| (p[0], p[1])).collect(),
        mesh.triangles.clone(),
    )
    .unwrap();
    use fs_plate::shell::arching::{
        ArchingParameters, ArchingSurface as PlateArch, ReferenceLoadState,
    };
    let arch = PlateArch::new(
        ArchingParameters::NurbsHeightMap(a.surface().clone()),
        ReferenceLoadState::Unloaded,
    )
    .unwrap();
    // 25 controls x 4 homogeneous coordinates x 5 evaluations: a node on a
    // seam corner of this degree-one chart needs four one-sided partials and
    // its exact point. The hole's axis points lie on the x = 0 and y = 0 seams.
    let budget = fs_plate::shell::profile::ProfileBudget {
        max_nodes: plate.nodes.len(),
        max_triangles: plate.tris.len(),
        max_feature_evaluations: 500 * plate.nodes.len(),
    };
    let lifted = arch.lift(&plate, budget, || false).unwrap();
    assert_eq!(lifted.mesh.tris, plate.tris);
    for ((x, y), p) in plate.nodes.iter().zip(&lifted.mesh.nodes) {
        assert_eq!(p[0].to_bits(), x.to_bits());
        assert_eq!(p[1].to_bits(), y.to_bits());
        assert!((p[2] - a.height_at([*x, *y]).unwrap().value.value).abs() < 1e-14);
    }
    assert!(arch.lift(&plate, budget, || true).is_err());
    assert!(
        arch.lift(
            &plate,
            fs_plate::shell::profile::ProfileBudget {
                max_feature_evaluations: budget.max_feature_evaluations - 1,
                ..budget
            },
            || false
        )
        .is_err()
    );
    // A chart whose x coordinate is 2u is not a physical-coordinate graph.
    let knots = fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let chart = fs_rep_nurbs::NurbsSurface::new(
        knots.clone(),
        knots,
        &[
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
        ],
        &[vec![1.0; 2], vec![1.0; 2]],
    )
    .unwrap();
    let invalid = PlateArch::new(
        ArchingParameters::NurbsHeightMap(std::sync::Arc::new(chart)),
        ReferenceLoadState::Unloaded,
    )
    .unwrap();
    assert!(invalid.sample([0.5, 0.5]).is_err());
}

#[test]
fn ag_012_generator_tables_survive_assembly_wire_and_revalidate() {
    use fs_scenario::assembly::{
        AssemblyDescription, AssemblyInput, DEFAULT_ASSEMBLY_BUDGET, PartDecl, PartKind,
    };
    use fs_scenario::{EntityCatalog, EntityDeclaration, Environment, Scenario};
    let scenario = Scenario::new(LABEL, 179, Environment::earth_lab());
    let mut catalog = EntityCatalog::new();
    let root = catalog.declare(EntityDeclaration::assembly(LABEL)).unwrap();
    let mut table = stations();
    table[0].0 = length("signed_zero", -0.0);
    let b = BoreProfile::new(LABEL.into(), table, BoreInterpolant::PiecewiseConical, B).unwrap();
    let outline = OutlineWithHoles::new(LABEL.into(), outer(), vec![circle_hole()], B).unwrap();
    let a = arch(2, |x, y| 0.02 + x * x + y * y);
    let sources = [
        b.geometry_source(),
        outline.geometry_source(),
        a.geometry_source(),
    ];
    let mut parts = Vec::new();
    for (i, source) in sources.iter().enumerate() {
        let entity = catalog
            .declare(EntityDeclaration::part(
                root,
                &format!("test-fixture:part{i}"),
            ))
            .unwrap();
        parts.push(PartDecl {
            entity,
            kind: PartKind::Shell,
            geometry: source.clone(),
            materials: vec![],
            uncertainty: vec![],
        });
    }
    let input = AssemblyInput {
        root,
        catalog,
        parts,
        interfaces: vec![],
        configurations: vec![],
        observers: vec![],
        reference_state: Some(ReferenceState::Unloaded),
    };
    let description =
        AssemblyDescription::new(&scenario, input.clone(), DEFAULT_ASSEMBLY_BUDGET).unwrap();
    let bytes = description.to_bytes(DEFAULT_ASSEMBLY_BUDGET).unwrap();
    let decoded =
        AssemblyDescription::from_bytes(&bytes, &scenario, DEFAULT_ASSEMBLY_BUDGET).unwrap();
    assert_eq!(bytes, decoded.to_bytes(DEFAULT_ASSEMBLY_BUDGET).unwrap());
    for (part, original) in decoded.input().parts.iter().zip(sources) {
        assert_eq!(part.geometry, original);
    }
    let Generator::BoreProfile(b) =
        Generator::from_source(&decoded.input().parts[0].geometry, B).unwrap()
    else {
        panic!();
    };
    assert_eq!(b.stations()[0].0.value.value.to_bits(), (-0.0f64).to_bits());
    let mut broken = input;
    if let GeometrySource::Parametric { params, .. } = &mut broken.parts[0].geometry {
        params.pop();
    }
    assert!(AssemblyDescription::new(&scenario, broken, DEFAULT_ASSEMBLY_BUDGET).is_err());
}
