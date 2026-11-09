//! Native source-temperature adjoints through filtered physical density design.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticAdjointSettings3, ThermoelasticDesignObjective3, ThermoelasticMotion3,
    ThermoelasticObjective3, ThermoelasticResponse3,
};
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_ivl::Interval;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{AffineIsotropicThermoelastic, IsothermalThermoelasticInput};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::LinearOp;
use fs_topopt::sdf3::thermal::{
    CutThermoelasticError3, CutThermoelasticEvaluation3, CutThermoelasticObjective3,
    ThermoelasticSdf3Elasticity,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3, Sdf3Elasticity};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let v = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(v, v)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn params() -> SimpParams {
    SimpParams {
        penal: 3.0,
        beta: 2.0,
        eta: 0.5,
        e_min: 0.03,
    }
}
fn settings() -> ThermoelasticAdjointSettings3 {
    ThermoelasticAdjointSettings3 {
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 5000,
        poll_iterations: 8,
        maximum_kinematic_strain: 0.1,
        maximum_points: 100_000,
    }
}
fn law() -> AffineIsotropicThermoelastic {
    let quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0, -1.8),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.3, 0.0001),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            0.0001,
            0.000001,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, value - 20.0 * slope),
                        (300.0, value),
                        (360.0, value + 60.0 * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", quantity, 280.0, 360.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored thermal density derivative fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermoelastic source".into(),
            phase: "solid".into(),
            process: "numerical derivative verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("T", quantity, 300.0)
        .unwrap();
    let elastic = resolve_material_state_point(
        &card,
        &query,
        &[
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY,
                Dims([0; 6]),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::Finite,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    AffineIsotropicThermoelastic::new(&card, &elastic, &expansion, [280.0, 360.0], 0.05).unwrap()
}
fn cartesian(delta: f64, formulation: ElasticityFormulation3) -> CutElasticity3 {
    let source = law();
    let field = |cell: HexCell| {
        let p = cell.lo();
        Ok(source
            .at_temperature(316.0 + 3.17 * p[0] + 2.31 * p[1] + 1.07 * p[2] + delta)
            .unwrap()
            .isothermal_elasticity(IsothermalThermoelasticInput {
                reference_density_kg_m3: 3.0 + p[0],
                mechanical_material_identity: ContentHash([73; 32]),
            })
            .unwrap())
    };
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Top,
        ElasticMaterialField3::ThermoelasticByCell(&field),
        &|p| p[0] == 0.0,
        &|_, n| n[2] > 0.5,
        ElasticityOptions3 {
            formulation,
            ghost_gamma: 2.7,
            ..Default::default()
        },
        Default::default(),
        Default::default(),
        &mut q,
    )
    .unwrap()
}
fn adaptive() -> AdaptiveElasticity3 {
    let source = law();
    let field = |cell: HexCell| {
        let p = cell.lo();
        Ok(source
            .at_temperature(316.0 + 3.17 * p[0] + 2.31 * p[1] + 1.07 * p[2])
            .unwrap()
            .isothermal_elasticity(IsothermalThermoelasticInput {
                reference_density_kg_m3: 3.0 + p[0],
                mechanical_material_identity: ContentHash([73; 32]),
            })
            .unwrap())
    };
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|cell| cell.index() == [0, 0, 1])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Top,
        ElasticMaterialField3::ThermoelasticByCell(&field),
        &|p| p[0] == 0.0,
        &|_, n| n[2] > 0.5,
        ElasticityOptions3 {
            formulation: ElasticityFormulation3::BBar,
            ghost_gamma: 2.7,
            ..Default::default()
        },
        Default::default(),
        Default::default(),
        &mut q,
    )
    .unwrap()
}
fn study<O: Sdf3Elasticity>(mut op: O) -> CutDensityStudy3<O> {
    op.set_scales(&vec![0.83; op.cells()]).unwrap();
    let body = op
        .prepare_material_body_load(&|p| [0.004 * (1.0 + p[1]), 0.001, -0.003], 100_000, poll)
        .unwrap();
    let mut regions = vec![PhysicalRegion3::Design; op.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = PhysicalRegion3::Void;
    CutDensityStudy3::new(op, 0.12, params())
        .with_material_body_load(body)
        .unwrap()
        .with_physical_regions(regions)
        .unwrap()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-10),
        "{a:e} != {b:e}"
    );
}
fn motion(p: [f64; 3], delta: f64) -> [f64; 3] {
    [
        0.0002 * p[0] + delta * 0.000001 * p[0],
        -0.0001 * p[0] * p[2],
        0.0001 * p[0],
    ]
}
fn objective(
    u: &[f64],
    points: &[ThermoelasticResponse3],
    rho: &[f64],
    force: &[f64],
) -> Result<CutThermoelasticObjective3, ElasticityError3> {
    let mut value = dot(u, force) + 0.5 * dot(u, u) + 1e-4 * force.iter().sum::<f64>();
    let displacement_derivative = u.iter().zip(force).map(|(u, f)| u + f).collect();
    let mut stress = Vec::with_capacity(points.len());
    let mut temperature = vec![0.0; rho.len()];
    let mut seen = vec![false; rho.len()];
    for point in points {
        let derivative =
            core::array::from_fn(|i| 0.0001 * point.mechanical.weight * point.mechanical.stress[i]);
        value += 0.5 * dot(&derivative, &point.mechanical.stress);
        stress.push(derivative);
        let cell = point.mechanical.cell;
        if !seen[cell] {
            let dt = point.temperature_k - 300.0;
            value += 0.5e-8 * dt * dt;
            temperature[cell] = 1e-8 * dt;
            seen[cell] = true;
        }
    }
    let mut scales = Vec::with_capacity(rho.len());
    for &r in rho {
        let scale = 0.03 + 0.97 * r * r * r;
        value += 0.5e-5 * r * r + 0.5e-5 * scale * scale;
        scales.push(1e-5 * scale);
    }
    Ok(CutThermoelasticObjective3 {
        native: ThermoelasticDesignObjective3 {
            objective: ThermoelasticObjective3 {
                value,
                displacement_derivative,
                stress_derivatives: Some(stress),
                explicit_temperature_derivative: Some(temperature),
            },
            explicit_scale_derivative: Some(scales),
        },
        explicit_occupancy_derivative: Some(rho.iter().map(|r| 1e-5 * r).collect()),
        explicit_force_derivative: Some(u.iter().map(|u| u + 1e-4).collect()),
    })
}
fn evaluate<O: ThermoelasticSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    rho: &[f64],
    force: &[f64],
    delta: f64,
) -> CutThermoelasticEvaluation3 {
    let g = |p, n| {
        let _ = n;
        motion(p, delta)
    };
    let gt = |p: [f64; 3], _| [0.000001 * p[0], 0.0, 0.0];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let calls = Cell::new(0);
    let result = study
        .evaluate_thermoelastic(
            rho,
            force,
            Some(ThermoelasticMotion3 {
                displacement: &g,
                temperature_derivative: Some(&gt),
            }),
            settings(),
            |u, p, r, f| {
                calls.set(calls.get() + 1);
                objective(u, p, r, f)
            },
            &mut control,
        )
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert!(result.work.linear_solves >= 3);
    assert!(result.sensitivity.thermoelastic.solution.iterations() > 0);
    assert!(study.operator().scales().iter().all(|s| *s == 0.83));
    result
}
fn check_chain<O: ThermoelasticSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
) -> (Vec<f64>, CutThermoelasticEvaluation3) {
    let rho: Vec<_> = (0..study.cells())
        .map(|i| 0.44 + 0.012 * (i % 8) as f64)
        .collect();
    let force = vec![0.0; study.operator().n()];
    let base = evaluate(study, &rho, &force, 0.0);
    assert_eq!(base.projected_rho[0], 1.0);
    assert_eq!(*base.projected_rho.last().unwrap(), 0.0);
    for (i, &value) in base.force_derivative.iter().enumerate() {
        if study.operator().fixed()[i / 3] {
            assert_eq!(value, 0.0)
        } else {
            close(
                value,
                base.sensitivity.thermoelastic.force_derivative[i]
                    + base.sensitivity.thermoelastic.solution.coefficients()[i]
                    + 1e-4,
                2e-13,
            )
        }
    }
    for direction in [
        vec![1.0; rho.len()],
        (0..rho.len()).map(|i| (i % 3) as f64 - 1.0).collect(),
        (0..rho.len())
            .map(|i| if i == 0 { 1.0 } else { 0.0 })
            .collect(),
    ] {
        let h = 2e-5;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = evaluate(study, &plus, &force, 0.0);
        let b = evaluate(study, &minus, &force, 0.0);
        close(
            (a.sensitivity.thermoelastic.value - b.sensitivity.thermoelastic.value) / (2.0 * h),
            dot(&base.gradient, &direction),
            5e-4,
        );
        close(
            (a.volume_fraction - b.volume_fraction) / (2.0 * h),
            dot(&base.volume_gradient, &direction),
            5e-5,
        );
    }
    // A raw control under protected material still influences adjacent design cells.
    assert!(base.volume_gradient[0] > 0.0);
    let direction: Vec<_> = (0..force.len())
        .map(|i| 0.0001 * ((i % 5) as f64 - 2.0))
        .collect();
    let h = 1e-4;
    let a = evaluate(
        study,
        &rho,
        &direction.iter().map(|v| h * v).collect::<Vec<_>>(),
        0.0,
    );
    let b = evaluate(
        study,
        &rho,
        &direction.iter().map(|v| -h * v).collect::<Vec<_>>(),
        0.0,
    );
    close(
        (a.sensitivity.thermoelastic.value - b.sensitivity.thermoelastic.value) / (2.0 * h),
        dot(&base.force_derivative, &direction),
        5e-4,
    );
    (rho, base)
}

#[test]
fn cartesian_density_temperature_and_authored_force_share_the_native_objective() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut base_study = study(cartesian(0.0, formulation));
        let (rho, base) = check_chain(&mut base_study);
        let h = 1e-4;
        let mut plus = study(cartesian(h, formulation));
        let mut minus = study(cartesian(-h, formulation));
        let force = vec![0.0; plus.operator().n()];
        let a = evaluate(&mut plus, &rho, &force, h);
        let b = evaluate(&mut minus, &rho, &force, -h);
        close(
            (a.sensitivity.thermoelastic.value - b.sensitivity.thermoelastic.value) / (2.0 * h),
            base.sensitivity
                .thermoelastic
                .temperature_derivative
                .iter()
                .sum(),
            5e-4,
        );
    }
}
#[test]
fn adaptive_hanging_field_uses_the_same_density_and_body_force_transposes() {
    let mut study = study(adaptive());
    check_chain(&mut study);
}
#[test]
fn failures_and_unwinding_callbacks_restore_the_incoming_scales() {
    let mut study = study(cartesian(0.0, ElasticityFormulation3::BBar));
    let rho = vec![0.5; study.cells()];
    let force = vec![0.0; study.operator().n()];
    let original = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    assert!(
        study
            .evaluate_thermoelastic(
                &rho,
                &force,
                None,
                settings(),
                |_, _, _, _| Err(ElasticityError3::Invalid("caller refuses objective")),
                &mut control
            )
            .is_err()
    );
    assert_eq!(study.operator().scales(), original);
    assert!(
        study
            .evaluate_thermoelastic(
                &rho,
                &force,
                None,
                settings(),
                |u, p, r, f| {
                    let mut value = objective(u, p, r, f)?;
                    value.explicit_occupancy_derivative = Some(vec![]);
                    Ok(value)
                },
                &mut control
            )
            .is_err()
    );
    assert_eq!(study.operator().scales(), original);
    let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = study.evaluate_thermoelastic(
            &rho,
            &force,
            None,
            settings(),
            |_, _, _, _| panic!("authored callback panic"),
            &mut control,
        );
    }));
    assert!(unwound.is_err());
    assert_eq!(study.operator().scales(), original);
    let cancelled = Cell::new(false);
    let mut callback = |_| {
        if cancelled.get() {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let result = study.evaluate_thermoelastic(
        &rho,
        &force,
        None,
        settings(),
        |u, p, r, f| {
            cancelled.set(true);
            objective(u, p, r, f)
        },
        &mut control,
    );
    assert!(matches!(
        result,
        Err(CutThermoelasticError3::Evaluation(
            EvaluationStop::Cancelled
        ))
    ));
    assert_eq!(study.operator().scales(), original);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let no_mechanics = ThermoelasticAdjointSettings3 {
        maximum_iterations: 0,
        ..settings()
    };
    assert!(
        study
            .evaluate_thermoelastic(&rho, &force, None, no_mechanics, objective, &mut control)
            .is_err()
    );
    assert_eq!(study.operator().scales(), original);
}
