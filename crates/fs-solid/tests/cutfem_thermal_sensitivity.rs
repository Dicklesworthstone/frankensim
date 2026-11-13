//! G0/G1/G3: temperature transposes of refreshed native CutFEM mechanics.
//! Authored source curves verify numerical mechanics, not material calibration.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::thermal::ThermalJacobianLimits3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3, MaterialRefreshLimits3,
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
use fs_material::elastic::thermal::{
    AffineIsotropicThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::op::LinearOp;

const POINTS: usize = 100_000;
const TOTAL_STRAIN: f64 = 0.1;

fn jacobian_limits() -> ThermalJacobianLimits3 {
    ThermalJacobianLimits3 {
        max_points: POINTS,
        max_workspace_components: 1_000_000,
        max_action_terms: 1_000_000,
    }
}

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

struct Top(f64);
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - self.0, hi[2] - self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}

fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    action(&mut control)
}

fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}

fn law(parameters: [f64; 5]) -> AffineIsotropicThermoelastic {
    let [poisson, young_slope, nu_slope, alpha, alpha_slope] = parameters;
    let temperature_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0, young_slope),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), poisson, nu_slope),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            alpha,
            alpha_slope,
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
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity,
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored CutFEM temperature-adjoint fixture".into(),
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
        .with_quantity("T", temperature_quantity, 300.0)
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

fn hot(law: &AffineIsotropicThermoelastic, temperature: f64) -> IsothermalThermoelasticState {
    law.at_temperature(temperature)
        .unwrap()
        .isothermal_elasticity(IsothermalThermoelasticInput {
            reference_density_kg_m3: 3.0,
            mechanical_material_identity: ContentHash([73; 32]),
        })
        .unwrap()
}

fn temperature(cell: HexCell) -> f64 {
    let lo = cell.lo();
    316.0 + 3.0 * lo[0] + 2.0 * lo[1] + lo[2]
}

fn cartesian(
    law: &AffineIsotropicThermoelastic,
    formulation: ElasticityFormulation3,
) -> CutElasticity3 {
    let materials = |cell| Ok(hot(law, temperature(cell)));
    control(|control| {
        CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            [2; 3],
            &Top(0.73),
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 {
                formulation,
                ghost_gamma: 2.7,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    })
}

fn fields(nodes: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>) {
    let u = nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.001 * x * x * (1.0 + y),
                -0.0008 * x * (z * z + 0.2 * y),
                0.0006 * x * y * y,
            ]
        })
        .collect();
    // Algebraic dual strains deliberately exceed the source's physical limit.
    let dual = nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                1.7 * x * (1.0 + 0.3 * x + z),
                -0.4 * (1.0 + x * y),
                0.9 * x * z * z,
            ]
        })
        .collect();
    (u, dual)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}",
    );
}

fn cartesian_residual(op: &CutElasticity3, u: &[f64]) -> Vec<f64> {
    let mut residual = vec![0.0; op.n()];
    op.apply(u, &mut residual);
    let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let free = op.eigenstrain_load(&source, poll).unwrap();
    for (value, force) in residual.iter_mut().zip(free) {
        *value -= force;
    }
    residual
}

fn support_motion([x, y, z]: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [0.0003 * x * (1.0 + y), -0.0002 * x * z, 0.0004 * x]
}

fn support_motion_slope([x, y, z]: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [0.00001 * x, -0.000006 * x * y, 0.000004 * x * z]
}

fn subtract_motion(
    op: &CutElasticity3,
    residual: &mut [f64],
    motion: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
) {
    for (value, load) in residual
        .iter_mut()
        .zip(op.prescribed_displacement_load(motion, poll).unwrap())
    {
        *value -= load;
    }
}

#[test]
fn full_temperature_pullback_matches_each_refreshed_cell_with_nitsche_and_ghosts() {
    for (formulation, parameters) in [
        (
            ElasticityFormulation3::Displacement,
            [0.2, -2.0, 0.0008, 0.0001, 0.000002],
        ),
        (
            ElasticityFormulation3::BBar,
            // Poisson slope 2^-16: the knots 0.49 - 20 s and 0.49 + 60 s
            // are then exact, so the path is affine in binary as the source
            // check requires (one ulp off the line is a kink; 2e-5 was).
            [0.49, -2.0, 1.525_878_906_25e-5, 0.0001, 0.000002],
        ),
    ] {
        for with_motion in [false, true] {
            let law = law(parameters);
            let at = hot(&law, 317.0);
            let c = at.acoustic().stiffness_mandel_pa();
            let dc = at
                .source()
                .temperature_derivatives()
                .unwrap()
                .stiffness_mandel_pa_k;
            // Changing Poisson response cannot be replaced by a scalar dE/E.
            assert!((dc[0][0] * c[3][3] - dc[3][3] * c[0][0]).abs() > 1.0);
            let mut op = cartesian(&law, formulation);
            let scales: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.04 * i as f64).collect();
            op.set_scales(&scales).unwrap();
            let (u, dual) = fields(op.nodes());
            let gradient = if with_motion {
                op.thermal_residual_temperature_pullback_with_motion(
                    &u,
                    &dual,
                    &support_motion,
                    None,
                    TOTAL_STRAIN,
                    POINTS,
                    poll,
                )
                .unwrap()
            } else {
                op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
                    .unwrap()
            };
            let lower: Vec<_> = op
                .cell_nodes()
                .iter()
                .map(|nodes| op.nodes()[nodes[0]])
                .collect();
            let original_mass = op.reference_cell_masses_kg().unwrap().to_vec();
            let h = 0.001;
            for (cell, &target) in lower.iter().enumerate() {
                let mut values = [0.0; 2];
                for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
                    let materials = |bounds: HexCell| {
                        let step = if bounds.lo() == target { sign * h } else { 0.0 };
                        Ok(hot(&law, temperature(bounds) + step))
                    };
                    op.refresh_material_field(
                        ElasticMaterialField3::ThermoelasticByCell(&materials),
                        MaterialRefreshLimits3::default(),
                        poll,
                    )
                    .unwrap();
                    let mut residual = cartesian_residual(&op, &u);
                    if with_motion {
                        subtract_motion(&op, &mut residual, &support_motion);
                    }
                    values[sample] = dot(&dual, &residual);
                    assert_eq!(op.scales(), scales.as_slice());
                    assert_eq!(
                        op.reference_cell_masses_kg().unwrap(),
                        original_mass.as_slice()
                    );
                }
                close(
                    gradient[cell],
                    (values[1] - values[0]) / (2.0 * h),
                    3e-6,
                    2e-8,
                );
            }
        }
    }
}

#[test]
fn adaptive_temperature_pullback_uses_the_original_hanging_map_on_both_fields() {
    for with_motion in [false, true] {
        let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
        let tree = Octree3::uniform(1, 4, 4096).unwrap();
        let tree = tree
            .refined(&[*tree.leaves().iter().next().unwrap()], poll)
            .unwrap();
        let materials = |cell| Ok(hot(&law, temperature(cell)));
        let mut op = control(|control| {
            AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
                domain(),
                &tree,
                &Top(0.73),
                ElasticMaterialField3::ThermoelasticByCell(&materials),
                &|p| p[0] == 0.0,
                &|_, normal| normal[2] > 0.5,
                ElasticityOptions3 {
                    ghost_gamma: 2.7,
                    ..Default::default()
                },
                Default::default(),
                Default::default(),
                control,
            )
            .unwrap()
        });
        assert!(op.physical_nodes().len() > op.nodes().len());
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.02 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let (u, dual) = fields(op.nodes());
        let gradient = if with_motion {
            op.thermal_residual_temperature_pullback_with_motion(
                &u,
                &dual,
                &support_motion,
                None,
                TOTAL_STRAIN,
                POINTS,
                poll,
            )
            .unwrap()
        } else {
            op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
                .unwrap()
        };
        let lower: Vec<_> = op
            .cell_nodes()
            .iter()
            .map(|nodes| op.physical_nodes()[nodes[0]])
            .collect();
        let h = 0.001;
        for (cell, &target) in lower.iter().enumerate() {
            let mut values = [0.0; 2];
            for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let materials = |bounds: HexCell| {
                    let step = if bounds.lo() == target { sign * h } else { 0.0 };
                    Ok(hot(&law, temperature(bounds) + step))
                };
                op.refresh_material_field(
                    ElasticMaterialField3::ThermoelasticByCell(&materials),
                    MaterialRefreshLimits3::default(),
                    poll,
                )
                .unwrap();
                let mut residual = vec![0.0; op.n()];
                op.apply(&u, &mut residual);
                let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
                for (value, force) in residual
                    .iter_mut()
                    .zip(op.eigenstrain_load(&source, poll).unwrap())
                {
                    *value -= force;
                }
                if with_motion {
                    for (value, load) in residual.iter_mut().zip(
                        op.prescribed_displacement_load(&support_motion, poll)
                            .unwrap(),
                    ) {
                        *value -= load;
                    }
                }
                values[sample] = dot(&dual, &residual);
            }
            close(
                gradient[cell],
                (values[1] - values[0]) / (2.0 * h),
                4e-6,
                3e-8,
            );
        }
    }
}

fn plain(state: &IsothermalThermoelasticState, counts: [usize; 3], top: f64) -> CutElasticity3 {
    control(|control| {
        CutElasticity3::build_with_material_field(
            domain(),
            counts,
            &Top(top),
            ElasticMaterialField3::ThermoelasticUniform(state),
            &|p| p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    })
}

#[test]
fn nonsmooth_native_stabilization_refuses_only_when_used() {
    let kink = hot(&law([0.0, 0.0, 0.001, 0.0001, 0.0]), 300.0);
    let full = plain(&kink, [1; 3], 1.2); // No Nitsche or ghost modulus.
    let (u, dual) = fields(full.nodes());
    full.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
        .unwrap();

    let cut = plain(&kink, [2; 3], 0.73);
    let (u, dual) = fields(cut.nodes());
    let error = cut
        .thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
        .unwrap_err();
    assert!(matches!(error, ElasticityError3::Invalid(message) if message.contains("nonsmooth")));

    // Zero C' and equal active row slopes are smooth, despite zero couplings.
    let smooth = hot(&law([0.0, 0.0, 0.0, 0.0001, 0.000002]), 300.0);
    let cut = plain(&smooth, [2; 3], 0.73);
    let (u, dual) = fields(cut.nodes());
    let derivative = cut
        .thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
        .unwrap();
    assert!(derivative.iter().any(|value| value.abs() > 1e-8));
}

#[test]
fn physical_domain_budgets_and_late_cancellation_leave_the_owner_unchanged() {
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let (u, dual) = fields(op.nodes());
    let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let before = cartesian_residual(&op, &u);
    let mut calls = 0usize;
    let expected = op
        .thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, || {
            calls += 1;
            poll()
        })
        .unwrap();
    assert!(calls > 10);
    for stop in [1, calls / 2, calls] {
        let mut at = 0usize;
        assert!(matches!(
            op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, || {
                at += 1;
                if at == stop {
                    ControlFlow::Break(())
                } else {
                    poll()
                }
            },),
            Err(ElasticityError3::Cancelled),
        ));
    }
    assert!(
        op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, 0, poll,)
            .is_err()
    );
    assert!(
        op.thermal_residual_temperature_pullback(&u, &dual, 1e-6, POINTS, poll,)
            .is_err()
    );
    let invalid: Vec<_> = u.iter().map(|value| value * 100.0).collect();
    assert!(
        op.thermal_residual_temperature_pullback(&invalid, &dual, 0.9, POINTS, poll,)
            .is_err()
    ); // Native total/elastic source bounds remain stricter.
    let large_dual: Vec<_> = dual.iter().map(|value| value * 1e4).collect();
    let actual = op
        .thermal_residual_temperature_pullback(&u, &large_dual, TOTAL_STRAIN, POINTS, poll)
        .unwrap();
    for (actual, expected) in actual.iter().zip(&expected) {
        close(*actual, 1e4 * expected, 3e-12, 1e-8);
    }
    op.validate_eigenstrain(&source).unwrap();
    assert_eq!(cartesian_residual(&op, &u), before);
    assert_eq!(
        op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll,)
            .unwrap(),
        expected
    );

    let state = hot(&law, 316.0);
    let cold = control(|control| {
        CutElasticity3::build_with_material_field(
            domain(),
            [1; 3],
            &Top(1.2),
            ElasticMaterialField3::Uniform(state.acoustic()),
            &|p| p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    });
    let (u, dual) = fields(cold.nodes());
    assert!(
        cold.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll,)
            .is_err()
    );
}

#[test]
fn temperature_dependent_motion_includes_the_original_lifting_product_rule() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
        let mut op = cartesian(&law, formulation);
        let at = 317.0;
        let reference = hot(&law, at);
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticUniform(&reference),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.03 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let (u, dual) = fields(op.nodes());
        let fixed = op
            .thermal_residual_temperature_pullback_with_motion(
                &u,
                &dual,
                &support_motion,
                None,
                TOTAL_STRAIN,
                POINTS,
                poll,
            )
            .unwrap();
        let moving = op
            .thermal_residual_temperature_pullback_with_motion(
                &u,
                &dual,
                &support_motion,
                Some(&support_motion_slope),
                TOTAL_STRAIN,
                POINTS,
                poll,
            )
            .unwrap();
        let motion_load = op
            .prescribed_displacement_load(&support_motion_slope, poll)
            .unwrap();
        close(
            moving.iter().sum::<f64>() - fixed.iter().sum::<f64>(),
            -dot(&dual, &motion_load),
            2e-11,
            2e-11,
        );
        assert!(dot(&dual, &motion_load).abs() > 1e-5);
        let h = 0.001;
        let mut observed = [0.0; 2];
        for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let source = hot(&law, at + sign * h);
            op.refresh_material_field(
                ElasticMaterialField3::ThermoelasticUniform(&source),
                MaterialRefreshLimits3::default(),
                poll,
            )
            .unwrap();
            let displaced = |position, normal| {
                let base = support_motion(position, normal);
                let slope = support_motion_slope(position, normal);
                std::array::from_fn(|c| base[c] + sign * h * slope[c])
            };
            let mut residual = cartesian_residual(&op, &u);
            subtract_motion(&op, &mut residual, &displaced);
            observed[sample] = dot(&dual, &residual);
        }
        close(
            moving.iter().sum(),
            (observed[1] - observed[0]) / (2.0 * h),
            4e-6,
            3e-8,
        );
    }
}

#[test]
fn motion_pullback_admission_and_callback_cancellation_publish_no_partial_state() {
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let (u, dual) = fields(op.nodes());
    let before = cartesian_residual(&op, &u);
    let preparation = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let calls = std::cell::Cell::new(0usize);
    let counted_motion = |p, n| {
        calls.set(calls.get() + 1);
        support_motion(p, n)
    };
    assert!(
        op.thermal_residual_temperature_pullback_with_motion(
            &u,
            &dual,
            &counted_motion,
            None,
            TOTAL_STRAIN,
            0,
            poll,
        )
        .is_err()
    );
    assert_eq!(
        calls.get(),
        0,
        "point admission precedes boundary callbacks"
    );
    let bad = |_, _| [f64::NAN, 0.0, 0.0];
    for (motion, derivative) in [
        (&bad as &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3], None),
        (
            &support_motion,
            Some(&bad as &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]),
        ),
    ] {
        assert!(
            op.thermal_residual_temperature_pullback_with_motion(
                &u,
                &dual,
                motion,
                derivative,
                TOTAL_STRAIN,
                POINTS,
                poll,
            )
            .is_err()
        );
    }
    calls.set(0);
    assert!(matches!(
        op.thermal_residual_temperature_pullback_with_motion(
            &u,
            &dual,
            &counted_motion,
            None,
            TOTAL_STRAIN,
            POINTS,
            || if calls.get() > 0 {
                ControlFlow::Break(())
            } else {
                poll()
            },
        ),
        Err(ElasticityError3::Cancelled),
    ));
    op.validate_eigenstrain(&preparation).unwrap();
    assert_eq!(cartesian_residual(&op, &u), before);
    let homogeneous = op
        .thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll)
        .unwrap();
    let zero = op
        .thermal_residual_temperature_pullback_with_motion(
            &u,
            &dual,
            &|_, _| [0.0; 3],
            None,
            TOTAL_STRAIN,
            POINTS,
            poll,
        )
        .unwrap();
    assert_eq!(zero, homogeneous);
    let state = hot(&law, 317.0);
    let no_embedded = plain(&state, [1; 3], 1.2);
    let (u, dual) = fields(no_embedded.nodes());
    assert!(matches!(
        no_embedded.thermal_residual_temperature_pullback_with_motion(
            &u, &dual, &support_motion, None, TOTAL_STRAIN, POINTS, poll,
        ),
        Err(ElasticityError3::Invalid(message)) if message.contains("embedded"),
    ));
}

fn adjoint_settings() -> fs_cutfem::elastic3::thermal::ThermoelasticAdjointSettings3 {
    fs_cutfem::elastic3::thermal::ThermoelasticAdjointSettings3 {
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 5000,
        poll_iterations: 8,
        maximum_kinematic_strain: TOTAL_STRAIN,
        maximum_points: POINTS,
    }
}

// A nonlinear objective with independent displacement, physical stress and
// explicit temperature dependence. Coefficients define the fixture's units.
fn thermal_objective(
    u: &[f64],
    points: &[fs_cutfem::elastic3::thermal::ThermoelasticResponse3],
    cells: usize,
) -> Result<fs_cutfem::elastic3::thermal::ThermoelasticObjective3, ElasticityError3> {
    let displacement: Vec<_> = u
        .iter()
        .enumerate()
        .map(|(i, &u)| (1.0 + 0.3 * (i % 3) as f64) * u)
        .collect();
    let mut value = 0.5 * dot(u, &displacement);
    let mut stress = Vec::with_capacity(points.len());
    let mut explicit = vec![0.0; cells];
    let mut seen = vec![false; cells];
    for point in points {
        let derivative = std::array::from_fn(|component| {
            1e-3 * point.mechanical.weight
                * (1.0 + 0.2 * component as f64)
                * point.mechanical.stress[component]
        });
        value += 0.5 * dot(&derivative, &point.mechanical.stress);
        stress.push(derivative);
        let cell = point.mechanical.cell;
        if !seen[cell] {
            let change = point.temperature_k - 300.0;
            explicit[cell] = 1e-6 * change;
            value += 0.5e-6 * change * change;
            seen[cell] = true;
        }
    }
    Ok(fs_cutfem::elastic3::thermal::ThermoelasticObjective3 {
        value,
        displacement_derivative: displacement,
        stress_derivatives: Some(stress),
        explicit_temperature_derivative: Some(explicit),
    })
}

#[test]
fn complete_thermal_objective_adjoint_matches_refreshed_equilibria_and_external_forces() {
    use fs_cutfem::elastic3::thermal::ThermoelasticMotion3;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
        let mut op = cartesian(&law, formulation);
        let baseline = hot(&law, 317.0);
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticUniform(&baseline),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.03 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
        let external: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|&[x, y, z]| [0.003 * x, -0.002 * x * y, 0.001 * x * z])
            .collect();
        let cells = op.cells();
        let settings = adjoint_settings();
        let result = op
            .solve_thermoelastic_adjoint(
                &external,
                Some(ThermoelasticMotion3 {
                    displacement: &support_motion,
                    temperature_derivative: Some(&support_motion_slope),
                }),
                settings,
                |u, points| thermal_objective(u, points, cells),
                |_| poll(),
            )
            .unwrap();
        assert!(result.solution.iterations() > 0);
        assert!(
            result.solution.iterations() + result.adjoint_iterations <= settings.maximum_iterations
        );
        assert!(
            result.adjoint_residual.euclidean().unwrap() < settings.relative_residual_tolerance
        );
        assert_eq!(op.scales(), scales);
        assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory);
        let h = 0.002;
        let mut values = [0.0; 2];
        for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let state = hot(&law, 317.0 + sign * h);
            op.refresh_material_field(
                ElasticMaterialField3::ThermoelasticUniform(&state),
                MaterialRefreshLimits3::default(),
                poll,
            )
            .unwrap();
            let moved = |p, n| {
                let u = support_motion(p, n);
                let slope = support_motion_slope(p, n);
                std::array::from_fn(|c| u[c] + sign * h * slope[c])
            };
            let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
            let mut force = op.eigenstrain_load(&prepared, poll).unwrap();
            let lift = op.prescribed_displacement_load(&moved, poll).unwrap();
            for ((force, external), lift) in force.iter_mut().zip(&external).zip(lift) {
                *force += external + lift;
            }
            let primal = op
                .solve_controlled(&force, 1e-11, 5000, 8, |_| poll())
                .unwrap();
            let points = op
                .thermoelastic_response(primal.coefficients(), TOTAL_STRAIN, POINTS, poll)
                .unwrap();
            values[sample] = thermal_objective(primal.coefficients(), &points, cells)
                .unwrap()
                .value;
        }
        close(
            result.temperature_derivative.iter().sum(),
            (values[1] - values[0]) / (2.0 * h),
            3e-5,
            2e-9,
        );

        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticUniform(&baseline),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let direction: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|&[x, y, z]| [x * (1.0 + y), -x * z, 0.4 * x])
            .collect();
        let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
        let free = op.eigenstrain_load(&prepared, poll).unwrap();
        let lift = op
            .prescribed_displacement_load(&support_motion, poll)
            .unwrap();
        for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let force: Vec<_> = external
                .iter()
                .zip(&free)
                .zip(&lift)
                .zip(&direction)
                .map(|(((external, free), lift), direction)| {
                    external + free + lift + sign * h * direction
                })
                .collect();
            let primal = op
                .solve_controlled(&force, 1e-11, 5000, 8, |_| poll())
                .unwrap();
            let points = op
                .thermoelastic_response(primal.coefficients(), TOTAL_STRAIN, POINTS, poll)
                .unwrap();
            values[sample] = thermal_objective(primal.coefficients(), &points, cells)
                .unwrap()
                .value;
        }
        close(
            dot(&result.force_derivative, &direction),
            (values[1] - values[0]) / (2.0 * h),
            3e-5,
            2e-9,
        );
        assert!(matches!(
            op.solve_thermoelastic_adjoint(
                &external,
                Some(ThermoelasticMotion3 {
                    displacement: &support_motion,
                    temperature_derivative: Some(&support_motion_slope),
                }),
                fs_cutfem::elastic3::thermal::ThermoelasticAdjointSettings3 {
                    maximum_iterations: result.solution.iterations(),
                    ..settings
                },
                |u, points| thermal_objective(u, points, cells),
                |_| poll(),
            ),
            Err(ElasticityError3::NotConverged { .. })
        ));
    }
}

#[test]
fn complete_thermal_adjoint_preserves_original_adaptive_hanging_reconstruction() {
    use fs_cutfem::elastic3::thermal::ThermoelasticMotion3;
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let materials = |bounds| Ok(hot(&law, temperature(bounds)));
    let mut op = control(|control| {
        AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            &tree,
            &Top(0.73),
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 {
                ghost_gamma: 2.7,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    });
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.02 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let external: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| [0.003 * x, -0.002 * x * y, 0.001 * x * z])
        .collect();
    let cells = op.cells();
    let result = op
        .solve_thermoelastic_adjoint(
            &external,
            Some(ThermoelasticMotion3 {
                displacement: &support_motion,
                temperature_derivative: None,
            }),
            adjoint_settings(),
            |u, points| thermal_objective(u, points, cells),
            |_| poll(),
        )
        .unwrap();
    let slope = |p: [f64; 3]| 0.4 + 0.5 * p[0] - 0.2 * p[1] + 0.3 * p[2];
    let directions: Vec<_> = op
        .cell_nodes()
        .iter()
        .map(|nodes| slope(op.physical_nodes()[nodes[0]]))
        .collect();
    let h = 0.002;
    let mut values = [0.0; 2];
    for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let changed = |bounds: HexCell| {
            Ok(hot(
                &law,
                temperature(bounds) + sign * h * slope(bounds.lo()),
            ))
        };
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticByCell(&changed),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
        let mut force = op.eigenstrain_load(&prepared, poll).unwrap();
        let lift = op
            .prescribed_displacement_load(&support_motion, poll)
            .unwrap();
        for ((force, external), lift) in force.iter_mut().zip(&external).zip(lift) {
            *force += external + lift;
        }
        let primal = op
            .solve_controlled(&force, 1e-11, 5000, 8, |_| poll())
            .unwrap();
        let points = op
            .thermoelastic_response(primal.coefficients(), TOTAL_STRAIN, POINTS, poll)
            .unwrap();
        values[sample] = thermal_objective(primal.coefficients(), &points, cells)
            .unwrap()
            .value;
    }
    close(
        dot(&result.temperature_derivative, &directions),
        (values[1] - values[0]) / (2.0 * h),
        5e-5,
        3e-9,
    );
}

#[test]
fn complete_thermal_adjoint_admits_budgets_and_cancels_after_the_objective() {
    use std::cell::Cell;
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let force = vec![0.0; op.n()];
    let called = Cell::new(false);
    let settings = adjoint_settings();
    let original = op.scales().to_vec();
    assert!(matches!(
        op.solve_thermoelastic_adjoint(
            &force,
            None,
            fs_cutfem::elastic3::thermal::ThermoelasticAdjointSettings3 {
                maximum_points: 0,
                ..settings
            },
            |u, points| {
                called.set(true);
                thermal_objective(u, points, op.cells())
            },
            |_| poll(),
        ),
        Err(ElasticityError3::Invalid(_))
    ));
    assert!(!called.get());
    assert!(matches!(
        op.solve_thermoelastic_adjoint(
            &force,
            None,
            settings,
            |u, points| {
                called.set(true);
                thermal_objective(u, points, op.cells())
            },
            |_| if called.get() {
                ControlFlow::Break(())
            } else {
                poll()
            },
        ),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(called.get());
    assert!(matches!(
        op.solve_thermoelastic_adjoint(
            &force,
            None,
            settings,
            |u, points| {
                let mut goal = thermal_objective(u, points, op.cells())?;
                goal.stress_derivatives.as_mut().unwrap()[0][0] = f64::NAN;
                Ok(goal)
            },
            |_| poll(),
        ),
        Err(ElasticityError3::Invalid(_))
    ));
    assert_eq!(op.scales(), original);
}

fn temperature_direction([x, y, z]: [f64; 3]) -> f64 {
    0.4 + 0.8 * x - 1.2 * y + 0.5 * z
}

fn local_direction(bounds: &[([f64; 3], [f64; 3])], point: [f64; 3]) -> f64 {
    let (lo, _) = bounds.iter().find(|(lo, hi)| {
        (0..3).all(|a| lo[a] <= point[a] && point[a] < hi[a])
    }).expect("selected quadrature point is interior to its owning leaf");
    temperature_direction(*lo)
}

#[test]
fn prepared_temperature_block_matches_full_residual_jvp_and_direct_vjp() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        for with_motion in [false, true] {
            let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
            let mut op = cartesian(&law, formulation);
            let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.03 * i as f64).collect();
            op.set_scales(&scales).unwrap();
            let (u, dual) = fields(op.nodes());
            let bounds: Vec<_> = op.cell_nodes().iter()
                .map(|nodes| (op.nodes()[nodes[0]], op.nodes()[nodes[7]])).collect();
            let direction: Vec<_> = bounds.iter().map(|(lo, _)| temperature_direction(*lo)).collect();
            let calls = std::cell::Cell::new(0usize);
            let counted_motion = |p, n| {
                calls.set(calls.get() + 1);
                support_motion(p, n)
            };
            let prepared = if with_motion {
                op.prepare_thermal_residual_jacobian_with_motion(&u, &counted_motion,
                    Some(&support_motion_slope), TOTAL_STRAIN, jacobian_limits(), poll).unwrap()
            } else {
                op.prepare_thermal_residual_jacobian(&u, TOTAL_STRAIN, jacobian_limits(), poll).unwrap()
            };
            assert_eq!(prepared.n(), op.n());
            assert_eq!(prepared.temperatures(), op.cells());
            let samples = calls.get();
            if with_motion { assert!(samples > 0); }
            let mut jvp = vec![91.0; op.n()];
            let mut vjp = vec![92.0; op.cells()];
            op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut jvp, poll).unwrap();
            op.thermal_residual_temperature_transpose_into(&prepared, &dual, &mut vjp, poll).unwrap();
            assert_eq!(calls.get(), samples, "reuse never re-queries the prescribed boundary");
            let direct = if with_motion {
                op.thermal_residual_temperature_pullback_with_motion(&u, &dual, &support_motion,
                    Some(&support_motion_slope), TOTAL_STRAIN, POINTS, poll).unwrap()
            } else {
                op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll).unwrap()
            };
            for (actual, expected) in vjp.iter().zip(&direct) {
                close(*actual, *expected, 3e-11, 1e-11);
            }
            close(dot(&dual, &jvp), dot(&direction, &vjp), 3e-12, 1e-12);
            for (node, &fixed) in op.fixed().iter().enumerate() {
                if fixed { assert_eq!(&jvp[3 * node..3 * node + 3], &[0.0; 3]); }
            }
            let h = 0.001;
            let mut residuals = Vec::new();
            for sign in [-1.0, 1.0] {
                let materials = |cell: HexCell| Ok(hot(&law,
                    temperature(cell) + sign * h * temperature_direction(cell.lo())));
                op.refresh_material_field(ElasticMaterialField3::ThermoelasticByCell(&materials),
                    MaterialRefreshLimits3::default(), poll).unwrap();
                let mut residual = cartesian_residual(&op, &u);
                if with_motion {
                    let displaced = |p, n| {
                        let base = support_motion(p, n);
                        let slope = support_motion_slope(p, n);
                        let delta = sign * h * local_direction(&bounds, p);
                        std::array::from_fn(|c| base[c] + delta * slope[c])
                    };
                    subtract_motion(&op, &mut residual, &displaced);
                }
                residuals.push(residual);
            }
            for i in 0..op.n() {
                close(jvp[i], (residuals[1][i] - residuals[0][i]) / (2.0 * h), 4e-6, 3e-8);
            }
        }
    }
}

#[test]
fn adaptive_prepared_temperature_block_preserves_hanging_transpose_and_moving_supports() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
        let tree = Octree3::uniform(1, 4, 4096).unwrap();
        let tree = tree.refined(&[*tree.leaves().iter().next().unwrap()], poll).unwrap();
        let materials = |cell| Ok(hot(&law, temperature(cell)));
        let mut op = control(|control| AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(), &tree, &Top(0.73), ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0, &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 { formulation, ghost_gamma: 2.7, ..Default::default() },
            Default::default(), Default::default(), control).unwrap());
        assert!(op.physical_nodes().len() > op.nodes().len());
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.02 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let (u, dual) = fields(op.nodes());
        let bounds: Vec<_> = op.cell_nodes().iter()
            .map(|nodes| (op.physical_nodes()[nodes[0]], op.physical_nodes()[nodes[7]])).collect();
        let direction: Vec<_> = bounds.iter().map(|(lo, _)| temperature_direction(*lo)).collect();
        let prepared = op.prepare_thermal_residual_jacobian_with_motion(&u, &support_motion,
            Some(&support_motion_slope), TOTAL_STRAIN, jacobian_limits(), poll).unwrap();
        let mut jvp = vec![0.0; op.n()];
        let mut vjp = vec![0.0; op.cells()];
        op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut jvp, poll).unwrap();
        op.thermal_residual_temperature_transpose_into(&prepared, &dual, &mut vjp, poll).unwrap();
        let direct = op.thermal_residual_temperature_pullback_with_motion(&u, &dual, &support_motion,
            Some(&support_motion_slope), TOTAL_STRAIN, POINTS, poll).unwrap();
        for (actual, expected) in vjp.iter().zip(&direct) {
            close(*actual, *expected, 4e-11, 2e-11);
        }
        close(dot(&dual, &jvp), dot(&direction, &vjp), 4e-12, 2e-12);
        let h = 0.001;
        let mut residuals = Vec::new();
        for sign in [-1.0, 1.0] {
            let materials = |cell: HexCell| Ok(hot(&law,
                temperature(cell) + sign * h * temperature_direction(cell.lo())));
            op.refresh_material_field(ElasticMaterialField3::ThermoelasticByCell(&materials),
                MaterialRefreshLimits3::default(), poll).unwrap();
            let mut residual = vec![0.0; op.n()];
            op.apply(&u, &mut residual);
            let eigenstrain = op.prepare_thermoelastic(POINTS, poll).unwrap();
            for (value, force) in residual.iter_mut().zip(op.eigenstrain_load(&eigenstrain, poll).unwrap()) {
                *value -= force;
            }
            let displaced = |p, n| {
                let base = support_motion(p, n);
                let slope = support_motion_slope(p, n);
                let delta = sign * h * local_direction(&bounds, p);
                std::array::from_fn(|c| base[c] + delta * slope[c])
            };
            for (value, force) in residual.iter_mut().zip(op.prescribed_displacement_load(&displaced, poll).unwrap()) {
                *value -= force;
            }
            residuals.push(residual);
        }
        for i in 0..op.n() {
            close(jvp[i], (residuals[1][i] - residuals[0][i]) / (2.0 * h), 5e-6, 4e-8);
        }
    }
}

#[test]
fn prepared_temperature_block_refuses_stale_owners_scales_and_invalid_actions() {
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let mut op = cartesian(&law, ElasticityFormulation3::Displacement);
    let (u, dual) = fields(op.nodes());
    let prepared = op.prepare_thermal_residual_jacobian(&u, TOTAL_STRAIN, jacobian_limits(), poll).unwrap();
    let direction = vec![0.5; op.cells()];
    let sentinel = vec![19.0; op.n()];
    let mut output = sentinel.clone();
    let another = cartesian(&law, ElasticityFormulation3::Displacement);
    assert!(another.thermal_residual_temperature_apply_into(&prepared, &direction, &mut output, poll).is_err());
    assert_eq!(output, sentinel);
    op.set_scales(&vec![0.7; op.cells()]).unwrap();
    assert!(op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut output, poll).is_err());
    assert_eq!(output, sentinel);
    op.set_scales(&vec![1.0; op.cells()]).unwrap(); // Exact restoration remains usable.
    let bad = vec![f64::NAN; op.cells()];
    assert!(op.thermal_residual_temperature_apply_into(&prepared, &bad, &mut output, poll).is_err());
    assert_eq!(output, sentinel);
    assert!(op.thermal_residual_temperature_apply_into(&prepared, &direction[..direction.len()-1], &mut output, poll).is_err());
    assert_eq!(output, sentinel);
    assert!(matches!(op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut output,
        || ControlFlow::Break(())), Err(ElasticityError3::Cancelled)));
    assert_eq!(output, sentinel);
    op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut output, poll).unwrap();
    assert_ne!(output, sentinel);
    let mut transposed = vec![0.0; op.cells()];
    op.thermal_residual_temperature_transpose_into(&prepared, &dual, &mut transposed, poll).unwrap();
    close(dot(&dual, &output), dot(&direction, &transposed), 3e-12, 1e-12);
    let refreshed = hot(&law, 318.0);
    op.refresh_material_field(ElasticMaterialField3::ThermoelasticUniform(&refreshed),
        MaterialRefreshLimits3::default(), poll).unwrap();
    output.copy_from_slice(&sentinel);
    assert!(op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut output, poll).is_err());
    assert_eq!(output, sentinel);
}

#[test]
fn prepared_temperature_admission_precedes_callbacks_and_keeps_the_original_owner() {
    let material_law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let op = cartesian(&material_law, ElasticityFormulation3::Displacement);
    let (u, _) = fields(op.nodes());
    let before = cartesian_residual(&op, &u);
    let baseline = op.prepare_thermal_residual_jacobian(&u, TOTAL_STRAIN, jacobian_limits(), poll).unwrap();
    let calls = std::cell::Cell::new(0usize);
    let motion = |p, n| { calls.set(calls.get() + 1); support_motion(p, n) };
    for limits in [
        ThermalJacobianLimits3 { max_points: baseline.points() - 1, ..jacobian_limits() },
        ThermalJacobianLimits3 { max_workspace_components: baseline.workspace_components() - 1, ..jacobian_limits() },
        ThermalJacobianLimits3 { max_action_terms: baseline.action_terms() - 1, ..jacobian_limits() },
    ] {
        assert!(op.prepare_thermal_residual_jacobian_with_motion(&u, &motion, None,
            TOTAL_STRAIN, limits, poll).is_err());
        assert_eq!(calls.get(), 0, "all budgets precede boundary evaluation");
    }
    assert!(op.prepare_thermal_residual_jacobian(&u, 1e-6, jacobian_limits(), poll).is_err());
    let bad = |_, _| [f64::NAN, 0.0, 0.0];
    assert!(op.prepare_thermal_residual_jacobian_with_motion(&u, &motion, Some(&bad),
        TOTAL_STRAIN, jacobian_limits(), poll).is_err());
    calls.set(0);
    assert!(matches!(op.prepare_thermal_residual_jacobian_with_motion(&u, &motion, None,
        TOTAL_STRAIN, jacobian_limits(), || if calls.get() > 0 { ControlFlow::Break(()) } else { poll() }),
        Err(ElasticityError3::Cancelled)));
    assert_eq!(cartesian_residual(&op, &u), before);
    let kink = hot(&law([0.0, 0.0, 0.001, 0.0001, 0.0]), 300.0);
    let cut = plain(&kink, [2; 3], 0.73);
    let (u, _) = fields(cut.nodes());
    assert!(matches!(cut.prepare_thermal_residual_jacobian(&u, TOTAL_STRAIN, jacobian_limits(), poll),
        Err(ElasticityError3::Invalid(message)) if message.contains("nonsmooth")));
}

#[test]
fn explicit_isotropic_temperature_actions_retain_weak_signed_modes_without_bbar() {
    for (poisson, hydrostatic) in [(-1.0 + 1e-9, true), (0.5 - 1e-9, false)] {
        let law = law([poisson, -2.0, 0.0, 0.0, 0.0]);
        let source = hot(&law, 317.0);
        let (bulk_slope, shear_slope) = source.source().temperature_derivatives().unwrap().isotropic_moduli_pa_k.unwrap();
        let op = control(|control| CutElasticity3::build_with_material_field(domain(), [1; 3], &Top(1.2),
            ElasticMaterialField3::ThermoelasticUniform(&source), &|_| false,
            ElasticityOptions3 { formulation: ElasticityFormulation3::Displacement, ..Default::default() }, control).unwrap());
        let amplitude = 0.001;
        let dual: Vec<_> = op.nodes().iter().flat_map(|&[x, y, z]| {
            if hydrostatic { [x, y, z] } else { [x, -y, 0.0] }
        }).collect();
        let u: Vec<_> = dual.iter().map(|value| amplitude * value).collect();
        let expected = if hydrostatic { 9.0 * bulk_slope * amplitude } else { 4.0 * shear_slope * amplitude };
        let direct = op.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll).unwrap();
        close(direct[0], expected, 2e-11, 1e-14);
        let prepared = op.prepare_thermal_residual_jacobian(&u, TOTAL_STRAIN, jacobian_limits(), poll).unwrap();
        let mut force = vec![0.0; op.n()];
        op.thermal_residual_temperature_apply_into(&prepared, &[1.0], &mut force, poll).unwrap();
        close(dot(&dual, &force), expected, 2e-11, 1e-14);
    }
}

fn thermal_design_objective(
    u: &[f64],
    points: &[fs_cutfem::elastic3::thermal::ThermoelasticResponse3],
    scales: &[f64],
) -> Result<fs_cutfem::elastic3::thermal::ThermoelasticDesignObjective3, ElasticityError3> {
    let mut objective = thermal_objective(u, points, scales.len())?;
    let mut explicit = Vec::with_capacity(scales.len());
    for (cell, &scale) in scales.iter().enumerate() {
        let weight = 2e-3 * (1.0 + 0.1 * cell as f64);
        objective.value += 0.5 * weight * scale * scale;
        explicit.push(weight * scale);
    }
    Ok(
        fs_cutfem::elastic3::thermal::ThermoelasticDesignObjective3 {
            objective,
            explicit_scale_derivative: Some(explicit),
        },
    )
}

#[test]
fn thermoelastic_design_gradient_includes_thermal_load_stress_and_support_motion() {
    use fs_cutfem::elastic3::thermal::ThermoelasticMotion3;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
        let mut op = cartesian(&law, formulation);
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.55 + 0.035 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let external: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|&[x, y, z]| [0.003 * x, -0.002 * x * y, 0.001 * x * z])
            .collect();
        let motion = Some(ThermoelasticMotion3 {
            displacement: &support_motion,
            temperature_derivative: Some(&support_motion_slope),
        });
        let calls = std::cell::Cell::new(0);
        let result = op
            .solve_thermoelastic_design_adjoint(
                &external,
                motion,
                adjoint_settings(),
                |u, points| {
                    calls.set(calls.get() + 1);
                    thermal_design_objective(u, points, &scales)
                },
                |_| poll(),
            )
            .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(op.scales(), scales);
        assert_eq!(result.scale_derivative.len(), op.cells());
        // Adding design derivatives reuses exactly the same primal and adjoint.
        let temperature_only = op
            .solve_thermoelastic_adjoint(
                &external,
                motion,
                adjoint_settings(),
                |u, points| Ok(thermal_design_objective(u, points, &scales)?.objective),
                |_| poll(),
            )
            .unwrap();
        assert_eq!(
            result.thermoelastic.solution.coefficients(),
            temperature_only.solution.coefficients()
        );
        assert_eq!(
            result.thermoelastic.temperature_derivative,
            temperature_only.temperature_derivative
        );
        assert_eq!(
            result.thermoelastic.force_derivative,
            temperature_only.force_derivative
        );
        assert_eq!(
            result.thermoelastic.solution.iterations(),
            temperature_only.solution.iterations()
        );
        assert_eq!(
            result.thermoelastic.adjoint_iterations,
            temperature_only.adjoint_iterations
        );
        let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
        let h = 2e-5;
        for cell in 0..op.cells() {
            let mut values = [0.0; 2];
            for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut varied = scales.clone();
                varied[cell] += sign * h;
                op.set_scales(&varied).unwrap();
                let mut rhs = op.eigenstrain_load(&prepared, poll).unwrap();
                let lift = op
                    .prescribed_displacement_load(&support_motion, poll)
                    .unwrap();
                for ((rhs, force), lift) in rhs.iter_mut().zip(&external).zip(lift) {
                    *rhs += force + lift;
                }
                let primal = op
                    .solve_controlled(&rhs, 1e-11, 5000, 8, |_| poll())
                    .unwrap();
                let points = op
                    .thermoelastic_response(primal.coefficients(), TOTAL_STRAIN, POINTS, poll)
                    .unwrap();
                values[sample] = thermal_design_objective(primal.coefficients(), &points, &varied)
                    .unwrap()
                    .objective
                    .value;
            }
            close(
                result.scale_derivative[cell],
                (values[1] - values[0]) / (2.0 * h),
                3e-5,
                2e-9,
            );
        }
        op.set_scales(&scales).unwrap();
    }
}

#[test]
fn adaptive_thermoelastic_design_gradient_matches_joint_material_load_and_scale_variation() {
    use fs_cutfem::elastic3::thermal::ThermoelasticMotion3;
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let materials = |bounds| Ok(hot(&law, temperature(bounds)));
    let mut op = control(|control| {
        AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            &tree,
            &Top(0.73),
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, n| n[2] > 0.5,
            ElasticityOptions3 {
                ghost_gamma: 2.7,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    });
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.02 * i as f64).collect();
    let scale_direction: Vec<_> = (0..op.cells()).map(|i| 0.15 - 0.017 * i as f64).collect();
    let temperature_direction = |p: [f64; 3]| 0.4 + 0.5 * p[0] - 0.2 * p[1] + 0.3 * p[2];
    let temperature_directions: Vec<_> = op
        .cell_nodes()
        .iter()
        .map(|nodes| temperature_direction(op.physical_nodes()[nodes[0]]))
        .collect();
    let external: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| [0.003 * x, -0.002 * x * y, 0.001 * x * z])
        .collect();
    let force_direction: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| [0.002 * x, 0.001 * x * z, -0.003 * x * y])
        .collect();
    op.set_scales(&scales).unwrap();
    let result = op
        .solve_thermoelastic_design_adjoint(
            &external,
            Some(ThermoelasticMotion3 {
                displacement: &support_motion,
                temperature_derivative: None,
            }),
            adjoint_settings(),
            |u, points| thermal_design_objective(u, points, &scales),
            |_| poll(),
        )
        .unwrap();
    let derivative = dot(&result.scale_derivative, &scale_direction)
        + dot(
            &result.thermoelastic.temperature_derivative,
            &temperature_directions,
        )
        + dot(&result.thermoelastic.force_derivative, &force_direction);
    let mut values = [0.0; 2];
    let h = 2e-4;
    for (sample, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let materials = |bounds: HexCell| {
            Ok(hot(
                &law,
                temperature(bounds) + sign * h * temperature_direction(bounds.lo()),
            ))
        };
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let varied: Vec<_> = scales
            .iter()
            .zip(&scale_direction)
            .map(|(scale, direction)| scale + sign * h * direction)
            .collect();
        op.set_scales(&varied).unwrap();
        let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
        let mut rhs = op.eigenstrain_load(&prepared, poll).unwrap();
        let lift = op
            .prescribed_displacement_load(&support_motion, poll)
            .unwrap();
        for (((rhs, force), direction), lift) in rhs
            .iter_mut()
            .zip(&external)
            .zip(&force_direction)
            .zip(lift)
        {
            *rhs += force + sign * h * direction + lift;
        }
        let primal = op
            .solve_controlled(&rhs, 1e-11, 5000, 8, |_| poll())
            .unwrap();
        let points = op
            .thermoelastic_response(primal.coefficients(), TOTAL_STRAIN, POINTS, poll)
            .unwrap();
        values[sample] = thermal_design_objective(primal.coefficients(), &points, &varied)
            .unwrap()
            .objective
            .value;
    }
    close(derivative, (values[1] - values[0]) / (2.0 * h), 3e-5, 3e-9);
}

#[test]
fn thermoelastic_design_partial_admission_and_callback_cancellation_are_atomic() {
    use fs_cutfem::elastic3::thermal::ThermoelasticDesignObjective3;
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let external = vec![0.0; op.n()];
    let scales = op.scales().to_vec();
    for explicit in [vec![0.0; op.cells() + 1], vec![f64::NAN; op.cells()]] {
        assert!(matches!(
            op.solve_thermoelastic_design_adjoint(
                &external,
                None,
                adjoint_settings(),
                |u, points| Ok(ThermoelasticDesignObjective3 {
                    objective: thermal_objective(u, points, op.cells())?,
                    explicit_scale_derivative: Some(explicit),
                }),
                |_| poll(),
            ),
            Err(ElasticityError3::Invalid(_))
        ));
        assert_eq!(op.scales(), scales);
    }
    let cancel = std::cell::Cell::new(false);
    assert!(matches!(
        op.solve_thermoelastic_design_adjoint(
            &external,
            None,
            adjoint_settings(),
            |u, points| {
                cancel.set(true);
                thermal_design_objective(u, points, &scales)
            },
            |_| if cancel.get() {
                ControlFlow::Break(())
            } else {
                poll()
            },
        ),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(op.scales(), scales);
}
