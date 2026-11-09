//! G0/G1/G3: temperature transposes of refreshed native CutFEM mechanics.
//! Authored source curves verify numerical mechanics, not material calibration.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
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
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{
    AffineIsotropicThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::op::LinearOp;

const POINTS: usize = 100_000;
const TOTAL_STRAIN: f64 = 0.1;

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
        QuadratureOptions3 { depth: 0, ..Default::default() },
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
        claims.insert_claim(PropertyClaim {
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
                .with_quantity("T", temperature_quantity, 280.0, 360.0),
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
                YOUNG_MODULUS_PROPERTY, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY, Dims([0; 6]), ScalarAdmissibility::Finite,
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
    AffineIsotropicThermoelastic::new(
        &card, &elastic, &expansion, [280.0, 360.0], 0.05,
    )
    .unwrap()
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

fn cartesian(law: &AffineIsotropicThermoelastic, formulation: ElasticityFormulation3) -> CutElasticity3 {
    let materials = |cell| Ok(hot(law, temperature(cell)));
    control(|control| {
        CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            [2; 3],
            &Top(0.73),
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 { formulation, ghost_gamma: 2.7, ..Default::default() },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    })
}

fn fields(nodes: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>) {
    let u = nodes.iter().flat_map(|&[x, y, z]| [
        0.001 * x * x * (1.0 + y),
        -0.0008 * x * (z * z + 0.2 * y),
        0.0006 * x * y * y,
    ]).collect();
    // Algebraic dual strains deliberately exceed the source's physical limit.
    let dual = nodes.iter().flat_map(|&[x, y, z]| [
        1.7 * x * (1.0 + 0.3 * x + z),
        -0.4 * (1.0 + x * y),
        0.9 * x * z * z,
    ]).collect();
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

#[test]
fn full_temperature_pullback_matches_each_refreshed_cell_with_nitsche_and_ghosts() {
    for (formulation, parameters) in [
        (ElasticityFormulation3::Displacement, [0.2, -2.0, 0.0008, 0.0001, 0.000002]),
        (ElasticityFormulation3::BBar, [0.49, -2.0, 0.00002, 0.0001, 0.000002]),
    ] {
        let law = law(parameters);
        let at = hot(&law, 317.0);
        let c = at.acoustic().stiffness_mandel_pa();
        let dc = at.source().temperature_derivatives().unwrap().stiffness_mandel_pa_k;
        // Changing Poisson response cannot be replaced by a scalar dE/E.
        assert!((dc[0][0] * c[3][3] - dc[3][3] * c[0][0]).abs() > 1.0);
        let mut op = cartesian(&law, formulation);
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.04 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let (u, dual) = fields(op.nodes());
        let gradient = op.thermal_residual_temperature_pullback(
            &u, &dual, TOTAL_STRAIN, POINTS, poll,
        ).unwrap();
        let lower: Vec<_> = op.cell_nodes().iter().map(|nodes| op.nodes()[nodes[0]]).collect();
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
                ).unwrap();
                values[sample] = dot(&dual, &cartesian_residual(&op, &u));
                assert_eq!(op.scales(), scales.as_slice());
                assert_eq!(op.reference_cell_masses_kg().unwrap(), original_mass.as_slice());
            }
            close(gradient[cell], (values[1] - values[0]) / (2.0 * h), 3e-6, 2e-8);
        }
    }
}

#[test]
fn adaptive_temperature_pullback_uses_the_original_hanging_map_on_both_fields() {
    let law = law([0.2, -2.0, 0.0008, 0.0001, 0.000002]);
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree.refined(&[*tree.leaves().iter().next().unwrap()], poll).unwrap();
    let materials = |cell| Ok(hot(&law, temperature(cell)));
    let mut op = control(|control| {
        AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            &tree,
            &Top(0.73),
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 { ghost_gamma: 2.7, ..Default::default() },
            Default::default(),
            Default::default(),
            control,
        ).unwrap()
    });
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.02 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let (u, dual) = fields(op.nodes());
    let gradient = op.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS, poll,
    ).unwrap();
    let lower: Vec<_> = op.cell_nodes().iter()
        .map(|nodes| op.physical_nodes()[nodes[0]]).collect();
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
            ).unwrap();
            let mut residual = vec![0.0; op.n()];
            op.apply(&u, &mut residual);
            let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
            for (value, force) in residual.iter_mut()
                .zip(op.eigenstrain_load(&source, poll).unwrap())
            {
                *value -= force;
            }
            values[sample] = dot(&dual, &residual);
        }
        close(gradient[cell], (values[1] - values[0]) / (2.0 * h), 4e-6, 3e-8);
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
        ).unwrap()
    })
}

#[test]
fn nonsmooth_native_stabilization_refuses_only_when_used() {
    let kink = hot(&law([0.0, 0.0, 0.001, 0.0001, 0.0]), 300.0);
    let full = plain(&kink, [1; 3], 1.2); // No Nitsche or ghost modulus.
    let (u, dual) = fields(full.nodes());
    full.thermal_residual_temperature_pullback(&u, &dual, TOTAL_STRAIN, POINTS, poll).unwrap();

    let cut = plain(&kink, [2; 3], 0.73);
    let (u, dual) = fields(cut.nodes());
    let error = cut.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS, poll,
    ).unwrap_err();
    assert!(matches!(error, ElasticityError3::Invalid(message) if message.contains("nonsmooth")));

    // Zero C' and equal active row slopes are smooth, despite zero couplings.
    let smooth = hot(&law([0.0, 0.0, 0.0, 0.0001, 0.000002]), 300.0);
    let cut = plain(&smooth, [2; 3], 0.73);
    let (u, dual) = fields(cut.nodes());
    let derivative = cut.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS, poll,
    ).unwrap();
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
    let expected = op.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS,
        || { calls += 1; poll() },
    ).unwrap();
    assert!(calls > 10);
    for stop in [1, calls / 2, calls] {
        let mut at = 0usize;
        assert!(matches!(
            op.thermal_residual_temperature_pullback(
                &u, &dual, TOTAL_STRAIN, POINTS,
                || {
                    at += 1;
                    if at == stop { ControlFlow::Break(()) } else { poll() }
                },
            ),
            Err(ElasticityError3::Cancelled),
        ));
    }
    assert!(op.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, 0, poll,
    ).is_err());
    assert!(op.thermal_residual_temperature_pullback(
        &u, &dual, 1e-6, POINTS, poll,
    ).is_err());
    let invalid: Vec<_> = u.iter().map(|value| value * 100.0).collect();
    assert!(op.thermal_residual_temperature_pullback(
        &invalid, &dual, 0.9, POINTS, poll,
    ).is_err()); // Native total/elastic source bounds remain stricter.
    let large_dual: Vec<_> = dual.iter().map(|value| value * 1e4).collect();
    let actual = op.thermal_residual_temperature_pullback(
        &u, &large_dual, TOTAL_STRAIN, POINTS, poll,
    ).unwrap();
    for (actual, expected) in actual.iter().zip(&expected) {
        close(*actual, 1e4 * expected, 3e-12, 1e-8);
    }
    op.validate_eigenstrain(&source).unwrap();
    assert_eq!(cartesian_residual(&op, &u), before);
    assert_eq!(op.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS, poll,
    ).unwrap(), expected);

    let state = hot(&law, 316.0);
    let cold = control(|control| {
        CutElasticity3::build_with_material_field(
            domain(), [1; 3], &Top(1.2),
            ElasticMaterialField3::Uniform(state.acoustic()),
            &|p| p[0] == 0.0, Default::default(), control,
        ).unwrap()
    });
    let (u, dual) = fields(cold.nodes());
    assert!(cold.thermal_residual_temperature_pullback(
        &u, &dual, TOTAL_STRAIN, POINTS, poll,
    ).is_err());
}
