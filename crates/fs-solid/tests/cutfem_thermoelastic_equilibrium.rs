//! Atomic source-temperature equilibrium with original inertia and downstream modes.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticEquilibriumSettings3,
};
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
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
use fs_solid::cutfront::dynamics3::CutElasticityRef3;
use fs_solid::cutfront::prestress3::{CutPrestressModeOptions3, CutPrestressModes3};
use fs_solver::LinearOp;

const POINTS: usize = 100_000;
fn proceed() -> ControlFlow<()> { ControlFlow::Continue(()) }
struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
enum Model {
    Cartesian(CutElasticity3),
    Adaptive(AdaptiveElasticity3),
}
macro_rules! native {
    ($owner:expr, $method:ident $(, $arg:expr)* $(,)?) => {
        match $owner {
            Model::Cartesian(op) => op.$method($($arg),*),
            Model::Adaptive(op) => op.$method($($arg),*),
        }
    };
}
impl Model {
    fn operator(&self) -> &dyn LinearOp {
        match self { Self::Cartesian(op) => op, Self::Adaptive(op) => op }
    }
    fn source(&self) -> CutElasticityRef3<'_> {
        match self {
            Self::Cartesian(op) => CutElasticityRef3::Cartesian(op),
            Self::Adaptive(op) => CutElasticityRef3::Adaptive(op),
        }
    }
}
fn build(adaptive: bool, source: &IsothermalThermoelasticState, embedded: bool) -> Model {
    let domain = HexCell::try_new([0.0; 3], [1.0; 3]).unwrap();
    let mut callback = |_| proceed();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 { depth: 0, ..Default::default() }, &mut callback,
    ).unwrap();
    let options = ElasticityOptions3 {
        formulation: ElasticityFormulation3::BBar, ..Default::default()
    };
    let field = ElasticMaterialField3::ThermoelasticUniform(source);
    let clamp = |p: [f64; 3]| p[0] == 0.0;
    if adaptive {
        let tree = Octree3::uniform(1, 4, 4096).unwrap();
        let tree = tree.refined(&[*tree.leaves().iter().next().unwrap()], proceed).unwrap();
        Model::Adaptive(if embedded {
            AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
                domain, &tree, &Top, field, &clamp, &|_, _| true, options,
                Default::default(), Default::default(), &mut control,
            ).unwrap()
        } else {
            AdaptiveElasticity3::build_with_material_field(
                domain, &tree, &Top, field, &clamp, options, &mut control,
            ).unwrap()
        })
    } else {
        let counts = if embedded { [1; 3] } else { [2, 1, 2] };
        Model::Cartesian(if embedded {
            CutElasticity3::build_with_material_field_and_embedded_dirichlet(
                domain, counts, &Top, field, &clamp, &|_, _| true, options,
                Default::default(), Default::default(), &mut control,
            ).unwrap()
        } else {
            CutElasticity3::build_with_material_field(
                domain, counts, &Top, field, &clamp, options, &mut control,
            ).unwrap()
        })
    }
}
fn settings(n: usize) -> ThermoelasticEquilibriumSettings3 {
    ThermoelasticEquilibriumSettings3 {
        refresh: Default::default(),
        maximum_native_dofs: 8 * n,
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 4000,
        poll_iterations: 32,
        maximum_free_residual_n: 2e-8,
        maximum_kinematic_strain: 0.08,
        maximum_points: POINTS,
        maximum_mass_terms: POINTS,
        maximum_prestress_entries: Some(POINTS),
        moment_origin_m: [0.2, 0.1, -0.1],
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite() &&
        (a-b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0), "{a:e} != {b:e}");
}
fn same(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) { close(a, b, tolerance); }
}
fn action(model: &Model, u: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; model.operator().n()];
    model.operator().apply(u, &mut out);
    out
}

fn law(strain_limit: f64) -> AffineIsotropicThermoelastic {
    let [poisson, young_slope, nu_slope, alpha, alpha_slope] =
        [0.2, -2.0, 0.0008, 0.0001, 0.000002];
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
                source: "authored CutFEM native thermal-state fixture".into(),
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
        &card, &elastic, &expansion, [280.0, 360.0], strain_limit,
    )
    .unwrap()
}


fn frozen(law: &AffineIsotropicThermoelastic, temperature: f64) -> IsothermalThermoelasticState {
    law.at_temperature(temperature).unwrap()
        .isothermal_elasticity(IsothermalThermoelasticInput {
            // Explicit authored reference inventory, independent of T.
            reference_density_kg_m3: 3.0,
            mechanical_material_identity: ContentHash([82; 32]),
        }).unwrap()
}

#[test]
fn changed_temperature_equilibrates_on_both_spaces_and_reaches_original_vibration() {
    let law = law(0.05);
    let cold = frozen(&law, 300.0);
    let hot = |bounds: HexCell| Ok(frozen(&law, 325.0 + 12.0 * bounds.lo()[0]
        + 8.0 * bounds.lo()[2]));
    let body = |_| [0.0, -0.02, 0.01];
    for adapted in [false, true] {
        let mut model = build(adapted, &cold, false);
        let n = model.operator().n();
        let cells = native!(&model, cells);
        let scales: Vec<_> = (0..cells).map(|i| 0.6 + 0.01 * (i % 9) as f64).collect();
        native!(&mut model, set_scales, &scales).unwrap();
        let nodes = native!(&model, nodes).as_ptr();
        let inventory = native!(&model, reference_cell_masses_kg).unwrap().to_vec();
        let old_mass = native!(&model, prepare_material_mass, POINTS, proceed).unwrap();
        let old_strain = native!(&model, prepare_thermoelastic, POINTS, proceed).unwrap();
        let occupancy = vec![0.61; cells];
        let original_mass = old_mass.diagonal(&occupancy, proceed).unwrap();
        let fixed = native!(&model, fixed).iter().position(|fixed| *fixed).unwrap();
        let mut nodal = vec![0.0; n];
        nodal[3 * fixed + 1] = 0.25;
        let free = native!(&model, fixed).iter().rposition(|fixed| !*fixed).unwrap();
        nodal[3 * free + 2] = -0.01;
        let calls = Cell::new(0);
        let query = |bounds| { calls.set(calls.get() + 1); hot(bounds) };
        if adapted {
            // Native adaptive load assembly may use the original physical
            // rows even though the equilibrium solution uses fewer masters.
            let Model::Adaptive(op) = &model else { unreachable!() };
            assert!(3 * op.physical_nodes().len() > n);
            assert!(native!(&mut model, refresh_thermoelastic_equilibrium,
                ElasticMaterialField3::ThermoelasticByCell(&query), ReferenceLoad3::body(&body),
                Some(nodal.as_slice()), None,
                ThermoelasticEquilibriumSettings3 { maximum_native_dofs: n, ..settings(n) },
                |_| proceed()).is_err());
            assert_eq!(calls.get(), 0, "physical DOF budget must precede source callbacks");
        }
        let accepted = native!(&mut model, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticByCell(&query), ReferenceLoad3::body(&body),
            Some(nodal.as_slice()), None, settings(n), |_| proceed()).unwrap();
        assert_eq!(calls.get(), cells);
        assert_eq!(native!(&model, nodes).as_ptr(), nodes);
        assert_eq!(native!(&model, reference_cell_masses_kg).unwrap(), inventory);
        assert_eq!(accepted.mass.diagonal(&occupancy, proceed).unwrap(), original_mass);
        native!(&model, validate_material_mass, &accepted.mass).unwrap();
        native!(&model, validate_eigenstrain, &accepted.eigenstrain).unwrap();
        native!(&model, validate_prestress, accepted.prestress.as_ref().unwrap()).unwrap();
        assert!(native!(&model, validate_material_mass, &old_mass).is_err());
        assert!(native!(&model, validate_eigenstrain, &old_strain).is_err());
        assert!(accepted.solution.coefficients().iter().any(|u| u.abs() > 1e-5));
        assert!(accepted.equilibrium.free_residual_norm <= settings(n).maximum_free_residual_n);
        // Full supported external loads must survive strong-row elimination.
        close(accepted.equilibrium.support_force[1], -(0.25 - 0.02 * 0.73), 3e-8);

        let mut fresh = build(adapted, &cold, false);
        native!(&mut fresh, set_scales, &scales).unwrap();
        native!(&mut fresh, refresh_material_field,
            ElasticMaterialField3::ThermoelasticByCell(&hot),
            Default::default(), proceed).unwrap();
        let free_strain = native!(&fresh, prepare_thermoelastic, POINTS, proceed).unwrap();
        let mut rhs = native!(&fresh, eigenstrain_equilibrium_load,
            &free_strain, ReferenceLoad3::body(&body), None, proceed).unwrap();
        for (f, extra) in rhs.iter_mut().zip(&nodal) { *f += extra; }
        let expected = native!(&fresh, solve_controlled, &rhs,
            settings(n).relative_residual_tolerance, 4000, 32, |_| proceed()).unwrap();
        same(accepted.solution.coefficients(), expected.coefficients(), 2e-10);
        if !adapted {
            // Consume the returned G and original-density mass directly.
            let gate = CancelGate::new_clock_free();
            let modes = CutPrestressModes3::new(model.source(),
                accepted.prestress.as_ref().unwrap(), &accepted.mass, &occupancy,
                Default::default(), &gate).unwrap();
            let report = modes.lowest_modes(CutPrestressModeOptions3 {
                modes: 1, max_ticks: 768, spectral_tolerance: 2e-7,
                residual_absolute_tolerance: 2e-7, ..Default::default()
            }).unwrap();
            assert_eq!(report.modes.len(), 1);
            assert!(report.modes[0].eigenvalue.is_finite() && report.modes[0].eigenvalue > 0.0);
        }
    }
}

#[test]
fn solve_output_budget_domain_and_final_cancellation_restore_the_accepted_owner() {
    let law = law(0.05);
    let cold = frozen(&law, 300.0);
    let hot = frozen(&law, 335.0);
    for adapted in [false, true] {
        let mut model = build(adapted, &cold, false);
        let n = model.operator().n();
        let previous = native!(&mut model, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticUniform(&cold), ReferenceLoad3::default(),
            None, None, settings(n), |_| proceed()).unwrap();
        let sample: Vec<_> = (0..n).map(|i| 0.0001 * (i % 7) as f64).collect();
        let before = action(&model, &sample);
        let inventory = native!(&model, reference_cell_masses_kg).unwrap().to_vec();
        for refused in [
            ThermoelasticEquilibriumSettings3 { maximum_iterations: 0, ..settings(n) },
            ThermoelasticEquilibriumSettings3 { maximum_mass_terms: 0, ..settings(n) },
            ThermoelasticEquilibriumSettings3 { maximum_kinematic_strain: 1e-9, ..settings(n) },
        ] {
            assert!(native!(&mut model, refresh_thermoelastic_equilibrium,
                ElasticMaterialField3::ThermoelasticUniform(&hot), ReferenceLoad3::default(),
                None, None, refused, |_| proceed()).is_err());
            assert_eq!(action(&model, &sample), before);
            assert_eq!(native!(&model, reference_cell_masses_kg).unwrap(), inventory);
            native!(&model, validate_material_mass, &previous.mass).unwrap();
            native!(&model, validate_eigenstrain, &previous.eigenstrain).unwrap();
            native!(&model, validate_prestress, previous.prestress.as_ref().unwrap()).unwrap();
        }
        // Measure a complete matching successful call; cancel the actual final
        // publication checkpoint after the next candidate payload is complete.
        let mut probe = build(adapted, &cold, false);
        let mut total = 0usize;
        native!(&mut probe, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticUniform(&hot), ReferenceLoad3::default(),
            None, None, settings(n), |_| { total += 1; proceed() }).unwrap();
        let mut calls = 0usize;
        let result = native!(&mut model, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticUniform(&hot), ReferenceLoad3::default(),
            None, None, settings(n), |_| {
                calls += 1;
                if calls == total { ControlFlow::Break(()) } else { proceed() }
            });
        assert!(matches!(result, Err(ElasticityError3::Cancelled)));
        assert_eq!(calls, total);
        assert_eq!(action(&model, &sample), before);
        native!(&model, validate_material_mass, &previous.mass).unwrap();
        native!(&model, validate_eigenstrain, &previous.eigenstrain).unwrap();
        native!(&model, validate_prestress, previous.prestress.as_ref().unwrap()).unwrap();
        let retry = native!(&mut model, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticUniform(&hot), ReferenceLoad3::default(),
            None, None, settings(n), |_| proceed()).unwrap();
        native!(&model, validate_material_mass, &retry.mass).unwrap();
        assert!(native!(&model, validate_material_mass, &previous.mass).is_err());
    }
}

#[test]
fn selected_support_strain_is_admitted_without_allocating_prestress() {
    let law = law(0.05);
    let cold = frozen(&law, 295.0);
    let target = frozen(&law, 300.0); // Integrated free strain is exactly zero.
    let mut model = build(false, &cold, true);
    let fresh = build(false, &target, true);
    let n = model.operator().n();
    let desired: Vec<_> = native!(&fresh, nodes).iter()
        .flat_map(|p| [0.01 * p[0] * p[2], 0.0, 0.0]).collect();
    let motion = |p: [f64; 3], _: [f64; 3]| [0.01 * p[0] * p[2], 0.0, 0.0];
    let lifting = native!(&fresh, prescribed_displacement_load, &motion, proceed).unwrap();
    let mut force = action(&fresh, &desired);
    for (f, b) in force.iter_mut().zip(lifting) { *f -= b; }
    let limited = ThermoelasticEquilibriumSettings3 {
        maximum_kinematic_strain: 0.0094, maximum_prestress_entries: None,
        ..settings(n)
    };
    // For u_x=.01*x*z, the retained bulk points have total norm below .0094,
    // while selected z=.73 support points exceed it. This checks the surface
    // gate independently of whether G is requested.
    native!(&fresh, admit_thermoelastic, &desired, limited.maximum_kinematic_strain,
        POINTS, proceed).unwrap();
    let original_mass = native!(&model, prepare_material_mass, POINTS, proceed).unwrap();
    let original_strain = native!(&model, prepare_thermoelastic, POINTS, proceed).unwrap();
    let before = action(&model, &desired);
    assert!(native!(&mut model, refresh_thermoelastic_equilibrium,
        ElasticMaterialField3::ThermoelasticUniform(&target), ReferenceLoad3::default(),
        Some(force.as_slice()), Some(&motion), limited, |_| proceed()).is_err());
    assert_eq!(action(&model, &desired), before);
    native!(&model, validate_material_mass, &original_mass).unwrap();
    native!(&model, validate_eigenstrain, &original_strain).unwrap();
    let accepted = native!(&mut model, refresh_thermoelastic_equilibrium,
        ElasticMaterialField3::ThermoelasticUniform(&target), ReferenceLoad3::default(),
        Some(force.as_slice()), Some(&motion),
        ThermoelasticEquilibriumSettings3 { maximum_kinematic_strain: 0.012, ..limited },
        |_| proceed()).unwrap();
    assert!(accepted.prestress.is_none());
    same(accepted.solution.coefficients(), &desired, 2e-9);
    assert!(accepted.equilibrium.free_residual_norm < limited.maximum_free_residual_n);
}

#[path = "cutfem_thermoelastic_equilibrium/dynamics.rs"]
mod dynamics;

#[path = "cutfem_thermoelastic_equilibrium/midpoint.rs"]
mod midpoint;
