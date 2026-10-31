//! Independent nonlinear slab/series oracles at picofarad and picocoulomb scales.
//! Authored constitutive fixtures; these are not measured material predictions.
use fs_em::electrostatic::{
    ElectrostaticControls, PlanarElectrostatics,
    nonlinear::{
        DielectricMaterial, ElectrodeControl, NonlinearElectrostaticControls,
        PlanarNonlinearElectrostatics, tangent::DifferentialCapacitanceControls,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        dielectric::reversible::{
            ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS,
            REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{Dims, QuantitySpec};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

const GAP: f64 = 1e-3;
const WIDTH: f64 = 2e-2;
const THICKNESS: f64 = 5e-3;
const AREA: f64 = WIDTH * THICKNESS;
const E0: f64 = 1e5;
const D0: f64 = 1e-6;
const EPSILON0: f64 = D0 / E0;
const V0: f64 = E0 * GAP;
const Q0: f64 = AREA * D0;
const U0: f64 = AREA * GAP * E0 * D0;

fn near(actual: f64, expected: f64, physical_scale: f64) {
    let tolerance = 5e-10 * physical_scale.max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance={tolerance:e}"
    );
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 32,
        max_dofs: 64,
        max_triplets: 288,
    }
}
fn linear_controls() -> CheckedCgConfig {
    CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 128,
        max_corrections: 2,
    }
}
fn controls() -> NonlinearElectrostaticControls {
    NonlinearElectrostaticControls {
        absolute_residual_tolerance_c: 1e-22,
        relative_residual_tolerance: 1e-12,
        charge_tolerance_c: 1e-20,
        energy_tolerance_j: 1e-16,
        max_newton_iterations: 32,
        max_backtracks: 24,
        max_operator_applications: 2048,
        linear: linear_controls(),
    }
}
fn differential_controls() -> DifferentialCapacitanceControls {
    DifferentialCapacitanceControls {
        maximum_terminals: 8,
        maximum_matrix_entries: 64,
        max_operator_applications: 1024,
        coefficient_tolerance_f: 1e-22,
        linear: linear_controls(),
    }
}
fn card(
    name: &str,
    dims: Dims,
    value: PropertyValue,
    validity: ValidityDomain,
    interpolation: InterpolationPolicy,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored nonlinear electrostatic fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: name.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reversible(knots: &[(f64, f64)]) -> ReversibleDielectricMaterial {
    let source = card(
        REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ELECTRIC_DISPLACEMENT_DIMS,
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: knots.iter().map(|&(e, d)| (e * E0, d * D0)).collect(),
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        ValidityDomain::unconstrained().with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
            0.0,
            4.0 * E0,
        ),
        InterpolationPolicy::LinearInside,
    );
    let point = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    ReversibleDielectricMaterial::resolve(
        &source,
        &point,
        "E",
        4.0 * E0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn nonlinear() -> ReversibleDielectricMaterial {
    reversible(&[(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)])
}
fn scalar(name: &str, dims: Dims, value: f64) -> ResolvedMaterialStatePoint {
    let source = card(
        name,
        dims,
        PropertyValue::Scalar { value, dims },
        ValidityDomain::unconstrained(),
        InterpolationPolicy::ConstantWithinValidity,
    );
    resolve_material_state_point(
        &source,
        &QueryPoint::new(),
        &[
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn linear(epsilon: f64) -> ResolvedMaterialStatePoint {
    scalar(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        epsilon,
    )
}
fn slab_mesh() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "nonlinear-capacitor-slab",
        vec![
            [0.0, 0.0],
            [GAP / 2.0, 0.0],
            [GAP, 0.0],
            [0.0, WIDTH],
            [GAP / 2.0, WIDTH],
            [GAP, WIDTH],
        ],
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::planar(THICKNESS).unwrap(),
    )
    .unwrap()
}
fn split_mesh() -> TriComplex2 {
    // Distinct dielectric components meet the same conductor. The duplicated
    // interface vertices are joined electrically only by its terminal trace.
    TriComplex2::from_indexed_triangles(
        "two-dielectrics-one-floating-conductor",
        vec![
            [0.0, 0.0],
            [GAP, 0.0],
            [GAP, WIDTH],
            [0.0, WIDTH],
            [GAP, 0.0],
            [2.0 * GAP, 0.0],
            [2.0 * GAP, WIDTH],
            [GAP, WIDTH],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(THICKNESS).unwrap(),
    )
    .unwrap()
}
fn series_problem() -> PlanarNonlinearElectrostatics {
    let nonlinear = DielectricMaterial::Reversible(nonlinear());
    let linear = DielectricMaterial::isotropic(linear(EPSILON0));
    PlanarNonlinearElectrostatics::new(
        split_mesh(),
        vec![nonlinear.clone(), nonlinear, linear.clone(), linear],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn nonlinear_picofarad_slab_has_integral_storage_and_distinct_differential_capacitance() {
    let law = nonlinear();
    let identity = law.material().card_identity();
    let problem = PlanarNonlinearElectrostatics::new(
        slab_mesh(),
        vec![DielectricMaterial::Reversible(law); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let solve = |voltage| {
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(voltage),
                ],
                &[0.0, voltage / 2.0, voltage, 0.0, voltage / 2.0, voltage],
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    let field = solve(2.0 * V0);
    // At E/E0=2, D/D0=2.5. Independently integrating the three explicit
    // affine segments gives W/(E0 D0)=1.75 and W*/(E0 D0)=3.25.
    near(field.terminal_charge_c()[0], -2.5 * Q0, Q0);
    near(field.terminal_charge_c()[1], 2.5 * Q0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    assert!((field.energy_j() - 0.5 * 2.0 * V0 * 2.5 * Q0).abs() > 0.5 * U0);
    for (electric, displacement) in field
        .electric_field_v_m()
        .iter()
        .zip(field.displacement_c_m2())
    {
        near(electric[0], -2.0 * E0, E0);
        near(electric[1], 0.0, E0);
        near(displacement[0], -2.5 * D0, D0);
        near(displacement[1], 0.0, D0);
    }
    for material in field.problem().materials() {
        assert_eq!(material.material().card_identity(), identity);
    }
    let low = solve(2.0 * V0 - 1.0);
    let high = solve(2.0 * V0 + 1.0);
    let differential = (high.terminal_charge_c()[1] - low.terminal_charge_c()[1]) / 2.0;
    let secant = field.terminal_charge_c()[1] / (2.0 * V0);
    near(differential, 0.5 * EPSILON0 * AREA / GAP, 1e-12);
    near(secant, 1.25 * EPSILON0 * AREA / GAP, 1e-12);
    assert!(secant > 2.0 * differential);
    assert!(field.energy_identity_residual_j().abs() <= controls().energy_tolerance_j);
    assert!(field.maximum_free_charge_residual_c() <= controls().charge_tolerance_c);
    // An accepted-state differential matrix follows the actual branch slope,
    // independently of both secant capacitance and integrated stored energy.
    for (voltage, slope) in [(0.5 * V0, 2.0), (2.0 * V0, 0.5), (3.5 * V0, 4.0)] {
        let biased = solve(voltage);
        let tangent = biased
            .differential_capacitance(differential_controls(), &mut || false)
            .unwrap();
        assert!(std::ptr::eq(tangent.field(), &biased));
        assert_eq!(tangent.terminals(), &[0, 1]);
        let c = slope * EPSILON0 * AREA / GAP;
        for (actual, expected) in tangent.matrix_f().iter().zip([c, -c, -c, c]) {
            near(*actual, expected, 1e-12);
        }
        assert!(tangent.maximum_residual_f() <= differential_controls().coefficient_tolerance_f);
    }
}

#[test]
fn charged_common_electrode_solves_the_actual_nonlinear_series_equation() {
    let problem = series_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2, 4, 7], &[5, 6]], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(0.5 * Q0),
                ElectrodeControl::Voltage(4.0 * V0),
            ],
            &[0.0; 8],
            controls(),
            &mut || false,
        )
        .unwrap();
    // D_left/D0=1.5+0.5*v; D_right/D0=4-v; their difference is 0.5.
    // Thus the independent root is v=2, including the prescribed conductor charge.
    for (actual, expected) in field
        .terminal_voltage_v()
        .iter()
        .zip([0.0, 2.0 * V0, 4.0 * V0])
    {
        near(*actual, expected, V0);
    }
    for (actual, expected) in field
        .terminal_charge_c()
        .iter()
        .zip([-2.5 * Q0, 0.5 * Q0, 2.0 * Q0])
    {
        near(*actual, expected, Q0);
    }
    near(field.energy_j(), 3.75 * U0, U0);
    near(field.coenergy_j(), 5.25 * U0, U0);
    for node in [1, 2, 4, 7] {
        near(field.relative_voltage_v()[node], 2.0 * V0, V0);
    }
    assert!(field.report().newton_iterations > 0);
    assert!(field.maximum_floating_charge_residual_c() <= controls().charge_tolerance_c);
    assert!(field.charge_balance_residual_c().abs() <= controls().charge_tolerance_c);
    let tangent = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(tangent.terminals(), &[0, 2]);
    // The middle charge is fixed while its voltage responds. Independent
    // series condensation gives (0.5 pF * 1 pF)/(0.5 pF + 1 pF) = 1/3 pF.
    let series_c = EPSILON0 * AREA / GAP / 3.0;
    for (actual, expected) in tangent
        .matrix_f()
        .iter()
        .zip([series_c, -series_c, -series_c, series_c])
    {
        near(*actual, expected, 1e-12);
    }
    let perturbed = |increment| {
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(0.5 * Q0),
                    ElectrodeControl::Voltage(4.0 * V0 + increment),
                ],
                field.relative_voltage_v(),
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    let increment = 0.25;
    let low = perturbed(-increment);
    let high = perturbed(increment);
    for (terminal, expected) in [(0, -series_c), (1, 0.0), (2, series_c)] {
        near(
            (high.terminal_charge_c()[terminal] - low.terminal_charge_c()[terminal])
                / (2.0 * increment),
            expected,
            1e-12,
        );
    }
    near(
        (high.terminal_voltage_v()[1] - low.terminal_voltage_v()[1]) / (2.0 * increment),
        2.0 / 3.0,
        1.0,
    );
}

#[test]
fn differential_extraction_refuses_genuine_kinks_and_preserves_accepted_fields_on_abort() {
    // Binary geometry places the represented field exactly on its source
    // knot; decimal element inversion can put it one ULP to either side.
    let gap = 1.0 / 1024.0;
    let width = 1.0 / 64.0;
    let thickness = 1.0 / 256.0;
    let voltage_scale = E0 * gap;
    let energy_scale = width * thickness * gap * E0 * D0;
    let mesh = || {
        TriComplex2::from_indexed_triangles(
            "exact-knot-dielectric",
            vec![
                [0.0, 0.0],
                [gap / 2.0, 0.0],
                [gap, 0.0],
                [0.0, width],
                [gap / 2.0, width],
                [gap, width],
            ],
            vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
            Metric2::planar(thickness).unwrap(),
        )
        .unwrap()
    };
    let problem = PlanarNonlinearElectrostatics::new(
        mesh(),
        vec![DielectricMaterial::Reversible(nonlinear()); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let solve = |voltage| {
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(voltage),
                ],
                &[0.0, voltage / 2.0, voltage, 0.0, voltage / 2.0, voltage],
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    for voltage in [voltage_scale, 3.0 * voltage_scale] {
        let kink = solve(voltage);
        for electric in kink.electric_field_v_m() {
            assert_eq!(electric[0].hypot(electric[1]), voltage / gap);
        }
        let error = kink
            .differential_capacitance(differential_controls(), &mut || false)
            .err()
            .unwrap_or_else(|| {
                panic!(
                    "distinct one-sided slopes need refusal: voltage={voltage:e}, fields={:?}",
                    kink.electric_field_v_m()
                )
            });
        assert!(error.to_string().contains("not unique"), "{error}");
    }
    let field = solve(2.0 * voltage_scale);
    let mut polls = 0;
    let admitted = field
        .differential_capacitance(differential_controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(
            field
                .differential_capacitance(differential_controls(), &mut || {
                    seen += 1;
                    seen >= stop
                })
                .is_err()
        );
        assert_eq!(seen, stop);
    }
    let mut bounded = differential_controls();
    bounded.maximum_matrix_entries = 3;
    assert!(
        field
            .differential_capacitance(bounded, &mut || false)
            .is_err()
    );
    bounded = differential_controls();
    bounded.max_operator_applications = 0;
    assert!(
        field
            .differential_capacitance(bounded, &mut || false)
            .is_err()
    );
    let retried = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(retried.matrix_f(), admitted.matrix_f());
    near(field.energy_j(), 1.75 * energy_scale, energy_scale);

    // A source breakpoint with equal slopes is still differentiable.
    let smooth = PlanarNonlinearElectrostatics::new(
        mesh(),
        vec![
            DielectricMaterial::Reversible(reversible(&[
                (0.0, 0.0),
                (1.0, 2.0),
                (2.0, 4.0),
                (4.0, 8.0),
            ]));
            4
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let smooth_electrodes = smooth
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let smooth_field = smooth_electrodes
        .solve(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(voltage_scale),
            ],
            &[
                0.0,
                voltage_scale / 2.0,
                voltage_scale,
                0.0,
                voltage_scale / 2.0,
                voltage_scale,
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
    let smooth_tangent = smooth_field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    near(
        smooth_tangent.matrix_f()[3],
        2.0 * EPSILON0 * width * thickness / gap,
        1e-12,
    );
}

#[test]
fn one_voltage_reference_anchors_components_joined_by_a_shared_floating_conductor() {
    let problem = series_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2, 4, 7], &[5, 6]], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(0.5 * Q0),
                ElectrodeControl::Charge(2.0 * Q0),
            ],
            &[0.0; 8],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.terminal_voltage_v()[1], 2.0 * V0, V0);
    near(field.terminal_voltage_v()[2], 4.0 * V0, V0);
    near(field.energy_j(), 3.75 * U0, U0);
    // Splitting that common conductor leaves the right component without a
    // voltage reference. Compatible total charges alone cannot fix its gauge.
    let separated = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2], &[4, 7], &[5, 6]], &mut || false)
        .unwrap();
    assert!(
        separated
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(2.5 * Q0),
                    ElectrodeControl::Charge(-2.0 * Q0),
                    ElectrodeControl::Charge(2.0 * Q0)
                ],
                &[0.0; 8],
                controls(),
                &mut || false,
            )
            .is_err()
    );
    let all_floating = [
        ElectrodeControl::Charge(-2.5 * Q0),
        ElectrodeControl::Charge(0.5 * Q0),
        ElectrodeControl::Charge(2.0 * Q0),
    ];
    assert!(
        electrodes
            .solve(&[0.0; 4], &all_floating, &[0.0; 8], controls(), &mut || {
                false
            })
            .is_err()
    );
}

#[test]
fn original_linear_source_and_reversible_linear_limit_match_with_volume_charge() {
    let material = linear(2.0 * EPSILON0);
    let baseline = PlanarElectrostatics::new(
        slab_mesh(),
        vec![material.clone(); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let rho = [1e-3; 4];
    let original = baseline
        .solve(
            &rho,
            &[(0, 0.0), (3, 0.0), (2, 2.0 * V0), (5, 2.0 * V0)],
            ElectrostaticControls {
                linear: linear_controls(),
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    for model in [
        DielectricMaterial::isotropic(material.clone()),
        DielectricMaterial::Reversible(reversible(&[(0.0, 0.0), (4.0, 8.0)])),
    ] {
        let identity = model.material().card_identity();
        let problem =
            PlanarNonlinearElectrostatics::new(slab_mesh(), vec![model; 4], limits(), &mut || {
                false
            })
            .unwrap();
        let electrodes = problem
            .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
            .unwrap();
        let field = electrodes
            .solve(
                &rho,
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(2.0 * V0),
                ],
                &[0.0, V0, 2.0 * V0, 0.0, V0, 2.0 * V0],
                controls(),
                &mut || false,
            )
            .unwrap();
        for (actual, expected) in field
            .electric_field_v_m()
            .iter()
            .zip(original.electric_field_v_m())
        {
            for axis in 0..2 {
                near(actual[axis], expected[axis], E0);
            }
        }
        near(field.energy_j(), original.energy_j(), U0);
        near(field.coenergy_j(), field.energy_j(), U0);
        near(field.volume_charge_c(), 1e-3 * AREA * GAP, Q0);
        near(
            field.terminal_charge_c().iter().sum::<f64>() + field.volume_charge_c(),
            0.0,
            Q0,
        );
        for source in field.problem().materials() {
            assert_eq!(source.material().card_identity(), identity);
        }
    }
    // A conductivity coefficient carries no dielectric storage authority.
    let conductor = scalar("electrical_conductivity", ELECTRICAL_CONDUCTIVITY_DIMS, 2.0);
    assert!(
        PlanarNonlinearElectrostatics::new(
            slab_mesh(),
            vec![DielectricMaterial::isotropic(conductor); 4],
            limits(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn common_bias_and_bounded_refusals_preserve_the_same_prepared_physical_problem() {
    let problem = series_problem();
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[1, 2, 4, 7], &[5, 6]], &mut || false)
        .unwrap();
    let electrical = [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Charge(0.5 * Q0),
        ElectrodeControl::Voltage(4.0 * V0),
    ];
    let mut polls = 0;
    let accepted = electrodes
        .solve(&[0.0; 4], &electrical, &[0.0; 8], controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let bias = 1e12;
    let shifted = electrodes
        .solve(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(bias),
                ElectrodeControl::Charge(0.5 * Q0),
                ElectrodeControl::Voltage(bias + 4.0 * V0),
            ],
            &[0.0; 8],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(shifted.reference_terminal(), 0);
    assert_eq!(shifted.reference_voltage_v(), bias);
    assert_eq!(shifted.relative_voltage_v(), accepted.relative_voltage_v());
    assert_eq!(shifted.terminal_charge_c(), accepted.terminal_charge_c());
    assert_eq!(shifted.energy_j(), accepted.energy_j());
    let mut no_updates = controls();
    no_updates.max_newton_iterations = 0;
    assert!(
        electrodes
            .solve(&[0.0; 4], &electrical, &[0.0; 8], no_updates, &mut || false)
            .is_err()
    );
    let mut no_operator = controls();
    no_operator.max_operator_applications = 0;
    assert!(
        electrodes
            .solve(&[0.0; 4], &electrical, &[0.0; 8], no_operator, &mut || {
                false
            })
            .is_err()
    );
    assert!(
        electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(0.5 * Q0),
                    ElectrodeControl::Voltage(20.0 * V0)
                ],
                &[0.0; 8],
                controls(),
                &mut || false,
            )
            .is_err()
    );
    let mut inconsistent = [0.0; 8];
    inconsistent[1] = V0;
    assert!(
        electrodes
            .solve(
                &[0.0; 4],
                &electrical,
                &inconsistent,
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(polls > 3);
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(
            electrodes
                .solve(&[0.0; 4], &electrical, &[0.0; 8], controls(), &mut || {
                    seen += 1;
                    seen >= stop
                })
                .is_err()
        );
        assert_eq!(seen, stop);
    }
    let retried = electrodes
        .solve(&[0.0; 4], &electrical, &[0.0; 8], controls(), &mut || false)
        .unwrap();
    assert_eq!(retried.relative_voltage_v(), accepted.relative_voltage_v());
    assert_eq!(retried.terminal_charge_c(), accepted.terminal_charge_c());
    assert_eq!(retried.energy_j(), accepted.energy_j());
}
