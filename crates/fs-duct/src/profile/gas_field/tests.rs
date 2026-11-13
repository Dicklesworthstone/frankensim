use super::*;
use crate::compose::{ComposeOptions, compose_scattering};
use crate::{RadiusProfile, input_impedance};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::gas_field::PrimitiveGasState;
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_mms::{LadderSide, OrderGate, RefinementLadder};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    mixture::FrozenIdealGasMixtureModelV1, required_nasa9_parameter_dims,
};

const LENGTH: f64 = 0.57;
const RADIUS: f64 = 0.0073;
const TAU: f64 = core::f64::consts::TAU;

#[test]
fn nonuniform_gas_closed_chimney_keeps_resonant_limits() {
    // G3: a half-wave sealed chimney disappears at its local gas state, even
    // when the downstream bore crosses a temperature/composition gradient.
    let material = field([285.0, 335.0], [0.8, 0.4], 101325.0);
    let gas = *material.state_at(0.0).unwrap().state();
    let h = 0.04;
    let mut duct = cylinder();
    duct.segments.insert(0, Segment::ToneHole {
        hole_radius: 0.0025,
        chimney_height: h,
        bore_radius: RADIUS,
        state: crate::HoleState::Closed,
    });
    let with_hole = GasFieldDuct::new(&duct, &material, RiccatiSteps(64)).unwrap();
    let bare = GasFieldDuct::new(&cylinder(), &material, RiccatiSteps(64)).unwrap();
    let options = ComposeOptions::default();
    for phase in [core::f64::consts::PI, 0.5 * core::f64::consts::PI] {
        let omega = phase * gas.sound_speed / h;
        let chain = with_hole.chain_controlled(omega, LossModel::Lossless,
            &mut || false).unwrap();
        let (actual, _) = compose_scattering(&chain.blocks, options, &mut || false).unwrap();
        if phase == core::f64::consts::PI {
            let reference = bare.chain_controlled(omega, LossModel::Lossless,
                &mut || false).unwrap();
            let (expected, _) = compose_scattering(&reference.blocks, options,
                &mut || false).unwrap();
            for (a, e) in actual.scattering.into_iter().zip(expected.scattering) {
                assert!((a - e).abs() < 2e-12);
            }
        } else {
            assert!((actual.scattering[0] + C64::ONE).abs() < 2e-12);
            assert!(actual.scattering[1].abs() < 2e-12);
            assert!(actual.scattering[2].abs() < 2e-12);
        }
    }
}

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic constant-cp field convergence fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
// A calorically perfect NASA-9 fixture gives an independent perturbation
// oracle. These are declared analytic coefficients, not measured air data.
fn species(
    name: &str,
    mass: f64,
    cp_over_r: f64,
    viscosity: f64,
    conductivity: f64,
) -> (Nasa9StandardStateModelV1, SpeciesTransportMaterial) {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_over_r;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: vec![],
        provenance: provenance(),
    };
    let id = SpeciesId::new(name).unwrap();
    let species = Nasa9StandardStateModelV1::new(
        id.clone(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("analytic").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, viscosity),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            conductivity,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (1000.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        200.0,
                        1000.0,
                    )
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: id.to_string(),
            phase: "gas".into(),
            process: "analytic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (
        species,
        SpeciesTransportMaterial {
            species: id,
            material,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    )
}

fn gas_model() -> ThermochemicalGasModel {
    let (a, transport_a) = species("analytic-field-a", 0.028, 3.5, 1.8e-5, 0.026);
    let (b, transport_b) = species("analytic-field-b", 0.044, 4.5, 1.45e-5, 0.016);
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        )
        .unwrap(),
        vec![transport_a, transport_b],
        None,
    )
    .unwrap()
}

fn primitive(temperature_k: f64, fraction: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k,
        composition: Composition::new(
            CompositionBasis::MoleFraction,
            vec![fraction, 1.0 - fraction],
        )
        .unwrap(),
    }
}
pub(crate) fn field(temperatures: [f64; 2], fractions: [f64; 2], pressure: f64) -> AxialGasField {
    AxialGasField::new(
        vec![0.0, LENGTH],
        pressure,
        vec![
            primitive(temperatures[0], fractions[0]),
            primitive(temperatures[1], fractions[1]),
        ],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}
fn cylinder() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: RADIUS,
            length: LENGTH,
        }],
    }
}

fn profile() -> Duct {
    Duct {
        segments: vec![Segment::Profile {
            profile: std::sync::Arc::new(RadiusProfile::exponential(RADIUS, 2.0, LENGTH).unwrap()),
            steps: RiccatiSteps(64),
        }],
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}
fn matrix(field: &GasFieldDuct, omega: f64, loss: LossModel) -> [C64; 4] {
    let chain = field.chain_controlled(omega, loss, &mut || false).unwrap();
    compose_scattering(&chain.blocks, ComposeOptions::default(), &mut || false)
        .unwrap()
        .0
        .abcd()
        .unwrap()
}
fn matrix_error(a: [C64; 4], b: [C64; 4], zref: f64) -> f64 {
    let scales = [1.0, 1.0 / zref, zref, 1.0];
    let numerator = (0..4)
        .map(|i| (a[i] - b[i]).abs() * scales[i])
        .fold(0.0_f64, f64::max);
    let denominator = (0..4)
        .map(|i| b[i].abs() * scales[i])
        .fold(0.0_f64, f64::max);
    numerator / denominator
}

#[test]
fn material_wall_uniform_field_retains_original_chain_exactly() {
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let wall = WallPin { surface_density: 0.2, stiffness_per_area: 4e7, resistance: 300.0 };
    for geometry in [cylinder(), profile()] {
        let field = GasFieldDuct::new(&geometry, &source, RiccatiSteps(16)).unwrap();
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            let got = field.response_wall(TAU * 300.0, loss,
                Termination::FlangedOpen, None, Some(&wall)).unwrap();
            let exact = crate::input_impedance_wall(&geometry, field.inlet_gas(), TAU * 300.0,
                loss, Termination::FlangedOpen, Some(&wall)).unwrap();
            assert_eq!(got.impedance, exact.impedance);
            assert_eq!(got.u_mouth_over_u_in, exact.u_mouth_over_u_in);
            assert_eq!(got.p_mouth_over_u_in, exact.p_mouth_over_u_in);
        }
    }
}

#[test]
fn material_wall_and_composition_field_share_fourth_order_propagation() {
    let source = field([290.0, 390.0], [0.99, 0.45], 101325.0);
    let wall = WallPin { surface_density: 0.2, stiffness_per_area: 4e7, resistance: 300.0 };
    let geometry = profile();
    let omega = TAU * 430.0;
    let bind = |n| GasFieldDuct::new(&geometry, &source, RiccatiSteps(n)).unwrap();
    let solve = |field: &GasFieldDuct, wall| field.response_wall(
        omega, LossModel::Bessel, Termination::FlangedOpen, None, wall).unwrap();
    let fine = bind(512);
    let reference = solve(&fine, Some(&wall));
    // The walled bore's input impedance is about 17 times smaller than the rigid bore's, and
    // its relative error at each step count about 10 times larger: 6.8e-7 at 64 steps (the
    // rigid bore has 7.1e-8 there), 4.2e-8 at 128. The 2e-7 bound applies at 128 steps.
    let errors = [8, 16, 32, 64, 128].map(|n| relative(solve(&bind(n), Some(&wall)).impedance,
        reference.impedance));
    for pair in errors.windows(2) {
        let observed = det::ln(pair[0] / pair[1]) / det::ln(2.0);
        assert!((observed - 4.0).abs() < 0.2, "{errors:?}, order={observed}");
    }
    assert!(errors[4] < 2e-7, "{errors:?}");
    let rigid = solve(&fine, None);
    let stiffer = WallPin { stiffness_per_area: 2.0 * wall.stiffness_per_area, ..wall };
    let stiff = solve(&fine, Some(&stiffer));
    assert!(relative(rigid.impedance, reference.impedance) > 1e-3);
    assert!(relative(stiff.impedance, reference.impedance) > 1e-3);
    assert!(reference.impedance.re > 0.0);
    let reflection = fine.reflection_wall(omega, LossModel::Bessel,
        Termination::FlangedOpen, None, Some(&wall)).unwrap()
        .referenced_to(fine.inlet_impedance()).unwrap();
    assert!(reflection.abs() <= 1.0 + 1e-12);
    let expected = (reference.impedance - C64::from_re(fine.inlet_impedance()))
        / (reference.impedance + C64::from_re(fine.inlet_impedance()));
    assert!((reflection - expected).abs() < 1e-12);
    let different = field([290.0, 390.0], [0.99, 0.9], 101325.0);
    let alternate = GasFieldDuct::new(&geometry, &different, RiccatiSteps(64)).unwrap();
    assert!(relative(solve(&alternate, Some(&wall)).impedance, reference.impedance) > 1e-3);
    let invalid = WallPin { resistance: f64::NAN, ..wall };
    assert!(fine.response_wall(omega, LossModel::Bessel,
        Termination::FlangedOpen, None, Some(&invalid)).is_err());
    let mut polls = 0;
    assert!(matches!(fine.chain_wall_controlled(omega, LossModel::Bessel, Some(&wall),
        &mut || { polls += 1; polls > 3 }), Err(DuctError::Cancelled { completed_steps: 3 })));
}

// Independent analytic Frobenius solution of the lossless linear-T cylinder:
// (1+epsilon*s) p'' + epsilon*p' + (k0*L)^2 p = 0, s=x/L.
// a[n+2] = -(epsilon*(n+1)^2*a[n+1] + (k0*L)^2*a[n])/((n+1)*(n+2)).
// This convergent series is the same J0/Y0-in-sqrt(T) solution expressed at
// the inlet; it needs neither a marching integrator nor Bessel approximations.
// Our oracle domain |epsilon|<=0.1, k0*L<=4.3 lies inside its convergence disk.
pub(crate) fn linear_temperature_exact(gas: &GasState, delta_t: f64, omega: f64) -> [C64; 4] {
    let epsilon = delta_t / gas.temperature;
    let kl = omega * LENGTH / gas.sound_speed;
    let solution = |p: f64, derivative: f64| {
        let mut coefficients = vec![p, derivative];
        for n in 0..126 {
            let j = (n + 1) as f64;
            coefficients.push(
                -(epsilon * j * j * coefficients[n + 1] + kl * kl * coefficients[n])
                    / (j * (j + 1.0)),
            );
        }
        let p: f64 = coefficients.iter().sum();
        let derivative: f64 = coefficients
            .iter()
            .enumerate()
            .map(|(n, a)| n as f64 * a)
            .sum();
        (p, derivative)
    };
    let (f, df) = solution(1.0, 0.0);
    let (g, dg) = solution(0.0, 1.0);
    let zref = gas.density * gas.sound_speed / (core::f64::consts::PI * RADIUS * RADIUS);
    [
        C64::from_re((1.0 + epsilon) * dg),
        C64::new(0.0, -kl * g * zref),
        C64::new(0.0, (1.0 + epsilon) * df / (kl * zref)),
        C64::from_re(f),
    ]
}

#[test]
fn uniform_continuous_field_is_bitwise_original_smooth_duct() {
    let geometry = profile();
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(8)).unwrap();
    assert_eq!(bound.steps_used(), 0);
    assert_eq!(bound.field(), &source);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for termination in [
            Termination::IdealOpen,
            Termination::Closed,
            Termination::FlangedOpen,
        ] {
            let original =
                input_impedance(&geometry, bound.inlet_gas(), TAU * 610.0, loss, termination)
                    .unwrap();
            let got = bound
                .response(TAU * 610.0, loss, termination, None)
                .unwrap();
            for (a, b) in [
                (got.impedance, original.impedance),
                (got.u_mouth_over_u_in, original.u_mouth_over_u_in),
                (got.p_mouth_over_u_in, original.p_mouth_over_u_in),
            ] {
                assert_eq!(a.re.to_bits(), b.re.to_bits());
                assert_eq!(a.im.to_bits(), b.im.to_bits());
            }
            assert_eq!(got.mouth_ka.to_bits(), original.mouth_ka.to_bits());
            assert_eq!(
                got.min_shear_number.to_bits(),
                original.min_shear_number.to_bits()
            );
            assert_eq!(got.composition, original.composition);
        }
    }
}

#[test]
fn linear_temperature_field_matches_independent_analytic_solution() {
    let geometry = cylinder();
    for delta in [0.0, 10.0, 30.0] {
        let source = field([300.0, 300.0 + delta], [1.0, 1.0], 101325.0);
        let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(256)).unwrap();
        for kl in [0.4, 1.2, 2.4, 4.3] {
            let omega = kl * bound.inlet_gas().sound_speed / LENGTH;
            let expected = linear_temperature_exact(bound.inlet_gas(), delta, omega);
            let error = matrix_error(
                matrix(&bound, omega, LossModel::Lossless),
                expected,
                bound.inlet_impedance(),
            );
            assert!(error < 1e-9, "delta={delta} kL={kl}: {error:e}");
            for termination in [Termination::Closed, Termination::IdealOpen] {
                let got = bound
                    .response(omega, LossModel::Lossless, termination, None)
                    .unwrap();
                let z = if matches!(termination, Termination::Closed) {
                    expected[0] / expected[2]
                } else {
                    expected[1] / expected[3]
                };
                assert!(relative(got.impedance, z) < 1e-9);
                let reflected = bound
                    .reflection(omega, LossModel::Lossless, termination, None)
                    .unwrap()
                    .reflection;
                let reference = C64::from_re(bound.inlet_impedance());
                assert!(relative(reflected, (z - reference) / (z + reference)) < 1e-9);
            }
        }
    }
}

#[test]
fn smooth_temperature_composition_and_radius_have_fourth_order_step_ladder() {
    let geometry = profile();
    let source = field([290.0, 390.0], [0.99, 0.45], 101325.0);
    let omega = TAU * 1350.0;
    let fine = GasFieldDuct::new(&geometry, &source, RiccatiSteps(1024))
        .unwrap()
        .response(omega, LossModel::Bessel, Termination::FlangedOpen, None)
        .unwrap();
    let ns = [16, 32, 64, 128];
    let errors: Vec<_> = ns
        .iter()
        .map(|&n| {
            let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(n)).unwrap();
            assert_eq!(bound.steps_used(), n);
            let got = bound
                .response(omega, LossModel::Bessel, Termination::FlangedOpen, None)
                .unwrap();
            relative(got.impedance, fine.impedance)
        })
        .collect();
    assert!(
        errors[0] > 1e-8,
        "field/geometry discretization must be live: {errors:?}"
    );
    let ladder =
        RefinementLadder::new(ns.iter().map(|&n| 1.0 / n as f64).collect(), errors).unwrap();
    let verdict = OrderGate { theoretical: 4.0 }
        .check("continuous_gas_and_geometry", LadderSide::Primal, &ladder)
        .unwrap();
    eprintln!("{}", verdict.json_line(true));
}

#[test]
fn shared_material_pressure_and_mixture_changes_reach_physical_response() {
    let geometry = profile();
    let base = field([300.0, 340.0], [0.9, 0.6], 100000.0);
    let pressure = field([300.0, 340.0], [0.9, 0.6], 200000.0);
    let mixture = field([300.0, 340.0], [0.6, 0.1], 100000.0);
    let bind =
        |source: &AxialGasField| GasFieldDuct::new(&geometry, source, RiccatiSteps(64)).unwrap();
    let (a, b, c) = (bind(&base), bind(&pressure), bind(&mixture));
    let observe = |bound: &GasFieldDuct| {
        bound
            .response(
                TAU * 650.0,
                LossModel::Lossless,
                Termination::FlangedOpen,
                None,
            )
            .unwrap()
    };
    let (ra, rb, rc) = (observe(&a), observe(&b), observe(&c));
    assert!(relative(rb.impedance, ra.impedance.scale(2.0)) < 1e-12);
    assert!(relative(rb.p_mouth_over_u_in, ra.p_mouth_over_u_in.scale(2.0)) < 1e-12);
    assert!(relative(rb.u_mouth_over_u_in, ra.u_mouth_over_u_in) < 1e-12);
    assert_eq!(a.one_way_time_s(), b.one_way_time_s());
    assert!(relative(rc.impedance, ra.impedance) > 0.01);
    assert!((c.one_way_time_s() / a.one_way_time_s() - 1.0).abs() > 0.01);
    assert_eq!(*c.inlet_gas(), *mixture.state_at(0.0).unwrap().state());
    assert_eq!(*c.outlet_gas(), *mixture.state_at(LENGTH).unwrap().state());
    assert!(
        (rc.mouth_ka
            - TAU * 650.0 * geometry.segments[0].outlet_radius() / c.outlet_gas().sound_speed)
            .abs()
            < 1e-15
    );
}

#[test]
fn field_quadrature_reevaluates_primitives_and_splits_authored_knots() {
    let source = field([300.0, 390.0], [0.99, 0.45], 101325.0);
    let geometry = profile();
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(1)).unwrap();
    let step = &bound.quadrature[0];
    for (index, coefficient) in [0.5 - det::sqrt(3.0) / 6.0, 0.5 + det::sqrt(3.0) / 6.0]
        .into_iter()
        .enumerate()
    {
        let x = LENGTH + coefficient * (-LENGTH);
        assert_eq!(step.gases[index], *source.state_at(x).unwrap().state());
        let Segment::Profile { profile, .. } = &geometry.segments[0] else {
            unreachable!()
        };
        assert_eq!(step.radii[index], profile.radius(x).unwrap());
    }
    let station = 0.21;
    let piecewise = AxialGasField::new(
        vec![0.0, station, LENGTH],
        101325.0,
        vec![
            primitive(300.0, 0.9),
            primitive(370.0, 0.8),
            primitive(320.0, 0.4),
        ],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let spline = Duct {
        segments: vec![Segment::Profile {
            profile: std::sync::Arc::new(
                RadiusProfile::spline(vec![(0.0, RADIUS), (0.35, 0.009), (LENGTH, 0.012)], None)
                    .unwrap(),
            ),
            steps: RiccatiSteps(8),
        }],
    };
    let bound = GasFieldDuct::new(&spline, &piecewise, RiccatiSteps(4)).unwrap();
    assert_eq!(
        bound.steps_used(),
        12,
        "gas and geometry knots split three spans"
    );
}

#[test]
fn continuous_field_refuses_missing_geometry_bad_resolution_and_partial_work() {
    let geometry = cylinder();
    let uniform = field([300.0, 300.0], [1.0, 1.0], 101325.0);
    let gradient = field([300.0, 330.0], [1.0, 1.0], 101325.0);
    assert!(GasFieldDuct::new(&geometry, &uniform, RiccatiSteps(MAX_RICCATI_STEPS)).is_ok());
    for n in [0, MAX_RICCATI_STEPS + 1] {
        assert!(GasFieldDuct::new(&geometry, &gradient, RiccatiSteps(n)).is_err());
    }
    assert!(GasFieldDuct::new(&Duct { segments: vec![] }, &gradient, RiccatiSteps(4)).is_err());
    let mut short = geometry.clone();
    short.segments[0] = Segment::Cylinder {
        radius: RADIUS,
        length: LENGTH * 0.9,
    };
    assert!(GasFieldDuct::new(&short, &gradient, RiccatiSteps(4)).is_err());
    let split = AxialGasField::new(
        vec![0.0, LENGTH / 2.0, LENGTH],
        101325.0,
        vec![
            primitive(300.0, 1.0),
            primitive(320.0, 1.0),
            primitive(330.0, 1.0),
        ],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    assert!(GasFieldDuct::new(&geometry, &split, RiccatiSteps(MAX_RICCATI_STEPS / 2 + 1)).is_err());
    let mut polls = 0;
    let refused = GasFieldDuct::new_controlled(&geometry, &gradient, RiccatiSteps(4), &mut || {
        polls += 1;
        polls == 3
    });
    assert!(matches!(
        refused,
        Err(GasFieldDuctError::Duct(DuctError::Cancelled {
            completed_steps: 1
        },))
    ));
    let bound = GasFieldDuct::new(&geometry, &gradient, RiccatiSteps(4)).unwrap();
    let mut polls = 0;
    assert!(matches!(
        bound.chain_controlled(TAU * 300.0, LossModel::Bessel, &mut || {
            polls += 1;
            polls == 3
        }),
        Err(DuctError::Cancelled { completed_steps: 2 })
    ));
    let mut stale = *bound.outlet_gas();
    stale.density *= 1.1;
    assert!(
        bound
            .response(
                TAU * 300.0,
                LossModel::Bessel,
                Termination::FlangedOpen,
                Some(&stale)
            )
            .is_err()
    );
}

fn tone_hole(state: crate::HoleState) -> Segment {
    Segment::ToneHole {
        hole_radius: 0.002,
        chimney_height: 0.005,
        bore_radius: RADIUS,
        state,
    }
}

fn fingered_bore(states: [crate::HoleState; 2]) -> Duct {
    let axial = Segment::Cylinder {
        radius: RADIUS,
        length: LENGTH / 3.0,
    };
    Duct {
        segments: vec![
            axial.clone(),
            tone_hole(states[0]),
            axial.clone(),
            tone_hole(states[1]),
            axial,
        ],
    }
}

#[test]
fn uniform_fingered_field_preserves_original_bits() {
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let geometry = fingered_bore([crate::HoleState::Open, crate::HoleState::Vent(0.4)]);
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(8)).unwrap();
    assert_eq!(bound.steps_used(), 0);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for frequency in [210.0, 470.0, 900.0] {
            let expected = input_impedance(
                &geometry,
                bound.inlet_gas(),
                TAU * frequency,
                loss,
                Termination::FlangedOpen,
            )
            .unwrap();
            let actual = bound
                .response(TAU * frequency, loss, Termination::FlangedOpen, None)
                .unwrap();
            assert_eq!(actual.impedance, expected.impedance);
            assert_eq!(actual.u_mouth_over_u_in, expected.u_mouth_over_u_in);
            assert_eq!(actual.p_mouth_over_u_in, expected.p_mouth_over_u_in);
            let (observed, holes) = bound.response_with_wall_and_holes(
                TAU * frequency, loss, Termination::FlangedOpen, None, WallModel::Rigid,
            ).unwrap();
            let (original, original_holes) = crate::input_impedance_with_wall_and_holes(
                &geometry, bound.inlet_gas(), TAU * frequency, loss,
                Termination::FlangedOpen, WallModel::Rigid,
            ).unwrap();
            assert_eq!(observed.impedance, original.impedance);
            assert_eq!(holes, original_holes);
            assert_eq!(holes.len(), bound.open_hole_sites().len());
        }
    }
}

#[test]
fn warmed_entry_preserves_the_complete_hole_cluster_in_an_independent_exterior() {
    let length = LENGTH / 3.0;
    let source = AxialGasField::new(
        vec![0.0, length, LENGTH],
        101325.0,
        vec![
            primitive(330.0, 0.8),
            primitive(300.0, 0.5),
            primitive(300.0, 0.5),
        ],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let geometry = fingered_bore([crate::HoleState::Open, crate::HoleState::Vent(0.4)]);
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(32)).unwrap();
    let prefix_field = AxialGasField::new(
        vec![0.0, length],
        101325.0,
        vec![primitive(330.0, 0.8), primitive(300.0, 0.5)],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let prefix = GasFieldDuct::new(
        &Duct {
            segments: vec![geometry.segments[0].clone()],
        },
        &prefix_field,
        RiccatiSteps(32),
    )
    .unwrap();
    let tail = Duct {
        segments: geometry.segments[1..].to_vec(),
    };
    let tail_gas = source.state_at(length).unwrap();
    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 4e7,
        resistance: 300.0,
    };
    // Independent cascade: solve the entire existing uniform hole cluster as
    // one physical load, then apply only the warm prefix. This retains both
    // mutual corrections and the original lossy chimney/wall laws.
    let outside = source.state_at(0.0).unwrap();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for wall in [None, Some(&wall)] {
          for exterior in [None, Some(outside.state())] {
            let omega = TAU * 430.0;
            let (tail_response, tail_holes) = crate::exterior::input_impedance_and_holes(
                &tail,
                tail_gas.state(),
                exterior.unwrap_or(tail_gas.state()),
                omega,
                loss,
                Termination::FlangedOpen,
                WallModel::from_pin(wall),
            )
            .unwrap();
            let prefix_chain = prefix
                .chain_wall_controlled(omega, loss, wall, &mut || false)
                .unwrap();
            let expected = compose_chain(
                &prefix_chain.blocks,
                Some(tail_response.impedance),
                ComposeOptions::default(),
            )
            .unwrap();
            let actual = bound
                .response_wall(omega, loss, Termination::FlangedOpen, exterior, wall)
                .unwrap();
            assert!(relative(actual.impedance, expected.z_in) < 2e-11);
            assert!(
                relative(
                    actual.u_mouth_over_u_in,
                    expected.u_out_over_u_in * tail_response.u_mouth_over_u_in
                ) < 2e-11
            );
            assert!(
                relative(
                    actual.p_mouth_over_u_in,
                    expected.u_out_over_u_in * tail_response.p_mouth_over_u_in
                ) < 2e-11
            );
            let reflected = bound.reflection_wall(omega, loss, Termination::FlangedOpen,
                exterior, wall).unwrap().referenced_to(bound.inlet_impedance()).unwrap();
            let zref = C64::from_re(bound.inlet_impedance());
            assert!(relative(reflected, (actual.impedance-zref)/(actual.impedance+zref))<2e-11);
            let (observed, holes) = bound.response_with_wall_and_holes(
                omega, loss, Termination::FlangedOpen, exterior, WallModel::from_pin(wall),
            ).unwrap();
            assert_eq!(observed.impedance, actual.impedance);
            assert_eq!(observed.u_mouth_over_u_in, actual.u_mouth_over_u_in);
            assert_eq!(holes.len(), 2, "a field cannot silence its open chimneys");
            for (hole, tail_hole) in holes.iter().zip(&tail_holes) {
                assert!((hole.axial_position_m - length - tail_hole.axial_position_m).abs() < 1e-15);
                assert!(relative(hole.u_hole_over_u_in,
                    expected.u_out_over_u_in * tail_hole.u_hole_over_u_in) < 2e-11);
                assert!(relative(hole.u_chimney_in_over_u_in,
                    expected.u_out_over_u_in * tail_hole.u_chimney_in_over_u_in) < 2e-11);
            }
          }
        }
    }
}

#[test]
fn hole_travel_uses_the_actual_temperature_field_and_local_chimney_gas() {
    let source = field([300.0, 390.0], [0.8, 0.8], 101325.0);
    let geometry = fingered_bore([crate::HoleState::Closed, crate::HoleState::Vent(0.4)]);
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(32)).unwrap();
    let sites = bound.open_hole_sites();
    assert_eq!(sites.len(), 1);
    let site = sites[0];
    let x = 2.0 * LENGTH / 3.0;
    let local = source.state_at(x).unwrap();
    // This fixture has constant cp/composition, hence c is proportional to
    // sqrt(T). Integrating 1/c over a linear temperature field is analytic.
    let travel = 2.0 * x / (bound.inlet_gas().sound_speed + local.state().sound_speed);
    assert!((site.bore_one_way_time_s - travel).abs() < 1e-13);
    let chimney = 0.005 + fs_phs::side_hole_inner_length(0.002, RADIUS);
    assert!((site.one_way_time_s - travel - chimney / local.state().sound_speed).abs() < 1e-13);
    assert_eq!(site.effective_radius_m, 0.002 * 0.4_f64.sqrt());
    let (_, holes) = bound.response_with_wall_and_holes(
        TAU * 430.0, LossModel::Bessel, Termination::FlangedOpen, None, WallModel::Rigid,
    ).unwrap();
    assert_eq!(holes.len(), 1);
    assert_eq!(holes[0].axial_position_m, site.axial_position_m);
    assert!(holes[0].u_hole_over_u_in.abs() > 1e-5);
}

#[test]
fn single_open_hole_with_closed_fingerings_retains_fourth_order_field_convergence() {
    let source = field([300.0, 390.0], [0.95, 0.5], 101325.0);
    let geometry = fingered_bore([crate::HoleState::Closed, crate::HoleState::Open]);
    let response = |steps| {
        GasFieldDuct::new(&geometry, &source, RiccatiSteps(steps))
            .unwrap()
            .response(
                TAU * 900.0,
                LossModel::Bessel,
                Termination::FlangedOpen,
                None,
            )
            .unwrap()
            .impedance
    };
    let reference = response(256);
    let errors: Vec<_> = [2, 4, 8, 16]
        .into_iter()
        .map(|steps| relative(response(steps), reference))
        .collect();
    assert!(
        errors[0] > 1e-7,
        "field resolution must affect the fingered bore: {errors:?}"
    );
    for pair in errors.windows(2) {
        let order = det::ln(pair[0] / pair[1]) / det::ln(2.0);
        assert!((order - 4.0).abs() < 0.4, "field order={order}: {errors:?}");
    }
    let uniform = field([300.0, 300.0], [0.95, 0.95], 101325.0);
    let baseline = GasFieldDuct::new(&geometry, &uniform, RiccatiSteps(16))
        .unwrap()
        .response(
            TAU * 900.0,
            LossModel::Bessel,
            Termination::FlangedOpen,
            None,
        )
        .unwrap();
    assert!(
        relative(reference, baseline.impedance) > 1e-3,
        "the primitive field must change the physical response"
    );
}

/// Three holes (open, vented, open) between four equal cylinders.
fn fingered_bore3(states: [crate::HoleState; 3]) -> Duct {
    let axial = Segment::Cylinder {
        radius: RADIUS,
        length: LENGTH / 4.0,
    };
    Duct {
        segments: vec![
            axial.clone(),
            tone_hole(states[0]),
            axial.clone(),
            tone_hole(states[1]),
            axial.clone(),
            tone_hole(states[2]),
            axial,
        ],
    }
}

/// Two open holes 12 mm apart (about 1.6 bore radii), where the mutual term is live.
fn close_open_pair() -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder { radius: RADIUS, length: 0.3 },
            tone_hole(crate::HoleState::Open),
            Segment::Cylinder { radius: RADIUS, length: 0.012 },
            tone_hole(crate::HoleState::Open),
            Segment::Cylinder { radius: RADIUS, length: LENGTH - 0.312 },
        ],
    }
}

/// The mutual series length of every open hole, summed exactly as the binding sums it:
/// `fs_phs::side_hole_mutual_length` over all open holes in site order, from zero.
fn expected_mutual_lengths(geometry: &Duct) -> Vec<f64> {
    let mut sites = Vec::new();
    let mut position = 0.0;
    for segment in &geometry.segments {
        if let Segment::ToneHole { hole_radius, bore_radius, state, .. } = *segment {
            let sigma = crate::hole_sigma(state);
            if sigma > 0.0 {
                sites.push((position, hole_radius * sigma.sqrt(), bore_radius));
            }
        }
        position += segment.axial_length();
    }
    sites
        .iter()
        .map(|&(xi, bi, ai)| {
            let mut total = 0.0;
            for &(xj, bj, _) in &sites {
                total += fs_phs::side_hole_mutual_length(bi, bj, ai, (xj - xi).abs());
            }
            total
        })
        .collect()
}

#[test]
fn uniform_field_with_three_open_holes_is_bitwise_the_uniform_owner() {
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let geometry = fingered_bore3([
        crate::HoleState::Open,
        crate::HoleState::Vent(0.4),
        crate::HoleState::Open,
    ]);
    let bound = GasFieldDuct::new(&geometry, &source, RiccatiSteps(8)).unwrap();
    assert_eq!(bound.steps_used(), 0);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for frequency in [210.0, 470.0, 900.0] {
            let expected = input_impedance(
                &geometry,
                bound.inlet_gas(),
                TAU * frequency,
                loss,
                Termination::FlangedOpen,
            )
            .unwrap();
            let actual = bound
                .response(TAU * frequency, loss, Termination::FlangedOpen, None)
                .unwrap();
            assert_eq!(actual.impedance, expected.impedance);
            assert_eq!(actual.u_mouth_over_u_in, expected.u_mouth_over_u_in);
            assert_eq!(actual.p_mouth_over_u_in, expected.p_mouth_over_u_in);
        }
    }
}

#[test]
fn mutual_length_across_a_gradient_is_geometric_and_independent_of_the_gas() {
    // The law: t_m = -0.37 b_i b_j / a exp(-s/a), geometry only; each hole adds its sum in its own gas.
    let geometry = close_open_pair();
    let expected = expected_mutual_lengths(&geometry);
    assert_eq!(expected.len(), 2);
    assert!(
        expected.iter().all(|t| *t < -3.0e-5 && *t > -5.0e-5),
        "the term must be live (about -0.04 mm) on this pair: {expected:?}"
    );
    for gradient in [
        field([300.0, 330.0], [0.8, 0.8], 101325.0),
        field([360.0, 300.0], [0.8, 0.8], 101325.0),
        field([300.0, 330.0], [0.8, 0.3], 101325.0),
    ] {
        let bound = GasFieldDuct::new(&geometry, &gradient, RiccatiSteps(8)).unwrap();
        assert_eq!(bound.holes.len(), 2);
        for (hole, want) in bound.holes.iter().zip(&expected) {
            assert_eq!(hole.extra_series_m, *want, "t_m must not depend on the gas");
        }
        // each hole is still built in its own local gas
        for hole in &bound.holes {
            assert_eq!(hole.gas, *gradient.state_at(hole.axial_position_m).unwrap().state());
        }
    }
}

#[test]
fn a_vanishing_gradient_across_open_holes_is_continuous_with_the_uniform_owner() {
    // The nonuniform path must keep the mutual correction: dropping t_m (-0.04 mm here) would
    // move the response by about 1e-4 relative or more, far above this bound.
    let geometry = close_open_pair();
    let nearly = field([300.0, 300.0 + 1.0e-9], [0.8, 0.8], 101325.0);
    assert!(nearly.uniform_state().is_none());
    let bound = GasFieldDuct::new(&geometry, &nearly, RiccatiSteps(8)).unwrap();
    assert!(bound.steps_used() > 0);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for frequency in [210.0, 470.0, 900.0] {
            let uniform = input_impedance(
                &geometry,
                bound.inlet_gas(),
                TAU * frequency,
                loss,
                Termination::FlangedOpen,
            )
            .unwrap();
            let actual = bound
                .response(TAU * frequency, loss, Termination::FlangedOpen, None)
                .unwrap();
            assert!(
                relative(actual.impedance, uniform.impedance) < 1.0e-6,
                "{loss:?} {frequency} Hz: {}",
                relative(actual.impedance, uniform.impedance)
            );
        }
    }
}

#[test]
fn a_gradient_across_open_holes_changes_the_response() {
    let geometry = fingered_bore([crate::HoleState::Open, crate::HoleState::Vent(0.4)]);
    let gradient = field([300.0, 330.0], [0.8, 0.8], 101325.0);
    let bound = GasFieldDuct::new(&geometry, &gradient, RiccatiSteps(8)).unwrap();
    let uniform = input_impedance(
        &geometry,
        bound.inlet_gas(),
        TAU * 430.0,
        LossModel::Bessel,
        Termination::FlangedOpen,
    )
    .unwrap();
    let actual = bound
        .response(TAU * 430.0, LossModel::Bessel, Termination::FlangedOpen, None)
        .unwrap();
    assert!(
        relative(actual.impedance, uniform.impedance) > 1.0e-3,
        "the gradient must reach the physical response"
    );
}

fn watery_model() -> ThermochemicalGasModel {
    // Species A is declared the water species so that the liquid-water limit can be checked.
    let (a, transport_a) = species("analytic-field-a", 0.028, 3.5, 1.8e-5, 0.026);
    let (b, transport_b) = species("analytic-field-b", 0.044, 4.5, 1.45e-5, 0.016);
    let water = SpeciesId::new("analytic-field-a").unwrap();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        )
        .unwrap(),
        vec![transport_a, transport_b],
        Some(&water),
    )
    .unwrap()
}

/// 322 K with 0.10 water and 290 K with 0.018 water are both below the Buck saturation fractions
/// (0.115 and 0.0189 at 101325 Pa), but the linear chord at 306 K (0.059) is above 0.049.
fn chord_field(phase: ThermochemicalPhase) -> AxialGasField {
    AxialGasField::new(
        vec![0.0, LENGTH],
        101325.0,
        vec![primitive(322.0, 0.10), primitive(290.0, 0.018)],
        watery_model(),
        QueryPoint::new(),
        phase,
    )
    .unwrap()
}

#[test]
fn a_supersaturated_chord_across_open_holes_still_refuses() {
    let geometry = fingered_bore([crate::HoleState::Open, crate::HoleState::Vent(0.4)]);
    let refused = GasFieldDuct::new(
        &geometry,
        &chord_field(ThermochemicalPhase::LiquidWaterLimit),
        RiccatiSteps(8),
    )
    .unwrap_err();
    assert!(matches!(refused, GasFieldDuctError::Material(_)), "{refused}");
    assert!(refused.to_string().contains("supersaturated"), "{refused}");
    // the same stations without the phase check are admitted, so the refusal is the phase check
    assert!(
        GasFieldDuct::new(&geometry, &chord_field(ThermochemicalPhase::FrozenGas), RiccatiSteps(8)).is_ok()
    );
}

#[test]
fn fractions_that_do_not_sum_to_exact_one_still_refuse() {
    let refused = match Composition::new(CompositionBasis::MoleFraction, vec![0.8, 0.2 + 1.0e-12]) {
        Err(error) => format!("{error:?}"),
        Ok(composition) => format!(
            "{:?}",
            AxialGasField::new(
                vec![0.0, LENGTH],
                101325.0,
                vec![
                    PrimitiveGasState { temperature_k: 300.0, composition },
                    primitive(300.0, 0.8),
                ],
                gas_model(),
                QueryPoint::new(),
                ThermochemicalPhase::FrozenGas,
            )
            .unwrap_err()
        ),
    };
    assert!(refused.contains("exact") || refused.contains("sum"), "{refused}");
}

#[test]
fn varying_media_across_open_holes_are_admitted_and_closed_pads_stay_local_shunts() {
    let geometry = fingered_bore([crate::HoleState::Open, crate::HoleState::Vent(0.4)]);
    let gradient = field([300.0, 330.0], [0.8, 0.8], 101325.0);
    assert!(GasFieldDuct::new(&geometry, &gradient, RiccatiSteps(8)).is_ok());
    let hidden_gradient = AxialGasField::new(
        vec![0.0, LENGTH / 3.0, LENGTH / 2.0, 2.0 * LENGTH / 3.0, LENGTH],
        101325.0,
        vec![
            primitive(300.0, 0.8),
            primitive(300.0, 0.8),
            primitive(330.0, 0.8),
            primitive(300.0, 0.8),
            primitive(300.0, 0.8),
        ],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    // A bump between equal endpoint states is admitted too: the mutual term is geometric.
    assert!(GasFieldDuct::new(&geometry, &hidden_gradient, RiccatiSteps(8)).is_ok());
    // Closed pads have no mutual open-hole term and remain local shunts.
    let closed = fingered_bore([crate::HoleState::Closed, crate::HoleState::Closed]);
    assert!(GasFieldDuct::new(&closed, &hidden_gradient, RiccatiSteps(8)).is_ok());
    // A hole consumes a cancellable cascade block despite zero axial length.
    let bound = GasFieldDuct::new(&closed, &gradient, RiccatiSteps(1)).unwrap();
    assert_eq!(bound.steps_used(), 5);
    let mut polls = 0;
    assert!(matches!(
        bound.chain_controlled(TAU * 430.0, LossModel::Bessel, &mut || {
            polls += 1;
            polls == 2
        }),
        Err(DuctError::Cancelled { completed_steps: 1 })
    ));
}

// ---- Opt-in radiation interaction between the open holes and the bell (TMMI) ----

fn interaction(epsilon: f64) -> crate::OpenHoleInteraction {
    crate::OpenHoleInteraction::new(epsilon).unwrap()
}

/// The independent owner and the same binding with the interaction declared.
fn owners(geometry: &Duct, source: &AxialGasField, epsilon: f64) -> (GasFieldDuct, GasFieldDuct) {
    let plain = GasFieldDuct::new(geometry, source, RiccatiSteps(16)).unwrap();
    let coupled = plain.clone().with_open_hole_interaction(interaction(epsilon));
    (plain, coupled)
}

/// Frequency of the |Z| maximum inside [low, high] Hz, in cents re `reference` Hz: the best point of a
/// 1-cent grid refined by the parabola through it and its two neighbours in ln |Z|.
fn peak_cents(owner: &GasFieldDuct, low: f64, high: f64, reference: f64) -> f64 {
    let exterior = *owner.outlet_gas();
    let (lo, hi) = (1200.0 * (low / reference).log2(), 1200.0 * (high / reference).log2());
    let n = (hi - lo).round() as usize;
    let magnitude = |cents: f64| {
        owner
            .response(TAU * reference * (cents / 1200.0).exp2(), LossModel::Bessel,
                Termination::FlangedOpen, Some(&exterior))
            .unwrap()
            .impedance
            .abs()
            .ln()
    };
    let values: Vec<f64> = (0..=n).map(|i| magnitude(lo + i as f64)).collect();
    let best = (1..n).max_by(|&a, &b| values[a].total_cmp(&values[b])).unwrap();
    let (y0, y1, y2) = (values[best - 1], values[best], values[best + 1]);
    lo + best as f64 + 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2)
}

#[test]
fn a_declared_interaction_on_a_bore_with_no_open_hole_changes_nothing() {
    let source = field([300.0, 340.0], [0.8, 0.4], 101325.0);
    let geometry = fingered_bore([crate::HoleState::Closed, crate::HoleState::Closed]);
    let (plain, coupled) = owners(&geometry, &source, 2.0);
    assert_eq!(coupled.open_hole_interaction().unwrap().epsilon(), 2.0);
    assert!(plain.open_hole_interaction().is_none());
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for frequency in [210.0, 470.0, 900.0] {
            let omega = TAU * frequency;
            let a = plain.response(omega, loss, Termination::FlangedOpen, None).unwrap();
            let b = coupled.response(omega, loss, Termination::FlangedOpen, None).unwrap();
            assert_eq!(a.impedance, b.impedance);
            assert_eq!(a.u_mouth_over_u_in, b.u_mouth_over_u_in);
            let ra = plain.reflection(omega, loss, Termination::FlangedOpen, None).unwrap();
            let rb = coupled.reflection(omega, loss, Termination::FlangedOpen, None).unwrap();
            assert_eq!(ra.reflection, rb.reflection);
        }
    }
}

#[test]
fn negligible_coupling_reproduces_the_independent_response_reflection_and_hole_flows() {
    // With eps -> infinity the mutual impedance vanishes: the coupled assembly (chimney matrices, self
    // radiation, series halves, admittance two-ports, nodal solve) must close to the independent owner.
    let geometry = fingered_bore3([crate::HoleState::Open, crate::HoleState::Closed, crate::HoleState::Open]);
    for source in [
        field([300.0, 340.0], [0.8, 0.4], 101325.0),
        field([300.0, 300.0], [0.8, 0.8], 101325.0),
    ] {
        let (plain, coupled) = owners(&geometry, &source, 1.0e15);
        let exterior = *plain.outlet_gas();
        for termination in [Termination::FlangedOpen, Termination::Closed] {
            for frequency in [150.0, 310.0, 470.0, 900.0, 2400.0] {
                let omega = TAU * frequency;
                let (a, holes_a) = plain.response_with_wall_and_holes(
                    omega, LossModel::Bessel, termination, Some(&exterior), WallModel::Rigid).unwrap();
                let (b, holes_b) = coupled.response_with_wall_and_holes(
                    omega, LossModel::Bessel, termination, Some(&exterior), WallModel::Rigid).unwrap();
                assert!(relative(b.impedance, a.impedance) < 1.0e-8, "{frequency} {:?} {:?}", b.impedance, a.impedance);
                assert!((b.u_mouth_over_u_in - a.u_mouth_over_u_in).abs() < 1.0e-8 * (1.0 + a.u_mouth_over_u_in.abs()));
                assert!((b.p_mouth_over_u_in - a.p_mouth_over_u_in).abs() < 1.0e-8 * (1.0 + a.p_mouth_over_u_in.abs()));
                assert_eq!(holes_a.len(), 2);
                for (x, y) in holes_a.iter().zip(&holes_b) {
                    assert!(relative(y.u_hole_over_u_in, x.u_hole_over_u_in) < 1.0e-7, "{frequency}");
                    assert!(relative(y.u_chimney_in_over_u_in, x.u_chimney_in_over_u_in) < 1.0e-7);
                }
                let ra = plain.reflection(omega, LossModel::Bessel, termination, Some(&exterior)).unwrap();
                let rb = coupled.reflection(omega, LossModel::Bessel, termination, Some(&exterior)).unwrap();
                assert!((rb.reflection - ra.reflection).abs() < 1.0e-8, "{frequency}");
                assert_eq!(ra.reference_impedance, rb.reference_impedance);
            }
        }
    }
}

#[test]
fn the_interaction_flattens_a_close_open_pair_and_the_full_space_factor_halves_it() {
    // Two open holes 12 mm apart: the lowest resonance falls (Lefebvre, Scavone and Kergomard 2013:
    // "reduce slightly the lower resonance frequencies"); to first order the shift goes as 1 / eps.
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let (plain, half_space) = owners(&close_open_pair(), &source, 2.0);
    let (_, full_space) = owners(&close_open_pair(), &source, 4.0);
    let base = peak_cents(&plain, 150.0, 450.0, 200.0);
    let shift2 = peak_cents(&half_space, 150.0, 450.0, 200.0) - base;
    let shift4 = peak_cents(&full_space, 150.0, 450.0, 200.0) - base;
    assert!(shift2 < -0.5 && shift2 > -40.0, "{shift2}");
    assert!(shift4 < 0.0 && (shift4 / shift2 > 0.35) && (shift4 / shift2 < 0.65), "{shift2} {shift4}");
}

#[test]
fn the_interaction_refuses_vents_an_ideal_open_end_and_the_per_block_chain() {
    let source = field([300.0, 300.0], [0.8, 0.8], 101325.0);
    let omega = TAU * 300.0;
    let (_, vented) = owners(
        &fingered_bore3([crate::HoleState::Open, crate::HoleState::Vent(0.4), crate::HoleState::Open]),
        &source, 2.0);
    assert!(vented.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).is_err());
    let (plain, coupled) = owners(
        &fingered_bore3([crate::HoleState::Open, crate::HoleState::Closed, crate::HoleState::Open]),
        &source, 2.0);
    assert!(coupled.response(omega, LossModel::Bessel, Termination::IdealOpen, None).is_err());
    assert!(coupled.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).is_ok());
    // A coupled radiation matrix has no per-block chain; the independent owner keeps it.
    assert!(coupled.chain_controlled(omega, LossModel::Bessel, &mut || false).is_err());
    assert!(plain.chain_controlled(omega, LossModel::Bessel, &mut || false).is_ok());
}

/// The posted clarinet card's bore and holes (`segments` rows), with the card's hole state as open or closed.
fn card_duct(text: &str) -> Duct {
    let mut segments = Vec::new();
    for line in text.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let v = |i: usize| w[i].parse::<f64>().unwrap();
        match w.first().copied() {
            Some("cylinder") => segments.push(Segment::Cylinder { radius: v(1), length: v(2) }),
            Some("cone") => segments.push(Segment::Cone { inlet_radius: v(1), outlet_radius: v(2), length: v(3) }),
            Some("hole") => segments.push(Segment::ToneHole {
                hole_radius: v(1),
                chimney_height: v(2),
                bore_radius: v(3),
                state: if v(4) > 0.5 { crate::HoleState::Open } else { crate::HoleState::Closed },
            }),
            _ => {}
        }
    }
    Duct { segments }
}

fn clarinet_card(note: &str) -> Duct {
    let path = format!(
        "{}/../../data/instruments/descriptions/clarinet/performances/BbClar.mf.{note}.performance",
        env!("CARGO_MANIFEST_DIR"));
    card_duct(&std::fs::read_to_string(path).unwrap())
}

fn clarinet_field(duct: &Duct) -> AxialGasField {
    let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
    AxialGasField::new(
        vec![0.0, length],
        101325.0,
        vec![primitive(300.0, 0.8), primitive(300.0, 0.8)],
        gas_model(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

#[test]
fn the_posted_clarinet_cards_shift_as_the_independent_transfer_matrix_model_predicts() {
    // Passive |Z| peak with and without the TMMI coupling (eps = 2) on the posted F4 and G4 cards (12 and 14
    // open holes plus the bell), cents re the nominal pitch. An independent plane-wave transfer-matrix model
    // (Zwikker-Kosten losses, the same chimney and radiation laws, nodal TMMI solve) gives -8.8 and -9.1
    // cent on the profile card; the bands below leave room for the different test gas. D3 has no open
    // hole and must not move at all.
    let d3 = clarinet_card("D3");
    let (plain, coupled) = owners(&d3, &clarinet_field(&d3), 2.0);
    let exterior = *plain.outlet_gas();
    for frequency in [140.0, 148.0, 155.0, 450.0, 3000.0] {
        let omega = TAU * frequency;
        let a = plain.response(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior)).unwrap();
        let b = coupled.response(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior)).unwrap();
        assert_eq!(a.impedance, b.impedance);
    }
    for (note, nominal, band) in [("F4", 349.23, (-11.0, -6.5)), ("G4", 392.0, (-11.5, -6.5))] {
        let duct = clarinet_card(note);
        let (plain, coupled) = owners(&duct, &clarinet_field(&duct), 2.0);
        let (lo, hi) = (nominal * (-100.0_f64 / 1200.0).exp2(), nominal * (250.0_f64 / 1200.0).exp2());
        let (p_plain, p_coupled) = (peak_cents(&plain, lo, hi, nominal), peak_cents(&coupled, lo, hi, nominal));
        let shift = p_coupled - p_plain;
        eprintln!("TMMI {note}: passive peak {p_plain:+.3} cent, coupled {p_coupled:+.3}, shift {shift:+.3}; gas {:?}", plain.inlet_gas());
        assert!(shift > band.0 && shift < band.1, "{note}: {shift}");
    }
}

#[test]
fn the_coupled_clarinet_response_is_admitted_at_every_frequency_of_the_line_realization() {
    // The reed line is realized from the reflection on a uniform grid up to Nyquist (24 kHz at 48 kHz).
    let duct = clarinet_card("G4");
    let (_, coupled) = owners(&duct, &clarinet_field(&duct), 2.0);
    let exterior = *coupled.outlet_gas();
    let mut worst: f64 = 0.0;
    for step in 1..=480 {
        let omega = TAU * 50.0 * f64::from(step);
        let r = coupled.reflection(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior))
            .unwrap_or_else(|error| panic!("{} Hz: {error:?}", 50 * step));
        worst = worst.max(r.reflection.abs());
    }
    assert!(worst <= 1.0 + 1.0e-9, "|r| reaches {worst}");
}

// ---- Opt-in printed thick-wall hole lengths (Lefebvre and Scavone 2012) ----

/// The posted card's owner with the printed lengths, and with the TMMI coupling too when `epsilon` is given.
fn printed_owner(duct: &Duct, epsilon: Option<f64>) -> GasFieldDuct {
    let owner = GasFieldDuct::new(duct, &clarinet_field(duct), RiccatiSteps(16))
        .unwrap()
        .with_printed_hole_lengths();
    match epsilon {
        Some(epsilon) => owner.with_open_hole_interaction(interaction(epsilon)),
        None => owner,
    }
}

#[test]
fn printed_lengths_are_off_unless_declared_change_nothing_without_a_hole_and_reach_a_holed_bore() {
    let omega = TAU * 470.0;
    let exterior_of = |owner: &GasFieldDuct| *owner.outlet_gas();
    for source in [field([290.0, 390.0], [0.99, 0.45], 101325.0), field([300.0, 300.0], [0.8, 0.8], 101325.0)] {
        for geometry in [cylinder(), profile()] {
            let plain = GasFieldDuct::new(&geometry, &source, RiccatiSteps(16)).unwrap();
            let flagged = plain.clone().with_printed_hole_lengths();
            assert!(!plain.printed_hole_lengths() && flagged.printed_hole_lengths());
            let a = plain.response(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior_of(&plain))).unwrap();
            let b = flagged.response(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior_of(&plain))).unwrap();
            assert_eq!(a.impedance, b.impedance);
        }
        let holed = fingered_bore3([crate::HoleState::Open, crate::HoleState::Closed, crate::HoleState::Open]);
        let plain = GasFieldDuct::new(&holed, &source, RiccatiSteps(16)).unwrap();
        let flagged = plain.clone().with_printed_hole_lengths();
        let a = plain.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).unwrap();
        let b = flagged.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).unwrap();
        assert!(relative(b.impedance, a.impedance) > 1.0e-6, "{:?} {:?}", a.impedance, b.impedance);
    }
}

#[test]
fn the_printed_lengths_shift_the_posted_cards_as_the_independent_transfer_matrix_model_predicts() {
    // Passive |Z| peak, cents re the nominal pitch, on the posted D3 (all holes closed), F4 and G4 cards in this
    // test's gas. The independent plane-wave transfer-matrix model (same chimney and radiation laws, Zwikker-Kosten
    // losses, the printed equations coded separately in Python from Lefebvre and Scavone 2012) predicts, in this
    // exact gas: printed set alone D3 +14.17, F4 +1.84, G4 +0.97 cent; printed set with TMMI (eps 2) against
    // the independent engine F4 -6.81, G4 -7.95 cent. D3 moves through the closed-hole series length t_a(c)
    // (+15.43 cent alone) less the closed-hole matching length t_m (-1.28 cent).
    let tolerance = 0.3;
    for (note, nominal, alone, together) in
        [("D3", 146.83, 14.169, None), ("F4", 349.23, 1.838, Some(-6.813)), ("G4", 392.0, 0.965, Some(-7.952))]
    {
        let duct = clarinet_card(note);
        let plain = GasFieldDuct::new(&duct, &clarinet_field(&duct), RiccatiSteps(16)).unwrap();
        let (lo, hi) = (nominal * (-150.0_f64 / 1200.0).exp2(), nominal * (250.0_f64 / 1200.0).exp2());
        let base = peak_cents(&plain, lo, hi, nominal);
        let shift = peak_cents(&printed_owner(&duct, None), lo, hi, nominal) - base;
        eprintln!("PRINTED {note}: passive peak {base:+.3} cent, printed set shift {shift:+.3}");
        assert!((shift - alone).abs() < tolerance, "{note}: printed set shifts {shift}, model {alone}");
        if let Some(together) = together {
            let both = peak_cents(&printed_owner(&duct, Some(2.0)), lo, hi, nominal) - base;
            eprintln!("PRINTED {note}: with TMMI {both:+.3}");
            assert!((both - together).abs() < tolerance, "{note}: printed set with TMMI shifts {both}, model {together}");
        }
    }
}

#[test]
fn negligible_coupling_with_the_printed_lengths_reproduces_the_independent_printed_response() {
    // The coupled port must carry the same printed chimney the independent assembly does (inner, matching and
    // radiation lengths), and the series halves of the printed set.
    let geometry = fingered_bore3([crate::HoleState::Open, crate::HoleState::Closed, crate::HoleState::Open]);
    for source in [field([300.0, 340.0], [0.8, 0.4], 101325.0), field([300.0, 300.0], [0.8, 0.8], 101325.0)] {
        let plain = GasFieldDuct::new(&geometry, &source, RiccatiSteps(16)).unwrap().with_printed_hole_lengths();
        let coupled = plain.clone().with_open_hole_interaction(interaction(1.0e15));
        let exterior = *plain.outlet_gas();
        for frequency in [150.0, 470.0, 2400.0] {
            let omega = TAU * frequency;
            let (a, holes_a) = plain.response_with_wall_and_holes(
                omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior), WallModel::Rigid).unwrap();
            let (b, holes_b) = coupled.response_with_wall_and_holes(
                omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior), WallModel::Rigid).unwrap();
            assert!(relative(b.impedance, a.impedance) < 1.0e-8, "{frequency} {:?} {:?}", b.impedance, a.impedance);
            for (x, y) in holes_a.iter().zip(&holes_b) {
                assert!(relative(y.u_hole_over_u_in, x.u_hole_over_u_in) < 1.0e-7, "{frequency}");
            }
            let ra = plain.reflection(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior)).unwrap();
            let rb = coupled.reflection(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior)).unwrap();
            assert!((rb.reflection - ra.reflection).abs() < 1.0e-8, "{frequency}");
        }
    }
}

#[test]
fn the_printed_set_with_the_coupling_is_admitted_at_every_frequency_of_the_line_realization() {
    let duct = clarinet_card("G4");
    let owner = printed_owner(&duct, Some(2.0));
    let exterior = *owner.outlet_gas();
    let mut worst: f64 = 0.0;
    for step in 1..=480 {
        let omega = TAU * 50.0 * f64::from(step);
        let r = owner.reflection(omega, LossModel::Bessel, Termination::FlangedOpen, Some(&exterior))
            .unwrap_or_else(|error| panic!("{} Hz: {error:?}", 50 * step));
        worst = worst.max(r.reflection.abs());
    }
    assert!(worst <= 1.0 + 1.0e-9, "|r| reaches {worst}");
}
