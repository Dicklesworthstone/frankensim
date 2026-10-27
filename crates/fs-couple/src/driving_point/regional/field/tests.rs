use super::*;
use fs_duct::{LossModel, Termination};
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
fn continuous_gas_fields_keep_open_holes_in_physical_microphone_output() {
    use crate::pcm_wav::baffled::CircularOutletReceiver;
    use fs_duct::{HoleState, RiccatiSteps, WallModel};

    let radius = 0.0015;
    let geometry = Duct { segments: vec![
        Segment::Cylinder { radius, length: LENGTH / 2.0 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.003,
            bore_radius: radius, state: HoleState::Open },
        Segment::Cylinder { radius, length: LENGTH / 2.0 },
    ] };
    let source = field(360.0);
    let outside = *source.state_at(0.0).unwrap().state();
    let duct = RegionalDuct::from_continuous_gas_field(&geometry, &source, RiccatiSteps(16))
        .unwrap().with_exterior_gas(outside).unwrap();
    let field_owner = duct.continuous_gas_field().unwrap();
    let (response, holes) = duct.response_with_wall_and_holes(
        TAU * 430.0, LossModel::Bessel, Termination::FlangedOpen, WallModel::Rigid,
    ).unwrap();
    let (direct, direct_holes) = field_owner.response_with_wall_and_holes(
        TAU * 430.0, LossModel::Bessel, Termination::FlangedOpen, Some(&outside), WallModel::Rigid,
    ).unwrap();
    assert_eq!(response.impedance, direct.impedance);
    assert_eq!(holes, direct_holes);
    assert_eq!(holes.len(), 1);

    let receiver = CircularOutletReceiver { position_m: [0.0, 0.0, 1.0],
        radial_rings: 2, angular_points: 8, maximum_frequency_hz: 700.0 };
    let (mut bell, mut microphones, _) = duct.baffled_outlet_and_open_holes_and_shunt_delayed_at(
        8000, 256, receiver, WallModel::Rigid, outside.sound_speed,
    ).unwrap();
    assert_eq!(microphones.len(), 1, "changing to a continuous field must not drop a source");
    let info = microphones[0].info();
    assert_eq!(info.gas, outside);
    assert_eq!(info.radius_m, 0.0005);
    assert_eq!(info.receiver.position_m[0], LENGTH / 2.0);
    assert_eq!(info.axial_delay_samples,
        (field_owner.open_hole_sites()[0].one_way_time_s * 8000.0).floor() as usize);
    let mut hole_energy = 0.0;
    let mut difference = 0.0;
    for n in 0..1024 {
        let input = if n == 0 { 1.0 } else { 0.0 };
        let bell_pressure = bell.step(input).unwrap();
        let hole_pressure = microphones[0].step(input).unwrap();
        assert!(bell_pressure.is_finite() && hole_pressure.is_finite());
        hole_energy += hole_pressure * hole_pressure;
        difference += ((bell_pressure + hole_pressure) - bell_pressure).powi(2);
    }
    assert!(hole_energy > 0.0 && difference > 0.0);
}

#[test]
fn default_continuous_exterior_is_shared_by_feedback_holes_baked_load_and_microphones() {
    use crate::pcm_wav::baffled::CircularOutletReceiver;
    use fs_duct::{HoleState, RiccatiSteps, WallModel};
    use fs_duct::compose::{ComposeOptions, compose_chain_reflection};
    use fs_math::c64::C64;

    let radius = 0.0015;
    let geometry = Duct { segments: vec![
        Segment::Cylinder { radius, length: LENGTH / 2.0 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.003,
            bore_radius: radius, state: HoleState::Open },
        Segment::Cylinder { radius, length: LENGTH / 2.0 },
    ] };
    let source = field(360.0);
    let regional = RegionalDuct::from_continuous_gas_field(&geometry, &source, RiccatiSteps(16))
        .unwrap();
    let exterior = *source.state_at(LENGTH).unwrap().state();
    assert_ne!(source.state_at(LENGTH / 2.0).unwrap().state(), &exterior);
    let explicit = regional.clone().with_exterior_gas(exterior).unwrap();
    let owner = regional.continuous_gas_field().unwrap();
    let omega = TAU * 430.0;
    for termination in [Termination::Closed, Termination::FlangedOpen] {
        let actual = regional.response_with_wall_and_holes(
            omega, LossModel::Bessel, termination, WallModel::Rigid).unwrap();
        let expected = owner.response_with_wall_and_holes(
            omega, LossModel::Bessel, termination, Some(&exterior), WallModel::Rigid).unwrap();
        assert_eq!(actual.0.impedance, expected.0.impedance);
        assert_eq!(actual.1, expected.1);
        assert_eq!(regional.response(omega, LossModel::Bessel, termination, None)
            .unwrap().impedance, expected.0.impedance);
        assert_eq!(regional.reflection_with_wall(omega, LossModel::Bessel,
            termination, WallModel::Rigid).unwrap(), owner.reflection_with_wall(
            omega, LossModel::Bessel, termination, Some(&exterior), WallModel::Rigid)
            .unwrap().referenced_to(regional.inlet_impedance()).unwrap());
        // The lower-level owner deliberately retains its local-gas convention.
        let local = owner.response_with_wall_and_holes(
            omega, LossModel::Bessel, termination, None, WallModel::Rigid).unwrap();
        assert_ne!(local.1[0].u_hole_over_u_in, actual.1[0].u_hole_over_u_in);
    }
    let load = C64::from_re(300000.0);
    let chain = owner.chain_with_wall_exterior_controlled(omega, LossModel::Bessel,
        Some(&exterior), WallModel::Rigid, &mut || false).unwrap();
    let expected = compose_chain_reflection(&chain.blocks, Some(load), ComposeOptions::default())
        .unwrap().referenced_to(regional.inlet_impedance()).unwrap();
    assert_eq!(regional.loaded_reflection(omega, load, WallModel::Rigid).unwrap(), expected);

    let receiver = CircularOutletReceiver { position_m: [0.0, 0.0, 1.0],
        radial_rings: 2, angular_points: 8, maximum_frequency_hz: 700.0 };
    let microphones = |duct: &RegionalDuct| duct.baffled_outlet_and_open_holes_delayed_at(
        8000, 256, receiver, WallModel::Rigid, exterior.sound_speed).unwrap();
    let (mut default_bell, mut default_holes) = microphones(&regional);
    let (mut explicit_bell, mut explicit_holes) = microphones(&explicit);
    assert_eq!(default_holes.len(), 1);
    assert_eq!(default_holes[0].info().gas, exterior);
    for sample in 0..1024 {
        let input = if sample == 0 { 1.0 } else { 0.0 };
        assert_eq!(default_bell.step(input).unwrap().to_bits(),
            explicit_bell.step(input).unwrap().to_bits());
        assert_eq!(default_holes[0].step(input).unwrap().to_bits(),
            explicit_holes[0].step(input).unwrap().to_bits());
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
fn gas_model() -> ThermochemicalGasModel {
    gas_model_with_viscosity([1.8e-5, 1.8e-5])
}
// Viscosity knots at 200 K and 1000 K; the default above is flat in T.
fn gas_model_with_viscosity(viscosity_pa_s: [f64; 2]) -> ThermochemicalGasModel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
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
    let id = SpeciesId::new("analytic-field-gas").unwrap();
    let species = Nasa9StandardStateModelV1::new(
        id.clone(),
        MolarMass::new(0.028),
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
    for (key, dims, [low, high]) in [
        ("dynamic_viscosity", DynViscosity::DIMS, viscosity_pa_s),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.026; 2]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, low), (1000.0, high)],
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
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![SpeciesTransportMaterial {
            species: id,
            material,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        None,
    )
    .unwrap()
}
pub(crate) fn field(outlet_temperature: f64) -> AxialGasField {
    field_with_model(outlet_temperature, gas_model())
}
/// The same two-station field with a temperature-dependent viscosity claim.
pub(crate) fn field_with_viscosity(outlet_temperature: f64, viscosity_pa_s: [f64; 2])
    -> AxialGasField {
    field_with_model(outlet_temperature, gas_model_with_viscosity(viscosity_pa_s))
}
fn field_with_model(outlet_temperature: f64, model: ThermochemicalGasModel) -> AxialGasField {
    let state = |temperature_k| PrimitiveGasState {
        temperature_k,
        composition: Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    };
    AxialGasField::new(
        vec![0.0, LENGTH],
        101325.0,
        vec![state(300.0), state(outlet_temperature)],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}
pub(crate) fn cylinder() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: RADIUS,
            length: LENGTH,
        }],
    }
}
fn resonance(duct: &RegionalDuct, harmonic: usize) -> f64 {
    let nominal = harmonic as f64 * duct.inlet_gas().sound_speed / (4.0 * LENGTH);
    let inverse = |f| {
        duct.response(TAU * f, LossModel::Lossless, Termination::IdealOpen, None)
            // Re(1/Z)=0; its imaginary part crosses zero at an impedance peak.
            .map(|r| -1.0 / r.impedance.im)
    };
    let mut lo = nominal * 0.94;
    let mut hi = nominal * 1.02;
    let mut flo = inverse(lo).unwrap();
    assert!(flo * inverse(hi).unwrap() < 0.0);
    // 32 fixed bisections resolve cents far below the spatial ladder error.
    for _ in 0..32 {
        let mid = 0.5 * (lo + hi);
        // A midpoint on the exact lossless pole is refused as singular: the
        // bisection has landed on the peak itself.
        let fmid = match inverse(mid) {
            Ok(value) => value,
            Err(fs_duct::DuctError::Singular) => return mid,
            Err(error) => panic!("resonance bisection at {mid} Hz: {error}"),
        };
        if flo * fmid > 0.0 {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn uniform_field_keeps_original_duct_response_bitwise() {
    let field = field(300.0);
    let gas = *field.state_at(0.0).unwrap().state();
    for geometry in [
        cylinder(),
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: RADIUS,
                outlet_radius: 0.012,
                length: LENGTH,
            }],
        },
    ] {
        let original = RegionalDuct::new(&geometry, &[gas]).unwrap();
        let sampled = RegionalDuct::from_gas_field(&geometry, &field, 64).unwrap();
        assert_eq!(sampled.segment_gases(), &[gas]);
        assert_eq!(
            sampled.one_way_time_s().to_bits(),
            original.one_way_time_s().to_bits()
        );
        for frequency in [100.0, 350.0, 1000.0] {
            for loss in [LossModel::Lossless, LossModel::Bessel] {
                let a = original
                    .response(TAU * frequency, loss, Termination::UnflangedOpen, None)
                    .unwrap();
                let b = sampled
                    .response(TAU * frequency, loss, Termination::UnflangedOpen, None)
                    .unwrap();
                assert_eq!(a.impedance, b.impedance);
                assert_eq!(a.u_mouth_over_u_in, b.u_mouth_over_u_in);
                assert_eq!(a.p_mouth_over_u_in, b.p_mouth_over_u_in);
            }
        }
    }
}

#[test]
fn primitive_field_changes_register_intervals_with_second_order_convergence() {
    let field = field(299.0);
    let mut intervals = Vec::new();
    let mut fundamentals = Vec::new();
    for cells in [4, 8, 16, 32, 64] {
        let duct = RegionalDuct::from_gas_field(&cylinder(), &field, cells).unwrap();
        assert_eq!(duct.inlet_gas(), field.state_at(0.0).unwrap().state());
        assert_eq!(duct.exterior_gas(), field.state_at(LENGTH).unwrap().state());
        assert_ne!(duct.inlet_gas(), &duct.segment_gases()[0]);
        let f1 = resonance(&duct, 1);
        let f3 = resonance(&duct, 3);
        let interval = 1200.0 * (f3 / (3.0 * f1)).log2();
        intervals.push(interval);
        fundamentals.push(f1);
        let lossy = duct
            .response(
                TAU * (0.9 * f1),
                LossModel::Bessel,
                Termination::IdealOpen,
                None,
            )
            .unwrap();
        assert!(lossy.impedance.re >= 0.0);
        assert!(lossy.impedance.re.is_finite() && lossy.impedance.im.is_finite());
        eprintln!(
            "{{\"test\":\"gas_field_register\",\"cells\":{cells},\"f1_hz\":{f1},\"f3_hz\":{f3},\"interval_cents\":{interval}}}"
        );
    }
    let n = intervals.len();
    let coarse = (fundamentals[n - 3] - fundamentals[n - 2]).abs();
    let fine = (fundamentals[n - 2] - fundamentals[n - 1]).abs();
    let order = (coarse / fine).log2();
    assert!((1.7..2.3).contains(&order), "order {order}; {intervals:?}");
    assert!((intervals[0] - intervals[n - 1]).abs() > 1e-4);
    let extrapolated = intervals[n - 1] + (intervals[n - 1] - intervals[n - 2]) / 3.0;
    // First-order energy perturbation for T(x)=300K-x/L K at fixed K=gamma*p:
    // <x/L>_(u_q^2)=1/2+2/(q*pi)^2 for quarter-wave q=1,3.
    // delta log(f3/3f1)=(deltaT/300)*8/(9*pi^2).
    let expected =
        1200.0 / core::f64::consts::LN_2 / 300.0 * 8.0 / (9.0 * core::f64::consts::PI.powi(2));
    assert!(extrapolated > 0.0);
    assert!((extrapolated - expected).abs() <= (0.02 * expected).max(0.02));
    // A uniform temperature can match f1 but has exactly odd harmonic ratios;
    // this residual interval survives one global tuning action.
    let interval_error = (intervals[n - 2] - intervals[n - 1]).abs() / 3.0;
    assert!(extrapolated.abs() > 10.0 * interval_error);
    eprintln!(
        "{{\"test\":\"gas_field_register_order\",\"order\":{order},\"extrapolated_cents\":{extrapolated},\"perturbation_cents\":{expected}}}"
    );
}

#[test]
fn nonuniform_cone_keeps_geometry_and_exact_boundary_media() {
    let field = field(280.0);
    let geometry = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.003,
            outlet_radius: 0.012,
            length: LENGTH,
        }],
    };
    let duct = RegionalDuct::from_gas_field(&geometry, &field, 8).unwrap();
    assert_eq!(duct.segment_gases().len(), 8);
    assert_eq!(duct.inlet_gas().temperature, 300.0);
    assert_eq!(duct.exterior_gas().temperature, 280.0);
    let expected =
        duct.inlet_gas().characteristic_impedance / (core::f64::consts::PI * 0.003f64.powi(2));
    assert_eq!(duct.inlet_impedance(), expected);
    let response = duct
        .response(
            TAU * 400.0,
            LossModel::Bessel,
            Termination::UnflangedOpen,
            None,
        )
        .unwrap();
    assert!(
        (response.mouth_ka - TAU * 400.0 * 0.012 / duct.exterior_gas().sound_speed).abs() < 1e-14
    );
    assert!(response.impedance.re > 0.0);
}

#[test]
fn field_subdivision_and_geometry_refuse_before_rendering() {
    let field = field(299.0);
    for count in [0, 65, usize::MAX] {
        assert!(RegionalDuct::from_gas_field(&cylinder(), &field, count).is_err());
    }
    assert!(RegionalDuct::from_gas_field(&cylinder(), &field, 64).is_ok());
    let mismatch = Duct {
        segments: vec![Segment::Cylinder {
            radius: RADIUS,
            length: 0.56,
        }],
    };
    assert!(RegionalDuct::from_gas_field(&mismatch, &field, 8).is_err());
    let two = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH / 2.0
            };
            2
        ],
    };
    assert!(RegionalDuct::from_gas_field(&two, &field, 32).is_ok());
    assert!(RegionalDuct::from_gas_field(&two, &field, 33).is_err());
    let profile = fs_duct::RadiusProfile::exponential(RADIUS, 0.0, LENGTH).unwrap();
    let unsupported = Duct {
        segments: vec![Segment::Profile {
            profile: std::sync::Arc::new(profile),
            steps: fs_duct::RiccatiSteps(8),
        }],
    };
    assert!(RegionalDuct::from_gas_field(&unsupported, &field, 8).is_err());
}

#[test]
fn continuous_field_reaches_smooth_profile_and_exact_physical_ports() {
    let source = field(275.0);
    let geometry = Duct {
        segments: vec![Segment::Profile {
            profile: std::sync::Arc::new(
                fs_duct::RadiusProfile::exponential(RADIUS, 0.8, LENGTH).unwrap(),
            ),
            steps: fs_duct::RiccatiSteps(64),
        }],
    };
    let regional = RegionalDuct::from_continuous_gas_field(
        &geometry, &source, fs_duct::RiccatiSteps(32),
    ).unwrap();
    let owner = regional.continuous_gas_field().unwrap();
    assert_eq!(owner.steps_used(), 32);
    assert_eq!(regional.inlet_gas(), source.state_at(0.0).unwrap().state());
    assert_eq!(regional.outlet_gas(), source.state_at(LENGTH).unwrap().state());
    assert_eq!(regional.exterior_gas(), regional.outlet_gas());
    assert_eq!(regional.inlet_impedance(), owner.inlet_impedance());
    assert_eq!(regional.one_way_time_s(), owner.one_way_time_s());
    for termination in [Termination::IdealOpen, Termination::FlangedOpen] {
        let regional_response = regional.response(TAU * 510.0, LossModel::Bessel, termination, None)
            .unwrap();
        let direct = owner.response(TAU * 510.0, LossModel::Bessel, termination, None).unwrap();
        assert_eq!(regional_response.impedance, direct.impedance);
        assert_eq!(regional_response.u_mouth_over_u_in, direct.u_mouth_over_u_in);
        assert_eq!(regional_response.p_mouth_over_u_in, direct.p_mouth_over_u_in);
    }
}

#[test]
fn midpoint_field_is_a_second_order_rung_of_continuous_propagation() {
    let source = field(270.0);
    let continuous = RegionalDuct::from_continuous_gas_field(
        &cylinder(), &source, fs_duct::RiccatiSteps(256),
    ).unwrap();
    let omega = TAU * 700.0;
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let reference = continuous.response(omega, loss, Termination::IdealOpen, None)
            .unwrap().impedance;
        let errors: Vec<_> = [8, 16, 32, 64].into_iter().map(|n| {
            let midpoint = RegionalDuct::from_gas_field(&cylinder(), &source, n).unwrap();
            let got = midpoint.response(omega, loss, Termination::IdealOpen, None)
                .unwrap().impedance;
            (got - reference).abs() / reference.abs()
        }).collect();
        assert!(errors[0] > 1e-6, "the piecewise material approximation must be live: {errors:?}");
        for pair in errors.windows(2) {
            let order = fs_math::det::ln(pair[0] / pair[1]) / fs_math::det::ln(2.0);
            assert!((order - 2.0).abs() < 0.2, "midpoint material order={order}: {errors:?}");
        }
    }
}

#[test]
fn continuous_regional_field_requires_axial_end_segments() {
    let mut geometry = cylinder();
    geometry.segments.insert(0, Segment::ToneHole {
        hole_radius: 0.002,
        chimney_height: 0.003,
        bore_radius: RADIUS,
        state: fs_duct::HoleState::Open,
    });
    assert!(RegionalDuct::from_continuous_gas_field(
        &geometry, &field(280.0), fs_duct::RiccatiSteps(16),
    ).is_err());
}

#[test]
fn continuous_field_routes_local_tone_holes_into_the_characteristic_line() {
    let axial = Segment::Cylinder {
        radius: RADIUS,
        length: LENGTH / 3.0,
    };
    let hole = |state| Segment::ToneHole {
        hole_radius: 0.002,
        chimney_height: 0.005,
        bore_radius: RADIUS,
        state,
    };
    let geometry = Duct {
        segments: vec![
            axial.clone(),
            hole(fs_duct::HoleState::Closed),
            axial.clone(),
            hole(fs_duct::HoleState::Open),
            axial,
        ],
    };
    let source = field(280.0);
    let regional =
        RegionalDuct::from_continuous_gas_field(&geometry, &source, fs_duct::RiccatiSteps(16))
            .unwrap();
    let owner = regional.continuous_gas_field().unwrap();
    assert_eq!(owner.steps_used(), 50);
    assert_eq!(
        regional.segment_gases()[1],
        *source.state_at(LENGTH / 3.0).unwrap().state()
    );
    assert_eq!(
        regional.segment_gases()[3],
        *source.state_at(2.0 * LENGTH / 3.0).unwrap().state()
    );
    let response = regional
        .response(
            TAU * 430.0,
            LossModel::Bessel,
            Termination::UnflangedOpen,
            None,
        )
        .unwrap();
    let direct = owner
        .response(
            TAU * 430.0,
            LossModel::Bessel,
            Termination::UnflangedOpen,
            Some(regional.exterior_gas()),
        )
        .unwrap();
    assert_eq!(response.impedance, direct.impedance);
    assert_eq!(response.u_mouth_over_u_in, direct.u_mouth_over_u_in);
    let mut line = regional
        .characteristic_line(Termination::UnflangedOpen, 8000, 2048, None, Some(512))
        .unwrap();
    let mut energy = 0.0;
    for sample in 0..2048 {
        let incoming = line.push(if sample == 0 { 1.0 } else { 0.0 }).unwrap();
        assert!(incoming.is_finite());
        energy += incoming * incoming;
    }
    assert!(
        energy > 0.0,
        "the fingered field must actually feed pressure back"
    );
}

#[test]
fn material_wall_with_continuous_gas_reaches_feedback_and_outlet_pressure() {
    let source = field(270.0);
    let regional = RegionalDuct::from_continuous_gas_field(
        &cylinder(), &source, fs_duct::RiccatiSteps(32),
    ).unwrap();
    let wall = fs_phs::WallPin {
        surface_density: 0.01, stiffness_per_area: 1e8, resistance: 300.0,
    };
    let owner = regional.continuous_gas_field().unwrap();
    let actual = regional.response(TAU * 430.0, LossModel::Bessel,
        Termination::FlangedOpen, Some(&wall)).unwrap();
    let direct = owner.response_wall(TAU * 430.0, LossModel::Bessel,
        Termination::FlangedOpen, None, Some(&wall)).unwrap();
    assert_eq!(actual.impedance, direct.impedance);
    assert_eq!(actual.p_mouth_over_u_in, direct.p_mouth_over_u_in);

    let mut compliant = regional.characteristic_line(Termination::FlangedOpen,
        8000, 4801, Some(&wall), Some(512)).unwrap();
    let mut rigid = regional.characteristic_line(Termination::FlangedOpen,
        8000, 4801, None, Some(512)).unwrap();
    let receiver = crate::pcm_wav::baffled::CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0], radial_rings: 2, angular_points: 8,
        // The single-rate circular receiver admits up to 0.1 x the rate.
        maximum_frequency_hz: 800.0,
    };
    let mut compliant_mic = regional.baffled_outlet(8000, 512, receiver, Some(&wall)).unwrap();
    let mut rigid_mic = regional.baffled_outlet(8000, 512, receiver, None).unwrap();
    let first = compliant_mic.info().axial_delay_samples
        + compliant_mic.info().exterior_delay_samples.0;
    let (mut feedback_change, mut pressure_change, mut pressure_energy) = (0.0, 0.0, 0.0);
    for n in 0..1024 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let a = compliant.push(input).unwrap();
        let b = rigid.push(input).unwrap();
        let pa = compliant_mic.step(input).unwrap();
        let pb = rigid_mic.step(input).unwrap();
        assert!([a, b, pa, pb].iter().all(|x| x.is_finite()));
        if n < first { assert_eq!(pa, 0.0); }
        feedback_change += (a - b) * (a - b);
        pressure_change += (pa - pb) * (pa - pb);
        pressure_energy += pa * pa;
    }
    assert!(feedback_change > 1e-5);
    assert!(pressure_energy > 1e-12 && pressure_change > 1e-4 * pressure_energy);
}

mod steady;
