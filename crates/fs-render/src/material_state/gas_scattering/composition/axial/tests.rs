//! Synthetic original-owner controls, not measured-gas accuracy claims.
use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::gas_field::PrimitiveGasState;
use fs_material::gas_scattering::{
    MolecularRayleigh,
    composition::MolecularCauchyDispersion,
    thermochemical::{SpeciesMolecularRayleigh, ThermochemicalRayleighState},
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_qty::{DynViscosity, MolarMass, Pressure};
use fs_thermochem::{
    ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
    StandardStateConventionV1, StandardStatePhaseV1, mixture::FrozenIdealGasMixtureModelV1,
    required_nasa9_parameter_dims,
};

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic axial gas optics control".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn original_model(water: bool) -> ThermochemicalGasModel {
    let mut species = Vec::new();
    let mut transport = Vec::new();
    for (name, mass) in [("A", 0.028), ("B", 0.018)] {
        let mut parameters = std::collections::BTreeMap::new();
        for (i, (key, dims)) in required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .enumerate()
        {
            parameters.insert(
                key.to_owned(),
                LawParameter {
                    value: if i == 2 { 3.5 } else { 0.0 },
                    dims,
                },
            );
        }
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: 100_000.0,
                dims: Pressure::DIMS,
            },
        );
        let card = ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 1600.0),
            sources: vec![],
            provenance: provenance(),
        };
        species.push(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new(name).unwrap(),
                MolarMass::new(mass),
                StandardStateConventionV1::new(
                    StandardStatePhaseV1::Gas,
                    ReferenceEquationOfStateV1::IdealGas,
                    Pressure::new(100_000.0),
                    ElementalReferenceIdV1::new("synthetic-reference").unwrap(),
                )
                .unwrap(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap(),
        );
        let mut claims = ClaimSet::new();
        for (property, dims, value) in [
            ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(property, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        transport.push(SpeciesTransportMaterial {
            species: SpeciesId::new(name).unwrap(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: name.into(),
                    phase: "gas".into(),
                    process: "synthetic".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
    }
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        transport,
        water.then(|| SpeciesId::new("B").unwrap()),
    )
    .unwrap()
}
fn field(stations: &[(f64, f64)], cellwise: bool, phase: ThermochemicalPhase) -> AxialGasField {
    let positions = (0..stations.len()).map(|i| i as f64).collect();
    let states = stations
        .iter()
        .map(|&(temperature_k, a)| PrimitiveGasState {
            temperature_k,
            composition: Composition::new(CompositionBasis::MoleFraction, vec![a, 1.0 - a])
                .unwrap(),
        })
        .collect();
    let constructor = if cellwise {
        AxialGasField::new_piecewise_constant
    } else {
        AxialGasField::new
    };
    constructor(
        positions,
        101_325.0,
        states,
        original_model(phase != ThermochemicalPhase::FrozenGas),
        QueryPoint::new(),
        phase,
    )
    .unwrap()
}
fn optics(field: &AxialGasField, position: f64, dispersive: bool) -> MolecularRefractionModel {
    let state = field.state_at(position).unwrap();
    let rows: Vec<_> = state
        .thermodynamics()
        .receipt()
        .species()
        .iter()
        .map(|species| {
            let (r, delta, tag) = if species.as_str() == "A" {
                (0.00015, 0.03, 112)
            } else {
                (0.00065, 0.18, 113)
            };
            SpeciesMolecularRayleigh {
                species: species.clone(),
                optical: MolecularRayleigh::try_new(
                    r,
                    2.5e25,
                    550.0,
                    delta,
                    ContentHash([tag; 32]),
                )
                .unwrap(),
            }
        })
        .collect();
    let model = MolecularRefractionModel::from_thermochemical(
        &ThermochemicalRayleighState::try_new(state, &rows).unwrap(),
    )
    .unwrap();
    if !dispersive {
        return model;
    }
    let spectra: Vec<_> = model
        .components()
        .iter()
        .map(|component| {
            MolecularCauchyDispersion::try_new(
                component.optical(),
                2e-6,
                1e-7,
                [380.0, 780.0],
                ContentHash([114; 32]),
            )
            .unwrap()
        })
        .collect();
    model.with_dispersion(&spectra).unwrap()
}
fn bind<'a>(
    field: &'a AxialGasField,
    model: &'a MolecularRefractionModel,
    span: usize,
) -> AxialGasCompositionField<'a> {
    AxialGasCompositionField::try_new(
        field,
        model,
        span,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        ContentHash([115; 32]),
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_original_temperature_composition_and_eos_supply_complete_partial_density_jets() {
    let original = field(
        &[(300.0, 1.0), (900.0, 0.0)],
        false,
        ThermochemicalPhase::FrozenGas,
    );
    let model = optics(&original, 0.5, false);
    let provider = bind(&original, &model, 0);
    assert!(std::ptr::eq(provider.source(), &original));
    assert_eq!(provider.axial_bounds_m(), [0.0, 1.0]);
    for x in [0.0, 0.125, 0.5, 0.875, 1.0] {
        let point = Point3::new(x, 0.0, 0.0);
        let state = original.state_at(x).unwrap();
        assert_eq!(provider.sample_state(point).unwrap(), state);
        assert_eq!(
            provider.temperature_k(point).unwrap(),
            state.state().temperature
        );
        let mut jets = [UNWRITTEN_JET; 2];
        provider
            .sample_partial_density_jets(point, &mut jets)
            .unwrap();
        let number = AVOGADRO_PER_MOL / state.equation_of_state().molar_volume().value();
        let t = 300.0 + 600.0 * x;
        for (i, (fraction, slope)) in [(1.0 - x, -1.0), (x, 1.0)].into_iter().enumerate() {
            let gradient = number * (slope - fraction * 600.0 / t);
            close(jets[i].sample.density, number * fraction, number * 2e-15);
            close(jets[i].sample.gradient_per_m.x, gradient, number * 2e-14);
            close(
                jets[i].hessian_diagonal_per_m2[0],
                -1200.0 / t * gradient,
                number * 8e-14,
            );
            assert_eq!(jets[i].sample.gradient_per_m.y, 0.0);
            assert_eq!(jets[i].hessian_diagonal_per_m2[1..], [0.0, 0.0]);
            assert!(jets[i].sample.density <= provider.component_number_density_upper_bounds()[i]);
        }
        if x == 0.0 {
            assert_eq!(jets[1].sample.density, 0.0);
            assert!(
                jets[1].sample.gradient_per_m.x > 0.0,
                "absent species can enter the original span"
            );
        }
    }
}

#[test]
fn g1_original_span_feeds_shared_spectral_refraction_and_molecular_scattering() {
    let original = field(
        &[(300.0, 0.75), (900.0, 0.25)],
        false,
        ThermochemicalPhase::FrozenGas,
    );
    let model = optics(&original, 0.5, true);
    let provider = bind(&original, &model, 0);
    let refractive = provider.bind_refraction().unwrap();
    let scattering = refractive
        .bind_scattering_field(provider.component_number_density_upper_bounds())
        .unwrap();
    for wavelength in [380.0, 550.0, 780.0] {
        let prepared = refractive.at_wavelength(wavelength).unwrap();
        let spectral = model
            .at_wavelength(wavelength, PartialDensityBasis::MoleculesPerM3)
            .unwrap();
        let majorant = scattering
            .scattering_majorant_optics(wavelength)
            .unwrap()
            .scattering_per_m;
        for i in 0..=32 {
            let x = f64::from(i) / 32.0;
            let p = Point3::new(x, 0.0, 0.0);
            let state = original.state_at(x).unwrap();
            let number = AVOGADRO_PER_MOL / state.equation_of_state().molar_volume().value();
            let partials: Vec<_> = state
                .thermodynamics()
                .receipt()
                .mole_fraction_bits()
                .iter()
                .map(|bits| number * f64::from_bits(*bits))
                .collect();
            let expected = spectral.refractive_index(&partials).unwrap();
            let sample = prepared.sample(p, wavelength).unwrap();
            assert_eq!(sample.phase_index(), expected.phase_index);
            assert_eq!(sample.group_index(), Some(expected.group_index));
            assert!(sample.group_index().unwrap() > sample.phase_index());
            assert_eq!(sample.absorption_per_m(), 0.0);
            assert_eq!(
                prepared
                    .sample_second_order(p, wavelength)
                    .unwrap()
                    .optical(),
                sample
            );
            let actual = scattering.sample_scattering_optics(p, wavelength).unwrap();
            assert_eq!(actual, spectral.scattering(&partials).unwrap());
            close(
                actual.mass_density_kg_m3,
                state.state().density,
                state.state().density * 5e-14,
            );
            assert!(actual.scattering_per_m <= majorant);
        }
    }
}

#[test]
fn g1_native_curved_ray_bundles_focus_with_the_same_analytic_thermochemical_jet() {
    let original = field(
        &[(300.0, 0.75), (900.0, 0.25)],
        false,
        ThermochemicalPhase::FrozenGas,
    );
    let model = optics(&original, 0.5, true);
    let provider = bind(&original, &model, 0);
    let refractive = provider.bind_refraction().unwrap();
    let ray = Ray {
        origin: Point3::new(0.5, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    };
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let options = GradientIndexOptions {
        max_step_m: 0.25,
        relative_tolerance: 1e-10,
        ..Default::default()
    };
    let length = 30.0;
    for wavelength in [380.0, 780.0] {
        let prepared = refractive.at_wavelength(wavelength).unwrap();
        let bundle = propagate_bundle(
            &prepared,
            ray,
            axis,
            wavelength,
            length,
            BundleOptions {
                integration: options,
                ..Default::default()
            },
        )
        .unwrap();
        let path = bundle.propagation();
        assert!(
            path.ray.origin.x < 0.49,
            "the physical thermal/composition gradient must bend the ray"
        );
        assert!(path.group_optical_path_length_m.unwrap() > path.optical_path_length_m);
        assert_eq!(path.optical_depth, 0.0);
        let projection = bundle
            .project_at_endpoint(path.ray.dir, ProjectionOptions::default())
            .unwrap();
        assert!((projection.area_per_solid_angle_m2() - length * length).abs() > 1.0);
        // Neighboring complete native solves independently check the transported
        // bundle variation; finite differences are used only in this control.
        let epsilon = 1e-5;
        let endpoints: Vec<_> = [-epsilon, epsilon]
            .into_iter()
            .map(|angle| {
                propagate(
                    &prepared,
                    Ray {
                        origin: ray.origin,
                        dir: Vec3::new(angle, 0.0, 1.0),
                    },
                    axis,
                    wavelength,
                    length,
                    options,
                )
                .unwrap()
                .ray
                .origin
            })
            .collect();
        let derivative = endpoints[1].delta_from(endpoints[0]).scale(0.5 / epsilon);
        let actual = bundle.angular_derivatives()[0].position_m_per_rad;
        close(actual.x, derivative.x, 3e-5);
        close(actual.y, derivative.y, 3e-5);
        close(actual.z, derivative.z, 3e-5);
    }
}

#[test]
fn g0_original_knots_are_explicit_domains_and_unsupported_physics_refuses() {
    let original = field(
        &[(300.0, 0.75), (900.0, 0.25), (300.0, 0.75)],
        false,
        ThermochemicalPhase::FrozenGas,
    );
    let model = optics(&original, 0.5, false);
    let left = bind(&original, &model, 0);
    let right = bind(&original, &model, 1);
    let a = left.bind_refraction().unwrap();
    let b = right.bind_refraction().unwrap();
    let knot = Point3::new(1.0, 0.0, 0.0);
    let sa = a.sample(knot, 550.0).unwrap();
    let sb = b.sample(knot, 550.0).unwrap();
    assert_eq!(sa.phase_index(), sb.phase_index());
    assert_ne!(sa.phase_gradient_per_m(), sb.phase_gradient_per_m());
    assert_eq!(
        a.sample(Point3::new(1.1, 0.0, 0.0), 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert_eq!(
        right.temperature_k(Point3::new(0.9, 0.0, 0.0)),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert_ne!(left.content_identity(), right.content_identity());
    let point = Point3::new(0.0, 0.0, 0.0);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let refusal = |source, index, direction, identity| {
        AxialGasCompositionField::try_new(source, &model, index, point, direction, identity).err()
    };
    assert!(matches!(
        refusal(&original, 2, axis, ContentHash([1; 32])),
        Some(AxialGasCompositionError::InvalidSpan)
    ));
    assert!(matches!(
        refusal(&original, 0, axis, ContentHash([0; 32])),
        Some(AxialGasCompositionError::MissingSnapshotIdentity)
    ));
    assert!(matches!(
        refusal(&original, 0, Vec3::new(0.0, 0.0, 0.0), ContentHash([1; 32])),
        Some(AxialGasCompositionError::InvalidEmbedding)
    ));
    let cellwise = field(
        &[(300.0, 0.75), (900.0, 0.25)],
        true,
        ThermochemicalPhase::FrozenGas,
    );
    assert!(matches!(
        refusal(&cellwise, 0, axis, ContentHash([1; 32])),
        Some(AxialGasCompositionError::UnsupportedReconstruction)
    ));
    let phase = field(
        &[(300.0, 1.0), (310.0, 1.0)],
        false,
        ThermochemicalPhase::LiquidWaterLimit,
    );
    assert!(matches!(
        refusal(&phase, 0, axis, ContentHash([1; 32])),
        Some(AxialGasCompositionError::UnsupportedPhase)
    ));
    let pure = field(
        &[(300.0, 1.0), (900.0, 0.0)],
        false,
        ThermochemicalPhase::FrozenGas,
    );
    let incomplete = optics(&pure, 0.0, false);
    assert!(matches!(
        AxialGasCompositionField::try_new(&pure, &incomplete, 0, point, axis, ContentHash([1; 32])),
        Err(AxialGasCompositionError::ModelMismatch)
    ));
}
