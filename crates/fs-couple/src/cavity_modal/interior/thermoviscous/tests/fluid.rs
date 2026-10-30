use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::FluidAcousticParameter;
use fs_material::fluid_acoustics::thermoviscous::{
    ResolvedThermoviscousLiquidState, resolve_thermoviscous_liquid_state,
    resolve_thermoviscous_liquid_state_with_expansion,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, VOLUMETRIC_EXPANSION_DIMS,
    VOLUMETRIC_EXPANSION_PROPERTY,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Velocity};

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
        .with("pressure", 101325.0)
        .unwrap()
}

// Manufactured complete source cards, never measured water/rosin constants.
fn liquid(
    alpha: Option<f64>,
    viscosity: f64,
    cv_scale: f64,
) -> (MaterialCard, ResolvedThermoviscousLiquidState) {
    let beta_s = 1.0 / (1000.0 * 1400.0_f64.powi(2));
    let beta_t = beta_s + alpha.unwrap_or(2e-4).powi(2) * 300.0 / (1000.0 * 4000.0);
    let cv = 4000.0 * beta_s / beta_t * cv_scale;
    let mut rows = vec![
        (DENSITY_PROPERTY, Density::DIMS, 1000.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            4000.0,
        ),
        ("specific_heat_capacity_cv", SPECIFIC_HEAT_CAPACITY_DIMS, cv),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.6,
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            viscosity,
        ),
        ("sound_speed", Velocity::DIMS, 1400.0),
    ];
    if let Some(alpha) = alpha {
        rows.push((
            VOLUMETRIC_EXPANSION_PROPERTY,
            VOLUMETRIC_EXPANSION_DIMS,
            alpha,
        ));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in rows {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        300.0,
                        300.0,
                    )
                    .with("pressure", 101325.0, 101325.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured complete fluid EOS".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured fluid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let parameter = FluidAcousticParameter::sound_speed().unwrap();
    let source = if alpha.is_some() {
        resolve_thermoviscous_liquid_state_with_expansion(
            &card,
            &point(),
            &parameter,
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
    } else {
        resolve_thermoviscous_liquid_state(
            &card,
            &point(),
            &parameter,
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
    }
    .unwrap();
    (card, source)
}

fn build_fluid<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    source: &ResolvedThermoviscousLiquidState,
) -> ThermoviscousFe<'c> {
    let (medium, thermodynamics) = source_fluid_medium(source, 101325.0, 0.002).unwrap();
    let geometry = HelmholtzFe::new(
        mesh,
        points,
        2,
        Medium {
            density_kg_m3: source.state().acoustic().density_kg_m3(),
            sound_speed_m_s: source.state().acoustic().sound_speed_m_s(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let end = |x: f64| {
        faces
            .iter()
            .copied()
            .filter(|face| {
                face.iter()
                    .all(|&v| (points[v as usize][0] - x).abs() < 1e-14)
            })
            .collect::<Vec<_>>()
    };
    let left = end(0.0);
    let right = end(LENGTH);
    let walls = faces
        .into_iter()
        .filter(|f| !left.contains(f) && !right.contains(f))
        .collect::<Vec<_>>();
    ThermoviscousFe::new_with_thermodynamics(
        mesh,
        points,
        medium,
        thermodynamics,
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: walls,
            ports: vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn original_liquid_card_reaches_resolved_fields_work_and_native_duct_scattering() {
    let (card, source) = liquid(Some(2e-4), 0.001, 1.0);
    assert_eq!(source.material().properties().len(), 7);
    assert_eq!(source.material().card_identity(), card.content_hash());
    assert_eq!(
        source
            .thermal_transport_state()
            .unwrap()
            .volumetric_expansion_per_k(),
        Some(2e-4)
    );
    let (mesh, xyz) = kuhn_cube(2);
    let xyz = points(&xyz);
    let model = ThermoviscousInterior::from_material(
        build_fluid(&mesh, &xyz, &source),
        source.clone(),
        BAND,
    )
    .unwrap();
    assert_eq!(model.material_source().identity(), source.identity());
    let drive = [C64::new(1.0, 0.2), C64::new(-0.5, 0.1)];
    let f = 400.0;
    let field = model.response(f, &drive, &mut || false).unwrap();
    assert!(field.power.viscous_w > 0.0);
    assert!(field.power.conductive_w >= 0.0);
    assert!(field.power.balance_residual_w().abs() < 2e-8 * field.power.dissipated_w());
    let resistance = source.state().acoustic().characteristic_impedance_pa_s_m() / (HEIGHT * WIDTH);
    let block = model.duct_two_port(f, resistance, &mut || false).unwrap();
    let native = model
        .inner()
        .scattering_with_heat_capacity(f, C64::from_re(4000.0), &[resistance; 2], &mut || false)
        .unwrap();
    for (a, b) in block.scattering.iter().zip(native) {
        assert_eq!(*a, b.conj());
    }
    let (_, viscous) = liquid(Some(2e-4), 0.008, 1.0);
    let changed =
        ThermoviscousInterior::from_material(build_fluid(&mesh, &xyz, &viscous), viscous, BAND)
            .unwrap();
    let other = changed.response(f, &drive, &mut || false).unwrap();
    let relative = (field.inward_flow_m3_s[0] - other.inward_flow_m3_s[0]).abs()
        / field.inward_flow_m3_s[0].abs();
    assert!(
        relative > 1e-3,
        "viscosity data must reach the actual traction/flow solve: {relative}"
    );
    assert!(model.response(f, &drive, &mut || true).is_err());
    assert!(
        model
            .response(BAND[1] + 1.0, &drive, &mut || false)
            .is_err()
    );
    assert!(
        source
            .at_query_point(&card, &point().with("pressure", 2e5).unwrap())
            .is_err()
    );
    assert_eq!(source.material().query(), point());
}

#[test]
fn source_fluid_refuses_unknown_expansion_inconsistent_calorics_and_wrong_fe() {
    let (_, unknown) = liquid(None, 0.001, 1.0);
    assert!(source_fluid_medium(&unknown, 101325.0, 0.0).is_err());
    let (_, inconsistent) = liquid(Some(2e-4), 0.001, 0.95);
    assert!(source_fluid_medium(&inconsistent, 101325.0, 0.0).is_err());
    for alpha in [0.0, -2e-4] {
        let (_, source) = liquid(Some(alpha), 0.001, 1.0);
        let (medium, thermo) = source_fluid_medium(&source, 101325.0, 0.0).unwrap();
        assert_eq!(thermo.isobaric_expansion_per_k, alpha);
        assert!((thermo.sound_speed_m_s(medium).unwrap() / 1400.0 - 1.0).abs() < 2e-15);
    }
    let (_, source) = liquid(Some(2e-4), 0.001, 1.0);
    let (_, wrong) = liquid(Some(2e-4), 0.002, 1.0);
    let (mesh, xyz) = kuhn_cube(1);
    let xyz = points(&xyz);
    assert!(
        ThermoviscousInterior::from_material(build_fluid(&mesh, &xyz, &source), wrong, BAND)
            .is_err()
    );
    let original = gas();
    let (medium, thermo) = source_fluid_medium(&original, original.pressure, 0.0).unwrap();
    assert!((thermo.sound_speed_m_s(medium).unwrap() / original.sound_speed - 1.0).abs() < 2e-14);
}
