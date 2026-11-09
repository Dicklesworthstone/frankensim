//! G1/G4 original material publication and independent nozzle acceptance.
use super::super::nozzle::{
    EquilibriumNozzleControl,
    material::{
        MaterialEquilibriumNozzleControl, MaterialEquilibriumNozzleError,
        equilibrium_material_nozzle_flow,
    },
};
use super::*;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::{
    SpeciesTransportMaterial, ThermochemicalGasModel,
    energy::equilibrium::pressure::GasEquilibriumEnthalpySolve,
};
use fs_qty::{
    Area, DynViscosity, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1,
    transport::kinetic::ThermalConductivityQuantityV1,
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            T0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), P0)
        .unwrap()
        .with("sample_id", 7.0)
        .unwrap()
}
fn bank(
    chemistry: &IdealGasEquilibriumV1,
    transport_minimum_temperature: f64,
) -> ThermochemicalGasModel {
    let sources = chemistry.models();
    assert!(sources.len() == 1 || sources.len() == 2);
    let transport = sources
        .iter()
        .map(|source| {
            let mut claims = ClaimSet::new();
            for (name, dims, value) in [
                ("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
                (
                    "thermal_conductivity",
                    ThermalConductivityQuantityV1::DIMS,
                    0.025,
                ),
            ] {
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(name, dims),
                        value: PropertyValue::Scalar { value, dims },
                        validity: ValidityDomain::unconstrained()
                            .with_quantity(
                                "temperature",
                                quantity(QuantityKind::AbsoluteTemperature),
                                transport_minimum_temperature,
                                3000.0,
                            )
                            .with_quantity(
                                "pressure",
                                quantity(QuantityKind::Pressure),
                                100.0,
                                1e8,
                            ),
                        uncertainty: UncertaintyModel::Unstated,
                        interpolation: InterpolationPolicy::ConstantWithinValidity,
                        observations: Vec::new(),
                        provenance: Provenance {
                            source: "manufactured reacting nozzle material".into(),
                            license: "test fixture".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
            SpeciesTransportMaterial {
                species: source.species().clone(),
                material: MaterialCard::assemble(
                    MaterialStateId {
                        chemistry: source.species().as_str().into(),
                        phase: "gas".into(),
                        process: "reacting-nozzle-conformance".into(),
                        revision: 0,
                    },
                    claims,
                    Vec::new(),
                )
                .unwrap(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            sources.to_vec(),
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 / sources.len() as f64; sources.len()],
            )
            .unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap()
}
fn material_controls() -> MaterialEquilibriumNozzleControl {
    let mut native = controls();
    native.isentropic.entropy_tolerance_j_k = 1e-9;
    native.enthalpy.enthalpy_tolerance_j = 1e-9;
    MaterialEquilibriumNozzleControl {
        nozzle: EquilibriumNozzleControl::new(native.isentropic, Pressure::new(100.0)),
        material: GasEquilibriumEnthalpySolve {
            equilibrium: native.enthalpy,
            mass_tolerance_kg: 1e-12,
            enthalpy_tolerance_j: 1e-7,
        },
        relative_endpoint_tolerance: 1e-8,
    }
}

#[test]
fn g1_reacting_material_flow_retains_cards_mass_and_reacting_sonic_closure() {
    let chemistry = dissociation();
    let bank = bank(&chemistry, 200.0);
    let before = bank.clone();
    let context = context();
    let feed = [(id("A2"), Amount::new(1.0))];
    for back in [0.9 * P0, 0.0] {
        let flow = equilibrium_material_nozzle_flow(
            &bank,
            &chemistry,
            &context,
            &feed,
            Pressure::new(back),
            Area::new(0.002),
            material_controls(),
            || false,
        )
        .unwrap();
        let parcel = flow.outlet.parcel();
        assert_eq!(parcel.mass_kg(), flow.inlet.parcel().mass_kg());
        close(parcel.mass_kg(), 0.002, 1e-16);
        assert_eq!(parcel.point().axes()["sample_id"], 7.0);
        assert_eq!(
            parcel.model().thermodynamics().components(),
            chemistry.models()
        );
        assert_eq!(
            parcel.model().transport_materials(),
            bank.transport_materials()
        );
        close(
            flow.mass_flow.value(),
            0.002 * parcel.state().state().density * flow.velocity_m_s,
            1e-14,
        );
        assert!(
            flow.balance.stagnation_enthalpy_residual_j_kg.abs()
                <= material_controls().nozzle.absolute_energy_tolerance_j_kg
        );
        assert!(
            flow.balance.entropy_residual_j_k.abs()
                <= material_controls().nozzle.isentropic.entropy_tolerance_j_k
        );
        assert!(
            flow.balance.relative_mass_flux_residual.abs()
                <= material_controls().relative_endpoint_tolerance
        );
        assert!(
            flow.balance.relative_density_residual.abs()
                <= material_controls().relative_endpoint_tolerance
        );
        assert!(
            flow.outlet
                .equilibrium()
                .equilibrium()
                .is_bound_to(&chemistry)
        );
        check_products(flow.outlet.equilibrium(), &chemistry);
        let response = flow.outlet.native_sample().response;
        assert!(
            response.heat_capacity_p_j_k
                > parcel.mass_kg() * parcel.state().state().specific_heat_cp
        );
        assert_eq!(flow.material_pressure_evaluations, 2);
        assert!(flow.composition_iterations >= flow.isentropic.composition_iterations);
        if flow.isentropic.choked {
            close(flow.mach_number, 1.0, 2e-8);
            assert!((response.sound_speed_m_s - parcel.state().state().sound_speed).abs() > 1.0);
        } else {
            assert!(flow.mach_number < 1.0);
        }
    }
    assert_eq!(bank, before);
}

#[test]
fn g4_material_transport_and_loose_hp_do_not_bypass_original_flow_acceptance() {
    let chemistry = dissociation();
    let limited = bank(&chemistry, 900.0);
    let complete = bank(&chemistry, 200.0);
    let context = context();
    let feed = [(id("A2"), Amount::new(1.0))];
    let call = |bank, c| {
        equilibrium_material_nozzle_flow(
            bank,
            &chemistry,
            &context,
            &feed,
            Pressure::new(0.0),
            Area::new(0.002),
            c,
            || false,
        )
    };
    assert!(matches!(
        call(&limited, material_controls()),
        Err(MaterialEquilibriumNozzleError::Material(_))
    ));
    let mut c = material_controls();
    c.material.equilibrium.equilibrium.maximum_iterations -= 1;
    assert!(matches!(
        call(&complete, c),
        Err(MaterialEquilibriumNozzleError::Invalid(_))
    ));
    c = material_controls();
    c.material.equilibrium.enthalpy_tolerance_j = 1e6;
    c.material.enthalpy_tolerance_j = 1e6;
    // A loose HP endpoint can close its own declared tolerance while violating
    // the unchanged chemical isentrope, reacting sonic state and flow rate.
    assert!(matches!(
        call(&complete, c),
        Err(MaterialEquilibriumNozzleError::Endpoint(_))
    ));
}

#[test]
fn g0_material_zero_flow_keeps_the_exact_reference_without_an_hp_exit_query() {
    let chemistry = dissociation();
    let bank = bank(&chemistry, 200.0);
    for (back, area) in [(P0, 0.001), (0.0, 0.0)] {
        let flow = equilibrium_material_nozzle_flow(
            &bank,
            &chemistry,
            &context(),
            &[(id("A2"), Amount::new(1.0))],
            Pressure::new(back),
            Area::new(area),
            material_controls(),
            || false,
        )
        .unwrap();
        assert_eq!(flow.inlet, flow.outlet);
        assert_eq!(flow.mass_flow.value(), 0.0);
        assert_eq!(flow.material_pressure_evaluations, 1);
        assert_eq!(flow.isentropic.pressure_evaluations, 1);
        assert_eq!(flow.balance.entropy_residual_j_k, 0.0);
        assert_eq!(flow.balance.stagnation_enthalpy_residual_j_kg, 0.0);
    }
}

#[test]
fn g4_final_material_nozzle_cancellation_preserves_borrowed_sources_and_context() {
    let chemistry = dissociation();
    let bank = bank(&chemistry, 200.0);
    let bank_before = bank.clone();
    let context = context();
    let context_before = context.clone();
    let feed = [(id("A2"), Amount::new(1.0))];
    let call = |cancelled: &mut dyn FnMut() -> bool| {
        equilibrium_material_nozzle_flow(
            &bank,
            &chemistry,
            &context,
            &feed,
            Pressure::new(0.0),
            Area::new(0.002),
            material_controls(),
            cancelled,
        )
    };
    let mut total = 0;
    let expected = call(&mut || {
        total += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        call(&mut || {
            count += 1;
            count == total
        }),
        Err(MaterialEquilibriumNozzleError::Cancelled)
    ));
    assert_eq!(count, total);
    assert_eq!(call(&mut || false).unwrap(), expected);
    assert_eq!(bank, bank_before);
    assert_eq!(context, context_before);
}
