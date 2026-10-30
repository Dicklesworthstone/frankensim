//! G0/G3: physical acoustic scaling uses one material's sound and transport.
use fs_evidence::ValidityDomain;
use fs_lbm::lattice3::{AcousticScalingPlan, acoustic_scaling};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{
    FluidAcousticParameter, FluidAcousticState, FluidAcoustics, resolve_acoustic_liquid_state,
};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Velocity};

fn card(sound_speed: f64, viscosity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, quantity, value) in [
        ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
        (
            "specific_heat_capacity",
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            4000.0,
        ),
        (
            "thermal_conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            0.6,
        ),
        (
            "dynamic_viscosity",
            QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS),
            viscosity,
        ),
        (
            "sound_speed",
            QuantitySpec::dimensional(Velocity::DIMS),
            sound_speed,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured acoustic liquid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn source_liquid_swap_changes_physical_timestep_and_relaxation() {
    let mut plans = Vec::new();
    for (speed, viscosity) in [(1500.0, 0.001), (1000.0, 0.004)] {
        let source = card(speed, viscosity);
        let liquid = resolve_acoustic_liquid_state(
            &source,
            &QueryPoint::new(),
            &FluidAcousticParameter::sound_speed().unwrap(),
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(liquid.material().card_identity(), source.content_hash());
        let plan = AcousticScalingPlan::from_material(&liquid, 1e-4, 1.0, 0.05).unwrap();
        assert_eq!(plan.sound_speed_m_s, speed);
        assert_eq!(plan.nu_m2_s, viscosity / 1000.0);
        assert!((plan.dt_s * speed / (fs_lbm::CS2.sqrt() * 1e-4) - 1.0).abs() < 1e-14);
        plans.push(plan);
    }
    assert!((plans[1].dt_s / plans[0].dt_s - 1.5).abs() < 1e-14);
    assert!((plans[1].lattice.tau_margin / plans[0].lattice.tau_margin - 6.0).abs() < 1e-13);
}

#[test]
fn gas_entry_reuses_coherent_material_admission_and_preserves_accepted_bits() {
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let expected = acoustic_scaling(
        gas.dynamic_viscosity / gas.density,
        gas.sound_speed,
        1e-4,
        1.0,
        0.05,
    )
    .unwrap();
    let plan = AcousticScalingPlan::from_gas(&gas, 1e-4, 1.0, 0.05).unwrap();
    assert_eq!(plan.dt_s.to_bits(), expected.dt_s.to_bits());
    assert_eq!(plan.lattice.tau.to_bits(), expected.lattice.tau.to_bits());
    assert_eq!(
        plan.lattice.viscosity.to_bits(),
        expected.lattice.viscosity.to_bits()
    );
    gas.prandtl *= 1.01;
    assert!(AcousticScalingPlan::from_gas(&gas, 1e-4, 1.0, 0.05).is_err());
}

#[test]
fn consumer_refuses_mismatched_acoustic_and_transport_density() {
    struct Inconsistent;
    impl FluidAcoustics for Inconsistent {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            FluidAcousticState::try_new(1000.0, 1500.0)
        }
    }
    impl FluidTransport for Inconsistent {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            FluidTransportState::try_new(999.0, 4000.0, 0.6, 0.001)
        }
    }
    assert!(AcousticScalingPlan::from_material(&Inconsistent, 1e-4, 1.0, 0.05).is_err());
}
