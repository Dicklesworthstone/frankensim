//! G1/G3: original liquid laws reach native lossy pressure/flow propagation.
//! Source curves are manufactured numerical fixtures, not measured water data.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::FluidAcousticParameter;
use fs_material::fluid_acoustics::thermoviscous::{
    ResolvedThermoviscousLiquidState, ThermoviscousFluid, ThermoviscousFluidState,
    resolve_thermoviscous_liquid_state_with_requirements,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, LiquidPropertyRequirements};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Velocity};

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("phase-liquid", 1.0)
        .unwrap()
}

fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}

fn source_liquid() -> (MaterialCard, ResolvedThermoviscousLiquidState) {
    let heat = QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS);
    let conductivity = semantic(QuantityKind::ThermalConductivity);
    let mut claims = ClaimSet::new();
    for (name, quantity, values) in [
        (
            "rho",
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, 920.0],
        ),
        ("cp", heat, [4000.0, 4100.0]),
        ("cv", heat, [3900.0, 3700.0]),
        ("k", conductivity, [0.6, 1.2]),
        (
            "mu",
            QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS),
            [0.001, 0.004],
        ),
        (
            "c",
            QuantitySpec::dimensional(Velocity::DIMS),
            [1500.0, 1300.0],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        semantic(QuantityKind::AbsoluteTemperature),
                        300.0,
                        400.0,
                    )
                    .with("pressure", 100_000.0, 100_000.0)
                    .with("phase-liquid", 1.0, 1.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured duct liquid source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical fluid".into(),
            phase: "liquid".into(),
            process: "test fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirements = LiquidPropertyRequirements::try_new(
        requirement("rho", Density::DIMS),
        requirement("cp", SPECIFIC_HEAT_CAPACITY_DIMS),
        ScalarPropertyRequirement::try_with_quantity(
            "k",
            conductivity,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        requirement("mu", DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let material = resolve_thermoviscous_liquid_state_with_requirements(
        &card,
        &point(300.0),
        &FluidAcousticParameter::SoundSpeed(requirement("c", Velocity::DIMS)),
        &requirements,
        &requirement("cv", SPECIFIC_HEAT_CAPACITY_DIMS),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, material)
}

fn lossy<M: ThermoviscousFluid + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    loss: LossModel,
    load: Option<C64>,
    wall: WallModel<'_>,
) -> ChainResult {
    fluid::input_impedance_load_with_loss(
        duct,
        material,
        omega,
        loss,
        load,
        wall,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap()
}

#[test]
fn original_liquid_temperature_and_transport_change_actual_pipe_transmission() {
    let (card, cold) = source_liquid();
    let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
    let original = cold.clone();
    let (radius, length, omega) = (0.004_f64, 0.7_f64, 800.0_f64);
    let area = core::f64::consts::PI * radius * radius;
    let duct = cylinder(radius, length);
    let mut transmissions = Vec::new();
    for source in [&cold, &hot] {
        let state = source.thermoviscous_state().unwrap();
        let rho = state.acoustic().density_kg_m3();
        let c = state.acoustic().sound_speed_m_s();
        let transport = state.transport();
        let mu = transport.dynamic_viscosity_pa_s();
        let pr = mu * transport.specific_heat_j_kg_k() / transport.thermal_conductivity_w_m_k();
        let layer = (mu / (2.0 * rho * omega)).sqrt() / radius;
        let thermal = (state.heat_capacity_ratio() - 1.0) / pr.sqrt();
        // Independent Kirchhoff first-order line solution; no native wave or
        // transfer helper supplies these reference coefficients.
        let delta_k = layer * (1.0 + thermal);
        let delta_z = layer * (1.0 - thermal);
        let k = C64::new(1.0 + delta_k, delta_k).scale(omega / c);
        let z = C64::new(1.0 + delta_z, delta_z).scale(rho * c / area);
        let phase = k.re * length;
        let expected = C64::new(phase.cos(), phase.sin()).scale((-k.im * length).exp());
        let actual = lossy(
            &duct,
            source,
            omega,
            LossModel::WideTube,
            Some(z),
            WallModel::Rigid,
        );
        near(actual.z_in, z, z.abs(), 4e-12);
        near(actual.u_out_over_u_in, expected, 1.0, 4e-12);
        near(actual.p_out_over_u_in, z * expected, z.abs(), 4e-12);
        let absorbed = actual.z_in.re - z.re * actual.u_out_over_u_in.abs().powi(2);
        let expected_absorbed = z.re * (1.0 - (-2.0 * k.im * length).exp());
        assert!(absorbed > 0.0);
        assert!((absorbed - expected_absorbed).abs() < 5e-12 * z.re);
        transmissions.push(actual.u_out_over_u_in);
        for property in source.material().properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                cold.material()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    assert!((transmissions[0] - transmissions[1]).abs() > 0.04);
    assert!(transmissions[1].abs() < transmissions[0].abs());
    assert!(cold.at_query_point(&card, &point(401.0)).is_err());
    assert_eq!(cold, original);
}

#[test]
fn gas_and_original_liquid_share_poiseuille_inertia_and_isothermal_compliance_limits() {
    let (card, liquid) = source_liquid();
    let hot = liquid.at_query_point(&card, &point(400.0)).unwrap();
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    for material in [&air as &dyn ThermoviscousFluid, &liquid, &hot] {
        let state = material.thermoviscous_state().unwrap();
        let rho = state.acoustic().density_kg_m3();
        let c = state.acoustic().sound_speed_m_s();
        let mu = state.transport().dynamic_viscosity_pa_s();
        let radius = 0.002;
        let area = core::f64::consts::PI * radius * radius;
        for rv in [0.001_f64, 0.01_f64] {
            let omega = rv * rv * mu / (rho * radius * radius);
            let wave =
                fluid::segment_wave_with_loss(material, radius, omega, LossModel::Bessel).unwrap();
            let series = C64::new(0.0, -1.0) * wave.wavenumber * wave.characteristic_impedance;
            let shunt = C64::new(0.0, -1.0) * wave.wavenumber / wave.characteristic_impedance;
            let resistance = 8.0 * mu / (area * radius * radius);
            let inertia = (4.0 / 3.0) * omega * rho / area;
            let compliance = omega * area * state.heat_capacity_ratio() / (rho * c * c);
            assert!((series.re / resistance - 1.0).abs() < 2e-8);
            assert!((series.im / -inertia - 1.0).abs() < 3e-7);
            assert!((shunt.im / -compliance - 1.0).abs() < 2e-8);
            assert!(shunt.re >= 0.0);
            let duct = cylinder(radius, 0.05);
            let response = lossy(
                &duct,
                material,
                omega,
                LossModel::Bessel,
                Some(C64::ZERO),
                WallModel::Rigid,
            );
            assert!((response.z_in.re / (resistance * 0.05) - 1.0).abs() < 2e-7);
        }
    }
}

#[test]
fn low_prandtl_refuses_wide_thermal_layer_and_zero_thermal_coupling_remains_exact() {
    let acoustic = FluidAcousticState::try_new(7300.0, 2000.0).unwrap();
    let transport = FluidTransportState::try_new(7300.0, 250.0, 100.0, 0.002).unwrap();
    let coupled = ThermoviscousFluidState::try_new(acoustic, transport, 200.0).unwrap();
    let decoupled = ThermoviscousFluidState::try_new(acoustic, transport, 250.0).unwrap();
    let radius = 0.001;
    let omega = 30.0 * 30.0 * transport.dynamic_viscosity_pa_s()
        / (transport.density_kg_m3() * radius * radius);
    assert!(matches!(
        fluid::segment_wave_with_loss(&coupled, radius, omega, LossModel::WideTube),
        Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
    ));
    assert!(fluid::segment_wave_with_loss(&coupled, radius, omega, LossModel::Bessel).is_ok());
    assert!(fluid::segment_wave_with_loss(&decoupled, radius, omega, LossModel::WideTube).is_ok());
    let wave = fluid::segment_wave_with_loss(&decoupled, radius, omega, LossModel::Bessel).unwrap();
    let shunt = C64::new(0.0, -1.0) * wave.wavenumber / wave.characteristic_impedance;
    let compliance =
        omega * core::f64::consts::PI * radius * radius / acoustic.adiabatic_bulk_modulus_pa();
    near(shunt, C64::new(0.0, -compliance), compliance, 3e-13);
}

#[test]
fn original_gas_and_generic_fluid_share_the_same_thermal_validity_boundary() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let radius = 0.01;
    let omega = 11.0 * 11.0 * gas.dynamic_viscosity / (gas.density * radius * radius);
    // The old gas-only path admitted this state solely because rv > 10.
    // Its thermal layer has rt < 10 and cannot use that asymptotic model.
    assert!(11.0 * gas.prandtl.sqrt() < fs_duct::MIN_SHEAR_NUMBER);
    assert!(matches!(
        fs_duct::segment_wave(&gas, radius, omega, LossModel::WideTube),
        Err(DuctError::BadParameter { what }) if what.contains("thermal shear")
    ));
    assert!(matches!(
        fluid::segment_wave_with_loss(&gas, radius, omega, LossModel::WideTube),
        Err(FluidDuctError::Duct(DuctError::BadParameter { what })) if what.contains("thermal shear")
    ));
    let native = fs_duct::segment_wave(&gas, radius, omega, LossModel::Bessel).unwrap();
    let shared = fluid::segment_wave_with_loss(&gas, radius, omega, LossModel::Bessel).unwrap();
    near(
        shared.wavenumber,
        native.wavenumber,
        native.wavenumber.abs(),
        3e-13,
    );
    near(
        shared.characteristic_impedance,
        native.characteristic_impedance,
        native.characteristic_impedance.abs(),
        3e-13,
    );
}

#[test]
fn all_loss_models_retain_gas_pressure_flow_across_cylinders_cones_profiles_and_walls() {
    let air = GasState::try_new_moist_air(310.0, 150_000.0, 0.5).unwrap();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.1,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.016,
                length: 0.2,
            },
            Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(0.016, 1.0, 0.2).unwrap()),
                steps: RiccatiSteps(24),
            },
        ],
    };
    let wall = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 4e7,
        resistance: 300.0,
    };
    for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
        for wall in [WallModel::Rigid, WallModel::Pin(&wall)] {
            let omega = 1700.0;
            let load = C64::new(1e6, -2e5);
            let old = fs_duct::input_impedance_load_with_wall(&duct, &air, omega, loss, load, wall)
                .unwrap();
            let new = lossy(&duct, &air, omega, loss, Some(load), wall);
            near(new.z_in, old.impedance, 1.0, 4e-11);
            near(new.u_out_over_u_in, old.u_mouth_over_u_in, 1.0, 4e-11);
            near(new.p_out_over_u_in, old.p_mouth_over_u_in, 1.0, 4e-11);
            let chain = fluid::duct_chain_with_loss(
                &duct,
                &air,
                omega,
                loss,
                wall,
                ConeSlices::Auto,
                options(),
                || false,
            )
            .unwrap();
            assert_eq!(
                chain.min_shear_number.is_finite(),
                loss != LossModel::Lossless
            );
        }
    }
}

#[test]
fn lossy_material_admission_cancellation_and_work_refusals_publish_no_partial_chain() {
    let (_, liquid) = source_liquid();
    let original = liquid.clone();
    let duct = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.001,
            outlet_radius: 0.004,
            length: 0.1,
        }],
    };
    let mut polls = 0;
    assert!(matches!(
        fluid::duct_chain_with_loss(
            &duct,
            &liquid,
            500.0,
            LossModel::Bessel,
            WallModel::Rigid,
            ConeSlices::Fixed(16),
            options(),
            || {
                polls += 1;
                polls == 5
            }
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let capped = ComposeOptions {
        maximum_blocks: 3,
        ..options()
    };
    assert!(
        fluid::duct_chain_with_loss(
            &duct,
            &liquid,
            500.0,
            LossModel::Bessel,
            WallModel::Rigid,
            ConeSlices::Fixed(16),
            capped,
            || false
        )
        .is_err()
    );
    let mut stale = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    stale.specific_heat_cp *= 2.0;
    assert!(matches!(
        fluid::segment_wave_with_loss(&stale, 0.01, 1000.0, LossModel::Bessel),
        Err(FluidDuctError::Material(_))
    ));
    let retry = || {
        fluid::duct_chain_with_loss(
            &duct,
            &liquid,
            500.0,
            LossModel::Bessel,
            WallModel::Rigid,
            ConeSlices::Fixed(16),
            options(),
            || false,
        )
        .unwrap()
    };
    let first =
        compose::compose_chain_controlled(&retry().blocks, Some(C64::ZERO), options(), &mut || {
            false
        })
        .unwrap();
    let second =
        compose::compose_chain_controlled(&retry().blocks, Some(C64::ZERO), options(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(first.z_in, second.z_in);
    assert_eq!(first.u_out_over_u_in, second.u_out_over_u_in);
    assert_eq!(liquid, original);
}

#[path = "thermoviscous/porous.rs"]
mod porous;
