//! Actual coexistence compositions drive existing flow and component exchange.

use super::*;
use fs_flux::capillary::material::miscibility::{
    step_binary_liquid_horizontal, step_binary_liquid_planar_horizontal,
};
use fs_material::liquid::resolve_liquid_state;
use fs_material::liquid::solution::gradient::BinaryLiquidSquareGradientLaw;
use fs_material::liquid::solution::gradient::interface::ResolvedBinaryLiquidPlanarInterface;
use fs_material::liquid::solution::{BinaryLiquidSolutionLaw, ResolvedBinaryLiquidEquilibrium};
use fs_material::state_point::{
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{
    CompositionBasis, FrequencyConvention, QuantityKind, SemanticType, ValueForm,
};
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;
use fs_thermochem::solution::gradient::interface::BinaryPlanarInterfaceOptionsV1;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn equilibrium(temperature: f64) -> ResolvedBinaryLiquidEquilibrium {
    equilibrium_with_source_context(temperature, None)
}

fn equilibrium_with_source_context(
    temperature: f64,
    context: Option<(&str, QuantitySpec, bool)>,
) -> ResolvedBinaryLiquidEquilibrium {
    with_liquid_law(context, [800.0, 1200.0], |law| {
        law.equilibrate(
            Temperature::new(temperature),
            [4.0, 6.0],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap()
    })
}

fn with_liquid_law<T>(
    context: Option<(&str, QuantitySpec, bool)>,
    density: [f64; 2],
    use_law: impl FnOnce(&BinaryLiquidSolutionLaw<'_>) -> T,
) -> T {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [2000.0, 3000.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            [0.2, 0.6],
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [0.001, 0.005],
        ),
    ] {
        let mut source = claim(name, dims, values[0]);
        if name == THERMAL_CONDUCTIVITY_PROPERTY {
            source.key =
                PropertyKey::with_quantity(name, quantity(QuantityKind::ThermalConductivity));
        }
        source.value = PropertyValue::Curve {
            abscissa: "x_a".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, values[0]), (1.0, values[1])],
            dims,
        };
        source.interpolation = InterpolationPolicy::LinearInside;
        source.validity = ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
                600.0,
            )
            .with_quantity(
                "x_a",
                quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
                0.0,
                1.0,
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                100_000.0,
                100_000.0,
            );
        if let Some((property, response_quantity, response_curve)) = context
            && property == name
        {
            source.validity =
                source
                    .validity
                    .with_quantity("response_context", response_quantity, 0.0, 10.0);
            if response_curve {
                source.value = PropertyValue::Curve {
                    abscissa: "response_context".into(),
                    abscissa_dims: response_quantity.dims(),
                    knots: vec![(0.0, values[0]), (10.0, values[1])],
                    dims,
                };
            }
        }
        claims.insert_claim(source).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic binary capillary solution".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity(
            "x_a",
            quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
            0.4,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    if let Some((_, response_quantity, _)) = context {
        point = point
            .with_quantity("response_context", response_quantity, 1.0)
            .unwrap();
    }
    let reference =
        resolve_liquid_state(&card, &point, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let model = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        900.0,
        Temperature::new(300.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let law = BinaryLiquidSolutionLaw::try_new(
        model,
        [MolarMass::new(0.018), MolarMass::new(0.046)],
        &card,
        reference,
        "T",
        "x_a",
    )
    .unwrap();
    use_law(&law)
}

fn planar_interface() -> ResolvedBinaryLiquidPlanarInterface {
    // rho = c * (x_A M_A + x_B M_B), matching the same fixed site volume.
    with_liquid_law(None, [10_000.0 * 0.046, 10_000.0 * 0.018], |law| {
        let gradient = BinarySquareGradientV1::new(
            law.thermodynamics().clone(),
            Temperature::new(300.0),
            AmountConcentration::new(10_000.0),
            1e-9,
        )
        .unwrap();
        BinaryLiquidSquareGradientLaw::try_new(law, gradient, 1e-12)
            .unwrap()
            .equilibrium_interface(
                [4.0, 6.0],
                BinaryPlanarInterfaceOptionsV1::default(),
                || false,
            )
            .unwrap()
            .unwrap()
    })
}

#[test]
fn derived_planar_tension_preserves_the_original_flow_sources_and_inventory() {
    let planar = planar_interface();
    let before = planar.clone();
    let sigma = planar.interface().surface_tension_n_m();
    assert!(sigma > 0.0);
    for invading in [0, 1] {
        for pressure in [1000.0, -1000.0] {
            let input = MaterialMeniscusInput {
                wet_length_m: 0.03,
                applied_pressure_pa: pressure,
                ..input()
            };
            let step =
                step_binary_liquid_planar_horizontal(&planar, invading, 0.4, 0.01, input).unwrap();
            let explicit =
                step_binary_liquid_horizontal(planar.equilibrium(), invading, sigma, 0.4, input)
                    .unwrap();
            assert_eq!(step.meniscus, explicit);
            assert!(std::ptr::eq(step.planar_interface, &planar));
            assert!(std::ptr::eq(
                step.meniscus.equilibrium,
                planar.equilibrium()
            ));
            assert_eq!(step.maximum_interface_width_radius_ratio, 0.01);
            close(
                step.interface_width_radius_ratio,
                planar.interface().central_slope_width_m() / input.radius_m,
            );
            assert!(step.interface_width_radius_ratio < 0.01);
            assert_eq!(
                step.meniscus
                    .discrete
                    .meniscus
                    .transferred_volume_m3
                    .signum(),
                pressure.signum(),
            );
        }
    }
    assert_eq!(planar, before);
}

#[test]
fn diffuse_interface_and_invalid_wall_angle_refuse_sharp_capillary_flow() {
    let planar = planar_interface();
    let before = planar.clone();
    let ratio = planar.interface().central_slope_width_m() / input().radius_m;
    for allowance in [0.5 * ratio, 0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            step_binary_liquid_planar_horizontal(
                &planar,
                0,
                0.4,
                allowance,
                MaterialMeniscusInput {
                    duration_s: f64::NAN,
                    ..input()
                },
            ),
            Err(MaterialCapillaryError::Unsupported(_))
        ));
    }
    for cosine in [-1.01, 1.01, f64::NAN] {
        assert!(matches!(
            step_binary_liquid_planar_horizontal(&planar, 0, cosine, 0.01, input()),
            Err(MaterialCapillaryError::Flow(_))
        ));
    }
    assert!(matches!(
        step_binary_liquid_planar_horizontal(&planar, 2, 0.4, 0.01, input()),
        Err(MaterialCapillaryError::Unsupported(_))
    ));
    assert_eq!(planar, before);
}

#[test]
fn coexisting_sources_set_exact_flow_work_and_signed_component_exchange() {
    let equilibrium = equilibrium(300.0);
    let before = equilibrium.clone();
    for invading in [0, 1] {
        for pressure in [1000.0, -1000.0] {
            let input = MaterialMeniscusInput {
                wet_length_m: 0.03,
                applied_pressure_pa: pressure,
                ..input()
            };
            let step =
                step_binary_liquid_horizontal(&equilibrium, invading, 0.02, 0.4, input).unwrap();
            let inlet = &equilibrium.phases()[invading];
            let outlet = &equilibrium.phases()[1 - invading];
            let m = step.discrete.meniscus;
            let distance = m.wet_length_m - input.wet_length_m;
            let mu_in = inlet.liquid().dynamic_viscosity_pa_s();
            let mu_out = outlet.liquid().dynamic_viscosity_pa_s();
            let resistance_integral = distance
                * (mu_in * input.wet_length_m
                    + mu_out * (input.tube_length_m - input.wet_length_m)
                    + 0.5 * (mu_in - mu_out) * distance);
            let drive = pressure + 2.0 * 0.02 * 0.4 / input.radius_m;
            close(
                resistance_integral,
                input.radius_m.powi(2) * drive * m.elapsed_time_s / 8.0,
            );
            close(
                m.pressure_work_j - m.wall_surface_energy_change_j,
                m.viscous_dissipation_j,
            );
            close(
                step.discrete.liquid_viscous_dissipation_j
                    + step.discrete.displaced_viscous_dissipation_j,
                m.viscous_dissipation_j,
            );
            assert!(m.viscous_dissipation_j > 0.0);
            assert_eq!(m.transferred_volume_m3.signum(), pressure.signum());
            assert!(std::ptr::eq(step.equilibrium, &equilibrium));
            for (phase, transfers) in [
                (inlet, step.invading_component_inflow_mol),
                (outlet, step.displaced_component_outflow_mol),
            ] {
                let x = phase.mole_fractions();
                let mean_molar_mass = x[0] * 0.018 + x[1] * 0.046;
                let amount =
                    m.transferred_volume_m3 * phase.liquid().density_kg_m3() / mean_molar_mass;
                for component in 0..2 {
                    close(transfers[component], x[component] * amount);
                    assert_eq!(transfers[component].signum(), pressure.signum());
                }
            }
            close(
                step.invading_mass_inflow_kg,
                step.invading_component_inflow_mol[0] * 0.018
                    + step.invading_component_inflow_mol[1] * 0.046,
            );
            close(
                step.displaced_mass_outflow_kg,
                step.displaced_component_outflow_mol[0] * 0.018
                    + step.displaced_component_outflow_mol[1] * 0.046,
            );
            for component in 0..2 {
                close(
                    step.tube_component_change_mol[component],
                    step.invading_component_inflow_mol[component]
                        - step.displaced_component_outflow_mol[component],
                );
            }
        }
    }
    let end = step_binary_liquid_horizontal(
        &equilibrium,
        0,
        0.02,
        0.4,
        MaterialMeniscusInput {
            wet_length_m: 0.03,
            duration_s: 1e6,
            ..input()
        },
    )
    .unwrap();
    assert_eq!(end.discrete.meniscus.end, CircularMeniscusEnd::Outlet);
    close(
        end.discrete.meniscus.transferred_volume_m3,
        PI * input().radius_m.powi(2) * 0.07,
    );
    assert!(end.discrete.meniscus.elapsed_time_s < 1e6);
    assert_eq!(equilibrium, before);
}

#[test]
fn miscible_solution_refuses_a_meniscus_and_pair_parameters_stay_explicit() {
    let mixed = equilibrium(600.0);
    assert_eq!(mixed.phases().len(), 1);
    // Invalid geometry proves the phase refusal precedes the flow operator.
    assert!(matches!(
        step_binary_liquid_horizontal(
            &mixed,
            0,
            0.02,
            0.4,
            MaterialMeniscusInput {
                radius_m: f64::NAN,
                ..input()
            }
        ),
        Err(MaterialCapillaryError::Unsupported(_))
    ));
    let split = equilibrium(300.0);
    assert!(matches!(
        step_binary_liquid_horizontal(&split, 2, 0.02, 0.4, input()),
        Err(MaterialCapillaryError::Unsupported(_))
    ));
    for (sigma, cosine) in [(-0.01, 0.4), (f64::NAN, 0.4), (0.02, 1.01)] {
        assert!(matches!(
            step_binary_liquid_horizontal(&split, 0, sigma, cosine, input()),
            Err(MaterialCapillaryError::Flow(_))
        ));
    }
    let zero = step_binary_liquid_horizontal(&split, 0, 0.0, 0.4, input()).unwrap();
    assert_eq!(zero.discrete.meniscus.transferred_volume_m3, 0.0);
    assert_eq!(zero.discrete.meniscus.viscous_dissipation_j, 0.0);
    assert_eq!(zero.invading_component_inflow_mol, [0.0; 2]);
    assert_eq!(zero.displaced_component_outflow_mol, [0.0; 2]);
    assert_eq!(zero.tube_component_change_mol, [0.0; 2]);
}

#[test]
fn frozen_snapshots_do_not_infer_newtonian_density_or_viscosity() {
    for property in [DENSITY_PROPERTY, DYNAMIC_VISCOSITY_PROPERTY] {
        for response_quantity in [
            quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            QuantitySpec::dimensional(Frequency::DIMS),
        ] {
            for response_curve in [false, true] {
                let split = equilibrium_with_source_context(
                    300.0,
                    Some((property, response_quantity, response_curve)),
                );
                let before = split.clone();
                assert_eq!(split.phases().len(), 2);
                assert!(matches!(
                    step_binary_liquid_horizontal(
                        &split,
                        0,
                        0.02,
                        0.4,
                        MaterialMeniscusInput {
                            radius_m: f64::NAN,
                            ..input()
                        }
                    ),
                    Err(MaterialCapillaryError::Unsupported(_))
                ));
                assert_eq!(split, before);
            }
        }
    }
}
