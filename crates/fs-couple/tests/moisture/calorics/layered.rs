//! Unlike-material heat storage, native curved inventories and one shared gas.
//! Constitutive sources are authored numerical verification data.
#[path = "layered/sdirk.rs"]
mod sdirk;

use super::*;
use fs_couple::moisture::{
    layered::{MoistureActivityTransfer, MoistureLayer, MoistureLayeredSlab},
    radial::{
        MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry,
        layered::{MoistureLayeredRadial, MoistureRadialLayer},
    },
};
use fs_material::state_point::MonotoneMoistureSorption;

const SLOPE: [f64; 2] = [0.2, 0.3];
const LENGTH: [f64; 2] = [0.2, 0.1];
const RHO: f64 = 500.0;
const AREA: f64 = 0.001;
const DIFFUSION: f64 = 0.01;
const CONTACT: f64 = 0.005; // Whole slab interface conductance [kg/s].
const CONSTANT_Q: [[f64; 2]; 2] = [[2e6; 2], [3e6; 2]];
const LOADING_Q: [[f64; 2]; 2] = [[2e6, 2.6e6], [3.2e6, 2.3e6]];

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn material(slope: f64, q: [f64; 2]) -> MaterialCard {
    let original = source(q);
    let mut claims = original.claims().clone();
    let mut sorption = original.claims_for("moisture_diffusivity")[0].1.clone();
    sorption.key = PropertyKey::with_quantity(
        "sorption",
        qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
    );
    sorption.value = PropertyValue::Curve {
        abscissa: "relative_humidity".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(0.0, 0.0), (1.0, slope)],
        dims: Dims::NONE,
    };
    sorption.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
            293.15,
        )
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.0,
            1.0,
        );
    sorption.provenance.source = "authored layered caloric isotherm".into();
    claims.insert_claim(sorption).unwrap();
    MaterialCard::assemble(original.id().clone(), claims, original.models().to_vec()).unwrap()
}

fn isotherm(card: &MaterialCard) -> MonotoneMoistureSorption {
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    let binding = resolve_material_state_point(
        card,
        &point,
        &[ScalarPropertyRequirement::try_with_quantity(
            "sorption",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MonotoneMoistureSorption::from_binding(card, &binding, "sorption", [0.0, 1.0]).unwrap()
}

fn stack(
    q: [[f64; 2]; 2],
    activities: [Vec<f64>; 2],
    cx: &Cx<'_>,
) -> (MoistureLayeredSlab, [MoistureSorptionEnergy; 2]) {
    let cards = std::array::from_fn::<_, 2, _>(|i| material(SLOPE[i], q[i]));
    let laws = std::array::from_fn(|i| law(&cards[i], [0.0, 0.3]));
    let layers = cards
        .into_iter()
        .zip(activities)
        .enumerate()
        .map(|(i, (card, activity))| {
            let sorption = isotherm(&card);
            let field = MoistureSlab::new_kirchhoff(
                card,
                point(0.1),
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                MoistureSlabGeometry {
                    length_m: LENGTH[i],
                    cross_section_m2: AREA,
                    maximum_moisture: 0.3,
                    maximum_cells: 16,
                },
                activity.into_iter().map(|a| SLOPE[i] * a).collect(),
                nonlinear(),
                cx,
            )
            .unwrap();
            MoistureLayer::new(field, sorption, cx).unwrap()
        })
        .collect();
    let body = MoistureLayeredSlab::new(layers, 32, cx)
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::IdealContact; 2],
            vec![MoistureActivityTransfer::ConductanceKgM2S(CONTACT / AREA)],
            cx,
        )
        .unwrap();
    (body, laws)
}

fn coated(
    symmetry: RadialMoistureSymmetry,
    cx: &Cx<'_>,
) -> (MoistureLayeredRadial, [MoistureSorptionEnergy; 2]) {
    let cards = std::array::from_fn::<_, 2, _>(|i| material(SLOPE[i], LOADING_Q[i]));
    let laws = std::array::from_fn(|i| law(&cards[i], [0.0, 0.3]));
    let layers = cards
        .into_iter()
        .enumerate()
        .map(|(i, card)| {
            let sorption = isotherm(&card);
            let field = MoistureRadial::new(
                card,
                point(0.1),
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                MoistureRadialGeometry {
                    symmetry,
                    faces_m: if i == 0 {
                        vec![0.0, 0.004, 0.01]
                    } else {
                        vec![0.01, 0.02, 0.05]
                    },
                    maximum_moisture: 0.3,
                    maximum_cells: 16,
                },
                if i == 0 {
                    vec![0.8 * SLOPE[i], 0.3 * SLOPE[i]]
                } else {
                    vec![0.2 * SLOPE[i], 0.5 * SLOPE[i]]
                },
                cx,
            )
            .unwrap();
            MoistureRadialLayer::new(field, sorption, cx).unwrap()
        })
        .collect();
    let body = MoistureLayeredRadial::new(layers, 32, cx)
        .unwrap()
        .with_transfer(
            [
                MoistureActivityTransfer::IdealContact,
                MoistureActivityTransfer::ConductanceKgM2S(0.04),
            ],
            vec![MoistureActivityTransfer::ConductanceKgM2S(0.03)],
            cx,
        )
        .unwrap();
    (body, laws)
}

fn check_layer(
    report: &SorptionLayerEnergy,
    before: &[f64],
    after: &[f64],
    dry: &[f64],
    q: [f64; 2],
    rt: f64,
) -> f64 {
    assert_eq!(report.cell_binding_energy_change_j.len(), before.len());
    assert_eq!(
        report.cell_bound_internal_energy_change_j.len(),
        before.len()
    );
    let mut expected_heat = 0.0;
    let mut water = 0.0;
    for (i, ((&a, &b), &mass)) in before.iter().zip(after).zip(dry).enumerate() {
        // Independent primitive for q(M)=q0+(q1-q0)*M/0.3.
        let dm = mass * (b - a);
        let heat = dm * (rt - q[0] - (q[1] - q[0]) * (a + b) / 0.6);
        near(report.cell_binding_energy_change_j[i], heat, 3e-9);
        expected_heat += heat;
        water += dm;
    }
    near(report.water_change_kg, water, 2e-14);
    near(report.binding_energy_change_j, expected_heat, 3e-9);
    near(
        report.cell_bound_internal_energy_change_j.iter().sum(),
        report.bound_internal_energy_change_j,
        3e-9,
    );
    expected_heat
}

fn check_whole(energy: &LayeredSorptionEnergyStep, expected_heat: f64, gas_change: f64, dt: f64) {
    assert_eq!(energy.interval_s, [0.0, dt]);
    assert_eq!(energy.gas_internal_energy_change_j, gas_change);
    near(energy.thermostat_heat_j, expected_heat, 8e-9);
    near(energy.mean_thermostat_power_w, expected_heat / dt, 8e-7);
    near(
        energy.bound_internal_energy_change_j,
        energy
            .layers
            .iter()
            .map(|layer| layer.bound_internal_energy_change_j)
            .sum(),
        8e-9,
    );
    near(
        energy.bound_internal_energy_change_j + gas_change - energy.thermostat_heat_j,
        0.0,
        energy_controls().energy_balance_tolerance_j,
    );
    assert!(energy.energy_balance_residual_j.abs() <= energy_controls().energy_balance_tolerance_j);
    assert_eq!(
        energy.source_segments,
        energy
            .layers
            .iter()
            .map(|layer| layer.source_segments)
            .sum::<usize>()
    );
}

#[test]
fn sealed_unlike_layers_release_binding_difference_with_zero_gas_exchange() {
    with_cx(|cx| {
        let capacity = [
            RHO * AREA * LENGTH[0] * SLOPE[0],
            RHO * AREA * LENGTH[1] * SLOPE[1],
        ];
        let conductance = 1.0
            / (0.5 * LENGTH[0] / (RHO * DIFFUSION * SLOPE[0] * AREA)
                + 1.0 / CONTACT
                + 0.5 * LENGTH[1] / (RHO * DIFFUSION * SLOPE[1] * AREA));
        let dt = 0.5;
        for initial in [[0.8, 0.15], [0.15, 0.8]] {
            let (old, laws) = stack(CONSTANT_Q, initial.map(|a| vec![a]), cx);
            let atmosphere = air::atmosphere(0.5, 1.0);
            let result = atmosphere
                .advance_layered_with_sorption_energy(
                    &old,
                    dt,
                    [false; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    &laws,
                    energy_controls(),
                    cx,
                )
                .unwrap();
            // Independent backward-Euler solution of two original water stores.
            let exchange = dt * conductance * (initial[0] - initial[1])
                / (1.0 + dt * conductance * (1.0 / capacity[0] + 1.0 / capacity[1]));
            let expected = [
                initial[0] - exchange / capacity[0],
                initial[1] + exchange / capacity[1],
            ];
            for (i, (&actual, target)) in result
                .native()
                .moisture()
                .activity()
                .iter()
                .zip(expected)
                .enumerate()
            {
                near(actual, target, 2e-12);
                near(
                    result.energy().layers[i].water_change_kg,
                    if i == 0 { -exchange } else { exchange },
                    3e-14,
                );
                assert_eq!(
                    result.energy().layers[i].sorption_energy_source,
                    laws[i].identity()
                );
            }
            near(
                result.native().step().transport.interfaces[0].water_left_to_right_kg,
                exchange,
                3e-14,
            );
            let heat = (CONSTANT_Q[0][0] - CONSTANT_Q[1][0]) * exchange;
            assert!(heat.abs() > 500.0);
            check_whole(
                result.energy(),
                heat,
                result.native().step().gas_internal_energy_change_j,
                dt,
            );
            assert_eq!(result.native().atmosphere(), &atmosphere);
            near(result.native().step().water_to_solid_kg, 0.0, 2e-14);
            near(
                result.native().moisture().water_mass_kg(),
                old.water_mass_kg(),
                2e-14,
            );
            assert_eq!(result.energy().source_segments, 2);
            assert_eq!(old.time_s(), 0.0);
        }
    });
}

#[test]
fn finite_air_layered_uptake_and_desorption_integrate_each_loading_dependent_source() {
    with_cx(|cx| {
        let (old, laws) = stack(LOADING_Q, [vec![0.2, 0.5], vec![0.1, 0.3, 0.6]], cx);
        for rh in [0.85, 0.05] {
            let atmosphere = air::atmosphere(rh, 1.0);
            let dt = 0.5;
            let expected = atmosphere
                .advance_layered(
                    &old,
                    dt,
                    [true; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    cx,
                )
                .unwrap();
            let result = atmosphere
                .advance_layered_with_sorption_energy(
                    &old,
                    dt,
                    [true; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    &laws,
                    energy_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                result.native().moisture().activity(),
                expected.moisture().activity()
            );
            assert_eq!(result.native().atmosphere(), expected.atmosphere());
            let mut expected_heat = 0.0;
            for i in 0..2 {
                let before = old.layers()[i].slab();
                let after = result.native().moisture().layers()[i].slab();
                let dry = vec![
                    before.dry_mass_kg() / before.moisture().len() as f64;
                    before.moisture().len()
                ];
                expected_heat += check_layer(
                    &result.energy().layers[i],
                    before.moisture(),
                    after.moisture(),
                    &dry,
                    LOADING_Q[i],
                    rt(&atmosphere),
                );
                assert_eq!(
                    result.energy().layers[i].sorption_energy_source,
                    laws[i].identity()
                );
            }
            check_whole(
                result.energy(),
                expected_heat,
                expected.step().gas_internal_energy_change_j,
                dt,
            );
            let change = result.native().step().water_to_solid_kg;
            assert_eq!(change.signum(), (rh - 0.5_f64).signum());
            near(
                result
                    .energy()
                    .layers
                    .iter()
                    .map(|x| x.water_change_kg)
                    .sum(),
                change,
                2e-13,
            );
            near(
                result.native().atmosphere().vapor_mass_kg().unwrap()
                    + result.native().moisture().water_mass_kg(),
                atmosphere.vapor_mass_kg().unwrap() + old.water_mass_kg(),
                2e-12,
            );
            assert_eq!(result.energy().source_segments, 5);
        }
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn coated_spheres_and_cylinders_use_original_unequal_shells_and_one_gas_datum() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let (old, laws) = coated(symmetry, cx);
            assert!(
                old.layers()[1].radial().dry_masses_kg()[1]
                    > 40.0 * old.layers()[0].radial().dry_masses_kg()[0]
            );
            for (rh, exposed) in [
                (0.8, [false, true]),
                (0.05, [false, true]),
                (0.5, [false; 2]),
            ] {
                let atmosphere = air::atmosphere(rh, 1.0);
                let dt = 0.002;
                let expected = atmosphere
                    .advance_layered_radial(
                        &old,
                        dt,
                        exposed,
                        config(),
                        nonlinear(),
                        air::controls(),
                        cx,
                    )
                    .unwrap();
                let result = atmosphere
                    .advance_layered_radial_with_sorption_energy(
                        &old,
                        dt,
                        exposed,
                        config(),
                        nonlinear(),
                        air::controls(),
                        &laws,
                        energy_controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(
                    result.native().moisture().activity(),
                    expected.moisture().activity()
                );
                assert_eq!(result.native().atmosphere(), expected.atmosphere());
                let mut expected_heat = 0.0;
                for i in 0..2 {
                    let before = old.layers()[i].radial();
                    let after = result.native().moisture().layers()[i].radial();
                    expected_heat += check_layer(
                        &result.energy().layers[i],
                        before.moisture(),
                        after.moisture(),
                        before.dry_masses_kg(),
                        LOADING_Q[i],
                        rt(&atmosphere),
                    );
                    assert_eq!(
                        result.energy().layers[i].sorption_energy_source,
                        laws[i].identity()
                    );
                }
                check_whole(
                    result.energy(),
                    expected_heat,
                    expected.step().gas_internal_energy_change_j,
                    dt,
                );
                assert_eq!(result.energy().source_segments, 4);
                if exposed == [false; 2] {
                    assert_eq!(result.native().atmosphere(), &atmosphere);
                    near(result.native().step().water_to_solid_kg, 0.0, 2e-13);
                    assert!(expected_heat.abs() > 1.0);
                }
                assert_eq!(
                    result.energy().vapor_calorics.species(),
                    atmosphere.source_catalog().water_species().unwrap()
                );
            }
            assert_eq!(old.time_s(), 0.0);
        }
    });
}

#[test]
fn layered_partition_equilibrium_needs_no_source_integration_work() {
    with_cx(|cx| {
        let (old, laws) = stack(CONSTANT_Q, [vec![0.5; 2], vec![0.5; 3]], cx);
        assert_ne!(
            old.layers()[0].slab().moisture()[0],
            old.layers()[1].slab().moisture()[0]
        );
        let atmosphere = air::atmosphere(0.5, 1.0);
        let controls = SorptionEnergyControls {
            maximum_source_segments: 0,
            ..energy_controls()
        };
        let result = atmosphere
            .advance_layered_with_sorption_energy(
                &old,
                1.0,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                &laws,
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(result.native().moisture().activity(), old.activity());
        assert_eq!(result.energy().source_segments, 0);
        assert_eq!(result.energy().thermostat_heat_j, 0.0);
        for layer in &result.energy().layers {
            assert_eq!(layer.water_change_kg, 0.0);
            assert_eq!(layer.binding_energy_change_j, 0.0);
        }
    });
}

#[test]
fn layered_caloric_refusals_leave_every_inventory_and_clock_unchanged() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (old, laws) = stack(CONSTANT_Q, [vec![0.8], vec![0.15]], &cx);
        let saved_activity = old.activity().to_vec();
        let atmosphere = air::atmosphere(0.5, 1.0);
        let saved_air = atmosphere.clone();
        let advance = |laws: &[MoistureSorptionEnergy], controls| {
            atmosphere.advance_layered_with_sorption_energy(
                &old,
                0.5,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                laws,
                controls,
                &cx,
            )
        };
        assert!(matches!(
            advance(&[], energy_controls()),
            Err(SorptionEnergyError::Invalid(_))
        ));
        let reversed = [laws[1].clone(), laws[0].clone()];
        assert!(matches!(
            advance(&reversed, energy_controls()),
            Err(SorptionEnergyError::Invalid(_))
        ));
        let narrow = [
            laws[0].clone(),
            law(&material(SLOPE[1], CONSTANT_Q[1]), [0.0, 0.045000000001]),
        ];
        assert!(matches!(
            advance(&narrow, energy_controls()),
            Err(SorptionEnergyError::Moisture(MoistureError::Material(_)))
        ));
        assert!(matches!(
            advance(
                &laws,
                SorptionEnergyControls {
                    maximum_source_segments: 1,
                    ..energy_controls()
                }
            ),
            Err(SorptionEnergyError::SourceBudget {
                required: 2,
                maximum: 1
            })
        ));
        // The first material's full integral fits; the shared cap refuses on
        // the second material rather than silently renewing the allowance.
        let valid = advance(&laws, energy_controls()).unwrap();
        assert_eq!(valid.energy().source_segments, 2);
        gate.request();
        assert!(matches!(
            advance(&laws, energy_controls()),
            Err(SorptionEnergyError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(atmosphere, saved_air);
        assert_eq!(old.activity(), saved_activity);
        assert_eq!(old.time_s(), 0.0);
    });
}
