//! Original phase emissivity closes native heat and survives for observation.

use super::{input, source_fixture, thermal_source_fixture};
use fs_conduction::lumped::LumpedPhaseThermalTransport;
use fs_conduction::lumped::endpoint::LumpedPhaseThermalEnvironment;
use fs_conduction::radiation::{MaterialGraySurfaceLaw, STEFAN_BOLTZMANN_W_M2_K4};
use fs_couple::incandescent::PreparedIncandescentSource;
use fs_couple::incandescent::history::{
    IncandescentDriveHistory, IncandescentHistoryAccuracy, IncandescentHistoryControls,
    advance_history_adaptive,
};
use fs_matdb::{ClaimSelection, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::conductor::{ElectricalPropertyBasis, OhmicDrive, graph::MaterialOhmicNode};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QtyAny};
use fs_scenario::TimeSignal;

fn laws<'a>(
    node: &'a MaterialOhmicNode,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
) -> (PreparedPhaseConductivity<'a>, MaterialGraySurfaceLaw<'a>) {
    let point = ThermalMaterialPoint {
        material: node.source().card(),
        context,
        temperature_axis: "T",
        liquid_fraction_axis: Some("f"),
    };
    let bounds = [
        curve.knots()[0].specific_enthalpy_j_kg,
        curve.knots().last().unwrap().specific_enthalpy_j_kg,
    ];
    let conductivity = PreparedPhaseConductivity::new(
        curve,
        point,
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap();
    let radiation = MaterialGraySurfaceLaw::new(
        "filament",
        curve,
        point,
        MaterialPropertySelection::SingleClaimOnly,
        bounds,
    )
    .unwrap();
    (conductivity, radiation)
}

fn environment<'a>(
    transport: &'a LumpedPhaseThermalTransport<'a>,
) -> LumpedPhaseThermalEnvironment<'a> {
    LumpedPhaseThermalEnvironment {
        temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        convection_w_per_m2_k: 0.0,
        transport,
        enthalpy_tolerance_j_kg: 1e-9,
        maximum_thermal_residual_j: 1e-9,
    }
}

#[test]
fn g1_phase_radiation_and_electrical_work_close_at_the_same_latent_state() {
    let (node, curve) = source_fixture(ElectricalPropertyBasis::Resistivity, true);
    let context = node.source().original().material().query();
    let (conductivity, radiation) = laws(&node, &curve, &context);
    let transport = LumpedPhaseThermalTransport::new(&conductivity, &radiation, || Ok(())).unwrap();
    let bath = environment(&transport);
    let source = PreparedIncandescentSource::try_new_with_caloric_curve(
        &node,
        OhmicDrive::Current(2.0),
        &curve,
    )
    .unwrap();
    let proposal = input(&curve, 4000.0, 0.1);
    let black_power = STEFAN_BOLTZMANN_W_M2_K4
        * proposal.surface_area_m2
        * (600.0_f64.powi(4) - 300.0_f64.powi(4));
    // At h=4000: epsilon=0.5, Joule power=8 W. Along the plateau,
    // dP_J/dh=0.004 and dP_rad,in/dh=0.0003*black_power.
    let slope = 0.004 + 0.0003 * black_power;
    let limit = source.maximum_step_s(&bath, proposal, || Ok(())).unwrap();
    let exact_limit = proposal.mass_kg / (2.0 * slope);
    assert!(limit <= exact_limit && limit > exact_limit * (1.0 - 1e-12));
    let delta_h = proposal.duration_s * (8.0 - 0.5 * black_power)
        / (proposal.mass_kg - proposal.duration_s * slope);
    let step = source.advance(&bath, proposal, || Ok(())).unwrap();
    assert_eq!(step.thermal.state.temperature_k(), 600.0);
    assert!((step.thermal.state.specific_enthalpy_j_kg() - 4000.0 - delta_h).abs() < 1e-8);
    let observed = step.thermal_surface.as_ref().unwrap();
    let epsilon = 0.8 - 0.6 * step.thermal.state.liquid_mass_fraction();
    assert_eq!(observed.phase, step.thermal.state);
    assert_eq!(
        observed.material.query().axes()["f"],
        step.thermal.state.liquid_mass_fraction()
    );
    assert!((observed.emissivity.value() - epsilon).abs() < 1e-14);
    assert!(
        (step.thermal.endpoint.radiation_heat_j + proposal.duration_s * epsilon * black_power)
            .abs()
            < 1e-9
    );
    assert_eq!(
        step.thermal.additional_internal_heat_j,
        step.electrical.joule_power_w * proposal.duration_s
    );
    assert!(step.thermal.endpoint.residual_j.abs() <= 1e-9);
    assert!(
        source
            .advance(&bath, input(&curve, 4000.0, limit * 1.01), || Ok(()))
            .is_err()
    );

    let signal = TimeSignal::Constant(QtyAny::new(2.0, Dims([0, 0, 0, 0, 1, 0])));
    let history = advance_history_adaptive(
        &node,
        &bath,
        proposal,
        IncandescentDriveHistory::Current(&signal),
        IncandescentHistoryControls {
            start_s: 0.0,
            maximum_step_s: 0.1,
            maximum_steps: 100,
        },
        IncandescentHistoryAccuracy {
            maximum_local_enthalpy_error_j_kg: 0.05,
        },
        || Ok(()),
    )
    .unwrap();
    let heat: f64 = history
        .steps
        .iter()
        .map(|row| {
            assert!(row.step.thermal_surface.is_some());
            row.step.thermal.additional_internal_heat_j + row.step.thermal.endpoint.radiation_heat_j
        })
        .sum();
    let final_h = history
        .steps
        .last()
        .unwrap()
        .step
        .thermal
        .state
        .specific_enthalpy_j_kg();
    assert!((proposal.mass_kg * (final_h - 4000.0) - heat).abs() < 1e-8);
    assert!(history.maximum_accepted_error_j_kg <= 0.05);
}

#[test]
fn g4_liquid_conductivity_and_actual_source_conditions_remain_load_bearing() {
    for liquid_k in [1e6, 0.1] {
        let (node, curve) =
            thermal_source_fixture(ElectricalPropertyBasis::Resistivity, true, liquid_k);
        for pressure in [1e5, 2e5] {
            let context = node
                .source()
                .original()
                .material()
                .query()
                .with("pressure", pressure)
                .unwrap();
            let (conductivity, radiation) = laws(&node, &curve, &context);
            let transport =
                LumpedPhaseThermalTransport::new(&conductivity, &radiation, || Ok(())).unwrap();
            let source = PreparedIncandescentSource::try_new_with_caloric_curve(
                &node,
                OhmicDrive::Current(0.0),
                &curve,
            )
            .unwrap();
            // Initially solid still needs whole-chart conductivity admission.
            let accepted = source.advance(
                &environment(&transport),
                input(&curve, 0.0, 0.01),
                || Ok(()),
            );
            assert_eq!(accepted.is_ok(), liquid_k == 1e6 && pressure == 1e5);
        }
    }
}
