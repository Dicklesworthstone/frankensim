//! G1/G3: contact and uncovered baths share both finite-body endpoints.
use super::*;
use fs_conduction::interface::lumped_pair::advance_contact_pair_with_environments;
use fs_conduction::lumped::LumpedThermalEnvironment;

fn environment(temperature_k: f64, h: f64, emissivity: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k,
        radiation_temperature_k: temperature_k,
        convection_w_per_m2_k: h,
        transport: LumpedThermalTransport::try_declared(1000.0, emissivity).unwrap(),
        enthalpy_tolerance_j_kg: 1e-12,
        maximum_thermal_residual_j: 1e-12,
    }
}

fn boundary_controls() -> ContactIterationConfig {
    ContactIterationConfig {
        max_iterations: 128,
        relaxation: 0.5,
        residual_tolerance_j: 1e-13,
    }
}

#[test]
fn g1_two_baths_and_contact_match_the_joint_backward_euler_linear_system() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let mut inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    inputs[1].internal_heat_j = 1e-5;
    let (source, id) = temperature_interface(&cards, 0.0);
    let context = point();
    let interface = InterfaceResistance::from_card_pinned(
        "contact",
        &source,
        &binding(&cards, &source, &context)
            .at(inputs.map(|i| i.initial))
            .unwrap(),
        id,
    )
    .unwrap();
    let baths = [environment(320.0, 10.0, 0.0), environment(300.0, 20.0, 0.0)];
    let area = 2e-4;
    let dt = inputs[0].duration_s;
    let capacities = [2e-5, 4e-5];
    let conductance = area / 1e-4;
    let exposed = inputs.map(|i| i.surface_area_m2 - area);
    let bath_g = [0, 1].map(|side| baths[side].convection_w_per_m2_k * exposed[side]);
    let diagonal = [0, 1].map(|side| capacities[side] + dt * (conductance + bath_g[side]));
    let cross = dt * conductance;
    let rhs = [0, 1].map(|side| {
        capacities[side] * inputs[side].initial.temperature_k()
            + inputs[side].internal_heat_j
            + dt * bath_g[side] * baths[side].temperature_k
    });
    let determinant = diagonal[0] * diagonal[1] - cross * cross;
    let expected = [
        (diagonal[1] * rhs[0] + cross * rhs[1]) / determinant,
        (cross * rhs[0] + diagonal[0] * rhs[1]) / determinant,
    ];
    support::with_cx(|cx| {
        let pair = advance_contact_pair_with_environments(
            cx,
            inputs,
            [&baths[0], &baths[1]],
            &interface,
            area,
            1e-12,
            boundary_controls(),
        )
        .unwrap();
        let report = &pair[0].report;
        assert!(report.iterations > 1);
        assert_eq!(report.exposed_area_m2, exposed);
        for side in 0..2 {
            let temperature = pair[side].state.temperature_k();
            assert!((temperature - expected[side]).abs() < 2e-7);
            let boundary = dt * bath_g[side] * (baths[side].temperature_k - temperature);
            assert!((report.convection_heat_j[side] - boundary).abs() < 1e-17);
            assert_eq!(report.radiation_heat_j[side], 0.0);
            let contact = if side == 0 {
                -report.heat_a_to_b_j
            } else {
                report.heat_a_to_b_j
            };
            assert!((pair[side].external_heat_j - boundary - contact).abs() < 1e-17);
            assert!(report.energy_residual_j[side].abs() <= 1e-12);
        }
        assert!(
            (report.heat_a_to_b_j
                - dt * conductance
                    * (pair[0].state.temperature_k() - pair[1].state.temperature_k()))
            .abs()
                < 1e-12
        );
        assert!(report.total_energy_residual_j.abs() <= 1e-12);
    });
}

#[test]
fn g3_material_contact_radiation_coverage_and_separation_keep_original_owners() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let (source, id) = temperature_interface(&cards, 0.1);
    let context = point();
    let law = MaterialContactLaw::new(
        binding(&cards, &source, &context),
        inputs.map(|i| i.initial),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let baths = [environment(320.0, 10.0, 0.8), environment(305.0, 20.0, 0.4)];
    let area = 2e-4;
    support::with_cx(|cx| {
        let pair = law
            .advance_with_environments(
                cx,
                inputs,
                [&baths[0], &baths[1]],
                area,
                1e-14,
                contact_controls(),
                boundary_controls(),
            )
            .unwrap();
        let report = &pair[0].report;
        assert_eq!(
            report.contact.endpoint.properties()[0]
                .answer()
                .receipt
                .selected,
            id
        );
        let temperatures = pair.each_ref().map(|p| p.state.temperature_k());
        let resistance = 1e-4 * (1.0 + 0.1 * (0.5 * temperatures.iter().sum::<f64>() - 300.0));
        assert!(
            (report.heat_a_to_b_j
                - inputs[0].duration_s * area / resistance * (temperatures[0] - temperatures[1]))
                .abs()
                < 1e-13
        );
        for side in 0..2 {
            let radiation = [0.8, 0.4][side]
                * 5.670_374_419e-8
                * (inputs[side].surface_area_m2 - area)
                * inputs[side].duration_s
                * (baths[side].radiation_temperature_k.powi(4) - temperatures[side].powi(4));
            assert!((report.radiation_heat_j[side] - radiation).abs() < 1e-17);
            assert!(report.energy_residual_j[side].abs() <= 1e-12);
        }
        // A fully covered surface receives none of either bath's heat.
        let covered = law
            .advance_with_environments(
                cx,
                inputs,
                [&baths[0], &baths[1]],
                inputs[0].surface_area_m2,
                1e-14,
                contact_controls(),
                boundary_controls(),
            )
            .unwrap();
        let isolated = law
            .advance(
                cx,
                inputs,
                [&baths[0].transport, &baths[1].transport],
                inputs[0].surface_area_m2,
                1e-14,
                contact_controls(),
            )
            .unwrap();
        for side in 0..2 {
            assert_eq!(covered[side].state, isolated[side].state);
            assert_eq!(
                covered[side].external_heat_j,
                isolated[side].external_heat_j
            );
            assert_eq!(covered[side].report.convection_heat_j, [0.0; 2]);
            assert_eq!(covered[side].report.radiation_heat_j, [0.0; 2]);
        }
        // Separation restores the complete exposed areas and both native baths.
        let open = law
            .advance_with_environments(
                cx,
                inputs,
                [&baths[0], &baths[1]],
                0.0,
                1e-14,
                contact_controls(),
                boundary_controls(),
            )
            .unwrap();
        assert_eq!(open[0].report.heat_a_to_b_j, 0.0);
        assert!(open[0].report.contact.applied.contact.is_none());
        for side in 0..2 {
            let native = baths[side].advance(cx, inputs[side]).unwrap();
            assert!((open[side].state.temperature_k() - native.state.temperature_k()).abs() < 2e-7);
        }
    });
}

#[test]
fn g4_joint_biot_work_and_cancellation_refusals_preserve_the_pair() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let initial = inputs.map(|i| i.initial);
    let (source, id) = temperature_interface(&cards, 0.0);
    let context = point();
    let interface = InterfaceResistance::from_card_pinned(
        "contact",
        &source,
        &binding(&cards, &source, &context).at(initial).unwrap(),
        id,
    )
    .unwrap();
    let baths = [environment(320.0, 10.0, 0.8), environment(305.0, 20.0, 0.4)];
    support::with_cx(|cx| {
        let mut exhausted = boundary_controls();
        exhausted.max_iterations = 1;
        assert!(
            advance_contact_pair_with_environments(
                cx,
                inputs,
                [&baths[0], &baths[1]],
                &interface,
                2e-4,
                1e-12,
                exhausted
            )
            .is_err()
        );
        let mut fast = environment(310.0, 4500.0, 0.0);
        fast.transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
        let short = inputs.map(|input| UniformEnthalpyStepInput {
            duration_s: 1e-8,
            ..input
        });
        // Each path alone admits Bi < 0.1, but their sum must refuse.
        fast.prepare_endpoint(short[0]).unwrap();
        let contact =
            advance_contact_pair(cx, short, [&fast.transport; 2], &interface, 2e-4, 1e-12).unwrap();
        assert!(contact[0].report.contact.as_ref().unwrap().maximum_biot < 0.1);
        assert!(
            advance_contact_pair_with_environments(
                cx,
                short,
                [&fast; 2],
                &interface,
                2e-4,
                1e-12,
                boundary_controls()
            )
            .is_err()
        );
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            advance_contact_pair_with_environments(
                cx,
                inputs,
                [&baths[0], &baths[1]],
                &interface,
                2e-4,
                1e-12,
                boundary_controls()
            ),
            Err(fs_conduction::ConductionError::Cancelled { .. })
        ));
    });
    assert_eq!(inputs.map(|i| i.initial), initial);
}
