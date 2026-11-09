//! Moving deposition, local radiation and original fa/fb contact feedback.
use super::*;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_evidence::ValidityDomain;

#[test]
fn g1_pointwise_scan_keeps_both_latent_contact_fractions_and_external_energy() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        // Two area-1/2 exterior faces: [0,1,3] and [4,5,7]. Their shared
        // contact-edge nodes make radiation change the live interface fractions.
        let boundary = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh)
                .region("cold exterior", |face| face.centroid[1] == 0.0 && face.centroid[2] < 0.0,
                    ThermalBc::neumann(0.0).unwrap()).unwrap()
                .region("hot exterior", |face| face.centroid[1] == 0.0 && face.centroid[2] > 0.0,
                    ThermalBc::neumann(0.0).unwrap()).unwrap()
                .adiabatic_remainder(),
            &Environment::earth_lab(), Vec::new(), 0.0,
        ).unwrap();
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source_card = source(&cards, 1.0, false, false, 1.0);
        let query_point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let contact_laws = contacts(
            cx, &mesh, &initial, &cards, &curves, &source_card, &query_point,
            sides(&mesh, 4, 2, 0.0), 0, false,
        );

        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar { value: 0.2, dims: EMISSIVITY_DIMS },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance { source: "synthetic latent scan radiation".into(),
                license: "internal-test-use".into(), artifact: None },
        }).unwrap();
        let surface = MaterialCard::assemble(MaterialStateId {
            chemistry: "latent scan surface".into(), phase: "declared gray".into(),
            process: "fixture".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let patches = ["cold exterior", "hot exterior"].into_iter().map(|region| {
            RadiationRegionHistory {
                region: region.into(),
                emissivity: SurfaceEmissivity::from_card(region, &surface, 400.0,
                    SelectionPolicy::SingleClaimOnly).unwrap(),
                ambient_temperature: TimeSignal::Constant(QtyAny::new(300.0, TEMPERATURE_DIMS)),
            }
        }).collect();
        let path = TimeSignal::Ramp {
            t_start: 0.0, t_end: 0.02,
            from: QtyAny::new(-0.25, LENGTH), to: QtyAny::new(0.25, LENGTH),
        };
        let mut scan = ScanningEnthalpyRun::new(
            cx, scan_base(cx, &boundary, &conductivity, &curves, 2),
            vec![heater(path)], 1, None,
        ).unwrap()
            .with_contacts_newton(cx, contact_laws, config()).unwrap()
            .with_pointwise_ambient_radiation(cx, patches).unwrap();

        let radiation_w = [400_f64, 500_f64].map(|temperature| {
            0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5
                * (temperature.powi(4) - 300_f64.powi(4))
        });
        let mut old = OLD;
        let mut previous = 0.0;
        let mut accumulated_j = 0.0;
        for (ordinal, end) in [0.005, 0.01].into_iter().enumerate() {
            let dt = end - previous;
            let center_z = -0.25 + 25.0 * end;
            // Each constant plateau face has load Q/3 at its three nodes.
            // Their masses are 1/24 and 1/12 kg. Absorb this known external
            // impulse into the old h supplied to the existing independent
            // P1-source + nonlinear contact quadratic, never a thermal solver.
            let mut radiated = old;
            for vertex in [0, 1, 3] { radiated[vertex] -= dt * 8.0 * radiation_w[0]; }
            for vertex in [4, 5, 7] { radiated[vertex] -= dt * 4.0 * radiation_w[1]; }
            let (expected, expected_heat) = endpoint(&radiated, dt, center_z);
            let step = scan.advance_to(cx, end, config().step).unwrap();
            let ScanTransport::Contact(contact) = &step.transport else {
                panic!("phase-dependent contact was omitted from the radiating scan");
            };
            let solved = step.solution();
            for i in 0..8 {
                assert!((0.0..100.0).contains(&expected[i]));
                near(solved.specific_enthalpy_j_kg[i], expected[i], 2e-8);
                near(solved.liquid_mass_fraction[i], expected[i] / 100.0, 2e-10);
            }
            assert_eq!(solved.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
            near(step.deposition[0].center_m[2], center_z, 1e-14);
            near(step.deposition[0].delivered_power_w, POWER_W, 1e-12);
            near(solved.source_w, POWER_W, 1e-12);
            let radiation = solved.pointwise_radiation.as_ref().unwrap();
            assert_eq!(radiation.patches.len(), 2);
            near(radiation.outward_heat_w, radiation_w.iter().sum(), 1e-9);
            near(contact.physical.radiation_out_w, radiation_w.iter().sum(), 1e-9);
            assert!(contact.radiation.is_none() && contact.combined_boundary.is_none());

            let traces = contact.physical_contacts.material_trace_states();
            assert!(!traces.is_empty());
            for trace in traces {
                near(trace.liquid_mass_fractions[0], expected[..3].iter().sum::<f64>() / 300.0, 2e-10);
                near(trace.liquid_mass_fractions[1], expected[4..7].iter().sum::<f64>() / 300.0, 2e-10);
            }
            let flux = contact.physical_contacts.frozen_interfaces().fluxes(&solved.temperature).unwrap();
            near(-flux[0].heat_rate_a_to_b_w, expected_heat, 2e-8);
            // Reusing old fractions would leave the same plateau temperatures
            // but produce the wrong contact heat and all three coupled loads.
            let frozen_heat = 50.0 / (1.0
                + 0.02 * old[..3].iter().sum::<f64>() / 3.0
                + 0.01 * old[4..7].iter().sum::<f64>() / 3.0);
            assert!((expected_heat - frozen_heat).abs() > 1e-3);
            let input_j = dt * (POWER_W - radiation_w.iter().sum::<f64>());
            let stored_j = old.iter().zip(&solved.specific_enthalpy_j_kg).enumerate()
                .map(|(i, (&before, &after))| (after - before) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 })
                .sum::<f64>();
            near(stored_j, input_j, 1e-9);
            near(solved.stored_energy_change_j, input_j, 1e-9);
            near(contact.physical.energy_residual_j, 0.0, config().step.energy_tolerance_j);
            assert!(contact.physical.residual_norm_j <= config().residual_tolerance_j);
            accumulated_j += input_j;
            near(solved.specific_enthalpy_j_kg.iter().zip(OLD).enumerate()
                .map(|(i, (&after, before))| (after - before) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 })
                .sum(), accumulated_j, 2e-9);
            assert_eq!(scan.conduction().accepted_steps(), ordinal + 1);
            old.copy_from_slice(&solved.specific_enthalpy_j_kg);
            previous = end;
        }
    });
}
