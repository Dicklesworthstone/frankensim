//! G1/G3/G4: actual simultaneous T/f source data in the native implicit graph.
use super::*;

fn alloy(scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(KP, KD),
            // k=s*(2+.02*(T-400)+4f), not a derived mixture rule.
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "f".into()],
                axis_dims: [Temperature::DIMS, Dims::NONE],
                samples: vec![
                    ([300.0, 0.0], 0.0),
                    ([800.0, 0.0], 10.0 * scale),
                    ([800.0, 1.0], 14.0 * scale),
                    ([300.0, 1.0], 4.0 * scale),
                ],
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                dims: KD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured alloy".into(),
            phase: "solidus/liquidus interval".into(),
            process: "fixed thermal inventory".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 450.0, 1.0),
            (450.0, 800.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap()
}
// Eliminate each finite hot body from BE. With center h=25+Q,
// k=a+bQ, T=412.5+.5Q and each hot T=500-Q/n:
// Q(k+2)=4*dt*n*k*(87.5-(.5+1/n)*Q).
fn heat(scale: f64, n: usize, dt: f64) -> f64 {
    let a = 3.25 * scale;
    let b = 0.05 * scale;
    let w = 4.0 * dt * n as f64;
    let d = 0.5 + 1.0 / n as f64;
    let aa = (1.0 + w * d) * b;
    let bb = (1.0 + w * d) * a + 2.0 - w * 87.5 * b;
    let cc = w * 87.5 * a;
    2.0 * cc / (bb + (bb * bb + 4.0 * aa * cc).sqrt())
}

#[test]
fn source_data_drives_finite_branched_alloy_heat_without_apparent_heat_capacity() {
    let mut endpoints = Vec::new();
    for scale in [0.5, 2.0] {
        let cards = [alloy(scale), card("hot finite source", 2.0, false)];
        let curves = [chart(&cards[0]), caloric(&cards[1])];
        let query = context();
        let laws = [
            material_law(&curves[0], &cards[0], &query, [0.0, 100.0]),
            material_law(&curves[1], &cards[1], &query, [-100.0, 500.0]),
        ];
        let old = initial(&curves, 3);
        let bulk = conduction(&old, &laws, 3);
        let bindings = [
            (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &laws[1] as &dyn ImplicitConductionLaw),
        ];
        let graph = ImplicitConductionGraph::new(
            &old,
            &bulk,
            &bindings,
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(4).unwrap();
        let (next, report) = graph
            .advance_to(&old, 0.1, controls(0.1), &mut workspace, || Ok(()))
            .unwrap();
        let q = heat(scale, 3, 0.1);
        assert!(0.0 < q && q < 75.0);
        near(next.cells()[&CENTER].specific_h_j_kg, 25.0 + q, 2e-8);
        near(next.temperature_k(CENTER).unwrap(), 412.5 + 0.5 * q, 2e-8);
        for hot in HOT {
            near(next.cells()[&hot].specific_h_j_kg, 200.0 - q / 3.0, 2e-8);
        }
        let response = laws[0].state_at_specific_enthalpy(25.0 + q).unwrap();
        near(response.liquid_mass_fraction(), (25.0 + q) / 100.0, 2e-13);
        near(
            response.conductivity_enthalpy_derivative_kg_w_m_k_j(),
            0.05 * scale,
            2e-13,
        );
        near(response.temperature_enthalpy_derivative_kg_k_j(), 0.5, 0.0);
        near(next.energy_j(), old.energy_j(), 2e-8);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.energy.external_heat_j, 0.0);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        endpoints.push(next.cells()[&CENTER].specific_h_j_kg);

        let mut calls = 0;
        graph
            .advance_to(&old, 0.1, controls(0.1), &mut workspace, || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        let mut retry_calls = 0;
        assert!(
            graph
                .advance_to(&old, 0.1, controls(0.1), &mut workspace, || {
                    retry_calls += 1;
                    if retry_calls == calls {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CENTER].specific_h_j_kg, 25.0);
        let (retry, _) = graph
            .advance_to(&old, 0.1, controls(0.1), &mut workspace, || Ok(()))
            .unwrap();
        for (&id, cell) in next.cells() {
            assert_eq!(
                cell.specific_h_j_kg.to_bits(),
                retry.cells()[&id].specific_h_j_kg.to_bits()
            );
        }
    }
    assert!(endpoints[1] > endpoints[0]);
}
