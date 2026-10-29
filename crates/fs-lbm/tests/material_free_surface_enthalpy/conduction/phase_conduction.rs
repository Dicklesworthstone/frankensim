//! The original declared phase law reaches heat, momentum and atomic admission.
use super::*;
use fs_material::phase::conductivity::LiquidMassFractionConductivity;

#[test]
fn g1_declared_phase_conductivity_matches_independent_finite_heat_balance() {
    let card = conducting_card(1.0);
    for multiplier in [0.25, 1.0, 2.0] {
        with_liquid(&card, false, |thermal, gate| {
            let law = LiquidMassFractionConductivity::declared(4.0, multiplier).unwrap();
            *thermal = thermal.clone().with_phase_conductivity(law, gate).unwrap();
            assert_eq!(thermal.phase_conductivity_law(), Some(law));
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_conduction(&[face()], &[], controls(), gate)
                .unwrap();

            // Two 1 mg cells, Cp=2000 J/(kg K). Apply the declared material
            // multiplier to each physical half-cell resistance independently.
            let mut t = [350.0, 450.0];
            let mut last = 0.0;
            for _ in 0..4 {
                let k = t.map(|t| multiplier * (t - 200.0));
                last = 1e-6 / (0.0005 / k[0] + 0.0005 / k[1]);
                let q = 0.00025 * last * (t[1] - t[0]);
                t[0] += q / 0.002;
                t[1] -= q / 0.002;
            }
            assert_eq!(report.thermal.substeps, 4);
            near(report.last_applied[0].conductance_w_k, last, 2e-13);
            for i in 0..2 {
                let endpoint = &report.thermal_endpoint[0];
                let raw = &endpoint.materials[i];
                let applied = endpoint.conductivities[i];
                assert_eq!(raw.material().card_identity(), card.content_hash());
                near(temperature(raw), t[i], 2e-13);
                near(raw.thermal_conductivity_w_m_k(), t[i] - 200.0, 2e-13);
                near(applied.source_w_m_k, t[i] - 200.0, 2e-13);
                assert_eq!(applied.phase_multiplier, multiplier);
                near(applied.effective_w_m_k, multiplier * (t[i] - 200.0), 2e-13);
                near(
                    thermal
                        .conductivity_at(PAIR[i])
                        .unwrap()
                        .unwrap()
                        .effective_w_m_k,
                    multiplier * (t[i] - 200.0),
                    2e-12,
                );
                near(
                    thermal.specific_enthalpy_at(PAIR[i]).unwrap(),
                    (t[i] - 300.0) * 2000.0,
                    2e-12,
                );
            }
            assert_eq!(report.flow.external_heat_j, 0.0);
            near(thermal.total_enthalpy_j().unwrap(), before, 2e-12);
            assert!(report.flow.balance_defect_j.abs() < 1e-9);
        });
    }
}

#[test]
fn g3_phase_data_changes_native_flow_and_identity_retains_original_bits() {
    let card = conducting_card(1.0);
    with_liquid(&card, true, |thermal, gate| {
        let mut original = thermal.clone();
        let mut identity = thermal
            .clone()
            .with_phase_conductivity(
                LiquidMassFractionConductivity::declared(1.0, 1.0).unwrap(),
                gate,
            )
            .unwrap();
        original
            .step_with_conduction(&[face()], &[], controls(), gate)
            .unwrap();
        identity
            .step_with_conduction(&[face()], &[], controls(), gate)
            .unwrap();
        assert_eq!(
            original.flow().lattice().grid().state_bits(),
            identity.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            original.total_enthalpy_j().unwrap(),
            identity.total_enthalpy_j().unwrap()
        );
        assert_eq!(
            original.specific_enthalpy_at(PAIR[0]),
            identity.specific_enthalpy_at(PAIR[0])
        );
    });
    let run = |multiplier| {
        with_liquid(&card, true, |thermal, gate| {
            *thermal = thermal
                .clone()
                .with_phase_conductivity(
                    LiquidMassFractionConductivity::declared(4.0, multiplier).unwrap(),
                    gate,
                )
                .unwrap();
            let energy = thermal.total_enthalpy_j().unwrap();
            for _ in 0..3 {
                thermal
                    .step_with_conduction(&[face()], &[], controls(), gate)
                    .unwrap();
            }
            near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s,
                temperature(thermal.material_at(PAIR[0]).unwrap()),
                thermal
                    .material_at(PAIR[0])
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
            )
        })
    };
    let slow = run(0.25);
    let fast = run(2.0);
    assert_ne!(slow.0, fast.0);
    assert!((slow.1[1] - fast.1[1]).abs() > 1e-10);
    assert!(fast.2 > slow.2);
    assert!(fast.3 < slow.3);
}

#[test]
fn g4_phase_conductivity_overflow_and_cancellation_preserve_native_history() {
    let card = conducting_card(1.0);
    with_liquid(&card, false, |thermal, gate| {
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        assert!(
            thermal
                .clone()
                .with_phase_conductivity(
                    LiquidMassFractionConductivity::declared(1.0, f64::MAX).unwrap(),
                    gate,
                )
                .is_err()
        );
        assert_eq!(thermal.phase_conductivity_law(), None);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);

        // The actual source k is initially <=250. Heating the hot cell by
        // 0.06 J gives 480 K and k=280, overflowing only the declared product.
        // The original material and caloric domains still support that state.
        let law = LiquidMassFractionConductivity::declared(1.0, f64::MAX / 260.0).unwrap();
        *thermal = thermal.clone().with_phase_conductivity(law, gate).unwrap();
        let mut clean = thermal.clone();
        let coefficient = thermal.conductivity_at(PAIR[1]).unwrap();
        let source = thermal.material_at(PAIR[1]).unwrap().clone();
        let unchanged = |candidate: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
            assert_eq!(candidate.flow().lattice().grid().state_bits(), bits);
            assert_eq!(candidate.total_enthalpy_j().unwrap(), energy);
            assert_eq!(candidate.specific_enthalpy_at(PAIR[1]), Some(300000.0));
            assert_eq!(candidate.material_at(PAIR[1]).unwrap(), &source);
            assert_eq!(candidate.conductivity_at(PAIR[1]).unwrap(), coefficient);
            assert_eq!(candidate.accepted_steps(), 0);
            assert_eq!(candidate.physical_time_s(), 0.0);
        };
        assert!(matches!(thermal.step_with_heat(&[FreeSurfaceCellHeat {
            coordinates: PAIR[1], energy_j: 0.06,
        }], gate), Err(FreeSurfaceEnthalpyError::Invalid(message)) if message.contains("conductivity")));
        unchanged(thermal);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step(&cancelled).is_err());
        unchanged(thermal);
        clean.step(gate).unwrap();
        thermal.step(gate).unwrap();
        assert_eq!(
            thermal.flow().lattice().grid().state_bits(),
            clean.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            thermal.total_enthalpy_j().unwrap(),
            clean.total_enthalpy_j().unwrap()
        );
        assert_eq!(
            thermal.conductivity_at(PAIR[1]).unwrap(),
            clean.conductivity_at(PAIR[1]).unwrap()
        );
        assert_eq!(thermal.accepted_steps(), 1);
    });
}
