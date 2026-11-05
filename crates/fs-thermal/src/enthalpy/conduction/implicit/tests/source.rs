use super::*;

struct Feedback<'a> {
    caloric: &'a dyn ThermalLaw,
    coefficients: [f64; 3],
    bounds: [f64; 2],
}

impl ImplicitHeatSourceLaw for Feedback<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        if !(self.bounds[0]..=self.bounds[1]).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        self.caloric.temperature_k(h)?;
        let x = h - 300.0;
        let [a, b, c] = self.coefficients;
        Ok(ImplicitHeatSourceState {
            power_w: a + b * x + c * x * x,
            power_enthalpy_derivative_w_kg_j: b + 2.0 * c * x,
        })
    }
}

fn feedback(law: &dyn ThermalLaw, coefficients: [f64; 3]) -> Feedback<'_> {
    Feedback {
        caloric: law,
        coefficients,
        bounds: [100.0, 1000.0],
    }
}

#[test]
fn endpoint_generation_and_authored_power_share_the_original_conductive_balance() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let authored = [3.0, -1.0]
        .into_iter()
        .enumerate()
        .map(|(i, source_w)| TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(
                CellId(i as u64),
                Boundary {
                    source_w,
                    ..boundary(None)
                },
            ),
        })
        .collect();
    let (mut thermal, observer) = setup(&law, &[300.0, 320.0], None, authored);
    let masses = [2.0, 3.0];
    for (i, mass) in masses.into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    let heating = feedback(&law, [10.0, 0.2, 0.0]);
    let cooling = feedback(&law, [-5.0, -0.1, 0.0]);
    let sources: &[(CellId, &dyn ImplicitHeatSourceLaw)] =
        &[(CellId(0), &heating), (CellId(1), &cooling)];
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, sources, || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_heat_sources(&thermal, sources, || Ok(()))
    .unwrap();
    let dt = 0.5;
    // Independent two-cell balance in x = T - 300. Series Fourier G = 1.
    let a = masses[0] + 0.8 * dt;
    let d = masses[1] + 1.1 * dt;
    let b0 = 13.0 * dt;
    let b1 = masses[1] * 20.0 - 6.0 * dt;
    let determinant = a * d - dt * dt;
    let expected = [
        (d * b0 + dt * b1) / determinant,
        (dt * b0 + a * b1) / determinant,
    ];
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let chain_result = chain
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let graph_result = graph
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    for (next, report) in [chain_result, graph_result] {
        for (i, x) in expected.into_iter().enumerate() {
            near(
                next.temperature_k(CellId(i as u64)).unwrap(),
                300.0 + x,
                1e-9,
            );
        }
        let energy_change = masses[0] * expected[0] + masses[1] * (expected[1] - 20.0);
        let endpoint_heat = dt * (13.0 + 0.2 * expected[0] - 6.0 - 0.1 * expected[1]);
        near(report.energy.external_heat_j, endpoint_heat, 1e-9);
        near(report.energy.stored_energy_change_j, energy_change, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        assert_eq!(report.energy.substeps, 1);
        assert!(report.entropy.has_unclassified_sources);
        // The endpoint source is transient solver input, not authored history.
        near(next.cells()[&CellId(0)].boundary.source_w, 3.0, 0.0);
    }
    near(thermal.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
}

#[test]
fn nonlinear_source_and_fourier_tangents_match_the_original_joule_residual() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (mut thermal, observer) = setup(&law, &[400.0, 350.0], Some(300.0), vec![]);
    for (i, mass) in [1.5, 2.5].into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    let heating = feedback(&law, [10.0, 0.2, 0.003]);
    let cooling = feedback(&law, [-5.0, -0.1, -0.002]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(
                &thermal,
                &[(CellId(0), &heating), (CellId(1), &cooling)],
                || Ok(()),
            )
            .unwrap();
    let dt = 0.7;
    let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
    let energy = [-10.0, 3.0];
    let mut report = empty_report();
    let trial = chain
        .evaluate(
            &step,
            &energy,
            &BTreeMap::new(),
            &[],
            false,
            controls(1.0),
            &mut report,
            &mut || Ok(()),
        )
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    chain
        .jacobian(
            &thermal,
            trial.payload(),
            dt,
            &mut workspace,
            &mut || Ok(()),
        )
        .unwrap();
    for column in 0..2 {
        let mut plus = energy;
        let mut minus = energy;
        plus[column] += 1e-3;
        minus[column] -= 1e-3;
        let a = chain
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        let b = chain
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        for row in 0..2 {
            let analytic = if row == column {
                workspace.diagonal[row]
            } else if row > column {
                workspace.lower[column]
            } else {
                workspace.upper[row]
            };
            near(
                analytic,
                (a.residual_j()[row] - b.residual_j()[row]) / 0.002,
                1e-8,
            );
        }
    }
}

#[test]
fn latent_phase_generation_uses_stored_energy_even_with_constant_temperature() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let (thermal, observer) = setup(&law, &[350.0], Some(300.0), vec![]);
    let heating = feedback(&law, [80.0, -0.2, 0.0]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, &[(CellId(0), &heating)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let expected_h = 490.0 / 1.2;
    near(next.cells()[&CellId(0)].specific_h_j_kg, expected_h, 1e-9);
    near(next.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
    near(report.energy.external_heat_j, expected_h - 350.0, 1e-9);
    near(report.energy.residual_j, 0.0, 1e-9);
}

#[test]
fn endpoint_source_survives_both_sdirk_stages_and_adaptive_step_doubling() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(300.0), vec![]);
    // Source cooling rate 3 plus Fourier bath rate 1 gives x' = -4 x.
    let cooling = feedback(&law, [0.0, -3.0, 0.0]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, &[(CellId(0), &cooling)], || Ok(()))
            .unwrap();
    let exact = 300.0 + 100.0 * (-4.0_f64).exp();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (next, report) = chain
            .advance_sdirk2_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let t = next.temperature_k(CellId(0)).unwrap();
        errors.push((t - exact).abs());
        near(report.energy.external_heat_j, t - 400.0, 1e-8);
        near(report.energy.residual_j, 0.0, 1e-8);
    }
    assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
    let adaptive = adaptive::AdaptiveImplicitConductionControls {
        step: controls(0.5),
        initial_dt_s: 0.5,
        enthalpy_absolute_tolerance_j_kg: 1e-4,
        temperature_absolute_tolerance_k: 1e-4,
        relative_tolerance: 1e-4,
        ..Default::default()
    };
    let (next, _) = chain
        .advance_adaptive_to(&thermal, 1.0, adaptive, &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), exact, 0.02);
}

#[test]
fn source_attachment_and_late_domain_refusal_leave_the_original_unchanged() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let foreign = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0, 400.0], None, vec![]);
    let chain = || {
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
    };
    let source = feedback(&law, [1.0, 0.0, 0.0]);
    let wrong = feedback(&foreign, [1.0, 0.0, 0.0]);
    let nonfinite = feedback(&law, [f64::NAN, 0.0, 0.0]);
    assert!(chain().with_heat_sources(&thermal, &[], || Ok(())).is_err());
    assert!(
        chain()
            .with_heat_sources(&thermal, &[(CellId(2), &source)], || Ok(()))
            .is_err()
    );
    assert!(
        chain()
            .with_heat_sources(
                &thermal,
                &[(CellId(0), &source), (CellId(0), &source)],
                || Ok(())
            )
            .is_err()
    );
    assert!(
        chain()
            .with_heat_sources(&thermal, &[(CellId(0), &wrong)], || Ok(()))
            .is_err()
    );
    assert!(
        chain()
            .with_heat_sources(&thermal, &[(CellId(0), &nonfinite)], || Ok(()))
            .is_err()
    );
    let bounded = Feedback {
        bounds: [350.0, 1000.0],
        ..feedback(&law, [-1000.0, 0.0, 0.0])
    };
    let chain = chain()
        .with_heat_sources(&thermal, &[(CellId(0), &bounded)], || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let mut control = controls(1.0);
    control.maximum_halvings = 0;
    assert!(
        chain
            .advance_to(&thermal, 1.0, control, &mut workspace, || Ok(()))
            .is_err()
    );
    near(thermal.time_s(), 0.0, 0.0);
    near(thermal.temperature_k(CellId(0)).unwrap(), 400.0, 0.0);
    assert!(
        chain
            .advance_to(&thermal, 0.01, controls(0.01), &mut workspace, || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    near(thermal.temperature_k(CellId(1)).unwrap(), 400.0, 0.0);
}
