//! Native goal adaptation must re-query the physical contact on both fields.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

fn solved(f: &Fixture, ordinal: usize) -> (J, J, J) {
    let (source, contact, ledger) = prepare(f, ordinal, None);
    let output = if f
        .project
        .cooling
        .as_ref()
        .unwrap()
        .conduction
        .as_ref()
        .unwrap()
        .radiation
        .is_some()
    {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        fs_cli::run(vec![
            "--json".into(),
            "solve".into(),
            source.display().to_string(),
            ledger.display().to_string(),
            "--materials".into(),
            f.dir.join("solid.fsmcdpk").display().to_string(),
            "--materials".into(),
            root.join("data/reference-project/gray-surface.fsmcdpk")
                .display()
                .to_string(),
            "--interfaces".into(),
            contact.display().to_string(),
        ])
    } else {
        run(f, &source, &contact, &ledger)
    };
    assert_eq!(
        output.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        output.stdout,
        output.stderr
    );
    let output = J::parse(&output.stdout).unwrap();
    let receipt = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stages = receipt.get("stages").unwrap().as_array().unwrap();
    let stage = |name| {
        let row = stages
            .iter()
            .find(|row| row.str_field("stage") == Some(name))
            .unwrap();
        artifact(&ledger, row.str_field("receipt").unwrap())
    };
    let conduction = stage("conduction");
    let field = artifact(&ledger, conduction.str_field("solution_artifact").unwrap());
    (conduction, field, stage("qoi"))
}

fn discretization(qoi: &J) -> &J {
    qoi.get("budget").unwrap().as_array().unwrap()[0]
        .get("terms")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|term| term.str_field("kind") == Some("discretization"))
        .unwrap()
}

fn maximum(field: &J) -> f64 {
    field
        .get("temperature")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .fold(f64::NEG_INFINITY, f64::max)
}

#[test]
fn live_contact_adaptive_goals_measure_the_full_residual_and_publish_discretization() {
    // Exercise a fixed reservoir and prescribed Robin exchange, with K'(T)
    // present in the latter. Both use the same live R''(T) = 0.001 T source.
    for nonlinear in [false, true] {
        let mut f = Fixture::new(nonlinear, false);
        if nonlinear {
            f.project
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .boundaries[0]
                .condition = B::Convection {
                coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
            };
        }
        f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
        f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
        f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
        let (receipt, field, qoi) = solved(&f, 80);
        let adaptive = receipt
            .get("adaptive")
            .expect("complete live-contact goal producer executed");
        assert_eq!(adaptive.str_field("status"), Some("observed-tolerance-met"));
        assert!(adaptive.f64_field("solved_meshes").unwrap() > 1.0);
        assert_eq!(adaptive.get("continuum_error_bound"), Some(&J::Bool(false)));
        assert_eq!(
            receipt.get("solver_control"),
            Some(&J::Null),
            "a nonlinear contact comparison must not claim a frozen-operator inverse certificate"
        );
        let history = adaptive.get("history").unwrap().as_array().unwrap();
        assert!(!history.is_empty());
        for row in history {
            assert!(row.f64_field("interface_pairs").unwrap() > 0.0);
            assert_eq!(
                row.f64_field("probe_interface_pairs").unwrap(),
                4.0 * row.f64_field("interface_pairs").unwrap()
            );
            assert_eq!(
                row.get("uses_nonlinear_jacobian"),
                Some(&J::Bool(true)),
                "temperature-dependent contact requires its complete tangent even for constant bulk K"
            );
            let balance = row.get("algebraic_balance").unwrap();
            assert_eq!(balance.str_field("status"), Some("unsupported-model"));
            assert_eq!(
                balance.str_field("comparison_model"),
                Some("source-contact-full-transpose")
            );
            assert_eq!(
                balance.get("contact_temperature_feedback"),
                Some(&J::Bool(true))
            );
            assert!(row.f64_field("probe_tets").unwrap() > row.f64_field("tets").unwrap());
            assert!(row.f64_field("absolute_contribution_sum_k").unwrap() > 0.0);
            let linear = row.f64_field("signed_linear_change_k").unwrap();
            let remainder = row.f64_field("linearization_remainder_k").unwrap();
            let ownership = row.f64_field("maximum_remainder_k").unwrap();
            let measured = row.f64_field("measured_change_k").unwrap();
            assert!((linear + remainder + ownership - measured).abs() < 1e-9);
            assert!(row.f64_field("estimated_change_k").unwrap() + 1e-12 >= measured.abs());
            assert!(row.f64_field("dual_residual").unwrap() <= 1e-9);
        }
        assert!(
            history
                .iter()
                .any(|row| row.f64_field("linearization_remainder_k").unwrap().abs() > 1e-12),
            "the finite contact/material remainder must not disappear into a frozen-resistance estimate"
        );
        let last = history.last().unwrap();
        assert!((maximum(&field) - last.f64_field("t_max_k").unwrap()).abs() < 1e-9);
        let term = discretization(&qoi);
        assert_eq!(term.str_field("state"), Some("interval"));
        assert_eq!(
            term.get("derivation").unwrap().str_field("method"),
            Some("adaptive-enrichment-order-one")
        );
        let expected = 2.0
            * last
                .f64_field("estimated_change_k")
                .unwrap()
                .abs()
                .max(last.f64_field("measured_change_k").unwrap().abs());
        assert_eq!(
            term.f64_field("upper_kelvin").unwrap().to_bits(),
            expected.to_bits()
        );

        // A separately requested uniform ladder must reproduce the adaptive
        // probe on its identical first refined mesh and physical source law.
        f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
        let (uniform, _, _) = solved(&f, 81);
        let rungs = uniform
            .get("ladder")
            .unwrap()
            .get("rungs")
            .unwrap()
            .as_array()
            .unwrap();
        assert_eq!(rungs.len(), 3);
        let first = &history[0];
        assert_eq!(first.f64_field("probe_tets"), rungs[1].f64_field("tets"));
        assert!(
            (first.f64_field("probe_t_max_k").unwrap() - rungs[1].f64_field("t_max_k").unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (rungs[2].f64_field("t_max_k").unwrap() - maximum(&field)).abs()
                <= term.f64_field("upper_kelvin").unwrap() + 1e-9,
            "this independently resolved finer fixture must fit its observed discretization allowance"
        );

        f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
        let (replayed, replayed_field, _) = solved(&f, 82);
        assert_eq!(receipt.get("adaptive"), replayed.get("adaptive"));
        assert_eq!(field.get("temperature"), replayed_field.get("temperature"));
    }
}

#[test]
fn live_contact_adaptive_memory_limit_retains_the_last_probed_mesh_without_a_bound() {
    let mut f = Fixture::new(false, false);
    f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
    f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
    let (baseline, _, _) = solved(&f, 83);
    let first = &baseline
        .get("adaptive")
        .unwrap()
        .get("history")
        .unwrap()
        .as_array()
        .unwrap()[0];
    // Admit exactly the imported mesh's first required enrichment. Every
    // subsequent refinement then exceeds the next-probe memory gate before
    // local shape quality can become the stopping reason.
    let probe_tets = first.f64_field("probe_tets").unwrap() as u64;
    f.project.budgets.as_mut().unwrap().accuracy_rel = 1e-12;
    f.project.budgets.as_mut().unwrap().memory_bytes = probe_tets.checked_mul(256).unwrap();
    let (receipt, field, qoi) = solved(&f, 86);
    let adaptive = receipt.get("adaptive").unwrap();
    assert_eq!(adaptive.str_field("status"), Some("unresolved"));
    assert_eq!(adaptive.str_field("stop"), Some("memory-budget"));
    let history = adaptive.get("history").unwrap().as_array().unwrap();
    let last = history
        .last()
        .expect("the memory declaration admits a complete enriched comparison");
    assert!((maximum(&field) - last.f64_field("t_max_k").unwrap()).abs() < 1e-9);
    assert_eq!(
        receipt.get("interfaces").unwrap().f64_field("pair_count"),
        last.f64_field("interface_pairs")
    );
    assert_eq!(
        adaptive.f64_field("last_estimated_change_k"),
        last.f64_field("estimated_change_k")
    );
    assert_eq!(discretization(&qoi).str_field("state"), Some("no-data"));
}

#[test]
fn live_contact_adaptive_source_corner_refuses_an_invented_smooth_goal_tangent() {
    let mut f = Fixture::new(false, false);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = B::FixedTemperature {
        temperature: QtyAny::new(300.0, dims::TEMPERATURE),
    };
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 300.0;
    envelope.ambient_hi.value = 300.0;
    let source_law =
        || temperature_pack_from_knots(vec![(200.0, 0.2), (300.0, 0.25), (450.0, 0.55)]);
    let (source, contact, ledger) = prepare_with_pack(&f, 84, None, source_law());
    let base = run(&f, &source, &contact, &ledger);
    assert_eq!(
        base.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        base.stdout,
        base.stderr
    );
    f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
    let (source, contact, ledger) = prepare_with_pack(&f, 85, None, source_law());
    let refused = run(&f, &source, &contact, &ledger);
    assert_ne!(refused.exit_code, fs_cli::exit::SUCCESS);
    let diagnostic = format!("{}{}", refused.stdout, refused.stderr);
    assert!(
        diagnostic.contains("derivative") || diagnostic.contains("slope"),
        "{diagnostic}"
    );
    assert!(
        !diagnostic.contains("contact-composition"),
        "the source itself must refuse, rather than the former blanket adaptive-composition gate"
    );
}

#[test]
fn live_contact_adaptive_natural_air_and_mixed_goals_retain_complete_boundary_feedback() {
    for (natural, air, nonlinear) in [
        (true, false, false),
        (false, true, true),
        (true, true, true),
    ] {
        let mut f = Fixture::new(nonlinear, air);
        let setup = f
            .project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap();
        setup.radiation = None;
        if natural {
            setup.boundaries[usize::from(air)].condition = B::NaturalConvection {
                characteristic_length: QtyAny::new(1.0, dims::LENGTH),
                ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".into(),
            };
        } else {
            // The downstream hot wall must see heat from the upstream cold
            // wall through the original ordered air-reference map.
            let mut downstream = setup.boundaries[0].condition.clone();
            let B::AirflowConvection { order, .. } = &mut downstream else {
                panic!("air fixture");
            };
            *order = 1;
            setup.boundaries[1].condition = downstream;
        }
        let envelope = f.project.envelope.as_mut().unwrap();
        envelope.ambient_lo.value = 297.0;
        envelope.ambient_hi.value = 297.0;
        f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
        f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
        f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
        let (receipt, field, qoi) = solved(&f, 90);
        let adaptive = receipt.get("adaptive").unwrap();
        assert_eq!(adaptive.str_field("status"), Some("observed-tolerance-met"));
        assert!(adaptive.f64_field("solved_meshes").unwrap() > 1.0);
        let history = adaptive.get("history").unwrap().as_array().unwrap();
        assert!(!history.is_empty());
        for row in history {
            assert_eq!(
                row.get("natural_coefficient_feedback"),
                Some(&J::Bool(natural))
            );
            assert_eq!(row.get("uses_nonlinear_jacobian"), Some(&J::Bool(true)));
            let balance = row.get("algebraic_balance").unwrap();
            assert_eq!(
                balance.str_field("comparison_model"),
                Some("source-contact-full-transpose")
            );
            assert_eq!(
                balance.get("contact_temperature_feedback"),
                Some(&J::Bool(true))
            );
            assert_eq!(balance.get("air_reference_feedback"), Some(&J::Bool(air)));
            assert_eq!(
                row.get("coupled_adjoint"),
                Some(&J::Null),
                "the complete FGMRES transpose must not be reported as the legacy interface iteration"
            );
            let linear = row.f64_field("signed_linear_change_k").unwrap();
            let nonlinear = row.f64_field("linearization_remainder_k").unwrap();
            let ownership = row.f64_field("maximum_remainder_k").unwrap();
            let measured = row.f64_field("measured_change_k").unwrap();
            assert!((linear + nonlinear + ownership - measured).abs() < 1e-9);
            assert!(row.f64_field("estimated_change_k").unwrap() + 1e-12 >= measured.abs());
            assert!(row.f64_field("dual_residual").unwrap() <= 1e-9);
        }
        // The receipt carries `natural_convection` only when a natural law
        // runs; without one the key is absent, not null (fs-cli solve.rs).
        assert_eq!(
            receipt.get("natural_convection").is_some(),
            natural
        );
        assert_eq!(receipt.get("conjugate").unwrap() != &J::Null, air);
        let contact = receipt.get("temperature_contact").unwrap();
        if air {
            let exchange = receipt.get("conjugate").unwrap();
            let segments = exchange.get("segments").unwrap().as_array().unwrap();
            assert_eq!(segments.len(), if natural { 1 } else { 2 });
            assert!(
                (exchange.f64_field("solid_total_w").unwrap()
                    - exchange.f64_field("air_total_w").unwrap())
                .abs()
                    <= exchange.f64_field("balance_tolerance_w").unwrap()
            );
            // These heat rates are independently computed by the solid and
            // physical air march. Their absolute discrepancies must fit the
            // contact's watt gate, without cancellation between segments.
            let imbalance: f64 = segments
                .iter()
                .map(|segment| segment.f64_field("imbalance_w").unwrap().abs())
                .sum();
            assert!(
                imbalance <= contact.f64_field("residual_tolerance_w").unwrap(),
                "the accepted exchange must close at the live-contact physical tolerance"
            );
        }
        assert!(
            contact.f64_field("residual_norm_w").unwrap()
                <= contact.f64_field("residual_tolerance_w").unwrap()
        );
        assert!(
            contact.f64_field("energy_closure_w").unwrap().abs()
                <= contact.f64_field("energy_tolerance_w").unwrap()
        );
        for face in contact.get("face_states").unwrap().as_array().unwrap() {
            let query = face.f64_field("query_temperature_k").unwrap();
            assert!((face.f64_field("resistance_m2_k_w").unwrap() - 0.001 * query).abs() < 1e-14);
        }
        assert!(
            (maximum(&field) - history.last().unwrap().f64_field("t_max_k").unwrap()).abs() < 1e-9
        );
        assert_eq!(discretization(&qoi).str_field("state"), Some("interval"));
        if natural && air {
            let (replayed, replayed_field, _) = solved(&f, 91);
            assert_eq!(receipt.get("adaptive"), replayed.get("adaptive"));
            assert_eq!(field.get("temperature"), replayed_field.get("temperature"));
        }
        // An independently requested two-rung ladder resolves the same first
        // probe. Its ordinary coupling stop is looser than the adaptive stop;
        // a microkelvin allows that declared difference without hiding a
        // frozen contact, natural coefficient or downstream air reference.
        let first = &history[0];
        f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
        f.project.budgets.as_mut().unwrap().memory_bytes = (first.f64_field("probe_tets").unwrap()
            as u64)
            .checked_mul(256)
            .unwrap();
        let (uniform, _, _) = solved(&f, 92);
        let rungs = uniform
            .get("ladder")
            .unwrap()
            .get("rungs")
            .unwrap()
            .as_array()
            .unwrap();
        assert_eq!(rungs.len(), 2);
        assert_eq!(first.f64_field("probe_tets"), rungs[1].f64_field("tets"));
        assert!(
            (first.f64_field("probe_t_max_k").unwrap() - rungs[1].f64_field("t_max_k").unwrap())
                .abs()
                < 1e-6,
            "natural={natural} air={air}: the independent physical uniform solve disagrees with the adaptive probe"
        );
    }
}

#[test]
fn live_contact_adaptive_radiation_keeps_vacuum_natural_and_air_heat_physical() {
    let number = |row: &J, key: &str| row.f64_field(key).unwrap();
    // The same imported bodies and source curve exercise real radiation-only
    // cooling and each original convection owner, with nonlinear bulk K(T).
    for mode in ["vacuum", "robin", "natural", "air"] {
        let mut f = Fixture::new(true, true);
        let setup = f
            .project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap();
        let bath = setup.radiation.as_ref().unwrap().surfaces[0]
            .reservoir_temperature
            .value;
        setup.boundaries[0].condition = match mode {
            "vacuum" => B::HeatFlux {
                outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
            },
            "robin" => B::Convection {
                coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
            },
            "natural" => B::NaturalConvection {
                characteristic_length: QtyAny::new(1.0, dims::LENGTH),
                ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".into(),
            },
            "air" => setup.boundaries[0].condition.clone(),
            _ => unreachable!(),
        };
        // One fixed ambient, as in radiation.rs's fixture: the project refuses
        // a natural-convection ambient or an air inlet outside the envelope's
        // ambient range, and that ambient also sets the air density. Where the
        // convection has its own ambient it is that one; otherwise the bath.
        let ambient = match &setup.boundaries[0].condition {
            B::NaturalConvection { ambient_temperature, .. } => ambient_temperature.value,
            B::AirflowConvection { inlet_temperature, .. } => inlet_temperature.value,
            _ => bath,
        };
        let envelope = f.project.envelope.as_mut().unwrap();
        envelope.ambient_lo.value = ambient;
        envelope.ambient_hi.value = ambient;
        f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
        f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
        f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
        let (receipt, field, qoi) = solved(&f, 100);
        let adaptive = receipt.get("adaptive").unwrap();
        assert_eq!(
            adaptive.str_field("status"),
            Some("observed-tolerance-met"),
            "{mode}"
        );
        assert!(number(adaptive, "solved_meshes") > 1.0);
        assert_eq!(adaptive.get("continuum_error_bound"), Some(&J::Bool(false)));
        assert_eq!(receipt.get("solver_control"), Some(&J::Null));
        let history = adaptive.get("history").unwrap().as_array().unwrap();
        assert!(!history.is_empty());
        for row in history {
            let balance = row.get("algebraic_balance").unwrap();
            assert_eq!(
                balance.str_field("comparison_model"),
                Some("source-contact-radiation-full-transpose")
            );
            assert_eq!(balance.str_field("status"), Some("unsupported-model"));
            assert_eq!(
                balance.get("contact_temperature_feedback"),
                Some(&J::Bool(true))
            );
            assert_eq!(balance.get("radiation_feedback"), Some(&J::Bool(true)));
            assert_eq!(
                balance.get("air_reference_feedback"),
                Some(&J::Bool(mode == "air"))
            );
            assert_eq!(
                row.get("natural_coefficient_feedback"),
                Some(&J::Bool(mode == "natural"))
            );
            assert_eq!(row.get("uses_nonlinear_jacobian"), Some(&J::Bool(true)));
            assert_eq!(row.get("coupled_adjoint"), Some(&J::Null));
            assert!(number(row, "interface_pairs") > 0.0);
            assert_eq!(
                number(row, "probe_interface_pairs"),
                4.0 * number(row, "interface_pairs")
            );
            assert!(number(row, "absolute_contribution_sum_k") > 0.0);
            let measured = number(row, "measured_change_k");
            assert!(
                (number(row, "signed_linear_change_k")
                    + number(row, "linearization_remainder_k")
                    + number(row, "maximum_remainder_k")
                    - measured)
                    .abs()
                    < 1e-9
            );
            assert!(number(row, "estimated_change_k") + 1e-12 >= measured.abs());
            assert!(number(row, "dual_residual") <= 1e-9);
            if mode == "vacuum" {
                // With no convection or prescribed bath, the radiative
                // reservoir defines the rise; absolute kelvin is not an
                // admissible shortcut to a very loose accuracy target.
                let tolerance = 0.2 * (number(row, "t_max_k") - bath).abs();
                assert_eq!(number(row, "tolerance_k").to_bits(), tolerance.to_bits());
            }
        }
        let last = history.last().unwrap();
        assert!((maximum(&field) - number(last, "t_max_k")).abs() < 1e-9);
        assert_eq!(
            receipt.get("natural_convection").is_some(),
            mode == "natural"
        );
        assert_eq!(receipt.get("conjugate").unwrap() != &J::Null, mode == "air");
        let contact = receipt.get("temperature_contact").unwrap();
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        assert_eq!(
            contact
                .get("ambient_radiation")
                .unwrap()
                .get("physical_endpoint_reassembled"),
            Some(&J::Bool(true))
        );
        let (_, _, selected) = temperature_pack();
        for face in contact.get("face_states").unwrap().as_array().unwrap() {
            let mean = 0.5 * number(face, "side_a_temperature_k")
                + 0.5 * number(face, "side_b_temperature_k");
            assert_eq!(number(face, "query_temperature_k"), mean);
            assert!((number(face, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            assert_eq!(face.str_field("selected_claim"), Some(selected.as_str()));
            assert_eq!(face.str_field("usage_receipt").unwrap().len(), 64);
        }
        let radiation = receipt.get("radiation").unwrap();
        let convective = number(radiation, "convective_out_w");
        let radiative = number(radiation, "nonlinear_radiative_out_w");
        let energy = receipt.get("energy").unwrap();
        let source = number(energy, "source_w");
        assert!((source - 4.44).abs() < 1e-9);
        assert!((source - convective - radiative).abs() <= number(contact, "energy_tolerance_w"));
        assert!(
            radiative.abs() > 0.1,
            "{mode}: radiation must materially affect this heat balance"
        );
        if mode == "vacuum" {
            assert_eq!(
                convective, 0.0,
                "a radiation-only trace must not acquire a convection bath"
            );
        } else if mode == "air" {
            let exchange = receipt.get("conjugate").unwrap();
            assert!((number(exchange, "solid_total_w") - convective).abs() < 1e-9);
            assert!(
                (number(exchange, "air_total_w") - convective).abs()
                    <= number(contact, "residual_tolerance_w")
            );
            assert!(
                (number(exchange, "air_total_w") - (convective + radiative)).abs() > 0.1,
                "radiation exchanges with its independent reservoir, not transported air"
            );
        }
        let term = discretization(&qoi);
        assert_eq!(qoi.str_field("authority"), Some("estimated-candidate"));
        assert_eq!(term.str_field("state"), Some("interval"));
        assert_eq!(
            term.get("derivation").unwrap().str_field("method"),
            Some("adaptive-enrichment-order-one")
        );
        let half_width = 2.0
            * number(last, "estimated_change_k")
                .abs()
                .max(number(last, "measured_change_k").abs());
        assert_eq!(number(term, "upper_kelvin").to_bits(), half_width.to_bits());

        if mode == "vacuum" {
            let (replayed, replayed_field, _) = solved(&f, 101);
            assert_eq!(receipt.get("adaptive"), replayed.get("adaptive"));
            assert_eq!(field.get("temperature"), replayed_field.get("temperature"));
            // Resolve the identical first enrichment independently, without
            // paying for an unrelated third uniform rung in every composition.
            let first = &history[0];
            f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
            f.project.budgets.as_mut().unwrap().memory_bytes = (number(first, "probe_tets") as u64)
                .checked_mul(256)
                .unwrap();
            let (uniform, _, _) = solved(&f, 102);
            let rungs = uniform
                .get("ladder")
                .unwrap()
                .get("rungs")
                .unwrap()
                .as_array()
                .unwrap();
            assert_eq!(rungs.len(), 2);
            assert_eq!(first.f64_field("probe_tets"), rungs[1].f64_field("tets"));
            assert!(
                (number(first, "probe_t_max_k") - number(&rungs[1], "t_max_k")).abs() < 1e-6,
                "the independent uniform radiative/contact solve must reproduce the enriched field"
            );
        }
    }
}

#[test]
fn fixed_contact_adaptive_vacuum_uses_the_same_physical_radiation_comparison() {
    let mut f = Fixture::new(true, true);
    let setup = f
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    let bath = setup.radiation.as_ref().unwrap().surfaces[0]
        .reservoir_temperature
        .value;
    setup.boundaries[0].condition = B::HeatFlux {
        outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
    };
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = bath;
    envelope.ambient_hi.value = bath;
    f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
    f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
    // This uses the original scalar resistance pack, not temperature_pack:
    // vacuum adaptation must also work without a live contact source.
    let (receipt, field, _) = f.solve(0.13, None, 110);
    assert!(receipt.get("temperature_contact").is_none());
    assert_eq!(receipt.get("conjugate"), Some(&J::Null));
    let adaptive = receipt.get("adaptive").unwrap();
    assert_eq!(adaptive.str_field("status"), Some("observed-tolerance-met"));
    assert!(adaptive.f64_field("solved_meshes").unwrap() > 1.0);
    assert_eq!(adaptive.get("continuum_error_bound"), Some(&J::Bool(false)));
    let history = adaptive.get("history").unwrap().as_array().unwrap();
    assert!(!history.is_empty());
    for row in history {
        let balance = row.get("algebraic_balance").unwrap();
        assert_eq!(
            balance.str_field("comparison_model"),
            Some("area-mean-radiation-full-transpose")
        );
        assert_eq!(balance.get("air_reference_feedback"), Some(&J::Bool(false)));
        assert!(row.f64_field("interface_pairs").unwrap() > 0.0);
        assert_eq!(
            row.f64_field("probe_interface_pairs").unwrap(),
            4.0 * row.f64_field("interface_pairs").unwrap()
        );
        assert!(row.f64_field("absolute_contribution_sum_k").unwrap() > 0.0);
        assert!(row.f64_field("dual_residual").unwrap() <= 1e-9);
        let tolerance = 0.2 * (row.f64_field("t_max_k").unwrap() - bath).abs();
        assert_eq!(
            row.f64_field("tolerance_k").unwrap().to_bits(),
            tolerance.to_bits()
        );
    }
    assert!((maximum(&field) - history.last().unwrap().f64_field("t_max_k").unwrap()).abs() < 1e-9);
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(radiation.f64_field("convective_out_w"), Some(0.0));
    assert!((radiation.f64_field("nonlinear_radiative_out_w").unwrap() - 4.44).abs() < 1e-7);

    let first = &history[0];
    f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
    f.project.budgets.as_mut().unwrap().memory_bytes = (first.f64_field("probe_tets").unwrap()
        as u64)
        .checked_mul(256)
        .unwrap();
    let (uniform, _, _) = f.solve(0.13, None, 111);
    let rungs = uniform
        .get("ladder")
        .unwrap()
        .get("rungs")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(rungs.len(), 2);
    assert_eq!(first.f64_field("probe_tets"), rungs[1].f64_field("tets"));
    assert!(
        (first.f64_field("probe_t_max_k").unwrap() - rungs[1].f64_field("t_max_k").unwrap()).abs()
            < 1e-6
    );
}
