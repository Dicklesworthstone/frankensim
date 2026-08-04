//! End-to-end contract checks for deterministic Euler-disc campaign records.

use std::collections::BTreeMap;
use std::process::Command;

use fs_blake3::hash_domain;

#[test]
fn closed_campaign_emits_deterministic_controlled_trajectory_output() {
    let binary = env!("CARGO_BIN_EXE_euler_disc_campaign");
    let first = Command::new(binary)
        .arg("--closed-only")
        .output()
        .expect("closed campaign launches");
    let second = Command::new(binary)
        .arg("--closed-only")
        .output()
        .expect("closed campaign replay launches");
    assert!(
        first.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert_eq!(first.stdout, second.stdout);
    let output = String::from_utf8(first.stdout).expect("utf8");
    assert!(!output.contains("CONTOUR_FORCE_PER_NORMAL_FORCE"));
    assert!(output.contains("closed-time-evolving-profile-native-reduced-euler-disc"));
    assert!(output.contains("channel_work_j"));
    assert!(output.contains("\"relative_defect\""));
    assert!(output.contains("\"reimpact_count\""));
    assert!(output.contains("last_step_channel_work_j"));
    assert!(output.contains("precession_acceleration_rad_per_s2"));
    assert!(output.contains("\"defect_j\""));
    assert!(output.contains("model_disagreement"));
    assert!(output.contains("higher-fidelity transactional adapters"));
    assert!(!output.contains("shape_inertia_factor"));
    assert!(!output.contains("cone-inertia-cylinder-contact-surrogate"));
    assert!(!output.contains("\"duration_s\""));
    let v3_records: Vec<_> = output
        .lines()
        .filter(|line| line.contains("euler-disc-campaign-jsonl-v3"))
        .map(|line| JsonParser::parse(line).expect("v3 record is valid JSON"))
        .collect();
    let v3_roots: Vec<_> = v3_records.iter().map(object).collect();
    let closed: Vec<_> = v3_roots
        .iter()
        .copied()
        .filter(|fields| {
            fields.get("model").is_some_and(|_| {
                string(fields, "model") == "fs-mbd-profile-native-reduced-coupled-runner"
            })
        })
        .collect();
    assert_eq!(closed.len(), 11);
    assert!(
        closed
            .iter()
            .any(|fields| { string(fields, "scenario") == "closed-reduced-ring-equal-mass" })
    );
    assert!(
        closed
            .iter()
            .any(|fields| string(fields, "scenario") == "closed-reduced-solid-no-gas")
    );
    assert!(
        closed
            .iter()
            .any(|fields| { string(fields, "scenario") == "closed-reduced-solid-no-rolling" })
    );
    assert!(
        closed
            .iter()
            .any(|fields| { string(fields, "scenario") == "closed-reduced-symmetric-tapered" })
    );
    assert!(
        closed
            .iter()
            .any(|fields| { string(fields, "scenario") == "closed-reduced-solid-fillet-1mm" })
    );
    let inputs = object(closed[0].get("inputs").expect("closed inputs"));
    for key in [
        "mass_kg",
        "input_units",
        "radius_m",
        "thickness_m",
        "density_kg_m3",
        "gravity_m_per_s2",
        "transverse_inertia_kg_m2",
        "axial_inertia_kg_m2",
        "timestep_s",
        "maximum_steps",
        "initial_horizon_s",
        "maximum_horizon_s",
        "declared_final_horizon_s",
        "continuation_count",
        "terminal_inclination_rad",
        "reimpact_limit",
        "initial_inclination_rad",
        "initial_precession_rad_per_s",
        "initial_spin_rad_per_s",
        "sliding_friction_coefficient",
        "rolling_resistance_m",
        "base_effective_mass_kg",
        "base_stiffness_n_per_m",
        "base_damping_n_s_per_m",
        "contact_stiffness_n_per_m",
        "contact_damping_n_s_per_m",
        "gas_rotational_damping_n_m_s",
        "gas_translation_damping_n_s_per_m",
    ] {
        assert!(inputs.contains_key(key), "missing v3 input {key}");
    }
    let solid = closed
        .iter()
        .find(|fields| string(fields, "scenario") == "closed-reduced-solid")
        .expect("solid controlled case");
    let equal_mass_ring = closed
        .iter()
        .find(|fields| string(fields, "scenario") == "closed-reduced-ring-equal-mass")
        .expect("equal-mass ring controlled case");
    let solid_inputs = object(solid.get("inputs").expect("solid inputs"));
    let ring_inputs = object(equal_mass_ring.get("inputs").expect("ring inputs"));
    let solid_mass = number(solid_inputs, "mass_kg");
    let ring_mass = number(ring_inputs, "mass_kg");
    assert!((solid_mass - ring_mass).abs() / solid_mass < 1.0e-12);
    assert_ne!(
        number(solid_inputs, "axial_inertia_kg_m2"),
        number(ring_inputs, "axial_inertia_kg_m2")
    );
    let solid_outcome = object(solid.get("outcome").expect("solid outcome"));
    let ring_outcome = object(equal_mass_ring.get("outcome").expect("ring outcome"));
    assert_eq!(string(solid_outcome, "kind"), "right-censored");
    assert_eq!(
        string(ring_outcome, "kind"),
        "physical-terminal-inclination"
    );
    assert!(number(ring_outcome, "retained_time_s") < number(solid_outcome, "retained_time_s"));
    for fields in &closed {
        let profile = object(fields.get("profile").expect("resolved profile"));
        assert!(boolean(profile, "mass_and_support_same_chart"));
        let outcome = object(fields.get("outcome").expect("typed outcome"));
        let outcome_kind = string(outcome, "kind");
        assert!(matches!(
            outcome_kind,
            "physical-terminal-inclination" | "right-censored" | "numerical-refusal"
        ));
        assert_eq!(
            boolean(outcome, "observed_physical_terminal"),
            outcome_kind == "physical-terminal-inclination"
        );
        assert!(outcome.contains_key("numerical_refusal_reason"));
        let energy = object(fields.get("energy").expect("closed energy accounting"));
        assert!(number(energy, "initial_total_j") > 0.0);
        assert!(number(energy, "final_total_j") >= 0.0);
        let relative_defect = number(energy, "relative_defect");
        assert!(
            (0.0..=0.01).contains(&relative_defect),
            "unexpected relative energy defect for {}: {relative_defect}",
            string(fields, "scenario")
        );
    }
    let convergence = v3_roots
        .iter()
        .find(|fields| string(fields, "scenario") == "closed-reduced-solid-timestep-convergence")
        .expect("h/h2/h4 convergence record");
    assert_eq!(
        string(convergence, "observed_order"),
        "withheld-eventful-mode"
    );
    assert!(convergence.contains_key("fine_reference_qoi_linf"));
    assert!(convergence.contains_key("fine_reference_qoi"));
    let ranking = v3_roots
        .iter()
        .find(|fields| string(fields, "scenario") == "equal-mass-ring-vs-solid-ranking-convergence")
        .expect("censor-aware ranking convergence record");
    assert!(!boolean(ranking, "ordering_agreement"));
    assert!(!boolean(
        ranking,
        "ring_shorter_than_solid_bound_proven_at_all_rungs"
    ));
    assert!(!boolean(ranking, "ranking_numerically_supported"));
    let ranking_rungs = object(ranking.get("rungs").expect("ranking rungs"));
    let h2 = object(ranking_rungs.get("h2").expect("h2 ranking rung"));
    assert_eq!(
        string(h2, "censor_aware_ordering"),
        "ring-numerical-refusal"
    );
    assert_eq!(
        string(h2, "ring_numerical_refusal_reason"),
        "reimpact-limit-exceeded"
    );
    let calibration = v3_roots
        .iter()
        .find(|fields| string(fields, "scenario") == "physical-calibration-readiness")
        .expect("calibration readiness record");
    assert_eq!(string(calibration, "terminal"), "no-data");
    assert!(!boolean(calibration, "synthetic_substitution"));
    let manifest = v3_roots
        .iter()
        .find(|fields| string(fields, "scenario") == "campaign-complete")
        .expect("v3 manifest");
    assert_eq!(number(manifest, "record_count"), (closed.len() + 3) as f64);
}

#[derive(Debug)]
enum JsonValue {
    Number(f64),
    String(String),
    Boolean(bool),
    Object(BTreeMap<String, JsonValue>),
}

/// Strict JSON subset parser for this ASCII-only numeric producer; it is not a
/// general JSON implementation.
struct JsonParser<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> JsonParser<'a> {
    fn parse(input: &'a str) -> Result<JsonValue, String> {
        let mut parser = Self {
            bytes: input.as_bytes(),
            position: 0,
        };
        let value = parser.value()?;
        parser.whitespace();
        if parser.position != parser.bytes.len() {
            return Err("trailing JSON bytes".to_owned());
        }
        Ok(value)
    }

    fn value(&mut self) -> Result<JsonValue, String> {
        self.whitespace();
        match self.peek() {
            Some(b'{') => self.object(),
            Some(b'"') => self.string().map(JsonValue::String),
            Some(b'-' | b'0'..=b'9') => self.number().map(JsonValue::Number),
            Some(b't') => self.literal(b"true", JsonValue::Boolean(true)),
            Some(b'f') => self.literal(b"false", JsonValue::Boolean(false)),
            _ => Err("expected JSON object, string, or number".to_owned()),
        }
    }

    fn object(&mut self) -> Result<JsonValue, String> {
        self.expect(b'{')?;
        self.whitespace();
        let mut fields = BTreeMap::new();
        if self.consume(b'}') {
            return Ok(JsonValue::Object(fields));
        }
        loop {
            self.whitespace();
            let key = self.string()?;
            self.whitespace();
            self.expect(b':')?;
            let value = self.value()?;
            if fields.insert(key, value).is_some() {
                return Err("duplicate JSON key".to_owned());
            }
            self.whitespace();
            if self.consume(b'}') {
                return Ok(JsonValue::Object(fields));
            }
            self.expect(b',')?;
        }
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut output = String::new();
        while let Some(byte) = self.next() {
            match byte {
                b'"' => return Ok(output),
                b'\\' => match self.next() {
                    Some(b'"') => output.push('"'),
                    Some(b'\\') => output.push('\\'),
                    Some(b'/') => output.push('/'),
                    Some(b'b') => output.push('\u{0008}'),
                    Some(b'f') => output.push('\u{000c}'),
                    Some(b'n') => output.push('\n'),
                    Some(b'r') => output.push('\r'),
                    Some(b't') => output.push('\t'),
                    _ => return Err("unsupported JSON escape".to_owned()),
                },
                0..=0x1f => return Err("unescaped control character".to_owned()),
                byte => output.push(byte as char),
            }
        }
        Err("unterminated JSON string".to_owned())
    }

    fn number(&mut self) -> Result<f64, String> {
        let start = self.position;
        while matches!(
            self.peek(),
            Some(b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9')
        ) {
            self.position += 1;
        }
        std::str::from_utf8(&self.bytes[start..self.position])
            .map_err(|_| "non-UTF8 JSON number".to_owned())?
            .parse()
            .map_err(|_| "invalid JSON number".to_owned())
    }

    fn literal(&mut self, expected: &[u8], value: JsonValue) -> Result<JsonValue, String> {
        let end = self.position.saturating_add(expected.len());
        if self.bytes.get(self.position..end) == Some(expected) {
            self.position = end;
            Ok(value)
        } else {
            Err("invalid JSON literal".to_owned())
        }
    }

    fn whitespace(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.position += 1;
        }
    }

    fn expect(&mut self, expected: u8) -> Result<(), String> {
        if self.next() == Some(expected) {
            Ok(())
        } else {
            Err(format!("expected JSON byte {}", expected as char))
        }
    }

    fn consume(&mut self, expected: u8) -> bool {
        if self.peek() == Some(expected) {
            self.position += 1;
            true
        } else {
            false
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.position).copied()
    }

    fn next(&mut self) -> Option<u8> {
        let byte = self.peek()?;
        self.position += 1;
        Some(byte)
    }
}

fn object(value: &JsonValue) -> &BTreeMap<String, JsonValue> {
    let JsonValue::Object(fields) = value else {
        panic!("expected JSON object");
    };
    fields
}

fn string<'a>(fields: &'a BTreeMap<String, JsonValue>, key: &str) -> &'a str {
    let Some(JsonValue::String(value)) = fields.get(key) else {
        panic!("expected JSON string {key}");
    };
    value
}

fn number(fields: &BTreeMap<String, JsonValue>, key: &str) -> f64 {
    let Some(JsonValue::Number(value)) = fields.get(key) else {
        panic!("expected JSON number {key}");
    };
    *value
}

fn boolean(fields: &BTreeMap<String, JsonValue>, key: &str) -> bool {
    let Some(JsonValue::Boolean(value)) = fields.get(key) else {
        panic!("expected JSON boolean {key}");
    };
    *value
}

#[test]
fn campaign_executable_emits_deterministic_substantive_production_records() {
    let binary = env!("CARGO_BIN_EXE_euler_disc_campaign");
    let first = Command::new(binary).output().expect("campaign launches");
    assert!(
        first.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    let second = Command::new(binary).output().expect("campaign relaunches");
    assert!(second.status.success());
    assert_eq!(
        first.stdout, second.stdout,
        "campaign output is deterministic"
    );

    let stdout = String::from_utf8(first.stdout).expect("JSONL is UTF-8");
    let lines: Vec<_> = stdout.lines().collect();
    let records: Vec<_> = lines
        .iter()
        .map(|line| JsonParser::parse(line).expect("record is valid JSON"))
        .collect();
    assert!(
        records.len() >= 17,
        "legacy records plus closed trajectory records"
    );
    let roots: Vec<_> = records[..10].iter().map(object).collect();
    let scenarios: Vec<_> = roots
        .iter()
        .map(|fields| string(fields, "scenario"))
        .collect();
    assert_eq!(
        scenarios,
        [
            "geometry-sharp-squat-disc",
            "geometry-filleted-squat-disc",
            "conservative-steady-oracle",
            "dynamic-unilateral-contact",
            "reduced-flexible-base",
            "contour-only-decay",
            "boundary-layer-only-decay",
            "combined-decay",
            "reduced-exterior-wrench-passivity",
            "campaign-complete",
        ]
    );
    for fields in &roots {
        assert_eq!(string(fields, "schema"), "euler-disc-campaign-jsonl-v1");
        for key in [
            "model",
            "source",
            "authority",
            "units",
            "budget",
            "terminal",
            "residual",
            "no_claim",
        ] {
            assert!(fields.contains_key(key), "missing {key}");
        }
        assert!(string(fields, "no_claim").contains("no-physical-validation"));
        if string(fields, "scenario") == "campaign-complete" {
            assert_eq!(
                string(fields, "campaign_seed_u64_dec"),
                "4995983222254027085"
            );
        } else {
            assert_eq!(
                string(fields, "campaign_seed_manifest_ref"),
                "campaign-complete"
            );
        }
    }

    let sharp_inputs = object(roots[0].get("inputs").expect("sharp inputs"));
    let filleted_inputs = object(roots[1].get("inputs").expect("filleted inputs"));
    assert_eq!(number(sharp_inputs, "radius_m"), 0.038);
    assert_eq!(number(sharp_inputs, "thickness_m"), 0.006);
    assert_eq!(number(sharp_inputs, "edge_radius_m"), 0.0);
    assert_eq!(number(filleted_inputs, "edge_radius_m"), 0.001);
    let sharp_budget = object(roots[0].get("budget").expect("sharp budget"));
    assert_eq!(string(sharp_budget, "seed_u64_dec"), "4995983222254027085");
    assert!(
        number(roots[0], "mass_kg") > number(roots[1], "mass_kg"),
        "the actual fillet changes the geometry-derived mass"
    );

    let oracle_inputs = object(roots[2].get("inputs").expect("oracle inputs"));
    assert!(number(oracle_inputs, "thickness_m") < 0.01 * number(oracle_inputs, "radius_m"));
    assert_ne!(
        number(oracle_inputs, "thickness_m"),
        number(sharp_inputs, "thickness_m")
    );
    let oracle_residual = object(roots[2].get("residual").expect("oracle residual"));
    assert!(number(oracle_residual, "precession_s_inv2").abs() < 1.0e-6);

    let contact_inputs = object(roots[3].get("inputs").expect("contact inputs"));
    let contact_residual = object(roots[3].get("residual").expect("contact residual"));
    let contact_loads = object(roots[3].get("loads_n").expect("contact loads"));
    assert_eq!(
        string(contact_inputs, "primary_case"),
        "fillet-fixed-density"
    );
    assert_eq!(number(contact_inputs, "static_friction_coefficient"), 100.0);
    assert_eq!(
        string(contact_inputs, "interface_system_id"),
        "campaign/squat-disc->plane"
    );
    assert_eq!(
        string(contact_inputs, "interface_history_id"),
        "campaign/shared-matched-profile-contact-history-v1"
    );
    assert_eq!(
        string(contact_inputs, "dry_interface_authority"),
        "caller-declared-dry-interface"
    );
    assert_eq!(string(roots[3], "authority"), "numerical-reference-only");
    let contact_normal_reaction_n = number(contact_loads, "normal_reaction_mean_n");
    assert!(contact_normal_reaction_n > 0.0);
    assert!(
        number(contact_residual, "energy_j").abs()
            <= number(contact_residual, "energy_residual_limit_j")
    );
    assert!(
        number(contact_residual, "sum_abs_energy_j")
            <= number(contact_residual, "sum_abs_energy_limit_j")
    );
    assert!(
        number(contact_residual, "max_abs_energy_j")
            <= number(contact_residual, "max_abs_energy_limit_j")
    );
    assert!(
        number(contact_residual, "refinement_position_m").abs()
            <= number(contact_residual, "refinement_position_limit_m")
    );
    assert!(
        number(contact_residual, "refinement_energy_j").abs()
            <= number(contact_residual, "refinement_energy_limit_j")
    );
    let top_energy_scale_j = number(contact_residual, "energy_scale_j");
    assert_eq!(
        number(contact_residual, "energy_residual_limit_j"),
        top_energy_scale_j * 1.0e-4
    );
    assert_eq!(
        number(contact_residual, "sum_abs_energy_limit_j"),
        top_energy_scale_j * 1.0e-3
    );
    assert_eq!(
        number(contact_residual, "max_abs_energy_limit_j"),
        top_energy_scale_j * 1.0e-4
    );
    assert_eq!(
        number(contact_residual, "refinement_position_limit_m"),
        number(contact_inputs, "radius_m") * 1.0e-1
    );
    assert_eq!(
        number(contact_residual, "refinement_energy_limit_j"),
        top_energy_scale_j * 1.0e-1
    );
    let cases = object(roots[3].get("cases").expect("contact cases"));
    let sharp_case = object(cases.get("sharp_fixed_density").expect("sharp case"));
    let fillet_case = object(cases.get("fillet_fixed_density").expect("fillet case"));
    let equal_mass_case = object(cases.get("fillet_equal_mass").expect("equal mass case"));
    assert_eq!(
        contact_normal_reaction_n,
        number(fillet_case, "mean_reaction_n"),
        "the propagated contact load is the primary case mean reaction"
    );
    assert_eq!(
        number(contact_residual, "energy_j"),
        number(fillet_case, "summed_mechanical_balance_residual_j")
    );
    assert_eq!(
        number(contact_residual, "sum_abs_energy_j"),
        number(fillet_case, "sum_abs_mechanical_balance_residual_j")
    );
    assert_eq!(
        number(contact_residual, "max_abs_energy_j"),
        number(fillet_case, "max_abs_mechanical_balance_residual_j")
    );
    assert_eq!(
        number(contact_residual, "energy_scale_j"),
        number(fillet_case, "mechanical_energy_scale_j")
    );
    for case_fields in [sharp_case, fillet_case, equal_mass_case] {
        assert_eq!(string(case_fields, "terminal"), "horizon-reached");
        assert!(number(case_fields, "mean_reaction_n") > 0.0);
        assert!(number(case_fields, "final_reaction_n") > 0.0);
        assert!(number(case_fields, "peak_reaction_n") > 0.0);
        assert!(number(case_fields, "max_required_static_mu") >= 0.0);
        assert!(number(case_fields, "initial_material_contact_speed_m_per_s") < 1.0e-9);
        assert!(number(case_fields, "inertia_transverse_kg_m2") > 0.0);
        assert!(number(case_fields, "inertia_axial_kg_m2") > 0.0);
        assert!(number(case_fields, "center_height_m") > 0.0);
        assert!(number(case_fields, "mechanical_energy_scale_j") > 0.0);
        assert!(number(case_fields, "sum_abs_mechanical_balance_residual_j") >= 0.0);
        assert!(number(case_fields, "max_abs_mechanical_balance_residual_j") >= 0.0);
        let energy_scale_j = number(case_fields, "mechanical_energy_scale_j");
        assert!(
            number(case_fields, "summed_mechanical_balance_residual_j").abs()
                <= number(case_fields, "energy_residual_limit_j")
        );
        assert!(
            number(case_fields, "sum_abs_mechanical_balance_residual_j")
                <= number(case_fields, "sum_abs_energy_limit_j")
        );
        assert!(
            number(case_fields, "max_abs_mechanical_balance_residual_j")
                <= number(case_fields, "max_abs_energy_limit_j")
        );
        assert_eq!(
            number(case_fields, "energy_residual_limit_j"),
            energy_scale_j * 1.0e-4
        );
        assert_eq!(
            number(case_fields, "sum_abs_energy_limit_j"),
            energy_scale_j * 1.0e-3
        );
        assert_eq!(
            number(case_fields, "max_abs_energy_limit_j"),
            energy_scale_j * 1.0e-4
        );
        let refinement = object(
            case_fields
                .get("quarter_step_refinement")
                .expect("quarter-step refinement receipt"),
        );
        for (coarse, fine, flag) in [
            (
                "coarse_position_error_m",
                "fine_position_error_m",
                "position_refinement_improved",
            ),
            (
                "coarse_linear_momentum_error_kg_m_per_s",
                "fine_linear_momentum_error_kg_m_per_s",
                "linear_momentum_refinement_improved",
            ),
            (
                "coarse_angular_momentum_error_kg_m2_per_s",
                "fine_angular_momentum_error_kg_m2_per_s",
                "angular_momentum_refinement_improved",
            ),
            (
                "coarse_orientation_error_rad",
                "fine_orientation_error_rad",
                "orientation_refinement_improved",
            ),
            (
                "coarse_energy_error_j",
                "fine_energy_error_j",
                "energy_refinement_improved",
            ),
        ] {
            assert!(number(refinement, coarse).is_finite());
            assert!(number(refinement, fine).is_finite());
            assert_eq!(
                boolean(refinement, flag),
                number(refinement, fine) <= number(refinement, coarse)
            );
            assert!(
                boolean(refinement, flag),
                "quarter-step {flag} did not improve"
            );
        }
    }
    assert!(number(sharp_case, "mass_kg") > number(fillet_case, "mass_kg"));
    assert!((number(sharp_case, "mass_kg") - number(equal_mass_case, "mass_kg")).abs() < 1.0e-12);
    assert!(number(equal_mass_case, "density_kg_m3") > number(fillet_case, "density_kg_m3"));
    let sharp_support = object(
        sharp_case
            .get("support_vector_world_m")
            .expect("sharp support vector"),
    );
    let fillet_support = object(
        fillet_case
            .get("support_vector_world_m")
            .expect("fillet support vector"),
    );
    let support_difference = ["x", "y", "z"]
        .into_iter()
        .map(|axis| (number(sharp_support, axis) - number(fillet_support, axis)).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        support_difference > 1.0e-9,
        "sharp and fillet supports differ"
    );
    let contact_deltas = object(roots[3].get("deltas").expect("contact deltas"));
    for (delta, left, right) in [
        (
            object(
                contact_deltas
                    .get("fillet_fixed_density_minus_sharp_fixed_density")
                    .expect("fixed density delta"),
            ),
            fillet_case,
            sharp_case,
        ),
        (
            object(
                contact_deltas
                    .get("fillet_equal_mass_minus_sharp_fixed_density")
                    .expect("equal mass delta"),
            ),
            equal_mass_case,
            sharp_case,
        ),
    ] {
        for key in [
            "mean_reaction_n",
            "final_reaction_n",
            "peak_reaction_n",
            "max_required_static_mu",
            "summed_mechanical_balance_residual_j",
            "sum_abs_mechanical_balance_residual_j",
            "max_abs_mechanical_balance_residual_j",
            "mechanical_energy_scale_j",
            "refinement_position_m",
            "refinement_energy_j",
        ] {
            assert!(number(delta, key).is_finite());
            assert!(
                (number(delta, key) - (number(left, key) - number(right, key))).abs() < 1.0e-12,
                "delta {key} did not match its two retained cases"
            );
        }
    }
    let contact_no_claim = string(roots[3], "no_claim");
    assert!(contact_no_claim.contains("one-way-not-closed-coupling"));
    assert!(contact_no_claim.contains("not-outcome-ranking"));
    assert_eq!(
        string(roots[3], "source"),
        "frep/matched-profile-contact-comparison"
    );

    let base_inputs = object(roots[4].get("inputs").expect("base inputs"));
    let base_residual = object(roots[4].get("residual").expect("base residual"));
    let base_loads = object(roots[4].get("loads_n").expect("base loads"));
    let base_powers = object(roots[4].get("powers_w").expect("base powers"));
    assert_eq!(
        number(base_inputs, "normal_force_n"),
        contact_normal_reaction_n
    );
    assert!(number(base_residual, "energy_j").abs() < 2.0e-6);
    assert!(number(base_residual, "refinement_displacement_m") < 1.0e-6);
    assert!(number(base_loads, "modal_force_n").is_finite());
    assert!(base_powers.is_empty(), "a force is not a power channel");

    let contour_work = object(roots[5].get("work_j").expect("contour work"));
    let boundary_work = object(roots[6].get("work_j").expect("boundary work"));
    let combined_work = object(roots[7].get("work_j").expect("combined work"));
    assert!(number(contour_work, "dry_contour") > 0.0);
    assert_eq!(number(contour_work, "bildsten_boundary_layer"), 0.0);
    assert_eq!(number(boundary_work, "dry_contour"), 0.0);
    assert!(number(boundary_work, "bildsten_boundary_layer") > 0.0);
    assert!(number(combined_work, "dry_contour") > 0.0);
    assert!(number(combined_work, "bildsten_boundary_layer") > 0.0);

    for record in &roots[5..8] {
        let inputs = object(record.get("inputs").expect("decay inputs"));
        let final_state = object(record.get("final").expect("decay final state"));
        let crossover = object(record.get("crossover").expect("decay crossover"));
        assert_eq!(string(record, "terminal"), "validity-cutoff");
        assert!(number(final_state, "time_s") > 0.0);
        assert!(number(final_state, "theta_rad") > 0.0);
        assert!(number(final_state, "omega_rad_s") > 0.0);
        assert!(number(final_state, "energy_j") > 0.0);
        match string(crossover, "class") {
            "encoded-power-law-crossover" => {
                assert!(number(crossover, "theta_rad") > 0.0);
                if crossover.contains_key("time_s") {
                    assert!(number(crossover, "time_s") >= 0.0);
                } else {
                    assert_eq!(
                        string(crossover, "time_status"),
                        "outside-retained-trajectory"
                    );
                }
            }
            "none" | "not-comparable" => {
                assert!(!crossover.contains_key("theta_rad"));
                assert!(!crossover.contains_key("time_s"));
                assert!(crossover.contains_key("theta_status"));
                assert!(crossover.contains_key("time_status"));
            }
            other => panic!("unexpected crossover class {other}"),
        }
        assert_eq!(number(inputs, "mass_kg"), number(roots[1], "mass_kg"));
        let residual = object(record.get("residual").expect("decay residual"));
        assert!(number(residual, "energy_scale_j") > 0.0);
        assert!(number(residual, "energy_j").abs() <= number(residual, "energy_residual_limit_j"));
        assert!(
            number(residual, "refinement_time_s").abs()
                <= number(residual, "refinement_time_limit_s")
        );
        assert!(
            number(residual, "refinement_work_j").abs()
                <= number(residual, "refinement_work_limit_j")
        );
    }
    let contour_inputs = object(roots[5].get("inputs").expect("contour inputs"));
    let combined_inputs = object(roots[7].get("inputs").expect("combined inputs"));
    assert_eq!(
        number(contour_inputs, "dry_normal_force_n"),
        contact_normal_reaction_n
    );
    assert_eq!(
        number(combined_inputs, "dry_normal_force_n"),
        contact_normal_reaction_n
    );
    assert_eq!(
        number(contour_inputs, "contour_force_n"),
        contact_normal_reaction_n * number(contour_inputs, "contour_force_per_normal_force")
    );
    assert_eq!(
        number(combined_inputs, "contour_force_n"),
        contact_normal_reaction_n * number(combined_inputs, "contour_force_per_normal_force")
    );

    let combined_sources = object(roots[7].get("sources").expect("combined sources"));
    assert_ne!(string(combined_sources, "dry"), "none");
    assert_ne!(string(combined_sources, "bildsten"), "none");
    let flux_powers = object(roots[8].get("powers_w").expect("flux powers"));
    let flux_residual = object(roots[8].get("residual").expect("flux residual"));
    assert!(number(flux_powers, "relative_dissipation") > 0.0);
    assert!(number(flux_residual, "passivity_w") <= 1.0e-10);
    let flux_inputs = object(roots[8].get("inputs").expect("flux inputs"));
    assert_eq!(number(flux_inputs, "radius_m"), 0.038);
    assert_eq!(number(flux_inputs, "thickness_m"), 0.006);
    assert!(number(flux_inputs, "snapshot_angular_speed_rad_s") > 0.0);
    assert_eq!(string(flux_inputs, "gas_source_id"), "campaign/air-card-v1");
    assert_eq!(
        string(flux_inputs, "roughness_source_id"),
        "campaign/exterior-roughness-card-v1"
    );
    assert_eq!(
        string(flux_inputs, "correlation_authority"),
        "caller-declared-screening"
    );
    assert_eq!(
        string(roots[8], "source"),
        "campaign/caller-declared-screening-correlation-v1"
    );

    let manifest_budget = object(roots[9].get("budget").expect("manifest budget"));
    assert_eq!(number(manifest_budget, "record_count"), 9.0);
    assert_eq!(number(manifest_budget, "declared_cx_poll_quota"), 10_000.0);
    assert_eq!(number(manifest_budget, "declared_cx_cost_quota"), 100_000.0);
    assert_eq!(
        string(roots[9], "digest_domain"),
        "org.frankensim.euler-disc-campaign-jsonl.v1"
    );
    assert_eq!(
        string(roots[9], "digest_scope"),
        "preceding-data-records-LF-joined-no-trailing-LF"
    );
    assert_eq!(string(roots[9], "digest_blake3").len(), 64);
    let exact_payload = lines[..9].join("\n");
    assert_eq!(
        string(roots[9], "digest_blake3"),
        hash_domain(
            "org.frankensim.euler-disc-campaign-jsonl.v1",
            exact_payload.as_bytes()
        )
        .to_hex()
    );
}
