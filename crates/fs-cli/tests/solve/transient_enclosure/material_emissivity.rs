//! Original source emissivity changes emission and reflected enclosure exchange.
//! Exact triangle polynomial integrals independently audit the physical reports.

use super::adjoint::evaluate;
use super::*;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};

const SOURCE_KNOTS: [(f64, f64); 3] = [(200.0, 0.11), (280.0, 0.35), (450.0, 0.86)];
const OPEN_VIEW: f64 = 0.35;
const OPEN_BATH_K: f64 = 330.0;
const AREA: f64 = 0.5;

struct MaterialFixture {
    cards: CardPackSet,
    card: String,
    state: String,
    claim: String,
}

fn material(knots: &[(f64, f64)], specimen: &str) -> MaterialFixture {
    let raw = format!("synthetic original material; specimen={specimen}; epsilon(T)={knots:?}");
    let source = hash_bytes(raw.as_bytes());
    let provenance = || Provenance {
        source: "synthetic combined conductivity and emissivity; not experimental data".into(),
        license: "CC-BY-4.0; redistribution permitted".into(),
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: specimen.into(),
            method: "authored scalar conductivity and original emissivity ordinates".into(),
            artifact: source,
            caveats: "numerical verification only; no measured material qualification".into(),
            provenance: provenance(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal-conductivity", CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 6.0,
                dims: CONDUCTIVITY_DIMS,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let claim = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::SURFACE_EMISSIVITY_PROPERTY,
                fs_qty::Dims::NONE,
            ),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: dims::TEMPERATURE,
                knots: knots.to_vec(),
                dims: fs_qty::Dims::NONE,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with(
                "T",
                knots.first().unwrap().0,
                knots.last().unwrap().0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap()
        .0
        .to_hex();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "native-enclosure-emissivity".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        NormalizedPack::new(
            "native-enclosure-emissivity",
            "synthetic-v1",
            source,
            "CC-BY-4.0; redistribution permitted",
            claims,
            Vec::new(),
            Vec::new(),
        )
        .unwrap(),
    )
    .unwrap();
    let cards = CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "synthetic-original-material.fsmcdpk",
        pack.to_bytes(),
    )])
    .unwrap();
    let card = cards.materials()[0].card().to_hex();
    let state = cards.materials()[0].identity().to_string();
    MaterialFixture {
        cards,
        card,
        state,
        claim,
    }
}

fn fixture(material: &MaterialFixture, open: bool) -> ProjectSpec {
    let mut spec = project(&material.card);
    rebind_to(&mut spec, &material.cards);
    let binding = &mut spec.materials.as_mut().unwrap()[0];
    binding.claim = None;
    binding.source =
        "synthetic original storage material with evolving enclosure emissivity".into();
    close(binding.temp_lo.value, 233.15, 0.0);
    close(binding.temp_hi.value, 398.15, 0.0);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    for (index, boundary) in setup.boundaries.iter_mut().enumerate() {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(
                if open {
                    0.0
                } else if index == 0 {
                    -12.0
                } else {
                    2.0
                },
                dims::HEAT_FLUX,
            ),
        };
    }
    let declaration = setup.enclosure_radiation.as_mut().unwrap();
    declaration.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    declaration.emissivity_response =
        Some(fs_project::RadiationEmissivityResponse::LocalTemperature);
    declaration.max_iterations = 1;
    for surface in &mut declaration.surfaces {
        surface.claim = Some(material.claim.clone());
    }
    if open {
        declaration.ambient_temperature = Some(QtyAny::new(OPEN_BATH_K, dims::TEMPERATURE));
        declaration.view_factors = vec![vec![0.0, OPEN_VIEW], vec![OPEN_VIEW, 0.0]];
    }
    let policy = time(&mut spec);
    policy.initial_temperature = fs_project::TransientInitialTemperature::FromRegionalEnthalpy;
    policy.enthalpy.as_mut().unwrap().materials[0].material_card_identity =
        fs_blake3::ContentHash::from_hex(&material.card).unwrap();
    spec
}

fn epsilon(temperature: f64) -> f64 {
    let pair = SOURCE_KNOTS
        .windows(2)
        .find(|pair| pair[0].0 <= temperature && temperature <= pair[1].0)
        .expect("independent emissivity oracle refuses extrapolation");
    pair[0].1 + (temperature - pair[0].0) * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)
}

/// Exact barycentric polynomial averages on a triangle. The normalized
/// integral of lambda_0^a lambda_1^b lambda_2^c is 2 a! b! c!/(a+b+c+2)!.
/// Expanding T_h^4 leaves 15 monomials with weight 1/15; multiplying by the
/// independent P1 epsilon field gives the weighted sum below with weight 1/105.
/// This uses neither the production quadrature nor an emission-at-mean-T law.
fn emission_moments(temperatures: [f64; 3], emissivities: [f64; 3]) -> (f64, f64) {
    let mut black = 0.0;
    let mut emitted = 0.0;
    for a in 0..=4 {
        for b in 0..=4 - a {
            let c = 4 - a - b;
            let term = temperatures[0].powi(a) * temperatures[1].powi(b) * temperatures[2].powi(c);
            black += term / 15.0;
            emitted += term
                * (emissivities[0] * f64::from(a + 1)
                    + emissivities[1] * f64::from(b + 1)
                    + emissivities[2] * f64::from(c + 1))
                / 105.0;
        }
    }
    (SIGMA * black, SIGMA * emitted)
}

fn check_exchange(report: &J, material: &MaterialFixture, curve: &str, open: bool) -> [f64; 2] {
    assert_eq!(
        report.str_field("model"),
        Some(if open {
            "pointwise-p1-material-open-gray-diffuse-enclosure"
        } else {
            "pointwise-p1-material-closed-gray-diffuse-enclosure"
        })
    );
    assert_eq!(report.str_field("spatial_model"), Some("pointwise"));
    assert_eq!(
        report.str_field("emissivity_response"),
        Some("local-temperature")
    );
    assert!(
        fs_blake3::ContentHash::from_hex(report.str_field("view_factor_identity").unwrap())
            .is_some()
    );
    let surfaces = report.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 2);
    let mut means = [0.0; 2];
    let mut emissions = [0.0; 2];
    for (index, surface) in surfaces.iter().enumerate() {
        assert_eq!(surface.str_field("card"), Some(material.card.as_str()));
        assert_eq!(
            surface.str_field("material_state"),
            Some(material.state.as_str())
        );
        assert_eq!(
            surface.str_field("selected_claim"),
            Some(material.claim.as_str())
        );
        close(n(surface, "area_m2"), AREA, 1e-14);
        close(n(surface, "reference_emissivity"), epsilon(300.0), 1e-15);
        let states = surface.get("material_states").unwrap().as_array().unwrap();
        assert_eq!(
            states.len(),
            3,
            "each imported patch is one original triangle"
        );
        let mut temperatures = [0.0; 3];
        let mut emissivities = [0.0; 3];
        for (node, state) in states.iter().enumerate() {
            let temperature = n(state, "temperature_k");
            temperatures[node] = temperature;
            emissivities[node] = epsilon(temperature);
            close(n(state, "emissivity"), emissivities[node], 1e-12);
            close(
                n(state, "specific_enthalpy_j_kg"),
                CP * (temperature - INITIAL_K),
                3e-7,
            );
            close(n(state, "liquid_mass_fraction"), 0.0, 0.0);
            assert_eq!(state.str_field("phase_curve_identity"), Some(curve));
            assert_eq!(
                state.str_field("selected_claim"),
                Some(material.claim.as_str())
            );
            for key in ["material_state_identity", "usage_receipt"] {
                assert!(fs_blake3::ContentHash::from_hex(state.str_field(key).unwrap()).is_some());
            }
        }
        means[index] = emissivities.iter().sum::<f64>() / 3.0;
        let (black, emitted) = emission_moments(temperatures, emissivities);
        emissions[index] = emitted;
        close(
            n(surface, "mean_temperature_k"),
            temperatures.iter().sum::<f64>() / 3.0,
            2e-11,
        );
        close(n(surface, "material_mean_emissivity"), means[index], 2e-14);
        close(n(surface, "mean_blackbody_emission_w_m2"), black, 2e-8);
    }
    // J_i = <epsilon B>_i + (1-<epsilon>_i) [F_ij J_j + F_i,b B_b].
    // The two-by-two inverse includes all reflected paths and rebuilds its
    // matrix from the current source state, independently of the native owner.
    let view = if open { OPEN_VIEW } else { 1.0 };
    let bath_black = if open {
        SIGMA * OPEN_BATH_K.powi(4)
    } else {
        0.0
    };
    let reflected = means.map(|value| (1.0 - value) * view);
    let rhs =
        [0, 1].map(|index| emissions[index] + (1.0 - means[index]) * (1.0 - view) * bath_black);
    let determinant = 1.0 - reflected[0] * reflected[1];
    let radiosities = [
        (rhs[0] + reflected[0] * rhs[1]) / determinant,
        (rhs[1] + reflected[1] * rhs[0]) / determinant,
    ];
    let mut heat = [0.0; 2];
    for index in 0..2 {
        let irradiation = view * radiosities[1 - index] + (1.0 - view) * bath_black;
        heat[index] = AREA * (emissions[index] - means[index] * irradiation);
        close(
            n(&surfaces[index], "radiosity_w_m2"),
            radiosities[index],
            2e-8,
        );
        close(n(&surfaces[index], "irradiation_w_m2"), irradiation, 2e-8);
        close(n(&surfaces[index], "outward_heat_w"), heat[index], 2e-8);
    }
    close(n(report, "radiative_out_w"), heat.iter().sum(), 5e-8);
    close(
        n(report, "absolute_patch_heat_w"),
        heat.iter().map(|x| x.abs()).sum(),
        5e-8,
    );
    close(n(report, "convective_out_w"), 0.0, 0.0);
    close(n(report, "enclosure_energy_closure_w"), 0.0, 1e-6);
    if open {
        // Reciprocal closure gives one black row with equal view to both
        // equal-area patches. Its independent heat is outside solid storage.
        let bath_area = 2.0 * AREA * (1.0 - view);
        let bath_irradiation = 0.5 * (radiosities[0] + radiosities[1]);
        let bath_absorbed = bath_area * (bath_irradiation - bath_black);
        close(n(report, "ambient_temperature_k"), OPEN_BATH_K, 0.0);
        close(n(report, "ambient_effective_area_m2"), bath_area, 1e-14);
        close(n(report, "ambient_radiosity_w_m2"), bath_black, 2e-8);
        close(
            n(report, "ambient_irradiation_w_m2"),
            bath_irradiation,
            2e-8,
        );
        close(n(report, "ambient_absorbed_w"), bath_absorbed, 5e-8);
        close(heat.iter().sum(), bath_absorbed, 1e-6);
    } else {
        assert!(report.get("ambient_temperature_k").is_none());
        close(heat.iter().sum(), 0.0, 1e-6);
    }
    heat
}

#[test]
fn g1_local_enclosure_emissivity_rebuilds_reflection_and_closes_physical_heat_accounts() {
    let material = material(&SOURCE_KNOTS, "original-enclosure-specimen");
    for open in [false, true] {
        let spec = fixture(&material, open);
        assert!(
            print_sexpr(&spec)
                .unwrap()
                .contains(":emissivity-response local-temperature")
        );
        let solved = run(&spec, &material.cards);
        let receipt = J::parse(&solved.conduction).unwrap();
        let history = receipt.get("transient").unwrap();
        let source = &history.get("materials").unwrap().as_array().unwrap()[0];
        assert_eq!(
            source.str_field("material_card_identity"),
            Some(material.card.as_str())
        );
        let curve = source.str_field("phase_curve_identity").unwrap();
        let mut nonzero_exchange = false;
        let mut nonuniform_source = false;
        let mut fine_stored = 0.0;
        for grid in ["coarse", "fine"] {
            for row in history.get(grid).unwrap().as_array().unwrap() {
                let radiation = row.get("radiation").unwrap();
                let exchange = radiation.get("exchange").unwrap();
                let heat = check_exchange(exchange, &material, curve, open);
                nonzero_exchange |= heat.iter().any(|q| q.abs() > 1e-5);
                for surface in exchange.get("surfaces").unwrap().as_array().unwrap() {
                    nonuniform_source |=
                        n(surface, "maximum_emissivity") - n(surface, "minimum_emissivity") > 1e-6;
                }
                let radiative: f64 = heat.iter().sum();
                let neumann = if open { 0.0 } else { -5.0 };
                let net = n(row, "source_w") - neumann - radiative;
                close(n(radiation, "solid_solves"), 1.0, 0.0);
                close(n(radiation, "prescribed_neumann_out_w"), neumann, 1e-12);
                close(n(radiation, "convective_out_w"), 0.0, 0.0);
                close(n(radiation, "physical_dirichlet_in_w"), 0.0, 0.0);
                close(n(radiation, "nonlinear_radiative_out_w"), radiative, 5e-8);
                close(n(row, "net_input_w"), net, 5e-8);
                close(n(row, "stored_energy_change_j"), n(row, "dt_s") * net, 2e-6);
                assert!(
                    n(radiation, "physical_residual_j")
                        <= n(radiation, "physical_residual_threshold_j")
                );
                assert!(
                    n(radiation, "physical_energy_residual_j").abs()
                        <= n(history, "energy_tolerance_j")
                );
                if grid == "fine" {
                    fine_stored += n(row, "stored_energy_change_j");
                }
            }
        }
        assert!(
            nonzero_exchange,
            "the fixture must exchange actual reflected heat"
        );
        assert!(
            nonuniform_source,
            "actual nodal emissivity must vary within an imported face"
        );
        close(
            n(history.get("energy").unwrap(), "stored_change_j"),
            fine_stored,
            2e-6,
        );
        if !open {
            close(fine_stored, 14.0, 1e-5);
        }
        let report = receipt.get("radiation").unwrap();
        check_exchange(report, &material, curve, open);
        let field = values(&J::parse(&solved.field).unwrap(), "temperature");
        for surface in report.get("surfaces").unwrap().as_array().unwrap() {
            for state in surface.get("material_states").unwrap().as_array().unwrap() {
                close(
                    n(state, "temperature_k"),
                    field[n(state, "vertex") as usize],
                    0.0,
                );
            }
        }
        assert!(solved.html.contains("local-temperature"));
    }
}

#[test]
fn g1_omitted_enclosure_response_keeps_original_fixed_query_behavior() {
    let material = material(&SOURCE_KNOTS, "original-enclosure-specimen");
    let mut spec = fixture(&material, true);
    let (local, local_field) = evaluate(&spec, &material.cards, None);
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .enclosure_radiation
        .as_mut()
        .unwrap()
        .emissivity_response = None;
    assert!(!print_sexpr(&spec).unwrap().contains(":emissivity-response"));
    let (fixed, fixed_field) = evaluate(&spec, &material.cards, None);
    let report = fixed.get("radiation").unwrap();
    assert_eq!(
        report.str_field("model"),
        Some("pointwise-p1-open-gray-diffuse-enclosure")
    );
    for surface in report.get("surfaces").unwrap().as_array().unwrap() {
        close(n(surface, "emissivity"), epsilon(300.0), 1e-15);
    }
    assert!(
        local_field
            .iter()
            .zip(&fixed_field)
            .any(|(a, b)| (a - b).abs() > 1e-5)
    );
    assert!(
        (n(local.get("radiation").unwrap(), "radiative_out_w") - n(report, "radiative_out_w"))
            .abs()
            > 1e-5
    );
}

fn refused_before_solution(spec: &ProjectSpec, cards: &CardPackSet) -> SolveRefusal {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let refusal = run_solve(
        &ledger,
        &CancelGate::new_clock_free(),
        &mut benign_clock(),
        &decoded,
        cards,
        &mut Vec::new(),
    )
    .unwrap_err();
    assert_eq!(
        refusal.code, "cli-solve-conduction-transient",
        "{refusal:?}"
    );
    assert_eq!(refusal.stage, Some("conduction"), "{refusal:?}");
    assert_eq!(
        stage_receipt_hashes(&ledger, refusal.run.as_ref().unwrap()).len(),
        4,
        "no failed material state may publish a conduction field, QoI or report"
    );
    refusal
}

#[test]
fn g0_local_enclosure_refuses_foreign_source_owners_and_uncovered_initial_states() {
    for open in [false, true] {
        let original = material(&SOURCE_KNOTS, "original-enclosure-specimen");
        let other = material(&SOURCE_KNOTS, "same-law-different-original-card");
        let cards = CardPackSet::admit(
            original
                .cards
                .iter()
                .chain(other.cards.iter())
                .map(|entry| {
                    raw_pack(
                        entry.kind(),
                        "independent-material-owner",
                        entry.bytes().to_vec(),
                    )
                })
                .collect(),
        )
        .unwrap();
        let mut spec = fixture(&original, open);
        let surface = &mut spec
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .enclosure_radiation
            .as_mut()
            .unwrap()
            .surfaces[1];
        surface.card = other.card;
        surface.claim = Some(other.claim);
        let refusal = refused_before_solution(&spec, &cards);
        assert!(
            refusal.what.contains("same original material card"),
            "{refusal:?}"
        );

        // The declared 300 K query passes, but initial 293.15 K has no emissivity
        // source support. Preserve the original complete material/chart range.
        let narrow = material(&[(295.0, 0.4), (450.0, 0.8)], "uncovered-initial-source");
        let refusal = refused_before_solution(&fixture(&narrow, open), &narrow.cards);
        assert!(refusal.what.contains("emissivity"), "{refusal:?}");
        assert!(
            refusal.what.contains("does not cover") || refusal.what.contains("extrapolat"),
            "{refusal:?}"
        );
    }
}

#[path = "material_emissivity/adjoint.rs"]
mod adjoint;
