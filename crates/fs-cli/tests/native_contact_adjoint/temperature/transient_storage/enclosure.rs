//! Equal-area gray exchange shares the physical capacity/contact history.
//! These reciprocal factors are a manufactured algebraic fixture, not a claim
//! that the imported contact pair geometrically forms a closed enclosure.

use super::*;
use super::radiation::{Cooling, radiating_fixture, reservoir_heat, solve_radiating, solve_radiating_output};
use fs_project::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};

const START_K: f64 = 300.0;
const EPSILON: f64 = 0.85; // Existing immutable gray-surface fixture card.
const SIGMA: f64 = 5.670_374_419e-8;
const OPEN_BATH_K: f64 = 290.0;
const MUTUAL_VIEW: f64 = 0.25;

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap_or_else(|| panic!("missing numeric {key}: {row:?}"))
}

fn enclosed_fixture(cooling: Cooling) -> Fixture {
    let mut f = radiating_fixture(cooling, INITIAL_K);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let ambient = setup.radiation.take().unwrap();
    setup.enclosure_radiation = Some(ConductionEnclosure {
        emissivity_response: None,
        spatial_model: None,
        ambient_temperature: None,
        surfaces: ["cold", "hot"].into_iter().map(|target| EnclosureRadiatingSurface {
            name: format!("{target}-radiator"), target: target.into(),
            card: ambient.surfaces[0].card.clone(), claim: ambient.surfaces[0].claim.clone(),
            query_temperature: QtyAny::new(START_K, dims::TEMPERATURE),
        }).collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "manufactured equal-area reciprocal exchange; not geometric visibility evidence".into(),
        },
        closure_tolerance: 1e-12, reciprocity_relative_tolerance: 1e-12,
        max_iterations: 256,
        temperature_tolerance: QtyAny::new(1e-11, dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-10, dims::POWER), relaxation: 0.5,
    });
    if cooling == Cooling::Prescribed {
        for boundary in &mut setup.boundaries {
            boundary.condition = ThermalBoundaryCondition::HeatFlux {
                outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
            };
        }
    }
    f
}

fn open_fixture(cooling: Cooling) -> Fixture {
    let mut f = enclosed_fixture(cooling);
    let declaration = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap();
    declaration.ambient_temperature = Some(QtyAny::new(OPEN_BATH_K, dims::TEMPERATURE));
    declaration.view_factors = vec![vec![0.0, MUTUAL_VIEW], vec![MUTUAL_VIEW, 0.0]];
    f
}

fn mesh_masses() -> Vec<f64> {
    let mut f = fixture(true, false);
    time(&mut f).max_step.value = 1.0;
    let (receipt, _) = solve(&f, 200);
    // Only the immutable geometry/density weights are reused. Pointwise
    // enthalpy radiation is not a trajectory oracle for this area-mean law.
    let masses = values(receipt.get("transient").unwrap().get("endpoint").unwrap(),
        "reference_nodal_mass_kg");
    close(masses.iter().sum(), 20.0, 1e-10);
    masses
}

fn exchange(row: &J) -> f64 {
    let radiation = row.get("radiation").unwrap();
    let open = radiation.str_field("model") == Some("open-gray-diffuse-enclosure");
    if !open { assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure")); }
    assert_eq!(radiation.str_field("storage"), Some("volumetric-heat-capacity"));
    let patches = radiation.get("patches").unwrap().as_array().unwrap();
    assert_eq!(patches.len(), 2);
    let patch = |target| patches.iter().find(|patch| patch.str_field("region") == Some(target)).unwrap();
    let (cold, hot) = (patch("cold"), patch("hot"));
    close(n(cold, "area_m2"), 5.0, 1e-12);
    close(n(hot, "area_m2"), 5.0, 1e-12);
    let (a, b) = (n(cold, "mean_temperature_k"), n(hot, "mean_temperature_k"));
    let (cold_out, hot_out) = if open {
        // Independently solve the reflected 2x2 equations relative to the
        // black bath, then recover each original solid patch's outward heat.
        let emission = |t: f64| EPSILON * SIGMA * (t - OPEN_BATH_K) * (t + OPEN_BATH_K)
            * (t * t + OPEN_BATH_K * OPEN_BATH_K);
        let r = (1.0 - EPSILON) * MUTUAL_VIEW;
        let ja = (emission(a) + r * emission(b)) / (1.0 - r * r);
        let jb = (emission(b) + r * emission(a)) / (1.0 - r * r);
        (5.0 * (ja - MUTUAL_VIEW * jb), 5.0 * (jb - MUTUAL_VIEW * ja))
    } else {
        // Independent closed two-surface gray resistance network.
        let q = SIGMA * 5.0 * (a - b) * (a + b) * (a * a + b * b)
            / (2.0 / EPSILON - 1.0);
        (q, -q)
    };
    assert!(cold_out.abs() > 1e-5, "the fixture must exchange internal radiation");
    close(n(cold, "nonlinear_heat_w"), cold_out, 2e-8);
    close(n(hot, "nonlinear_heat_w"), hot_out, 2e-8);
    for key in ["applied_radiation_out_w", "nonlinear_radiation_out_w"] {
        close(n(radiation, key), cold_out + hot_out, 2e-8);
    }
    for key in ["prescribed_neumann_out_w", "physical_dirichlet_in_w"] {
        close(n(radiation, key), 0.0, 2e-8);
    }
    if open {
        assert!(cold_out + hot_out > 0.1, "the black bath must remove real stored heat");
        close(n(radiation, "ambient_temperature_k"), OPEN_BATH_K, 0.0);
        close(n(radiation, "ambient_absorbed_w"), cold_out + hot_out, 2e-8);
        close(n(radiation, "enclosure_energy_closure_w"), 0.0, 2e-8);
    }
    assert!(n(radiation, "physical_residual_j") <= n(radiation, "physical_residual_threshold_j"));
    assert!(n(radiation, "physical_energy_residual_j").abs() <= ENERGY_TOLERANCE_J);
    if let Some(contact) = row.get("temperature_contact").filter(|value| value.get("method").is_some()) {
        close(n(contact, "physical_residual_norm_j"), n(radiation, "physical_residual_j"), 0.0);
        close(n(contact, "physical_energy_residual_j"), n(radiation, "physical_energy_residual_j"), 0.0);
    }
    cold_out
}

fn external_heat(row: &J, cooling: Cooling, internal: f64) -> f64 {
    let radiative = n(row.get("radiation").unwrap(), "nonlinear_radiation_out_w");
    let natural = match cooling {
        Cooling::Natural => reservoir_heat(row, "cold"),
        Cooling::Mixed => reservoir_heat(row, "hot"),
        _ => 0.0,
    };
    let air = if matches!(cooling, Cooling::Air | Cooling::Mixed) {
        let coupled = row.get("conjugate").unwrap();
        let actual = coupled.get("exchange").unwrap();
        assert_eq!(actual.str_field("path"), Some("vent:cold"));
        close(n(actual, "inlet_k"), INITIAL_K, 0.0);
        close(n(actual.get("air_properties").unwrap(), "specific_heat_j_kg_k"), 1007.0, 0.0);
        let heat = 1007.0 * n(actual, "mass_flow_kg_s") * (n(actual, "outlet_k") - INITIAL_K);
        assert!(heat > 0.1);
        close(n(actual, "air_total_w"), heat, 3e-6);
        close(n(coupled, "air_heat_gain_w"), heat, 3e-6);
        close(n(coupled, "off_path_convective_out_w"), natural, 3e-6);
        close(n(coupled, "radiative_out_w"), radiative, 2e-8);
        assert!(n(coupled, "coupled_energy_residual_j").abs() <= ENERGY_TOLERANCE_J);
        heat
    } else { 0.0 };
    let outgoing = n(row.get("radiation").unwrap(), "convective_out_w");
    close(outgoing, air + natural, 3e-6);
    assert!((outgoing - (air + natural + internal)).abs() > 1e-5,
        "reflected radiation must not be added to convective or air-stream heat");
    if cooling != Cooling::Prescribed {
        let owner = row.get(if cooling == Cooling::Air { "conjugate" } else { "natural" }).unwrap();
        let contact = row.get("temperature_contact").unwrap();
        assert!(n(owner, "physical_residual_norm_j") <= n(owner, "physical_residual_tolerance_j"));
        assert!(n(owner, "coupled_energy_residual_j").abs() <= ENERGY_TOLERANCE_J);
        for key in ["physical_residual_norm_j", "physical_energy_residual_j"] {
            close(n(contact, key), n(owner, key), 0.0);
        }
        close(n(contact, "thermal_responses"), n(owner, "solid_solves"), 0.0);
        assert!(n(contact, "thermal_responses") >= n(owner, "radiation_trials"));
    }
    air + natural
}

fn audit(receipt: &J, field: &J, cooling: Cooling, card: &str, claim: &str, masses: &[f64], surface: bool) {
    let trajectory = receipt.get("transient").unwrap();
    let open = receipt.get("radiation").unwrap().str_field("model") == Some("open-gray-diffuse-enclosure");
    assert_eq!(trajectory.get("radiation").unwrap().str_field("method"), Some(if open {
        "implicit-open-gray-diffuse-enclosure"
    } else { "implicit-closed-gray-diffuse-enclosure" }));
    let contact_work = trajectory.get("temperature_contact").unwrap();
    let method = if surface { "source-backed-contact-fixed-point" } else { "source-backed-contact-newton" };
    assert_eq!(contact_work.str_field("method"), Some(method));
    let (mut fine_stored, mut min_query, mut max_query, mut responses) = (0.0, f64::INFINITY, f64::NEG_INFINITY, 0.0);
    for (grid, ends) in [("coarse", vec![0.25, 0.75, 1.0]),
        ("fine", vec![0.125, 0.25, 0.5, 0.75, 0.875, 1.0])] {
        let rows = trajectory.get(grid).unwrap().as_array().unwrap();
        assert_eq!(rows.len(), ends.len());
        let (mut start, mut supplied, mut stored, mut grid_responses) = (0.0, 0.0, 0.0, 0.0);
        for (row, end) in rows.iter().zip(ends) {
            let dt = end - start;
            close(n(row, "time_s"), end, 0.0); close(n(row, "dt_s"), dt, 0.0);
            let query = check_contact(row, claim);
            min_query = min_query.min(query); max_query = max_query.max(query);
            let contact = row.get("temperature_contact").unwrap();
            assert_eq!(contact.str_field("method"), Some(method));
            for face in contact.get("face_states").unwrap().as_array().unwrap() {
                assert_eq!(face.str_field("interface_card"), Some(card));
            }
            let internal = exchange(row);
            let external = external_heat(row, cooling, internal);
            let radiation = row.get("radiation").unwrap();
            assert!(n(radiation, "solid_solves") >= n(radiation, "radiation_trials"));
            if surface { assert!(n(radiation, "solid_solves") > n(radiation, "radiation_trials")); }
            if cooling == Cooling::Prescribed {
                close(n(contact, "thermal_responses"), n(radiation, "solid_solves"), 0.0);
            }
            let power = if start < 0.25 { 20.0 } else if start < 0.75 { 5.0 } else { 10.0 };
            let net = power - external - n(radiation, "nonlinear_radiation_out_w");
            close(n(row, "source_w"), power, 1e-10);
            close(n(row, "net_input_w"), net, 2.0 * ENERGY_TOLERANCE_J / dt + 1e-8);
            close(n(row, "stored_energy_change_j"), dt * net, 2.0 * ENERGY_TOLERANCE_J + 1e-8);
            supplied += dt * power; stored += dt * net;
            grid_responses += n(contact, "thermal_responses"); start = end;
        }
        close(supplied, 10.0, 0.0);
        if cooling == Cooling::Prescribed && !open { close(stored, 10.0, 2e-8); }
        else { assert!(stored < 0.0, "external cooling must remove real stored heat"); }
        close(n(contact_work, &format!("{grid}_thermal_responses")), grid_responses, 0.0);
        responses += grid_responses;
        if grid == "fine" { fine_stored = stored; }
    }
    close(n(contact_work, "total_thermal_responses"), responses, 0.0);
    assert!(max_query - min_query > 1e-4, "the contact source must follow the evolving field");
    let temperature = values(field, "temperature");
    assert_eq!(temperature.len(), masses.len());
    let stored: f64 = temperature.iter().zip(masses).map(|(t, m)| m * CP * (t - START_K)).sum();
    close(stored, fine_stored, 8.0 * ENERGY_TOLERANCE_J);
    let energy = trajectory.get("energy").unwrap();
    close(n(energy, "stored_change_j"), stored, 1e-8);
    close(n(energy, "integrated_net_input_j"), fine_stored, 4e-6);
    let last = trajectory.get("fine").unwrap().as_array().unwrap().last().unwrap();
    let contact_heat: f64 = last.get("temperature_contact").unwrap().get("face_states").unwrap()
        .as_array().unwrap().iter().map(|face| 0.5 * (n(face, "side_a_temperature_k")
            - n(face, "side_b_temperature_k")) / n(face, "resistance_m2_k_w")).sum();
    assert!(contact_heat.abs() > 1e-5);
    close(interface_heat(receipt), contact_heat, 1e-8);
    // The published enclosure is the accepted final physical endpoint.
    let final_radiation = receipt.get("radiation").unwrap();
    assert_eq!(final_radiation.str_field("time_integration"), Some("backward-euler-capacity"));
    let last_radiation = last.get("radiation").unwrap();
    for key in ["solid_solves", "radiation_trials"] {
        close(n(final_radiation, key), n(last_radiation, key), 0.0);
    }
    assert!(n(final_radiation, "solid_solves") >= n(final_radiation, "radiation_trials"));
    if surface { assert!(n(final_radiation, "solid_solves") > n(final_radiation, "radiation_trials")); }
    let patches = last_radiation.get("patches").unwrap().as_array().unwrap();
    for surface in final_radiation.get("surfaces").unwrap().as_array().unwrap() {
        close(n(surface, "emissivity"), EPSILON, 0.0);
        let patch = patches.iter().find(|patch| patch.str_field("region") == surface.str_field("target")).unwrap();
        close(n(surface, "mean_temperature_k"), n(patch, "mean_temperature_k"), 0.0);
        close(n(surface, "nonlinear_heat_w"), n(patch, "nonlinear_heat_w"), 0.0);
    }
    if open {
        close(n(final_radiation, "ambient_temperature_k"), OPEN_BATH_K, 0.0);
        close(n(final_radiation, "ambient_effective_area_m2"), 10.0 * (1.0 - MUTUAL_VIEW), 1e-12);
        close(n(final_radiation, "ambient_absorbed_w"), n(last_radiation, "ambient_absorbed_w"), 0.0);
        close(n(final_radiation, "nonlinear_radiative_out_w"), n(last_radiation, "nonlinear_radiation_out_w"), 0.0);
        close(n(final_radiation, "enclosure_energy_closure_w"), 0.0, 2e-8);
    }
}

#[test]
fn native_capacity_contact_closed_enclosure_redistributes_heat_and_stores_exact_workload() {
    let masses = mesh_masses();
    let mut f = enclosed_fixture(Cooling::Prescribed);
    let enclosure = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().enclosure_radiation.take();
    let (control, control_field) = solve(&f, 201);
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().enclosure_radiation = enclosure;
    let mut curve_temperature = Vec::new();
    for surface in [false, true] {
        let pack = if surface { temperature_surface_pack() } else { temperature_pack() };
        let (card, claim) = (pack.1.clone(), pack.2.clone());
        let (receipt, field) = solve_radiating(&f, 202 + usize::from(surface), pack);
        audit(&receipt, &field, Cooling::Prescribed, &card, &claim, &masses, surface);
        assert!(n(receipt.get("temperature").unwrap(), "max") < n(control.get("temperature").unwrap(), "max") - 1e-5,
            "reflected exchange must cool the hotter body, not merely appear in a receipt");
        let temperature = values(&field, "temperature");
        assert!(temperature.iter().zip(values(&control_field, "temperature")).any(|(a, b)| (a - b).abs() > 1e-5));
        if surface {
            for (a, b) in temperature.iter().zip(&curve_temperature) { close(*a, *b, 3e-6); }
        } else { curve_temperature = temperature; }
    }
}

#[test]
fn native_capacity_contact_enclosure_keeps_natural_and_air_heat_in_their_own_reservoirs() {
    let masses = mesh_masses();
    for (case, cooling) in [Cooling::Natural, Cooling::Air, Cooling::Mixed].into_iter().enumerate() {
        let f = enclosed_fixture(cooling);
        let mut curve_temperature = Vec::new();
        for surface in [false, true] {
            let pack = if surface { temperature_surface_pack() } else { temperature_pack() };
            let (card, claim) = (pack.1.clone(), pack.2.clone());
            let (receipt, field) = solve_radiating(&f, 204 + 2 * case + usize::from(surface), pack);
            audit(&receipt, &field, cooling, &card, &claim, &masses, surface);
            let temperature = values(&field, "temperature");
            if surface {
                for (a, b) in temperature.iter().zip(&curve_temperature) { close(*a, *b, 3e-6); }
            } else { curve_temperature = temperature; }
        }
    }
}

#[test]
fn native_capacity_contact_open_enclosure_preserves_bath_storage_and_all_cooling_paths() {
    let masses = mesh_masses();
    for (case, cooling) in [Cooling::Prescribed, Cooling::Natural, Cooling::Air, Cooling::Mixed].into_iter().enumerate() {
        let f = open_fixture(cooling);
        let mut curve_temperature = Vec::new();
        for surface in [false, true] {
            let pack = if surface { temperature_surface_pack() } else { temperature_pack() };
            let (card, claim) = (pack.1.clone(), pack.2.clone());
            let (receipt, field) = solve_radiating(&f, 210 + 2 * case + usize::from(surface), pack);
            audit(&receipt, &field, cooling, &card, &claim, &masses, surface);
            let temperature = values(&field, "temperature");
            if surface {
                for (a, b) in temperature.iter().zip(&curve_temperature) { close(*a, *b, 3e-6); }
            } else { curve_temperature = temperature; }
        }
    }
}

#[test]
fn native_capacity_open_enclosure_uses_the_same_owner_with_static_contact() {
    let f = open_fixture(Cooling::Prescribed);
    let masses = mesh_masses();
    let (receipt, field, _) = f.solve(0.3, None, 218);
    let trajectory = receipt.get("transient").unwrap();
    let mut integrated = 0.0;
    for row in trajectory.get("fine").unwrap().as_array().unwrap() {
        exchange(row);
        let radiation = row.get("radiation").unwrap();
        close(n(radiation, "convective_out_w"), 0.0, 0.0);
        let net = n(row, "source_w") - n(radiation, "nonlinear_radiation_out_w");
        let dt = n(row, "dt_s");
        close(n(row, "stored_energy_change_j"), dt * net, 2.0 * ENERGY_TOLERANCE_J);
        integrated += dt * net;
    }
    let stored: f64 = values(&field, "temperature").iter().zip(&masses)
        .map(|(t, mass)| mass * CP * (t - START_K)).sum();
    close(stored, integrated, 8.0 * ENERGY_TOLERANCE_J);
    assert!(stored < 0.0);
    let final_radiation = receipt.get("radiation").unwrap();
    assert_eq!(final_radiation.str_field("time_integration"), Some("backward-euler-capacity"));
    close(n(final_radiation, "ambient_absorbed_w"), n(final_radiation, "nonlinear_radiative_out_w"), 2e-8);
}

#[test]
fn native_capacity_open_contact_enclosure_refuses_unimplemented_material_controls() {
    let f = open_fixture(Cooling::Prescribed);
    for (ordinal, output) in [(219, "temperature-final-material-adjoint"), (220, "temperature-sampled-peak-material-adjoint")] {
        let (source, contact, ledger) = prepare_with_pack(&f, ordinal, Some(output), temperature_pack());
        let output = super::radiation::run_radiating(&f, &source, &contact, &ledger);
        assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
        let refusal = format!("{}{}", output.stdout, output.stderr);
        assert!(refusal.contains("radiation") || refusal.contains("enclosure"), "{refusal}");
        assert!(!refusal.contains("run_receipt"), "unsupported derivative published history: {refusal}");
    }
}

fn adjoint(receipt: &J) -> &J { receipt.get("transient").unwrap().get("adjoint").unwrap() }

fn same_history(plain: &J, plain_field: &J, differentiated: &J, field: &J) {
    assert_eq!(plain_field.get("temperature"), field.get("temperature"));
    for key in ["energy", "interfaces", "radiation", "conjugate", "natural"] {
        assert_eq!(
            forward_section(plain, key),
            forward_section(differentiated, key),
            "changed forward {key}"
        );
    }
    for key in ["coarse", "fine", "energy", "temperature_contact", "sampled_peak"] {
        assert_eq!(plain.get("transient").unwrap().get(key),
            differentiated.get("transient").unwrap().get(key), "changed trajectory {key}");
    }
    let derivative = adjoint(differentiated);
    assert_eq!(n(derivative, "reconstructed_solid_endpoints"), 0.0);
    assert_eq!(n(derivative, "reconstruction_krylov_iterations"), 0.0);
    assert_eq!(n(derivative, "revalidated_solid_endpoints"), n(derivative, "state_index"));
    assert!(n(derivative, "maximum_dual_relative_residual") < 1e-12);
}

fn bath_derivative(receipt: &J) -> f64 {
    let radiation = adjoint(receipt).get("radiation").unwrap();
    assert_eq!(radiation.str_field("method"), Some("complete-open-gray-diffuse-capacity-tangent"));
    let bath = radiation.get("black_environment").unwrap();
    close(n(bath, "temperature_k"), OPEN_BATH_K, 0.0);
    n(bath, "dtemperature_dambient_temperature_k_per_k")
}

fn agrees(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-4 * expected.abs().max(0.01),
        "complete trajectory derivative {actual:e} versus physical histories {expected:e}");
}

#[test]
fn native_capacity_enclosure_contact_adjoints_preserve_schedules_and_full_feedback() {
    for open in [false, true] {
        let mut f = if open { open_fixture(Cooling::Prescribed) } else { enclosed_fixture(Cooling::Prescribed) };
        let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries.iter_mut().find(|row| row.target == "cold").unwrap().condition =
            ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(2.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
            };
        time(&mut f).boundary_schedules = vec![TransientBoundarySchedule {
            target: "cold".into(), kind: TransientBoundaryTarget::ConvectionTemperature,
            source: "synthetic held reservoirs".into(),
            steps: [(0.25, 297.0), (0.75, 302.0), (1.0, 298.0)].into_iter()
                .map(|(until, reference)| TransientBoundaryStep {
                    until: QtyAny::new(until, dims::TIME), value: QtyAny::new(reference, dims::TEMPERATURE),
                }).collect(),
        }];
        let pack = || if open { temperature_surface_pack() } else { temperature_pack() };
        let (plain, plain_field) = solve_radiating(&f, 230, pack());
        let (receipt, field) = solve_radiating_output(&f, 231, Some("temperature-final-adjoint"), pack());
        same_history(&plain, &plain_field, &receipt, &field);
        let derivative = adjoint(&receipt);
        let vertex = n(derivative, "active_vertex") as usize;
        close(n(derivative, "value_k"), values(&field, "temperature")[vertex], 0.0);
        let radiation = derivative.get("radiation").unwrap();
        assert_eq!(n(radiation, "patch_count"), 2.0);
        assert_eq!(n(radiation, "control_steps"), n(derivative, "state_index"));
        assert_eq!(radiation.get("black_environment").is_some(), open);
        for row in radiation.get("surfaces").unwrap().as_array().unwrap() {
            assert_eq!(row.str_field("card").unwrap().len(), 64);
            assert!(n(row, "dtemperature_demissivity_k").is_finite());
        }
        let robin = &derivative.get("robin").unwrap().as_array().unwrap()[0];
        assert_eq!(robin.get("reference_history_held_fixed"), Some(&J::Bool(true)));
        close(n(robin, "dtemperature_dreference_k_per_k"), 0.0, 0.0);
        let hot = derivative.get("regional_power").unwrap().as_array().unwrap().iter()
            .find(|row| row.str_field("region") == Some("hot")).unwrap();
        let derivatives = [n(hot, "dtemperature_duniform_added_watts_k_per_w"),
            n(derivative, "dtemperature_dcapacity_multiplier_k"),
            n(derivative, "dtemperature_duniform_initial_k_per_k"),
            if open { bath_derivative(&receipt) } else { 0.0 }];
        for (coordinate, epsilon) in [0.01, 0.001, 0.01, 0.01].into_iter().enumerate() {
            if coordinate == 3 && !open { continue; }
            let mut probes = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                match coordinate {
                    0 => for step in &mut time(&mut changed).power_schedules[0].steps {
                        step.watts.value += sign * epsilon;
                    },
                    1 => for row in &mut time(&mut changed).capacities {
                        row.volumetric_heat_capacity.value *= 1.0 + sign * epsilon;
                    },
                    2 => time(&mut changed).initial_temperature = TransientInitialTemperature::Uniform(
                        QtyAny::new(START_K + sign * epsilon, dims::TEMPERATURE)),
                    _ => changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap().value += sign * epsilon,
                }
                let (_, field) = solve_radiating(&changed, 232 + 2 * coordinate + side, pack());
                probes.push(values(&field, "temperature")[vertex]);
            }
            agrees(derivatives[coordinate], (probes[1] - probes[0]) / (2.0 * epsilon));
        }
    }
}

#[test]
fn native_open_capacity_contact_peaks_reverse_only_the_selected_prefix() {
    for initial in [false, true] {
        let mut f = open_fixture(Cooling::Prescribed);
        if !initial {
            for (index, step) in time(&mut f).power_schedules[0].steps.iter_mut().enumerate() {
                step.watts.value = if index == 0 { 1000.0 } else { 0.0 };
            }
        }
        let (plain, plain_field) = solve_radiating(&f, 245, temperature_pack());
        let (receipt, field) = solve_radiating_output(&f, 246,
            Some("temperature-sampled-peak-adjoint"), temperature_pack());
        same_history(&plain, &plain_field, &receipt, &field);
        let derivative = adjoint(&receipt);
        let selected = n(derivative, "state_index") as usize;
        let intervals = derivative.get("intervals").unwrap().as_array().unwrap();
        if initial {
            assert_eq!(selected, 0);
            close(bath_derivative(&receipt), 0.0, 0.0);
            close(n(derivative, "dtemperature_dcapacity_multiplier_k"), 0.0, 0.0);
            close(n(derivative, "dtemperature_duniform_initial_k_per_k"), 1.0, 0.0);
            for surface in derivative.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap() {
                close(n(surface, "dtemperature_demissivity_k"), 0.0, 0.0);
            }
        } else {
            assert!(selected > 0 && selected < intervals.len());
            let epsilon = 0.01;
            let mut peaks = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap().value += sign * epsilon;
                let (probe, _) = solve_radiating(&changed, 247 + side, temperature_pack());
                let peak = probe.get("transient").unwrap().get("sampled_peak").unwrap();
                close(n(peak, "time_s"), n(derivative, "time_s"), 0.0);
                peaks.push(n(peak, "value_k"));
            }
            agrees(bath_derivative(&receipt), (peaks[1] - peaks[0]) / (2.0 * epsilon));
        }
        for interval in &intervals[selected..] {
            close(n(interval, "dtemperature_dsource_multiplier_k"), 0.0, 0.0);
            for region in interval.get("regional_power").unwrap().as_array().unwrap() {
                close(n(region, "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
            }
        }
    }
}

#[test]
fn native_open_capacity_static_contact_keeps_existing_natural_air_derivatives() {
    let f = open_fixture(Cooling::Mixed);
    let (plain, plain_field, _) = f.solve(0.3, None, 250);
    let (receipt, field, _) = f.solve(0.3, Some("temperature-final-adjoint"), 251);
    same_history(&plain, &plain_field, &receipt, &field);
    let vertex = n(adjoint(&receipt), "active_vertex") as usize;
    let epsilon = 0.01;
    let mut probes = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
        changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap().value += sign * epsilon;
        let (_, field, _) = changed.solve(0.3, None, 252 + side);
        probes.push(values(&field, "temperature")[vertex]);
    }
    agrees(bath_derivative(&receipt), (probes[1] - probes[0]) / (2.0 * epsilon));
}
