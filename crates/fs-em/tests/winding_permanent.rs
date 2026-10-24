//! Fixed-remanence affine fields and native circuits. Source cards here are
//! manufactured analytic models, not measurements of a named magnet or wire.
mod common;
use common::{config, ends, layered, limits, near};
use fs_em::{
    force::MaxwellForceLimits,
    planar::{
        MagneticError, PlanarMagnetostatics,
        permanent::{MagneticEnergyReference, PlanarMagneticMaterial, PlanarRecoilMaterial},
    },
    transfer::{CellSample, PlanarFrame},
    winding::{
        PlanarWindingBasis, WindingConfig, WindingError,
        circuit::WindingCircuit,
        ohmic::{
            OhmicWindingCircuit, OhmicWindingControls, OhmicWindingError, WindingConductor,
            WindingMaterialRebindControls,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph};

const WIDTH: f64 = 0.01;
const HEIGHT: f64 = 0.02;
const DEPTH: f64 = 0.03;
const PROFILE: f64 = 3e4;
const BOUNDARY: [usize; 4] = [0, 3, 2, 5];

fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn material_card(
    br: Option<f64>,
    transverse: f64,
    mu: f64,
) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
    let mut data = vec![(
        ElectricalPropertyBasis::Resistivity.property_name(),
        ElectricalPropertyBasis::Resistivity.dimensions(),
        2e-4,
        ScalarAdmissibility::StrictlyPositive,
    )];
    if let Some(br) = br {
        data.extend([
            (
                RECOIL_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                mu,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                REMANENT_FLUX_DENSITY_PROPERTY,
                FLUX_DENSITY_DIMS,
                br,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                -1e6,
                ScalarAdmissibility::Finite,
            ),
            (
                RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                1e6,
                ScalarAdmissibility::Finite,
            ),
            (
                RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                transverse,
                ScalarAdmissibility::NonNegative,
            ),
        ]);
    } else {
        data.push((
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            mu,
            ScalarAdmissibility::StrictlyPositive,
        ));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value, admissibility) in data {
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, admissibility).unwrap());
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured affine winding material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic affine winding".into(),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        requirements,
    )
}
fn state(
    card: &MaterialCard,
    requirements: &[ScalarPropertyRequirement],
) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(300.0),
        requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn problem(br: Option<f64>, transverse: f64) -> PlanarMagnetostatics {
    let (left, requirements) = material_card(br, transverse, 2e-6);
    let left = state(&left, &requirements);
    let left = if br.is_some() {
        PlanarMagneticMaterial::Recoil(
            PlanarRecoilMaterial::new(
                ResolvedRecoilMaterial::try_from_material(left).unwrap(),
                [0.0, 1.0],
            )
            .unwrap(),
        )
    } else {
        PlanarMagneticMaterial::Linear(left)
    };
    let (right, requirements) = material_card(None, transverse, 1e-6);
    let right = PlanarMagneticMaterial::Linear(state(&right, &requirements));
    PlanarMagnetostatics::new_with_material_models(
        layered(WIDTH, HEIGHT, DEPTH),
        vec![left.clone(), left, right.clone(), right],
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn winding_config() -> WindingConfig {
    WindingConfig {
        max_windings: 4,
        max_basis_values: 1000,
        pairing_tolerance_h: 1e-10,
        current_balance_tolerance_turns: 1e-9,
        inverse_residual_tolerance: 1e-10,
    }
}
fn profile(p: f64) -> Vec<f64> {
    vec![0.0, 0.0, p, -p]
}
fn basis<'a>(
    problem: &'a PlanarMagnetostatics,
    profiles: Vec<Vec<f64>>,
    ends: Vec<f64>,
) -> PlanarWindingBasis<'a> {
    PlanarWindingBasis::new_with_fixed_remanence(
        problem,
        profiles,
        &BOUNDARY,
        ends,
        winding_config(),
        1e-9,
        config(),
        &mut || false,
    )
    .unwrap()
}
fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: 1.0 }),
            (2, 0, Branch::Resistor { ohms: 2.0 }),
        ],
        transformers: vec![],
    }
}
fn conductors(circuit: &WindingCircuit<'_>, triangles: &[usize]) -> Vec<WindingConductor> {
    triangles
        .iter()
        .map(|&triangle| WindingConductor {
            triangle,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                circuit.basis().problem().materials()[triangle].clone(),
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap(),
        })
        .collect()
}
fn heat_controls() -> OhmicWindingControls {
    OhmicWindingControls {
        power_tolerance_w: 1e-9,
        energy_tolerance_j: 1e-9,
    }
}
fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-8 * expected.abs().max(1e-12),
        "{actual} != {expected}"
    );
}
fn active_inductance(p: f64) -> f64 {
    // Two free interface vertices. Symmetric/antisymmetric Ritz directions
    // have stiffness a and a+2b; the balanced right-side profile loads only
    // the upper vertex by -p*V/3. This is an independent closed-form P1 oracle.
    let volume = WIDTH * HEIGHT * DEPTH / 2.0;
    let a = volume * (1.0 / 2e-6 + 1.0 / 1e-6) / (WIDTH * WIDTH);
    let b = volume * (1.0 / 2e-6 + 1.0 / 1e-6) / (HEIGHT * HEIGHT);
    (p * volume / 3.0).powi(2) * (a + b) / (a * (a + 2.0 * b))
}

#[test]
fn affine_linkage_and_shifted_storage_match_layered_magnet_oracles_and_direct_fields() {
    let problem = problem(Some(0.8), 1e6);
    let winding = basis(&problem, vec![profile(PROFILE)], vec![0.04]);
    let magnet_b = 0.8 * 1e-6 / (2e-6 + 1e-6);
    let expected_bias = PROFILE * HEIGHT * DEPTH * WIDTH * WIDTH * magnet_b / 6.0;
    relative(winding.bias_flux_linkages_wb()[0], expected_bias);
    relative(winding.active_inductance_h()[0], active_inductance(PROFILE));
    relative(winding.inductance_h()[0], 0.04 + active_inductance(PROFILE));
    let h = -0.8 / (2e-6 + 1e-6);
    relative(
        winding.bias_energy_j(),
        0.5 * WIDTH * HEIGHT * DEPTH * (2e-6 + 1e-6) * h * h,
    );
    assert!(winding.maximum_bias_stationarity_defect_wb() <= 1e-9);
    assert!(winding.unit_fields()[0].is_incremental_response());
    assert_eq!(
        winding.unit_fields()[0].energy_reference(),
        MagneticEnergyReference::ZeroFlux
    );
    assert!(
        winding.unit_fields()[0]
            .recoil_operating_points()
            .iter()
            .all(Option::is_none)
    );
    for current in [-0.4, 0.0, 0.7] {
        let total = winding.flux_linkages_wb(&[current]).unwrap();
        relative(
            total[0],
            expected_bias + winding.inductance_h()[0] * current,
        );
        near(winding.currents_a(&total).unwrap()[0], current);
        near(
            winding
                .currents_from_incremental_flux_wb(
                    &winding.incremental_flux_linkages_wb(&[current]).unwrap(),
                )
                .unwrap()[0],
            current,
        );
        let field = winding
            .field_for_currents(&[current], &mut || false)
            .unwrap();
        assert!(!field.is_incremental_response());
        let direct = problem
            .solve(
                &profile(PROFILE * current),
                &ends(0.0),
                config(),
                &mut || false,
            )
            .unwrap();
        for (got, expected) in field
            .magnetic_flux_density_t()
            .iter()
            .zip(direct.magnetic_flux_density_t())
        {
            near(got[0], expected[0]);
            near(got[1], expected[1]);
        }
        near(
            field.energy_j(),
            winding.bias_energy_j() + 0.5 * winding.active_inductance_h()[0] * current * current,
        );
        near(
            winding
                .pair_flux_linkages_wb(field.potential_wb_m(), &mut || false)
                .unwrap()[0]
                + 0.04 * current,
            total[0],
        );
    }
}

#[test]
fn generic_entry_accepts_linear_materials_and_data_only_magnet_swaps_do_not_create_free_current() {
    let linear = problem(None, 1e6);
    let plain = basis(&linear, vec![profile(PROFILE)], vec![0.04]);
    assert!(plain.bias_field().is_none());
    assert_eq!(plain.bias_flux_linkages_wb(), &[0.0]);
    assert_eq!(plain.bias_energy_j(), 0.0);
    let mut previous = None;
    for br in [0.8, 1.2] {
        let problem = problem(Some(br), 1e6);
        let winding = basis(&problem, vec![profile(PROFILE)], vec![0.04]);
        let values = (
            winding.inductance_h()[0],
            winding.bias_flux_linkages_wb()[0],
            winding.bias_energy_j(),
        );
        if let Some((l, flux, energy)) = previous {
            assert_eq!(values.0, l);
            relative(values.1, 1.5 * flux);
            relative(values.2, 2.25 * energy);
        }
        previous = Some(values);
        let circuit = winding
            .into_circuit(graph(), vec![1], &mut || false)
            .unwrap();
        let initial = circuit
            .consistent_initial_state(&[0.0], &[], &[0.0], &mut || false)
            .unwrap();
        assert_eq!(initial[circuit.flux_coordinates()[0]], 0.0);
        relative(
            circuit.total_flux_linkages_wb(&initial).unwrap()[0],
            values.1,
        );
        assert_eq!(circuit.dae().system.hamiltonian(&initial), 0.0);
        relative(
            circuit.stored_energy_with_bias_j(&initial).unwrap(),
            values.2,
        );
        let (step, defect) = circuit
            .step_audited(&initial, &[0.0], 0.001, &mut || false)
            .unwrap();
        near(circuit.currents_a(&step.x).unwrap()[0], 0.0);
        near(step.supplied, 0.0);
        near(step.dissipated, 0.0);
        near(defect, 0.0);
    }
}

#[test]
fn actual_affine_circuit_and_ohmic_heat_follow_the_native_rl_step() {
    let problem = problem(Some(0.8), 1e6);
    let circuit = basis(&problem, vec![profile(PROFILE)], vec![0.04])
        .into_circuit(graph(), vec![1], &mut || false)
        .unwrap();
    let initial = circuit
        .consistent_initial_state(&[0.2], &[], &[1.4], &mut || false)
        .unwrap();
    let l = circuit.basis().inductance_h()[0];
    let dt = 0.001;
    let (step, defect) = circuit
        .step_audited(&initial, &[1.4], dt, &mut || false)
        .unwrap();
    let expected = ((l - dt) * 0.2 + dt * 1.4) / (l + dt);
    near(circuit.currents_a(&step.x).unwrap()[0], expected);
    near(step.delta_h, 0.5 * l * (expected * expected - 0.2 * 0.2));
    near(
        step.dissipated,
        dt * 2.0 * f64::midpoint(0.2, expected).powi(2),
    );
    assert!(defect < 1e-9);
    let ohmic = OhmicWindingCircuit::new(
        &circuit,
        conductors(&circuit, &[2, 3]),
        0.0,
        &[0.2],
        &[],
        &[1.4],
        &mut || false,
    )
    .unwrap();
    relative(
        ohmic.resistance_ohm()[0],
        2e-4 * PROFILE * PROFILE * WIDTH * HEIGHT * DEPTH,
    );
    let candidate = ohmic
        .step_to_candidate(
            dt,
            ohmic.conductors().to_vec(),
            &[1.4],
            heat_controls(),
            &mut || false,
        )
        .unwrap();
    let r = 2.0 + ohmic.resistance_ohm()[0];
    let expected = ((l - 0.5 * dt * r) * 0.2 + dt * 1.4) / (l + 0.5 * dt * r);
    near(candidate.circuit().currents_a().unwrap()[0], expected);
    near(
        candidate.step().heat_j,
        dt * r * f64::midpoint(0.2, expected).powi(2),
    );
    assert!(candidate.step().energy_residual_j.abs() < 1e-9);
    near(
        candidate.step().native.supplied,
        candidate.step().native.delta_h + candidate.step().heat_j,
    );
    assert_eq!(ohmic.time_s(), 0.0);
    near(
        circuit
            .total_flux_linkages_wb(candidate.circuit().state())
            .unwrap()[0],
        circuit.basis().bias_flux_linkages_wb()[0] + l * expected,
    );
}

#[test]
fn reversed_turns_and_branch_mapping_preserve_affine_mutual_currents() {
    let problem = problem(Some(0.8), 1e6);
    let p = 1000.0 * PROFILE;
    let winding = basis(
        &problem,
        vec![profile(p), profile(-2.0 * p)],
        vec![0.04, 0.06],
    );
    relative(
        winding.bias_flux_linkages_wb()[1],
        -2.0 * winding.bias_flux_linkages_wb()[0],
    );
    let active = active_inductance(p);
    for (&got, want) in winding.active_inductance_h().iter().zip([
        active,
        -2.0 * active,
        -2.0 * active,
        4.0 * active,
    ]) {
        relative(got, want);
    }
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::Resistor { ohms: 2.0 }),
            (2, 0, Branch::Inductor { henries: 1.0 }),
            (2, 0, Branch::Resistor { ohms: 3.0 }),
        ],
        transformers: vec![],
    };
    let circuit = winding
        .into_circuit(graph, vec![2, 0], &mut || false)
        .unwrap();
    assert_eq!(circuit.flux_coordinates(), &[1, 0]);
    let currents = [0.3, -0.2];
    let initial = circuit
        .consistent_initial_state(&currents, &[], &[], &mut || false)
        .unwrap();
    let l = circuit.basis().inductance_h();
    let dt = 0.002;
    let a = l[0] + 1.5 * dt;
    let b = l[1];
    let c = l[2];
    let d = l[3] + dt;
    let rhs = [
        (l[0] - 1.5 * dt) * currents[0] + l[1] * currents[1],
        l[2] * currents[0] + (l[3] - dt) * currents[1],
    ];
    let expected = [
        (d * rhs[0] - b * rhs[1]) / (a * d - b * c),
        (a * rhs[1] - c * rhs[0]) / (a * d - b * c),
    ];
    let (step, defect) = circuit
        .step_audited(&initial, &[], dt, &mut || false)
        .unwrap();
    for (got, want) in circuit.currents_a(&step.x).unwrap().iter().zip(expected) {
        near(*got, want);
    }
    assert!(defect < 1e-9);
}

#[test]
fn only_actual_fields_are_source_admitted_and_rejected_candidates_preserve_history() {
    // One ampere would exceed this transverse limit, although the response
    // derivative is well defined. Valid small currents must remain usable.
    let problem = problem(Some(0.8), 0.5);
    assert!(matches!(
        problem.solve(&profile(PROFILE), &ends(0.0), config(), &mut || false),
        Err(MagneticError::Material(_))
    ));
    let winding = basis(&problem, vec![profile(PROFILE)], vec![0.04]);
    assert!(winding.admit_currents(&[0.01], &mut || false).is_ok());
    assert!(winding.admit_currents(&[0.1], &mut || false).is_err());
    let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    let response = &winding.unit_fields()[0];
    assert!(response.is_incremental_response());
    assert!(response.impressed_current_density_a_m2().is_some());
    assert!(matches!(
        response.maxwell_wrench(
            &[0],
            &frame,
            [0.0; 3],
            MaxwellForceLimits {
                max_cells: 1000,
                max_edges: 1000,
            },
            &mut || false,
        ),
        Err(MagneticError::Invalid(
            "Maxwell force requires the actual operating field including permanent-magnet bias"
        ))
    ));
    assert!(
        winding
            .sample_for_currents(
                &[0.1],
                &frame,
                &[CellSample {
                    triangle: 2,
                    position_m: [5.0 * WIDTH / 3.0, HEIGHT / 3.0, DEPTH / 2.0]
                }],
                1e-12,
                1,
                &mut || false
            )
            .is_err()
    );
    let circuit = winding
        .into_circuit(graph(), vec![1], &mut || false)
        .unwrap();
    assert!(
        circuit
            .consistent_initial_state(&[0.1], &[], &[0.0], &mut || false)
            .is_err()
    );
    let initial = circuit
        .consistent_initial_state(&[0.01], &[], &[0.0], &mut || false)
        .unwrap();
    let original = initial.clone();
    assert!(matches!(
        circuit.step_audited(&initial, &[10.0], 0.001, &mut || false),
        Err(WindingError::Magnetic(MagneticError::Material(_)))
    ));
    assert_eq!(initial, original);
    let ohmic = OhmicWindingCircuit::new(
        &circuit,
        conductors(&circuit, &[2, 3]),
        0.0,
        &[0.01],
        &[],
        &[0.0],
        &mut || false,
    )
    .unwrap();
    let original = ohmic.state().to_vec();
    assert!(matches!(
        ohmic.step_to_candidate(
            0.001,
            ohmic.conductors().to_vec(),
            &[10.0],
            heat_controls(),
            &mut || false
        ),
        Err(OhmicWindingError::Winding(WindingError::Magnetic(
            MagneticError::Material(_)
        )))
    ));
    assert_eq!(ohmic.state(), original);
    assert_eq!(ohmic.time_s(), 0.0);
    assert!(
        circuit
            .step_audited(&initial, &[0.0], 0.001, &mut || true)
            .is_err()
    );
}

#[test]
fn a_conducting_magnet_uses_the_original_complete_recoil_source_without_generic_mu() {
    let problem = problem(Some(0.8), 1e6);
    assert!(
        problem.materials()[0]
            .property(MAGNETIC_PERMEABILITY_PROPERTY)
            .is_none()
    );
    let winding = basis(
        &problem,
        vec![vec![PROFILE, -PROFILE, 0.0, 0.0]],
        vec![0.04],
    );
    let circuit = winding
        .into_circuit(graph(), vec![1], &mut || false)
        .unwrap();
    let ohmic = OhmicWindingCircuit::new(
        &circuit,
        conductors(&circuit, &[0, 1]),
        0.0,
        &[0.01],
        &[],
        &[0.0],
        &mut || false,
    )
    .unwrap();
    assert!(
        ohmic
            .step_to_candidate(
                0.001,
                ohmic.conductors().to_vec(),
                &[0.0],
                heat_controls(),
                &mut || false
            )
            .is_ok()
    );
    let (card, _) = material_card(Some(0.8), 1e6, 2e-6);
    let mut changed = ohmic.conductors().to_vec();
    for source in &mut changed {
        let hot = source
            .conductor
            .material()
            .at_query_point(&card, &point(400.0))
            .unwrap();
        source.conductor = ResolvedConductorState::try_from_material_with_basis(
            hot,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap();
    }
    assert!(matches!(
        ohmic.step_to_candidate(0.001, changed, &[0.0], heat_controls(), &mut || false),
        Err(OhmicWindingError::Invalid(
            "conducting recoil material requires its exact original fixed-state source bundle"
        ))
    ));
    let (foreign, requirements) = material_card(Some(0.9), 1e6, 2e-6);
    let mut changed = ohmic.conductors().to_vec();
    changed[0].conductor = ResolvedConductorState::try_from_material_with_basis(
        state(&foreign, &requirements),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    assert!(
        ohmic
            .step_to_candidate(0.001, changed, &[0.0], heat_controls(), &mut || false)
            .is_err()
    );
    assert_eq!(ohmic.time_s(), 0.0);
}

#[test]
fn affine_bias_storage_is_charged_to_the_existing_capacity_and_has_explicit_admission() {
    let problem = problem(Some(0.8), 1e6);
    let mut budget = winding_config();
    let vertices = problem.mesh().vertices().len();
    let cells = problem.mesh().faces().len();
    let response_values = 2 * vertices + 7 * cells;
    let retained_values = response_values + 2 * vertices + 12 * cells + 1;
    for capacity in [response_values, retained_values - 1, retained_values] {
        budget.max_basis_values = capacity;
        assert_eq!(
            PlanarWindingBasis::new_with_fixed_remanence(
                &problem,
                vec![profile(PROFILE)],
                &BOUNDARY,
                vec![0.04],
                budget,
                1e-9,
                config(),
                &mut || false
            )
            .is_ok(),
            capacity == retained_values
        );
    }
    for allowance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PlanarWindingBasis::new_with_fixed_remanence(
                &problem,
                vec![profile(PROFILE)],
                &BOUNDARY,
                vec![0.04],
                winding_config(),
                allowance,
                config(),
                &mut || false
            )
            .is_err()
        );
    }
}

#[test]
fn material_replacement_requires_unbiased_old_and_new_flux_coordinates() {
    for (old_br, new_br) in [
        (None, None),
        (Some(0.8), None),
        (None, Some(0.8)),
        (Some(0.8), Some(0.9)),
    ] {
        let old_problem = problem(old_br, 1e6);
        let new_problem = problem(new_br, 1e6);
        let old = basis(&old_problem, vec![profile(PROFILE)], vec![0.04])
            .into_circuit(graph(), vec![1], &mut || false)
            .unwrap();
        let new = basis(&new_problem, vec![profile(PROFILE)], vec![0.04])
            .into_circuit(graph(), vec![1], &mut || false)
            .unwrap();
        let ohmic = OhmicWindingCircuit::new(
            &old,
            conductors(&old, &[2, 3]),
            0.0,
            &[0.01],
            &[],
            &[0.0],
            &mut || false,
        )
        .unwrap();
        let original = ohmic.state().to_vec();
        let candidate = ohmic.rebind_material_candidate(
            &new,
            conductors(&new, &[2, 3]),
            &[0.0],
            WindingMaterialRebindControls {
                maximum_absolute_work_j: 0.0,
                energy_tolerance_j: 1e-9,
            },
            &mut || false,
        );
        if old_br.is_some() || new_br.is_some() {
            assert!(matches!(
                candidate,
                Err(OhmicWindingError::Invalid(
                    "fixed-remanence material replacement requires total-flux offset and bias-storage work ownership"
                ))
            ));
        } else {
            let (replacement, report) = candidate.unwrap();
            assert_eq!(report.supplied_work_j, 0.0);
            near(replacement.currents_a().unwrap()[0], 0.01);
            assert_eq!(replacement.time_s(), 0.0);
        }
        assert_eq!(ohmic.state(), original);
        assert_eq!(ohmic.time_s(), 0.0);
    }
}
