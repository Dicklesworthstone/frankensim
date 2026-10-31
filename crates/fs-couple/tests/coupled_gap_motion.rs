//! Shared-gap motion against an independent 2x2 magnetic solve and RK4 dynamics.
//! Constitutive numbers are synthetic test fixtures, not measured actuator data.
#[path = "electrical_contacts/source.rs"]
#[allow(dead_code)]
mod source;

use fs_couple::material_circuit::motion::FieldMotionInput;
use fs_couple::material_circuit::motion::coupled_gap::{
    CommonMagneticGap, CoupledGapInitial, CoupledGapMotion, GapWinding,
};
use fs_couple::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, MaterialResistor};
use fs_ga::Vec3;
use fs_matdb::QueryPoint;
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::field_storage::{
    MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::Density;

const SCALE: f64 = 0.4; // Actual fixture permeability 0.8 times area 0.5.
const LEAKAGE: [f64; 2] = [0.3, 0.4];
const TURNS: [f64; 2] = [2.0, -3.0];
const MASS: f64 = 2.5;
const RESISTANCE: [f64; 2] = [0.7, 1.1];
const PORTS: [f64; 2] = [1.2, -0.9];
const FORCE: f64 = 0.4;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs().max(1.0), "{a} != {b}");
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 2e-8,
        power_tolerance_w: 2e-8,
    }
}
fn body(mass: f64) -> ResolvedBody {
    let card = source::scalar_card("density", Density::DIMS, mass);
    resolve_body(
        &[RegionBinding {
            key: "guided common pole",
            material: &card,
            point: &QueryPoint::new(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        }],
        1,
        || false,
    )
    .unwrap()
}
fn medium(gap: f64, permeability: f64) -> ResolvedInductor {
    let card = source::scalar_card(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        permeability,
    );
    ResolvedInductor::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.5,
        gap,
        1,
    )
    .unwrap()
}
fn winding(reverse: bool) -> Vec<GapWinding> {
    // Reverse input order: the owner must retain original circuit branch order.
    vec![
        GapWinding {
            branch: 1,
            signed_turns: if reverse { 3 } else { -3 },
            leakage_h: LEAKAGE[1],
        },
        GapWinding {
            branch: 0,
            signed_turns: 2,
            leakage_h: LEAKAGE[0],
        },
    ]
}
fn resistor(branch: usize, ohms: f64) -> MaterialResistor {
    let card = source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, ohms);
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn graph(reverse: bool) -> CircuitGraph {
    CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::Inductor { henries: f64::NAN }),
            if reverse {
                (0, 2, Branch::Inductor { henries: f64::NAN })
            } else {
                (2, 0, Branch::Inductor { henries: f64::NAN })
            },
            (3, 1, Branch::Resistor { ohms: f64::NAN }),
            (4, 2, Branch::Resistor { ohms: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 0 }),
            (4, 0, Branch::VoltageSource { port: 1 }),
        ],
        transformers: vec![],
    }
}
fn model(gap: f64, reverse: bool, permeability: f64, mass: f64) -> CoupledGapMotion {
    let gap =
        CommonMagneticGap::new(medium(gap, permeability), winding(reverse), [0.2, 1.2]).unwrap();
    CoupledGapMotion::new(
        graph(reverse),
        vec![resistor(2, RESISTANCE[0]), resistor(3, RESISTANCE[1])],
        gap,
        body(mass),
        CoupledGapInitial {
            fluxes_wb: vec![0.9, if reverse { 0.7 } else { -0.7 }],
            charges_c: vec![],
            electrical_ports: PORTS.to_vec(),
            velocity_m_s: 0.25,
        },
    )
    .unwrap()
}
fn input(dt: f64) -> FieldMotionInput<'static> {
    FieldMotionInput {
        electrical_ports: &PORTS,
        external_force_n: FORCE,
        duration_s: dt,
    }
}

// Independent direct inverse, rather than the production rank-one solve.
fn magnetic(gap: f64, linkage: [f64; 2]) -> ([f64; 2], f64, f64) {
    let l1 = LEAKAGE[0] + SCALE / gap * TURNS[0] * TURNS[0];
    let l2 = LEAKAGE[1] + SCALE / gap * TURNS[1] * TURNS[1];
    let m = SCALE / gap * TURNS[0] * TURNS[1];
    let determinant = l1 * l2 - m * m;
    let currents = [
        (l2 * linkage[0] - m * linkage[1]) / determinant,
        (l1 * linkage[1] - m * linkage[0]) / determinant,
    ];
    let flux = SCALE / gap * (TURNS[0] * currents[0] + TURNS[1] * currents[1]);
    let force = -flux * flux / (2.0 * SCALE);
    let energy = 0.5 * (linkage[0] * currents[0] + linkage[1] * currents[1]);
    (currents, force, energy)
}

#[test]
fn source_geometry_full_matrix_currents_and_fixed_linkage_force_agree() {
    let model = model(0.7, false, 0.8, MASS);
    let (current, force, energy) = magnetic(0.7, [0.9, -0.7]);
    for (&actual, expected) in model.winding_currents_a().unwrap().iter().zip(current) {
        close(actual, expected, 2e-14);
    }
    close(model.magnetic_energy_j().unwrap(), energy, 2e-14);
    close(model.field_force_n().unwrap(), force, 2e-14);
    close(
        model.common_flux_wb().unwrap(),
        SCALE / 0.7 * (TURNS[0] * current[0] + TURNS[1] * current[1]),
        2e-14,
    );
    for (branch, winding) in model.gap().windings().iter().enumerate() {
        assert_eq!(winding.branch, branch);
        let Branch::Inductor { henries } = model.circuit().graph().branches[branch].2 else {
            panic!()
        };
        close(
            henries,
            LEAKAGE[branch] + SCALE / 0.7 * TURNS[branch].powi(2),
            2e-14,
        );
    }
    assert_eq!(model.circuit().mutual_inductance()[0].branches, (0, 1));
    close(
        model.circuit().mutual_inductance()[0].henries,
        SCALE / 0.7 * TURNS[0] * TURNS[1],
        2e-14,
    );
    let delta = 1e-5;
    let force_difference = -(magnetic(0.7 + delta, [0.9, -0.7]).2
        - magnetic(0.7 - delta, [0.9, -0.7]).2)
        / (2.0 * delta);
    close(force, force_difference, 1e-10);
    assert_eq!(model.gap().medium().material(), medium(0.7, 0.8).material());
    assert_eq!(model.body(), &body(MASS));
}

fn derivative(x: [f64; 4]) -> [f64; 4] {
    let (current, force, _) = magnetic(x[2], [x[0], x[1]]);
    [
        PORTS[0] - RESISTANCE[0] * current[0],
        PORTS[1] - RESISTANCE[1] * current[1],
        x[3] / MASS,
        force + FORCE,
    ]
}
fn reference(duration: f64) -> [f64; 4] {
    let mut x = [0.9, -0.7, 0.7, MASS * 0.25];
    let dt = duration / 2048.0;
    for _ in 0..2048 {
        let a = derivative(x);
        let b = derivative(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| x[i] + dt * c[i]));
        x = std::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}
fn run(steps: usize) -> [f64; 4] {
    let mut model = model(0.7, false, 0.8, MASS);
    let initial = model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2);
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..steps {
        let candidate = model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                input(0.2 / steps as f64),
                controls(),
                || false,
            )
            .unwrap();
        let report = candidate.step();
        close(
            report.electrical_energy_j[0],
            model.magnetic_energy_j().unwrap(),
            2e-12,
        );
        close(
            report.electrical_energy_j[1],
            candidate.model().magnetic_energy_j().unwrap(),
            2e-12,
        );
        close(
            report.mechanical_supply_j,
            FORCE * (candidate.model().gap_m() - model.gap_m()),
            2e-12,
        );
        close(
            report.field_work_on_body_j,
            report.electrical_supply_j
                - report.heat_j
                - (report.electrical_energy_j[1] - report.electrical_energy_j[0]),
            2e-8,
        );
        assert!(report.heat_j > 0.0);
        close(report.heat_j, report.native.dissipated, 2e-10);
        close(report.energy_residual_j, 0.0, 2e-10);
        heat += report.heat_j;
        work += report.electrical_supply_j + report.mechanical_supply_j;
        let (next, _) = candidate.into_parts();
        assert_eq!(
            next.gap().medium().material(),
            model.gap().medium().material()
        );
        assert_eq!(next.body(), model.body());
        model = next;
    }
    close(
        model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2) + heat
            - work,
        initial,
        2e-9,
    );
    [
        model.circuit().fluxes_wb()[0],
        model.circuit().fluxes_wb()[1],
        model.gap_m(),
        model.momentum_kg_m_s(),
    ]
}

#[test]
fn loaded_mutual_motion_converges_to_independent_dynamics_and_closes_work_heat() {
    let reference = reference(0.2);
    let error = |actual: [f64; 4]| {
        actual
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let coarse = error(run(8));
    let fine = error(run(16));
    assert!(
        coarse > 3.5 * fine && coarse < 4.5 * fine,
        "coarse={coarse}, fine={fine}"
    );
    assert!(fine < 1e-4, "fine={fine}");
}

#[test]
fn winding_orientation_changes_linkage_sign_but_not_physical_motion_or_heat() {
    let a = model(0.7, false, 0.8, MASS);
    let b = model(0.7, true, 0.8, MASS);
    let a = a
        .step_candidate(
            a.circuit().resistors().to_vec(),
            input(0.03),
            controls(),
            || false,
        )
        .unwrap();
    let b = b
        .step_candidate(
            b.circuit().resistors().to_vec(),
            input(0.03),
            controls(),
            || false,
        )
        .unwrap();
    close(
        a.model().circuit().fluxes_wb()[0],
        b.model().circuit().fluxes_wb()[0],
        1e-10,
    );
    close(
        a.model().circuit().fluxes_wb()[1],
        -b.model().circuit().fluxes_wb()[1],
        1e-10,
    );
    close(a.model().gap_m(), b.model().gap_m(), 1e-10);
    close(a.model().velocity_m_s(), b.model().velocity_m_s(), 1e-10);
    close(a.step().heat_j, b.step().heat_j, 1e-10);
    close(
        a.step().electrical_supply_j,
        b.step().electrical_supply_j,
        1e-10,
    );
    close(
        a.model().common_flux_wb().unwrap(),
        b.model().common_flux_wb().unwrap(),
        1e-10,
    );
}

#[test]
fn actual_permeability_and_body_density_change_the_coupled_response() {
    let advance = |permeability, mass| {
        let model = model(0.7, false, permeability, mass);
        model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                input(0.03),
                controls(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0
    };
    let original = advance(0.8, MASS);
    let changed_field = advance(0.4, MASS);
    let changed_body = advance(0.8, 2.0 * MASS);
    assert!(
        (original.field_force_n().unwrap() - changed_field.field_force_n().unwrap()).abs() > 1e-3
    );
    assert!(
        (original.circuit().fluxes_wb()[0] - changed_field.circuit().fluxes_wb()[0]).abs() > 1e-5
    );
    assert!((original.velocity_m_s() - changed_body.velocity_m_s()).abs() > 1e-4);
    close(changed_body.generalized_mass_kg(), 2.0 * MASS, 1e-14);
}

#[test]
fn noncontiguous_windings_keep_other_storage_and_contact_heat() {
    let mut circuit = graph(false);
    // Insert an unrelated RLC island before the second common winding. This
    // makes original branch indices differ from compressed flux-slot indices.
    circuit.node_count = 6;
    circuit
        .branches
        .insert(1, (5, 0, Branch::Inductor { henries: 0.8 }));
    circuit
        .branches
        .push((5, 0, Branch::Capacitor { farads: 0.4 }));
    circuit
        .branches
        .push((5, 0, Branch::Resistor { ohms: f64::NAN }));
    let mut windings = winding(false);
    windings[0].branch = 2;
    let gap = CommonMagneticGap::new(medium(0.7, 0.8), windings, [0.2, 1.2]).unwrap();
    let a = source::scalar_card("density", Density::DIMS, 1.0);
    let b = source::scalar_card("density", Density::DIMS, 2.0);
    let interface = source::card(&a, &b, 0.2, 0.0);
    let contact = source::contact(19, 3, &interface, 1.0);
    let model = CoupledGapMotion::new_with_contacts_and_support(
        circuit,
        vec![resistor(3, 0.7), resistor(4, 1.1), resistor(8, 0.6)],
        vec![contact.clone()],
        gap,
        body(MASS),
        None,
        CoupledGapInitial {
            fluxes_wb: vec![0.9, 0.25, -0.7],
            charges_c: vec![0.2],
            electrical_ports: PORTS.to_vec(),
            velocity_m_s: 0.25,
        },
    )
    .unwrap();
    assert!(
        model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                input(0.01),
                controls(),
                || false
            )
            .is_err()
    );
    let candidate = model
        .step_with_contacts_candidate(
            model.circuit().resistors().to_vec(),
            vec![contact],
            input(0.01),
            controls(),
            || false,
        )
        .unwrap();
    let other_energy = |model: &CoupledGapMotion| {
        0.5 * model.circuit().fluxes_wb()[1].powi(2) / 0.8
            + 0.5 * model.circuit().charges_c()[0].powi(2) / 0.4
    };
    close(
        candidate.step().electrical_energy_j[0],
        model.magnetic_energy_j().unwrap() + other_energy(&model),
        2e-12,
    );
    close(
        candidate.step().electrical_energy_j[1],
        candidate.model().magnetic_energy_j().unwrap() + other_energy(candidate.model()),
        2e-12,
    );
    assert!(candidate.step().contacts[0].heat_j > 0.0);
    close(
        candidate.step().heat_j,
        candidate
            .step()
            .resistors
            .iter()
            .map(|h| h.heat_j)
            .sum::<f64>()
            + candidate.step().contacts[0].heat_j,
        2e-12,
    );
    assert_eq!(
        candidate.model().circuit().mutual_inductance()[0].branches,
        (0, 2)
    );
    assert_eq!(
        candidate.model().circuit().contacts(),
        model.circuit().contacts()
    );
}

#[test]
fn leakage_domain_branch_and_final_cancellation_fail_without_mutating_history() {
    for leakage in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut w = winding(false);
        w[0].leakage_h = leakage;
        assert!(CommonMagneticGap::new(medium(0.7, 0.8), w, [0.2, 1.2]).is_err());
    }
    let mut duplicate = winding(false);
    duplicate[0].branch = duplicate[1].branch;
    assert!(CommonMagneticGap::new(medium(0.7, 0.8), duplicate, [0.2, 1.2]).is_err());
    assert!(CommonMagneticGap::new(medium(0.7, 0.8), winding(false), [0.7, 1.2]).is_err());
    let mut wrong_branch = winding(false);
    wrong_branch[0].branch = 2;
    assert!(
        CoupledGapMotion::new(
            graph(false),
            vec![resistor(2, 0.7), resistor(3, 1.1)],
            CommonMagneticGap::new(medium(0.7, 0.8), wrong_branch, [0.2, 1.2]).unwrap(),
            body(MASS),
            CoupledGapInitial {
                fluxes_wb: vec![0.9, -0.7],
                charges_c: vec![],
                electrical_ports: PORTS.to_vec(),
                velocity_m_s: 0.0
            }
        )
        .is_err()
    );
    let model = model(0.7, false, 0.8, MASS);
    let saved = model.clone();
    let mut calls = 0;
    let first = model
        .step_candidate(
            model.circuit().resistors().to_vec(),
            input(0.01),
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let repeated = model
        .step_candidate(
            model.circuit().resistors().to_vec(),
            input(0.01),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(first.model(), repeated.model());
    let mut actual = 0;
    assert!(matches!(
        model.step_candidate(
            model.circuit().resistors().to_vec(),
            input(0.01),
            controls(),
            || {
                actual += 1;
                actual == calls
            }
        ),
        Err(MaterialCircuitError::Cancelled)
    ));
    let outside = FieldMotionInput {
        external_force_n: 1e8,
        ..input(0.01)
    };
    assert!(
        model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                outside,
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(model, saved);
}
