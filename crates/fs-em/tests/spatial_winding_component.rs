//! fs-em's 3D saturating winding circuit as one Component: its ports, a step
//! that is the owner's own step at standstill bit for bit, heat that leaves
//! through the ports as the owner books it, the stored energy the state
//! carries against the owner's sums, the kernel's audit, conformance,
//! refusals and identity; then the same with conductor and resistor laws at
//! the temperatures the heat ports receive. The magnetic lane is in
//! `magnetic.rs`, and, with `magnetic-hysteresis`, the play-memory owner in
//! `hysteretic.rs`.
#![cfg(feature = "component")]

#[cfg(feature = "magnetic-hysteresis")]
#[path = "spatial_winding_component/hysteretic.rs"]
mod hysteretic;
#[path = "spatial_winding_component/magnetic.rs"]
mod magnetic;
#[path = "nonlinear_eddy/common.rs"]
mod material;

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::spatial::nonlinear::SpatialNonlinearMaterial;
use fs_em::spatial::vector_potential::induction::nonlinear::SpatialNonlinearInduction;
use fs_em::spatial::vector_potential::induction::nonlinear::SpatialNonlinearInductionControls;
use fs_em::spatial::vector_potential::induction::nonlinear::winding::{
    SpatialNonlinearWindingCandidate, SpatialNonlinearWindingCircuit,
    SpatialNonlinearWindingConfig, SpatialNonlinearWindingStep, SpatialNonlinearWindings,
    SpatialWindingConductor,
};
use fs_em::spatial::vector_potential::induction::{
    SpatialInductionControls, SpatialInductionLimits,
};
use fs_em::spatial::vector_potential::nonlinear::{
    NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics,
};
use fs_em::spatial::vector_potential::{
    VectorMagneticBoundary, VectorMagneticControls, VectorPotentialFace,
};
use fs_em::winding::component::ResistorLaw;
use fs_em::winding::component::spatial::{
    SpatialConductorLaw, SpatialHeatReceiver, SpatialLawElement, SpatialWinding,
    SpatialWindingError, SpatialWindingState,
};
use fs_em::winding::nonlinear::{
    NonlinearWindingControls, NonlinearWindingInitial, resistors::GraphResistorMaterial,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductor, ResolvedConductorState};
use fs_material::field_storage::{
    nonlinear::REVERSIBLE_FIELD_STRENGTH_PROPERTY, permanent::FIELD_STRENGTH_DIMS,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Carrier, Class, Component, EndOfStep, PortLayout, PortRole, Proposal, audit_step,
};
use fs_qty::{Density, Dims};
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};

const DT: f64 = 0.1;
const VOLTS: f64 = 2.0;
const GRAPH_OHM: f64 = 0.3;
const FILL: f64 = 0.6;
const TOLERANCE_J: f64 = 1e-9;
/// The cube's assembly bounds.
const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 1000,
    max_dofs: 4000,
    max_triplets: 100_000,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The card: a saturating reversible law, and a resistivity constant up to
/// 300 K that doubles by 500 K, where the card's validity ends.
fn card(sigma: f64) -> MaterialCard {
    material::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], sigma)
}

/// The unit cube in 48 tetrahedra of 1/48 m³, each with the card's
/// reversible law.
fn cube(card: &MaterialCard) -> SpatialNonlinearVectorMagnetostatics {
    cube_of(SpatialNonlinearMaterial::Reversible(material::law(card)))
}

/// The unit cube in 48 tetrahedra of 1/48 m³, each with `material`.
fn cube_of(material: SpatialNonlinearMaterial) -> SpatialNonlinearVectorMagnetostatics {
    let mut positions = vec![];
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let id = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = vec![];
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i, j + 1, k),
                    id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1),
                    id(i + 1, j + 1, k + 1),
                );
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    let count = tets.len();
    SpatialNonlinearVectorMagnetostatics::new(
        positions,
        tets,
        vec![material; count],
        LIMITS,
        &mut || false,
    )
    .unwrap()
}

/// The winding's unit-current profile, `2 curl(N_center e_z)`: a closed
/// loop around the centre vertex, carried by the 24 cells that touch it.
fn profile(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<[f64; 3]> {
    let center = p
        .positions_m()
        .iter()
        .position(|x| x.map(f64::to_bits) == [0.5_f64.to_bits(); 3])
        .unwrap() as u32;
    p.mesh()
        .tets
        .iter()
        .map(|tet| {
            if !tet.contains(&center) {
                return [0.0; 3];
            }
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != center).collect();
            let point = p.positions_m()[face[0] as usize];
            let normal = cross(
                sub(p.positions_m()[face[1] as usize], point),
                sub(p.positions_m()[face[2] as usize], point),
            );
            let height = dot(normal, sub(p.positions_m()[center as usize], point));
            cross(normal.map(|v| v / height), [0.0, 0.0, 2.0])
        })
        .collect()
}

/// The cells the winding crosses (`true`) or does not.
fn cells(p: &SpatialNonlinearVectorMagnetostatics, crossed: bool) -> Vec<usize> {
    profile(p)
        .iter()
        .enumerate()
        .filter(|(_, j)| j.iter().any(|&v| v != 0.0) == crossed)
        .map(|(cell, _)| cell)
        .collect()
}

fn controls() -> NonlinearWindingControls<SpatialNonlinearInductionControls> {
    NonlinearWindingControls {
        field: SpatialNonlinearInductionControls {
            induction: SpatialInductionControls {
                relative_tolerance: 2e-11,
                absolute_residual_a: 1e-11,
                restart: 64,
                maximum_cycles: 100,
                maximum_krylov_scalars: 2_000_000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 2e-9,
                energy_tolerance_j: 2e-8,
                maximum_cell_magnetic_peclet: 100.0,
            },
            maximum_newton_iterations: 50,
            maximum_backtracks: 35,
        },
        maximum_field_trials: 200,
        maximum_newton_iterations: 25,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 2e-7,
        power_tolerance_w: 2e-7,
        energy_tolerance_j: 3e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}

/// The field at rest with these core conductors.
fn history(
    p: &SpatialNonlinearVectorMagnetostatics,
    core: Vec<Option<ResolvedConductorState>>,
) -> SpatialNonlinearInduction<'_> {
    let cg = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    let boundary = VectorMagneticBoundary {
        potential: p
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    };
    let rest = p
        .solve(
            &vec![[0.0; 3]; p.mesh().tets.len()],
            &boundary,
            &vec![0.0; p.mesh().edges.len()],
            NonlinearVectorControls {
                field: VectorMagneticControls {
                    field: cg,
                    gauge: cg,
                    maximum_gauge_product_terms: 2_000_000,
                    current_tolerance_a: 2e-8,
                    flux_tolerance_wb: 2e-9,
                    energy_tolerance_j: 2e-8,
                },
                relative_residual_tolerance: 2e-11,
                maximum_newton_iterations: 50,
                maximum_backtracks: 35,
            },
            &mut || false,
        )
        .unwrap();
    SpatialNonlinearInduction::new(
        &rest,
        core,
        SpatialInductionLimits {
            maximum_unknowns: 20000,
            maximum_conducting_cells: 1000,
            maximum_assembly_terms: 1_000_000,
        },
        2e-9,
        &mut || false,
    )
    .unwrap()
}

/// The card's conductor at `kelvin`, refused rather than unwrapped outside
/// the card's validity.
fn resolve(card: &MaterialCard, kelvin: f64) -> Result<ResolvedConductorState, String> {
    let basis = ElectricalPropertyBasis::Resistivity;
    let scalar = |name: &str, dims: Dims, admissibility| {
        ScalarPropertyRequirement::try_new(name, dims, admissibility).map_err(|e| format!("{e:?}"))
    };
    let requirements = [
        basis.requirement().map_err(|e| format!("{e:?}"))?,
        scalar(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
        scalar(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        scalar(
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        scalar(
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let bundle = resolve_material_state_point(
        card,
        &material::point(kelvin),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| format!("{e:?}"))?;
    ResolvedConductorState::try_from_material_with_basis(bundle, basis)
        .map_err(|e| format!("{e:?}"))
}

/// Core conductors at `kelvin` in the 24 cells the winding does not cross.
fn core_at(
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    kelvin: f64,
) -> Vec<Option<ResolvedConductorState>> {
    let conductor = resolve(card, kelvin).unwrap();
    let mut core = vec![None; p.mesh().tets.len()];
    for cell in cells(p, false) {
        core[cell] = Some(conductor.clone());
    }
    core
}

/// The winding's strands at `kelvin`, filling `fill` of each cell it
/// crosses.
fn strands_at(
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    kelvin: f64,
    fill: f64,
) -> Vec<SpatialWindingConductor> {
    let conductor = resolve(card, kelvin).unwrap();
    cells(p, true)
        .into_iter()
        .map(|cell| SpatialWindingConductor {
            winding: 0,
            cell,
            fill_fraction: fill,
            conductor: conductor.clone(),
        })
        .collect()
}

/// A voltage source driving the winding into a graph resistor (branch 2)
/// with `across` beside it. The winding has an end inductance and the cells
/// it does not cross carry core conductors at `core_k`, so every storage
/// term and every kind of heat is present.
fn series_on<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    core_k: f64,
    across: Branch,
    fill: f64,
) -> SpatialNonlinearWindingCircuit<'m> {
    series_of(
        p,
        core_at(p, card, core_k),
        strands_at(p, card, 300.0, fill),
        across,
    )
}

/// [`series_on`] with these core conductors and strands.
fn series_of(
    p: &SpatialNonlinearVectorMagnetostatics,
    core: Vec<Option<ResolvedConductorState>>,
    strands: Vec<SpatialWindingConductor>,
    across: Branch,
) -> SpatialNonlinearWindingCircuit<'_> {
    SpatialNonlinearWindingCircuit::new(
        history(p, core),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: f64::NAN }),
                (2, 0, Branch::Resistor { ohms: GRAPH_OHM }),
                (2, 0, across),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile(p)],
            end_inductance_h: vec![0.03],
        },
        strands,
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: if matches!(across, Branch::Capacitor { .. }) {
                vec![0.0]
            } else {
                vec![]
            },
            external_inputs: vec![0.0],
        },
        SpatialNonlinearWindingConfig {
            max_windings: 4,
            max_profile_values: 100_000,
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}

/// [`series_on`] with a capacitor beside the graph resistor.
fn series<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialNonlinearWindingCircuit<'m> {
    series_on(p, card, 300.0, Branch::Capacitor { farads: 0.05 }, FILL)
}

/// The strands, the core and the graph resistor.
fn receivers(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<SpatialHeatReceiver> {
    vec![
        SpatialHeatReceiver {
            cells: cells(p, true),
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: cells(p, false),
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: vec![],
            resistors: vec![2],
        },
    ]
}

fn component<'m>(
    circuit: &SpatialNonlinearWindingCircuit<'m>,
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialWinding<'m> {
    SpatialWinding::new(
        circuit,
        core_at(p, card, 300.0),
        strands_at(p, card, 300.0, FILL),
        receivers(p),
        controls(),
        TOLERANCE_J,
    )
    .unwrap()
}

/// The owner's own step from `owner` at standstill.
fn owner_step<'m>(
    owner: &SpatialNonlinearWindingCircuit<'m>,
    core: Vec<Option<ResolvedConductorState>>,
    strands: Vec<SpatialWindingConductor>,
) -> SpatialNonlinearWindingCandidate<'m> {
    owner
        .step_with_conductors_candidate(
            owner.time_s() + DT,
            &[[0.0; 3]; 27],
            core,
            strands,
            &[VOLTS],
            controls(),
            &mut || false,
        )
        .unwrap()
}

/// The owner's heat records, booked to the three receivers by hand.
fn booked(step: &SpatialNonlinearWindingStep) -> [f64; 3] {
    let mut heat = [0.0; 3];
    for cell in &step.cells {
        heat[0] += cell.heat_j;
    }
    for cell in &step.field.cells {
        heat[1] += cell.heat_j;
    }
    for resistor in &step.graph_resistors {
        assert_eq!(resistor.branch, 2);
        heat[2] += resistor.heat_j;
    }
    heat
}

fn state(circuit: SpatialNonlinearWindingCircuit<'_>) -> SpatialWindingState<'_> {
    SpatialWindingState::new(circuit).unwrap()
}

#[test]
fn ports_follow_the_native_graph_and_no_port_tells_its_end() {
    let card = card(18.0);
    let p = cube(&card);
    let winding = component(&series(&p, &card), &p, &card);
    let ports = winding.ports();
    let names: Vec<_> = ports.iter().map(|p| p.id().as_str().to_owned()).collect();
    assert_eq!(names, ["source0", "heat0", "heat1", "heat2"]);
    assert_eq!(
        ports[0].carrier(),
        Carrier::Power(PortKind::ElectricalVoltageCurrent)
    );
    assert_eq!(ports[0].role(), PortRole::EffortIn);
    for port in &ports[1..] {
        assert_eq!(port.carrier(), Carrier::Heat);
        assert_eq!(port.role(), PortRole::EffortIn);
    }
    assert_eq!(winding.class(), Class::Smooth);
    for port in 0..ports.len() {
        assert_eq!(winding.end_of_step(port), EndOfStep::Untold);
    }
    assert_eq!(winding.receivers(), receivers(&p));
}

#[test]
fn a_step_is_the_owners_own_step_and_its_heat_leaves_through_the_ports() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let winding = component(&circuit, &p, &card);
    let layout = PortLayout::new(winding.ports()).unwrap();
    let kelvin = [300.0, 310.0, 320.0];
    let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
    let mut accepted = state(circuit.clone());
    assert!(accepted.step().is_none());
    let mut owner = circuit;
    for n in 1..=3 {
        let mut next = accepted.clone();
        let mut outputs = [0.0; 4];
        let proposal =
            with_cx(|cx| winding.propose(&accepted, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let expected = owner_step(
            &owner,
            core_at(&p, &card, 300.0),
            strands_at(&p, &card, 300.0, FILL),
        );
        let step = expected.step();
        let rebuilt = state(expected.circuit().clone());
        assert_eq!(next, rebuilt);
        // The state keeps the owner's own step.
        let kept = next.step().unwrap();
        assert_eq!(kept.heat_j.to_bits(), step.heat_j.to_bits());
        assert_eq!(kept.stored_after_j.to_bits(), step.stored_after_j.to_bits());
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        let heat = booked(step);
        for (k, q) in heat.iter().enumerate() {
            assert!(*q > 0.0, "receiver {k} received no heat at step {n}");
            assert_eq!(outputs[1 + k].to_bits(), (-q / DT).to_bits());
        }
        let returned: f64 = outputs[1..].iter().map(|rate| -rate * DT).sum();
        assert!((returned - step.heat_j).abs() <= receipt.energy_tolerance_j);
        assert_eq!(
            receipt.stored_change_j.to_bits(),
            step.stored_change_j.to_bits()
        );
        assert_eq!(
            receipt.dissipated_j.to_bits(),
            heat.iter().sum::<f64>().to_bits()
        );
        let entropy: f64 = heat.iter().zip(kelvin).map(|(q, t)| q / t).sum();
        assert_eq!(
            receipt.entropy_generated_j_per_k.to_bits(),
            entropy.to_bits()
        );
        assert!(
            (receipt.numerical_loss_j - step.numerical_dissipation_j).abs()
                <= receipt.energy_tolerance_j + TOLERANCE_J
        );
        assert!(step.numerical_dissipation_j > 0.0);
        assert!(receipt.residual <= receipt.tolerance);
        // The stored energy the states carry, bit for bit: the old state's
        // is the owner's start of step, and the new state's (the owner's end
        // of step) is the same sum recomputed from the new field.
        let old_j = winding.stored_energy_j(&accepted).unwrap();
        let new_j = winding.stored_energy_j(&next).unwrap();
        assert_eq!(new_j.to_bits(), step.stored_after_j.to_bits());
        let start = old_j - step.stored_before_j;
        let recomputed = new_j - rebuilt.stored_energy_j();
        assert_eq!(
            old_j.to_bits(),
            step.stored_before_j.to_bits(),
            "start {start:e} J"
        );
        assert_eq!(
            new_j.to_bits(),
            rebuilt.stored_energy_j().to_bits(),
            "recomputed {recomputed:e} J"
        );
        let audit = audit_step(
            winding.ports(),
            &layout,
            DT,
            &inputs,
            &outputs,
            &receipt,
            old_j,
            new_j,
        )
        .unwrap();
        println!(
            "component against owner: step {n}, current {:.9} A, heat {heat:?} J, stored \
             {new_j:.6e} J, numerical {:.6e} J; stored against the owner's start {start:.1e} \
             J, against the recomputed end {recomputed:.1e} J, floor {:.1e} J; balance defect \
             {:.1e} J, endpoint defect {:.1e} J",
            next.circuit().currents_a()[0],
            receipt.numerical_loss_j,
            receipt.energy_tolerance_j,
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        owner = expected.into_parts().0;
        accepted = next;
    }
}

#[test]
fn the_component_passes_the_conformance_kit() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let winding = component(&circuit, &p, &card);
    let state = state(circuit);
    let inputs = [VOLTS, 300.0, 310.0, 320.0];
    let mut stale = state.clone();
    let mut outputs = [0.0; 4];
    with_cx(|cx| winding.propose(&state, DT, &inputs, &mut stale, &mut outputs, cx)).unwrap();
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: DT,
        inputs: &inputs,
    };
    let observed =
        with_cx(|cx| check_component(&winding, &case, cx)).expect("component conformance");
    let Observed::Stepped { outputs, .. } = observed else {
        panic!("the winding shortened the step");
    };
    assert!(
        outputs[0] > 0.0 && outputs[1..].iter().all(|rate| *rate < 0.0),
        "{outputs:?}"
    );
}

#[test]
fn construction_and_steps_refuse_what_the_component_cannot_book() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let build = |receivers: Vec<SpatialHeatReceiver>, tolerance: f64| {
        SpatialWinding::new(
            &circuit,
            core_at(&p, &card, 300.0),
            strands_at(&p, &card, 300.0, FILL),
            receivers,
            controls(),
            tolerance,
        )
    };
    let refused = |result: Result<SpatialWinding<'_>, SpatialWindingError>, what: &str| match result
    {
        Err(SpatialWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    let mut missing = receivers(&p);
    missing.truncate(2);
    refused(
        build(missing, TOLERANCE_J),
        "a dissipating element is mapped to no receiver",
    );
    let core = cells(&p, false);
    let mut unmapped = receivers(&p);
    unmapped[1].cells.retain(|&cell| cell != core[0]);
    refused(
        build(unmapped, TOLERANCE_J),
        "a dissipating element is mapped to no receiver",
    );
    let mut twice = receivers(&p);
    twice[0].cells.push(core[0]);
    refused(
        build(twice, TOLERANCE_J),
        "a dissipating element is mapped to two receivers",
    );
    let mut empty = receivers(&p);
    empty.push(SpatialHeatReceiver::default());
    refused(build(empty, TOLERANCE_J), "a heat receiver maps no element");
    let mut inductor = receivers(&p);
    inductor[2].resistors.push(1);
    refused(
        build(inductor, TOLERANCE_J),
        "a receiver names a branch that is not a resistor",
    );
    let mut foreign = receivers(&p);
    foreign[0].cells.push(48);
    refused(
        build(foreign, TOLERANCE_J),
        "a receiver names a cell that carries no conductor",
    );
    refused(
        build(receivers(&p), 0.0),
        "energy tolerance must be positive and finite",
    );

    let winding = component(&circuit, &p, &card);
    let accepted = state(circuit.clone());
    let step = |inputs: &[f64], from: &SpatialWindingState<'_>| {
        let mut next = from.clone();
        let mut outputs = vec![0.0; inputs.len()];
        with_cx(|cx| winding.propose(from, DT, inputs, &mut next, &mut outputs, cx))
    };
    let invalid = |result: Result<Proposal, SpatialWindingError>, what: &str| match result {
        Err(SpatialWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    let ports = "positive duration, complete finite ports and positive temperatures required";
    invalid(step(&[VOLTS, 300.0, 0.0, 320.0], &accepted), ports);
    invalid(step(&[VOLTS, 300.0, 310.0], &accepted), ports);
    invalid(step(&[f64::NAN, 300.0, 310.0, 320.0], &accepted), ports);
    let topology = "the state's circuit is not this Component's topology";
    let other = state(series_on(
        &p,
        &card,
        300.0,
        Branch::Resistor { ohms: 1.0 },
        FILL,
    ));
    invalid(step(&[VOLTS, 300.0, 310.0, 320.0], &other), topology);
    // With fixed conductors, the graph resistor's value is the Component's.
    let hot = circuit
        .with_graph_resistors(vec![path_at(&card, 400.0)], &mut || false)
        .unwrap();
    invalid(step(&[VOLTS, 300.0, 310.0, 320.0], &state(hot)), topology);
    // A refused step leaves the accepted state as it was, and a state whose
    // core conductors alone differ is another state.
    assert_eq!(accepted, state(circuit));
    let warm = state(series_on(
        &p,
        &card,
        400.0,
        Branch::Capacitor { farads: 0.05 },
        FILL,
    ));
    assert_eq!(warm.circuit().currents_a(), accepted.circuit().currents_a());
    assert_ne!(warm, accepted);
}

#[test]
fn identity_words_tell_conductors_receivers_resistors_and_tolerances_apart() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let words = |winding: &SpatialWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let build = |core, strands, receivers, controls, tolerance| {
        words(
            &SpatialWinding::new(&circuit, core, strands, receivers, controls, tolerance).unwrap(),
        )
    };
    let base = words(&component(&circuit, &p, &card));
    assert_eq!(base, words(&component(&circuit, &p, &card)));
    let fixed_core = || core_at(&p, &card, 300.0);
    let fixed_strands = || strands_at(&p, &card, 300.0, FILL);
    let hot_core = build(
        core_at(&p, &card, 400.0),
        fixed_strands(),
        receivers(&p),
        controls(),
        TOLERANCE_J,
    );
    let hot_strands = build(
        fixed_core(),
        strands_at(&p, &card, 400.0, FILL),
        receivers(&p),
        controls(),
        TOLERANCE_J,
    );
    let thin = build(
        fixed_core(),
        strands_at(&p, &card, 300.0, 0.5),
        receivers(&p),
        controls(),
        TOLERANCE_J,
    );
    let mut swapped = receivers(&p);
    swapped.swap(0, 1);
    let moved = build(
        fixed_core(),
        fixed_strands(),
        swapped,
        controls(),
        TOLERANCE_J,
    );
    let mut looser = controls();
    looser.power_tolerance_w *= 2.0;
    let loose = build(
        fixed_core(),
        fixed_strands(),
        receivers(&p),
        looser,
        TOLERANCE_J,
    );
    let wider = build(
        fixed_core(),
        fixed_strands(),
        receivers(&p),
        controls(),
        2.0 * TOLERANCE_J,
    );
    for other in [&hot_core, &hot_strands, &thin, &moved, &loose, &wider] {
        assert_ne!(&base, other);
    }
    // Two graphs that differ in one resistor's value only, and two receiver
    // maps that differ only in which receiver heats each graph resistor.
    let split = |ohms: f64, branches: [usize; 2]| {
        let parallel = series_on(&p, &card, 300.0, Branch::Resistor { ohms }, FILL);
        let mut receivers = receivers(&p);
        receivers[1].resistors = vec![branches[0]];
        receivers[2].resistors = vec![branches[1]];
        words(
            &SpatialWinding::new(
                &parallel,
                fixed_core(),
                fixed_strands(),
                receivers,
                controls(),
                TOLERANCE_J,
            )
            .unwrap(),
        )
    };
    assert_ne!(split(1.0, [2, 3]), split(2.0, [2, 3]));
    assert_ne!(split(1.0, [2, 3]), split(1.0, [3, 2]));
}

/// Every cell's conductor from one card, at the port's temperature or, when
/// frozen, at that fixed temperature.
#[derive(Debug)]
struct CardLaw {
    card: MaterialCard,
    frozen: Option<f64>,
}

impl CardLaw {
    fn at(&self, kelvin: f64) -> Result<ResolvedConductorState, String> {
        resolve(&self.card, self.frozen.unwrap_or(kelvin))
    }
}

impl SpatialConductorLaw<ResolvedConductorState> for CardLaw {
    fn core(&self, _cell: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        self.at(kelvin)
    }
    fn strand(
        &self,
        _winding: usize,
        _cell: usize,
        kelvin: f64,
    ) -> Result<ResolvedConductorState, String> {
        self.at(kelvin)
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        let state = resolve(&self.card, 300.0).unwrap();
        words.push(state.resistivity_ohm_m().to_bits());
        for chunk in state.material().card_identity().0.chunks(8) {
            words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        words.push(self.frozen.map_or(0, f64::to_bits));
    }
}

/// The length of a unit-area path of the card with `GRAPH_OHM` at 300 K.
fn path_length(card: &MaterialCard) -> f64 {
    GRAPH_OHM / resolve(card, 300.0).unwrap().resistivity_ohm_m()
}

/// A uniform unit-area path of one card.
#[derive(Debug)]
struct PathLaw {
    law: CardLaw,
    length_m: f64,
}

impl ResistorLaw for PathLaw {
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String> {
        self.law
            .at(kelvin)?
            .into_uniform_conductor(self.length_m, 1.0)
            .map_err(|e| format!("{e:?}"))
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        self.law.identity_words(words);
        words.push(self.length_m.to_bits());
    }
}

/// Conductor laws by receiver, and resistor laws by branch.
type Laws = (
    Vec<Arc<dyn SpatialConductorLaw<ResolvedConductorState>>>,
    BTreeMap<usize, Arc<dyn ResistorLaw>>,
);

/// One card law per receiver and a path law for the graph resistor.
fn laws(card: &MaterialCard, frozen: Option<f64>) -> Laws {
    let cells: Arc<dyn SpatialConductorLaw<ResolvedConductorState>> = Arc::new(CardLaw {
        card: card.clone(),
        frozen,
    });
    let path: Arc<dyn ResistorLaw> = Arc::new(PathLaw {
        law: CardLaw {
            card: card.clone(),
            frozen,
        },
        length_m: path_length(card),
    });
    (
        vec![cells.clone(), cells.clone(), cells],
        BTreeMap::from([(2, path)]),
    )
}

fn following<'m>(
    circuit: &SpatialNonlinearWindingCircuit<'m>,
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialWinding<'m> {
    let (cells, paths) = laws(card, None);
    SpatialWinding::at_port_temperatures(
        circuit,
        cells,
        paths,
        receivers(p),
        controls(),
        TOLERANCE_J,
    )
    .unwrap()
}

/// The graph resistor of `card`'s path at `kelvin`.
fn path_at(card: &MaterialCard, kelvin: f64) -> GraphResistorMaterial {
    GraphResistorMaterial {
        branch: 2,
        conductor: resolve(card, kelvin)
            .unwrap()
            .into_uniform_conductor(path_length(card), 1.0)
            .unwrap(),
    }
}

#[test]
fn laws_give_each_element_at_its_ports_temperature_as_the_owner_would() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    // The resistivity moves by more than 5 % across 300-400 K.
    let ratio = resolve(&card, 300.0).unwrap().resistivity_ohm_m()
        / resolve(&card, 400.0).unwrap().resistivity_ohm_m();
    println!("laws: conductivity at 400 K is {ratio} of its value at 300 K");
    assert!(ratio < 0.95);
    let winding = following(&circuit, &p, &card);
    let layout = PortLayout::new(winding.ports()).unwrap();
    let mut accepted = state(circuit.clone());
    let mut owner = circuit;
    for (n, kelvin) in [
        [300.0, 350.0, 400.0],
        [330.0, 390.0, 310.0],
        [400.0, 300.0, 360.0],
    ]
    .into_iter()
    .enumerate()
    {
        let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
        let mut next = accepted.clone();
        let mut outputs = [0.0; 4];
        let proposal =
            with_cx(|cx| winding.propose(&accepted, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let path = path_at(&card, kelvin[2]);
        let ohms = path.conductor.resistance_ohm();
        let driven = owner
            .with_graph_resistors(vec![path], &mut || false)
            .unwrap();
        let expected = owner_step(
            &driven,
            core_at(&p, &card, kelvin[1]),
            strands_at(&p, &card, kelvin[0], FILL),
        );
        let step = expected.step();
        assert_eq!(next, state(expected.circuit().clone()));
        assert_eq!(
            next.circuit().graph().branches[2].2,
            Branch::Resistor { ohms }
        );
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        let heat = booked(step);
        for (k, q) in heat.iter().enumerate() {
            assert_eq!(outputs[1 + k].to_bits(), (-q / DT).to_bits());
        }
        assert_eq!(
            receipt.stored_change_j.to_bits(),
            step.stored_change_j.to_bits()
        );
        let entropy: f64 = heat.iter().zip(kelvin).map(|(q, t)| q / t).sum();
        assert_eq!(
            receipt.entropy_generated_j_per_k.to_bits(),
            entropy.to_bits()
        );
        let old_j = winding.stored_energy_j(&accepted).unwrap();
        let new_j = winding.stored_energy_j(&next).unwrap();
        let audit = audit_step(
            winding.ports(),
            &layout,
            DT,
            &inputs,
            &outputs,
            &receipt,
            old_j,
            new_j,
        )
        .unwrap();
        println!(
            "laws against owner: step {}, ports at {kelvin:?} K, current {:.9} A, graph \
             resistor {ohms:.9} ohm, heat {heat:?} J, balance defect {:.1e} J, endpoint \
             defect {:.1e} J",
            n + 1,
            next.circuit().currents_a()[0],
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        owner = expected.into_parts().0;
        accepted = next;
    }
}

#[test]
fn a_law_that_ignores_temperature_steps_as_fixed_conductors_do() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let (cells, paths) = laws(&card, Some(300.0));
    let frozen = SpatialWinding::at_port_temperatures(
        &circuit,
        cells,
        paths,
        receivers(&p),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    // The fixed Component steps the circuit whose graph resistor already
    // carries the law's resistance.
    let bound = circuit
        .with_graph_resistors(vec![path_at(&card, 300.0)], &mut || false)
        .unwrap();
    let fixed = component(&bound, &p, &card);
    let mut a = state(circuit);
    let mut b = state(bound);
    for kelvin in [[300.0, 350.0, 400.0], [400.0, 330.0, 310.0]] {
        let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
        let (mut na, mut nb) = (a.clone(), b.clone());
        let (mut oa, mut ob) = ([0.0; 4], [0.0; 4]);
        with_cx(|cx| frozen.propose(&a, DT, &inputs, &mut na, &mut oa, cx)).unwrap();
        with_cx(|cx| fixed.propose(&b, DT, &inputs, &mut nb, &mut ob, cx)).unwrap();
        assert_eq!(na, nb);
        assert_eq!(oa.map(f64::to_bits), ob.map(f64::to_bits));
        (a, b) = (na, nb);
    }
}

#[test]
fn the_component_with_laws_passes_the_conformance_kit() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let winding = following(&circuit, &p, &card);
    let state = state(circuit);
    let inputs = [VOLTS, 300.0, 350.0, 400.0];
    let mut stale = state.clone();
    let mut outputs = [0.0; 4];
    with_cx(|cx| winding.propose(&state, DT, &inputs, &mut stale, &mut outputs, cx)).unwrap();
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: DT,
        inputs: &inputs,
    };
    with_cx(|cx| check_component(&winding, &case, cx)).expect("component conformance");
}

#[test]
fn laws_refuse_what_they_cannot_answer_and_constructions_they_cannot_drive() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let (conductors, paths) = laws(&card, None);
    let build = |circuit: &SpatialNonlinearWindingCircuit<'_>,
                 conductors: Vec<Arc<dyn SpatialConductorLaw<ResolvedConductorState>>>,
                 paths: BTreeMap<usize, Arc<dyn ResistorLaw>>,
                 receivers: Vec<SpatialHeatReceiver>| {
        SpatialWinding::at_port_temperatures(
            circuit,
            conductors,
            paths,
            receivers,
            controls(),
            TOLERANCE_J,
        )
        .map(|_| ())
    };
    let refused = |result: Result<(), SpatialWindingError>, what: &str| match result {
        Err(SpatialWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    refused(
        build(
            &circuit,
            conductors[..2].to_vec(),
            paths.clone(),
            receivers(&p),
        ),
        "one conductor law per heat receiver required",
    );
    let mut inductor = paths.clone();
    inductor.insert(1, paths[&2].clone());
    refused(
        build(&circuit, conductors.clone(), inductor, receivers(&p)),
        "a resistor law names a branch that is not a resistor",
    );
    // The owner rebinds every graph resistor or none.
    let parallel = series_on(&p, &card, 300.0, Branch::Resistor { ohms: 1.0 }, FILL);
    let mut both = receivers(&p);
    both[2].resistors.push(3);
    refused(
        build(&parallel, conductors.clone(), paths.clone(), both),
        "every graph resistor or none has a law",
    );

    // The card's validity ends at 500 K. A law asked above it refuses with
    // the element, its receiver and the temperature; nothing is clamped.
    let winding = following(&circuit, &p, &card);
    let accepted = state(circuit.clone());
    let step = |winding: &SpatialWinding<'_>, inputs: &[f64]| {
        let mut next = accepted.clone();
        let mut outputs = vec![0.0; inputs.len()];
        with_cx(|cx| winding.propose(&accepted, DT, inputs, &mut next, &mut outputs, cx))
    };
    let strand = cells(&p, true)[0];
    let core = cells(&p, false)[0];
    for (inputs, receiver, element) in [
        (
            [VOLTS, 520.0, 300.0, 310.0],
            0,
            SpatialLawElement::Strand {
                winding: 0,
                cell: strand,
            },
        ),
        (
            [VOLTS, 300.0, 520.0, 310.0],
            1,
            SpatialLawElement::Core(core),
        ),
        (
            [VOLTS, 300.0, 310.0, 520.0],
            2,
            SpatialLawElement::Resistor(2),
        ),
    ] {
        match step(&winding, &inputs) {
            Err(SpatialWindingError::Law {
                receiver: r,
                element: e,
                kelvin,
                ..
            }) if r == receiver && e == element && kelvin.to_bits() == 520.0_f64.to_bits() => {}
            other => panic!("expected the refusal of {element:?}, got {other:?}"),
        }
    }
    // A law of another card meets the owner's own continuity refusal.
    let (other, _) = laws(&self::card(9.0), None);
    let foreign = SpatialWinding::at_port_temperatures(
        &circuit,
        other,
        paths,
        receivers(&p),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    match step(&foreign, &[VOLTS, 300.0, 310.0, 320.0]) {
        Err(SpatialWindingError::Owner(error)) => println!("laws: another card: {error}"),
        other => panic!("expected the owner's refusal, got {other:?}"),
    }
    assert_eq!(accepted, state(circuit));
}

#[test]
fn identity_words_tell_laws_and_cells_apart() {
    let card = card(18.0);
    let p = cube(&card);
    let circuit = series(&p, &card);
    let words = |winding: &SpatialWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let at = |circuit: &SpatialNonlinearWindingCircuit<'_>, (cells, paths): Laws| {
        words(
            &SpatialWinding::at_port_temperatures(
                circuit,
                cells,
                paths,
                receivers(&p),
                controls(),
                TOLERANCE_J,
            )
            .unwrap(),
        )
    };
    let base = at(&circuit, laws(&card, None));
    assert_eq!(base, at(&circuit, laws(&card, None)));
    assert_ne!(base, words(&component(&circuit, &p, &card)));
    assert_ne!(base, at(&circuit, laws(&card, Some(300.0))));
    let other = self::card(9.0);
    let (cells, paths) = laws(&card, None);
    let (foreign_cells, foreign_paths) = laws(&other, None);
    assert_ne!(base, at(&circuit, (foreign_cells, paths)));
    assert_ne!(base, at(&circuit, (cells, foreign_paths)));
    // The same graph with thinner strands.
    let thin = series_on(&p, &card, 300.0, Branch::Capacitor { farads: 0.05 }, 0.5);
    assert_ne!(base, at(&thin, laws(&card, None)));
}
