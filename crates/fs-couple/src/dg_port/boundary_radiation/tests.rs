use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::bake::boundary::tests::baked_fixture;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_vfit::mimo::load::PreparedMimoLoad;

#[path = "../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(11);

fn duct(source: &BoundaryMimoLoad, gas: &GasState, area_factor: f64) -> DgDuct {
    let area = source.record().aperture_support_areas_m2()[1] * area_factor;
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: det::sqrt(area / core::f64::consts::PI),
                length: 0.05,
            }],
        },
        gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(INLET),
        DgBoundary::Port(OUTLET),
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: SIDE,
        downstream_element: 1,
    }])
    .unwrap()
}

fn bindings(source: &BoundaryMimoLoad) -> Vec<(PortId, SourceId)> {
    // The physical matrix order deliberately differs from DG's inlet/outlet/side.
    vec![
        (SIDE, source.record().sources()[0].clone()),
        (OUTLET, source.record().sources()[1].clone()),
    ]
}

pub(super) fn load() -> DgPortLoad {
    let source = baked_fixture();
    DgPortLoad::new_with_boundary_radiation(
        duct(&source, source.record().gas(), 1.0),
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        bindings(&source),
        source,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn material_self_and_mutual_reaction_advance_the_original_work_conjugate_loop() {
    let source = baked_fixture();
    let original = duct(&source, source.record().gas(), 1.0);
    let dt = source.load().time_step_s();
    let identity = source.record().realization_identity();
    let mut reference = DgPortLoad::new_with_terminal_ports(
        original.clone(),
        dt,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET],
        source.load().clone(),
    )
    .unwrap();
    // An independent self-only control exposes accidental loss of every mutual
    // term. Diagonal positive-real channels retain the original scalar owners.
    let mut self_only = source.load().source_model().clone();
    for index in [1, 2] {
        self_only.d[index] = 0.0;
        self_only.e[index] = 0.0;
        for residue in &mut self_only.residues {
            residue[index] = C64::ZERO;
        }
    }
    let diagonal = PreparedMimoLoad::new(&self_only, dt, *source.load().limits()).unwrap();
    let mut uncoupled = DgPortLoad::new_with_terminal_ports(
        original.clone(),
        dt,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET],
        diagonal,
    )
    .unwrap();
    let mut actual = DgPortLoad::new_with_boundary_radiation(
        original,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        bindings(&source),
        source,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(actual.terminal_ports(), &[SIDE, OUTLET]);
    assert!(actual.radiation_record().is_none());
    let (mut inlet_work, mut exterior_work, mut gap, mut scale) = (0.0, 0.0, 0.0, 0.0);
    for sample in 0..96 {
        let flow = if sample < 12 { 1e-9 } else { 0.0 };
        let candidate = actual.preview_step(flow).unwrap();
        inlet_work += candidate.inlet_work_j();
        for frame in candidate.terminal_frames() {
            exterior_work += frame.load_work_j;
            let work_scale = frame.port_work_j.abs() + frame.load_work_j.abs();
            assert!(frame.port_work_residual_j.abs() <= 2e-10 * work_scale.max(1e-300));
        }
        assert_eq!(
            candidate.accept(),
            reference.preview_step(flow).unwrap().accept()
        );
        uncoupled.preview_step(flow).unwrap().accept();
        assert_eq!(actual.state(), reference.state());
        assert_eq!(
            actual.terminal().unwrap().pole_states(),
            reference.terminal().unwrap().pole_states()
        );
        for (full, diagonal) in actual
            .last_terminal_frames()
            .iter()
            .zip(uncoupled.last_terminal_frames())
        {
            gap += (full.outward_flow_m3_s - diagonal.outward_flow_m3_s).powi(2);
            scale += full.outward_flow_m3_s.powi(2);
        }
    }
    assert!(
        gap > 1e-8 * scale,
        "native material geometry must retain audible-loop mutual reaction"
    );
    let balance =
        actual.stored_energy_j() + actual.cumulative_dissipated_j() + exterior_work - inlet_work;
    assert!(balance.abs() <= 2e-9 * (inlet_work.abs() + exterior_work.abs()).max(1e-300));
    assert_eq!(
        actual
            .boundary_radiation_record()
            .unwrap()
            .realization_identity(),
        identity
    );
}

#[test]
fn original_aperture_names_area_pressure_and_all_ports_are_required() {
    let source = baked_fixture();
    let original = duct(&source, source.record().gas(), 1.0);
    let valid = bindings(&source);
    let mut reversed_names = valid.clone();
    reversed_names[0].1 = valid[1].1.clone();
    reversed_names[1].1 = valid[0].1.clone();
    for invalid_bindings in [
        vec![valid[0].clone()],
        vec![valid[0].clone(), valid[0].clone()],
        vec![(INLET, valid[0].1.clone()), valid[1].clone()],
        vec![(PortId(999), valid[0].1.clone()), valid[1].clone()],
        reversed_names,
    ] {
        assert!(
            DgPortLoad::new_with_boundary_radiation(
                original.clone(),
                1e-11,
                MAX_DG_PROPAGATOR_DOFS,
                invalid_bindings,
                source.clone(),
                &CancelGate::new_clock_free(),
            )
            .is_err()
        );
    }
    let different_pressure = GasState::try_new_moist_air(293.15, 90_000.0, 0.0).unwrap();
    for mismatched in [
        duct(&source, source.record().gas(), 1.2),
        duct(&source, &different_pressure, 1.0),
    ] {
        assert!(
            DgPortLoad::new_with_boundary_radiation(
                mismatched,
                1e-11,
                MAX_DG_PROPAGATOR_DOFS,
                valid.clone(),
                source.clone(),
                &CancelGate::new_clock_free(),
            )
            .is_err()
        );
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        DgPortLoad::new_with_boundary_radiation(
            original,
            1e-11,
            MAX_DG_PROPAGATOR_DOFS,
            valid,
            source,
            &gate,
        ),
        Err(AcousticRealizeError::Cancelled)
    ));
}

#[test]
fn finite_branch_radiates_its_outer_flow_and_keeps_original_storage() {
    let source = baked_fixture();
    let radius = det::sqrt(source.record().aperture_support_areas_m2()[0] / core::f64::consts::PI);
    let branch = DgBranch {
        port: SIDE,
        duct: DgDuct::new(
            &Duct {
                segments: vec![Segment::Cylinder {
                    radius,
                    length: 0.03,
                }],
            },
            source.record().gas(),
            DgMesh {
                elements: 1,
                order: 2,
            },
            DgBoundary::Port(PortId(0)),
            DgBoundary::Port(PortId(1)),
        )
        .unwrap(),
    };
    let original = duct(&source, source.record().gas(), 1.0);
    let mut actual = DgPortLoad::new_with_boundary_radiation_branches(
        original,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        bindings(&source),
        source,
        vec![branch],
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let (mut junction_difference, mut exterior_scale) = (0.0, 0.0);
    for sample in 0..48 {
        actual
            .preview_step(if sample < 4 { 1e-9 } else { 0.0 })
            .unwrap()
            .accept();
        let branch = &actual.last_branch_frames()[0];
        let exterior = &actual.last_terminal_frames()[0];
        assert_eq!(exterior.port, SIDE);
        // The native branch's distal work is the negative exterior work;
        // substituting its inlet flow would miss retained compressive storage.
        let distal_work = branch.duct.ledger.port_work_step_j - branch.branch_work_j;
        let scale = branch.duct.ledger.port_work_step_j.abs()
            + branch.branch_work_j.abs()
            + exterior.port_work_j.abs();
        assert!((distal_work + exterior.port_work_j).abs() <= 2e-10 * scale.max(1e-300));
        junction_difference += (branch.junction_flow_m3_s - exterior.outward_flow_m3_s).powi(2);
        exterior_scale += exterior.outward_flow_m3_s.powi(2);
    }
    assert!(junction_difference > 1e-8 * exterior_scale);
    assert!(actual.branches()[0].state().step_count() > 0);
}

#[test]
fn continuing_interior_gas_cannot_redefine_the_retained_diffraction_exterior() {
    use crate::dg_port::DgGasRefresh;
    use fs_matdb::QueryPoint;
    use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
    use fs_material::gas_species::input::ThermochemicalPhase;
    use fs_material::thermochemical_gas::packs::from_packs;
    use fs_thermochem::{Composition, CompositionBasis};

    let field = |pressure| {
        let source = gas_support::packs("boundary-loop-gas", 0.028, 3.5, 200.0);
        let composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
        let model = from_packs(&[source.binding()], composition.clone(), None).unwrap();
        AxialGasField::new(
            vec![0.0, 0.05],
            pressure,
            vec![
                PrimitiveGasState {
                    temperature_k: 330.0,
                    composition
                };
                2
            ],
            model,
            QueryPoint::new(),
            ThermochemicalPhase::FrozenGas,
        )
        .unwrap()
    };
    let mut actual = load();
    actual.preview_step(1e-9).unwrap().accept();
    let before = actual.clone();
    assert!(
        actual
            .rebind_gas_field_controlled(&field(90_000.0), DgGasRefresh::Lossless, 1.0, &mut || {
                false
            },)
            .is_err()
    );
    assert_eq!(actual.state(), before.state());
    assert_eq!(
        actual.duct().operator_identity(),
        before.duct().operator_identity()
    );
    assert_eq!(
        actual.terminal().unwrap().pole_states(),
        before.terminal().unwrap().pole_states()
    );
    assert_eq!(actual.last_terminal_frames(), before.last_terminal_frames());
    let pressure = before.boundary_radiation_record().unwrap().gas().pressure;
    actual
        .rebind_gas_field_controlled(&field(pressure), DgGasRefresh::Lossless, 1.0, &mut || false)
        .unwrap();
    assert_eq!(actual.state().step_count(), before.state().step_count());
    assert_eq!(
        actual.terminal().unwrap().pole_states(),
        before.terminal().unwrap().pole_states()
    );
    assert_eq!(
        actual
            .boundary_radiation_record()
            .unwrap()
            .realization_identity(),
        before
            .boundary_radiation_record()
            .unwrap()
            .realization_identity()
    );
}
