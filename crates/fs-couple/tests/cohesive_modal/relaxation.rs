//! G2 native histories, full mechanical response and transactional refusals.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::cohesive::relaxation::{
    ViscoelasticBondModalFrame, ViscoelasticBondModalSystem,
};
use fs_matdb::QueryPoint;
use fs_material::cohesive::relaxation::source::{
    ResolvedViscoelasticBond, resolve_viscoelastic_bond,
};
use fs_material::state_point::MaterialPropertySelection;

#[path = "../viscoelastic_bond_source/common.rs"]
mod material;

#[path = "relaxation/temperature.rs"]
mod temperature;

fn source() -> ResolvedViscoelasticBond {
    resolve_viscoelastic_bond(
        &material::card(false, false),
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        material::limits(),
    )
    .unwrap()
}
fn build(
    source: &ResolvedViscoelasticBond,
    relative: [f64; 3],
    cross: bool,
) -> ViscoelasticBondModalSystem<'_> {
    ViscoelasticBondModalSystem::from_source(
        network(relative, 1.0, cross),
        source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        AREA,
        columns(),
        config(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn native_states(system: &ViscoelasticBondModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn native_momentum(system: &ViscoelasticBondModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        det::sqrt(MA) * system.components()[0].states()[i].velocity_m_sqrt_kg_per_s
            + det::sqrt(MB) * system.components()[1].states()[i].velocity_m_sqrt_kg_per_s
    })
}
fn point_observation(system: &ViscoelasticBondModalSystem<'_>) -> ([f64; 3], [f64; 3], f64, f64) {
    (
        system.point().jump_m(),
        system.point().traction_pa(),
        system.point().stored_energy_j_m2(),
        system.point().elapsed_s(),
    )
}

#[test]
fn original_two_body_bond_relaxes_with_conserved_momenta_and_complete_energy() {
    let source = source();
    let mut system = build(&source, [0.002, -0.001, 0.0005], false);
    assert!(core::ptr::eq(system.source(), &source));
    let initial = system.total_energy_j().unwrap();
    let momenta = native_momentum(&system);
    let mut heat = 0.0;
    let mut native_loss = 0.0;
    let mut numerical_loss = 0.0;
    let mut different_mean = false;
    for _ in 0..1000 {
        let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
        heat += frame.relaxation_heat_j;
        native_loss += frame.network_dissipation_j;
        numerical_loss += frame.numerical_dissipation_j;
        different_mean |=
            (frame.mean_force_n[0] - AREA * frame.endpoint_traction_pa[0]).abs() > 1e-6;
        assert!(frame.relaxation_heat_j >= 0.0);
        assert_eq!(frame.external_work_j, 0.0);
        assert_eq!(frame.prescribed_work_j, 0.0);
        assert!(frame.bond_work_residual_j.abs() <= frame.bond_work_tolerance_j);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        for i in 0..3 {
            close(native_momentum(&system)[i], momenta[i], 3e-13);
        }
        close(
            system.total_energy_j().unwrap() + heat + native_loss + numerical_loss,
            initial,
            3e-11,
        );
    }
    assert!(heat > 1e-8 && different_mean);
    assert_eq!(system.samples_rendered(), 1000);
    close(system.point().elapsed_s(), 1000.0 / f64::from(RATE), 1e-13);
}

#[test]
fn hereditary_mean_traction_drives_full_cross_axis_native_motion() {
    let source = source();
    let mut system = build(&source, [0.0; 3], true);
    let compliance = system.compliance_m_per_n();
    assert!(compliance[0][1].abs() > 1e-10 && compliance[0][2].abs() > 1e-10);
    let target = [1e-5, -5e-6, 2e-6];
    let prepared = source.law().prepare(system.sample_period_s()).unwrap();
    let expected = system.point().trial(&prepared, target, 1e-8).unwrap();
    let force = expected.mean_traction_pa.map(|v| AREA * v);
    let rate = core::array::from_fn(|i| {
        (target[i] + (0..3).map(|j| compliance[i][j] * force[j]).sum::<f64>())
            / system.sample_period_s()
    });
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0; 6], rate).unwrap();
    for i in 0..3 {
        close(frame.jump_m[i], target[i], 5e-13);
        close(frame.mean_force_n[i], force[i], 3e-9);
    }
    close(
        frame.relaxation_heat_j,
        AREA * expected.viscous_dissipated_j_m2,
        1e-13,
    );
    close(
        system.total_energy_j().unwrap() - before
            + frame.network_dissipation_j
            + frame.numerical_dissipation_j
            + frame.relaxation_heat_j,
        frame.prescribed_work_j,
        3e-12,
    );
}

#[test]
fn receiver_refusal_late_cancellation_and_material_bounds_preserve_branch_memory() {
    let source = source();
    let mut system = build(&source, [0.001, -0.0005, 0.0003], false);
    let mut replay = build(&source, [0.001, -0.0005, 0.0003], false);
    for _ in 0..20 {
        for model in [&mut system, &mut replay] {
            model.step(&[0.0; 6], [0.0; 3]).unwrap();
        }
    }
    let state = native_states(&system);
    let point = point_observation(&system);
    let last = system.last_frame().copied();
    let sample = system.samples_rendered();
    let gate = CancelGate::new();
    let refused: Result<(), CohesiveModalError> = system.step_with_admission(
        &[0.0; 6],
        [0.002, -0.001, 0.0005],
        &gate,
        |frame, candidate| {
            assert_eq!(frame.sample, sample + 1);
            assert_ne!(candidate[0].states(), &state[..3]);
            Err(CohesiveModalError::Network(
                ModalCouplingError::EnergyBalance {
                    residual_j: 1.0,
                    tolerance_j: 0.0,
                },
            ))
        },
    );
    assert!(refused.is_err());
    let cancelled: Result<ViscoelasticBondModalFrame, CohesiveModalError> = system
        .step_with_admission(&[0.0; 6], [0.002, -0.001, 0.0005], &gate, |frame, _| {
            gate.request();
            Ok(*frame)
        });
    assert!(cancelled.is_err());
    // The original point's rate / separation bounds refuse this proposed path.
    assert!(system.step(&[0.0; 6], [10.0, 0.0, 0.0]).is_err());
    assert_eq!(native_states(&system), state);
    assert_eq!(point_observation(&system), point);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.samples_rendered(), sample);
    // A full continuation proves hidden Prony histories survived every refusal.
    for _ in 0..50 {
        assert_eq!(
            *system.step(&[0.0; 6], [0.0; 3]).unwrap(),
            *replay.step(&[0.0; 6], [0.0; 3]).unwrap()
        );
        assert_eq!(native_states(&system), native_states(&replay));
    }
    assert!(
        ViscoelasticBondModalSystem::from_source(
            network([0.0; 3], 1.0, false),
            &source,
            [
                &source.surfaces()[1].material,
                &source.surfaces()[0].material
            ],
            AREA,
            columns(),
            config(),
            &CancelGate::new(),
        )
        .is_err()
    );
}
