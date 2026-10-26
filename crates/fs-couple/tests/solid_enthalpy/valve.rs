//! G1/G3/G4: finite local heat reaches the actual tangent and pressure port.
use super::*;
use fs_couple::solid_enthalpy::valve::{HeatedValveError, HeatedValveReduction};
use fs_scenario::reed_source::BeatingReedSource;
use fs_solid::reduce::thermal::ThermalValveEquilibrium;
use fs_solid::reduce::{OrificeSpec, SingleModeValveAdmission};

#[path = "valve/motion.rs"]
mod motion;

#[path = "valve/gas_motion.rs"]
mod gas_motion;

const CUBE: [[f64; 3]; 8] = [
    [0., 0., 0.],
    [0., 0., 1.],
    [0., 1., 0.],
    [0., 1., 1.],
    [1., 0., 0.],
    [1., 0., 1.],
    [1., 1., 0.],
    [1., 1., 1.],
];
const CUBE_TETS: [[usize; 4]; 6] = [
    [0, 4, 6, 7],
    [0, 6, 2, 7],
    [0, 2, 3, 7],
    [0, 3, 1, 7],
    [0, 1, 5, 7],
    [0, 5, 4, 7],
];
const FACE: [usize; 4] = [4, 5, 6, 7];
const FORCE: [(usize, f64); 1] = [(21, -0.1)];
fn reduction<'a>(fixed: &'a [usize], pressure: &'a [[f64; 3]]) -> HeatedValveReduction<'a> {
    HeatedValveReduction {
        fixed_dofs: fixed,
        orifice: OrificeSpec {
            face_nodes: &FACE,
            opening_axis: [1., 0., 0.],
            width_axis: [0., 1., 0.],
            opposing_plane_offset_m: 0.5,
        },
        equilibrium_nodal_forces_n: &FORCE,
        nodal_load_per_pa_m2: pressure,
        window_hz: (0.01, 100.0),
        retained_compliance_floor: 0.999,
        loss_factor: 0.02,
        source_id: "manufactured finite-heat cube",
        static_solve: TetStaticSolveConfig::default(),
        maximum_nodal_displacement_m: 0.1,
    }
}

#[test]
fn finite_heat_capacity_and_expansion_change_the_analytic_valve_without_changing_mass() {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..24).filter(|&d| d != 21).collect();
        let mut pressure = [[0.0; 3]; 8];
        pressure[7][0] = 1.0;
        let input = reduction(&fixed, &pressure);
        let context = point(300.0);
        let mut identities = Vec::new();
        for (alpha, cp) in [(1e-5, 10.0), (2e-5, 10.0), (1e-5, 20.0)] {
            let (card, caloric) = material(1000., alpha, cp, 0., false);
            let cells: Vec<_> = (0..6).map(|id| (id, 200.0)).collect();
            let thermal = initial(caloric.curve(), &cells, 1.0);
            let source = SolidEnthalpy::try_new(
                &thermal,
                &CUBE,
                &CUBE_TETS,
                (0..6)
                    .map(|id| super::source(id, &card, caloric.curve(), &context))
                    .collect(),
                TetAssemblyBudget::standard(),
                0.1,
                config().mesh_frame,
                cx,
            )
            .unwrap();
            let next = source
                .advance_pressure_valve_to(&thermal, 1.0, thermal_config(), &input, cx)
                .unwrap();
            let rise = 200.0 / cp;
            let strain = alpha * (rise + 0.005 * rise * rise);
            let e = 1000.0 * (1.0 - 0.002 * rise);
            // Direct integrals on the six affine tetrahedra, one free x DOF:
            // K_e=2E/3, M=rho/10, f_thermal=E*epsilon/3,
            // K_g=E*(u-3*epsilon)/3, mean aperture translation=u/4.
            let u = 0.5 * strain + FORCE[0].1 / (2.0 * e / 3.0);
            let stiffness = 2.0 * e / 3.0 + e * (u - 3.0 * strain) / 3.0;
            let card = next.valve.valve().card();
            near(next.thermal_audit.external_heat_j, 1200., 1e-8);
            near(next.thermal_audit.stored_energy_change_j, 1200., 1e-8);
            near(
                next.thermal.temperature_k(CellId(0)).unwrap(),
                300. + rise,
                1e-10,
            );
            near(next.valve.equilibrium().displacement_m()[7][0], u, 1e-10);
            near(card.total_mass_kg, 6.0, 1e-12);
            near(
                card.modes[0].frequency_hz,
                (stiffness / 0.6).sqrt() / core::f64::consts::TAU,
                1e-8,
            );
            near(card.modes[0].mass_kg, 9.6, 1e-8);
            near(card.modes[0].stiffness_n_m, 16.0 * stiffness, 1e-5);
            near(card.rest_gap_m, 0.5 + u.min(0.0), 1e-12);
            near(
                next.valve
                    .valve()
                    .pressure_compliance()
                    .exact_volume_m3_per_pa(),
                1.0 / stiffness,
                1e-11,
            );
            near(
                next.valve
                    .valve()
                    .pressure_compliance()
                    .exact_opening_m_per_pa(),
                0.25 / stiffness,
                1e-11,
            );
            let ThermalValveEquilibrium::Tet4 { tangent, .. } = next.valve.equilibrium() else {
                panic!("P1 inventory must keep P1 equilibrium");
            };
            assert_eq!(tangent.pencil().total_mass_kg, card.total_mass_kg);
            for sample in next.state.samples() {
                assert_eq!(
                    sample.elastic.card_identity(),
                    sample.reference_elastic.card_identity()
                );
                near(sample.reference_mass_kg, 1., 1e-12);
            }
            let reed = next
                .valve
                .compile_beating_reed(
                    0.01,
                    0.001,
                    &SingleModeValveAdmission {
                        minimum_opening_compliance_fraction: 0.99,
                        minimum_pressure_compliance_fraction: 0.99,
                        maximum_static_opening_relative_error: 1e-8,
                    },
                )
                .unwrap();
            near(
                reed.closing_pressure_pa,
                reed.stiffness_n_m * reed.rest_opening_m / 4.0,
                1e-8,
            );
            let repeat = source
                .advance_pressure_valve_to(&thermal, 1.0, thermal_config(), &input, cx)
                .unwrap();
            assert_eq!(repeat.valve.identity(), next.valve.identity());
            identities.push(next.valve.identity());
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 0.0);
        }
        assert_ne!(identities[0], identities[1]);
        assert_ne!(identities[0], identities[2]);
    });
}

#[test]
fn failed_aperture_or_material_path_does_not_publish_the_thermal_debit() {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..24).filter(|&d| d != 21).collect();
        let mut pressure = [[0.0; 3]; 8];
        pressure[7][0] = 1.0;
        let context = point(300.0);
        for split in [false, true] {
            let (card, caloric) = material(1000., 1e-5, 10., 0., split);
            let cells: Vec<_> = (0..6).map(|id| (id, 1000.0)).collect();
            let thermal = initial(caloric.curve(), &cells, 1.0);
            let source = SolidEnthalpy::try_new(
                &thermal,
                &CUBE,
                &CUBE_TETS,
                (0..6)
                    .map(|id| super::source(id, &card, caloric.curve(), &context))
                    .collect(),
                TetAssemblyBudget::standard(),
                0.1,
                config().mesh_frame,
                cx,
            )
            .unwrap();
            let mut input = reduction(&fixed, &pressure);
            if !split {
                input.orifice.opposing_plane_offset_m = 2.0;
            }
            let failed =
                source.advance_pressure_valve_to(&thermal, 1.0, thermal_config(), &input, cx);
            if split {
                assert!(matches!(
                    failed,
                    Err(HeatedValveError::Solid(SolidEnthalpyError::Material(_)))
                ));
            } else {
                assert!(matches!(failed, Err(HeatedValveError::Reduction(_))));
            }
            assert_eq!(thermal.time_s(), 0.0);
            assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
            input = reduction(&fixed, &pressure);
            assert!(
                source
                    .advance_pressure_valve_to(&thermal, 0.125, thermal_config(), &input, cx)
                    .is_ok()
            );
        }
    });
}

#[test]
fn cancelled_heated_valve_keeps_the_accepted_inventory() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (card, caloric) = material(1000., 1e-5, 10., 0., false);
        let thermal = initial(
            caloric.curve(),
            &(0..6).map(|id| (id, 200.0)).collect::<Vec<_>>(),
            1.0,
        );
        let context = point(300.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &CUBE,
            &CUBE_TETS,
            (0..6)
                .map(|id| super::source(id, &card, caloric.curve(), &context))
                .collect(),
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            &cx,
        )
        .unwrap();
        let fixed: Vec<_> = (0..24).filter(|&d| d != 21).collect();
        let mut pressure = [[0.0; 3]; 8];
        pressure[7][0] = 1.0;
        gate.request();
        assert!(matches!(
            source.advance_pressure_valve_to(
                &thermal,
                1.0,
                thermal_config(),
                &reduction(&fixed, &pressure),
                &cx
            ),
            Err(HeatedValveError::Solid(SolidEnthalpyError::Cancelled))
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    });
}
