//! G1 native membrane/heat oracle, G3 layered thermal reversal and G4 atomicity.
//! These authored material numbers exercise the model, not measured calibration.
#[path = "shell/dynamics.rs"]
mod dynamics;
#[path = "shell/creep.rs"]
mod creep;
#[path = "shell/contact.rs"]
mod contact;
#[path = "shell/hill.rs"]
mod hill;

use super::*;
use fs_couple::shell_plastic_enthalpy::{
    ShellThermoplasticDiscretization, ShellThermoplasticError, ShellThermoplasticLayer,
    ThermoplasticShell,
};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{PlaneStressLayer, SectionBudget, ThicknessQuadrature};
use fs_material::plastic::thermal::ThermoplasticSource;
use fs_plate::shell::ShellMesh;
use fs_solid::plastic_shell::PlasticShellBudget;

const THICKNESS: f64 = 0.001;
const AREA: f64 = L * L / 2.0;
const VOLUME: f64 = AREA * THICKNESS;

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [L, 0.0, 0.0], [0.0, L, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn discretization() -> ShellThermoplasticDiscretization {
    ShellThermoplasticDiscretization {
        coordinate_length_m: L,
        rotation_limit_rad: 0.04,
        drilling_alpha: 0.02,
        quadrature: ThicknessQuadrature::Gauss3,
        reduction: ReductionConfig {
            max_iterations: 32,
            stress_tolerance_pa: 1e-4,
            strain_limit: 0.04,
        },
        section_budget: SectionBudget {
            maximum_layers: 4,
            maximum_quadrature_points: 48,
        },
        shell_budget: PlasticShellBudget::default(),
    }
}

fn thermal<'a>(fixture: &'a Fixture, powers: &[f64], conductance: f64) -> DepositionSimulation<'a> {
    let volume = VOLUME / powers.len() as f64;
    let mut events: Vec<_> = powers
        .iter()
        .enumerate()
        .map(|(i, &source_w)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(100 + i as u64),
                material_id: MaterialId(0),
                segment_id: i as u64,
                body_id: 7,
                mass_kg: RHO * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L, L, THICKNESS],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        })
        .collect();
    if conductance > 0.0 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 0,
                a: CellId(100),
                b: CellId(101),
                area_m2: AREA,
                conductance_w_k: conductance,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn bind_shell<'a>(
    fixture: &'a Fixture,
    mesh: &'a ShellMesh,
    thermal: &DepositionSimulation<'_>,
    count: usize,
    cx: &Cx<'_>,
) -> ThermoplasticShell<'a> {
    let step = THICKNESS / count as f64;
    let layers = (0..count)
        .map(|i| ShellThermoplasticLayer {
            cell: CellId(100 + i as u64),
            layer: PlaneStressLayer {
                bottom_z_m: -THICKNESS / 2.0 + step * i as f64,
                top_z_m: -THICKNESS / 2.0 + step * (i + 1) as f64,
                subdivisions: 2,
                law: ThermoplasticSource::try_new(fixture.query(), 0.0).unwrap(),
            },
        })
        .collect();
    let mut d = discretization();
    // Reference preparation must not consume the free-coordinate solve budget.
    // This membrane has exactly two free coordinates, even though it has 18
    // physical coordinates in its original mass and kinematic operators.
    d.shell_budget.maximum_free_dofs = 2;
    ThermoplasticShell::try_new(thermal, mesh, vec![layers], d, cx).unwrap()
}

fn membrane_supports() -> Vec<(usize, f64)> {
    // Free axial and lateral contraction: epsilon_yy follows sigma_yy=0.
    (0..18)
        .filter(|&d| d != 6 && d != 13)
        .map(|d| (d, 0.0))
        .collect()
}

#[test]
fn g1_loaded_shell_softens_yields_heats_and_retains_permanent_set() {
    context(|cx| {
        let cp = 100.0;
        let dt = 0.001;
        let sigma = 3.4e6;
        let fixture = Fixture::new(cp);
        let mesh = mesh();
        let heat = RHO * VOLUME * cp * 6.0;
        let thermal = thermal(&fixture, &[heat / dt], 0.0);
        let shell = bind_shell(&fixture, &mesh, &thermal, 1, cx);
        let old_solution = shell.solution().clone();
        let old_cells = thermal.cells().clone();
        let loads = [(6, sigma * VOLUME / L)];
        let supports = membrane_supports();
        let step = shell
            .advance(
                &thermal,
                &(),
                dt,
                &supports,
                &loads,
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        // Independent uniaxial J2+caloric oracle. Hardening is stored, and the
        // endpoint integration loss is not deposited in the caloric cell.
        let mut lo = 300.0;
        let mut hi = 319.0;
        for _ in 0..80 {
            let t = 0.5 * (lo + hi);
            let yield_pa = YIELD - SOFTENING * (t - 300.0);
            let p = (sigma - yield_pa) / H;
            let residual = RHO * cp * (t - 300.0) - heat / VOLUME - yield_pa * p;
            if residual > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let expected_t = 0.5 * (lo + hi);
        let yield_pa = YIELD - SOFTENING * (expected_t - 300.0);
        let p = (sigma - yield_pa) / H;
        close(
            step.thermal().temperature_k(CellId(100)).unwrap(),
            expected_t,
            3e-8,
        );
        close(
            step.shell().solution().coordinates[1][0],
            L * (sigma / E + p),
            3e-11,
        );
        close(
            step.shell().solution().coordinates[2][1],
            L * (-0.2 * sigma / E - 0.5 * p),
            3e-11,
        );
        for point in step
            .shell()
            .solution()
            .section_states
            .iter()
            .flat_map(|s| &s.points)
        {
            close(point.history.alpha, p, 2e-9);
            close(point.history.plastic_strain[2], -0.5 * p, 2e-9);
        }
        close(
            step.report().mechanical.physical_heat_j,
            VOLUME * yield_pa * p,
            2e-9,
        );
        assert!(step.report().mechanical.material_numerical_loss_j > 0.0);
        close(step.report().mechanical.reversible_heat_j, 0.0, 1e-12);
        assert!(step.report().combined_residual_j.abs() <= config().energy_tolerance_j);
        assert_eq!(shell.solution(), &old_solution);
        assert_eq!(thermal.cells(), &old_cells);

        // The next thermal interval still has authored external heat. Releasing
        // force must retain the already accumulated plastic set and add no heat
        // merely because elastic strain and numerical work are released.
        let unloaded = step
            .shell()
            .advance(
                step.thermal(),
                &(),
                2.0 * dt,
                &supports,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(unloaded.shell().solution().coordinates[1][0], L * p, 3e-11);
        close(unloaded.report().mechanical.physical_heat_j, 0.0, 1e-12);
        assert!(unloaded.report().combined_residual_j.abs() <= config().energy_tolerance_j);
    });
}

#[test]
fn g3_layer_heat_and_original_conduction_reverse_the_thermal_bending_moment() {
    context(|cx| {
        let cp = 100.0;
        let dt = 0.001;
        let fixture = Fixture::with_expansion(cp, Some(1e-5));
        let mesh = mesh();
        let power = RHO * VOLUME * cp * 3.0 / dt;
        let fixed: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let mut moments = Vec::new();
        for powers in [[power, 0.0], [0.0, power]] {
            let thermal = thermal(&fixture, &powers, 0.002);
            let shell = bind_shell(&fixture, &mesh, &thermal, 2, cx);
            let step = shell
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let t0 = step.thermal().temperature_k(CellId(100)).unwrap();
            let t1 = step.thermal().temperature_k(CellId(101)).unwrap();
            let a = -0.5 * E / (1.0 - 0.2) * 1e-5 * THICKNESS;
            // Both planes are constrained. Native thermal expansion therefore
            // produces sigma_xx=sigma_yy=-E alpha dT/(1-nu), sigma_zz=0.
            let expected_n = a * (t0 + t1 - 600.0);
            let expected_m = a * THICKNESS / 4.0 * (t1 - t0);
            let resultants = step.shell().solution().section_resultants[0];
            close(resultants[0], expected_n, 2e-6);
            close(resultants[3], expected_m, 2e-9);
            moments.push(resultants[3]);
            assert_eq!(step.report().cell_heat_j.len(), 2);
            close(step.report().mechanical.physical_heat_j, 0.0, 1e-12);
            assert!(step.report().mechanical.reversible_heat_j != 0.0);
            assert!(step.report().combined_residual_j.abs() <= config().energy_tolerance_j);
            assert!(t0.min(t1) > 300.0); // Original graph has transferred heat.
        }
        assert!(moments[0] * moments[1] < 0.0);
        close(moments[0], -moments[1], 2e-9);
    });
}

#[derive(Clone)]
struct RefuseAfterHeat;
impl DepositionObserver for RefuseAfterHeat {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Err(DepositionError::InvalidInput("synthetic late refusal"))
    }
}

#[test]
fn g4_late_refusal_and_mismatched_thermal_checkpoint_leave_shell_and_heat_unchanged() {
    context(|cx| {
        let fixture = Fixture::new(100.0);
        let mesh = mesh();
        let thermal = thermal(&fixture, &[10.0], 0.0);
        let shell = bind_shell(&fixture, &mesh, &thermal, 1, cx);
        let old_solution = shell.solution().clone();
        let old_cells = thermal.cells().clone();
        let result = shell.advance(
            &thermal,
            &RefuseAfterHeat,
            0.001,
            &membrane_supports(),
            &[(6, 3.4e6 * VOLUME / L)],
            IncrementSettings::default(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(result, Err(ShellThermoplasticError::Observer(_))));
        let advanced = thermal
            .advance_to(0.001, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert!(
            shell
                .advance(
                    &advanced,
                    &(),
                    0.002,
                    &membrane_supports(),
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                )
                .is_err()
        );
        assert_eq!(shell.solution(), &old_solution);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(shell.time_s(), 0.0);
    });
}
