//! G1 finite-body heat partition and G4 atomic sustained shell contact.
//! Material and friction coefficients are authored test fixtures.

use super::*;
use fs_contact::sustained::NaturalMapScales;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_couple::shell_plastic_enthalpy::contact::{
    ShellContactHeatBinding, ShellContactHeatControls, ShellContactHeatError,
    ShellContactHeatSource, ShellContactInterfaceSource,
};
use fs_solid::plastic_shell::contact::{
    ShellContactControls, ShellContactStep, ShellContactWitness,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;

const INTERFACE: u64 = 41;
const CONTACT: u64 = 73;
const OBSTACLE: CellId = CellId(200);
const CONDUCTIVITY: f64 = 20.0;
const CP_CONTACT: f64 = 100.0;
const DT_CONTACT: f64 = 0.001;
const STRETCH: f64 = 0.008;
const INDENT: f64 = 2e-5;
const SLIDE: f64 = 1e-4;

fn thermal_pair(fixture: &Fixture) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    // Original orthogonal prism centrelines run between opposite face centres.
    // The top skin is the actual shared thermal face, not an arbitrary point.
    for (i, id) in [CellId(100), OBSTACLE].into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                material_id: MaterialId(0),
                segment_id: i as u64,
                body_id: 7 + i as u64,
                mass_kg: RHO * VOLUME,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [L / 3.0, L / 3.0, (i as f64 - 0.5) * THICKNESS],
                    end_m: [L / 3.0, L / 3.0, (i as f64 + 0.5) * THICKNESS],
                    volume_m3: VOLUME,
                    surface_area_m2: 2.0 * AREA + (2.0 + 2.0_f64.sqrt()) * L * THICKNESS,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: INTERFACE,
            a: CellId(100),
            b: OBSTACLE,
            area_m2: AREA,
            conductance_w_k: CONDUCTIVITY * AREA / THICKNESS,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn witness() -> ShellContactWitness {
    ShellContactWitness {
        id: CONTACT,
        triangle: 0,
        barycentric: [1.0 / 3.0; 3],
        offset_m: [0.0, 0.0, THICKNESS / 2.0],
        reference_obstacle_point_m: [L / 3.0, L / 3.0, THICKNESS / 2.0],
        axes_world: [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        area_m2: AREA,
    }
}

fn supports() -> Vec<(usize, f64)> {
    // A single physical transverse coordinate responds to the indenter. All
    // membrane coordinates are prescribed, giving a controlled plastic stretch.
    (0..18)
        .filter(|&d| d != 14)
        .map(|d| (d, if d == 6 { STRETCH * L } else { 0.0 }))
        .collect()
}

fn interval() -> ShellContactStep {
    ShellContactStep {
        obstacle_start_displacement_m: [0.0; 3],
        obstacle_end_displacement_m: [SLIDE, 0.0, -INDENT],
        // The coupled source replaces this value at the same thermal trial.
        friction_coefficient: 0.55,
        initial_force_n: [0.0; 3],
        scales: NaturalMapScales {
            normal_n_per_m: 1e5,
            tangent_n_per_m: 1e5,
        },
    }
}

fn source(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<ShellContactHeatSource<f64>, &'static str> {
    assert_eq!(ids, &[INTERFACE]);
    assert_eq!(temperatures.len(), 1);
    let mu = 0.55 - 0.02 * (temperatures[0] - 300.0);
    Ok(ShellContactHeatSource {
        interfaces: vec![ShellContactInterfaceSource {
            link_id: INTERFACE,
            friction_coefficient: mu,
            conductivity_w_m_k: [CONDUCTIVITY; 2],
        }],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: mu,
        },
    })
}

fn controls() -> ShellContactHeatControls {
    let mut mechanical = ShellContactControls::default();
    mechanical.solve.maximum_contacts = 1;
    mechanical.work_tolerance_j = 1e-8;
    ShellContactHeatControls {
        mechanical,
        maximum_cells: 2,
        maximum_unknowns: 3,
        maximum_trials: 1024,
        maximum_newton_iterations: 32,
        difference_step_j: 1e-7,
        temperature_scale_j_per_k: RHO * VOLUME * CP_CONTACT,
        heat: FrictionHeatTolerances {
            cell_j: 1e-8,
            thermal_j: 1e-8,
            interface_j: 1e-8,
            total_j: 1e-8,
        },
        material_energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        combined_energy_tolerance_j: 2e-8,
        geometry_tolerance_m: 1e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-10,
            linear_restart: 3,
            forcing_minimum: 1e-5,
            ..Default::default()
        },
        audit_entropy: false,
    }
}

#[test]
fn g1_yielding_shell_and_sliding_obstacle_share_only_physical_heat() {
    context(|cx| {
        let fixture = Fixture::new(CP_CONTACT);
        let mesh = mesh();
        let thermal = thermal_pair(&fixture);
        let shell = bind_shell(&fixture, &mesh, &thermal, 1, cx);
        let accepted = shell.solution().clone();
        let accepted_cells = thermal.cells().clone();
        let contacts = shell.prepare_contacts(&[witness()], 1, cx).unwrap();
        let route = FrictionHeatRoute::from_link(&thermal, INTERFACE, [CONDUCTIVITY; 2]).unwrap();
        let bindings = [ShellContactHeatBinding {
            contact_id: CONTACT,
            route: &route,
            shell_side: 0,
            temperature_guess_k: 300.0,
        }];
        let step = shell
            .advance_contact_heat(
                &thermal,
                DT_CONTACT,
                &supports(),
                &[],
                &contacts,
                &[interval()],
                &bindings,
                None,
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                source,
                cx,
            )
            .unwrap();
        assert!(step.motion.is_none());
        assert_eq!(step.contacts.len(), 1);
        let contact = &step.contacts[0];
        assert_eq!(contact.id, CONTACT);
        assert!(contact.force_n[0] > 0.0);
        close(contact.initial_gap_m, 0.0, 1e-14);
        assert!(contact.final_gap_m >= 0.0);
        close(contact.final_gap_m, 0.0, 1e-10);
        close(
            step.shell.solution().coordinates[2][2],
            -3.0 * INDENT,
            3e-10,
        );
        let slip = STRETCH * L / 3.0 - SLIDE;
        close(contact.slip_increment_m[0], slip, 1e-12);
        close(contact.slip_increment_m[1], 0.0, 1e-12);
        close(contact.force_n[1], step.source * contact.force_n[0], 1e-7);
        close(contact.normal_pressure_pa, contact.force_n[0] / AREA, 1e-6);
        let friction_j = -step.source * contact.force_n[0] * slip;
        close(contact.friction_loss_j, friction_j, 1e-9);
        close(contact.normal_numerical_loss_j, 0.0, 1e-8);
        close(
            contact.obstacle_work_j - contact.shell_work_j,
            friction_j,
            1e-8,
        );

        let bulk_j = step.mechanical.physical_heat_j + step.mechanical.reversible_heat_j;
        assert!(step.mechanical.physical_heat_j > 0.0);
        assert!(step.mechanical.material_numerical_loss_j > 1e-8);
        close(step.mechanical.reversible_heat_j, 0.0, 1e-12);
        assert!(
            step.shell
                .solution()
                .section_states
                .iter()
                .flat_map(|s| &s.points)
                .any(|point| point.history.alpha > 0.0)
        );
        let shell_t = step.thermal.temperature_k(CellId(100)).unwrap();
        let obstacle_t = step.thermal.temperature_k(OBSTACLE).unwrap();
        assert!(shell_t > obstacle_t && obstacle_t > 300.0);
        let capacity = RHO * VOLUME * CP_CONTACT;
        let half_conductance = 2.0 * CONDUCTIVITY * AREA / THICKNESS;
        // Independent two-cell backward-Euler equations. Equal side
        // conductances do not imply equal heat shares when only the shell also
        // dissipates plastically. Their difference includes real conduction.
        close(
            capacity * (shell_t + obstacle_t - 600.0),
            bulk_j + friction_j,
            2e-8,
        );
        close(
            (capacity + half_conductance * DT_CONTACT) * (shell_t - obstacle_t),
            bulk_j,
            2e-8,
        );
        let face = &step.heat.interfaces[0];
        let face_t =
            (shell_t + obstacle_t) / 2.0 + friction_j / (2.0 * half_conductance * DT_CONTACT);
        close(face.interface_temperature_k, face_t, 3e-6);
        close(step.source, 0.55 - 0.02 * (face_t - 300.0), 1e-7);
        assert!(step.source < 0.55);
        close(step.heat.dissipated_j, friction_j, 1e-9);
        close(
            step.heat.thermal.stored_energy_change_j,
            bulk_j + friction_j,
            2e-8,
        );
        close(step.external_heat_j, 0.0, 1e-12);
        assert!(step.combined_residual_j.abs() <= controls().combined_energy_tolerance_j);
        assert_eq!(shell.solution(), &accepted);
        assert_eq!(thermal.cells(), &accepted_cells);
    });
}

#[test]
fn g1_dynamic_contact_retains_consistent_mass_and_excludes_inertial_loss_from_heat() {
    context(|cx| {
        let dt = 1e-5;
        let fixture = Fixture::new(CP_CONTACT);
        let mesh = mesh();
        let thermal = thermal_pair(&fixture);
        let shell = bind_shell(&fixture, &mesh, &thermal, 1, cx);
        let mass = shell.prepare_reference_inertia(cx).unwrap();
        let contacts = shell.prepare_contacts(&[witness()], 1, cx).unwrap();
        let route = FrictionHeatRoute::from_link(&thermal, INTERFACE, [CONDUCTIVITY; 2]).unwrap();
        let velocity = vec![[0.0; 6]; mesh.nodes.len()];
        let step = shell
            .advance_contact_heat(
                &thermal,
                dt,
                &supports(),
                &[],
                &contacts,
                &[interval()],
                &[ShellContactHeatBinding {
                    contact_id: CONTACT,
                    route: &route,
                    shell_side: 0,
                    temperature_guess_k: 300.0,
                }],
                Some((&velocity, &mass)),
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                source,
                cx,
            )
            .unwrap();
        let motion = step.motion.as_ref().unwrap();
        // Integral N_i^2 dA=A/6; the two nonzero nodal translations are
        // orthogonal, so their consistent-mass kinetic energies add directly.
        let kinetic_j =
            RHO * VOLUME / 12.0 * ((STRETCH * L / dt).powi(2) + (3.0 * INDENT / dt).powi(2));
        close(motion.velocity[1][0], STRETCH * L / dt, 1e-10);
        close(motion.velocity[2][2], -3.0 * INDENT / dt, 3e-7);
        close(motion.inertia.physical.kinetic_after_j, kinetic_j, 1e-10);
        close(motion.inertia.physical.numerical_loss_j, kinetic_j, 1e-10);
        close(motion.inertia.drilling_gauge.kinetic_after_j, 0.0, 1e-14);
        close(motion.inertia.total.inertial_work_j, 2.0 * kinetic_j, 2e-10);
        assert!(kinetic_j > 1e-5);
        assert!(step.mechanical.physical_heat_j > 0.0);
        assert!(step.contacts[0].friction_loss_j > 0.0);
        let physical_heat_j = step.mechanical.physical_heat_j
            + step.mechanical.reversible_heat_j
            + step.contacts[0].friction_loss_j;
        close(
            step.heat.thermal.stored_energy_change_j,
            physical_heat_j,
            2e-8,
        );
        close(step.contacts[0].normal_numerical_loss_j, 0.0, 1e-8);
        assert!(step.combined_residual_j.abs() <= controls().combined_energy_tolerance_j);
    });
}

#[test]
fn g4_late_source_refusal_does_not_publish_shell_contact_or_cell_history() {
    context(|cx| {
        let fixture = Fixture::new(CP_CONTACT);
        let mesh = mesh();
        let thermal = thermal_pair(&fixture);
        let shell = bind_shell(&fixture, &mesh, &thermal, 1, cx);
        let accepted = shell.solution().clone();
        let accepted_cells = thermal.cells().clone();
        let contacts = shell.prepare_contacts(&[witness()], 1, cx).unwrap();
        let route = FrictionHeatRoute::from_link(&thermal, INTERFACE, [CONDUCTIVITY; 2]).unwrap();
        let calls = std::cell::Cell::new(0usize);
        let result = shell.advance_contact_heat(
            &thermal,
            DT_CONTACT,
            &supports(),
            &[],
            &contacts,
            &[interval()],
            &[ShellContactHeatBinding {
                contact_id: CONTACT,
                route: &route,
                shell_side: 0,
                temperature_guess_k: 300.0,
            }],
            None,
            IncrementSettings::default(),
            controls(),
            &[0.0; 2],
            |state, ids, temperatures| {
                calls.set(calls.get() + 1);
                if calls.get() >= 2 {
                    Err("synthetic late source refusal")
                } else {
                    source(state, ids, temperatures)
                }
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(ShellContactHeatError::Source(
                "synthetic late source refusal"
            ))
        ));
        assert!(calls.get() >= 2);
        assert_eq!(shell.solution(), &accepted);
        assert_eq!(thermal.cells(), &accepted_cells);
        assert_eq!(shell.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
