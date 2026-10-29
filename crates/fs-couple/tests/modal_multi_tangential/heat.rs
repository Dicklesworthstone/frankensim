//! Independent shared-capacity and shared-inertia references on original owners.
use super::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::multiple::heat::{
    ModalContactHeatSource, ModalFrictionHeatSolve, ModalTangentialHeatRoute,
    MultiTangentialHeatError, MultiTangentialHeatSolve, MultiTangentialHeatSource,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{DryTangentialState, TangentialStep, TangentialStepInput, TriboError};
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1.0 / 48000.0;
const CAP: [f64; 3] = [1e-6, 2e-6, 3e-6];
const K: [f64; 3] = [0.1, 0.2, 0.15];
const MASS: [f64; 3] = [0.002, 0.006, 0.01];
const SPEED: [f64; 3] = [1.0, 0.25, 0.0];
const MU: [f64; 2] = [0.8, 0.6];
const BETA: [f64; 2] = [0.04, 0.03];
const ROTATION: [[f64; 2]; 2] = [[0.8, -0.6], [0.6, 0.8]];
const LINK: [u64; 2] = [10, 20];
// The middle caloric cell belongs to the body shared by both modal contacts.
const CELLS: [[usize; 2]; 2] = [[1, 0], [1, 2]];
const BODY: [usize; 3] = [1, 0, 2];

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 {
            Ok(300.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        from: f64,
        to: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial = self.temperature_k(from)?;
        self.temperature_k(to)?;
        Ok(Some(((to - from) / initial).ln_1p()))
    }
}
static CALORIC: Caloric = Caloric;

fn thermal(temperatures: [f64; 3]) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: BODY[i] as u64,
                material_id: MaterialId(0),
                mass_kg: CAP[i],
                specific_h_j_kg: temperatures[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
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
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINK[i],
                a: CellId(CELLS[i][0] as u64),
                b: CellId(CELLS[i][1] as u64),
                area_m2: 1.0,
                conductance_w_k: series_conductance_w_k(
                    1.0,
                    [0.5; 2],
                    CELLS[i].map(|cell| K[cell]),
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn routes(thermal: &DepositionSimulation<'_>) -> [FrictionHeatRoute; 2] {
    std::array::from_fn(|i| {
        FrictionHeatRoute::from_link(thermal, LINK[i], CELLS[i].map(|cell| K[cell])).unwrap()
    })
}

fn bindings(
    routes: &[FrictionHeatRoute; 2],
    order: [usize; 2],
) -> [ModalTangentialHeatRoute<'_>; 2] {
    // Deliberately reverse canonical thermal order, independently of mechanics.
    [1, 0].map(|site| ModalTangentialHeatRoute {
        contact_index: order.iter().position(|&i| i == site).unwrap(),
        route: &routes[site],
        temperature_guess_k: 300.0,
    })
}

fn controls() -> MultiTangentialHeatSolve {
    MultiTangentialHeatSolve {
        solve: ModalFrictionHeatSolve {
            maximum_cells: 3,
            maximum_trials: 512,
            maximum_newton_iterations: 24,
            difference_step_j: 1e-8,
            temperature_scale_j_per_k: 1e-6,
            heat: FrictionHeatTolerances {
                cell_j: 1e-12,
                thermal_j: 1e-12,
                interface_j: 1e-12,
                total_j: 1e-12,
            },
            combined_energy_tolerance_j: 2e-12,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-13,
                relative_tolerance: 1e-10,
                linear_restart: 5,
                forcing_minimum: 1e-6,
                ..Default::default()
            },
        },
        maximum_interfaces: 2,
        maximum_unknowns: 5,
        audit_entropy: true,
    }
}

fn source(
    _: &DepositionSimulation<'_>,
    links: &[u64],
    temperatures: &[f64],
) -> Result<MultiTangentialHeatSource<Vec<f64>>, DepositionError> {
    assert_eq!(links, &LINK);
    Ok(MultiTangentialHeatSource {
        interfaces: (0..2)
            .map(|i| ModalContactHeatSource {
                link_id: LINK[i],
                conductivity_w_m_k: CELLS[i].map(|cell| K[cell]),
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperatures.to_vec(),
        },
    })
}

fn input(speed: [f64; 2]) -> PlanarTangentialModalDrive {
    PlanarTangentialModalDrive {
        contact_temperature_k: None,
        ..drive(speed)
    }
}

// Synthetic affine temperature dependence defines a reference, not calibration.
struct SofteningCoulomb {
    site: usize,
    beta: f64,
}
impl TangentialLaw for SofteningCoulomb {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        "test-shared-modal-thermal-coulomb"
    }
    fn initial_state(&self) -> Self::State {
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: MU[self.site],
        }
        .initial_state()
    }
    fn stored_energy_j(&self, _: &Self::State) -> f64 {
        0.0
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let temperature = input
            .contact_temperature_k
            .ok_or(TriboError::InvalidInput {
                field: "temperature",
            })?;
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: MU[self.site] - self.beta * (temperature - 300.0),
        }
        .advance(state, input)
    }
}

fn free_system(
    order: [usize; 2],
    initial: [f64; 3],
    beta: [f64; 2],
) -> MultiTangentialModalSystem<SofteningCoulomb> {
    let models = (0..3)
        .map(|body| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                48000,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO
                    };
                    2
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: initial[body] * MASS[body].sqrt(),
                    },
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: 0.0,
                    },
                ])
                .unwrap();
            model
        })
        .collect();
    let network =
        CoupledModalSystem::new(models, vec![], network_config(), &CancelGate::new()).unwrap();
    let attachment = |body| {
        std::array::from_fn(|axis| ModalAttachment {
            component: body,
            shapes: (0..2)
                .map(|mode| ROTATION[axis][mode] / (MASS[body] * (mode + 1) as f64).sqrt())
                .collect(),
        })
    };
    let contacts = order
        .into_iter()
        .map(|site| ModalTangentialContact {
            left: attachment(0),
            right: Some(attachment(site + 1)),
            law: SofteningCoulomb {
                site,
                beta: beta[site],
            },
            limits: limits(),
        })
        .collect();
    MultiTangentialModalSystem::new(network, contacts, config(), &CancelGate::new()).unwrap()
}

// Eliminate the three finite cells and solve the remaining two physical
// quadratic equations with their analytic 2x2 Jacobian. This does not call the
// native thermal owner, constitutive adapter or Newton/Krylov implementation.
fn thermal_softening_reference() -> ([f64; 2], [f64; 3], [f64; 2], [f64; 2]) {
    let g = K.map(|k| 2.0 * k);
    let middle = CAP[1] + 2.0 * DT * g[1];
    let coupling = DT * DT * g[1] * g[1] / middle;
    let admittance = [
        [
            DT * g[0] * CAP[0] / (CAP[0] + DT * g[0]) + DT * g[1] - coupling,
            -coupling,
        ],
        [
            -coupling,
            DT * g[2] * CAP[2] / (CAP[2] + DT * g[2]) + DT * g[1] - coupling,
        ],
    ];
    let c = std::array::from_fn::<_, 2, _>(|i| {
        std::array::from_fn::<_, 2, _>(|j| {
            0.5 * DT * DT * (MASS[0].recip() + if i == j { MASS[i + 1].recip() } else { 0.0 })
        })
    });
    let distance = [DT * (SPEED[0] - SPEED[1]), DT * (SPEED[0] - SPEED[2])];
    let mut rise = [0.0; 2];
    for _ in 0..20 {
        let mu = std::array::from_fn::<_, 2, _>(|i| MU[i] - BETA[i] * rise[i]);
        let slip =
            std::array::from_fn::<_, 2, _>(|i| distance[i] - c[i][0] * mu[0] - c[i][1] * mu[1]);
        let residual = std::array::from_fn::<_, 2, _>(|i| {
            admittance[i][0] * rise[0] + admittance[i][1] * rise[1] - mu[i] * slip[i]
        });
        if residual.iter().all(|x| x.abs() < 1e-20) {
            break;
        }
        let j = std::array::from_fn::<_, 2, _>(|i| {
            std::array::from_fn::<_, 2, _>(|k| {
                admittance[i][k] - mu[i] * c[i][k] * BETA[k]
                    + if i == k { BETA[i] * slip[i] } else { 0.0 }
            })
        });
        let determinant = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        rise[0] -= (j[1][1] * residual[0] - j[0][1] * residual[1]) / determinant;
        rise[1] -= (j[0][0] * residual[1] - j[1][0] * residual[0]) / determinant;
    }
    let mu = std::array::from_fn(|i| MU[i] - BETA[i] * rise[i]);
    let heat = std::array::from_fn(|i| mu[i] * (distance[i] - c[i][0] * mu[0] - c[i][1] * mu[1]));
    for i in 0..2 {
        close(
            admittance[i][0] * rise[0] + admittance[i][1] * rise[1],
            heat[i],
            1e-19,
        );
    }
    let body = [
        300.0 + DT * g[0] * rise[0] / (CAP[0] + DT * g[0]),
        300.0 + DT * g[1] * (rise[0] + rise[1]) / middle,
        300.0 + DT * g[2] * rise[1] / (CAP[2] + DT * g[2]),
    ];
    (rise.map(|x| 300.0 + x), body, mu, heat)
}

#[test]
fn both_interfaces_share_actual_capacity_and_change_all_three_body_recoils() {
    let mut owner = free_system([0, 1], SPEED, BETA);
    let mut cells = thermal([300.0; 3]);
    let routes = routes(&cells);
    let before = owner.total_energy_j().unwrap();
    let result = owner
        .step_with_friction_heat(
            &mut cells,
            &bindings(&routes, [0, 1]),
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let (temperature, body_temperature, mu, heat) = thermal_softening_reference();
    assert!((temperature[0] - temperature[1]).abs() > 0.1);
    assert!(owner.compliance_m_per_n()[2].abs() > 0.2 * owner.compliance_m_per_n()[0]);
    assert_eq!(result.contact_indices, [0, 1]);
    for i in 0..2 {
        assert!(mu[i] < MU[i] - 0.05);
        close(
            result.heat.interfaces[i].interface_temperature_k,
            temperature[i],
            2e-7,
        );
        close(result.heat.interfaces[i].dissipated_j, heat[i], 2e-12);
        close(result.source[i], temperature[i], 2e-7);
        assert_eq!(
            owner.states()[i].contact_temperature_k(),
            Some(result.source[i])
        );
        for axis in 0..2 {
            close(
                result.mechanical.contacts[i].traction_n[axis],
                -ROTATION[axis][0] * mu[i],
                2e-8,
            );
        }
    }
    let force = [-mu[0] - mu[1], mu[0], mu[1]];
    for body in 0..3 {
        for axis in 0..2 {
            let velocity: f64 = (0..2)
                .map(|mode| {
                    ROTATION[axis][mode]
                        * owner.network().components()[body].states()[mode].velocity_m_sqrt_kg_per_s
                        / (MASS[body] * (mode + 1) as f64).sqrt()
                })
                .sum();
            close(
                velocity,
                ROTATION[axis][0] * (SPEED[body] + DT * force[body] / MASS[body]),
                2e-10,
            );
        }
        close(
            cells.temperature_k(CellId(body as u64)).unwrap(),
            body_temperature[body],
            2e-7,
        );
    }
    let deposited_on_shared =
        result.heat.interfaces[0].body_heat_j[0] + result.heat.interfaces[1].body_heat_j[0];
    close(
        deposited_on_shared,
        CAP[1] * (body_temperature[1] - 300.0),
        2e-12,
    );
    close(cells.energy_j(), heat.iter().sum(), 2e-12);
    close(
        result.heat.dissipated_j,
        result.mechanical.friction_dissipation_j,
        1e-16,
    );
    close(
        owner.total_energy_j().unwrap() + cells.energy_j(),
        before,
        2e-12,
    );
    assert!(result.combined_energy_residual_j.abs() < 2e-12);
    let entropy = result.heat.entropy.unwrap();
    close(
        entropy.friction_generation_j_k,
        heat[0] / temperature[0] + heat[1] / temperature[1],
        2e-14,
    );
    close(
        entropy.thermal_storage_change_j_k,
        (0..3)
            .map(|i| CAP[i] * ((body_temperature[i] - 300.0) / 300.0).ln_1p())
            .sum(),
        2e-14,
    );
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k >= 0.0);
    assert!(!entropy.has_unclassified_sources);
}

#[test]
fn contact_and_route_permutations_preserve_the_same_shared_physical_solution() {
    let mut reference = free_system([0, 1], SPEED, BETA);
    let mut permuted = free_system([1, 0], SPEED, BETA);
    let mut cells = thermal([300.0; 3]);
    let mut other_cells = thermal([300.0; 3]);
    let routes = routes(&cells);
    let original = reference
        .step_with_friction_heat(
            &mut cells,
            &bindings(&routes, [0, 1]),
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let mut other_bindings = bindings(&routes, [1, 0]);
    other_bindings.reverse();
    let other = permuted
        .step_with_friction_heat(
            &mut other_cells,
            &other_bindings,
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(other.contact_indices, [1, 0]);
    assert_eq!(
        other
            .heat
            .interfaces
            .iter()
            .map(|x| x.link_id)
            .collect::<Vec<_>>(),
        LINK
    );
    for i in 0..2 {
        close(
            other.heat.interfaces[i].interface_temperature_k,
            original.heat.interfaces[i].interface_temperature_k,
            2e-7,
        );
        for axis in 0..2 {
            close(
                other.mechanical.contacts[1 - i].traction_n[axis],
                original.mechanical.contacts[i].traction_n[axis],
                2e-8,
            );
            close(
                permuted.states()[1 - i].slip_m()[axis],
                reference.states()[i].slip_m()[axis],
                1e-13,
            );
        }
    }
    for (a, b) in states(&reference).iter().zip(states(&permuted)) {
        close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 1e-13);
        close(
            a.velocity_m_sqrt_kg_per_s,
            b.velocity_m_sqrt_kg_per_s,
            2e-10,
        );
    }
    for i in 0..3 {
        close(
            cells.temperature_k(CellId(i)).unwrap(),
            other_cells.temperature_k(CellId(i)).unwrap(),
            2e-7,
        );
    }
}

#[test]
fn zero_slip_keeps_original_three_cell_conduction_and_entropy() {
    let mut owner = free_system([0, 1], [0.0; 3], [0.0; 2]);
    let initial = [320.0, 295.0, 280.0];
    let mut cells = thermal(initial);
    let routes = routes(&cells);
    let before = cells.energy_j();
    let result = owner
        .step_with_friction_heat(
            &mut cells,
            &bindings(&routes, [0, 1]),
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    // Independent tridiagonal BE system after eliminating massless interfaces.
    let g = [
        DT * series_conductance_w_k(1.0, [0.5; 2], [K[0], K[1]]).unwrap(),
        DT * series_conductance_w_k(1.0, [0.5; 2], [K[1], K[2]]).unwrap(),
    ];
    let left = CAP[0] + g[0];
    let right = CAP[2] + g[1];
    let middle = CAP[1] + g[0] + g[1] - g[0] * g[0] / left - g[1] * g[1] / right;
    let t1 = (CAP[1] * initial[1]
        + g[0] * CAP[0] * initial[0] / left
        + g[1] * CAP[2] * initial[2] / right)
        / middle;
    let expected = [
        (CAP[0] * initial[0] + g[0] * t1) / left,
        t1,
        (CAP[2] * initial[2] + g[1] * t1) / right,
    ];
    for i in 0..3 {
        close(
            cells.temperature_k(CellId(i as u64)).unwrap(),
            expected[i],
            2e-7,
        );
    }
    for i in 0..2 {
        let t = (K[CELLS[i][0]] * expected[CELLS[i][0]] + K[CELLS[i][1]] * expected[CELLS[i][1]])
            / (K[CELLS[i][0]] + K[CELLS[i][1]]);
        close(result.heat.interfaces[i].interface_temperature_k, t, 2e-7);
        assert_eq!(result.mechanical.contacts[i].traction_n, [0.0; 2]);
        assert_eq!(result.mechanical.contacts[i].slip_distance_m, [0.0; 2]);
    }
    assert_eq!(result.heat.dissipated_j, 0.0);
    assert!(result.heat.interfaces[0].body_heat_j[1] < 0.0);
    assert!(result.heat.interfaces[1].body_heat_j[1] > 0.0);
    close(cells.energy_j(), before, 2e-12);
    let entropy = result.heat.entropy.unwrap();
    assert_eq!(entropy.friction_generation_j_k, 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    close(
        entropy.thermal_storage_change_j_k,
        entropy.interface_conduction_generation_j_k
            + entropy.temporal_defect_j_k
            + entropy.interface_residual_entropy_j_k,
        2e-14,
    );
}

#[test]
fn duplicate_missing_or_misidentified_heat_owners_refuse_before_publication() {
    let mut owner = free_system([0, 1], SPEED, BETA);
    let mut cells = thermal([300.0; 3]);
    let routes = routes(&cells);
    let bindings = bindings(&routes, [0, 1]);
    let history = owner.states().to_vec();
    let bodies = states(&owner);
    let caloric = cells.cells().clone();
    let mut duplicate_contact = bindings;
    duplicate_contact[1].contact_index = duplicate_contact[0].contact_index;
    let mut duplicate_route = bindings;
    duplicate_route[1].route = duplicate_route[0].route;
    let mut out_of_range = bindings;
    out_of_range[0].contact_index = 2;
    for invalid in [
        duplicate_contact.as_slice(),
        duplicate_route.as_slice(),
        out_of_range.as_slice(),
        &bindings[..1],
    ] {
        assert!(matches!(
            owner.step_with_friction_heat(
                &mut cells,
                invalid,
                &[0.0; 6],
                &[input([0.0; 2]); 2],
                controls(),
                &[0.0; 3],
                source,
                &CancelGate::new(),
            ),
            Err(MultiTangentialHeatError::Invalid(_))
        ));
    }
    let wrong_source = |thermal: &DepositionSimulation<'_>, links: &[u64], temperatures: &[f64]| {
        let mut material = source(thermal, links, temperatures)?;
        material.interfaces.reverse();
        Ok::<_, DepositionError>(material)
    };
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &bindings,
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            controls(),
            &[0.0; 3],
            wrong_source,
            &CancelGate::new(),
        ),
        Err(MultiTangentialHeatError::Invalid(_))
    ));
    let mut small = controls();
    small.maximum_unknowns = 4;
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &bindings,
            &[0.0; 6],
            &[input([0.0; 2]); 2],
            small,
            &[0.0; 3],
            source,
            &CancelGate::new(),
        ),
        Err(MultiTangentialHeatError::Invalid(_))
    ));
    assert_eq!(owner.states(), history);
    assert_eq!(states(&owner), bodies);
    assert!(owner.last_frame().is_none());
    assert_eq!(owner.network().samples_rendered(), 0);
    assert_eq!(cells.cells(), &caloric);
    assert_eq!(cells.time_s(), 0.0);
}

#[test]
fn all_loaded_vector_histories_survive_late_thermal_refusal_cancel_and_exact_retry() {
    let mut owner = layered(true, 100.0);
    let mut reference = layered(true, 100.0);
    let mut cells = thermal([300.0; 3]);
    let mut reference_cells = thermal([300.0; 3]);
    let routes = routes(&cells);
    let bindings = bindings(&routes, [0, 1]);
    let initial = [input([0.01, 0.02]), input([0.015, -0.005])];
    for (body, heat) in [
        (&mut owner, &mut cells),
        (&mut reference, &mut reference_cells),
    ] {
        let result = body
            .step_with_friction_heat(
                heat,
                &bindings,
                &[0.0; 6],
                &initial,
                controls(),
                &[0.0; 3],
                source,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(result.heat.dissipated_j, 0.0);
        assert_eq!(heat.energy_j(), 0.0);
        assert!(result.mechanical.contact_energy_j > 0.0);
        assert!(
            body.states()
                .iter()
                .all(|state| state.elastic_slip_m().iter().all(|x| x.abs() > 1e-10))
        );
    }
    let history = owner.states().to_vec();
    let bodies = states(&owner);
    let frame = owner.last_frame().cloned();
    let caloric = cells.cells().clone();
    let next = [input([-0.01, 0.005]), input([-0.005, 0.01])];
    let mut exhausted = controls();
    exhausted.solve.maximum_trials = 1;
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &bindings,
            &[0.0; 6],
            &next,
            exhausted,
            &[0.0; 3],
            source,
            &CancelGate::new(),
        ),
        Err(MultiTangentialHeatError::TrialBudget)
    ));
    let calls = Cell::new(0usize);
    let refuse = |thermal: &DepositionSimulation<'_>, links: &[u64], temperatures: &[f64]| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            Err(DepositionError::MaterialDomain)
        } else {
            source(thermal, links, temperatures)
        }
    };
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &bindings,
            &[0.0; 6],
            &next,
            controls(),
            &[0.0; 3],
            refuse,
            &CancelGate::new(),
        ),
        Err(MultiTangentialHeatError::Source(
            DepositionError::MaterialDomain
        ))
    ));
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);
    let cancel = |thermal: &DepositionSimulation<'_>, links: &[u64], temperatures: &[f64]| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            gate.request();
        }
        source(thermal, links, temperatures)
    };
    assert!(
        owner
            .step_with_friction_heat(
                &mut cells,
                &bindings,
                &[0.0; 6],
                &next,
                controls(),
                &[0.0; 3],
                cancel,
                &gate,
            )
            .is_err()
    );
    assert_eq!(owner.states(), history);
    assert_eq!(states(&owner), bodies);
    assert_eq!(owner.last_frame(), frame.as_ref());
    assert_eq!(owner.network().samples_rendered(), 1);
    assert_eq!(cells.cells(), &caloric);
    assert_eq!(cells.time_s(), DT);
    let actual = owner
        .step_with_friction_heat(
            &mut cells,
            &bindings,
            &[0.0; 6],
            &next,
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let expected = reference
        .step_with_friction_heat(
            &mut reference_cells,
            &bindings,
            &[0.0; 6],
            &next,
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual.mechanical, expected.mechanical);
    assert_eq!(owner.states(), reference.states());
    assert_eq!(states(&owner), states(&reference));
    assert_eq!(cells.cells(), reference_cells.cells());
}
