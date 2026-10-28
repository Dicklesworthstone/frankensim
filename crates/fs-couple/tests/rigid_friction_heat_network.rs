//! G1: independent three-body/two-interface mechanics and finite caloric network.
//! Analytical coefficients are fixtures, not claims about a sourced material.
use fs_contact::tangential::nonsmooth::{
    NonsmoothContact, NonsmoothControls, NonsmoothTangentialState,
};
use fs_couple::friction_heat::network::*;
use fs_couple::friction_heat::rigid::network::*;
use fs_couple::friction_heat::rigid::{RigidFrictionHeatControls, RigidFrictionHeatError};
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::*;
use std::collections::BTreeMap;

const DT: f64 = 0.01;
const A: f64 = 0.01;
const AREA: f64 = core::f64::consts::PI * 0.02 * 0.02;
const MASS: [f64; 3] = [2.0, 3.0, 4.0];
const CAP: [f64; 3] = [0.02, 0.03, 0.04];
const K: [f64; 3] = [40.0, 60.0, 100.0];
const STIFFNESS: f64 = 10000.0;
const PATCH: [&str; 2] = ["z-left", "a-right"];

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
}
static CALORIC: Caloric = Caloric;
fn law(i: usize, t: f64) -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        ["fixture/left", "fixture/right"][i],
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: [0.4, 0.35][i] - [0.01, 0.015][i] * (t - 300.0),
            tangential_stiffness_n_per_m: STIFFNESS,
            torsional_stiffness_nm_per_rad: 1.0,
            torsional_capacity_factor: A / 0.0002_f64.sqrt(),
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .unwrap()
}
fn thermal(temperatures: [f64; 3]) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAP[i],
                specific_h_j_kg: temperatures[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, 2.0 * A * (i as f64 - 1.0)],
                    end_m: [0.0, 0.0, 2.0 * A * i as f64],
                    volume_m3: AREA * 2.0 * A,
                    surface_area_m2: 2.0 * AREA + 2.0 * core::f64::consts::PI * 0.02 * 2.0 * A,
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
                id: 10 * (i as u64 + 1),
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(AREA, [A; 2], [K[i], K[i + 1]]).unwrap(),
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
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn mechanics() -> ([RigidBodyState; 3], [MassProperties; 3]) {
    let properties = MASS.map(|m| {
        MassProperties::new(
            m,
            Vec3::ZERO,
            Vec3::new(2.0 * m * A * A, 2.0 * m * A * A, m * A * A),
        )
        .unwrap()
    });
    let bodies = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                Vec3::new(0.0, 0.0, (2.0 * i as f64 - 1.0) * A),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(MASS[i] * [8.0, 3.0, 0.0][i], 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    (bodies, properties)
}
fn controls() -> RigidFrictionHeatControls {
    let newton = NewtonKrylovConfig {
        absolute_tolerance: 1e-11,
        relative_tolerance: 1e-11,
        linear_restart: 6,
        forcing_minimum: 1e-6,
        ..Default::default()
    };
    RigidFrictionHeatControls {
        mechanical: fs_mbd::tangential_friction::RigidFrictionControls {
            maximum_bodies: 3,
            maximum_contacts: 2,
            geometry_tolerance_m: 1e-12,
            velocity_tolerance_mps: 1e-9,
            linear_momentum_tolerance_ns: 1e-10,
            angular_momentum_tolerance_nms: 1e-10,
            energy_tolerance_j: 1e-9,
            contact: NonsmoothControls {
                maximum_evaluations: 4000,
                maximum_newton_iterations: 40,
                maximum_mechanical_dofs: 18,
                difference_step_mps: 1e-5,
                velocity_tolerance_mps: 1e-10,
                cone_tolerance: 1e-10,
                energy_tolerance_j: 1e-9,
                newton,
            },
        },
        maximum_cells: 3,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 0.05,
        heat: FrictionHeatTolerances {
            cell_j: 1e-9,
            thermal_j: 1e-9,
            interface_j: 1e-9,
            total_j: 1e-9,
        },
        combined_energy_tolerance_j: 2e-9,
        newton,
    }
}
fn source(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<RigidFrictionHeatNetworkSource<()>, &'static str> {
    assert_eq!(ids, &[10, 20]);
    Ok(RigidFrictionHeatNetworkSource {
        interfaces: (0..2)
            .map(|i| RigidInterfaceHeatSource {
                link_id: ids[i],
                law: law(i, temperatures[i]),
                conductivity_w_m_k: [K[i], K[i + 1]],
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    })
}
fn run(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut NonsmoothTangentialState,
    thermal: &mut DepositionSimulation<'static>,
    reverse: bool,
    controls: RigidFrictionHeatControls,
    resolve: impl Fn(
        &DepositionSimulation<'static>,
        &[u64],
        &[f64],
    ) -> Result<RigidFrictionHeatNetworkSource<()>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidFrictionNetworkReport<()>, RigidFrictionHeatError<&'static str>> {
    let laws = [law(0, 300.0), law(1, 300.0)];
    let patches = std::array::from_fn::<_, 2, _>(|i| {
        NormalPatchView::new(
            PATCH[i],
            "fixture/normal",
            "fixture/normal-source",
            NormalPatchAuthority::SyntheticFixture,
            10.0,
            0.02,
            0.02,
            0.0002,
        )
        .unwrap()
    });
    let interfaces = std::array::from_fn::<_, 2, _>(|i| {
        PartialSlipInterface::new(
            PATCH[i],
            "retained-two-region",
            "fixture/interface",
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap()
    });
    let routes = std::array::from_fn::<_, 2, _>(|i| {
        FrictionHeatRoute::from_link(thermal, 10 * (i as u64 + 1), [K[i], K[i + 1]]).unwrap()
    });
    let interval = format!("interval/{}", history.version());
    let contacts = std::array::from_fn::<_, 2, _>(|i| NonsmoothContact {
        law: &laws[i],
        patch: &patches[i],
        interface: &interfaces[i],
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        common_spin_rad_per_s: 0.0,
        torsional_length_m: A,
        creepage_reference_speed_mps: 1.0,
        work: GeneralizedWorkOwnership::new(PATCH[i], &interval, "x", "y", "spin").unwrap(),
    });
    let order = if reverse { [1, 0] } else { [0, 1] };
    let inputs = order.map(|i| HeatedRigidContact {
        contact: &contacts[i],
        route: &routes[i],
        temperature_guess_k: 300.0,
        binding: RigidContactBinding {
            body_a: i,
            body_b: i + 1,
            arm_a_body: Vec3::new(0.0, 0.0, A),
            arm_b_body: Vec3::new(0.0, 0.0, -A),
        },
    });
    apply_rigid_friction_heat_network(
        bodies, properties, history, &inputs, thermal, DT, controls, &[0.0; 3], resolve, gate,
    )
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}");
}

// Formula-disjoint reference: eliminate all THREE caloric temperatures,
// solve only two interface temperatures by an analytical 2x2 Newton system.
// Native shared-body coupling W01=-1/(2*m_middle) is retained explicitly.
fn reference() -> ([f64; 2], [f64; 3], [f64; 2], [f64; 2]) {
    let d = K.map(|k| DT * k * AREA / A);
    let cross = -d[1] * d[1] / (CAP[1] + 2.0 * d[1]);
    let matrix = [
        [d[0] * CAP[0] / (CAP[0] + d[0]) + d[1] + cross, cross],
        [cross, d[2] * CAP[2] / (CAP[2] + d[2]) + d[1] + cross],
    ];
    let mobility = [[1.25, -1.0 / 6.0], [-1.0 / 6.0, 0.875]];
    let regional = 0.3_f64.powi(2) / 0.6 + 0.7_f64.powi(2) / 0.4;
    let slopes = [0.1, 0.15];
    let mut theta = [0.0; 2];
    let mut heat = [0.0; 2];
    for _ in 0..12 {
        let c = [4.0 - slopes[0] * theta[0], 3.5 - slopes[1] * theta[1]];
        let mut residual = [0.0; 2];
        let mut jac = matrix;
        for i in 0..2 {
            let j = 1 - i;
            let self_coefficient = DT * DT * mobility[i][i] + regional / STIFFNESS;
            heat[i] = c[i] * DT * [5.0, 3.0][i]
                - self_coefficient * c[i] * c[i]
                - DT * DT * mobility[i][j] * c[i] * c[j];
            residual[i] = matrix[i][0] * theta[0] + matrix[i][1] * theta[1] - heat[i];
            jac[i][i] += slopes[i]
                * (DT * [5.0, 3.0][i]
                    - 2.0 * self_coefficient * c[i]
                    - DT * DT * mobility[i][j] * c[j]);
            jac[i][j] -= DT * DT * mobility[i][j] * c[i] * slopes[j];
        }
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        theta[0] -= (jac[1][1] * residual[0] - jac[0][1] * residual[1]) / det;
        theta[1] -= (jac[0][0] * residual[1] - jac[1][0] * residual[0]) / det;
    }
    let cell = [
        d[0] * theta[0] / (CAP[0] + d[0]),
        d[1] * (theta[0] + theta[1]) / (CAP[1] + 2.0 * d[1]),
        d[2] * theta[1] / (CAP[2] + d[2]),
    ];
    (
        theta.map(|x| 300.0 + x),
        cell.map(|x| 300.0 + x),
        heat,
        [4.0 - slopes[0] * theta[0], 3.5 - slopes[1] * theta[1]],
    )
}

#[test]
fn g1_shared_rigid_body_and_finite_cell_match_independent_joint_reference() {
    let (mut bodies, properties) = mechanics();
    let before = bodies;
    let mut network = thermal([300.0; 3]);
    let mut history = NonsmoothTangentialState::new(2, 16).unwrap();
    let report = run(
        &mut bodies,
        &properties,
        &mut history,
        &mut network,
        false,
        controls(),
        source,
        &CancelGate::default(),
    )
    .unwrap();
    let (interface, cell, heat, force) = reference();
    assert_eq!(report.mechanical.contacts[0].patch_id, "a-right");
    for i in 0..2 {
        close(
            report.heat.interfaces[i].interface_temperature_k,
            interface[i],
            2e-7,
        );
        close(report.heat.interfaces[i].dissipated_j, heat[i], 1e-9);
        close(report.heat.interfaces[i].interface_residual_j, 0.0, 1e-9);
        close(
            report.interfaces[i].law.parameters().kinetic_mu,
            force[i] / 10.0,
            2e-9,
        );
        let mechanical = &report.mechanical.contacts[1 - i];
        assert_eq!(mechanical.law.state, PartialSlipStateKind::GrossSlide);
        close(mechanical.impulse_ns[0], -DT * force[i], 1e-10);
    }
    for i in 0..3 {
        close(
            network.temperature_k(CellId(i as u64)).unwrap(),
            cell[i],
            2e-7,
        );
    }
    close(
        bodies[0].linear_momentum_world().x,
        before[0].linear_momentum_world().x - DT * force[0],
        1e-10,
    );
    close(
        bodies[1].linear_momentum_world().x,
        before[1].linear_momentum_world().x + DT * (force[0] - force[1]),
        1e-10,
    );
    close(
        bodies[2].linear_momentum_world().x,
        before[2].linear_momentum_world().x + DT * force[1],
        1e-10,
    );
    close(
        bodies[1].angular_momentum_body().y,
        -A * DT * (force[0] + force[1]),
        1e-12,
    );
    close(
        report.heat.thermal.stored_energy_change_j,
        heat[0] + heat[1],
        1e-9,
    );
    close(
        CAP[1] * (cell[1] - 300.0),
        report.heat.interfaces[0].body_heat_j[1] + report.heat.interfaces[1].body_heat_j[0],
        1e-9,
    );
    assert!(
        report.mechanical.kinetic_before_j - report.mechanical.kinetic_after_j > heat[0] + heat[1]
    );
    close(report.combined_energy_residual_j, 0.0, 2e-9);
    assert_eq!(history.time_s(), network.time_s());
}

#[test]
fn g5_input_permutation_and_retained_second_step_are_identical() {
    let (mut a, properties) = mechanics();
    let mut b = a;
    let mut ta = thermal([300.0; 3]);
    let mut tb = ta.clone();
    let mut ha = NonsmoothTangentialState::new(2, 16).unwrap();
    let mut hb = ha.clone();
    for _ in 0..2 {
        let ra = run(
            &mut a,
            &properties,
            &mut ha,
            &mut ta,
            false,
            controls(),
            source,
            &CancelGate::default(),
        )
        .unwrap();
        let rb = run(
            &mut b,
            &properties,
            &mut hb,
            &mut tb,
            true,
            controls(),
            source,
            &CancelGate::default(),
        )
        .unwrap();
        assert_eq!(a, b);
        assert_eq!(ha, hb);
        assert_eq!(ra.mechanical.contacts, rb.mechanical.contacts);
        for i in 0..3 {
            assert_eq!(
                ta.temperature_k(CellId(i)).unwrap(),
                tb.temperature_k(CellId(i)).unwrap()
            );
        }
        assert!(ra.heat.dissipated_j > 0.0);
    }
    assert_eq!(ha.version(), 2);
}

#[test]
fn g0_late_source_budget_and_competing_owner_refusals_publish_nothing() {
    let (original, properties) = mechanics();
    let accepted = thermal([300.0; 3]);
    let original_history = NonsmoothTangentialState::new(2, 16).unwrap();
    for refusal in 0..5 {
        let mut bodies = original;
        let mut network = accepted.clone();
        let mut history = original_history.clone();
        let mut control = controls();
        if refusal == 0 {
            control.maximum_trials = 1;
        }
        if refusal == 1 {
            control.mechanical.maximum_contacts = 1;
        }
        let gate = CancelGate::default();
        let result = run(
            &mut bodies,
            &properties,
            &mut history,
            &mut network,
            false,
            control,
            |state, ids, t| {
                let mut input = source(state, ids, t)?;
                if refusal == 2 {
                    input.transport.conductances_w_k.insert(20, 1.0);
                }
                if refusal == 3 && t[0] > 300.001 {
                    return Err("source-domain");
                }
                if refusal == 4 {
                    gate.request();
                }
                Ok(input)
            },
            &gate,
        );
        assert!(result.is_err());
        assert_eq!(bodies, original);
        assert_eq!(history, original_history);
        assert_eq!(network.time_s(), accepted.time_s());
        for i in 0..3 {
            assert_eq!(network.temperature_k(CellId(i)).unwrap(), 300.0);
        }
        let retry = run(
            &mut bodies,
            &properties,
            &mut history,
            &mut network,
            false,
            controls(),
            source,
            &CancelGate::default(),
        )
        .unwrap();
        close(
            retry.heat.interfaces[0].interface_temperature_k,
            reference().0[0],
            2e-7,
        );
    }
}

#[test]
fn g1_zero_friction_shared_cell_conduction_and_duplicate_route_refusal() {
    let accepted = thermal([330.0, 300.0, 280.0]);
    let routes = [
        FrictionHeatRoute::from_link(&accepted, 10, [K[0], K[1]]).unwrap(),
        FrictionHeatRoute::from_link(&accepted, 20, [K[1], K[2]]).unwrap(),
    ];
    assert!(
        FrictionHeatNetworkStep::new(&accepted, DT, &[&routes[0], &routes[0]], 2, || Ok(()))
            .is_err()
    );
    let step = FrictionHeatNetworkStep::new(&accepted, DT, &[&routes[1], &routes[0]], 2, || Ok(()))
        .unwrap();
    // Independent tridiagonal backward Euler with the original series links.
    let g = [
        DT * series_conductance_w_k(AREA, [A; 2], [K[0], K[1]]).unwrap(),
        DT * series_conductance_w_k(AREA, [A; 2], [K[1], K[2]]).unwrap(),
    ];
    let left = CAP[0] + g[0];
    let right = CAP[2] + g[1];
    let middle = CAP[1] + g[0] + g[1] - g[0] * g[0] / left - g[1] * g[1] / right;
    let t1 =
        (CAP[1] * 300.0 + g[0] * CAP[0] * 330.0 / left + g[1] * CAP[2] * 280.0 / right) / middle;
    let t = [
        (CAP[0] * 330.0 + g[0] * t1) / left,
        t1,
        (CAP[2] * 280.0 + g[1] * t1) / right,
    ];
    let temperatures = [
        (K[0] * t[0] + K[1] * t[1]) / (K[0] + K[1]),
        (K[1] * t[1] + K[2] * t[2]) / (K[1] + K[2]),
    ];
    let energies = std::array::from_fn::<_, 3, _>(|i| CAP[i] * (t[i] - [330.0, 300.0, 280.0][i]));
    let trial = step
        .evaluate(
            &energies,
            &temperatures,
            |_, _| {
                Ok::<_, DepositionError>(FrictionHeatNetworkInput {
                    interfaces: (0..2)
                        .map(|i| InterfaceHeatInput {
                            dissipated_j: 0.0,
                            conductivity_w_m_k: [K[i], K[i + 1]],
                        })
                        .collect(),
                    transport: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    },
                })
            },
            || Ok(()),
        )
        .unwrap();
    let (_, _, report) = trial.accept(controls().heat, || Ok(())).unwrap();
    assert!(report.interfaces[0].body_heat_j[0] < 0.0);
    assert!(report.interfaces[1].body_heat_j[1] > 0.0);
    close(report.thermal.stored_energy_change_j, 0.0, 1e-10);
    assert_eq!(report.dissipated_j, 0.0);
}
