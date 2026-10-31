//! Work-supplied molecular heating, source feedback and atomic acceptance.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::{MaterialFlowError, ViscousHeatingAudit};
use fs_lbm::conjugate::{
    ChtError, FluidProperties, ViscousHeatControls, creeping_flow_with_viscous_heat,
};

#[path = "viscous/forcing.rs"]
mod forcing;
#[path = "viscous/unsteady.rs"]
mod unsteady;
#[path = "viscous/obstacles.rs"]
mod obstacles;
#[path = "viscous/rheology.rs"]
mod rheology;
#[path = "viscous/porous.rs"]
mod porous;

fn driven_box() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [1.0, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-3;
    config
}

fn controls() -> ViscousHeatControls {
    ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-7,
    }
}

fn thermal_step() -> StepConfig {
    StepConfig {
        // The shared Advected owner further tightens max_dt harmonically.
        // This proposal leaves the explicit 0.5 s advance endpoint in control.
        max_dt_s: 10.0,
        energy_tolerance_j: 1e-8,
    }
}

fn source(scale: f64, thinning: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            100.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.01,
        ),
    ] {
        claims.insert_claim(scalar(name, dims, value)).unwrap();
    }
    let mut mu = scalar(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, scale);
    if thinning {
        mu.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, scale), (310.0, 0.1 * scale), (810.0, 0.01 * scale)],
            dims: DYNAMIC_VISCOSITY_DIMS,
        };
        mu.interpolation = InterpolationPolicy::LinearInside;
    }
    claims.insert_claim(mu).unwrap();
    card(claims)
}

#[test]
fn native_stokes_wall_traction_and_heat_match_an_independent_discrete_solution() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let mut energies = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu, false);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let ids = [
            Some(CellId(0)),
            Some(CellId(1)),
            Some(CellId(2)),
            Some(CellId(3)),
        ];
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_viscous_heating(controls())
            .unwrap();
        let (next, accepted, thermal) = old
            .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .unwrap();
        let epoch = accepted.applied().unwrap();
        let power = epoch.viscous_power.as_ref().unwrap();
        // The four divergence-free free MAC faces have one amplitude a.
        // Direct assembly gives D = mu*(24a²-4Ua+2U²), P=2mu*U*(U-a).
        // The independent Stokes balance gives a=U/12 and P=D=11mu*U²/6.
        // This is the exact discrete box, not a continuum corner-singularity claim.
        near(epoch.flow.field.flux_x(1, 0, 0), -1.0 / 12.0);
        near(epoch.flow.field.flux_x(1, 1, 0), 1.0 / 12.0);
        near(epoch.flow.field.flux_y(0, 1, 0), 1.0 / 12.0);
        near(epoch.flow.field.flux_y(1, 1, 0), -1.0 / 12.0);
        near(power.heat_w, 11.0 * mu / 6.0);
        near(power.wall_tangential_force_n[3][0], 11.0 * mu / 6.0);
        near(power.supplied_wall_power_w, 11.0 * mu / 6.0);
        near(power.wall_power_w[3], power.heat_w);
        near(power.pressure_constraint_power_w, 0.0);
        near(power.residual_w, 0.0);
        for c in 0..4 {
            let expected = if c < 2 { mu / 24.0 } else { 7.0 * mu / 8.0 };
            near(power.cell_heat_w[c], expected);
            // Initial h is uniform, so this first step has no advected heat.
            near(
                next.cells()[&CellId(c as u64)].specific_h_j_kg,
                0.5 * expected,
            );
            assert_eq!(next.cells()[&CellId(c as u64)].boundary.source_w, 0.0);
        }
        let work = accepted.viscous_heating_audit().unwrap();
        assert_eq!(work.substeps, 1);
        near(work.duration_s, 0.5);
        near(work.generated_heat_j, 0.5 * power.heat_w);
        near(work.supplied_wall_work_j, work.generated_heat_j);
        near(next.energy_j() - old.energy_j(), work.supplied_wall_work_j);
        near(thermal.external_heat_j, work.generated_heat_j);
        near(thermal.residual_j, 0.0);
        near(accepted.audit().external_enthalpy_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(
            observer.viscous_heating_audit(),
            Some(ViscousHeatingAudit::default())
        );
        energies.push(next.energy_j());
    }
    near(energies[1], 2.0 * energies[0]);
}

#[test]
fn graded_heterogeneous_stress_includes_reciprocal_shear_and_shared_heat_partition() {
    let x = [1.0, 2.0];
    let y = [0.5, 1.0];
    let depth = 0.75;
    let domain = VoxelDomain::graded([x.to_vec(), y.to_vec(), vec![depth]]).unwrap();
    let mu = [100.0, 200.0, 300.0, 400.0];
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 100.0,
    };
    let (flow, power) = creeping_flow_with_viscous_heat(
        &domain,
        &fluid,
        &mu,
        &driven_box(),
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // On a 2x2 closed MAC grid, exact divergence leaves ONE volumetric flux q.
    // Assemble its scalar dissipation independently: D=Aq²-2Bq+C. For unit
    // top-wall speed, q=B/A and P=C-Bq. Unlike the square box, gamma_ad != 0.
    let lx = 0.5 * (x[0] + x[1]);
    let ly = 0.5 * (y[0] + y[1]);
    let resistance: f64 = (0..4).map(|c| x[c % 2] * y[c / 2] / mu[c]).sum();
    let edge_mu = (x[0] + x[1]) * (y[0] + y[1]) / resistance;
    let gamma_per_q = ((1.0 / y[0] + 1.0 / y[1]) / ly - (1.0 / x[0] + 1.0 / x[1]) / lx) / depth;
    let shear_a = edge_mu * lx * ly * depth * gamma_per_q * gamma_per_q;
    assert!(shear_a > 0.0);
    let normal_a: f64 = (0..4)
        .map(|c| 4.0 * mu[c] / (x[c % 2] * y[c / 2] * depth))
        .sum();
    let wall_x = [
        depth * (x[0] * mu[0] + x[1] * mu[1]) / y[0],
        depth * (x[0] * mu[2] + x[1] * mu[3]) / y[1],
    ];
    let wall_y = [
        depth * (y[0] * mu[0] + y[1] * mu[2]) / x[0],
        depth * (y[0] * mu[1] + y[1] * mu[3]) / x[1],
    ];
    let a = normal_a
        + shear_a
        + (0..2)
            .map(|i| wall_x[i] / (y[i] * depth).powi(2) + wall_y[i] / (x[i] * depth).powi(2))
            .sum::<f64>();
    let b = wall_x[1] / (y[1] * depth);
    let q = b / a;
    near(flow.field.flux_x(1, 0, 0), -q);
    near(flow.field.flux_x(1, 1, 0), q);
    near(flow.field.flux_y(0, 1, 0), q);
    near(flow.field.flux_y(1, 1, 0), -q);
    near(power.supplied_wall_power_w, wall_x[1] - b * q);
    near(power.heat_w, wall_x[1] - b * q);
    let shear_w = shear_a * q * q;
    for c in 0..4 {
        let i = c % 2;
        let j = c / 2;
        let normal = 4.0 * mu[c] * q * q / (x[i] * y[j] * depth);
        let shear = shear_w * (x[i] * y[j] / mu[c]) / resistance;
        let u = if j == 0 {
            -q / (y[j] * depth)
        } else {
            q / (y[j] * depth)
        };
        let wall = if j == 0 { 0.0 } else { 1.0 };
        let along_x = depth * x[i] * mu[c] / y[j] * (wall - u).powi(2);
        let along_y = depth * y[j] * mu[c] / x[i] * (q / (x[i] * depth)).powi(2);
        near(power.cell_heat_w[c], normal + shear + along_x + along_y);
    }
    near(power.residual_w, 0.0);
}

#[test]
fn source_temperature_alters_next_traction_power_and_accepted_enthalpy() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let mut outcomes = Vec::new();
    for thinning in [false, true] {
        let card = source(100.0, thinning);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_viscous_heating(controls())
            .unwrap();
        let (warm, warm_observer, _) = old
            .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .unwrap();
        let first_work = warm_observer.viscous_heating_audit().unwrap();
        assert_eq!(first_work.substeps, 1);
        near(first_work.generated_heat_j, 0.5 * 1100.0 / 6.0);
        let (next, accepted, thermal) = warm
            .advance_observed(&warm_observer, 1.0, thermal_step(), || Ok(()))
            .unwrap();
        assert_eq!(accepted.momentum_solves(), 1);
        let applied = accepted.applied().unwrap();
        near(applied.time_s, 0.5);
        for sample in &applied.samples {
            let first_heat = if sample.voxel < 2 { 100.0 / 24.0 } else { 87.5 };
            let expected_t = 300.0 + 0.5 * first_heat / 100.0;
            near(sample.temperature_k, expected_t);
            let expected_mu = if thinning {
                100.0 * (1.0 - 0.09 * (expected_t - 300.0))
            } else {
                100.0
            };
            near(sample.transport.dynamic_viscosity_pa_s(), expected_mu);
            let receipt = &sample
                .material
                .material()
                .property(DYNAMIC_VISCOSITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt;
            assert_eq!(
                receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
        let power = applied.viscous_power.as_ref().unwrap();
        let work = accepted.viscous_heating_audit().unwrap();
        near(work.generated_heat_j, 0.5 * power.heat_w);
        near(next.energy_j() - warm.energy_j(), work.supplied_wall_work_j);
        near(
            next.energy_j() - old.energy_j(),
            first_work.supplied_wall_work_j + work.supplied_wall_work_j,
        );
        near(thermal.external_heat_j, work.generated_heat_j);
        near(work.residual_j, 0.0);
        near(power.wall_tangential_force_n[3][0], power.wall_power_w[3]);
        for sample in accepted.endpoint() {
            near(
                sample.temperature_k,
                next.temperature_k(sample.cell).unwrap(),
            );
            near(
                sample.specific_enthalpy_j_kg,
                next.cells()[&sample.cell].specific_h_j_kg,
            );
        }
        outcomes.push((
            power.wall_tangential_force_n[3][0],
            power.heat_w,
            next.energy_j(),
        ));
    }
    assert!(outcomes[1].0 < 0.99 * outcomes[0].0);
    assert!(outcomes[1].1 < 0.99 * outcomes[0].1);
    assert!(outcomes[1].2 < outcomes[0].2);
}

#[derive(Clone, Default, Debug, PartialEq)]
struct History {
    steps: usize,
    reject: Option<usize>,
}
impl DepositionObserver for History {
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
        self.steps += 1;
        if self.reject == Some(self.steps) {
            return Err(DepositionError::InvalidInput(
                "independent endpoint refusal",
            ));
        }
        Ok(())
    }
}

#[test]
fn failed_late_candidate_and_retry_do_not_repeat_viscous_heat() {
    let card = source(100.0, true);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let old = initial(&domain, 1.0, &caloric);
    let history = History {
        steps: 0,
        reject: Some(2),
    };
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, history.clone(), &gate, 10)
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    assert!(matches!(
        old.advance_observed(
            &observer,
            1.0,
            StepConfig {
                max_dt_s: 0.5,
                ..thermal_step()
            },
            || Ok(())
        ),
        Err(MaterialFlowError::Inner(DepositionError::InvalidInput(
            "independent endpoint refusal"
        )))
    ));
    assert_eq!(observer.inner(), &history);
    assert_eq!(
        observer.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
    assert!(observer.applied().is_none());
    near(old.energy_j(), 0.0);
    near(old.time_s(), 0.0);
    let (retried, accepted, _) = old
        .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    let clean = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, History::default(), &gate, 10)
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
    let (fresh, fresh_observer, _) = old
        .advance_observed(&clean, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(retried.cells(), fresh.cells());
    assert_eq!(
        accepted.viscous_heating_audit(),
        fresh_observer.viscous_heating_audit()
    );
    assert_eq!(accepted.inner().steps, 1);
    near(retried.energy_j(), 0.5 * 1100.0 / 6.0);
    gate.request();
    assert!(matches!(
        old.advance_observed(&observer, 0.5, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
    ));
    near(old.energy_j(), 0.0);
    assert_eq!(
        observer.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
}

#[test]
fn unsupported_open_driven_or_high_reynolds_regimes_refuse_and_default_adds_no_heat() {
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = driven_box();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let old = initial(&domain, 1.0, &caloric);
    let ordinary = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap();
    let (next, accepted, _) = old
        .advance_observed(&ordinary, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    near(next.energy_j(), 0.0);
    assert!(accepted.viscous_heating_audit().is_none());
    assert!(accepted.applied().unwrap().viscous_power.is_none());
    let open = throughflow();
    let unclosed = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &open, (), &gate, 10)
        .unwrap();
    assert!(unclosed.with_viscous_heating(controls()).is_err());
    let mut unresolved = config.clone();
    unresolved.faces[3] = FvBoundary::Wall {
        velocity: [0.0, 0.0, 1.0],
    };
    let no_dof = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &unresolved, (), &gate, 10)
        .unwrap();
    assert!(no_dof.with_viscous_heating(controls()).is_err());
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 1.0,
    };
    assert!(matches!(
        creeping_flow_with_viscous_heat(&domain, &fluid, &[1.0; 4], &config, controls(), &gate),
        Err(ChtError::InvalidInput {
            field: "simple.viscous_heat",
            ..
        })
    ));
}

#[path = "viscous/pressure_work.rs"]
mod pressure_work;
