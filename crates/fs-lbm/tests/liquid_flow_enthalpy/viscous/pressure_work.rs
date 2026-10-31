//! G1/G3/G4: pressure work drives source-resolved flow, heat and stream enthalpy.
//! Synthetic channel/porous data verify the stated normal-traction and zero
//! tangential-velocity mouth model, not a general unconstrained reservoir mouth.
use super::*;
use fs_lbm::conjugate::{FlowField, creeping_flow_with_pressure_work_heat};
use fs_material::porous::permeability::INTRINSIC_PERMEABILITY_DIMS;
use fs_material::porous::permeability::temperature::PermeabilityTemperatureLaw;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

fn channel(plane: [usize; 3], drop: f64, gauge: f64, walls: bool) -> SimpleConfig {
    let [a, b, _] = plane;
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * a] = FvBoundary::Pressure {
        pressure_pa: gauge + drop,
    };
    faces[2 * a + 1] = FvBoundary::Pressure { pressure_pa: gauge };
    if walls {
        faces[2 * b] = FvBoundary::wall();
        faces[2 * b + 1] = FvBoundary::wall();
    }
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-10;
    config
}

fn flux(field: &FlowField, a: usize, p: [usize; 3]) -> f64 {
    match a {
        0 => field.flux_x(p[0], p[1], p[2]),
        1 => field.flux_y(p[0], p[1], p[2]),
        _ => field.flux_z(p[0], p[1], p[2]),
    }
}

fn fluid(mu: f64) -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: mu,
    }
}

// Independent quadratic in physical face velocities. It has no access to
// SIMPLE matrix entries or the production heat report. Its variation under
// mass-conserving flux perturbations tests the mouth's reciprocal stress rows.
fn channel_dissipation(x: &[f64], y: &[f64], mu: &[f64], u: &[Vec<f64>], v: &[Vec<f64>]) -> f64 {
    let nx = x.len();
    let ny = y.len();
    let m = |i: usize, j: usize| mu[j * nx + i];
    let mut d = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            d += 2.0
                * m(i, j)
                * x[i]
                * y[j]
                * (((u[i + 1][j] - u[i][j]) / x[i]).powi(2)
                    + ((v[i][j + 1] - v[i][j]) / y[j]).powi(2));
        }
    }
    for i in 1..nx {
        for j in 1..ny {
            let lx = 0.5 * (x[i - 1] + x[i]);
            let ly = 0.5 * (y[j - 1] + y[j]);
            let resistance: f64 = [i - 1, i]
                .into_iter()
                .flat_map(|a| [j - 1, j].into_iter().map(move |b| x[a] * y[b] / m(a, b)))
                .sum();
            let edge_mu = 4.0 * lx * ly / resistance;
            let shear = (u[i][j] - u[i][j - 1]) / ly + (v[i][j] - v[i - 1][j]) / lx;
            d += edge_mu * lx * ly * shear * shear;
        }
    }
    // No-slip transverse walls. Each real adjacent cell supplies half its
    // streamwise width; an opening has one real half, never a mirrored volume.
    for i in 0..=nx {
        for j in [0, ny - 1] {
            let weight: f64 = [i.checked_sub(1), (i < nx).then_some(i)]
                .into_iter()
                .flatten()
                .map(|a| x[a] * m(a, j) / y[j])
                .sum();
            d += weight * u[i][j] * u[i][j];
        }
    }
    // At the mouths, the tangential trace is exactly zero. The half-edge
    // shear includes BOTH the normal-velocity variation and the interior
    // tangential velocity's centre-to-mouth derivative.
    for (face, cell, sign) in [(0, 0, 1.0), (nx, nx - 1, -1.0)] {
        for j in 1..ny {
            let half = 0.5 * x[cell];
            let ly = 0.5 * (y[j - 1] + y[j]);
            let edge_mu = (y[j - 1] + y[j]) / (y[j - 1] / m(cell, j - 1) + y[j] / m(cell, j));
            let shear = (u[face][j] - u[face][j - 1]) / ly + sign * v[cell][j] / half;
            d += edge_mu * half * ly * shear * shear;
        }
    }
    d
}

#[test]
fn heterogeneous_graded_channel_satisfies_independent_constrained_work_stationarity() {
    let x = [0.5, 1.0, 1.5];
    let y = [0.75, 1.25];
    let mu = [100.0, 180.0, 75.0, 260.0, 90.0, 210.0];
    let domain = VoxelDomain::graded([x.to_vec(), y.to_vec(), vec![1.0]]).unwrap();
    let config = channel([0, 1, 2], 30.0, 1e5, true);
    let (flow, power) = creeping_flow_with_pressure_work_heat(
        &domain,
        &fluid(100.0),
        &mu,
        &config,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let u: Vec<Vec<_>> = (0..=3)
        .map(|i| (0..2).map(|j| flow.field.flux_x(i, j, 0) / y[j]).collect())
        .collect();
    let v: Vec<Vec<_>> = (0..3)
        .map(|i| (0..=2).map(|j| flow.field.flux_y(i, j, 0) / x[i]).collect())
        .collect();
    assert!(
        v.iter().any(|row| row[1].abs() > 1e-8),
        "fixture must excite reciprocal mouth shear"
    );
    let q: f64 = (0..2).map(|j| flow.field.flux_x(0, j, 0)).sum();
    let d = channel_dissipation(&x, &y, &mu, &u, &v);
    near(d, 30.0 * q);
    near(power.heat_w, d);
    assert!(power.cell_heat_w.iter().all(|heat| *heat >= 0.0));

    let epsilon = 1e-4;
    for through in [false, true] {
        let energies = [-1.0, 1.0].map(|sign| {
            let mut trial_u = u.clone();
            let mut trial_v = v.clone();
            if through {
                // Equal flux perturbation on every x face in one row: div=0
                // and the independently prescribed pressure work changes.
                for row in &mut trial_u {
                    row[0] += sign * epsilon / y[0];
                }
            } else {
                // Four-face circulation across two adjacent cells in x/y.
                // Every affected cell still has exactly zero net flux.
                trial_u[1][0] += sign * epsilon / y[0];
                trial_u[1][1] -= sign * epsilon / y[1];
                trial_v[0][1] -= sign * epsilon / x[0];
                trial_v[1][1] += sign * epsilon / x[1];
            }
            channel_dissipation(&x, &y, &mu, &trial_u, &trial_v)
        });
        // Stationarity of 1/2 D(u)-W(u), restricted to incompressible fields.
        near(
            (energies[1] - energies[0]) / (4.0 * epsilon),
            if through { 30.0 } else { 0.0 },
        );
    }
}

#[test]
fn normal_traction_channel_matches_discrete_poiseuille_in_each_axis_and_gauge() {
    let gate = CancelGate::new();
    for plane in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
        let [a, b, c] = plane;
        let mut n = [0; 3];
        n[a] = 3;
        n[b] = 4;
        n[c] = 1;
        let domain = VoxelDomain::new(n[0], n[1], n[2], 1.0).unwrap();
        for mu in [100.0, 200.0] {
            for (sign, gauge) in [(1.0, 0.0), (-1.0, 1e6)] {
                let config = channel(plane, sign * 30.0, gauge, true);
                let (flow, power) = creeping_flow_with_pressure_work_heat(
                    &domain,
                    &fluid(mu),
                    &vec![mu; domain.cell_count()],
                    &config,
                    controls(),
                    &gate,
                )
                .unwrap();
                // Independent centred FV channel solution, dy=depth=1:
                // u_j=G/(2mu) [y_j(H-y_j)+dy²/4], H=4, L=3.
                // The dy²/4 term is the exact half-cell no-slip wall stencil.
                let gradient = sign * 10.0;
                let mut q = 0.0;
                for j in 0..4 {
                    let y = j as f64 + 0.5;
                    let u = gradient / (2.0 * mu) * (y * (4.0 - y) + 0.25);
                    q += u;
                    for i in 0..=3 {
                        let mut p = [0; 3];
                        p[a] = i;
                        p[b] = j;
                        near(flux(&flow.field, a, p), u);
                    }
                }
                near(
                    q,
                    gradient * 4.0_f64.powi(3) / (12.0 * mu) * (1.0 + 2.0 / 16.0),
                );
                near(power.supplied_pressure_power_w, sign * 30.0 * q);
                near(
                    power.pressure_power_w.iter().sum(),
                    power.supplied_pressure_power_w,
                );
                near(power.molecular_heat_w, sign * 30.0 * q);
                near(power.heat_w, power.supplied_pressure_power_w);
                near(power.cell_heat_w.iter().sum(), power.heat_w);
                near(power.supplied_wall_power_w, 0.0);
                near(power.supplied_body_power_w, 0.0);
                near(power.pressure_constraint_power_w, 0.0);
                near(power.residual_w, 0.0);
                assert!(power.heat_w > 0.0);
            }
        }
    }
}

#[test]
fn original_viscosity_feedback_changes_pressure_flow_and_stream_heat_once() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = channel([0, 1, 2], 30.0, 1e5, true);
    let gate = CancelGate::new();
    let card = source(100.0, true);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..domain.cell_count())
        .map(|c| Some(CellId(c as u64)))
        .collect();
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_pressure_work_heating(controls())
        .unwrap();
    let dt = 0.5;
    let (first, first_owner, first_audit) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    // Two transverse cells have equal u=G/(2mu)=.05 and equal heat .5 W.
    let initial_power = 3.0;
    for c in first.cells().values() {
        near(c.specific_h_j_kg, 0.25);
        assert_eq!(c.boundary.source_w, 0.0);
    }
    near(first_audit.external_heat_j, dt * initial_power);
    near(first_owner.audit().external_enthalpy_j, 0.0);
    near(
        first_owner
            .viscous_heating_audit()
            .unwrap()
            .supplied_pressure_work_j,
        first_audit.external_heat_j,
    );

    let mu = 100.0 * (1.0 - 0.09 * (300.0025 - 300.0));
    let speed = 10.0 / (2.0 * mu);
    let q = 2.0 * speed;
    let mechanical = 30.0 * q;
    // h_in=0; the actual starting outlet h=.25 leaves with mass rate rho*Q.
    let stream = -q * 0.25;
    let (second, owner, audit) = first
        .advance_observed(&first_owner, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let applied = owner.applied().unwrap();
    near(applied.flow.field.flux_x(0, 0, 0), speed);
    near(applied.viscous_power.as_ref().unwrap().heat_w, mechanical);
    assert!(speed > 0.05);
    near(owner.audit().external_enthalpy_j, dt * stream);
    near(
        owner.viscous_heating_audit().unwrap().generated_heat_j,
        dt * mechanical,
    );
    near(
        owner
            .viscous_heating_audit()
            .unwrap()
            .supplied_pressure_work_j,
        dt * mechanical,
    );
    near(
        second.energy_j() - first.energy_j(),
        dt * (mechanical + stream),
    );
    near(audit.external_heat_j, dt * (mechanical + stream));
    near(audit.residual_j, 0.0);
    near(second.mass_kg(), old.mass_kg());
    for sample in &applied.samples {
        near(sample.temperature_k, 300.0025);
        near(sample.transport.dynamic_viscosity_pa_s(), mu);
    }
    let (replay, replay_owner, _) = first
        .advance_observed(&first_owner, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    assert_eq!(replay.cells(), second.cells());
    assert_eq!(
        replay_owner.viscous_heating_audit(),
        owner.viscous_heating_audit()
    );
    assert!(observer.applied().is_none());
}

const K: [&str; 3] = ["pressure_k1", "pressure_k2", "pressure_k3"];
const BETA: [&str; 3] = ["pressure_beta1", "pressure_beta2", "pressure_beta3"];
const BETA_DIMS: Dims = Dims([-1, 0, 0, 0, 0, 0]);

fn porous_card(beta: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (names, dims, values) in [
        (K, INTRINSIC_PERMEABILITY_DIMS, [0.5, 1.0, 2.0]),
        (BETA, BETA_DIMS, [beta; 3]),
    ] {
        for (name, value) in names.into_iter().zip(values) {
            claims.insert_claim(scalar(name, dims, value)).unwrap();
        }
    }
    card(claims)
}

fn porous_law(card: &MaterialCard) -> PermeabilityTemperatureLaw<'_> {
    let keys = K.map(|name| PropertyKey::new(name, INTRINSIC_PERMEABILITY_DIMS));
    let beta = BETA.map(|name| PropertyKey::new(name, BETA_DIMS));
    let requirements: Vec<_> = keys
        .iter()
        .chain(&beta)
        .map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        })
        .collect();
    let resolved = resolve_material_state_point(
        card,
        &point(300.0, 1e5),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PermeabilityTemperatureLaw::try_new(
        card,
        resolved,
        keys,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "T",
    )
    .unwrap()
    .with_forchheimer(beta, 1.0)
    .unwrap()
}

#[test]
fn material_porous_plug_converts_actual_pressure_work_to_linear_and_cubic_heat() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..domain.cell_count())
        .map(|c| Some(CellId(c as u64)))
        .collect();
    for beta in [1000.0, 2000.0] {
        let skeleton = porous_card(beta);
        for sign in [-1.0, 1.0] {
            let config = channel([0, 1, 2], sign * 30.0, 1e5, false);
            let observer = input(&ids, &liquid, &caloric)
                .feedback(&old, &domain, &config, (), &gate, 10)
                .unwrap()
                .with_porous_materials(vec![Some(porous_law(&skeleton)); domain.cell_count()])
                .unwrap()
                .with_pressure_work_heating(controls())
                .unwrap()
                .with_body_force_heating(
                    vec![[sign * 2.0, 0.0, 0.0]; domain.cell_count()],
                    controls(),
                )
                .unwrap();
            let (next, owner, audit) = old
                .advance_observed(&observer, 0.125, thermal_step(), || Ok(()))
                .unwrap();
            // Pressure gradient 10 and explicit body force 2 act together:
            // (mu/K)u + rho beta |u|u = sign*12. Both supply accounts stay distinct.
            let linear = 200.0;
            let u = sign * 24.0 / (linear + (linear * linear + 48.0 * beta).sqrt());
            let epoch = owner.applied().unwrap();
            let power = epoch.viscous_power.as_ref().unwrap();
            for i in 0..=3 {
                for j in 0..2 {
                    near(epoch.flow.field.flux_x(i, j, 0), u);
                }
            }
            near(power.molecular_heat_w, 0.0);
            near(power.darcy_heat_w, 6.0 * linear * u * u);
            near(power.forchheimer_heat_w, 6.0 * beta * u.abs().powi(3));
            near(power.supplied_pressure_power_w, sign * 30.0 * 2.0 * u);
            near(power.supplied_body_power_w, 6.0 * sign * 2.0 * u);
            near(
                power.heat_w,
                power.supplied_pressure_power_w + power.supplied_body_power_w,
            );
            near(power.residual_w, 0.0);
            for c in 0..domain.cell_count() {
                near(
                    power.cell_heat_w[c],
                    linear * u * u + beta * u.abs().powi(3),
                );
                assert_eq!(
                    epoch.porous_samples[c].material.material().card_identity(),
                    skeleton.content_hash()
                );
            }
            near(next.energy_j() - old.energy_j(), 0.125 * power.heat_w);
            near(audit.external_heat_j, 0.125 * power.heat_w);
            near(
                owner
                    .viscous_heating_audit()
                    .unwrap()
                    .supplied_pressure_work_j,
                0.125 * power.supplied_pressure_power_w,
            );
        }
    }
}

fn narrow_source() -> MaterialCard {
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
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 100.0),
    ] {
        let mut claim = scalar(name, dims, value);
        if name == DYNAMIC_VISCOSITY_PROPERTY {
            claim.validity = claim.validity.with("T", 300.0, 300.001);
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

#[test]
fn final_source_exit_and_late_cancel_preserve_open_flow_and_thermal_history() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = channel([0, 1, 2], 30.0, 0.0, true);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..domain.cell_count())
        .map(|c| Some(CellId(c as u64)))
        .collect();
    for narrow in [true, false] {
        let card = if narrow {
            narrow_source()
        } else {
            source(100.0, false)
        };
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_pressure_work_heating(controls())
            .unwrap();
        let saved = old.cells().clone();
        if narrow {
            // .5 W/cell for .5 s makes T=300.0025, outside the original mu claim.
            assert!(
                old.advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
                    .is_err()
            );
        } else {
            let mut calls = 0;
            old.advance_observed(&observer, 0.5, thermal_step(), || {
                calls += 1;
                Ok(())
            })
            .unwrap();
            assert!(calls > 4);
            let mut at = 0;
            assert!(
                old.advance_observed(&observer, 0.5, thermal_step(), || {
                    at += 1;
                    if at == calls - 1 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
            );
        }
        assert_eq!(old.cells(), &saved);
        near(old.time_s(), 0.0);
        assert!(observer.applied().is_none());
        assert_eq!(
            observer.viscous_heating_audit(),
            Some(ViscousHeatingAudit::default())
        );
        let (short, _, _) = old
            .advance_observed(&observer, 1.0 / 1024.0, thermal_step(), || Ok(()))
            .unwrap();
        assert!(short.energy_j() > old.energy_j());
    }
}
