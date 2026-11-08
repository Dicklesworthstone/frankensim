//! G1/G4 analytic thermal expansion, primary/secondary pyroelectricity and
//! electrocaloric exchange. Constant synthetic SI coefficients make the whole
//! finite-body energy balance a quadratic in temperature, independently of the
//! constitutive derivative and Newton implementations.
use super::*;
use fs_couple::piezoelectric::spatial::transient::acoustic::{
    PlaneWavePortControls, TransientAcousticPort,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

const REFERENCE_K: f64 = 300.0;
const ALPHA: f64 = 1e-5;
const PYRO: f64 = 1e-4;
const CURVATURE: f64 = STIFFNESS[2] * ALPHA * ALPHA + PYRO * PYRO / EPSILON[2];

fn scalar(name: &str, dims: fs_qty::Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            absolute_temperature(),
            280.0,
            400.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored constant pyroelectric/thermal-expansion fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn source() -> (MaterialCard, ResolvedHeatCapacityCurve, TensorPiezoelectric) {
    let original = material([0.0, 0.0, E33, 0.0, 0.0, 0.0]);
    let mut claims = ClaimSet::new();
    for property in original.material().properties() {
        let requirement = property.requirement();
        claims
            .insert_claim(scalar(
                requirement.name(),
                requirement.dims(),
                property.value_si(),
            ))
            .unwrap();
    }
    for requirement in ThermalTensorPiezoelectric::thermal_bias_requirements().unwrap() {
        let value = match requirement.name() {
            "piezoelectric_thermal_expansion_mandel_3" => ALPHA,
            "piezoelectric_pyroelectric_clamped_3" => PYRO,
            _ => 0.0,
        };
        claims
            .insert_claim(scalar(requirement.name(), requirement.dims(), value))
            .unwrap();
    }
    // This caloric chart is the stress-free, zero-field reference. Clamping and
    // charge constraints contribute their own field heat capacity, not a second
    // independently specified clamped heat capacity.
    claims
        .insert_claim(scalar(
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            CP,
        ))
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermally biased reciprocal ceramic".into(),
            phase: "solid".into(),
            process: "stress-free zero-field caloric reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    let tensor = TensorPiezoelectric::resolve(
        &card,
        &point(REFERENCE_K),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    (card, curve, tensor)
}

fn biased_law(card: &MaterialCard, tensor: &TensorPiezoelectric) -> ThermalTensorPiezoelectric {
    law(card, tensor)
        .with_thermal_bias(REFERENCE_K, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

// The physically continuous root of a*dT - b*dT²/2 = q. Rationalization
// avoids subtracting nearly equal capacities for a small heat input.
fn temperature_increment(a: f64, q: f64) -> f64 {
    2.0 * q / (a + (a * a - 2.0 * CURVATURE * q).sqrt())
}

// At S=0 the independent energy U=psi+T*entropy simplifies to this expression;
// p*Tref*D/epsilon is the reversible electrocaloric term, not electrical loss.
fn clamped_internal_energy(temperature: f64, displacement: f64) -> f64 {
    VOLUME
        * (displacement * displacement / (2.0 * EPSILON[2])
            + PYRO * REFERENCE_K * displacement / EPSILON[2]
            - 0.5 * CURVATURE * (temperature * temperature - REFERENCE_K * REFERENCE_K))
}

#[test]
fn g1_clamped_open_circuit_heating_creates_pyroelectric_voltage_and_thermal_reaction() {
    with_cx(|cx, _| {
        for refinement in [1, 2] {
            let geometry = cube(refinement);
            let (card, curve, tensor) = source();
            let materials = [tensor.clone()];
            let model = Model::assemble(
                &geometry.nodes,
                &geometry.tets,
                &materials,
                FRAME,
                budget(),
                cx,
            )
            .unwrap();
            let n = geometry.nodes.len();
            let zero = vec![[0.0; 3]; n];
            let charges = vec![0.0; n];
            let supports = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
            let ordinary = Body::new(
                &model,
                &supports,
                &electrodes(&geometry, 0.0, Condition::Charge(0.0)),
                &zero,
                &zero,
                &charges,
                field_controls(),
                cx,
            )
            .unwrap();
            let cards = [card];
            let curves = [curve];
            let context = point(REFERENCE_K);
            let bindings = cells(&cards, &curves, &context);
            let network = MaterialNetwork {
                cells: &bindings,
                resistors: &[],
                links: &[],
            };
            let power = 0.015;
            let dt = 2e-4;
            let thermal = caloric(&curves, power);
            let prepared = network.prepare(&thermal, || false).unwrap();
            let body = ThermalTensorPiezoelectricBody::new(
                &ordinary,
                vec![biased_law(&cards[0], &tensor)],
                CELLS[0],
                0.002,
                cx,
            )
            .unwrap();

            // Heat/V = (rho*cp - Tref*b)*dT - b*dT²/2, with
            // b=C*alpha²+p²/epsilon. Both mechanical clamping and the open
            // electrical circuit change the temperature reached by the heat.
            let delta =
                temperature_increment(RHO * CP - REFERENCE_K * CURVATURE, power * dt / VOLUME);
            let temperature = REFERENCE_K + delta;
            let electric_field = -PYRO * delta / EPSILON[2];
            let stress = -STIFFNESS[2] * ALPHA * delta - E33 * electric_field;
            let internal = clamped_internal_energy(temperature, 0.0);
            let candidate = body
                .step_candidate(
                    &thermal,
                    dt,
                    &prepared,
                    &[],
                    &[Drive::Voltage(0.0), Drive::Current(0.0)],
                    &zero,
                    &charges,
                    None,
                    field_controls(),
                    thermal_controls(),
                    &[0.0],
                    cx,
                )
                .unwrap();
            check_report(candidate.report());
            let report = &candidate.report().coupled;
            let state = candidate.body().state();
            near(
                temperatures(candidate.thermal(), &curves)[0],
                temperature,
                2e-8,
                0.0,
            );
            near(
                state.electrodes()[1].voltage_v,
                -LENGTH * electric_field,
                2e-7,
                3e-8,
            );
            near(state.electrodes()[1].charge_c, 0.0, 1e-18, 0.0);
            near(
                state.stored_energy_j(),
                0.5 * VOLUME * CURVATURE * delta * delta,
                2e-16,
                3e-8,
            );
            for cell in state.cells() {
                near(
                    cell.response.electric_field_v_m[2],
                    electric_field,
                    2e-4,
                    3e-8,
                );
                near(cell.response.displacement_c_m2[2], 0.0, 1e-12, 0.0);
                near(cell.response.stress_mandel_pa[2], stress, 0.03, 3e-8);
            }
            let top_force: f64 = geometry
                .top
                .iter()
                .map(|&node| report.electrical.midpoint_support_force_n[node][2])
                .sum();
            near(top_force, AREA * stress, 3e-8, 3e-8);
            near(
                report.source.field_internal_energy_j[1],
                internal,
                2e-13,
                3e-8,
            );
            near(report.source.reversible_heat_j, -internal, 2e-13, 3e-8);
            near(
                report.thermal.stored_energy_change_j,
                RHO * VOLUME * CP * delta,
                2e-13,
                3e-8,
            );
            near(candidate.report().external_heat_j, power * dt, 1e-16, 0.0);
            near(report.electrical.native.mechanical_work_j, 0.0, 1e-16, 0.0);
            near(
                report.electrical.loads.electrical_source_work_j,
                0.0,
                1e-16,
                0.0,
            );
            assert_eq!(report.deposited_joule_j, 0.0);
            assert_eq!(state.displacement_m(), zero);
            assert_eq!(state.velocity_m_s(), zero);
            assert_eq!(body.state().time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
        }
    });
}

#[test]
fn g1_g4_signed_charge_drive_exchanges_reversible_heat_and_preserves_rejected_history() {
    with_cx(|cx, gate| {
        let geometry = cube(1);
        let (card, curve, tensor) = source();
        let materials = [tensor.clone()];
        let model = Model::assemble(
            &geometry.nodes,
            &geometry.tets,
            &materials,
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let n = geometry.nodes.len();
        let zero = vec![[0.0; 3]; n];
        let charges = vec![0.0; n];
        let supports = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let ordinary = Body::new(
            &model,
            &supports,
            &electrodes(&geometry, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charges,
            field_controls(),
            cx,
        )
        .unwrap();
        let cards = [card];
        let curves = [curve];
        let context = point(REFERENCE_K);
        let bindings = cells(&cards, &curves, &context);
        let network = MaterialNetwork {
            cells: &bindings,
            resistors: &[],
            links: &[],
        };
        let thermal = caloric(&curves, 0.0);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let body = ThermalTensorPiezoelectricBody::new(
            &ordinary,
            vec![biased_law(&cards[0], &tensor)],
            CELLS[0],
            0.002,
            cx,
        )
        .unwrap();
        let dt = 1e-4;
        for displacement in [-D_INITIAL, D_INITIAL] {
            let charge = -AREA * displacement;
            let drives = [Drive::Voltage(0.0), Drive::Current(charge / dt)];
            // The endpoint-temperature midpoint voltage pays electrical work
            // V*(D²/2-p*dT*D)/epsilon. Combining that actual work with U and
            // the finite caloric chart leaves this independent quadratic.
            let coupling = PYRO * displacement / EPSILON[2];
            let delta = temperature_increment(
                RHO * CP - REFERENCE_K * CURVATURE + coupling,
                -REFERENCE_K * coupling,
            );
            let temperature = REFERENCE_K + delta;
            let work = VOLUME * (0.5 * displacement * displacement - PYRO * delta * displacement)
                / EPSILON[2];
            let internal = clamped_internal_energy(temperature, displacement);
            let candidate = body
                .step_candidate(
                    &thermal,
                    dt,
                    &prepared,
                    &[],
                    &drives,
                    &zero,
                    &charges,
                    None,
                    field_controls(),
                    thermal_controls(),
                    &[0.0],
                    cx,
                )
                .unwrap();
            check_report(candidate.report());
            let report = &candidate.report().coupled;
            near(
                temperatures(candidate.thermal(), &curves)[0],
                temperature,
                2e-8,
                0.0,
            );
            assert!(
                delta * displacement < 0.0,
                "the pyroelectric sign sets heating versus cooling"
            );
            near(
                candidate.body().state().electrodes()[1].charge_c,
                charge,
                1e-18,
                3e-8,
            );
            near(
                candidate.body().state().electrodes()[1].voltage_v,
                -LENGTH * (displacement - PYRO * delta) / EPSILON[2],
                2e-6,
                3e-8,
            );
            near(
                report.electrical.loads.electrical_source_work_j,
                work,
                2e-15,
                3e-8,
            );
            near(
                report.source.field_internal_energy_j[1],
                internal,
                2e-13,
                3e-8,
            );
            near(
                report.source.reversible_heat_j,
                RHO * VOLUME * CP * delta,
                2e-13,
                3e-8,
            );
            near(internal + RHO * VOLUME * CP * delta, work, 2e-15, 3e-8);
            near(candidate.report().external_heat_j, 0.0, 1e-16, 0.0);
            assert_eq!(report.deposited_joule_j, 0.0);
        }
        let drives = [Drive::Voltage(0.0), Drive::Current(-AREA * D_INITIAL / dt)];
        let mut limited = field_controls();
        limited.maximum_electric_field_v_m = 10.0;
        assert!(
            body.step_candidate(
                &thermal,
                dt,
                &prepared,
                &[],
                &drives,
                &zero,
                &charges,
                None,
                limited,
                thermal_controls(),
                &[0.0],
                cx
            )
            .is_err()
        );
        let retry = body
            .step_candidate(
                &thermal,
                dt,
                &prepared,
                &[],
                &drives,
                &zero,
                &charges,
                None,
                field_controls(),
                thermal_controls(),
                &[0.0],
                cx,
            )
            .unwrap();
        assert!(temperatures(retry.thermal(), &curves)[0] < REFERENCE_K);
        gate.request();
        assert!(matches!(
            body.step_candidate(
                &thermal,
                dt,
                &prepared,
                &[],
                &drives,
                &zero,
                &charges,
                None,
                field_controls(),
                thermal_controls(),
                &[0.0],
                cx
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().voltage_v(), ordinary.state().voltage_v());
        assert_eq!(body.state().displacement_m(), zero);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(temperatures(&thermal, &curves), vec![REFERENCE_K]);
    });
}

#[test]
fn g1_short_circuit_free_expansion_includes_secondary_charge_and_paid_acoustic_work() {
    with_cx(|cx, _| {
        let geometry = cube(1);
        let (card, curve, tensor) = source();
        let materials = [tensor.clone()];
        let model = Model::assemble(
            &geometry.nodes,
            &geometry.tets,
            &materials,
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let n = geometry.nodes.len();
        let zero = vec![[0.0; 3]; n];
        let charges = vec![0.0; n];
        let supports = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let ordinary = Body::new(
            &model,
            &supports,
            &electrodes(&geometry, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charges,
            field_controls(),
            cx,
        )
        .unwrap();
        let cards = [card];
        let curves = [curve];
        let context = point(REFERENCE_K);
        let bindings = cells(&cards, &curves, &context);
        let network = MaterialNetwork {
            cells: &bindings,
            resistors: &[],
            links: &[],
        };
        let dt = 1e-4;
        let delta = 2.0;
        let strain = ALPHA * delta;
        let rate = 2.0 * strain / dt;
        let kinetic = RHO * VOLUME * LENGTH * LENGTH * rate * rate / 6.0;
        let frozen_elastic = 0.5 * VOLUME * STIFFNESS[2] * strain * strain;
        // Prescribe an exact stress-free endpoint, paying all midpoint support
        // work and consistent inertia. This is not a quasistatic assumption
        // about a freely accelerating body. U_final=0, W_support=K-psi_frozen.
        let heat = RHO * VOLUME * CP * delta + frozen_elastic;
        let thermal = caloric(&curves, heat / dt);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let body = ThermalTensorPiezoelectricBody::new(
            &ordinary,
            vec![biased_law(&cards[0], &tensor)],
            CELLS[0],
            0.002,
            cx,
        )
        .unwrap();
        let velocities = geometry
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(i, p)| {
                let v = [0.0, 0.0, rate * p[2]];
                (0..3).map(move |axis| (3 * i + axis, v[axis]))
            })
            .collect::<Vec<_>>();
        let candidate = body
            .step_candidate(
                &thermal,
                dt,
                &prepared,
                &[],
                &[Drive::Voltage(0.0), Drive::Voltage(0.0)],
                &zero,
                &charges,
                Some(&velocities),
                field_controls(),
                thermal_controls(),
                &[0.0],
                cx,
            )
            .unwrap();
        check_report(candidate.report());
        let state = candidate.body().state();
        let report = &candidate.report().coupled;
        near(
            temperatures(candidate.thermal(), &curves)[0],
            REFERENCE_K + delta,
            2e-8,
            0.0,
        );
        near(
            state.electrodes()[1].charge_c,
            -AREA * (PYRO + E33 * ALPHA) * delta,
            1e-18,
            3e-8,
        );
        near(state.electrodes()[1].voltage_v, 0.0, 1e-10, 0.0);
        for cell in state.cells() {
            near(cell.strain_mandel[2], strain, 1e-12, 3e-8);
            near(cell.response.stress_mandel_pa[2], 0.0, 0.03, 0.0);
            near(cell.response.electric_field_v_m[2], 0.0, 2e-4, 0.0);
            near(
                cell.response.displacement_c_m2[2],
                (PYRO + E33 * ALPHA) * delta,
                1e-12,
                3e-8,
            );
        }
        near(state.stored_energy_j(), 0.0, 1e-16, 0.0);
        near(state.kinetic_energy_j(), kinetic, 1e-16, 3e-8);
        near(
            report.electrical.native.support_work_j,
            kinetic - frozen_elastic,
            2e-16,
            3e-8,
        );
        near(
            report.electrical.loads.electrical_source_work_j,
            0.0,
            1e-16,
            0.0,
        );
        near(report.source.field_internal_energy_j[1], 0.0, 2e-13, 0.0);
        near(
            report.source.reversible_heat_j,
            -frozen_elastic,
            2e-13,
            3e-8,
        );
        near(
            report.thermal.stored_energy_change_j,
            RHO * VOLUME * CP * delta,
            2e-13,
            3e-8,
        );
        near(candidate.report().external_heat_j, heat, 1e-16, 0.0);
        assert_eq!(report.deposited_joule_j, 0.0);

        let mut faces = BTreeSet::new();
        for tet in &geometry.tets {
            for omit in 0..4 {
                let face = (0..4)
                    .filter(|&i| i != omit)
                    .map(|i| tet[i])
                    .collect::<Vec<_>>();
                if face.iter().all(|node| geometry.top.contains(node)) {
                    let mut face = [face[0], face[1], face[2]];
                    face.sort_unstable();
                    faces.insert(face);
                }
            }
        }
        let faces = faces.into_iter().collect::<Vec<_>>();
        let maximum_wavenumber = core::f64::consts::PI / LENGTH;
        // The trace belongs to the ORIGINAL model. Thermal trials must retain
        // its geometry while rebinding it to their changed material operators.
        let port = model
            .plane_wave_acoustic_port(
                &faces,
                AREA,
                PlaneWavePortControls {
                    maximum_surface_dofs: 64,
                    coplanarity_tolerance: 1e-10,
                    maximum_wavenumber_rad_m: maximum_wavenumber,
                    maximum_relative_amplitude: 1e-2,
                },
                cx,
            )
            .unwrap();
        let (density, sound_speed) = (1000.0, 1480.0);
        let medium = FluidAcousticState::try_new(density, sound_speed).unwrap();
        let wave = TransientAcousticPort {
            port: &port,
            medium: &medium,
            midpoint_incident_pressure_pa: 0.0,
            maximum_angular_frequency_rad_s: 0.5 * sound_speed * maximum_wavenumber,
        };
        // Uniform top velocity gives W_sound=dt*rho*c*A*v_mid² independently
        // of the native trace weights. Supports pay this outgoing wave work;
        // it cannot heat the ceramic or change the required external heat.
        let midpoint_velocity = 0.5 * rate * LENGTH;
        let wave_work = dt * density * sound_speed * AREA * midpoint_velocity * midpoint_velocity;
        for inputs in [&[][..], std::slice::from_ref(&wave)] {
            let acoustic = body
                .step_with_acoustic_ports_candidate(
                    &thermal,
                    dt,
                    inputs,
                    &prepared,
                    &[],
                    &[Drive::Voltage(0.0), Drive::Voltage(0.0)],
                    &zero,
                    &charges,
                    Some(&velocities),
                    field_controls(),
                    thermal_controls(),
                    &[0.0],
                    cx,
                )
                .unwrap();
            check_report(acoustic.report());
            let actual = &acoustic.report().coupled;
            let outward = if inputs.is_empty() { 0.0 } else { wave_work };
            near(
                temperatures(acoustic.thermal(), &curves)[0],
                REFERENCE_K + delta,
                2e-8,
                0.0,
            );
            near(
                acoustic.body().state().electrodes()[1].charge_c,
                -AREA * (PYRO + E33 * ALPHA) * delta,
                1e-18,
                3e-8,
            );
            assert_eq!(
                acoustic.body().state().displacement_m(),
                state.displacement_m()
            );
            assert_eq!(acoustic.body().state().velocity_m_s(), state.velocity_m_s());
            near(
                actual.electrical.native.support_work_j,
                kinetic - frozen_elastic + outward,
                2e-16,
                3e-8,
            );
            near(
                actual.electrical.native.acoustic_outward_work_j,
                outward,
                2e-18,
                3e-8,
            );
            near(
                actual.source.reversible_heat_j,
                -frozen_elastic,
                2e-13,
                3e-8,
            );
            near(acoustic.report().external_heat_j, heat, 1e-16, 0.0);
            assert_eq!(actual.deposited_joule_j, 0.0);
            assert_eq!(actual.electrical.acoustic_ports.len(), inputs.len());
            if let Some(sound) = actual.electrical.acoustic_ports.first() {
                near(
                    sound.midpoint_volume_velocity_m3_s,
                    AREA * midpoint_velocity,
                    1e-18,
                    3e-8,
                );
                near(sound.outgoing_energy_j, wave_work, 2e-18, 3e-8);
                near(sound.net_outward_work_j, wave_work, 2e-18, 3e-8);
                assert_eq!(sound.incident_energy_j, 0.0);
            } else {
                assert_eq!(
                    acoustic.thermal().cells()[&CELLS[0]].specific_h_j_kg,
                    candidate.thermal().cells()[&CELLS[0]].specific_h_j_kg
                );
            }
        }
    });
}
