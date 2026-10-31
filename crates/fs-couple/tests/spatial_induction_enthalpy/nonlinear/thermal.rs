//! G1/G3/G4: sourced H(B,T), actual 3D currents and complete finite heat.
//! Independent scalar cube and two-temperature balances qualify implementation
//! accuracy using manufactured sources, not measured-material fidelity.
use super::*;
use fs_couple::circuit_enthalpy::magnetic::caloric::MagneticConstitutiveThermalLaw;
use fs_material::field_storage::nonlinear::thermal::ThermalReversibleMagneticMaterial;
#[path = "../../../../fs-material/tests/nonlinear_magnetic_thermal/common.rs"]
mod scale_source;

const CELL: CellId = CellId(7);
const SINK: CellId = CellId(8);
const DT: f64 = 0.05;
const A0: f64 = 1.2;
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;

fn with_scale(base: MaterialCard, slope: f64) -> MaterialCard {
    let source = scale_source::card(slope);
    let claim = source
        .claims()
        .claims_ordered()
        .find(|(_, claim)| claim.key.name() == scale_source::SCALE)
        .unwrap()
        .1;
    let mut claims = base.claims().clone();
    claims.insert_claim(claim.clone()).unwrap();
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn law(card: &MaterialCard) -> ThermalReversibleMagneticMaterial {
    let original = conductor(card, 300.0, BASIS);
    let mut requirements = original
        .material()
        .properties()
        .iter()
        .map(|p| p.requirement().clone())
        .collect::<Vec<_>>();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            scale_source::SCALE,
            fs_qty::Dims([0; 6]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        card,
        &context(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let reference = ReversibleMagneticMaterial::resolve(
        card,
        &context(300.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ThermalReversibleMagneticMaterial::new(
        &reference,
        &bundle,
        &scale_source::scale_key(),
        "T",
        [290.0, 360.0],
    )
    .unwrap()
}

fn problem(
    card: &MaterialCard,
    law: &ThermalReversibleMagneticMaterial,
) -> SpatialNonlinearVectorMagnetostatics {
    let geometry = super::problem(card, BASIS);
    let material = law.at_temperature(300.0).unwrap().magnetic().clone();
    SpatialNonlinearVectorMagnetostatics::new(
        geometry.positions_m().to_vec(),
        geometry.mesh().tets.clone(),
        vec![SpatialNonlinearMaterial::Reversible(material); 6],
        assembly(),
        &mut || false,
    )
    .unwrap()
}

fn bindings(
    law: &ThermalReversibleMagneticMaterial,
) -> Vec<SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>> {
    (0..6)
        .map(|tetrahedron| SpatialMagneticThermalLaw {
            tetrahedron,
            law: MagneticConstitutiveThermalLaw::Separable(law.clone()),
            maximum_flux_density_t: 3.0,
        })
        .collect()
}

fn initial<'m>(
    problem: &'m SpatialNonlinearVectorMagnetostatics,
    law: &ThermalReversibleMagneticMaterial,
) -> SpatialNonlinearInduction<'m> {
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        law.at_temperature(300.0).unwrap().material().clone(),
        BASIS,
    )
    .unwrap();
    initial_with_conductor(problem, conductor, A0, true)
}

// Unit-cube Whitney mass after charge elimination is 3 sigma/40. Each caloric
// cell has mass 1 kg and cp=4 J/(kg K). For reference energy F(a) and affine
// scale s=1+b(T-300), U=(1-300b)F and Qrev=b*T1*(F1-F0).
fn reference(t0: [f64; 2], a0: f64, slope: f64) -> ([f64; 2], f64, f64, f64, f64, f64) {
    let amplitude = |t| {
        bracket(0.0, a0, |a| {
            let mass = 3.0 * sigma(t, 0.2, 0.009, BASIS) / 40.0;
            mass * (a - a0) / DT + 2.0 * (1.0 + slope * (t - 300.0)) * h(2.0 * a, true)
        })
    };
    let sink = |t| {
        bracket(290.0, 359.0, |other| {
            4.0 * (other - t0[1]) + DT * heat::transfer(t, other, 0.4, 0.2)
        })
    };
    let t = bracket(290.0, 359.0, |t| {
        let a = amplitude(t);
        let q = 3.0 * sigma(t, 0.2, 0.009, BASIS) / 40.0 * (a - a0).powi(2) / DT;
        let reversible = slope * t * (energy(a, true) - energy(a0, true));
        4.0 * (t - t0[0]) - q - reversible - DT * heat::transfer(t, sink(t), 0.4, 0.2)
    });
    let a = amplitude(t);
    let scale = 1.0 + slope * (t - 300.0);
    let q = 3.0 * sigma(t, 0.2, 0.009, BASIS) / 40.0 * (a - a0).powi(2) / DT;
    let reversible = slope * t * (energy(a, true) - energy(a0, true));
    let numerical =
        scale * (2.0 * h(2.0 * a, true) * (a - a0) - energy(a, true) + energy(a0, true));
    let parameter = slope * (t - t0[0]) * energy(a0, true);
    ([t, sink(t)], a, q, reversible, numerical, parameter)
}

#[test]
fn sourced_temperature_and_complete_spatial_heat_match_independent_balance() {
    let mut endpoints = Vec::new();
    for slope in [0.0, 0.001, -0.001] {
        endpoints.push(heat::with_heat_material(
            true,
            0.4,
            0.2,
            500.0,
            None,
            |base| with_scale(base, slope),
            |cx, _, assembly, mut thermal, card| {
                let law = law(card);
                let problem = problem(card, &law);
                let mut field = initial(&problem, &law);
                let bindings = bindings(&law);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                for step in 1..=2 {
                    let a0 = field.potential_wb()[diagonal(&problem)];
                    let expected = reference(
                        [
                            thermal.temperature_k(CELL).unwrap(),
                            thermal.temperature_k(SINK).unwrap(),
                        ],
                        a0,
                        slope,
                    );
                    let advance = || {
                        prepared
                            .step_spatial_nonlinear_induction_thermal_candidate(
                                cx,
                                &thermal,
                                step as f64 * DT,
                                &field,
                                &[[0.0; 3]; 6],
                                &[[0.0; 3]; 8],
                                &regions(),
                                1e-12,
                                &bindings,
                                thermal_controls(),
                                field_controls(),
                                &[0.0; 2],
                            )
                            .unwrap()
                    };
                    let next = advance();
                    let replay = advance();
                    assert_eq!(next.field().potential_wb(), replay.field().potential_wb());
                    assert_eq!(
                        next.thermal().energy_j().to_bits(),
                        replay.thermal().energy_j().to_bits()
                    );
                    near(
                        next.field().potential_wb()[diagonal(&problem)],
                        expected.1,
                        3e-9,
                    );
                    for (cell, t) in [CELL, SINK].into_iter().zip(expected.0) {
                        near(next.thermal().temperature_k(cell).unwrap(), t, 3e-8);
                    }
                    let report = next.report();
                    let native = &report.coupled.electrical;
                    near(native.heat_j, expected.2, 3e-8);
                    near(report.coupled.deposited_joule_j, expected.2, 3e-8);
                    near(report.coupled.source.reversible_heat_j, expected.3, 3e-8);
                    near(native.numerical_dissipation_j, expected.4, 3e-8);
                    near(native.parameter_work_j, expected.5, 3e-8);
                    near(
                        report.coupled.thermal.stored_energy_change_j,
                        expected.2 + expected.3,
                        4e-8,
                    );
                    near(
                        report.coupled.source.field_internal_energy_j[0],
                        (1.0 - 300.0 * slope) * energy(a0, true),
                        3e-8,
                    );
                    near(
                        report.coupled.source.field_internal_energy_j[1],
                        (1.0 - 300.0 * slope) * energy(expected.1, true),
                        3e-8,
                    );
                    near(report.external_heat_j, 0.0, 3e-8);
                    near(report.total_energy_residual_j, 0.0, 3e-8);
                    assert_eq!(report.coupled.source.fields.len(), 6);
                    assert!(native.numerical_dissipation_j > 1e-5);
                    let sources = &report.coupled.source.network;
                    assert_eq!(sources.thermal_contacts.len(), 1);
                    assert_eq!(sources.radiation.as_ref().unwrap().surfaces().len(), 2);
                    for source in next.field().conductors().iter().flatten() {
                        near(source.material().query().axes()["T"], expected.0[0], 3e-8);
                        assert!(source.material().property(scale_source::SCALE).is_some());
                    }
                    let (accepted_thermal, accepted_field, _) = next.into_parts();
                    thermal = accepted_thermal;
                    field = accepted_field;
                }
                [
                    thermal.temperature_k(CELL).unwrap(),
                    field.potential_wb()[diagonal(&problem)],
                ]
            },
        ));
    }
    for changed in &endpoints[1..] {
        assert!((changed[0] - endpoints[0][0]).abs() > 1e-5);
        assert!((changed[1] - endpoints[0][1]).abs() > 1e-7);
    }
}

#[test]
fn bare_spatial_entropy_owner_refuses_missing_laws_and_late_cancellation() {
    heat::with_heat_material(
        true,
        0.4,
        0.2,
        500.0,
        None,
        |base| with_scale(base, 0.001),
        |_, _, assembly, thermal, card| {
            let law = law(card);
            let problem = problem(card, &law);
            let field = initial(&problem, &law);
            let bindings = bindings(&law);
            let regions = regions();
            let network = SpatialMagneticThermalNetwork {
                cells: assembly.electrical.bulk.cells,
                links: assembly.electrical.bulk.links,
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let calls = std::cell::Cell::new(0usize);
            let next = network
                .step_nonlinear_thermal_candidate(
                    &thermal,
                    DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &bindings,
                    thermal_controls(),
                    field_controls(),
                    &[0.0; 2],
                    || {
                        calls.set(calls.get() + 1);
                        false
                    },
                )
                .unwrap();
            assert!(calls.get() > 10);
            near(next.report().total_energy_residual_j, 0.0, 3e-8);
            assert!(next.report().coupled.source.reversible_heat_j.abs() > 1e-5);
            let advance = |laws: &[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]| {
                network.step_nonlinear_thermal_candidate(
                    &thermal,
                    DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    laws,
                    thermal_controls(),
                    field_controls(),
                    &[0.0; 2],
                    || false,
                )
            };
            assert!(advance(&[]).is_err());
            let mut outside = bindings.clone();
            outside[0].maximum_flux_density_t = 0.5;
            assert!(advance(&outside).is_err());
            assert!(
                network
                    .step_nonlinear_candidate(
                        &thermal,
                        DT,
                        &field,
                        &[[0.0; 3]; 6],
                        &[[0.0; 3]; 8],
                        thermal_controls(),
                        field_controls(),
                        &[0.0; 2],
                        || false
                    )
                    .is_err()
            );
            let stop = calls.get() * 3 / 4;
            calls.set(0);
            let cancelled = network.step_nonlinear_thermal_candidate(
                &thermal,
                DT,
                &field,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &bindings,
                thermal_controls(),
                field_controls(),
                &[0.0; 2],
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= stop
                },
            );
            assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(field.time_s(), 0.0);
            near(field.potential_wb()[diagonal(&problem)], A0, 2e-10);
            near(thermal.temperature_k(CELL).unwrap(), 300.0, 1e-12);
            near(thermal.temperature_k(SINK).unwrap(), 340.0, 1e-12);
            assert_eq!(
                field.conductors()[0]
                    .as_ref()
                    .unwrap()
                    .material()
                    .query()
                    .axes()["T"],
                300.0
            );
        },
    );
}
