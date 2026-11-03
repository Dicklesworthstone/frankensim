//! Manufactured original sources and independent conservative physical limits.
use super::*;
use crate::conjugate::domain::Voxel;
use crate::conjugate::flow::{
    CoexistenceVolumeControls, DensityProjectionControls, MaterialAdvectionControls,
    ThermalVolumeControls,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};

fn card(axis: &str, values: &[(f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("k", KD),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: if axis == "T" {
                    Temperature::DIMS
                } else {
                    Dims::NONE
                },
                knots: values.to_vec(),
                dims: KD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured conductive density-volume fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermal-volume material".into(),
            phase: "declared original caloric chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn curve(card: &MaterialCard, values: &[(f64, f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    let knots = values
        .iter()
        .map(|&(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect();
    if values.iter().all(|value| value.2 == 1.0) {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Liquid,
            knots,
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots).unwrap()
    }
}

fn source<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

fn controls(dt: f64) -> ConductiveThermalVolumeControls {
    ConductiveThermalVolumeControls {
        volume: PiecewiseThermalVolumeControls {
            thermal: ThermalVolumeControls {
                volume: CoexistenceVolumeControls {
                    transport: MaterialAdvectionControls {
                        projection: DensityProjectionControls {
                            duration_s: dt,
                            relative_tolerance: 1e-12,
                            maximum_volume_defect_m3_s: 1e-10,
                            maximum_energy_defect_j: 1e-8,
                        },
                        max_cells: 100,
                        max_connections: 1000,
                        max_iterations: 1000,
                        relative_tolerance: 1e-12,
                        maximum_mass_defect_kg: 1e-10,
                        maximum_enthalpy_defect_j: 1e-9,
                    },
                    maximum_cell_volume_defect_m3: 1e-10,
                },
                max_volume_iterations: 64,
                source_relaxation: 1.0,
            },
            max_mixing_intervals: 10_000,
        },
        max_conductive_faces: 300,
    }
}

fn input<'a>(
    mass: &'a [f64],
    states: &'a [Option<EquilibriumPhaseState>],
    increments: &'a [f64],
) -> CoexistenceVolumeInput<'a> {
    CoexistenceVolumeInput {
        mass_kg: mass,
        phase_states: states,
        net_enthalpy_increment_j: increments,
        boundary_inflow: [None; 6],
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}",
    );
}

fn invalid(result: Result<ConductiveThermalVolumeStep, ChtError>, reason: &str) {
    match result {
        Err(ChtError::InvalidInput { reason: actual, .. }) => {
            assert!(
                actual.contains(reason),
                "{actual:?} does not contain {reason:?}"
            );
        }
        other => panic!("expected {reason:?}, got {other:?}"),
    }
}

#[test]
fn sourced_heat_drives_exact_sealed_density_flow_and_paired_donor_enthalpy() {
    // Independent model: T=300+h, k=T-299, v=1+h/8. At h=[1,3],
    // G=1/(1/4+1/8)=8/3. dt=3/32 gives Q=1/2 J and volume flux
    // q=(Q/8)/dt=2/3 m3/s. The left donor exports (1/16)*M_left_new.
    // This determines M=[128/153,1312/1683] and h=[25/16,371/164].
    let material = card("T", &[(300.0, 1.0), (304.0, 5.0)]);
    for axis in 0..3 {
        for datum in [0.0, 1024.0] {
            let caloric = curve(
                &material,
                &[
                    (datum, 300.0, 1.0, 1.0),
                    (datum + 4.0, 304.0, 1.0, 2.0 / 3.0),
                ],
            );
            let conductivity = source(&material, &caloric, [datum, datum + 4.0]);
            let mut dims = [1; 3];
            dims[axis] = 2;
            let domain = VoxelDomain::new(dims[0], dims[1], dims[2], 1.0).unwrap();
            let field = FlowField::quiescent(&domain);
            let mass = [8.0 / 9.0, 8.0 / 11.0];
            let states =
                [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(datum + h).unwrap()));
            let result = field
                .advance_isotropic_conductive_thermal_volume(
                    &domain,
                    [FlowFace::Wall; 6],
                    &conductivity,
                    input(&mass, &states, &[0.0; 2]),
                    controls(3.0 / 32.0),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            assert_eq!(
                result.report.conductivity_source_identity,
                conductivity.identity()
            );
            assert_eq!(result.report.conductive_faces, 1);
            close(result.report.exchanged_heat_j, 0.5);
            close(result.report.maximum_enthalpy_mixing_weight, 11.0 / 32.0);
            close(result.conductive_enthalpy_increment_j[0], 0.5);
            close(result.conductive_enthalpy_increment_j[1], -0.5);
            let volume = &result.thermal.volume;
            let transport = &volume.transport;
            for (actual, expected) in transport
                .mass_kg
                .iter()
                .zip([128.0 / 153.0, 1312.0 / 1683.0])
            {
                close(*actual, expected);
            }
            for (c, expected) in [25.0 / 16.0, 371.0 / 164.0].into_iter().enumerate() {
                close(transport.specific_enthalpy_j_kg[c] - datum, expected);
                let state = volume.phase_states[c].unwrap();
                close(transport.mass_kg[c] / state.bulk_density_kg_m3(), 1.0);
                assert_eq!(state.phase_curve_identity(), caloric.identity());
            }
            close(
                transport.projection.field.outward(&domain, 0, 2 * axis + 1),
                2.0 / 3.0,
            );
            close(volume.volume_source_m3_s[0], 2.0 / 3.0);
            close(volume.volume_source_m3_s[1], -2.0 / 3.0);
            // The native projection freezes inertia at the actual original
            // cell inventories; the new EOS volume is a continuity source.
            let dual_mass = 0.5 * (8.0 / 9.0 + 8.0 / 11.0);
            let pressure = -(2.0 / 3.0) * dual_mass / (3.0 / 32.0);
            close(transport.projection.pressure_correction_pa[1], pressure);
            close(
                volume.report.required_source_pressure_work_j,
                -(3.0 / 32.0) * pressure * (2.0 / 3.0),
            );
            assert!(volume.report.required_source_pressure_work_j > 0.0);
            assert!(transport.projection.report.projection_loss_j > 0.0);
            close(
                transport.report.mass_before_kg,
                transport.report.mass_after_kg,
            );
            close(transport.report.carried_enthalpy_j, 0.0);
            close(volume.report.material_enthalpy_change_j, 0.0);
            close(volume.report.requested_enthalpy_increment_j, 0.0);
            assert_eq!(result.thermal.source_iterations, 1);
            assert_eq!(field, FlowField::quiescent(&domain));
            assert_eq!(mass, [8.0 / 9.0, 8.0 / 11.0]);
        }
    }
}

#[test]
fn graded_half_cell_resistance_and_actual_mass_recover_conservative_conduction() {
    // With widths [1,2] and face area 6, actual half-cell resistance is
    // .5/2+1/8=3/8, hence G=16 W/K. Q=.01*16*40=6.4 J. Masses are
    // rho*V=[12,24], so the two specific increments are Q/12 and -Q/24.
    let material = card("T", &[(300.0, 2.0), (340.0, 8.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (4.0, 340.0, 1.0, 2.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    for axis in 0..3 {
        let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let states = [0.0, 4.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
        let result = FlowField::quiescent(&domain)
            .advance_isotropic_conductive_thermal_volume(
                &domain,
                [FlowFace::Wall; 6],
                &conductivity,
                input(&[12.0, 24.0], &states, &[0.0; 2]),
                controls(0.01),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(result.conductive_enthalpy_increment_j[0], 6.4);
        close(result.conductive_enthalpy_increment_j[1], -6.4);
        close(result.report.maximum_enthalpy_mixing_weight, 2.0 / 15.0);
        let volume = &result.thermal.volume;
        close(volume.transport.specific_enthalpy_j_kg[0], 6.4 / 12.0);
        close(volume.transport.specific_enthalpy_j_kg[1], 4.0 - 6.4 / 24.0);
        assert_eq!(volume.transport.mass_kg, vec![12.0, 24.0]);
        assert_eq!(volume.volume_source_m3_s, vec![0.0; 2]);
        close(volume.report.required_source_pressure_work_j, 0.0);
        close(volume.report.material_enthalpy_change_j, 0.0);
    }
}

#[test]
fn liquid_cooling_and_latent_heating_keep_one_original_density_and_latent_account() {
    let material = card("f", &[(0.0, 2.0), (1.0, 6.0)]);
    let caloric = curve(
        &material,
        &[
            (-2.0, 300.0, 0.0, 2.0),
            (0.0, 400.0, 0.0, 2.0),
            (4.0, 400.0, 1.0, 1.0),
            (6.0, 500.0, 1.0, 1.0),
        ],
    );
    let conductivity = source(&material, &caloric, [-2.0, 6.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let states = [1.0, 5.5].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    // Cold k=3, hot k=6 give G=4. Q=.01*4*75=3, h*=[23/8,5/2].
    // Original M=[8/5,1] changes occupied volume by [3/8,-3/16].
    // The left-wall constraint fixes internal flux 37.5 and outlet 18.75.
    let result = FlowField::quiescent(&domain)
        .advance_isotropic_conductive_thermal_volume(
            &domain,
            faces,
            &conductivity,
            input(&[8.0 / 5.0, 1.0], &states, &[0.0; 2]),
            controls(0.01),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(result.conductive_enthalpy_increment_j[0], 3.0);
    close(result.conductive_enthalpy_increment_j[1], -3.0);
    close(result.report.maximum_enthalpy_mixing_weight, 2.0 / 3.0);
    let volume = &result.thermal.volume;
    close(volume.volume_source_m3_s[0], 37.5);
    close(volume.volume_source_m3_s[1], -18.75);
    close(
        volume.transport.projection.field.outward(&domain, 0, 1),
        37.5,
    );
    close(
        volume.transport.projection.field.outward(&domain, 1, 1),
        18.75,
    );
    for (c, (mass, h)) in [(64.0 / 55.0, 23.0 / 8.0), (1264.0 / 1045.0, 413.0 / 158.0)]
        .into_iter()
        .enumerate()
    {
        close(volume.transport.mass_kg[c], mass);
        close(volume.transport.specific_enthalpy_j_kg[c], h);
        let state = volume.phase_states[c].unwrap();
        close(state.temperature_k(), 400.0);
        close(state.liquid_mass_fraction(), h / 4.0);
        close(mass / state.bulk_density_kg_m3(), 1.0);
    }
    close(volume.transport.report.outgoing_mass_kg, 237.0 / 1045.0);
    close(volume.transport.report.carried_enthalpy_j, -1239.0 / 2090.0);
    close(volume.report.material_enthalpy_change_j, -1239.0 / 2090.0);
    assert_eq!(result.thermal.source_iterations, 1);

    // At the accepted endpoint both cells lie on the true isothermal plateau.
    // Re-evaluating their original phase conductivities produces no heat or
    // latent change even though their h, k, phase fractions and densities differ.
    // A Free outlet remains a possible inlet: projecting the prior outflow
    // to zero can leave a signed roundoff flux. Author the same accepted
    // outlet material for backflow; no flux is clipped or treated as heat.
    let mut continuation = input(&volume.transport.mass_kg, &volume.phase_states, &[0.0; 2]);
    continuation.boundary_inflow[1] = volume.phase_states[1];
    let next = volume
        .transport
        .projection
        .field
        .advance_isotropic_conductive_thermal_volume(
            &domain,
            faces,
            &conductivity,
            continuation,
            controls(0.01),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(next.conductive_enthalpy_increment_j, vec![0.0; 2]);
    close(next.report.exchanged_heat_j, 0.0);
    close(next.report.maximum_enthalpy_mixing_weight, 0.0);
    for c in 0..2 {
        close(
            next.thermal.volume.transport.mass_kg[c],
            volume.transport.mass_kg[c],
        );
        close(
            next.thermal.volume.transport.specific_enthalpy_j_kg[c],
            volume.transport.specific_enthalpy_j_kg[c],
        );
    }
}

#[test]
fn authored_net_heat_is_added_once_to_original_state_conduction() {
    let material = card("T", &[(300.0, 2.0), (304.0, 2.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (4.0, 304.0, 1.0, 2.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let states = [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let result = FlowField::quiescent(&domain)
        .advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            input(&[2.0; 2], &states, &[0.25, -0.125]),
            controls(0.125),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.conductive_enthalpy_increment_j, vec![0.5, -0.5]);
    let volume = &result.thermal.volume;
    close(volume.transport.specific_enthalpy_j_kg[0], 1.375);
    close(volume.transport.specific_enthalpy_j_kg[1], 2.6875);
    close(volume.report.material_enthalpy_change_j, 0.125);
    close(volume.report.requested_enthalpy_increment_j, 0.125);
    close(result.report.increment_assembly_defect_j, 0.0);
}

#[test]
fn monotonicity_limits_sum_all_neighbor_weights_without_latent_heat_capacity() {
    let material = card("T", &[(300.0, 2.0), (304.0, 2.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 1.0), (4.0, 304.0, 1.0, 1.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [3.0, 1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let data = input(&[1.0; 3], &states, &[0.0; 3]);
    // Each face alone has weight .6; the middle cell's sum 1.2 must refuse.
    invalid(
        field.advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            data,
            controls(0.3),
            &CancelGate::new_clock_free(),
        ),
        "mixing weight exceeds one",
    );
    // Exactly one is admitted, without a hidden safety factor or h clipping.
    let result = field
        .advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            data,
            controls(0.25),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        result.thermal.volume.transport.specific_enthalpy_j_kg,
        vec![2.0, 3.0, 2.0]
    );
    close(result.report.maximum_enthalpy_mixing_weight, 1.0);
    assert_eq!(field, FlowField::quiescent(&domain));
}

#[test]
fn zero_source_conductivity_and_solid_faces_insulate_without_material_heat() {
    for zero_source in [true, false] {
        let k = if zero_source { 0.0 } else { 2.0 };
        let material = card("T", &[(300.0, k), (304.0, k)]);
        let caloric = curve(
            &material,
            &[(0.0, 300.0, 1.0, 1.0), (4.0, 304.0, 1.0, 2.0 / 3.0)],
        );
        let conductivity = source(&material, &caloric, [0.0, 4.0]);
        let mut domain = VoxelDomain::new(if zero_source { 2 } else { 3 }, 1, 1, 1.0).unwrap();
        let mut states = vec![
            Some(caloric.state_at_specific_enthalpy(1.0).unwrap()),
            Some(caloric.state_at_specific_enthalpy(3.0).unwrap()),
        ];
        let mut mass = vec![8.0 / 9.0, 8.0 / 11.0];
        if !zero_source {
            domain.set(1, 0, 0, Voxel::Solid(0));
            states.insert(1, None);
            mass.insert(1, 0.0);
        }
        let mut limits = controls(10.0);
        limits.max_conductive_faces = usize::from(zero_source);
        let result = FlowField::quiescent(&domain)
            .advance_isotropic_conductive_thermal_volume(
                &domain,
                [FlowFace::Wall; 6],
                &conductivity,
                input(&mass, &states, &vec![0.0; mass.len()]),
                limits,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.report.conductive_faces, usize::from(zero_source));
        assert_eq!(
            result.conductive_enthalpy_increment_j,
            vec![0.0; mass.len()]
        );
        assert_eq!(result.thermal.volume.transport.mass_kg, mass);
        assert_eq!(result.thermal.volume.phase_states, states);
        close(result.report.exchanged_heat_j, 0.0);
    }
}

#[test]
fn source_chart_budget_and_cancel_refusals_leave_the_original_state_untouched() {
    let material = card("T", &[(300.0, 2.0), (304.0, 2.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 1.0), (4.0, 304.0, 1.0, 1.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    let narrow = source(&material, &caloric, [1.0, 3.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let data = input(&[1.0; 2], &states, &[0.0; 2]);
    let gate = CancelGate::new_clock_free();
    let call = |source, data, limits| {
        field.advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            source,
            data,
            limits,
            &gate,
        )
    };
    let mut limits = controls(0.125);
    limits.max_conductive_faces = 0;
    invalid(call(&conductivity, data, limits), "face-work budget");
    limits = controls(0.125);
    limits.volume.thermal.volume.transport.max_cells = 1;
    invalid(call(&conductivity, data, limits), "cell budget");
    invalid(
        call(
            &conductivity,
            CoexistenceVolumeInput {
                mass_kg: &[],
                ..data
            },
            controls(0.125),
        ),
        "shape",
    );
    invalid(
        call(
            &narrow,
            CoexistenceVolumeInput {
                net_enthalpy_increment_j: &[3.0, 0.0],
                ..data
            },
            controls(0.125),
        ),
        "calorically advanced",
    );
    let mut incoming = data;
    incoming.boundary_inflow[0] = Some(caloric.state_at_specific_enthalpy(0.0).unwrap());
    invalid(call(&narrow, incoming, controls(0.125)), "source/chart");
    let other = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (4.0, 304.0, 1.0, 2.0)]);
    let foreign = [1.0, 3.0].map(|h| Some(other.state_at_specific_enthalpy(h).unwrap()));
    invalid(
        call(
            &conductivity,
            CoexistenceVolumeInput {
                phase_states: &foreign,
                ..data
            },
            controls(0.125),
        ),
        "source/chart",
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        field.advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            data,
            controls(0.125),
            &cancelled,
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(
        states,
        [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()))
    );
}

#[test]
fn small_conductive_heat_lost_beside_authored_heat_consumes_the_shared_budget() {
    let material = card("T", &[(300.0, 2.0), (500.0, 2.0)]);
    let caloric = curve(
        &material,
        &[
            (0.0, 300.0, 1.0, 1.0),
            (4.0, 304.0, 1.0, 1.0),
            (18_014_398_509_481_984.0, 500.0, 1.0, 1.0),
        ],
    );
    let conductivity = source(&material, &caloric, [0.0, 18_014_398_509_481_984.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    // The physical face heat is 0.5 J, but 2^53+0.5 rounds to 2^53.
    // All enthalpies remain in-domain. The addition residual itself must
    // refuse before native transport can treat the lost heat as a valid source.
    invalid(
        field.advance_isotropic_conductive_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            input(&[1.0; 2], &states, &[9_007_199_254_740_992.0, 0.0]),
            controls(0.125),
            &CancelGate::new_clock_free(),
        ),
        "heat assembly exceeds",
    );
    assert_eq!(field, FlowField::quiescent(&domain));
}
