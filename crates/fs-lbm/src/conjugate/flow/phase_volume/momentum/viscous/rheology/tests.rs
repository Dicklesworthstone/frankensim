//! Real source queries and native strain on the evolving material inventory.

use super::super::super::tests::{close, controls as momentum_controls};
use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::{CoexistenceVolumeInput, MaterialMomentumViscosityInput};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::liquid::{KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY as NU_KEY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::phase::HeatCapacitySourceQuery;
use fs_material::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, resolve_phase_transport_state_with_shear_for_volume_change,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_KEY, THERMAL_CONDUCTIVITY_PROPERTY as K_KEY,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Frequency, QuantitySpec, Temperature};

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), t)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 3.0)
        .unwrap()
}

fn fixture(changing_density: bool, shear_dependent: bool) -> (MaterialCard, ResolvedMeltingCurve) {
    let cp = QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = QuantitySpec::dimensional(Density::DIMS);
    let latent = QuantitySpec::dimensional(Dims([2, 0, -2, 0, 0, 0]));
    let nu = QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS);
    let k = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let densities = [2.0, if changing_density { 1.0 } else { 2.0 }];
    let viscosities = if shear_dependent {
        [0.25, 0.75]
    } else {
        [0.5; 2]
    };
    let mut claims = ClaimSet::new();
    for (name, quantity, axis, range, values) in [
        ("solid_cp", cp, "T", [299.0, 300.0], [1.0; 2]),
        ("solid_density", rho, "T", [299.0, 300.0], [4.0; 2]),
        (CP_KEY, cp, "T", [300.0, 302.0], [1.0; 2]),
        ("density", rho, "T", [300.0, 302.0], densities),
        (NU_KEY, nu, "rate", [0.0, 4.0], viscosities),
        (K_KEY, k, "T", [299.0, 302.0], [1.0; 2]),
        ("melting_point", absolute(), "T", [299.0, 302.0], [300.0; 2]),
        ("fusion", latent, "T", [299.0, 302.0], [1.0; 2]),
    ] {
        let temperatures = if axis == "T" { range } else { [299.0, 302.0] };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: if axis == "T" {
                        Temperature::DIMS
                    } else {
                        Frequency::DIMS
                    },
                    knots: vec![(range[0], values[0]), (range[1], values[1])],
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), temperatures[0], temperatures[1])
                    .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 4.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED caloric/rheology coupling fixture, not measurements"
                        .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured pure substance".into(),
            phase: "solid/liquid".into(),
            process: "density and strain integration test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let solid_cp = PropertyKey::with_quantity("solid_cp", cp);
    let solid_rho = PropertyKey::with_quantity("solid_density", rho);
    let liquid_cp = PropertyKey::with_quantity(CP_KEY, cp);
    let liquid_rho = PropertyKey::with_quantity("density", rho);
    let tm = PropertyKey::with_quantity("melting_point", absolute());
    let fusion = PropertyKey::with_quantity("fusion", latent);
    let points = [point(299.0), point(300.0), point(302.0)];
    let melting = ResolvedMeltingCurve::try_new(
        &card,
        MeltingSourceQuery {
            solid: HeatCapacitySourceQuery {
                heat_capacity: &solid_cp,
                density: &solid_rho,
                temperature_axis: "T",
                lower: &points[0],
                upper: &points[1],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            liquid: HeatCapacitySourceQuery {
                heat_capacity: &liquid_cp,
                density: &liquid_rho,
                temperature_axis: "T",
                lower: &points[1],
                upper: &points[2],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            melting_point: &tm,
            latent_heat_fusion: &fusion,
            transition: &points[1],
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: None,
        },
        -2.0,
        1e-12,
    )
    .unwrap();
    (card, melting)
}

fn resolve(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    rate: f64,
) -> Result<ResolvedPhaseTransportState, ChtError> {
    resolve_phase_transport_state_with_shear_for_volume_change(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material: card,
                context: &point(301.5),
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            viscosity_basis: LiquidViscosityBasis::Kinematic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        "rate",
        rate,
    )
    .map_err(|e| refuse(e.to_string()))
}
fn controls() -> ThermalVolumeMomentumControls {
    let mut c = momentum_controls();
    c.mechanics.projection.duration_s = 0.5;
    c.volume.thermal.volume.transport.projection = c.mechanics.projection;
    c
}
fn viscous_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 100,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-9,
        maximum_angular_impulse_defect_kg_m2_s: 1e-9,
        maximum_energy_defect_j: 1e-8,
    }
}
fn input<'a>(
    mass: &'a [f64],
    phases: &'a [Option<EquilibriumPhaseState>],
    momentum: &'a [[f64; 3]],
    heat: &'a [f64],
    forces: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: phases,
            net_enthalpy_increment_j: heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: forces,
        boundary_velocity_m_s: [None; 6],
    }
}
fn mouths() -> ([FlowFace; 6], [FvBoundary; 6]) {
    let mut faces = [FlowFace::Wall; 6];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    faces[..2].fill(FlowFace::Free);
    boundaries[..2].fill(FvBoundary::Outlet);
    boundaries[2..4].fill(FvBoundary::wall());
    (faces, boundaries)
}
fn pressure_run(
    field: &FlowField,
    domain: &VoxelDomain,
    curve: &EquilibriumEnthalpyPhaseCurve,
    input: ThermalVolumeMomentumInput<'_>,
    gate: &CancelGate,
    resolve: impl FnMut(
        usize,
        &EquilibriumPhaseState,
        f64,
    ) -> Result<ResolvedPhaseTransportState, ChtError>,
) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
    let x = [Some(1e12); 2];
    let sides = [None; 2];
    let (faces, boundaries) = mouths();
    field.advance_piecewise_thermal_volume_momentum_with_source_rheology_and_pressure_profile(
        domain,
        faces,
        StaggeredBoundaryPressure {
            x_pa: &x,
            y_pa: &sides,
            z_pa: &sides,
        },
        curve,
        input,
        boundaries,
        controls(),
        viscous_controls(),
        gate,
        resolve,
    )
}

#[test]
fn original_native_strain_and_current_density_resolve_one_kinematic_source() {
    let (card, melting) = fixture(true, true);
    let curve = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let phases = [Some(state)];
    let mass = [2.0];
    let momentum = [[2.0, 0.0, 0.0]];
    let mut input = input(&mass, &phases, &momentum, &[0.0], &[[0.0; 3]]);
    input.volume.boundary_inflow[0] = Some(state);
    input.boundary_velocity_m_s[0] = Some([0.5, 0.0, 0.0]);
    let mut calls = Vec::new();
    let result = pressure_run(
        &field,
        &domain,
        curve,
        input,
        &CancelGate::new_clock_free(),
        |cell, phase, rate| {
            calls.push((cell, *phase, rate));
            resolve(&card, &law, phase.specific_enthalpy_j_kg(), rate)
        },
    )
    .unwrap();
    assert_eq!(calls.len(), 1);
    assert_eq!((calls[0].0, calls[0].1), (0, state));
    close(calls[0].2, 2.0);
    close(result.original_shear_rate_per_s[0], 2.0);
    let source = result.original_sources[0].as_ref().unwrap();
    assert_eq!(source.phase_state(), &state);
    let transport = source.transport_state().unwrap();
    close(transport.density_kg_m3(), 2.0);
    close(transport.dynamic_viscosity_pa_s(), 1.0);
    let property = source.material().property(NU_KEY).unwrap();
    close(property.value_si(), 0.5);
    for property in source.material().properties() {
        let point = &property.answer().receipt.query_point;
        assert!(point.contains(&("T".into(), 300.0)));
        assert!(point.contains(&("rate".into(), 2.0)));
    }
    // Unit native stress has gamma²=4. The resolved mu=1 gives u*=1/2,
    // physical heat=1/2 and separate BE loss=1/4 over dt=1/2.
    close(result.evolution.viscosity.molecular_heat_j, 0.5);
    close(result.evolution.viscosity.backward_euler_loss_j, 0.25);
    let volume = &result.evolution.evolution.thermal.volume;
    assert!(volume.volume_source_m3_s[0] > 0.0);
    let after = volume.phase_states[0].unwrap();
    assert!(after.bulk_density_kg_m3() < 2.0);
    let mass = volume.transport.mass_kg[0];
    close(mass / after.bulk_density_kg_m3(), 1.0);
}

#[test]
fn zero_shear_and_constant_sources_reduce_to_the_original_closed_owner() {
    let (card, melting) = fixture(false, false);
    let curve = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phases = [Some(curve.state_at_specific_enthalpy(0.5).unwrap()); 2];
    let mass = [2.0; 2];
    let boundaries = [FvBoundary::wall(); 6];
    for momentum in [[[0.0; 3]; 2], [[2.0, 0.0, 0.0], [0.0; 3]]] {
        let input = input(&mass, &phases, &momentum, &[0.0; 2], &[[0.0; 3]; 2]);
        let mut calls = [0; 2];
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_source_rheology(
                &domain,
                curve,
                input,
                boundaries,
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
                |cell, phase, rate| {
                    calls[cell] += 1;
                    resolve(&card, &law, phase.specific_enthalpy_j_kg(), rate)
                },
            )
            .unwrap();
        assert_eq!(calls, [1; 2]);
        if momentum == [[0.0; 3]; 2] {
            assert_eq!(result.original_shear_rate_per_s, [0.0; 2]);
        }
        let sources: Vec<Option<&dyn FluidTransport>> = result
            .original_sources
            .iter()
            .map(|source| source.as_ref().map(|s| s as &dyn FluidTransport))
            .collect();
        let raw = field
            .advance_piecewise_thermal_volume_momentum_with_viscosity(
                &domain,
                curve,
                input,
                MaterialMomentumViscosityInput {
                    sources: &sources,
                    boundaries,
                },
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.evolution, raw);
    }
}

#[test]
fn source_mismatch_phase_exit_and_cancellation_preserve_input_and_retry() {
    let (card, melting) = fixture(false, true);
    let curve = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let phases = [Some(state)];
    let mass = [2.0];
    let momentum = [[2.0, 0.0, 0.0]];
    let mut input = input(&mass, &phases, &momentum, &[0.0], &[[0.0; 3]]);
    input.volume.boundary_inflow[0] = Some(state);
    input.boundary_velocity_m_s[0] = Some([0.5, 0.0, 0.0]);
    let run = |input, wrong: bool, gate: &CancelGate, calls: &mut usize| {
        pressure_run(&field, &domain, curve, input, gate, |_, phase, rate| {
            *calls += 1;
            let h = phase.specific_enthalpy_j_kg() + if wrong { 0.1 } else { 0.0 };
            resolve(&card, &law, h, rate)
        })
    };
    let gate = CancelGate::new_clock_free();
    let expected = run(input, false, &gate, &mut 0).unwrap();
    assert!(run(input, true, &gate, &mut 0).is_err());
    let mut mixed_inflow = input;
    mixed_inflow.volume.boundary_inflow[0] = Some(curve.state_at_specific_enthalpy(-0.5).unwrap());
    let mut incoming_calls = 0;
    assert!(
        run(mixed_inflow, false, &gate, &mut incoming_calls)
            .unwrap_err()
            .to_string()
            .contains("fully liquid incoming material")
    );
    assert_eq!(incoming_calls, 0);
    let mut cooling = input;
    cooling.volume.net_enthalpy_increment_j = &[-2.0];
    assert!(run(cooling, false, &gate, &mut 0).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let mut calls = 0;
    assert!(matches!(
        run(input, false, &cancelled, &mut calls),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(calls, 0);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(phases, [Some(state)]);
    assert_eq!(mass, [2.0]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]]);
    assert_eq!(run(input, false, &gate, &mut 0).unwrap(), expected);
}
