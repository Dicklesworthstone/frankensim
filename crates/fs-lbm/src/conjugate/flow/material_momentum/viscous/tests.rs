//! Source viscosity, actual inertia, molecular heat and atomic split acceptance.

use super::*;
use crate::conjugate::flow::{
    DensityForceControls, DensityMomentumReconstructionControls, DensityProjectionControls,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, ResolvedLiquidState, resolve_liquid_state,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

fn controls() -> MaterialMomentumAdvanceControls {
    let projection = DensityProjectionControls {
        duration_s: 0.5,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-7,
    };
    MaterialMomentumAdvanceControls {
        mechanics: MaterialDensityMomentumControls {
            reconstruction: DensityMomentumReconstructionControls {
                max_cells: 100,
                maximum_mass_defect_kg: 1e-10,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            forcing: DensityForceControls {
                max_cells: 100,
                maximum_impulse_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            projection,
            maximum_feedback_momentum_defect_kg_m_s: 1e-8,
            maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-8,
            maximum_feedback_energy_defect_j: 1e-7,
        },
        material: MaterialAdvectionControls {
            projection,
            max_cells: 100,
            max_connections: 1000,
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-9,
            maximum_enthalpy_defect_j: 1e-7,
        },
        momentum: MaterialMomentumControls {
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_momentum_defect_kg_m_s: 1e-8,
            maximum_energy_defect_j: 1e-7,
        },
    }
}

fn viscous_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 200,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-8,
        maximum_angular_impulse_defect_kg_m2_s: 1e-8,
        maximum_energy_defect_j: 1e-7,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

/// Synthetic source fixture, not measured liquid data. Source rho is
/// deliberately unrelated to the actual evolving M/V used by mechanics.
fn liquid(mu: f64) -> ResolvedLiquidState {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1000.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            1.0,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic viscous momentum fixture, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_liquid_state(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn input<'a>(
    mass: &'a [f64; 2],
    h: &'a [f64; 2],
    p: &'a [[f64; 3]; 2],
) -> MaterialMomentumAdvanceInput<'a> {
    MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: mass,
            specific_enthalpy_j_kg: h,
            volume_sources: &[MaterialVolumeSource::None; 2],
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: p,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 2],
        source_velocity_m_s: &[None; 2],
        boundary_velocity_m_s: [None; 6],
    }
}

#[test]
fn sourced_viscosity_uses_actual_graded_inertia_and_heats_only_molecular_loss() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0]; // actual rho=[2,6], unlike source rho=1000
    let h = [10.0, 30.0];
    let momentum = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let slow = liquid(1.0);
    let fast = liquid(3.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&slow), Some(&fast)];
    let run = |sources| {
        field
            .advance_material_momentum_with_viscosity(
                &domain,
                input(&mass, &h, &momentum),
                MaterialMomentumViscosityInput {
                    sources,
                    boundaries: [FvBoundary::Symmetry; 6],
                },
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let next = run(&sources);
    // One active x face: M_face=(4+36)/2=20 kg; normal-stress
    // stiffness=2*A*(mu0/dx0+mu1/dx1)=8 kg/s.
    // Thus BE u*=20/(20+0.5*8)=5/6, independently of source rho.
    let u = 5.0 / 6.0;
    let impulse = 20.0 * (u - 1.0);
    close(next.viscosity.face_impulse_kg_m_s[0], impulse);
    close(next.viscosity.cell_momentum_change_kg_m_s[0], impulse);
    close(
        next.viscosity.face_angular_impulse_kg_m2_s[1],
        0.5 * impulse,
    );
    close(next.viscosity.face_angular_impulse_kg_m2_s[2], -impulse);
    close(
        next.viscosity.cell_angular_momentum_change_kg_m2_s[1],
        0.5 * impulse,
    );
    close(
        next.viscosity.cell_angular_momentum_change_kg_m2_s[2],
        -impulse,
    );
    close(next.viscosity.native_kinetic_energy_change_j, -55.0 / 18.0);
    close(next.viscosity.backward_euler_loss_j, 5.0 / 18.0);
    close(next.viscosity.numerical_scatter_loss_j, 5.0 / 36.0);
    close(next.viscosity.cell_kinetic_energy_change_j, -115.0 / 36.0);
    close(next.viscosity.representation_work_j, 0.0);
    close(next.viscosity.wall_work_j, 0.0);
    close(next.viscosity.energy_defect_j, 0.0);
    for c in 0..2 {
        // dt*2*mu*V*(u*/dx)^2, not total kinetic decay.
        close(next.viscosity.cell_molecular_heat_j[c], 25.0 / 18.0);
        close(next.advance.transport.mass_kg[c], mass[c]);
        close(
            next.advance.transport.specific_enthalpy_j_kg[c],
            h[c] + (25.0 / 18.0) / mass[c],
        );
        close(next.advance.momentum.momentum_kg_m_s[c][0], 0.0);
    }
    close(next.represented_molecular_heat_j, 25.0 / 9.0);
    assert!(next.viscosity.backward_euler_loss_j + next.viscosity.numerical_scatter_loss_j > 0.4);
    assert_eq!(
        next.advance.momentum.report.represented_impulse_kg_m_s,
        [0.0; 3]
    );
    // The original closed pressure owner removes its current x throughflow
    // after viscosity; no force, pressure, or heat was replayed by transport.
    close(next.advance.transport.projection.field.flux_x(1, 0, 0), 0.0);

    // Changing an original liquid source changes stress, without changing
    // actual inventory or substituting reference-density face masses.
    let changed_sources: [Option<&dyn FluidTransport>; 2] = [Some(&fast), Some(&fast)];
    let changed = run(&changed_sources);
    let changed_u = 20.0 / (20.0 + 0.5 * 16.0);
    close(
        changed.viscosity.face_impulse_kg_m_s[0],
        20.0 * (changed_u - 1.0),
    );
    assert!(changed.viscosity.face_impulse_kg_m_s[0] < impulse);
    assert_eq!(mass, [4.0, 36.0]);
    assert_eq!(h, [10.0, 30.0]);
}

#[test]
fn moving_wall_work_drives_actual_inertia_and_has_its_own_physical_heat() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 6.0];
    let h = [0.0; 2];
    let p = [[0.0; 3]; 2];
    let source = FluidTransportState::try_new(999.0, 2.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    boundaries[3] = FvBoundary::Wall {
        velocity: [1.0, 0.0, 0.0],
    };
    let result = field
        .advance_material_momentum_with_viscosity(
            &domain,
            input(&mass, &h, &p),
            MaterialMomentumViscosityInput {
                sources: &sources,
                boundaries,
            },
            controls(),
            viscous_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // One native DOF: M=4, normal stiffness=4, wall stiffness=2.
    // 4*u/0.5+6*u=2 gives u=1/7. Physical wall work=6/7 J,
    // physical heat=38/49 J, native storage and BE loss each=2/49 J.
    let report = &result.viscosity;
    close(report.face_impulse_kg_m_s[0], 4.0 / 7.0);
    close(report.wall_work_by_face_j[3], 6.0 / 7.0);
    close(report.wall_work_j, 6.0 / 7.0);
    close(report.molecular_heat_j, 38.0 / 49.0);
    close(report.native_kinetic_energy_change_j, 2.0 / 49.0);
    close(report.backward_euler_loss_j, 2.0 / 49.0);
    close(report.cell_kinetic_energy_change_j, 1.0 / 49.0);
    close(report.numerical_scatter_loss_j, 1.0 / 49.0);
    close(report.energy_defect_j, 0.0);
    for c in 0..2 {
        close(report.cell_molecular_heat_j[c], 19.0 / 49.0);
        close(
            result.advance.transport.specific_enthalpy_j_kg[c],
            (19.0 / 49.0) / mass[c],
        );
    }
    close(result.represented_molecular_heat_j, 38.0 / 49.0);
}

#[test]
fn repeated_step_reconstructs_current_momentum_and_quiescence_preserves_old_advance() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 6.0];
    let h = [10.0, 30.0];
    let zero = [[0.0; 3]; 2];
    let source = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let material = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Symmetry; 6],
    };
    let still = field
        .advance_material_momentum_with_viscosity(
            &domain,
            input(&mass, &h, &zero),
            material,
            controls(),
            viscous_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let original = field
        .advance_material_momentum(
            &domain,
            [FlowFace::Wall; 6],
            input(&mass, &h, &zero),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(still.advance, original);
    assert_eq!(still.viscosity.molecular_heat_j, 0.0);
    let unrelated = FlowField::from_face_velocity(&domain, [FlowFace::Wall; 6], |axis, _| {
        if axis == 0 { 17.0 } else { 0.0 }
    })
    .unwrap();
    let repeated = unrelated
        .advance_material_momentum_with_viscosity(
            &domain,
            MaterialMomentumAdvanceInput {
                material: MaterialAdvectionInput {
                    mass_kg: &still.advance.transport.mass_kg,
                    specific_enthalpy_j_kg: &still.advance.transport.specific_enthalpy_j_kg,
                    ..input(&mass, &h, &zero).material
                },
                momentum_kg_m_s: &still.advance.momentum.momentum_kg_m_s,
                ..input(&mass, &h, &zero)
            },
            material,
            controls(),
            viscous_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(repeated, still);
}

#[test]
fn source_scope_work_cancellation_and_unrepresentable_heat_refuse_atomically() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original_field = field.clone();
    let mass = [2.0, 6.0];
    let h = [10.0, 30.0];
    let p = [[2.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    let source = liquid(1.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let material = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Symmetry; 6],
    };
    let run = |data, viscosity, controls, gate: &CancelGate| {
        field.advance_material_momentum_with_viscosity(
            &domain,
            data,
            viscosity,
            self::controls(),
            controls,
            gate,
        )
    };
    let accepted = run(
        input(&mass, &h, &p),
        material,
        viscous_controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let mut no_cells = viscous_controls();
    no_cells.max_cells = 1;
    let mut no_sweeps = viscous_controls();
    no_sweeps.max_sweeps = 0;
    let mut one_sweep = viscous_controls();
    one_sweep.max_sweeps = 1;
    for limited in [no_cells, no_sweeps, one_sweep] {
        assert!(
            run(
                input(&mass, &h, &p),
                material,
                limited,
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
    let missing: [Option<&dyn FluidTransport>; 2] = [Some(&source), None];
    assert!(
        run(
            input(&mass, &h, &p),
            MaterialMomentumViscosityInput {
                sources: &missing,
                ..material
            },
            viscous_controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut open = material;
    open.boundaries[0] = FvBoundary::Outlet;
    assert!(
        run(
            input(&mass, &h, &p),
            open,
            viscous_controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    // Native stress/heat succeeds but h cannot represent the positive heat.
    // No enthalpy or momentum candidate escapes this late physical gate.
    let huge_h = [1e30; 2];
    assert!(
        run(
            input(&mass, &huge_h, &p),
            material,
            viscous_controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        run(input(&mass, &h, &p), material, viscous_controls(), &gate),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, original_field);
    assert_eq!(mass, [2.0, 6.0]);
    assert_eq!(h, [10.0, 30.0]);
    assert_eq!(p, [[2.0, 0.0, 0.0], [6.0, 0.0, 0.0]]);
    assert_eq!(
        run(
            input(&mass, &h, &p),
            material,
            viscous_controls(),
            &CancelGate::new_clock_free()
        )
        .unwrap(),
        accepted
    );
}
