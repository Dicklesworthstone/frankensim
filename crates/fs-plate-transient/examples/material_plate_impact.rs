//! Material-backed two-way impact on a spatial DKT plate.
//!
//! Run: cargo run -p fs-plate-transient --example material_plate_impact
//! Both cases are synthetic fixtures, not named real materials or calibration.
//! The plate's E, nu and density come from the same resolved surface-B bulk
//! state retained by the contact source. Contact K and chi come from that source.
//! The declared half-space depth and layer thickness equal the plate thickness.
//! The fixed quadrature transfers the resultant and samples mean plate motion;
//! it does not resolve the changing Hertz pressure distribution.
//! Physical energy and instantaneous dissipation are reported without claiming
//! exact generalized-alpha step energy conservation.

use fs_contact::interface_binding::{
    ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, NORMAL_CONTACT_LAW_VERSION,
    NORMAL_CONTACT_STATE_SCHEMA_VERSION, NORMAL_HUNT_CROSSLEY_POINT_LAW_ID,
    NormalContactModelSelection, bind_dry_interface_system_card, bind_isotropic_elastic_interface,
    bind_normal_contact_model, bind_normal_interface_state,
};
use fs_contact::normal_patch::{
    InputUncertainty,
    source::{PointContactEnvelope, PointContactGeometry, SourcePointContact},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterfaceSystemCard, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateMesh, PlateModel, PlateSection, assemble,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact,
    PlateImpactBudget,
};
use fs_qty::{Density, Dims, Pressure};
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState,
    initialization::{InitialSolveConfig, initial_workspace_components, second_order_acceleration},
    second_order_adjoint::trajectory::StructuralTrajectoryModel,
};
use fs_tribo::InputAuthority;
use std::{collections::BTreeMap, error::Error};

/// Actual section thickness and source depth/layer dimensions.
pub const THICKNESS_M: f64 = 0.008;
/// Initial clearance before contact closes.
pub const GAP_M: f64 = 1.0e-6;
/// Finite normal translating mass, supplied by the geometry owner.
pub const BODY_MASS_KG: f64 = 0.001;
/// Fixed integration step, inside the original source envelope.
pub const STEP_S: f64 = 5.0e-5;

fn provenance(label: &str) -> Provenance {
    Provenance {
        source: format!("synthetic material-plate example: {label}"),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance(name),
    }
}
fn solid(name: &str, young_pa: f64, density: f64) -> Result<MaterialCard, Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        ("density", density, Density::DIMS),
        ("young_modulus", young_pa, Pressure::DIMS),
        ("poisson_ratio", 0.25, Dims::NONE),
        ("yield_stress", 1.0e10, Pressure::DIMS),
    ] {
        claims.insert_claim(claim(key, value, dims))?;
    }
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic example".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

/// Resolve explicit synthetic bulk and interface cards through their real owners.
pub fn synthetic_source(
    plate_name: &str,
    young_pa: f64,
    density: f64,
    chi: f64,
) -> Result<SourcePointContact, Box<dyn Error>> {
    let body = solid("synthetic indenter", 1.0e8, 1000.0)?;
    let plate = solid(plate_name, young_pa, density)?;
    let mut parameters = BTreeMap::new();
    for (name, value, dims) in [
        ("characteristic-rate", 1.0, Dims([1, 0, -1, 0, 0, 0])),
        ("max-patch-to-radius", 0.1, Dims::NONE),
        ("max-strain", 0.01, Dims::NONE),
        ("max-patch-to-depth", 0.1, Dims::NONE),
        ("max-patch-to-layer", 0.1, Dims::NONE),
        ("max-pressure-to-yield", 0.2, Dims::NONE),
        ("max-rate-ratio", 1.0, Dims::NONE),
        ("dissipation", chi, Dims([-1, 0, 1, 0, 0, 0])),
    ] {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    let model = ConstitutiveModelCard {
        law: LawId(NORMAL_HUNT_CROSSLEY_POINT_LAW_ID.into()),
        law_version: NORMAL_CONTACT_LAW_VERSION,
        parameters,
        state_schema_version: NORMAL_CONTACT_STATE_SCHEMA_VERSION,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        sources: vec![],
        provenance: provenance("Hunt-Crossley interface"),
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim(ADHESION_ENERGY_PROPERTY, 0.0, ADHESION_ENERGY_DIMS))?;
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: body.id().clone(),
            texture_frame: "body surface".into(),
        },
        SurfaceSpec {
            material: plate.id().clone(),
            texture_frame: "plate surface".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic example".into(),
            history: "initial".into(),
        },
        claims,
        vec![model],
    )?;
    let query = QueryPoint::new().with("T", 293.15)?;
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let body_state = resolve_isotropic_solid_state_point(&body, &query, selection.clone())?;
    let plate_state = resolve_isotropic_solid_state_point(&plate, &query, selection.clone())?;
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture)?;
    let elastic = bind_isotropic_elastic_interface(&dry, &body_state, &plate_state)?;
    let interface = resolve_interface_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?],
        selection,
    )?;
    let normal = bind_normal_interface_state(&elastic, &interface)?;
    let model = bind_normal_contact_model(&normal, NormalContactModelSelection::SingleSupported)?;
    Ok(SourcePointContact::new(
        model,
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.01,
        },
        PointContactEnvelope {
            maximum_indentation_m: 2.0e-5,
            maximum_absolute_rate_m_per_s: 0.1,
            maximum_step_s: 1.0e-4,
            half_space_depth_m: THICKNESS_M,
            layer_thickness_m: THICKNESS_M,
            uncertainty: InputUncertainty {
                radius_relative: 0.01,
                modulus_relative: 0.02,
                load_relative: 0.03,
            },
        },
    )?)
}

/// Assemble every spatial DOF from the exact retained surface-B bulk state.
pub fn build_plate(source: &SourcePointContact) -> Result<(PlateMesh, PlateModel), Box<dyn Error>> {
    let bulk = source.model().normal_state().elastic().surface_b_state();
    let section = PlateSection::isotropic(
        bulk.young_modulus_pa(),
        bulk.poisson_ratio(),
        THICKNESS_M,
        bulk.density_kg_m3(),
    )?;
    let mesh = PlateMesh::rectangle(0.4, 0.3, 3, 3);
    let plate = assemble(
        &mesh,
        &section,
        &PlateMesh::rectangle_boundary(3, 3),
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )?;
    Ok((mesh, plate))
}

/// Bind source and spatial plate through one fixed total-force transfer footprint.
pub fn build_impact<'a>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    source: SourcePointContact,
    load: &'a NoPlateLoad,
) -> Result<PlateImpact<'a, NoPlateLoad>, Box<dyn Error>> {
    let points = [
        PlatePatchPoint {
            point: [0.19985, 0.14985],
            area_m2: 9.0e-8,
        },
        PlatePatchPoint {
            point: [0.20015, 0.14985],
            area_m2: 9.0e-8,
        },
        PlatePatchPoint {
            point: [0.19985, 0.15015],
            area_m2: 9.0e-8,
        },
        PlatePatchPoint {
            point: [0.20015, 0.15015],
            area_m2: 9.0e-8,
        },
    ];
    let patch = PlatePatchStencil::locate(
        mesh,
        plate,
        &points,
        PlatePatchBudget {
            max_points: points.len(),
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    )?;
    let dynamics = PlateDynamics::new(
        plate,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: 1.0,
            mass_damping_per_s: 0.0,
            stiffness_damping_s: 0.0,
        },
        load,
        PlateDynamicsBudget {
            max_dofs: plate.free,
            max_nonzeros: plate.k.nnz() + plate.m.nnz(),
            max_full_dofs: plate.dof_map.len(),
            max_load_parameters: 0,
        },
        &mut || false,
    )?;
    Ok(PlateImpact::from_source(
        dynamics,
        patch,
        BODY_MASS_KG,
        source,
        GAP_M,
        PlateImpactBudget {
            max_free_dofs: plate.free + 1,
            max_contact_points: 1,
        },
    )?)
}

/// Bounded generalized-alpha solver used by the experiment.
pub fn method(n: usize, step_s: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        step_s,
        1.0,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1.0e-12,
                relative_tolerance: 1.0e-11,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_minimum: 1.0e-13,
                forcing_maximum: 1.0e-11,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 16,
        },
    )
}

/// Production consistent initialization; the indenter approaches at 5 mm/s.
pub fn initial_state(
    impact: &PlateImpact<'_, NoPlateLoad>,
) -> Result<SecondOrderState, Box<dyn Error>> {
    let n = impact.dimension();
    let q = vec![0.0; n];
    let mut v = vec![0.0; n];
    v[impact.plate_dofs()] = 0.005;
    let mut force = vec![0.0; n];
    impact.forcing(0.0, &mut force)?;
    let config = InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1.0e-12,
    };
    let a = second_order_acceleration(
        impact,
        &q,
        &v,
        &force,
        &impact.mass_preconditioner(),
        config,
        initial_workspace_components(n, 0, config).ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    Ok(SecondOrderState::new(0.0, &q, &v, &a.value))
}

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "case,young_pa,density_kg_m3,contact_k_n_per_m_3_2,chi_s_per_m,peak_force_n,peak_plate_displacement_m,final_body_velocity_m_per_s,final_energy_j"
    );
    for (name, young, density, chi) in [
        ("synthetic_A", 1.0e7, 900.0, 0.4),
        ("synthetic_B", 3.0e7, 1800.0, 0.8),
    ] {
        let source = synthetic_source(name, young, density, chi)?;
        let k = source.stiffness_n_per_m_3_2();
        let selected_chi = source.dissipation_s_per_m();
        let (mesh, plate) = build_plate(&source)?;
        let no_load = NoPlateLoad;
        let impact = build_impact(&mesh, &plate, source, &no_load)?;
        let mut state = initial_state(&impact)?;
        let integrator = method(impact.dimension(), STEP_S);
        let mut force = vec![0.0; impact.dimension()];
        let mut peak_force: f64 = 0.0;
        let mut peak_displacement: f64 = 0.0;
        for _ in 0..100 {
            impact.forcing(integrator.forcing_time(state.t)?, &mut force)?;
            integrator.step(&mut state, &impact, &force)?;
            let sample = impact.probe(&state.q, &state.v)?;
            peak_force = peak_force.max(sample.normal_force_n);
            peak_displacement = peak_displacement.max(sample.patch.displacement.abs());
        }
        let sample = impact.probe(&state.q, &state.v)?;
        println!(
            "{name},{young:.8e},{density:.8e},{k:.8e},{selected_chi:.8e},{peak_force:.8e},{peak_displacement:.8e},{:.8e},{:.8e}",
            sample.indenter_velocity_m_per_s, sample.total_energy_j,
        );
    }
    Ok(())
}
