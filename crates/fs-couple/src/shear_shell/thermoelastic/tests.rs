//! G1/G3 physical shell strain, constitutive response and exact work/heat checks.
//! All source data are synthetic; these are discrete-model, not experimental tests.
use super::{ShellThermalProperties, ShellThermoelasticOptions, source};
use super::source::tests::{fixture, fixture_with};
use crate::shear_shell::{
    MaterialShellImage, ResolvedShearShell, ShellProjectionBudget, ShellReductionOptions,
    compile_material_shear_shell,
};
use crate::shear_shell::connections::{
    MaterialShellConnection, ShellInterfacePoint, ShellInterfaceSide,
    stiffness_per_area_dims, viscosity_per_area_dims,
};
use crate::thin_plate::PlateMaterialModel;
use crate::modal_acoustic_time::ModalAcousticTimeBudget;
use crate::render::schedule::force::coupled::{
    ModalCouplingConfig, PreparedModalBudget,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point,
};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;
use fs_math::{c64::C64, det};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
    connections::ShellConnectionBudget, strain::{PreparedElementStrain, integration_rule},
};

fn options() -> ShellThermoelasticOptions {
    ShellThermoelasticOptions {
        band_rad_s: (200.0, 20_000.0), normalized_tolerance: 1e-4,
        max_thermal_terms: 16, max_memory_states: 64,
        max_projection_terms: 100_000_000, max_factor_entries: 20_000,
        reference_stiffness_n_m: 1e4,
    }
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8, max_connections: 64, max_setup_terms: 10_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5, maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    }
}
fn preparation() -> PreparedModalBudget {
    PreparedModalBudget { max_modes: 8, max_dense_bytes: 128*1024*1024,
        max_multiply_terms: 20_000_000_000 }
}

fn with_interface(source: ResolvedShearShell, band: (f64,f64), gate: &CancelGate) -> ResolvedShearShell {
    let material = &source.regions()[0].input().region.material;
    let left = material.material().clone();
    let mut right = left.clone();
    right.process = "synthetic fixed counterface".into();
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("spring", stiffness_per_area_dims(), 1e7),
        ("dashpot", viscosity_per_area_dims(), 500.0)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("omega", band.0, band.1),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "synthetic independent shell attachment".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::NonNegative).unwrap());
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec { material: left, texture_frame: "shell point".into() },
        SurfaceSpec { material: right.clone(), texture_frame: "fixed point".into() },
        SystemContext { medium: "synthetic".into(), third_body: None,
            environment: "held".into(), history: "linear bilateral".into() },
        claims, vec![],
    ).unwrap();
    let point = material.query().with("omega", 0.5*(band.0+band.1)).unwrap();
    let interface = resolve_interface_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    source.with_interfaces(vec![MaterialShellConnection {
        name: "independent support".into(), interface,
        stiffness_property: "spring".into(), viscosity_property: "dashpot".into(),
        area_m2: 1e-4, direction: [0.0, 0.0, 1.0],
        left: ShellInterfaceSide::Region { name: "source".into(),
            points: vec![ShellInterfacePoint { node: 8, weight: 1.0, offset_m: [0.0; 3] }] },
        right: ShellInterfaceSide::Fixed { material: right },
    }], ShellConnectionBudget { max_connections: 1, max_terms: 1, max_work: 100_000 }, gate).unwrap()
}

fn reduced(curved: bool, zero_alpha: bool, independent_loss: bool)
    -> (MaterialShellImage, ShellThermalProperties)
{
    let gate = CancelGate::new_clock_free();
    let (base, _) = fixture(PlateMaterialModel::Orthotropic12);
    let point = base.regions()[0].input().region.material.query();
    let (mut source, properties) = fixture_with(PlateMaterialModel::Orthotropic12, point, |claims| {
        for p in claims {
            if zero_alpha && p.key.name().starts_with("alpha_") {
                p.value = PropertyValue::Scalar { value: 0.0,
                    dims: fs_material::state_point::INVERSE_TEMPERATURE_DIMS };
            }
            if independent_loss && PLANE_STRESS_VISCOSITY_PROPERTIES[..3].contains(&p.key.name()) {
                p.value = PropertyValue::Scalar { value: 1000.0, dims: fs_qty::DynViscosity::DIMS };
            }
        }
    });
    if curved {
        let mut mesh = source.mesh().clone();
        for p in &mut mesh.nodes { p[2] = 0.03*p[0]+0.08*p[0]*p[1]; }
        let inputs = source.regions().iter().map(|r| r.input().clone()).collect();
        source = compile_material_shear_shell(mesh, inputs,
            MitcBudget { max_nodes: 9, max_triangles: 8 }, &gate).unwrap();
    }
    if independent_loss { source = with_interface(source, (50.0,1e7), &gate); }
    (reduce_source(&source, curved, &gate), properties)
}

fn reduce_source(source: &ResolvedShearShell, curved: bool, gate: &CancelGate) -> MaterialShellImage {
    let normals = source.mesh().nodes.iter().map(|p| {
        if !curved { return [0.0, 0.0, 1.0]; }
        let v = [-0.03-0.08*p[1], -0.08*p[0], 1.0];
        let norm = det::sqrt(v.iter().map(|x| x*x).sum());
        v.map(|x| x/norm)
    }).collect();
    let directors = NodalDirectors::new(source.mesh(), normals, 0.9).unwrap();
    let reduction = ShellReductionOptions {
        constraints: [0, 3, 6].map(|node| ShellConstraint {
            node, mask: [true; 6], frame: ConstraintFrame::Global,
        }).to_vec(),
        mass: ShellMass::Consistent, drilling_alpha: 1e-3,
        eigenvalue_window: (200.0_f64.powi(2), 20_000.0_f64.powi(2)), n_modes: 3,
        assembly_budget: MitcBudget { max_nodes: 9, max_triangles: 8 },
        projection_budget: ShellProjectionBudget { max_modes: 3,
            max_projection_terms: 100_000_000, reference_resistance_n_s_m: 1.0 },
    };
    let mut ports = vec![vec![0.0; 54]; 2];
    ports[0][6*8+2] = 1.0;
    ports[1][6*2+2] = 1.0;
    source.reduce(&directors, &reduction, &ports, gate).unwrap()
}

/// Independent integral of actual original MITC strain polynomials. Transform
/// strain into the principal frame, rather than rotating the production factor.
fn original_thermal_grams(original: &MaterialShellImage, properties: &[ShellThermalProperties],
    region_filter: Option<usize>)
    -> (Vec<f64>, Vec<f64>)
{
    let bound = source::bind_regions(original, properties,
        options().band_rad_s, &CancelGate::new_clock_free()).unwrap();
    let n = original.physical_modes().len();
    let mut membrane = vec![0.0; n*n];
    let mut bending = vec![0.0; n*n];
    let mesh = original.source().mesh();
    let rule = integration_rule();
    for (element, tri) in mesh.tris.iter().enumerate() {
        let index = original.source().element_regions[element];
        if region_filter.is_some_and(|r| r!=index) { continue; }
        let region = &original.source().regions()[index];
        let f = bound[index].law.relaxation_factor();
        let scale = det::sqrt(12.0)/region.thickness_m();
        let facet = mesh.facet(element).unwrap();
        let angle = fs_plate::shell::survey::material_axis_angle(
            &facet, region.input().axis1_global).unwrap()+region.input().region.material_angle_rad;
        let (s,c) = (det::sin(angle), det::cos(angle));
        let contract = |e: [f64; 3]| {
            let local = [c*c*e[0]+s*s*e[1]+c*s*e[2],
                s*s*e[0]+c*c*e[1]-c*s*e[2],
                2.0*c*s*(e[1]-e[0])+(c*c-s*s)*e[2]];
            (0..3).map(|k| f[k]*local[k]).sum::<f64>()
        };
        let observer = PreparedElementStrain::new(mesh, element, original.directors(),
            &original.source().sections()[element], &mut || false).unwrap();
        let fields: Vec<_> = original.physical_modes().iter().map(|phi| {
            let displacement = core::array::from_fn(|k| phi[6*tri[k/6]+k%6]);
            observer.real_field(&displacement, &mut || false).unwrap()
        }).collect();
        for &(lambda, weight) in &rule {
            let strain: Vec<_> = fields.iter().map(|f| f.at_barycentric(lambda).unwrap()).collect();
            let a: Vec<_> = strain.iter().map(|e| scale*contract(e.membrane)).collect();
            let d: Vec<_> = strain.iter().map(|e| contract(e.gradient_per_m)).collect();
            for i in 0..n { for j in 0..n {
                membrane[n*i+j] += facet.area_m2*weight*a[i]*a[j];
                bending[n*i+j] += facet.area_m2*weight*d[i]*d[j];
            } }
        }
    }
    (membrane, bending)
}

fn max_complex(z: C64) -> f64 { z.re.abs().max(z.im.abs()) }
fn quadratic(matrix: &[f64], q: &[C64]) -> f64 {
    let n = q.len();
    (0..n).map(|i| (0..n).map(|j|
        matrix[n*i+j]*(q[i].re*q[j].re+q[i].im*q[j].im)).sum::<f64>()).sum()
}

#[test]
fn g1_curved_anisotropic_shell_matches_original_strain_and_complex_dynamic_stiffness() {
    let gate = CancelGate::new_clock_free();
    let (original, properties) = reduced(true, false, true);
    let image = original.with_thermoelastic_memory(std::slice::from_ref(&properties), options(), &gate).unwrap();
    assert_eq!(image.original().physical_modes(), original.physical_modes());
    assert_eq!(image.original().force_ports(), original.force_ports());
    assert_eq!(image.original().volume_shapes_m2_per_sqrt_kg(), original.volume_shapes_m2_per_sqrt_kg());
    assert_eq!(image.original().cavity_volume_shapes_m2_per_sqrt_kg(),
        original.cavity_volume_shapes_m2_per_sqrt_kg());
    let (membrane, bending) = original_thermal_grams(&original, std::slice::from_ref(&properties), None);
    let n = original.physical_modes().len();
    let weights: f64 = image.spectra()[0].terms().iter().map(|t| t.weight()).sum();
    for (expected, is_membrane) in [(&membrane, true), (&bending, false)] {
        let norm = expected.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(norm > 0.0);
        for i in 0..n { for j in 0..n {
            let actual: f64 = if is_membrane {
                image.membrane_columns().iter().map(|c|
                    image.reference_stiffness_n_m()*c[i]*c[j]).sum()
            } else {
                image.relaxation_branches().iter().map(|b|
                    b.stiffness*b.projection[2*i]*b.projection[2*j]).sum::<f64>()/weights
            };
            assert!((actual-expected[n*i+j]).abs() <= 5e-10*norm);
        } }
    }
    let off_diagonal = (0..n).flat_map(|i| (0..n).filter(move |&j| i!=j)
        .map(move |j| (i,j))).map(|(i,j)| bending[n*i+j].abs()).fold(0.0_f64, f64::max);
    assert!(off_diagonal > 1e-8*bending.iter().map(|v| v.abs()).fold(0.0_f64, f64::max));
    let system = image.prepare_coupled(96_000, ModalAcousticTimeBudget::audible_reference(),
        config(), preparation(), &gate).unwrap();
    let omega = 1300.0;
    let force = [C64::new(0.2,0.05), C64::new(-0.3,0.1), C64::new(0.4,-0.02)];
    let response = system.harmonic_response(omega, &force, 1_000_000, &gate).unwrap();
    let shape = image.spectra()[0].terms().iter().fold(C64::ZERO, |sum,t| {
        let x = omega*t.relaxation_time_s();
        sum+C64::new(x*x/(1.0+x*x),-x/(1.0+x*x)).scale(t.weight())
    });
    let mut damping = original.damping().matrix_per_s().to_vec();
    for (law, column) in original.source().interfaces().iter().zip(original.interface_columns()) {
        for i in 0..n { for j in 0..n {
            damping[n*i+j] += law.law().damping_n_s_m*column[i]*column[j];
        } }
    }
    let q = &response.displacement_m_sqrt_kg;
    for i in 0..n {
        let mut actual = C64::ZERO;
        for j in 0..n {
            let elastic = if i==j { original.radiators()[i].omega.powi(2)-omega*omega } else { 0.0 };
            let stiffness = C64::new(elastic+membrane[n*i+j], -omega*damping[n*i+j])
                + shape.scale(bending[n*i+j]);
            actual = actual+stiffness*q[j];
        }
        assert!(max_complex(actual-force[i]) < 2e-9);
    }
    let thermal_heat = 0.5*omega*(-shape.im)*quadratic(&bending,q);
    let viscous_heat = 0.5*omega*omega*quadratic(&damping,q);
    assert!(thermal_heat>0.0 && viscous_heat>0.0);
    assert!((response.material_dissipation_w/thermal_heat-1.0).abs() < 2e-9);
    assert!((response.connection_dissipation_w/viscous_heat-1.0).abs() < 2e-9);
    assert!(response.power_residual_w.abs() < 1e-10*response.external_power_w.abs()+1e-18);
    assert_eq!(system.samples_rendered(),0,"harmonic observation cannot advance histories");
}


#[test]
fn g3_regions_share_only_identical_kernels_and_preserve_every_owned_gram() {
    let gate = CancelGate::new_clock_free();
    let (original, property) = reduced(true,false,false);
    let single = original.with_thermoelastic_memory(
        std::slice::from_ref(&property),options(),&gate).unwrap();
    let mut inputs = vec![original.source().regions()[0].input().clone();2];
    for (r,input) in inputs.iter_mut().enumerate() {
        input.region.region.name = format!("region-{r}");
        input.region.region.triangle_indices = (4*r..4*r+4).collect();
    }
    let properties = vec![property.clone();2];
    let compile = |inputs| compile_material_shear_shell(original.source().mesh().clone(),
        inputs,MitcBudget {max_nodes:9,max_triangles:8},&gate).unwrap();
    let split = reduce_source(&compile(inputs.clone()),true,&gate);
    let shared = split.with_thermoelastic_memory(&properties,options(),&gate).unwrap();
    assert_eq!(split.physical_modes(),original.physical_modes());
    assert_eq!(shared.relaxation_branches().len(),single.relaxation_branches().len());
    let n=split.physical_modes().len();
    let (a,b)=(shared.membrane_columns(),single.membrane_columns());
    for i in 0..n { for j in 0..n {
        let left:f64=a.iter().map(|row|row[i]*row[j]).sum();
        let right:f64=b.iter().map(|row|row[i]*row[j]).sum();
        assert!((left-right).abs() <= 5e-11*left.abs().max(right.abs()).max(1e-20));
    }}
    for i in 0..n { for j in 0..n {
        let gram=|image:&super::ThermoelasticShellImage| image.relaxation_branches().iter()
            .map(|b|b.stiffness*b.projection[2*i]*b.projection[2*j]).sum::<f64>();
        let (a,b)=(gram(&shared),gram(&single));
        assert!((a-b).abs() <= 5e-11*a.abs().max(b.abs()).max(1e-20));
    }}
    let (changed,_) = fixture_with(PlateMaterialModel::Orthotropic12,
        original.source().regions()[0].input().region.material.query(),|claims| {
            let p=claims.iter_mut().find(|p|p.key.name()=="k_normal").unwrap();
            p.value=PropertyValue::Scalar {value:200.0,
                dims:fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS};
        });
    inputs[1].region.material=changed.regions()[0].input().region.material.clone();
    let split=reduce_source(&compile(inputs),true,&gate);
    let mut allowance=options();allowance.max_memory_states=128;
    let distinct=split.with_thermoelastic_memory(&properties,allowance,&gate).unwrap();
    assert_eq!(split.physical_modes(),original.physical_modes());
    let times:Vec<_>=distinct.spectra().iter().map(|s|s.terms().iter()
        .map(|t|t.relaxation_time_s().to_bits()).collect::<Vec<_>>()).collect();
    assert!(times[0].iter().all(|t|!times[1].contains(t)));
    assert!((distinct.spectra()[1].terms()[0].relaxation_time_s()
        /distinct.spectra()[0].terms()[0].relaxation_time_s()-0.5).abs()<1e-14);
    for r in 0..2 {
        let (_,gram)=original_thermal_grams(&split,&properties,Some(r));
        let norm=gram.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
        for term in distinct.spectra()[r].terms() {
            let arms:Vec<_>=distinct.relaxation_branches().iter().filter(|b|
                b.relaxation_time_s.to_bits()==term.relaxation_time_s().to_bits()).collect();
            assert!(!arms.is_empty());
            for i in 0..n { for j in 0..n {
                let actual:f64=arms.iter().map(|b|
                    b.stiffness*b.projection[2*i]*b.projection[2*j]).sum();
                assert!((actual-term.weight()*gram[n*i+j]).abs()<=5e-10*term.weight()*norm);
            }}
        }
    }
    assert_eq!(distinct.relaxation_branches().iter().filter(|b|
        times[0].contains(&b.relaxation_time_s.to_bits())
        ||times[1].contains(&b.relaxation_time_s.to_bits())).count(),
        distinct.relaxation_branches().len());
}

#[test]
fn g1_pulse_decay_closes_work_heat_and_cancel_preserves_every_state() {
    let gate = CancelGate::new_clock_free();
    let (original, properties) = reduced(true, false, true);
    let image = original.with_thermoelastic_memory(&[properties], options(), &gate).unwrap();
    let mut system = image.prepare_coupled(96_000, ModalAcousticTimeBudget::audible_reference(),
        config(), preparation(), &gate).unwrap();
    assert_eq!(system.material_memory().len(),image.relaxation_branches().len());
    let mut forces = [0.0; 3];
    let (mut work, mut heat, mut material_heat, mut viscous_heat) = (0.0,0.0,0.0,0.0);
    for sample in 0..256 {
        image.fill_force(&[if sample<24 {0.1} else {0.0},0.0], &mut forces).unwrap();
        let step = system.step_under_gate(&forces,&gate).unwrap();
        assert!(step.material_dissipation_j>=0.0 && step.connection_dissipation_j>=0.0);
        assert_eq!(step.numerical_dissipation_j,0.0);
        work += step.external_work_j;
        heat += step.total_dissipation_j();
        material_heat += step.material_dissipation_j;
        viscous_heat += step.connection_dissipation_j;
    }
    let stored = system.total_energy_j().unwrap();
    assert!(material_heat>0.0 && viscous_heat>0.0 && stored>0.0);
    assert!((work-stored-heat).abs() < 2e-10*work.abs().max(stored)+1e-15);
    assert!(system.material_memory().iter().any(|z| *z!=0.0));
    for _ in 0..128 { system.step(&[0.0;3]).unwrap(); }
    assert!(system.total_energy_j().unwrap()<stored);
    let before_states = system.components()[0].states().to_vec();
    let before_memory = system.material_memory().to_vec();
    let before_frame = system.last_frame().cloned();
    let stop = CancelGate::new_clock_free(); stop.request();
    assert!(system.step_under_gate(&[0.1;3],&stop).is_err());
    assert_eq!(system.components()[0].states(),before_states);
    assert_eq!(system.material_memory(),before_memory);
    assert_eq!(system.last_frame(),before_frame.as_ref());
}

#[test]
fn g0_null_channels_and_transactional_work_band_capacity_refusals() {
    let gate = CancelGate::new_clock_free();
    let (original, properties) = reduced(false, true, false);
    let zero = original.with_thermoelastic_memory(&[properties], options(), &gate).unwrap();
    assert!(zero.relaxation_branches().is_empty() && zero.membrane_columns().is_empty());
    assert_eq!(zero.original().physical_modes(),original.physical_modes());
    let (original, properties) = reduced(false, false, false);
    for change in 0..5 {
        let mut limit = options();
        match change {
            0 => limit.max_projection_terms=0,
            1 => limit.max_factor_entries=0,
            2 => limit.max_memory_states=0,
            3 => limit.band_rad_s=(10_000.0,20_000.0),
            _ => {limit.max_thermal_terms=1;limit.normalized_tolerance=1e-12;},
        }
        assert!(original.with_thermoelastic_memory(std::slice::from_ref(&properties),limit,&gate).is_err());
    }
    let stop=CancelGate::new_clock_free();stop.request();
    assert!(original.with_thermoelastic_memory(std::slice::from_ref(&properties),options(),&stop).is_err());
    let attached=with_interface(original.source().clone(),(50.0,1e7),&gate);
    let attached=reduce_source(&attached,false,&gate);
    // Same physical attachment and pencil; only its source applicability edge changes.
    let maximum=attached.radiators().iter().map(|r|r.omega).fold(0.0_f64,f64::max);
    let narrow=with_interface(attached.source().clone(),(50.0,maximum*(1.0+1e-9)),&gate);
    let narrow=reduce_source(&narrow,false,&gate);
    assert!(narrow.with_thermoelastic_memory(std::slice::from_ref(&properties),options(),&gate).is_err());
    let image=original.with_thermoelastic_memory(&[properties],options(),&gate).unwrap();
    let mut out=[7.0;3];
    assert!(image.fill_force(&[0.1,f64::NAN],&mut out).is_err());
    assert_eq!(out,[7.0;3]);
    assert!(image.prepare_coupled(1,ModalAcousticTimeBudget::audible_reference(),
        config(),preparation(),&gate).is_err());
}
