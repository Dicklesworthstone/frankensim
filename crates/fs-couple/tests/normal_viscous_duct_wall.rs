//! Manufactured normal-plane materials; no measured instrument qualification.
use fs_couple::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall, WallMaterialLoss};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const BAND: [f64; 2] = [20.0, 4000.0];
const H: f64 = 0.0002;
const R: f64 = 0.006;
const LENGTH: f64 = 0.3;
const NAMES: [&str; 3] = ["eta_hh", "eta_ha", "eta_aa"];
const E: [f64; 3] = [4e6, 1e6, 3e6];
const NU: [f64; 3] = [0.2, 0.1, 0.25];
type Datum = (&'static str, Dims, f64);

fn values(d: [f64; 3]) -> Vec<Datum> {
    vec![("density", Density::DIMS, 1000.0),
        ("young_modulus_1", Pressure::DIMS, E[0]),
        ("young_modulus_2", Pressure::DIMS, E[1]),
        ("young_modulus_3", Pressure::DIMS, E[2]),
        ("poisson_ratio_12", Dims::NONE, NU[0]),
        ("poisson_ratio_13", Dims::NONE, NU[1]),
        ("poisson_ratio_23", Dims::NONE, NU[2]),
        (NAMES[0], DynViscosity::DIMS, d[0]),
        (NAMES[1], DynViscosity::DIMS, d[1]),
        (NAMES[2], DynViscosity::DIMS, d[2])]
}
fn state(data: &[Datum], band: [f64; 2]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in data {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0)
                .with("omega", band[0], band[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "synthetic normal-plane viscosity, not measured wood".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured-wall".into(), phase: "solid".into(),
        process: "independent tensor regression".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn wall(d: [f64; 3], axes: [usize; 2], axial: AxialConstraint) -> LocalMaterialWall {
    LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
        state(&values(d), [0.0, 1e6]), H, axial, axes, NAMES, BAND, 1.0, 0.8,
    ).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1e-20), "{a:e} != {b:e}");
}
// Independently invert the dimensional reciprocal normal compliance.
fn elastic(axes: [usize; 2]) -> [f64; 3] {
    let [h, a] = axes;
    let i = h.min(a); let j = h.max(a);
    let nu = NU[match (i,j) { (0,1) => 0, (0,2) => 1, _ => 2 }];
    let (hh, ha, aa) = (1.0/E[h], -nu/E[i], 1.0/E[a]);
    let determinant = hh*aa-ha*ha;
    [aa/determinant, -ha/determinant, hh/determinant]
}

#[test]
fn complete_complex_constitutive_schur_matches_the_same_native_wall_impedance() {
    let d = [400.0, -60.0, 90.0];
    for axes in [[0,1], [1,0], [0,2], [2,0], [1,2], [2,1]] {
        let [hh,ha,aa] = elastic(axes);
        for axial in [AxialConstraint::StressFree, AxialConstraint::StrainFixed] {
            let wall = wall(d, axes, axial);
            for f in [20.0, 80.0, 300.0, 1200.0, 4000.0] {
                let omega = core::f64::consts::TAU * f;
                // Direct division of (C_ha-i*w*D_ha)^2 by C_aa-i*w*D_aa;
                // this oracle has no Maxwell decomposition or fitted spectrum.
                let (re, im) = if axial == AxialConstraint::StressFree {
                    let (ar,ai) = (aa, -omega*d[2]);
                    let (nr,ni) = (ha*ha-omega*omega*d[1]*d[1], -2.0*ha*omega*d[1]);
                    let norm = ar*ar+ai*ai;
                    (hh-(nr*ar+ni*ai)/norm, -omega*d[0]-(ni*ar-nr*ai)/norm)
                } else { (hh, -omega*d[0]) };
                let actual = wall.specific_impedance(R, omega).unwrap();
                close(actual.re, -im/omega * H/(R*R), 2e-12);
                close(actual.im, re/omega * H/(R*R) - omega*1000.0*H, 2e-12);
            }
        }
    }
}

#[test]
fn zero_rank_deficient_and_proportional_viscosities_keep_their_actual_limits() {
    let zero = wall([0.0;3], [0,1], AxialConstraint::StressFree);
    assert!(matches!(zero.loss(), WallMaterialLoss::KelvinVoigt));
    assert_eq!(zero.base_patch(R,LENGTH).unwrap().wall.resistance, 0.0);
    let no_axial_rate = wall([400.0,0.0,0.0], [0,1], AxialConstraint::StressFree);
    close(no_axial_rate.base_patch(R,LENGTH).unwrap().wall.resistance, 400.0*H/(R*R), 2e-15);
    let rank_one = wall([400.0,200.0,100.0], [0,1], AxialConstraint::StressFree);
    assert!(matches!(rank_one.loss(), WallMaterialLoss::Maxwell(_)));
    assert_eq!(rank_one.base_patch(R,LENGTH).unwrap().wall.resistance, 0.0);
    let fixed = wall([400.0,-60.0,90.0], [0,1], AxialConstraint::StrainFixed);
    assert!(matches!(fixed.loss(), WallMaterialLoss::KelvinVoigt));
    close(fixed.base_patch(R,LENGTH).unwrap().wall.resistance, 400.0*H/(R*R), 2e-15);

    let c = elastic([0,1]);
    let tau = 1e-4;
    let proportional = wall(c.map(|x| tau*x), [0,1], AxialConstraint::StressFree);
    // Even if finite arithmetic retains a tiny positive arm, its complete
    // response must recover the proportional law, not another fitted damping.
    for f in [20.0,400.0,4000.0] {
        let omega = core::f64::consts::TAU*f;
        let z = proportional.specific_impedance(R,omega).unwrap();
        close(z.re, tau*E[0]*H/(R*R), 1e-12);
        close(z.im, E[0]*H/(R*R)/omega-omega*1000.0*H, 1e-12);
    }
}

#[test]
fn source_gaps_dimensions_passivity_bands_and_limits_refuse_without_defaults() {
    for d in [[-1.0,0.0,1.0],[1.0,0.0,-1.0],[1.0,2.0,1.0],[1.0,0.1,0.0],[0.0,-0.1,1.0]] {
        assert!(LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
            state(&values(d), [0.0,1e6]), H, AxialConstraint::StressFree,
            [0,1], NAMES, BAND, 1.0, 0.8).is_err());
    }
    let original = values([400.0,-60.0,90.0]);
    for missing in NAMES {
        let partial:Vec<_> = original.iter().copied().filter(|p|p.0!=missing).collect();
        assert!(LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
            state(&partial, [0.0,1e6]), H, AxialConstraint::StressFree,
            [0,1], NAMES, BAND, 1.0, 0.8).is_err());
    }
    let mut wrong = original.clone();
    wrong.last_mut().unwrap().1 = Pressure::DIMS;
    for source in [state(&wrong,[0.0,1e6]), state(&original,[0.0,2000.0])] {
        assert!(LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
            source,H,AxialConstraint::StressFree,[0,1],NAMES,BAND,1.0,0.8).is_err());
    }
    assert!(LocalMaterialWall::orthotropic_normal_viscosity_requirements(
        [0,1],["density",NAMES[1],NAMES[2]]).is_err());
    let source = state(&original,[0.0,1e6]);
    for (dt,angle) in [(0.0,0.8),(2.1,0.8),(1.0,0.0),(1.0,1.1)] {
        assert!(LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
            source.clone(),H,AxialConstraint::StressFree,[0,1],NAMES,BAND,dt,angle).is_err());
    }
}

fn derivative(y: [f64;3], c: [f64;3], d: [f64;3], beta: f64) -> [f64;3] {
    let [x,v,a] = y;
    let axial_rate = -(c[2]*a+c[1]*x/R+d[1]*v/R)/d[2];
    let hoop_stress = c[0]*x/R+c[1]*a+d[0]*v/R+d[1]*axial_rate;
    [v,(1.0-beta*x-H/R*hoop_stress)/(1000.0*H),axial_rate]
}
fn reference_step(y: [f64;3], c: [f64;3], d: [f64;3], beta: f64, dt: f64) -> [f64;3] {
    let add = |y: [f64;3], k: [f64;3], scale| core::array::from_fn(|i|y[i]+scale*k[i]);
    let k1=derivative(y,c,d,beta);
    let k2=derivative(add(y,k1,0.5*dt),c,d,beta);
    let k3=derivative(add(y,k2,0.5*dt),c,d,beta);
    let k4=derivative(add(y,k3,dt),c,d,beta);
    core::array::from_fn(|i|y[i]+dt/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]))
}

#[test]
fn original_two_strain_dynamics_match_native_gas_pressure_storage_and_loss() {
    let c=elastic([0,1]);
    let d=[400.0,-60.0,90.0];
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let bulk=gas.density*gas.sound_speed.powi(2);
    let beta=2.0*bulk/R;
    let geometry=fs_duct::Duct {segments:vec![fs_duct::Segment::Cylinder {radius:R,length:LENGTH}]};
    let duct=DgDuct::new(&geometry,&gas,DgMesh {elements:2,order:2},DgBoundary::Closed,DgBoundary::Closed).unwrap();
    let bound=wall(d,[0,1],AxialConstraint::StressFree).bind_dg(duct).unwrap();
    let duct=bound.duct();
    let propagator=duct.propagator(1e-6,1e-13).unwrap();
    let mut current=duct.project_initial_state(|_|1.0,|_|0.0);
    let initial=current.stored_energy_j();
    let mut previous=initial;
    let mut reference=[0.0;3];
    for _ in 0..128 {
        // Independent full [radial displacement, velocity, axial strain] law,
        // integrated on 32 finer substeps; not the condensed Maxwell equation.
        for _ in 0..32 { reference=reference_step(reference,c,d,beta,1e-6/32.0); }
        propagator.step(&mut current,&[]).unwrap();
        propagator.sync_nodal(&mut current).unwrap();
        let motion=duct.compliant_wall_state(&current).unwrap();
        assert_eq!(motion.relaxation_pressure_pa.len(),1);
        let [x,v,a]=reference;
        let pressure=1.0-beta*x;
        assert!((current.pressure()[0]-pressure).abs()<2e-8);
        assert!((motion.displacement_m[0]-x).abs()<1e-15);
        let strain=x/R;
        let wall_energy=core::f64::consts::TAU*R*LENGTH*H*0.5*
            (1000.0*v*v+c[0]*strain*strain+2.0*c[1]*strain*a+c[2]*a*a);
        close(motion.stored_energy_j,wall_energy,2e-7);
        let energy=core::f64::consts::PI*R*R*LENGTH*pressure*pressure/(2.0*bulk)+wall_energy;
        close(current.stored_energy_j(),energy,2e-8);
        let da=derivative(reference,c,d,beta)[2];
        let dh=v/R;
        let power=core::f64::consts::TAU*R*LENGTH*H*(d[0]*dh*dh+2.0*d[1]*dh*da+d[2]*da*da);
        close(motion.dissipated_power_w,power,2e-7);
        assert!(current.stored_energy_j()<=previous+1e-10*initial);
        previous=current.stored_energy_j();
    }
    assert!(current.stored_energy_j()<initial);
}

#[test]
fn cross_viscosity_changes_real_constitutive_memory_not_just_an_effective_decay() {
    let negative=wall([400.0,-60.0,90.0],[0,1],AxialConstraint::StressFree);
    let positive=wall([400.0,60.0,90.0],[0,1],AxialConstraint::StressFree);
    assert_eq!(negative.base_patch(R,LENGTH).unwrap().wall,positive.base_patch(R,LENGTH).unwrap().wall);
    let law = |wall: &LocalMaterialWall| match wall.loss() {
        WallMaterialLoss::Maxwell(spectrum)=>spectrum.material().terms.clone(),
        _=>panic!("independent normal coupling must retain its physical memory"),
    };
    assert_ne!(law(&negative),law(&positive));
    let omega=core::f64::consts::TAU*1000.0;
    assert_ne!(negative.specific_impedance(R,omega).unwrap(),positive.specific_impedance(R,omega).unwrap());
}
