//! Real material and circuit owners; synthetic coefficients, not measured fits.
use fs_couple::material_circuit::*;
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ResolvedConductor, resolve_uniform_conductor};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel, InterpolationPolicy};
use fs_evidence::ValidityDomain;
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance, TransformerLink};

fn conductor(r: f64) -> ResolvedConductor {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS),
        value:PropertyValue::Scalar { value:r, dims:ELECTRICAL_RESISTIVITY_DIMS },
        validity:ValidityDomain::unconstrained(), uncertainty:UncertaintyModel::Unstated,
        // A scalar claim is a plateau; matdb refuses LinearInside without knots.
        interpolation:InterpolationPolicy::ConstantWithinValidity, observations:vec![],
        provenance:Provenance { source:"synthetic circuit resistivity".into(),license:"CC0-1.0".into(),artifact:None },
    }).unwrap();
    let card=MaterialCard::assemble(MaterialStateId { chemistry:"synthetic circuit".into(),
        phase:"solid".into(),process:"authored numerical test".into(),revision:0 },claims,vec![]).unwrap();
    resolve_uniform_conductor(&card,&QueryPoint::new(),MaterialPropertySelection::SingleClaimOnly,1.0,1.0).unwrap()
}
fn material(branch:usize,r:f64) -> MaterialResistor { MaterialResistor {branch,conductor:conductor(r)} }
fn controls() -> MaterialCircuitSolve { MaterialCircuitSolve {energy_tolerance_j:1e-8,power_tolerance_w:1e-9} }
fn rc(r:f64) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:2, branches:vec![
        (1,0,Branch::Capacitor {farads:0.002}),
        (1,0,Branch::Resistor {ohms:999.0}),
    ],transformers:vec![]},vec![material(1,r)],vec![],vec![0.2],&[]).unwrap()
}
fn close(a:f64,b:f64) { assert!((a-b).abs()<=1e-8*b.abs().max(1.0),"{a} != {b}"); }

#[test]
fn finite_capacitor_discharge_matches_analytic_midpoint_and_branch_heat() {
    let mut circuit=rc(10.0);
    let initial=circuit.clone();
    let mut total_heat=0.0;
    let dt=0.0005_f64;
    let factor=(1.0-dt/(2.0*10.0*0.002))/(1.0+dt/(2.0*10.0*0.002));
    for step in 1..=64 {
        let proposal=circuit.step_candidate(circuit.resistors().to_vec(),&[],dt,controls(),||false).unwrap();
        let report=proposal.step();
        close(proposal.circuit().charges_c()[0],0.2*factor.powi(step));
        assert!(report.heat_j>0.0 && report.stored_after_j<report.stored_before_j);
        assert_eq!(report.native.supplied,0.0);
        assert_eq!(report.resistors[0].material.conductor.material(),circuit.resistors()[0].conductor.material());
        close(report.heat_j,report.stored_before_j-report.stored_after_j);
        close(report.resistors[0].current_a,report.resistors[0].voltage_v/10.0);
        total_heat+=report.heat_j;
        close(total_heat+report.stored_after_j,10.0);
        circuit=proposal.into_parts().0;
    }
    assert_eq!(initial.charges_c(),&[0.2]);
}

#[test]
fn every_resistor_gets_its_own_material_power_without_duplicating_supply() {
    let graph=CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:0.002}),
        (1,2,Branch::Resistor {ohms:1.0}),
        (2,0,Branch::Resistor {ohms:1.0}),
    ],transformers:vec![]};
    // Bind out of order deliberately; report remains in original branch order.
    let circuit=MaterialCircuit::new(graph,vec![material(2,30.0),material(1,10.0)],vec![],vec![0.2],&[]).unwrap();
    let candidate=circuit.step_candidate(circuit.resistors().to_vec(),&[],0.001,controls(),||false).unwrap();
    let heat=&candidate.step().resistors;
    assert_eq!(heat[0].material.branch,1); assert_eq!(heat[1].material.branch,2);
    close(heat[0].current_a,heat[1].current_a);
    close(heat[1].heat_j/heat[0].heat_j,3.0);
    close(candidate.step().heat_j,heat[0].heat_j+heat[1].heat_j);
    close(candidate.step().heat_j,candidate.step().native.dissipated);
}

#[test]
fn actual_transformer_and_driven_ports_retain_native_kirchhoff_work() {
    let graph=CircuitGraph {node_count:4,branches:vec![
        (1,0,Branch::VoltageSource {port:0}),
        (1,2,Branch::Resistor {ohms:1.0}),
        (3,0,Branch::Resistor {ohms:1.0}),
    ],transformers:vec![TransformerLink {primary:(2,0),secondary:(3,0),ratio:2.0}]};
    let circuit=MaterialCircuit::new(graph,vec![material(1,10.0),material(2,5.0)],vec![],vec![],&[0.0]).unwrap();
    let candidate=circuit.step_candidate(circuit.resistors().to_vec(),&[12.0],0.001,controls(),||false).unwrap();
    let report=candidate.step();
    close(report.native.y[0],12.0/(10.0+4.0*5.0));
    close(report.resistors[1].current_a,2.0*report.resistors[0].current_a);
    close(report.heat_j,report.native.supplied);
    assert_eq!(report.stored_after_j,0.0);
}

#[test]
fn resistance_only_swap_preserves_charge_but_changes_current_and_drain_rate() {
    let circuit=rc(10.0);
    let saved=circuit.clone();
    let a=circuit.step_candidate(vec![material(1,10.0)],&[],0.001,controls(),||false).unwrap();
    let b=circuit.step_candidate(vec![material(1,20.0)],&[],0.001,controls(),||false).unwrap();
    assert_eq!(a.step().stored_before_j,b.step().stored_before_j);
    assert!(b.circuit().charges_c()[0]>a.circuit().charges_c()[0]);
    assert!(b.step().heat_j<a.step().heat_j);
    assert!(b.step().resistors[0].current_a<a.step().resistors[0].current_a);
    assert_eq!(circuit,saved);
}

#[test]
fn inductor_storage_and_capacitor_charge_share_the_original_circuit() {
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Inductor {henries:0.05}),
        (1,0,Branch::Capacitor {farads:0.002}),
        (1,0,Branch::Resistor {ohms:10.0}),
    ],transformers:vec![]},vec![material(2,10.0)],vec![0.01],vec![0.02],&[]).unwrap();
    let candidate=circuit.step_candidate(circuit.resistors().to_vec(),&[],0.0005,controls(),||false).unwrap();
    assert_ne!(candidate.circuit().fluxes_wb(),circuit.fluxes_wb());
    assert_ne!(candidate.circuit().charges_c(),circuit.charges_c());
    close(candidate.step().stored_before_j,0.5*0.01_f64.powi(2)/0.05+0.5*0.02_f64.powi(2)/0.002);
    close(candidate.step().stored_after_j+candidate.step().heat_j,candidate.step().stored_before_j);
}

#[test]
fn bindings_nonfinite_ports_capacity_and_late_cancellation_refuse_atomically() {
    let circuit=rc(10.0); let saved=circuit.clone();
    for sources in [vec![],vec![material(0,10.0)],vec![material(99,10.0)],
        vec![material(1,10.0),material(1,20.0)]] {
        assert!(circuit.step_candidate(sources,&[],0.001,controls(),||false).is_err());
        assert_eq!(circuit,saved);
    }
    for dt in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(circuit.step_candidate(circuit.resistors().to_vec(),&[],dt,controls(),||false).is_err());
    }
    assert!(circuit.step_candidate(circuit.resistors().to_vec(),&[f64::NAN],0.001,controls(),||false).is_err());
    let mut count=0;
    circuit.step_candidate(circuit.resistors().to_vec(),&[],0.001,controls(),||{count+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(circuit.step_candidate(circuit.resistors().to_vec(),&[],0.001,controls(),||{
        actual+=1;actual==count
    }),Err(MaterialCircuitError::Cancelled)));
    assert_eq!(circuit,saved);
    let too_many=CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::VoltageSource {port:usize::MAX}),
    ],transformers:vec![]};
    assert!(MaterialCircuit::new(too_many,vec![],vec![],vec![],&[]).is_err());
}

#[test]
fn mutual_coils_retain_cross_energy_and_material_heat_when_resistance_changes() {
    // Two electrically isolated loops exchange energy through their common
    // field. The independent reference solves (L + dt R/2)i1=(L - dt R/2)i0.
    for m in [-0.75_f64, 0.75] {
        let graph = CircuitGraph { node_count:3, branches:vec![
            (1,0,Branch::Inductor {henries:2.0}),
            (1,0,Branch::Resistor {ohms:999.0}),
            (2,0,Branch::Resistor {ohms:999.0}),
            (2,0,Branch::Inductor {henries:3.0}),
        ], transformers:vec![] };
        let links = vec![MutualInductance { branches:(0,3), henries:m }];
        let mut current = [2.0_f64,-0.5_f64];
        let flux = |i:[f64;2]| vec![2.0*i[0]+m*i[1],m*i[0]+3.0*i[1]];
        let energy = |i:[f64;2]| 0.5*(2.0*i[0]*i[0]+2.0*m*i[0]*i[1]+3.0*i[1]*i[1]);
        let mut circuit = MaterialCircuit::new_with_mutual_inductance(graph,
            vec![material(1,1.0),material(2,3.0)], links.clone(), flux(current),vec![],&[]).unwrap();
        let initial_energy = energy(current);
        let mut total_heat = 0.0;
        for k in 0..12 {
            // Reconfiguration must retain M even when changing only contacts.
            circuit = circuit.rebind_contacts_candidate(vec![],&[]).unwrap();
            let r = [if k < 6 { 1.0 } else { 2.0 },3.0];
            let dt = 0.1;
            let a = 2.0+0.5*dt*r[0];
            let d = 3.0+0.5*dt*r[1];
            let rhs = [(2.0-0.5*dt*r[0])*current[0]+m*current[1],
                m*current[0]+(3.0-0.5*dt*r[1])*current[1]];
            let next = [(d*rhs[0]-m*rhs[1])/(a*d-m*m),
                (a*rhs[1]-m*rhs[0])/(a*d-m*m)];
            let original = circuit.clone();
            let proposal = circuit.step_candidate(vec![material(1,r[0]),material(2,r[1])],
                &[],dt,controls(),||false).unwrap();
            assert_eq!(proposal.circuit().mutual_inductance(),links);
            for (actual,expected) in proposal.circuit().fluxes_wb().iter().zip(flux(next)) {
                close(*actual,expected);
            }
            close(proposal.step().stored_before_j,energy(current));
            close(proposal.step().stored_after_j,energy(next));
            for j in 0..2 {
                let midpoint = 0.5*(current[j]+next[j]);
                close(proposal.step().resistors[j].current_a,-midpoint);
                close(proposal.step().resistors[j].heat_j,dt*r[j]*midpoint*midpoint);
            }
            total_heat += proposal.step().heat_j;
            close(energy(next)+total_heat,initial_energy);
            assert_eq!(circuit,original);
            circuit = proposal.into_parts().0;
            current = next;
        }
        let original = circuit.clone();
        assert!(matches!(circuit.step_candidate(circuit.resistors().to_vec(),&[],0.1,
            controls(),||true),Err(MaterialCircuitError::Cancelled)));
        assert_eq!(circuit,original);
    }
}

#[test]
fn material_circuit_rejects_an_indefinite_multiwinding_energy() {
    // Every |M| is below sqrt(Li Lj), but the three-winding matrix is indefinite.
    let graph = CircuitGraph { node_count:4, branches:vec![
        (1,0,Branch::Inductor {henries:1.0}),
        (2,0,Branch::Inductor {henries:1.0}),
        (3,0,Branch::Inductor {henries:1.0}),
    ], transformers:vec![] };
    let mutual = [(0,1),(0,2),(1,2)].into_iter().map(|branches|
        MutualInductance {branches,henries:-0.75}).collect();
    assert!(MaterialCircuit::new_with_mutual_inductance(graph,vec![],mutual,
        vec![0.0;3],vec![],&[]).is_err());
}
