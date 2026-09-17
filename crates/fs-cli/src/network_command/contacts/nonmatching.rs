//! Explicit contact-side declarations for independent volume meshes.
//! Geometry admission has finite pair/overlap budgets; it never searches for
//! unspecified partners. It runs during ordinary request parsing, before the
//! cooling solve's wall-time watchdog, just like existing mesh admission.
use super::*;
use fs_conduction::interface::{NonmatchingOptions, NonmatchingSurface};

pub(super) fn parse(value:&J,row:&Declaration,mesh:&ConductionMesh,
    slots:&BTreeMap<[u32;3],usize>,owned:&mut BTreeSet<[u32;3]>) -> Result<NonmatchingSurface> {
    object(value,&["side_a_faces","side_b_faces","plane_tolerance_m","coverage_relative_tolerance",
        "max_pair_tests","max_overlap_triangles"],"contact.nonmatching")?;
    let mut side=|key:&str|->Result<Vec<usize>> {
        let mut result=Vec::new();
        for item in array(get(value,key)?,key,mesh.boundary().len())? {
            let mut face=indices::<3>(item,key,mesh.vertex_count())?;
            face.sort_unstable();
            let slot=*slots.get(&face).ok_or_else(||bad("nonmatching contact side is not an exterior trace triangle"))?;
            if !owned.insert(face) {return Err(bad("nonmatching contact face is repeated or already externally owned"));}
            result.push(slot);
        }
        Ok(result)
    };
    let a=side("side_a_faces")?;let b=side("side_b_faces")?;
    let options=NonmatchingOptions {
        plane_tolerance_m:number(get(value,"plane_tolerance_m")?,"contact plane_tolerance_m")?,
        coverage_relative_tolerance:positive(get(value,"coverage_relative_tolerance")?,"contact coverage_relative_tolerance")?,
        max_pair_tests:count(get(value,"max_pair_tests")?,"contact max_pair_tests",1_000_000)?,
        max_overlap_triangles:count(get(value,"max_overlap_triangles")?,"contact max_overlap_triangles",100_000)?,
    };
    NonmatchingSurface::new(row.name.clone(),a,b,resistance(row)?,options).map_err(producer)
}

pub(super) fn bind(mesh:&ConductionMesh,boundary:&ThermalBoundary,matching:Vec<InterfaceSurface>,
    surfaces:Vec<NonmatchingSurface>)->Result<ThermalInterfaces> {
    if surfaces.is_empty(){return ThermalInterfaces::new(mesh,boundary,matching).map_err(producer);}
    // No random work is performed. This context owns bounded preprocessing,
    // not a fresh time allowance for any thermal solve or iterative producer.
    context(|cx|ThermalInterfaces::with_nonmatching(cx,mesh,boundary,matching,surfaces).map_err(producer))
}

pub(super) fn gradient(interfaces:&ThermalInterfaces,name:&str,temperature:&[f64],adjoint:&[f64])->Result<f64> {
    context(|cx|interfaces.nonmatching_log_resistance_pullback(cx,name,temperature,adjoint)
        .map_err(producer)?.ok_or_else(||bad("missing admitted nonmatching contact derivative")))
}

fn context<T>(run:impl FnOnce(&Cx<'_>)->Result<T>)->Result<T> {
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey{seed:0,kernel_id:717,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        run(&cx)
    })
}
