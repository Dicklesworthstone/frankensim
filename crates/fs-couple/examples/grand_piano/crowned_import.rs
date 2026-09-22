//! Preserve selected OBJ heights while reusing every existing part/material/
//! unit/frame/support/bridge admission. Flattening is only a temporary chart
//! for topology and point location; actual heights go to the structural shell.
use super::{Imported, Spec, MAX_OBJ_BYTES, MAX_NODES, MAX_TRIS, cross, dot, lower};
use super::super::crowned_board;
use fs_io::obj::read_obj_document;
use std::collections::{BTreeMap,BTreeSet};

pub fn import(obj:&str,specification:&str)->Result<Imported,String> {
    if obj.len()>MAX_OBJ_BYTES {return Err("OBJ exceeds 32 MiB".into());}
    let spec=Spec::read(specification)?;
    let mut document=read_obj_document(obj).map_err(|e|e.to_string())?;
    let mut selected=BTreeSet::new();let mut triangles=0usize;
    for region in &document.regions {
        if !region.has_label(&spec.part) {continue;}
        triangles=triangles.checked_add(region.triangles.len()).ok_or("triangle count overflow")?;
        if triangles>MAX_TRIS {return Err("selected crowned panel exceeds its triangle budget".into());}
        for face in region.triangles.clone() {
            for id in document.soup.triangles[face] {
                selected.insert(id as usize);
                if selected.len()>MAX_NODES {return Err("selected crowned panel exceeds its node budget".into());}
            }
        }
    }
    let normal=cross(spec.u,spec.v);
    let mut heights=BTreeMap::new();
    for id in selected {
        let p=&mut document.soup.positions[id];
        let d=[p.x-spec.origin[0],p.y-spec.origin[1],p.z-spec.origin[2]];
        let z=dot(d,normal)*spec.scale;
        if !z.is_finite() || z.abs()>0.05 {return Err("selected OBJ crown exceeds 50 mm from its declared reference plane".into());}
        heights.insert(id+1,z);
        let (u,v)=(dot(d,spec.u),dot(d,spec.v));
        let plane:[f64;3]=std::array::from_fn(|i|spec.origin[i]+u*spec.u[i]+v*spec.v[i]);
        p.x=plane[0];p.y=plane[1];p.z=plane[2];
    }
    // This preserves label assignments and polygon triangulation exactly. The
    // original importer still validates the numerical reference-plane chart;
    // its flatness tolerance is NOT used to erase or accept structural crown.
    let mut out=lower(&document,&spec)?;
    let z=out.source_vertices.iter().map(|id|heights.get(id).copied()
        .ok_or_else(||"selected crown vertex missing after chart remap".to_owned()))
        .collect::<Result<Vec<_>,_>>()?;
    out.fsb=crowned_board::elevate(&out.fsb,&z,"source OBJ heights retained in the declared SI frame; rectangular rib/bridge sections")?;
    out.max_projection_m=z.iter().map(|x|x.abs()).fold(0.,f64::max);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    const OBJ:&str="o ignored\nv 90 90 90\no soundboard\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0.5 0.5 0.015\nusemtl spruce\nf 2 3 6\nf 3 4 6\nf 4 5 6\nf 5 2 6\n";
    const SPEC:&str="frankensim-obj-board-v1\nsource,estimated,synthetic crown regression\npart,soundboard\nunits,1\nframe,0,0,0,1,0,0,0,1,0\nflatness,1e-8\nmaterial,spruce,0.008,450,1e10,8e8,0.3,6e8,0.2\nsupport,clamped\nboundary,all\ndamping,0.01\npretension,0\nbridge,69,0.5,0.5\n";
    #[test]
    fn source_crown_reaches_shell_instead_of_disappearing_in_projection() {
        assert!(super::super::import(OBJ,SPEC).is_err());
        let imported=import(OBJ,SPEC).unwrap();
        assert_eq!(imported.source_vertices,vec![2,3,4,5,6]);
        assert_eq!(imported.max_projection_m,0.015);
        let board=crowned_board::CrownedBoard::read(&imported.fsb).unwrap();
        assert_eq!(board.max_height_m,0.015);
        assert!(board.mass_kg>3.6);
        assert!(!board.prepare(&[69],400.).unwrap().modes.is_empty());
    }
    #[test]
    fn millimetres_and_rotated_asset_frames_preserve_structural_heights() {
        let obj="o soundboard\nv 10 20 30\nv 10 1020 30\nv 10 1020 1030\nv 10 20 1030\nv 25 520 530\nusemtl spruce\nf 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n";
        let spec=SPEC.replace("units,1\n","units,0.001\n")
            .replace("frame,0,0,0,1,0,0,0,1,0","frame,10,20,30,0,1,0,0,0,1");
        let a=import(OBJ,SPEC).unwrap();let b=import(obj,&spec).unwrap();
        assert_eq!(a.fsb,b.fsb);
    }
    #[test]
    fn shell_import_does_not_guess_missing_materials_or_downbearing() {
        assert!(import(OBJ,&SPEC.replace("material,spruce","material,maple")).is_err());
        assert!(import(OBJ,&SPEC.replace("pretension,0","pretension,100")).is_err());
        assert!(import(&OBJ.replace("0.015","0.06"),SPEC).is_err());
        assert!(import(OBJ,&SPEC.replace("support,clamped","support,simply_supported")).is_err());
    }
}
