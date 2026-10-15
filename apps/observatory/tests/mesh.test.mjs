import {test} from 'node:test';
import assert from 'node:assert/strict';
import {MESH_SCHEMA,validateMeshSnapshot,frameAt,fieldRange,probe,joinFrames} from '../mesh-data.mjs';
import {parseVTU,exterior} from '../import-vtu.mjs';
import {spawnSync} from 'node:child_process';
import {mkdtemp,writeFile,readFile,rm} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
export function fixture({type=10, points=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],cells=[[0,1,2,3]],field='310 320 330 340'}={}) {
  return `<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
<UnstructuredGrid><Piece NumberOfPoints="${points.length}" NumberOfCells="${cells.length}">
<Points><DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">${points.flat().join(' ')}</DataArray></Points>
<Cells><DataArray type="Int64" Name="connectivity" format="ascii">${cells.flat().join(' ')}</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">${cells.map((c,i)=>cells.slice(0,i+1).reduce((n,c)=>n+c.length,0)).join(' ')}</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">${cells.map(()=>type).join(' ')}</DataArray></Cells>
<PointData><DataArray type="Float64" Name="Temperature" NumberOfComponents="1" format="ascii" Unit="K">${field}</DataArray>
<DataArray type="Float64" Name="Heat &amp; flux" NumberOfComponents="3" format="ascii" Unit="W/m²">${points.map(()=>[3,4,0]).flat().join(' ')}</DataArray></PointData>
<CellData><DataArray type="Int32" Name="MaterialIndex" NumberOfComponents="1" format="ascii">${cells.map((_,i)=>i+7).join(' ')}</DataArray></CellData>
</Piece></UnstructuredGrid></VTKFile>`;
}
const imported=()=>parseVTU(fixture(),{lengthUnit:'m'});
test('read native writer layout: retain units, vectors, cell ownership and exact values',()=>{
 const s=imported(); assert.equal(s.schema,MESH_SCHEMA);assert.equal(s.mesh.triangles.length,4);
 assert.deepEqual(s.mesh.triangle_cells,[0,0,0,0]); assert.deepEqual(s.mesh.fields[0].values,[[310,320,330,340]]);
 assert.equal(s.mesh.fields[1].label,'Heat & flux'); assert.equal(s.mesh.fields[1].unit,'W/m²');
 assert.equal(s.mesh.fields[2].unit,'unspecified'); assert.deepEqual(fieldRange(s.mesh.fields[1]),[5,5]);
 assert.deepEqual(probe(s.mesh,1,2,'1').series[0].values,[4]);
});
test('coordinate unit is mandatory and only coordinates are converted',()=>{
 assert.throws(()=>parseVTU(fixture()),/Declare/);
 const s=parseVTU(fixture(),{lengthUnit:'mm'});assert.equal(s.mesh.positions_m[1][0],.001);
 assert.equal(s.mesh.fields[0].values[0][0],310);
});
test('joined tetrahedra remove only their common face and retain original cell indexes',()=>{
 const points=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]], cells=[[0,1,2,3],[0,2,1,4]];
 const s=parseVTU(fixture({points,cells,field:'1 2 3 4 5'}),{lengthUnit:'m'});
 assert.equal(s.mesh.triangles.length,6);assert.deepEqual(s.mesh.triangle_cells,[0,0,0,1,1,1]);
 assert.deepEqual(s.mesh.fields[2].values,[[7,8]]);
});
test('convex hexahedron and planar quad triangulate with source-cell mapping',()=>{
 const p=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
 assert.equal(exterior(p,[[0,1,2,3,4,5,6,7]],[12]).triangles.length,12);
 assert.deepEqual(exterior(p,[[0,1,2,3]],[9]).triangle_cells,[0,0]);
 const q=structuredClone(p);q[3][2]=.2;assert.throws(()=>exterior(q,[[0,1,2,3]],[9]),/Nonplanar/);
 assert.throws(()=>exterior(p,[[0,2,1,3]],[9]),/Nonconvex/);
});
test('explicit sequence clocks retain real frame values without interpolation',()=>{
 const a=imported(),b=imported();b.mesh.fields[0].values[0][0]=400;
 const s=joinFrames([a,b],[0,.3]);assert.deepEqual(probe(s.mesh,0,0).series[0].values,[310,400]);
 assert.equal(frameAt(s.mesh.time_s,.299),0);assert.equal(frameAt(s.mesh.time_s,.3),1);
 assert.equal(frameAt(s.mesh.time_s,-1),0);assert.equal(frameAt(s.mesh.time_s,10),1);
 assert.deepEqual(fieldRange(s.mesh.fields[0]),[310,400]);assert.deepEqual(a.mesh.time_s,[0]);
 assert.throws(()=>joinFrames([a,b],[0,0]),/strictly increasing/);
 b.mesh.positions_m[1][0]=2;assert.throws(()=>joinFrames([a,b],[0,1]),/changed geometry/);
});
test('explicit source units and association gate deformation',()=>{
 const s=imported();s.mesh.displacement_field=1;assert.throws(()=>validateMeshSnapshot(s),/metres/);
 s.mesh.fields[1].unit='m';assert.equal(validateMeshSnapshot(s).mesh.displacement_field,1);
});
for (const [name, change] of [
 ['nonfinite',x=>x.replace('310 320','NaN 320')],['truncated array',x=>x.replace('310 320','310')],
 ['binary',x=>x.replace('format="ascii"','format="binary"')],['unknown entity',x=>x.replace('&amp;','&remote;')],
 ['DTD',x=>'<!DOCTYPE VTKFile [<!ENTITY remote SYSTEM "file:///etc/passwd">]>'+x],
 ['duplicate attributes',x=>x.replace('Unit="K"','Unit="K" Unit="Pa"')],['bad XML',x=>x.replace('</Cells>','</PointData>')],
 ['huge count',x=>x.replace('NumberOfPoints="4"','NumberOfPoints="900000"')],
 ['unknown cell type',x=>x.replace('format="ascii">10<','format="ascii">42<')],
 ['bad connectivity',x=>x.replace('>0 1 2 3<','>0 1 2 999<')],
 ['unsafe integer',x=>x.replace('format="ascii">7<','format="ascii">9007199254740993<')]
])test(`VTU rejects ${name} atomically`,()=>assert.throws(()=>parseVTU(change(fixture()),{lengthUnit:'m'})));
test('invalid mesh fields, indexes and recorded topology refuse rather than repair',()=>{
 for(const mutate of [s=>s.mesh.triangles[0][0]=99,s=>s.mesh.fields[0].values[0].pop(),
   s=>s.mesh.time_s.push(1),s=>s.mesh.triangles.push(s.mesh.triangles[0]),s=>s.mesh.cell_count=0]){
   const s=imported();mutate(s);assert.throws(()=>validateMeshSnapshot(s));
 }
});
test('converter runs as shipped, preserves native values and never overwrites outputs',async()=>{
 const dir=await mkdtemp(join(tmpdir(),'observatory-vtu-test-'));
 try {
  const input=join(dir,'field.vtu'),output=join(dir,'result.json');await writeFile(input,fixture());
  const command=new URL('../tools/vtu-to-observatory.mjs',import.meta.url);
  const result=spawnSync(process.execPath,[command.pathname,'--units','m','--output',output,input],{encoding:'utf8'});
  assert.equal(result.status,0,result.stderr);assert.deepEqual(JSON.parse(await readFile(output,'utf8')).mesh,imported().mesh);
  assert.notEqual(spawnSync(process.execPath,[command.pathname,'--units','m','--output',output,input]).status,0);
 }finally{await rm(dir,{recursive:true,force:true});}
});
