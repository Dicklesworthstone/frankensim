import {test} from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp,writeFile,readFile} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {parseVTK,parseVTU} from '../import-vtu.mjs';
import {joinFrames,probe} from '../mesh-data.mjs';
// Exact cooling-cht field/coordinate structure, with synthetic values. This is
// exporter compatibility evidence, not execution of the Rust cooling solver.
export function fixture({axes=[[0,1,3],[0,2,5],[0,4,9]],origin=[0,0,0],units=true}={}) {
  const extent=axes.flatMap((a,i)=>[origin[i],origin[i]+a.length-1]).join(' ');
  const cells=axes.reduce((n,a)=>n*Math.max(1,a.length-1),1),points=axes.reduce((n,a)=>n*a.length,1);
  const arr=(name,values,components=1,type='Float64',unit='')=>`<DataArray type="${type}" Name="${name}" NumberOfComponents="${components}"${units&&unit?` Unit="${unit}"`:''} format="ascii">${values.join(' ')}</DataArray>`;
  return `<?xml version="1.0"?>\n<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
<RectilinearGrid WholeExtent="${extent}"><Piece Extent="${extent}">
<PointData>${arr('displacement',Array.from({length:points*3},(_,i)=>i*.001),3,'Float64','m')}</PointData>
<CellData Scalars="temperature_k" Vectors="velocity_m_s">
${arr('temperature_k',Array.from({length:cells},(_,i)=>300+i),1,'Float64','K')}
${arr('velocity_m_s',Array.from({length:cells*3},(_,i)=>i-7),3,'Float64','m/s')}
${arr('pressure_pa',Array.from({length:cells},(_,i)=>i*100),1,'Float64','Pa')}
${arr('material',Array.from({length:cells},(_,i)=>i%3),1,'Int32')}
</CellData><Coordinates>${axes.map((a,i)=>arr('xyz'[i],a)).join('')}</Coordinates>
</Piece></RectilinearGrid></VTKFile>`;
}
test('VTR lowers nonuniform axes to original x-fastest cells with exact associated data',()=>{
  const s=parseVTK(fixture(),{lengthUnit:'m'}),m=s.mesh;
  assert.equal(m.positions_m.length,27);assert.equal(m.cell_count,8);assert.equal(m.triangles.length,48);
  assert.deepEqual(m.positions_m[14],[3,2,4]);
  assert.deepEqual(m.source_cells[0],[0,1,4,3,9,10,13,12]);
  assert.deepEqual(m.source_cells[7],[13,14,17,16,22,23,26,25]);
  assert.deepEqual(m.fields.map(f=>[f.label,f.unit,f.association]),[
    ['displacement','m','point'],['temperature_k','K','cell'],['velocity_m_s','m/s','cell'],['pressure_pa','Pa','cell'],['material','unspecified','cell']]);
  assert.equal(probe(m,1,7).series[0].values[0],307);
  assert.equal(probe(m,2,7,'2').series[0].values[0],16);
  assert.equal(m.displacement_field,undefined); // A field name never enables deformation.
  assert.deepEqual(parseVTU(fixture(),{lengthUnit:'m'}),s); // Legacy API stays compatible.
});
test('nonzero and negative extent origins do not shift or permute coordinates',()=>{
  const a=parseVTK(fixture(),{lengthUnit:'m'}),b=parseVTK(fixture({origin:[-12,15,6]}),{lengthUnit:'m'});
  assert.deepEqual(a.mesh,b.mesh);
});
test('coordinate scaling never rescales field values or invents missing field units',()=>{
  const a=parseVTK(fixture({units:false}),{lengthUnit:'mm'}).mesh;
  assert.deepEqual(a.positions_m.at(-1),[3,5,9].map(v=>v*.001));
  assert.equal(a.fields[1].values[0][7],307);assert.equal(a.fields[1].unit,'unspecified');
  assert.throws(()=>parseVTK(fixture()),/coordinate unit/);
});
for(const fixed of [0,1,2])test(`planar grid on axis ${fixed} retains all quad source cells`,()=>{
  const axes=[[0,1,3],[0,2,5],[0,4,9]];axes[fixed]=[7];
  const m=parseVTK(fixture({axes}),{lengthUnit:'cm'}).mesh;
  assert.equal(m.cell_count,4);assert.equal(m.triangles.length,8);
  assert.ok(m.source_types.every(t=>t===9));assert.ok(m.positions_m.every(p=>p[fixed]===.07));
});
test('recorded sequence keeps topology and uses supplied physical times',()=>{
  const a=parseVTK(fixture(),{lengthUnit:'m'}),b=structuredClone(a);b.mesh.fields[1].values[0][0]=900;
  const s=joinFrames([a,b],[.1,.5]);
  assert.deepEqual(probe(s.mesh,1,0).series[0].values,[300,900]);
  b.mesh.positions_m[0][0]=-.1;assert.throws(()=>joinFrames([a,b],[0,1]),/changed geometry/);
});
for(const [name,change,pattern] of [
  ['multiple pieces',s=>s.replace('</RectilinearGrid>','<Piece Extent="0 2 0 2 0 2"/></RectilinearGrid>'),/one Piece/],
  ['partial pieces',s=>s.replace('Piece Extent="0 2','Piece Extent="0 1'),/Partial/],
  ['reversed extent',s=>s.replaceAll('0 2 0 2 0 2','2 0 0 2 0 2'),/extent/],
  ['grid budget',s=>s.replaceAll('0 2 0 2 0 2','0 49999 0 49999 0 49999'),/budget/],
  ['coordinate order',s=>s.replace('>0 1 3</DataArray>','>0 3 1</DataArray>'),/strictly increasing/],
  ['duplicate axis',s=>s.replace('>0 1 3</DataArray>','>0 1 1</DataArray>'),/strictly increasing/],
  ['coordinate arity',s=>s.replace('>0 1 3</DataArray>','>0 1</DataArray>'),/expected 3/],
  ['field arity',s=>s.replace('300 301 302 303 304 305 306 307','300 301'),/expected 8/],
  ['compressed',s=>s.replace('type="RectilinearGrid"','type="RectilinearGrid" compressor="vtkZLibDataCompressor"'),/uncompressed/],
  ['binary',s=>s.replace('format="ascii"','format="binary"'),/ASCII/],
  ['nonfinite field',s=>s.replace('300 301','NaN 301'),/numeric/],
  ['external entities',s=>'<!DOCTYPE x [<!ENTITY x SYSTEM "file:///etc/passwd">]>'+s,/DTD/],
])test(`VTR refuses ${name} rather than publishing a partial mesh`,()=>assert.throws(()=>parseVTK(change(fixture()),{lengthUnit:'m'}),pattern));
test('actual converter reads VTR sequences and never overwrites an output',async()=>{
  const dir=await mkdtemp(join(tmpdir(),'observatory-vtr-')),a=join(dir,'a.vtr'),b=join(dir,'b.vtr'),out=join(dir,'sequence.json');
  await writeFile(a,fixture());await writeFile(b,fixture().replace('300 301','310 311'));
  const tool=fileURLToPath(new URL('../tools/vtu-to-observatory.mjs',import.meta.url));
  const args=[tool,'--units','m','--output',out,'--times','0,2',a,b];
  const run=spawnSync(process.execPath,args,{encoding:'utf8'});assert.equal(run.status,0,run.stderr);
  const s=JSON.parse(await readFile(out,'utf8'));assert.deepEqual(s.mesh.time_s,[0,2]);
  assert.deepEqual(probe(s.mesh,1,0).series[0].values,[300,310]);
  assert.notEqual(spawnSync(process.execPath,args).status,0);
});
