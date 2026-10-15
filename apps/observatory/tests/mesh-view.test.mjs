import {test} from 'node:test';
import assert from 'node:assert/strict';
import {projectPoint,pickTriangle} from '../mesh-view.mjs';
import {MESH_SCHEMA,validateMeshSnapshot,joinFrames} from '../mesh-data.mjs';
const snapshot=()=>({schema:MESH_SCHEMA,title:'Explicit test fixture',source:{model:'Test',material:'Synthetic',notes:['No solver result.']},mesh:{positions_m:[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],triangles:[[0,1,2],[0,1,3],[0,2,3],[1,2,3]],triangle_cells:[0,0,0,0],cell_count:1,source_cells:[[0,1,2,3]],source_types:[10],time_s:[0],fields:[{label:'MaterialIndex',unit:'1',association:'cell',components:1,values:[[3]]}]}});
test('camera projection retains independent XYZ, viewport aspect, zoom and pan',()=>{
 const camera={yaw:0,pitch:0,zoom:2,aspect:2,pan:[.1,.2]};
 assert.deepEqual(projectPoint([1,2,3],camera),[1.1,4.2,3]);
 const p=projectPoint([1,0,0],{...camera,yaw:Math.PI/2,zoom:1,aspect:1,pan:[0,0]});
 assert.ok(Math.abs(p[0])<1e-15);assert.equal(p[2],-1);
});
test('picking chooses the visible front triangle and a real source point, not interpolated values',()=>{
 const p=[[-1,-1,0],[1,-1,0],[0,1,0],[-1,-1,1],[1,-1,1],[0,1,1]],tri=[[0,1,2],[3,4,5]];
 assert.deepEqual(pickTriangle(p,tri,0,.8),{triangle:1,point:5});
 assert.equal(pickTriangle(p,tri,2,2),null);
 assert.equal(pickTriangle([[0,0,0],[0,0,0],[0,0,0]],[[0,1,2]],0,0),null);
});
test('cell identity survives serialization and sequence checks include interior connectivity',()=>{
 const a=snapshot(),b=snapshot();
 assert.deepEqual(validateMeshSnapshot(a).mesh.source_cells,a.mesh.source_cells);
 assert.deepEqual(joinFrames([a,b],[0,.1]).mesh.source_types,[10]);
 b.mesh.source_cells[0]=[0,2,1,3];
 assert.throws(()=>joinFrames([a,b],[0,.1]),/changed geometry/);
 delete a.mesh.source_cells;delete a.mesh.source_types;
 assert.throws(()=>joinFrames([a,a],[0,1]),/connectivity/);
});
test('surface/source mapping cannot label a triangle as belonging to an unrelated cell',()=>{
 const a=snapshot();a.mesh.cell_count=2;a.mesh.source_types=[5,10];a.mesh.source_cells=[[0,1,2],[0,1,2,3]];
 a.mesh.fields[0].values=[[1,2]];
 assert.throws(()=>validateMeshSnapshot(a),/does not belong/);
});

import {rasterize,ramp} from '../mesh-raster.mjs';
function triangle(z,value){return new Float32Array([-.8,-.8,z,value,1,0,0,0, .8,-.8,z,value,0,1,0,0, 0,.8,z,value,0,0,1,0]);}
test('software depth test preserves the front field regardless of draw order',()=>{
 const front=triangle(1,1),back=triangle(0,0),options={wire:false};
 const a=rasterize(new Float32Array([...front,...back]),p=>p,30,30,options);
 const b=rasterize(new Float32Array([...back,...front]),p=>p,30,30,options);
 assert.deepEqual(a,b);assert.deepEqual([...a.slice((15*30+15)*4,(15*30+15)*4+4)],[217,76,46,255]);
});
test('software cell colors remain discontinuous and outside geometry stays transparent',()=>{
 const a=rasterize(triangle(0,.5),p=>p,30,30,{wire:false});
 assert.deepEqual([...a.slice(0,4)],[0,0,0,0]);
 assert.deepEqual([...a.slice((15*30+15)*4,(15*30+15)*4+4)],new Uint8ClampedArray([...ramp(.5).map(v=>v*255),255]).reduce((a,v)=>(a.push(v),a),[]));
});
test('software renderer refuses an oversized target or excessive overlap before publication',()=>{
 assert.throws(()=>rasterize(triangle(0,0),p=>p,1000,1000),/pixel budget/);
 const values=new Float32Array(24*200);for(let i=0;i<200;i++)values.set(triangle(0,0),i*24);
 assert.throws(()=>rasterize(values,p=>p,600,500),/overlap/);
});
