import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createObservatory} from '../server.mjs';
test('actual loopback HTTP serves the complete mesh module graph with strict headers',async()=>{
 const app=createObservatory({engine:'/nonexistent/native-fixture'});
 await new Promise(resolve=>app.server.listen(0,'127.0.0.1',resolve));
 const base=`http://127.0.0.1:${app.server.address().port}`;
 try{
  for(const asset of ['app.mjs','protocol.mjs','visuals.mjs','mesh-data.mjs','import-vtu.mjs','mesh-ui.mjs','mesh-view.mjs','mesh-raster.mjs','mesh.css']){
   const response=await fetch(`${base}/${asset}`);assert.equal(response.status,200,asset);
   assert.match(response.headers.get('content-type'),asset.endsWith('.css')?/text\/css/:/text\/javascript/);
   assert.match(response.headers.get('content-security-policy'),/script-src 'self'/);
   assert.ok((await response.text()).length>100);
  }
  assert.equal((await fetch(base+'/tools/vtu-to-observatory.mjs')).status,404);
  assert.equal((await fetch(base+'/mesh-data.mjs',{headers:{Origin:'https://example.com'}})).status,403);
 }finally{app.shutdown();}
});
