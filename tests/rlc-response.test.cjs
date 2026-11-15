const test=require('node:test'),assert=require('node:assert/strict');
const {readFileSync}=require('node:fs'),{join}=require('node:path'),vm=require('node:vm');
const html=readFileSync(join(__dirname,'..','rlc_capacitor_response_static_slider.html'),'utf8');
const code=html.match(/<script id="rlc-math">([\s\S]*?)<\/script>/)[1];
const ctx=vm.createContext({});vm.runInContext(code,ctx);const M=ctx.RLCMath;
const base={f:1.5,R:3,L:.4,C:.015,V0:.45,I0:.08,output:'capacitor'};
const make=s=>({...base,...s});
const close=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<=tol*Math.max(1,Math.abs(b)),`${a} != ${b}`);
const outputs=['capacitor','resistor','inductor'];
const divide=(a,b)=>({re:(a.re*b.re+a.im*b.im)/(b.re*b.re+b.im*b.im),im:(a.im*b.re-a.re*b.im)/(b.re*b.re+b.im*b.im)});

test('all transfer functions match independently computed voltage-divider impedances',()=>{
 for(const f of [.2,1.5,2.05468148,12]) {
  const s=make({f}),w=2*Math.PI*f,zc={re:0,im:-1/(w*s.C)},zl={re:0,im:w*s.L},zr={re:s.R,im:0};
  const sum={re:s.R,im:zl.im+zc.im},p=M.parameters(s);
  for(const [out,z] of Object.entries({capacitor:zc,resistor:zr,inductor:zl})) {
   const expected=divide(z,sum);close(p.phasors[out].re,expected.re);close(p.phasors[out].im,expected.im);
   const model=M.model({...s,output:out});
   for(const t of [0,.1,.7])close(model.values(t).steady,expected.re*Math.sin(w*t)+expected.im*Math.cos(w*t));
  }
 }
});
test('selected outputs have the correct physical initial values',()=>{
 for(const V0 of [-2,.45,2])for(const I0 of [-.3,0,.3]) {
  const s=make({V0,I0});
  close(M.model({...s,output:'capacitor'}).values(0).total,V0);
  close(M.model({...s,output:'resistor'}).values(0).total,s.R*I0);
  close(M.model({...s,output:'inductor'}).values(0).total,-V0-s.R*I0);
 }
});
test('Kirchhoff voltage law holds for total and separate response contributions',()=>{
 for(const zeta of [.3,1,2]) {
  const s=make({R:2*Math.sqrt(base.L/base.C)*zeta});
  const models=outputs.map(output=>M.model({...s,output}));
  for(const t of [0,.001,.03,.2,1,10]) {
   const v=models.map(m=>m.values(t));
   close(v.reduce((a,b)=>a+b.total,0),Math.sin(2*Math.PI*s.f*t));
   close(v.reduce((a,b)=>a+b.steady,0),Math.sin(2*Math.PI*s.f*t));
   close(v.reduce((a,b)=>a+b.transient,0),0);
   close(v.reduce((a,b)=>a+b.initial,0),0);
  }
 }
});
test('both circuit differential equations hold in all damping regimes',()=>{
 for(const zeta of [.3,1,2]) {
  const s=make({R:2*Math.sqrt(base.L/base.C)*zeta});
  const c=M.model(s),r=M.model({...s,output:'resistor'});
  for(const t of [.03,.2,.7]) {
   const dt=1e-6,vc=c.values(t).total,i=r.values(t).total/s.R;
   const dv=(c.values(t+dt).total-c.values(t-dt).total)/(2*dt);
   const di=(r.values(t+dt).total-r.values(t-dt).total)/(2*dt*s.R);
   close(s.C*dv,i,2e-7);close(s.L*di+s.R*i+vc,Math.sin(2*Math.PI*s.f*t),2e-7);
  }
 }
});
test('matching both initial states cancels the transient for every output',()=>{
 for(const zeta of [.02,.3,1,2,20])for(const output of outputs) {
  const s=make({R:2*Math.sqrt(base.L/base.C)*zeta,output}),p=M.parameters(s);
  const m=M.model({...s,V0:p.matchV,I0:p.matchI});
  for(const t of [0,.001,.1,.6,5,15]) {
   const v=m.values(t);assert.equal(v.total,v.steady);close(v.transient+v.initial,0);
  }
 }
});
test('near-critical damping is continuous and finite on either side of repeated poles',()=>{
 const s=make({R:2*Math.sqrt(base.L/base.C)});
 for(const output of outputs) {
  const critical=M.model({...s,output});
  for(const factor of [1-1e-8,1+1e-8]) {
   const model=M.model({...s,output,R:s.R*factor});
   for(const t of [0,.01,.2,1,15])close(model.values(t).total,critical.values(t).total,1e-6);
  }
 }
});
test('matching states near resonance is not clipped to the manual voltage/current ranges',()=>{
 const s=make({f:2.9,R:.1,L:1.5,C:.002}),p=M.parameters(s);
 assert.ok(Math.abs(p.matchV)>2);assert.ok(Math.abs(p.matchI)>.3);
 for(const output of outputs) {
  const model=M.model({...s,output,V0:p.matchV,I0:p.matchI});
  for(const t of [0,.1,1,4])assert.equal(model.values(t).total,model.values(t).steady);
 }
});
