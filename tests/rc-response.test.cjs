const test = require('node:test');
const assert = require('node:assert/strict');
const {readFileSync} = require('node:fs');
const {join} = require('node:path');
const vm = require('node:vm');
const html=readFileSync(join(__dirname,'..','rc_sinusoidal_response_static_slider.html'),'utf8');
const source=html.match(/<script id="rc-math">([\s\S]*?)<\/script>/)[1];
const context=vm.createContext({});vm.runInContext(source,context);const M=context.RCMath;
const close=(a,b,tolerance=1e-11)=>assert.ok(Math.abs(a-b)<=tolerance*Math.max(1,Math.abs(b)),`${a} != ${b}`);

test('q is angular frequency times the RC time constant',()=>{
  close(M.parameters(1/(2*Math.PI*.4),.4).q,1);
  close(M.parameters(1.55,.46).q,2*Math.PI*1.55*.46);
});
test('the q=1 case requires an initial capacitor voltage of -0.5 V',()=>{
  const f=1/(2*Math.PI*.1),p=M.parameters(f,.1);
  close(p.cancellationVoltage,-.5);
  close(M.values(0,f,.1,-.5).steady,-.5);
});
test('manual initial voltage and the first-order circuit equation are preserved',()=>{
  for(const [f,tau] of [[.1,.02],[1.6,.1],[6,1]])for(const vc of [-2,0,1,2]) {
    close(M.values(0,f,tau,vc).total,vc);
    for(const t of [.03,.1,.7]) {
      const dt=1e-6;
      const derivative=(M.values(t+dt,f,tau,vc).total-M.values(t-dt,f,tau,vc).total)/(2*dt);
      close(tau*derivative+M.values(t,f,tau,vc).total,Math.sin(2*Math.PI*f*t),2e-7);
    }
  }
});
test('the exact cancellation voltage removes the net transient across parameter limits',()=>{
  for(const f of [.1,1.55,1.6,6])for(const tau of [.02,.1,.46,1]) {
    const p=M.parameters(f,tau);
    assert.ok(p.cancellationVoltage>=-.5&&p.cancellationVoltage<=0);
    for(const t of [0,.0001,.1,.8,5]) {
      const v=M.values(t,f,tau,p.cancellationVoltage);
      close(v.transient+v.initial,0);
      assert.equal(v.total,v.steady);
    }
  }
});
