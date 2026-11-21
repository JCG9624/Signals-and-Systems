const assert = require('node:assert/strict');
const { test } = require('node:test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const vm = require('node:vm');

const html = readFileSync(join(__dirname, '..', 'Bode_Pole_Zero_Explorer.html'), 'utf8');
const context = {};
vm.runInNewContext(html.match(/<script id="bode-math">([\s\S]*?)<\/script>/)[1], context);
const B = context.BodeMath;
const real = (id, kind, re) => ({ id, kind, form: 'real', re, im: 0 });
const pair = (id, kind, re, im) => ({ id, kind, form: 'pair', re, im });

function near(actual, expected, message) {
  if (!Number.isFinite(expected)) {
    assert.ok(Object.is(actual, expected), `${message}: ${actual} versus ${expected}`);
  } else {
    assert.ok(Math.abs(actual - expected) < 1e-9, `${message}: ${actual} versus ${expected}`);
  }
}

test('all non-origin factors have unit DC gain, including right-half-plane roots', () => {
  for (const kind of ['zero', 'pole']) {
    for (const root of [real(1, kind, -10), real(1, kind, 10), pair(1, kind, -.6, 2.9), pair(1, kind, .6, 2.9), pair(1, kind, 0, 3)]) {
      const f = B.normalisedFactor(root, 0);
      for (const key of ['mag', 'ma', 'phase', 'pa']) near(f[key], 0, `${kind}, ${root.re}, ${key} at DC`);
    }
  }
});

test('normalised first-order magnitude and phase match 1 ± s/a at the corner', () => {
  for (const kind of ['zero', 'pole']) {
    for (const re of [-10, 10]) {
      const sign = kind === 'zero' ? 1 : -1;
      const f = B.normalisedFactor(real(1, kind, re), 10);
      near(f.mag, sign * 10 * Math.log10(2), 'corner magnitude');
      near(f.ma, 0, 'corner magnitude asymptote');
      near(f.phase, sign * (re < 0 ? 45 : -45), 'corner phase');
      near(f.pa, f.phase, 'corner phase asymptote');
    }
  }
});

test('normalised quadratic factors match 1 + 2ζ(s/ωn) + (s/ωn)²', () => {
  for (const kind of ['zero', 'pole']) {
    for (const re of [-.6, .6]) {
      const root = pair(1, kind, re, Math.sqrt(8.64));
      const sign = kind === 'zero' ? 1 : -1;
      const f = B.normalisedFactor(root, 3);
      near(f.mag, sign * 20 * Math.log10(.4), 'quadratic magnitude at ωn');
      near(f.ma, 0, 'quadratic magnitude asymptote at ωn');
      near(f.phase, sign * (re < 0 ? 90 : -90), 'quadratic phase at ωn');
      near(f.pa, f.phase, 'quadratic phase asymptote at ωn');
    }
  }
});

test('normalised gain and contributions sum to the original exact and asymptotic totals', () => {
  const models = [
    [real(1, 'zero', -2), real(2, 'pole', -10), real(3, 'pole', -100)],
    [real(1, 'zero', 2), real(2, 'pole', 10), pair(3, 'pole', -.6, 2.9)],
    [real(1, 'zero', 2), real(2, 'zero', 3), pair(3, 'pole', .6, 2.9)],
    [real(1, 'zero', 0), real(2, 'pole', -10), pair(3, 'pole', .6, 2.9)],
    [real(1, 'pole', 0), real(2, 'pole', -10), pair(3, 'zero', 0, 3)],
    [real(1, 'zero', 2), real(2, 'pole', 2), real(3, 'pole', -10)],
    [real(1, 'zero', 0), real(2, 'pole', 0), real(3, 'pole', -10)],
    [],
  ];
  for (const roots of models) for (const gain of [-5000, 0, 5000]) {
    for (const omega of [.001, .1, 1, 2, 10, 100, 10000]) {
      const total = B.total(B.reduce(roots).roots, gain, omega);
      const sum = B.normalisedGain(roots, gain);
      for (const root of roots) {
        const f = B.normalisedFactor(root, omega);
        for (const key of ['mag', 'ma', 'phase', 'pa']) sum[key] += f[key];
      }
      for (const key of ['mag', 'ma', 'phase', 'pa']) near(sum[key], total[key], `sum ${key} at ${omega}`);
    }
  }
});

test('leading gain equals signed H(0) when there are no origin factors', () => {
  const roots = [real(1, 'zero', 2), real(2, 'pole', -10), pair(3, 'pole', -.6, Math.sqrt(8.64))];
  for (const gain of [-5000, 0, 5000]) near(B.normalisedGain(roots, gain).value, B.dcGain(roots, gain).value, 'signed DC gain');
  const twoRhpZeros = [real(1, 'zero', 2), real(2, 'zero', 3)];
  near(B.normalisedGain(twoRhpZeros, 1).phase, 360, 'continuous phase branch is retained');
});

test('origin factors and imaginary-axis singularities remain explicit', () => {
  for (const kind of ['pole', 'zero']) {
    const root = real(1, kind, 0);
    for (const omega of [.001, 1, 100]) for (const key of ['mag', 'ma', 'phase', 'pa']) {
      near(B.normalisedFactor(root, omega)[key], B.factor(root, omega)[key], 'origin factor unchanged');
    }
    near(B.normalisedGain([root, real(2, 'pole', -10)], 2).value, .2, 'origin excluded from leading constant');
    near(B.normalisedFactor(pair(2, kind, 0, 3), 3).phase, NaN, 'phase undefined at imaginary root');
  }
});

test('normalised contributions stay finite when the leading scalar overflows', () => {
  const roots = Array.from({ length: 20 }, (_, id) => pair(id + 1, 'pole', -1e-9, 1e-8));
  const gain = B.normalisedGain(roots, 1e100);
  assert.equal(gain.value, Infinity);
  assert.ok(Number.isFinite(gain.mag));
  assert.ok(Number.isFinite(gain.log10Abs));
});
